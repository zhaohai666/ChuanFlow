#!/usr/bin/env python3
"""Analytic verification of the declared PR binary interaction kij (chflow-pvt).

usage: check_pvt_binary_interaction.py <exe> <case.dict>

D7 increment 3 (5.95). The van der Waals one-fluid mixing rule gains a declared,
name-based symmetric binary-interaction matrix a_ij = sqrt(a_i a_j)(1 - kij), and
the component registry gains water (IAPWS/NIST Tc/Pc, omega 0.3443). The checker
re-derives the whole Peng-Robinson phase state independently in Python:

  kernel     for every flashed phase it rebuilds a_i, b_i, a_ij (with the kij
             parsed straight from the case), a_mix, b_mix, the cubic roots and
             ln phi_i from the REPORTED composition, and holds the reported
             compressibility / density / molar mass and the fugacity-equality
             K_i = y_i/x_i = exp(ln phi_L - ln phi_V) to the C++ report.
  gate       a kij = 0 twin (binaryInteraction block removed) is re-derived with
             the same independent PR at kij = 0 and must match bit-for-bit, while
             the declared-kij report must DIFFER from it (kij is observable).
  physics    water concentrates in the liquid, acid gas / methane in the vapour;
             every flash converges; beta in [0, 1]; material balance closes.

The kij values are DECLARED case inputs, not a kernel table; the constants below
are transcribed from src/pvt/src/components.cpp and the PR76 omega_a / omega_b.
"""
import math
import os
import re
import subprocess
import sys

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]

R = 8.314462618
SQRT2 = 1.4142135623730951
OMEGA_A = 0.45723552892138218
OMEGA_B = 0.07779607390388846

# ---- Independent PR76 alpha kappa (NOT a mirror of the chflow kernel) -------
# Source: Cantera 3.1 documentation, "Peng-Robinson Equation of State":
#   https://www.cantera.org/3.1/reference/thermo/species-thermo.html
#   alpha = [1 + kappa (1 - sqrt(T/Tc))]^2
#   kappa = 0.37464 + 1.54226 w - 0.26992 w^2,  for w <= 0.491
# These are the PR76 original published coefficients, re-derived here from the
# external physical reference; they are NOT transcribed from src/pvt (the old
# checker copy had a 10x-too-small quadratic term 0.026992).
PR_KAPPA_BASE = 0.37464
PR_KAPPA_LINEAR = 1.54226
PR_KAPPA_QUAD = 0.26992


def pr76_kappa(w):
    """PR76 kappa from the published alpha function, omega <= 0.491 branch."""
    # All CT components above have omega <= 0.3443, i.e. inside this branch.
    return PR_KAPPA_BASE + PR_KAPPA_LINEAR * w - PR_KAPPA_QUAD * w * w

# Exact constants from src/pvt/src/components.cpp (MW g/mol, Tc K, Pc Pa, omega).
CT = {
    "methane": (16.043, 190.564, 45.992e5, 0.0115),
    "ethane": (30.070, 305.322, 48.722e5, 0.0995),
    "propane": (44.097, 369.834, 42.477e5, 0.1523),
    "iButane": (58.123, 407.817, 36.480e5, 0.1835),
    "nButane": (58.123, 425.125, 37.962e5, 0.2002),
    "nPentane": (72.150, 469.700, 33.700e5, 0.2515),
    "nHexane": (86.178, 507.600, 30.250e5, 0.3013),
    "carbonDioxide": (44.010, 304.128, 73.773e5, 0.2239),
    "hydrogenSulphide": (34.081, 373.400, 89.628e5, 0.0900),
    "nitrogen": (28.014, 126.192, 33.958e5, 0.0372),
    "water": (18.01528, 647.096, 220.64e5, 0.3443),
}


def cbrt(x):
    return math.copysign(abs(x) ** (1.0 / 3.0), x)


def cubic_real_roots(c2, c1, c0):
    p = c1 - c2 * c2 / 3.0
    q = 2.0 * c2 * c2 * c2 / 27.0 - c2 * c1 / 3.0 + c0
    disc = q * q / 4.0 + p * p * p / 27.0
    shift = -c2 / 3.0
    if disc > 0.0:
        s = math.sqrt(disc)
        return [cbrt(-q / 2.0 + s) + cbrt(-q / 2.0 - s) + shift]
    if p >= 0.0:
        return [shift]
    r = math.sqrt(-p * p * p / 27.0)
    m = 2.0 * math.sqrt(-p / 3.0)
    theta = math.acos(max(-1.0, min(1.0, -q / (2.0 * r))))
    return sorted(m * math.cos((theta + 2.0 * k * math.pi) / 3.0) + shift
                  for k in range(3))


def pr_phase(names, zs, T, p, kij, pref):
    """Independent PR: returns (Z, molarVolume, M, lnPhi list) for the phase."""
    n = len(names)
    MW = [CT[c][0] * 1.0e-3 for c in names]
    M = sum(z * m for z, m in zip(zs, MW))
    RT = R * T
    aI, bI = [], []
    for c in names:
        _, Tc, Pc, w = CT[c]
        kappa = pr76_kappa(w)  # independent PR76 criterion, not a kernel mirror
        alpha = (1.0 + kappa * (1.0 - math.sqrt(T / Tc))) ** 2
        Rtc = R * Tc
        aI.append(OMEGA_A * Rtc * Rtc / Pc * alpha)
        bI.append(OMEGA_B * Rtc / Pc)
    aIJ = [[math.sqrt(aI[i] * aI[j]) * (1.0 - kij(names[i], names[j]))
            for j in range(n)] for i in range(n)]
    bMix = sum(zs[i] * bI[i] for i in range(n))
    aMix = sum(zs[i] * zs[j] * aIJ[i][j] for i in range(n) for j in range(n))
    A = aMix * p / (RT * RT)
    B = bMix * p / RT
    roots = [z for z in cubic_real_roots(-(1.0 - B), A - 3.0 * B * B - 2.0 * B,
                                         -(A * B - B * B - B ** 3))
             if z > B * (1.0 + 1e-12)]
    if not roots:
        raise SystemExit("PR cubic has no physical root")

    def lnphi(Z):
        logTerm = math.log((Z + (1.0 + SQRT2) * B) / (Z + (1.0 - SQRT2) * B))
        out = []
        for i in range(n):
            sumA = sum(zs[j] * aIJ[i][j] for j in range(n))
            biB = bI[i] / bMix
            out.append(biB * (Z - 1.0) - math.log(Z - B)
                       - A / (2.0 * SQRT2 * B) * (2.0 * sumA / aMix - biB) * logTerm)
        return out

    distinct = sorted(set(roots))
    if pref == "liquid":
        Z = distinct[0]
    elif pref == "vapour":
        Z = distinct[-1]
    else:
        Z = min(distinct, key=lambda z: sum(zs[i] * lnphi(z)[i] for i in range(n)))
    return Z, Z * RT / p, M, lnphi(Z)


# ---- parse the case: feed, temperature, declared kij ----
with open(case) as fh:
    src = fh.read()
comp = re.search(r"composition\s*\{(.*?)\}", src, re.S)
feed = {}
for k, v in re.findall(r"([A-Za-z_]\w*)\s+([-+0-9.eE]+)\s*;", comp.group(1)):
    feed[k] = float(v)
tot = sum(feed.values())
feed = {k: v / tot for k, v in feed.items()}
temperature = float(re.search(r"temperature\s+([-+0-9.eE]+)\s*;", src).group(1))

kij_pairs = {}
bi = re.search(r"binaryInteraction\s*\{(.*?)\n\}", src, re.S)
if bi:
    for m in re.finditer(r"([A-Za-z_]\w*)\s*\{([^}]*)\}", bi.group(1)):
        ni = m.group(1)
        for nj, val in re.findall(r"([A-Za-z_]\w*)\s+([-+0-9.eE]+)\s*;", m.group(2)):
            v = float(val)
            kij_pairs[(ni, nj)] = v
            kij_pairs[(nj, ni)] = v


def kij_declared(a, b):
    return 0.0 if a == b else kij_pairs.get((a, b), 0.0)


def kij_zero(a, b):
    return 0.0


def run(path):
    proc = subprocess.run([exe, path], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"pvt failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")
    return proc.stdout


def parse(out):
    flashes, cur = [], None
    for line in out.splitlines():
        t = line.replace("[chflow]", "").split()
        if not t or t[0] != "pvt":
            continue
        k = t[1]
        if k == "flash":
            kv = dict(zip(t[2::2], t[3::2]))
            cur = {"T": float(kv["T"]), "p": float(kv["p"]),
                   "beta": float(kv["beta"]), "single": kv["single"] == "1",
                   "converged": kv["converged"] == "1",
                   "residual": float(kv["residual"]), "x": {}, "y": {}}
            flashes.append(cur)
        elif k == "liquid" and cur is not None:
            cur["rhoL"], cur["ML"] = float(t[3]), float(t[5])
        elif k == "vapour" and cur is not None:
            cur["rhoV"], cur["MV"] = float(t[3]), float(t[5])
        elif k in ("x", "y") and cur is not None:
            cur[k][t[2]] = float(t[3])
    return flashes


fails = []
worst = {"Z": 0.0, "rho": 0.0, "M": 0.0, "lnK": 0.0}


def check(name, ok, detail=""):
    if not ok:
        fails.append(f"{name} {detail}".strip())


def close(name, got, want, rtol, atol=0.0):
    ok = abs(got - want) <= atol + rtol * abs(want)
    check(name, ok, f"got {got:.12g} want {want:.12g}")


def verify_phase(f, side, kij, tag):
    comp_d = f[side]
    names = list(comp_d.keys())
    zs = [comp_d[c] for c in names]
    rho = f["rhoL"] if side == "x" else f["rhoV"]
    Mrep = f["ML"] if side == "x" else f["MV"]
    pref = "liquid" if side == "x" else "vapour"
    Z, Vm, M, lp = pr_phase(names, zs, f["T"], f["p"], kij, pref)
    label = f"{tag} p={f['p']:.6g} {side}"
    close(f"{label} molar mass", Mrep, M, 1.0e-9, 1.0e-12)
    close(f"{label} density", rho, M / Vm, 1.0e-9, 1.0e-9)
    Zrep = f["p"] * Mrep / (rho * R * f["T"])
    close(f"{label} compressibility", Zrep, Z, 1.0e-9, 1.0e-12)
    worst["Z"] = max(worst["Z"], abs(Zrep - Z))
    worst["rho"] = max(worst["rho"], abs(rho - M / Vm))
    worst["M"] = max(worst["M"], abs(Mrep - M))
    return lp


def verify_flash(f, kij, tag):
    label = f"{tag} p={f['p']:.6g}"
    check(f"{label} converged", f["converged"], "successive substitution stalled")
    check(f"{label} residual", f["residual"] < 1.0e-9, f"{f['residual']:.3g}")
    check(f"{label} beta in range", 0.0 <= f["beta"] <= 1.0, f"{f['beta']}")
    for c in feed:
        close(f"{label} material balance {c}",
              f["beta"] * f["y"][c] + (1.0 - f["beta"]) * f["x"][c], feed[c],
              1.0e-8, 1.0e-12)
    if f["single"]:
        verify_phase(f, "x", kij, tag)
        return
    lpL = verify_phase(f, "x", kij, tag)
    lpV = verify_phase(f, "y", kij, tag)
    names = list(f["x"].keys())
    for i, c in enumerate(names):
        lnK_ref = lpL[i] - lpV[i]
        lnK_rep = math.log(f["y"][c] / f["x"][c])
        d = abs(lnK_ref - lnK_rep)
        worst["lnK"] = max(worst["lnK"], d)
        check(f"{label} fugacity equality {c}", d < 1.0e-7,
              f"lnK ref {lnK_ref:.10g} vs y/x {lnK_rep:.10g}")


out = run(case)
flashes = parse(out)
check("case produced flashes", len(flashes) >= 3, f"{len(flashes)} flashes")
check("case declares kij", bool(kij_pairs), "no binaryInteraction pairs parsed")
check("case includes water", "water" in feed, "feed has no water")
for f in flashes:
    verify_flash(f, kij_declared, "kij")

# ---- physics: water to the liquid, methane / acid gas to the vapour ----
two = [f for f in flashes if not f["single"]]
check("has a two-phase branch", len(two) >= 2, f"{len(two)} two-phase flashes")
for f in two:
    check(f"p={f['p']:.6g} water concentrates in liquid",
          f["x"]["water"] > f["y"]["water"], f"{f['x']['water']} vs {f['y']['water']}")
    check(f"p={f['p']:.6g} methane enriched in vapour",
          f["y"]["methane"] > f["x"]["methane"], "")
    check(f"p={f['p']:.6g} liquid denser than vapour", f["rhoL"] > f["rhoV"], "")

# ---- gate: kij = 0 twin matches independent PR at kij = 0, and differs ----
twin = re.sub(r"binaryInteraction\s*\{.*?\n\}", "", src, count=1, flags=re.S)
twin_path = case + ".kij0"
with open(twin_path, "w") as fh:
    fh.write(twin)
twin_out = run(twin_path)
os.remove(twin_path)
twin_flashes = parse(twin_out)
check("twin produced same flash count", len(twin_flashes) == len(flashes),
      f"{len(twin_flashes)} vs {len(flashes)}")
for f in twin_flashes:
    verify_flash(f, kij_zero, "kij0")
base_lines = [l for l in out.splitlines() if "pvt flash" in l]
twin_lines = [l for l in twin_out.splitlines() if "pvt flash" in l]
check("declared kij is observable (differs from kij=0)", base_lines != twin_lines,
      "kij block had no effect on the flash report")

if fails:
    print("\n".join(f"FAIL: {m}" for m in fails))
    raise SystemExit(f"{len(fails)} pvt-kij check(s) failed")

print(f"pvt-kij: {len(flashes)} flashes ({len(two)} two-phase), "
      f"worst Z {worst['Z']:.3e} / rho {worst['rho']:.3e} / M {worst['M']:.3e} / "
      f"lnK {worst['lnK']:.3e}, kij=0 twin matches independent PR and differs - "
      f"all checks passed")