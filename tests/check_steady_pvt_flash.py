#!/usr/bin/env python3
"""Analytic verification of per-station true-flash PVT march in 1D (D10/D7, 5.121).

usage: check_steady_pvt_flash.py <exe> <case.dict>

The `compositional{ liveFlash true; sigmaFromP true; }` opt-in makes the 1D
PVT march flash the feed at EVERY station (PR76 successive-substitution TP
flash + homogeneous mixture density + mass-weighted LBC viscosity + Weinaug-
Katz surface tension) instead of reading a precomputed rho/mu table.

This checker never trusts the kernel's property functions. It:
  flash    re-derives each station's (rho, mu, sigma) with an independent
           Python transcription of PR76 TP flash, LBC and Weinaug-Katz,
           and compares to the printed per-station readout to < 1e-8;
  march    replays the 24-segment energy balance with those independent
           properties and checks the reported mass rate, outlet velocity,
           Reynolds, friction factor and downstream pressure;
  sigma    checks the WK sigma really varies along the march (it is a
           per-station quantity, not a constant) and is the value the march
           reported at the outlet;
  converge runs plain rho(p)/mu(p) twins built from the SAME independent
           flash on a coarse and a fine pressure grid and checks the fine
           grid converges to the live-flash result (interpolation error,
           not a different physics);
  gate     a no-declaration case and an explicit `liveFlash false;
           sigmaFromP false;` case are byte-identical, and neither leaks a
           flashStation/sigmaOutlet line (zero-declaration iron law);
  refuse   malformed declarations fail loudly.

Pure Python standard library only.
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]

R = 8.314462618
SQRT2 = 1.4142135623730951
OMEGA_A = 0.45723552892138218
OMEGA_B = 0.07779607390388846
ZC = 0.3074
K_CRIT_OVER_B = 0.30740 / 0.07779607390388846
ATM = 101325.0

# Exact component constants from src/pvt/src/components.cpp:
# name -> (MW g/mol, Tc K, Pc Pa, acentric).
CT = {
    "methane": (16.043, 190.564, 45.992e5, 0.0115),
    "ethane": (30.070, 305.322, 48.722e5, 0.0995),
    "propane": (44.097, 369.834, 42.477e5, 0.1523),
    "nButane": (58.123, 425.125, 37.962e5, 0.2002),
}
# Weinaug & Katz (1943) cgs parachors from src/pvt/src/componentData.cpp.
PARA_CGS = {
    "methane": 77.0, "ethane": 108.0, "propane": 150.3, "nButane": 189.9,
}

fails = []


def check(name, ok, detail=""):
    if not ok:
        fails.append(f"{name} {detail}".strip())


def close(name, got, want, rtol, atol=0.0):
    ok = abs(got - want) <= atol + rtol * abs(want)
    check(name, ok, f"got {got:.12g} want {want:.12g}")


def cbrt(x):
    return math.copysign(abs(x) ** (1.0 / 3.0), x)


def cubic_real_roots(c2, c1, c0):
    """Analytic roots, termwise identical to pengRobinson.cpp's solver."""
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


def wilson_k(c, T, p):
    _, Tc, Pc, w = CT[c]
    return Pc / p * math.exp(5.373 * (1.0 + w) * (1.0 - Tc / T))


def rr(z, K, beta):
    f = 0.0
    for zi, ki in zip(z, K):
        d = ki - 1.0
        f += zi * d / (1.0 + beta * d)
    return f


def solve_rr(z, K):
    lo, hi, beta = 0.0, 1.0, 0.5
    for _ in range(100):
        f = rr(z, K, beta)
        if f > 0.0:
            lo = beta
        else:
            hi = beta
        df = 0.0
        for zi, ki in zip(z, K):
            d = ki - 1.0
            w = 1.0 + beta * d
            df -= zi * d * d / (w * w)
        nxt = beta - f / df if df != 0.0 else 0.5 * (lo + hi)
        if not (lo < nxt < hi):
            nxt = 0.5 * (lo + hi)
        if abs(nxt - beta) < 1e-15:
            beta = nxt
            break
        beta = nxt
    return beta


def normalise(v):
    s = sum(v)
    return [x / s for x in v]


def pr_state(names, zs, T, p, pref):
    """PR state of one phase. pref = 'liquid' | 'vapour' | 'stable'."""
    n = len(names)
    MW = [CT[c][0] * 1.0e-3 for c in names]
    M = sum(z * m for z, m in zip(zs, MW))
    RT = R * T
    aI, bI = [], []
    for c in names:
        _, Tc, Pc, w = CT[c]
        kappa = 0.37464 + 1.54226 * w - 0.026992 * w * w
        alpha = (1.0 + kappa * (1.0 - math.sqrt(T / Tc))) ** 2
        Rtc = R * Tc
        aI.append(OMEGA_A * Rtc * Rtc / Pc * alpha)
        bI.append(OMEGA_B * Rtc / Pc)
    bMix = sum(zs[i] * bI[i] for i in range(n))
    aMix = sum(zs[i] * zs[j] * math.sqrt(aI[i] * aI[j])
               for i in range(n) for j in range(n))
    A = aMix * p / (RT * RT)
    B = bMix * p / RT
    roots = sorted(set(z for z in cubic_real_roots(
        -(1.0 - B), A - 3.0 * B * B - 2.0 * B, -(A * B - B * B - B ** 3))
        if z > B * (1.0 + 1e-12)))
    if not roots:
        raise SystemExit("PR cubic has no physical root")

    def build(Z):
        Vm = Z * RT / p
        logTerm = math.log((Z + (1.0 + SQRT2) * B) / (Z + (1.0 - SQRT2) * B))
        lnphi, gibbs = [], 0.0
        for i in range(n):
            sumA = sum(zs[j] * math.sqrt(aI[i] * aI[j]) for j in range(n))
            biB = bI[i] / bMix
            lphi = (biB * (Z - 1.0) - math.log(Z - B)
                    - A / (2.0 * SQRT2 * B) * (2.0 * sumA / aMix - biB)
                    * logTerm)
            lnphi.append(lphi)
            gibbs += zs[i] * lphi
        return {"Z": Z, "Vm": Vm, "rho": M / Vm, "overB": Vm / bMix,
                "lnphi": lnphi, "gibbs": gibbs, "M": M}

    if pref == "liquid":
        pick = min(roots)
    elif pref == "vapour":
        pick = max(roots)
    else:
        pick = min(roots, key=lambda Z: build(Z)["gibbs"])
    return build(pick)


def tp_flash(names, zs, T, p, tol=1.0e-12, maxIter=1000):
    """Successive-substitution TP flash, mirrored from src/pvt/src/flash.cpp."""
    n = len(names)
    K = [wilson_k(c, T, p) for c in names]
    for _ in range(1, maxIter + 1):
        f0, f1 = rr(zs, K, 0.0), rr(zs, K, 1.0)
        if f0 <= 0.0:
            beta = 0.0
            x = list(zs)
            y = normalise([K[i] * zs[i] for i in range(n)])
        elif f1 >= 0.0:
            beta = 1.0
            y = list(zs)
            x = normalise([zs[i] / K[i] for i in range(n)])
        else:
            beta = solve_rr(zs, K)
            x = normalise([zs[i] / (1.0 + beta * (K[i] - 1.0))
                           for i in range(n)])
            y = normalise([K[i] * x[i] for i in range(n)])
        sl = pr_state(names, x, T, p, "liquid")
        sv = pr_state(names, y, T, p, "vapour")
        newK = [math.exp(sl["lnphi"][i] - sv["lnphi"][i]) for i in range(n)]
        change = max(abs(math.log(nk / k)) for nk, k in zip(newK, K))
        K = newK
        if change < tol:
            break

    f0, f1 = rr(zs, K, 0.0), rr(zs, K, 1.0)
    spread = max(abs(k - 1.0) for k in K)
    two = f0 > 0.0 and f1 < 0.0 and spread > 1e-6
    if two:
        beta = solve_rr(zs, K)
        if 0.0 < beta < 1.0:
            x = normalise([zs[i] / (1.0 + beta * (K[i] - 1.0))
                           for i in range(n)])
            y = normalise([K[i] * x[i] for i in range(n)])
            distinct = max(abs(y[i] - x[i]) for i in range(n))
            two = distinct > 1e-8
    if not two:
        x = list(zs)
        y = list(zs)
        one = pr_state(names, zs, T, p, "stable")
        beta = 0.0 if one["overB"] < K_CRIT_OVER_B else 1.0
        liq = vap = one
    else:
        liq = pr_state(names, x, T, p, "liquid")
        vap = pr_state(names, y, T, p, "vapour")

    ML = sum(x[i] * CT[names[i]][0] * 1.0e-3 for i in range(n))
    MV = sum(y[i] * CT[names[i]][0] * 1.0e-3 for i in range(n))
    return {"single": not two, "beta": beta, "x": x, "y": y,
            "rhoL": liq["rho"], "rhoV": vap["rho"],
            "Vml": liq["Vm"], "Vmg": vap["Vm"], "ML": ML, "MV": MV}


def xi(Tc, Pc, MW):
    return Tc ** (1.0 / 6.0) / (math.sqrt(MW) * (Pc / ATM) ** (2.0 / 3.0))


def stiel_thodos(T, Tc, Pc, MW):
    Tr = T / Tc
    x = xi(Tc, Pc, MW)
    cp = (17.78e-5 * (4.58 * Tr - 1.67) ** 0.625 / x if Tr > 1.5
          else 34.0e-5 * Tr ** 0.94 / x)
    return cp * 1.0e-3


def lbc(T, Vm, zs, names):
    """Low-level LBC exactly as src/pvt/src/viscosity.cpp (Vm in m3/mol)."""
    MWs = [CT[c][0] for c in names]
    Tcs = [CT[c][1] for c in names]
    Pcs = [CT[c][2] for c in names]
    Vcs = [ZC * R * tc / pc for tc, pc in zip(Tcs, Pcs)]
    Tc = sum(z * t for z, t in zip(zs, Tcs))
    MW = sum(z * m for z, m in zip(zs, MWs))
    Pc = sum(z * p for z, p in zip(zs, Pcs))
    Vc = sum(z * v for z, v in zip(zs, Vcs))
    rr_ = Vc / Vm
    num = sum(z * stiel_thodos(T, t, p, m) * math.sqrt(m)
              for z, t, p, m in zip(zs, Tcs, Pcs, MWs))
    den = sum(z * math.sqrt(m) for z, m in zip(zs, MWs))
    mu0_cp = num / den * 1.0e3
    poly = (0.1023 + 0.023364 * rr_ + 0.058533 * rr_ ** 2
            - 0.040758 * rr_ ** 3 + 0.0093324 * rr_ ** 4)
    return (mu0_cp + (poly ** 4 - 1.0e-4) / xi(Tc, Pc, MW)) * 1.0e-3


def weinaug_katz(parachors, Vml, Vmg, xs, ys):
    tot = sum(P * (x / Vml - y / Vmg) for P, x, y in zip(parachors, xs, ys))
    return tot ** 4


def props(names, zs, T, p, parachors):
    """Independent per-station rho/mu/sigma from ONE flash."""
    r = tp_flash(names, zs, T, p)
    beta = r["beta"]
    mL = (1.0 - beta) * r["ML"]
    mV = beta * r["MV"]
    mt = mL + mV
    wL, wV = mL / mt, mV / mt
    rho = 1.0 / (wL / r["rhoL"] + wV / r["rhoV"])
    muL = lbc(T, r["Vml"], r["x"], names)
    muV = lbc(T, r["Vmg"], r["y"], names)
    mu = wL * muL + wV * muV
    sigma = 0.0 if r["single"] else weinaug_katz(
        parachors, r["Vml"], r["Vmg"], r["x"], r["y"])
    return {"rho": rho, "mu": mu, "sigma": sigma,
            "beta": beta, "single": r["single"],
            "rhoL": r["rhoL"], "rhoV": r["rhoV"]}


def churchill(re, eps):
    re = max(re, 1e-12)
    a1 = (7.0 / re) ** 0.9 + 0.27 * eps
    a = (2.457 * math.log(1.0 / a1)) ** 16
    b = (37530.0 / re) ** 16
    return 8.0 * ((8.0 / re) ** 12 + (a + b) ** (-1.5)) ** (1.0 / 12.0)


# ---- parse the case ----
with open(case) as fh:
    src = fh.read()
grid = [float(x) for x in re.search(r"pressure\s*\(([^)]*)\)", src).group(1).split()]
T = float(re.search(r"temperature\s+([-+0-9.eE]+)", src).group(1))
comp_block = re.search(r"composition\s*\{([^}]*)\}", src).group(1)
names = []
zs = []
for tok in re.finditer(r"(\w+)\s+([-+0-9.eE]+)\s*;", comp_block):
    names.append(tok.group(1))
    zs.append(float(tok.group(2)))
check("composition parsed", len(names) == 2, f"{names} {zs}")
check("parachors fromTable declared", "parachors fromTable" in src)
parachors = []
for c in names:
    check(f"{c} has a built-in parachor", c in PARA_CGS)
    parachors.append(PARA_CGS[c] * (10.0 ** -6.75))  # cgs -> SI


def run(path):
    return subprocess.run([exe, path], capture_output=True, text=True)


proc = run(case)
check("liveFlash case returns 0", proc.returncode == 0,
      f"rc={proc.returncode} {proc.stderr[:200]}")
out = proc.stdout

# ---- flash: per-station rho/mu/sigma vs independent re-derivation ----
stations = []
for m in re.finditer(r"flashStation\s+(\S+)\s+i\s+(\d+)\s+p\s+([-+0-9.eE]+)"
                     r"\s+rho\s+([-+0-9.eE]+)\s+mu\s+([-+0-9.eE]+)"
                     r"\s+sigma\s+([-+0-9.eE]+)", out):
    stations.append({"pipe": m.group(1), "i": int(m.group(2)),
                     "p": float(m.group(3)), "rho": float(m.group(4)),
                     "mu": float(m.group(5)), "sigma": float(m.group(6))})
check("24 flash stations printed", len(stations) == 24, f"{len(stations)}")
check("all stations belong to Line",
      all(s["pipe"] == "Line" for s in stations))

py_props = []
for s in stations:
    pr = props(names, zs, T, s["p"], parachors)
    py_props.append(pr)
    tag = f"station {s['i']} p={s['p']:.6g}"
    close(f"{tag} rho independent flash", s["rho"], pr["rho"], 1e-8)
    close(f"{tag} mu independent LBC", s["mu"], pr["mu"], 1e-8)
    close(f"{tag} sigma independent WK", s["sigma"], pr["sigma"],
          1e-8, 1e-15)
    check(f"{tag} two-phase on this branch", not pr["single"])

# ---- sigma: per-station and varying ----
sigmas = [s["sigma"] for s in stations]
pressures = [s["p"] for s in stations]
check("station pressures descend", all(b < a for a, b in zip(pressures, pressures[1:])))
check("sigma varies along the march",
      max(sigmas) > 1.01 * min(sigmas),
      f"min {min(sigmas):.6g} max {max(sigmas):.6g}")
check("WK sigma rises as pressure falls (more vapour)",
      all(b > a for a, b in zip(sigmas, sigmas[1:])),
      f"{sigmas[0]:.6g} -> {sigmas[-1]:.6g}")
m = re.search(r"pipe Line .*sigmaOutlet\s+([-+0-9.eE]+)", out)
check("pipe line reports sigmaOutlet", m is not None)
if m:
    close("sigmaOutlet == last station sigma", float(m.group(1)),
          sigmas[-1], 1e-12)

# ---- march: replay the energy balance with the independent properties ----
node_w = float(re.search(r"node W p\s+([-+0-9.eE]+)", out).group(1))
node_t = float(re.search(r"node T p\s+([-+0-9.eE]+)", out).group(1))
m = re.search(r"pipe Line Q\s+([-+0-9.eE]+)\s+v\s+([-+0-9.eE]+)"
              r"\s+Re\s+([-+0-9.eE]+)\s+f\s+([-+0-9.eE]+)"
              r"\s+dp\s+([-+0-9.eE]+)", out)
check("pipe Line parsed", m is not None)
if m and len(py_props) == 24:
    Q, v, Re, f, dp = (float(x) for x in m.groups())
    L, D, eps = 1000.0, 0.03, 1.0e-5
    area = math.pi * D * D / 4.0
    ds = L / 24
    close("station 0 pressure == upstream node", pressures[0], node_w, 1e-11)
    mass = py_props[0]["rho"] * Q
    for i in range(24):
        rho = py_props[i]["rho"]
        mu = py_props[i]["mu"]
        vv = mass / (rho * area)
        reyn = rho * vv * D / mu
        ff = churchill(reyn, eps / D)
        dpseg = (ff * ds / D) * rho * vv * vv / 2.0
        p_next = pressures[i] - dpseg
        if i < 23:
            close(f"segment {i} lands on next station p", p_next,
                  pressures[i + 1], 1e-8)
        else:
            close("last segment lands on downstream node", p_next,
                  node_t, 1e-8)
        if i == 23:
            close("replayed outlet velocity", vv, v, 1e-8)
            close("replayed outlet Reynolds", reyn, Re, 1e-6)
            close("replayed outlet f", ff, f, 1e-6)
    close("reported dp consistent with node pressures", dp,
          node_w - node_t, 1e-9)


# ---- convergence: live flash vs precomputed rho/mu tables on fine grids ----
def table_twin(n, path):
    """Explicit pvt case with rho/mu from the independent flash on n nodes."""
    lo, hi = pressures[-1] * 0.99, pressures[0] * 1.01
    pg = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
    rows = [props(names, zs, T, p, parachors) for p in pg]
    plist = " ".join(repr(p) for p in pg)
    rholist = " ".join(repr(r["rho"]) for r in rows)
    mulist = " ".join(repr(r["mu"]) for r in rows)
    text = src
    # Strip // comments, replace the whole compositional{...} with rho/mu.
    lines = [l.split("//")[0] for l in text.splitlines()]
    clean = "\n".join(lines)
    start = clean.index("compositional")
    brace = clean.index("{", start)
    depth = 0
    end = brace
    for j in range(brace, len(clean)):
        if clean[j] == "{":
            depth += 1
        elif clean[j] == "}":
            depth -= 1
            if depth == 0:
                end = j + 1
                break
    twin = (clean[:start]
            + f"rho ( {rholist} );\n        mu ( {mulist} );\n"
            + clean[end:])
    # Inject the new pressure grid (the source grid is only the envelope).
    twin = re.sub(r"pressure\s*\([^)]*\)", f"pressure ( {plist} );", twin,
                  count=1)
    assert "compositional" not in twin and "liveFlash" not in twin
    with open(path, "w") as fh:
        fh.write(twin)
    return run(path)


with tempfile.TemporaryDirectory(prefix="cf_flash_conv_") as d:
    results = []
    for n in (75, 600):
        rp = table_twin(n, os.path.join(d, f"twin{n}.dict"))
        check(f"fine-grid twin n={n} returns 0", rp.returncode == 0,
              f"rc={rp.returncode} {rp.stderr[:160]}")
        mm = re.search(r"node W p\s+([-+0-9.eE]+)", rp.stdout)
        if mm:
            results.append((n, float(mm.group(1))))
    check("two convergence twins ran", len(results) == 2)
    if len(results) == 2:
        coarse_err = abs(results[0][1] - node_w)
        fine_err = abs(results[1][1] - node_w)
        check("fine grid is closer to live flash than coarse grid",
              fine_err < coarse_err,
              f"coarse err {coarse_err:.6g} fine err {fine_err:.6g}")
        check("fine-grid twin converged to live-flash physics (rel 1e-3)",
              fine_err <= 1e-3 * abs(node_w),
              f"live W {node_w:.12g} fine W {results[1][1]:.12g}")
        print(f"  convergence: live W={node_w:.9g}, "
              f"coarse(75) W={results[0][1]:.9g}, "
              f"fine(600) W={results[1][1]:.9g}")


# ---- gate: absence / false declaration of the opt-in is byte-identical ----
def without(src, needles):
    text = src
    for nd in needles:
        text = text.replace(nd, "")
    return text


no_decl = without(src, ["liveFlash true;", "sigmaFromP true;",
                        "parachors fromTable;"])
false_decl = (src.replace("liveFlash true;", "liveFlash false;")
                 .replace("sigmaFromP true;", "sigmaFromP false;"))
with tempfile.TemporaryDirectory(prefix="cf_flash_gate_") as d:
    p1 = os.path.join(d, "none.dict")
    p2 = os.path.join(d, "false.dict")
    with open(p1, "w") as fh:
        fh.write(no_decl)
    with open(p2, "w") as fh:
        fh.write(false_decl)
    r1, r2 = run(p1), run(p2)
    check("no-declaration case returns 0", r1.returncode == 0,
          f"rc={r1.returncode} {r1.stderr[:160]}")
    check("false-declaration case returns 0", r2.returncode == 0,
          f"rc={r2.returncode} {r2.stderr[:160]}")
    check("no-declaration case produces no flashStation line",
          "flashStation" not in r1.stdout)
    check("false-declaration case produces no flashStation line",
          "flashStation" not in r2.stdout)
    check("abstain / false opt-in outputs byte-identical",
          r1.stdout == r2.stdout,
          f"none:\n{r1.stdout}\nfalse:\n{r2.stdout}")


# ---- refusals ----
def refuse(label, text, needle):
    with tempfile.TemporaryDirectory(prefix="cf_flash_refuse_") as d:
        p = os.path.join(d, "case.dict")
        with open(p, "w") as fh:
            fh.write(text)
        r = run(p)
        check(f"refuse {label}",
              r.returncode != 0 and needle in (r.stdout + r.stderr),
              f"rc={r.returncode} out={r.stdout[:120]} err={r.stderr[:160]}")


refuse("missing composition",
       src.replace("composition { methane 0.5; nButane 0.5; }", ""),
       "needs a composition block")
refuse("missing temperature",
       src.replace("temperature 320.0;", ""),
       "needs temperature")
refuse("sigmaFromP without parachors",
       src.replace("parachors fromTable;", ""),
       "sigmaFromP needs parachors")
refuse("parachors wrong size",
       src.replace("parachors fromTable;", "parachors ( 1.0e-5 );"),
       "one value per component")
refuse("rho with compositional",
       src.replace("pressure ( 2.0e6",
                   "rho ( 1 2 3 4 5 6 7 );\n        pressure ( 2.0e6"),
       "cannot declare rho/mu")

if fails:
    print(f"FAIL {len(fails)} checks:")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print(f"PASS {len(stations)} stations: independent PR flash/LBC/WK < 1e-8, "
      f"march replay, per-station sigma, fine-grid convergence, "
      f"zero-declaration gate, 5 refusals")
