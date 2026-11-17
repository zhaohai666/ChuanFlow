#!/usr/bin/env python3
"""Analytic verification of the LBC compositional viscosity (chflow-pvt).

usage: check_pvt_lbc.py <exe> <case.dict>

D7 increment 1 (5.93). The Lohrenz-Bray-Clark (1964, SPE 915) residual
viscosity is checked two independent ways:

  formula    a pure-Python transcription of the LBC equations reproduces the
             Whitson & Brule (2000) worked example as published by the
             `chemicals` library to 1e-12 -- this pins the coefficients
             (a4 = 0.0093724), the xi units, the Stiel-Thodos dilute-gas term,
             the Herning-Zipperer mixing and rr = Vc/Vm.
  kernel     for every flashed phase the checker re-derives mu from the
             reported phase state (Vm = M/rho, PR critical volume
             Vc = 0.3074 R Tc/Pc, Kay-averaged constants) and weighs it against
             the C++ report.
  physics    mu > 0; the dilute-gas limit at the lowest pressure; liquid more
             viscous than vapour on the two-phase branch; the Herning-Zipperer
             mu0 bounded by the pure-component viscosities.
  gate       removing the `viscosity` block drops the report and leaves every
             flash line bit-identical (untouched path).
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
ATM = 101325.0
ZC = 0.3074  # PR76 universal critical compressibility

# Exact critical constants from src/pvt/src/components.cpp (MW g/mol, Tc K, Pc Pa).
CT = {
    "methane": (16.043, 190.564, 45.992e5),
    "ethane": (30.070, 305.322, 48.722e5),
    "propane": (44.097, 369.834, 42.477e5),
    "iButane": (58.123, 407.817, 36.480e5),
    "nButane": (58.123, 425.125, 37.962e5),
    "nPentane": (72.150, 469.700, 33.700e5),
    "nHexane": (86.178, 507.600, 30.250e5),
    "carbonDioxide": (44.010, 304.128, 73.773e5),
    "hydrogenSulphide": (34.081, 373.400, 89.628e5),
    "nitrogen": (28.014, 126.192, 33.958e5),
}


def xi(Tc, Pc, MW):
    return Tc ** (1.0 / 6.0) / (math.sqrt(MW) * (Pc / ATM) ** (2.0 / 3.0))


def stiel_thodos(T, Tc, Pc, MW):
    Tr = T / Tc
    x = xi(Tc, Pc, MW)
    cp = (17.78e-5 * (4.58 * Tr - 1.67) ** 0.625 / x if Tr > 1.5
          else 34.0e-5 * Tr ** 0.94 / x)
    return cp * 1.0e-3  # Pa s


def lbc(T, Vm, zs, MWs, Tcs, Pcs, Vcs):
    Tc = sum(z * t for z, t in zip(zs, Tcs))
    MW = sum(z * m for z, m in zip(zs, MWs))
    Pc = sum(z * p for z, p in zip(zs, Pcs))
    Vc = sum(z * v for z, v in zip(zs, Vcs))
    rr = Vc / Vm
    num = sum(z * stiel_thodos(T, t, p, m) * math.sqrt(m)
              for z, t, p, m in zip(zs, Tcs, Pcs, MWs))
    den = sum(z * math.sqrt(m) for z, m in zip(zs, MWs))
    mu0_cp = num / den * 1.0e3
    poly = (0.1023 + 0.023364 * rr + 0.058533 * rr ** 2
            - 0.040758 * rr ** 3 + 0.0093724 * rr ** 4)
    return (mu0_cp + (poly ** 4 - 1.0e-4) / xi(Tc, Pc, MW)) * 1.0e-3


def mu0_mix(T, zs, MWs, Tcs, Pcs):
    num = sum(z * stiel_thodos(T, t, p, m) * math.sqrt(m)
              for z, t, p, m in zip(zs, Tcs, Pcs, MWs))
    den = sum(z * math.sqrt(m) for z, m in zip(zs, MWs))
    return num / den


fails = []


def check(name, ok, detail=""):
    if not ok:
        fails.append(f"{name} {detail}".strip())


def close(name, got, want, rtol, atol=0.0):
    ok = abs(got - want) <= atol + rtol * abs(want)
    check(name, ok, f"got {got:.12g} want {want:.12g}")


# ---- formula self-test against the published worked example ----
example = lbc(300.0, 0.0023025, [.4, .3, .3],
              [16.04246, 30.06904, 44.09562],
              [190.564, 305.32, 369.83],
              [4599000.0, 4872000.0, 4248000.0],
              [9.86e-05, 0.0001455, 0.0002])
close("LBC reproduces the Whitson-Brule published example",
      example, 9.925488160761484e-06, 0.0, 1.0e-12)


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
                   "x": {}, "y": {}}
            flashes.append(cur)
        elif k == "liquid" and cur is not None:
            cur["rhoL"], cur["ML"] = float(t[3]), float(t[5])
        elif k == "vapour" and cur is not None:
            cur["rhoV"], cur["MV"] = float(t[3]), float(t[5])
        elif k in ("x", "y") and cur is not None:
            cur[k][t[2]] = float(t[3])
        elif k == "viscosity" and cur is not None:
            cur["muL"], cur["muV"] = float(t[3]), float(t[5])
    return flashes


def phase_mu(T, comp, M, rho):
    names = list(comp.keys())
    zs = [comp[n] for n in names]
    MWs = [CT[n][0] for n in names]
    Tcs = [CT[n][1] for n in names]
    Pcs = [CT[n][2] for n in names]
    Vcs = [ZC * R * CT[n][1] / CT[n][2] for n in names]
    return lbc(T, M / rho, zs, MWs, Tcs, Pcs, Vcs), (T, zs, MWs, Tcs, Pcs)


out = run(case)
flashes = parse(out)
check("case produced flashes", len(flashes) >= 3, f"{len(flashes)} flashes")
check("case reports viscosity", all("muL" in f for f in flashes),
      "no pvt viscosity lines (is the viscosity block present?)")

for f in flashes:
    label = f"p={f['p']:.6g}"
    if "muL" not in f:
        continue
    mL, argsL = phase_mu(f["T"], f["x"], f["ML"], f["rhoL"])
    mV, argsV = phase_mu(f["T"], f["y"], f["MV"], f["rhoV"])
    close(f"{label} liquid LBC vs independent reference", f["muL"], mL, 1.0e-9, 1.0e-15)
    close(f"{label} vapour LBC vs independent reference", f["muV"], mV, 1.0e-9, 1.0e-15)
    check(f"{label} liquid viscosity positive", f["muL"] > 0.0, f"{f['muL']}")
    check(f"{label} vapour viscosity positive", f["muV"] > 0.0, f"{f['muV']}")
    # Herning-Zipperer mu0 bounded by the pure-component dilute-gas viscosities
    T, zs, MWs, Tcs, Pcs = argsV
    pure = [stiel_thodos(T, t, p, m) for t, p, m in zip(Tcs, Pcs, MWs)]
    mu0 = mu0_mix(T, zs, MWs, Tcs, Pcs)
    check(f"{label} vapour mu0 within pure-component bounds",
          min(pure) - 1.0e-15 <= mu0 <= max(pure) + 1.0e-15,
          f"mu0 {mu0:.6g} bounds [{min(pure):.6g}, {max(pure):.6g}]")
    if not f["single"]:
        check(f"{label} liquid more viscous than vapour", f["muL"] > f["muV"],
              f"{f['muL']:.6g} vs {f['muV']:.6g}")

# ---- dilute-gas limit: the lowest-pressure vapour is nearly mu0 ----
low = min(flashes, key=lambda f: f["p"])
if "muV" in low and low["rhoV"] > 0.0:
    _, args = phase_mu(low["T"], low["y"], low["MV"], low["rhoV"])
    T, zs, MWs, Tcs, Pcs = args
    mu0 = mu0_mix(T, zs, MWs, Tcs, Pcs)
    # at low density the residual (poly^4 - 1e-4)/xi is a small correction
    check("low-pressure vapour near the dilute-gas limit",
          abs(low["muV"] - mu0) < 0.10 * mu0,
          f"mu {low['muV']:.6g} mu0 {mu0:.6g}")

# ---- gate: removing the viscosity block drops the report, flashes unchanged ----
with open(case) as fh:
    src = fh.read()
twin = re.sub(r"viscosity\s*\{[^}]*\}", "", src, count=1)
twin_path = case + ".novisc"
with open(twin_path, "w") as fh:
    fh.write(twin)
twin_out = run(twin_path)
os.remove(twin_path)
check("untouched case reports no viscosity", "pvt viscosity" not in twin_out,
      "viscosity leaked into the gated-off report")
base_flash = [l for l in out.splitlines() if "pvt flash" in l]
twin_flash = [l for l in twin_out.splitlines() if "pvt flash" in l]
check("gated-off flash report bit-identical", base_flash == twin_flash,
      "flash lines changed when viscosity was removed")

if fails:
    print("\n".join(f"FAIL: {m}" for m in fails))
    raise SystemExit(f"{len(fails)} pvt-lbc check(s) failed")

print(f"pvt-lbc: {len(flashes)} flashes, published example matched, "
      f"kernel == independent LBC reference, gated report clean - all checks passed")