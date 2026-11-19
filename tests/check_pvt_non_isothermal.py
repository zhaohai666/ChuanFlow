#!/usr/bin/env python3
"""D13 increment 5 verification: non-isothermal (PH) flash with declared cp(T).

usage: check_pvt_non_isothermal.py <exe> <case.dict>

Nothing is trusted from the app's property functions. The checker independently
transcribes PR76 (cubic roots, fugacity, successive-substitution TP flash) and
re-derives the molar enthalpy at each sample T:

  h_ig_i(T) - h_ig(Tref) = cpA_i (T-Tref) + cpB_i/2 (T^2-Tref^2)
  da_i/dT               = -OmegaA (R Tc)^2/Pc * kappa f / sqrt(T Tc)
  H^R_phase             = R T (Z-1) + (T daMix/dT - aMix)/(2 sqrt2 bMix)
                          * ln[(Z+(1+sqrt2)B)/(Z+(1-sqrt2)B)]
  H_mix                 = (1-beta)(sum x_i(h_ig,i+H^R_L)) + beta(sum y_i(h_ig,i+H^R_V))

cpA/cpB/Tref/hTarget are read from the case. It checks every enthalpyAt sample,
that phFlash solves back to the T whose H equals hTarget, monotonicity, a
zero-declaration byte-identical twin, and loud rejections. Pure stdlib.
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
CT = {  # name -> (MW g/mol, Tc K, Pc Pa, acentric), from components.cpp
    "methane": (16.043, 190.564, 45.992e5, 0.0115),
    "ethane": (30.070, 305.322, 48.722e5, 0.0995),
    "propane": (44.097, 369.834, 42.477e5, 0.1523),
    "nButane": (58.123, 425.125, 37.962e5, 0.2002),
}

with open(case) as fh:
    source = fh.read()


def block(name):
    m = re.search(name + r"\s*\{(.*?)\}", source, re.S)
    if m is None:
        raise SystemExit(f"case has no {name} block")
    return m.group(1)


def scalar(text, key):
    m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)\s*;", text)
    if m is None:
        raise SystemExit(f"missing scalar {key}")
    return float(m.group(1))


def numbers(text, key):
    m = re.search(r"\b" + key + r"\s*\(([^)]*)\)", text)
    if m is None:
        raise SystemExit(f"missing list {key}")
    return [float(x) for x in m.group(1).split()]


comp = block("composition")
feed = {k: float(v) for k, v in re.findall(
    r"([A-Za-z_][A-Za-z0-9_]*)\s+([-+0-9.eE]+)\s*;", comp)}
names = list(feed)
zs = [feed[k] for k in names]
s = sum(zs)
zs = [z / s for z in zs]

ni = block("nonIsothermal")
cpA, cpB = numbers(ni, "cpA"), numbers(ni, "cpB")
Tref, pIso = scalar(ni, "referenceT"), scalar(ni, "pressure")
sampleT = numbers(ni, "enthalpyT")
hTarget = scalar(ni, "hTarget")


# ---- independent PR transcription ----
def cbrt(x):
    return math.copysign(abs(x) ** (1.0 / 3.0), x)


def cubic_real_roots(c2, c1, c0):
    p = c1 - c2 * c2 / 3.0
    q = 2.0 * c2 ** 3 / 27.0 - c2 * c1 / 3.0 + c0
    disc = q * q / 4.0 + p ** 3 / 27.0
    shift = -c2 / 3.0
    if disc > 0.0:
        sq = math.sqrt(disc)
        return [cbrt(-q / 2.0 + sq) + cbrt(-q / 2.0 - sq) + shift]
    if p >= 0.0:
        return [shift]
    r = math.sqrt(-p ** 3 / 27.0)
    m = 2.0 * math.sqrt(-p / 3.0)
    theta = math.acos(max(-1.0, min(1.0, -q / (2.0 * r))))
    return sorted(m * math.cos((theta + 2.0 * k * math.pi) / 3.0) + shift
                  for k in range(3))


def wilson_k(name, T, p):
    _, Tc, Pc, w = CT[name]
    return Pc / p * math.exp(5.373 * (1.0 + w) * (1.0 - Tc / T))


def rr(z, K, beta):
    return sum(zi * (ki - 1.0) / (1.0 + beta * (ki - 1.0))
               for zi, ki in zip(z, K))


def solve_rr(z, K):
    lo, hi, beta = 0.0, 1.0, 0.5
    for _ in range(200):
        f = rr(z, K, beta)
        if f > 0.0:
            lo = beta
        else:
            hi = beta
        df = -sum(zi * (ki - 1.0) ** 2 / (1.0 + beta * (ki - 1.0)) ** 2
                  for zi, ki in zip(z, K))
        nxt = beta - f / df if df != 0.0 else 0.5 * (lo + hi)
        if not (lo < nxt < hi):
            nxt = 0.5 * (lo + hi)
        if abs(nxt - beta) < 1e-15:
            return nxt
        beta = nxt
    return beta


def normalise(v):
    t = sum(v)
    return [x / t for x in v]


def pr_params(z, T, p):
    n = len(names)
    aI, bI = [], []
    for name in names:
        _, Tc, Pc, w = CT[name]
        kappa = 0.37464 + 1.54226 * w - 0.026992 * w * w
        alpha = (1.0 + kappa * (1.0 - math.sqrt(T / Tc))) ** 2
        Rtc = R * Tc
        aI.append(OMEGA_A * Rtc * Rtc / Pc * alpha)
        bI.append(OMEGA_B * Rtc / Pc)
    aMix = sum(z[i] * z[j] * math.sqrt(aI[i] * aI[j])
               for i in range(n) for j in range(n))
    bMix = sum(z[i] * bI[i] for i in range(n))
    A = aMix * p / (R * T) ** 2
    B = bMix * p / (R * T)
    return aI, bI, aMix, bMix, A, B


def pr_root(z, T, p, pref):
    aI, bI, aMix, bMix, A, B = pr_params(z, T, p)
    roots = sorted(set(x for x in cubic_real_roots(
        -(1.0 - B), A - 3.0 * B * B - 2.0 * B, -(A * B - B * B - B ** 3))
        if x > B * (1.0 + 1e-12)))
    if not roots:
        raise SystemExit("PR cubic has no physical root")
    Z = min(roots) if pref == "L" else max(roots)
    return aI, aMix, bMix, B, Z


def da_dt(name, T):
    _, Tc, Pc, w = CT[name]
    kappa = 0.37464 + 1.54226 * w - 0.026992 * w * w
    f = 1.0 + kappa * (1.0 - math.sqrt(T / Tc))
    return -OMEGA_A * (R * Tc) ** 2 / Pc * kappa * f / math.sqrt(T * Tc)


def daMix_dt(z, aI, T):
    n = len(names)
    return sum(z[i] * z[j] * 0.5 * (
        da_dt(names[i], T) * math.sqrt(aI[j] / aI[i])
        + da_dt(names[j], T) * math.sqrt(aI[i] / aI[j]))
        for i in range(n) for j in range(n))


def h_residual(z, T, p, pref):
    aI, aMix, bMix, B, Z = pr_root(z, T, p, pref)
    dadT = daMix_dt(z, aI, T)
    logTerm = math.log((Z + (1.0 + SQRT2) * B) / (Z + (1.0 - SQRT2) * B))
    return R * T * (Z - 1.0) + (T * dadT - aMix) / (2.0 * SQRT2 * bMix) * logTerm


def h_ideal(i, T):
    return cpA[i] * (T - Tref) + 0.5 * cpB[i] * (T * T - Tref * Tref)


def h_mix(T, p):
    n = len(names)
    K = [wilson_k(names[i], T, p) for i in range(n)]
    f0, f1 = rr(zs, K, 0.0), rr(zs, K, 1.0)
    if f0 <= 0.0:
        beta, x, y = 0.0, list(zs), normalise([K[i] * zs[i] for i in range(n)])
    elif f1 >= 0.0:
        beta, x, y = 1.0, normalise([zs[i] / K[i] for i in range(n)]), list(zs)
    else:
        beta = solve_rr(zs, K)
        x = normalise([zs[i] / (1.0 + beta * (K[i] - 1.0)) for i in range(n)])
        y = normalise([K[i] * x[i] for i in range(n)])
    hL = sum(x[i] * (h_ideal(i, T) + h_residual(x, T, p, "L")) for i in range(n))
    hV = sum(y[i] * (h_ideal(i, T) + h_residual(y, T, p, "V")) for i in range(n))
    return (1.0 - beta) * hL + beta * hV, beta


# ---- parse app output ----
proc = subprocess.run([exe, case], capture_output=True, text=True)
if proc.returncode != 0:
    raise SystemExit(f"pvt failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")
out = proc.stdout
h_lines = [ln for ln in out.splitlines() if "[chflow] pvt enthalpyAt" in ln]
assert len(h_lines) == len(sampleT), (
    f"expected {len(sampleT)} enthalpyAt lines, got {len(h_lines)}")

fails = []


def close(name, got, want, rtol=1e-6, atol=1e-6):
    if abs(got - want) > atol + rtol * abs(want):
        fails.append(f"{name}: got {got:.10g} want {want:.10g}")


for ln, Tcase in zip(h_lines, sampleT):
    parts = ln.replace("[chflow]", "").split()
    # pvt enthalpyAt T <T> p <p> h <h>
    Tgot, pgot, hgot = float(parts[3]), float(parts[5]), float(parts[7])
    assert abs(Tgot - Tcase) < 1e-9 and abs(pgot - pIso) < 1e-3
    hWant, beta = h_mix(Tcase, pIso)
    close(f"H(T={Tcase})", hgot, hWant, rtol=1e-6, atol=1e-4)
print(f"enthalpyAt: {len(sampleT)} points match independent PR departure enthalpy")

hs = [float(ln.replace("[chflow]", "").split()[7]) for ln in h_lines]
assert all(hs[i] < hs[i + 1] for i in range(len(hs) - 1)), "H must rise with T"
print("enthalpy: monotone increasing with T")

ph_lines = [ln for ln in out.splitlines() if "[chflow] pvt phFlash" in ln]
assert len(ph_lines) == 1, "expected one phFlash line"
p = ph_lines[0].replace("[chflow]", "").split()
kv = dict(zip(p[2::2], p[3::2]))
Tsol, hsol, res, conv = (float(kv["T"]), float(kv["h"]),
                         float(kv["residual"]), int(kv["converged"]))
assert conv == 1, "phFlash did not converge"
assert res < 1e-4, f"phFlash residual too large: {res}"
close("phFlash h == hTarget", hsol, hTarget, rtol=1e-6, atol=1e-3)
hAtTsol, _ = h_mix(Tsol, pIso)
close("independent H(Tsol) == hTarget", hAtTsol, hTarget, rtol=1e-5, atol=1.0)
# The solved T must reproduce the same h the independent model gives.
assert Tsol > 0
print(f"phFlash: solved T={Tsol:.6f} K, H={hsol:.6g}, residual={res:.3g}, "
      f"independent H(T)={hAtTsol:.6g}")


# ---- rejections ----
def run_text(text):
    fd, path = tempfile.mkstemp(suffix=".dict", prefix="pvt_ni_")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        return subprocess.run([exe, path], capture_output=True, text=True)
    finally:
        os.unlink(path)


rejects = [
    ("cpA wrong length", source.replace("cpA (20 40);", "cpA (20);")),
    ("cpB wrong length", source.replace("cpB (30 50);", "cpB (30 50 70);")),
    ("Tmin >= Tmax", source.replace("Tmin 250;", "Tmin 500;")),
    ("hTarget out of range", source.replace("hTarget 242384.567025;", "hTarget 9e9;")),
    ("negative pressure", source.replace("pressure 5.0e6;", "pressure -5.0e6;")),
]
for label, bad in rejects:
    r = run_text(bad)
    assert r.returncode != 0, f"{label} was not rejected"
print(f"rejections: {len(rejects)} invalid declarations all refused")

# ---- zero-declaration twin ----
twin = re.sub(r"nonIsothermal\s*\{.*?\}", "", source, flags=re.S)
r = run_text(twin)
assert r.returncode == 0, f"twin failed: {r.stderr}"
base_flash = [ln for ln in out.splitlines() if "[chflow] pvt flash" in ln]
twin_flash = [ln for ln in r.stdout.splitlines() if "[chflow] pvt flash" in ln]
assert base_flash == twin_flash, "flash report changed when nonIsothermal absent"
assert "enthalpyAt" not in r.stdout and "phFlash" not in r.stdout
print("zero-declaration twin: flash report byte-identical, no new tokens")

if fails:
    raise SystemExit("FAIL\n" + "\n".join(fails))
print("PASS: D13 increment 5 (non-isothermal PH flash)")
