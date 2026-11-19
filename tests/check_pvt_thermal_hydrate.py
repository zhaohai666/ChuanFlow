#!/usr/bin/env python3
"""D13 increment 4 verification: Prandtl thermal conductivity + hydrate
Clausius–Clapeyron equilibrium (pvt layer, declared-input closures).

usage: check_pvt_thermal_hydrate.py <exe> <case.dict>

Both reported quantities are restated with plain Python math, independent of
the C++ implementation and of any tabulated data:

  k(T)      = mu * cp / Pr
  p_eq(T)   = pRef * exp( -(dH/R) (1/T - 1/Tref) )

Constants (mu, cp, Pr, Tref, pRef, dH) are read from the case file, so the
checker never trusts a number the app printed. Also checks the monotonicity /
reference-point identities, six loud rejections, and a zero-declaration twin
whose flash report stays byte-identical when both new blocks are removed.
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
R = 8.31446261815324  # kGasConstant, must match cf/pvt/Components.h

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


tc = block("thermalConductivity")
mu, cp, Pr = scalar(tc, "mu"), scalar(tc, "cp"), scalar(tc, "Pr")
he = block("hydrateEquilibrium")
Tref, pRef, dH = scalar(he, "Tref"), scalar(he, "pRef"), scalar(he, "dH")
temps = [float(x) for x in re.search(r"temperatures\s*\(([^)]*)\)", he).group(1).split()]

proc = subprocess.run([exe, case], capture_output=True, text=True)
if proc.returncode != 0:
    raise SystemExit(f"pvt failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")
out = proc.stdout

# ---- thermal conductivity ----
k_lines = [ln for ln in out.splitlines() if "[chflow] pvt thermalConductivity" in ln]
assert len(k_lines) == 1, f"expected one thermalConductivity line, got {len(k_lines)}"
kv = dict(zip(k_lines[0].split()[3::2], k_lines[0].split()[4::2]))
k_got = float(kv["k"])
k_want = mu * cp / Pr
assert abs(k_got - k_want) <= 1e-6 * abs(k_want), (
    f"k={k_got} != mu*cp/Pr={k_want}")
print(f"thermal conductivity k={k_got:.10g} W/(m.K) matches mu*cp/Pr={k_want:.10g}")

# ---- hydrate equilibrium ----
h_lines = [ln for ln in out.splitlines() if "[chflow] pvt hydrateEquilibrium" in ln]
assert len(h_lines) == len(temps), (
    f"expected {len(temps)} hydrate lines, got {len(h_lines)}")
got_T, got_p = [], []
for ln in h_lines:
    parts = ln.replace("[chflow]", "").split()
    # pvt hydrateEquilibrium T <T> p <p>
    got_T.append(float(parts[3]))
    got_p.append(float(parts[5]))
assert [round(t, 9) for t in got_T] == [round(t, 9) for t in temps], (
    "hydrate temperatures not in case order")

for T, p in zip(got_T, got_p):
    ln_ratio = -(dH / R) * (1.0 / T - 1.0 / Tref)
    want = pRef * math.exp(ln_ratio)
    assert abs(p - want) <= 1e-6 * abs(want), (
        f"p_eq({T})={p} != {want}")
# Reference-point identity: at Tref the pressure is exactly pRef.
anchor = min(range(len(temps)), key=lambda i: abs(temps[i] - Tref))
assert abs(got_p[anchor] - pRef) <= 1e-9 * pRef, "p_eq(Tref) != pRef"
# Monotonic rise with T (positive dissociation enthalpy).
assert all(got_p[i] < got_p[i + 1] for i in range(len(got_p) - 1)), (
    "p_eq must increase with T")
print(f"hydrate equilibrium: {len(temps)} points match Clausius–Clapeyron; "
      f"p(Tref)={got_p[anchor]:.6g} == pRef; monotonic")

# ---- loud rejections ----
def run_case_text(text):
    fd, path = tempfile.mkstemp(suffix=".dict", prefix="pvt_th_")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        return subprocess.run([exe, path], capture_output=True, text=True)
    finally:
        os.unlink(path)


rejects = [
    ("negative mu", source.replace("mu  0.00089", "mu -0.00089")),
    ("zero cp", source.replace("cp  4180.0", "cp 0")),
    ("negative Pr", source.replace("Pr  6.14", "Pr -6.14")),
    ("negative Tref", source.replace("Tref 285.25", "Tref -285.25")),
    ("negative pRef", source.replace("pRef 3.0e6", "pRef -3.0e6")),
    ("zero dH", source.replace("dH   5.1868e4", "dH 0")),
]
for label, bad in rejects:
    r = run_case_text(bad)
    assert r.returncode != 0, f"{label} was not rejected"
print(f"rejections: {len(rejects)} invalid declarations all refused")

# ---- zero-declaration twin: flash lines byte-identical ----
twin = re.sub(r"thermalConductivity\s*\{.*?\}", "", source, flags=re.S)
twin = re.sub(r"hydrateEquilibrium\s*\{.*?\}", "", twin, flags=re.S)
r = run_case_text(twin)
assert r.returncode == 0, f"twin case failed: {r.stderr}"
twin_lines = [ln for ln in r.stdout.splitlines() if "[chflow] pvt flash" in ln]
base_lines = [ln for ln in out.splitlines() if "[chflow] pvt flash" in ln]
assert twin_lines == base_lines, "flash report changed when new blocks absent"
assert "[chflow] pvt thermalConductivity" not in r.stdout
assert "[chflow] pvt hydrateEquilibrium" not in r.stdout
print("zero-declaration twin: flash report byte-identical, no new tokens")

print("PASS: D13 increment 4 (thermal conductivity + hydrate equilibrium)")
