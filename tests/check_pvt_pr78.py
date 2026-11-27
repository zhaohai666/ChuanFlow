#!/usr/bin/env python3
"""5.230 (H3): independent verification of the PR1978 alpha branch.

The prAlpha{} echo is checked against the two published closed forms,
re-derived here in Python with the kernel's own operation order:

  omega <= 0.491 (PR76):  kappa = 0.37464 + 1.54226 w - 0.26992 w^2
  omega >   0.491 (PR78): kappa = 0.379642 + 1.48503 w - 0.164423 w^2
                                 + 0.016666 w^3
  alpha = (1 + kappa (1 - sqrt(T/Tc)))^2
  a_i   = 0.45724 (R Tc)^2 / Pc * alpha,  b_i = 0.07780 R Tc / Pc
  da_i/dT = 0.45724 (R Tc)^2 / Pc * (-kappa f / sqrt(T Tc)), f = 1+kappa(1-sqrtTr)

The PR76 arm (omega = 0.4) is compared BITWISE -- the operation order is
reproduced exactly, so equality is the "unchanged since before 5.230" witness.
The PR78 arm (omega = 0.5) is compared at 1e-14. dadT must carry the SAME
branch (the enthalpy departure's consistency criterion). The flash/enthalpy
twin (omega moved 0.5 -> 0.491) must stay continuous: the branch switch is a
0.4% kappa step, so every reported output moves by far less than that.

usage: check_pvt_pr78.py <pvtExe> <caseDict>
"""
import math
import os
import re
import subprocess
import sys

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

PVT_EXE, CASE = sys.argv[1:3]
CASE = os.path.abspath(CASE)

R_GAS = 8.314462618
OMEGA_A = 0.45723552892138218   # kOmegaA, PR(1976) exact
OMEGA_B = 0.07779607390388846   # kOmegaB, PR(1976) exact

fails = []
NCHECK = [0]


def check(name, cond, detail=""):
    NCHECK[0] += 1
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def near(got, want, rtol=1e-14, atol=0.0):
    return abs(got - want) <= rtol * abs(want) + atol


def run(case):
    r = subprocess.run([PVT_EXE, case], capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def pr_kappa(w):
    if w <= 0.491:
        return 0.37464 + 1.54226 * w - 0.26992 * w * w
    return 0.379642 + 1.48503 * w - 0.164423 * w * w + 0.016666 * w * w * w


src = open(CASE).read()

# the declared custom components
comps = {}
for m in re.finditer(r"(\w+)\s*\{\s*fraction ([-+0-9.eE]+);\s*Tc ([-+0-9.eE]+);"
                     r"\s*Pc ([-+0-9.eE]+);\s*omega ([-+0-9.eE]+);"
                     r"\s*molarMass ([-+0-9.eE]+);", src):
    comps[m.group(1)] = dict(fraction=float(m.group(2)), Tc=float(m.group(3)),
                             Pc=float(m.group(4)), omega=float(m.group(5)),
                             M=float(m.group(6)))
T = float(re.search(r"prAlpha\s*\{\s*temperature ([-+0-9.eE]+)", src).group(1))

rc, out = run(CASE)
check("app ran", rc == 0)

echo = {}
for m in re.finditer(r"prAlpha comp (\w+) T ([-+0-9.eE]+) omega ([-+0-9.eE]+)"
                     r" kappa ([-+0-9.eE]+) alpha ([-+0-9.eE]+) a ([-+0-9.eE]+)"
                     r" b ([-+0-9.eE]+) dadT ([-+0-9.eE]+)", out):
    echo[m.group(1)] = dict(
        T=float(m.group(2)), omega=float(m.group(3)), kappa=float(m.group(4)),
        alpha=float(m.group(5)), a=float(m.group(6)), b=float(m.group(7)),
        dadT=float(m.group(8)))

check("prAlpha echoed every component",
      set(echo) == set(comps), f"{sorted(echo)} vs {sorted(comps)}")

for name, c in comps.items():
    e = echo.get(name)
    if e is None:
        continue
    k = pr_kappa(c["omega"])
    sqrtTr = math.sqrt(T / c["Tc"])
    f = 1.0 + k * (1.0 - sqrtTr)
    alpha = f * f
    Rtc = R_GAS * c["Tc"]
    a = OMEGA_A * Rtc * Rtc / c["Pc"] * alpha
    b = OMEGA_B * Rtc / c["Pc"]
    dadT = OMEGA_A * Rtc * Rtc / c["Pc"] * (-k * f / math.sqrt(T * c["Tc"]))
    check(f"{name} omega echoed", e["omega"] == c["omega"])
    if c["omega"] <= 0.491:
        # PR76 arm: bit-identical to the pre-5.230 kernel (same expression,
        # same operation order) -- equality, not tolerance.
        check(f"{name} PR76 kappa bitwise", e["kappa"] == k,
              f"{e['kappa']!r} vs {k!r}")
        check(f"{name} PR76 alpha bitwise", e["alpha"] == alpha)
        check(f"{name} PR76 a bitwise", e["a"] == a)
        check(f"{name} PR76 dadT bitwise", e["dadT"] == dadT)
    else:
        check(f"{name} PR78 kappa (1e-14)", near(e["kappa"], k),
              f"{e['kappa']!r} vs {k!r}")
        check(f"{name} PR78 alpha (1e-14)", near(e["alpha"], alpha))
        check(f"{name} PR78 a (1e-14)", near(e["a"], a))
        check(f"{name} PR78 b (1e-14)", near(e["b"], b))
        check(f"{name} PR78 dadT (1e-14) -- same branch", near(e["dadT"], dadT),
              f"{e['dadT']!r} vs {dadT!r}")

# the PR78 kappa really differs from PR76 at omega=0.5 (the branch is live)
k76 = 0.37464 + 1.54226 * 0.5 - 0.26992 * 0.5 * 0.5
k78 = pr_kappa(0.5)
check("PR78 kappa(0.5) != PR76 kappa(0.5)", k76 != k78,
      f"{k76!r} vs {k78!r}")

# --- flash / enthalpy twin continuity: omega 0.5 -> 0.491 ---
def outputs(text):
    vals = {}
    for m in re.finditer(r"flash T ([-+0-9.eE]+) p ([-+0-9.eE]+) beta ([-+0-9.eE]+)",
                         text):
        vals["beta@" + m.group(2)] = float(m.group(3))
    for m in re.finditer(r"liquid rho ([-+0-9.eE]+)", text):
        vals.setdefault("rhoL", []).append(float(m.group(1)))
    for m in re.finditer(r"enthalpyAt T ([-+0-9.eE]+) p [-+0-9.eE]+ h ([-+0-9.eE]+)",
                         text):
        vals["h@" + m.group(1)] = float(m.group(2))
    return vals

base = outputs(out)
twin = src.replace("omega 0.5;", "omega 0.491;")
with open("/tmp/_cf_pr78_twin.dict", "w") as fh:
    fh.write(twin)
rc2, out2 = run("/tmp/_cf_pr78_twin.dict")
check("twin (omega 0.491) ran", rc2 == 0)
t = outputs(out2)
cont = True
for k in base:
    if isinstance(base[k], list) or k not in t:
        continue
    a, b = base[k], t[k]
    if abs(a - b) > 0.05 * max(abs(a), abs(b)):
        cont = False
        print(f"   discontinuity at {k}: {a!r} vs {b!r}")
check("flash/enthalpy outputs continuous across the branch", cont)

# --- custom-component loader rejections ---
bad = src.replace("Tc 500.0;", "Tc -1.0;")
with open("/tmp/_cf_pr78_badTc.dict", "w") as fh:
    fh.write(bad)
rc3, err3 = run("/tmp/_cf_pr78_badTc.dict")
check("reject non-positive Tc", rc3 != 0 and "Tc" in err3, err3[-120:])

bad = src.replace("fraction 0.3;", "fraction 0.0;")
with open("/tmp/_cf_pr78_badf.dict", "w") as fh:
    fh.write(bad)
rc4, err4 = run("/tmp/_cf_pr78_badf.dict")
check("reject non-positive fraction", rc4 != 0 and "fraction" in err4,
      err4[-120:])

print(f"\n{NCHECK[0] - len(fails)}/{NCHECK[0]} checks passed")
if fails:
    raise SystemExit("FAILED: " + ", ".join(fails))
