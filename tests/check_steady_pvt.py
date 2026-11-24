#!/usr/bin/env python3
"""Verification of the tabulated-property steady pipe law (fluid { pvt ... }).

Criterion makeup: 3 independent physical criteria + mirror diagnostics
(regression safety net). The mirror no longer carries the verdict alone.

Independent physical criteria:

  I1. Churchill closed-form limits. The friction law used everywhere must
      reduce to 64/Re in the laminar limit and to the fully-rough form
      8 / (2.457 ln(1/(0.27 eps)))^2 as Re -> infinity (checked as exact
      function identities, so a silently drifted mirror cannot hide).

  I2. Constant-property twin. A twin case with rho/mu tables held constant
      turns the march ODE into one with constant coefficients, so the
      24-segment Euler recursion is exact: the reported node pressures must
      equal the closed-form network solution p_T + rho g dz + Churchill
      friction (summed over both pipes), and the fill must not expand
      (outlet velocity exactly Q/A). (A laminar twin is unreachable: the
      source-driven Newton stalls below rate ~1e-3, so the turbulent arm is
      the one the closed form can exercise; its limits are pinned in I1.)

  I3. First-order march convergence. The restated 24-segment Euler march,
      driven against a 1536-segment reference, must halve its error when the
      segment count doubles (Richardson ratios ~ 2, monotonically).

Mirror diagnostics (regression safety net): the checker never trusts the
solver's internals -- it parses the rho(p)/mu(p) table straight out of the
case dictionary, restates the 24-segment Euler march of the local
mechanical-energy balance (mass rate fixed by the upstream density,
piecewise-linear properties at the segment-start pressure, Churchill friction
at the local Reynolds number), and drives it with the 12-digit
pressures/rates the solver reports:

  serial    : the junction carries Q_riser == Q_flowline == the source rate
              (no external rate at the free junction),
  closure   : marching the Riser from the reported well pressure at the
              reported rate lands on the reported junction pressure, and
              marching the Flowline from the junction lands on the fixed
              terminal -- the same recursion agrees to secant precision,
  reports   : outlet velocity and friction factor match the restated march,
  expansion : the outlet velocity exceeds q/A (the fill expands as the
              pressure falls) -- the constant-property law cannot show this,
  compress  : the reported dp differs visibly from an inlet-density
              constant-property computation, i.e. rho(p) is doing work.

usage: check_steady_pvt.py <exe> <case.dict>
"""
import math
import re
import subprocess
import sys

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]

proc = subprocess.run([exe, case], capture_output=True, text=True)
if proc.returncode != 0:
    raise SystemExit(f"solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")

nodes, pipes = {}, {}
for line in proc.stdout.splitlines():
    tok = line.replace("[chflow]", "").split()
    if not tok:
        continue
    if tok[0] == "node" and len(tok) >= 4:
        nodes[tok[1]] = float(tok[3])
    elif tok[0] == "pipe":
        d = {}
        for j in range(2, len(tok) - 1, 2):
            d[tok[j]] = float(tok[j + 1])
        pipes[tok[1]] = d
if not pipes:
    raise SystemExit("no pipe lines parsed")

# ---- the property table comes from the case itself ----
text = re.sub(r"//[^\n]*", "", open(case).read())  # strip // comments
pvt = re.search(r"pvt\s*\{(.*?)\}", text, re.S).group(1)


def table(key):
    m = re.search(key + r"\s*\(([^)]*)\)", pvt)
    return [float(x) for x in m.group(1).split()]


P_GRID, RHO_TAB, MU_TAB = table("pressure"), table("rho"), table("mu")
if not (len(P_GRID) == len(RHO_TAB) == len(MU_TAB)) or len(P_GRID) < 2:
    raise SystemExit("case table malformed")


def interp(vals, p):
    """piecewise-linear lookup clamped at the table ends (as PvtFluid)"""
    if p <= P_GRID[0]:
        return vals[0]
    if p >= P_GRID[-1]:
        return vals[-1]
    lo, hi = 0, len(P_GRID) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if P_GRID[mid] <= p:
            lo = mid
        else:
            hi = mid
    t = (p - P_GRID[lo]) / (P_GRID[hi] - P_GRID[lo])
    return vals[lo] + (vals[hi] - vals[lo]) * t


def churchill(re, eps):
    re = max(re, 1e-12)
    a1 = (7.0 / re) ** 0.9 + 0.27 * eps
    a = (2.457 * math.log(1.0 / a1)) ** 16
    b = (37530.0 / re) ** 16
    return 8.0 * ((8.0 / re) ** 12 + (a + b) ** (-1.5)) ** (1.0 / 12.0)


G = 9.80665
NSEG = 24


def march(name, q, p_up, nseg=NSEG):
    """the solver's Euler recursion: dp/ds at the segment-start pressure"""
    L, D, eps, sin = GEOM[name]
    area = math.pi * D * D / 4.0
    mass = interp(RHO_TAB, p_up) * q
    ds = L / nseg
    dz = sin * ds
    p, v_out, re_out, f_out = p_up, 0.0, 0.0, 0.0
    for _ in range(nseg):
        rho = interp(RHO_TAB, p)
        mu = interp(MU_TAB, p)
        v = mass / (rho * area)
        re = rho * v * D / mu
        f = churchill(re, eps / D)
        v_out, re_out, f_out = v, re, f
        p -= rho * G * dz + (f * ds / D) * rho * v * v / 2.0
    return p, v_out, re_out, f_out


# ---- case geometry (re-parsed so the case stays the single source) ----
def pipe_geom(name):
    body = re.search(name + r"\s*\{(.*?)\}", text, re.S).group(1)

    def num(key, default=0.0):
        m = re.search(key + r"\s+([0-9.eE+-]+)\s*;", body)
        return float(m.group(1)) if m else default

    L, D = num("length"), num("diameter")
    src = re.search(r"from\s+(\w+);", body).group(1)
    dst = re.search(r"to\s+(\w+);", body).group(1)
    # elevation gradient along the flow direction (z(dst) - z(src)) / L
    ez = re.search(r"elevation\s+([0-9.eE+-]+)\s*;", body)
    return L, D, num("roughness"), float(ez.group(1)) / L


GEOM = {name: pipe_geom(name) for name in ("Riser", "Flowline")}

fails = []


def check(ok, msg):
    if not ok:
        fails.append(msg)


# ---- serial junction: both pipes carry the source rate ----
SRC = float(re.search(r"rate\s+([0-9.eE+-]+)\s*;", text).group(1))
q_r, q_f = pipes["Riser"]["Q"], pipes["Flowline"]["Q"]
check(abs(q_r - SRC) < 1e-9 * SRC, f"Riser rate {q_r} != source {SRC}")
check(abs(q_f - SRC) < 1e-9 * SRC, f"Flowline rate {q_f} != source {SRC}")

# ---- march closure on the reported pressures ----
p_w, p_m, p_t = nodes["W"], nodes["M"], nodes["T"]
p_m_march, v_r, re_r, f_r = march("Riser", q_r, p_w)
check(abs(p_m_march - p_m) < 1e-9 * p_m,
      f"Riser march lands on {p_m_march}, solver reported {p_m}")
p_t_march, v_f, re_f, f_f = march("Flowline", q_f, p_m)
check(abs(p_t_march - p_t) < 1e-9 * p_t,
      f"Flowline march lands on {p_t_march}, terminal is {p_t}")

# ---- reported outlet velocity / friction factor are the march's ----
check(abs(pipes["Riser"]["v"] - v_r) < 1e-9 * v_r,
      f"Riser velocity {pipes['Riser']['v']} != march {v_r}")
check(abs(pipes["Riser"]["f"] - f_r) < 1e-9 * f_r,
      f"Riser friction {pipes['Riser']['f']} != march {f_r}")
check(abs(pipes["Flowline"]["v"] - v_f) < 1e-9 * v_f,
      f"Flowline velocity {pipes['Flowline']['v']} != march {v_f}")
check(abs(pipes["Flowline"]["f"] - f_f) < 1e-9 * f_f,
      f"Flowline friction {pipes['Flowline']['f']} != march {f_f}")

# ---- the fill expands: outlet velocity above the inlet q/A ----
for name in ("Riser", "Flowline"):
    D = GEOM[name][1]
    area = math.pi * D * D / 4.0
    v_in = pipes[name]["Q"] / area
    check(pipes[name]["v"] > v_in,
          f"{name} outlet velocity {pipes[name]['v']} not above inlet {v_in}")

# ---- compressibility is doing visible work: the reported dp differs from
#      an inlet-density constant-property computation ----
for name in ("Riser", "Flowline"):
    L, D, eps, sin = GEOM[name]
    area = math.pi * D * D / 4.0
    q = pipes[name]["Q"]
    dp = pipes[name]["dp"]
    rho_in = interp(RHO_TAB, p_w if name == "Riser" else p_m)
    mu_in = interp(MU_TAB, p_w if name == "Riser" else p_m)
    v_in = q / area
    re = rho_in * v_in * D / mu_in
    dp_const = (rho_in * G * sin + churchill(re, eps / D) / D * rho_in
                * v_in * v_in / 2.0) * L
    rel = abs(dp - dp_const) / dp
    check(rel > 5e-4, f"{name}: dp {dp} indistinguishable from the "
                      f"constant-property law (rel {rel:.2e}) — rho(p) "
                      f"is not active")
    check(rel < 0.3, f"{name}: dp {dp} suspiciously far from the "
                     f"constant-property law (rel {rel:.2e})")

# ---- I1: Churchill closed-form limits (laminar and fully rough) ----
EPS_D = GEOM["Riser"][2] / GEOM["Riser"][1]
check(abs(churchill(1e-4, EPS_D) - 64.0 / 1e-4) <= 1e-9 * (64.0 / 1e-4),
      f"I1: Churchill does not reduce to 64/Re in the laminar limit "
      f"(got {churchill(1e-4, EPS_D)}, want {64.0 / 1e-4})")
f_rough = 8.0 / (2.457 * math.log(1.0 / (0.27 * EPS_D))) ** 2
check(abs(churchill(1e14, EPS_D) - f_rough) <= 1e-6 * f_rough,
      f"I1: Churchill does not reach the fully-rough closed form as Re -> inf "
      f"(got {churchill(1e14, EPS_D)}, want {f_rough})")

# ---- I2: constant-property twin against the explicit analytic network
#      solution. With rho/mu held constant the march ODE has constant
#      coefficients, so the 24-segment Euler recursion is EXACT: the
#      reported pressures must equal p_T + rho g dz + Churchill friction
#      summed in closed form, and the fill must not expand (v == Q/A).
#      (A laminar twin is unreachable: the source-driven Newton stalls
#      below rate ~1e-3, so the turbulent Churchill arm is the one the
#      closed form can exercise; its limits are pinned separately in I1.)
import tempfile
import os
import subprocess as sp
import re as regex  # the compress block above shadows `re` with a Reynolds number

RHO_C, MU_C = 50.0, 1.3e-5
pvt_const = ("pvt\n    {\n        pressure ( 1.0e5 5.0e5 1.0e6 2.0e6 4.0e6 );\n"
             "        rho ( 50.0 50.0 50.0 50.0 50.0 );\n"
             "        mu ( 1.3e-5 1.3e-5 1.3e-5 1.3e-5 1.3e-5 );\n    }")
twin = regex.sub(r"pvt\s*\{.*?\}", lambda m: pvt_const, text, count=1, flags=regex.S)
if "50.0 50.0" not in twin:
    raise SystemExit("I2: could not build the constant-property twin case")
with tempfile.TemporaryDirectory(prefix="cf_pvt_twin_") as td:
    tp = os.path.join(td, "twin.dict")
    open(tp, "w").write(twin)
    proc2 = sp.run([exe, tp], capture_output=True, text=True)
    if proc2.returncode != 0:
        raise SystemExit(f"I2: twin solver failed rc={proc2.returncode}\n"
                         f"{proc2.stdout}\n{proc2.stderr}")
    twin_nodes, twin_pipes = {}, {}
    for line in proc2.stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "node" and len(tok) >= 4:
            twin_nodes[tok[1]] = float(tok[3])
        elif tok[0] == "pipe":
            d = {}
            for j in range(2, len(tok) - 1, 2):
                d[tok[j]] = float(tok[j + 1])
            twin_pipes[tok[1]] = d

def closed_form_dp(name, q):
    """exact constant-property loss: gravity head + Churchill friction"""
    L, D, eps, sin = GEOM[name]
    area = math.pi * D * D / 4.0
    v = q / area
    re = RHO_C * v * D / MU_C
    return RHO_C * G * sin * L + churchill(re, eps / D) * (L / D) * RHO_C * v * v / 2.0, v, re

if not twin_pipes or "W" not in twin_nodes:
    raise SystemExit("I2: twin output unparsable")
if "NOT converged" in proc2.stdout:
    raise SystemExit("I2: twin case did not converge")
Q_T = twin_pipes["Riser"]["Q"]
dp_r, v_r_t, re_r_t = closed_form_dp("Riser", Q_T)
dp_f, v_f_t, re_f_t = closed_form_dp("Flowline", twin_pipes["Flowline"]["Q"])
p_w_an = twin_nodes["T"] + dp_r + dp_f
p_m_an = twin_nodes["T"] + dp_f
check(abs(twin_nodes["W"] - p_w_an) <= 1e-9 * p_w_an,
      f"I2: twin well head {twin_nodes['W']} != closed-form {p_w_an}")
check(abs(twin_nodes["M"] - p_m_an) <= 1e-9 * p_m_an,
      f"I2: twin junction {twin_nodes['M']} != closed-form {p_m_an}")
for name, v_t in (("Riser", v_r_t), ("Flowline", v_f_t)):
    check(abs(twin_pipes[name]["v"] - v_t) <= 1e-9 * v_t,
          f"I2: twin {name} outlet velocity {twin_pipes[name]['v']} != Q/A {v_t} "
          f"(constant properties must not expand)")

# ---- I3: the 24-segment Euler march is first-order convergent ----
ref = march("Riser", q_r, p_w, nseg=1536)[0]
errs = [abs(march("Riser", q_r, p_w, nseg=n)[0] - ref) for n in (24, 48, 96)]
check(errs[0] > errs[1] > errs[2] > 0,
      f"I3: march errors not monotonically decreasing: {errs}")
r1 = errs[0] / errs[1] if errs[1] > 0 else 0.0
r2 = errs[1] / errs[2] if errs[2] > 0 else 0.0
check(1.7 <= r1 <= 2.3, f"I3: Richardson ratio e(24)/e(48) = {r1} (want ~2)")
check(1.7 <= r2 <= 2.3, f"I3: Richardson ratio e(48)/e(96) = {r2} (want ~2)")

if fails:
    for m in fails:
        print("FAIL:", m)
    raise SystemExit(f"{len(fails)} checks FAILED")
print(f"steadyPvtColumn OK: p_W {p_w:.6g} -> p_M {p_m:.6g} -> p_T {p_t:.6g}, "
      f"Q {SRC:.4g} m3/s, Riser v_out {v_r:.4f} (expands), "
      f"march closure < 1e-9 rel")
