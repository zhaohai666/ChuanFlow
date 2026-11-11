#!/usr/bin/env python3
"""steadyPowerLawRheology: a power-law fill in the steady network, checked
against the exact laminar closed form rather than against a correlation.

The kernel makes one substitution -- the Metzner-Reed generalized Reynolds
number in place of rho*u*d/mu inside the same Churchill call:

    Re_MR = rho u^(2-n) d^n / (K' 8^(n-1)),   K' = K ((3n+1)/(4n))^n

and the claim is that Churchill's laminar branch f = 64/Re then *is* the
solution for tau = K gdot^n, not an approximation of it:

    dp = (L/d) K (3n+1)/n (8u/d)^n
     Q = (pi d^3 / 32) [ n dp d / ((3n+1) K L) ]^(1/n)

which at n = 1, K = mu is 128 mu L Q = pi d^4 dp. Both ends of every pipe are
fixed-pressure nodes, so each pipe carries the imposed dp alone and each of the
nine statements below is independent of the others and of the Newton solve.

Eight families of judgement:
  1. the declaration read back (nine pipes, one drop, two files)
  2. Q against the closed form, every pipe, both driving differences
  3. the equivalence control: n = 1, K = mu is bit-identical to undeclared
  4. the laminar branch is the one in force (f * Re == 64, Re well inside)
  5. the three structural exponents: dp -> 1/n, d -> (3n+1)/n, L -> -1
  6. the wall stress is not a material property (same dp, d, L -> same tauw)
  7. the parse-time refusal that is still a refusal, plus the coupled-edge
     declaration 5.63 lifted: accepted, but the 1D trace is bit-identical to a
     twin without the three keys (the declaration belongs to the 3D partner)
"""

import math
import os
import re
import subprocess
import sys
import tempfile

BINARY = sys.argv[1]
CASE = sys.argv[2]
CASE_DIR = CASE.rsplit("/", 1)[0]
DP_BASE = 1.0e4           # p(A) 3e4 - p(B) 2e4
DP_DOUBLE = 2.0e4         # case_dpDouble: p(A) 4e4

RHO = 850.0
MU = 0.1

# name: (length, diameter, n, K); n = 0 means the pipe declares no rheology
PIPES = {
    "Newt":     (20.0, 0.02, 0.0, 0.0),
    "PowN1":    (20.0, 0.02, 1.0, MU),
    "Thin05":   (20.0, 0.02, 0.5, 1.0),
    "Thin07":   (20.0, 0.02, 0.7, 2.0),
    "Thick13":  (20.0, 0.02, 1.3, 0.02),
    "D2xNewt":  (20.0, 0.04, 0.0, 0.0),
    "D2xThin":  (20.0, 0.04, 0.5, 1.0),
    "LongNewt": (40.0, 0.02, 0.0, 0.0),
    "LongThin": (40.0, 0.02, 0.5, 1.0),
}
NEWTONIAN = ("Newt", "PowN1", "D2xNewt", "LongNewt")
THIN = ("Thin05", "Thin07", "Thick13", "D2xThin", "LongThin")

failures = []


def require(cond, what):
    if cond:
        print("ok:  ", what)
    else:
        print("FAIL:", what)
        failures.append(what)


def run(case_path):
    p = subprocess.run([BINARY, case_path], capture_output=True, text=True,
                       cwd=CASE_DIR)
    return p


def parse(stdout):
    """`pipe <name> Q <q> v <u> Re <re> f <f> dp <dp> ...`, token-paired: the
    field order is the app's, and reading it as key/value pairs means a new
    field printed after dp cannot silently shift the numbers."""
    out = {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) < 2 or tok[0] != "pipe" or tok[1] not in PIPES:
            continue
        name = tok[1]
        fields = {}
        i = 2
        while i + 1 < len(tok):
            key, val = tok[i], tok[i + 1]
            if key in ("Q", "v", "Re", "f", "dp"):
                try:
                    fields[key] = float(val)
                except ValueError:
                    break
            i += 2
        if not {"Q", "v", "Re", "f", "dp"} <= set(fields) or name in out:
            continue
        out[name] = fields
    return out


def q_closed_form(dp, length, diam, n, k):
    """Exact laminar discharge of a power-law fluid in a circular pipe, from
    the integrated profile: with tau_w = dp d/(4L),
        Q = pi R^3 n/(3n+1) (tau_w/K)^(1/n)
    which at n = 1, k = mu is 128 mu L Q = pi d^4 dp. The Metzner-Reed
    constant with the exponent 1 on (3n+1)/(4n) -- the turbulent-collapse
    convention -- is what 5.45 originally checked against; the profile
    integral is the exact one and the two differ by (4n/(3n+1))^(1-1/n),
    25% in the rate at n = 0.5 (measured in 5.47)."""
    if n == 0.0:
        return math.pi * diam ** 4 * dp / (128.0 * MU * length)
    tau_w = dp * diam / (4.0 * length)
    return (math.pi * diam ** 3 / 8.0) * n / (3.0 * n + 1.0) * (tau_w / k) ** (1.0 / n)


def re_mr(u, diam, n, k):
    if n == 0.0:
        return RHO * u * diam / MU
    kp = k * ((3.0 * n + 1.0) / (4.0 * n)) ** n
    return RHO * u ** (2.0 - n) * diam ** n / (kp * 8.0 ** (n - 1.0))


def kp_of(name):
    """K' = K((3n+1)/4n)^n, the consistency the generalized number is written
    with: Rabinowitsch's wall rate cubed... rather, taken to the n."""
    n, k = PIPES[name][2], PIPES[name][3]
    return k * ((3.0 * n + 1.0) / (4.0 * n)) ** n


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


base = run(f"{CASE_DIR}/case.dict")
dbl = run(f"{CASE_DIR}/case_dpDouble.dict")
require(base.returncode == 0 and dbl.returncode == 0,
        f"both driving differences solve and converge (rc {base.returncode},"
        f" {dbl.returncode})")
B = parse(base.stdout)
D = parse(dbl.stdout)

# --- 1. the declaration read back -------------------------------------------
require(len(B) == len(PIPES) and len(D) == len(PIPES),
        f"all nine pipes report a rate, a velocity, a Reynolds number and a"
        f" drop in both files ({len(B)} and {len(D)} parsed)")
require(all(abs(B[n]["dp"] - DP_BASE) < 1e-9 for n in PIPES)
        and all(abs(D[n]["dp"] - DP_DOUBLE) < 1e-9 for n in PIPES),
        f"every pipe carries exactly the imposed difference -- the nodes are"
        f" fixed, so dp is {DP_BASE:.0f} Pa and {DP_DOUBLE:.0f} Pa and the"
        f" closed form needs no Newton solve to be exact")

# --- 2. against the closed form ---------------------------------------------
worst_q = max(rel(B[n]["Q"], q_closed_form(DP_BASE, *PIPES[n][:2],
                                           *[v for v in PIPES[n][2:]]))
              for n in PIPES)
worst_q2 = max(rel(D[n]["Q"], q_closed_form(DP_DOUBLE, PIPES[n][0],
                                            PIPES[n][1], PIPES[n][2],
                                            PIPES[n][3]))
               for n in PIPES)
require(worst_q < 1e-9 and worst_q2 < 1e-9,
        f"all nine rates are the exact power-law solution, to "
        f"{worst_q:.3e} relative at dp={DP_BASE:.0f} and {worst_q2:.3e} at"
        f" {DP_DOUBLE:.0f} -- the rheology is not a fitted correction here")
worst_v = max(rel(B[n]["v"], B[n]["Q"] / (math.pi * PIPES[n][1] ** 2 / 4))
              for n in PIPES)
require(worst_v < 1e-10,
        f"and the reported velocity is the rate over the bore area to "
        f"{worst_v:.3e} -- that residual is the printed precision, not the"
        f" physics, and the diameter the case declares is the one the kernel"
        f" used")

# --- 3. the equivalence control ---------------------------------------------
same = all(B["Newt"][k] == B["PowN1"][k] and D["Newt"][k] == D["PowN1"][k]
           for k in ("Q", "v", "Re", "f", "dp"))
require(same,
        "declaring n = 1 with K = mu is bit-identical to declaring nothing at"
        " all, in both files and in every field -- Re_MR collapses to rho u d /"
        " mu term for term (4/4, 8^0, x^1), so the generalized number costs"
        " nothing and the other 98 tests cannot be moved by it")

# --- 4. the laminar branch is what is in force -------------------------------
worst_f = max(abs(B[n]["f"] * B[n]["Re"] - 64.0) for n in PIPES)
worst_f2 = max(abs(D[n]["f"] * D[n]["Re"] - 64.0) for n in PIPES)
worst_re = max(max(B[n]["Re"] for n in PIPES), max(D[n]["Re"] for n in PIPES))
worst_f_all = max(worst_f, worst_f2)
require(worst_f_all < 1e-9,
        f"f * Re_MR == 64 on all eighteen states (worst {worst_f_all:.3e})"
        f" -- the closed form above is the friction law the"
        f" kernel actually used, not one the checker assumed")
require(worst_re < 1000.0,
        f"and the whole run stays in that branch: the largest generalized"
        f" Reynolds number anywhere is {worst_re:.4g} (the doubled drop takes"
        f" the wide Newtonian bore to it), an order of magnitude below the"
        f" turnover, so no pipe leans on the turbulent half of Churchill --"
        f" which is why the assertions above are exact rather than"
        f" tolerance-based")
worst_re_calc = max(rel(B[n]["Re"], re_mr(B[n]["v"], PIPES[n][1], PIPES[n][2],
                                          PIPES[n][3]))
                    for n in PIPES)
require(worst_re_calc < 1e-11,
        f"the reported Re is the generalized number itself, to "
        f"{worst_re_calc:.3e} (the twelve printed digits floor any such "
        f"comparison near 5e-13, and pow() of the corrected constant adds a "
        f"couple of them) -- the log says which Reynolds number was used, and "
        f"it is the generalized one (K' = K((3n+1)/4n)^n), not the Newtonian "
        f"one with a relabelled viscosity")

# --- 5. the three structural exponents ---------------------------------------
ratio_d_newt = B["D2xNewt"]["Q"] / B["Newt"]["Q"]
ratio_d_thin = B["D2xThin"]["Q"] / B["Thin05"]["Q"]
require(abs(ratio_d_newt - 16.0) < 1e-9 and abs(ratio_d_thin - 32.0) < 1e-9,
        f"doubling the bore multiplies the rate by {ratio_d_newt:.6f} for the"
        f" Newtonian and {ratio_d_thin:.6f} at n = 0.5: the exponent is"
        f" (3n+1)/n, so 4 turns into 5 -- a shear-thinning line is sized by a"
        f" different power of the diameter than a Newtonian one")
ratio_l_newt = B["LongNewt"]["Q"] / B["Newt"]["Q"]
ratio_l_thin = B["LongThin"]["Q"] / B["Thin05"]["Q"]
require(abs(ratio_l_newt - 0.5) < 1e-9 and abs(ratio_l_thin - 0.25) < 1e-9,
        f"doubling the length divides it by {1/ratio_l_newt:.6f} Newtonian and"
        f" {1/ratio_l_thin:.6f} at n = 0.5 -- the length exponent stays 1, so"
        f" what the rheology moves is the dp exponent and the d exponent, not"
        f" the dependence on L")
worst_dp = max(rel(D[n]["Q"] / B[n]["Q"],
                   (DP_DOUBLE / DP_BASE) ** (1.0 / (1.0 if PIPES[n][2] == 0
                                                    else PIPES[n][2])))
               for n in PIPES)
require(worst_dp < 1e-9,
        f"doubling the drop multiplies the rate by 2^(1/n) for every one of the"
        f" nine, worst {worst_dp:.3e} relative: "
        + ", ".join(f"n={PIPES[n][2] or 1:g} -> "
                    f"{D[n]['Q']/B[n]['Q']:.4f}" for n in NEWTONIAN[:1] + THIN))

# --- 6. the wall stress is not a material property --------------------------
SAME_BORE = ("Newt", "PowN1", "Thin05", "Thin07", "Thick13")
taus = {n: B[n]["dp"] * PIPES[n][1] / (4.0 * PIPES[n][0]) for n in PIPES}
require(max(abs(taus[n] - taus["Thin05"]) for n in SAME_BORE) == 0.0,
        f"the five pipes that share a bore and a length ask for the *same* wall"
        f" stress {taus['Thin05']:.6g} Pa at this drop -- K and n do not enter"
        f" tau_w = dp d / 4L, they enter the shear rate needed to hold it")
worst_tk = max(rel(kp_of(n) * (8.0 * B[n]["v"] / PIPES[n][1]) ** PIPES[n][2],
                   taus[n])
               for n in THIN)
require(worst_tk < 1e-9,
        f"and each rate does hold its own stress: K'(8u/d)^n reproduces"
        f" tau_w = dp d/4L on all five rheology pipes to {worst_tk:.3e}, so the"
        f" fluid is free to choose the shear rate but not the wall stress --"
        f" which is why the shear-thinning pair {taus['Thin05']:.6g} vs"
        f" {taus['D2xThin']:.6g} Pa differs by geometry alone")

# --- 7. refusals ------------------------------------------------------------
r1 = run(f"{CASE_DIR}/case_rejectNoConsistency.dict")
require(r1.returncode != 0
        and "flow index" in (r1.stdout + r1.stderr).lower()
        and "consistency" in (r1.stdout + r1.stderr).lower(),
        "a flow index without a consistency is refused at parse time, not"
        " defaulted to mu: at n != 1 the fluid's viscosity is a different unit"
        " (Pa.s^n), so a default would be a silent wrong number")

# 5.63 opened this cell: a coupled edge MAY declare the fluid of its 3D
# partner region. The 1D parser moves the declaration to
# NetworkPipe::partnerRheology and clears its own wall-law fields, so the 1D
# trace must be bit-identical to a twin case without the three keys -- the
# no-double-counting rule 5.43 established for the shape account still holds.
src = open(f"{CASE_DIR}/case_rejectCoupled.dict").read()
control = re.sub(r"\s*(flowIndex|consistency|yieldStress) [^;]+;", "", src)
with tempfile.TemporaryDirectory() as td:
    ctrl = os.path.join(td, "case_coupledControl.dict")
    open(ctrl, "w").write(control)
    r2 = run(f"{CASE_DIR}/case_rejectCoupled.dict")
    r2c = run(ctrl)
require(r2.returncode == 0 and r2c.returncode == 0,
        "a coupled element that declares a rheology is accepted: the keys now"
        " describe the fluid inside the partner region, not a wall law the 1D"
        " row evaluates")
require(r2.stdout == r2c.stdout,
        "and the 1D trace is bit-identical to a twin without the declaration:"
        " the coupling edge still charges only R*q, so no friction term was"
        " silently added by the declaration (5.43 / 5.63)")

print()
print(f"{len(failures)} failure(s)")
sys.exit(1 if failures else 0)
