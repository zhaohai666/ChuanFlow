#!/usr/bin/env python3
"""transientPowerLawRheology: the wall law 5.45 - 5.47 proved in the steady
network is now the number the transient momentum row reads, checked by
re-deriving every reported step.

The row for a pipe between two fixed-pressure nodes is one equation,

    rho L/A * (Q_new - Q_prev)/dt = dp - dp_wall(Q_new)

and dp_wall is the shared kernel function (pipeFrictionPa): Churchill's f on the
Metzner-Reed number, whose laminar branch is the exact solution of the declared
tau = K gdot^n law. So the checker does not compare this run against a remembered
number -- it solves the implicit-Euler recursion itself, with the law written out
from the case file, and walks every reported state of all six pipes
against it. Two consequences carry the increment:

  * relaxed to rest, the row lands on the profile integral 5.47 corrected 5.45
    with, Q = (pi d^3/8) n/(3n+1) (dp d/(4KL))^(1/n), which is the same closed
    form the steady kernel is checked against. Before 5.48 the row computed
    `rho u d / mu` by hand, so a declared `flowIndex` was not refused but
    *ignored*: all six traces below would have been the Newtonian one.
  * the damping the row applies is the law's own slope, n dp/Q, not the chord
    dp/Q and not the Newtonian floor 0.2 * 128 mu L/(pi d^4). The fit of the
    measured tail separates the three: at n = 0.5 the chord answer is 2x off,
    which is a factor an assertion can see.

Nine families of judgement:
  1. the declaration is reported -- which is not the same as being used
  2. n = 1, K = mu is inert: the declared and undeclared traces are identical
     strings, step by step
  3. every terminal rate against the closed form of its own law
  4. the sign of the exponent, and how far the old silent-Newtonian answer was
  5. every reported state of every pipe against the implicit-Euler recursion
  6. the tail fit: the row damps with n dp/Q, not with dp/Q and not with the floor
  7. the shape of the approach, including the odd continuation from a negative
     start rate
  8. the laminar branch really is the one in force in every state
  9. the three declaration-face refusals that are still refusals
 10. the fourth refusal -- the one 5.48 wrote for a yield stress -- lifted:
     the same declarations now solve, and land on 5.46's and 5.47's closed
     forms, with the sub-threshold one at exactly zero in every state
"""

import math
import subprocess
import sys

BINARY = sys.argv[1]
CASE = sys.argv[2]
CASE_DIR = CASE.rsplit("/", 1)[0]

RHO = 850.0
MU = 0.1
DP = 1.0e4
DT = 0.002
END = 4.0
NS = int(round(END / DT)) + 1

# name: (length, diameter, n, K, initialFlow); n = 0 means nothing declared, so
# the row must read the Newtonian number -- the same number n = 1, K = mu
# resolves to term for term, which is what family 2 tests as strings
PIPES = {
    "Newt":      (20.0, 0.02, 0.0, 0.0, 0.0),
    "Decl1":     (20.0, 0.02, 1.0, 0.1, 0.0),
    "Thin05":    (20.0, 0.02, 0.5, 1.0, 0.0),
    "Thin07":    (20.0, 0.02, 0.7, 2.0, 0.0),
    "Thick13":   (20.0, 0.02, 1.3, 0.02, 0.0),
    "Thin05Rev": (20.0, 0.02, 0.5, 1.0, -2.0e-6),
}
NONNEWTON = [k for k, p in PIPES.items() if p[2] > 0 and abs(p[2] - 1.0) > 1e-12]
# rho L / A, the inertance every pipe's row carries (no partner measures it here)
INERTIA = RHO * PIPES["Newt"][0] / (0.25 * math.pi * PIPES["Newt"][1] ** 2)

failures = []


def require(cond, what):
    if cond:
        print("ok:  ", what)
    else:
        print("FAIL:", what)
        failures.append(what)


def rel(a, b):
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale > 0 else 0.0


def area_of(diam):
    return 0.25 * math.pi * diam * diam


def run(case_path):
    return subprocess.run([BINARY, case_path], capture_output=True, text=True,
                          cwd=CASE_DIR)


def parse_trace(stdout):
    """The reported wall laws, and the Q trace of every pipe both as floats and
    as the strings the app actually printed."""
    laws = {}
    times = []
    trace = {name: [] for name in PIPES}
    raw = {name: [] for name in PIPES}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 7 and tok[0] == "pipe" and tok[2] == "wall-law":
            laws[tok[1]] = (float(tok[4]), float(tok[6]))
            continue
        if len(tok) >= 2 and tok[0] == "t":
            times.append(float(tok[1]))
            continue
        if len(tok) == 4 and tok[0] == "pipe" and tok[2] == "Q" and tok[1] in PIPES:
            trace[tok[1]].append(float(tok[3]))
            raw[tok[1]].append(tok[3])
    return laws, times, trace, raw


# --- the law, written out from the case file -------------------------------
# pipeFrictionPa is Churchill's Darcy f on the Metzner-Reed number, so on its
# laminar branch it is the exact solution of tau = K gdot^n.

def churchill_f(re, rough=0.0):
    re = max(re, 1e-12)
    a1 = (7.0 / re) ** 0.9 + 0.27 * rough
    a = (2.457 * math.log(1.0 / a1)) ** 16.0
    b = (37530.0 / re) ** 16.0
    return 8.0 * ((8.0 / re) ** 12.0 + (a + b) ** -1.5) ** (1.0 / 12.0)


def re_generalized(speed, diam, n, k):
    if not n > 0:
        return RHO * speed * diam / MU
    kp = k * ((3.0 * n + 1.0) / (4.0 * n)) ** n
    return RHO * speed ** (2.0 - n) * diam ** n / (kp * 8.0 ** (n - 1.0))


def dp_wall(q, diam, length, n, k):
    """Signed friction drop of the declared law, in Pa: the kernel's own
    parentheses, made odd in Q."""
    speed = abs(q) / area_of(diam)
    if speed == 0.0:
        return 0.0
    f = churchill_f(re_generalized(speed, diam, n, k))
    dp = (f * length / diam) * RHO * speed * speed / 2.0
    return dp if q > 0 else -dp


def q_closed_form(diam, length, n, k, dp):
    """The exact laminar discharge of the declared law, by integrating the
    profile: with tau_w = dp d/(4L),
        Q = (pi d^3 / 8) * n/(3n+1) * (tau_w/K)^(1/n)
    which at n = 1, k = mu is 128 mu L Q = pi d^4 dp. The constant with exponent
    1 on (3n+1)/(4n) -- what 5.45 first checked against -- is ((3n+1)/(4n))^(1-1/n)
    off this, 25% in the rate at n = 0.5; family 3 is where that correction is
    now visible on the transient side too."""
    if not n > 0:
        return math.pi * diam ** 4 * dp / (128.0 * MU * length)
    tau_w = dp * diam / (4.0 * length)
    return (math.pi * diam ** 3 / 8.0) * n / (3.0 * n + 1.0) \
        * (tau_w / k) ** (1.0 / n)


def solve_step(q_prev, diam, length, n, k):
    """One implicit-Euler step of the single row: monotone in Q, so bisect."""
    i_dt = INERTIA / DT

    def g(q):
        return i_dt * (q - q_prev) + dp_wall(q, diam, length, n, k) - DP

    step = DP / i_dt + abs(q_prev) + 1.0
    lo, hi = q_prev - step, q_prev + step
    for _ in range(200):
        if g(lo) <= 0.0:
            break
        lo -= step
    for _ in range(200):
        if g(hi) >= 0.0:
            break
        hi += step
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo <= 1e-16 * max(abs(mid), 1e-30):
            break
    return 0.5 * (lo + hi)


def tail_slope(qs, q_ss):
    """The slope (Pa.s/m3) the row is damping with, read off the approach where
    it is neither still nonlinear nor already at the printed resolution:
    |delta|/Q_ss in [3e-5, 5e-2] -- the lower bound is what lets the fastest pipe
    (n = 0.7, whose time constant is 3.8 dt) still hand over twenty samples. Implicit Euler on the linearised row gives
    delta_{i+1}/delta_i = 1/(1 + dt s/I) exactly, so s comes back without a log
    approximation; the window's endpoints give the ratio as a geometric mean.
    Returns None when the window is too short to fit."""
    idx = [i for i, q in enumerate(qs)
           if q_ss > 0 and 3e-5 <= abs(q_ss - q) / q_ss <= 5e-2]
    if len(idx) < 8:
        return None
    first, last = idx[0], idx[-1]
    if last == first:
        return None
    r = ((q_ss - qs[last]) / (q_ss - qs[first])) ** (1.0 / (last - first))
    if not (0.0 < r < 1.0):
        return None
    return (INERTIA / DT) * (1.0 / r - 1.0)


proc = run(CASE)
if proc.returncode != 0:
    raise SystemExit(f"solver failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")

laws, times, trace, raw = parse_trace(proc.stdout)

DECLARED = {k for k, p in PIPES.items() if p[2] > 0}

# --- 1. the declaration is reported ---------------------------------------
require(set(laws) == DECLARED,
        f"the app reports a wall law for exactly the {len(DECLARED)} declared "
        f"pipes and not for the undeclared one: {sorted(laws)}")
require(all(laws[k][0] == PIPES[k][2] and laws[k][1] == PIPES[k][3] for k in laws),
        "the reported n and K are the case file's own numbers read back without a "
        "conversion -- which is the whole of what 5.45 to 5.47 could only claim "
        "about the steady solver")
require(all(len(trace[k]) == NS for k in PIPES),
        f"every pipe reports {NS} states (t = 0 plus {NS - 1} steps)")

# --- 2. declaring the Newtonian number changes no string -------------------
same = sum(1 for a, b in zip(raw["Newt"], raw["Decl1"]) if a == b)
require(same == NS,
        f"n = 1 with K = mu is inert in the printed trace: {same}/{NS} states of "
        f"Newt and Decl1 are the same string")

# --- 3. terminal rate against the closed form of its own law --------------
q_ss = {}
worst_q, worst_q_name = 0.0, ""
for name, (length, diam, n, k, _) in PIPES.items():
    q_ss[name] = q_closed_form(diam, length, n, k, DP)
    e = rel(trace[name][-1], q_ss[name])
    if e > worst_q:
        worst_q, worst_q_name = e, name
    require(e < 1e-9,
            f"{name}: relaxed to rest, the transient row lands on n="
            f"{'-' if n == 0 else n}'s closed form, Q = {trace[name][-1]:.10e} "
            f"vs {q_ss[name]:.10e} (rel {e:.3e})")
require(worst_q < 1e-8,
        f"the same closed form 5.47 measures on the steady side is met on the "
        f"transient side to {worst_q:.3e} (worst {worst_q_name}): one law, two "
        f"kernels, one number")
worst_root = max(rel(dp_wall(q_ss[name], PIPES[name][1], PIPES[name][0],
                             PIPES[name][2], PIPES[name][3]), DP)
                 for name in PIPES)
require(worst_root < 1e-12,
        f"and the gap above is the solver's stop, not the law's: the closed "
        f"form dropped back into the kernel's own friction returns the asked "
        f"1.0e4 Pa to {worst_root:.3e}, so Churchill's laminar branch and the "
        f"profile integral are the same statement here -- the leftover is each "
        f"step's Newton exit (1e-11 on a residual scaled by rho L/(A dt)) "
        f"settling as far below the fixed point as the relaxation lets it pile "
        f"up")

# --- 4. the sign of the exponent, and the size of the old silent answer ---
ratios = {name: q_ss[name] / q_ss["Newt"] for name in PIPES}
require(ratios["Thin07"] < ratios["Thin05"] < 1.0 < ratios["Thick13"],
        f"at one drop and one geometry the steady rates order by the exponent: "
        f"n = 0.7 gives {ratios['Thin07']:.4f} of the Newtonian rate, n = 0.5 "
        f"{ratios['Thin05']:.4f}, n = 1.3 {ratios['Thick13']:.4f}")
min_gap = min(abs(ratios[name] - 1.0) for name in NONNEWTON)
require(min_gap > 0.25,
        f"a declared law is not a decoration: the nearest non-Newtonian terminal "
        f"rate is {min_gap * 100:.1f}% away from the Newtonian one, which is how "
        f"far the six traces were forced apart -- the row before 5.48 would have "
        f"printed one answer for all of them")

# --- 5. every reported state against the recursion ------------------------
worst_trace, worst_trace_name, worst_step = 0.0, "", 0
for name, (length, diam, n, k, q0) in PIPES.items():
    qs = [q0]
    for _ in range(1, len(trace[name])):
        qs.append(solve_step(qs[-1], diam, length, n, k))
    e_pipe = max(rel(m, t) for m, t in zip(qs, trace[name]))
    if e_pipe > worst_trace:
        worst_trace, worst_trace_name = e_pipe, name
        worst_step = max(i for i, (m, t) in enumerate(zip(qs, trace[name]))
                         if rel(m, t) == e_pipe)
    require(e_pipe < 1e-7,
            f"{name}: all {len(qs)} states reproduce the implicit-Euler row whose "
            f"friction is the declared law (worst rel {e_pipe:.3e})")
require(worst_trace < 1e-7,
        f"the whole trajectory of {len(PIPES)} pipes x {NS} states is the row's "
        f"own algebra, not a shape that happens to end in the right place "
        f"(worst {worst_trace:.3e} at step {worst_step} of {worst_trace_name})")

# --- 6. the damping is the law's slope, not the chord and not the floor ----
FLOOR = 0.2 * 128.0 * MU * PIPES["Newt"][0] / (math.pi * 0.02 ** 4)
for name in ("Newt", "Thin05", "Thin07", "Thick13", "Thin05Rev"):
    n = PIPES[name][2]
    s_meas = tail_slope(trace[name], q_ss[name])
    s_law = (n if n > 0 else 1.0) * DP / q_ss[name]
    if s_meas is None:
        require(False, f"{name}: the tail of the approach is long enough to fit")
        continue
    require(rel(s_meas, s_law) < 0.05,
            f"{name}: the row damps with slope {s_meas:.6e} Pa.s/m3 against the "
            f"law's own n dp/Q = {s_law:.6e} for n={'-' if n == 0 else n} "
            f"(rel {rel(s_meas, s_law):.3e}, tau = "
            f"{INERTIA / s_meas:.4f} s)")
    if name in NONNEWTON:
        s_chord = DP / q_ss[name]
        require(rel(s_meas, s_chord) > 0.10 and rel(s_meas, FLOOR) > 0.10,
                f"{name}: the chord dp/Q = {s_chord:.6e} and the Newtonian floor "
                f"0.2*128 mu L/(pi d^4) = {FLOOR:.6e} both miss the measured "
                f"{s_meas:.6e} by more than 10% (chord by "
                f"{rel(s_meas, s_chord) * 100:.1f}%), so neither is what the row "
                f"actually used")

# --- 7. the shape of the approach ----------------------------------------
from_rest = [trace[k] for k in PIPES if PIPES[k][4] == 0.0]
require(all(qs[i + 1] >= qs[i] for i in range(len(qs) - 1) for qs in from_rest),
        f"the {len(from_rest)} runs started at rest approach their rates without "
        f"ever stepping back -- the row is a relaxation, not an oscillation")
require(all(max(qs) <= q_ss[name] * (1.0 + 1e-9)
            for name in PIPES if PIPES[name][4] == 0.0 for qs in [trace[name]]),
        "and none overshoots: no reported state exceeds its law's own steady "
        "rate, which is all the imposed drop can pay")
signs = [1 if q > 0 else (-1 if q < 0 else 0) for q in trace["Thin05Rev"]]
crossings = sum(1 for a, b in zip(signs, signs[1:]) if a and b and a != b)
require(crossings == 1 and signs[-1] == 1,
        f"Thin05Rev is released flowing the wrong way and turns around exactly "
        f"once ({crossings} sign change): the law continues oddly in Q, "
        f"dp_wall(-q) = -dp_wall(q), and only a marching row can show it")
require(rel(trace["Thin05Rev"][-1], trace["Thin05"][-1]) < 1e-9,
        f"and it ends where the same law started from rest ends, "
        f"{trace['Thin05Rev'][-1]:.10e} vs {trace['Thin05'][-1]:.10e}, so the "
        f"row's fixed point belongs to the law and not to the initial condition")

# --- 8. the laminar branch is the one in force ----------------------------
worst_re, worst_re_name = 0.0, ""
for name, (length, diam, n, k, _) in PIPES.items():
    for q in trace[name]:
        re = re_generalized(abs(q) / area_of(diam), diam, n, k)
        if re > worst_re:
            worst_re, worst_re_name = re, name
require(worst_re < 2100.0,
        f"every reported state of every pipe is on Churchill's laminar branch "
        f"(largest generalized number {worst_re:.4f}, on {worst_re_name}), so "
        f"family 3's closed form is the exact solution rather than a fit")

# --- 9. the declaration face refuses what this kernel does not own --------
REFUSALS = [
    ("case_rejectNegativeIndex.dict", "negative",
     "a negative index is no law at all (0 means nothing was declared)"),
    ("case_rejectNoConsistency.dict", "without a positive consistency",
     "at n != 1 the fluid's mu is not the same quantity as K"),
    ("case_rejectCoupledRheology.dict", "cannot carry a rheology",
     "a coupled element's loss is already its partner's statement"),
]
for fname, needle, why in REFUSALS:
    p = run(f"{CASE_DIR}/{fname}")
    msg = (p.stdout + p.stderr).replace("[chflow] ", "")
    line = next((ln for ln in msg.splitlines() if needle in ln), "")
    require(p.returncode != 0 and needle in msg,
            f"{fname}: refused, naming {why} -- '{line}'")

# --- 10. the refusal 5.48 wrote down for a threshold is lifted by 5.49 -----
# The fourth file of family 9 used to be case_rejectYieldStress.dict; the static
# branch it refused for is now in the row, so the same declarations solve and
# land on the steady kernel's own yielded closed forms. Only the third statement
# is new physics: below the threshold the answer is *exactly* zero, which is the
# shape 5.46 has and a marching row could not.
YIELD_PIPES = {
    "Bing": (0.02, 20.0, 0.0, 0.0, 0.5),   # tau_0 only -> Buckingham (5.46)
    "Hb":   (0.02, 20.0, 0.5, 1.0, 0.5),   # with an index -> profile integral
    "Stuck": (0.02, 20.0, 0.0, 0.0, 2.6),  # dp_static 1.04e4 > the asked drop
}


def g_buckingham(xi):
    return 1.0 - (4.0 / 3.0) * xi + xi ** 4 / 3.0


def q_yield_closed_form(diam, length, n, k, tau0, dp):
    """The exact laminar discharge of the yielded law at an imposed drop, with
    tau_w = dp d/(4L): Buckingham at n = 1, k = mu and 5.47's profile integral
    otherwise, both written as the kernel writes them."""
    tw = dp * diam / (4.0 * length)
    if tw <= tau0:
        return 0.0
    if not n > 0:
        xi = tau0 / tw
        return (math.pi * diam ** 3 / (32.0 * MU)) * tw * g_buckingham(xi)
    d = tw - tau0
    r = 0.5 * diam
    return (math.pi * r ** 3 * tw ** -3 * k ** (-1.0 / n) * n * d ** (1.0 + 1.0 / n)
            * (tau0 * tau0 / (n + 1.0) + 2.0 * tau0 * d / (2.0 * n + 1.0)
               + d * d / (3.0 * n + 1.0)))


p = run(f"{CASE_DIR}/case_yieldNowSolves.dict")
require(p.returncode == 0,
        "case_yieldNowSolves.dict: 5.48's refusal is lifted -- the same "
        "declarations the row used to reject now solve, with a threshold on all "
        "three pipes")

ylaws, ytimes, ytrace, yraw = {}, [], {}, {}
for line in p.stdout.splitlines():
    tok = line.replace("[chflow] ", "").split()
    if len(tok) == 4 and tok[0] == "pipe" and tok[2] == "Q" and tok[1] in YIELD_PIPES:
        ytrace[tok[1]] = ytrace.get(tok[1], []) + [float(tok[3])]
    elif len(tok) >= 7 and tok[0] == "pipe" and tok[2] == "wall-law":
        ylaws[tok[1]] = (float(tok[4]), float(tok[6]))
require(all(len(v) == NS for v in ytrace.values()) and len(ytrace) == 3,
        f"and all three report the same {NS} states as the six power-law traces")
worst_yield = 0.0
for name, (diam, length, n, k, tau0) in YIELD_PIPES.items():
    want = q_yield_closed_form(diam, length, n, k, tau0, DP)
    e = rel(ytrace[name][-1], want)
    worst_yield = max(worst_yield, e)
    law = ("Buckingham" if not n > 0 else "the HB integral")
    at = f"Q = {ytrace[name][-1]:.10e} vs {want:.10e} (rel {e:.3e})"
    if want == 0.0:
        at = f"its wall needs 1.04e4 Pa against the 1.0e4 Pa on offer, so the " \
             f"closed form's own dead branch says 0 and the trace says 0"
    require(e < 1e-8,
            f"{name}: ended on {law} of the same 1.0e4 Pa -- {at}")
require(ytrace["Stuck"] == [0.0] * NS,
        f"the sub-threshold pipe is exactly zero in all {NS} states -- not small, "
        f"not below a tolerance: the static branch is the row, so this is the "
        f"sentence 5.48 could only write down as a refusal")
require(ylaws.get("Hb") == (0.5, 1.0) and "Bing" not in ylaws,
        f"the lifted declaration still reads back its own n and K, and only where "
        f"one was declared ({sorted(ylaws)})")

print()
print(f"transientPowerLawRheology: terminal rate vs closed form worst "
      f"{worst_q:.3e} (the closed form into the kernel's own law returns the "
      f"drop to {worst_root:.3e}), whole-trace reproduction worst "
      f"{worst_trace:.3e}, largest generalized Re {worst_re:.4f}, "
      f"declared-vs-Newtonian gap at least {min_gap * 100:.1f}%, "
      f"yielded-terminal-vs-steady-closed-form worst {worst_yield:.3e}")

if failures:
    print(f"FAILED: {len(failures)} assertion(s)")
    sys.exit(1)
print("all assertions passed")
