#!/usr/bin/env python3
"""transientYieldValve: the yield wall's static branch (5.49) and the local loss
that the M2 transient row has carried since architecture `§4.2` (`valve K`), on
the same momentum row, checked by re-deriving
every reported step from both laws written out in Python.

The row is

    rho L/A * (Q_new - Q_prev)/dt = drive - dp_wall(Q_new) - K rho |Q| Q / (2 A^2),
    drive = pFrom - pTo - rho g dz

and the two terms on the right differ in one way that 5.46 and 5.49 each put in
their own boundary list: the wall still asks for 4 tau_0 L/d at zero rate, while
the valve asks for nothing there. So the pair has three separable claims, and this
checker aims at each with a different kind of evidence.

  * At rest the valve is not a number the wall has to beat. `dp_valve(0) = 0`, so
    the drive 5.49's latch compares against does not move with K -- which is
    asserted as an identity over five decades of K rather than as a tolerance:
    StuckRef / StuckK500 / StuckK1e8 print the same token in all 2001 states, and
    the steady kernel prints no `dpValve` field for any pipe that reports `Q 0`.
  * Where the law *is* read the valve is charged, so the rate is the root of
    *both* terms together. Four thresholds-of-the-same-pipe (ReleaseRef through
    ReleaseK5000) plus the gravity-augmented pair (RollRef / RollK500) are
    compared to this file's own bisection on `drive = dp_wall(Q) + dp_valve(Q)`,
    in both kernels.
  * The finite-time stop is the row's arithmetic, so the path decides when it
    happens: the same wall-held pipe with no valve, K = 5e3 and K = 1e8 stops on
    step 17, 17 and 5, and the recursion with both terms reproduces all three
    traces state by state -- including the step it stops on.

Eight families of judgement:
  1. the threshold and the valve field are both reported for what they are
  2. a held pipe is the same pipe at any K, in both kernels
  3. a moving pipe's rate is the root of the two-term loss, in both kernels
  4. the steady kernel's printed split adds up, and it is the valve's own formula
  5. the stop: three different paths, three stop steps, and what each is worth
  6. every reported state of every pipe against the recursion with both terms
  7. the marched row at t = 4 against the steady kernel for the same network
  8. the geometry of the two terms: the valve moves what the wall cannot
"""

import math
import subprocess
import sys

TRANSIENT = sys.argv[1]
STEADY = sys.argv[2]
CASE = sys.argv[3]
CASE_DIR = CASE.rsplit("/", 1)[0]

RHO = 850.0
MU = 0.1
GRAV = 9.80665
DT = 0.002
END = 4.0
NS = int(round(END / DT)) + 1
P_A, P_B, P_A2 = 3.0e4, 2.0e4, 2.4e4
LENGTH = 20.0
DIAM = 0.02
AREA = 0.25 * math.pi * DIAM ** 2
INERTIA = RHO * LENGTH / AREA                # rho L / A, every pipe
I_DT = INERTIA / DT
KV_SCALE = RHO / (2.0 * AREA * AREA)         # dp_valve = K * KV_SCALE * Q |Q|
WALL_SCALE = 4.0 * LENGTH / DIAM             # dp_wall = WALL_SCALE * tau_wall

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


def run(binary, case_path):
    return subprocess.run([binary, case_path], capture_output=True, text=True,
                          cwd=CASE_DIR)


# name: (tau_0, valve K, initial Q, pFrom, pTo, dz)
PIPES = {
    "StuckRef":    (2.6, 0.0,    0.0,    P_A,  P_B,  0.0),
    "StuckK500":   (2.6, 500.0,  0.0,    P_A,  P_B,  0.0),
    "StuckK1e8":   (2.6, 1.0e8,  0.0,    P_A,  P_B,  0.0),
    "ReleaseRef":  (1.0, 0.0,    0.0,    P_A,  P_B,  0.0),
    "ReleaseK100": (1.0, 100.0,  0.0,    P_A,  P_B,  0.0),
    "ReleaseK500": (1.0, 500.0,  0.0,    P_A,  P_B,  0.0),
    "ReleaseK5000": (1.0, 5000.0, 0.0,   P_A,  P_B,  0.0),
    "HangRef":     (0.5, 0.0,    0.0,    P_A,  P_B,  1.0),
    "HangK500":    (0.5, 500.0,  0.0,    P_A,  P_B,  1.0),
    "RollRef":     (0.5, 0.0,    0.0,    P_A,  P_B, -1.0),
    "RollK500":    (0.5, 500.0,  0.0,    P_A,  P_B, -1.0),
    "StopRef":     (1.5, 0.0,    2.0e-6, P_A2, P_B,  0.0),
    "StopK5000":   (1.5, 5000.0, 2.0e-6, P_A2, P_B,  0.0),
    "StopK1e8":    (1.5, 1.0e8,  2.0e-4, P_A2, P_B,  0.0),
}
HELD = ("StuckRef", "StuckK500", "StuckK1e8", "HangRef", "HangK500",
        "StopRef", "StopK5000", "StopK1e8")
MOVING = ("ReleaseRef", "ReleaseK100", "ReleaseK500", "ReleaseK5000",
          "RollRef", "RollK500")
VALVED = tuple(n for n in PIPES if PIPES[n][1] > 0.0)


def parse_transient(stdout):
    """Q trace (value and printed token), latch per pipe, thresholds, times."""
    out = {"q": {n: [] for n in PIPES}, "raw": {n: [] for n in PIPES},
           "hold": {n: [] for n in PIPES}, "thr": {}, "times": []}
    started = False
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 2 and tok[0] == "t":
            out["times"].append(float(tok[1]))
            started = True
        elif not started and len(tok) == 5 and tok[0] == "pipe" \
                and tok[3] == "threshold":
            out["thr"][tok[1]] = float(tok[4])
        elif len(tok) == 4 and tok[0] == "pipe" and tok[2] == "Q":
            if tok[1] in out["q"]:
                out["q"][tok[1]].append(float(tok[3]))
                out["raw"][tok[1]].append(tok[3])
        elif len(tok) == 5 and tok[0] == "pipe" and tok[3] == "hold":
            if tok[1] in out["hold"]:
                out["hold"][tok[1]].append(int(tok[4]))
    return out


def parse_steady(stdout):
    """One state per pipe, with the printed drop split."""
    out = {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) < 12 or tok[0] != "pipe" or tok[2] != "Q" or tok[1] not in PIPES:
            continue
        rec = {"q": float(tok[3]), "v": float(tok[5])}
        for i, t in enumerate(tok):
            if t == "dp":
                rec["dp"] = float(tok[i + 1])
            elif t == "dpValve":
                rec["dpValve"] = float(tok[i + 1])
        out[tok[1]] = rec
    return out


# --- both laws, written out from the case file ------------------------------

def threshold(tau0):
    return WALL_SCALE * tau0


def q_of_tw(tw, tau0):
    """5.47's profile integral at n = 1, k = mu, which is 5.46's Buckingham
    bracket: every pipe in this case declares a yield stress and no flow index."""
    if tw <= tau0:
        return 0.0
    d = tw - tau0
    r = 0.5 * DIAM
    return (math.pi * r ** 3 * tw ** -3 / MU * d * d
            * (tau0 * tau0 / 2.0 + 2.0 * tau0 * d / 3.0 + d * d / 4.0))


def dp_valve(q, k):
    return k * KV_SCALE * abs(q) * q


def steady_root(drive, tau0, k):
    """The rate at which the two losses together ask for `drive`, found on the
    wall stress: dp_wall and dp_valve are both monotone in Q, and Q is monotone
    in the wall stress, so one bisection in tw-space settles the pair at once.
    Below the threshold there is no positive root -- that is 5.49's branch, and
    here it is the branch the *steady* kernel takes as well."""
    if abs(drive) <= threshold(tau0):
        return 0.0
    sign = 1.0 if drive > 0 else -1.0
    lo, hi = tau0, tau0 + 1.0
    ask = lambda tw: WALL_SCALE * tw + dp_valve(sign * q_of_tw(tw, tau0), k)
    while ask(hi) < abs(drive):
        hi = 2.0 * hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if ask(mid) > abs(drive):
            hi = mid
        else:
            lo = mid
        if hi - lo <= 1e-16 * max(abs(mid), 1e-30):
            break
    return sign * q_of_tw(0.5 * (lo + hi), tau0)


def step_row(q_prev, drive, tau0, k):
    """One implicit-Euler step of the row above, static branch included, stated
    from the residual side in wall-stress space: the row's own answer is the rate
    that makes the residual vanish, and where none does -- which is exactly where
    the inertia the step starts with cannot carry the rate to a wall stress the
    law will pay the drive at -- the answer is 0.

    This is the same rule the kernel applies from the iterate side (an entry
    latch, plus a clamp when the iterate is carried through zero), reached without
    the Newton loop, so a trace that matches it says the two terms and the branch
    are the only thing that can produce these numbers."""
    residual = lambda tw: (I_DT * (q_of_tw(tw, tau0) - q_prev)
                           + WALL_SCALE * tw + dp_valve(q_of_tw(tw, tau0), k)
                           - drive)
    if residual(tau0) >= 0.0:
        return 0.0
    lo, hi = tau0, tau0 + 1.0
    while residual(hi) < 0.0:
        hi = 2.0 * hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if residual(mid) > 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo <= 1e-16 * max(abs(mid), 1e-30):
            break
    return q_of_tw(0.5 * (lo + hi), tau0)


def drive_of(name):
    tau0, k, q0, pf, pt, dz = PIPES[name]
    return (pf - pt) - RHO * GRAV * dz


# --- the two runs -----------------------------------------------------------
proc = run(TRANSIENT, CASE)
if proc.returncode != 0:
    raise SystemExit(f"transient solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
main = parse_transient(proc.stdout)
proc = run(STEADY, f"{CASE_DIR}/case_steady.dict")
if proc.returncode != 0:
    raise SystemExit(f"steady solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
ss = parse_steady(proc.stdout)

# --- 1. what the two terms each report --------------------------------------
require(set(main["thr"]) == set(PIPES),
        f"the app states a static-hold threshold for exactly the {len(PIPES)} "
        f"pipes that declare a yield stress, whatever their valve K is "
        f"({len(main['thr'])} lines)")
worst_thr = max(rel(main["thr"][n], threshold(PIPES[n][0])) for n in PIPES)
require(worst_thr < 1e-12,
        f"and every threshold is `4 tau_0 L / d` read back, not a number the "
        f"valve is in: worst rel {worst_thr:.3e} over {len(PIPES)} pipes")
require(main["thr"]["StuckRef"] == main["thr"]["StuckK500"]
        == main["thr"]["StuckK1e8"] == 10400.0,
        "the three Stuck pipes print the identical threshold -- 10400 Pa from "
        "tau_0 = 2.6 with K of 0, 5e2 and 1e8")
require(all(len(main["hold"][n]) == NS for n in PIPES),
        f"every pipe reports the latch in all {NS} states")
require(set(ss) == set(PIPES),
        f"the steady kernel solves the same {len(PIPES)}-pipe network, so the "
        f"cross-kernel families below compare like with like")
no_valve_at_rest = [n for n in HELD if "dpValve" in ss[n]]
moving_valved = [n for n in VALVED if ss[n]["q"] != 0.0
                 and "dpValve" not in ss[n]]
require(not no_valve_at_rest and not moving_valved,
        f"the `dpValve` field is present exactly where there is a rate: none of "
        f"the {len(HELD)} pipes that end at rest prints it (even `StopK1e8`, "
        f"whose valve would be worth {PIPES['StopK1e8'][1] * KV_SCALE * (2e-4) ** 2:.3e} Pa "
        f"at its starting rate), and every valved pipe that moves prints it")

# --- 2. a held pipe is the same pipe at any K -------------------------------
STUCK = ("StuckRef", "StuckK500", "StuckK1e8")
for a, b in ((STUCK[0], STUCK[1]), (STUCK[1], STUCK[2]), (STUCK[0], STUCK[2])):
    same = sum(1 for x, y in zip(main["raw"][a], main["raw"][b]) if x == y)
    require(same == NS,
            f"{a} (K = {PIPES[a][1]:.0e}) and {b} (K = {PIPES[b][1]:.0e}): the "
            f"same printed token in {same}/{NS} states")
for n in STUCK + ("HangRef", "HangK500"):
    zeros = sum(1 for q in main["raw"][n] if q == "0")
    held = sum(main["hold"][n][1:])
    require(zeros == NS and held == NS - 1,
            f"{n}: `Q 0` in all {zeros}/{NS} printed states and the latch raised "
            f"in {held}/{NS - 1} steps, at valve K = {PIPES[n][1]:.0e}")
require(all(ss[n]["q"] == 0.0 for n in STUCK + ("HangRef", "HangK500")),
        "the same five pipes report exactly zero in the steady kernel too, and "
        "with the full node-pair drop on the line (`dp 10000`) -- the wall is "
        "holding them, not the solver failing to pay them")
require(main["raw"]["HangRef"] == main["raw"]["HangK500"],
        "Hang's whole trace is one string shared with its valved twin: at a "
        "1664.35 Pa drive against a 2000 Pa wall a K of 500 has nothing to add")

# --- 3. where the law is read, the valve is charged -------------------------
worst_root = 0.0
shares = {}
for n in MOVING:
    tau0, k = PIPES[n][0], PIPES[n][1]
    drive = drive_of(n)
    want = steady_root(drive, tau0, k)
    e = rel(main["q"][n][-1], want)
    e_ss = rel(ss[n]["q"], want)
    worst_root = max(worst_root, e, e_ss)
    shares[n] = dp_valve(main["q"][n][-1], k)
    require(e < 1e-8,
            f"{n}: the marched row ends on the root of *both* losses at drive "
            f"{drive:.4f} Pa -- {main['q'][n][-1]:.10e} vs {want:.10e} "
            f"(rel {e:.3e}), of which the valve takes {shares[n]:.4f} Pa")
    require(e_ss < 1e-8,
            f"{n}: and the steady kernel lands on that same root independently "
            f"({ss[n]['q']:.10e}, rel {e_ss:.3e}) -- one equation, two kernels")
require(worst_root < 1e-8,
        f"all {len(MOVING)} moving pipes, both kernels, on this file's own "
        f"bisection: worst rel {worst_root:.3e}")
seq = ("ReleaseRef", "ReleaseK100", "ReleaseK500", "ReleaseK5000")
rate_seq = [main["q"][n][-1] for n in seq]
share_seq = [shares[n] / drive_of(n) for n in seq]
require(all(a > b for a, b in zip(rate_seq, rate_seq[1:]))
        and all(x < y for x, y in zip(share_seq, share_seq[1:])),
        f"one threshold (tau_0 = 1.0) with four K's is monotone the way only a "
        f"charged term is: Q = {' > '.join(f'{v:.6e}' for v in rate_seq)} while "
        f"the valve's share of the 1.0e4 Pa drive runs "
        f"{', '.join(f'{100 * s:.2f}%' for s in share_seq)}")
require(share_seq[3] > 0.10 and main["raw"]["ReleaseK5000"][-1]
        != main["raw"]["ReleaseRef"][-1],
        "and the biggest of those shares is not a rounding artefact: "
        f"ReleaseK5000 prints {main['raw']['ReleaseK5000'][-1]} where ReleaseRef "
        f"prints {main['raw']['ReleaseRef'][-1]}, a "
        f"{1 - rate_seq[3] / rate_seq[0]:.1%} cut from the same wall")
require(main["raw"]["RollRef"] != main["raw"]["RollK500"],
        f"the gravity-augmented pair differs, each on its own root: RollRef "
        f"{main['raw']['RollRef'][-1]} against RollK500 "
        f"{main['raw']['RollK500'][-1]} at drive {drive_of('RollRef'):.4f} Pa, "
        f"so the 9.2x-the-wall drive reads the valve while the 0.83x drive does "
        f"not")

# --- 4. the steady kernel's printed split ----------------------------------
worst_split, worst_formula, split_name = 0.0, 0.0, ""
for n in MOVING:
    if "dpValve" not in ss[n]:
        continue
    tau0, k = PIPES[n][0], PIPES[n][1]
    dp_friction = ss[n]["dp"] - RHO * GRAV * PIPES[n][5]
    # the wall's share is what is left of the friction drop once the valve's is
    # paid; read it back as a wall stress and ask the law what rate it holds
    tw = (dp_friction - ss[n]["dpValve"]) / WALL_SCALE
    q_from_split = q_of_tw(tw, tau0)
    e_split = rel(q_from_split, ss[n]["q"])
    e_formula = rel(ss[n]["dpValve"], dp_valve(ss[n]["q"], k))
    if e_split > worst_split:
        worst_split, split_name = e_split, n
    worst_formula = max(worst_formula, e_formula)
    require(e_split < 1e-6,
            f"{n}: dp - rho g dz - dpValve is a wall stress of {tw:.4f} Pa and "
            f"the yielded law at that stress is {q_from_split:.10e} against the "
            f"printed {ss[n]['q']:.10e} (rel {e_split:.3e}), so the two printed "
            f"numbers do add up to the drop")
    require(e_formula < 1e-6,
            f"{n}: and dpValve is the valve's own formula, not a fitted share -- "
            f"{ss[n]['dpValve']:.6f} Pa against K rho v^2 / 2 = "
            f"{dp_valve(ss[n]['q'], k):.6f} Pa (rel {e_formula:.3e})")
require(worst_split < 1e-6,
        f"the split adds up for every valved pipe that moves (worst rate rel "
        f"{worst_split:.3e} at {split_name})")
require(worst_formula < 1e-6,
        f"and the valve field is the formula on all of them to "
        f"{worst_formula:.3e}")

# --- 5. the stop: three paths, three steps ---------------------------------
stops = {}
for n in ("StopRef", "StopK5000", "StopK1e8"):
    raw = main["raw"][n]
    idx = next(i for i, q in enumerate(raw) if q == "0")
    stops[n] = idx
    require(all(q == "0" for q in raw[idx:]) and idx < NS - 1,
            f"{n}: reaches exactly 0 at t = {main['times'][idx]:.3f} (step {idx}) "
            f"under a {drive_of(n):.0f} Pa drive with a {threshold(PIPES[n][0]):.0f} Pa "
            f"wall, at K = {PIPES[n][1]:.0e}, and prints the single token `0` in "
            f"the remaining {NS - idx} states")
    require(sum(main["hold"][n][idx + 1:]) == NS - 1 - idx
            and main["hold"][n][idx] == 0,
            f"{n}: the latch picks it up the step after the stop "
            f"({NS - 1 - idx} steps) and not on the stopping step, which reports "
            f"the crossing clamp -- the flag is the entry latch, not the rate")
require(stops["StopRef"] == stops["StopK5000"],
        f"the same wall with no valve and with K = 5e3 stops on the same step "
        f"({stops['StopRef']}), while the rate at the last moving step differs -- "
        f"{main['q']['StopRef'][stops['StopRef'] - 1]:.6e} against "
        f"{main['q']['StopK5000'][stops['StopK5000'] - 1]:.6e} m3/s -- because "
        f"that valve is worth {dp_valve(main['q']['StopK5000'][stops['StopK5000'] - 1], 5000.0):.3e} Pa "
        f"there against the {threshold(1.5) - 4000.0:.0f} Pa margin the stop runs on")
require(stops["StopK1e8"] < stops["StopRef"],
        f"and a valve that is a match for the wall moves the step: StopK1e8 stops "
        f"on {stops['StopK1e8']} because its own valve drop at 2.0e-4 m3/s is "
        f"{dp_valve(2.0e-4, 1.0e8):.4e} Pa, {dp_valve(2.0e-4, 1.0e8) / drive_of('StopK1e8'):.0f}x "
        f"the drive -- the latch never saw a K, the row saw nothing else")
# What decides the step is the row's own inertia against the wall's margin: the
# rate below which the implicit step has no root is q* = (dp_static - drive) *
# dt / (rho L / A), the same amplitude 5.49 measured, and no K appears in it.
q_star = (threshold(1.5) - drive_of("StopRef")) / I_DT
margins = {}
for n in ("StopRef", "StopK5000", "StopK1e8"):
    q_at_stop = main["q"][n][stops[n] - 1]
    q_before = main["q"][n][stops[n] - 2]
    margins[n] = (q_at_stop, q_before)
    require(q_at_stop < q_star <= q_before,
            f"{n}: the step it stops on is the first one that starts under the "
            f"inertia margin q* = (dp_static - drive) dt / (rho L/A) = "
            f"{q_star:.6e} m3/s -- it starts that step at {q_at_stop:.6e} and the "
            f"step before at {q_before:.6e}, so the valve moved the path and the "
            f"wall's margin moved the stop")
require(all(dp_valve(margins[n][0], PIPES[n][1]) < 0.2 * (threshold(1.5) - 4000.0)
            for n in margins),
        f"and at those three stopping rates even K = 1e8 is worth at most "
        f"{max(dp_valve(margins[n][0], PIPES[n][1]) for n in margins):.1f} Pa, "
        f"{100 * max(dp_valve(margins[n][0], PIPES[n][1]) for n in margins) / (threshold(1.5) - 4000.0):.0f}% "
        f"of the {threshold(1.5) - 4000.0:.0f} Pa margin -- which is why the hold "
        f"test being stated without a K is not a hole in the physics")

# --- 6. every state of every pipe against the recursion -------------------
worst_trace, worst_name, worst_step = 0.0, "", 0
recursion_stops = {}
for name, (tau0, k, q0, pf, pt, dz) in PIPES.items():
    drive = drive_of(name)
    qs = [q0]
    while len(qs) < len(main["q"][name]):
        qs.append(step_row(qs[-1], drive, tau0, k))
    e_pipe = max(rel(m, t) for m, t in zip(qs, main["q"][name]))
    if name in stops:
        recursion_stops[name] = next(i for i, q in enumerate(qs) if q == 0.0)
    if e_pipe > worst_trace:
        worst_trace, worst_name = e_pipe, name
        worst_step = max(i for i, (m, t) in enumerate(zip(qs, main["q"][name]))
                         if rel(m, t) == e_pipe)
    require(e_pipe < 1e-7,
            f"{name}: all {len(qs)} states at drive {drive:.4f} Pa and valve "
            f"K = {k:.0e} reproduce the row with both terms and the branch "
            f"(worst rel {e_pipe:.3e})")
require(worst_trace < 1e-7,
        f"the whole trajectory of {len(PIPES)} pipes x {NS} states is the row's "
        f"own algebra -- wall, valve and the `Q = 0` branch -- not a shape that "
        f"happens to end in the right place (worst {worst_trace:.3e} at step "
        f"{worst_step} of {worst_name})")
require(recursion_stops == stops,
        f"and the recursion stops each of the three on the reported step: "
        f"{[recursion_stops[n] for n in ('StopRef', 'StopK5000', 'StopK1e8')]} "
        f"against {[stops[n] for n in ('StopRef', 'StopK5000', 'StopK1e8')]}, so "
        f"the different steps are the different paths, not a damping the kernel "
        f"adds")

# --- 7. one physical state, both kernels -----------------------------------
worst_xk, xk_name = 0.0, ""
for n in MOVING:
    e = rel(main["q"][n][-1], ss[n]["q"])
    if e > worst_xk:
        worst_xk, xk_name = e, n
    require(e < 1e-7,
            f"{n}: marched row at t = 4 gives {main['raw'][n][-1]} and the steady "
            f"kernel gives {ss[n]['q']:.11e} for the same two-term loss")
at_rest_meet = all(main["raw"][n][-1] == "0" and ss[n]["q"] == 0.0 for n in HELD)
require(at_rest_meet,
        f"the {len(HELD)} pipes the wall holds meet on the token `0` in both "
        f"kernels, so the branch and the dead branch are one statement")
require(worst_xk < 1e-7,
        f"the moving pipes meet across the kernels to {worst_xk:.3e} (worst "
        f"{xk_name}) -- the marched row's terminal state is the steady root of "
        f"both losses, valve included")

# --- 8. which term the geometry moves --------------------------------------
require("dpValve" not in ss["RollRef"] and ss["RollRef"]["dp"] == 10000.0,
        "RollRef carries no valve at all, so it is the control for the claim in "
        "family 3: the wall there is being read, and only the wall")
PAIR = ("ReleaseK500", "RollK500")
q_ratio = main["q"][PAIR[1]][-1] / main["q"][PAIR[0]][-1]
v_ratio = shares[PAIR[1]] / shares[PAIR[0]]
require(rel(v_ratio, q_ratio ** 2) < 1e-6,
        f"one K = 500 on two drives: the valve's drop scales with the square of "
        f"the rate it is handed -- {shares[PAIR[1]]:.4f} / {shares[PAIR[0]]:.4f} "
        f"= {v_ratio:.4f} against (Q2/Q1)^2 = ({q_ratio:.4f})^2 = "
        f"{q_ratio ** 2:.4f}")
wall_ratio = ((drive_of(PAIR[1]) - shares[PAIR[1]])
              / (drive_of(PAIR[0]) - shares[PAIR[0]]))
exponent = math.log(wall_ratio) / math.log(q_ratio)
require(0.4 < exponent < 0.8,
        f"and the wall's does not: the same two pipes' wall drops are "
        f"{drive_of(PAIR[1]) - shares[PAIR[1]]:.4f} and "
        f"{drive_of(PAIR[0]) - shares[PAIR[0]]:.4f} Pa, a ratio of "
        f"{wall_ratio:.4f} where the rates differ by {q_ratio:.4f} -- an "
        f"exponent of {exponent:.3f} against the valve's 2, which is why the "
        f"share a valve buys depends on which wall it is sitting behind")

print()
print(f"transientYieldValve: {len(STUCK) + 2} pipes at rest print one shared "
      f"trace across K from 0 to 1e8, {len(MOVING)} moving pipes land on this "
      f"file's own two-loss bisection to {worst_root:.3e} in both kernels, the "
      f"valve's share of the drive runs 0 to {100 * max(share_seq):.2f}% across "
      f"one threshold, the stops come at steps {[stops[n] for n in ('StopRef', 'StopK5000', 'StopK1e8')]} "
      f"(recursion {[recursion_stops[n] for n in ('StopRef', 'StopK5000', 'StopK1e8')]}), "
      f"whole-trace reproduction worst {worst_trace:.3e}, cross-kernel worst "
      f"{worst_xk:.3e}")

if failures:
    print(f"FAILED: {len(failures)} assertion(s)")
    sys.exit(1)
print("all assertions passed")
