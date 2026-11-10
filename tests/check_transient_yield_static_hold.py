#!/usr/bin/env python3
"""transientYieldStaticHold: the yield threshold has a branch of its own in the
transient momentum row, checked by re-deriving every reported step and by
reading the same physical state off the steady kernel.

A yielded wall does not stop holding when the rate stops: dp_wall(Q) -> dp_static
= 4 tau_0 L / d as Q -> 0+, and *below* that number the moving law has no root at
all. The steady kernel sees the asked drop, so it can decide that before
bracketing (5.46, 5.47). A marching row solves for the drop, so it needs the
branch: while the fluid is at rest and the drive it is offered does not reach the
threshold, the row *is* `Q = 0`.

For a pipe between two fixed-pressure nodes that row is one equation,

    rho L/A * (Q_new - Q_prev)/dt = drive - dp_wall(Q_new),
    drive = pFrom - pTo - rho g dz

and dp_wall is the shared kernel function (pipeFrictionPa) made odd in Q. So the
checker does not compare this run against remembered numbers: it solves the
recursion itself, with the law and the branch written out from the case file, and
walks every reported state. Where the branch is the whole statement it says so on
the printed string instead of on a tolerance -- `Q 0` in all 2001 samples is not a
small number, it is the row.

Ten families of judgement:
  1. the threshold is reported, per pipe, as `4 tau_0 L / d` and nowhere else
  2. the static branch: exactly zero, in every state, for four different reasons
  3. gravity is inside the drive the wall resists (one threshold, both sides of
     it, same pipe, same node pair)
  4. released from rest, the row lands on 5.46's and 5.47's closed forms, and
     n = 1 with K = mu is bit-identical to declaring the threshold alone
  5. the finite-time stop, and what the moving law alone could have printed
  6. every reported state of every pipe against the recursion *with* the branch
  7. a held pipe seals the node it feeds: the balance left over is the solver's
     own exit residual, not a leak through the wall
  8. the network itself can let go: the latch flips at the step whose entry drive
     crosses the threshold, which is the printed node pressure of the previous
     sample
  9. the terminal states of the released pair, and of the two stuck pairs, read
     off the *steady* kernel for the same topology
 10. the two declaration-face refusals this branch still refuses
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

# name: (diameter, length, n, K, tau_0, initialFlow, pFrom, pTo, dz)
PIPES = {
    "Stuck":        (0.02, 20.0, 0.0, 0.0, 2.6,  0.0,     P_A,  P_B,  0.0),
    "StuckHB":      (0.02, 20.0, 0.7, 2.0, 2.6,  0.0,     P_A,  P_B,  0.0),
    "ReleaseB05":   (0.02, 20.0, 0.0, 0.0, 0.5,  0.0,     P_A,  P_B,  0.0),
    "ReleaseB10":   (0.02, 20.0, 0.0, 0.0, 1.0,  0.0,     P_A,  P_B,  0.0),
    "ReleaseN1":    (0.02, 20.0, 1.0, 0.1, 1.0,  0.0,     P_A,  P_B,  0.0),
    "ReleaseHB05":  (0.02, 20.0, 0.5, 1.0, 0.5,  0.0,     P_A,  P_B,  0.0),
    "Hang":         (0.02, 20.0, 0.0, 0.0, 0.5,  0.0,     P_A,  P_B,  1.0),
    "Roll":         (0.02, 20.0, 0.0, 0.0, 0.5,  0.0,     P_A,  P_B, -1.0),
    "HeldBelow":    (0.02, 20.0, 0.0, 0.0, 1.5,  0.0,     P_A2, P_B,  0.0),
    "StopBelow":    (0.02, 20.0, 0.0, 0.0, 1.5,  2.0e-6,  P_A2, P_B,  0.0),
}
INERTIA = RHO * 20.0 / (0.25 * math.pi * 0.02 ** 2)   # rho L / A, every pipe
I_DT = INERTIA / DT
# the Newtonian resistance of the moving half of every series pair
R_NEWTON = 128.0 * MU * 20.0 / (math.pi * 0.02 ** 4)

# case_series.dict: one stuck pipe and one moving pipe on a free node, twice
# mirrored, plus the pair the network lets go of.
SERIES = {
    "HoldUp":      (0.02, 20.0, 0.0, 0.0, 2.6),
    "BleedDown":   (0.02, 20.0, 0.0, 0.0, 0.0),
    "HoldDown":    (0.02, 20.0, 0.0, 0.0, 2.6),
    "FeedUp":      (0.02, 20.0, 0.0, 0.0, 0.0),
    "ReleaseLate": (0.02, 20.0, 0.0, 0.0, 2.0),
    "Bleed":       (0.02, 20.0, 0.0, 0.0, 0.0),
}
NODE_VOLUME = {"C": 1.0, "D": 1.0, "F": 0.1}
NODE_BULK = 1.0e8
NODE_P0 = {"A": P_A, "B": P_B, "C": 2.5e4, "D": 2.1e4, "F": 2.5e4}

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


def parse_transient(stdout, names):
    """Q trace and static-hold latch per pipe, node pressure traces, thresholds."""
    out = {"q": {n: [] for n in names}, "raw": {n: [] for n in names},
           "hold": {n: [] for n in names}, "p": {n: [] for n in NODE_P0},
           "thr": {}, "times": []}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 2 and tok[0] == "t":
            out["times"].append(float(tok[1]))
        elif len(tok) >= 4 and tok[0] == "node" and tok[2] == "p":
            if tok[1] in out["p"]:
                out["p"][tok[1]].append(float(tok[3]))
        elif len(tok) == 4 and tok[0] == "pipe" and tok[2] == "Q":
            if tok[1] in out["q"]:
                out["q"][tok[1]].append(float(tok[3]))
                out["raw"][tok[1]].append(tok[3])
        elif len(tok) == 5 and tok[0] == "pipe" and tok[3] == "hold":
            if tok[1] in out["hold"]:
                out["hold"][tok[1]].append(int(tok[4]))
        elif len(tok) == 5 and tok[0] == "pipe" and tok[3] == "threshold":
            out["thr"][tok[1]] = float(tok[4])
    return out


def parse_steady(stdout, names):
    """One state per pipe: the steady app prints Q in a longer line."""
    q, p = {}, {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 4 and tok[0] == "pipe" and tok[2] == "Q" and tok[1] in names:
            q[tok[1]] = float(tok[3])
        elif len(tok) >= 4 and tok[0] == "node" and tok[2] == "p":
            p[tok[1]] = float(tok[3])
    return q, p


# --- the wall law and the branch, written out from the case file -------------

def threshold(tau0, length, diam):
    return 4.0 * tau0 * length / diam


def q_of_tw(tw, diam, length, n, k, tau0):
    """The exact laminar discharge a yielded wall stress holds: 5.47's profile
    integral, which is 5.46's Buckingham bracket at n = 1, k = mu."""
    if tw <= tau0:
        return 0.0
    nn = n if n > 0 else 1.0
    kk = k if n > 0 else MU
    d = tw - tau0
    r = 0.5 * diam
    return (math.pi * r ** 3 * tw ** -3 * kk ** (-1.0 / nn) * nn
            * d ** (1.0 + 1.0 / nn)
            * (tau0 * tau0 / (nn + 1.0) + 2.0 * tau0 * d / (2.0 * nn + 1.0)
               + d * d / (3.0 * nn + 1.0)))


def dp_wall_mag(q, diam, length, n, k, tau0):
    """The moving law's drop at a rate magnitude, by inverting the profile
    integral for the wall stress: monotone, so a bisection lands on it."""
    if q == 0.0:
        return threshold(tau0, length, diam)
    lo, hi = tau0, tau0 + 1e4
    for _ in range(300):
        if q_of_tw(hi, diam, length, n, k, tau0) > q:
            break
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if q_of_tw(mid, diam, length, n, k, tau0) > q:
            hi = mid
        else:
            lo = mid
        if hi - lo <= 1e-15 * max(abs(mid), 1e-30):
            break
    return 4.0 * length * 0.5 * (lo + hi) / diam


def dp_wall(q, diam, length, n, k, tau0):
    """Its odd continuation -- and the jump at Q = 0 is the whole difficulty:
    dp_wall(0+) is the threshold, not nothing."""
    if q == 0.0:
        return 0.0
    dp = dp_wall_mag(abs(q), diam, length, n, k, tau0)
    return dp if q > 0 else -dp


def step_row(q_prev, drive, diam, length, n, k, tau0):
    """One implicit-Euler step of a single row, branch included, stated from the
    residual side: solve for the rate the row asks for, and where the moving law
    cannot pay the drive at all -- which happens exactly when |drive| does not
    reach the threshold -- the row has no root and its answer is 0.

    This is the same rule the kernel applies from the iterate side (a latch at
    step entry, plus a clamp when the Newton iterate is carried through zero),
    reached without the Newton loop, so a trace that matches it says the branch
    is the only thing that can produce these numbers."""
    eps = 1e-14
    scale = max(abs(q_prev), abs(drive) / I_DT, 1e-6)
    for lo, hi in ((eps, 4.0 * scale), (-4.0 * scale, -eps)):
        a, b = lo, hi
        fa = I_DT * (a - q_prev) + dp_wall(a, diam, length, n, k, tau0) - drive
        fb = I_DT * (b - q_prev) + dp_wall(b, diam, length, n, k, tau0) - drive
        if fa * fb > 0.0:
            continue
        for _ in range(140):
            m = 0.5 * (a + b)
            fm = I_DT * (m - q_prev) + dp_wall(m, diam, length, n, k, tau0) - drive
            if fa * fm <= 0.0:
                b, fb = m, fm
            else:
                a, fa = m, fm
            if b - a <= 1e-16 * max(abs(m), 1e-30):
                break
        return 0.5 * (a + b)
    return 0.0


# --- case.dict: the ten single-row pipes ------------------------------------
proc = run(TRANSIENT, CASE)
if proc.returncode != 0:
    raise SystemExit(f"solver failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")
main = parse_transient(proc.stdout, PIPES)

# --- 1. the threshold is reported, per pipe, and nowhere else ---------------
require(set(main["thr"]) == set(PIPES),
        f"the app states a static-hold threshold for exactly the {len(PIPES)} "
        f"pipes that declare a yield stress, and no line for a pipe that does "
        f"not ({sorted(main['thr'])})")
require(all(len(main["hold"][k]) == NS for k in PIPES),
        f"and every one of them reports the latch in all {NS} states")
worst_thr = max(rel(main["thr"][k], threshold(PIPES[k][4], PIPES[k][1], PIPES[k][0]))
                for k in PIPES)
require(worst_thr < 1e-12 and main["thr"]["Stuck"] == main["thr"]["StuckHB"],
        f"the reported threshold is `4 tau_0 L / d` read back, not the fluid's "
        f"number: worst rel {worst_thr:.3e}, and the same "
        f"{main['thr']['Stuck']:.6g} Pa comes out of the undeclared `Stuck` and "
        f"of `StuckHB` with n = 0.7 in front of it")

# --- 2. the static branch, four different reasons ---------------------------
HELD = {
    "Stuck": "the wall's 1.04e4 Pa is more than the 1.0e4 Pa on offer",
    "StuckHB": "and the threshold does not know the flow index, so the declared "
               "n = 0.7 changes nothing about it",
    "Hang": "the 1.0e4 Pa node pair is 1664.35 Pa once the pipe's own 1 m rise "
            "is paid for, against a 2000 Pa wall",
    "HeldBelow": "at rest under a 4000 Pa drive with a 6000 Pa wall, there is no "
                 "rate the moving law can ask for",
}
for name, why in HELD.items():
    zeros = sum(1 for q in main["raw"][name] if q == "0")
    held = sum(main["hold"][name][1:])
    require(zeros == NS and held == NS - 1,
            f"{name}: `Q 0` in all {zeros}/{NS} printed states and the latch "
            f"raised in {held}/{NS - 1} steps -- {why}")

# --- 3. gravity belongs to the drive the wall resists -----------------------
drive_hang = (P_A - P_B) - RHO * GRAV * 1.0
drive_roll = (P_A - P_B) + RHO * GRAV * 1.0
require(main["thr"]["Hang"] == main["thr"]["Roll"] == 2000.0,
        "Hang and Roll are one pipe with one threshold (2000 Pa) and one node "
        "pair, differing only in the sign of dz")
require(rel(drive_hang, 1664.3475) < 1e-6 and rel(drive_roll, 18335.6525) < 1e-6,
        f"so the drive is {drive_hang:.4f} Pa hanging up and {drive_roll:.4f} Pa "
        f"rolling down, and rho g dz = {RHO * GRAV:.4f} Pa is what separates "
        f"them: below the threshold on one side, 9.2x above it on the other")
q_roll_ss = q_of_tw(drive_roll * 0.02 / (4.0 * 20.0), 0.02, 20.0, 0.0, 0.0, 0.5)
e_roll = rel(main["q"]["Roll"][-1], q_roll_ss)
require(e_roll < 1e-8,
        f"Roll relaxes onto Buckingham at the *gravity-augmented* drop: "
        f"{main['q']['Roll'][-1]:.10e} vs {q_roll_ss:.10e} (rel {e_roll:.3e}) -- "
        f"tw = {drive_roll * 0.02 / 80.0:.4f} Pa, not the 2.5 Pa the node pair "
        f"alone would give")
require(max(abs(q) for q in main["q"]["Hang"]) == 0.0,
        "and a branch that compared the node pressures alone would have moved "
        "Hang: it does not, because what the wall resists is the head minus the "
        "rise")

# --- 4. released from rest lands on the steady closed forms -----------------
released = ("ReleaseB05", "ReleaseB10", "ReleaseN1", "ReleaseHB05")
worst_rel = 0.0
for name in released:
    diam, length, n, k, tau0 = PIPES[name][:5]
    tw = (P_A - P_B) * diam / (4.0 * length)
    want = q_of_tw(tw, diam, length, n, k, tau0)
    e = rel(main["q"][name][-1], want)
    worst_rel = max(worst_rel, e)
    require(e < 1e-8,
            f"{name}: the row ends on {'Buckingham' if not n > 0 else 'HB'} at "
            f"tw = {tw:.4f} Pa (xi = {tau0 / tw:.4f}, D = {tw - tau0:.4f} Pa), "
            f"Q = {main['q'][name][-1]:.10e} vs {want:.10e} (rel {e:.3e})")
require(worst_rel < 1e-8,
        f"the same yielded closed forms 5.46 and 5.47 measure in the steady "
        f"network are met by the marched row to {worst_rel:.3e} on all "
        f"{len(released)} released pipes -- one law, two kernels")
same = sum(1 for a, b in zip(main["raw"]["ReleaseB10"], main["raw"]["ReleaseN1"])
           if a == b)
require(same == NS,
        f"n = 1 with K = mu is inert on a yielded row too: {same}/{NS} states of "
        f"ReleaseB10 and ReleaseN1 are the same printed string, so 5.47's "
        f"containment of 5.46 now holds in time as well as in the bracket")

# --- 5. the stop, in finite time --------------------------------------------
stop, stop_raw = main["q"]["StopBelow"], main["raw"]["StopBelow"]
first_zero = next(i for i, q in enumerate(stop) if q == 0.0)
require(all(q == 0.0 for q in stop[first_zero:])
        and set(stop_raw[first_zero:]) == {"0"},
        f"StopBelow starts flowing backwards-free at 2.0e-6 m3/s under a 4000 Pa "
        f"drive with a 6000 Pa wall: it reaches exactly 0 at t = "
        f"{main['times'][first_zero]:.3f} ({first_zero} steps) and prints the "
        f"single token `0` in every one of the remaining {NS - first_zero} states")
require(sum(main["hold"]["StopBelow"][first_zero + 1:]) == NS - 1 - first_zero
        and main["hold"]["StopBelow"][first_zero] == 0,
        f"the latch picks it up the step after the one that stopped it: "
        f"{NS - 1 - first_zero} steps report it, and the stopping step itself "
        f"reports the crossing clamp instead -- the flag is the entry latch, not "
        f"the rate")
naive = (threshold(1.5, 20.0, 0.02) - 4000.0) * DT / INERTIA
require(naive > 0.0 and stop_raw[-1] == "0" and stop[first_zero - 1] < naive,
        f"without the branch this row has nothing to converge to: the odd "
        f"continuation of the moving law leaves the row stepping by "
        f"(dp_static - drive) dt / (rho L/A) = {naive:.3e} m3/s around zero, "
        f"while the last rate before the reported stop is "
        f"{stop[first_zero - 1]:.3e} m3/s -- inside one such step of the token "
        f"`0` the branch prints, so what separates the two readings is a branch "
        f"and not a tolerance")
require(stop[first_zero - 1] > 0.0,
        f"and it stops *through* zero: the last rate before the stop is still "
        f"{stop[first_zero - 1]:.3e} m3/s going forward, and the step that "
        f"follows it is the row `Q = 0` rather than a small reverse flow")

# --- 6. every state of every pipe against the recursion ---------------------
worst_trace, worst_trace_name, worst_step = 0.0, "", 0
recursion_zero = {}
for name, (diam, length, n, k, tau0, q0, pf, pt, dz) in PIPES.items():
    drive = (pf - pt) - RHO * GRAV * dz
    qs = [q0]
    while len(qs) < len(main["q"][name]):
        qs.append(step_row(qs[-1], drive, diam, length, n, k, tau0))
    e_pipe = max(rel(m, t) for m, t in zip(qs, main["q"][name]))
    if name == "StopBelow":
        recursion_zero["stop"] = next(i for i, q in enumerate(qs) if q == 0.0)
    if e_pipe > worst_trace:
        worst_trace, worst_trace_name = e_pipe, name
        worst_step = max(i for i, (m, t) in enumerate(zip(qs, main["q"][name]))
                         if rel(m, t) == e_pipe)
    require(e_pipe < 1e-7,
            f"{name}: all {len(qs)} states at drive {drive:.4f} Pa reproduce the "
            f"row with the branch (worst rel {e_pipe:.3e})")
require(worst_trace < 1e-7,
        f"the whole trajectory of {len(PIPES)} pipes x {NS} states is the row's "
        f"own algebra including `Q = 0`, not a shape that happens to end in the "
        f"right place (worst {worst_trace:.3e} at step {worst_step} of "
        f"{worst_trace_name})")
require(recursion_zero.get("stop") == first_zero,
        f"and the recursion stops that pipe at the reported step: "
        f"{recursion_zero.get('stop')} there against {first_zero} here, so the "
        f"finite-time stop is the branch's arithmetic rather than a damping the "
        f"kernel adds")

# --- case_series.dict -------------------------------------------------------
proc = run(TRANSIENT, f"{CASE_DIR}/case_series.dict")
if proc.returncode != 0:
    raise SystemExit(f"series solver failed\n{proc.stdout}\n{proc.stderr}")
ser = parse_transient(proc.stdout, SERIES)

# --- 7. a held pipe seals the node it feeds ---------------------------------
# The node row the kernel assembles is
#     volume/(bulk*dt) * (p_i - p_{i-1})  -  (Q_in - Q_out) = 0
# and that row's own exit scale is |stor| + |sumIn| + 1.0 (transient.cpp), which
# at these throughputs is 1 m3/s -- so the solver is allowed to leave an
# *absolute* residue of tol*1 = 1e-11 m3/s in it. That is the floor this family
# measures against: not machine zero, and not a leak either, because the held
# pipe enters the balance as the token `0`.
NEWTON_TOL = 1.0e-11
worst_seal, seal_name, worst_leak, thru = 0.0, "", 0.0, 0.0
for node, inflow, outflow in (("C", "HoldUp", "BleedDown"),
                              ("D", "FeedUp", "HoldDown")):
    stuck = inflow if inflow.startswith("Hold") else outflow
    moving = outflow if stuck == inflow else inflow
    require(ser["raw"][stuck] == ["0"] * NS
            and sum(ser["hold"][stuck][1:]) == NS - 1,
            f"{stuck}: exactly zero and latched in every state while its far end "
            f"is a node pressure the solve has to find (node {node})")
    compliance = NODE_VOLUME[node] / NODE_BULK
    leaks, scaled = [], []
    for i in range(1, NS):
        stor = compliance * (ser["p"][node][i] - ser["p"][node][i - 1]) / DT
        sum_in = ser["q"][inflow][i] - ser["q"][outflow][i]
        leaks.append(abs(stor - sum_in))
        scaled.append(abs(stor - sum_in) / (abs(stor) + abs(sum_in) + 1.0))
    e_scaled, e_leak = max(scaled), max(leaks)
    thru = max(thru, max(abs(v) for v in ser["q"][moving]))
    if e_scaled > worst_seal:
        worst_seal, seal_name = e_scaled, node
    worst_leak = max(worst_leak, e_leak)
    require(e_scaled < 5.0 * NEWTON_TOL,
            f"node {node}: the balance the kernel writes closes to its own exit "
            f"criterion, {e_scaled:.3e} against the {NEWTON_TOL:.0e} it stops on "
            f"(worst absolute leak {e_leak:.3e} m3/s), so the seal is the "
            f"solver's residual and not a wall letting fluid through")
require(worst_seal < 5.0 * NEWTON_TOL,
        f"both floating nodes are sealed to the exit residual, not to a "
        f"tolerance chosen for the check (worst {worst_seal:.3e} at {seal_name})")
require(worst_leak / thru < 1e-5,
        f"and that residual is small enough to be a solver number, not a "
        f"physical one: {worst_leak:.3e} m3/s is {worst_leak / thru:.3e} of the "
        f"largest rate any node actually moves ({thru:.3e}) -- the +1.0 in the "
        f"row's scale showing through -- and as a leak it would need "
        f"{NODE_VOLUME['C'] * 5000.0 / NODE_BULK / worst_leak:.3e} s to drain "
        f"node C the 5000 Pa its own solve covers in {END:.1f} s")

# --- 8. the network can let go ----------------------------------------------
hold_rl = ser["hold"]["ReleaseLate"]
rel_idx = next(i for i in range(1, NS) if hold_rl[i] == 0)
drive_entry = P_A - ser["p"]["F"][rel_idx - 1]
drive_before = P_A - ser["p"]["F"][rel_idx - 2]
require(all(hold_rl[j] == 1 for j in range(1, rel_idx))
        and ser["raw"]["ReleaseLate"][:rel_idx] == ["0"] * rel_idx,
        f"ReleaseLate (tau_0 = 2.0, threshold 8000 Pa) starts held at a 5000 Pa "
        f"drive, prints exactly 0 in all {rel_idx} states before it moves, and is "
        f"released at t = {ser['times'][rel_idx]:.3f}")
require(drive_entry > 8000.0 and drive_before <= 8000.0,
        f"the latch flips on the step whose *entry* drive crosses the threshold: "
        f"node F had fallen to {ser['p']['F'][rel_idx - 1]:.4f} Pa, which leaves "
        f"{drive_entry:.4f} Pa against 8000 Pa, and one sample earlier it was "
        f"{drive_before:.4f} Pa -- so being held is a state of the row, not a "
        f"property of the pipe")
require(all(h == 0 for h in hold_rl[rel_idx:]),
        f"and it does not go back: no later step re-latches, because the node it "
        f"drains ends at {ser['p']['F'][-1]:.4f} Pa, a drive of "
        f"{P_A - ser['p']['F'][-1]:.4f} Pa its own wall must then resist")

# --- 9. one physical state, both kernels ------------------------------------
proc = run(STEADY, f"{CASE_DIR}/case_steadySeries.dict")
if proc.returncode != 0:
    raise SystemExit(f"steady reference failed\n{proc.stdout}\n{proc.stderr}")
ss, ss_p = parse_steady(proc.stdout, SERIES)
e_q = rel(ser["q"]["ReleaseLate"][-1], ss["ReleaseLate"])
e_p = rel(ser["p"]["F"][-1], ss_p["F"])
require(e_q < 1e-6,
        f"the pair the network let go of ends on the rate the *steady* kernel "
        f"computes for the same topology: transient "
        f"{ser['q']['ReleaseLate'][-1]:.10e} vs steady {ss['ReleaseLate']:.10e} "
        f"(rel {e_q:.3e}), and its node lands on {ser['p']['F'][-1]:.6f} vs "
        f"{ss_p['F']:.6f} Pa (rel {e_p:.3e})")
require(ss["HoldUp"] == 0.0 and ss["HoldDown"] == 0.0,
        "the two stuck pairs report exactly zero in the steady kernel too, which "
        "is the cross-kernel half of family 2")
tau_c = (NODE_VOLUME["C"] / NODE_BULK) * R_NEWTON
require(ser["p"]["C"][-1] > ss_p["C"] and ser["p"]["D"][-1] < ss_p["D"]
        and rel(ser["p"]["C"][-1], ss_p["C"]) > 0.05,
        f"and the two sealed nodes are *not* at their steady pressures yet "
        f"({ser['p']['C'][-1]:.1f} vs {ss_p['C']:.1f} Pa, {ser['p']['D'][-1]:.1f} "
        f"vs {ss_p['D']:.1f} Pa), each on the side it is heading: the compliance "
        f"and the moving pipe's resistance give a relaxation time of "
        f"(volume/bulk) * 128 mu L/(pi d^4) = {tau_c:.2f} s, so 4 s is "
        f"{END / tau_c:.2f} of them -- the run agrees with 5.46 about the rate "
        f"long before it agrees about the floating pressure")

# --- 10. the declaration face still refuses ---------------------------------
REFUSALS = [
    ("case_rejectNegativeYield.dict", "negative yield",
     "a threshold below zero empties the held set and then pays every rate"),
    ("case_rejectCoupledYield.dict", "cannot carry a threshold",
     "a coupled element's loss is its partner's, so no wall is left to hold"),
]
for fname, needle, why in REFUSALS:
    p = run(TRANSIENT, f"{CASE_DIR}/{fname}")
    msg = (p.stdout + p.stderr).replace("[chflow] ", "")
    line = next((ln for ln in msg.splitlines() if needle in ln), "")
    require(p.returncode != 0 and needle in msg,
            f"{fname}: refused, naming {why} -- '{line}'")

print()
print(f"transientYieldStaticHold: {len(HELD)} pipes exactly zero in all {NS} "
      f"states, {len(released)} released pipes onto the yielded closed forms to "
      f"{worst_rel:.3e}, stop at step {first_zero} (recursion "
      f"{recursion_zero.get('stop')}) against an odd-law amplitude of "
      f"{naive:.3e}, whole-trace reproduction worst {worst_trace:.3e}, node seal "
      f"worst {worst_seal:.3e} against the {NEWTON_TOL:.0e} the solver stops on "
      f"({worst_leak:.3e} m3/s absolute), released at step {rel_idx} and ending "
      f"{e_q:.3e} off the steady kernel's own rate")

if failures:
    print(f"FAILED: {len(failures)} assertion(s)")
    sys.exit(1)
print("all assertions passed")
