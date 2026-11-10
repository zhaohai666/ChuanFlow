#!/usr/bin/env python3
"""transientYieldValveSchedule: the valve *timetable* (M2's `valveClose` /
`valveOpen`, architecture `§4.2`) crossed with the yield wall's static branch
(5.49), which is the combination 5.50 put in its own boundary ③.

The row, per pipe, both ends pinned:

    rho L/A * (Q_new - Q_prev)/dt = drive - dp_wall(Q_new) - K(t) rho|Q|Q/(2A^2)
    drive = pFrom - pTo - rho g dz,   K(t) = valveKAt(t + dt)

Two switches live on that row and they are read from different places, so the
question is not "can a big K stop a yielded pipe" but *which switch a schedule
can reach*:

  * the latch (`transient.cpp:40`-`:44`) compares `drive` with `4 tau_0 L/d` and
    `drive` contains no K (`:52`-`:56`);
  * the schedule (`:119`-`:133`, sampled at `t + dt` by `:293`-`:300`) only
    changes the term the *moving* row pays.

Seven families of judgement:
  1. above the threshold a valve is a discount, not a refusal (no pipe is ever
     held, even when the closed valve is worth 9.8x the whole drive)
  2. below it the wall does the stopping and the schedule only sets the step
     (184 / 72 / 54, monotone in the ramp, and a ramp that arrives late is
     invisible to the trace)
  3. a schedule cannot wake what the wall holds: the yielded pipe stays on the
     single token `0` through the whole reopening while the Newtonian twin with
     the identical schedule returns to its own no-valve rate
  4. a pipe already at rest under a sub-threshold drive does not answer a
     schedule at all (string-for-string identity with the same pipe, no valve)
  5. every state of every pipe against this file's own recursion of both laws
     with the sampled K -- the relative gate on the moving states, an absolute
     one over all of them, and the step each pipe stops on
  6. the K really is the one at `t + dt`, not the one at `t`: for a ramp that
     overlaps motion the two conventions disagree, and the kernel's stop step is
     the `t + dt` one; for a ramp that never overlaps motion they are the same
     trace exactly, which is the same invisibility family 4 found
  7. the settled row against the steady kernel with every valve frozen at its
     end-of-schedule value, and no `dpValve` field on a pipe reporting `Q 0`
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
INERTIA = RHO * LENGTH / AREA
I_DT = INERTIA / DT
KV_SCALE = RHO / (2.0 * AREA * AREA)
WALL_SCALE = 4.0 * LENGTH / DIAM
CLOSED_K = 1.0e6

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


# name: (tau_0, K_open, initial Q, pFrom, pTo, closeStart, closeEnd, openStart,
#        openEnd)
PIPES = {
    "CreepOpen":  (1.5, 0.0, 0.0, P_A, P_B, -1.0, -1.0, -1.0, -1.0),
    "CreepCyc":   (1.5, 0.0, 0.0, P_A, P_B, 0.5, 0.65, 1.5, 1.65),
    "CreepShut":  (1.5, 0.0, 0.0, P_A, P_B, 0.5, 0.65, -1.0, -1.0),
    "CreepK1e6":  (1.5, CLOSED_K, 0.0, P_A, P_B, -1.0, -1.0, -1.0, -1.0),
    "StopPlain":  (1.5, 0.0, 2.0e-4, P_A2, P_B, -1.0, -1.0, -1.0, -1.0),
    "StopEarly":  (1.5, 0.0, 2.0e-4, P_A2, P_B, 0.05, 0.15, 1.5, 1.65),
    "StopMid":    (1.5, 0.0, 2.0e-4, P_A2, P_B, 0.05, 0.6, 1.5, 1.65),
    "StopLate":   (1.5, 0.0, 2.0e-4, P_A2, P_B, 0.5, 0.65, 1.5, 1.65),
    "HeldPlain":  (1.5, 0.0, 0.0, P_A2, P_B, -1.0, -1.0, -1.0, -1.0),
    "HeldCyc":    (1.5, 0.0, 0.0, P_A2, P_B, 0.05, 0.15, 1.5, 1.65),
    "NewtPlain":  (0.0, 0.0, 2.0e-4, P_A2, P_B, -1.0, -1.0, -1.0, -1.0),
    "NewtEarly":  (0.0, 0.0, 2.0e-4, P_A2, P_B, 0.05, 0.15, 1.5, 1.65),
}
YIELDED = tuple(n for n in PIPES if PIPES[n][0] > 0)
ABOVE = ("CreepOpen", "CreepCyc", "CreepShut", "CreepK1e6")
STOP = ("StopPlain", "StopEarly", "StopMid", "StopLate")
HELD_REST = ("HeldPlain", "HeldCyc")
SCHEDULED = tuple(n for n in PIPES if PIPES[n][5] >= 0)

# the valve frozen at the value its schedule holds at t = endTime, for the
# steady cross-check
K_END = {n: (CLOSED_K if PIPES[n][7] < 0 and PIPES[n][5] >= 0
             else PIPES[n][1]) for n in PIPES}


def valve_k(name, t):
    """transient.cpp:119-133 verbatim."""
    tau0, kop, q0, pf, pt, cs, ce, os, oe = PIPES[name]
    if cs < 0:
        return kop
    if t < cs:
        return kop
    if t < ce:
        return kop + (t - cs) / (ce - cs) * (CLOSED_K - kop)
    if os < 0 or t < os:
        return CLOSED_K
    if t < oe:
        return CLOSED_K + (t - os) / (oe - os) * (kop - CLOSED_K)
    return kop


def threshold(tau0):
    return WALL_SCALE * tau0


def q_of_tw(tw, tau0):
    if tw <= tau0:
        return 0.0
    d = tw - tau0
    r = 0.5 * DIAM
    return (math.pi * r ** 3 * tw ** -3 / MU * d * d
            * (tau0 * tau0 / 2.0 + 2.0 * tau0 * d / 3.0 + d * d / 4.0))


def dp_valve(q, k):
    return k * KV_SCALE * abs(q) * q


def drive_of(name):
    pf, pt = PIPES[name][3], PIPES[name][4]
    return (pf - pt) - RHO * GRAV * 0.0


def steady_root(drive, tau0, k):
    if abs(drive) <= threshold(tau0):
        return 0.0
    lo, hi = tau0, tau0 + 1.0
    ask = lambda tw: WALL_SCALE * tw + dp_valve(q_of_tw(tw, tau0), k)
    while ask(hi) < abs(drive):
        hi = 2.0 * hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if ask(mid) > abs(drive):
            hi = mid
        else:
            lo = mid
    return q_of_tw(0.5 * (lo + hi), tau0)


def step_row(q_prev, drive, tau0, k):
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
    return q_of_tw(0.5 * (lo + hi), tau0)


def recurse(name, k_at_t=False):
    """The whole trace, from this file's own copy of both laws and the sampled K."""
    tau0 = PIPES[name][0]
    q = PIPES[name][2]
    drive = drive_of(name)
    out = [q]
    for i in range(1, NS):
        t = i * DT
        k = valve_k(name, (i - 1) * DT if k_at_t else t)
        q = step_row(q, drive, tau0, k)
        out.append(q)
    return out


def stop_step(trace):
    """The first state a trace prints 0 on after it was last moving."""
    for i in range(1, len(trace)):
        if trace[i] == 0.0 and trace[i - 1] != 0.0:
            return i
    return None


def parse_transient(stdout):
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
    out = {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) < 8 or tok[0] != "pipe" or tok[2] != "Q" or tok[1] not in PIPES:
            continue
        rec = {"q": float(tok[3]), "raw": tok}
        for i, t in enumerate(tok):
            if t == "dp":
                rec["dp"] = float(tok[i + 1])
            elif t == "dpValve":
                rec["dpValve"] = float(tok[i + 1])
        out[tok[1]] = rec
    return out


proc = run(TRANSIENT, CASE)
if proc.returncode != 0:
    raise SystemExit(f"transient solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
main = parse_transient(proc.stdout)
proc = run(STEADY, f"{CASE_DIR}/case_steady.dict")
if proc.returncode != 0:
    raise SystemExit(f"steady solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
steady = parse_steady(proc.stdout)

require(len(main["times"]) == NS, f"{NS} states printed")
for n in PIPES:
    require(len(main["q"][n]) == NS, f"pipe {n}: {NS} Q states")

# --- 1. above the threshold a valve is a discount, not a refusal ------------
thr_above = threshold(1.5)
drive_above = drive_of("CreepOpen")
open_root = steady_root(drive_above, 1.5, 0.0)
bill_at_open = dp_valve(open_root, CLOSED_K)
require(abs(thr_above - 6000.0) < 1e-9 and drive_above > thr_above,
        f"the above-threshold group is on drive {drive_above:.1f} Pa against a "
        f"{thr_above:.1f} Pa wall, so the latch can never be its")
require(bill_at_open > 9.0 * drive_above,
        f"the closed valve is worth {bill_at_open:.1f} Pa at the open-valve rate, "
        f"{bill_at_open / drive_above:.1f}x the whole drive")
for n in ABOVE:
    held = sum(main["hold"][n])
    req0 = [i for i, q in enumerate(main["q"][n]) if q == 0.0]
    require(held == 0, f"{n}: held in 0/{NS} states, not one, across a "
                       f"close/hold/open cycle to K = {CLOSED_K:.0e}")
    require(len(req0) <= 1, f"{n}: exactly {len(req0)} zero-rate state(s) -- the "
                            f"initial one, since the wall is not holding it")
    require(min(q for q in main["q"][n][1:]) > 0.0,
            f"{n}: every moving state is strictly nonzero, min "
            f"{min(q for q in main['q'][n][1:]):.6e} m3/s")
require(main["raw"]["CreepShut"] != main["raw"]["CreepK1e6"],
        "a valve that closes and stays closed is not the same trace as one that "
        "was written closed from the start -- the ramp is in the row")
tail_gap = rel(main["q"]["CreepShut"][-1], main["q"]["CreepK1e6"][-1])
require(tail_gap < 1e-11,
        f"but both end on the same root (rel {tail_gap:.3e}): the schedule "
        f"disappears once the ramp is over")

# --- 2. below it: the wall stops it, the schedule only says when ------------
stops = {}
for n in STOP:
    nz = [i for i, q in enumerate(main["q"][n]) if q != 0.0]
    stops[n] = max(nz) + 1 if nz else 0
    require(sum(1 for q in main["q"][n][stops[n]:] if q != 0.0) == 0,
            f"{n}: prints exactly 0 from state {stops[n]} (t = "
            f"{stops[n] * DT:.3f} s) to the end of the run")
require(stops["StopEarly"] < stops["StopMid"] < stops["StopPlain"],
        f"the stop step is the schedule's, monotonically: "
        f"{stops['StopEarly']} / {stops['StopMid']} / {stops['StopPlain']} for a "
        f"ramp ending at 0.15 / 0.6 / never")
require(main["raw"]["StopLate"] == main["raw"]["StopPlain"],
        f"a ramp that starts after the wall has already stopped the pipe is "
        f"invisible in all {NS} printed tokens")
require([bool(h) for h in main["hold"]["StopLate"]]
        == [bool(h) for h in main["hold"]["StopPlain"]],
        "and invisible in the latch flags too, state for state")

# --- 3. no schedule wakes what the wall holds ------------------------------
reopen = int(round(1.65 / DT))
for n in ("StopEarly", "StopMid", "StopLate"):
    after = main["q"][n][reopen:]
    require(all(q == 0.0 for q in after),
            f"{n}: after its valve has finished reopening, all {len(after)} "
            f"states are still exactly 0 -- the K schedule cannot restart it")
    require(sum(main["hold"][n][reopen:]) == len(after),
            f"{n}: and the latch says so in all {len(after)} of them")
nt, np_ = main["q"]["NewtEarly"], main["q"]["NewtPlain"]
dip = min(nt)
require(dip > 0.0,
        f"the Newtonian twin of the same schedule is never at rest: its floor is "
        f"{dip:.6e} m3/s against {np_[0]:.6e} to start")
tail_newt = rel(nt[-1], np_[-1])
require(tail_newt < 1e-9,
        f"and it comes back to its no-valve rate once reopened (tail rel "
        f"{tail_newt:.3e}) -- the asymmetry with the yielded pipe is the "
        f"threshold, not the valve")
require(nt[-1] != np_[-1],
        f"the return is to within {tail_newt:.3e}, not bit-identical: the "
        f"Newtonian row is still relaxing at t = 4 s")

# --- 4. a held pipe does not answer a schedule at all ----------------------
require(main["raw"]["HeldCyc"] == main["raw"]["HeldPlain"],
        f"HeldCyc runs the full close/hold/open cycle on a pipe the wall already "
        f"holds and prints the same token as HeldPlain in all {NS} states")
require(all(q == 0.0 for q in main["q"]["HeldCyc"]),
        "and every one of those tokens is 0")
for n in SCHEDULED:
    require(valve_k(n, 3.0) == K_END[n],
            f"{n}: its schedule is over at t = 3.0 s, where K sits at "
            f"{valve_k(n, 3.0):.1e} -- the value family 7 freezes it at")
for t, want in ((0.075, 2.5e5), (0.1, 5.0e5), (0.125, 7.5e5)):
    got = valve_k("StopEarly", t)
    e = rel(got, want)
    require(e < 1e-12, f"StopEarly: K({t} s) = {got:.6e} is the straight line "
                       f"from 0 to 1e6 over [0.05, 0.15] s (rel {e:.3e} against "
                       f"{want:.1e})")
require(valve_k("StopEarly", 0.05) == 0.0
        and valve_k("StopEarly", 0.15) == CLOSED_K,
        "and the ramp starts on the open K and ends on the closed one, exactly")

# --- 5. the whole 12 x 2001 against the recursion --------------------------
# Two gates. The kernel's Newton quits on a residual scaled by what the row
# itself carries (drive + wall, ~1e4 Pa), so on the states that sit just before
# a stop -- where the rate is ~1e-10 of that scale -- the relative gap is the
# solver's exit precision and not a disagreement about the law. The absolute
# gate is what the two statements actually agree on.
RATE_FLOOR = 1e-8
worst_rel = {}
worst_abs = {}
stop_pred = {}
for n in PIPES:
    pred = recurse(n)
    obs = main["q"][n]
    wr, wa = 0.0, 0.0
    for i in range(NS):
        a = abs(pred[i] - obs[i])
        wa = max(wa, a)
        scale = max(abs(pred[i]), abs(obs[i]))
        if scale >= RATE_FLOOR:
            wr = max(wr, a / scale)
    worst_rel[n] = wr
    worst_abs[n] = wa
    stop_pred[n] = stop_step(pred)
for n in PIPES:
    require(worst_rel[n] < 5e-9, f"{n}: every state above {RATE_FLOOR:.0e} m3/s "
                                f"reproduced by the two laws plus the sampled K, "
                                f"worst rel {worst_rel[n]:.3e}")
gap_worst = max(worst_abs, key=lambda k: worst_abs[k])
require(max(worst_abs.values()) < 5e-15,
        f"and all {len(PIPES) * NS} states agree to "
        f"{max(worst_abs.values()):.3e} m3/s in absolute units (worst is "
        f"{gap_worst}) -- the relative tail is the exit criterion, not the law")
for n in HELD_REST:
    require(worst_abs[n] == 0.0 and worst_rel[n] == 0.0,
            f"{n}: the trace the wall holds at rest is reproduced bit for bit "
            f"(0.0 in both units)")
for n in STOP:
    require(stop_pred[n] == stops[n],
            f"{n}: the recursion stops it on state {stop_pred[n]}, the app prints "
            f"its first 0 on {stops[n]}")

# --- 6. the sampled K is the one at t + dt ---------------------------------
MOVING_RAMP = ("CreepCyc", "CreepShut", "StopEarly", "StopMid", "NewtEarly")
INVISIBLE_RAMP = ("StopLate", "HeldCyc")
for n in MOVING_RAMP:
    lag = 0.0
    ref = recurse(n, k_at_t=True)
    for i, q in enumerate(ref):
        scale = max(abs(q), abs(main["q"][n][i]))
        if scale >= RATE_FLOOR:
            lag = max(lag, abs(q - main["q"][n][i]) / scale)
    require(lag > 100.0 * worst_rel[n],
            f"{n}: sampling K at t instead of t + dt misses the app at rel "
            f"{lag:.3e}, where the t + dt recursion gets {worst_rel[n]:.3e}")
for n in INVISIBLE_RAMP:
    pair = max(rel(a, b) for a, b in zip(recurse(n), recurse(n, k_at_t=True)))
    require(pair == 0.0,
            f"{n}: where the ramp never overlaps a moving state the two sampling "
            f"conventions are one and the same trace ({pair:.3e} apart)")
for n in ("StopEarly", "StopMid"):
    require(stop_step(recurse(n, k_at_t=True)) == stops[n] + 1,
            f"{n}: the one-step reading of that choice -- the app stops it on "
            f"state {stops[n]}, a K sampled at t leaves it moving one step longer "
            f"(state {stop_step(recurse(n, k_at_t=True))})")

# --- 7. the settled row against the steady kernel --------------------------
for n in ("CreepOpen", "CreepCyc", "CreepShut", "CreepK1e6", "NewtPlain",
          "NewtEarly"):
    drive = drive_of(n)
    s = steady[n]["q"]
    e = rel(s, main["q"][n][-1])
    require(e < 2e-9, f"{n}: the marched tail {main['q'][n][-1]:.7e} meets the "
                      f"steady root {s:.7e} at rel {e:.3e}")
    tau0 = PIPES[n][0]
    require(abs(s - steady_root(drive, tau0, K_END[n]))
            <= 2e-9 * max(abs(s), 1e-30),
            f"{n}: and that root is this file's own bisection on drive = "
            f"dp_wall + K rho Q^2/2A^2 with K frozen at {K_END[n]:.1e}")
for n in STOP + HELD_REST:
    rec = steady[n]
    require(rec["q"] == 0.0, f"{n}: the steady kernel reports Q 0 for the same "
                             f"sub-threshold drive {drive_of(n):.0f} Pa")
    require("dpValve" not in rec,
            f"{n}: and prints no dpValve field at all, only dp {rec['dp']:.0f}")
    require(abs(rec["dp"] - drive_of(n)) < 1e-6,
            f"{n}: the drop it does print is the drive, which is what a static "
            f"wall supplies")
require("dpValve" in steady["CreepShut"],
        "the one steady pipe whose valve is closed *and* which is moving does "
        "print dpValve")
v = steady["CreepShut"]
bill = rel(v["dpValve"], dp_valve(v["q"], CLOSED_K))
require(bill < 1e-10, f"CreepShut: dpValve {v['dpValve']:.6f} Pa is the valve's "
                      f"own formula at its own rate (rel {bill:.3e})")
wall_only = v["dp"] - v["dpValve"]
require(wall_only > threshold(1.5),
        f"and what is left for the wall, {wall_only:.3f} Pa, is above its "
        f"{threshold(1.5):.0f} Pa threshold -- which is why this pipe is not held")

# --- the threshold is the same number for every yielded pipe, K or no K ----
require(len(main["thr"]) == len(YIELDED),
        f"{len(main['thr'])} static-hold thresholds printed, one per yielded "
        f"pipe ({len(YIELDED)}), and the two Newtonian pipes print none")
require(all(abs(v - 6000.0) < 1e-9 for v in main["thr"].values()),
        "and every one of them is 4 tau_0 L/d = 6000 Pa, the same number for the "
        "scheduled and the unscheduled pipes alike")

print(f"\n{len(failures)} failure(s)")
if failures:
    for f in failures:
        print("  -", f)
    sys.exit(1)
