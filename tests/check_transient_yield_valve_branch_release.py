#!/usr/bin/env python3
"""transientYieldValveBranchRelease: two yielded walls on ONE free node, a
scheduled drain, and the two things 5.53 could only declare structurally
unreachable -- its boundary (2), a node with more than one pair, and its
boundary (6), a latch that closes again after it opened.

5.53 found the release is one-way and said why: a single wall's own flow is the
only thing that can lift its node back to the pressure it needs, and that flow
goes to zero exactly there, so the node arrives from below and stops. This file
changes the desktop and not the kernel. Put a second wall on the node and the
node acquires a ceiling it did not have: while the *weaker* wall flows, the
network's limit is that wall's own release pressure, and every wall stronger
than it passes under that limit on the way up. The rows per arm are

    rho L/A (Qk - Qk')/dt = (pA - pF) - 4 L tw_k / d - Kk rho|Qk|Qk/(2A^2)
    rho L/A (Qd - Qd')/dt = (pF - pB) - R Qd - Kd(t) rho|Qd|Qd/(2A^2)
    V/BC (pF - pF')        = sum_k Qk - Qd,   Kd(t) = valveKAt(t + dt)

with the latch rule untouched -- the entry latch on `Q_prev == 0` and a drive
under 4 tau_0 L / d (transient.cpp:40-:43) and the in-step clamp for an iterate
carried back through zero (transient.cpp:353-:355) -- and the node row already
degree-generic (transient.cpp:395-:416). Three walls therefore need no code; they
need the case file to stop putting one wall per node.

What the run says, all of it measured here:

  * The latch answers one number in both directions. Eleven walls released, in
    fourteen events, and three of them were put back to sleep; every one of the
    seventeen flips is the step whose entry drive crosses its own 4 tau_0 L / d:
    the releases are caught 1.0 - 16.5 Pa *over* the threshold, the re-latches
    3.1 - 3.6 Pa *under* it. There is no hysteresis to fit, because there is
    nothing between the two branches but the threshold.
  * Which wall can be put back to sleep is decided by the ordering of the
    thresholds, not by the schedule. On the six arms whose node climbs back, a
    wall re-latches iff some *other* wall on the same node has a strictly smaller
    4 tau_0 L / d -- the weakest wall's release pressure is the node's ceiling,
    and a ceiling cannot be crossed from below. Two identical walls therefore
    re-latch neither (arm `i`), and one wall alone re-latches nothing (arm `u`).
    The inequality predicts two walls more than the run delivers, and the run
    gives a different reason for each: arm `t`'s strong wall is never asked (its
    node bottoms 593.53 Pa above the pressure that wall would need), arm `e`'s
    strong wall releases into a node whose steady destination is 620.27 Pa below
    its own release pressure.
  * A held wall's row is silent in the strongest sense the print allows. Arm `t`
    and arm `s` are the same node with and without a permanently-held wall on it,
    and the two print the same numbers for the node and for every flowing pipe.
  * The seal that arrives while every inflow is held takes the node nowhere. Arm
    `e` closes 0.4 s before the strong wall's first flip, and its node then walks
    *down* on the drain's trickle; the re-latch is not an independent event, it
    is downstream of the release it cancels.

Nine families of judgement:
  1. three thresholds, thirteen yielded pipes, no Newtonian pipe, and each flip
     is the step whose entry drive crosses the threshold of its own wall
  2. no hysteresis: the release side and the re-latch side bracket the same
     number, and the gap each leaves is one step's worth of node motion
  3. which walls re-latch, arm by arm, against the rule "some other wall on this
     node is strictly weaker"; arms `i` and `u` are the controls, and the two
     walls the rule predicts but the run refuses are named with their reasons
  4. the node's ceiling is the weakest inflow's release pressure (and pA itself
     when an inflow is Newtonian); every sealed arm rises monotonically and stays
     under it, and every re-latch happens strictly under it
  5. the silent row: arm `t` = arm `s` + one wall that never wakes, and the two
     print the same node and the same three flowing pipes
  6. the long seal: 7001 steps under `1e8` separate "not yet" from "not ever" for
     the weakest wall, and give the strong wall 7295 steps of sleep; the long run
     and the 4 s run share their flip steps
  7. the cancelled seal: arm `e` never lifts its node at all, releases once
     instead of three times, and still settles where arm `a` settles
  8. every state of every arm, in the 4 s and the 18 s run, against the
     checker's own copy of the three rows, the node row and the sampled K
  9. the node balance closes on the printed tokens, and the steady kernel holds
     exactly the walls the transient run never woke
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
DT = 0.002
END = 4.0
NS = int(round(END / DT)) + 1
LONG_END = 18.0
NL = int(round(LONG_END / DT)) + 1
P_A, P_B = 3.0e4, 2.0e4
LENGTH = 20.0
DIAM = 0.02
AREA = 0.25 * math.pi * DIAM ** 2
INERTIA = RHO * LENGTH / AREA
I_DT = INERTIA / DT
KV_SCALE = RHO / (2.0 * AREA * AREA)
WALL_SCALE = 4.0 * LENGTH / DIAM          # dp = WALL_SCALE * tw
R_LAM = 128.0 * MU * LENGTH / (math.pi * DIAM ** 4)
NODE_VOLUME = 0.1
NODE_BULK = 1.0e8
CC = NODE_VOLUME / NODE_BULK              # 1e-9 m3/Pa
P_NODE0 = 2.7e4
RATE_FLOOR = 1.0e-8
NEWTON_TOL = 1.0e-11

CLOSE_STD = (0.8, 1.0, 1.0e8, 2.5, 2.7)
CLOSE_EARLY = (0.2, 0.4, 1.0e8, 2.5, 2.7)
CLOSE_LONG = (0.8, 1.0, 1.0e8, 15.0, 15.2)

# arm: (node, [(inflow, tau_0)], (drain, schedule))
# tau_0 = 0 is not a wall at all -- arm `p`'s Newtonian feed.
ARMS = {
    "a": ("Fa", [("Sa", 2.0), ("Wa", 1.5)], ("Da", CLOSE_STD)),
    "i": ("Fi", [("Sia", 2.0), ("Sib", 2.0)], ("Di", CLOSE_STD)),
    "t": ("Ft", [("St", 2.0), ("Wt", 1.5), ("Tt", 1.0)], ("Dt", CLOSE_STD)),
    "s": ("Fs", [("Ws", 1.5), ("Ts", 1.0)], ("Ds", CLOSE_STD)),
    "u": ("Fu", [("Su", 2.0)], ("Du", CLOSE_STD)),
    "p": ("Fp", [("Sp", 2.0), ("Np", 0.0)], ("Dp", CLOSE_STD)),
    "e": ("Fe", [("Se", 2.0), ("We", 1.5)], ("De", CLOSE_EARLY)),
}
WALLS = {arm: [n for n, t in ARMS[arm][1] if t > 0] for arm in ARMS}
THRESH = {n: WALL_SCALE * t for arm in ARMS for n, t in ARMS[arm][1] if t > 0}
FEEDS = {n: t for arm in ARMS for n, t in ARMS[arm][1]}
DRAINS = {arm: ARMS[arm][2][0] for arm in ARMS}
YIELDED = sorted(THRESH)
NEWTONIAN = sorted(n for arm in ARMS for n, t in ARMS[arm][1] if t == 0)
ALL_PIPE_NAMES = sorted(set(YIELDED) | set(NEWTONIAN) | set(DRAINS.values()))
ALL_NODES = sorted(a[0] for a in ARMS.values())

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


def threshold(tau0):
    return WALL_SCALE * tau0


def valve_k(sched, t):
    """transient.cpp:125-:133 verbatim: linear ramps, closedK held between the
    ramps. No arm here keeps a valve shut at the end, so openStart is never < 0."""
    cs, ce, kc, os_, oe = sched
    if t < cs:
        return 0.0
    if t < ce:
        return 0.0 + (t - cs) / (ce - cs) * (kc - 0.0)
    if t < os_:
        return kc
    if t < oe:
        return kc + (t - os_) / (oe - os_) * (0.0 - kc)
    return 0.0


def sched_of(pipe, long_run=False):
    """The schedule written on this pipe in case.dict / case_long.dict, read by
    nothing here -- the tables above are the same words, written out so a failure
    can be read against the file."""
    for arm, (node, inflows, (drain, sched)) in ARMS.items():
        if pipe == drain:
            return CLOSE_LONG if long_run and arm == "a" else sched
    return None


# --- the rows, stated from this file's copy of the laws ----------------------


def q_of_tw(tw, tau0):
    """The exact laminar discharge a yielded wall stress holds (5.46 at n = 1,
    5.47's profile integral), and with tau0 = 0 it is Hagen-Poiseuille."""
    if tw <= tau0:
        return 0.0
    d = tw - tau0
    r = 0.5 * DIAM
    return (math.pi * r ** 3 * tw ** -3 / MU * d * d
            * (tau0 * tau0 / 2.0 + 2.0 * tau0 * d / 3.0 + d * d / 4.0))


def row_yield(q_prev, drive, kv, tau0):
    """One implicit row of a yielded pipe, static branch included, in wall-stress
    space where the residual is monotone. g(tau0) >= 0 is the in-step clamp of
    transient.cpp:353-:355: the row has no root above the wall's own threshold,
    because the momentum term is already worth more than the drive can pay."""
    thr = threshold(tau0)

    def g(tw):
        q = q_of_tw(tw, tau0)
        return I_DT * (q - q_prev) + WALL_SCALE * tw + kv * KV_SCALE * q * q - drive
    if g(tau0) >= 0.0:
        return 0.0
    lo, hi = tau0, tau0 + 1.0
    while g(hi) < 0.0:
        hi = 2.0 * hi
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo <= 1e-14 * max(mid, 1e-30):
            break
    return q_of_tw(0.5 * (lo + hi), tau0)


def row_laminar(q_prev, drive, kv):
    """The same row for a Newtonian pipe, whose residual is a quadratic in Q with
    a known sign, so the root is written out instead of searched for."""
    b = I_DT + R_LAM
    c = kv * KV_SCALE
    rhs = I_DT * q_prev + drive
    if drive >= 0.0 and rhs >= 0.0:
        if c > 0.0:
            return (-b + math.sqrt(b * b + 4.0 * c * rhs)) / (2.0 * c)
        return rhs / b
    return -row_laminar(-q_prev, -drive, kv)


def inflow_row(name, q_prev, drive, kv):
    tau0 = FEEDS[name]
    if tau0 > 0:
        return row_yield(q_prev, drive, kv, tau0)
    return row_laminar(q_prev, drive, kv)


def march(arm, nsteps, long_run=False, npass=4):
    """The whole trace of one arm: every inflow row, the drain row, the node row
    and the sampled K. The node row is the sum, so an arm with three walls carries
    three terms in it and an arm with a held wall carries that wall as zero."""
    node, inflows, (drain, sched) = ARMS[arm]
    if long_run:
        sched = CLOSE_LONG
    q = {n: 0.0 for n, _ in inflows}
    qd = 0.0
    p = P_NODE0
    logs = {n: [0.0] for n, _ in inflows}
    logs[drain] = [0.0]
    flags = {n: [0] for n in WALLS[arm]}
    plist = [p]
    for i in range(1, nsteps + 1):
        tn = i * DT
        kd = valve_k(sched, tn)
        drive0 = P_A - p
        held = {n: (t > 0 and q[n] == 0.0 and abs(drive0) <= threshold(t))
                for n, t in inflows}
        pn = p
        for _ in range(npass):
            tot = 0.0
            for n, t in inflows:
                qn = 0.0 if held[n] else inflow_row(n, q[n], P_A - pn, 0.0)
                tot += qn
            qdn = row_laminar(qd, pn - P_B, kd)
            pn = p + (DT / CC) * (tot - qdn)
        newq = {}
        for n, t in inflows:
            newq[n] = 0.0 if held[n] else inflow_row(n, q[n], P_A - pn, 0.0)
        qdn = row_laminar(qd, pn - P_B, kd)
        for n in q:
            q[n] = newq[n]
        qd = qdn
        p = pn
        for n, _ in inflows:
            logs[n].append(newq[n])
        for n in WALLS[arm]:
            flags[n].append(1 if held[n] else 0)
        logs[drain].append(qdn)
        plist.append(p)
    return logs, flags, plist


def parse_transient(stdout, pipes, nodes):
    out = {"q": {n: [] for n in pipes}, "raw": {n: [] for n in pipes},
           "hold": {n: [] for n in pipes}, "p": {n: [] for n in nodes},
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


def parse_steady(stdout):
    flows, nodes = {}, {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 4 and tok[0] == "pipe" and tok[2] == "Q":
            rec = {"q": float(tok[3]), "raw": tok[3]}
            for i, t in enumerate(tok):
                if t == "dp":
                    rec["dp"] = float(tok[i + 1])
            flows[tok[1]] = rec
        elif len(tok) >= 4 and tok[0] == "node" and tok[2] == "p":
            nodes[tok[1]] = float(tok[3])
    return flows, nodes


def events(flags, ps):
    """Every latch flip after the first step, with the drive the step started on
    and the drive the step before it started on."""
    out = []
    prev = None
    for i, v in enumerate(flags):
        if v != prev and i >= 2:
            out.append({"step": i, "hold": v, "entry": P_A - ps[i - 1],
                        "prev": P_A - ps[i - 2]})
        prev = v
    return out


# --- runs --------------------------------------------------------------------
proc = run(TRANSIENT, CASE)
if proc.returncode != 0:
    raise SystemExit(f"transient solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
main = parse_transient(proc.stdout, ALL_PIPE_NAMES, ALL_NODES)
proc = run(TRANSIENT, f"{CASE_DIR}/case_long.dict")
if proc.returncode != 0:
    raise SystemExit(f"long solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
longn = parse_transient(proc.stdout, ["Sa", "Wa", "Da"], ["Fa"])
proc = run(STEADY, f"{CASE_DIR}/case_steady.dict")
if proc.returncode != 0:
    raise SystemExit(f"steady solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
steady, steady_p = parse_steady(proc.stdout)

EV = {}
for arm in ARMS:
    node = ARMS[arm][0]
    EV[arm] = {n: events(main["hold"][n], main["p"][node]) for n in WALLS[arm]}
RELATCH = {n: [e for e in EV[arm][n] if e["hold"] == 1]
           for arm in ARMS for n in WALLS[arm]}
RELEASING = {}
for arm in ARMS:
    for n in WALLS[arm]:
        ev = [e for e in EV[arm][n] if e["hold"] == 0]
        RELEASING[n] = ev[0] if ev else None
EVERY_FLIP = [e for arm in ARMS for n in WALLS[arm] for e in EV[arm][n]]
LONG_EV = {n: events(longn["hold"][n], longn["p"]["Fa"]) for n in ("Sa", "Wa")}

# --- 1. three thresholds, thirteen walls, and the flip that is theirs --------
thr_by_value = {}
for n, v in main["thr"].items():
    thr_by_value.setdefault(v, []).append(n)
require(len(main["thr"]) == len(WALLS["a"]) + len(WALLS["i"]) + len(WALLS["t"])
        + len(WALLS["s"]) + len(WALLS["u"]) + len(WALLS["p"]) + len(WALLS["e"])
        and set(main["thr"]) == {n for arm in ARMS for n in WALLS[arm]},
        f"a static-hold threshold is printed for exactly the "
        f"{sum(len(WALLS[a]) for a in ARMS)} yielded pipes of the seven arms and "
        f"for none of the {len(DRAINS) + 1} Newtonian pipes "
        f"({sorted(main['thr'])})")
require(sorted(thr_by_value) == [threshold(1.0), threshold(1.5), threshold(2.0)],
        f"the {len(main['thr'])} thresholds take exactly three values, "
        f"4 tau_0 L / d = {threshold(2.0):.0f} / {threshold(1.5):.0f} / "
        f"{threshold(1.0):.0f} Pa for tau_0 = 2.0 / 1.5 / 1.0, in "
        + ", ".join(f"{v:.0f}x{len(thr_by_value[v])}" for v in sorted(thr_by_value)))
require(all(thr_by_value[threshold(2.0)][i] in
            ("Sa", "Sia", "Sib", "St", "Su", "Sp", "Se")
            for i in range(len(thr_by_value[threshold(2.0)]))),
        "and each wall's threshold is its own tau_0, not its neighbour's: the "
        f"seven 8000 Pa lines are {sorted(thr_by_value[threshold(2.0)])}, the four "
        f"6000 Pa lines {sorted(thr_by_value[threshold(1.5)])} and the two 4000 Pa "
        f"lines {sorted(thr_by_value[threshold(1.0)])}")
require(all(e["hold"] == 0 and e["prev"] <= THRESH[n] < e["entry"]
            for n in RELEASING if RELEASING[n] for e in [RELEASING[n]])
        and all(e["entry"] <= THRESH[n] < e["prev"]
                for arm in ARMS for n in WALLS[arm] for e in RELATCH[n]),
        f"every flip is the step whose entry drive crosses its own wall's "
        f"threshold, from the correct side: {sum(1 for n in RELEASING if RELEASING[n])} "
        f"walls released at least once -- each one's *first* release has the "
        f"previous sample at or under 4 tau_0 L / d and this one over it -- and "
        f"{sum(len(v) for v in RELATCH.values())} of them were put back to sleep, "
        f"the mirror image, for {len(EVERY_FLIP)} flips in the file")
never = sorted(n for n in RELEASING if RELEASING[n] is None)
require(never == ["Sp", "St"],
        f"two walls are never released inside 4 s and the run says so: {never}, "
        f"which print the token 0 in all {NS} states")

# --- 2. one number, both directions ------------------------------------------
ALL_REL = [(n, e) for arm in ARMS for n in WALLS[arm] for e in EV[arm][n]
           if e["hold"] == 0]
ALL_RL = [(n, e) for arm in ARMS for n in WALLS[arm] for e in RELATCH[n]]
gap_over = [e["entry"] - THRESH[n] for n, e in ALL_REL]
gap_under = [THRESH[n] - e["entry"] for n, e in ALL_RL]
require(ALL_REL and max(gap_over) < 20.0 and min(gap_over) > 0.0,
        f"the {len(gap_over)} releases start above the wall's own threshold by at "
        f"most {max(gap_over):.2f} Pa -- one step's worth of node -- and by "
        f"{min(gap_over):.2f} Pa at the tightest")
require(ALL_RL and max(gap_under) < 5.0 and min(gap_under) > 0.0,
        f"the {len(gap_under)} re-latches start *under* the same number, by "
        f"{max(gap_under):.2f} Pa at most and {min(gap_under):.2f} Pa at least")
require(all(abs(e["prev"] - THRESH[n]) < 20.0
            for n, e in ALL_REL + ALL_RL),
        "and the sample one step before is on the other side of the same number "
        f"by the same one step, i.e. within "
        f"{max(abs(e['prev'] - THRESH[n]) for n, e in ALL_REL + ALL_RL):.2f} Pa")
require(all(abs((P_A - e["entry"]) - (P_A - THRESH[n])) < 20.0
            for n, e in ALL_REL + ALL_RL),
        "so the latch is read off one pressure in both directions and there is no "
        "hysteresis to fit: every flip happens within "
        f"{max(max(gap_over), max(gap_under)):.2f} Pa of pA - 4 tau_0 L / d, the "
        f"three re-latches at "
        + ", ".join(f"{n} {P_A - e['entry']:.2f} Pa" for n, e in ALL_RL)
        + f" against release pressures {P_A - THRESH['Sa']:.0f} / "
        f"{P_A - THRESH['Wa']:.0f} / {P_A - THRESH['Tt']:.0f} Pa")

# --- 3. who may go back to sleep ---------------------------------------------
law = {n: any(0 < THRESH[m] < THRESH[n] for m in WALLS[arm])
       for arm in ARMS for n in WALLS[arm]}
require(sorted(RELEASING) == sorted(law),
        f"the rule is stated for all {len(law)} walls of the file")
did = {n: bool(RELEASING[n]) and bool(RELATCH[n]) for n in RELEASING}
# The inequality is a statement about the ceiling a node climbs back to, so it
# can only be tested on an arm whose node does climb back: the six arms whose
# seal ends while the file is still running. Arm `e` is not one of them and the
# next require says why in its own numbers.
predicted = sorted(n for arm in ("a", "i", "t", "s", "u", "p")
                   for n in WALLS[arm] if RELEASING[n] and law[n])
actual = sorted(n for n in RELEASING if did[n])
require(predicted == actual,
        "a wall re-latches iff it released and some *other* wall on its node has "
        f"a strictly smaller threshold -- on the six arms that climb back that is "
        f"exactly the set that does, {actual}: "
        + ", ".join(f"{n}={'yes' if did[n] else 'no'}" for n in sorted(did)
                    if RELEASING[n]))
odd = sorted(n for arm in ARMS for n in WALLS[arm] if law[n] and not did.get(n, False))
require(odd == ["Se", "St"],
        "and the two walls the inequality alone would put back to sleep are the "
        "two this file refuses, each for a different reason: St never released "
        f"(its node bottoms at {min(main['p']['Ft']):.2f} Pa, "
        f"{min(main['p']['Ft']) - (P_A - THRESH['St']):.2f} Pa above the pressure "
        "Sa and Wa needed), and Se released once at step "
        f"{RELEASING['Se']['step']} into a node that never gets back under it -- "
        f"the largest pressure Fe reaches after that release is "
        f"{max(main['p']['Fe'][RELEASING['Se']['step'] + 1:]):.2f} Pa, "
        f"{P_A - THRESH['Se'] - max(main['p']['Fe'][RELEASING['Se']['step'] + 1:]):.2f} Pa "
        f"short of its own {P_A - THRESH['Se']:.0f} Pa, and its steady destination "
        f"{steady_p['Fe']:.2f} Pa is short by "
        f"{P_A - THRESH['Se'] - steady_p['Fe']:.2f} Pa")
require(not any(RELATCH[n] for n in ("Wa", "Tt", "Ts")),
        "the weakest wall of each arm is the one that never re-latches (Wa 6000 "
        "among 8000/6000, Tt and Ts 4000 among 8000/6000/4000), and it is not a "
        "matter of degree: its release pressure is the node's ceiling")
require(not any(RELATCH[n] for n in ("Sia", "Sib")),
        "the ordering control: arm `i` has two walls of the same threshold, so "
        "neither is strictly weaker and neither re-latches -- they release on the "
        f"same step {RELEASING['Sia']['step']} and sleep through the same seal "
        "arm `a` wakes Sa in")
require(all(len([e for e in EV[arm][n] if e["hold"] == 0]) == 1
            and not RELATCH[n] for arm, n in (("u", "Su"), ("e", "We"))),
        "and the arms whose weakest flowing wall has nothing weaker beside it "
        f"release once and never re-latch: Su at {RELEASING['Su']['step']}, We at "
        f"{RELEASING['We']['step']} -- We is the interesting one, because it is "
        "the wall that keeps arm `e`'s node off Se's release pressure")
sa_steps = [e["step"] for e in EV["a"]["Sa"]]
require([e["hold"] for e in EV["a"]["Sa"]] == [0, 1, 0],
        f"the reference arm `a` is one wall waking, sleeping and waking again: "
        f"steps {sa_steps} at t = "
        + ", ".join(f"{main['times'][s]:.3f}" for s in sa_steps) + " s")
require(sum(main["hold"]["Sa"]) == (sa_steps[0] - 1) + (sa_steps[2] - sa_steps[1]),
        f"and its sleep is a counted interval, not a story: {sum(main['hold']['Sa'])} "
        f"of {NS} states held = {sa_steps[0] - 1} before the first release plus "
        f"{sa_steps[2] - sa_steps[1]} between the re-latch and the second release")

# --- 4. the ceiling ----------------------------------------------------------
ceil = {arm: (P_A - min(threshold(t) for n, t in ARMS[arm][1] if t > 0)
              if all(t > 0 for _, t in ARMS[arm][1]) else P_A)
        for arm in ARMS}
seal = {arm: (int(CLOSE_STD[1] / DT), int(CLOSE_STD[3] / DT) + 1)
        for arm in ARMS}
seal["e"] = (int(CLOSE_EARLY[1] / DT), int(CLOSE_EARLY[3] / DT) + 1)
for arm in ("a", "i", "t", "s", "u", "p"):
    lo, hi = seal[arm]
    ps = main["p"][ARMS[arm][0]][lo:hi]
    rising = all(b > a for a, b in zip(ps, ps[1:]))
    require(rising and max(ps) < ceil[arm],
            f"arm {arm}: with the drain shut the node rises monotonically from "
            f"{ps[0]:.2f} to {ps[-1]:.2f} Pa and never reaches the ceiling "
            f"{ceil[arm]:.0f} Pa implied by its weakest inflow "
            f"(deficit {ceil[arm] - max(ps):.2f} Pa)")
require(ceil["p"] == P_A,
        "arm `p`'s ceiling is pA itself, not a wall pressure, because its second "
        f"inflow is Newtonian: its node climbs to {max(main['p']['Fp'][seal['p'][0]:seal['p'][1]]):.2f} "
        f"Pa, and the wall on it is never released -- the closest its drive gets "
        f"to 4 tau_0 L / d = 8000 Pa is {P_A - min(main['p']['Fp']):.2f} Pa")
rl_under = [(n, P_A - e["entry"], ceil[arm]) for arm in ARMS for n in WALLS[arm]
            for e in RELATCH[n]]
require(all(p < c for _, p, c in rl_under),
        "and every re-latch happens strictly under its arm's ceiling, which is "
        "the point: the node has room to climb past the wall's own release "
        "pressure only because something weaker is still flowing -- "
        + ", ".join(f"{n} {p:.2f} < {c:.0f}" for n, p, c in rl_under))
require(min(main["p"]["Ft"]) > P_A - THRESH["St"],
        f"arm `t` never asks its strongest wall at all: its node bottoms at "
        f"{min(main['p']['Ft']):.2f} Pa, "
        f"{min(main['p']['Ft']) - (P_A - THRESH['St']):.2f} Pa above the "
        f"{P_A - THRESH['St']:.0f} Pa Sa/Wa would need, because two weaker walls "
        "keep feeding it")

# --- 5. the silent row -------------------------------------------------------
require(set(main["raw"]["St"]) == {"0"} and set(main["raw"]["Sp"]) == {"0"},
        f"the two walls that never wake print the single token `0` in all {NS} "
        "states -- the held row is the forced zero, not a small number")
same_p = max(abs(a - b) for a, b in zip(main["p"]["Ft"], main["p"]["Fs"]))
require(same_p == 0.0,
        f"arm `t` and arm `s` are the same node with and without that silent row, "
        f"and the two node traces are the same numbers in every state "
        f"(largest difference {same_p:.1e} Pa)")
require(all(main["raw"][a] == main["raw"][b]
            for a, b in (("Wt", "Ws"), ("Tt", "Ts"), ("Dt", "Ds"))),
        "and so is every flowing pipe token for token: Wt/Ws, Tt/Ts and Dt/Ds "
        f"print identical strings in all {NS} states, while St's own row still "
        "occupies a column in the node matrix")
require(abs(steady_p["Ft"] - steady_p["Fs"]) < 1e-9
        and steady["St"]["q"] == 0.0 and steady["Sp"]["q"] == 0.0,
        f"the steady kernel says the same thing without a latch: Ft and Fs settle "
        f"to {steady_p['Ft']:.4f} / {steady_p['Fs']:.4f} Pa and the silent wall "
        "reports exactly zero rate")

# --- 6. the long seal -------------------------------------------------------
long_ev_sa = LONG_EV["Sa"]
require([e["step"] for e in long_ev_sa] == [333, 507, 7802],
        f"the 18 s run gives arm `a`'s strong wall three flips, steps "
        f"{[e['step'] for e in long_ev_sa]} at t = "
        + ", ".join(f"{longn['times'][e['step']]:.3f}" for e in long_ev_sa)
        + " s, and its sleep is 7295 of them")
require(sum(longn["hold"]["Sa"]) == 7627 and sum(longn["hold"]["Wa"]) == 179,
        f"held-state counts: Sa {sum(longn['hold']['Sa'])} of {NL}, Wa only the "
        f"first {sum(longn['hold']['Wa'])} -- the weakest wall of this arm is "
        "asleep at the start and awake to the end")
lp = longn["p"]["Fa"]
release_weak = P_A - threshold(1.5)
ceiling = [release_weak - p for p in lp]
deficit = [(i, ceiling[i]) for i in (1000, 2000, 3000, 4000, 5000, 6000, 7000)]
long_seal = (int(CLOSE_STD[1] / DT), int(CLOSE_LONG[3] / DT) + 1)
require(all(d > 0 for _, d in deficit)
        and all(b < a for (_, a), (_, b) in zip(deficit, deficit[1:]))
        and long_seal[0] < deficit[0][0] and deficit[-1][0] < long_seal[1],
        f"the node climbs all {long_seal[1] - long_seal[0]} steps of the long "
        f"seal toward {release_weak:.0f} Pa and never arrives -- the deficit "
        f"samples sit inside that window, from {deficit[0][1]:.2f} Pa at step "
        f"{deficit[0][0]} to {deficit[-1][1]:.2f} Pa at step {deficit[-1][0]}")
require(all(b > a for a, b in zip(lp[long_seal[0]:long_seal[1]],
                                  lp[long_seal[0] + 1:long_seal[1]])),
        "monotonically over the whole shut window, with no flip for the 6000 Pa "
        f"wall to make: the last sample before the reopen is "
        f"{lp[long_seal[1] - 1]:.4f} Pa, {ceiling[long_seal[1] - 1]:.2f} Pa under "
        "its own release pressure")
shrink = [deficit[k + 1][1] / deficit[k][1] for k in range(len(deficit) - 1)]
require(all(s < 1.0 for s in shrink)
        and all(b > a for a, b in zip(shrink, shrink[1:]))
        and shrink[-1] > 0.99,
        "and the remaining gap stops closing like a schedule and starts closing "
        f"like an asymptote: the per-1000-step shrink of the deficit runs "
        + ", ".join(f"{s:.5f}" for s in shrink) + f" -- the first thousand steps "
        f"take it {deficit[0][1]:.1f} -> {deficit[1][1]:.1f} Pa, the last only "
        f"{deficit[-2][1]:.1f} -> {deficit[-1][1]:.1f} Pa, so `not yet` has "
        "become `not ever`")
same_flips = ([e["step"] for e in LONG_EV["Wa"]] == [180]
              and RELEASING["Wa"]["step"] == 180
              and LONG_EV["Sa"][0]["step"] == RELEASING["Sa"]["step"]
              and LONG_EV["Sa"][1]["step"] == RELATCH["Sa"][0]["step"])
require(same_flips,
        "the long run and the 4 s run agree on the flips they share: Wa releases "
        f"at 180 in both, Sa releases at {LONG_EV['Sa'][0]['step']} and re-latches "
        f"at {LONG_EV['Sa'][1]['step']} in both")
dl = max(abs(a - b) for a, b in zip(lp[:seal["a"][1]], main["p"]["Fa"][:seal["a"][1]]))
require(dl < 1e-6,
        f"and the two runs agree on the whole shared window as numbers, though not "
        f"as tokens: the largest node difference over the first {seal['a'][1]} "
        f"states is {dl:.2e} Pa (relative {dl / P_A:.1e}), because 21 pipes and 3 "
        "pipes do not sum the same dense rows in the same order")

# --- 7. the seal that arrives too early -------------------------------------
elo, ehi = seal["e"]
eps = main["p"]["Fe"][elo:ehi]
require(all(b < a for a, b in zip(eps, eps[1:])),
        f"arm `e` closes 0.4 s after the start, i.e. while *both* its walls are "
        f"held, and its node then falls for the whole {ehi - elo} steps of the "
        f"seal: {eps[0]:.4f} -> {eps[-1]:.4f} Pa")
drops = [a - b for a, b in zip(eps, eps[1:])]
qd_seal = main["q"][ARMS["e"][2][0]][elo:ehi]
per_step = DT * sum(qd_seal) / len(qd_seal) / CC
require(all(d > 0 for d in drops) and max(drops) < 0.23,
        f"at a rate set by the drain's own trickle through `1e8` and not by any "
        f"wall: {min(drops):.4f} - {max(drops):.4f} Pa per step, i.e. "
        f"{min(drops) / DT:.1f} - {max(drops) / DT:.1f} Pa/s, and dt * Q_drain / "
        f"Cc for this file's drained mean {sum(qd_seal) / len(qd_seal):.4e} m3/s "
        f"is {per_step:.4f} Pa -- the node row with both inflows at zero")
require(RELEASING["Se"]["step"] == 1624 and not RELATCH["Se"],
        f"so the strong wall of arm `e` releases once, at step "
        f"{RELEASING['Se']['step']} (t = {main['times'][RELEASING['Se']['step']]:.3f} s), "
        "and never re-latches: the second flip of arm `a` is not an independent "
        "event, it is downstream of the release this seal cancelled")
require(abs(steady_p["Fe"] - steady_p["Fa"]) < 1e-9
        and abs(main["p"]["Fe"][-1] - main["p"]["Fa"][-1]) < 6.0,
        f"and both arms end at the same place once the schedules freeze: steady "
        f"gives Fa/Fe {steady_p['Fa']:.4f} / {steady_p['Fe']:.4f} Pa and the two "
        f"marched tails are at {main['p']['Fa'][-1]:.4f} / {main['p']['Fe'][-1]:.4f} Pa")

# --- 8. every state, from this file's own rows ------------------------------
worst_q = 0.0
worst_p = 0.0
worst_node_row = 0.0
flag_mismatch = []
flip_mismatch = []
for arm in ARMS:
    node = ARMS[arm][0]
    logs, flags, plist = march(arm, NS - 1)
    for i in range(len(plist)):
        worst_p = max(worst_p, rel(plist[i], main["p"][node][i]))
    for n, _ in ARMS[arm][1] + [ARMS[arm][2]]:
        for a, b in zip(logs[n], main["q"][n]):
            if max(abs(a), abs(b)) > RATE_FLOOR:
                worst_q = max(worst_q, rel(a, b))
    for n in WALLS[arm]:
        for i, (a, b) in enumerate(zip(flags[n], main["hold"][n])):
            if a != b and i >= 1:
                flag_mismatch.append((n, i, a, b))
        mine = [i for i in range(2, len(flags[n]))
                if flags[n][i] != flags[n][i - 1]]
        theirs = [e["step"] for e in EV[arm][n]]
        if mine != theirs:
            flip_mismatch.append((n, mine[:4], theirs[:4]))
logs, flags, plist = march("a", NL - 1, long_run=True)
for i in range(len(plist)):
    worst_p = max(worst_p, rel(plist[i], longn["p"]["Fa"][i]))
for n in ("Sa", "Wa", "Da"):
    for a, b in zip(logs[n], longn["q"][n]):
        if max(abs(a), abs(b)) > RATE_FLOOR:
            worst_q = max(worst_q, rel(a, b))
for n in ("Sa", "Wa"):
    mine = [i for i in range(2, len(flags[n])) if flags[n][i] != flags[n][i - 1]]
    theirs = [e["step"] for e in LONG_EV[n]]
    if mine != theirs:
        flip_mismatch.append((n + " long", mine[:4], theirs[:4]))
require(not flag_mismatch and not flip_mismatch,
        f"the checker's own three rows, node row and sampled K reproduce every "
        f"latch state of {(NS - 1) * len(ARMS) + (NL - 1)} steps: no printed flag "
        f"differs and the flip lists agree step for step "
        f"({len(flip_mismatch)} mismatches)")
require(worst_q < 5e-9,
        f"with the largest rate difference in the file {worst_q:.3e} relative "
        f"({(NS - 1) * len(ARMS) + (NL - 1)} states x 21 pipes, |Q| > "
        f"{RATE_FLOOR:.0e} m3/s)")
require(worst_p < 2e-11,
        f"and the largest node difference {worst_p:.3e} relative -- including the "
        "two re-latch steps, where the row that is solved changes shape")

# --- 9. continuity and the steady destination --------------------------------
abs_res = 0.0
rel_res = 0.0
for arm in ARMS:
    node, inflows, drain = ARMS[arm][0], [n for n, _ in ARMS[arm][1]], ARMS[arm][2][0]
    ps = main["p"][node]
    for i in range(1, len(ps)):
        tin = sum(main["q"][n][i] for n in inflows)
        out = main["q"][drain][i]
        r = CC / DT * (ps[i] - ps[i - 1]) - (tin - out)
        abs_res = max(abs_res, abs(r))
        if max(abs(tin), abs(out)) > RATE_FLOOR:
            rel_res = max(rel_res, abs(r) / max(abs(tin), abs(out)))
require(abs_res < 1e-12,
        f"the node rows close on the printed tokens alone, with the held walls "
        f"entering as the zeros they print: largest continuity residual "
        f"{abs_res:.3e} m3/s against through-flows of 1e-6 m3/s "
        f"(relative {rel_res:.3e}, printed precision)")
def held_pairing(parsed):
    """(how many states the print flags as held, the tokens those states show,
    and any held state whose rate is not exactly zero)."""
    idx = [(n, i) for n in parsed["hold"]
           for i, f in enumerate(parsed["hold"][n]) if f]
    return len(idx), {parsed["raw"][n][i] for n, i in idx}, \
        [(n, i, parsed["raw"][n][i]) for n, i in idx if parsed["q"][n][i] != 0.0]


n_main, tok_main, odd_main = held_pairing(main)
n_long, tok_long, odd_long = held_pairing(longn)
require(not odd_main and not odd_long and tok_main == {"0"} and tok_long == {"0"},
        "which is a statement about volume, not about rounding: every one of the "
        f"{n_main} held states of the 4 s run and {n_long} of the 18 s run prints "
        "the single token 0, so a wall the print calls held leaves the node row "
        "with the zero it actually marched, not with the row of the law it refused")
require(sum(main["hold"]["St"]) == sum(main["hold"]["Sp"]) == NS - 1
        and all(v == 0.0 for v in main["q"]["St"] + main["q"]["Sp"]),
        f"the two never-woke walls are the whole 4 s file minus the t 0 sample "
        f"({NS - 1} of {NS} states held each) and they never once showed a rate, "
        f"while the re-latch step {RELATCH['Sa'][0]['step']} of arm `a` prints the "
        f"pair side by side: {main['raw']['Sa'][RELATCH['Sa'][0]['step']]} for the "
        "wall going back to sleep next to "
        f"{main['raw']['Wa'][RELATCH['Sa'][0]['step']]} for the neighbour that is "
        "still carrying it")
steady_zero = sorted(n for n in THRESH if steady[n]["q"] == 0.0)
never_woken = sorted(n for n in RELEASING if RELEASING[n] is None)
require(steady_zero == never_woken,
        f"the steady kernel refuses to pay exactly the walls the 4 s window never "
        f"woke: it reports 0 for {steady_zero} at settled drives "
        + ", ".join(f"{n} {steady[n]['dp']:.2f}" for n in steady_zero)
        + f" against thresholds {', '.join(f'{n} {THRESH[n]:.0f}' for n in steady_zero)}")
require(all(steady[n]["q"] > 0.0 for n in THRESH if n not in steady_zero),
        f"and a nonzero rate for the other {len(THRESH) - len(steady_zero)} walls, "
        f"the smallest being Sa's {steady['Sa']['q']:.4e} m3/s at "
        f"{steady['Sa']['dp']:.2f} Pa")
tails = [(arm, main["p"][ARMS[arm][0]][-1], steady_p[ARMS[arm][0]])
         for arm in ("a", "i", "t", "s", "u")]
require(all(abs(m - s) < 900.0 for _, m, s in tails),
        "and the marched tails are heading there: "
        + ", ".join(f"{arm} {m:.2f} vs {s:.2f}" for arm, m, s in tails)
        + " Pa (arm p and arm e are named by their own families above; p is "
        f"mid-flight toward its settled {steady_p['Fp']:.2f} Pa node and e is "
        f"{abs(main['p']['Fe'][-1] - main['p']['Fa'][-1]):.2f} Pa from arm a's)")

print()
if failures:
    raise SystemExit(f"{len(failures)} check(s) failed")
print(f"transientYieldValveBranchRelease: {len(ARMS)} arms, "
      f"{len(THRESH)} walls, {(NS - 1) * len(ARMS) + (NL - 1)} marched steps, "
      f"{len(EVERY_FLIP)} latch flips of the 4 s file and "
      f"{sum(len(v) for v in LONG_EV.values())} of the 18 s file verified")
