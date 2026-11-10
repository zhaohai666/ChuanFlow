#!/usr/bin/env python3
"""5.55: the branch on the outflow side -- a yield wall whose drive is `pF - pB`.

Usage: check_transient_yield_valve_split_release.py TRANSIENT STEADY CASE

Everything 5.46 through 5.54 put a yielded pipe on the *inflow* side of a free
node: its drive had one end pinned at the high-pressure source, it could only
RAISE that node, and the pressure it stops at -- `pA - 4 tau_0 L / d` -- is a
ceiling it approaches from below. That single fact is what made 5.53 call the
release one-way and what made 5.54 need a SECOND wall to break through it. This
file turns the arm around. A drain wall can only LOWER its node, and it stops at
`pB + 4 tau_0 L / d`, a FLOOR approached from above, so every claim about the
latch has a mirror image to be checked against rather than assumed:

  * the re-latch of 5.54's ceiling law, read from the floor side (`Wz`),
  * the same event produced by a Newtonian neighbour instead of a weaker wall
    (`Wb`), which 5.54's arm `p` only ever showed in the release direction,
  * the arms where nothing crosses the limit at all and the wall therefore never
  gets to report a held state (`Ef`, `Tz`, `Em`),
  * the dead band: a chain of two walls whose flow windows do not overlap, which
    is a fixed point of the whole solve and an *exact* steady solution,
  * a node with two free pressures and the schedule sitting on the row between
    them, so a row whose drive has two unknowns is in two node rows at once.

The kernel is not touched to get any of this. `holdsStatic` compares the
ABSOLUTE drive against the geometry-only threshold (transient.cpp:40-:43),
`drivePa` is signed by the row (transient.cpp:52-:56), `frictionDp` copies the
drive's sign at Q = 0 (transient.cpp:72-:81), and the node row sums whatever the
topology hands it (transient.cpp:395-:416).

The families below are stated as measured numbers. The one claim that is not a
measurement -- which wall a network is even allowed to put back to sleep -- is
computed from the topology and the thresholds alone, before the trace is read,
and then compared with the events the run reports; family 3 prints both sets and
the row of the carrier that each wall is matched with.
"""

import math
import re
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
SOURCE_SCALE = 6.0e-5                     # the `sources` block in both steady files

CLOSE_STD = (0.8, 1.0, 1.0e8, 2.5, 2.7)
CLOSE_LONG = (0.8, 1.0, 1.0e8, 15.0, 15.2)

# arm: ([free nodes], [(pipe, from, to, tau_0, schedule)])
# tau_0 = 0 is not a wall at all. `to` or `from` equal to a free node name is
# where the row pays the node row; A and B are the fixed pressures.
ARMS = {
    "a": (["Fa"], [("Sa", "A", "Fa", 2.0, None),
                   ("Wa", "A", "Fa", 1.5, None),
                   ("Da", "Fa", "B", 0.0, CLOSE_STD)]),
    "b": (["Fb"], [("Sb", "A", "Fb", 2.0, None),
                   ("Wb", "Fb", "B", 1.5, None),
                   ("Db", "Fb", "B", 0.0, CLOSE_STD)]),
    "c": (["Fc"], [("Sc", "A", "Fc", 2.0, None),
                   ("Dc", "Fc", "B", 0.0, CLOSE_STD),
                   ("Nc", "Fc", "B", 0.0, None)]),
    "d": (["Fd"], [("Sd", "A", "Fd", 2.0, None),
                   ("Ed", "Fd", "B", 2.0, None)]),
    "f": (["Ff"], [("Sf", "A", "Ff", 2.0, None),
                   ("Ef", "Ff", "B", 1.5, None)]),
    "z": (["Fz"], [("Sz", "A", "Fz", 2.0, None),
                   ("Wz", "Fz", "B", 1.5, None),
                   ("Tz", "Fz", "B", 1.0, None)]),
    "m": (["Fm1", "Fm2"], [("Sm", "A", "Fm1", 2.0, None),
                           ("Mm", "Fm1", "Fm2", 0.0, CLOSE_STD),
                           ("Em", "Fm2", "B", 1.5, None)]),
}
PIPES = [n for arm in ARMS for (n, _, _, _, _) in ARMS[arm][1]]
YIELDED = [(n, t, f, to) for arm in ARMS for (n, f, to, t, _) in ARMS[arm][1] if t > 0]
WALLS = sorted(n for n, _, _, _ in YIELDED)
THRESH = {n: WALL_SCALE * t for n, t, _, _ in YIELDED}
TAU = {n: t for (n, t, f, to) in YIELDED}
# Which side of its free node a wall sits on: `in` if it brings flow to the
# node from the fixed source, `out` if it takes flow away to the fixed sink.
SIDE = {n: ("in" if f == "A" else "out") for n, t, f, to in YIELDED}
NODE_OF = {n: (f if f.startswith("F") else to) for n, t, f, to in YIELDED}
LIMIT = {n: (P_A - THRESH[n]) if SIDE[n] == "in" else (P_B + THRESH[n])
         for n in WALLS}
DRIVER = {n: (P_A - LIMIT[n]) if SIDE[n] == "in" else (LIMIT[n] - P_B)
          for n in WALLS}     # the drive each wall needs to reach its own limit
SCHEDULED = {n: s for arm in ARMS for (n, _, _, _, s) in ARMS[arm][1] if s}
NEWTONIAN = sorted(n for n in PIPES if THRESH.get(n, 0) == 0)
ALL_NODES = sorted({x for arm in ARMS for x in ARMS[arm][0]}) + ["A", "B"]
LONG_ARMS = ("b", "f", "z")

failures = []
CHECKS = 0


def require(cond, what):
    global CHECKS
    CHECKS += 1
    print("ok:  " if cond else "FAIL:", what)
    if not cond:
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
    """transient.cpp:125-:133 verbatim: linear ramps, closedK held between them."""
    cs, ce, kc, os_, oe = sched
    if t < cs:
        return 0.0
    if t < ce:
        return (t - cs) / (ce - cs) * kc
    if t < os_:
        return kc
    if t < oe:
        return kc + (t - os_) / (oe - os_) * (0.0 - kc)
    return 0.0


def q_of_tw(tw, tau0):
    """The exact laminar discharge a yielded wall stress holds (5.46 at n = 1):
    with tau0 = 0 this is Hagen-Poiseuille."""
    if tw <= tau0:
        return 0.0
    d = tw - tau0
    r = 0.5 * DIAM
    return (math.pi * r ** 3 * tw ** -3 / MU * d * d
            * (tau0 * tau0 / 2.0 + 2.0 * tau0 * d / 3.0 + d * d / 4.0))


def row_yield(q_prev, drive, kv, tau0):
    """One implicit row of a yielded pipe in wall-stress space, where the
    residual is monotone. g(tau0) >= 0 is the in-step clamp of
    transient.cpp:353-:355: the row has no root once the wall's own threshold
    outruns what the drive can pay. Signed, because a drain wall's drive is read
    the other way across the node and a wall can be asked to run backwards."""
    if drive < 0.0:
        return -row_yield(-q_prev, -drive, kv, tau0)
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
    """The same row without a wall: a quadratic in Q with a known sign."""
    b = I_DT + R_LAM
    c = kv * KV_SCALE
    rhs = I_DT * q_prev + drive
    if drive >= 0.0 and rhs >= 0.0:
        if c > 0.0:
            return (-b + math.sqrt(b * b + 4.0 * c * rhs)) / (2.0 * c)
        return rhs / b
    return -row_laminar(-q_prev, -drive, kv)


def row_of(name, q_prev, drive, kv):
    tau0 = TAU.get(name, 0.0)
    if tau0 > 0:
        return row_yield(q_prev, drive, kv, tau0)
    return row_laminar(q_prev, drive, kv)


def march(arm, nsteps, sched=None, npass=16):
    """One arm's whole trace from this file's copy of the three laws. A node with
    two free pressures (arm `m`) is the same loop with two node rows, so the mid
    row appears in both."""
    nodes, rows = ARMS[arm]
    q = {n: 0.0 for (n, _, _, _, _) in rows}
    p = {nd: P_NODE0 for nd in nodes}
    logs = {n: [0.0] for (n, _, _, _, _) in rows}
    flags = {n: [0] for (n, f, to, t, s) in rows if t > 0}
    plist = {nd: [p[nd]] for nd in nodes}

    ROW_ENDS = {n: (f, to) for (n, f, to, t, s) in rows}

    def ends(nm):
        f, to = ROW_ENDS[nm]
        pf = P_A if f == "A" else p[f]
        pt = P_B if to == "B" else p[to]
        return pf, pt

    for i in range(1, nsteps + 1):
        tn = i * DT
        # the entry latch: this step's decision, made from the PREVIOUS state's
        # rate and the previous state's node pressures (transient.cpp:310-:316)
        held = {}
        for (n, f, to, t, s) in rows:
            if t <= 0:
                continue
            pf, pt = ends(n)
            held[n] = (q[n] == 0.0 and abs(pf - pt) <= threshold(t))
        # fixed point on the node pressures; the sampled K is at t + dt
        trial = dict(p)
        for _ in range(npass):
            for nd in nodes:
                tot = 0.0
                for (n, f, to, t, s) in rows:
                    kv = valve_k(sched if sched else SCHEDULED.get(n), tn) if s else 0.0
                    pf = P_A if f == "A" else trial[f]
                    pt = P_B if to == "B" else trial[to]
                    qn = 0.0 if held.get(n) else row_of(n, q[n], pf - pt, kv)
                    if f == nd:
                        tot -= qn
                    if to == nd:
                        tot += qn
                trial[nd] = p[nd] + (DT / CC) * tot
        p = trial
        newq = {}
        for (n, f, to, t, s) in rows:
            kv = valve_k(SCHEDULED[n], tn) if s else 0.0
            pf = P_A if f == "A" else p[f]
            pt = P_B if to == "B" else p[to]
            newq[n] = 0.0 if held.get(n) else row_of(n, q[n], pf - pt, kv)
        for n in q:
            q[n] = newq[n]
        for (n, f, to, t, s) in rows:
            logs[n].append(newq[n])
        for n in flags:
            flags[n].append(1 if held[n] else 0)
        for nd in nodes:
            plist[nd].append(p[nd])
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
    flows, nodes, banner = {}, {}, None
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
        elif tok[:2] == ["steady", "solved:"]:
            banner = {"iter": int(tok[2].split("=")[1]),
                      "res": float(tok[3].split("=")[1].split()[0]),
                      "conv": "NOT" not in line,
                      "raw": line.replace("[chflow] ", "")}
    return flows, nodes, banner


def drive_at(name, main):
    """The drive the latch read on the step that made each decision, taken from
    the printed node pressure of the state BEFORE it, on the side this wall sits
    on."""
    nd = NODE_OF[name]
    ps = main["p"][nd]
    if SIDE[name] == "in":
        return [P_A - x for x in ps]
    return [x - P_B for x in ps]


def events(name, main):
    """Every latch flip AFTER the first step, with the drive the step started on
    and the drive the step before it started on. The first state is skipped on
    purpose: the app prints `static hold 0` at t = 0, before it has made a
    decision, so the 0 -> 1 that every held wall shows at step 1 is the first
    report of the latch, not an event the network caused."""
    flags = main["hold"][name]
    dv = drive_at(name, main)
    out = []
    for i in range(2, len(flags)):
        if flags[i] != flags[i - 1]:
            out.append({"step": i, "hold": flags[i], "entry": dv[i - 1],
                        "prev": dv[i - 2]})
    return out


proc = run(TRANSIENT, f"{CASE_DIR}/case.dict")
if proc.returncode != 0:
    raise SystemExit(f"transient solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
main = parse_transient(proc.stdout, PIPES, ALL_NODES)

proc_long = run(TRANSIENT, f"{CASE_DIR}/case_long.dict")
if proc_long.returncode != 0:
    raise SystemExit(f"long solver failed rc={proc_long.returncode}\n"
                     f"{proc_long.stdout}\n{proc_long.stderr}")
lp = [n for arm in LONG_ARMS for (n, _, _, _, _) in ARMS[arm][1]]
lnodes = sorted({x for arm in LONG_ARMS for x in ARMS[arm][0]})
longn = parse_transient(proc_long.stdout, lp, lnodes)

proc_st = run(STEADY, f"{CASE_DIR}/case_steady.dict")
if proc_st.returncode not in (0, 1):
    raise SystemExit(f"steady solver failed rc={proc_st.returncode}")
stq, stp, stbanner = parse_steady(proc_st.stdout)
proc_sd = run(STEADY, f"{CASE_DIR}/case_steady_dead.dict")
sdq, sdp, sdbanner = parse_steady(proc_sd.stdout)
proc_ss = run(STEADY, f"{CASE_DIR}/case_steady_slow.dict")
ssq, ssp, ssbanner = parse_steady(proc_ss.stdout)

EV = {n: events(n, main) for n in WALLS}
REL = {n: [e for e in EV[n] if e["hold"] == 0] for n in WALLS}
RL = {n: [e for e in EV[n] if e["hold"] == 1] for n in WALLS}
FLIPS = {n: EV[n] for n in WALLS}
NFLAG = sum(len(EV[n]) for n in WALLS)

# --- family 1: the roster, and that the kernel does not know which side ----
require(len(main["p"]["A"]) == NS and len(main["times"]) == NS,
        f"the 4 s file prints {NS} states on every row and every node of the "
        f"7-arm split network")
require(sorted(main["thr"]) == WALLS and len(WALLS) == 14,
        "14 yielded rows print their threshold and they are exactly the 14 walls "
        "of the case file -- 8 bringing flow to a free node, 6 taking it away")
require(sorted(set(main["thr"].values())) == [threshold(1.0), threshold(1.5),
                                              threshold(2.0)],
        "the thresholds printed are 4 tau_0 L / d = 4000 / 6000 / 8000 Pa and no "
        "fourth number appears, so inflow walls and drain walls share one line")
require(all(main["thr"][n] == THRESH[n] for n in WALLS),
        "each wall's printed threshold is 4 tau_0 L / d for its own tau_0, "
        "whatever side of the node it sits on")
require(all(p not in main["thr"] for p in NEWTONIAN),
        f"{len(NEWTONIAN)} Newtonian rows (drains, the mid row, the always-open "
        f"neighbour) print no threshold and no hold flag")
require(all(len(main["hold"][n]) == NS for n in WALLS),
        "every wall prints one hold flag per state, so the latch is reported on "
        "the outflow side exactly as on the inflow side")

# --- family 2: the two directions of the latch, counted -------------------
released = sorted(n for n in WALLS if REL[n])
relatched = sorted(n for n in WALLS if RL[n])
never_wakes = sorted(n for n in WALLS if main["hold"][n][1] == 1 and not REL[n])
asleep_at_first = sorted(n for n in WALLS if main["hold"][n][1] == 1)
awake_at_first = sorted(n for n in WALLS if main["hold"][n][1] == 0)
flowing_all_through = sorted(n for n in WALLS if sum(main["hold"][n]) == 0)
require(len(asleep_at_first) == 9 and len(awake_at_first) == 5,
        f"{len(asleep_at_first)} of the 14 walls are asleep at the first decision "
        f"({' '.join(asleep_at_first)}: all 8 walls that bring flow from A, whose "
        f"drive is 3000 Pa, plus the one 8000 Pa drain wall whose floor the "
        f"27000 Pa start does not reach), and {len(awake_at_first)} start awake "
        f"({' '.join(awake_at_first)}): the drain side of a node at 27000 Pa is "
        f"driven by 7000 Pa")
require(released == ["Sa", "Sb", "Sc", "Wa"],
        f"4 of the 9 sleeping walls wake ({len(sum(REL.values(), []))} wake events) "
        f"and all four are on the INFLOW side: {' '.join(released)}")
require(relatched == ["Sa", "Wb", "Wz"],
        f"3 walls are put back to sleep ({len(sum(RL.values(), []))} re-latch "
        f"events): Sa by a weaker inflow wall as in 5.54, Wz and Wb by the "
        f"mirror mechanism on the DRAIN side -- and both of those two started "
        f"awake, so the drain side of this file only ever closes the latch: "
        f"{' '.join(relatched)}")
require(NFLAG == 8 and max(len(EV[n]) for n in WALLS) == 3,
        f"{NFLAG} latch flips in the whole file: "
        f"{len(sum(REL.values(), []))} wake events and "
        f"{len(sum(RL.values(), []))} re-latches; only Sa flips more than once "
        f"(3 times: {' '.join(str(e['step']) for e in EV['Sa'])})")
require(never_wakes == ["Ed", "Sd", "Sf", "Sm", "Sz"],
        f"{len(never_wakes)} walls go to sleep at step 1 and never leave: "
        f"{' '.join(never_wakes)} -- the 5 whose own limit the node never crosses "
        f"(family 4 prints every distance), each printing a single token")
require(flowing_all_through == ["Ef", "Em", "Tz"],
        f"3 drain walls flow in every one of the {NS} states and never get to "
        f"report a held state: {' '.join(flowing_all_through)}")

# --- family 3: both directions read the same number ----------------------
gaps_over = [(n, e["entry"] - THRESH[n]) for n in WALLS for e in REL[n]]
gaps_under = [(n, THRESH[n] - e["entry"]) for n in WALLS for e in RL[n]]
require(gaps_over and all(0.0 < g < 20.0 for _, g in gaps_over),
        f"every wake event starts from a drive {min(g for _, g in gaps_over):.2f} - "
        f"{max(g for _, g in gaps_over):.2f} Pa OVER its own 4 tau_0 L / d, the "
        f"same one-step crossing 5.54 measured on the inflow side")
rl_onestep = [(n, g) for n, g in gaps_under if g < 20.0]
require(len(rl_onestep) == 2 and sorted(n for n, _ in rl_onestep) == ["Sa", "Wz"],
        f"two of the three re-latches are one-step crossings too: "
        f"{' '.join(f'{n} {g:.2f}' for n, g in rl_onestep)} Pa UNDER the same "
        f"number the release read OVER")
require(abs(THRESH["Wb"] - EV["Wb"][0]["entry"]) > 50.0,
        f"Wb's re-latch is not a one-step crossing: the drive at the step that "
        f"holds it is {THRESH['Wb'] - EV['Wb'][0]['entry']:.2f} Pa UNDER its "
        f"threshold, because the Newtonian drain keeps decelerating the column "
        f"after the drive has already gone inside the wall")
require(all(abs(e["prev"] - THRESH[n]) < 120.0 for n in WALLS for e in EV[n]),
        "the step before every flip is on the other side of the same number "
        "(|drive - 4 tau_0 L / d| < 120 Pa on all 8 flips), so there is no "
        "hysteresis band to fit")
# the law: which wall a network is even allowed to put back to sleep. Predicted
# from ARMS and the thresholds alone -- no trace is read -- and then compared with
# the three events family 2 counted.
ROWS_AT = {}
for arm in ARMS:
    for (n, f, to, t, _) in ARMS[arm][1]:
        for nd in (f, to):
            if nd.startswith("F"):
                ROWS_AT.setdefault(nd, []).append(
                    (n, t, "supply" if to == nd else "drain"))


def carrier(n):
    """The rows that can carry n's node past n's own limit. Two conditions, and
    they are the same word on both sides of the node: the row has to push the node
    the way n's limit lies -- a supplier of the node can only raise it, a drainer
    can only lower it -- and its own limit has to be FURTHER away than n's, which
    for the ceiling `pA - 4 tau_0 L / d` and the floor `pB + 4 tau_0 L / d` is a
    SMALLER threshold either way. A Newtonian row has threshold 0, so it always
    qualifies; two walls of equal threshold never carry each other, which is
    5.54's arm `i` and this file's arm `d`. A wall sitting alone on its side of
    its node has no carrier whatever else the node is connected to."""
    nd = NODE_OF[n]
    mine = [d for (m, t, d) in ROWS_AT[nd] if m == n][0]
    return [m for (m, t, d) in ROWS_AT[nd]
            if m != n and d == mine and t < TAU[n]]


predicted = sorted(n for n in WALLS if carrier(n))
require(predicted == relatched,
        f"the ceiling-or-floor law, read before the trace, names exactly the "
        f"{len(predicted)} walls that go back to sleep: "
        f"{' '.join(f'{n}<-{carrier(n)[0]}' for n in predicted)}")
require(all(len(carrier(n)) == 1 for n in predicted) and
        all(TAU.get(carrier(n)[0], 0.0) < TAU[n] for n in predicted),
        "and each of the three has one carrier, not a choice: Sa by Wa (6000 "
        "under 8000), Wz by Tz (4000 under 6000), Wb by the Newtonian drain Db "
        "(no threshold at all), so the mechanism is the same event mirrored")
require(sorted(set(WALLS) - set(predicted)) ==
        sorted(never_wakes + flowing_all_through + ["Sb", "Sc", "Wa"]),
        f"the other {len(WALLS) - len(predicted)} walls have no carrier on their "
        f"own side of their node and none of them re-latches")

# --- family 4: the floor, measured from both sides -----------------------
z_floor = P_B + THRESH["Tz"]
f_floor = P_B + THRESH["Ef"]
gap_f = min(main["p"]["Ff"]) - f_floor
gap_z = min(main["p"]["Fz"]) - z_floor
require(0.0 < gap_f < 400.0,
        f"arm f: nobody crosses the 6000 drain's floor, so its node falls to "
        f"{f_floor:.1f} Pa and stops {gap_f:.2f} Pa above it -- the rate never "
        f"reaches zero, so the latch is never asked")
require(0.0 < gap_z < 400.0,
        f"arm z: the node is left {gap_z:.2f} Pa above the WEAKEST drain's floor "
        f"{z_floor:.1f} Pa, not above the 6000 drain's {f_floor:.1f} Pa, which is "
        f"what let Wz be re-latched at step {RL['Wz'][0]['step']}")
require(min(main["p"]["Fz"]) < P_B + THRESH["Wz"] and
        min(main["p"]["Fb"]) < P_B + THRESH["Wb"],
        f"the two re-latched drain walls both had their floor crossed: Fz to "
        f"{min(main['p']['Fz']):.2f} against {P_B + THRESH['Wz']:.1f} and Fb to "
        f"{min(main['p']['Fb']):.2f} against {P_B + THRESH['Wb']:.1f}")
require(max(drive_at("Sz", main)) < THRESH["Sz"],
        f"arm z's inflow wall is never even asked: its drive reaches at most "
        f"{max(drive_at('Sz', main)):.2f} Pa against a threshold of "
        f"{THRESH['Sz']:.1f}, so it stays held while two drains empty the node")
require(max(drive_at("Sf", main)) < THRESH["Sf"] and
        max(drive_at("Sm", main)) < THRESH["Sm"] and
        max(drive_at("Sd", main)) < THRESH["Sd"],
        "the same for Sf, Sm and Sd: the node never gets low enough for the "
        "inflow wall's own release to be on offer, so a chain of walls can leave "
        "one of its members permanently asleep")

# --- family 5: the dead band is a fixed point of the WHOLE solve --------
require(set(main["raw"]["Sd"]) == {"0"} and set(main["raw"]["Ed"]) == {"0"},
        "arm d (8000 in, 8000 out) prints the single token 0 on BOTH rows in all "
        "2001 states -- the two walls' flow windows do not overlap")
require(len(set(main["p"]["Fd"])) == 1 and main["p"]["Fd"][0] == P_NODE0,
        f"and its node prints one token too ({P_NODE0:.1f} Pa): with no row "
        f"moving, the node row is 0 = 0 and the march cannot leave the state it "
        f"started in")
require(sum(main["hold"]["Sd"]) == NS - 1 and sum(main["hold"]["Ed"]) == NS - 1,
        f"both flags read 1 from the second state on ({NS - 1} of {NS}): the "
        f"latch condition is satisfied by the state the latch itself freezes")
require(min(drive_at("Ed", main)) == max(drive_at("Ed", main)) and
        max(drive_at("Ed", main)) < THRESH["Ed"],
        f"the drain wall of the dead band sits "
        f"{THRESH['Ed'] - max(drive_at('Ed', main)):.2f} Pa inside its own "
        f"threshold and its drive never changes there ({P_NODE0 - P_B:.0f} Pa "
        f"against {THRESH['Ed']:.0f}): the node it hangs on is the one node the "
        f"file never moves, so nothing in it is ever asked")

# --- family 6: silence, and what a split costs --------------------------
require(set(main["raw"]["Sz"]) == {"0"} and set(main["raw"]["Sm"]) == {"0"} and
        set(main["raw"]["Sf"]) == {"0"},
        "the three never-asked inflow walls print one token each, the same "
        "silence 5.54 measured for a held row on the inflow side")
n_open = int(0.8 / DT)
same_before = sum(1 for a, b in zip(main["raw"]["Dc"][:n_open + 1],
                                   main["raw"]["Nc"][:n_open + 1]) if a == b)
require(same_before == n_open + 1,
        f"arm c's two geometrically identical drains print token-identical rates "
        f"for all {n_open + 1} states before the schedule moves ({n_open + 1} of "
        f"{n_open + 1}) -- the split ratio is 1 and it is not a fit")
diff_after = [i for i in range(n_open + 1, NS)
              if main["raw"]["Dc"][i] != main["raw"]["Nc"][i]]
require(len(diff_after) == NS - n_open - 1 and diff_after[0] == n_open + 1,
        f"they differ in EVERY state after the schedule starts ({len(diff_after)} "
        f"of the {NS - n_open - 1} states after t={CLOSE_STD[0]} s, first at "
        f"t={diff_after[0] * DT:.3f} s), which is the same "
        f"K-is-in-the-row-but-not-in-the-drive asymmetry 5.52 found")
i_shut = int(CLOSE_STD[1] / DT)
ratio_shut = main["q"]["Dc"][i_shut] / main["q"]["Nc"][i_shut]
require(ratio_shut < 0.1 and main["q"]["Dc"][i_shut] > 0.0,
        f"and the split ratio at closedK is {ratio_shut:.4f}, not zero: the row "
        f"the timetable shut still carries {main['q']['Dc'][i_shut]:.4e} against "
        f"its open twin's {main['q']['Nc'][i_shut]:.4e} m3/s, because K sits in "
        f"the momentum row and the row has to pay the inertia it already has")
i_re = int(round((CLOSE_STD[3] + DT) / DT))
rel_after = [rel(main["q"]["Dc"][i], main["q"]["Nc"][i]) for i in range(i_re, NS)]
require(all(rel_after[i] >= rel_after[i + 1] for i in range(len(rel_after) - 1)),
        f"from the first step of the reopen ({i_re * DT:.3f} s) the gap only "
        f"shrinks, state by state ({rel_after[0]:.3f} -> {rel_after[-1]:.3e} "
        f"relative over {len(rel_after)} states), and at the last one the twins "
        f"agree to {main['q']['Dc'][-1]:.6e} / {main['q']['Nc'][-1]:.6e} -- what "
        f"is left of the difference is the column's memory, not the schedule")
require(len(RL["Sc"] if "Sc" in RL else []) == 0 and len(REL["Sc"]) == 1,
        f"arm c's wall releases once at step {REL['Sc'][0]['step']} and is never "
        f"re-latched although its node has degree 3: the two extra rows are on "
        f"the DRAIN side and Sc is the node's only supplier, so the same-side "
        f"law of family 3 gives it no carrier -- what buys a second direction is "
        f"a weaker row on the wall's OWN side (arm b has one and re-latches), "
        f"not a large degree")

# --- family 7: the 4 s trace reproduced by this file's own rows ---------
worst_q = 0.0
worst_p = 0.0
flag_bad = []
step_bad = []
for arm in ARMS:
    logs, flags, plist = march(arm, NS - 1)
    nodes = ARMS[arm][0]
    for (n, _, _, _, _) in ARMS[arm][1]:
        got = main["q"][n]
        for i in range(NS):
            r = rel(logs[n][i], got[i])
            if logs[n][i] != got[i] and r > worst_q and max(abs(got[i]), abs(logs[n][i])) >= RATE_FLOOR:
                worst_q = r
        if got and got[0] == 0.0 and logs[n][0] != 0.0:
            raise SystemExit("first state is not zero")
    for nd in nodes:
        got = main["p"][nd]
        for i in range(NS):
            worst_p = max(worst_p, rel(plist[nd][i], got[i]))
    for n in flags:
        if flags[n] != main["hold"][n]:
            bad = [i for i in range(NS) if flags[n][i] != main["hold"][n][i]]
            flag_bad.append((n, len(bad), bad[:3]))
        ev = [i for i in range(2, NS) if flags[n][i] != flags[n][i - 1]]
        if ev != [e["step"] for e in EV[n]]:
            step_bad.append((n, ev, [e["step"] for e in EV[n]]))
require(worst_q < 2e-9,
        f"{sum(len(ARMS[arm][1]) for arm in ARMS)} x {NS} = "
        f"{sum(len(ARMS[arm][1]) for arm in ARMS) * NS} printed rates are "
        f"reproduced by this file's own three rows to {worst_q:.3e} relative on "
        f"every state whose rate is at least {RATE_FLOOR:.0e} m3/s")
require(worst_p < 2e-11,
        f"and the {len(ALL_NODES) - 2} free node pressures -- including the two "
        f"of arm m's chain -- to {worst_p:.3e} relative")
require(not flag_bad, f"zero latch-flag mismatches over {14 * NS} printed flags "
                      f"({flag_bad})")
require(not step_bad, f"and the flip step list is reproduced state by state for "
                      f"all 14 walls ({step_bad})")

# --- family 8: continuity from the printed tokens, both branch types ----
worst_abs = 0.0
worst_rel = 0.0
bal_where = None
rel_where = None
QUANTUM = 1e-7 * CC / DT              # one unit of the last printed pressure digit
for arm in ARMS:
    for nd in ARMS[arm][0]:
        for i in range(1, NS):
            tot = 0.0
            for (n, f, to, t, s) in ARMS[arm][1]:
                v = main["q"][n][i]
                if to == nd:
                    tot += v
                if f == nd:
                    tot -= v
            moved = (main["p"][nd][i] - main["p"][nd][i - 1]) * CC / DT
            scale = max(abs(tot), abs(moved), 1e-30)
            if abs(tot - moved) > worst_abs:
                worst_abs = abs(tot - moved)
                bal_where = (arm, nd, i, tot, moved)
            if abs(tot - moved) / scale > worst_rel:
                worst_rel = abs(tot - moved) / scale
                rel_where = (arm, nd, i, tot, moved)
require(worst_abs <= QUANTUM,
        f"all {len(ALL_NODES) - 2} node rows x {NS - 1} steps close on the "
        f"printed tokens alone to {worst_abs:.3e} m3/s, which is one unit of the "
        f"last printed pressure digit ({QUANTUM:.1e} = 1e-07 Pa x V/BC/dt); the "
        f"largest term anywhere near it is {abs(bal_where[3]):.4e} at "
        f"{bal_where[1]} step {bal_where[2]}")
require(worst_rel < 2e-3,
        f"the relative mismatch only grows to {worst_rel:.3e}, at "
        f"{rel_where[1]} step {rel_where[2]}, where both sides of the row are "
        f"{abs(rel_where[3]):.1e} and {abs(rel_where[4]):.1e} m3/s -- four steps "
        f"into a node whose whole trace lives in the last digits the print has "
        f"room for; a held wall enters the row as the zero it prints and Fd, the "
        f"node nothing moves, closes on 0 = 0")

# --- family 9: the calendar question (18 s) ----------------------------
require([e["step"] for e in events("Wb", longn)] == [83] and
        [e["step"] for e in events("Wz", longn)] == [193] and
        [e["step"] for e in events("Sb", longn)] == [317],
        "the long run's three flips land on the same steps as the 4 s run "
        "(Wb 83, Wz 193, Sb 317) even though its reopen is 12 s later")
require(sum(longn["hold"]["Wb"]) == NL - 83 and sum(longn["hold"]["Wz"]) == NL - 193,
        f"and they stay asleep: Wb is held in {sum(longn['hold']['Wb'])} of "
        f"{NL} states and Wz in {sum(longn['hold']['Wz'])}")
i_open = int(CLOSE_LONG[3] / DT)
require(set(longn["hold"]["Wb"][i_open:]) == {1},
        f"reopening the drain at {CLOSE_LONG[3]} s does not wake the wall it put "
        f"to sleep ({i_open} to {NL - 1} all read 1) -- the seal pulls the node "
        f"under a floor, the reopen only lets the node fall further")
gaps_ef = [(s, longn["p"]["Ff"][s] - f_floor) for s in (1000, 3000, 5000, 7000, 9000)]
gaps_tz = [(s, longn["p"]["Fz"][s] - z_floor) for s in (1000, 3000, 5000, 7000, 9000)]
require(all(g > 0 for _, g in gaps_ef + gaps_tz) and
        all(b[1] < a[1] for a, b in zip(gaps_ef, gaps_ef[1:])) and
        all(b[1] < a[1] for a, b in zip(gaps_tz, gaps_tz[1:])),
        f"the two floors are approached and never reached: Ff {gaps_ef[0][1]:.2f} "
        f"-> {gaps_ef[-1][1]:.2f} Pa and Fz {gaps_tz[0][1]:.2f} -> "
        f"{gaps_tz[-1][1]:.2f} Pa over steps 1000 - 9000")
shr_per_k = [(gaps_tz[j + 2][1] / gaps_tz[j][1]) ** 0.5 for j in (0, 1, 2)]
require(all(s < 1.0 for s in shr_per_k),
        f"and arm z's gap shrinks by {shr_per_k[0]:.5f} -> {shr_per_k[-1]:.5f} "
        f"per thousand steps, below one everywhere: 18 s of clock turns "
        f"'not yet' into 'not reachable', the same verdict 5.54 read on the "
        f"ceiling side")
require(min(abs(x) for x in longn["q"]["Ef"] if x != 0) > 0 and
        sum(longn["hold"]["Ef"]) == 0 and sum(longn["hold"]["Tz"]) == 0,
        f"the two walls that were never allowed to cross never report a held "
        f"state in {NL - 1} marched steps either")
shared = 0
maxdiff = 0.0
for nd in ("Ff", "Fz"):
    w = min(NS, NL)
    for i in range(w):
        shared += 1
        maxdiff = max(maxdiff, abs(main["p"][nd][i] - longn["p"][nd][i]))
require(maxdiff < 5e-7,
        f"the two arms with no schedule in them agree between the runs over "
        f"{shared} states to {maxdiff:.3e} Pa, not token for token: the 19-pipe "
        f"and the 8-pipe systems sum different dense rows")
w_open = int(CLOSE_STD[0] / DT) + 1
d_fb = max(abs(main["p"]["Fb"][i] - longn["p"]["Fb"][i]) for i in range(w_open))
require(d_fb < 5e-7,
        f"arm b agrees over the window where BOTH files' schedules are the same "
        f"({w_open} states, {d_fb:.3e} Pa) and diverges after it by up to "
        f"{max(abs(main['p']['Fb'][i] - longn['p']['Fb'][i]) for i in range(NS)):.1f} Pa, "
        f"which is the reopen the long file has not performed yet")

# --- family 10: the steady kernel, three budgets and two failure kinds --
require(sdbanner["conv"] and sdbanner["iter"] == 0 and sdbanner["res"] == 0.0,
        f"the dead-band arms alone are an EXACT steady solution: the banner is "
        f"'{sdbanner['raw']}' -- Newton stops at its first look because 0 = 0")
DEAD_ROWS = ("Sd", "Ed", "Sf", "Ef", "Sm", "Mm", "Em")
require(all(sdq[n]["q"] == 0.0 for n in DEAD_ROWS),
        f"and every one of its {len(DEAD_ROWS)} rows ({' '.join(DEAD_ROWS)}) "
        f"reports no rate: the arm `d` pair, the arm `f` pair and the whole "
        f"two-node chain of arm `m`")
require(all((sdq[n]["raw"], sdq[n]["dp"]) == (stq[n]["raw"], stq[n]["dp"])
            for n in DEAD_ROWS) and
        all(sdp[nd] == stp[nd] == 0.5 * (P_A + P_B)
            for nd in ("Fd", "Ff", "Fm1", "Fm2")),
        f"the rank deficiency is LOCAL: those {len(DEAD_ROWS)} rows print "
        f"token-identical rate AND drop in the 7-arm file and the 3-arm file "
        f"(Sd/Ed/Sf/Ef/Sm/Em at dp 5000 = the midpoint the guess puts them, Mm at "
        f"0), and their 4 nodes print the same guess in both, so what breaks "
        f"Newton in the full file is the two arms that DO flow (a and z), not the "
        f"three that print nothing")
require(all(sdp[n] == 0.5 * (P_A + P_B) for n in ("Fd", "Ff", "Fm1", "Fm2")),
        f"the undetermined nodes print the solver's own guess "
        f"{0.5 * (P_A + P_B):.1f} Pa = (pA + pB) / 2 -- steady has no way to say "
        f"which element of the dead band the march parked in")
require(main["p"]["Fd"][-1] != sdp["Fd"] and main["p"]["Ff"][-1] != sdp["Ff"],
        f"and it is not the same number: the transient left Fd at "
        f"{main['p']['Fd'][-1]:.1f} and Ff at {main['p']['Ff'][-1]:.2f} Pa, both "
        f"inside the band steady prints as one point")
require(not stbanner["conv"], f"the full seven-arm steady file does NOT converge "
                              f"at maxIter 200: '{stbanner['raw']}'")
imb_a = stq["Da"]["q"] - stq["Sa"]["q"] - stq["Wa"]["q"]
require(abs(stbanner["res"] * SOURCE_SCALE - imb_a) < 1e-9 * abs(imb_a) + 1e-12,
        f"and the residual it stops at IS arm a's node imbalance: {imb_a:.6e} "
        f"m3/s against the banner's {stbanner['res']:.10f} x the sources' scale "
        f"{SOURCE_SCALE:.0e} = {stbanner['res'] * SOURCE_SCALE:.6e}")
imb_z = stq["Tz"]["q"]
require(abs(stq["Sz"]["q"]) == 0.0 and abs(stq["Wz"]["q"]) == 0.0,
        f"arm z's two stronger walls report no rate in steady, so the only flow "
        f"at that node is the 4000 drain's own {imb_z:.4e} m3/s, which continuity "
        f"refuses")
require(abs(ssp["Fa"] - 21379.7317) < 1e-3 and ssbanner["res"] < stbanner["res"],
        f"at maxIter 50000 the same file creeps to Fa = {ssp['Fa']:.7f} Pa -- the "
        f"landing point 5.54 published for this arm -- and the residual falls to "
        f"{ssbanner['res']:.4e}, so this is a rate of convergence and not a "
        f"refusal")
require(abs(stp["Fa"] - ssp["Fa"]) > 100.0 and
        abs(stp["Fb"] - ssp["Fb"]) < 1e-9 and abs(stp["Fc"] - ssp["Fc"]) < 1e-9,
        f"the two stalling arms move with the budget (Fa "
        f"{stp['Fa']:.4f} -> {ssp['Fa']:.4f}) while the two arms that converged "
        f"print the SAME token at both budgets (Fb {stp['Fb']:.4f}, "
        f"Fc {stp['Fc']:.4f})")
require(stp["Fb"] > P_B and abs(ssp["Fz"] - z_floor) < 1.0,
        f"steady's own approach to the floor is the same limit the march fell "
        f"toward: Fz settles {ssp['Fz']:.6f} Pa, {ssp['Fz'] - z_floor:.4f} Pa "
        f"above the 4000 drain's floor")
asleep_st = sorted(n for n in WALLS if n in stq and stq[n]["q"] == 0.0)
march_asleep = sorted(n for n in WALLS if main["hold"][n][-1] == 1)
require(set(march_asleep) <= set(asleep_st),
        f"steady refuses every one of the {len(march_asleep)} walls the march "
        f"leaves asleep at its last state ({' '.join(march_asleep)})")
steady_extra = sorted(set(asleep_st) - set(march_asleep))
require(asleep_st == sorted(march_asleep + steady_extra) and
        steady_extra == ["Ef", "Em"],
        f"steady reports {len(asleep_st)} walls at exactly zero, {len(steady_extra)} "
        f"more than the march: {' '.join(steady_extra)} -- the drain wall of each "
        f"of the two arms whose node has no steady balance at all (supply ceiling "
        f"22000 below drain floor 26000), which the march nevertheless keeps "
        f"flowing because it only ever approaches that floor")
require(all(main["hold"][n][-1] == 0 and drive_at(n, main)[-1] > THRESH[n]
            for n in steady_extra),
        f"and both of the two are flowing hard at the last marched state: Ef at "
        f"{drive_at('Ef', main)[-1]:.2f} and Em at {drive_at('Em', main)[-1]:.2f} "
        f"Pa of drive against a 6000 Pa threshold, awake in the march and asleep "
        f"in steady -- the two files answer two different questions")
require(abs(stq["Sf"]["q"]) == 0.0 and abs(stq["Sm"]["q"]) == 0.0,
        f"the two supply walls of those arms report no rate in steady either, so "
        f"arms f and m are dead bands for the steady kernel while the march runs "
        f"them as ordinary flowing arms")
require(abs(stq["Db"]["q"] - stq["Sb"]["q"]) < 1e-12 and stq["Wb"]["q"] == 0.0,
        f"arm b's settled state is the SERIES pair Sb -> node -> Db with the "
        f"drain wall out of it: {stq['Sb']['q']:.11e} on both rows, because a "
        f"wall asleep below its floor contributes exactly what the transient "
        f"prints for it")
require(abs(stq["Dc"]["q"] - stq["Nc"]["q"]) < 1e-15,
        f"steady gives arm c's two identical drains the same rate to the last "
        f"printed digit ({stq['Dc']['raw']} / {stq['Nc']['raw']}), which is the "
        f"split ratio the transient also reports once the valve is open")

# --- family 11: the tail, the smallest rates, the arithmetic of sleep ----
sa_steps = [e["step"] for e in EV["Sa"]]
require(sa_steps == [333, 507, 1530],
        f"arm a's strong wall wakes, sleeps and wakes again at {sa_steps} -- the "
        f"same three steps 5.54 published for the same three-pipe arm, reached "
        f"here through a 19-pipe network")
require(abs(main["p"]["Fa"][-1] - 21375.816193) < 1e-4,
        f"and its last marched node pressure is {main['p']['Fa'][-1]:.6f} Pa, "
        f"the number 5.54 printed as this arm's tail")
held_states = sum(sum(main["hold"][n]) for n in WALLS)
tok_zero = sum(1 for n in WALLS for x in main["raw"][n] if x == "0")
quiet = [(n, i) for n in WALLS for i, x in enumerate(main["raw"][n])
         if x == "0" and main["hold"][n][i] == 0 and i > 0]
require(tok_zero - held_states == len(quiet) + 14,
        f"{held_states} printed held states and {tok_zero} printed zero tokens "
        f"over 14 walls: the 14 rows of the t=0 state print 0 with flag 0, before "
        f"the solver has decided anything, and the remaining {len(quiet)} are "
        f"reported next")
require(sorted(quiet) == [("Sa", 506), ("Wb", 82), ("Wz", 192)] and
        all([e["step"] for e in RL[n]][0] == i + 1 for n, i in quiet),
        f"the {len(quiet)} zero rates printed with the flag still open are the "
        f"clamp steps: " + " ".join(f"{n} at {i}" for n, i in quiet) +
        ", and each is exactly the state before that wall's re-latch. Q = 0 with "
        "flag 0 is what makes `qPrev == 0` true on the next step, so the latch is "
        "not a guess about the drive, it is a consequence of the rate the row "
        "just reported")
nz = [(n, min(abs(x) for x in main["q"][n] if x != 0)) for n in WALLS
      if any(x != 0 for x in main["q"][n])]
smallest = min(nz, key=lambda t: t[1])
i_sm = next(i for i, x in enumerate(main["q"][smallest[0]])
            if abs(x) == smallest[1])
require(smallest[0] == "Sa" and smallest[1] < 1e-9 and i_sm == 505,
        f"the smallest non-zero rate any wall ever prints is {smallest[1]:.4e} "
        f"m3/s on {smallest[0]} at step {i_sm} -- the step before the clamp "
        f"carries it to exactly zero, and the row is a 12-digit print of a "
        f"momentum balance whose inertial term is rho L / A / dt * dq, so a rate "
        f"that small is the row saying it has run out of column to move")

print()
print(f"{len(sys.argv)} argv, {len(ARMS)} arms, {len(WALLS)} walls, {len(PIPES)} "
      f"pipes, {NS} + {NL} states, {NFLAG} latch flips "
      f"({len(sum(REL.values(), []))} wakes, {len(sum(RL.values(), []))} "
      f"re-latches) verified against 3 steady budgets in {CHECKS} checks")
if failures:
    print(f"{len(failures)} FAILURES")
    for f in failures:
        print("  -", f)
    sys.exit(1)
print("all checks passed")
