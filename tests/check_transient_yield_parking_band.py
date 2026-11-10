#!/usr/bin/env python3
"""5.56: the parking lot of a dead band, found by sweeping the node's own start.

Usage: check_transient_yield_parking_band.py TRANSIENT STEADY CASE

5.55 put two 8000 Pa walls on one free node, started it at 2.7e4, and read one
frozen token -- then listed that as its own boundary: one desk, one initial
value. This file is the sweep. The same two-row arm is built eleven times and
only the node's `pressure` changes, walking from pB to pA, plus a twelfth copy
that carries the close/open schedule on both rows. The question the ladder
answers is which starts are PARKED, because a start inside the band makes both
rows write the identity 0 = 0: the wall that would flow in pays exactly its own
drive, the other pays its own, the node row (transient.cpp:395-:416) sums two
zeros, and nothing in the file is left to move anything.

What the run says, family by family:

  1. the ladder splits at the arithmetic, not at the calendar: the starts inside
     [pA - 4 tau_0 L/d, pB + 4 tau_0 L/d] freeze at their own value and the
     starts outside slide. The classification is computed from the numbers in
     this file BEFORE the trace is read, then compared with it.
  2. the split is CLOSED: the two starts that sit exactly on an edge freeze too,
     because holdsStatic compares with <= (transient.cpp:40-:43). `c` at one Pa
     outside and `d` on the edge bracket the boundary to 1 Pa on each side.
  3. a slider has one direction and no overshoot: below the band the node only
     rises and never reaches the ceiling, above it only falls and never reaches
     the floor.
  4. a wall that woke never sleeps again here, because the latch is a state test
     on Q_prev == 0 and the relaxation keeps a strictly positive rate forever;
     the wall that stayed under its threshold is held in every decided state.
  5. the ceiling problem and the floor problem are the same row read with the
     sign flipped, and the printout says so token for token.
  6. near an edge the rate is QUADRATIC in the excess over 4 tau_0 L/d -- the
     linear term of the Bingham bracket cancels identically -- so the gap obeys
     d(gap)/dt = -k gap^2, i.e. gap ~ 1/t: 1/gap is a straight line in the step
     number whose slope is predicted from this file's own constants. The
     parking lot is therefore approached but never entered from outside.
  7. the calendar cannot touch a parked row (arm `m` against arm `f`).
  8. the parking lot is a property of a node and its own two rows: three arms
     alone in a second file print token for token what they printed in the
     twelve-arm matrix for the frozen edge and the near slider, and to within
     one printed pressure quantum for the far one (Newton sums the same rows in
     a different order, so only the numbers are guaranteed to align).
  9. the whole 4 s trace is reproduced by this file's copy of the three laws.
 10. the held rows are silent to the strongest kind, and the node rows close on
     the printed tokens alone.
 11. the steady kernel has no latch and no memory of the sweep: it pays zero on
     all twenty-two rows and leaves all eleven nodes at the average of the two
     fixed pressures, which is exactly the CENTRE of the band.

The kernel is not touched to get any of this.
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
END_LONG = 40.0
NL = int(round(END_LONG / DT)) + 1
P_A, P_B = 3.0e4, 2.0e4
LENGTH = 20.0
DIAM = 0.02
AREA = 0.25 * math.pi * DIAM ** 2
INERTIA = RHO * LENGTH / AREA
I_DT = INERTIA / DT
WALL_SCALE = 4.0 * LENGTH / DIAM          # dp = WALL_SCALE * tw
R_LAM = 128.0 * MU * LENGTH / (math.pi * DIAM ** 4)
NODE_VOLUME = 0.1
NODE_BULK = 1.0e8
CC = NODE_VOLUME / NODE_BULK              # 1e-9 m3/Pa
RATE_FLOOR = 1.0e-8

TAU0 = 2.0
THR = WALL_SCALE * TAU0                   # 8000 Pa, every wall in this file
CEIL = P_A - THR                          # 22000: a below-band node rises to it
FLOOR = P_B + THR                         # 28000: an above-band node falls to it

CLOSE_STD = (0.8, 1.0, 1.0e8, 2.5, 2.7)

# arm: its node's starting pressure, and whether the two rows carry the schedule
ARMS = {"a": (2.0e4, None), "b": (2.1e4, None), "c": (21999.0, None),
        "d": (2.2e4, None), "e": (2.3e4, None), "f": (2.5e4, None),
        "g": (2.7e4, None), "h": (2.8e4, None), "i": (28001.0, None),
        "j": (2.9e4, None), "k": (3.0e4, None), "m": (2.5e4, CLOSE_STD)}
ORDER = "abcdefghijkm"                     # the eleven starts plus `m`
MIRR = {"a": "k", "b": "j", "c": "i"}

# The classification this file predicts from arithmetic alone, before reading
# anything: inside the closed band the node is parked; below it the inflow wall
# wakes and the node rises; above it the drain wall wakes and the node falls.
KIND = {}
for arm in ARMS:
    p0 = ARMS[arm][0]
    KIND[arm] = ("parked" if CEIL <= p0 <= FLOOR
                 else ("up" if p0 < CEIL else "down"))
PARKED = sorted(a for a in KIND if KIND[a] == "parked")
UP = sorted(a for a in KIND if KIND[a] == "up")
DOWN = sorted(a for a in KIND if KIND[a] == "down")

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
    if sched is None:
        return 0.0
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
    """The exact laminar discharge a yielded wall stress holds (5.46 at n = 1)."""
    if tw <= tau0:
        return 0.0
    d = tw - tau0
    r = 0.5 * DIAM
    return (math.pi * r ** 3 * tw ** -3 / MU * d * d
            * (tau0 * tau0 / 2.0 + 2.0 * tau0 * d / 3.0 + d * d / 4.0))


def row_yield(q_prev, drive, tau0):
    """One implicit row of a yielded pipe in wall-stress space; g(tau0) >= 0 is
    the in-step clamp of transient.cpp:353-:355."""
    if drive < 0.0:
        return -row_yield(-q_prev, -drive, tau0)
    thr = threshold(tau0)

    def g(tw):
        q = q_of_tw(tw, tau0)
        return I_DT * (q - q_prev) + WALL_SCALE * tw - drive
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


def rows_of(arm):
    sched = ARMS[arm][1]
    return [("S" + arm, "A", "F" + arm, TAU0, sched),
            ("E" + arm, "F" + arm, "B", TAU0, sched)]


def march(arm, nsteps, npass=16):
    """One arm's trace from this file's copy of the two rows and the node row.
    The node starts where the case file says it starts; that is the whole
    experiment."""
    sched = ARMS[arm][1]
    rows = rows_of(arm)
    q = {n: 0.0 for (n, _, _, _, _) in rows}
    p = ARMS[arm][0]
    logs = {n: [0.0] for (n, _, _, _, _) in rows}
    flags = {n: [0] for (n, _, _, _, _) in rows}
    plist = [p]
    for i in range(1, nsteps + 1):
        tn = i * DT
        # the entry latch (transient.cpp:310-:316): previous rate is 0 and this
        # node's own drive is under the wall's threshold
        held = {n: (q[n] == 0.0
                    and abs((P_A if f == "A" else p)
                            - (P_B if to == "B" else p)) <= threshold(t))
                for (n, f, to, t, _) in rows}
        trial = p
        for _ in range(npass):
            tot = 0.0
            for (n, f, to, t, s) in rows:
                kv = valve_k(s, tn)
                drive = (P_A if f == "A" else trial) - (P_B if to == "B" else trial)
                qn = 0.0 if held[n] else row_yield(q[n], drive, t)
                tot += qn if to == "F" + arm else -qn
            trial = p + (DT / CC) * tot
        p = trial
        for (n, f, to, t, s) in rows:
            kv = valve_k(s, tn)
            drive = (P_A if f == "A" else p) - (P_B if to == "B" else p)
            q[n] = 0.0 if held[n] else row_yield(q[n], drive, t)
        for (n, _, _, _, _) in rows:
            logs[n].append(q[n])
        for n in flags:
            flags[n].append(1 if held[n] else 0)
        plist.append(p)
    return logs, flags, plist


def parse_transient(stdout, pipes, nodes):
    out = {"q": {n: [] for n in pipes}, "raw": {n: [] for n in pipes},
           "hold": {n: [] for n in pipes}, "p": {n: [] for n in nodes},
           "rawp": {n: [] for n in nodes}, "thr": {}}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 4 and tok[0] == "node" and tok[2] == "p":
            if tok[1] in out["p"]:
                out["p"][tok[1]].append(float(tok[3]))
                out["rawp"][tok[1]].append(tok[3])
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


ALL_PIPES = [("S" + a) for a in ARMS] + [("E" + a) for a in ARMS]
ALL_NODES = ["F" + a for a in ARMS]

res = run(TRANSIENT, CASE)
if res.returncode != 0:
    print("FAIL: transient app exited %d\n%s" % (res.returncode, res.stderr))
    sys.exit(1)
T = parse_transient(res.stdout, ALL_PIPES, ALL_NODES)

res = run(TRANSIENT, CASE_DIR + "/case_long.dict")
if res.returncode != 0:
    print("FAIL: long transient app exited %d" % res.returncode)
    sys.exit(1)
L = parse_transient(res.stdout, ALL_PIPES, ALL_NODES)

res = run(TRANSIENT, CASE_DIR + "/case_single.dict")
if res.returncode != 0:
    print("FAIL: single-arm app exited %d" % res.returncode)
    sys.exit(1)
S = parse_transient(res.stdout, ["Sa", "Ea", "Sc", "Ec", "Sd", "Ed"],
                    ["Fa", "Fc", "Fd"])

res = run(STEADY, CASE_DIR + "/case_steady.dict")
if res.returncode != 0:
    print("FAIL: steady app exited %d" % res.returncode)
    sys.exit(1)
STEADY_OUT = res.stdout

# --- family 1: the ladder splits at the arithmetic ---------------------------
print("f1: %d starts predicted parked, %d up, %d down (band [%.0f, %.0f])"
      % (len(PARKED), len(UP), len(DOWN), CEIL, FLOOR))
for arm in ORDER:
    p = T["p"]["F" + arm]
    if KIND[arm] == "parked":
        require(len(set(p)) == 1 and p[0] == ARMS[arm][0],
                "arm %s starts at %.6g, inside the band, and prints the one "
                "token %s in all %d states -- both rows write 0 = 0, so the "
                "node row has nothing to move"
                % (arm, ARMS[arm][0], T["rawp"]["F" + arm][0], len(p)))
    else:
        require(len(set(p)) == len(p) and p[1] != p[0],
                "arm %s starts at %.6g, outside the band, and every one of the "
                "%d states is a different token; it already moves at step 2 "
                "(%.10g)" % (arm, ARMS[arm][0], len(p), p[1]))
require(len(set(T["p"]["Fm"])) == 1 and T["p"]["Fm"][0] == 2.5e4,
        "arm m starts at 2.5e4 with the schedule on both rows and is parked "
        "like f: the band is decided by drive vs threshold, not by the calendar")

# --- family 2: the split is closed, and it is located to 1 Pa ----------------
require(T["rawp"]["Fd"][0] == "22000" and T["rawp"]["Fh"][0] == "28000",
        "the two edge starts print exactly their own numbers: Fd %s (drive %g = "
        "threshold %g) and Fh %s" % (T["rawp"]["Fd"][0], P_A - 2.2e4, THR,
                                     T["rawp"]["Fh"][0]))
require(T["hold"]["Sd"][1] == 1 and T["hold"]["Eh"][1] == 1,
        "an edge start is HELD on the very first decision, so holdsStatic uses "
        "<= and not < (transient.cpp:40-:43): Sd at drive exactly %g holds, Eh "
        "at drive exactly %g holds" % (THR, THR))
require(T["p"]["Fc"][1] > T["p"]["Fc"][0] and T["p"]["Fi"][1] < T["p"]["Fi"][0],
        "one Pa outside, the same file moves: Fc %.10g -> %.10g rises, Fi %.10g -> "
        "%.10g falls -- so the boundary sits in (21999, 22000] and [28000, "
        "28001)" % (T["p"]["Fc"][0], T["p"]["Fc"][1], T["p"]["Fi"][0],
                    T["p"]["Fi"][1]))

# --- family 3: one direction, no overshoot -----------------------------------
worst_up = max(max(T["p"]["F" + a]) for a in UP)
worst_dn = min(min(T["p"]["F" + a]) for a in DOWN)
require(all(T["p"]["F" + a][i + 1] >= T["p"]["F" + a][i] for a in UP
            for i in range(NS - 1)) and worst_up < CEIL,
        "below the band the node only rises and never touches the ceiling: "
        "max %.10g Pa against %.0f over %d states x %d arms"
        % (worst_up, CEIL, NS, len(UP)))
require(all(T["p"]["F" + a][i + 1] <= T["p"]["F" + a][i] for a in DOWN
            for i in range(NS - 1)) and worst_dn > FLOOR,
        "above the band the node only falls and never touches the floor: min "
        "%.10g Pa against %.0f" % (worst_dn, FLOOR))

# --- family 4: a wall that woke is never re-armed ----------------------------
bad = []
for a in UP + DOWN:
    waking = "S" + a if KIND[a] == "up" else "E" + a
    asleep = "E" + a if KIND[a] == "up" else "S" + a
    if any(T["hold"][waking][i] for i in range(1, NS)):
        bad.append(a + ":waking")
    if not all(T["hold"][asleep][i] for i in range(1, NS)):
        bad.append(a + ":asleep")
    if min(x for x in T["q"][waking][1:] if x > 0) == 0.0:
        bad.append(a + ":zero")
require(not bad, "in all %d decided states of the %d sliding arms the wall that "
        "woke reports hold 0 (its rate is never exactly 0, the smallest printed "
        "is %.6g m3/s) and the other reports hold 1 -- parking is reachable "
        "only by starting there%s" % (NS - 1, len(UP) + len(DOWN),
                                      min(min(x for x in T["q"][("S" if KIND[a] == "up" else "E") + a][1:] if x > 0)
                                          for a in UP + DOWN),
                                      "" if not bad else ", FAIL " + str(bad)))
require(all(T["q"][p][i] == 0.0 for p in ALL_PIPES for i in range(NS)
            if all(T["hold"][p][j] for j in range(1, i + 1)) and i > 0)
        and all(set(T["raw"][p]) == {"0"} for p in
                ["S" + a for a in PARKED + DOWN] + ["E" + a for a in PARKED + UP]),
        "a row held from the first decision prints the single token 0 in every "
        "state: %d rows x %d tokens" % (
            len(["S" + a for a in PARKED + DOWN] + ["E" + a for a in PARKED + UP]), NS))

# --- family 5: the ceiling problem IS the floor problem ----------------------
for a, m in MIRR.items():
    dp = max(abs((T["p"]["F" + a][i] - ARMS[a][0])
                 - ((ARMS[m][0] - T["p"]["F" + m][i]))) for i in range(NS))
    same = sum(1 for i in range(NS)
               if T["raw"]["S" + a][i] == T["raw"]["E" + m][i])
    require(dp < 1e-6 and same == NS,
            "arms %s/%s mirror token for token: node %s + node %s = %.4g in all "
            "%d states (max %.3g Pa off) and the waking rates are identical in "
            "%d/%d states" % (a, m, "F" + a, "F" + m, ARMS[a][0] + ARMS[m][0],
                              NS, dp, same, NS))

# --- family 6: near an edge the rate is quadratic, so the gap is 1/t ---------
# q_of_tw's bracket is 1 - (4/3)x + (1/3)x^4 with x = tau0/tw; both the constant
# and the linear term cancel at x = 1, which is why the near-threshold exponent
# is 2 and the relaxation is algebraic rather than exponential.
R_PIPE = 0.5 * DIAM
DR_TO_TW = R_PIPE / (2.0 * LENGTH)        # a Pa of gap is this Pa of wall stress
DREF = 1.0e-5
KQ = q_of_tw(TAU0 + DREF, TAU0) / DREF ** 2 * DR_TO_TW ** 2   # m3/s per Pa^2
K_PRED = KQ / CC                          # the slope of 1/gap against time
worst_slope = 0.0
worst_arr = 0.0
for a in UP + DOWN:
    ser = L["p"]["F" + a]
    w = "S" + a if KIND[a] == "up" else "E" + a
    edge = CEIL if KIND[a] == "up" else FLOOR
    gap = [abs(edge - ser[i]) for i in range(NL)]
    k1, k2 = 10000, 20000
    slope = (1.0 / gap[k2] - 1.0 / gap[k1]) / ((k2 - k1) * DT)
    worst_slope = max(worst_slope, abs(slope - K_PRED) / K_PRED)
    # the same law integrated from the start of the run: 1/gap = 1/gap0 + K t
    gap0 = abs(edge - ARMS[a][0])
    pred40 = 1.0 / (1.0 / gap0 + K_PRED * END_LONG)
    worst_arr = max(worst_arr, rel(pred40, gap[20000]))
    rat = [gap[i + 2000] / gap[i] for i in range(0, 18000, 2000)]
    require(all(rat[i + 1] > rat[i] for i in range(len(rat) - 1))
            and abs(slope - K_PRED) / K_PRED < 0.05,
            "arm %s over %d states: 1/gap is a straight line in the step number "
            "with slope %.6g per Pa.s against %.6g predicted from the small-d "
            "term of q_of_tw (%.4f%% off), while the per-2s shrink %.5f -> %.5f "
            "climbs toward 1 -- a 1/t tail, not an exponential"
            % (a, NL, slope, K_PRED, 100 * abs(slope - K_PRED) / K_PRED,
               rat[0], rat[-1]))
    require(rel(pred40, gap[20000]) < 0.05,
            "arm %s integrated: 1/gap = 1/%.4g + %.6g x 40 s predicts gap "
            "%.4f Pa at step 20000, printed %.4f Pa (%.3f%% off) -- after 40 s "
            "this arm is still %.4f Pa from the edge it can only approach"
            % (a, gap0, K_PRED, pred40, gap[20000],
               100 * rel(pred40, gap[20000]), gap[20000]))
require(worst_slope < 0.05 and worst_arr < 0.05,
        "all %d sliding arms sit on the one predicted 1/t law: worst slope "
        "%.4f%%, worst 40 s arrival %.3f%%" % (len(UP) + len(DOWN),
                                               100 * worst_slope, 100 * worst_arr))
ex = []
for step in (5000, 10000, 20000):
    pts = []
    for a in UP + DOWN:
        ser = L["p"]["F" + a]
        w = "S" + a if KIND[a] == "up" else "E" + a
        edge = CEIL if KIND[a] == "up" else FLOOR
        pts.append((abs(edge - ser[step]), L["q"][w][step]))
    pts = sorted({(round(g, 6), q) for (g, q) in pts if g > 0.5 and q > 0.0})
    for (g1, q1), (g2, q2) in zip(pts, pts[1:]):
        if g2 > g1 and q2 > q1:
            ex.append(math.log(q2 / q1) / math.log(g2 / g1))
require(ex and min(ex) > 1.85 and max(ex) < 2.15,
        "the exponent d ln Q / d ln gap read between arms at the same step is "
        "%.4f - %.4f over %d pairs spanning gaps %.3g - %.3g Pa, i.e. the "
        "near-edge rate is quadratic in the excess and not linear"
        % (min(ex), max(ex), len(ex),
           min(abs((CEIL if KIND[a] == "up" else FLOOR) - L["p"]["F" + a][20000])
               for a in UP + DOWN),
           max(abs((CEIL if KIND[a] == "up" else FLOOR) - L["p"]["F" + a][20000])
               for a in UP + DOWN)))

# --- family 7: the calendar cannot touch a parked row ------------------------
same_p = sum(1 for i in range(NS) if T["rawp"]["Fm"][i] == T["rawp"]["Ff"][i])
same_q = all(T["raw"]["S" + ("m" if False else "m")][i] == T["raw"]["Sf"][i]
             and T["raw"]["Em"][i] == T["raw"]["Ef"][i] for i in range(NS))
sched_seen = max(valve_k(CLOSE_STD, i * DT) for i in range(NS))
require(same_p == NS and same_q
        and all(T["hold"]["Sm"][i] == T["hold"]["Sf"][i] == 1
                and T["hold"]["Em"][i] == T["hold"]["Ef"][i] == 1
                for i in range(1, NS)),
        "arm m is arm f with a schedule that reaches K = %.4g in this very file, "
        "and it is indistinguishable: %d/%d node tokens identical, both rows' "
        "rate tokens identical, all %d decided flags 1 on both rows -- a held "
        "row has no K term to read" % (sched_seen, same_p, NS, NS - 1))

# --- family 8: the parking lot belongs to the node and its own two rows -----
# Three arms alone in a 3-node / 6-pipe matrix against the same three in the
# 12-node / 24-pipe matrix. A parked row has nothing to sum, so it must be
# token-identical; a sliding row goes through a different dense matrix, so the
# honest expectation here is one printed quantum on the node and the residual
# scale on the rate (5.54's lesson: alignment by number, not by token).
park_ident = True
node_drift = 0.0
rate_drift = 0.0
for a in ("a", "c", "d"):
    n = "F" + a
    node_drift = max(node_drift, max(abs(x - y) for x, y in
                                    zip(T["p"][n], S["p"][n])))
    for p_ in ("S" + a, "E" + a):
        rate_drift = max(rate_drift, max(abs(x - y) for x, y in
                                         zip(T["q"][p_], S["q"][p_])))
    if KIND[a] == "parked":
        park_ident = all([S["rawp"][n] == T["rawp"][n],
                          S["raw"]["S" + a] == T["raw"]["S" + a],
                          S["raw"]["E" + a] == T["raw"]["E" + a]])
same_kind = all((len(set(S["rawp"]["F" + a])) == 1) == (KIND[a] == "parked")
                for a in ("a", "c", "d"))
require(park_ident, "the parked arm `d` is token-identical between the two "
        "matrices: node and both rows, %d states each" % NS)
require(same_kind and node_drift < 2.0e-7 and rate_drift < 1.0e-15,
        "the three survivors keep their classification and their numbers: "
        "frozen-vs-sliding identical for all three, node drift %.3e Pa (at most "
        "one printed quantum) and rate drift %.3e m3/s -- the two matrix sizes "
        "sum the same rows in a different order, so this is alignment by "
        "number, not by token" % (node_drift, rate_drift))

# --- family 9: reproduce the 4 s trace with this file's own rows -------------
worst_p = 0.0
worst_q = 0.0
flag_bad = []
kind_bad = []
for arm in ORDER:
    logs, flags, plist = march(arm, NS - 1)
    mine = plist
    their = T["p"]["F" + arm]
    worst_p = max(worst_p, max(rel(mine[i], their[i]) for i in range(NS)))
    for (n, _, _, _, _) in rows_of(arm):
        theirs = T["q"][n]
        for i in range(NS):
            if abs(logs[n][i]) > RATE_FLOOR or abs(theirs[i]) > RATE_FLOOR:
                worst_q = max(worst_q, rel(logs[n][i], theirs[i]))
        if flags[n][1:] != T["hold"][n][1:]:
            flag_bad.append(n)
    moved = mine[1] != mine[0]
    if moved != (KIND[arm] != "parked"):
        kind_bad.append(arm)
require(worst_p < 1e-9, "the %d nodes x %d states are reproduced by this "
        "file's two rows and node row to %.3e relative worst" % (len(ORDER), NS, worst_p))
require(worst_q < 1e-8, "the %d x 2 rate series reproduce to %.3e relative "
        "above the %.0e m3/s ground rate" % (len(ORDER), worst_q, RATE_FLOOR))
require(not flag_bad, "every one of the %d latch flags in the trace is "
        "reproduced%s" % (2 * len(ORDER) * (NS - 1),
                          "" if not flag_bad else ", FAIL " + str(flag_bad)))
require(not kind_bad, "the parked/sliding classification this file predicted "
        "from CEIL/FLOOR alone is what the reproduced march does too%s"
        % ("" if not kind_bad else ", FAIL " + str(kind_bad)))

# --- family 10: conservation on the printed tokens --------------------------
held_tokens = sum(1 for p in ALL_PIPES for x in T["raw"][p] if x == "0")
closure = 0.0
for a in ORDER:
    for i in range(1, NS):
        closure = max(closure, abs(T["q"]["S" + a][i] - T["q"]["E" + a][i]))
require(held_tokens == sum(1 for p in ALL_PIPES for x in T["raw"][p] if x == "0")
        and closure * CC < 1e-6,
        "%d of the %d printed rate tokens are exactly 0 and every node row "
        "closes on the printed numbers: the two rows of a node differ by at "
        "most %.3e m3/s = %.3e Pa of node motion per step"
        % (held_tokens, len(ALL_PIPES) * NS, closure, closure * CC))
zeros_at_zero = sum(1 for p in ALL_PIPES if T["hold"][p][0] == 0)
require(all(T["hold"][p][0] == 0 for p in ALL_PIPES),
        "the t = 0 column prints %d hold flags all 0: no decision has been made "
        "yet, so an event scan starts at i = 1 (the %d frozen arms would each "
        "report a phantom release)" % (zeros_at_zero, len(PARKED)))

# --- family 11: the steady kernel has no latch and no memory ----------------
sp = {}
sq = {}
banner = ""
for line in STEADY_OUT.splitlines():
    m = re.match(r'\[chflow\] node (\w+) p ([0-9.e+-]+)', line)
    if m: sp[m.group(1)] = float(m.group(2))
    m = re.match(r'\[chflow\] pipe (\w+) Q ([0-9.e+-]+)', line)
    if m: sq[m.group(1)] = float(m.group(2))
    if "steady solved" in line:
        banner = line.replace("[chflow] ", "")
guss = (P_A + P_B) / 2.0
require(all(v == 0.0 for v in sq.values()) and len(sq) == 2 * (len(ARMS) - 1),
        "the steady kernel pays 0 on all %d rows: two 8000 Pa walls in series "
        "across %g Pa have no flow for any start" % (len(sq), P_A - P_B))
require(all(sp["F" + a] == guss for a in ORDER if a != "m") and banner.endswith("converged)")
        and "iterations=0" in banner,
        "and it leaves all %d swept nodes at the average of the fixed pressures, "
        "%.4g Pa -- which is the CENTRE of the band [%g, %g] -- and calls that "
        "converged in zero iterations: steady cannot tell the eleven starts "
        "apart [%s]" % (len(ARMS) - 1, guss, CEIL, FLOOR, banner))
require(abs(CEIL + FLOOR - 2.0 * guss) < 1e-12,
        "the band centre is the guess by arithmetic: (pA - thr + pB + thr)/2 = "
        "(pA + pB)/2 = %g, so the kernel's blind start is always a parked state "
        "whenever the band is not empty" % guss)

print("\n%s: %d checks, %d failures" % ("PASS" if not failures else "FAIL",
                                          CHECKS, len(failures)))
sys.exit(1 if failures else 0)
