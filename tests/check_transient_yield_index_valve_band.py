#!/usr/bin/env python3
"""5.58: what can a calendar touch once the band itself is index-blind?

5.57 split the parking band: the threshold `4 tau_0 L / d` and the `<=` that
latches on it read no rheology, so n = 0.5 / 1 / 2 all park in the same
[22000, 28000] Pa, while everything downstream of the latch reads the index.
That left 5.56's boundary 4 open for HB walls -- can a valve schedule reach into
the dead band?  The schedule lives in the SAME row as the wall law

    (inertia/dt) dQ + 4L/d * tau_w + K(t) rho |Q| Q / (2 A^2) = drive
        (transient.cpp:363-:382)

and the held row is replaced by `Q = 0` before `kv` is ever read (:353-:355), so
the two halves of this file's answer are predicted by those two excerpts and
measured here.  Ten arms per index, one 4 s clock, three windows inside it
(open until 0.5 s, closing 0.5-1.0 s, closed 1.0-2.0 s, reopening 2.0-2.5 s,
open after):

  1. the desk and the latch never see the calendar: a parked arm with the full
     schedule on BOTH rows is token-for-token the parked arm without it, at every
     index, and the edge arm whose drive is exactly 8000 is still HELD on the
     very first decision with K(t) sweeping to 1e8 behind it;
  2. the schedule's first byte is the schedule's own arithmetic: an armed row is
     identical to its unarmed twin until the first sample strictly past
     closeStart, and the armed row's wall stress -- recovered with the K term
     moved off the drive -- is this file's HB section integral from the FIRST
     decided step;
  3. what the closed window costs, per index: the throttle on the rate, the Pa
     of approach the arm loses and still has lost a full second after the valve
     is open again, and the payback -- after reopening EVERY state puts the
     armed arm's rate above its unarmed twin's;
  4. the calendar re-prices the drive, not the wall: the valve's own
     rho K |Q| Q / 2A^2 peaks at a fraction of the 8000 Pa threshold that is
     monotone in n, while the wall stress the armed row implies stays within
     ~11% of its twin's -- and on a 1 Pa start the same schedule prints nothing;
  5. the calendar respects the mirror: `|Q|Q` is odd, so the same schedule on the
     waking row of an UP arm and of its DOWN mirror still prints mirrored
     tokens;
  6. the band is still not crossable and not enterable by a calendar, at any
     index, and the throttled rows never re-arm because their rate is never
     exactly 0;
  7. the payback: with K back to exactly 0 the armed arm's debt to its twin
     shrinks strictly and monotonically and is repaid by an index-dependent
     fraction -- at m = 3 the bill is nearly permanent;
  8. every node row still closes on the printed tokens alone, scheduled rows
     included: the calendar moves money inside one row, not across a node.
"""
import math
import os
import re
import subprocess
import sys

TRANSIENT = sys.argv[1]
CASE = sys.argv[2]
CASE_DIR = CASE.rsplit("/", 1)[0]

RHO, MU, L, D = 850.0, 0.1, 20.0, 0.02
R_PIPE, AREA = D / 2.0, math.pi * D * D / 4.0
CC = 0.1 / 1.0e8
WALL_SCALE = 4.0 * L / D
THR = WALL_SCALE * 2.0
DT = 0.002
INERTIA = RHO * L / AREA
I_DT = INERTIA / DT
P_A, P_B, TAU0 = 3.0e4, 2.0e4, 2.0
SLOPE = D / (4.0 * L)
VK = RHO / (2.0 * AREA * AREA)
CS, CE, KC, OS_, OE = 0.5, 1.0, 1.0e8, 2.0, 2.5
FIRST = int(round(CS / DT)) + 1          # first sample strictly past closeStart
CEILING, FLOOR = P_A - THR, P_B + THR     # 22000 / 28000

CASES = [("5", 0.5, 1.0), ("1", 1.0, 0.1), ("2", 2.0, 0.02)]
PARKED = "ahbc"
SCHED = {"h": ("S", "E"), "b": ("S", "E"), "e": ("S",), "i": ("S",), "f": ("E",)}
UP = "degi"          # below the ceiling, waking row is the S row
DOWN = "jf"          # above the floor, waking row is the E row
ARMS = PARKED + UP + DOWN

failures = []
CHECKS = 0


def require(cond, what):
    global CHECKS
    CHECKS += 1
    print("ok:  " if cond else "FAIL:", what)
    if not cond:
        failures.append(what)


def run(case):
    out = subprocess.run([TRANSIENT, case], cwd=CASE_DIR, capture_output=True,
                         text=True)
    require(out.returncode == 0, "%s exited %d: %s"
            % (case, out.returncode, out.stderr.strip()[:300]))
    N, P, F = {}, {}, {}
    for line in out.stdout.splitlines():
        m = re.match(r"\[chflow\] node (\S+) p (\S+) ", line)
        if m:
            N.setdefault(m.group(1), []).append(m.group(2))
            continue
        m = re.match(r"\[chflow\] pipe (\S+) Q (\S+)", line)
        if m:
            P.setdefault(m.group(1), []).append(m.group(2))
            continue
        m = re.match(r"\[chflow\] pipe (\S+) static hold (\d)", line)
        if m:
            F.setdefault(m.group(1), []).append(int(m.group(2)))
    return N, P, F


def valve_k(t):
    """transient.cpp:125-:133 verbatim: linear ramps, closedK held between them,
    and the open value is 0 here, so an armed row is unarmed before cs."""
    if t < CS:
        return 0.0
    if t < CE:
        return (t - CS) / (CE - CS) * KC
    if t < OS_:
        return KC
    if t < OE:
        return KC + (t - OS_) / (OE - OS_) * (0.0 - KC)
    return 0.0


def q_hb(tw, n, K):
    if tw <= TAU0:
        return 0.0
    d = tw - TAU0
    p = 1.0 / n
    return math.pi * R_PIPE ** 3 * tw ** -3 * K ** -p * (
        d ** (p + 3.0) / (p + 3.0) + 2.0 * TAU0 * d ** (p + 2.0) / (p + 2.0)
        + TAU0 * TAU0 * d ** (p + 1.0) / (p + 1.0))


def waking(up):
    return "S" if up else "E"


def drive_of(up, pn, i):
    return (P_A - pn[i]) if up else (pn[i] - P_B)


def tw_row(up, pn, q, i, n, K, with_valve=True):
    """Wall stress the row itself implies (transient.cpp:363-:382), with the
    valve's rho K|Q|Q/2A^2 either moved off the drive or ignored."""
    drive = drive_of(up, pn, i)
    v = float(q[i])
    acc = I_DT * (v - float(q[i - 1]))
    valve = valve_k(i * DT) * VK * v * abs(v) if with_valve else 0.0
    return (drive - acc - valve) / WALL_SCALE


R = {}
for tag, n, K in CASES:
    R[tag] = run("case_%s.dict" % tag)

# --- family 1: the calendar cannot reach into the desk or the latch -----------
for tag, n, K in CASES:
    N, P, F = R[tag]
    ns = len(N["F" + tag + "a"])
    frozen = {a: set(N["F" + tag + a]) for a in PARKED}
    require(frozen == {"a": {"25000"}, "h": {"25000"}, "b": {"22000"},
                       "c": {"28000"}},
            "n=%s: the four parked starts print one token each over %d states "
            "(25000 / 25000 / 22000 / 28000) -- the threshold is 4 tau_0 L/d "
            "with no n, no K and no valve in it (transient.cpp:28-:31)" % (n, ns))
    same = (N["F" + tag + "a"] == N["F" + tag + "h"]
            and P["S" + tag + "a"] == P["S" + tag + "h"]
            and P["E" + tag + "a"] == P["E" + tag + "h"])
    require(same,
            "n=%s: arm h IS arm a token for token -- %d node tokens and %d+%d "
            "rate tokens identical -- even though both of h's rows carry the "
            "close-to-%.0e / reopen schedule, because a held row becomes Q = 0 "
            "before kv is read (transient.cpp:353-:355)"
            % (n, ns, ns, ns, KC))
    require(set(F["S" + tag + "h"][1:]) == {1}
            and set(F["E" + tag + "h"][1:]) == {1}
            and all(v == "0" for v in P["S" + tag + "h"]) 
            and all(v == "0" for v in P["E" + tag + "h"]),
            "n=%s: both of h's rows are HELD in every one of the %d decided "
            "states and print a single token 0 through the closed window, so the "
            "calendar has no term to act on" % (n, ns - 1))
    for a, row in (("b", "S"), ("c", "E")):
        require(F[row + tag + a][1] == 1 and set(P[row + tag + a]) == {"0"},
                "n=%s: the %s edge arm %s is still HELD on the very first "
                "decision (drive exactly %.0f = threshold %.0f) and never wakes "
                "while K(t) sweeps to %.0e -- the `<=` is calendar-blind too "
                "(transient.cpp:40-:43)" % (n, "lower" if a == "b" else "upper",
                                            a, THR, THR, KC))
    require(len(F) == 20 and all(v[0] == 0 for v in F.values()),
            "n=%s: t = 0 prints %d hold flags, every one 0, so an event scan "
            "must still start at i = 1" % (n, len(F)))

# --- family 2: the armed row's wall stress, valve term and all ---------------
worst_w = {}
worst_wo = {}
first_diff = {}
for tag, n, K in CASES:
    N, P, F = R[tag]
    ns = len(N["F" + tag + "a"])
    for a in UP + DOWN:
        up = a in UP
        row = waking(up)
        pn = [float(x) for x in N["F" + tag + a]]
        q = P[row + tag + a]
        armed = a in SCHED
        w = wo = 0.0
        for i in range(1, ns):
            v = float(q[i])
            if v == 0.0:
                continue
            tw = tw_row(up, pn, q, i, n, K, armed)
            if tw > TAU0:
                w = max(w, abs(q_hb(tw, n, K) - abs(v)) / abs(v))
            tw2 = tw_row(up, pn, q, i, n, K, False)
            if tw2 > TAU0:
                wo = max(wo, abs(q_hb(tw2, n, K) - abs(v)) / abs(v))
        worst_w[(tag, a)] = w
        worst_wo[(tag, a)] = wo
        require(w < 2e-6,
                "n=%s arm %s: %s row's printed rate IS this file's HB section "
                "integral at the wall stress the row implies -- drive minus "
                "(inertia/dt) dQ%s, divided by 4L/d -- worst relative %.3e over "
                "%d decided states, from the FIRST decided step"
                % (n, a, row + tag + a,
                   " minus K(t) rho |Q|Q/2A^2" if armed else " (no valve term: "
                   "this row declares none)", w, ns - 1))
    for a, twin in (("e", "d"), ("i", "g"), ("f", "j")):
        diffs = [i for i in range(ns)
                 if N["F" + tag + a][i] != N["F" + tag + twin][i]]
        first_diff[(tag, a)] = diffs[0] if diffs else None
        require(all(N["F" + tag + a][i] == N["F" + tag + twin][i]
                    for i in range(min(FIRST, ns))),
                "n=%s: arms %s/%s are token for token through the whole open "
                "window (steps 0-%d, K(t) = 0 there) -- the open value is "
                "`valve 0`, so the schedule's bytes enter the row only when the "
                "ramp starts" % (n, a, twin, FIRST - 1))
    armed = [first_diff[(tag, a)] for a in ("e", "i", "f")]
    require(all(x is None or x >= FIRST for x in armed)
            and armed[1] is not None and armed[1] == FIRST,
            "n=%s: the three armed sliders first differ from their unarmed "
            "twins at steps %s -- never before step %d, the first sample "
            "strictly past closeStart %.3f s (K there %.3e), and the far slider "
            "does move off its twin on exactly that step"
            % (n, [str(x) for x in armed], FIRST, CS, valve_k(FIRST * DT)))
    require(worst_w[(tag, "i")] < 2e-6 and worst_wo[(tag, "i")] > 1e-1,
            "n=%s: on the armed far slider the SAME reconstruction read without "
            "the valve term is off by %.3e relative while with it it is off by "
            "%.3e -- the K(t) |Q|Q term lives in the same row as the wall law, "
            "and a checker that forgets it would report a wrong rheology "
            "(transient.cpp:363-:382)" % (n, worst_wo[(tag, "i")],
                                          worst_w[(tag, "i")]))
    require(abs(worst_w[(tag, "g")] - worst_wo[(tag, "g")]) < 1e-12,
            "n=%s: and on its unarmed twin g the two readings are the same "
            "number (%.3e), because that row declares no valve at all"
            % (n, worst_w[(tag, "g")]))
for tag, n, K in CASES:
    late = {a: first_diff[(tag, a)] for a in ("e", "f", "i")}
    if n == 0.5:
        require(late["e"] is None and late["f"] is None,
                "n=0.5: the 1 Pa armed arms never differ from their unarmed "
                "twins in ANY of the 2001 printed states -- at m=3 the schedule "
                "is real but sub-quantum, exactly the blind spot 5.57 found for "
                "the tail it was measured against (far slider first differs at "
                "step %s)" % late["i"])
    else:
        require(late["e"] is not None and late["e"] > FIRST
                and late["i"] is not None and late["i"] == FIRST,
                "n=%s: the far slider separates from its twin at step %d, the "
                "first sample past closeStart (K %.3e), but the 1 Pa armed arm "
                "takes until step %d (t = %.3f s, %.0f steps into the closed "
                "window) -- one calendar, two different printout latencies, set "
                "by the exponent" % (n, late["i"], valve_k(late["i"] * DT),
                                     late["e"], late["e"] * DT,
                                     late["e"] - FIRST))

# --- family 3: what the closed window costs the trajectory, per index --------
T_CLOSED = int(round(1.5 / DT))          # K = closedK held, ramp long finished
T_OPEN = int(round(3.5 / DT))            # K = 0 again, a full second after OE
ARMD = int(round(OE / DT)) + 1
cost = {}
throttle = {}
rate_closed = {}
for tag, n, K in CASES:
    N, P, F = R[tag]
    ns = len(N["F" + tag + "a"])
    pn = [float(x) for x in N["F" + tag + "i"]]
    pt = [float(x) for x in N["F" + tag + "g"]]
    qa = [float(x) for x in P["S" + tag + "i"]]
    qt = [float(x) for x in P["S" + tag + "g"]]
    lost = pt[-1] - pn[-1]                       # armed still below its twin
    moved = pt[-1] - pt[0]                       # the twin's whole approach
    throttle[tag] = (qt[T_CLOSED] / qa[T_CLOSED], qt[T_OPEN] / qa[T_OPEN])
    cost[tag] = (lost, moved, 100.0 * lost / moved)
    rate_closed[tag] = (qa[T_CLOSED], qt[T_CLOSED])
    require(qa[T_CLOSED] < qt[T_CLOSED] and qa[T_OPEN] > qt[T_OPEN],
            "n=%s: the far slider's rate at t = %.1f s (K = %.0e held) is %.6e "
            "against its unarmed twin's %.6e -- throttled by %.2fx -- and at "
            "t = %.1f s, a full second after the valve is wide open again, it is "
            "%.6e against %.6e: the armed arm now flows %.2fx its twin"
            % (n, T_CLOSED * DT, KC, qa[T_CLOSED], qt[T_CLOSED],
               qt[T_CLOSED] / qa[T_CLOSED], T_OPEN * DT, qa[T_OPEN], qt[T_OPEN],
               qa[T_OPEN] / qt[T_OPEN]))
    require(lost > 0.0 and moved > lost,
            "n=%s: the closed valve costs this arm %.4f Pa of the %.4f Pa of "
            "approach it would otherwise have made (%.1f%% of the way), and the "
            "%.4f Pa is still missing at t = %.3f s with K back to exactly 0 -- "
            "a calendar moves a trajectory, it does not move a law"
            % (n, lost, moved, 100.0 * lost / moved, lost, (ns - 1) * DT))
    require(all(qa[i] > qt[i] for i in range(ARMD, ns)),
            "n=%s: after the open window every one of the %d states puts the "
            "armed arm's rate ABOVE its unarmed twin's -- the arm was left %.4f "
            "Pa further from the ceiling, so the same wall law is handed more "
            "drive, and it pays the schedule back at %.2f..%.2f x until the end "
            "of the clock"
            % (n, ns - ARMD, cost[tag][0],
               min(qa[i] / qt[i] for i in range(ARMD, ns)),
               max(qa[i] / qt[i] for i in range(ARMD, ns))))
ordr = [(throttle[tag][0], throttle[tag][1]) for tag, _, _ in CASES]
require(all(x > 1.0 for x, _ in ordr)
        and [x for x, _ in ordr] == sorted(x for x, _ in ordr)
        and all(y < 1.0 for _, y in ordr),
        "the throttle at t = 1.5 s is %.4fx / %.4fx / %.4fx and the rate ratio "
        "at t = 3.5 s is %.4fx / %.4fx / %.4fx (n = 0.5 / 1 / 2, twin over "
        "armed in both cases): the index decides how much a closed valve takes "
        "away, and every index pays it back after the valve is open again -- one "
        "schedule, three different prices"
        % tuple([x for pair in zip(*ordr) for x in pair]))
frac = {tag: cost[tag][2] for tag, _, _ in CASES}
require(min(frac.values()) > 10.0 and max(frac.values()) < 25.0,
        "the same schedule costs the three indices %.1f%% / %.1f%% / %.1f%% of "
        "their own approach even though the approaches themselves are %.1f / "
        "%.1f / %.1f Pa -- the fraction is the calendar's, the distance is the "
        "index's" % (frac["5"], frac["1"], frac["2"],
                     *[cost[tag][1] for tag, _, _ in CASES]))

# --- family 4: the calendar re-prices the drive, not the wall ----------------
vmax = {}
twmax = {}
for tag, n, K in CASES:
    N, P, F = R[tag]
    ns = len(N["F" + tag + "a"])
    pn = [float(x) for x in N["F" + tag + "i"]]
    pt = [float(x) for x in N["F" + tag + "g"]]
    qa = [float(x) for x in P["S" + tag + "i"]]
    qt = [float(x) for x in P["S" + tag + "g"]]
    vpk = twm = 0.0
    pk_at = (0, 1.0)
    for i in range(1, ns):
        k = valve_k(i * DT)
        if k == 0.0:
            continue
        vpa = k * VK * qa[i] * abs(qa[i])
        twa = (P_A - pn[i] - I_DT * (qa[i] - qa[i - 1]) - vpa) / WALL_SCALE
        twt = (P_A - pt[i] - I_DT * (qt[i] - qt[i - 1])) / WALL_SCALE
        if vpa > vpk:
            vpk = vpa
            pk_at = (i, P_A - pn[i])
        twm = max(twm, abs(twa - twt) / twt)
    vmax[tag] = (vpk, vpk / THR, 100.0 * vpk / pk_at[1], pk_at[0])
    twmax[tag] = twm
    require(twm < 0.12,
            "n=%s: with the K(t) rho|Q|Q/2A^2 term moved off the drive, the "
            "armed row's own wall stress never differs from its unarmed twin's "
            "by more than %.4f%% over the %d states where K > 0 (at t = %.1f s: "
            "%.4f Pa against %.4f Pa) -- the schedule did not change which law "
            "the row obeys, it only changed how much of the drive reaches it"
            % (n, 100.0 * twm, sum(1 for i in range(1, ns)
                                   if valve_k(i * DT) > 0.0),
               T_CLOSED * DT,
               (P_A - pn[T_CLOSED] - I_DT * (qa[T_CLOSED] - qa[T_CLOSED - 1])
                - valve_k(T_CLOSED * DT) * VK * qa[T_CLOSED] * qa[T_CLOSED])
               / WALL_SCALE,
               (P_A - pt[T_CLOSED] - I_DT * (qt[T_CLOSED] - qt[T_CLOSED - 1]))
               / WALL_SCALE))
    require(vpk > 300.0,
            "n=%s: the valve's own term peaks at %.1f Pa of drive (step %d, "
            "t = %.3f s) -- %.1f%% of the %.0f Pa threshold the latch reads and "
            "%.1f%% of the drive on the armed row right there -- and none of "
            "that money is allowed to buy a step of the band "
            "(transient.cpp:28-:31)"
            % (n, vpk, vmax[tag][3], vmax[tag][3] * DT, 100.0 * vpk / THR,
               THR, vmax[tag][2]))
share = [vmax[tag][1] for tag, _, _ in CASES]
require(share == sorted(share) and share[0] < 0.1 and share[-1] > 0.4,
        "the peak valve term is %.1f%% / %.1f%% / %.1f%% of the 8000 Pa "
        "threshold at n = 0.5 / 1 / 2, monotone in the index (the term is "
        "K rho |Q| Q / 2A^2 and Q itself is index-dependent): at n = 2 the "
        "calendar spends nearly half the threshold on flow and still cannot "
        "touch the number the latch compares against"
        % tuple(100.0 * s for s in share))
for tag, n, K in CASES:
    N, P, F = R[tag]
    pn = [float(x) for x in N["F" + tag + "e"]]
    pt = [float(x) for x in N["F" + tag + "d"]]
    qa = [float(x) for x in P["S" + tag + "e"]]
    qt = [float(x) for x in P["S" + tag + "d"]]
    r = [qt[T_CLOSED] / qa[T_CLOSED], qt[T_OPEN] / qa[T_OPEN]]
    sep = first_diff[(tag, "e")]
    require(max(abs(r[0] - 1.0), abs(r[1] - 1.0)) < 0.01,
            "n=%s: on the 1 Pa armed arm the same schedule is worth %.4f%% of "
            "the rate at t = 1.5 s and %.4f%% at t = 3.5 s, where the far arm on "
            "the same clock is throttled %.2fx and pays back %.2fx -- the valve "
            "taxes rho K |Q|Q/2A^2, and a near-threshold wall at m = %.2f hardly "
            "moves at all. Its node token leaves its unarmed twin's %s"
            % (n, 100.0 * abs(r[0] - 1.0), 100.0 * abs(r[1] - 1.0),
               throttle[tag][0], throttle[tag][1], (n + 1.0) / n,
               "nowhere in the %d printed states" % len(pn)
               if sep is None else "at step %d" % sep))

# --- family 5: the calendar respects the mirror ------------------------------
for tag, n, K in CASES:
    N, P, F = R[tag]
    ns = len(N["F" + tag + "a"])
    off = max(abs(float(x) + float(y) - (P_A + P_B))
              for x, y in zip(N["F" + tag + "e"], N["F" + tag + "f"]))
    same = sum(1 for x, y in zip(P["S" + tag + "e"], P["E" + tag + "f"])
               if x == y)
    require(off <= 1e-7 and same == ns,
            "n=%s: arms e (21999, valve on its S row) and f (28001, valve on its "
            "E row) mirror token for token -- node sums to 50000 within %.1e and "
            "the waking rates are identical in %d/%d states -- because the loss "
            "term is rho K |Q| Q / 2A^2, odd in Q, so the same schedule closes "
            "the same amount on both sides" % (n, off, same, ns))

# --- family 6: the band is still not enterable and not crossable -------------
for tag, n, K in CASES:
    N, P, F = R[tag]
    ns = len(N["F" + tag + "a"])
    for a in UP:
        pn = [float(x) for x in N["F" + tag + a]]
        require(max(pn) < CEILING,
                "n=%s: arm %s (unscheduled or scheduled, 1 Pa or 2000 Pa outside) "
                "never reaches the ceiling 22000 -- max %.5f Pa, still %.5f Pa "
                "short" % (n, a, max(pn), CEILING - max(pn)))
    for a in DOWN:
        pn = [float(x) for x in N["F" + tag + a]]
        require(min(pn) > FLOOR,
                "n=%s: arm %s never reaches the floor 28000 -- min %.5f Pa"
                % (n, a, min(pn)))
    qa = [float(x) for x in P["S" + tag + "i"]]
    qf = [float(x) for x in P["E" + tag + "f"]]
    window = qa[int(CS / DT):int(OE / DT) + 1] + qf[int(CS / DT):int(OE / DT) + 1]
    require(all(F["S" + tag + a][1:] and set(F["S" + tag + a][1:]) == {0}
                for a in ("e", "i"))
            and all(F["E" + tag + f][1:] and set(F["E" + tag + f][1:]) == {0}
                    for f in DOWN),
            "n=%s: the armed sliders' waking rows are never asked to hold again, "
            "not in the %d decided states of any window; the smallest rate either "
            "of them printed while the valve was doing its work is %.4e m3/s "
            "(%.4e at K = %.0e held) -- never exactly 0, so the latch has no "
            "handle: a calendar throttles a wall, it does not park it"
            % (n, ns - 1, min(window), rate_closed[tag][0], KC))

# --- family 7: the payback curve -- who repays the schedule, and how fast ----
DEBT = [2.5, 2.7, 3.0, 3.5, (len(R["5"][0]["F5a"]) - 1) * DT]
debt = {}
for tag, n, K in CASES:
    N, P, F = R[tag]
    pn = [float(x) for x in N["F" + tag + "i"]]
    pt = [float(x) for x in N["F" + tag + "g"]]
    samples = [pt[int(round(t / DT))] - pn[int(round(t / DT))] for t in DEBT]
    debt[tag] = samples
    require(all(b > 0.0 for b in samples)
            and samples == sorted(samples, reverse=True),
            "n=%s: the armed arm's debt to its unarmed twin is %.4f / %.4f / "
            "%.4f / %.4f / %.4f Pa at t = %.1f / %.1f / %.1f / %.1f / %.3f s -- "
            "strictly shrinking once K is back to 0 and never repaid in full "
            "inside this clock, and it is a debt in the trajectory only: K = 0 "
            "in every one of these states" % (n, *samples, *DEBT))
    repaid = 100.0 * (1.0 - samples[-1] / samples[0])
    require(0.0 < repaid < 100.0,
            "n=%s: %.1f%% of the debt is repaid in the %.1f s of open valve left "
            "on the clock (the remaining %.4f Pa is the part the schedule took "
            "out of the 1/t tail, and a 1/t tail does not catch up)"
            % (n, repaid, DEBT[-1] - DEBT[0], samples[-1]))
rep = {tag: 100.0 * (1.0 - debt[tag][-1] / debt[tag][0]) for tag, _, _ in CASES}
require(rep["5"] < rep["1"] < rep["2"],
        "the same open valve repays %.1f%% / %.1f%% / %.1f%% of the debt at "
        "n = 0.5 / 1 / 2 -- at m = 3 the arm is still %.4f Pa behind after 1.5 s "
        "of free flow while at m = 1.5 it has nearly doubled its catch-up rate. "
        "The calendar's bill is settled by the index" %
        (rep["5"], rep["1"], rep["2"], debt["5"][-1]))

# --- family 8: the armed rows still obey the node row, tokens only -----------
for tag, n, K in CASES:
    N, P, F = R[tag]
    ns = len(N["F" + tag + "a"])
    worst = 0.0
    who = None
    for a in ARMS:
        pn = [float(x) for x in N["F" + tag + a]]
        qs = [float(x) for x in P["S" + tag + a]]
        qe = [float(x) for x in P["E" + tag + a]]
        for i in range(1, ns):
            r = abs((qs[i] - qe[i]) - CC * (pn[i] - pn[i - 1]) / DT)
            if r > worst:
                worst, who = r, (a, i)
    require(worst < 1e-12,
            "n=%s: %d free nodes x %d states -- every node row closes on the "
            "PRINTED tokens alone to %.4e m3/s of volume imbalance (worst on arm "
            "%s step %d, %.2e Pa of node motion per step), scheduled rows "
            "included: the valve moved money inside one row, it did not create "
            "or destroy flux at the node (transient.cpp:395-:416)"
            % (n, len(ARMS), ns, worst, who[0], who[1], worst * DT / CC))

print("\n%s: %d checks, %d failures"
      % ("PASS" if not failures else "FAIL", CHECKS, len(failures)))
for f in failures:
    print("  -", f)
sys.exit(1 if failures else 0)
