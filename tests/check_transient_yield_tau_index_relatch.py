#!/usr/bin/env python3
"""5.60: the tau_0 x flowIndex cross on the two-way latch.

5.59 swept the index at ONE yield stress and found the re-latch present at
n = 1 / 2, absent at n = 0.5 -- but it could not say whether the absence belonged
to the index or to that one threshold. 5.57 swept the index on the parking band
and left "tau_0 x index" explicitly untested. This file sweeps both axes on the
same desk: four arms per tier, strong wall tau_0 = 1.75 / 2.0 / 2.25 / 2.5 against
a tau_0 = 1.5 neighbour, three tiers (n, K) = (0.5, 1.0) / (1.0, 0.1) / (2.0, 0.02).

    threshold = 4 tau_0 L / d = 7000 / 8000 / 9000 / 10000 Pa   (index-blind)
    release   = pA - threshold = 23000 / 22000 / 21000 / 20000 Pa

No kernel change: transient.cpp:28-:31 (threshold), :40-:43 (entry desk),
:353-:355 (in-step clamp, before kv is read), :395-:416 (degree-generic node row).

  1. the desk is blind to BOTH axes: all twelve weak walls release at step 180 in
     all three tiers, and the four arms of a tier are the same strings until the
     first strong wall wakes -- a held neighbour contributes nothing;
  2. tau_0 enters only as 4 tau_0 L / d, and every wake is the first entry drive
     above that number;
  3. the re-latch is a WINDOW in tau_0, not a half-line: 4 of 12 cells sleep,
     2 wake and are never asked, 1 wakes only after the seal broke, 5 never wake;
  4. the two axes buy different things -- tau_0 moves the release pressure by
     3000 Pa and the refill by < 60 Pa, n moves the refill by > 2000 Pa at fixed
     release pressure;
  5. a positive margin does NOT buy a sleep: at n = 0.5 arm c tops out 397 Pa
     ABOVE its release pressure and still never sleeps, because it woke 129 steps
     after the seal broke;
  6. the coast between "drive under threshold" and "flag flips" is the law's and
     the threshold's: lag 0 / 0 at n = 1, 5 / 7 at n = 2;
  7. every unheld step with a root is the HB section integral of its own implied
     wall stress, and the ONLY sub-yield unheld step in a sleeping cell is the one
     immediately before the clamp;
  8. arm b at n = 1, K = mu is 5.59's arm a in a new file, the node rows close on
     printed tokens alone, and the latch stays event-quiet.
"""
import math
import os
import re
import subprocess
import sys

TRANSIENT = sys.argv[1]
CASE = sys.argv[2]
CASE_DIR = CASE.rsplit("/", 1)[0]
SIB = os.path.join(os.path.dirname(os.path.abspath(CASE_DIR)),
                   "transientYieldIndexValveRelatch")

RHO, L, D = 850.0, 20.0, 0.02
R_PIPE, AREA = D / 2.0, math.pi * D * D / 4.0
VOL, BULK = 0.1, 1.0e8
CC = VOL / BULK
WALL_SCALE = 4.0 * L / D
DT = 0.002
I_DT = RHO * L / AREA / DT
P_A, P_B = 3.0e4, 2.0e4
CE, OE = 1.0, 2.7
CLOSED = int(round(CE / DT))                 # 500: first sample with closedK
OPENED = int(round(OE / DT))                 # 1350: first sample with K back to 0
CASES = [("5", 0.5, 1.0), ("1", 1.0, 0.1), ("2", 2.0, 0.02)]
ARMS = [("a", 1.75), ("b", 2.0), ("c", 2.25), ("d", 2.5)]
TAU_W = 1.5
TH = {a: WALL_SCALE * ts for a, ts in ARMS}
REL = {a: P_A - TH[a] for a, _ in ARMS}
FREE = tuple("F" + a for a, _ in ARMS)

failures = []
CHECKS = 0


def require(cond, what):
    global CHECKS
    CHECKS += 1
    print("ok:  " if cond else "FAIL:", what)
    if not cond:
        failures.append(what)


def run(case_dir, name):
    out = subprocess.run([TRANSIENT, name], cwd=case_dir, capture_output=True,
                         text=True)
    require(out.returncode == 0, "%s/%s exited %d: %s"
            % (os.path.basename(case_dir), name, out.returncode,
               out.stderr.strip()[:300]))
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


def q_hb(tw, tau0, n, K):
    """HB volumetric flow at wall stress tw -- the bracket 5.47/5.57/5.59 used."""
    if tw <= tau0:
        return 0.0
    d, p = tw - tau0, 1.0 / n
    return math.pi * R_PIPE ** 3 * tw ** -3 * K ** -p * (
        d ** (p + 3.0) / (p + 3.0) + 2.0 * tau0 * d ** (p + 2.0) / (p + 2.0)
        + tau0 * tau0 * d ** (p + 1.0) / (p + 1.0))


def events(flags):
    return [i for i in range(1, len(flags)) if flags[i] != flags[i - 1]]


R = {}
for tag, n, K in CASES:
    R[tag] = run(CASE_DIR, "case_%s.dict" % tag)
NS = len(R["1"][0]["Fa"])
FL = {(t, a): [float(x) for x in R[t][0]["F" + a]] for t, _, _ in CASES for a, _ in ARMS}
QL = {(t, w): [float(x) for x in R[t][1][w]] for t, _, _ in CASES for w in
      [p + a for a, _ in ARMS for p in ("S", "W", "D")]}
EV = {(t, a): events(R[t][2]["S" + a]) for t, _, _ in CASES for a, _ in ARMS}


def entry(t, a, i):
    """The drive step i's latch read: the node as step i-1 left it."""
    return P_A - FL[(t, a)][i - 1]


def classify(t, a):
    ev = EV[(t, a)]
    if len(ev) == 1:
        return "NEVER"
    if ev[1] > OPENED:
        return "LATE"
    return "SLEPT" if len(ev) > 2 else "AWAKE"


KIND = {(t, a): classify(t, a) for t, _, _ in CASES for a, _ in ARMS}
BEST = {(t, a): max(FL[(t, a)][CLOSED:OPENED + 1]) for t, _, _ in CASES for a, _ in ARMS}

# --- family 1: the desk is blind to both axes, and a held neighbour is silent ---
for tag, n, K in CASES:
    N, P, F = R[tag]
    weak = [events(F["W" + a]) for a, _ in ARMS]
    require(all(w == [1, 180] for w in weak),
            "n=%s: all FOUR weak walls (tau_0 = 1.5, threshold 6000 Pa) release at "
            "step 180 -- %s -- whatever their strong neighbour's tau_0 is. Before "
            "any wall flows the node is the drain's arithmetic alone, so the desk "
            "(:40-:43) reads neither n nor the neighbour's yield stress"
            % (n, weak))
    wakes = [EV[(tag, a)][1] if len(EV[(tag, a)]) > 1 else NS for a, _ in ARMS]
    first = min(wakes)
    same_q = all(P["W" + a][:first] == P["W" + ARMS[0][0]][:first] for a, _ in ARMS)
    same_p = all(N["F" + a][:first] == N["F" + ARMS[0][0]][:first] for a, _ in ARMS)
    require(same_q and same_p,
            "n=%s: the four arms are the SAME STRINGS for the first %d states -- "
            "weak-wall rates and node pressures alike -- because every strong wall "
            "is still held there and a held row is x=0, r=0, J=1 (:353-:355). The "
            "prefix ends exactly at the earliest strong wake (%s)"
            % (n, first, wakes))
    div = [i for i in range(NS) if len({N["F" + a][i] for a, _ in ARMS}) > 1]
    require(not div or div[0] >= first,
            "n=%s: the node tokens first split at %s, never before the earliest "
            "strong wake %d -- the four arms differ only once a strong row stops "
            "being zero" % (n, div[0] if div else "never", first))
    held = [a for a, _ in ARMS if KIND[(tag, a)] == "NEVER"]
    for a in held:
        other = "d" if a == "c" else "c"
        upto = EV[(tag, other)][1] if len(EV[(tag, other)]) > 1 else NS
        require(N["F" + a][:upto] == N["F" + other][:upto]
                and P["W" + a][:upto] == P["W" + other][:upto],
                "n=%s: arm %s (tau_0 = %.2f, never wakes) and arm %s are the same "
                "strings for %d states -- two different yield stresses that the "
                "desk never asks about are indistinguishable in print"
                % (n, a, dict(ARMS)[a], other, upto))

# --- family 2: tau_0 is exactly 4 tau_0 L/d, and every wake is the first crossing
for tag, n, K in CASES:
    for a, ts in ARMS:
        ev = EV[(tag, a)]
        require(TH[a] == 4000.0 * ts and REL[a] == P_A - TH[a],
                "n=%s arm %s: threshold %.0f Pa = 4 tau_0 L/d with tau_0 = %.2f, "
                "release pressure %.0f Pa -- pure geometry, no n anywhere in it "
                "(transient.cpp:28-:31)" % (n, a, TH[a], ts, REL[a]))
        if len(ev) == 1:
            closest = max(entry(tag, a, i) for i in range(1, NS))
            require(closest <= TH[a],
                    "n=%s arm %s: NEVER woken -- the largest entry drive in %d "
                    "states is %.4f Pa, %.4f Pa short of %.0f, so :40-:43 holds it "
                    "from the first decision to the last"
                    % (n, a, NS - 1, closest, TH[a] - closest, TH[a]))
            continue
        i = ev[1]
        require(entry(tag, a, i) > TH[a] >= entry(tag, a, i - 1),
                "n=%s arm %s: wakes at step %d on the FIRST entry drive above "
                "%.0f Pa (%.4f there, %.4f one step earlier, overshoot %.4f Pa = "
                "%.3f%% of the threshold)"
                % (n, a, i, TH[a], entry(tag, a, i), entry(tag, a, i - 1),
                   entry(tag, a, i) - TH[a],
                   100.0 * (entry(tag, a, i) - TH[a]) / TH[a]))

# --- family 3: the frontier is a window, and the census is exact ---------------
EXPECT = {("5", "a"): ("AWAKE", [1, 239]), ("5", "b"): ("AWAKE", [1, 318]),
          ("5", "c"): ("LATE", [1, 1479]), ("5", "d"): ("NEVER", [1]),
          ("1", "a"): ("SLEPT", [1, 240, 963, 1403]),
          ("1", "b"): ("SLEPT", [1, 333, 507, 1530]),
          ("1", "c"): ("NEVER", [1]), ("1", "d"): ("NEVER", [1]),
          ("2", "a"): ("SLEPT", [1, 241, 709, 1439]),
          ("2", "b"): ("SLEPT", [1, 343, 472, 1568]),
          ("2", "c"): ("NEVER", [1]), ("2", "d"): ("NEVER", [1])}
for key, (kind, ev) in sorted(EXPECT.items()):
    require(KIND[key] == kind and EV[key] == ev,
            "n=%s arm %s (tau_0 = %.2f): %s, flag ladder %s -- measured %s %s"
            % (key[0], key[1], dict(ARMS)[key[1]], kind, ev, KIND[key], EV[key]))
counts = {k: sum(1 for v in KIND.values() if v == k)
          for k in ("SLEPT", "AWAKE", "LATE", "NEVER")}
require(counts == {"SLEPT": 4, "AWAKE": 2, "LATE": 1, "NEVER": 5},
        "the 12-cell census is %s: the re-latch is a WINDOW in tau_0 (1.75 and 2.0 "
        "sleep at n = 1 / 2, 2.25 and 2.5 never even wake), not the half-line a "
        "single-threshold sweep could suggest" % counts)
for tag, n, K in CASES:
    slept = [a for a, _ in ARMS if KIND[(tag, a)] == "SLEPT"]
    require(slept == (["a", "b"] if tag != "5" else []),
            "n=%s: the sleeping arms are %s -- at fixed index the sleep is lost on "
            "BOTH sides of the window, by never waking (tau_0 >= 2.25) and by "
            "never being asked (n = 0.5)" % (n, slept or "none"))

# --- family 4: the two axes buy different things -------------------------------
for tag, n, K in CASES:
    vals = [BEST[(tag, a)] for a, _ in ARMS]
    spread = max(vals) - min(vals)
    require(spread < 60.0,
            "n=%s: across tau_0 = 1.75..2.5 the sealed window's best node moves "
            "only %.5f Pa (%s) while the release pressure moves %.0f Pa -- the "
            "refill barely knows the strong wall's yield stress, the ceiling is "
            "made of it" % (n, spread, " / ".join("%.5f" % v for v in vals),
                           REL["a"] - REL["d"]))
NSPREAD = {}
for a, ts in ARMS:
    vals = [BEST[(t, a)] for t, _, _ in CASES]
    spread = max(vals) - min(vals)
    NSPREAD[a] = spread
    require(spread > 2000.0,
            "tau_0 = %.2f: across n = 0.5 / 1 / 2 the best node moves %.5f Pa "
            "(%s) at a FIXED release pressure %.0f Pa -- the index is what buys "
            "the refill%s"
            % (ts, spread, " / ".join("%.5f" % v for v in vals), REL[a],
               "" if a in ("a", "b") else
               ", and here the strong wall never wakes at %d of the three "
               "indices, so all of it is the weak wall's law"
               % sum(1 for t, _, _ in CASES if KIND[(t, a)] == "NEVER")))
require(max(NSPREAD.values()) - min(NSPREAD.values()) < 10.0,
        "and the index buys the SAME refill at every tau_0: the four spreads are "
        "%s Pa, within %.5f Pa of one another -- the strong wall's yield stress "
        "changes what the node must reach, not what the law delivers"
        % (" / ".join("%.5f" % NSPREAD[a] for a, _ in ARMS),
           max(NSPREAD.values()) - min(NSPREAD.values())))
sa = [BEST[(t, "a")] for t, _, _ in CASES]
tau_spread = max(BEST[("5", a)] for a, _ in ARMS) - min(BEST[("5", a)] for a, _ in ARMS)
require((sa[2] - sa[0]) / tau_spread > 10.0,
        "the decomposition is quantitative: n moves the refill %.5f Pa along arm "
        "a while tau_0 moves it %.5f Pa at n = 0.5 -- a ratio of %.1f, against a "
        "release pressure that tau_0 moves by %.0f Pa and n does not move at all"
        % (sa[2] - sa[0], tau_spread, (sa[2] - sa[0]) / tau_spread,
           REL["a"] - REL["d"]))
for tag, n, K in CASES:
    require(BEST[(tag, "a")] > BEST[(tag, "b")] > BEST[(tag, "c")]
            >= BEST[(tag, "d")],
            "n=%s: the best node is monotone DECREASING in the strong wall's "
            "tau_0 (%.5f > %.5f > %.5f >= %.5f) -- a stronger neighbour supplies "
            "less, so it refills the node less, even though it also needs less"
            % (n, BEST[(tag, "a")], BEST[(tag, "b")], BEST[(tag, "c")],
               BEST[(tag, "d")]))

# --- family 5: a positive margin is not a sleep --------------------------------
t, a = "5", "c"
margin = BEST[(t, a)] - REL[a]
wake = EV[(t, a)][1]
require(margin > 0.0 and KIND[(t, a)] == "LATE" and wake > OPENED,
        "n=0.5 arm c: the sealed window's best node %.5f Pa is %.3f Pa ABOVE its "
        "release pressure %.0f, and it STILL never sleeps -- because it woke at "
        "step %d, %d steps AFTER the seal broke at %d. A margin is necessary, not "
        "sufficient: the wall has to be awake inside the window for the desk to "
        "have anything to hold"
        % (BEST[(t, a)], margin, REL[a], wake, wake - OPENED, OPENED))
require(all(R[t][2]["S" + a][i] == 1 for i in range(CLOSED, OPENED + 1)),
        "n=0.5 arm c: all %d states of the sealed window print flag 1 -- the row "
        "was held the whole time, so there was no moving wall to put back to sleep"
        % (OPENED - CLOSED + 1))
for tag, n, K in CASES:
    for a, ts in ARMS:
        if KIND[(tag, a)] != "AWAKE":
            continue
        mind = min(entry(tag, a, i) for i in range(EV[(tag, a)][1], NS))
        require(mind > TH[a],
                "n=%s arm %s: awake from %d and never asked -- the smallest entry "
                "drive after the wake is %.4f Pa, still %.4f Pa above %.0f, so "
                ":40-:43 is never offered a sleep (best node %.5f, margin %+.3f)"
                % (n, a, EV[(tag, a)][1], mind, mind - TH[a], TH[a],
                   BEST[(tag, a)], BEST[(tag, a)] - REL[a]))

# --- family 6: the coast, the sleep window, the second wake --------------------
for tag, n, K in CASES:
    for a, ts in ARMS:
        if KIND[(tag, a)] != "SLEPT":
            continue
        ev = EV[(tag, a)]
        wake, sleep, wake2 = ev[1], ev[2], ev[3]
        cross = next(i for i in range(wake, NS) if entry(tag, a, i) <= TH[a])
        lag = sleep - cross
        q = QL[(tag, "S" + a)]
        require(lag >= 0 and R[tag][2]["S" + a][sleep] == 1,
                "n=%s arm %s: the entry drive falls under %.0f Pa at step %d "
                "(%.4f) and the flag flips %d step(s) later at %d -- the row still "
                "prints %.6g m3/s at the crossing, so under-threshold drive alone "
                "does not sleep a moving wall (:353-:355 also wants the iterate "
                "carried through zero)"
                % (n, a, TH[a], cross, entry(tag, a, cross), lag, sleep,
                   q[cross]))
        require(all(x == 0.0 for x in q[sleep:wake2])
                and all(R[tag][2]["S" + a][i] == 1 for i in range(sleep, wake2)),
                "n=%s arm %s: the sleep window [%d, %d) is %d states of exactly 0 "
                "and flag 1, one contiguous block with no third decision inside"
                % (n, a, sleep, wake2, wake2 - sleep))
        require(entry(tag, a, wake2) > TH[a] and wake2 > OPENED,
                "n=%s arm %s: the second wake at %d is %d steps (%.3f s) after K "
                "returned to 0 at %d, on entry drive %.4f Pa > %.0f -- the "
                "reopened drain buys it, not the rheology"
                % (n, a, wake2, wake2 - OPENED, (wake2 - OPENED) * DT, OPENED,
                   entry(tag, a, wake2), TH[a]))
require([EV[("1", a)][2] - next(i for i in range(EV[("1", a)][1], NS)
                               if entry("1", a, i) <= TH[a]) for a in ("a", "b")]
        == [0, 0]
        and [EV[("2", a)][2] - next(i for i in range(EV[("2", a)][1], NS)
                                   if entry("2", a, i) <= TH[a]) for a in ("a", "b")]
        == [5, 7],
        "the coast is the law's AND the threshold's: at n = 1 both arms flip on the "
        "crossing step itself (lag 0 / 0), at n = 2 tau_0 = 1.75 coasts 5 steps and "
        "tau_0 = 2.0 coasts 7 -- one line of C++ with no n and no tau_0 in it, four "
        "different lags")

# --- family 7: the rows close until the law runs out of root -------------------
SUB = {}
for tag, n, K in CASES:
    N, P, F = R[tag]
    for a, ts in ARMS:
        for w, t0 in (("S" + a, ts), ("W" + a, TAU_W)):
            pn, q, fl = FL[(tag, a)], QL[(tag, w)], F[w]
            mx, worst_at, sub, used = 0.0, None, [], 0
            for i in range(1, NS):
                if fl[i] or fl[i - 1]:
                    continue
                tw = (P_A - pn[i] - I_DT * (q[i] - q[i - 1])) / WALL_SCALE
                if tw <= t0:
                    sub.append(i)
                    continue
                rq = q_hb(tw, t0, n, K)
                rel = abs(rq - q[i]) / rq
                used += 1
                if rel > mx:
                    mx, worst_at = rel, i
            SUB[(tag, w)] = sub
            if used == 0:
                require(all(fl[1:]),
                        "n=%s %s (tau_0 = %.2f): the bracket is offered nothing -- "
                        "all %d decided states print flag 1, so there is no unheld "
                        "step whose row could be compared to the section integral"
                        % (n, w, t0, NS - 1))
                continue
            require(mx < 1.0e-5,
                    "n=%s %s (tau_0 = %.2f): the printed rate IS the HB section "
                    "integral of the wall stress its own row implies at every "
                    "unheld step with a root -- worst %.3e relative at step %s "
                    "over %d of %d decided states"
                    % (n, w, t0, mx, worst_at, used, NS - 1))
for tag, n, K in CASES:
    for a, ts in ARMS:
        key, sub = (tag, a), SUB[(tag, "S" + a)]
        if KIND[key] == "SLEPT":
            require(sub == [EV[key][2] - 1],
                    "n=%s arm %s: exactly ONE unheld step implies a sub-yield wall "
                    "stress -- %d, the last before the clamp at %d. The row runs "
                    "out of root first and :353-:355 follows one step later"
                    % (n, a, sub[0], EV[key][2]))
        else:
            require(sub == [],
                    "n=%s arm %s (%s): no unheld step implies a sub-yield stress "
                    "at all -- the row never ran out of root, so there was nothing "
                    "for the clamp to do" % (n, a, KIND[key]))

# --- family 8: 5.59's arm a in a new file, node closure, event quietness -------
if os.path.isdir(SIB):
    for tag, n, K in CASES:
        NB, PB, FB = run(SIB, "case_%s.dict" % tag)
        N, P, F = R[tag]
        pairs = [("Sb", "Sa", True), ("Wb", "Wa", True), ("Db", "Da", False)]
        for new, old, has_flag in pairs:
            same = sum(1 for x, y in zip(P[new], PB[old]) if x == y)
            md = max(abs(float(x) - float(y)) / max(abs(float(y)), 1e-30)
                     for x, y in zip(P[new], PB[old]))
            fl_ok = (F[new] == FB[old]) if has_flag else True
            fl_txt = ("hold flags equal for all %d states" % NS) if has_flag \
                else "a drain row, no hold flag in either file"
            require(fl_ok and md < 1.0e-7,
                    "n=%s: arm b's %s is 5.59's %s -- %s, %d of %d rate tokens the "
                    "same string, worst %.3e relative. Same physics, new file: only "
                    "the summation order of a 12-pipe network against a 10-pipe one "
                    "pays a decimal"
                    % (n, new, old, fl_txt, same, NS, md))
        same = sum(1 for x, y in zip(N["Fb"], NB["Fa"]) if x == y)
        md = max(abs(float(x) - float(y)) for x, y in zip(N["Fb"], NB["Fa"]))
        require(md <= 2.0e-07 and events(F["Sb"]) == events(FB["Sa"]),
                "n=%s: node Fb is 5.59's Fa within %.3e Pa (%d of %d tokens the "
                "same string) and the ladders are the same integers %s"
                % (n, md, same, NS, events(F["Sb"])))
for tag, n, K in CASES:
    N, P, F = R[tag]
    worst = 0.0
    for a, _ in ARMS:
        pn = [float(x) for x in N["F" + a]]
        for i in range(1, NS):
            bal = (QL[(tag, "S" + a)][i] + QL[(tag, "W" + a)][i]
                   - QL[(tag, "D" + a)][i]) - CC * (pn[i] - pn[i - 1]) / DT
            worst = max(worst, abs(bal))
    require(worst < 1.0e-12,
            "n=%s: all four node rows close on the PRINTED TOKENS alone to %.4e "
            "m3/s = %.2e Pa of node motion per step -- the tau_0 sweep moves money "
            "between rows, it does not create or destroy flux at a node"
            % (n, worst, worst * DT / CC))
    require(all(F[w][0] == 0 for w in F),
            "n=%s: all %d t = 0 hold flags print 0 -- the decision has not happened "
            "yet, so an event scan from i = 0 would invent %d extra flips"
            % (n, len(F), len(F)))
    gaps = [EV[(tag, a)][i + 1] - EV[(tag, a)][i]
            for a, _ in ARMS for i in range(len(EV[(tag, a)]) - 1)]
    require(min(gaps) >= 50,
            "n=%s: the smallest gap between two decisions on any strong wall is %d "
            "steps (%.3f s) over %d gaps -- no dt-scale chatter anywhere in the "
            "12-cell sweep" % (n, min(gaps), min(gaps) * DT, len(gaps)))

print()
print("%d checks, %d failures" % (CHECKS, len(failures)))
for f in failures:
    print("FAILED:", f)
sys.exit(1 if failures else 0)