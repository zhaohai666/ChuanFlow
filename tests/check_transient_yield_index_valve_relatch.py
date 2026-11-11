#!/usr/bin/env python3
"""5.59: the re-latch moved onto an HB wall.

5.54 proved a wall that woke can be put back to sleep by its own node, and left
boundary 3 open: none of its thirteen walls declared a `flowIndex`.  5.57 proved
the desk is index-blind -- threshold `4 tau_0 L / d` (transient.cpp:28-:31), the
`<=` inside `holdsStatic` (:40-:43) -- and 5.58 proved a calendar only reaches
the row downstream of the latch.  Those two together make a sharper question than
"does HB change the band": a sleep needs the NODE to come back up, and the node
is moved by the wall law, which does read `n`.

The kernel has two static-branch shapes and only the second one can sleep a
moving wall (:340-:360):

    held(i)     = qPrev == 0 and |drive| <= 4 tau_0 L / d      (the desk, :40-:43)
    crossing(i) = qPrev != 0 and |drive| <= 4 tau_0 L / d
                  and the Newton iterate carried x[k] through zero
    -> x[k] = 0, r[k] = 0, J[k][k] = 1                          (:353-:355)

so the drive being under the threshold is NECESSARY but not sufficient, and the
half that is missing is the wall's own inertia.  Arm `a` is 5.54's two-way latch
with both walls rewritten as Herschel-Bulkley:

    (inertia/dt) dQ_k + 4L/d * tw_k          = pA - pF     (Sa: tau_0 2.0, Wa: 1.5)
    (inertia/dt) dQ_d + Kd(t) rho|Q_d|Q_d/2A^2 = pF - pB
    (V/BC) dpF   = Q_Sa + Q_Wa - Q_d

with (n, K) = (0.5, 1.0) / (1.0, 0.1) / (2.0, 0.02) in three files, the same
close 0.8 -> 1.0 s to closedK = 1e8 and reopen 2.5 -> 2.7 s, plus `s` (arm a
minus its strong row), `u` (the strong row alone) and `p` (a Newtonian
neighbour).  No kernel change.

  1. what no wall has moved yet is index-blind to the last digit: while every
     wall prints exactly 0 the node is the drain's arithmetic alone, so the weak
     wall releases at step 180 and the lone wall of `u` at 318 in ALL THREE
     tiers, and their node tokens are the same strings until the step whose rate
     stops being zero;
  2. the second release is not: arm a's 8000 Pa wall wakes at 318 / 333 / 343 --
     one threshold, three steps, because between 180 and those steps the node is
     being emptied by the weak wall's HB law;
  3. the sleep needs both halves of :340-:360 and the coast between them is the
     law's: at n = 2 the entry drive falls under 8000 at step 465 and the flag
     only flips at 472 (7 steps, 1.4e-2 s), at n = 1 both land on 507;
  4. at n = 0.5 the sleep is never offered: the node's best rise across the sealed
     window tops out 599 Pa below the release pressure and the entry drive never
     falls under the threshold, while n = 1 / n = 2 overshoot by 1192 / 1544 Pa;
  5. the wall rows close on the HB section integral at every unheld step, and the
     EXACT step before each sleep is the only unheld step whose row implies a
     sub-yield wall stress -- the row runs out of root first, the clamp follows;
  6. `n = 1, K = mu` is 5.54's wall: identical flags, identical event steps
     (1 / 333 / 507 / 1530), node tokens within 1e-07 Pa and rates within
     7.2e-09 relative -- and bit-identical on the arm whose wall never wakes;
  7. a held wall contributes nothing (Wa and Ws are the same strings until the
     strong row first flows), the four node rows close on printed tokens alone,
     and the latch stays event-quiet: no chatter at the dt scale.
"""
import math
import os
import re
import subprocess
import sys

TRANSIENT = sys.argv[1]
CASE = sys.argv[2]
CASE_DIR = CASE.rsplit("/", 1)[0]
BLINK = os.path.join(os.path.dirname(os.path.abspath(CASE_DIR)),
                     "transientYieldValveBranchRelease", "case.dict")

RHO, L, D = 850.0, 20.0, 0.02
R_PIPE, AREA = D / 2.0, math.pi * D * D / 4.0
VOL, BULK = 0.1, 1.0e8
CC = VOL / BULK
WALL_SCALE = 4.0 * L / D
DT = 0.002
I_DT = RHO * L / AREA / DT
P_A, P_B = 3.0e4, 2.0e4
CS, CE, KC, OS_, OE = 0.8, 1.0, 1.0e8, 2.5, 2.7
CLOSED = int(round(CE / DT))                 # first sample with closedK in force
OPENED = int(round(OE / DT))                 # first sample with K back to 0
CEIL_8 = P_A - WALL_SCALE * 2.0              # 22000 Pa: tau_0 = 2.0's release pressure
CEIL_6 = P_A - WALL_SCALE * 1.5              # 24000 Pa
CASES = [("5", 0.5, 1.0), ("1", 1.0, 0.1), ("2", 2.0, 0.02)]
TAU = {"Sa": 2.0, "Wa": 1.5, "Ws": 1.5, "Su": 2.0, "Sp": 2.0}
NODE = {"Sa": "Fa", "Wa": "Fa", "Ws": "Fs", "Su": "Fu", "Sp": "Fp"}
WALLS = tuple(TAU)
FREE = ("Fa", "Fs", "Fu", "Fp")
IN = {"Fa": ("Sa", "Wa"), "Fs": ("Ws",), "Fu": ("Su",), "Fp": ("Sp", "Np")}
OUT = {"Fa": "Da", "Fs": "Ds", "Fu": "Du", "Fp": "Dp"}

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
    """HB volumetric flow at wall stress tw -- the same bracket 5.47/5.57 used."""
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
RB = run(os.path.dirname(BLINK), os.path.basename(BLINK)) \
    if os.path.exists(BLINK) else None
NS = len(R["1"][0]["Fa"])
FL = {(t, w): [float(x) for x in R[t][0][NODE[w]]] for t, _, _ in CASES for w in WALLS}
QL = {(t, w): [float(x) for x in R[t][1][w]] for t, _, _ in CASES for w in WALLS}


def entry(tag, w, i):
    """The drive step i's latch read: the node as step i-1 left it."""
    return P_A - FL[(tag, w)][i - 1]


# --- family 1: what no wall has moved yet is index-blind to the last digit ----
for tag, n, K in CASES:
    N, P, F = R[tag]
    require(events(F["Wa"]) == [1, 180] and events(F["Su"]) == [1, 318]
            and events(F["Sp"]) == [1],
            "n=%s: Wa's flag changes at 1 and 180, arm u's single wall at 1 and "
            "318, arm p's at 1 only -- all three decided while every wall prints "
            "0, i.e. by the Newtonian drain alone, and the %d states of p are a "
            "desk the node never leaves (transient.cpp:40-:43)"
            % (n, NS - 1))
    require(set(P["Wa"][:180]) == {"0"} and set(P["Su"][:318]) == {"0"},
            "n=%s: %d then %d held states print one token 0 before their release "
            "-- the drive the latch compares against 8000 / 6000 Pa contains no "
            "wall term at all in that stretch"
            % (n, 180, 318))
require(all(events(R[t][2]["Wa"]) == events(R["5"][2]["Wa"]) == [1, 180]
            and events(R[t][2]["Su"]) == events(R["5"][2]["Su"]) == [1, 318]
            for t in ("1", "2")),
        "the release STEPS 180 and 318 are the same integers in all three tiers "
        "-- with nothing but a Newtonian drain moving the node the desk is blind "
        "to `n` twice over, in the number and in the step it is reached on")
for nd, upto in (("Fa", 180), ("Fu", 318)):
    cols = [R[t][0][nd] for t in ("5", "1", "2")]
    head = all(c[:upto] == cols[0][:upto] for c in cols[1:])
    diffs = [i for i in range(len(cols[0])) if len({c[i] for c in cols}) > 1]
    require(head and diffs and diffs[0] >= upto,
            "%s: the node tokens are the SAME STRINGS in the three tiers for the "
            "first %d states and first split at %d -- that is the step whose rate "
            "stops being zero, the first step a wall law gets to move the node"
            % (nd, upto, diffs[0] if diffs else -1))
    margin = max(abs(float(cols[a][i]) - float(cols[b][i]))
                 for a in range(3) for b in range(3)
                 for i in range(diffs[0], min(diffs[0] + 40, len(cols[0]))))
    require(margin < 20.0,
            "%s: and the split is small where it starts -- over the 40 states "
            "after it the three tiers stay within %.3f Pa of one another, %.2f%% "
            "of the 8000 Pa the latch compares against, so the law's first print "
            "is a decimal, not a new trajectory"
            % (nd, margin, 100.0 * margin / 8000.0))

# --- family 2: the second release is index-dependent, the threshold is not ----
wake = {}
for tag, n, K in CASES:
    ev = events(R[tag][2]["Sa"])
    wake[tag] = ev[1]
    i = ev[1]
    require(entry(tag, "Sa", i) > 8000.0 >= entry(tag, "Sa", i - 1)
            and QL[(tag, "Sa")][i] > 0.0,
            "n=%s: arm a's 8000 Pa wall wakes at %d on the FIRST entry drive "
            "above 4 tau_0 L/d (%.4f there, %.4f one step earlier) and prints "
            "%.6g m3/s -- by then the node had been emptied since step 180 by the "
            "WEAK wall's HB law" % (n, i, entry(tag, "Sa", i),
                                    entry(tag, "Sa", i - 1), QL[(tag, "Sa")][i]))
require(wake["5"] < wake["1"] < wake["2"],
        "one 8000 Pa crossed at step %s / %s / %s for n = 0.5 / 1 / 2 -- the "
        "threshold (transient.cpp:28-:31) is index-blind and the ARRIVAL is not: "
        "%s s, a spread of %d steps = %.3f s on one number none of them reads"
        % (wake["5"], wake["1"], wake["2"],
           " / ".join("%.3f" % (wake[t] * DT) for t in ("5", "1", "2")),
           wake["2"] - wake["5"], (wake["2"] - wake["5"]) * DT))
for tag, n, K in CASES:
    i = wake[tag]
    d = entry(tag, "Sa", i) - 8000.0
    require(0.0 < d < 80.0,
            "n=%s: the drive at the wake is 8000 + %.4f Pa, i.e. %.3f%% past the "
            "threshold -- three tiers wake within %.1f Pa of one number they "
            "cannot see the inside of" % (n, d, 100.0 * d / 8000.0, 80.0))

# --- family 3: the sleep needs the iterate too, and the coast is the law's ----
seq, cross_step, sleep_step, wake2 = {}, {}, {}, {}
for tag, n, K in CASES:
    ev = events(R[tag][2]["Sa"])
    seq[tag] = ev
    if len(ev) < 3:
        continue
    sleep_step[tag] = ev[2]
    wake2[tag] = ev[3] if len(ev) > 3 else None
    cross_step[tag] = next(i for i in range(wake[tag], NS)
                           if entry(tag, "Sa", i) <= 8000.0)
    lag = sleep_step[tag] - cross_step[tag]
    require(lag >= 0 and R[tag][2]["Sa"][sleep_step[tag]] == 1
            and (lag == 0 or QL[(tag, "Sa")][cross_step[tag]] > 0.0),
            "n=%s: the entry drive falls under 8000 Pa at step %d (%.4f) and the "
            "wall is STILL AWAKE there, printing %.6g m3/s; the flag flips %d "
            "step(s) later at %d -- under-threshold drive alone does not sleep a "
            "moving wall, :353-:355 also wants the iterate carried through zero"
            % (n, cross_step[tag], entry(tag, "Sa", cross_step[tag]),
               QL[(tag, "Sa")][cross_step[tag]], lag, sleep_step[tag]))
require(sleep_step["1"] - cross_step["1"] == 0
        and sleep_step["2"] - cross_step["2"] == 7,
        "the coast is 0 steps at n = 1 and 7 steps (1.4e-02 s) at n = 2 -- one "
        "line of C++ (`x[k] * pd.qPrev <= 0.0`) with no `n` in it, timed entirely "
        "by how much momentum each law had put into the column")
for tag, n, K in CASES:
    if tag not in sleep_step:
        continue
    i, j = sleep_step[tag], wake2[tag]
    q, fl = QL[(tag, "Sa")], R[tag][2]["Sa"]
    require(set(q[i:j]) == {0.0} and sum(fl[i:j]) == j - i,
            "n=%s: the slept wall prints exactly 0 in all %d states of [%d, %d) "
            "-- the sleep is absolute, there is no sub-yield creep to print"
            % (n, j - i, i, j))
    require(entry(tag, "Sa", j) > 8000.0 >= entry(tag, "Sa", j - 1)
            and j > OPENED,
            "n=%s: the same wall wakes again at %d (entry drive %.4f over 8000, "
            "previous %.4f) -- %.2f s after K returned to exactly 0 at step %d, "
            "so the reopened drain bought the second release, not the law"
            % (n, j, entry(tag, "Sa", j), entry(tag, "Sa", j - 1),
               (j - OPENED) * DT, OPENED))
    require(all(fl[k] == 1 for k in range(i, j)),
            "n=%s: the hold is one continuous block of %d states between the "
            "sleep at %d and the second wake at %d -- no intermediate decision "
            "re-opens it" % (n, j - i, i, j))
require(seq["1"] == [1, 333, 507, 1530] and seq["2"] == [1, 343, 472, 1568],
        "the event ladders are 1/333/507/1530 at n = 1 and 1/343/472/1568 at "
        "n = 2 -- the shear-thickening wall wakes 10 steps later, sleeps 35 steps "
        "earlier and wakes again 38 steps later, all four on the same 8000 Pa")
require(seq["5"] == [1, 318],
        "n=0.5: the flag changes twice and never again -- n = 0.5 has no sleep in "
        "this window at all, which is the subject of the next family")

# --- family 4: at n = 0.5 the sleep is never offered, and the node says why ------
best, gap = {}, {}
for tag, n, K in CASES:
    N, P, F = R[tag]
    pn = [float(x) for x in N["Fa"]]
    best[tag] = max(pn[CLOSED:OPENED + 1])
    gap[tag] = CEIL_8 - best[tag]
    sup = sum(float(P[w][CLOSED]) for w in ("Sa", "Wa"))
    leak = float(P["Da"][CLOSED])
    require(sup > leak and sup / leak > 1.5,
            "n=%s: at the first sealed state (step %d, closedK = %.0e) the two "
            "walls supply %.6g m3/s against the valve's %.6g -- ratio %.3f, so "
            "every tier starts winning the race and the node rises from %.5f Pa"
            % (n, CLOSED, KC, sup, leak, sup / leak, pn[CLOSED]))
require(gap["5"] > 0.0 > gap["1"] > gap["2"] - 400.0,
        "the sealed window's best node value is %.5f / %.5f / %.5f Pa against the "
        "release pressure 22000 -- n = 0.5 tops out %.3f Pa SHORT and is never "
        "offered a sleep, n = 1 overshoots by %.3f Pa and n = 2 by %.3f Pa "
        "(n = 2 overshot so early it slept before the seal completed)"
        % (best["5"], best["1"], best["2"], gap["5"], -gap["1"], -gap["2"]))
drives5 = [entry("5", "Sa", i) for i in range(wake["5"], NS)]
require(min(drives5) > 8000.0,
        "n=0.5: the entry drive never once falls to 4 tau_0 L/d after the wake -- "
        "its smallest is 8000 + %.4f Pa (%.3f%% of the threshold), so the desk is "
        "not declining a sleep, it is never asked; the %.3f Pa the wall law misses "
        "the release pressure by is all there is"
        % (min(drives5) - 8000.0, (min(drives5) - 8000.0) / 8000.0 * 100.0, gap["5"]))
for tag, n, K in CASES:
    pn = [float(x) for x in R[tag][0]["Fa"]]
    rise = best[tag] - pn[CLOSED]
    require(rise > 0.0,
            "n=%s: from the first sealed state the node climbs %.4f Pa on the "
            "walls' supply before the reopened drain takes over -- the refill "
            "exists at every index, only its size knows `n`" % (n, rise))
q5 = max(QL[("5", "Wa")][CLOSED:OPENED])
q1 = max(QL[("1", "Wa")][CLOSED:OPENED])
q2 = max(QL[("2", "Wa")][CLOSED:OPENED])
require(q5 < q1 and q5 < q2,
        "the weak wall alone peaks at %.6g m3/s across the sealed window at "
        "n = 0.5 against %.6g (n = 1) and %.6g (n = 2) -- one order of magnitude "
        "of supply is the difference between a sleep and none" % (q5, q1, q2))

# --- family 5: the rows close until the law runs out of root, then the clamp ---
for tag, n, K in CASES:
    N, P, F = R[tag]
    for w in WALLS:
        tau0 = TAU[w]
        pn, q, fl = FL[(tag, w)], QL[(tag, w)], F[w]
        mx, worst_at, sub, used = 0.0, None, [], 0
        for i in range(1, NS):
            if fl[i] or fl[i - 1]:
                continue
            tw = (P_A - pn[i] - I_DT * (q[i] - q[i - 1])) / WALL_SCALE
            if tw <= tau0:
                sub.append(i)
                continue
            rq = q_hb(tw, tau0, n, K)
            rel = abs(rq - q[i]) / rq
            used += 1
            if rel > mx:
                mx, worst_at = rel, i
        if used == 0:
            require(all(fl[1:]),
                    "n=%s: %s offers the bracket nothing at all -- all %d decided "
                    "states print flag 1, so there is no unheld step after t = 0 "
                    "whose row could be compared to the section integral. Arm p is "
                    "the desk 5.57 predicted and family 1 measured "
                    "(transient.cpp:40-:43)" % (n, w, NS - 1))
            continue
        require(mx < 1.0e-5,
                "n=%s: %s's printed rate IS the HB section integral of the wall "
                "stress its own row implies (drive minus (inertia/dt) dQ, no "
                "valve on this row) at every unheld step with a root -- worst "
                "%.3e relative at step %s over %d of %d decided states; the held "
                "steps are excluded because there the row was replaced by x=0"
                % (n, w, mx, worst_at, used, NS - 1))
        if w == "Sa" and tag in sleep_step:
            require(sub == [sleep_step[tag] - 1],
                    "n=%s: exactly ONE unheld step in the whole run implies a "
                    "sub-yield wall stress -- step %d, the last one before the "
                    "clamp at %d. The row runs out of root first and :353-:355 "
                    "follows one step later; holding is not a residual"
                    % (n, sub[0], sleep_step[tag]))
        elif w == "Sa":
            require(sub == [],
                    "n=%s: no unheld step of Sa implies a sub-yield stress at all "
                    "-- the row never ran out of root, because the drive never "
                    "came back down under 8000 Pa (family 4)" % n)

# --- family 6: n = 1, K = mu is 5.54's wall, not a look-alike ------------------
if RB:
    NB, PB, FB = RB
    for w in ("Sa", "Wa", "Su", "Sp", "Np", "Da", "Du", "Dp"):
        if w not in PB or w not in R["1"][1]:
            continue
        eq = sum(1 for a, b in zip(PB[w], R["1"][1][w]) if a == b)
        md = max(abs(float(a) - float(b)) for a, b in zip(PB[w], R["1"][1][w]))
        rd = max(abs(float(a) - float(b)) / max(abs(float(a)), 1e-300)
                 for a, b in zip(PB[w], R["1"][1][w]) if a != b) if md else 0.0
        have = w in FB and w in R["1"][2]
        fl_eq = (not have) or FB[w] == R["1"][2][w]
        first = next((i for i, (a, b) in enumerate(zip(PB[w], R["1"][1][w]))
                      if a != b), None)
        fl_txt = ("hold flags equal for all %d states" % NS) if have else \
            ("%s is a drain row, no hold flag in either file" % w)
        require(fl_eq and (md == 0.0 or rd < 1.0e-7),
                "5.54's %s vs this file's n=1 %s: %s, %d/%d rate tokens identical, "
                "worst absolute diff %.1e m3/s (%.1e relative) first showing at "
                "step %s -- (n, K) = (1, mu) is the same row, and only the last "
                "printed decimal knows it was written twice"
                % (w, w, fl_txt, eq, NS, md, rd, first))
    for nd in ("Fa", "Fu", "Fp"):
        md = max(abs(float(a) - float(b)) for a, b in zip(NB[nd], R["1"][0][nd]))
        eq = sum(1 for a, b in zip(NB[nd], R["1"][0][nd]) if a == b)
        wl = {"Fa": "Sa", "Fu": "Su", "Fp": "Sp"}[nd]
        require(md <= 2.0e-07 and events(FB[wl]) == events(R["1"][2][wl]),
                "5.54's node %s and this file's n=1 %s agree to %.1e Pa over %d "
                "states (%d of %d tokens are the same string) and their %s ladders "
                "are the same integers %s -- one printed decimal of a node "
                "pressure is 1e-07 Pa, 1.2e-09 of the 8000 Pa the latch reads"
                % (nd, nd, md, NS, eq, NS, wl, events(FB[wl])))
    require(events(FB["Sa"]) == events(R["1"][2]["Sa"]) == [1, 333, 507, 1530]
            and events(FB["Su"]) == events(R["1"][2]["Su"]),
            "5.54's Bingham two-way latch and this file's HB (n=1, K=mu) latch "
            "flip at identical steps (%s), and so does the single-wall arm (%s) "
            "-- the whole event ladder is one formula's two spellings"
            % (events(FB["Sa"]), events(FB["Su"])))
    require(events(FB["Sp"]) == events(R["1"][2]["Sp"])
            and set(PB["Sp"]) == {"0"} == set(R["1"][1]["Sp"]),
            "on arm p, whose wall never wakes, 5.54's and this file's tokens are "
            "bit-identical for all %d states (max diff 0.0) -- the arithmetic "
            "that touches no law reproduces exactly; only flowing states pay the "
            "decimal" % NS)
else:
    require(True, "5.54's case file was not found next to this one (%s); the "
                  "cross-increment identity family was skipped" % BLINK)

# --- family 7: a held wall is silent, the nodes close, the latch is quiet ------
for tag, n, K in CASES:
    N, P, F = R[tag]
    upto = wake[tag]
    ew = sum(1 for a, b in zip(P["Wa"][:upto], P["Ws"][:upto]) if a == b)
    en = sum(1 for a, b in zip(N["Fa"][:upto], N["Fs"][:upto]) if a == b)
    first = next((i for i in range(NS) if P["Wa"][i] != P["Ws"][i]), None)
    require(ew == upto and en == upto,
            "n=%s: `s` is `a` with the strong row taken out, and until that row "
            "first flows (step %d) the weak wall's %d rate tokens and the node's "
            "%d tokens are the same strings (they part at %s) -- a held wall "
            "contributes exactly nothing to node continuity"
            % (n, upto, upto, upto, first))
worst = 0.0
for tag, n, K in CASES:
    N, P, F = R[tag]
    for nd in FREE:
        pn = [float(x) for x in N[nd]]
        for i in range(1, NS):
            r = (sum(float(P[w][i]) for w in IN[nd])
                 - float(P[OUT[nd]][i]) - CC * (pn[i] - pn[i - 1]) / DT)
            worst = max(worst, abs(r))
require(worst < 1.0e-12,
        "all twelve node rows (4 free nodes x 3 tiers) close on the printed "
        "tokens alone across every release, sleep and second release -- worst "
        "%.4e m3/s, i.e. %.2e Pa of node motion per step "
        "(transient.cpp:395-:416)" % (worst, worst / CC * DT))
for tag, n, K in CASES:
    N, P, F = R[tag]
    ev = events(F["Sa"])
    q, fl = QL[(tag, "Sa")], F["Sa"]
    quiet = all(ev[i] - ev[i - 1] >= 50 for i in range(2, len(ev)))
    awake = [i for i in range(ev[1] + 1, len(q)) if not fl[i]]
    zeros = [i for i in awake if q[i] == 0.0]
    nz = min(abs(q[i]) for i in awake if q[i] != 0.0)
    want = [sleep_step[tag] - 1] if tag in sleep_step else []
    require(quiet and zeros == want,
            "n=%s: no two latch decisions after the first are closer than 50 "
            "steps (0.1 s), and the ONLY awake state printing exactly 0 is %s "
            "-- the desk and the printer agree one step apart: the row runs out "
            "of root, its rate falls under the print quantum (smallest nonzero "
            "elsewhere %.3e m3/s), then the clamp writes the same 0 with a flag "
            "on it" % (n, zeros if zeros else "nowhere in %d awake states"
                       % len(awake), nz))
    has_sleep = len(ev) > 3
    ok_tail = all(fl[i] == 0 for i in range(ev[-1], NS))
    ok_sleep = all(q[i] == 0.0 for i in range(ev[2], ev[3])) if has_sleep else True
    require(ok_tail and ok_sleep,
            ("n=%s: after the last decision (%d) the wall stays awake to t = %.1f s "
             "and the sleep window [%d, %d) is all zeros -- the flag and the token "
             "say the same thing in both directions"
             % (n, ev[-1], (NS - 1) * DT, ev[2], ev[3])) if has_sleep else
            ("n=%s: after its single decision at %d the wall stays awake for the "
             "remaining %d states, and every one of them carries flag 0"
             % (n, ev[-1], NS - 1 - ev[-1])))
require(all(R[tag][2][w][0] == 0 for tag, _, _ in CASES for w in WALLS),
        "t = 0 prints a hold flag of 0 on all %d yielded walls of all %d tiers "
        "(%d flags), so an event scan still has to start at i = 1"
        % (len(WALLS), len(CASES), len(WALLS) * len(CASES)))
for tag, n, K in CASES:
    N, P, F = R[tag]
    pn = [float(x) for x in N["Fp"]]
    require(max(pn) > CEIL_6 and all(float(x) > CEIL_8 for x in pn),
            "n=%s: arm p's node never even reaches 22000 Pa (its minimum is "
            "%.5f) because the Newtonian neighbour's flow does not vanish at the "
            "wall's own release pressure -- the ceiling is a property of the "
            "neighbour's law, and it holds at every index" % (n, min(pn)))

print("\n%d checks, %d failures" % (CHECKS, len(failures)))
sys.exit(1 if failures else 0)
