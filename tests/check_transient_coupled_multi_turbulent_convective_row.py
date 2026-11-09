#!/usr/bin/env python3
"""5.43: the shape account moves from the block's chord into the network's own
momentum row -- and is charged once, not twice.

    Src -> A -> { X1 (3D kOmegaSST duct, parallel, rate free)
                  || By (1D pipe) } -> X2 (3D RANS duct, series, q pinned)
         -> X2 -> B -> Snk

Four arms, one variable: `convectiveToRow` says per arm and per edge whether the
region's own changing profile is paid by the 1D row (`charge`) or only reported
(`report`, which is what every case until now did).

    c_off  -- report on both edges (= 5.41's `p_net`, and the inertness control)
    c_both -- charge on both edges
    c_pin  -- charge on X2 only, the edge whose rate the topology pins
    c_free -- charge on X1 only, the edge whose rate the split relaxes

What the kernel promises, and what this script checks:

  1. a charged edge's row is
         rowDrop = R*q + M*dq/dt + dpConv
     and the interface now iterates R against
         chord   = dpFriction - dpConv
     so the two sides still add up to the block's own friction half: the account
     moved, it did not double. A reported edge's rowBeta is exactly 0.
  2. nothing is charged by default: `c_off` reproduces 5.41's `p_net` edge line at
     the start and at the end of the trajectory, in the printed fields and in the
     coefficients, and it needs the same number of interface passes at every one of
     the 8 levels -- with no shape term anywhere, which is what makes the whole
     thing a test of the switch rather than of the case.
  3. the relabeling is computable: R_charged / R_off == 1 - dpConv/dpFriction, so
     the coefficient the drawing used to declare stops being the wall's number by
     9.5..90.5% on X1 and 12.5..39.4% on X2 -- the edge whose shape account is
     closest to its own friction is the one whose resistance moves most.
  4. on the pinned edge the charge is a relabeling and nothing else: the rate is
     6.0e-5 to the last digit in all four arms at all eight levels, and the only
     fields that move are the three the split touches (R, dpNetwork, rowBeta).
     Everything else that moves moves by 3.5e-2 Pa against an account of 283 Pa.
  5. the choice on one edge does not reach the other: c_pin's X1 equals c_off's
     X1 in all 17 fields and c_both's X1 equals c_free's X1 in all 17, so an
     edge's row answers to its own declaration while its neighbour's is free to
     disagree -- the same per-edge falsifiability 5.35 bought for attribution.
  6. on the free edge the charge is a different trajectory: X1's rate moves
     1.42e-4 relative (max, at level 8), its block's drop 6.05e-4 and the split
     of the coefficient 90.5%, while the row's TOTAL moves 1.25e-1 Pa against
     the 1.489e2 Pa moved between its two sides -- 8.4e-4 of the account.
  7. 5.40's two brackets survive the split unchanged in structure:
         rowError*dpBlock == (dpNetwork - chord) + (M*dq/dt - dpInertia)
     closes to 6.8e-13 over all 64 states, with the chord bracket capped by the
     interface tolerance at 5.9e-1 Pa and the attribution bracket at 754.24 Pa
     (X1) and 1903.42 Pa (X2) -- three and a half to four orders apart, whatever
     the case decided about the shape term.
  8. the account itself is a property of the region and not of who pays it:
     dpConv/q^2 agrees across all four arms to 2.4e-5 relative on X1 and to
     1.0e-4 on X2, at every level.

usage: check_transient_coupled_multi_turbulent_convective_row.py <solver> <case>
"""
import math
import os
import re
import subprocess
import sys

RHO = 1000.0
Q_IN = 6.0e-5
DT, N_STEPS = 0.2, 8
TOL = 1e-3
ARMS = ("c_off", "c_both", "c_pin", "c_free")
CHARGE = {"c_off": ("report", "report"), "c_both": ("charge", "charge"),
          "c_pin": ("report", "charge"), "c_free": ("charge", "report")}
NAME = ("X1", "X2")
SEED = {"X1": 1.4e7, "X2": 1.36e7}
FIELDS = ("q", "dpBlock", "dpNetwork", "R", "residual", "M", "dpFriction",
          "dpInertia", "qPrev", "dq", "rowDrop", "rowError", "betaIn", "betaOut",
          "pRisePa", "dpConv", "rowBeta")
# 5.41's `p_net` level-1 states, measured in that case's own log and re-asserted
# there; quoted here so the inertness of the new switch is graded against the
# shipped number rather than against a copy of this run.
PINET_L1 = {"X1": {"dpBlock": 851.08258956984685, "dpFriction": 579.40211914132692,
                   "dpNetwork": 579.40721946887868, "dpConv": -55.239173285473818,
                   "rowError": 0.88621650000126295},
            "X2": {"dpBlock": 2868.6479781972207, "dpFriction": 2225.5865729289981,
                   "dpNetwork": 2225.5627280310264, "dpConv": -277.29484678803749,
                   "rowError": 0.66351600258052723}}
# ... and the same at the end of the trajectory, so the inertness is graded on where
# the run starts and where it lands rather than on one level alone.
PINET_L8 = {"X1": {"q": 3.1118825189579614e-05, "dpBlock": 177.93872006380388,
                   "dpNetwork": 164.74362239216259, "dpFriction": 164.77917883517469,
                   "R": 5294018.0546188587, "rowDrop": 212.60041875660994,
                   "dpConv": -148.97387627505432, "rowError": 0.19479570652400632},
            "X2": {"q": 6.0000000000000002e-05, "dpBlock": 720.44658376591372,
                   "dpNetwork": 720.4249438241302, "dpFriction": 720.44848007364146,
                   "R": 12007082.397068836, "rowDrop": 720.4249438241302,
                   "dpConv": -283.2722617002043, "rowError": -3.0036844189617738e-05}}
# The interface's pass sequence for that same 5.41 run -- how hard it worked at each
# level, which is the part of the behaviour a printed field cannot show.
PINET_PASSES = (3, 6, 5, 5, 4, 4, 4, 4)


def value(line, key):
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([-+0-9.eE]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return float(m.group(1))


def word(line, key):
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([A-Za-z0-9_]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return m.group(1)


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("transient shape-charge solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    arms, cur = [], None
    for ln in lines:
        t = ln.split()
        if ln.startswith("transient start ") and " accelerate " in ln:
            cur = {"arm": t[2].rstrip(":"), "conv": {}, "shape": {}, "inertia": {},
                   "seed": {}, "edge": {}, "pass": {}, "pipe": {}, "solves": 0}
            arms.append(cur)
        elif cur is None:
            continue
        if ln.startswith("SIMPLE iter 1:"):
            cur["solves"] += 1
        elif ln.startswith("transient arm ") and " convective " in ln:
            cur["conv"][t[4]] = t[6]
        elif ln.startswith("transient arm ") and " profile " in ln:
            cur["shape"][t[4]] = t[6]
        elif ln.startswith("transient arm ") and " inertia " in ln:
            cur["inertia"][t[4]] = t[6]
        elif ln.startswith("transient start ") and " resistance " in ln:
            cur["seed"][word(ln, "edge")] = value(ln, "resistance")
        elif ln.startswith("transient step ") and " pass " in ln:
            cur["pass"].setdefault(int(t[2]), []).append(ln)
        elif ln.startswith("transient step ") and " edge " in ln:
            cur["edge"].setdefault(int(t[2]), {})[t[6]] = ln
        elif ln.startswith("pipes step "):
            cur["pipe"][int(t[2].rstrip(":"))] = ln

    by = {c["arm"]: c for c in arms}

    def E(a, n, i):
        return by[a]["edge"][n][NAME[i]]

    charged = lambda a, i: by[a]["conv"][NAME[i]] == "charge"
    ALL = [(a, n, i) for a in ARMS for n in range(1, N_STEPS + 1) for i in (0, 1)]

    # ===== 1. one axis: four arms, everything else held
    header = [ln for ln in lines if ln.startswith("transient network+blocks")]
    require(len(header) == 1 and "2 coupled edges, 4 starting sets" in header[0]
            and "shape charged on some edges" in header[0],
            f"one network, two coupled edges, four starting sets, and the header "
            f"announces that the shape account is charged somewhere")
    require([c["arm"] for c in arms] == list(ARMS),
            f"the four arms ran in the order the case declares them (got "
            f"{[c['arm'] for c in arms]})")
    require(all(c["seed"] == SEED for c in arms),
            f"every arm seeds both edges with the same pair ({SEED['X1']:.3e}, "
            f"{SEED['X2']:.3e}), so the charge cannot come in through the starting "
            "coefficients")
    require(all(c["shape"] == {"X1": "parabolic", "X2": "parabolic"} for c in arms),
            "every arm is shaped on every edge: the shape axis is held, not crossed")
    require(all(c["inertia"] == {"X1": "network", "X2": "network"} for c in arms),
            "every arm keeps the drawing's rho*L/A: attribution is held still while "
            "the shape account moves")
    require([tuple(by[a]["conv"][n] for n in NAME) for a in ARMS]
            == [CHARGE[a] for a in ARMS],
            "the charge table is read back per arm and per edge from the log, and "
            f"matches the case: {[CHARGE[a] for a in ARMS]}")
    require(all(int(value(E(a, n, i), "unsteady")) == 1 for a in ARMS
                for n in range(1, N_STEPS + 1) for i in (0, 1)),
            "all 64 (arm x level x edge) states are stepping: a quasi-steady partner "
            "has no friction half to take the account out of, so this case does not "
            "ask one")

    # ===== 2. the row carries it, and carries it exactly
    beta_bad = [(a, n, NAME[i], value(E(a, n, i), "rowBeta"),
                 value(E(a, n, i), "dpConv"))
                for (a, n, i) in ALL
                if abs(value(E(a, n, i), "rowBeta")
                       - (value(E(a, n, i), "dpConv") if charged(a, i) else 0.0))
                > 1e-9]
    require(not beta_bad,
            f"rowBeta is the whole account on a charged edge and exactly 0 on a "
            f"reported one in all 64 states ({len([1 for x in ALL if charged(x[0], x[2])])} "
            f"charged, {len([1 for x in ALL if not charged(x[0], x[2])])} reported)"
            + (f"; first mismatch {beta_bad[0]}" if beta_bad else ""))
    row_bad = max(abs(value(E(a, n, i), "rowDrop")
                      - (value(E(a, n, i), "dpNetwork")
                         + value(E(a, n, i), "M") * value(E(a, n, i), "dq") / DT
                         + value(E(a, n, i), "rowBeta"))) for (a, n, i) in ALL)
    require(row_bad == 0.0,
            f"the printed row closes as dpNetwork + M*dq/dt + rowBeta with no slack "
            f"at all in any state (worst {row_bad:.1e} Pa): the third term is in the "
            "equation, not next to it")
    two_term = all(abs(value(E(a, n, i), "rowDrop")
                       - (value(E(a, n, i), "dpNetwork")
                          + value(E(a, n, i), "M") * value(E(a, n, i), "dq") / DT))
                   <= 1e-9 for (a, n, i) in ALL if not charged(a, i))
    require(two_term,
            "and on a reported edge the same row is still the two-term row 5.41 "
            "solved -- both statements come from the same field of the log")

    # ===== 3. the chord gave it up: the interface's own criterion moved with it
    # residual is what the kernel iterates on, so re-deriving it from the printed
    # fields against the REDUCED chord is the double-count test: had the chord kept
    # the account, this gap would be the account itself (283 Pa on X2, 39% of it).
    gap_bad, gaps = [], {}
    for (a, n, i) in ALL:
        l = E(a, n, i)
        target = value(l, "dpFriction") - value(l, "rowBeta")
        rel = abs(value(l, "dpNetwork") - target) / abs(target)
        gaps[(a, i)] = max(gaps.get((a, i), 0.0), rel)
        if abs(rel - value(l, "residual")) > 1e-12:
            gap_bad.append((a, n, NAME[i], rel, value(l, "residual")))
    require(not gap_bad,
            "the kernel's own residual is measured against dpFriction - rowBeta in "
            "every state, charged or not, so the shape term cannot be paid by both "
            "the resistance and the row")
    require(all(v <= TOL * (1.0 + 1e-9) for v in gaps.values()),
            f"every arm's worst relative chord gap is "
            f"{max(gaps.values()):.3e}, inside the case's own tolerance {TOL:.0e} "
            f"(c_off X1 {gaps[('c_off', 0)]:.3e}, c_both X1 {gaps[('c_both', 0)]:.3e}, "
            f"c_both X2 {gaps[('c_both', 1)]:.3e}) -- a charged edge converges on a "
            "smaller chord, so the same absolute agreement is harder to reach")
    absGap = {(a, i): max(abs(value(E(a, n, i), "dpNetwork")
                               - (value(E(a, n, i), "dpFriction")
                                  - value(E(a, n, i), "rowBeta")))
                          for n in range(1, N_STEPS + 1))
              for (a, i) in gaps}
    acct = {(a, i): max(abs(value(E(a, n, i), "rowBeta"))
                        for n in range(1, N_STEPS + 1))
            for (a, i) in gaps if charged(a, i)}
    require(all(absGap[k] < 0.01 * acct[k] for k in acct),
            f"per charged edge the leftover is at most {max(absGap[k] for k in acct):.3f} Pa "
            f"against an account of {min(acct.values()):.1f} Pa moved between the two "
            f"sides ({max(absGap[k] / acct[k] for k in acct) * 100:.2f}%): the split is "
            "reported as a relabeling and graded as one")

    # ===== 4. behaviour inertness: c_off IS 5.41's p_net
    l1 = [value(E("c_off", 1, i), "dpBlock") for i in (0, 1)]
    bad541 = [(NAME[i], f, value(E("c_off", 1, i), f), PINET_L1[NAME[i]][f])
              for i in (0, 1) for f in PINET_L1["X1"]
              if abs(value(E("c_off", 1, i), f) - PINET_L1[NAME[i]][f]) > 1e-12]
    require(not bad541,
            f"c_off level 1 reproduces 5.41's p_net on both edges to the last printed "
            f"digit ({l1[0]:.8f} against the shipped {PINET_L1['X1']['dpBlock']:.8f} "
            f"and {l1[1]:.8f} against {PINET_L1['X2']['dpBlock']:.8f}, and the "
            "coefficients, frictions and row errors with it)")
    bad8 = [(NAME[i], f, value(E("c_off", 8, i), f), PINET_L8[NAME[i]][f])
            for i in (0, 1) for f in PINET_L8["X1"]
            if abs(value(E("c_off", 8, i), f) - PINET_L8[NAME[i]][f]) > 1e-12]
    require(not bad8,
            "and the same at the far end of the trajectory: level 8 lands on 5.41's "
            f"X1 rate {value(E('c_off', 8, 0), 'q'):.6e} and row totals "
            f"{value(E('c_off', 8, 0), 'rowDrop'):.4f} / "
            f"{value(E('c_off', 8, 1), 'rowDrop'):.4f} Pa, so nothing drifted on the "
            "way either")
    pOff = [len(by["c_off"]["pass"][n]) for n in range(1, N_STEPS + 1)]
    require(tuple(pOff) == PINET_PASSES,
            f"the interface worked exactly as hard as it did there: passes per level "
            f"{pOff} against the shipped {list(PINET_PASSES)} -- a switch that was "
            "installed but never read cannot change how many passes a level needs")
    require(all(value(E("c_off", n, i), "rowBeta") == 0.0 for n in range(1, 9)
                for i in (0, 1)),
            "and it carries no shape term at all, which is what makes the last "
            "assertion a test of the switch rather than of the case")

    # ===== 5. the relabeling is computable from the control
    relR, shift = {}, {}
    for a in ("c_both", "c_pin", "c_free"):
        for i in (0, 1):
            if not charged(a, i):
                continue
            w = max(abs(value(E(a, n, i), "R") / value(E("c_off", n, i), "R")
                        - (1.0 - value(E("c_off", n, i), "dpConv")
                           / value(E("c_off", n, i), "dpFriction")))
                    for n in range(1, N_STEPS + 1))
            relR[(a, i)] = w
            shift[(a, i)] = (min(value(E(a, n, i), "R") / value(E("c_off", n, i), "R")
                                 for n in range(1, 9)),
                             max(value(E(a, n, i), "R") / value(E("c_off", n, i), "R")
                                 for n in range(1, 9)))
    require(max(relR.values()) < 2e-3,
            f"R_charged/R_report equals 1 - dpConv/dpFriction at every level and on "
            f"both edges (worst {max(relR.values()):.1e}), so the number the drawing "
            "declared is now the wall's friction alone")
    require(shift[("c_free", 0)][1] > 1.5 and shift[("c_pin", 1)][1] < 1.5,
            f"and the two edges move by very different amounts: X1's coefficient "
            f"grows {100*(shift[('c_free',0)][0]-1):.1f}% to {100*(shift[('c_free',0)][1]-1):.1f}% "
            f"while X2's grows {100*(shift[('c_pin',1)][0]-1):.1f}% to "
            f"{100*(shift[('c_pin',1)][1]-1):.1f}% -- the shape account is 90% of X1's "
            "friction and 39% of X2's")

    # ===== 6. the pinned edge, and where the charge does not reach
    moved = {f: max(abs(value(E("c_pin", n, 1), f) - value(E("c_off", n, 1), f))
                    for n in range(1, N_STEPS + 1)) for f in FIELDS}
    triple = {"R", "dpNetwork", "rowBeta"}
    phys = [f for f in FIELDS if f not in triple | {"residual"}]
    worstPhys = max((moved[f], f) for f in phys)
    require(len({value(E(a, n, 1), "q") for a in ARMS for n in range(1, 9)}) == 1
            and abs(value(E("c_off", 1, 1), "q") - Q_IN) < 1e-20,
            f"X2's rate is the injected 6.0e-5 to the last digit in all four arms at "
            "all eight levels: the topology pins it, so the charge cannot move it")
    require(worstPhys[0] < 0.05,
            f"of the fields that are not the split's own three, the largest movement "
            f"on the pinned edge is {worstPhys[1]} at {moved[worstPhys[1]]:.3e} -- "
            f"against rowBeta {moved['rowBeta']:.1f} and R {moved['R']:.3e}: the rest "
            "is the interface's stopping noise, not the account")
    require(moved["rowDrop"] < 0.05,
            f"the pinned edge's row total moves {moved['rowDrop']:.4f} Pa out of the "
            f"{moved['rowBeta']:.1f} Pa its row took on: the destination is where it was")
    other = {f: max(abs(value(E("c_pin", n, 0), f) - value(E("c_off", n, 0), f))
                    for n in range(1, N_STEPS + 1)) for f in FIELDS}
    require(all(v == 0.0 for v in other.values()),
            f"an edge's row answers only to its own declaration: charging X2 leaves the "
            f"untouched X1 bit-identical to the reporting arm in all {len(FIELDS)} "
            "printed fields at all 8 levels")
    same = {f: max(abs(value(E("c_both", n, 0), f) - value(E("c_free", n, 0), f))
                   for n in range(1, N_STEPS + 1)) for f in FIELDS}
    require(all(v == 0.0 for v in same.values()),
            f"and the same holds with the shared edge already charged: X1, charged in "
            f"c_both and in c_free alike, is bit-identical across them too -- the table "
            "is read per edge, never per arm")

    # ===== 7. the free edge: a real trajectory, a small one
    qMove = [abs(value(E("c_free", n, 0), "q") / value(E("c_off", n, 0), "q") - 1.0)
             for n in range(1, N_STEPS + 1)]
    dpMove = [abs(value(E("c_free", n, 0), "rowDrop") - value(E("c_off", n, 0), "rowDrop"))
              for n in range(1, N_STEPS + 1)]
    acctX1 = max(abs(value(E("c_free", n, 0), "rowBeta")) for n in range(1, 9))
    require(max(qMove) > 1e-5,
            f"charging the FREE edge does move it: X1's rate differs from c_off by up "
            f"to {100*max(qMove):.4f}% (level {qMove.index(max(qMove))+1}), because the "
            "row it re-solves now has a term quadratic in the rate")
    require(max(dpMove) < 1e-3 * acctX1,
            f"and the row's total moves {max(dpMove):.4f} Pa against the {acctX1:.1f} Pa "
            f"({100*max(dpMove)/acctX1:.3f}% of it) it took from the chord: the split "
            "buys a different path, not a different bill")
    qBy = lambda a, n: value(by[a]["pipe"][n], "By")
    by_split = max(abs(qBy(a, 8) - qBy("c_off", 8)) for a in ("c_free", "c_both"))
    require(1e-12 < by_split < 1e-8,
            f"the branch takes the rest: at level 8 By carries {qBy('c_free', 8):.6e} "
            f"where the reporting arm carries {qBy('c_off', 8):.6e} (difference "
            f"{by_split:.3e} m3/s) -- the 1D pipe gets the remainder without being "
            "asked, which is why its rate moves at all")
    closed = max(abs(value(by[a]["pipe"][n], "X1") + qBy(a, n) - Q_IN)
                 for a in ARMS for n in range(1, N_STEPS + 1))
    require(closed < 1e-19,
            f"X1 + By is exactly the injected 6.0e-5 at every level of every arm "
            f"(worst {closed:.2e} m3/s): the shape account moves where the rate goes, "
            "never how much of it there is")

    # ===== 8. 5.40's two brackets survive, with the chord's target moved
    idw, chordMax, attribMax = 0.0, {}, {}
    for (a, n, i) in ALL:
        l = E(a, n, i)
        c1 = value(l, "dpNetwork") - (value(l, "dpFriction") - value(l, "rowBeta"))
        c2 = value(l, "M") * value(l, "dq") / DT - value(l, "dpInertia")
        idw = max(idw, abs(value(l, "rowError") * value(l, "dpBlock") - (c1 + c2)))
        chordMax[(a, i)] = max(chordMax.get((a, i), 0.0), abs(c1))
        attribMax[(a, i)] = max(attribMax.get((a, i), 0.0), abs(c2))
    require(idw < 1e-11,
            f"the identity rowError*dpBlock = chord bracket + attribution bracket "
            f"closes to {idw:.1e} Pa over all 64 states with the chord's target now "
            "reduced by the shape term: the structure 5.40 found is unchanged by who "
            "pays for the shape")
    require(max(chordMax.values()) < 1.0 and attribMax[("c_off", 1)] > 1900.0
            and attribMax[("c_off", 0)] > 750.0,
            f"and the two brackets are still {min(attribMax[k]/chordMax[k] for k in attribMax):.0f}x "
            f"to {max(attribMax[k]/chordMax[k] for k in attribMax):.0f}x apart: chord at "
            f"most {max(chordMax.values()):.3f} Pa, attribution {attribMax[('c_off',0)]:.2f} Pa "
            f"on X1 and {attribMax[('c_off',1)]:.2f} Pa on X2 -- moving the shape account "
            "did not turn it into a third bill")

    # ===== 9. the account is the region's, not the row's
    for i in (0, 1):
        coef = lambda a, n: abs(value(E(a, n, i), "dpConv")) / value(E(a, n, i), "q") ** 2
        sp = [max(coef(a, n) for a in ARMS) / min(coef(a, n) for a in ARMS) - 1.0
              for n in range(1, N_STEPS + 1)]
        require(max(sp) < 2e-4,
                f"on {NAME[i]} the coefficient dpConv/q^2 agrees across all four arms "
                f"to {max(sp):.1e} relative at every level (it is worth "
                f"{min(coef(a, n) for a in ARMS for n in (1, 8)):.4e} to "
                f"{max(coef(a, n) for a in ARMS for n in (1, 8)):.4e} Pa.s2/m6): who "
                "pays for the shape does not change what the shape is worth")
    gradeX2 = max(abs(value(E("c_pin", n, 1), "rowError")
                      - value(E("c_off", n, 1), "rowError")) for n in range(1, 9))
    gradeX1 = max(abs(value(E("c_free", n, 0), "rowError")
                      - value(E("c_off", n, 0), "rowError")) for n in range(1, 9))
    require(gradeX2 < 2e-4 and gradeX1 < 1e-2,
            f"the grade the two sides are given on does not care about the split "
            f"either: rowError moves by {gradeX2:.1e} on the pinned edge and "
            f"{gradeX1:.1e} on the free one")

    # ===== 10. what it costs
    passes = {a: {n: len(v) for n, v in by[a]["pass"].items()} for a in ARMS}
    tot = {a: sum(p.values()) for a, p in passes.items()}
    cap = max(max(p.values()) for p in passes.values())
    require(all(max(p.values()) <= 7 for p in passes.values()),
            f"no arm came near the 12-pass cap: the worst level anywhere needed "
            f"{cap} passes ({list(passes['c_off'].values())} for c_off against "
            f"{list(passes['c_free'].values())} for c_free)")
    require(tot["c_free"] <= tot["c_off"] + 1 and tot["c_pin"] == tot["c_off"],
            f"the whole interface costs {tot['c_off']} passes reported and "
            f"{tot['c_free']} charged on the free edge ({tot['c_both']} on both), with "
            f"{by['c_off']['solves']} against {by['c_free']['solves']} block solves: "
            "one extra pass at level 2 is paid back at levels 6..8, and charging the "
            "pinned edge costs nothing at all")

    print(f"\n{len(fails)} failure(s)")
    for f in fails:
        print("  " + f)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
