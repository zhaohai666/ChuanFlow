#!/usr/bin/env python3
"""5.44: the two accounts the 1D row owes are paid on the same edge at the same
time -- and neither lever touches the other's bill.

    Src -> A -> { X1 (3D kOmegaSST duct, parallel, rate free)
                  || By (1D pipe) } -> X2 (3D RANS duct, series, q pinned)
         -> X2 -> B -> Snk

Six arms, a full 3x2 factorial on X2 alone. X1 keeps the drawing's rho*L/A and
carries no shape charge in any arm, so every difference in the run comes from one
of two declarations on one edge:

    inertiaFrom.X2     -- whose rho*L/A^2 the row installs as M: `network` (the
                          drawing), `partner` (the region the block solves in),
                          `slab` (the fluid the block's drop actually crosses)
    convectiveToRow.X2 -- whether the region's changing profile is paid by the
                          block's chord (`report`) or written into the row
                          (`charge`)

                     n_rep   n_chr   g_rep   g_chr   s_rep   s_chr
    X2 rung          network network partner partner slab    slab
    X2 shape account report  charge  report  charge  report  charge

`n_rep` and `n_chr` are shipped cells -- they are 5.43's `c_off` and `c_pin` -- so
the crossing is graded against published numbers and not only against itself.

What the kernel promises, and what this script checks:

  1. both tables are read back per arm and per edge from the log, and they are the
     only things that differ: same seed, same shape, same scheme, same dt.
  2. each rung installs the closed form it claims: rho*L/A for the drawing,
     rho*V_region/A_region^2 for the region, rho*V_probe/A_region^2 for the slab.
  3. the attribution staircase falls rung by rung with the shape axis held still,
     and its top step is 5.41's published 1903.4177 Pa on the pinned edge at the
     level where that edge's rate first moves.
  4. the charge does not bend the staircase: within a rung the installed inertance
     is bit-identical with and without it at every level, and the block's own
     dpInertia agrees to 2.3e-3 Pa out of 643.
  5. the rung does not bend the shape account either: dpConv on X2 is bit-identical
     across all six arms at levels 1..7 and agrees to 2.8e-6 relative at level 8,
     while the inertance installed beside it differs by 3.88x.
  6. nothing is charged twice: rowBeta is dpConv on the 24 charged states, exactly 0
     on the 72 reported ones (every one of which is an X1 state or a `_rep` arm's
     X2), and rowDrop closes as dpNetwork + M*dq/dt + rowBeta with zero slack.
  7. 5.40's decomposition survives both levers at once, and the two brackets it
     splits are individually addressable: the chord bracket stays inside the
     interface tolerance in every arm no matter which rung is installed, the bill
     bracket falls 1903.42 -> 1606.94 -> 13.19 Pa at the declaration, and the sum
     falls 144x when both levers are used -- while X1, declared on neither, keeps
     its 754.24 Pa bill all six arms through.
  8. neither declaration reaches the neighbour: X1 is bit-identical to `n_rep` in
     all 17 fields at all 8 levels in the other five arms.
  9. `n_rep` reproduces 5.41's `p_net` state for state and pass for pass, and all
     six arms need the same eight pass counts and the same 70 block solves, so the
     two tables are declarations the interface reads, not extra work.

usage: check_transient_coupled_multi_turbulent_shape_and_inertia_split.py <solver> <case>
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
NAME = ("X1", "X2")
ARMS = ("n_rep", "n_chr", "g_rep", "g_chr", "s_rep", "s_chr")
RUNGS = ("network", "partner", "slab")
PAIRS = (("n_rep", "n_chr"), ("g_rep", "g_chr"), ("s_rep", "s_chr"))
RUNG = {"n_rep": "network", "n_chr": "network", "g_rep": "partner",
        "g_chr": "partner", "s_rep": "slab", "s_chr": "slab"}
CHARGE = {"n_rep": "report", "n_chr": "charge", "g_rep": "report",
          "g_chr": "charge", "s_rep": "report", "s_chr": "charge"}
SEED = {"X1": 1.4e7, "X2": 1.36e7}
FIELDS = ("q", "dpBlock", "dpNetwork", "R", "residual", "M", "dpFriction",
          "dpInertia", "qPrev", "dq", "rowDrop", "rowError", "betaIn", "betaOut",
          "pRisePa", "dpConv", "rowBeta")
# The drawing's two coupled pipes and the two blocks that solve inside them: X2 is
# a 0.24 m pipe drawn at diameter 0.006, solved in as a 0.24 x 0.004 x 0.008 m
# region whose probe is the last 0.07 m; X1 likewise at 0.2 x 0.005 x 0.005 with a
# 0.07 m probe. The three rungs are then three closed forms, per edge.
GEOM = {"X1": {"L": 0.2, "d": 0.005, "A": 0.005 * 0.005, "Lp": 0.07},
        "X2": {"L": 0.24, "d": 0.006, "A": 0.004 * 0.008, "Lp": 0.07}}


def M_of(edge, source):
    g = GEOM[edge]
    if source == "network":
        return RHO * g["L"] / (math.pi * g["d"] ** 2 / 4.0)
    if source == "partner":
        return RHO * g["L"] / g["A"]
    return RHO * g["Lp"] / g["A"]


M_RUNG = {"X2": {w: M_of("X2", w) for w in RUNGS}}
DQDT_L1 = Q_IN / DT
# 5.41's `p_net` level-1 and level-8 states, which are this case's `n_rep`.
PINET_L1 = {"X1": {"dpBlock": 851.08258956984685, "dpFriction": 579.40211914132692,
                   "dpNetwork": 579.40721946887868, "dpConv": -55.239173285473818,
                   "rowError": 0.88621650000126295},
            "X2": {"dpBlock": 2868.6479781972207, "dpFriction": 2225.5865729289981,
                   "dpNetwork": 2225.5627280310264, "dpConv": -277.29484678803749,
                   "rowError": 0.66351600258052723}}
PINET_L8 = {"X1": {"q": 3.1118825189579614e-05, "dpBlock": 177.93872006380388,
                   "dpNetwork": 164.74362239216259, "dpFriction": 164.77917883517469,
                   "R": 5294018.0546188587, "rowDrop": 212.60041875660994,
                   "dpConv": -148.97387627505432, "rowError": 0.19479570652400632},
            "X2": {"q": 6.0000000000000002e-05, "dpBlock": 720.44658376591372,
                   "dpNetwork": 720.4249438241302, "dpFriction": 720.44848007364146,
                   "R": 12007082.397068836, "rowDrop": 720.4249438241302,
                   "dpConv": -283.2722617002043, "rowError": -3.0036844189617738e-05}}
PINET_PASSES = (3, 6, 5, 5, 4, 4, 4, 4)
# 5.41's published parabolic bill on the pinned edge at level 1, rebuilt from the
# four shipped fields it decomposes into rather than quoted from a memory of it.
BILL_TOP_L1 = (PINET_L1["X2"]["rowError"] * PINET_L1["X2"]["dpBlock"]
               + PINET_L1["X2"]["dpFriction"] - PINET_L1["X2"]["dpNetwork"])


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
        raise SystemExit("transient shape-and-inertia-split solver failed")

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
    rung = lambda a, e: by[a]["inertia"][e]
    ALL = [(a, n, i) for a in ARMS for n in range(1, N_STEPS + 1) for i in (0, 1)]
    bill = lambda a, n, i: (value(E(a, n, i), "M") * value(E(a, n, i), "dq") / DT
                            - value(E(a, n, i), "dpInertia"))
    chord = lambda a, n, i: (value(E(a, n, i), "dpFriction")
                             - value(E(a, n, i), "rowBeta"))
    gap = lambda a, n, i: value(E(a, n, i), "dpNetwork") - chord(a, n, i)

    # ===== 1. two tables, one edge, everything else held
    header = [ln for ln in lines if ln.startswith("transient network+blocks")]
    require(len(header) == 1 and "2 coupled edges, 6 starting sets" in header[0]
            and "per-edge on some arms" in header[0]
            and "shape charged on some edges" in header[0],
            "one network, two coupled edges, six starting sets, and the header "
            "announces that both who-pays tables are per-edge somewhere")
    require([c["arm"] for c in arms] == list(ARMS),
            f"the six arms ran in the order the case declares them (got "
            f"{[c['arm'] for c in arms]})")
    require(all(c["seed"] == SEED for c in arms),
            f"every arm seeds both edges with the same pair ({SEED['X1']:.3e}, "
            f"{SEED['X2']:.3e}), so neither table can come in through the "
            "starting coefficients")
    require(all(c["shape"] == {"X1": "parabolic", "X2": "parabolic"} for c in arms),
            "every arm is parabolic on every edge: the shape stays put, which is "
            "what keeps dpConv nonzero and negative on all 96 states")
    require([tuple(rung(a, e) for e in NAME) for a in ARMS]
            == [("network", RUNG[a]) for a in ARMS],
            "the attribution table is read back per arm and per edge from the log: "
            f"X1 keeps the drawing everywhere while X2 walks {RUNGS[0]} -> "
            f"{RUNGS[1]} -> {RUNGS[2]}")
    require([tuple(by[a]["conv"][n] for n in NAME) for a in ARMS]
            == [("report", CHARGE[a]) for a in ARMS],
            "and the charge table is read back the same way, moved only on X2: "
            f"{[CHARGE[a] for a in ARMS]}")
    require(all(int(value(E(a, n, i), "unsteady")) == 1 for (a, n, i) in ALL),
            "all 96 (arm x level x edge) states are stepping: the account the two "
            "levers divide is the momentum budget of a block that carries its own "
            "du/dt, so a quasi-steady partner would leave this case asking nothing "
            "of either table")

    # ===== 2. each rung installs the closed form it claims
    bad_m = [(a, NAME[i], rung(a, NAME[i]), value(E(a, n, i), "M"),
              M_of(NAME[i], rung(a, NAME[i])))
             for (a, n, i) in ALL
             if abs(value(E(a, n, i), "M") - M_of(NAME[i], rung(a, NAME[i])))
             > 1e-9 * M_of(NAME[i], rung(a, NAME[i]))]
    require(not bad_m,
            "every arm installs the inertance its rung defines, at every level and "
            "on both edges: rho*L/A for the drawing, rho*L/A_region for the region "
            "(= rho*V/A_region^2), rho*L_probe/A_region for the slab -- on X2 that "
            "is " + ", ".join(f"{w} {M_RUNG['X2'][w]:.6e}" for w in RUNGS)
            + " kg/m4, and X1 keeps the drawing's "
            f"{M_of('X1', 'network'):.6e} in all six"
            + (f"; first mismatch {bad_m[0]}" if bad_m else ""))
    require(abs(bill("n_rep", 1, 1) - BILL_TOP_L1) < 1e-9,
            f"the top step of the staircase is a shipped number, not a fresh one: "
            f"{bill('n_rep', 1, 1):.4f} Pa on the pinned edge at level 1, which is "
            f"{BILL_TOP_L1:.4f} Pa read out of 5.41's published rowError, dpBlock, "
            "dpNetwork and dpFriction for this rung at this rate change")

    # ===== 3. the staircase falls, and the charge does not bend it
    stairs = {a: bill(a, 1, 1) for a in ARMS}
    require(stairs["s_rep"] < stairs["g_rep"] < stairs["n_rep"],
            "with the shape axis held at `report` the attribution axis still does "
            f"the one job it was bought for: {stairs['n_rep']:.4f} Pa at the "
            f"drawing -> {stairs['g_rep']:.4f} at the region -> "
            f"{stairs['s_rep']:.4f} at the slab, a fall of "
            f"{stairs['n_rep'] / stairs['s_rep']:.1f}x")
    require(all(abs(stairs[a] - stairs[b]) < 1e-2 for a, b in PAIRS),
            "and the same three heights with the shape account written into the "
            f"row: {stairs['n_chr']:.4f} / {stairs['g_chr']:.4f} / "
            f"{stairs['s_chr']:.4f} Pa, each within "
            f"{max(abs(stairs[a] - stairs[b]) for a, b in PAIRS):.1e} Pa of its "
            "reported twin -- so the two levers have different addresses")
    d_m = {p: max(abs(value(E(p[0], n, 1), "M") - value(E(p[1], n, 1), "M"))
                  for n in range(1, N_STEPS + 1)) for p in PAIRS}
    require(all(v == 0.0 for v in d_m.values()),
            "within a rung the installed inertance is bit-identical with and "
            f"without the charge at every level ({d_m}) -- the charge moves a term "
            "in the row, not the mass the row is made of")
    d_i = {p: max(abs(value(E(p[0], n, 1), "dpInertia")
                      - value(E(p[1], n, 1), "dpInertia"))
                   for n in range(1, N_STEPS + 1)) for p in PAIRS}
    worst_i = max(d_i.values())
    require(worst_i < 1e-2,
            "and so is the block's own report of what its region accelerated: "
            f"the rung pairs agree to {d_i[('g_rep', 'g_chr')]:.1e} and "
            f"{d_i[('s_rep', 's_chr')]:.1e} Pa exactly and the drawing's pair to "
            f"{worst_i:.3e} Pa out of {value(E('n_rep', 1, 1), 'dpInertia'):.1f} "
            f"({worst_i / value(E('n_rep', 1, 1), 'dpInertia'):.1e} relative) -- "
            "the charge is invisible to the side that measured it")

    # ===== 4. the rung does not bend the shape account either
    q_pin = [(a, n, value(E(a, n, 1), "q")) for a in ARMS
             for n in range(1, N_STEPS + 1) if value(E(a, n, 1), "q") != Q_IN]
    require(not q_pin,
            "X2's rate is the injected 6.0e-5 to the last bit in all six arms at "
            "all eight levels, so a rung and a charge there cannot be confounded "
            "with a different trajectory"
            + (f"; first exception {q_pin[0]}" if q_pin else ""))
    bit_conv = [n for n in range(1, 8)
                if len({value(E(a, n, 1), "dpConv") for a in ARMS}) == 1]
    spread = {n: max(abs(value(E(a, n, 1), "dpConv")) for a in ARMS)
              / min(abs(value(E(a, n, 1), "dpConv")) for a in ARMS) - 1.0
              for n in range(1, N_STEPS + 1)}
    require(bit_conv == list(range(1, 8)) and spread[8] < 1e-5,
            f"the shape account is the same number whoever installs the inertance: "
            f"dpConv on X2 is bit-identical across all six arms at levels "
            f"{bit_conv[0]}..{bit_conv[-1]} and agrees to {spread[8]:.1e} relative "
            f"at level 8, while the inertance installed beside it differs by "
            f"{value(E('n_rep', 1, 1), 'M') / value(E('s_rep', 1, 1), 'M'):.4f}x "
            "-- the account is a property of the region, and the run measured that "
            "rather than inferring it")
    coef = {n: max(abs(value(E(a, n, 1), "dpConv")) / value(E(a, n, 1), "q") ** 2
                   for a in ARMS)
            / min(abs(value(E(a, n, 1), "dpConv")) / value(E(a, n, 1), "q") ** 2
                  for a in ARMS) - 1.0 for n in range(1, N_STEPS + 1)}
    require(max(coef.values()) < 1e-5,
            f"and so is its coefficient per q^2, worst {max(coef.values()):.1e} "
            f"(it is worth {abs(value(E('n_rep', 1, 1), 'dpConv')) / Q_IN ** 2:.4e} "
            f"to {abs(value(E('n_rep', 8, 1), 'dpConv')) / Q_IN ** 2:.4e} "
            "Pa.s2/m6 up and down the ramp) -- 5.43 found the charge does not "
            "change what the shape is worth; this finds the rung does not either")

    # ===== 5. one bill, two payers, no third
    n_ch = len([1 for x in ALL if charged(x[0], x[2])])
    beta_bad = [(a, n, NAME[i], value(E(a, n, i), "rowBeta"))
                for (a, n, i) in ALL
                if abs(value(E(a, n, i), "rowBeta")
                       - (value(E(a, n, i), "dpConv") if charged(a, i) else 0.0))
                > 1e-9]
    require(not beta_bad,
            f"rowBeta is the whole account on the {n_ch} charged states and exactly "
            f"0 on the {len(ALL) - n_ch} reported ones, and every reported state is "
            "either an X1 state or a `_rep` arm's X2 -- X1 was never asked to pay"
            + (f"; first mismatch {beta_bad[0]}" if beta_bad else ""))
    row_bad = max(abs(value(E(a, n, i), "rowDrop")
                      - (value(E(a, n, i), "dpNetwork")
                         + value(E(a, n, i), "M") * value(E(a, n, i), "dq") / DT
                         + value(E(a, n, i), "rowBeta"))) for (a, n, i) in ALL)
    require(row_bad == 0.0,
            f"the printed row closes as dpNetwork + M*dq/dt + rowBeta with no slack "
            f"at all in any of the {len(ALL)} states (worst {row_bad:.1e} Pa): the "
            "two tables split one statement, they do not add a second bill")

    # ===== 6. 5.40's identity, with both brackets under declaration
    worst_id = max(abs(value(E(a, n, i), "rowError") * value(E(a, n, i), "dpBlock")
                       - (gap(a, n, i) + bill(a, n, i))) for (a, n, i) in ALL)
    require(worst_id < 1e-9,
            f"the decomposition 5.40 found survives both levers: rowError*dpBlock "
            "== (dpNetwork - (dpFriction - rowBeta)) + (M*dq/dt - dpInertia) in all "
            f"{len(ALL)} states, worst {worst_id:.3e} Pa -- the charge changes what "
            "the chord bracket is measured against, the rung changes what the bill "
            "bracket holds, and neither changes the sum being split")
    rel_gap = max(abs(gap(a, n, i)) / abs(chord(a, n, i)) for (a, n, i) in ALL)
    abs_gap = {a: max(abs(gap(a, n, i)) for n in range(1, N_STEPS + 1)
                      for i in (0, 1)) for a in ARMS}
    require(rel_gap <= TOL,
            f"the chord bracket stays inside the interface's own tolerance in every "
            f"arm of the factorial (worst {rel_gap:.4e} relative against the "
            f"{TOL:.0e} the negotiation is held to; {min(abs_gap.values()):.4f} to "
            f"{max(abs_gap.values()):.4f} Pa absolute) -- installing a different "
            "rung buys the row a smaller inertial term, not a looser chord")
    worst_bill = {(a, i): max(abs(bill(a, n, i)) for n in range(1, N_STEPS + 1))
                  for a in ARMS for i in (0, 1)}
    bill_edge = {}
    for a in ARMS:
        for i in (0, 1):
            bill_edge[(i, rung(a, NAME[i]))] = max(
                bill_edge.get((i, rung(a, NAME[i])), 0.0), worst_bill[(a, i)])
    require(abs(bill_edge[(0, "network")] - 754.2383) < 5e-3
            and abs(bill_edge[(1, "network")] - 1903.4177) < 5e-3,
            "and the bill bracket is a property of where the declaration was made, "
            f"not of the run: X1, declared on neither table, keeps its "
            f"{bill_edge[(0, 'network')]:.4f} Pa in all six arms, while X2's own "
            f"worst falls {bill_edge[(1, 'network')]:.4f} -> "
            f"{bill_edge[(1, 'partner')]:.4f} -> {bill_edge[(1, 'slab')]:.4f} Pa "
            "rung by rung")
    both = {a: abs(gap(a, 1, 1)) + abs(bill(a, 1, 1)) for a in ARMS}
    re_l1 = {a: value(E(a, 1, 1), "rowError") for a in ARMS}
    require(both["s_chr"] * 100.0 < both["n_rep"]
            and re_l1["s_rep"] < re_l1["g_rep"] < re_l1["n_rep"],
            f"the two brackets sum to {both['n_rep']:.2f} Pa at the drawing's rung "
            f"reported and {both['s_chr']:.2f} Pa at the slab rung charged, "
            f"{both['n_rep'] / both['s_chr']:.1f}x apart, and the row's own error "
            f"follows it: {re_l1['n_rep']:.4f} -> {re_l1['g_rep']:.4f} -> "
            f"{re_l1['s_rep']:.4f} -- the first state in this chain whose row error "
            "is bounded from both sides at once")
    require(all(abs(re_l1[a] - re_l1[b]) < 1e-5 for a, b in PAIRS),
            "and the charge leaves that error where it found it: per rung the pair "
            f"agrees to {max(abs(re_l1[a] - re_l1[b]) for a, b in PAIRS):.1e} "
            f"({re_l1['n_rep']:.6f} / {re_l1['g_rep']:.6f} / {re_l1['s_rep']:.6f} "
            "reported against "
            f"{re_l1['n_chr']:.6f} / {re_l1['g_chr']:.6f} / {re_l1['s_chr']:.6f} "
            "charged) -- paying the shape half does not pay the inertial half, "
            "which is exactly what 5.43's boundary left open")

    # ===== 7. neither declaration reaches the neighbour
    for a in ARMS[1:]:
        diff = {f: max(abs(value(E(a, n, 0), f) - value(E("n_rep", n, 0), f))
                       for n in range(1, N_STEPS + 1)) for f in FIELDS}
        require(all(v == 0.0 for v in diff.values()),
                f"{a}'s X1 is bit-identical to n_rep's in all {len(FIELDS)} fields "
                f"at all {N_STEPS} levels: a rung on a neighbour and a charge on a "
                "neighbour are both invisible here, the way 5.43 found for the "
                "charge alone")

    # ===== 8. the controls, and what the run cost
    bad1 = [(NAME[i], f, value(E("n_rep", 1, i), f), PINET_L1[NAME[i]][f])
            for i in (0, 1) for f in PINET_L1[NAME[i]]
            if abs(value(E("n_rep", 1, i), f) - PINET_L1[NAME[i]][f]) > 1e-12]
    require(not bad1,
            "n_rep reproduces 5.41's `p_net` level-1 states field for field -- the "
            "shipped numbers, not a copy of this run"
            + (f"; first mismatch {bad1[0]}" if bad1 else ""))
    bad8 = [(NAME[i], f, value(E("n_rep", N_STEPS, i), f), PINET_L8[NAME[i]][f])
            for i in (0, 1) for f in PINET_L8[NAME[i]]
            if abs(value(E("n_rep", N_STEPS, i), f) - PINET_L8[NAME[i]][f]) > 1e-12]
    require(not bad8,
            "and its level-8 states too, including the coefficient the interface "
            "converged to and the row's total"
            + (f"; first mismatch {bad8[0]}" if bad8 else ""))
    seq = {a: [len(by[a]["pass"][n]) for n in range(1, N_STEPS + 1)] for a in ARMS}
    require(tuple(seq["n_rep"]) == PINET_PASSES,
            f"and the interface worked just as hard level by level (got "
            f"{seq['n_rep']} against the shipped {list(PINET_PASSES)}), which is the "
            "part no printed field can show")
    require(all(v == seq["n_rep"] for v in seq.values()),
            f"and all six arms of the factorial work exactly as hard as each other "
            f"({seq['n_rep']} everywhere, worst level {max(max(v) for v in seq.values())} "
            f"against a cap of 12): neither table buys extra negotiation, because "
            "neither table is an iteration -- it is a coefficient the row is written "
            "with")
    solves = {a: by[a]["solves"] for a in ARMS}
    require(len(set(solves.values())) == 1,
            f"block solves are the same in every arm ({list(solves.values())} -- "
            f"{max(solves.values()) * len(ARMS)} in total for the factorial): the "
            "cost of the crossing is the crossing itself, not the physics")
    closed = max(abs(value(by[a]["pipe"][n], "X1")
                     + value(by[a]["pipe"][n], "By") - Q_IN)
                 for a in ARMS for n in range(1, N_STEPS + 1))
    require(closed < 1e-19,
            f"and the branch still closes: X1 + By is the injected 6.0e-5 at every "
            f"level of every arm (worst {closed:.2e} m3/s)")

    print(f"\n{len(fails)} failure(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
