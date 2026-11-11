#!/usr/bin/env python3
"""B5 (5.67): the three per-edge tables AND the per-edge profile land in ONE
case -- the full coupling grammar in a single file, checked end to end.

    Src -> A -> { X1 (3D kOmegaSST duct, parallel, rate free)
                  || By (1D pipe) } -> X2 (3D RANS duct, series, q pinned)
         -> B -> Snk

Every cell of the (ddtFrom, inertiaFrom, convectiveToRow, profile) board had
been visited somewhere -- 5.40/5.41 ran ddtFrom + inertiaFrom + profile,
5.43/5.44 ran inertiaFrom + convectiveToRow + profile -- but no single case
ever carried all three tables at once. Five arms on 5.41's topology, meshes,
probes, dt and seeds, so published twins exist:

              ddtFrom        inertiaFrom        convectiveToRow
    base      Euler, Euler   network, network   report, report    (= 5.41 p_net)
    mix       Euler, steady  network, slab      report, report    (= 5.41 p_mix)
    mix_charge  "            "                  report, CHARGE    <- THE new cell
    pin       steady, Euler  network, slab      report, report    (= 5.41 p_pin)
    both_chr  Euler, Euler   network, slab      report, CHARGE    (= 5.44 s_chr)

What the kernel promises, and what this script checks:

  1. all three tables and the shape are read back per arm and per edge from
     the log -- scheme from the `unsteady` flag of every state, inertia source
     and shape from the arm lines, the charge account from rowBeta -- and the
     arms run in the order the case declares.
  2. each rung installs the closed form it claims, at every level: rho*L/A for
     the drawing, rho*L_probe/A_region for the slab.
  3. the quasi edges ddtFrom buys are real: unsteady 0, rowError 0 by
     definition, dpInertia 0 and dpFriction == dpBlock (the whole drop, nothing
     graded) on all 24 quasi states.
  4. the new cell works: mix_charge's X2 is quasi AND slab AND charged at once
     -- its row pays M*dq/dt with the slab's number and carries rowBeta ==
     dpConv != 0 at all eight levels; the mix twin, identical but for the
     charge word, has rowBeta 0 and M/dpBlock/dpConv/dpFriction/dpInertia
     bit-identical, so the charge moved exactly one term of one row.
  5. the charge does not reach the neighbour: mix's X1 is bit-identical to
     mix_charge's X1 in all 17 fields at all 8 levels.
  6. every row closes as dpNetwork + M*dq/dt + rowBeta with zero slack in all
     80 states, and rowBeta is the whole account on the 16 charged states and
     exactly 0 on the 24 reported ones.
  7. the chord bracket stays inside the interface tolerance in every arm; the
     base arm reproduces 5.41's published p_net level-1 and level-8 states and
     its pass sequence; X2's rate is the injected 6.0e-5 everywhere; the branch
     still closes; the charge buys no extra negotiation.

usage: check_transient_coupled_multi_turbulent_triple_table.py <solver> <case>
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
ARMS = ("base", "mix", "mix_charge", "pin", "both_chr")
SEED = {"X1": 1.4e7, "X2": 1.36e7}
FIELDS = ("q", "dpBlock", "dpNetwork", "R", "residual", "M", "dpFriction",
          "dpInertia", "qPrev", "dq", "rowDrop", "rowError", "betaIn", "betaOut",
          "pRisePa", "dpConv", "rowBeta")
GEOM = {"X1": {"L": 0.2, "d": 0.005, "A": 0.005 * 0.005, "Lp": 0.07},
        "X2": {"L": 0.24, "d": 0.006, "A": 0.004 * 0.008, "Lp": 0.07}}
# The quasi cells ddtFrom buys, and the charge cells convectiveToRow buys.
QUASI = {("mix", "X2"), ("mix_charge", "X2"), ("pin", "X1")}
CHARGED = {("mix_charge", "X2"), ("both_chr", "X2")}
# 5.41's published `p_net` states, which are this case's `base` -- shipped
# numbers, not a copy of this run.
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


def M_of(edge, source):
    g = GEOM[edge]
    if source == "network":
        return RHO * g["L"] / (math.pi * g["d"] ** 2 / 4.0)
    return RHO * g["Lp"] / g["A"]


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
        raise SystemExit("transient triple-table solver failed")

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
            cur = {"arm": t[2].rstrip(":"), "edge": {}, "pass": {}, "seed": {},
                   "solves": 0}
            arms.append(cur)
        elif cur is None:
            continue
        if ln.startswith("SIMPLE iter 1:"):
            cur["solves"] += 1
        elif ln.startswith("transient step ") and " pass " in ln:
            cur["pass"].setdefault(int(t[2]), []).append(ln)
        elif ln.startswith("transient step ") and " edge " in ln:
            cur["edge"].setdefault(int(t[2]), {})[t[6]] = ln
        elif ln.startswith("transient start ") and " resistance " in ln:
            cur["seed"][word(ln, "edge")] = value(ln, "resistance")

    by = {c["arm"]: c for c in arms}

    def E(a, n, i):
        return by[a]["edge"][n][NAME[i]]

    ALL = [(a, n, i) for a in ARMS for n in range(1, N_STEPS + 1) for i in (0, 1)]

    # the three tables, read back from their announcements
    inertia, conv, shape = {}, {}, {}
    for ln in lines:
        m = re.match(r"transient arm (\w+) edge (\w+) inertia (\w+)", ln)
        if m: inertia[(m.group(1), m.group(2))] = m.group(3)
        m = re.match(r"transient arm (\w+) edge (\w+) convective (\w+)", ln)
        if m: conv[(m.group(1), m.group(2))] = m.group(3)
        m = re.match(r"transient arm (\w+) edge (\w+) profile (\w+)", ln)
        if m: shape[(m.group(1), m.group(2))] = m.group(3)

    # ===== 1. one case, three tables, read back per arm and per edge
    header = [ln for ln in lines if ln.startswith("transient network+blocks")]
    require(len(header) == 1
            and f"2 coupled edges, {len(ARMS)} starting sets" in header[0]
            and "per-edge on some arms" in header[0]
            and "shape charged on some edges" in header[0],
            f"one network, two coupled edges, {len(ARMS)} starting sets, and the "
            "header announces both a per-edge inertia table and a charge table")
    require([c["arm"] for c in arms] == list(ARMS),
            f"the five arms ran in the order the case declares them (got "
            f"{[c['arm'] for c in arms]})")
    require(all(c["seed"] == SEED for c in arms),
            f"every arm seeds both edges with the same pair ({SEED['X1']:.3e}, "
            f"{SEED['X2']:.3e}), so no table can come in through the starting "
            "coefficients")
    require(all(shape[(a, e)] == "parabolic" for a in ARMS for e in NAME),
            "the shape axis is parabolic on every edge of every arm: the profile "
            "is a per-edge field of the same case, held still while the three "
            "tables move")
    require({k: v for k, v in inertia.items() if k[0] in ARMS}
            == {(a, "X1"): "network" for a in ARMS}
            | {("base", "X2"): "network"} | {("mix", "X2"): "slab"}
            | {("mix_charge", "X2"): "slab"} | {("pin", "X2"): "slab"}
            | {("both_chr", "X2"): "slab"},
            "the inertia table is read back per arm and per edge: X1 keeps the "
            "drawing in all five arms while X2 walks network -> slab x4")
    require({k: v for k, v in conv.items() if k[0] in ARMS}
            == {(a, "X1"): "report" for a in ARMS}
            | {("mix_charge", "X2"): "charge"} | {("both_chr", "X2"): "charge"}
            | {("base", "X2"): "report"} | {("mix", "X2"): "report"}
            | {("pin", "X2"): "report"},
            "and the charge table with it: report everywhere except the two "
            "arms that say charge on X2 -- three tables, one grammar, one file")
    got_quasi = {(a, NAME[i]) for (a, n, i) in ALL
                 if value(E(a, n, i), "unsteady") == 0}
    require(got_quasi == QUASI,
            "the ddtFrom table is read back from the states themselves: exactly "
            f"{sorted(QUASI)} are quasi (unsteady 0) in all 8 levels and every "
            "other (arm, edge) steps -- the third table, read back the same way")

    # ===== 2. each rung installs the closed form it claims
    bad_m = [(a, NAME[i], inertia[(a, NAME[i])])
             for (a, n, i) in ALL
             if abs(value(E(a, n, i), "M") - M_of(NAME[i], inertia[(a, NAME[i])]))
             > 1e-9 * M_of(NAME[i], inertia[(a, NAME[i])])]
    require(not bad_m,
            f"every arm installs the inertance its rung defines at every level "
            f"(all {len(ALL)} states): rho*L/A for the drawing "
            f"({M_of('X1', 'network'):.6e} / {M_of('X2', 'network'):.6e} kg/m4), "
            f"rho*L_probe/A_region for the slab ({M_of('X1', 'slab'):.6e} / "
            f"{M_of('X2', 'slab'):.6e} kg/m4)"
            + (f"; first mismatch {bad_m[0]}" if bad_m else ""))

    # ===== 3. the quasi cells are real quasi cells
    bad_q = [(a, NAME[i]) for (a, n, i) in ALL
             if (a, NAME[i]) in QUASI
             and (value(E(a, n, i), "rowError") != 0.0
                  or value(E(a, n, i), "dpInertia") != 0.0
                  or value(E(a, n, i), "dpFriction") != value(E(a, n, i), "dpBlock"))]
    require(not bad_q,
            "every quasi state is the real thing: rowError 0 by definition, "
            "dpInertia 0, and dpFriction == dpBlock -- the partner's whole drop "
            "with nothing graded on it, on all 24 quasi states"
            + (f"; first mismatch {bad_q[0]}" if bad_q else ""))

    # ===== 4. the new cell: three tables stacked on one edge
    stack_ok = all(
        value(E("mix_charge", n, 1), "unsteady") == 0
        and abs(value(E("mix_charge", n, 1), "M") - M_of("X2", "slab"))
        < 1e-9 * M_of("X2", "slab")
        and value(E("mix_charge", n, 1), "rowBeta")
        == value(E("mix_charge", n, 1), "dpConv")
        and value(E("mix_charge", n, 1), "dpConv") != 0.0
        for n in range(1, N_STEPS + 1))
    require(stack_ok,
            "THE new cell: mix_charge's X2 is quasi AND slab AND charged at once "
            "-- unsteady 0, M = the slab's "
            f"{M_of('X2', 'slab'):.6e}, and rowBeta == dpConv != 0 at all "
            f"{N_STEPS} levels (smallest {min(abs(value(E('mix_charge', n, 1), 'dpConv')) for n in range(1, N_STEPS + 1)):.4f} "
            "Pa) -- three tables on one edge of one arm of one case")
    twin = all(value(E("mix", n, 1), f) == value(E("mix_charge", n, 1), f)
               for n in range(1, N_STEPS + 1)
               for f in ("dpBlock", "dpConv", "dpFriction", "dpInertia", "M"))
    twin &= value(E("mix", 1, 1), "rowBeta") == 0.0
    require(twin,
            "and the mix twin, identical but for one word, has rowBeta 0 and "
            "dpBlock/dpConv/dpFriction/dpInertia/M bit-identical at every level: "
            "the charge moved exactly one term of one row and nothing else")

    # ===== 5. the charge does not reach the neighbour
    diff = {f: max(abs(value(E("mix_charge", n, 0), f) - value(E("mix", n, 0), f))
                   for n in range(1, N_STEPS + 1)) for f in FIELDS}
    require(all(v == 0.0 for v in diff.values()),
            f"mix_charge's X1 is bit-identical to mix's X1 in all {len(FIELDS)} "
            "fields at all 8 levels: the charge on the neighbour changes what "
            "X2's own row is written with and nothing about X1")

    # ===== 6. one row, one bill, no slack
    worst_row = max(abs(value(E(a, n, i), "rowDrop")
                        - (value(E(a, n, i), "dpNetwork")
                           + value(E(a, n, i), "M") * value(E(a, n, i), "dq") / DT
                           + value(E(a, n, i), "rowBeta")))
                    for (a, n, i) in ALL)
    require(worst_row == 0.0,
            f"the printed row closes as dpNetwork + M*dq/dt + rowBeta with zero "
            f"slack in all {len(ALL)} states -- the same identity 5.44 verified "
            "for two tables, now with the ddtFrom table live in three arms")
    bad_beta = [(a, NAME[i], value(E(a, n, i), "rowBeta"))
                for (a, n, i) in ALL
                if abs(value(E(a, n, i), "rowBeta")
                       - (value(E(a, n, i), "dpConv")
                          if (a, NAME[i]) in CHARGED else 0.0)) > 1e-9]
    n_ch = len([1 for x in ALL if (x[0], NAME[x[2]]) in CHARGED])
    require(not bad_beta,
            f"rowBeta is the whole account on the {n_ch} charged states and "
            f"exactly 0 on the {len(ALL) - n_ch} reported ones -- the account "
            "follows the table, including on the two quasi cells where the row "
            "pays M*dq/dt with a grade of 0"
            + (f"; first mismatch {bad_beta[0]}" if bad_beta else ""))

    # ===== 7. the controls, and what the run cost
    rel_gap = max(abs(value(E(a, n, i), "dpNetwork")
                      - (value(E(a, n, i), "dpFriction")
                         - value(E(a, n, i), "rowBeta")))
                  / abs(value(E(a, n, i), "dpFriction")
                        - value(E(a, n, i), "rowBeta"))
                  for (a, n, i) in ALL)
    require(rel_gap <= TOL,
            f"the chord bracket stays inside the interface's own tolerance in "
            f"every arm (worst {rel_gap:.4e} relative against the {TOL:.0e}) -- "
            "a quasi edge, a slab rung and a charged row together do not loosen "
            "the negotiation")
    bad1 = [(NAME[i], f, value(E("base", 1, i), f), PINET_L1[NAME[i]][f])
            for i in (0, 1) for f in PINET_L1[NAME[i]]
            if abs(value(E("base", 1, i), f) - PINET_L1[NAME[i]][f]) > 1e-12]
    require(not bad1,
            "base reproduces 5.41's published p_net level-1 states field for "
            "field -- the shipped numbers, so the three-table case starts where "
            "the two-table literature ended"
            + (f"; first mismatch {bad1[0]}" if bad1 else ""))
    bad8 = [(NAME[i], f, value(E("base", N_STEPS, i), f), PINET_L8[NAME[i]][f])
            for i in (0, 1) for f in PINET_L8[NAME[i]]
            if abs(value(E("base", N_STEPS, i), f) - PINET_L8[NAME[i]][f]) > 1e-12]
    require(not bad8,
            "and its level-8 states too, including the coefficient the "
            "interface converged to"
            + (f"; first mismatch {bad8[0]}" if bad8 else ""))
    seq = [len(by["base"]["pass"][n]) for n in range(1, N_STEPS + 1)]
    require(tuple(seq) == PINET_PASSES,
            f"and base works exactly as hard as the published arm level by "
            f"level ({seq}): the extra tables cost nothing when they say "
            "nothing")
    q_pin = [(a, n, value(E(a, n, 1), "q")) for a in ARMS
             for n in range(1, N_STEPS + 1)
             if abs(value(E(a, n, 1), "q") / Q_IN - 1.0) > 1e-12]
    require(not q_pin,
            "X2's rate is the injected 6.0e-5 to 1e-12 in all five arms at all "
            "eight levels -- the series edge is pinned by topology in every "
            "cell of the board (pin's level-1 split lands one ulp off, which "
            "is what a quasi free edge does to the split)"
            + (f"; first exception {q_pin[0]}" if q_pin else ""))
    pipe_lines = {}
    lo = None
    for k, ln in enumerate(lines):
        m = re.match(r"transient start (\w+): inertia", ln)
        if m:
            lo = m.group(1)
        m = re.match(r"pipes step (\d+):", ln)
        if m and lo:
            t = ln.split()
            pipe_lines.setdefault(lo, {})[int(m.group(1))] = dict(
                zip(t[3::2], t[4::2]))
    closed = max(abs(float(pipe_lines[a][n]["X1"]) + float(pipe_lines[a][n]["By"])
                     - Q_IN)
                 for a in ARMS for n in range(1, N_STEPS + 1))
    require(closed < 1e-18,
            f"and the branch still closes: X1 + By is the injected 6.0e-5 at "
            f"every level of every arm (worst {closed:.2e} m3/s)")
    solves = {a: by[a]["solves"] for a in ARMS}
    require(solves["mix"] == solves["mix_charge"]
            and solves["base"] == solves["both_chr"],
            f"the charge buys no extra block work ({solves}) -- mix and "
            f"mix_charge solve the same {solves['mix_charge']}, base and "
            f"both_chr the same {solves['both_chr']}; pin's {solves['pin']} is "
            "smaller because its quasi X1 ends the negotiation sooner, which "
            "is what a quasi partner is for")

    print(f"\n{len(fails)} failure(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
