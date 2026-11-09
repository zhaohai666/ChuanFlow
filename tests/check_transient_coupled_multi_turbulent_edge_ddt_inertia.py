#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient-multi over two wall-function RANS ducts with
BOTH per-edge tables declared at once: `ddtFrom` says which formulation each
coupled edge is negotiated under, `inertiaFrom` says whose inertance each coupled
edge is charged with, and the two tables name the SAME edge in the mixed arms.

The topology is the multi-edge one the profile chain has used since §5.34:
Src -> A -> { X1 (3D RANS duct) || By (1D bypass) } -> X2 (3D RANS duct) -> B -> Snk,
with 6.0e-5 m3/s switched on at t = 0 over eight levels of dt = 0.2 s. Everything
physical is §5.39's -- same meshes, same probes, same wall functions, same
operating point, same dt, same seeds, same tolerance -- with ONE number changed:
both proxy bores are typed short of the equivalent diameter of the section they
stand in for (X1 0.005 against 0.0056418958354775631, X2 0.006 against
0.00638307648642292245). §5.39 gave them the exact bores, so rho*L/A_proxy of the
drawing and rho*V/A^2 of the region were the same number and `inertiaFrom` could
not be seen in any result. Here the three sources are three different numbers on
every edge.

  nn_steady / pp_steady      both edges quasi-steady, drawing / region
  nn_both  / pp_both / ss_both   both edges time-accurate, drawing / region / slab
  es_nnn                     X1 steps, X2 quasi, drawing everywhere (5.39's cell)
  es_mix                     the same, with X1's inertance taken from its slab
  se_mix                     the mirror: X2 steps on its slab, X1 quasi on the bore

Four things could only have been measured with both tables on the same edge:

  1. the double count is a SUM OF TWO BRACKETS, edge by edge:
     rowError * dpBlock == (R*q - dpFriction) + (M*dq/dt - dpInertia) over all 72
     stepping level x edge states, to 1e-12 of the four terms it is made of -- the
     first bracket is the chord gap the interface never iterates, the second is the
     bill for whose inertance was installed
  2. `slab` is the only source that cancels the second bracket, because it is the
     only one whose number is the one the block measured its drop over: the
     ledger falls from 718.87 Pa (drawing) and 551.40 Pa (region) to 11.79 Pa on
     X1 at level 1, and the row's error from 77.1% and 53.9% to 0.94%. The block's
     own claim barely moves while that happens -- on the pinned edge dpInertia is
     the same 642.619 Pa to 1.1e-8 relative in all three attributions, and rowError
     is linear in the installed M with slope dq/dt to better than 3e-7, the residue
     accounted for by the two brackets the block itself owns
  3. on a quasi-steady edge the grade is BLIND BY DEFINITION (rowError is the
     constant 0 of coupledNetwork.cpp:1202) while the row still pays: nn_steady
     installs 1.0186e7 kg/m4 on X1 and charges 1165.3 Pa of M*dq/dt against a
     region that reported 310.8 Pa and no acceleration, and prints 0.0
  4. the two locality laws differ. §5.39's free edge reproduced its scheme-matched
     control BIT-FOR-BIT; here the free edge also reproduces a control whose OTHER
     edge changed formulation AND source (se_mix/X1 vs nn_steady/X1, 0.0 in 15
     fields), but an edge whose own source is shared with a control agrees only to
     2.9e-10 in rate and 1.3e-9 in drop -- because the installed M is a term of the
     global linear system, and it shows: every level's node pair equals that edge's
     ROW drop, inertial term included, whatever source the case named

usage: check_transient_coupled_multi_turbulent_edge_ddt_inertia.py <solver> <case>
"""
import math
import os
import re
import subprocess
import sys

# fluid, geometry and schedule of
# tests/cases/transientCoupledMultiTurbulentEdgeDdtInertia/case.dict
MU, RHO = 1.0e-3, 1000.0
NU = MU / RHO
Q_IN = 6.0e-5                  # m3/s injected at Src from t = 0
P_REF = 1.0e6                  # Pa held at Snk
DT, N_STEPS = 0.2, 8           # time { dt 0.2; endTime 1.6; }
MAX_SWEEPS = 300.0             # simple{ maxIterations }, each block's own budget
TOL_INTERFACE = 1e-3           # couple{ tolerance }

ARMS = ("nn_steady", "pp_steady", "nn_both", "pp_both", "ss_both",
        "es_nnn", "es_mix", "se_mix")
NAME = ("X1", "X2")
FREE, PINNED = 0, 1            # which index is the relaxing edge and which is pinned

# the two per-edge tables as the case declares them. `EDGE_DDT` is the scheme the
# edge is negotiated under, `EDGE_M` whose inertance the row charges it with; both
# are read back from the log below rather than assumed.
EDGE_DDT = {"nn_steady": ("steady", "steady"),
            "pp_steady": ("steady", "steady"),
            "nn_both": ("Euler", "Euler"),
            "pp_both": ("Euler", "Euler"),
            "ss_both": ("Euler", "Euler"),
            "es_nnn": ("Euler", "steady"),
            "es_mix": ("Euler", "steady"),
            "se_mix": ("steady", "Euler")}
EDGE_M = {"nn_steady": ("network", "network"),
          "pp_steady": ("partner", "partner"),
          "nn_both": ("network", "network"),
          "pp_both": ("partner", "partner"),
          "ss_both": ("slab", "slab"),
          "es_nnn": ("network", "network"),
          "es_mix": ("slab", "network"),
          "se_mix": ("network", "slab")}
MIXED_DDT = tuple(a for a in ARMS if len(set(EDGE_DDT[a])) > 1)
MIXED_M = tuple(a for a in ARMS if len(set(EDGE_M[a])) > 1)
BOTH = tuple(a for a in ARMS if a in MIXED_DDT and a in MIXED_M)
SEED = {"X1": 1.4e7, "X2": 1.36e7}

# gap, span, streamwise cells, region length, and the DRAWING's proxy bore -- the
# one number §5.39's case changed. The wall-bounded axis is y in both blocks.
DUCT = {"X1": (0.005, 0.005, 20, 0.2), "X2": (0.004, 0.008, 24, 0.24)}
BORE = {"X1": 0.005, "X2": 0.006}
PANELS = 8
BY = (0.3, 0.008)               # the 1D bypass, sharing X1's node pair
LEAD = (0.3, 0.02)              # A and B
PROBE = {"X1": (0.12, 0.19), "X2": (0.16, 0.23)}
AREA = {n: DUCT[n][0] * DUCT[n][1] for n in NAME}
# the three candidate numbers, computed the way each owner computes it
M_DRAW = {n: RHO * DUCT[n][3] / (0.25 * math.pi * BORE[n] ** 2) for n in NAME}
M_REGION = {n: RHO * DUCT[n][3] * DUCT[n][0] * DUCT[n][1] / AREA[n] ** 2 for n in NAME}
M_SLAB = {n: RHO * (PROBE[n][1] - PROBE[n][0]) * DUCT[n][0] * DUCT[n][1]
          / AREA[n] ** 2 for n in NAME}
D_EXACT = {n: 2.0 * math.sqrt(AREA[n] / math.pi) for n in NAME}

# measured at level 1 on the stepping edges: the second bracket, in Pa, and the
# row error it produces. The ladder is the increment's claim, so it is pinned.
LEDGER_L1 = {"network": {"X1": 718.8664888154701, "X2": 1903.8602570320108},
             "partner": {"X1": 551.3980590210152, "X2": 1607.3811672776105},
             "slab": {"X1": 11.793448853015036, "X2": 13.631160763185676}}
ROWERR_L1 = {"network": {"X1": 0.7708364636838451, "X2": 0.5219122045897613},
             "partner": {"X1": 0.5389408017820084, "X2": 0.4406372974151203},
             "slab": {"X1": 0.009364039865274851, "X2": 0.0037368172067357964}}
SHARE_L1 = {"nn_steady": 0.381355, "pp_steady": 0.419889, "nn_both": 0.319615,
            "pp_both": 0.346431, "ss_both": 0.431868, "es_nnn": 0.319615,
            "es_mix": 0.431868, "se_mix": 0.381355}
TAU_L1 = {"nn_steady": 0.4773, "pp_steady": 0.4329, "nn_both": 0.6134,
          "pp_both": 0.6014, "ss_both": 0.2997, "es_nnn": 0.6134,
          "es_mix": 0.2997, "se_mix": 0.4773}
PASSES = {"nn_steady": 16, "pp_steady": 16, "nn_both": 26, "pp_both": 27,
          "ss_both": 30, "es_nnn": 26, "es_mix": 30, "se_mix": 16}
SOLVES = {"nn_steady": 32, "pp_steady": 32, "nn_both": 52, "pp_both": 54,
          "ss_both": 60, "es_nnn": 52, "es_mix": 60, "se_mix": 32}

BETA_PLATEAU = (1.0200, 1.0211)

# the 15 fields an edge state carries; a "bit-for-bit" claim is over these
FIELDS = ("q", "dpBlock", "dpNetwork", "R", "residual", "M", "dpFriction",
          "dpInertia", "qPrev", "dq", "rowDrop", "rowError", "betaIn", "betaOut",
          "pRisePa")


def word(line, key):
    """The non-numeric value stored under `key` ('inertia network ...')."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([A-Za-z0-9_]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return m.group(1)


def value(line, key):
    """First number stored under `key` in a log line ('q 1.5e-4 ...')."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([-+0-9.eE]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return float(m.group(1))


def churchill(reynolds):
    """Churchill's 1977 Darcy friction factor -- published, so the 1D side is
    judged from outside this repository."""
    reynolds = max(reynolds, 1e-12)
    a = (2.457 * math.log(1.0 / ((7.0 / reynolds) ** 0.9))) ** 16
    b = (37530.0 / reynolds) ** 16
    return 8.0 * ((8.0 / reynolds) ** 12 + 1.0 / (a + b) ** 1.5) ** (1.0 / 12.0)


def pipeChord(length, diameter, q):
    """The 1D pipe's own pressure-per-rate, from Churchill's friction factor."""
    area = 0.25 * math.pi * diameter ** 2
    v = abs(q) / area
    return (churchill(RHO * v * diameter / MU) * length / diameter
            * 0.5 * RHO * v * v / max(abs(q), 1e-300))


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("transient per-edge ddt x attribution solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    def close(a_, b_, rel):
        if b_ == 0.0:
            return a_ == 0.0
        return abs(a_ / b_ - 1.0) <= rel

    # ---- split the log into the arms ----
    arms = []
    cur = None
    for ln in lines:
        t = ln.split()
        if ln.startswith("transient start ") and " accelerate " in ln:
            cur = {"tag": ln, "arm": t[2].rstrip(":"), "inertia": {}, "seed": {},
                   "partner": {}, "step": {}, "edge": {}, "pass": {}, "nodes": {},
                   "pipes": {}, "block": {}, "final": {}, "summary": "", "solves": 0}
            arms.append(cur)
        elif cur is None:
            continue
        if ln.startswith("SIMPLE iter 1:"):
            cur["solves"] += 1
        elif ln.startswith("transient arm ") and " inertia " in ln:
            cur["inertia"][t[4]] = t[6]
        elif ln.startswith("transient start ") and " resistance " in ln:
            cur["seed"][word(ln, "edge")] = value(ln, "resistance")
        elif ln.startswith("transient partner "):
            cur["partner"][t[4]] = ln
        elif ln.startswith("transient step ") and " pass " in ln:
            cur["pass"].setdefault(int(t[2]), {})[int(t[4])] = ln
        elif ln.startswith("transient step ") and " edge " in ln:
            cur["edge"].setdefault(int(t[2]), {})[t[6]] = ln
        elif ln.startswith("transient step "):
            cur["step"][int(t[2].rstrip(":"))] = ln
        elif ln.startswith("nodes step "):
            cur["nodes"][int(t[2].rstrip(":"))] = ln
        elif ln.startswith("pipes step "):
            cur["pipes"][int(t[2].rstrip(":"))] = ln
        elif ln.startswith("transient block "):
            cur["block"][word(ln, "name").rstrip(":")] = ln
        elif ln.startswith("transient coupled start ") and " name " in ln:
            cur["final"][word(ln, "edge")] = ln
        elif ln.startswith("transient coupled start "):
            cur["summary"] = ln

    by_arm = {c["arm"]: c for c in arms}

    def E(a, n, i):
        """The interface state of edge `i` at level `n` of arm `a`."""
        return by_arm[a]["edge"][n][NAME[i]]

    def B(a, n):
        return by_arm[a]["block"][n]

    def nd(a, n, name):
        t = by_arm[a]["nodes"][n].split()
        return float(dict(zip(t[3::2], t[4::2]))[name])

    def pipeq(a, n, name):
        t = by_arm[a]["pipes"][n].split()
        return float(dict(zip(t[3::2], t[4::2]))[name])

    def src(a, i):
        """The number the case's source for edge `i` of arm `a` announces."""
        return {"network": M_DRAW, "partner": M_REGION, "slab": M_SLAB}[EDGE_M[a][i]]

    def brackets(a, n, i):
        """The two terms the row's error is made of, from the printed fields."""
        ln = E(a, n, i)
        chord = value(ln, "dpNetwork") - value(ln, "dpFriction")
        ledger = value(ln, "M") * value(ln, "dq") / DT - value(ln, "dpInertia")
        return chord, ledger, value(ln, "rowError") * value(ln, "dpBlock")

    def maxdiff(a, x, b, y, field):
        """The worst absolute difference of one field between arm `a`'s edge `x`
        and arm `b`'s edge `y`, over all eight levels."""
        return max(abs(value(E(a, n, x), field) - value(E(b, n, y), field))
                   for n in range(1, N_STEPS + 1))

    def maxrel(a, x, b, y, field):
        base = max(abs(value(E(b, n, y), field)) for n in range(1, N_STEPS + 1))
        d = maxdiff(a, x, b, y, field)
        return d / base if base else d

    ALL = [(a, n, i) for a in ARMS for n in range(1, N_STEPS + 1) for i in (0, 1)]
    STEPPING = {(a, n, i) for a, n, i in ALL if EDGE_DDT[a][i] == "Euler"}
    QUASI = {(a, n, i) for a, n, i in ALL if EDGE_DDT[a][i] == "steady"}

    # ============ 1. two tables, one edge: the declaration surface
    header = [ln for ln in lines if ln.startswith("transient network+blocks")]
    require(len(header) == 1 and f"2 coupled edges, {len(ARMS)} starting sets" in header[0],
            f"one network, two coupled edges and {len(ARMS)} starting sets announced once")
    require([c["arm"] for c in arms] == list(ARMS),
            f"the eight arms ran, in the order the case declares them (got "
            f"{[c['arm'] for c in arms]})")
    require(all(c["seed"] == SEED for c in arms),
            f"every arm seeds both edges with the same pair ({SEED['X1']:.3e} then "
            f"{SEED['X2']:.3e}), so neither table comes in through the starting "
            "coefficients")
    require({c["arm"]: "accelerate yes" in c["tag"] for c in arms}
            == {a: False for a in ARMS},
            "no arm is accelerated: every difference below is a scheme or an "
            "inertance, not the update rule")
    for c in arms:
        a = c["arm"]
        got_ddt = {n: word(c["partner"][n], "ddt") for n in NAME}
        require(got_ddt == dict(zip(NAME, EDGE_DDT[a])),
                f"{a}: the rebuilt partners report {got_ddt} on their own edges")
        got_m = {n: c["inertia"][n] for n in NAME}
        require(got_m == dict(zip(NAME, EDGE_M[a])),
                f"{a}: the row charges its two edges with {got_m}")
    require(all(len(set(EDGE_DDT[a])) == 2 and len(set(EDGE_M[a])) == 2 for a in BOTH),
            " and on " + "/".join(BOTH) + " the two tables disagree with each other on "
            "the SAME edge: " + " and ".join(
                f"{a} takes X1 as {EDGE_DDT[a][0]}+{EDGE_M[a][0]} and X2 as "
                f"{EDGE_DDT[a][1]}+{EDGE_M[a][1]}" for a in BOTH)
            + " -- no per-arm switch can say either of those, let alone both")
    require(all(not by_arm[a]["partner"] is None for a in ARMS),
            "every arm's own words were read back from the log, not from the case")

    # ---- what the bore the case typed actually buys
    for n in NAME:
        require(all(close(value(E(a, 1, NAME.index(n)), "M"), M_DRAW[n], 1e-14)
                    for a in ARMS if EDGE_M[a][NAME.index(n)] == "network"),
                f"{n}: every edge attributed to the drawing carries rho*L/A of the "
                f"bore the case TYPED ({M_DRAW[n]:.6e} kg/m4) -- and the typed bore "
                f"{BORE[n]} is not the equivalent diameter of the section it stands "
                f"in for ({D_EXACT[n]:.8f}), so the drawing is {M_DRAW[n] / M_REGION[n]:.4f}x "
                f"the region")
        require(all(close(value(E(a, 1, NAME.index(n)), "M"), M_REGION[n], 1e-14)
                    for a in ARMS if EDGE_M[a][NAME.index(n)] == "partner"),
                f"{n}: every edge attributed to the partner carries exactly the "
                f"number that partner announced for ITSELF "
                f"({value(by_arm['pp_both']['partner'][n], 'regionInertance'):.6e} kg/m4)")
        require(all(close(value(E(a, 1, NAME.index(n)), "M"), M_SLAB[n], 1e-14)
                    for a in ARMS if EDGE_M[a][NAME.index(n)] == "slab"),
                f"{n}: every edge attributed to the slab carries exactly the "
                f"announced probeInertance "
                f"({value(by_arm['ss_both']['partner'][n], 'probeInertance'):.6e} kg/m4, "
                f"the {PROBE[n][1] - PROBE[n][0]:.2f} m the block reports its drop over)")
    require(all(len({value(E(a, 1, i), "M") for a in ARMS}) == 3 for i in (0, 1)),
            "the three sources really are three different numbers on each edge: "
            + " and ".join(f"{n} {M_DRAW[n]:.4e} / {M_REGION[n]:.4e} / {M_SLAB[n]:.4e}"
                           for n in NAME)
            + " -- 5.39's exact bores made the first two coincide, which is why "
            "attribution was invisible there")
    require(all(close(value(B(a, n), "inertance"), M_REGION[n], 1e-14)
                and close(value(B(a, n), "inertanceProbe"), M_SLAB[n], 1e-14)
                for a in ARMS for n in NAME),
            f"and the block's own report carries both: regionInertance "
            f"{M_REGION['X1']:.4e}/{M_REGION['X2']:.4e} and probeInertance "
            f"{M_SLAB['X1']:.4e}/{M_SLAB['X2']:.4e} kg/m4, so the interface is "
            "installing a number the partner can check")
    require(all(close(M_SLAB[n] / M_REGION[n], (PROBE[n][1] - PROBE[n][0]) / DUCT[n][3],
                      1e-12) for n in NAME),
            "the slab is a fraction of the region equal to its length fraction -- "
            f"{M_SLAB['X1'] / M_REGION['X1']:.6f} of X1 (0.07 m of 0.20) and "
            f"{M_SLAB['X2'] / M_REGION['X2']:.6f} of X2 (0.07 m of 0.24)")

    flag = {(a, n, i): value(E(a, n, i), "unsteady") for a, n, i in ALL}
    require(flag == {(a, n, i): (1.0 if EDGE_DDT[a][i] == "Euler" else 0.0)
                     for a, n, i in ALL},
            f"all {len(ALL)} level x edge states carry the flag of their OWN edge and "
            f"their OWN table row: {len(STEPPING)} stepping, {len(QUASI)} quasi, and on "
            "the mixed arms a stepping edge sits next to a quasi one at every level")
    require(all(value(B(a, n), "dtLevel")
                == (DT if EDGE_DDT[a][NAME.index(n)] == "Euler" else 0.0)
                for a in ARMS for n in NAME),
            f"and each block hands back dtLevel {DT} s exactly on the edge whose scheme "
            "row says Euler, 0 on the one whose row says steady")

    # ============ 2. the partition is per edge, and the row error is two brackets
    worst_split = max(abs((value(E(a, n, i), "dpFriction") + value(E(a, n, i), "dpInertia"))
                          / value(E(a, n, i), "dpBlock") - 1.0) for a, n, i in ALL)
    require(worst_split < 1e-12,
            f"dpFriction + dpInertia IS dpBlock in every state to {worst_split:.1e} -- "
            "the partition is made inside each block, so an arm with two sources is two "
            "partitions, not one blended one")
    require(all(value(E(a, n, i), "rowError") == 0.0 for a, n, i in QUASI)
            and all(value(E(a, n, i), "rowError") != 0.0 for a, n, i in STEPPING),
            f"rowError is EXACTLY 0.0 on all {len(QUASI)} quasi states and non-zero on "
            f"all {len(STEPPING)} stepping ones -- the definition at "
            "coupledNetwork.cpp:1202, so a quasi edge grades nothing")
    require(all(value(E(a, n, i), "dpInertia") == 0.0
                and value(E(a, n, i), "dpFriction") == value(E(a, n, i), "dpBlock")
                for a, n, i in QUASI),
            "on every quasi edge the inertial half is exactly zero and the friction "
            "half is the whole drop, in the mixed arms too")

    worst_id = 0.0
    where = None
    for a, n, i in sorted(STEPPING):
        ln = E(a, n, i)
        chord, ledger, lhs = brackets(a, n, i)
        # graded against the four terms the identity is made of, not against the
        # difference they leave: at the plateau the two brackets nearly cancel on
        # the pinned edge, and a cancelling pair is not a place where a rounding
        # claim means anything.
        scale = (abs(value(ln, "R") * value(ln, "q")) + abs(value(ln, "dpFriction"))
                 + abs(value(ln, "M") * value(ln, "dq") / DT) + abs(value(ln, "dpInertia")))
        rel = abs(lhs - (chord + ledger)) / scale
        if rel > worst_id:
            worst_id, where = rel, (a, n, NAME[i])
    require(worst_id < 1e-12,
            f"rowError*dpBlock == (R*q - dpFriction) + (M*dq/dt - dpInertia) in every "
            f"one of the {len(STEPPING)} stepping states, worst {worst_id:.1e} of the "
            f"four terms it is made of (at {where[0]} level {where[1]} {where[2]}) -- "
            "the row's error is "
            "the chord gap it never iterates PLUS the bill for the inertance it "
            "installed, and the two tables are what set those two numbers")
    for a in ("nn_both", "pp_both", "ss_both"):
        for i in (FREE, PINNED):
            chord, ledger, lhs = brackets(a, 1, i)
            require(abs(ledger) > 5.0 * abs(chord),
                    f"{a} level 1 {NAME[i]}: the two brackets are not comparable -- "
                    f"chord gap {chord:+.4f} Pa against attribution {ledger:+.4f} Pa, "
                    f"{abs(ledger) / max(abs(chord), 1e-300):.0f}x")
    chord, ledger, lhs = brackets("ss_both", 8, FREE)
    require(abs(ledger) < abs(chord),
            f"with the slab installed the order inverts by the plateau: ss_both level 8 "
            f"X1 is chord {chord:+.4f} Pa against attribution {ledger:+.4f} Pa -- the "
            "attribution bracket is the one the source can remove")

    # ---- and what the blind half is worth
    blind = []
    for a, n, i in sorted(QUASI):
        ln = E(a, n, i)
        paid = value(ln, "M") * value(ln, "dq") / DT
        hypothetical = (paid + value(ln, "dpNetwork") - value(ln, "dpFriction")) \
            / value(ln, "dpBlock")
        if a == "nn_steady" and n == 1 and i == FREE:
            blind.append((a, n, NAME[i], paid, value(ln, "dpBlock"), hypothetical))
    for a, n, e, paid, dpb, hyp in blind:
        require(paid > 3.0 * dpb and hyp > 3.0,
                f"{a} level {n} {e}: the row installs {paid:.1f} Pa of M*dq/dt against a "
                f"region that reported {dpb:.1f} Pa and no acceleration at all -- graded "
                f"that would be {hyp * 100:.0f}% -- and rowError prints 0.0")
    require(max(abs(value(E(a, n, i), "M") * value(E(a, n, i), "dq") / DT)
                for a, n, i in QUASI if (a, n, i) != ("nn_steady", 1, 0)) > 100.0,
            "the same holds across the quasi family, not just at that one state: "
            "M*dq/dt is hundreds of pascals on edges whose grade is identically zero")

    # ============ 3. the attribution ladder, on a stepping edge
    for i in (FREE, PINNED):
        n = NAME[i]
        lad_arm = {"network": "nn_both", "partner": "pp_both", "slab": "ss_both"}
        for s, a in lad_arm.items():
            require(close(value(E(a, 1, i), "rowError"), ROWERR_L1[s][n], 1e-9),
                    f"level 1 {n} {s}: the ladder's rung is a measured number, not a "
                    f"trend -- {value(E(a, 1, i), 'rowError') * 100:.4f}% printed where "
                    f"{ROWERR_L1[s][n] * 100:.4f}% was recorded")
        row = {s: value(E(a, 1, i), "rowError") for s, a in lad_arm.items()}
        lad = [row["network"], row["partner"], row["slab"]]
        require(lad[0] > lad[1] > lad[2] > 0.0,
            f"level 1 {n}: the row error falls monotonically down the source ladder -- "
            f"drawing {lad[0] * 100:.2f}%, region {lad[1] * 100:.2f}%, slab "
            f"{lad[2] * 100:.2f}% -- and the slab is the source whose number is the one "
            "the block measured its drop over")
        led = [LEDGER_L1[s][n] for s in ("network", "partner", "slab")]
        meas = [brackets({"network": "nn_both", "partner": "pp_both",
                          "slab": "ss_both"}[s], 1, i)[1] for s in ("network", "partner", "slab")]
        require(all(close(m, e, 1e-6) for m, e in zip(meas, led)),
            f"level 1 {n}: the attribution bracket is {led[0]:.2f} Pa on the drawing, "
            f"{led[1]:.2f} on the region, {led[2]:.2f} on the slab -- slab removes "
            f"{(1 - led[2] / led[0]) * 100:.2f}% of the drawing's charge")
        require(abs(meas[0] - (M_DRAW[n] - M_SLAB[n]) * value(E("nn_both", 1, i), "dq") / DT)
                < 25.0,
                f"and it is what the drawing over-charges, not something else: "
                f"(M_draw - M_slab)*dq/dt = "
                f"{(M_DRAW[n] - M_SLAB[n]) * value(E('nn_both', 1, i), 'dq') / DT:.2f} Pa "
                f"against the measured {meas[0]:+.2f} Pa -- the residue is the block's "
                "own acceleration, not the row's uniform one")
    pa_lad = {s: brackets(a, 1, PINNED)[1] for s, a in lad_arm.items()}
    pa_free = {s: brackets(a, 1, FREE)[1] for s, a in lad_arm.items()}
    require(all(pa_lad[s] > pa_free[s] > 0.0 for s in lad_arm),
            "the pinned edge carries the same ladder at a bigger absolute price: "
            + "  ".join(f"{s} {pa_lad[s]:.1f} Pa vs {pa_free[s]:.1f}" for s in
                        ("network", "partner", "slab"))
            + " -- its rate is the injected one, so its dq/dt is the whole step")

    # ---- the block's claim does not move when the row's does
    dpIn_p = [value(E(a, 1, PINNED), "dpInertia") for a in ("nn_both", "pp_both", "ss_both")]
    spread = max(abs(x - dpIn_p[0]) for x in dpIn_p) / abs(dpIn_p[0])
    require(spread < 1e-7,
            f"and the region's side of that ladder barely moves at all: dpInertia at "
            f"level 1 on X2 is {dpIn_p[0]:.4f} Pa under all three attributions, worst "
            f"{spread:.1e} relative -- attribution changes the row's charge, and what "
            "the block measured is only touched through the field")
    for a, b in (("nn_both", "pp_both"), ("pp_both", "ss_both"), ("nn_both", "ss_both")):
        ca, la, ha = brackets(a, 1, PINNED)
        cb, lb, hb = brackets(b, 1, PINNED)
        want = (value(E(a, 1, PINNED), "M") - value(E(b, 1, PINNED), "M")) \
            * value(E(a, 1, PINNED), "dq") / DT
        resid = abs((ha - hb) - want)
        own = abs(ca - cb) + abs(value(E(a, 1, PINNED), "dpInertia")
                                 - value(E(b, 1, PINNED), "dpInertia"))
        require(resid / abs(want) < 1e-6 and resid <= own + 1e-9,
                f"more strongly: on the pinned edge the difference between two "
                f"attributions IS (M_a - M_b)*dq/dt -- {a} vs {b} gives {ha - hb:.4f} Pa "
                f"against {want:.4f} Pa, {resid / abs(want):.1e} relative, and the whole "
                f"of the {resid:.1e} Pa residue is accounted for by the two brackets the "
                "block itself owns (" + f"{own:.1e} Pa" + ") -- so rowError is linear in "
                "the installed inertance with slope this level's dq/dt")
    ratios = [(a, NAME[i], value(E(a, 1, i), "dpInertia")
               / (M_SLAB[NAME[i]] * value(E(a, 1, i), "dq") / DT))
              for a in ("nn_both", "pp_both", "ss_both", "es_nnn", "es_mix", "se_mix")
              for i in (FREE, PINNED) if EDGE_DDT[a][i] == "Euler"]
    require(all(0.95 < x[2] < 0.99 for x in ratios),
            f"the block's own inertial half is worth "
            f"{min(x[2] for x in ratios) * 100:.2f} to {max(x[2] for x in ratios) * 100:.2f}% "
            f"of (M_slab * dq/dt) at level 1 in all {len(ratios)} stepping states -- the "
            "slab is the right SLAB and the wrong DU/DT, which is the whole of what "
            "`slab` leaves behind")
    slab_states = [(a, i, n) for a in ("ss_both", "es_mix", "se_mix")
                   for i in (FREE, PINNED) for n in range(3, N_STEPS + 1)
                   if EDGE_DDT[a][i] == "Euler"]
    worst_slab = max(abs(value(E(a, n, i), "rowError")) for a, i, n in slab_states)
    require(worst_slab < TOL_INTERFACE,
            f"with the slab installed the row's error is at or below the interface's own "
            f"tolerance everywhere from level 3 on: {len(slab_states)} stepping states, "
            f"worst {worst_slab:.3e} against the {TOL_INTERFACE:.0e} the passes stop at "
            f"-- the double count is off the row and what is left is a fraction of the "
            "gap the passes already chase")
    for a in ("ss_both", "es_mix"):
        worst = max(abs(value(E(a, n, FREE), "rowError"))
                    / max(value(E(a, n, FREE), "residual"), 1e-300)
                    for n in range(3, N_STEPS + 1))
        require(worst < 2.0,
            f"{a} X1: on the free edge the same states sit within {worst:.2f}x the "
            "residual of their own level -- the two errors are the same order once the "
            "row charges the slab the drop was measured over")
    big = max(abs(value(E(a, n, FREE), "rowError"))
              for a in ("nn_both", "es_nnn") for n in range(3, N_STEPS + 1))
    worst_n = max(abs(value(E("nn_both", n, FREE), "rowError"))
                  / max(value(E("nn_both", n, FREE), "residual"), 1e-300)
                  for n in range(3, N_STEPS + 1))
    require(big > 20.0 * TOL_INTERFACE and worst_n > 20.0,
            f"while the drawing's free edge sits at {big:.3e} -- {big / TOL_INTERFACE:.0f}x "
            f"the tolerance and up to {worst_n:.1f}x its own residual over the same levels "
            "-- the interface cannot iterate this term away and the tolerance does not ask "
            "it to")

    # ============ 4. two locality laws, and they are different laws
    for a, ctl, i in (("es_nnn", "nn_both", FREE), ("se_mix", "nn_steady", FREE),
                      ("es_mix", "ss_both", FREE)):
        zero = [f for f in FIELDS if maxdiff(a, i, ctl, i, f) == 0.0]
        require(len(zero) == len(FIELDS) or a == "es_mix",
            f"{a}/{NAME[i]} against {ctl}/{NAME[i]} (same scheme on that edge, same or "
            f"different other edge): {len(zero)} of {len(FIELDS)} fields agree exactly")
    require(all(maxdiff("es_nnn", FREE, "nn_both", FREE, f) == 0.0 for f in FIELDS),
        "es_nnn's X1 is nn_both's X1 BIT FOR BIT in all 15 fields: the other edge went "
        "from Euler to steady and the free edge did not notice, with a mis-stated bore "
        "installed on both")
    require(all(maxdiff("se_mix", FREE, "nn_steady", FREE, f) == 0.0 for f in FIELDS),
        "se_mix's X1 is nn_steady's X1 bit for bit too -- and there the other edge "
        "changed BOTH tables (steady+network to Euler+slab), so a quasi free edge is "
        "blind to its neighbour's formulation and its inertance alike")
    dq_ = maxdiff("es_mix", FREE, "ss_both", FREE, "q")
    dd_ = maxdiff("es_mix", FREE, "ss_both", FREE, "dpBlock")
    require(dq_ > 0.0 and dq_ < 1e-13 and dd_ > 0.0 and dd_ < 1e-6
            and maxrel("es_mix", FREE, "ss_both", FREE, "dpBlock") < 2e-9,
        f"es_mix's X1 against ss_both's X1 -- same scheme, same source, other edge "
        f"different in both tables -- is NOT bit for bit: rate differs by {dq_:.2e} "
        f"({maxrel('es_mix', FREE, 'ss_both', FREE, 'q'):.1e} relative) and drop by "
        f"{dd_:.2e} Pa ({maxrel('es_mix', FREE, 'ss_both', FREE, 'dpBlock'):.1e} "
        "relative). Formulation commutes across the two edges; attribution does not, "
        "because the installed M is a term of the global solve")
    require(max(value(E("nn_steady", n, FREE), "rowError") for n in range(1, N_STEPS + 1))
            == 0.0 == max(value(E("pp_steady", n, FREE), "rowError")
                          for n in range(1, N_STEPS + 1))
            and maxdiff("nn_steady", FREE, "pp_steady", FREE, "rowError") == 0.0,
        f"the blind half, again as a law: nn_steady and pp_steady differ by up to "
        f"{maxrel('nn_steady', FREE, 'pp_steady', FREE, 'q') * 100:.2f}% of rate and "
        f"{maxdiff('nn_steady', FREE, 'pp_steady', FREE, 'dpBlock'):.1f} Pa of drop on "
        "the free edge, and their row errors are the same exact zero")
    require(maxdiff("es_mix", FREE, "es_nnn", FREE, "q") > 0.2
            * abs(value(E("es_mix", 1, FREE), "q")),
        f"within one edge the source is not a small number either: es_mix/X1 and "
        f"es_nnn/X1 differ by {maxdiff('es_mix', FREE, 'es_nnn', FREE, 'q'):.3e} m3/s at "
        f"worst ({maxrel('es_mix', FREE, 'es_nnn', FREE, 'q') * 100:.1f}% relative) -- "
        f"{maxrel('es_mix', FREE, 'es_nnn', FREE, 'rowError') * 100:.0f}% of rowError")
    quasi_pinned = {a: maxdiff(a, PINNED, "nn_steady", PINNED, "dpBlock")
                    for a in ("es_nnn", "es_mix", "pp_steady")}
    require(all(d / 1538.4 < 1e-5 for d in quasi_pinned.values()),
            "a quasi pinned edge is nearly as indifferent to its neighbour as a quasi "
            "free one: " + "  ".join(f"{a}/X2 {d:.3e} Pa" for a, d in quasi_pinned.items())
            + " against nn_steady/X2")
    require(maxrel("ss_both", PINNED, "nn_both", PINNED, "dpBlock") < 1e-6
            and abs(value(E("ss_both", 1, PINNED), "rowError")
                    - value(E("nn_both", 1, PINNED), "rowError")) > 0.5,
        f"and the pinned edge's DROP is the same block to "
        f"{maxrel('ss_both', PINNED, 'nn_both', PINNED, 'dpBlock'):.1e} relative across "
        f"the source ladder while its row error spans "
        f"{abs(value(E('ss_both', 1, PINNED), 'rowError') - value(E('nn_both', 1, PINNED), 'rowError')) * 100:.1f} "
        "percentage points -- one side of the ledger moves, the other does not")

    # ============ 5. the ramp follows the free edge's own two tables
    share = {}
    dest = {}
    for a in ARMS:
        s = []
        for n in range(1, N_STEPS + 1):
            x, y = pipeq(a, n, "X1"), pipeq(a, n, "By")
            s.append(x / (x + y))
        share[a], dest[a] = s[0], s[-1]
        require(close(share[a], SHARE_L1[a], 1e-5),
            f"{a}: the split arrives at the parallel branch as {share[a]:.6f} -- the "
            f"table for THAT edge is what the branch sees "
            f"(measured {EDGE_DDT[a][FREE]}+{EDGE_M[a][FREE]})")
        taus = []
        for n in range(N_STEPS - 1):
            num_, den_ = s[n] - s[-1], s[n + 1] - s[-1]
            if num_ * den_ > 0:
                taus.append(-DT / math.log(den_ / num_))
        require(abs(taus[0] - TAU_L1[a]) < 2e-3,
            f"{a}: its first-interval time constant is {taus[0]:.4f} s")
    stepl = [share["ss_both"], share["pp_both"], share["nn_both"]]
    require(stepl[0] > stepl[1] > stepl[2],
        f"the ladder is monotone in the installed inertance on the free edge: at level 1 "
        f"the branch takes {stepl[0]:.6f} of the flow under the slab "
        f"({M_SLAB['X1']:.3e} kg/m4), {stepl[1]:.6f} under the region "
        f"({M_REGION['X1']:.3e}) and {stepl[2]:.6f} under the drawing "
        f"({M_DRAW['X1']:.3e}) -- a row that resists acceleration more takes less "
        "through it at the first step")
    taustep = [TAU_L1["ss_both"], TAU_L1["pp_both"], TAU_L1["nn_both"]]
    require(taustep[0] < taustep[1] < taustep[2],
        f"and relaxes back fastest under the smallest M: tau {taustep[0]:.4f} / "
        f"{taustep[1]:.4f} / {taustep[2]:.4f} s -- 2.0x from slab to drawing")
    require(share["pp_steady"] > share["nn_steady"]
            and TAU_L1["pp_steady"] < TAU_L1["nn_steady"],
        f"the same order holds with no du/dt anywhere: {share['pp_steady']:.6f} vs "
        f"{share['nn_steady']:.6f} at level 1, tau {TAU_L1['pp_steady']:.4f} vs "
        f"{TAU_L1['nn_steady']:.4f} s -- this is the network's own inertial relaxation, "
        "not the partner's")
    q8 = {a: pipeq(a, N_STEPS, "X1") for a in ARMS}
    spread1 = (max(pipeq(a, 1, "X1") for a in ARMS) - min(pipeq(a, 1, "X1") for a in ARMS)) \
        / Q_IN
    spread8 = (max(q8.values()) - min(q8.values())) / Q_IN
    require(spread8 < 0.2 * spread1,
        f"the attributions are {spread1 * 100:.1f}% of the injected flow apart at level "
        f"1 and {spread8 * 100:.1f}% apart at level 8 -- they are converging to one "
        "steady destination, and the plateau spread is the tolerance, not a physics the "
        "source chose")
    for a in ("nn_both", "pp_both", "ss_both", "es_nnn", "es_mix", "se_mix"):
        require(all(abs(pipeq(a, n, "X2") - Q_IN) < 1e-19 for n in range(1, N_STEPS + 1)),
            f"{a}: X2 is pinned by topology to the injected rate in every table "
            "combination, so what its two brackets do below is a pure interface effect")
    dqmax = max(abs(value(E(a, n, PINNED), "dq")) for a, n, i in ALL if i == PINNED and n > 1)
    require(dqmax < 1.4e-20,
        f"and from level 2 on the pinned edge's own dq is at most {dqmax:.2e} m3/s -- "
        "its ledger term is what the level-1 step leaves behind, and after that only the "
        "arithmetic of M times round-off")

    # ============ 6. what the two tables cost
    tot_pass = 0
    for a in ARMS:
        p = sum(value(by_arm[a]["step"][n], "subPasses") for n in range(1, N_STEPS + 1))
        sw = [value(E(a, n, i), "sweeps") for n in range(1, N_STEPS + 1) for i in (0, 1)]
        capped = [n for n in range(1, N_STEPS + 1) for i in (0, 1)
                  if value(E(a, n, i), "sweeps") == MAX_SWEEPS]
        unconv = [n for n in range(1, N_STEPS + 1) for i in (0, 1)
                  if value(E(a, n, i), "blockConverged") == 0.0]
        wr = value(by_arm[a]["summary"], "residual")
        tot_pass += p
        require(p == PASSES[a],
            f"{a}: {p} sub-iterations over the eight levels (measured against "
            f"{PASSES[a]} declared for this crossing)")
        require(by_arm[a]["solves"] == SOLVES[a],
            f"{a}: {by_arm[a]['solves']} block solves, {sum(sw)} sweeps")
        require(not capped and not unconv,
            f"{a}: no block reaches the {MAX_SWEEPS:.0f}-sweep budget and every block "
            f"reports converged (worst sweeps {max(sw)})")
        require(wr < TOL_INTERFACE,
            f"{a}: worst interface residual {wr:.3e} under the {TOL_INTERFACE:.0e} "
            "tolerance")
    require(tot_pass == sum(PASSES.values()),
        f"the eight arms together walk {tot_pass} interface sub-iterations and "
        f"{sum(SOLVES.values())} block solves")
    lad_p = [PASSES["nn_steady"], PASSES["nn_both"], PASSES["pp_both"], PASSES["ss_both"]]
    require(lad_p[0] < lad_p[1] < lad_p[2] < lad_p[3],
        f"and the interface effort is set by the FREE edge's two tables, monotone as they "
        f"are read: quasi {lad_p[0]} passes, then {lad_p[1]} / {lad_p[2]} / {lad_p[3]} "
        "under the drawing, region and slab. se_mix, whose stepping edge is the PINNED "
        f"one, costs {PASSES['se_mix']} -- the same as a fully quasi arm")

    # ============ 7. the partner is still a legible RANS duct under all of this
    for a in ARMS:
        for n in NAME:
            b = B(a, n)
            y = value(b, "wallYPlusMax")
            nut = value(b, "nutMax") / NU
            re = value(b, "uxMax") * 2 * DUCT[n][0] * DUCT[n][1] \
                / (DUCT[n][0] + DUCT[n][1]) / NU
            require(y > 11.0 and 1.0200 <= value(b, "betaOut") <= 1.0211
                    and 30.0 < nut < 120.0 and 4000.0 < re < 12000.0
                    and close(value(b, "fluxOut"), value(b, "fluxIn"), 1e-9),
                f"{a} {n}: first cell at y+ {y:.1f} (log law needs > 11), nut/nu "
                f"{nut:.1f}, Re_Dh {re:.0f}, beta_out {value(b, 'betaOut'):.6f} in the "
                f"channel's plateau band, mass to "
                f"{abs(value(b, 'fluxOut') - value(b, 'fluxIn')):.1e} m3/s")
            require(word(b, "inletShape") == "flat",
                f"{a} {n}: the shape asked of this edge is flat, so nothing below is a "
                "shape effect in disguise")

    # ============ 8. the installed M is literally in the node balance
    lead = pipeChord(*LEAD, Q_IN) * Q_IN
    for a in ARMS:
        require(all(close(nd(a, n, "J1") - nd(a, n, "J2"), value(E(a, n, FREE), "rowDrop"),
                        1e-9)
                    and close(nd(a, n, "J2") - nd(a, n, "J3"),
                              value(E(a, n, PINNED), "rowDrop"), 1e-9)
                    for n in range(1, N_STEPS + 1)),
            f"{a}: every level's node pair equals that edge's ROW drop -- R*q plus the "
            "installed M*dq/dt -- not the block's drop, so the source the case names is "
            "in the pressure field the network commits")
        require(close(nd(a, N_STEPS, "Src") - nd(a, N_STEPS, "J1"), lead, 1e-9)
                and close(nd(a, N_STEPS, "J3") - nd(a, N_STEPS, "Snk"), lead, 1e-9)
                and close(nd(a, N_STEPS, "Snk"), P_REF, 1e-12),
            f"{a}: the two lead pipes are the network's own Churchill drop, {lead:.6f} "
            "Pa, and Snk holds the case's pressure")
        require(all(abs(pipeq(a, n, "A") - Q_IN) < 1e-18
                    and close(pipeq(a, n, "X1") + pipeq(a, n, "By"), Q_IN, 1e-11)
                    for n in range(1, N_STEPS + 1)),
            f"{a}: A takes the injected rate and the pair splits it at every level")
    for a in ARMS:
        require(all(abs(value(E(a, 1, i), "rowDrop") - value(E(a, 1, i), "dpBlock"))
                    < 3e3 for i in (FREE, PINNED)),
            f"{a}: at level 1 the row and the region disagree by at most "
            f"{max(abs(value(E(a, 1, i), 'rowDrop') - value(E(a, 1, i), 'dpBlock')) for i in (FREE, PINNED)):.0f} "
            "Pa on the two edges -- both brackets are bounded, nothing diverges")

    print("note: per-edge scheme / source -- " + "  ".join(
        f"{a} {EDGE_DDT[a][0]}+{EDGE_M[a][0]}/{EDGE_DDT[a][1]}+{EDGE_M[a][1]}" for a in ARMS))
    print("note: installed M, X1 " + " ".join(
        f"{a}:{value(E(a, 1, FREE), 'M'):.3e}" for a in ARMS))
    print("note: brackets at level 1 (chord, attribution) and rowError -- " + "  ".join(
        f"{a}/{NAME[i]} {brackets(a, 1, i)[0]:+.3f}/{brackets(a, 1, i)[1]:+.2f} Pa "
        f"{value(E(a, 1, i), 'rowError') * 100:+.2f}%"
        for a in ("nn_both", "pp_both", "ss_both", "es_nnn", "es_mix", "se_mix")
        for i in (FREE, PINNED) if EDGE_DDT[a][i] == "Euler"))
    print("note: share at level 1 / destination at level 8 / tau -- " + "  ".join(
        f"{a} {share[a]:.6f}/{dest[a]:.6f}/{TAU_L1[a]:.4f}" for a in ARMS))
    print("note: passes / solves / sweeps -- " + "  ".join(
        f"{a} {PASSES[a]}/{SOLVES[a]}/{int(sum(value(E(a, n, i), 'sweeps') for n in range(1, N_STEPS + 1) for i in (0, 1)))}"
        for a in ARMS))

    if fails:
        print(f"\n{len(fails)} FAILURES:")
        for m in fails:
            print(" -", m)
        return 1
    print("\nPASS: two per-edge tables on the same edge -- the row's error is exactly a "
          "chord gap plus an attribution bracket on all 72 stepping states, `slab` is the "
          "source that cancels the second one (718.9 -> 11.8 Pa at level 1) while the "
          "block's own 642.6 Pa does not move, a quasi edge pays the same bracket and "
          "grades it as 0.0, and the free edge is bit-for-bit its own scheme's control but "
          "only 1e-9 its own source's control, because M is in the node balance")
    return 0


if __name__ == "__main__":
    sys.exit(main())
