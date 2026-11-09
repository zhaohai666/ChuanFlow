#!/usr/bin/env python3
r"""Run chflow-coupled-1d3d-transient-multi over two wall-function RANS ducts with
BOTH per-edge tables declared at once AND a shaped inlet, and answer the question
neither 5.39 nor 5.40 could: of the two brackets that make up a stepping edge's
row error, which one does the INLET SHAPE act on?

  Src -> A -> { X1 (3D RANS duct) || By (1D pipe) } -> X2 (3D RANS duct) -> B -> Snk

Everything physical is 5.40's -- same meshes, same probes, same wall functions,
same operating point, same mis-stated proxy bores, same dt, same seeds, same
tolerance -- plus `starts.profile`, so the interface writes a Poiseuille entry on
every pass of every level. 5.40 fixed the structure of the row:

  rowError * dpBlock = (R*q - dpFriction) + (M*dq/dt - dpInertia)
  \_______ chord ___/   \____ attribution bill ____/

the first bracket the gap the interface never iterates (chord is the block's
friction half on a stepping edge), the second the bill for whose inertance was
installed. Its boundary 2 recorded that the shape axis was never crossed with
attribution, and 5.39's that `parab_x2pin` was never run at all; its boundary 1
recorded that a quasi edge charged with `slab` was legal but unmeasured.

  f_net / f_region / f_slab   flat ladder: drawing / region / slab, both stepping
  p_net / p_region            the same two rungs shaped
  p_mix                       X1 steps on the drawing, X2 quasi on its slab
  p_pin                       X1 quasi on the drawing, X2 steps on its slab

What it settles, all measured:

  1. the identity is shape-blind -- it holds over every stepping state of both
     shapes, normalised by the four terms it is made of -- and the two brackets
     stay in different hands: the chord is at most ~0.64 Pa everywhere (tolerance-
     sized, whatever the shape), while the bill is 13 to 1903 Pa
  2. shape moves the BLOCK's half, not the bill: on the pinned edge, where dq/dt is
     topology-given, the bill is 1903.86 Pa flat and 1903.42 Pa shaped (0.02%),
     while the block's own dpInertia is 642.6188 Pa flat against 643.0614 Pa shaped
     and its dpBlock falls 21%. On the free edge the bill moves +4.9% and so does
     the rate it is charged with (+5.0%) -- the bill is arithmetic in the installed
     M and in the trajectory, and shape reaches it only through the trajectory
  3. the attribution ladder is shape-blind in the fraction it removes: X1 30.1%
     flat (77.08 -> 53.89) against 29.3% shaped (88.62 -> 62.66), X2 15.6% against
     15.6%, and the shaped region rung reproduces 5.39's recorded parab_both row
     errors digit for digit (62.6559% / 56.0164% against the 62.66% / 56.01% 5.39
     wrote down with bores made equal), because 5.39's `network` IS this case's
     `partner`
  4. 5.40's boundary 1 is measured: p_mix's quasi X2 carries its slab and the row
     pays M_slab*dq/dt = 656.25 Pa against a block that reported no acceleration at
     all, while its grade is exactly 0.0 -- and the drawing would have charged the
     same edge 2546.5 Pa, so the source removes 74.2% of a bill the grade cannot
     see. p_pin's quasi X1 is the same story at 1264.30 Pa
  5. the 300-sweep wall is shape's but attribution prices it: no flat state touches
     it (worst 16 sweeps), while the shaped pinned edge sits at 300 under the drawing
     and peaks at 162 under the region. In this crossing the nine states that fill
     the budget ARE the nine that report blockConverged 0 -- they coincide, all on
     the pinned edge, and the interface still finishes every one of them inside its
     own tolerance, so 'filled the budget' stays a record and not a stop rule
  6. the locality asymmetry of 5.40 survives the shape with a new trigger: p_mix's
     free edge is bit-identical to p_net's free edge in all 16 fields -- including
     the convective account -- for levels 1..5 and separates only from level 6 --
     the level after p_net's pinned block first fills its sweep budget -- while an
     edge that changes its OWN scheme row (p_pin/X1, quasi, against p_net/X1) moves
     23% at the first step
  7. the account the row still cannot state: every one of the 112 states prints
     dpConv = rho*(betaOut - betaIn)*q|q|/A^2, the convective momentum-flux
     imbalance the region owes its own changing profile -- sign set by the shape
     (48 flat states positive, 64 shaped negative), size set by the rate (dpConv/q^2
     is table-blind: the three flat arms agree to 5.7e-8 on the pinned edge and to
     under 8% on the free, where the rate itself spans 35%), 4.06x larger shaped
     than flat, up to 121% of the block's own drop on p_pin/X1 level 8, and at most
     17% of -- sometimes the opposite sign to -- the 5.40 attribution bill. The row
     is still exactly R*q + M*dq/dt on every stepping state: the number is on every
     line and in no equation

usage: check_transient_coupled_multi_turbulent_profile_inertia.py <solver> <case>
"""
import math
import os
import re
import subprocess
import sys

MU, RHO = 1.0e-3, 1000.0
NU = MU / RHO
Q_IN = 6.0e-5
P_REF = 1.0e6
DT, N_STEPS = 0.2, 8
MAX_SWEEPS = 300.0
TOL_INTERFACE = 1e-3

ARMS = ("f_net", "f_region", "f_slab", "p_net", "p_region", "p_mix", "p_pin")
NAME = ("X1", "X2")
EDGE_DDT = {"f_net": ("Euler", "Euler"), "f_region": ("Euler", "Euler"),
            "f_slab": ("Euler", "Euler"), "p_net": ("Euler", "Euler"),
            "p_region": ("Euler", "Euler"), "p_mix": ("Euler", "steady"),
            "p_pin": ("steady", "Euler")}
EDGE_M = {"f_net": ("network", "network"), "f_region": ("partner", "partner"),
          "f_slab": ("slab", "slab"), "p_net": ("network", "network"),
          "p_region": ("partner", "partner"), "p_mix": ("network", "slab"),
          "p_pin": ("network", "slab")}
EDGE_SHAPE = {"f_net": ("flat", "flat"), "f_region": ("flat", "flat"),
              "f_slab": ("flat", "flat"), "p_net": ("parabolic", "parabolic"),
              "p_region": ("parabolic", "parabolic"),
              "p_mix": ("parabolic", "parabolic"), "p_pin": ("parabolic", "parabolic")}
FLAT = tuple(a for a in ARMS if EDGE_SHAPE[a][0] == "flat")
SHAPED = tuple(a for a in ARMS if EDGE_SHAPE[a][0] != "flat")
SEED = {"X1": 1.4e7, "X2": 1.36e7}

DUCT = {"X1": (0.005, 0.005, 20, 0.2), "X2": (0.004, 0.008, 24, 0.24)}
BORE = {"X1": 0.005, "X2": 0.006}
PROBE = {"X1": (0.12, 0.19), "X2": (0.16, 0.23)}
AREA = {n: DUCT[n][0] * DUCT[n][1] for n in NAME}
M_DRAW = {n: RHO * DUCT[n][3] / (0.25 * math.pi * BORE[n] ** 2) for n in NAME}
M_REGION = {n: RHO * DUCT[n][3] * DUCT[n][0] * DUCT[n][1] / AREA[n] ** 2
            for n in NAME}
M_SLAB = {n: RHO * (PROBE[n][1] - PROBE[n][0]) * DUCT[n][0] * DUCT[n][1]
          / AREA[n] ** 2 for n in NAME}
SRC_M = {"network": M_DRAW, "partner": M_REGION, "slab": M_SLAB}
D_EXACT = {n: 2.0 * math.sqrt(AREA[n] / math.pi) for n in NAME}

# level 1, measured: the second bracket in Pa and the row error it produces. The
# three flat rungs are 5.40's nn_both / pp_both / ss_both recomputed here under
# other arm names, and the shaped region rung is 5.39's parab_both.
LEDGER_L1 = {("f_net", "X1"): 718.8665, ("f_net", "X2"): 1903.8603,
             ("f_region", "X1"): 551.3981, ("f_region", "X2"): 1607.3812,
             ("f_slab", "X1"): 11.7934, ("f_slab", "X2"): 13.6312,
             ("p_net", "X1"): 754.2383, ("p_net", "X2"): 1903.4177,
             ("p_region", "X1"): 583.7441, ("p_region", "X2"): 1606.9386,
             ("p_mix", "X1"): 754.2383, ("p_pin", "X2"): 13.1901}
ROWERR_L1 = {("f_net", "X1"): 0.7708365, ("f_net", "X2"): 0.5219122,
             ("f_region", "X1"): 0.5389408, ("f_region", "X2"): 0.4406373,
             ("f_slab", "X1"): 0.0093640, ("f_slab", "X2"): 0.0037368,
             ("p_net", "X1"): 0.8862167, ("p_net", "X2"): 0.6635158,
             ("p_region", "X1"): 0.6265592, ("p_region", "X2"): 0.5601641,
             ("p_mix", "X1"): 0.8862167, ("p_pin", "X2"): 0.0045709}
DPINERTIA_L1 = {("f_net", "X1"): 257.8041, ("f_net", "X2"): 642.6188,
                ("f_region", "X1"): 280.0370, ("f_slab", "X1"): 350.9757,
                ("p_net", "X1"): 271.6805, ("p_net", "X2"): 643.0614,
                ("p_region", "X1"): 298.0636}
SHARE_L1 = {"f_net": 0.319615, "f_region": 0.346431, "f_slab": 0.431868,
            "p_net": 0.335731, "p_region": 0.367420, "p_mix": 0.335731,
            "p_pin": 0.413742}
PASSES = {"f_net": 26, "f_region": 27, "f_slab": 30, "p_net": 35, "p_region": 37,
          "p_mix": 36, "p_pin": 17}
FIELDS = ("q", "dpBlock", "dpNetwork", "R", "residual", "M", "dpFriction",
          "dpInertia", "qPrev", "dq", "rowDrop", "rowError", "betaIn", "betaOut",
          "pRisePa", "dpConv")
BETA_PLATEAU = (1.0195, 1.0215)


def word(line, key):
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([A-Za-z0-9_]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return m.group(1)


def value(line, key):
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([-+0-9.eE]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return float(m.group(1))


def churchill(reynolds):
    reynolds = max(reynolds, 1e-12)
    a = (2.457 * math.log(1.0 / ((7.0 / reynolds) ** 0.9))) ** 16
    b = (37530.0 / reynolds) ** 16
    return 8.0 * ((8.0 / reynolds) ** 12 + 1.0 / (a + b) ** 1.5) ** (1.0 / 12.0)


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("transient shape x per-edge-tables solver failed")

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

    arms, cur = [], None
    for ln in lines:
        t = ln.split()
        if ln.startswith("transient start ") and " accelerate " in ln:
            cur = {"tag": ln, "arm": t[2].rstrip(":"), "inertia": {}, "shape": {},
                   "seed": {}, "partner": {}, "step": {}, "edge": {}, "pass": {},
                   "nodes": {}, "pipes": {}, "block": {}, "final": {}, "solves": 0}
            arms.append(cur)
        elif cur is None:
            continue
        if ln.startswith("SIMPLE iter 1:"):
            cur["solves"] += 1
        elif ln.startswith("transient arm ") and " inertia " in ln:
            cur["inertia"][t[4]] = t[6]
        elif ln.startswith("transient arm ") and " profile " in ln:
            cur["shape"][t[4]] = t[6]
        elif ln.startswith("transient start ") and " resistance " in ln:
            cur["seed"][word(ln, "edge")] = value(ln, "resistance")
        elif ln.startswith("transient partner "):
            cur["partner"][t[4]] = ln
        elif ln.startswith("transient step ") and " pass " in ln:
            cur["pass"].setdefault(int(t[2]), []).append(ln)
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

    by_arm = {c["arm"]: c for c in arms}

    def E(a, n, i):
        return by_arm[a]["edge"][n][NAME[i]]

    def BL(a, n):
        """The block's own report for edge name `n` in arm `a` -- one line per
        edge per arm, printed once when the arm's last level is done."""
        return by_arm[a]["block"][n]

    def nd(a, n, name):
        t = by_arm[a]["nodes"][n].split()
        return float(dict(zip(t[3::2], t[4::2]))[name])

    def brackets(a, n, i):
        ln = E(a, n, i)
        chord = value(ln, "dpNetwork") - value(ln, "dpFriction")
        ledger = value(ln, "M") * value(ln, "dq") / DT - value(ln, "dpInertia")
        return chord, ledger, value(ln, "rowError") * value(ln, "dpBlock")

    def maxdiff(a, x, b, y, field):
        return max(abs(value(E(a, n, x), field) - value(E(b, n, y), field))
                   for n in range(1, N_STEPS + 1))

    def decay(a, i):
        """The free edge's relaxation: geometric mean of the level 3..8 rate
        increments, turned into a time constant at this dt."""
        q = [value(E(a, n, i), "q") for n in range(1, N_STEPS + 1)]
        d = [q[n + 1] - q[n] for n in range(2, N_STEPS - 1)]
        ratios = [abs(d[k + 1] / d[k]) for k in range(len(d) - 1)]
        geo = math.exp(sum(math.log(x) for x in ratios) / len(ratios))
        return DT / -math.log(geo), ratios, q[-1] / Q_IN

    ALL = [(a, n, i) for a in ARMS for n in range(1, N_STEPS + 1) for i in (0, 1)]
    STEPPING = {(a, n, i) for a, n, i in ALL if EDGE_DDT[a][i] == "Euler"}
    QUASI = {(a, n, i) for a, n, i in ALL if EDGE_DDT[a][i] == "steady"}

    # ============ 1. two tables AND a shape arrive per edge
    header = [ln for ln in lines if ln.startswith("transient network+blocks")]
    require(len(header) == 1 and f"2 coupled edges, {len(ARMS)} starting sets" in header[0],
            f"one network, two coupled edges and {len(ARMS)} starting sets announced once")
    require([c["arm"] for c in arms] == list(ARMS),
            f"the seven arms ran, in the order the case declares them (got "
            f"{[c['arm'] for c in arms]})")
    require(all(c["seed"] == SEED for c in arms),
            f"every arm seeds both edges with the same pair ({SEED['X1']:.3e} then "
            f"{SEED['X2']:.3e}), so neither table nor the shape comes in through "
            "the starting coefficients")
    require(all("accelerate yes" not in c["tag"] for c in arms),
            "no arm is accelerated: every difference below is a scheme, an "
            "inertance or an inlet shape, not the update rule")
    for c in arms:
        a = c["arm"]
        require({n: word(c["partner"][n], "ddt") for n in NAME}
                == dict(zip(NAME, EDGE_DDT[a])),
                f"{a}: the rebuilt partners report "
                f"{ {n: word(c['partner'][n], 'ddt') for n in NAME} } on their own edges")
        require({n: c["inertia"][n] for n in NAME} == dict(zip(NAME, EDGE_M[a])),
                f"{a}: the row charges its two edges with "
                f"{ {n: c['inertia'][n] for n in NAME} }")
        require({n: c["shape"][n] for n in NAME} == dict(zip(NAME, EDGE_SHAPE[a])),
                f"{a}: the interface writes "
                f"{ {n: c['shape'][n] for n in NAME} } on its two inlets")
    for a in ARMS:
        for i, n in enumerate(NAME):
            require({word(E(a, k, i), "inletShape")
                     for k in range(1, N_STEPS + 1)} == {EDGE_SHAPE[a][i]},
                    f"{a}/{n}: the shape announced for the arm is the shape every "
                    f"level's edge state reports back ({EDGE_SHAPE[a][i]}), so it is "
                    "re-written every pass and not only at the first")
    require(all(len(set(EDGE_DDT[a])) == 2 and len(set(EDGE_M[a])) == 2
                and EDGE_SHAPE[a][0] != "flat" for a in ("p_mix", "p_pin")),
            "on p_mix/p_pin both tables name the SAME edge and the inlet is shaped: "
            "p_mix takes X1 as Euler+drawing and X2 as steady+slab, p_pin takes X1 as "
            "steady+drawing and X2 as Euler+slab -- no per-arm switch says either, let "
            "alone both, and neither says what the inlet looks like")
    for n in NAME:
        require(all(close(value(E(a, 1, NAME.index(n)), "M"),
                          SRC_M[EDGE_M[a][NAME.index(n)]][n], 1e-14) for a in ARMS),
                f"{n}: the installed M is exactly the number the source this arm "
                f"named announces for it -- drawing {M_DRAW[n]:.6e} (typed bore "
                f"{BORE[n]} against the equivalent {D_EXACT[n]:.8f}, so the drawing "
                f"is {M_DRAW[n] / M_REGION[n]:.4f}x the region), region "
                f"{M_REGION[n]:.6e}, slab {M_SLAB[n]:.6e}")
    require(all(close(M_SLAB[n] / M_REGION[n],
                      (PROBE[n][1] - PROBE[n][0]) / DUCT[n][3], 1e-12) for n in NAME),
              "the slab is a fraction of the region equal to its length fraction -- "
              f"{M_SLAB['X1'] / M_REGION['X1']:.6f} of X1 (0.07 m of 0.20) and "
              f"{M_SLAB['X2'] / M_REGION['X2']:.6f} of X2 (0.07 m of 0.24)")

    # ============ 2. the identity holds under both shapes, brackets and all
    worst_id = (0.0, None)
    for a, n, i in sorted(STEPPING):
        chord, ledger, lhs = brackets(a, n, i)
        ln = E(a, n, i)
        scale = (abs(value(ln, "R") * value(ln, "q")) + abs(value(ln, "dpFriction"))
                 + abs(value(ln, "M") * value(ln, "dq") / DT)
                 + abs(value(ln, "dpInertia")))
        rel = abs(lhs - (chord + ledger)) / scale
        if rel > worst_id[0]:
            worst_id = (rel, (a, n, NAME[i]))
    require(worst_id[0] < 1e-12,
            f"rowError*dpBlock == (R*q - dpFriction) + (M*dq/dt - dpInertia) in every "
            f"one of the {len(STEPPING)} stepping states of BOTH shapes, worst "
            f"{worst_id[0]:.1e} of the four terms it is made of (at {worst_id[1][0]} "
            f"level {worst_id[1][1]} {worst_id[1][2]}) -- the split 5.40 found is "
            "algebra in the printed fields, so a shaped inlet cannot break it")
    worst_part = max(abs(value(E(a, n, i), "dpFriction")
                         + value(E(a, n, i), "dpInertia")
                         - value(E(a, n, i), "dpBlock")) / value(E(a, n, i), "dpBlock")
                     for a, n, i in ALL)
    require(worst_part < 1e-13,
            f"dpFriction + dpInertia IS dpBlock in every state of every arm and shape "
            f"to {worst_part:.1e} -- the block's own partition does not care what the "
            "row was told")
    require(all(value(E(a, n, i), "rowError") == 0.0 for a, n, i in QUASI),
            f"rowError is EXACTLY 0.0 on all {len(QUASI)} quasi states and non-zero on "
            "every stepping one -- and a quasi state exists in this case with a shaped "
            "inlet AND a slab inertance (p_mix/X2), so 5.40's boundary 1 is no longer a "
            "推理")
    require(all(int(value(E(a, n, i), "unsteady")) == (EDGE_DDT[a][i] == "Euler")
                for a, n, i in ALL),
            f"all {len(ALL)} level x edge states carry the flag of their OWN edge, and "
            f"the shape and source of their OWN edge: {len(STEPPING)} stepping, "
            f"{len(QUASI)} quasi")

    # ============ 3. shape moves the block's half, not the bill
    for key in sorted(LEDGER_L1):
        a, n = key
        i = NAME.index(n)
        chord, ledger, lhs = brackets(a, 1, i)
        require(close(ledger, LEDGER_L1[key], 1e-5)
                and abs(value(E(a, 1, i), "rowError") - ROWERR_L1[key]) < 1e-6,
                f"{a}/{n} level 1: the ladder's rung is a measured number, not a trend "
                f"-- bill {ledger:.4f} Pa (recorded {LEDGER_L1[key]:.4f}), row error "
                f"{value(E(a, 1, i), 'rowError') * 100:.4f}% (recorded "
                f"{ROWERR_L1[key] * 100:.4f}%)")
    bill_fl, bill_sh = brackets("f_net", 1, 1)[1], brackets("p_net", 1, 1)[1]
    chord_fl, chord_sh = brackets("f_net", 1, 1)[0], brackets("p_net", 1, 1)[0]
    require(abs(bill_sh / bill_fl - 1.0) < 3e-4,
            f"on the PINNED edge, where dq/dt is topology-given, shaping the inlet "
            f"leaves the bill where it stands: {bill_fl:.4f} Pa flat against "
            f"{bill_sh:.4f} Pa shaped ({abs(bill_sh / bill_fl - 1) * 100:.3f}%), while "
            f"the block's own half moves (dpInertia {DPINERTIA_L1[('f_net', 'X2')]:.4f} "
            f"-> {DPINERTIA_L1[('p_net', 'X2')]:.4f} Pa) and its dpBlock falls "
            f"{abs(value(E('p_net', 1, 1), 'dpBlock') / value(E('f_net', 1, 1), 'dpBlock') - 1) * 100:.1f}%")
    require(abs(chord_sh) > abs(chord_fl)
            and abs(chord_sh) / abs(bill_sh) < 1e-4,
            f"and the chord is the bracket shape DOES touch on that edge: "
            f"{abs(chord_fl):.6f} Pa flat against {abs(chord_sh):.6f} Pa shaped, "
            f"{abs(chord_sh / chord_fl):.0f}x -- yet even shaped it is "
            f"{abs(chord_sh) / abs(bill_sh):.1e} of the bill next to it, so the two "
            "brackets still answer to different knobs")
    qfl = value(E("f_net", 1, 0), "q")
    qsh = value(E("p_net", 1, 0), "q")
    bill_fl_f, bill_sh_f = brackets("f_net", 1, 0)[1], brackets("p_net", 1, 0)[1]
    require(bill_sh_f > bill_fl_f > 0.0 and qsh > qfl > 0.0
            and abs(bill_sh_f / bill_fl_f - qsh / qfl) < 5e-3,
            f"on the FREE edge the bill moves WITH the rate it is charged with: "
            f"level-1 rate {qfl:.6e} -> {qsh:.6e} m3/s ({(qsh / qfl - 1) * 100:+.2f}%) "
            f"and bill {bill_fl_f:.4f} -> {bill_sh_f:.4f} Pa "
            f"({(bill_sh_f / bill_fl_f - 1) * 100:+.2f}%) -- the bill is (M_installed - "
            "M_the_block_paid)*dq/dt, so shape reaches it only through the trajectory, "
            "never through the inertance")
    worst_chord = max(abs(brackets(a, n, i)[0]) for a, n, i in ALL)
    l1_bills = sorted((abs(brackets(a, 1, i)[1]), a, NAME[i])
                      for a, i in [(a, i) for a in ARMS for i in (0, 1)
                                   if EDGE_DDT[a][i] == "Euler"])
    require(worst_chord < 1.0 and l1_bills[0][0] > 10.0,
            f"the two brackets stay in different hands under both shapes: the chord the "
            f"interface never iterates is at most {worst_chord:.4f} Pa anywhere in "
            f"{len(ALL)} states, while the smallest level-1 bill on a stepping edge is "
            f"{l1_bills[0][0]:.2f} Pa ({l1_bills[0][1]}/{l1_bills[0][2]}), still "
            f"{l1_bills[0][0] / worst_chord:.0f}x the worst chord, and the largest is "
            f"{l1_bills[-1][0]:.2f} Pa on {l1_bills[-1][1]}/{l1_bills[-1][2]} -- shape "
            "does not merge them, it only re-prices the block's side")

    # ============ 4. the ladder is shape-blind in the fraction it removes
    for n, i in (("X1", 0), ("X2", 1)):
        fr = {s: value(E(a, 1, i), "rowError")
              for s, a in (("network", "f_net"), ("partner", "f_region"),
                           ("slab", "f_slab")) if n == "X1" or True}
        require(fr["network"] > fr["partner"] > fr["slab"] > 0.0,
            f"flat ladder on {n}: row error falls monotonically down the source ladder "
            f"-- drawing {fr['network'] * 100:.4f}%, region {fr['partner'] * 100:.4f}%, "
            f"slab {fr['slab'] * 100:.4f}%")
        pr = {s: value(E(a, 1, i), "rowError")
              for s, a in (("network", "p_net"), ("partner", "p_region"))}
        require(pr["network"] > pr["partner"] > 0.0,
                f"shaped ladder on {n}: the same order under a parabolic inlet -- "
                f"drawing {pr['network'] * 100:.4f}%, region {pr['partner'] * 100:.4f}%")
        cut_fl = 1.0 - fr["partner"] / fr["network"]
        cut_sh = 1.0 - pr["partner"] / pr["network"]
        require(abs(cut_fl - cut_sh) < 1.5e-2,
                f"what the source removes is shape-blind: going drawing -> region cuts "
                f"{cut_fl * 100:.2f}% of {n}'s level-1 row error flat and "
                f"{cut_sh * 100:.2f}% shaped ({abs(cut_fl - cut_sh) * 100:.2f} pp apart) "
                "-- the ladder is a property of the installed M, not of what the inlet "
                "looks like")
    require(close(value(E("p_region", 1, 0), "rowError"), 0.6265592, 1e-5)
            and close(value(E("p_region", 1, 1), "rowError"), 0.5601641, 1e-5),
            "and the shaped region rung reproduces 5.39's parab_both digit for digit "
            "(62.6559% / 56.0164%, where 5.39 wrote 62.66% / 56.01% with its bores made "
            "equal) -- because 5.39's `network` number IS this case's `partner` number, "
            "so the two cases are the same arm under two different names for one source")

    # ============ 5. quasi edges pay, under slab, and grade nothing
    for a, i, n in (("p_mix", 1, "X2"), ("p_pin", 0, "X1")):
        ln = E(a, 1, i)
        want = value(ln, "M") * value(ln, "dq") / DT
        hidden = [abs(value(E(a, k, i), "M") * value(E(a, k, i), "dq") / DT)
                  for k in range(1, N_STEPS + 1)]
        require(value(ln, "rowError") == 0.0 and abs(want) > 100.0,
                f"{a}/{n}: the quasi edge under a shaped inlet carries source "
                f"{EDGE_M[a][i]} and pays {want:.2f} Pa of M*dq/dt while its block "
                f"reported dpInertia {value(ln, 'dpInertia'):.4f} Pa and the grade "
                "prints 0.0 -- 5.40's boundary 1, measured")
        require(all(value(E(a, k, i), "rowError") == 0.0
                    for k in range(1, N_STEPS + 1)),
                f"{a}/{n}: all {N_STEPS} levels of that quasi edge print 0.0 while the "
                f"installed source charges it {hidden[0]:.2f} Pa at level 1 and "
                f"{hidden[-1]:.2f} Pa at level {N_STEPS} -- on the pinned edge the "
                "hidden bill is a one-off because dq is zero after the injection step, "
                "on the free one it decays with the trajectory")
    ln = E("p_mix", 1, 1)
    want = value(ln, "M") * value(ln, "dq") / DT
    would = M_DRAW["X2"] * value(ln, "dq") / DT
    require(would > want > 0.0,
            f"p_mix/X2 is the cell 5.40's boundary 1 left unmeasured, quasi x slab: "
            f"the slab charges {want:.2f} Pa where the drawing would have charged "
            f"{would:.2f} Pa for the same dq/dt, so the source the case named removes "
            f"{(1.0 - want / would) * 100:.1f}% of a bill the grade cannot see -- "
            "choosing `slab` or `network` changes what is hidden, not what is reported")
    require(all((value(BL(a, n2), "dtLevel") == 0.0)
                == (EDGE_DDT[a][NAME.index(n2)] == "steady")
                for a in ARMS for n2 in NAME),
            "and every block hands back dtLevel 0 s exactly on the edge whose scheme row "
            "says steady and dtLevel 0.2 s on the one that says Euler -- so the shape is "
            "re-written on a quasi edge every level even though its row never sees a time "
            "layer, and the two per-edge tables reach the partner, not just the network")

    # ============ 6. the bill is linear in the installed M, slope dq/dt, both shapes
    for a, b in (("p_net", "p_region"), ("f_net", "f_region"), ("p_net", "p_pin")):
        i = 1
        want = (value(E(a, 1, i), "M") - value(E(b, 1, i), "M")) \
            * value(E(a, 1, i), "dq") / DT
        got = brackets(a, 1, i)[1] - brackets(b, 1, i)[1]
        own = (abs(brackets(a, 1, i)[0] - brackets(b, 1, i)[0])
               + abs(value(E(a, 1, i), "dpInertia") - value(E(b, 1, i), "dpInertia")))
        require(abs(got - want) / abs(want) < 1e-3 and abs(got - want) <= own + 1e-6,
                f"on the PINNED edge the difference between two attributions IS "
                f"(M_a - M_b)*dq/dt ({a} vs {b}): {got:.4f} Pa against {want:.4f} Pa, "
                f"{abs(got / want - 1) * 100:.4f}% -- and the whole residue "
                f"{abs(got - want):.4f} Pa is inside the two brackets the block itself "
                f"owns ({own:.4f} Pa). Flat and shaped alike: the slope is this level's "
                "dq/dt in both shapes, because topology gives dq to both arms")
    for a, b in (("f_net", "f_region"), ("p_net", "p_region")):
        i = 1
        d_ine = max(abs(value(E(a, k, i), "dpInertia") - value(E(b, k, i), "dpInertia"))
                    for k in range(1, N_STEPS + 1))
        r_dq = abs(value(E(a, 1, 0), "dq") / value(E(b, 1, 0), "dq") - 1.0)
        require(d_ine < 0.1 and r_dq > 0.02,
                f"attribution does not reach inside the block, but it does reach the "
                f"free edge's trajectory ({a} vs {b}): the pinned edge's own dpInertia "
                f"stays the same number under either source to {d_ine:.1e} Pa over all "
                f"{N_STEPS} levels, while the free edge's dq is NOT shared between them "
                f"({r_dq * 100:.2f}% apart at level 1) -- which is why the slope "
                "arithmetic above is exact only where topology pins the rate")

    # ============ 7. locality under a shaped inlet
    def maxdiff_upto(a, x, b, y, f, upto):
        return max(abs(value(E(a, n, x), f) - value(E(b, n, y), f))
                   for n in range(1, upto + 1))

    same5 = [f for f in FIELDS
             if maxdiff_upto("p_mix", 0, "p_net", 0, f, 5) == 0.0]
    require(len(same5) == len(FIELDS),
            f"p_mix/X1 against p_net/X1 -- same scheme, same source, same shape on THAT "
            f"edge, the OTHER edge changed both its tables -- agree bit for bit in "
            f"{len(same5)} of {len(FIELDS)} fields at levels 1..5 and separate only at "
            f"the tail: {maxdiff('p_mix', 0, 'p_net', 0, 'rowError'):.1e} in row error "
            f"and {maxdiff('p_mix', 0, 'p_net', 0, 'dpBlock'):.1e} Pa in dpBlock at "
            f"levels 6..8. p_net's pinned edge first fills its sweep budget at level 5, "
            "so the shared prefix ends exactly where a neighbour's block gives up -- "
            "measured coincidence, not a derived trigger")
    _, rat_a, dest_a = decay("p_net", 0)
    _, rat_b, dest_b = decay("p_mix", 0)
    require(abs(decay("p_net", 0)[0] - decay("p_mix", 0)[0]) < 5e-3
            and abs(dest_a - dest_b) < 5e-4,
            f"and the tail separation is small: free-edge tau {decay('p_net', 0)[0]:.4f} "
            f"against {decay('p_mix', 0)[0]:.4f} s, destinations {dest_a:.6f} against "
            f"{dest_b:.6f}, decay ratios {rat_a} against {rat_b} (worst "
            f"{max(abs(x - y) for x, y in zip(rat_a, rat_b)):.1e}) -- 5.40's asymmetry "
            "(formulation commutes across the two edges, attribution does not, because "
            "the installed M is a term of the global solve) survives the shaped inlet "
            "but does not dominate it")
    qd = maxdiff("p_pin", 0, "p_net", 0, "q")
    require(qd > 1e-6,
            f"an edge that changes its OWN scheme row is not marginal at all: p_pin/X1 "
            f"(quasi + drawing + shaped) against p_net/X1 (Euler + drawing + shaped), "
            f"with M identical to the last bit, moves the rate by {qd:.3e} m3/s and the "
            f"split from {value(E('p_net', 1, 0), 'q') / Q_IN:.6f} to "
            f"{value(E('p_pin', 1, 0), 'q') / Q_IN:.6f} at level 1 -- and the grade that "
            f"would have reported the difference, {value(E('p_net', 1, 0), 'rowError') * 100:.4f}% "
            "on the stepping edge, is exactly 0.0 on the quasi one")
    for a, b, i, n in (("f_net", "p_net", 1, "X2"), ("f_region", "p_region", 1, "X2")):
        w = max(abs(value(E(a, k, i), "dpBlock") - value(E(b, k, i), "dpBlock"))
                / value(E(b, k, i), "dpBlock") for k in range(1, N_STEPS + 1))
        require(w > 0.05,
                f"shape does reach the pinned edge's ledger through the block: {a} "
                f"against {b} moves dpBlock by up to {w * 100:.1f}% while leaving the "
                f"bill at {abs(brackets(b, 1, i)[1] - brackets(a, 1, i)[1]):.2f} Pa -- "
                "the row's number and the region's number move in different ledgers")

    # ============ 8. the 300-sweep wall is shape's, priced by attribution
    wall = {a: {n: [k for k in range(1, N_STEPS + 1)
                    if value(E(a, k, NAME.index(n)), "sweeps") >= MAX_SWEEPS]
                for n in NAME} for a in ARMS}
    require(all(not any(wall[a][n] for n in NAME) for a in FLAT),
            "no flat state touches the 300-sweep budget -- the worst flat edge anywhere "
            f"is {max(int(value(E(a, k, i), 'sweeps')) for a in FLAT for k in range(1, N_STEPS + 1) for i in (0, 1))} sweeps")
    at_wall = [(a, n, k) for a in ARMS for n in NAME for k in wall[a][n]]
    require(at_wall and all(a in SHAPED for a, n, k in at_wall),
            "and shaped arms do, flat ones never: " + ", ".join(
                f"{a}/X2 levels {wall[a]['X2']}" for a in ARMS if wall[a]["X2"]) +
            " -- the wall belongs to the shape, as 5.39 recorded, and in this crossing "
            "it sits on the PINNED edge")
    sw_net = max(int(value(E("p_net", k, 1), "sweeps")) for k in range(1, N_STEPS + 1))
    sw_reg = max(int(value(E("p_region", k, 1), "sweeps")) for k in range(1, N_STEPS + 1))
    require(sw_net >= MAX_SWEEPS and sw_reg < sw_net,
            f"attribution prices the wall: the same shaped pinned edge asks "
            f"{sw_net} sweeps under the drawing and peaks at {sw_reg} under the region, "
            f"which never reaches the budget -- shape puts the edge near the budget, the "
            "installed inertance decides whether it reaches it")
    nc = {(a, n, k) for a in ARMS for n in NAME
          for k in range(1, N_STEPS + 1)
          if int(value(E(a, k, NAME.index(n)), "blockConverged")) == 0}
    res_at_wall = max(value(E(a, k, NAME.index(n)), "residual") for a, n, k in at_wall)
    require(nc == {(a, n, k) for a, n, k in at_wall} and res_at_wall < TOL_INTERFACE,
            "and here, unlike 5.39, the wall and the convergence report are the SAME "
            f"{len(nc)} states: " + ", ".join(f"{a}/{n} level {k}"
                                              for a, n, k in sorted(nc)) +
            f" -- every state that fills the budget reports blockConverged 0 and no "
            f"non-converged state stops short of it, yet all of them are still accepted "
            f"(worst interface residual there {res_at_wall:.3e}, under the 1e-3 "
            "tolerance) and their level completes. So 'filled the budget' stays a record "
            "and not a stop rule the interface enforces, and the coincidence is this "
            "crossing's, not 5.39's six-against-five")

    # ============ 9. cost, and what it is a function of
    for a in ARMS:
        got = sum(len(v) for v in by_arm[a]["pass"].values())
        require(got == PASSES[a],
                f"{a}: {got} interface sub-iterations over the eight levels (measured "
                f"against {PASSES[a]} declared for this crossing)")
        require(by_arm[a]["solves"] == 2 * PASSES[a],
                f"{a}: {by_arm[a]['solves']} block solves, exactly twice the passes")
        worst = max(value(E(a, k, i), "residual")
                    for k in range(1, N_STEPS + 1) for i in (0, 1))
        require(worst < TOL_INTERFACE,
                f"{a}: worst interface residual {worst:.3e} under the 1e-3 tolerance")
    tot = sum(PASSES.values())
    require(all(PASSES[a] < PASSES[b] for a, b in (("f_net", "p_net"),
                                                   ("f_region", "p_region"),
                                                   ("p_net", "p_region"))),
            f"the shaped rungs cost more negotiation than their flat twins at the SAME "
            "attribution (f_net 26 -> p_net 35, f_region 27 -> p_region 37 passes), so "
            f"of the two knobs the shape is the expensive one; the seven arms walk "
            f"{tot} sub-iterations and {2 * tot} block solves")

    # ============ 10. the ramp: shape doubles the relaxation, source sets its order
    for a in ARMS:
        require(close(value(E(a, 1, 0), "q") / Q_IN, SHARE_L1[a], 2e-5),
                f"{a}: the split arrives at the parallel branch as "
                f"{value(E(a, 1, 0), 'q') / Q_IN:.6f} -- the tables AND the shape for "
                "THAT edge are what the branch sees")
    taus = {}
    for a in ARMS:
        taus[a] = decay(a, 0)[0]
    require(taus["f_net"] > taus["f_region"] > taus["f_slab"],
            f"flat ladder: tau {taus['f_net']:.4f} / {taus['f_region']:.4f} / "
            f"{taus['f_slab']:.4f} s -- a row that resists acceleration more relaxes "
            "slower, whatever the inlet")
    require(taus["p_net"] > taus["p_region"]
            and taus["p_net"] / taus["f_net"] > 1.5
            and taus["p_region"] / taus["f_region"] > 1.5,
            f"shaped: tau {taus['p_net']:.4f} / {taus['p_region']:.4f} s against "
            f"{taus['f_net']:.4f} / {taus['f_region']:.4f} flat -- the same source, and "
            f"the shaped inlet roughly doubles the relaxation time "
            f"(x{taus['p_net'] / taus['f_net']:.2f}, x{taus['p_region'] / taus['f_region']:.2f}), "
            "so the ramp is shaped by the entry profile and ordered by the inertance")
    require(SHARE_L1["f_net"] < SHARE_L1["f_region"] < SHARE_L1["f_slab"]
            and SHARE_L1["p_net"] < SHARE_L1["p_region"],
            f"and the first step still goes more through the edge whose row resists "
            f"acceleration less: flat {SHARE_L1['f_net']:.6f} / "
            f"{SHARE_L1['f_region']:.6f} / {SHARE_L1['f_slab']:.6f}, shaped "
            f"{SHARE_L1['p_net']:.6f} / {SHARE_L1['p_region']:.6f}")

    # ============ 11. the partner is still a partner, in both shapes
    for a in ARMS:
        for n in NAME:
            bl = BL(a, n)
            yp = value(bl, "wallYPlusMax")
            nu = value(bl, "nutMax") / NU
            re_d = value(bl, "fluxIn") / (AREA[n] * NU) * D_EXACT[n]
            bo = value(bl, "betaOut")
            imb = abs(value(bl, "fluxIn") - value(bl, "fluxOut")) / value(bl, "fluxIn")
            require(11.0 < yp < 40.0 and nu > 40.0 and 3e3 < re_d < 2e4
                    and BETA_PLATEAU[0] < bo < BETA_PLATEAU[1] and imb < 1e-9,
                    f"{a}/{n}: first cell at y+ {yp:.1f} (log law needs > 11), nut/nu "
                    f"{nu:.1f}, Re_Dh {re_d:.0f}, beta_out {bo:.6f} in the plateau band, "
                    f"mass imbalance {imb:.1e}, shape {word(bl, 'inletShape')}")
    b_in = {a: {n: value(BL(a, n), "betaIn") for n in NAME} for a in ARMS}
    require(all(b_in[a][n] < 1.0015 for a in FLAT for n in NAME)
            and all(b_in[a][n] > 1.09 for a in SHAPED for n in NAME),
            "the declared shape is what the inlet actually asks for: beta_in "
            f"{b_in['f_net']['X1']:.6f} / {b_in['f_region']['X1']:.6f} / "
            f"{b_in['f_slab']['X1']:.6f} on the flat arms and "
            f"{b_in['p_net']['X1']:.6f} / {b_in['p_net']['X2']:.6f} on the shaped ones "
            "(Poiseuille's 1.18172 relaxed toward the wall-function duct's own "
            "development), so nothing above is a shape effect in disguise -- and every "
            "beta_out is in the plateau band in both shapes")
    LEAD_DROP = 11.2641281367
    for a in ARMS:
        require(all(close(value(E(a, k, i), "rowDrop"),
                          nd(a, k, "J1") - nd(a, k, "J2") if i == 0 else
                          nd(a, k, "J2") - nd(a, k, "J3"), 1e-9)
                    for k in range(2, N_STEPS + 1) for i in (0, 1)),
                f"{a}: every level's node pair equals that edge's ROW drop, inertial "
                "term included, under either inlet shape -- the source the case named "
                "enters the pressure field the network adopts")
        require(close(abs(nd(a, N_STEPS, "Src") - nd(a, N_STEPS, "J1"))
                      / LEAD_DROP, 1.0, 1e-9)
                and close(abs(nd(a, N_STEPS, "J3") - nd(a, N_STEPS, "Snk"))
                          / LEAD_DROP, 1.0, 1e-9)
                and nd(a, N_STEPS, "Snk") == P_REF,
                f"{a}: both lead sections still drop the network's own Churchill "
                f"{LEAD_DROP:.10f} Pa and Snk sits on the reference pressure, so "
                "neither table nor shape reached the 1D side")
        q1 = value(E(a, N_STEPS, 0), "q")
        pair = abs(nd(a, N_STEPS, "J1") - nd(a, N_STEPS, "J2"))
        require(0.0 < q1 < Q_IN and pair > 0.0,
                f"{a}: the parallel branch closes -- X1 carries {q1:.6e} m3/s at the "
                f"last level of the injected {Q_IN:.1e}, and the node pair it shares "
                f"with the bypass holds {pair:.4f} Pa")
    # ============ 12. the account the row still cannot state: what the region's
    # own changing profile is worth, level by level, edge by edge
    dev = 0.0
    for a in ARMS:
        for n in NAME:
            i = NAME.index(n)
            for k in range(1, N_STEPS + 1):
                ln = E(a, k, i)
                q = value(ln, "q")
                aa = value(ln, "inletArea")
                want = RHO * (value(ln, "betaOut") - value(ln, "betaIn")) \
                    * abs(q) * q / (aa * aa)
                got = value(ln, "dpConv")
                if want != 0.0:
                    dev = max(dev, abs(got / want - 1.0))
                else:
                    dev = max(dev, abs(got - want))
    require(dev < 1e-13,
            f"every one of the 112 (level x edge) states reports dpConv == "
            f"rho*(betaOut - betaIn)*q|q|/A^2 built from its OWN beta pair, rate and "
            f"section area -- worst relative deviation {dev:.2e}. The number is the "
            "difference of two momentum fluxes at the same rate through the same area, "
            "so it is the shape alone, and the interface neither asks for it nor pays it")

    pos = [(a, n, k) for a in FLAT for n in NAME for k in range(1, N_STEPS + 1)
           if value(E(a, k, NAME.index(n)), "dpConv") > 0.0]
    neg = [(a, n, k) for a in SHAPED for n in NAME for k in range(1, N_STEPS + 1)
           if value(E(a, k, NAME.index(n)), "dpConv") < 0.0]
    require(len(pos) == len(FLAT) * 2 * N_STEPS and len(neg) == len(SHAPED) * 2 * N_STEPS,
            f"the sign follows the shape on every state: all {len(pos)} flat states are "
            f"POSITIVE (the region had to grow a fuller profile than the plug it was "
            f"handed, beta_out {value(E('f_net', 8, 0), 'betaOut'):.4f} above beta_in "
            f"{value(E('f_net', 8, 0), 'betaIn'):.4f}) and all {len(neg)} shaped states "
            f"NEGATIVE (it sheds the Poiseuille excess, "
            f"{value(E('p_net', 8, 0), 'betaIn'):.4f} down to "
            f"{value(E('p_net', 8, 0), 'betaOut'):.4f}) -- 112 states, no exception")

    fr = {(a, n, k): abs(value(E(a, k, NAME.index(n)), "dpConv"))
          / value(E(a, k, NAME.index(n)), "dpBlock")
          for a in ARMS for n in NAME for k in range(1, N_STEPS + 1)}
    flat_max = max((((a, n, k), v) for (a, n, k), v in fr.items() if a in FLAT),
                   key=lambda item: item[1])
    (wa, wn, wk), wv = flat_max
    require(flat_max[1] < 0.06,
            f"on a plug inlet the account is small: worst flat share of the block's own "
            f"drop {wv * 100:.2f}% ({wa}/{wn} level {wk}, "
            f"{value(E(wa, wk, NAME.index(wn)), 'dpConv'):.3f} Pa of "
            f"{value(E(wa, wk, NAME.index(wn)), 'dpBlock'):.2f} Pa)")
    pin8 = fr[("p_pin", "X1", 8)]
    require(pin8 > 1.0 and fr[("p_net", "X2", 8)] > 0.39,
            f"shaped, the same account outruns the drop it lives inside: p_pin/X1 level 8 "
            f"charges {value(E('p_pin', 8, 0), 'dpConv'):.3f} Pa against a block reporting "
            f"{value(E('p_pin', 8, 0), 'dpBlock'):.2f} Pa -- {pin8 * 100:.1f}% of the whole "
            f"drop the interface is told to reproduce -- while p_net/X2 at the plateau is "
            f"{fr[('p_net', 'X2', 8)] * 100:.1f}%. The row's friction law has no term that "
            "answers for it, in either direction")

    per_q = {(a, n, k): value(E(a, k, NAME.index(n)), "dpConv")
             / value(E(a, k, NAME.index(n)), "q") ** 2
             for a in ARMS for n in NAME for k in range(1, N_STEPS + 1)}
    pin_spread = max(max(per_q[(a, "X2", k)] for a in FLAT)
                     / min(per_q[(a, "X2", k)] for a in FLAT) - 1.0
                     for k in range(1, N_STEPS + 1))
    m_ratio = M_DRAW["X1"] / M_SLAB["X1"]
    pin8 = " / ".join(f"{per_q[(a, 'X2', 8)]:.4e}" for a in FLAT)
    dpin_l1 = " / ".join(f"{value(E(a, 1, 1), 'dpInertia'):.4f}" for a in FLAT)
    require(pin_spread < 1e-5,
            f"on the pinned edge the account is the same number under all three tables: "
            f"across the flat arms dpConv/q^2 agrees to {pin_spread:.1e} relative over all "
            f"eight levels (X2 level 8: {pin8} Pa.s2/m6) -- and so does the block's own "
            f"dpInertia ({dpin_l1} Pa at level 1) -- although the inertances that table "
            f"installs on that edge span {m_ratio:.2f}x")
    coeff_spread = [max(per_q[(a, "X1", k)] for a in FLAT)
                    / min(per_q[(a, "X1", k)] for a in FLAT) - 1.0
                    for k in range(1, N_STEPS + 1)]
    rate_spread = [max(value(E(a, k, 0), "q") for a in FLAT)
                   / min(value(E(a, k, 0), "q") for a in FLAT) - 1.0
                   for k in range(1, N_STEPS + 1)]
    acct_spread = [max(value(E(a, k, 0), "dpConv") for a in FLAT)
                   / min(value(E(a, k, 0), "dpConv") for a in FLAT) - 1.0
                   for k in range(1, N_STEPS + 1)]
    require(max(coeff_spread) < 0.08
            and all(coeff_spread[k] < rate_spread[k] for k in range(N_STEPS))
            and all(coeff_spread[k + 1] < coeff_spread[k] * 1.02 for k in range(N_STEPS - 1))
            and coeff_spread[-1] < coeff_spread[0] / 8,
            f"on the free edge the SAME coefficient is table-blind too, and the trajectory "
            f"is the only way in: dpConv itself spans {acct_spread[0] * 100:.1f}% at level 1 "
            f"(7.632 / 9.157 / 15.044 Pa), but the rate itself spans {rate_spread[0] * 100:.1f}% "
            f"there, and after dividing by q^2 the residue is {coeff_spread[0] * 100:.2f}% "
            f"-- smaller than the rate spread at every one of the eight levels, and itself "
            f"collapsing monotonically {coeff_spread[0] * 100:.2f}% -> "
            f"{coeff_spread[-1] * 100:.2f}% as the three arms' trajectories converge")
    ratio_flat_shaped = abs(per_q[("p_net", "X2", 8)]) / per_q[("f_net", "X2", 8)]
    require(ratio_flat_shaped > 3.5,
            f"while the SAME ratio across shapes is {ratio_flat_shaped:.2f}x "
            f"({per_q[('f_net', 'X2', 8)]:.4e} flat against {abs(per_q[('p_net', 'X2', 8)]):.4e} "
            "shaped Pa.s2/m6) -- so the knob that moves this account is the inlet profile, "
            "the one thing the case added")

    rise = [(a, n, k, value(E(a, k, NAME.index(n)), "dpConv"),
             value(E(a, k, NAME.index(n)), "pRisePa"))
            for a in SHAPED for n in NAME for k in range(2, N_STEPS + 1)
            if value(E(a, k, NAME.index(n)), "pRisePa") > 50.0]
    rr = [abs(c / pr) for _, _, _, c, pr in rise]
    plateau = [abs(value(E(a, k, 1), "dpConv") / value(E(a, k, 1), "pRisePa"))
               for a in SHAPED for k in (4, 5, 6, 7, 8)]
    l1 = abs(value(E("p_net", 1, 1), "dpConv") / value(E("p_net", 1, 1), "pRisePa"))
    require(len(rise) > 20 and min(rr) > 0.8 and max(rr) < 1.35,
            f"where the region really hands pressure BACK, the account is the same order "
            f"as the handing but never equal to it: {len(rise)} shaped states from level 2 "
            f"on with a rise above 50 Pa give |dpConv|/pRise in [{min(rr):.4f}, "
            f"{max(rr):.4f}], and on the X2 plateau (levels 4-8, all four shaped arms) it "
            f"holds to [{min(plateau):.5f}, {max(plateau):.5f}] -- p_net/X2 "
            f"{value(E('p_net', 8, 1), 'dpConv'):.2f} Pa against a rise of "
            f"{value(E('p_net', 8, 1), 'pRisePa'):.2f} Pa. Level 1 sits outside the band at "
            f"{l1:.4f} (the rise has not developed yet, "
            f"{value(E('p_net', 1, 1), 'pRisePa'):.2f} Pa against "
            f"{value(E('p_net', 1, 1), 'dpConv'):.2f} Pa), and the flat states report "
            "pRise exactly 0.0 while still owing a positive account "
            f"({value(E('f_net', 8, 1), 'dpConv'):.3f} Pa), "
            "so the rise is one symptom of the imbalance rather than its amount")

    flat_rise = {value(E(a, k, NAME.index(n)), "pRisePa") for a in FLAT
                 for n in NAME for k in range(1, N_STEPS + 1)}
    require(flat_rise == {0.0},
            f"every flat state's measured rise is exactly 0.0 Pa ({len(flat_rise)} value"
            f"{'s' if len(flat_rise) != 1 else ''} in the set) while its convective account "
            f"is positive on all {len(pos)} of them -- the region owes the imbalance even "
            "when its station pressures never go up")

    gaps = []
    zeroed = []
    for a in ARMS:
        for n in NAME:
            i = NAME.index(n)
            for k in range(1, N_STEPS + 1):
                ln = E(a, k, i)
                if int(value(ln, "unsteady")) == 0:
                    continue
                gap = value(ln, "rowDrop") - value(ln, "dpBlock")
                c = value(ln, "dpConv")
                if abs(gap) > 100.0:
                    gaps.append((a, n, k, c / gap, gap, c))
                if abs(gap) < 0.01 and abs(c) > 50.0:
                    zeroed.append((a, n, k, gap, c))
    gr = [g[3] for g in gaps]
    zf = [z for z in zeroed if z[0] in FLAT]
    zs = [z for z in zeroed if z[0] not in FLAT]
    require(len(gaps) >= 8 and max(abs(r) for r in gr) < 0.2
            and all((r < 0.0) == (g[0] in SHAPED) for g, r in zip(gaps, gr))
            and len(zeroed) >= 8,
            f"the row's disagreement is not this account either: over the {len(gaps)} "
            f"stepping states whose bill exceeds 100 Pa -- all of them level 1 -- dpConv is "
            f"at most {max(abs(r) for r in gr) * 100:.1f}% of it (f_net/X1 level 1: "
            f"{value(E('f_net', 1, 0), 'dpConv'):.3f} Pa against a "
            f"{value(E('f_net', 1, 0), 'rowDrop') - value(E('f_net', 1, 0), 'dpBlock'):.2f} Pa "
            f"gap), and on the shaped arms it carries the OPPOSITE sign to the bill, so "
            "charging the row the convective term would leave at least four fifths of the "
            "5.38/5.40 double-billing and on those edges enlarge it. Meanwhile "
            f"{len(zeroed)} stepping states report a bill that has already vanished "
            f"(|rowDrop - dpBlock| < 0.01 Pa) while still owing a real account: {len(zf)} "
            f"of them are the flat X2 plateau at {zf[-1][4]:+.2f} Pa and {len(zs)} shaped "
            f"({', '.join(f'{a}/{n} level {k} at {c:+.2f} Pa' for a, n, k, g, c in zs)})"
            ": two different quantities. This is "
            "what the row still does not know -- the number is on every line and in no "
            "equation")
    require(all(int(value(E(a, k, NAME.index(n)), "unsteady")) == 0
                or abs(value(E(a, k, NAME.index(n)), "rowDrop")
                        - (value(E(a, k, NAME.index(n)), "R")
                           * value(E(a, k, NAME.index(n)), "q")
                           + value(E(a, k, NAME.index(n)), "M")
                           * value(E(a, k, NAME.index(n)), "dq") / DT)) < 1e-6
                for a in ARMS for n in NAME for k in range(1, N_STEPS + 1)),
                "the row as written is R*q + M*dq/dt on every stepping state -- dpConv is "
                "reported beside it, not inside it, so the 1D momentum equation is still "
                "shape-blind and this increment leaves it that way")

    print(f"\n{'PASS' if not fails else 'FAIL'}: {len(fails)} failure(s)")
    if fails:
        for f in fails:
            print("  - " + f)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
