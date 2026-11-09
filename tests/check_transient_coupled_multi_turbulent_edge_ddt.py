#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient-multi over two wall-function RANS ducts whose
time formulation is chosen EDGE by EDGE, so one time level holds a partner that
answers du/dt next to a partner that answers a chord.

The topology is the multi-edge one the profile chain has used since §5.34:
Src -> A -> { X1 (3D RANS duct) || By (1D bypass) } -> X2 (3D RANS duct) -> B -> Snk,
with a constant 6.0e-5 m3/s switched on at t = 0 so every rate starts at zero and
ramps over eight levels of dt = 0.2 s. Everything physical is §5.38's, unchanged:
same meshes, same bores, same wall functions, same operating point, same dt, same
seeds, same tolerance. What is new is that `couple{ ddtFrom }` names a scheme per
coupled edge inside an arm, which the per-arm `ddts` table cannot express:

  lo_steady    both edges quasi-steady               (5.38's control)
  lo_both      both edges time-accurate              (5.38's crossing)
  lo_x1free    the edge whose rate RELAXES steps, the pinned one does not
  lo_x2pin     the mirror: the pinned edge steps, the free one does not
and the shaped family repeats three of those rows with parabolic inlets, so the
mixed arm's edges can be judged against controls that agree with them on BOTH
that edge's scheme and that edge's shape:

  parab_steady / parab_both / parab_x1free

The kernel already took `chord` block by block (coupledNetwork.cpp:1182) and already
knew each block's own formulation (coupledNetwork.cpp:253), so nothing about the
interface had to change: this increment is the first one that can ASK for a mixed
level. Four things could only have been measured here:

  1. the declaration reaches the edge, not the arm -- one arm prints `ddt Euler` and
     `ddt steady` on its two partners, and the per-edge `unsteady` flag follows it at
     all 112 level x edge states
  2. an untouched edge IS the control's computation -- each mixed arm's edge is
     bit-for-bit the arm that matches it on that edge's scheme and shape in every
     field (rate, both halves of the drop, coefficient, residual, sweeps, row error,
     beta, recovery, flag: 0.0 over 8 levels), and the pinned edge it left
     quasi-steady differs from its own control by less than a pascal
  3. what does not compose is the attribution -- the double count 5.38 measured
     appears on exactly the edge that steps and is identically 0.0 on every
     quasi-steady one, in the same arm and the same level, while the row keeps its
     own inertance on both
  4. the ramp belongs to the edge whose rate is free -- first-level split,
     destination, time constant, exchange schedule and the block sweep counts of a
     mixed arm are the ones of whichever control has the same scheme AND the same
     shape on the PARALLEL edge; stepping the pinned edge buys nothing there

usage: check_transient_coupled_multi_turbulent_edge_ddt.py <solver> <case>
"""
import math
import os
import re
import subprocess
import sys

# fluid, geometry and schedule of
# tests/cases/transientCoupledMultiTurbulentEdgeDdt/case.dict
MU, RHO = 1.0e-3, 1000.0
NU = MU / RHO
Q_IN = 6.0e-5                  # m3/s injected at Src from t = 0
P_REF = 1.0e6                  # Pa held at Snk
DT, N_STEPS = 0.2, 8           # time { dt 0.2; endTime 1.6; }
MAX_SWEEPS = 300.0             # simple{ maxIterations }, each block's own budget
TOL_INTERFACE = 1e-3           # couple{ tolerance }

ARMS = ("lo_steady", "lo_both", "lo_x1free", "lo_x2pin",
        "parab_steady", "parab_both", "parab_x1free")
NAME = ("X1", "X2")
FREE, PINNED = 0, 1            # which index is the relaxing edge and which is pinned

# the per-edge scheme the case's `ddtFrom` table asks for -- the fact this whole
# checker is organised around, since a per-arm switch cannot express rows 3, 4 and 7
EDGE_DDT = {"lo_steady": ("steady", "steady"),
            "lo_both": ("Euler", "Euler"),
            "lo_x1free": ("Euler", "steady"),
            "lo_x2pin": ("steady", "Euler"),
            "parab_steady": ("steady", "steady"),
            "parab_both": ("Euler", "Euler"),
            "parab_x1free": ("Euler", "steady")}
SHAPE = {"lo_steady": ("flat", "flat"),
         "lo_both": ("flat", "flat"),
         "lo_x1free": ("flat", "flat"),
         "lo_x2pin": ("flat", "flat"),
         "parab_steady": ("parabolic", "parabolic"),
         "parab_both": ("parabolic", "parabolic"),
         "parab_x1free": ("parabolic", "parabolic")}
MIXED = tuple(a for a in ARMS if len(set(EDGE_DDT[a])) > 1)   # the three mixed arms
SHAPED = tuple(a for a in ARMS if SHAPE[a] == ("parabolic", "parabolic"))
# for each mixed arm and each of its edges, the arm that agrees with it on THAT
# edge's scheme AND on that edge's shape -- the control a mixed arm's edge is
# supposed to reproduce exactly, and the reason the shaped family is in the case
CONTROL = {"lo_x1free": {FREE: "lo_both", PINNED: "lo_steady"},
           "lo_x2pin": {FREE: "lo_steady", PINNED: "lo_both"},
           "parab_x1free": {FREE: "parab_both", PINNED: "parab_steady"}}
SEED = {"X1": 1.4e7, "X2": 1.36e7}

# gap, span, streamwise cells, region length. The wall-bounded axis is y in both
# blocks and the z faces are symmetry planes, so a section is gap x span.
DUCT = {"X1": (0.005, 0.005, 20, 0.2), "X2": (0.004, 0.008, 24, 0.24)}
PANELS = 8                     # wall-normal panels the inlet patch is made of
BY = (0.3, 0.008)               # the 1D bypass, sharing X1's node pair
LEAD = (0.3, 0.02)              # A and B
AREA = {n: DUCT[n][0] * DUCT[n][1] for n in NAME}
PROBE_INERTANCE = {"X1": 2.8e6, "X2": 2.1875e6}   # rho * 0.07 / A, the slab announced

CONTINUUM = {"flat": 1.0, "parabolic": 6.0 / 5.0,
             "power1over7": (7.0 / 9.0) / (7.0 / 8.0) ** 2}
# the plateau this wall-function channel develops to
BETA_PLATEAU = (1.0200, 1.0211)


def word(line, key):
    """The non-numeric value stored under `key` ('inletShape flat ...')."""
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


def weight(shape, eta):
    """The declared profile at a section coordinate eta in [-1, 1]."""
    if shape == "flat":
        return 1.0
    if shape == "parabolic":
        return 1.0 - eta * eta
    if shape == "power1over7":
        return (1.0 - abs(eta)) ** (1.0 / 7.0)
    raise SystemExit(f"unknown shape '{shape}'")


def betaDeclared(shape, panels):
    """<u^2>/<u>^2 of the declared shape as the mesh can see it, on `panels`
    equal panels across the wall-bounded gap."""
    w = [weight(shape, (2.0 * (j + 0.5) / panels - 1.0)) for j in range(panels)]
    mean = sum(w) / len(w)
    return sum(x * x for x in w) / len(w) / (mean * mean)


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


def inertance(length, area):
    return RHO * length / area


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("transient per-edge ddt solver failed")

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

    def same(a_, b_, x_, y_, field):
        """The worst difference of one field between arm `a_`'s edge and arm `b_`'s,
        over all eight levels."""
        return max(abs(value(by_arm[a_]["edge"][n][x_], field)
                       - value(by_arm[b_]["edge"][n][y_], field))
                  for n in range(1, N_STEPS + 1))

    # ---- split the log into the arms ----
    arms = []
    cur = None
    for ln in lines:
        t = ln.split()
        if ln.startswith("transient start ") and " accelerate " in ln:
            cur = {"tag": ln, "arm": t[2].rstrip(":"), "shape": {}, "seed": {},
                   "partner": {}, "step": {}, "edge": {}, "pass": {}, "nodes": {},
                   "pipes": {}, "block": {}, "final": {}, "summary": "", "solves": 0}
            arms.append(cur)
        elif cur is None:
            continue
        if ln.startswith("SIMPLE iter 1:"):
            cur["solves"] += 1
        elif ln.startswith("transient arm ") and " profile " in ln:
            cur["shape"][t[4]] = t[6]
        elif ln.startswith("transient start ") and " resistance " in ln:
            cur["seed"][word(ln, "edge")] = value(ln, "resistance")
        elif ln.startswith("transient partner "):
            cur["partner"][t[4]] = ln
        elif ln.startswith("transient step ") and " pass " in ln:
            cur["pass"].setdefault(int(t[2]), {})[int(t[4])] = ln
        elif ln.startswith("transient step ") and " edge " in ln:
            cur["edge"].setdefault(int(t[2]), {})[int(t[4])] = ln
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
        return by_arm[a]["edge"][n][i]

    def B(a, n):
        return by_arm[a]["block"][n]

    def nd(a, n, name):
        t = by_arm[a]["nodes"][n].split()
        return float(dict(zip(t[3::2], t[4::2]))[name])

    def pipeq(a, n, name):
        t = by_arm[a]["pipes"][n].split()
        return float(dict(zip(t[3::2], t[4::2]))[name])

    ALL = [(a, n, i) for a in ARMS for n in range(1, N_STEPS + 1) for i in (0, 1)]
    STEPPING = {(a, n, i) for a, n, i in ALL if EDGE_DDT[a][i] == "Euler"}
    QUASI = {(a, n, i) for a, n, i in ALL if EDGE_DDT[a][i] == "steady"}

    # ------------- 1. the declaration reaches the EDGE, and it beats the arm
    header = [ln for ln in lines if ln.startswith("transient network+blocks")]
    require(len(header) == 1 and f"2 coupled edges, {len(ARMS)} starting sets" in header[0],
            f"one network, two coupled edges and {len(ARMS)} starting sets announced once")
    require([c["arm"] for c in arms] == list(ARMS),
            f"the seven arms ran, in the order the case declares them (got "
            f"{[c['arm'] for c in arms]})")
    require(all(c["seed"] == SEED for c in arms),
            f"every arm seeds both edges with the same pair ({SEED['X1']:.3e} then "
            f"{SEED['X2']:.3e}), so no arm starts at its own shaped resistance")
    require({c["arm"]: "accelerate yes" in c["tag"] for c in arms}
            == {a: False for a in ARMS},
            "no arm is accelerated: every difference below is the scheme on one edge "
            "or the shape, not the update rule")
    for c in arms:
        require(c["shape"] == dict(zip(NAME, SHAPE[c["arm"]])),
                f"{c['arm']}: the `profile` reached its two edges as declared "
                f"({c['shape']['X1']} on both)")
        got = {n: word(c["partner"][n], "ddt") for n in NAME}
        require(got == dict(zip(NAME, EDGE_DDT[c["arm"]])),
                f"{c['arm']}: the rebuilt partners report {got} -- one scheme per "
                "edge, which no per-arm switch can say")
    mixed_words = {a: {n: word(by_arm[a]["partner"][n], "ddt") for n in NAME}
                   for a in MIXED}
    require(all(len(set(v.values())) == 2 for v in mixed_words.values()),
            "and on all three mixed arms the two words are DIFFERENT "
            + " / ".join(f"{a}: {v['X1']} on the free edge, {v['X2']} on the pinned "
                         "one" for a, v in mixed_words.items())
            + " -- the case's per-arm word is `Euler` on arms 3, 4 and 7, so the table "
            "beat it on one edge and let it stand on the other")
    require(all(word(by_arm["lo_steady"]["partner"][n], "ddt") == "steady"
                for n in NAME),
            "the control's per-arm word `steady` beat the blocks' own "
            "`schemes{ ddt Euler; }` on both edges, so the override is exercised in "
            "the other direction too")
    flag = {(a, n, i): value(E(a, n, i), "unsteady") for a, n, i in ALL}
    require(flag == {(a, n, i): (1.0 if EDGE_DDT[a][i] == "Euler" else 0.0)
                     for a, n, i in ALL},
            f"all {len(ALL)} level x edge states carry the flag of their OWN edge: "
            f"{len(STEPPING)} say the region answered du/dt at that level and "
            f"{len(QUASI)} say it did not, and on the two mixed arms both kinds sit "
            "in the same level")
    require(all(value(B(a, n), "unsteady")
                == (1.0 if EDGE_DDT[a][NAME.index(n)] == "Euler" else 0.0)
                and value(B(a, n), "dtLevel")
                == (DT if EDGE_DDT[a][NAME.index(n)] == "Euler" else 0.0)
                for a in ARMS for n in NAME),
            f"the block's own report agrees edge by edge and hands back the level it "
            f"solved: dtLevel {DT} s where its edge stepped, 0 where it did not")
    for a in ARMS:
        M_ = {i: value(E(a, N_STEPS, i), "M") for i in (FREE, PINNED)}
        require(all(close(M_[i], inertance(DUCT[n][3], AREA[n]), 1e-14)
                    and close(value(B(a, n), "inertance"), M_[i], 1e-14)
                    for i, n in enumerate(NAME)),
                f"{a}: each edge's row inertance is rho*L/A of the exact bore of "
                f"THAT edge's mesh ({M_[0]:.6g} for X1, {M_[1]:.6g} for X2), and the "
                "block measures the same number on its own grid")
        require({n: value(B(a, n), "inertanceProbe") for n in NAME}
                == {n: value(by_arm[a]["partner"][n], "probeInertance")
                    for n in NAME},
                f"{a}: and the slab each edge bills its acceleration over is the one "
                f"it announced at the start ({value(B(a, 'X1'), 'inertanceProbe'):.6g} "
                f"for X1, {value(B(a, 'X2'), 'inertanceProbe'):.6g} for X2)")

    # -------- 2. two attributions in one level: the halves are per edge, not per arm
    split = max(abs((value(E(a, n, i), "dpFriction") + value(E(a, n, i), "dpInertia"))
                    / value(E(a, n, i), "dpBlock") - 1.0) for a, n, i in ALL)
    require(split < 1e-12,
            f"dpFriction + dpInertia IS dpBlock in every state to {split:.1e} -- the "
            "partition is made inside each block, so a mixed level is two partitions, "
            "not one blended one")
    require(all(value(E(a, n, i), "dpInertia") == 0.0
                and value(E(a, n, i), "dpFriction") == value(E(a, n, i), "dpBlock")
                and value(E(a, n, i), "rowError") == 0.0 for a, n, i in QUASI),
            f"on every quasi-steady edge, in every arm including the mixed ones, the "
            f"inertial half is EXACTLY zero, the friction half is the whole drop, and "
            f"the row error is exactly zero ({len(QUASI)} states)")
    require(all(abs(value(E(a, n, i), "dpInertia")) > 0.0 for a, n, i in STEPPING),
            f"and on every stepping edge it is not: {len(STEPPING)} states all carry a "
            "non-zero inertial half, including the pinned ones whose rate is held")
    for a in MIXED:
        require(all((value(E(a, n, FREE), "dpInertia") == 0.0)
                    != (value(E(a, n, PINNED), "dpInertia") == 0.0)
                    for n in range(1, N_STEPS + 1)),
                f"{a}: at every one of the eight levels exactly one of the two edges "
                "has an inertial half -- the arm is mixed level by level, not "
                "time-averaged into one formulation")
    rowerr = max(abs(value(E(a, n, i), "rowDrop")
                     - (value(E(a, n, i), "dpNetwork")
                        + value(E(a, n, i), "M") * value(E(a, n, i), "dq") / DT))
                 for a, n, i in ALL)
    require(rowerr == 0.0,
            f"the node pair's claim is the friction it negotiated plus rho*L/A*dq/dt "
            f"in all {len(ALL)} states to {rowerr:.1e} Pa -- the row's inertance stays "
            "in on BOTH edges whether or not the block on that edge is paying its own")

    # ------- 3. which half the interface stops on is read block by block as well
    chord_bad = 0.0
    for a, n, i in ALL:
        l = E(a, n, i)
        chord = (value(l, "dpFriction") if EDGE_DDT[a][i] == "Euler"
                 else value(l, "dpBlock"))
        chord_bad = max(chord_bad, abs(value(l, "residual")
                                       - abs(value(l, "R") * value(l, "q") - chord)
                                       / max(abs(chord), 1e-300)))
    require(chord_bad < 1e-14,
            f"the stop is graded on the friction half of the stepping edge and on the "
            f"whole drop of the quasi-steady one IN THE SAME ARM: recomputing "
            f"|R*q - chord|/chord per edge reproduces the printed residual to "
            f"{chord_bad:.1e}")
    require(all(value(E(a, n, i), "residual") < TOL_INTERFACE for a, n, i in ALL),
            f"every state of every arm ends inside the {TOL_INTERFACE:.0e} tolerance, "
            "so a mixed level is a converged level")
    whole = {(a, i): abs(value(E(a, 1, i), "R") * value(E(a, 1, i), "q")
                         - value(E(a, 1, i), "dpBlock"))
             / value(E(a, 1, i), "dpBlock")
             for a in ARMS for i in (FREE, PINNED)}
    graded = {(a, i): value(E(a, 1, i), "residual")
              for a in ARMS for i in (FREE, PINNED)}
    stepped1 = [(a, i) for a in ARMS for i in (FREE, PINNED)
                if EDGE_DDT[a][i] == "Euler"]
    require(all(whole[(a, i)] > TOL_INTERFACE > graded[(a, i)] for a, i in stepped1),
            f"the choice is not cosmetic: at level 1 the {len(stepped1)} stepped edges "
            "are graded "
            + " / ".join(f"{graded[(a, i)] * 100:.4f}%" for a, i in stepped1)
            + " on their friction half, and the same states read "
            + " / ".join(f"{whole[(a, i)] * 100:.2f}%" for a, i in stepped1)
            + f" on their whole drop -- two to three orders above the "
            f"{TOL_INTERFACE:.0e} the run is asked to stop inside")
    quasi_same = max(abs(whole[(a, i)] - graded[(a, i)])
                     for a in ARMS for i in (FREE, PINNED)
                     if EDGE_DDT[a][i] == "steady")
    require(quasi_same < 1e-12,
            f"and on every quasi-steady edge the two readings are the SAME number "
            f"(worst {quasi_same:.1e}), because there the block's friction half IS its "
            "whole drop")

    # ---------------- 4. the edge you left alone is the control's computation
    # Each mixed arm has, on each of its two edges, a control that agrees with it on
    # THAT edge's scheme AND on that edge's shape -- lo_x1free against lo_both and
    # lo_steady, lo_x2pin against lo_steady and lo_both, parab_x1free against
    # parab_both and parab_steady. That is why the shaped family is in this case:
    # without parab_both the only thing to compare a shaped mixed edge against would
    # be a flat one, and a difference could not be attributed to the neighbour.
    fields = ("q", "dpBlock", "dpFriction", "dpInertia", "R", "residual", "sweeps",
              "rowError", "rowDrop", "dq", "betaOut", "pRisePa", "pRiseX",
              "fluxImposed", "unsteady")
    pin_rel = {}
    for a in MIXED:
        c0 = CONTROL[a][FREE]
        worst = {f: same(a, c0, FREE, FREE, f) for f in fields}
        require(all(v == 0.0 for v in worst.values()),
                f"{a} on X1 is not approximately {c0}, it is {c0}'s computation: "
                f"{len(fields)} fields x 8 levels all agree to exactly 0.0 -- rate, both "
                "halves of the drop, coefficient, residual, sweeps, row error, beta, "
                "recovery, flag")
        c1 = CONTROL[a][PINNED]
        dp = same(a, c1, PINNED, PINNED, "dpBlock")
        dr = same(a, c1, PINNED, PINNED, "R")
        rerr = same(a, c1, PINNED, PINNED, "rowError")
        dq = same(a, c1, PINNED, PINNED, "q")
        pin_dp = value(E(c1, N_STEPS, PINNED), "dpBlock")
        require(dq < 2e-20 and dp / pin_dp < 1e-3 and rerr < 1e-4,
                f"{a} on X2 reproduces {c1}'s computation on that edge to "
                f"{dp:.3e} Pa of {pin_dp:.1f} Pa ({dp / pin_dp:.1e} relative), the same "
                f"rate to {dq:.1e} and a row error of {rerr:.1e} -- hold an edge's own "
                "scheme and shape fixed and its neighbour's formulation does not reach it")
        rel_dp = dp / pin_dp
        rel_r = dr / value(E(c1, N_STEPS, PINNED), "R")
        pin_rel[a] = rel_dp
        require(rel_r < 1e-3 and dr > 1.0,
                f"but the coefficient the interface leaves on that untouched edge moves "
                f"by {dr:.1f} Pa.s/m3 out of {value(E(c1, N_STEPS, PINNED), 'R'):.4g} "
                f"({rel_r:.1e} relative) while its drop moves {rel_dp:.1e}: what "
                "composes is the edge's physics, not the price the interface ends on")
    wall = same("parab_x1free", "parab_steady", PINNED, PINNED, "sweeps")
    require(wall > 100.0,
            f"one thing does NOT compose, and it is the block budget on the shaped "
            f"pinned edge: the two arms differ by {wall:.0f} sweeps because the mixed "
            f"arm spends {sum(1 for k in range(1, N_STEPS + 1) if value(E('parab_x1free', k, PINNED), 'sweeps') == MAX_SWEEPS)} "
            f"of eight levels at the {MAX_SWEEPS:.0f}-sweep wall and its control only "
            f"{sum(1 for k in range(1, N_STEPS + 1) if value(E('parab_steady', k, PINNED), 'sweeps') == MAX_SWEEPS)}. "
            "That is 5.37's tolerance floor showing through: at a rate the network "
            "holds, a shaped inlet is solved to somewhere near the budget, so whether a "
            f"level reaches it is a coin the neighbour's trajectory can flip -- while "
            f"the drop it hands back still agrees to {pin_rel['parab_x1free']:.1e} "
            "relative")
    switched = {(i, tag): same(x, y, i, i, "dpBlock")
                for i, (tag, x, y) in
                ((FREE, ("flat", "lo_both", "lo_steady")),
                 (PINNED, ("flat", "lo_x1free", "lo_both")),
                 (FREE, ("shaped", "parab_both", "parab_steady")),
                 (PINNED, ("shaped", "parab_x1free", "parab_both")))}
    require(all(v > 600.0 for v in switched.values()),
            "for scale: the same edge with its OWN scheme flipped moves by "
            + "  ".join(f"{tag} {NAME[i]} {v:.1f} Pa"
                        for (i, tag), v in switched.items())
            + f" -- three to six orders more than what a mixed neighbour leaves on the "
            f"edge it does not name (worst {max(pin_rel.values()):.1e} relative of the "
            "same drop)")
    pin_both = same("lo_x2pin", "lo_both", PINNED, PINNED, "dpBlock")
    require(pin_both < 1e-3 and same("lo_x2pin", "lo_both", PINNED, PINNED, "sweeps")
            == 0.0,
            f"and on the mirror arm the stepped pinned edge is {pin_both:.1e} Pa from "
            "the all-Euler arm's stepped pinned edge with the same sweep count at every "
            "level, so the equality is with the edge's own scheme and not with a "
            "particular arm")

    # -------- 5. the double count is born on the edge that steps and stops there
    lvl1 = {(a, i): value(E(a, 1, i), "rowError") for a, i in
            ((a, i) for a in ARMS for i in (FREE, PINNED))
            if EDGE_DDT[a][i] == "Euler"}
    require(all(0.44 < abs(v) < 0.63 for v in lvl1.values()),
            "the first level is where the double count lives, and it is per edge: "
            + "  ".join(f"{a}/{NAME[i]} {v * 100:.2f}%" for (a, i), v in sorted(lvl1.items())))
    require(all(value(E(a, n, i), "rowError") == 0.0 for a, n, i in QUASI),
            "while the same eight levels of the same mixed arms report exactly 0.0 for "
            "the edge they left quasi-steady -- the double count does not leak across "
            "the node pair it shares with its neighbour")
    lvl8 = {(a, i): value(E(a, N_STEPS, i), "rowError") for a in ARMS
            for i in (FREE, PINNED) if EDGE_DDT[a][i] == "Euler"}
    require(all(1e-2 < abs(lvl8[(a, FREE)]) < 0.13 for a in ARMS
                if EDGE_DDT[a][FREE] == "Euler"),
            "the free edge still disagrees at the last level: "
            + "  ".join(f"{a}/X1 {lvl8[(a, FREE)] * 100:.3f}%" for a in ARMS
                        if EDGE_DDT[a][FREE] == "Euler"))
    require(all(abs(lvl8[(a, PINNED)]) < 1e-4 for a in ARMS
                if EDGE_DDT[a][PINNED] == "Euler"),
            "and the pinned one is down to "
            + "  ".join(f"{a}/X2 {abs(lvl8[(a, PINNED)]) * 100:.6f}%" for a in ARMS
                        if EDGE_DDT[a][PINNED] == "Euler")
            + " -- the term that doubled it is the one that stopped moving")
    pinned_states = [(a, n) for a in ARMS for n in range(2, N_STEPS + 1)]
    dqpin = max(abs(value(E(a, n, PINNED), "dq")) for a, n in pinned_states)
    dpi2 = {a: value(E(a, 2, PINNED), "dpInertia") for a in ARMS
            if EDGE_DDT[a][PINNED] == "Euler"}
    require(dqpin < 2e-20 and all(8.0 < v < 8.5 for v in dpi2.values()),
            f"yet the stepped pinned edge is a laboratory: its rate is 6.0e-5 at all "
            f"{len(pinned_states)} states from level 2 on (|dq| <= {dqpin:.2e}) and it "
            "still reports "
            + "  ".join(f"{a}/X2 dpInertia {v:.3f} Pa at level 2"
                        for a, v in dpi2.items())
            + " -- the jump out of rest at level 1 left a claim the block is still "
            "carrying, at a rate that has not moved since, and the shaped block pays a "
            "different one because its field is a different field")

    # ---- 6. the ramp belongs to the edge whose rate is free: split, tau, schedule
    m_by = inertance(BY[0], 0.25 * math.pi * BY[1] ** 2)
    inert_share = m_by / (inertance(DUCT["X1"][3], AREA["X1"]) + m_by)
    share = {a: value(E(a, 1, FREE), "q") / Q_IN for a in ARMS}
    require(share["lo_x2pin"] == share["lo_steady"]
            and share["lo_x1free"] == share["lo_both"]
            and share["parab_x1free"] == share["parab_both"],
            f"the first-level split of a mixed arm is EXACTLY the split of whichever "
            f"control has the same scheme AND the same shape on the FREE edge: "
            f"lo_steady and lo_x2pin both {share['lo_steady']:.6f}, lo_both and "
            f"lo_x1free both {share['lo_x1free']:.6f}, parab_both and parab_x1free both "
            f"{share['parab_x1free']:.6f} -- the edge the table does not name changes "
            "the split by 0.0")
    require(abs(share["lo_steady"] / inert_share - 1.0) < 0.03
            and abs(share["parab_steady"] / inert_share - 1.0) < 0.12
            and all(-0.195 < share[a] / inert_share - 1.0 < -0.13
                    for a in ("lo_both", "lo_x1free", "parab_both", "parab_x1free")),
            f"against the inertial split M_By/(M_X1+M_By) = {inert_share:.6f}: the "
            f"quasi-steady free edges sit {100 * (share['lo_steady'] / inert_share - 1):+.2f}% "
            f"(flat) and {100 * (share['parab_steady'] / inert_share - 1):+.2f}% (shaped) "
            "and the stepping ones below the flat control by "
            + " / ".join(f"{100 * (share[a] / inert_share - 1):+.2f}%"
                         for a in ("lo_both", "lo_x1free", "parab_both", "parab_x1free"))
            + " -- the whole-region acceleration the row does not credit is what "
            "short-changes the bypass at the first foot, and the shape moves the whole "
            "family up")
    dest = {a: value(E(a, N_STEPS, FREE), "q") / Q_IN for a in ARMS}
    require(dest["lo_x2pin"] == dest["lo_steady"]
            and dest["lo_x1free"] == dest["lo_both"]
            and dest["parab_x1free"] == dest["parab_both"],
            f"and so is where the ramp ends: {dest['lo_steady']:.6f} for the flat "
            f"quasi-steady pair, {dest['lo_both']:.6f} for the flat stepped pair, "
            f"{dest['parab_both']:.6f} for the shaped stepped pair, again to exactly 0.0")
    require(dest["parab_x1free"] > dest["lo_x1free"],
            f"within the mixed formulation the shape still buys the parallel branch a "
            f"bigger share of the destination ({dest['parab_x1free']:.6f} against "
            f"{dest['lo_x1free']:.6f})")

    def tau(a):
        d = [abs(value(E(a, n, FREE), "dq")) for n in range(1, N_STEPS + 1)]
        rat = [d[n] / d[n - 1] for n in range(3, N_STEPS)]
        g = sum(rat) / len(rat)
        return g, -DT / math.log(g)

    taus = {a: tau(a) for a in ARMS}
    require(taus["lo_x2pin"][0] == taus["lo_steady"][0]
            and taus["lo_x1free"][0] == taus["lo_both"][0]
            and taus["parab_x1free"][0] == taus["parab_both"][0],
            f"the tail decay ratio is the same number, not a similar one: "
            f"{taus['lo_steady'][0]:.5f} for lo_steady and lo_x2pin, "
            f"{taus['lo_both'][0]:.5f} for lo_both and lo_x1free, "
            f"{taus['parab_both'][0]:.5f} for parab_both and parab_x1free, so tau is "
            + " / ".join(f"{a} {taus[a][1]:.4f} s" for a in ARMS))
    require(1.20 < taus["lo_both"][1] / taus["lo_steady"][1] < 1.40
            and 1.15 < taus["parab_both"][1] / taus["parab_steady"][1] < 1.25,
            f"so the free edge's own du/dt lengthens the ramp by "
            f"+{100 * (taus['lo_both'][1] / taus['lo_steady'][1] - 1):.1f}% flat and "
            f"+{100 * (taus['parab_both'][1] / taus['parab_steady'][1] - 1):.1f}% shaped, "
            "while the pinned edge's adds nothing (the pairs are identical)")
    require(taus["parab_x1free"][1] / taus["lo_x1free"][1] > 1.8,
            f"and the shape still multiplies it by {taus['parab_x1free'][1] / taus['lo_x1free'][1]:.2f} "
            f"inside a mixed arm ({taus['parab_x1free'][1]:.4f} s against "
            f"{taus['lo_x1free'][1]:.4f} s) -- the same lever 5.37 and 5.38 measured, "
            "now on an arm whose second partner never left the chord")
    passes = {a: [int(val_) for val_ in
                  (value(by_arm[a]["step"][n], "passes") for n in range(1, N_STEPS + 1))]
              for a in ARMS}
    require(passes["lo_x2pin"] == passes["lo_steady"]
            and passes["lo_x1free"] == passes["lo_both"]
            and passes["parab_x1free"] == passes["parab_both"],
            "the exchange schedule too, level by level: "
            + f"{passes['lo_steady']} == {passes['lo_x2pin']}, "
            f"{passes['lo_both']} == {passes['lo_x1free']}, "
            f"{passes['parab_both']} == {passes['parab_x1free']}")
    tot = {a: sum(passes[a]) for a in ARMS}
    require(tot["lo_steady"] == 16 and tot["lo_x2pin"] == 16
            and tot["lo_both"] == 27 and tot["lo_x1free"] == 27
            and tot["parab_both"] == 37 and tot["parab_x1free"] == 37
            and tot["parab_steady"] == 18,
            f"totals " + "/".join(str(tot[a]) for a in ARMS) + " passes for the seven "
            f"arms: du/dt costs {tot['lo_both'] / tot['lo_steady']:.2f}x when it is put "
            f"on the edge that relaxes and 1.00x when it is put on the one that does "
            f"not, and the same holds shaped ({tot['parab_both'] / tot['parab_steady']:.2f}x)")
    require(all(by_arm[a]["solves"] == 2 * tot[a] for a in ARMS),
            f"and SIMPLE is still called exactly twice per exchange "
            f"({'/'.join(str(by_arm[a]['solves']) for a in ARMS)} solves), so the "
            "schedule is interface work, not block work")

    # ------------- 7. the acceleration the row is credited with is still shape-blind
    chi = {}
    for a in ARMS:
        if EDGE_DDT[a][FREE] != "Euler":
            continue
        l = E(a, 1, FREE)
        chi[a] = value(l, "dpInertia") / (PROBE_INERTANCE["X1"] * value(l, "dq") / DT)
    require(all(0.95 < v < 0.98 for v in chi.values())
            and max(chi.values()) - min(chi.values()) < 0.005,
            "on the mixed arm the level-1 acceleration of the shaped edge is "
            f"{chi['parab_x1free']:.6f} of rho*V_probe/A^2*dq/dt and of the flat one "
            f"{chi['lo_x1free']:.6f} -- {(max(chi.values()) - min(chi.values())) * 100:.2f} "
            "apart, the same blind 5.38 measured, and identical within each shape "
            "whether or not the other partner is in the time level -- taking the second "
            "partner out of the level did not teach the first one the shape")
    fric = {n: (value(E("parab_x1free", n, FREE), "dpFriction")
                / value(E("lo_x1free", n, FREE), "dpFriction") - 1.0) for n in (1, N_STEPS)}
    require(-0.18 < fric[1] < -0.12 and fric[N_STEPS] < -0.45,
            f"while the SAME switch on the SAME edge moves the friction half "
            f"{fric[1] * 100:+.1f}% at level 1 and {fric[N_STEPS] * 100:+.1f}% at the "
            "last: the shape is paid on the chord even in an arm whose other partner is "
            "quasi-steady")
    beta_err = max(abs(value(E(a, n, i), "betaImposed")
                       / betaDeclared(SHAPE[a][i], PANELS) - 1.0) for a, n, i in ALL)
    require(beta_err < 1e-12,
            f"the imposed factor is still the declared shape's discrete value on "
            f"{PANELS} panels to {beta_err:.1e}, on a stepping edge and a quasi-steady "
            "one alike")
    flux_err = max(abs(value(E(a, n, i), "fluxImposed") / value(E(a, n, i), "q") - 1.0)
                   for a, n, i in ALL)
    require(flux_err < 1e-12,
            f"and the section integral of the condition the inlet received is still "
            f"that level's rate to {flux_err:.1e} on both kinds of edge")

    # ---------------- 8. the plateau and the recovery are per-edge facts as well
    for a in ARMS:
        for i in (FREE, PINNED):
            v = [value(E(a, n, i), "betaOut") for n in range(1, N_STEPS + 1)]
            if EDGE_DDT[a][i] == "steady":
                require(max(v) - min(v) < 1.2e-4 and all(BETA_PLATEAU[0] < x
                                                         < BETA_PLATEAU[1] for x in v),
                        f"{a}/{NAME[i]}: a quasi-steady edge is AT its plateau from "
                        f"level 1 -- {v[0]:.8f} to {v[-1]:.8f}, spread "
                        f"{max(v) - min(v):.1e}, a number rather than a trajectory")
            else:
                require(v[0] < BETA_PLATEAU[0] and BETA_PLATEAU[0] < v[-1]
                        < BETA_PLATEAU[1],
                        f"{a}/{NAME[i]}: a stepping edge starts BELOW the plateau "
                        f"({v[0]:.8f}) and arrives at it ({v[-1]:.8f}) -- the region is "
                        "behind its own partner by an amount of level, not of shape")
    filled = {(a, n, i): value(E(a, n, i), "pRisePa") > 0.0 for a, n, i in ALL}
    above = {(a, n, i): value(E(a, n, i), "betaIn") > value(E(a, n, i), "betaOut")
             for a, n, i in ALL}
    require(sum(filled.values()) == len(SHAPED) * 2 * N_STEPS
            and filled == above,
            f"the sign law is per edge and per level too: exactly "
            f"{sum(filled.values())} of {len(ALL)} states report a pressure rise, they "
            f"are the {len(SHAPED)} shaped arms' two edges at all eight levels -- both "
            "formulations and both positions -- and they are the same states whose "
            "imposed beta is above what the region reads back")
    rise_par = {}
    for a in SHAPED:
        for i in (FREE, PINNED):
            v = [value(E(a, n, i), "pRisePa") for n in range(1, N_STEPS + 1)]
            sx = sorted({round(value(E(a, n, i), "pRiseX"), 4)
                         for n in range(1, N_STEPS + 1)})
            rise_par[(a, i)] = (v[0], v[-1], max(v) - min(v), sx)
    f_, p_ = rise_par[("parab_x1free", FREE)], rise_par[("parab_x1free", PINNED)]
    require(p_[2] < 2e-6 and len(p_[3]) == 1,
            f"on the quasi-steady pinned edge of the shaped mixed arm the recovery is "
            f"{p_[1]:.4f} Pa at one station {p_[3][0]} m at all eight levels, moving "
            f"{p_[2]:.1e} Pa between them -- a pinned rate through a fixed shape has no "
            "ramp to sign it")
    require(f_[2] > 100.0 and len(f_[3]) == 3,
            f"on the stepping free edge of the same arm it is a curve "
            f"{f_[0]:.4f} -> {f_[1]:.4f} Pa whose peak station walks "
            f"{' -> '.join(f'{x} m' for x in f_[3])} as the region develops: the two "
            "edges of one shaped arm differ only in whether their own du/dt is in the "
            "equation that produced their drop")
    rb = rise_par[("parab_both", PINNED)]
    require(len(rb[3]) == 2 and rb[2] > 100.0,
            f"and the mirror of that is the shaped pinned edge when it DOES step "
            f"(parab_both): its peak walks {rb[3]} and swings {rb[2]:.1f} Pa, so the "
            "station the recovery peaks on is that edge's own formulation talking, not "
            "the shape's and not the position's")
    rise_diff = same("parab_x1free", "parab_steady", PINNED, PINNED, "pRisePa")
    require(rise_diff < 1e-5,
            f"and the recovery of the quasi-steady pinned edge is separable too: "
            f"{rise_diff:.1e} Pa between the mixed arm and its shaped control, against "
            f"{f_[2]:.1f} Pa of curve on the edge that steps")

    # -------------------- 9. cost inside the block, and where the sweep cap bites
    sw = {(a, n): [int(value(E(a, k, NAME.index(n)), "sweeps"))
                   for k in range(1, N_STEPS + 1)] for a in ARMS for n in NAME}
    require(sw[("lo_x1free", "X2")] == sw[("lo_steady", "X2")]
            and sw[("lo_x2pin", "X2")] == sw[("lo_both", "X2")]
            and sw[("parab_x1free", "X1")] == sw[("parab_both", "X1")],
            f"sweep by sweep too: the quasi-steady pinned edge takes "
            f"{sw[('lo_steady', 'X2')][0]} at every level in both arms that leave it "
            f"steady, the stepping one {sw[('lo_both', 'X2')][0]} in both arms that step "
            f"it, and the shaped free edge {sw[('parab_both', 'X1')]} in both arms that "
            "step it -- the block budgets are per edge as well")
    capped = {(a, n): sum(1 for k in range(1, N_STEPS + 1)
                          if value(E(a, k, NAME.index(n)), "sweeps") == MAX_SWEEPS)
              for a in ARMS for n in NAME}
    unconv = {(a, n): sum(1 for k in range(1, N_STEPS + 1)
                          if value(E(a, k, NAME.index(n)), "blockConverged") == 0)
              for a in ARMS for n in NAME}
    require(capped[("parab_x1free", "X2")] == 6 and unconv[("parab_x1free", "X2")] == 5,
            f"the wall is still under the shaped pinned block: it exhausts all "
            f"{MAX_SWEEPS:.0f} sweeps on {capped[('parab_x1free', 'X2')]} of eight "
            f"levels and reports blockConverged 0 on {unconv[('parab_x1free', 'X2')]} "
            "of those six -- here the cap and the failure are NOT quite the same event, "
            "one level reaches the budget and still meets its tolerance")
    at_wall = {(a, n) for (a, n), v in capped.items() if v > 0}
    require(at_wall == {(a, "X2") for a in SHAPED}
            and capped[("parab_steady", "X2")] == 3
            and capped[("parab_both", "X2")] == 1,
            "and the wall is the shape's, in every formulation: "
            + "  ".join(f"{a}/X2 {capped[(a, 'X2')]} of eight levels at "
                        f"{MAX_SWEEPS:.0f} sweeps" for a in SHAPED)
            + f", while every flat arm and every X1 block never reaches it "
            f"(worst of the other eleven: {max(max(sw[(a, n)]) for (a, n) in sw if (a, n) not in at_wall)}"
            f" sweeps, and the mixed flat arms never exceed "
            f"{max(max(sw[(a, 'X1')]) for a in MIXED if a in MIXED and SHAPE[a] == ('flat', 'flat'))})")
    require(max(value(E(a, k, i), "residual") for a in ARMS
                for k in range(1, N_STEPS + 1) for i in (FREE, PINNED)) < TOL_INTERFACE,
            f"the worst residual over all {len(ALL)} states is "
            f"{max(value(E(a, k, i), 'residual') for a in ARMS for k in range(1, N_STEPS + 1) for i in (FREE, PINNED)):.3e}"
            ", so the capped levels are a block-budget fact and not a coupling failure")
    newt = max(value(by_arm[a]["step"][n], "newton") for a in ARMS
               for n in range(1, N_STEPS + 1))
    require(newt < 1e-11,
            f"the network's own linear algebra closes every level to {newt:.2e}; the "
            f"worst of each arm is "
            + " ".join(f"{a} {max(value(by_arm[a]['step'][k], 'newton') for k in range(1, N_STEPS + 1)):.1e}"
                       for a in ARMS)
            + " -- iteration is spent on the interface, never on the matrix")

    # ------------------ 10. the partner stays legal and the network still closes
    yp = {(a, n): value(B(a, n), "wallYPlusMax") for a in ARMS for n in NAME}
    require(all(20.0 < v < 35.0 for v in yp.values()),
            "the first cell centres stay in the log-law band on both formulations of "
            "both edges " + " ".join(f"{a}/{n} {v:.1f}" for (a, n), v in sorted(yp.items())))
    nr = {(a, n): value(B(a, n), "nutMax") / NU for a in ARMS for n in NAME}
    require(all(40.0 < v < 150.0 for v in nr.values()),
            "the eddy viscosity runs "
            + " ".join(f"{a}/{n} {v:.0f}x" for (a, n), v in sorted(nr.items()))
            + " the molecular one, so the wall functions are the closure being used on "
            "every edge of every arm")
    for a in ARMS:
        for n in NAME:
            gap_, span_, _, length = DUCT[n]
            dh = 2.0 * gap_ * span_ / (gap_ + span_)
            plug = value(B(a, n), "plugSpeed")
            require(RHO * plug * dh / MU > 4000.0,
                    f"{a} {n}: Re_Dh {RHO * plug * dh / MU:.0f} at the rate the ramp "
                    "ended on")
            require(abs(value(B(a, n), "fluxOut") / value(B(a, n), "fluxIn") - 1.0)
                    < 5e-11,
                    f"{a} {n}: what enters the region leaves it to "
                    f"{abs(value(B(a, n), 'fluxOut') / value(B(a, n), 'fluxIn') - 1):.1e}"
                    " -- neither a shape nor a time level creates mass")
            require(close(value(by_arm[a]["final"][n], "dp"),
                          value(B(a, n), "dpPa"), 1e-12),
                    f"{a} {n}: the drop the interface ended with is the drop the "
                    "block's own field measures")
            plain = abs(value(B(a, n), "probeSlope")) * gap_ / (2.0 * RHO
                                                               * value(B(a, n), "wallUTau") ** 2) - 1.0
            acc = RHO * value(E(a, N_STEPS, NAME.index(n)), "dq") / DT / gap_ / span_
            credited = (abs(value(B(a, n), "probeSlope")) * gap_
                        / (2.0 * (RHO * value(B(a, n), "wallUTau") ** 2
                                  + gap_ / 2.0 * acc)) - 1.0)
            require(-0.06 < plain < 0.19,
                    f"{a} {n}: the developed window's gradient balances the model's "
                    f"own wall shear to {plain * 100:+.2f}% without the level's "
                    "acceleration")
            i = NAME.index(n)
            if i == PINNED:
                require(abs(credited - plain) < 1e-10,
                        f"{a} {n}: on the pinned edge the credit is the SAME number "
                        f"({credited * 100:+.2f}% against {plain * 100:+.2f}%): its rate "
                        "has stopped moving, so there is nothing to add whether or not "
                        "the block carries du/dt")
            if i == FREE and EDGE_DDT[a][i] == "steady":
                require(abs(credited - plain) > 1e-3,
                        f"{a} {n}: and on the FREE edge of the same arm the credit moves "
                        f"even though the flag says the region paid no du/dt "
                        f"({plain * 100:+.2f}% -> {credited * 100:+.2f}%, "
                        f"{acc:+.1f} Pa/m) -- that term is the ROW's, from the network's "
                        "own rate change, and the block's formulation does not own it")
            if i == FREE and a in SHAPED and EDGE_DDT[a][i] == "Euler":
                require(abs(credited) < abs(plain),
                        f"{a} {n}: for the shaped stepping edge the credit runs toward "
                        f"zero, {plain * 100:+.2f}% -> {credited * 100:+.2f}%, from an "
                        f"acceleration of {acc:+.1f} Pa/m -- the term is real, it is a "
                        "few percent of the shear, and it lives only on the edge that "
                        "pays for it")
            if i == FREE and a not in SHAPED and EDGE_DDT[a][i] == "Euler":
                require(abs(credited) > abs(plain),
                        f"{a} {n}: on the flat stepping edge the same credit runs AWAY "
                        f"from zero, {plain * 100:+.2f}% -> {credited * 100:+.2f}%, from "
                        f"{acc:+.1f} Pa/m -- the sign of the correction is the region's "
                        "own history, not the formulation's")
    lead = pipeChord(*LEAD, Q_IN) * Q_IN
    for a in ARMS:
        require(close(nd(a, N_STEPS, "Src") - nd(a, N_STEPS, "J1"), lead, 1e-9)
                and close(nd(a, N_STEPS, "J3") - nd(a, N_STEPS, "Snk"), lead, 1e-9),
                f"{a}: the two lead pipes are the network's own Churchill drop, "
                f"{lead:.6f} Pa at the injected rate, and no shape or time level "
                "enters it")
        require(close(nd(a, N_STEPS, "Snk"), P_REF, 1e-12),
                f"{a}: Snk is held at the case's fixed pressure")
        require(all(close(nd(a, n, "J1") - nd(a, n, "J2"),
                          value(E(a, n, FREE), "rowDrop"), 1e-9)
                    and close(nd(a, n, "J2") - nd(a, n, "J3"),
                              value(E(a, n, PINNED), "rowDrop"), 1e-9)
                    for n in range(1, N_STEPS + 1)),
                f"{a}: every level's node pair is every edge's momentum row, inertial "
                "term included, on a mixed arm and on a uniform one alike")
        require(all(abs(pipeq(a, n, "A") - Q_IN) < 1e-18
                    and close(pipeq(a, n, "X1") + pipeq(a, n, "By"), Q_IN, 1e-11)
                    for n in range(1, N_STEPS + 1)),
                f"{a}: A takes the injected rate and the pair splits it at every level")

    print("note: per-edge scheme " + "  ".join(
        f"{a} {EDGE_DDT[a][0]}/{EDGE_DDT[a][1]}" for a in ARMS))
    print("note: split at level 1 / destination at level 8 -- " + "  ".join(
        f"{a} {share[a]:.6f}/{dest[a]:.6f}" for a in ARMS))
    print("note: tau from the free edge's tail -- " + "  ".join(
        f"{a} {taus[a][1]:.4f} s (ratio {taus[a][0]:.5f})" for a in ARMS))
    print("note: rowError at level 1, X1/X2 -- " + "  ".join(
        f"{a} {value(E(a, 1, FREE), 'rowError') * 100:+.2f}%/"
        f"{value(E(a, 1, PINNED), 'rowError') * 100:+.2f}%" for a in ARMS))
    print("note: passes per arm " + " ".join(f"{a} {tot[a]}" for a in ARMS))
    print("note: mixed arm vs its shape-matched control, worst dpBlock over 8 levels "
          "-- " + "  ".join(
              f"{a}/X1 {same(a, CONTROL[a][FREE], FREE, FREE, 'dpBlock'):.1e} "
              f"{a}/X2 {same(a, CONTROL[a][PINNED], PINNED, PINNED, 'dpBlock'):.3f}"
              for a in MIXED))
    print("note: shaped pinned block at the 300-sweep wall -- " + "  ".join(
        f"{a} {capped[(a, 'X2')]}/8 capped, {unconv[(a, 'X2')]} of them not converged"
        for a in SHAPED))

    if fails:
        print(f"\n{len(fails)} FAILURES:")
        for m in fails:
            print(" -", m)
        return 1
    print("\nPASS: a scheme declared per edge -- one level holds two formulations, the "
          "free edge of a mixed arm is its control's computation to exactly 0.0 on "
          "flat and shaped alike, its other edge agrees with its own control to under "
          "1e-3 relative, and what the mix changes is only which half the interface "
          "stops on, where the double count lives, and the price it ends on")
    return 0


if __name__ == "__main__":
    sys.exit(main())
