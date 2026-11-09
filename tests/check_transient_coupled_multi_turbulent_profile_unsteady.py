#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient-multi over two wall-function RANS ducts whose
inlet is handed a SHAPE and whose own fluid answers du/dt, and check what a shape
is worth once the partner it is handed to is not only nonlinear but unsteady.

The topology is the multi-edge one the profile chain has used since §5.34:
Src -> A -> { X1 (3D RANS duct) || By (1D bypass) } -> X2 (3D RANS duct) -> B -> Snk,
with a constant 6.0e-5 m3/s switched on at t = 0 so every rate starts at zero and
ramps over eight levels of dt = 0.2 s. Four arms differ only in the `profile`
table and `ddts` table: `lo_steady` flat with a quasi-steady partner (the control),
`parab_steady` parabolic with the same quasi-steady partner, then `lo_euler`,
`parab_euler` and `p17_euler` repeat all three shapes with `schemes{ ddt Euler; }`
honoured inside the block. So the case is a 2x2 with a third shape tacked on:
shape x (partner steps in time or does not).

§5.36 put a shape in front of a RANS partner at one operating point, §5.37 put a
shape inside a time level against a partner that carried du/dt, and §5.33 put a
du/dt-carrying partner in a time level without a shape. Only here are all three
present, and the two things that could only have been settled here are:

  1. the declaration, and which half of the interface it reaches - five arms, and
     each level reports both the shape word and the partner's ddt
  2. the trade still closes on the rate - the section integral of the condition
     the inlet received IS that level's rate, for the stepped partner too
  3. the block's two halves add up - dpFriction + dpInertia == dpBlock exactly,
     and a quasi-steady partner reports an acceleration of identically zero
  4. the row and the block BOTH bill the acceleration - the row's inertance is
     never removed when the partner starts answering du/dt, so the interface
     converges on the friction half while the inertial half disagrees by 44-63%
     at the first level and 1-12% at the last
  5. the acceleration is shape-blind - three shapes, one du/dt: the level's
     measured inertia agrees across them to 0.4 percentage points while the same
     three shapes move the friction half by more than 40%
  6. the pressure recovery is negotiated away, not carried by the row - the
     exact identity rowError*dpBlock = (R*q - dpFriction) + (M*dq/dt - dpInertia),
     with the recovery absent from every term
  7. the sign law and the plateau, now as functions of time - beta_out at the
     first level is 68-71% of the plateau gap when the partner steps and 100%
     when it does not, and the rise appears at a station that migrates
     downstream as the ramp runs
  8. the ramp: what the partner's own inertia costs - tau grows 19-24% when the
     block starts answering du/dt, the shape multiplies it by ~2 in both
     formulations, and the textbook secant formula becomes ACCURATE for the
     stepped partner (within 3%) where it over-counts the quasi-steady one by 23%
  9. the pinned edge is a clean laboratory - its rate is 6.0e-5 at all 40 states,
     so a region there accelerates (8 Pa of dpInertia) at zero rate change, and
     the formulation changes its drop by <= 0.17% while it changes the parallel
     edge's rate by 1.6%
 10. schedule, cost, and the walls - du/dt costs 1.7-2.1x more exchanges, the
     300-sweep cap and blockConverged 0 are the same event, and the dynamic
     inertance read off a pinned edge is not a number worth asserting
 11. the partner stays legal, and the network closes - log-law y+, the wall
     balance credited with the level's own acceleration, and the node pair still
     equal to the row

usage: check_transient_coupled_multi_turbulent_profile_unsteady.py <solver> <case>
"""
import math
import os
import re
import subprocess
import sys

# fluid, geometry and schedule of
# tests/cases/transientCoupledMultiTurbulentProfileUnsteady/case.dict
MU, RHO = 1.0e-3, 1000.0
NU = MU / RHO
Q_IN = 6.0e-5                  # m3/s injected at Src from t = 0
P_REF = 1.0e6                  # Pa held at Snk
DT, N_STEPS = 0.2, 8           # time { dt 0.2; endTime 1.6; }
MAX_PASSES = 12                # couple{ maxSubIterations }
TOL_INTERFACE = 1e-3           # couple{ tolerance }
MAX_SWEEPS = 300.0             # simple{ maxIterations }, each block's own budget

ARMS = ("lo_steady", "parab_steady", "lo_euler", "parab_euler", "p17_euler")
NAME = ("X1", "X2")
SHAPE = {"lo_steady": ("flat", "flat"),
         "parab_steady": ("parabolic", "parabolic"),
         "lo_euler": ("flat", "flat"),
         "parab_euler": ("parabolic", "parabolic"),
         "p17_euler": ("power1over7", "power1over7")}
DDT = {"lo_steady": "steady", "parab_steady": "steady",
       "lo_euler": "Euler", "parab_euler": "Euler", "p17_euler": "Euler"}
STEP = {"lo_euler", "parab_euler", "p17_euler"}          # the arms whose partner moves
SHAPED = {"parab_steady", "parab_euler"}                # the arms fed a parabola
SEED = {"X1": 1.4e7, "X2": 1.36e7}

# gap, span, streamwise cells, region length. The wall-bounded axis is y in both
# blocks and the z faces are symmetry planes, so a section is gap x span.
DUCT = {"X1": (0.005, 0.005, 20, 0.2), "X2": (0.004, 0.008, 24, 0.24)}
PANELS = 8                     # wall-normal panels the inlet patch is made of
BY = (0.3, 0.008)               # the 1D bypass, sharing X1's node pair
LEAD = (0.3, 0.02)              # A and B
AREA = {n: DUCT[n][0] * DUCT[n][1] for n in NAME}

CONTINUUM = {"flat": 1.0, "parabolic": 6.0 / 5.0,
             "power1over7": (7.0 / 9.0) / (7.0 / 8.0) ** 2}
# the plateau this wall-function channel develops to, measured over all 80 states
BETA_PLATEAU = (1.0200, 1.0211)
# rho*(beta_in - beta_out)*U^2 measured against the rise the region reports
RISE_ENTRY_QUASI = (0.970, 1.000)   # §5.36's band: a quasi-steady partner is there
RISE_ENTRY_STEPPED = (0.060, 0.960)  # and a stepped one is still building it
RISE_SERIES = (1.15, 1.20)


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


def stations(line):
    """The (x, y) pairs of a `transient betaProfile`/`pressureProfile` line."""
    tail = re.search(r": (.*)$", line).group(1).split()
    return [(float(tail[k]), float(tail[k + 1])) for k in range(0, len(tail), 2)]


def pass_nums(line, key):
    """The run of numbers stored under `key` in a pass line, up to the next word."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+((?:[-+0-9.eE]+\s*)+)", line)
    if m is None:
        raise SystemExit(f"pass line has no '{key}': {line}")
    return [float(v) for v in m.group(1).split()]


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
        raise SystemExit("transient shape x block du/dt solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    def check(cond, msg):
        if not cond:
            require(False, msg)

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
            cur = {"tag": ln, "arm": t[2].rstrip(":"), "shape": {}, "seed": {},
                   "partner": {}, "step": {}, "edge": {}, "pass": {}, "nodes": {},
                   "pipes": {}, "block": {}, "beta": {}, "press": {}, "final": {},
                   "summary": "", "solves": 0}
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
        elif ln.startswith("transient betaProfile "):
            cur["beta"][word(ln, "edge").rstrip(":")] = ln
        elif ln.startswith("transient pressureProfile "):
            cur["press"][word(ln, "edge").rstrip(":")] = ln
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

    # ---------------------------- 1. the declaration, and which half it reaches
    header = [ln for ln in lines if ln.startswith("transient network+blocks")]
    require(len(header) == 1 and "2 coupled edges, 5 starting sets" in header[0],
            "one network, two coupled edges and five starting sets announced once")
    require([c["arm"] for c in arms] == list(ARMS),
            f"the five arms ran, in the order the case declares them (got "
            f"{[c['arm'] for c in arms]})")
    require(all(c["seed"] == SEED for c in arms),
            f"every arm seeds both edges with the same pair ({SEED['X1']:.3e} then "
            f"{SEED['X2']:.3e}), so no arm starts at its own shaped resistance")
    require({c["arm"]: "accelerate yes" in c["tag"] for c in arms}
            == {a: False for a in ARMS},
            "no arm is accelerated: every difference below is the shape or the "
            "partner's own time level, not the update rule")
    for c in arms:
        require(c["shape"] == dict(zip(NAME, SHAPE[c["arm"]])),
                f"{c['arm']}: the case's `profile` reached its two edges as declared "
                f"({c['shape']['X1']} on both)")
        require({n: word(c["partner"][n], "ddt") for n in NAME}
                == {n: DDT[c["arm"]] for n in NAME},
                f"{c['arm']}: both partners were rebuilt under the table's word "
                f"`{DDT[c['arm']]}` -- the ddt the arm is given, not the case default")
    declared = {(a, n, i): word(E(a, n, i), "inletShape") for a, n, i in ALL}
    require(declared == {(a, n, i): SHAPE[a][i] for a, n, i in ALL},
            f"all {len(declared)} level x edge states report the shape that arm names:"
            " the interface re-hands the shape at every level, whether or not the"
            " region it lands in is stepping")
    flag = {(a, n, i): value(E(a, n, i), "unsteady") for a, n, i in ALL}
    require(flag == {(a, n, i): (1.0 if a in STEP else 0.0) for a, n, i in ALL},
            "and the same 80 states say whether the block carried du/dt at that "
            f"level: 48 with `unsteady 1`, {len(ALL) - 48} with `unsteady 0`, one "
            "flag per arm with no level left out")
    require(all(value(B(a, n), "unsteady") == (1.0 if a in STEP else 0.0)
                and value(B(a, n), "dtLevel") == (DT if a in STEP else 0.0)
                for a in ARMS for n in NAME),
            f"the block's own report agrees, and hands back the level it solved: "
            f"dtLevel {DT} s where it stepped and 0 where it did not")
    for a in ARMS:
        M_ = {i: value(E(a, N_STEPS, i), "M") for i in (0, 1)}
        require(all(close(M_[i], inertance(DUCT[n][3], AREA[n]), 1e-14)
                    and close(value(B(a, n), "inertance"), M_[i], 1e-14)
                    for i, n in enumerate(NAME)),
                f"{a}: the row's inertance is rho*L/A of the exact bores the blocks "
                f"mesh ({M_[0]:.6g} for X1, {M_[1]:.6g} for X2) and the block's own "
                "mesh measures the same numbers")
        require({n: value(B(a, n), "inertanceProbe") for n in NAME}
                == {n: value(by_arm[a]["partner"][n], "probeInertance") for n in NAME},
                f"{a}: and the probe slab the block bills its acceleration over is "
                f"the one it announced at the start "
                f"({value(B(a, 'X1'), 'inertanceProbe'):.6g} for X1, "
                f"{value(B(a, 'X2'), 'inertanceProbe'):.6g} for X2), "
                f"= {(value(B(a, 'X1'), 'inertanceProbe')
                     / value(E(a, N_STEPS, 0), 'M')) * 100:.1f}% of X1's region")

    # ------------------------------- 2. what the trade preserves: the rate, and the factor
    flux_err = max(abs(value(E(a, n, i), "fluxImposed")
                       / value(E(a, n, i), "q") - 1.0) for a, n, i in ALL)
    require(flux_err < 1e-12,
            f"the section integral of the condition the inlet received IS that "
            f"level's rate, for a stepping region too: worst {flux_err:.2e} over 80 "
            "states -- a shape is bought with a shape, not with a different amount")
    beta_err = max(abs(value(E(a, n, i), "betaImposed")
                       / betaDeclared(SHAPE[a][i], PANELS) - 1.0) for a, n, i in ALL)
    require(beta_err < 1e-12,
            f"and the factor read back off the boundary faces is the declared "
            f"profile's discrete value on {PANELS} panels to {beta_err:.1e} -- "
            "re-normalised per level, and identical under either partner")
    for s in ("flat", "parabolic", "power1over7"):
        disc = betaDeclared(s, PANELS)
        require(abs(disc / CONTINUUM[s] - 1.0) < 0.02,
                f"the {PANELS}-panel value of `{s}` is {disc:.6f} against the "
                f"continuum {CONTINUUM[s]:.6f}")
    mfp_err = max(abs(value(B(a, n), "momentumFluxPa")
                      / (RHO * value(B(a, n), "betaImposed")
                         * value(B(a, n), "fluxImposed") ** 2
                         / value(B(a, n), "inletArea")) - 1.0)
                  for a in ARMS for n in NAME)
    require(mfp_err < 1e-12,
            f"the ledger closes on the block's own faces: momentumFluxPa is "
            f"rho*beta*q^2/A to {mfp_err:.1e}")

    # --------------------------- 3. the block's two halves, and what a steady one says
    split = max(abs((value(E(a, n, i), "dpFriction") + value(E(a, n, i), "dpInertia"))
                    / value(E(a, n, i), "dpBlock") - 1.0) for a, n, i in ALL)
    require(split < 1e-12,
            f"dpFriction + dpInertia IS dpBlock, in every state, to {split:.1e} -- "
            "the split the interface reports is a partition of the drop the network "
            "was handed, not a third quantity")
    for a in ("lo_steady", "parab_steady"):
        require(all(value(E(a, n, i), "dpInertia") == 0.0
                    and value(E(a, n, i), "rowError") == 0.0
                    and value(E(a, n, i), "dpFriction") == value(E(a, n, i), "dpBlock")
                    for n in range(1, N_STEPS + 1) for i in (0, 1)),
                f"{a}: a quasi-steady partner reports exactly zero acceleration, "
                "exactly zero row error, and its whole drop in the friction half -- "
                "the control is not approximately the other case, it is the other "
                "case with the term removed")
    for a in STEP:
        require(all(abs(value(E(a, n, i), "dpInertia")) > 0.0
                    for n in range(1, N_STEPS + 1) for i in (0, 1)),
                f"{a}: a stepping partner has a non-zero inertial half at every "
                "level of both edges, including the pinned one")

    # --------------------- 4. the row and the block both bill the acceleration
    for a in ARMS:
        for n in range(1, N_STEPS + 1):
            for i in (0, 1):
                l = E(a, n, i)
                require(close(value(l, "rowDrop"),
                              value(l, "dpNetwork")
                              + value(l, "M") * value(l, "dq") / DT, 0.0),
                          f"{a} level {n} edge {i}: the node pair's claim is the "
                          "friction it negotiated plus rho*L/A*dq/dt -- the row's "
                          "inertance stayed in even after the block started paying")
    gap = {(a, n, i): value(E(a, n, i), "rowDrop") - value(E(a, n, i), "dpBlock")
           for a, n, i in ALL}
    inert_claim = {(a, n, i): value(E(a, n, i), "M") * value(E(a, n, i), "dq") / DT
                   for a, n, i in ALL}
    for a in ("lo_steady", "parab_steady"):
        worst = max(abs(gap[(a, n, i)] - (inert_claim[(a, n, i)]
                                          - (value(E(a, n, i), "dpBlock")
                                             - value(E(a, n, i), "dpNetwork"))))
                    for n in range(1, N_STEPS + 1) for i in (0, 1))
        require(worst < 1e-9,
                f"{a}: with a quasi-steady partner the row's whole inertial claim is "
                f"unbilled -- the gap between its drop and the block's is "
                f"M*dq/dt to {worst:.1e} Pa, and rowError is printed as 0 for it")
    lvl1 = {(a, i): value(E(a, 1, i), "rowError") for a in STEP for i in (0, 1)}
    lvl8 = {(a, i): value(E(a, N_STEPS, i), "rowError") for a in STEP for i in (0, 1)}
    require(all(0.44 < v < 0.63 for v in lvl1.values()),
            "the first level is where the double count lives: the row's claim and "
            "the block's disagree by "
            + " ".join(f"{a}/{NAME[i]} {v * 100:.1f}%" for (a, i), v in sorted(lvl1.items())))
    require(all(v < 0.125 for v in lvl8.values()),
            "and it drains along the ramp: "
            + " ".join(f"{a}/{NAME[i]} {value(E(a, N_STEPS, i), 'rowError') * 100:.3f}%"
                       for a in STEP for i in (0, 1)))
    require(all(1e-2 < lvl8[(a, 0)] < 0.125 for a in STEP)
            and all(lvl8[(a, 1)] < 1e-4 for a in STEP),
            f"the two halves drain at different speeds: X1 still disagrees by "
            f"{lvl8[('parab_euler', 0)] * 100:.1f}% / {lvl8[('p17_euler', 0)] * 100:.1f}%"
            f" / {lvl8[('lo_euler', 0)] * 100:.1f}% where the rate is its own, while the "
            f"pinned edge is down to {max(lvl8[(a, 1)] for a in STEP) * 100:.5f}% -- the "
            "term that doubled it is the one that stopped moving")
    require(abs(lvl1[("parab_euler", 0)] - lvl1[("lo_euler", 0)]) > 0.08
            and abs(lvl1[("p17_euler", 0)] - lvl1[("lo_euler", 0)]) < 0.01,
            f"the shape sets how much is double-billed only when it is a big shape: "
            f"parab {lvl1[('parab_euler', 0)] * 100:.1f}% against lo "
            f"{lvl1[('lo_euler', 0)] * 100:.1f}% and p17 {lvl1[('p17_euler', 0)] * 100:.1f}%"
            " at level 1 -- a 1/7 inlet is within a point of a plug here, because the "
            "double count is the row's inertance term, and the shape only moves the "
            "rate the level puts through it")
    require(all(abs(value(E(a, n, i), "residual")) < TOL_INTERFACE for a, n, i in ALL),
            f"yet the interface still stopped: every one of the 80 states is inside "
            f"{TOL_INTERFACE}")
    # coupledNetwork.cpp:1182 -- the stop is graded on the friction half alone, and
    # for a quasi-steady partner that half is the whole drop
    chord_of = lambda l: (value(l, "dpFriction") if value(l, "unsteady")
                          else value(l, "dpBlock"))
    res_err = max(abs(value(E(a, n, i), "residual")
                      - abs(value(E(a, n, i), "R") * value(E(a, n, i), "q")
                            - chord_of(E(a, n, i)))
                      / chord_of(E(a, n, i)))
                  for a, n, i in ALL)
    require(res_err < 1e-15,
            f"and what it stops on is |R*q - chord|/chord with the chord the "
            f"friction half: the identity holds to {res_err:.1e} over all 80 states, "
            f"while the same states' inertial halves disagree by up to "
            f"{max(lvl1.values()) * 100:.1f}% -- the loop converges on one half of the "
            "momentum and never looks at the other")

    # --------------------- 5. the acceleration is blind to the shape; friction is not
    probe = {n: value(B(a, n), "inertanceProbe") for a in ARMS for n in NAME}
    chi = {(a, n, i): (value(E(a, n, i), "dpInertia")
                       / (probe[NAME[i]] * value(E(a, n, i), "dq") / DT)
                       if abs(value(E(a, n, i), "dq")) > 1e-15 else None)
           for a in STEP for n in range(1, N_STEPS + 1) for i in (0, 1)}
    first = {i: [chi[(a, 1, i)] for a in STEP] for i in (0, 1)}
    spread = {i: max(first[i]) - min(first[i]) for i in (0, 1)}
    require(spread[0] < 0.005 and spread[1] < 0.002,
            f"at the level the ramp actually moves (level 1, dq = 6e-5 on X2) the "
            f"measured inertia is {first[0][0]:.4f} .. {first[0][2]:.4f} of "
            f"rho*V_probe/A^2*dq/dt on X1 for the three shapes -- a spread of "
            f"{spread[0] * 100:.2f} percentage points -- and "
            f"{first[1][0]:.4f} .. {first[1][2]:.4f} on X2, {spread[1] * 100:.2f} pp")
    fric1 = {i: (value(E("parab_euler", 1, i), "dpFriction")
                 / value(E("lo_euler", 1, i), "dpFriction") - 1.0) for i in (0, 1)}
    fric8 = {i: (value(E("parab_euler", N_STEPS, i), "dpFriction")
                 / value(E("lo_euler", N_STEPS, i), "dpFriction") - 1.0) for i in (0, 1)}
    require(fric1[0] < -0.14 and fric8[0] < -0.50 and fric8[1] < -0.52,
            f"on the same states the friction half moves "
            f"{fric1[0] * 100:+.1f}% at level 1 and {fric8[0] * 100:+.1f}% at level "
            f"{N_STEPS} on X1 ({fric8[1] * 100:+.1f}% on X2) -- the shape is worth "
            "more than half the drop and none of the inertia")
    late = {i: max(abs(chi[(a, n, i)] - 1.0) for a in STEP
                   for n in range(3, N_STEPS + 1) if chi[(a, n, i)] is not None)
            for i in (0,)}
    require(all(v < 0.025 for v in late.values()),
              f"once the ramp is on its own tail (levels 3-{N_STEPS}) X1's measured "
              f"inertia is within {max(late.values()) * 100:.2f}% of the probe slab's "
              "own rho*V/A^2*dq/dt in all three arms -- du/dt is a linear term on a "
              "section that barely changes shape between levels, and the pinned edge "
              "is left out here because its dq is zero and the ratio undefined")
    par = [abs(chi[("parab_euler", n, 0)] - 1.0) for n in range(3, N_STEPS + 1)]
    lo_ = [abs(chi[("lo_euler", n, 0)] - 1.0) for n in range(3, N_STEPS + 1)]
    require(all(p < l for p, l in zip(par, lo_)),
            "and the parabola is the closest of the three to that prediction at "
            "every one of those levels: "
            + " ".join(f"L{n} {chi[('parab_euler', n, 0)]:.4f} vs "
                       f"{chi[('lo_euler', n, 0)]:.4f}" for n in range(3, N_STEPS + 1)))
    dyn = {a: value(B(a, "X1"), "inertanceDynamic")
           / value(B(a, "X1"), "inertanceProbe") for a in STEP}
    require(all(0.999 < v < 1.006 for v in dyn.values()),
            "the block's own dynamic-over-probe inertance at the last solve is "
            + " ".join(f"{a} {v:.5f}" for a, v in dyn.items())
            + " -- inside 0.6% of the constant it announced, where §5.33's "
              "quasi-steady reading was allowed 2.8%")

    # --------------------- 6. the recovery is inside the friction half, not the row
    ident = max(abs(value(E(a, n, i), "rowError") * value(E(a, n, i), "dpBlock")
                    - ((value(E(a, n, i), "R") * value(E(a, n, i), "q")
                        - value(E(a, n, i), "dpFriction"))
                       + (value(E(a, n, i), "M") * value(E(a, n, i), "dq") / DT
                          - value(E(a, n, i), "dpInertia"))))
                for a in STEP for n in range(1, N_STEPS + 1) for i in (0, 1))
    require(ident < 1e-11,
            f"the row error is exactly the attribution arithmetic, and nothing else: "
            f"rowError*dpBlock == (R*q - dpFriction) + (M*dq/dt - dpInertia) to "
            f"{ident:.1e} Pa over 48 states")
    ins = {(a, n, i): value(E(a, n, i), "betaIn") for a, n, i in ALL}
    outs = {(a, n, i): value(E(a, n, i), "betaOut") for a, n, i in ALL}
    full = {(a, n, i): ins[(a, n, i)] - outs[(a, n, i)] for a, n, i in ALL}
    rise = {(a, n, i): value(E(a, n, i), "pRisePa") for a, n, i in ALL}
    k = ("parab_euler", N_STEPS, 0)
    terms = ((value(E(*k), "R") * value(E(*k), "q") - value(E(*k), "dpFriction")),
             (value(E(*k), "M") * value(E(*k), "dq") / DT - value(E(*k), "dpInertia")))
    require(abs(terms[0]) + abs(terms[1]) < 0.16 * rise[k],
            f"and the shape's pressure recovery is nowhere in it: {k[0]}'s parallel "
            f"edge ends with {rise[k]:.2f} Pa of recovery on a "
            f"{value(E(*k), 'dpBlock'):.2f} Pa drop "
            f"({rise[k] / value(E(*k), 'dpBlock') * 100:.1f}% of it), yet the row and "
            f"the block differ by only {terms[0]:.2f} + {terms[1]:.2f} = "
            f"{sum(terms):.2f} Pa -- the chord absorbs the recovery into the "
            "friction half it negotiates, so the row never sees the shape at all")

    # --------------------- 7. the sign law and the plateau, now in time
    over = {kk for kk, v in full.items() if v > 0}
    under = {kk for kk, v in full.items() if v <= 0}
    require(len(over) == 32 and len(under) == 48
            and over == {(a, n, i) for a in SHAPED
                         for n in range(1, N_STEPS + 1) for i in (0, 1)},
            f"the inlet is fuller than the plateau in exactly the {len(over)} "
            f"parabolic-fed states and poorer in {len(under)}: the plug and 1/7 "
            "inlets never hand over more than the region develops to, at any level, "
            "under either partner")
    require(all((rise[kk] > 0.0) == (kk in over) for kk in ALL),
            f"and pRisePa is strictly positive in exactly those {len(over)} states "
            f"and exactly 0.0 in the other {len(under)} -- the sign law survived "
            "being asked 80 times, twice per level and once per formulation")
    QUASI = ("lo_steady", "parab_steady")
    qvals = {kk: v for kk, v in outs.items() if kk[0] in QUASI}
    require(all(BETA_PLATEAU[0] <= v <= BETA_PLATEAU[1] for v in qvals.values()),
            f"with a quasi-steady partner every level is already on the plateau: all "
            f"{len(qvals)} of {QUASI}'s beta_out read inside {BETA_PLATEAU} "
            f"({min(qvals.values()):.8f} .. {max(qvals.values()):.8f}) -- the region's "
            "own end state, reached at the first level and held there")
    for a in QUASI:
        for i, n in enumerate(NAME):
            v = [outs[(a, m, i)] for m in range(1, N_STEPS + 1)]
            require(max(v) - min(v) < 1e-4,
                    f"{a} {n}: and it is one number across the ramp, not a trajectory "
                    f"({max(v):.8f} .. {min(v):.8f}, spread {max(v) - min(v):.1e}) -- "
                    "a partner with no time level of its own is never behind")
    svals = {kk: v for kk, v in outs.items() if kk[0] in STEP}
    first = [v for kk, v in svals.items() if kk[1] == 1]
    require(all(v < BETA_PLATEAU[0] for v in first),
            f"the stepping arms' first level is below that band by design "
            f"({min(first):.8f} .. {max(first):.8f}): the region was handed the inlet "
            "and has not yet relaxed it, so what level 1 reports is development, not "
            "the plateau")
    require(max(svals.values()) <= max(qvals.values()) + 1e-5,
            f"but by the last levels a stepping region is on the same plateau as its "
            f"quasi-steady twin, not above it ({max(svals.values()):.8f} against "
            f"{max(qvals.values()):.8f}) -- the partner's formulation decides WHEN the "
            "end state is reached, never WHICH end state")
    dev = {(a, i): (outs[(a, 1, i)] - 1.0) / (outs[(a, N_STEPS, i)] - 1.0)
           for a in STEP for i in (0, 1)}
    require(all(0.66 < dev[(a, 0)] < 0.72 for a in STEP)
            and all(0.82 < dev[(a, 1)] < 0.85 for a in STEP),
            "a stepping region has not finished developing at the first level: "
            f"beta_out is {dev[('lo_euler', 0)] * 100:.1f}% "
            f"{dev[('parab_euler', 0)] * 100:.1f}% "
            f"{dev[('p17_euler', 0)] * 100:.1f}% of the plateau gap on X1 and "
            f"{dev[('lo_euler', 1)] * 100:.1f}% on X2")
    require(all(abs(outs[(a, 1, 0)] - outs[(b, 1, 0)]) < 1e-4
                for a in STEP for b in STEP)
            and all(abs(outs[(a, 1, 1)] - outs[(b, 1, 1)]) < 3e-4 for a in STEP for b in STEP),
            f"and that unfinished state is the same whatever the inlet was handed: "
            f"X1's beta_out at level 1 is {outs[('lo_euler', 1, 0)]:.6f} / "
            f"{outs[('parab_euler', 1, 0)]:.6f} / {outs[('p17_euler', 1, 0)]:.6f} "
            "across a 0.18 spread of imposed beta -- the level, not the shape, owns "
            "the development when the region answers du/dt")
    require(all(outs[(a, 1, i)] > outs[(b, 1, i)] for a in ("lo_steady", "parab_steady")
                for b in STEP for i in (0, 1))
            and all(dev2 > 0.999 for dev2 in
                    [(outs[(a, 1, i)] - 1.0) / (outs[(a, N_STEPS, i)] - 1.0)
                     for a in ("lo_steady", "parab_steady") for i in (0, 1)]),
            "the quasi-steady control finishes that development inside its first "
            "level (>99.9% of the gap): the same dt, the same rate, the same inlet, "
            "and the difference is one term in the block")
    for a in STEP:
        for i, n in enumerate(NAME):
            if rise[(a, N_STEPS, i)] > 0.0:
                rat = [rise[(a, m, i)] / (RHO * full[(a, m, i)]
                                          * (value(E(a, m, i), "q")
                                             / value(E(a, m, i), "inletArea")) ** 2)
                       for m in range(1, N_STEPS + 1) if full[(a, m, i)] > 0]
                band = RISE_SERIES if i else RISE_ENTRY_STEPPED
                if i:
                    require(rat[0] < band[0],
                            f"{a} {n}: and its very first level is the one exception "
                            f"-- {rat[0]:.4f}, short of {band[0]}, because a region "
                            "that only started answering du/dt one level ago has not "
                            "yet developed the plateau the shape is measured against")
                    rat = rat[1:]
                require(all(band[0] < x < band[1] for x in rat),
                        f"{a} {n}: the rise is {min(rat):.4f} .. {max(rat):.4f} of "
                        f"rho*(beta_in-beta_out)*U^2 over the {len(rat)} levels it is "
                        "checked on -- the entry region relaxing the inlet toward its "
                        "own shape, one level at a time")
                if i == 0:
                    require(rat == sorted(rat) and rat[0] < 0.10,
                            f"{a} {n}: and it is a curve, not a constant -- the "
                            f"fraction climbs from {rat[0]:.4f} at the first level "
                            f"to {rat[-1]:.4f} at the last, because at 0.2 s the "
                            "region has only begun to relax what it was handed")
                else:
                    require(max(rat) - min(rat) < 0.01,
                            f"{a} {n}: the pinned edge sits on one ratio from level "
                            f"2 on ({min(rat):.4f} .. {max(rat):.4f}) -- its "
                            "rate never moved, so nothing about it was allowed to")
    rat_q = [rise[("parab_steady", m, 0)] / (RHO * full[("parab_steady", m, 0)]
                                             * (value(E("parab_steady", m, 0), "q")
                                                / value(E("parab_steady", m, 0),
                                                        "inletArea")) ** 2)
             for m in range(1, N_STEPS + 1)]
    require(all(RISE_ENTRY_QUASI[0] < x < RISE_ENTRY_QUASI[1] for x in rat_q),
            f"parab_steady X1: the same instrument on the quasi-steady partner reads "
            f"{min(rat_q):.4f} .. {max(rat_q):.4f} -- §5.36's band, so the reason the "
            "stepped arm falls short of it is the ramp, not the law")
    px = {"parab_euler": (0.015, 0.035, 0.045), "parab_steady": (0.045, 0.045, 0.045)}
    same_x = lambda x_, y_: abs(x_ - y_) < 1e-9
    for a, seq in px.items():
        seen = [value(E(a, m, 0), "pRiseX") for m in range(1, N_STEPS + 1)]
        require(same_x(seen[0], seq[0]) and same_x(seen[1], seq[1])
                and same_x(seen[2], seq[1]) and all(same_x(x, seq[2]) for x in seen[3:]),
                f"{a} X1: the station the rise peaks at walks downstream as the "
                f"region develops -- {seq[0]} m at level 1, {seq[1]} m at levels "
                f"2-3, {seq[2]} m from level 4 on; the quasi-steady arm sits at "
                f"{px['parab_steady'][0]} m for all eight")
    for a in ARMS:
        for i, n in enumerate(NAME):
            if rise[(a, N_STEPS, i)] == 0.0:
                require(value(E(a, N_STEPS, i), "pRiseX")
                        == stations(by_arm[a]["press"][n])[0][0],
                        f"{a} {n}: where nothing is handed back the reported station "
                        "is the region's own first one -- the instrument cannot "
                        "invent a rise")
                continue
            ps = stations(by_arm[a]["press"][n])
            k8 = max(range(len(ps)), key=lambda j: ps[j][1])
            require(close(rise[(a, N_STEPS, i)], ps[k8][1] - ps[0][1], 0.0)
                    and ps[k8][0] == value(E(a, N_STEPS, i), "pRiseX"),
                    f"{a} {n}: level {N_STEPS}'s rise and peak station are the "
                    "last solve's own pressure table, bit for bit")
        for n in NAME:
            st = stations(by_arm[a]["beta"][n])
            bs = [b for _, b in st]
            gap_, span_, ncell, length = DUCT[n]
            require(len(st) == ncell
                    and all(abs(xs - (j + 0.5) * length / ncell) < 1e-9
                            for j, (xs, _) in enumerate(st)),
                    f"{a} {n}: the region reports one shape station per cross-section "
                    f"({len(st)} over {length} m), each at its own cell centre")
            require(bs[-1] == outs[(a, N_STEPS, NAME.index(n))],
                    f"{a} {n}: the last station IS the betaOut the interface logged "
                    "at the level that ended the arm, bit for bit")
            if SHAPE[a][0] == "parabolic":
                require(bs[0] > bs[-1]
                        and max(bs[j] - bs[j + 1] for j in range(len(bs) - 1)) > 0.005,
                        f"{a} {n}: fed a parabola, the region relaxes it -- "
                        f"{bs[0]:.6f} -> {bs[-1]:.6f} along x")
            else:
                require(min(bs[j + 1] - bs[j] for j in range(len(bs) - 1)) > -1e-4
                        and bs[-1] - bs[0] > 0.005,
                        f"{a} {n}: a {SHAPE[a][0]} inlet is built by the region -- "
                        f"{bs[0]:.6f} -> {bs[-1]:.6f} with no fall bigger than 1e-4")

    # --------------------- 8. the ramp: what the partner's inertia costs
    def tau_fit(a, i=0):
        dq = [abs(value(E(a, n, i), "dq")) for n in range(1, N_STEPS + 1)]
        rt = [dq[j + 1] / dq[j] for j in range(2, 6)]   # L4/L3 .. L7/L6
        mean = sum(rt) / len(rt)
        return mean, -DT / math.log(mean), max(abs(x / mean - 1.0) for x in rt)
    taum = {}
    for a in ARMS:
        r, t, dev_ = tau_fit(a)
        taum[a] = t
        require(dev_ < 0.08,
                f"{a}: the tail of the ramp is geometric -- |dq| decays by "
                f"{r:.5f} +/- {dev_ * 100:.1f}% over levels 3-{N_STEPS}, so one time "
                f"constant {t:.4f} s describes it")
    require(1.90 < taum["parab_euler"] / taum["lo_euler"] < 2.04
            and 1.97 < taum["parab_steady"] / taum["lo_steady"] < 2.10,
            f"the shape multiplies the time constant by about the same factor under "
            f"either partner: {taum['parab_euler'] / taum['lo_euler']:.4f} when the "
            f"block steps ({taum['lo_euler']:.4f} -> {taum['parab_euler']:.4f} s) "
            f"and {taum['parab_steady'] / taum['lo_steady']:.4f} when it does not "
            f"({taum['lo_steady']:.4f} -> {taum['parab_steady']:.4f} s)")
    require(1.03 < taum["p17_euler"] / taum["lo_euler"] < 1.07,
            f"a 1/7 inlet moves it {taum['p17_euler'] / taum['lo_euler'] - 1:+.2%} "
            f"({taum['p17_euler']:.4f} s), "
            f"{(taum['p17_euler'] / taum['lo_euler'] - 1) / (taum['parab_euler'] / taum['lo_euler'] - 1) * 100:.0f}%"
            " of the parabola's move, which is where a shape only "
            f"{abs(full[('p17_euler', N_STEPS, 0)]):.4f} of beta below the plateau belongs")
    up = {s: taum[s + "_euler"] / taum[s + "_steady"] for s in ("lo", "parab")}
    require(all(1.17 < v < 1.26 for v in up.values()),
            f"and the block's own du/dt lengthens the ramp "
            f"{(up['lo'] - 1) * 100:+.1f}% for a plug and {(up['parab'] - 1) * 100:+.1f}% "
            f"for a parabola -- an almost formulation-blind surcharge, where the "
            "shape is a factor")
    m_x1 = inertance(DUCT["X1"][3], AREA["X1"])
    m_by = inertance(BY[0], 0.25 * math.pi * BY[1] ** 2)
    form = {}
    for a in ARMS:
        r1 = value(E(a, N_STEPS, 0), "R")
        r_by = (nd(a, N_STEPS, "J1") - nd(a, N_STEPS, "J2")) / pipeq(a, N_STEPS, "By")
        form[a] = (m_x1 + m_by) / (r1 + r_by)
    for a in ARMS:
        require(0.97 < form[a] / taum[a] < 1.30,
                f"{a}: the textbook (M_X1+M_By)/(R_X1+R_By) = {form[a]:.4f} s is "
                f"{form[a] / taum[a]:.4f}x the measured {taum[a]:.4f} s")
    require(all(abs(form[a] / taum[a] - 1.0) < 0.03 for a in STEP)
            and all(form[a] / taum[a] > 1.22 for a in ("lo_steady", "parab_steady")),
            "which is the finding: with the partner quasi-steady the secant formula "
            f"over-counts by {form['lo_steady'] / taum['lo_steady'] - 1:+.1%} and "
            f"{form['parab_steady'] / taum['parab_steady'] - 1:+.1%}, but once the "
            "block carries the acceleration it is right to "
            + " ".join(f"{form[a] / taum[a] - 1:+.2%}" for a in STEP)
            + " -- the inertance the row charged to itself was the term that was missing")
    share = {a: value(E(a, 1, 0), "q") / Q_IN for a in ARMS}
    inert_share = m_by / (m_x1 + m_by)
    require(all(abs(share[a] / inert_share - 1.0) < 0.09
                for a in ("lo_steady", "parab_steady")),
            f"the first accepted level is close to the pure inertance division "
            f"M_By/(M_X1+M_By) = {inert_share:.6f} only where the partner is "
            "quasi-steady: "
            + " ".join(f"{a} {share[a]:.6f}"
                       for a in ("lo_steady", "parab_steady")))
    require(all(-0.195 < share[a] / inert_share - 1.0 < -0.13
                for a in STEP),
            "and it is systematically BELOW it once the block carries its own du/dt: "
            + " ".join(f"{a} {share[a] / inert_share - 1:+.2%}" for a in STEP)
            + f" against {inert_share:.6f} -- at dt = {DT} s, short of these arms' own "
            "tau, the region's inertia is not the only mass the first kick has to "
            "move, so the bypass takes a larger cut than the row's inertance claim "
            "would give it")
    dest = {a: value(E(a, N_STEPS, 0), "q") / Q_IN for a in ARMS}
    require(dest["parab_steady"] > dest["lo_steady"]
            and dest["parab_euler"] > dest["p17_euler"] > dest["lo_euler"],
            "the share the ramp relaxes toward keeps the shape's order inside each "
            f"formulation: parab {dest['parab_steady']:.6f} > plug "
            f"{dest['lo_steady']:.6f} quasi-steady, and parab "
            f"{dest['parab_euler']:.6f} > 1/7 {dest['p17_euler']:.6f} > plug "
            f"{dest['lo_euler']:.6f} stepped")
    require(dest["parab_euler"] < dest["parab_steady"]
            and dest["lo_euler"] < dest["lo_steady"],
            "and every stepped arm is still travelling at the last level, behind its "
            f"quasi-steady twin: parab {dest['parab_euler']:.6f} < "
            f"{dest['parab_steady']:.6f}, plug {dest['lo_euler']:.6f} < "
            f"{dest['lo_steady']:.6f} -- 1.6 s is roughly three tau, not the end of "
            "the ramp")
    shift = {s: (value(E("parab_" + s, 1, 0), "q")
                 / value(E("lo_" + s, 1, 0), "q") - 1.0)
             for s in ("steady", "euler")}
    form_shift = (value(E("lo_euler", 1, 0), "q") / value(E("lo_steady", 1, 0), "q")
                  - 1.0, value(E("parab_euler", 1, 0), "q")
                  / value(E("parab_steady", 1, 0), "q") - 1.0)
    require(abs(form_shift[0]) > abs(shift["steady"])
            and abs(form_shift[1]) > abs(shift["euler"])
            and all(v < -0.15 for v in form_shift)
            and all(0.0 < v < 0.12 for v in shift.values()),
            f"at that first level the partner's formulation moves the parallel rate "
            f"more than the shape does: du/dt costs plug {form_shift[0]:+.1%} and "
            f"parabola {form_shift[1]:+.1%} while a parabola inside one formulation "
            f"is worth {shift['steady']:+.1%} quasi-steady and "
            f"{shift['euler']:+.1%} stepped -- the double count is already "
            "steering, the shape is not yet")

    # --------------------- 9. the pinned edge, where only the shape is allowed to move
    for a in ARMS:
        for n in range(1, N_STEPS + 1):
            require(close(value(E(a, n, 1), "q"), Q_IN, 1e-11),
                    f"{a} level {n}: the series edge carries the whole injected rate "
                    "however its inlet is shaped and however its block steps")
        dq2 = [abs(value(E(a, n, 1), "dq")) for n in range(2, N_STEPS + 1)]
        require(max(dq2) < 1e-18,
                f"{a}: and from level 2 on its rate change is {max(dq2):.2e} m3/s -- "
                "machine zero, so any inertia it reports was bought by the profile, "
                "not by the flow")
    dpI2 = {a: value(E(a, 2, 1), "dpInertia") for a in STEP}
    require(all(7.0 < v < 9.0 for v in dpI2.values()),
            "a region can accelerate with no rate change at all: at level 2 the "
            "pinned edge reports "
            + " ".join(f"{a} {v:.3f} Pa" for a, v in dpI2.items())
            + " of rho*du/dt over the probe slab while dq is 7e-21 -- the mass "
            "weighted speed moved because the profile did")
    require(all(value(B(a, "X2"), "inertanceDynamic") > 1e10 for a in STEP),
            f"and its dynamic inertance at the last solve is "
            + " ".join(f"{a} {value(B(a, 'X2'), 'inertanceDynamic'):.1e}" for a in STEP)
            + " -- rate over zero, which is why this checker asserts no criterion on "
            "it: the number is only meaningful on an edge whose rate moved")
    for f in ("lo", "parab"):
        worst = max(abs(value(E(f + "_euler", n, 1), "dpBlock")
                        / value(E(f + "_steady", n, 1), "dpBlock") - 1.0)
                    for n in range(3, N_STEPS + 1))
        require(worst < 0.0025,
                f"{f}: on the pinned edge the two formulations end within "
                f"{worst * 100:.2f}% of each other's drop from level 3 on -- the "
                "acceleration term is real and it is negligible there, because the "
                "rate that drives it is zero")
    dx1 = {f: (value(E(f + "_euler", N_STEPS, 0), "q")
               / value(E(f + "_steady", N_STEPS, 0), "q") - 1.0)
           for f in ("lo", "parab")}
    require(dx1["lo"] < -0.01 and dx1["parab"] < -0.01,
            "on the parallel edge the same switch moves the rate "
            f"plug {dx1['lo'] * 100:.2f}% and parabola {dx1['parab'] * 100:.2f}% "
            "at the last level -- where the rate is free the inertial double "
            "count is not negligible, it is the trajectory")
    require(dx1["parab"] / dx1["lo"] > 4.0,
            "and it is the shaped arm that feels it hardest, by "
            f"{dx1['parab'] / dx1['lo']:.2f}x: {dx1['parab'] * 100:.2f}% against "
            f"{dx1['lo'] * 100:.2f}% -- the shape is paid on the chord, and the "
            "chord is what sets the rate")
    dr8 = {a: value(B(a, "X2"), "dpPa") for a in ARMS}
    require(dr8["lo_euler"] > dr8["p17_euler"] > dr8["parab_euler"]
            and dr8["lo_steady"] > dr8["parab_steady"],
            f"X2's drop at the end keeps the shape's order in both formulations: "
            f"plug {dr8['lo_steady']:.2f} / {dr8['lo_euler']:.2f} Pa, 1/7 "
            f"{dr8['p17_euler']:.2f}, parabola {dr8['parab_steady']:.2f} / "
            f"{dr8['parab_euler']:.2f} -- the shape is "
            f"{(1 - dr8['parab_euler'] / dr8['lo_euler']) * 100:.1f}% of a stepped "
            "edge's drop and it is paid entirely in the friction half")

    # --------------------- 10. the schedule, the cost, and the walls
    pass_tot = {}
    for a in ARMS:
        passes = sum(len(v) for v in by_arm[a]["pass"].values())
        pass_tot[a] = passes
        s = by_arm[a]["summary"]
        require(value(s, "steps") == N_STEPS
                and close(value(s, "t"), N_STEPS * DT, 1e-12),
                f"{a}: {value(s, 'steps'):.0f} levels to t = {value(s, 't'):.4f} s")
        require(value(s, "residual") < TOL_INTERFACE,
                f"{a}: the interface closed, worst residual "
                f"{value(s, 'residual'):.3e} inside {TOL_INTERFACE}")
        require(by_arm[a]["solves"] == 2 * passes,
                f"{a}: {by_arm[a]['solves']} SIMPLE solves for {passes} passes -- one "
                "network solve and one solve per block, per pass, per level, whether "
                "or not the block carries du/dt")
        for n in range(1, N_STEPS + 1):
            p = by_arm[a]["pass"][n]
            last = p[max(p)]
            require(1 <= len(p) <= MAX_PASSES,
                    f"{a} level {n}: {len(p)} passes, inside the {MAX_PASSES} the "
                    "case allows")
            res = pass_nums(last, "residuals")
            rate = pass_nums(last, "rates")
            require(len(res) == 2 and len(rate) == 2
                    and all(close(rate[i], value(E(a, n, i), "q"), 1e-12)
                            and abs(res[i] - value(E(a, n, i), "residual")) < 1e-15
                            and max(res) < TOL_INTERFACE for i in (0, 1)),
                    f"{a} level {n}: the level's state IS its last pass")
            worst_i = int(re.search(r" worst (\d+)", last).group(1))
            require(worst_i in (0, 1) and abs(res[worst_i] - max(res)) < 1e-15,
                    f"{a} level {n}: the stop is taken on edge {worst_i}, whose "
                    "residual is the largest in the same line")
        require(all(value(by_arm[a]["step"][n], "newton") < 1e-11
                    for n in range(1, N_STEPS + 1)),
                f"{a}: the network solve itself is at machine precision at every "
                "level -- the interface, not the linear algebra, is what iterates")
    for f in ("lo", "parab"):
        require(pass_tot[f + "_euler"] > pass_tot[f + "_steady"],
            f"du/dt is not free: {f} goes from {pass_tot[f + '_steady']} exchanges "
            f"with a quasi-steady partner to {pass_tot[f + '_euler']} with its own "
            f"time level ({pass_tot[f + '_euler'] / pass_tot[f + '_steady']:.2f}x), "
            "because the coefficient the interface negotiates on now wanders with "
            "the level")
    require(pass_tot["parab_euler"] > pass_tot["lo_euler"]
            and pass_tot["p17_euler"] > pass_tot["lo_euler"],
            f"within the stepped family the shaped arms cost more: lo "
            f"{pass_tot['lo_euler']}, p17 {pass_tot['p17_euler']}, parab "
            f"{pass_tot['parab_euler']} passes")
    capped = {(a, n): sum(1 for k in range(1, N_STEPS + 1)
                          if value(E(a, k, NAME.index(n)), "sweeps") == MAX_SWEEPS)
              for a in ARMS for n in NAME}
    unconv = {(a, n): sum(1 for k in range(1, N_STEPS + 1)
                          if value(E(a, k, NAME.index(n)), "blockConverged") == 0)
              for a in ARMS for n in NAME}
    require(capped[("parab_steady", "X2")] == 3
            and unconv[("parab_steady", "X2")] == 3
            and capped[("parab_euler", "X2")] == 1
            and unconv[("parab_euler", "X2")] == 1,
            f"the honest wall, and it is the same event as blockConverged 0: the "
            f"parabolic pinned block spends all {MAX_SWEEPS:.0f} sweeps on "
            f"{capped[('parab_steady', 'X2')]}/8 levels with a quasi-steady partner "
            f"and {capped[('parab_euler', 'X2')]}/8 with its own time level")
    require(all(capped[(a, n)] == 0 and unconv[(a, n)] == 0
                for a in ARMS for n in NAME
                if (a, n) not in (("parab_steady", "X2"), ("parab_euler", "X2"))),
            f"no other edge comes close: a plug or 1/7 inlet on the same pinned rate "
            f"resolves in {value(E('lo_euler', N_STEPS, 1), 'sweeps'):.0f} sweeps, "
            f"and the shaped parallel edge in "
            f"{value(E('parab_euler', N_STEPS, 0), 'sweeps'):.0f}")
    require(max(value(E(a, k, i), "sweeps") for a in ("lo_euler", "p17_euler")
                for k in range(1, N_STEPS + 1) for i in (0, 1)) <= 17,
            "stepping costs the parallel block nothing in sweeps: the worst of the "
            f"two flat/1-7 Euler arms is "
            f"{max(value(E(a, k, i), 'sweeps') for a in ('lo_euler', 'p17_euler') for k in range(1, N_STEPS + 1) for i in (0, 1)):.0f}"
            " -- the extra exchanges are the interface's, not the block's")

    # --------------------- 11. the partner stays legal, and the network closes
    yp = {(a, n): value(B(a, n), "wallYPlusMax") for a in ARMS for n in NAME}
    require(all(20.0 < v < 35.0 for v in yp.values()),
            "the first cell centres stay in the log-law band under both partners "
            + " ".join(f"{a}/{n} {v:.1f}" for (a, n), v in sorted(yp.items())))
    nr = {(a, n): value(B(a, n), "nutMax") / NU for a in ARMS for n in NAME}
    require(all(40.0 < v < 150.0 for v in nr.values()),
            "the eddy viscosity runs "
            + " ".join(f"{a}/{n} {v:.0f}x" for (a, n), v in sorted(nr.items()))
            + " the molecular one, so the wall functions are the closure being used")
    for a in ARMS:
        for n in NAME:
            gap_, span_, _, length = DUCT[n]
            dh = 2.0 * gap_ * span_ / (gap_ + span_)
            plug = value(B(a, n), "plugSpeed")
            require(RHO * plug * dh / MU > 4000.0,
                    f"{a} {n}: Re_Dh {RHO * plug * dh / MU:.0f} at the rate the ramp "
                    "ended on")
            require(close(value(B(a, n), "inletArea"), gap_ * span_, 1e-12),
                    f"{a} {n}: the inlet patch is the section the shape was sampled on")
            require(abs(value(B(a, n), "fluxOut") / value(B(a, n), "fluxIn") - 1.0)
                    < 1e-10,
                    f"{a} {n}: what enters the region leaves it to "
                    f"{abs(value(B(a, n), 'fluxOut') / value(B(a, n), 'fluxIn') - 1):.1e}"
                    " -- the shape is redistributed across the section, and stepping "
                    "it does not create mass")
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
            if a in STEP and NAME.index(n) == 0:
                require(abs(credited - plain) < 0.11 and abs(credited) < 0.20,
                        f"{a} {n}: crediting the entry edge's own rho*du/dt "
                        f"({acc:+.1f} Pa/m) moves the window balance by "
                        f"{(credited - plain) * 100:+.2f} pp to "
                        f"{credited * 100:+.2f}%, from {plain * 100:+.2f}% -- the term "
                        "is real, it is a few percent of the shear, and it does not "
                        "rescue or wreck the balance either way")
            if a in SHAPED and NAME.index(n) == 0:
                require(abs(credited) < abs(plain),
                        f"{a} {n}: for the arms handed a parabola the credit runs "
                        f"toward zero, {plain * 100:+.2f}% -> "
                        f"{credited * 100:+.2f}%, because their parallel rate is still "
                        f"climbing at the last level ({acc:+.1f} Pa/m) -- the window "
                        "balance is a momentum balance and the acceleration belongs in "
                        "it")
            if a in STEP - SHAPED and NAME.index(n) == 0:
                require(abs(plain) < 0.04,
                        f"{a} {n}: the plug-fed stepped arm needs no credit to sit "
                        f"inside 4% ({plain * 100:+.2f}%, {credited * 100:+.2f}% with "
                        "it) -- the unbalanced half of the parabola's window is the "
                        "shape's own doing, not the time level's")
    lead = pipeChord(*LEAD, Q_IN) * Q_IN
    for a in ARMS:
        require(close(nd(a, N_STEPS, "Src") - nd(a, N_STEPS, "J1"), lead, 1e-9)
                and close(nd(a, N_STEPS, "J3") - nd(a, N_STEPS, "Snk"), lead, 1e-9),
                f"{a}: the two lead pipes are the network's own Churchill drop, "
                f"{lead:.6f} Pa at the injected rate, and no shape or ddt enters it")
        require(close(nd(a, N_STEPS, "Snk"), P_REF, 1e-12),
                f"{a}: Snk is held at the case's fixed pressure")
        require(all(close(nd(a, n, "J1") - nd(a, n, "J2"),
                          value(E(a, n, 0), "rowDrop"), 1e-9)
                    and close(nd(a, n, "J2") - nd(a, n, "J3"),
                              value(E(a, n, 1), "rowDrop"), 1e-9)
                    for n in range(1, N_STEPS + 1)),
                f"{a}: every level's node pair is every edge's momentum row, "
                "inertial term and all")
        require(all(abs(pipeq(a, n, "A") - Q_IN) < 1e-18
                    and close(pipeq(a, n, "X1") + pipeq(a, n, "By"), Q_IN, 1e-11)
                    for n in range(1, N_STEPS + 1)),
                f"{a}: A takes the injected rate and the pair splits it at every level")

    print("note: rowError*dpBlock (Pa) at levels 1/4/8, X1 -- " + "  ".join(
        f"{a} {value(E(a, 1, 0), 'rowError') * value(E(a, 1, 0), 'dpBlock'):.2f}/"
        f"{value(E(a, 4, 0), 'rowError') * value(E(a, 4, 0), 'dpBlock'):.2f}/"
        f"{value(E(a, N_STEPS, 0), 'rowError') * value(E(a, N_STEPS, 0), 'dpBlock'):.2f}"
        for a in STEP))
    print("note: chi (measured/probe inertance) level 1 X1 " + "  ".join(
        f"{a} {chi[(a, 1, 0)]:.4f}" for a in STEP)
        + "  |  dpFriction parab/lo " + "  ".join(
            f"L{n} {value(E('parab_euler', n, 0), 'dpFriction') / value(E('lo_euler', n, 0), 'dpFriction') - 1:+.1%}"
            for n in (1, N_STEPS)))
    print("note: tau " + "  ".join(f"{a} {taum[a]:.4f} s (secant {form[a]:.4f})"
                                    for a in ARMS))
    print("note: beta_out level 1 -> 8 " + "  ".join(
        f"{a}/{n} {outs[(a, 1, NAME.index(n))]:.6f}->{outs[(a, N_STEPS, NAME.index(n))]:.6f}"
        for a in STEP for n in NAME))
    print("note: pRisePa level 1 -> 8 on X1 " + "  ".join(
        f"{a} {rise[(a, 1, 0)]:.3f}->{rise[(a, N_STEPS, 0)]:.3f} Pa at "
        f"{value(E(a, 1, 0), 'pRiseX'):.3f}->{value(E(a, N_STEPS, 0), 'pRiseX'):.3f} m"
        for a in ("parab_steady", "parab_euler")))
    print("note: interface passes per arm " + " ".join(
        f"{a} {pass_tot[a]}" for a in ARMS))

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: a shape and a block that answers du/dt together -- the acceleration "
          "is blind to the shape and the friction is not, the recovery stays inside "
          "the friction half the row negotiates, and the row's inertance now bills "
          "the same term twice where the textbook formula finally adds up")


if __name__ == "__main__":
    sys.exit(main())
