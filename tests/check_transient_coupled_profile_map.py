#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient-multi and verify that the 1D<->3D interface
hands each detailed region an inlet SHAPE inside every time step, and that a
shape has two things to say once it is inside a ramp.

The case is the multi-edge network with a rate switched on at t=0, fed ten times
the un-shaped rate so that the entrance region sits inside each duct:
Src -> A -> { X1 (3D duct) || By (1D bypass) } -> X2 (3D duct) -> B -> Snk, with
five arms - `lo` flat on both edges (the control), `parab` and `p17` one fuller
profile on both, `mixed` parabolic on X1 and flat on X2, and `parab_qn` the
parabolic arm through the accelerated operator.

Why a ramp and not a steady solve: the resistance a shape settles is what the
parallel branch relaxes TOWARD, and it is one of the two numbers in the ramp's
own time constant tau = (M_X1 + M_By)/(R_X1 + R_By). So the same declaration
moves the destination and the speed, and only a case that keeps re-applying the
shape at every level can show either.

Ten families of statements:

  1. the declaration is per edge, per arm, per level - what the case asked for is
     what each edge reports at all twelve levels and what each block reports
  2. re-normalised every level - the section integral of the condition the inlet
     received is that level's rate, and the declared shape's own factor is a
     level-independent constant
  3. the beta ledger - the shape handed over, the shape the region's first plane
     carries and the shape its last one carries: developing arms climb, arms fed
     a full profile are edited downward, and the outlet is one value for all five
  4. the destination moves - the branch share of the switched-on rate lands where
     the shaped resistance says, while the series edge never moves its rate
  5. per-edge ownership - an edge's whole trajectory is decided by its own shape
     and by nothing its neighbour is given
  6. the 1D network still closes - every level, every arm, every edge's own
     inertance, the sudden-start kick included
  7. the interface schedule - one schedule for the whole network, accepted only
     when the worst edge agrees, one block solve per edge per pass
  8. the ramp is the discrete inertial solve - first value, decay per step and
     destination, arm by arm, and the time constants differ because the shapes do
  9. duct laws and each block's own field - what a shaped resistance is worth
     against planar Poiseuille, conservation, the developed gradient, and the
     scatter inside the window the drop is measured over
 10. laminar, and the momentum-flux ledger as measured - the imposed inlet flux is
     exactly rho*beta*q^2/A and is worth 1e-4 of the drop at this Reynolds number

usage: check_transient_coupled_profile_map.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# geometry and fluid of tests/cases/transientCoupledProfileMap/case.dict
MU, RHO = 0.1, 800.0
Q_IN = 3.0e-4                  # m3/s injected at Src at t = 0
P_REF = 1.0e6                  # Pa held at Snk
DT, N_STEPS = 0.05, 12         # time { dt 0.05; endTime 0.6; }
MAX_PASSES = 8                 # couple{ maxSubIterations }
TOL_INTERFACE = 1e-4           # couple{ tolerance }

ARMS = ("lo", "parab", "p17", "mixed", "parab_qn")
DECLARED = {
    "lo": {"X1": "flat", "X2": "flat"},
    "parab": {"X1": "parabolic", "X2": "parabolic"},
    "p17": {"X1": "power1over7", "X2": "power1over7"},
    "mixed": {"X1": "parabolic", "X2": "flat"},
    "parab_qn": {"X1": "parabolic", "X2": "parabolic"},
}
ACCELERATED = {"parab_qn"}
# what each declared shape is worth as a momentum-flux factor, integrated over
# the 20 wall-normal panels the inlet patch is made of
BETA_ANCHOR = {"flat": 1.0,
               "parabolic": 1.1970121617640874,
               "power1over7": 1.0134305810983724}
BETA_POISEUILLE = 6.0 / 5.0     # the continuum value the last plane approaches

EDGES = {"X1": (0.5, 0.02, 0.02), "X2": (0.4, 0.01, 0.03)}   # length, gap, width
BY = (1.5, 0.03)               # the 1D bypass, sharing X1's node pair
LEAD = (5.0, 0.1)              # A and B, carrying everything
PIPE_D = {"X1": 0.02, "X2": 0.015}
NAME = ("X1", "X2")


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


def values(line, key):
    """The run of numbers stored under `key`, stopping at the next word."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+((?:[-+0-9.eE]+\s*)+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return [float(v) for v in m.group(1).split()]


def planarDuct(length, gap, width):
    """Resistance of a plane channel: q = gap^3 width dp / (12 mu L)."""
    return 12.0 * MU * length / (gap ** 3 * width)


def hagenPoiseuille(length, diameter):
    return 128.0 * MU * length / (math.pi * diameter ** 4)


def stations(line):
    """The (x, beta) pairs of a `transient betaProfile` line."""
    tail = re.search(r": (.*)$", line).group(1).split()
    return [(float(tail[k]), float(tail[k + 1])) for k in range(0, len(tail) - 1, 2)]


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("transient profile-mapping solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    def check(cond, msg):
        if not cond:
            require(False, msg)

    # ---- split the log into the arms ----
    arms = []
    cur = None
    for ln in lines:
        if ln.startswith("transient start ") and " accelerate " in ln:
            cur = {"tag": ln, "arm": word(ln, "start").rstrip(":"), "shape": {},
                   "seed": {}, "step": {}, "pass": {}, "edge": {},
                   "nodes": {}, "pipes": {}, "block": {}, "beta": {},
                   "summary": [], "solves": 0}
            arms.append(cur)
        elif cur is None:
            continue
        if ln.startswith("SIMPLE iter 1:"):
            cur["solves"] += 1
        elif ln.startswith("transient arm ") and " profile " in ln:
            cur["shape"][word(ln, "edge")] = word(ln, "profile")
        elif ln.startswith("transient start ") and " resistance " in ln:
            cur["seed"][word(ln, "edge")] = value(ln, "resistance")
        elif ln.startswith("transient step ") and " pass " in ln:
            cur["pass"].setdefault(int(ln.split()[2]), {})[int(ln.split()[4])] = ln
        elif ln.startswith("transient step ") and " edge " in ln:
            cur["edge"].setdefault(int(ln.split()[2]), {})[int(ln.split()[4])] = ln
        elif ln.startswith("transient step "):
            cur["step"][int(ln.split()[2].rstrip(":"))] = ln
        elif ln.startswith("nodes step "):
            cur["nodes"][int(ln.split()[2].rstrip(":"))] = ln
        elif ln.startswith("pipes step "):
            cur["pipes"][int(ln.split()[2].rstrip(":"))] = ln
        elif ln.startswith("transient block "):
            cur["block"][word(ln, "name").rstrip(":")] = ln
        elif ln.startswith("transient betaProfile "):
            cur["beta"][word(ln, "edge").rstrip(":")] = ln
        elif ln.startswith("transient coupled start ") and " name " in ln:
            cur["summary"].append(ln)

    def E(a, n, i):
        """The interface state of edge `i` at level `n` of arm `a`."""
        return by_arm[a]["edge"][n][i]

    by_arm = {c["arm"]: c for c in arms}

    header = [ln for ln in lines if ln.startswith("transient network+blocks")]
    require(len(header) == 1 and "2 coupled edges, 5 starting sets" in header[0],
            "one network, two coupled edges and five arms announced once")
    require([c["arm"] for c in arms] == list(ARMS),
            f"the five arms ran, in the order the case declares them (got "
            f"{[c['arm'] for c in arms]})")
    require(all(sorted(c["seed"]) == list(NAME) for c in arms),
            "every arm seeds both edges, and with the same pair the un-shaped"
            " transient case quotes: 1.25e6 then 5.3e6")
    require(abs(arms[0]["seed"]["X1"] / arms[0]["seed"]["X2"]
                - 1.25e6 / 5.3e6) <= 1e-12,
            "the seeds are one third under each duct value, so no arm starts at"
            " its own resistance")
    accel = {c["arm"]: "accelerate yes" in c["tag"] for c in arms}
    require(accel == {a: (a in ACCELERATED) for a in ARMS},
            "only parab_qn runs through the accelerated operator; the other four"
            " are plain Gauss-Seidel, so a difference between them is a shape and"
            " not an integrator")

    # ---- (1) the declaration reaches every edge at every level ----
    for c in arms:
        require(c["shape"] == DECLARED[c["arm"]],
                f"{c['arm']}: the case's `profile` reached its two edges as declared"
                f" ({c['shape']['X1']} on X1, {c['shape']['X2']} on X2)")
    shape_levels = {(a, n, i): word(E(a, n, i), "inletShape")
                    for a in ARMS for n in range(1, N_STEPS + 1) for i in (0, 1)}
    require(all(shape_levels[k] == DECLARED[k[0]][NAME[k[2]]]
                for k in shape_levels),
            f"all {len(shape_levels)} arm x level x edge states report the shape that"
                " arm declared for that edge - the interface re-writes it at every"
                " pass of every time level, not once at the start of a run")
    for a in ARMS:
        for i, name in enumerate(NAME):
            b = by_arm[a]["block"][name]
            check(word(b, "inletShape") == DECLARED[a][name],
                  f"{a} {name}: the block's own last solve carries that shape")
            check(abs(value(b, "inletArea")
                      - EDGES[name][1] * EDGES[name][2]) <= 1e-15,
                  f"{a} {name}: the patch it was spread over is the section the"
                  f" inertance is referred to ({value(b, 'inletArea'):.6g} m2)")

    # ---- (2) what the level asked for is what the inlet received ----
    worstFlux = 0.0
    for a in ARMS:
        for n in range(1, N_STEPS + 1):
            for i, name in enumerate(NAME):
                e = E(a, n, i)
                q, fi = value(e, "q"), value(e, "fluxImposed")
                check(q > 0, f"{a} step {n} {name}: the level delivers a positive"
                             " rate to shape")
                worstFlux = max(worstFlux, abs(fi / q - 1.0))
    require(worstFlux <= 1e-12,
            f"at every one of the {N_STEPS * 2 * len(ARMS)} levels the section integral"
            f" of the condition the inlet got equals the rate that level solved for"
            f" (worst {worstFlux:.2e} relative) - a shape was traded for a shape, not"
            " for a different amount of fluid")
    for a in ARMS:
        for i, name in enumerate(NAME):
            shape = DECLARED[a][name]
            lev = [value(E(a, n, i), "betaImposed") for n in range(1, N_STEPS + 1)]
            anchor = BETA_ANCHOR[shape]
            require(max(abs(v / anchor - 1.0) for v in lev) <= 1e-12,
                    f"{a} {name}: every level reads back the factor of `{shape}`,"
                    f" {anchor:.16g}, on the faces the interface wrote (worst"
                    f" {max(abs(v / anchor - 1.0) for v in lev):.2e})")
            require(max(lev) - min(lev) <= 1e-12 * max(1.0, abs(anchor)),
                    f"{a} {name}: and it is the same number at level 1 and level"
                    f" {N_STEPS} (spread {max(lev) - min(lev):.2e}) - the shape is a"
                    " declaration, only its normalisation follows the rate")

    # ---- (3) the beta ledger: handed over vs carried ----
    for a in ARMS:
        for i, name in enumerate(NAME):
            e = E(a, N_STEPS, i)
            bi, bIn, bOut = (value(e, "betaImposed"), value(e, "betaIn"),
                             value(e, "betaOut"))
            if BETA_ANCHOR[DECLARED[a][name]] < bOut:
                require(bi < bIn < bOut,
                        f"{a} {name}: a {DECLARED[a][name]} inlet arrives at {bi:.6g},"
                        f" its first solved plane reads {bIn:.6g} and its last"
                        f" {bOut:.6g} - the profile handed over is not yet the region's"
                        " own, so the region still has to develop it")
            else:
                require(bIn > bOut and abs(bi / bIn - 1.0) < 0.01,
                        f"{a} {name}: a full profile arrives at {bi:.6g}, is read"
                        f" back at {bIn:.6g} on the first plane and settles to"
                        f" {bOut:.6g} at the last - the region is edited down"
                        " towards its own developed shape instead of climbing to it")
    for i, name in enumerate(NAME):
        outs = [value(E(a, N_STEPS, i), "betaOut") for a in ARMS]
        ins = [value(E(a, N_STEPS, i), "betaIn") for a in ARMS]
        require(max(outs) - min(outs) <= 1e-5,
                f"{name}: the last plane reads one shape factor whichever inlet was"
                f" handed it ({min(outs):.7g} to {max(outs):.7g}, spread"
                f" {max(outs) - min(outs):.2e}) - the shape is visible at the inlet"
                " and gone at the outlet")
        inSpread, outSpread = max(ins) - min(ins), max(outs) - min(outs)
        require(inSpread >= 0.05 and inSpread >= 1000 * outSpread,
                f"{name}: the first plane reads {min(ins):.6g} to {max(ins):.6g} -"
                f" {inSpread:.3f} apart, {inSpread / min(ins) * 100:.1f}% of the value"
                f" and {inSpread / outSpread:.0f}x the spread at the outlet - the same"
                " five arms are indistinguishable at one end of the region and obvious"
                " at the other")
        require(abs(outs[0] / BETA_POISEUILLE - 1.0) <= 0.01,
                f"{name}: that one value is {outs[0]:.6g} against the continuum"
                f" plane-Poiseuille factor 6/5 = {BETA_POISEUILLE:.6g}"
                f" ({abs(outs[0] / BETA_POISEUILLE - 1) * 100:.2f}% off, which is what"
                " a finite duct on this mesh leaves)")
    flatArm, parabArm = by_arm["lo"], by_arm["parab"]
    spread = {}
    for tag, c in (("lo", flatArm), ("parab", parabArm)):
        v = [value(c["edge"][n][0], "betaIn") for n in range(1, N_STEPS + 1)]
        spread[tag] = max(v) - min(v)
    require(spread["lo"] >= 20 * spread["parab"],
            f"the ledger is a per-level number, not a per-run constant: on the flat"
            f" arm X1's first-plane factor moves {spread['lo']:.4g} between level 1"
            f" and level {N_STEPS} as the ramp carries more rate through a fixed"
            f" geometry, while on the parabolic arm it moves {spread['parab']:.4g}"
            f" ({spread['lo'] / spread['parab']:.0f}x less) - a profile that arrives"
            " already shaped has nothing left to develop")
    for a in ARMS:
        for i, name in enumerate(NAME):
            st = stations(by_arm[a]["beta"][name])
            betas = [b for _, b in st]
            first, last = st[0], st[-1]
            if BETA_ANCHOR[DECLARED[a][name]] < last[1]:
                require(min(betas) == first[1] and max(betas) == last[1]
                        and sum(1 for u, v in zip(betas, betas[1:]) if v > u)
                        >= 0.8 * (len(betas) - 1),
                        f"{a} {name}: {len(st)} cross-sections, the factor climbing"
                        f" from {first[1]:.6g} at x={first[0]:.4g} to {last[1]:.6g} at"
                        f" x={last[0]:.4g} in almost every step - a single imposed rate"
                        " could not say where in the duct that happens")
            else:
                argmax = betas.index(max(betas))
                argmin = betas.index(min(betas))
                require(argmax == 0 and 0 < argmin < len(betas) - 1,
                        f"{a} {name}: its highest factor is the first plane"
                        f" ({max(betas):.6g} at x={first[0]:.4g}) and its lowest is"
                        f" inside the duct ({min(betas):.6g} at x={st[argmin][0]:.4g}),"
                        " the profile handed over being relaxed rather than grown")

    # ---- (4) the destination a shape settles to ----
    qEnd = {a: value(E(a, N_STEPS, 0), "q") for a in ARMS}
    duct = {n: planarDuct(*EDGES[n]) for n in NAME}
    R_BY = hagenPoiseuille(*BY)
    for a in ("parab", "p17"):
        shift = qEnd[a] / qEnd["lo"] - 1.0
        require(0.0 <= shift <= 0.04,
                f"{a}: the branch share the ramp ends at is {shift * 100:+.4f}% against"
                f" the flat control ({qEnd[a]:.10g} vs {qEnd['lo']:.10g} m3/s) - the"
                " fuller inlet costs the region less, so the bypass yields flow")
    require(qEnd["lo"] < qEnd["p17"] < qEnd["parab"],
            f"and the three declared shapes order themselves by how full they are:"
            f" flat {qEnd['lo']:.6g}, 1/7 power {qEnd['p17']:.6g}, parabolic"
            f" {qEnd['parab']:.6g} m3/s - the same rate in three distributions is"
            " three fixed points, not one fixed point described three ways")
    for a in ARMS:
        for n in range(1, N_STEPS + 1):
            check(abs(value(E(a, n, 1), "q") / Q_IN - 1.0) <= 1e-12,
                  f"{a} step {n}: X2 carries exactly the injected rate")
    R2 = {a: value(E(a, N_STEPS, 1), "R") for a in ARMS}
    require(abs(R2["parab"] / R2["lo"] - 1.0) <= 0.07
            and abs(R2["p17"] / R2["lo"] - 1.0) <= 0.02
            and abs(R2["mixed"] / R2["lo"] - 1.0) <= 1e-4,
            f"while the series edge keeps its rate and changes its drop:"
            f" X2's resistance ends at {R2['lo']:.6g} flat, {R2['p17']:.6g} under the"
            f" 1/7-power inlet ({(R2['p17'] / R2['lo'] - 1) * 100:+.2f}%) and"
            f" {R2['parab']:.6g} under the parabolic one"
            f" ({(R2['parab'] / R2['lo'] - 1) * 100:+.2f}%), with `mixed` back on the"
            " flat value")

    # ---- (5) an edge answers only its own declaration ----
    for n in range(1, N_STEPS + 1):
        for key in ("q", "R", "dpBlock", "betaIn", "betaOut", "sweeps"):
            check(value(E("mixed", n, 0), key) == value(E("parab", n, 0), key),
                  f"level {n} {key}: mixed's X1 is parab's X1 to the digit")
    require(True, "mixed's X1 equals the both-parabolic arm's X1 on every key at every"
                  " one of the 12 levels, digit for digit - X1 and By form a current"
                  " divider whose split no downstream region can move")
    moved = max(abs(value(E("mixed", n, 0), "q") / value(E("lo", n, 0), "q") - 1.0)
                for n in range(1, N_STEPS + 1))
    require(moved >= 0.01,
            f"and that identity is not vacuous: the same arm's X1 is"
            f" {moved * 100:.2f}% away from the flat control's at the same levels")
    drift2 = max(abs(value(E("mixed", n, 1), "R") / value(E("lo", n, 1), "R") - 1.0)
                 for n in range(1, N_STEPS + 1))
    away2 = max(abs(value(E("mixed", n, 1), "R") / value(E("parab", n, 1), "R") - 1.0)
                for n in range(1, N_STEPS + 1))
    require(drift2 <= 5e-6 and away2 >= 0.05,
            f"mixed's X2 follows the FLAT arm to {drift2:.2e} relative while sitting"
            f" {away2 * 100:.1f}% away from the arm that shaped it like its own X1 -"
            " each edge takes its own shape. Not to the digit, because the interface"
            " stops at a tolerance and each level's warm start depends on how many"
            " passes the run has spent")

    # ---- (6) the 1D network closes around both shapes ----
    M_BY = RHO * BY[0] / (0.25 * math.pi * BY[1] ** 2)
    R_LEAD = hagenPoiseuille(*LEAD)
    M_LEAD = RHO * LEAD[0] / (0.25 * math.pi * LEAD[1] ** 2)
    worst = {a: {k: 0.0 for k in ("X1", "X2", "By", "A", "B")} for a in ARMS}
    for a in ARMS:
        prev = {}
        for n in range(1, N_STEPS + 1):
            nd, pp = by_arm[a]["nodes"][n], by_arm[a]["pipes"][n]
            p = {k: value(nd, k) for k in ("Src", "J1", "J2", "J3", "Snk")}
            f = {k: value(pp, k) for k in ("A", "X1", "By", "X2", "B")}
            check(abs(f["A"] - Q_IN) <= 1e-12 * Q_IN
                  and abs(f["B"] - Q_IN) <= 1e-12 * Q_IN,
                  f"{a} step {n}: the lead pipes carry the whole injected rate")
            check(abs(f["X2"] - Q_IN) <= 1e-12 * Q_IN,
                  f"{a} step {n}: the series region takes it exactly, shape or no shape")
            check(abs(f["A"] - f["X1"] - f["By"]) <= 1e-12 * Q_IN
                  and abs(f["X2"] - f["X1"] - f["By"]) <= 1e-12 * Q_IN,
                  f"{a} step {n}: the branch split adds up at both ends")
            check(abs(p["Snk"] - P_REF) <= 1e-9, f"{a} step {n}: the sink holds its"
                                                 " reference pressure")
            for key, drop, resist, inert in (
                    ("X1", p["J1"] - p["J2"], value(E(a, n, 0), "R"),
                     value(E(a, n, 0), "M")),
                    ("X2", p["J2"] - p["J3"], value(E(a, n, 1), "R"),
                     value(E(a, n, 1), "M")),
                    ("By", p["J1"] - p["J2"], R_BY, M_BY),
                    ("A", p["Src"] - p["J1"], R_LEAD, M_LEAD),
                    ("B", p["J3"] - p["Snk"], R_LEAD, M_LEAD)):
                dev = abs(drop - (resist * f[key]
                                  + inert * (f[key] - prev.get(key, 0.0)) / DT)) / abs(drop)
                worst[a][key] = max(worst[a][key], dev)
            prev.update(f)
    for a in ARMS:
        require(max(worst[a].values()) <= 1e-8,
                f"{a}: all five pipes' node drops close as R*q + M*dq/dt at every"
                f" level, with each edge's own inertance (worst {max(worst[a].values()):.2e})"
                " - declaring a shape did not move the momentum balance the network"
                " solves")
    for i, name in enumerate(NAME):
        want = RHO * EDGES[name][0] / (0.25 * math.pi * PIPE_D[name] ** 2)
        require(all(abs(value(E(a, n, i), "M") / want - 1.0) <= 1e-9
                    for a in ARMS for n in range(1, N_STEPS + 1)),
                f"{name}: every level of every arm carries the coupled pipe's"
                f" rho*L/A = {want:.6g} kg/m4 - the inertance stays where"
                " `inertias ( network )` puts it, whatever the inlet looks like")
    for a in ARMS:
        nd = by_arm[a]["nodes"][1]
        p = {k: value(nd, k) for k in ("J2", "J3")}
        check(abs(p["J2"] - p["J3"] - value(E(a, 1, 1), "rowDrop"))
              <= 1e-8 * abs(p["J2"] - p["J3"]),
              f"{a} step 1: the series edge charges the sudden start it never negotiates"
              f" - {p['J2'] - p['J3']:.4f} Pa against the row's R*q + M*q/dt ="
              f" {value(E(a, 1, 1), 'rowDrop'):.4f} Pa")
        check(value(E(a, N_STEPS, 1), "unsteady") == 0
              and value(E(a, N_STEPS, 0), "unsteady") == 0,
              f"{a}: both partners are quasi-steady blocks here (`ddt steady`), so the"
              " interface has no du/dt half of the drop to split off and `rowError`"
              " grades nothing - this case is about the shape, not about the region's"
              " own history")

    # ---- (7) one schedule for the whole network ----
    totalPass = {}
    for a in ARMS:
        c = by_arm[a]
        nPass = 0
        for n in sorted(c["pass"]):
            passes = c["pass"][n]
            last = passes[max(passes)]
            rel = values(last, "residuals")
            check(max(rel) <= TOL_INTERFACE,
                  f"{a} step {n}: the accepted pass has BOTH edges under the tolerance"
                  f" (worst {max(rel):.3e})")
            check(int(value(last, "worst")) == rel.index(max(rel)),
                  f"{a} step {n}: the pass reports as its worst edge the one whose"
                  " residual really is largest")
            check(1 <= len(passes) <= MAX_PASSES,
                  f"{a} step {n}: agreed inside the sub-iteration budget")
            check(int(value(c["step"][n], "passes")) == len(passes)
                  == int(value(c["step"][n], "subPasses")),
                  f"{a} step {n}: {len(passes)} passes recorded, which is what the step"
                  " reports it spent")
            for i, name in enumerate(NAME):
                e = c["edge"][n][i]
                x, res = value(e, "q"), value(e, "residual")
                check(abs(res - abs(value(e, "dpNetwork") - value(e, "dpBlock"))
                          / value(e, "dpBlock")) <= 1e-9,
                      f"{a} step {n} {name}: the residual is the gap between the two"
                          " solvers' drops at the rate this level delivered")
                check(value(e, "blockConverged") == 1,
                      f"{a} step {n} {name}: its block converged at the state the level"
                          " accepted")
            nPass += len(passes)
        totalPass[a] = nPass
        require(c["solves"] == nPass * 2,
                f"{a}: {c['solves']} block solves for {nPass} passes x 2 edges - every"
                " pass refreshes both regions, including the one that already agreed")
    require(totalPass["lo"] == 25 and totalPass["parab"] == 17
            and totalPass["p17"] == 25 and totalPass["mixed"] == 17
            and totalPass["parab_qn"] == 17,
            f"the shaped arms agreed in {totalPass['parab']} interface passes over the"
            f" 12 levels against the flat control's {totalPass['lo']} - the shape"
            " changes the schedule as well as the answer, and `mixed` pays exactly what"
            " `parab` pays because its own agreement is dominated by the same shaped"
            " duct")
    require(totalPass["parab_qn"] == totalPass["parab"],
            f"the accelerated parabolic arm spent the same {totalPass['parab']} passes as"
            " the plain one: at this tolerance the shaped interface was already agreeing"
            " in one or two passes per level, so the operator had nothing to buy here"
            " and only the destination it did not touch")

    # ---- (8) the ramp is the discrete inertial solve, arm by arm ----
    M_X1 = RHO * EDGES["X1"][0] / (0.25 * math.pi * PIPE_D["X1"] ** 2)
    ramp = {}
    for a in ARMS:
        q = [value(E(a, n, 0), "q") for n in range(1, N_STEPS + 1)]
        r = [value(E(a, n, 0), "R") for n in range(1, N_STEPS + 1)]
        q1 = (Q_IN * (R_BY + M_BY / DT)
              / ((M_X1 + M_BY) / DT + r[0] + R_BY))
        require(abs(q[0] / q1 - 1.0) <= 1e-9,
                f"{a}: the first level splits the switched-on rate by inertance alone"
                f" - {q[0]:.10g} m3/s against"
                f" (R_By + M_By/dt)*Q/((M_X1 + M_By)/dt + R_X1 + R_By) = {q1:.10g}")
        u, v = q[1:-1], q[2:]
        k = len(u)
        su, sv = sum(u), sum(v)
        B = (k * sum(x * y for x, y in zip(u, v)) - su * sv) / (
            k * sum(x * x for x in u) - su * su)
        destination = (sv - B * su) / k / (1.0 - B)
        tau = (M_X1 + M_BY) / (r[-1] + R_BY)
        decay = 1.0 / (1.0 + DT / tau)
        qInf = Q_IN * R_BY / (r[-1] + R_BY)
        dev = max(abs(q[n] - (qInf + (q[0] - qInf) * decay ** n)) / qInf
                  for n in range(N_STEPS))
        require(0.5 < B < 1.0,
                f"{a}: the gap to the resistance share decays geometrically"
                f" (fitted ratio {B:.6f}, {1.0 / (1.0 - B):.1f} steps per decay)")
        require(abs(B / decay - 1.0) <= 6e-3,
                f"{a}: the measured decay per step is that arm's own discrete inertial"
                f" solve 1/(1 + dt/tau) = {decay:.6f}, measured {B:.6f}"
                f" ({(B / decay - 1) * 100:+.3f}%), from tau ="
                f" (M_X1 + M_By)/(R_X1 + R_By) = {tau:.6g} s built on the resistance"
                f" its inlet shape settled ({r[-1]:.6g})")
        require(abs(destination / qInf - 1.0) <= 1.5e-3,
                f"{a}: and it extrapolates to the resistance share"
                f" Q*R_By/(R_X1 + R_By) = {qInf:.10g} m3/s, fitted destination"
                f" {destination:.10g} ({abs(destination / qInf - 1) * 100:.4f}%)")
        require(dev <= 1.5e-3,
                f"{a}: all {N_STEPS} levels follow that discrete ramp from their own"
                f" first value, worst {dev * 100:.4f}% of the operating point")
        ramp[a] = (tau, B, qInf, q)
    tauLo, tauPar = ramp["lo"][0], ramp["parab"][0]
    require(1.01 < tauPar / tauLo <= 1.06,
            f"the shape is inside the time constant: X1's resistance ends at"
            f" {value(E('lo', N_STEPS, 0), 'R'):.6g} on the flat control and"
            f" {value(E('parab', N_STEPS, 0), 'R'):.6g} on the parabolic one, so the"
            f" ramp's own constant is {tauLo:.6g} s against {tauPar:.6g} s"
            f" ({(tauPar / tauLo - 1) * 100:+.2f}%) - the fuller inlet settles more"
            " slowly, and it also has further to go")
    require(ramp["parab"][3][-1] / ramp["parab"][2] < ramp["lo"][3][-1] / ramp["lo"][2],
            f"the parabolic arm ends {100 * (1 - ramp['parab'][3][-1] / ramp['parab'][2]):.2f}%"
            f" short of its destination where the flat one is"
            f" {100 * (1 - ramp['lo'][3][-1] / ramp['lo'][2]):.2f}% short - a shape"
            " declares both the answer and how long the network takes to get there")
    require(abs(ramp["parab_qn"][3][-1] / ramp["parab"][3][-1] - 1.0) <= 1e-6,
            f"parab_qn lands on parab's trajectory (last level"
            f" {ramp['parab_qn'][3][-1]:.10g} against {ramp['parab'][3][-1]:.10g} m3/s,"
            f" worst {max(abs(x / y - 1.0) for x, y in zip(ramp['parab_qn'][3], ramp['parab'][3])) * 100:.4f}%"
            " apart) - the operator moved nothing that the shape had decided")

    # ---- (9) duct laws, and each block's own field ----
    for a in ARMS:
        for i, name in enumerate(NAME):
            b = by_arm[a]["block"][name]
            L, H, W = EDGES[name]
            f = value(b, "fluxIn")
            plug = value(b, "plugSpeed")
            q = value(E(a, N_STEPS, i), "q")
            check(abs(f / q - 1.0) <= 1e-9,
                  f"{a} {name}: the region carries the rate its edge ended with"
                  f" ({f:.6g} m3/s)")
            check(abs(plug / (q / (H * W)) - 1.0) <= 1e-9,
                  f"{a} {name}: and its mean speed is that rate over the section the"
                  " shape was spread across")
            check(abs(value(b, "fluxOut") - f) <= 1e-7 * f,
                  f"{a} {name}: what enters the region leaves it")
            check(value(b, "planeFluxMax") - value(b, "planeFluxMin") <= 1e-7 * f,
                  f"{a} {name}: every cross-section inside carries the same rate, so a"
                  " shaped inlet did not buy mass")
            gradExact = 12.0 * MU * plug / (H * H)
            check(abs(abs(value(b, "probeSlope")) - gradExact) <= 0.04 * gradExact,
                  f"{a} {name}: the developed gradient in the window the drop is taken"
                  f" over is {abs(value(b, 'probeSlope')):.5g} Pa/m against planar"
                  f" Poiseuille {gradExact:.5g} ({abs(abs(value(b, 'probeSlope')) / gradExact - 1) * 100:.2f}%)")
            check(abs(value(b, "dpPa") - value(E(a, N_STEPS, i), "dpBlock")) <= 1e-12,
                  f"{a} {name}: the drop the interface ended with is the drop the"
                              " block's own field measures")
            check(abs(value(b, "momentumFluxPa")
                      - RHO * value(b, "betaImposed") * f * f / value(b, "inletArea"))
                  <= 1e-9 * value(b, "momentumFluxPa"),
                  f"{a} {name}: the inlet's own momentum flux is rho*beta*q^2/A ="
                  f" {RHO * value(b, 'betaImposed') * f * f / value(b, 'inletArea'):.6g} Pa")
    print(f"settled resistances against the duct law:"
          f" X1 flat {value(E('lo', N_STEPS, 0), 'R') / duct['X1']:.4f}x,"
          f" 1/7-power {value(E('p17', N_STEPS, 0), 'R') / duct['X1']:.4f}x,"
          f" parabolic {value(E('parab', N_STEPS, 0), 'R') / duct['X1']:.4f}x;"
          f" X2 flat {value(E('lo', N_STEPS, 1), 'R') / duct['X2']:.4f}x,"
          f" parabolic {value(E('parab', N_STEPS, 1), 'R') / duct['X2']:.4f}x")
    require(value(E("lo", N_STEPS, 0), "R") > value(E("p17", N_STEPS, 0), "R")
            > value(E("parab", N_STEPS, 0), "R"),
            "the resistance a region is worth falls as the inlet it is handed fills"
            " out - flat over the duct law (the plug still has to be developed and the"
            " chord pays for it), 1/7 power less so, parabolic under it")
    flatSpread = value(by_arm["lo"]["block"]["X1"], "probeSpread")
    parSpread = value(by_arm["parab"]["block"]["X1"], "probeSpread")
    require(parSpread <= 0.1 * flatSpread,
            f"X1: the gradient inside the window its drop is taken over is"
            f" {flatSpread:.4g} Pa wide on the flat arm and {parSpread:.4g} Pa on the"
            f" parabolic one ({flatSpread / parSpread:.0f}x flatter) - the drop a shaped"
            " inlet reports is a straighter line, which is what makes its resistance a"
            " number instead of an average over development")
    for name in NAME:
        for a in ("lo", "parab"):
            b = by_arm[a]["block"][name]
            dp = value(b, "dpPa")
            require(value(b, "probeSpread") <= 1e-3 * dp,
                    f"{a} {name}: its whole measurement window is"
                    f" {value(b, 'probeSpread') / dp * 100:.4f}% wide against its own"
                    f" drop ({dp:.5g} Pa), so the resistance the interface read is not"
                    " a chord across a curved line")
    for i, name in enumerate(NAME):
        lev = [value(E("lo", n, i), "unsteady") for n in range(1, N_STEPS + 1)]
        require(set(lev) == {0.0},
                f"{name}: no level of the control arm asked for an acceleration term"
                " in its drop, so the resistance carries the whole chord here")

    # ---- (10) laminar, and the ledger as measured ----
    for name in NAME:
        for a in ARMS:
            b = by_arm[a]["block"][name]
            require(value(b, "nutMax") == 0 and value(b, "wallYPlusMax") == 0,
                    f"{a} {name}: stayed laminar, which is what makes the shaped"
                    " resistances comparable to the duct law at all")
    nu = MU / RHO
    for name in NAME:
        b = by_arm["lo"]["block"][name]
        H = EDGES[name][1]
        print(f"{name} at Re {value(b, 'plugSpeed') * H / nu:.1f} on the gap,"
              f" {value(b, 'plugSpeed') * 2 * H / nu:.1f} on the hydraulic diameter")
    share = {}
    for a in ("lo", "parab"):
        for i, name in enumerate(NAME):
            b = by_arm[a]["block"][name]
            share[(a, name)] = (value(b, "momentumFluxPa"), value(b, "dpPa"))
    for (a, name), (mf, dp) in share.items():
        require(mf / dp <= 2e-4,
                f"{a} {name}: the momentum flux handed over is {mf:.4g} Pa worth against"
                f" the {dp:.5g} Pa the region reports ({mf / dp:.2e} of it) - the ledger"
                " is readable at every level, and at Re 160 it is not what the drop is"
                " made of")
    ratio = (share[("parab", "X2")][0] / share[("lo", "X2")][0])
    qRatio = value(E("parab", N_STEPS, 1), "q") / value(E("lo", N_STEPS, 1), "q")
    require(abs(ratio / (qRatio ** 2 * BETA_ANCHOR["parabolic"]) - 1.0) <= 1e-9,
            f"X2 carries the same rate in both arms, so the whole change in what its"
            f" inlet momentum flux is worth ({share[('parab', 'X2')][0]:.6g} against"
            f" {share[('lo', 'X2')][0]:.6g} Pa) is the declared shape: ratio {ratio:.7g}"
            f" = beta(parabolic) {BETA_ANCHOR['parabolic']:.7g} to 1e-9")

    nPass = sum(totalPass.values())
    nSolve = sum(c["solves"] for c in arms)
    print(f"cost: {nSolve} block solves and {nPass} interface passes for"
          f" {N_STEPS * len(ARMS)} agreed time levels on two edges and five arms;"
          f" the run is {len([ln for ln in lines if ln.startswith('transient step ') and ' subPasses ' in ln])}"
          f" reported levels long")

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: the interface hands each region its declared inlet shape at every"
          " time level, and inside a ramp that shape moves both the destination the"
          " branch relaxes to and the constant it relaxes with")


if __name__ == "__main__":
    sys.exit(main())
