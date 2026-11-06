#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient-multi and verify a time-stepping network
with TWO detailed regions inside it.

The case is the steady multi-edge network driven by a rate switched on at t=0:
Src -> A -> { X1 (3D duct) || By (1D bypass) } -> X2 (3D duct) -> B -> Snk. The
two edges are asked to agree in different ways, and that pairing is what makes
this one loop rather than two loops that happen to share a binary:

  * X1 is a parallel branch, so its rate is a result of the momentum balance -
    it starts at the inertial split and ramps to the resistance split with
    tau = (M_X1 + M_By)/(R_X1 + R_By).
  * X2 is in series, so an incompressible network hands it the whole injected
    rate from the first step: its interface negotiates the drop only, while its
    own inertance pays for the sudden start inside the momentum row.

Ten families of statements:

  1. one state per edge     - every time level reports both edges, each with the
                              rate its own pipe carries, R times that rate as the
                              drop the network charges, and the residual as the
                              gap between the two solvers
  2. the 1D network         - the series edges carry exactly the injected rate, the
                              branch split adds up, and every node drop closes as
                              R*Q + M*dQ/dt with THAT edge's installed inertance,
                              the series one including its sudden-start kick
  3. the interface schedule - the pass a step accepted has every edge under the
                              tolerance, the pass reports as its worst edge the one
                              whose residual is largest, and the passes a step spent
                              are the ones the slowest edge needed
  4. one pass, one solve    - the blocks were solved exactly passes x edges times,
                              so the edge that agreed first is still refreshed with
                              its neighbour rather than left on a half state
  5. edges do not lean      - the two sets seed X1 three times apart and X2
                              identically, so X2's first answer is the same digits in
                              both while X1's coefficient starts elsewhere
  6. the settled state      - both sets land on one per-edge rate, drop and
                              resistance, and that resistance is the duct's own
                              planar Poiseuille value to within its discretisation
  7. the first step         - the parallel split is by inertance alone
  8. the ramp               - tau = (M_X1 + M_By)/(R_X1 + R_By), and the whole
                              measured trajectory follows the discrete inertial solve
                              built from it, step by step
  9. each block's field     - conservation, one rate per plane, developed gradient
                              and parabolic profile, in each block's own duct
 10. laminar               - both ducts stayed laminar, which is what makes the
                              drops linear in the rate and family 6 valid

usage: check_transient_coupled_multi.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# geometry and fluid of tests/cases/transientCoupledMulti1d3d/case.dict
MU, RHO = 0.1, 800.0
Q_IN = 3.0e-5                  # m3/s injected at Src at t = 0
P_REF = 1.0e6                  # Pa held at Snk
DT, N_STEPS = 0.05, 12         # time { dt 0.05; endTime 0.6; }
MAX_PASSES = 8                 # couple{ maxSubIterations }
TOL_INTERFACE = 1e-4           # couple{ tolerance }

EDGES = {
    # edge: (length, gap, width) of the duct the 3D block stands in for
    "X1": (0.5, 0.02, 0.02),
    "X2": (0.4, 0.01, 0.03),
}
BY = (1.5, 0.03)               # the 1D bypass, sharing X1's node pair
LEAD = (5.0, 0.1)              # A and B, carrying everything
PIPE_D = {"X1": 0.02, "X2": 0.015}   # the coupled pipe's bore, for rho*L/A and Re
NAME = ("X1", "X2")


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


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("transient multi-edge coupled solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    def check(cond, msg):
        if not cond:
            require(False, msg)

    # ---- split the log into the starting sets ----
    sets = []
    cur = None
    for ln in lines:
        if ln.startswith("transient start ") and " inertia " in ln:
            cur = {"tag": ln, "seed": {}, "step": {}, "pass": {}, "edge": {},
                   "nodes": {}, "pipes": {}, "block": [], "summary": [], "solves": 0}
            sets.append(cur)
        elif cur is None:
            continue
        if ln.startswith("SIMPLE iter 1:"):
            cur["solves"] += 1
        elif ln.startswith("transient start ") and " index " in ln:
            cur["seed"][ln.split()[4]] = ln
        elif ln.startswith("transient step ") and " pass " in ln:
            n, p = int(ln.split()[2]), int(ln.split()[4])
            cur["pass"].setdefault(n, {})[p] = ln
        elif ln.startswith("transient step ") and " edge " in ln:
            cur["edge"].setdefault(int(ln.split()[2]), {})[int(ln.split()[4])] = ln
        elif ln.startswith("transient step "):
            cur["step"][int(ln.split()[2].rstrip(":"))] = ln
        elif ln.startswith("nodes step "):
            cur["nodes"][int(ln.split()[2].rstrip(":"))] = ln
        elif ln.startswith("pipes step "):
            cur["pipes"][int(ln.split()[2].rstrip(":"))] = ln
        elif ln.startswith("transient block "):
            cur["block"].append(ln)
        elif ln.startswith("transient coupled start ") and " name " in ln:
            cur["summary"].append(ln)

    header = [ln for ln in lines if ln.startswith("transient network+blocks")]
    require(len(header) == 1 and "2 coupled edges, 2 starting sets" in header[0]
            and "region inertia network" in header[0],
            "one network, two coupled edges, two starting sets and one inertia"
            " source, announced once")
    require(len(sets) == 2, f"the case ran both starting sets (got {len(sets)})")
    require(all(sorted(c["seed"]) == list(NAME) for c in sets),
            "each starting set quotes a coefficient for both edges at once")
    require(all("accelerate no" in c["tag"] for c in sets),
            "both sets run through the plain interface step - what is under test is two"
            " edges in one loop, not the accelerator")

    seed = [{k: value(c["seed"][k], "resistance") for k in NAME} for c in sets]

    for a, c in enumerate(sets):
        require(sorted(c["step"]) == list(range(1, N_STEPS + 1)),
                f"set {a}: all {N_STEPS} steps reported a state (got {len(c['step'])})")
        require(sorted(c["nodes"]) == sorted(c["step"])
                and sorted(c["pipes"]) == sorted(c["step"]),
                f"set {a}: every step reports its whole 1D state")
        require(len(c["block"]) == 2 and len(c["summary"]) == 2,
                f"set {a}: both blocks and both per-edge interface states reported once")

    # ---- (1) one state per edge, at every time level ----
    q = [{}, {}]      # per set, per edge: the rate the step delivers through it
    R = [{}, {}]      # the coefficient the network solved that step with
    M = [{}, {}]      # the inertance installed in that edge's momentum row
    dpB = [{}, {}]    # and the drop its block reported
    for a, c in enumerate(sets):
        for n in sorted(c["step"]):
            s = c["step"][n]
            check(abs(value(s, "t") - n * DT) <= 1e-15,
                  f"set {a} step {n}: t is {n} steps of {DT} s")
            check(int(value(s, "passes")) == len(c["pass"][n])
                  == int(value(s, "subPasses")),
                  f"set {a} step {n}: {len(c['pass'][n])} passes recorded, which is what"
                  " the step reports it spent")
            check(value(s, "newton") < 1e-11,
                  f"set {a} step {n}: the 1D system solved to its own floor")
            check(sorted(c["edge"][n]) == [0, 1],
                  f"set {a} step {n}: both edges report, not only the one that moved")
            for i, name in enumerate(NAME):
                e = c["edge"][n][i]
                check(e.split()[6] == name, f"set {a} step {n}: edge {i} is {name}")
                x, res = value(e, "q"), value(e, "residual")
                check(abs(x - value(c["pipes"][n], name)) <= 1e-12 * Q_IN,
                      f"set {a} step {n} {name}: the reported rate is that pipe's flow")
                check(abs(value(e, "dpNetwork") - value(e, "R") * x)
                      <= 1e-12 * abs(value(e, "dpNetwork")),
                      f"set {a} step {n} {name}: the network charges its coefficient"
                      " times the rate that step delivers")
                check(abs(res - abs(value(e, "dpNetwork") - value(e, "dpBlock"))
                          / value(e, "dpBlock")) <= 1e-9,
                      f"set {a} step {n} {name}: the residual is the gap between two"
                      " solvers, not a number the case was given")
                check(res <= TOL_INTERFACE,
                      f"set {a} step {n} {name}: interface residual {res:.3e} under the"
                      f" {TOL_INTERFACE:g} tolerance")
                check(value(e, "blockConverged") == 1,
                      f"set {a} step {n} {name}: its block converged at the state the"
                      " step accepted")
                q[a].setdefault(i, []).append(x)
                R[a].setdefault(i, []).append(value(e, "R"))
                M[a].setdefault(i, []).append(value(e, "M"))
                dpB[a].setdefault(i, []).append(value(e, "dpBlock"))

    # while the case says `inertias ( network )` the row carries the bore's rho*L/A
    for i, name in enumerate(NAME):
        L = EDGES[name][0]
        want = RHO * L / (0.25 * math.pi * PIPE_D[name] ** 2)
        check(all(abs(m / want - 1.0) <= 1e-12 for m in M[a][i] for a in (0, 1)),
              f"{name}: every time level carries the coupled pipe's rho*L/A ="
              f" {want:.6g} kg/m4")

    # ---- (2) the network closes with each edge's own inertance ----
    R_BY = hagenPoiseuille(*BY)
    R_LEAD = hagenPoiseuille(*LEAD)
    M_BY = RHO * BY[0] / (0.25 * math.pi * BY[1] ** 2)
    M_LEAD = RHO * LEAD[0] / (0.25 * math.pi * LEAD[1] ** 2)
    worst = [{k: 0.0 for k in ("X1", "X2", "By", "A", "B")} for _ in sets]
    seriesKick = []
    for a, c in enumerate(sets):
        prev = {}
        for n in sorted(c["step"]):
            nd, pp = c["nodes"][n], c["pipes"][n]
            p = {k: value(nd, k) for k in ("Src", "J1", "J2", "J3", "Snk")}
            f = {k: value(pp, k) for k in ("A", "X1", "By", "X2", "B")}
            x1, x2 = q[a][0][n - 1], q[a][1][n - 1]

            check(abs(f["A"] - Q_IN) <= 1e-12 * Q_IN,
                  f"set {a} step {n}: the lead carries the whole injected rate")
            check(abs(f["X2"] - Q_IN) <= 1e-12 * Q_IN,
                  f"set {a} step {n}: the series region takes exactly the injected rate"
                  " - topology decides it, the interface does not")
            check(abs(f["B"] - Q_IN) <= 1e-12 * Q_IN,
                  f"set {a} step {n}: nothing is lost between the two regions")
            check(abs(f["A"] - x1 - f["By"]) <= 1e-12 * Q_IN
                  and abs(f["X2"] - x1 - f["By"]) <= 1e-12 * Q_IN,
                  f"set {a} step {n}: the branch split adds up at both ends of X1")
            check(abs(f["X1"] - x1) <= 1e-12 * Q_IN
                  and abs(f["X2"] - x2) <= 1e-12 * Q_IN,
                  f"set {a} step {n}: each block is solved at the rate its own edge"
                  " carries")
            check(abs(p["Snk"] - P_REF) <= 1e-9,
                  f"set {a} step {n}: the sink holds its reference pressure")

            for key, drop, flow, resist, inert in (
                    ("X1", p["J1"] - p["J2"], x1, R[a][0][n - 1], M[a][0][n - 1]),
                    ("X2", p["J2"] - p["J3"], x2, R[a][1][n - 1], M[a][1][n - 1]),
                    ("By", p["J1"] - p["J2"], f["By"], R_BY, M_BY),
                    ("A", p["Src"] - p["J1"], f["A"], R_LEAD, M_LEAD),
                    ("B", p["J3"] - p["Snk"], f["B"], R_LEAD, M_LEAD)):
                dev = abs(drop - (resist * flow
                                  + inert * (flow - prev.get(key, 0.0)) / DT)) / abs(drop)
                worst[a][key] = max(worst[a][key], dev)
            prev.update(f)
            if n == 1:
                seriesKick.append((p["J2"] - p["J3"], R[a][1][0] * x2
                                   + M[a][1][0] * x2 / DT, M[a][1][0] * x2 / DT))

    for a, w in enumerate(worst):
        require(w["X1"] <= 1e-8 and w["X2"] <= 1e-8,
                f"set {a}: both detailed regions' node drops close as R*Q + M*dQ/dt with"
                f" the inertance installed in that edge (worst {w['X1']:.2e}/"
                f"{w['X2']:.2e})")
        require(w["By"] <= 1e-9 and w["A"] <= 1e-9 and w["B"] <= 1e-9,
                f"set {a}: the 1D edges close the same way with their own inertia"
                f" (worst {w['By']:.1e}/{w['A']:.1e}/{w['B']:.1e})")

    for a, (drop, want, kick) in enumerate(seriesKick):
        require(abs(drop - want) <= 1e-8 * want,
                f"set {a}: the series edge has nothing to negotiate about its rate and"
                f" still charges an inertial first drop - {drop:.4f} Pa against R*q +"
                f" M*q/dt = {want:.4f} Pa, of which {kick:.1f} Pa is the sudden start")

    # ---- (3) the schedule each step ran through ----
    lone = 0
    for a, c in enumerate(sets):
        over = 0
        for n in sorted(c["pass"]):
            passes = c["pass"][n]
            last = passes[max(passes)]
            rel = values(last, "residuals")
            check(max(rel) <= TOL_INTERFACE,
                  f"set {a} step {n}: the accepted pass has BOTH edges under the"
                  f" tolerance (worst {max(rel):.3e})")
            check(int(value(last, "worst")) == rel.index(max(rel)),
                  f"set {a} step {n}: the pass reports as its worst edge the one whose"
                  " residual really is largest")
            check(1 <= len(passes) <= MAX_PASSES,
                  f"set {a} step {n}: agreed inside the sub-iteration budget")
            for k in sorted(passes)[:-1]:
                r = values(passes[k], "residuals")
                if max(r) > TOL_INTERFACE:
                    over += 1
                if min(r) <= TOL_INTERFACE:
                    lone += 1
        nPass = sum(len(c["pass"][n]) for n in c["pass"])
        require(over == nPass - N_STEPS,
                f"set {a}: {over} of the {nPass} passes in the run ended over tolerance"
                f" and were re-run, and exactly one per step ({N_STEPS}) was accepted"
                " - the step stops when the WORST edge agrees, not the first one")
    require(lone >= 1,
            f"{lone} of those passes were spent on one edge while its neighbour was"
            " already under the tolerance: the schedule is one number for the whole"
            " network, so an edge that agrees early is refreshed with the one that"
            " does not")

    # ---- (4) one pass = one network solve + one solve per block ----
    for a, c in enumerate(sets):
        nPass = sum(int(value(c["step"][n], "passes")) for n in c["step"])
        require(c["solves"] == nPass * 2,
                f"set {a}: {c['solves']} block solves for {nPass} passes x 2 edges, so"
                " a second region is refreshed by every pass rather than run in its own"
                " loop")

    # ---- (5) an edge that starts far off does not move its neighbour ----
    require(abs(seed[1]["X1"] / seed[0]["X1"] - 3.0) <= 1e-12
            and seed[0]["X2"] == seed[1]["X2"],
            f"the two sets differ by one coefficient: X1 starts {seed[0]['X1']:.3g}"
            f" then {seed[1]['X1']:.3g}, X2 at {seed[0]['X2']:.3g} in both")
    first = [values(c["pass"][1][1], "coeff") for c in sets]
    require(abs(first[1][0] / first[0][0] - 3.0) <= 1e-12,
            f"the first pass of each set is still its own seed: X1's coefficient starts"
            f" {first[1][0] / first[0][0]:.1f}x apart")
    require(abs(first[1][1] / first[0][1] - 1.0) <= 1e-12,
            "while X2 starts on the same number in both and is fed the same rate")
    e2 = [value(c["edge"][1][1], "R") for c in sets]
    require(abs(e2[0] / e2[1] - 1.0) <= 1e-5,
            f"so the series edge answers as good as identically however its neighbour is"
            f" seeded: step 1 resistance {e2[0]:.10g} and {e2[1]:.10g}"
            f" ({abs(e2[0] / e2[1] - 1):.1e} apart - what is left is the number of"
            " passes its own block was warm-started through, not the rate it was fed")
    e1 = [value(c["edge"][1][0], "q") for c in sets]
    require(abs(e1[0] / e1[1] - 1.0) <= 1e-6,
            f"and the parallel edge still delivers the same agreed rate on the same"
            f" step ({e1[0]:.6g} and {e1[1]:.6g} m3/s)")
    p1 = [values(c["pass"][1][1], "residuals") for c in sets]
    lo, hi = min(p1[1]), max(p1[1])
    require(hi > TOL_INTERFACE and lo < 0.1 * hi,
            f"in the second set X1 was already agreed to {lo:.3e} while X2 was still"
            f" {hi:.3e} off on the same first pass - one edge cannot report the"
            " interface state on its own")

    # ---- (6) both sets land on one state ----
    for i, name in enumerate(NAME):
        for what, hist in (("rate", q), ("resistance", R), ("drop", dpB)):
            x, y = hist[0][i][-1], hist[1][i][-1]
            require(abs(x / y - 1.0) <= 2e-4,
                    f"{name} ends on one {what} whichever edge was seeded far:"
                    f" {x:.10g} and {y:.10g} ({abs(x / y - 1) * 100:.4f}%)")
    duct = [planarDuct(*EDGES[n]) for n in NAME]
    print(f"settled resistances: X1 {R[0][0][-1]:.6g} (duct {duct[0]:.6g},"
          f" {R[0][0][-1] / duct[0]:.4f}x), X2 {R[0][1][-1]:.6g} (duct {duct[1]:.6g},"
          f" {R[0][1][-1] / duct[1]:.4f}x)")
    for i, name in enumerate(NAME):
        require(abs(R[0][i][-1] / duct[i] - 1.0) <= 0.03,
                f"{name}: the block's chord resistance is its duct's planar Poiseuille"
                f" value to within {abs(R[0][i][-1] / duct[i] - 1) * 100:.2f}% - the"
                " deficit is entrance length and discretisation, not coupling")

    # ---- (7) the first step divides the sudden rate by inertance ----
    x1, r1 = q[0][0][0], R[0][0][0]
    q1 = Q_IN * (R_BY + M_BY / DT) / ((M[0][0][0] + M_BY) / DT + r1 + R_BY)
    require(abs(x1 / q1 - 1.0) <= 1e-3,
            f"the first step splits by inertance: X1 gets {x1:.6g} against"
            f" (R_By + M_By/dt)*Q/((M_X1 + M_By)/dt + R_X1 + R_By) = {q1:.6g} m3/s"
            f" ({abs(x1 / q1 - 1) * 100:.3f}%)")
    require(x1 < 0.95 * q[0][0][-1],
            f"and the split really ramps: X1 starts {x1 / q[0][0][-1]:.3f} of the way to"
            " its resistance share")

    # ---- (8) the ramp is the discrete inertial solve, step by step ----
    xs, ys = q[0][0][1:-1], q[0][0][2:]
    k = len(xs)
    sx, sy = sum(xs), sum(ys)
    B = (k * sum(u * v for u, v in zip(xs, ys)) - sx * sy) / (
        k * sum(u * u for u in xs) - sx * sx)
    A = (sy - B * sx) / k
    qInfFit = A / (1.0 - B)
    tau = (M[0][0][-1] + M_BY) / (R[0][0][-1] + R_BY)
    qInf = Q_IN * R_BY / (R[0][0][-1] + R_BY)
    decay = 1.0 / (1.0 + DT / tau)
    tauDisc = -DT / math.log(B)
    print(f"ramp: tau {tau:.6g} s closed form, decay per step {decay:.6f} predicted /"
          f" {B:.6f} fitted, that decay's own time constant {tauDisc:.6g} s; q_inf"
          f" {qInf:.6g} predicted / {qInfFit:.6g} m3/s fitted")
    require(0.5 < B < 1.0,
            f"the gap to the resistance share decays geometrically, {1.0 / (1.0 - B):.1f}"
            f" steps per decay - {N_STEPS} steps run {N_STEPS * DT / tau:.2f} of the"
            " closed-form time constants")
    require(abs(B / decay - 1.0) <= 1e-2,
            f"the decay per step is the discrete inertial solve's own"
            f" 1/(1 + dt/tau) = {decay:.6f}, measured {B:.6f}"
            f" ({(B / decay - 1) * 100:+.3f}%): one implicit Euler step of {DT} s against"
            f" a {tau:.3g} s ramp is dt/tau = {DT / tau:.2f} of a time constant, so the"
            f" scheme's effective constant is the longer {tauDisc:.6g} s"
            f" ({tauDisc / tau - 1:+.2%}) rather than the continuous {tau:.6g}")
    require(abs(qInfFit / qInf - 1.0) <= 5e-3,
            f"the ramp extrapolates to the resistance share Q*R_By/(R_X1+R_By) ="
            f" {qInf:.6g} m3/s, fitted destination {qInfFit:.6g}"
            f" ({abs(qInfFit / qInf - 1) * 100:.3f}%)")
    dev = max(abs(q[a][0][n] - (qInf + (q[a][0][0] - qInf) * decay ** n)) / qInf
              for a in (0, 1) for n in range(N_STEPS))
    require(dev <= 3e-3,
            f"both sets follow that discrete ramp from their own first value, every step"
            f" and within {dev * 100:.3f}% of the operating point")
    require(abs(q[0][1][-1] / Q_IN - 1.0) <= 1e-12,
            "and the series edge never moved: its rate is decided by topology, so all"
            " its interface ever negotiated was its drop")

    # ---- (9) each block's own field ----
    nu = MU / RHO
    for a, c in enumerate(sets):
        for i, name in enumerate(NAME):
            fin = c["block"][i]
            L, H, W = EDGES[name]
            area = H * W
            f = value(fin, "fluxIn")
            plug = value(fin, "plugSpeed")
            require(abs(f / q[a][i][-1] - 1.0) <= 1e-6,
                    f"set {a} {name}: the block carries the rate its edge ended with"
                    f" ({f:.6g} m3/s)")
            require(abs(value(fin, "inletArea") - area) <= 1e-12 * area,
                    f"set {a} {name}: the inlet patch is the section the inertance was"
                    f" referred to ({value(fin, 'inletArea'):.6g} m2)")
            require(abs(value(fin, "fluxOut") - f) <= 1e-7 * f,
                    f"set {a} {name}: what enters the region leaves it")
            require(value(fin, "planeFluxMax") - value(fin, "planeFluxMin")
                    <= 1e-7 * f,
                    f"set {a} {name}: every cross-section plane inside carries the same"
                    " rate")
            gradExact = 12.0 * MU * plug / (H * H)
            require(abs(abs(value(fin, "probeSlope")) - gradExact) <= 0.03 * gradExact,
                    f"set {a} {name}: developed dp/dx"
                    f" {abs(value(fin, 'probeSlope')):.4g} against planar Poiseuille"
                    f" {gradExact:.4g}")
            require(abs(value(fin, "uxCentre") / plug - 1.5) <= 0.03,
                    f"set {a} {name}: centreline speed is 3/2 the mean")
            require(abs(value(fin, "inertance") - RHO * L / area)
                    <= 1e-9 * RHO * L / area,
                    f"set {a} {name}: the block measures rho*V/A^2 ="
                    f" {RHO * L / area:.6g} kg/m4 on its own mesh"
                    f" ({value(fin, 'inertance'):.10g})")
            require(abs(dpB[a][i][-1] - value(fin, "dpPa")) <= 1e-12,
                    f"set {a} {name}: the drop the interface ended with is the drop the"
                    " block's own field measures")
            if a == 0:
                require(value(fin, "nutMax") == 0 and value(fin, "wallYPlusMax") == 0,
                        f"{name}: stayed laminar (Re {plug * PIPE_D[name] / nu:.1f} on"
                        f" the pipe bore, {plug * H / nu:.1f} on the gap), which is what"
                        " makes both drops linear in the rate")

    nPass = sum(int(value(c["step"][n], "passes")) for c in sets for n in c["step"])
    print(f"cost: {sum(c['solves'] for c in sets)} block solves and {nPass} interface"
          f" passes for {N_STEPS * len(sets)} agreed time levels on two edges")

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: two detailed regions agree with one network inside every time step, on"
          " one schedule that stops on the worst of them")


if __name__ == "__main__":
    sys.exit(main())
