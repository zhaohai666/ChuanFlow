#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient-multi with `inertias ( network partner )`
and verify that TWO detailed regions get their inertia from two different places
inside one time-stepping network.

The network is the multi-edge case's -- Src -> A -> { X1 (3D duct) || By (1D
bypass) } -> X2 (3D duct) -> B -> Snk, rate switched on at t = 0 -- and it is run
twice from the SAME seeds. The only difference between the arms is where the M in
the 1D momentum row comes from:

  * `bore`   keeps the coupled pipe's rho*L/A, i.e. the nominal bore of the round
             pipe standing in for a rectangular region,
  * `region` takes each block's own rho*V/Area^2 measured on the region's mesh.

What is under test is that this is a PER-EDGE number. The two regions are
rectangles with different aspect ratios, so the bore mis-states them by different
amounts (X1 0.02 x 0.02 -> 27.3%, X2 0.01 x 0.03 -> 69.8%), and the arms have to
disagree by exactly those two numbers, in one run, on the same schedule.

Nine families of statements:

  1. one difference       - the arms quote the same coefficient for both edges
                            and differ only in the attribution word
  2. where M comes from   - every time level of every edge carries the number that
                            arm's source gives, to 1e-12, and the two edges are
                            mis-stated by different amounts
  3. the block's own mesh - the measurement is identical in both arms because it is
                            a mesh property; installed == measured only where the
                            region was asked
  4. the mapping untouched- the coefficient the interface negotiates, and the number
                            of passes each step spends, are the same in both arms:
                            the inertance goes to the momentum row and nowhere else
  5. the momentum row     - each edge closes as R*q + M*dq/dt with its own arm's M,
                            including the series edge's sudden-start kick
  6. the first step       - the parallel split follows that arm's own inertance, and
                            the series edge's drop differs by exactly dM*Q/dt
  7. the ramp             - tau = (M_X1 + M_By)/(R_X1 + R_By) per arm, the ratio of
                            the two taus is the inertance ratio diluted by the
                            bypass, and each arm follows its own discrete decay
  8. the destination      - the arms end on one state, and its per-edge resistance
                            is the same duct deficit the single-attribution case
                            measured
  9. each block's field   - conservation, developed gradient, parabolic profile,
                            laminar in both ducts

usage: check_transient_coupled_multi_inertia.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# geometry and fluid of tests/cases/transientCoupledMultiInertia/case.dict
MU, RHO = 0.1, 800.0
Q_IN = 3.0e-5                  # m3/s injected at Src at t = 0
DT, N_STEPS = 0.05, 12         # time { dt 0.05; endTime 0.6; }
MAX_PASSES = 8                 # couple{ maxSubIterations }
TOL_INTERFACE = 1e-4           # couple{ tolerance }

ARMS = ("network", "partner")
NAME = ("X1", "X2")
EDGES = {"X1": (0.5, 0.02, 0.02), "X2": (0.4, 0.01, 0.03)}   # L, gap, width
PIPE_D = {"X1": 0.02, "X2": 0.015}
BY = (1.5, 0.03)               # the 1D bypass, sharing X1's node pair
LEAD = (5.0, 0.1)              # A and B, carrying everything
SEED = {"X1": 3.75e6, "X2": 1.2e7}   # the one starting set both arms share


def value(line, key):
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([-+0-9.eE]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return float(m.group(1))


def values(line, key):
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+((?:[-+0-9.eE]+\s*)+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return [float(v) for v in m.group(1).split()]


def planarDuct(length, gap, width):
    return 12.0 * MU * length / (gap ** 3 * width)


def hagenPoiseuille(length, diameter):
    return 128.0 * MU * length / (math.pi * diameter ** 4)


def boreInertance(length, diameter):
    """rho*L/A on the nominal round bore the 1D side was given."""
    return RHO * length / (0.25 * math.pi * diameter ** 2)


def regionInertance(length, gap, width):
    """rho*V/A^2 of the region; V = L*A, so this is rho*L/A_section."""
    return RHO * length / (gap * width)


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("transient multi-edge inertia-attribution solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    def check(cond, msg):
        if not cond:
            require(False, msg)

    # ---- split the log into the two attribution arms ----
    arms = []
    cur = None
    for ln in lines:
        if ln.startswith("transient start ") and " inertia " in ln:
            cur = {"tag": ln, "inertia": ln.split()[4], "seed": {}, "step": {},
                   "pass": {}, "edge": {}, "nodes": {}, "pipes": {}, "block": [],
                   "summary": [], "solves": 0}
            arms.append(cur)
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

    # ---- (1) the arms differ by attribution and nothing else ----
    header = [ln for ln in lines if ln.startswith("transient network+blocks")]
    require(len(header) == 1 and "2 coupled edges, 2 starting sets" in header[0]
            and "region inertia network (arms)" in header[0],
            "one network, two coupled edges, two starting sets and an attribution that"
            " is quoted per arm (the header names the first and says 'arms')")
    require(len(arms) == 2 and [c["inertia"] for c in arms] == list(ARMS),
            f"both arms ran, in the order the case asked for {list(ARMS)}"
            f" (got {[c['inertia'] for c in arms]})")
    require(all(sorted(c["seed"]) == list(NAME) for c in arms),
            "each arm quotes a coefficient for both edges at once")
    seeds = [{k: value(c["seed"][k], "resistance") for k in NAME} for c in arms]
    require(seeds[0] == seeds[1]
            and all(abs(seeds[a][k] - SEED[k]) <= 1e-12 for a in (0, 1) for k in NAME),
            f"the two arms start on the same two numbers (X1 {SEED['X1']:.3g},"
            f" X2 {SEED['X2']:.3g}) - the only thing the case changes between them is"
            " where M comes from")
    require(all("accelerate no" in c["tag"] for c in arms),
            "neither arm is accelerated: what is under test is the inertance's"
            " attribution, not the interface operator")

    for a, c in enumerate(arms):
        require(sorted(c["step"]) == list(range(1, N_STEPS + 1)),
                f"{ARMS[a]}: all {N_STEPS} steps reported a state")
        require(len(c["block"]) == 2 and len(c["summary"]) == 2,
                f"{ARMS[a]}: both blocks and both per-edge end states reported once")

    # ---- (2) where each edge's M comes from, per edge and per arm ----
    M = [{}, {}]
    R = [{}, {}]
    q = [{}, {}]
    dpB = [{}, {}]
    dpN = [{}, {}]
    for a, c in enumerate(arms):
        for i, name in enumerate(NAME):
            M[a][i] = [value(c["edge"][n][i], "M") for n in sorted(c["step"])]
            R[a][i] = [value(c["edge"][n][i], "R") for n in sorted(c["step"])]
            q[a][i] = [value(c["edge"][n][i], "q") for n in sorted(c["step"])]
            dpB[a][i] = [value(c["edge"][n][i], "dpBlock") for n in sorted(c["step"])]
            dpN[a][i] = [value(c["edge"][n][i], "dpNetwork") for n in sorted(c["step"])]

    for a in (0, 1):
        for i, name in enumerate(NAME):
            dev = max(abs(dpN[a][i][n] - R[a][i][n] * q[a][i][n]) for n in range(N_STEPS))
            require(dev <= 1e-9,
                    f"{name} / {ARMS[a]}: what the network charges its coupled edge is"
                    f" its own coefficient times the rate that step delivers (worst"
                    f" {dev:.1e} Pa over {N_STEPS} levels) - the interface state the"
                    " inertance was installed into is still self-consistent")

    for i, name in enumerate(NAME):
        L, H, W = EDGES[name]
        mb, mr = boreInertance(L, PIPE_D[name]), regionInertance(L, H, W)
        for a in (0, 1):
            w = mb if a == 0 else mr
            require(all(abs(m / w - 1.0) <= 1e-12 for m in M[a][i]),
                    f"{name} / {ARMS[a]}: all {N_STEPS} time levels carry"
                    f" {'rho*L/A on the bore' if a == 0 else 'rho*V/A^2 of the region'} ="
                    f" {w:.6g} kg/m4, one number for the whole arm")
        # the two edges are mis-stated by DIFFERENT amounts - the point of doing
        # this per edge instead of once for the network
        require(abs(mb / mr - (H * W) / (0.25 * math.pi * PIPE_D[name] ** 2))
                <= 1e-12,
                f"{name}: the bore over-states the region by {mb / mr:.6f}x, which is"
                f" the real section over the inscribed circle ({H * W:.6g} m2 against"
                f" pi*d^2/4 = {0.25 * math.pi * PIPE_D[name] ** 2:.6g} m2,"
                f" {0.25 * math.pi * PIPE_D[name] ** 2 / (H * W):.6f} of it)")
        print(f"{name}: bore {mb:.7f} vs region {mr:.7f} kg/m4,"
              f" {(mb / mr - 1) * 100:+.2f}%")
    span = [M[0][i][0] / M[1][i][0] for i in (0, 1)]
    require(abs(span[0] - span[1]) > 0.3,
            f"the two edges do not disagree by the same amount: the bore is"
            f" {(span[0] - 1) * 100:.1f}% high on X1 and {(span[1] - 1) * 100:.1f}% on"
            " X2, so one global correction factor could not have produced both")

    # ---- (3) the block's own measurement is a mesh property ----
    for i, name in enumerate(NAME):
        L, H, W = EDGES[name]
        measured = [value(c["block"][i], "inertance") for c in arms]
        mr = regionInertance(L, H, W)
        require(abs(measured[0] - measured[1]) <= 1e-12 * mr,
                f"{name}: both arms read the same {measured[0]:.7f} kg/m4 off the mesh"
                " - the measurement does not depend on what the arm then does with it")
        require(abs(measured[0] - mr) <= 1e-9 * mr,
                f"{name}: that number is rho*V/A^2 = {mr:.6g} kg/m4 for a"
                f" {H} x {W} m section {L} m long")
        require(abs(M[1][i][0] / measured[0] - 1.0) <= 1e-9,
                f"{name} / partner: what the row carries is what the block measured")
        require(abs(M[0][i][0] / measured[0] - boreInertance(L, PIPE_D[name]) / mr)
                <= 1e-9,
                f"{name} / network: what the row carries is the bore's, not the"
                f" region's - {M[0][i][0] / measured[0]:.4f}x the measurement")

    # ---- (4) attribution moves the momentum row, not the interface mapping ----
    for i, name in enumerate(NAME):
        worst = max(abs(R[0][i][n] / R[1][i][n] - 1.0) for n in range(N_STEPS))
        require(worst <= 2e-3,
                f"{name}: the coefficient the two arms negotiate is the same number"
                f" every step (worst {worst:.2e} apart, inside twice the"
                f" {TOL_INTERFACE:g} interface tolerance that bounds how far a step's R"
                " may sit from its block's dp/q) - the partner still reports only the"
                " drop it resists at the rate it is given, so the fixed point the"
                " interface iterates is unchanged by where M came from")
    passes = {a: [int(value(c["step"][n], "passes")) for n in sorted(c["step"])]
              for a, c in enumerate(arms)}
    require(max(passes[0] + passes[1]) <= 2,
            f"neither arm ever needed more than 2 of the {MAX_PASSES} passes a step may"
            f" spend (bore {passes[0]}, region {passes[1]}): the interface is hunting the"
            " same fixed point in both arms, so a different source for M buys it no new"
            " work - the passes are set by how far the partner's drop is from R*q")
    require(abs(sum(passes[0]) - sum(passes[1])) <= 4,
            f"and the two schedules stay in the same band even though the trajectories"
            f" do not: {sum(passes[0])} passes for the bore arm against"
            f" {sum(passes[1])} for the region arm"
            f" ({sum(passes[0]) - sum(passes[1])} apart over {N_STEPS} levels)")
    for a, c in enumerate(arms):
        for n in sorted(c["pass"]):
            last = c["pass"][n][max(c["pass"][n])]
            require(max(values(last, "residuals")) <= TOL_INTERFACE,
                    f"{ARMS[a]} step {n}: the accepted pass has both edges under the"
                    f" {TOL_INTERFACE:g} tolerance")
            require(1 <= len(c["pass"][n]) <= MAX_PASSES,
                    f"{ARMS[a]} step {n}: agreed inside the sub-iteration budget")
    require(all(c["solves"] == sum(int(value(c["step"][n], "passes"))
                                   for n in c["step"]) * 2 for c in arms),
            f"block solves are passes x 2 edges in both arms"
            f" ({arms[0]['solves']}/{arms[1]['solves']}), so the second region is still"
            " refreshed by every pass")

    # ---- (5) the momentum row closes with each arm's own M ----
    R_BY = hagenPoiseuille(*BY)
    M_BY = boreInertance(*BY)
    R_LEAD = hagenPoiseuille(*LEAD)
    M_LEAD = boreInertance(*LEAD)
    worst = [{k: 0.0 for k in ("X1", "X2", "By", "A", "B")} for _ in arms]
    seriesDrop = {}
    for a, c in enumerate(arms):
        prev = {}
        for n in sorted(c["step"]):
            p = {k: value(c["nodes"][n], k) for k in ("Src", "J1", "J2", "J3", "Snk")}
            f = {k: value(c["pipes"][n], k) for k in ("A", "X1", "By", "X2", "B")}
            if n == 1:
                seriesDrop[a] = p["J2"] - p["J3"]
            require(abs(f["X2"] - Q_IN) <= 1e-12 * Q_IN,
                    f"{ARMS[a]} step {n}: the series edge still takes exactly the"
                    " injected rate - attribution cannot change topology")
            require(abs(f["A"] - q[a][0][n - 1] - f["By"]) <= 1e-12 * Q_IN,
                    f"{ARMS[a]} step {n}: the branch split adds up")
            for key, drop, flow, resist, inert in (
                    ("X1", p["J1"] - p["J2"], q[a][0][n - 1], R[a][0][n - 1],
                     M[a][0][n - 1]),
                    ("X2", p["J2"] - p["J3"], q[a][1][n - 1], R[a][1][n - 1],
                     M[a][1][n - 1]),
                    ("By", p["J1"] - p["J2"], f["By"], R_BY, M_BY),
                    ("A", p["Src"] - p["J1"], f["A"], R_LEAD, M_LEAD),
                    ("B", p["J3"] - p["Snk"], f["B"], R_LEAD, M_LEAD)):
                dev = abs(drop - (resist * flow
                                  + inert * (flow - prev.get(key, 0.0)) / DT)) / abs(drop)
                worst[a][key] = max(worst[a][key], dev)
            prev.update(f)
    for a, w in enumerate(worst):
        require(w["X1"] <= 1e-8 and w["X2"] <= 1e-8,
                f"{ARMS[a]}: both regions' node drops close as R*q + M*dq/dt with that"
                f" arm's inertance (worst {w['X1']:.2e}/{w['X2']:.2e})")
        require(w["By"] <= 1e-9 and w["A"] <= 1e-9 and w["B"] <= 1e-9,
                f"{ARMS[a]}: the 1D edges close the same way with their own inertia"
                f" (worst {w['By']:.1e}/{w['A']:.1e}/{w['B']:.1e})")

    # ---- (6) the first step, both ways one number can enter the row ----
    for a in (0, 1):
        r1 = R[a][0][0]
        q1 = Q_IN * (R_BY + M_BY / DT) / ((M[a][0][0] + M_BY) / DT + r1 + R_BY)
        require(abs(q[a][0][0] / q1 - 1.0) <= 1e-8,
                f"{ARMS[a]} step 1: the switched-on rate divides by that arm's own"
                f" inertance - X1 gets {q[a][0][0]:.6g} against {q1:.6g} m3/s,"
                f" the same closed form read with a different M ({abs(q[a][0][0] / q1 - 1):.1e}"
                f" apart: the row's inertial term, not a tuned number)")
    require(q[1][0][0] > q[0][0][0],
            f"and the arms separate at once: the lighter region starts the parallel"
            f" branch higher, {q[1][0][0]:.6g} against {q[0][0][0]:.6g} m3/s"
            f" ({q[1][0][0] / q[0][0][0] - 1:+.2%})")
    dDrop = seriesDrop[0] - seriesDrop[1]
    dM = M[0][1][0] - M[1][1][0]
    want = dM * Q_IN / DT
    require(abs(dDrop - want) <= 5e-3,
            f"the series edge's rate never moved, so all its arms' difference is the"
            f" sudden start: first-step node drop {seriesDrop[0]:.4f} Pa (bore) against"
            f" {seriesDrop[1]:.4f} Pa (region), {(dDrop / want - 1) * 100:+.3f}% off the"
            f" arithmetic dM*Q/dt = {want:.4f} Pa")
    print(f"X2 first step: bore pays {M[0][1][0] * Q_IN / DT:.4f} Pa of sudden start,"
          f" region {M[1][1][0] * Q_IN / DT:.4f} Pa of the same {Q_IN:.1e} m3/s in"
          f" {DT} s")

    # ---- (7) each arm's ramp is its own, and the two taus ratio to the inertances ----
    inf = {}
    for a, c in enumerate(arms):
        xs, ys = q[a][0][1:-1], q[a][0][2:]
        k = len(xs)
        sx, sy = sum(xs), sum(ys)
        B = (k * sum(u * v for u, v in zip(xs, ys)) - sx * sy) / (
            k * sum(u * u for u in xs) - sx * sx)
        A = (sy - B * sx) / k
        inf[a] = (B, A / (1.0 - B))
        tau = (M[a][0][-1] + M_BY) / (R[a][0][-1] + R_BY)
        tauDisc = -DT / math.log(B)
        qInf = Q_IN * R_BY / (R[a][0][-1] + R_BY)
        require(0.5 < B < 1.0,
                f"{ARMS[a]}: the gap to the resistance share decays geometrically"
                f" ({B:.6f} per step, {tauDisc:.6g} s as a time constant against"
                f" dt/ln(1+dt/tau) = {DT / math.log(1 + DT / tau):.6g} s for that tau)")
        require(abs(B * (1 + DT / tau) - 1.0) <= 1e-2,
                f"{ARMS[a]}: the decay is 1/(1 + dt/tau) for THAT arm's tau ="
                f" {tau:.6g} s (1/(1+dt/tau) = {1 / (1 + DT / tau):.6f}, fitted {B:.6f}"
                f", {(B / (1 / (1 + DT / tau)) - 1) * 100:+.3f}%)")
        require(abs(inf[a][1] / qInf - 1.0) <= 5e-3,
                f"{ARMS[a]}: the ramp extrapolates to the resistance share"
                f" {qInf:.6g} m3/s, fitted destination {inf[a][1]:.6g}"
                f" ({abs(inf[a][1] / qInf - 1) * 100:.3f}%)")
        dev = max(abs(q[a][0][n] - (qInf + (q[a][0][0] - qInf)
                                    * (1.0 / (1.0 + DT / tau)) ** n)) / qInf
                  for n in range(N_STEPS))
        require(dev <= 3e-3,
                f"{ARMS[a]}: the whole {N_STEPS}-step trajectory follows that discrete"
                f" ramp within {dev * 100:.3f}% of the operating point")
    tau = {a: (M[a][0][-1] + M_BY) / (R[a][0][-1] + R_BY) for a in (0, 1)}
    L, H, W = EDGES["X1"]
    expect = (regionInertance(L, H, W) + M_BY) / (boreInertance(L, PIPE_D["X1"]) + M_BY)
    print(f"tau: bore {tau[0]:.6g} s, region {tau[1]:.6g} s, ratio"
          f" {tau[1] / tau[0]:.6f} (closed form {expect:.6f})")
    require(abs(tau[1] / tau[0] - expect) <= 2e-3,
            f"the two taus ratio to the two inertances, diluted by the bypass: X1's M is"
            f" {regionInertance(L, H, W) / boreInertance(L, PIPE_D['X1']):.4f}x the"
            f" bore's but tau only moves {(tau[1] / tau[0] - 1) * 100:.2f}%, because"
            f" M_By = {M_BY:.6g} kg/m4 is already {M_BY / (M_BY + boreInertance(L, PIPE_D['X1'])):.3f}"
            " of the branch's inertia")
    require(abs(inf[1][0] / inf[0][0] - 1.0) >= 1e-3,
            f"and the two arms really walk different ramps: decay per step {inf[0][0]:.6f}"
            f" (bore) against {inf[1][0]:.6f} (region)")

    # ---- (8) one destination, whichever source M came from ----
    # The arms are still 1/3 of a time constant short of that destination at t = 0.6,
    # so the END of the run is not where they have to agree - the extrapolation is.
    # What the last step does have to agree on is the two things that carry no memory
    # of the path: the coefficient, and the series edge's rate.
    x, y = inf[0][1], inf[1][1]
    require(abs(x / y - 1.0) <= 2e-4,
            f"X1's ramp extrapolates to one operating point whichever source its"
            f" inertance came from: {x:.10g} and {y:.10g} m3/s"
            f" ({abs(x / y - 1) * 100:.4f}%)")
    require(q[1][0][-1] > q[0][0][-1],
            f"while the paths stay apart at the end of the run, in the direction the"
            f" inertances say: X1 is at {q[0][0][-1]:.6g} (bore) and {q[1][0][-1]:.6g}"
            f" (region) m3/s, {(q[1][0][-1] / q[0][0][-1] - 1) * 100:.2f}% apart, because"
            f" the lighter region's ramp is the faster of the two")
    for i, name in enumerate(NAME):
        require(abs(R[0][i][-1] / R[1][i][-1] - 1.0) <= 2e-3,
                f"{name} ends on one resistance whichever source the inertance came"
                f" from: {R[0][i][-1]:.10g} and {R[1][i][-1]:.10g}")
    require(abs(q[0][1][-1] / q[1][1][-1] - 1.0) <= 1e-12,
            f"and the series edge was on one rate the whole time"
            f" ({q[0][1][-1]:.6g} m3/s in both): attribution cannot move a flow that"
            " topology already decided")
    duct = {name: planarDuct(*EDGES[name]) for name in NAME}
    for i, name in enumerate(NAME):
        for a in (0, 1):
            require(abs(R[a][i][-1] / duct[name] - 1.0) <= 0.03,
                    f"{name} / {ARMS[a]}: settled resistance is {R[a][i][-1]:.6g} against"
                    f" the duct's planar Poiseuille {duct[name]:.6g} Pa*s/m3, short by"
                    f" {abs(R[a][i][-1] / duct[name] - 1) * 100:.2f}% - the same deficit"
                    " the single-attribution case measured on these two blocks, so"
                    " attribution did not move the partition's answer")
    print(f"settled: X1 {R[0][0][-1] / duct['X1']:.4f}x and X2"
          f" {R[0][1][-1] / duct['X2']:.4f}x of planar Poiseuille in both arms")

    # ---- (9) each block's own field, in both arms ----
    nu = MU / RHO
    for a, c in enumerate(arms):
        for i, name in enumerate(NAME):
            fin = c["block"][i]
            L, H, W = EDGES[name]
            f = value(fin, "fluxIn")
            plug = value(fin, "plugSpeed")
            require(abs(f / q[a][i][-1] - 1.0) <= 1e-6,
                    f"{ARMS[a]} {name}: the block carries the rate its edge ended with")
            require(abs(value(fin, "inletArea") - H * W) <= 1e-12 * H * W,
                    f"{ARMS[a]} {name}: the inlet patch is the {H * W:.6g} m2 section the"
                    " arm's inertance was referred to")
            require(abs(value(fin, "fluxOut") - f) <= 1e-7 * f,
                    f"{ARMS[a]} {name}: what enters the region leaves it")
            gradExact = 12.0 * MU * plug / (H * H)
            require(abs(abs(value(fin, "probeSlope")) - gradExact) <= 0.03 * gradExact,
                    f"{ARMS[a]} {name}: developed dp/dx"
                    f" {abs(value(fin, 'probeSlope')):.4g} against planar Poiseuille"
                    f" {gradExact:.4g}")
            require(abs(value(fin, "uxCentre") / plug - 1.5) <= 0.03,
                    f"{ARMS[a]} {name}: centreline speed is 3/2 the mean")
            require(value(fin, "nutMax") == 0 and value(fin, "wallYPlusMax") == 0,
                    f"{ARMS[a]} {name}: stayed laminar (Re {plug * H / nu:.1f} on the"
                    " gap), which is what makes the drops linear in the rate")

    nPass = sum(int(value(c["step"][n], "passes")) for c in arms for n in c["step"])
    print(f"cost: {sum(c['solves'] for c in arms)} block solves and {nPass} interface"
          f" passes for {N_STEPS * 2} agreed time levels on two edges")

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: two detailed regions take their inertia from two different places, one"
          " edge at a time, inside one schedule that never notices the difference")


if __name__ == "__main__":
    sys.exit(main())
