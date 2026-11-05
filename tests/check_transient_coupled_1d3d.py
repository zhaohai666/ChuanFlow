#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient and verify the partitioned 1D<->3D interface
at every time level.

The case is the steady coupled network put in motion: Src -> A -> { X (3D block)
|| By (1D bypass) } -> B -> Snk, at rest until an injected rate is switched on at
t = 0. Seven families of statements are checked:

  1. one state per step     - the block carries the rate the step delivers, and
                              the residual is the gap between two solvers rather
                              than a number compared with itself
  2. the discrete 1D law    - every pipe's node drop is its partner resistance (or
                              its own friction) times its NEW rate plus its
                              inertance times the rate change: the momentum row
                              the network solved, closed with the block's answer
  3. continuity             - the leads carry the whole rate at every time level
                              and the branch split adds up
  4. the interface schedule - the starting coefficient is corrected on the first
                              step and every later step agrees within a couple of
                              partner runs, at the tolerance the case asks for
  5. the ramp               - the sudden rate first splits by INERTANCE and then
                              relaxes with tau = (M_X + M_By)/(R_X + R_By), so a
                              system carrying the region's inertance twice reads
                              there (the region holds 43% of the loop's inertance)
  6. the settled end        - the series lands on the operating point a laminar
                              duct defines, and the block's field at that rate is
                              developed Poiseuille flow
  7. block continuity       - what enters the detailed region leaves it, and every
                              cross-section plane inside carries the same rate

usage: check_transient_coupled_1d3d.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# geometry and fluid of tests/cases/transientCoupled1d3d/case.dict
MU, RHO = 0.1, 800.0
H, W, L = 0.02, 0.02, 0.5      # duct height (wall-normal), width, length
AREA = H * W
Q_IN = 2.4e-4                  # m3/s switched on at Src at t = 0
P_REF = 1.0e6                  # Pa held at the sink
DT, N_STEPS = 0.02, 60         # time { dt 0.02; endTime 1.2; }
MAX_PASSES = 8                 # couple{ maxSubIterations }
TOL_INTERFACE = 1e-4           # couple{ tolerance } the case asks for

# the 1D side of the loop: the coupled region keeps its inertance here, so these
# are the numbers the discrete closed form is built from
SEC_X = 0.25 * math.pi * 0.02 ** 2     # pipe X: length 0.5, diameter 0.02
SEC_BY = 0.25 * math.pi * 0.03 ** 2    # pipe By: length 1.5, diameter 0.03
SEC_LEAD = 0.25 * math.pi * 0.1 ** 2   # pipes A and B: length 5, diameter 0.1
M_X = RHO * 0.5 / SEC_X
M_BY = RHO * 1.5 / SEC_BY
M_LEAD = RHO * 5.0 / SEC_LEAD


def value(line, key):
    """First number stored under `key` in a log line ('q 1.5e-4 ...')."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([-+0-9.eE]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return float(m.group(1))


def hagenPoiseuille(length, diameter):
    return 128.0 * MU * length / (math.pi * diameter ** 4)


def read_vtu_cells(path):
    """Cell centroids (x, y) and the streamwise/transverse velocity components."""
    txt = open(path).read()
    pts = re.search(
        r'<DataArray type="Float64" NumberOfComponents="3" format="ascii">\s*(.*?)</DataArray>',
        txt, re.S).group(1)
    coords = [float(v) for v in pts.split()]
    points = list(zip(coords[0::3], coords[1::3], coords[2::3]))

    conn = re.search(r'Name="connectivity" format="ascii">\s*(.*?)</DataArray>', txt, re.S).group(1)
    offs = re.search(r'Name="offsets" format="ascii">\s*(.*?)</DataArray>', txt, re.S).group(1)
    uvals = re.search(
        r'<DataArray type="Float64" Name="U" NumberOfComponents="3" format="ascii">\s*(.*?)</DataArray>',
        txt, re.S).group(1)
    uv = [float(v) for v in uvals.split()]

    ci = [int(v) for v in conn.split()]
    oi = [int(v) for v in offs.split()]
    cells, prev, k = [], 0, 0
    for o in oi:
        vids = ci[prev:o]
        prev = o
        n = len(vids)
        cx = sum(points[v][0] for v in vids) / n
        cy = sum(points[v][1] for v in vids) / n
        cells.append((cx, cy, uv[3 * k], uv[3 * k + 1]))
        k += 1
    return cells


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("transient coupled solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    # the same statement, for the assertions that hold at every one of the 60
    # time levels: a passing step has nothing to say
    def check(cond, msg):
        if not cond:
            require(False, msg)

    # ---- collect the series ----
    step, nodes, pipes = {}, {}, {}
    for ln in lines:
        if ln.startswith("transient step "):
            step[int(ln.split()[2].rstrip(":"))] = ln
        elif ln.startswith("nodes step "):
            nodes[int(ln.split()[2].rstrip(":"))] = ln
        elif ln.startswith("pipes step "):
            pipes[int(ln.split()[2].rstrip(":"))] = ln

    require(sorted(step) == list(range(1, N_STEPS + 1)),
            f"every one of the {N_STEPS} steps the case asks for reported a state"
            f" (got {len(step)})")
    require(sorted(nodes) == sorted(step) and sorted(pipes) == sorted(step),
            "every step reports its whole 1D state")

    blk = [ln for ln in lines if ln.startswith("transient block ")]
    require(len(blk) == 1, "the block reports the field it ended on")
    summary = [ln for ln in lines if ln.startswith("transient coupled seed ")]
    require(len(summary) == 1, "one interface series ran")

    # ---- per step: one state, the discrete 1D law, continuity ----
    qX, qBy, qA, qB = [], [], [], []
    worst = {"X": 0.0, "By": 0.0, "A": 0.0, "B": 0.0}
    for n in sorted(step):
        s, nd, pp = step[n], nodes[n], pipes[n]
        t, q = value(s, "t"), value(s, "q")
        R, dpNet, dpBlock = value(s, "R"), value(s, "dpNetwork"), value(s, "dpBlock")
        res = value(s, "residual")
        x, by, a, b = value(pp, "X"), value(pp, "By"), value(pp, "A"), value(pp, "B")
        pSrc, pJ1, pJ2, pSnk = (value(nd, k) for k in ("Src", "J1", "J2", "Snk"))
        xPrev = qX[-1] if qX else 0.0
        byPrev = qBy[-1] if qBy else 0.0
        aPrev = qA[-1] if qA else 0.0
        bPrev = qB[-1] if qB else 0.0

        # the clock, and one rate: what the interface reports is what the
        # network's own state says
        check(abs(t - n * DT) <= 1e-12, f"step {n}: t is {n} steps of {DT} s")
        check(abs(q - x) <= 1e-12 * q, f"step {n}: the reported rate is the coupled pipe's")

        # (1) the residual is the gap between the two solvers, and it is what it says
        check(abs(res - abs(dpNet - dpBlock) / dpBlock) <= 1e-9,
              f"step {n}: residual is the cross-solver gap")
        check(res <= TOL_INTERFACE,
              f"step {n}: interface residual {res:.3e} below the case tolerance")
        check(abs(dpNet - R * q) <= 1e-12 * abs(dpNet),
              f"step {n}: the network charges its coefficient times the rate")
        check(value(s, "blockConverged") == 1,
              f"step {n}: the block solve converged at the reported state")
        check(value(s, "newton") < 1e-11,
              f"step {n}: the 1D step solved its own system to its floor")

        # (3) node continuity at this time level
        check(abs(a - Q_IN) <= 1e-12 * Q_IN,
              f"step {n}: the lead carries the whole injected rate")
        check(abs(a - x - by) <= 1e-12 * Q_IN,
              f"step {n}: upstream node: in = coupled + bypass")
        check(abs(b - x - by) <= 1e-12 * Q_IN,
              f"step {n}: downstream node: out = coupled + bypass")

        # (2) the momentum rows the step solved, closed with the block's answer:
        #     the node drop is the loss at the NEW rate plus inertance times the
        #     rate change. The coupled pipe's loss is the interface coefficient,
        #     so this is the statement the two solvers shake hands on.
        dp12 = pJ1 - pJ2
        for key, drop, flow, prev, inert, resist in (
                ("X", dp12, x, xPrev, M_X, R),
                ("By", dp12, by, byPrev, M_BY, hagenPoiseuille(1.5, 0.03)),
                ("A", pSrc - pJ1, a, aPrev, M_LEAD, hagenPoiseuille(5.0, 0.1)),
                ("B", pJ2 - pSnk, b, bPrev, M_LEAD, hagenPoiseuille(5.0, 0.1))):
            dev = abs(drop - (resist * flow + inert * (flow - prev) / DT)) / abs(drop)
            worst[key] = max(worst[key], dev)
        check(abs(pSnk - P_REF) <= 1e-9,
              f"step {n}: the sink holds its reference pressure")
        qX.append(x); qBy.append(by); qA.append(a); qB.append(b)

    require(worst["X"] <= 1e-8,
            f"the coupled region's node drop is R*Q + M*dQ/dt at every step"
            f" (worst {worst['X']:.2e})")
    require(worst["By"] <= 1e-9,
            f"the bypass's laminar friction plus its inertance closes at every step"
            f" (worst {worst['By']:.2e})")
    require(worst["A"] <= 1e-9 and worst["B"] <= 1e-9,
            f"the leads close too, inertial kick on the first step included"
            f" (worst {worst['A']:.2e}, {worst['B']:.2e})")

    # ---- (4) the interface schedule ----
    passes = [int(value(step[n], "passes")) for n in sorted(step)]
    require(all(1 <= p <= MAX_PASSES for p in passes),
            "every step agrees inside the allowed sub-iterations")
    require(passes[0] == 2,
            "the starting coefficient is corrected on the first step")
    require(sum(passes) <= 2 * N_STEPS,
            f"the series costs {sum(passes)} partner runs for {N_STEPS} steps"
            " - a laminar region's coefficient is rate-independent")

    # the coefficient the interface converged to, and the operating point it
    # implies: the closed form below is built from these, not from inputs
    R_X = value(step[N_STEPS], "R")
    R_BY = hagenPoiseuille(1.5, 0.03)
    qInf = R_BY * Q_IN / (R_X + R_BY)
    tau = (M_X + M_BY) / (R_X + R_BY)

    # ---- (5) the ramp ----
    # the first step: the network was at rest, so continuity did not hold at the
    # previous time level and the sudden rate divides by INERTANCE
    q1Inertial = Q_IN * (R_BY + M_BY / DT) / ((M_X + M_BY) / DT + R_X + R_BY)
    dev1 = abs(qX[0] - q1Inertial) / q1Inertial
    require(dev1 <= 1e-3,
            f"the first step splits the sudden rate by inertance"
            f" ({qX[0]:.6g} vs {q1Inertial:.6g} m3/s, {dev1:.2e})")
    require(qX[0] < 0.95 * qInf,
            f"the series really ramps: it starts {qX[0] / qInf:.3f} of the way to"
            " its operating point, not on it")
    require(all(qX[i] < qX[i + 1] for i in range(N_STEPS - 1))
            and all(qBy[i] > qBy[i + 1] for i in range(N_STEPS - 1)),
            "the rate moves toward the resistive split monotonically, one branch"
            " up and the other down")

    # every later step satisfies q^n = A + B q^(n-1); recovering A and B from the
    # series gives the time constant and the asymptote the two solvers produced
    xs, ys = qX[1:-1], qX[2:]
    k = len(xs)
    sx, sy = sum(xs), sum(ys)
    B = (k * sum(p * f for p, f in zip(xs, ys)) - sx * sy) / (k * sum(v * v for v in xs) - sx * sx)
    A = (sy - B * sx) / k
    tauFit = -DT / math.log(B)
    qInfFit = A / (1.0 - B)
    print(f"measured ramp: tau {tauFit:.6g} s (closed form {tau:.6g}),"
          f" q_inf {qInfFit:.6g} m3/s (closed form {qInf:.6g}), decay/step {B:.6g}")
    require(0.9 < B < 1.0, f"one time constant spans {1.0 / (1.0 - B):.1f} steps")
    require(abs(qInfFit - qInf) <= 0.01 * qInf,
            "the series relaxes to the resistive split of the two branches")
    require(abs(tauFit - tau) <= 0.05 * tau,
            f"the ramp's time constant is (M_X + M_By)/(R_X + R_By): counting the"
            f" region's inertance twice would read {tau * (1 + M_X / (M_X + M_BY)):.4g} s")

    # ---- (6) the settled end ----
    gap = (qInf - qX[-1]) / qInf
    require(0.0 < gap <= 0.005,
            f"the series ends {gap * 100:.2f}% short of the operating point, the"
            " exponential tail its own time constant predicts")
    rBlock = 12.0 * MU * L / (H ** 3 * W)
    print(f"settled state: block R {R_X:.6g} Pa.s/m3 (planar duct {rBlock:.6g}),"
          f" q_X {qX[-1]:.6g} m3/s, dp {value(step[N_STEPS], 'dpBlock'):.6g} Pa")
    require(1.0 <= R_X / rBlock <= 1.06,
            "the region's settled resistance is the laminar duct value plus its"
            " entrance length")
    reX = RHO * (qX[-1] / SEC_X) * 0.02 / MU     # the 1D kernel's own convention
    require(20.0 <= reX <= 2000.0,
            f"the detailed region stayed laminar at Re {reX:.1f}, as the closed"
            " forms above assume")
    require(abs(value(step[N_STEPS], "dpBlock") - value(blk[0], "dpPa")) <= 1e-12,
            "the drop the interface ended with is the drop the block's own field measures")

    # ---- (7) the block at the settled rate: continuity inside the region ----
    fin = blk[0]
    q, plug = value(fin, "fluxIn"), value(fin, "plugSpeed")
    dpBlock = value(fin, "dpPa")
    slope = value(fin, "probeSlope")
    require(abs(q - qX[-1]) <= 1e-10 * q, "the block carries the settled branch rate")
    require(abs(value(fin, "inletArea") * plug - q) <= 1e-12 * q,
            "plug speed times inlet area is the carried rate")
    require(abs(value(fin, "fluxOut") - q) <= 1e-7 * q,
            "what enters the detailed region leaves it")
    require(value(fin, "planeFluxMax") - value(fin, "planeFluxMin") <= 1e-7 * q,
            "every cross-section plane inside carries the same rate")

    gradExact = 12.0 * MU * plug / (H * H)
    require(abs(abs(slope) - gradExact) <= 0.03 * gradExact,
            f"developed dp/dx {abs(slope):.4g} vs planar Poiseuille {gradExact:.4g}")
    require(value(fin, "probeSections") >= 20, "enough sections inside the probe window")
    require(value(fin, "probeSpread") <= 0.002 * abs(dpBlock),
            f"section means hug the developed line (spread {value(fin, 'probeSpread'):.3g} Pa)")
    require(abs(value(fin, "uxCentre") / plug - 1.5) <= 0.03,
            "centreline speed is 3/2 the mean")

    cells = read_vtu_cells(os.path.join(workdir, "out_transient_coupled_block.vtu"))
    dev = [c for c in cells if 0.15 <= c[0] <= 0.45]
    if len(dev) < 20 * 20:
        raise SystemExit(f"expected developed columns in the .vtu, found {len(dev)} cells")
    peak = 1.5 * plug

    def profile(y):
        return peak * (1.0 - ((y - H / 2.0) / (H / 2.0)) ** 2)

    profErr = max(abs(c[2] - profile(c[1])) for c in dev)
    vMax = max(abs(c[3]) for c in dev)
    # a column of cells spans dy = H/20 and the full duct width
    columnArea = (H / 20.0) * W
    fluxCol = min(sum(c[2] for c in dev if abs(c[0] - x) < 1e-9) * columnArea
                  for x in sorted({round(c[0], 6) for c in dev}))
    require(profErr <= 0.05 * peak, f"velocity profile is parabolic (err {profErr:.4g})")
    require(vMax <= 0.02 * plug, f"no cross-stream velocity ({vMax:.3g} m/s)")
    require(abs(fluxCol - q) <= 0.02 * q,
            "the written field carries the branch rate column by column")

    # ---- the summary line describes the same state as the series ----
    su = summary[0]
    require(abs(value(su, "q") - qX[-1]) <= 1e-12 * q and
            abs(value(su, "R") - R_X) <= 1e-12 * R_X and
            int(value(su, "steps")) == N_STEPS,
            "the summary reports the last step of the series")

    if fails:
        raise SystemExit(f"{len(fails)} transient coupled-interface checks FAILED")
    print("transient coupled 1D<->3D verification passed")


if __name__ == "__main__":
    main()
