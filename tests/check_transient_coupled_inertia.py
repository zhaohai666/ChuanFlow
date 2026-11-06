#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient twice on one case and check where the
detailed region's inertia comes from.

The network is the transient coupled one - Src -> A -> { X (3D block) || By (1D
bypass) } -> B -> Snk, at rest until a rate is switched on - and the two arms
differ in a single source: `inertias ( network partner )` lets arm 0 keep the
coupled pipe's own rho*L/A while arm 1 installs the rho*V/A^2 the block measures
on the region's mesh. The region is a square duct and the pipe standing for it is
circular, so the two numbers differ by 27% and everything else is identical.

Seven families of statements:

  1. where the number came from - each arm carries one inertance at every time
     level, arm 0's is the bore's rho*L/A and arm 1's the region's rho*V/A^2,
     and the block reports the measurement it handed over
  2. the row that number sits in - each arm's momentum rows close against its
     OWN M, so the installed value is the one the network solved with
  3. the interface - every time level agrees, in both arms, inside the budget
  4. the first step - a network at rest divides a sudden rate by INERTANCE, so
     the arms part company immediately, each at its own split
  5. the ramp - tau = (M_X + M_By)/(R_X + R_By) per arm, and the ratio of the
     two taus is the bore's shape error diluted by the bypass's inertia
  6. the partition - the arms differ during the ramp and reconverge: the
     operating point holds no inertia, so the same fixed point is owed either way
  7. the block - continuity inside the region and the developed-duct closures,
     unchanged by who states the inertia

usage: check_transient_coupled_inertia.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# geometry and fluid of tests/cases/transientCoupledInertia/case.dict
MU, RHO = 0.1, 800.0
H, W, L = 0.02, 0.02, 0.5      # the detailed duct: height, width, length
AREA = H * W                   # m2, the section every cross-stream plane shows
VOL = L * H * W                # m3, the region's own volume
Q_IN = 2.4e-4                  # m3/s switched on at Src at t = 0
P_REF = 1.0e6                  # Pa held at the sink
DT, N_STEPS = 0.02, 60         # time { dt 0.02; endTime 1.2; }
MAX_PASSES = 8                 # couple{ maxSubIterations }
TOL_INTERFACE = 1e-4           # couple{ tolerance } the case asks for
SEED_R = 3.75e6                # couple{ seeds }, the same in both arms

# the two candidate inertances, rho*L/A of the pipe and rho*V/A^2 of the region
SEC_X = 0.25 * math.pi * 0.02 ** 2     # pipe X as the case declares it
SEC_BY = 0.25 * math.pi * 0.03 ** 2
SEC_LEAD = 0.25 * math.pi * 0.1 ** 2
M_BORE = RHO * L / SEC_X               # what `inertias ( network )` keeps
M_REGION = RHO * VOL / (AREA * AREA)   # what `inertias ( partner )` installs
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

    def check(cond, msg):
        if not cond:
            require(False, msg)

    # ---- split the log into the two arms ----
    arms = []
    cur = None
    for ln in lines:
        if ln.startswith("transient seed "):
            cur = {"tag": ln, "step": {}, "nodes": {}, "pipes": {}, "block": [],
                   "summary": []}
            arms.append(cur)
        elif cur is None:
            continue
        elif ln.startswith("transient step "):
            cur["step"][int(ln.split()[2].rstrip(":"))] = ln
        elif ln.startswith("nodes step "):
            cur["nodes"][int(ln.split()[2].rstrip(":"))] = ln
        elif ln.startswith("pipes step "):
            cur["pipes"][int(ln.split()[2].rstrip(":"))] = ln
        elif ln.startswith("transient block "):
            cur["block"].append(ln)
        elif ln.startswith("transient coupled seed "):
            cur["summary"].append(ln)

    require(len(arms) == 2, "the case ran both inertia conventions as arms of one network")
    require("inertia network" in arms[0]["tag"]
            and "inertia partner" in arms[1]["tag"],
            "arm 0 keeps the bore's inertance, arm 1 takes the region's")
    header = [ln for ln in lines if ln.startswith("transient network+block")]
    require(len(header) == 1 and "region inertia network (arms)" in header[0],
            "one block, one convention list, both arms announced")

    series = []
    for a, arm in enumerate(arms):
        require(sorted(arm["step"]) == list(range(1, N_STEPS + 1)),
                f"arm {a}: all {N_STEPS} steps reported a state"
                f" (got {len(arm['step'])})")
        require(sorted(arm["nodes"]) == sorted(arm["step"])
                and sorted(arm["pipes"]) == sorted(arm["step"]),
                f"arm {a}: every step reports its whole 1D state")
        require(len(arm["block"]) == 1 and len(arm["summary"]) == 1,
                f"arm {a}: the block and the interface each reported once")
        series.append(arm)

    # ---- (1) where each arm's inertance came from ----
    ms, qX, qBy, qA, qB, Rsettled = [], [], [], [], [], []
    for a, arm in enumerate(series):
        each = [value(arm["step"][n], "M") for n in sorted(arm["step"])]
        check(len(set(each)) == 1,
              f"arm {a}: one inertance at every time level, not one per step")
        m = each[0]
        ms.append(m)
        want = M_BORE if a == 0 else M_REGION
        check(abs(m / want - 1.0) <= 1e-12,
              f"arm {a}: the row carries "
              + ("the coupled pipe's own rho*L/A" if a == 0
                 else "the region's rho*V/A^2")
              + f" = {want:.6g} kg/m4 (logged {m:.6g}, {abs(m / want - 1):.1e})")
        fin = arm["block"][0]
        check(abs(value(fin, "inertance") - M_REGION) <= 1e-9 * M_REGION,
              f"arm {a}: the block measures rho*V/A^2 = {M_REGION:.6g} kg/m4 on its"
              f" own mesh ({value(fin, 'inertance'):.10g})")
        # the same measurement in both arms: it is a property of the mesh, and
        # reading it must not depend on whether anyone asked for it
        check(abs(value(fin, "inertance") / M_REGION - 1.0)
              <= abs(value(series[1 - a]["block"][0], "inertance") / M_REGION - 1.0) * 2
              + 1e-15,
              f"arm {a}: the block's measurement is the same in either convention")

    require(abs(ms[0] / ms[1] - AREA / SEC_X) <= 1e-12,
            "the two arms differ by exactly the bore's shape error: rho*L/A of the"
            f" circle over rho*V/A^2 of the square is {ms[0] / ms[1]:.6f} against the"
            f" inverse area ratio {AREA / SEC_X:.6f}")
    require(abs(M_BORE / M_REGION - 1.2732395447351628) <= 1e-12,
            f"that error is 27.3%: a circle inscribed in the region's square section"
            f" holds {(1 - SEC_X / AREA) * 100:.1f}% less fluid for the same width")
    print(f"inertances: bore {ms[0]:.6g}, region {ms[1]:.6g}, bypass {M_BY:.6g} kg/m4")

    # ---- (2) the rows the number sits in, closed per arm ----
    worst = {}
    for a, arm in enumerate(series):
        w = {"X": 0.0, "By": 0.0, "A": 0.0, "B": 0.0}
        xq, bq, aq, bqq = [], [], [], []
        for n in sorted(arm["step"]):
            s, nd, pp = arm["step"][n], arm["nodes"][n], arm["pipes"][n]
            q = value(s, "q")
            x, by, la, lb = (value(pp, k) for k in ("X", "By", "A", "B"))
            pSrc, pJ1, pJ2, pSnk = (value(nd, k) for k in ("Src", "J1", "J2", "Snk"))
            R, dpNet, res = (value(s, k) for k in ("R", "dpNetwork", "residual"))
            xPrev, byPrev = (xq[-1] if xq else 0.0, bq[-1] if bq else 0.0)
            aPrev, bPrev = (aq[-1] if aq else 0.0, bqq[-1] if bqq else 0.0)

            check(abs(value(s, "t") - n * DT) <= 1e-12,
                  f"arm {a} step {n}: t is {n} steps of {DT} s")
            check(abs(q - x) <= 1e-12 * q,
                  f"arm {a} step {n}: the reported rate is the coupled pipe's")

            # (3) the interface at this time level
            check(abs(res - abs(dpNet - value(s, "dpBlock")) / value(s, "dpBlock"))
                  <= 1e-9,
                  f"arm {a} step {n}: residual is the cross-solver gap")
            check(res <= TOL_INTERFACE,
                  f"arm {a} step {n}: interface residual {res:.3e} under the tolerance")
            check(abs(dpNet - R * q) <= 1e-12 * abs(dpNet),
                  f"arm {a} step {n}: the network charges its coefficient times the rate")
            check(value(s, "blockConverged") == 1,
                  f"arm {a} step {n}: the block converged at the reported state")
            check(value(s, "newton") < 1e-11,
                  f"arm {a} step {n}: the 1D step solved its own system to its floor")
            check(1 <= value(s, "passes") <= MAX_PASSES,
                  f"arm {a} step {n}: agreed inside the sub-iteration budget")

            # (2) continuity and the momentum rows, with THIS arm's M
            check(abs(la - Q_IN) <= 1e-12 * Q_IN,
                  f"arm {a} step {n}: the lead carries the whole injected rate")
            check(abs(la - x - by) <= 1e-12 * Q_IN,
                  f"arm {a} step {n}: upstream node: in = coupled + bypass")
            check(abs(lb - x - by) <= 1e-12 * Q_IN,
                  f"arm {a} step {n}: downstream node: out = coupled + bypass")
            check(abs(pSnk - P_REF) <= 1e-9,
                  f"arm {a} step {n}: the sink holds its reference pressure")
            dp12 = pJ1 - pJ2
            for key, drop, flow, prev, inert, resist in (
                    ("X", dp12, x, xPrev, ms[a], R),
                    ("By", dp12, by, byPrev, M_BY, hagenPoiseuille(1.5, 0.03)),
                    ("A", pSrc - pJ1, la, aPrev, M_LEAD, hagenPoiseuille(5.0, 0.1)),
                    ("B", pJ2 - pSnk, lb, bPrev, M_LEAD, hagenPoiseuille(5.0, 0.1))):
                dev = abs(drop - (resist * flow + inert * (flow - prev) / DT)) / abs(drop)
                w[key] = max(w[key], dev)
            xq.append(x); bq.append(by); aq.append(la); bqq.append(lb)
        qX.append(xq); qBy.append(bq); qA.append(aq); qB.append(bqq)
        worst[a] = w
        Rsettled.append(value(arm["step"][N_STEPS], "R"))

    for a, w in worst.items():
        require(w["X"] <= 1e-8,
                f"arm {a}: the coupled region's node drop is R*Q + M*dQ/dt with the"
                f" inertance that arm installed (worst {w['X']:.2e})")
        require(w["By"] <= 1e-9 and w["A"] <= 1e-9 and w["B"] <= 1e-9,
                f"arm {a}: the bypass and the leads close the same way, inertial kick"
                f" on the first step included (worst {w['By']:.1e}/{w['A']:.1e}/{w['B']:.1e})")

    R_BY = hagenPoiseuille(1.5, 0.03)
    fits = []
    for a, arm in enumerate(series):
        q = qX[a]
        # ---- (4) the first step splits the sudden rate by inertance ----
        q1 = Q_IN * (R_BY + M_BY / DT) / ((ms[a] + M_BY) / DT
                                           + value(arm["step"][1], "R") + R_BY)
        dev1 = abs(q[0] - q1) / q1
        check(dev1 <= 1e-3,
              f"arm {a}: the first step divides the rate by inertance"
              f" ({q[0]:.6g} vs {q1:.6g} m3/s, {dev1:.2e})")
        check(q[0] < 0.99 * q[-1],
              f"arm {a}: the series really ramps, starting {q[0] / q[-1]:.3f} of the way"
              " to its operating point")

        # ---- (5) the ramp: q^n = A + B q^(n-1) ----
        xs, ys = q[1:-1], q[2:]
        k = len(xs)
        sx, sy = sum(xs), sum(ys)
        B = (k * sum(p * f for p, f in zip(xs, ys)) - sx * sy) / (
            k * sum(v * v for v in xs) - sx * sx)
        A = (sy - B * sx) / k
        tauFit = -DT / math.log(B)
        qInfFit = A / (1.0 - B)
        tau = (ms[a] + M_BY) / (Rsettled[a] + R_BY)
        qInf = Q_IN * R_BY / (Rsettled[a] + R_BY)
        print(f"arm {a} ramp: tau {tauFit:.6g} s (closed form {tau:.6g}),"
              f" q_inf {qInfFit:.6g} m3/s, decay/step {B:.6g}")
        require(0.9 < B < 1.0,
                f"arm {a}: one time constant spans {1.0 / (1.0 - B):.1f} steps")
        require(abs(tauFit - tau) <= 0.05 * tau,
                f"arm {a}: tau is (M_X + M_By)/(R_X + R_By) with M_X read from "
                + ("the bore" if a == 0 else "the region") + ":"
                f" {tau:.6g} s against the fitted {tauFit:.6g} s"
                f" ({(tauFit / tau - 1) * 100:+.2f}%)")
        fits.append({"B": B, "tau": tauFit, "tauPred": tau, "qInf": qInfFit})

    # ---- (5b) the ratio of the two taus is the bore's shape error ----
    ratioFit = fits[0]["tau"] / fits[1]["tau"]
    ratioPred = (ms[0] + M_BY) / (ms[1] + M_BY)
    require(abs(ratioFit / ratioPred - 1.0) <= 0.02,
            f"reading the inertia off the bore instead of the region stretches the"
            f" ramp by {ratioPred:.5f}: measured {ratioFit:.5f}"
            f" ({(ratioFit / ratioPred - 1) * 100:+.2f}%)")
    require(fits[0]["tau"] > fits[1]["tau"],
            f"the arm with the larger inertance ramps slower"
            f" ({fits[0]['tau']:.4f} s against {fits[1]['tau']:.4f} s)")

    # ---- (6) the arms part company on the way and meet at the operating point ----
    d = [abs(x - y) for x, y in zip(qX[0], qX[1])]
    iMax = d.index(max(d))
    qInf = 0.5 * (fits[0]["qInf"] + fits[1]["qInf"])
    require(iMax == 0,
            f"the two conventions differ most on the very first step, when the whole"
            f" split is inertia (step {iMax + 1})")
    require(qX[0][0] < qX[1][0],
            f"and the lag is the inertia's: the arm that overstates it takes"
            f" {qX[0][0] / qX[1][0]:.4f} of the rate the other one delivers at the same"
            " instant")
    require(max(d) >= 0.05 * qInf,
            f"and the bore's 27% inertia error moves the rate by {max(d) / qInf * 100:.1f}%"
            f" of the operating point ({max(d):.3e} m3/s)")
    require(d[-1] <= 0.02 * max(d),
            f"an inertia the loop has settled out decays away: the gap narrows to"
            f" {d[-1]:.3e} m3/s, {d[-1] / max(d) * 100:.2f}% of its peak")
    require(abs(qX[0][-1] / qX[1][-1] - 1.0) <= 5e-3,
            f"both arms end within {abs(qX[0][-1] / qX[1][-1] - 1) * 100:.3f}% of the same"
            f" rate ({qX[0][-1]:.6g} and {qX[1][-1]:.6g} m3/s)")
    require(abs(fits[0]["qInf"] / fits[1]["qInf"] - 1.0) <= 1e-3,
            f"and extrapolate to the same operating point to within"
            f" {abs(fits[0]['qInf'] / fits[1]['qInf'] - 1) * 1e2:.4f}% - inertia sets the"
            " path, resistance sets the destination")

    # ---- (7) the detailed region, in each arm ----
    for a, arm in enumerate(series):
        fin = arm["block"][0]
        q = value(fin, "fluxIn")
        plug = value(fin, "plugSpeed")
        dpPa = value(fin, "dpPa")
        require(abs(q / qX[a][-1] - 1.0) <= 1e-6,
                f"arm {a}: the block carries the rate its arm ended with"
                f" ({q:.6g} m3/s)")
        require(abs(value(fin, "inletArea") * plug - q) <= 1e-12 * q,
                f"arm {a}: plug speed times inlet area is the carried rate")
        require(abs(value(fin, "fluxOut") - q) <= 1e-7 * q,
                f"arm {a}: what enters the detailed region leaves it")
        require(value(fin, "planeFluxMax") - value(fin, "planeFluxMin") <= 1e-7 * q,
                f"arm {a}: every cross-section plane inside carries the same rate")
        require(abs(value(fin, "inletArea") - AREA) <= 1e-12 * AREA,
                f"arm {a}: the inlet patch is the section the inertance was referred to"
                f" ({value(fin, 'inletArea'):.6g} m2)")
        gradExact = 12.0 * MU * plug / (H * H)
        require(abs(abs(value(fin, "probeSlope")) - gradExact) <= 0.03 * gradExact,
                f"arm {a}: developed dp/dx {abs(value(fin, 'probeSlope')):.4g} vs planar"
                f" Poiseuille {gradExact:.4g}")
        require(abs(value(fin, "uxCentre") / plug - 1.5) <= 0.03,
                f"arm {a}: centreline speed is 3/2 the mean")
        require(abs(value(arm["step"][N_STEPS], "dpBlock") - dpPa) <= 1e-12,
                f"arm {a}: the drop the interface ended with is the drop the block's"
                " own field measures")

    print(f"interface: {N_STEPS * 2} time levels over two conventions, worst residual"
          f" {max(value(arm['step'][n], 'residual') for arm in series for n in arm['step']):.3e}"
          f" under {TOL_INTERFACE:g}")
    total = sum(int(value(arm["step"][n], "sweeps")) for arm in series for n in arm["step"])
    passes = sum(int(value(arm["step"][n], "passes")) for arm in series for n in arm["step"])
    print(f"cost: {total} block sweeps and {passes} interface passes for"
          f" {N_STEPS * 2} agreed time levels")

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: the region's inertia is one number with two sources, and the network"
          " solves with whichever one the case names")


if __name__ == "__main__":
    sys.exit(main())
