#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-multi and verify a network with TWO detailed regions.

The case is Src -> A -> { X1 (3D duct) || By (1D bypass) } -> X2 (3D duct) -> B ->
Snk. X1 sits in a parallel branch, so the node split decides its rate; X2 sits in
series, so topology pins its rate to the injected one. That pairing lets the
checker read some statements off the network alone and others off each block's own
field. Eight families are checked:

  1. one state per edge     - each block carries the rate its edge delivers, the
                              drop the network charged is R times that rate, and
                              the residual is the gap between the two solvers
  2. the 1D network         - leads carry the whole rate, the branch split adds up,
                              the two parallel branches share their drop, and the
                              node pressures chain across both interfaces
  3. the interface schedule - every start begins far off, settles inside the case
                              tolerance and inside the outer budget, and one pass
                              costs exactly one network solve plus one solve per
                              block
  4. edges do not lean on each other - X2 starts and is fed identically in all
                              three runs, so its first answer must be the same
                              digits, while X1's rate wanders by 2x between them
  5. the settled state      - all three starting sets land on the same per-edge
                              rate, drop and resistance
  6. the operating point    - the series-parallel closed form built from each
                              duct's own planar Poiseuille resistance and the
                              bypass's Hagen-Poiseuille one
  7. each block's field     - conservation, one rate per plane, developed gradient
                              and the parabolic profile of that block's own duct
  8. laminar              - both ducts stayed laminar, which is what makes the
                              drops linear in the rate and family 6 valid

usage: check_coupled_1d3d_multi.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# geometry and fluid of tests/cases/coupled1d3dMulti/case.dict
MU, RHO = 0.1, 800.0
Q_IN = 3.0e-5                  # m3/s injected at Src
P_REF = 1.0e6                  # Pa held at Snk
TOL_INTERFACE = 1e-4           # couple{ tolerance }
MAX_OUTER = 8                  # couple{ maxOuterIterations }

EDGES = {
    # edge: (length, gap, width) of the duct that replaces the 1D law there
    "X1": (0.5, 0.02, 0.02),
    "X2": (0.4, 0.01, 0.03),
}
# the 1D side: the bypass shares a node pair with X1, the leads carry everything
BY = (1.5, 0.03)
LEAD = (5.0, 0.1)
# the coupled pipe's own diameter, used only for the reported Reynolds number
PIPE_D = {"X1": 0.02, "X2": 0.015}


def value(line, key):
    """First number stored under `key` in a log line ('q 1.5e-4 ...')."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([-+0-9.eE]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return float(m.group(1))


def run(binary, case):
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("coupled solver failed")
    return r.stdout, workdir


def planarDuct(length, gap, width):
    """Resistance of a plane channel: q = gap^3 width dp / (12 mu L)."""
    return 12.0 * MU * length / (gap ** 3 * width)


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
    out, workdir = run(binary, case)
    lines = [ln.replace("[chflow] ", "") for ln in out.splitlines()]

    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    # ---- collect the per-edge records ----
    starts = {}                      # start name -> its records
    for ln in lines:
        if ln.startswith("interface "):
            t = ln.split()
            starts.setdefault(t[1], {}).setdefault("passes", []).append(
                {"iter": int(value(ln, "iter")), "edge": t[5], "R": value(ln, "R"),
                 "q": value(ln, "q"), "dpBlock": value(ln, "dpBlock"),
                 "sweeps": value(ln, "sweeps"),
                 "blockConverged": value(ln, "blockConverged"),
                 "residual": value(ln, "residual")})
        elif ln.startswith("settled "):
            t = ln.split()
            starts.setdefault(t[1], {}).setdefault("settled", {})[t[3]] = ln
        elif ln.startswith("coupled ") and ": iterations" in ln:
            t = ln.split()
            starts.setdefault(t[1].rstrip(":"), {})["final"] = ln
        elif ln.startswith("coupledEdges "):
            t = ln.split()
            name = t[1]
            pairs = {}
            for i in range(4, len(t), 3):
                pairs[t[i]] = {"dp": float(t[i + 1]), "q": float(t[i + 2])}
            starts.setdefault(name, {})["network"] = pairs
        elif ln.startswith("nodes "):
            t = ln.split()
            starts.setdefault(t[1].rstrip(":"), {})["nodes"] = {
                t[i]: float(t[i + 1]) for i in range(2, len(t), 2)}
        elif ln.startswith("pipes "):
            t = ln.split()
            starts.setdefault(t[1].rstrip(":"), {})["pipes"] = {
                t[i]: float(t[i + 1]) for i in range(2, len(t), 2)}
        elif ln.startswith("block ") and " edge " in ln and ": cells" in ln:
            t = ln.split()
            starts.setdefault(t[1], {}).setdefault("blocks", {})[t[3].rstrip(":")] = ln

    require(set(starts) == {"lo", "mix", "high"},
            f"three starting sets ran (got {sorted(starts)})")
    require(all(set(s["settled"]) == set(EDGES) for s in starts.values()),
            "every run reports a settled state for both coupled edges")

    # ---- 1 & 3: per edge, the interface is one exchange of two numbers ----
    for name in sorted(starts):
        s = starts[name]
        iters = value(s["final"], "iterations")
        for p in s["passes"]:
            e = p["edge"]
            k = f"{name} pass {p['iter']} {e}"
            # the two sides of one pass, as the log quotes them: the coefficient
            # the network was charged with, the rate it delivered, and the drop
            # the partner reported at that rate
            require(abs(p["residual"] - abs(p["R"] * p["q"] - p["dpBlock"]) / p["dpBlock"])
                    <= 1e-9 * p["residual"] + 1e-15,
                    f"{k}: residual is the cross-solver gap")
            if p["iter"] == iters:
                net = s["network"][e]
                require(abs(p["q"] - net["q"]) <= 1e-12 * Q_IN,
                        f"{k}: the exchange uses the network's own rate")
                require(abs(p["R"] * p["q"] - net["dp"]) <= 1e-9 * abs(net["dp"]),
                        f"{k}: the charged drop is resistance times rate")
                require(abs(p["residual"] - abs(net["dp"] - p["dpBlock"]) / p["dpBlock"])
                        <= 1e-9, f"{k}: the accepted gap is between the two solves")
                require(abs(value(s["blocks"][e], "fluxIn") - p["q"]) <= 1e-12 * p["q"],
                        f"{k}: the block was measured at the rate being reported")
        fin = s["final"]
        require("converged" in fin and "NOT" not in fin,
                f"{name}: the multi-edge interface converged")
        require(1 <= iters <= MAX_OUTER, f"{name}: settled in {int(iters)} passes")
        last = [p for p in s["passes"] if p["iter"] == iters]
        worst = value(fin, "worstResidual")
        require(abs(worst - max(p["residual"] for p in last)) <= 1e-18,
                f"{name}: the pass stops on the worst edge, not on one of them")
        require(max(p["residual"] for p in last) <= TOL_INTERFACE,
                f"{name}: every edge is inside the case tolerance")
        for p in last:
            st = s["settled"][p["edge"]]
            require(abs(value(st, "residual") - p["residual"]) <= 1e-15 * p["residual"]
                    and value(st, "R") == p["R"] and value(st, "q") == p["q"],
                    f"{name} {p['edge']}: the settled line repeats the accepted pass")

    # one pass = one solve per block, so the SIMPLE runs are counted exactly
    blockSolves = sum(1 for ln in lines if ln.startswith("SIMPLE iter 1:"))
    passes = sum(len(s["passes"]) for s in starts.values()) // len(EDGES)
    require(blockSolves == passes * len(EDGES),
            f"{blockSolves} block solves for {passes} passes over {len(EDGES)} edges")

    # ---- 2: the 1D network on its own terms ----
    for name in sorted(starts):
        s = starts[name]
        pr, q = s["pipes"], s["network"]
        require(abs(pr["A"] - Q_IN) <= 1e-8 * Q_IN and abs(pr["B"] - Q_IN) <= 1e-8 * Q_IN,
                f"{name}: both leads carry the injected rate")
        require(abs(q["X2"]["q"] - Q_IN) <= 1e-8 * Q_IN,
                f"{name}: the series edge carries the whole rate")
        require(abs(pr["X1"] + pr["By"] - Q_IN) <= 1e-8 * Q_IN,
                f"{name}: the parallel branch split adds up")
        rBy = hagenPoiseuille(*BY)
        dpBy = rBy * pr["By"]
        require(abs(q["X1"]["dp"] - dpBy) <= 1e-6 * dpBy,
                f"{name}: the detailed branch and the bypass share their drop")
        n = s["nodes"]
        require(abs(n["Snk"] - P_REF) < 1e-9, f"{name}: the sink holds its pressure")
        for e, (a, b) in (("X1", ("J1", "J2")), ("X2", ("J2", "J3"))):
            require(abs((n[a] - n[b]) - q[e]["dp"]) <= 1e-9 * q[e]["dp"],
                    f"{name}: {e}'s drop is its own node difference")
        chain = (hagenPoiseuille(*LEAD) * Q_IN + q["X1"]["dp"] + q["X2"]["dp"]
                 + hagenPoiseuille(*LEAD) * Q_IN)
        require(abs(chain - (n["Src"] - n["Snk"])) <= 1e-9 * chain,
                f"{name}: the pressure chain closes across both interfaces")

    # ---- 4: the edges are separate problems, not one exchange twice ----
    first1 = {name: [p for p in starts[name]["passes"] if p["iter"] == 1]
              for name in starts}
    x2 = {name: [p for p in lst if p["edge"] == "X2"][0] for name, lst in first1.items()}
    require(len({x2[n]["dpBlock"] for n in x2}) == 1,
            "X2's first answer is the same digits in all three runs: it never saw "
            "what X1 was doing")
    x1q = {[p for p in first1[n] if p["edge"] == "X1"][0]["q"] for n in first1}
    require(max(x1q) > 2.0 * min(x1q),
            "X1's rate really moves between the starts, so the two edges are not "
            "the same exchange twice")
    require(all(p["R"] == 5.3e6 and abs(p["q"] - Q_IN) <= 1e-12 * Q_IN
                for p in x2.values()),
            "every run hands X2 the same starting coefficient and the same rate")

    # ---- 5: the fixed point does not belong to the starting set ----
    for e in EDGES:
        for field in ("q", "dp", "R"):
            vals = [value(starts[n]["settled"][e], field) for n in starts]
            span = (max(vals) - min(vals)) / max(vals)
            require(span <= 1e-3,
                    f"the three starts agree on {e}'s {field} to {span:.2e}")

    # ---- 6: closed form of the network plus both ducts ----
    rDuct = {e: planarDuct(*EDGES[e]) for e in EDGES}
    rBy = hagenPoiseuille(*BY)
    rLead = hagenPoiseuille(*LEAD)
    rBranch = 1.0 / (1.0 / rDuct["X1"] + 1.0 / rBy)
    dpExact = Q_IN * (2.0 * rLead + rBranch + rDuct["X2"])
    qX1Exact = Q_IN * rBy / (rBy + rDuct["X1"])
    last = starts["high"]
    print(f"analytic: dp {dpExact:.4g} Pa, q_X1 {qX1Exact:.6g} m3/s, "
          f"R_X1 {rDuct['X1']:.4g}, R_X2 {rDuct['X2']:.4g} Pa.s/m3")
    dpMeasured = last["nodes"]["Src"] - last["nodes"]["Snk"]
    print(f"coupled : dp {dpMeasured:.4g} Pa, q_X1 "
          f"{last['pipes']['X1']:.6g} m3/s, q_X2 {last['pipes']['X2']:.6g} m3/s")
    require(abs(dpMeasured - dpExact) <= 0.05 * dpExact,
            "end-to-end pressure within 5% of the series-parallel closed form")
    require(abs(last["pipes"]["X1"] - qX1Exact) <= 0.03 * qX1Exact,
            "the split across the detailed branch is the closed-form one")
    for e in EDGES:
        st = last["settled"][e]
        r = value(st, "R")
        require(0.94 <= r / rDuct[e] <= 1.06,
                f"{e}: its block reports its own duct's resistance "
                f"({r / rDuct[e]:.4f} of it)")

    # ---- 7: each block's own field, at each block's own rate ----
    for e in EDGES:
        L, H, W = EDGES[e]
        ln = last["blocks"][e]
        q, plug = value(ln, "fluxIn"), value(ln, "plugSpeed")
        require(value(ln, "blockConverged") == 1, f"{e}: its block solve converged")
        require(abs(value(ln, "inletArea") * plug - q) <= 1e-9 * q,
                f"{e}: plug speed times inlet area is the carried rate")
        require(abs(value(ln, "fluxOut") - q) <= 1e-3 * q,
                f"{e}: what enters the region leaves it")
        spread = value(ln, "planeFluxMax") - value(ln, "planeFluxMin")
        require(spread <= 1e-3 * q, f"{e}: every plane inside carries one rate")
        slope, dpBlock = abs(value(ln, "probeSlope")), value(ln, "dpPa")
        require(abs(slope * L - dpBlock) <= 0.02 * dpBlock,
                f"{e}: the developed gradient integrates to its own drop")
        require(spread <= 1e-3 * q and value(ln, "probeSections") >= 20,
                f"{e}: enough developed sections inside the probe window")
        require(abs(value(ln, "uxCentre") / plug - 1.5) <= 0.04,
                f"{e}: centreline speed is 3/2 the mean")

        cells = read_vtu_cells(os.path.join(workdir, f"out_multi_{e}.vtu"))
        dev = [c for c in cells if 0.35 * L <= c[0] <= 0.95 * L]
        require(len(dev) >= 20 * H, f"{e}: developed columns present in the .vtu")
        peak = 1.5 * plug

        def profile(y, H=H, peak=peak):
            return peak * (1.0 - ((y - H / 2.0) / (H / 2.0)) ** 2)

        profErr = max(abs(c[2] - profile(c[1])) for c in dev)
        vMax = max(abs(c[3]) for c in dev)
        # a column of cells tiles the cross-section, whatever its own split is
        columns = {}
        for c in dev:
            columns.setdefault(round(c[0], 6), []).append(c[2])
        fluxCol = min(sum(us) / len(us) * H * W for us in columns.values())
        require(profErr <= 0.06 * peak,
                f"{e}: profile is parabolic in its own gap (err {profErr:.4g})")
        require(vMax <= 0.03 * plug, f"{e}: no cross-stream velocity ({vMax:.3g} m/s)")
        require(abs(fluxCol - q) <= 0.02 * q,
                f"{e}: the written field carries its rate column by column")

        # ---- 8: laminar, which is what makes 6 linear ----
        re_ = 4.0 * RHO * q / (math.pi * MU * PIPE_D[e])
        require(5.0 <= re_ <= 200.0, f"{e}: Re {re_:.1f} stayed laminar")

    if fails:
        raise SystemExit(f"{len(fails)} multi-edge coupled checks FAILED")
    print("multi-edge coupled 1D<->3D verification passed")


if __name__ == "__main__":
    main()
