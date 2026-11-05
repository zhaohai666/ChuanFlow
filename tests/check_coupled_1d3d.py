#!/usr/bin/env python3
"""Run chflow-coupled-1d3d and verify the partitioned 1D<->3D interface.

The case is a steady network, Src -> A -> { X (3D block) || By (1D bypass) } ->
B -> Snk, driven by an injected rate. Six families of statements are checked,
per seed:

  1. interface identity   - the rate the network delivers is the rate the block
                            was measured at, and the residual is the gap
                            between two solvers rather than a number compared
                            with itself
  2. node continuity      - the leads carry the whole rate and the branch split
                            adds up
  3. pressure-chain       - the drops around the loop sum to the end-to-end
                            difference, and the parallel branches share theirs
  4. block continuity     - what enters the detailed region leaves it, and every
                            cross-section plane inside carries the same rate
  5. developed gradient   - inside the block, dp/dx is the planar Poiseuille
                            value and the profile is the parabola
  6. operating point      - the closed-form network solution with analytic
                            laminar resistances reproduces the coupled result

plus seed independence: three starting coefficients, one answer.

usage: check_coupled_1d3d.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# geometry and fluid of tests/cases/coupled1d3d/case.dict
MU, RHO = 0.1, 800.0
H, W, L = 0.02, 0.02, 0.5      # duct height (wall-normal), width, length
AREA = H * W
Q_IN = 2.4e-4                  # m3/s injected at Src
P_REF = 1.0e6                  # Pa held at the sink
TOL_INTERFACE = 1e-4           # couple{} tolerance the case asks for


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

    seeds = {}
    for ln in lines:
        if ln.startswith("coupled seed "):
            key = ln.split()[2].rstrip(":")
            seeds[key] = {"final": ln}
        elif ln.startswith("block seed "):
            seeds[ln.split()[2].rstrip(":")]["block"] = ln
        elif ln.startswith("network seed "):
            seeds[ln.split()[2].rstrip(":")]["network"] = ln

    require(len(seeds) == 3, f"three interface seeds ran (got {len(seeds)})")

    nodes = [float(ln.split()[-1]) for ln in lines if ln.startswith("node final p ")]
    pipes = {}
    for ln in lines:
        if ln.startswith("pipe final "):
            t = ln.split()
            pipes[t[2]] = {"Q": value(ln, "Q"), "dp": value(ln, "dp"),
                           "Re": value(ln, "Re")}
    require(len(nodes) == 4, "four node pressures reported")
    require(set(pipes) == {"A", "X", "By", "B"}, "all four pipes reported")

    # ---- per-seed: interface identity, block continuity, convergence ----
    for name in sorted(seeds, key=int):
        s = seeds[name]
        fin, blk, net = s["final"], s["block"], s["network"]
        require("converged" in fin and "NOT" not in fin, f"seed {name}: interface converged")
        q, dpBlock, R = value(fin, "q"), value(fin, "dp"), value(fin, "R")
        res = value(fin, "residual")
        iters = value(fin, "iterations")
        dpNet, qNet = value(net, "dpNetwork"), value(net, "qNetwork")
        fluxIn, fluxOut = value(blk, "fluxIn"), value(blk, "fluxOut")
        pMin, pMax = value(blk, "planeFluxMin"), value(blk, "planeFluxMax")

        require(value(blk, "blockConverged") == 1,
                f"seed {name}: the block solve converged at the reported state")
        require(value(net, "nodeResidual") < 1e-9,
                f"seed {name}: the 1D network solved tight")
        # (1) one state: the block measured the rate the network delivers
        require(abs(qNet - fluxIn) <= 1e-12 * q,
                f"seed {name}: the block carries exactly the network's branch rate")
        require(abs(qNet - q) <= 1e-12 * q, f"seed {name}: reported rate is consistent")
        # the coupled element's own law, as the network solved it
        require(abs(dpNet - R * qNet) <= 1e-9 * abs(dpNet),
                f"seed {name}: dpNetwork equals resistance times rate")
        # the residual is a gap between two solvers, and it is what it says
        require(abs(res - abs(dpNet - dpBlock) / dpBlock) <= 1e-9,
                f"seed {name}: residual is the cross-solver gap")
        require(res <= TOL_INTERFACE,
                f"seed {name}: interface residual {res:.3e} below the case tolerance")
        require(1 <= iters <= 6, f"seed {name}: converged in {int(iters)} passes")
        # (4) what enters leaves, and every plane inside carries the same rate
        require(abs(fluxOut - fluxIn) <= 1e-7 * fluxIn,
                f"seed {name}: the block conserves volume")
        require(pMax - pMin <= 1e-7 * q,
                f"seed {name}: every cross-section plane inside carries one rate")

    # ---- (2) node continuity and (3) pressure chain, from the final state ----
    pSrc, pJ1, pJ2, pSnk = nodes
    require(abs(pipes["A"]["Q"] - Q_IN) <= 1e-9 * Q_IN,
            "the lead carries the injected rate")
    require(abs(pipes["A"]["Q"] - pipes["X"]["Q"] - pipes["By"]["Q"])
            <= 1e-9 * Q_IN, "upstream node: in = coupled + bypass")
    require(abs(pipes["B"]["Q"] - pipes["X"]["Q"] - pipes["By"]["Q"])
            <= 1e-9 * Q_IN, "downstream node: out = coupled + bypass")
    chain = pipes["A"]["dp"] + pipes["X"]["dp"] + pipes["B"]["dp"]
    require(abs(chain - (pSrc - pSnk)) <= 1e-9 * abs(pSrc - pSnk),
            "the pressure chain closes across the interface")
    require(abs(pipes["X"]["dp"] - pipes["By"]["dp"]) <= 1e-9 * pipes["X"]["dp"],
            "parallel branches share their drop")
    require(abs((pJ1 - pJ2) - pipes["X"]["dp"]) <= 1e-9 * pipes["X"]["dp"],
            "the branch drop is the interface node difference")
    require(abs(pSnk - P_REF) < 1e-9, "the sink holds its reference pressure")

    # the seed that wrote the field is the one the final state describes
    lastBlk = seeds[max(seeds, key=int)]["block"]
    plug = value(lastBlk, "plugSpeed")
    dpBlock = value(lastBlk, "dpPa")
    slope = value(lastBlk, "probeSlope")
    spread = value(lastBlk, "probeSpread")
    qBlock = value(lastBlk, "fluxIn")

    # ---- (5) developed gradient and profile ----
    gradExact = 12.0 * MU * plug / (H * H)
    require(abs(abs(slope) - gradExact) <= 0.03 * gradExact,
            f"developed dp/dx {abs(slope):.4g} vs planar Poiseuille {gradExact:.4g}")
    require(abs(value(lastBlk, "inletArea") * plug - qBlock) <= 1e-12 * qBlock,
            "plug speed times inlet area is the carried rate")
    require(spread <= 0.002 * abs(dpBlock),
            f"section means hug the developed line (spread {spread:.3g} Pa)")
    require(value(lastBlk, "probeSections") >= 20, "enough sections inside the probe window")
    require(abs(value(lastBlk, "uxCentre") / plug - 1.5) <= 0.03,
            "centreline speed is 3/2 the mean")

    cells = read_vtu_cells(os.path.join(workdir, "out_coupled_block.vtu"))
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
    require(abs(fluxCol - qBlock) <= 0.02 * qBlock,
            "the written field carries the branch rate column by column")

    # ---- (6) analytic operating point of the whole coupled system ----
    rBlock = 12.0 * MU * L / (H ** 3 * W)
    rBy, rLead = hagenPoiseuille(1.5, 0.03), hagenPoiseuille(5.0, 0.1)
    rBranch = 1.0 / (1.0 / rBlock + 1.0 / rBy)
    dpExact = Q_IN * (2.0 * rLead + rBranch)
    qBlockExact = Q_IN * rBy / (rBlock + rBy)
    measuredR = dpBlock / qBlock
    print(f"analytic: dp {dpExact:.4g} Pa, coupled block Q {qBlockExact:.6g} m3/s, "
          f"block R {rBlock:.4g} Pa.s/m3")
    print(f"coupled : dp {pSrc - pSnk:.4g} Pa, block Q {pipes['X']['Q']:.6g} m3/s, "
          f"block R {measuredR:.4g} Pa.s/m3")
    require(abs((pSrc - pSnk) - dpExact) <= 0.04 * dpExact,
            "end-to-end pressure within 4% of the closed form")
    require(abs(pipes["X"]["Q"] - qBlockExact) <= 0.04 * qBlockExact,
            "branch split within 4% of the closed form")
    require(1.0 <= measuredR / rBlock <= 1.06,
            "the block's resistance is the laminar duct value plus its entrance")
    require(20.0 <= pipes["X"]["Re"] <= 200.0,
            "the detailed region stayed laminar, as the analytic check assumes")

    # ---- (7) the fixed point does not belong to the starting coefficient ----
    qs = [value(seeds[k]["final"], "q") for k in seeds]
    ds = [value(seeds[k]["final"], "dp") for k in seeds]
    rs = [value(seeds[k]["final"], "R") for k in seeds]
    for label, vals in (("rate", qs), ("drop", ds), ("resistance", rs)):
        span = (max(vals) - min(vals)) / max(vals)
        require(span <= 1e-4, f"seeds agree on the interface {label} to {span:.2e}")

    if fails:
        raise SystemExit(f"{len(fails)} coupled-interface checks FAILED")
    print("coupled 1D<->3D verification passed")


if __name__ == "__main__":
    main()
