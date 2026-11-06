#!/usr/bin/env python3
"""Run chflow-coupled-1d3d against a RANS detailed region and verify the
partitioned interface whose partner is genuinely nonlinear.

The case is a steady network, Src -> A -> { X (3D block) || By (1D bypass) } ->
B -> Snk, driven by an injected rate, where X is a plane channel solved with
k-omega SST. The same starting coefficient is run twice: once with the plain
interface Gauss-Seidel and once with the quasi-Newton operator, so the two pass
counts differ only in how the next coefficient was chosen.

Six families of statements:

  1. interface identity   - the rate the network delivers is the rate the block
                            was measured at, and the residual is the gap
                            between the two solvers' drops
  2. partner nonlinearity - the block's own dp(q) pairs, read straight off
                            the interface passes, fit a power law with exponent
                            ~1.74; the same fit over the laminar case's passes
                            gives 1.06, and its dp/q moves by 1.06x rather than
                            the 2.1x measured here
  3. wall layer           - the first cell sits in the wall-function band the
                            model is valid in, and the eddy viscosity swamps
                            the molecular one
  4. inside the block     - the developed pressure gradient balances the wall
                            shear the model reports and every cross-section
                            carries the rate that entered: identities of the
                            block's own fields, no 1D side
  5. friction law         - the bulk speed the block carries matches Pope's
                            smooth-wall law integrated with the block's own
                            u_tau, and the profile is flatter than a parabola
  6. the two arms         - same seed, same operating point; the accelerated arm
                            gets there in no more passes, each pass being one
                            block solve

usage: check_coupled_turbulent.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# geometry and fluid of tests/cases/coupledTurbulent/case.dict
MU, RHO = 1.0e-3, 1000.0
NU = MU / RHO
H, W, L = 0.02, 0.02, 0.6       # channel gap (wall-normal), span, length
AREA = H * W
DH = 2.0 * H                    # plane channel: the two walls carry the shear
Q_IN = 1.0e-3                   # m3/s injected at Src
TOL_INTERFACE = 1e-3            # couple{} tolerance the case asks for
ARMS = ("0", "1")               # seed indices: plain then accelerated


def value(line, key):
    """First number stored under `key` in a log line ('q 1.5e-4 ...')."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([-+0-9.eE]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return float(m.group(1))


def word(line, key):
    """The token stored under `key` in a log line ('turb kOmegaSST ...')."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+(\S+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return m.group(1)


def run(binary, case):
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("turbulent coupled solver failed")
    return r.stdout, workdir


def popeBulkSpeed(uTau, nu, dh):
    """Pope's smooth-wall friction law, 1/sqrt(f) = 2 log10(Re_b sqrt(f)) - 0.6,
    solved for the bulk speed that carries a given wall shear: fixed point on
    U = u_tau sqrt(8/f)."""
    u = uTau * 20.0
    for _ in range(200):
        f = 8.0 * (uTau / u) ** 2
        inv = 2.0 * math.log10(u * dh / nu * math.sqrt(f)) - 0.6
        u = uTau * math.sqrt(8.0) * inv
    return u


def powerLawFit(pairs):
    """Least-squares exponent of dp = c q^m over (q, dp) pairs."""
    xs = [math.log(q) for q, dp in pairs]
    ys = [math.log(dp) for q, dp in pairs]
    n = float(len(xs))
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def main():
    binary, case = sys.argv[1], sys.argv[2]
    out, workdir = run(binary, case)
    lines = [ln.replace("[chflow] ", "") for ln in out.splitlines()]

    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    arms = {}
    for ln in lines:
        if ln.startswith("coupled seed "):
            arms[ln.split()[2].rstrip(":")] = {"final": ln}
        elif ln.startswith("block seed "):
            arms[ln.split()[2].rstrip(":")]["block"] = ln
        elif ln.startswith("network seed "):
            arms[ln.split()[2].rstrip(":")]["network"] = ln
        elif ln.startswith("interface seed "):
            key = ln.split()[2]
            arms.setdefault(key, {}).setdefault("passes", []).append(ln)

    require(set(arms) == set(ARMS), f"both interface arms ran (got {sorted(arms)})")

    pipes = {}
    for ln in lines:
        if ln.startswith("pipe final "):
            pipes[ln.split()[2]] = ln
    require(set(pipes) == {"A", "X", "By", "B"}, "all four pipes reported")

    finals = {}
    for name in ARMS:
        a = arms[name]
        fin, blk, net = a["final"], a["block"], a["network"]
        converged = "converged" in fin and "NOT" not in fin
        iters = value(fin, "iterations")
        q, dp, R, res = value(fin, "q"), value(fin, "dp"), value(fin, "R"), value(fin, "residual")
        acc = value(fin, "accelerate")
        require(converged, f"arm {name}: interface converged")
        require(acc == float(name), f"arm {name}: the arm ran with the operator asked for")

        # ---- 1. interface identity --------------------------------------
        require(abs(value(net, "qNetwork") - q) <= 1e-12 * abs(q),
                f"arm {name}: the rate the network delivered is the rate measured")
        require(abs(value(net, "dpNetwork") - R * q) <= 1e-9 * abs(dp),
                f"arm {name}: what the network charges is its coefficient "
                f"times the rate it delivered")
        last = a["passes"][-1]
        require(abs(value(last, "dpBlock") - dp) <= 1e-12 * abs(dp),
                f"arm {name}: the drop reported is the partner's own, at that rate")
        gap = abs(R * q - value(last, "dpBlock")) / abs(value(last, "dpBlock"))
        require(abs(gap - res) <= 1e-12 * max(res, 1.0),
                f"arm {name}: residual is the gap between the two solvers")
        require(res <= TOL_INTERFACE,
                f"arm {name}: residual {res:.3e} at or below the asked {TOL_INTERFACE:.0e}")
        require(value(net, "nodeResidual") < 1e-9,
                f"arm {name}: the 1D side solved tight, so its rate is not the "
                f"loose end of the comparison")
        require(all(value(p, "residual") >= 0 for p in a["passes"]),
                f"arm {name}: every pass reports a residual")

        # ---- block-side numbers -----------------------------------------
        require(word(blk, "turb") == "kOmegaSST",
                f"arm {name}: the block reports the model it solved with")
        blkq = {"fluxIn": value(blk, "fluxIn"), "fluxOut": value(blk, "fluxOut"),
                "planeMin": value(blk, "planeFluxMin"), "planeMax": value(blk, "planeFluxMax"),
                "plug": value(blk, "plugSpeed"),
                "dpPa": value(blk, "dpPa"), "slope": value(blk, "probeSlope"),
                "uxCentre": value(blk, "uxCentre"), "nutMax": value(blk, "nutMax"),
                "yPlus": value(blk, "wallYPlusMax"), "uTau": value(blk, "wallUTau"),
                "sweeps": value(blk, "sweeps")}
        require("blockConverged 1" in blk,
                f"arm {name}: the block is converged at the rate the interface "
                f"settled on")
        finals[name] = (q, dp, R, res, iters, blkq, a["passes"])

    # ---- 2. the partner is nonlinear, and the passes show it -------------
    for name in ARMS:
        q, dp, R, res, iters, blkq, passes = finals[name]
        pairs = [(value(p, "q"), value(p, "dpBlock")) for p in passes]
        qs = [qq for qq, _ in pairs]
        require(len(pairs) >= 4 and max(qs) / min(qs) > 1.5,
                f"arm {name}: the passes visit rates spread over "
                f"{max(qs) / min(qs):.2f}x, enough to measure dp(q)")
        m = powerLawFit(pairs)
        require(1.40 < m < 2.15,
                f"arm {name}: the block's own dp = c q^m fit reads m = {m:.3f} "
                f"(the same fit over the laminar case's passes gives 1.06)")
        apparent = [dpBlock / qq for qq, dpBlock in pairs]
        spread = max(apparent) / min(apparent)
        require(spread > 1.8,
                f"arm {name}: the drop-per-rate the partner reports moves by "
                f"{spread:.2f}x across the rates it was asked for (the laminar "
                f"block's own passes move it by 1.06x)")

    # ---- 3. wall layer and eddy viscosity -------------------------------
    for name in ARMS:
        blkq = finals[name][5]
        require(30.0 < blkq["yPlus"] < 300.0,
                f"arm {name}: first cell y+ = {blkq['yPlus']:.1f} is inside the "
                f"wall-function band")
        require(blkq["nutMax"] / NU > 50.0,
                f"arm {name}: largest eddy viscosity is {blkq['nutMax'] / NU:.1f}x "
                f"the molecular one")
        require(blkq["uxCentre"] / blkq["plug"] < 1.25,
                f"arm {name}: centreline speed is {blkq['uxCentre'] / blkq['plug']:.3f}x "
                f"the bulk speed, flatter than the 1.5 a parabola gives")

    # ---- 4. momentum closure inside the block ---------------------------
    for name in ARMS:
        blkq = finals[name][5]
        tauBalance = abs(blkq["slope"]) * H / 2.0
        tauModel = RHO * blkq["uTau"] ** 2
        require(abs(tauBalance - tauModel) / tauModel < 0.02,
                f"arm {name}: developed gradient times H/2 balances the reported "
                f"wall shear ({tauBalance:.4f} vs {tauModel:.4f} Pa)")
        fluxIn, fluxOut = blkq["fluxIn"], blkq["fluxOut"]
        require(abs(fluxOut - fluxIn) / abs(fluxIn) < 5e-3,
                f"arm {name}: what enters the region leaves it to "
                f"{abs(fluxOut - fluxIn) / abs(fluxIn):.1e}")
        planeWorst = max(abs(blkq["planeMin"] - fluxIn),
                         abs(blkq["planeMax"] - fluxIn)) / abs(fluxIn)
        require(planeWorst < 5e-3,
                f"arm {name}: every cross-section plane inside carries the same "
                f"rate to {planeWorst:.1e}")

    # ---- 5. friction law ------------------------------------------------
    for name in ARMS:
        blkq = finals[name][5]
        uPred = popeBulkSpeed(blkq["uTau"], NU, DH)
        rel = (blkq["plug"] - uPred) / uPred
        require(abs(rel) < 0.08,
                f"arm {name}: bulk speed {blkq['plug']:.4f} vs Pope's "
                f"{uPred:.4f} m/s from the block's own u_tau ({100 * rel:+.1f}%)")

    # ---- 6. the two arms, and the cost of each -------------------------
    q0, dp0, R0 = finals["0"][0], finals["0"][1], finals["0"][2]
    q1, dp1, R1 = finals["1"][0], finals["1"][1], finals["1"][2]
    it0, it1 = finals["0"][4], finals["1"][4]
    require(abs(R1 - R0) / R0 < 3.0 * TOL_INTERFACE,
            f"same seed, same interface resistance from both arms "
            f"({R0:.6e} vs {R1:.6e})")
    require(abs(q1 - q0) / q0 < 3.0 * TOL_INTERFACE,
            f"both arms carry the same rate through the block "
            f"({q0:.6e} vs {q1:.6e})")
    require(it1 <= it0,
            f"the accelerated arm needs no more passes than the plain one "
            f"({int(it1)} vs {int(it0)})")
    for name in ARMS:
        passes = finals[name][6]
        require(all(value(p, "sweeps") >= 1 for p in passes),
                f"arm {name}: every pass actually solved the block")
        require(abs(finals[name][5]["dpPa"] - value(passes[-1], "dpBlock"))
                <= 1e-12 * abs(value(passes[-1], "dpBlock")),
                f"arm {name}: the block's own drop is the number the interface used")
    # one pass is one block solve: the accelerated arm wins on fewer of them, not
    # on more, so the total number of SIMPLE starts is the sum of the two budgets
    solves = sum(1 for ln in lines if ln.startswith("SIMPLE iter 1:"))
    require(solves == sum(int(finals[n][4]) for n in ARMS),
            f"block solves equal the passes of both arms and nothing else "
            f"({solves} == {sum(int(finals[n][4]) for n in ARMS)})")

    # ---- the 1D side around the block ----------------------------------
    qX = value(pipes["X"], "Q")
    qBy = value(pipes["By"], "Q")
    require(abs(value(pipes["A"], "Q") - Q_IN) <= 1e-9 * Q_IN,
            "the lead carries the whole injected rate")
    require(abs(qX + qBy - Q_IN) <= 1e-9 * Q_IN, "the branch split adds up to the rate")
    require(abs(value(pipes["X"], "dp") - value(pipes["By"], "dp")) <= 1e-9 * dp1,
            "the parallel branches share their drop")
    require(value(pipes["By"], "Re") > 4000.0,
            f"the bypass is turbulent too (Re {value(pipes['By'], 'Re'):.0f}), so its "
            f"1D law is not a laminar shortcut")

    print()
    if fails:
        print(f"{len(fails)} FAILED CHECKS:")
        for f in fails:
            print("  - " + f)
        raise SystemExit(1)
    print("turbulent-partner coupled 1D<->3D verification passed")


if __name__ == "__main__":
    main()
