#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient-multi against TWO RANS detailed regions and
verify a time-stepping network whose interface is walked once plain and once
accelerated.

The case puts on one table the three ingredients that each had their own
increment: a nonlinear partner (§5.23), two edges inside the same time level
(§5.25) and the stepping network. Src -> A -> { X1 (3D k-omega SST duct) || By
(1D bypass) } -> X2 (3D k-omega SST duct) -> B -> Snk, from rest under a rate
switched on at t = 0.

What a second nonlinear edge could have done but does not: X1 is the parallel
branch, so its rate is a result of the momentum balance and the ramp toward the
resistive split is set by the sum of the two partners' SLOPES; X2 is in series,
so the incompressible network hands it the whole injected rate and there is no
rate left to negotiate - once its field is warm its coefficient stops moving
too. Ten families:

  1. two operators, one case  - both arms quote the same starting coefficient for
                                both edges and differ only in `accelerates`
  2. one state per edge       - every step reports both edges, each with the rate its
                                own pipe carries, R times that rate as the drop the
                                network charges, and the residual as the solvers' gap
  3. the proxy bore is exact  - d = 2*sqrt(H*W/pi) makes the installed rho*L/A the
                                number the block measures on its own mesh, to 4e-15,
                                so attribution is not a variable in this case
  4. the 1D network           - every node drop closes as R*Q + M*dQ/dt with the block's
                                own turbulent chord, the bypass against published
                                Churchill, the series edge including its sudden-start kick
  5. both partners nonlinear  - the parallel block's dp(q) across the ramp fits
                                exponent 1.68, and both chords are 5.8x and 7.8x the
                                planar Poiseuille their geometry would give - the same
                                two blocks read 0.98x laminar in §5.25
  6. the series edge is free  - its rate is pinned to the injected value to 1e-16, its
                                coefficient stops moving after step 1, its residual runs
                                2-5 orders under the tolerance, and the interface's
                                `worst` edge is X1 on every single pass
  7. the schedule             - one pass refreshes both blocks, one pass per time level is
                                accepted, block solves == passes x 2, and the budget is
                                never reached
  8. the ramp is slopes       - (M_X1 + M_By)/(m_X R_X + m_By R_By) against the fitted
                                time constant, and the chord formula - exact for the
                                laminar pair - reading it 64% high
  9. two operators, one answer - plain Gauss-Seidel and quasi-Newton land on one state per
                                edge, and on the pinned edge they report the same digits
 10. each block's field       - conservation, wall resolution, shear-against-gradient
                                closure and Pope's friction law, in each block's own duct

usage: check_transient_coupled_multi_turbulent.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# fluid and geometry of tests/cases/transientCoupledMultiTurbulent1d3d/case.dict
MU, RHO = 1.0e-3, 1000.0
NU = MU / RHO
Q_IN = 1.2e-4                   # m3/s injected at Src at t = 0
P_REF = 1.0e6                   # Pa held at Snk
DT, N_STEPS = 0.05, 12          # time { dt 0.05; endTime 0.6; }
MAX_PASSES = 6                  # couple{ maxSubIterations }
TOL_INTERFACE = 2e-4            # couple{ tolerance }
ARMS = ("gs", "iqn")            # couple{ starts } with accelerates ( 0 1 )
SEED = {"X1": 1.0e7, "X2": 1.5e7}

EDGES = {"X1": (0.2, 0.005, 0.005),       # length, wall-normal gap, span
         "X2": (0.12, 0.004, 0.008)}
BORE = {"X1": 0.0056418958354775631,      # 2*sqrt(H*W/pi) of EDGES["X1"]
        "X2": 0.00638307648642292245}     # 2*sqrt(H*W/pi) of EDGES["X2"]
BY = (0.3, 0.008)               # the 1D bypass, sharing X1's node pair
LEAD = (0.3, 0.02)              # A and B, carrying everything
NAME = ("X1", "X2")
YPLUS_LOG_LAW = 11.81           # the wall functions are valid above this


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


def word(line, key):
    """The token stored under `key` in a log line ('turb kOmegaSST ...')."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+(\S+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return m.group(1)


def churchillFriction(re, roughness=0.0):
    """Churchill's 1977 Darcy friction factor, all three flow regimes - the
    published correlation the 1D bypass's momentum row is closed with, so the
    network's own friction law is graded against the literature rather than
    against itself."""
    re = max(re, 1e-12)
    a = (2.457 * math.log(1.0 / ((7.0 / re) ** 0.9 + 0.27 * roughness))) ** 16
    b = (37530.0 / re) ** 16
    return 8.0 * ((8.0 / re) ** 12 + 1.0 / (a + b) ** 1.5) ** (1.0 / 12.0)


def pipeChord(length, diameter, q):
    """dp/q of a 1D pipe carrying q, straight from Darcy-Weisbach."""
    area = 0.25 * math.pi * diameter ** 2
    v = q / area
    re = RHO * v * diameter / MU
    return churchillFriction(re) * length / diameter * 0.5 * RHO * v * v / q


def planarDuct(length, gap, width):
    """Laminar resistance of a plane channel: q = gap^3 width dp / (12 mu L)."""
    return 12.0 * MU * length / (gap ** 3 * width)


def popeBulkSpeed(uTau, nu, dh):
    """Pope's smooth-wall friction law, 1/sqrt(f) = 2 log10(Re_b sqrt(f)) - 0.6,
    solved for the bulk speed a given wall shear carries."""
    u = uTau * 20.0
    for _ in range(200):
        f = 8.0 * (uTau / u) ** 2
        u = uTau * math.sqrt(8.0) * (2.0 * math.log10(u * dh / nu * math.sqrt(f)) - 0.6)
    return u


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("transient multi-edge turbulent coupled solver failed")

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
    arms, cur = [], None
    for ln in lines:
        if ln.startswith("transient start ") and " inertia " in ln:
            cur = {"tag": ln, "seed": {}, "step": {}, "pass": {}, "edge": {},
                   "nodes": {}, "pipes": {}, "block": [], "summary": [], "solves": 0}
            arms.append(cur)
        elif cur is None:
            continue
        if ln.startswith("SIMPLE iter 1:"):
            cur["solves"] += 1
        elif ln.startswith("transient start ") and " index " in ln:
            cur["seed"][ln.split()[4]] = ln
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
            cur["block"].append(ln)
        elif ln.startswith("transient coupled start ") and " name " in ln:
            cur["summary"].append(ln)

    # ---- 1. two operators, one network, one seed pair ----
    require(len(arms) == len(ARMS), f"both arms ran (got {len(arms)})")
    tags = [c["tag"].split()[2].rstrip(":") for c in arms]
    require(tags == list(ARMS), f"the arms announce themselves in order: {' '.join(tags)}")
    require([("accelerate yes" in c["tag"]) for c in arms] == [False, True],
            "the two runs differ in one thing and one thing only: the second interface"
            " is the quasi-Newton, the first is plain block Gauss-Seidel")
    require(all("inertia network" in c["tag"] for c in arms),
            "both arms read the region's inertia from the 1D side, which family 3 makes"
            " the same number the block measures on its own mesh")
    seed = [{k: value(c["seed"][k], "resistance") for k in NAME} for c in arms]
    require(all(s == SEED for s in seed),
            f"both arms seed the same coefficient for both edges at once: X1"
            f" {SEED['X1']:.3g}, X2 {SEED['X2']:.3g}")
    for a, c in enumerate(arms):
        require(sorted(c["step"]) == list(range(1, N_STEPS + 1)),
                f"{ARMS[a]}: all {N_STEPS} time levels reported a state")
        require(sorted(c["nodes"]) == sorted(c["step"])
                and sorted(c["pipes"]) == sorted(c["step"]),
                f"{ARMS[a]}: every step reports its whole 1D state")
        require(len(c["block"]) == 2 and len(c["summary"]) == 2,
                f"{ARMS[a]}: both blocks and both per-edge interface states reported once")
        require(abs(value(c["step"][N_STEPS], "t") - 0.6) <= 1e-12,
                f"{ARMS[a]}: the series ends at the endTime the case asks for")
        require(all(word(c["block"][i], "turb") == "kOmegaSST" for i in (0, 1)),
                f"{ARMS[a]}: both partners that faced the interface were RANS blocks")

    # ---- (2) one state per edge, and (3) what inertia it is referred to ----
    q = [{} for _ in arms]
    R = [{} for _ in arms]
    M = [{} for _ in arms]
    dpB = [{} for _ in arms]
    dpN = [{} for _ in arms]
    res = [{} for _ in arms]
    sw = [{} for _ in arms]
    w = [{"chord": 0.0, "res": 0.0, "tol": 0.0, "rate": 0.0, "newton": 0.0,
          "inertia": 0.0, "blockInertia": 0.0} for _ in arms]
    for a, c in enumerate(arms):
        for i, name in enumerate(NAME):
            steps = sorted(c["edge"])
            q[a][name] = [value(c["edge"][n][i], "q") for n in steps]
            R[a][name] = [value(c["edge"][n][i], "R") for n in steps]
            M[a][name] = [value(c["edge"][n][i], "M") for n in steps]
            dpB[a][name] = [value(c["edge"][n][i], "dpBlock") for n in steps]
            dpN[a][name] = [value(c["edge"][n][i], "dpNetwork") for n in steps]
            res[a][name] = [value(c["edge"][n][i], "residual") for n in steps]
            sw[a][name] = [int(value(c["edge"][n][i], "sweeps")) for n in steps]
        for n in sorted(c["step"]):
            check(abs(value(c["step"][n], "t") - n * DT) <= 1e-12,
                  f"{ARMS[a]} step {n}: t is {n} steps of {DT} s")
            check(value(c["step"][n], "newton") < 1e-11,
                  f"{ARMS[a]} step {n}: the 1D step solved its own system")
            w[a]["newton"] = max(w[a]["newton"], value(c["step"][n], "newton"))
            for i, name in enumerate(NAME):
                e = c["edge"][n][i]
                x, rr = value(e, "q"), value(e, "R")
                check(value(e, "blockConverged") == 1,
                      f"{ARMS[a]} {name} step {n}: the block converged")
                check(abs(value(c["pipes"][n], name) - x) <= 1e-16 * Q_IN,
                      f"{ARMS[a]} {name} step {n}: the edge reports the rate its own"
                      " pipe carries")
                check(abs(value(e, "dpNetwork") - rr * x) <= 1e-12 * value(e, "dpNetwork"),
                      f"{ARMS[a]} {name} step {n}: the network charges R times the rate")
                gap = abs(value(e, "dpNetwork") - value(e, "dpBlock")) / value(e, "dpBlock")
                check(abs(gap - value(e, "residual")) <= 1e-12,
                      f"{ARMS[a]} {name} step {n}: residual is the gap between two solvers")
                check(value(e, "residual") <= TOL_INTERFACE,
                      f"{ARMS[a]} {name} step {n}: interface residual"
                      f" {value(e, 'residual'):.3e} under the case tolerance")
                w[a]["chord"] = max(w[a]["chord"],
                                    abs(value(e, "dpNetwork") - rr * x)
                                    / value(e, "dpNetwork"))
                w[a]["res"] = max(w[a]["res"], abs(gap - value(e, "residual")))
                w[a]["tol"] = max(w[a]["tol"], value(e, "residual"))
                w[a]["rate"] = max(w[a]["rate"], abs(value(c["pipes"][n], name) - x) / x)
                L, H, W = EDGES[name]
                bore = RHO * L / (0.25 * math.pi * BORE[name] ** 2)
                region = RHO * L / (H * W)
                w[a]["inertia"] = max(w[a]["inertia"],
                                      max(abs(value(e, "M") / bore - 1.0),
                                          abs(value(e, "M") / region - 1.0)))
                w[a]["blockInertia"] = max(w[a]["blockInertia"],
                                           abs(value(c["block"][i], "inertance") / region - 1.0))

    for a in range(len(arms)):
        require(w[a]["chord"] <= 1e-12 and w[a]["res"] <= 1e-12
                and w[a]["rate"] <= 1e-15,
                f"{ARMS[a]}: across all 24 edge-steps the interface is still a chord -"
                f" the drop the network charges is R*Q to {w[a]['chord']:.1e}, the"
                f" residual is the two solvers' gap to {w[a]['res']:.1e}, and each block"
                f" is fed its own pipe's rate to {w[a]['rate']:.1e}")
        require(w[a]["tol"] <= TOL_INTERFACE,
                f"{ARMS[a]}: every accepted state sits under the tolerance the case asks"
                f" for (worst {w[a]['tol']:.3e})")
        require(w[a]["newton"] < 1e-11,
                f"{ARMS[a]}: every 1D time level solved its own network to"
                f" {w[a]['newton']:.1e}")
        require(w[a]["inertia"] <= 1e-14 and w[a]["blockInertia"] <= 1e-14,
                f"{ARMS[a]}: one bore for both regions: the installed inertance is rho*L/(pi d^2/4)"
                f" and rho*L/(H*W) to {w[a]['inertia']:.1e}, and the block's own"
                f" rho*V/A^2 is the same number to {w[a]['blockInertia']:.1e} - the two"
                " sources of §5.24 cannot be told apart here, so what is left to measure"
                " is the partner's law")

    # ---- 4. the 1D rows, closed with what each partner actually reports ----
    M_BY = RHO * BY[0] / (0.25 * math.pi * BY[1] ** 2)
    M_LEAD = RHO * LEAD[0] / (0.25 * math.pi * LEAD[1] ** 2)
    qBy = [{} for _ in arms]
    for a, c in enumerate(arms):
        worst = {k: 0.0 for k in ("X1", "X2", "By", "A", "B")}
        kick = 0.0
        prev = {}
        for n in sorted(c["step"]):
            nd, pp = c["nodes"][n], c["pipes"][n]
            p = {k: value(nd, k) for k in ("Src", "J1", "J2", "J3", "Snk")}
            f = {k: value(pp, k) for k in ("A", "X1", "By", "X2", "B")}
            qBy[a][n] = f["By"]
            check(abs(p["Snk"] - P_REF) <= 1e-9,
                  f"{ARMS[a]} step {n}: the sink holds its reference pressure")
            check(abs(f["A"] - Q_IN) <= 1e-15 * Q_IN and abs(f["B"] - Q_IN) <= 1e-15 * Q_IN,
                  f"{ARMS[a]} step {n}: both leads carry exactly the injected rate")
            check(abs(f["A"] - f["X1"] - f["By"]) <= 1e-15 * Q_IN
                  and abs(f["X2"] - f["X1"] - f["By"]) <= 1e-15 * Q_IN,
                  f"{ARMS[a]} step {n}: the branch split adds up at both ends of X1")
            for key, drop, flow, resist, inert in (
                    ("X1", p["J1"] - p["J2"], f["X1"], R[a]["X1"][n - 1],
                     M[a]["X1"][n - 1]),
                    ("X2", p["J2"] - p["J3"], f["X2"], R[a]["X2"][n - 1],
                     M[a]["X2"][n - 1]),
                    ("By", p["J1"] - p["J2"], f["By"], pipeChord(*BY, f["By"]), M_BY),
                    ("A", p["Src"] - p["J1"], f["A"], pipeChord(*LEAD, f["A"]), M_LEAD),
                    ("B", p["J3"] - p["Snk"], f["B"], pipeChord(*LEAD, f["B"]), M_LEAD)):
                worst[key] = max(worst[key],
                                 abs(drop - (resist * flow + inert
                                             * (flow - prev.get(key, 0.0)) / DT))
                                 / abs(drop))
            if n == 1:
                want = R[a]["X2"][0] * f["X2"] + M[a]["X2"][0] * f["X2"] / DT
                kick = max(kick, abs((p["J2"] - p["J3"]) - want) / want)
            prev.update(f)
        require(worst["X1"] <= 1e-9 and worst["X2"] <= 1e-12,
                f"{ARMS[a]}: both detailed regions' node drops close as R*Q + M*dQ/dt"
                f" with the RANS chord the block reported on the same pass (worst"
                f" {worst['X1']:.2e} / {worst['X2']:.2e})")
        require(worst["By"] <= 1e-12 and worst["A"] <= 1e-11 and worst["B"] <= 1e-11,
                f"{ARMS[a]}: the bypass and both leads close the same row against"
                f" published Churchill at Re {RHO * qBy[a][N_STEPS] / (0.25 * math.pi * BY[1] ** 2) * BY[1] / MU:.0f}"
                f" - the 1D side is graded against the literature, not itself (worst"
                f" {worst['By']:.1e}/{worst['A']:.1e}/{worst['B']:.1e})")
        require(kick <= 1e-14,
                f"{ARMS[a]}: the series edge negotiates no rate and still pays the sudden"
                f" start in its own momentum row - step 1's node drop is R*Q + M*Q/dt to"
                f" {kick:.1e}, of which {M[a]['X2'][0] * Q_IN / DT:.1f} Pa is the inertial"
                f" kick against {R[a]['X2'][0] * Q_IN:.1f} Pa of friction")

    # ---- 5. both partners are nonlinear, and not by the same amount ----
    lam = {n: planarDuct(*EDGES[n]) for n in NAME}
    for a, c in enumerate(arms):
        mX = math.log(dpB[a]["X1"][-1] / dpB[a]["X1"][0]) \
            / math.log(q[a]["X1"][-1] / q[a]["X1"][0])
        spread = (q[a]["X1"][0] - q[a]["X1"][-1]) / q[a]["X1"][-1]
        require(1.5 < mX < 2.0,
                f"{ARMS[a]}: the parallel block's own dp(q) answers across a"
                f" {spread * 100:.1f}% rate move fit dp = c*q^{mX:.3f} - the laminar"
                " partner of §5.25, measured the same way, gives 1.06")
        require(all(dpB[a]["X1"][i] > dpB[a]["X1"][i + 1] for i in range(N_STEPS - 1))
                and all(q[a]["X1"][i] > q[a]["X1"][i + 1] for i in range(N_STEPS - 1)),
                f"{ARMS[a]}: the negotiated branch relaxes downward in rate and drop"
                " together, as a monotone partner law requires")
        dqX2 = max(abs(v - Q_IN) / Q_IN for v in q[a]["X2"])
        ddpX2 = max(dpB[a]["X2"]) / min(dpB[a]["X2"]) - 1.0
        require(dqX2 <= 1e-15 and ddpX2 <= 1e-5,
                f"{ARMS[a]}: the series block never gets to show a law - its rate holds"
                f" the injected value to {dqX2:.1e} and its drop spans {ddpX2:.1e}"
                " over the whole ramp, so no exponent can be read off this trajectory")
    print("settled chords against the laminar duct of the same geometry:"
          + "".join(f" {n} {R[0][n][-1] / lam[n]:.4f}x" for n in NAME))
    for n in NAME:
        require(R[0][n][-1] / lam[n] >= 4.0,
                f"{n}: the RANS chord is {R[0][n][-1] / lam[n]:.4f}x the planar"
                f" Poiseuille resistance of its own duct ({lam[n]:.6g} vs"
                f" {R[0][n][-1]:.6g}) - §5.25's laminar pair sat at 0.98x, so this is"
                " the friction regime talking, not the coupling")
    require(abs(R[0]["X1"][-1] / lam["X1"] / (R[0]["X2"][-1] / lam["X2"]) - 1.0) >= 0.2,
            f"and the two edges are not the same problem twice: {R[0]['X1'][-1] / lam['X1']:.3f}x"
            f" against {R[0]['X2'][-1] / lam['X2']:.3f}x, because their aspect ratios"
            " differ")

    # ---- 6. the series edge is not on the critical path ----
    allWorst = set()
    lone = 0
    for a, c in enumerate(arms):
        rLate = [res[a]["X2"][n] for n in range(1, N_STEPS)]
        require(max(abs(v - Q_IN) / Q_IN for v in q[a]["X2"]) <= 1e-15,
                f"{ARMS[a]}: the series edge's rate is decided by topology at every step,"
                f" not negotiated (worst {max(abs(v - Q_IN) / Q_IN for v in q[a]['X2']):.1e})")
        require(max(R[a]["X2"][1:]) / min(R[a]["X2"][1:]) - 1.0 <= 1e-6,
                f"{ARMS[a]}: and once its field is warm its coefficient stops moving too"
                f" - steps 2..{N_STEPS} span {max(R[a]['X2'][1:]) / min(R[a]['X2'][1:]) - 1:.1e}")
        require(max(rLate) <= TOL_INTERFACE / 100.0,
                f"{ARMS[a]}: the pinned edge agrees to {max(rLate):.2e} on every step but"
                f" the first - {TOL_INTERFACE / max(rLate):.0f}x inside what the interface"
                " asked for, because there is nothing left to iterate")
        require(min(sw[a]["X2"]) == max(sw[a]["X2"]),
                f"{ARMS[a]}: its block takes the same {sw[a]['X2'][0]} SIMPLE sweeps at"
                f" every time level, while the negotiated one takes"
                f" {min(sw[a]['X1'])}-{max(sw[a]['X1'])}")
        for n in sorted(c["pass"]):
            passes = c["pass"][n]
            for p in sorted(passes):
                rel = values(passes[p], "residuals")
                allWorst.add(int(value(passes[p], "worst")))
                if max(rel) > TOL_INTERFACE and min(rel) <= TOL_INTERFACE:
                    lone += 1
    require(allWorst == {0},
            f"on every pass of both arms the interface names edge 0 (X1) as its worst -"
            f" the schedule of a two-edge network with two RANS partners is the"
            f" negotiated edge's alone")
    require(lone >= 20,
            f"{lone} of the passes spent so far ran while the series edge was already"
            " under the tolerance: one number decides the whole network, so an edge that"
            " agrees early is refreshed with the one that does not")

    # ---- 7. the schedule, and what it cost ----
    nPassTot = []
    for a, c in enumerate(arms):
        nPass = sum(int(value(c["step"][n], "passes")) for n in c["step"])
        nPassTot.append(nPass)
        over = sum(1 for n in c["pass"] for p in c["pass"][n]
                   if max(values(c["pass"][n][p], "residuals")) > TOL_INTERFACE)
        require(over == nPass - N_STEPS,
                f"{ARMS[a]}: {over} of its {nPass} passes ended over tolerance and were"
                f" re-run, and exactly one per time level ({N_STEPS}) was accepted - the"
                " step stops when the WORST edge agrees")
        require(all(1 <= len(c["pass"][n]) <= MAX_PASSES for n in c["pass"]),
                f"{ARMS[a]}: every time level agreed inside the sub-iteration budget")
        require(max(len(c["pass"][n]) for n in c["pass"]) < MAX_PASSES,
                f"{ARMS[a]}: the budget was never reached - the worst step needed"
                f" {max(len(c['pass'][n]) for n in c['pass'])} of {MAX_PASSES} passes")
        require(c["solves"] == nPass * 2,
                f"{ARMS[a]}: {c['solves']} block solves for {nPass} passes x 2 edges, so"
                " the pinned region is refreshed by every pass rather than left on a"
                " half state")
    require(nPassTot[1] <= nPassTot[0],
            f"the accelerated arm needed {nPassTot[1]} passes against the plain"
            f" {nPassTot[0]} - on a network this easy the quasi-Newton buys"
            f" {nPassTot[0] - nPassTot[1]} pass{'es' if nPassTot[0] - nPassTot[1] != 1 else ''}"
            " out of a 53-pass run, and the first time level is where it shows")

    # ---- 8. the ramp of the negotiated branch: slopes, not chords ----
    ramp = []
    for a, c in enumerate(arms):
        x, Rv = q[a]["X1"], R[a]["X1"]
        mX = math.log(dpB[a]["X1"][-1] / dpB[a]["X1"][0]) \
            / math.log(x[-1] / x[0])
        mBy = 1.0 + math.log(pipeChord(*BY, qBy[a][N_STEPS]) / pipeChord(*BY, qBy[a][1])) \
                  / math.log(qBy[a][N_STEPS] / qBy[a][1])
        RX, RBy = Rv[-1], pipeChord(*BY, qBy[a][N_STEPS])
        Mtot = M[a]["X1"][-1] + M_BY
        slope = mX * RX + mBy * RBy
        chordR = RX + RBy
        xs, ys = x[1:-1], x[2:]
        k = float(len(xs))
        sx, sy = sum(xs), sum(ys)
        B = (k * math.fsum(u * v for u, v in zip(xs, ys)) - sx * sy) / (
            k * math.fsum(u * u for u in xs) - sx * sx)
        Acoef = (sy - B * sx) / k
        tauFit = -DT / math.log(B)
        qInfFit = Acoef / (1.0 - B)
        predict = lambda restore: DT / math.log(1.0 + DT * restore / Mtot)
        tauSlope, tauChord = predict(slope), predict(chordR)
        Bslope = 1.0 / (1.0 + DT * slope / Mtot)
        qInf = Q_IN * RBy / chordR
        q1 = Q_IN * (RBy + M_BY / DT) / (Mtot / DT + Rv[0] + RBy)
        gap = (x[-1] - qInfFit) / qInfFit
        tail = (x[0] - qInfFit) / qInfFit * B ** (N_STEPS - 1)
        dev = max(abs(x[n] - (qInf + (x[0] - qInf) * Bslope ** n)) / qInf
                  for n in range(N_STEPS))
        print(f"{ARMS[a]} ramp: tau {tauSlope:.6f} s from slopes against {tauFit:.6f} s"
              f" fitted ({(tauSlope / tauFit - 1) * 100:+.3f}%), chord formula"
              f" {tauChord:.6f} s ({(tauChord / tauFit - 1) * 100:+.2f}%); m_X {mX:.4f}"
              f" m_By {mBy:.4f},"
              f" weighted {slope / chordR:.4f}")
        require(0.5 < B < 0.95,
                f"{ARMS[a]}: the gap to the resistive share decays geometrically,"
                f" {1.0 / (1.0 - B):.1f} steps per decay")
        require(abs(tauSlope / tauFit - 1.0) <= 0.01,
                f"{ARMS[a]}: the ramp's time constant is the partners' SLOPES,"
                f" (M_X1+M_By)/(m_X R_X + m_By R_By) = {tauSlope:.4f} s against the"
                f" fitted {tauFit:.4f} s ({(tauSlope / tauFit - 1) * 100:+.2f}%) - the"
                " same closed form the laminar two-edge case satisfies with m = 1")
        require(abs(B / Bslope - 1.0) <= 0.005,
                f"{ARMS[a]}: per step that is the implicit Euler decay"
                f" 1/(1 + dt m R/M) = {Bslope:.6f} against the measured {B:.6f}"
                f" ({(B / Bslope - 1) * 100:+.3f}%)")
        require(tauChord / tauFit - 1.0 >= 0.5,
                f"{ARMS[a]}: the CHORD formula that the laminar pair satisfies to 1.6%"
                f" reads {tauChord:.4f} s here, {(tauChord / tauFit - 1) * 100:.0f}% high"
                " - the same miss §5.23 measured on a single turbulent edge, now with a"
                " second nonlinear partner in the loop that does not enter it")
        require(abs((tauChord / tauSlope) / (slope / chordR) - 1.0) <= 0.06,
                f"{ARMS[a]}: and the miss is the resistance-weighted exponent - chord/slope"
                f" {tauChord / tauSlope:.3f} against m = {slope / chordR:.3f}, the rest"
                " being the implicit step's own ln(1+x)")
        require(abs(x[0] / q1 - 1.0) <= 5e-3,
                f"{ARMS[a]}: the first step splits the sudden rate by inertance"
                f" ({x[0]:.6g} against {q1:.6g} m3/s, {(x[0] / q1 - 1) * 100:+.3f}%)")
        require(abs(x[0] / qInf - 1.0) >= 0.02,
                f"{ARMS[a]}: that inertial split puts the detailed branch"
                f" {abs(x[0] / qInf - 1) * 100:.1f}% away from the resistive share it"
                " relaxes to, not on it")
        require(abs(qInfFit / qInf - 1.0) <= 5e-3,
                f"{ARMS[a]}: the fitted destination {qInfFit:.6g} is the resistive split"
                f" of the two settled chords, {qInf:.6g} m3/s"
                f" ({(qInfFit / qInf - 1) * 100:+.3f}%)")
        require(0.0 < gap <= 4e-3 and abs(gap / tail - 1.0) <= 0.15,
                f"{ARMS[a]}: the series ends {gap * 100:.4f}% short of its own"
                f" asymptote, {tail * 100:.4f}% being the exponential tail its measured"
                " decay still has to burn at the 12-step cut")
        require(dev <= 4e-3,
                f"{ARMS[a]}: the whole 12-step trajectory follows the discrete inertial"
                f" solve built from the slopes, every step within {dev * 100:.3f}% of the"
                " operating point")
        ramp.append((tauFit, mX))

    # ---- 9. two operators, one answer ----
    for i, name in enumerate(NAME):
        for what, hist in (("rate", q), ("resistance", R), ("drop", dpB)):
            x, y = hist[0][name][-1], hist[1][name][-1]
            require(abs(x / y - 1.0) <= 5e-6,
                    f"{name} ends on one {what} whichever operator walked the interface:"
                    f" {x:.10g} and {y:.10g} ({abs(x / y - 1) * 100:.4f}%)")
    require(abs(ramp[0][0] / ramp[1][0] - 1.0) <= 5e-3,
            f"and the same ramp comes out of both: fitted time constant {ramp[0][0]:.6f}"
            f" s against {ramp[1][0]:.6f} s ({abs(ramp[0][0] / ramp[1][0] - 1) * 100:.3f}%),"
            f" partner exponent {ramp[0][1]:.4f} against {ramp[1][1]:.4f}")
    step1X2 = [res[a]["X2"][0] for a in range(len(arms))]
    require(step1X2[1] <= TOL_INTERFACE and step1X2[1] >= 5.0 * step1X2[0],
            f"on the first time level the accelerated arm accepted the pinned edge at"
            f" {step1X2[1]:.2e} against the plain one's {step1X2[0]:.2e} - the quasi-Newton"
            " mixes both edges' corrections, so the step where it moves the negotiated"
            f" edge hardest lets its neighbour in with the larger gap, still"
            f" {TOL_INTERFACE / step1X2[1]:.0f}x inside the tolerance")

    # ---- 10. each block's own field ----
    nu = NU
    for a, c in enumerate(arms):
        for i, name in enumerate(NAME):
            fin = c["block"][i]
            L, H, W = EDGES[name]
            area = H * W
            dh = 2.0 * H
            f, plug = value(fin, "fluxIn"), value(fin, "plugSpeed")
            uTau, yp = value(fin, "wallUTau"), value(fin, "wallYPlusMax")
            slope, dp = abs(value(fin, "probeSlope")), value(fin, "dpPa")
            require(abs(f / q[a][name][-1] - 1.0) <= 1e-9,
                  f"{ARMS[a]} {name}: the block carries the rate its edge ended with"
                  f" ({f:.6g} m3/s)")
            require(abs(value(fin, "inletArea") - area) <= 1e-12 * area,
                  f"{ARMS[a]} {name}: the inlet patch is the section the inertance was"
                  f" referred to ({value(fin, 'inletArea'):.6g} m2)")
            require(abs(value(fin, "fluxOut") - f) <= 1e-11 * f,
                  f"{ARMS[a]} {name}: what enters the region leaves it")
            require(value(fin, "planeFluxMax") - value(fin, "planeFluxMin")
                    <= 1e-11 * f,
                  f"{ARMS[a]} {name}: every cross-section plane inside carries the same"
                  " rate")
            require(abs(value(fin, "inletArea") * plug - f) <= 1e-12 * f,
                  f"{ARMS[a]} {name}: plug speed times inlet area is the carried rate")
            require(RHO * plug * dh / MU > 4000.0,
                  f"{ARMS[a]} {name}: the region is turbulent at Re_bulk"
                  f" {RHO * plug * dh / MU:.0f}, which is what the partner model is")
            require(value(fin, "nutMax") / nu > 20.0,
                  f"{ARMS[a]} {name}: eddy viscosity {value(fin, 'nutMax') / nu:.1f}x the"
                  " molecular one - the RANS term is doing the work, not the laminar one")
            require(yp >= YPLUS_LOG_LAW and yp < 300.0,
                  f"{ARMS[a]} {name}: the wall cells sit at y+ {yp:.1f}, inside the"
                  " log-law band the wall functions assume")
            closure = slope * H / 2.0
            require(abs(closure - RHO * uTau ** 2)
                    <= (0.01 if name == "X1" else 0.08) * closure,
                    f"{ARMS[a]} {name}: the developed gradient balances the wall shear"
                    f" the model reports - |dp/dx|*H/2 {closure:.6g} against rho*u_tau^2"
                    f" {RHO * uTau ** 2:.6g} Pa ({closure / (RHO * uTau ** 2) - 1:+.2%})"
                    + (" on a 0.2 m duct of 20 x-cells" if name == "X1" else
                       " on a 0.12 m duct of 12 x-cells, half of which is still entrance"))
            require(0.94 <= dp / L / slope <= 1.06,
                  f"{ARMS[a]} {name}: the whole drop is {dp / L / slope:.4f} of the"
                  " developed line's - the entry region of a duct that starts from a plug"
                  " profile, not a mismatch")
            require(value(fin, "probeSpread") <= 0.02 * dp,
                  f"{ARMS[a]} {name}: section means hug the developed line to"
                  f" {value(fin, 'probeSpread') / dp * 100:.2f}% of the drop over"
                  f" {value(fin, 'probeSections'):.0f} sections")
            require(abs(value(fin, "uxMax") / plug - 1.0) <= 0.15,
                  f"{ARMS[a]} {name}: a turbulent profile, u_max/u_bulk"
                  f" {value(fin, 'uxMax') / plug:.4f} against 1.5 for the laminar"
                  " parabola these same two ducts gave §5.25")
            uPope = popeBulkSpeed(uTau, nu, dh)
            require(abs(plug / uPope - 1.0) < 0.08,
                  f"{ARMS[a]} {name}: Pope's smooth-wall law put at the block's own u_tau"
                  f" gives the bulk speed it carries ({plug:.4f} against {uPope:.4f} m/s,"
                  f" {(plug / uPope - 1) * 100:+.1f}%)")
            require(abs(dpB[a][name][-1] - dp) <= 1e-9 * dp,
                  f"{ARMS[a]} {name}: the drop the interface ended with is the drop the"
                  " block's own field measures")

    tot = sum(c["solves"] for c in arms)
    print(f"cost: {tot} block solves and {sum(nPassTot)} interface passes for"
          f" {N_STEPS * len(arms)} agreed time levels on two RANS edges")

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: two RANS detailed regions agree with one network inside every time"
          " step, on a schedule the negotiated edge alone owns")


if __name__ == "__main__":
    sys.exit(main())
