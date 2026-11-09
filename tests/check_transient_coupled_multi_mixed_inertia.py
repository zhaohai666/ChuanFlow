#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient-multi over two RANS detailed regions whose
inertance is sourced PER EDGE in one of the three arms, and verify that each
edge answers to its own source while its neighbour's choice is free to disagree.

transientCoupledMultiTurbulentInertia could only state the attribution per
whole interface: one arm took both regions' inertance from the proxy bore, the
next from both blocks' meshes. A network whose regions are described by drawings
of uneven quality needs the choice per edge, so `couple.inertiaFrom` now carries
one source per coupled edge for the arms that declare it. Here both bores are
mis-stated - X1 by 1.273240, X2 by 1.696946 - and both partners are nonlinear, so
the mixed arm (X1 from its block, X2 from its bore) is the first run in this
repository that installs two different sources at the same time.

What that buys, and what this checker refuses to approximate:

  1. three arms, one table  - the runs differ only in `inertias` / `inertiaFrom`,
                              same seeds, both partners k-omega SST, none accelerated
  2. what got installed     - every edge carries the inertance ITS arm named for
                              ITSE, every step, and the block's own measurement is
                              the same pair of numbers in all three arms
  3. superposition, not average - the mixed arm's parallel edge is bit-identical to
                              the region arm's at every step (rate, coefficient,
                              both drops, residual) while its series edge carries
                              the bore arm's number, and the bore arm differs from
                              both on the parallel edge
  4. the interface          - dpNetwork == R*Q exactly in all three arms, residual
                              == the gap between the two solvers, under tolerance
  5. the discrete 1D rows   - each edge's node drop closes with that EDGE's own
                              inertance, in the arm that mixed them; the bypass and
                              leads close with published Churchill
  6. the starting kick      - the series edge's first step pays M*Q/dt out of the
                              inertance its arm installed: the mixed arm pays the
                              bore's 15272.5 Pa, not its neighbour's 9000 Pa
  7. the ramp               - tau is (M_X1 + M_By) over the sum of SLOPES with THIS
                              edge's M, and the mixed arm's ramp is the region
                              arm's ramp because its X1 is the region arm's X1
  8. the destination        - the resistive split holds no M, so all three arms
                              relax to the same rate while their paths differ
  9. the partner's law      - the same exponent read from all three arms: which
                              inertance the neighbour was given cannot touch it
 10. schedule, cost, fields - who pays for the interface, and each block's own
                              wall-resolved field

usage: check_transient_coupled_multi_mixed_inertia.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# fluid and geometry of tests/cases/transientCoupledMultiMixedInertia/case.dict
MU, RHO = 1.0e-3, 1000.0
Q_IN = 1.2e-4
P_REF = 1.0e6
DT, N_STEPS = 0.05, 12
MAX_PASSES = 6
TOL = 2e-4
NAME = ("X1", "X2")
ARMS = ("bore", "region", "mixed")
DECLARED = {"bore": ("network", "network"),
            "region": ("partner", "partner"),
            "mixed": ("partner", "network")}
EDGES = {"X1": (0.2, 0.005, 0.005), "X2": (0.12, 0.004, 0.008)}
BORE = {"X1": 0.005, "X2": 0.0049}
BY = (0.3, 0.008)
LEAD = (0.3, 0.02)
YPLUS_LOG_LAW = 11.81
M_BY = RHO * BY[0] / (0.25 * math.pi * BY[1] ** 2)
M_LEAD = RHO * LEAD[0] / (0.25 * math.pi * LEAD[1] ** 2)


def value(line, key):
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([-+0-9.eE]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return float(m.group(1))


def word(line, key):
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+(\S+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return m.group(1)


def churchill(re, roughness=0.0):
    """Churchill's 1977 Darcy friction factor - the published correlation the 1D
    side is closed with, so the network is judged from outside this repository."""
    re = max(re, 1e-12)
    a = (2.457 * math.log(1.0 / ((7.0 / re) ** 0.9 + 0.27 * roughness))) ** 16
    b = (37530.0 / re) ** 16
    return 8.0 * ((8.0 / re) ** 12 + 1.0 / (a + b) ** 1.5) ** (1.0 / 12.0)


def pipeChord(length, diameter, q):
    area = 0.25 * math.pi * diameter ** 2
    v = q / area
    return (churchill(RHO * v * diameter / MU) * length / diameter
            * 0.5 * RHO * v * v / q)


def ductChord(length, height, width, q):
    return 12.0 * MU * length / (height ** 3 * width)


def popeBulkSpeed(uTau, dh):
    u = uTau * 20.0
    for _ in range(200):
        f = 8.0 * (uTau / u) ** 2
        u = uTau * math.sqrt(8.0) * (2.0 * math.log10(u * dh / NU * math.sqrt(f)) - 0.6)
    return u


def parse(lines):
    """One dict per arm: the per-edge source it was given, its accepted states,
    its passes, and its block reports."""
    arms = []
    cur = None
    for ln in lines:
        if ln.startswith("transient start ") and ln.count(":") == 1 and " edge " not in ln:
            t = ln.split()
            cur = dict(name=t[2].rstrip(":"), inertia=t[4], accelerate=t[6],
                       source={}, steps={}, passes=[], blocks=[], summary=None,
                       solves=0)
            arms.append(cur)
        elif cur is None:
            continue
        elif ln.startswith("SIMPLE iter 1:"):
            cur["solves"] += 1
        elif ln.startswith("transient arm ") and " inertia " in ln:
            t = ln.split()
            cur["source"][t[4]] = word(ln, "inertia")
        elif ln.startswith("transient step ") and " pass " in ln:
            cur["passes"].append(ln)
        elif ln.startswith("transient step ") and " edge " in ln:
            t = ln.split()
            cur["steps"].setdefault(int(t[2]), {})[int(t[4])] = ln
        elif ln.startswith("transient step ") and " newton " in ln:
            cur["steps"].setdefault(int(ln.split()[2].rstrip(":")), {})["s"] = ln
        elif ln.startswith("nodes step "):
            cur["steps"][int(ln.split()[2].rstrip(":"))]["nd"] = ln
        elif ln.startswith("pipes step "):
            cur["steps"][int(ln.split()[2].rstrip(":"))]["pp"] = ln
        elif ln.startswith("transient block "):
            cur["blocks"].append(ln)
        elif ln.startswith("transient coupled start ") and " edge " not in ln:
            cur["summary"] = ln
    return arms


def main():
    binary, case = sys.argv[1], sys.argv[2]
    r = subprocess.run([binary, os.path.basename(case)],
                       capture_output=True, text=True,
                       cwd=os.path.dirname(os.path.abspath(case)))
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("transient multi mixed-inertia solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    arms = parse(lines)

    # ---- 1. three arms, one table ----
    require([a["name"] for a in arms] == list(ARMS),
            f"all three arms ran, in the order declared (got {[a['name'] for a in arms]})")
    for a in arms:
        require(sorted(a["steps"]) == list(range(1, N_STEPS + 1)),
                f"{a['name']}: all {N_STEPS} time levels reported an accepted state")
        require(sorted(a["source"]) == list(NAME),
                f"{a['name']}: the run reported a source for each coupled edge")
        require(tuple(a["source"][n] for n in NAME) == DECLARED[a["name"]],
                f"{a['name']}: {a['source']['X1']} on X1, {a['source']['X2']} on X2"
                f" - the case says {DECLARED[a['name']]}")
        require(a["accelerate"] == "no",
                f"{a['name']}: plain block Gauss-Seidel, so an arm-to-arm difference"
                " can only be the inertance")
        require(all(word(b, "turb") == "kOmegaSST" for b in a["blocks"]),
                f"{a['name']}: both partners faced the interface as RANS blocks")
    require(len({tuple(a["source"][n] for n in NAME) for a in arms}) == 3,
            "the three arms state three different per-edge tables - the mixed one is"
            " not a re-run of either whole-interface choice")

    # ---- 2. what each edge was given ----
    M = [{i: [value(a["steps"][n][i], "M") for n in range(1, N_STEPS + 1)]
          for i in (0, 1)} for a in arms]
    inertSelf = [[value(b, "inertance") for b in a["blocks"]] for a in arms]
    for a_i, a in enumerate(arms):
        for i, n in enumerate(NAME):
            L, H, W = EDGES[n]
            exact = RHO * L / (H * W)
            bore = RHO * L / (0.25 * math.pi * BORE[n] ** 2)
            want = bore if a["source"][n] == "network" else inertSelf[a_i][i]
            require(max(abs(x - want) for x in M[a_i][i]) <= 1e-9 * abs(want),
                    f"{a['name']} {n}: its momentum row carries"
                    f" {'the bore rho*L/A = %.10g' % bore if a['source'][n] == 'network' else 'the block-measured rho*V/A^2 = %.10g' % exact}"
                    " kg/m4, every step")
            require(abs(inertSelf[a_i][i] - exact) <= 1e-9 * exact,
                    f"{a['name']} {n}: its block measures {inertSelf[a_i][i]:.10g} kg/m4"
                    f" on its own mesh, equal to rho*L/(H*W) = {exact:.10g}")
    for i, n in enumerate(NAME):
        L, H, W = EDGES[n]
        ratio = (RHO * L / (0.25 * math.pi * BORE[n] ** 2)) / (RHO * L / (H * W))
        require(abs(M[0][i][0] / M[1][i][0] - ratio) <= 1e-9,
                f"{n}: the bore overstates this region by {ratio:.6f}"
                f" (H/W {W / H:.1f}), so the two sources name different numbers here")
        mixedOnBore = DECLARED["mixed"][i] == "network"
        say = (f"the same bore as the bore arm, so the two differ by nothing here"
               if mixedOnBore else
               f"the mesh, so the bore arm sits {ratio:.6f}x above it")
        require(abs(M[0][i][0] / M[2][i][0] - (1.0 if mixedOnBore else ratio)) <= 1e-9,
                f"{n}: the mixed arm names {say} ({M[0][i][0]:.10g} against"
                f" {M[2][i][0]:.10g} kg/m4)")
        require({inertSelf[k][i] for k in range(len(arms))} == {inertSelf[0][i]},
                f"{n}: the region's own inertance is one number whatever the arms"
                f" decide to do with it ({inertSelf[0][i]:.10g})")

    # ---- 3. superposition, not average ----
    reg, mix, bore = 1, 2, 0
    for i, n in enumerate(NAME):
        same = all(arms[reg]["steps"][s][i] == arms[mix]["steps"][s][i]
                   for s in range(1, N_STEPS + 1))
        wantSame = DECLARED["mixed"][i] == DECLARED["region"][i]
        require(same == wantSame,
                f"{n}: the mixed arm's accepted states are"
                f" {'BIT-IDENTICAL to the region arm at all ' + str(N_STEPS) + ' steps' if wantSame else 'NOT those of the region arm - it was given a different source'}"
                f" ({'same edge, same source, so the same numbers' if wantSame else 'a different inertance in the row'})")
    require(M[mix][0] == M[reg][0] and M[mix][1] == M[bore][1],
            f"the mixed arm installs the region's number on X1 ({M[mix][0][0]:.10g}) and"
            f" the bore's on X2 ({M[mix][1][0]:.10g}) in the SAME run - two sources,"
            " one momentum system")
    require(M[mix][1] != M[reg][1],
            f"and the two really are different numbers in one run:"
            f" X2 reads {M[mix][1][0]:.10g} against the region arm's {M[reg][1][0]:.10g}")
    qDiff = abs(value(arms[bore]["steps"][1][0], "q") - value(arms[reg]["steps"][1][0], "q"))
    require(qDiff > 5e-6,
            f"the whole-interface choices are not the mixed one either: X1's first step"
            f" is {value(arms[bore]['steps'][1][0], 'q'):.6g} m3/s under the bore against"
            f" {value(arms[reg]['steps'][1][0], 'q'):.6g} under the mesh"
            f" ({qDiff / Q_IN * 100:.2f}% of the injected rate apart)")

    # ---- 4. the interface, and 5. the rows it solved ----
    series = []
    for a_i, a in enumerate(arms):
        st = a["steps"]
        q = {i: [value(st[n][i], "q") for n in range(1, N_STEPS + 1)] for i in (0, 1)}
        R = {i: [value(st[n][i], "R") for n in range(1, N_STEPS + 1)] for i in (0, 1)}
        dp = {i: [value(st[n][i], "dpBlock") for n in range(1, N_STEPS + 1)]
              for i in (0, 1)}
        res = {i: [value(st[n][i], "residual") for n in range(1, N_STEPS + 1)]
               for i in (0, 1)}
        pas = [int(value(st[n]["s"], "passes")) for n in range(1, N_STEPS + 1)]
        nd = {k: [value(st[n]["nd"], k) for n in range(1, N_STEPS + 1)]
              for k in ("Src", "J1", "J2", "J3", "Snk")}
        pp = {k: [value(st[n]["pp"], k) for n in range(1, N_STEPS + 1)]
              for k in ("A", "X1", "By", "X2", "B")}

        w_id = max(abs(value(st[n][i], "dpNetwork") - R[i][n - 1] * q[i][n - 1])
                   for i in (0, 1) for n in range(1, N_STEPS + 1))
        require(w_id == 0.0,
                f"{a['name']}: the network still charges R times the rate, under a"
                f" mixed source table and a mis-stated bore both (worst {w_id:.1e} Pa)")
        w_gap = max(abs(res[i][n - 1] - abs(value(st[n][i], "dpNetwork") - dp[i][n - 1])
                        / dp[i][n - 1]) for i in (0, 1) for n in range(1, N_STEPS + 1))
        require(w_gap == 0.0,
                f"{a['name']}: the residual is the gap between the two solvers"
                f" (worst {w_gap:.1e})")
        require(max(max(res[i]) for i in (0, 1)) <= TOL,
                f"{a['name']}: every accepted level is under the case tolerance"
                f" (worst {max(max(res[i]) for i in (0, 1)):.3e})")
        require(max(value(st[n]["s"], "newton") for n in range(1, N_STEPS + 1)) < 1e-11,
                f"{a['name']}: the 1D step solved its own system at every level")
        require(all(abs(q[i][n - 1] - pp[NAME[i]][n - 1]) <= 1e-15 * Q_IN
                    for i in (0, 1) for n in range(1, N_STEPS + 1))
                and all(abs(pp["A"][n - 1] - Q_IN) <= 1e-15 * Q_IN
                        and abs(pp["A"][n - 1] - pp["X1"][n - 1] - pp["By"][n - 1])
                        <= 1e-12 * Q_IN for n in range(1, N_STEPS + 1)),
                f"{a['name']}: the negotiated rates are the pipes' own and the branch"
                " split adds up to the lead")
        require(all(abs(nd["Snk"][n - 1] - P_REF) <= 1e-9 for n in range(1, N_STEPS + 1)),
                f"{a['name']}: the sink holds its reference pressure")

        dp12 = [nd["J1"][n - 1] - nd["J2"][n - 1] for n in range(1, N_STEPS + 1)]
        dp23 = [nd["J2"][n - 1] - nd["J3"][n - 1] for n in range(1, N_STEPS + 1)]

        def row(drop, flow, i):
            return max(abs(drop[n] - (R[i][n] * flow[n]
                                      + M[a_i][i][n] * (flow[n]
                                                        - (flow[n - 1] if n else 0.0))
                                        / DT)) / abs(drop[n])
                       for n in range(N_STEPS))

        rX1, rX2 = row(dp12, q[0], 0), row(dp23, q[1], 1)
        require(rX1 <= 1e-9 and rX2 <= 1e-12,
                f"{a['name']}: each region's node drop closes with THAT EDGE's own"
                f" inertance ({rX1:.1e} parallel, {rX2:.1e} series) - X1 from"
                f" {a['source']['X1']}, X2 from {a['source']['X2']}")
        rBy = max(abs(dp12[n] - (pipeChord(*BY, pp["By"][n]) * pp["By"][n]
                                 + M_BY * (pp["By"][n]
                                           - (pp["By"][n - 1] if n else 0.0)) / DT))
                  / dp12[n] for n in range(N_STEPS))
        rLead = max(min([
            abs((nd["Src"][n] - nd["J1"][n])
                - (pipeChord(*LEAD, pp["A"][n]) * pp["A"][n]
                   + M_LEAD * (pp["A"][n] - (pp["A"][n - 1] if n else 0.0)) / DT))
            / abs(nd["Src"][n] - nd["J1"][n]),
            abs((nd["J3"][n] - nd["Snk"][n])
                - (pipeChord(*LEAD, pp["B"][n]) * pp["B"][n]
                   + M_LEAD * (pp["B"][n] - (pp["B"][n - 1] if n else 0.0)) / DT))
            / abs(nd["J3"][n] - nd["Snk"][n])]
            for n in range(N_STEPS)))
        require(rBy <= 1e-12 and rLead <= 1e-12,
                f"{a['name']}: the bypass and the leads still close against published"
                f" Churchill chords plus their own inertance ({rBy:.1e}, {rLead:.1e})")

        series.append(dict(q=q, R=R, dp=dp, res=res, pas=pas, pp=pp, nd=nd,
                           dp12=dp12, dp23=dp23))

    # ---- 6. the starting kick follows the edge, not the arm ----
    kick = [M[i][1][0] * Q_IN / DT for i in range(len(arms))]
    for i, a in enumerate(series):
        err = abs(a["dp23"][0] - (a["R"][1][0] * Q_IN + kick[i])) / a["dp23"][0]
        require(err <= 1e-12,
                f"{arms[i]['name']}: the series edge's first node drop is friction at"
                f" Q_IN plus its OWN inertance's kick ({a['dp23'][0]:.6f} ="
                f" {a['R'][1][0] * Q_IN:.4f} + {kick[i]:.4f}, {err:.1e})")
    require(arms[mix]["source"]["X2"] == "network"
            and abs(series[mix]["dp23"][0] - series[bore]["dp23"][0])
            <= 1e-3 * abs(series[bore]["dp23"][0]),
            f"the mixed arm pays the BORE's kick on X2 ({series[mix]['dp23'][0]:.6f} Pa"
            f" against the bore arm's {series[bore]['dp23'][0]:.6f} Pa),"
            f" not the mesh arm's {series[reg]['dp23'][0]:.6f} Pa - the difference it"
            f" keeps is {series[mix]['dp23'][0] - series[reg]['dp23'][0]:.6f} Pa against"
            f" dM*Q/dt = {kick[bore] - kick[reg]:.6f} Pa")
    require(series[mix]["q"] == series[reg]["q"]
            and abs(series[mix]["dp23"][0] - series[reg]["dp23"][0]) > 1e3,
            f"and the mixed arm is therefore neither of the other two: its X1 is the"
            f" region arm's X1 to the last bit while its X2 drop is"
            f" {series[mix]['dp23'][0] - series[reg]['dp23'][0]:.6f} Pa above the region"
            " arm's")

    # ---- 7. the ramp, per arm ----
    for a_i, a in enumerate(arms):
        s = series[a_i]
        q, R, pp = s["q"], s["R"], s["pp"]
        xs, ys = q[0][1:-1], q[0][2:]
        k = float(len(xs))
        sx, sy = sum(xs), sum(ys)
        Bc = (k * math.fsum(u * w for u, w in zip(xs, ys)) - sx * sy) / (
            k * math.fsum(u * u for u in xs) - sx * sx)
        qInf = ((sy - Bc * sx) / k) / (1.0 - Bc)
        tauFit = -DT / math.log(Bc)
        s["B"], s["tau"], s["qInf"] = Bc, tauFit, qInf
        Mtot = M[a_i][0][0] + M_BY
        RX, RBy = R[0][-1], pipeChord(*BY, pp["By"][-1])
        mX = math.log(s["dp"][0][-1] / s["dp"][0][0]) / math.log(q[0][-1] / q[0][0])
        mBy = (1.0 + math.log(pipeChord(*BY, pp["By"][-1])
                              / pipeChord(*BY, pp["By"][0]))
               / math.log(pp["By"][-1] / pp["By"][0]))
        s["mX"] = mX
        slope = mX * RX + mBy * RBy
        tauS = DT / math.log(1.0 + DT * slope / Mtot)
        require(abs(tauS / tauFit - 1.0) <= 0.01,
                f"{a['name']}: the ramp's time constant is (M_X1+M_By) over the partners'"
                f" SLOPES = {tauS:.6f} s against the fitted {tauFit:.6f} s"
                f" ({(tauS / tauFit - 1) * 100:+.3f}%), with X1 on {a['source']['X1']}"
                f" ({M[a_i][0][0]:.6g} kg/m4) and X2 on {a['source']['X2']}")
        require(abs(qInf / (Q_IN * RBy / (RBy + RX)) - 1.0) <= 0.02,
                f"{a['name']}: and it relaxes towards the resistive split of the settled"
                f" chords, {qInf:.6g} m3/s")
    require(series[mix]["tau"] == series[reg]["tau"],
            f"the mixed arm's ramp IS the region arm's ramp - tau {series[mix]['tau']:.6f}"
            f" s against {series[reg]['tau']:.6f} s, bit for bit, because the edge that"
            " does the ramping was given the same inertance")
    tauBore = series[bore]["tau"]
    mRatio = (M[bore][0][0] + M_BY) / (M[reg][0][0] + M_BY)
    require(abs(tauBore / series[reg]["tau"] / mRatio - 1.0) <= 0.02,
            f"while the bore arm's is the same ratio of inertances away:"
            f" {tauBore:.6f} s / {series[reg]['tau']:.6f} s = {tauBore / series[reg]['tau']:.5f}"
            f" against {mRatio:.5f}")

    # ---- 8. the destination holds no M ----
    spread = max(s["qInf"] for s in series) - min(s["qInf"] for s in series)
    require(spread / min(s["qInf"] for s in series) <= 2e-4,
            f"all three arms relax to one rate: fitted asymptotes"
            f" {' / '.join('%.6g' % s['qInf'] for s in series)} m3/s, spread"
            f" {spread / min(s['qInf'] for s in series) * 100:.4f}% - while their end"
            f" rates still differ by"
            f" {(max(s['q'][0][-1] for s in series) - min(s['q'][0][-1] for s in series)) / Q_IN * 100:.3f}%"
            " of the injected rate")

    # ---- 9. the partner's law is not the neighbour's business ----
    ex = [s["mX"] for s in series]
    require(max(ex) - min(ex) <= 0.01,
            f"X1's block answers with one law under every table:"
            f" dp = c*q^{ex[0]:.4f} / {ex[1]:.4f} / {ex[2]:.4f} (bore / region / mixed)"
            " - attribution does not touch the slope")
    require(all(max(abs(x - Q_IN) for x in s["q"][1]) / Q_IN <= 1e-14 for s in series),
            "the series edge's rate is pinned in all three arms, which is why its"
            " nonlinearity has no exponent to read and only its inertance is per-arm")

    # ---- 10. schedule, cost, and each block's own field ----
    for a_i, a in enumerate(arms):
        s = series[a_i]
        idx = [int(value(p, "worst")) for p in a["passes"]]
        require(set(idx) == {0},
                f"{a['name']}: the schedule stops on the negotiated edge in all"
                f" {len(idx)} passes")
        require(sum(s["pas"]) == len(a["passes"]),
                f"{a['name']}: one pass is one network solve - {sum(s['pas'])} accepted"
                f" of {len(a['passes'])} run, {len(a['passes']) - N_STEPS} over tolerance")
        require(all(1 <= p <= MAX_PASSES for p in s["pas"]),
                f"{a['name']}: every level agrees inside the allowed sub-iterations"
                f" (max {max(s['pas'])} of {MAX_PASSES})")
    nPass = [len(a["passes"]) for a in arms]
    require(all(a["solves"] == 2 * len(a["passes"]) for a in arms),
            f"cost: {' / '.join('%d passes' % p for p in nPass)} interface passes for"
            f" {3 * N_STEPS} agreed time levels and {sum(a['solves'] for a in arms)}"
            " block solves - every pass runs every partner once, under whichever"
            " source that edge was given")
    for a_i, a in enumerate(arms):
        s = series[a_i]
        for i, n in enumerate(NAME):
            fin = a["blocks"][i]
            L, H, W = EDGES[n]
            f = value(fin, "fluxIn")
            plug = value(fin, "plugSpeed")
            require(abs(f / s["q"][i][-1] - 1.0) <= 1e-4,
                    f"{a['name']} {n}: the block carries the rate its edge ended with"
                    f" ({f:.6g} m3/s)")
            require(abs(value(fin, "inletArea") - H * W) <= 1e-12 * H * W,
                    f"{a['name']} {n}: the inlet patch is the section its inertance was"
                    " referred to")
            require(abs(value(fin, "fluxOut") - f) <= 1e-9 * f,
                    f"{a['name']} {n}: what enters the region leaves it")
            require(RHO * plug * (2.0 * H) / MU > 4000.0,
                    f"{a['name']} {n}: turbulent at Re_bulk"
                    f" {RHO * plug * 2.0 * H / MU:.0f}")
            yp = value(fin, "wallYPlusMax")
            require(YPLUS_LOG_LAW <= yp < 300.0,
                    f"{a['name']} {n}: wall cells at y+ {yp:.1f}, inside the log-law band")
            closure = abs(value(fin, "probeSlope")) * H / 2.0
            require(abs(closure - RHO * value(fin, "wallUTau") ** 2) <= 0.08 * closure,
                    f"{a['name']} {n}: developed gradient balances wall shear to"
                    f" {(closure / (RHO * value(fin, 'wallUTau') ** 2) - 1) * 100:+.2f}%")
            require(abs(value(fin, "dpPa") - s["dp"][i][-1]) <= 1e-6 * s["dp"][i][-1],
                    f"{a['name']} {n}: the drop the interface ended with is the drop the"
                    " block's own field measures")
    print(f"note: X1 reads {' / '.join(arms[a_i]['source']['X1'] for a_i in range(3))}"
          f" across the arms, X2 reads"
          f" {' / '.join(arms[a_i]['source']['X2'] for a_i in range(3))}")

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: one run installs two different inertance sources, and each edge"
          " answers to its own")


if __name__ == "__main__":
    sys.exit(main())
