#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient-multi over two RANS detailed regions as a
complete 2x2: which source each region's inertia comes from, and which operator
negotiates the interface.

Every earlier run in this chain varied one of those two axes. Attribution got
two arms and a plain loop (transientCoupledMultiTurbulentInertia); the interface
operator got two arms with the attribution degree of freedom deliberately
neutralised (transientCoupledMultiTurbulent1d3d). What was never on the table is
the crossing, and the crossing is where the two accounts could have leaked into
each other: quasi-Newton changes the UPDATE the interface makes, and the update
is the only thing between a block's reported drop and the momentum row that
carries the region's mass.

So the claim this checker tests is that the axes are independent, in both
directions and with no approximations:

  1. four arms, two axes    - the runs differ only in `inertias` and
                              `accelerates`; both partners are k-omega SST; all
                              four seeds identical; the (source, operator) pairs
                              are the complete factorial, not four near-misses
  2. what got installed     - each edge carries the inertance its source names,
                              every step, and the OPERATOR does not enter that
                              table: within a source pair the two arms install
                              the same numbers to the last bit
  3. the interface          - dpNetwork == R*Q exactly, in all four arms: the
                              chord is a property of the exchange, not of the
                              loop that iterates it
  4. the discrete 1D rows   - every arm's node drop closes with its own edge's
                              inertance, so an accepted quasi-Newton state is
                              still a state of the same momentum system
  5. the operator axis      - within a source pair the two arms end on the same
                              rate, the same asymptote, the same tau and the same
                              first-step node drop; the operator buys passes, not
                              physics
  6. the source axis        - tau still is (M_X1 + M_By) over the sum of SLOPES,
                              the tau ratio still equals the inertance ratio, and
                              the destination still holds no M - four arms, one
                              destination
  7. the starting kick      - the pinned series edge pays M*Q/dt out of the
                              source the arm named, and pays the same kick under
                              either operator
  8. the residual floors    - quasi-Newton acts on the interface as a whole, so
                              the edge it never needed to move rides along with a
                              residual an order of magnitude above the plain
                              loop's; both stay inside the case tolerance
  9. schedule, cost         - who pays for the interface, and that one pass is
                              one network solve under either operator
 10. the blocks' own fields - wall-resolved closure, quoted from outside the
                              repository, in every arm

usage: check_transient_coupled_multi_turbulent_inertia_newton.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# fluid and geometry of tests/cases/transientCoupledMultiTurbulentInertiaNewton/case.dict
MU, RHO = 1.0e-3, 1000.0
Q_IN = 1.2e-4
P_REF = 1.0e6
DT, N_STEPS = 0.05, 12
MAX_PASSES = 6
TOL = 2e-4
NAME = ("X1", "X2")
ARMS = ("bore_plain", "bore_newton", "region_plain", "region_newton")
SOURCE = {"bore_plain": "network", "bore_newton": "network",
          "region_plain": "partner", "region_newton": "partner"}
OPERATOR = {"bore_plain": "no", "bore_newton": "yes",
            "region_plain": "no", "region_newton": "yes"}
PAIRS = (("bore_plain", "bore_newton"), ("region_plain", "region_newton"))
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


def parse(lines):
    """One dict per arm: the source and operator it was given, its accepted
    states, its passes, and its block reports."""
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
        raise SystemExit("transient multi inertia-x-Newton solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    arms = parse(lines)
    at = {a["name"]: k for k, a in enumerate(arms)}

    # ---- 1. four arms, two axes ----
    require([a["name"] for a in arms] == list(ARMS),
            f"all four arms ran, in the order declared (got {[a['name'] for a in arms]})")
    for a in arms:
        require(sorted(a["steps"]) == list(range(1, N_STEPS + 1)),
                f"{a['name']}: all {N_STEPS} time levels reported an accepted state")
        require(tuple(a["source"][n] for n in NAME) == (SOURCE[a["name"]],) * 2,
                f"{a['name']}: {a['source']['X1']} on X1, {a['source']['X2']} on X2"
                f" - the case says {SOURCE[a['name']]} on both")
        require(a["accelerate"] == OPERATOR[a["name"]],
                f"{a['name']}: interface operator is"
                f" {'quasi-Newton' if OPERATOR[a['name']] == 'yes' else 'plain block Gauss-Seidel'},"
                " as its name says")
        require(a["inertia"] == SOURCE[a["name"]],
                f"{a['name']}: the arm header and the per-edge table agree on the source")
        require(all(word(b, "turb") == "kOmegaSST" for b in a["blocks"]),
                f"{a['name']}: both partners faced the interface as RANS blocks")
    combos = {(SOURCE[a["name"]], OPERATOR[a["name"]]) for a in arms}
    require(combos == {("network", "no"), ("network", "yes"),
                       ("partner", "no"), ("partner", "yes")},
            "the four arms are the complete 2x2 - two sources times two operators,"
            " so every axis has a partner to be compared with")

    # ---- 2. what got installed, and that the operator is not part of it ----
    M = [{i: [value(a["steps"][n][i], "M") for n in range(1, N_STEPS + 1)]
          for i in (0, 1)} for a in arms]
    inertSelf = [[value(b, "inertance") for b in a["blocks"]] for a in arms]
    for a_i, a in enumerate(arms):
        for i, n in enumerate(NAME):
            L, H, W = EDGES[n]
            exact = RHO * L / (H * W)
            bore = RHO * L / (0.25 * math.pi * BORE[n] ** 2)
            want = bore if SOURCE[a["name"]] == "network" else inertSelf[a_i][i]
            say = ("the bore rho*L/A" if SOURCE[a["name"]] == "network"
                   else "the block-measured rho*V/A^2")
            require(max(abs(x - want) for x in M[a_i][i]) <= 1e-9 * abs(want),
                    f"{a['name']} {n}: its momentum row carries {say} ="
                    f" {want:.10g} kg/m4, every step")
            require(abs(inertSelf[a_i][i] - exact) <= 1e-9 * exact,
                    f"{a['name']} {n}: its block measures {inertSelf[a_i][i]:.10g} kg/m4"
                    f" on its own mesh, equal to rho*L/(H*W) = {exact:.10g}")
    for i, n in enumerate(NAME):
        L, H, W = EDGES[n]
        ratio = (RHO * L / (0.25 * math.pi * BORE[n] ** 2)) / (RHO * L / (H * W))
        require(abs(M[0][i][0] / M[2][i][0] - ratio) <= 1e-9,
                f"{n}: the bore overstates this region by {ratio:.6f}"
                f" (H/W {W / H:.1f}), so the two sources name different numbers here")
        require({inertSelf[k][i] for k in range(len(arms))} == {inertSelf[0][i]},
                f"{n}: the region's own inertance is one number under either operator"
                f" and either source ({inertSelf[0][i]:.10g})")
    for p, q in PAIRS:
        require(M[at[p]] == M[at[q]],
                f"{p} / {q}: the two operators install the SAME inertance table"
                f" ({M[at[p]][0][0]:.10g}, {M[at[p]][1][0]:.10g} kg/m4) - choosing how"
                " to iterate the interface does not choose where the mass comes from")

    # ---- 3. the interface, and 4. the rows it solved ----
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
                f"{a['name']}: the network still charges R times the rate under"
                f" {'quasi-Newton' if a['accelerate'] == 'yes' else 'the plain loop'}"
                f" (worst {w_id:.1e} Pa)")
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
                f" inertance ({rX1:.1e} parallel, {rX2:.1e} series) - {SOURCE[a['name']]},"
                f" {'quasi-Newton' if a['accelerate'] == 'yes' else 'plain loop'}")
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

        # the ramp fitted from the parallel edge's own trajectory
        xs, ys = q[0][1:N_STEPS - 1], q[0][2:N_STEPS]
        k = float(len(xs))
        sx, sy = sum(xs), sum(ys)
        Bc = (k * math.fsum(u * w for u, w in zip(xs, ys)) - sx * sy) / (
            k * math.fsum(u * u for u in xs) - sx * sx)
        qInf = ((sy - Bc * sx) / k) / (1.0 - Bc)
        RX, RBy_ = R[0][-1], pipeChord(*BY, pp["By"][-1])
        mX = math.log(dp[0][-1] / dp[0][0]) / math.log(q[0][-1] / q[0][0])
        mBy = (1.0 + math.log(pipeChord(*BY, pp["By"][-1])
                              / pipeChord(*BY, pp["By"][0]))
               / math.log(pp["By"][-1] / pp["By"][0]))
        slope = mX * RX + mBy * RBy_
        series.append(dict(q=q, R=R, dp=dp, res=res, pas=pas, pp=pp, nd=nd,
                           dp12=dp12, dp23=dp23, B=Bc, tau=-DT / math.log(Bc),
                           qInf=qInf, mX=mX, mBy=mBy,
                           tauS=DT / math.log(1.0 + DT * slope
                                              / (M[a_i][0][0] + M_BY))))

    # ---- 5. the operator axis buys passes, not physics ----
    for p, q in PAIRS:
        sp, sq_ = series[at[p]], series[at[q]]
        dq = abs(sp["q"][0][-1] - sq_["q"][0][-1]) / sp["q"][0][-1]
        require(dq <= 1e-4,
                f"{p} / {q}: the two operators end on the same parallel-branch rate"
                f" ({sp['q'][0][-1]:.6g} against {sq_['q'][0][-1]:.6g} m3/s, {dq:.1e}"
                " relative) - the interface has one fixed point")
        dInf = abs(sp["qInf"] - sq_["qInf"]) / sp["qInf"]
        dTau = abs(sp["tau"] - sq_["tau"]) / sp["tau"]
        require(dInf <= 1e-4 and dTau <= 5e-3,
                f"{p} / {q}: and they fit the same asymptote and the same time constant"
                f" ({sp['qInf']:.6g} against {sq_['qInf']:.6g} m3/s, {dInf:.1e};"
                f" tau {sp['tau']:.6f} against {sq_['tau']:.6f} s, {dTau:.1e})")
        n1 = len(arms[at[p]]["passes"])
        n2 = len(arms[at[q]]["passes"])
        require(n2 <= n1,
                f"{p} / {q}: the accelerated arm needs no more passes than the plain"
                f" one ({n2} against {n1}) - what the operator buys is the schedule")

    # ---- 6. the source axis still owns its own accounts ----
    for a_i, a in enumerate(arms):
        s = series[a_i]
        require(abs(s["tauS"] / s["tau"] - 1.0) <= 0.01,
                f"{a['name']}: the ramp's time constant is (M_X1+M_By) over the partners'"
                f" SLOPES = {s['tauS']:.6f} s against the fitted {s['tau']:.6f} s"
                f" ({(s['tauS'] / s['tau'] - 1) * 100:+.3f}%), with X1 on"
                f" {SOURCE[a['name']]} ({M[a_i][0][0]:.6g} kg/m4)")
        RX = s["R"][0][-1]
        RBy_ = pipeChord(*BY, s["pp"]["By"][-1])
        require(abs(s["qInf"] / (Q_IN * RBy_ / (RBy_ + RX)) - 1.0) <= 0.02,
                f"{a['name']}: and it relaxes towards the resistive split of the settled"
                f" chords, {s['qInf']:.6g} m3/s")
    tauBore = series[at["bore_plain"]]["tau"]
    tauRegion = series[at["region_plain"]]["tau"]
    mRatio = (M[at["bore_plain"]][0][0] + M_BY) / (M[at["region_plain"]][0][0] + M_BY)
    require(abs(tauBore / tauRegion / mRatio - 1.0) <= 0.02,
            f"the source axis moves tau by exactly the inertance ratio: {tauBore:.6f} s /"
            f" {tauRegion:.6f} s = {tauBore / tauRegion:.5f} against {mRatio:.5f},"
            " and the same holds for the accelerated pair"
            f" ({series[at['bore_newton']]['tau']:.6f} /"
            f" {series[at['region_newton']]['tau']:.6f} s)")
    spread = max(s["qInf"] for s in series) - min(s["qInf"] for s in series)
    require(spread / min(s["qInf"] for s in series) <= 2e-4,
            f"all four arms relax to one rate: fitted asymptotes"
            f" {' / '.join('%.6g' % s['qInf'] for s in series)} m3/s, spread"
            f" {spread / min(s['qInf'] for s in series) * 100:.4f}% - the destination"
            " holds neither M nor the operator")

    # ---- 7. the starting kick follows the source, not the operator ----
    kick = {k: M[at[k]][1][0] * Q_IN / DT for k in ARMS}
    for a_i, a in enumerate(arms):
        s = series[a_i]
        err = abs(s["dp23"][0] - (s["R"][1][0] * Q_IN + kick[a["name"]])) / s["dp23"][0]
        require(err <= 1e-12,
                f"{a['name']}: the series edge's first node drop is friction at Q_IN plus"
                f" its OWN inertance's kick ({s['dp23'][0]:.6f} ="
                f" {s['R'][1][0] * Q_IN:.4f} + {kick[a['name']]:.4f}, {err:.1e})")
    for p, q in PAIRS:
        d = abs(series[at[p]]["dp23"][0] - series[at[q]]["dp23"][0])
        require(d <= 1e-5 * abs(series[at[p]]["dp23"][0]),
                f"{p} / {q}: both operators pay the same kick on the pinned edge"
                f" ({series[at[p]]['dp23'][0]:.6f} against {series[at[q]]['dp23'][0]:.6f}"
                f" Pa, {d:.3e} apart) - the rate there is never an unknown")
    dAcross = series[at["bore_plain"]]["dp23"][0] - series[at["region_plain"]]["dp23"][0]
    require(abs(dAcross - (kick["bore_plain"] - kick["region_plain"]))
            <= 1e-3 * abs(dAcross),
            f"while the source decides it: the bore pair kicks {series[at['bore_plain']]['dp23'][0]:.6f} Pa"
            f" and the mesh pair {series[at['region_plain']]['dp23'][0]:.6f} Pa, a difference of"
            f" {dAcross:.6f} Pa against dM*Q/dt ="
            f" {kick['bore_plain'] - kick['region_plain']:.6f} Pa")
    require(all(max(abs(x - Q_IN) for x in s["q"][1]) / Q_IN <= 1e-14 for s in series),
            "the series edge's rate is pinned in all four arms, which is why its"
            " nonlinearity has no exponent to read and only its inertance is per-arm")

    # ---- 8. what quasi-Newton does to the edge it was not asked about ----
    for p, q in PAIRS:
        rp, rq = series[at[p]]["res"], series[at[q]]["res"]
        fp, fq = max(rp[1]), max(rq[1])
        require(fq > 10.0 * fp and fq <= TOL,
                f"{p} / {q}: quasi-Newton negotiates the interface as a WHOLE, so the"
                f" pinned edge rides along with a worst residual of {fq:.3e} against the"
                f" plain loop's {fp:.3e} ({fq / fp:.0f}x) - still inside the case"
                f" tolerance {TOL:.0e}, and the edge that is negotiated is unchanged"
                f" ({max(rp[0]):.3e} against {max(rq[0]):.3e})")
    ex = [s["mX"] for s in series]
    require(max(ex) - min(ex) <= 0.01,
            f"X1's block answers with one law under every arm: dp = c*q^"
            f"{' / '.join('%.4f' % e for e in ex)} - neither the source nor the operator"
            " touches the slope")

    # ---- 9. schedule and cost ----
    for a in arms:
        idx = {int(value(p, "worst")) for p in a["passes"]}
        require(idx == {0},
                f"{a['name']}: the schedule stops on the negotiated edge in all"
                f" {len(a['passes'])} passes")
        s = series[at[a["name"]]]
        require(sum(s["pas"]) == len(a["passes"]),
                f"{a['name']}: one pass is one network solve - {sum(s['pas'])} accepted"
                f" of {len(a['passes'])} run, {len(a['passes']) - N_STEPS} over tolerance")
        require(all(1 <= p <= MAX_PASSES for p in s["pas"]),
                f"{a['name']}: every level agrees inside the allowed sub-iterations"
                f" (max {max(s['pas'])} of {MAX_PASSES})")
    nPass = [len(a["passes"]) for a in arms]
    require(all(a["solves"] == 2 * len(a["passes"]) for a in arms),
            f"cost: {' / '.join('%d passes' % p for p in nPass)} interface passes for"
            f" {4 * N_STEPS} agreed time levels and {sum(a['solves'] for a in arms)}"
            " block solves - every pass runs every partner once, whichever source and"
            " whichever operator that arm was given")

    # ---- 10. each block's own field, quoted from outside ----
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
    print(f"note: sources read {' / '.join(SOURCE[n] for n in ARMS)} and operators"
          f" {' / '.join(OPERATOR[n] for n in ARMS)} across the four arms")

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: inertia attribution and the interface operator are independent axes"
          " of one transient 1D-3D interface")


if __name__ == "__main__":
    sys.exit(main())
