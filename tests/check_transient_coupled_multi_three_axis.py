#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient-multi over two RANS detailed regions with
all three degrees of freedom this chain has been turning one at a time present at
once: where each edge's region inertia comes from (per arm AND per edge), and
which operator negotiates the interface.

transientCoupledMultiTurbulentInertiaNewton crossed the two WHOLE-INTERFACE
sources with the operator and got a clean 2x2, but it left the per-edge table out
of the accelerated half: there the claim "this edge trusts the block, that edge
trusts the drawing, under quasi-Newton" was an inference from two results that
had never been run together. This case makes it a measurement - 8 arms = 4
inertia tables (both bore, both region, X1 from the region with X2 from the
drawing, and the mirror) x 2 interface operators, on one topology, with identical
seeds.

What has to survive the crossing:

  1. eight arms, three axes   - the runs differ only in `inertias`,
                                 `inertiaFrom{}` and `accelerates`; both partners
                                 are k-omega SST; all eight seeds identical
  2. what got installed       - each edge carries the inertance ITS source names,
                                 every step, in every arm: the per-edge table is
                                 installed verbatim under the quasi-Newton map
                                 too, and the operator never enters it
  3. the interface            - dpNetwork == R*Q exactly in all eight arms, and
                                 every arm's momentum rows close with THAT
                                 EDGE's inertance - an accepted accelerated state
                                 is still a state of the momentum system the
                                 table says it is
  4. per-edge attribution     - under Gauss-Seidel an edge's whole accepted
                                 history is the history of the arm that gave that
                                 edge the same source, to the last bit, whatever
                                 the other edge was given; under quasi-Newton the
                                 edges negotiate simultaneously, so separability
                                 becomes ordering rather than identity, and this
                                 checker measures how far
  5. the operator axis        - one table, two operators: the same end rate, the
                                 same asymptote, the same tau; passes, not physics
  6. the source axis          - tau still is (M_X1 + M_By) over the sum of SLOPES,
                                 the tau ratio still equals the inertance ratio,
                                 the destination still holds no M - eight arms
  7. the starting kick        - the pinned series edge pays M*Q/dt out of the
                                 source ITS edge was given, under either operator
  8. floors and laws          - the residual floor the pinned edge rides along
                                 with, and the partner's own exponent, which
                                 neither axis touches
  9. schedule, cost           - who pays for the interface, in all eight arms
 10. the blocks' own fields   - wall-resolved closure, quoted from outside this
                                 repository, arm by arm

usage: check_transient_coupled_multi_three_axis.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# fluid and geometry of tests/cases/transientCoupledMultiThreeAxis/case.dict
MU, RHO = 1.0e-3, 1000.0
Q_IN = 1.2e-4
P_REF = 1.0e6
DT, N_STEPS = 0.05, 12
MAX_PASSES = 6
TOL = 2e-4
NAME = ("X1", "X2")
ARMS = ("bore_gs", "region_gs", "mixX1_gs", "mixX2_gs",
        "bore_qn", "region_qn", "mixX1_qn", "mixX2_qn")
# what each arm's per-edge table says, which is what the interface installs
TABLE = {"bore_gs": ("network", "network"),
         "region_gs": ("partner", "partner"),
         "mixX1_gs": ("partner", "network"),
         "mixX2_gs": ("network", "partner"),
         "bore_qn": ("network", "network"),
         "region_qn": ("partner", "partner"),
         "mixX1_qn": ("partner", "network"),
         "mixX2_qn": ("network", "partner")}
OPERATOR = dict((a, "no" if a.endswith("_gs") else "yes") for a in ARMS)
# the whole-interface default each arm declares in `inertias`, which a mixed
# arm's per-edge table then overrides edge by edge
DEFAULT = {"bore_gs": "network", "region_gs": "partner",
           "mixX1_gs": "network", "mixX2_gs": "network",
           "bore_qn": "network", "region_qn": "partner",
           "mixX1_qn": "network", "mixX2_qn": "network"}
GS = ("bore_gs", "region_gs", "mixX1_gs", "mixX2_gs")
QN = ("bore_qn", "region_qn", "mixX1_qn", "mixX2_qn")
# inside one operator group, the other arm that gave edge i the same source
TWIN = {}
for _a in ARMS:
    for _i in (0, 1):
        _c = [b for b in ARMS if b != _a and OPERATOR[b] == OPERATOR[_a]
              and TABLE[b][_i] == TABLE[_a][_i]]
        assert len(_c) == 1, (_a, _i, _c)
        TWIN[(_a, _i)] = _c[0]
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
    """One dict per arm: the table and operator it was given, its accepted
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
        raise SystemExit("transient multi three-axis solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    arms = parse(lines)
    at = {a["name"]: k for k, a in enumerate(arms)}

    # ---- 1. eight arms, three axes ----
    require([a["name"] for a in arms] == list(ARMS),
            f"all eight arms ran, in the order declared (got {[a['name'] for a in arms]})")
    for a in arms:
        require(sorted(a["steps"]) == list(range(1, N_STEPS + 1)),
                f"{a['name']}: all {N_STEPS} time levels reported an accepted state")
        require(sorted(a["source"]) == list(NAME),
                f"{a['name']}: the run reported a source for each coupled edge")
        require(tuple(a["source"][n] for n in NAME) == TABLE[a["name"]],
                f"{a['name']}: installed {a['source']['X1']} on X1 and"
                f" {a['source']['X2']} on X2, which is the table this arm declares")
        require(a["accelerate"] == OPERATOR[a["name"]],
                f"{a['name']}: negotiated by the"
                f" {'quasi-Newton' if a['accelerate'] == 'yes' else 'plain Gauss-Seidel'}"
                " map this arm declares")
        require(all(word(b, "turb") == "kOmegaSST" for b in a["blocks"]),
                f"{a['name']}: both partners faced the interface as RANS blocks")
        require(a["inertia"] == DEFAULT[a["name"]],
                f"{a['name']}: the header's single `inertia {a['inertia']}` word is"
                " the arm's whole-interface default, which for a mixed arm is not"
                f" what X1 was given ({a['source']['X1']}) - the per-edge table is"
                " what got installed")
    require(len({TABLE[a] for a in ARMS}) == 4,
            "the eight arms carry four different per-edge tables - mixing is not a"
            " re-run of either whole-interface choice")
    require(len({OPERATOR[a] for a in ARMS}) == 2,
            "and every table is run under both operators - a complete 4x2, not eight"
            " near-misses")
    require(all(TABLE[a][0] != TABLE[a][1] for a in ARMS if a.startswith("mix")),
            "the mixed arms really do name two different sources in one momentum"
            " system")

    # ---- 2. what each edge was given, under either operator ----
    M = [{i: [value(a["steps"][n][i], "M") for n in range(1, N_STEPS + 1)]
          for i in (0, 1)} for a in arms]
    inertSelf = [[value(b, "inertance") for b in a["blocks"]] for a in arms]
    for a_i, a in enumerate(arms):
        for i, n in enumerate(NAME):
            L, H, W = EDGES[n]
            exact = RHO * L / (H * W)
            bore = RHO * L / (0.25 * math.pi * BORE[n] ** 2)
            want = bore if TABLE[a["name"]][i] == "network" else inertSelf[a_i][i]
            require(max(abs(x - want) for x in M[a_i][i]) <= 1e-9 * abs(want),
                    f"{a['name']} {n}: its momentum row carries"
                    f" {'the bore rho*L/A = %.10g' % bore if TABLE[a['name']][i] == 'network' else 'the block-measured rho*V/A^2 = %.10g' % exact}"
                    " kg/m4, every step")
            require(abs(inertSelf[a_i][i] - exact) <= 1e-9 * exact,
                    f"{a['name']} {n}: its block measures {inertSelf[a_i][i]:.10g} kg/m4"
                    f" on its own mesh, equal to rho*L/(H*W) = {exact:.10g}")
    for i, n in enumerate(NAME):
        L, H, W = EDGES[n]
        ratio = (RHO * L / (0.25 * math.pi * BORE[n] ** 2)) / (RHO * L / (H * W))
        require(abs(M[at["bore_gs"]][i][0] / M[at["region_gs"]][i][0] - ratio) <= 1e-9,
                f"{n}: the bore overstates this region by {ratio:.6f}"
                f" (H/W {W / H:.1f}), so the two sources name different numbers here")
        require(abs(ratio - 1.0) > 1e-3,
                f"{n}: the two sources are genuinely apart here ({ratio:.6f}),"
                " so naming them edge by edge is not naming one number twice")
    for i, n in enumerate(NAME):
        bySource = {}
        for a_i, a in enumerate(arms):
            bySource.setdefault(TABLE[a["name"]][i], set()).add(M[a_i][i][0])
        one = dict((k, next(iter(v))) for k, v in bySource.items())
        require(all(len(v) == 1 for v in bySource.values()),
                f"{n}: one source, one number - all eight arms install"
                f" {one['network']:.10g} for network and {one['partner']:.10g}"
                " kg/m4 for partner, whatever the operator or the other edge")

    # ---- 3. the interface, and the rows it solved ----
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
        opWord = "the quasi-Newton map" if a["accelerate"] == "yes" else "a plain loop"
        require(w_id == 0.0,
                f"{a['name']}: the network still charges R times the rate, under a"
                f" per-edge table and {opWord} both (worst {w_id:.1e} Pa)")
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

    # ---- 4. per-edge attribution under each operator ----
    def sep(a, b, i):
        """Relative distance between two arms on edge i, in the quantity THAT edge
        actually moves: its rate where the topology lets it negotiate one, its node
        drop where its rate is pinned by the series position."""
        if i == 0:
            return (max(abs(series[at[a]]["q"][0][k] - series[at[b]]["q"][0][k])
                        for k in range(N_STEPS)) / Q_IN)
        return (max(abs(series[at[a]]["dp23"][k] - series[at[b]]["dp23"][k])
                    for k in range(N_STEPS)) / abs(series[at[a]]["dp23"][0]))

    for grp, op in ((GS, "the plain loop"), (QN, "quasi-Newton")):
        for a in grp:
            for i, n in enumerate(NAME):
                twin = TWIN[(a, i)]
                bit = all(arms[at[a]]["steps"][s][i] == arms[at[twin]]["steps"][s][i]
                          for s in range(1, N_STEPS + 1))
                dq = sep(a, twin, i)
                rivals = [b for b in grp if TABLE[b][i] != TABLE[a][i]]
                dr = min(sep(a, b, i) for b in rivals)
                what = "rate" if i == 0 else "node drop"
                require(dq <= 1e-4 and dr >= 1e-3 and dr / max(dq, 1e-16) >= 1e3,
                        f"{a} {n}: under {op} the arm that gave {n} the SAME source"
                        f" ({twin}) is {dq:.1e} of the {what} away, the nearest arm"
                        f" that gave it the other source is {dr:.3e} - attribution"
                        f" separates by {dr / max(dq, 1e-16):.1e}x"
                        f" ({'and the two states are bit-identical lines' if bit else 'not bit-identical lines'})")
    for a, w1, w2 in (("mixX1_gs", "region_gs", "bore_gs"),
                      ("mixX2_gs", "bore_gs", "region_gs"),
                      ("mixX1_qn", "region_qn", "bore_qn"),
                      ("mixX2_qn", "bore_qn", "region_qn")):
        d1, d2 = sep(a, w1, 0), sep(a, w2, 1)
        require(d1 <= 1e-4 and d2 <= 1e-4,
                f"{a} IS the pair pasted together: X1 is {w1}'s X1 to {d1:.1e} of the"
                f" rate and X2 is {w2}'s X2 to {d2:.1e} of the drop, in ONE run - the"
                " two rows of one momentum system written from two different books")
    d1, d2 = sep("mixX1_gs", "mixX2_gs", 0), sep("mixX1_gs", "mixX2_gs", 1)
    require(d1 >= 1e-3 and d2 >= 1e-3,
            f"and the two mirror tables are far apart on BOTH edges ({d1:.3e} of X1's"
            f" rate, {d2:.3e} of X2's drop) - neither mixed arm is a compromise"
            " between the two whole-interface choices")

    # ---- 5. the operator axis: one table, two maps ----
    for a in GS:
        b = a[:-2] + "qn"
        p, s = series[at[a]], series[at[b]]
        dq = max(abs(p["q"][i][k] - s["q"][i][k])
                 for i in (0, 1) for k in range(N_STEPS)) / Q_IN
        require(dq <= 2e-4,
                f"{a} against {b}: the same table under the accelerated map ends on the"
                f" same rates ({dq:.3e} of the injected rate over all {N_STEPS} levels"
                " and both edges)")
        require(sum(s["pas"]) <= sum(p["pas"]),
                f"{b} buys the interface no more expensive than {a}"
                f" ({sum(s['pas'])} against {sum(p['pas'])} passes)")

    # ---- 6. the source axis, per arm ----
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
                f" ({M[a_i][0][0]:.6g} kg/m4), X2 on {a['source']['X2']} and the"
                f" {'quasi-Newton' if a['accelerate'] == 'yes' else 'plain'} map")
        require(abs(qInf / (Q_IN * RBy / (RBy + RX)) - 1.0) <= 0.02,
                f"{a['name']}: and it relaxes towards the resistive split of the settled"
                f" chords, {qInf:.6g} m3/s")
    for op in (GS, QN):
        require(abs(series[at[op[0]]]["tau"] - series[at[op[1]]]["tau"]) > 1e-6,
                f"{op[0]}/{op[1]}: the two whole-interface sources still give two"
                f" different ramps under the {'plain' if op is GS else 'accelerated'} map"
                f" ({series[at[op[0]]]['tau']:.6f} against {series[at[op[1]]]['tau']:.6f} s)")
        mRatio = (M[at[op[0]]][0][0] + M_BY) / (M[at[op[1]]][0][0] + M_BY)
        require(abs(series[at[op[0]]]["tau"] / series[at[op[1]]]["tau"] / mRatio
                    - 1.0) <= 0.02,
                f"and their tau ratio is the inertance ratio in that group too:"
                f" {series[at[op[0]]]['tau'] / series[at[op[1]]]['tau']:.5f} against"
                f" {mRatio:.5f}")
        for mix, whole in ((op[2], op[1]), (op[3], op[0])):
            tM, tW = series[at[mix]]["tau"], series[at[whole]]["tau"]
            rel = abs(tM / tW - 1.0)
            say = ("at the plain loop's own floor" if op is GS
                   else "only as far as the simultaneous solve reaches")
            require(rel <= (1e-5 if op is GS else 5e-3),
                    f"{mix}: its X1 is {whole}'s X1, so its ramp is {whole}'s ramp"
                    f" {say} - tau {rel:.1e} relative apart"
                    f" ({tM:.6f} against {tW:.6f} s). The edge doing the ramping"
                    " decides; the other edge is worth"
                    f" {abs(series[at[op[0]]]['tau'] / series[at[op[1]]]['tau'] - 1.0):.3f}"
                    " of that")
    spread = max(s["qInf"] for s in series) - min(s["qInf"] for s in series)
    require(spread / min(s["qInf"] for s in series) <= 2e-4,
            f"all eight arms relax to one rate: fitted asymptotes"
            f" {min(s['qInf'] for s in series):.6g} to {max(s['qInf'] for s in series):.6g}"
            f" m3/s, spread {spread / min(s['qInf'] for s in series) * 100:.4f}% - the"
            " destination holds no inertance, whichever axis moved")

    # ---- 7. the starting kick follows the edge, not the arm ----
    kick = [M[i][1][0] * Q_IN / DT for i in range(len(arms))]
    for i, a in enumerate(series):
        err = abs(a["dp23"][0] - (a["R"][1][0] * Q_IN + kick[i])) / a["dp23"][0]
        require(err <= 1e-12,
                f"{arms[i]['name']}: the series edge's first node drop is friction at"
                f" Q_IN plus its OWN inertance's kick ({a['dp23'][0]:.6f} ="
                f" {a['R'][1][0] * Q_IN:.4f} + {kick[i]:.4f}, {err:.1e})")
    for op in (GS, QN):
        g, r = series[at[op[0]]], series[at[op[1]]]
        m1, m2 = series[at[op[2]]], series[at[op[3]]]
        require(abs(m1["dp23"][0] - g["dp23"][0]) <= 1e-3 * abs(g["dp23"][0]),
                f"{op[2]} pays the BORE's kick on X2 under the"
                f" {'plain' if op is GS else 'accelerated'} map"
                f" ({m1['dp23'][0]:.6f} Pa against the bore arm's {g['dp23'][0]:.6f})"
                f" even though it gave X1 to the region ({r['dp23'][0]:.6f} Pa there)")
        require(abs(m2["dp23"][0] - r["dp23"][0]) <= 1e-3 * abs(r["dp23"][0]),
                f"{op[3]} pays the REGION's kick on X2"
                f" ({m2['dp23'][0]:.6f} Pa against {r['dp23'][0]:.6f}) even though it"
                f" gave X1 to the drawing ({g['dp23'][0]:.6f} Pa there)")
        require(abs(g["dp23"][0] - r["dp23"][0]) > 1e3
                and abs(m1["dp23"][0] - m2["dp23"][0]
                        - (kick[at[op[2]]] - kick[at[op[3]]])) <= 1e-3 * abs(kick[at[op[2]]]),
                f"and the two mixed arms really are apart: {m1['dp23'][0]:.6f} against"
                f" {m2['dp23'][0]:.6f} Pa, a difference of"
                f" {m1['dp23'][0] - m2['dp23'][0]:.6f} Pa against dM*Q/dt ="
                f" {kick[at[op[2]]] - kick[at[op[3]]]:.6f} Pa")
    require(all(abs(s["q"][1][k] - Q_IN) <= 1e-14 for s in series for k in range(N_STEPS)),
            "the series edge's rate is pinned in all eight arms, which is why its"
            " nonlinearity has no exponent to read and only its inertance is per-edge")

    # ---- 8. floors the operator rides along with, and the partner's law ----
    for p, s_ in ((0, 4), (1, 5), (2, 6), (3, 7)):
        fp = max(max(series[p]["res"][i]) for i in (0, 1))
        fq = max(max(series[s_]["res"][i]) for i in (0, 1))
        pinned = [max(series[p]["res"][1]), max(series[s_]["res"][1])]
        talked = [max(series[p]["res"][0]), max(series[s_]["res"][0])]
        require(pinned[1] > 5.0 * pinned[0] and fq <= TOL,
                f"TABLE {ARMS[p]}/{ARMS[s_]}: quasi-Newton acts on the interface as a"
                f" whole, so the pinned series edge's floor goes"
                f" {pinned[0]:.3e} to {pinned[1]:.3e} ({pinned[1] / pinned[0]:.0f}x)"
                f" while the edge actually negotiated stays at"
                f" {talked[0]:.3e} / {talked[1]:.3e} - both inside {TOL:.0e}")
    ex = [s["mX"] for s in series]
    require(max(ex) - min(ex) <= 0.01,
            f"X1's block answers with one law under all eight arms:"
            f" dp = c*q^{min(ex):.4f}..{max(ex):.4f} - neither axis touches the slope")

    # ---- 9. schedule and cost ----
    for a_i, a in enumerate(arms):
        s = series[a_i]
        idx = [int(value(p, "worst")) for p in a["passes"]]
        require(set(idx) == {0},
                f"{a['name']}: the schedule stops on the negotiated edge in all"
                f" {len(a['passes'])} passes")
        require(sum(s["pas"]) == len(a["passes"]),
                f"{a['name']}: one pass is one network solve - {sum(s['pas'])} accepted"
                f" of {len(a['passes'])} run, {len(a['passes']) - N_STEPS} over tolerance")
        require(all(1 <= p <= MAX_PASSES for p in s["pas"]),
                f"{a['name']}: every level agrees inside the allowed sub-iterations"
                f" (max {max(s['pas'])} of {MAX_PASSES})")
        require(a["solves"] == 2 * len(a["passes"]),
                f"{a['name']}: every pass runs both partners once")
    nPass = [len(a["passes"]) for a in arms]
    require(sum(a["solves"] for a in arms) == 2 * sum(nPass),
            f"cost: {' / '.join('%d' % p for p in nPass)} passes per arm = {sum(nPass)}"
            f" interface passes for {8 * N_STEPS} agreed time levels and"
            f" {sum(a['solves'] for a in arms)} block solves - under either operator and"
            " any per-edge table")

    # ---- 10. each block's own field ----
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
    print(f"note: the eight arms install"
          f" {sorted({(TABLE[a][0], TABLE[a][1], OPERATOR[a]) for a in ARMS})}")

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: per-edge attribution and the interface operator cross with the two"
          " whole-interface sources - eight arms, one momentum system per table")


if __name__ == "__main__":
    sys.exit(main())
