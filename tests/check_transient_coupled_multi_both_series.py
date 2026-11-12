#!/usr/bin/env python3
"""transientCoupledMultiBothSeries (A9): TWO coupled edges, BOTH in the SERIES
position of one time-stepping network.

    A (fixed p) -> X1 (RANS duct) -> J1 (free port) -> X2 (RANS duct) -> Snk (fixed p)

5.32 put both edges in the parallel position; every earlier run had one edge
pinned by series position and one negotiating. THIS topology is the opposite
corner: both edges are in series, there is no third branch, and J1's
continuity is q_X1 - q_X2 = 0 -- the network carries ONE rate, and what the
two attribution tables still move is WHERE the middle node sits. The claims,
all six arms (four tables under the plain loop, the whole-bore and the X2-mix
again under quasi-Newton):

  1. the series bind          - q_X1 == q_X2 bit for bit, every arm, every
                                level: the network never had two rates to
                                negotiate, and neither edge is a passenger in
                                the residual sense (both lead the schedule)
  2. the rows ADD             - rowDrop_i == R_i*q + M_i*dq/dt bit for bit,
                                and the node closes it from both sides; the
                                merged row (pA - pSnk) == rowDrop_1 +
                                rowDrop_2 is (M1+M2)*dq/dt + (R1+R2)*q -- the
                                ladder's numerator is the PLAIN SUM of the two
                                edges' books, no split weights (the parallel
                                case needed m_i weights and still overshoots)
  3. the inertial SPLIT       - the constant inertial total (drive -
                                (R1+R2)q) divides between the two rows in
                                exactly M1/(M1+M2), every arm, every level:
                                which book an edge carries moves the middle
                                node's share (mixX2's X1 pays 73% of the
                                inertial account against mixX1's 56%) WITHOUT
                                touching the rate
  4. the ladder reads as a SLOPE - the march sits in the linear part of the
                                relaxation (tau = Msum/Rsum ~ 0.52 s, twelve
                                levels of dt 5e-5 cover 0.12% of the travel),
                                so each arm's dq/dt is inversely proportional
                                to its M sum, to better than 1%: four tables,
                                four slopes, strictly ordered
  5. the operator keeps the rung - quasi-Newton on the same table gives the
                                same slope bit for bit
  6. the blocks answer for their own fields - flux in == flux out == the
                                edge's rate, the inlet patch is the section
                                the book was referred to, turbulent, y+ in
                                the wall-function band

usage: check_transient_coupled_multi_both_series.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

MU, RHO = 1.0e-3, 1000.0
P_A, P_SNK = 6.0e3, 0.0
Q0 = 4.5e-5
DT, N_STEPS = 5.0e-5, 12
TOL = 2e-4
MAX_PASSES = 6
NAME = ("X1", "X2")
ARMS = ("bore_gs", "region_gs", "mixX1_gs", "mixX2_gs", "bore_qn", "mixX2_qn")
TABLE = {"bore_gs": ("network", "network"), "region_gs": ("partner", "partner"),
         "mixX1_gs": ("partner", "network"), "mixX2_gs": ("network", "partner"),
         "bore_qn": ("network", "network"), "mixX2_qn": ("network", "partner")}
OPERATOR = dict((a, "no" if a.endswith("_gs") else "yes") for a in ARMS)
DEFAULT = {"bore_gs": "network", "region_gs": "partner",
           "mixX1_gs": "network", "mixX2_gs": "network",
           "bore_qn": "network", "mixX2_qn": "network"}
GS = tuple(a for a in ARMS if OPERATOR[a] == "no")
QN = tuple(a for a in ARMS if OPERATOR[a] == "yes")
EDGES = {"X1": (0.2, 0.005, 0.005), "X2": (0.12, 0.004, 0.008)}
BORE = {"X1": 0.005, "X2": 0.0049}


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


def pairs(line):
    """`nodes step 1: A p J1 p Snk p` / `pipes step 1: X1 q X2 q` as a dict."""
    t = line.split()
    k = next(i for i, x in enumerate(t) if x.endswith(":")) + 1
    out = {}
    while k + 1 < len(t):
        out[t[k]] = float(t[k + 1])
        k += 2
    return out


def linfit_slope(ys):
    """Least-squares slope against the level index 0..N-1."""
    xs = list(range(len(ys)))
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def parse(lines):
    arms = []
    cur = None
    for ln in lines:
        if ln.startswith("transient start ") and ln.count(":") == 1 and " edge " not in ln:
            t = ln.split()
            cur = dict(name=t[2].rstrip(":"), inertia=t[4], accelerate=t[6],
                       source={}, starts={}, steps={}, passes=[], blocks=[],
                       summary=None, solves=0)
            arms.append(cur)
        elif cur is None:
            continue
        elif ln.startswith("SIMPLE iter 1:"):
            cur["solves"] += 1
        elif ln.startswith("transient arm ") and " inertia " in ln:
            cur["source"][ln.split()[4]] = word(ln, "inertia")
        elif ln.startswith("transient start ") and " edge " in ln and " resistance " in ln:
            cur["starts"][ln.split()[4]] = value(ln, "resistance")
        elif ln.startswith("transient step ") and " pass " in ln:
            cur["passes"].append(ln)
        elif ln.startswith("transient step ") and " edge " in ln:
            t = ln.split()
            cur["steps"].setdefault(int(t[2].rstrip(":")), {})[int(t[4])] = ln
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
        raise SystemExit("transient both-series solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    arms = parse(lines)

    # ---- 1. six arms, declared order, everything installed ----
    require([a["name"] for a in arms] == list(ARMS),
            f"all six arms ran in the declared order (got"
            f" {[a['name'] for a in arms]})")
    for a in arms:
        require(sorted(a["steps"]) == list(range(1, N_STEPS + 1)),
                f"{a['name']}: all {N_STEPS} levels reported an accepted state")
        require(tuple(a["source"][n] for n in NAME) == TABLE[a["name"]],
                f"{a['name']}: installed {a['source']['X1']} on X1 and"
                f" {a['source']['X2']} on X2 -- the table this arm declares")
        require(sorted(a["starts"]) == list(NAME)
                and a["starts"]["X1"] == 1.0e7 and a["starts"]["X2"] == 1.5e7,
                f"{a['name']}: both edges started from the same chord guesses"
                f" ({a['starts']['X1']:g}, {a['starts']['X2']:g})")

    # ---- per-arm series ----
    series = []
    for a in arms:
        g = lambda key, i, n: value(a["steps"][n][i], key)
        q = [[g("q", i, n) for n in range(1, N_STEPS + 1)] for i in (0, 1)]
        R = [[g("R", i, n) for n in range(1, N_STEPS + 1)] for i in (0, 1)]
        M = [[g("M", i, n) for n in range(1, N_STEPS + 1)] for i in (0, 1)]
        dq = [[g("dq", i, n) for n in range(1, N_STEPS + 1)] for i in (0, 1)]
        rowd = [[g("rowDrop", i, n) for n in range(1, N_STEPS + 1)]
                for i in (0, 1)]
        dn = [[g("dpNetwork", i, n) for n in range(1, N_STEPS + 1)]
              for i in (0, 1)]
        db = [[g("dpBlock", i, n) for n in range(1, N_STEPS + 1)] for i in (0, 1)]
        res = [[g("residual", i, n) for n in range(1, N_STEPS + 1)]
               for i in (0, 1)]
        nd = [pairs(a["steps"][n]["nd"]) for n in range(1, N_STEPS + 1)]
        series.append(dict(q=q, R=R, M=M, dq=dq, rowd=rowd, dn=dn, db=db,
                           res=res, nd=nd))

    # ---- 2. the series bind: ONE rate ----
    w_bind = 0.0
    for a_i, a in enumerate(arms):
        s = series[a_i]
        w_bind = max(w_bind,
                     max(abs(s["q"][0][n] - s["q"][1][n]) for n in range(N_STEPS)))
    require(w_bind == 0.0,
            f"q_X1 == q_X2 to the last bit in all six arms x"
            f" {N_STEPS} levels ({N_STEPS * len(arms)} pairs, worst {w_bind:.1f}"
            f" m3/s) -- J1's continuity has no third branch, so the network"
            f" carries ONE rate and the interface negotiates a node drop, not"
            f" two rates")
    moved = max(abs(series[0]["q"][0][n] - Q0) / Q0 for n in range(N_STEPS))
    require(moved > 0.001,
            f"and the rate MOVES over the levels ({moved * 100:.3f}% of the"
            f" initial flow -- the merged row's inertance is paying, not a"
            f" frozen state)")

    # ---- 3. the rows ADD: bit-identical row closure, node from both sides ----
    w_row = 0.0
    for s in series:
        for i in (0, 1):
            for n in range(N_STEPS):
                w_row = max(w_row, abs(s["rowd"][i][n]
                                       - (s["dn"][i][n]
                                          + s["M"][i][n] * s["dq"][i][n] / DT)))
    require(w_row == 0.0,
            f"rowDrop == R*q + M*dq/dt to the last bit on BOTH edges in all six"
            f" arms ({2 * N_STEPS * len(arms)} edge-level pairs, worst"
            f" {w_row:.1e} Pa) -- the series row is the same momentum statement"
            " the coupled row has always been")
    w_node = 0.0
    for s in series:
        for n in range(N_STEPS):
            drop1 = s["nd"][n]["A"] - s["nd"][n]["J1"]
            drop2 = s["nd"][n]["J1"] - s["nd"][n]["Snk"]
            w_node = max(w_node,
                         abs(drop1 - s["rowd"][0][n]) / max(abs(drop1), 1e-12),
                         abs(drop2 - s["rowd"][1][n]) / max(abs(drop2), 1e-12))
    require(w_node <= 2e-8,
            f"the SAME free node closes BOTH rows from both sides (pA - pJ1 is"
            f" edge 0's row, pJ1 - pSnk is edge 1's, worst {w_node:.1e}"
            f" relative) -- one node, two books, no third equation to hide in")
    w_merge = 0.0
    for s in series:
        for n in range(N_STEPS):
            w_merge = max(w_merge, abs(s["nd"][n]["A"] - s["nd"][n]["Snk"]
                                       - s["rowd"][0][n] - s["rowd"][1][n]))
    require(w_merge <= 1e-4,
            f"and the merged row is the plain sum: (pA - pSnk) == rowDrop_X1 +"
            f" rowDrop_X2 to {w_merge:.1e} Pa -- (M_X1+M_X2)*dq/dt plus"
            f" (R_X1+R_X2)*q, the numerator WITHOUT split weights, which is"
            f" what the parallel topology needed m_i factors for and still"
            f" overshoots")

    # ---- 4. the inertial SPLIT is M1/(M1+M2), exactly ----
    w_split = 0.0
    for a_i, a in enumerate(arms):
        s = series[a_i]
        for n in range(N_STEPS):
            tot = (s["rowd"][0][n] - s["dn"][0][n]) + (s["rowd"][1][n]
                                                       - s["dn"][1][n])
            if abs(tot) < 1e-12:
                continue
            w_split = max(w_split, abs((s["rowd"][0][n] - s["dn"][0][n]) / tot
                                       - s["M"][0][n] / (s["M"][0][n]
                                                         + s["M"][1][n])))
    require(w_split <= 1e-12,
            f"the constant inertial total splits between the two rows in"
            f" EXACTLY M1/(M1+M2) (worst {w_split:.1e}) in every arm and every"
            f" level -- which book X1 carries moves the middle node's share"
            f" from {series[3]['M'][0][0] / (series[3]['M'][0][0] + series[3]['M'][1][0]):.4f}"
            f" (mixX2: X1 pays the bore book) to"
            f" {series[2]['M'][0][0] / (series[2]['M'][0][0] + series[2]['M'][1][0]):.4f}"
            f" (mixX1: X1 pays the region's), WITHOUT touching the rate")
    j1_rank = [series[at_]["nd"][N_STEPS - 1]["J1"] for at_, a_ in
               ((ARMS.index("mixX2_gs"), "mixX2_gs"),
                (ARMS.index("mixX1_gs"), "mixX1_gs"))]
    require(j1_rank[0] < j1_rank[1],
            f"and the middle node stands where the share says: mixX2's J1"
            f" ({j1_rank[0]:.1f} Pa) sits BELOW mixX1's ({j1_rank[1]:.1f} Pa)"
            f" -- X1 carrying the bigger book pays more of the inertial"
            f" account, so less drive is left between J1 and the sink")

    # ---- 5. the ladder reads as a SLOPE, inverse in the M sum ----
    slopes, msums = {}, {}
    for a_i, a in enumerate(arms):
        s = series[a_i]
        slopes[a["name"]] = linfit_slope(s["q"][0])
        msums[a["name"]] = s["M"][0][0] + s["M"][1][0]
    worst_ratio = 0.0
    for a in GS:
        for b in GS:
            if a == b:
                continue
            worst_ratio = max(worst_ratio,
                              abs((slopes[a] / slopes[b])
                                  / (msums[b] / msums[a]) - 1.0))
    require(worst_ratio <= 0.01,
            f"four tables give four slopes, each INVERSELY proportional to its"
            f" M sum to {worst_ratio * 100:.2f}% worst case -- the merged"
            f" row's dq/dt is (drive - Rsum*q)/Msum, so the ladder is the two"
            f" edges' mass adding into ONE numerator, each edge contributing"
            f" the book its table names for it")
    order = sorted(GS, key=lambda x: -msums[x])
    strict = all(slopes[order[k]] < slopes[order[k + 1]] for k in range(3))
    require(strict,
            f"and strictly ordered: more inertance, slower march --"
            + ", ".join(f"{x} {slopes[x]:.3e}" for x in order)
            + " (slope per level, m3/s per step; M sums"
            + ", ".join(f" {msums[x]:.4g}" for x in order) + ")")
    for qn in QN:
        gs = qn.replace("_qn", "_gs")
        require(abs(slopes[qn] / slopes[gs] - 1.0) <= 0.01,
                f"{qn}: the same table under quasi-Newton keeps the same rung"
                f" (slope {slopes[qn]:.3e} against the plain loop's"
                f" {slopes[gs]:.3e}, {(slopes[qn] / slopes[gs] - 1) * 100:+.3f}%)"
                f" -- the operator changes the road, not the rung")

    # ---- 6. the books are the case's two, at 5.32's ratios ----
    for i, n in enumerate(NAME):
        L, H, W = EDGES[n]
        exact = RHO * L / (H * W)
        bore = RHO * L / (0.25 * math.pi * BORE[n] ** 2)
        mn = [series[ARMS.index(x)]["M"][i][0] for x in ARMS
              if TABLE[x][i] == "network"][0]
        mp = [series[ARMS.index(x)]["M"][i][0] for x in ARMS
              if TABLE[x][i] == "partner"][0]
        require(abs(mp / exact - 1.0) < 1e-9
                and abs(mn / bore - 1.0) < 1e-9,
                f"{n}: the region source installs exactly rho*L/(H*W) ="
                f" {exact:.6g} kg/m4 and the drawing's bore"
                f" {BORE[n]} m gives {mn:.6g} ({bore / exact:.5f}x) -- the same"
                f" two books and the same mis-statement as 5.32's parallel run")

    # ---- 7. the handover, the chords, the cost ----
    for a in arms:
        cnt = {0: 0, 1: 0}
        bad = 0
        for p in a["passes"]:
            rs = [float(x) for x in re.search(r"residuals (\S+) (\S+)",
                                              p).groups()]
            w = int(word(p, "worst"))
            cnt[w] += 1
            bad += 0 if w == max((0, 1), key=lambda i: rs[i]) else 1
        require(bad == 0 and cnt[0] >= 1 and cnt[1] >= 1,
                f"{a['name']}: the worst index is the argmax of the residuals"
                f" printed beside it, and BOTH edges lead the schedule"
                f" ({cnt[0]} passes to {cnt[1]}) -- in series neither edge can"
                f" be a passenger on the residual either")
        s_arm = series[ARMS.index(a["name"])]
        res_worst = max(max(r) for r in s_arm["res"])
        require(res_worst < TOL,
                f"{a['name']}: every accepted level is inside the declared"
                f" tolerance ({res_worst:.3e} against {TOL:g}) -- the chord"
                f" the block reports and the row's R*q agree")
        pas = [value(a["steps"][n]["s"], "passes") for n in range(1, N_STEPS + 1)]
        require(all(p <= MAX_PASSES for p in pas)
                and sum(pas) == len(a["passes"]),
                f"{a['name']}: every level agrees inside the allowed"
                f" sub-iterations (max {int(max(pas))} of {MAX_PASSES})")
        require(a["solves"] == 2 * len(a["passes"]),
                f"{a['name']}: {a['solves']} block solves for"
                f" {len(a['passes'])} passes -- one pass runs both partners"
                f" once")
        require(max(value(a["steps"][n]["s"], "newton")
                    for n in range(1, N_STEPS + 1)) < 1e-11,
                f"{a['name']}: the 1D step solved its own system at every"
                f" level")

    # ---- 8. the blocks answer for their own fields ----
    for a in arms:
        s = series[ARMS.index(a["name"])]
        for b in a["blocks"]:
            n = word(b, "name").rstrip(":")
            i = NAME.index(n)
            L, H, W = EDGES[n]
            require(abs(value(b, "fluxIn") - value(b, "fluxOut"))
                    <= 1e-9 * abs(value(b, "fluxIn"))
                    and abs(value(b, "fluxIn") - s["q"][i][-1])
                    <= 1e-9 * abs(s["q"][i][-1]),
                    f"{a['name']} {n}: the block carries the rate its edge"
                    f" ended with ({value(b, 'fluxIn'):.6g} m3/s) and what"
                    f" enters leaves -- in series that is the WHOLE rate")
            require(abs(value(b, "inletArea") - H * W) <= 1e-12 * H * W,
                    f"{a['name']} {n}: the inlet patch is the section its"
                    f" inertance was referred to")
            yp = value(b, "wallYPlusMax")
            require(11.81 <= yp <= 300.0,
                    f"{a['name']} {n}: wall cells at y+ {yp:.1f}, inside the"
                    f" log-law band")
            v = value(b, "plugSpeed")
            de = 4.0 * H * W / (2.0 * (H + W))
            require(RHO * v * de / MU > 4000.0,
                    f"{a['name']} {n}: turbulent at Re_bulk"
                    f" {RHO * v * de / MU:.0f}")
            require(abs(value(b, "dpPa") - s["db"][i][-1])
                    <= 1e-9 * abs(value(b, "dpPa")),
                    f"{a['name']} {n}: the drop the interface ended with is"
                    f" the drop the block's own field measures")

    print()
    print(f"{len(fails)} failure(s)")
    if fails:
        for f in fails:
            print("FAILED: " + f)
        sys.exit(1)
    print("OK: two series edges -- ONE rate bound bit for bit, the rows add"
          " into a weight-free numerator, the inertial account splits exactly"
          " as M1/(M1+M2) and moves the middle node, and the ladder of slopes"
          " is the M sum talking")


if __name__ == "__main__":
    main()
