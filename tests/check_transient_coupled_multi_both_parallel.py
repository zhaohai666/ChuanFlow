#!/usr/bin/env python3
"""Run the 1D3D transient multi-edge coupling with BOTH coupled edges sitting in
the parallel position, and find out what per-edge attribution is worth when
neither edge is a passenger.

Every coupled run in this chain has had one edge that negotiates a rate and one
whose rate is fixed by the series position. That asymmetry is what made §5.31's
per-edge measurement two different quantities - a rate on one edge, a node drop
on the other - and it is what left two things untested:

  * the residual's `worst` edge never changed hands in 224 passes there, because
    one edge had nothing to negotiate. If the schedule is driven by who is
    furthest from its chord, handover has to appear the moment both edges are
    live; if it does not, the schedule is being driven by something else.
  * the ramp's numerator had only one coupled inertance in it. Here BOTH
    branches' mass is in the same time constant, so a per-edge table has to
    contribute two different books to ONE numerator - and a mixed table is then
    predicted to sit BETWEEN the two whole-interface tables rather than copying
    one of them, which is exactly what §5.31 found it did on the pinned
    topology.

So the claims, in both directions and with no approximations:

  1. eight arms, three axes    - four tables x two operators on one topology,
                                 identical seeds, and neither edge pinned: both
                                 rates move and both are strictly inside the
                                 injected rate
  2. what got installed        - one source, one number per edge, in all eight
                                 arms, and four distinct SUMS of the two edges'
                                 inertances - the mixed tables are neither of
                                 the whole ones
  3. the interface             - dpNetwork == R*Q exactly, on both edges, every
                                 level, under either operator
  4. two rows, one node drop   - the SAME J1-J2 closes against edge 0's row with
                                 edge 0's M and against edge 1's row with edge
                                 1's M, at every level, in machine arithmetic
  5. the split is predictable  - the pair of rows plus continuity reproduces
                                 every level's split from (R, M, qPrev) alone;
                                 put the other book's M into either row and the
                                 reconstruction breaks by orders of magnitude,
                                 which is what makes the per-edge source an
                                 identified property and not a relabelling
  6. the ladder              - ONE time constant now carries BOTH branches' mass,
                               so four tables give four strictly ordered ramps
                               whose tau ratios are the ratios of M_X1 + M_X2, the
                               mixed tables sit BETWEEN the whole-interface ones,
                               and the operator moves tau far less than the ladder
                               spans. The analytic tau = (M_X1+M_X2)/sum of secant
                               slopes is 1.48x to 1.53x the measured rate here
                               (the inner passes of one level are their own discrete
                               map), so the ladder is claimed as a ratio, not as a
                               formula
  7. the destination           - eight arms relax to one split, and that split is
                                 the resistive split of the two blocks' own
                                 reported power laws - no inertance in it
  8. the handover              - the worst edge really does change hands, in
                                 every arm, and the index printed is the argmax
                                 of the two residuals printed with it
  9. no passenger              - each edge's own residual floor stays inside a
                                 quarter of the tolerance under either operator and
                                 the two operators differ on it by at most 40x, in
                                 BOTH directions - the §5.31 ride-along inflation
                                 (73x..87x, one-signed) belonged to the pinned
                                 series position, not to quasi-Newton
 10. schedule, cost, fields    - who pays for the interface, one pass = one
                                 network solve, and both wall-resolved blocks
                                 close on their own fields in every arm

usage: check_transient_coupled_multi_both_parallel.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# fluid and geometry of tests/cases/transientCoupledMultiBothParallel/case.dict
MU, RHO = 1.0e-3, 1000.0
Q_IN, P_REF = 1.2e-4, 1.0e6
DT, N_STEPS = 0.05, 12
MAX_PASSES = 6
TOL = 2e-4
NAME = ("X1", "X2")
ARMS = ("bore_gs", "region_gs", "mixX1_gs", "mixX2_gs",
        "bore_qn", "region_qn", "mixX1_qn", "mixX2_qn")
TABLE = {"bore_gs": ("network", "network"), "region_gs": ("partner", "partner"),
         "mixX1_gs": ("partner", "network"), "mixX2_gs": ("network", "partner"),
         "bore_qn": ("network", "network"), "region_qn": ("partner", "partner"),
         "mixX1_qn": ("partner", "network"), "mixX2_qn": ("network", "partner")}
OPERATOR = dict((a, "no" if a.endswith("_gs") else "yes") for a in ARMS)
DEFAULT = {"bore_gs": "network", "region_gs": "partner",
           "mixX1_gs": "network", "mixX2_gs": "network",
           "bore_qn": "network", "region_qn": "partner",
           "mixX1_qn": "network", "mixX2_qn": "network"}
GS = tuple(a for a in ARMS if OPERATOR[a] == "no")
QN = tuple(a for a in ARMS if OPERATOR[a] == "yes")
PAIRS = tuple(zip(GS, QN))
# (length, gap, span) of each detailed region, and the bore the drawing gives it
EDGES = {"X1": (0.2, 0.005, 0.005), "X2": (0.12, 0.004, 0.008)}
BORE = {"X1": 0.005, "X2": 0.0049}
LEADS = {"A": (0.3, 0.02), "B": (0.3, 0.02)}
M_LEAD = RHO * LEADS["A"][0] / (0.25 * math.pi * LEADS["A"][1] ** 2)


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
    """`nodes step 1: Src p J1 p ...` and `pipes step 1: A q X1 q ...` as a dict
    of name to number, so the 1D side is read straight off the network's own
    report."""
    t = line.split()
    k = next(i for i, x in enumerate(t) if x.endswith(":")) + 1
    out = {}
    while k + 1 < len(t):
        out[t[k]] = float(t[k + 1])
        k += 2
    return out


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


def powerLaw(qs, dps):
    """Least squares on log dp = log C + m log q, so each block's own reported
    drop-versus-rate law is read from the outside."""
    xs = [math.log(x) for x in qs]
    ys = [math.log(x) for x in dps]
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    m = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    return m, math.exp((sy - m * sx) / n)


def resistiveSplit(c1, m1, c2, m2, total):
    lo, hi = 1e-9 * total, total * (1.0 - 1e-9)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if c1 * mid ** m1 - c2 * (total - mid) ** m2 > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi) / total


def relax(sp, dt, nSteps):
    """One relaxation rate and the state it relaxes toward, from a split
    fraction history: the rate from consecutive-difference ratios (a geometric
    decay has a constant ratio), the destination from a linear fit at that rate.

    The window band is diagnostic only - each window divides two differences of
    numbers that have already relaxed most of the way, so it amplifies the last
    digits.  What says ONE rate carries the split is the residual of a single
    exponential through all the levels, returned here as `dev`: the worst
    distance of the history from dest + amp*exp(-t/tau), as a fraction of the
    total travel."""
    taus = []
    for i in range(2, nSteps - 1):
        d1 = sp[i] - sp[i - 1]
        d2 = sp[i + 1] - sp[i]
        if d1 * d2 > 0.0 and abs(d1) > 1e-14:
            taus.append(-dt / math.log(abs(d2 / d1)))
    if len(taus) < 3:
        raise SystemExit("the split history does not decay monotonically")
    tau = sum(taus[-4:]) / len(taus[-4:])
    last = taus[-4:]
    spread = (max(last) - min(last)) / tau
    xs = [math.exp(-(i + 1) * dt / tau) for i in range(nSteps)]
    n = len(xs)
    sx, sy = sum(xs), sum(sp)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, sp))
    amp = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    dest = (sy - amp * sx) / n
    dev = max(abs(dest + amp * xs[i] - sp[i]) for i in range(n)) / abs(dest - sp[0])
    return tau, dest, spread, dev


def parse(lines):
    """One dict per arm: the table and operator it was given, its accepted
    states, its passes, and its block reports."""
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
            t = ln.split()
            cur["source"][t[4]] = word(ln, "inertia")
        elif ln.startswith("transient start ") and " edge " in ln and " resistance " in ln:
            cur["starts"][ln.split()[4]] = value(ln, "resistance")
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
        raise SystemExit("transient both-parallel solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    arms = parse(lines)
    at = {a["name"]: k for k, a in enumerate(arms)}

    # ---- 1. eight arms, three axes, and nothing pinned ----
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
        require(sorted(a["starts"]) == list(NAME)
                and abs(a["starts"]["X1"] / 1.0e7 - 1.0) < 1e-12
                and abs(a["starts"]["X2"] / 1.5e7 - 1.0) < 1e-12,
                f"{a['name']}: both edges started from the same chord guesses as"
                f" every other arm ({a['starts']['X1']:g}, {a['starts']['X2']:g})")
    require(len({TABLE[a] for a in ARMS}) == 4,
            "the eight arms carry four different per-edge tables - mixing is not a"
            " re-run of either whole-interface choice")
    require(len({OPERATOR[a] for a in ARMS}) == 2,
            "and every table is run under both operators - a complete 4x2, not eight"
            " near-misses")
    require(all(TABLE[a][0] != TABLE[a][1] for a in ARMS if a.startswith("mix")),
            "the mixed arms really do name two different sources in one momentum"
            " system")

    # ---- the per-arm series, and the two edges' freedom ----
    series = []
    for a_i, a in enumerate(arms):
        q = [[value(a["steps"][n][i], "q") for n in range(1, N_STEPS + 1)] for i in (0, 1)]
        R = [[value(a["steps"][n][i], "R") for n in range(1, N_STEPS + 1)] for i in (0, 1)]
        M = [[value(a["steps"][n][i], "M") for n in range(1, N_STEPS + 1)] for i in (0, 1)]
        dpn = [[value(a["steps"][n][i], "dpNetwork") for n in range(1, N_STEPS + 1)]
               for i in (0, 1)]
        dpb = [[value(a["steps"][n][i], "dpBlock") for n in range(1, N_STEPS + 1)]
               for i in (0, 1)]
        res = [[value(a["steps"][n][i], "residual") for n in range(1, N_STEPS + 1)]
               for i in (0, 1)]
        nd = [pairs(a["steps"][n]["nd"]) for n in range(1, N_STEPS + 1)]
        pp = [pairs(a["steps"][n]["pp"]) for n in range(1, N_STEPS + 1)]
        dp12 = [nd[n]["J1"] - nd[n]["J2"] for n in range(N_STEPS)]
        series.append(dict(q=q, R=R, M=M, dpn=dpn, dpb=dpb, res=res, nd=nd, pp=pp,
                           dp12=dp12))

    for a_i, a in enumerate(arms):
        s = series[a_i]
        far = [max(abs(x - Q_IN) for x in s["q"][i]) / Q_IN for i in (0, 1)]
        moved = [max(abs(s["q"][i][n] - s["q"][i][0]) for n in range(N_STEPS)) / Q_IN
                 for i in (0, 1)]
        require(min(far) > 0.1 and min(moved) > 0.01,
                f"{a['name']}: NEITHER edge is pinned - both rates stay inside the"
                f" injected rate (furthest {far[0]:.3f} / {far[1]:.3f} of it) and"
                f" both travel over the levels ({moved[0]:.4f} / {moved[1]:.4f}) -"
                " this is the topology §5.31 could not put two rates into")

    # ---- 2. one source, one number, per edge, in all eight arms ----
    for i, n in enumerate(NAME):
        L, H, W = EDGES[n]
        exact = RHO * L / (H * W)
        bore = RHO * L / (0.25 * math.pi * BORE[n] ** 2)
        for a_i, a in enumerate(arms):
            require(all(v == series[a_i]["M"][i][0] for v in series[a_i]["M"][i]),
                    f"{a['name']} {n}: the inertance installed on this edge does not"
                    f" move over the {N_STEPS} levels ({series[a_i]['M'][i][0]:.6g})")
        net = set(series[at[a]]["M"][i][0] for a in ARMS if TABLE[a][i] == "network")
        par = set(series[at[a]]["M"][i][0] for a in ARMS if TABLE[a][i] == "partner")
        require(len(net) == 1 and len(par) == 1,
                f"{n}: all four bore arms install one number"
                f" ({net.pop():.6g} kg/m4) and all four region arms another"
                f" ({par.pop():.6g}) - neither the operator nor the other edge's"
                " choice enters this table")
    for i, n in enumerate(NAME):
        L, H, W = EDGES[n]
        exact = RHO * L / (H * W)
        bore = RHO * L / (0.25 * math.pi * BORE[n] ** 2)
        mn = [series[at[a]]["M"][i][0] for a in ARMS if TABLE[a][i] == "network"][0]
        mp = [series[at[a]]["M"][i][0] for a in ARMS if TABLE[a][i] == "partner"][0]
        require(abs(mp / exact - 1.0) < 1e-9,
                f"{n}: the region source installs exactly rho*L/(H*W) ="
                f" {exact:.6g} kg/m4 as measured off its own mesh")
        require(abs(mn / exact - bore / exact) < 1e-9 and bore / exact > 1.05,
                f"{n}: the drawing's bore {BORE[n]} m gives rho*L/(pi d^2/4) ="
                f" {mn:.6g}, i.e. {bore / exact:.6f}x the region's own - the two"
                " books really disagree on this edge")
    sums = {}
    for a in ARMS:
        sM = series[at[a]]["M"][0][0] + series[at[a]]["M"][1][0]
        sums.setdefault(sM, []).append(a)
    require(len(sums) == 4,
            "the eight arms carry four different SUMS of the two edges' inertances"
            f" ({', '.join(f'{k:.6g}' for k in sorted(sums))}) - here both books"
            " land in ONE numerator, so a mixed table cannot be a copy of either"
            " whole one")

    # ---- 3. the interface is still a chord, on both edges, under either operator ----
    w_id = max(abs(s["dpn"][i][n] - s["R"][i][n] * s["q"][i][n])
               for s in series for i in (0, 1) for n in range(N_STEPS))
    require(w_id == 0.0,
            f"dpNetwork == R*q to the last bit on BOTH edges in all eight arms"
            f" ({N_STEPS * 2 * 8} edge-by-level pairs, worst {w_id:.1f} Pa) - the"
            " chord is a property of the exchange, not of the loop that iterates it"
            " or of which book each edge's mass came from")

    # ---- 4. one node drop, two momentum rows, two different M ----
    for a_i, a in enumerate(arms):
        s = series[a_i]
        errs = []
        for i in (0, 1):
            errs.append(max(abs(s["dp12"][n] - (s["R"][i][n] * s["q"][i][n]
                                                + s["M"][i][n]
                                                * (s["q"][i][n]
                                                   - (s["q"][i][n - 1] if n else 0.0))
                                                / DT)) / abs(s["dp12"][n])
                            for n in range(N_STEPS)))
        require(errs[0] <= 1e-9 and errs[1] <= 1e-9,
                f"{a['name']}: the SAME J1-J2 drop closes against edge 0's row with"
                f" edge 0's M ({errs[0]:.1e}) and against edge 1's row with edge 1's"
                f" M ({errs[1]:.1e}) - X1 from {a['source']['X1']}, X2 from"
                " {0}, two books, one momentum system".format(a["source"]["X2"]))
        require(all(s["nd"][n]["Snk"] - P_REF <= 1e-9 for n in range(N_STEPS))
                and all(abs(s["pp"][n]["X1"] + s["pp"][n]["X2"] - Q_IN) <= 1e-12 * Q_IN
                        and abs(s["pp"][n]["A"] - Q_IN) <= 1e-15 * Q_IN
                        and abs(s["pp"][n]["B"] - Q_IN) <= 1e-15 * Q_IN
                        for n in range(N_STEPS)),
                f"{a['name']}: the sink holds its reference and the two negotiated"
                " branches are the whole injected rate - there is no third branch to"
                " hide an imbalance in")
        rLead = max(min([
            abs((s["nd"][n]["Src"] - s["nd"][n]["J1"])
                - (pipeChord(*LEADS["A"], s["pp"][n]["A"]) * s["pp"][n]["A"]
                   + M_LEAD * (s["pp"][n]["A"] - (s["pp"][n - 1]["A"] if n else 0.0))
                   / DT)) / abs(s["nd"][n]["Src"] - s["nd"][n]["J1"]),
            abs((s["nd"][n]["J2"] - s["nd"][n]["Snk"])
                - (pipeChord(*LEADS["B"], s["pp"][n]["B"]) * s["pp"][n]["B"]
                   + M_LEAD * (s["pp"][n]["B"] - (s["pp"][n - 1]["B"] if n else 0.0))
                   / DT)) / abs(s["nd"][n]["J2"] - s["nd"][n]["Snk"])]
            for n in range(N_STEPS)))
        require(rLead <= 1e-12,
                f"{a['name']}: both leads still close against published Churchill"
                f" chords plus their own inertance ({rLead:.1e})")

    # ---- 5. the split is reconstructed from the two rows, and only from those ----
    def reconstruct(a_i, override=None):
        """Solve the pair of rows (equal drop, rates summing to the injected one)
        from the arm's own reported R, M and previous rates. `override` names an
        edge whose M is swapped for the other book's, to see if the history cares."""
        s = series[a_i]
        worst = 0.0
        alt = {"network": "partner", "partner": "network"}
        for n in range(N_STEPS):
            mm = [s["M"][0][n], s["M"][1][n]]
            if override is not None:
                src = arms[a_i]["source"][NAME[override]]
                L, H, W = EDGES[NAME[override]]
                mm[override] = (RHO * L / (0.25 * math.pi * BORE[NAME[override]] ** 2)
                                if alt[src] == "network" else RHO * L / (H * W))
            k = [s["R"][0][n] + mm[0] / DT, s["R"][1][n] + mm[1] / DT]
            qp = [(s["q"][0][n - 1] if n else 0.0), (s["q"][1][n - 1] if n else 0.0)]
            dp = ((Q_IN - mm[0] * qp[0] / DT / k[0] - mm[1] * qp[1] / DT / k[1])
                  / (1.0 / k[0] + 1.0 / k[1]))
            pred = [(dp + mm[0] * qp[0] / DT) / k[0], (dp + mm[1] * qp[1] / DT) / k[1]]
            worst = max(worst, abs(pred[0] - s["q"][0][n]) / Q_IN,
                        abs(pred[1] - s["q"][1][n]) / Q_IN)
        return worst

    for a_i, a in enumerate(arms):
        exact = reconstruct(a_i)
        wrong0 = reconstruct(a_i, 0)
        wrong1 = reconstruct(a_i, 1)
        require(exact <= 1e-10,
                f"{a['name']}: the two rows, each with the M THIS EDGE was given,"
                f" reproduce every level's split ({exact:.1e} of the injected rate)")
        require(wrong0 >= 1e-2 and wrong1 >= 1e-2,
                f"{a['name']}: swap in the other book on X1 and the reconstruction"
                f" misses by {wrong0:.3f}; on X2 by {wrong1:.3f} - {math.log10(wrong0 / max(exact, 1e-18)):.0f}x"
                f" / {math.log10(wrong1 / max(exact, 1e-18)):.0f}x worse, so each"
                " edge's source is an identified property of this run and not a"
                " relabelling")

    # ---- 6. the ladder: four tables, four time constants, strictly ordered ----
    lad = {}
    for a_i, a in enumerate(arms):
        s = series[a_i]
        sp = [s["q"][0][n] / Q_IN for n in range(N_STEPS)]
        tau, dest, spread, dev = relax(sp, DT, N_STEPS)
        m = [powerLaw(s["q"][i], s["dpb"][i]) for i in (0, 1)]
        lad[a["name"]] = dict(tau=tau, dest=dest, spread=spread, dev=dev,
                              m=[m[0][0], m[1][0]],
                              C=[m[0][1], m[1][1]],
                              sumM=s["M"][0][0] + s["M"][1][0])
    for a_i, a in enumerate(arms):
        d = lad[a["name"]]
        require(d["dev"] <= 5e-3,
                f"{a['name']}: one relaxation rate carries the split - a single"
                f" exponential through all {N_STEPS} levels (tau {d['tau']:.5f} s)"
                f" never misses the history by more than {d['dev'] * 100:.3f}% of the"
                f" travel, and the four tail windows agree to"
                f" {d['spread'] * 100:.1f}%")
        require(1.5 <= d["m"][0] <= 1.85 and 1.5 <= d["m"][1] <= 1.85,
                f"{a['name']}: BOTH edges now report their own law - dp ~ q^"
                f"{d['m'][0]:.4f} on X1 and q^{d['m'][1]:.4f} on X2 - where the"
                " pinned edge of 5.31 had no exponent to read")
    for gs, qn in PAIRS:
        d1, d2 = lad[gs], lad[qn]
        require(abs(d1["tau"] / d2["tau"] - 1.0) <= 0.01
                and d1["sumM"] == d2["sumM"],
                f"{gs} against {qn}: the same table under the accelerated map keeps"
                f" the same ladder rung (tau {d1['tau']:.5f} against {d2['tau']:.5f} s,"
                f" {(d1['tau'] / d2['tau'] - 1) * 100:+.3f}%) and the same numerator"
                f" ({d1['sumM']:.6g} kg/m4)")
    for grp in (GS, QN):
        srt = sorted(grp, key=lambda x: lad[x]["sumM"])
        for lo, hi in zip(srt, srt[1:]):
            require(lad[hi]["tau"] > lad[lo]["tau"],
                    f"{lo} against {hi} under the"
                    f" {'plain loop' if grp is GS else 'quasi-Newton'}: more"
                    f" inertance in the numerator, slower ramp - tau"
                    f" {lad[lo]['tau']:.5f} to {lad[hi]['tau']:.5f} s for a numerator"
                    f" ratio of {lad[hi]['sumM'] / lad[lo]['sumM']:.4f}")
        w = max(abs((lad[b]["tau"] / lad[a]["tau"]) / (lad[b]["sumM"] / lad[a]["sumM"])
                    - 1.0) for a in srt for b in srt if a != b)
        if grp is GS:
            require(w <= 0.05,
                    f"the ladder is the numerator talking: across the four plain-loop"
                    f" tables every tau ratio matches the corresponding ratio of"
                    f" M_X1+M_X2 to within {w * 100:.2f}% worst case - the two edges'"
                    " mass adds into ONE time constant, each edge contributing the"
                    " book its table names for it")
    srt = sorted(GS, key=lambda x: lad[x]["sumM"])
    fast, slow = srt[0], srt[-1]
    for a in ("mixX1_gs", "mixX2_gs"):
        require(lad[fast]["tau"] < lad[a]["tau"] < lad[slow]["tau"],
                f"{a} is BETWEEN {fast} ({lad[fast]['tau']:.5f} s) and {slow}"
                f" ({lad[slow]['tau']:.5f} s) at {lad[a]['tau']:.5f} s - a mixed"
                " table is neither of the whole ones here, it takes one M from each")
        other = "mixX2_gs" if a == "mixX1_gs" else "mixX1_gs"
        require((lad[other]["tau"] > lad[a]["tau"])
                == (lad[other]["sumM"] > lad[a]["sumM"]),
                f"{a} against {other}: the two mirror tables land where their"
                f" numerators are ({lad[a]['sumM']:.6g} against"
                f" {lad[other]['sumM']:.6g} kg/m4, tau {lad[a]['tau']:.5f} against"
                f" {lad[other]['tau']:.5f} s) - swapping which edge keeps which book"
                " moves the ramp by the amount the sum says, in the direction the sum"
                " says")

    # ---- 7. the destination holds no inertance ----
    dest = [lad[a]["dest"] for a in ARMS]
    require((max(dest) / min(dest) - 1.0) <= 1e-3,
            f"eight arms, four numerators, two operators - all relax toward one"
            f" split {(max(dest) + min(dest)) / 2:.6f} of the injected rate"
            f" (spread {(max(dest) / min(dest) - 1) * 100:.4f}%)")
    for a in ARMS:
        d = lad[a]
        r = resistiveSplit(d["C"][0], d["m"][0], d["C"][1], d["m"][1], Q_IN)
        require(abs(d["dest"] / r - 1.0) <= 1e-3,
                f"{a}: that destination is the resistive split of the two blocks'"
                f" own reported laws ({r:.6f} predicted against {d['dest']:.6f}"
                f" fitted, {(d['dest'] / r - 1) * 100:+.4f}%) - no M in it")
    rr = [resistiveSplit(lad[a]["C"][0], lad[a]["m"][0], lad[a]["C"][1], lad[a]["m"][1],
                         Q_IN) for a in ARMS]
    require((max(rr) / min(rr) - 1.0) <= 5e-3,
            f"and the resistive split itself barely moves across the eight arms"
            f" ({min(rr):.6f} to {max(rr):.6f}, {(max(rr) / min(rr) - 1) * 100:.4f}%)"
            " - attribution changes the road, the arithmetic of the destination is"
            " the blocks'")

    # ---- 8. the handover ----
    counts = []
    for a in arms:
        cnt = {0: 0, 1: 0}
        bad = 0
        for p in a["passes"]:
            rs = [float(x) for x in re.search(r"residuals (\S+) (\S+)", p).groups()]
            w = int(word(p, "worst"))
            cnt[w] += 1
            bad += 0 if w == max((0, 1), key=lambda i: rs[i]) else 1
        counts.append(cnt)
        require(bad == 0,
                f"{a['name']}: the worst index printed for all {len(a['passes'])}"
                " passes is the argmax of the two residuals printed beside it")
        require(cnt[0] >= 4 and cnt[1] >= 4,
                f"{a['name']}: the handover happens - edge 0 is the worst edge in"
                f" {cnt[0]} passes and edge 1 in {cnt[1]}")
    require(all(c[1] >= 4 for c in counts),
            f"both edges lead the schedule in every one of the eight arms"
            f" (edge 1 worst {min(c[1] for c in counts)} to {max(c[1] for c in counts)}"
            f" times out of {sum(counts[0].values())}) - in 5.31 that index was 0 in"
            " all 224 passes because one edge had nothing to negotiate")

    # ---- 9. no passenger: the floors are the edges' own ----
    ratios = []
    for tbl, gs, qn in [("bore", GS[0], QN[0]), ("region", GS[1], QN[1]),
                        ("mixX1", GS[2], QN[2]), ("mixX2", GS[3], QN[3])]:
        sg, sq = series[at[gs]], series[at[qn]]
        for i, n in enumerate(NAME):
            fg = min(sg["res"][i])
            fq = min(sq["res"][i])
            rat = max(fg / fq, fq / fg)
            ratios.append(fg / fq)
            require(rat <= 40.0 and fg <= 0.25 * TOL and fq <= 0.25 * TOL,
                f"{tbl} table, {n}: both operators bring THIS edge's own residual"
                f" floor to {max(fg, fq):.3e} - no worse than"
                f" {max(fg, fq) / TOL * 100:.1f}% of the declared tolerance - and"
                f" they differ by {rat:.1f}x ({fg:.3e} plain against {fq:.3e}"
                " accelerated). 5.31's pinned edge sat at the same place in every arm"
                " while its floor moved 73x to 87x between the operators, because it"
                " had nothing to negotiate and was reading the other edge's crosstalk;"
                " here neither edge can do that.")
        require(max(max(s["res"][i]) for s in (sg, sq) for i in (0, 1)) < TOL,
                f"{tbl} table: every accepted level on both edges is inside the"
                f" declared tolerance ({max(max(s['res'][i]) for s in (sg, sq) for i in (0, 1)):.3e}"
                f" against {TOL:g})")
    require(min(ratios) < 1.0 / 3.0 and max(ratios) > 3.0,
            f"the operator comparison has no consistent direction across the eight"
            f" edge-table pairs (plain/accelerated floor ratio from {min(ratios):.3f}"
            f" to {max(ratios):.1f}, {sum(1 for r in ratios if r > 1.0)} of"
            f" {len(ratios)} above 1) - a passenger edge inflates one way only, as it"
            " did in 5.31, while here each edge answers for itself")

    # ---- 10. schedule, cost, and the blocks' own fields ----
    tot = 0
    for a in arms:
        pas = [value(a["steps"][n]["s"], "passes") for n in range(1, N_STEPS + 1)]
        tot += len(a["passes"])
        require(all(p <= MAX_PASSES for p in pas)
                and sum(pas) == len(a["passes"]),
                f"{a['name']}: every level agrees inside the allowed sub-iterations"
                f" (max {int(max(pas))} of {MAX_PASSES}) and one pass is one"
                " interface exchange")
        require(a["solves"] == 2 * len(a["passes"]),
                f"{a['name']}: {a['solves']} block solves for {len(a['passes'])}"
                " passes - each pass runs both partners once, and the network solve"
                " is not in that count")
        require(max(value(a["steps"][n]["s"], "newton")
                    for n in range(1, N_STEPS + 1)) < 1e-11,
                f"{a['name']}: the 1D step solved its own system at every level")
        require(all(word(b, "blockConverged") == "1" for b in a["blocks"]),
                f"{a['name']}: both blocks closed their own fields before reporting")
    require(tot == 279,
            f"cost: {tot} interface passes for {N_STEPS * 8} agreed time levels and"
            f" {2 * tot} block solves - with both edges live the plain loop and the"
            " accelerated map pay within a few passes of each other")
    for a in arms:
        for b in a["blocks"]:
            n = word(b, "name").rstrip(":")
            L, H, W = EDGES[n]
            i = NAME.index(n)
            s = series[at[a["name"]]]
            require(abs(value(b, "fluxIn") - value(b, "fluxOut")) <= 1e-9 * Q_IN
                    and abs(value(b, "fluxIn") - s["q"][i][-1]) <= 1e-9 * Q_IN,
                    f"{a['name']} {n}: the block carries the rate its edge ended"
                    f" with ({value(b, 'fluxIn'):.6g} m3/s) and what enters leaves")
            require(abs(value(b, "inletArea") - H * W) <= 1e-12 * H * W,
                    f"{a['name']} {n}: the inlet patch is the section its inertance"
                    " was referred to")
            yp = value(b, "wallYPlusMax")
            require(11.81 <= yp <= 300.0,
                    f"{a['name']} {n}: wall cells at y+ {yp:.1f}, inside the log-law"
                    " band")
            v = value(b, "plugSpeed")
            de = 4.0 * H * W / (2.0 * (H + W))
            require(RHO * v * de / MU > 4000.0,
                    f"{a['name']} {n}: turbulent at Re_bulk"
                    f" {RHO * v * de / MU:.0f}")
            closure = abs(value(b, "probeSlope")) * H / 2.0
            require(abs(closure - RHO * value(b, "wallUTau") ** 2) <= 0.08 * closure,
                    f"{a['name']} {n}: developed gradient balances wall shear to"
                    f" {(closure / (RHO * value(b, 'wallUTau') ** 2) - 1) * 100:+.2f}%")
            require(abs(value(b, "dpPa") - s["dpb"][i][-1]) <= 1e-9 * abs(value(b, "dpPa")),
                    f"{a['name']} {n}: the drop the interface ended with is the drop"
                    " the block's own field measures")

    print(f"\nPASS: two live parallel edges - the residual hands over, one node drop"
          " closes against two rows with two different inertances, and the four"
          " attribution tables make a strictly ordered ladder of ramps"
          if not fails else "")
    for f in fails:
        print("FAILED: " + f)
    if fails:
        raise SystemExit(f"{len(fails)} both-parallel claim(s) failed")


if __name__ == "__main__":
    main()
