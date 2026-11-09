#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient-multi over two RANS detailed regions whose
proxy bores are deliberately NOT the regions' equivalent sections, and verify the
partitioned interface inside every time level for both sources of the region's
inertia.

The case is transientCoupledMultiTurbulent1d3d's topology with the one degree of
freedom it neutralised put back on the table: §5.27 measured partner nonlinearity
on two edges with `d = 2*sqrt(H*W/pi)` so that the bore arm and the region arm
installed the same inertance, and §5.26 measured attribution with two laminar
ducts. Here both bores are mis-stated - X1 by 1.273240, X2 by 1.696946 - and both
partners answer with a nonlinear law, so the question is whether the two accounts
stay separate. Ten families of statements:

  1. two arms, one knob   - the runs differ only in `inertias`, same seeds, both
                             partners k-omega SST
  2. what got installed   - each edge carries the inertance its arm names, every
                             step, and the block's own measurement is the same
                             pair of numbers in both arms
  3. the interface        - dpNetwork == R*Q exactly, residual == the gap between
                             the two solvers, under tolerance at every level
  4. the discrete 1D rows - each pipe's node drop is its loss at the NEW rate plus
                             its OWN arm's inertance times the rate change; the
                             bypass and leads close with published Churchill
  5. partner laws         - both blocks fit dp = c*q^1.68 across their rate travel
                             (a laminar partner gives 1.06), and X2's series
                             position leaves no exponent to read at all
  6. who pays for the schedule - the negotiated edge alone: worst index 0 on every
                             pass, the pinned edge constant and deep inside tol
  7. attribution moves mass, not resistance - the pinned edge is bit-identical
                             across arms while its first-step kick follows that
                             arm's M; the difference between arms is delta*M*Q/dt
  8. the ramp, per arm    - the restoring term is the sum of SLOPES for both arms
                             (chords read ~64% high in both), and the whole
                             trajectory follows the discrete inertial solve
  9. the two ramps differ by M alone - tau's ratio is the inertance ratio diluted
                             by the bypass, and the destination is the same
  10. each block's own field - wall functions, shear closure, conservation

usage: check_transient_coupled_multi_turbulent_inertia.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# fluid and geometry of tests/cases/transientCoupledMultiTurbulentInertia/case.dict
MU, RHO = 1.0e-3, 1000.0
NU = MU / RHO
Q_IN = 1.2e-4                    # m3/s switched on at Src at t = 0
P_REF = 1.0e6                    # Pa held at the sink
DT, N_STEPS = 0.05, 12           # time { dt 0.05; endTime 0.6; }
MAX_PASSES = 6                   # couple{ maxSubIterations }
TOL = 2e-4                       # couple{ tolerance }
NAME = ("X1", "X2")              # edge 0 is parallel, edge 1 is in series
ARMS = ("bore", "region")        # couple{ inertias ( network partner ) }
INERTIA = ("network", "partner")
EDGES = {"X1": (0.2, 0.005, 0.005), "X2": (0.12, 0.004, 0.008)}
BORE = {"X1": 0.005, "X2": 0.0049}      # the nominal circles, both too small
BY = (0.3, 0.008)
LEAD = (0.3, 0.02)
YPLUS_LOG_LAW = 11.81
M_BY = RHO * BY[0] / (0.25 * math.pi * BY[1] ** 2)
M_LEAD = RHO * LEAD[0] / (0.25 * math.pi * LEAD[1] ** 2)


def value(line, key):
    """First number stored under `key` in a log line ('q 1.5e-4 ...')."""
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
    bypass and leads are closed with, so the network side is judged by a formula
    from outside this repository."""
    re = max(re, 1e-12)
    a = (2.457 * math.log(1.0 / ((7.0 / re) ** 0.9 + 0.27 * roughness))) ** 16
    b = (37530.0 / re) ** 16
    return 8.0 * ((8.0 / re) ** 12 + 1.0 / (a + b) ** 1.5) ** (1.0 / 12.0)


def pipeChord(length, diameter, q):
    """dp/q of a 1D pipe carrying q, straight from Darcy-Weisbach."""
    area = 0.25 * math.pi * diameter ** 2
    v = q / area
    return (churchill(RHO * v * diameter / MU) * length / diameter
            * 0.5 * RHO * v * v / q)


def ductChord(name, q):
    """The resistance the block itself would publish for rate q is not known
    here; only its endpoints are read off the log. This is the analytic laminar
    reference the turbulent coefficient is compared against."""
    L, H, W = EDGES[name]
    return 12.0 * MU * L / (H ** 3 * W)


def popeBulkSpeed(uTau, dh):
    """Pope's smooth-wall friction law, 1/sqrt(f) = 2 log10(Re_b sqrt(f)) - 0.6,
    solved for the bulk speed a given wall shear carries."""
    u = uTau * 20.0
    for _ in range(200):
        f = 8.0 * (uTau / u) ** 2
        u = uTau * math.sqrt(8.0) * (2.0 * math.log10(u * dh / NU * math.sqrt(f)) - 0.6)
    return u


def parse(lines):
    """One dict per arm: accepted state per step, every interface pass, blocks."""
    arms = []
    cur = None
    for ln in lines:
        if (ln.startswith("transient start ") and ln.count(":") == 1
                and " edge " not in ln):
            t = ln.split()
            cur = dict(name=t[2].rstrip(":"), inertia=t[4], accelerate=t[6],
                       steps={}, passes=[], blocks=[], summary=None)
            arms.append(cur)
        elif cur is None:
            continue
        elif ln.startswith("transient step ") and " pass " in ln:
            cur["passes"].append(ln)
        elif ln.startswith("transient step ") and " edge " in ln:
            t = ln.split()
            cur["steps"].setdefault(int(t[2]), {})[int(t[4])] = ln
        elif ln.startswith("transient step ") and " newton " in ln:
            cur["steps"].setdefault(
                int(ln.split()[2].rstrip(":")), {})["s"] = ln
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
        raise SystemExit("transient multi turbulent-inertia solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    def check(cond, msg):
        if not cond:
            require(False, msg)

    arms = parse(lines)

    # ---- 1. two arms, one knob ----
    require([a["name"] for a in arms] == list(ARMS),
            f"both arms ran, in the order declared (got {[a['name'] for a in arms]})")
    for a in arms:
        require(sorted(a["steps"]) == list(range(1, N_STEPS + 1)),
                f"{a['name']}: all {N_STEPS} time levels reported an accepted state")
        require(all(all(i in a["steps"][n] for i in (0, 1))
                    for n in a["steps"]),
                f"{a['name']}: every step negotiated both edges")
    require([a["inertia"] for a in arms] == list(INERTIA),
            f"the two runs differ only in where the region's inertia comes from"
            f" ({[a['inertia'] for a in arms]})")
    require(all(a["accelerate"] == "no" for a in arms),
            "neither arm is accelerated: the schedule is the same operator twice,"
            " so an arm-to-arm difference can only be the inertance")
    for a in arms:
        for i, n in enumerate(NAME):
            e = a["steps"][1][i]
            check(word(e, "name") == n, f"{a['name']}: edge {i} is {n}")
        require(all(word(b, "turb") == "kOmegaSST" for b in a["blocks"]),
                f"{a['name']}: both partners faced the interface as RANS blocks")

    # ---- 2. what each arm installed ----
    M = [{i: [value(a["steps"][n][i], "M") for n in range(1, N_STEPS + 1)]
          for i in (0, 1)} for a in arms]
    inertSelf = [[value(b, "inertance") for b in a["blocks"]] for a in arms]
    for a_i, a in enumerate(arms):
        for i, n in enumerate(NAME):
            L, H, W = EDGES[n]
            section = H * W
            exact = RHO * L / section
            bore = RHO * L / (0.25 * math.pi * BORE[n] ** 2)
            want = bore if a["inertia"] == "network" else inertSelf[a_i][i]
            require(max(abs(x - want) for x in M[a_i][i]) <= 1e-9 * abs(want),
                    f"{a['name']} {n}: the inertance in its momentum row is"
                    f" {'the proxy bore rho*L/A' if a['inertia'] == 'network' else 'the one the block measured'}"
                    f" = {want:.10g} kg/m4, every step")
            require(abs(value(a["blocks"][i], "inertance") - exact)
                    <= 1e-9 * exact,
                    f"{a['name']} {n}: the block measures rho*V/A^2 = rho*L/(H*W) ="
                    f" {exact:.10g} kg/m4 on its own mesh"
                    f" ({value(a['blocks'][i], 'inertance'):.10g})")
    for i, n in enumerate(NAME):
        L, H, W = EDGES[n]
        ratio = (RHO * L / (0.25 * math.pi * BORE[n] ** 2)) / (RHO * L / (H * W))
        require(abs(M[0][i][0] / M[1][i][0] - ratio) <= 1e-9,
                f"{n}: the bore arm overstates this region's inertance by"
                f" {ratio:.6f} = the real section over the nominal circle"
                f" (H/W {W / H:.1f})")
        require(inertSelf[1][i] == inertSelf[0][i],
                f"{n}: the region's own inertance is the same number in both arms"
                f" ({inertSelf[0][i]:.10g}) - attribution is the arm's choice, not"
                " the block's")
    require(M[0][0][0] / M[1][0][0] != M[0][1][0] / M[1][1][0],
            f"the two edges are mis-stated by DIFFERENT amounts"
            f" ({M[0][0][0] / M[1][0][0]:.6f} vs {M[0][1][0] / M[1][1][0]:.6f}), so no"
            " single correction factor reproduces both arms")

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
        sw = {i: [int(value(st[n][i], "sweeps")) for n in range(1, N_STEPS + 1)]
              for i in (0, 1)}
        pas = [int(value(st[n]["s"], "passes")) for n in range(1, N_STEPS + 1)]
        nd = {k: [value(st[n]["nd"], k) for n in range(1, N_STEPS + 1)]
              for k in ("Src", "J1", "J2", "J3", "Snk")}
        pp = {k: [value(st[n]["pp"], k) for n in range(1, N_STEPS + 1)]
              for k in ("A", "X1", "By", "X2", "B")}

        w_id = max(abs(value(st[n][i], "dpNetwork") - R[i][n - 1] * q[i][n - 1])
                   for i in (0, 1) for n in range(1, N_STEPS + 1))
        require(w_id == 0.0,
                f"{a['name']}: the network still charges R times the rate, RANS partners"
                f" and a mis-stated bore both included (worst {w_id:.1e} Pa)")
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
        require(all(abs(value(st[n]["s"], "t") - n * DT) <= 1e-12
                    for n in range(1, N_STEPS + 1)),
                f"{a['name']}: step n is at t = n*dt")
        require(all(value(st[n][i], "blockConverged") == 1
                    for i in (0, 1) for n in range(1, N_STEPS + 1)),
                f"{a['name']}: every block solve that reached the interface converged")
        require(all(abs(q[i][n - 1] - pp[NAME[i]][n - 1]) <= 1e-15 * Q_IN
                    for i in (0, 1) for n in range(1, N_STEPS + 1))
                and all(abs(pp["A"][n - 1] - Q_IN) <= 1e-15 * Q_IN
                        and abs(pp["A"][n - 1] - pp["X1"][n - 1] - pp["By"][n - 1])
                        <= 1e-12 * Q_IN for n in range(1, N_STEPS + 1)),
                f"{a['name']}: the rates the interface negotiated are the pipes' own, and"
                " the branch split adds up to the lead")
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
                f"{a['name']}: both regions' node drops close with THIS arm's"
                f" inertance ({rX1:.1e} parallel, {rX2:.1e} series)")
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
                f"{a['name']}: the bypass and the two leads close against published"
                f" Churchill chords plus their own inertance ({rBy:.1e}, {rLead:.1e})")

        series.append(dict(q=q, R=R, dp=dp, res=res, sw=sw, pas=pas, pp=pp, nd=nd,
                           dp12=dp12, dp23=dp23, rX1=rX1, rX2=rX2))

    # ---- 5. the partners' laws ----
    for a_i, a in enumerate(arms):
        s = series[a_i]
        q, dp, R, pp = s["q"], s["dp"], s["R"], s["pp"]
        mX = math.log(dp[0][-1] / dp[0][0]) / math.log(q[0][-1] / q[0][0])
        mBy = (1.0 + math.log(pipeChord(*BY, pp["By"][-1])
                              / pipeChord(*BY, pp["By"][0]))
               / math.log(pp["By"][-1] / pp["By"][0]))
        s["mX"], s["mBy"] = mX, mBy
        travel = abs(q[0][-1] - q[0][0]) / q[0][-1]
        require(1.5 < mX < 2.0,
                f"{a['name']}: the parallel block's own answers fit dp = c*q^{mX:.4f}"
                f" across a {travel * 100:.1f}% rate move - a laminar partner measured"
                " the same way gives 1.06")
        require(1.6 < mBy < 1.85,
                f"{a['name']}: the 1D bypass is turbulent too, exponent {mBy:.4f} at Re"
                f" {RHO * pp['By'][-1] / (0.25 * math.pi * BY[1] ** 2) * BY[1] / MU:.0f}")
        s["flatX2"] = (max(abs(x - Q_IN) for x in q[1]) / Q_IN,
                       max(abs(x - dp[1][0]) for x in dp[1]) / dp[1][0],
                       max(abs(x - R[1][0]) for x in R[1]) / R[1][0])
        require(s["flatX2"][0] <= 1e-14 and s["flatX2"][1] <= 1e-6,
                f"{a['name']}: the series edge has no exponent to read - the network"
                f" hands it Q_IN to {s['flatX2'][0]:.1e} and its drop moves only"
                f" {s['flatX2'][1]:.1e} over the whole ramp")
        s["lam"] = [R[i][-1] / ductChord(NAME[i], 1.0) for i in (0, 1)]
        require(all(x >= 4.0 for x in s["lam"]) and abs(s["lam"][0] - s["lam"][1]) >= 0.2,
                f"{a['name']}: the settled coefficients are {s['lam'][0]:.4f}x and"
                f" {s['lam'][1]:.4f} the planar-Poiseuille resistance of the same"
                " ducts - turbulence, and two different problems")

    # ---- 6. who pays for the schedule ----
    for a_i, a in enumerate(arms):
        s = series[a_i]
        idx = [int(value(p, "worst")) for p in a["passes"]]
        r0 = [float(p.split()[8]) for p in a["passes"]]
        r1 = [float(p.split()[9]) for p in a["passes"]]
        require(set(idx) == {0},
                f"{a['name']}: the schedule stops on the negotiated edge in all"
                f" {len(idx)} passes - the pinned one never holds the interface up")
        require(sum(1 for x, y in zip(r0, r1) if x > TOL >= y) >= 10,
                f"{a['name']}: in {sum(1 for x, y in zip(r0, r1) if x > TOL >= y)}"
                f" of the {len(r0)} passes the series edge was already inside tolerance"
                " and was refreshed anyway - a multi-edge schedule is one number for all")
        require(all(1 <= p <= MAX_PASSES for p in s["pas"]),
                f"{a['name']}: every time level agrees inside the allowed"
                f" sub-iterations (max {max(s['pas'])} of {MAX_PASSES})")
        require(sum(s["pas"]) - N_STEPS == len(a["passes"]) - N_STEPS,
                f"{a['name']}: one pass is one network solve - {sum(s['pas'])} accepted"
                f" of {len(a['passes'])} run, {len(a['passes']) - N_STEPS} over tolerance")
        require(all(x == 6 for x in s["sw"][1])
                and min(s["sw"][0]) >= 10 and max(s["sw"][0]) <= 14,
                f"{a['name']}: the pinned block needs {s['sw'][1][0]} SIMPLE sweeps every"
                f" step while the negotiated one takes {min(s['sw'][0])}-{max(s['sw'][0])}")
        require(max(s["res"][1]) <= TOL / 100,
                f"{a['name']}: the series edge sits {TOL / max(s['res'][1]):.0f}x inside"
                " the tolerance the schedule is chasing")

    # ---- 7. attribution moves mass, not resistance ----
    a, b = series
    dX2 = (abs(a["q"][1][-1] - b["q"][1][-1]),
           abs(a["R"][1][-1] / b["R"][1][-1] - 1.0),
           abs(a["dp"][1][-1] / b["dp"][1][-1] - 1.0))
    require(dX2[0] <= 1e-19 and dX2[1] <= 1e-9,
            f"the pinned edge ends on the same state in both arms - rate to one unit in"
            f" the last place ({dX2[0]:.1e} m3/s), coefficient to {dX2[1]:.1e} - because"
            " its rate was never a question")
    err = [abs(s["dp23"][0] - (s["R"][1][0] * Q_IN + M[i][1][0] * Q_IN / DT))
           / s["dp23"][0] for i, s in ((0, a), (1, b))]
    require(all(e <= 1e-12 for e in err),
            f"each arm's first-step series drop is its friction at Q_IN plus its OWN"
            f" inertance's kick ({a['dp23'][0]:.6f} = {a['R'][1][0] * Q_IN:.4f}"
            f" + {M[0][1][0] * Q_IN / DT:.4f}; {b['dp23'][0]:.6f} ="
            f" {M[1][1][0] * Q_IN / DT:.4f} inertial), worst {max(err):.1e}")
    dKick = (M[0][1][0] - M[1][1][0]) * Q_IN / DT
    require(abs((a["dp23"][0] - b["dp23"][0]) - dKick) <= 1e-6 * abs(dKick),
            f"the whole difference between the arms on that edge is the inertance they"
            f" disagree about: {a['dp23'][0] - b['dp23'][0]:.6f} Pa against"
            f" dM*Q/dt = {dKick:.6f} Pa")
    require(a["res"][1][0] <= TOL and b["res"][1][0] <= TOL,
            "both arms still agreed the pinned edge on the first pass of the first"
            " level, kick and all")

    # ---- 8. the ramp, per arm ----
    for a_i, arm in enumerate(arms):
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
        MX1 = M[a_i][0][0]
        Mtot = MX1 + M_BY
        RX, RBy = R[0][-1], pipeChord(*BY, pp["By"][-1])
        slope = s["mX"] * RX + s["mBy"] * RBy
        chord = RX + RBy
        tauS = DT / math.log(1.0 + DT * slope / Mtot)
        tauC = DT / math.log(1.0 + DT * chord / Mtot)
        Bpred = 1.0 / (1.0 + DT * slope / Mtot)
        s["tauS"], s["tauC"], s["mBar"] = tauS, tauC, slope / chord
        require(0.5 < Bc < 0.95, f"{arm['name']}: one time constant spans"
                                 f" {1.0 / (1.0 - Bc):.1f} steps")
        require(abs(tauS / tauFit - 1.0) <= 0.01,
                f"{arm['name']}: the ramp's time constant is the partners' SLOPES,"
                f" (M_X1+M_By)/(m_X R_X + m_By R_By) = {tauS:.6f} s against the fitted"
                f" {tauFit:.6f} s ({(tauS / tauFit - 1) * 100:+.3f}%)"
                f" - with this arm's inertance {MX1:.6g} kg/m4")
        require(abs(Bc / Bpred - 1.0) <= 0.005,
                f"{arm['name']}: the same closed form per step, B {Bpred:.6f} against"
                f" the measured {Bc:.6f} ({(Bc / Bpred - 1) * 100:+.3f}%)")
        require(tauC / tauFit - 1.0 >= 0.5,
                f"{arm['name']}: the CHORD formula reads {tauC:.6f} s here,"
                f" {(tauC / tauFit - 1) * 100:.1f}% off the fitted ramp")
        require(abs((tauC / tauS) / s["mBar"] - 1.0) <= 0.06,
                f"{arm['name']}: the miss is the resistance-weighted exponent -"
                f" chord/slope {tauC / tauS:.4f} against m = {s['mBar']:.4f}")
        q1 = Q_IN * (RBy + M_BY / DT) / (Mtot / DT + R[0][0] + pipeChord(*BY, pp["By"][0]))
        dev1 = abs(q[0][0] - q1) / q1
        require(dev1 <= 5e-3,
                f"{arm['name']}: the first step splits the sudden rate by this arm's"
                f" inertance ({q[0][0]:.6g} against {q1:.6g} m3/s, {dev1:.1e})")
        resistive = Q_IN * RBy / (RBy + RX)
        require(abs(q[0][0] / resistive - 1.0) >= 0.02,
                f"{arm['name']}: that split is {abs(q[0][0] / resistive - 1) * 100:.1f}%"
                " away from the resistive one it relaxes towards")
        require(abs(qInf / resistive - 1.0) <= 0.02,
                f"{arm['name']}: the fitted asymptote {qInf:.6g} m3/s is the resistive"
                f" split of the settled chords, {resistive:.6g}"
                f" ({(qInf / resistive - 1) * 100:+.3f}%)")
        dev = max(abs(q[0][n] - (qInf + (q[0][0] - qInf) * Bpred ** n)) / qInf
                  for n in range(N_STEPS))
        require(dev <= 5e-3,
                f"{arm['name']}: the discrete slope-form ramp reproduces every step of"
                f" the trajectory to {dev * 100:.4f}% of the operating point")
        gap = (q[0][-1] - qInf) / qInf
        tail = Bpred ** N_STEPS * (q[0][0] - qInf) / q[0][-1]
        require(abs(gap) <= 0.02 and abs(abs(gap) / abs(tail) - 1.0) <= 0.25,
                f"{arm['name']}: the series ends {gap * 100:+.4f}% short of its"
                f" asymptote, against the {tail * 100:+.4f}% its own decay leaves"
                f" unburnt after {N_STEPS} steps ({N_STEPS * DT / tauFit:.2f} constants)")

    # ---- 9. the two ramps differ by M alone ----
    require(abs(a["R"][1][-1] / b["R"][1][-1] - 1.0) <= 1e-9
            and max(abs(a["R"][0][n] - b["R"][0][n]) / a["R"][0][n]
                    for n in range(N_STEPS)) <= 0.2,
            f"attribution does not move the interface's answer: the pinned coefficient"
            f" agrees to {abs(a['R'][1][-1] / b['R'][1][-1] - 1.0):.1e}, and the"
            f" negotiated one only follows the rate its arm gave it"
            f" (worst {max(abs(a['R'][0][n] - b['R'][0][n]) / a['R'][0][n] for n in range(N_STEPS)) * 100:.2f}%)")
    tauRatio = a["tau"] / b["tau"]
    mRatio = (M[0][0][0] + M_BY) / (M[1][0][0] + M_BY)
    require(abs(tauRatio / mRatio - 1.0) <= 0.02,
            f"the ramps' ratio IS the inertance ratio diluted by the bypass,"
            f" {tauRatio:.5f} against {mRatio:.5f} ({(tauRatio / mRatio - 1) * 100:+.2f}%)"
            f" - tau {a['tau']:.6f} s (bore) vs {b['tau']:.6f} s (region), while the"
            f" resistive destination differs by only"
            f" {abs(a['qInf'] / b['qInf'] - 1) * 100:.4f}%")
    require(abs(a["qInf"] / b["qInf"] - 1.0) <= 0.01
            and abs(a["q"][0][-1] / b["q"][0][-1] - 1.0) >= 0.005,
            f"the destination is the same and the arms are not there yet:"
            f" fitted asymptotes {a['qInf']:.6g} / {b['qInf']:.6g} m3/s"
            f" ({abs(a['qInf'] / b['qInf'] - 1) * 100:.4f}% apart) but the rates after"
            f" {N_STEPS} steps are {a['q'][0][-1]:.6g} / {b['q'][0][-1]:.6g}"
            f" ({abs(a['q'][0][-1] / b['q'][0][-1] - 1) * 100:.3f}% apart)")
    require(abs((a["q"][0][-1] - a["qInf"]) / a["qInf"])
            > abs((b["q"][0][-1] - b["qInf"]) / b["qInf"]),
            "and the arm that overstates the inertance is the one still further from"
            f" that shared destination ({(a['q'][0][-1] / a['qInf'] - 1) * 100:+.4f}% vs"
            f" {(b['q'][0][-1] / b['qInf'] - 1) * 100:+.4f}%) - the heavier region"
            " takes longer to reach the rate the resistances agree on")
    require(abs(a["mX"] / b["mX"] - 1.0) <= 0.01,
            f"the partner's exponent is the same law read from either arm"
            f" ({a['mX']:.4f} vs {b['mX']:.4f}) - inertance does not touch the slope")

    # ---- 10. each block's own field ----
    for a_i, arm in enumerate(arms):
        s = series[a_i]
        for i, n in enumerate(NAME):
            fin = arm["blocks"][i]
            L, H, W = EDGES[n]
            area = H * W
            f = value(fin, "fluxIn")
            plug = value(fin, "plugSpeed")
            require(abs(f / s["q"][i][-1] - 1.0) <= 1e-4,
                    f"{arm['name']} {n}: the block carries the rate its edge ended with"
                    f" ({f:.6g} m3/s)")
            require(abs(value(fin, "inletArea") - area) <= 1e-12 * area,
                    f"{arm['name']} {n}: the inlet patch is the section the inertance was"
                    f" referred to ({value(fin, 'inletArea'):.6g} m2)")
            require(abs(value(fin, "fluxOut") - f) <= 1e-9 * f
                    and value(fin, "planeFluxMax") - value(fin, "planeFluxMin")
                    <= 1e-9 * f,
                    f"{arm['name']} {n}: what enters the region leaves it, and every"
                    " cross-section plane inside carries the same rate")
            require(RHO * plug * (2.0 * H) / MU > 4000.0,
                    f"{arm['name']} {n}: turbulent at Re_bulk"
                    f" {RHO * plug * 2.0 * H / MU:.0f}")
            require(value(fin, "nutMax") / NU > 20.0,
                    f"{arm['name']} {n}: eddy viscosity is {value(fin, 'nutMax') / NU:.1f}"
                    " times the molecular one")
            yp = value(fin, "wallYPlusMax")
            require(YPLUS_LOG_LAW <= yp < 300.0,
                    f"{arm['name']} {n}: wall cells at y+ {yp:.1f}, inside the log-law"
                    " band the wall functions assume")
            grad = abs(value(fin, "probeSlope"))
            closure = grad * H / 2.0
            uTau = value(fin, "wallUTau")
            require(abs(closure - RHO * uTau ** 2) <= 0.08 * closure,
                    f"{arm['name']} {n}: developed gradient balances wall shear,"
                    f" |dp/dx|*H/2 {closure:.6g} against rho*u_tau^2"
                    f" {RHO * uTau ** 2:.6g} Pa ({(closure / (RHO * uTau ** 2) - 1) * 100:+.2f}%)")
            require(0.94 <= value(fin, "dpPa") / (grad * L) <= 1.06,
                    f"{arm['name']} {n}: the whole block's drop is"
                    f" {value(fin, 'dpPa') / (grad * L):.4f} of its developed gradient"
                    " over its length")
            ubulk = popeBulkSpeed(uTau, 2.0 * H)
            require(abs(ubulk / plug - 1.0) <= 0.08,
                    f"{arm['name']} {n}: Pope's friction law from the block's own u_tau"
                    f" gives {ubulk:.6g} m3/s-area-averaged against the measured"
                    f" {plug:.6g} ({(ubulk / plug - 1) * 100:+.1f}%)")
            require(1.0 < value(fin, "uxMax") / plug <= 1.25,
                    f"{arm['name']} {n}: a flat turbulent profile,"
                    f" u_max/ubar = {value(fin, 'uxMax') / plug:.3f}")
            require(abs(value(fin, "dpPa") - s["dp"][i][-1]) <= 1e-6 * s["dp"][i][-1],
                    f"{arm['name']} {n}: the drop the interface ended with is the drop"
                    " the block's own field measures")

    nPass = len(arms[0]["passes"]) + len(arms[1]["passes"])
    nAcc = sum(sum(s["pas"]) for s in series)
    print(f"cost: {nPass * 2} block solves and {nPass} interface passes for"
          f" {2 * N_STEPS} agreed time levels on two RANS edges,"
          f" {nAcc} of them accepted")

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: a nonlinear partner and a mis-stated region inertance stay two"
          " separate accounts inside every time level")


if __name__ == "__main__":
    sys.exit(main())
