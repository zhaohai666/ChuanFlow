#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient against a RANS detailed region and verify
the partitioned interface inside every time level.

The case is the steady turbulent coupled network put in motion: Src -> A ->
{ X (3D k-omega SST channel) || By (1D bypass) } -> B -> Snk, at rest until an
injected rate is switched on at t = 0. Each time level agrees the interface
again, so the series shows what a non-linear partner does to the network's own
dynamics. Seven families of statements:

  1. one state per step     - the block carries the rate the step delivers, and
                              the residual is the gap between two solvers
  2. the discrete 1D law    - every pipe's node drop is its loss at the NEW rate
                              plus its inertance times the rate change; for the
                              bypass the loss is the published Churchill chord,
                              so the 1D side closes without the code under test
  3. partner nonlinearity   - the block's own dp(q) answers across the ramp fit
                              dp = c q^1.68, and the bypass's own law gives 1.73;
                              a laminar partner measured the same way gives 1.06
  4. the ramp               - the sudden rate first divides by INERTANCE, then
                              relaxes with the sum of the partners' SLOPES,
                              m_i R_i, as the closed form's restoring term; the
                              chord formula R_X + R_By - exact for the laminar
                              case, where both m_i are 1 - misreads this ramp
  5. the interface schedule - every time level agrees within a couple of partner
                              runs, and each of those is exactly one block solve
  6. two seeds, one answer  - starting off the fixed point and next to it land on
                              the same state and the same time constant
  7. the settled block      - the region ends in a wall-function-resolved
                              turbulent channel whose shear balances its own
                              developed gradient

usage: check_transient_coupled_turbulent.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# fluid and geometry of tests/cases/transientCoupledTurbulent/case.dict
MU, RHO = 1.0e-3, 1000.0
NU = MU / RHO
H, W, LB = 0.005, 0.005, 0.2     # block gap (wall-normal), span, length
AREA = H * W                     # 2.5e-5 m2, the region's own cross-section
DH = 2.0 * H                     # plane channel: the two walls carry the shear

# the 1D side: X is the coupled proxy pipe, By the bypass, A and B the leads
DXPROXY, LX = 0.005641895835, 0.2
DBY, LBY = 0.008, 0.3
DLEAD, LLEAD = 0.02, 0.3
SEC_X = 0.25 * math.pi * DXPROXY ** 2
SEC_BY = 0.25 * math.pi * DBY ** 2
SEC_LEAD = 0.25 * math.pi * DLEAD ** 2
M_X = RHO * LX / SEC_X           # the region's inertance, carried by the proxy
M_BY = RHO * LBY / SEC_BY
M_LEAD = RHO * LLEAD / SEC_LEAD

Q_IN = 1.2e-4                    # m3/s switched on at Src at t = 0
P_REF = 1.0e6                    # Pa held at the sink
DT, N_STEPS = 0.06, 27           # time { dt 0.06; endTime 1.62; }
MAX_PASSES = 6                   # couple{ maxSubIterations }
TOL_INTERFACE = 2e-4             # couple{ tolerance } the case asks for
SEEDS = ("0", "1")               # starting coefficients 1.0e7 and 2.2e7
YPLUS_LOG_LAW = 11.81            # the wall functions are valid above this


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


def churchillFriction(re, roughness=0.0):
    """Churchill's 1977 Darcy friction factor, all three flow regimes. This is
    the published correlation, written from the literature: it is what the 1D
    bypass's own momentum row is closed with, so the network's friction law is
    checked against an outside formula rather than against itself."""
    re = max(re, 1e-12)
    a = (2.457 * math.log(1.0 / ((7.0 / re) ** 0.9 + 0.27 * roughness))) ** 16
    b = (37530.0 / re) ** 16
    return 8.0 * ((8.0 / re) ** 12 + 1.0 / (a + b) ** 1.5) ** (1.0 / 12.0)


def pipeChord(length, diameter, q):
    """dp/q of a 1D pipe carrying q, straight from Darcy-Weisbach."""
    area = 0.25 * math.pi * diameter ** 2
    v = q / area
    re = RHO * v * diameter / MU
    dp = churchillFriction(re) * length / diameter * 0.5 * RHO * v * v
    return dp / q


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
        raise SystemExit("transient turbulent coupled solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    # the same statement for all 27 time levels: a passing step has nothing to say
    def check(cond, msg):
        if not cond:
            require(False, msg)

    # ---- collect one entry per time level, per seed ----
    runs = []
    for ln in lines:
        if ln.startswith("transient step "):
            n = int(ln.split()[2].rstrip(":"))
            if not runs or n in runs[-1]:     # the same step number starts again
                runs.append({})              # only when the next seed begins
            runs[-1][n] = [ln, None, None]
        elif ln.startswith("nodes step "):
            runs[-1][int(ln.split()[2].rstrip(":"))][1] = ln
        elif ln.startswith("pipes step "):
            runs[-1][int(ln.split()[2].rstrip(":"))][2] = ln

    blocks, summaries = {}, {}
    for ln in lines:
        if ln.startswith("transient block "):
            blocks[ln.split()[2].rstrip(":")] = ln
        elif ln.startswith("transient coupled seed "):
            summaries[ln.split()[3].rstrip(":")] = ln

    require(len(runs) == len(SEEDS),
            f"both starting coefficients ran their own series (got {len(runs)})")
    require(all(sorted(s) == list(range(1, N_STEPS + 1)) for s in runs),
            f"every one of the {N_STEPS} steps the case asks for reported a state")
    require(all(all(p is not None for p in v.values()) for v in runs),
            "every step reports its whole 1D state")
    require(set(blocks) == set(SEEDS) and set(summaries) == set(SEEDS),
            "each seed reported its block and its interface summary")

    series = []
    for seed, steps in zip(SEEDS, runs):
        q, qBy, qA, qB, dpB, R, pas, sw = ([] for _ in range(8))
        dp12, pSrc, pJ1, pJ2, pSnk = ([] for _ in range(5))
        for n in sorted(steps):
            s, nd, pp = steps[n]
            t = value(s, "t")
            check(abs(t - n * DT) <= 1e-12, f"seed {seed} step {n}: t is {n} steps")
            check(value(s, "blockConverged") == 1,
                  f"seed {seed} step {n}: the block solve converged")
            check(value(s, "newton") < 1e-11,
                  f"seed {seed} step {n}: the 1D step solved its own system")
            x = value(s, "q")
            check(abs(x - value(pp, "X")) <= 1e-12 * x,
                  f"seed {seed} step {n}: the reported rate is the coupled pipe's")
            check(abs(value(s, "dpNetwork") - value(s, "R") * x)
                  <= 1e-12 * value(s, "dpNetwork"),
                  f"seed {seed} step {n}: the network charges R times the rate")
            gap = abs(value(s, "dpNetwork") - value(s, "dpBlock")) / value(s, "dpBlock")
            check(abs(gap - value(s, "residual")) <= 1e-9,
                  f"seed {seed} step {n}: residual is the cross-solver gap")
            check(value(s, "residual") <= TOL_INTERFACE,
                  f"seed {seed} step {n}: interface residual "
                  f"{value(s, 'residual'):.2e} under the case tolerance")
            check(abs(value(pp, "A") - Q_IN) <= 1e-15 * Q_IN,
                  f"seed {seed} step {n}: the lead carries the whole injected rate")
            check(abs(value(pp, "A") - x - value(pp, "By")) <= 1e-12 * Q_IN,
                  f"seed {seed} step {n}: the branch split adds up to the lead")

            q.append(x); qBy.append(value(pp, "By")); qA.append(value(pp, "A"))
            qB.append(value(pp, "B"))
            dpB.append(value(s, "dpBlock"))
            R.append(value(s, "R"))
            pas.append(int(value(s, "passes"))); sw.append(int(value(s, "sweeps")))
            pSrc.append(value(nd, "Src")); pJ1.append(value(nd, "J1"))
            pJ2.append(value(nd, "J2")); pSnk.append(value(nd, "Snk"))
            dp12.append(pJ1[-1] - pJ2[-1])

        # ---- 2. the discrete momentum rows the step solved ----
        #     the node drop is the loss at the NEW rate plus inertance times the
        #     rate change; the network starts at rest, so the first step's
        #     previous rate is zero and its inertial kick is the whole drop
        worstX = worstBy = worstLead = 0.0
        for i in range(N_STEPS):
            prev = lambda a: a[i - 1] if i else 0.0
            rowX = abs(dp12[i] - (R[i] * q[i] + M_X * (q[i] - prev(q)) / DT)) / dp12[i]
            chord = pipeChord(LBY, DBY, qBy[i])
            rowBy = abs(dp12[i] - (chord * qBy[i] + M_BY * (qBy[i] - prev(qBy)) / DT)) / dp12[i]
            worstX = max(worstX, rowX); worstBy = max(worstBy, rowBy)
            for drop, flow, prevFlow, length, bore in (
                    (pSrc[i] - pJ1[i], qA[i], prev(qA), LLEAD, DLEAD),
                    (pJ2[i] - pSnk[i], qB[i], prev(qB), LLEAD, DLEAD)):
                worstLead = max(worstLead, abs(drop - (pipeChord(length, bore, flow) * flow
                                                       + M_LEAD * (flow - prevFlow) / DT))
                                / abs(drop))
            check(abs(pSnk[i] - P_REF) <= 1e-9,
                  f"seed {seed} step {i+1}: the sink holds its reference pressure")

        series.append(dict(q=q, qBy=qBy, dp=dpB, R=R, pas=pas, sw=sw,
                           worstX=worstX, worstBy=worstBy, worstLead=worstLead))

    # one number per family, so 27 steps x 12 statements do not bury the result
    require(max(s["worstX"] for s in series) <= 1e-9,
            f"the coupled region's node drop is R*Q + M*dQ/dt at every step of both"
            f" seeds (worst {max(s['worstX'] for s in series):.2e})")
    require(max(s["worstBy"] for s in series) <= 1e-12,
            f"the bypass's row closes with Churchill's own chord plus its inertance"
            f" (worst {max(s['worstBy'] for s in series):.2e})")
    require(max(s["worstLead"] for s in series) <= 1e-9,
            f"the two leads close the same way, inertial kick on the first step"
            f" included (worst {max(s['worstLead'] for s in series):.2e})")

    # ---- 3. the two partners' laws, read off the ramp they produced ----
    for seed, s in zip(SEEDS, series):
        q, dp, qBy = s["q"], s["dp"], s["qBy"]
        mX = math.log(dp[-1] / dp[0]) / math.log(q[-1] / q[0])
        mBy = 1.0 + math.log(pipeChord(LBY, DBY, qBy[-1]) / pipeChord(LBY, DBY, qBy[0])) \
                  / math.log(qBy[-1] / qBy[0])
        s["mX"], s["mBy"] = mX, mBy
        spread = (q[0] - q[-1]) / q[-1]
        require(1.5 < mX < 2.0,
                f"seed {seed}: the block's own answers across a {spread*100:.1f}% rate"
                f" move fit dp = c*q^{mX:.3f} - a laminar partner measured the same"
                " way gives 1.06")
        require(1.6 < mBy < 1.85,
                f"seed {seed}: the 1D bypass is turbulent too, q^{mBy:.3f} at Re"
                f" {RHO * qBy[-1] / SEC_BY * DBY / MU:.0f}")
        require(all(dp[i] > dp[i + 1] for i in range(N_STEPS - 1))
                and all(q[i] > q[i + 1] for i in range(N_STEPS - 1)),
                f"seed {seed}: the detailed branch relaxes downward in rate and drop"
                " together, as a monotone partner law requires")

    # ---- 4. the ramp: slopes, not chords ----
    for seed, s in zip(SEEDS, series):
        q, R = s["q"], s["R"]
        xs, ys = q[1:-1], q[2:]
        k = float(len(xs))
        sx, sy = sum(xs), sum(ys)
        B = (k * math.fsum(a * b for a, b in zip(xs, ys)) - sx * sy) / (
            k * math.fsum(a * a for a in xs) - sx * sx)
        Acoef = (sy - B * sx) / k
        tauFit = -DT / math.log(B)
        qInf = Acoef / (1.0 - B)
        s["B"], s["tau"], s["qInf"] = B, tauFit, qInf

        RX, RBy = R[-1], pipeChord(LBY, DBY, s["qBy"][-1])
        Mtot = M_X + M_BY
        slope = s["mX"] * RX + s["mBy"] * RBy
        chord = RX + RBy
        predict = lambda restore: DT / math.log(1.0 + DT * restore / Mtot)
        tauSlope, tauChord = predict(slope), predict(chord)
        Bslope = 1.0 / (1.0 + DT * slope / Mtot)
        mBar = slope / chord
        s["tauSlope"], s["tauChord"], s["mBar"] = tauSlope, tauChord, mBar
        s["RX"], s["RBy"] = RX, RBy

        # the first step: the network was at rest, so the rate divides by inertance
        q1 = Q_IN * (RBy + M_BY / DT) / (Mtot / DT + R[0] + RBy)
        dev1 = abs(q[0] - q1) / q1
        require(dev1 <= 5e-3,
                f"seed {seed}: the first step splits the sudden rate by inertance"
                f" ({q[0]:.6g} vs {q1:.6g} m3/s, {dev1:.1e})")
        require(abs(q[0] / (RBy / chord * Q_IN) - 1.0) >= 0.02,
                f"seed {seed}: the inertial split puts the detailed branch"
                f" {abs(q[0] / (RBy / chord * Q_IN) - 1.0) * 100:.1f}% away from the"
                " resistive split it relaxes to, not on it")

        require(0.5 < B < 0.95,
                f"seed {seed}: one time constant spans {1.0 / (1.0 - B):.1f} steps")
        require(abs(tauSlope / tauFit - 1.0) <= 0.01,
                f"seed {seed}: the ramp's time constant is the partners' SLOPES,"
                f" (M_X+M_By)/(m_X R_X + m_By R_By) = {tauSlope:.4f} s against the"
                f" fitted {tauFit:.4f} s ({(tauSlope / tauFit - 1) * 100:+.2f}%)")
        require(abs(B / Bslope - 1.0) <= 0.005,
                f"seed {seed}: the same closed form per step, B {Bslope:.4f} against"
                f" the measured {B:.4f}")
        require(tauChord / tauFit - 1.0 >= 0.5,
                f"seed {seed}: the CHORD formula that the laminar case satisfies to"
                f" 1.6% reads {tauChord:.4f} s here, {(tauChord / tauFit - 1) * 100:.0f}% off")
        require(abs((tauChord / tauSlope) / mBar - 1.0) <= 0.06,
                f"seed {seed}: the miss is the resistance-weighted exponent -"
                f" chord/slope {tauChord / tauSlope:.3f} against m = {mBar:.3f},"
                " the rest being the implicit step's own ln(1+x)")
        gap = (q[-1] - qInf) / qInf
        tail = math.exp(-N_STEPS * DT / tauFit) * (q[0] - qInf) / q[-1]
        require(0.0 < gap <= 4e-4,
                f"seed {seed}: the series ends {gap*100:.4f}% short of its asymptote,"
                f" the {tail*100:.4f}% exponential tail its own time constant predicts")
        require(abs(qInf - RBy / chord * Q_IN) <= 0.005 * Q_IN,
                f"seed {seed}: the fitted asymptote {qInf:.6g} is the resistive split"
                f" of the two settled chords, {RBy / chord * Q_IN:.6g} m3/s")

    # ---- 5. the interface schedule, and what it cost ----
    allPas = [p for s in series for p in s["pas"]]
    require(all(1 <= p <= MAX_PASSES for p in allPas),
            "every time level agrees inside the allowed sub-iterations")
    require(all(s["pas"][0] == 3 for s in series),
            "both starting coefficients are corrected within three partner runs")
    require(max(allPas) == 3,
            f"no step ever needed more than {max(allPas)} partner runs,"
            f" the accelerated interface holding the ramp to {sum(allPas)/len(allPas):.2f} on average")
    require(sum(allPas) <= 2.0 * len(allPas),
            f"the two series cost {sum(allPas)} partner runs for {2 * N_STEPS} steps")
    require(min(s["sw"][0] for s in series) >= 8,
            "the first step's block solve is a cold field and says so in its sweeps")
    solves = sum(1 for ln in lines if ln.startswith("SIMPLE iter 1:"))
    require(solves == sum(allPas),
            f"one pass is exactly one block solve: {solves} SIMPLE starts for"
            f" {sum(allPas)} interface passes")

    # ---- 6. two seeds, one answer ----
    a, b = series
    dq = abs(a["q"][-1] - b["q"][-1]) / a["q"][-1]
    dR = abs(a["RX"] - b["RX"]) / a["RX"]
    dtau = abs(a["tau"] / b["tau"] - 1.0)
    dm = abs(a["mX"] / b["mX"] - 1.0)
    require(dq <= 1e-4 and dR <= 1e-4,
            f"seeds 1.0e7 and 2.2e7 end on one state: rate to {dq:.1e},"
            f" coefficient to {dR:.1e}")
    require(dtau <= 5e-3,
            f"the same ramp time constant comes out of both seeds"
            f" ({a['tau']:.4f} s, {b['tau']:.4f} s, {dtau:.1e})")
    require(dm <= 0.01,
            f"and the same partner exponent ({a['mX']:.3f} vs {b['mX']:.3f})")

    # ---- 7. the settled region, from its own fields ----
    for seed, s in zip(SEEDS, series):
        fin = blocks[seed]
        su = summaries[seed]
        q, plug = value(fin, "fluxIn"), value(fin, "plugSpeed")
        dp, slope = value(fin, "dpPa"), abs(value(fin, "probeSlope"))
        uTau, yp = value(fin, "wallUTau"), value(fin, "wallYPlusMax")
        nut = value(fin, "nutMax")
        require(word(fin, "turb") == "kOmegaSST",
                f"seed {seed}: the partner that faced the interface was the RANS block")
        require(abs(q - s["q"][-1]) <= 1e-7 * q,
                f"seed {seed}: the block carries the rate the series settled on")
        require(abs(dp - value(su, "dp")) <= 1e-12 * dp
                and abs(value(su, "R") - s["RX"]) <= 1e-12 * s["RX"]
                and int(value(su, "steps")) == N_STEPS,
                f"seed {seed}: the summary reports the last step of its series")
        require(abs(value(fin, "inletArea") * plug - q) <= 1e-12 * q,
                f"seed {seed}: plug speed times inlet area is the carried rate")
        require(abs(value(fin, "fluxOut") - q) <= 1e-6 * q,
                f"seed {seed}: what enters the region leaves it")
        require(value(fin, "planeFluxMax") - value(fin, "planeFluxMin") <= 1e-6 * q,
                f"seed {seed}: every cross-section plane inside carries the same rate")
        require(RHO * plug * DH / MU > 4000.0,
                f"seed {seed}: the region is turbulent at Re_bulk"
                f" {RHO * plug * DH / MU:.0f}, which is what the partner model is")
        require(yp >= YPLUS_LOG_LAW and yp < 300.0,
                f"seed {seed}: the wall cells sit at y+ {yp:.1f}, inside the log-law"
                " band the wall functions assume (first cell centre y+ ="
                f" {(H / 16.0) * uTau / NU:.1f})")
        closure = slope * H / 2.0
        require(abs(closure - RHO * uTau ** 2) <= 0.01 * closure,
                f"seed {seed}: the developed gradient balances the wall shear the"
                f" model reports: |dp/dx|*H/2 {closure:.6g} vs rho u_tau^2"
                f" {RHO * uTau ** 2:.6g} Pa ({closure / (RHO * uTau ** 2) - 1:+.2%})")
        require(1.01 <= dp / LB / slope <= 1.06,
                f"seed {seed}: the whole drop is {dp / LB / slope:.4f} of the developed"
                " line - the entry region, not a mismatch")
        require(value(fin, "probeSpread") <= 0.02 * dp,
                f"seed {seed}: section means hug the developed line to"
                f" {value(fin, 'probeSpread') / dp * 100:.2f}% of the drop")
        require(value(fin, "probeSections") >= 8, "enough sections inside the probe window")
        require(nut / NU > 20.0,
                f"seed {seed}: eddy viscosity {nut / NU:.1f}x the molecular one")
        require(value(fin, "uxMax") / plug <= 1.2,
                f"seed {seed}: the profile is flat, u_max/u_bulk"
                f" {value(fin, 'uxMax') / plug:.3f} against 1.5 for a parabola")
        u_pope = popeBulkSpeed(uTau, NU, DH)
        rel = (plug - u_pope) / u_pope
        require(abs(rel) < 0.08,
                f"seed {seed}: Pope's smooth-wall law put at the block's own u_tau"
                f" gives the bulk speed it carries ({plug:.4f} vs {u_pope:.4f} m/s,"
                f" {100 * rel:+.1f}%)")

    if fails:
        raise SystemExit(f"{len(fails)} transient turbulent coupled checks FAILED")
    print("transient turbulent 1D<->3D verification passed")


if __name__ == "__main__":
    main()
