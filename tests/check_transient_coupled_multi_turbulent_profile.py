#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-transient-multi over two wall-function RANS ducts whose
inlet is handed a SHAPE inside every time step, and check what a shape says once
the partner it is handed to is nonlinear AND the question is being re-asked at
eight levels instead of once.

The topology is the multi-edge one the profile chain has used since §5.34:
Src -> A -> { X1 (3D RANS duct) || By (1D bypass) } -> X2 (3D RANS duct) -> B -> Snk,
with a constant 6.0e-5 m3/s switched on at t = 0 so every rate starts at zero and
ramps. Five arms differ only in the `profile` table: `lo` flat on both edges (the
control), `parab` and `p17` one fuller shape on both, `mixed` parabolic on the
parallel edge and flat on the series one, `parab_qn` the parabolic arm through the
accelerated operator.

Three earlier runs each held two of the three axes. §5.35 put a shape inside a
time level against a laminar duct law, §5.36 put a shape in front of a RANS
partner at one operating point, and §5.30 put a RANS partner inside a time level
without a shape. Only here are all three present, and the two things that could
only have been settled here are:

  1. the declaration, level by level   - five arms, and the shape the case names is
                                        the shape each edge reports at all eight
                                        levels and each block solves with
  2. the trade is rate-preserving      - the section integral of the condition the
                                        inlet received IS that level's rate, and the
                                        factor read back off the boundary faces is
                                        the declared shape's own discrete value
  3. the sign law survives time        - the plateau the region develops to is one
                                        number over 80 states, and the sign of
                                        (beta_in - beta_out) still sorts the arms
  4. the recovery has a time signature - `pRisePa` per level: 40 states hand back
                                        exactly 0, 40 hand back a rise, and the rise
                                        is the region's rho*(beta_in - beta_out)*U^2
                                        to under 1.2% on the edge that is free to
                                        move and constant on the one that is not
  5. the 1D rows close                 - the momentum row's drop is the node pair's
                                        drop, the chord identity is exact, and the
                                        attribution is still the network's own
  6. the level-1 split is inertial     - the first accepted level is the pure
                                        inertance division to 0.6%, and a shape
                                        moves it by 2.6 percentage points
  7. the shape moves the ramp, not 3%  - tau doubles (2.06x) where §5.35's laminar
                                        partner moved it by 3.29%, and the secant
                                        formula keeps the RATIO while over-counting
                                        the absolute by the same 1.5x in all five arms
  8. separability against a RANS ramp  - the parallel edge's whole trajectory is
                                        bit-identical under either partner, the
                                        pinned edge's is identical to 1.6e-5
  9. schedule and cost                 - one schedule for the network, one solve per
                                        block per pass, and the worst edge names
                                        which shape is holding the loop up
 10. the partner stays legal           - log-law y+, eddy viscosity two orders over
                                        molecular, and the wall balance the model
                                        itself implies
 11. the price, and the floor under it - friction is unchanged by the shape; the
                                        price is the secant, and the arm that pays
                                        it most visibly is the one whose pinned
                                        region never resolves inside its own budget

usage: check_transient_coupled_multi_turbulent_profile.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# fluid, geometry and schedule of
# tests/cases/transientCoupledMultiTurbulentProfile/case.dict
MU, RHO = 1.0e-3, 1000.0
NU = MU / RHO
Q_IN = 6.0e-5                  # m3/s injected at Src from t = 0
P_REF = 1.0e6                  # Pa held at Snk
DT, N_STEPS = 0.05, 8          # time { dt 0.05; endTime 0.4; }
MAX_PASSES = 12                # couple{ maxSubIterations }
TOL_INTERFACE = 1e-3           # couple{ tolerance }
MAX_SWEEPS = 300.0             # simple{ maxIterations }, each block's own budget

ARMS = ("lo", "parab", "p17", "mixed", "parab_qn")
NAME = ("X1", "X2")
SHAPE = {"lo": ("flat", "flat"),
         "parab": ("parabolic", "parabolic"),
         "p17": ("power1over7", "power1over7"),
         "mixed": ("parabolic", "flat"),
         "parab_qn": ("parabolic", "parabolic")}
ACCELERATED = {"parab_qn"}
SEED = {"X1": 1.4e7, "X2": 1.36e7}

# gap, span, streamwise cells, region length. The wall-bounded axis is y in both
# blocks and the z faces are symmetry planes, so a section is gap x span.
DUCT = {"X1": (0.005, 0.005, 20, 0.2), "X2": (0.004, 0.008, 24, 0.24)}
PANELS = 8                     # wall-normal panels the inlet patch is made of
BY = (0.3, 0.008)              # the 1D bypass, sharing X1's node pair
LEAD = (0.3, 0.02)             # A and B
NODE_PAIR = {"X1": ("J1", "J2"), "X2": ("J2", "J3")}

CONTINUUM = {"flat": 1.0, "parabolic": 6.0 / 5.0,
             "power1over7": (7.0 / 9.0) / (7.0 / 8.0) ** 2}
BETA_LAMINAR = 6.0 / 5.0       # what a laminar duct hands back; this one must not
# the plateau this wall-function channel develops to, measured over all 80 states
BETA_PLATEAU = (1.0200, 1.0213)
# rho*(beta_in - beta_out)*U^2 measured against the rise the region reports
RISE_ENTRY = (0.985, 1.005)    # the edge whose rate the shape is allowed to move
RISE_SERIES = (1.15, 1.20)     # the edge whose rate topology pins


def word(line, key):
    """The non-numeric value stored under `key` ('inletShape flat ...')."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([A-Za-z0-9_]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return m.group(1)


def value(line, key):
    """First number stored under `key` in a log line ('q 1.5e-4 ...')."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([-+0-9.eE]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return float(m.group(1))


def stations(line):
    """The (x, y) pairs of a `transient betaProfile`/`pressureProfile` line."""
    tail = re.search(r": (.*)$", line).group(1).split()
    return [(float(tail[k]), float(tail[k + 1])) for k in range(0, len(tail), 2)]


def pass_nums(line, key):
    """The run of numbers stored under `key` in a pass line, up to the next word."""
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+((?:[-+0-9.eE]+\s*)+)", line)
    if m is None:
        raise SystemExit(f"pass line has no '{key}': {line}")
    return [float(v) for v in m.group(1).split()]


def weight(shape, eta):
    """The declared profile at a section coordinate eta in [-1, 1] across the
    wall-bounded gap: the same three words the case may name."""
    if shape == "flat":
        return 1.0
    if shape == "parabolic":
        return 1.0 - eta * eta
    if shape == "power1over7":
        return (1.0 - abs(eta)) ** (1.0 / 7.0)
    raise SystemExit(f"unknown shape '{shape}'")


def betaDeclared(shape, panels):
    """<u^2>/<u>^2 of the declared shape as the mesh can see it: the profile
    sampled at the middle of each of `panels` equal panels across the gap, which
    is where a cell-centred field puts its weight."""
    w = [weight(shape, (2.0 * (j + 0.5) / panels - 1.0)) for j in range(panels)]
    mean = sum(w) / len(w)
    return sum(x * x for x in w) / len(w) / (mean * mean)


def developedGradient(gap, uTau):
    """dp/dx implied by the model's own friction velocity. Only the two y-walls
    carry shear, so a control volume of area H*W balances pressure against
    2*tau_w*W and the gradient is 2*rho*u_tau^2/H."""
    return 2.0 * RHO * uTau * uTau / gap


def churchill(reynolds):
    """Churchill's 1977 Darcy friction factor -- published, so the 1D side of the
    pair is judged from outside this repository."""
    reynolds = max(reynolds, 1e-12)
    a = (2.457 * math.log(1.0 / ((7.0 / reynolds) ** 0.9))) ** 16
    b = (37530.0 / reynolds) ** 16
    return 8.0 * ((8.0 / reynolds) ** 12 + 1.0 / (a + b) ** 1.5) ** (1.0 / 12.0)


def pipeChord(length, diameter, q):
    """The 1D pipe's own pressure-per-rate, from Churchill's friction factor."""
    area = 0.25 * math.pi * diameter ** 2
    v = abs(q) / area
    return (churchill(RHO * v * diameter / MU) * length / diameter
            * 0.5 * RHO * v * v / max(abs(q), 1e-300))


def inertance(length, area):
    return RHO * length / area


def main():
    binary, case = sys.argv[1], sys.argv[2]
    workdir = os.path.dirname(case)
    r = subprocess.run([binary, case], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("transient turbulent profile-mapping solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    def check(cond, msg):
        if not cond:
            require(False, msg)

    def close(a, b, rel):
        if b == 0.0:
            return a == 0.0
        return abs(a / b - 1.0) <= rel

    # ---- split the log into the arms ----
    arms = []
    cur = None
    for ln in lines:
        t = ln.split()
        if ln.startswith("transient start ") and " accelerate " in ln:
            cur = {"tag": ln, "arm": word(ln, "start").rstrip(":"), "shape": {},
                   "seed": {}, "partner": {}, "step": {}, "edge": {}, "pass": {},
                   "nodes": {}, "pipes": {}, "block": {}, "beta": {}, "press": {},
                   "final": {}, "summary": "", "solves": 0}
            arms.append(cur)
        elif cur is None:
            continue
        if ln.startswith("SIMPLE iter 1:"):
            cur["solves"] += 1
        elif ln.startswith("transient arm ") and " profile " in ln:
            cur["shape"][t[4]] = t[6]
        elif ln.startswith("transient start ") and " resistance " in ln:
            cur["seed"][word(ln, "edge")] = value(ln, "resistance")
        elif ln.startswith("transient partner "):
            cur["partner"][t[4]] = ln
        elif ln.startswith("transient step ") and " pass " in ln:
            cur["pass"].setdefault(int(t[2]), {})[int(t[4])] = ln
        elif ln.startswith("transient step ") and " edge " in ln:
            cur["edge"].setdefault(int(t[2]), {})[int(t[4])] = ln
        elif ln.startswith("transient step "):
            cur["step"][int(t[2].rstrip(":"))] = ln
        elif ln.startswith("nodes step "):
            cur["nodes"][int(t[2].rstrip(":"))] = ln
        elif ln.startswith("pipes step "):
            cur["pipes"][int(t[2].rstrip(":"))] = ln
        elif ln.startswith("transient block "):
            cur["block"][word(ln, "name").rstrip(":")] = ln
        elif ln.startswith("transient betaProfile "):
            cur["beta"][word(ln, "edge").rstrip(":")] = ln
        elif ln.startswith("transient pressureProfile "):
            cur["press"][word(ln, "edge").rstrip(":")] = ln
        elif ln.startswith("transient coupled start ") and " name " in ln:
            cur["final"][word(ln, "edge")] = ln
        elif ln.startswith("transient coupled start "):
            cur["summary"] = ln

    by_arm = {c["arm"]: c for c in arms}

    def E(a, n, i):
        """The interface state of edge `i` at level `n` of arm `a`."""
        return by_arm[a]["edge"][n][i]

    def B(a, n):
        return by_arm[a]["block"][n]

    def nd(a, n, name):
        t = by_arm[a]["nodes"][n].split()
        return float(dict(zip(t[3::2], t[4::2]))[name])

    def pipeq(a, n, name):
        t = by_arm[a]["pipes"][n].split()
        return float(dict(zip(t[3::2], t[4::2]))[name])

    # ------------------------------------------------ 1. the declaration, level by level
    header = [ln for ln in lines if ln.startswith("transient network+blocks")]
    require(len(header) == 1 and "2 coupled edges, 5 starting sets" in header[0],
            "one network, two coupled edges and five starting sets announced once")
    require([c["arm"] for c in arms] == list(ARMS),
            f"the five arms ran, in the order the case declares them (got "
            f"{[c['arm'] for c in arms]})")
    require(all(c["seed"] == SEED for c in arms),
            "every arm seeds both edges with the same pair, so no arm starts at its"
            f" own shaped resistance (seeds {SEED['X1']:.3e} then {SEED['X2']:.3e})")
    require({c["arm"]: "accelerate yes" in c["tag"] for c in arms}
            == {a: (a in ACCELERATED) for a in ARMS},
            "only parab_qn runs through the accelerated operator, so a difference "
            "between it and parab is the update rule and not the shape")
    for c in arms:
        require(c["shape"] == dict(zip(NAME, SHAPE[c["arm"]])),
                f"{c['arm']}: the case's `profile` reached its two edges as declared"
                f" ({c['shape']['X1']} on X1, {c['shape']['X2']} on X2)")
    declared = {(a, n, i): word(E(a, n, i), "inletShape")
                for a in ARMS for n in range(1, N_STEPS + 1) for i in (0, 1)}
    require(declared == {(a, n, i): SHAPE[a][i]
                         for a in ARMS for n in range(1, N_STEPS + 1) for i in (0, 1)},
            f"all {len(declared)} level x edge states report the shape that arm names:"
            " the interface re-hands the shape at every level, it does not settle "
            "once at the start")
    require(all(word(B(a, n), "inletShape") == SHAPE[a][NAME.index(n)]
                for a in ARMS for n in NAME),
            "and each block's own report carries the same word it was handed")

    # --------------------------------- 2. what the trade preserves: the rate, and the factor
    flux_err = max(abs(value(E(a, n, i), "fluxImposed") / value(E(a, n, i), "q") - 1.0)
                   for a in ARMS for n in range(1, N_STEPS + 1) for i in (0, 1))
    require(flux_err < 1e-12,
            f"the section integral of the condition the inlet received IS that "
            f"level's rate: worst {flux_err:.2e} over 80 states. A shape is bought "
            "with a shape, not with a different amount of fluid")
    beta_err = {(a, n): max(abs(value(E(a, n, i), "betaImposed")
                                / betaDeclared(SHAPE[a][i], PANELS) - 1.0)
                             for i in (0, 1)) for a in ARMS
                for n in range(1, N_STEPS + 1)}
    require(max(beta_err.values()) < 1e-12,
            f"and the factor read back off the boundary faces is the declared "
            f"profile's own value on {PANELS} panels, to {max(beta_err.values()):.1e}"
            " -- the mapping is re-normalised per level, not per run")
    for a in ARMS:
        for i, n in enumerate(NAME):
            vals = [value(E(a, k, i), "betaImposed") for k in range(1, N_STEPS + 1)]
            spread = (max(vals) - min(vals)) / max(vals)
            require(spread < 1e-12,
                    f"{a} {n}: betaImposed is one number across all {N_STEPS} levels "
                    f"({vals[0]:.16g}, spread {spread:.1e}) -- it grades the "
                    "mapping, not the flow")
    for s in ("flat", "parabolic", "power1over7"):
        disc = betaDeclared(s, PANELS)
        require(abs(disc / CONTINUUM[s] - 1.0) < 0.02,
                f"the {PANELS}-panel value of `{s}` is {disc:.6f} against the "
                f"continuum {CONTINUUM[s]:.6f}: the shape does not care who owns "
                "the walls, only how the mesh samples it")
    mfp_err = max(abs(value(B(a, n), "momentumFluxPa")
                      / (RHO * value(B(a, n), "betaImposed")
                         * value(B(a, n), "fluxImposed") ** 2
                         / value(B(a, n), "inletArea")) - 1.0)
                  for a in ARMS for n in NAME)
    require(mfp_err < 1e-12,
            f"the ledger closes on the block's own faces: momentumFluxPa is "
            f"rho*beta*q^2/A to {mfp_err:.1e}")

    # ---------------------------------- 3. the sign law, measured 80 times
    outs = [value(E(a, n, i), "betaOut") for a in ARMS
            for n in range(1, N_STEPS + 1) for i in (0, 1)]
    ins = {(a, n, i): value(E(a, n, i), "betaIn") for a in ARMS
           for n in range(1, N_STEPS + 1) for i in (0, 1)}
    require(min(outs) >= BETA_PLATEAU[0] and max(outs) <= BETA_PLATEAU[1],
            f"every level's outlet reads inside {BETA_PLATEAU} -- "
            f"{min(outs):.8f} .. {max(outs):.8f} over 80 states: the plateau is the "
            "region's, so it does not move with the level or with the shape")
    require(abs(sum(outs) / len(outs) / BETA_LAMINAR - 1.0) > 0.10,
            f"and it is {abs(sum(outs) / len(outs) / BETA_LAMINAR - 1) * 100:.1f}% "
            f"below the 6/5 a laminar duct hands back ({sum(outs) / len(outs):.6f}): "
            "a wall-function core is flat, and that is what the partner develops to")
    full = {(a, n, i): ins[(a, n, i)] - value(E(a, n, i), "betaOut")
            for a in ARMS for n in range(1, N_STEPS + 1) for i in (0, 1)}
    over = {k: v for k, v in full.items() if v > 0}
    under = {k: v for k, v in full.items() if v <= 0}
    require(len(over) == 40 and len(under) == 40,
            f"the inlet is fuller than the plateau in {len(over)} states and poorer "
            f"in {len(under)} -- the five arms split 40/40, and no state sits on the "
            "boundary between them")
    require(all((k[0], k[2]) in {("parab", 0), ("parab", 1), ("mixed", 0),
                                 ("parab_qn", 0), ("parab_qn", 1)} for k in over),
            "and exactly the parabolic-fed states are the full ones: `lo` and `p17` "
            "never hand over more than the region develops to, at any level")
    for a in ARMS:
        for i, n in enumerate(NAME):
            if full[(a, 1, i)] > 0:
                require(all(full[(a, k, i)] > 0 for k in range(1, N_STEPS + 1)),
                        f"{a} {n}: an inlet fuller than the plateau stays fuller at "
                        f"every level (+{min(full[(a, k, i)] for k in range(1, 9)):.4f}"
                        f" .. +{max(full[(a, k, i)] for k in range(1, 9)):.4f} of beta)")
            else:
                require(all(full[(a, k, i)] <= 0 for k in range(1, N_STEPS + 1)),
                        f"{a} {n}: a plug or 1/7 inlet stays under-full at every "
                        f"level ({min(full[(a, k, i)] for k in range(1, 9)):.4f} .. "
                        f"{max(full[(a, k, i)] for k in range(1, 9)):.4f}) -- the "
                        "region has to build the profile itself")
    # the developing arms climb along x, on the block's own stations
    for a in ("lo", "p17"):
        for n in NAME:
            st = stations(by_arm[a]["beta"][n])
            bs = [b for _, b in st]
            gap, span, ncell, length = DUCT[n]
            require(len(st) == ncell
                    and all(abs(xs - (k + 0.5) * length / ncell) < 1e-9
                            for k, (xs, _) in enumerate(st)),
                    f"{a} {n}: the region reports one shape station per cross-section "
                    f"({len(st)} over {length} m), each at its own cell centre")
            require(min(bs[k + 1] - bs[k] for k in range(len(bs) - 1)) > -1e-4
                    and bs[-1] - bs[0] > (0.015 if a == "lo" else 0.005),
                    f"{a} {n}: the profile the region builds for itself climbs "
                    f"{bs[0]:.6f} -> {bs[-1]:.6f} along x with no fall bigger than "
                    "1e-4 -- development, not a boundary artefact")
            require(bs[-1] == value(E(a, N_STEPS, NAME.index(n)), "betaOut"),
                    f"{a} {n}: the last station IS the betaOut the interface logged "
                    "at the level that ended the arm, bit for bit")

    # ------------------------------ 4. the recovery, level by level: pRise
    rise = {(a, n, i): value(E(a, n, i), "pRisePa") for a in ARMS
            for n in range(1, N_STEPS + 1) for i in (0, 1)}
    zero = [k for k, v in rise.items() if v == 0.0]
    pos = [k for k, v in rise.items() if v > 0.0]
    require(len(zero) == 40 and len(pos) == 40,
            f"the rise is exactly 0 in {len(zero)} states and strictly positive in "
            f"{len(pos)} -- the same 40/40 cut as the sign of (beta_in - beta_out), "
            "so the two are the same event measured twice")
    require({k for k in pos} == set(over) and {k for k in zero} == set(under),
            "state for state: a region hands pressure back if and only if it was "
            "fed fuller than it develops to")
    require(all(value(E(a, n, i), "pRiseX") == stations(by_arm[a]["press"][NAME[i]])[0][0]
                for (a, n, i) in zero),
            f"where nothing is handed back the reported station is the region's own "
            f"first one ({stations(by_arm['lo']['press']['X1'])[0][0]} m for X1, "
            f"{stations(by_arm['lo']['press']['X2'])[0][0]} m for X2): the maximum of "
            "a monotone decline is its start, so the instrument cannot invent a rise")
    peak_x = {"X1": 0.045, "X2": 0.035}
    require(all(abs(value(E(a, k, i), "pRiseX") - peak_x[n]) < 1e-6
                for a in ARMS for i, n in enumerate(NAME)
                for k in range(1, N_STEPS + 1)
                if rise[(a, k, i)] > 0),
            f"where something is, it is always the same station -- "
            f"X1 at {peak_x['X1']} m and X2 at {peak_x['X2']} m, at every level "
            "the rate was raised through")
    # the momentum balance the rise is supposed to be
    for a in ARMS:
        for i, n in enumerate(NAME):
            if rise[(a, N_STEPS, i)] == 0.0:
                continue
            q = value(E(a, N_STEPS, i), "q")
            area = value(E(a, N_STEPS, i), "inletArea")
            mom = RHO * (ins[(a, N_STEPS, i)] - value(E(a, N_STEPS, i), "betaOut")) \
                * (q / area) ** 2
            ratios = [rise[(a, k, i)]
                      / (RHO * (ins[(a, k, i)] - value(E(a, k, i), "betaOut"))
                         * (value(E(a, k, i), "q") / value(E(a, k, i), "inletArea")) ** 2)
                      for k in range(1, N_STEPS + 1)]
            band = RISE_SERIES if i else RISE_ENTRY
            require(all(band[0] < x < band[1] for x in ratios),
                    f"{a} {n}: the rise is {min(ratios):.4f} .. {max(ratios):.4f} of "
                    f"rho*(beta_in-beta_out)*U^2 ({mom:.4f} Pa against the "
                    f"{rise[(a, N_STEPS, i)]:.4f} measured) -- the entry region "
                    "relaxing the inlet toward its own shape, level by level")
            if i == 0:
                require(max(ratios) - min(ratios) < 0.02,
                        f"{a} {n}: and the ratio does not wander as the ramp runs "
                        f"({min(ratios):.4f} .. {max(ratios):.4f}), so the rise is "
                        "the same physics at every level, not a start-up artefact")
            else:
                spread = max(ratios) - min(ratios)
                require(spread < 1e-2,
                        f"{a} {n}: the pinned edge sits on one ratio to {spread:.1e} "
                        f"({rise[(a, 1, i)]:.4f} Pa at level 1 against "
                        f"{rise[(a, N_STEPS, i)]:.4f} at level {N_STEPS}) -- its rate "
                        "never moved, so nothing about it was allowed to")
    for a in ARMS:
        if rise[(a, N_STEPS, 0)] == 0.0:
            continue
        r1, r8 = rise[(a, 1, 0)], rise[(a, N_STEPS, 0)]
        q1, q8 = value(E(a, 1, 0), "q"), value(E(a, N_STEPS, 0), "q")
        expo = math.log(r8 / r1) / math.log(q8 / q1)
        require(r8 > r1 and 2.05 < expo < 2.30,
                f"{a}: the recovery grows with the ramp -- {r1:.4f} Pa at the first "
                f"accepted level to {r8:.4f} Pa at the last, over a rate that grows "
                f"{q8 / q1:.4f}x, an exponent of {expo:.4f} on q")
    # the same instrument, read off the block's own station table
    for a in ARMS:
        for i, n in enumerate(NAME):
            ps = stations(by_arm[a]["press"][n])
            k8 = max(range(len(ps)), key=lambda k: ps[k][1])
            require(close(rise[(a, N_STEPS, i)], ps[k8][1] - ps[0][1], 0.0)
                    and ps[k8][0] == value(E(a, N_STEPS, i), "pRiseX"),
                    f"{a} {n}: level {N_STEPS}'s rise and peak station are the "
                    "last solve's own pressure table, bit for bit -- the number is "
                    "captured with the solve it belongs to, not left over from the "
                    "arm's ending")

    # ----------------------------------- 5. the 1D side still closes with these edges
    for a in ARMS:
        for n in range(1, N_STEPS + 1):
            require(close(nd(a, n, "J1") - nd(a, n, "J2"), value(E(a, n, 0), "rowDrop"), 1e-9)
                    and close(nd(a, n, "J2") - nd(a, n, "J3"), value(E(a, n, 1), "rowDrop"), 1e-9),
                    f"{a} level {n}: each edge's momentum row is closed against the "
                    "node pair it sits across, and a shape does not enter that row")
            require(close(value(E(a, n, 0), "rowDrop"),
                          value(E(a, n, 0), "dpNetwork")
                          + value(E(a, n, 0), "M") * value(E(a, n, 0), "dq") / DT, 0.0),
                    f"{a} level {n}: rowDrop == dpNetwork + M*dq/dt exactly, in Pa, "
                    "with the inertance the arm installed")
            for i in (0, 1):
                require(abs(value(E(a, n, i), "residual")) < TOL_INTERFACE
                        and close(value(E(a, n, i), "residual"),
                                  abs(value(E(a, n, i), "dpNetwork")
                                      - value(E(a, n, i), "dpBlock"))
                                  / value(E(a, n, i), "dpNetwork"), 1e-3),
                        f"{a} level {n} edge {i}: the residual is the two drops' "
                        "disagreement relatively, and inside the case tolerance")
                require(close(value(E(a, n, i), "dpFriction"),
                              value(E(a, n, i), "dpBlock"), 1e-12)
                        and value(E(a, n, i), "dpInertia") == 0,
                        f"{a} level {n} edge {i}: a steady partner's whole drop is "
                        "the friction half the row negotiates on")
            require(value(E(a, n, 0), "unsteady") == 0
                    and value(E(a, n, 1), "unsteady") == 0
                    and value(E(a, n, 0), "dpInertia") == 0
                    and value(E(a, n, 1), "dpInertia") == 0
                    and value(E(a, n, 0), "rowError") == 0
                    and value(E(a, n, 1), "rowError") == 0,
                    f"{a} level {n}: nothing was attributed to the block -- both "
                    "partners solve steady, so the shape moves the friction half "
                    "and the acceleration stays the network's")
        require(all(close(value(E(a, N_STEPS, i), "M"),
                          inertance(DUCT[n][3], DUCT[n][0] * DUCT[n][1]), 1e-14)
                    for i, n in enumerate(NAME)),
            f"{a}: the two inertances are rho*L/A of the exact bores the blocks "
            f"mesh ({value(E(a, N_STEPS, 0), 'M'):.6g}, "
            f"{value(E(a, N_STEPS, 1), 'M'):.6g}) -- §5.26's attribution is a "
            "third axis and this case holds it at the network")
        require(all(close(value(by_arm[a]["partner"][n], "regionInertance"),
                          value(E(a, N_STEPS, i), "M"), 1e-14)
                    for i, n in enumerate(NAME)),
            f"{a}: and the block's own mesh measures the same rho*V/Area^2 it is "
            "charged with, so the arms would agree even if the bore had been guessed")
    lead = pipeChord(LEAD[0], LEAD[1], Q_IN) * Q_IN
    for a in ARMS:
        require(close(nd(a, N_STEPS, "Src") - nd(a, N_STEPS, "J1"), lead, 1e-9)
                and close(nd(a, N_STEPS, "J3") - nd(a, N_STEPS, "Snk"), lead, 1e-9),
                f"{a}: the two lead pipes are the network's own, and their "
                f"Churchill drop over A and B is {lead:.6f} Pa at the injected rate "
                "-- an external number, and no shape of any kind enters it")
        require(close(nd(a, N_STEPS, "Snk"), P_REF, 1e-12),
                f"{a}: Snk is held at the case's fixed pressure "
                f"({nd(a, N_STEPS, 'Snk'):.4f} Pa)")
        check(all(abs(pipeq(a, k, "A") - Q_IN) < 1e-18
                  and abs(pipeq(a, k, "B") - Q_IN) < 1e-18
                  for k in range(1, N_STEPS + 1)),
              f"{a}: every level, A and B carry the injected rate exactly")

    # ------------------------------------------- 6. the first level is the inertial split
    m_by = inertance(BY[0], 0.25 * math.pi * BY[1] ** 2)
    m_x1 = inertance(DUCT["X1"][3], DUCT["X1"][0] * DUCT["X1"][1])
    share_inertial = m_by / (m_x1 + m_by)
    for a in ARMS:
        q1 = value(E(a, 1, 0), "q")
        q2 = value(E(a, 1, 1), "q")
        require(close(q1 + pipeq(a, 1, "By"), Q_IN, 1e-11)
                and close(q2, Q_IN, 1e-11),
                f"{a} level 1: the pair takes the whole injected rate "
                f"({q1 + pipeq(a, 1, 'By'):.17g}) and the series edge carries it "
                f"({q2:.17g}) -- topology, before any shape has been negotiated")
    lvl1 = {a: value(E(a, 1, 0), "q") / Q_IN for a in ARMS}
    require(all(abs(lvl1[a] / share_inertial - 1.0) < 0.03 for a in ARMS),
            f"the first accepted level is the pure inertance division "
            f"M_By/(M_X1+M_By) = {share_inertial:.6f} to within 3%: "
            + " ".join(f"{a} {lvl1[a]:.6f} ({(lvl1[a] / share_inertial - 1) * 100:+.3f}%)"
                       for a in ("lo", "parab", "p17")))
    require(lvl1["parab"] > lvl1["p17"] > lvl1["lo"]
            and abs(lvl1["mixed"] - lvl1["parab"]) < 1e-12,
            f"and the ordering is the shapes': parab {lvl1['parab']:.6f} > p17 "
            f"{lvl1['p17']:.6f} > lo {lvl1['lo']:.6f} -- a fuller inlet costs less "
            f"drop, so at the same inertial kick the shaped branch is handed more "
            "of the rate")

    # --------------------------------- 7. destination and time constant
    dest = {a: value(E(a, N_STEPS, 0), "q") / Q_IN for a in ARMS}
    require(dest["parab"] > dest["lo"] and dest["p17"] > dest["lo"],
            "the branch share the ramp relaxes TOWARD moves with the shape too: "
            + " ".join(f"{a} {dest[a]:.6f}" for a in ARMS))
    g1 = {a: lvl1[a] / lvl1["lo"] - 1.0 for a in ARMS}
    g8 = {a: dest[a] / dest["lo"] - 1.0 for a in ARMS}
    require(0.025 < g1["parab"] < 0.040 and 0.17 < g8["parab"] < 0.24,
            f"a parabolic inlet is worth {g1['parab'] * 100:+.4f}% at the first level "
            f"and {g8['parab'] * 100:+.4f}% at the last -- the shape's price is paid "
            "over the ramp, not at its start")
    require(0.002 < g1["p17"] < 0.008 and 0.015 < g8["p17"] < 0.030,
            f"and a 1/7 inlet {g1['p17'] * 100:+.4f}% then {g8['p17'] * 100:+.4f}%"
            f", which is what a shape only {abs(full[('p17', N_STEPS, 0)]):.4f} of "
            "beta below the plateau should do")
    tau = {}
    for a in ("lo", "parab", "p17"):
        dq = [abs(value(E(a, n, 0), "dq")) for n in range(1, N_STEPS + 1)]
        ratios = [dq[n + 1] / dq[n] for n in range(1, 6)]
        rr = sum(ratios) / len(ratios)
        tau[a] = -DT / math.log(rr)
        require(max(abs(x / rr - 1.0) for x in ratios) < 0.05,
                f"{a}: the ramp is geometric -- successive |dq| decay by "
                f"{min(ratios):.6f} .. {max(ratios):.6f} (mean {rr:.6f}), so one "
                "time constant describes the whole series")
    require(1.85 < tau["parab"] / tau["lo"] < 2.25,
            f"tau is {tau['lo']:.4f} s on a plug and {tau['parab']:.4f} s on a "
            f"parabola, a ratio of {tau['parab'] / tau['lo']:.4f} -- the shape moves "
            "the SPEED of the ramp by a factor, where the laminar partner of §5.35 "
            "moved it by 3.29%")
    require(1.004 < tau["p17"] / tau["lo"] < 1.02,
            f"and a 1/7 inlet by {tau['p17'] / tau['lo'] - 1:+.2%} "
            f"({tau['p17']:.4f} s) -- {abs(tau['p17'] / tau['lo'] - 1) / abs(tau['parab'] / tau['lo'] - 1) * 100:.1f}% "
            "of the parabola's move on the same axis, which is where a shape that "
            f"only owes {abs(full[('p17', N_STEPS, 0)]):.4f} of beta to the plateau "
            "belongs")
    # the network's own secant formula, arm by arm
    form = {}
    for a in tau:
        r1 = value(E(a, N_STEPS, 0), "R")
        r_by = abs(nd(a, N_STEPS, "J1") - nd(a, N_STEPS, "J2")) / pipeq(a, N_STEPS, "By")
        form[a] = (m_x1 + m_by) / (r1 + r_by)
        require(1.40 < form[a] / tau[a] < 1.62,
                f"{a}: the textbook (M_X1+M_By)/(R_X1+R_By) = {form[a]:.4f} s is "
                f"{form[a] / tau[a]:.4f}x the measured {tau[a]:.4f} s -- the chord of "
                "a region whose inlet is developing is not the slope it settles on")
    require(abs((form["parab"] / form["lo"]) / (tau["parab"] / tau["lo"]) - 1.0) < 0.02,
            f"but the RATIO survives: formula {form['parab'] / form['lo']:.4f} "
            f"against measured {tau['parab'] / tau['lo']:.4f}, because the same "
            f"~1.52x over-count divides out of both arms")

    # ------------------------------------------- 8. separability, inside the ramp
    fields = ("q", "dpBlock", "R", "pRisePa", "betaIn", "betaOut", "rowDrop")
    worst_pm = max(abs(value(E("parab", n, 0), f) - value(E("mixed", n, 0), f))
                   for n in range(1, N_STEPS + 1) for f in fields)
    require(worst_pm == 0.0,
            "the parallel edge's whole trajectory is bit-identical between parab and "
            f"mixed ({worst_pm:.1e} on {len(fields)} fields x {N_STEPS} levels) -- "
            "what the series edge is given cannot reach it, not with a RANS partner "
            "on both sides")
    worst_ml = max(abs(value(E("mixed", n, 1), f) / value(E("lo", n, 1), f) - 1.0)
                   for n in range(1, N_STEPS + 1) for f in ("q", "dpBlock"))
    rise_ml = max(abs(value(E("mixed", n, 1), "pRisePa"))
                  for n in range(1, N_STEPS + 1))
    require(worst_ml < 5e-5 and rise_ml == 0.0,
            f"and the flat half of mixed is lo's series edge to {worst_ml:.2e} "
            f"in rate and drop, with its pressure rise at exactly {rise_ml:.1f} Pa "
            "on all 8 levels: the neighbour's parabola costs it nothing beyond the "
            "block's own sweep history")
    worst_pq = max(abs(value(E("parab_qn", n, 0), f) / value(E("parab", n, 0), f) - 1.0)
                   for n in range(1, N_STEPS + 1) for f in ("q", "pRisePa", "betaOut"))
    require(worst_pq < 1e-4,
            f"the operator changes the path, not the physics: parab and parab_qn "
            f"agree on the shaped edge to {worst_pq:.2e} in rate, rise and shape")
    for a in ARMS:
        for n in range(1, N_STEPS + 1):
            require(close(value(E(a, n, 1), "q"), Q_IN, 1e-11),
                    f"{a} level {n}: the series edge carries the whole injected rate "
                    "however its inlet is shaped -- a shape there moves the drop and "
                    "can never move the rate")
    require(abs(value(E("mixed", N_STEPS, 1), "betaImposed") - 1.0) < 1e-12
            and abs(value(E("mixed", N_STEPS, 0), "betaImposed")
                    - betaDeclared("parabolic", PANELS)) < 1e-12,
            "mixed is a real per-edge table inside the time loop: one arm, one run, "
            "one edge parabolic and its neighbour flat, at every level")

    # --------------------------------------- 9. the schedule and what it costs
    for a in ARMS:
        passes = sum(len(v) for v in by_arm[a]["pass"].values())
        s = by_arm[a]["summary"]
        require(value(s, "steps") == N_STEPS
                and close(value(s, "t"), (N_STEPS - 1) * DT + DT, 1e-12),
                f"{a}: {value(s, 'steps'):.0f} levels to t = {value(s, 't'):.4f} s")
        require(value(s, "residual") < TOL_INTERFACE,
                f"{a}: the interface closed, worst residual "
                f"{value(s, 'residual'):.3e} inside {TOL_INTERFACE}")
        require(by_arm[a]["solves"] == 2 * passes,
                f"{a}: {by_arm[a]['solves']} SIMPLE solves for {passes} passes -- one "
                "network solve and one solve per block, per pass, per level")
        for n in range(1, N_STEPS + 1):
            p = by_arm[a]["pass"][n]
            last = p[max(p)]
            require(1 <= len(p) <= MAX_PASSES,
                    f"{a} level {n}: {len(p)} passes, inside the {MAX_PASSES} the "
                    "case allows")
            res = pass_nums(last, "residuals")
            rate = pass_nums(last, "rates")
            require(len(res) == 2 and len(rate) == 2
                    and all(close(rate[i], value(E(a, n, i), "q"), 1e-12)
                            and abs(res[i] - value(E(a, n, i), "residual")) < 1e-15
                            for i in (0, 1))
                    and max(res) < TOL_INTERFACE,
                    f"{a} level {n}: the level's state IS its last pass -- the two "
                    f"rates and residuals recorded there are what the interface "
                    f"settled on, and the worst of them ({max(res):.3e}) is inside "
                    f"{TOL_INTERFACE}")
            worst_i = int(re.search(r" worst (\d+)", last).group(1))
            require(worst_i in (0, 1) and abs(res[worst_i] - max(res)) < 1e-15,
                    f"{a} level {n}: the stop is taken on edge {worst_i}, and that "
                    "is the edge whose residual is the largest in the same line")
        require(all(value(by_arm[a]["step"][n], "newton") < 1e-11
                    for n in range(1, N_STEPS + 1)),
                f"{a}: the network solve itself is at machine precision at every "
                "level -- the interface, not the linear algebra, is what iterates")
    mixed_worst = {int(re.search(r" worst (\d+)", p).group(1))
                   for v in by_arm["mixed"]["pass"].values() for p in v.values()}
    require(mixed_worst == {0},
            f"and in `mixed` the worst edge is the shaped one in every pass of every "
            f"level ({sorted(mixed_worst)}): the flat series edge is a passenger "
            "even though it is negotiating its own drop")

    # ---------------------------------------- 10. the partner stays a legal RANS solve
    yplus = {(a, n): value(B(a, n), "wallYPlusMax") for a in ARMS for n in NAME}
    require(all(20.0 < v < 35.0 for v in yplus.values()),
            "the first cell centres stay in the log-law band at every level's shape "
            + " ".join(f"{a}/{n} {v:.1f}" for (a, n), v in sorted(yplus.items())))
    nu_ratio = {(a, n): value(B(a, n), "nutMax") / NU for a in ARMS for n in NAME}
    require(all(40.0 < v < 150.0 for v in nu_ratio.values()),
            "the eddy viscosity runs 47-134x the molecular one, so the wall "
            "functions are the closure being used: "
            + " ".join(f"{a}/{n} {v:.0f}x" for (a, n), v in sorted(nu_ratio.items())))
    for a in ARMS:
        for n in NAME:
            gap, span, _, _ = DUCT[n]
            plug = value(B(a, n), "plugSpeed")
            dh = 2.0 * gap * span / (gap + span)
            require(RHO * plug * dh / MU > 4000.0,
                    f"{a} {n}: Re_Dh {RHO * plug * dh / MU:.0f} at the rate the ramp "
                    "ended on")
            require(close(value(B(a, n), "inletArea"), gap * span, 1e-12),
                    f"{a} {n}: the inlet patch is the section the shape was sampled on")
            require(abs(value(B(a, n), "fluxOut") / value(B(a, n), "fluxIn") - 1.0)
                    < 1e-10,
                    f"{a} {n}: what enters the region leaves it to "
                    f"{abs(value(B(a, n), 'fluxOut') / value(B(a, n), 'fluxIn') - 1):.1e}"
                    " -- the shape is redistributed across the section, not created")
            spread = (value(B(a, n), "planeFluxMax")
                      - value(B(a, n), "planeFluxMin")) / value(B(a, n), "fluxIn")
            check(spread < 1e-10,
                  f"{a} {n}: every cross-section of the last solve already carries "
                  f"the rate ({spread:.1e} spread)")
    # the wall balance, per shape
    for a in ARMS:
        for n in NAME:
            gap = DUCT[n][0]
            closure = abs(value(B(a, n), "probeSlope")) * gap / 2.0
            tauw = RHO * value(B(a, n), "wallUTau") ** 2
            require(0.88 < closure / tauw < 1.12,
                    f"{a} {n}: the developed window's gradient balances the model's "
                    f"own wall shear to {(closure / tauw - 1) * 100:+.2f}% "
                    f"({value(B(a, n), 'probeSlope'):.1f} Pa/m against 2*rho*u_tau^2/H)")

    # --------------------------- 11. the price, and the floor under the measurement
    dp8 = {a: value(B(a, "X2"), "dpPa") for a in ARMS}
    require(dp8["lo"] > dp8["p17"] > dp8["parab"],
            f"X2, the pinned edge, at level {N_STEPS}: plug {dp8['lo']:.4f} Pa, "
            f"1/7 {dp8['p17']:.4f} Pa, parabola {dp8['parab']:.4f} Pa -- the shape "
            f"is {(1 - dp8['parab'] / dp8['lo']) * 100:.1f}% of the drop, and it is "
            "paid by an edge whose rate never changed")
    for a in ARMS:
        require(all(close(value(by_arm[a]["final"][n], "dp"), value(B(a, n), "dpPa"), 1e-12)
                    for n in NAME),
                f"{a}: the drop the interface ended with is the drop the block's own "
                "field measures -- the last level's ledger is the run's answer")
    # friction is not the price: what a shaped drop is short by is the entry region
    for a in ARMS:
        for i, n in enumerate(NAME):
            gap, span, _, length = DUCT[n]
            extr = developedGradient(gap, value(B(a, n), "wallUTau")) * length
            dp = value(B(a, n), "dpPa")
            mom = RHO * max(0.0, full[(a, N_STEPS, i)]) \
                * (value(E(a, N_STEPS, i), "q")
                   / value(E(a, N_STEPS, i), "inletArea")) ** 2
            if SHAPE[a][i] == "parabolic":
                require(dp / extr < 0.60,
                        f"{a} {n}: fed a parabola, the region reports {dp:.2f} Pa "
                        f"across {length} m -- {dp / extr:.4f} of the developed "
                        f"extrapolation ({extr:.2f} Pa)")
                require((extr - dp) > 1.5 * mom,
                        f"{a} {n}: the {extr - dp:.2f} Pa it is short by is "
                        f"{(extr - dp) / mom:.2f}x the rho*(beta_in-beta_out)*U^2 = "
                        f"{mom:.2f} Pa the inlet relaxation alone accounts for -- "
                        "the rest is the profile it never had to rebuild")
            else:
                require(0.97 < dp / extr < 1.03,
                        f"{a} {n}: a plug or 1/7 inlet reports {dp / extr:.4f} of "
                        f"the developed drop ({dp:.2f} against {extr:.2f} Pa) -- "
                        "within 3% of its own fully developed self, at the end of "
                        "a ramp")
    r_x1 = {a: value(E(a, N_STEPS, 0), "R") for a in ARMS}
    require(r_x1["lo"] > r_x1["p17"] > r_x1["parab"],
            "the order the laminar chain found survives a RANS partner inside a "
            f"ramp: R_X1 plug {r_x1['lo']:.6e} > 1/7 {r_x1['p17']:.6e} > parabola "
            f"{r_x1['parab']:.6e}, and the parabola is "
            f"{r_x1['lo'] / r_x1['parab']:.1f}x cheaper, not 20%")
    capped = {(a, n): sum(1 for k in range(1, N_STEPS + 1)
                          if value(E(a, k, NAME.index(n)), "sweeps") == MAX_SWEEPS)
              for a in ARMS for n in NAME}
    unconv = {(a, n): sum(1 for k in range(1, N_STEPS + 1)
                          if value(E(a, k, NAME.index(n)), "blockConverged") == 0)
              for a in ARMS for n in NAME}
    require(capped[("parab", "X2")] >= 3 and unconv[("parab", "X2")] == capped[("parab", "X2")]
            and capped[("parab_qn", "X2")] >= 5
            and unconv[("parab_qn", "X2")] == capped[("parab_qn", "X2")],
            f"and the honest wall: the parabolic SERIES block spends all "
            f"{MAX_SWEEPS:.0f} sweeps on {capped[('parab', 'X2')]}/"
            f"{capped[('parab_qn', 'X2')]} of {N_STEPS} levels under the plain loop "
            f"and the accelerated one ({unconv[('parab', 'X2')]}/"
            f"{unconv[('parab_qn', 'X2')]} of them reporting blockConverged 0)")
    require(capped[("lo", "X2")] == capped[("p17", "X2")] == capped[("mixed", "X2")] == 0
            and capped[("lo", "X1")] == capped[("parab", "X1")] == 0,
            "no other edge comes close: a plug or 1/7 inlet on the same pinned rate "
            f"resolves in {value(E('lo', N_STEPS, 1), 'sweeps'):.0f} sweeps, and the "
            f"shaped parallel edge in {value(E('parab', N_STEPS, 0), 'sweeps'):.0f}")
    require(value(E("mixed", N_STEPS, 1), "sweeps") == value(E("lo", N_STEPS, 1), "sweeps")
            and capped[("mixed", "X2")] == 0,
            "mixed's series edge is flat, so it is fast even with a parabolic "
            "neighbour -- the wall belongs to the shape a pinned rate is fed, not "
            "to the series topology")
    require(max(value(E(a, k, i), "residual") for a in ("parab", "parab_qn")
                for k in range(1, N_STEPS + 1) for i in (0, 1)) < TOL_INTERFACE,
            "the interface still closed every one of those levels: the network owns "
            "that rate, so an unresolved partner shows up as a wandering drop, not "
            "as a failure -- which is why the tolerance is 1e-3 and not smaller")

    print("note: beta_out level 1 / 8 per arm+edge " + "  ".join(
        f"{a}/{n} {value(E(a, 1, NAME.index(n)), 'betaOut'):.6f}->"
        f"{value(E(a, N_STEPS, NAME.index(n)), 'betaOut'):.6f}"
        for a in ARMS for n in NAME))
    print("note: pRisePa level 1 -> 8 on the shaped entry edges " + "  ".join(
        f"{a}/X1 {value(E(a, 1, 0), 'pRisePa'):.4f}->{value(E(a, N_STEPS, 0), 'pRisePa'):.4f}"
        for a in ("parab", "mixed", "parab_qn")))
    print(f"note: tau lo {tau['lo']:.4f} s / parab {tau['parab']:.4f} s / "
          f"p17 {tau['p17']:.4f} s, secant formula "
          f"{form['lo']:.4f}/{form['parab']:.4f}/{form['p17']:.4f} s")
    print("note: interface passes per arm " + " ".join(
        f"{a} {sum(len(v) for v in by_arm[a]['pass'].values())}" for a in ARMS))

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: a shape, a wall-function partner and a time level together -- the "
          "sign law and the plateau hold at every level, the entry region's pressure "
          "recovery becomes a curve in time, and the shape moves the ramp by a "
          "factor rather than by percent")


if __name__ == "__main__":
    sys.exit(main())
