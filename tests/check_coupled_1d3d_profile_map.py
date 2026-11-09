#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-multi over two detailed ducts whose inlet is handed a
SHAPE rather than a flat rate, and check that the shape is what the run says it
is -- from the section integral up to the wall shear the region ends with.

Every earlier case in this chain maps one number across the interface: the
network asks for a volume rate, the block spreads it evenly over the inlet
section. That is a legitimate statement about a pipe long enough to forget its
entrance, and a false one about a region that is not: the same rate pushed in as
a plug and in as a Poiseuille profile develops over different lengths, so it
costs different pressure, and the operating point the two solvers agree on is
different too. Nothing in the interface could say which of those a region was
being given.

So the claim this checker tests is that a rate distributed over a section
arrives as a distribution, in ten separate senses:

  1. what got handed over   - six arms, each naming a shape per edge; the words
                              are flat / parabolic / power1over7 and `mixed` is
                              one run with two different shapes on two edges
  2. the section integral   - the profile carries EXACTLY the rate the network
                              asked for: imposed flux == network rate == the
                              flux through the patch, no leakage in the mapping
  3. the declared shape     - the momentum-flux factor the block reads back off
                              its own boundary faces equals the midpoint sum the
                              checker does over the declared profile on the
                              declared panels, and sits near the continuum value
  4. the momentum flux      - rho*beta*q^2/A, the number the 1D row would need,
                              reproduced from the block's own beta and q
  5. development along x    - a plug inlet starts flat and climbs to the
                              duct's own plateau; a parabolic one is there from
                              the first plane; 1/7 is between the two
  6. the plateau is the region's - the developed beta is the same number under
                              all three inlets, to one part in 10^4
  7. resistance and its law - the drop over the measured span lands on the
                              Poiseuille value when the inlet is already the
                              Poiseuille shape, and 20% above it when it is a
                              plug; probe-window scatter drops by two orders
  8. the operating point    - the shape moves the parallel split by 8% and the
                              series edge by nothing, and the bypass's own 1D
                              law closes on the block's drop at that split
  9. invariances            - starting coefficient, interface operator, and the
                              other edge's shape do not enter an edge's answer
 10. schedule and cost      - every settled solve converged, one pass is one
                              network solve plus both blocks, and the shape now
                              reaches the predictor flux, not only the corrector

usage: check_coupled_1d3d_profile_map.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# fluid and geometry of tests/cases/coupled1d3dProfileMap/case.dict
MU, RHO = 0.1, 800.0
Q_IN = 8.0e-4
P_REF = 1.0e6
NAME = ("X1", "X2")
ARMS = ("flat_lo", "flat_high", "parab", "p17", "mixed", "parab_qn")
SHAPE = {"flat_lo": ("flat", "flat"), "flat_high": ("flat", "flat"),
         "parab": ("parabolic", "parabolic"),
         "p17": ("power1over7", "power1over7"),
         "mixed": ("parabolic", "flat"),
         "parab_qn": ("parabolic", "parabolic")}
TOL = 1e-4
MAX_PASSES = 12
# gap, span, streamwise cells, wall-normal panels, region length, probe window
DUCT = {"X1": (0.02, 0.02, 75, 20, 0.5), "X2": (0.01, 0.03, 50, 16, 0.4)}
CONTINUUM = {"flat": 1.0, "parabolic": 6.0 / 5.0,
             "power1over7": (7.0 / 9.0) / (7.0 / 8.0) ** 2}
BY = (1.5, 0.03)
LEAD = (5.0, 0.1)


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
    """Churchill's 1977 Darcy friction factor -- published, so the 1D side is
    judged from outside this repository."""
    re = max(re, 1e-12)
    a = (2.457 * math.log(1.0 / ((7.0 / re) ** 0.9 + 0.27 * roughness))) ** 16
    b = (37530.0 / re) ** 16
    return 8.0 * ((8.0 / re) ** 12 + 1.0 / (a + b) ** 1.5) ** (1.0 / 12.0)


def pipeChord(length, diameter, q):
    area = 0.25 * math.pi * diameter ** 2
    v = abs(q) / area
    return (churchill(RHO * v * diameter / MU) * length / diameter
            * 0.5 * RHO * v * v / max(abs(q), 1e-300))


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


def parse(lines):
    """One dict per arm: the shapes it was handed, its passes, the network it
    ended with, each block's own report, the beta stations along each region,
    and the continuity floor every one of its SIMPLE solves stopped at."""
    arms = []
    shapes = {}
    cur = None

    def close_solve():
        if cur["run"]:
            cur["floors"].append(cur["run"][-1])
            cur["floorLens"].append(len(cur["run"]))
        cur["run"] = []

    for ln in lines:
        # the case names each edge's shape before the arm opens, so these are
        # filed by arm and handed over when the arm is built
        if ln.startswith("start ") and " profile " in ln:
            t = ln.split()
            shapes.setdefault(t[1], {})[t[3]] = t[7]
            continue
        if ln.startswith("start ") and " interface " in ln:
            if cur is not None:
                close_solve()
            cur = dict(name=ln.split()[1], shape=dict(shapes.get(ln.split()[1], {})),
                       states={}, passes=[], blocks={}, summary=None, solves=0,
                       stations={}, edges={}, nodes={}, pipes={}, floors=[],
                       floorLens=[], run=[])
            arms.append(cur)
        elif cur is None:
            continue
        elif ln.startswith("SIMPLE iter 1:"):
            close_solve()
            cur["solves"] += 1
        elif ln.startswith("SIMPLE iter "):
            cur["run"].append(value(ln, "continuity"))
        elif ln.startswith("interface ") and " iter " in ln:
            t = ln.split()
            cur["passes"].append(ln)
            cur["states"].setdefault(int(t[3]), {})[t[5]] = ln
        elif ln.startswith("settled "):
            t = ln.split()
            cur["states"].setdefault(-1, {})[t[3]] = ln
        elif ln.startswith("coupledEdges "):
            t = ln.split()
            cur["edges"] = {t[i]: (float(t[i + 1]), float(t[i + 2]))
                            for i in range(4, len(t), 3)}
        elif ln.startswith("nodes ") or ln.startswith("pipes "):
            t = ln.split()
            cur["nodes" if ln.startswith("nodes ") else "pipes"] = {
                t[i]: float(t[i + 1]) for i in range(2, len(t), 2)}
        elif ln.startswith("block ") and " edge " in ln:
            cur["blocks"][ln.split()[3].rstrip(":")] = ln
        elif ln.startswith("betaProfile "):
            t = ln.split()
            nums = [float(x) for x in t[4:]]
            cur["stations"][t[3].rstrip(":")] = list(zip(nums[0::2], nums[1::2]))
        elif ln.startswith("coupled ") and ": iterations " in ln:
            cur["summary"] = ln
    if cur is not None:
        close_solve()
    return arms


def main():
    binary, case = sys.argv[1], sys.argv[2]
    r = subprocess.run([binary, os.path.basename(case)],
                       capture_output=True, text=True,
                       cwd=os.path.dirname(os.path.abspath(case)))
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        raise SystemExit("profile-mapped multi solver failed")

    lines = [ln.replace("[chflow] ", "") for ln in r.stdout.splitlines()]
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    def close(a, b, rel, abs_=0.0):
        return abs(a - b) <= max(rel * abs(b), abs_)

    arms = {a["name"]: a for a in parse(lines)}
    require(sorted(arms) == sorted(ARMS), f"arms run: {sorted(arms)}")

    # --------------------------------------------------------------- 1. shapes
    for a in ARMS:
        for n, want in zip(NAME, SHAPE[a]):
            got = arms[a]["shape"].get(n)
            require(got == want, f"{a} {n}: handed shape {got}, case says {want}")
    require(arms["mixed"]["shape"] == {"X1": "parabolic", "X2": "flat"},
            "mixed: the shape is a per-edge choice, one arm with two different "
            "inlets in one run")

    blocks = {a: arms[a]["blocks"] for a in ARMS}
    for a in ARMS:
        for n in NAME:
            b = blocks[a][n]
            require(word(b, "inletShape") == SHAPE[a][NAME.index(n)],
                    f"{a} {n}: the block solved with the shape it was given")

    # --------------------------------------------- 2. the rate survives the map
    for a in ARMS:
        for n in NAME:
            b = blocks[a][n]
            q_net = arms[a]["edges"][n][1]
            imp, fin, fout = value(b, "fluxImposed"), value(b, "fluxIn"), value(b, "fluxOut")
            require(close(imp, q_net, 1e-12),
                    f"{a} {n}: the profile carries the asked rate to {imp / q_net - 1:+.2e}"
                    " of it")
            require(close(imp, fin, 1e-12),
                    f"{a} {n}: the same section integral the patch reports")
            lo, hi = value(b, "planeFluxMin"), value(b, "planeFluxMax")
            require(lo <= imp * (1 + 2e-5) and hi >= imp * (1 - 2e-5)
                    and close(hi, lo, 2e-5),
                    f"{a} {n}: every cross-section holds the rate "
                    f"({lo:.6e} .. {hi:.6e} against {imp:.6e})")
    require(close(arms["parab"]["pipes"]["A"], Q_IN, 1e-9)
            and close(arms["parab"]["nodes"]["Snk"], P_REF, 1e-9),
            "the injected rate and the reference pressure are where the case put them")

    # ------------------------------------------- 3. the shape is the shape asked
    for a in ARMS:
        for n in NAME:
            b = blocks[a][n]
            shape, panels = SHAPE[a][NAME.index(n)], DUCT[n][3]
            got, want = value(b, "betaImposed"), betaDeclared(shape, panels)
            require(close(got, want, 1e-9),
                    f"{a} {n}: betaImposed {got:.7f} against the midpoint sum of "
                    f"{shape} over {panels} panels, {want:.7f}")
            cont = CONTINUUM[shape]
            require(close(got, cont, 0.012),
                    f"{a} {n}: {shape} is within {abs(got / cont - 1) * 100:.2f}% of "
                    f"its continuum factor {cont:.6f}")
            if shape == "flat":
                require(close(got, 1.0, 1e-12), f"{a} {n}: a plug is beta = 1 exactly")

    # ------------------------------------------------- 4. momentum flux ledger
    for a in ARMS:
        for n in NAME:
            b = blocks[a][n]
            q, area = value(b, "fluxImposed"), value(b, "inletArea")
            want = RHO * value(b, "betaImposed") * q * q / area
            require(close(value(b, "momentumFluxPa"), want, 1e-12),
                    f"{a} {n}: rho*beta*q^2/A = {want:.7f} Pa is what the block reports")

    # ------------------------------------------ 5. how far it has left to run
    for a in ARMS:
        for n in NAME:
            st = arms[a]["stations"][n]
            gap, span, nxc, panels, length = DUCT[n]
            dx = length / nxc
            xs = [x for x, _ in st]
            require(len(st) >= 8 and all(xs[i + 1] > xs[i] for i in range(len(xs) - 1))
                    and close(xs[0], dx / 2, 1e-9)
                    and close(xs[-1], length - dx / 2, 1e-9),
                    f"{a} {n}: {len(st)} stations march from the first cell plane "
                    f"{xs[0]:.6f} to the last {xs[-1]:.6f} of a {length} m region")
            plateau = st[-1][1]
            b = blocks[a][n]
            shape = SHAPE[a][NAME.index(n)]
            if shape == "flat":
                require(st[0][1] < 1.10,
                        f"{a} {n}: a plug starts at beta {st[0][1]:.4f}, flat as advertised")
                run = [x for _, x in st]
                require(all(run[i + 1] >= run[i] - 5e-4 for i in range(len(run) - 1)),
                        f"{a} {n}: the plug inlet only gains momentum-flux factor as it runs")
                require(st[0][1] < plateau - 0.10,
                        f"{a} {n}: the plug is {plateau - st[0][1]:.4f} of beta away from "
                        "the duct's own profile at the inlet plane")
            elif shape == "parabolic":
                interior = [x for _, x in st[1:]]
                drift = max(abs(x - plateau) for x in interior) / plateau
                require(drift < 2e-3,
                        f"{a} {n}: a Poiseuille inlet is translation invariant -- "
                        f"past the plane that carries the boundary, beta never leaves "
                        f"{plateau:.6f} by more than {drift * 100:.2f}%")
                require(close(st[0][1], value(b, "betaImposed"), 0.006),
                        f"{a} {n}: the boundary plane reads the shape it was handed "
                        f"({st[0][1]:.6f} against {value(b, 'betaImposed'):.6f}) plus "
                        "the entrance cell's own bias")
            else:
                require(1.03 < st[0][1] < 1.12,
                        f"{a} {n}: 1/7 enters between the plug and the solution at "
                        f"{st[0][1]:.4f}")
                run = [x for _, x in st]
                require(all(run[i + 1] >= run[i] - 5e-4 for i in range(len(run) - 1)),
                        f"{a} {n}: 1/7 develops the same way a plug does, from closer in")

    # -------------------------------------------- 6. the plateau is the region's
    for n in NAME:
        vals = [blocks[a][n] and value(blocks[a][n], "betaOut") for a in
                ("flat_lo", "parab", "p17")]
        require((max(vals) - min(vals)) / max(vals) < 2e-4,
                f"{n}: developed beta is one number under all three inlets: "
                + " / ".join(f"{v:.6f}" for v in vals))
        require(1.19 < vals[1] < 1.20,
                f"{n}: the developed factor sits by the continuum Poiseuille 1.2 "
                f"({vals[1]:.6f})")

    # ---------------------------------------- 7. resistance and its own physics
    for n in NAME:
        gap, span, nxc, _, _ = DUCT[n]
        st = arms["parab"]["stations"][n]
        l_meas = st[-1][0] - st[0][0]
        r_analytic = 12.0 * MU * l_meas / (gap ** 3 * span)
        rs = {a: value(blocks[a][n], "dpPa") / arms[a]["edges"][n][1] for a in ARMS}
        require(rs["flat_lo"] > rs["p17"] > rs["parab"],
                f"{n}: R flat {rs['flat_lo']:.6e} > R 1/7 {rs['p17']:.6e} > R parab "
                f"{rs['parab']:.6e} -- the fuller the inlet, the less it has to earn")
        require(close(rs["parab"], r_analytic, 0.04),
                f"{n}: a Poiseuille inlet over L_meas {l_meas:.6f} m gives "
                f"{rs['parab']:.6e} against the duct law {r_analytic:.6e} "
                f"({(rs['parab'] / r_analytic - 1) * 100:+.2f}%)")
        excess = rs["flat_lo"] / r_analytic - 1.0
        require(0.10 < excess < 0.40,
                f"{n}: the same region fed a plug costs {(excess) * 100:+.1f}% more "
                "than its own developed law")
        spread = {a: value(blocks[a][n], "probeSpread") for a in ("flat_lo", "parab", "p17")}
        require(spread["parab"] < 0.1 * spread["flat_lo"],
                f"{n}: the probe window is already a straight line under a shaped inlet "
                f"(scatter {spread['parab']:.4f} Pa against {spread['flat_lo']:.4f} Pa)")
        require(spread["p17"] < spread["flat_lo"],
                f"{n}: and 1/7 sits between them ({spread['p17']:.4f} Pa)")
    # the block's drop is measured between the first and last planes it reports
    for a in ARMS:
        for n in NAME:
            st = arms[a]["stations"][n]
            require(close(value(blocks[a][n], "dpPa") / (st[-1][0] - st[0][0]),
                          abs(value(blocks[a][n], "probeSlope")),
                          0.35 if SHAPE[a][NAME.index(n)] != "parabolic" else 0.03),
                    f"{a} {n}: the drop per measured metre {value(blocks[a][n], 'dpPa') / (st[-1][0] - st[0][0]):.1f} Pa/m "
                    f"against the developed slope {value(blocks[a][n], 'probeSlope'):.1f} Pa/m")

    # -------------------------------------- 8. where the network ends up sitting
    q_flat = arms["flat_lo"]["edges"]["X1"][1]
    q_parab = arms["parab"]["edges"]["X1"][1]
    shift = q_parab / q_flat - 1.0
    require(0.04 < shift < 0.15,
            f"the shape alone moves the parallel split by {shift * 100:+.2f}% "
            f"({q_flat:.6e} -> {q_parab:.6e} m3/s)")
    for a in ARMS:
        require(close(arms[a]["edges"]["X2"][1], Q_IN, 1e-12),
                f"{a}: the series edge carries the whole injected rate whatever "
                "its inlet looks like")
        pp = arms[a]["pipes"]
        require(close(pp["A"], Q_IN, 1e-9)
                and close(pp["X1"] + pp["By"], pp["A"], 1e-12),
                f"{a}: node J1 balances -- {pp['X1']:.9e} + {pp['By']:.9e} "
                f"against {pp['A']:.9e}")
        # the bypass is a 1D pipe: its own law has to close on the block's drop
        dp_by = pipeChord(BY[0], BY[1], pp["By"]) * pp["By"]
        require(close(arms[a]["edges"]["X1"][0], dp_by, 2e-3),
                f"{a}: the 3D drop between J1 and J2 is the 1D bypass's own "
                f"{dp_by:.2f} Pa at the split rate it forces")
        dp_lead = pipeChord(LEAD[0], LEAD[1], Q_IN) * Q_IN
        require(close(arms[a]["nodes"]["Src"] - arms[a]["nodes"]["J1"], dp_lead, 2e-3),
                f"{a}: the lead pipe's Churchill drop is the head the source pays")
    # a channel between walls a gap apart and a span wide: Dh = 2*gap, u = q/(gap*span)
    re_bulk = {(a, n): 2.0 * RHO * abs(arms[a]["edges"][n][1]) / (MU * DUCT[n][1])
               for a in ARMS for n in NAME}
    require(all(200 < v < 900 for v in re_bulk.values()),
            "both regions stay laminar -- "
            + " ".join(f"{a}/{n} Re {v:.0f}" for (a, n), v in sorted(re_bulk.items()))
            + " -- so the difference between the arms is the inlet and nothing else")

    # ------------------------------------------------------- 9. what cannot move it
    for n in NAME:
        lo, hi = (value(blocks["flat_lo"][n], "dpPa") / arms["flat_lo"]["edges"][n][1],
                  value(blocks["flat_high"][n], "dpPa") / arms["flat_high"]["edges"][n][1])
        require(close(lo, hi, 1e-4),
                f"{n}: two starting coefficients three times apart land on one "
                f"resistance ({lo:.6e} / {hi:.6e})")
    require(close(value(blocks["mixed"]["X1"], "dpPa"),
                  value(blocks["parab"]["X1"], "dpPa"), 1e-12)
            and close(arms["mixed"]["edges"]["X1"][1],
                      arms["parab"]["edges"]["X1"][1], 1e-12),
            "the edge fed a parabola does not care that its partner's inlet changed: "
            "mixed X1 is parab X1 to the last bit")
    require(close(arms["mixed"]["edges"]["X2"][1], arms["flat_lo"]["edges"]["X2"][1], 1e-9)
            and close(value(blocks["mixed"]["X2"], "betaImposed"), 1.0, 1e-12),
            "and the edge fed a plug is the flat arm's edge, with its own shape")
    require(close(value(blocks["parab_qn"]["X1"], "dpPa"),
                  value(blocks["parab"]["X1"], "dpPa"), 1e-5),
            "plain sweep and quasi-Newton agree on the parabola's drop to "
            f"{abs(value(blocks['parab_qn']['X1'], 'dpPa') / value(blocks['parab']['X1'], 'dpPa') - 1):.1e}")

    # ------------------------------------------------- 10. who paid for the map
    for a in ARMS:
        s = arms[a]["summary"]
        require("(converged)" in s, f"{a}: the interface did not converge")
        it, worst = value(s, "iterations"), value(s, "worstResidual")
        require(2 <= it <= MAX_PASSES, f"{a}: {it} interface passes")
        require(worst < TOL, f"{a}: worst edge residual {worst:.3e} inside {TOL}")
        require(value(s, "blockSolves") == 2 * it,
                f"{a}: one pass is one network solve and one solve per block")
        require(arms[a]["solves"] == 2 * it,
                f"{a}: {arms[a]['solves']} SIMPLE solves against {2 * it} block solves")
        last = max(arms[a]["states"])
        # the block solves the pass in edge order, so the last floor is X2's and
        # the one before is X1's -- and their sweep counts say so out loud
        floors = arms[a]["floors"]
        require(len(floors) == 2 * it, f"{a}: {len(floors)} floors for {2 * it} solves")
        for slot, n in zip((-2, -1), NAME):
            b = blocks[a][n]
            fin = arms[a]["states"][-1][n]
            pass_ = arms[a]["states"][last][n]
            require(value(fin, "sweeps") == value(b, "sweeps")
                    and close(value(fin, "residual"), value(pass_, "residual"), 1e-12),
                    f"{a} {n}: the settled state is the block's own last answer")
            require(word(b, "blockConverged") == "1" and value(b, "sweeps") < 250,
                    f"{a} {n}: the block stopped on its own gate at "
                    f"{value(b, 'sweeps'):.0f} sweeps")
            require(arms[a]["floorLens"][slot] + 1 == value(b, "sweeps"),
                    f"{a} {n}: the {value(b, 'sweeps'):.0f}-sweep report is the "
                    f"{arms[a]['floorLens'][slot]}-sweep log solve it belongs to")
            # the imposed profile must reach the PREDICTOR flux, not only the
            # correction: an HbyA that still reads one number over the patch
            # leaves the boundary column out of balance by the shape's own size,
            # and here a parabolic inlet measured 6.46e-2 before the faces were
            # adopted and 2.0e-3 / 3.6e-3 after.
            floor = floors[slot]
            require(floor < 7e-2, f"{a} {n}: predictor floor {floor:.3e} stays bounded")
            if SHAPE[a][NAME.index(n)] == "parabolic":
                require(floor < 5e-3,
                        f"{a} {n}: the predictor floor under a shaped inlet is "
                        f"{floor:.3e}, where the same solve used to sit at 6.5e-2 "
                        "while the response still read one number over the patch")
                plug = arms["flat_lo"]["floors"][slot]
                require(floor < plug / 10.0,
                        f"{a} {n}: the same region fed a plug stops at {plug:.3e}, "
                        f"{plug / floor:.1f} times higher, so what the shaped inlet "
                        "gained is the predictor, not the pressure")

    print(f"note: shapes {['/'.join(SHAPE[a]) for a in ARMS]}")
    print(f"note: q_X1 flat {q_flat:.6e} / 1-7 {arms['p17']['edges']['X1'][1]:.6e} / "
          f"parab {q_parab:.6e}, q_X2 always {Q_IN:.3e}")
    print(f"note: inlet shape shifts the parallel split by {shift * 100:+.2f}%")
    print("note: final-pass predictor floors X1 / X2 "
          + "  ".join(f"{a} {arms[a]['floors'][-2]:.2e}/{arms[a]['floors'][-1]:.2e}"
                      for a in ("flat_lo", "parab", "p17")))

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: the interface maps a rate AND the shape it arrives in, and the "
          "region answers with the pressure that shape costs")


if __name__ == "__main__":
    sys.exit(main())
