#!/usr/bin/env python3
"""Run chflow-coupled-1d3d-multi over two wall-resolved RANS ducts whose inlet is
handed a SHAPE, and check what the shape is worth when the partner it is handed
to is not a duct law.

The laminar version of this case (§5.34/§5.35) could ask a single question: a
declared shape is judged by how far the region still has to climb to reach the
developed profile, and against a laminar duct that profile is the continuum
Poiseuille one, beta = 6/5. The fuller the inlet, the cheaper the region. Here
the partner is kOmegaSST with wall functions, and it develops to something else
entirely -- a flat core with thin wall layers, whose factor this run measures at
1.0206 -- so the two axes can finally meet:

  1. what got handed over   - five arms naming a shape per edge, and the block
                              solving with the shape it was given
  2. the section integral   - the profile carries exactly the rate the network
                              asked for, through every cross-section
  3. the declared shape     - betaImposed recomputed as the midpoint sum over the
                              declared panels, and within 2% of its continuum
                              value -- the shape does not care who owns the walls
  4. the momentum flux      - rho*beta*q^2/A reproduced from the block's own q
  5. development along x    - the SIGN of (beta_in - beta_out) is what the region
                              does with the shape: a plug and 1/7 are under-full
                              and climb; a Poiseuille inlet is over-full for a
                              turbulent core and decays
  6. the plateau is the region's - developed beta is one number over all ten
                              arm x edge pairs, and it is 15% away from the 6/5
                              a laminar duct hands back
  7. the shape's price, paid in pressure - a constant-area duct can only have its
                              pressure pulled down by wall shear, so the parabola
                              arms show a station pair where it RISES: the entry
                              region handing part of the drop back. The plug and
                              1/7 arms decline monotonically.
  8. friction is not the price - over the developed window every arm's gradient
                              sits on 2*rho*u_tau^2/H, the model's own wall
                              balance, to within 10%, for every inlet shape. What
                              the shape moves is the SECANT the block reports
                              across the whole region: 1.01x the developed
                              extrapolation for a plug, 0.28x for a parabola.
  9. the operating point    - that secant is what the parallel bypass sees, so a
                              shape moves the split by 54% (parabola) and 3%
                              (1/7) while the series edge never moves at all
 10. invariances           - per-edge ownership and the interface operator, still
                              exactly themselves against a nonlinear partner
 11. cost and its floor    - every arm converged, one pass is one network solve
                              plus both blocks, and the one arm that did NOT
                              resolve inside its own sweep budget is the reason
                              the interface tolerance is 1e-3 rather than smaller

usage: check_coupled_1d3d_turbulent_profile.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

# fluid and geometry of tests/cases/coupled1d3dTurbulentProfile/case.dict
MU, RHO = 1.0e-3, 1000.0
Q_IN = 6.0e-5
P_REF = 1.0e6
NAME = ("X1", "X2")
ARMS = ("lo", "parab", "p17", "mixed", "parab_qn")
SHAPE = {"lo": ("flat", "flat"),
         "parab": ("parabolic", "parabolic"),
         "p17": ("power1over7", "power1over7"),
         "mixed": ("parabolic", "flat"),
         "parab_qn": ("parabolic", "parabolic")}
TOL = 1e-3
MAX_PASSES = 12
# gap, span, streamwise cells, wall-normal panels, region length -- the wall-
# bounded axis is y in both blocks and the z faces are symmetry planes
DUCT = {"X1": (0.005, 0.005, 20, 8, 0.2), "X2": (0.004, 0.008, 24, 8, 0.24)}
CONTINUUM = {"flat": 1.0, "parabolic": 6.0 / 5.0,
             "power1over7": (7.0 / 9.0) / (7.0 / 8.0) ** 2}
BY = (0.3, 0.008)
LEAD = (0.3, 0.02)
# the number the laminar chain measured and this one must not reproduce
BETA_LAMINAR = 6.0 / 5.0
# the block's own SIMPLE budget, from the case
MAX_SWEEPS = 300.0


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


def developedGradient(gap, uTau):
    """dp/dx implied by the turbulence model's own friction velocity. Only the
    two y-walls carry shear -- the z faces are symmetry planes -- so a control
    volume of area H*W balances pressure against 2*tau_w*W, and
    2*rho*u_tau^2*W/(H*W) = 2*rho*u_tau^2/H."""
    return 2.0 * RHO * uTau * uTau / gap


def parse(lines):
    """One dict per arm: the shapes it was handed, its passes, the network it
    ended with, each block's own report, the beta and pressure stations along
    each region, every per-pass drop a block gave back, and the continuity floor
    each SIMPLE solve started from."""
    arms = []
    shapes = {}
    cur = None

    def close_solve():
        if cur["run"]:
            cur["floors"].append(cur["run"][-1])
            cur["floorLens"].append(len(cur["run"]))
        cur["run"] = []

    for ln in lines:
        if ln.startswith("start ") and " profile " in ln:
            t = ln.split()
            shapes.setdefault(t[1], {})[t[3]] = t[7]
            continue
        if ln.startswith("start ") and " interface " in ln:
            if cur is not None:
                close_solve()
            cur = dict(name=ln.split()[1], shape=dict(shapes.get(ln.split()[1], {})),
                       states={}, passes=[], blocks={}, summary=None, solves=0,
                       stations={}, pstations={}, edges={}, nodes={}, pipes={},
                       floors=[], floorLens=[], run=[], dpBlocks={})
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
            edge = t[5]
            cur["passes"].append(ln)
            cur["states"].setdefault(int(t[3]), {})[edge] = ln
            cur["dpBlocks"].setdefault(edge, []).append(value(ln, "dpBlock"))
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
        elif ln.startswith("pressureProfile "):
            t = ln.split()
            nums = [float(x) for x in t[4:]]
            cur["pstations"][t[3].rstrip(":")] = list(zip(nums[0::2], nums[1::2]))
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
        raise SystemExit("turbulent profile-mapped multi solver failed")

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

    def block(a, n):
        return arms[a]["blocks"][n]

    # ----------------------------------------------------------- 1. the shapes
    for a in ARMS:
        for n, want in zip(NAME, SHAPE[a]):
            got = arms[a]["shape"].get(n)
            require(got == want, f"{a} {n}: handed shape {got}, case says {want}")
    require(arms["mixed"]["shape"] == {"X1": "parabolic", "X2": "flat"},
            "mixed: one run, two different shapes on two edges, and the RANS "
            "partner is on the edge that changed")
    for a in ARMS:
        for n in NAME:
            require(word(block(a, n), "inletShape") == SHAPE[a][NAME.index(n)],
                    f"{a} {n}: the block solved with the shape it was given")
            require(word(block(a, n), "turb") == "kOmegaSST",
                    f"{a} {n}: the partner is the nonlinear one, not a duct law")

    # ------------------------------------- 2. the rate still survives the map
    for a in ARMS:
        for n in NAME:
            b = block(a, n)
            q_net = arms[a]["edges"][n][1]
            imp, fin = value(b, "fluxImposed"), value(b, "fluxIn")
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
    require(close(arms["lo"]["pipes"]["A"], Q_IN, 1e-9)
            and close(arms["lo"]["nodes"]["Snk"], P_REF, 1e-9),
            "the injected rate and the reference pressure are where the case put them")

    # ----------------------------------------- 3. the shape is the shape asked
    for a in ARMS:
        for n in NAME:
            b = block(a, n)
            shape, panels = SHAPE[a][NAME.index(n)], DUCT[n][3]
            got, want = value(b, "betaImposed"), betaDeclared(shape, panels)
            require(close(got, want, 1e-9),
                    f"{a} {n}: betaImposed {got:.7f} against the midpoint sum of "
                    f"{shape} over {panels} panels, {want:.7f}")
            cont = CONTINUUM[shape]
            require(close(got, cont, 0.02),
                    f"{a} {n}: {shape} is within {abs(got / cont - 1) * 100:.2f}% of "
                    f"its continuum factor {cont:.6f} -- the declaration is the "
                    "declaration whatever the partner develops to")
            if shape == "flat":
                require(close(got, 1.0, 1e-12), f"{a} {n}: a plug is beta = 1 exactly")

    # --------------------------------------------- 4. the momentum flux ledger
    for a in ARMS:
        for n in NAME:
            b = block(a, n)
            q, area = value(b, "fluxImposed"), value(b, "inletArea")
            want = RHO * value(b, "betaImposed") * q * q / area
            require(close(value(b, "momentumFluxPa"), want, 1e-12),
                    f"{a} {n}: rho*beta*q^2/A = {want:.7f} Pa is what the block reports")

    # ------------------------------- 5. what the region does with the shape
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
            b = block(a, n)
            shape = SHAPE[a][NAME.index(n)]
            beta_in, beta_out = value(b, "betaIn"), value(b, "betaOut")
            run = [x for _, x in st]
            if shape == "parabolic":
                require(value(b, "betaImposed") > beta_in > beta_out,
                        f"{a} {n}: the inlet was handed OVER-full for a turbulent "
                        f"core ({value(b, 'betaImposed'):.6f} -> {beta_in:.6f} -> "
                        f"{beta_out:.6f}) and the region bleeds it off")
                require(all(run[i + 1] <= run[i] + 5e-4 for i in range(len(run) - 1)),
                        f"{a} {n}: the decay is monotone down the region")
                tail = max(abs(x - beta_out) for x in run[-4:])
                require(tail / beta_out < 2e-4,
                        f"{a} {n}: and it does arrive -- the last four stations are "
                        f"{tail:.2e} off the plateau, {(value(b, 'betaImposed') - beta_out) / tail:.0f}x "
                        "inside the gap the inlet was handed with")
            else:
                require(value(b, "betaImposed") < beta_out and beta_in < beta_out,
                        f"{a} {n}: {shape} enters UNDER-full -- {value(b, 'betaImposed'):.6f} "
                        f"handed, {beta_in:.6f} at the first plane, both below the "
                        f"plateau {beta_out:.6f} -- so the region still has "
                        f"{beta_out - beta_in:.4f} of beta to climb")
                require(all(run[i + 1] >= run[i] - 5e-4 for i in range(len(run) - 1)),
                        f"{a} {n}: the {shape} inlet only gains momentum-flux factor "
                        "as it runs")
    # the sign of (beta_in - beta_out) is the whole story, and it flips here
    over = [(a, n) for a in ARMS for n in NAME
            if value(block(a, n), "betaIn") > value(block(a, n), "betaOut")]
    under = [(a, n) for a in ARMS for n in NAME if (a, n) not in over]
    require(len(over) == 5 and len(under) == 5,
            f"five of the ten arm x edge pairs are fed over-full ({over}) and five "
            "under-full -- against a laminar duct the parabolic arms were the "
            "under-full ones")

    # -------------------------------- 6. the plateau belongs to the region
    betas = {(a, n): value(block(a, n), "betaOut") for a in ARMS for n in NAME}
    spread = max(betas.values()) - min(betas.values())
    require(spread < 1e-3,
            f"developed beta is one number over all ten arm x edge pairs: "
            f"{min(betas.values()):.6f} .. {max(betas.values()):.6f}, spread "
            f"{spread:.2e}")
    require(all(1.005 < v < 1.035 for v in betas.values()),
            "and that number is a turbulent core: "
            + " ".join(f"{a}/{n} {v:.5f}" for (a, n), v in sorted(betas.items())))
    ref = betas["lo", "X1"]
    require(abs(ref - BETA_LAMINAR) / BETA_LAMINAR > 0.10,
            f"the shape relaxes toward {ref:.5f}, {abs(ref / BETA_LAMINAR - 1) * 100:.1f}% "
            f"below the {BETA_LAMINAR:.2f} a laminar duct develops to -- so the "
            "laminar reading of a declared shape does not carry over as a number")

    # --------------------------- 7. the shape's price, paid in pressure
    for a in ARMS:
        for n in NAME:
            ps = arms[a]["pstations"][n]
            st = arms[a]["stations"][n]
            require(len(ps) == len(st) and all(x == y for (x, _), (y, _) in zip(ps, st)),
                    f"{a} {n}: the pressure is reported on the same {len(ps)} planes "
                    "as the shape, so the two curves can be read against each other")
            p = [v for _, v in ps]
            rise = max(p[:8]) - p[0]
            kmax = max(range(len(p)), key=lambda i: p[i])
            shape = SHAPE[a][NAME.index(n)]
            if shape == "parabolic":
                require(rise > 150.0 and kmax <= 8,
                        f"{a} {n}: the entry region HANDS BACK {rise:.2f} Pa -- the "
                        f"pressure peaks at x = {ps[kmax][0]:.4f} m, {ps[kmax][1]:.1f} Pa "
                        "above the inlet plane, in a duct whose walls can only pull it down")
            else:
                require(all(p[i + 1] <= p[i] for i in range(len(p) - 1)),
                        f"{a} {n}: the {shape} arm's pressure falls station by station "
                        "all the way to the outlet, as a wall-bounded duct must")

    # -------------------------------- 8. friction is not what the shape moved
    for a in ARMS:
        for n in NAME:
            b = block(a, n)
            gap, span, nxc, panels, length = DUCT[n]
            g_wall = developedGradient(gap, value(b, "wallUTau"))
            ratio = abs(value(b, "probeSlope")) / g_wall
            require(0.90 < ratio < 1.15,
                    f"{a} {n}: over the developed window the gradient is "
                    f"{abs(value(b, 'probeSlope')):.1f} Pa/m against the wall balance "
                    f"2*rho*u_tau^2/H = {g_wall:.1f} Pa/m ({ratio:.4f}) -- the inlet "
                    "shape costs nothing to the friction itself")
            extr = g_wall * length
            dp = value(b, "dpPa")
            shape = SHAPE[a][NAME.index(n)]
            if shape == "parabolic":
                require(dp / extr < 0.60,
                        f"{a} {n}: yet the secant the block reports across the whole "
                        f"{length} m is {dp / extr:.4f} of the developed drop "
                        f"({dp:.2f} against {extr:.2f} Pa)")
                mom = RHO * (value(b, "betaIn") - value(b, "betaOut")) \
                    * (value(b, "fluxImposed") / value(b, "inletArea")) ** 2
                require((extr - dp) > 1.5 * mom,
                        f"{a} {n}: the {extr - dp:.2f} Pa the region gave back is "
                        f"{(extr - dp) / mom:.2f}x the rho*(beta_in-beta_out)*U^2 = "
                        f"{mom:.2f} Pa the momentum balance alone accounts for")
            else:
                require(0.97 < dp / extr < 1.03,
                        f"{a} {n}: a plug or 1/7 inlet reports {dp / extr:.4f} of the "
                        f"developed drop ({dp:.2f} against {extr:.2f} Pa) -- within "
                        "3% of its own fully developed self")
    for n in NAME:
        rs = {a: value(block(a, n), "dpPa") / arms[a]["edges"][n][1] for a in ARMS}
        require(rs["lo"] > rs["p17"] > rs["parab"],
                f"{n}: R plug {rs['lo']:.6e} > R 1/7 {rs['p17']:.6e} > R parab "
                f"{rs['parab']:.6e} -- the ORDER the laminar case found survives, "
                f"the price is not: the parabola is {rs['lo'] / rs['parab']:.1f}x "
                "cheaper, not 20%")
        require(value(block("p17", n), "probeSpread")
                < value(block("lo", n), "probeSpread"),
                f"{n}: 1/7 crosses the probe window straighter than a plug does "
                f"({value(block('p17', n), 'probeSpread'):.4f} against "
                f"{value(block('lo', n), 'probeSpread'):.4f} Pa)")

    # ------------------------------ 9. where the network ends up sitting
    q_flat = arms["lo"]["edges"]["X1"][1]
    for a in ARMS:
        require(close(arms[a]["edges"]["X2"][1], Q_IN, 1e-12),
                f"{a}: the series edge carries the whole injected rate whatever "
                "its inlet looks like")
        pp = arms[a]["pipes"]
        require(close(pp["A"], Q_IN, 1e-9)
                and close(pp["X1"] + pp["By"], pp["A"], 1e-11),
                f"{a}: node J1 balances to {(pp['X1'] + pp['By']) / pp['A'] - 1:+.2e} "
                f"-- {pp['X1']:.9e} + {pp['By']:.9e} against {pp['A']:.9e}")
        dp_by = pipeChord(BY[0], BY[1], pp["By"]) * pp["By"]
        require(close(arms[a]["edges"]["X1"][0], dp_by, 2e-3),
                f"{a}: the 3D drop between J1 and J2 is the 1D bypass's own "
                f"{dp_by:.2f} Pa at the split rate it forces")
        dp_lead = pipeChord(LEAD[0], LEAD[1], Q_IN) * Q_IN
        require(close(arms[a]["nodes"]["Src"] - arms[a]["nodes"]["J1"], dp_lead, 2e-3),
                f"{a}: the lead pipe's Churchill drop is the head the source pays")
    shift_parab = arms["parab"]["edges"]["X1"][1] / q_flat - 1.0
    shift_p17 = arms["p17"]["edges"]["X1"][1] / q_flat - 1.0
    require(0.45 < shift_parab < 0.65,
            f"a parabolic inlet alone moves the parallel split by {shift_parab * 100:+.2f}% "
            f"({q_flat:.6e} -> {arms['parab']['edges']['X1'][1]:.6e} m3/s)")
    require(0.02 < shift_p17 < 0.05,
            f"and a 1/7 inlet by {shift_p17 * 100:+.2f}% -- an order of magnitude "
            f"less ({shift_parab / shift_p17:.1f}x), which is what a shape only "
            f"{value(block('p17', 'X1'), 'betaOut') - value(block('p17', 'X1'), 'betaIn'):.4f} "
            "of beta below the plateau should do")
    j3 = [arms[a]["nodes"]["J3"] for a in ARMS]
    require(max(j3) - min(j3) < 1e-6,
            f"every arm holds J3 at {j3[0]:.4f} Pa: everything downstream of the "
            "series edge sees the same rate, so only J1 and J2 are allowed to move")

    # ------------------------------------------ 10. what cannot move an edge
    require(close(arms["mixed"]["edges"]["X1"][1],
                  arms["parab"]["edges"]["X1"][1], 1e-9)
            and close(value(block("mixed", "X1"), "dpPa"),
                      value(block("parab", "X1"), "dpPa"), 1e-3),
            "the edge fed a parabola does not care that its partner's inlet changed: "
            f"mixed X1 is parab X1 to {abs(arms['mixed']['edges']['X1'][1] / arms['parab']['edges']['X1'][1] - 1):.1e} "
            f"in rate and {abs(value(block('mixed', 'X1'), 'dpPa') / value(block('parab', 'X1'), 'dpPa') - 1):.1e} "
            "in drop -- the last digit is the block's own sweep count, 29 against 32")
    require(close(arms["mixed"]["edges"]["X2"][1], arms["lo"]["edges"]["X2"][1], 1e-9)
            and close(value(block("mixed", "X2"), "betaImposed"), 1.0, 1e-12)
            and close(value(block("mixed", "X2"), "dpPa"),
                      value(block("lo", "X2"), "dpPa"), 1e-5),
            "and the edge fed a plug is the plug arm's edge, drop and all, with its "
            "own shape")
    for n in NAME:
        require(close(value(block("parab_qn", n), "dpPa"),
                      value(block("parab", n), "dpPa"), 1e-3),
                f"plain sweep and quasi-Newton agree on {n}'s parabolic drop to "
                f"{abs(value(block('parab_qn', n), 'dpPa') / value(block('parab', n), 'dpPa') - 1):.1e}")

    # ----------------------------------------- 11. cost, and the floor on it
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
        require(max(arms[a]["floors"]) < 2e-2,
                f"{a}: every solve starts from a bounded predictor flux "
                f"(worst first-sweep continuity {max(arms[a]['floors']):.3e})")
        last = max(arms[a]["states"])
        for slot, n in zip((-2, -1), NAME):
            b = block(a, n)
            fin = arms[a]["states"][-1][n]
            pass_ = arms[a]["states"][last][n]
            require(value(fin, "sweeps") == value(b, "sweeps")
                    and close(value(fin, "residual"), value(pass_, "residual"), 1e-12),
                    f"{a} {n}: the settled state is the block's own last answer")
            require(arms[a]["floorLens"][slot] + 1 == value(b, "sweeps"),
                    f"{a} {n}: the {value(b, 'sweeps'):.0f}-sweep report is the "
                    f"{arms[a]['floorLens'][slot]}-sweep log solve it belongs to")
    for n in NAME:
        for a in ("lo", "p17", "mixed"):
            require(word(block(a, n), "blockConverged") == "1"
                    and value(block(a, n), "sweeps") < MAX_SWEEPS,
                    f"{a} {n}: a laminar-shaped or plug inlet resolves inside the "
                    f"block's own budget at {value(block(a, n), 'sweeps'):.0f} sweeps")
    # the measured reason the interface tolerance is 1e-3
    wander = arms["parab"]["dpBlocks"]["X2"]
    require(word(block("parab", "X2"), "blockConverged") == "0"
            and value(block("parab", "X2"), "sweeps") == MAX_SWEEPS,
            f"the parabolic X2 is the one solve that never resolves on its own: it "
            f"spends all {MAX_SWEEPS:.0f} sweeps on every pass")
    require((max(wander) - min(wander)) / max(wander) < 2e-3,
            f"and at a rate pinned by topology it still wanders {(max(wander) - min(wander)) / max(wander) * 100:.2f}% "
            f"between passes ({min(wander):.2f} .. {max(wander):.2f} Pa) -- asking "
            f"the interface for less than the partner's own resolution would only "
            "burn passes")

    # ------------------------------------------ the wall treatment stays legal
    yplus = {(a, n): value(block(a, n), "wallYPlusMax") for a in ARMS for n in NAME}
    require(all(11.0 < v < 60.0 for v in yplus.values()),
            "first cell centres sit in the log-law band, so the wall functions are "
            "the ones being used -- "
            + " ".join(f"{a}/{n} {v:.1f}" for (a, n), v in sorted(yplus.items())))
    nu = {(a, n): value(block(a, n), "nutMax") / (MU / RHO) for a in ARMS for n in NAME}
    require(all(v > 10.0 for v in nu.values()),
            "and the eddy viscosity dominates the molecular one everywhere: "
            + " ".join(f"{a}/{n} {v:.0f}x" for (a, n), v in sorted(nu.items())))
    re_bulk = {(a, n): RHO * (abs(arms[a]["edges"][n][1]) / DUCT[n][0] / DUCT[n][1])
               * (2.0 * DUCT[n][0] * DUCT[n][1] / (DUCT[n][0] + DUCT[n][1])) / MU
               for a in ARMS for n in NAME}
    require(all(v > 2000.0 for v in re_bulk.values()),
            "every region is fed well past the laminar limit -- "
            + " ".join(f"{a}/{n} Re {v:.0f}" for (a, n), v in sorted(re_bulk.items())))

    print("note: beta in/out per arm " + "  ".join(
        f"{a} {value(block(a, n), 'betaIn'):.5f}->{value(block(a, n), 'betaOut'):.5f}"
        for a in ARMS for n in NAME))
    print(f"note: q_X1 plug {q_flat:.6e} / 1-7 {arms['p17']['edges']['X1'][1]:.6e} / "
          f"parab {arms['parab']['edges']['X1'][1]:.6e}, q_X2 always {Q_IN:.3e}")
    print(f"note: split shift parab {shift_parab * 100:+.2f}% / 1-7 {shift_p17 * 100:+.2f}%")
    print("note: dp / developed extrapolation " + "  ".join(
        f"{a}/{n} {value(block(a, n), 'dpPa') / (developedGradient(DUCT[n][0], value(block(a, n), 'wallUTau')) * DUCT[n][4]):.4f}"
        for a in ARMS for n in NAME))
    print("note: interface passes " + " ".join(
        f"{a} {value(arms[a]['summary'], 'iterations'):.0f}" for a in ARMS))

    if fails:
        raise SystemExit(f"{len(fails)} ASSERTION(S) FAILED")
    print("PASS: the interface hands a SHAPE to a nonlinear partner, the partner "
          "answers with the profile it actually develops to, and the price of the "
          "shape is paid in the secant and the operating point -- not in friction")


if __name__ == "__main__":
    sys.exit(main())
