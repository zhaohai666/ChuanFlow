#!/usr/bin/env python3
r"""steadyRectDuct (A10 -> C3 / 5.80): the rectangular duct enters the steady
network. The laminar anchor is the exact Fourier-series solution of
grad^2 u = -G/mu on a W x H rectangle with no-slip walls (self-derived, not
transcribed): bracket(aspect) = 1 - (192 aspect/pi^5) sum_{n odd}
tanh(n pi/(2 aspect))/n^5, Q = (dp/L) long short^3 bracket/(12 mu), and
f.Re_Dh = C_rect = 96/((1+aspect)^2 bracket). Churchill swaps the laminar
member for C_rect (96 at aspect->0 parallel plates, 56.908 for a square). The
circular branch is bit-still (the whole pre-existing suite is the witness).

The checker re-sums the series independently (20001 odd terms vs the kernel's
4001) AND cross-checks the square bracket with a pure-Python 2D finite-difference
Poisson solve (SOR, no scipy), Richardson-extrapolated -- a second, series-free
numerical anchor.

usage: check_steady_rect_duct.py <solver> <case>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

RHO, MU, L = 850.0, 0.05, 10.0
W = 0.04
DP_LAM, DP_TURB = 1000.0, 2.0e6


def bracket(aspect, terms=20001):
    """Independent re-sum of the rectangular-duct laminar bracket."""
    s = 0.0
    for n in range(1, 2 * terms, 2):
        s += math.tanh(n * math.pi / (2.0 * aspect)) / n**5
    return 1.0 - (192.0 * aspect / math.pi**5) * s


def geo(h):
    lo, sh = max(W, h), min(W, h)
    a = sh / lo
    br = bracket(a)
    area = W * h
    dh = 2.0 * W * h / (W + h)
    cond = lo * sh**3 * br / (12.0 * MU * L)
    c = 96.0 / ((1.0 + a)**2 * br)
    return a, br, area, dh, cond, c


def fd_unit_square(n):
    """Solve grad^2 u = -1 on the unit square (u=0 on the boundary) by SOR,
    return the area integral of u. The exact value is bracket(1)/12 -- a
    series-free cross-check of the kernel's Fourier anchor."""
    h = 1.0 / n
    h2 = h * h
    u = [[0.0] * (n + 1) for _ in range(n + 1)]
    omega = 2.0 / (1.0 + math.sin(math.pi / n))
    for _ in range(20000):
        mx = 0.0
        for i in range(1, n):
            ui, up, dn = u[i], u[i - 1], u[i + 1]
            for j in range(1, n):
                gs = 0.25 * (up[j] + dn[j] + ui[j - 1] + ui[j + 1] + h2)
                d = omega * (gs - ui[j])
                ui[j] += d
                if abs(d) > mx:
                    mx = abs(d)
        if mx < 1e-12:
            break
    q = 0.0
    for i in range(1, n):
        q += sum(u[i][1:n])
    return q * h2


def value(line, key):
    m = re.search(r"(?:^|\s)" + re.escape(key) + r"\s+([-+0-9.eE]+)", line)
    if m is None:
        raise SystemExit(f"log line has no '{key}': {line}")
    return float(m.group(1))


def main():
    binary, case = sys.argv[1], sys.argv[2]
    r = subprocess.run([binary, case], capture_output=True, text=True,
                       cwd=os.path.dirname(case))
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:])
        raise SystemExit("steady rectangular-duct solver failed")
    fails = []

    def require(cond, msg):
        if not cond:
            fails.append(msg)
        print(("ok:   " if cond else "FAIL: ") + msg)

    pipes = {}
    for ln in r.stdout.splitlines():
        if ln.startswith("[chflow] pipe "):
            t = ln.split()
            pipes[t[2]] = ln
    head = [ln for ln in r.stdout.splitlines()
            if "nodes" in ln and "pipes" in ln]
    require(bool(head) and "4 nodes, 9 pipes" in head[0],
            f"the case is the shipped 4-node / 9-pipe network ({head})")

    # --- the laminar ladder (dp = 1000): f.Re IS C_rect, Q IS conductance.dp
    ladder = (("RectHundredth", 0.0004), ("RectTenth", 0.004),
              ("RectQuarter", 0.01), ("RectHalf", 0.02), ("RectSquare", 0.04))
    consts = []
    for name, h in ladder:
        ln = pipes[name]
        a, br, area, dh, cond, c = geo(h)
        consts.append(c)
        fr = value(ln, "f") * value(ln, "Re")
        require(abs(fr - c) <= 1e-9 * c,
                f"{name}: f.Re = {fr:.6f} equals the self-derived C_rect = "
                f"{c:.6f} (aspect {a:g}, bracket {br:.6f}) to 1e-9 -- the "
                "all-regime blend's laminar member IS the rectangle solution")
        q = value(ln, "Q")
        require(abs(q - cond * DP_LAM) <= 1e-9 * cond * DP_LAM,
                f"{name}: Q = {q:.6e} is the closed conductance.dp "
                f"({cond * DP_LAM:.6e}) to 1e-9 -- the fallback and law agree")

    require(all(consts[i] > consts[i + 1] for i in range(len(consts) - 1)),
            f"the laminar constant falls monotonically as the duct opens: "
            + " > ".join(f"{c:.3f}" for c in consts)
            + " (aspect 0.01 -> 1)")
    require(abs(consts[0] - 96.0) < 1.5 and consts[0] < 96.0,
            f"the thin-duct member approaches the parallel-plate 96: "
            f"C(0.01) = {consts[0]:.4f} < 96 (the aspect -> 0 limit)")
    require(abs(consts[-1] - 56.908) < 0.01,
            f"the square member is the textbook 56.91: C(1) = "
            f"{consts[-1]:.4f} (bracket = {geo(0.04)[1]:.6f})")

    # --- the rotated twin: a 0.02 x 0.04 duct IS the 0.04 x 0.02 duct
    lr, lh = pipes["RectRotated"], pipes["RectHalf"]
    same = all(abs(value(lr, k) - value(lh, k)) <= 1e-9 * max(abs(value(lh, k)),
                                                             1e-30)
               for k in ("Q", "v", "Re", "f"))
    require(same,
            "RectRotated (W=0.02,H=0.04) matches RectHalf (W=0.04,H=0.02) on "
            "Q/v/Re/f to 1e-9 -- the bore is the same rectangle, the long/short "
            "swap is symmetric")

    # --- the circular twin is bit-still: f.Re = 64, Q = Hagen-Poiseuille.dp
    ln = pipes["CircTwin"]
    q_poi = math.pi * W**4 / (128.0 * MU * L) * DP_LAM
    require(abs(value(ln, "f") * value(ln, "Re") - 64.0) <= 1e-9 * 64.0
            and abs(value(ln, "Q") - q_poi) <= 1e-9 * q_poi,
            "CircTwin: f.Re = 64 and Q = Hagen-Poiseuille.dp to 1e-9 -- the "
            "circular branch never moved (C_rect's aspect->0 member is 96, "
            "not 64; the circle keeps its own 64)")

    # --- the series-free anchor: 2D FD Poisson on the square, Richardson
    exact = bracket(1.0) / 12.0
    q48, q96 = fd_unit_square(48), fd_unit_square(96)
    rich = q96 + (q96 - q48) / 3.0
    require(abs(q96 - exact) / exact < 1e-3
            and abs(q96 - exact) < abs(q48 - exact),
            f"the FD Poisson integral converges toward the series: "
            f"n=48 rel {abs(q48 - exact) / exact:.2e}, n=96 rel "
            f"{abs(q96 - exact) / exact:.2e} (exact bracket/12 = {exact:.8f})")
    require(abs(rich - exact) / exact < 1e-5,
            f"Richardson(48,96) = {rich:.10f} confirms the square bracket to "
            f"{abs(rich - exact) / exact:.2e} -- a second, series-free numerical "
            "method agrees with the kernel's Fourier anchor")

    # --- the wiring, everywhere: Re = rho v Dh/mu, dp = f L rho v^2/2Dh,
    #     v = Q/(W H). Laminar ladder + circular twin + turbulent tier.
    bad = []
    arms = [(n, h, DP_LAM) for n, h in ladder]
    arms += [("CircTwin", 0.0, DP_LAM), ("TurbRectSquare", 0.04, DP_TURB),
             ("TurbCirc", 0.0, DP_TURB)]
    for name, h, dp in arms:
        ln = pipes[name]
        if h > 0:
            area, dh = W * h, 2.0 * W * h / (W + h)
        else:
            area, dh = math.pi * 0.25 * W**2, W
        q, v, re, f = (value(ln, k) for k in ("Q", "v", "Re", "f"))
        if abs(v - q / area) > 1e-9 * max(v, 1e-30) \
                or abs(re - RHO * v * dh / MU) > 1e-9 * max(re, 1e-30) \
                or abs(f * L / dh * RHO * v * v / 2.0 - dp) > 1e-9 * dp:
            bad.append(name)
    require(not bad,
            f"the Dh wiring holds at 1e-9 on every tier (v = Q/A, "
            f"Re = rho v Dh/mu, dp = f L rho v^2/2Dh): mismatches {bad}")

    # --- the turbulent tier is genuinely turbulent and the geometries split
    require(value(pipes["TurbRectSquare"], "Re") > 4000
            and value(pipes["TurbCirc"], "Re") > 4000,
            "the turbulent tier is genuinely turbulent on both geometries "
            f"(square Re = {value(pipes['TurbRectSquare'], 'Re'):.0f}, "
            f"circle Re = {value(pipes['TurbCirc'], 'Re'):.0f})")
    require(abs(value(pipes["TurbRectSquare"], "Q")
                - value(pipes["TurbCirc"], "Q"))
            > 1e-3 * value(pipes["TurbCirc"], "Q"),
            "the rectangular duct and its circular twin carry different rates "
            "at the same turbulent dp -- the bore shape is live, not a circle")

    # --- the rejections: loud, not silent circles (probes in a temp dir)
    text = open(case).read()
    probes = [
        ("rect + annulus",
         text.replace("ductHeight 0.04; }", "ductHeight 0.04; dInner 0.02; }", 1),
         "annulus"),
        ("rect + yieldStress",
         text.replace("ductHeight 0.04; }",
                      "ductHeight 0.04; yieldStress 1.0; }", 1),
         "yield"),
        ("rect + power-law",
         text.replace("ductHeight 0.04; }",
                      "ductHeight 0.04; flowIndex 0.5; consistency 30; }", 1),
         "power-law"),
        ("rect + coupled edge",
         text.replace("ductHeight 0.04; }",
                      "ductHeight 0.04; coupledResistance 1e8; }", 1),
         "coupled"),
        ("negative ductHeight",
         text.replace("ductHeight 0.04; }", "ductHeight -0.01; }", 1),
         "negative ductHeight"),
    ]
    with tempfile.TemporaryDirectory() as td:
        for what, txt, needle in probes:
            p = os.path.join(td, "probe.dict")
            open(p, "w").write(txt)
            rr = subprocess.run([binary, p], capture_output=True, text=True,
                                cwd=td)
            require(rr.returncode != 0 and needle in rr.stdout + rr.stderr,
                    f"rejecting {what}: non-zero exit naming the reason "
                    f"('{needle}')")

    print()
    if fails:
        print(f"FAILED: {len(fails)}")
        for f in fails:
            print("  -", f)
        raise SystemExit(1)
    print("OK: the rectangular section carries the self-derived Fourier-series "
          "closed form as its laminar anchor (f.Re = C_rect to 1e-9, Q to "
          "1e-9), a monotone aspect ladder from the parallel-plate 96 to the "
          "square 56.91, a series-free 2D FD Poisson cross-check (Richardson "
          "to 1e-5), a symmetric rotated twin, Dh-consistent reports at 1e-9 "
          "on every tier, a bit-still circular branch, and loud rejections "
          "outside the declaration scope.")


if __name__ == "__main__":
    main()