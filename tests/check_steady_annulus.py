#!/usr/bin/env python3
r"""steadyAnnulus (A10 -> A6 / 5.73): the annulus cross-section enters the
steady network. FB 1958 laminar anchor, Churchill with the laminar member
swapped for C(kappa), Dodge-Metzner Re_GP on Dh without the Rabinowitsch
factor, and the FB power-law profile (lambda by wall matching). The circular
branch is bit-still (the whole pre-existing suite is the witness).

usage: check_steady_annulus.py <solver> <case>
"""
import math
import os
import re
import subprocess
import sys

RHO, MU, L = 850.0, 0.05, 10.0
DP_LAM, DP_TURB = 1.0e4, 4.98e6
D_OUT = 0.06
N, K_POW = 0.5, 30.0


def geo(d_in):
    ro, ri = 0.5 * D_OUT, 0.5 * d_in
    area = math.pi * (ro * ro - ri * ri)
    dh = D_OUT - d_in
    g = ro**4 - ri**4 - (ro * ro - ri * ri) ** 2 / math.log(ro / ri)
    return ro, ri, area, dh, g


def c_lam(d_in):
    _, _, _, dh, g = geo(d_in)
    return 16.0 * dh * dh * (0.25 * D_OUT * D_OUT - 0.25 * d_in * d_in) / g


def fb_integrals(d_in, n, strips=2048):
    ro, ri = 0.5 * D_OUT, 0.5 * d_in
    def parts(lam):
        li = ri2 = 0.0
        inner = outer = 0.0
        h1, h2 = (lam - ri) / strips, (ro - lam) / strips
        for i in range(strips + 1):
            s1 = ri + i * h1
            w = 1.0 if i in (0, strips) else (4.0 if i % 2 else 2.0)
            v1 = max(0.0, lam * lam / s1 - s1) ** (1.0 / n)
            li += w * v1
            inner += w * v1 * (lam * lam - s1 * s1)
            s2 = lam + i * h2
            v2 = max(0.0, s2 - lam * lam / s2) ** (1.0 / n)
            ri2 += w * v2
            outer += w * v2 * (s2 * s2 - lam * lam)
        return (li * h1 - ri2 * h2) / 3.0, \
               0.5 * (inner * h1 + outer * h2) / 3.0
    lo, hi = ri * (1 + 1e-12), ro * (1 - 1e-12)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if parts(mid)[0] < 0.0:
            lo = mid
        else:
            hi = mid
    lam = 0.5 * (lo + hi)
    return lam, parts(lam)[1]


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
        raise SystemExit("steady annulus solver failed")
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
    require(bool(head) and "5 nodes, 12 pipes" in head[0],
            f"the case is the shipped 5-node / 12-pipe network ({head})")

    # --- the exact deep-laminar tier (dp = 1 Pa): f.Re IS the closed constant
    for name, d_in in (("AnnLamEx", 0.04), ("KappaSmallEx", 0.006),
                       ("KappaTinyEx", 6e-6)):
        ln = pipes[name]
        fr = value(ln, "f") * value(ln, "Re")
        require(abs(fr - c_lam(d_in)) <= 1e-9 * c_lam(d_in),
                f"{name}: f.Re = {fr:.6f} equals the FB closed constant "
                f"C({d_in / D_OUT:g}) = {c_lam(d_in):.6f} to 1e-9 -- the "
                "all-regime blend's laminar member IS the annulus solution")
        q = value(ln, "Q")
        q_fb = 1.0 * math.pi * geo(d_in)[4] / (8.0 * MU * L)
        require(abs(q - q_fb) <= 1e-9 * q_fb,
                f"{name}: Q = {q:.6e} is the FB closed form "
                f"({q_fb:.6e}) -- the conductance fallback and the law agree")
    ln = pipes["CircLamEx"]
    q_poi = 1.0 * math.pi * D_OUT ** 4 / (128.0 * MU * L)
    require(abs(value(ln, "f") * value(ln, "Re") - 64.0) <= 1e-9 * 64.0
            and abs(value(ln, "Q") - q_poi) <= 1e-9 * q_poi,
            "CircLamEx: f.Re = 64 and Q = Hagen-Poiseuille to 1e-9 -- the "
            "kappa -> 0 member is the circular pipe, bit for bit")
    ladder = [c_lam(d) for d in (0.04, 0.006, 6e-6)]
    require(ladder[0] > ladder[1] > ladder[2] > 64.0,
            f"the laminar constant falls monotonically toward 64: "
            f"{ladder[0]:.2f} > {ladder[1]:.2f} > {ladder[2]:.2f} > 64 "
            "(the logarithmic kappa -> 0 limit)")

    # --- the power-law tier (dp = 1e4): the FB profile, and C_pl
    lam, dim_q = fb_integrals(0.04, N)
    require(abs(lam - 0.0235) < 0.005 and 0.022 < lam < 0.027,
            f"the FB zero-shear radius sits between the walls: "
            f"lambda = {lam * 1000:.3f} mm in (20, 30)")
    q_fb = (DP_LAM / (2.0 * K_POW * L)) ** (1.0 / N) * 2.0 * math.pi * dim_q
    ln = pipes["AnnPow"]
    require(abs(value(ln, "Q") - q_fb) <= 5e-6 * q_fb,
            f"AnnPow: Q = {value(ln, 'Q'):.6e} matches the independent "
            f"FB profile integral ({q_fb:.6e}, 2048-strip re-solve) to 5e-6")
    a2 = geo(0.04)[2]
    c_pl = 4.0 * geo(0.04)[3] ** (N + 1.0) * a2 ** N \
        / ((2.0 * math.pi * dim_q) ** N * 8.0 ** (N - 1.0))
    fr = value(ln, "f") * value(ln, "Re")
    require(abs(fr - c_pl) <= 1e-9 * c_pl,
            f"AnnPow: f.Re_GP = {fr:.4f} equals the closed C_pl = "
            f"{c_pl:.4f} (K cancels) -- the blend's laminar asymptote is "
            "the FB profile, not a 64/Re line")
    ln2 = pipes["KappaSmallPow"]
    lam2, dq2 = fb_integrals(0.006, N)
    a3 = geo(0.006)[2]
    c_pl2 = 4.0 * geo(0.006)[3] ** (N + 1.0) * a3 ** N \
        / ((2.0 * math.pi * dq2) ** N * 8.0 ** (N - 1.0))
    fr2 = value(ln2, "f") * value(ln2, "Re")
    dm_limit = 64.0 * ((3.0 * N + 1.0) / (4.0 * N)) ** N
    require(abs(fr2 - c_pl2) <= 1e-9 * c_pl2
            and c_pl > c_pl2 > dm_limit,
            f"KappaSmallPow: f.Re_GP = {fr2:.4f} == C_pl = {c_pl2:.4f}, and "
            f"the power-law ladder falls toward the kappa -> 0 member "
            f"{dm_limit:.2f} = 64((3n+1)/4n)^n: {c_pl:.2f} > {c_pl2:.2f}")
    ln3 = pipes["CircPow"]
    q_mr = (math.pi * 0.03 ** 3 * N / (3.0 * N + 1.0)
            * (DP_LAM * 0.03 / (2.0 * L * K_POW)) ** (1.0 / N))
    require(abs(value(ln3, "Q") - q_mr) <= 1e-9 * q_mr
            and abs(value(ln3, "f") * value(ln3, "Re") - 64.0)
            <= 1e-6 * 64.0,
            "CircPow: Q is the Metzner-Reed closed form and f.Re = 64 to "
            "1e-6 -- the circular power-law branch is untouched")

    # --- the wiring, everywhere: Re = rho v Dh / mu and dp = f L rho v^2/2Dh
    bad = []
    for name, dp in (("AnnLam", DP_LAM), ("CircLam", DP_LAM),
                     ("KappaSmall", DP_LAM), ("TurbAnn", DP_TURB),
                     ("TurbCirc", DP_TURB)):
        ln = pipes[name]
        d_in = {"AnnLam": 0.04, "CircLam": 0.0, "KappaSmall": 0.006,
                "TurbAnn": 0.04, "TurbCirc": 0.0}[name]
        _, _, a, dh, _ = geo(d_in) if d_in > 0 else \
            (0, 0, math.pi * 0.25 * D_OUT ** 2, D_OUT, 0)
        q, v, re, f = (value(ln, k) for k in ("Q", "v", "Re", "f"))
        if abs(re - RHO * v * dh / MU) > 1e-9 * max(re, 1e-30) \
                or abs(f * L / dh * RHO * v * v / 2.0 - dp) > 1e-9 * dp:
            bad.append(name)
    require(not bad,
            f"the Dh wiring holds at 1e-9 on every tier (Re = rho v Dh/mu, "
            f"dp = f L rho v^2 / 2Dh): mismatches {bad}")
    require(value(pipes["TurbAnn"], "Re") > 4000
            and value(pipes["TurbCirc"], "Re") > 20000
            and value(pipes["TurbAnn"], "f") > value(pipes["TurbCirc"], "f"),
            "the turbulent tier is genuinely turbulent on both geometries "
            "and the annulus pays the larger friction factor at the same "
            "pressure")

    # --- the rejections: loud, not silent circles (probes mutate this case)
    text = open(case).read()
    probes = [
        ("dInner >= diameter",
         text.replace("dInner 0.04; }", "dInner 0.06; }", 1),
         "annulusInner"),
        ("annulus + yieldStress",
         text.replace("dInner 0.04; }",
                      "dInner 0.04; yieldStress 1.0; }", 1),
         "yield"),
        ("annulus + coupled edge",
         text.replace("dInner 0.04; }",
                      "dInner 0.04; coupledResistance 1e8; }", 1),
         "coupled"),
        ("negative dInner", text.replace("dInner 0.04; }",
                                         "dInner -0.01; }", 1),
         "negative dInner"),
    ]
    for what, txt, needle in probes:
        p = os.path.join(os.path.dirname(case), "probe.dict")
        open(p, "w").write(txt)
        rr = subprocess.run([binary, p], capture_output=True, text=True,
                            cwd=os.path.dirname(case))
        os.remove(p)
        require(rr.returncode != 0 and needle in rr.stdout + rr.stderr,
                f"rejecting {what}: non-zero exit naming the reason "
                f"('{needle}')")

    print()
    if fails:
        print(f"FAILED: {len(fails)}")
        for f in fails:
            print("  -", f)
        raise SystemExit(1)
    print("OK: the annulus section carries the FB 1958 closed form as its "
          "laminar anchor (f.Re = C(kappa) to 1e-9, Q to 1e-9), a "
          "monotone kappa ladder falling to the circular 64, the power-law "
          "FB profile and its own C_pl ladder toward 64((3n+1)/4n)^n, "
          "Dh-consistent reports at 1e-9 on every tier, and loud "
          "rejections outside the declaration scope.")


if __name__ == "__main__":
    main()
