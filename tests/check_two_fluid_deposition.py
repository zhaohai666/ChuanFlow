#!/usr/bin/env python3
"""D2 increment 4 (5.99): asphaltene onset + CaCO3 scale deposition,
independently re-derived.

Both species share waxStep's implicit-Euler relaxation
    w_new = (w/dt + k*w_eq) / (1/dt + k),  k = growthRate if w < w_eq else
                                                dissolutionRate
with SPECIES-SPECIFIC declared equilibria:

asphaltene  w_eq = z_a * clamp(drive, 0, 1),  drive = (pAOP - p)/onsetWidth
            (onsetWidth 0 = sharp step). p is the local cell-centre
            HYDRAULIC pressure. For the pressure-driven laminar case the
            profile is exactly linear: loss per cell = dp/nsec with
            dp = 32 mu v L / D^2 and p(cell i) = p_in - (i+0.5)*dp/nsec.
            AOP is a DECLARED experimental input - no kernel AOP model.

scale       w_eq = z_s * clamp(Omega - 1, 0, 1),
            Omega = a(Ca2+) * a(CO3 2-) / Ksp (the standard saturation
            ratio; SI = log10(Omega)). Activities and Ksp are DECLARED; the
            default Ksp cites calcite at 25 C, log Ksp = -8.480 +/- 0.020
            (Plummer & Busenberg 1982, USGS). Only the drive Omega vs 1 is
            the verified physics; the linear yield over Omega in [1, 2] is a
            declared pipeline-geometry model.

The checker reproduces every frame from the previous frame and the reported
vm/dp (asphaltene) or the declared chemistry (scale), never calling a kernel
constitutive function. It also verifies the onset front, the laminar closed
forms, the dissolution arm, and loud rejections.

Usage: check_two_fluid_deposition.py <exe> <case.dict> <asphaltene|scale>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

FAILS = []


def check(name, ok, detail=""):
    if not ok:
        FAILS.append(f"{name} {detail}".strip())


def near(name, got, want, rtol=1e-9, atol=0.0):
    ok = abs(got - want) <= atol + rtol * abs(want)
    check(name, ok, f"got {got:.12g} want {want:.12g}")


def scalar(src, name, default=None):
    m = re.search(r"\b" + re.escape(name) + r"\s+([-\d.eE+]+)\s*;", src)
    return float(m.group(1)) if m else default


def block(src, name):
    m = re.search(r"\b" + name + r"\s*\{(.*?)\n\}", src, re.S)
    return m.group(1) if m else ""


def run(exe, case):
    r = subprocess.run([exe, case], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"case failed ({r.returncode}): "
                           + (r.stdout + r.stderr)[-400:])
    frames, cur = [], None
    for line in r.stdout.splitlines():
        line = line.replace("[chflow] ", "")
        tok = line.split()
        if not tok:
            continue
        if tok[0] == "t" and len(tok) >= 6:
            cur = {"t": float(tok[1]), "vm": float(tok[tok.index("vm") + 1]),
                   "dp": float(tok[tok.index("dp") + 1]), "cells": {}}
            frames.append(cur)
        elif tok[0] == "cell" and cur is not None:
            i = int(tok[1])
            d = {}
            for j in range(2, len(tok) - 1):
                if tok[j] in ("asph", "scale", "alpha"):
                    d[tok[j]] = float(tok[j + 1])
            cur["cells"][i] = d
    return frames


def set_key(text, key, value):
    """Replace the numeric value of `key` with `value` (verbatim text)."""
    pat = re.compile(r"(\b" + re.escape(key) + r"\s+)[-\d.eE+]+(\s*;)")
    out, n = pat.subn(lambda m: m.group(1) + value + m.group(2), text, 1)
    if n != 1:
        raise RuntimeError(f"set_key: key not found: {key}")
    return out


def run_reject(exe, src, what, reps, needle):
    text = src
    for key, value in reps:
        text = set_key(text, key, value)
    fd, p = tempfile.mkstemp(suffix=".dict", prefix="cf_d2i4_reject_")
    os.write(fd, text.encode())
    os.close(fd)
    try:
        rr = subprocess.run([exe, p], capture_output=True, text=True)
        ok = rr.returncode != 0 and needle in rr.stdout + rr.stderr
        check(f"reject: {what}", ok, f"rc={rr.returncode}")
    finally:
        os.remove(p)


def relax(w, dt, k, wEq, cap):
    return min(max((w / dt + k * wEq) / (1.0 / dt + k), 0.0), cap)


def check_asphaltene(exe, case):
    """5.99: onset front in a pressure-driven laminar line."""
    src = open(case).read()
    pb = block(src, "pipe")
    L = float(re.search(r"\blength\s+([-\d.eE+]+)", pb).group(1))
    D = float(re.search(r"\bdiameter\s+([-\d.eE+]+)", pb).group(1))
    nsec = int(float(re.search(r"\bsections\s+([-\d.eE+]+)", pb).group(1)))
    dx = L / nsec
    mu = scalar(src, "muL")
    bb = block(src, "boundaries")
    pIn = float(re.search(r"\binlet\s*\{.*?\bp\s+([-\d.eE+]+)", bb, re.S)
                .group(1))
    pOut = float(re.search(r"\boutlet\s*\{.*?\bp\s+([-\d.eE+]+)", bb, re.S)
                 .group(1))
    ab = block(src, "asphaltene")
    AOP = scalar(ab, "onsetPressure")
    width = scalar(ab, "onsetWidth", 0.0)
    za = scalar(ab, "asphalteneFraction")
    kg = scalar(ab, "growthRate")
    kd = scalar(ab, "dissolutionRate")
    w0 = scalar(ab, "initial", 0.0)
    dt = scalar(src, "dt")
    check("asphaltene case parsed",
          None not in (L, D, mu, pIn, pOut, AOP, za, kg, kd, dt)
          and nsec > 0 and dx > 0,
          f"D={D} L={L} nsec={nsec} pIn={pIn} pOut={pOut} AOP={AOP}")

    frames = run(exe, case)
    check("frames emitted", len(frames) >= 2, f"{len(frames)} frames")
    check("every frame has all cells",
          all(len(f["cells"]) == nsec for f in frames), "cell count")
    check("frame 0 starts at the declared initial deposit",
          all(abs(frames[0]["cells"][i]["asph"] - w0) < 1e-12
              for i in range(nsec)), "frame 0")
    check("no scale token in an asphaltene-only case",
          all("scale" not in f["cells"][i] for f in frames
              for i in range(nsec)), "scale token leaked")

    # steady laminar closed forms: v = dp D^2 / (32 mu L), uniform loss
    last = frames[-1]
    vWant = last["dp"] * D * D / (32.0 * mu * L)
    near("steady vm = dp D^2 / (32 mu L)", last["vm"], vWant)
    near("dp = pIn - pOut", last["dp"], pIn - pOut)
    loss = last["dp"] / nsec
    near("per-cell loss = 32 mu v dx / D^2",
         loss, 32.0 * mu * last["vm"] * dx / (D * D))

    # onset front: cell centres below the AOP hold deposit, above stay clean
    def centreP(n, i):
        return pIn - (i + 0.5) * frames[n]["dp"] / nsec

    below = [i for i in range(nsec) if centreP(len(frames) - 1, i) < AOP]
    above = [i for i in range(nsec) if centreP(len(frames) - 1, i) >= AOP]
    check("onset front splits the line", 0 < len(below) < nsec,
          f"below={below[:3]}..{below[-1:]} n_below={len(below)}")
    check("cells above onset stay clean",
          all(abs(last["cells"][i]["asph"] - w0) < 1e-12 for i in above),
          "deposit above onset")
    check("cells below onset deposit",
          all(last["cells"][i]["asph"] > 0.9 * za for i in below),
          "missing deposit below onset")

    # frame-by-frame implicit-Euler recursion from the reported dp/vm:
    # asphalteneStep used this step's nodePressure, which is exactly the
    # reported profile (uniform loss, fixedPressure inlet march).
    worst = 0.0
    for n in range(1, len(frames)):
        lossN = frames[n]["dp"] / nsec
        for i in range(nsec):
            p = pIn - (i + 0.5) * lossN
            # the kernel's step drive: width 0 = (AOP - p > 0) ? 1 : 0
            if width == 0:
                drive = 1.0 if AOP - p > 0 else 0.0
            else:
                drive = (AOP - p) / width
            aEq = za * min(max(drive, 0.0), 1.0)
            wPrev = frames[n - 1]["cells"][i]["asph"]
            k = kg if wPrev < aEq else kd
            pred = relax(wPrev, dt, k, aEq, za)
            worst = max(worst,
                        abs(frames[n]["cells"][i]["asph"] - pred))
    check("asphaltene recursion reproduced (1e-9)", worst < 1e-9,
          f"worst {worst:.3e}")

    # dissolution arm: drop the AOP below every cell centre and start from a
    # full deposit - everything must redissolve to zero.
    pMin = min(centreP(len(frames) - 1, i) for i in range(nsec))
    dis_src = set_key(set_key(src, "onsetPressure",
                              repr(max(pMin - 1000.0, 1.0))),
                      "initial", repr(za))
    fd, dis = tempfile.mkstemp(suffix=".dict", prefix="cf_d2i4_dissolve_")
    os.write(fd, dis_src.encode())
    os.close(fd)
    try:
        fr = run(exe, dis)
        mx = max(fr[-1]["cells"][i]["asph"] for i in range(nsec))
        # closed form: w decays by 1/(1 + kd*dt) per step from the initial
        steps = len(fr) - 1
        want = za / (1.0 + kd * dt) ** steps
        near("above-onset sibling decays on the closed form", mx, want)
        check("above-onset sibling decays", mx < 0.05 * za,
              f"max {mx:.3e}")
    finally:
        os.remove(dis)

    # loud rejections
    run_reject(exe, src, "missing onsetPressure",
               [("onsetPressure", "0")], "positive onsetPressure")
    run_reject(exe, src.replace("onsetPressure", "onsetWidth -1;\n    onsetPressure", 1),
               "negative onsetWidth", [("onsetPressure", repr(AOP))],
               "negative asphaltene onsetWidth")
    run_reject(exe, src, "asphalteneFraction out of range",
               [("asphalteneFraction", "1.5")],
               "asphalteneFraction outside [0,1]")
    run_reject(exe, src, "negative asphaltene rate",
               [("growthRate", "-1")], "negative asphaltene rate")

    if FAILS:
        print("\n".join("FAIL: " + m for m in FAILS))
        raise SystemExit(f"{len(FAILS)} asphaltene check(s) failed")
    print(f"asphaltene: {len(frames)} frames x {nsec} cells, "
          f"{len(below)} below onset / {len(above)} clean, "
          f"dp {last['dp']:.3f} Pa, worst recursion {worst:.3e} "
          f"- all checks passed")


def check_scale(exe, case):
    """5.99: Omega-driven scale in a sealed static line."""
    src = open(case).read()
    sb = block(src, "scale")
    aCa = scalar(sb, "calciumActivity")
    aCo = scalar(sb, "carbonateActivity")
    ksp = scalar(sb, "solubilityProduct")
    zs = scalar(sb, "scaleFraction")
    kg = scalar(sb, "growthRate")
    kd = scalar(sb, "dissolutionRate")
    w0 = scalar(sb, "initial", 0.0)
    dt = scalar(src, "dt")
    pb = block(src, "pipe")
    nsec = int(float(re.search(r"\bsections\s+([-\d.eE+]+)", pb).group(1)))
    check("scale case parsed",
          None not in (aCa, aCo, ksp, zs, kg, kd, dt) and ksp > 0
          and nsec > 0, f"aCa={aCa} aCo={aCo} Ksp={ksp}")

    omega = aCa * aCo / ksp
    sEq = zs * min(max(omega - 1.0, 0.0), 1.0)
    print(f"  Omega = {omega:.9f}, SI = {math.log10(omega):+.6f}, "
          f"s_eq = {sEq:.9f}")

    frames = run(exe, case)
    check("frames emitted", len(frames) >= 2, f"{len(frames)} frames")
    check("every frame has all cells",
          all(len(f["cells"]) == nsec for f in frames), "cell count")
    check("no asph token in a scale-only case",
          all("asph" not in f["cells"][i] for f in frames
              for i in range(nsec)), "asph token leaked")

    # sealed static line: every cell identical, uniform chemistry
    for f in frames:
        vals = [f["cells"][i]["scale"] for i in range(nsec)]
        check(f"frame {f['t']} cells uniform", max(vals) - min(vals) < 1e-12,
              f"spread {max(vals) - min(vals):.3e}")

    # frame-by-frame implicit-Euler recursion from the declared chemistry
    worst = 0.0
    predLast = frames[0]["cells"][0]["scale"]
    for n in range(1, len(frames)):
        wPrev = frames[n - 1]["cells"][0]["scale"]
        k = kg if wPrev < sEq else kd
        predLast = relax(wPrev, dt, k, sEq, zs)
        worst = max(worst,
                    abs(frames[n]["cells"][0]["scale"] - predLast))
    check("scale recursion reproduced (1e-9)", worst < 1e-9,
          f"worst {worst:.3e}")
    near("final scale matches the closed-form sequence",
         frames[-1]["cells"][0]["scale"], predLast)
    check("scale bounded by z_s",
          all(-1e-12 <= f["cells"][i]["scale"] <= zs + 1e-12
              for f in frames for i in range(nsec)), "bound")

    # undersaturated sibling: scale the activities so Omega < 1 - the line
    # must stay at the initial deposit (no drive either way at w0=0).
    under_src = set_key(src, "calciumActivity",
                        repr(aCa * 0.3 / omega))
    fd, under = tempfile.mkstemp(suffix=".dict", prefix="cf_d2i4_under_")
    os.write(fd, under_src.encode())
    os.close(fd)
    try:
        fr = run(exe, under)
        mx = max(fr[-1]["cells"][i]["scale"] for i in range(nsec))
        check("undersaturated sibling stays clean", mx < 1e-12,
              f"max {mx:.3e}")
    finally:
        os.remove(under)

    # dissolution arm: start from a full deposit in the same undersaturated
    # brine - every cell must decay toward zero via dissolutionRate.
    dis2_src = set_key(under_src, "initial", repr(zs))
    fd, dis = tempfile.mkstemp(suffix=".dict", prefix="cf_d2i4_dis2_")
    os.write(fd, dis2_src.encode())
    os.close(fd)
    try:
        fr = run(exe, dis)
        first, lastf = fr[1]["cells"][0]["scale"], fr[-1]["cells"][0]["scale"]
        check("oversaturated-deposit sibling dissolves",
              lastf < first and lastf < 0.05 * zs,
              f"first {first:.4f} last {lastf:.4f}")
        worstD = 0.0
        sEqU = zs * min(max(0.3 - 1.0, 0.0), 1.0)  # scaled brine: Omega=0.3
        for n in range(1, len(fr)):
            wPrev = fr[n - 1]["cells"][0]["scale"]
            pred = relax(wPrev, dt, kd, sEqU, zs)
            worstD = max(worstD, abs(fr[n]["cells"][0]["scale"] - pred))
        check("dissolution recursion reproduced (1e-9)", worstD < 1e-9,
              f"worst {worstD:.3e}")
    finally:
        os.remove(dis)

    # loud rejections
    run_reject(exe, src, "negative ion activity",
               [("calciumActivity", "-1e-4")],
               "negative scale ion activity")
    run_reject(exe, src, "zero solubilityProduct",
               [("solubilityProduct", "0")],
               "positive solubilityProduct")
    run_reject(exe, src, "scaleFraction out of range",
               [("scaleFraction", "1.5")],
               "scaleFraction outside [0,1]")
    run_reject(exe, src, "negative scale rate",
               [("growthRate", "-1")], "negative scale rate")

    if FAILS:
        print("\n".join("FAIL: " + m for m in FAILS))
        raise SystemExit(f"{len(FAILS)} scale check(s) failed")
    print(f"scale: {len(frames)} frames x {nsec} cells, "
          f"worst recursion {worst:.3e} - all checks passed")


def main():
    if len(sys.argv) != 4:
        raise SystemExit("usage: check_two_fluid_deposition.py <exe> "
                         "<case.dict> <asphaltene|scale>")
    exe, case, mode = sys.argv[1], sys.argv[2], sys.argv[3]
    if mode == "asphaltene":
        return check_asphaltene(exe, case)
    if mode == "scale":
        return check_scale(exe, case)
    raise SystemExit(f"unknown deposition mode {mode}")


if __name__ == "__main__":
    main()
