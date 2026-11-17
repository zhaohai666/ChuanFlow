#!/usr/bin/env python3
"""D2 increment 3 (5.98): thin-layer deposit feedback on wall shear,
independently re-derived.

The kernel's wall law normally uses the declared bore D and relative
roughness rough/D. With wax{ feedback yes; } it substitutes, per section,

    d_eff = D sqrt(1 - w)          (exact remaining circular diameter when
                                   the deposit occupies area fraction w)
    rr    = depositRoughness / d_eff

into the SAME Churchill law; area(), transport faces and the thermal area
stay at D. In the laminar Newtonian branch f = 64/Re is roughness-free, so
the whole pressure loss has a closed form with no constitutive helper:

    dp = sum_i 32 mu v dx / d_i^2 = 32 mu v L / (D^2 (1 - w))

for a horizontal, alpha = 1, no-slip, single-phase liquid line. A no-feedback
twin must read back 32 mu v L / D^2, so the ratio is exactly 1/(1-w).
This checker parses D/L/mu/v/w from the case, runs both twins, and re-derives
both losses. It also verifies every reported d_eff against the area identity,
checks that the deposit roughness is observable in the turbulent branch, and
checks loud rejections. It never calls a kernel constitutive function.

Usage: check_two_fluid_deposit_feedback.py <exe> <case.dict>
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


def scalar(src, name):
    m = re.search(r"\b" + re.escape(name) + r"\s+([-\d.eE+]+)\s*;", src)
    return float(m.group(1)) if m else None


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
            cur = {"t": float(tok[1]), "dp": float(tok[tok.index("dp") + 1]),
                   "cells": {}}
            frames.append(cur)
        elif tok[0] == "cell" and cur is not None:
            i = int(tok[1])
            d = {}
            for j in range(2, len(tok) - 1):
                if tok[j] in ("wax", "dEff", "vsl"):
                    d[tok[j]] = float(tok[j + 1])
            cur["cells"][i] = d
    return frames


def write_variant(src, replaces, tag):
    text = src
    for old, new in replaces:
        if old not in text:
            raise RuntimeError(f"variant {tag}: marker not found: {old}")
        text = text.replace(old, new, 1)
    fd, p = tempfile.mkstemp(suffix=".dict", prefix="cf_d2i3_" + tag + "_")
    os.write(fd, text.encode())
    os.close(fd)
    return p


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: check_two_fluid_deposit_feedback.py "
                         "<exe> <case.dict>")
    exe, case = sys.argv[1], sys.argv[2]
    src = open(case).read()

    pb = re.search(r"\bpipe\s*\{(.*?)\n\}", src, re.S).group(1)
    L = float(re.search(r"\blength\s+([-\d.eE+]+)", pb).group(1))
    D = float(re.search(r"\bdiameter\s+([-\d.eE+]+)", pb).group(1))
    nsec = int(float(re.search(r"\bsections\s+([-\d.eE+]+)", pb).group(1)))
    dx = L / nsec
    mu = scalar(src, "muL")
    wb = re.search(r"\bwax\s*\{(.*?)\n\}", src, re.S).group(1)
    w0 = scalar(wb, "initial")
    bb = re.search(r"\binlet\s*\{(.*?)\}", src, re.S).group(1)
    vsl = scalar(bb, "vsl")
    check("case parsed", None not in (D, L, mu, w0, vsl) and all(
        x > 0 for x in (D, L, mu, vsl)) and 0 <= w0 < 1,
        f"D={D} L={L} dx={dx} mu={mu} v={vsl} w={w0}")

    dEffWant = D * math.sqrt(1.0 - w0)
    dpPlainWant = 32.0 * mu * vsl * L / (D * D)
    dpFedWant = 32.0 * mu * vsl * L / (dEffWant * dEffWant)

    frames = run(exe, case)
    twin_path = os.path.join(os.path.dirname(case), "case_nofeedback.dict")
    check("no-feedback twin exists", os.path.exists(twin_path), twin_path)
    plain = run(exe, twin_path) if os.path.exists(twin_path) else []

    check("frames emitted", len(frames) >= 2 and len(plain) >= 2,
          f"{len(frames)} / {len(plain)}")
    check("every frame has all cells",
          all(len(f["cells"]) == nsec for f in frames + plain), "cell count")

    # every reported d_eff is the exact bore-area identity from the reported
    # deposit fraction -- not a stored constant and not the declared bore.
    worstDe = 0.0
    for f in frames:
        for i in range(nsec):
            c = f["cells"][i]
            want = D * math.sqrt(1.0 - c["wax"])
            worstDe = max(worstDe, abs(c["dEff"] - want))
    check("reported d_eff = D sqrt(1-w) (1e-9)", worstDe < 1e-9,
          f"worst {worstDe:.3e}")
    check("d_eff is strictly below the declared bore",
          all(f["cells"][i]["dEff"] < D
              for f in frames for i in range(nsec)), f"D={D}")

    # laminar closed-form pressure loss for both twins, from case constants
    # only; the ratio is 1/(1-w) exactly.
    near("no-feedback dp = 32 mu v L / D^2", plain[-1]["dp"], dpPlainWant)
    near("feedback dp = 32 mu v L / (D^2 (1-w))",
         frames[-1]["dp"], dpFedWant)
    near("dp ratio = 1/(1-w)",
         frames[-1]["dp"] / plain[-1]["dp"], 1.0 / (1.0 - w0))
    check("feedback raises pressure loss",
          frames[-1]["dp"] > plain[-1]["dp"] * 1.5,
          f"{frames[-1]['dp']} vs {plain[-1]['dp']}")

    # no-feedback path is untouched: no dEff token, identical vsl/wax fields
    check("no-feedback twin reports no dEff token",
          all("dEff" not in c for f in plain for c in f["cells"].values()),
          "dEff appeared without feedback")
    worstF = max(abs(f["cells"][i]["vsl"] - plain[n]["cells"][i]["vsl"])
                 for n, f in enumerate(frames) for i in range(nsec))
    check("feedback leaves the flow field unchanged", worstF < 1e-12,
          f"worst vsl diff {worstF:.3e}")

    # deposit roughness must be observable in the TURBULENT branch (the
    # laminar f = 64/Re is roughness-free). Raise the velocity, then compare
    # rough=0 against rough=1e-4 at the same deposit.
    turb = write_variant(src, [("vsl 0.05;", "vsl 5;"),
                               ("dt 10;", "dt 0.5;"),
                               ("endTime 10;", "endTime 1;")], "turb")
    rough = write_variant(open(turb).read(),
                          [("depositRoughness 0;", "depositRoughness 1e-4;")],
                          "rough")
    try:
        ft = run(exe, turb)[-1]["dp"]
        fr = run(exe, rough)[-1]["dp"]
        check("deposit roughness observable in turbulent branch",
              fr > ft * 1.001, f"rough {fr:.3f} vs smooth {ft:.3f}")
    finally:
        os.remove(turb)
        os.remove(rough)

    # loud rejections: negative deposit roughness, and a non-boolean switch.
    rejects = [
        ("negative depositRoughness",
         [("depositRoughness 0;", "depositRoughness -1;")],
         "negative wax depositRoughness"),
        ("non-boolean feedback",
         [("feedback yes;", "feedback maybe;")], "not boolean"),
    ]
    for what, reps, needle in rejects:
        p = write_variant(src, reps, "reject")
        try:
            rr = subprocess.run([exe, p], capture_output=True, text=True)
            ok = rr.returncode != 0 and needle in rr.stdout + rr.stderr
            check(f"reject: {what}", ok, f"rc={rr.returncode}")
        finally:
            os.remove(p)

    if FAILS:
        print("\n".join("FAIL: " + m for m in FAILS))
        raise SystemExit(f"{len(FAILS)} deposit-feedback check(s) failed")
    print(f"deposit-feedback: {len(frames)} frames x {nsec} cells, "
          f"d_eff {dEffWant:.9f}, dp {dpPlainWant:.3f} -> {dpFedWant:.3f} Pa, "
          f"worst d_eff {worstDe:.3e} - all checks passed")


if __name__ == "__main__":
    main()
