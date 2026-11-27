#!/usr/bin/env python3
"""5.232: verify the 3D Euler-Euler (two-fluid) kernel.

The checker derives every arm from cases/3dEulerEuler/case.dict (the drag
arm) by rewriting the declared blocks in its tmp directory, runs the app,
and re-derives the physics in Python without trusting the solver:

  * uniform arm  - spatially uniform alphaG and UG = UL: every reported
    quantity stays at its initial value to 1e-14 (the upwind M-matrix and
    the equal-velocity drag pair must be identities, not approximations);
  * drag arm     - uniform alpha, UG != UL, no gravity, no pressure: every
    step is replayed against the closed form of the coupled implicit pair
        UG' = ((1+b) UG + a UL) / (1+a+b),  UL' = (b UG + (1+a) UL) / (1+a+b)
        a = K dt / (alphaG rhoG), b = K dt / (alphaL rhoL)
    i.e. each phase's own response time tau_k = alpha_k rho_k / K (the
    implicit pairing identity), to 1e-9 per step; the mixture momentum
    alphaG rhoG UG + alphaL rhoL UL is invariant to 1e-12;
  * slip arm     - open column, gravity + declared dpdz = -rho_mix g: the
    terminal slip matches the discretization's closed form
        uG - uL = alphaG alphaL (rhoL - rhoG) g / K,   dp/dz = -rho_mix g
    with the mixture momentum conserved;
  * schillerNaumann arm - the same per-step replay with
        K = (3/4) Cd(Re) rhoL alphaG |s| / dp,  Cd = 24/Re (1+0.15 Re^0.687)
        (Cd = 0.44 above Re = 1000), Re = rhoL |s| dp / muL;
  * every frame, every arm: alphaG in [0,1], phase masses conserved to 1e-10;
  * refusals: non-positive rho/mu, alphaG0 outside [0,1], steps < 1,
    dt <= 0, unknown drag model, schillerNaumann without a diameter.

usage: check_3d_euler_euler.py <binary> <case.dict> <tmpdir>
"""
import math
import os
import re
import subprocess
import sys

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


def close(a, b, rel, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


def scalar(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", txt)
    if not m:
        raise RuntimeError("no scalar " + key)
    return float(m.group(1))


def vec(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(([^)]*)\)\s*;", txt)
    if not m:
        raise RuntimeError("no vector " + key)
    return [float(x) for x in m.group(1).split()]


EE_RE = re.compile(
    r"EE (" + NUM + r") alphaGmax (" + NUM + r") UGz (" + NUM + r") ULz (" + NUM +
    r") dpdz (" + NUM + r") alphaGmin (" + NUM + r") massG (" + NUM + r") massL (" + NUM +
    r") cont (" + NUM + r")")


def run(binary, text, tmpdir, tag):
    path = os.path.join(tmpdir, tag + ".dict")
    with open(path, "w") as fh:
        fh.write(text)
    r = subprocess.run([binary, path], capture_output=True, text=True, cwd=tmpdir)
    return r, path


def frames(stdout):
    out = []
    for line in stdout.splitlines():
        m = EE_RE.search(line)
        if m:
            out.append({"t": float(m.group(1)), "amax": float(m.group(2)),
                        "ug": float(m.group(3)), "ul": float(m.group(4)),
                        "dpdz": float(m.group(5)), "amin": float(m.group(6)),
                        "massG": float(m.group(7)), "massL": float(m.group(8))})
    return out


def common_checks(fr, tag):
    check(len(fr) >= 2, tag + ": no EE frames")
    if not fr:
        return
    m0, l0 = fr[0]["massG"], fr[0]["massL"]
    for i, f in enumerate(fr):
        check(-1e-14 <= f["amin"] and f["amax"] <= 1.0 + 1e-14,
              "%s frame %d: alphaG outside [0,1]: [%g, %g]" % (tag, i, f["amin"], f["amax"]))
        check(abs(f["massG"] - m0) <= 1e-10 * max(abs(m0), 1e-300)
              and abs(f["massL"] - l0) <= 1e-10 * max(abs(l0), 1e-300),
              "%s frame %d: phase mass not conserved: %g/%g vs %g/%g"
              % (tag, i, f["massG"], f["massL"], m0, l0))


if len(sys.argv) != 4:
    raise SystemExit("usage: check_3d_euler_euler.py <binary> <case.dict> <tmpdir>")
BINARY, CASE, TMP = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(TMP, exist_ok=True)
BASE = open(CASE).read()

rhoG = scalar(BASE, "rhoG")
rhoL = scalar(BASE, "rhoL")
muL = scalar(BASE, "muL")
alphaG0 = scalar(BASE, "alphaG0")
beta = scalar(BASE, "beta")
dt = scalar(BASE, "dt")
steps = int(scalar(BASE, "steps"))
aG, aL = alphaG0, 1.0 - alphaG0
mG, mL = aG * rhoG, aL * rhoL  # phase mass per unit volume

# ---- arm 1: drag. The base case IS the drag arm --------------------------
r, _ = run(BINARY, BASE, TMP, "drag")
check(r.returncode == 0, "drag arm returned %d: %s" % (r.returncode, (r.stdout + r.stderr)[-500:]))
fr = frames(r.stdout)
common_checks(fr, "drag")
if len(fr) == steps + 1:
    K = beta
    a, b = K * dt / mG, K * dt / mL
    D = 1.0 + a + b
    P0 = mG * fr[0]["ug"] + mL * fr[0]["ul"]
    s0 = abs(fr[0]["ug"] - fr[0]["ul"])
    pairTol = 1e-9 * K * s0
    for i in range(1, len(fr)):
        ug0, ul0 = fr[i - 1]["ug"], fr[i - 1]["ul"]
        ug1 = ((1.0 + b) * ug0 + a * ul0) / D
        ul1 = (b * ug0 + (1.0 + a) * ul0) / D
        check(close(fr[i]["ug"], ug1, 1e-9, 1e-9),
              "drag step %d: UG %g != implicit-pair closed form %g" % (i, fr[i]["ug"], ug1))
        check(close(fr[i]["ul"], ul1, 1e-9, 1e-9),
              "drag step %d: UL %g != implicit-pair closed form %g" % (i, fr[i]["ul"], ul1))
        P = mG * fr[i]["ug"] + mL * fr[i]["ul"]
        check(abs(P - P0) <= 1e-12 * max(abs(P0), 1.0),
              "drag step %d: mixture momentum %g != invariant %g" % (i, P, P0))
        # the pairing identity per phase: alpha_k rho_k (Uk' - Uk)/dt = K (Uother' - Uk')
        check(abs(mG * (fr[i]["ug"] - ug0) / dt - K * (fr[i]["ul"] - fr[i]["ug"]))
              <= pairTol,
              "drag step %d: gas implicit pairing tau_G identity broken" % i)
        check(abs(mL * (fr[i]["ul"] - ul0) / dt - K * (fr[i]["ug"] - fr[i]["ul"]))
              <= pairTol,
              "drag step %d: liquid implicit pairing tau_L identity broken" % i)
    s0, sN = fr[0]["ug"] - fr[0]["ul"], fr[-1]["ug"] - fr[-1]["ul"]
    check(abs(sN) < 1e-6 * max(abs(s0), 1e-300), "drag arm: slip did not relax to zero")
    check(abs(fr[-1]["ug"] - P0 / (mG + mL)) < 1e-6,
          "drag arm: did not relax to the common (momentum) velocity")
else:
    check(False, "drag arm: %d frames != steps+1" % len(fr))

# ---- arm 2: uniform fields stay exactly constant --------------------------
uni = BASE.replace("UG (0 0 0.3);", "UG (0 0 0.1);")
uni = uni.replace("UL (0 0 -0.000240480961923848);", "UL (0 0 0.1);")
uni = uni.replace("steps 200;", "steps 60;")
uni = uni.replace("linearTol 1e-14;", "linearTol 1e-15;")
r, _ = run(BINARY, uni, TMP, "uniform")
check(r.returncode == 0, "uniform arm returned %d: %s" % (r.returncode, (r.stdout + r.stderr)[-500:]))
fr = frames(r.stdout)
common_checks(fr, "uniform")
for i, f in enumerate(fr):
    check(abs(f["amax"] - aG) <= 1e-14 and abs(f["amin"] - aG) <= 1e-14,
          "uniform frame %d: alphaG moved: [%g, %g]" % (i, f["amin"], f["amax"]))
    check(abs(f["ug"] - 0.1) <= 1e-14 and abs(f["ul"] - 0.1) <= 1e-14,
          "uniform frame %d: velocity moved: %g / %g" % (i, f["ug"], f["ul"]))
    check(abs(f["dpdz"]) <= 1e-14, "uniform frame %d: dpdz %g != 0" % (i, f["dpdz"]))

# ---- arm 3: buoyancy slip: open column, gravity + declared dpdz -----------
# The declared gradient dpdz = -rho_mix g supports the mixture weight, so the
# uniform-slip state is an EXACT fixed point of the discretization: each
# phase's buoyancy force is balanced by the drag,
#   uG - uL = alphaG alphaL (rhoL - rhoG) g / K.
g = 9.81
rho_mix = aG * rhoG + aL * rhoL
slip = BASE
slip = slip.replace("gravity (0 0 0);", "gravity (0 0 -%r);" % g)
slip = slip.replace("dpdz 0;", "dpdz %r;" % (-rho_mix * g))
slip = slip.replace("UG (0 0 0.3);", "UG (0 0 0);")
slip = slip.replace("UL (0 0 -0.000240480961923848);", "UL (0 0 0);")
slip = slip.replace("steps 200;", "steps 200;")
slip = slip.replace("reportEvery 1;", "reportEvery 10;")
r, _ = run(BINARY, slip, TMP, "slip")
check(r.returncode == 0, "slip arm returned %d: %s" % (r.returncode, (r.stdout + r.stderr)[-500:]))
fr = frames(r.stdout)
common_checks(fr, "slip")
if fr:
    s_term = aG * aL * (rhoL - rhoG) * g / beta
    dpdz_term = -rho_mix * g
    s_end = fr[-1]["ug"] - fr[-1]["ul"]
    check(close(s_end, s_term, 1e-6, 1e-9),
          "slip arm: terminal slip %g != closed form alphaG alphaL drho g / K = %g" % (s_end, s_term))
    check(close(fr[-1]["dpdz"], dpdz_term, 1e-6, 1e-9),
          "slip arm: dpdz %g != -rho_mix g = %g" % (fr[-1]["dpdz"], dpdz_term))
    # the mixture momentum stays at its initial value (the declared gradient
    # cancels the mixture weight and the drag conserves): with zero initial
    # momentum the heavy liquid barely moves while the gas rises at ~s*
    Pmix = mG * fr[-1]["ug"] + mL * fr[-1]["ul"]
    check(abs(Pmix) <= 1e-6 * (mG + mL) * abs(s_term),
          "slip arm: mixture momentum %g != 0 at terminal slip" % Pmix)

# ---- arm 4: Schiller-Naumann per-step replay ------------------------------
sn = BASE.replace("model constant;\n    beta  1200.0;",
                  "model schillerNaumann;\n    diameter 0.002;")
sn = sn.replace("steps 200;", "steps 200;")
r, _ = run(BINARY, sn, TMP, "schillerNaumann")
check(r.returncode == 0, "schillerNaumann arm returned %d: %s" % (r.returncode, (r.stdout + r.stderr)[-500:]))
fr = frames(r.stdout)
common_checks(fr, "schillerNaumann")
if len(fr) == steps + 1:
    dp = 0.002
    P0 = mG * fr[0]["ug"] + mL * fr[0]["ul"]
    for i in range(1, len(fr)):
        ug0, ul0 = fr[i - 1]["ug"], fr[i - 1]["ul"]
        s = abs(ug0 - ul0)
        if s > 0:
            Re = rhoL * s * dp / muL
            cd = 0.44 if Re > 1000.0 else (24.0 / Re) * (1.0 + 0.15 * Re ** 0.687)
            K = 0.75 * cd * rhoL * aG * s / dp
        else:
            K = 0.0
        a, b = K * dt / mG, K * dt / mL
        D = 1.0 + a + b
        ug1 = ((1.0 + b) * ug0 + a * ul0) / D
        ul1 = (b * ug0 + (1.0 + a) * ul0) / D
        check(close(fr[i]["ug"], ug1, 1e-9, 1e-9),
              "schillerNaumann step %d: UG %g != replay %g" % (i, fr[i]["ug"], ug1))
        check(close(fr[i]["ul"], ul1, 1e-9, 1e-9),
              "schillerNaumann step %d: UL %g != replay %g" % (i, fr[i]["ul"], ul1))
        P = mG * fr[i]["ug"] + mL * fr[i]["ul"]
        check(abs(P - P0) <= 1e-12 * max(abs(P0), 1.0),
              "schillerNaumann step %d: mixture momentum not invariant" % i)
else:
    check(False, "schillerNaumann arm: %d frames != steps+1" % len(fr))

# ---- refusals --------------------------------------------------------------
refusals = [
    (BASE.replace("rhoG 1.2;", "rhoG -1.2;"), "rhoG and rhoL must be finite and positive", "neg_rhoG"),
    (BASE.replace("rhoL 998.0;", "rhoL 0.0;"), "rhoG and rhoL must be finite and positive", "zero_rhoL"),
    (BASE.replace("muG  1.8e-5;", "muG  -1;"), "muG and muL must be finite and positive", "neg_muG"),
    (BASE.replace("muL  1.0e-3;", "muL  0.0;"), "muG and muL must be finite and positive", "zero_muL"),
    (BASE.replace("alphaG0 0.4;", "alphaG0 1.5;"), "alphaG0 must lie in [0, 1]", "alpha_gt1"),
    (BASE.replace("alphaG0 0.4;", "alphaG0 -0.2;"), "alphaG0 must lie in [0, 1]", "alpha_lt0"),
    (BASE.replace("steps 200;", "steps 0;"), "run steps must be at least 1", "zero_steps"),
    (BASE.replace("dt 0.001;", "dt 0;"), "dt must be finite and positive", "zero_dt"),
    (BASE.replace("model constant;", "model stokes;"), "drag model 'stokes' is unknown", "bad_model"),
    (BASE.replace("model constant;\n    beta  1200.0;", "model schillerNaumann;"),
     "needs a declared positive", "sn_no_diameter"),
]
for text, needle, tag in refusals:
    r, _ = run(BINARY, text, TMP, "refuse_" + tag)
    check(r.returncode != 0 and needle in (r.stdout + r.stderr),
          "refusal '%s' not raised loudly (rc=%d)" % (needle, r.returncode))

if fails:
    print("eulerEuler FAILED %d checks:" % len(fails))
    for f in fails[:30]:
        print("  -", f)
    sys.exit(1)

print("eulerEuler OK: uniform identity 1e-14, implicit-pair per-step closed form "
      "1e-9 with momentum invariant 1e-12, terminal slip = alphaG alphaL drho g/K, "
      "Schiller-Naumann replay, alpha in [0,1] and mass conserved every frame, "
      "10 loud refusals")
