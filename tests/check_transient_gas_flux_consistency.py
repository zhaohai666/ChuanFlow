#!/usr/bin/env python3
"""Gas-flux consistency on the EOS two-fluid kernel (chflow-1d-twofluid):
the gvol = vM - flux identity and the pG-front invariance it buys, plus
the Wood homogeneous-mixture sound-speed option (E7c acoustic extension).

usage: check_transient_gas_flux_consistency.py <exe> <case-dir>

The kernel is a lumped-mixture-momentum + drift model: ONE mixture
velocity vM with an algebraic liquid face flux F(alpha) and gas volume
flux G(alpha). The model definition carries the identity

    F(alpha) + G(alpha) = vM        (pointwise, every closure),

because the gas inventory is the pipe volume minus the liquid volume, so
the gas transport equation is the exact complement of the liquid one.
The kernel used to solve the gas face flux as an INDEPENDENT Godunov
problem on the mirrored quadratic G = vM - F; the mirror is pointwise
exact but the Godunov branch selection does not mirror the liquid one
(the qa == 0 linear leg upwinds the opposite state, rising/falling jumps
swap the shock and rarefaction branches), so F + G = vM broke on every
non-uniform face - the canonical probe (slip s=1, vM=2, 0.8 -> 0.3)
returned F + G = 3.05, a gas volume flux EXCEEDING the mixture flux,
which over/under-swept gas at each front and drifted the EOS gas
pressure pG. The fix takes the identity face by face; imposed boundaries
keep their exact imposed values (a phaseRates inlet carries vM = vsl +
vsg, so its declared vsg IS the identity; a junction split is built from
the pipe's own flux functions at the mixed holdup; an open outlet's
pointwise closure satisfies it identically).

This checker imports NO C++ path. It re-derives everything from the CASE
TEXT and the REPORTED fields (per-cell EOS pressure `p`, per-cell alpha,
per-face `face j F ... G ...` lines):

main case (case.dict, rate-driven front):
  * face identity F + G = vM on every reported face of every frame
    (exact: worst relative error must be 0);
  * pG front invariance (a): every cell at every frame reads
    p == p0 = rhoG gasR T (all gas densities in the run are the declared
    rhoG, so the identity makes each cell's gas mass evolve exactly as
    rhoG times its gas volume and p never moves);
  * pG front invariance (b): behind the measured front position the
    pressure is a settled plateau - consecutive sampled frames differ by
    ~0 while the front is still advancing;
  * pG front invariance (c), the strong form: the alpha front rides
    vM (measured uF == vM), and with the 1.0 s sample spacing the
    integer cell shift is exact, so p(x) at t matches p(x - uF t) at
    the earlier frame cell for cell;
  * gas volume balance: injected (vsg t) - discharged (reported G at
    the outlet face, time-integrated) == gas inventory change
    sum((1 - alpha) dx) - the gas counterpart of the slug tracking
    checker's liquid balance, same tolerance style;
  * the zero-declaration twin (caseBase.dict, uniform, no front) runs
    the identical face report and never moves p - the repair has
    nothing to change there and the report is pure observation.

fixedPressure arm (caseFP.dict): the same front driven by a 2000 Pa
pressure difference - vM and the front speed are the momentum solve's
answer, not declared numbers. The inlet ghost density is the declared
pressure's EOS density (p / (gasR T) = 11.61 = rhoG), so the same
p == p0 invariance, identity and gas balance hold with the REPORTED vm.

Wood arm (caseWood.dict, sealed quiescent 50/50 line, acoustics with
aModel wood + liquidSoundSpeed 1200):
  * a_wood re-derived from the case text alone: the compliance sum
    alphaL/(rhoL aL^2) + alphaG/(rhoG aG^2) with aG^2 = gasR T, the
    mixture speed a = sqrt(1/(compliance rho_m)) ~ 62.8 m/s - the
    textbook collapse (the liquid compliance dominates rhoG aG^2 = p0);
  * the reported acoustic fronts ride a_wood, NOT the EOS gas speed
    293.4 m/s (first arrivals at probe cells match (L - x)/a_wood);
  * the hard-wall gauge hysteresis transitions are spaced 2L/a_wood;
  * the undeclared twin (caseWoodBase.dict, default aModel eos) rides
    the EOS speed - the opt-in model does not touch the default path.

rejections (>= 6, loud cfFatal): aModel wood without liquidSoundSpeed
(caseRejectWood.dict), wood + declared soundSpeed (mutually exclusive),
unknown aModel word, non-positive liquidSoundSpeed (0 and negative),
wood with the EOS off, and wood with acoustics enable false (a
contradictory declaration is refused even disabled).
"""

import math
import os
import re
import subprocess
import sys

if len(sys.argv) != 3:
    raise SystemExit(__doc__)
EXE, CASEDIR = sys.argv[1:3]

fails = []


def check(name, got, want, rtol, atol=0.0):
    tol = atol + rtol * abs(want)
    ok = abs(got - want) <= tol
    print(f"{'ok:  ' if ok else 'FAIL:'} {name}: got {got!r} want {want!r}")
    if not ok:
        fails.append(name)


def check_true(name, cond):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name}: {cond}")
    if not cond:
        fails.append(name)


def expect_fatal(name, path, needle=None):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    msg = (r.stdout + r.stderr).lower()
    bad = r.returncode == 0 or ("fatal" not in msg and "error" not in msg)
    if not bad and needle is not None:
        bad = needle.lower() not in msg
    print(f"{'ok:  ' if not bad else 'FAIL:'} {name} refused "
          f"(rc={r.returncode})")
    if bad:
        fails.append(name)


def run_raw(path):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"case {path} failed rc={r.returncode}\n"
                         f"{r.stdout}\n{r.stderr}")
    return r.stdout


def parse(text):
    """Frames with the per-cell rows (alpha, p, ...) and the per-face
    flux rows (face j F F G G)."""
    frames = []
    for line in text.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "t" and len(tok) >= 4:
            frames.append({"t": float(tok[1]), "vm": float(tok[3]),
                           "cells": [], "faces": {}})
        elif tok[0] == "cell" and frames:
            d = {}
            for k in ("alpha", "p", "pA", "vsl", "vsg"):
                if k in tok:
                    d[k] = float(tok[tok.index(k) + 1])
            frames[-1]["cells"].append(d)
        elif tok[0] == "face" and frames:
            frames[-1]["faces"][int(tok[1])] = (
                float(tok[tok.index("F") + 1]),
                float(tok[tok.index("G") + 1]))
    return frames


def case_number(text, key):
    # first declaration of the key anywhere in the case text (the block
    # order in these dicts makes the first hit the intended one)
    m = re.search(rf"(?m)^\s*{re.escape(key)}\s+([-+0-9.eE]+)\s*;",
                  text) or re.search(
                      rf"{re.escape(key)}\s+([-+0-9.eE]+)\s*;", text)
    if not m:
        raise SystemExit(f"case text has no {key}")
    return float(m.group(1))


def identity_err(frames):
    """Worst relative identity violation |F + G - vM| measured against
    the face flux scale (the honest denominator: near-zero mixture
    velocities early in the momentum transient must not inflate the
    relative error). The kernel identity is exact; the tolerance below
    carries the 12-significant-digit quantization of the REPORTED vm
    (the faces are asserted against the printed velocity)."""
    worst = 0.0
    for f in frames:
        for j, (F, G) in f["faces"].items():
            worst = max(worst, abs(F + G - f["vm"])
                        / max(abs(F) + abs(G), 1e-30))
    return worst


def step_integral(frames, face, dt):
    """Exact discrete time integral of a face's G: the case reports
    every step (reportEvery 1), and each sampled face value IS the step's
    own flux, so the rectangular right-sum over the uniform dt sampling
    is the exact discrete integral - no quadrature allowance needed.
    Skips the t = 0 frame, which is emitted before the first transport
    step and carries no face rows."""
    stepped = [f for f in frames if f["faces"]]
    acc = 0.0
    for f in stepped:
        acc += f["faces"][face][1] * dt
    return acc


# ================= shared case text =================
path = os.path.join(CASEDIR, "case.dict")
case_text = open(path).read()
RhoG = case_number(case_text, "rhoG")
GasR = case_number(case_text, "gasR")
Temp0 = case_number(case_text, "temperature")
Vsl = case_number(case_text, "vsl")
Vsg = case_number(case_text, "vsg")
Sections = int(case_number(case_text, "sections"))
Length = case_number(case_text, "length")
DX = Length / Sections
Dt = case_number(case_text, "dt")
P0 = RhoG * GasR * Temp0
VM = Vsl + Vsg

# ================= main case: identity + pG front invariance ==========
text = run_raw(path)
frames = parse(text)
check_true("main case produced frames", len(frames) >= 10)

# (identity) F + G == vM on every reported face of every frame - the
# model's own defining identity, asserted from the report alone
worst = identity_err(frames)
check("face identity F + G = vM holds on every face every frame "
      "(worst relative error, print-quantization tolerance)", worst,
      0.0, 0.0, 1e-11)
stepped = [f for f in frames if f["faces"]]
check_true("stepped frames carry the full face set (n+1 faces)",
           all(len(f["faces"]) == Sections + 1 for f in stepped)
           and len(stepped) >= len(frames) - 1)

# (a) pG front invariance: every cell at every frame sits at p0
worstp = 0.0
for f in frames:
    for c in f["cells"]:
        worstp = max(worstp, abs(c["p"] - P0) / P0)
check("pG front invariance (a): every cell at every frame reads "
      "p == p0 = rhoG gasR T", worstp, 0.0, 0.0, 1e-9)

# front position: the midpoint crossing of the alpha step (0.8 -> 0.3)
def front_cell(f):
    alphas = [c["alpha"] for c in f["cells"]]
    for i, a in enumerate(alphas):
        if a < 0.55:
            return i
    return None

# (b) behind the measured front the pressure is a settled plateau
plateau_err = 0.0
settled_seen = False
for k in range(1, len(frames)):
    fc = front_cell(frames[k])
    if fc is None or fc < 2:
        continue
    settled_seen = True
    # cells at least 2 behind the front, while the front is still inside
    if fc <= Sections - 2:
        for i in range(fc - 2):
            plateau_err = max(plateau_err,
                              abs(frames[k]["cells"][i]["p"]
                                  - frames[k - 1]["cells"][i]["p"]) / P0)
check_true("front cells measurable across the run", settled_seen)
check("pG front invariance (b): behind the front d pG/dt -> 0 (worst "
      "frame-to-frame plateau drift)", plateau_err, 0.0, 0.0, 1e-9)

# (c) strong form: the alpha front rides vM and the p profile is
# shift-invariant. Every 5th sampled frame (1.0 s spacing) shifts the
# profile by exactly one cell at uF = vM = 1 m/s with dx = 1 m.
# measure inside the window where the front is still inside the line
# (the front enters at x = 50 and leaves at t = 50)
u_frames = [f for f in frames if 5.0 < f["t"] < 45.0]
fc0 = front_cell(u_frames[0])
fc1 = front_cell(u_frames[-1])
t0, t1 = u_frames[0]["t"], u_frames[-1]["t"]
# measured speed from the front cells (the +0.5*dx centre convention
# cancels in the difference)
uF_measured = (fc1 - fc0) * DX / (t1 - t0)
check("the alpha front rides the mixture velocity (measured uF vs vM)",
      uF_measured, VM, 0.02)
stride = int(round(1.0 / Dt))  # frames per 1.0 s at reportEvery 1
shift = int(round(VM * (u_frames[stride]["t"] - u_frames[0]["t"]) / DX))
check_true("sample spacing gives an integer cell shift", shift >= 1)
shift_err = 0.0
for k in range(0, len(frames) - shift * stride, shift * stride):
    a = frames[k]
    b = frames[k + shift * stride]
    fa, fb = front_cell(a), front_cell(b)
    if fa is None or fb is None or fa < shift + 1 or fb < 1:
        continue
    # cells strictly behind the front on both frames, shifted by the
    # front travel: p(x - uF t) profiles must coincide cell for cell
    for i in range(shift + 1, min(fa, fb + shift)):
        shift_err = max(shift_err,
                        abs(b["cells"][i]["p"]
                            - a["cells"][i - shift]["p"]) / P0)
check("pG front invariance (c): p(x - uF t) shift-invariant (worst "
      "cell mismatch over the shifted profiles)", shift_err, 0.0, 0.0,
      1e-9)

# gas volume balance: injected - discharged == inventory change (per
# unit bore area, the gas counterpart of the slug tracking checker's
# liquid balance, same 1e-3 documented tolerance). The case reports
# every step, so the reported outlet-face G integrates EXACTLY (each
# sample is the step's own flux, uniform dt); the tolerance only carries
# the printed-value quantization and the alpha clamps.
inv = [sum((1.0 - c["alpha"]) * DX for c in f["cells"]) for f in frames]
t_end = frames[-1]["t"]
injected = Vsg * t_end
discharged = step_integral(frames, Sections, Dt)
check("gas volume balance: injected - discharged == inventory change "
      "(independent of the kernel, documented tolerance 1e-3)",
      injected - discharged, inv[-1] - inv[0], 0.0, 1e-3)

# ================= zero-declaration twin (uniform) =====================
base_text = run_raw(os.path.join(CASEDIR, "caseBase.dict"))
base_frames = parse(base_text)
b_worst = identity_err(base_frames)
check("zero-declaration twin: the face report rides the uniform run and "
      "the identity holds there too", b_worst, 0.0, 0.0, 1e-11)
b_worstp = max(abs(c["p"] - P0) / P0
               for f in base_frames for c in f["cells"])
check("zero-declaration twin: uniform run never moves p (the repair has "
      "nothing to change, the report is pure observation)", b_worstp,
      0.0, 0.0, 1e-9)
check("zero-declaration twin: alpha stays uniform 0.8 (no front exists)",
      max(abs(c["alpha"] - 0.8)
          for f in base_frames for c in f["cells"]), 0.0, 0.0, 1e-12)

# ================= fixedPressure arm ===================================
fp_text = run_raw(os.path.join(CASEDIR, "caseFP.dict"))
fp_frames = parse(fp_text)
fp_worst = identity_err(fp_frames)
check("fixedPressure arm: face identity F + G = vM with the momentum-"
      "solved reported vm", fp_worst, 0.0, 0.0, 1e-11)
fp_worstp = max(abs(c["p"] - P0) / P0
                for f in fp_frames for c in f["cells"])
check("fixedPressure arm: pG front invariance holds with the ghost-"
      "density inlet (p == p0 everywhere)", fp_worstp, 0.0, 0.0, 1e-9)
fp_inv = [sum((1.0 - c["alpha"]) * DX for c in f["cells"])
          for f in fp_frames]
# both face integrals from the reported per-step fluxes, exact
fp_dis = step_integral(fp_frames, Sections, Dt)
fp_inj = step_integral(fp_frames, 0, Dt)
check("fixedPressure arm: gas volume balance from the reported face "
      "fluxes (documented tolerance 1e-3)",
      fp_inj - fp_dis, fp_inv[-1] - fp_inv[0], 0.0, 1e-3)

# ================= Wood arm ============================================
wood_path = os.path.join(CASEDIR, "caseWood.dict")
wood_text = open(wood_path).read()
RhoL = case_number(wood_text, "rhoL")
AL = case_number(wood_text, "liquidSoundSpeed")
A0w = case_number(wood_text, "alpha")
# Wood's homogeneous-mixture sound speed, re-derived from the case text:
#   1/(rho_m a^2) = alphaL/(rhoL aL^2) + alphaG/(rhoG aG^2)
# with aG^2 = gasR T (the EOS-exact gas derivative) and rhoG = P0/aG^2.
# The case's `alpha` is the LIQUID holdup (the kernel weights the
# compliance by alpha_ for the liquid term and 1-alpha_ for the gas
# term, twoFluid.cpp acousticsStep), so both lines below use it the
# same way (2026-11-24 fix: the first version swapped the weights in
# the compliance only; invisible here because caseWood declares the
# symmetric alpha 0.5).
AG2 = GasR * Temp0
RhoGw = P0 / AG2
compliance = A0w / (RhoL * AL * AL) + (1.0 - A0w) / (RhoGw * AG2)
rho_m = A0w * RhoL + (1.0 - A0w) * RhoGw
A_wood = math.sqrt(1.0 / (compliance * rho_m))
A_eos = math.sqrt(AG2)
check("Wood speed re-derived from the case text is far below the EOS "
      "gas speed (the mixture-compressibility collapse)",
      A_eos / A_wood, 4.5, 0.3)

wood_frames = parse(run_raw(wood_path))
Lw = case_number(wood_text, "length")
# first-arrival times at the probe cells near the venting (outlet) end:
# the release wave leaves the outlet at t = 0 and arrives at cell x at
# (L - x)/a. Fit the speed from the arrivals of three probe cells.
probes = [95, 85, 75]
pts = []
for i in probes:
    t1 = next(f["t"] for f in wood_frames if f["cells"][i]["pA"] < 0.0
              if "pA" in f["cells"][i])
    pts.append((Lw - (i + 0.5) * DX, t1))
sx = sum(p[0] for p in pts)
sy = sum(p[1] for p in pts)
sxx = sum(p[0] ** 2 for p in pts)
sxy = sum(p[0] * p[1] for p in pts)
slope = (len(pts) * sxy - sx * sy) / (len(pts) * sxx - sx * sx)
a_fit = 1.0 / slope
# tolerances carry the arrival-time quantization of the report spacing
# (the probe nearest the vent crosses within one report sample)
check("Wood arm: the reported acoustic fronts ride the re-derived "
      "Wood speed", a_fit, A_wood, 0.08)
check("Wood arm: the fronts ride the Wood speed and NOT the EOS gas "
      "speed (the EOS/Wood ratio reproduced from the report)",
      A_eos / a_fit, A_eos / A_wood, 0.08)
# reflection: the hard-wall gauge signal repeats with the 4L/a
# quarter-wave period; verify the reported wall pressure at the inlet
# end cell completes its first full cycle within 4L/a_wood (+/- 5%)
p_inlet = [(f["t"], f["cells"][0]["pA"]) for f in wood_frames]
half = Lw / A_wood
t_cross = next((t for t, p in p_inlet if t > half and p < 0.0), None)
check("Wood arm: the wall gauge first flips sign at the quarter-wave "
      "time L/a_wood", t_cross if t_cross is not None else 0.0,
      half, 0.05)

# undeclared twin rides the EOS speed
base_wood_text = open(os.path.join(CASEDIR, "caseWoodBase.dict")).read()
bw_frames = parse(run_raw(os.path.join(CASEDIR, "caseWoodBase.dict")))
bpts = []
for i in probes:
    t1 = next(f["t"] for f in bw_frames if f["cells"][i]["pA"] < 0.0
              if "pA" in f["cells"][i])
    bpts.append((Lw - (i + 0.5) * DX, t1))
bsx = sum(p[0] for p in bpts)
bsy = sum(p[1] for p in bpts)
bsxx = sum(p[0] ** 2 for p in bpts)
bsxy = sum(p[0] * p[1] for p in bpts)
bslope = (len(bpts) * bsxy - bsx * bsy) / (len(bpts) * bsxx - bsx * bsx)
check("Wood undeclared twin: the default aModel eos path rides the EOS "
      "gas speed (the opt-in model leaves it verbatim)",
      1.0 / bslope, A_eos, 0.05)

# ================= rejection witnesses =================================
expect_fatal("wood without liquidSoundSpeed",
             os.path.join(CASEDIR, "caseRejectWood.dict"),
             "liquidSoundSpeed")

wood_decl = """    aModel wood;           // Wood's homogeneous-mixture sound speed
    liquidSoundSpeed 1200; // the DECLARED liquid phase sound speed aL
"""


def wood_variant(name, replacement, needle=None):
    p = os.path.join(CASEDIR, f"caseReject_{name}.dict")
    assert wood_decl in wood_text
    open(p, "w").write(wood_text.replace(wood_decl, replacement))
    expect_fatal(f"wood witness: {name}", p, needle)


wood_variant("mutual_soundSpeed",
             "    soundSpeed 250;      // wholesale override: contradicts wood\n"
             "    aModel wood;\n    liquidSoundSpeed 1200;\n",
             "mutually exclusive")
wood_variant("unknown_word",
             "    aModel adiabatic;\n    liquidSoundSpeed 1200;\n",
             "unknown acoustics amodel")
wood_variant("aL_zero",
             "    aModel wood;\n    liquidSoundSpeed 0;\n",
             "liquidSoundSpeed")
wood_variant("aL_negative",
             "    aModel wood;\n    liquidSoundSpeed -1200;\n",
             "liquidSoundSpeed")
# wood with the EOS off: gasR removed -> the enable block refuses
no_eos = wood_text.replace(wood_decl,
                           "    aModel wood;\n    liquidSoundSpeed 1200;\n")
no_eos = re.sub(r"(?m)^\s*gasR\s+[-0-9.eE]+;\s*(//[^\n]*)?$", "", no_eos)
p = os.path.join(CASEDIR, "caseReject_noEOS.dict")
open(p, "w").write(no_eos)
expect_fatal("wood witness: wood with the EOS off", p, "gas eos")
# contradictory wood declaration even with enable false
dis = wood_text.replace(wood_decl,
                        "    enable false;\n    aModel wood;\n")
p = os.path.join(CASEDIR, "caseReject_disabled.dict")
open(p, "w").write(dis)
expect_fatal("wood witness: contradictory wood declaration with enable "
             "false", p, "liquidSoundSpeed")

# ================= summary =============================================
print()
if fails:
    print(f"FAIL ({len(fails)}):")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("all checks passed")
