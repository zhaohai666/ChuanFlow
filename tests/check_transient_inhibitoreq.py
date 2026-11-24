#!/usr/bin/env python3
"""O4 (5.163) rigorous inhibitor phase equilibrium for the 1D two-fluid
kernel (chflow-1d-twofluid).

usage: check_transient_inhibitoreq.py <exe> <case.dict> <caseBase.dict>
        <caseDisabled.dict> <casePure.dict> <caseFlow.dict> <caseReject.dict>

inhibitorEquilibrium{} replaces the Hammerschmidt offset with the rigorous
root of the declared curve at the water activity of the transported
inhibitor: aW = 1/(1 + activityK w), T_eq = curveB/(ln p + ln aW - curveA).

Sealed case: an inclined line whose geothermal wall falls across the
rigorous equilibrium line, so the cool half grows hydrate and the warm half
dissociates it (spatial dTsub sign flip). The checker replicates the full
per-cell recursion (activity root, implicit growth/decay, 1:1 water move,
stoichiometric gas sink, implicit latent thermal step) at 1e-9 and checks
the reported Teq/dTsub columns against an independent formula on the
reported p/inh/T columns.

Flowing case: an inhibitor front (inhibitorFraction 0.4) advects through a
hydrate-bearing line: growth ahead of the front, dissociation behind it
(temporal flip). Full coupled replication (Godunov transport of alpha /
water / inhibitor / hydrate / gas mass + rigorous kinetics + implicit
enthalpy solve) at 1e-9.

Invariance: enable false is byte-identical to the undeclared twin; a pure
line (inhibitor 0) is trajectory-identical to its block-stripped twin (the
declaration line, which names the block, is filtered); the undeclared run
reports no Teq token at all. K monotonicity siblings and a basis-water
sibling (with the documented c/water clamp exercised) close the physics;
eight illegal declarations are refused.
"""

import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 8:
    raise SystemExit(__doc__)
EXE, CASE, CASE_BASE, CASE_DISABLED, CASE_PURE, CASE_FLOW, CASE_REJECT = \
    sys.argv[1:8]

fails = []


def check(name, got, want, rtol, atol=0.0):
    tol = atol + rtol * abs(want)
    ok = abs(got - want) <= tol
    print(f"{'ok:  ' if ok else 'FAIL:'} {name}: got {got!r} want {want!r}")
    if not ok:
        fails.append(name)


def ok(name):
    print(f"ok:   {name}")


def parse(text):
    frames = []
    for line in text.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "t" and len(tok) >= 6:
            frames.append({"t": float(tok[1]), "vm": float(tok[3]),
                           "dp": float(tok[5]), "cells": []})
        elif tok[0] == "cell" and frames:
            c = {}
            key = None
            for w in tok[2:]:
                if key is None:
                    key = w
                else:
                    c[key] = float(w)
                    key = None
            frames[-1]["cells"].append(c)
    return frames


def run(path):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"case {path} failed rc={r.returncode}\n"
                         f"{r.stdout}\n{r.stderr}")
    return r.stdout, parse(r.stdout)


def run_text(text):
    with tempfile.NamedTemporaryFile("w", suffix=".dict", delete=False) as f:
        f.write(text)
        path = f.name
    try:
        return run(path)
    finally:
        os.unlink(path)


def expect_fatal(name, text):
    with tempfile.NamedTemporaryFile("w", suffix=".dict", delete=False) as f:
        f.write(text)
        path = f.name
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    os.unlink(path)
    bad = r.returncode == 0
    msg = (r.stdout + r.stderr).lower()
    if not bad and "fatal" not in msg and "error" not in msg:
        bad = True
    print(f"{'ok:  ' if not bad else 'FAIL:'} {name} refused "
          f"(rc={r.returncode})")
    if bad:
        fails.append(name)


def godunov(qa, qb, qc, aL, aR):
    # port of the kernel's exact face flux for F(a) = qa a^2 + qb a + qc
    fL = qc + aL * (qb + qa * aL)
    if aL == aR:
        return fL
    fR = qc + aR * (qb + qa * aR)
    shock = (aL > aR) if qa > 0 else (aL < aR)
    if shock:
        s = (fL - fR) / (aL - aR)
        return fL if s >= 0 else fR
    if qa == 0:
        return fL if qb >= 0 else fR
    sonic = -qb / (2.0 * qa)
    lo, hi = min(aL, aR), max(aL, aR)
    a = min(max(sonic, lo), hi)
    return qc + a * (qb + qa * a)


def solve_dense(mat, rhs):
    # port of src/1d/src/DenseSolve.h (Gauss elimination with partial
    # pivoting + back substitution)
    n = len(rhs)
    for k in range(n):
        piv, best = k, abs(mat[k][k])
        for i in range(k + 1, n):
            if abs(mat[i][k]) > best:
                best, piv = abs(mat[i][k]), i
        if piv != k:
            mat[k], mat[piv] = mat[piv], mat[k]
            rhs[k], rhs[piv] = rhs[piv], rhs[k]
        for i in range(k + 1, n):
            m = mat[i][k] / mat[k][k]
            if m == 0:
                continue
            for j in range(k, n):
                mat[i][j] -= m * mat[k][j]
            rhs[i] -= m * rhs[k]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = rhs[i] - sum(mat[i][j] * x[j] for j in range(i + 1, n))
        x[i] = s / mat[i][i]
    return x


# =====================================================================
# shared case constants (mirrors of the case files)
# =====================================================================
RHO_L, RHO_G = 1000.0, 200.0
CL, CG, GAS_R = 4180.0, 2150.0, 287.0
CURVE_A, CURVE_B = 35.2, -5684.0
K_G, K_D = 0.01, 0.01
GAS_CONTENT, LATENT = 122.0, 4.2e8
KW = 50.0
DI = 0.02
AREA = math.pi * DI * DI / 4.0


def teq_rig(p, w, ak):
    aW = 1.0 / (1.0 + ak * w)
    return CURVE_B / (math.log(p) + math.log(aW) - CURVE_A)


CASE_TEXT = open(CASE).read()

# =====================================================================
# sealed line: full per-cell replication of the rigorous recursion
# =====================================================================
N, DT, DX = 10, 0.5, 2.0
A0, W0, T0 = 0.5, 0.5, 300.0
H0, C0, AK = 0.04, 0.3, 8.0
STEPS_PER_FRAME = 40
# geothermal wall: depth_i = inletDepth - (i + 0.5) dx sin(30 deg)
TW = [172.6 + 1.2 * (100.0 - (i + 0.5)) for i in range(N)]

stdout_sealed, frames = run(CASE)
if "twoFluid inhibitorEquilibrium: activityK 8 basis liquid" \
        not in stdout_sealed:
    fails.append("declaration line missing")
else:
    ok("declaration line names activityK 8 basis liquid")


def replicate_sealed(basis_water=False, ak=AK, c0=C0):
    """Exact port of the sealed rigorous recursion; returns the per-frame
    replicated states."""
    alpha = [A0] * N
    water = [W0] * N
    hyd = [H0] * N
    temp = [T0] * N
    inh = [c0] * N
    mg = [RHO_G * (1.0 - A0) * DX] * N
    steps_done = 0
    out = []
    for k in range(len(frames)):
        want_steps = 0 if k == 0 else k * STEPS_PER_FRAME
        while steps_done < want_steps:
            for i in range(N):
                ag = 1.0 - alpha[i]
                p_i = mg[i] * GAS_R * temp[i] / (ag * DX)
                if basis_water:
                    w = min(inh[i] / water[i], 1.0) if water[i] > 0 else 0.0
                else:
                    w = inh[i]
                dT = teq_rig(p_i, w, ak) - temp[i]
                total = water[i] + hyd[i]
                if dT >= 0:
                    r = K_G * dT
                    h_new = (hyd[i] / DT + r * total) / (1.0 / DT + r)
                    trace = 1e-9 * RHO_G * DX
                    dh_max = (mg[i] - trace) / (GAS_CONTENT * alpha[i] * DX)
                    h_new = min(h_new, hyd[i] + max(dh_max, 0.0))
                else:
                    h_new = hyd[i] / (1.0 + DT * K_D * (-dT))
                h_new = min(max(h_new, 0.0), total)
                dh = h_new - hyd[i]
                water[i] = min(max(total - h_new, 0.0), 1.0)
                hyd[i] = h_new
                mg[i] -= GAS_CONTENT * alpha[i] * dh * DX
                cap = AREA * DX * (alpha[i] * RHO_L * CL
                                   + ag * (mg[i] / (ag * DX)) * CG)
                kw = KW * DX
                temp[i] = (cap / DT * temp[i] + kw * TW[i]
                           + LATENT * alpha[i] * dh * AREA * DX / DT) \
                          / (cap / DT + kw)
            steps_done += 1
        out.append((list(alpha), list(water), list(hyd), list(temp),
                    list(inh), list(mg)))
    return out


def compare_sealed(frames_ref, rep):
    bad = None
    for k, fr in enumerate(frames_ref):
        alpha, water, hyd, temp, inh, mg = rep[k]
        for i in range(N):
            got = fr["cells"][i]
            ag = 1.0 - alpha[i]
            want_p = mg[i] * GAS_R * temp[i] / (ag * DX)
            want_teq = teq_rig(want_p, inh[i], AK)
            if (abs(got["water"] - water[i]) > 1e-9 + 1e-9 * water[i]
                    or abs(got["h"] - hyd[i]) > 1e-9 + 1e-9 * hyd[i]
                    or abs(got["T"] - temp[i]) > 1e-9 + 1e-9 * temp[i]
                    or abs(got["p"] - want_p) > 1e-9 * want_p
                    or abs(got["inh"] - inh[i]) > 1e-9
                    or abs(got["Teq"] - want_teq) > 1e-9 + 1e-9 * want_teq
                    or abs(got["dTsub"] - (want_teq - temp[i])) > 1e-8):
                bad = (k, fr["t"], i, got)
                break
        if bad:
            break
    return bad


rep = replicate_sealed()
bad = compare_sealed(frames, rep)
if bad:
    fails.append(f"sealed replication mismatch frame {bad[0]} t={bad[1]} "
                 f"cell {bad[2]}: {bad[3]}")
else:
    for k in range(0, len(frames), 2):
        ok(f"sealed rigorous recursion replicated at t={frames[k]['t']} "
           f"(w/h/T/p/inh/Teq/dTsub, 1e-9)")

# reported Teq/dTsub columns vs the independent formula on the reported p
worst_teq = worst_dts = 0.0
for fr in frames:
    for got in fr["cells"]:
        want = teq_rig(got["p"], got["inh"], AK)
        worst_teq = max(worst_teq, abs(got["Teq"] - want))
        worst_dts = max(worst_dts, abs(got["dTsub"] - (got["Teq"] - got["T"])))
check("Teq column = curveB/(ln p + ln aW - curveA) from reported columns",
      worst_teq, 0.0, 0.0, 1e-8)
check("dTsub column = Teq - T from reported columns", worst_dts, 0.0, 0.0,
      1e-8)

# spatial sign flip: cool half grows, warm half dissociates
fr1 = frames[1]
grow = [i for i, c in enumerate(fr1["cells"]) if c["dTsub"] > 0.3]
shrink = [i for i, c in enumerate(fr1["cells"]) if c["dTsub"] < -0.3]
if grow and shrink:
    ok(f"frame t=20: {len(shrink)} cells dTsub < -0.3, {len(grow)} cells "
       f"dTsub > +0.3 (spatial sign flip)")
else:
    fails.append(f"no spatial sign flip at t=20: grow={grow} "
                 f"shrink={shrink}")
flip_ok = True
for k in range(2, len(frames)):
    for i in range(N):
        d = frames[k]["cells"][i]["h"] - frames[k - 1]["cells"][i]["h"]
        ds = frames[k]["cells"][i]["dTsub"]
        if ds > 0.3 and d <= 0.0:
            flip_ok = False
        if ds < -0.3 and frames[k - 1]["cells"][i]["h"] > 1e-6 and d >= 1e-15:
            flip_ok = False
if flip_ok:
    ok("every frame (from t=40, past the transient): dTsub > +0.3 cells "
       "gain h, dTsub < -0.3 cells with hydrate lose h")
else:
    fails.append("dTsub sign vs delta-h inconsistent somewhere")

# conservation from reported columns only: sealed 1:1 conversion keeps
# water + hydrate at the initial total
cons = max(abs(c["water"] + c["h"] - (W0 + H0))
           for fr in frames for c in fr["cells"])
check("per-cell water + h = 0.54 (1:1 conversion, no transport)", cons, 0.0,
      0.0, 1e-9)
# gas inventory: emitted EOS mass vs the stoichiometric sink on h
mg_emit = sum(c["p"] * (1.0 - c["alpha"]) * DX / (GAS_R * c["T"])
              for c in frames[-1]["cells"])
mg_want = sum(RHO_G * (1.0 - A0) * DX
              - GAS_CONTENT * A0 * (c["h"] - H0) * DX
              for c in frames[-1]["cells"])
check("gas inventory = initial - gasContent * alpha * dh * dx",
      mg_emit, mg_want, 1e-6)
# latent witness: while growth runs (t=100, past the initial transient)
# the growth cell sits well above its own wall; the fully inhibited K40
# sibling (dissociation everywhere, no latent release) lands ON the wall.
stdout_pure, frames_pure = run(CASE_PURE)

# =====================================================================
# invariance twins
# =====================================================================
stdout_base, _ = run(CASE_BASE)
stdout_dis, _ = run(CASE_DISABLED)
if stdout_dis == stdout_base:
    ok("enable false is byte-identical to the undeclared twin")
else:
    fails.append("disabled twin differs from the undeclared twin")
if "Teq" not in stdout_base and "inhibitorEquilibrium" not in stdout_base:
    ok("undeclared run reports no Teq / inhibitorEquilibrium token")
else:
    fails.append("undeclared run leaked the rigorous tokens")


def strip_block(text):
    return re.sub(r"inhibitorEquilibrium\s*\{[^}]*\}\s*", "", text)


def filt(text):
    # the declaration line naming the block and the gated Teq/dTsub columns
    # are expected output of the declared path, not trajectory: strip them,
    # everything else must match byte for byte
    out = []
    for l in text.splitlines():
        if "inhibitorEquilibrium" in l:
            continue
        out.append(re.sub(r" Teq \S+ dTsub \S+", "", l))
    return "\n".join(out)


stdout_stripped, _ = run_text(strip_block(open(CASE_PURE).read()))
if filt(stdout_pure) == filt(stdout_stripped):
    ok("pure line (inhibitor 0) trajectory byte-identical to the "
       "block-stripped twin")
else:
    fails.append("pure twin differs from the stripped run")
# at w = 0 the rigorous root IS the declared base line
worst = max(abs(c["Teq"] - CURVE_B / (math.log(c["p"]) - CURVE_A))
            for fr in frames_pure for c in fr["cells"])
check("pure line Teq = declared curve verbatim (log(1) = 0)", worst, 0.0,
      0.0, 1e-8)

# =====================================================================
# K monotonicity siblings
# =====================================================================
_, frames_k10 = run_text(CASE_TEXT.replace("activityK 8;", "activityK 10;"))
_, frames_k40 = run_text(CASE_TEXT.replace("activityK 8;", "activityK 40;"))
c8, c10, c40 = (frames[-1]["cells"][N - 1], frames_k10[-1]["cells"][N - 1],
                frames_k40[-1]["cells"][N - 1])
if c40["Teq"] < c10["Teq"] < c8["Teq"]:
    ok("final Teq strictly lower at higher activityK")
else:
    fails.append(f"K monotonicity of Teq broken: {c8['Teq']} {c10['Teq']} "
                 f"{c40['Teq']}")
if c40["h"] < c10["h"] < c8["h"]:
    ok("final hydrate inventory strictly lower at higher activityK")
else:
    fails.append(f"K monotonicity of h broken: {c8['h']} {c10['h']} "
                 f"{c40['h']}")
# the temperature witness is mid-run (t=100): once conversion completes
# every cell relaxes back onto its wall regardless of how much hydrate
# formed, so only running growth keeps a cell above the wall; stronger
# inhibition (higher K) means less growth and less latent release
m8, m10, m40 = (frames[5]["cells"][N - 1], frames_k10[5]["cells"][N - 1],
                frames_k40[5]["cells"][N - 1])
if m8["T"] > m10["T"] > m40["T"]:
    ok("t=100 temperature strictly lower at higher activityK (less latent "
       "release)")
else:
    fails.append(f"K monotonicity of T broken: {m8['T']} {m10['T']} "
                 f"{m40['T']}")
warm = m8["T"] - TW[N - 1]
flat = m40["T"] - TW[N - 1]
if warm > 1.0 and abs(flat) < 0.1:
    ok(f"latent release: growth cell at t=100 sits {warm:.2f} K above its "
       f"wall, the dissociating K40 sibling {flat:.2e}")
else:
    fails.append(f"latent witness failed: warm={warm} flat={flat}")

# =====================================================================
# basis water sibling (activityK 2 so growth runs and the c/water ratio
# crosses 1, exercising the documented clamp)
# =====================================================================
BW_AK = 2.0
bw_text = (CASE_TEXT.replace("basis liquid;", "basis water;")
           .replace("activityK 8;", f"activityK {BW_AK};"))
_, frames_bw = run_text(bw_text)
rep_bw = replicate_sealed(basis_water=True, ak=BW_AK)
bad_bw = None
for k, fr in enumerate(frames_bw):
    alpha, water, hyd, temp, inh, mg = rep_bw[k]
    for i in range(N):
        got = fr["cells"][i]
        ag = 1.0 - alpha[i]
        want_p = mg[i] * GAS_R * temp[i] / (ag * DX)
        w = min(inh[i] / water[i], 1.0) if water[i] > 0 else 0.0
        want_teq = teq_rig(want_p, w, BW_AK)
        if (abs(got["h"] - hyd[i]) > 1e-9 + 1e-9 * hyd[i]
                or abs(got["T"] - temp[i]) > 1e-9 + 1e-9 * temp[i]
                or abs(got["water"] - water[i]) > 1e-9 + 1e-9 * water[i]
                or abs(got["Teq"] - want_teq) > 1e-9 + 1e-9 * want_teq):
            bad_bw = (k, i, got)
            break
    if bad_bw:
        break
if bad_bw:
    fails.append(f"basis-water replication mismatch frame {bad_bw[0]} "
                 f"cell {bad_bw[1]}: {bad_bw[2]}")
else:
    ok("basis-water sealed recursion replicated at 1e-9 (all frames)")
clamped = [i for i, c in enumerate(frames_bw[-1]["cells"])
           if c["water"] > 0 and c["inh"] / c["water"] > 1.0]
if clamped:
    ok(f"c/water clamp exercised at t=400 in cells {clamped} (w pinned "
       f"at 1)")
else:
    fails.append("basis-water sibling never exercised the clamp")

# =====================================================================
# flowing line: inhibitor front, growth -> dissociation temporal flip
# =====================================================================
FN, FDT, FDX = 20, 0.2, 0.5
FVM, FUS, VSL, VSG = 0.3, 0.3, 0.05, 0.25
FA0, FW0, FH0, FC0 = 0.5, 0.5, 0.02, 0.0
FT0, FIN, FOUT_P, FOUT_T = 295.0, 290.0, 1.7e7, 290.0
FWALL, FAK, FINH_IN = 295.0, 20.0, 0.4
FSTEPS = 20
stdout_flow, fframes = run(CASE_FLOW)
if "twoFluid inhibitorEquilibrium: activityK 20 basis liquid" \
        not in stdout_flow:
    fails.append("flow declaration line missing")
else:
    ok("flow declaration line names activityK 20 basis liquid")

alpha = [FA0] * FN
water = [FW0] * FN
hyd = [FH0] * FN
temp = [FT0] * FN
inh = [FC0] * FN
mg = [RHO_G * (1.0 - FA0) * FDX] * FN
steps_done = 0
badf = None
worst_te = 0.0
for k, fr in enumerate(fframes):
    want_steps = 0 if k == 0 else k * FSTEPS
    while steps_done < want_steps:
        # ---- transport (wallis, imposed vM) ----
        flux = [VSL]
        for j in range(1, FN):
            flux.append(godunov(FUS, FVM - FUS, 0.0, alpha[j - 1],
                                alpha[j]))
        a_end = alpha[FN - 1]
        flux.append(FUS * a_end * a_end + (FVM - FUS) * a_end)
        # gas face flux = the identity G = vm - F: an independent Godunov
        # solve of the mirrored quadratic does not mirror the liquid
        # branch selection and breaks F + G = vm on non-uniform faces
        gvol = [VSG]
        for j in range(1, FN):
            gvol.append(FVM - flux[j])
        gvol.append(FVM + (FUS - FVM) * a_end - FUS * a_end * a_end)

        def rho_up(j, q):
            if q > 0:
                if j == 0:
                    return RHO_G
                if j == FN:
                    return mg[FN - 1] / ((1.0 - alpha[FN - 1]) * FDX)
                return mg[j - 1] / ((1.0 - alpha[j - 1]) * FDX)
            if q < 0:
                if j == FN:
                    return FOUT_P / (GAS_R * FOUT_T)
                return mg[j] / ((1.0 - alpha[j]) * FDX)
            return 0.0

        wflux = [0.0] * (FN + 1)
        hflux = [0.0] * (FN + 1)
        iflux = [0.0] * (FN + 1)
        for j in range(FN + 1):
            if flux[j] == 0:
                continue
            if j == 0:
                wflux[j] = FW0 * flux[j]
                iflux[j] = FINH_IN * flux[j]
            elif j == FN:
                wflux[j] = water[FN - 1] * flux[j]
                hflux[j] = hyd[FN - 1] * flux[j]
                iflux[j] = inh[FN - 1] * flux[j]
            else:
                wu = water[j - 1] if flux[j] >= 0 else water[j]
                hu = hyd[j - 1] if flux[j] >= 0 else hyd[j]
                iu = inh[j - 1] if flux[j] >= 0 else inh[j]
                wflux[j] = wu * flux[j]
                hflux[j] = hu * flux[j]
                iflux[j] = iu * flux[j]
        gmflux = [rho_up(j, gvol[j]) * gvol[j] for j in range(FN + 1)]
        for i in range(FN):
            a_new = alpha[i] + FDT / FDX * (flux[i] - flux[i + 1])
            a_cap = max(a_new, 0.0)
            af = alpha[i] * water[i] + FDT / FDX * (wflux[i] - wflux[i + 1])
            ah = alpha[i] * hyd[i] + FDT / FDX * (hflux[i] - hflux[i + 1])
            ai = alpha[i] * inh[i] + FDT / FDX * (iflux[i] - iflux[i + 1])
            alpha[i] = min(max(a_new, 0.0), 1.0)
            af = min(max(af, 0.0), a_cap)
            water[i] = af / a_cap if a_cap > 1e-9 else 0.0
            ah = min(max(ah, 0.0), a_cap)
            hyd[i] = ah / a_cap if a_cap > 1e-9 else 0.0
            ai = min(max(ai, 0.0), a_cap)
            inh[i] = ai / a_cap if a_cap > 1e-9 else 0.0
            mg[i] = max(mg[i] - FDT * (gmflux[i + 1] - gmflux[i]), 0.0)
        # ---- rigorous kinetics ----
        h_delta = [0.0] * FN
        for i in range(FN):
            ag = 1.0 - alpha[i]
            p_i = mg[i] * GAS_R * temp[i] / (ag * FDX) if ag > 1e-12 else 0.0
            dT = teq_rig(p_i, inh[i], FAK) - temp[i]
            total = water[i] + hyd[i]
            if dT >= 0:
                r = K_G * dT
                h_new = (hyd[i] / FDT + r * total) / (1.0 / FDT + r)
                trace = 1e-9 * RHO_G * FDX
                dh_max = (mg[i] - trace) / (GAS_CONTENT * alpha[i] * FDX)
                h_new = min(h_new, hyd[i] + max(dh_max, 0.0))
            else:
                h_new = hyd[i] / (1.0 + FDT * K_D * (-dT))
            h_new = min(max(h_new, 0.0), total)
            dh = h_new - hyd[i]
            water[i] = min(max(total - h_new, 0.0), 1.0)
            hyd[i] = h_new
            mg[i] -= GAS_CONTENT * alpha[i] * dh * FDX
            h_delta[i] = alpha[i] * dh
        # ---- implicit thermal on the step's phase fluxes ----
        mat = [[0.0] * FN for _ in range(FN)]
        rhs = [0.0] * FN
        for i in range(FN):
            cap = AREA * FDX * (alpha[i] * RHO_L * CL
                                + (1.0 - alpha[i])
                                * (mg[i] / ((1.0 - alpha[i]) * FDX)) * CG)
            kw = KW * FDX
            mat[i][i] += cap / FDT + kw
            rhs[i] += cap / FDT * temp[i] + kw * FWALL
            rhs[i] += LATENT * h_delta[i] * AREA * FDX / FDT
        for j in range(FN + 1):
            kl = AREA * RHO_L * CL * abs(flux[j])
            kg = AREA * rho_up(j, gvol[j]) * CG * abs(gvol[j])
            for q_f, k_f in ((flux[j], kl), (gvol[j], kg)):
                if q_f == 0:
                    continue
                up = j - 1 if q_f > 0 else j
                dn = j if q_f > 0 else j - 1
                if 0 <= dn < FN:
                    if 0 <= up < FN:
                        mat[dn][up] -= k_f
                    else:
                        rhs[dn] += k_f * (FIN if j == 0 else FOUT_T)
                    mat[dn][dn] += k_f
        temp = solve_dense(mat, rhs)
        steps_done += 1
    for i in range(FN):
        got = fr["cells"][i]
        ag_i = 1.0 - alpha[i]
        want_p = (mg[i] * GAS_R * temp[i] / (ag_i * FDX)
                  if ag_i > 1e-12 else 0.0)
        want_teq = teq_rig(want_p, inh[i], FAK)
        worst_te = max(worst_te, abs(got["Teq"] - want_teq))
        if (abs(got["alpha"] - alpha[i]) > 1e-9 + 1e-9 * alpha[i]
                or abs(got["water"] - water[i]) > 1e-9 + 1e-9 * water[i]
                or abs(got["inh"] - inh[i])
                > 1e-9 + 1e-9 * max(inh[i], 1e-3)
                or abs(got["h"] - hyd[i]) > 1e-9 + 1e-9 * hyd[i]
                or abs(got["T"] - temp[i]) > 1e-9 + 1e-9 * temp[i]
                or abs(got["p"] - want_p) > 1e-9 * max(want_p, 1.0)):
            badf = (k, fr["t"], i, got)
            break
    if badf:
        break
if badf:
    fails.append(f"flow replication mismatch frame {badf[0]} t={badf[1]} "
                 f"cell {badf[2]}: {badf[3]}")
else:
    for k in range(0, len(fframes), 5):
        ok(f"flow transport + rigorous kinetics + thermal replicated at "
           f"t={fframes[k]['t']} (1e-9)")
check("flow Teq column = rigorous formula from reported p/inh", worst_te,
      0.0, 0.0, 1e-6)

# front advection: the inh > 0.2 tip moves downstream, never back, and the
# front clears the line
front = []
for fr in fframes:
    xs = [i for i, c in enumerate(fr["cells"]) if c["inh"] > 0.2]
    front.append(xs[-1] if xs else -1)
mono = all(front[k] <= front[k + 1] for k in range(len(front) - 1))
if mono and front[0] == -1 and front[-1] == FN - 1:
    ok("inhibitor front tip advects monotonically downstream and clears "
       "the line")
else:
    fails.append(f"front not monotone/crossing: {front}")
check("inlet cell saturated at the injected fraction 0.4",
      fframes[-1]["cells"][0]["inh"], 0.4, 1e-6)
var_mid = max(
    (lambda v: sum((x - sum(v) / len(v)) ** 2 for x in v) / len(v))
    ([c["inh"] for c in fr["cells"]]) for fr in fframes[1:-1])
if var_mid > 1e-4:
    ok("front carries real heterogeneity (max variance of inh > 1e-4)")
else:
    fails.append("inhibitor field stayed uniform")
# temporal flip on a mid-line cell: growth while clean, decay behind
grow_f = decay_f = None
ci = 10
for k in range(1, len(fframes)):
    d = fframes[k]["cells"][ci]["h"] - fframes[k - 1]["cells"][ci]["h"]
    ds = fframes[k]["cells"][ci]["dTsub"]
    if d > 1e-6 and ds > 0.3 and grow_f is None:
        grow_f = fframes[k]["t"]
    if d < -1e-6 and ds < -0.3:
        decay_f = fframes[k]["t"]
        break
if grow_f is not None and decay_f is not None and decay_f > grow_f:
    ok(f"cell 10 grows at t={grow_f} (dTsub > 0.3), dissociates by "
       f"t={decay_f} (dTsub < -0.3) behind the front")
else:
    fails.append(f"no growth->decay flip in cell 10: {grow_f} {decay_f}")

# =====================================================================
# rejections
# =====================================================================
expect_fatal("inhibitorEquilibrium without inhibitor{}",
             open(CASE_REJECT).read())
no_hy = re.sub(r"hydrate\s*\{[^}]*\}\s*", "", CASE_TEXT)
expect_fatal("rigorous chain without hydrate{} (inhibitor{} rejects "
             "first)", no_hy)
expect_fatal("missing activityK",
             CASE_TEXT.replace("    activityK 8;\n", ""))
expect_fatal("activityK 0", CASE_TEXT.replace("activityK 8;", "activityK 0;"))
expect_fatal("activityK -1",
             CASE_TEXT.replace("activityK 8;", "activityK -1;"))
expect_fatal("unknown basis",
             CASE_TEXT.replace("basis liquid;", "basis salt;"))
expect_fatal("unknown enable word",
             CASE_TEXT.replace("activityK 8;", "enable yes; activityK 8;"))
expect_fatal("inhibitor initial out of [0,1]",
             CASE_TEXT.replace("initial 0.3;", "initial 1.5;"))
# a disabled block is inert: it runs even with the parents declared
run_text(CASE_TEXT.replace("activityK 8;", "enable false; activityK 8;"))
ok("enable false with declared parents runs (inert)")

print()
if fails:
    print(f"FAILED ({len(fails)}):")
    for f in fails:
        print(f"  - {f}")
    raise SystemExit(1)
print("all rigorous-inhibitor-equilibrium checks passed")
