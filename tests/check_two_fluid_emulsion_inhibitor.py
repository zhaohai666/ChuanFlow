#!/usr/bin/env python3
"""G8 emulsion mixing viscosity + G5 thermodynamic hydrate inhibition
for the distributed two-fluid kernel (chflow-1d-twofluid).

usage: check_two_fluid_emulsion_inhibitor.py <exe> <emulsion.case.dict>
         <inhibitor.case.dict>

Emulsion (G8): a rate-driven line with a uniform water cut prices the
Newtonian liquid share of the wall law at the declared mixing law on the
transported cut. The holdup kinematics are viscosity-blind, so the steady
holdup is the undeclared twin's alpha* and the law shows itself exactly in
the muEff profile and the Churchill pressure drop:
  brinkman:   mu_e = mu_oil (1 + 2.5 phi)
  richardson: mu_e = mu_oil^(1-phi) mu_water^phi
The dry line (phi = 0) lands on the emulsion-free run bit for bit, and a
power-law liquid alongside emulsion{} is refused.

Inhibitor (G5): the sealed subcooled line's kinetics see the
Hammerschmidt-offset equilibrium line dT = T_eq(p) - k c - T with the
transported fraction c (still line: c stays at its initial). The checker
replicates the offset recursion cell by cell (1e-9), witnesses full
suppression at k = 30 (h stays exactly 0, p pinned at p0) and the
undeclared twin trajectory at c = 0.
"""

import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 4:
    raise SystemExit(__doc__)
EXE, EMUL_CASE, INH_CASE = sys.argv[1:4]

fails = []


def check(name, got, want, rtol, atol=0.0):
    tol = atol + rtol * abs(want)
    ok = abs(got - want) <= tol
    print(f"{'ok:  ' if ok else 'FAIL:'} {name}: got {got!r} want {want!r}")
    if not ok:
        fails.append(name)


def run(path):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"case {path} failed rc={r.returncode}\n"
                         f"{r.stdout}\n{r.stderr}")
    frames = []
    for line in r.stdout.splitlines():
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


def expect_fatal(name, path):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    bad = r.returncode == 0
    msg = (r.stdout + r.stderr).lower()
    if not bad and "fatal" not in msg and "error" not in msg:
        bad = True
    print(f"{'ok:  ' if not bad else 'FAIL:'} {name} refused "
          f"(rc={r.returncode})")
    if bad:
        fails.append(name)


def churchill(re_, relrough):
    re_ = max(re_, 1e-12)
    a1 = (7.0 / re_) ** 0.9 + 0.27 * relrough
    a = (2.457 * math.log(1.0 / a1)) ** 16
    b = (37530.0 / re_) ** 16
    return 8.0 * ((8.0 / re_) ** 12 + (a + b) ** (-1.5)) ** (1.0 / 12.0)


def godunov(qa, qb, qc, aL, aR):
    # port of tests/check_two_fluid.py godunov (the kernel's exact face
    # flux for F(a) = qa a^2 + qb a + qc)
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


def read(path):
    return open(path).read()


# =====================================================================
# G8 emulsion
# =====================================================================
RHO_L, RHO_G = 1000.0, 20.0
MU_L, MU_G = 0.01, 1e-5
MU_W = 0.001
N, DI, DT = 20, 0.02, 0.1
DX = 1.0
VM, US, VSL, VSG, PHI = 0.4, 0.3, 0.15, 0.25, 0.4
ROUGH = 0.0
STEPS_PER_FRAME = 100

emul_text = read(EMUL_CASE)

# steady holdup root of F(a) = us a^2 + (vm - us) a = vsl
A_STAR = (-(VM - US) + math.sqrt((VM - US) ** 2 + 4.0 * US * VSL)) / (2.0 * US)


def emul_dp(cells, mu_of):
    dp = 0.0
    for c in cells:
        a = c["alpha"]
        rho = a * RHO_L + (1.0 - a) * RHO_G
        mu = a * mu_of(c["water"]) + (1.0 - a) * MU_G
        re_ = rho * VM * DI / mu
        dp += churchill(re_, ROUGH / DI) * (DX / DI) * 0.5 * rho * VM ** 2
    return dp


frames = run(EMUL_CASE)
MU_BK = MU_L * (1.0 + 2.5 * PHI)
alpha_e = [0.5] * N
steps_done = 0
for k, fr in enumerate(frames):
    want_steps = 0 if k == 0 else k * STEPS_PER_FRAME
    # replicate the holdup transport (Godunov, inlet rate imposed)
    while steps_done < want_steps:
        flux = [VSL]
        for j in range(1, N):
            flux.append(godunov(US, VM - US, 0.0, alpha_e[j - 1],
                                alpha_e[j]))
        a_end = alpha_e[N - 1]
        flux.append(US * a_end * a_end + (VM - US) * a_end)
        for i in range(N):
            na = alpha_e[i] + DT / DX * (flux[i] - flux[i + 1])
            alpha_e[i] = min(max(na, 0.0), 1.0)
        steps_done += 1
    bad = [i for i in range(N)
           if abs(fr["cells"][i]["alpha"] - alpha_e[i]) > 1e-9]
    if bad:
        fails.append(f"emulsion alpha cells {bad[:4]} (t={fr['t']})")
    mu_bad = [i for i in range(N)
              if abs(fr["cells"][i].get("muEff", -1.0) - MU_BK) > 1e-15]
    if mu_bad:
        fails.append(f"muEff != brinkman mu_oil (1+2.5 phi) cells "
                     f"{mu_bad[:4]} (t={fr['t']})")
    check(f"emulsion dp = Churchill sum with brinkman liquid share "
          f"(t={fr['t']})", fr["dp"], emul_dp(fr["cells"], lambda w: MU_BK),
          1e-9)
if not any("emulsion alpha" in f for f in fails):
    print("ok:   viscosity-blind holdup transport replicated (1e-9)")
if not any("muEff" in f for f in fails):
    print(f"ok:   muEff = mu_oil (1 + 2.5 * {PHI}) = {MU_BK} in every cell "
          f"(1e-15)")
last = frames[-1]["cells"]
bad = [i for i in range(N) if abs(last[i]["alpha"] - A_STAR) > 1e-9]
check("steady holdup lands on the undeclared alpha*", float(not bad), 1, 0)

# ---- richardson sibling ----
with tempfile.TemporaryDirectory(prefix="cf_emul_") as td:
    p = os.path.join(td, "rich.dict")
    open(p, "w").write(emul_text.replace("model brinkman",
                                         "model richardson"))
    fr_r = run(p)
    MU_R = MU_L ** (1.0 - PHI) * MU_W ** PHI
    bad = [i for i in range(N)
           if abs(fr_r[-1]["cells"][i]["muEff"] - MU_R) > 1e-11 * MU_R]
    check("richardson muEff = mu_oil^0.6 mu_water^0.4", float(not bad), 1, 0)
    check("richardson dp = Churchill sum with the product-rule share",
          fr_r[-1]["dp"],
          emul_dp(fr_r[-1]["cells"], lambda w: MU_R), 1e-9)

    # ---- dry-line identity: phi = 0 lands on the emulsion-free run ----
    dry = emul_text.replace("waterFraction 0.4", "waterFraction 0.0")
    dry = dry.replace("water 0.4;", "water 0.0;")
    pd = os.path.join(td, "dry.dict")
    open(pd, "w").write(dry)
    fr_d = run(pd)
    bad = [i for i in range(N)
           if abs(fr_d[-1]["cells"][i]["muEff"] - MU_L) > 1e-18]
    check("dry line muEff == mu_oil exactly", float(not bad), 1, 0)
    off = re.sub(r"emulsion\s*\{[^}]*\}", "", dry)
    po = os.path.join(td, "off.dict")
    open(po, "w").write(off)
    fr_o = run(po)
    # the dry run still prints the muEff column (emulsion is declared);
    # the physics identity is every OTHER field bit-for-bit
    strip = [{k: v for k, v in c.items() if k != "muEff"}
             for c in fr_d[-1]["cells"]]
    same = (strip == fr_o[-1]["cells"]
            and fr_d[-1]["dp"] == fr_o[-1]["dp"]
            and fr_d[-1]["vm"] == fr_o[-1]["vm"])
    print(f"{'ok:  ' if same else 'FAIL:'} dry emulsion run bit-identical "
          f"to the emulsion-free run")
    if not same:
        fails.append("dry emulsion not bit-identical")

    # ---- rejections ----
    pr = os.path.join(td, "rej_flowindex.dict")
    open(pr, "w").write(emul_text.replace("muG  1e-5;",
                                          "muG  1e-5;\n    flowIndex 0.6;"))
    expect_fatal("emulsion{} with a power-law liquid", pr)
    pr = os.path.join(td, "rej_model.dict")
    open(pr, "w").write(emul_text.replace("model brinkman", "model pal-rhodes"))
    expect_fatal("unknown emulsion model word", pr)
    pr = os.path.join(td, "rej_muwater.dict")
    open(pr, "w").write(emul_text.replace("muWater 0.001", "muWater 0"))
    expect_fatal("non-positive emulsion muWater", pr)


# =====================================================================
# G5 inhibitor
# =====================================================================
N_I, DX_I, DI_I = 10, 2.0, 0.02
AREA_I = math.pi * DI_I * DI_I / 4.0
A0, W0, T0 = 0.5, 0.5, 300.0
RHO_G_I = 200.0
GAS_R, KW, T_WALL = 287.0, 50.0, 280.0
CL, CG = 4180.0, 2150.0
CURVE_A, CURVE_B = 35.2, -5684.0
K_G, K_D = 0.01, 0.01
GAS_CONTENT, LATENT = 122.0, 4.2e8
DT_I, STEPS_PER_FRAME_I = 0.5, 160

inh_text = read(INH_CASE)


def replicate_sealed(k_offset):
    """the sealed per-cell recursion with the equilibrium line depressed by
    k_offset (a constant offset for the whole run on a still line)"""
    alpha = [A0] * N_I
    water = [W0] * N_I
    hyd = [0.0] * N_I
    temp = [T0] * N_I
    mg = [RHO_G_I * (1.0 - A0) * DX_I] * N_I
    frames_i = []
    for step in range(0, 1601):
        if step % STEPS_PER_FRAME_I == 0:
            frames_i.append({"t": step * DT_I,
                             "w": list(water), "h": list(hyd),
                             "T": list(temp),
                             "p": [mg[i] * GAS_R * temp[i]
                                   / ((1.0 - alpha[i]) * DX_I)
                                   for i in range(N_I)]})
        if step == 1600:
            break
        for i in range(N_I):
            ag = 1.0 - alpha[i]
            p_i = mg[i] * GAS_R * temp[i] / (ag * DX_I)
            dT = CURVE_B / (math.log(p_i) - CURVE_A) - k_offset - temp[i]
            total = water[i] + hyd[i]
            if dT >= 0:
                r = K_G * dT
                h_new = (hyd[i] / DT_I + r * total) / (1.0 / DT_I + r)
                trace = 1e-9 * RHO_G_I * DX_I
                dh_max = (mg[i] - trace) / (GAS_CONTENT * alpha[i] * DX_I)
                h_new = min(h_new, hyd[i] + max(dh_max, 0.0))
            else:
                h_new = hyd[i] / (1.0 + DT_I * K_D * (-dT))
            h_new = min(max(h_new, 0.0), total)
            dh = h_new - hyd[i]
            water[i] = min(max(total - h_new, 0.0), 1.0)
            hyd[i] = h_new
            mg[i] -= GAS_CONTENT * alpha[i] * dh * DX_I
            cap = AREA_I * DX_I * (alpha[i] * RHO_L * CL
                                   + ag * (mg[i] / (ag * DX_I)) * CG)
            kw = KW * DX_I
            temp[i] = (cap / DT_I * temp[i] + kw * T_WALL
                       + LATENT * alpha[i] * dh * AREA_I * DX_I / DT_I) \
                      / (cap / DT_I + kw)
    return frames_i


frames_i = run(INH_CASE)
rep = replicate_sealed(5.0)
ok_all = True
for k, fr in enumerate(frames_i):
    r = rep[k]
    for i in range(N_I):
        got = fr["cells"][i]
        if (abs(got["water"] - r["w"][i]) > 1e-9 + 1e-9 * r["w"][i]
                or abs(got.get("h", 0.0) - r["h"][i]) > 1e-9
                or abs(got.get("T", 0.0) - r["T"][i]) > 1e-9
                    + 1e-9 * r["T"][i]
                or abs(got.get("p", 0.0) - r["p"][i]) > 1e-9 * r["p"][i]):
            fails.append(f"inhibitor cell {i} (t={fr['t']}): got={got} "
                         f"want w={r['w'][i]} h={r['h'][i]}")
            ok_all = False
            break
    if not ok_all:
        break
    if abs(got.get("inh", -1.0) - 1.0) > 1e-15:
        fails.append(f"inhibitor fraction drifted on a still line "
                     f"(t={fr['t']})")
        ok_all = False
        break
if ok_all:
    print("ok:   offset kinetics recursion dT = T_eq(p) - 5 c - T "
          "replicated (1e-9), c pinned at 1")

with tempfile.TemporaryDirectory(prefix="cf_inh_") as td:
    # c = 0 sibling: the undeclared twin's trajectory verbatim
    p0 = os.path.join(td, "c0.dict")
    open(p0, "w").write(inh_text.replace("initial 1;", "initial 0;"))
    fr0 = run(p0)
    rep0 = replicate_sealed(0.0)
    ok = True
    for k, fr in enumerate(fr0):
        r = rep0[k]
        for i in range(N_I):
            got = fr["cells"][i]
            if (abs(got.get("h", 0.0) - r["h"][i]) > 1e-9
                    or abs(got.get("T", 0.0) - r["T"][i]) > 1e-9
                        + 1e-9 * r["T"][i]):
                ok = False
                break
        if not ok:
            break
    print(f"{'ok:  ' if ok else 'FAIL:'} c = 0 lands on the undeclared "
          f"twin's trajectory")
    if not ok:
        fails.append("c=0 twin mismatch")

    # k = 30 sibling: the drive is negative for the whole run -> no growth,
    # p pinned at p0, T a pure wall relaxation
    p30 = os.path.join(td, "k30.dict")
    open(p30, "w").write(inh_text.replace("kDepression 5;",
                                          "kDepression 30;"))
    fr30 = run(p30)
    h_bad = [i for i in range(N_I)
             if fr30[-1]["cells"][i].get("h", 0.0) != 0.0]
    check("suppressed line keeps h exactly 0", float(not h_bad), 1, 0)
    w_bad = [i for i in range(N_I)
             if abs(fr30[-1]["cells"][i]["water"] - W0) > 1e-15]
    check("suppressed line keeps the free water untouched",
          float(not w_bad), 1, 0)
    # pure wall relaxation, closed form: T_n = Tw + (T0 - Tw) g^n with
    # g = (cap/dt)/(cap/dt + kw), cap constant (no gas consumed)
    cap0 = AREA_I * DX_I * (A0 * RHO_L * CL
                            + (1.0 - A0) * RHO_G_I * CG)
    g = (cap0 / DT_I) / (cap0 / DT_I + KW * DX_I)
    t_end = T_WALL + (T0 - T_WALL) * g ** 1600
    check("suppressed temperature = pure wall relaxation",
          fr30[-1]["cells"][0].get("T", 0.0), t_end, 1e-9)
    # the gas inventory is untouched, so the EOS pressure tracks the
    # reported temperature exactly: p = rhoG0 gasR T (no latent spikes)
    p_bad = [i for i in range(N_I)
             if abs(fr30[-1]["cells"][i].get("p", 0.0)
                    - RHO_G_I * GAS_R * fr30[-1]["cells"][i]["T"])
             > 1e-9 * RHO_G_I * GAS_R * fr30[-1]["cells"][i]["T"]]
    check("suppressed line's p = rhoG0 gasR T (inventory untouched)",
          float(not p_bad), 1, 0)

    # ---- rejections ----
    ph = os.path.join(td, "nohydrate.dict")
    open(ph, "w").write(re.sub(r"hydrate\s*\{[^}]*\}", "", inh_text))
    expect_fatal("inhibitor{} without hydrate{}", ph)
    pk = os.path.join(td, "k0.dict")
    open(pk, "w").write(inh_text.replace("kDepression 5;",
                                         "kDepression 0;"))
    expect_fatal("non-positive kDepression", pk)
    pi = os.path.join(td, "ini.dict")
    open(pi, "w").write(inh_text.replace("initial 1;", "initial 1.5;"))
    expect_fatal("inhibitor initial outside [0,1]", pi)

print()
if fails:
    print(f"FAILED ({len(fails)}):")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("all emulsion/inhibitor judgments passed")
