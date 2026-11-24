#!/usr/bin/env python3
"""G10 1D sand / cuttings transport for the distributed two-fluid kernel
(chflow-1d-twofluid).

usage: check_two_fluid_sand_transport.py <exe> <case.dict>

Criterion makeup: 3 independent physical criteria + mirror diagnostics
(regression safety net). The mirror no longer carries the verdict alone.

Independent physical criteria:

  I1. Inventory conservation closure. On a dense twin (reportEvery 1) the
      total inventory sum_c (alpha*sand_c + bed_c) must change, every single
      step, by exactly the upwind through-flow: dI = dt/dx * q_L * (s_in -
      s_last) -- settling/erosion only exchanges between suspension and bed,
      so the balance is closed by the boundary fluxes alone. Checked on the
      settling frame (s_in = 0.2) and the erosion twin (s_in = 0).

  I2. Thin-layer monotonicity. With feedback on, dEff = D*sqrt(1 - B): the
      bore thins monotonically wherever the bed grows (settling run) and
      thickens monotonically wherever the bed decays (erosion run).

  I3. Critical-velocity band. Two siblings straddling uCrit -- vM = 0.99
      uCrit and vM = 1.01 uCrit, same sand everywhere -- must land on
      opposite branches: the bed grows below uCrit and decays above it, by
      margins the erode/settle rates make unambiguous.

Mirror diagnostics (regression safety net):

  M1. Settling branch (main frame, vM = 0.6 < uCrit = 1.0): the suspended
      fraction rides the upwind transport from the inlet composition and the
      bed grows by exactly the removed liquid share, so the checker
      replicates the kernel's per-step recursion -- transport (explicit Euler
      on the constant liquid flux) then the implicit-Euler settling
      S -> S / (1 + dt vt/(D/2)) with B += alpha (S - S_new). Cell 0
      therefore relaxes to the settling/replenishment steady state S = 0.15
      (0.6 (0.2 - S) = 0.2 S), not a pure decay. The thin-layer law
      dEff = D sqrt(1 - B) prices the bore.
  M2. Erosion sibling (vM = 1.2 >= uCrit): the bed decays
      B_n = B_0 (1 + dt erodeRate (vM/uCrit - 1))^-n while the eroded sand
      rides the upwind transport.
  M3. Rejections are loud.
"""

import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)
EXE, CASE = sys.argv[1:3]

fails = []


def check(name, got, want, rtol, atol=0.0):
    tol = atol + rtol * abs(want)
    ok = abs(got - want) <= tol
    print(f"{'ok:  ' if ok else 'FAIL:'} {name}: got {got!r} want {want!r}")
    if not ok:
        fails.append(name)


def expect_fatal(name, path):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    msg = (r.stdout + r.stderr).lower()
    bad = r.returncode == 0 or ("fatal" not in msg and "error" not in msg)
    print(f"{'ok:  ' if not bad else 'FAIL:'} {name} refused "
          f"(rc={r.returncode})")
    if bad:
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
                           "cells": []})
        elif tok[0] == "cell" and frames:
            c, key = {}, None
            for w in tok[2:]:
                if key is None:
                    key = w
                else:
                    c[key] = float(w)
                    key = None
            frames[-1]["cells"].append(c)
    return frames


N, DX, DI = 10, 1.0, 0.1
DT = 0.5
SPF = 20
A0, S0 = 0.5, 0.2
VSL, VSG = 0.3, 0.3
VT, UCRIT, ERODE = 0.01, 1.0, 0.05
K_SETTLE = VT / (0.5 * DI)

text = open(CASE).read()
frames = run(CASE)

# ---- settling recursion + thin-layer bore ----
S_IN = 0.2  # inlet sandFraction (phaseRates injection)
FLUX = VSL  # uniform alpha keeps the liquid flux constant along the line
s = [S0] * N
b = [0.0] * N
ok = True
step = 0
for k, fr in enumerate(frames):
    want_steps = k * SPF
    while step < want_steps:
        # transport first (upwind on the constant liquid flux, inlet
        # composition S_IN), then the settling exchange, kernel order
        ns = list(s)
        for i in range(N):
            fin = FLUX * (S_IN if i == 0 else s[i - 1])
            fout = FLUX * s[i]
            content = A0 * s[i] + DT / DX * (fin - fout)
            ns[i] = max(content, 0.0) / A0
        kb = 1.0 + DT * K_SETTLE
        for i in range(N):
            snew = ns[i] / kb
            b[i] += A0 * (ns[i] - snew)
            ns[i] = snew
        s = ns
        step += 1
    for i, c in enumerate(fr["cells"]):
        if (abs(c["sand"] - s[i]) > 1e-9 * max(s[i], 1e-12) + 1e-12
                or abs(c["bed"] - b[i]) > 1e-12):
            fails.append(f"settling: t={fr['t']} cell {i}: got "
                         f"({c['sand']}, {c['bed']}) want "
                         f"({s[i]}, {b[i]})")
            ok = False
            break
        want_d = DI * math.sqrt(max(1.0 - b[i], 0.0))
        if abs(c["dEff"] - want_d) > 1e-9 * want_d:
            fails.append(f"thin-layer bore: t={fr['t']} cell {i}: "
                         f"{c['dEff']} != {want_d}")
            ok = False
            break
    if not ok:
        break
if ok:
    print(f"ok:   settling transport+decay recursion replicated, bed grows by "
          f"the removed content, dEff = D sqrt(1 - B) on all {len(frames)} "
          f"frames (cell 0 -> steady 0.15)")

# ---- erosion sibling: vM = 1.2 > uCrit, bedInitial 0.1 ----
with tempfile.TemporaryDirectory(prefix="cf_sand_") as td:
    ero = text.replace("vsl 0.3; vsg 0.3", "vsl 0.6; vsg 0.6")
    ero = ero.replace("initial 0.2;", "initial 0;")
    ero = ero.replace("bedInitial 0;", "bedInitial 0.1;")
    ero = ero.replace("sandFraction 0.2", "sandFraction 0.0")
    p = os.path.join(td, "ero.dict")
    open(p, "w").write(ero)
    fr_e = run(p)
    VM_E = 1.2
    K_ER = ERODE * (VM_E / UCRIT - 1.0)
    FLUX = 0.6  # liquid flux (vsl imposed; uniform alpha keeps it constant)
    s = [0.0] * N
    b = [0.1] * N
    ok = True
    step = 0
    for k, fr in enumerate(fr_e):
        want_steps = k * SPF
        while step < want_steps:
            # transport first (upwind on the constant liquid flux, inlet
            # composition 0), then the erosion exchange
            ns = list(s)
            for i in range(N):
                fin = FLUX * (0.0 if i == 0 else s[i - 1])
                fout = FLUX * s[i]
                content = A0 * s[i] + DT / DX * (fin - fout)
                ns[i] = max(content, 0.0) / A0
            kb = 1.0 + DT * K_ER
            for i in range(N):
                bnew = b[i] / kb
                ns[i] = min(ns[i] + (b[i] - bnew) / A0, 1.0)
                b[i] = bnew
            s = ns
            step += 1
        for i, c in enumerate(fr["cells"]):
            if (abs(c["sand"] - s[i]) > 1e-9 * max(s[i], 1e-12)
                    + 1e-12 or abs(c["bed"] - b[i]) > 1e-12):
                fails.append(f"erosion: t={fr['t']} cell {i}: got "
                             f"({c['sand']}, {c['bed']}) want "
                             f"({s[i]}, {b[i]})")
                ok = False
                break
        if not ok:
            break
    if ok:
        print("ok:   erosion bed decay B_n = 0.1 (1+dt erodeRate (vM/uCrit-1))"
              "^-n advected by the upwind transport replicated (1e-9)")
    if b[0] > 0.09:
        fails.append("erosion barely moved the bed")

# ---- I1: inventory conservation closure on dense twins (every step) ----
DENSE = text.replace("reportEvery 20;", "reportEvery 1;")
FLUX_S = VSL  # settling twin keeps the original vsl 0.3 line
with tempfile.TemporaryDirectory(prefix="cf_sand_dense_") as td:
    # settling twin: inlet injects S_IN = 0.2 with the liquid flux
    p = os.path.join(td, "dense_settle.dict")
    open(p, "w").write(DENSE)
    fr_d = run(p)
    ok = True
    for k in range(len(fr_d) - 1):
        inv0 = sum(A0 * c["sand"] + c["bed"] for c in fr_d[k]["cells"])
        inv1 = sum(A0 * c["sand"] + c["bed"] for c in fr_d[k + 1]["cells"])
        s_last = fr_d[k]["cells"][-1]["sand"]
        want = DT / DX * FLUX_S * (S_IN - s_last)
        got = inv1 - inv0
        if abs(got - want) > 1e-9 * max(abs(want), 1e-3) + 1e-12:
            fails.append(f"I1 settling closure: t={fr_d[k + 1]['t']}: "
                         f"dInventory {got!r} != through-flow {want!r}")
            ok = False
            break
    if ok:
        print(f"ok:   I1 settling inventory closure dI = dt/dx qL (sIn - sLast) "
              f"on all {len(fr_d) - 1} dense steps")

    # erosion twin: inlet sand 0, initial 0, bedInitial 0.1 (mirror recipe)
    ero_d = DENSE.replace("vsl 0.3; vsg 0.3", "vsl 0.6; vsg 0.6")
    ero_d = ero_d.replace("initial 0.2;", "initial 0;")
    ero_d = ero_d.replace("bedInitial 0;", "bedInitial 0.1;")
    ero_d = ero_d.replace("sandFraction 0.2", "sandFraction 0.0")
    p = os.path.join(td, "dense_erode.dict")
    open(p, "w").write(ero_d)
    fr_de = run(p)
    FLUX_E = 0.6
    ok = True
    for k in range(len(fr_de) - 1):
        inv0 = sum(A0 * c["sand"] + c["bed"] for c in fr_de[k]["cells"])
        inv1 = sum(A0 * c["sand"] + c["bed"] for c in fr_de[k + 1]["cells"])
        s_last = fr_de[k]["cells"][-1]["sand"]
        want = -DT / DX * FLUX_E * s_last
        got = inv1 - inv0
        if abs(got - want) > 1e-9 * max(abs(want), 1e-3) + 1e-12:
            fails.append(f"I1 erosion closure: t={fr_de[k + 1]['t']}: "
                         f"dInventory {got!r} != outflow {-want!r}")
            ok = False
            break
    if ok:
        print(f"ok:   I1 erosion inventory closure dI = -dt/dx qL sLast "
              f"on all {len(fr_de) - 1} dense steps")

# ---- I2: thin-layer dEff monotonicity wherever the bed moves ----
bad = []
for k in range(len(frames) - 1):
    for i, (c0, c1) in enumerate(zip(frames[k]["cells"], frames[k + 1]["cells"])):
        if c1["dEff"] > c0["dEff"] + 1e-12:
            bad.append(f"settling t={frames[k + 1]['t']} cell {i}: "
                       f"dEff rose {c0['dEff']} -> {c1['dEff']}")
            break
    if bad:
        break
if bad:
    fails.extend(bad)
else:
    print("ok:   I2 settling dEff = D sqrt(1-B) monotone thinning on all frames")
bad = []
for k in range(len(fr_e) - 1):
    for i, (c0, c1) in enumerate(zip(fr_e[k]["cells"], fr_e[k + 1]["cells"])):
        if c1["dEff"] < c0["dEff"] - 1e-12:
            bad.append(f"erosion t={fr_e[k + 1]['t']} cell {i}: "
                       f"dEff fell {c0['dEff']} -> {c1['dEff']}")
            break
    if bad:
        break
if bad:
    fails.extend(bad)
else:
    print("ok:   I2 erosion dEff monotone thickening on all frames")

# ---- I3: critical-velocity band -- 0.99 uCrit settles, 1.01 uCrit erodes ----
with tempfile.TemporaryDirectory(prefix="cf_sand_band_") as td:
    below = text.replace("vsl 0.3; vsg 0.3", "vsl 0.495; vsg 0.495")
    above = text.replace("vsl 0.3; vsg 0.3", "vsl 0.505; vsg 0.505")
    above = above.replace("bedInitial 0;", "bedInitial 0.1;")
    p = os.path.join(td, "below.dict")
    open(p, "w").write(below)
    fr_b = run(p)
    beds = [c["bed"] for c in fr_b[-1]["cells"]]
    if not all(b > 1e-3 for b in beds):
        fails.append(f"I3: vM = 0.99 uCrit did not grow the bed everywhere "
                     f"(final beds {beds})")
    else:
        print(f"ok:   I3 vM = 0.99 uCrit: bed grows everywhere "
              f"(final beds {['%.3f' % b for b in beds]})")
    p = os.path.join(td, "above.dict")
    open(p, "w").write(above)
    fr_a = run(p)
    beds = [c["bed"] for c in fr_a[-1]["cells"]]
    if not all(b < 0.1 - 1e-3 for b in beds):
        fails.append(f"I3: vM = 1.01 uCrit did not erode the bed everywhere "
                     f"(final beds {beds})")
    else:
        print(f"ok:   I3 vM = 1.01 uCrit: bed erodes everywhere "
              f"(final beds {['%.3f' % b for b in beds]})")

    # ---- rejections ----
    p = os.path.join(td, "r1.dict")
    open(p, "w").write(text.replace("terminalVelocity 0.01",
                                    "terminalVelocity 0"))
    expect_fatal("non-positive terminalVelocity", p)
    p = os.path.join(td, "r2.dict")
    open(p, "w").write(text.replace("criticalVelocity 1.0",
                                    "criticalVelocity 0"))
    expect_fatal("non-positive criticalVelocity", p)
    p = os.path.join(td, "r3.dict")
    open(p, "w").write(text.replace("initial 0.2;", "initial 1.5;"))
    expect_fatal("suspended initial outside [0,1]", p)
    p = os.path.join(td, "r4.dict")
    open(p, "w").write(text.replace("bedInitial 0;", "bedInitial 1;"))
    expect_fatal("bed initial outside [0,1)", p)
    p = os.path.join(td, "r5.dict")
    open(p, "w").write(text.replace("feedback true", "feedback maybe"))
    expect_fatal("sand feedback word outside true/false", p)

print()
if fails:
    print(f"FAILED ({len(fails)}):")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("all sand-transport judgments passed")
