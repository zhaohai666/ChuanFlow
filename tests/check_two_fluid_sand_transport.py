#!/usr/bin/env python3
"""G10 1D sand / cuttings transport for the distributed two-fluid kernel
(chflow-1d-twofluid).

usage: check_two_fluid_sand_transport.py <exe> <case.dict>

Settling branch (main frame, vM = 0.6 < uCrit = 1.0): the suspended
fraction rides the upwind transport from the inlet composition and the
bed grows by exactly the removed liquid share, so the checker replicates
the kernel's per-step recursion -- transport (explicit Euler on the
constant liquid flux) then the implicit-Euler settling S -> S / (1 + dt
vt/(D/2)) with B += alpha (S - S_new). Cell 0 therefore relaxes to the
settling/replenishment steady state S = 0.15 (0.6 (0.2 - S) = 0.2 S), not
a pure decay. The thin-layer law dEff = D sqrt(1 - B) prices the bore.
Erosion sibling (vM = 1.2 >= uCrit): the bed decays
B_n = B_0 (1 + dt erodeRate (vM/uCrit - 1))^-n while the eroded sand
rides the upwind transport. Rejections are loud.
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
