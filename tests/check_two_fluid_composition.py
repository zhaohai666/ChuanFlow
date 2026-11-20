#!/usr/bin/env python3
"""G1 transient composition transport + per-cell PR flash diagnostics for
the distributed two-fluid kernel (chflow-1d-twofluid).

usage: check_two_fluid_composition.py <exe> <case.dict>

A horizontal pressure-driven line (fixedPressure 12 MPa inlet, 2 MPa
outlet, wallis slip 0 -> uniform stationary holdup, liquid flux = alpha *
vM) whose initial gas inventory sits at the outlet pressure: the dense
inlet ghost injects rhoG = p_in/(gasR T), so a compression front and the
liquid composition front (zInlet heavy feed vs the declared light feed
z0) travel together at vM. The checker

  * replicates the kernel's per-step upwind recursion for the N - 1
    transported mole fractions (content alpha*z on the flux alpha*vM read
    from each frame's own vm, inlet injection zInlet, last fraction closed
    by 1 - sum) and compares every reported z column;
  * applies the Wilson-K single-phase clamps independently on every cell
    of every frame (K_i = pc_i/p * exp(5.373 (1+w_i)(1 - Tc_i/T)) with the
    case's component table hardcoded here): f(0) = sum z(K-1) < 0 forces
    beta == 0, f(1) = sum z(K-1)/K > 0 forces beta == 1;
  * witnesses both regimes in the run (a high-pressure all-liquid inlet
    section with beta == 0 and a low-pressure flashing section with
    0 < beta < 1);
  * proves beta monotone (non-increasing) in the reported cell pressure
    along the known monotone p profile;
  * proves the dry twin: removing the composition block leaves the whole
    report bit-identical (modulo the comp/z columns and the declaration
    line the wet run adds);
  * and that the feed validation is loud (negative z, sum != 1, single
    component, > 8 components, unknown name, wrong zInlet size, no flash
    pressure).
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
    return r.stdout


def parse(stdout):
    frames = []
    for line in stdout.splitlines():
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


# ---- case constants (mirrored from case.dict) ----
N, DX = 10, 10.0
DT = 0.1
A0 = 0.5
T0 = 300.0
NC = 3
Z0_RAW = [0.5, 0.3, 0.2]
ZIN_RAW = [0.1, 0.3, 0.6]
# Wilson inputs hardcoded from cf::pvt component data for the declared
# feed: (Tc [K], pc [Pa], omega) -- methane, propane, nHexane
CRIT = [("methane", 190.564, 45.992e5, 0.0115),
        ("propane", 369.834, 42.477e5, 0.1523),
        ("nHexane", 507.600, 30.250e5, 0.3013)]


def normalise(z):
    s = sum(z)
    return [v / s for v in z]


Z0, ZIN = normalise(Z0_RAW), normalise(ZIN_RAW)


def wilsonK(p):
    return [pc / p * math.exp(5.373 * (1.0 + w) * (1.0 - tc / T0))
            for (_n, tc, pc, w) in CRIT]


text = open(CASE).read()
stdout = run(CASE)
frames = parse(stdout)

# the declaration line names the feed and its z0
decl = [ln for ln in stdout.splitlines() if "twoFluid composition:" in ln]
if (len(decl) != 1
        or not all(name in decl[0] for name, *_ in CRIT)):
    fails.append(f"composition declaration line missing/incomplete: {decl}")
else:
    print("ok:   declaration line lists the 3 components and their z0")

# ---- z transport recursion (kernel order: transport, per component, on
# the frame's own flux alpha * vm; content clamped to [0, alpha]) ----
z = [list(Z0) for _ in range(N)]
step = 0
ok = True
seen_liquid = seen_two_phase = False
for k, fr in enumerate(frames):
    while step < k:  # reportEvery 1: frame k is produced by step k
        vm = frames[step + 1]["vm"]  # the vm printed with the step's frame
        flux = A0 * vm
        ns = [list(c) for c in z]
        for i in range(N):
            for kk in range(NC - 1):
                fin = flux * (ZIN[kk] if i == 0 else z[i - 1][kk])
                fout = flux * z[i][kk]
                content = A0 * z[i][kk] + DT / DX * (fin - fout)
                content = min(max(content, 0.0), A0)
                ns[i][kk] = content / A0
        for i in range(N):
            ns[i][NC - 1] = 1.0 - ns[i][0] - ns[i][1]
        z = ns
        step += 1
    for i, c in enumerate(fr["cells"]):
        if abs(c["alpha"] - A0) > 1e-12:
            fails.append(f"holdup moved: t={fr['t']} cell {i} "
                         f"alpha {c['alpha']}")
            ok = False
            break
        for kk in range(NC):
            if abs(c[f"z{kk}"] - z[i][kk]) > 1e-9 * max(abs(z[i][kk]), 1e-12) \
                    + 1e-12:
                fails.append(f"z recursion: t={fr['t']} cell {i} comp {kk}: "
                             f"got {c[f'z{kk}']} want {z[i][kk]}")
                ok = False
                break
        if not ok:
            break
    if not ok:
        break
if ok:
    print(f"ok:   upwind z recursion (N-1 transported + closure) replicated "
          f"on all {len(frames)} frames, holdup stationary")

# ---- Wilson-K clamps + monotonicity on every cell of every frame ----
ok = True
n_clamp_liq = n_two = 0
for fr in frames:
    cells = fr["cells"]
    for i, c in enumerate(cells):
        zz = [c[f"z{kk}"] for kk in range(NC)]
        K = wilsonK(c["p"])
        f0 = sum(zz[j] * (K[j] - 1.0) for j in range(NC))
        f1 = sum(zz[j] * (K[j] - 1.0) / K[j] for j in range(NC))
        beta = c["comp"]
        if f0 < -1e-9 and beta > 1e-12:
            fails.append(f"Wilson f(0)<0 but beta>0: t={fr['t']} cell {i} "
                         f"f0={f0} beta={beta}")
            ok = False
        if f1 > 1e-9 and beta < 1.0 - 1e-12:
            fails.append(f"Wilson f(1)>0 but beta<1: t={fr['t']} cell {i} "
                         f"f1={f1} beta={beta}")
            ok = False
        if beta > 1e-12 and beta < 1.0 - 1e-12:
            n_two += 1
            if not 0.0 < beta < 1.0:
                fails.append(f"beta out of (0,1): {beta}")
                ok = False
        if beta <= 1e-12 and f0 < -1e-9:
            n_clamp_liq += 1
        if not ok:
            break
    # p profile monotone non-increasing along the line, beta monotone
    # non-increasing in p (non-decreasing along the line)
    for i in range(N - 1):
        if cells[i + 1]["p"] > cells[i]["p"] + 1e-6:
            fails.append(f"p profile rose downstream: t={fr['t']} cell {i}")
            ok = False
        if cells[i + 1]["comp"] < cells[i]["comp"] - 1e-12:
            fails.append(f"beta not monotone in p: t={fr['t']} cell {i}: "
                         f"p {cells[i]['p']}->{cells[i+1]['p']} but beta "
                         f"{cells[i]['comp']}->{cells[i+1]['comp']}")
            ok = False
    if not ok:
        break
if ok:
    print(f"ok:   Wilson-K clamps hold on every cell (all-liquid witnesses "
          f"{n_clamp_liq}, two-phase witnesses {n_two}), beta monotone in p")

# ---- both regimes present: high-pressure all-liquid inlet section AND a
# low-pressure flashing section ----
liq = [(fr["t"], i) for fr in frames for i, c in enumerate(fr["cells"])
       if c["comp"] <= 1e-12 and c["p"] > 8e6]
flash = [(fr["t"], i) for fr in frames for i, c in enumerate(fr["cells"])
         if 1e-6 < c["comp"] < 1.0 - 1e-6 and c["p"] < 3e6]
same_frame = any(
    any(c["comp"] <= 1e-12 and c["p"] > 8e6 for c in fr["cells"])
    and any(1e-6 < c["comp"] < 1.0 - 1e-6 for c in fr["cells"])
    for fr in frames)
if liq and flash and same_frame:
    print(f"ok:   both regimes witnessed ({len(liq)} high-p all-liquid "
          f"cells, {len(flash)} low-p flashing cells, coexisting in one "
          f"frame)")
else:
    fails.append(f"regimes missing: liq={len(liq)} flash={len(flash)} "
                 f"same_frame={same_frame}")

# ---- dry twin: without the composition block the report is bit-identical
# (wet output minus the comp/z columns and the declaration line) ----
dry_text = re.sub(r"composition\s*\{[^}]*\}", "", text)
with tempfile.TemporaryDirectory(prefix="cf_comp_") as td:
    p = os.path.join(td, "dry.dict")
    open(p, "w").write(dry_text)
    dry = run(p)
    wet_norm = []
    for ln in stdout.splitlines():
        if "twoFluid composition:" in ln:
            continue
        wet_norm.append(re.sub(r" comp \S+(?: z\d \S+)+", "", ln))
    if wet_norm == dry.splitlines():
        print("ok:   dry twin (composition block removed) bit-identical")
    else:
        fails.append("dry twin diverged from the wet run")
        a, b = wet_norm, dry.splitlines()
        for j in range(max(len(a), len(b))):
            x = a[j] if j < len(a) else "<eof>"
            y = b[j] if j < len(b) else "<eof>"
            if x != y:
                print(f"first diff line {j}:\n  wet: {x}\n  dry: {y}")
                break

    # ---- loud rejections ----
    def write(name, t):
        q = os.path.join(td, name)
        open(q, "w").write(t)
        return q

    expect_fatal("negative mole fraction", write(
        "r1.dict", re.sub(r"methane 0\.5;", "methane -0.1;", text)))
    expect_fatal("sum != 1", write(
        "r2.dict", re.sub(r"nHexane 0\.2;", "nHexane 0.1;", text)))
    expect_fatal("single component", write(
        "r3.dict", re.sub(r"composition\s*\{[^}]*\}",
                          "composition { methane 1; }", text)))
    expect_fatal("more than 8 components", write(
        "r4.dict", re.sub(
            r"composition\s*\{[^}]*\}",
            "composition { methane 0.2; ethane 0.1; propane 0.1; iButane"
            " 0.1; nButane 0.1; nPentane 0.1; nHexane 0.1; carbonDioxide"
            " 0.1; nitrogen 0.1; }", text)))
    expect_fatal("unknown component name", write(
        "r5.dict", re.sub(r"methane 0\.5;", "methanol 0.5;", text)))
    expect_fatal("zInlet wrong size", write(
        "r6.dict", text.replace("zInlet (0.1 0.3 0.6)", "zInlet (0.1 0.9)")))
    norel = text.replace("    gasR 500;\n", "")
    norel = re.sub(r"inlet\s*\{[^}]*\}",
                   "inlet { type phaseRates; vsl 0.5; vsg 0.5; }", norel)
    norel = re.sub(r"outlet\s*\{[^}]*\}", "outlet { type sealed; }", norel)
    expect_fatal("no flash pressure (no EOS, no fixedPressure)", write(
        "r7.dict", norel))

print()
if fails:
    print(f"FAILED ({len(fails)}):")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("all composition-transport judgments passed")
