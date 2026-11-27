#!/usr/bin/env python3
"""5.235 independent regression for the 3D level-symmetric DO arm.

Nothing here trusts a number printed by the solver: every expected value is
re-derived from the case constants (kappa_a, kappa_s, sigmaSB, wall
temperatures/emissivities, medium profile, direction-set data) with Python.

Arms (the base case is the degenerate one; the others are declared mutations):

  A1 degenerate: nx x 1 x 1 box, symmetry/empty transverse faces, black
     walls, kappa_s = 0. The 3D upwind sweep then reduces EXACTLY to the
     1D slab discrete recursion
         I_i = (alpha I_{i-1} + kappa_a I_b) / (beta + alpha),
         alpha = |s_x| / dx
     for every ordinate of the declared S2 set; the checker rebuilds G, q,
     S cell by cell from that recursion and demands 1e-10 agreement. This
     is the "extension does not change the physics" arm: the same sweep on
     a 1-cell-thick transverse box IS the 1D slab scheme.

  A2 Rosseland: kappa_a >> kappa_s (beta L ~ 400), linear medium profile,
     walls matched to the medium ends. Interior cell flux must match the
     Rosseland diffusion closed form
         q = -(4 sigma / (3 beta)) d(T^4)/dx
     at the cell centre. The wall fluxes themselves are NOT comparable:
     the O(1/beta) boundary layer is unresolved by design, and the
     comparison is declared on the interior cell at 5%.

  A3 isothermal cavity + A4 energy conservation: all six walls grey at one
     temperature, uniform medium at the same temperature. The source
     balance I = I_b is an exact fixed point, so every cell must read
     G = 4 sigma T^4 (this also proves sum(w) = 4 pi of the set), every
     net flux and every volumetric source vanish, and the six-wall flux
     integral minus the volume source integral (the printed balance)
     closes to 1e-10 of the sigma T^4 scale.

  A5 refusals: order not 2/4, emissivity outside [0,1], negative kappa_a,
     negative kappa_s -- each must fail loudly.

usage: check_3d_radiation_do.py <binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
fails = []
nchecks = 0


def check(cond, msg):
    global nchecks
    nchecks += 1
    if not cond:
        fails.append(msg)


def close(a, b, rel, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


def scalar(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", txt)
    if not m:
        raise RuntimeError("no scalar " + key)
    return float(m.group(1))


def block(txt, name):
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", txt)
    if not m:
        raise RuntimeError("no " + name + " block")
    start = m.end() - 1
    depth = 0
    for j in range(start, len(txt)):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
            if depth == 0:
                return txt[start + 1:j]
    raise RuntimeError("unbalanced " + name + " block")


def vectors(block_text, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(\s*(" + NUM + r")\s+("
                  + NUM + r")\s+(" + NUM + r")\s*\)", block_text)
    if not m:
        raise RuntimeError("no vector " + key)
    return tuple(float(m.group(i)) for i in (1, 2, 3))


def ordinates(order):
    """The declared level-symmetric set, rebuilt from published geometry."""
    out = []
    if order == 2:
        c = 1.0 / math.sqrt(3.0)
        for sx in (-1, 1):
            for sy in (-1, 1):
                for sz in (-1, 1):
                    out.append(((sx * c, sy * c, sz * c), math.pi / 2.0))
    elif order == 4:
        n = 1.0 / math.sqrt(6.0)
        for p in ((2, 1, 1), (1, 2, 1), (1, 1, 2)):
            for sx in (-1, 1):
                for sy in (-1, 1):
                    for sz in (-1, 1):
                        out.append(((sx * p[0] * n, sy * p[1] * n, sz * p[2] * n),
                                    math.pi / 6.0))
    else:
        raise RuntimeError("unsupported order %s" % order)
    return out


def parse_cells(output):
    cells = {}
    for line in output.splitlines():
        m = re.search(r"RADDO3D cell (\d+) (\d+) (\d+) x (" + NUM + r") y ("
                      + NUM + r") z (" + NUM + r") G (" + NUM + r") q ("
                      + NUM + r") (" + NUM + r") (" + NUM + r") S (" + NUM + r")", line)
        if m:
            cells[(int(m.group(1)), int(m.group(2)), int(m.group(3)))] = {
                "G": float(m.group(7)),
                "q": (float(m.group(8)), float(m.group(9)), float(m.group(10))),
                "S": float(m.group(11)),
            }
    return cells


def parse_summary(output):
    m = re.search(r"RADDO3D summary qX0 (" + NUM + r") qX1 (" + NUM + r") qY0 ("
                  + NUM + r") qY1 (" + NUM + r") qZ0 (" + NUM + r") qZ1 ("
                  + NUM + r") volSrc (" + NUM + r") balance (" + NUM + r") Gmean ("
                  + NUM + r") Gmin (" + NUM + r") Gmax (" + NUM + r")", output)
    if not m:
        return None
    keys = ["qX0", "qX1", "qY0", "qY1", "qZ0", "qZ1", "volSrc", "balance",
            "Gmean", "Gmin", "Gmax"]
    return {k: float(m.group(i + 1)) for i, k in enumerate(keys)}


def parse_config(output):
    m = re.search(r"RADDO3D config order (\d+) dirs (\d+) cells (\d+) (\d+) (\d+) "
                  r"length (" + NUM + r") (" + NUM + r") (" + NUM + r") iterations (\d+) "
                  r"residual (" + NUM + r")", output)
    if not m:
        return None
    return {"order": int(m.group(1)), "dirs": int(m.group(2)),
            "cells": (int(m.group(3)), int(m.group(4)), int(m.group(5))),
            "length": (float(m.group(6)), float(m.group(7)), float(m.group(8))),
            "iterations": int(m.group(9)), "residual": float(m.group(10))}


def run(text, td, tag):
    p = os.path.join(td, tag + ".dict")
    open(p, "w").write(text)
    return subprocess.run([BINARY, p], capture_output=True, text=True, cwd=td)


def mutate(base, pairs, tag):
    txt = base
    for old, new in pairs:
        check(old in txt, "%s: pattern not found in case: %r" % (tag, old))
        txt = txt.replace(old, new, 1)
    return txt


if len(sys.argv) != 3:
    raise SystemExit("usage: check_3d_radiation_do.py <binary> <case.dict>")
BINARY, CASE = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
BASE = "\n".join(l[:l.find("//")] if "//" in l else l
                 for l in open(CASE).read().splitlines())

rad = block(BASE, "radiation3d")
ka = scalar(rad, "absorption")
ks = scalar(rad, "scattering")
sigma = scalar(rad, "sigmaSB")
T1 = scalar(rad, "mediumT1")
T2 = scalar(rad, "mediumT2")
order = int(re.search(r"\border (\d+)\s*;", rad).group(1))
blk = block(BASE, "block")
mn = vectors(blk, "min")
mx = vectors(blk, "max")
nc = vectors(blk, "cells")
Lx, Ly, Lz = mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]
Nx, Ny, Nz = int(nc[0]), int(nc[1]), int(nc[2])
patches = block(BASE, "patches")
TL = scalar(block(patches, "left"), "wallTemperature")
TR = scalar(block(patches, "right"), "wallTemperature")

ords = ordinates(order)
check(abs(sum(w for _, w in ords) - 4.0 * math.pi) < 1e-12,
      "direction set weights must sum to 4 pi, got %.15g" % sum(w for _, w in ords))

# ------------------------------------------------------------------ A1 degenerate
r = subprocess.run([BINARY, CASE], capture_output=True, text=True,
                   cwd=os.path.dirname(os.path.abspath(CASE)))
check(r.returncode == 0, "A1 solver returned %d: %s" %
      (r.returncode, (r.stdout + r.stderr)[-500:]))
cfg = parse_config(r.stdout)
cells = parse_cells(r.stdout)
summ = parse_summary(r.stdout)
check(cfg is not None and summ is not None, "A1: missing RADDO3D lines")
if cfg:
    check(cfg["order"] == order and cfg["dirs"] == len(ords)
          and cfg["cells"] == (Nx, Ny, Nz),
          "A1: config line disagrees with the case: %s" % cfg)
check(len(cells) == Nx * Ny * Nz, "A1: parsed %d cells, expected %d"
      % (len(cells), Nx * Ny * Nz))

# rebuild the 1D upwind recursion per ordinate (kappa_s = 0 => direct)
beta = ka + ks
dx = Lx / Nx
IL = sigma * TL ** 4 / math.pi
IR = sigma * TR ** 4 / math.pi
Ib = sigma * T1 ** 4 / math.pi
check(ks == 0.0, "A1 requires scattering 0 in the base case")
G_a = [0.0] * Nx
q_a = [0.0] * Nx
for (s, w) in ords:
    a = s[0]
    alpha = abs(a) / dx
    I = [0.0] * Nx
    prev = IL if a > 0 else IR
    rng = range(Nx) if a > 0 else range(Nx - 1, -1, -1)
    for i in rng:
        I[i] = (alpha * prev + ka * Ib) / (beta + alpha)
        prev = I[i]
    for i in range(Nx):
        G_a[i] += w * I[i]
        q_a[i] += w * a * I[i]
S_a = [ka * (g - 4.0 * math.pi * Ib) for g in G_a]

worst = 0.0
for i in range(Nx):
    c = cells.get((i, 0, 0))
    if c is None:
        check(False, "A1: cell %d missing" % i)
        continue
    for key, exact in (("G", G_a[i]), ("S", S_a[i])):
        rel = abs(c[key] - exact) / max(abs(exact), 1.0)
        worst = max(worst, rel)
        check(close(c[key], exact, 1e-10, 1e-9),
              "A1 cell %d %s %g != 1D recursion %g (rel %.2e)"
              % (i, key, c[key], exact, rel))
    relq = abs(c["q"][0] - q_a[i]) / max(abs(q_a[i]), 1.0)
    worst = max(worst, relq)
    check(close(c["q"][0], q_a[i], 1e-10, 1e-9),
          "A1 cell %d qx %g != %g" % (i, c["q"][0], q_a[i]))
    check(abs(c["q"][1]) == 0.0 and abs(c["q"][2]) == 0.0,
          "A1 cell %d has transverse flux on the degenerate box" % i)

# the discrete balance of the same recursion: q(0)-q(L) = integral S dx
if summ:
    scale = max(abs(summ["qX0"]), abs(summ["qX1"]), abs(4 * sigma * TL ** 4), 1.0)
    # A1 box: only the x walls are open, so the balance is
    # qX0*LyLz + qX1*LyLz - integral S dV
    check(abs(summ["balance"]) / (scale * Ly * Lz) <= 1e-10,
          "A1 discrete balance does not close: %g (scale %g)"
          % (summ["balance"], scale * Ly * Lz))
    check(summ["qY0"] == 0 and summ["qY1"] == 0
          and summ["qZ0"] == 0 and summ["qZ1"] == 0,
          "A1 symmetry/empty faces must report exactly zero net flux")

# ------------------------------------------------------------------ A2 Rosseland
ross = mutate(BASE, [
    ("cells (20 1 1);", "cells (40 1 1);"),
    ("absorption 1.0;", "absorption 400;"),
    ("scattering 0;", "scattering 1;"),
    ("mediumT2 1200;", "mediumT2 900;"),
    ("left  { wallTemperature 1500; emissivity 1.0; }",
     "left  { wallTemperature 1200; emissivity 1.0; }"),
], "A2")
with tempfile.TemporaryDirectory(prefix="cf_do3d_") as td:
    r2 = run(ross, td, "a2")
    check(r2.returncode == 0, "A2 solver returned %d: %s" %
          (r2.returncode, (r2.stdout + r2.stderr)[-500:]))
    c2 = parse_cells(r2.stdout)
    s2 = parse_summary(r2.stdout)
    check(len(c2) == 40, "A2: expected 40 cells, parsed %d" % len(c2))
    if c2 and s2:
        beta2 = 400.0 + 1.0
        # local Rosseland flux at the centre cell: T linear 1200 -> 900
        Tc = 0.5 * (1200.0 + 900.0)
        dTdx = (900.0 - 1200.0) / 1.0
        q_ross = -(4.0 * sigma / (3.0 * beta2)) * 4.0 * Tc ** 3 * dTdx
        qc = 0.5 * (c2[(19, 0, 0)]["q"][0] + c2[(20, 0, 0)]["q"][0])
        check(abs(qc - q_ross) <= 0.05 * abs(q_ross),
              "A2 centre flux %g is not within 5%% of Rosseland %g" % (qc, q_ross))
        check(qc > 0, "A2 flux must run hot (x=0) to cold (x=L), got %g" % qc)
        # monotone decreasing G from the hot end (interior, past the layers)
        gs = [c2[(i, 0, 0)]["G"] for i in range(4, 36)]
        check(all(gs[i] > gs[i + 1] for i in range(len(gs) - 1)),
              "A2 G is not monotonically decreasing in the interior")

# ------------------------------------------------- A3/A4 isothermal cavity
iso = mutate(BASE, [
    ("ymin  { name bottom; type symmetry; }", "ymin  { name bottom; type wall; }"),
    ("ymax  { name top;    type symmetry; }", "ymax  { name top;    type wall; }"),
    ("zmin  { name front;  type empty;    }", "zmin  { name front;  type wall;    }"),
    ("zmax  { name back;   type empty;    }", "zmax  { name back;   type wall;    }"),
    ("left  { wallTemperature 1500; emissivity 1.0; }",
     "left  { wallTemperature 1000; emissivity 0.7; }\n"
     "    right { wallTemperature 1000; emissivity 0.7; }\n"
     "    bottom { wallTemperature 1000; emissivity 0.7; }\n"
     "    top { wallTemperature 1000; emissivity 0.7; }\n"
     "    front { wallTemperature 1000; emissivity 0.7; }\n"
     "    back { wallTemperature 1000; emissivity 0.7; }"),
    ("right { wallTemperature 900;  emissivity 1.0; }", ""),
    ("scattering 0;", "scattering 1.5;"),
    ("mediumT1 1200;", "mediumT1 1000;"),
    ("mediumT2 1200;", "mediumT2 1000;"),
], "A3")
with tempfile.TemporaryDirectory(prefix="cf_do3d_") as td:
    r3 = run(iso, td, "a3")
    check(r3.returncode == 0, "A3 solver returned %d: %s" %
          (r3.returncode, (r3.stdout + r3.stderr)[-500:]))
    c3 = parse_cells(r3.stdout)
    s3 = parse_summary(r3.stdout)
    if c3 and s3:
        G_exact = 4.0 * sigma * 1000.0 ** 4
        # every cell, every direction: I == Ib => G == 4 sigma T^4 exactly
        check(close(s3["Gmean"], G_exact, 1e-10, 1e-6),
              "A3 Gmean %g != 4 sigma T^4 %g" % (s3["Gmean"], G_exact))
        check(close(s3["Gmin"], s3["Gmax"], 1e-12, 1e-9),
              "A3 isothermal cavity must read uniform G: %g .. %g"
              % (s3["Gmin"], s3["Gmax"]))
        worstS = max(abs(c["S"]) for c in c3.values())
        worstq = max(max(abs(v) for v in c["q"]) for c in c3.values())
        check(worstS / G_exact <= 1e-10,
              "A3 per-cell source not zero: %g (G scale %g)" % (worstS, G_exact))
        check(worstq / G_exact <= 1e-10,
              "A3 per-cell net flux divergence not zero: max |q| %g" % worstq)
        wallsum = sum(abs(s3[k]) for k in ("qX0", "qX1", "qY0", "qY1", "qZ0", "qZ1"))
        check(wallsum / G_exact <= 1e-10,
              "A3 closed isothermal cavity wall fluxes do not vanish: %g" % wallsum)
        area = max(Ly * Lz, Lx * Lz, Lx * Ly)
        check(abs(s3["balance"]) / (G_exact * area) <= 1e-10,
              "A4 energy balance does not close: %g" % s3["balance"])

# ------------------------------------------------------------------ A5 refusals
refusals = [
    (mutate(BASE, [("order 2;", "order 3;")], "r1"), "order must be 2 or 4"),
    (mutate(BASE, [("order 2;", "order 5;")], "r2"), "order must be 2 or 4"),
    (mutate(BASE, [("emissivity 1.0; }", "emissivity 1.5; }")], "r3"),
     "emissivity must lie in [0,1]"),
    (mutate(BASE, [("emissivity 1.0; }", "emissivity -0.1; }")], "r4"),
     "emissivity must lie in [0,1]"),
    (mutate(BASE, [("absorption 1.0;", "absorption -1.0;")], "r5"),
     "absorption kappa_a cannot be negative"),
    (mutate(BASE, [("scattering 0;", "scattering -2.0;")], "r6"),
     "scattering kappa_s cannot be negative"),
]
with tempfile.TemporaryDirectory(prefix="cf_do3d_rej_") as td:
    for txt, needle in refusals:
        rr = run(txt, td, "rej")
        check(rr.returncode != 0 and needle in (rr.stdout + rr.stderr),
              "refusal '%s' not raised loudly (rc=%d)" % (needle, rr.returncode))

# ------------------------------------------------------------------ verdict
if fails:
    print("radiation DO 3D FAILED %d checks:" % len(fails))
    for f in fails[:40]:
        print("  -", f)
    sys.exit(1)

print("radiation DO 3D OK: %d assertions; degenerate S%d box == 1D slab "
      "recursion (worst %.2e rel), Rosseland centre flux within 5%%, "
      "isothermal cavity G == 4 sigma T^4 with zero fluxes (<1e-10), "
      "closed-cavity balance < 1e-10, %d loud refusals"
      % (nchecks, order, worst, len(refusals)))
