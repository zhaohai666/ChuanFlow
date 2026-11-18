#!/usr/bin/env python3
"""D8 increment 4: verify the steady P1 radiation solver against a 1D grey slab
whose answer is a published closed form, not a number the solver was given.

The case is a slab between two black walls at different temperatures with a
uniform participating medium. The steady P1 equation

    gamma G'' - kappa_a G + 4 kappa_a sigma T_m^4 = 0,   gamma = 1/(3 beta)

has the exact solution  G(x) = 4 sigma T_m^4 + C1 cosh(m x) + C2 sinh(m x),
m^2 = kappa_a/gamma = 3 beta kappa_a, with C1, C2 fixed by the two black-wall
values G(0) = 4 sigma T_L^4, G(L) = 4 sigma T_R^4. Every constant is rebuilt
from the case file; the solver is never told the closed form.

  * analytic profile  - cell-centre G against the double-exponential, loose
       tolerance because a cell-centred FV scheme has O(dx^2) error against a
       curved continuous profile (the tolerance is a fraction of the wall jump);
  * black walls       - the printed wall G equals 4 sigma T_w^4 to 1e-9;
  * exact balance     - integral(S_rad) + surface(q_r.n) = 0 to 1e-6 relative:
       this is a discrete identity of the assembled operator, so it is sharp;
  * flux law          - q_r = -gamma grad G rebuilt by central differences;
  * net source        - S_rad = kappa_a (G - 4 sigma T_m^4), tight;
  * monotonicity      - G decreases from the hot wall to the cold wall;
  * refusals          - negative absorption / scattering, missing block, zero
       medium temperature fail loudly.

usage: check_p1_slab.py <solver-binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

SIGMA = 5.670374419e-8
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


def block(txt, name):
    """Body of a named { } block, balancing nested braces."""
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


def run(path, tag, td):
    p = os.path.join(td, tag + ".dict")
    open(p, "w").write(path)
    return subprocess.run([BINARY, p], capture_output=True, text=True)


def parse_cells(output):
    cells = {}
    for line in output.splitlines():
        m = re.search(r"p1 cell (\d+): G (" + NUM + r") q (" + NUM + r") ("
                      + NUM + r") (" + NUM + r") S (" + NUM + r")", line)
        if m:
            cells[int(m.group(1))] = {
                "G": float(m.group(2)),
                "q": (float(m.group(3)), float(m.group(4)), float(m.group(5))),
                "S": float(m.group(6)),
            }
    return cells


def parse_balance(output):
    m = re.search(r"p1 balance: volumeSource (" + NUM + r") boundaryFlux ("
                  + NUM + r") sum (" + NUM + r")", output)
    return (float(m.group(1)), float(m.group(2)), float(m.group(3))) if m else None


def parse_walls(output):
    w = {}
    for line in output.splitlines():
        m = re.search(r"p1 wall (\w+): T (" + NUM + r") K -> G (" + NUM + r")", line)
        if m:
            w[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    return w


if len(sys.argv) != 3:
    raise SystemExit("usage: check_p1_slab.py <solver-binary> <case.dict>")
BINARY = sys.argv[1]
CASE = sys.argv[2]
BASE = open(CASE).read()

rad = block(BASE, "radiation")
ka = scalar(rad, "absorption")
ks = scalar(rad, "scattering")
Tm = scalar(rad, "mediumTemperature")
blk = block(BASE, "block")
mn = vectors(blk, "min")
mx = vectors(blk, "max")
nc = vectors(blk, "cells")
Lx = mx[0] - mn[0]
Nx = int(nc[0])
dx = Lx / Nx
beta = ka + ks
gamma = 1.0 / (3.0 * beta)
m2 = ka / gamma
mm = math.sqrt(m2)

patches = block(BASE, "patches")
wls = {}
for name in ("left", "right"):
    pb = block(patches, name)
    wls[name] = scalar(pb, "wallTemperature")
TL, TR = wls["left"], wls["right"]

Gm = 4 * SIGMA * Tm ** 4
G0 = 4 * SIGMA * TL ** 4
GL = 4 * SIGMA * TR ** 4

r = subprocess.run([BINARY, CASE], capture_output=True, text=True)
check(r.returncode == 0, "solver returned %d: %s" % (r.returncode, (r.stdout + r.stderr)[-400:]))
out = r.stdout
cells = parse_cells(out)
bal = parse_balance(out)
walls = parse_walls(out)

check(len(cells) == Nx, "expected %d cells, parsed %d" % (Nx, len(cells)))
check(bal is not None, "no balance line")

# ---- black walls -------------------------------------------------------
for name, (tw, gw) in walls.items():
    check(close(gw, 4 * SIGMA * tw ** 4, 1e-9),
          "wall %s: G %g != 4 sigma T^4 %g" % (name, gw, 4 * SIGMA * tw ** 4))
check(set(walls.keys()) == {"left", "right"}, "expected left+right walls, got %s" % list(walls))

# ---- analytic double-exponential ---------------------------------------
C1 = G0 - Gm
C2 = (GL - Gm - C1 * math.cosh(mm * Lx)) / math.sinh(mm * Lx)
jump = abs(G0 - GL)
tol_abs = 3e-3 * jump  # O(dx^2) discretisation against a curved continuous profile
worst_rel = 0.0
for i in sorted(cells):
    xc = mn[0] + (i + 0.5) * dx
    Ga = Gm + C1 * math.cosh(mm * xc) + C2 * math.sinh(mm * xc)
    Gn = cells[i]["G"]
    err = abs(Gn - Ga)
    worst_rel = max(worst_rel, err / max(Ga, 1.0))
    check(err <= tol_abs, "cell %d: G %g vs analytic %g (err %g > %g)"
          % (i, Gn, Ga, err, tol_abs))

# ---- monotonic decrease ------------------------------------------------
gseq = [cells[i]["G"] for i in sorted(cells)]
check(all(gseq[i] > gseq[i + 1] for i in range(len(gseq) - 1)),
      "G is not monotonically decreasing from hot to cold wall")

# ---- net source law ----------------------------------------------------
for i, c in cells.items():
    S_expected = ka * (c["G"] - Gm)
    check(close(c["S"], S_expected, 1e-9, 1.0),
          "cell %d: S %g != kappa(G-4sigmaT4) %g" % (i, c["S"], S_expected))

# ---- radiative flux law (central difference on the printed field) ------
# q_x = -gamma dG/dx; interior cells only (one-sided would add O(dx) error)
idx = sorted(cells)
for j in range(1, len(idx) - 1):
    i = idx[j]
    dGdx = (cells[i + 1]["G"] - cells[i - 1]["G"]) / (2 * dx)
    qx = -gamma * dGdx
    check(close(cells[i]["q"][0], qx, 8e-3, 200.0),
          "cell %d: qx %g != -gamma gradG %g" % (i, cells[i]["q"][0], qx))
    check(abs(cells[i]["q"][1]) < 1e-4 and abs(cells[i]["q"][2]) < 1e-4,
          "cell %d: transverse flux not zero on a 1D slab" % i)

# ---- exact discrete balance --------------------------------------------
vol_src, bnd_flux, total = bal
scale = max(abs(vol_src), abs(bnd_flux), 1.0)
check(abs(total) / scale <= 1e-6,
      "balance does not close: vol %g + bnd %g = %g" % (vol_src, bnd_flux, total))
check(vol_src * bnd_flux < 0, "volume source and boundary flux must have opposite signs")

# ---- refusals -----------------------------------------------------------
neg_ka = BASE.replace("absorption 1.0;", "absorption -1.0;")
neg_ks = BASE.replace("scattering 0.5;", "scattering -0.5;")
zero_T = BASE.replace("mediumTemperature 1200;", "mediumTemperature 0;")
no_block = re.sub(r"radiation\s*\{.*?^\}", "", BASE, count=1, flags=re.M | re.S)
for txt, needle, tag in [
    (neg_ka, "absorption kappa_a must be a finite positive", "neg_ka"),
    (neg_ks, "scattering kappa_s cannot be negative", "neg_ks"),
    (zero_T, "mediumTemperature must be a finite positive", "zero_T"),
    (no_block, "needs a radiation", "no_block"),
]:
    with tempfile.TemporaryDirectory(prefix="cf_p1_") as td:
        rr = run(txt, tag, td)
        check(rr.returncode != 0 and needle in (rr.stdout + rr.stderr),
              "refusal '%s' not raised loudly (rc=%d)" % (needle, rr.returncode))

if fails:
    print("p1_slab FAILED %d checks:" % len(fails))
    for f in fails[:30]:
        print("  -", f)
    sys.exit(1)

print("p1_slab OK: double-exponential closed form (worst %.2e relative), black-wall "
      "G=4sigmaT^4, exact discrete balance (%.1e), flux/source laws, monotonic, "
      "4 loud refusals" % (worst_rel, abs(total) / scale))
