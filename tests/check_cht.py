#!/usr/bin/env python3
"""D8 increment 7: verify two-region conjugate heat transfer against the
one-dimensional two-layer series-resistance closed form.

The two regions have different conductivities and outer wall temperatures; the
interface temperature and heat flux are found by a Dirichlet--Neumann iteration,
not prescribed. From the case constants alone the checker rebuilds:

    R   = L_s/k_s + L_f/k_f,   q = (T_L - T_R)/R,
    T_i = T_L - q L_s/k_s,

then the linear profile in each layer at the cell centres, and requires the two
interface fluxes to be equal and opposite (flux continuity) to a tight tolerance.

  * interface T and flux - against the closed form, tight;
  * flux continuity      - |q_s + q_f| / q small (the CHT coupling residual);
  * linear profiles      - every cell-centre T against the closed form;
  * outer flux           - each region's outer face carries the same |q|;
  * few iterations       - the secant converges in a handful of sweeps;
  * refusals             - negative k, equal outer T, bad relax, missing block.

usage: check_cht.py <solver-binary> <case.dict>
"""
import os
import re
import subprocess
import sys
import tempfile

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
    return tuple(float(m.group(i)) for i in (1, 2, 3))


if len(sys.argv) != 3:
    raise SystemExit("usage: check_cht.py <solver-binary> <case.dict>")
BINARY, CASE = sys.argv[1], sys.argv[2]
BASE = open(CASE).read()
cd = block(BASE, "cht")
sd = block(cd, "solid")
fd = block(cd, "fluid")

ks, kf = scalar(sd, "k"), scalar(fd, "k")
TL, TR = scalar(sd, "outerT"), scalar(fd, "outerT")
sb = block(sd, "block")
fb = block(fd, "block")
smin, smax, snc = vectors(sb, "min"), vectors(sb, "max"), vectors(sb, "cells")
fmin, fmax, fnc = vectors(fb, "min"), vectors(fb, "max"), vectors(fb, "cells")
Ls = smax[0] - smin[0]
Lf = fmax[0] - fmin[0]
Ns, Nf = int(snc[0]), int(fnc[0])

R = Ls / ks + Lf / kf
q = (TL - TR) / R
Ti = TL - q * Ls / ks

r = subprocess.run([BINARY, CASE], capture_output=True, text=True)
check(r.returncode == 0, "solver returned %d: %s" % (r.returncode, (r.stdout + r.stderr)[-400:]))
out = r.stdout

mi = re.search(r"cht interface: T (" + NUM + r") fluxMismatch (" + NUM
               + r") outerIters (\d+)", out)
check(mi is not None, "no interface line")
mf = re.search(r"cht flux: solid outer (" + NUM + r") solid iface (" + NUM
               + r") fluid iface (" + NUM + r") fluid outer (" + NUM + r")", out)
check(mf is not None, "no flux line")
if not (mi and mf):
    print("cht FAILED: missing output lines")
    sys.exit(1)

Ti_num = float(mi.group(1))
mismatch = float(mi.group(2))
outer_iters = int(mi.group(3))
so, si, fi, fo = (float(mf.group(i)) for i in (1, 2, 3, 4))

# ---- interface closed form ---------------------------------------------
check(close(Ti_num, Ti, 1e-8, 1e-6), "interface T %g != closed form %g" % (Ti_num, Ti))
check(mismatch / abs(q) <= 1e-6, "flux mismatch %g too large (q=%g)" % (mismatch, q))
check(outer_iters <= 10, "secant took %d outer iterations" % outer_iters)

# ---- interface flux continuity -----------------------------------------
check(close(si, q, 1e-6, 1e-4), "solid iface flux %g != q %g" % (si, q))
check(close(fi, -q, 1e-6, 1e-4), "fluid iface flux %g != -q %g" % (fi, -q))
check(abs(si + fi) <= 1e-6 * abs(q), "interface flux not continuous: %g + %g" % (si, fi))

# ---- outer flux carries the same |q| -----------------------------------
check(close(abs(so), abs(q), 1e-5, 1e-3), "solid outer |flux| %g != |q| %g" % (abs(so), abs(q)))
check(close(abs(fo), abs(q), 1e-5, 1e-3), "fluid outer |flux| %g != |q| %g" % (abs(fo), abs(q)))
check(so * si < 0 and fo * fi < 0, "outer and interface flux of a layer must oppose")

# ---- linear profiles at cell centres -----------------------------------
solidT, fluidT = {}, {}
for line in out.splitlines():
    a = re.search(r"cht solid cell (\d+): T (" + NUM + r")", line)
    if a:
        solidT[int(a.group(1))] = float(a.group(2))
    b = re.search(r"cht fluid cell (\d+): T (" + NUM + r")", line)
    if b:
        fluidT[int(b.group(1))] = float(b.group(2))
check(len(solidT) == Ns and len(fluidT) == Nf,
      "cell records %d/%d != %d/%d" % (len(solidT), len(fluidT), Ns, Nf))

dxS = Ls / Ns
for i, T in solidT.items():
    xc = smin[0] + (i + 0.5) * dxS
    Ta = TL - q * xc / ks
    check(abs(T - Ta) <= 1e-6, "solid cell %d: T %g != %g" % (i, T, Ta))
dxF = Lf / Nf
for i, T in fluidT.items():
    xc = fmin[0] + (i + 0.5) * dxF
    Ta = Ti - q * (xc - fmin[0]) / kf
    check(abs(T - Ta) <= 1e-6, "fluid cell %d: T %g != %g" % (i, T, Ta))

# ---- refusals -----------------------------------------------------------
refusals = [
    (BASE.replace("k 2.0;", "k -2.0;", 1), "solid k must be finite and positive", "neg_k"),
    (BASE.replace("outerT 500.0;", "outerT 300.0;", 1),
     "two outer temperatures must differ", "equal_T"),
    (BASE.replace("relax 1.0;", "relax 1.5;", 1), "relax must be in", "bad_relax"),
    (re.sub(r"cht\s*\{.*?^\}", "", BASE, count=1, flags=re.M | re.S),
     "needs a cht", "no_block"),
]
for txt, needle, tag in refusals:
    with tempfile.TemporaryDirectory(prefix="cf_cht_") as td:
        p = os.path.join(td, tag + ".dict")
        open(p, "w").write(txt)
        rr = subprocess.run([BINARY, p], capture_output=True, text=True)
        check(rr.returncode != 0 and needle in (rr.stdout + rr.stderr),
              "refusal '%s' not raised loudly (rc=%d): %s"
              % (needle, rr.returncode, (rr.stdout + rr.stderr)[-200:]))

if fails:
    print("cht FAILED %d checks:" % len(fails))
    for f in fails[:30]:
        print("  -", f)
    sys.exit(1)

print("cht OK: series-resistance interface T and flux, flux continuity, linear "
      "profiles in both regions, outer flux, %d outer iters, 4 loud refusals" % outer_iters)
