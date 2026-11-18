#!/usr/bin/env python3
"""D8 increment 5: verify steady species transport + one-step Arrhenius
combustion against a case whose fuel profile is a published closed form.

The case pins the oxidiser and is stagnant, so the fuel obeys
    D Y'' - k Yo Y = 0,   k = A exp(-Ea/(R T)),   p = sqrt(k Yo / D),
with Y(0)=fuelInlet and Y(L)=0 giving
    Y(x) = Y0 sinh(p (L-x)) / sinh(p L).
Every constant is rebuilt from the case with an INDEPENDENT math.exp, so the
Arrhenius temperature law and the reaction-diffusion balance are both graded
without trusting the solver's own rate.

  * arrhenius     - the three printed k(T) match math.exp exactly (1e-9 rel);
  * analytic fuel - cell-centre Yf against the sinh profile (O(dx^2) slack);
  * rate law      - printed rate = k Yf Yo, tight;
  * pinned ox     - Yo stays at fixedOxidiser everywhere;
  * fuel balance  - inflow - outflow = integral of the sink, sharp (discrete
                    identity); the pinned oxidiser is an infinite reservoir and
                    is deliberately NOT required to balance;
  * refusals      - negative diffusivity / preExp / activationE / no block.

usage: check_combustion.py <solver-binary> <case.dict>
"""
import math
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
    raise SystemExit("usage: check_combustion.py <solver-binary> <case.dict>")
BINARY, CASE = sys.argv[1], sys.argv[2]
BASE = open(CASE).read()

cb = block(BASE, "combustion")
D = scalar(cb, "diffusivity")
A = scalar(cb, "preExp")
Ea = scalar(cb, "activationE")
R = scalar(cb, "gasConstant")
T = scalar(cb, "temperature")
Y0 = scalar(cb, "fuelInlet")
Yo = scalar(cb, "fixedOxidiser")

blk = block(BASE, "block")
mn = vectors(blk, "min")
mx = vectors(blk, "max")
nc = vectors(blk, "cells")
Lx = mx[0] - mn[0]
Nx = int(nc[0])
dx = Lx / Nx

r = subprocess.run([BINARY, CASE], capture_output=True, text=True)
check(r.returncode == 0, "solver returned %d: %s" % (r.returncode, (r.stdout + r.stderr)[-400:]))
out = r.stdout

cells = {}
for line in out.splitlines():
    m = re.search(r"combustion cell (\d+): Yf (" + NUM + r") Yo (" + NUM
                  + r") rate (" + NUM + r")", line)
    if m:
        cells[int(m.group(1))] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
check(len(cells) == Nx, "expected %d cells, parsed %d" % (Nx, len(cells)))

# ---- Arrhenius independently -------------------------------------------
k = A * math.exp(-Ea / (R * T))
scan = {}
for line in out.splitlines():
    m = re.search(r"combustion arrhenius: T (" + NUM + r") k (" + NUM + r")", line)
    if m:
        scan[float(m.group(1))] = float(m.group(2))
check(len(scan) == 3, "expected 3 arrhenius scan points, got %d" % len(scan))
for Tp, kp in scan.items():
    kexp = A * math.exp(-Ea / (R * Tp))
    check(close(kp, kexp, 1e-9), "arrhenius k(T=%g)=%g != %g" % (Tp, kp, kexp))

# ---- analytic fuel profile ---------------------------------------------
p = math.sqrt(k * Yo / D)
tol = 3e-3 * Y0
for i, (yf, yo, rate) in cells.items():
    xc = mn[0] + (i + 0.5) * dx
    ya = Y0 * math.sinh(p * (Lx - xc)) / math.sinh(p * Lx)
    check(abs(yf - ya) <= tol, "cell %d: Yf %g vs analytic %g (err %g > %g)"
          % (i, yf, ya, abs(yf - ya), tol))

# ---- rate law and pinned oxidiser --------------------------------------
for i, (yf, yo, rate) in cells.items():
    check(close(rate, k * yf * yo, 1e-9, 1e-12),
          "cell %d: rate %g != k Yf Yo %g" % (i, rate, k * yf * yo))
    check(abs(yo - Yo) < 1e-12, "cell %d: Yo %g != pinned %g" % (i, yo, Yo))

# ---- fuel conservation -------------------------------------------------
mb = re.search(r"combustion balance: fuel (" + NUM + r")", out)
check(mb is not None, "no fuel balance line")
if mb:
    resid = float(mb.group(1))
    # scale by the analytic sink magnitude
    analytic_sink = sum(k * Y0 * math.sinh(p * (Lx - (mn[0] + (i + 0.5) * dx))) /
                        math.sinh(p * Lx) * dx for i in range(Nx)) * Yo
    check(abs(resid) <= 1e-6 * max(analytic_sink, 1e-3),
          "fuel balance residual %g too large (sink ~%g)" % (resid, analytic_sink))

# ---- refusals -----------------------------------------------------------
repls = [
    (BASE.replace("diffusivity 1.0;", "diffusivity -1.0;"),
     "diffusivity must be finite and positive", "neg_D"),
    (BASE.replace("preExp 100.0;", "preExp -100.0;"),
     "preExp A must be finite and positive", "neg_A"),
    (BASE.replace("activationE 50000.0;", "activationE -1.0;"),
     "activationE Ea cannot be negative", "neg_Ea"),
    (re.sub(r"combustion\s*\{.*?^\}", "", BASE, count=1, flags=re.M | re.S),
     "needs a combustion", "no_block"),
]
for txt, needle, tag in repls:
    with tempfile.TemporaryDirectory(prefix="cf_comb_") as td:
        p = os.path.join(td, tag + ".dict")
        open(p, "w").write(txt)
        rr = subprocess.run([BINARY, p], capture_output=True, text=True)
        check(rr.returncode != 0 and needle in (rr.stdout + rr.stderr),
              "refusal '%s' not raised loudly (rc=%d)" % (needle, rr.returncode))

if fails:
    print("combustion FAILED %d checks:" % len(fails))
    for f in fails[:30]:
        print("  -", f)
    sys.exit(1)

print("combustion OK: independent Arrhenius k(T), sinh closed form, rate law, "
      "pinned oxidiser, sharp fuel balance, 4 loud refusals")
