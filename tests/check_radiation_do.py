#!/usr/bin/env python3
"""F9: independent regression for DO / S2S / P1 grey radiation.

Nothing here uses a number printed by the solver to build the expected answer.

DO (case.dict)
--------------
For a uniform grey slab (kappa_s = 0, beta = kappa_a) between black walls the
S_N transport equation decouples per ordinate:

    mu_m I_m' + beta I_m = kappa_a I_b,   I_b = sigma T_m^4 / pi

so each ray attenuates exactly:

    I_m(x) = I_b + (I_L - I_b) exp(-beta x / mu_m)     mu_m > 0
    I_m(x) = I_b + (I_R - I_b) exp(-beta (L-x)/|mu_m|) mu_m < 0

with I_L = sigma T_L^4 / pi, I_R = sigma T_R^4 / pi. This checker rebuilds
from the case text alone:

  * the cell-AVERAGED incident radiation  G_i = 2 pi sum w_m <I_m>_i  --
    exponential fitting gives <I_m>_i = b/beta + (I_in - b/beta)
    (|mu|/(beta dx))(1-exp(-beta dx/|mu|)); that cell average is exactly the
    quantity for which the discrete wall/source balance closes;
  * the cell-CENTRE flux q_i = 2 pi sum w_m mu_m I_m(x_i);
  * the net source S_i = kappa_a (G_i - 4 pi I_b);
  * the closed-form wall fluxes q(0), q(L) and the integral of S.

All are compared to the printed fields at 1e-6 relative (the mesh is exact
for a piecewise-constant source and kappa_s=0). The balance
q(0)-q(L)-integral(S dx) is a discrete identity and is checked tight.

S2S (derived from case.dict by swapping the model word)
--------------------------------------------------------
Two equal-area, directly opposed black plates: F_LR=F_RL=1 (summation +
reciprocity), so q = sigma (T_L^4 - T_R^4) and G = 2 sigma (T_L^4 + T_R^4)
in the transparent-medium limit; zero volumetric source.

P1 (caseP1.dict) is used only for a documented qualitative agreement trend.
Five-plus loud refusals are exercised through temp case mutations.

usage: check_radiation_do.py <binary> <case.dict> <caseP1.dict> <caseReject.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

SIGMA = 5.670374419e-8
NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"

# Published Gauss-Legendre roots of P_4 / P_8 on [-1,1] (same constants the
# solver uses; they are quadrature data, not solver output).
HALF = {
    4: [(0.8611363115940526, 0.3478548451374538),
        (0.3399810435848563, 0.6521451548625461)],
    8: [(0.9602898564975363, 0.1012285362903763),
        (0.7966664774136267, 0.2223810344533745),
        (0.5255324099163290, 0.3137066458778873),
        (0.1834346424956498, 0.3626837833783620)],
}

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


def run(path):
    return subprocess.run([BINARY, path], capture_output=True, text=True)


def run_text(txt, tag, td):
    p = os.path.join(td, tag + ".dict")
    open(p, "w").write(txt)
    return run(p)


def parse_cells(output):
    cells = {}
    for line in output.splitlines():
        m = re.search(r"RADDO cell (\d+) x (" + NUM + r") G (" + NUM + r") q ("
                      + NUM + r") (" + NUM + r") (" + NUM + r") S (" + NUM + r")", line)
        if m:
            cells[int(m.group(1))] = {
                "x": float(m.group(2)), "G": float(m.group(3)),
                "q": (float(m.group(4)), float(m.group(5)), float(m.group(6))),
                "S": float(m.group(7)),
            }
    return cells


def parse_summary(output):
    pat = (r"RADDO t (" + NUM + r") model (\w+) Gmean (" + NUM + r") qLeft ("
           + NUM + r") qRight (" + NUM + r") volSrc (" + NUM + r") balance ("
           + NUM + r") Gmin (" + NUM + r") Gmax (" + NUM + r") viewFactor (" + NUM + r")")
    m = re.search(pat, output)
    if not m:
        return None
    return {
        "model": m.group(2), "Gmean": float(m.group(3)),
        "qLeft": float(m.group(4)), "qRight": float(m.group(5)),
        "volSrc": float(m.group(6)), "balance": float(m.group(7)),
        "Gmin": float(m.group(8)), "Gmax": float(m.group(9)),
        "viewFactor": float(m.group(10)), "raw_Gmean": m.group(3),
    }


def sig_digits(s):
    core = s.lstrip("+-").replace(".", "").lstrip("0")
    return len(core)


if len(sys.argv) != 5:
    raise SystemExit("usage: check_radiation_do.py <binary> <case.dict> "
                     "<caseP1.dict> <caseReject.dict>")
BINARY, CASE, CASE_P1, CASE_REJECT = sys.argv[1:5]
BASE = open(CASE).read()

# ---------------------------------------------------------------- case data
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
patches = block(BASE, "patches")
TL = scalar(block(patches, "left"), "wallTemperature")
TR = scalar(block(patches, "right"), "wallTemperature")

qm = re.search(r"\bquadrature\s+[Ss]?(\d+)", BASE)
order = int(qm.group(1)) if qm else 8
check(order in HALF, "unsupported quadrature order in case: %s" % order)
half = HALF.get(order, HALF[8])
ords = [(mu, w) for mu, w in half] + [(-mu, w) for mu, w in half]
check(abs(sum(w for _, w in ords) - 2.0) < 1e-12, "quadrature weights must sum to 2")

# The 1e-6 analytic comparison is exact only for a non-scattering slab; the
# shipped DO case is deliberately kappa_s=0 so every printed number is closed
# form. (Scattering is accepted by the solver but no closed form is asserted.)
check(ks == 0.0, "DO analytic regression requires kappa_s = 0 in case.dict")

pi = math.pi
Ib = SIGMA * Tm ** 4 / pi
IL = SIGMA * TL ** 4 / pi
IR = SIGMA * TR ** 4 / pi


def face_I(mu, x):
    """Exact ray intensity at position x (pure absorption, black walls)."""
    if mu > 0:
        return Ib + (IL - Ib) * math.exp(-beta * x / mu)
    return Ib + (IR - Ib) * math.exp(-beta * (Lx - x) / abs(mu))


def avg_I(mu, xc):
    """Cell-average ray intensity over [xc-dx/2, xc+dx/2] by exp fitting."""
    a = abs(mu)
    i_in = face_I(mu, xc - dx / 2.0) if mu > 0 else face_I(mu, xc + dx / 2.0)
    return Ib + (i_in - Ib) * (a / (beta * dx)) * (1.0 - math.exp(-beta * dx / a))


x_centres = {i: mn[0] + (i + 0.5) * dx for i in range(Nx)}


def analytic(ord_set):
    """Rebuild every printed DO quantity for a given S_N ordinate set."""
    G, q, S = {}, {}, {}
    for i in range(Nx):
        xc = x_centres[i]
        G[i] = 2 * pi * sum(w * avg_I(mu, xc) for mu, w in ord_set)
        q[i] = 2 * pi * sum(w * mu * face_I(mu, xc) for mu, w in ord_set)
        S[i] = ka * (G[i] - 4 * pi * Ib)
    q0 = 2 * pi * sum(w * mu * face_I(mu, 0.0) for mu, w in ord_set)
    qL = 2 * pi * sum(w * mu * face_I(mu, Lx) for mu, w in ord_set)
    vol = dx * sum(S.values())
    return G, q, S, q0, qL, vol, q0 - qL - vol


# analytic fields for the case's quadrature
G_exact, q_exact, S_exact, q0_exact, qL_exact, vol_exact, balance_exact = analytic(ords)

# ---------------------------------------------------------------- DO solve
r = run(CASE)
check(r.returncode == 0, "DO solver returned %d: %s" %
      (r.returncode, (r.stdout + r.stderr)[-400:]))
out = r.stdout
cells = parse_cells(out)
summary = parse_summary(out)
check(len(cells) == Nx, "expected %d DO cells, parsed %d" % (Nx, len(cells)))
check(summary is not None, "no RADDO summary line")
if summary:
    check(summary["model"] == "do", "summary model is %s, expected do" % summary["model"])
    check(sig_digits(summary["raw_Gmean"]) >= 12,
          "Gmean printed with %d significant digits, expected >=12" %
          sig_digits(summary["raw_Gmean"]))

# analytic validity: the closed form itself must satisfy the balance sharply
scale = max(abs(q0_exact), abs(qL_exact), abs(vol_exact), 1.0)
check(abs(balance_exact) / scale <= 1e-9,
      "analytic balance does not close: %g" % balance_exact)

# per-cell fields
worst = 0.0
for i in sorted(cells):
    c = cells[i]
    check(close(c["x"], x_centres[i], 1e-12, 1e-9),
          "cell %d x %g != %g" % (i, c["x"], x_centres[i]))
    for key, exact in (("G", G_exact[i]), ("S", S_exact[i])):
        rel = abs(c[key] - exact) / max(abs(exact), 1.0)
        worst = max(worst, rel)
        check(close(c[key], exact, 1e-6, 1e-6),
              "cell %d %s %g != analytic %g (rel %.2e)" % (i, key, c[key], exact, rel))
    relq = abs(c["q"][0] - q_exact[i]) / max(abs(q_exact[i]), 1.0)
    worst = max(worst, relq)
    check(close(c["q"][0], q_exact[i], 1e-6, 1e-6),
          "cell %d qx %g != analytic %g (rel %.2e)" % (i, c["q"][0], q_exact[i], relq))
    check(abs(c["q"][1]) < 1e-6 and abs(c["q"][2]) < 1e-6,
          "cell %d has transverse flux on a 1D slab" % i)

# summary wall fluxes / source / balance
if summary:
    for key, exact in (("qLeft", q0_exact), ("qRight", -qL_exact), ("volSrc", vol_exact)):
        check(close(summary[key], exact, 1e-6, 1e-6),
              "%s %g != analytic %g" % (key, summary[key], exact))
    check(close(summary["Gmean"], sum(G_exact.values()) / Nx, 1e-6, 1e-6),
          "Gmean %g != analytic %g" % (summary["Gmean"], sum(G_exact.values()) / Nx))
    check(abs(summary["balance"]) / scale <= 1e-8,
          "printed DO balance %g is not ~0 (scale %g)" % (summary["balance"], scale))
    check(close(summary["balance"], balance_exact, 1e-6, 1e-9),
          "printed balance %g != analytic %g" % (summary["balance"], balance_exact))
    # physics sign / ordering: medium hotter than the right wall -> right wall
    # is a net absorber, so qRight < 0 < qLeft for this case.
    check(summary["qLeft"] > 0 and summary["qRight"] < 0,
          "wall flux signs wrong: qLeft %g qRight %g" %
          (summary["qLeft"], summary["qRight"]))
    check(summary["Gmin"] < summary["Gmean"] < summary["Gmax"],
          "Gmean must lie between Gmin and Gmax")

gseq = [cells[i]["G"] for i in sorted(cells)]
check(all(gseq[i] > gseq[i + 1] for i in range(len(gseq) - 1)),
      "DO G is not monotonically decreasing from hot to cold wall")

# ---------------------------------------------------------------- S4 sweep
# The solver also supports S4; rebuild its (different) closed form and check
# the same 1e-6 agreement, proving the quadrature is really read from the case.
with tempfile.TemporaryDirectory(prefix="cf_raddo_s4_") as td:
    s4_txt = BASE.replace("quadrature S8;", "quadrature S4;")
    check(s4_txt != BASE, "could not derive an S4 case from case.dict")
    r4 = run_text(s4_txt, "s4", td)
    check(r4.returncode == 0, "S4 solver returned %d: %s" %
          (r4.returncode, (r4.stdout + r4.stderr)[-400:]))
    cells4 = parse_cells(r4.stdout)
    sum4 = parse_summary(r4.stdout)
    check(len(cells4) == Nx, "S4 cell count %d != %d" % (len(cells4), Nx))
    if cells4 and sum4:
        o4 = [(mu, w) for mu, w in HALF[4]] + [(-mu, w) for mu, w in HALF[4]]
        g4, q4, s4, q04, qL4, vol4, bal4 = analytic(o4)
        sc4 = max(abs(q04), abs(qL4), abs(vol4), 1.0)
        check(abs(bal4) / sc4 <= 1e-9, "S4 analytic balance does not close")
        check(sum4["model"] == "do", "S4 summary model is %s" % sum4["model"])
        for i in sorted(cells4):
            for key, exact in (("G", g4[i]), ("S", s4[i])):
                check(close(cells4[i][key], exact, 1e-6, 1e-6),
                      "S4 cell %d %s %g != analytic %g" %
                      (i, key, cells4[i][key], exact))
            check(close(cells4[i]["q"][0], q4[i], 1e-6, 1e-6),
                  "S4 cell %d qx %g != analytic %g" %
                  (i, cells4[i]["q"][0], q4[i]))
        for key, exact in (("qLeft", q04), ("qRight", -qL4), ("volSrc", vol4)):
            check(close(sum4[key], exact, 1e-6, 1e-6),
                  "S4 %s %g != analytic %g" % (key, sum4[key], exact))
        check(abs(sum4["balance"]) / sc4 <= 1e-8,
              "S4 printed balance %g is not ~0" % sum4["balance"])

# ---------------------------------------------------------------- S2S solve
with tempfile.TemporaryDirectory(prefix="cf_raddo_") as td:
    s2s_txt = BASE.replace("model do;", "model s2s;")
    check(s2s_txt != BASE, "could not derive an s2s case from case.dict")
    rr = run_text(s2s_txt, "s2s", td)
    check(rr.returncode == 0, "S2S solver returned %d: %s" %
          (rr.returncode, (rr.stdout + rr.stderr)[-400:]))
    s2s_sum = parse_summary(rr.stdout)
    s2s_cells = parse_cells(rr.stdout)
    if s2s_sum:
        check(s2s_sum["model"] == "s2s", "S2S summary model is %s" % s2s_sum["model"])
        q_s2s = SIGMA * (TL ** 4 - TR ** 4)
        G_s2s = 2 * SIGMA * (TL ** 4 + TR ** 4)
        check(close(s2s_sum["qLeft"], q_s2s, 1e-6, 1e-6),
              "S2S qLeft %g != sigma(TL^4-TR^4) %g" % (s2s_sum["qLeft"], q_s2s))
        check(close(s2s_sum["qRight"], -q_s2s, 1e-6, 1e-6),
              "S2S qRight %g != -sigma(TL^4-TR^4) %g" % (s2s_sum["qRight"], -q_s2s))
        check(close(s2s_sum["Gmean"], G_s2s, 1e-6, 1e-6),
              "S2S Gmean %g != 2 sigma(TL^4+TR^4) %g" % (s2s_sum["Gmean"], G_s2s))
        check(close(s2s_sum["Gmin"], s2s_sum["Gmax"], 1e-9, 1e-9),
              "S2S G must be uniform across the transparent gap")
        check(abs(s2s_sum["volSrc"]) < 1e-6 * max(abs(q_s2s), 1.0),
              "S2S volumetric source must be zero")
        check(abs(s2s_sum["balance"]) / max(abs(q_s2s), 1.0) <= 1e-9,
              "S2S balance must close: %g" % s2s_sum["balance"])
        check(close(s2s_sum["viewFactor"], 1.0, 1e-12, 1e-12),
              "S2S view factor F_LR must be 1 for parallel plates")
        check(len(s2s_cells) == Nx, "S2S cell count %d != %d" % (len(s2s_cells), Nx))
        for i, c in s2s_cells.items():
            check(close(c["q"][0], q_s2s, 1e-6, 1e-6) and abs(c["q"][1]) < 1e-6
                  and abs(c["q"][2]) < 1e-6,
                  "S2S cell %d flux must be uniform sigma dT^4" % i)
            check(close(c["G"], G_s2s, 1e-6, 1e-6),
                  "S2S cell %d G must be uniform" % i)
            check(abs(c["S"]) < 1e-6, "S2S cell %d source must vanish" % i)

# ---------------------------------------------------------------- P1 trend
rp1 = run(CASE_P1)
check(rp1.returncode == 0, "P1 solver returned %d: %s" %
      (rp1.returncode, (rp1.stdout + rp1.stderr)[-400:]))
p1_sum = parse_summary(rp1.stdout)
if p1_sum and summary:
    check(p1_sum["model"] == "p1", "P1 summary model is %s" % p1_sum["model"])
    # P1 is a diffusion approximation: on this tau=1 slab its Gmean is within
    # ~15% of the exact transport value and the wall fluxes keep their signs.
    check(close(p1_sum["Gmean"], summary["Gmean"], 0.25),
          "P1/DO Gmean disagree: P1 %g vs DO %g" %
          (p1_sum["Gmean"], summary["Gmean"]))
    check(p1_sum["qLeft"] > 0 and p1_sum["qRight"] < 0,
          "P1/DO wall-flux signs disagree: %g %g" %
          (p1_sum["qLeft"], p1_sum["qRight"]))
    check(p1_sum["Gmin"] < p1_sum["Gmax"],
          "P1 G must vary across the slab")

# ---------------------------------------------------------------- refusals
refusals = [
    (open(CASE_REJECT).read(), "unknown radiation model", "bad_model"),
    (BASE.replace("absorption 1.0;", "absorption -1.0;"),
     "absorption kappa_a must be a finite positive", "neg_absorption"),
    (BASE.replace("scattering 0;", "scattering -0.5;"),
     "scattering kappa_s cannot be negative", "neg_scattering"),
    (BASE.replace("mediumTemperature 1200;", "mediumTemperature 0;"),
     "mediumTemperature must be a finite positive", "zero_medium_T"),
    (BASE.replace("wallTemperature 1500;", "wallTemperature 0;"),
     "wallTemperature must be finite and positive", "zero_wall_T"),
    (BASE.replace("quadrature S8;", "quadrature S16;"),
     "unknown quadrature", "bad_quadrature"),
    (re.sub(r"radiation\s*\{.*?^\}", "", BASE, count=1, flags=re.M | re.S),
     "needs a radiation", "no_radiation"),
]
with tempfile.TemporaryDirectory(prefix="cf_raddo_rej_") as td:
    for txt, needle, tag in refusals:
        rr = run_text(txt, tag, td)
        check(rr.returncode != 0 and needle in (rr.stdout + rr.stderr),
              "refusal '%s' not raised loudly (rc=%d)" % (needle, rr.returncode))

# ---------------------------------------------------------------- verdict
if fails:
    print("radiationDoS2s FAILED %d checks:" % len(fails))
    for f in fails[:40]:
        print("  -", f)
    sys.exit(1)

print("radiationDoS2s OK: %d assertions; analytic S%d exponential slab "
      "(worst %.2e relative), closed-form wall flux, exact discrete balance "
      "(%.1e), S2S view-factor limit, P1/DO trend, %d loud refusals" %
      (nchecks, order, worst, abs(balance_exact) / scale, len(refusals)))
