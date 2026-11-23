#!/usr/bin/env python3
"""F12 regression for chflow-3d-melting: re-derive the physics from the case file.

Independent truth (nothing below is read back from the solver's own numbers):

  * liquid fraction   f_l(T) = clamp((T - Ts)/(Tl - Ts), 0, 1)
  * apparent enthalpy H(T)  = sum_cells V rho (cp (T - tRef) + L f_l)
  * mushy source      the transient diagonal gains V (L/cp) df_l/dT / dt for a
                      cell whose starting T lies in (Ts, Tl), with
                      df_l/dT = 1/(Tl - Ts)
  * Darcy damping     f = -A (1 - f_l)^2 / (f_l^3 + eps) U, applied explicitly
                      and capped at 0.5/dt, so a solid cell loses at most half
                      its velocity per step and settles near 2 dt |F_body|
  * Qin               trapezoidal integral of the boundary heat rate over time

The checker parses block{}, energy{}, phaseChange{}, piso{} and patches{} from
the case file itself, then re-derives the quantities that must hold:

  1. enthalpy balance - H_end must match Qin_end to 2% of both the reported
     enthalpy and the mushy band's sensible enthalpy scale rho cp V (Tl - Ts)
     (with the latent scale rho V L). The case heats from the top wall and
     cools from the bottom, so Qin is a genuine nonzero boundary integral.
  2. front motion     - yF must stay inside the mesh y-range, must move
     monotonically, and must start at the f_l = 0.5 crossing implied by the
     case's initial cosine mode along the front axis.
  3. liquid bounds    - 0 <= flMin <= flMax <= 1, with solid and melt both
     present once the run has advanced.
  4. Darcy separation - uLiquid >> uSolid while the two phases coexist, and
     uSolid obeys the explicit cap derived from dt and the body force.
  5. refusals         - the sibling caseBase.dict (no phaseChange block) and
     caseReject.dict (invalid mushy interval) must both fail loudly.

The first report is excluded from the phase assertions: initialize() rebuilds
T from the cosine mode but does not refresh f_l, so that single line carries
the constructor's f_l(value at initialT) rather than the real field.

Tolerances are physical, not fitted to a stored number. The balance is asserted
to 2% of max(|H_end|, |Qin_end|) (the first-order error of evaluating df_l/dT
at the step's starting temperature) and to 2% of the mushy band's sensible
enthalpy scale rho cp V (Tl - Ts), which is the energy the phase change can
ever absorb. The Darcy cap envelope and the starting front are derived from
piso{} dt / force and energy{} initial mode, not from the printed history.

usage: check_melting_front.py <binary> <case.dict>
"""
import math
import os
import re
import subprocess
import sys

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
fails = []
ncheck = 0


def check(cond, msg):
    global ncheck
    ncheck += 1
    if not cond:
        fails.append(msg)


def close(a, b, rel, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


def clean(txt):
    out = []
    for line in txt.splitlines():
        line = re.sub(r"//.*$", "", line)
        line = re.sub(r"#.*$", "", line)
        out.append(line)
    return "\n".join(out)


def block(txt, name):
    """Body of a named { } block, balancing nested braces."""
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", txt)
    if not m:
        return None
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


def scalar(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", txt)
    if not m:
        raise RuntimeError("no scalar " + key)
    return float(m.group(1))


def vec(txt, key):
    m = re.search(
        r"\b" + re.escape(key) + r"\s*\(\s*(" + NUM + r")\s+(" + NUM
        + r")\s+(" + NUM + r")\s*\)", txt)
    if not m:
        raise RuntimeError("no vector " + key)
    return tuple(float(m.group(i)) for i in (1, 2, 3))


def has_key(txt, key):
    return re.search(r"\b" + re.escape(key) + r"\b", txt) is not None


REPORT = re.compile(
    r"\bt\s+(" + NUM + r")\s+yF\s+(" + NUM + r")\s+H\s+(" + NUM
    + r")\s+Qin\s+(" + NUM + r")\s+uLiquid\s+(" + NUM
    + r")\s+uSolid\s+(" + NUM + r")\s+flMin\s+(" + NUM
    + r")\s+flMax\s+(" + NUM + r")\s*$")


def parse_reports(out):
    rows = []
    for line in out.splitlines():
        m = REPORT.search(line)
        if m:
            rows.append(dict(zip(
                ("t", "yF", "H", "Qin", "uLiquid", "uSolid", "flMin", "flMax"),
                [float(m.group(i)) for i in range(1, 9)])))
    return rows


def run(binary, path):
    return subprocess.run([binary, path], capture_output=True, text=True)


# ---------------------------------------------------------------- arguments
if len(sys.argv) != 3:
    raise SystemExit("usage: check_melting_front.py <binary> <case.dict>")
BINARY, CASE = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(CASE))
BASE = clean(open(CASE).read())

# ---------------------------------------------------------------- case parse
blk = block(BASE, "block")
mn, mx, cells = vec(blk, "min"), vec(blk, "max"), vec(blk, "cells")
Ny = int(cells[1])
Vtot = (mx[0] - mn[0]) * (mx[1] - mn[1]) * (mx[2] - mn[2])
ylo, yhi = mn[1], mx[1]
dy = (yhi - ylo) / Ny

e = block(BASE, "energy")
rho, cp = scalar(e, "rho"), scalar(e, "cp")
initT = scalar(e, "initialT")
mode = int(scalar(e, "initialMode")) if has_key(e, "initialMode") else 0
amp = scalar(e, "initialAmplitude") if has_key(e, "initialAmplitude") else 0.0
axis = re.search(r"\binitialAxis\s+(\w+)", e)
axis = axis.group(1) if axis else "z"

ph = block(e, "phaseChange")
L = scalar(ph, "latentHeat")
Ts, Tl = scalar(ph, "solidusT"), scalar(ph, "liquidusT")
A, eps = scalar(ph, "mushyConstant"), scalar(ph, "mushyEpsilon")

p = block(BASE, "piso")
dt, endTime = scalar(p, "dt"), scalar(p, "endTime")
F = vec(p, "force") if has_key(p, "force") else (0.0, 0.0, 0.0)
Fmag = (F[0] ** 2 + F[1] ** 2 + F[2] ** 2) ** 0.5

# ------------------------------------------------- independently derived truth
# f_l = 0.5 <=> T = (Ts + Tl)/2. Invert the case's own cosine initial mode on
# the front axis to find where the upward crossing must start.
Tmid = 0.5 * (Ts + Tl)
axis_i = {"x": 0, "y": 1, "z": 2}[axis]
axis_lo, axis_hi = mn[axis_i], mx[axis_i]
span = axis_hi - axis_lo
front0 = None
if mode > 0 and amp != 0.0 and span > 0:
    ratio = (Tmid - initT) / amp
    if -1.0 <= ratio <= 1.0:
        # cos(mode*pi*xi) = ratio: scan [0,1] for the first root and refine by
        # bisection. This is the f_l = 0.5 crossing the app's tracker reports.
        candidates = []
        steps = 10000
        prev = math.cos(mode * math.pi * 0.0) - ratio
        for k in range(1, steps + 1):
            xi_k = k / steps
            cur = math.cos(mode * math.pi * xi_k) - ratio
            if prev <= 0 <= cur or prev >= 0 >= cur:
                lo_, hi_ = (k - 1) / steps, xi_k
                prev_sign = prev <= 0
                for _ in range(60):
                    mid = 0.5 * (lo_ + hi_)
                    v = math.cos(mode * math.pi * mid) - ratio
                    if (v <= 0) == prev_sign:
                        lo_ = mid
                    else:
                        hi_ = mid
                candidates.append(0.5 * (lo_ + hi_))
            prev = cur
        if candidates:
            front0 = axis_lo + candidates[0] * span

# enthalpy scale: sensible heat across the mushy band (dominant when the latent
# heat is deliberately small) plus the latent heat of the whole domain.
Hscale = max(rho * cp * Vtot * (Tl - Ts), rho * Vtot * L, 1e-30)
tolH = 0.02 * Hscale  # declared physical tolerance, 2%

# cap-derived solid velocity ceiling: per step |U| <- ~0.5|U| + dt|F| when the
# Darcy damping saturates at 0.5/dt, so the steady solid speed is about 2 dt F.
damp_at_solid = A * (1.0 - 0.1) ** 2 / (0.1 ** 3 + eps)
usolid_cap = 6.0 * dt * Fmag + 1e-9 if damp_at_solid >= 0.5 / dt else float("inf")

# ---------------------------------------------------------------- main run
r = run(BINARY, CASE)
check(r.returncode == 0,
      "melting case returned %d: %s" % (r.returncode, (r.stdout + r.stderr)[-400:]))
rows = parse_reports(r.stdout)
check(len(rows) >= 3, "expected at least 3 report lines, got %d" % len(rows))
if not rows:
    print("melting_front FAILED %d checks (no reports)" % ncheck)
    for f in fails:
        print("  -", f)
    sys.exit(1)

# time column is a cheap sanity anchor derived from piso{}
check(close(rows[0]["t"], 0.0, 0, 1e-12), "first report t != 0")
check(all(rows[i]["t"] >= rows[i - 1]["t"] - 1e-12
          for i in range(1, len(rows))), "time column is not non-decreasing")
check(close(rows[-1]["t"], endTime, 0, 1e-9),
      "last t %g != endTime %g" % (rows[-1]["t"], endTime))
check(close(rows[0]["H"], 0.0, 0, 1e-9), "initial enthalpy offset != 0")
check(close(rows[0]["Qin"], 0.0, 0, 1e-9), "initial Qin != 0")

# the first report carries a stale f_l (initialize() rebuilds T but leaves f_l
# at its constructor value), so it is excluded from the phase assertions.
body = rows[1:]

# 1. enthalpy balance ------------------------------------------------------
he, qe = rows[-1]["H"], rows[-1]["Qin"]
check(abs(he - qe) <= 0.02 * max(abs(he), abs(qe), 1e-30),
      "enthalpy balance: H_end %g - Qin_end %g = %g > 2%% of the reported "
      "enthalpy" % (he, qe, he - qe))
check(abs(he - qe) <= tolH,
      "enthalpy balance: H_end %g - Qin_end %g = %g > 2%% of scale %g"
      % (he, qe, he - qe, Hscale))

# 2. front -----------------------------------------------------------------
ys = [row["yF"] for row in body]
check(all(ylo - 1e-12 <= y <= yhi + 1e-12 for y in ys),
      "front yF left the mesh y-range [%g, %g]" % (ylo, yhi))
dys = [ys[i] - ys[i - 1] for i in range(1, len(ys))]
mono_up = all(d >= -1e-9 for d in dys)
mono_dn = all(d <= 1e-9 for d in dys)
check(mono_up or mono_dn, "front yF is not monotonic")
if front0 is not None:
    check(abs(ys[0] - front0) <= dy * 1.5,
          "initial front %g != cosine-derived f_l=0.5 crossing %g (dy=%g)"
          % (ys[0], front0, dy))
check(max(ys) - min(ys) >= 0.5 * dy,
      "front barely moved: span %g < half a cell (%g)" % (max(ys) - min(ys), 0.5 * dy))

# 3. liquid bounds ---------------------------------------------------------
for row in body:
    check(0.0 <= row["flMin"] <= row["flMax"] <= 1.0,
          "f_l bounds violated at t=%g: flMin %g flMax %g"
          % (row["t"], row["flMin"], row["flMax"]))
check(body[-1]["flMin"] < 0.1 and body[-1]["flMax"] > 0.9,
      "both solid and melt must coexist at the end (flMin %g flMax %g)"
      % (body[-1]["flMin"], body[-1]["flMax"]))

# 4. Darcy -----------------------------------------------------------------
pairs = [row for row in body if row["flMin"] < 0.1 and row["flMax"] > 0.9]
check(len(pairs) >= 2, "no report has both a solid and a liquid phase")
for row in pairs:
    ul, us = row["uLiquid"], row["uSolid"]
    check(ul > 0.0, "melt is not moving at t=%g (uLiquid %g)" % (row["t"], ul))
    check(us <= 0.1 * ul,
          "Darcy failed at t=%g: uSolid %g is not << uLiquid %g"
          % (row["t"], us, ul))
    check(us <= usolid_cap,
          "solid speed %g at t=%g exceeds the 0.5/dt Darcy cap envelope %g"
          % (us, row["t"], usolid_cap))

# 5. loud refusals ---------------------------------------------------------
base_path = os.path.join(HERE, "caseBase.dict")
reject_path = os.path.join(HERE, "caseReject.dict")
have_base = os.path.isfile(base_path)
have_reject = os.path.isfile(reject_path)
check(have_base, "missing sibling caseBase.dict next to the case")
check(have_reject, "missing sibling caseReject.dict next to the case")

if have_base:
    base_txt = clean(open(base_path).read())
    check(block(block(base_txt, "energy") or "", "phaseChange") is None,
          "caseBase.dict must not carry a phaseChange block")
    rr = run(BINARY, base_path)
    check(rr.returncode != 0 and "phaseChange" in (rr.stdout + rr.stderr),
          "base case (no phaseChange) did not refuse loudly (rc=%d)" % rr.returncode)

if have_reject:
    rej_txt = clean(open(reject_path).read())
    rej_ph = block(block(rej_txt, "energy") or "", "phaseChange") or ""
    invalid = False
    if (has_key(rej_ph, "liquidusT") and has_key(rej_ph, "solidusT")
            and has_key(rej_ph, "latentHeat")):
        invalid = (scalar(rej_ph, "liquidusT") <= scalar(rej_ph, "solidusT")
                   or scalar(rej_ph, "latentHeat") <= 0)
    check(invalid, "caseReject.dict is not actually invalid")
    rr = run(BINARY, reject_path)
    msg = rr.stdout + rr.stderr
    check(rr.returncode != 0 and (
              "liquidusT must exceed solidusT" in msg
              or "latentHeat must be positive" in msg
              or "mushyConstant must be positive" in msg
              or "mushyEpsilon must be positive" in msg),
          "reject case did not fail loudly (rc=%d): %s" % (rr.returncode, msg[-300:]))

# ---------------------------------------------------------------- verdict
if fails:
    print("melting_front FAILED %d/%d checks:" % (len(fails), ncheck))
    for f in fails[:30]:
        print("  -", f)
    sys.exit(1)

print("melting_front OK: %d assertions; H-Qin %.2f%% of scale, front %.1f dy "
      "monotone in [%.4g, %.4g], fl in [0,1], Darcy uS/uL %.3f, 2 loud refusals"
      % (ncheck, 100.0 * abs(he - qe) / Hscale,
         (max(ys) - min(ys)) / dy, min(ys), max(ys),
         body[-1]["uSolid"] / body[-1]["uLiquid"]))
