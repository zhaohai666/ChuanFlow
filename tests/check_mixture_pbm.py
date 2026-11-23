#!/usr/bin/env python3
"""Batch D item F7: verify chflow-3d-mixture-pbm independently.

The solver's printed trajectory is never used as truth. The checker parses the
case text itself, re-derives the documented monodisperse moment ODEs and the
drift-flux slip relation, then grades the solver against them:

  1. replay the SAME explicit-Euler rule of the monodisperse moment ODEs
       dM0/dt = J + (nFrag-1) b M0
       dM1/dt = J v                      (breakage conserves volume)
       dM2/dt = J v^2 - (1-1/nFrag) b v^2 M0
       dM3/dt = J v^3 - (1-1/nFrag^2) b v^3 M0
     with v = pi d0^3/6, from M_0 = (alpha0/v, alpha0, alpha0 v, alpha0 v^2),
     and require every reported M0..M3 and alpha/uG/uL to match tightly;
  2. closed forms: M1(t) = alpha0 + J v t is exact (linear), M0(t) is the
     exponential (nFrag-1)b growth with nucleation, uG-uL = slip identically;
  3. physicality: every moment positive and finite, alpha inside (0,1),
     growth in M0/M1 with positive rates;
  4. slip changes phase distribution: re-run with slip=0 and require its holdup
     to equal M1 (no redistribution) and exceed the slip>0 main holdup;
  5. no-growth twin (caseBase.dict, J=breakage=0) keeps every moment constant;
  6. loud refusals: invalid rho (caseReject.dict), negative J, negative
     breakage, alpha0 too large, and a removed mixturePbm{} block.

usage: check_mixture_pbm.py <solver-binary> <case.dict>
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

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


def parse_case(path):
    d = block(open(path).read(), "mixturePbm")
    return {
        "alpha0": scalar(d, "alpha0"), "d0": scalar(d, "d0"),
        "rhoL": scalar(d, "rhoL"), "rhoG": scalar(d, "rhoG"),
        "slip": scalar(d, "slip"), "j": scalar(d, "j"),
        "J": scalar(d, "J"), "breakage": scalar(d, "breakage"),
        "nFrag": scalar(d, "nFrag"), "dt": scalar(d, "dt"),
        "steps": int(scalar(d, "steps")),
        "reportEvery": int(scalar(d, "reportEvery")),
    }


def derive(cfg):
    """Independently re-derive every closure constant from the case text."""
    d0 = cfg["d0"]
    v = math.pi * d0 * d0 * d0 / 6.0
    k = cfg["j"] / (cfg["j"] + cfg["slip"])
    m0 = cfg["alpha0"] / v
    m = [m0, cfg["alpha0"], m0 * v * v, m0 * v * v * v]
    return {"v": v, "k": k, "m0": m}


def source(cfg, d, m):
    v = d["v"]
    n = cfg["nFrag"]
    f2 = 1.0 - 1.0 / n
    f3 = 1.0 - 1.0 / (n * n)
    b = cfg["breakage"]
    return [
        cfg["J"] + (n - 1.0) * b * m[0],
        cfg["J"] * v,
        cfg["J"] * v * v - f2 * b * v * v * m[0],
        cfg["J"] * v * v * v - f3 * b * v * v * v * m[0],
    ]


def replay(cfg, d):
    """Explicit-Euler replay from the case constants, independent of the app."""
    m = list(d["m0"])
    out = {0: list(m)}
    for step in range(1, cfg["steps"] + 1):
        dm = source(cfg, d, m)
        m = [m[i] + cfg["dt"] * dm[i] for i in range(4)]
        out[step] = list(m)
    return out


def parse_steps(output):
    rec = {}
    pat = (r"MIXTUREPBM step (\d+) t (" + NUM + r") M0 (" + NUM + r") M1 ("
           + NUM + r") M2 (" + NUM + r") M3 (" + NUM + r") alpha (" + NUM
           + r") uG (" + NUM + r") uL (" + NUM + r")")
    for line in output.splitlines():
        m = re.search(pat, line)
        if m:
            k = int(m.group(1))
            rec[k] = {
                "t": float(m.group(2)),
                "M": [float(m.group(i)) for i in (3, 4, 5, 6)],
                "alpha": float(m.group(7)),
                "uG": float(m.group(8)), "uL": float(m.group(9)),
            }
    return rec


def write_case(txt, td, tag):
    p = os.path.join(td, tag + ".dict")
    with open(p, "w") as f:
        f.write(txt)
    return p


def run(binary, path, td, tag):
    with open(path) as f:
        txt = f.read()
    p = write_case(txt, td, tag)
    return subprocess.run([binary, p], capture_output=True, text=True)


if len(sys.argv) != 3:
    raise SystemExit("usage: check_mixture_pbm.py <solver-binary> <case.dict>")
BINARY, CASE = sys.argv[1], sys.argv[2]
here = os.path.dirname(CASE)
BASE = open(CASE).read()
cfg = parse_case(CASE)
d = derive(cfg)

td = tempfile.mkdtemp(prefix="cf_mixture_pbm_")
main = run(BINARY, CASE, td, "main")
twin = run(BINARY, os.path.join(here, "caseBase.dict"), td, "base")
reject = run(BINARY, os.path.join(here, "caseReject.dict"), td, "reject")

check(main.returncode == 0, "main case exited %d: %s"
      % (main.returncode, (main.stdout + main.stderr)[-300:]))
check(twin.returncode == 0, "no-growth twin exited %d" % twin.returncode)
check(reject.returncode != 0, "invalid rho must be refused loudly")
check("rhoG must be finite and positive" in (reject.stdout + reject.stderr),
      "invalid-rho refusal must name rhoG")

steps = parse_steps(main.stdout)
expected = [k for k in range(0, cfg["steps"] + 1, cfg["reportEvery"])]
if cfg["steps"] not in expected:
    expected.append(cfg["steps"])
check(set(steps) == set(expected), "step records %s != expected %s"
      % (sorted(steps), expected))

# ---- 1. explicit-Euler replay is the primary independent truth ----------
hist = replay(cfg, d)
for step, rec in sorted(steps.items()):
    check(step in hist, "step %d not replayable" % step)
    if step not in hist:
        continue
    check(close(rec["t"], step * cfg["dt"], 1e-9, 1e-12),
          "step %d t %g != %g" % (step, rec["t"], step * cfg["dt"]))
    for i in range(4):
        check(close(rec["M"][i], hist[step][i], 1e-8, 1e-30),
              "step %d M%d %g != Euler replay %g"
              % (step, i, rec["M"][i], hist[step][i]))
    alpha = hist[step][1] * d["k"]
    check(close(rec["alpha"], alpha, 1e-8, 1e-30),
          "step %d alpha %g != M1*k %g" % (step, rec["alpha"], alpha))
    uG = cfg["j"] + (1.0 - alpha) * cfg["slip"]
    uL = cfg["j"] - alpha * cfg["slip"]
    check(close(rec["uG"], uG, 1e-8, 1e-12),
          "step %d uG %g != %g" % (step, rec["uG"], uG))
    check(close(rec["uL"], uL, 1e-8, 1e-12),
          "step %d uL %g != %g" % (step, rec["uL"], uL))

# ---- 2. closed forms (independent of the step replay) -------------------
t = cfg["steps"] * cfg["dt"]
last = steps[cfg["steps"]]
m1_cf = cfg["alpha0"] + cfg["J"] * d["v"] * t
check(close(last["M"][1], m1_cf, 1e-8, 1e-15),
      "final M1 %g != linear closed form %g" % (last["M"][1], m1_cf))
A = (cfg["nFrag"] - 1.0) * cfg["breakage"]
if A > 0:
    m0_cf = (d["m0"][0] + cfg["J"] / A) * math.exp(A * t) - cfg["J"] / A
    check(close(last["M"][0], m0_cf, 2e-2, 1e-12),
          "final M0 %g != exponential closed form %g" % (last["M"][0], m0_cf))
check(close(last["uG"] - last["uL"], cfg["slip"], 1e-8, 1e-12),
      "uG-uL %g != slip %g" % (last["uG"] - last["uL"], cfg["slip"]))

# ---- 3. physicality and growth ------------------------------------------
for step, rec in sorted(steps.items()):
    for i in range(4):
        check(math.isfinite(rec["M"][i]) and rec["M"][i] > 0,
              "step %d M%d non-physical: %g" % (step, i, rec["M"][i]))
    check(0.0 < rec["alpha"] < 1.0,
          "step %d alpha outside (0,1): %g" % (step, rec["alpha"]))
check(last["M"][0] > steps[0]["M"][0], "M0 must grow with nucleation+breakage")
check(last["M"][1] > steps[0]["M"][1], "M1 must grow with nucleation")
check(last["M"][2] > 0 and last["M"][3] > 0, "M2/M3 must stay positive")

# ---- 4. slip changes phase distribution ---------------------------------
p0_path = write_case(BASE.replace("slip 0.2;", "slip 0.0;"), td, "slip0")
p0 = subprocess.run([BINARY, p0_path], capture_output=True, text=True)
check(p0.returncode == 0, "slip=0 twin exited %d" % p0.returncode)
s0 = parse_steps(p0.stdout)
if 0 in s0 and 0 in steps:
    check(close(s0[0]["alpha"], cfg["alpha0"], 1e-8, 1e-12),
          "slip=0 holdup at t=0 must equal alpha0: %g" % s0[0]["alpha"])
    check(s0[0]["alpha"] > steps[0]["alpha"],
          "removing slip must raise holdup: %g vs %g"
          % (s0[0]["alpha"], steps[0]["alpha"]))
    check(close(s0[0]["uG"] - s0[0]["uL"], 0.0, 1e-8, 1e-12),
          "slip=0 twin must have uG==uL")
check(close(steps[0]["alpha"], cfg["alpha0"] * d["k"], 1e-8, 1e-12),
      "main holdup != alpha0*j/(j+slip): slip has no effect")

# ---- 5. no-growth twin: rates zero reproduce the constant analytic state --
bt = parse_steps(twin.stdout)
bd = derive(parse_case(os.path.join(here, "caseBase.dict")))
for step, rec in sorted(bt.items()):
    for i in range(4):
        check(close(rec["M"][i], bd["m0"][i], 1e-8, 1e-30),
              "no-growth step %d M%d %g != constant %g"
              % (step, i, rec["M"][i], bd["m0"][i]))
    check(close(rec["alpha"], bd["m0"][1] * bd["k"], 1e-8, 1e-30),
          "no-growth step %d alpha drifted" % step)

# ---- 6. loud refusals ----------------------------------------------------
refusals = [
    (BASE.replace("J 1000000;", "J -1000000;"),
     "J (nucleation rate) must be finite and non-negative", "neg_J"),
    (BASE.replace("breakage 0.5;", "breakage -0.5;"),
     "breakage frequency must be finite and non-negative", "neg_breakage"),
    (BASE.replace("alpha0 0.05;", "alpha0 1.5;"),
     "alpha0 too large", "large_alpha"),
    (BASE.replace("rhoL 1000;", "rhoL -1000;"),
     "rhoL must be finite and positive", "neg_rhoL"),
    (re.sub(r"mixturePbm\s*\{.*?^\}", "", BASE, count=1, flags=re.M | re.S),
     "needs a mixturePbm", "no_block"),
]
for txt, needle, tag in refusals:
    p = write_case(txt, td, tag)
    r = subprocess.run([BINARY, p], capture_output=True, text=True)
    check(r.returncode != 0 and needle in (r.stdout + r.stderr),
          "refusal '%s' not raised loudly (rc=%d)" % (needle, r.returncode))

shutil.rmtree(td, ignore_errors=True)
if fails:
    print("mixture_pbm FAILED %d checks:" % len(fails))
    for f in fails[:30]:
        print("  -", f)
    sys.exit(1)
print("mixture_pbm OK: %d assertions; Euler replay of monodisperse moments, "
      "closed forms, physical moments, slip redistributes holdup, no-growth "
      "twin constant, 5 loud refusals" % ncheck)
