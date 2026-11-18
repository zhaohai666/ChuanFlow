#!/usr/bin/env python3
"""D8 increment 6: verify two-node QMOM constant-kernel coagulation.

Two independent anchors, neither trusting the solver:

  * moment trajectory  - replay the SAME published explicit-Euler rule of the
    Smoluchowski constant-kernel moment ODEs step by step and require every
    printed M_k to match tightly; also compare the final state to the closed
    form (M0 = M00/(1+beta M00 t/2), M1 const, M2 = M20+beta M10^2 t,
    M3 = M30 + 3 beta M10 M20 t + (3/2) beta^2 M10^2 t^2) with O(dt) slack;
  * quadrature         - the printed two nodes/weights must reproduce all four
    moments, sum_i w_i v_i^k = M_k (k=0..3), to 1e-8 relative. That identity is
    the DEFINITION of the Gaussian closure, so it is graded independently of
    how the solver reconstructed the nodes.

Plus M_1 conservation (mass), nodes positive and moving right as particles
grow, and loud refusals of non-positive kernel/moments, a moment set that is
not a positive measure, a missing block and zero steps.

usage: check_qmom.py <solver-binary> <case.dict>
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


if len(sys.argv) != 3:
    raise SystemExit("usage: check_qmom.py <solver-binary> <case.dict>")
BINARY, CASE = sys.argv[1], sys.argv[2]
BASE = open(CASE).read()
qd = block(BASE, "qmom")
beta = scalar(qd, "beta")
m0 = scalar(qd, "m0")
m1 = scalar(qd, "m1")
m2 = scalar(qd, "m2")
m3 = scalar(qd, "m3")
dt = scalar(qd, "dt")
steps = int(scalar(qd, "steps"))
every = int(scalar(qd, "reportEvery"))

r = subprocess.run([BINARY, CASE], capture_output=True, text=True)
check(r.returncode == 0, "solver returned %d: %s" % (r.returncode, (r.stdout + r.stderr)[-400:]))
out = r.stdout

records = {}
quad = {}
for line in out.splitlines():
    m = re.search(r"qmom step (\d+): t (" + NUM + r") M (" + NUM + r") ("
                  + NUM + r") (" + NUM + r") (" + NUM + r") nodes ("
                  + NUM + r") (" + NUM + r") weights (" + NUM + r") ("
                  + NUM + r") massErr (" + NUM + r")", line)
    if m:
        k = int(m.group(1))
        records[k] = {
            "t": float(m.group(2)),
            "m": [float(m.group(i)) for i in (3, 4, 5, 6)],
            "v": [float(m.group(i)) for i in (7, 8)],
            "w": [float(m.group(i)) for i in (9, 10)],
            "massErr": float(m.group(11)),
        }
    q = re.search(r"qmom moment (\d+) (\d): (" + NUM + r")", line)
    if q:
        quad[(int(q.group(1)), int(q.group(2)))] = float(q.group(3))

expected_steps = [k for k in range(0, steps + 1, every)]
if steps not in expected_steps:
    expected_steps.append(steps)
check(set(records.keys()) == set(expected_steps),
      "step records %s != expected %s" % (sorted(records), expected_steps))

# ---- replay the explicit Euler of the published moment ODEs ------------
def source(m):
    return (-0.5 * beta * m[0] ** 2,
            0.0,
            beta * m[1] ** 2,
            3.0 * beta * m[1] * m[2])


state = [m0, m1, m2, m3]
replay = {0: list(state)}
for step in range(1, steps + 1):
    dm = source(state)
    state = [state[k] + dt * dm[k] for k in range(4)]
    replay[step] = list(state)

for step, rec in records.items():
    for k in range(4):
        check(close(rec["m"][k], replay[step][k], 1e-9, 1e-9),
              "step %d M%d %g != Euler replay %g" % (step, k, rec["m"][k], replay[step][k]))

# ---- closed form at the final step (O(dt) slack for explicit Euler) ----
final = records[steps]
t = steps * dt
cf = [m0 / (1.0 + beta * m0 * t / 2.0),
      m1,
      m2 + beta * m1 ** 2 * t,
      m3 + 3.0 * beta * m1 * m2 * t + 1.5 * beta ** 2 * m1 ** 2 * t ** 2]
for k in range(4):
    check(close(final["m"][k], cf[k], 2e-2, 1e-6),
          "final M%d %g != closed form %g" % (k, final["m"][k], cf[k]))

# ---- quadrature reproduces the moments ---------------------------------
for step, rec in records.items():
    w, v = rec["w"], rec["v"]
    for k in range(4):
        rebuilt = w[0] * v[0] ** k + w[1] * v[1] ** k
        check(close(rebuilt, rec["m"][k], 1e-8, 1e-8),
              "step %d quadrature moment %d: %g != M %g" % (step, k, rebuilt, rec["m"][k]))
        qline = quad.get((step, k))
        check(qline is not None and close(qline, rebuilt, 1e-9, 1e-9),
              "step %d printed moment %d %s != rebuilt %g" % (step, k, qline, rebuilt))

# ---- mass conservation, node positivity and rightward drift ------------
for step, rec in records.items():
    check(rec["massErr"] < 1e-12, "step %d massErr %g" % (step, rec["massErr"]))
    check(rec["v"][0] > 0 and rec["v"][1] > 0, "step %d non-positive node %s" % (step, rec["v"]))
    check(rec["w"][0] >= 0 and rec["w"][1] >= 0, "step %d negative weight %s" % (step, rec["w"]))
first, last = records[0], records[steps]
check(last["v"][1] > first["v"][1], "largest node did not grow with coagulation")
check(last["m"][0] < first["m"][0], "number density did not fall with coagulation")

# ---- refusals -----------------------------------------------------------
refusals = [
    (BASE.replace("beta 0.05;", "beta -0.05;"), "beta must be a finite positive", "neg_beta"),
    (BASE.replace("m1 1.0;", "m1 -1.0;"), "M_1 must be finite and positive", "neg_m1"),
    (BASE.replace("m2 1.2;", "m2 0.5;"), "violate M_0 M_2 >= M_1^2", "bad_measure"),
    (BASE.replace("steps 20;", "steps 0;"), "steps must be at least 1", "zero_steps"),
    (re.sub(r"qmom\s*\{.*?^\}", "", BASE, count=1, flags=re.M | re.S),
     "needs a qmom", "no_block"),
]
for txt, needle, tag in refusals:
    with tempfile.TemporaryDirectory(prefix="cf_qmom_") as td:
        p = os.path.join(td, tag + ".dict")
        open(p, "w").write(txt)
        rr = subprocess.run([BINARY, p], capture_output=True, text=True)
        check(rr.returncode != 0 and needle in (rr.stdout + rr.stderr),
              "refusal '%s' not raised loudly (rc=%d)" % (needle, rr.returncode))

if fails:
    print("qmom FAILED %d checks:" % len(fails))
    for f in fails[:30]:
        print("  -", f)
    sys.exit(1)

print("qmom OK: Euler replay of constant-kernel Smoluchowski moments, closed form, "
      "two-node quadrature reproduces M0..M3, mass conserved, nodes grow, 5 loud refusals")
