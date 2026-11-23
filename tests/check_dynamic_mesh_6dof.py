#!/usr/bin/env python3
"""F11: verify GCL-corrected relative-flux transport on a translating mesh and
the rigid 6DOF (translation + rotation) update, independently of the C++ code.

  * Uniform-scalar invariance: a spatially constant phi is transported by the
    relative face flux phi_rel = (U - U_mesh).Sf with the GCL correction, so a
    rigid translation leaves scalarMean exactly phi0 at every step.
  * Mass balance: the per-cell GCL residual dV/dt - sum_f U_mesh.Sf must be at
    roundoff for a rigid translation (no geometric truncation).
  * Rigid body: the translational ODE m v' = F - k v is integrated in closed
    form (v -> F/k), and the rotational ODE I theta'' + c theta' + K theta = 0
    is integrated as a damped harmonic oscillator. The app's RK4 solution is
    compared against both closed forms.
  * Static twin: meshVelocity = 0 and bodyForce = 0 must reduce to the
    fixed-mesh answer (phi unchanged, GCL zero, body at rest).
  * Loud refusals of negative mass, drag and inertia.

usage: check_dynamic_mesh_6dof.py <solver> <case> <static-case> <reject-case>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
fails = []
count = 0


def check(cond, msg):
    global count
    count += 1
    if not cond:
        fails.append(msg)


def close(a, b, rel, abs_=0.0):
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


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


def scalar(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(" + NUM + r")\s*;", txt)
    if not m:
        raise RuntimeError("no scalar " + key)
    return float(m.group(1))


def vectors(txt, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(\s*(" + NUM + r")\s+("
                  + NUM + r")\s+(" + NUM + r")\s*\)", txt)
    if not m:
        raise RuntimeError("no vector " + key)
    return tuple(float(m.group(i)) for i in (1, 2, 3))


def set_scalar(txt, key, val):
    return re.sub(r"(\b" + re.escape(key) + r"\s+)" + NUM + r"(\s*;)",
                  lambda m: m.group(1) + str(val) + m.group(2), txt, count=1)


def run(txt, tag):
    with tempfile.TemporaryDirectory(prefix="cf_dynmesh_") as td:
        p = os.path.join(td, tag + ".dict")
        with open(p, "w") as fh:
            fh.write(txt)
        return subprocess.run([BINARY, p], capture_output=True, text=True)


def parse(out):
    keys = ("t", "x", "v", "theta", "omega", "scalarMean", "massBalance")
    pat = re.compile(r"DYNMESH\s+" + r"\s+".join("(" + NUM + ")" for _ in keys))
    recs = []
    for line in out.splitlines():
        m = pat.search(line)
        if m:
            recs.append({k: float(m.group(i + 1)) for i, k in enumerate(keys)})
    return recs


def body_exact(t, m, k, F, I, c, K, theta0, omega0):
    """Closed forms, derived independently of the C++ integrator."""
    if k > 0:
        tau = m / k
        vt = F[0] / k
        decay = math.exp(-t / tau)
        v = vt * (1.0 - decay)
        x = vt * (t - tau * (1.0 - decay))
    else:
        v = F[0] * t / m
        x = F[0] * t * t / (2.0 * m)

    if K > 0:
        alpha = c / (2.0 * I)
        disc = K / I - alpha * alpha
        if disc > 0:              # underdamped (the case's regime)
            wd = math.sqrt(disc)
            A = theta0
            B = (omega0 + alpha * theta0) / wd
            env = math.exp(-alpha * t)
            th = env * (A * math.cos(wd * t) + B * math.sin(wd * t))
            om = env * (omega0 * math.cos(wd * t) - (alpha * B + A * wd) * math.sin(wd * t))
        elif disc == 0:           # critically damped
            w = alpha
            A = theta0
            B = omega0 + w * theta0
            env = math.exp(-w * t)
            th = env * (A + B * t)
            om = env * (B - w * (A + B * t))
        else:                     # overdamped
            r = math.sqrt(-disc)
            l1, l2 = -alpha + r, -alpha - r
            B = (omega0 - l1 * theta0) / (l2 - l1)
            A = theta0 - B
            th = A * math.exp(l1 * t) + B * math.exp(l2 * t)
            om = A * l1 * math.exp(l1 * t) + B * l2 * math.exp(l2 * t)
    elif c > 0:                   # free rotation with drag
        om = omega0 * math.exp(-c * t / I)
        th = theta0 + omega0 * I / c * (1.0 - math.exp(-c * t / I))
    else:
        om = omega0
        th = theta0 + omega0 * t
    return x, v, th, om


if len(sys.argv) != 5:
    raise SystemExit("usage: check_dynamic_mesh_6dof.py <solver> <case> <static-case> <reject-case>")
BINARY, CASE, STATIC, REJECT = sys.argv[1:5]

# ---------------------------------------------------------------- main case
text = open(CASE).read()
dm = block(text, "dynamicMesh")
meshV = vectors(dm, "meshVelocity")
F = vectors(dm, "bodyForce")
m = scalar(dm, "mass")
k = scalar(dm, "drag")
I = scalar(dm, "inertia")
c = scalar(dm, "rotDrag")
K = scalar(dm, "rotStiffness")
theta0 = scalar(dm, "theta0")
omega0 = scalar(dm, "omega0")
phi0 = scalar(dm, "initialValue")
dt = scalar(dm, "dt")
steps = int(scalar(dm, "steps"))
reportEvery = int(scalar(dm, "reportEvery"))

r = subprocess.run([BINARY, CASE], capture_output=True, text=True)
check(r.returncode == 0, "solver rc=%d: %s" % (r.returncode, (r.stdout + r.stderr)[-200:]))
recs = parse(r.stdout)
expected = steps + 1 if reportEvery == 1 else 1 + (steps + reportEvery - 1) // reportEvery
check(len(recs) == expected, "DYNMESH records %d != %d" % (len(recs), expected))

maxGcl = 0.0
for rec in recs:
    t = rec["t"]
    maxGcl = max(maxGcl, rec["massBalance"])
    # uniform scalar invariance under rigid translation
    check(close(rec["scalarMean"], phi0, 0.0, 1e-9),
          "t=%g: scalarMean %g != uniform %g" % (t, rec["scalarMean"], phi0))
    check(rec["massBalance"] < 1e-9,
          "t=%g: GCL massBalance %g not tight" % (t, rec["massBalance"]))
    # rigid-body closed forms
    ex, ev, eth, eom = body_exact(t, m, k, F, I, c, K, theta0, omega0)
    check(close(rec["x"], ex, 1e-6, 1e-7),
          "t=%g: x %g != exact %g" % (t, rec["x"], ex))
    check(close(rec["v"], ev, 1e-6, 1e-7),
          "t=%g: v %g != exact %g" % (t, rec["v"], ev))
    check(close(rec["theta"], eth, 1e-6, 1e-7),
          "t=%g: theta %g != exact %g" % (t, rec["theta"], eth))
    check(close(rec["omega"], eom, 1e-6, 1e-7),
          "t=%g: omega %g != exact %g" % (t, rec["omega"], eom))

check(meshV == (0.5, 0.0, 0.0), "main case meshVelocity fixture changed")

# ---------------------------------------------------------------- static twin
sr = subprocess.run([BINARY, STATIC], capture_output=True, text=True)
check(sr.returncode == 0, "static solver rc=%d" % sr.returncode)
srecs = parse(sr.stdout)
stext = open(STATIC).read()
sdm = block(stext, "dynamicMesh")
ssteps = int(scalar(sdm, "steps"))
sphi0 = scalar(sdm, "initialValue")
check(len(srecs) == ssteps + 1,
      "static DYNMESH records %d != %d" % (len(srecs), ssteps + 1))
for rec in srecs:
    t = rec["t"]
    check(close(rec["x"], 0.0, 0.0, 1e-12), "static t=%g: x %g != 0" % (t, rec["x"]))
    check(close(rec["v"], 0.0, 0.0, 1e-12), "static t=%g: v %g != 0" % (t, rec["v"]))
    check(close(rec["theta"], 0.0, 0.0, 1e-12), "static t=%g: theta %g != 0" % (t, rec["theta"]))
    check(close(rec["omega"], 0.0, 0.0, 1e-12), "static t=%g: omega %g != 0" % (t, rec["omega"]))
    check(close(rec["scalarMean"], sphi0, 0.0, 1e-9),
          "static t=%g: scalarMean %g != %g" % (t, rec["scalarMean"], sphi0))
    check(rec["massBalance"] < 1e-12, "static t=%g: massBalance %g" % (t, rec["massBalance"]))

# ---------------------------------------------------------------- refusals
refusals = [
    (open(REJECT).read(), "mass must be finite and positive", "neg_mass"),
    (set_scalar(text, "drag", -1), "drag must be finite and non-negative", "neg_drag"),
    (set_scalar(text, "inertia", -0.5), "inertia must be finite and positive", "neg_inertia"),
]
for txt, needle, tag in refusals:
    rr = run(txt, tag)
    check(rr.returncode != 0 and needle in (rr.stdout + rr.stderr),
          "refusal '%s' not raised loudly (rc=%d): %s"
          % (needle, rr.returncode, (rr.stdout + rr.stderr)[-160:]))

if fails:
    print("FAIL: dynamicMesh6Dof (%d/%d assertions failed)" % (len(fails), count))
    for f in fails[:40]:
        print("  -", f)
    sys.exit(1)

print("PASS: dynamicMesh6Dof — %d assertions (uniform-scalar invariance, "
      "GCL mass balance max %.1e, closed-form translation/rotation, static "
      "twin, 3 loud rejections)" % (count, maxGcl))
