#!/usr/bin/env python3
"""5.236 (H6): verify the opt-in fluid-stress load on the 6DOF body,
independently of the C++ code.

  * Closed-surface integration: on the block box (V = 0.01 m3) with the
    hydrostatic closure p = pRef + rho g.x the net pressure force is exactly
    the buoyancy -rho g V, the viscous force vanishes (no gradient declared),
    the net moment is r_c(t) x Fp about the origin, and the divergence-theorem
    volume equals the box volume.
  * The integrated load REPLACES bodyForce / drag: the translational DOF stays
    at rest (Fp.x = 0) while the rotational DOF keeps its damped oscillation
    (moment.z = 0), checked against the closed forms.
  * GCL non-regression: scalarMean stays 7 and the mass balance stays tight,
    exactly as in the F11 reference.
  * Passive twin: stripping the fluidStress block must reproduce the F11
    baseline stdout byte-for-byte.
  * Coupled arm: gravity along x drives x(t) = F/(2m) t^2 and, with the solid
    closures zeroed, theta(t) = M_z/(2I) t^2 in closed form.
  * Shear arm: a constant velocity gradient adds no net force and no net
    moment on a closed surface (Fv = 0, M unchanged).
  * Loud refusals: an open named surface, a negative density and an unknown
    patch name.

usage: check_3d_dynamicMesh6DofStress.py <solver> <case> <baseline-case> <reject-case>
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


def block_span(txt, name):
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", txt)
    if not m:
        raise RuntimeError("no " + name + " block")
    start = m.start()
    depth = 0
    for j in range(m.end() - 1, len(txt)):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
            if depth == 0:
                return start, j + 1
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
    pat = r"(\b" + re.escape(key) + r"\s+)" + NUM + r"(\s*;)"
    if not re.search(pat, txt):
        raise RuntimeError("cannot set scalar " + key)
    return re.sub(pat, lambda m: m.group(1) + str(val) + m.group(2), txt, count=1)


def set_vector(txt, key, comps):
    out = re.sub(r"(\b" + re.escape(key) + r"\s*)\([^)]*\)",
                 lambda m: m.group(1) + "(" + " ".join(str(c) for c in comps) + ")",
                 txt, count=1)
    if out == txt:
        raise RuntimeError("cannot set vector " + key)
    return out


def run(txt, tag):
    with tempfile.TemporaryDirectory(prefix="cf_stress_") as td:
        p = os.path.join(td, tag + ".dict")
        with open(p, "w") as fh:
            fh.write(txt)
        return subprocess.run([BINARY, p], capture_output=True, text=True)


def parse_dynmesh(out):
    keys = ("t", "x", "v", "theta", "omega", "scalarMean", "massBalance")
    pat = re.compile(r"DYNMESH\s+" + r"\s+".join("(" + NUM + ")" for _ in keys))
    recs = []
    for line in out.splitlines():
        m = pat.search(line)
        if m:
            recs.append({k: float(m.group(i + 1)) for i, k in enumerate(keys)})
    return recs


def parse_stress(out):
    keys = ("t", "fpx", "fpy", "fpz", "fvx", "fvy", "fvz", "mx", "my", "mz", "vol")
    pat = re.compile(r"STRESS\s+" + r"\s+".join("(" + NUM + ")" for _ in keys))
    recs = []
    for line in out.splitlines():
        m = pat.search(line)
        if m:
            recs.append({k: float(m.group(i + 1)) for i, k in enumerate(keys)})
    return recs


def body_exact(t, m, k, F, I, c, K, theta0, omega0):
    """Closed forms of the F11 checker (translation with drag k, damped
    rotational oscillator); used with F = 0 for the stress-loaded arms."""
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
        if disc > 0:
            wd = math.sqrt(disc)
            A = theta0
            B = (omega0 + alpha * theta0) / wd
            env = math.exp(-alpha * t)
            th = env * (A * math.cos(wd * t) + B * math.sin(wd * t))
            om = env * (omega0 * math.cos(wd * t) - (alpha * B + A * wd) * math.sin(wd * t))
        elif disc == 0:
            w = alpha
            A = theta0
            B = omega0 + w * theta0
            env = math.exp(-w * t)
            th = env * (A + B * t)
            om = env * (B - w * (A + B * t))
        else:
            r = math.sqrt(-disc)
            l1, l2 = -alpha + r, -alpha - r
            B = (omega0 - l1 * theta0) / (l2 - l1)
            A = theta0 - B
            th = A * math.exp(l1 * t) + B * math.exp(l2 * t)
            om = A * l1 * math.exp(l1 * t) + B * l2 * math.exp(l2 * t)
    elif c > 0:
        om = omega0 * math.exp(-c * t / I)
        th = theta0 + omega0 * I / c * (1.0 - math.exp(-c * t / I))
    else:
        om = omega0
        th = theta0 + omega0 * t
    return x, v, th, om


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


if len(sys.argv) != 5:
    raise SystemExit("usage: check_3d_dynamicMesh6DofStress.py <solver> <case> "
                     "<baseline-case> <reject-case>")
BINARY, CASE, BASELINE, REJECT = sys.argv[1:5]

# ---------------------------------------------------------------- main case
text = open(CASE).read()
blk = block(text, "block")
bmin = vectors(blk, "min")
bmax = vectors(blk, "max")
dm = block(text, "dynamicMesh")
fs = block(text, "fluidStress")
meshV = vectors(dm, "meshVelocity")
rho = scalar(fs, "density")
grav = vectors(fs, "gravity")
pRef = scalar(fs, "pressureRef")
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

Vbox = (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[2] - bmin[2])
gmag = math.sqrt(sum(g * g for g in grav))
# buoyancy of the hydrostatic closure: Fp = -rho g V
FpWant = tuple(-rho * g * Vbox for g in grav)
# the box centre follows the prescribed mesh translation
def rc(t):
    return (0.5 * (bmin[0] + bmax[0]) + meshV[0] * t,
            0.5 * (bmin[1] + bmax[1]) + meshV[1] * t,
            0.5 * (bmin[2] + bmax[2]) + meshV[2] * t)

r = subprocess.run([BINARY, CASE], capture_output=True, text=True)
check(r.returncode == 0, "solver rc=%d: %s" % (r.returncode, (r.stdout + r.stderr)[-200:]))
dyn = parse_dynmesh(r.stdout)
str_ = parse_stress(r.stdout)
check(len(dyn) == steps + 1, "DYNMESH records %d != %d" % (len(dyn), steps + 1))
check(len(str_) == steps + 1, "STRESS records %d != %d" % (len(str_), steps + 1))

for rec in str_:
    t = rec["t"]
    got = (rec["fpx"], rec["fpy"], rec["fpz"])
    for i, name in enumerate("xyz"):
        check(close(got[i], FpWant[i], 0.0, 1e-9),
              "STRESS t=%g: Fp.%s %g != buoyancy %g" % (t, name, got[i], FpWant[i]))
        check(close(rec["fv" + name], 0.0, 0.0, 1e-12),
              "STRESS t=%g: Fv.%s %g != 0 (no gradient)" % (t, name, rec["fv" + name]))
    mw = cross(rc(t), FpWant)
    for i, name in enumerate("xyz"):
        check(close(rec["m" + name], mw[i], 1e-9, 1e-9),
              "STRESS t=%g: M.%s %g != rc x Fp %g" % (t, name, rec["m" + name], mw[i]))
    check(close(rec["vol"], Vbox, 0.0, 1e-12),
          "STRESS t=%g: volume %g != box volume %g" % (t, rec["vol"], Vbox))

for rec in dyn:
    t = rec["t"]
    ex, ev, eth, eom = body_exact(t, m, k, (0.0, 0.0, 0.0), I, c, K, theta0, omega0)
    check(close(rec["x"], 0.0, 0.0, 1e-12),
          "t=%g: x %g != 0 (stress replaces bodyForce, Fp.x=0)" % (t, rec["x"]))
    check(close(rec["v"], 0.0, 0.0, 1e-12), "t=%g: v %g != 0" % (t, rec["v"]))
    check(close(rec["theta"], eth, 1e-6, 1e-7),
          "t=%g: theta %g != damped oscillator %g" % (t, rec["theta"], eth))
    check(close(rec["omega"], eom, 1e-6, 1e-7),
          "t=%g: omega %g != %g" % (t, rec["omega"], eom))
    check(close(rec["scalarMean"], phi0, 0.0, 1e-9),
          "t=%g: scalarMean %g != uniform %g (GCL regression)" % (t, rec["scalarMean"], phi0))
    check(rec["massBalance"] < 1e-9, "t=%g: GCL massBalance %g not tight" % (t, rec["massBalance"]))

# ------------------------------------------------------- passive byte twin
start, end = block_span(text, "fluidStress")
stripped = text[:start].rstrip() + "\n" + text[end:].lstrip()
check("fluidStress" not in stripped, "twin: fluidStress block not stripped")
twin = run(stripped, "twin")
base = subprocess.run([BINARY, BASELINE], capture_output=True, text=True)
check(twin.returncode == 0 and base.returncode == 0, "twin run failed")
check(twin.stdout == base.stdout,
      "twin: stdout differs from the F11 baseline (not byte-identical)")
check("STRESS" not in twin.stdout, "twin: undeclared stress must print no STRESS line")

# ---------------------------------------------------------------- coupled arm
coup = set_vector(text, "gravity", (-9.81, 0, 0))
coup = set_scalar(coup, "rotDrag", 0)
coup = set_scalar(coup, "rotStiffness", 0)
coup = set_scalar(coup, "theta0", 0)
coup = set_scalar(coup, "omega0", 0)
cr = run(coup, "coupled")
check(cr.returncode == 0, "coupled rc=%d: %s" % (cr.returncode, (cr.stdout + cr.stderr)[-200:]))
cdyn, cstr = parse_dynmesh(cr.stdout), parse_stress(cr.stdout)
check(len(cdyn) == steps + 1 and len(cstr) == steps + 1, "coupled: record counts")
Fx = -rho * (-9.81) * Vbox  # = 98.1
for rec in cstr:
    check(close(rec["fpx"], Fx, 1e-9, 1e-9), "coupled: Fp.x %g != %g" % (rec["fpx"], Fx))
    check(close(rec["fpy"], 0.0, 0.0, 1e-9) and close(rec["fpz"], 0.0, 0.0, 1e-9),
          "coupled: Fp.y/z must vanish")
    mw = cross(rc(rec["t"]), (Fx, 0.0, 0.0))
    check(close(rec["mz"], mw[2], 1e-9, 1e-9), "coupled: M.z %g != %g" % (rec["mz"], mw[2]))
Mz = -0.05 * Fx  # = -4.905, constant
for rec in cdyn:
    t = rec["t"]
    check(close(rec["x"], Fx / (2.0 * m) * t * t, 1e-9, 1e-9),
          "coupled t=%g: x %g != F/(2m) t^2" % (t, rec["x"]))
    check(close(rec["v"], Fx / m * t, 1e-9, 1e-9), "coupled t=%g: v %g" % (t, rec["v"]))
    check(close(rec["theta"], Mz / (2.0 * I) * t * t, 1e-9, 1e-9),
          "coupled t=%g: theta %g != Mz/(2I) t^2" % (t, rec["theta"]))
    check(close(rec["omega"], Mz / I * t, 1e-9, 1e-9), "coupled t=%g: omega %g" % (t, rec["omega"]))

# ------------------------------------------------------------------ shear arm
shear = text.replace("pressureRef 0;",
                     "pressureRef 0;\n    velocityGradient (0 1 0 0 0 0 0 0 0);", 1)
check(shear != text, "shear: could not inject velocityGradient")
sr = run(shear, "shear")
check(sr.returncode == 0, "shear rc=%d: %s" % (sr.returncode, (sr.stdout + sr.stderr)[-200:]))
sdyn, sstr = parse_dynmesh(sr.stdout), parse_stress(sr.stdout)
check(len(sdyn) == steps + 1 and len(sstr) == steps + 1, "shear: record counts")
for rec in sstr:
    for name in "xyz":
        check(close(rec["fv" + name], 0.0, 0.0, 1e-12),
              "shear t=%g: Fv.%s %g != 0 (constant gradient on a closed surface)"
              % (rec["t"], name, rec["fv" + name]))
    mw = cross(rc(rec["t"]), FpWant)
    for i, name in enumerate("xyz"):
        check(close(rec["m" + name], mw[i], 1e-9, 1e-9),
              "shear t=%g: M.%s %g != pressure-only moment %g"
              % (rec["t"], name, rec["m" + name], mw[i]))
for a, b in zip(sdyn, dyn):
    for key in ("x", "v", "theta", "omega"):
        check(close(a[key], b[key], 1e-9, 1e-9),
              "shear t=%g: body %s changed by the shear arm" % (a["t"], key))

# ----------------------------------------------------------------- refusals
refusals = [
    (open(REJECT).read(), "not closed", "open_surface"),
    (set_scalar(text, "density", -1), "density must be finite and positive", "neg_density"),
    (text.replace("pressureRef 0;", "pressureRef 0;\n    surface (nosuchpatch);", 1),
     "no patch named", "unknown_patch"),
]
for txt, needle, tag in refusals:
    rr = run(txt, tag)
    check(rr.returncode != 0 and needle in (rr.stdout + rr.stderr),
          "refusal '%s' not raised loudly (rc=%d): %s"
          % (needle, rr.returncode, (rr.stdout + rr.stderr)[-160:]))

if fails:
    print("FAIL: dynamicMesh6DofStress (%d/%d assertions failed)" % (len(fails), count))
    for f in fails[:40]:
        print("  -", f)
    sys.exit(1)

print("PASS: dynamicMesh6DofStress — %d assertions (buoyancy -rho g V, moment "
      "rc x Fp, zero viscous load, stress-driven body ODE, GCL non-regression, "
      "byte-identical passive twin, coupled + shear arms, 3 loud rejections)" % count)
sys.exit(0)
