#!/usr/bin/env python3
"""D8 increment 8: verify moving-mesh geometry and the Geometric Conservation
Law with two motions whose discrete answers are known:

  translate - every vertex moves at a constant velocity, so each cell's face
      area vectors are unchanged and sum(U_m.Sf) over a closed cell is exactly
      zero: dV/dt = 0 must hold to ROUNDOFF (>1e-12 fails). The mesh must
      translate rigidly (xMax = xMax0 + v*t) with constant total volume.

  expand    - U_m = rate*x makes every vertex position multiply by (1+rate*dt)
      each step, so after N steps xMax = xMax0 (1+rate dt)^N and the volume
      (all three coordinates stretch) is V0 (1+rate dt)^(3N). The forward-Euler
      GCL residual dV/dt - sum(U_m.Sf) is O(dt); the checker runs three dt with
      the same end time and requires the final max|GCL| to halve as dt halves.

usage: check_ale.py <solver-binary> <translate-case> <expand-case>
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


def run(txt, tag):
    with tempfile.TemporaryDirectory(prefix="cf_ale_") as td:
        p = os.path.join(td, tag + ".dict")
        open(p, "w").write(txt)
        return subprocess.run([BINARY, p], capture_output=True, text=True)


def steps_of(out):
    recs = []
    for line in out.splitlines():
        m = re.search(r"ale step (\d+): t (" + NUM + r") totalVolume (" + NUM
                      + r") maxGcl (" + NUM + r") sumGcl (" + NUM + r") xMin ("
                      + NUM + r") xMax (" + NUM + r")", line)
        if m:
            recs.append({
                "step": int(m.group(1)), "t": float(m.group(2)),
                "vol": float(m.group(3)), "maxGcl": float(m.group(4)),
                "sumGcl": float(m.group(5)), "xMin": float(m.group(6)),
                "xMax": float(m.group(7)),
            })
    return recs


if len(sys.argv) != 4:
    raise SystemExit("usage: check_ale.py <solver> <translate-case> <expand-case>")
BINARY, TCASE, ECASE = sys.argv[1], sys.argv[2], sys.argv[3]

# ---------------- translation -------------------------------------------
tb = open(TCASE).read()
tad = block(tb, "ale")
tv = vectors(tad, "translation")
tdt = scalar(tad, "dt")
tsteps = int(scalar(tad, "steps"))
tr = subprocess.run([BINARY, TCASE], capture_output=True, text=True)
check(tr.returncode == 0, "translate solver rc=%d" % tr.returncode)
mi = re.search(r"ale initial: totalVolume (" + NUM + r")", tr.stdout)
V0 = float(mi.group(1)) if mi else 0
trecs = steps_of(tr.stdout)
check(len(trecs) == tsteps, "translate records %d != %d" % (len(trecs), tsteps))
mmax = 0
x0 = None
for rec in trecs:
    mmax = max(mmax, rec["maxGcl"])
    if x0 is None:
        x0 = rec["xMax"] - tv[0] * rec["t"]  # xMax0
    check(rec["maxGcl"] < 1e-12,
          "translate step %d: maxGcl %g not roundoff" % (rec["step"], rec["maxGcl"]))
    check(close(rec["vol"], V0, 1e-12, 1e-15),
          "translate step %d: volume %g changed from %g" % (rec["step"], rec["vol"], V0))
    check(close(rec["xMax"], x0 + tv[0] * rec["t"], 1e-9, 1e-9),
          "translate step %d: xMax %g != rigid %g" % (rec["step"], rec["xMax"], x0 + tv[0] * rec["t"]))
    check(close(rec["sumGcl"], 0.0, 0, 1e-12),
          "translate step %d: sumGcl %g != 0" % (rec["step"], rec["sumGcl"]))

# ---------------- expansion ---------------------------------------------
eb = open(ECASE).read()
ead = block(eb, "ale")
rate = scalar(ead, "expansionRate")
edt = scalar(ead, "dt")
esteps = int(scalar(ead, "steps"))
eblk = block(eb, "block")
emin, emax = vectors(eblk, "min"), vectors(eblk, "max")
Lx0 = emax[0] - emin[0]
Ly0 = emax[1] - emin[1]
Lz0 = emax[2] - emin[2]
V0e = Lx0 * Ly0 * Lz0 * 1.0  # cells(20 1 1) -> volume = bbox

er = subprocess.run([BINARY, ECASE], capture_output=True, text=True)
check(er.returncode == 0, "expand solver rc=%d" % er.returncode)
erecs = steps_of(er.stdout)
check(len(erecs) == esteps, "expand records %d != %d" % (len(erecs), esteps))
for rec in erecs:
    factor = (1.0 + rate * edt) ** rec["step"]
    check(close(rec["vol"], V0e * factor ** 3, 1e-9, 1e-12),
          "expand step %d: volume %g != %g" % (rec["step"], rec["vol"], V0e * factor ** 3))
    check(close(rec["xMax"], emax[0] * factor, 1e-9, 1e-9),
          "expand step %d: xMax %g != %g" % (rec["step"], rec["xMax"], emax[0] * factor))
    check(rec["maxGcl"] < 1e-3, "expand step %d: maxGcl %g too large" % (rec["step"], rec["maxGcl"]))

# ---- GCL error is first order in dt (same end time) ---------------------
def final_gcl(dt, steps):
    txt = eb.replace("dt %s;" % edt, "dt %s;" % dt).replace("steps %d;" % esteps,
                                                            "steps %d;" % steps)
    rr = run(txt, "dt%s" % dt)
    recs = steps_of(rr.stdout)
    return recs[-1]["maxGcl"] if recs else float("inf")

g1 = final_gcl("0.01", 10)   # t = 0.10
g2 = final_gcl("0.005", 20)  # t = 0.10
g4 = final_gcl("0.0025", 40) # t = 0.10
ratio1 = g1 / g2
ratio2 = g2 / g4
check(1.7 < ratio1 < 2.3, "GCL error not ~halved from dt 0.01->0.005 (ratio %g)" % ratio1)
check(1.7 < ratio2 < 2.3, "GCL error not ~halved from dt 0.005->0.0025 (ratio %g)" % ratio2)

# ---- refusals -----------------------------------------------------------
refusals = [
    (tb.replace("dt 0.1;", "dt -0.1;"), "dt must be finite and positive", "neg_dt"),
    (tb.replace("steps 10;", "steps 0;"), "steps must be at least 1", "zero_steps"),
    (tb.replace("motion translate;", "motion spin;"),
     "motion must be 'translate' or 'expand'", "bad_motion"),
    (eb.replace("expansionRate 1.0;", "expansionRate 0;"),
     "expansionRate must be finite and positive", "zero_rate"),
    (re.sub(r"ale\s*\{.*?^\}", "", tb, count=1, flags=re.M | re.S),
     "needs an ale", "no_block"),
]
for txt, needle, tag in refusals:
    rr = run(txt, tag)
    check(rr.returncode != 0 and needle in (rr.stdout + rr.stderr),
          "refusal '%s' not raised loudly (rc=%d): %s"
          % (needle, rr.returncode, (rr.stdout + rr.stderr)[-160:]))

if fails:
    print("ale FAILED %d checks:" % len(fails))
    for f in fails[:30]:
        print("  -", f)
    sys.exit(1)

print("ale OK: translation GCL at roundoff (max %.1e), rigid motion, expansion "
      "volume law, GCL error first order (ratios %.2f / %.2f), 5 loud refusals"
      % (mmax, ratio1, ratio2))
