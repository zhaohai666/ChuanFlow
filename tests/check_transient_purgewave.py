"""O3 (5.162): verify the distributed riser purge (blowdown) wave.

The physics is re-derived INDEPENDENTLY in Python - the reference model
here is the documented closure and the documented wave kinematics, not a
copy of the C++:

  * the wallis liquid flux F(a) = us a^2 + (vM - us) a with the mixture
    velocity vM = vsl + vsg imposed by the rate inlet;
  * the refill root alpha_root (F(a) = vsl) and the fill-shock speed
    s = (vsl - F(alpha_res)) / (alpha_root - alpha_res), which set the
    accumulation time of the cycle;
  * the cycle period T = L/uWave + L/s (traverse + refill), which must
    land in the severe-slugging frequency band;
  * the trigger inequality cap < factor rhoL g H on the reported values;
  * the implicit-Euler cap relaxation cap_k = (cap_{k-1}/dt +
    pbase_k/tau)/(1/dt + 1/tau) from the reported pbase;
  * the closed-mode front speed
    u = sqrt(2 max(factor rhoL g H_trig - cap, 0) / rhoL)
    and the step-by-step integration x_k = x_{k-1} + u_k dt;
  * the blowdown accounting: every cell the front passed sits at the
    residual holdup exactly, the purged volume of a cycle matches the
    head drop up to the bounded transport inflow, and the hydrostatic
    dp across the pipe collapses during the purge and recovers after it.

The undeclared path is witnessed twice: caseDisabled.dict (enable false)
must be byte-identical to the twin caseBase.dict without the block, and
a declared-but-never-triggering variant (triggerFactor 1e-6) must match
caseBase.dict byte-for-byte once its purge state lines are stripped.

Usage: check_transient_purgewave.py <chflow-1d-twofluid> <casedir>
"""

import math
import os
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)
EXE, CASEDIR = sys.argv[1:3]

fails = []


def check(name, got, want, rtol, atol=0.0):
    tol = atol + rtol * abs(want)
    ok = abs(got - want) <= tol
    print(f"{'ok:  ' if ok else 'FAIL:'} {name}: got {got!r} want {want!r}")
    if not ok:
        fails.append(name)


def expect_fatal(name, path, needle=None):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    msg = (r.stdout + r.stderr).lower()
    bad = r.returncode == 0 or ("fatal" not in msg and "error" not in msg)
    if not bad and needle is not None:
        bad = needle.lower() not in msg
    print(f"{'ok:  ' if not bad else 'FAIL:'} {name} refused "
          f"(rc={r.returncode})")
    if bad:
        fails.append(name)


# ---- case constants (mirrors case.dict: 60 sections x 1 m, vertical) ----
N = 60
DX = 1.0
L = N * DX
DT = 0.05
RHO_L = 1000.0
G = 9.80665
US = 0.5
VSL = 0.5
VSG = 0.05
VM = VSL + VSG
FACTOR = 1.0
CONT = 0.5
UW = 2.0
TAU = 20.0
ARES = 0.05


def flux(a):
    """Wallis drift liquid flux at the imposed mixture velocity."""
    return US * a * a + (VM - US) * a


# independent reference model: the steady fill root and the fill shock
ALPHA_ROOT = (-(VM - US) + math.sqrt((VM - US) ** 2 + 4 * US * VSL)) / (2 * US)
S_SHOCK = (VSL - flux(ARES)) / (ALPHA_ROOT - ARES)
T_REF = L / UW + L / S_SHOCK


def run_raw(path):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"case {path} failed rc={r.returncode}\n"
                         f"{r.stdout}\n{r.stderr}")
    return r.stdout


def kv(rest):
    d, key = {}, None
    for w in rest:
        if key is None:
            key = w
        else:
            d[key] = w
            key = None
    return d


def parse(text):
    frames = []
    for line in text.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "t" and len(tok) >= 6:
            frames.append({"t": float(tok[1]), "vm": float(tok[3]),
                           "dp": float(tok[5]), "cells": [],
                           "purge": None, "events": []})
        elif tok[0] == "cell" and frames:
            d = kv(tok[2:])
            c = {"cell": float(tok[1])}
            for k, v in d.items():
                c[k] = float(v)
            frames[-1]["cells"].append(c)
        elif tok[0] == "purge" and len(tok) > 4 and tok[3] == "phase":
            d = kv(tok[5:])
            p = {"t": float(tok[2]), "phase": tok[4]}
            for k in ("x_wave", "u_wave", "cap", "pbase", "head", "out"):
                p[k] = float(d[k])
            p["cycles"] = int(float(d["cycles"]))
            frames[-1]["purge"] = p
        elif tok[0] == "purge" and frames:
            d = kv(tok[3:])
            e = {"t": float(tok[2]), "phase": d["phase"]}
            for k in ("x_wave", "head", "cap", "out", "duration", "period"):
                if k in d:
                    e[k] = float(d[k])
            if "cycle" in d:
                e["cycle"] = int(float(d["cycle"]))
            frames[-1]["events"].append(e)
    return frames


def strip_purge(text):
    return "".join(line + "\n" for line in text.splitlines()
                   if not line.replace("[chflow]", "").lstrip()
                   .startswith("purge "))


ACTIVE = ("TRIGGER", "PROPAGATE", "COMPLETE")

# ================= constant-mode full-cycle case =================
text = run_raw(os.path.join(CASEDIR, "case.dict"))
frames = parse(text)
check("case.dict: every frame carries a purge state line",
      sum(1 for f in frames if f["purge"] is None), 0, 0.0)
check("case.dict: 6001 report frames (reportEvery 1, 300 s)",
      len(frames), 6001, 0.0)
f0 = frames[0]["purge"]
check("case.dict: t=0 idle (phase NONE, zeroed state)",
      (f0["phase"] == "NONE" and f0["x_wave"] == 0.0
       and f0["u_wave"] == 0.0 and f0["cycles"] == 0
       and f0["out"] == 0.0), True, 0.0)

trig = [e for f in frames for e in f["events"] if e["phase"] == "TRIGGER"]
comp = [e for f in frames for e in f["events"] if e["phase"] == "COMPLETE"]
check("case.dict: three blowdowns triggered in 300 s (the third is still "
      "refilling at the end)", len(trig), 3, 0.0)
check("case.dict: two blowdowns completed (>= 2 cycles)", len(comp), 2, 0.0)
check("case.dict: cycle numbers consecutive from 1",
      [e["cycle"] for e in comp] == list(range(1, len(comp) + 1)), True, 0.0)

# trigger physics: the cap sits below the declared fraction of the static
# head of a continuous column, on the reported values
check("case.dict: every trigger fires with cap < factor rhoL g H",
      min(FACTOR * RHO_L * G * e["head"] - e["cap"] for e in trig) > 0,
      True, 0.0)
check("case.dict: every trigger sees a continuous column (H >= 0.5 L)",
      min(e["head"] for e in trig) >= CONT * L, True, 0.0)
check("case.dict: first trigger fires at the very first step",
      trig[0]["t"], DT, 1e-12)

# cycle accounting re-derived from the event fields
check("case.dict: durations match trigger-to-complete",
      max(abs(c["duration"] - (c["t"] - t["t"]))
          for c, t in zip(comp, trig)), 0.0, 0.0, 1e-9)
check("case.dict: first period is zero", comp[0]["period"], 0.0, 0.0, 1e-12)
check("case.dict: period matches consecutive completions",
      abs(comp[1]["period"] - (comp[1]["t"] - comp[0]["t"])), 0.0, 0.0, 1e-9)
check("case.dict: period matches the independent traverse+fill model",
      comp[1]["period"], T_REF, 0.15)
check("case.dict: period inside the severe-slugging band (20 s..500 s)",
      20.0 <= comp[1]["period"] <= 500.0, True, 0.0)
check("case.dict: cycle frequency inside 2 mHz..50 mHz",
      0.002 <= 1.0 / comp[1]["period"] <= 0.05, True, 0.0)

# purged volume per cycle vs the head drop, bounded by the transport inflow
vol_ok = True
for c, t in zip(comp, trig):
    dI = t["head"] - c["head"]
    if not (dI - 0.5 <= c["out"] <= dI + VM * c["duration"] + 0.5):
        vol_ok = False
check("case.dict: per-cycle out matches the head drop within the "
      "transport bound", vol_ok, True, 0.0)

outs = [f["purge"]["out"] for f in frames]
check("case.dict: out is zero while the column only accumulates",
      max(o for f, o in zip(frames, outs) if f["t"] < trig[0]["t"] - 1e-12),
      0.0, 0.0, 1e-12)
check("case.dict: cumulative out grows across completed cycles",
      all(comp[i]["out"] > 0.5 * L * (1.0 - CONT)
          and outs[-1] >= sum(c["out"] for c in comp) - 1e-6
          for i in range(len(comp))), True, 0.0)
# the state-line out minus the trigger-frame out equals the event out
out_ok = True
for c, t in zip(comp, trig):
    o_t = next(f["purge"]["out"] for f in frames if f["t"] == t["t"])
    o_c = next(f["purge"]["out"] for f in frames if f["t"] == c["t"])
    out_ok = out_ok and abs((o_c - o_t) - c["out"]) <= 1e-6
check("case.dict: event out equals the cumulative-out difference",
      out_ok, True, 0.0)

# the drained state and the front footprint on the field
drained_ok = True
below_ok = True
for f in frames:
    p = f["purge"]
    if p["phase"] == "COMPLETE":
        if max(abs(c["alpha"] - ARES) for c in f["cells"]) > 1e-12:
            drained_ok = False
    if p["phase"] in ACTIVE:
        # strict below: the cell the front sits on is the one-cell-thick
        # transitional front cell, its value is not pinned to alpha_res
        for c in f["cells"]:
            centre = (c["cell"] + 0.5) * DX
            if centre < p["x_wave"] and abs(c["alpha"] - ARES) > 1e-12:
                below_ok = False
check("case.dict: at every COMPLETE the riser mean alpha equals the "
      "residual holdup (<= alpha_res + tol)", drained_ok, True, 0.0)
check("case.dict: every cell behind the front sits at alpha_res exactly",
      below_ok, True, 0.0)

# front kinematics: declared constant speed, step-by-step integration
check("case.dict: u_wave equals the declared constant while active",
      max(f["purge"]["u_wave"] for f in frames
          if f["purge"]["phase"] in ACTIVE), UW, 1e-12)
check("case.dict: u_wave is zero while idle",
      max((f["purge"]["u_wave"] for f in frames
           if f["purge"]["phase"] == "NONE"), default=0.0), 0.0, 0.0, 1e-12)
xerr = 0.0
for a, b in zip(frames, frames[1:]):
    pa, pb = a["purge"], b["purge"]
    if pb["phase"] in ACTIVE and pa["phase"] in ACTIVE:
        xerr = max(xerr, abs(pb["x_wave"] - min(pa["x_wave"]
                                                + pb["u_wave"] * DT, L)))
    elif pb["phase"] == "TRIGGER":
        xerr = max(xerr, abs(pb["x_wave"] - min(pb["u_wave"] * DT, L)))
check("case.dict: front position integrates u_wave step by step",
      xerr, 0.0, 0.0, 1e-9)
check("case.dict: front clamps at the riser top on completion",
      min(c["x_wave"] for c in comp), L, 1e-12)

# cap relaxation: implicit Euler on the reported pbase, seeded at step 1
check("case.dict: cap seeded from the base EOS pressure at the first step",
      frames[1]["purge"]["cap"], frames[1]["purge"]["pbase"], 1e-12)
cap_err = 0.0
for a, b in zip(frames[1:], frames[2:]):
    want = (a["purge"]["cap"] / DT + b["purge"]["pbase"] / TAU) \
        / (1.0 / DT + 1.0 / TAU)
    cap_err = max(cap_err, abs(b["purge"]["cap"] - want))
check("case.dict: cap follows the implicit-Euler relaxation of pbase",
      cap_err, 0.0, 0.0, 1e-4)
# the cap flushes out during the blowdown and re-charges while accumulating
f30 = next(f for f in frames if f["t"] == 30.0)
f60 = next(f for f in frames if f["t"] == 60.0)
check("case.dict: the cap depressurises during the purge",
      f30["purge"]["cap"] < frames[1]["purge"]["cap"], True, 0.0)
check("case.dict: the cap re-charges towards the EOS pin while refilling",
      f60["purge"]["cap"] > f30["purge"]["cap"], True, 0.0)

# phase sequence NONE -> TRIGGER -> PROPAGATE... -> COMPLETE -> NONE
seq_ok = True
prev = "NONE"
for f in frames[1:]:
    s = f["purge"]["phase"]
    ok = ((s == "TRIGGER" and prev in ("NONE",))
          or (s == "PROPAGATE" and prev in ("TRIGGER", "PROPAGATE"))
          or (s == "COMPLETE" and prev in ("TRIGGER", "PROPAGATE"))
          or (s == "NONE" and prev in ("COMPLETE", "NONE")))
    if not ok:
        seq_ok = False
    prev = s
check("case.dict: phase machine TRIGGER->PROPAGATE->COMPLETE->NONE holds",
      seq_ok, True, 0.0)

# feedback: the hydrostatic head collapses during the purge, recovers after
times = [f["t"] for f in frames]


def dp_at(tq):
    i = min(range(len(times)), key=lambda k: abs(times[k] - tq))
    return frames[i]["dp"]


check("case.dict: dp collapses below a quarter of the trigger dp",
      max(dp_at(c["t"]) / dp_at(t["t"]) for c, t in zip(comp, trig)) < 0.25,
      True, 0.0)
check("case.dict: dp recovers above 80% of the previous trigger dp before "
      "the next blowdown",
      dp_at(trig[1]["t"] - 1.0) / dp_at(trig[0]["t"]) > 0.8, True, 0.0)
check("case.dict: mixture velocity pinned by the rate inlet",
      max(abs(f["vm"] - VM) for f in frames), 0.0, 0.0, 1e-12)
check("case.dict: final cycle counter equals the completions",
      frames[-1]["purge"]["cycles"], len(comp), 0.0)


def finite_all(fl):
    for f in fl:
        vals = [f["t"], f["vm"], f["dp"]] + \
            [v for v in f["purge"].values() if isinstance(v, (int, float))]
        if any(isinstance(v, float) and not math.isfinite(v) for v in vals):
            return False
        for c in f["cells"]:
            if any(not math.isfinite(v) for v in c.values()):
                return False
    return True


check("case.dict: every reported number is finite", finite_all(frames),
      True, 0.0)

# ================= closed-mode formula case =================
cframes = parse(run_raw(os.path.join(CASEDIR, "caseClosed.dict")))
ctrig = [e for f in cframes for e in f["events"] if e["phase"] == "TRIGGER"]
ccomp = [e for f in cframes for e in f["events"] if e["phase"] == "COMPLETE"]
check("caseClosed.dict: at least two closed-mode cycles complete",
      len(ccomp) >= 2, True, 0.0)
u_err = 0.0
x_err = 0.0
for a, b in zip(cframes, cframes[1:]):
    pa, pb = a["purge"], b["purge"]
    if pb["phase"] in ACTIVE:
        prior = [e for e in ctrig if e["t"] <= pb["t"] + 1e-9]
        h_trig = prior[-1]["head"] if prior else 0.0
        want = math.sqrt(2.0 * max(FACTOR * RHO_L * G * h_trig
                                   - pb["cap"], 0.0) / RHO_L)
        u_err = max(u_err, abs(pb["u_wave"] - want))
        if pa["phase"] in ACTIVE:
            x_err = max(x_err, abs(pb["x_wave"] - min(pa["x_wave"]
                                                      + pb["u_wave"] * DT,
                                                      L)))
        elif pb["phase"] == "TRIGGER":
            x_err = max(x_err, abs(pb["x_wave"] - min(pb["u_wave"] * DT, L)))
check("caseClosed.dict: u_wave equals the energy-balance closed form "
      "(trigger head, live cap) every step", u_err, 0.0, 0.0, 1e-8)
check("caseClosed.dict: closed-mode front integrates step by step",
      x_err, 0.0, 0.0, 1e-8)
check("caseClosed.dict: closed-mode front never stalls (x strictly grows "
      "while propagating)",
      all(b["purge"]["x_wave"] > a["purge"]["x_wave"]
          for a, b in zip(cframes, cframes[1:])
          if a["purge"]["phase"] == "PROPAGATE"
          and b["purge"]["phase"] in ("PROPAGATE", "COMPLETE")), True, 0.0)
check("caseClosed.dict: drained riser sits at alpha_res at every COMPLETE",
      all(max(abs(c["alpha"] - ARES) for c in f["cells"]) <= 1e-12
          for f in cframes if f["purge"]["phase"] == "COMPLETE"), True, 0.0)
check("caseClosed.dict: closed-mode triggers respect cap < factor rhoL g H",
      min(FACTOR * RHO_L * G * e["head"] - e["cap"] for e in ctrig) > 0,
      True, 0.0)
check("caseClosed.dict: every reported number is finite",
      finite_all(cframes), True, 0.0)

# ================= undeclared-path byte identity =================
base_text = run_raw(os.path.join(CASEDIR, "caseBase.dict"))
dis_text = run_raw(os.path.join(CASEDIR, "caseDisabled.dict"))
check("caseDisabled.dict (enable false) is byte-identical to caseBase.dict",
      dis_text == base_text, True, 0.0)

case_text = open(os.path.join(CASEDIR, "case.dict")).read()
with tempfile.TemporaryDirectory() as tmp:
    def variant(name, old, new):
        path = os.path.join(tmp, name + ".dict")
        with open(path, "w") as f:
            f.write(case_text.replace(old, new))
        return path

    idle_path = variant("idle", "triggerFactor 1.0", "triggerFactor 1e-6")
    idle_text = run_raw(idle_path)
    iframes = parse(idle_text)
    check("idle variant: the layer never triggers",
          len([e for f in iframes for e in f["events"]]), 0, 0.0)
    check("idle variant: the layer still runs (cap relaxes towards the "
          "EOS pin)", iframes[-1]["purge"]["cap"] > 0.0, True, 0.0)
    check("idle variant: stripping the purge lines reproduces caseBase "
          "byte-for-byte", strip_purge(idle_text) == base_text, True, 0.0)

    # ---- parameter rejections (loud cfFatal, no silent clamping) ----
    expect_fatal("unknown waveMode",
                 variant("mode", "waveMode constant", "waveMode travelling"),
                 "unknown purgeWave waveMode")
    expect_fatal("enable outside true/false",
                 variant("en", "enable true", "enable maybe"),
                 "enable must be true or false")
    expect_fatal("triggerFactor 0",
                 variant("tf", "triggerFactor 1.0", "triggerFactor 0"),
                 "triggerFactor must be positive")
    expect_fatal("continuityAlpha 1.5 outside (0, 1)",
                 variant("ca", "continuityAlpha 0.5", "continuityAlpha 1.5"),
                 "continuityAlpha must be in (0, 1)")
    expect_fatal("continuityAlpha 0 outside (0, 1)",
                 variant("ca0", "continuityAlpha 0.5", "continuityAlpha 0"),
                 "continuityAlpha must be in (0, 1)")
    expect_fatal("alphaResidual above continuityAlpha",
                 variant("ar", "alphaResidual 0.05", "alphaResidual 0.6"),
                 "below continuityAlpha")
    expect_fatal("negative alphaResidual",
                 variant("arn", "alphaResidual 0.05", "alphaResidual -0.1"),
                 "alphaResidual must be in [0, 1)")
    expect_fatal("constant mode without a positive uWave",
                 variant("uw0", "uWave 2.0", "uWave 0"),
                 "positive uWave")
    expect_fatal("negative uWave",
                 variant("uwn", "uWave 2.0", "uWave -2"),
                 "negative purgeWave uWave")
    expect_fatal("gasCapTime 0",
                 variant("tau0", "gasCapTime 20", "gasCapTime 0"),
                 "positive gasCapTime")
    expect_fatal("negative maxCycles",
                 variant("mc", "maxCycles 0", "maxCycles -1"),
                 "negative purgeWave maxCycles")
    expect_fatal("riserEnd not above riserStart",
                 variant("win", "maxCycles 0",
                         "maxCycles 0;\n    riserStart 10;\n    riserEnd 5"),
                 "riserEnd must exceed riserStart")
    expect_fatal("negative riserStart",
                 variant("rs", "maxCycles 0",
                         "maxCycles 0;\n    riserStart -1"),
                 "negative purgeWave riserStart")

expect_fatal("purgeWave without the EOS gasR (caseReject.dict)",
             os.path.join(CASEDIR, "caseReject.dict"),
             "no gasR is declared")

print()
if fails:
    print(f"FAILED {len(fails)} check(s):")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("transient purge wave (5.162): all checks passed")
