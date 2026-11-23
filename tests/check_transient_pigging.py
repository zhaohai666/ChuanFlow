#!/usr/bin/env python3
"""E6 (5.176): reduced 1D pigging (pipeline pig / scraper) on the two-fluid
kernel (chflow-1d-twofluid).

usage: check_transient_pigging.py <exe> <case-dir>

The kernel maintains a single Lagrangian pig beside the Eulerian holdup
field. It launches at the declared arc, advances at
    v = u0 (1 - beta clamp(H/C, 0, 1))
(a growing piled front chokes the pig), sweeps the field over [x_old, x_new]
and piles `pushEfficiency` of that liquid into `holdupFront`, and releases
the whole pile as `out` at ARRIVE (x = endX).

This checker imports NO C++ path. It re-derives every position, velocity and
piled volume independently from the REPORTED per-cell alpha field and the
case declarations, then verifies:
  * every frame's pig state line and the ARRIVE journal match to 1e-7;
  * the front-inventory increments equal eta times the field sweep of that
    step (liquid inventory conservation of the piled slug);
  * positions rise monotonically, v starts at u0 and falls on a growing
    front, H rises monotonically, `out` is zero until the single ARRIVE;
  * the undeclared twin (no pigging block) reproduces the declared run's
    whole report byte-for-byte once the `pig` lines are stripped, and the
    enable-false twin is byte-identical outright;
  * the blockage-0 arm cruises at exactly u0 and arrives in endX/u0 seconds
    (one knob moves only the velocity law);
  * the global field liquid balance still holds (the pig never touches it);
  * out-of-range / inconsistent declarations are refused loudly (cfFatal).

The model boundary is honest: a declared kinematic reduced law, no bypass
flow, seal failure or bend resistance, and no write-back to holdup.
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


def check_true(name, cond):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name}: {cond}")
    if not cond:
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


# ---- case constants (mirror case.dict) ----
DX = 1.0
N = 20
TOTAL = N * DX
DT = 0.05
U0 = 2.0
BETA = 0.5
ETA = 0.6
CAP = 8.0
START = 0.0
END = 20.0


def run_raw(path):
    r = subprocess.run([EXE, path], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"case {path} failed rc={r.returncode}\n"
                         f"{r.stdout}\n{r.stderr}")
    return r.stdout


def parse(text):
    frames = []
    for line in text.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "t" and len(tok) >= 6:
            frames.append({"t": float(tok[1]), "cells": [], "states": [],
                           "arrives": []})
        elif tok[0] == "cell" and frames:
            c = {}
            i = tok.index
            c["alpha"] = float(tok[i("alpha") + 1])
            c["vsl"] = float(tok[i("vsl") + 1])
            frames[-1]["cells"].append(c)
        elif tok[0] == "pig" and frames:
            i = tok.index
            ev = {
                "t": float(tok[i("t") + 1]),
                "x": float(tok[i("x") + 1]),
                "v": float(tok[i("v") + 1]),
                "H": float(tok[i("holdupFront") + 1]),
                "out": float(tok[i("outVol") + 1]),
                "event": tok[i("event") + 1],
            }
            if ev["event"] == "ARRIVE":
                frames[-1]["arrives"].append(ev)
            else:
                frames[-1]["states"].append(ev)
    return frames


def swept(fr, x0, x1):
    """Independently integrate reported alpha over [x0, x1] at this frame."""
    lo, hi = max(0.0, x0), min(TOTAL, x1)
    if hi <= lo:
        return 0.0
    inv = 0.0
    for i in range(N):
        o = min(hi, (i + 1) * DX) - max(lo, i * DX)
        if o > 0.0:
            inv += fr["cells"][i]["alpha"] * o
    return inv


def replay(frames, u0, beta, eta, cap, end):
    """Re-derive every state line and ARRIVE event from the reported field."""
    x, H, out = START, 0.0, 0.0
    arrived = False
    states, arrives = [], []
    incs = []
    for n, fr in enumerate(frames):
        if n == 0:
            states.append(dict(x=x, v=u0, H=H, out=out, event="RUN"))
            continue
        if arrived:
            states.append(dict(x=end, v=0.0, H=H, out=out, event="NONE"))
            continue
        choke = min(1.0, max(0.0, H / cap)) if cap > 0.0 else 0.0
        v = u0 * (1.0 - beta * choke)
        x1 = min(x + v * DT, end)
        dH = eta * swept(fr, x, x1)
        incs.append((fr["t"], x, x1, dH))
        H += dH
        x = x1
        if x >= end:
            out = H
            arrived = True
            states.append(dict(x=end, v=0.0, H=H, out=out, event="NONE"))
            arrives.append(dict(t=fr["t"], x=end, v=v, H=H, out=out,
                                event="ARRIVE"))
        else:
            states.append(dict(x=x, v=v, H=H, out=out, event="RUN"))
    return states, arrives, incs


def compare_states(frames, ref):
    """Compare reported pig state lines to the independent replay."""
    max_e = dict(x=0.0, v=0.0, H=0.0, out=0.0)
    n = 0
    for fr, r in zip(frames, ref):
        st = fr["states"]
        if len(st) != 1:
            return f"frame {fr['t']}: {len(st)} state lines (want 1)", False
        s = st[0]
        if s["event"] != r["event"]:
            return f"frame {fr['t']}: event {s['event']} vs {r['event']}", False
        for k in ("x", "v", "H", "out"):
            e = abs(s[k] - r[k])
            max_e[k] = max(max_e[k], e)
            if e > 1e-7 * max(1.0, abs(r[k])):
                return f"frame {fr['t']} field {k}: {s[k]} vs {r[k]}", False
        n += 1
    return "", True, max_e, n


def compare_arrives(frames, ref):
    rep = [e for fr in frames for e in fr["arrives"]]
    if len(rep) != len(ref):
        return f"ARRIVE count {len(rep)} vs {len(ref)}", False
    for a, b in zip(rep, ref):
        if a["event"] != b["event"]:
            return "ARRIVE event token mismatch", False
        for k in ("t", "x", "v", "H", "out"):
            if abs(a[k] - b[k]) > 1e-7 * max(1.0, abs(b[k])):
                return f"ARRIVE field {k}: {a[k]} vs {b[k]}", False
    return "", True


# ================= declared run =================
text = run_raw(os.path.join(CASEDIR, "case.dict"))
frames = parse(text)
check("frames reported (dt 0.05, endTime 20, reportEvery 1)",
      len(frames), 401, 0)
check_true("every frame has 20 reported cells",
           all(len(fr["cells"]) == N for fr in frames))

states, arrives, incs = replay(frames, U0, BETA, ETA, CAP, END)
bad, ok, max_e, n = compare_states(frames, states)
check("every pig state line re-derived independently from the field",
      float(ok), 1.0, 0)
if not ok:
    fails.append("state lines: " + bad)
else:
    print(f"       max state errors x/v/H/out = "
          f"{max_e['x']:.2e}/{max_e['v']:.2e}/"
          f"{max_e['H']:.2e}/{max_e['out']:.2e} over {n} frames")
bad_a, ok_a = compare_arrives(frames, arrives)
check("the single ARRIVE event (time, arc, velocity, released volume)",
      float(ok_a), 1.0, 0)
if not ok_a:
    fails.append("ARRIVE: " + bad_a)

# ---- front-inventory conservation: dH = eta * field sweep ----
rep_states = [fr["states"][0] for fr in frames if fr["states"]]
worst_inc = 0.0
for e in incs:
    t, x0, x1, dH = e
    idx = int(round(t / DT))
    if idx <= 0 or idx >= len(rep_states):
        continue
    dH_rep = rep_states[idx]["H"] - rep_states[idx - 1]["H"]
    worst_inc = max(worst_inc, abs(dH_rep - dH))
check("piled-front increments equal eta x the reported field sweep",
      worst_inc, 0.0, 0, 1e-7)
check("H and position rise monotonically",
      float(all(rep_states[i]["H"] >= rep_states[i - 1]["H"] - 1e-12
                and rep_states[i]["x"] >= rep_states[i - 1]["x"] - 1e-12
                for i in range(1, len(rep_states)))), 1.0, 0)
v_run = [s["v"] for s in rep_states if s["event"] == "RUN"]
check("the pig launches at exactly the declared u0", v_run[0], U0, 0, 1e-9)
check("a growing front slows the pig below u0",
      float(any(v < U0 - 0.05 for v in v_run)), 1.0, 0)
check("velocity stays positive (blockage < 1, no stall)",
      float(all(v > 0 for v in v_run)), 1.0, 0)
check("out is zero until ARRIVE",
      float(all(s["out"] == 0.0 for s in rep_states
                if s["event"] == "RUN")), 1.0, 0)

# ---- release conservation: outVol == final piled H == sum of increments ----
check("released volume equals the final piled front inventory",
      float(len(arrives) == 1 and abs(arrives[0]["out"]
                                      - rep_states[-1]["H"]) < 1e-7), 1.0, 0)
check("released volume equals the sum of all field-sweep increments",
      float(len(arrives) == 1
            and abs(arrives[0]["out"] - sum(e[3] for e in incs)) < 1e-7),
      1.0, 0)
check("the pig arrives at the declared end arc",
      float(len(arrives) == 1 and abs(arrives[0]["x"] - END) < 1e-9),
      1.0, 0)

# ---- global field liquid balance (pig is read-only) ----
def imposed_vsl(t_mid):
    return 0.8 if 3.0 <= t_mid < 7.0 else 0.5


injected = sum(imposed_vsl(n * DT - 0.5 * DT) * DT
               for n in range(1, len(frames)))
discharged = sum(fr["cells"][N - 1]["vsl"] * DT
                 for fr in frames[1:])
inv_end = sum(frames[-1]["cells"][i]["alpha"] * DX for i in range(N))
inv_0 = sum(frames[0]["cells"][i]["alpha"] * DX for i in range(N))
# The Godunov interface flux differs from the reported cell-centre closure
# flux while a holdup pulse crosses a face, so the telescoped balance is
# exact only in the steady limit; 1e-2 m per m2 of bore over 20 s (about 1
# percent of the injected liquid) is the documented tolerance for this
# pulsed case. The byte-identity witness below is the strict proof that the
# pig never touches the field.
check("global field liquid balance (injected - discharged == inventory "
      "change, pig untouched, documented tolerance 1e-2)",
      injected - discharged, inv_end - inv_0, 0, 1e-2)

# ---- undeclared path byte-identity ----
base_text = run_raw(os.path.join(CASEDIR, "caseBase.dict"))
stripped = "\n".join(l for l in text.splitlines()
                     if not l.replace("[chflow]", "").lstrip().startswith("pig "))
check("declared run minus the pig lines is byte-identical to the "
      "undeclared twin", float(stripped == base_text.rstrip("\n")), 1.0, 0)
dis_text = run_raw(os.path.join(CASEDIR, "caseDisabled.dict"))
check("enable-false twin is byte-identical to the undeclared twin",
      float(dis_text.rstrip("\n") == base_text.rstrip("\n")), 1.0, 0)

# ================= one-knob arm: blockage 0 =================
nb_text = run_raw(os.path.join(CASEDIR, "caseNoBlock.dict"))
nb_frames = parse(nb_text)
nb_states, nb_arrives, nb_incs = replay(nb_frames, U0, 0.0, ETA, CAP, END)
bad, ok, _, _ = compare_states(nb_frames, nb_states)
check("blockage-0 arm: every state line re-derived (v = u0 throughout)",
      float(ok), 1.0, 0)
if not ok:
    fails.append("no-block states: " + bad)
nb_v = [s["v"] for s in nb_states if s["event"] == "RUN"]
check("blockage-0 arm: the pig cruises at exactly u0",
      float(all(abs(v - U0) < 1e-12 for v in nb_v)), 1.0, 0)
check("blockage-0 arm: arrival time is exactly endX/u0 = 10 s",
      float(len(nb_arrives) == 1 and abs(nb_arrives[0]["t"] - 10.0) < DT),
      1.0, 0)
check("blockage-0 arm: the front still piles liquid",
      float(len(nb_arrives) == 1 and nb_arrives[0]["out"] > 0), 1.0, 0)

# ================= rejections =================
with open(os.path.join(CASEDIR, "case.dict")) as f:
    case_text = f.read()
with tempfile.TemporaryDirectory(prefix="cf_pig_") as td:
    def variant(name, old, new):
        assert old in case_text, name
        p = os.path.join(td, name + ".dict")
        with open(p, "w") as f:
            f.write(case_text.replace(old, new))
        return p

    expect_fatal("blockage 1.5 outside [0, 1) (caseReject.dict)",
                 os.path.join(CASEDIR, "caseReject.dict"),
                 "blockage must be in [0, 1)")
    expect_fatal("negative blockage",
                 variant("bneg", "blockage 0.5;", "blockage -0.1;"),
                 "blockage must be in [0, 1)")
    expect_fatal("negative pigVelocity",
                 variant("uneg", "pigVelocity 2.0;", "pigVelocity -1;"),
                 "negative pigging pigVelocity")
    expect_fatal("zero pigVelocity while enabled",
                 variant("u0", "pigVelocity 2.0;", "pigVelocity 0;"),
                 "positive pigVelocity")
    expect_fatal("pushEfficiency > 1",
                 variant("eta", "pushEfficiency 0.6;",
                         "pushEfficiency 1.5;"),
                 "pushEfficiency must be in [0, 1]")
    expect_fatal("non-positive frontCapacity",
                 variant("cap", "frontCapacity 8;", "frontCapacity 0;"),
                 "positive frontCapacity")
    expect_fatal("negative startX",
                 variant("sx", "startX 0;", "startX -1;"),
                 "negative pigging startX")
    expect_fatal("negative endX",
                 variant("ex", "endX 20;", "endX -1;"),
                 "negative pigging endX")
    expect_fatal("endX not exceeding startX",
                 variant("ord", "startX 0;\n    endX 20;",
                         "startX 10;\n    endX 5;"),
                 "endX must exceed startX")
    expect_fatal("enable outside true/false",
                 variant("en", "enable true;", "enable maybe;"),
                 "true or false")

print()
if fails:
    print(f"FAILED {len(fails)} check(s):")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("transient pigging (5.176): all checks passed")
