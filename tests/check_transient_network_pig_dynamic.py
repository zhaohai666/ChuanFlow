#!/usr/bin/env python3
"""5.235 (E6c): the DYNAMIC pig on the 1D network transient row.

mass > 0 switches the kinematic point (5.190) to a momentum ODE integrated over
the committed step. The claim, re-derivable from the case text and the trace:

  1. the discrete momentum law is what the report prints, frame by frame
     (dt == reportEvery, so every reported frame is one step):
         Fdrive = driveFraction (p_from - p_to) A      [the frame's OWN pressures]
         Fslug  = kSlug H_n A                          [the front at step start]
         v_{n+1}= max(0, (m v_n/dt + Fdrive - f0 - Fslug) / (m/dt + f1))
         x_{n+1}= min(x_n + v_{n+1} dt, end)
         H_{n+1}= H_n + eta lambda (x_{n+1} - x_n)     [the liquid the PIG sweeps]
     u0 is the launch speed v_0; the front slug grows by the pig's swept distance
     (so a stalled pig stops piling liquid) and resists it as kSlug H A;
  2. arrival releases the whole front (out = H), freezes the pig (v = 0, x = end,
     H and out constant) and journals ARRIVE once, recording the velocity that
     carried it there;
  3. NOTHING IS WRITTEN BACK: the undeclared twin (pigging{} deleted) prints
     byte-identical output once the pig / pigEvent lines are stripped.

Also checked: a stalled-pig arm (huge slugResistance) decelerates to v = 0, never
arrives, and its front freezes; the developing trajectory accelerates then the
growing slug slows it; six dynamic declarations are refused by name.

usage: check_transient_network_pig_dynamic.py <transient-exe> <case.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) < 3:
    raise SystemExit(__doc__)

EXE, CASE = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def run(path):
    p = subprocess.run([EXE, path], capture_output=True, text=True,
                       cwd=os.path.dirname(path) or ".")
    return p.returncode, p.stdout + p.stderr


raw = open(CASE).read()
CASE_TEXT = re.sub(r"//[^\n]*", "", raw)


def block_of(name, text=CASE_TEXT):
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", text)
    if not m:
        return None
    i, depth = m.end(), 1
    while i < len(text) and depth:
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
        i += 1
    return text[m.end():i - 1]


def scalar(key, text, default=None):
    m = re.search(r"\b" + re.escape(key) + r"\s+([0-9.eE+-]+)\s*;", text)
    return float(m.group(1)) if m else default


pig = block_of("pigging")
pipe = block_of("LINE")
d = float(scalar("diameter", pipe))
L = float(scalar("length", pipe))
dt = float(scalar("dt", block_of("time")))
U0 = scalar("u0", pig)
ETA = scalar("pushEfficiency", pig, 0.0)
LAM = scalar("liquidFraction", pig, 1.0)
X0 = scalar("startX", pig, 0.0)
END = scalar("endX", pig, 0.0) or L
MASS = scalar("mass", pig, 0.0)
F0 = scalar("frictionForce", pig, 0.0)
F1 = scalar("frictionViscous", pig, 0.0)
KSLUG = scalar("slugResistance", pig, 0.0)
DRIVE = scalar("driveFraction", pig, 1.0)
AREA = 0.25 * math.pi * d * d

rc, out = run(CASE)
check("the case solves", rc == 0, f"(exit {rc})")
check("the case declares the dynamic mode (mass > 0)", MASS > 0, f"(mass {MASS})")


def parse(text):
    frames, events, rates, tank, outp = [], [], [], [], []
    for line in text.splitlines():
        m = re.match(r"\[chflow\] pig LINE t (\S+) x (\S+) v (\S+) front (\S+)"
                     r" out (\S+) event (\w+)$", line)
        if m:
            frames.append((float(m.group(1)), float(m.group(2)),
                           float(m.group(3)), float(m.group(4)),
                           float(m.group(5)), m.group(6)))
            continue
        m = re.match(r"\[chflow\] pigEvent LINE ARRIVE t (\S+) x (\S+) v (\S+)"
                     r" front (\S+) out (\S+)$", line)
        if m:
            events.append(tuple(float(m.group(i)) for i in range(1, 6)))
            continue
        m = re.match(r"\[chflow\] pipe LINE Q (\S+)$", line)
        if m:
            rates.append(float(m.group(1)))
            continue
        m = re.match(r"\[chflow\] node TANK p (\S+) ", line)
        if m:
            tank.append(float(m.group(1)))
            continue
        m = re.match(r"\[chflow\] node OUT p (\S+) ", line)
        if m:
            outp.append(float(m.group(1)))
    return frames, events, rates, tank, outp


frames, events, rates, tank, outp = parse(out)
check("a frame per step", len(frames) == len(rates) == len(tank) == len(outp),
      f"({len(frames)}/{len(rates)}/{len(tank)}/{len(outp)})")
check("the pig starts where it was launched",
      frames and abs(frames[0][1] - X0) < 1e-12 and frames[0][3] == 0.0
      and frames[0][4] == 0.0, str(frames[0] if frames else None))
check("at t=0 it runs at the launch speed u0",
      frames and abs(frames[0][2] - U0) < 1e-12,
      f"({frames[0][2] if frames else 0} vs {U0})")

# ---- 1. the discrete momentum law, re-derived frame by frame ----
# v reads the frame's OWN committed pressure drop and the front at the step's
# start; the front then grows by the distance the pig actually sweeps this step.
if frames and len(frames) == len(tank):
    worst_x = worst_v = worst_h = worst_o = 0.0
    x, H, vprev = frames[0][1], frames[0][3], frames[0][2]
    arrived_at = None
    vmax = vprev
    for n in range(1, len(frames)):
        t, fx, fv, fh, fout, ev = frames[n]
        if arrived_at is not None:
            worst_x = max(worst_x, abs(fx - END))
            worst_h = max(worst_h, abs(fh - H))
            worst_o = max(worst_o, abs(fout - H))
            worst_v = max(worst_v, abs(fv - 0.0))
            continue
        dp = tank[n] - outp[n]                 # the frame's own committed drop
        fdrive = DRIVE * dp * AREA
        fslug = KSLUG * H * AREA               # front at the step's start
        v = (MASS * vprev / dt + fdrive - F0 - fslug) / (MASS / dt + F1)
        if v < 0:
            v = 0.0
        x_next = min(x + v * dt, END)
        H_next = H + ETA * LAM * (x_next - x)  # the liquid the pig sweeps
        worst_x = max(worst_x, abs(fx - x_next))
        worst_v = max(worst_v, abs(fv - (v if x_next < END else 0.0)))
        worst_h = max(worst_h, abs(fh - H_next))
        worst_o = max(worst_o, abs(fout - (H_next if x_next >= END else 0.0)))
        vmax = max(vmax, v)
        if x_next >= END:
            arrived_at = n
        x, H, vprev = x_next, H_next, v
    check("x matched the re-derived arc position (abs 1e-9)", worst_x < 1e-9,
          f"(worst {worst_x:.3e})")
    check("v matched the implicit-Euler momentum ODE / the freeze (abs 1e-9)",
          worst_v < 1e-9, f"(worst {worst_v:.3e})")
    check("front matched eta lambda * swept distance (abs 1e-9)", worst_h < 1e-9,
          f"(worst {worst_h:.3e})")
    check("out matched the released front (abs 1e-9)", worst_o < 1e-9,
          f"(worst {worst_o:.3e})")
    check("the pig arrives in this run", arrived_at is not None,
          f"(frame {arrived_at})")
    check("the pig accelerates off the launch (peak v > u0)", vmax > U0 * 1.05,
          f"(peak {vmax:.4g} vs u0 {U0})")
    # the growing slug slows it: the arrival speed is below the peak
    if arrived_at is not None and arrived_at >= 2:
        check("the growing slug decelerates it (arrival v < peak v)",
              events and events[0][2] < vmax * 0.999,
              f"(arrival {events[0][2] if events else 0:.4g} vs peak {vmax:.4g})")

# ---- 2. the arrival statement ----
check("one ARRIVE journal entry", len(events) == 1, f"({len(events)})")
arr = [f for f in frames if f[5] == "ARRIVE"]
check("every ARRIVE frame sits exactly at the end arc",
      all(abs(f[1] - END) < 1e-12 for f in arr))
check("out stays 0 before arrival",
      all(f[4] == 0.0 for f in frames[:next((i for i, f in enumerate(frames)
                                             if f[5] == "ARRIVE"), len(frames))]))
if events and arr:
    t_ev, x_ev, v_ev, h_ev, o_ev = events[0]
    check("ARRIVE time is the frame that arrived",
          abs(t_ev - arr[0][0]) < 1e-12, f"({t_ev} vs {arr[0][0]})")
    check("ARRIVE releases the whole front", abs(o_ev - h_ev) < 1e-12
          and abs(o_ev - arr[0][3]) < 1e-9 and abs(h_ev - arr[0][4]) < 1e-9,
          f"({o_ev} / {h_ev})")
    check("ARRIVE records a positive carrying velocity", v_ev > 0, f"({v_ev})")
    after = arr[1:]
    check("the pig freezes after arrival",
          all(f[2] == 0.0 and abs(f[1] - END) < 1e-12 for f in after)
          and all(f[3] == arr[0][3] and f[4] == arr[0][4] for f in after),
          f"({len(after)} frames)")
    # the released front is exactly the swept liquid over the whole arc
    check("the released front is eta lambda (end - startX)",
          abs(o_ev - ETA * LAM * (END - X0)) < 1e-9,
          f"({o_ev} vs {ETA * LAM * (END - X0)})")

# ---- 3. the strong invariant: no write-back ----
def strip_pig(text):
    return "\n".join(l for l in text.splitlines()
                     if not re.match(r"\[chflow\] (pig|pigEvent) ", l))


twin_text = re.sub(r"\n\s*pigging\s*\{[^}]*\}", "", raw)
check("the twin deleted exactly one block", "pigging" not in twin_text)
path = os.path.join(os.path.dirname(CASE), ".twin_no_pig_dyn.dict")
with open(path, "w") as fh:
    fh.write(twin_text)
rc2, out2 = run(path)
os.unlink(path)
check("the twin solves", rc2 == 0, f"(exit {rc2})")
check("everything except the pig lines is byte-identical",
      strip_pig(out).splitlines() == out2.splitlines())
check("the twin prints no pig line at all", "pig " not in out2 and "pigEvent" not in out2)

# ---- the stalled-pig arm: a huge slug resistance chokes it to rest ----
stall = re.sub(r"slugResistance [0-9.eE+-]+;", "slugResistance 5000;", raw)
path = os.path.join(tempfile.mkdtemp(), "pigStall.dict")
with open(path, "w") as fh:
    fh.write(stall)
rc3, out3 = run(path)
sf, se, sr, st, so = parse(out3)
check("stalled arm solves", rc3 == 0, f"(exit {rc3})")
check("stalled arm never arrives", len(se) == 0
      and all(f[5] == "RUN" for f in sf), f"({len(se)} events)")
# once v hits 0 it stays 0 and the front freezes (the pig sweeps no more liquid)
zero_i = next((i for i, f in enumerate(sf) if f[2] == 0.0), None)
check("stalled arm decelerates to rest", zero_i is not None and zero_i > 0,
      f"(first v=0 at frame {zero_i})")
if zero_i is not None:
    rest = sf[zero_i:]
    check("a stalled pig stays at rest (v == 0)", all(f[2] == 0.0 for f in rest))
    check("a stalled pig's front freezes (swept-distance law)",
          all(abs(f[3] - rest[0][3]) < 1e-12 for f in rest)
          and all(abs(f[1] - rest[0][1]) < 1e-12 for f in rest),
          f"(front {rest[0][3]:.4g}, x {rest[0][1]:.4g})")
    check("a stalled pig stops short of the end arc", rest[0][1] < END - 1e-9,
          f"(x {rest[0][1]:.4g} vs end {END})")

# ---- refusals: the dynamic declarations are validated by name ----
probes = [
    ("negative mass", "mass -5;", "finite non-negative kg"),
    ("negative frictionForce", "frictionForce -1;", "non-negative N"),
    ("negative frictionViscous", "frictionViscous -1;", "non-negative N.s/m"),
    ("negative slugResistance", "slugResistance -1;", "non-negative Pa/m"),
    ("driveFraction > 1", "driveFraction 1.5;", "fraction in [0,1]"),
    ("driveFraction < 0", "driveFraction -0.1;", "fraction in [0,1]"),
]
for label, repl, needle in probes:
    key = repl.split()[0]
    text = re.sub(r"\b" + key + r" [0-9.eE+-]+;", repl, raw, count=1)
    path = os.path.join(tempfile.mkdtemp(), "pigDynRefuse.dict")
    with open(path, "w") as fh:
        fh.write(text)
    rc4, out4 = run(path)
    check(f"refused: {label}", rc4 != 0 and needle in out4,
          f"(exit {rc4}, needle {'hit' if needle in out4 else 'missing'})")

print(f"\n{len(fails)} failed" if fails else
      "\nnetworkPigDynamic: all checks passed")
sys.exit(1 if fails else 0)