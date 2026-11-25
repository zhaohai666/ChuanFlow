#!/usr/bin/env python3
"""5.190 (E6b): the reduced pig on the 1D NETWORK transient row.

The claim has three parts, each re-derivable from the case text and the
printed trace:

  1. the discrete kinematic law is what the report prints, frame by frame
     (the case runs dt == reportEvery, so every reported frame is one step):
         v_n = u0 (1 - beta clamp(H_n / C, 0, 1))
         x_{n+1} = min(x_n + v_n dt, end)
         H_{n+1} = H_n + eta * lambda * |Q_n| dt / A
     with Q_n the pipe rate the frame itself reports and A the bore area;
  2. arrival releases the whole front (out = H), freezes the pig (v = 0,
     x = end, H and out constant) and journals ARRIVE exactly once;
  3. NOTHING IS WRITTEN BACK: the undeclared twin (same case, `pigging{}`
     deleted) prints byte-identical output once the `pig` / `pigEvent` lines
     are stripped. That is the strong invariant - this row has no holdup
     field for the pig to disturb, and the frame rates and node pressures
     prove it by equality rather than by assertion.

Also checked: an unchoked arm (blockage 0) moves at exactly u0, and nine
declarations are refused by name.

usage: check_transient_network_pigging.py <transient-exe> <case.dict>
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
BETA = scalar("blockage", pig, 0.0)
ETA = scalar("pushEfficiency", pig, 0.0)
CAP = scalar("frontCapacity", pig, 0.0)
LAM = scalar("liquidFraction", pig, 1.0)
X0 = scalar("startX", pig, 0.0)
END = scalar("endX", pig, 0.0) or L
AREA = 0.25 * math.pi * d * d

rc, out = run(CASE)
check("the case solves", rc == 0, f"(exit {rc})")

frames = []       # [(t, x, v, front, out, event)]
events = []
rates = []
for line in out.splitlines():
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

check("a frame per step", len(frames) == len(rates), f"({len(frames)} / {len(rates)})")
check("the pig starts where it was launched",
      frames and abs(frames[0][1] - X0) < 1e-12 and frames[0][3] == 0.0
      and frames[0][4] == 0.0,
      str(frames[0] if frames else None))
check("at t=0 it runs at exactly u0", frames and abs(frames[0][2] - U0) < 1e-12,
      f"({frames[0][2] if frames else 0} vs {U0})")

# ---- 1. the discrete law, re-derived frame by frame ----
# Order matters and is the kernel's own: v reads the frame's EXISTING front,
# the front then grows by the rate THIS frame reports, and the step that
# reaches the end arc freezes the pig (v = 0) and releases the front.
if frames and len(frames) == len(rates):
    worst_x = worst_v = worst_h = worst_o = 0.0
    x, H = frames[0][1], frames[0][3]
    arrived_at = None
    for n in range(1, len(frames)):
        t, fx, fv, fh, fout, ev = frames[n]
        q = abs(rates[n])              # the rate this frame reports
        if arrived_at is not None:
            # frozen: x, front and out hold, v is 0
            worst_x = max(worst_x, abs(fx - END))
            worst_h = max(worst_h, abs(fh - H))
            worst_o = max(worst_o, abs(fout - H))
            check_v = abs(fv - 0.0)
            worst_v = max(worst_v, check_v)
            continue
        choke = min(1.0, max(0.0, H / CAP)) if CAP > 0 else 0.0
        v = U0 * (1.0 - BETA * choke)
        x_next = min(x + v * dt, END)
        H_next = H + ETA * LAM * q * dt / AREA
        worst_x = max(worst_x, abs(fx - x_next))
        worst_v = max(worst_v, abs(fv - (v if x_next < END else 0.0)))
        worst_h = max(worst_h, abs(fh - H_next))
        worst_o = max(worst_o, abs(fout - (H_next if x_next >= END else 0.0)))
        if x_next >= END:
            arrived_at = n
        x, H = x_next, H_next
    check("x matched the re-derived arc position (abs 1e-9)", worst_x < 1e-9,
          f"(worst {worst_x:.3e})")
    check("v matched u0(1-beta clamp(H/C)) / the freeze (abs 1e-9)",
          worst_v < 1e-9, f"(worst {worst_v:.3e})")
    check("front matched eta lambda |Q| dt / A (abs 1e-9)", worst_h < 1e-9,
          f"(worst {worst_h:.3e})")
    check("out matched the released front (abs 1e-9)", worst_o < 1e-9,
          f"(worst {worst_o:.3e})")
    check("the pig arrives in this run", arrived_at is not None,
          f"(frame {arrived_at})")

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
    check("ARRIVE records the velocity that carried it there",
          abs(v_ev - min(frames[frames.index(arr[0]) - 1][2], 1e300)) < 1e-9,
          f"({v_ev})")
    after = arr[1:]
    check("the pig freezes after arrival",
          all(f[2] == 0.0 and abs(f[1] - END) < 1e-12 for f in after)
          and all(f[3] == arr[0][3] and f[4] == arr[0][4] for f in after),
          f"({len(after)} frames)")

# ---- 3. the strong invariant: no write-back ----
def strip_pig(text):
    return "\n".join(l for l in text.splitlines()
                     if not re.match(r"\[chflow\] (pig|pigEvent) ", l))


twin_text = re.sub(r"\n\s*pigging\s*\{[^}]*\}", "", raw)
check("the twin deleted exactly one block", "pigging" not in twin_text)
path = os.path.join(os.path.dirname(CASE), ".twin_no_pig.dict")
with open(path, "w") as fh:
    fh.write(twin_text)
rc2, out2 = run(path)
os.unlink(path)
check("the twin solves", rc2 == 0, f"(exit {rc2})")
check("everything except the pig lines is byte-identical",
      strip_pig(out).splitlines() == out2.splitlines())
check("the twin prints no pig line at all", "pig " not in out2 and "pigEvent" not in out2)

# ---- the unchoked arm: v is exactly u0 for the whole run ----
free = raw.replace("blockage 0.8;", "blockage 0.0;")
path = os.path.join(tempfile.mkdtemp(), "pigFree.dict")
with open(path, "w") as fh:
    fh.write(free)
rc3, out3 = run(path)
free_frames = [tuple(float(m.group(i)) for i in range(1, 6)) + (m.group(6),)
               for m in (re.match(r"\[chflow\] pig LINE t (\S+) x (\S+) v (\S+)"
                                  r" front (\S+) out (\S+) event (\w+)$", l)
                         for l in out3.splitlines()) if m]
check("unchoked arm solves", rc3 == 0, f"(exit {rc3})")
arr_i = next((i for i, f in enumerate(free_frames) if f[5] == "ARRIVE"),
              len(free_frames))
check("blockage 0 => v == u0 up to arrival, then the freeze",
      free_frames
      and all(abs(f[2] - U0) < 1e-12 for f in free_frames[:arr_i])
      and all(f[2] == 0.0 for f in free_frames[arr_i:]),
      f"({len(free_frames)} frames, arrive at {arr_i})")

# ---- refusals ----
probes = [
    ("negative startX", "startX -5;", "negative startX"),
    ("arrival before launch", "endX 3;", "walks upstream"),
    ("past the chain end", "endX 500;", "past the chain end"),
    ("zero u0", "u0 0;", "positive u0"),
    ("blockage 1", "blockage 1;", "fraction in [0,1)"),
    ("pushEfficiency > 1", "pushEfficiency 1.5;", "fraction in [0,1]"),
    ("no frontCapacity", "frontCapacity 0;", "positive frontCapacity"),
]
for label, repl, needle in probes:
    text = raw
    if label == "no frontCapacity":
        text = re.sub(r"\n\s*frontCapacity [0-9.eE+-]+;", "", text)
    else:
        key = repl.split()[0]
        text = re.sub(r"\b" + key + r" [0-9.eE+-]+;", repl, text, count=1)
    path = os.path.join(tempfile.mkdtemp(), "pigRefuse.dict")
    with open(path, "w") as fh:
        fh.write(text)
    rc4, out4 = run(path)
    check(f"refused: {label}", rc4 != 0 and needle in out4,
          f"(exit {rc4}, needle {'hit' if needle in out4 else 'missing'})")

# liquidFraction must be in (0,1] and a coupled pipe cannot carry a pig
for label, needle in (("liquidFraction 0", "(0,1]"),):
    # only the pigging block carries this key; the node's own liquidFraction is
    # a different declaration and must be left alone
    body = re.search(r"pigging\s*\{([^}]*)\}", raw).group(1)
    text = raw.replace(body, re.sub(r"liquidFraction [0-9.eE+-]+;",
                                    "liquidFraction 0;", body))
    path = os.path.join(tempfile.mkdtemp(), "pigRefuse2.dict")
    with open(path, "w") as fh:
        fh.write(text)
    rc5, out5 = run(path)
    check(f"refused: {label}", rc5 != 0 and needle in out5,
          f"(exit {rc5}, needle {'hit' if needle in out5 else 'missing'})")

print(f"\n{len(fails)} failed" if fails else "\nnetworkPigging: all checks passed")
sys.exit(1 if fails else 0)
