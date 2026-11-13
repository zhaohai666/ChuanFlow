#!/usr/bin/env python3
"""Closed-form verification of the multi-stage ESP + affinity laws
(chflow-1d-steady, 5.81 / D1 increment 1) on laminar single-phase pipes.

The pump curve is H(q) = stages * (s^2 h0 - (h0-h1)(q/q1)^2) clamped at 0,
with s = pumpSpeed (relative speed) and stages = pumpStages. Series stages add
head at the same flow; the Euler affinity laws (Q ~ s, H ~ s^2) give
H_s(q) = s^2 H_1(q/s), whose quadratic term is speed-invariant. At stages 1 /
speed 1 this is the 5.33 pump law bit for bit.

Every arm is laminar, so the pipe law is the explicit linear system curve
R q (R = 128 mu L/(pi D^4)) and the checker never iterates:
  * rate-fed arms impose q = 2e-4, so dpPump = -rho g H(q) is read at a KNOWN
    flow -- the stage ladder is exactly 1:2:3 and the speed arms isolate the
    s^2 shutoff scaling pointwise (not just at one operating point);
  * closed loops (both ends 1e5) put the operating point at the positive root
    of rho g H(q) = R q, where the pump lift equals friction and net dp is 0;
  * the constant-head arm (pumpRate 0) gives H = stages s^2 h0 exactly and a
    linear q = rho g H / R.

usage: check_steady_esp.py <exe> <case.dict>
"""
import math
import os
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]


def run(exe_path, case_path):
    proc = subprocess.run([exe_path, case_path], capture_output=True, text=True)
    return proc


def parse(stdout):
    nodes, pipes = {}, {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "node" and len(tok) >= 4:
            nodes[tok[1]] = float(tok[3])
        elif tok[0] == "pipe":
            d = {}
            for j in range(2, len(tok) - 1, 2):
                d[tok[j]] = float(tok[j + 1])
            pipes[tok[1]] = d
    return nodes, pipes


proc = run(exe, case)
if proc.returncode != 0:
    raise SystemExit(f"solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
nodes, pipes = parse(proc.stdout)
if not pipes:
    raise SystemExit("no pipe lines parsed")

fails = []
ncrit = 0


def check(ok, msg):
    global ncrit
    ncrit += 1
    if not ok:
        fails.append(msg)


RHO, MU, G = 850.0, 0.1, 9.80665
L, D = 1000.0, 0.05
A = math.pi * D * D / 4.0
R = 128.0 * MU * L / (math.pi * D ** 4)   # laminar system curve d(q) = R q
RG = RHO * G
QFED = 2.0e-4
H0, H1, Q1 = 40.0, 25.0, 3.0e-4
PTERM = 1.0e6   # rate-fed terminals at 10 bar keep every source node positive


def head(q, stages, speed, h0=H0, h1=H1, q1=Q1):
    """H(q) = stages * (s^2 h0 - (h0-h1)(q/q1)^2), clamped at 0."""
    if q1 <= 0:
        return max(stages * speed * speed * h0, 0.0)
    s2 = speed * speed
    return max(stages * (s2 * h0 - (h0 - h1) * q * q / (q1 * q1)), 0.0)


def quad_root(a, b, c):
    """positive root of a q^2 + b q + c = 0 (a > 0)."""
    return (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)


# ---- head count -------------------------------------------------------------
check("22 nodes, 11 pipes" in proc.stdout,
      f"expected 22 nodes / 11 pipes, got: {proc.stdout.splitlines()[0] if proc.stdout else ''}")

# ---- rate-fed arms: pointwise curve at q = QFED ------------------------------
fed = [
    ("ESP_FedS1", 1, 1.0, "FS1", "FT1"),
    ("ESP_FedS2", 2, 1.0, "FS2", "FT2"),
    ("ESP_FedS3", 3, 1.0, "FS3", "FT3"),
    ("ESP_FedV08", 1, 0.8, "FV08S", "FV08T"),
    ("ESP_FedV12", 1, 1.2, "FV12S", "FV12T"),
    ("ESP_FedComb", 2, 1.1, "FCS", "FCT"),
]
fed_dpp = {}
for name, stages, speed, src, term in fed:
    p = pipes[name]
    H = head(QFED, stages, speed)
    dpp = -RG * H
    fed_dpp[name] = p["dpPump"]
    check(abs(p["Q"] - QFED) < 1e-12,
          f"{name}: rate {p['Q']} != imposed {QFED}")
    check(abs(p["dpPump"] - dpp) < 1e-6 * abs(dpp),
          f"{name}: dpPump {p['dpPump']} != -rho g H {dpp}")
    # the fixed terminal is honored, and the free source node follows from it:
    # p_src = p_term + R q - rho g H (dpTotal = friction - lift). The terminal
    # sits at PTERM so the lift never drives the source node negative.
    check(abs(nodes[term] - PTERM) < 1e-6,
          f"{name}: terminal {nodes[term]} != declared {PTERM}")
    p_src = PTERM + R * QFED - RG * H
    check(abs(nodes[src] - p_src) < 1e-6,
          f"{name}: source node {nodes[src]} != closed form {p_src}")
    check(40 < p["Re"] < 50, f"{name}: Re {p['Re']} left the laminar band")

# pointwise stage stacking at the SAME flow: H_N = N H_1 (speed 1) -> 1:2:3
check(abs(fed_dpp["ESP_FedS2"] - 2 * fed_dpp["ESP_FedS1"])
      < 1e-9 * abs(fed_dpp["ESP_FedS1"]),
      "stage ladder not 1:2 at fixed flow")
check(abs(fed_dpp["ESP_FedS3"] - 3 * fed_dpp["ESP_FedS1"])
      < 1e-9 * abs(fed_dpp["ESP_FedS1"]),
      "stage ladder not 1:3 at fixed flow")
# affinity: the quadratic term is speed-invariant, so two speeds differ only
# by the shutoff scaling: H(s2)-H(s1) = (s2^2-s1^2) h0 at the same q
d_aff = fed_dpp["ESP_FedV12"] - fed_dpp["ESP_FedV08"]
check(abs(d_aff - (-RG * (1.2 ** 2 - 0.8 ** 2) * H0)) < 1e-6 * abs(d_aff),
      f"affinity shutoff scaling wrong: {d_aff}")

# ---- closed loops: operating point = root of rho g H(q) = R q ----------------
loops = [
    ("ESP_LoopBase", 1, 1.0),
    ("ESP_Loop2Stage", 2, 1.0),
    ("ESP_LoopSpeed", 1, 0.8),
    ("ESP_LoopComb", 2, 1.1),
]
for name, stages, speed in loops:
    p = pipes[name]
    # rho g stages (s^2 h0 - (h0-h1)(q/q1)^2) = R q
    a = RG * stages * (H0 - H1) / (Q1 * Q1)
    c = RG * stages * speed * speed * H0
    q = quad_root(a, R, -c)
    H = head(q, stages, speed)
    check(abs(p["Q"] - q) < 1e-9 * q,
          f"{name}: rate {p['Q']} != operating root {q}")
    check(abs(p["dpPump"] + RG * H) < 1e-6 * RG * H,
          f"{name}: dpPump {p['dpPump']} != -rho g H(q*) {-RG * H}")
    check(abs(R * q + p["dpPump"]) < 1e-6 * abs(p["dpPump"]),
          f"{name}: pump lift != friction drop (closed loop)")
    check(abs(p["dp"]) < 1e-6, f"{name}: closed loop reports net dp {p['dp']}")
    check(40 < p["Re"] < 400, f"{name}: Re {p['Re']} left the laminar band")

# the baseline reproduces the 5.33 documented anchor
check(abs(pipes["ESP_LoopBase"]["Q"] - 3.085629e-4) < 1e-9,
      f"baseline operating point {pipes['ESP_LoopBase']['Q']} != 5.33 anchor "
      f"3.085629e-4 (zero-declaration drift)")

# ---- constant head: H = stages s^2 h0, q linear ------------------------------
ch = pipes["ESP_ConstShutoff"]
Hch = 3 * 0.9 ** 2 * 40.0          # = 97.2 m exactly
qch = RG * Hch / R
check(abs(Hch - 97.2) < 1e-12, "constant-head shutoff law not N s^2 h0")
check(abs(ch["Q"] - qch) < 1e-9 * qch,
      f"const-head rate {ch['Q']} != rho g H / R {qch}")
check(abs(ch["dpPump"] + RG * Hch) < 1e-6 * RG * Hch,
      f"const-head dpPump {ch['dpPump']} != -rho g * 97.2 {-RG * Hch}")
check(abs(R * qch + ch["dpPump"]) < 1e-6 * abs(ch["dpPump"]),
      "const-head pump lift != friction")
check(abs(ch["dp"]) < 1e-6, "const-head closed loop reports net dp")

# ---- loud rejections outside the declaration scope --------------------------
with open(case) as f:
    base_text = f.read()


def probe(mutate, needle, label):
    with tempfile.TemporaryDirectory() as td:
        text = mutate(base_text)
        cp = os.path.join(td, "case.dict")
        with open(cp, "w") as f:
            f.write(text)
        pr = run(exe, cp)
        check(pr.returncode != 0 and needle in (pr.stdout + pr.stderr),
              f"rejection '{label}' not loud (needle '{needle}')")


# stages/speed without a pump: strip pumpHead from the baseline loop arm
probe(lambda t: t.replace(
    "ESP_LoopBase { from LBS; to LBD; length 1000; diameter 0.05; roughness 0;\n"
    "                   pumpHead 25; pumpRate 3.0e-4; pumpShutoffHead 40; }",
    "ESP_LoopBase { from LBS; to LBD; length 1000; diameter 0.05; roughness 0;\n"
    "                   pumpStages 2; }"),
    "without a pump", "stages without pump")
# non-positive speed
probe(lambda t: t.replace("pumpSpeed 0.8; }", "pumpSpeed 0; }", 1),
      "non-positive", "zero speed")
# stages < 1
probe(lambda t: t.replace("pumpStages 2; pumpSpeed 1.0; }",
                          "pumpStages 0; pumpSpeed 1.0; }", 1),
      "pumpStages < 1", "zero stages")

if fails:
    for m in fails:
        print("FAIL:", m)
    raise SystemExit(f"esp: {len(fails)} of {ncrit} checks FAILED")
print(f"esp: all {ncrit} checks passed")