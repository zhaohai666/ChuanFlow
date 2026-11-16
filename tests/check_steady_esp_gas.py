#!/usr/bin/env python3
"""Closed-form verification of ESP gas interference (chflow-1d-steady,
5.86 / D1 increment 2) against the published Turpin (1986) correlation, on
laminar single-phase rate-fed pump arms.

The kernel multiplies the 5.81 curve by
    H_tp = H_sp * exp(-a1 * (qG/qL)),
    a1   = (346430/p_i^2)(qG/qL) - 410/p_i,     p_i in psia
(Turpin 1986, verbatim in SPE 81910-PA Eqs 5-6 and Energies 2021 Table 2),
where qG/qL = g/(1-g) is the free-gas/liquid ratio at intake and p_i is the
declared ABSOLUTE intake pressure converted with 1 psi = 6894.757293168 Pa.

This checker re-derives that multiplier INDEPENDENTLY of the kernel (same
published formula, its own arithmetic) and holds dpPump to it -- the claim is
weighed here, not taken from the log. It also pins the two structural
properties: g absent/0 is bit-identical to the 5.81 single-phase head, and the
multiplier is independent of stages (so a 2-stage degraded head is exactly 2x
the 1-stage degraded head).

usage: check_steady_esp_gas.py <exe> <case.dict>
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
    return subprocess.run([exe_path, case_path], capture_output=True, text=True)


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
    raise SystemExit(f"solver failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")
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


# ---- published constants / curve (independent of the kernel) ----------------
RHO, MU, G = 850.0, 0.1, 9.80665
L, D = 1000.0, 0.05
A = math.pi * D * D / 4.0
R = 128.0 * MU * L / (math.pi * D ** 4)   # laminar system curve d(q) = R q
RG = RHO * G
QFED = 2.0e-4
H0, H1, Q1 = 40.0, 25.0, 3.0e-4
PASCAL_PER_PSI = 6894.757293168           # 1 psi exactly, Turpin's native unit


def head_sp(q, stages=1.0, speed=1.0, h0=H0, h1=H1, q1=Q1):
    """The 5.81 multi-stage / affinity curve, before gas interference."""
    s2 = speed * speed
    if q1 <= 0:
        return max(stages * s2 * h0, 0.0)
    return max(stages * (s2 * h0 - (h0 - h1) * q * q / (q1 * q1)), 0.0)


def turb_f(g, p_pa):
    """Turpin (1986) two-phase multiplier, re-derived here."""
    r = g / (1.0 - g)                     # q_G/q_L at intake
    pi = p_pa / PASCAL_PER_PSI            # psia
    a1 = (346430.0 / (pi * pi)) * r - 410.0 / pi
    return math.exp(-a1 * r), a1, r, pi


def r_cross(pi_psia):
    """The correlation's own crossover ratio: a1 = 0 <=> r = p_i/845."""
    return pi_psia / 845.0


# ---- head count -------------------------------------------------------------
check("14 nodes, 7 pipes" in proc.stdout,
      "expected 14 nodes / 7 pipes; got first line: "
      + (proc.stdout.splitlines()[0] if proc.stdout else "<empty>"))

# ---- 1. zero-declaration is bit-identical to the 5.81 single-phase head -----
hs_base = head_sp(QFED, 1.0, 1.0)                       # = 33.3333... m
dpp_base = -RG * hs_base
check(abs(pipes["GAS_Base"]["dpPump"] - dpp_base) < 1e-9 * abs(dpp_base),
      f"GAS_Base dpPump {pipes['GAS_Base']['dpPump']} != 5.81 head "
      f"{dpp_base}")
# g absent vs g == 0 declared: identical to the last printed bit
check(pipes["GAS_Base"]["dpPump"] == pipes["GAS_G0"]["dpPump"],
      f"g-absent {pipes['GAS_Base']['dpPump']} != g==0 {pipes['GAS_G0']['dpPump']} "
      f"(short-circuit must be bit-identical)")
check(abs(pipes["GAS_G0"]["dpPump"] - dpp_base) < 1e-9 * abs(dpp_base),
      "GAS_G0 (declared g=0) left the single-phase head")

# ---- 2. every degraded arm reproduces the published multiplier at 1e-9 ------
arms = [
    ("GAS_G20",       1.0, 0.2, 1.0e6),
    ("GAS_G30",       1.0, 0.3, 1.0e6),
    ("GAS_G20HiP",    1.0, 0.2, 3.0e6),
    ("GAS_2Stage_G20", 2.0, 0.2, 1.0e6),
]
for name, stages, g, p_pa in arms:
    p = pipes[name]
    check(abs(p["Q"] - QFED) < 1e-12, f"{name}: rate {p['Q']} != imposed {QFED}")
    f, a1, r, pi = turb_f(g, p_pa)
    Htp = head_sp(QFED, stages, 1.0) * f
    expect = -RG * Htp
    check(abs(p["dpPump"] - expect) < 1e-6 * abs(expect),
          f"{name}: dpPump {p['dpPump']} != H_sp*F {expect} (F={f}, a1={a1}, "
          f"r={r}, p_i={pi} psia)")
    check(40 < p["Re"] < 50, f"{name}: Re {p['Re']} left the laminar band")

# ---- 3. degradation deepens monotonically with gas fraction at fixed p_i ----
mag = lambda n: abs(pipes[n]["dpPump"])
check(mag("GAS_G20") < mag("GAS_Base") < 1e12,
      "g=0.2 at 145 psia did not degrade below the single-phase head")
check(mag("GAS_G30") < mag("GAS_G20"),
      "g=0.3 degraded LESS than g=0.2 at the same intake pressure")

# ---- 4. pressure dependence + the correlation's own crossover sign ----------
# F's sign relative to 1 is set by whether r is above/below r_cross = p_i/845:
#   a1 = 0 at r = p_i/845;  a1 > 0 (degrade) for r above, a1 < 0 (gain) below.
for name, g, p_pa in [("GAS_G20", 0.2, 1.0e6), ("GAS_G20HiP", 0.2, 3.0e6)]:
    f, a1, r, pi = turb_f(g, p_pa)
    above = r > r_cross(pi)
    # kernel must agree with the crossover prediction of the sign of (1 - F):
    degraded = mag(name) < mag("GAS_Base")
    check(degraded == above,
          f"{name}: F-1 sign {degraded} disagrees with crossover "
          f"(r={r}, r_cross={r_cross(pi)})")
    # and the kernel's own dpPump must be on the same side as the re-derived F:
    check((mag(name) < mag("GAS_Base")) == (f < 1.0),
          f"{name}: re-derived F {f} disagrees with the observed head side")

# at the SAME gas fraction, a higher intake pressure keeps the gas compressed:
# its multiplier is strictly larger (here the HiP arm even reads above 1).
f_lo, *_ = turb_f(0.2, 1.0e6)
f_hi, *_ = turb_f(0.2, 3.0e6)
check(f_hi > f_lo,
      f"higher intake pressure degraded more: F(3e6)={f_hi} !> F(1e6)={f_lo}")
check(mag("GAS_G20HiP") > mag("GAS_G20"),
      "the 3e6-Pa arm did not hold MORE head than the 1e6-Pa arm at the same g")

# ---- 5. the multiplier is independent of stages (5.81 interaction) ----------
check(abs(pipes["GAS_2Stage_G20"]["dpPump"] - 2.0 * pipes["GAS_G20"]["dpPump"])
      < 1e-9 * abs(pipes["GAS_G20"]["dpPump"]),
      "2-stage degraded head is not exactly 2x the 1-stage degraded head")

# ---- 6. constant head (pumpRate 0) + gas: H = h0 * F, q = rho g H / R --------
ch = pipes["GAS_Const"]
f_c, *_ = turb_f(0.2, 1.0e6)
Hch = H0 * f_c
qch = RG * Hch / R
check(abs(ch["dpPump"] + RG * Hch) < 1e-6 * RG * Hch,
      f"const-head dpPump {ch['dpPump']} != -rho g h0 F {-RG * Hch}")
check(abs(ch["Q"] - qch) < 1e-9 * qch,
      f"const-head rate {ch['Q']} != rho g H / R {qch}")
check(abs(ch["dp"]) < 1e-6, f"const-head closed loop reports net dp {ch['dp']}")

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


# g >= 1 (a volume fraction must be strictly below 1)
probe(lambda t: t.replace("pumpGasFraction 0.3;", "pumpGasFraction 1.0;", 1),
      "outside [0,1)", "g == 1")
# g < 0
probe(lambda t: t.replace("pumpGasFraction 0.3;", "pumpGasFraction -0.1;", 1),
      "outside [0,1)", "g negative")
# g without a pump (strip pumpHead from the G20 arm, keep gas)
probe(lambda t: t.replace(
    "GAS_G20  { from G20S; to G20T; length 1000; diameter 0.05; roughness 0;\n"
    "               pumpHead 25; pumpRate 3.0e-4; pumpShutoffHead 40;\n"
    "               pumpGasFraction 0.2; pumpIntakePressure 1.0e6; }",
    "GAS_G20  { from G20S; to G20T; length 1000; diameter 0.05; roughness 0;\n"
    "               pumpGasFraction 0.2; pumpIntakePressure 1.0e6; }"),
      "without a pump", "gas without pump")
# g without an intake pressure (the correlation cannot price p_i)
probe(lambda t: t.replace(
    "pumpGasFraction 0.2; pumpIntakePressure 1.0e6; }",
    "pumpGasFraction 0.2; }", 1),
      "positive pumpIntakePressure", "gas without intake pressure")

if fails:
    for m in fails:
        print("FAIL:", m)
    raise SystemExit(f"esp-gas: {len(fails)} of {ncrit} checks FAILED")
print(f"esp-gas: all {ncrit} checks passed")
