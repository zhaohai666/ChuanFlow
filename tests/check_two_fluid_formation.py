#!/usr/bin/env python3
"""Analytic verification of the radial formation conduction grid coupled to
the two-fluid temperature equation (chflow-1d-twofluid, D3 increment 1).

usage: check_two_fluid_formation.py <exe> <case.dict>

The case is steady single-phase liquid flow through a horizontal pipe whose
wall is a finite rock mass conducting radially from the bore wall r_w = D/2
out to a far field r_e held at T_e. Two self-derived anchors (cylindrical
Fourier conduction, no transcription of any published correlation):

  steady fluid   the film k'dx and the formation resistance
                 R_form = ln(r_e/r_w)/(2 pi k dx) act in series, so the fluid
                 sees g_eff = 1/(1/(k'dx) + R_form) to T_e and the upwind
                 implicit steady state is the exchanger recurrence
                   T_i - T_e = (T_in - T_e) r^(i+1),  r = mcd/(mcd + g_eff).
  steady rock    the face conductance g_j = 2 pi k dx / ln(rho_{j+1}/rho_j)
                 makes the node profile the exact logarithm
                   T_j = T_w + (T_e - T_w) ln(rho_j/r_w)/ln(r_e/r_w),
                 which on geometric nodes rho_j = r_w(r_e/r_w)^(j/M) is
                 LINEAR in j: T_j = T_w + (T_e - T_w) j/M. The wall node sits
                 at the series divider T_w = T_e + (T_i - T_e) R_form/(R_film
                 + R_form).
  energy         every frame closes intake = storage + far-field rejection
                 exactly (the internal faces telescope); at steady state the
                 storage vanishes and the fluid enthalpy loss equals the heat
                 the rock rejects to the far field.
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

exe, case = sys.argv[1], sys.argv[2]


def run(text):
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "case.dict")
        with open(p, "w") as fh:
            fh.write(text)
        return subprocess.run([exe, p], capture_output=True, text=True)


with open(case) as fh:
    original = fh.read()

proc = run(original)
if proc.returncode != 0:
    raise SystemExit(f"solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")

# ---- parse the time series ----
frames = []
cur = None


def flush():
    global cur
    if cur is not None:
        frames.append(cur)
        cur = None


for line in proc.stdout.splitlines():
    tok = line.replace("[chflow]", "").split()
    if not tok:
        continue
    if tok[0] == "t":
        flush()
        cur = {"t": float(tok[1]), "T": [], "form": {},
               "thermal": None, "formation": None}
    elif tok[0] == "cell":
        cur["T"].append(float(tok[11]))
    elif tok[0] == "thermal:":
        cur["thermal"] = (float(tok[2]), float(tok[4]))
    elif tok[0] == "formation:":
        cur["formation"] = (float(tok[2]), float(tok[4]), float(tok[6]))
    elif tok[0] == "form":
        cur["form"][int(tok[1])] = [float(x) for x in tok[2:]]
flush()
if not frames:
    raise SystemExit("no time series parsed")

fails = []
ncrit = 0


def check(ok, msg):
    global ncrit
    ncrit += 1
    if not ok:
        fails.append(msg)


# ---- case parameters (mirrored from the dict) ----
RHO_L, C_L = 1000.0, 4180.0
LENGTH, DIAMETER, SECTIONS = 2.0, 0.02, 20
K_WALL, T_E, T_IN, V_M = 20.0, 280.0, 330.0, 0.02
K_FORM, R_E, RHOC, M = 2.0, 0.2, 10000.0, 10

DX = LENGTH / SECTIONS
AREA = math.pi * DIAMETER * DIAMETER / 4.0
RW = DIAMETER / 2.0
MCD = AREA * V_M * RHO_L * C_L
KW = K_WALL * DX                                   # film conductance, W/K
R_FILM = 1.0 / KW
R_FORM = math.log(R_E / RW) / (2.0 * math.pi * K_FORM * DX)
G_FORM = 1.0 / R_FORM
G_EFF = 1.0 / (R_FILM + R_FORM)
R_REC = MCD / (MCD + G_EFF)
# geometric node radii and the linear-in-j log profile fraction
RHO = [RW * (R_E / RW) ** (j / M) for j in range(M + 1)]
FRAC = [math.log(RHO[j] / RW) / math.log(R_E / RW) for j in range(M + 1)]

last = frames[-1]
Ts = last["T"]
if len(Ts) != SECTIONS:
    raise SystemExit(f"expected {SECTIONS} cell temperatures, got {len(Ts)}")
if len(last["form"]) != SECTIONS:
    raise SystemExit(f"expected {SECTIONS} formation profiles, "
                     f"got {len(last['form'])}")

# ---- steady fluid recurrence with the series effective conductance ----
for i, T in enumerate(Ts):
    exact = T_E + (T_IN - T_E) * R_REC ** (i + 1)
    check(abs(T - exact) < 1e-9,
          f"cell {i} T {T} != series exchanger recurrence {exact}")

# ---- steady formation: wall divider + exact log profile (linear in j) ----
wall_div = R_FORM / (R_FILM + R_FORM)
for i in range(SECTIONS):
    fp = last["form"][i]
    if len(fp) != M + 1:
        raise SystemExit(f"section {i}: expected {M+1} nodes, got {len(fp)}")
    tw_exact = T_E + (Ts[i] - T_E) * wall_div
    check(abs(fp[0] - tw_exact) < 1e-9,
          f"section {i} wall node {fp[0]} != series divider {tw_exact}")
    tw = fp[0]
    worst = max(abs(fp[j] - (tw + (T_E - tw) * FRAC[j]))
                for j in range(M + 1))
    check(worst < 1e-9,
          f"section {i} log profile max residual {worst} >= 1e-9")
    # far field is the Dirichlet boundary
    check(abs(fp[M] - T_E) < 1e-9,
          f"section {i} far node {fp[M]} != T_e {T_E}")
    # monotone decreasing outward (hot fluid, cold far field)
    check(all(fp[j] > fp[j + 1] for j in range(M)),
          f"section {i} radial profile not monotone out to the far field")

# axial fluid profile monotone down the flow
check(all(Ts[i] > Ts[i + 1] for i in range(SECTIONS - 1)),
      "fluid profile not monotone along the flow")

# ---- steady energy: fluid loss = far-field rejection, storage vanishes ----
qw, qi = last["thermal"]
qf, qin, qs = last["formation"]
check(abs(qs) < 1e-9, f"steady storage {qs} did not vanish")
check(abs(qw + qf) < 1e-9, f"steady qWall {qw} != -qFar {qf}")
check(abs(qin - qf) < 1e-9, f"steady intake {qin} != far rejection {qf}")
stream_loss = MCD * (T_IN - Ts[-1])
check(abs(stream_loss - qf) < 1e-9 * max(abs(qf), 1.0),
      f"stream enthalpy loss {stream_loss} != far-field heat {qf}")
# the summed formation conductance is 2 pi k L / ln(r_e/r_w)
wall_sum = sum(last["form"][i][0] - T_E for i in range(SECTIONS))
check(abs(qf - G_FORM * wall_sum) < 1e-9 * max(abs(qf), 1.0),
      f"qFar {qf} != G_FORM * sum(T_w - T_e) {G_FORM * wall_sum}")

# ---- transient energy balance closes on EVERY frame (intake = store + far) ----
nbal = 0
for fr in frames:
    if fr["t"] == 0.0:
        continue
    qf, qin, qs = fr["formation"]
    check(abs(qin - (qs + qf)) < 1e-9 * max(abs(qin), 1.0),
          f"t={fr['t']} formation balance broken: intake {qin} != "
          f"storage {qs} + far {qf}")
    nbal += 1
# the first reported frame must carry real storage (finite rock capacity)
first = next(fr for fr in frames if fr["t"] > 0.0)
check(first["formation"][2] > 1.0,
      f"first frame storage {first['formation'][2]} W is not finite-capacity")

# ---- refusals: each malformed formation block is named loudly ----
refusals = [
    (re.sub(r"thermal\s*\{[^}]*\}", "", original), "needs a thermal"),
    (original.replace("conductivity 2;", "conductivity 0;"),
     "positive conductivity"),
    (original.replace("outerRadius 0.2;", "outerRadius 0.005;"),
     "exceed the bore"),
    (original.replace("rhoC 10000;", "rhoC 0;"), "positive rhoC"),
    (original.replace("shells 10;", "shells 0;"), "at least one"),
]
for text, needle in refusals:
    r = run(text)
    check(r.returncode != 0 and needle in (r.stdout + r.stderr),
          f"refusal '{needle}' not raised loudly (rc={r.returncode})")

if fails:
    for m in fails:
        print("FAIL:", m)
    raise SystemExit(f"formation: {len(fails)} of {ncrit} checks FAILED")
print(f"formation: all {ncrit} checks passed "
      f"(T_out {Ts[-1]:.6f} K, r {R_REC:.6f}, g_eff {G_EFF:.6f} W/K, "
      f"R_form {R_FORM:.6f} K/W, qFar {qf:.4f} W, {nbal} transient balances)")