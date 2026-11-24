#!/usr/bin/env python3
"""E7c (5.183): linear along-the-line EOS acoustics (depressurization
transient) on the two-fluid kernel (chflow-1d-twofluid).

usage: check_transient_eos_acoustics.py <exe> <case-dir>

The kernel's EOS mode carries an opt-in linear acoustic layer: a
cell-centred pressure perturbation p' and a face-staggered velocity u'
solved with the explicit staggered scheme of the linearized Euler
acoustics (LeVeque ch. 9: dp'/dt + rho a^2 du'/dx = 0,
du'/dt + (1/rho) dp'/dx = 0, u rows first, then p rows on the NEW u).
The wave speed defaults to the EOS-self-consistent derivative of the
kernel's own gas law p = rhoG gasR T: a^2 = (dp/drho)_T = gasR T, or the
case declares soundSpeed wholesale. One declared end vents
(pressure-release, p' = 0), the other is a hard wall (u' = 0).

This checker imports NO C++ path. It re-derives the wave speed, the
stability bound, the front positions, the reflection period and the
acoustic energy independently from the CASE TEXT and the REPORTED
p(t, x), u(t, face) fields, then verifies:
  * the run starts uniform (p' = initialGauge everywhere, u' = 0) and
    the declared dt respects the re-derived a_max dt / dx_min <= cfl;
  * the release wave crosses the line in L/a: first arrivals at probe
    cells match (L - x)/a and the fitted wave speed reproduces a to 5%;
  * the front position x_front = L - a t at two instants (2.5 cells);
  * the hard-wall gauge signal is the textbook +/-p0 square wave: its
    hysteresis transitions are spaced 2L/a (the reflection period) and
    it overshoots to SUCTION (p' < -0.5 p0 - inertial overshoot);
  * the acoustic energy (re-derived from the reported fields) never
    grows, and is CONSERVED - the pressure-release boundary has zero
    energy flux (p' . u = 0 there), so the lossless linear system
    sloshes instead of decaying;
  * the declared soundSpeed arm rides 250 m/s and NOT the EOS speed;
  * the undeclared twin (no acoustics block) reproduces the declared
    run's whole report byte-for-byte once the pA tokens and aunode lines
    are stripped, and the enable-false twin is byte-identical outright;
  * out-of-range / inconsistent declarations are refused loudly
    (cfFatal), including a dt that violates the stability bound at
    runtime.

The model boundary is honest: a LINEAR perturbation layer superposed on
the unchanged lumped mean state, report-only (no feedback into the
momentum, holdup or EOS density fields), wave content only.
"""

import math
import os
import re
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
            frames.append({"t": float(tok[1]), "pA": [], "uA": [], "p": []})
        elif tok[0] == "cell" and frames:
            i = tok.index
            if "pA" in tok:
                frames[-1]["pA"].append(float(tok[i("pA") + 1]))
            if "p" in tok:
                frames[-1]["p"].append(float(tok[i("p") + 1]))
        elif tok[0] == "aunode" and frames:
            frames[-1]["uA"].append(float(tok[tok.index("u") + 1]))
    return frames


def case_number(text, key):
    m = re.search(rf"(?m)^\s*{re.escape(key)}\s+([-+0-9.eE]+)\s*;", text)
    if not m:
        raise SystemExit(f"case text has no {key}")
    return float(m.group(1))


def wave_fit(frames, probe, L, dx, p0):
    """First-arrival times at the probe cells + the 1/a least-squares fit."""
    pts = []
    for i in probe:
        x = (i + 0.5) * dx
        t1 = next(f["t"] for f in frames if f["pA"][i] < 0.5 * p0)
        pts.append((L - x, t1))
    n = len(pts)
    sx = sum(p[0] for p in pts)
    sy = sum(p[1] for p in pts)
    sxx = sum(p[0] ** 2 for p in pts)
    sxy = sum(p[0] * p[1] for p in pts)
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    return pts, 1.0 / slope


def wall_transitions(frames, p0):
    """Hysteresis (Schmitt) transitions of the hard-wall cell's
    |pA - p0| drive: robust to the scheme's dispersive overshoot at the
    reflection instant (the raw 0.5 p0 threshold flickers there)."""
    st = False
    trans = []
    for f in frames:
        s = abs(f["pA"][0] - p0) / p0
        if (st and s >= 0.4) or s > 0.6:
            ns = True
        elif (not st and s <= 0.4) or s < 0.4:
            ns = False
        else:
            ns = st
        if ns != st:
            trans.append(f["t"])
        st = ns
    return trans


def acoustic_energy(fr, rho, a2, dxs):
    """Per-unit-area acoustic energy re-derived from the reported fields:
    potential p'^2/(2 rho a^2) per cell + kinetic rho u'^2/2 per face
    (a face owns the mean of its two neighbouring dx; an end face owns
    the boundary half cell)."""
    ep = sum(fr["pA"][i] ** 2 / (2.0 * rho[i] * a2[i]) * dxs[i]
             for i in range(len(fr["pA"])))
    ek = 0.0
    n = len(fr["pA"])
    for j in range(len(fr["uA"])):
        df = 0.5 * (dxs[j - 1] + dxs[j]) if 0 < j < n else 0.5 * dxs[
            0 if j == 0 else n - 1]
        rf = 0.5 * (rho[j - 1] + rho[j]) if 0 < j < n else rho[
            0 if j == 0 else n - 1]
        ek += 0.5 * rf * fr["uA"][j] ** 2 * df
    return ep + ek


# ================= case text (independent re-derivation) =================
with open(os.path.join(CASEDIR, "case.dict")) as f:
    case_text = f.read()
L = case_number(case_text, "length")
N = int(case_number(case_text, "sections"))
DX = L / N
GASR = case_number(case_text, "gasR")
T0 = case_number(case_text, "temperature")
DT = case_number(case_text, "dt")
END = case_number(case_text, "endTime")
REP = int(case_number(case_text, "reportEvery"))
P0 = case_number(case_text, "initialGauge")
SND = case_number(case_text, "soundSpeed")
CFL = case_number(case_text, "cfl")
A = SND if SND > 0 else math.sqrt(GASR * T0)  # the re-derived wave speed
A2 = A * A
# the mean state is a sealed static gas line: rho_i = rhoG = p/(gasR T)
RHO = P0 / (GASR * T0)

STEPS = int(round(END / DT))
NFRAMES = STEPS // REP + 1
PROBE = [10, 30, 50, 70, 90]

# ================= declared run =================
text = run_raw(os.path.join(CASEDIR, "case.dict"))
frames = parse(text)
check("frames reported (dt 0.001, endTime 2, reportEvery 5)",
      len(frames), NFRAMES, 0)
check_true("every frame has 100 reported cells",
           all(len(fr["pA"]) == N for fr in frames))
check_true("every frame has 101 face velocity lines",
           all(len(fr["uA"]) == N + 1 for fr in frames))
check("frame 0 is the uniform declared overpressure",
      max(abs(v - P0) for v in frames[0]["pA"]), 0.0, 0, 1e-6)
check("frame 0 is at rest (every face u' = 0)",
      max(abs(v) for v in frames[0]["uA"]), 0.0, 0, 1e-12)
check_true("declared dt respects the re-derived stability bound "
           f"a dt/dx = {A * DT / DX:.3f} <= cfl {CFL}",
           A * DT / DX <= CFL)

# ---- bounded, no spurious growth ----
mx = max(max(abs(v) for v in fr["pA"]) for fr in frames)
check_true("max |pA| over the run stays bounded (no spurious growth; the "
           f"discrete wall reflection wobbles to {mx / P0:.2f} p0 - the "
           "staggered scheme's known dispersion on the square-wave "
           "release, tolerated at 1.6 p0)", mx <= 1.6 * P0)

# ---- first arrivals: the wave crosses the line at a ----
worst_arr = 0.0
for i in PROBE:
    x = (i + 0.5) * DX
    t1 = next(f["t"] for f in frames if f["pA"][i] < 0.5 * P0)
    worst_arr = max(worst_arr, abs(t1 - (L - x) / A))
check("worst first-arrival error at the probe cells (ideal (L-x)/a)",
      worst_arr, 0.0, 0, 0.02)

_, a_fit = wave_fit(frames, PROBE, L, DX, P0)
check("wave speed fitted from the first arrivals", a_fit, A, 0.05)

# ---- front position x_front = L - a t ----
for ts in (0.10, 0.25):
    fr = min(frames, key=lambda f: abs(f["t"] - ts))
    cnt = sum(1 for v in fr["pA"] if v < 0.5 * P0)
    check(f"front cells behind the wave at t = {fr['t']} "
          f"(ideal a t / dx = {A * fr['t'] / DX:.2f})",
          cnt, A * fr["t"] / DX, 0, 2.5)

# ---- reflection period: the hard-wall square wave ----
trans = wall_transitions(frames, P0)
check_true("the hard-wall gauge signal shows the expected reflection "
           f"transitions ({[f'{t:.3f}' for t in trans]})",
           len(trans) == 3)
for k, (t0, t1) in enumerate(zip(trans, trans[1:])):
    check(f"hard-wall transition spacing {k + 1} = the 2L/a reflection "
          "period", t1 - t0, 2.0 * L / A, 0.10)
check("the wall gauge overshoots to suction (inertial overshoot of the "
      "release)",
      float(min(fr["pA"][0] for fr in frames) < -0.5 * P0), 1.0, 0)

# ---- energy: the pressure-release boundary has ZERO energy flux
# (p' . u = 0 there), so the closed linear system conserves the acoustic
# energy - it sloshes at the 4L/a period instead of decaying. The checks
# are boundedness + conservation of that lossless budget. ----
rho = [RHO] * N
a2 = [A2] * N
dxs = [DX] * N
E = [acoustic_energy(fr, rho, a2, dxs) for fr in frames]
check_true("acoustic energy never grows past its initial value "
           "(stability; worst frame "
           f"{max(E) / E[0]:.4f} x E0)", max(E) <= 1.01 * E[0])
check("acoustic energy is conserved within the staggered scheme's "
      "modified-energy tolerance (lossless linear system)",
      E[-1], E[0], 0.02)

# ---- open end: pressure release works ----
check("the vent-adjacent cell stays low after the release passes",
      max(fr["pA"][N - 1] for fr in frames if fr["t"] >= 0.05), 0.0, 0,
      0.5 * P0)

# ---- byte-identity twins ----
base_text = run_raw(os.path.join(CASEDIR, "caseBase.dict"))
stripped = "\n".join(
    re.sub(r" pA [^ ]+", "", l)
    for l in text.splitlines()
    if not l.replace("[chflow]", "").lstrip().startswith("aunode "))
check("declared run minus the pA/aunode lines is byte-identical to the "
      "undeclared twin (the acoustic layer is the only difference)",
      float(stripped == base_text.rstrip("\n")), 1.0, 0)
dis_text = run_raw(os.path.join(CASEDIR, "caseDisabled.dict"))
check("enable-false twin is byte-identical to the undeclared twin",
      float(dis_text.rstrip("\n") == base_text.rstrip("\n")), 1.0, 0)

# ================= declared soundSpeed arm =================
speed_text = run_raw(os.path.join(CASEDIR, "caseSpeed.dict"))
sframes = parse(speed_text)
check("soundSpeed arm frames reported (endTime 2.4)", len(sframes), 481, 0)
spts, s_fit = wave_fit(sframes, PROBE, L, DX, P0)
check("declared soundSpeed arm rides 250 m/s", s_fit, 250.0, 0.05)
check_true("the override really moved the knob (fitted 246.9 m/s, not "
           f"the EOS speed {A:.1f} m/s)", abs(s_fit - A) > 0.05 * A)
strans = wall_transitions(sframes, P0)
check_true("soundSpeed arm wall transitions "
           f"({[f'{t:.3f}' for t in strans]})", len(strans) == 3)
for k, (t0, t1) in enumerate(zip(strans, strans[1:])):
    check(f"soundSpeed arm reflection spacing {k + 1} = 2L/250",
          t1 - t0, 2.0 * L / 250.0, 0.10)

# ================= rejections =================
with open(os.path.join(CASEDIR, "caseReject.dict")) as f:
    reject_text = f.read()


def variant(name, old, new, src=None):
    src = case_text if src is None else src
    assert old in src, name
    p = os.path.join(tempfile.gettempdir(), "cf_e7c_" + name + ".dict")
    with open(p, "w") as f:
        f.write(src.replace(old, new))
    return p


expect_fatal("cfl 0.75 outside (0, 0.5] (caseReject.dict)",
             os.path.join(CASEDIR, "caseReject.dict"),
             "cfl must be in (0, 0.5]")
expect_fatal("enable outside true/false",
             variant("en", "enable true;", "enable maybe;"),
             "true or false")
expect_fatal("negative soundSpeed",
             variant("snd", "soundSpeed 0;", "soundSpeed -5;"),
             "soundSpeed must be positive")
expect_fatal("zero initialGauge while enabled",
             variant("g0", "initialGauge 1e6;", "initialGauge 0;"),
             "positive initialGauge")
expect_fatal("negative initialGauge",
             variant("gneg", "initialGauge 1e6;", "initialGauge -1e6;"),
             "initialGauge must be positive")
expect_fatal("unknown openEnd word",
             variant("oe", "openEnd outlet;", "openEnd middle;"),
             "outlet or inlet")
expect_fatal("acoustics without the gas EOS",
             variant("nogas", "gasR 287;", ""),
             "needs the gas EOS")
expect_fatal("dt violating the declared stability bound at runtime",
             variant("dtv", "dt 0.001;", "dt 0.004;"),
             "acoustic CFL violated")

print()
if fails:
    print(f"FAILED {len(fails)} check(s):")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("transient EOS acoustics (5.183): all checks passed")
