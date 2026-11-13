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

D3 increment 2 (5.84) adds ONE declaration, `formation{ semiInfinite yes; }`,
and re-verifies the same shipped case plus five derived arms (all written to
a system temp dir). The outer field then answers the recorded per-step
far-field flux by a line-source Duhamel superposition,
  dT_e(t_n) = sum_m (q'_m - q'_{m-1}) F(t_n - t_m),
  F(s) = E1(r_e^2/(4 alpha s))/(4 pi k),  F(0) = 0,
so the current step's own flux cannot feed back into its own boundary. The
E1 probes are weighed against an INDEPENDENT adaptive-Simpson quadrature of
int_0^inf exp(-x exp(w)) dw, and every rise against a reconstruction from the
printed flux history -- the C++ series/continued-fraction kernel is never
trusted. A geothermal arm proves the rise is global while each section keeps
its own depth-dependent far-field offset; a scope twin proves the word is
ignored outside formation{}; the shipped case remains the bit-identical
zero-declaration witness.
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

# ---- parse the time series (shared by every derived arm) ----
def parse_series(text):
    frames = []
    cur = None

    def flush():
        nonlocal cur
        if cur is not None:
            frames.append(cur)
            cur = None

    for line in text.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "t":
            flush()
            cur = {"t": float(tok[1]), "T": [], "form": {},
                   "thermal": None, "formation": None,
                   "farRise": None, "e1": None}
        elif tok[0] == "cell":
            cur["T"].append(float(tok[11]))
        elif tok[0] == "thermal:":
            cur["thermal"] = (float(tok[2]), float(tok[4]))
        elif tok[0] == "formation:":
            cur["formation"] = (float(tok[2]), float(tok[4]), float(tok[6]))
            cur["farRise"] = float(tok[8]) if len(tok) > 8 else 0.0
        elif tok[0] == "formation" and len(tok) > 2 and tok[1] == "e1:":
            cur["e1"] = {float(tok[2 + 2 * j]): float(tok[3 + 2 * j])
                         for j in range((len(tok) - 2) // 2)}
        elif tok[0] == "form":
            cur["form"][int(tok[1])] = [float(x) for x in tok[2:]]
    flush()
    if not frames:
        raise SystemExit("no time series parsed")
    return frames


frames = parse_series(proc.stdout)

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

# =====================================================================
# 5.84 (D3 increment 2): the semi-infinite line-source outer field.
# All arms are derived from the ONE shipped case into a system temp dir;
# the shipped case itself stays the zero-declaration witness above.
# =====================================================================
TIME = "time { dt 1; endTime 2000; reportEvery 100; }"
SEMI = original.replace("shells 10;", "shells 10;\n    semiInfinite yes;")
FINE = "time { dt 1; endTime 30; reportEvery 1; }"
SEMI_FINE = SEMI.replace(TIME, FINE)
PLAIN_FINE = original.replace(TIME, FINE)
GEO = (SEMI.replace("angleDeg 0;", "angleDeg 30;")
           .replace("wallModel constant;",
                    "wallModel geothermal; surfaceTemperature 280; "
                    "gradient 0.03; inletDepth 1000;")
           .replace("wallTemperature 280;", ""))
# The same word, but declared in thermal{} instead of formation{}: not a
# partner declaration, so the whole run must be bit-identical to the base.
TOPLEVEL = original.replace("thermal\n{", "thermal\n{\n    semiInfinite yes;")


def series(text, tag):
    rr = run(text)
    if rr.returncode != 0:
        raise SystemExit(f"{tag} arm failed rc={rr.returncode}\n"
                         f"{rr.stdout}\n{rr.stderr}")
    return rr, parse_series(rr.stdout)


# ---- E1 by an INDEPENDENT adaptive-Simpson quadrature -----------------
# E1(x) = int_x^inf exp(-u)/u du = int_0^inf exp(-x exp(w)) dw (u = x e^w).
# The C++ kernel is a series / continued fraction; nothing here shares it.
def adaptive_simpson(f, a, b, rel=1e-13, depth=40):
    def rec(lo, hi, flo, fmid, fhi, whole, eps, d):
        mid = 0.5 * (lo + hi)
        lmid, rmid = 0.5 * (lo + mid), 0.5 * (mid + hi)
        flm, frm = f(lmid), f(rmid)
        left = (mid - lo) / 6.0 * (flo + 4.0 * flm + fmid)
        right = (hi - mid) / 6.0 * (fmid + 4.0 * frm + fhi)
        delta = left + right - whole
        if d <= 0 or abs(delta) <= 15.0 * eps:
            return left + right + delta / 15.0
        return (rec(lo, mid, flo, flm, fmid, left, eps / 2, d - 1)
                + rec(mid, hi, fmid, frm, fhi, right, eps / 2, d - 1))

    mid = 0.5 * (a + b)
    fa, fm, fb = f(a), f(mid), f(b)
    whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb)
    # RELATIVE tolerance: E1(50) is 3.8e-24, so an absolute 1e-14 would
    # accept the first (very wrong) Simpson panel without refining at all.
    eps = rel * max(abs(whole), 1e-300)
    return rec(a, b, fa, fm, fb, whole, eps, depth)


def E1(x):
    wmax = math.log(80.0 / x) + 1.0
    return adaptive_simpson(lambda w: math.exp(-x * math.exp(w)), 0.0, wmax)


ALPHA = K_FORM / RHOC
X_SCALE = R_E * R_E / (4.0 * ALPHA)  # 50 s for this case


def kernel(s):
    """F(s) = E1(r_e^2/(4 alpha s))/(4 pi k); F(0) = 0 exactly."""
    if s <= 0.0:
        return 0.0
    return E1(X_SCALE / s) / (4.0 * math.pi * K_FORM)


r_semi, semi_frames = series(SEMI, "semi")
r_fine, fine_frames = series(SEMI_FINE, "semiFine")
_, plain_fine = series(PLAIN_FINE, "plainFine")
r_geo, geo_frames = series(GEO, "geo")
r_top, top_frames = series(TOPLEVEL, "toplevel")

# (1) the four printed E1 probes vs an independent quadrature
probe_frame = next(fr for fr in semi_frames if fr["e1"] is not None)
for x, printed in sorted(probe_frame["e1"].items()):
    exact = E1(x)
    check(abs(printed - exact) <= 1e-9 * exact,
          f"E1({x}) printed {printed} != independent quadrature {exact}")

# (2) zero-declaration discipline: no rise, and a foreign scope is ignored
check(all(fr["farRise"] == 0.0 for fr in frames),
      "zero-declaration run developed a non-zero far-field rise")
check(r_top.stdout == proc.stdout,
      "semiInfinite declared in thermal{} changed the run: the declaration "
      "scoped to formation{} only")
check(all(fr["formation"] == tf["formation"]
          for fr, tf in zip(top_frames, frames)),
      "toplevel twin's formation rates differ from the base")

# (3) Duhamel: rebuild every rise from the PRINTED per-step far-field flux
q = {fr["t"]: fr["formation"][0] / LENGTH for fr in fine_frames}
for n in range(1, 31):
    rise = sum((q.get(m, 0.0) - q.get(m - 1, 0.0)) * kernel(n - m)
               for m in range(1, n))
    got = fine_frames[n]["farRise"]
    check(abs(got - rise) <= 1e-10 + 1e-8 * abs(rise),
          f"t={n} farRise {got} != Duhamel reconstruction {rise}")
# the first few steps see an effectively impermeable outer field (x >= 10)
for t in range(1, 6):
    check(abs(fine_frames[t]["farRise"]) < 1e-6,
          f"t={t} line-source rise {fine_frames[t]['farRise']} is not small")
    a, b = fine_frames[t]["formation"][0], plain_fine[t]["formation"][0]
    check(abs(a - b) <= 1e-9 + 1e-6 * abs(b),
          f"t={t} semi-infinite qFar {a} != Dirichlet twin {b} early on")

# (4) the outer field is really unbounded: no steady state, no fixed T_e
rises = [fr["farRise"] for fr in semi_frames if fr["t"] > 0.0]
check(all(b > a for a, b in zip(rises, rises[1:])),
      "line-source far-field rise is not strictly increasing")
check(semi_frames[-1]["farRise"] > 10.0,
      f"late rise {semi_frames[-1]['farRise']} K is too small to be unbounded")
check(semi_frames[-1]["formation"][0] < 0.9 * last["formation"][0],
      "semi-infinite late qFar did not fall below the Dirichlet steady value")
for tag, fr in (("first", next(f for f in semi_frames if f["t"] > 0)),
                ("last", semi_frames[-1])):
    for i in range(SECTIONS):
        check(abs(fr["form"][i][M] - (T_E + fr["farRise"])) < 1e-9,
              f"{tag} frame section {i} far node {fr['form'][i][M]} != "
              f"T_e + rise {T_E + fr['farRise']}")

# (5) energy balance still closes exactly on the semi-infinite path
for tag, fs in (("semi", semi_frames), ("semiFine", fine_frames)):
    for fr in fs:
        if fr["t"] == 0.0:
            continue
        qf, qin, qs = fr["formation"]
        check(abs(qin - (qs + qf)) <= 1e-9 * max(abs(qin), 1.0),
              f"{tag} t={fr['t']} balance broken: intake {qin} != "
              f"storage {qs} + far {qf}")

# (6) geothermal far field: the rise is global, the depth profile survives
glast = geo_frames[-1]
geo_sin = math.sin(math.radians(30.0))


def te_geo(i):
    drop = (i + 0.5) * DX * geo_sin
    return 280.0 + 0.03 * (1000.0 - drop)


for i in range(SECTIONS):
    check(abs(glast["form"][i][M] - (te_geo(i) + glast["farRise"])) < 1e-9,
          f"geothermal section {i} far node != Te({i}) + shared rise")
check(abs((glast["form"][0][M] - glast["form"][SECTIONS - 1][M])
          - 0.03 * (SECTIONS - 1) * DX * geo_sin) < 1e-9,
      "geothermal far-field spread is not the declared depth gradient")
check(glast["form"][0][M] > glast["form"][SECTIONS - 1][M],
      "geothermal far field is not warmer at the deeper inlet")

# (7) loud refusals on the semi-infinite declaration face
semi_refusals = [
    (original.replace("shells 10;", "shells 10;\n    semiInfinite maybe;"),
     "must be 'yes' or 'no'"),
    (re.sub(r"thermal\s*\{[^}]*\}", "", SEMI), "needs a thermal"),
    (SEMI.replace("conductivity 2;", "conductivity 0;"),
     "positive conductivity"),
    (SEMI.replace("outerRadius 0.2;", "outerRadius 0.005;"),
     "exceed the bore"),
    (SEMI.replace("rhoC 10000;", "rhoC 0;"), "positive rhoC"),
    (SEMI.replace("shells 10;", "shells 0;"), "at least one"),
]
for text, needle in semi_refusals:
    rr = run(text)
    check(rr.returncode != 0 and needle in (rr.stdout + rr.stderr),
          f"semi refusal '{needle}' not raised loudly (rc={rr.returncode})")

if fails:
    for m in fails:
        print("FAIL:", m)
    raise SystemExit(f"formation: {len(fails)} of {ncrit} checks FAILED")
print(f"formation: all {ncrit} checks passed "
      f"(T_out {Ts[-1]:.6f} K, r {R_REC:.6f}, g_eff {G_EFF:.6f} W/K, "
      f"R_form {R_FORM:.6f} K/W, steady qFar {last['formation'][0]:.4f} W, "
      f"{nbal} Dirichlet transient balances; 5.84 semi-infinite: "
      f"rise {semi_frames[-1]['farRise']:.6f} K, "
      f"late qFar {semi_frames[-1]['formation'][0]:.4f} W, "
      f"{len(fine_frames) - 1} Duhamel steps, 4 E1 probes)")