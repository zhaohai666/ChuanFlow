#!/usr/bin/env python3
"""5.188 / G6b: the steady thermal profile with a LAYERED wall.

The claim has three parts, all re-derivable from the case text alone:

  1. the layered stack prices the per-unit-length coefficient by the series
     conductance of the cylindrical shells (plus the optional outer film),
       U = 1 / ( sum_i ln(r_{i+1}/r_i)/(2 pi k_i) + 1/(2 pi r_out h_o) ),
     and that is the number the march actually used (`thermal <pipe> u`),
  2. the per-segment P/T profile the report prints is the march itself:
     T[k+1] = T_w + (T[k] + jt (p[k+1]-p[k]) - T_w) exp(-U ds/(m c_p)),
     with T[0] = the fluid inlet temperature, so the printed outlet T on the
     pipe line is the last segment's exit;
  3. declaring the wall changes nothing hydraulically: an isothermal twin of
     the same network is byte-identical once every `thermal` line and every
     ` T <num>` token is stripped, and prints no `thermal` line at all.

Refusals: a layer list whose thickness/conductivity lengths disagree, a
non-positive thickness, a non-positive conductivity, a missing innerRadius, a
negative outer film, and `insulation` together with a declared
`wallHeatTransfer`.

usage: check_steady_thermal_profile.py <steady-exe> <case.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) < 3:
    raise SystemExit(__doc__)

STEADY, CASE = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
CASE_DIR = os.path.dirname(CASE)
TWIN = os.path.join(CASE_DIR, "caseIsothermal.dict")

fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def run(path):
    p = subprocess.run([STEADY, path], capture_output=True, text=True,
                       cwd=os.path.dirname(path) or ".")
    return p.returncode, p.stdout + p.stderr


def strip_comments(text):
    return re.sub(r"//[^\n]*", "", text)


CASE_TEXT = strip_comments(open(CASE).read())


def block_of(name, text=CASE_TEXT):
    """The `name { ... }` body, brace-balanced."""
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", text)
    if not m:
        return None
    i = m.end()
    depth = 1
    while i < len(text) and depth:
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
        i += 1
    return text[m.end():i - 1]


def scalar(body, key, default=None):
    m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)\s*;", body)
    return float(m.group(1)) if m else default


def numlist(body, key):
    m = re.search(r"\b" + key + r"\s*\(([^)]*)\)", body)
    return [float(x) for x in m.group(1).split()] if m else []


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


pipes_body = block_of("pipes")
p1_body = block_of("P1", pipes_body)
p2_body = block_of("P2", pipes_body)
ins_body = block_of("insulation", p1_body)

# ---- the declaration, re-read -------------------------------------------
P1_LEN = scalar(p1_body, "length")
P1_JT = scalar(p1_body, "jouleThomson", 0.0)
P1_CP = scalar(p1_body, "heatCapacity", 0.0)
P1_TW = scalar(p1_body, "wallTemperature")
P2_LEN = scalar(p2_body, "length")
P2_JT = scalar(p2_body, "jouleThomson", 0.0)
INNER = scalar(ins_body, "innerRadius")
THICK = numlist(ins_body, "thickness")
COND = numlist(ins_body, "conductivity")
FILM = scalar(ins_body, "outerCoefficient", 0.0)
SEGMENTS = 24

# the coefficient the stack prices, independently of the kernel
r = INNER
series = 0.0
for d, k in zip(THICK, COND):
    r_next = r + d
    series += math.log(r_next / r) / (2.0 * math.pi * k)
    r = r_next
ROUT = r
if FILM > 0:
    series += 1.0 / (2.0 * math.pi * r * FILM)
U_LAYERS = 1.0 / series
check("the stack's own series conductance",
      abs(U_LAYERS - 0.299793859276) < 1e-9,
      f"U = {U_LAYERS:.12f} W/(m K), r_out = {ROUT:.6f} m")

RC, OUT = run(CASE)
check("case runs", RC == 0, f"rc={RC}")

PIPE = {}
for line in OUT.splitlines():
    m = re.match(r"\[chflow\] pipe (\S+) Q (?:.* )?m (\S+)", line)
    if m:
        PIPE[m.group(1)] = float(m.group(2))

THERM = {}
POINTS = {}
for line in OUT.splitlines():
    m = re.match(r"\[chflow\] thermal (\S+) pt (\d+) x (\S+) p (\S+) T (\S+)",
                 line)
    if m:
        POINTS.setdefault(m.group(1), []).append(
            (int(m.group(2)), float(m.group(3)), float(m.group(4)),
             float(m.group(5))))
        continue
    m = re.match(r"\[chflow\] thermal (\S+) u (\S+) rout (\S+) pts (\d+) Tout (\S+)",
                 line)
    if m:
        THERM[m.group(1)] = dict(u=float(m.group(2)), rout=float(m.group(3)),
                                 pts=int(m.group(4)), tout=float(m.group(5)))

check("the two thermal pipes report a profile",
      sorted(THERM) == ["P1", "P2"], f"{sorted(THERM)}")
check("the undeclared pipe reports no thermal line",
      "P3" not in THERM and "P3" not in POINTS)
check("one point per march segment plus the exit",
      all(THERM[k]["pts"] == SEGMENTS
          and len(POINTS[k]) == SEGMENTS + 1 for k in THERM))

# ---- 1. the number the march used is the stack's ------------------------
check("layered U reaches the march", rel(THERM["P1"]["u"], U_LAYERS) < 1e-12,
      f"printed {THERM['P1']['u']:.12f} vs stack {U_LAYERS:.12f}")
check("the outer radius is the stack's own",
      rel(THERM["P1"]["rout"], ROUT) < 1e-12,
      f"{THERM['P1']['rout']} vs {ROUT}")
CHECK_OK = True

# ---- 2. the printed profile IS the march ---------------------------------
NODE_P = {}
for line in OUT.splitlines():
    m = re.match(r"\[chflow\] node (\S+) p (\S+)", line)
    if m:
        NODE_P[m.group(1)] = float(m.group(2))

P1 = POINTS["P1"]
P2 = POINTS["P2"]
m1 = PIPE["P1"]
m2 = PIPE["P2"]
ds1 = P1_LEN / SEGMENTS
ds2 = P2_LEN / SEGMENTS


def wall(t, tw, u, ds, mass_rate, cp):
    if u <= 0 or mass_rate <= 0:
        return t
    return tw + (t - tw) * math.exp(-u * ds / (mass_rate * cp))


# P1: the wall law on every segment, then the exit on the march's own p
t = P1[0][3]
rebuilt1 = [t]
for k in range(SEGMENTS - 1):
    t = t + P1_JT * (P1[k + 1][2] - P1[k][2])
    t = wall(t, P1_TW, U_LAYERS, ds1, m1, P1_CP)
    rebuilt1.append(t)
worst1 = max(rel(a[3], b) for a, b in zip(P1, rebuilt1))
check("P1 profile follows the wall law segment by segment",
      worst1 < 1e-9,
      f"worst relative deviation {worst1:.3e} (m = {m1:.6f} kg/s, "
      f"ds = {ds1} m, U = {U_LAYERS:.9f})")

t_exit = t + P1_JT * (P1[SEGMENTS][2] - P1[SEGMENTS - 1][2])
t_exit = wall(t_exit, P1_TW, U_LAYERS, ds1, m1, P1_CP)
check("P1 outlet T is the march's own last exit",
      rel(P1[SEGMENTS][3], t_exit) < 1e-12 and
      rel(THERM["P1"]["tout"], t_exit) < 1e-12,
      f"printed {P1[SEGMENTS][3]:.9f} vs rebuilt {t_exit:.9f} "
      f"(last dp {P1[SEGMENTS][2] - P1[SEGMENTS - 1][2]:.3f} Pa)")
check("P1 exit pressure is the last step's own",
      P1[SEGMENTS][2] < P1[SEGMENTS - 1][2],
      f"p {P1[SEGMENTS - 1][2]:.3f} -> {P1[SEGMENTS][2]:.3f} Pa")

# P2: pure Joule-Thomson on the printed pressures (no wall to relax toward)
rebuilt2 = [P2[0][3] + P2_JT * (p[2] - P2[0][2]) for p in P2]
worst2 = max(rel(a[3], b) for a, b in zip(P2, rebuilt2))
check("P2 is the pure J-T telescope on the printed pressures",
      worst2 < 1e-11, f"worst relative deviation {worst2:.3e}")
check("P2 prices no wall", THERM["P2"]["u"] == 0.0)

# the profile is monotone in x and lands on the pipe length
for name, pts, length in (("P1", P1, P1_LEN), ("P2", P2, P2_LEN)):
    dx = length / SEGMENTS
    check(f"{name} profile grid is x = k ds",
          all(rel(p[1], k * dx) < 1e-12 for k, p in enumerate(pts))
          and all(p[0] == k for k, p in enumerate(pts)))

# ---- 3. the twin ----------------------------------------------------------
if os.path.exists(TWIN):
    RT, OUT_T = run(TWIN)
    check("twin runs", RT == 0, f"rc={RT}")
    check("twin prints no thermal line", "thermal " not in OUT_T)

    def strip(text):
        text = re.sub(r"\[chflow\] thermal .*\n", "", text)
        text = re.sub(r" T [-+0-9.eE]+", "", text)
        return text

    check("hydraulics are the twin's, once the thermal readout is stripped",
          strip(OUT) == strip(OUT_T),
          "the pipes, nodes and solve lines are byte-identical")
else:
    check("twin present", False, f"missing {TWIN}")

# ---- 4. refusals ----------------------------------------------------------
base_p1 = p1_body


def refuse(name, ins_replacement=None, extra="", needle=""):
    """Rewrite P1's insulation block (or append a key) and expect a refusal."""
    body = base_p1
    if ins_replacement is not None:
        body = re.sub(r"insulation\s*\{[^}]*\}", ins_replacement, body,
                      flags=re.S)
    body = body + extra
    mutated = CASE_TEXT.replace(base_p1, body)
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "caseReject.dict")
        open(path, "w").write(mutated)
        rc, out = run(path)
    check(name, rc != 0 and needle in out,
          out.strip().splitlines()[-1][:110] if out.strip() else "no output")


refuse("refuses mismatched layer lists",
       "insulation { innerRadius 0.04; thickness (0.03 0.02); "
       "conductivity (0.05); }",
       needle="one conductivity per thickness entry")
refuse("refuses a zero shell",
       "insulation { innerRadius 0.04; thickness (0.0 0.02); "
       "conductivity (0.05 0.03); }",
       needle="thickness must be positive")
refuse("refuses a non-positive conductivity",
       "insulation { innerRadius 0.04; thickness (0.03 0.02); "
       "conductivity (-0.05 0.03); }",
       needle="conductivity must be positive")
refuse("refuses a missing innerRadius",
       "insulation { thickness (0.03); conductivity (0.05); }",
       needle="needs a positive innerRadius")
refuse("refuses a negative outer film",
       "insulation { innerRadius 0.04; thickness (0.03); "
       "conductivity (0.05); outerCoefficient -8; }",
       needle="negative insulation outerCoefficient")
refuse("refuses insulation and a declared wall coefficient",
       extra=" wallHeatTransfer 2.0;",
       needle="declares both insulation and a wallHeatTransfer")

print()
if fails:
    print(f"{len(fails)} ASSERTION(S) FAILED: {fails}")
    sys.exit(1)
print("all checks passed")
