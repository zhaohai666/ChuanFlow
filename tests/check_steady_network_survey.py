#!/usr/bin/env python3
"""Verification of the WELL-TRAJECTORY SURVEY (5.193, G9) on the steady network.

The increment: a pipe could state one slope for its whole arc, so a build
section and a horizontal lateral were both priced as one straight line. A
`survey { md (...); inclination (...); }` states the inclination against
measured depth, the elevation becomes the survey's own integral, and the march
prices gravity segment by segment.

Criterion makeup. Every number below is re-derived from the case text and the
fluid table alone -- the checker never trusts a solver internal. The five
independent families are:

  S1. THE ARC IS THE SURVEY'S. The printed `dz` is the closed-form integral
      `sum_i (cos th_i - cos th_{i+1})/(th_{i+1} - th_i) * ds_i` of the
      inclination list in the case, to 1e-12; the printed `sinmin`/`sinmax`
      are the min/max of the 24 printed segment samples; every `surveypt`
      sits at a segment MIDPOINT and carries sin(theta(md)) interpolated
      linearly in md. (A survey read as a straight mean would fail here.)

  S2. THE MARCH IS PRICED ON THOSE SEGMENTS. For each pvt pipe the checker
      marches the case's own 24 segments -- per-segment sin from S1, the
      clamped linear rho(p)/mu(p) table, Churchill at the local Reynolds
      number, mass rate fixed at the upstream density -- and bisects the
      rate that lands the outlet on the fixed low-side pressure. It must
      reproduce the printed Q to 1e-9 relative. This is the criterion that
      makes the segmentation, not the declaration, falsifiable.

  S3. THE INVOLUTION IS ACCEPTED AND THE ASYMMETRY IS NOT. P1 walks the arc
      70 -> 50 -> 20 -> 0 degrees and P2 walks 0 -> 20 -> 50 -> 70: the
      closed-form integral is symmetric in its endpoints, so their `dz` must
      agree to 1e-15 while their rates must DIFFER -- the steep part sits at
      the dense end for P1 and at the sparse end for P2. P5 is P1's arc
      walked downhill (negative inclinations): its `dz` must be exactly
      -dz(P1) and its rate must exceed P2's, because gravity then assists.

  S4. THE SINGLE-SLOPE IDENTITY. P3 states one inclination as a survey and P4
      declares the survey's own integral as `elevation` with no survey at
      all: Q, velocity, Reynolds and friction must agree to 1e-12. The only
      difference between them is whether the march is handed one sin or 24
      identical ones, so a plumbing error in the segmentation shows here and
      nowhere else. (The same pair exists on the black-oil march as O2/O4,
      where the holdup law reads the inclination a second time.)

  S5. THE GRAMMAR IS REFUSED. Ten malformed declarations are each run and
      must fail with the named message: unequal lists, md not from 0, md
      not to the pipe's length, md not monotone, |inclination| > 90, an
      elevation that disagrees with the integral, an unknown key inside
      survey{}, a single station, a survey on a coupled edge, and a survey
      beside a slugPrediction trajectory.

On the black-oil case the checker verifies the same dz/sinmin/surveypt
structure and the two twin pairs, plus one asymmetry: the constant 12-degree
twin O2 must NOT reproduce the varying arc O3, which is what makes "the arc is
the survey, not its mean" falsifiable on the Beggs-Brill march too.

usage: check_steady_network_survey.py <exe> <case.dict> <caseBlackOil.dict>
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

if len(sys.argv) != 4:
    raise SystemExit(__doc__)

EXE, CASE, CASE_BO = sys.argv[1], sys.argv[2], sys.argv[3]
SEGMENTS = 24
G = 9.80665

fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def run(path, cwd=None):
    p = subprocess.run([EXE, path], capture_output=True, text=True,
                       cwd=cwd or (os.path.dirname(path) or "."))
    return p.returncode, p.stdout + p.stderr


def strip_comments(text):
    return re.sub(r"//[^\n]*", "", text)


def block_of(name, text):
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


def scalar(body, key, default=None):
    m = re.search(r"\b" + key + r"\s+([-+0-9.eE]+)\s*;", body)
    return float(m.group(1)) if m else default


def numlist(body, key):
    m = re.search(r"\b" + key + r"\s*\(([^)]*)\)", body)
    return [float(x) for x in m.group(1).split()] if m else []


def z_integral(md, incl):
    """closed form of integral sin(theta) ds for a theta linear in md"""
    tot = 0.0
    for i in range(len(md) - 1):
        ds = md[i + 1] - md[i]
        a, b = math.radians(incl[i]), math.radians(incl[i + 1])
        tot += (math.sin(a) * ds if abs(b - a) < 1e-12
                else (math.cos(a) - math.cos(b)) / (b - a) * ds)
    return tot


def seg_sin(md, incl, length, nseg=SEGMENTS):
    """sin(theta) at each segment midpoint, theta piecewise-linear in md"""
    out = []
    ds = length / nseg
    k = 0
    for i in range(nseg):
        mid = (i + 0.5) * ds
        while k + 2 < len(md) and mid >= md[k + 1]:
            k += 1
        f = ((mid - md[k]) / (md[k + 1] - md[k])) if md[k + 1] > md[k] else 0.0
        th = math.radians(incl[k] + f * (incl[k + 1] - incl[k]))
        out.append(math.sin(th))
    return out


# ---------------------------------------------------------------- parsing
def parse_run(text):
    """pipe / survey / surveypt tables out of one app run"""
    pipes, survey, pts = {}, {}, {}
    for line in text.splitlines():
        body = re.sub(r"^\S+chflow\S*\]\s*", "", line.strip()) \
            if line.startswith("[chflow]") else None
        if body is None:
            continue
        tok = body.split()
        if not tok:
            continue
        if tok[0] in ("pipe", "survey", "surveypt", "node"):
            key = tok[0]
            name = tok[1]
            vals = tok[2:]
            if key == "survey":
                if name not in survey:
                    survey[name] = {"dz": float(vals[vals.index("dz") + 1]),
                                    "sinmin": float(vals[vals.index("sinmin") + 1]),
                                    "sinmax": float(vals[vals.index("sinmax") + 1]),
                                    "st": int(vals[vals.index("st") + 1]),
                                    "seg": int(vals[vals.index("seg") + 1])}
                continue
            if key == "surveypt":
                # surveypt P1 k 0 md 25 sin 0.932
                pts.setdefault(name, []).append(
                    (int(tok[3]), float(tok[5]), float(tok[7])))
                continue
            if key == "node":
                pipes.setdefault("_nodes", {})[name] = float(tok[3])
                continue
            d = pipes.setdefault(name, {})
            for i in range(0, len(vals) - 1, 2):
                try:
                    d[vals[i]] = float(vals[i + 1])
                except ValueError:
                    d[vals[i]] = float("nan")
    return pipes, survey, pts


CASE_TEXT = strip_comments(open(CASE).read())
RC, OUT = run(CASE)
check("the pvt case runs", RC == 0, f"rc={RC}")
if RC != 0:
    print(OUT)
    raise SystemExit(1)

PIPES_BODY = block_of("pipes", CASE_TEXT)
FLUID = block_of("pvt", block_of("fluid", CASE_TEXT))
P_TABLE = numlist(FLUID, "pressure")
RHO_TABLE = numlist(FLUID, "rho")
MU_TABLE = numlist(FLUID, "mu")


def lookup(p):
    if p <= P_TABLE[0]:
        return RHO_TABLE[0], MU_TABLE[0]
    if p >= P_TABLE[-1]:
        return RHO_TABLE[-1], MU_TABLE[-1]
    for i in range(len(P_TABLE) - 1):
        if P_TABLE[i] <= p <= P_TABLE[i + 1]:
            f = (p - P_TABLE[i]) / (P_TABLE[i + 1] - P_TABLE[i])
            return (RHO_TABLE[i] + f * (RHO_TABLE[i + 1] - RHO_TABLE[i]),
                    MU_TABLE[i] + f * (MU_TABLE[i + 1] - MU_TABLE[i]))
    raise AssertionError


def churchill(re, eps_rel):
    re = max(re, 1e-12)
    a = (2.457 * math.log(1.0 / ((7.0 / re) ** 0.9 + 0.27 * eps_rel))) ** 16.0
    b = (37530.0 / re) ** 16.0
    return 8.0 * ((8.0 / re) ** 12.0 + (a + b) ** -1.5) ** (1.0 / 12.0)


def march(length, diameter, roughness, q, p_up, sin_list, nseg=SEGMENTS):
    area = 0.25 * math.pi * diameter ** 2
    mass = lookup(p_up)[0] * q
    ds = length / nseg
    p = p_up
    floor = max(0.02 * p_up, 1000.0)
    for i in range(nseg):
        r, m = lookup(p)
        v = mass / (r * area)
        f = churchill(r * v * diameter / m, roughness / diameter)
        dp = (r * G * sin_list[i] * ds
              + f * ds / diameter * r * v * v / 2.0)
        p -= dp
        if not p > floor:
            return floor
    return p


def rate_for_drop(length, diameter, roughness, p_up, p_low, sin_list):
    def resid(q):
        return march(length, diameter, roughness, q, p_up, sin_list) - p_low
    lo, hi = 1e-6, 10.0
    f_lo, f_hi = resid(lo), resid(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f_mid = resid(mid)
        if f_mid > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def pipe_decl(name):
    body = block_of(name, PIPES_BODY)
    sv = block_of("survey", body)
    return dict(length=scalar(body, "length"),
                diameter=scalar(body, "diameter"),
                roughness=scalar(body, "roughness", 0.0),
                elevation=scalar(body, "elevation", None),
                md=numlist(sv, "md") if sv else [],
                incl=numlist(sv, "inclination") if sv else [])


PIPE, SURVEY, PT = parse_run(OUT)

# ------------------------------------------------------------ S1: the arc
SURVEYED = ["P1", "P2", "P3", "P5"]
check("the surveyed pipes report a survey line, the twin does not",
      sorted(SURVEY) == sorted(SURVEYED) and "P4" not in PT
      and "P4" not in SURVEY,
      f"{sorted(SURVEY)}")

for name in SURVEYED:
    d = pipe_decl(name)
    want = z_integral(d["md"], d["incl"])
    check(f"{name}: dz is the survey's closed-form integral",
          rel(SURVEY[name]["dz"], want) < 1e-12,
          f"printed {SURVEY[name]['dz']:.12f} vs integral {want:.12f} m")
    check(f"{name}: one sample per march segment",
          SURVEY[name]["seg"] == SEGMENTS
          and len(PT[name]) == SEGMENTS, f"seg={SURVEY[name]['seg']}")
    check(f"{name}: stations are the declaration's",
          SURVEY[name]["st"] == len(d["md"]), f"st={SURVEY[name]['st']}")
    want_sin = seg_sin(d["md"], d["incl"], d["length"])
    worst = max(rel(want_sin[k], PT[name][k][2])
                for k in range(SEGMENTS))
    # the printed tokens carry 12 significant digits, so this bound is the
    # print's own resolution, not a physics tolerance
    check(f"{name}: every segment sin is sin(theta(md)) at the midpoint",
          worst < 1e-10, f"worst rel {worst:.3e}")
    md_worst = max(rel((k + 0.5) * d["length"] / SEGMENTS, PT[name][k][1])
                   for k in range(SEGMENTS))
    check(f"{name}: the samples sit at segment midpoints",
          md_worst < 1e-12, f"worst rel {md_worst:.3e}")
    lo = min(s for _, _, s in PT[name])
    hi = max(s for _, _, s in PT[name])
    check(f"{name}: sinmin/sinmax are the samples' own bounds",
          rel(SURVEY[name]["sinmin"], lo) < 1e-12
          and rel(SURVEY[name]["sinmax"], hi) < 1e-12,
          f"{SURVEY[name]['sinmin']:.6f}/{SURVEY[name]['sinmax']:.6f} vs "
          f"{lo:.6f}/{hi:.6f}")

# ----------------------------------------------------- S2: the march's own
for name in SURVEYED + ["P4"]:
    d = pipe_decl(name)
    sins = seg_sin(d["md"], d["incl"], d["length"]) if d["md"] else \
        [d["elevation"] / d["length"]] * SEGMENTS
    q = rate_for_drop(d["length"], d["diameter"], d["roughness"],
                      3.0e6, 1.0e6, sins)
    check(f"{name}: the restated 24-segment march lands on the printed rate",
          rel(q, PIPE[name]["Q"]) < 1e-9,
          f"march {q:.12f} vs printed {PIPE[name]['Q']:.12f} m3/s"
          f" (rel {rel(q, PIPE[name]['Q']):.2e})")

# ------------------------------- S3: same drop, different bill, and downhill
check("P1 and P2 state the same vertical drop",
      rel(SURVEY["P1"]["dz"], SURVEY["P2"]["dz"]) < 1e-15,
      f"{SURVEY['P1']['dz']:.12f} vs {SURVEY['P2']['dz']:.12f} m")
check("P1 (steep at the dense end) is slower than P2 (steep at the sparse end)",
      PIPE["P1"]["Q"] < PIPE["P2"]["Q"],
      f"{PIPE['P1']['Q']:.10f} < {PIPE['P2']['Q']:.10f} m3/s")
margin = (PIPE["P2"]["Q"] - PIPE["P1"]["Q"]) / PIPE["P2"]["Q"]
check("the reversal costs ~2.4% of the rate, not a rounding",
      0.015 < margin < 0.035, f"{margin * 100:.3f}%")
check("P5's arc is exactly P1's, negated",
      rel(SURVEY["P5"]["dz"], -SURVEY["P1"]["dz"]) < 1e-15,
      f"{SURVEY['P5']['dz']:.12f} m")
check("downhill beats both uphill arcs",
      PIPE["P5"]["Q"] > PIPE["P2"]["Q"] > PIPE["P1"]["Q"],
      f"{PIPE['P5']['Q']:.10f} > {PIPE['P2']['Q']:.10f} > "
      f"{PIPE['P1']['Q']:.10f}")

# -------------------------------------- S4: one slope, two ways of stating it
for token in ("Q", "v", "Re", "f", "dp"):
    check(f"P3 (survey) and P4 (elevation) agree on {token}",
          rel(PIPE["P3"][token], PIPE["P4"][token]) < 1e-12,
          f"{PIPE['P3'][token]:.12g} vs {PIPE['P4'][token]:.12g}")
check("the twin's survey integral is the elevation it declared",
      rel(SURVEY["P3"]["dz"], pipe_decl("P4")["elevation"]) < 1e-12,
      f"{SURVEY['P3']['dz']:.12f} vs {pipe_decl('P4')['elevation']:.12f} m")

# ------------------------------------------------------------ S5: refusals
def refuse(label, mutate, needle):
    text = mutate(CASE_TEXT) if mutate else CASE_TEXT
    d = tempfile.mkdtemp()
    try:
        path = os.path.join(d, "case.dict")
        open(path, "w").write(text)
        rc, out = run(path, cwd=d)
        ok = rc != 0 and needle in out
        check(f"refused: {label}", ok,
              f"rc={rc}" + ("" if ok else " :: " + out.strip().splitlines()[-1]))
    finally:
        shutil.rmtree(d)


def sub1(pattern, repl, text):
    new, n = re.subn(pattern, repl, text, count=1)
    if n != 1:
        raise AssertionError(f"pattern not found: {pattern}")
    return new


P1_SURVEY = 'md ( 0 400 800 1200 );\n            inclination ( 70 50 20 0 );'

refuse("unequal md and inclination lists",
       lambda t: sub1(r"inclination \( 70 50 20 0 \);",
                      "inclination ( 70 50 20 );", t),
       "one angle per station")
refuse("md not starting at 0",
       lambda t: sub1(r"md \( 0 400 800 1200 \);", "md ( 5 400 800 1200 );", t),
       "starts at 0")
refuse("md not ending at the pipe's length",
       lambda t: sub1(r"md \( 0 400 800 1200 \);", "md ( 0 400 800 1100 );", t),
       "the survey lays the pipe's own arc out")
refuse("md not monotone",
       lambda t: sub1(r"md \( 0 400 800 1200 \);", "md ( 0 800 400 1200 );", t),
       "measured depth only grows")
refuse("inclination outside [-90, 90]",
       lambda t: sub1(r"inclination \( 70 50 20 0 \);",
                      "inclination ( 70 50 20 120 );", t),
       "outside [-90, 90]")
refuse("a declared elevation that disagrees with the integral",
       lambda t: sub1(r"roughness 4\.6e-5;\n        survey",
                      "roughness 4.6e-5; elevation 700;\n        survey", t),
       "the survey is the arc")
refuse("an unknown key inside survey{}",
       lambda t: sub1(r"inclination \( 70 50 20 0 \);",
                      "inclination ( 70 50 20 0 ); azimuth ( 0 0 0 0 );", t),
       "a survey is md and inclination")
refuse("a single station",
       lambda t: sub1(r"md \( 0 400 800 1200 \);\n            inclination "
                      r"\( 70 50 20 0 \);",
                      "md ( 0 1200 );\n            inclination ( 30 30 );", t)
       if False else sub1(r"inclination \( 70 50 20 0 \);", "inclination ( 70 );", t),
       "at least two stations")
refuse("a survey on a coupled edge (the coupled-element rule fires first, so "
       "the declaration never reaches the march)",
       lambda t: sub1(r"roughness 4\.6e-5;\n        survey",
                      "roughness 4.6e-5; coupledResistance 1.0e6;\n        survey", t),
       "coupled element")

# ---------------------------------------------- the black-oil march, same test
BO_TEXT = strip_comments(open(CASE_BO).read())
RC2, OUT2 = run(CASE_BO)
check("the black-oil case runs", RC2 == 0, f"rc={RC2}")
if RC2 != 0:
    print(OUT2)
    raise SystemExit(1)

P2, S2, T2 = parse_run(OUT2)
BO_PIPES = block_of("pipes", BO_TEXT)

check("O1, O2 and O3 report a survey, O4 does not",
      sorted(S2) == ["O1", "O2", "O3"] and "O4" not in T2, f"{sorted(S2)}")

for name in ("O1", "O3"):
    body = block_of(name, BO_PIPES)
    sv = block_of("survey", body)
    md, incl = numlist(sv, "md"), numlist(sv, "inclination")
    length = scalar(body, "length")
    want = z_integral(md, incl)
    check(f"{name}: the closed-form integral of the declared arc",
          rel(S2[name]["dz"], want) < 1e-12,
          f"printed {S2[name]['dz']:.12f} vs integral {want:.12f} m")
    want_sin = seg_sin(md, incl, length)
    worst = max(rel(want_sin[k], T2[name][k][2]) for k in range(SEGMENTS))
    check(f"{name}: segment sins are the interpolated inclinations",
          worst < 1e-10, f"worst rel {worst:.3e}")
    diffs = [T2[name][k + 1][2] - T2[name][k][2] for k in range(SEGMENTS - 1)]
    signed = 1 if all(d > 0 for d in diffs) else -1
    check(f"{name}: sin(theta) moves monotonically along the arc",
          all(d > 0 for d in diffs) or all(d < 0 for d in diffs),
          f"{'rising' if signed > 0 else 'falling'}")
    check(f"{name}: exactly one sample per segment",
          len(T2[name]) == SEGMENTS, f"{len(T2[name])}")

check("O1 and O3 state the same vertical drop",
      rel(S2["O1"]["dz"], S2["O3"]["dz"]) < 1e-15,
      f"{S2['O1']['dz']:.12f} m")
check("the reversed arc is faster (holdup sees the inclination too)",
      P2["O1"]["Q"] < P2["O3"]["Q"],
      f"{P2['O1']['Q']:.10f} < {P2['O3']['Q']:.10f} m3/s")

for token in ("Q", "v", "HL", "lambda", "Nfr"):
    check(f"O2 (survey) and O4 (elevation) agree on {token}",
          rel(P2["O2"][token], P2["O4"][token]) < 1e-12,
          f"{P2['O2'][token]:.12g} vs {P2['O4'][token]:.12g}")

mw = abs(P2["O2"]["Q"] - P2["O3"]["Q"]) / P2["O3"]["Q"]
check("the mean-inclination twin does NOT stand in for the varying arc",
      mw > 5e-3, f"|Q(O2) - Q(O3)|/Q(O3) = {mw:.4f}")

refuse_bo = tempfile.mkdtemp()
try:
    text = sub1(r"roughness 4\.6e-5;\n        survey",
                "roughness 4.6e-5;\n        slugPrediction { trajectory "
                "{ s1 { length 600; angleDeg 0; } s2 { length 600; angleDeg 0; } } }\n        survey",
                BO_TEXT)
    path = os.path.join(refuse_bo, "case.dict")
    open(path, "w").write(text)
    rc, out = run(path, cwd=refuse_bo)
    check("refused: a survey beside a slugPrediction trajectory",
          rc != 0 and "two statements of one arc" in out,
          f"rc={rc}" + ("" if rc != 0 else ""))
finally:
    shutil.rmtree(refuse_bo)

print()
if fails:
    print(f"{len(fails)} criterion(s) failed")
    raise SystemExit(1)
print("all criteria passed")
