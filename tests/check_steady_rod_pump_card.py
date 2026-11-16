#!/usr/bin/env python3
"""Verification of the subsurface rod pump and its card (5.90, D1 increment 4).

A rod-pump row is the one row in this kernel whose rate is not an answer. It
takes its volume from the machine -- barrel area, stroke, speed -- and hands the
network back the pressure its two ports differ by. The card is what that
exchange looks like drawn on an axis: four statements, two of them statics and
one of them Hooke's law,

  rod weight   Wr  = (rho_r - rho_l) g A_r L        (weight in air minus buoyancy)
  fluid load   Fo  = A_p (p_above - p_below)        (the two solved node pressures)
  stretch      df   = Fo L / (A_r E)                (Hooke, on the load's own part)
  stroke       Seff = S - df, and the card is the parallelogram
                 (0,Wr) (df,Wr+Fo) (S,Wr+Fo) (S-df,Wr), area Fo Seff.

The checker reads no solver internals. It takes the pump's seven numbers, the
two column bores and the four boundary pressures out of the case dictionary, the
solved node pressures and every card field out of the 12-digit report, and
restates all four in closed form -- including the differential the plunger sees,
which for a well arm is the classic sum

  dp = (p_wh - p_cas) + rho g (L - h) + 128 mu L_t Q_t/(pi d_t^4)
                          - 128 mu L_a Q_a/(pi d_a^4)

so the level, the casing charge, the tubing's weight, the tubing's friction and
the annulus's friction are all priced twice over: once by the network's march and
once by the card's own load, and the two must agree.

What the arms are then made to prove against each other, from the printed rows
alone:

  the card is the load, not the pump  : RP_BUOY and RP_DEPR are the same well
                          with a string of a different density, so Fo, stretch,
                          Seff, EV, rate, work and power are the SAME numbers and
                          only the rod-weight datum moves, by exactly
                          (rho_1 - rho_2) g A_r L;
  a stiffer string lifts more         : RP_GAUGE's 7/8 in rod carries a heavier
                          string and a longer effective stroke on the same load,
                          and A_r df / Fo is the same compliance, L/E, for both;
  speed is not a load                 : RP_FAST reaches 12/8 of the rate while
                          its efficiency moves by the friction the extra rate
                          costs, not by anything in the force balance;
  a flooded annulus unloads the plunger: RP_FLOOD's whole hydrostatic column is
                          taken off underneath, so its card is a sliver and its
                          upthrust is the largest of the set;
  a charged annulus flats it out      : RP_IDLE runs its full swept volume at
                          zero load, zero work and zero power, with dp the wrong
                          way round -- and the volume it delivers is NOT paid for
                          by the rods;
  a pump can be stroked and stalled   : RP_STALL's string stretches past its
                          stroke, so it delivers 0.0 and every number still says
                          the load is there;
  an imposed rate finds its own load  : RP_RIG closes A_p Seff N = 2.8e-4, so the
                          stretch the solve needed is S - Q/(A_p N) and the
                          differential is that stretch times A_r E/(A_p L);
  the card absorbs a sub-plunger friction: RP_ANNFRI runs with its whole rate down
                          a narrow annulus, and the 1.125 MPa that costs is in
                          the load on the plunger, not in a correction;
  the identity that ties them together: power = dp x Q on every arm that lifts,
                          and the parallelogram's area Fo Seff is that dp times
                          the volume actually swept;
  and the sign of upthrust is reachable: the flooded well, re-run against a
                          14 MPa header, holds MORE than its frictionless column
                          and reports a negative upthrust -- the ceiling stops
                          being one when the separator outweighs the lift.

usage: check_steady_rod_pump_card.py <steady-exe> <case.dict> \
        [<transient-exe> <twofluid-exe> <nodal-exe>]
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) < 3:
    raise SystemExit(__doc__)


def _exe(arg):
    return os.path.abspath(arg) if arg else None


STEADY, CASE = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
TRANSIENT, TWOFLUID, NODAL = (_exe(sys.argv[i]) if len(sys.argv) > i else None
                              for i in (3, 4, 5))
HERE = os.path.dirname(os.path.abspath(__file__))

G = 9.80665
TIGHT = 1e-10   # a printed number against a closed form of its own declarations
SOLVE = 1e-7    # a printed number against something the Newton solve had to find

fails = []
err = {}


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def near(got, want, tol):
    return rel(got, want) <= tol


def worst(key, val):
    err[key] = max(err.get(key, 0.0), val)


def rows_from(stdout):
    """The report's node pressures and pipe rows, as printed numbers."""
    node_p, rows, raws = {}, {}, {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "node" and len(tok) >= 4:
            node_p[tok[1]] = float(tok[3])
        elif tok[0] == "pipe":
            f, raw = {}, {}
            j = 2
            while j < len(tok) - 1:
                if tok[j] in ("rod", "valve"):
                    f[tok[j]] = tok[j + 1]
                    j += 2
                elif tok[j] in ("SONIC", "CHOKED", "FORCED"):
                    f[tok[j]] = 1.0
                    j += 1
                else:
                    f[tok[j]] = float(tok[j + 1])
                    raw[tok[j]] = tok[j + 1]
                    j += 2
            rows[tok[1]] = f
            raws[tok[1]] = raw
    return node_p, rows, raws


# ---------------------------------------------------------------- run the case
proc = subprocess.run([STEADY, CASE], capture_output=True, text=True)
if proc.returncode != 0:
    raise SystemExit(f"solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")

NODE_P, ROW, RAW = rows_from(proc.stdout)
if not ROW:
    raise SystemExit("no pipe lines parsed")

last = proc.stdout.strip().splitlines()[-1]
check("the network solved", "converged" in last and "NOT converged" not in last,
      last.replace("[chflow] steady solved: ", ""))
iters = float(re.search(r"iterations=(\S+)", last).group(1))
relres = float(re.search(r"residual=(\S+)", last).group(1))
check("and a positive-displacement row did not need the whole budget to say so",
      0 < iters < 200, f"{iters:.0f} Newton iterations, residual {relres:.3g}")

# ---------------------------------------------------------------- the case text
src = open(CASE).read()
text = re.sub(r"//[^\n]*", "", src)


def block(body, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", body)
    if not m:
        raise SystemExit(f"case has no {key}{{}} block")
    start, depth = m.end() - 1, 0
    for j in range(start, len(body)):
        if body[j] == "{":
            depth += 1
        elif body[j] == "}":
            depth -= 1
            if depth == 0:
                return body[start + 1:j]
    raise SystemExit(f"unterminated {key}{{}} block")


def entries(body):
    out = {}
    for m in re.finditer(r"(\w+)\s*\{", body):
        start, depth = m.end() - 1, 0
        for j in range(start, len(body)):
            if body[j] == "{":
                depth += 1
            elif body[j] == "}":
                depth -= 1
                if depth == 0:
                    out[m.group(1)] = body[start + 1:j]
                    break
    return out


def num(body, key, default=None):
    m = re.search(r"\b" + re.escape(key) + r"\s+(-?[0-9.eE+-]+)\s*;", body)
    return float(m.group(1)) if m else default


def has(body, key):
    return re.search(r"\b" + re.escape(key) + r"\s+-?[0-9.eE+-]+\s*;", body) is not None


def word(body, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(\w+)\s*;", body)
    return m.group(1) if m else None


FILL = block(text, "fluid")
RHO, MU = num(FILL, "rho"), num(FILL, "mu")
if not (RHO and MU and RHO > 0 and MU > 0):
    raise SystemExit("case fill malformed")

NODES = entries(block(text, "nodes"))
PIPES = entries(block(text, "pipes"))
SOURCES = entries(block(text, "sources"))

FIXED = {n: num(b, "pressure") for n, b in NODES.items() if has(b, "pressure")}
FREE = sorted(set(NODES) - set(FIXED))
FEED = {}
for name, body in SOURCES.items():
    node = word(body, "node")
    FEED[node] = FEED.get(node, 0.0) + num(body, "rate", 0.0)

GEOM = {}
for name, body in PIPES.items():
    GEOM[name] = dict(
        src=word(body, "from"), dst=word(body, "to"),
        length=num(body, "length", 0.0),
        diameter=num(body, "diameter", 0.0),
        elevation=num(body, "elevation", 0.0),
        ap=math.pi / 4.0 * (num(body, "rodPlungerDiameter", 0.0) or 0.0) ** 2,
        ar=math.pi / 4.0 * (num(body, "rodDiameter", 0.0) or 0.0) ** 2,
        stroke=num(body, "rodStrokeLength", 0.0),
        spm=num(body, "rodStrokesPerMinute", 0.0),
        depth=num(body, "rodLength", 0.0),
        young=num(body, "rodYoungModulus", 0.0),
        rodrho=num(body, "rodDensity", 0.0))
PUMPS = sorted(n for n, g in GEOM.items() if g["ap"] > 0)
MARCH = sorted(n for n, g in GEOM.items() if g["ap"] <= 0)
check("the case is a pump case", len(PUMPS) >= 8 and len(MARCH) >= 8,
      f"{len(PUMPS)} pump rows, {len(MARCH)} marching rows")
check("every pump row is a point element",
      all(GEOM[n]["length"] == 0.0 and GEOM[n]["diameter"] == 0.0
          and GEOM[n]["elevation"] == 0.0 for n in PUMPS),
      "length 0; diameter 0; and no elevation")
check("and no marching row pretends to be one",
      all(has(PIPES[n], "rodPlungerDiameter") is False for n in MARCH))
check("every pump row is in the report", all(n in ROW and "rod" in ROW[n]
                                            for n in PUMPS))

ARMS = ["RP_FLOOD", "RP_DEPR", "RP_GAUGE", "RP_DEPTH", "RP_STALL", "RP_IDLE",
        "RP_BUOY", "RP_FAST", "RP_RIG", "RP_ANNFRI"]
check("all ten arms are here", all(a in PUMPS for a in ARMS),
      " ".join(sorted(set(ARMS) - set(PUMPS))) or "as documented")

# ------------------------------------------------- the card, from the case alone
CARD = {}
for name in PUMPS:
    g = GEOM[name]
    row = ROW[name]
    # dpTotal = p(from) - p(to); the plunger works against the negation of it.
    rise = -row["dp"]
    lift = max(rise, 0.0)
    fo = lift * g["ap"]
    stretch = fo * g["depth"] / (g["ar"] * g["young"])
    seff = max(g["stroke"] - stretch, 0.0)
    ns = g["spm"] / 60.0
    wr = (g["rodrho"] - RHO) * G * g["ar"] * g["depth"]
    state = ("NO-LOAD" if lift <= 0.0 else
             "STALLED" if seff <= 0.0 else "FULL")
    CARD[name] = dict(
        rise=rise, fo=fo, fomax=RHO * G * g["depth"] * g["ap"],
        wr=wr, wair=g["rodrho"] * G * g["ar"] * g["depth"],
        buoy=RHO * G * g["ar"] * g["depth"],
        stretch=stretch, seff=seff, ev=seff / g["stroke"], disp=g["ap"] * seff * ns,
        work=fo * seff, power=fo * seff * ns, state=state, ns=ns,
        wmax=wr + fo, wmin=wr, cbe=wr + fo / 2.0,
        ap=g["ap"], ar=g["ar"], depth=g["depth"], spm=g["spm"],
        young=g["young"], rodrho=g["rodrho"], stroke=g["stroke"])

FIELDS = dict(wmax="wmax", wmin="wmin", dmax="fo", dmin="zero", wr="wr",
              fo="fo", fomax="fomax", upthrust="upthrust", stretch="stretch",
              stroke="stroke", seff="seff", ev="ev", disp="disp", m="mass",
              work="work", power="power", cbe="cbe")
seen = {}
for name in PUMPS:
    c, row = CARD[name], ROW[name]
    want = dict(
        wmax=c["wmax"], wmin=c["wmin"], dmax=c["fo"], dmin=0.0, wr=c["wr"],
        fo=c["fo"], fomax=c["fomax"], upthrust=c["fomax"] - c["fo"],
        stretch=c["stretch"], stroke=GEOM[name]["stroke"], seff=c["seff"],
        ev=c["ev"], disp=c["disp"], m=RHO * c["disp"], work=c["work"],
        power=c["power"], cbe=c["cbe"])
    exact = {"wmax", "wmin", "wr", "fomax", "stroke", "cbe", "dmin"}
    for key, val in want.items():
        got = row[key]
        tol = TIGHT if key in exact else SOLVE
        e = rel(got, val)
        worst(key, e)
        seen[key] = seen.get(key, 0) + (e <= tol)
        if e > tol:
            check(f"{name}.{key}", False, f"printed {got!r}, closed form {val!r}")
    check(f"{name} prints the card the case implies",
          all(rel(row[k] if k != "dmin" else 0.0, v) <=
              (TIGHT if k in exact else SOLVE)
              for k, v in want.items() if k != "dmin"),
          f"{c['state']}, Fo {c['fo']:.1f} N, df {c['stretch']:.4f} m")
    check(f"{name}'s state is the one the arithmetic says",
          row["rod"] == c["state"], f"printed {row['rod']}, arithmetic {c['state']}")
    check(f"{name}'s swept volume is its row rate",
          near(row["Q"], c["disp"], SOLVE),
          f"Q {row['Q']:.6e}, A_p Seff N {c['disp']:.6e} m3/s")
    check(f"{name}'s face velocity is that rate over the barrel",
          near(row["v"], c["disp"] / c["ap"], SOLVE),
          f"v {row['v']:.6f} m/s")
    check(f"{name}'s whole drop is the machine's",
          row["dp"] == row["dpPump"], f"dp and dpPump both {row['dpPump']:.6e} Pa")
    check(f"{name}'s card sits on the rod weight, not below it",
          near(row["wmax"] - row["wmin"], row["fo"], TIGHT),
          f"Wmax - Wmin = {row['wmax'] - row['wmin']:.3f} = Fo")
    check(f"{name}'s downhole card is the load squared into a rectangle",
          near(row["dmax"], row["fo"], TIGHT) and row["dmin"] == 0.0,
          f"downhole {row['dmin']:.1f} .. {row['dmax']:.1f} N")
    check(f"{name}'s counterbalance is the card's mid-height",
          near(row["cbe"], 0.5 * (row["wmax"] + row["wmin"]), TIGHT),
          f"{row['cbe']:.1f} N")
    check(f"{name}'s buoyed weight is its weight in air minus its buoyancy",
          near(row["wr"], c["wair"] - c["buoy"], TIGHT),
          f"{c['wair']:.1f} - {c['buoy']:.1f} = {row['wr']:.1f} N")
    check(f"{name}'s load is the barrel area times the differential it solved",
          near(row["fo"], c["ap"] * max(c["rise"], 0.0), SOLVE),
          f"A_p {c['ap']:.6e} m2 x dp {max(c['rise'], 0.0):.6e} Pa")
    check(f"{name}'s stretch is Hooke on that load",
          near(row["stretch"], row["fo"] * c["depth"] / (c["ar"] * c["young"]),
               SOLVE),
          f"df = Fo L/(A_r E) = {row['stretch']:.6f} m on A_r {c['ar']:.6e} m2")
check("every card field was checked on every arm",
      all(seen.get(k, 0) == len(PUMPS) for k in FIELDS),
      ", ".join(f"{k} {seen.get(k, 0)}" for k in FIELDS))

# ------------------------------------------------------ the four states, out loud
check("six arms are full, one is stalled, one is idle",
      sum(1 for n in PUMPS if ROW[n]["rod"] == "FULL") >= 6
      and sum(1 for n in PUMPS if ROW[n]["rod"] == "STALLED") == 1
      and sum(1 for n in PUMPS if ROW[n]["rod"] == "NO-LOAD") == 1,
      " ".join(f"{n}:{ROW[n]['rod']}" for n in ARMS))
STALL = [n for n in PUMPS if ROW[n]["rod"] == "STALLED"][0]
IDLE = [n for n in PUMPS if ROW[n]["rod"] == "NO-LOAD"][0]
check(f"{STALL} never leaves bottom dead centre",
      ROW[STALL]["seff"] == 0.0 and ROW[STALL]["ev"] == 0.0
      and ROW[STALL]["Q"] == 0.0 and ROW[STALL]["m"] == 0.0
      and ROW[STALL]["work"] == 0.0 and ROW[STALL]["power"] == 0.0
      and ROW[STALL]["stretch"] > GEOM[STALL]["stroke"],
      f"df {ROW[STALL]['stretch']:.4f} m past an S of {GEOM[STALL]['stroke']} m, "
      f"and the load is still there: Fo {ROW[STALL]['fo']:.0f} N")
check(f"{IDLE} strokes without lifting",
      ROW[IDLE]["fo"] == 0.0 and ROW[IDLE]["stretch"] == 0.0
      and ROW[IDLE]["seff"] == GEOM[IDLE]["stroke"]
      and ROW[IDLE]["ev"] == 1.0 and ROW[IDLE]["work"] == 0.0
      and ROW[IDLE]["power"] == 0.0
      and ROW[IDLE]["wmax"] == ROW[IDLE]["wmin"] == ROW[IDLE]["wr"]
      and ROW[IDLE]["upthrust"] == ROW[IDLE]["fomax"],
      f"a flat line at {ROW[IDLE]['wr']:.0f} N, and no load on it")
check(f"but {IDLE} still delivers its swept volume",
      near(ROW[IDLE]["Q"], CARD[IDLE]["ap"] * GEOM[IDLE]["stroke"]
           * GEOM[IDLE]["spm"] / 60.0, TIGHT),
      f"{ROW[IDLE]['Q']:.6e} m3/s = A_p S N")
check(f"and that volume is not paid for by the rods",
      CARD[IDLE]["rise"] < 0.0
      and CARD[IDLE]["rise"] * ROW[IDLE]["Q"] < 0.0
      and ROW[IDLE]["power"] == 0.0,
      f"dp across the plunger {CARD[IDLE]['rise']:.6e} Pa the wrong way round: "
      f"the annulus charge is what moves it, not the stroke")

# ------------------------------------------- the power identity, arm by arm
for name in PUMPS:
    c = CARD[name]
    if name == IDLE:
        continue
    flows = c["rise"] * ROW[name]["Q"]
    check(f"{name}: the card's area is the pressure rise times the swept volume",
          near(ROW[name]["power"], flows, SOLVE),
          f"dp x Q {flows:.4f} W, Fo Seff N {ROW[name]['power']:.4f} W")
    check(f"{name}: work per stroke times strokes per second is its power",
          near(ROW[name]["power"], ROW[name]["work"] * c["ns"], TIGHT),
          f"{ROW[name]['work']:.3f} J x {c['spm']:.0f} spm")
    check(f"{name}: and the swept volume per stroke is A_p x the effective stroke",
          near(ROW[name]["work"], max(c["rise"], 0.0)
               * (ROW[name]["disp"] / c["ns"]), SOLVE),
          f"dp x A_p Seff = {ROW[name]['work']:.3f} J")

# ------------------------------------------------ the network around each pump
def march_law(g, q):
    """the laminar closed form a constant-property row of this case obeys"""
    return RHO * G * g["elevation"] + 128.0 * MU * g["length"] * q \
        / (math.pi * g["diameter"] ** 4)


for name in MARCH:
    g, row = GEOM[name], ROW[name]
    area = math.pi / 4.0 * g["diameter"] ** 2
    q = row["Q"]
    check(f"{name} re-marches on its own printed rate",
          near(row["dp"], march_law(g, q), SOLVE),
          f"dp {row['dp']:.6e}, rho g dz + 128 mu L Q/(pi d^4) {march_law(g, q):.6e}")
    check(f"{name}'s velocity is that rate over its bore",
          near(row["v"], q / area, SOLVE), f"v {row['v']:.6f} m/s")
    reyn = RHO * abs(row["v"]) * g["diameter"] / MU
    if abs(q) < 1e-15:
        # A stalled pump leaves its produced column static. The report still
        # prints Churchill's 64/Re at the solver's own Re floor, which is a
        # number about the code, not about the fluid, so what is checked here is
        # the thing that is physically true: no rate, no friction, pure weight.
        check(f"{name} carries nothing, so its drop is its own column's weight",
              row["Re"] == 0.0 and row["v"] == 0.0
              and near(row["dp"], RHO * G * g["elevation"], TIGHT),
              f"Re {row['Re']:.1f}, dp {row['dp']:.6e} = rho g dz")
    else:
        check(f"{name} is laminar, so f is 64/Re and not a correlation",
              reyn < 2000.0 and near(row["f"], 64.0 / reyn, SOLVE),
              f"Re {reyn:.1f}, f {row['f']:.6f}")

for name in PUMPS:
    g, c = GEOM[name], CARD[name]
    pfrom, pto = NODE_P[g["src"]], NODE_P[g["dst"]]
    check(f"{name}'s printed dp is its two solved node pressures",
          near(row_dp := ROW[name]["dp"], pfrom - pto, SOLVE),
          f"p({g['src']}) {pfrom:.6e} - p({g['dst']}) {pto:.6e}")
    worst("dpnode", rel(row_dp, pfrom - pto))

# node closure: mass in = mass out at every free node
for node in FREE:
    if node not in NODE_P:
        continue
    into = sum(ROW[n]["Q"] for n, g in GEOM.items() if g["dst"] == node)
    outof = sum(ROW[n]["Q"] for n, g in GEOM.items() if g["src"] == node)
    scale = max(abs(into), abs(outof), abs(FEED.get(node, 0.0)), 1e-300)
    bal = into - outof + FEED.get(node, 0.0)
    check(f"node {node} closes on mass", abs(bal) <= SOLVE * scale,
          f"residual {bal:.3e} m3/s against {scale:.3e}")

# ------------------------------------------------ the classic load, term by term
TUB_OF = {}
for name in MARCH:
    g = GEOM[name]
    if g["elevation"] > 0 and g["src"] in {GEOM[p]["dst"] for p in PUMPS}:
        TUB_OF[g["src"]] = name
PMP_OF = {GEOM[p]["src"]: p for p in PUMPS}

for node, pump in sorted(PMP_OF.items()):
    g = GEOM[pump]
    tub = TUB_OF.get(g["dst"])
    if tub is None:
        continue
    ann = [n for n in MARCH if GEOM[n]["src"] == node]
    if len(ann) != 1:
        continue
    a = GEOM[ann[0]]
    c = CARD[pump]

    def fric(row_name):
        gg = GEOM[row_name]
        return 128.0 * MU * gg["length"] * ROW[row_name]["Q"] \
            / (math.pi * gg["diameter"] ** 4)

    fric_tub, fric_ann = fric(tub), fric(ann[0])
    p_cas = FIXED[a["dst"]]
    p_wh = FIXED[GEOM[tub]["dst"]]
    l_minus_h = GEOM[tub]["elevation"] - a["elevation"]
    classic = (p_wh - p_cas) + RHO * G * l_minus_h + fric_tub - fric_ann
    check(f"{pump}: the load is the level, the heads and both frictions",
          near(c["rise"], classic, SOLVE),
          f"dp {c['rise']:.6e} vs (p_wh-p_cas) {p_wh - p_cas:.4e} "
          f"+ rho g (L-h) {RHO * G * l_minus_h:.4e} + tubing fric {fric_tub:.4e} "
          f"- annulus fric {fric_ann:.4e}")
    check(f"{pump}: upthrust is exactly the part of the ceiling the level keeps",
          near(ROW[pump]["upthrust"], c["fomax"] - c["fo"], SOLVE)
          and near(c["fomax"], RHO * G * GEOM[pump]["depth"] * c["ap"], TIGHT),
          f"Fo_max {c['fomax']:.1f} N, Fo {c['fo']:.1f} N, "
          f"upthrust {ROW[pump]['upthrust']:.1f} N")

# ------------------------------------------------------ what the arms say to each
DEPR, GAUGE, FAST, BUOY, FLOOD = ("RP_DEPR", "RP_GAUGE", "RP_FAST", "RP_BUOY",
                                  "RP_FLOOD")
if all(n in PUMPS for n in (DEPR, GAUGE, FAST, BUOY, FLOOD)):
    check("a heavier string does not change the load, only the stretch",
          CARD[GAUGE]["fo"] > 0 and CARD[DEPR]["fo"] > 0
          and abs(ROW[GAUGE]["fo"] / ROW[DEPR]["fo"] - 1.0) < 0.02
          and ROW[GAUGE]["stretch"] < ROW[DEPR]["stretch"]
          and ROW[GAUGE]["seff"] > ROW[DEPR]["seff"]
          and ROW[GAUGE]["ev"] > ROW[DEPR]["ev"],
          f"7/8 in rod: Fo {ROW[GAUGE]['fo']:.0f} vs {ROW[DEPR]['fo']:.0f} N, "
          f"df {ROW[GAUGE]['stretch']:.4f} vs {ROW[DEPR]['stretch']:.4f} m, "
          f"EV {ROW[GAUGE]['ev']:.4f} vs {ROW[DEPR]['ev']:.4f}")
    check("and df x A_r over Fo is one compliance, L/E, for both strings",
          near(ROW[GAUGE]["stretch"] * CARD[GAUGE]["ar"] / ROW[GAUGE]["fo"],
               ROW[DEPR]["stretch"] * CARD[DEPR]["ar"] / ROW[DEPR]["fo"], SOLVE),
          f"L/E = {GEOM[DEPR]['depth'] / GEOM[DEPR]['young']:.3e} m^3/N "
          f"on either gauge")
    check("a string lighter than the fluid moves the datum, not the card",
          near(ROW[BUOY]["fo"], ROW[DEPR]["fo"], SOLVE)
          and near(ROW[BUOY]["stretch"], ROW[DEPR]["stretch"], SOLVE)
          and near(ROW[BUOY]["seff"], ROW[DEPR]["seff"], SOLVE)
          and near(ROW[BUOY]["ev"], ROW[DEPR]["ev"], SOLVE)
          and near(ROW[BUOY]["disp"], ROW[DEPR]["disp"], SOLVE)
          and near(ROW[BUOY]["work"], ROW[DEPR]["work"], SOLVE),
          f"every load-derived number identical: Fo {ROW[BUOY]['fo']:.1f} N, "
          f"EV {ROW[BUOY]['ev']:.4f}")
    shift = (GEOM[BUOY]["rodrho"] - GEOM[DEPR]["rodrho"]) * G \
        * CARD[DEPR]["ar"] * GEOM[DEPR]["depth"]
    check("and the rod weight moves by exactly the density difference it declared",
          near(ROW[BUOY]["wr"] - ROW[DEPR]["wr"], shift, TIGHT)
          and near(ROW[BUOY]["wmin"] - ROW[DEPR]["wmin"], shift, TIGHT)
          and near(ROW[BUOY]["cbe"] - ROW[DEPR]["cbe"], shift, TIGHT),
          f"(700-850) g A_r L = {shift:.2f} N, and Wmin goes negative: "
          f"{ROW[BUOY]['wmin']:.1f} N")
    check("a buoyant string puts the card across the zero-load axis",
          ROW[BUOY]["wmin"] < 0.0 < ROW[BUOY]["wmax"],
          f"{ROW[BUOY]['wmin']:.1f} .. {ROW[BUOY]['wmax']:.1f} N")
    ratio = GEOM[FAST]["spm"] / GEOM[DEPR]["spm"]
    check("speed scales the volume, and buys its own friction",
          near(ROW[FAST]["disp"] / ROW[DEPR]["disp"], ratio, 5e-3)
          and abs(ROW[FAST]["ev"] - ROW[DEPR]["ev"]) < 5e-3
          and ROW[FAST]["dp"] < ROW[DEPR]["dp"]
          and near(ROW[FAST]["power"], ROW[FAST]["work"]
                   * GEOM[FAST]["spm"] / 60.0, TIGHT),
          "%.0f/%.0f spm gives %.4f of the rate at EV %.4f vs %.4f"
          % (GEOM[FAST]["spm"], GEOM[DEPR]["spm"],
             ROW[FAST]["disp"] / ROW[DEPR]["disp"], ROW[FAST]["ev"],
             ROW[DEPR]["ev"]))
    check("a flooded annulus takes the whole column off the plunger's underside",
          ROW[FLOOD]["fo"] < 0.3 * ROW[DEPR]["fo"]
          and ROW[FLOOD]["upthrust"] > ROW[DEPR]["upthrust"]
          and near(ROW[FLOOD]["fomax"], ROW[DEPR]["fomax"], TIGHT),
          f"Fo {ROW[FLOOD]['fo']:.0f} vs {ROW[DEPR]['fo']:.0f} N against the same "
          f"ceiling {ROW[FLOOD]['fomax']:.0f} N")

RIG = "RP_RIG"
if RIG in PUMPS:
    c = CARD[RIG]
    q = FEED.get(GEOM[RIG]["src"], 0.0)
    want_df = GEOM[RIG]["stroke"] - q / (c["ap"] * c["ns"])
    check("the rig: an imposed rate decides its own stretch",
          q > 0 and near(ROW[RIG]["stretch"], want_df, SOLVE),
          f"S - Q/(A_p N) = {want_df:.6f} m, printed {ROW[RIG]['stretch']:.6f} m")
    check("and that stretch is a differential, through the string's stiffness",
          near(c["rise"], ROW[RIG]["stretch"] * c["ar"] * c["young"]
               / (c["ap"] * c["depth"]), SOLVE),
          f"dp = df A_r E/(A_p L) = {c['rise']:.6e} Pa")
    check("so the row rate is the feed rate, to the digits the solve agreed on",
          near(ROW[RIG]["Q"], q, 1e-6),
          f"Q {ROW[RIG]['Q']:.6e} against a {q:.6e} m3/s feed")

ANNFRI = "RP_ANNFRI"
if ANNFRI in PUMPS and DEPR in PUMPS:
    ann = [n for n in MARCH if GEOM[n]["src"] == GEOM[ANNFRI]["src"]]
    if ann:
        a = GEOM[ann[0]]
        fr = 128.0 * MU * a["length"] * abs(ROW[ann[0]]["Q"]) \
            / (math.pi * a["diameter"] ** 4)
        check("the narrow annulus is carrying the whole rate, and it shows in the load",
          ROW[ann[0]]["Q"] < 0.0 and fr > 5e5,
          f"Q_ann {ROW[ann[0]]['Q']:.6e} m3/s, friction {fr:.4e} Pa")
        check("and the card absorbs it as a term of the plunger's differential",
              near(ROW[ANNFRI]["fo"] - ROW[DEPR]["fo"],
                   (CARD[ANNFRI]["rise"] - CARD[DEPR]["rise"]) * CARD[ANNFRI]["ap"],
                   TIGHT)
              and CARD[ANNFRI]["rise"] > CARD[DEPR]["rise"],
              f"dp {CARD[ANNFRI]['rise']:.6e} vs {CARD[DEPR]['rise']:.6e} Pa "
              f"for the same well with a wide annulus")

# --------------------------------------------------------------- the card drawn
# The report does not spell the parallelogram out: it is fixed by the four
# printed numbers (stroke, stretch, Wmin, Wmax), so the checker draws it and
# measures it. The shoelace on those four vertices is the work per stroke, and
# that is the identity the increment claims -- the card's area is not a figure of
# speech for the energy.
for name in PUMPS:
    r, c = ROW[name], CARD[name]
    df = r["stretch"]
    if r["rod"] == "STALLED":
        # The four vertices would need more travel than the cranks give: the
        # load transfer is not finished before the rod reaches top dead centre,
        # so there is no closed card to measure and the machine delivers nothing.
        check(f"{name}'s load transfer never finishes inside its stroke",
              df > r["stroke"] and r["work"] == 0.0 and r["seff"] == 0.0,
              f"df {df:.4f} m past an S of {r['stroke']} m, area zero")
        continue
    verts = ((0.0, r["wmin"]), (df, r["wmax"]),
             (r["stroke"], r["wmax"]), (r["stroke"] - df, r["wmin"]))
    area = 0.5 * abs(sum(verts[i][0] * verts[(i + 1) % 4][1]
                         - verts[(i + 1) % 4][0] * verts[i][1]
                         for i in range(4)))
    check(f"{name}'s parallelogram closes on the four printed numbers",
          near(area, r["work"], TIGHT)
          and near(verts[2][0] - verts[1][0], r["seff"], TIGHT),
          f"({df:.4f}, {r['wmax']:.1f}) N/m vertices, area {area:.3f} J = "
          f"Fo x Seff = {r['work']:.3f} J")
    check(f"{name}'s effective stroke is that parallelogram's width",
          near(r["stroke"] - df, r["seff"], TIGHT),
          f"S - df = {r['stroke'] - df:.6f} = Seff")

# --------------------------------------------------------------- bit stillness
check("no marching row reported a card",
      all("rod" not in ROW[n] for n in MARCH),
      f"{len(MARCH)} marched rows print Re and f as they always did")
check("and the pump rows report no wall law",
      all("Re" not in ROW[n] and "f" not in ROW[n] for n in PUMPS))

# ------------------------------------------------------------- loud refusals
def run_on(exe, body):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "r.dict")
        with open(path, "w") as fh:
            fh.write(body)
        rr = subprocess.run([exe, path], capture_output=True, text=True, cwd=tmp)
    return rr.returncode, (rr.stdout + rr.stderr).strip()


def refuses(label, body, want, exe=None):
    code, out = run_on(exe or STEADY, body)
    check(f"refused: {label}", code != 0 and want in out,
          out.splitlines()[-1][:150] if out else "no message")


def block_span(body, name):
    m = re.search(r"\b" + re.escape(name) + r"\s*\{", body)
    if not m:
        raise SystemExit(f"probe: the case has no {name} block")
    start, depth = m.end() - 1, 0
    for j in range(start, len(body)):
        if body[j] == "{":
            depth += 1
        elif body[j] == "}":
            depth -= 1
            if depth == 0:
                return start, j + 1
    raise SystemExit(f"probe: unterminated {name} block")


PUMP_HEAD = "        from P2; to D2; length 0; diameter 0;"
GEOM_SET = "rodPlungerDiameter 0.044450; rodStrokeLength 1.5;"
STRING_SET = "rodStrokesPerMinute 8; rodDiameter 0.019050;"
DEEP_SET = "rodLength 1500; rodYoungModulus 2.0e11; rodDensity 7850;"
for anchor in (PUMP_HEAD, GEOM_SET, STRING_SET, DEEP_SET):
    lo, hi = block_span(src, DEPR)
    if anchor not in src[lo:hi]:
        raise SystemExit(f"probe: {anchor!r} is no longer in the {DEPR} block")


def patched(replace=None, extra="", scope=DEPR):
    lo, hi = block_span(src, scope)
    head = src[lo:hi]
    if extra:
        new = head.replace(PUMP_HEAD, PUMP_HEAD + "\n        " + extra, 1)
        if new == head:
            raise SystemExit("probe: the pump anchor is not in the scope")
        head = new
    if replace:
        old, new = replace
        cut = head.replace(old, new, 1)
        if cut == head:
            raise SystemExit(f"probe: {old!r} is not in the {scope} scope")
        head = cut
    return src[:lo] + head + src[hi:]


refuses("a rod pump with a length",
        patched((PUMP_HEAD, "        from P2; to D2; length 900; diameter 0;")),
        "rod pump with a length")
refuses("a rod pump with a bore",
        patched((PUMP_HEAD, "        from P2; to D2; length 0; diameter 0.062;")),
        "rod pump with a diameter")
refuses("a rod pump with an elevation", patched(extra="elevation 10;"),
        "rod pump with an elevation")
refuses("a pump with no stroke", patched((GEOM_SET, "rodPlungerDiameter 0.044450;")),
        "rod pump without a stroke")
refuses("a pump with no speed", patched((STRING_SET, "rodDiameter 0.019050;")),
        "rod pump without a speed")
refuses("a pump with no string", patched((STRING_SET, "rodStrokesPerMinute 8;")),
        "rod pump without a rod")
refuses("a pump with no setting depth",
        patched((DEEP_SET, "rodYoungModulus 2.0e11; rodDensity 7850;")),
        "rod pump without a setting depth")
refuses("a pump with no modulus",
        patched((DEEP_SET, "rodLength 1500; rodDensity 7850;")),
        "rod pump without a modulus")
refuses("a pump whose string has no density",
        patched((DEEP_SET, "rodLength 1500; rodYoungModulus 2.0e11;")),
        "non-positive rod density")
refuses("a pump with a negative rod density",
        patched((DEEP_SET, DEEP_SET.replace("7850", "-7850"))),
        "non-positive rod density")
refuses("a rod wider than its barrel",
        patched((STRING_SET, STRING_SET.replace("0.019050", "0.044450"))),
        "rod at least as wide as its plunger")
refuses("a negative barrel bore",
        patched((GEOM_SET, GEOM_SET.replace("0.044450", "-0.044450"))),
        "negative rodPlungerDiameter")
refuses("a pump parameter with no barrel to switch it on",
        patched((GEOM_SET, "rodStrokeLength 1.5;")),
        "rod-pump parameter without rodPlungerDiameter")

# The pump's own scope refuses anything that would price a wall loss, another
# machine or a two-phase dispatch on a point element. Most of these keys are
# refused one gate earlier by the law that owns them on any pipe at all -- an
# annulus bore needs a bore to be inside, a rheology index needs its units, a
# stage count needs a curve, a dome needs a gas table -- and the probe asserts
# the message the user actually gets rather than the one this increment added.
CO_DECL = (
    ("roughness 1e-5;", "rod pump that also declares"),
    ("valve 0.5;", "rod pump that also declares"),
    ("pumpHead 10;", "rod pump that also declares"),
    ("pumpIntakePressure 1e6;", "rod pump that also declares"),
    ("yieldStress 10;", "rod pump that also declares"),
    ("ductHeight 0.01;", "rod pump that also declares"),
    ("coupledResistance 1e6;", "rod pump that also declares"),
    ("dInner 0.05;", "needs 0 < annulusInner < diameter"),
    ("eccentricity 0.01;", "eccentricity without an inner bore"),
    ("flowIndex 0.6;", "flow index without a positive consistency"),
    ("consistency 0.5;", "consistency without a flowIndex"),
    ("flowPattern 2;", "flowPattern 2 on a single-phase network"),
    ("pumpStages 2;", "pumpStages/pumpSpeed without a pump"),
    ("pumpGasFraction 0.1;", "pumpGasFraction without a pump"),
    ("liftDomePressure 7.0e6;", "gas-lift valve on a constant-property network"))
for decl, want in CO_DECL:
    refuses(f"a pump that also declares {decl[:-1]}", patched(extra=decl), want)

FILL_END = text.index("nodes")
FILL_START = text.index("fluid")
for fill, want, tag in (
        ("fluid { pvt { pressure ( 1.0e6 1.0e7 ); rho ( 800.0 900.0 );"
         " mu ( 0.1 0.1 ); } }", "property-varying fill", "tabulated"),
        ("fluid { blackoil { oilDensity 850; } }", "black-oil network",
         "black-oil")):
    refuses(f"a pump on a {tag} fill",
            text[:FILL_START] + fill + "\n\n" + text[FILL_END:], want)
refuses("a zero-length row that is not a pump",
        src.replace("TUB2 { from D2; to W2; length 1500;",
                    "TUB2 { from D2; to W2; length 0;", 1),
        "needs positive length and diameter")

# ------------------------------------------------- a signed upthrust, if asked
# The element states that upthrust is signed: a level below the surface and a
# casing pressure take weight off the plunger's underside (positive), tubing
# friction and a wellhead backpressure put weight on it (negative). No arm of
# the case reaches the negative side -- every well runs to a 2 MPa header and the
# column it lifts weighs 12.5 MPa -- so the sign gets its witness here, on the
# flooded well with its separator wound up past the whole column. The ceiling
# then stops being a ceiling: the load the string carries is more than the
# frictionless column, and the card says so instead of clamping.
BACK = src.replace("W1 { pressure 2.0e6; }", "W1 { pressure 1.4e7; }", 1)
if BACK == src:
    raise SystemExit("probe: the case no longer fixes W1 at 2 MPa")
code, out = run_on(STEADY, BACK)
check("a flooded well against a 14 MPa header still solves",
      code == 0 and "converged" in out.splitlines()[-1],
      out.splitlines()[-1].replace("[chflow] ", ""))
bp, br, _ = rows_from(out)
bg, b = GEOM[FLOOD], br[FLOOD]
fric = lambda r: br[r]["dp"] - RHO * G * GEOM[r]["elevation"]
ceiling = RHO * G * bg["depth"] * bg["ap"]
# The annulus is full to the surface here, so the two columns' weights cancel
# and the plunger's differential is the header over the casing plus the two
# frictions it has to beat.
expected = (bp["W1"] - bp["H1"]) + fric("TUB1") - fric("ANN1")
check("and the backpressure is what the plunger now holds",
      near(-b["dp"], expected, SOLVE) and near(b["fomax"], ceiling, TIGHT),
      f"dp {-b['dp']:.6e} vs (p_wh - p_cas) {bp['W1'] - bp['H1']:.6e}"
      f" + tubing fric {fric('TUB1'):.6e} - annulus fric {fric('ANN1'):.6e}"
      f" Pa, against a {ceiling:.1f} N ceiling")
check("so the load passes its ceiling and upthrust turns negative",
      b["fo"] > ceiling and b["upthrust"] < 0.0
      and near(b["upthrust"], ceiling - b["fo"], TIGHT),
      f"fo {b['fo']:.1f} N > fomax {ceiling:.1f} N, upthrust {b['upthrust']:.2f} N")
check("and the card is still the same four statements under it",
      b["rod"] == "FULL"
      and near(b["fo"], bg["ap"] * -b["dp"], TIGHT)
      and near(b["stretch"], b["fo"] * bg["depth"] / (bg["ar"] * bg["young"]), TIGHT)
      and near(b["seff"], bg["stroke"] - b["stretch"], TIGHT),
      f"{b['rod']}: df {b['stretch']:.4f} m of an S {bg['stroke']} m, fo"
      f" {b['fo']:.1f} N, power {b['power']:.1f} W")

for exe, case_name, tag in ((TRANSIENT, "transientBetaMomentum", "transient: pipe"),
                            (TWOFLUID, "twoFluidPressureLaminar", "twoFluid: pipe")):
    if not exe:
        continue
    other = open(os.path.join(HERE, "cases", case_name, "case.dict")).read()
    i = other.index("diameter")
    j = other.index(";", i) + 1
    code, out = run_on(exe, other[:j] + " rodPlungerDiameter 0.044450;" + other[j:])
    check(f"refused: a rod pump on the {case_name} row",
          code != 0 and tag in out and "rod pump" in out,
          out.splitlines()[-1][:150] if out else "no message")
if NODAL:
    nd = open(os.path.join(HERE, "cases", "nodalVogel", "case.dict")).read()
    key = "    elevation 2200;    // m, wellhead above the pay"
    if key not in nd:
        raise SystemExit("probe: the nodalVogel well block moved")
    code, out = run_on(NODAL, nd.replace(key,
                                         key + "\n    rodPlungerDiameter 0.044450;"))
    check("refused: a rod pump on the nodal VLP well",
          code != 0 and "nodal: well declares a rod pump" in out,
          out.splitlines()[-1][:150] if out else "no message")

if fails:
    raise SystemExit(f"rod-pump card verification FAILED: {len(fails)}: {fails}")
print(f"the card is four statements and one parallelogram: Wr, Fo = A_p dp, "
      f"df = Fo L/(A_r E), Seff = S - df restated from the case on {len(PUMPS)} "
      f"arms (worst {max(err.values()):.2g} rel over "
      f"{', '.join(f'{k} {v:.2g}' for k, v in sorted(err.items()))}), the whole "
      f"load re-summed from the level, both heads and both column frictions, "
      f"{sum(1 for n in PUMPS if ROW[n]['rod'] == 'FULL')} full / "
      f"1 stalled / 1 idle, one compliance L/E across two string gauges, a "
      f"density change that moves only the datum, an imposed rate that finds its "
      f"own differential, the card area = dp x swept volume on every arm that "
      f"lifts, and a separator that outweighs the whole column priced as a "
      f"negative upthrust. Untouched paths are bit-still -- the standing suite "
      f"is the witness.")
