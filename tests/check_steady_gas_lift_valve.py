#!/usr/bin/env python3
"""Verification of the IPO gas-lift valve (5.89, D1 increment 3).

A valve pipe is not a pipe with a valve in it: it is a bellows force balance
that decides whether gas may pass at all, and an isentropic nozzle that prices
how much does. The checker never reads the solver's internals. It takes the
valve's settings and the manifold's property table out of the case dictionary,
the solved node pressures out of the 12-digit report, and restates both halves
of the law in closed form:

  force balance : Pvo = Pd/(1-R) + St - R/(1-R) Pt, Pvc = Pd + St(1-R), and the
                  spread, the production effect and the test-rack opening that
                  go with them - re-priced at the tubing pressure the arm's own
                  flow produced;
  nozzle        : mdot = Cd Ap sqrt(2 rho1 p1 k/(k-1) [r^(2/k) - r^((k+1)/k)]),
                  sonic for r <= r* = (2/(k+1))^(1/(k-1)), on the density the
                  table gives at the manifold pressure;
  disposition   : which side of its own thresholds the solved manifold pressure
                  came in on, including the third answer a steady solve can
                  give - inside the spread, where there is no state to inherit
                  and the arm is priced shut.

What the arms are then made to prove against each other, from the printed rows
alone:

  a choked port is deaf : the choked arms into shoes a megapascal apart carry
                          the SAME mass rate, to the last digit printed, while
                          their drops differ by that megapascal;
  the bellows is not the port: the 0.05-ratio valve and the 0.10-ratio valve on
                          the same dome, port and shoe print the same rate and
                          different opening pressures - and the same closing
                          pressure, because R never reaches Pvc;
  a spring moves the band: the springless arm is open, the 0.5 MPa arm is inside
                          its spread, the 1.5 MPa arm is shut - one knob, three
                          answers, and Pvc rises by exactly St(1-R) each time;
  a dead arm is dead    : an arm with nothing to inherit holds its tubing at the
                          shoe pressure exactly and passes nothing;
  the flux is Bernoulli with the expansion priced in: at r ~ 0.98 the port is
                          within 2% of sqrt(2 rho dp) and under it;
  the arms hang on a real network: every marching tubing re-marched from the
                          reported upstream pressure at the reported rate lands
                          on the reported downstream pressure, and the manifold
                          balance closes against the tie-in plus the compressor.

usage: check_steady_gas_lift_valve.py <steady-exe> <case.dict> \
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

fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def near(got, want, rtol):
    return rel(got, want) <= rtol


# ---------------------------------------------------------------- run the case
proc = subprocess.run([STEADY, CASE], capture_output=True, text=True)
if proc.returncode != 0:
    raise SystemExit(f"solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")

NODE_P, ROW, RAW = {}, {}, {}
for line in proc.stdout.splitlines():
    tok = line.replace("[chflow]", "").split()
    if not tok:
        continue
    if tok[0] == "node" and len(tok) >= 4:
        NODE_P[tok[1]] = float(tok[3])
    elif tok[0] == "pipe":
        f, raw = {}, {}
        j = 2
        while j < len(tok) - 1:
            if tok[j] == "valve":
                f["valve"] = tok[j + 1]
                j += 2
            elif tok[j] == "SONIC":
                f["SONIC"] = 1.0
                j += 1
            else:
                f[tok[j]] = float(tok[j + 1])
                raw[tok[j]] = tok[j + 1]
                j += 2
        ROW[tok[1]] = f
        RAW[tok[1]] = raw
if not ROW:
    raise SystemExit("no pipe lines parsed")

last = proc.stdout.strip().splitlines()[-1]
check("the network solved", "converged" in last and "NOT converged" not in last,
      last.replace("[chflow] steady solved: ", ""))
iters = float(re.search(r"iterations=(\S+)", last).group(1))
relres = float(re.search(r"residual=(\S+)", last).group(1))
check("and it did not need the whole budget to say so", 0 < iters < 200,
      f"{iters:.0f} Newton iterations, relative residual {relres:.3g}")

# ---------------------------------------------------------------- the case text
src = open(CASE).read()
text = re.sub(r"//[^\n]*", "", src)  # comment-stripped, for parsing


def block(body, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", body)
    if not m:
        raise SystemExit(f"case has no {key}{{}} block")
    start = m.end() - 1
    depth = 0
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
        start = m.end() - 1
        depth = 0
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


def listof(body, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(([^)]*)\)", body)
    return [float(x) for x in m.group(1).split()] if m else []


PVT = block(block(text, "fluid"), "pvt")
P_GRID, RHO_TAB, MU_TAB = listof(PVT, "pressure"), listof(PVT, "rho"), \
    listof(PVT, "mu")
if not (len(P_GRID) == len(RHO_TAB) == len(MU_TAB)) or len(P_GRID) < 2:
    raise SystemExit("case table malformed")

NODES = entries(block(text, "nodes"))
PIPES = entries(block(text, "pipes"))
SOURCES = entries(block(text, "sources"))

GEOM = {}
for name, body in PIPES.items():
    GEOM[name] = dict(
        src=word(body, "from"), dst=word(body, "to"),
        length=num(body, "length", 0.0), diameter=num(body, "diameter"),
        roughness=num(body, "roughness", 0.0),
        elevation=num(body, "elevation", 0.0),
        dome=num(body, "liftDomePressure", 0.0),
        ratio=num(body, "liftAreaRatio", 0.0),
        spring=num(body, "liftSpringPressure", 0.0),
        cd=num(body, "liftDischargeCoeff", 0.0),
        kappa=num(body, "liftKappa", 0.0))
VALVES = sorted(n for n, g in GEOM.items() if g["dome"] > 0)
MARCHES = sorted(n for n, g in GEOM.items() if g["dome"] <= 0)
check("the case is a valve case", len(VALVES) >= 6 and len(MARCHES) >= 6,
      f"{len(VALVES)} valve rows, {len(MARCHES)} marching rows")


def interp(vals, p):
    """piecewise-linear lookup clamped at the table ends (as PvtFluid)"""
    if p <= P_GRID[0]:
        return vals[0]
    if p >= P_GRID[-1]:
        return vals[-1]
    lo, hi = 0, len(P_GRID) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if P_GRID[mid] <= p:
            lo = mid
        else:
            hi = mid
    t = (p - P_GRID[lo]) / (P_GRID[hi] - P_GRID[lo])
    return vals[lo] + (vals[hi] - vals[lo]) * t


def rho(p):
    return interp(RHO_TAB, p)


# The case claims its table lies on one straight line through the origin, so the
# fill is an ideal gas at a fixed Rs T and every density below is a closed form
# rather than a lookup only the solver could have checked.
RS_T = [P_GRID[i] / RHO_TAB[i] for i in range(len(P_GRID))]
MID_RHO = [rho(0.5 * (P_GRID[i] + P_GRID[i + 1])) for i in range(len(P_GRID) - 1)]
MID_EXACT = [0.5 * (P_GRID[i] + P_GRID[i + 1]) / RS_T[0]
             for i in range(len(P_GRID) - 1)]
check("the table is one straight line through the origin (rho = P/RsT)",
      max(rel(v, RS_T[0]) for v in RS_T) < 1e-10
      and max(rel(a, b) for a, b in zip(MID_RHO, MID_EXACT)) < 1e-12,
      f"Rs T = {RS_T[0]:.6f} J/kg, mu flat at {MU_TAB[0]:.1e} Pa s")

# ------------------------------------------------------------ the valve's halves
ATM = 101325.0
G = 9.80665
NSEG = 24


def bore_area(g):
    return 0.25 * math.pi * g["diameter"] * g["diameter"]


def thresholds(g, p_tubing):
    """the bellows' own settings: absolute pressures"""
    r = g["ratio"]
    pvo = g["dome"] / (1.0 - r) + g["spring"] - (r / (1.0 - r)) * p_tubing
    pvc = g["dome"] + g["spring"] * (1.0 - r)
    return dict(pvo=pvo, pvc=pvc, spread=pvo - pvc,
                ppe=(r / (1.0 - r)) * p_tubing,
                tro=g["dome"] / (1.0 - r) + g["spring"]
                    - (r / (1.0 - r)) * ATM)


def flux_bracket(r, k):
    """the isentropic flux's own shape factor, unclamped"""
    return r ** (2.0 / k) - r ** ((k + 1.0) / k)


_RSTAR = {}


def rstar(k):
    """the ratio the flux brackets at, located by searching the bracket

    The throat chokes where the mass flux is at its maximum, so the checker
    finds that ratio by ternary search on the flux expression it is already
    using rather than by quoting a critical-ratio formula - which is what makes
    these rate checks a test of the solver's threshold, not a restatement of it.
    """
    if k not in _RSTAR:
        lo, hi = 1e-12, 1.0
        for _ in range(300):
            m1 = lo + (hi - lo) / 3.0
            m2 = hi - (hi - lo) / 3.0
            if flux_bracket(m1, k) < flux_bracket(m2, k):
                lo = m1
            else:
                hi = m2
        _RSTAR[k] = 0.5 * (lo + hi)
    return _RSTAR[k]


def nozzle(g, p_up, p_down):
    """the isentropic port: mass flux on the upstream state, sonic at r*"""
    k = g["kappa"]
    rc = rstar(k)
    out = dict(massRate=0.0, chedzy=0.0, choked=False, r=0.0, rstar=rc)
    if p_down >= p_up:
        return out
    rho1 = rho(p_up)
    out["r"] = p_down / p_up
    out["choked"] = out["r"] <= rc
    out["chedzy"] = (g["cd"] * bore_area(g)
                     * math.sqrt(2.0 * rho1 * (p_up - p_down)))
    # Past the critical ratio the throat no longer feels the downstream, so the
    # flux holds at its maximum G(r*) rather than climbing the unphysical branch
    # the algebra alone would give for r < r*.
    rs = max(out["r"], rc)
    bracket = flux_bracket(rs, k)
    out["massRate"] = (g["cd"] * bore_area(g)
                       * math.sqrt(max(2.0 * rho1 * p_up * (k / (k - 1.0))
                                       * bracket, 0.0)))
    return out


DISPOSITION = {1: "OPEN", 0: "CLOSED", 2: "IN-SPREAD"}

# ------------------------------------------------- every valve row, from scratch
err = {}


def note(key, got, want):
    err[key] = max(err.get(key, 0.0), rel(got, want))


for name in VALVES:
    g = GEOM[name]
    row, raw = ROW.get(name), RAW.get(name)
    if row is None or "valve" not in row:
        check(f"{name}: reported as a valve row", False, "no valve token")
        continue
    p_up, p_dn = NODE_P[g["src"]], NODE_P[g["dst"]]
    th, fl = thresholds(g, p_dn), nozzle(g, p_up, p_dn)
    disp = 1 if p_up >= th["pvo"] else (0 if p_up <= th["pvc"] else 2)

    check(f"{name}: disposition", row["valve"] == DISPOSITION[disp],
          f"{row['valve']}, with the manifold at {p_up:.1f} against Pvo "
          f"{th['pvo']:.1f} and Pvc {th['pvc']:.1f}")
    for key, want in (("pvo", th["pvo"]), ("pvc", th["pvc"]),
                      ("spread", th["spread"]), ("ppe", th["ppe"]),
                      ("tro", th["tro"])):
        check(f"{name}: {key} is the force balance at this tubing pressure",
              near(row[key], want, 1e-9),
              f"reported {row[key]:.8g}, restated {want:.8g}")
        note(key, row[key], want)
    check(f"{name}: the test rack is the same valve with the tubing side open",
          near(row["tro"] - row["pvo"],
               (g["ratio"] / (1.0 - g["ratio"])) * (p_dn - ATM), 1e-9)
          and row["tro"] >= row["pvo"] >= row["pvo"] - 1e9,
          f"tro - pvo = {row['tro'] - row['pvo']:.6g} Pa")
    check(f"{name}: R, St, Cd and kappa are all read from the case",
          all(has(PIPES[name], k) for k in
              ("liftDomePressure", "liftAreaRatio", "liftDischargeCoeff",
               "liftKappa")),
          f"R = {g['ratio']}, St = {g['spring']:.1e}, Cd = {g['cd']}, "
          f"k = {g['kappa']}")
    check(f"{name}: no wall law is priced on a valve row",
          "Re" not in row and "f" not in row and "HL" not in row
          and "He" not in row and "mapRegime" not in row,
          " ".join(sorted(row.keys())))
    check(f"{name}: the whole drop is the port's",
          near(row["dp"], p_up - p_dn, 1e-9)
          and raw["dpValve"] == raw["dp"],
          f"dp = dpValve = {row['dp']:.6g} Pa")
    if disp == 1:
        check(f"{name}: the throat's verdict",
              ("SONIC" in row) == fl["choked"],
              f"r = {fl['r']:.6f} against r* = {fl['rstar']:.6f}")
    else:
        check(f"{name}: a shut port has no throat to report",
              "SONIC" not in row,
              f"{row['valve']}, so r = {fl['r']:.6f} against r* = "
              f"{fl['rstar']:.6f} is a question the arm never has to answer")

    want_m = fl["massRate"] if disp == 1 else 0.0
    check(f"{name}: mass rate through the port",
          near(row["m"], want_m, 1e-9),
          f"reported {row['m']:.10g} kg/s, nozzle {want_m:.10g}")
    note("m", row["m"], want_m)
    if disp != 1:
        check(f"{name}: priced shut passes nothing",
              row["m"] == 0 and row["Q"] == 0 and row["v"] == 0,
              f"{row['valve']} at a {row['dp']:.6g} Pa drop")
        continue
    rho_up = rho(p_up)
    check(f"{name}: the volume rate is priced at the manifold density",
          near(row["Q"], row["m"] / rho_up, 1e-11),
          f"Q {row['Q']:.10g} = m/rho({p_up:.0f}) = {row['m']:.6g}/{rho_up:.6g}")
    check(f"{name}: the bore velocity is that volume rate over the port",
          near(row["v"], row["Q"] / bore_area(g), 1e-11),
          f"v {row['v']:.6g} m/s over {bore_area(g) * 1e6:.4f} mm^2")
    check(f"{name}: the compressible flux never beats Bernoulli",
          fl["massRate"] <= fl["chedzy"],
          f"m/m_Chedzy = {fl['massRate'] / fl['chedzy']:.6f} at r = {fl['r']:.4f}")

# ------------------------------------------- what the arms were built to show
DESIGN = {
    "VLV_SUB": ("OPEN", False), "VLV_CHOK": ("OPEN", True),
    "VLV_TWIN": ("OPEN", True), "VLV_UP": ("OPEN", True),
    "VLV_R05": ("OPEN", False), "VLV_WIDE": ("OPEN", False),
    "VLV_SPR": ("IN-SPREAD", None),
    "VLV_HARD": ("CLOSED", None), "VLV_CLOSE": ("CLOSED", None),
    "VLV_LO": ("CLOSED", None),
}
missing = [n for n in DESIGN if n not in ROW]
check("every designed arm is still in the case", not missing, str(missing))
for name, (want_disp, want_sonic) in DESIGN.items():
    if name in missing:
        continue
    row = ROW[name]
    check(f"{name}: still demonstrating {want_disp}"
          + ("" if want_sonic is None else f", sonic={want_sonic}"),
          row["valve"] == want_disp
          and (want_sonic is None or ("SONIC" in row) == want_sonic),
          f"{row['valve']}{' SONIC' if 'SONIC' in row else ''}")
check("the arms cover all three dispositions and both throat states",
      {ROW[n]["valve"] for n in VALVES} == {"OPEN", "CLOSED", "IN-SPREAD"}
      and any("SONIC" in ROW[n] for n in VALVES)
      and any("SONIC" not in ROW[n] for n in VALVES
              if ROW[n]["valve"] == "OPEN"),
      ", ".join(sorted({ROW[n]["valve"] for n in VALVES})))

# A choked port is deaf, and the report has to show it: three arms off one
# manifold, shoes a megapascal apart, differ in drop by that megapascal and in
# nothing else.
deaf = [n for n in ("VLV_CHOK", "VLV_TWIN", "VLV_UP") if n in ROW]
if len(deaf) == 3:
    check("a choked port cannot hear its downstream",
          len({RAW[n]["m"] for n in deaf}) == 1
          and len({RAW[n]["Q"] for n in deaf}) == 1
          and len({RAW[n]["v"] for n in deaf}) == 1
          and len({RAW[n]["dp"] for n in deaf}) == 3,
          f"m = {RAW['VLV_CHOK']['m']} kg/s on all three, drops "
          + " / ".join(RAW[n]["dp"] for n in deaf))
    check("and it is the throat, not the tolerance, that says so",
          all(nozzle(GEOM[n], NODE_P[GEOM[n]["src"]],
                     NODE_P[GEOM[n]["dst"]])["choked"] for n in deaf)
          and max(rel(ROW[n]["m"], ROW[deaf[0]]["m"]) for n in deaf) == 0.0,
          f"every one of them at r < "
          f"{nozzle(GEOM['VLV_CHOK'], NODE_P['ANN'], NODE_P['T2'])['rstar']:.6f}")
    open_sub = [n for n in ("VLV_SUB", "VLV_R05") if "SONIC" not in ROW.get(n, {})]
    unchoked = [n for n in open_sub
                if ROW[n]["m"] < ROW["VLV_CHOK"]["m"] * (1.0 - 1e-9)]
    check("while an unchoked port at the same manifold passes less",
          len(unchoked) == len(open_sub) > 0,
          ", ".join(f"{n} {ROW[n]['m']:.6g}" for n in open_sub)
          + f" against {ROW['VLV_CHOK']['m']:.6g} kg/s")

# R is a bellows statement, not a port statement.
if {"VLV_SUB", "VLV_R05"} <= set(ROW):
    check("the area ratio never enters the nozzle",
          RAW["VLV_R05"]["m"] == RAW["VLV_SUB"]["m"]
          and RAW["VLV_R05"]["Q"] == RAW["VLV_SUB"]["Q"],
          f"R = 0.10 and R = 0.05 both pass m = {RAW['VLV_SUB']['m']} kg/s")
    check("while R moves Pvo and the spread and leaves Pvc alone",
          RAW["VLV_R05"]["pvc"] == RAW["VLV_SUB"]["pvc"]
          and ROW["VLV_R05"]["pvo"] < ROW["VLV_SUB"]["pvo"]
          and ROW["VLV_R05"]["spread"] < ROW["VLV_SUB"]["spread"] > 0,
          f"Pvc {RAW['VLV_SUB']['pvc']} shared; Pvo {ROW['VLV_SUB']['pvo']:.1f} "
          f"against {ROW['VLV_R05']['pvo']:.1f}")
    for a, b in (("VLV_SUB", 0.10), ("VLV_R05", 0.05)):
        check(f"{a}: the case still declares R = {b}",
              near(GEOM[a]["ratio"], b, 1e-12))

# The spring knob: one dome, one shoe, three answers.
knob = [n for n in ("VLV_CHOK", "VLV_SPR", "VLV_HARD") if n in ROW]
if len(knob) == 3:
    for a, b in zip(knob, knob[1:]):
        d_st = GEOM[b]["spring"] - GEOM[a]["spring"]
        r_b = GEOM[b]["ratio"]
        check(f"{a} -> {b}: the spring raises Pvc by St(1-R)",
              near(GEOM[a]["dome"], GEOM[b]["dome"], 1e-12)
              and near(ROW[b]["pvc"] - ROW[a]["pvc"], d_st * (1.0 - r_b), 1e-9),
              f"dPvc {ROW[b]['pvc'] - ROW[a]['pvc']:.1f} Pa for a "
              f"dSt(1-R) {d_st * (1.0 - r_b):.1f} Pa")
    check("and the manifold lands in a different band in each arm",
          ROW["VLV_CHOK"]["pvo"] < NODE_P["ANN"]
          and ROW["VLV_SPR"]["pvc"] < NODE_P["ANN"] < ROW["VLV_SPR"]["pvo"]
          and NODE_P["ANN"] < ROW["VLV_HARD"]["pvc"],
          f"manifold {NODE_P['ANN']:.1f} Pa: open above "
          f"{ROW['VLV_CHOK']['pvo']:.0f}, in-band up to "
          f"{ROW['VLV_SPR']['pvo']:.0f}, shut above "
          f"{ROW['VLV_HARD']['pvc']:.0f}")

# An arm with nothing to inherit is exactly its shoe: no flow, no drop, and a
# tubing pressure that is a boundary condition rather than a solved number.
for name in VALVES:
    row = ROW[name]
    if row["valve"] == "OPEN":
        continue
    dn = GEOM[name]["dst"]
    shoe = num(NODES[dn], "pressure")
    runs = [t for t in MARCHES if GEOM[t]["src"] == dn]
    if shoe is None or not runs:
        continue
    for t in runs:
        check(f"{t}: a shut arm holds its tubing at the shoe",
              ROW[t]["Q"] == 0 and ROW[t]["dp"] == 0
              and near(NODE_P[dn], shoe, 1e-12),
              f"{dn} at {NODE_P[dn]:.3f} = the shoe, Q = {RAW[t]['Q']}, "
              f"dp = {RAW[t]['dp']}")

# And the dead marching row's wall-law diagnostic is a number: Churchill at the
# Reynolds floor is 64/Re_min, which is the limit the laminar branch tends to,
# not a divide by zero dressed up as a friction factor.
for t in MARCHES:
    if "Re" not in ROW[t] or ROW[t]["Q"] != 0:
        continue
    check(f"{t}: a dead marching row reports the laminar floor 64/Re_min",
          near(ROW[t]["f"], 64.0 / 1e-12, 1e-12), f"f = {RAW[t]['f']}")

# ------------------------------------------ the compressible flux is Bernoulli
wide = "VLV_WIDE"
if wide in ROW:
    g = GEOM[wide]
    fl = nozzle(g, NODE_P[g["src"]], NODE_P[g["dst"]])
    ratio = fl["massRate"] / fl["chedzy"]
    check(f"{wide}: r ~ 1, so the port is Bernoulli to within the expansion",
          0.98 < fl["r"] < 0.99 and 0.97 < ratio < 1.0,
          f"r = {fl['r']:.6f}, m/m_Chedzy = {ratio:.6f}")
    check(f"{wide}: and the exact flux, not the Chedzy one, is what reported",
          near(ROW[wide]["m"], fl["massRate"], 1e-9)
          and not near(ROW[wide]["m"], fl["chedzy"], 1e-3),
          f"reported {ROW[wide]['m']:.8g} against Bernoulli's {fl['chedzy']:.8g}")
    check(f"{wide}: inverted spread, and running open",
          ROW[wide]["spread"] < 0 and ROW[wide]["valve"] == "OPEN",
          f"spread {ROW[wide]['spread']:.1f} Pa: the tubing holds the port off "
          f"its seat, so there is no band to sit in")

# The flux function's own anchors, restated where no arm can reach.
k0 = GEOM[VALVES[0]]["kappa"]
RC = rstar(k0)
RC_FORM = (2.0 / (k0 + 1.0)) ** (k0 / (k0 - 1.0))


def bracket(r):
    return flux_bracket(max(r, RC), k0)


check("the searched throat ratio is the closed form",
      near(RC_FORM, RC, 1e-6),
      f"r* = {RC:.9f}: (2/(k+1))^(k/(k-1)) for kappa = {k0}, located by search "
      f"to the width of the peak's flat top ({abs(RC_FORM - RC):.1e} apart)")
# The one place a valve law can be quietly wrong is which ratio it chokes at:
# quoting the exponent 1/(k-1) - the throat DENSITY ratio's exponent - still
# gives a continuous flux with a maximum that merely sits in the wrong place,
# and every passing arm is mis-priced by the bracket between the two.
WRONG = (2.0 / (k0 + 1.0)) ** (1.0 / (k0 - 1.0))
check("and the wrong exponent is not the same law at a different ratio",
      abs(flux_bracket(WRONG, k0) - flux_bracket(RC, k0))
      / flux_bracket(RC, k0) > 1e-3,
      f"the k-1 exponent would choke at {WRONG:.6f} and mis-price a passing "
      f"arm by {100.0 * (1.0 - math.sqrt(flux_bracket(WRONG, k0) / flux_bracket(RC, k0))):.2f}%")
edge = 1.0 - 1e-7
check("the bracket opens as Chedzy: k/(k-1) [ ] -> (1-r) as r -> 1",
      abs(k0 / (k0 - 1.0) * bracket(edge) / (1.0 - edge) - 1.0) < 1e-6,
      f"at r = {edge:.8f} the ratio is "
      f"{k0 / (k0 - 1.0) * bracket(edge) / (1.0 - edge):.9f}")
check("the bracket is maximal, and continuous, at the searched ratio",
      all(bracket(RC) > bracket(RC + d)
          and bracket(RC) >= bracket(RC - d)
          for d in (1e-6, 1e-3, 0.05))
      and abs(bracket(RC + 1e-9) - bracket(RC - 1e-9)) < 1e-12,
      f"the flux holds flat below r* = {RC:.6f}, so the law has no step")
check("and it falls monotonically as the downstream is felt harder",
      all(bracket(r) > bracket(r + 0.01)
          for r in (RC + 1e-9, 0.7, 0.8, 0.9, 0.99)))

# ------------------------------------------------------------- the marching arms
def churchill(re, eps):
    re = max(re, 1e-12)
    a1 = (7.0 / re) ** 0.9 + 0.27 * eps
    a = (2.457 * math.log(1.0 / a1)) ** 16
    b = (37530.0 / re) ** 16
    return 8.0 * ((8.0 / re) ** 12 + (a + b) ** (-1.5)) ** (1.0 / 12.0)


def march(name, q, p_up):
    """the solver's Euler recursion: dp/ds at the segment-start pressure"""
    g = GEOM[name]
    area = bore_area(g)
    mass = rho(p_up) * q
    ds = g["length"] / NSEG
    dz = (g["elevation"] / g["length"]) * ds
    p, out = p_up, (0.0, 0.0, 0.0)
    for _ in range(NSEG):
        r_ = rho(p)
        v = mass / (r_ * area)
        re = r_ * v * g["diameter"] / interp(MU_TAB, p)
        f = churchill(re, g["roughness"] / g["diameter"])
        out = (v, re, f)
        p -= r_ * G * dz + (f * ds / g["diameter"]) * r_ * v * v / 2.0
    return p, out


worst_close = 0.0
for name in MARCHES:
    row = ROW[name]
    a, b = GEOM[name]["src"], GEOM[name]["dst"]
    sign = 1.0
    if NODE_P[a] < NODE_P[b]:
        a, b, sign = b, a, -1.0
    q = abs(row["Q"])
    if q == 0.0:
        continue
    p_end, (v_out, re_out, f_out) = march(name, q, NODE_P[a])
    worst_close = max(worst_close, rel(p_end, NODE_P[b]))
    check(f"{name}: re-marched from {a} at the reported rate it lands on {b}",
          rel(p_end, NODE_P[b]) < 1e-9,
          f"{p_end:.3f} against the reported {NODE_P[b]:.3f} Pa over a "
          f"{row['dp']:.6g} Pa drop")
    check(f"{name}: and reports the march's own outlet state",
          near(row["Re"], re_out, 1e-8) and near(row["f"], f_out, 1e-8)
          and near(row["v"], sign * v_out, 1e-8),
          f"Re {row['Re']:.6g}, f {row['f']:.6g}, v {row['v']:.6g} "
          f"({'>' if sign < 0 else '<'} running backward against the declaration)")

# The network, not only the valves.
FREE = [n for n, b in NODES.items() if num(b, "pressure") is None]
check("the manifold is the only free node the valves hang on",
      "ANN" in FREE and all(GEOM[v]["src"] == "ANN" for v in VALVES),
      f"{len(FREE)} free nodes, {len(VALVES)} ports")
worst_node = 0.0
for n in FREE:
    terms = []
    for pname, g in GEOM.items():
        if g["dst"] == n:
            terms.append(ROW[pname]["Q"])
        elif g["src"] == n:
            terms.append(-ROW[pname]["Q"])
    for sbody in SOURCES.values():
        if word(sbody, "node") == n:
            terms.append(num(sbody, "rate"))
    scale = max(abs(t) for t in terms)
    imb = abs(sum(terms)) / max(scale, 1e-300)
    worst_node = max(worst_node, imb)
    check(f"node {n}: the reported rows balance", imb < 1e-9,
          f"{sum(terms):.3e} m3/s left over at a {scale:.3e} m3/s node")
check("every free node closes on the printed rows", worst_node < 1e-9,
      f"worst relative imbalance {worst_node:.3e}")

ports = sum(ROW[v]["Q"] for v in VALVES)
tie = -ROW["MFL"]["Q"]
gas = num(SOURCES["gas"], "rate")
check("the manifold is fed by the tie-in and the compressor, and by nothing else",
      near(ports, tie + gas, 1e-9),
      f"{ports:.6g} m3/s out of the ports against {tie:.6g} from the header "
      f"plus {gas:.6g} from the compressor")
check("the tie-in carries its gas against the direction it is declared in",
      ROW["MFL"]["Q"] < 0 and NODE_P["ANN"] < NODE_P["SUPPLY"],
      f"manifold {NODE_P['ANN']:.1f} Pa under the {NODE_P['SUPPLY']:.1f} Pa header")
check("every shoe is a boundary, printed as itself",
      all(near(NODE_P[n], num(NODES[n], "pressure"), 1e-15)
          for n in NODES if num(NODES[n], "pressure") is not None),
      ", ".join(f"{n} {NODE_P[n]:.0f}" for n in NODES
                if num(NODES[n], "pressure") is not None))

# ------------------------------------------------- declaration scope: refusals
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


SUB_PORT = "        length 0; diameter 0.0047625;   // a 3/16 in port"
SUB_SET = "liftDomePressure 7.4e6; liftAreaRatio 0.10;"
SUB_NOZZLE = "        liftDischargeCoeff 0.665; liftKappa 1.3;"
SUB_SPRING = "liftSpringPressure 0.5e6;"


def block_span(body, name):
    """the character range of one named entry, braces included"""
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


for anchor, why in ((SUB_PORT, "VLV_SUB's port line"), (SUB_SET, "VLV_SUB's dome"),
                    (SUB_NOZZLE, "VLV_SUB's nozzle")):
    lo, hi = block_span(src, "VLV_SUB")
    if anchor not in src[lo:hi]:
        raise SystemExit(f"probe: {why} is no longer in the VLV_SUB block")


def patched(replace=None, extra="", scope="VLV_SUB"):
    """rewrite once inside one pipe block, so a shared knob stays arm-specific"""
    lo, hi = block_span(src, scope)
    head = src[lo:hi]
    if extra:
        new = head.replace(SUB_PORT, SUB_PORT + "\n        " + extra, 1)
        if new == head:
            raise SystemExit("probe: the port anchor is not in the scope")
        head = new
    if replace:
        old, new = replace
        cut = head.replace(old, new, 1)
        if cut == head:
            raise SystemExit(f"probe: {old} is not in the {scope} scope")
        head = cut
    return src[:lo] + head + src[hi:]


refuses("a valve with a length",
        patched((SUB_PORT, "        length 5; diameter 0.0047625;  // a 3/16 in port")),
        "gas-lift valve with a length")
refuses("a valve with an elevation", patched(extra="elevation 10;"),
        "gas-lift valve with an elevation")
refuses("a valve with no port",
        patched((SUB_PORT, "        length 0; diameter 0;   // a 3/16 in port")),
        "gas-lift valve without a port")
# The valve's own scope refuses anything that would price a wall loss on a point
# element. Two of these keys are refused one gate earlier, by the law that owns
# them on any single-phase pipe at all, so each probe asserts the message the
# user would actually get rather than the one this increment added.
for decl, want in (("roughness 1e-5;", "gas-lift valve that also declares"),
                   ("valve 0.5;", "gas-lift valve that also declares"),
                   ("pumpHead 10;", "gas-lift valve that also declares"),
                   ("dInner 0.001;", "gas-lift valve that also declares"),
                   ("ductHeight 0.001;", "gas-lift valve that also declares"),
                   ("yieldStress 10;", "gas-lift valve that also declares"),
                   ("coupledResistance 1e6;", "is a coupled element"),
                   ("flowPattern 2;", "on a single-phase network")):
    refuses(f"a valve that also declares {decl[:-1]}", patched(extra=decl), want)
refuses("a valve whose R is 0", patched((SUB_SET, SUB_SET.replace("0.10", "0.0"))),
        "lift-valve area ratio outside (0,1)")
refuses("a valve whose R is 1", patched((SUB_SET, SUB_SET.replace("0.10", "1.0"))),
        "lift-valve area ratio outside (0,1)")
refuses("a negative spring",
        patched((SUB_SPRING, SUB_SPRING.replace("0.5e6", "-0.5e6")), scope="VLV_SPR"),
        "negative lift-valve spring pressure")
refuses("a discharge coefficient above 1",
        patched((SUB_NOZZLE, SUB_NOZZLE.replace("0.665", "1.5"))),
        "discharge coefficient outside (0,1]")
refuses("a nozzle with no discharge coefficient declared",
        patched((SUB_NOZZLE, "        liftKappa 1.3;")),
        "discharge coefficient outside (0,1]")
refuses("kappa at 1", patched((SUB_NOZZLE, SUB_NOZZLE.replace("1.3", "1.0"))),
        "lift-valve kappa <= 1")
refuses("a negative dome charge",
        patched((SUB_SET, SUB_SET.replace("7.4e6", "-7.4e6"))),
        "negative liftDomePressure")
refuses("lift parameters with no dome to switch them on",
        patched((SUB_SET, "liftAreaRatio 0.10;")),
        "lift-valve parameter without liftDomePressure")
refuses("a zero-length pipe that is not a valve",
        src.replace("TBG1 { from T1; to WH_HI; length 200;",
                    "TBG1 { from T1; to WH_HI; length 0;", 1),
        "needs positive length and diameter")
FILL_END = text.index("nodes")
for fill, want in (("fluid { rho 45; mu 1.2e-5; }", "constant-property network"),
                   ("fluid { blackoil { oilDensity 850; } }", "black-oil network")):
    refuses(f"a valve on a {want.split()[0]} fill",
            text[:text.index("fluid")] + fill + "\n\n" + text[FILL_END:], want)

# the key on the rows that have no injection string to price
for exe, case_name, tag in ((TRANSIENT, "transientBetaMomentum", "transient: pipe"),
                            (TWOFLUID, "twoFluidPressureLaminar", "twoFluid: pipe")):
    if not exe:
        continue
    other = open(os.path.join(HERE, "cases", case_name, "case.dict")).read()
    i = other.index("diameter")
    j = other.index(";", i) + 1
    code, out = run_on(exe, other[:j] + " liftDomePressure 7.0e6;" + other[j:])
    check(f"refused: a lift valve on the {case_name} row",
          code != 0 and tag in out and "lift valve" in out,
          out.splitlines()[-1][:150] if out else "no message")
if NODAL:
    nd = open(os.path.join(HERE, "cases", "nodalVogel", "case.dict")).read()
    key = "    elevation 2200;    // m, wellhead above the pay"
    if key not in nd:
        raise SystemExit("probe: the nodalVogel well block moved")
    code, out = run_on(NODAL, nd.replace(key, key + "\n    liftDomePressure 7.0e6;"))
    check("refused: a lift valve on the nodal VLP well",
          code != 0 and "nodal: well declares a lift valve" in out,
          out.splitlines()[-1][:150] if out else "no message")

if fails:
    raise SystemExit(f"gas-lift valve verification FAILED: {len(fails)}: {fails}")
print(f"the valve is a force balance and a nozzle: Pvo/Pvc/spread/PPE/TRO "
      f"restated at the solved tubing pressure (worst "
      f"{max(err['pvo'], err['pvc'], err['spread'], err['ppe'], err['tro']):.2g} "
      f"rel), the isentropic port re-integrated off the table's own density "
      f"(worst {err['m']:.2g} rel), {len(VALVES)} arms giving all three "
      f"dispositions and both throat states, choked arms bit-identical across a "
      f"megapascal of back pressure, and the marching arms closed to "
      f"{worst_close:.2g} rel")
