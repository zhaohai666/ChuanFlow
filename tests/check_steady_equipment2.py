#!/usr/bin/env python3
"""P10 (5.156): independent verification of the separator and the compressor
as 1D steady point elements, and of the compressor as a network pipe key on
the tabulated-PVT march.

The checker re-derives every number from the case text alone, with no solver
internals:

  * the separator's flash allocation through the shipped Standing/VB/DAK
    restatement (port of check_pvt_separator.py), the per-component mass
    balance with entrainment, and the directional pair eff=1 vs eff=0.9;
  * both compressor modes' closed forms, including the anchor
    W = m Z R_spec T1 (k/(k-1)) (PR^((k-1)/k) - 1), the head-mode round
    trip H -> PR -> H, and the eta=1 vs eta=0.9 direction;
  * the network march on the same 24-segment grid with the same Churchill
    friction and the same clamped linear rho(p)/mu(p) lookups, the boost
    riding inside the same Illinois bracket, for the forward curve arm,
    the constant-ratio reverse arm and the bare twin;
  * the full rejection set on both grammars.

usage: check_steady_equipment2.py <pointExe> <pointCase> <steadyExe> <steadyCase>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 5:
    raise SystemExit(__doc__)

POINT_EXE, POINT_CASE, STEADY_EXE, STEADY_CASE = sys.argv[1:5]
POINT_CASE = os.path.abspath(POINT_CASE)
STEADY_CASE = os.path.abspath(STEADY_CASE)

fails = []
NCHECK = [0]


def check(name, cond, detail=""):
    NCHECK[0] += 1
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def near(got, want, rtol=1e-9, atol=1e-12):
    return abs(got - want) <= rtol * abs(want) + atol


def block_body(src, key):
    """Body of `key { ... }` at depth 0 of src, or None."""
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", src)
    if not m:
        return None
    start = src.index("{", m.start())
    depth = 0
    for j in range(start, len(src)):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                return src[start + 1:j]
    raise SystemExit(f"case: unbalanced {key}{{}} block")


def top_names(src):
    """Names of the `NAME { ... }` blocks at depth 0 of src."""
    names = []
    depth = 0
    for m in re.finditer(r"[{}]|([A-Za-z_][A-Za-z0-9_]*)\s*\{", src):
        tok = m.group(0)
        if tok == "{":
            depth += 1
        elif tok == "}":
            depth -= 1
        else:
            if depth == 0:
                names.append(m.group(1))
            depth += 1
    return names


def scalar(src, key, default=None):
    m = re.search(r"\b" + re.escape(key) + r"\s+([-+0-9.eE]+)\s*;", src)
    if m:
        return float(m.group(1))
    if default is None:
        raise SystemExit(f"case has no scalar '{key}'")
    return default


def numlist(src, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\(([^)]*)\)", src)
    if not m:
        return None
    return [float(v) for v in m.group(1).replace(";", " ").split()]


# ---- physical constants, transcribed from src/pvt/src/blackOil.cpp -------
PSI = 6894.757293168361
SCF_STB = 5.614583333333333
M_AIR = 28.9647e-3
R_GAS = 8.314462618
P_STD, T_STD = 101325.0, 288.15
RHO_WATER_STD = 999.014


def standing_pb(gg, api_, tf, rs_scf):
    a = 0.00091 * tf - 0.0125 * api_
    return max(18.2 * ((rs_scf / gg) ** 0.83 * 10.0 ** a - 1.4), 0.0)


def vb_rs(p_psia, gg, api_, tf):
    if api_ <= 30.0:
        c1, c2, c3 = 0.0362, 1.0937, 25.724
    else:
        c1, c2, c3 = 0.0178, 1.187, 23.931
    tr = tf + 459.67
    return c1 * gg * max(p_psia, 0.0) ** c2 * math.exp(c3 * api_ / tr)


def vb_bo(p_psia, rs_scf, pb_psia, gg, api_, tf):
    heavy = api_ > 30.0
    c1 = 4.670e-4 if heavy else 4.677e-4
    c2 = 1.100e-5 if heavy else 1.751e-5
    c3 = 1.337e-9 if heavy else -1.811e-8
    t_term = (tf - 60.0) * api_ / gg
    bob = 1.0 + c1 * rs_scf + c2 * t_term + c3 * rs_scf * t_term
    if p_psia <= pb_psia:
        return bob
    co = (-1433.0 + 5.0 * rs_scf + 17.2 * tf - 1180.0 * gg + 12.61 * api_) \
         / (1.0e5 * p_psia)
    return bob * math.exp(co * (pb_psia - p_psia))


def standing_pc(gg):
    if gg < 0.75:
        return (677.0 + 15.0 * gg - 37.5 * gg * gg,
                168.0 + 325.0 * gg - 12.5 * gg * gg)
    return (706.0 + 51.0 * gg - 11.1 * gg * gg,
            187.0 + 330.0 * gg - 71.5 * gg * gg)


def dak_z(p_psia, tr_rankine, ppc, tpc):
    A = (0.3265, -1.0700, -0.5339, 0.01569, -0.05165,
         0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210)
    pr = max(p_psia, 1e-6) / ppc
    tr = tr_rankine / tpc
    t1 = (A[0] + A[1] / tr + A[2] / tr ** 3 + A[3] / tr ** 4 + A[4] / tr ** 5)
    t2 = A[5] + A[6] / tr + A[7] / (tr * tr)

    def z_of_rho(rho):
        r2 = rho * rho
        return (1.0 + t1 * rho + t2 * r2
                - A[8] * (A[6] / tr + A[7] / (tr * tr)) * rho ** 5
                + A[9] * (1.0 + A[10] * r2) * (r2 / tr ** 3)
                  * math.exp(-A[10] * r2))

    def residual(z):
        return z_of_rho(0.27 * pr / (z * tr)) - z

    lo, hi = 0.05, 3.0
    if residual(lo) <= 0.0:
        return lo
    if residual(hi) >= 0.0:
        return hi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if residual(mid) > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-13 * hi:
            break
    return 0.5 * (lo + hi)


def flash(gg, gor, t_k, p_pa):
    """rs, bo, bg, z at (p, T) for a stream whose std GOR is gor."""
    gamma_o = scalar_oil_density[0] / RHO_WATER_STD
    api_ = 141.5 / gamma_o - 131.5
    p_psia = max(p_pa, 1.0) / PSI
    tf = t_k * 1.8 - 459.67
    tr = t_k * 1.8
    gor_scf = gor * SCF_STB
    pb_psia = standing_pb(gg, api_, tf, gor_scf)
    if p_pa <= pb_psia * PSI:
        rs = min(vb_rs(p_psia, gg, api_, tf) / SCF_STB, gor)
        bo = vb_bo(p_psia, rs * SCF_STB, pb_psia, gg, api_, tf)
    else:
        rs = gor
        bo = vb_bo(p_psia, gor_scf, pb_psia, gg, api_, tf)
    ppc, tpc = standing_pc(gg)
    z = dak_z(p_psia, tr, ppc, tpc)
    bg = (P_STD / p_pa) * (t_k / T_STD) * z
    return rs, bo, bg, z


# ============================================================================
# Part A: separators
# ============================================================================
point_text = re.sub(r"//[^\n]*", "", open(POINT_CASE).read())
seps_block = block_body(point_text, "separators")
if seps_block is None:
    raise SystemExit("case has no separators{} block")

scalar_oil_density = [0.0]  # flash() reads it per arm below

sep_arms = {}
for name in top_names(seps_block):
    body = block_body(seps_block, name)
    bd = block_body(body, "blackoil") or ""
    sep_arms[name] = {
        "feedRate": scalar(body, "feedRate"),
        "pressure": scalar(body, "pressure"),
        "temperature": scalar(body, "temperature"),
        "efficiency": scalar(body, "efficiency", 1.0),
        "dP": scalar(body, "dP", 0.0),
        "gasGravity": scalar(bd, "gasGravity", 0.65),
        "oilDensity": scalar(bd, "oilDensity", 850.0),
        "gor": scalar(bd, "gor", 0.0),
        "waterCut": scalar(bd, "waterCut", 0.0),
        "waterDensity": scalar(bd, "waterDensity", 1000.0),
    }

proc = subprocess.run([POINT_EXE, POINT_CASE], capture_output=True, text=True)
if proc.returncode != 0:
    raise SystemExit(f"point exe failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")


def parse_kv(line, prefix):
    tok = line.replace("[chflow]", "").split()
    if not tok or tok[0] != prefix:
        return None
    d = {"name": tok[1]}
    for a, b in zip(tok[2::2], tok[3::2]):
        try:
            d[a] = float(b)
        except ValueError:
            d[a] = b
    return d


sep_lines = {}
cmp_lines = {}
for line in proc.stdout.splitlines():
    s = parse_kv(line, "sep")
    if s:
        sep_lines[s["name"]] = s
    c = parse_kv(line, "cmp")
    if c:
        cmp_lines[c["name"]] = c

check("all separator arms printed", len(sep_lines) == len(sep_arms),
      f"got {len(sep_lines)} want {len(sep_arms)}")

for name, arm in sep_arms.items():
    if name not in sep_lines:
        check(f"sep {name} printed", False)
        continue
    got = sep_lines[name]
    GG = arm["gasGravity"]
    RHO_STD = arm["oilDensity"]
    GOR = arm["gor"]
    WC = arm["waterCut"]
    WATER_DENSITY = arm["waterDensity"]
    scalar_oil_density[0] = RHO_STD
    wor = WC / (1.0 - WC) if 0.0 < WC < 1.0 else 0.0
    rho_gas_std = P_STD * GG * M_AIR / (R_GAS * T_STD)
    rs, bo, bg, z = flash(GG, GOR, arm["temperature"], arm["pressure"])
    glib = GOR - rs
    q = arm["feedRate"]
    mf = q * (RHO_STD + GOR * rho_gas_std + wor * WATER_DENSITY)
    mli = q * (RHO_STD + rs * rho_gas_std + wor * WATER_DENSITY)
    ment = (1.0 - arm["efficiency"]) * mli
    ml = arm["efficiency"] * mli
    mg = q * glib * rho_gas_std + ment

    check(f"sep {name}: stage state and flash allocation match independently",
          near(got.get("p", -1), arm["pressure"])
          and near(got.get("T", -1), arm["temperature"])
          and near(got.get("eff", -1), arm["efficiency"])
          and near(got.get("pout", -1), arm["pressure"] - arm["dP"])
          and near(got.get("rs", -1), rs)
          and near(got.get("glib", -1), glib)
          and near(got.get("bo", -1), bo)
          and near(got.get("bg", -1), bg)
          and near(got.get("z", -1), z)
          and near(got.get("rhogs", -1), rho_gas_std),
          f"rs got {got.get('rs')} want {rs}")
    check(f"sep {name}: per-component mass balance closes",
          near(got.get("mf", -1), mf)
          and near(got.get("mli", -1), mli)
          and near(got.get("ment", -1), ment)
          and near(got.get("ml", -1), ml)
          and near(got.get("mg", -1), mg)
          and abs(got.get("closure", 1e9)) <= 1e-9 * mf
          and glib >= -1e-12,
          f"closure={got.get('closure')}")

if "SEP_IDEAL" in sep_lines and "SEP_ENTRAIN" in sep_lines:
    i = sep_lines["SEP_IDEAL"]
    e = sep_lines["SEP_ENTRAIN"]
    check("directional: eff=0.9 loses liquid mass and gains it on the gas "
          "outlet, entrainment positive only for eff < 1",
          e.get("ml", 0) < i.get("ml", 0)
          and e.get("mg", 0) > i.get("mg", 0)
          and e.get("ment", 0) > 0
          and i.get("ment", -1) == 0)
if "SEP_DEAD" in sep_lines:
    d = sep_lines["SEP_DEAD"]
    check("dead oil: nothing liberates, the gas outlet is entrainment alone "
          "and the liquid outlet carries the whole feed",
          near(d.get("glib", -1), 0.0, 1e-12)
          and near(d.get("mg", -1), d.get("ment", -1), 1e-9)
          and near(d.get("ml", -1), d.get("mf", -1), 1e-9))


# ============================================================================
# Part B: compressors
# ============================================================================
cmps_block = block_body(point_text, "compressors")
if cmps_block is None:
    raise SystemExit("case has no compressors{} block")

cmp_arms = {}
for name in top_names(cmps_block):
    body = block_body(cmps_block, name)
    cmp_arms[name] = {
        "pressureRatio": scalar(body, "pressureRatio", 0.0),
        "shutoffRatio": scalar(body, "shutoffRatio", 0.0),
        "ratedRate": scalar(body, "ratedRate", 0.0),
        "polytropicHead": scalar(body, "polytropicHead", 0.0),
        "polytropicEfficiency": scalar(body, "polytropicEfficiency", 1.0),
        "kappa": scalar(body, "kappa", 0.0),
        "specificGasConstant": scalar(body, "specificGasConstant", 0.0),
        "compressibility": scalar(body, "compressibility", 1.0),
        "suctionPressure": scalar(body, "suctionPressure"),
        "suctionTemperature": scalar(body, "suctionTemperature"),
        "rate": scalar(body, "rate"),
    }

check("all compressor arms printed", len(cmp_lines) == len(cmp_arms),
      f"got {len(cmp_lines)} want {len(cmp_arms)}")

for name, a in cmp_arms.items():
    if name not in cmp_lines:
        check(f"cmp {name} printed", False)
        continue
    got = cmp_lines[name]
    ex = (a["kappa"] - 1.0) / a["kappa"]
    ZR = a["compressibility"] * a["specificGasConstant"]
    if a["pressureRatio"] > 0.0:
        shutoff = a["shutoffRatio"] if a["shutoffRatio"] > 0 \
            else a["pressureRatio"]
        x = a["rate"] / a["ratedRate"]
        pr = shutoff - (shutoff - a["pressureRatio"]) * x * x
        mode = "curve"
        h = ZR * a["suctionTemperature"] / ex * (pr ** ex - 1.0)
    else:
        h = a["polytropicHead"]
        pr = (1.0 + h * ex / (ZR * a["suctionTemperature"])) ** (1.0 / ex)
        mode = "head"
    p2 = a["suctionPressure"] * pr
    mdot = a["suctionPressure"] * a["rate"] / (ZR * a["suctionTemperature"])
    w0 = mdot * h
    w = w0 / a["polytropicEfficiency"]
    t2 = a["suctionTemperature"] * pr ** (ex / a["polytropicEfficiency"])

    ok_state = (got.get("mode") == mode
                and near(got.get("p1", -1), a["suctionPressure"])
                and near(got.get("t1", -1), a["suctionTemperature"])
                and near(got.get("q", -1), a["rate"])
                and near(got.get("pr", -1), pr)
                and near(got.get("p2", -1), p2)
                and near(got.get("m", -1), mdot))
    ok_power = (near(got.get("h", -1), h)
                and near(got.get("w0", -1), w0)
                and near(got.get("w", -1), w)
                and near(got.get("t2", -1), t2))
    check(f"cmp {name}: mode, suction state, PR, p2 and mass rate match",
          ok_state, f"pr got {got.get('pr')} want {pr}")
    check(f"cmp {name}: W0 = m Z R_spec T1 (k/(k-1)) (PR^((k-1)/k) - 1), "
          f"W = W0/eta, T2 = T1 PR^((k-1)/(k eta))", ok_power)
    if a["pressureRatio"] <= 0.0:
        # head-mode round trip: the anchor re-priced at the derived PR must
        # return the declared head
        h_back = ZR * a["suctionTemperature"] / ex * (pr ** ex - 1.0)
        check(f"cmp {name}: head round trip H -> PR -> H",
              near(h_back, a["polytropicHead"], 1e-9))

if "CMP_HEAD_IDEAL" in cmp_lines and "CMP_HEAD_LOSS" in cmp_lines:
    i = cmp_lines["CMP_HEAD_IDEAL"]
    l = cmp_lines["CMP_HEAD_LOSS"]
    check("directional: eta=1 shaft power equals the gas power, eta=0.9 "
          "costs more and discharges hotter at an identical declared-head "
          "ratio",
          near(i.get("w", -1), i.get("w0", -1), 1e-12)
          and l.get("w", 0) > i.get("w", 0)
          and l.get("t2", 0) > i.get("t2", 0)
          and near(l.get("pr", -1), i.get("pr", -1), 1e-12))
if "CMP_CURVE_RATED" in cmp_lines and "CMP_CURVE_HALF" in cmp_lines:
    check("curve: rated point reads pressureRatio exactly, half rate reads "
          "2.4 - 0.4/4 = 2.3, lower rate boosts harder",
          near(cmp_lines["CMP_CURVE_RATED"].get("pr", -1), 2.0, 1e-12)
          and near(cmp_lines["CMP_CURVE_HALF"].get("pr", -1), 2.3, 1e-12)
          and cmp_lines["CMP_CURVE_HALF"].get("pr", 0)
              > cmp_lines["CMP_CURVE_RATED"].get("pr", 0))


# ============================================================================
# Part C: network march re-derivation (port of marchPvtPipe/pvtRateForDrop)
# ============================================================================
steady_text = re.sub(r"//[^\n]*", "", open(STEADY_CASE).read())
pvt_block = block_body(block_body(steady_text, "fluid"), "pvt")
p_table = numlist(pvt_block, "pressure")
rho_table = numlist(pvt_block, "rho")
mu_table = numlist(pvt_block, "mu")
T_FLUID = scalar(pvt_block, "temperature")

nodes = {}
nb = block_body(steady_text, "nodes")
for n in top_names(nb):
    nodes[n] = scalar(block_body(nb, n), "pressure")

pipes = {}
pb = block_body(steady_text, "pipes")
for n in top_names(pb):
    body = block_body(pb, n)
    pipes[n] = {
        "from": re.search(r"\bfrom\s+([A-Za-z_][A-Za-z0-9_]*)\s*;", body).group(1),
        "to": re.search(r"\bto\s+([A-Za-z_][A-Za-z0-9_]*)\s*;", body).group(1),
        "length": scalar(body, "length"),
        "diameter": scalar(body, "diameter"),
        "roughness": scalar(body, "roughness", 0.0),
        "elevation": scalar(body, "elevation", 0.0),
        "pressureRatio": scalar(body, "pressureRatio", 0.0),
        "shutoffRatio": scalar(body, "shutoffRatio", 0.0),
        "ratedRate": scalar(body, "ratedRate", 0.0),
        "polytropicEfficiency": scalar(body, "polytropicEfficiency", 1.0),
        "kappa": scalar(body, "kappa", 0.0),
        "specificGasConstant": scalar(body, "specificGasConstant", 0.0),
        "compressibility": scalar(body, "compressibility", 1.0),
    }


def lookup(p, tp, tv):
    if p <= tp[0]:
        return tv[0]
    if p >= tp[-1]:
        return tv[-1]
    for i in range(len(tp) - 1):
        if p <= tp[i + 1]:
            t = (p - tp[i]) / (tp[i + 1] - tp[i])
            return tv[i] + t * (tv[i + 1] - tv[i])
    return tv[-1]


def churchill(re, eps):
    re = max(re, 1e-12)
    a1 = (7.0 / re) ** 0.9 + 0.27 * eps
    a = (2.457 * math.log(1.0 / a1)) ** 16.0
    b = (37530.0 / re) ** 16.0
    return 8.0 * ((8.0 / re) ** 12.0 + (a + b) ** -1.5) ** (1.0 / 12.0)


def march(pipe, q, p_up, sin_flow, p_start=0.0, nseg=24):
    area = 0.25 * math.pi * pipe["diameter"] ** 2
    mass_rate = lookup(p_up, p_table, rho_table) * q
    ds = pipe["length"] / nseg
    dz = sin_flow * ds
    p_floor = max(0.02 * p_up, 1000.0)
    p = p_start if p_start > 0.0 else p_up
    for _ in range(nseg):
        r = lookup(p, p_table, rho_table)
        m = lookup(p, p_table, mu_table)
        v = mass_rate / (r * area)
        f = churchill(r * v * pipe["diameter"] / m,
                      pipe["roughness"] / pipe["diameter"])
        dp = r * 9.80665 * dz + f * ds / pipe["diameter"] * r * v * v / 2.0
        p -= dp
        if not p > p_floor:
            return p_floor, False
    return p, True


def pr_of(pipe, q):
    if pipe["pressureRatio"] <= 0.0:
        return 1.0
    shutoff = pipe["shutoffRatio"] if pipe["shutoffRatio"] > 0 \
        else pipe["pressureRatio"]
    x = q / pipe["ratedRate"]
    return shutoff - (shutoff - pipe["pressureRatio"]) * x * x


def rate_for_drop(pipe, p_up, p_low, sin_flow, boost_up):
    p_floor = max(0.02 * p_up, 1000.0)

    def g(q):
        pr = pr_of(pipe, q)
        start = p_up * pr if boost_up else 0.0
        target = p_low if boost_up else p_low * pr
        po, ok = march(pipe, q, p_up, sin_flow, start)
        return (po if ok else p_floor) - target

    area = 0.25 * math.pi * pipe["diameter"] ** 2
    q_seed = max(1e-3 * area, 1e-8)
    g0 = g(0.0)
    if g0 <= 0.0:
        return 0.0
    q_lo, f_lo = 0.0, g0
    q_hi, f_hi = q_seed, g(q_seed)
    for _ in range(60):
        if f_hi <= 0.0:
            break
        q_lo, f_lo = q_hi, f_hi
        q_hi *= 4.0
        f_hi = g(q_hi)
    q = q_lo + (q_hi - q_lo) * f_lo / (f_lo - f_hi)
    for _ in range(100):
        fq = g(q)
        if abs(fq) < 1e-12 * p_up:
            break
        if fq > 0.0:
            q_lo, f_lo = q, fq
            f_lo *= 0.5
        else:
            q_hi, f_hi = q, fq
            f_hi *= 0.5
        q_new = q_lo + (q_hi - q_lo) * f_lo / (f_lo - f_hi)
        q = q_new if q_lo < q_new < q_hi else 0.5 * (q_lo + q_hi)
        if abs(q_hi - q_lo) < 1e-12 * max(q_hi, 1e-30):
            break
    return q


proc2 = subprocess.run([STEADY_EXE, STEADY_CASE], capture_output=True,
                       text=True)
if proc2.returncode != 0:
    raise SystemExit(f"steady exe failed rc={proc2.returncode}\n"
                     f"{proc2.stdout}\n{proc2.stderr}")

pipe_lines = {}
for line in proc2.stdout.splitlines():
    d = parse_kv(line, "pipe")
    if d:
        pipe_lines[d["name"]] = d

check("all pipes printed", len(pipe_lines) == len(pipes),
      f"got {len(pipe_lines)} want {len(pipes)}")

for name, pipe in pipes.items():
    if name not in pipe_lines:
        check(f"pipe {name} printed", False)
        continue
    got = pipe_lines[name]
    p_from = nodes[pipe["from"]]
    p_to = nodes[pipe["to"]]
    pr0 = pr_of(pipe, 0.0)
    forward = p_from * pr0 >= p_to
    p_up = p_from if forward else p_to
    p_low = p_to if forward else p_from
    boost_up = pipe["pressureRatio"] > 0.0 and forward
    sin_flow = (1.0 if forward else -1.0) * pipe["elevation"] / pipe["length"]
    q = rate_for_drop(pipe, p_up, p_low, sin_flow, boost_up)
    q_signed = q if forward else -q
    check(f"pipe {name}: Q matches the independent boost-aware march and dp "
          "echoes the fixed node pair",
          near(got.get("Q", 1e9), q_signed, 1e-6)
          and near(got.get("dp", 1e9), p_from - p_to, 1e-9)
          and ("cmpRatio" not in got) == (pipe["pressureRatio"] <= 0.0),
          f"Q got {got.get('Q')} want {q_signed}")
    if pipe["pressureRatio"] > 0.0:
        pr = pr_of(pipe, q)
        ex = (pipe["kappa"] - 1.0) / pipe["kappa"]
        h = pipe["compressibility"] * pipe["specificGasConstant"] \
          * T_FLUID / ex * (pr ** ex - 1.0)
        w = lookup(p_up, p_table, rho_table) * q * h \
          / pipe["polytropicEfficiency"]
        t2 = T_FLUID * pr ** (ex / pipe["polytropicEfficiency"])
        po, okm = march(pipe, q, p_up, sin_flow,
                        p_up * pr if boost_up else 0.0)
        target = p_low if boost_up else p_low * pr
        check(f"pipe {name}: curve ratio, boost, m H / eta power, passive "
              "T2, and the boosted march lands on the far end",
              near(got.get("cmpRatio", -1), pr)
              and near(got.get("cmpBoost", -1), p_from * (pr - 1.0))
              and near(got.get("cmpPower", -1), w, 1e-6)
              and near(got.get("cmpT2", -1), t2, 1e-9)
              and okm and abs(po - target) <= 1e-6 * p_up,
              f"landing residual {po - target}")

if "T" in pipe_lines and "K" in pipe_lines and "R" in pipe_lines:
    check("directional: the bare twin runs backwards 9 -> 5 MPa, the boosted "
          "line runs forwards against the same gradient, and the "
          "constant-ratio machine with an insufficient boost still runs "
          "backwards",
          pipe_lines["T"].get("Q", 0) < 0
          and pipe_lines["K"].get("Q", 0) > 0
          and pipe_lines["R"].get("Q", 0) < 0)
node_ok = True
for nm in ("A", "C", "Y"):
    m = re.search(r"node " + nm + r" p ([-+0-9.eE]+)", proc2.stdout)
    node_ok = node_ok and m is not None \
        and near(float(m.group(1)), nodes[nm], 1e-12)
check("every fixed node keeps its declared pressure", node_ok)


# ============================================================================
# Part D: rejections
# ============================================================================
def probe(exe, text, needle, label):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "probe.dict")
        open(path, "w").write(text)
        rr = subprocess.run([exe, path], capture_output=True, text=True)
        ok = rr.returncode != 0 and needle in rr.stdout + rr.stderr
        check(f"rejecting {label}", ok, f"needle '{needle}'")


def sub1(pattern, repl, text):
    out = re.sub(pattern, repl, text, count=1)
    if out == text:
        raise SystemExit(f"probe did not mutate the case: {pattern}")
    return out


# separator refusals
probe(POINT_EXE,
      sub1(r"\befficiency\s+[-+0-9.eE]+\s*;", "efficiency 1.2;", point_text),
      "efficiency", "sep efficiency above 1")
probe(POINT_EXE,
      sub1(r"\befficiency\s+[-+0-9.eE]+\s*;", "efficiency -0.1;", point_text),
      "efficiency", "sep efficiency below 0")
probe(POINT_EXE,
      sub1(r"\bfeedRate\s+[-+0-9.eE]+\s*;", "feedRate 0.0;", point_text),
      "feed rate", "sep zero feed rate")
probe(POINT_EXE,
      sub1(r"\bdP\s+[-+0-9.eE]+\s*;", "dP -1.0;", point_text),
      "dP", "sep negative dP")
probe(POINT_EXE,
      sub1(r"\bpressure\s+5\.0e5\s*;", "pressure -1.0;", point_text),
      "absolute pressure", "sep negative stage pressure")
probe(POINT_EXE,
      sub1(r"\btemperature\s+320\.0\s*;", "temperature -5.0;", point_text),
      "kelvin", "sep negative stage temperature")
probe(POINT_EXE,
      sub1(r"\bwaterCut\s+0\.1\s*;", "waterCut 1.2;", point_text),
      "waterCut", "sep waterCut above 1")

# compressor point refusals
probe(POINT_EXE,
      sub1(r"\bpolytropicEfficiency\s+[-+0-9.eE]+\s*;",
           "polytropicEfficiency 0.0;", point_text),
      "efficiency", "cmp eta 0")
probe(POINT_EXE,
      sub1(r"\bpolytropicEfficiency\s+[-+0-9.eE]+\s*;",
           "polytropicEfficiency 1.5;", point_text),
      "efficiency", "cmp eta above 1")
probe(POINT_EXE,
      sub1(r"\bkappa\s+[-+0-9.eE]+\s*;", "kappa 1.0;", point_text),
      "kappa", "cmp kappa 1")
probe(POINT_EXE,
      sub1(r"\brate\s+0\.05\s*;", "rate 0.0;", point_text),
      "rate", "cmp zero rate")
probe(POINT_EXE,
      sub1(r"\bsuctionPressure\s+[-+0-9.eE]+\s*;",
           "suctionPressure -1.0;", point_text),
      "suction pressure", "cmp negative suction")
probe(POINT_EXE,
      sub1(r"\bpressureRatio 2\.0;", "pressureRatio 2.0; polytropicHead 1.0e5;",
           point_text),
      "both", "cmp curve and head together")
probe(POINT_EXE,
      sub1(r"\bpressureRatio 2\.0;", "", point_text),
      "neither", "cmp curve loses its ratio")
probe(POINT_EXE,
      sub1(r"\bratedRate 0\.04;", "ratedRate 0.0;", point_text),
      "ratedRate", "cmp curve without a rated point")
probe(POINT_EXE,
      sub1(r"\bshutoffRatio 2\.4;", "shutoffRatio 1.5;", point_text),
      "rises with rate", "cmp curve shutoff below rated")
probe(POINT_EXE,
      sub1(r"\bspecificGasConstant 518\.28;", "specificGasConstant 0.0;",
           point_text),
      "specific gas constant", "cmp zero R_spec")

# network pipe-key refusals
probe(STEADY_EXE,
      sub1(r"\bshutoffRatio 2\.4;", "shutoffRatio 1.5;", steady_text),
      "rises with rate", "net curve shutoff below rated")
probe(STEADY_EXE,
      sub1(r"\bpolytropicEfficiency 0\.9;", "polytropicEfficiency 1.5;",
           steady_text),
      "polytropicEfficiency", "net eta above 1")
probe(STEADY_EXE,
      sub1(r"\bkappa 1\.3;", "kappa 0.9;", steady_text),
      "kappa", "net kappa below 1")
probe(STEADY_EXE,
      sub1(r"\bpressureRatio 2\.0;", "pressureRatio 0.5;", steady_text),
      "unity", "net ratio below 1")
probe(STEADY_EXE,
      sub1(r"\bpressureRatio 2\.0;", "pressureRatio -1.0;", steady_text),
      "negative pressureRatio", "net negative ratio")
probe(STEADY_EXE,
      sub1(r"\bpressureRatio 2\.0; shutoffRatio 2\.4;",
           "valve 100; pressureRatio 2.0; shutoffRatio 2.4;", steady_text),
      "one machine per line", "net valve clash on the compressor line")
probe(STEADY_EXE,
      sub1(r"\bpressureRatio 2\.0; shutoffRatio 2\.4;",
           "jouleThomson 1.0e-6; pressureRatio 2.0; shutoffRatio 2.4;",
           steady_text),
      "temperature readout", "net thermal clash on the compressor line")
probe(STEADY_EXE,
      sub1(r"\bpressureRatio 2\.0; shutoffRatio 2\.4;", "", steady_text),
      "without pressureRatio", "net compressor parameters without the switch")
probe(STEADY_EXE,
      sub1(r"\bratedRate 0\.05;", "ratedRate 0.0;", steady_text),
      "ratedRate", "net curve without a rated rate")
probe(STEADY_EXE,
      sub1(r"pvt\s*\{[^}]*\}", "rho 850.0; mu 0.1;", steady_text),
      "tabulated", "net compressor on a constant-property fill")
probe(STEADY_EXE,
      sub1(r"\btemperature 300\.0;", "", steady_text),
      "temperature", "net compressor without the table temperature")


# ============================================================================
if fails:
    raise SystemExit(f"steadyEquipment2 FAILED ({len(fails)}/{NCHECK[0]}): "
                     f"{fails}")
print(f"steadyEquipment2 OK: {NCHECK[0]} checks passed -- separator flash "
      "allocation, compressor anchors, boosted network march, all refusals")
