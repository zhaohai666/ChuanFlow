#!/usr/bin/env python3
"""M6: independent verification of the black-oil separator cascade.

The C++ model reuses evaluateBlackOil at every stage; this checker re-derives
the whole physics from the case dictionary alone, with no solver internals:

  Standing bubble point, Vazquez-Beggs Rs and Bo, Standing pseudo-criticals,
  Dranchuk-Abou-Kassem Z (bisection), the rhoOil mass balance of one STB, the
  standard gas density, and the stage-by-stage cascade (each stage re-flashes
  the liquid actually arriving, i.e. its remaining solution GOR).

It asserts for every stage:

  * rs, bo, bg, z, rhoOil-derived liveApi, gasMass and liquidMass match the
    independent restatement;
  * gas_i = rs_in - rs_out >= 0 and rs is non-increasing across the cascade;
  * liquidMass_i + cumGasMass_i == blackOilMassPerOilSm3(fluid) (mass closure);
  * total gas = gor - final rs, massIn == massOut, reported imbalance tiny;

plus parse-time refusals, including the gated twin with the separator{} block
deleted (must be refused, not silently defaulted).

usage: check_pvt_separator.py <chflow-separator> <case.dict>
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
case = os.path.abspath(case)
text = re.sub(r"//[^\n]*", "", open(case).read())

fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def near(got, want, rtol=1e-9, atol=1e-12):
    return abs(got - want) <= rtol * abs(want) + atol


def get_block(src, key):
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
LBFT3_TO_SI = 16.0184634

sep = get_block(text, "separator")
if sep is None:
    raise SystemExit("shipped case has no separator{} block")
pressures = numlist(sep, "pressure")
if not pressures:
    raise SystemExit("case has no pressure (...) list")
temperatures = numlist(sep, "temperatures")
if temperatures is None:
    temperatures = [scalar(sep, "temperature")] * len(pressures)
if len(temperatures) != len(pressures):
    raise SystemExit("temperatures must give one value per stage")

fl = get_block(get_block(sep, "fluid") or "", "blackoil")
if fl is None:
    raise SystemExit("case has no fluid{ blackoil{...} } block")
GG = scalar(fl, "gasGravity", 0.65)
RHO_STD = scalar(fl, "oilDensity", 850.0)
GOR = scalar(fl, "gor", 0.0)
WC = scalar(fl, "waterCut", 0.0)
WATER_DENSITY = scalar(fl, "waterDensity", 1000.0)
gamma_o = RHO_STD / RHO_WATER_STD
api = 141.5 / gamma_o - 131.5
wor = WC / (1.0 - WC) if 0.0 < WC < 1.0 else 0.0
rho_gas_std = P_STD * GG * M_AIR / (R_GAS * T_STD)
MASS_IN = RHO_STD + GOR * rho_gas_std + wor * WATER_DENSITY
STOCK_API = 141.5 / gamma_o - 131.5


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
         / (1.0e5 * max(pb_psia, 300.0))
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


def evaluate(fluid_gor, t_k, p_pa):
    p_psia = max(p_pa, 1.0) / PSI
    tf = t_k * 1.8 - 459.67
    tr = t_k * 1.8
    gor_scf = fluid_gor * SCF_STB
    pb_psia = standing_pb(GG, api, tf, gor_scf)
    if p_pa <= pb_psia * PSI:
        rs = min(vb_rs(p_psia, GG, api, tf) / SCF_STB, fluid_gor)
        bo = vb_bo(p_psia, rs * SCF_STB, pb_psia, GG, api, tf)
    else:
        rs = fluid_gor
        bo = vb_bo(p_psia, gor_scf, pb_psia, GG, api, tf)
    ppc, tpc = standing_pc(GG)
    z = dak_z(p_psia, tr, ppc, tpc)
    bg = (P_STD / p_pa) * (t_k / T_STD) * z
    rho_oil = (350.0 * gamma_o + 0.0764 * GG * rs * SCF_STB) \
              / (SCF_STB * bo) * LBFT3_TO_SI
    return rs, bo, bg, z, rho_oil


# ---- independent cascade ------------------------------------------------
expected = []
available = GOR
cum_gas = 0.0
total_gas = 0.0
for p, t in zip(pressures, temperatures):
    rs, bo, bg, z, rho_oil = evaluate(available, t, p)
    gas_lib = available - rs
    gas_mass = gas_lib * rho_gas_std
    cum_gas += gas_mass
    total_gas += gas_lib
    liquid_mass = RHO_STD + rs * rho_gas_std + wor * WATER_DENSITY
    live_api = 141.5 / (rho_oil / RHO_WATER_STD) - 131.5
    expected.append({
        "p": p, "T": t, "available": available, "rs": rs,
        "gasSm3": gas_lib, "gasMass": gas_mass, "liquidMass": liquid_mass,
        "cumGasMass": cum_gas, "bo": bo, "bg": bg, "z": z,
        "rhoOil": rho_oil, "liveApi": live_api,
    })
    available = rs
final_rs = available
MASS_OUT = expected[-1]["liquidMass"] + cum_gas


def run(path):
    proc = subprocess.run([exe, path], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"solver failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")
    stages, total = [], None
    for line in proc.stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok or tok[0] != "separator":
            continue
        if len(tok) > 1 and tok[1] == "stage":
            d = {"index": int(tok[2])}
            rest = tok[3:]
            for a, b in zip(rest[0::2], rest[1::2]):
                try:
                    d[a] = float(b)
                except ValueError:
                    d[a] = b
            stages.append(d)
        elif len(tok) > 1 and tok[1] == "total":
            total = {a: float(b) for a, b in zip(tok[2::2], tok[3::2])}
    return stages, total


stages, total = run(case)
check("one stage line per case stage", len(stages) == len(pressures),
      f"got {len(stages)} want {len(pressures)}")
check("total line present", total is not None)

prev_rs = float("inf")
for i, (got, want) in enumerate(zip(stages, expected)):
    check(f"stage {i}: p and T echo the case",
          near(got.get("p", -1), want["p"], 1e-12)
          and near(got.get("T", -1), want["T"], 1e-12))
    for key in ("available", "rs", "gasSm3", "gasMass", "liquidMass",
                "cumGasMass", "bo", "bg", "z", "liveApi"):
        check(f"stage {i}: {key} matches independent restatement",
              key in got and near(got[key], want[key], 1e-9, 1e-12),
              f"got {got.get(key)} want {want[key]}")
    check(f"stage {i}: gas release non-negative", want["gasSm3"] >= -1e-12,
          f"gasSm3={want['gasSm3']}")
    check(f"stage {i}: rs non-increasing", want["rs"] <= prev_rs + 1e-12,
          f"rs={want['rs']} prev={prev_rs}")
    prev_rs = want["rs"]
    check(f"stage {i}: mass closes liquid + cumulative gas = feed",
          near(want["liquidMass"] + want["cumGasMass"], MASS_IN, 1e-9),
          f"{want['liquidMass'] + want['cumGasMass']} vs {MASS_IN}")

check("total gas = gor - final rs",
      near(total.get("gasSm3", -1), GOR - final_rs, 1e-9),
      f"got {total.get('gasSm3')} want {GOR - final_rs}")
check("massIn is blackOilMassPerOilSm3", near(total.get("massIn", -1), MASS_IN, 1e-9),
      f"got {total.get('massIn')} want {MASS_IN}")
check("massOut is final liquid + all liberated gas",
      near(total.get("massOut", -1), MASS_OUT, 1e-9),
      f"got {total.get('massOut')} want {MASS_OUT}")
check("reported imbalance is negligible",
      abs(total.get("imbalance", 1e9)) <= 1e-9 * MASS_IN,
      f"imbalance={total.get('imbalance')}")
check("stockTankApi from oilDensityStd",
      near(total.get("stockTankApi", -1), STOCK_API, 1e-9),
      f"got {total.get('stockTankApi')} want {STOCK_API}")

# ---- parse-time refusals, including the gated twin ----------------------
mutations = []
# gated twin: delete the whole separator{} block
m = re.search(r"\bseparator\s*\{", text)
if m:
    start = text.index("{", m.start())
    depth = 0
    end = None
    for j in range(start, len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                end = j + 1
                break
    if end:
        mutations.append(("gated twin: separator{} deleted",
                          text[:m.start()] + text[end:], "separator{...} block"))
mutations += [
    ("negative stage pressure",
     re.sub(r"pressure\s*\(\s*[-+0-9.eE]+", "pressure ( -1.0", text, count=1),
     "positive"),
    ("waterCut >= 1",
     re.sub(r"\bwaterCut\s+[-+0-9.eE]+\s*;", "waterCut 1.2;", text, count=1),
     "waterCut"),
]
with tempfile.TemporaryDirectory() as tmp:
    path = os.path.join(tmp, "probe.dict")
    for what, mutated, needle in mutations:
        if mutated == text:
            raise SystemExit(f"probe '{what}' did not mutate the case")
        open(path, "w").write(mutated)
        rr = subprocess.run([exe, path], capture_output=True, text=True)
        ok = rr.returncode != 0 and needle in rr.stdout + rr.stderr
        check(f"rejecting {what}", ok, f"needle '{needle}'")

if fails:
    raise SystemExit(f"pvtSeparator FAILED: {fails}")
print("pvtSeparator OK: independent Standing/VB/DAK cascade matches every"
      " stage, mass closes per stage, gated twin refused")
