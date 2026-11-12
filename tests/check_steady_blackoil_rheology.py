#!/usr/bin/env python3
"""5.76 / A5 increment 2: a declared liquid wall law on the black-oil march.

The published Beggs-Brill chain keeps its holdup, flow-pattern and
inclination correlations EXACTLY as they were (none of them reads a
viscosity). Only the liquid share of the no-slip mixture viscosity is
repriced on the mixture wall rate gdot = 8 vm / D:

    muM = muL_app * lambda + muG * (1 - lambda),
    muL_app = K' gdot^(n-1) + tau_0 (1 - exp(-100 gdot)) / gdot,
    K' = K ((3n+1)/(4n))^n                 (Rabinowitsch).

Zero declaration returns blackOilLiquid()'s muL bit for bit.

The shipped case is DEAD OIL (gor 0, waterCut 0), horizontal, both ends
fixed. Then Rs = 0, pb = 0, freeGas = 0, lambda = HL = 1, exp(S) = 1,
Ek = 0, elevation = 0 and the gas correlations drop out. In the laminar
branch the mixture density cancels and

    dp/ds = 32 muM vm / D^2 = (4/D) [K' gdot^n + tau_0 (1-e^-100gdot)]

is the exact single-phase Metzner-Reed wall law. The only pressure
dependence left is Bo(p) in the in-situ velocity vm = q Bo(p)/A; over the
5000 Pa drop that is a 1.8e-5 relative effect, so this checker RE-MARCHES
the 24 segments instead of treating vm as a constant.

The checker re-derives Bo(p), rhoOil(p), muL and the surface mass from the
dictionary's blackoil{} block alone (no solver internals), restates the
whole 24-segment march at each reported stock-tank rate, and asserts:

  * convergence, no choke, iterations = 0 (no free node),
  * the re-march lands on the reported downstream pressure,
  * the reported outlet velocity IS the re-march's last vm,
  * HL == 1, the reported regime is the published fr map at that vm, the
    reported mass rate is q times the surface mass per sm3,
  * Re < 2040 in every segment (64/Re exact, no turbulent member),
  * tau_w = dp D/(4L) is reproduced by the declared wall law at the
    outlet rate (to the Bo(p) drift bound),
  * flowIndex 1 with K = muL is bit-identical to the undeclared pipe,
  * K' is the Rabinowitsch consistency, not the textbook K,
  * the five laws are visibly distinct,
  * four parse-time refusals, loud and naming their reason.

usage: check_steady_blackoil_rheology.py <exe> <case.dict>
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


def near(got, want, rtol, atol=0.0):
    return abs(got - want) <= rtol * abs(want) + atol


def block(src, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", src)
    if not m:
        raise SystemExit(f"case has no {key}{{}} block")
    start = m.end() - 1
    depth = 0
    for j in range(start, len(src)):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                return src[start + 1:j]
    raise SystemExit(f"unterminated {key}{{}} block")


def entries(src):
    out = {}
    for m in re.finditer(r"(\w+)\s*\{", src):
        start = m.end() - 1
        depth = 0
        for j in range(start, len(src)):
            if src[j] == "{":
                depth += 1
            elif src[j] == "}":
                depth -= 1
                if depth == 0:
                    out[m.group(1)] = src[start + 1:j]
                    break
    return out


def num(src, key, default=0.0):
    m = re.search(r"\b" + re.escape(key) + r"\s+([0-9.eE+-]+)\s*;", src)
    return float(m.group(1)) if m else default


# ---- the black-oil correlations, restated from the case dictionary ----
PSI = 6894.757293168361
SCF_STB = 5.614583333333333
M_AIR = 28.9647e-3
R_GAS = 8.314462618
P_STD, T_STD = 101325.0, 288.15
G = 9.80665
RHO_WATER_STD = 999.014
LBFT3_TO_SI = 16.0184634

fl = block(text, "blackoil")
GG = num(fl, "gasGravity", 0.65)
RHO_STD = num(fl, "oilDensity", 850.0)
GOR = num(fl, "gor", 0.0)
WC = num(fl, "waterCut", 0.0)
T = num(fl, "temperature", 300.0)
WATER_DENSITY = num(fl, "waterDensity", 1000.0)
WATER_MU = num(fl, "waterViscosity", 1e-3)

gamma_o = RHO_STD / RHO_WATER_STD
api = 141.5 / gamma_o - 131.5
TF = T * 1.8 - 459.67
heavy = api > 30.0
c2 = 1.100e-5 if heavy else 1.751e-5
rs_scf = GOR * SCF_STB
t_term = (TF - 60.0) * api / GG
c1 = 4.670e-4 if heavy else 4.677e-4
c3 = 1.337e-9 if heavy else -1.811e-8
bob = 1.0 + c1 * rs_scf + c2 * t_term + c3 * rs_scf * t_term
# dead oil: Standing's bubble point is clamped to 0, and VB's co floors
# the pb denominator at 300 psia
pb_psia = max(18.2 * ((rs_scf / GG) ** 0.83 * 10.0
                      ** (0.00091 * TF - 0.0125 * api) - 1.4), 0.0)
check("dead oil has a clamped zero bubble point", pb_psia == 0.0,
      f"pb={pb_psia} psia")
co = (-1433.0 + 5.0 * rs_scf + 17.2 * TF - 1180.0 * GG + 12.61 * api) \
     / (1.0e5 * max(pb_psia, 300.0))
mu_dead = 10.0 ** (10.0 ** (3.0324 - 0.02023 * api) * TF ** (-1.163)) - 1.0
mu_oil = 10.715 * (rs_scf + 100.0) ** (-0.515) \
         * mu_dead ** (5.44 * (rs_scf + 150.0) ** (-0.338)) * 1e-3
wor = WC / (1.0 - WC) if 0.0 < WC < 0.999 else 0.0


def bo_of(p):
    # above bubble point: Bo = Bob exp[co (pb - p)]
    return bob * math.exp(co * (pb_psia - p / PSI))


def rho_oil_of(p):
    return (350.0 * gamma_o + 0.0764 * GG * rs_scf) \
           / (SCF_STB * bo_of(p)) * LBFT3_TO_SI


def mu_liquid_of(p):
    v_o = bo_of(p)
    v_t = v_o + wor
    fo = v_o / v_t
    return math.exp(fo * math.log(mu_oil) + (1.0 - fo) * math.log(WATER_MU))


rho_gas_std = P_STD * GG * M_AIR / (R_GAS * T_STD)
MASS_PER_SM3 = RHO_STD + GOR * rho_gas_std + wor * WATER_DENSITY

REG = 100.0
SHEAR_FLOOR = 1e-6


def apparent_mu(pipe, gdot, mu_newtonian):
    """wallRateApparentViscosityPaS, identical to cf/common/Rheology.h."""
    n, K, t0 = pipe["n"], pipe["K"], pipe["tau0"]
    if n == 0.0 and t0 == 0.0:
        return mu_newtonian
    g = max(gdot, SHEAR_FLOOR)
    if n > 0.0:
        kp = K * ((3.0 * n + 1.0) / (4.0 * n)) ** n
        mu = kp * g ** (n - 1.0)
    else:
        mu = mu_newtonian
    if t0 > 0.0:
        mu += t0 * (-math.expm1(-REG * g)) / g
    return mu


def wall_stress(pipe, gdot, mu_newtonian):
    return apparent_mu(pipe, gdot, mu_newtonian) * max(gdot, SHEAR_FLOOR)


AB = (0.98, 0.0, 0.845, 1.065)
BB = (0.4846, 0.0, 0.5351, 0.5824)
CB = (0.0868, 0.0, 0.0173, 0.0609)
REGIME_NAME = ("segregated", "transition", "intermittent", "distributed")


def bb_gradient(pipe, vm, p):
    """beggsBrillDpDs for the dead-oil, horizontal, lambda = 1 case."""
    D = pipe["D"]
    lam = 1.0
    fr = vm * vm / (G * D)
    l1 = 316.0 * lam ** 0.302
    l2 = 0.0009252 * lam ** (-2.4684)
    l3 = 0.1 * lam ** (-1.4516)
    l4 = 0.5 * lam ** (-6.738)
    if (lam < 0.01 and fr < l1) or (lam >= 0.01 and fr < l2):
        regime = 0
    elif lam >= 0.01 and l2 <= fr <= l3:
        regime = 1
    elif (0.01 <= lam < 0.4 and l3 < fr <= l1) \
            or (lam >= 0.4 and l3 < fr <= l4):
        regime = 2
    else:
        regime = 3

    def hl0(r):
        return max(AB[r] * lam ** BB[r] * fr ** (-CB[r]), lam)

    if regime == 1:
        w = (l3 - fr) / (l3 - l2)
        hl = w * hl0(0) + (1.0 - w) * hl0(1)
    else:
        hl = hl0(regime)
    hl = min(max(hl, 0.0), 1.0)

    rho_l = rho_oil_of(p)
    rho_m = rho_l * lam           # + rhoG * (1 - lambda) == 0
    rho_s = rho_l * hl            # + rhoG * (1 - hl) == 0
    mu_l = mu_liquid_of(p)
    mu_app = apparent_mu(pipe, 8.0 * vm / D, mu_l)
    re = rho_m * vm * D / mu_app
    check_re = re < 2040.0
    f_ns = 64.0 / max(re, 1e-12)
    y = lam / (hl * hl)
    ly = math.log(y)
    if 1.0 < y < 1.2:
        s = math.log(2.2 * y - 1.2)
    else:
        s = ly / (ly * (ly * (0.01853 * ly * ly - 0.8725) + 3.182) - 0.0523)
    s = min(s, 7.0)
    f_tp = f_ns * math.exp(s)
    fric = f_tp * rho_m * vm * vm / (2.0 * D)
    ek = 0.0                      # vsg == 0
    return fric / max(1.0 - ek, 0.1), regime, hl, re, f_tp, s


NODE = {k: num(v, "pressure") for k, v in entries(block(text, "nodes")).items()}
PIPES = {}
for name, body in entries(block(text, "pipes")).items():
    PIPES[name] = {
        "from": re.search(r"\bfrom\s+(\w+)\s*;", body).group(1),
        "to": re.search(r"\bto\s+(\w+)\s*;", body).group(1),
        "L": num(body, "length"), "D": num(body, "diameter"),
        "eps": num(body, "roughness"), "z": num(body, "elevation"),
        "n": num(body, "flowIndex"), "K": num(body, "consistency"),
        "tau0": num(body, "yieldStress"),
    }


def run(path):
    proc = subprocess.run([exe, path], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"solver failed rc={proc.returncode}\n"
                         f"{proc.stdout}\n{proc.stderr}")
    nodes, pipes, raw, converged, iters = {}, {}, {}, False, None
    for line in proc.stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "node":
            nodes[tok[1]] = float(tok[3])
        elif tok[0] == "pipe":
            d, rest = {}, tok[2:]
            for a, b in zip(rest[0::2], rest[1::2]):
                d[a] = b
            for key in ("Q", "v", "HL", "m", "dp"):
                if key in d:
                    d[key] = float(d[key])
            pipes[tok[1]] = d
            raw[tok[1]] = line.strip()
        elif tok[0] == "steady":
            converged = "(converged)" in proc.stdout
            m = re.search(r"iterations=(\d+)", proc.stdout)
            iters = int(m.group(1)) if m else None
    return nodes, pipes, raw, converged, iters


nodes, pipes, raw, converged, iters = run(case)
check("steady solve converged with no free node (iterations = 0)",
      converged and iters == 0, f"iterations={iters}")
check("no pipe choked", all("CHOKED" not in line for line in raw.values()))
check("the two fixed pressures are the case's",
      near(nodes["A"], NODE["A"], 1e-14)
      and near(nodes["B"], NODE["B"], 1e-14))

for name, p in PIPES.items():
    if name not in pipes:
        raise SystemExit(f"missing pipe {name} in solver output")
    r = pipes[name]
    D, L, q = p["D"], p["L"], r["Q"]
    area = math.pi * D * D / 4.0
    p_up = max(NODE[p["from"]], NODE[p["to"]])
    p_down = min(NODE[p["from"]], NODE[p["to"]])
    ds = L / 24.0

    pcur, vm_last, hl_last, re_max, regime_max, s_max = p_up, 0.0, 0.0, 0.0, 0, 0.0
    all_laminar = True
    for _ in range(24):
        vm = q * bo_of(pcur) / area
        grad, regime, hl, re, f_tp, s = bb_gradient(p, vm, pcur)
        all_laminar = all_laminar and re < 2040.0
        pcur -= grad * ds
        vm_last, hl_last = vm, hl
        re_max, regime_max, s_max = max(re_max, re), max(regime_max, regime), max(s_max, abs(s))
    check(f"{name}: every segment is laminar (Re < 2040)", all_laminar,
          f"Re_max={re_max:.4g}")
    check(f"{name}: re-march lands on the downstream node",
          abs(pcur - p_down) <= 2e-7 * p_up,
          f"err={pcur - p_down:.4g} Pa (solver secant tol = {1e-7 * p_up:.4g} Pa)")
    check(f"{name}: reported outlet velocity is the re-march's last vm",
          near(r["v"], vm_last, 1e-9),
          f"got {r['v']:.10g} want {vm_last:.10g}")
    check(f"{name}: holdup is the published HL = 1",
          near(r["HL"], 1.0, 1e-12) and near(hl_last, 1.0, 1e-12))
    check(f"{name}: two-phase multiplier exp(S) is exactly 1",
          near(s_max, 0.0, 1e-12))
    check(f"{name}: regime is the published fr map (no viscosity input)",
          r["regime"] == REGIME_NAME[regime_max],
          f"reported {r['regime']}, restated {REGIME_NAME[regime_max]}")
    check(f"{name}: mass rate = q * surface mass per sm3",
          near(r["m"], q * MASS_PER_SM3, 1e-9),
          f"{r['m']:.9g} vs {q * MASS_PER_SM3:.9g}")
    check(f"{name}: dp is the imposed drop",
          near(r["dp"], p_up - p_down, 1e-12))

    # the wall stress is geometric; restate it at the OUTLET rate. Bo(p)
    # drifts < 2e-5 over the march, so this is a 1e-4 statement, not 1e-12.
    tau_w = (p_up - p_down) * D / (4.0 * L)
    tau_out = wall_stress(p, 8.0 * vm_last / D, mu_liquid_of(pcur))
    check(f"{name}: tau_w = dp D/(4L) is the declared wall law",
          near(tau_out, tau_w, 1e-4),
          f"tau_out {tau_out:.9g} vs geometric {tau_w:.9g}")

# ---- declared-but-Newtonian is bit-identical on the printed line -------
def payload(line):
    return " ".join(line.split()[3:])


check("flowIndex 1 with K = muL is bit-identical to the undeclared pipe",
      payload(raw["Newt"]) == payload(raw["PowN1"]),
      f"Newt [{payload(raw['Newt'])}] PowN1 [{payload(raw['PowN1'])}]")
check("a viscosity declaration at fixed rate cannot move HL or regime",
      pipes["Newt"]["HL"] == pipes["PowN1"]["HL"]
      and pipes["Newt"]["regime"] == pipes["PowN1"]["regime"])

# ---- K' is the Rabinowitsch consistency, not the textbook K -----------
thin = PIPES["Thin05"]
n, K = thin["n"], thin["K"]
kp = K * ((3.0 * n + 1.0) / (4.0 * n)) ** n
ratio = (K / kp) ** (1.0 / n)
check("Rabinowitsch K' = K ((3n+1)/(4n))^n is what the wall carries",
      near(ratio, 0.8, 1e-12), f"K'/K = {kp / K:.9g}")
newt_v, thin_v = pipes["Newt"]["v"], pipes["Thin05"]["v"]
tau_geom = (NODE["A"] - NODE["B"]) * thin["D"] / (4.0 * thin["L"])
v_kprime = (tau_geom / kp) ** (1.0 / n) * thin["D"] / 8.0
v_textbook = (tau_geom / K) ** (1.0 / n) * thin["D"] / 8.0
check("the undeclared pipe is the Newtonian wall law v = tau_w D/(8 muL)",
      near(newt_v, (tau_geom / mu_oil) * thin["D"] / 8.0, 1e-4),
      f"got {newt_v:.9g}")
check("the measured rate follows K', not the 1.25x textbook-K rate",
      near(thin_v, v_kprime, 1e-4)
      and abs(thin_v - v_textbook) > 0.1 * thin_v,
      f"v {thin_v:.9g}, K' {v_kprime:.9g}, textbook K {v_textbook:.9g}")

# ---- the five laws are distinct ---------------------------------------
qs = [pipes[k]["Q"] for k in ("Newt", "Thin05", "Yield", "HB")]
check("Newtonian / power-law / Bingham / Herschel-Bulkley are distinct",
      all(qs[i] != qs[j] for i in range(4) for j in range(i + 1, 4)),
      f"Q = {[f'{x:.6g}' for x in qs]}")
check("ordering follows the wall laws at tau_w = 0.125 Pa",
      qs[3] < qs[1] < qs[2] < qs[0])

# ---- parse-time refusals: loud, not silent defaults -------------------
KPOW = "0.044721359549995794"
probes = [
    ("flowIndex without a consistency",
     text.replace("flowIndex 0.5; consistency " + KPOW + ";",
                  "flowIndex 0.5;", 1), "positive consistency"),
    ("consistency without a flowIndex",
     text.replace("flowIndex 0.5; consistency " + KPOW + ";",
                  "consistency " + KPOW + ";", 1),
     "not a quantity the wall law reads"),
    ("negative flowIndex",
     text.replace("flowIndex 0.5;", "flowIndex -0.5;", 1), "negative"),
    ("negative yieldStress",
     text.replace("yieldStress 0.05; }", "yieldStress -0.05; }", 1),
     "negative yield"),
]
with tempfile.TemporaryDirectory() as tmp:
    path = os.path.join(tmp, "probe.dict")
    for what, mutated, needle in probes:
        if mutated == text:
            raise SystemExit(f"probe '{what}' did not mutate the case")
        open(path, "w").write(mutated)
        rr = subprocess.run([exe, path], capture_output=True, text=True)
        ok = rr.returncode != 0 and needle in rr.stdout + rr.stderr
        check(f"rejecting {what}", ok, f"needle '{needle}'")

if fails:
    raise SystemExit(f"steadyBlackOilRheology FAILED: {fails}")
print("steadyBlackOilRheology OK: BB liquid share repriced on the wall rate,"
      " published holdup/regime untouched; 24-segment re-march + K' + "
      "bit-identity + 4 refusals")
