#!/usr/bin/env python3
"""5.144 / G6: the PASSIVE steady temperature march of a marched pipe.

The claim under test is in the word PASSIVE: the Joule-Thomson and wall
terms price the outlet temperature along the pressure profile the
hydraulics already produced, so a thermal declaration cannot move a
single node pressure, rate, holdup or dp number - it can only append a
` T <K>` token. This checker re-derives the whole march from the case
dictionary alone - the black-oil PVT and Beggs-Brill gradient, the pvt
table march with its Churchill law, and the exact per-segment convective
relaxation - and asserts:

  * the standing multiphase case carries no T token anywhere (undeclared
    is the isothermal march),
  * a Joule-Thomson declaration on the standing case leaves every line
    byte-identical once the appended T token is stripped, and the token
    itself is T_in + jt * (p_out - p_in) on the solved node pressures -
    the per-segment steps telescope exactly,
  * a constant-wall declaration on a black-oil anchor reproduces the
    segment-by-segment relaxation T <- Tw + (T - Tw) exp(-U ds / (m c_p))
    on the re-derived march profile, landing on the terminal pressure,
  * a geothermal declaration on a pvt-table gas anchor reproduces the
    same law on the table march's own profile (constant table: the
    interpolation is exact by construction),
  * the wall profile is anchored at the inlet datum: a horizontal pipe
    sees surfaceTemperature everywhere;
  * refusals: thermal keys on a constant-property network, a
    wallHeatTransferSegments list of the wrong size, wall exchange
    without heatCapacity, thermal physics on a pvt fluid without a
    temperature, wall keys on a pure Joule-Thomson march, an unknown
    wallModel word, a negative wallHeatTransfer, and the keys on the
    transient row and the nodal VLP well.

usage: check_steady_thermal.py <steady-exe> <case.dict> \
        [<transient-exe> <nodal-exe>]
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) < 3:
    raise SystemExit(__doc__)

STEADY, CASE = sys.argv[1], os.path.abspath(sys.argv[2])
TRANSIENT = sys.argv[3] if len(sys.argv) > 3 else None
NODAL = sys.argv[4] if len(sys.argv) > 4 else None
HERE = os.path.dirname(os.path.abspath(__file__))
text = re.sub(r"//[^\n]*", "", open(CASE).read())

fails = []


def check(name, cond, detail=""):
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def near(got, want, rtol):
    return rel(got, want) <= rtol


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


def word(src, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(\w+)\s*;", src)
    return m.group(1) if m else None


# ---- constants ----
PSI = 6894.757293168361
SCF = 5.614583333333333
M_AIR = 28.9647e-3
R_GAS = 8.314462618
P_STD, T_STD = 101325.0, 288.15
G = 9.80665
RHO_WATER_STD = 999.014
LBFT3 = 16.0184634

fl = block(text, "blackoil")
GG = num(fl, "gasGravity", 0.65)
RHO_STD = num(fl, "oilDensity", 850.0)
GOR = num(fl, "gor", 0.0)
WC = num(fl, "waterCut", 0.0)
T_FLUID = num(fl, "temperature", 300.0)
SIG = num(fl, "surfaceTension", 0.025)
WATER_DENSITY = num(fl, "waterDensity", 1000.0)
WATER_MU = num(fl, "waterViscosity", 1e-3)

AB = (0.98, 0.0, 0.845, 1.065)
BB = (0.4846, 0.0, 0.5351, 0.5824)
CB = (0.0868, 0.0, 0.0173, 0.0609)
REGIME_NAME = ("segregated", "transition", "intermittent", "distributed")


def boundaries(lam):
    return (316.0 * lam ** 0.302, 0.0009252 * lam ** -2.4684,
            0.1 * lam ** -1.4516, 0.5 * lam ** -6.738)


def regime_of(lam, nfr):
    l1, l2, l3, l4 = boundaries(lam)
    if (lam < 0.01 and nfr < l1) or (lam >= 0.01 and nfr < l2):
        return 0
    if lam >= 0.01 and l2 <= nfr <= l3:
        return 1
    if (0.01 <= lam < 0.4 and l3 < nfr <= l1) or (lam >= 0.4 and l3 < nfr <= l4):
        return 2
    return 3


def hl0_of(lam, nfr, r):
    return max(AB[r] * lam ** BB[r] * nfr ** -CB[r], lam)


def hl_horizontal(lam, nfr, r):
    if r == 1:
        l1, l2, l3, l4 = boundaries(lam)
        w = (l3 - nfr) / (l3 - l2)
        return w * hl0_of(lam, nfr, 0) + (1.0 - w) * lam
    return hl0_of(lam, nfr, r)


UPHILL_R0 = (0.011, -3.768, 3.539, -1.614)
UPHILL_R2 = (2.96, 0.305, -0.4473, 0.0978)
DOWNHILL = (4.70, -0.3692, 0.1244, -0.5056)
NONE = "no-correction"


def incline_coeff(angle_deg, r):
    if abs(angle_deg) <= 1e-12:
        return NONE
    if angle_deg > 0.0:
        if r == 3:
            return NONE
        return UPHILL_R0 if r == 0 else UPHILL_R2
    return DOWNHILL


def hl_priced(lam, nfr, r, vsl, rho_l, angle_deg):
    hl = hl_horizontal(lam, nfr, r)
    coeff = incline_coeff(angle_deg, r)
    if coeff != NONE and lam > 1e-9:
        d, e, f, h = coeff
        nlv = vsl * math.sqrt(math.sqrt(rho_l / (G * SIG)))
        c = max((1.0 - lam) * math.log(d * lam ** e * nlv ** f * nfr ** h), 0.0)
        s1 = math.sin(1.8 * math.radians(angle_deg))
        hl = hl * (1.0 + c * (s1 - s1 ** 3 / 3.0))
    return min(max(hl, 0.0), 1.0)


def colebrook(re, ed):
    if re < 2040.0:
        return 64.0 / max(re, 1e-12)
    ed = max(ed, 0.0)
    y = max(1.14 - 2.0 * math.log10(ed + 21.25 / re ** 0.9), 1.0)
    for _ in range(100):
        yn = -2.0 * math.log10(ed / 3.7 + 2.51 * y / re)
        err = abs(yn - y) / yn
        y = yn
        if err < 1e-14:
            break
    return 1.0 / (y * y)


def churchill(re, ed):
    re = max(re, 1e-12)
    a1 = (7.0 / re) ** 0.9 + 0.27 * ed
    a = (2.457 * math.log(1.0 / a1)) ** 16.0
    b = (37530.0 / re) ** 16.0
    return 8.0 * ((8.0 / re) ** 12.0 + (a + b) ** -1.5) ** (1.0 / 12.0)


class BlackOil:
    """evaluateBlackOil + marchBeggsBrillPipe, restated; the march returns
    the per-segment pressure profile so the thermal relaxation can be
    re-derived on exactly the dp the kernel saw."""

    def __init__(self, gg, rho_std, gor, wc, temp, water_density, water_mu):
        self.gg, self.rho_std, self.gor, self.wc = gg, rho_std, gor, wc
        self.T, self.wd, self.wmu = temp, water_density, water_mu
        self.api = 141.5 / (rho_std / RHO_WATER_STD) - 131.5
        self.tf = temp * 1.8 - 459.67
        self.tr = temp * 1.8
        self.gor_scf = gor * SCF
        self.wor = wc / (1.0 - wc) if 0.0 < wc < 0.999 else 0.0
        rho_g_std = P_STD * gg * M_AIR / (R_GAS * T_STD)
        self.mass_per_sm3 = rho_std + gor * rho_g_std + self.wor * water_density

    def pb_psia(self):
        a = 0.00091 * self.tf - 0.0125 * self.api
        return max(18.2 * ((self.gor_scf / self.gg) ** 0.83 * 10 ** a - 1.4),
                   0.0)

    def rs_scf(self, p_psia):
        c1, c2, c3 = ((0.0362, 1.0937, 25.724) if self.api <= 30
                      else (0.0178, 1.187, 23.931))
        return (c1 * self.gg * max(p_psia, 0.0) ** c2
                * math.exp(c3 * self.api / (self.tf + 459.67)))

    def bo_of(self, p_psia, rs, pb):
        api, gg, tf = self.api, self.gg, self.tf
        c1, c2, c3 = ((4.670e-4, 1.100e-5, 1.337e-9) if api > 30
                      else (4.677e-4, 1.751e-5, -1.811e-8))
        tt = (tf - 60.0) * api / gg
        bob = 1 + c1 * rs + c2 * tt + c3 * rs * tt
        if p_psia <= pb:
            return bob
        co = ((-1433 + 5 * rs + 17.2 * tf - 1180 * gg + 12.61 * api)
              / (1e5 * p_psia))
        return bob * math.exp(co * (pb - p_psia))

    def z_dak(self, p_psia):
        gg, tr = self.gg, self.tr
        if gg < 0.75:
            tpc = 168 + 325 * gg - 12.5 * gg * gg
            ppc = 677 + 15 * gg - 37.5 * gg * gg
        else:
            tpc = 187 + 330 * gg - 71.5 * gg * gg
            ppc = 706 + 51 * gg - 11.1 * gg * gg
        A = [0.3265, -1.0700, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361,
             0.1844, 0.1056, 0.6134, 0.7210]
        pr, trr = max(p_psia, 1e-6) / ppc, tr / tpc
        t1 = A[0] + A[1] / trr + A[2] / trr ** 3 + A[3] / trr ** 4 + A[4] / trr ** 5
        t2 = A[5] + A[6] / trr + A[7] / trr ** 2

        def zofrho(rho):
            r2 = rho * rho
            return (1 + t1 * rho + t2 * r2
                    - A[8] * (A[6] / trr + A[7] / trr ** 2) * rho ** 5
                    + A[9] * (1 + A[10] * r2) * r2 / trr ** 3
                    * math.exp(-A[10] * r2))

        def res(z):
            return zofrho(0.27 * pr / (z * trr)) - z

        lo, hi = 0.05, 3.0
        if res(lo) <= 0:
            return lo
        if res(hi) >= 0:
            return hi
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if res(mid) > 0:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-13 * hi:
                break
        return 0.5 * (lo + hi)

    def state(self, p_pa):
        ppsia = max(p_pa, 1.0) / PSI
        pb = self.pb_psia()
        if p_pa <= pb * PSI:
            rs = min(self.rs_scf(ppsia) / SCF, self.gor)
            b = self.bo_of(ppsia, rs * SCF, pb)
            fg = self.gor - rs
        else:
            rs, fg = self.gor, 0.0
            b = self.bo_of(ppsia, self.gor_scf, pb)
        z = self.z_dak(ppsia)
        bg = (P_STD / p_pa) * (self.T / T_STD) * z
        rhog = p_pa * self.gg * M_AIR / (z * R_GAS * self.T)
        rho_oil = ((350 * (self.rho_std / RHO_WATER_STD)
                    + 0.0764 * self.gg * rs * SCF) / (SCF * b)) * LBFT3
        x = 10 ** (3.0324 - 0.02023 * self.api) * self.tf ** -1.163
        mu_dead = 10 ** x - 1
        mu_oil = (10.715 * (rs * SCF + 100) ** -0.515) * mu_dead ** (
            5.44 * (rs * SCF + 150) ** -0.338) * 1e-3
        mw = self.gg * 28.9647
        k = (9.4 + 0.02 * mw) * self.tr ** 1.5 / (209 + 19 * mw + self.tr)
        xx = 3.5 + 986 / self.tr + 0.01 * mw
        yy = 2.4 - 0.2 * xx
        mu_gas = 1e-4 * k * math.exp(xx * (rhog / 1000) ** yy) * 1e-3
        v_t = b + self.wor
        fo = b / v_t
        rho_l = (rho_oil * b + self.wd * self.wor) / v_t
        mu_l = math.exp(fo * math.log(mu_oil) + (1.0 - fo) * math.log(self.wmu))
        return dict(bo=b, fg=fg, bg=bg, rhog=rhog, rhol=rho_l, mul=mu_l,
                    mug=mu_gas, wor=self.wor)

    def march_profile(self, q, pin, L, D, rough, sintheta=0.0, nseg=24):
        area = 0.25 * math.pi * D * D
        ds = L / nseg
        p_floor = max(0.02 * pin, 1000.0)
        angle_deg = math.degrees(math.asin(max(-1.0, min(1.0, sintheta))))
        p = pin
        steps = []
        ok = True
        for _ in range(nseg):
            st = self.state(p)
            vsl = q * (st["bo"] + st["wor"]) / area
            vsg = q * st["fg"] * st["bg"] / area
            vm = vsl + vsg
            lam = vsl / vm
            nfr = vm * vm / (G * D)
            r = regime_of(lam, nfr)
            hl = hl_priced(lam, nfr, r, vsl, st["rhol"], angle_deg)
            rhom = st["rhol"] * lam + st["rhog"] * (1 - lam)
            rhos = st["rhol"] * hl + st["rhog"] * (1 - hl)
            mum = st["mul"] * lam + st["mug"] * (1 - lam)
            re = rhom * vm * D / mum
            fns = colebrook(re, rough / D)
            if lam < 1e-9 or hl < 1e-9:
                sm = 0.0
            else:
                yy = lam / (hl * hl)
                ly = math.log(yy)
                if 1.0 < yy < 1.2:
                    sm = math.log(2.2 * yy - 1.2)
                else:
                    sm = (ly / (ly * (ly * (0.01853 * ly * ly - 0.8725)
                                      + 3.182) - 0.0523))
                sm = min(sm, 7.0)
            ftp = fns * math.exp(sm)
            ek = vsg * vm * rhos / max(p, 1000.0)
            ele = G * math.sin(math.radians(angle_deg)) * rhos
            dpds = (ele + ftp * rhom * vm * vm / (2 * D)) / max(1 - ek, 0.1)
            pOld = p
            p -= dpds * ds
            if not p > p_floor:
                p = p_floor
                ok = False
            steps.append((pOld, p))
            if not ok:
                break
        return steps, p, ok


def thermal_steps(steps, t0, jt, u_of_seg, tw_of_seg, ds, mass_rate, cp):
    """The kernel's steadyThermalAdvance, segment by segment."""
    t = t0
    for i, (p_old, p_new) in enumerate(steps):
        t += jt * (p_new - p_old)
        u = u_of_seg(i)
        if u > 0.0 and mass_rate > 0.0:
            tw = tw_of_seg(i)
            t = tw + (t - tw) * math.exp(-u * ds / (mass_rate * cp))
    return t


# ---- run the solver and parse its readout ----
def parse_kv(toks):
    d, i = {}, 0
    while i < len(toks):
        t = toks[i]
        try:
            float(t)
            i += 1
            continue
        except ValueError:
            pass
        if t == "mapRegime":
            d["mapRegime"] = toks[i + 1]
            i += 3
            continue
        if t in REGIME_NAME or t == "FORCED":
            i += 1
            continue
        d[t] = toks[i + 1]
        i += 2
    return d


def run(case_path, exe=None):
    proc = subprocess.run([exe or STEADY, str(case_path)],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"solver failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")
    if "steady solved" not in proc.stdout or "(converged)" not in proc.stdout:
        raise SystemExit("solver did not report convergence")
    nodes, pipes, lines = {}, {}, {}
    for line in proc.stdout.splitlines():
        toks = line.replace("[chflow]", "").split()
        if not toks:
            continue
        if toks[0] == "node":
            nodes[toks[1]] = float(toks[3])
        elif toks[0] == "pipe":
            if "CHOKED" in toks:
                raise SystemExit(f"pipe {toks[1]} choked: {line}")
            pipes[toks[1]] = parse_kv(toks[2:])
            lines[toks[1]] = " ".join(toks)
    return nodes, pipes, lines


def run_stdout(body, exe=None):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "t.dict")
        with open(path, "w") as fh:
            fh.write(body)
        proc = subprocess.run([exe or STEADY, path], capture_output=True,
                              text=True, cwd=tmp)
    if proc.returncode != 0:
        raise SystemExit(f"solver failed\n{proc.stdout}\n{proc.stderr}")
    return proc.stdout


def run_dict(body):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "t.dict")
        with open(path, "w") as fh:
            fh.write(body)
        return run(path)


def strip_T(line):
    return re.sub(r" T [0-9.eE+-]+$", "", line)


FLUID = BlackOil(GG, RHO_STD, GOR, WC, T_FLUID, WATER_DENSITY, WATER_MU)
src = open(CASE).read()

# ---- 1. undeclared is the isothermal march: no T token anywhere ----
base_out = run_stdout(src)
check("the standing case carries no T token", " T " not in base_out)

# ---- 2. Joule-Thomson passivity on the standing case ----
pipe_entries = entries(block(text, "pipes"))
p3_name = p3_body = None
for name, body in pipe_entries.items():
    if name == "P3":
        p3_name, p3_body = name, body
if p3_body is None:
    raise SystemExit("probe: the standing case has no P3 pipe")
JT = 4.0e-7
patched = src.replace(p3_body, p3_body + f"jouleThomson {JT:.6g}; ")
jt_out = run_stdout(patched)
base_lines = base_out.splitlines()
jt_lines = jt_out.splitlines()
check("the JT declaration moves no other byte",
      len(base_lines) == len(jt_lines)
      and all(strip_T(a) == b for a, b in zip(jt_lines, base_lines)))
jt_nodes, jt_pipes, _ = run_dict(patched)
p_in, p_out = jt_nodes[word(p3_body, "from")], jt_nodes[word(p3_body, "to")]
want = T_FLUID + JT * (p_out - p_in)
got = float(jt_pipes[p3_name]["T"])
check("the JT token is T_in + jt (p_out - p_in) on the solved pressures",
      near(got, want, 1e-9), f"{got:.6f} vs {want:.6f}")

# ---- 3. constant-wall anchor on the black-oil march ----
OUT_P = jt_nodes["Outlet"]


def bo_dict(decl, rate=0.004):
    return f"""
fluid {{ blackoil {{ gasGravity {GG}; oilDensity {RHO_STD}; gor {GOR:g};
    waterCut {WC}; temperature {T_FLUID}; surfaceTension {SIG}; }} }}
nodes {{ W {{}}; Outlet {{ pressure {OUT_P:.1f}; }} }}
pipes {{ P {{ from W; to Outlet; length 800; diameter 0.1; roughness 4.6e-5;
    elevation 0;{decl} }} }}
sources {{ s {{ node W; rate {rate}; }} }}
solve {{ maxIter 100; tol 1e-10; }}
"""


U_W, TW_W, CP_W = 2.0, 380.0, 2000.0
wall_decl = (f" wallHeatTransfer {U_W}; wallTemperature {TW_W};"
             f" heatCapacity {CP_W};")
wn, wp, wl = run_dict(bo_dict(wall_decl))
bn, bp, bl = run_dict(bo_dict(""))
check("the wall declaration moves no hydraulic byte",
      strip_T(wl["P"]) == bl["P"],
      f"{strip_T(wl['P'])[-120:]}\n vs {bl['P'][-120:]}")
steps, p_end, ok = FLUID.march_profile(0.004, wn["W"], 800.0, 0.1, 4.6e-5)
check("wall anchor: the re-derived march lands on terminal",
      ok and near(p_end, wn["Outlet"], 1e-6), f"{p_end:.3f} vs {wn['Outlet']:.3f}")
mass_rate = 0.004 * FLUID.mass_per_sm3
t_wall = thermal_steps(steps, T_FLUID, 0.0, lambda i: U_W,
                       lambda i: TW_W, 800.0 / 24, mass_rate, CP_W)
check("wall anchor: the T token is the segment-by-segment relaxation",
      near(float(wp["P"]["T"]), t_wall, 1e-8),
      f"{float(wp['P']['T']):.6f} vs {t_wall:.6f}")
check("wall anchor: a hotter wall heats the outlet",
      float(wp["P"]["T"]) > T_FLUID, f"{float(wp['P']['T']):.3f} > {T_FLUID}")

# ---- 4. geothermal pvt-table anchor ----
PVT_RHO, PVT_MU, PVT_T = 55.0, 1.2e-5, 320.0
U_G, TS_G, GRAD_G, CP_G = 0.5, 300.0, 0.03, 2200.0


def pvt_dict(decl):
    return f"""
fluid {{ pvt {{ pressure ( 1.0e5 1.0e7 ); rho ( {PVT_RHO} {PVT_RHO} );
    mu ( {PVT_MU} {PVT_MU} ); temperature {PVT_T}; }} }}
nodes {{ In {{ pressure 1.0e7; }} Outlet {{ pressure 5.0e6; }} }}
pipes {{ P {{ from In; to Outlet; length 1000; diameter 0.1;
    roughness 4.6e-5; elevation 0;{decl} }} }}
solve {{ maxIter 100; tol 1e-10; }}
"""


geo_decl = (f" wallHeatTransfer {U_G}; wallModel geothermal;"
            f" surfaceTemperature {TS_G}; geothermalGradient {GRAD_G};"
            f" heatCapacity {CP_G};")
gn, gp, gl = run_dict(pvt_dict(geo_decl))
g0n, g0p, g0l = run_dict(pvt_dict(""))
check("the geothermal declaration moves no hydraulic byte",
      strip_T(gl["P"]) == g0l["P"])
Q_PVT = float(gp["P"]["Q"])


def pvt_profile(q, pin, L, D, rough):
    area = 0.25 * math.pi * D * D
    ds = L / 24
    mass = PVT_RHO * q
    p = pin
    steps = []
    ok = True
    for _ in range(24):
        v = mass / (PVT_RHO * area)
        re = PVT_RHO * v * D / PVT_MU
        f = churchill(re, rough / D)
        dp = (f * ds / D) * PVT_RHO * v * v / 2.0
        pOld = p
        p -= dp
        if not p > max(0.02 * pin, 1000.0):
            p = max(0.02 * pin, 1000.0)
            ok = False
        steps.append((pOld, p))
        if not ok:
            break
    return steps, p, ok


gsteps, gp_end, gok = pvt_profile(Q_PVT, gn["In"], 1000.0, 0.1, 4.6e-5)
check("pvt anchor: the re-derived table march lands on terminal",
      gok and near(gp_end, gn["Outlet"], 1e-6),
      f"{gp_end:.3f} vs {gn['Outlet']:.3f}")
# a horizontal pipe: every segment centre sits at the inlet datum, so the
# geothermal profile prices surfaceTemperature everywhere
tg = lambda i: TS_G + GRAD_G * 0.0
t_geo = thermal_steps(gsteps, PVT_T, 0.0, lambda i: U_G, tg,
                      1000.0 / 24, PVT_RHO * Q_PVT, CP_G)
check("pvt anchor: the T token is the relaxation on the table profile",
      near(float(gp["P"]["T"]), t_geo, 1e-8),
      f"{float(gp['P']['T']):.6f} vs {t_geo:.6f}")
check("pvt anchor: horizontal geothermal anchors at surfaceTemperature",
      PVT_T > TS_G and float(gp["P"]["T"]) < PVT_T
      and float(gp["P"]["T"]) > TS_G)

# ---- 5. refusals ----
def refuses(label, body, want, exe=None):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "r.dict")
        with open(path, "w") as fh:
            fh.write(body)
        rr = subprocess.run([exe or STEADY, path], capture_output=True,
                            text=True, cwd=tmp)
    out = (rr.stdout + rr.stderr).strip()
    check(f"refused: {label}", rr.returncode != 0 and want in out,
          out.splitlines()[-1][:130] if out else "no message")


refuses("thermal keys on a constant-property network", """
fluid { rho 850; mu 0.1; }
nodes { W {}; Outlet { pressure 4.0e6; } }
pipes { P { from W; to Outlet; length 500; diameter 0.08; roughness 4.6e-5;
    elevation 0; jouleThomson 1e-7; } }
sources { s { node W; rate 0.004; } }
solve { maxIter 100; tol 1e-10; }
""", "without a marched fluid")
refuses("a wallHeatTransferSegments list of the wrong size",
        bo_dict(" wallHeatTransferSegments (1 2 3); heatCapacity 2000;"
                " wallTemperature 350;"),
        "exactly 24")
refuses("wall exchange without heatCapacity",
        bo_dict(" wallHeatTransfer 2; wallTemperature 350;"),
        "positive heatCapacity")
refuses("thermal on a pvt fluid with no temperature", """
fluid { pvt { pressure (1.0e5 1.0e7); rho (55 55); mu (1.2e-5 1.2e-5); } }
nodes { W {}; Outlet { pressure 5.0e6; } }
pipes { P { from W; to Outlet; length 1000; diameter 0.1; roughness 4.6e-5;
    elevation 0; jouleThomson 1e-7; } }
sources { s { node W; rate 0.05; } }
solve { maxIter 100; tol 1e-10; }
""", "no temperature")
refuses("wall keys on a pure Joule-Thomson march",
        bo_dict(" jouleThomson 1e-7; heatCapacity 2000;"),
        "reads none of them")
refuses("an unknown wallModel word",
        bo_dict(" wallHeatTransfer 2; heatCapacity 2000; wallTemperature 350;"
                " wallModel spiral;"),
        "unknown wallModel")
refuses("a negative wallHeatTransfer",
        bo_dict(" wallHeatTransfer -2; heatCapacity 2000; wallTemperature 350;"),
        "negative wallHeatTransfer")

if TRANSIENT:
    other = open(os.path.join(HERE, "cases", "transientBetaMomentum",
                              "case.dict")).read()
    i = other.index("diameter")
    j = other.index(";", i) + 1
    refuses("thermal keys on the transient row",
            other[:j] + " jouleThomson 1e-7;" + other[j:],
            "transient", exe=TRANSIENT)

if NODAL:
    nd = open(os.path.join(HERE, "cases", "nodalVogel", "case.dict")).read()
    key = "    elevation 2200;    // m, wellhead above the pay"
    if key not in nd:
        raise SystemExit("probe: the nodalVogel well block moved")
    refuses("thermal keys on the nodal VLP well",
            nd.replace(key, key + "\n    wallHeatTransfer 2;"),
            "nodal: well declares a thermal parameter", exe=NODAL)

if fails:
    raise SystemExit(f"steady thermal march verification FAILED: {fails}")
print("the steady thermal march is passive: every declaration appends a T "
      "token and nothing else, the token is the re-derived segment-by-segment "
      "Joule-Thomson / wall relaxation on both marched paths, and every "
      "inert or out-of-scope declaration is refused")
