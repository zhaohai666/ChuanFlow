#!/usr/bin/env python3
"""5.85 / D4 increment 1: the Beggs-Brill flow-pattern map as a dispatch surface.

The shipped case is five horizontal black-oil pipes, one per regime of the
published map, all at the same stock-tank rate into a fixed-pressure terminal.
The diameter sets the outlet Froude number, so the five outlets land in
segregated / transition / intermittent / intermittent / distributed.

The claim under test is that the flow-pattern map is now a SINGLE dispatch
surface: beggsBrillDpDs reads its regime and the four L boundaries from
beggsBrillPattern, and the march reads out the outlet segment's map state.
This checker re-derives everything from the case dictionary alone -- the
black-oil PVT, the 24-segment march, and the published map -- and asserts:

  * the Newton solve converges, no pipe chokes,
  * the independent re-march from the reported upstream pressure lands on the
    terminal pressure, and its last segment IS the reported outlet state
    (lambda, vm, holdup, regime, L1..L4),
  * the reported Nfr is the kinematic identity vm^2/(g D),
  * the reported L1..L4 are the published boundaries at the reported lambda,
  * the reported outlet regime is the published precedence at (lambda, Nfr),
  * the along-pipe envelope brackets the outlet regime,
  * the transition arm's holdup is the published blend of the segregated
    holdup and lambda with weight (L3-Nfr)/(L3-L2),
  * the mass rate is the stock-tank rate times the surface mass per sm3,
  * all four regimes are visited,
  * a high-GOR twin (lambda < 0.01) takes the segregated branch Nfr < L1.

The dense standalone map (beggs_brill_pattern) is swept against an independent
Python restatement in tests/python/test_bindings.py; this script needs only
the steady executable.

usage: check_steady_flow_pattern_map.py <exe> <case.dict>
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


# ---- the black-oil correlations, restated from the case dictionary ----
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
T = num(fl, "temperature", 300.0)
SIG = num(fl, "surfaceTension", 0.025)
WATER_DENSITY = num(fl, "waterDensity", 1000.0)
WATER_MU = num(fl, "waterViscosity", 1e-3)

API = 141.5 / (RHO_STD / RHO_WATER_STD) - 131.5
TF = T * 1.8 - 459.67
TR = T * 1.8


class BlackOil:
    """Independent restatement of evaluateBlackOil + marchBeggsBrillPipe for
    one fluid (horizontal pipe, no declared liquid rheology)."""

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
        return max(18.2 * ((self.gor_scf / self.gg) ** 0.83 * 10 ** a - 1.4), 0.0)

    def rs_scf(self, p_psia):
        api = self.api
        c1, c2, c3 = (0.0362, 1.0937, 25.724) if api <= 30 else (0.0178, 1.187, 23.931)
        return c1 * self.gg * max(p_psia, 0.0) ** c2 * math.exp(c3 * api / (self.tf + 459.67))

    def bo_of(self, p_psia, rs, pb):
        api, gg, tf = self.api, self.gg, self.tf
        c1, c2, c3 = (4.670e-4, 1.100e-5, 1.337e-9) if api > 30 else (4.677e-4, 1.751e-5, -1.811e-8)
        tt = (tf - 60.0) * api / gg
        bob = 1 + c1 * rs + c2 * tt + c3 * rs * tt
        if p_psia <= pb:
            return bob
        co = (-1433 + 5 * rs + 17.2 * tf - 1180 * gg + 12.61 * api) / (1e5 * max(pb, 300.0))
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
        t1 = A[0] + A[1] / trr + A[2] / trr**3 + A[3] / trr**4 + A[4] / trr**5
        t2 = A[5] + A[6] / trr + A[7] / trr**2

        def zofrho(rho):
            r2 = rho * rho
            return (1 + t1 * rho + t2 * r2 - A[8] * (A[6] / trr + A[7] / trr**2) * rho**5
                    + A[9] * (1 + A[10] * r2) * r2 / trr**3 * math.exp(-A[10] * r2))

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
        gg, gor, T = self.gg, self.gor, self.T
        ppsia = max(p_pa, 1.0) / PSI
        pb = self.pb_psia()
        if p_pa <= pb * PSI:
            rs = min(self.rs_scf(ppsia) / SCF, gor)
            b = self.bo_of(ppsia, rs * SCF, pb)
            fg = gor - rs
        else:
            rs, fg = gor, 0.0
            b = self.bo_of(ppsia, self.gor_scf, pb)
        z = self.z_dak(ppsia)
        bg = (P_STD / p_pa) * (T / T_STD) * z
        rhog = p_pa * gg * M_AIR / (z * R_GAS * T)
        rho_oil = (350 * (self.rho_std / RHO_WATER_STD) + 0.0764 * gg * rs * SCF) / (SCF * b) * LBFT3
        x = 10 ** (3.0324 - 0.02023 * self.api) * self.tf ** -1.163
        mu_dead = 10 ** x - 1
        a = 10.715 * (rs * SCF + 100) ** -0.515
        bb = 5.44 * (rs * SCF + 150) ** -0.338
        mu_oil = a * mu_dead ** bb * 1e-3
        mw = gg * 28.9647
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

    def march(self, q, pin, L, D, rough, nseg=24):
        area = 0.25 * math.pi * D * D
        ds = L / nseg
        p_floor = max(0.02 * pin, 1000.0)
        p = pin
        last = None
        ok = True
        rmin, rmax = 9, -1
        for _ in range(nseg):
            st = self.state(p)
            vsl = q * (st["bo"] + st["wor"]) / area
            vsg = q * st["fg"] * st["bg"] / area
            vm = vsl + vsg
            lam = vsl / vm
            nfr = vm * vm / (G * D)
            r = regime_of(lam, nfr)
            l1, l2, l3, l4 = boundaries(lam)
            hl = holdup_of(lam, nfr)
            rhom = st["rhol"] * lam + st["rhog"] * (1 - lam)
            rhos = st["rhol"] * hl + st["rhog"] * (1 - hl)
            mum = st["mul"] * lam + st["mug"] * (1 - lam)
            re = rhom * vm * D / mum
            fns = colebrook(re, rough / D)
            if lam < 1e-9 or hl < 1e-9:
                s = 0.0
            else:
                yy = lam / (hl * hl)
                ly = math.log(yy)
                if 1.0 < yy < 1.2:
                    s = math.log(2.2 * yy - 1.2)
                else:
                    s = ly / (ly * (ly * (0.01853 * ly * ly - 0.8725) + 3.182) - 0.0523)
                s = min(s, 7.0)
            ftp = fns * math.exp(s)
            ek = vsg * vm * rhos / max(p, 1000.0)
            dpds = (ftp * rhom * vm * vm / (2 * D)) / max(1 - ek, 0.1)
            rmin, rmax = min(rmin, r), max(rmax, r)
            last = dict(p=p, lam=lam, nfr=nfr, r=r, hl=hl, vm=vm,
                        l1=l1, l2=l2, l3=l3, l4=l4, ftp=ftp)
            p -= dpds * ds
            if not p > p_floor:
                p = p_floor
                ok = False
                break
        return p, last, ok, rmin, rmax


fluid = BlackOil(GG, RHO_STD, GOR, WC, T, WATER_DENSITY, WATER_MU)
MASS_PER_SM3 = fluid.mass_per_sm3

# ---- the published flow-pattern map (independent restatement) ----
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


def holdup_of(lam, nfr):
    r = regime_of(lam, nfr)
    if r == 1:
        l1, l2, l3, l4 = boundaries(lam)
        w = (l3 - nfr) / (l3 - l2)
        hl = w * hl0_of(lam, nfr, 0) + (1.0 - w) * hl0_of(lam, nfr, 1)
    else:
        hl = hl0_of(lam, nfr, r)
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


# ---- run the solver and parse its readout ----
def run(case_path):
    proc = subprocess.run([exe, str(case_path)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"solver failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")
    if "steady solved" not in proc.stdout or "(converged)" not in proc.stdout:
        raise SystemExit("solver did not report convergence")
    nodes, pipes = {}, {}
    for line in proc.stdout.splitlines():
        toks = line.replace("[chflow]", "").split()
        if not toks:
            continue
        if toks[0] == "node":
            nodes[toks[1]] = float(toks[3])
        elif toks[0] == "pipe":
            if "CHOKED" in toks:
                raise SystemExit(f"pipe {toks[1]} choked: {line}")
            d = {}
            rest = toks[2:]
            for a, b in zip(rest[0::2], rest[1::2]):
                d[a] = b
            pipes[toks[1]] = d
    return nodes, pipes


nodes, pipes = run(case)

# geometry / sources from the dictionary
pipe_geo = {}
for name, body in entries(block(text, "pipes")).items():
    pipe_geo[name] = dict(D=num(body, "diameter"), L=num(body, "length"),
                          rough=num(body, "roughness"),
                          frm=word(body, "from"), to=word(body, "to"))
src_rate = {}
for name, body in entries(block(text, "sources")).items():
    src_rate[word(body, "node")] = num(body, "rate")

check("five pipes present", len(pipes) == 5, f"{sorted(pipes)}")
visited = set()
for name in sorted(pipes):
    f = pipes[name]
    geo = pipe_geo[name]
    D, L, rough = geo["D"], geo["L"], geo["rough"]
    p_up = nodes[geo["frm"]]
    p_dn = nodes[geo["to"]]
    q = src_rate[geo["frm"]]
    lam = float(f["lambda"]); nfr = float(f["Nfr"]); vm = float(f["v"])
    hl = float(f["HL"]); oreg = f["oRegime"]; reg = f["regime"]
    L1, L2, L3, L4 = (float(f[k]) for k in ("L1", "L2", "L3", "L4"))

    # independent re-march from the reported upstream pressure
    pout, last, ok, rmin, rmax = fluid.march(q, p_up, L, D, rough)
    check(f"{name}: re-march does not choke", ok)
    check(f"{name}: re-march lands on terminal", near(pout, p_dn, 1e-6),
          f"{pout:.3f} vs {p_dn:.3f}")
    check(f"{name}: outlet lambda is the re-march's", near(lam, last["lam"], 1e-8),
          f"{lam:.9f} vs {last['lam']:.9f}")
    check(f"{name}: outlet velocity is the re-march's", near(vm, last["vm"], 1e-8),
          f"{vm:.9f} vs {last['vm']:.9f}")
    check(f"{name}: outlet holdup is the re-march's", near(hl, last["hl"], 1e-8),
          f"{hl:.9f} vs {last['hl']:.9f}")
    check(f"{name}: Nfr = vm^2/(gD)", near(nfr, vm * vm / (G * D), 1e-9),
          f"{nfr:.9g}")
    check(f"{name}: Nfr is the re-march's", near(nfr, last["nfr"], 1e-8))
    # the readout is the published map at the reported outlet state
    b1, b2, b3, b4 = boundaries(lam)
    check(f"{name}: L1 = 316 lambda^0.302", near(L1, b1, 1e-9), f"{L1:.9g}")
    check(f"{name}: L2 = 0.0009252 lambda^-2.4684", near(L2, b2, 1e-9), f"{L2:.9g}")
    check(f"{name}: L3 = 0.1 lambda^-1.4516", near(L3, b3, 1e-9), f"{L3:.9g}")
    check(f"{name}: L4 = 0.5 lambda^-6.738", near(L4, b4, 1e-9), f"{L4:.9g}")
    check(f"{name}: L1..L4 are the re-march's",
          all(near(a, b, 1e-8) for a, b in zip((L1, L2, L3, L4),
              (last["l1"], last["l2"], last["l3"], last["l4"]))))
    r_map = regime_of(lam, nfr)
    check(f"{name}: oRegime is the map at (lambda,Nfr)", oreg == REGIME_NAME[r_map],
          f"{oreg} vs {REGIME_NAME[r_map]}")
    check(f"{name}: oRegime is the re-march's last regime", oreg == REGIME_NAME[last["r"]])
    check(f"{name}: envelope brackets outlet", reg in REGIME_NAME
          and REGIME_NAME.index(reg) >= r_map, f"max={reg} outlet={oreg}")
    check(f"{name}: lambda in (0,1)", 0.0 < lam < 1.0, f"{lam:.6f}")
    check(f"{name}: holdup in (0,1]", 0.0 < hl <= 1.0, f"{hl:.6f}")
    check(f"{name}: mass = q*surface mass", near(float(f["m"]), q * MASS_PER_SM3, 1e-6),
          f"{float(f['m']):.6f} vs {q * MASS_PER_SM3:.6f}")
    check(f"{name}: dp positive downstream", float(f["dp"]) > 0.0)
    if oreg == "transition":
        w = (L3 - nfr) / (L3 - L2)
        check(f"{name}: transition weight in (0,1)", 0.0 < w < 1.0, f"w={w:.6f}")
        hl_blend = w * hl0_of(lam, nfr, 0) + (1 - w) * lam
        check(f"{name}: transition holdup is the published blend",
              near(hl, min(max(hl_blend, 0), 1), 1e-8), f"{hl:.9f} vs {hl_blend:.9f}")
    visited.add(oreg)

check("all four regimes visited", visited == set(REGIME_NAME), f"{sorted(visited)}")

# ---- high-GOR twin: the lambda < 0.01 segregated branch ----
TWIN_GOR = 8000.0
twin_fluid = BlackOil(GG, RHO_STD, TWIN_GOR, 0.0, T, WATER_DENSITY, WATER_MU)
twin = f"""
fluid {{ blackoil {{ gasGravity {GG}; oilDensity {RHO_STD}; gor {TWIN_GOR:g}; waterCut 0;
    temperature {T}; surfaceTension {SIG}; }} }}
nodes {{ W {{}}; Outlet {{ pressure {nodes['Outlet']}; }} }}
pipes {{ P {{ from W; to Outlet; length 500; diameter 0.40; roughness 4.6e-5; elevation 0; }} }}
sources {{ s {{ node W; rate 0.004; }} }}
solve {{ maxIter 100; tol 1e-10; }}
"""
with tempfile.TemporaryDirectory() as tmp:
    twin_path = os.path.join(tmp, "twin.dict")
    with open(twin_path, "w") as fh:
        fh.write(twin)
    tn, tp = run(twin_path)
f = tp["P"]
lam = float(f["lambda"]); nfr = float(f["Nfr"])
L1 = float(f["L1"])
check("twin: lambda < 0.01", lam < 0.01, f"lambda={lam:.6f}")
check("twin: segregated via Nfr < L1", f["oRegime"] == "segregated" and nfr < L1,
      f"Nfr={nfr:.6g} L1={L1:.6g}")
check("twin: oRegime is the map", f["oRegime"] == REGIME_NAME[regime_of(lam, nfr)])
check("twin: L1 = 316 lambda^0.302", near(L1, 316 * lam ** 0.302, 1e-9))
pout, last, ok, _, _ = twin_fluid.march(0.004, tn["W"], 500, 0.40, 4.6e-5)
check("twin: re-march lands on terminal", ok and near(pout, tn["Outlet"], 1e-6))
check("twin: outlet lambda is the re-march's", near(lam, last["lam"], 1e-8))
check("twin: mass = q*surface mass", near(float(f["m"]), 0.004 * twin_fluid.mass_per_sm3, 1e-6))

if fails:
    raise SystemExit(f"flow-pattern map verification FAILED: {fails}")
print("flow-pattern map is the single dispatch surface: all regimes, readout,"
      " re-march and the lambda<0.01 branch agree")