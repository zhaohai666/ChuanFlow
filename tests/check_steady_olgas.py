#!/usr/bin/env python3
"""5.155: the OLGAS-steady holdup law as a correlation-family dispatch.

The claim under test is that olgasSteady replaces ONLY the horizontal +
inclination holdup block of the Beggs-Brill gradient: the flow-pattern map
stays evaluated and reported and its OWN verdict dispatches the mechanistic
branches (flowPattern is refused alongside the law), the friction
multiplier, the elevation and Ek terms and the 24-segment march read the
same numbers, and the undeclared path is bit for bit the shipped one. The
law itself is a transcription of the published mechanistic skeleton of the
steady OLGA family (Taitel & Dukler 1976 lineage), NOT the proprietary
model - the checker re-derives the SAME skeleton from the case dictionary,
so it proves the code matches the transcription, not the transcription
matches OLGA. Everything is re-derived from the case dictionary alone - the
black-oil PVT, the map, the segment-geometry / slug-content / drift holdup
branches and the march - and asserts:

  * the standing olgasSteady case with the word rewritten to
    `correlation beggsBrill;` is BYTE-identical to the same case with the
    declaration deleted, and the deleted twin is BYTE-identical to the
    shipped steadyNetworkMultiphase anchor (5.143),
  * on the network each marched pipe's reported outlet state (lambda, HL)
    IS the independent 24-segment march run under the olgasSteady law from
    the reported inlet pressure at the cut-derived rate, landing on the
    reported outlet node pressure without choking,
  * the map survives the swap: L1..L4 stay the published boundaries at the
    reported lambda and oRegime stays the map's own verdict (nothing was
    forced),
  * the swap moves the physics: HL and dp differ from the Beggs-Brill twin
    on every pipe,
  * single-pipe horizontal anchors at rates chosen so the map's verdict
    lands in each regime exercise every branch of the law - segregated
    (the circular-segment inversion, cross-checked in closed form),
    transition (the map-weighted blend), intermittent (the slug-content
    weight against the Taitel-Dukler locus), distributed (the Zuber-Hill
    drift line, cross-checked in closed form) - plus an uphill segregated
    anchor where the (1 + sin theta) factor enters,
  * an unknown correlation word, olgasSteady on a single-phase network,
    olgasSteady + flowPattern on one pipe, and the key on the transient
    row and the nodal VLP well are all refused.

usage: check_steady_olgas.py <steady-exe> <case.dict> \
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

STEADY = os.path.abspath(sys.argv[1])
CASE = os.path.abspath(sys.argv[2])
TRANSIENT = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else None
NODAL = os.path.abspath(sys.argv[4]) if len(sys.argv) > 4 else None
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


def num(src, key, default=0.0):
    m = re.search(r"\b" + re.escape(key) + r"\s+([0-9.eE+-]+)\s*;", src)
    return float(m.group(1)) if m else default


# ---- the black-oil correlations, restated from the case dictionary ----
PSI = 6894.757293168361
SCF = 5.614583333333333
M_AIR = 28.9647e-3
R_GAS = 8.314462618
P_STD, T_STD = 101325.0, 288.15
G = 9.80665
RHO_WATER_STD = 999.014
LBFT3 = 16.0184634
FT = 3.28083989501312

fl = block(text, "blackoil")
GG = num(fl, "gasGravity", 0.65)
RHO_STD = num(fl, "oilDensity", 850.0)
GOR = num(fl, "gor", 0.0)
WC = num(fl, "waterCut", 0.0)
T = num(fl, "temperature", 300.0)
SIG = num(fl, "surfaceTension", 0.025)
WATER_DENSITY = num(fl, "waterDensity", 1000.0)
WATER_MU = num(fl, "waterViscosity", 1e-3)

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


def blend_weight(lam, nfr):
    l1, l2, l3, l4 = boundaries(lam)
    return (l3 - nfr) / (l3 - l2)


def hl_horizontal(lam, nfr, r):
    if r == 1:
        w = blend_weight(lam, nfr)
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


def hl_lifted(hl, lam, nfr, vsl, rho_l, angle_deg, coeff):
    if coeff == NONE or lam <= 1e-9:
        return hl
    d, e, f, h = coeff
    nlv = vsl * math.sqrt(math.sqrt(rho_l / (G * SIG)))
    c = max((1.0 - lam) * math.log(d * lam ** e * nlv ** f * nfr ** h), 0.0)
    s1 = math.sin(1.8 * math.radians(angle_deg))
    return min(max(hl * (1.0 + c * (s1 - s1 ** 3 / 3.0)), 0.0), 1.0)


def hl_priced(lam, nfr, r, vsl, rho_l, angle_deg):
    return hl_lifted(hl_horizontal(lam, nfr, r), lam, nfr, vsl, rho_l,
                     angle_deg, incline_coeff(angle_deg, r))


# ---- 5.155: the OLGAS-steady mechanistic skeleton, restated exactly as
# the kernel transcribed it (Taitel & Dukler 1976 lineage; see the
# provenance warning in olgasSteadyHoldup's comment - this proves the
# code matches the transcription, not the transcription matches OLGA) ----
def segment_area_fraction(level):
    # area fraction below the depth level*D of a circular bore (closed form)
    x = min(max(level, 0.0), 1.0)
    c = 1.0 - 2.0 * x
    return (math.acos(c) - c * math.sqrt(max(4.0 * x * (1.0 - x), 0.0))) / math.pi


def hl_olgas(lam, vsl, vsg, rho_l, rho_g, d, sin_theta, regime, tw):
    if vsg <= 0.0:
        return 1.0
    if vsl <= 0.0:
        return 0.0
    vm = vsl + vsg
    drho = max(rho_l - rho_g, 1e-6)

    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if segment_area_fraction(mid) < lam:
            lo = mid
        else:
            hi = mid
    hl_seg = min(max(0.5 * (lo + hi) * (1.0 + sin_theta), 0.0), 1.0)

    ugs = vsg * math.sqrt(rho_g / (G * d * drho))
    usl = vsl * math.sqrt(rho_l / (G * d * drho))
    ugsb = 3.1 * math.sqrt(max(101.0 - 50.0 * usl, 0.0))
    beta = min(max(1.0 - ugs / ugsb, 0.0), 1.0) if ugsb > 0.0 else 0.0
    hl_int = min(max(beta + (1.0 - beta) * lam, 0.0), 1.0)

    hl_dist = min(max(1.0 - vsg / (1.2 * vm), lam), 1.0)

    if regime == 0:
        hl = hl_seg
    elif regime == 3:
        hl = hl_dist
    elif regime == 2:
        hl = hl_int
    else:
        hl = tw * hl_seg + (1.0 - tw) * hl_int
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


class BlackOil:
    """Independent restatement of evaluateBlackOil + marchBeggsBrillPipe,
    marching under a chosen holdup law exactly as the solver does."""

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

    def march(self, q, pin, L, D, rough, sintheta=0.0, forced=-1, nseg=24,
              corr=0):
        area = 0.25 * math.pi * D * D
        ds = L / nseg
        p_floor = max(0.02 * pin, 1000.0)
        angle_deg = math.degrees(math.asin(max(-1.0, min(1.0, sintheta))))
        p = pin
        last = None
        ok = True
        for _ in range(nseg):
            st = self.state(p)
            vsl = q * (st["bo"] + st["wor"]) / area
            vsg = q * st["fg"] * st["bg"] / area
            vm = vsl + vsg
            lam = vsl / vm
            nfr = vm * vm / (G * D)
            r_map = regime_of(lam, nfr)
            r = r_map if forced < 0 else forced
            l1, l2, l3, l4 = boundaries(lam)
            if corr == 0:
                hl = hl_priced(lam, nfr, r, vsl, st["rhol"], angle_deg)
            else:
                tw = blend_weight(lam, nfr) if r == 1 else 1.0
                hl = hl_olgas(lam, vsl, vsg, st["rhol"], st["rhog"], D,
                              math.sin(math.radians(angle_deg)), r, tw)
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
            ek = vsg * vm * rhom / max(p, 1000.0)
            ele = G * math.sin(math.radians(angle_deg)) * rhos
            dpds = (ele + ftp * rhom * vm * vm / (2 * D)) / max(1 - ek, 0.1)
            last = dict(p=p, lam=lam, nfr=nfr, r=r, r_map=r_map, hl=hl, vm=vm,
                        l1=l1, l2=l2, l3=l3, l4=l4)
            p -= dpds * ds
            if not p > p_floor:
                p = p_floor
                ok = False
                break
        return p, last, ok


FLUID = BlackOil(GG, RHO_STD, GOR, WC, T, WATER_DENSITY, WATER_MU)


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


def run_stdout(body):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "t.dict")
        with open(path, "w") as fh:
            fh.write(body)
        proc = subprocess.run([STEADY, path], capture_output=True, text=True,
                              cwd=tmp)
    if proc.returncode != 0:
        raise SystemExit(f"solver failed\n{proc.stdout}\n{proc.stderr}")
    return proc.stdout


def run_dict(body):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "t.dict")
        with open(path, "w") as fh:
            fh.write(body)
        return run(path)


def field(p, key):
    return float(p[key])


# ---- the standing case declares the law on every marched pipe ----
if text.count("correlation olgasSteady;") != 3:
    raise SystemExit("probe: the anchor case must declare olgasSteady on "
                     "exactly its three pipes")
undeclared = text.replace(" correlation olgasSteady;", "")
declared_bb = text.replace("correlation olgasSteady;", "correlation beggsBrill;")

# ---- the declared default is the shipped path, byte for byte ----
check("every pipe with `correlation beggsBrill;` is byte-identical to the "
      "same case with the declaration deleted",
      run_stdout(declared_bb) == run_stdout(undeclared))
shipped = open(os.path.join(HERE, "cases", "steadyNetworkMultiphase",
                            "case.dict")).read()
check("the undeclared twin is byte-identical to the shipped "
      "steadyNetworkMultiphase anchor (5.143)",
      run_stdout(undeclared) == run_stdout(shipped))

nodes, pipes, lines = run(CASE)
OUT_P = nodes["Outlet"]

# ---- the network geometry and rates, re-derived from the case ----
geo = {}
for m in re.finditer(r"(\w+)\s*\{([^{}]*)\}", block(text, "pipes")):
    g = {}
    for kv in re.finditer(r"(\w+)\s+([A-Za-z0-9._+-]+);", m.group(2)):
        g[kv.group(1)] = kv.group(2)
    geo[m.group(1)] = g
srcs = []
for m in re.finditer(r"(\w+)\s*\{([^{}]*)\}", block(text, "sources")):
    srcs.append((re.search(r"node\s+(\w+);", m.group(2)).group(1),
                 float(re.search(r"rate\s+([0-9.eE+-]+);", m.group(2)).group(1))))


def cut_rate(edge, fr):
    """the rate through a tree pipe = the sources on the from side when the
    edge is cut - derived from the case, not read off the solver"""
    adj = {}
    for pn, g in geo.items():
        adj.setdefault(g["from"], []).append((g["to"], pn))
        adj.setdefault(g["to"], []).append((g["from"], pn))
    seen, stack = {fr}, [fr]
    while stack:
        n = stack.pop()
        for nb, pn in adj.get(n, []):
            if pn != edge and nb not in seen:
                seen.add(nb)
                stack.append(nb)
    return sum(r for nd, r in srcs if nd in seen)


# ---- (a) + (c): each marched pipe IS the independent olgasSteady march,
# with the map still evaluated and reported ----
for name, g in sorted(geo.items()):
    fr, to = g["from"], g["to"]
    L, D = float(g["length"]), float(g["diameter"])
    rough, elev = float(g["roughness"]), float(g["elevation"])
    forward = nodes[fr] >= nodes[to]
    pin = nodes[fr] if forward else nodes[to]
    pto = nodes[to] if forward else nodes[fr]
    sintheta = (elev / L) if forward else (-elev / L)
    pout, last, ok = FLUID.march(cut_rate(name, fr), pin, L, D, rough,
                                 sintheta=sintheta, corr=3)
    f = pipes[name]
    lam = field(f, "lambda")
    check(f"{name}: the re-march lands on the reported outlet node "
          f"without choking", ok and near(pout, pto, 1e-6),
          f"{pout:.3f} vs {pto:.3f}")
    check(f"{name}: outlet holdup is the re-march's",
          near(field(f, "HL"), last["hl"], 1e-8),
          f"{field(f, 'HL'):.9f} vs {last['hl']:.9f}")
    check(f"{name}: outlet lambda is the re-march's",
          near(lam, last["lam"], 1e-8))
    check(f"{name}: L1..L4 stay the published boundaries at the reported "
          f"lambda",
          all(near(a, b, 1e-9) for a, b in
              zip((field(f, k) for k in ("L1", "L2", "L3", "L4")),
                  boundaries(lam))))
    check(f"{name}: oRegime is the map's own verdict (nothing was forced)",
          f["oRegime"] == REGIME_NAME[last["r_map"]] and "mapRegime" not in f)

# ---- (b): the swap moves the physics against the Beggs-Brill twin ----
bb_nodes, bb_pipes, bb_lines = run_dict(undeclared)
for name in sorted(geo):
    check(f"{name}: the holdup differs from the Beggs-Brill twin",
          abs(field(pipes[name], "HL") - field(bb_pipes[name], "HL")) > 1e-4,
          f"{field(pipes[name], 'HL'):.6f} vs "
          f"{field(bb_pipes[name], 'HL'):.6f}")
    check(f"{name}: ... and so does the pressure drop",
          abs(field(pipes[name], "dp") - field(bb_pipes[name], "dp")) > 100.0,
          f"{field(pipes[name], 'dp'):.1f} vs "
          f"{field(bb_pipes[name], 'dp'):.1f}")


# ---- the branch anchors: single horizontal pipes at rates that put the
# map's own verdict in each regime, plus an uphill segregated bore ----
def one_pipe_dict(decl, rate=0.004, elev=0.0):
    return f"""
fluid {{ blackoil {{ gasGravity {GG}; oilDensity {RHO_STD}; gor {GOR:g};
    waterCut {WC}; temperature {T}; surfaceTension {SIG}; }} }}
nodes {{ W {{}}; Outlet {{ pressure {OUT_P:.1f}; }} }}
pipes {{ P {{ from W; to Outlet; length 800; diameter 0.1; roughness 4.6e-5;
    elevation {elev};{decl} }} }}
sources {{ s {{ node W; rate {rate}; }} }}
solve {{ maxIter 100; tol 1e-10; }}
"""


anchors = (("segregated", "segregated", 0.0002, 0.0, True, 1e-6),
           # the uphill segregated bore: the BB twin at this near-stagnant
           # incline sits below the shipped march's own convergence knee,
           # so this anchor compares against the horizontal twin of the SAME
           # law instead - the (1 + sin theta) factor must back liquid up.
           # Near the static-head knee the network's inverse march carries a
           # larger landing offset, so the terminal tolerance is loosened to
           # 1e-5 and the dp check below allows exactly that offset.
           ("segregated-uphill", "segregated", 0.0002, 100.0, False, 1e-5),
           ("transition", "transition", 0.0005, 0.0, True, 1e-6),
           ("intermittent", "intermittent", 0.004, 0.0, True, 1e-6),
           ("distributed", "distributed", 0.05, 0.0, True, 1e-6))
anchor_pipes = {}
for tag, expect, rate, elev, twin, land_tol in anchors:
    nn_o, o_pipes, _ = run_dict(one_pipe_dict(" correlation olgasSteady;",
                                             rate, elev))
    if twin:
        nn_b, b_pipes, _ = run_dict(one_pipe_dict("", rate, elev))
    anchor_pipes[tag] = o_pipes["P"]
    f = o_pipes["P"]
    pout, last, ok = FLUID.march(rate, nn_o["W"], 800.0, 0.1, 4.6e-5,
                                 sintheta=elev / 800.0, corr=3)
    check(f"{tag} anchor: re-march lands on terminal unchoked",
          ok and near(pout, OUT_P, land_tol), f"{pout:.3f} vs {OUT_P:.3f}")
    check(f"{tag} anchor: outlet holdup is the re-march's",
          near(field(f, "HL"), last["hl"], 1e-8),
          f"{field(f, 'HL'):.9f} vs {last['hl']:.9f}")
    # the re-march's dp against the reported one: the reported dp is the
    # solver's own W - Outlet, so the only gap the re-march can read is its
    # landing offset at the terminal - allowed for exactly, plus a relative
    # slack the tiny-dp segregated anchors could not otherwise resolve
    check(f"{tag} anchor: dp is the re-march's",
          abs(field(f, "dp") - (nn_o["W"] - pout))
          <= max(1.0, 1e-5 * field(f, "dp")) + abs(pout - OUT_P),
          f"{field(f, 'dp'):.3f} vs {nn_o['W'] - pout:.3f}")
    check(f"{tag} anchor: the map's own verdict exercises this branch",
          f["oRegime"] == expect and "mapRegime" not in f,
          f"oRegime {f['oRegime']}")
    if twin:
        check(f"{tag} anchor: the holdup moves against the Beggs-Brill twin",
              abs(field(f, "HL") - field(b_pipes["P"], "HL")) > 1e-4,
              f"{field(f, 'HL'):.6f} vs {field(b_pipes['P'], 'HL'):.6f}")
    else:
        check(f"{tag} anchor: the inclination factor backs liquid up "
              f"against the horizontal twin of the same law",
              field(f, "HL") > field(anchor_pipes["segregated"], "HL"),
              f"{field(f, 'HL'):.6f} vs "
              f"{field(anchor_pipes['segregated'], 'HL'):.6f}")

# the two closed-form branches, cross-checked without mirroring the
# bisection: the segregated holdup IS the level whose circular-segment area
# equals lambda, and the distributed holdup IS the Zuber-Hill drift line
seg_f = anchor_pipes["segregated"]
check("segregated: the holdup inverts the circular-segment geometry of the "
      "bore", near(segment_area_fraction(field(seg_f, "HL")),
                    field(seg_f, "lambda"), 1e-6),
      f"f(HL) {segment_area_fraction(field(seg_f, 'HL')):.9f} vs lambda "
      f"{field(seg_f, 'lambda'):.9f}")
dist_f = anchor_pipes["distributed"]
lam_d = field(dist_f, "lambda")
check("distributed: the holdup is the drift line 1 - (1 - lambda)/1.2",
      near(field(dist_f, "HL"), 1.0 - (1.0 - lam_d) / 1.2, 1e-6),
      f"{field(dist_f, 'HL'):.9f} vs {1.0 - (1.0 - lam_d) / 1.2:.9f}")


# ---- declaration scope ----
def refuses(label, body, want):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "r.dict")
        with open(path, "w") as fh:
            fh.write(body)
        rr = subprocess.run([STEADY, path], capture_output=True, text=True,
                            cwd=tmp)
    out = (rr.stdout + rr.stderr).strip()
    check(f"refused: {label}", rr.returncode != 0 and want in out,
          out.splitlines()[-1][:130] if out else "no message")


refuses("an unknown correlation word",
        one_pipe_dict(" correlation olgas;"), "unknown correlation")
refuses("olgasSteady on a single-phase network", """
fluid { rho 850; mu 0.1; }
nodes { W {}; Outlet { pressure 4.0e6; } }
pipes { P { from W; to Outlet; length 500; diameter 0.08; roughness 4.6e-5;
    elevation 0; correlation olgasSteady; } }
sources { s { node W; rate 0.004; } }
solve { maxIter 100; tol 1e-10; }
""", "single-phase network")
refuses("olgasSteady + flowPattern on one pipe",
        one_pipe_dict(" correlation olgasSteady; flowPattern 2;"),
        "honour one dispatch, not both")

if TRANSIENT:
    other = open(os.path.join(HERE, "cases", "transientBetaMomentum",
                              "case.dict")).read()
    i = other.index("diameter")
    j = other.index(";", i) + 1
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "l.dict")
        with open(path, "w") as fh:
            fh.write(other[:j] + " correlation olgasSteady;" + other[j:])
        rr = subprocess.run([TRANSIENT, path], capture_output=True, text=True,
                            cwd=tmp)
    out = (rr.stdout + rr.stderr).strip()
    check("refused: olgasSteady on the transient row",
          rr.returncode != 0 and "transient" in out and "correlation" in out,
          out.splitlines()[-1][:130] if out else "no message")

if NODAL:
    nd = open(os.path.join(HERE, "cases", "nodalVogel", "case.dict")).read()
    key = "    elevation 2200;    // m, wellhead above the pay"
    if key not in nd:
        raise SystemExit("probe: the nodalVogel well block moved")
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "n.dict")
        with open(path, "w") as fh:
            fh.write(nd.replace(key, key + "\n    correlation olgasSteady;"))
        rr = subprocess.run([NODAL, path], capture_output=True, text=True,
                            cwd=tmp)
    out = (rr.stdout + rr.stderr).strip()
    check("refused: olgasSteady on the nodal VLP well",
          rr.returncode != 0 and "nodal: well declares correlation" in out,
          out.splitlines()[-1][:130] if out else "no message")

if fails:
    raise SystemExit(f"olgasSteady verification FAILED: {fails}")
print("the OLGAS-steady skeleton is a holdup dispatch: the undeclared path "
      "stays the shipped march byte for byte, every branch the map's verdict "
      "selects is the re-derived mechanistic law, the map itself is only "
      "read and never forced, and every out-of-scope declaration is refused")
