#!/usr/bin/env python3
"""5.143 / G4: the correlation family as a second holdup dispatch surface.

The claim under test is that hagedornBrown and dunsRos replace ONLY the
horizontal + inclination holdup block of the Beggs-Brill gradient: the
flow-pattern map stays evaluated and reported, the friction multiplier,
the elevation and Ek terms and the 24-segment march read the same numbers,
and the undeclared path is bit for bit the shipped one. This checker
re-derives everything from the case dictionary alone -- the black-oil PVT,
the map, the Hagedorn-Brown (Whitson et al. 2017, Appendix E) and Duns-Ros
(API Bulletin 34 skeleton) holdup laws and the march -- and asserts:

  * the standing multiphase case with `correlation beggsBrill;` written on
    every pipe is BYTE-identical to the undeclared case (the declared
    default is the shipped path),
  * a single-phase constant-property network accepts the word as a no-op
    and is byte-identical without it,
  * on a horizontal anchor each alternative law's reported outlet state
    (lambda, HL, Nfr, L1..L4, dp) IS the independent 24-segment march run
    with that law, landing on the terminal pressure without choking,
  * the map survives the swap: L1..L4 stay the published boundaries at the
    reported lambda and oRegime stays the map's own verdict (nothing was
    forced),
  * the swap moves the physics: HL and dp differ from the Beggs-Brill twin,
  * the Hagedorn-Brown holdup sits on its no-slip floor exactly where the
    restatement's lower bound prices it,
  * an unknown correlation word, a correlation on a single-phase network,
    correlation + flowPattern on one pipe, and the key on the transient
    row and the nodal VLP well are all refused.

usage: check_steady_correlations.py <steady-exe> <case.dict> \
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
        # the published transition blends the BRANCH holdups of the two
        # neighbouring regimes (not lambda); the march below prices both
        # branches' full gradients
        w = blend_weight(lam, nfr)
        hl = w * hl0_of(lam, nfr, 0) + (1.0 - w) * hl0_of(lam, nfr, 2)
        return min(max(hl, 0.0), 1.0)
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


# ---- 5.143 (G4): the two alternative holdup laws, restated from the
# published sources exactly as the kernel transcribed them ----
def hl_hlb(lam, vsl, vsg, rho_l, mu_l, d, p_pa):
    """Whitson, Brule & Stephenson, JPT 2017 (SPE 169350) Appendix E."""
    if vsg <= 0.0:
        return 1.0
    if vsl <= 0.0:
        return 0.0
    rs = rho_l / (G * SIG)
    nlv = vsl * rs ** 0.25
    ngv = vsg * rs ** 0.25
    nd = d * math.sqrt(G * rho_l / SIG)
    nmu = mu_l * (G / (rho_l * SIG ** 3)) ** 0.25
    x1 = min(max(nmu, 0.002), 0.5)
    nmuc = (9.0950e-2 * x1 ** 3 - 1.2183e-1 * x1 ** 2 + 5.6348e-2 * x1
            + 1.9122e-3)
    x2 = ngv * nmu ** 0.380 / nd ** 2.14
    psi = 1.0
    if x2 >= 0.01:
        x = min(x2, 0.088)
        num_ = (3.6917e-1 * x ** 3 + 1.0029e-2 * x ** 2 - 1.2769e-3 * x
                + 1.8089e-5)
        den_ = (x ** 4 - 4.894e-2 * x ** 3 + 2.7058e-2 * x ** 2
                - 1.2769e-3 * x + 1.9585e-5)
        psi = num_ / den_
    x3 = (nlv / ngv ** 0.575) * (p_pa / 1.0e5) ** 0.10 * (nmuc / nd)
    if x3 <= 4.6423e-5:
        xx = max(x3, 2.1e-6)
        y3 = (-3.7902e-5 * xx + 4.326e-11) / (xx * xx - 1.9323e-4 * xx
                                              - 3.752e-10)
    elif x3 <= 5.0e-3:
        y3 = (-18.53 * x3 * x3 + 1.147 * x3 + 6.4898e-5) / (x3 + 4.5319e-4)
    else:
        y3 = 0.98
    return max(psi * y3, lam)


def hl_dr(lam, vsl, vsg, rho_l, mu_l, d):
    """Duns & Ros region skeleton, API Bulletin 34 Vol I Section C
    (field units). Transcription, not re-validation - see the kernel."""
    if vsg <= 0.0:
        return 1.0
    if vsl <= 0.0:
        return 0.0
    mu_cp = mu_l * 1000.0
    rho_lb = rho_l / LBFT3
    sig_dyn = SIG * 1000.0
    vsl_ft, vsg_ft, d_ft = vsl * FT, vsg * FT, d * FT
    nl = 0.15726 * mu_cp / (rho_lb * sig_dyn ** 3) ** 0.25
    nlv = 1.938 * vsl_ft * (rho_lb / sig_dyn) ** 0.25
    nd = 120.872 * d_ft * math.sqrt(rho_lb / sig_dyn)
    l1 = 13.2 * (1.0 - lam)
    l2 = l1 + 0.42 * nlv ** 0.25 * nl ** 0.10 * nd ** 0.25
    rl1 = 1.0 + 0.3 * (vsl_ft + vsg_ft)
    if vsg_ft < l1:
        hl = lam * rl1
    elif vsg_ft <= l2:
        hl = lam * (rl1 + (vsg_ft - l1) / (l2 - l1) * (1.0 - rl1))
    else:
        hl = lam
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


def bo_vb(p_psia, rs_scf, pb_psia, gg, api, tf):
    """The published Vazquez-Beggs Bo with the PUBLISHED undersaturated
    compressibility: co's denominator is the pressure the state is
    evaluated at (Vazquez & Beggs 1980), so co is re-priced at every p."""
    c1, c2, c3 = ((4.670e-4, 1.100e-5, 1.337e-9) if api > 30
                  else (4.677e-4, 1.751e-5, -1.811e-8))
    tt = (tf - 60.0) * api / gg
    bob = 1 + c1 * rs_scf + c2 * tt + c3 * rs_scf * tt
    if p_psia <= pb_psia:
        return bob
    co = (-1433 + 5 * rs_scf + 17.2 * tf - 1180 * gg + 12.61 * api) / (1e5 * p_psia)
    return bob * math.exp(co * (pb_psia - p_psia))


# third-party anchor for the published co form, computed by hand:
#   gg = 0.75, 30 API, TF = 160 F, Rs = 350 scf/stb
#   pb = 1736.2810988821 psia, p = pb + 1000 = 2736.2810988821 psia
#   Bob = 1.208381; co = 2562.3 / (1e5 * 2736.2810988821) = 9.3641694965e-6/psi
#   Bo = 1.208381 * exp(-9.3641694965e-6 * 1000) = 1.1971183305720536
check("published VB undersaturated Bo anchor: co evaluated AT p",
      near(bo_vb(2736.2810988821, 350.0, 1736.2810988821, 0.75, 30.0, 160.0),
           1.1971183305720536, 1e-9))


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
        return bo_vb(p_psia, rs, pb, self.gg, self.api, self.tf)

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
        sin_t = math.sin(math.radians(angle_deg))
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
                hl_seg = hl_priced(lam, nfr, 0, vsl, st["rhol"], angle_deg) \
                    if r == 1 else None
                hl_int = hl_priced(lam, nfr, 2, vsl, st["rhol"], angle_deg) \
                    if r == 1 else None
            rhom = st["rhol"] * lam + st["rhog"] * (1 - lam)
            mum = st["mul"] * lam + st["mug"] * (1 - lam)
            re = rhom * vm * D / mum
            fns = colebrook(re, rough / D)

            def branch(hl_b):
                """One published branch's full gradient at holdup hl_b."""
                hl_b = min(max(hl_b, 0.0), 1.0)
                rhos_b = st["rhol"] * hl_b + st["rhog"] * (1 - hl_b)
                if lam < 1e-9 or hl_b < 1e-9:
                    sm_b = 0.0
                else:
                    yy = lam / (hl_b * hl_b)
                    ly = math.log(yy)
                    if 1.0 < yy < 1.2:
                        sm_b = math.log(2.2 * yy - 1.2)
                    else:
                        sm_b = ly / (ly * (ly * (0.01853 * ly * ly - 0.8725)
                                           + 3.182) - 0.0523)
                    sm_b = min(sm_b, 7.0)
                f_b = fns * math.exp(sm_b)
                ek_b = vsg * vm * rhos_b / max(p, 1000.0)
                ele_b = G * sin_t * rhos_b
                return f_b, (ele_b + f_b * rhom * vm * vm / (2 * D)) \
                    / max(1 - ek_b, 0.1)

            if corr == 0 and r == 1:
                # the published transition: BOTH neighbouring branches in
                # full, blended with the published weight
                w = (l3 - nfr) / (l3 - l2)
                hl = w * hl_seg + (1 - w) * hl_int
                f_seg, dp_seg = branch(hl_seg)
                f_int, dp_int = branch(hl_int)
                ftp = w * f_seg + (1 - w) * f_int
                dpds = w * dp_seg + (1 - w) * dp_int
            else:
                if corr == 0:
                    hl = hl_priced(lam, nfr, r, vsl, st["rhol"], angle_deg)
                elif corr == 1:
                    hl = min(max(hl_hlb(lam, vsl, vsg, st["rhol"], st["mul"],
                                        D, p), 0.0), 1.0)
                else:
                    hl = min(max(hl_dr(lam, vsl, vsg, st["rhol"], st["mul"],
                                       D), 0.0), 1.0)
                ftp, dpds = branch(hl)
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


def field(p, key):
    return float(p[key])


# ---- the declared default is the shipped path, byte for byte ----
src = open(CASE).read()
pipes_blk = block(re.sub(r"//[^\n]*", "", src), "pipes")
declared = src
for m in re.finditer(r"(\b\w+\s*\{[^}]*?\})", pipes_blk):
    body = m.group(1)
    if "correlation" in body:
        raise SystemExit("probe: the standing case already declares correlation")
    patched = body[:-1] + " correlation beggsBrill; }"
    declared = declared.replace(body, patched)
base_out = run_stdout(src)
twin_out = run_stdout(declared)
check("every pipe with `correlation beggsBrill;` is byte-identical to the "
      "undeclared case", base_out == twin_out)

nodes, pipes, lines = run(CASE)

# ---- the anchor: one horizontal pipe, each law in turn ----
OUT_P = nodes["Outlet"]


def one_pipe_dict(decl, rate=0.004):
    return f"""
fluid {{ blackoil {{ gasGravity {GG}; oilDensity {RHO_STD}; gor {GOR:g};
    waterCut {WC}; temperature {T}; surfaceTension {SIG}; }} }}
nodes {{ W {{}}; Outlet {{ pressure {OUT_P:.1f}; }} }}
pipes {{ P {{ from W; to Outlet; length 800; diameter 0.1; roughness 4.6e-5;
    elevation 0;{decl} }} }}
sources {{ s {{ node W; rate {rate}; }} }}
solve {{ maxIter 100; tol 1e-10; }}
"""


def run_dict(body):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "t.dict")
        with open(path, "w") as fh:
            fh.write(body)
        return run(path)


nn, bb_pipes, bb_lines = run_dict(one_pipe_dict(""))
nn1, hlb_pipes, hlb_lines = run_dict(one_pipe_dict(" correlation hagedornBrown;"))
nn2, dr_pipes, dr_lines = run_dict(one_pipe_dict(" correlation dunsRos;"))

for tag, corr, got, gnn in (("hagedornBrown", 1, hlb_pipes, nn1),
                            ("dunsRos", 2, dr_pipes, nn2)):
    f = got["P"]
    lam, nfr, vm, hl = (field(f, k) for k in ("lambda", "Nfr", "v", "HL"))
    pout, last, ok = FLUID.march(0.004, gnn["W"], 800.0, 0.1, 4.6e-5,
                                 corr=corr)
    check(f"{tag}: re-march does not choke", ok)
    check(f"{tag}: re-march lands on terminal", near(pout, gnn["Outlet"], 1e-6),
          f"{pout:.3f} vs {gnn['Outlet']:.3f}")
    check(f"{tag}: outlet holdup is the re-march's",
          near(hl, last["hl"], 1e-8), f"{hl:.9f} vs {last['hl']:.9f}")
    check(f"{tag}: outlet lambda is the re-march's",
          near(lam, last["lam"], 1e-8))
    check(f"{tag}: dp is the re-march's",
          near(field(f, "dp"), gnn["W"] - pout, 1e-5),
          f"{field(f, 'dp'):.3f} vs {gnn['W'] - pout:.3f}")
    # the map survives the swap: evaluated, reported, unforced
    check(f"{tag}: L1..L4 stay the published boundaries at the reported lambda",
          all(near(a, b, 1e-9) for a, b in
              zip((field(f, k) for k in ("L1", "L2", "L3", "L4")),
                  boundaries(lam))))
    check(f"{tag}: oRegime is the map's own verdict (nothing was forced)",
          f["oRegime"] == REGIME_NAME[last["r_map"]]
          and "mapRegime" not in f)
    # the swap moves the physics against the Beggs-Brill twin
    check(f"{tag}: the holdup differs from the Beggs-Brill twin",
          abs(hl - field(bb_pipes["P"], "HL")) > 1e-4,
          f"{hl:.6f} vs {field(bb_pipes['P'], 'HL'):.6f}")
    check(f"{tag}: ... and so does the pressure drop",
          abs(field(f, "dp") - field(bb_pipes["P"], "dp")) > 100.0,
          f"{field(f, 'dp'):.1f} vs {field(bb_pipes['P'], 'dp'):.1f}")

# the Hagedorn-Brown lower bound: where the restatement prices psi*y3 under
# the no-slip content, the holdup sits on lambda exactly
lam_h = field(hlb_pipes["P"], "lambda")
check("hagedornBrown: the holdup respects the no-slip floor",
      field(hlb_pipes["P"], "HL") >= lam_h - 1e-12,
      f"HL {field(hlb_pipes['P'], 'HL'):.9f} >= lambda {lam_h:.9f}")

# ---- the Beggs-Brill transition: independent criteria on the default law ----
# (a) boundary continuity, closed form: at Nfr = L2 the published weight is
# 1 and the transition IS the segregated branch; at Nfr = L3 it is 0 and the
# transition IS the intermittent branch -- holdup and full gradient alike.
def _branch_gradient(lam, nfr, vsl, vsg, rhol, rhog, mul, mug, D, rough,
                     p_pa, hl_b):
    hl_b = min(max(hl_b, 0.0), 1.0)
    rho_m = rhol * lam + rhog * (1 - lam)
    mu_m = mul * lam + mug * (1 - lam)
    vm_b = vsl + vsg
    re = rho_m * vm_b * D / mu_m
    fns = colebrook(re, rough / D)
    rho_s = rhol * hl_b + rhog * (1 - hl_b)
    yy = lam / (hl_b * hl_b)
    ly = math.log(yy)
    if 1.0 < yy < 1.2:
        s = math.log(2.2 * yy - 1.2)
    else:
        s = ly / (ly * (ly * (0.01853 * ly * ly - 0.8725) + 3.182) - 0.0523)
    s = min(s, 7.0)
    f_b = fns * math.exp(s)
    ek_b = vsg * vm_b * rho_s / max(p_pa, 1000.0)
    return f_b, (f_b * rho_m * vm_b * vm_b / (2 * D)) / max(1 - ek_b, 0.1)


for lam_c in (0.05, 0.1, 0.2, 0.4, 0.6, 0.9):
    l2c, l3c = boundaries(lam_c)[1], boundaries(lam_c)[2]
    check(f"continuity lam={lam_c}: w(L2) = 1, w(L3) = 0",
          abs(blend_weight(lam_c, l2c) - 1.0) < 1e-12
          and abs(blend_weight(lam_c, l3c)) < 1e-12)
    check(f"continuity lam={lam_c}: holdup(L2) = segregated law",
          abs(hl_horizontal(lam_c, l2c, 1)
              - min(max(hl0_of(lam_c, l2c, 0), 0.0), 1.0)) < 1e-12)
    check(f"continuity lam={lam_c}: holdup(L3) = intermittent law",
          abs(hl_horizontal(lam_c, l3c, 1)
              - min(max(hl0_of(lam_c, l3c, 2), 0.0), 1.0)) < 1e-12)
    vm_c = math.sqrt(0.05 * G * 0.05)
    common = (lam_c * vm_c, (1 - lam_c) * vm_c, 915.0, 2.67, 180e-6, 14e-6,
              0.05, 0.0, 1.0e7)
    w_l2, w_l3 = blend_weight(lam_c, l2c), blend_weight(lam_c, l3c)
    _, dp_s2 = _branch_gradient(lam_c, l2c, *common, hl0_of(lam_c, l2c, 0))
    _, dp_i3 = _branch_gradient(lam_c, l3c, *common, hl0_of(lam_c, l3c, 2))
    _, dp_b2 = _branch_gradient(lam_c, l2c, *common,
                                w_l2 * hl0_of(lam_c, l2c, 0)
                                + (1 - w_l2) * hl0_of(lam_c, l2c, 2))
    _, dp_b3 = _branch_gradient(lam_c, l3c, *common,
                                w_l3 * hl0_of(lam_c, l3c, 0)
                                + (1 - w_l3) * hl0_of(lam_c, l3c, 2))
    check(f"continuity lam={lam_c}: w=1 collapse reaches the segregated branch",
          abs(dp_b2 - dp_s2) < 1e-12)
    check(f"continuity lam={lam_c}: w=0 collapse reaches the intermittent branch",
          abs(dp_b3 - dp_i3) < 1e-12)

# (c) third-party numeric anchor, computed by hand from the paper's
# formulas (horizontal, D = 0.05 m, p = 1e7 Pa, roughness 0: rhoL 915,
# rhoG 2.67, muL 180e-6, muG 14e-6; lambda = 0.4, Nfr = 0.05 so
# vm = sqrt(Nfr g D) = 0.15657785603334848 m/s):
#   L2 = 0.008881999440997143, L3 = 0.37813747212622595
#   w      = (L3 - 0.05)/(L3 - L2)            = 0.8886461986331783
#   HL_seg = 0.98*0.4^0.4846*0.05^-0.0868     = 0.8152930705289793
#   HL_int = 0.845*0.4^0.5351*0.05^-0.0173    = 0.5450380588971202
#   HL_t   = w*HL_seg + (1-w)*HL_int          = 0.7851991476453372
#   fns = 0.022536993703394008 (Re 35794.98)
#   f_t    = w*f_seg + (1-w)*f_int            = 0.029805001869711356
#   dp_t   = w*dp_seg + (1-w)*dp_int          = 2.6861370122741586 Pa/m
vm_a = math.sqrt(0.05 * G * 0.05)
lam_a, nfr_a = 0.4, 0.05
w_a = blend_weight(lam_a, nfr_a)
hl_seg_a, hl_int_a = hl0_of(lam_a, nfr_a, 0), hl0_of(lam_a, nfr_a, 2)
hl_a = w_a * hl_seg_a + (1 - w_a) * hl_int_a
common_a = (lam_a * vm_a, (1 - lam_a) * vm_a, 915.0, 2.67, 180e-6, 14e-6,
            0.05, 0.0, 1.0e7)
f_seg_a, dp_seg_a = _branch_gradient(lam_a, nfr_a, *common_a, hl_seg_a)
f_int_a, dp_int_a = _branch_gradient(lam_a, nfr_a, *common_a, hl_int_a)
check("anchor: transition weight", near(w_a, 0.8886461986331783, 1e-12),
      f"{w_a!r}")
check("anchor: segregated branch holdup",
      near(hl_seg_a, 0.8152930705289793, 1e-12), f"{hl_seg_a!r}")
check("anchor: intermittent branch holdup",
      near(hl_int_a, 0.5450380588971202, 1e-12), f"{hl_int_a!r}")
check("anchor: blended transition holdup",
      near(hl_a, 0.7851991476453372, 1e-12), f"{hl_a!r}")
check("anchor: blended transition fTp",
      near(w_a * f_seg_a + (1 - w_a) * f_int_a, 0.029805001869711356, 1e-12))
check("anchor: blended transition gradient",
      near(w_a * dp_seg_a + (1 - w_a) * dp_int_a, 2.6861370122741586, 1e-12))

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
        one_pipe_dict(" correlation turpin;"), "unknown correlation")
refuses("correlation + flowPattern on one pipe",
        one_pipe_dict(" correlation hagedornBrown; flowPattern 2;"),
        "honour one dispatch, not both")
refuses("correlation on a single-phase network", """
fluid { rho 850; mu 0.1; }
nodes { W {}; Outlet { pressure 4.0e6; } }
pipes { P { from W; to Outlet; length 500; diameter 0.08; roughness 4.6e-5;
    elevation 0; correlation hagedornBrown; } }
sources { s { node W; rate 0.004; } }
solve { maxIter 100; tol 1e-10; }
""", "single-phase network")

# the word on a single-phase network is only refused when it would price
# something: the declared default is a no-op there, byte-identical
single = """
fluid { rho 850; mu 0.1; }
nodes { W {}; Outlet { pressure 4.0e6; } }
pipes { P { from W; to Outlet; length 500; diameter 0.08; roughness 4.6e-5;
    elevation 0;%s } }
sources { s { node W; rate 0.004; } }
solve { maxIter 100; tol 1e-10; }
"""
check("correlation beggsBrill on a single-phase network is a no-op",
      run_stdout(single % " correlation beggsBrill;") == run_stdout(single % ""))

# the key on the rows that have no Beggs-Brill holdup block to replace
if TRANSIENT:
    other = open(os.path.join(HERE, "cases", "transientBetaMomentum",
                              "case.dict")).read()
    i = other.index("diameter")
    j = other.index(";", i) + 1
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "l.dict")
        with open(path, "w") as fh:
            fh.write(other[:j] + " correlation dunsRos;" + other[j:])
        rr = subprocess.run([TRANSIENT, path], capture_output=True, text=True,
                            cwd=tmp)
    out = (rr.stdout + rr.stderr).strip()
    check("refused: correlation on the transient row",
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
            fh.write(nd.replace(key, key + "\n    correlation hagedornBrown;"))
        rr = subprocess.run([NODAL, path], capture_output=True, text=True,
                            cwd=tmp)
    out = (rr.stdout + rr.stderr).strip()
    check("refused: correlation on the nodal VLP well",
          rr.returncode != 0 and "nodal: well declares correlation" in out,
          out.splitlines()[-1][:130] if out else "no message")

if fails:
    raise SystemExit(f"correlation family verification FAILED: {fails}")
print("the correlation family is a holdup dispatch: the declared default is "
      "byte-identical to the shipped path, each alternative law is the "
      "re-derived march on its own published holdup while the map stays "
      "evaluated and reported, and every out-of-scope declaration is refused")
