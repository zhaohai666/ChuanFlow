#!/usr/bin/env python3
"""5.157 (P5): independent verification of nodal{} with lateral{} multi-branch
bores (chflow-1d-nodal).

A multi-lateral case declares laterals { L1 { ipr{} well{} } L2 { ... } }:
every lateral carries its own IPR (linear / Vogel / Jones, forms may differ
across branches and against the main ipr{}) and its own Beggs-Brill VLP bore,
and its top opens into the main bore's bottom hole. The junction pressure p_j
is the fixed point p_j = p_vlp_main(sum q_i(p_j)), q_i(p_j) being lateral i's
own IPR x VLP crossing at back-pressure p_j.

This checker re-derives the whole chain from the case text alone, in Python,
on the same published physics the kernel transcribes: the Standing /
Vazquez-Beggs / DAK black-oil state, the Beggs-Brill flow-pattern map with its
L1..L4 boundaries and holdup + inclination branches, the Churchill/Colebrook
friction multiplier, the 24-segment march, the closed-form IPRs and their
inverses, the expanding-bracket + Illinois-secant VLP inversion and nodal
crossing, and the junction fixed point. It asserts:

  * every flowing lateral's printed operating point satisfies its own IPR
    closed form recomputed here from the printed bottom-hole pressure,
  * every lateral's printed rate is reproduced INDEPENDENTLY by the Python
    IPR x VLP solve run at the printed junction pressure (a dead branch is
    independently confirmed dead),
  * the junction is a fixed point: the Python main-bore VLP at the printed
    total lands on the printed junction pressure, and for case2 the full
    Python junction iteration converges to the printed total and junction,
  * mass conservation at the junction: sum(q_i) == Q at print precision,
  * the allocation responds: adding a third lateral raises the total and
    LOWERS the two original branches (the branches compete, not proportional),
  * flow-pattern map verdicts: every reported regime is a named Beggs-Brill
    regime, a dead branch reports rate 0 with its bottom hole exactly on its
    static column at the junction,
  * heterogeneous IPR forms across branches are honoured,
  * golden anchors recorded at calibration,
  * the undeclared path: the standing nodal case prints no lateral line,
  * eight loud refusals: empty laterals{}, a lateral without ipr{} / well{},
    flowPattern / correlation / zero length / unknown ipr form on a lateral,
    and the main-bore flowPattern message staying byte-identical.

usage: check_nodal_lateral.py <chflow-1d-nodal> <case2.dict> <case3.dict>
       <caseHetero.dict> <caseDead.dict> <nodalVogel.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 7:
    raise SystemExit(__doc__)

EXE = os.path.abspath(sys.argv[1])
CASES = {name: os.path.abspath(sys.argv[i]) for i, name in
         enumerate(("case2", "case3", "caseHetero", "caseDead"), start=2)}
LEGACY = os.path.abspath(sys.argv[6])

fails = []
NCHECK = [0]


def check(name, cond, detail=""):
    NCHECK[0] += 1
    print(f"{'ok:  ' if cond else 'FAIL:'} {name} {detail}")
    if not cond:
        fails.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ---- case-text parsing ----------------------------------------------------
def strip_comments(src):
    return re.sub(r"//[^\n]*", "", src)


def _span(src, brace_idx):
    depth = 0
    for j in range(brace_idx, len(src)):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                return brace_idx, j
    raise SystemExit("case: unbalanced block")


def block_body(src, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", src)
    if not m:
        return None
    a, b = _span(src, src.index("{", m.start()))
    return src[a + 1:b]


def cut_block(src, key):
    """src with the first `key { ... }` chunk removed."""
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", src)
    a, b = _span(src, src.index("{", m.start()))
    return src[:m.start()] + src[b + 1:]


def named_children(src):
    """[(name, body)] of the NAME { ... } blocks at the top level of src."""
    out = []
    pat = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\s*\{")
    i, n = 0, len(src)
    while i < n:
        m = pat.match(src, i)
        if m:
            a, b = _span(src, src.index("{", m.end() - 1))
            out.append((m.group(0).split("{")[0].strip(), src[a + 1:b]))
            i = b + 1
        else:
            i += 1
    return out


def num(src, key, default=None):
    m = re.search(r"\b" + re.escape(key) + r"\s+([0-9.eE+-]+)\s*;", src)
    if m:
        return float(m.group(1))
    if default is None:
        raise SystemExit(f"case has no scalar '{key}'")
    return default


def word(src, key, default=None):
    m = re.search(r"\b" + re.escape(key) + r"\s+([A-Za-z_][A-Za-z0-9_]*)\s*;", src)
    return m.group(1) if m else default


# ---- physical constants, transcribed from src/pvt/src/blackOil.cpp --------
PSI = 6894.757293168361
SCF = 5.614583333333333
M_AIR = 28.9647e-3
R_GAS = 8.314462618
P_STD, T_STD = 101325.0, 288.15
G = 9.80665
RHO_WATER_STD = 999.014
LBFT3 = 16.0184634

# ---- the published Beggs-Brill flow-pattern map (the 5.85 machine) --------
AB = (0.98, 0.0, 0.845, 1.065)
BB_ = (0.4846, 0.0, 0.5351, 0.5824)
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
    return max(AB[r] * lam ** BB_[r] * nfr ** -CB[r], lam)


def blend_weight(lam, nfr):
    l1, l2, l3, _ = boundaries(lam)
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


def hl_lifted(hl, lam, nfr, vsl, rho_l, angle_deg, coeff, sig):
    if coeff == NONE or lam <= 1e-9:
        return hl
    d, e, f, h = coeff
    nlv = vsl * math.sqrt(math.sqrt(rho_l / (G * sig)))
    c = max((1.0 - lam) * math.log(d * lam ** e * nlv ** f * nfr ** h), 0.0)
    s1 = math.sin(1.8 * math.radians(angle_deg))
    return min(max(hl * (1.0 + c * (s1 - s1 ** 3 / 3.0)), 0.0), 1.0)


def hl_priced(lam, nfr, r, vsl, rho_l, angle_deg, sig):
    # the kernel clamps the final holdup to [0,1]; when the inclination
    # correction is NONE (horizontal or distributed regime uphill) the
    # unclamped hl_horizontal can exceed 1, so we clamp here
    hl = hl_lifted(hl_horizontal(lam, nfr, r), lam, nfr, vsl, rho_l,
                   angle_deg, incline_coeff(angle_deg, r), sig)
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
    """Independent restatement of evaluateBlackOil + marchBeggsBrillPipe."""

    def __init__(self, gg, rho_std, gor, wc, temp, sig):
        self.gg, self.rho_std, self.gor, self.wc = gg, rho_std, gor, wc
        self.T, self.sig = temp, sig
        self.api = 141.5 / (rho_std / RHO_WATER_STD) - 131.5
        self.tf = temp * 1.8 - 459.67
        self.tr = temp * 1.8
        self.gor_scf = gor * SCF
        self.wor = wc / (1.0 - wc) if 0.0 < wc < 0.999 else 0.0
        self.wd, self.wmu = 1000.0, 1e-3

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
              / (1e5 * max(pb, 300.0)))
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
            if hi - lo < 1e-12 * hi:
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
                    mug=mu_gas)

    def march(self, q, pin, well):
        """24-segment march; returns (pOutlet, ok, holdupOutlet, regimeMax)."""
        L, D = well["length"], well["diameter"]
        area = 0.25 * math.pi * D * D
        nseg = 24
        ds = L / nseg
        p_floor = max(0.02 * pin, 1000.0)
        sin_theta = well["elevation"] / L
        angle_deg = math.degrees(math.asin(max(-1.0, min(1.0, sin_theta))))
        p = pin
        hl_out, r_max = 0.0, -1
        ok = True
        for _ in range(nseg):
            st = self.state(p)
            vsl = q * (st["bo"] + self.wor) / area
            vsg = q * st["fg"] * st["bg"] / area
            vm = vsl + vsg
            lam = vsl / vm if vm > 0.0 else 0.0
            nfr = vm * vm / (G * D)
            r = regime_of(lam, nfr) if vm > 0.0 else 0
            if vm <= 0.0:  # the rate-0 static column: no-slip, no slip loss
                hl, rho_s = 0.0, st["rhol"]
            else:
                hl = hl_priced(lam, nfr, r, vsl, st["rhol"], angle_deg, self.sig)
                rho_s = st["rhol"] * hl + st["rhog"] * (1.0 - hl)
            rho_m = st["rhol"] * lam + st["rhog"] * (1.0 - lam)
            mu_m = st["mul"] * lam + st["mug"] * (1.0 - lam)
            re = rho_m * vm * D / mu_m if vm > 0.0 else 0.0
            fns = colebrook(re, well["roughness"] / D) if vm > 0.0 else 0.0
            if vm <= 0.0 or lam < 1e-9 or hl < 1e-9:
                sm = 0.0
            else:
                yv = lam / (hl * hl)
                ly = math.log(yv)
                if 1.0 < yv < 1.2:
                    sm = math.log(2.2 * yv - 1.2)
                else:
                    sm = (ly / (ly * (ly * (0.01853 * ly * ly - 0.8725)
                                      + 3.182) - 0.0523))
                sm = min(sm, 7.0)
            ftp = fns * math.exp(sm)
            ek = vsg * vm * rho_s / max(p, 1000.0)
            ele = G * math.sin(math.radians(angle_deg)) * rho_s
            dpds = (ele + ftp * rho_m * vm * vm / (2 * D)) / max(1 - ek, 0.1)
            hl_out, r_max = hl, max(r_max, r)
            p -= dpds * ds
            if not p > p_floor:
                p = p_floor
                ok = False
                break
        return p, ok, hl_out, r_max


def make_fluid(fl_body):
    return BlackOil(num(fl_body, "gasGravity", 0.65), num(fl_body, "oilDensity"),
                    num(fl_body, "gor", 0.0), num(fl_body, "waterCut", 0.0),
                    num(fl_body, "temperature", 300.0),
                    num(fl_body, "surfaceTension", 0.03))


def make_well(w_body):
    return dict(length=num(w_body, "length"), diameter=num(w_body, "diameter"),
                roughness=num(w_body, "roughness", 0.0),
                elevation=num(w_body, "elevation", 0.0))


# ---- the closed-form IPRs (mirrors of IpR in Nodal.h) --------------------
class Ipr:
    def __init__(self, body):
        self.pres = num(body, "reservoirPressure")
        self.form = word(body, "form", "vogel")
        if self.form == "linear":
            self.J = num(body, "productivityIndex")
        elif self.form == "vogel":
            self.qmax = num(body, "maxRate")
        elif self.form == "jones":
            self.C = num(body, "deliverability")
            self.n = num(body, "exponent", 1.0)
        else:
            raise SystemExit(f"case ipr form '{self.form}' unknown")

    def q_cap(self):
        if self.form == "linear":
            return self.J * self.pres
        if self.form == "vogel":
            return self.qmax
        return self.C * self.pres ** (2.0 * self.n)

    def pwf(self, q):
        q = max(0.0, min(q, self.q_cap()))
        if self.form == "linear":
            return self.pres - q / self.J
        if self.form == "vogel":
            x = (-0.2 + math.sqrt(3.24 - 3.2 * q / self.qmax)) / 1.6
            return self.pres * x
        arg = self.pres * self.pres - (q / self.C) ** (1.0 / self.n)
        return math.sqrt(max(arg, 0.0))

    def rate(self, p):
        p = max(0.0, min(p, self.pres))
        if self.form == "linear":
            return self.J * (self.pres - p)
        if self.form == "vogel":
            x = p / self.pres
            return self.qmax * (1.0 - 0.2 * x - 0.8 * x * x)
        return self.C * (self.pres * self.pres - p * p) ** self.n


# ---- the VLP inversion and the nodal crossing (mirrors of nodal.cpp) ------
def vlp(fl, well, q, pwh):
    def g(pin):
        po, ok, _, _ = fl.march(q, pin, well)
        return (po - pwh) if ok else max(0.02 * pin, 1000.0) - pwh

    p_lo = pwh
    f_lo = g(p_lo)
    if f_lo > -1e-9 * pwh:
        return p_lo
    p_hi = p_lo + max(1e4, 0.05 * pwh)
    f_hi = g(p_hi)
    while f_hi < 0.0:
        p_lo, f_lo = p_hi, f_hi
        p_hi *= 2.0
        f_hi = g(p_hi)
    p = p_lo + (p_hi - p_lo) * f_lo / (f_lo - f_hi)
    for _ in range(100):
        fp = g(p)
        if abs(fp) < 1e-9 * pwh:
            break
        if fp > 0.0:
            p_hi, f_hi = p, fp * 0.5
        else:
            p_lo, f_lo = p, fp * 0.5
        pn = p_lo + (p_hi - p_lo) * f_lo / (f_lo - f_hi)
        p = pn if p_lo < pn < p_hi else 0.5 * (p_lo + p_hi)
        if abs(p_hi - p_lo) < 1e-9 * p_hi:
            break
    return p


def solve_nodal(ipr, fl, well, pwh):
    """Returns (q, pwf, flows). The largest root of p_vlp(q) - p_ipr(q)."""
    static = vlp(fl, well, 0.0, pwh)
    qm = ipr.q_cap()

    def f(q):
        return vlp(fl, well, q, pwh) - ipr.pwf(q)

    q_hi, f_hi = qm, f(qm)
    q_lo, f_lo, br = 0.0, 0.0, False
    qp = q_hi
    for _ in range(60):
        qp *= 0.5
        fk = f(qp)
        if fk <= 0.0:
            q_lo, f_lo, br = qp, fk, True
            break
    if not br and f(0.0) <= 0.0:
        q_lo, f_lo, br = 0.0, f(0.0), True
    if not br:
        return 0.0, static, False
    p_scale = max(ipr.pres, pwh)
    q = q_lo + (q_hi - q_lo) * f_lo / (f_lo - f_hi)
    for _ in range(100):
        fq = f(q)
        if abs(fq) < 1e-6 * p_scale:
            break
        if fq > 0.0:
            q_hi, f_hi = q, fq * 0.5
        else:
            q_lo, f_lo = q, fq * 0.5
        qn = q_lo + (q_hi - q_lo) * f_lo / (f_lo - f_hi)
        q = qn if q_lo < qn < q_hi else 0.5 * (q_lo + q_hi)
        if abs(q_hi - q_lo) < 1e-9 * max(q_hi, 1e-30):
            break
    return q, vlp(fl, well, q, pwh), True


def junction_total(fl, main, pwh, laterals, p_j):
    return sum(solve_nodal(l["ipr"], fl, l["well"], p_j)[0] for l in laterals)


def junction_solve(fl, main, pwh, laterals):
    """The Python fixed point p_j = p_vlp_main(sum q_i(p_j))."""
    static = vlp(fl, main, 0.0, pwh)
    p_res_max = max(l["ipr"].pres for l in laterals)

    def f(p_j):
        return p_j - vlp(fl, main, junction_total(fl, main, pwh, laterals, p_j),
                         pwh)

    p_lo, f_lo = pwh, f(pwh)
    p_hi = max(2.0 * p_lo, 1.05 * static, 1.05 * p_res_max) + 1e4
    f_hi = f(p_hi)
    while f_hi <= 0.0:
        p_lo, f_lo = p_hi, f_hi
        p_hi *= 2.0
        f_hi = f(p_hi)
    p_scale = max(static, pwh, p_res_max)
    p = p_lo + (p_hi - p_lo) * f_lo / (f_lo - f_hi)
    for _ in range(60):
        fp = f(p)
        if abs(fp) < 1e-5 * p_scale:
            break
        if fp > 0.0:
            p_hi, f_hi = p, fp * 0.5
        else:
            p_lo, f_lo = p, fp * 0.5
        pn = p_lo + (p_hi - p_lo) * f_lo / (f_lo - f_hi)
        p = pn if p_lo < pn < p_hi else 0.5 * (p_lo + p_hi)
    return p, junction_total(fl, main, pwh, laterals, p)


# ---- run the solver and parse its readout ---------------------------------
def run_case(path):
    return subprocess.run([EXE, path], capture_output=True, text=True)


def parse_out(proc):
    out = {"laterals": []}
    for line in proc.stdout.splitlines():
        toks = line.replace("[chflow]", "").split()
        if not toks:
            continue
        if toks[0] == "column":
            out["static"] = float(toks[2])
        elif toks[0] == "lateral":
            out["laterals"].append(
                dict(name=toks[1], q=float(toks[3]), pwf=float(toks[5]),
                     holdup=float(toks[7]), regime=toks[9],
                     flows=toks[11] == "true"))
        elif toks[0] == "well" and toks[1] == "q":
            out["q"], out["pwf"], out["pth"] = (float(toks[2]), float(toks[4]),
                                                float(toks[6]))
        elif toks[0] == "well" and toks[1] == "holdup":
            out["holdup"], out["regime"] = float(toks[2]), toks[4]
        elif toks[0] == "nodal" and toks[1] == "solved:":
            out["iterations"] = int(toks[2].split("=")[1])
            out["flows"] = toks[3].split("=")[1] == "true"
            out["converged"] = "(converged)" in line and "NOT" not in line
    return out


def load_case(path):
    src = strip_comments(open(path).read())
    fl = make_fluid(block_body(src, "blackoil"))
    main = make_well(block_body(src, "well"))
    pwh = num(src, "wellheadPressure")
    laterals = []
    for name, body in named_children(block_body(src, "laterals")):
        laterals.append(dict(name=name, ipr=Ipr(block_body(body, "ipr")),
                             well=make_well(block_body(body, "well"))))
    return fl, main, pwh, laterals


# golden anchors recorded when this checker was calibrated; the solver is
# deterministic, so any drift means the physics moved.
GOLDEN_Q = {"case2": 0.00297038656391, "case3": 0.00363466635438,
            "caseHetero": 0.00274048090481, "caseDead": 0.00199074204004}

results, models = {}, {}
for name, path in CASES.items():
    proc = run_case(path)
    if proc.returncode != 0:
        raise SystemExit(f"solver failed on {name} rc={proc.returncode}\n"
                         f"{proc.stdout}\n{proc.stderr}")
    r = parse_out(proc)
    models[name] = load_case(path)
    check(f"{name}: the multi-lateral solve converges, the well flows and "
          f"every declared branch is reported",
          r["converged"] and r["flows"] and r["iterations"] > 0
          and len(r["laterals"]) == len(models[name][3]),
          f"iterations {r['iterations']}")
    results[name] = r

# ---- per-lateral physics ---------------------------------------------------
for name in CASES:
    r = results[name]
    fl, main, pwh, laterals = models[name]
    for lr, lat in zip(r["laterals"], laterals):
        ipr = lat["ipr"]
        if lr["flows"]:
            # the printed point sits on the branch's own closed-form inflow
            q_ipr = ipr.rate(lr["pwf"])
            check(f"{name}/{lr['name']}: the operating rate satisfies the "
                  f"{ipr.form} IPR identity recomputed here",
                  rel(lr["q"], q_ipr) < 1e-5,
                  f"q {lr['q']:.6g} vs {q_ipr:.6g}")
        # the branch rate IS its own IPR x VLP crossing at the junction
        # pressure - re-derived in Python, kernel internals untouched
        q_py, pwf_py, flows_py = solve_nodal(ipr, fl, lat["well"], r["pwf"])
        if lr["flows"]:
            check(f"{name}/{lr['name']}: the Python IPR x VLP solve at the "
                  f"junction pressure reproduces the branch rate and bottom "
                  f"hole",
                  flows_py and rel(q_py, lr["q"]) < 1e-4
                  and rel(pwf_py, lr["pwf"]) < 1e-5,
                  f"q {q_py:.6g} pwf {pwf_py:.6g}")
        else:
            # a dead branch: the Python solve agrees it is dead, and the
            # bottom hole sits exactly on its static column at the junction
            static_py = vlp(fl, lat["well"], 0.0, r["pwf"])
            check(f"{name}/{lr['name']}: the dead branch reports rate 0, the "
                  f"Python solve agrees it is dead and its bottom hole is "
                  f"its static column at the junction",
                  lr["q"] == 0.0 and not flows_py
                  and rel(lr["pwf"], static_py) < 1e-6,
                  f"pwf {lr['pwf']:.6g} vs column {static_py:.6g}")
    # the flow-pattern map verdicts across the whole completion
    verdicts = ([lr["regime"] for lr in r["laterals"] if lr["flows"]]
                + [r["regime"]])
    slips = ([lr["holdup"] for lr in r["laterals"] if lr["flows"]]
             + [r["holdup"]])
    check(f"{name}: every flowing bore reports a named Beggs-Brill regime "
          f"with a liquid holdup in (0,1]",
          all(v in REGIME_NAME for v in verdicts)
          and all(0.0 < h <= 1.0 for h in slips))

# ---- the junction: conservation, VLP membership, ordering ------------------
for name in CASES:
    r = results[name]
    fl, main, pwh, laterals = models[name]
    q_sum = sum(lr["q"] for lr in r["laterals"])
    check(f"{name}: mass conservation at the junction, sum(q_i) == Q",
          abs(q_sum - r["q"]) <= 1e-9 * max(r["q"], 1e-30),
          f"sum {q_sum:.12g} vs Q {r['q']:.12g}")
    p_vlp = vlp(fl, main, r["q"], pwh)
    check(f"{name}: the junction sits on the independently re-derived main "
          f"VLP at the total",
          rel(p_vlp, r["pwf"]) < 5e-6,
          f"p_vlp_main(Q) {p_vlp:.6g} vs p_wf {r['pwf']:.6g}")
    p_res_max = max(l["ipr"].pres for l in laterals)
    check(f"{name}: the junction lies between the wellhead and both the "
          f"largest branch reservoir and the main static column",
          pwh < r["pwf"] < min(p_res_max, r["static"]),
          f"{pwh:.3g} < {r['pwf']:.6g} < {min(p_res_max, r['static']):.6g}")

# ---- the full Python junction fixed point on case2 -------------------------
fl, main, pwh, laterals = models["case2"]
p_py, q_py = junction_solve(fl, main, pwh, laterals)
check("case2: the Python junction iteration converges to the printed total",
      rel(q_py, results["case2"]["q"]) < 1e-4,
      f"Q {q_py:.6g} vs {results['case2']['q']:.6g}")
check("case2: the Python junction iteration converges to the printed "
      "junction pressure",
      rel(p_py, results["case2"]["pwf"]) < 1e-4,
      f"p_j {p_py:.6g} vs {results['case2']['pwf']:.6g}")

# ---- the allocation responds to the completion -----------------------------
r2, r3 = results["case2"], results["case3"]
check("case3: the third branch raises the total and the junction together",
      r3["q"] > r2["q"] and r3["pwf"] > r2["pwf"],
      f"Q {r3['q']:.6g} > {r2['q']:.6g}, p_wf {r3['pwf']:.6g} > "
      f"{r2['pwf']:.6g}")
for nm in ("L1", "L2"):
    a = next(x for x in r2["laterals"] if x["name"] == nm)
    b = next(x for x in r3["laterals"] if x["name"] == nm)
    check(f"case3: adding L3 lowers {nm} (the branches compete, the "
          f"allocation is not proportional)",
          b["q"] < a["q"], f"{nm} {b['q']:.6g} < {a['q']:.6g}")

# ---- heterogeneous IPR forms are honoured ----------------------------------
rh = results["caseHetero"]
check("caseHetero: three different IPR forms (linear/jones/vogel) flow "
      "against the same junction",
      len(rh["laterals"]) == 3 and all(lr["flows"] for lr in rh["laterals"]))

# ---- the dead branch case ---------------------------------------------------
rd = results["caseDead"]
check("caseDead: the healthy branch carries the whole well",
      rd["laterals"][0]["q"] == rd["q"],
      f"L1 {rd['laterals'][0]['q']:.6g} == Q {rd['q']:.6g}")

# ---- golden anchors ----------------------------------------------------------
for name, ref in GOLDEN_Q.items():
    got = results[name]["q"]
    check(f"{name}: the total q {got:.12g} holds the golden anchor "
          f"{ref:.12g}", rel(got, ref) < 1e-9)

# ---- the undeclared path stays the single-well machine ----------------------
proc = run_case(LEGACY)
check("nodalVogel (undeclared twin): still runs and prints no lateral line",
      proc.returncode == 0 and "lateral" not in proc.stdout)

# ---- loud refusals -----------------------------------------------------------
def refuse(label, text, phrase):
    with tempfile.NamedTemporaryFile("w", suffix=".dict", delete=False) as fh:
        fh.write(text)
        path = fh.name
    try:
        proc = subprocess.run([EXE, path], capture_output=True, text=True)
        msg = proc.stdout + proc.stderr
        ok = proc.returncode != 0 and phrase in msg
        check(label, ok, "" if ok else msg[:200])
    finally:
        os.unlink(path)


base = strip_comments(open(CASES["case2"]).read())
lat_body = block_body(base, "laterals")
kids = dict(named_children(lat_body))
l1, l2 = kids["L1"], kids["L2"]
refuse("refuse: laterals{} with zero branches",
       base.replace(lat_body, " "), "declares no lateral")
refuse("refuse: a lateral without ipr{}",
       base.replace(l2, cut_block(l2, "ipr")),
       "lateral 'L2' has no ipr{} block")
refuse("refuse: a lateral without well{}",
       base.replace(l2, cut_block(l2, "well")),
       "lateral 'L2' has no well{} block")
l1w = block_body(l1, "well")
refuse("refuse: flowPattern on a lateral bore",
       base.replace(l1, l1.replace(l1w, l1w + " flowPattern 2;")),
       "lateral 'L1' well declares flowPattern")
refuse("refuse: correlation on a lateral bore",
       base.replace(l1, l1.replace(l1w, l1w + " correlation hagedornBrown;")),
       "lateral 'L1' well declares correlation")
l2w = block_body(l2, "well")
refuse("refuse: a lateral bore with zero length",
       base.replace(l2, l2.replace(l2w, l2w.replace("length    600;",
                                                    "length    0;"))),
       "lateral 'L2' well needs positive length")
l1i = block_body(l1, "ipr")
refuse("refuse: an unknown lateral ipr form",
       base.replace(l1, l1.replace(l1i, l1i.replace("form              "
                                                    "vogel;",
                                                    "form              "
                                                    "mueller;"))),
       "lateral 'L1' ipr form must be linear, vogel or jones")
mw = block_body(base, "well")
refuse("refuse: the main-bore flowPattern message stays byte-identical "
       "(no lateral prefix)",
       base.replace(mw, mw + " flowPattern 2;"),
       "nodal: well declares flowPattern")

print(f"\n{NCHECK[0]} nodal-lateral checks run")
if fails:
    print("FAILURES:")
    print("\n".join(fails))
    raise SystemExit(1)
print(f"all {NCHECK[0]} nodal-lateral checks passed")
