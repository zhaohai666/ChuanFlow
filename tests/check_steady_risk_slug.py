#!/usr/bin/env python3
"""5.158 / P18 + 5.159 / P22: the steady slug-start forecast and the
erosion / hydrate / corrosion risk report.

Both features are PURE REPORTS on the Beggs-Brill march's own per-segment
profile: an undeclared case must not gain a single byte, and a declared
case's report lines must reproduce what an independent restatement of the
march says. This checker re-derives the whole black-oil march from the case
dictionary alone (the PVT, the Beggs-Brill gradient, the per-segment
holdup / no-slip content / mixture density / regime / velocity / pressure /
temperature profile) and then asserts:

  * byte-identity: stripping the slugPrediction{}/riskReport{} blocks from
    the case leaves the full stdout exactly equal to the declared stdout
    minus the slug/risk lines, and the standing multiphase case carries no
    slug/risk token anywhere;
  * the dip geometry: the trajectory's strict downhill -> uphill junction,
    the upstream march segment's regime, the pooled liquid volume over the
    uphill run, the level criterion and the SLUG-START verdict all match the
    re-derived profile;
  * the API RP 14E erosion line Ve = C / sqrt(rho_m) in field units, read at
    the worst vm/Ve segment, with the over-limit count;
  * the hydrate line: the declared Clausius-Clapeyron dissociation
    temperature at each march segment's pressure against the segment's
    temperature (the fluid's isothermal one, or the passive march's own in
    the thermal variant), the worst subcooling and its K bands;
  * the corrosion bands: CO2 / H2S partial pressures at the highest march
    pressure, the declared Pa bands and the 0.5 psia sour line;
  * the thermal variant's passivity: identical hydraulic, slug, erosion and
    corrosion bytes, with the hydrate subcooling moved by the march
    temperature alone;
  * seventeen refusals of the declaration grammar.

usage: check_steady_risk_slug.py <steady-exe> <case.dict> <caseThermal.dict>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) < 3:
    raise SystemExit(__doc__)

STEADY, CASE = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
CASE2 = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 else None
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


def hl_priced(lam, nfr, r, vsl, rho_l, angle_deg):
    # the declared trajectory never enters the march: the hydraulics price
    # the pipe at its single lumped elevation slope, horizontal here
    hl = hl_horizontal(lam, nfr, r)
    if abs(angle_deg) > 1e-12 and lam > 1e-9 and r != 3:
        coeff = ((4.70, -0.3692, 0.1244, -0.5056) if angle_deg < 0.0
                 else ((0.011, -3.768, 3.539, -1.614) if r == 0
                       else (2.96, 0.305, -0.4473, 0.0978)))
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


class BlackOil:
    """evaluateBlackOil + marchBeggsBrillPipe, restated. The march returns
    the per-segment profile the 5.158/5.159 reports read: inlet pressure,
    holdup, no-slip content, mixture density, regime and mixture velocity."""

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
        """The march segment by segment: each entry is the segment's INLET
        pressure, holdup, no-slip content, mixture velocity, no-slip mixture
        density and regime -- exactly the vectors the kernel pushes."""
        area = 0.25 * math.pi * D * D
        ds = L / nseg
        p_floor = max(0.02 * pin, 1000.0)
        angle_deg = math.degrees(math.asin(max(-1.0, min(1.0, sintheta))))
        p = pin
        segs = []
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
            segs.append(dict(p=pOld, hl=hl, lam=lam, vm=vm, rhom=rhom,
                              r=r, pNew=p))
            if not ok:
                break
        return segs, p, ok


def teq_of(p_pa, tref, pref, dh):
    return 1.0 / (1.0 / tref - (R_GAS / dh) * math.log(p_pa / pref))


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
        if t in REGIME_NAME or t in ("FORCED", "SLUG-START", "NO-SLUG",
                                     "OVER", "OK", "NONE", "LOW", "MEDIUM",
                                     "HIGH", "SOUR", "SWEET"):
            i += 1
            continue
        d[t] = toks[i + 1]
        i += 2
    return d


def run_full(case_body, exe=None):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "t.dict")
        with open(path, "w") as fh:
            fh.write(case_body)
        proc = subprocess.run([exe or STEADY, path], capture_output=True,
                              text=True, cwd=tmp)
    if proc.returncode != 0:
        raise SystemExit(f"solver failed rc={proc.returncode}\n{proc.stdout}\n{proc.stderr}")
    if "steady solved" not in proc.stdout or "(converged)" not in proc.stdout:
        raise SystemExit("solver did not report convergence")
    nodes, pipes, slugs, risks, out_lines = {}, {}, [], [], []
    for line in proc.stdout.splitlines():
        toks = line.replace("[chflow]", "").split()
        if not toks:
            continue
        out_lines.append(" ".join(toks))
        if toks[0] == "node":
            nodes[toks[1]] = float(toks[3])
        elif toks[0] == "pipe":
            pipes[toks[1]] = parse_kv(toks[2:])
        elif toks[0] == "slug":
            slugs.append(toks)
        elif toks[0] == "risk":
            risks.append(toks)
    return dict(nodes=nodes, pipes=pipes, slugs=slugs, risks=risks,
                lines=out_lines, raw=proc.stdout)


def without_block(src, key):
    m = re.search(r"\b" + re.escape(key) + r"\s*\{", src)
    if not m:
        raise SystemExit(f"case has no {key}{{}} block to strip")
    start = m.start()
    depth = 0
    j = m.end() - 1
    while j < len(src):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                break
        j += 1
    return src[:start] + src[j + 1:]


FLUID = BlackOil(GG, RHO_STD, GOR, WC, T_FLUID, WATER_DENSITY, WATER_MU)

# ---- the declared run ----
base = run_full(text)
p1_name = "P1"
p1_body = entries(block(text, "pipes"))[p1_name]
L1 = num(p1_body, "length")
D1 = num(p1_body, "diameter")
ROUGH = num(p1_body, "roughness")
ELEV = num(p1_body, "elevation")
AREA = 0.25 * math.pi * D1 * D1

slug_lines = [s for s in base["slugs"] if len(s) > 1 and s[1] == p1_name]
risk_lines = [s for s in base["risks"] if len(s) > 1 and s[1] == p1_name]
check("the declared case prints the slug forecast", bool(slug_lines))
check("the declared case prints the risk report", bool(risk_lines))

# ---- 1. byte-identity: strip the declarations, the rest is exact ----
stripped = without_block(without_block(text, "slugPrediction"), "riskReport")
st = run_full(stripped)
kept = [l for l in base["lines"]
        if not (l.startswith("slug ") or l.startswith("risk "))]
check("stripping the declarations removes exactly the slug/risk lines",
      kept == st["lines"],
      f"{len(kept)} vs {len(st['lines'])} lines")
check("the undeclared variant prints no slug/risk line",
      not any(l.startswith(("slug ", "risk ")) for l in st["lines"]))
check("the control pipe P2 prints no report",
      not any(s[1] == "P2" for s in base["slugs"] + base["risks"]))
mp = open(os.path.join(HERE, "cases", "steadyNetworkMultiphase",
                       "case.dict")).read()
mp_out = run_full(mp)
check("the standing multiphase case carries no slug/risk token",
      "slug" not in mp_out["raw"] and "risk" not in mp_out["raw"])

# ---- 2. the dip: re-derive the march and the pooled liquid ----
q1 = float(base["pipes"][p1_name]["Q"])
segs, p_end, ok = FLUID.march_profile(q1, base["nodes"]["W"], L1, D1, ROUGH,
                                      ELEV / L1)
check("the re-derived march lands on the solved downstream node",
      ok and near(p_end, base["nodes"]["M"], 1e-6),
      f"{p_end:.3f} vs {base['nodes']['M']:.3f}")
n = len(segs)
ds = L1 / n

sp = block(p1_body, "slugPrediction")
FILL = num(sp, "fillThreshold", 1.0)
traj = [(num(b, "length"), num(b, "angleDeg"))
        for b in entries(block(sp, "trajectory")).values()]
cum = [0.0]
for ln, _ in traj:
    cum.append(cum[-1] + ln)

dips = []
for k in range(len(traj) - 1):
    if traj[k][1] < 0.0 and traj[k + 1][1] > 0.0:
        s_d = cum[k + 1]
        smid = 0.5 * (cum[k] + cum[k + 1])
        im = min(max(int(smid / ds), 0), n - 1)
        up_lo, up_hi = s_d, cum[k + 2]
        vol = sum(segs[i]["hl"] * AREA
                  * (min(up_hi, (i + 1) * ds) - max(up_lo, i * ds))
                  for i in range(n)
                  if min(up_hi, (i + 1) * ds) - max(up_lo, i * ds) > 0.0)
        crit = vol / (AREA * traj[k + 1][0])
        dips.append(dict(s=s_d, vol=vol, crit=crit, regime=segs[im]["r"],
                         start=segs[im]["r"] == 0 and crit >= FILL))
dip_lines = [s for s in slug_lines if s[2] == "dip"]
check("the dip count matches the re-derived trajectory",
      len(dip_lines) == len(dips), f"{len(dip_lines)} vs {len(dips)}")
for k, d in enumerate(dips):
    toks = next(s for s in dip_lines if s[3] == str(k))
    check(f"dip {k}: arc position", near(float(toks[5]), d["s"], 1e-6))
    check(f"dip {k}: pooled volume",
          near(float(toks[7]), d["vol"], 1e-6),
          f"{float(toks[7]):.6f} vs {d['vol']:.6f}")
    check(f"dip {k}: level criterion",
          near(float(toks[9]), d["crit"], 1e-6))
    check(f"dip {k}: upstream regime word",
          toks[11] == REGIME_NAME[d["regime"]], toks[11])
    check(f"dip {k}: verdict",
          toks[12] == ("SLUG-START" if d["start"] else "NO-SLUG"), toks[12])
dips_count = next(s for s in slug_lines if s[2] == "dips")
check("the dips summary counts the trajectory's dips",
      int(dips_count[3]) == len(dips))
check("the case exercises a segregated upstream run",
      all(d["regime"] == 0 for d in dips) and bool(dips))
check("the criterion is a liquid level fraction in (0, 1]",
      all(0.0 < d["crit"] <= 1.0 + 1e-12 for d in dips))
check("the global slug-risk line lists exactly the SLUG-START pipes",
      any(s[1:] == ["risk", "pipes", "P1"] for s in base["slugs"]))

# a threshold above the pooled level flips the verdict to NO-SLUG and the
# global line to none
hi_src = text.replace("fillThreshold 0.5;", "fillThreshold 0.9;")
hi_out = run_full(hi_src)
hi_dip = next(s for s in hi_out["slugs"] if s[2] == "dip")
check("a higher fill threshold flips the verdict to NO-SLUG",
      hi_dip[-1] == "NO-SLUG", hi_dip[-1])
check("with no slug start the global line reads none",
      any(s[1:] == ["risk", "pipes", "none"] for s in hi_out["slugs"]))

# ---- 3. erosion: API RP 14E on the march's own mixture density ----
C14E = num(block(p1_body, "riskReport"), "erosionC", 100.0)
ratios = [s["vm"] / (0.3048 * C14E / math.sqrt(s["rhom"] / LBFT3))
          for s in segs]
worst = max(range(n), key=lambda i: ratios[i])
er = next(s for s in risk_lines if s[2] == "erosion")
check("erosion: worst segment", int(er[4]) == worst, f"{er[4]} vs {worst}")
check("erosion: vm at the worst segment",
      near(float(er[6]), segs[worst]["vm"], 1e-6))
check("erosion: the 14E velocity limit",
      near(float(er[8]), 0.3048 * C14E / math.sqrt(
          segs[worst]["rhom"] / LBFT3), 1e-6),
      f"{float(er[8]):.6f}")
check("erosion: the ratio", near(float(er[10]), ratios[worst], 1e-6))
check("erosion: the OVER/OK word",
      er[11] == ("OVER" if ratios[worst] > 1.0 else "OK"), er[11])
check("erosion: the over-limit count",
      int(er[13]) == sum(1 for r in ratios if r > 1.0))
check("erosion: the case sits under the limit", ratios[worst] < 1.0)

# ---- 4. hydrate: the declared dissociation line, isothermal variant ----
rr = block(p1_body, "riskReport")
TREF, PREF, DH = num(rr, "hydrateTref"), num(rr, "hydratePRef"), num(rr, "hydrateDH")
LOWK, MEDK = num(rr, "hydrateLowK", 3.0), num(rr, "hydrateMediumK", 8.0)
dT = [teq_of(s["p"], TREF, PREF, DH) - T_FLUID for s in segs]
hworst = max(range(n), key=lambda i: dT[i])
hy = next(s for s in risk_lines if s[2] == "hydrate")
check("hydrate: worst (largest subcooling) segment",
      int(hy[4]) == hworst, f"{hy[4]} vs {hworst}")
check("hydrate: segment pressure", near(float(hy[6]), segs[hworst]["p"], 1e-6))
check("hydrate: isothermal T is the fluid temperature",
      near(float(hy[8]), T_FLUID, 1e-12))
check("hydrate: dissociation temperature",
      near(float(hy[10]), teq_of(segs[hworst]["p"], TREF, PREF, DH), 1e-6))
check("hydrate: subcooling", near(float(hy[12]), dT[hworst], 1e-6),
      f"{float(hy[12]):.6f} vs {dT[hworst]:.6f}")
grade = ("NONE" if dT[hworst] < 0.0 else
         "LOW" if dT[hworst] < LOWK else
         "MEDIUM" if dT[hworst] < MEDK else "HIGH")
check("hydrate: grade word", hy[13] == grade, f"{hy[13]} vs {grade}")
check("hydrate: the warm isothermal line is NONE", grade == "NONE")

# ---- 5. corrosion: partial pressures at the highest march pressure ----
co2y, h2sy = num(rr, "co2Fraction"), num(rr, "h2sFraction")
co2lo, co2hi = num(rr, "co2LowPa", 1e4), num(rr, "co2HighPa", 1e5)
pmax = max(s["p"] for s in segs)
co = next(s for s in risk_lines if s[2] == "corrosion")
check("corrosion: CO2 partial pressure",
      near(float(co[4]), co2y * pmax, 1e-6), f"{float(co[4]):.6f}")
band = ("LOW" if co2y * pmax < co2lo else
        "MEDIUM" if co2y * pmax < co2hi else "HIGH")
check("corrosion: CO2 band", co[6] == band, f"{co[6]} vs {band}")
check("corrosion: H2S partial pressure",
      near(float(co[8]), h2sy * pmax, 1e-6))
check("corrosion: the 0.5 psia sour line",
      co[9] == ("SOUR" if h2sy * pmax >= 3447.0 else "SWEET"), co[9])
check("corrosion: the highest march pressure is the inlet",
      near(pmax, base["nodes"]["W"], 1e-9))

# ---- 6. the thermal variant: the march temperature moves the hydrate ----
if CASE2:
    text2 = re.sub(r"//[^\n]*", "", open(CASE2).read())
    th = run_full(text2)
    p1b2 = entries(block(text2, "pipes"))["P1"]
    U2, TW2, CP2 = (num(p1b2, "wallHeatTransfer"), num(p1b2, "wallTemperature"),
                    num(p1b2, "heatCapacity"))
    m2 = float(th["pipes"]["P1"]["m"])
    t = T_FLUID
    temps = []
    for s in segs:
        temps.append(t)
        t = TW2 + (t - TW2) * math.exp(-U2 * ds / (m2 * CP2))
    dT2 = [teq_of(s["p"], TREF, PREF, DH) - temps[i]
           for i, s in enumerate(segs)]
    w2 = max(range(n), key=lambda i: dT2[i])
    hy2 = next(s for s in th["risks"] if s[2] == "hydrate")
    # the thermal march appends its ` T <K>` token to the pipe line (5.144
    # passivity); with that token removed the hydraulic bytes are identical
    def hyd(lines):
        return [re.sub(r" T [0-9][0-9.eE+-]*", "", l)
                for l in lines if l.startswith(("pipe", "node"))]
    check("thermal: the hydraulic lines are byte-identical (passivity)",
          hyd(th["lines"]) == hyd(base["lines"]))
    check("thermal: the slug forecast is untouched",
          [l for l in th["lines"] if l.startswith("slug")]
          == [l for l in base["lines"] if l.startswith("slug")])
    check("thermal: erosion unchanged",
          next(s for s in th["risks"] if s[2] == "erosion")
          == er)
    check("thermal: corrosion unchanged",
          next(s for s in th["risks"] if s[2] == "corrosion") == co)
    check("thermal: the hydrate T is the march temperature, not the fluid's",
          abs(float(hy2[8]) - temps[w2]) < 1e-9 * max(temps[w2], 1.0)
          and abs(float(hy2[8]) - T_FLUID) > 1.0,
          f"{float(hy2[8]):.3f}")
    check("thermal: the subcooling is the re-derived march profile",
          near(float(hy2[12]), dT2[w2], 1e-6),
          f"{float(hy2[12]):.6f} vs {dT2[w2]:.6f}")
    check("thermal: the worst segment moved", int(hy2[4]) == w2)
    grade2 = ("NONE" if dT2[w2] < 0.0 else
              "LOW" if dT2[w2] < LOWK else
              "MEDIUM" if dT2[w2] < MEDK else "HIGH")
    check("thermal: the grade word", hy2[13] == grade2, hy2[13])
    check("thermal: the cold wall raises the subcooling above the isothermal",
          dT2[w2] > dT[hworst] and grade2 != "NONE",
          f"{dT2[w2]:.3f} > {dT[hworst]:.3f}")

# ---- 7. refusals ----
def refuses(label, body, want, want2=None):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "r.dict")
        with open(path, "w") as fh:
            fh.write(body)
        rr2 = subprocess.run([STEADY, path], capture_output=True,
                             text=True, cwd=tmp)
    out = (rr2.stdout + rr2.stderr).strip()
    ok_msg = rr2.returncode != 0 and want in out and (want2 is None or want2 in out)
    check(f"refused: {label}", ok_msg,
          out.splitlines()[-1][:130] if out else "no message")


CONST = """
fluid { rho 850; mu 0.001; }
nodes { W {}; Outlet { pressure 4.0e6; } }
pipes { P { from W; to Outlet; length 500; diameter 0.1; roughness 4.6e-5;
    elevation 0;%s } }
sources { s { node W; rate 0.001; } }
solve { maxIter 100; tol 1e-10; }
"""
TRJ = " trajectory { D1 { length 250; angleDeg -3; } D2 { length 250; angleDeg 3; } } "
refuses("a pipe-level trajectory key",
        text.replace("slugPrediction", "trajectory { D1 { length 3000; "
                    "angleDeg -3; } D2 { length 3000; angleDeg 3; } } "
                    "slugPrediction", 1),
        "belongs inside slugPrediction")
refuses("an unknown slugPrediction key",
        text.replace("fillThreshold 0.5;", "maxLevel 0.9;"),
        "unknown key 'maxLevel'")
refuses("slugPrediction on a constant-property network",
        CONST % ("slugPrediction { trajectory { D1 { length 250; angleDeg -3; } "
                 "D2 { length 250; angleDeg 3; } } }"),
        "without a black-oil march")
refuses("a one-segment trajectory",
        re.sub(r"\s*D2 \{[^}]*\}", "", text).replace(
            "D1 { length 3000; angleDeg -3; }",
            "D1 { length 6000; angleDeg 0; }"),
        "at least two named segments")
refuses("trajectory lengths that do not sum to the pipe length",
        text.replace("D2 { length 3000; angleDeg 3; }",
                     "D2 { length 2900; angleDeg 3; }"),
        "sum to")
refuses("a trajectory that does not end where the pipe does",
        text.replace("D2 { length 3000; angleDeg 3; }",
                     "D2 { length 3000; angleDeg 2; }"),
        "rises to")
refuses("a trajectory segment with a non-positive length",
        text.replace("D1 { length 3000; angleDeg -3; }",
                     "D1 { length -5; angleDeg -3; }"),
        "positive length")
refuses("a trajectory angle outside [-90, 90]",
        text.replace("D2 { length 3000; angleDeg 3; }",
                     "D2 { length 3000; angleDeg 120; }"),
        "outside [-90, 90]")
refuses("a non-positive fill threshold",
        text.replace("fillThreshold 0.5;", "fillThreshold 0;"),
        "fillThreshold must be positive")
refuses("an unknown trajectory segment key",
        text.replace("D1 { length 3000; angleDeg -3; }",
                     "D1 { length 3000; angleDeg -3; rise 5; }"),
        "a segment is length and angleDeg")
refuses("an unknown riskReport key",
        text.replace("erosionC     100;", "sandLoadRate 5;"),
        "unknown key 'sandLoadRate'")
refuses("riskReport on a constant-property network",
        CONST % " riskReport { erosionC 100; } ",
        "without a black-oil march")
refuses("a non-positive erosionC",
        text.replace("erosionC     100;", "erosionC 0;"),
        "erosionC must be positive")
partial = re.sub(r"\s*hydratePRef\s+[^;]+;", "", text)
partial = re.sub(r"\s*hydrateDH\s+[^;]+;", "", partial)
refuses("a partial hydrate line", partial,
        "must be declared together", "missing hydratePRef")
refuses("descending hydrate bands",
        text.replace("hydrateLowK    3;", "hydrateLowK    9;"),
        "the bands must rise")
refuses("a sour-gas fraction outside [0,1]",
        text.replace("co2Fraction  0.02;", "co2Fraction  1.5;"),
        "outside [0,1]")
refuses("descending CO2 bands",
        text.replace("co2LowPa     1.0e4;", "co2LowPa     2.0e5;"),
        "the bands must rise with the partial pressure")

if fails:
    raise SystemExit(f"steady slug/risk verification FAILED: {fails}")
print("the slug forecast and the risk report are pure readouts: every number "
      "reproduces the independently re-derived march profile, the undeclared "
      "byte stream is exact, and the declaration grammar refuses every inert "
      "or impossible shape")
