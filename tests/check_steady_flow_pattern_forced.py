#!/usr/bin/env python3
"""5.88 / D4 increment 2: the flow-pattern map as a WRITABLE dispatch surface.

Twelve black-oil pipes on the 5.85 fluid, rate and terminal. Every family is
an auto pipe plus its forced twins on the SAME geometry, declaring
`flowPattern` 0 / 2 / 3 to name the segment law that prices the holdup while
the map keeps being evaluated and reported.

The claim under test is that the override writes to the one dispatch surface
5.85 read from: it selects the holdup law and the published inclination
coefficient set, it leaves the map's own verdict legible, and it changes
nothing about lambda, Nfr, L1..L4 or the auto path. This checker re-derives
everything from the case dictionary alone -- the black-oil PVT, the map, the
three laws, the transition blend, the >= lambda clamp, the per-regime
inclination correction and the 24-segment march -- and asserts:

  * every pipe's reported outlet state (lambda, vm, HL, Nfr, regime, L1..L4)
    IS the independent march run with that pipe's own declaration, and the
    re-march lands on the terminal pressure without choking,
  * the auto arms are the shipped 5.85 arms (same geometry, same answer),
  * an override that AGREES with the map is byte-identical to its auto twin
    (H_Map vs H_Auto apart from the name),
  * an override that DISAGREES carries `mapRegime <map> FORCED`, and only such
    a pipe does -- the auto pipes carry no token at all,
  * the map's numbers survive the override: L1..L4 stay the published
    boundaries at the reported lambda, and lambda stays the no-slip fraction,
  * the disagreement moves the physics (holdup and dp differ from the twin),
  * a forced law is really the forced law: at Nfr ~ 400 the segregated law
    falls under the >= lambda clamp, so the forced-segregated arm on the
    distributed pipe sits exactly on lambda,
  * on the transition arm the auto answer is the published blend and forcing
    0 or 2 collapses the blend onto that one law,
  * uphill, the forced regime also switches the inclination coefficient set
    (the reported holdup matches ITS OWN set and not the other arms'), and
    uphill-distributed takes no correction at all,
  * downhill, one set prices every regime: the sign of the angle, not the
    regime, picks it,
  * a high-GOR arm the map puts in the lambda < 0.01 / Nfr < L1 segregated
    branch can be forced to distributed, and its lambda stays the map's,
  * flowPattern 1 (the blend is not a law), > 3, < -1, a fractional index and
    the key on a single-phase network are refused, as is the key on the
    transient row, the two-fluid row and the nodal VLP well - none of which is
    a network pipe with a Beggs-Brill map to override.

usage: check_steady_flow_pattern_forced.py <steady-exe> <case.dict> \
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

STEADY, CASE = sys.argv[1], os.path.abspath(sys.argv[2])
TRANSIENT = sys.argv[3] if len(sys.argv) > 3 else None
TWOFLUID = sys.argv[4] if len(sys.argv) > 4 else None
NODAL = sys.argv[5] if len(sys.argv) > 5 else None
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


def ints(src, key):
    m = re.search(r"\b" + re.escape(key) + r"\s+(-?\d+)\s*;", src)
    return int(m.group(1)) if m else None


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
    """The published horizontal holdup law of regime r, with the >= lambda
    floor. r == 1 is not a law of its own: it is lambda itself."""
    return max(AB[r] * lam ** BB[r] * nfr ** -CB[r], lam)


def blend_weight(lam, nfr):
    l1, l2, l3, l4 = boundaries(lam)
    return (l3 - nfr) / (l3 - l2)


def hl_horizontal(lam, nfr, r):
    if r == 1:
        w = blend_weight(lam, nfr)
        return w * hl0_of(lam, nfr, 0) + (1.0 - w) * lam
    return hl0_of(lam, nfr, r)


# The published inclination coefficient sets. Uphill takes one per REGIME and
# none at all for distributed; downhill takes one set whatever the regime.
UPHILL_R0 = (0.011, -3.768, 3.539, -1.614)
UPHILL_R2 = (2.96, 0.305, -0.4473, 0.0978)
DOWNHILL = (4.70, -0.3692, 0.1244, -0.5056)
NONE = "no-correction"


def incline_coeff(angle_deg, r):
    if abs(angle_deg) <= 1e-12:
        return NONE
    if angle_deg > 0.0:
        if r == 3:
            return NONE  # uphill distributed: the published set is absent
        return UPHILL_R0 if r == 0 else UPHILL_R2
    return DOWNHILL


def hl_lifted(hl, lam, nfr, vsl, rho_l, angle_deg, coeff):
    """Apply one coefficient set to one horizontal holdup."""
    if coeff == NONE or lam <= 1e-9:
        return hl
    d, e, f, h = coeff
    nlv = vsl * math.sqrt(math.sqrt(rho_l / (G * SIG)))
    c = max((1.0 - lam) * math.log(d * lam ** e * nlv ** f * nfr ** h), 0.0)
    s1 = math.sin(1.8 * math.radians(angle_deg))
    return min(max(hl * (1.0 + c * (s1 - s1 ** 3 / 3.0)), 0.0), 1.0)


def hl_priced(lam, nfr, r, vsl, rho_l, angle_deg, coeff=None):
    """The whole chain: the horizontal law of the PRICED regime r, then the
    inclination correction with the set that regime selects."""
    return hl_lifted(hl_horizontal(lam, nfr, r), lam, nfr, vsl, rho_l, angle_deg,
                     incline_coeff(angle_deg, r) if coeff is None else coeff)


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
    """Independent restatement of evaluateBlackOil + marchBeggsBrillPipe for
    one fluid (horizontal or inclined, no declared liquid rheology), marching
    under a forced pattern exactly as the solver does."""

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
        c1, c2, c3 = (0.0362, 1.0937, 25.724) if self.api <= 30 else (0.0178, 1.187, 23.931)
        return c1 * self.gg * max(p_psia, 0.0) ** c2 * math.exp(c3 * self.api / (self.tf + 459.67))

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
        rho_oil = (350 * (self.rho_std / RHO_WATER_STD)
                   + 0.0764 * self.gg * rs * SCF) / (SCF * b) * LBFT3
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

    def march(self, q, pin, L, D, rough, sintheta=0.0, forced=-1, nseg=24):
        area = 0.25 * math.pi * D * D
        ds = L / nseg
        p_floor = max(0.02 * pin, 1000.0)
        angle_deg = math.degrees(math.asin(max(-1.0, min(1.0, sintheta))))
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
            r_map = regime_of(lam, nfr)
            r = r_map if forced < 0 else forced
            l1, l2, l3, l4 = boundaries(lam)
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
                    sm = ly / (ly * (ly * (0.01853 * ly * ly - 0.8725) + 3.182) - 0.0523)
                sm = min(sm, 7.0)
            ftp = fns * math.exp(sm)
            ek = vsg * vm * rhos / max(p, 1000.0)
            ele = G * math.sin(math.radians(angle_deg)) * rhos
            dpds = (ele + ftp * rhom * vm * vm / (2 * D)) / max(1 - ek, 0.1)
            rmin, rmax = min(rmin, r), max(rmax, r)
            last = dict(p=p, lam=lam, nfr=nfr, r=r, r_map=r_map, hl=hl, vm=vm,
                        vsl=vsl, rhol=st["rhol"], l1=l1, l2=l2, l3=l3, l4=l4,
                        ftp=ftp, angle=angle_deg)
            p -= dpds * ds
            if not p > p_floor:
                p = p_floor
                ok = False
                break
        return p, last, ok, rmin, rmax


FLUID = BlackOil(GG, RHO_STD, GOR, WC, T, WATER_DENSITY, WATER_MU)
MASS_PER_SM3 = FLUID.mass_per_sm3


# ---- run the solver and parse its readout ----
def parse_kv(toks):
    """The readout is `key value` pairs, but the value is a word for the
    regime keys and the FORCED marker is a third word, so a positional zip
    desynchronises the tail of the line. Walk it instead: a token is a key
    when it does not parse as a number."""
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
            d["mapRegime"] = toks[i + 1]  # ... <name> FORCED
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


def field(p, key):
    return float(p[key])


nodes, pipes, lines = run(CASE)

pipe_geo = {}
for name, body in entries(block(text, "pipes")).items():
    decl = ints(body, "flowPattern")
    pipe_geo[name] = dict(D=num(body, "diameter"), L=num(body, "length"),
                          rough=num(body, "roughness"),
                          elev=num(body, "elevation"),
                          forced=-1 if decl is None else decl,
                          frm=word(body, "from"), to=word(body, "to"))
src_rate = {}
for name, body in entries(block(text, "sources")).items():
    src_rate[word(body, "node")] = num(body, "rate")

check("twelve pipes present", len(pipes) == 12, f"{len(pipes)} pipes")

for name in sorted(pipes):
    f, geo = pipes[name], pipe_geo[name]
    D, L, rough, forced = geo["D"], geo["L"], geo["rough"], geo["forced"]
    p_up, p_dn = nodes[geo["frm"]], nodes[geo["to"]]
    q = src_rate[geo["frm"]]
    sintheta = geo["elev"] / L
    lam, nfr, vm = field(f, "lambda"), field(f, "Nfr"), field(f, "v")
    hl = field(f, "HL")
    oreg, envelope = f["oRegime"], f["regime"]
    L1, L2, L3, L4 = (field(f, k) for k in ("L1", "L2", "L3", "L4"))
    tokens = lines[name].split()

    pout, last, ok, rmin, rmax = FLUID.march(
        q, p_up, L, D, rough, sintheta=sintheta, forced=forced)
    r_map, r_priced = regime_of(lam, nfr), last["r"]

    # --- the reported state IS the march run with this pipe's declaration ---
    check(f"{name}: re-march does not choke", ok)
    check(f"{name}: re-march lands on terminal", near(pout, p_dn, 1e-6),
          f"{pout:.3f} vs {p_dn:.3f}")
    check(f"{name}: outlet lambda is the re-march's", near(lam, last["lam"], 1e-8),
          f"{lam:.9f} vs {last['lam']:.9f}")
    check(f"{name}: outlet velocity is the re-march's", near(vm, last["vm"], 1e-8))
    check(f"{name}: outlet holdup is the re-march's", near(hl, last["hl"], 1e-8),
          f"{hl:.9f} vs {last['hl']:.9f}")
    check(f"{name}: Nfr = vm^2/(gD)", near(nfr, vm * vm / (G * D), 1e-9))
    check(f"{name}: oRegime is the PRICED regime", oreg == REGIME_NAME[r_priced],
          f"{oreg} vs {REGIME_NAME[r_priced]}")
    check(f"{name}: regime is the max priced along the pipe",
          envelope == REGIME_NAME[rmax], f"{envelope} vs {REGIME_NAME[rmax]}")
    check(f"{name}: mass = q*surface mass", near(field(f, "m"), q * MASS_PER_SM3, 1e-6))
    check(f"{name}: dp positive downstream", field(f, "dp") > 0.0)
    check(f"{name}: lambda in (0,1)", 0.0 < lam < 1.0, f"{lam:.6f}")
    check(f"{name}: holdup in (0,1]", 0.0 < field(f, "HL") <= 1.0, f"{hl:.6f}")

    # --- the map is still evaluated and still reported under the override ---
    bb = boundaries(lam)
    check(f"{name}: L1..L4 stay the published boundaries at the reported lambda",
          all(near(a, b, 1e-9) for a, b in zip((L1, L2, L3, L4), bb)))
    check(f"{name}: lambda stays the no-slip fraction",
          near(lam, last["vsl"] / vm, 1e-9))

    # --- the FORCED token is exactly the disagreement footprint ---
    has_map = "mapRegime" in f
    want_token = r_priced != r_map
    check(f"{name}: mapRegime token iff the override disagrees",
          has_map == want_token,
          f"token={has_map} priced={REGIME_NAME[r_priced]} map={REGIME_NAME[r_map]}")
    if has_map:
        check(f"{name}: mapRegime names the map's verdict and FORCED follows it",
              f["mapRegime"] == REGIME_NAME[r_map] and forced >= 0
              and oreg == REGIME_NAME[forced] and "FORCED" in tokens,
              f"mapRegime {f['mapRegime']} ... {oreg}")
    else:
        check(f"{name}: no override footprint (auto, or it agrees with the map)",
              forced < 0 or forced == r_map)

# ---- the auto arms are 5.85's path untouched ----
for auto in ("H_Auto", "T_Auto", "D_Auto", "I_Auto"):
    f = pipes[auto]
    check(f"{auto}: declares nothing, takes the map's answer",
          pipe_geo[auto]["forced"] < 0 and "mapRegime" not in lines[auto].split()
          and f["oRegime"] == REGIME_NAME[regime_of(field(f, "lambda"),
                                                    field(f, "Nfr"))])

# The three horizontal auto arms share 5.85's P2 / P3 / P5 geometry exactly:
# the auto path of the writable surface is the shipped read-only one.
map_case = os.path.join(HERE, "cases", "steadyFlowPatternMap", "case.dict")
if os.path.exists(map_case):
    _, mp, ml = run(map_case)
    for twin, shipped in (("H_Auto", "P3"), ("T_Auto", "P2"), ("D_Auto", "P5")):
        same = all(near(field(pipes[twin], k), field(mp[shipped], k), 1e-9)
                   for k in ("lambda", "Nfr", "HL", "v", "L1", "L2", "L3", "L4", "dp"))
        check(f"{twin} is 5.85's {shipped} to 1e-9", same,
              f"HL {pipes[twin]['HL']} vs {mp[shipped]['HL']}, "
              f"dp {pipes[twin]['dp']} vs {mp[shipped]['dp']}")
else:
    raise SystemExit(f"the 5.85 witness case is missing: {map_case}")

# ---- an override that AGREES is byte-identical to its auto twin ----
def anonymise(line, name):
    return line.replace("pipe " + name, "pipe X")


anon_map = anonymise(lines["H_Map"], "H_Map")
anon_auto = anonymise(lines["H_Auto"], "H_Auto")
check("H_Map (flowPattern 2) is byte-identical to H_Auto", anon_map == anon_auto,
      "the map's own verdict, written down, changes nothing"
      if anon_map == anon_auto else f"{anon_map}\n vs {anon_auto}")

# ---- an override that DISAGREES moves the answer ----
for twin in ("H_Seg", "H_Dist"):
    check(f"{twin}: the forced law moves the holdup",
          pipes[twin]["HL"] != pipes["H_Auto"]["HL"],
          f"{field(pipes[twin], 'HL'):.6f} vs {field(pipes['H_Auto'], 'HL'):.6f}")
    check(f"{twin}: ... and the pressure drop",
          pipes[twin]["dp"] != pipes["H_Auto"]["dp"],
          f"{field(pipes[twin], 'dp'):.1f} vs {field(pipes['H_Auto'], 'dp'):.1f}")
    check(f"{twin}: ... but not lambda (the map's own number is untouched)",
          near(field(pipes[twin], "lambda"), field(pipes["H_Auto"], "lambda"), 1e-4))

# ---- a forced law is really the forced law ----
# At Nfr ~ 400 the segregated law 0.98 lam^0.4846 Nfr^-0.0868 falls under
# lambda, so the published >= lambda floor prices it: the forced-segregated
# arm on the distributed pipe sits EXACTLY on the no-slip holdup.
lam_d, nfr_d = field(pipes["D_Seg"], "lambda"), field(pipes["D_Seg"], "Nfr")
check("D_Seg: the segregated law alone is below lambda here",
      AB[0] * lam_d ** BB[0] * nfr_d ** -CB[0] < lam_d,
      f"raw {AB[0] * lam_d ** BB[0] * nfr_d ** -CB[0]:.9f} < lambda {lam_d:.9f}")
check("D_Seg: so the forced arm is clamped at lambda to the bit",
      abs(field(pipes["D_Seg"], "HL") - lam_d) <= 1e-12 * lam_d,
      f"HL {lam_d:.9f}")
check("D_Seg: the map still says distributed",
      pipes["D_Seg"]["oRegime"] == "segregated" and regime_of(lam_d, nfr_d) == 3
      and "mapRegime distributed FORCED" in lines["D_Seg"])

# ---- the transition arm: blend by default, one law when forced ----
lam_t, nfr_t = field(pipes["T_Auto"], "lambda"), field(pipes["T_Auto"], "Nfr")
w_t = blend_weight(lam_t, nfr_t)
check("T_Auto: the auto arm sits inside the map's L2..L3 band",
      0.0 < w_t < 1.0, f"w = {w_t:.6f}")
blend = w_t * hl0_of(lam_t, nfr_t, 0) + (1 - w_t) * lam_t
check("T_Auto: the auto holdup is the published blend",
      near(field(pipes["T_Auto"], "HL"), min(max(blend, 0.0), 1.0), 1e-8),
      f"{field(pipes['T_Auto'], 'HL'):.9f} vs {blend:.9f}")
for twin, law in (("T_Seg", 0), ("T_Int", 2)):
    lam_x, nfr_x = field(pipes[twin], "lambda"), field(pipes[twin], "Nfr")
    check(f"{twin}: forcing collapses the blend onto the r={law} law",
          near(field(pipes[twin], "HL"), hl0_of(lam_x, nfr_x, law), 1e-8),
          f"{field(pipes[twin], 'HL'):.9f} vs {hl0_of(lam_x, nfr_x, law):.9f}")
    check(f"{twin}: ... which is a different number than the blend",
          abs(field(pipes[twin], "HL") - field(pipes["T_Auto"], "HL")) > 1e-6)

# ---- uphill: the forced regime also switches the inclination set ----
for nm, r_want in (("I_Auto", 2), ("I_Seg", 0), ("I_Dist", 3)):
    geo, f = pipe_geo[nm], pipes[nm]
    pout, last, ok, _, _ = FLUID.march(src_rate[geo["frm"]], nodes[geo["frm"]],
                                       geo["L"], geo["D"], geo["rough"],
                                       sintheta=geo["elev"] / geo["L"],
                                       forced=geo["forced"])
    check(f"{nm}: priced with regime {r_want}'s own inclination set",
          last["r"] == r_want and near(field(f, "HL"), last["hl"], 1e-8),
          f"{field(f, 'HL'):.9f} vs {last['hl']:.9f}")
    others = [set_ for set_ in (UPHILL_R0, UPHILL_R2, DOWNHILL)
              if incline_coeff(last["angle"], r_want) != set_]
    if r_want == 3:
        plain = hl_horizontal(last["lam"], last["nfr"], 3)
        check(f"{nm}: uphill distributed takes NO inclination correction",
              near(field(f, "HL"), plain, 1e-9),
              f"lifted == the horizontal law {plain:.9f}")
    for set_ in others:
        alt = hl_priced(last["lam"], last["nfr"], r_want, last["vsl"],
                        last["rhol"], last["angle"], coeff=set_)
        check(f"{nm}: NOT priced with {set_}",
              abs(field(f, "HL") - alt) > 1e-6,
              f"{field(f, 'HL'):.9f} vs {alt:.9f}")
check("I_Seg: forcing the segregated law uphill lifts the holdup",
      field(pipes["I_Seg"], "HL") > field(pipes["I_Auto"], "HL"),
      f"{field(pipes['I_Seg'], 'HL'):.6f} > {field(pipes['I_Auto'], 'HL'):.6f}")
check("I_Auto: an uphill line is harder to push than the horizontal twin",
      field(pipes["I_Auto"], "dp") > field(pipes["H_Auto"], "dp"),
      f"{field(pipes['I_Auto'], 'dp'):.1f} > {field(pipes['H_Auto'], 'dp'):.1f}")

# ---- twins on their own dictionaries ----
OUT_P = nodes["Outlet"]


def one_pipe_dict(gor, wc, D, L, elev, decl, rate=0.004):
    return f"""
fluid {{ blackoil {{ gasGravity {GG}; oilDensity {RHO_STD}; gor {gor:g}; waterCut {wc};
    temperature {T}; surfaceTension {SIG}; }} }}
nodes {{ W {{}}; Outlet {{ pressure {OUT_P:.1f}; }} }}
pipes {{ P {{ from W; to Outlet; length {L}; diameter {D}; roughness 4.6e-5;
    elevation {elev};{decl} }} }}
sources {{ s {{ node W; rate {rate}; }} }}
solve {{ maxIter 100; tol 1e-10; }}
"""


def run_dict(body, exe=None):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "t.dict")
        with open(path, "w") as fh:
            fh.write(body)
        return run(path, exe)


# the lambda < 0.01 / Nfr < L1 branch, forced against the map
fluid_hi = BlackOil(GG, RHO_STD, 8000.0, 0.0, T, WATER_DENSITY, WATER_MU)
tn, tp, tl = run_dict(one_pipe_dict(8000.0, 0.0, 0.40, 500, 0, ""))
an, ap, al = run_dict(one_pipe_dict(8000.0, 0.0, 0.40, 500, 0, " flowPattern 3;"))
lam_g, nfr_g = field(tp["P"], "lambda"), field(tp["P"], "Nfr")
check("low-lambda twin: the map takes its lambda < 0.01 branch (Nfr < L1)",
      lam_g < 0.01 and tp["P"]["oRegime"] == "segregated"
      and nfr_g < field(tp["P"], "L1") and "mapRegime" not in tl["P"].split(),
      f"lambda={lam_g:.6f} Nfr={nfr_g:.6g} L1={field(tp['P'], 'L1'):.6g}")
lam_x, nfr_x = field(ap["P"], "lambda"), field(ap["P"], "Nfr")
check("low-lambda twin: the segregated branch can be forced to distributed",
      ap["P"]["oRegime"] == "distributed" and regime_of(lam_x, nfr_x) == 0
      and "mapRegime segregated FORCED" in al["P"], al["P"])
check("low-lambda twin: the forced arm prices the distributed law",
      near(field(ap["P"], "HL"), hl0_of(lam_x, nfr_x, 3), 1e-8),
      f"{field(ap['P'], 'HL'):.9f} vs {hl0_of(lam_x, nfr_x, 3):.9f}")
check("low-lambda twin: lambda stays the map's, not the declaration's",
      near(lam_x, lam_g, 1e-3), f"{lam_x:.6f} vs {lam_g:.6f}")
p2, l2s, ok2, _, _ = fluid_hi.march(0.004, an["W"], 500, 0.40, 4.6e-5, forced=3)
check("low-lambda twin: the re-march with the declaration lands on terminal",
      ok2 and near(p2, an["Outlet"], 1e-6)
      and near(field(ap["P"], "HL"), l2s["hl"], 1e-8))

# downhill: ONE coefficient set prices every regime. The drop is gentle
# (-20 m over 500 m) because a strongly assisting inclination takes the free
# node below the steady solve's convergence basin - that is a pre-existing
# property of the rate/drop secant, not of the dispatch surface, and the
# coefficient set only cares about the SIGN of the angle.
for decl, r_dn in (("", 2), (" flowPattern 0;", 0)):
    dn, dpp, dl = run_dict(one_pipe_dict(GOR, WC, 0.08, 500, -20, decl))
    pout, last, ok, _, _ = FLUID.march(0.004, dn["W"], 500.0, 0.08, 4.6e-5,
                                       sintheta=-20.0 / 500.0,
                                       forced=-1 if decl == "" else 0)
    check(f"downhill r={r_dn}: the re-march lands on terminal",
          ok and last["r"] == r_dn and last["angle"] < 0
          and near(pout, dn["Outlet"], 1e-6))
    check(f"downhill r={r_dn}: the published downhill set prices it",
          near(field(dpp["P"], "HL"), last["hl"], 1e-8),
          f"{field(dpp['P'], 'HL'):.9f} vs {last['hl']:.9f}")
    for set_ in (UPHILL_R0, UPHILL_R2):
        if r_dn == 0 and set_ == UPHILL_R0:
            continue
        alt = hl_priced(last["lam"], last["nfr"], r_dn, last["vsl"], last["rhol"],
                        last["angle"], coeff=set_)
        check(f"downhill r={r_dn}: NOT priced with the uphill set {set_}",
              abs(field(dpp["P"], "HL") - alt) > 1e-6,
              f"{field(dpp['P'], 'HL'):.9f} vs {alt:.9f}")
    check(f"downhill r={r_dn}: the map's verdict is untouched by the sign",
          (decl == "") == ("mapRegime" not in dl["P"].split()),
          f"map={REGIME_NAME[regime_of(last['lam'], last['nfr'])]}")

# ---- declaration scope: what the key refuses to be ----
src = open(CASE).read()
BASE = "    H_Auto { from WH1; to Outlet; length 500; diameter 0.08;  roughness 4.6e-5; elevation 0; }"
if BASE not in src:
    raise SystemExit("probe: the H_Auto arm moved out of the case")


def patched(decl):
    return src.replace(BASE, BASE[:-1] + " " + decl + " }")


def refuses(label, body, want):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "r.dict")
        with open(path, "w") as fh:
            fh.write(body)
        rr = subprocess.run([STEADY, path], capture_output=True, text=True, cwd=tmp)
    out = (rr.stdout + rr.stderr).strip()
    check(f"refused: {label}", rr.returncode != 0 and want in out,
          out.splitlines()[-1][:130] if out else "no message")


refuses("flowPattern 1 (the blend is not a law)", patched("flowPattern 1;"),
        "transition is not a holdup law")
refuses("flowPattern 4 (past the published regimes)", patched("flowPattern 4;"),
        "the published regimes are 0 segregated")
refuses("flowPattern -2 (below the auto sentinel)", patched("flowPattern -2;"),
        "the only negative value")
refuses("flowPattern 2.5 (an index, not a fraction)", patched("flowPattern 2.5;"),
        "fractional flowPattern")
refuses("flowPattern 2 on a single-phase network", """
fluid { rho 850; mu 0.1; }
nodes { W {}; Outlet { pressure 4.0e6; } }
pipes { P { from W; to Outlet; length 500; diameter 0.08; roughness 4.6e-5;
    elevation 0; flowPattern 2; } }
sources { s { node W; rate 0.004; } }
solve { maxIter 100; tol 1e-10; }
""", "single-phase")

# the key on the rows that have no Beggs-Brill map to override
for exe, case_name, tag in ((TRANSIENT, "transientBetaMomentum", "transient: pipe"),
                            (TWOFLUID, "twoFluidPressureLaminar", "twoFluid: pipe")):
    if not exe:
        continue
    other = open(os.path.join(HERE, "cases", case_name, "case.dict")).read()
    i = other.index("diameter")
    j = other.index(";", i) + 1
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "l.dict")
        with open(path, "w") as fh:
            fh.write(other[:j] + " flowPattern 2;" + other[j:])
        rr = subprocess.run([exe, path], capture_output=True, text=True, cwd=tmp)
    out = (rr.stdout + rr.stderr).strip()
    check(f"refused: flowPattern on the {case_name} row",
          rr.returncode != 0 and tag in out and "flowPattern" in out,
          out.splitlines()[-1][:130] if out else "no message")

# ... and on the nodal IPR-VLP well, whose march is not a network pipe
if NODAL:
    nd = open(os.path.join(HERE, "cases", "nodalVogel", "case.dict")).read()
    key = "    elevation 2200;    // m, wellhead above the pay"
    if key not in nd:
        raise SystemExit("probe: the nodalVogel well block moved")
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "n.dict")
        with open(path, "w") as fh:
            fh.write(nd.replace(key, key + "\n    flowPattern 2;"))
        rr = subprocess.run([NODAL, path], capture_output=True, text=True, cwd=tmp)
    out = (rr.stdout + rr.stderr).strip()
    check("refused: flowPattern on the nodal VLP well",
          rr.returncode != 0 and "nodal: well declares flowPattern" in out,
          out.splitlines()[-1][:130] if out else "no message")

if fails:
    raise SystemExit(f"forced flow-pattern dispatch verification FAILED: {fails}")
print("the flow-pattern map is writable: an agreed override is byte-identical "
      "to the auto answer and to 5.85, a disagreeing one reprices the holdup "
      "law and the inclination set while the map's own verdict stays on the "
      "readout")
