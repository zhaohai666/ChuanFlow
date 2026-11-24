#!/usr/bin/env python3
"""E8 verification: simplified van der Waals - Platteeuw hydrate equilibrium
(pvt layer, declared-input closure, gated on a `hydratePWHc{}` block).

usage: check_pvt_hydrate_pwhc.py <exe> <case.dict>

The checker re-derives every printed number independently in plain Python:

  PR flash       cubic roots, ln phi_i (PR76 omega_a / omega_b, van der Waals
                 one-fluid mixing, kij = 0), successive substitution on the
                 fugacity equality with Wilson K initialisation and the same
                 single-phase verdict rule as the kernel; nothing is taken on
                 trust from the C++ side.
  occupancy      theta_j = sum_i C_ij f_i / (1 + sum_i C_ij f_i), f_i the
                 flashed gas-side fugacity, C_ij the declared Langmuir
                 constants (case unit 1/MPa, converted to 1/Pa here).
  equilibrium    Delta_mu_w(T) + R T sum_j nu_j ln(1 - theta_j) = 0 with
                 Delta_mu_w(T) = muRef + muSlopeA (T - Tref) and nu_j the
                 unit-cell cages per water molecule (sI 46 H2O: 2 small + 6
                 large; sII 136 H2O: 16 small + 8 large). Bisection in p.

Physics gates: theta in (0,1) and strictly increasing in p; p_eq strictly
increasing in T; the Clausius-Clapeyron bound |d ln p/dT - dH/(R Tref^2)| <=
20% at the anchor. Twins: Langmuir -> 0 (C_ij x 1e-6) drives theta -> 0 at the
fixed anchor pressure and pushes the solved p_eq orders of magnitude right of
the base case (fugacity grows with p, so a smaller Langmuir constant still
finds a root below the search ceiling -- the ceiling is a domain guard, not an
expected outcome); a single-guest (pure methane) twin is compared against the
hand closed form theta = C f/(1 + C f);
a zero-declaration twin keeps the flash report byte-identical. Eight loud
rejections cover malformed declarations, including mutual exclusion with the
legacy `hydrateEquilibrium{}` block.
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

R = 8.314462618
SQRT2 = 1.4142135623730951
OMEGA_A = 0.45723552892138218
OMEGA_B = 0.07779607390388846

# PR76 alpha kappa (omega <= 0.491 branch): published PR76 standard constants,
# see Cantera 3.1 reference https://www.cantera.org/3.1/reference/thermo/species-thermo.html
PR_KAPPA_BASE = 0.37464
PR_KAPPA_LINEAR = 1.54226
PR_KAPPA_QUAD = 0.26992


def pr76_kappa(w):
    """PR76 standard kappa; independent of the chflow kernel constants."""
    return PR_KAPPA_BASE + PR_KAPPA_LINEAR * w - PR_KAPPA_QUAD * w * w


# Exact constants from src/pvt/src/components.cpp (MW g/mol, Tc K, Pc Pa, omega).
CT = {
    "methane": (16.043, 190.564, 45.992e5, 0.0115),
    "propane": (44.097, 369.834, 42.477e5, 0.1523),
}

# Declared-domain guard of the kernel (src/pvt/src/hydratePWHc.cpp).
P_FLOOR = 1.0e-3
P_CEILING = 1.0e10

checks = 0


def ok(name, cond, detail=""):
    global checks
    checks += 1
    assert cond, f"FAIL {name} {detail}"


def close(name, got, want, rtol):
    ok(name, abs(got - want) <= rtol * abs(want),
       f"got {got:.12g} want {want:.12g}")


def cbrt(x):
    return math.copysign(abs(x) ** (1.0 / 3.0), x)


def cubic_real_roots(c2, c1, c0):
    p = c1 - c2 * c2 / 3.0
    q = 2.0 * c2 * c2 * c2 / 27.0 - c2 * c1 / 3.0 + c0
    disc = q * q / 4.0 + p * p * p / 27.0
    shift = -c2 / 3.0
    if disc > 0.0:
        s = math.sqrt(disc)
        return [cbrt(-q / 2.0 + s) + cbrt(-q / 2.0 - s) + shift]
    if p >= 0.0:
        return [shift]
    r = math.sqrt(-p * p * p / 27.0)
    m = 2.0 * math.sqrt(-p / 3.0)
    theta = math.acos(max(-1.0, min(1.0, -q / (2.0 * r))))
    return sorted(m * math.cos((theta + 2.0 * k * math.pi) / 3.0) + shift
                  for k in range(3))


def pr_lnphi(names, zs, T, p, pref):
    """Independent PR state: ln phi list for the selected root."""
    n = len(names)
    RT = R * T
    aI, bI = [], []
    for c in names:
        _, Tc, Pc, w = CT[c]
        kappa = pr76_kappa(w)
        alpha = (1.0 + kappa * (1.0 - math.sqrt(T / Tc))) ** 2
        aI.append(OMEGA_A * (R * Tc) ** 2 / Pc * alpha)
        bI.append(OMEGA_B * R * Tc / Pc)
    aIJ = [[math.sqrt(aI[i] * aI[j]) for j in range(n)] for i in range(n)]
    bM = sum(zs[i] * bI[i] for i in range(n))
    aM = sum(zs[i] * zs[j] * aIJ[i][j] for i in range(n) for j in range(n))
    A = aM * p / (RT * RT)
    B = bM * p / RT
    roots = [z for z in cubic_real_roots(-(1.0 - B), A - 3.0 * B * B - 2.0 * B,
                                         -(A * B - B * B - B ** 3))
             if z > B * (1.0 + 1e-12)]
    if not roots:
        raise SystemExit(f"independent PR: no physical root at p={p} T={T}")

    def lnphi(Z):
        logTerm = math.log((Z + (1.0 + SQRT2) * B) / (Z + (1.0 - SQRT2) * B))
        out = []
        for i in range(n):
            sumA = sum(zs[j] * aIJ[i][j] for j in range(n))
            biB = bI[i] / bM
            out.append(biB * (Z - 1.0) - math.log(Z - B)
                       - A / (2.0 * SQRT2 * B) * (2.0 * sumA / aM - biB) * logTerm)
        return out

    distinct = sorted(set(roots))
    if pref == "liquid":
        Z = distinct[0]
    elif pref == "vapour":
        Z = distinct[-1]
    else:
        Z = min(distinct, key=lambda z: sum(zs[i] * lnphi(z)[i] for i in range(n)))
    return lnphi(Z)


def wilson_k(name, T, p):
    _, Tc, Pc, w = CT[name]
    return Pc / p * math.exp(5.373 * (1.0 + w) * (1.0 - Tc / T))


def flash_gas_fugacity(names, z, T, p, tol=1.0e-12, max_iter=500):
    """Independent flash (successive substitution, same verdict rule as
    src/pvt/src/flash.cpp). Returns (yOrFeed, lnPhi) of the gas side."""
    n = len(names)
    K = [wilson_k(names[i], T, p) for i in range(n)]

    def rr(beta):
        return sum(z[i] * (K[i] - 1.0) / (1.0 + beta * (K[i] - 1.0))
                   for i in range(n))

    def rr_solve():
        lo, hi, beta = 0.0, 1.0, 0.5
        for _ in range(100):
            f = rr(beta)
            if f > 0.0:
                lo = beta
            else:
                hi = beta
            df = -sum(z[i] * (K[i] - 1.0) ** 2 / (1.0 + beta * (K[i] - 1.0)) ** 2
                      for i in range(n))
            nxt = beta - f / df if df != 0.0 else 0.5 * (lo + hi)
            if not (lo < nxt < hi):
                nxt = 0.5 * (lo + hi)
            if abs(nxt - beta) < 1e-15:
                beta = nxt
                break
            beta = nxt
        return beta

    for _ in range(max_iter):
        f0, f1 = rr(0.0), rr(1.0)
        if f0 <= 0.0:
            beta = 0.0
            y = [K[i] * z[i] for i in range(n)]
            s = sum(y)
            x, y = z[:], [v / s for v in y]
        elif f1 >= 0.0:
            beta = 1.0
            y = z[:]
            s = sum(z[i] / K[i] for i in range(n))
            x = [z[i] / K[i] / s for i in range(n)]
        else:
            beta = rr_solve()
            x = [z[i] / (1.0 + beta * (K[i] - 1.0)) for i in range(n)]
            y = [K[i] * x[i] for i in range(n)]
            sx, sy = sum(x), sum(y)
            x, y = [v / sx for v in x], [v / sy for v in y]
        lnL = pr_lnphi(names, x, T, p, "liquid")
        lnV = pr_lnphi(names, y, T, p, "vapour")
        max_change = 0.0
        for i in range(n):
            new_k = math.exp(lnL[i] - lnV[i])
            max_change = max(max_change, abs(math.log(new_k / K[i])))
            K[i] = new_k
        if max_change < tol:
            break

    f0, f1 = rr(0.0), rr(1.0)
    spread = max(abs(k - 1.0) for k in K)
    two_phase = False
    if f0 > 0.0 and f1 < 0.0 and spread > 1e-6:
        beta = rr_solve()
        if 0.0 < beta < 1.0:
            x = [z[i] / (1.0 + beta * (K[i] - 1.0)) for i in range(n)]
            y = [K[i] * x[i] for i in range(n)]
            sx, sy = sum(x), sum(y)
            x, y = [v / sx for v in x], [v / sy for v in y]
            if max(abs(y[i] - x[i]) for i in range(n)) > 1e-8:
                two_phase = True
    if not two_phase:
        return z[:], pr_lnphi(names, z, T, p, "lowest")
    return y, pr_lnphi(names, y, T, p, "vapour")


def cages_per_water(structure, cage):
    ok("structure known", structure in ("sI", "sII"), structure)
    if structure == "sI":
        return 2.0 / 46.0 if cage == 0 else 6.0 / 46.0
    return 16.0 / 136.0 if cage == 0 else 8.0 / 136.0


def theta_at(names, z, T, p, langmuir_pa):
    y, lnV = flash_gas_fugacity(names, z, T, p)
    out = []
    for cage in (0, 1):
        sum_cf = sum(langmuir_pa[g][cage] * math.exp(lnV[g]) * y[g] * p
                     for g in range(len(names)))
        out.append(sum_cf / (1.0 + sum_cf))
    return out


def solve_p_eq(T, langmuir_pa, mu_ref, mu_slope, t_ref, nu, seed, feed_names,
               feed_z):
    """Independent bisection, same bracket rule as the kernel."""
    def g(p):
        return residual_at(T, p, langmuir_pa, mu_ref, mu_slope, t_ref, nu,
                           feed_names, feed_z)[0]

    lo, hi = P_FLOOR, seed
    if g(lo) <= 0.0:
        raise SystemExit(f"independent solve: hydrate stable at p->0, T={T}")
    while g(hi) > 0.0:
        hi *= 2.0
        if hi > P_CEILING:
            raise SystemExit(f"independent solve: no root below ceiling, T={T}")
    for _ in range(60):
        if hi - lo <= 1e-12 * hi:
            break
        mid = 0.5 * (lo + hi)
        if g(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    p = 0.5 * (lo + hi)
    return p, residual_at(T, p, langmuir_pa, mu_ref, mu_slope, t_ref, nu,
                          feed_names, feed_z)[1]


def residual_at(T, p, langmuir_pa, mu_ref, mu_slope, t_ref, nu, feed_names,
                feed_z):
    y, lnV = flash_gas_fugacity(feed_names, feed_z, T, p)
    th = []
    for cage in (0, 1):
        sum_cf = sum(langmuir_pa[g][cage] * math.exp(lnV[g]) * y[g] * p
                     for g in range(len(feed_names)))
        th.append(sum_cf / (1.0 + sum_cf))
    dmu = mu_ref + mu_slope * (T - t_ref)
    res = dmu + R * T * (nu[0] * math.log(1.0 - th[0])
                         + nu[1] * math.log(1.0 - th[1]))
    return res, th


# ---- parse the case ----
with open(case) as fh:
    source = fh.read()

comp = re.search(r"composition\s*\{(.*?)\}", source, re.S)
ok("composition block present", comp is not None)
feed_pairs = re.findall(r"([A-Za-z_]\w*)\s+([-+0-9.eE]+)\s*;", comp.group(1))
ok("composition non-empty", len(feed_pairs) > 0)
tot = sum(float(v) for _, v in feed_pairs)
FEED_NAMES = [k for k, _ in feed_pairs]
FEED_Z = [float(v) / tot for _, v in feed_pairs]

hp = re.search(r"hydratePWHc\s*\{(.*?)\n\}", source, re.S)
ok("hydratePWHc block present", hp is not None)
hp = hp.group(1)

structure = re.search(r"structure\s+(\w+)\s*;", hp).group(1)
guests = re.search(r"guests\s*\(([^)]*)\)", hp).group(1).split()
ok("guests non-empty", len(guests) > 0)
ok("guests declared in composition",
   all(g in FEED_NAMES for g in guests),
   f"guests {guests} vs feed {FEED_NAMES}")

lang_mpa = []
for g in guests:
    m = re.search(r"langmuir" + g + r"\s*\(([^)]*)\)", hp)
    ok(f"langmuir{g} present", m is not None)
    vals = [float(v) for v in m.group(1).split()]
    ok(f"langmuir{g} has two cage values", len(vals) == 2, str(vals))
    ok(f"langmuir{g} positive", all(v > 0.0 for v in vals))
    lang_mpa.append(vals)
# 1/MPa -> 1/Pa, exactly the conversion the app performs.
lang_pa = [[v * 1.0e-6 for v in row] for row in lang_mpa]

mu_ref = float(re.search(r"muRef\s+([-+0-9.eE]+)\s*;", hp).group(1))
mu_slope = float(re.search(r"muSlopeA\s+([-+0-9.eE]+)\s*;", hp).group(1))
t_ref = float(re.search(r"referenceTemperature\s+([-+0-9.eE]+)\s*;", hp).group(1))
p_ref = float(re.search(r"referencePressure\s+([-+0-9.eE]+)\s*;", hp).group(1))
d_h = float(re.search(r"enthalpyOfDissociation\s+([-+0-9.eE]+)\s*;", hp).group(1))
temps = [float(v) for v in
         re.search(r"temperatures\s*\(([^)]*)\)", hp).group(1).split()]
ok("temperatures non-empty", len(temps) > 0)
ok("temperatures positive", all(t > 0.0 for t in temps))

nu = [cages_per_water(structure, 0), cages_per_water(structure, 1)]

# ---- run the app on the declared case ----
proc = subprocess.run([exe, case], capture_output=True, text=True)
ok("app run succeeds", proc.returncode == 0, proc.stderr[-400:])
out = proc.stdout

h_lines = [ln for ln in out.splitlines() if "[chflow] pvt hydratePWHc" in ln]
hdr = [ln for ln in h_lines if "structure" in ln]
rows = [ln for ln in h_lines if " pEq " in ln]
mono_line = [ln for ln in h_lines if "monotonic" in ln]

ok("header line present", len(hdr) == 1, str(len(hdr)))
parts = hdr[0].replace("[chflow]", "").split()
ok("header structure echoes declaration", parts[3] == structure, parts[3])
got_cage_s = float(parts[parts.index("cagesSmall") + 1])
got_cage_l = float(parts[parts.index("cagesLarge") + 1])
close("cagesSmall matches unit-cell ratio", got_cage_s, nu[0], 1e-11)
close("cagesLarge matches unit-cell ratio", got_cage_l, nu[1], 1e-11)

ok("one equilibrium row per declared temperature", len(rows) == len(temps),
   f"{len(rows)} vs {len(temps)}")
ok("monotonic flag line present", len(mono_line) == 1)

got_T, got_p, got_ts, got_tl = [], [], [], []
for ln in rows:
    t = ln.replace("[chflow]", "").split()
    # pvt hydratePWHc T <T> pEq <p> thetaSmall <ts> thetaLarge <tl>
    got_T.append(float(t[3]))
    got_p.append(float(t[5]))
    got_ts.append(float(t[7]))
    got_tl.append(float(t[9]))

ok("temperatures echoed in case order",
   [round(t, 9) for t in got_T] == [round(t, 9) for t in temps])

# ---- independent re-derivation of every reported value ----
for i, T in enumerate(got_T):
    p_want, th_want = solve_p_eq(T, lang_pa, mu_ref, mu_slope, t_ref, nu,
                                 p_ref, FEED_NAMES, FEED_Z)
    close(f"T={T} p_eq independent", got_p[i], p_want, 1e-8)
    close(f"T={T} thetaSmall independent", got_ts[i], th_want[0], 1e-8)
    close(f"T={T} thetaLarge independent", got_tl[i], th_want[1], 1e-8)
    ok(f"T={T} thetaSmall in (0,1)", 0.0 < got_ts[i] < 1.0)
    ok(f"T={T} thetaLarge in (0,1)", 0.0 < got_tl[i] < 1.0)

# ---- physics: monotone curve in T, monotone theta in p ----
ok("p_eq strictly increases with T",
   all(got_p[i] < got_p[i + 1] for i in range(len(got_p) - 1)))
flag = int(mono_line[0].split()[-1])
ok("monotonic flag reported as 1", flag == 1, mono_line[0])

# Occupancy sweep at the anchor: theta must rise strictly with p, and the
# swept values must match the independent closed model.
p_anchor, th_anchor = solve_p_eq(t_ref, lang_pa, mu_ref, mu_slope, t_ref, nu,
                                 p_ref, FEED_NAMES, FEED_Z)
sweep = [0.5, 0.65, 0.8, 0.9, 1.0, 1.1, 1.25, 1.4, 1.6]
prev_s = prev_l = 0.0
for factor in sweep:
    p = factor * p_anchor
    th = theta_at(FEED_NAMES, FEED_Z, t_ref, p, lang_pa)
    ok(f"sweep p={p:.6g} thetaSmall in (0,1)", 0.0 < th[0] < 1.0)
    ok(f"sweep p={p:.6g} thetaLarge in (0,1)", 0.0 < th[1] < 1.0)
    ok(f"sweep p={p:.6g} thetaSmall rises with p", th[0] > prev_s)
    ok(f"sweep p={p:.6g} thetaLarge rises with p", th[1] > prev_l)
    prev_s, prev_l = th[0], th[1]

# ---- Clausius-Clapeyron physical bound on the slope at the anchor ----
p_lo, _ = solve_p_eq(t_ref - 0.5, lang_pa, mu_ref, mu_slope, t_ref, nu,
                     p_ref, FEED_NAMES, FEED_Z)
p_hi, _ = solve_p_eq(t_ref + 0.5, lang_pa, mu_ref, mu_slope, t_ref, nu,
                     p_ref, FEED_NAMES, FEED_Z)
slope = (math.log(p_hi) - math.log(p_lo)) / 1.0
target = d_h / (R * t_ref * t_ref)
ok("CC slope within 20% of dH/(R Tref^2)",
   abs(slope - target) <= 0.20 * target,
   f"slope {slope:.6g} target {target:.6g} ratio {slope / target:.4f}")


def run_case_text(text):
    fd, path = tempfile.mkstemp(suffix=".dict", prefix="pvt_pwhc_")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        return subprocess.run([exe, path], capture_output=True, text=True)
    finally:
        os.unlink(path)


def shrink_langmuir(text, factor):
    def repl(m):
        return m.group(1) + " ( " + " ".join(
            repr(float(v) * factor) for v in m.group(2).split()) + " );"
    return re.sub(r"(langmuir\w+)\s*\(([^)]*)\)\s*;", repl, text)


# ---- Langmuir -> 0 twin: theta -> 0 in closed form, p_eq diverges loudly ----
twin_lang = shrink_langmuir(source, 1.0e-6)
ok("twin case text rewrote the Langmuir lines", twin_lang != source)
lang_twin_pa = [[v * 1.0e-6 * 1.0e-6 for v in row] for row in lang_mpa]
th_twin = theta_at(FEED_NAMES, FEED_Z, t_ref, p_ref, lang_twin_pa)
y, lnV = flash_gas_fugacity(FEED_NAMES, FEED_Z, t_ref, p_ref)
for cage, name in ((0, "thetaSmall"), (1, "thetaLarge")):
    # Closed form with the flashed fugacity at the anchor point:
    sum_cf = sum(lang_twin_pa[g][cage] * math.exp(lnV[g]) * y[g] * p_ref
                 for g in range(len(FEED_NAMES)))
    theta_closed = sum_cf / (1.0 + sum_cf)
    close(f"Langmuir->0 twin {name} closed form", th_twin[cage], theta_closed,
          1e-8)
    ok(f"Langmuir->0 twin {name} theta -> 0", th_twin[cage] < 1.0e-5,
       f"{th_twin[cage]:.3g}")
r = run_case_text(twin_lang)
ok("Langmuir->0 twin runs", r.returncode == 0, r.stderr[-400:])
twin_rows = [ln for ln in r.stdout.splitlines()
             if "[chflow] pvt hydratePWHc" in ln and " pEq " in ln]
ok("Langmuir->0 twin has one row per temperature",
   len(twin_rows) == len(temps), str(len(twin_rows)))
for i, ln in enumerate(twin_rows):
    t = ln.replace("[chflow]", "").split()
    T_twin, p_twin = float(t[3]), float(t[5])
    p_want, _ = solve_p_eq(T_twin, lang_twin_pa, mu_ref, mu_slope, t_ref, nu,
                           p_ref, FEED_NAMES, FEED_Z)
    close(f"Langmuir->0 twin T={T_twin} p_eq independent",
          p_twin, p_want, 1e-8)
    ok(f"Langmuir->0 twin T={T_twin} p_eq far right of base",
       p_twin > 100.0 * got_p[i],
       f"p_eq {p_twin:.6g} vs base {got_p[i]:.6g}")
twin_p = [float(ln.replace("[chflow]", "").split()[5])
          for ln in twin_rows]
ok("Langmuir->0 twin p_eq strictly rises with T",
   all(twin_p[i] < twin_p[i + 1] for i in range(len(twin_p) - 1)))

# ---- single-guest degenerate twin: pure methane, closed-form occupancy ----
guest_case = f"""// E8 single-guest twin: pure methane, closed-form check.
temperature 280.0;
pressures (3.0e6);

composition
{{
    methane 1.0;
}}

solve {{ tol 1.0e-12; maxIter 500; }}

hydratePWHc
{{
    guests ( methane );
    structure {structure};
    langmuirmethane ( {lang_mpa[0][0]} {lang_mpa[0][1]} );
    muRef {mu_ref};
    muSlopeA {mu_slope};
    referenceTemperature {t_ref};
    referencePressure {p_ref};
    enthalpyOfDissociation {d_h};
    temperatures ( {' '.join(str(t) for t in temps)} );
}}
"""
r = run_case_text(guest_case)
ok("single-guest twin runs", r.returncode == 0, r.stderr[-400:])
rows_single = [ln for ln in r.stdout.splitlines()
               if "[chflow] pvt hydratePWHc" in ln and " pEq " in ln]
ok("single-guest twin has one row per temperature",
   len(rows_single) == len(temps), str(len(rows_single)))
for i, ln in enumerate(rows_single):
    t = ln.replace("[chflow]", "").split()
    T, p_rep = float(t[3]), float(t[5])
    ts_rep, tl_rep = float(t[7]), float(t[9])
    # Closed form: single gas, f = phi * p, theta = C f / (1 + C f).
    lnV = pr_lnphi(["methane"], [1.0], T, p_rep, "lowest")
    f_gas = math.exp(lnV[0]) * p_rep
    ts_want = lang_pa[0][0] * f_gas / (1.0 + lang_pa[0][0] * f_gas)
    tl_want = lang_pa[0][1] * f_gas / (1.0 + lang_pa[0][1] * f_gas)
    close(f"single-guest T={T} thetaSmall closed form", ts_rep, ts_want, 1e-8)
    close(f"single-guest T={T} thetaLarge closed form", tl_rep, tl_want, 1e-8)
    # p_eq from the independent single-guest solver.
    p_want, _ = solve_p_eq(T, lang_pa, mu_ref, mu_slope, t_ref, nu, p_ref,
                           ["methane"], [1.0])
    close(f"single-guest T={T} p_eq independent", p_rep, p_want, 1e-8)

# ---- loud rejections ----
rejects = [
    ("guest not declared in composition",
     source.replace("guests ( methane propane );",
                    "guests ( methane nitrogen );")),
    ("langmuir cage-count mismatch",
     source.replace("langmuirmethane ( 2.0 1.0 );",
                    "langmuirmethane ( 2.0 );")),
    ("missing langmuir key",
     source.replace("langmuirpropane ( 0.01 1.8 );  // 1/MPa", "")),
    ("structure typo",
     source.replace("structure sII;", "structure sIII;")),
    ("negative langmuir constant",
     source.replace("langmuirmethane ( 2.0 1.0 );",
                    "langmuirmethane ( -2.0 1.0 );")),
    ("empty temperature list",
     source.replace("temperatures ( 274.15 277.15 280.15 283.15 );",
                    "temperatures ( );")),
    ("mutual exclusion with hydrateEquilibrium",
     source.replace("hydratePWHc\n{",
                    "hydrateEquilibrium\n{\n"
                    "    Tref 285.25;\n    pRef 3.0e6;\n    dH 5.1868e4;\n"
                    "}\n\nhydratePWHc\n{")),
    ("negative muRef",
     source.replace("muRef 600.0;", "muRef -600.0;")),
    ("zero referencePressure",
     source.replace("referencePressure 2.6e6;", "referencePressure 0;")),
]
for label, bad in rejects:
    r = run_case_text(bad)
    ok(f"rejection: {label}", r.returncode != 0,
       f"rc={r.returncode} stderr={r.stderr[-200:]}")

# ---- zero-declaration twin: flash report byte-identical ----
twin_zero = re.sub(r"hydratePWHc\s*\{.*?\n\}", "", source, flags=re.S)
r = run_case_text(twin_zero)
ok("zero-declaration twin runs", r.returncode == 0, r.stderr[-400:])
twin_flash = [ln for ln in r.stdout.splitlines() if "[chflow] pvt flash" in ln]
base_flash = [ln for ln in out.splitlines() if "[chflow] pvt flash" in ln]
ok("flash report byte-identical without hydratePWHc",
   twin_flash == base_flash)
ok("no hydratePWHc tokens in zero-declaration twin",
   "hydratePWHc" not in r.stdout)
ok("declared hydratePWHc is observable (changes the report)",
   "[chflow] pvt hydratePWHc" in out)

print(f"PASS: E8 hydratePWHc ({checks} checks): "
      f"{len(temps)} equilibrium points re-derived to 1e-8, "
      f"p_eq(T={got_T[0]})={got_p[0]:.6g} Pa, "
      f"thetaSmall {got_ts[0]:.6f}->{got_ts[-1]:.6f}, "
      f"CC slope ratio {slope / target:.3f}, "
      f"{len(rejects)} rejections, "
      f"Langmuir->0 and single-guest twins green, "
      f"zero-declaration twin byte-identical")
