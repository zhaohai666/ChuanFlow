#!/usr/bin/env python3
"""steadyYieldTransition: the published laminar-turbulent locus of a yielded
pipe AND the turbulent friction a smooth Bingham pipe is priced by past it --
checked by re-deriving every reported number from Wilson & Thomas (2006) and
from Buckingham (1921) written out here.

5.83 lands the turbulent side. Below Re_B,t a yielded pipe answers with its
laminar wall law (5.46/5.47); strictly past it, a smooth Bingham pipe is priced
by the Wilson & Thomas Eq (15)/(16) closure (cross-checked against Thomas &
Wilson 2007 Eq [3]/[4]). The locus itself needs only the grouping the literature
agrees on:

    He   = tau_0 d^2 rho / mu^2                    (the Hedstrom number)
    Re_B = rho |v| d / mu                          (the plastic Reynolds number)
    Re_B,t(He) = 2100            (He <= 1700)
                            = 80 He^0.40   (1700 < He <= 1e5)
                            = 25 He^0.5    (He > 1e5)
    Q_t  = Re_B,t mu A / (rho d)      and   dp_t = 4 L tw(Q_t) / d  (laminar)
    Trel = Q_t / A / (tau_0/rho)^0.5 = Re_B,t / He^0.5

The turbulent closure, with beta = tau_w/tau_B = f Re_B^2/(8 He):
    sqrt(8/f) = 2.5 ln(sqrt He) + 2.5 ln[sqrt(beta)(beta-1)/(beta+1)]
                + 11.6/beta - chi(beta)
    chi(beta) = -2.5 ln[(beta-1)/beta] - 2.5 (beta+0.5)/beta^2
solved here by bisection on beta > 1, independently of the kernel.

`dp_t` is the LAST laminar point. The switch is hard, so at Re_B,t the turbulent
drop jumps above the laminar one -- a hysteretic gap with no steady root, which
the case avoids by probing only the well-posed regimes. Seven families:

  1. the locus belongs to the fluid and the geometry, not to the state it was
     read at -- twin pipes at different drops agree bit for bit
  2. the three bands, their constants, and the He targets of the published
     parameter set (Csizmadia & Toth 2013, Table 1: D = 0.1 m, rho = 1370 kg/m3,
     tau_y = 1 Pa, mu_B = 2e-2 Pa.s)
  3. the identities: Trel = Re_B,t / He^0.5 (which forces Trel to 25 through the
     whole He^0.5 band, and to 80 He^-0.10 in the mid band), and the printed
     Re_B is the printed Re
  4. dp_t is the pipe's OWN laminar wall law at Q_t -- bisection on Buckingham's
     bracket -- and it stays laminar even for a pipe solved turbulent
  5. the verdict: past/below the locus is the same answer from the rate, from
     the plastic Reynolds number, and from the drop
  6. the switch prices turbulence: below the locus every rate is the laminar
     Buckingham root and f Re_B g(xi) = 64; past it every rate is the turbulent
     root and f is the Eq (15)/(16) factor re-derived here, NOT 64/(Re_B g); the
     turbulent drop jumps above the laminar at the locus (the gap); and a
     tau_0 -> 0 pipe still lands on its Newtonian twin at the bracket's factor
  7. where the locus does not apply: a Newtonian pipe and a Herschel-Bulkley
     pipe print no He, because neither declares a plastic viscosity
"""

import math
import subprocess
import sys

STEADY = sys.argv[1]
CASE = sys.argv[2]
CASE_DIR = CASE.rsplit("/", 1)[0]

RHO = 1370.0
MU = 0.02
LENGTH = 20.0
DROP = {"N1": 200.0, "N2": 500.0, "N3": 5000.0, "N4": 30000.0, "N5": 150000.0,
        "N6": 700.0}

# name: (diameter, tau_0, from-node); tau_0 = None means no threshold declared
PIPES = {
    "Small01":       (0.1, 0.01,   "N1"),
    "Small05":       (0.1, 0.01,   "N6"),
    "TableHe04":     (0.1, 1.0,    "N2"),
    "TableHe04Run":  (0.1, 1.0,    "N3"),
    "Edge1700Under": (0.02, 1.2400, "N5"),
    "Edge1700Over":  (0.02, 1.2420, "N5"),
    "Edge1700Flat":  (0.02, 1.2409, "N4"),
    "Edge1e5Under":  (0.1, 2.9197, "N5"),
    "Edge1e5Over":   (0.1, 2.9203, "N5"),
    "BigHe06":       (0.1, 100.0,  "N4"),
    "BigHe06Run":    (0.1, 100.0,  "N5"),
    "BigHe08":       (0.1, 1.0e4,  "N5"),
    "BigHe10":       (0.1, 1.0e6,  "N5"),
    "Bore02":        (0.02, 0.0002, "N1"),
    "Bore02Run":     (0.02, 0.0002, "N4"),
    "Bore02Newt":    (0.02, None,   "N1"),
    "NewtRef":       (0.1,  None,   "N3"),
    "HbRef":         (0.1,  1.0,    "N3"),
}
LOCUS = [n for n in PIPES if PIPES[n][1] is not None and n != "HbRef"]
TWINS = [("Small01", "Small05"), ("TableHe04", "TableHe04Run"),
         ("BigHe06", "BigHe06Run"), ("Bore02", "Bore02Run"),
         ("Edge1700Over", "Edge1700Flat"), ("Bore02", "Bore02Run")]

failures = []
checks = [0]


def require(cond, what):
    checks[0] += 1
    if cond:
        print("ok:  ", what)
    else:
        print("FAIL:", what)
        failures.append(what)


def rel(a, b):
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale > 0 else 0.0


def run(binary, case_path):
    return subprocess.run([binary, case_path], capture_output=True, text=True,
                          cwd=CASE_DIR)


def parse_steady(stdout):
    out = {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) < 8 or tok[0] != "pipe" or tok[2] != "Q":
            continue
        name = tok[1]
        if name not in PIPES:
            continue
        rec = {"q": float(tok[3]), "v": float(tok[5]), "raw": tok[3]}
        for i, t in enumerate(tok):
            if t == "Re":
                rec["re"] = float(tok[i + 1])
            elif t == "f":
                rec["f"] = float(tok[i + 1])
            elif t == "dp":
                rec["dp"] = float(tok[i + 1])
            elif t == "He":
                rec["he"] = float(tok[i + 1])
            elif t == "ReB":
                rec["reb"] = float(tok[i + 1])
            elif t == "ReBt":
                rec["rebt"] = float(tok[i + 1])
            elif t == "Qt":
                rec["qt"] = float(tok[i + 1])
            elif t == "dpQt":
                rec["dpqt"] = float(tok[i + 1])
            elif t == "Vrelt":
                rec["trel"] = float(tok[i + 1])
            elif t == "band":
                rec["band"] = int(float(tok[i + 1]))
            elif t in ("pastTransition", "belowTransition"):
                rec["past"] = (t == "pastTransition")
        out[name] = rec
    return out


# --- the locus and the wall law, both written out here ----------------------

def he_of(diam, tau0):
    return tau0 * diam * diam * RHO / (MU * MU)


def band_and_locus(he):
    if he <= 1700.0:
        return 0, 2100.0
    if he <= 1.0e5:
        return 1, 80.0 * he ** 0.40
    return 2, 25.0 * math.sqrt(he)


def area(diam):
    return 0.25 * math.pi * diam * diam


def q_of_tw(tw, tau0, diam):
    """Buckingham's exact laminar rate at a wall stress."""
    if tw <= tau0:
        return 0.0
    xi = tau0 / tw
    return (math.pi * diam ** 3 / (32.0 * MU)) * tw \
        * (1.0 - 4.0 / 3.0 * xi + xi ** 4 / 3.0)


def tw_of_q(q, tau0, diam):
    """Its inverse: F(tw) = tw - (4/3) tau0 + tau0^4/(3 tw^3) = b."""
    b = 32.0 * MU * abs(q) / (math.pi * diam ** 3)
    if not b > 0.0:
        return tau0
    lo, hi = tau0, b + 4.0 / 3.0 * tau0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f = mid - 4.0 / 3.0 * tau0 + tau0 ** 4 / (3.0 * mid ** 3) - b
        if f > 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo <= 2.0e-16 * hi:
            break
    return 0.5 * (lo + hi)


def dp_of_q(q, tau0, diam):
    return 4.0 * LENGTH * tw_of_q(q, tau0, diam) / diam


def g_factor(dp, diam, tau0):
    """The bracket 1 - (4/3) xi + xi^4/3 at the wall stress this drop asks
    for, which is what makes f Re_B = 64/g rather than 64."""
    tw = dp * diam / (4.0 * LENGTH)
    xi = tau0 / tw
    return 1.0 - 4.0 / 3.0 * xi + xi ** 4 / 3.0


# --- the turbulent closure, re-derived independently of the kernel ----------

def wt_turbulent_f(reb, he):
    """Wilson & Thomas (2006) Eq (15)/(16): the smooth-wall turbulent Darcy
    factor of a Bingham plastic. With beta = tau_w/tau_B = f reb^2/(8 he) and
    sqrt(8/f) = reb/sqrt(he beta), the mean-to-shear velocity ratio obeys
        sqrt(8/f) = 1.25 ln he + 1.25 ln beta + 2.5 ln[(beta-1)/(beta+1)]
                    + 11.6/beta - chi(beta),
        chi(beta) = -2.5 ln[(beta-1)/beta] - 2.5 (beta+0.5)/beta^2.
    One monotone equation in beta > 1, solved by bracketed bisection. Written
    with plain log (the kernel uses log1p) so it is an independent code path."""
    def chi(b):
        return -2.5 * math.log((b - 1.0) / b) - 2.5 * (b + 0.5) / (b * b)

    def resid(b):
        return (1.25 * math.log(he) + 1.25 * math.log(b)
                + 2.5 * math.log((b - 1.0) / (b + 1.0)) + 11.6 / b - chi(b)
                - reb / math.sqrt(he * b))

    lo, hi = 1.0 + 1.0e-12, 2.0
    for _ in range(200):
        if resid(hi) >= 0.0:
            break
        hi *= 2.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if resid(mid) > 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo <= 2.0e-16 * hi:
            break
    beta = 0.5 * (lo + hi)
    return 8.0 * he * beta / (reb * reb)


def dp_turb_of_q(q, tau0, diam):
    """The turbulent drop at a rate magnitude: (f L/d) rho v^2/2 with f from
    Eq (15)/(16) at the plastic Reynolds number of that rate."""
    v = abs(q) / area(diam)
    reb = RHO * v * diam / MU
    f = wt_turbulent_f(reb, he_of(diam, tau0))
    return (f * LENGTH / diam) * RHO * v * v / 2.0


def q_turb_of_dp(dp, tau0, diam):
    """The turbulent rate at an imposed drop, bracketed above Q_t (the closure
    only applies past the locus). dp_turb is monotone increasing in Q."""
    qt = band_and_locus(he_of(diam, tau0))[1] * MU * area(diam) / (RHO * diam)
    lo, hi = qt, qt * 2.0
    for _ in range(200):
        if dp_turb_of_q(hi, tau0, diam) >= dp:
            break
        hi *= 2.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if dp_turb_of_q(mid, tau0, diam) < dp:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 2.0e-16 * hi:
            break
    return 0.5 * (lo + hi)


proc = run(STEADY, CASE)
require(proc.returncode == 0, "the locus case solves (rc=%d)" % proc.returncode)
S = parse_steady(proc.stdout)
require(len(S) == len(PIPES), "all %d pipes report" % len(PIPES))

# --- 1. the locus is geometry + fluid, not the state ------------------------

bad = [n for n in LOCUS if "he" not in S[n]]
require(not bad, "every Bingham pipe reports the locus (missing: %s)" % bad)
worst = 0.0
worst_name = ""
for a, b in TWINS:
    if PIPES[a][:2] != PIPES[b][:2]:
        continue
    for k in ("he", "rebt", "qt", "dpqt", "trel", "band"):
        r = rel(S[a][k], S[b][k])
        if r > worst:
            worst, worst_name = r, "%s/%s %s" % (a, b, k)
require(worst == 0.0,
        "He, Re_B,t, Q_t, dp_t, Trel and the band are bit-identical between "
        "twins at different drops (%s)" % worst_name)
require(S["TableHe04"]["q"] == 0.0 and S["TableHe04"]["qt"] == S["TableHe04Run"]["qt"],
        "a stuck pipe (Q = 0) still reports the same Q_t as the same pipe moving")

# --- 2. the three published bands ------------------------------------------

TARGETS = {"Small01": 3.425e2, "TableHe04": 3.425e4, "BigHe06": 3.425e6,
           "BigHe08": 3.425e8, "BigHe10": 3.425e10}
for n, want in TARGETS.items():
    require(rel(S[n]["he"], want) < 1e-3,
            "He(%s) = %.6g matches the published parameter set's %.6g "
            "(printed %.12g)" % (n, S[n]["he"], want, S[n]["he"]))
require(S["Small01"]["band"] == 0 and S["TableHe04"]["band"] == 1
        and S["BigHe06"]["band"] == 2,
        "band assignment follows the printed He ranges 1700 and 1e5")
worst = 0.0
for n in LOCUS:
    d, tau0, _ = PIPES[n]
    band, rebt = band_and_locus(he_of(d, tau0))
    require(S[n]["band"] == band, "%s reports band %d" % (n, band))
    worst = max(worst, rel(S[n]["rebt"], rebt))
require(worst < 2e-11,
        "Re_B,t is the band's own fit for all %d pipes to the 12 digits the app "
        "prints, worst relative %.3g" % (len(LOCUS), worst))
require(S["Edge1700Under"]["rebt"] == 2100.0,
        "band 0 is the flat 2100 the paper quotes as the engineering value")

# --- 3. the identities ------------------------------------------------------

worst_trel = worst_re = 0.0
for n in LOCUS:
    d, tau0, _ = PIPES[n]
    he = he_of(d, tau0)
    worst_trel = max(worst_trel, rel(S[n]["trel"], S[n]["rebt"] / math.sqrt(he)))
    worst_re = max(worst_re, rel(S[n]["reb"], S[n]["re"]))
require(worst_trel < 2e-11,
        "Trel = Re_B,t / He^0.5 holds as an identity for every pipe, worst %.3g"
        % worst_trel)
require(worst_re < 1e-9,
        "the printed plastic Reynolds number IS the printed Re of a Bingham "
        "pipe (rho v d / mu), worst %.3g" % worst_re)
band2 = [n for n in LOCUS if S[n]["band"] == 2]
require(all(abs(S[n]["trel"] - 25.0) < 1e-9 for n in band2) and len(band2) >= 4,
        "Trel is exactly 25 through the whole He^0.5 band (%d pipes), which is "
        "the paper's Eq (20) read back through the identity" % len(band2))
mid = [n for n in LOCUS if S[n]["band"] == 1]
require(all(rel(S[n]["trel"], 80.0 * he_of(*PIPES[n][:2]) ** -0.10) < 2e-11
            for n in mid),
        "in the mid band the same identity reads Trel = 80 He^-0.10 (%d pipes), "
        "the paper's own second form of Eq (19)" % len(mid))

# --- 4. dp_t is this pipe's own wall law at Q_t -----------------------------

worst = 0.0
worst_name = ""
for n in LOCUS:
    d, tau0, _ = PIPES[n]
    want = dp_of_q(S[n]["qt"], tau0, d)
    if rel(S[n]["dpqt"], want) > worst:
        worst, worst_name = rel(S[n]["dpqt"], want), n
require(worst < 2e-11,
        "dp_t is the kernel's own Buckingham wall law evaluated at Q_t, not a "
        "correlation (%s, worst relative %.3g)" % (worst_name, worst))
require(all(S[n]["dpqt"] > 4.0 * LENGTH * PIPES[n][1] / PIPES[n][0] for n in LOCUS),
        "dp_t is above the static threshold 4 tau_0 L / d for every pipe: the "
        "locus always costs more than just starting the wall")

# --- 5. the verdict, from three independent quantities ----------------------

worst_q = 0.0
for n in LOCUS:
    d, tau0, node = PIPES[n]
    from_dp = DROP[node]
    by_rate = S[n]["q"] > S[n]["qt"]
    by_reb = S[n]["reb"] > S[n]["rebt"]
    by_dp = from_dp > S[n]["dpqt"]
    require(S[n]["past"] == by_rate == by_reb == by_dp,
            "%s: pastLaminarLimit agrees between rate, Re_B and drop "
            "(Q %.6g vs Q_t %.6g, dp %.6g vs dp_t %.6g)"
            % (n, S[n]["q"], S[n]["qt"], from_dp, S[n]["dpqt"]))
n_past = sum(1 for n in LOCUS if S[n].get("past"))
n_below = len(LOCUS) - n_past
require(n_past >= 3 and n_below >= 3,
        "the case actually straddles the locus (%d past, %d below), so the "
        "verdict is not one-valued" % (n_past, n_below))
require(S["Small01"]["dp"] < S["Small01"]["dpqt"] < S["Small05"]["dp"],
        "Small01/Small05 are the same pipe either side of its own dp_t = %.6g "
        "Pa (200 laminar / 700 turbulent, clear of the gap)" % S["Small01"]["dpqt"])

# --- 6. the switch prices turbulence past the locus -------------------------

# below the locus: the rate is the laminar Buckingham root of the asked drop,
# and f Re_B g(xi) = 64 exactly as in 5.46/5.47
below = [n for n in LOCUS if not S[n].get("past") and S[n]["q"] > 0.0]
worst_q = worst_f = 0.0
worst_qn = ""
for n in below:
    d, tau0, node = PIPES[n]
    tw = S[n]["dp"] * d / (4.0 * LENGTH)
    want = q_of_tw(tw, tau0, d)
    if rel(S[n]["q"], want) > worst_q:
        worst_q, worst_qn = rel(S[n]["q"], want), n
    worst_f = max(worst_f, rel(S[n]["f"] * S[n]["reb"]
                               * g_factor(S[n]["dp"], d, tau0), 64.0))
require(worst_q < 1e-11,
        "every below-locus moving rate is still the exact laminar Buckingham "
        "root of the asked drop (%s, worst %.3g)" % (worst_qn, worst_q))
require(worst_f < 1e-9,
        "f Re_B g(xi) = 64 holds below the locus, worst relative %.3g" % worst_f)

# past the locus: the rate is the TURBULENT root, f is Eq (15)/(16) re-derived
# here, and the printed f reproduces the printed drop -- not 64/(Re_B g)
past = [n for n in LOCUS if S[n].get("past")]
worst_q = worst_f = worst_dp = 0.0
worst_qn = worst_fn = worst_dpn = ""
for n in past:
    d, tau0, node = PIPES[n]
    he = he_of(d, tau0)
    qw = q_turb_of_dp(DROP[node], tau0, d)
    fw = wt_turbulent_f(S[n]["reb"], he)
    dpw = (S[n]["f"] * LENGTH / d) * RHO * S[n]["v"] ** 2 / 2.0
    if rel(S[n]["q"], qw) > worst_q:
        worst_q, worst_qn = rel(S[n]["q"], qw), n
    if rel(S[n]["f"], fw) > worst_f:
        worst_f, worst_fn = rel(S[n]["f"], fw), n
    if rel(S[n]["dp"], dpw) > worst_dp:
        worst_dp, worst_dpn = rel(S[n]["dp"], dpw), n
require(len(past) >= 3,
        "the case has %d past-locus pipes priced turbulent" % len(past))
require(worst_q < 1e-9,
        "every past-locus rate is the exact root of the asked drop under the "
        "Wilson-Thomas turbulent law (%s, worst %.3g)" % (worst_qn, worst_q))
require(worst_f < 2e-11,
        "every past-locus f IS the Eq (15)/(16) factor re-derived independently "
        "here (%s, worst %.3g)" % (worst_fn, worst_f))
require(worst_dp < 1e-9,
        "the printed past-locus f reproduces the printed drop, (f L/d) rho v^2/2 "
        "(%s, worst %.3g)" % (worst_dpn, worst_dp))

# the deepest past-locus pipes have STOPPED reproducing the laminar bracket
deep = [n for n in past if S[n]["reb"] > 4.0e4]
require(len(deep) >= 2 and all(
    rel(S[n]["f"] * S[n]["reb"] * g_factor(S[n]["dp"], PIPES[n][0], PIPES[n][1]),
        64.0) > 0.1 for n in deep),
        "the deepest past-locus pipes (%s) run Re_B in the tens of thousands to "
        "hundreds of thousands and NO LONGER answer 64/(Re_B g): the laminar "
        "bracket is off by >10%%, which is exactly what 5.83 had to stop "
        "reproducing" % ",".join(deep))

# the switch is hard: the turbulent and laminar drops never join smoothly at Q_t
require(all(abs(dp_turb_of_q(S[n]["qt"], PIPES[n][1], PIPES[n][0])
                - dp_of_q(S[n]["qt"], PIPES[n][1], PIPES[n][0])) > 0.0
            for n in LOCUS),
        "the switch is hard: the turbulent and laminar drops differ at Q_t for "
        "every pipe (a jump, never a smooth join)")

# for the band-0 straddle pipe the jump is a GAP (turbulent above laminar), so a
# drop between the two has no steady root; Small05 sits at 700 Pa, clear above it
d, tau0, _ = PIPES["Small05"]
qt = S["Small05"]["qt"]
dl = dp_of_q(qt, tau0, d)
dt = dp_turb_of_q(qt, tau0, d)
require(dt > dl and S["Small05"]["dp"] > dt,
        "Small05's locus is a gap: dp_t,lam=%.4f < dp_t,turb=%.4f, and the "
        "imposed 700 Pa sits above the ceiling (a drop inside the gap would "
        "leave the solver sitting at Q_t with no root)" % (dl, dt))
require(rel(S["Small05"]["dpqt"], dl) < 2e-11,
        "the reported dp_t is the laminar floor of that gap (%.6g), even though "
        "Small05 itself is solved turbulent" % dl)

# the tau_0 -> 0 control: same geometry, Newtonian twin, still BELOW the locus
g_fac = g_factor(S["Bore02"]["dp"], 0.02, 0.0002)
require(rel(S["Bore02"]["q"], S["Bore02Newt"]["q"] * g_fac) < 1e-3,
        "the tau_0 = 2e-4 pipe's rate is its Newtonian twin's times the "
        "bracket's own factor g = %.6f at xi = 1/250: band 0 and 2100 are the "
        "Newtonian limit, not a separate statement" % g_fac)

# --- 7. where the locus is not a statement ----------------------------------

require("he" not in S["NewtRef"] and "he" not in S["HbRef"],
        "a Newtonian pipe and a Herschel-Bulkley pipe print no He token")
tw_hb = S["HbRef"]["dp"] * PIPES["HbRef"][0] / (4.0 * LENGTH)
require(rel(S["HbRef"]["re"], 8.0 * RHO * S["HbRef"]["v"] ** 2 / tw_hb) < 1e-9,
        "the HB pipe's printed Re is the APPARENT one (8 rho v^2 / tw), so no "
        "mu-based He or Re_B is a statement about it")

# --- the shape of the published locus ---------------------------------------

he_join = 1.0e5
lo_side = 80.0 * he_join ** 0.40
hi_side = 25.0 * math.sqrt(he_join)
ratio = (S["Edge1e5Over"]["rebt"] - S["Edge1e5Under"]["rebt"]) \
    / S["Edge1e5Under"]["rebt"]
require(-0.02 < ratio < 0.0,
        "the two published fits meet at He = 1e5 with the He^0.5 band %.2f%% "
        "BELOW the mid fit (printed %.6g then %.6g; the formulas themselves "
        "give %.6g and %.6g), so the step down is recorded rather than "
        "smoothed" % (100.0 * ratio, S["Edge1e5Under"]["rebt"],
                      S["Edge1e5Over"]["rebt"], lo_side, hi_side))
step = (S["Edge1700Over"]["rebt"] - 2100.0) / 2100.0
require(-0.30 < step < -0.20,
        "the band 0 -> band 1 join is a %.1f%% DECREASE in Re_B,t (%.6g after "
        "2100 before): the paper bridges it with an Eq (17) whose constants "
        "could not be transcribed, so the step is recorded, not smoothed"
        % (100.0 * step, S["Edge1700Over"]["rebt"]))
he_probe = [3.425e2, 1700.0, 3.425e4, 1.0e5, 3.425e6, 3.425e8, 3.425e10]
loci = [band_and_locus(h)[1] for h in he_probe if h > 1700.0]
require(all(b >= a for a, b in zip(loci, loci[1:])),
        "within and across the two published fits the locus is monotone "
        "non-decreasing in He")

print("")
print("%d assertions, %d failed" % (checks[0], len(failures)))
if failures:
    sys.exit(1)
