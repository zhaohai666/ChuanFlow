#!/usr/bin/env python3
"""steadyYieldTransition: the published laminar-turbulent locus of a yielded
pipe, read off the pipe's own wall law -- checked by re-deriving every reported
number from Wilson & Thomas (2006) and from Buckingham (1921) written out here.

What landed is a report, not a re-pricing. The row still answers a yielded pipe
with its laminar wall law (5.46/5.47), because the published turbulent-yield
FRICTION closure could not be transcribed with its constants intact. The locus,
though, needs only the grouping the literature agrees on and the laminar law
that is already in the kernel:

    He   = tau_0 d^2 rho / mu^2                    (the Hedstrom number)
    Re_B = rho |v| d / mu                          (the plastic Reynolds number)
    Re_B,t(He) = 2100            (He <= 1700)
                            = 80 He^0.40   (1700 < He <= 1e5)
                            = 25 He^0.5    (He > 1e5)
    Q_t  = Re_B,t mu A / (rho d)      and   dp_t = 4 L tw(Q_t) / d
    Trel = Q_t / A / (tau_0/rho)^0.5 = Re_B,t / He^0.5

`dp_t` is the useful one: it is how far the injection end has to move before the
statement the kernel is making stops being the one that applies. Seven families:

  1. the locus belongs to the fluid and the geometry, not to the state it was
     read at -- twin pipes at different drops agree bit for bit
  2. the three bands, their constants, and the He targets of the published
     parameter set (Csizmadia & Toth 2013, Table 1: D = 0.1 m, rho = 1370 kg/m3,
     tau_y = 1 Pa, mu_B = 2e-2 Pa.s)
  3. the identities: Trel = Re_B,t / He^0.5 (which forces Trel to 25 through the
     whole He^0.5 band, and to 80 He^-0.10 in the mid band), and the printed
     Re_B is the printed Re
  4. dp_t is the pipe's OWN wall law at Q_t -- bisection on Buckingham's bracket
  5. the verdict: past/below the locus is the same answer from the rate, from
     the plastic Reynolds number, and from the drop
  6. the report is inert: every rate is still the exact root of the asked drop,
     f Re_B g(xi) = 64 still holds past the locus, and a tau_0 -> 0 pipe lands
     on its Newtonian twin at the bracket's own factor
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
DROP = {"N1": 200.0, "N2": 500.0, "N3": 5000.0, "N4": 30000.0, "N5": 150000.0}

# name: (diameter, tau_0, from-node); tau_0 = None means no threshold declared
PIPES = {
    "Small01":       (0.1, 0.01,   "N1"),
    "Small05":       (0.1, 0.01,   "N2"),
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
        "Pa (200 / 500 Pa)" % S["Small01"]["dpqt"])

# --- 6. the report is inert -------------------------------------------------

worst = 0.0
worst_name = ""
for n in LOCUS:
    d, tau0, node = PIPES[n]
    tw = S[n]["dp"] * d / (4.0 * LENGTH)
    want = q_of_tw(tw, tau0, d)
    if want > 0 and rel(S[n]["q"], want) > worst:
        worst, worst_name = rel(S[n]["q"], want), n
require(worst < 1e-11,
        "every moving rate is still the exact root of the asked drop under the "
        "laminar law -- reading the locus charged nothing (%s, worst %.3g)"
        % (worst_name, worst))
worst = 0.0
for n in LOCUS:
    if S[n]["q"] <= 0.0:
        continue
    d, tau0, node = PIPES[n]
    worst = max(worst, rel(S[n]["f"] * S[n]["reb"] * g_factor(S[n]["dp"], d, tau0),
                           64.0))
require(worst < 1e-9,
        "f Re_B g(xi) = 64 still holds past the locus, worst relative %.3g -- "
        "the locus did not reach into the friction factor" % worst)
deep = [n for n in LOCUS if S[n].get("past") and S[n]["reb"] > 4.0e4]
require(len(deep) >= 2 and all(
    rel(S[n]["f"] * S[n]["reb"] * g_factor(S[n]["dp"], PIPES[n][0], PIPES[n][1]),
        64.0) < 1e-9 for n in deep),
        "the deepest past-locus pipes (%s) run Re_B in the tens of thousands to "
        "hundreds of thousands and still answer 64 / (Re_B g), i.e. the laminar "
        "bracket -- which is exactly the number the next increment has to stop "
        "reproducing" % ",".join(deep))

# the tau_0 -> 0 control: same geometry, Newtonian twin
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
