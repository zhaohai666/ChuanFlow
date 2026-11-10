#!/usr/bin/env python3
"""steadyHerschelBulkley: tau = tau_0 + K gdot^n in the steady network, checked
against the integrated velocity profile rather than against a convention.

One formula is claimed to contain both shipped cells: with tw = dp d/(4L) and
D = tw - tau_0,

    Q = pi R^3 tw^-3 K^(-1/n) n D^(1+1/n)
        [ tau0^2/(n+1) + 2 tau0 D/(2n+1) + D^2/(3n+1) ]

at n = 1, K = mu this is Buckingham's [1 - (4/3)xi + (1/3)xi^4] (5.46), and at
tau0 = 0 only the last term survives and it is the exact power-law discharge
(5.45). The second limit is what the increment also fixes: 5.45 wrote the
Metzner-Reed consistency as K (3n+1)/(4n) -- exponent 1, the turbulent-collapse
convention -- so its laminar line sat ((3n+1)/(4n))^(1-1/n) off the profile
integral, 25% in the rate at n = 0.5. The kernel now uses K ((3n+1)/(4n))^n, so
the M-R path and a vanishing yield stress land on the same number.

Every yielded pipe with an index reports the *apparent* Reynolds number
Re = 8 rho v^2 / tw, because no material viscosity enters its law; with
f = 8 tw/(rho v^2) the product f Re is then 64 by definition, so the load of
this test is carried by the closed form and by the two limits, and the
identity is asserted as a consistency of the two printed numbers.

Nine families of judgement:
  1. the declaration read back (twelve pipes, two driving differences)
  2. every yielded rate against the termwise closed form, then 2b an account of
     the two states that land soft enough for the solver's own stop to show
  3. the same rates against an independent quadrature of the profile integral
  4. limit n = 1: identical to the Bingham branch of 5.46, Re apart by g(xi)
  5. limit tau0 -> 0: identical to the corrected Metzner-Reed branch, and the
     size of the constant 5.45 had wrong
  6. the threshold is geometric and the dead branch is asked, not owned
  7. f Re == 64 and dp == f (L/d) rho v^2 / 2 on every yielded state
  8. the series pair: one free node splitting the drop between two yielded laws
  9. the declaration-face refusal that survives the lifting
"""

import math
import subprocess
import sys

BINARY = sys.argv[1]
CASE = sys.argv[2]
CASE_DIR = CASE.rsplit("/", 1)[0]
DP_BASE = 1.0e4
DP_DOUBLE = 2.0e4

RHO = 850.0
MU = 0.1

# name: (length, diameter, n, K, tau0); n = 0 means "no index declared" so the
# pipe runs the Bingham branch of 5.46, K = 0 means "no yield" so it runs the
# Metzner-Reed branch of 5.45, and tau0 = 1e-9 is the same pipe as HbZero with
# a threshold small enough to be a perturbation and large enough to be a branch
PIPES = {
    "HbN1":     (20.0, 0.02, 1.0, 0.1, 1.0),
    "Bing":     (20.0, 0.02, 0.0, 0.0, 1.0),
    "HbThin":   (20.0, 0.02, 0.5, 1.0, 0.5),
    "HbThin2":  (20.0, 0.02, 0.5, 1.0, 2.0),
    "HbThick":  (20.0, 0.02, 1.3, 0.02, 1.0),
    "HbNear":   (20.0, 0.02, 0.7, 2.0, 2.4),
    "HbStuck":  (20.0, 0.02, 0.7, 2.0, 2.6),
    "HbD2x":    (20.0, 0.04, 0.5, 1.0, 0.5),
    "HbZero":   (20.0, 0.02, 0.5, 1.0, 0.0),
    "HbTiny":   (20.0, 0.02, 0.5, 1.0, 1.0e-9),
    "HbS1":     (20.0, 0.02, 0.7, 2.0, 1.0),
    "HbS2":     (20.0, 0.02, 0.7, 2.0, 1.0),
    "Newt":     (20.0, 0.02, 0.0, 0.0, 0.0),
}
HB = [n for n, p in PIPES.items() if p[2] > 0 and p[4] > 0]
FIELDS = ("Q", "v", "Re", "f", "dp")

failures = []


def require(cond, what):
    if cond:
        print("ok:  ", what)
    else:
        print("FAIL:", what)
        failures.append(what)


def run(case_path):
    return subprocess.run([BINARY, case_path], capture_output=True, text=True,
                          cwd=CASE_DIR)


def parse(stdout):
    """`pipe <name> Q <q> v <u> Re <re> f <f> dp <dp>`, token-paired, plus the
    free node's pressure. Names come from the caller's table because the series
    pair is not in the parallel group."""
    out = {}
    nodes = {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 4 and tok[0] == "node" and tok[2] == "p":
            nodes[tok[1]] = float(tok[3])
            continue
        if len(tok) < 2 or tok[0] != "pipe" or tok[1] not in PIPES:
            continue
        name = tok[1]
        fields = {}
        i = 2
        while i + 1 < len(tok):
            key, val = tok[i], tok[i + 1]
            if key in FIELDS:
                try:
                    fields[key] = float(val)
                except ValueError:
                    break
            i += 2
        if not set(FIELDS) <= set(fields) or name in out:
            continue
        out[name] = fields
    return out, nodes


def tau_w(dp, length, diam):
    """The wall stress a drop asks for: geometry alone, no rheology."""
    return dp * diam / (4.0 * length)


def dp_static(length, diam, tau0):
    return 4.0 * tau0 * length / diam


def q_hb(dp, length, diam, n, k, tau0):
    """The claim itself: the termwise profile integral, read at this drop."""
    tw = tau_w(dp, length, diam)
    d = tw - tau0
    if d <= 0.0:
        return 0.0
    r = 0.5 * diam
    return (math.pi * r ** 3 * tw ** -3 * k ** (-1.0 / n) * n * d ** (1.0 + 1.0 / n)
            * (tau0 * tau0 / (n + 1.0) + 2.0 * tau0 * d / (2.0 * n + 1.0)
               + d * d / (3.0 * n + 1.0)))


def dp_of_q(q, length, diam, n, k, tau0):
    """The same closed form run backwards, by bisection on the monotone Q(tw):
    this is what the kernel's loss(q) answers with, so a drop rebuilt from a
    printed rate through here is the drop the kernel would have had to ask."""
    def q_at(tw):
        return q_hb(4.0 * tw * length / diam, length, diam, n, k, tau0)

    if q <= 0.0:
        return dp_static(length, diam, tau0)
    lo, hi = tau0, tau0 + 1.0
    while q_at(hi) < q:
        hi *= 2.0
    for _ in range(400):
        mid = 0.5 * (lo + hi)
        if q_at(mid) > q:
            hi = mid
        else:
            lo = mid
    return 4.0 * 0.5 * (lo + hi) * length / diam


def sens_hb(dp, length, diam, n, tau0):
    """d ln Q / d ln tw of the closed form, in closed form: with a = 1 + 1/n and
    A the bracket, Q ~ tw^-3 D^a A(D), so the logarithmic slope is
    -3 + (a + D A'/A) tw/D. It is the multiplier between a leftover drop error
    and the rate error it shows up as, and it is large exactly where the landing
    is soft -- which is why two pipes out of ten are off by 1e-9 and the rest by
    1e-12."""
    tw = tau_w(dp, length, diam)
    d = tw - tau0
    a = 1.0 + 1.0 / n
    aa = (tau0 * tau0 / (n + 1.0) + 2.0 * tau0 * d / (2.0 * n + 1.0)
          + d * d / (3.0 * n + 1.0))
    daa = 2.0 * tau0 / (2.0 * n + 1.0) + 2.0 * d / (3.0 * n + 1.0)
    return -3.0 + (a + d * daa / aa) * tw / d


def _gl(m):
    """Gauss-Legendre nodes and weights on [-1, 1]: Newton on P_m from the
    Christoffel-Darboux recurrence, with the asymptotic starting guess. numpy is
    not a test dependency, so the rule is built here and self-checked by the
    ladder in family 3 -- its nodes are the standard ones (2 gives +-1/sqrt(3)
    with weight 1) and its weights sum to 2.

    The recurrence is the point of the whole function: P_0 = 1, P_1 = z and
    j P_j = (2j-1) z P_{j-1} - (j-1) P_{j-2}, so after j runs 2..m the pair is
    (P_{m-1}, P_m) and P_m'(z) = m (z P_m - P_{m-1}) / (z^2 - 1) is built from
    the same two values. Off by one in that loop and the roots are of P_{m+1}
    while the weights still use m, which returns 20 well-formed-looking nodes
    whose weights do not sum to 2."""
    xs, ws = [], []
    for i in range(1, m + 1):
        z = math.cos(math.pi * (i - 0.25) / (m + 0.5))
        pp = 1.0
        for _ in range(100):
            p0, p1 = 1.0, z
            for j in range(2, m + 1):
                p0, p1 = p1, ((2.0 * j - 1.0) * z * p1 - (j - 1.0) * p0) / j
            pp = m * (z * p1 - p0) / (z * z - 1.0)
            step = p1 / pp
            z -= step
            if abs(step) < 1e-16:
                break
        xs.append(z)
        ws.append(2.0 / ((1.0 - z * z) * pp * pp))
    return xs, ws


_GL_CACHE = {}


def _panels(count):
    """[0, 1] cut geometrically, 2^-count as the first edge: the integrand
    behaves like s^(1/n), so its derivative at s = 0 is infinite whenever
    1/n < 1 and a single rule over the whole interval pays for that."""
    edges = [0.0] + [2.0 ** -j for j in range(count, -1, -1)]
    return list(zip(edges[:-1], edges[1:]))


def q_hb_quadrature(dp, length, diam, n, k, tau0, count=20, nodes=12):
    """The rate straight off the definition, Q = pi R^3 tw^-3 int tau^2
    gamma(tau) d tau with gamma = ((tau - tau0)/K)^(1/n), integrated over the
    profile instead of term by term. This never sees the three coefficients
    (n+1, 2n+1, 3n+1) that q_hb asserts, so a mistyped one of them shows up
    here as a disagreement rather than as a matching pair of typos."""
    tw = tau_w(dp, length, diam)
    span = tw - tau0
    if span <= 0.0:
        return 0.0
    if nodes not in _GL_CACHE:
        _GL_CACHE[nodes] = _gl(nodes)
    zs, ws = _GL_CACHE[nodes]
    total = 0.0
    for a, b in _panels(count):
        half = 0.5 * (b - a)
        mid = 0.5 * (a + b)
        for z, w in zip(zs, ws):
            s = mid + half * z
            total += w * half * (tau0 + span * s) ** 2 * s ** (1.0 / n)
    r = 0.5 * diam
    return (math.pi * r ** 3 * tw ** -3 * k ** (-1.0 / n)
            * span ** (1.0 + 1.0 / n) * total)


def q_buck(dp, length, diam, mu, tau0):
    """5.46's Buckingham closed form."""
    tw = tau_w(dp, length, diam)
    if tw <= tau0:
        return 0.0
    xi = tau0 / tw
    return math.pi * diam ** 3 / (32.0 * mu) * tw * (1.0 - 4.0 * xi / 3.0
                                                     + xi ** 4 / 3.0)


def q_profile(dp, length, diam, n, k):
    """The exact power-law discharge from the profile, tau0 = 0."""
    tw = tau_w(dp, length, diam)
    return (math.pi * diam ** 3 / 8.0) * n / (3.0 * n + 1.0) * (tw / k) ** (1.0 / n)


def q_metzner_old(dp, length, diam, n, k):
    """What 5.45 shipped: the same Churchill call with the exponent-1
    consistency K (3n+1)/(4n), i.e. the turbulent-collapse convention."""
    return (math.pi * diam ** 3 / 32.0) * (
        n * dp * diam / ((3.0 * n + 1.0) * k * length)) ** (1.0 / n)


def g_xi(tw, mu, tau0):
    xi = tau0 / tw
    return 1.0 - 4.0 * xi / 3.0 + xi ** 4 / 3.0


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


base = run(f"{CASE_DIR}/case.dict")
dbl = run(f"{CASE_DIR}/case_dp2.dict")
require(base.returncode == 0 and dbl.returncode == 0,
        f"both driving differences solve and converge (rc {base.returncode},"
        f" {dbl.returncode})")
B, NB = parse(base.stdout)
D, ND = parse(dbl.stdout)

# --- 1. the declaration read back --------------------------------------------
require(len(B) == len(PIPES) and len(D) == len(PIPES),
        f"all thirteen pipes report a rate, a velocity, a Reynolds number and a"
        f" drop in both files ({len(B)} and {len(D)} parsed)")
PARALLEL = [n for n in PIPES if n not in ("HbS1", "HbS2")]
require(all(abs(B[n]["dp"] - DP_BASE) < 1e-9 for n in PARALLEL)
        and all(abs(D[n]["dp"] - DP_DOUBLE) < 1e-9 for n in PARALLEL),
        f"every parallel pipe carries exactly the imposed difference"
        f" ({DP_BASE:.0f} Pa then {DP_DOUBLE:.0f} Pa), and only the series pair"
        f" splits it -- the nodes are fixed, so no Newton solve stands between"
        f" the case and the closed form")

# --- 2. the closed form, both files ----------------------------------------
# The series pair really carries 5000 Pa at the first driving and 10000 at the
# second, so it is read at its own drop; the others are told what the fixed
# nodes impose.
SERIES = ("HbS1", "HbS2")
HB_PARALLEL = [name for name in HB if name not in SERIES]
# The two states whose drop is barely past their own threshold: see the
# bracket accounting below, they are the only ones the solver stops early on.
SOFT = ("HbThin2", "HbNear")
HB_STEEP = [name for name in HB_PARALLEL if name not in SOFT]


def q_err(table, name, dp):
    return rel(table[name]["Q"], q_hb(dp, *PIPES[name][:2], *PIPES[name][2:]))


worst_q = max(q_err(B, name, DP_BASE) for name in HB_STEEP)
worst_q2 = max(q_err(D, name, DP_DOUBLE) for name in HB_STEEP)
worst_qs = max(q_err(B, name, 5000.0) for name in SERIES)
worst_qs2 = max(q_err(D, name, 10000.0) for name in SERIES)
require(worst_q < 1e-11 and worst_q2 < 1e-11 and worst_qs < 1e-11
        and worst_qs2 < 1e-11,
        f"the {len(HB_STEEP)} steep yielded states and the series pair are the"
        f" Herschel-Bulkley profile integral -- one of the steep ones is the"
        f" dead branch, where both sides are exactly zero -- to {worst_q:.3e}"
        f" relative at {DP_BASE:.0f} Pa, {worst_q2:.3e} at {DP_DOUBLE:.0f} Pa,"
        f" {worst_qs:.3e} on the series pair at its {5000.0:.0f} Pa and"
        f" {worst_qs2:.3e} at its doubled {10000.0:.0f} Pa: the law is the"
        f" integrated profile, not a fit to it")
worst_soft = max(q_err(B, name, DP_BASE) for name in SOFT)
worst_soft2 = max(q_err(D, name, DP_DOUBLE) for name in SOFT)
require(max(worst_soft, worst_soft2) < 1e-8,
        f"the two states family 2b calls soft land at"
        f" {max(worst_soft, worst_soft2):.3e} relative instead of at the"
        f" 1e-12 of the steep ones -- and 2b accounts for the whole of that gap"
        f" as the solver's stop, so the profile integral is not held to a"
        f" looser standard here than the other ten states are")

# --- 2b. what the leftover slack is made of ---------------------------------
soft_report = []
worst_slack = 0.0
for name in SOFT:
    length, diam, n, k, tau0 = PIPES[name]
    observed = q_err(B, name, DP_BASE)
    leftover = dp_of_q(B[name]["Q"], length, diam, n, k, tau0) - DP_BASE
    slope = sens_hb(DP_BASE, length, diam, n, tau0)
    predicted = abs(leftover) / DP_BASE * slope
    worst_slack = max(worst_slack, rel(observed, predicted))
    soft_report.append(f"{name}: {abs(leftover):.3e} Pa of drop left over x"
                       f" slope {slope:.3f} = {predicted:.3e}, observed"
                       f" {observed:.3e}")
require(worst_slack < 5e-3,
        f"the slack is fully accounted for, so nothing about the law is being"
        f" waived: rebuilt through the same closed form, the printed rates leave"
        f" {'; '.join(soft_report)}. equipmentFlow stops on a rate bracket 1e-9"
        f" wide relative (network.cpp:364) and the drop is the quantity the case"
        f" fixes, so what is left in the drop is multiplied by the landing's own"
        f" slope before it reaches the rate")

# --- 3. the same integral, not termwise ------------------------------------
QUAD = [(name, DP_BASE) for name in
        ("HbThin", "HbThick", "HbThin2", "HbNear")] + \
       [(name, 5000.0) for name in SERIES]
worst_quad = max(rel(B[name]["Q"],
                     q_hb_quadrature(dp, *PIPES[name][:2], *PIPES[name][2:]))
                 for name, dp in QUAD)
coarse = max(rel(q_hb_quadrature(dp, *PIPES[name][:2], *PIPES[name][2:],
                                 count=0),
                 q_hb_quadrature(dp, *PIPES[name][:2], *PIPES[name][2:],
                                 count=20))
             for name, dp in QUAD)
fine = max(rel(q_hb_quadrature(dp, *PIPES[name][:2], *PIPES[name][2:],
                               count=20),
               q_hb_quadrature(dp, *PIPES[name][:2], *PIPES[name][2:],
                               count=40))
           for name, dp in QUAD)
require(coarse > 1e-9 and fine < 1e-14,
        f"and the quadrature is resolved rather than merely different: the"
        f" same 12-node rule on one panel over s is {coarse:.3e} from the same"
        f" rule on the {len(_panels(20))} geometrically graded panels, because"
        f" s^(1/n) has an infinite derivative at the yield radius, but grading"
        f" twice as deep again moves nothing ({fine:.3e}), so what is compared"
        f" below is a converged integral")
require(worst_quad < 1e-8,
        f"and the closed form above is the right closed form: quadrating the"
        f" profile integral instead of integrating it term by term lands on the"
        f" same kernel numbers to {worst_quad:.3e} relative -- the two are"
        f" independent derivations of one statement, so a mistyped (3n+1) here"
        f" would not have a matching mistyped (3n+1) there")

# --- 4. limit n = 1: the branch 5.46 shipped -------------------------------
require(B["HbN1"]["Q"] == B["Bing"]["Q"] and B["HbN1"]["v"] == B["Bing"]["v"]
        and B["HbN1"]["f"] == B["Bing"]["f"],
        "n = 1 with K = mu declared as Herschel-Bulkley prints the *same* rate,"
        " velocity and friction factor as the Bingham branch: the two code"
        " paths -- Buckingham's Newton on a cubic and this one's bisection on"
        " the profile integral -- agree in the printed digits, which is what"
        " makes the new branch a generalisation of 5.46 rather than a rival")
tw25 = tau_w(DP_BASE, 20.0, 0.02)
g_ref = g_xi(tw25, MU, 1.0)
require(abs(B["HbN1"]["Re"] / (B["Bing"]["Re"] * g_ref) - 1.0) < 1e-9,
        f"the two Reynolds numbers that print beside them differ by exactly"
        f" g(xi) = {g_ref:.6f} ({B['HbN1']['Re']:.6g} vs"
        f" {B['Bing']['Re']:.6g}): the plastic number is rho v d / mu and the"
        f" apparent one is rho v d / (tw d / 8v), and the ratio between those"
        f" two viscosities is the rate defect of the yield -- a yield stress"
        f" with an index has no material viscosity to report, so the apparent"
        f" number is the one this branch prints")

# --- 5. limit tau0 = 0: the branch 5.45 shipped, corrected -----------------
require(abs(B["HbZero"]["Q"] - q_profile(DP_BASE, 20.0, 0.02, 0.5, 1.0))
        / q_profile(DP_BASE, 20.0, 0.02, 0.5, 1.0) < 1e-9,
        "the no-yield pipe at n = 0.5 is the profile integral too, to "
        f"{rel(B['HbZero']['Q'], q_profile(DP_BASE, 20.0, 0.02, 0.5, 1.0)):.3e}"
        " -- after the consistency constant is corrected, the Metzner-Reed"
        " laminar line *is* the exact power-law solution that 5.45 claimed it"
        " was")
require(abs(B["HbTiny"]["Q"] / B["HbZero"]["Q"] - 1.0) < 1e-8,
        f"and the same pipe with a 1e-9 Pa threshold, routed through the yield"
        f" branch instead, gives {B['HbTiny']['Q']:.12g} against"
        f" {B['HbZero']['Q']:.12g} -- {rel(B['HbTiny']['Q'], B['HbZero']['Q']):.3e}"
        " apart, which is the size of the threshold it declares (tau0 / tw ="
        f" {1.0e-9 / tw25:.1e}), not the two branches disagreeing")
ratio_old = q_profile(DP_BASE, 20.0, 0.02, 0.5, 1.0) / q_metzner_old(
    DP_BASE, 20.0, 0.02, 0.5, 1.0)
predict = (4.0 * 0.5 / (3.0 * 0.5 + 1.0)) ** (1.0 - 1.0 / 0.5)
require(abs(ratio_old - predict) < 1e-12 and abs(ratio_old - 1.25) < 1e-12,
        f"the size of the constant 5.45 had wrong is readable: the profile"
        f" integral is {ratio_old:.6f}x its closed form at n = 0.5, which is"
        f" (4n/(3n+1))^(1-1/n) = {predict:.6f} exactly -- 25% of rate, or 20%"
        f" of drop, is what the exponent-1 convention cost, and it vanishes at"
        f" n = 1 (where both are 1) which is why 5.45's own control, PowN1 vs"
        f" Newt, was bit-identical and never saw it")

# --- 6. the threshold is geometric, the branch is asked ---------------------
stk_dp = dp_static(20.0, 0.02, 2.6)
require(B["HbStuck"]["Q"] == 0.0 and B["HbStuck"]["v"] == 0.0
        and B["HbStuck"]["Re"] == 0.0 and B["HbStuck"]["f"] == 0.0
        and abs(B["HbStuck"]["dp"] - DP_BASE) < 1e-9 and stk_dp > DP_BASE,
        f"the n = 0.7 pipe whose threshold is {stk_dp:.0f} Pa prints zero rate,"
        f" zero velocity, zero Reynolds number and zero friction factor while"
        f" the full {DP_BASE:.0f} Pa sits across it: below dp_static = 4 tau0 L/d"
        f" the exact answer is no flow, and it takes everything asked of it")
require(D["HbStuck"]["Q"] > 0.0
        and rel(D["HbStuck"]["Q"], q_hb(DP_DOUBLE, 20.0, 0.02, 0.7, 2.0, 2.6))
        < 1e-9,
        f"the same pipe flows at the doubled difference"
        f" ({D['HbStuck']['Q']:.6g} m3/s, closed form to"
        f" {rel(D['HbStuck']['Q'], q_hb(DP_DOUBLE, 20.0, 0.02, 0.7, 2.0, 2.6)):.3e}):"
        f" the dead branch is decided by the drop asked across the pipe, not by"
        f" a property of the pipe -- 5.46's statement with n != 1")
require(B["HbThin"]["Q"] > B["HbThin2"]["Q"]
        and dp_static(20.0, 0.02, 2.0) > dp_static(20.0, 0.02, 0.5),
        f"at the same (n, K) a higher threshold is a lower rate"
        f" ({B['HbThin']['Q']:.6g} vs {B['HbThin2']['Q']:.6g} m3/s),"
        f" and the ordering follows the threshold that enters with the geometry"
        f" alone: 4 tau0 L/d is {dp_static(20.0, 0.02, 0.5):.0f} Pa then"
        f" {dp_static(20.0, 0.02, 2.0):.0f} Pa out of a {DP_BASE:.0f} Pa drop")

# --- 7. the two printed numbers are consistent with each other --------------
worst_fr = max(abs(B[n]["f"] * B[n]["Re"] - 64.0) for n in HB if n != "HbStuck")
worst_fr2 = max(abs(D[n]["f"] * D[n]["Re"] - 64.0) for n in HB)
require(max(worst_fr, worst_fr2) < 1e-6,
        f"f * Re == 64 on every flowing yielded state (worst"
        f" {max(worst_fr, worst_fr2):.3e}) -- but here that is an identity"
        f" between the two reports, not evidence: f is 8 tw/(rho v^2) and Re is"
        f" 8 rho v^2/tw, so the product is 64 for any tw whatever. It is"
        f" asserted because a mismatched pair would mean the log printed a"
        f" number from a different state than the drop beside it")
worst_dpf = max(rel(B[n]["dp"], B[n]["f"] * (PIPES[n][0] / PIPES[n][1])
                    * RHO * B[n]["v"] ** 2 / 2.0)
                for n in HB if n != "HbStuck")
require(worst_dpf < 1e-9,
        f"and the reported factor does reproduce the reported drop:"
        f" dp == f (L/d) rho v^2/2 on all {len(HB) - 1} yielded states to"
        f" {worst_dpf:.3e} -- the twelve-significant-digit print floor is what"
        f" this is measuring, so the friction number is the kernel's own")

# --- 8. the series pair -----------------------------------------------------
require(abs(NB["C"] - 25000.0) < 1e-9 and abs(ND["C"] - 30000.0) < 1e-9,
        f"the free node between two identical yielded pipes lands exactly"
        f" halfway ({NB['C']:.6g} Pa of 2.0e4 and 1.0e4, {ND['C']:.6g} Pa"
        f" doubled): symmetry is a property of the split, not of the driving")
require(B["HbS1"] == B["HbS2"],
        "and the two halves print bit-for-bit the same five fields, including"
        " the friction factor -- the network Newton is solving across a law"
        " whose loss does not vanish with its rate, and the dead branch and the"
        " profile integral come out symmetric rather than leaking")
require(D["HbS1"]["Q"] > B["HbS1"]["Q"],
        f"the pair carries more when the difference is doubled"
        f" ({D['HbS1']['Q']:.6g} vs {B['HbS1']['Q']:.6g} m3/s), and the index"
        f" is the yield law's own: each half goes from D = 0.25 Pa to 1.25 Pa"
        f" above its threshold, so the soft landing (Q ~ D^(1+1/n)) is being"
        f" read at two points rather than asserted at one")

# --- 9. the refusal that survives ------------------------------------------
r1 = run(f"{CASE_DIR}/case_rejectNoConsistency.dict")
require(r1.returncode != 0
        and "flow index" in (r1.stdout + r1.stderr).lower()
        and "consistency" in (r1.stdout + r1.stderr).lower(),
        "an index with a threshold and no consistency is still refused at parse"
        " time: 5.47 lifts the refusal that said Herschel-Bulkley had no"
        " inverted law, and keeps the one that says a K in Pa.s^n cannot"
        " default to a mu in Pa.s")

print()
print(f"{len(failures)} failure(s)")
sys.exit(1 if failures else 0)
