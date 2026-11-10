#!/usr/bin/env python3
"""steadyBinghamYield: a yield stress in the steady network, checked against
Buckingham's 1921 exact laminar solution and against the static threshold that
comes with it.

A power-law fill (5.45) only needed a substitution inside Churchill. A yield
stress cannot be smuggled in that way: the drop stops vanishing with the rate.
At q = 0 the wall still holds

    dpStatic = 4 tau_0 L / d

and below that the exact answer is no flow at all -- which is the one case the
Illinois bracket cannot see (g(0) already sits on the far side of the
crossing), so the kernel decides it before bracketing. Above the threshold the
law is the exact one, written for the wall stress tw = dp d / 4L with
xi = tau_0 / tw:

    Q = (pi d^3 / (32 mu)) tw [1 - (4/3) xi + (1/3) xi^4]
    f Re g(xi) = 64,   g(xi) in (0, 1],   Q = Q_Hagen-Poiseuille * g(xi)

Nine families of judgement:
  1. the declaration read back (twelve pipes, the imposed drop)
  2. Q against Buckingham, at the imposed drop and at half of it
  3. the sentinel: yieldStress 0 is bit-identical to declaring nothing
  4. the static threshold: q == 0 exactly, and Q > 0 iff dp > dpStatic
  5. f Re g(xi) == 64 -- the reported factor is the solution, not a correlation
  6. the yield multiplies the Newtonian rate by exactly g(xi)
  7. the wall stress is still geometry: tw does not contain tau_0
  8. the dead branch through the network Newton, plus the inverse problem
  9. the two parse-time refusals that are still refusals, and the one 5.47 lifted
"""

import math
import subprocess
import sys

BINARY = sys.argv[1]
CASE = sys.argv[2]
CASE_DIR = CASE.rsplit("/", 1)[0]

RHO = 850.0
MU = 0.1
DP_BASE = 1.0e4           # p(A) 3e4 - p(B) 2e4, every parallel pipe
DP_SERIES = 5.0e3         # Sr1: two identical yielded pipes split the drop
Q_DRIVEN = 5.0e-6         # case_driven: the rate is what is imposed

# name: (length, diameter, tau0); tau0 = 0 means no threshold, and for Newt the
# keyword itself is absent, which is what makes that pair a sentinel control.
PIPES = {
    "Newt":   (20.0, 0.02, 0.0),
    "T0zero": (20.0, 0.02, 0.0),
    "Xi02":   (20.0, 0.02, 0.5),
    "Xi04":   (20.0, 0.02, 1.0),
    "Xi06":   (20.0, 0.02, 1.5),
    "Xi08":   (20.0, 0.02, 2.0),
    "Xi096":  (20.0, 0.02, 2.4),
    "Stk":    (20.0, 0.02, 2.6),
    "Sr1a":   (20.0, 0.02, 1.0),
    "Sr1b":   (20.0, 0.02, 1.0),
    "Sr2a":   (20.0, 0.02, 2.6),
    "Sr2b":   (20.0, 0.02, 0.0),
}
DRIVEN = {
    "Yf": (20.0, 0.02, 1.0),
    "Yr": (20.0, 0.02, 1.0),
}
PAR = ("Newt", "T0zero", "Xi02", "Xi04", "Xi06", "Xi08", "Xi096", "Stk")
YLD = ("Xi02", "Xi04", "Xi06", "Xi08", "Xi096")     # yielded and flowing
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


def parse(stdout, known):
    """`pipe <name> Q <q> v <u> Re <re> f <f> dp <dp> ...` read as key/value
    pairs, so a field the app starts printing after dp cannot silently shift
    the numbers."""
    out = {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) < 2 or tok[0] != "pipe" or tok[1] not in known:
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
    return out


def parse_nodes(stdout):
    out = {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) == 4 and tok[0] == "node" and tok[2] == "p":
            try:
                out[tok[1]] = float(tok[3])
            except ValueError:
                pass
    return out


def q_hp(dp, length, diam):
    return math.pi * diam ** 4 * dp / (128.0 * MU * length)


def dp_static(length, diam, tau0):
    """The drop a stuck Bingham wall holds: 4 tau_0 L / d."""
    return 4.0 * tau0 * length / diam


def g_xi(xi):
    return 1.0 - 4.0 * xi / 3.0 + xi ** 4 / 3.0


def tauw_of(dp, length, diam):
    return dp * diam / (4.0 * length)


def dp_of(tauw, length, diam):
    return 4.0 * tauw * length / diam


def q_bingham(dp, length, diam, tau0):
    """Buckingham's discharge at an imposed drop; zero below the threshold,
    which is where the Hagen-Poiseuille limit of the same expression would have
    quietly kept flowing."""
    if tau0 == 0.0:
        return q_hp(dp, length, diam)
    if dp <= dp_static(length, diam, tau0):
        return 0.0
    tw = tauw_of(dp, length, diam)
    return (math.pi * diam ** 3 / (32.0 * MU)) * tw * g_xi(tau0 / tw)


def tauw_from_q(q, length, diam, tau0):
    """The same law inverted the other way: bisection on the monotone
    F(tw) = tw - 4 tau0/3 + tau0^4/(3 tw^3) - 32 mu Q/(pi d^3) over the
    analytic bracket [tau0, B + 4 tau0/3]. Deliberately not the kernel's
    Newton, so a wrong root here would not share a bug with the code."""
    b = 32.0 * MU * abs(q) / (math.pi * diam ** 3)
    if tau0 == 0.0:
        return b
    if b == 0.0:
        return tau0
    lo, hi = tau0, b + 4.0 * tau0 / 3.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f = mid - 4.0 * tau0 / 3.0 + tau0 ** 4 / (3.0 * mid ** 3) - b
        if f > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def dp_from_q(q, length, diam, tau0):
    return dp_of(tauw_from_q(q, length, diam, tau0), length, diam)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


base = run(f"{CASE_DIR}/case.dict")
drv = run(f"{CASE_DIR}/case_driven.dict")
require(base.returncode == 0 and drv.returncode == 0,
        f"both files solve and converge (rc {base.returncode},"
        f" {drv.returncode})")
B = parse(base.stdout, PIPES)
N = parse_nodes(base.stdout)
V = parse(drv.stdout, DRIVEN)
TAUS = {n: tauw_of(DP_BASE, PIPES[n][0], PIPES[n][1]) for n in PAR}

# --- 1. the declaration read back -------------------------------------------
require(len(B) == len(PIPES) and len(V) == len(DRIVEN),
        f"all {len(PIPES)} pipes report Q, v, Re, f and dp in case.dict and both"
        f" driven pipes do in case_driven ({len(B)} and {len(V)} parsed)")
worst_dp_par = max(abs(B[n]["dp"] - DP_BASE) for n in PAR)
require(worst_dp_par < 1e-9,
        f"the eight parallel pipes carry exactly the imposed {DP_BASE:.0f} Pa"
        f" (worst {worst_dp_par:.3e}): both ends are fixed-pressure nodes, so"
        f" each of those statements is the pipe law alone and not a Newton"
        f" artefact")

# --- 2. against Buckingham ---------------------------------------------------
worst_q = max(rel(B[n]["Q"], q_bingham(DP_BASE, *PIPES[n])) for n in YLD)
worst_q_series = max(rel(B[n]["Q"], q_bingham(DP_SERIES, *PIPES[n]))
                     for n in ("Sr1a", "Sr1b"))
worst_q_newt = max(rel(B[n]["Q"], q_hp(DP_BASE, PIPES[n][0], PIPES[n][1]))
                   for n in ("Newt", "T0zero"))
require(worst_q < 1e-9,
        f"the five flowing yielded rates are Buckingham's solution to "
        f"{worst_q:.3e} relative -- the kernel solved the wall-stress bracket,"
        f" it did not fit a correction onto the Newtonian law")
require(worst_q_series < 1e-9,
        f"and the same holds at half the drop, where the drop itself was an"
        f" unknown the Newton solve had to find: {worst_q_series:.3e} relative"
        f" on both series pipes")
require(worst_q_newt < 1e-9,
        f"the two no-threshold pipes stay Hagen-Poiseuille to {worst_q_newt:.3e}"
        f" (Churchill's laminar branch is 64/Re at this Reynolds number), which"
        f" is the baseline every g(xi) ratio below divides by")

# --- 3. the sentinel ---------------------------------------------------------
same = all(B["Newt"][k] == B["T0zero"][k] for k in FIELDS)
require(same,
        "yieldStress 0 is bit-identical to not declaring the keyword, in all"
        " five fields -- tau_0 -> 0 collapses F(tw) = tw onto the Newtonian"
        " slope term for term, so every other test in the suite cannot be"
        " moved by this branch and the case has a control rather than an"
        " assumption")

# --- 4. the static threshold -------------------------------------------------
stk_dp = dp_static(*PIPES["Stk"])
require(B["Stk"]["Q"] == 0.0 and B["Stk"]["v"] == 0.0
        and B["Stk"]["Re"] == 0.0 and B["Stk"]["f"] == 0.0,
        f"the stuck pipe reports exactly zero: Q, v, Re and f all print 0 while"
        f" dp is {B['Stk']['dp']:.0f} Pa against a threshold of {stk_dp:.0f} Pa"
        f" -- no residual rate left behind by a tolerance, and the bracket was"
        f" never entered")
require(B["Stk"]["dp"] < stk_dp and B["Stk"]["dp"] == DP_BASE,
        "and a stuck wall takes the drop it is asked for, not the threshold: dp"
        f" stays {DP_BASE:.0f} Pa while it could hold up to {stk_dp:.0f} Pa"
        " -- static friction, not a capped friction")
sign_ok = all((B[n]["Q"] > 0.0) == (DP_BASE > dp_static(*PIPES[n]))
              for n in YLD + ("Stk",))
thresholds = ", ".join(f"{n} {dp_static(*PIPES[n]):.0f}" for n in YLD + ("Stk",))
require(sign_ok,
        f"Q > 0 exactly when dp > dpStatic, over all six yielded parallel pipes"
        f" (thresholds in Pa: {thresholds}) -- the branch is decided by the"
        f" physics of the threshold, not by where the seed happened to land")
mono = all(B[a]["Q"] > B[b]["Q"] for a, b in
           (("Xi02", "Xi04"), ("Xi04", "Xi06"), ("Xi06", "Xi08"),
            ("Xi08", "Xi096")))
require(mono and B["Xi096"]["Q"] > 0.0,
        f"the rate falls as the threshold rises, and xi = 0.96 still carries"
        f" {B['Xi096']['Q']:.6g} m3/s -- the near-threshold branch is resolved"
        f" to the last digits rather than flushed to the left end of the"
        f" bracket, which is where a loose inner solve would hide")

# --- 5. the reported friction factor is the solution ------------------------
worst_id = 0.0
worst_dp_from_f = 0.0
for n in YLD:
    g = g_xi(PIPES[n][2] / TAUS[n])
    worst_id = max(worst_id, abs(B[n]["f"] * B[n]["Re"] * g - 64.0))
    dp_from_f = B[n]["f"] * PIPES[n][0] / PIPES[n][1] * RHO * B[n]["v"] ** 2 / 2
    worst_dp_from_f = max(worst_dp_from_f, rel(B[n]["dp"], dp_from_f))
require(worst_id < 1e-8,
        f"f * Re * g(xi) == 64 on all five flowing yielded pipes (worst "
        f"{worst_id:.3e} away from 64) -- the printed factor is the one"
        f" Buckingham's drop contains (f = 8 tw / rho v^2), not a Churchill"
        f" number pasted next to a yielded wall")
require(worst_dp_from_f < 1e-10,
        f"and it reproduces its own drop: dp == f (L/d) rho v^2 / 2 to "
        f"{worst_dp_from_f:.3e}, so the factor and the pressure on one log line"
        f" are a single statement rather than two that happen to be near")

# --- 6. the yield multiplies the Newtonian rate by g(xi) -------------------
ratios = {n: B[n]["Q"] / q_hp(DP_BASE, PIPES[n][0], PIPES[n][1]) for n in YLD}
worst_ratio = max(rel(ratios[n], g_xi(PIPES[n][2] / TAUS[n])) for n in YLD)
worst_g_bound = max(abs(g_xi(PIPES[n][2] / TAUS[n]) - 1.0) for n in YLD)
require(worst_ratio < 1e-9,
        f"Q / Q_Hagen-Poiseuille is exactly g(xi) on every yielded pipe (worst "
        f"{worst_ratio:.3e}): "
        + ", ".join(f"xi={PIPES[n][2] / TAUS[n]:g} -> {ratios[n]:.4f}"
                    for n in YLD)
        + " -- the threshold does not bend the parabolic profile, it deletes its"
          " core, and the deleted fraction is a function of xi alone")
require(all(0.0 < ratios[n] <= 1.0 for n in YLD) and worst_g_bound < 1.0,
        "every ratio sits in (0, 1], which is the same statement as 'a yield"
        " stress can only take flow away' -- above 1 would mean the threshold"
        " was reducing the drop, and exactly 1 would mean it did nothing")

# --- 7. the wall stress is geometry ----------------------------------------
require(max(abs(t - TAUS["Newt"]) for t in TAUS.values()) == 0.0,
        f"all eight parallel pipes ask for the same wall stress"
        f" {TAUS['Newt']:.3f} Pa -- tau_w = dp d / 4L has no tau_0 in it, so the"
        f" threshold moves the shear rate needed to hold the stress, never the"
        f" stress itself")
worst_tw = max(rel(TAUS[n], B[n]["f"] * RHO * B[n]["v"] ** 2 / 8.0)
               for n in YLD)
require(worst_tw < 1e-10,
        f"and each flowing rate really does hold that stress: 8 tw / (rho v^2)"
        f" rebuilds the reported f to {worst_tw:.3e}, so the fluid is free to"
        f" choose the rate but not the wall")

# --- 8. the dead branch through the Newton, and the inverse problem --------
require(N["C"] == 25000.0,
        f"the two identical yielded series pipes leave the free middle node at"
        f" exactly {N['C']:g} Pa -- symmetry, but symmetry the network Newton"
        f" had to earn through the yielded law and its conductance")
require(B["Sr1a"]["Q"] == B["Sr1b"]["Q"],
        f"and both carry the same rate {B['Sr1a']['Q']:.6g} m3/s to the last"
        f" printed bit: that is the continuity equation at the node, not a"
        f" tolerance")
require(N["D"] == 20000.0 and B["Sr2a"]["Q"] == 0.0 and B["Sr2b"]["Q"] == 0.0,
        f"a stuck pipe in series stops the whole line: the free node is pulled"
        f" onto the downstream fixed pressure ({N['D']:g} Pa) and both pipes"
        f" report exactly zero -- the plateau is what the iterate landed on,"
        f" not what it was handed")
require(B["Sr2a"]["dp"] == DP_BASE and B["Sr2a"]["dp"] < dp_static(*PIPES["Sr2a"]),
        f"and it stays stuck with the entire drop across it: dp"
        f" {B['Sr2a']['dp']:.0f} Pa against a threshold of"
        f" {dp_static(*PIPES['Sr2a']):.0f} Pa -- the series solution has no flow"
        f" left to leak through the bracket")
worst_inv = rel(V["Yf"]["dp"], dp_from_q(Q_DRIVEN, *DRIVEN["Yf"]))
require(worst_inv < 1e-10,
        f"the inverse problem agrees: impose {Q_DRIVEN:.0e} m3/s and the drop"
        f" the kernel solves for is {V['Yf']['dp']:.6g} Pa, matching the"
        f" checker's bisection inversion to {worst_inv:.3e} -- dp-known and"
        f" Q-known are one monotone law read from either end")
require(V["Yf"]["Q"] == -V["Yr"]["Q"] and V["Yf"]["dp"] == -V["Yr"]["dp"],
        f"and the law is even in the rate: the reversed pipe prints Q"
        f" {V['Yr']['Q']:.6g} at dp {V['Yr']['dp']:.6g} Pa, bit-for-bit the"
        f" negative of the forward one -- the dead branch is decided on the"
        f" sign of the chosen direction, not by luck in the forward case")

# --- 9. refusals, and the one that has been lifted --------------------------
r1 = run(f"{CASE_DIR}/case_rheologyNowSolves.dict")
Y1 = parse(r1.stdout, {"Bad"})
require(r1.returncode == 0 and Y1["Bad"]["Q"] > 0.0
        and Y1["Bad"]["dp"] == DP_BASE,
        f"the refusal this case used to carry is lifted by 5.47, and the"
        f" declaration is unchanged: a threshold together with an index and a"
        f" consistency now solves to {Y1['Bad']['Q']:.6g} m3/s across the full"
        f" {DP_BASE:.0f} Pa instead of being refused at parse time -- at n != 1"
        f" the bracket is no longer Buckingham's cubic but the profile integral,"
        f" which is monotone in tw for every n, so the case has its own inverted"
        f" law. What that number *is* belongs to 5.47's checker, not here")
r2 = run(f"{CASE_DIR}/case_rejectCoupled.dict")
r3 = run(f"{CASE_DIR}/case_rejectNegative.dict")
require(r2.returncode != 0 and "coupled" in (r2.stdout + r2.stderr).lower(),
        "a coupled element that also declares a yield stress is refused: its"
        " loss already is the partner model's statement, so there is no wall"
        " here for a threshold to hold")
require(r3.returncode != 0 and "negative" in (r3.stdout + r3.stderr).lower(),
        "a negative yield stress is refused instead of read as a threshold that"
        " helps the flow: tau_0 < 0 pushes g(xi) past 1, the one way this law"
        " could return a drop smaller than the Newtonian one")

print()
print(f"{len(failures)} failure(s)")
sys.exit(1 if failures else 0)
