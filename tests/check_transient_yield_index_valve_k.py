#!/usr/bin/env python3
"""transientYieldIndexValveK (C2): the Herschel-Bulkley wall (5.47) and the
bare local loss `valve K` (architecture §4.2, first measured with a Bingham
wall in 5.50) on the same momentum row, with NO schedule -- the one cell of the
(valve, schedule, index) table that no earlier case asked. Checked by
re-deriving every reported step from both terms written out in Python.

The row is

    rho L/A * (Q_new - Q_prev)/dt = drive - dp_wall(Q_new) - K rho |Q| Q / (2 A^2),
    drive = pFrom - pTo - rho g dz

with dp_wall the 5.47 profile integral inverted, Q = q_hb(tw, n, K_cons),
tw = dp_wall d / (4L). The valve term is the same one 5.50 charged a Bingham
wall with; what is new here is that the wall law carries an index, so the
valve's share of the asked drop is set by n and K_cons, not only by K:

  * at rest the valve is nothing (dp_valve(0) = 0): five decades of K on a held
    HB wall print one string, and the static-hold threshold carries no n, no
    K_cons and no K -- a tau_0 = 1.0 pipe prints 4000 Pa at n = 0.5 AND at n = 1.5;
  * where the law is read the valve is charged, and the same K = 500 that takes
    2.33 Pa of the n = 0.5 pipe's 10000 Pa takes 308.76 Pa of the n = 1.5 pipe's,
    because the shear-thickening wall flows 12x more and the quadratic term
    scales with Q^2;
  * the stop is the row's arithmetic, so the path decides the step: the checker
    reproduces all three Stop traces (including which step each stops on) from
    the recursion, and the steady kernel meets the marched row at t = 4 on the
    root of the same two-term loss, with the printed split dpWall + dpValve == dp.

Seven families of judgement:
  1. the wall law and the threshold are both reported, and the threshold has
     no n, no K_cons and no K in it
  2. a held pipe is the same pipe at any valve K, in both kernels
  3. every reported state of every flowing pipe against the recursion with both
     terms (HB profile integral + quadratic valve)
  4. the terminal rate is the root of the two-term loss; monotone in K
  5. the index moves the valve's share, not only the wall's rate
  6. the stop: three paths, three stop steps, the latch picks up after each
  7. the steady kernel: same root, printed split adds up through the inversion
"""

import math
import re
import subprocess
import sys

TRANSIENT = sys.argv[1]
STEADY = sys.argv[2]
CASE = sys.argv[3]
CASE_STEADY = sys.argv[4]
CASE_DIR = CASE.rsplit("/", 1)[0]

RHO = 850.0
MU = 0.1
GRAV = 9.80665
DT = 0.002
END = 4.0
NS = int(round(END / DT)) + 1
LENGTH = 20.0
DIAM = 0.02
AREA = 0.25 * math.pi * DIAM ** 2
R_PIPE = 0.5 * DIAM
INERTIA = RHO * LENGTH / AREA
I_DT = INERTIA / DT
KV_SCALE = RHO / (2.0 * AREA * AREA)     # dp_valve = K * KV_SCALE * Q|Q|
WALL_SCALE = 4.0 * LENGTH / DIAM         # dp_wall = WALL_SCALE * tw
P_A, P_B, P_A2 = 3.0e4, 2.0e4, 2.4e4

# name: (n, K_cons, tau0, valveK, q0, p_from, p_to, dz)
PIPES = {
    "StuckRef":      (0.5, 1.0, 2.6, 0.0, 0.0, P_A, P_B, 0.0),
    "StuckK500":     (0.5, 1.0, 2.6, 500.0, 0.0, P_A, P_B, 0.0),
    "StuckK1e8":     (0.5, 1.0, 2.6, 1.0e8, 0.0, P_A, P_B, 0.0),
    "RelHB05":       (0.5, 1.0, 1.0, 0.0, 0.0, P_A, P_B, 0.0),
    "RelHB05K100":   (0.5, 1.0, 1.0, 100.0, 0.0, P_A, P_B, 0.0),
    "RelHB05K500":   (0.5, 1.0, 1.0, 500.0, 0.0, P_A, P_B, 0.0),
    "RelHB05K5000":  (0.5, 1.0, 1.0, 5000.0, 0.0, P_A, P_B, 0.0),
    "RelHB15":       (1.5, 0.02, 1.0, 0.0, 0.0, P_A, P_B, 0.0),
    "RelHB15K500":   (1.5, 0.02, 1.0, 500.0, 0.0, P_A, P_B, 0.0),
    "HangHB05":      (0.5, 1.0, 0.5, 0.0, 0.0, P_A, P_B, 1.0),
    "HangHB05K500":  (0.5, 1.0, 0.5, 500.0, 0.0, P_A, P_B, 1.0),
    "RollHB05":      (0.5, 1.0, 0.5, 0.0, 0.0, P_A, P_B, -1.0),
    "RollHB05K500":  (0.5, 1.0, 0.5, 500.0, 0.0, P_A, P_B, -1.0),
    "StopRef":       (0.5, 1.0, 1.5, 0.0, 2.0e-6, P_A2, P_B, 0.0),
    "StopK5000":     (0.5, 1.0, 1.5, 5000.0, 2.0e-6, P_A2, P_B, 0.0),
    "StopK1e8":      (0.5, 1.0, 1.5, 1.0e8, 2.0e-4, P_A2, P_B, 0.0),
}
HELD = ("StuckRef", "StuckK500", "StuckK1e8", "HangHB05", "HangHB05K500")
FLOWING = ("RelHB05", "RelHB05K100", "RelHB05K500", "RelHB05K5000",
           "RelHB15", "RelHB15K500", "RollHB05", "RollHB05K500")
STOPPING = ("StopRef", "StopK5000", "StopK1e8")

failures = []


def require(cond, what):
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


def parse_transient(stdout):
    """Q trace (value and printed token), latch per pipe, thresholds, the
    wall-law declaration line, times."""
    out = {"q": {n: [] for n in PIPES}, "raw": {n: [] for n in PIPES},
           "hold": {n: [] for n in PIPES}, "thr": {}, "law": {}, "times": []}
    started = False
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 2 and tok[0] == "t":
            out["times"].append(float(tok[1]))
            started = True
        elif not started and len(tok) == 5 and tok[0] == "pipe" \
                and tok[2] == "static-hold" and tok[3] == "threshold":
            out["thr"][tok[1]] = float(tok[4])
        elif not started and len(tok) == 7 and tok[0] == "pipe" \
                and tok[2] == "wall-law" and tok[3] == "flowIndex":
            out["law"][tok[1]] = (float(tok[4]), float(tok[6]))
        elif len(tok) == 4 and tok[0] == "pipe" and tok[2] == "Q":
            if tok[1] in out["q"]:
                out["q"][tok[1]].append(float(tok[3]))
                out["raw"][tok[1]].append(tok[3])
        elif len(tok) == 5 and tok[0] == "pipe" and tok[3] == "hold":
            if tok[1] in out["hold"]:
                out["hold"][tok[1]].append(int(tok[4]))
    return out


def parse_steady(stdout):
    """One state per pipe, with the printed drop split."""
    out = {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) < 12 or tok[0] != "pipe" or tok[2] != "Q" \
                or tok[1] not in PIPES:
            continue
        rec = {"q": float(tok[3]), "v": float(tok[5])}
        for i, t in enumerate(tok):
            if t == "dp":
                rec["dp"] = float(tok[i + 1])
            elif t == "dpValve":
                rec["dpValve"] = float(tok[i + 1])
        out[tok[1]] = rec
    return out


# --- the two terms, written out from the case file --------------------------

def threshold(tau0):
    """4 tau_0 L / d -- no n, no K_cons, no valve K in it."""
    return WALL_SCALE * tau0


def q_hb(tw, n, k_cons, tau0):
    """5.47's laminar profile integral, the exact inverse wall of this case."""
    if tw <= tau0:
        return 0.0
    d = tw - tau0
    p = 1.0 / n
    return (math.pi * R_PIPE ** 3 * tw ** -3 * k_cons ** -p
            * (d ** (p + 3.0) / (p + 3.0)
               + 2.0 * tau0 * d ** (p + 2.0) / (p + 2.0)
               + tau0 * tau0 * d ** (p + 1.0) / (p + 1.0)))


def dp_valve(q, kv):
    return kv * KV_SCALE * abs(q) * q


def drive_of(name):
    n, kc, tau0, kv, q0, pf, pt, dz = PIPES[name]
    return (pf - pt) - RHO * GRAV * dz


def steady_root(name):
    """Rate at which the two losses together ask for `drive`; bisection in
    wall-stress space (dp_wall and dp_valve are both monotone in Q and Q is
    monotone in tw). Below the threshold there is no positive root."""
    n, kc, tau0, kv, q0, pf, pt, dz = PIPES[name]
    drive = drive_of(name)
    if abs(drive) <= threshold(tau0):
        return 0.0
    sign = 1.0 if drive > 0 else -1.0
    ask = lambda tw: (WALL_SCALE * tw
                      + dp_valve(sign * q_hb(tw, n, kc, tau0), kv))
    lo, hi = tau0, tau0 + 1.0
    while ask(hi) < abs(drive):
        hi = 2.0 * hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if ask(mid) > abs(drive):
            hi = mid
        else:
            lo = mid
        if hi - lo <= 1e-16 * max(abs(mid), 1e-30):
            break
    return sign * q_hb(0.5 * (lo + hi), n, kc, tau0)


def step_row(q_prev, name):
    """One implicit-Euler step, static branch included, stated from the
    residual side in wall-stress space -- the same rule the kernel applies
    (entry latch + clamp), reached without the Newton loop."""
    n, kc, tau0, kv, q0, pf, pt, dz = PIPES[name]
    drive = drive_of(name)
    flow = lambda tw: q_hb(tw, n, kc, tau0)
    residual = lambda tw: (I_DT * (flow(tw) - q_prev) + WALL_SCALE * tw
                           + dp_valve(flow(tw), kv) - drive)
    if residual(tau0) >= 0.0:
        return 0.0
    lo, hi = tau0, tau0 + 1.0
    while residual(hi) < 0.0:
        hi = 2.0 * hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if residual(mid) > 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo <= 1e-16 * max(abs(mid), 1e-30):
            break
    return flow(0.5 * (lo + hi))


# --- the two runs -----------------------------------------------------------

proc = run(TRANSIENT, CASE)
if proc.returncode != 0:
    raise SystemExit(f"transient solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
main = parse_transient(proc.stdout)
sproc = run(STEADY, CASE_STEADY)
if sproc.returncode != 0:
    raise SystemExit(f"steady solver failed rc={sproc.returncode}\n"
                     f"{sproc.stdout}\n{sproc.stderr}")
st = parse_steady(sproc.stdout)

# --- 1. the wall law and the threshold are reported, and the threshold ------
#     carries no n, no K_cons and no K
require(len(main["law"]) == len(PIPES),
        f"the app states the wall law for exactly the {len(PIPES)} pipes")
worst_law = max(rel(main["law"][n][0], PIPES[n][0])
                + rel(main["law"][n][1], PIPES[n][1]) for n in PIPES)
require(worst_law == 0.0,
        "every declared (flowIndex, consistency) is read back exactly -- "
        "the wall law this row reads is the one the dict names")
worst_thr = max(rel(main["thr"][n], threshold(PIPES[n][2])) for n in PIPES)
require(worst_thr < 1e-12,
        "and every threshold is `4 tau_0 L / d` read back, not a number the "
        "kernel could have taken from the index or the valve")
require(main["thr"]["RelHB05"] == main["thr"]["RelHB15"]
        and main["thr"]["StuckRef"] == main["thr"]["StuckK1e8"],
        "the SAME tau_0 prints the SAME threshold at n = 0.5 and n = 1.5, and "
        "five decades of valve K print one threshold -- neither the index nor "
        "the valve is in the latch")
require(all(len(main["q"][n]) == NS and len(main["hold"][n]) == NS
            for n in PIPES),
        f"every pipe reports Q and the latch in all {NS} states")

# --- 2. a held pipe is the same pipe at any valve K -------------------------
for a, b, c in (("StuckRef", "StuckK500", "StuckK1e8"),):
    require(main["raw"][a] == main["raw"][b] == main["raw"][c],
            f"{a} / {b} / {c}: identical token lists over {NS} states at "
            f"K = 0 / 500 / 1e8 -- dp_valve(0) = 0, so a held HB wall has no "
            "valve drop to report")
require(main["raw"]["HangHB05"] == main["raw"]["HangHB05K500"]
        and set(main["q"]["HangHB05"]) == {0.0},
        "HangHB05 / HangHB05K500: identical and exactly zero -- the gravity "
        f"drive {drive_of('HangHB05'):.2f} Pa is under the 2000 Pa threshold "
        "whatever K says")
require(all(set(main["q"][n]) == {0.0} for n in HELD),
        "every held pipe prints `Q 0` in all states: held by the wall, not "
        "by the solver")

# --- 3. every reported state of every flowing pipe vs the recursion ---------
worst = 0.0
worst_at = ""
for n in FLOWING + STOPPING:
    q_prev = PIPES[n][4]
    for i in range(1, NS):
        q_model = step_row(q_prev, n)
        err = rel(q_model, main["q"][n][i])
        if err > worst:
            worst, worst_at = err, f"{n} step {i}"
        q_prev = main["q"][n][i]
require(worst < 5e-9,
        f"all {NS} states of all {len(FLOWING) + len(STOPPING)} flowing or "
        f"stopping pipes reproduce the implicit-Euler row whose losses are the "
        f"HB profile integral and K rho|Q|Q/2A^2 (worst rel {worst:.3e} at "
        f"{worst_at})")

# --- 4. the terminal rate is the root of the two-term loss; monotone in K ---
for n in FLOWING:
    require(rel(main["q"][n][-1], steady_root(n)) < 1e-9,
            f"{n}: terminal rate {main['q'][n][-1]:.12g} is the two-term root "
            f"{steady_root(n):.12g} of drive = dp_wall_HB + K rho Q^2/(2A^2)")
k_seq = [main["q"][n][-1] for n in
         ("RelHB05", "RelHB05K100", "RelHB05K500", "RelHB05K5000")]
require(all(a > b for a, b in zip(k_seq, k_seq[1:])),
        "one HB wall (n = 0.5, tau_0 = 1.0) with four valve K's is strictly "
        "monotone in K, and each entry lands on its OWN root")

# --- 5. the index moves the valve's share, not only the wall's rate ---------
require(main["q"]["RelHB15"][-1] > 10.0 * main["q"]["RelHB05"][-1],
        "at the same tau_0 the shear-thickening wall (n = 1.5, K = 0.02) flows "
        f"{main['q']['RelHB15'][-1] / main['q']['RelHB05'][-1]:.1f}x the "
        "shear-thinning one -- the index is in the law the valve is charged "
        "against")
share05 = dp_valve(main["q"]["RelHB05K500"][-1], 500.0)
share15 = dp_valve(main["q"]["RelHB15K500"][-1], 500.0)
require(share15 > 100.0 * share05,
        f"the same K = 500 takes {share05:.2f} Pa of the n = 0.5 pipe's drop "
        f"but {share15:.2f} Pa of the n = 1.5 pipe's ({share15 / share05:.0f}x) "
        "-- the valve's share is set by the wall law's index, which is the "
        "cell 5.50 could not ask with a Bingham wall")

# --- 6. the stop: three paths, three stop steps -----------------------------
stops = {}
for n in STOPPING:
    q = main["q"][n]
    idx = next((i for i in range(1, NS) if q[i] == 0.0), None)
    stops[n] = idx
require(all(isinstance(stops[n], int) for n in STOPPING),
        f"all three Stop pipes stop in finite time, at steps {stops}")
require(stops["StopK1e8"] < stops["StopRef"],
        f"StopK1e8 stops on step {stops['StopK1e8']}, StopRef on "
        f"{stops['StopRef']} -- the huge valve's own drop at 2.0e-4 m3/s is "
        "1.72e10 Pa, four million times the 4000 Pa on offer")
require(main["hold"]["StopRef"][stops["StopRef"]] == 0
        and all(main["hold"]["StopRef"][i] == 1
                for i in range(stops["StopRef"] + 1, NS)),
        "the latch picks StopRef up the step after the stop and never lets go "
        "-- the flag is the entry latch, not the rate")
q_star = (threshold(1.5) - drive_of("StopRef")) / I_DT
margins = {n: dp_valve(q_star, PIPES[n][3]) for n in ("StopRef", "StopK5000")}
require(margins["StopK5000"] < 0.2 * (threshold(1.5) - drive_of("StopRef")),
        f"StopK5000's valve drop at the q* = (dp_static - drive) dt A / (rho L) "
        f"margin is {margins['StopK5000']:.2e} Pa against the 8000 Pa the stop "
        "runs on, which is why it stops when StopRef stops")

# --- 7. the steady kernel: same root, printed split adds up -----------------
for n in FLOWING:
    n_idx, kc, tau0, kv, q0, pf, pt, dz = PIPES[n]
    require(rel(st[n]["q"], steady_root(n)) < 1e-9,
            f"{n}: steady kernel Q = {st[n]['q']:.12g} is the same two-term "
            "root the marched row ends on")
    # the printed `dp` is the pressure drop only; gravity is a separate term,
    # so the wall half of the loss is dp - rho g dz - dpValve
    dp_fric = st[n]["dp"] - RHO * GRAV * dz
    if kv > 0.0:
        require(rel(st[n]["dpValve"], dp_valve(st[n]["q"], kv)) < 5e-9,
                f"{n}: the printed dpValve is the valve's own formula at the "
                "printed rate (1e-11 print truncation on Q)")
        tw = (dp_fric - st[n]["dpValve"]) / WALL_SCALE
    else:
        require("dpValve" not in st[n],
                f"{n}: no dpValve field without a valve -- the steady kernel "
                "does not invent a local loss")
        tw = dp_fric / WALL_SCALE
    require(rel(q_hb(tw, n_idx, kc, tau0), st[n]["q"]) < 1e-9,
            f"{n}: the wall half of the printed drop, run back through 5.47's "
            "inversion, returns the printed rate")
for n in HELD + STOPPING:
    require(st[n]["q"] == 0.0 and "dpValve" not in st[n],
            f"{n}: steady kernel prints Q 0, the asked drop, and NO dpValve "
            "field -- nothing to charge a valve against at rest")

print()
print(f"{len(failures)} failure(s)")
if failures:
    sys.exit(1)
print("OK: the HB wall and the bare valve K share one momentum row with no "
      "schedule in it, and the valve's share is set by the wall law's index")
