#!/usr/bin/env python3
"""5.62 band edge: pB lowered from 20000 to 19500 so arm d (n=0.5, tau_0=2.5)
wakes for the first time. Zero kernel change -- the same holdsStatic, in-step
clamp, and degree-generic node row as 5.49/5.54/5.60. The only move is the
sink pressure, which lets the node drop below the release floor pA - 4 tau_0
L/d = 20000 Pa that was asymptotically unreachable at pB = 20000.

verified by tests/check_transient_yield_tau_index_band_edge.py
(ctest verify.transientYieldTauIndexBandEdge).
"""
import math
import os
import re
import subprocess
import sys

TRANSIENT = sys.argv[1]
CASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "cases", "transientYieldTauIndexBandEdge")

P_A = 3.0e4
P_B = 1.95e4
DT = 0.002
RHO = 850
MU = 0.1
L = 20
D = 0.02
TAU0_WEAK = 1.5
THRESHOLD_WEAK = 4 * TAU0_WEAK * L / D  # 6000 Pa

CASES = [
    ("5", 0.5, 1.0),
    ("1", 1.0, 0.1),
    ("2", 2.0, 0.02),
]

ARMS = ["a", "b", "c", "d"]
TAU0_STRONG = {"a": 1.75, "b": 2.0, "c": 2.25, "d": 2.5}
THRESHOLD = {a: 4 * t * L / D for a, t in TAU0_STRONG.items()}
RELEASE = {a: P_A - th for a, th in THRESHOLD.items()}

failures = []
CHECKS = 0


def require(cond, what):
    global CHECKS
    CHECKS += 1
    print("ok:  " if cond else "FAIL:", what)
    if not cond:
        failures.append(what)


def run(case):
    out = subprocess.run([TRANSIENT, os.path.join(CASE_DIR, case)],
                         capture_output=True, text=True)
    require(out.returncode == 0, "%s exited %d: %s"
            % (case, out.returncode, out.stderr.strip()[:300]))
    N, P, F = {}, {}, {}
    for line in out.stdout.splitlines():
        m = re.match(r"\[chflow\] node (\S+) p (\S+) ", line)
        if m:
            N.setdefault(m.group(1), []).append(float(m.group(2)))
            continue
        m = re.match(r"\[chflow\] pipe (\S+) Q (\S+)", line)
        if m:
            P.setdefault(m.group(1), []).append(m.group(2))
            continue
        m = re.match(r"\[chflow\] pipe (\S+) static hold (\d)", line)
        if m:
            F.setdefault(m.group(1), []).append(int(m.group(2)))
    return N, P, F


def events(flags):
    return [i for i in range(len(flags)) if i == 0 or flags[i] != flags[i - 1]]


def drive_of(pn, i):
    return P_A - pn[i]


R = {}
for tag, n, K in CASES:
    R[tag] = run("case_%s.dict" % tag)

# --- family 1: arm d wakes in n=0.5, stays NEVER in n=1 and n=2 ----------
for tag, n, K in CASES:
    N, P, F = R[tag]
    ns = len(N["F" + "a"])
    require(ns == 2001, "%s: 2001 states, got %d" % (tag, ns))

    ev_d = events(F["S" + "d"])
    flowing_d = sum(1 for f in F["S" + "d"] if f == 0)
    if tag == "5":
        require(flowing_d > 100, "%s/d: arm d wakes (was NEVER), %d flowing"
                % (tag, flowing_d))
        require(len(ev_d) >= 2, "%s/d: at least release event" % tag)
    else:
        require(flowing_d <= 1, "%s/d: arm d stays NEVER (only t=0 boundary)"
                % tag)

# --- family 2: release drive above threshold at every first-release step ---
for tag, n, K in CASES:
    N, P, F = R[tag]
    for a in ARMS:
        flags = F["S" + a]
        rel_step = None
        for i in range(1, len(flags)):
            if flags[i - 1] == 1 and flags[i] == 0:
                rel_step = i
                break
        if rel_step is None:
            continue
        drive = drive_of(N["F" + a], rel_step - 1)
        th = THRESHOLD[a]
        require(drive > th,
                "%s/%s: drive at release %.4f Pa > threshold %.0f"
                % (tag, a, drive, th))
        if rel_step >= 2:
            drive_prev = drive_of(N["F" + a], rel_step - 2)
            require(drive_prev <= th,
                    "%s/%s: drive one step earlier %.4f Pa <= threshold %.0f"
                    % (tag, a, drive_prev, th))

# --- family 3: re-latch drive below threshold ----------------------------
for tag, n, K in CASES:
    N, P, F = R[tag]
    for a in ARMS:
        ev = events(F["S" + a])
        for j, step in enumerate(ev):
            if step == 0 or step == ev[0]:
                continue
            was_held = F["S" + a][step - 1] == 1
            now_free = F["S" + a][step] == 0
            was_free = F["S" + a][step - 1] == 0
            now_held = F["S" + a][step] == 1
            if was_free and now_held:
                drive = drive_of(N["F" + a], step - 1)
                th = THRESHOLD[a]
                require(drive <= th,
                        "%s/%s: re-latch at %d, drive %.4f <= %.0f"
                        % (tag, a, step, drive, th))

# --- family 4: t=0 flags all zero ----------------------------------------
for tag, n, K in CASES:
    N, P, F = R[tag]
    for a in ARMS:
        for w in ["S", "W"]:
            require(F[w + a][0] == 0,
                    "%s/%s%s: t=0 flag is 0" % (tag, w, a))

# --- family 5: node balance on printed tokens ----------------------------
V = 0.1
B = 1.0e8

for tag, n, K in CASES:
    N, P, F = R[tag]
    for a in ARMS:
        pF = N["F" + a]
        qS = [float(x) for x in P["S" + a]]
        qW = [float(x) for x in P["W" + a]]
        qD = [float(x) for x in P["D" + a]]
        fS = F["S" + a]
        fW = F["W" + a]
        worst = 0
        for i in range(1, len(pF)):
            qS_eff = 0.0 if fS[i] else qS[i]
            qW_eff = 0.0 if fW[i] else qW[i]
            dpdt = V / B * (pF[i] - pF[i - 1]) / DT
            residual = abs(qS_eff + qW_eff - qD[i] - dpdt)
            if residual > worst:
                worst = residual
        require(worst < 1e-13,
                "%s/%s: node balance worst %.4e m3/s"
                % (tag, a, worst))

# --- family 6: arm d release step is late (after drain reopen) -----------
for tag, n, K in CASES:
    N, P, F = R[tag]
    flags_d = F["S" + "d"]
    rel_step_d = None
    for i in range(1, len(flags_d)):
        if flags_d[i - 1] == 1 and flags_d[i] == 0:
            rel_step_d = i
            break
    if rel_step_d is not None:
        require(rel_step_d > 1350,
                "%s/d: arm d releases after drain reopen (step %d > 1350)"
                % (tag, rel_step_d))

# --- family 7: arm c in n=0.5 has re-latch cycle -------------------------
N5, P5, F5 = R["5"]
ev_c = events(F5["Sc"])
require(len(ev_c) >= 4, "5/c: arm c has re-latch + re-release (%d events)"
        % len(ev_c))

# --- family 8: weak walls release at same step across arms ---------------
for tag, n, K in CASES:
    N, P, F = R[tag]
    weak_rel = []
    for a in ARMS:
        ev = events(F["W" + a])
        if len(ev) >= 2:
            weak_rel.append(ev[1])
    require(len(set(weak_rel)) == 1,
            "%s: all weak walls release at same step %s"
            % (tag, weak_rel))

# --- summary -------------------------------------------------------------
print("\n%d checks, %d failures" % (CHECKS, len(failures)))
if failures:
    for f in failures:
        print("  FAIL:", f)
    sys.exit(1)
