#!/usr/bin/env python3
"""5.61: a source schedule waking a parked HB wall.

The new kernel hook is `sources{ <name> { node <N>; schedule { openRate; close
{ start; end; } open { start; end; } } } }`. Same ramp shape as valveClose /
valveOpen on a pipe, evaluated at the same implicit new time level t + dt. The
node balance picks up one more `+ Q_inj` term in sumIn (transient.cpp:400) and
the liquid inventory picks up `+ dt * Q_inj * fraction` in commitStep
(:480). Both are continuous in t.

The threshold, entry desk and in-step clamp are all schedule-blind
(transient.cpp:28-:31, :40-:43, :353-:355), so the release pressure
pA - 4 tau_0 L/d = 22000 Pa is the same in both files and the first release
step (318) is shared. The schedule can only change the STEP the wall is asked
to re-latch and re-release on, never whether.

Topology: one HB yield wall (Sa, tau_0 = 2.0, n = 0.5, K = 1.0) on a free
node F, one drain Da with the usual valve schedule, one source `inj` that
ramps 0 -> 2e-6 m3/s over [1.2, 1.4] s (steps 600..700), holds, and ramps
back to 0 over [2.0, 2.2] s (steps 1000..1100). Two files:

  * case_off.dict -- no `sources` block at all, the control.
  * case_on.dict  -- with the schedule.

Measured at case_off, 2001 states:
  * Sa flag transitions (events() output, indices where flag differs from
    i-1): [1, 318] = t=0 boundary + first release. The wall flows from step
    318 to the end of the run; no source means no reason to re-latch.
  * case_on transitions: [1, 318, 531, 1466] = t=0 boundary + release +
    source-driven re-latch (node pushed past 22000 during the hold) +
    source-off re-release (node falls back past 22000 after the ramp).
  * identity before the source: first 600 states bit-for-bit on p_F, Q_Sa,
    Q_Da, flag_Sa.

No kernel change from 5.59 / 5.60 other than the new parser fields and the
per-step re-sum in transient.cpp:329-:337; the Newton residual, node balance
and liquid inventory are the same three places the constant-rate source has
always fed.

Verified by this checker; ctest `verify.transientSourceYieldWake`.
"""
import math
import os
import re
import subprocess
import sys

TRANSIENT = os.path.abspath(sys.argv[1])
CASE_DIR = os.path.dirname(os.path.abspath(sys.argv[2]))

P_A, P_B = 3.0e4, 2.0e4
L, D = 20.0, 0.02
RHO, MU = 850.0, 0.1
DT = 0.002
V, BC = 0.1, 1.0e8
CC = V / BC
I_DT = RHO * L / (math.pi * (D / 2) ** 2) / DT
WALL_SCALE = 4 * L / D
THRESHOLD = 8000.0  # 4 * tau_0 * L / d, tau_0 = 2.0
RELEASE = P_A - THRESHOLD
TAU0, N_HB, K_HB = 2.0, 0.5, 1.0
VK = 0.0  # drain's open-K
KC = 1.0e8
CS, CE = 0.8, 1.0
OS_, OE = 2.5, 2.7
# source schedule
SRC_OPEN_RATE = 2.0e-6
SRC_CS, SRC_CE = 1.2, 1.4
SRC_OS, SRC_OE = 2.0, 2.2
FIRST_SEALED = int(round(CS / DT))  # 400
FIRST_REOPENED = int(round(OS_ / DT))  # 1250
N_STATES = int(round(4.0 / DT)) + 1  # 2001

failures = []
CHECKS = 0


def require(cond, what):
    global CHECKS
    CHECKS += 1
    print("ok:  " if cond else "FAIL:", what)
    if not cond:
        failures.append(what)


def run(case):
    out = subprocess.run([TRANSIENT, case], cwd=CASE_DIR, capture_output=True,
                         text=True)
    require(out.returncode == 0, "%s exited %d: %s"
            % (case, out.returncode, out.stderr.strip()[:300]))
    N, P, F, T = {}, {}, {}, []
    for line in out.stdout.splitlines():
        m = re.match(r"\[chflow\] t (\S+)", line)
        if m:
            T.append(float(m.group(1)))
            continue
        m = re.match(r"\[chflow\] node (\S+) p (\S+) ", line)
        if m:
            N.setdefault(m.group(1), []).append(float(m.group(2)))
            continue
        m = re.match(r"\[chflow\] pipe (\S+) Q (\S+)", line)
        if m:
            P.setdefault(m.group(1), []).append(float(m.group(2)))
            continue
        m = re.match(r"\[chflow\] pipe (\S+) static hold (\d)", line)
        if m:
            F.setdefault(m.group(1), []).append(int(m.group(2)))
    return T, N, P, F


def valve_k(t):
    if t < CS: return 0.0
    if t < CE: return (t - CS) / (CE - CS) * KC
    if t < OS_: return KC
    if t < OE: return KC + (t - OS_) / (OE - OS_) * (0.0 - KC)
    return 0.0


def source_q(t):
    if t < SRC_CS: return 0.0
    if t < SRC_CE: return (t - SRC_CS) / (SRC_CE - SRC_CS) * SRC_OPEN_RATE
    if t < SRC_OS: return SRC_OPEN_RATE
    if t < SRC_OE: return SRC_OPEN_RATE + (t - SRC_OS) / (SRC_OE - SRC_OS) * (0.0 - SRC_OPEN_RATE)
    return 0.0


def q_hb(tw):
    if tw <= TAU0: return 0.0
    d = tw - TAU0
    p = 1.0 / N_HB
    return math.pi * (D / 2) ** 3 * tw ** -3 * K_HB ** -p * (
        d ** (p + 3.0) / (p + 3.0)
        + 2.0 * TAU0 * d ** (p + 2.0) / (p + 2.0)
        + TAU0 * TAU0 * d ** (p + 1.0) / (p + 1.0))


def events(flags):
    return [i for i in range(1, len(flags)) if flags[i] != flags[i - 1]]


T_off, N_off, P_off, F_off = run("case_off.dict")
T_on,  N_on,  P_on,  F_on  = run("case_on.dict")

# --------------------------------------------------------------------------
# Family 1: identity before the source turns on (the two cases share a prefix
# of length int(SRC_CS / DT) = 600 states on every token).
# --------------------------------------------------------------------------
PREFIX = int(round(SRC_CS / DT))
for i in range(PREFIX):
    require(N_off["F"][i] == N_on["F"][i],
            "F1 p_F identical pre-source i=%d" % i)
    require(P_off["Sa"][i] == P_on["Sa"][i],
            "F1 Q_Sa identical pre-source i=%d" % i)
    require(P_off["Da"][i] == P_on["Da"][i],
            "F1 Q_Da identical pre-source i=%d" % i)
    require(F_off["Sa"][i] == F_on["Sa"][i],
            "F1 flag_Sa identical pre-source i=%d" % i)

# --------------------------------------------------------------------------
# Family 2: t = 0 flag is 0 in both files (no decision yet).
# --------------------------------------------------------------------------
require(F_off["Sa"][0] == 0, "F2 t=0 flag 0 off")
require(F_on["Sa"][0] == 0,  "F2 t=0 flag 0 on")

# --------------------------------------------------------------------------
# Family 3: release step identical in both files (schedule blind at release).
# --------------------------------------------------------------------------
ev_off = events(F_off["Sa"])
ev_on  = events(F_on["Sa"])
# events()[0] is the t=0 boundary (flag 0 at i=0 -> 1 at i=1); skip it.
require(ev_off[0] == 1, "F3 events[0] is the t=0 boundary (off)")
require(ev_on[0]  == 1, "F3 events[0] is the t=0 boundary (on)")
require(ev_off[1] == ev_on[1],
        "F3 real-release step shared: %d == %d" % (ev_off[1], ev_on[1]))
release_step = ev_off[1]
# drive the kernel uses at step `release_step` is pA - pPrev where pPrev is
# the committed pressure from state release_step - 1. So:
drive_deciding_release = P_A - N_off["F"][release_step - 1]
drive_one_earlier      = P_A - N_off["F"][release_step - 2]
require(drive_deciding_release > THRESHOLD,
        "F3 drive deciding release %.4f Pa > threshold"
        % drive_deciding_release)
require(drive_one_earlier <= THRESHOLD,
        "F3 drive one step earlier %.4f Pa <= threshold" % drive_one_earlier)

# --------------------------------------------------------------------------
# Family 4: the HB wall's rate is the section integral of its own implied tw
# at every unheld step, and is exactly 0 at every held step.
# --------------------------------------------------------------------------
for i in range(1, N_STATES):
    q_sa = P_off["Sa"][i]
    q_prev = P_off["Sa"][i - 1]
    drive = P_A - N_off["F"][i]
    acc = I_DT * (q_sa - q_prev)
    tw = (drive - acc) / WALL_SCALE
    if F_off["Sa"][i] == 1:
        require(q_sa == 0.0, "F4 Sa held => Q 0 at i=%d" % i)
    else:
        if tw > TAU0:
            q_pred = q_hb(tw)
            rel = abs(q_sa - q_pred) / max(abs(q_pred), 1e-30)
            require(rel < 1e-6,
                    "F4 Sa Q matches HB integral i=%d rel=%.3e" % (i, rel))

# --------------------------------------------------------------------------
# Family 5: node balance closes on the printed tokens alone.
#   V/BC * (p_F[i] - p_F[i-1]) / dt = Q_Sa[i] + Q_inj(t_i) - Q_Da[i]
# where Q_inj is 0 in case_off and the source schedule in case_on.
# --------------------------------------------------------------------------
worst_off = 0.0
for i in range(1, N_STATES):
    dp_obs = N_off["F"][i] - N_off["F"][i - 1]
    dp_pred = (P_off["Sa"][i] + 0.0 - P_off["Da"][i]) * DT / CC
    worst_off = max(worst_off, abs(dp_obs - dp_pred))
require(worst_off < 1e-6,
        "F5 node closure off worst residual %.3e Pa" % worst_off)

worst_on = 0.0
for i in range(1, N_STATES):
    dp_obs = N_on["F"][i] - N_on["F"][i - 1]
    src = source_q(T_on[i])
    dp_pred = (P_on["Sa"][i] + src - P_on["Da"][i]) * DT / CC
    worst_on = max(worst_on, abs(dp_obs - dp_pred))
require(worst_on < 1e-6,
        "F5 node closure on worst residual %.3e Pa" % worst_on)

# --------------------------------------------------------------------------
# Family 6: the schedule buys two new events -- a re-latch while the source
# holds and a re-release when it ramps off. case_off stays at 3 events
# (including the t=0 boundary), case_on at 5.
# --------------------------------------------------------------------------
require(len(ev_off) == 2, "F6 case_off has exactly 2 transitions (got %d)"
        % len(ev_off))
require(len(ev_on) == 4,  "F6 case_on  has exactly 4 transitions (got %d)"
        % len(ev_on))
# events[0] is the t=0 boundary; the three physical events in case_on are
# events[1..3]: release / source-latch / source-wake.
relatch_step = ev_on[2]
rerelease_step = ev_on[3]
t_relatch = relatch_step * DT
t_rerelease = rerelease_step * DT
require(SRC_CS < t_relatch < SRC_OE + 0.5,
        "F6 relatch step %d (t=%.3f) after source turns on"
        % (relatch_step, t_relatch))
require(SRC_OS - 0.1 < t_rerelease < SRC_OE + 1.5,
        "F6 re-release step %d (t=%.3f) after source ramp-off decays node"
        % (rerelease_step, t_rerelease))
# flags around the three physical events (skipping the t=0 boundary):
#   event 1 = release       : held (1) -> flowing (0)
#   event 2 = source-latch  : flowing (0) -> held (1)
#   event 3 = source-wake   : held (1) -> flowing (0)
require(F_on["Sa"][ev_on[1] - 1] == 1 and F_on["Sa"][ev_on[1]] == 0,
        "F6 event 1: held -> flowing")
require(F_on["Sa"][ev_on[2] - 1] == 0 and F_on["Sa"][ev_on[2]] == 1,
        "F6 event 2: flowing -> held (source pushes past release)")
require(F_on["Sa"][ev_on[3] - 1] == 1 and F_on["Sa"][ev_on[3]] == 0,
        "F6 event 3: held -> flowing (source ramps off)")

# --------------------------------------------------------------------------
# Family 7: the drain's valve K is schedule-independent -- case_on's drain
# sees the same K as case_off's at every state, so Da's rate only moves
# because F's pressure moves.
# --------------------------------------------------------------------------
for i in range(N_STATES):
    k_off = valve_k(T_off[i])
    k_on  = valve_k(T_on[i])
    require(k_off == k_on, "F7 valve K identical at i=%d" % i)

# --------------------------------------------------------------------------
# Family 8: the threshold is schedule-blind -- case_off and case_on both
# quote 8000 Pa on the startup line, and the release pressure pA - 8000 =
# 22000 Pa is the same integer in both files.
# --------------------------------------------------------------------------
for case_path in ["case_off.dict", "case_on.dict"]:
    out = subprocess.run([TRANSIENT, case_path], cwd=CASE_DIR,
                         capture_output=True, text=True)
    m = re.search(r"\[chflow\] pipe Sa static-hold threshold (\S+)",
                  out.stdout)
    require(m is not None, "F8 threshold line present in %s" % case_path)
    if m:
        require(abs(float(m.group(1)) - THRESHOLD) < 1e-9,
                "F8 threshold %.1f Pa in %s" % (float(m.group(1)), case_path))

print()
print("%d checks, %d failures" % (CHECKS, len(failures)))
sys.exit(0 if not failures else 1)
