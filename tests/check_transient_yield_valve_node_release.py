#!/usr/bin/env python3
"""transientYieldValveNodeRelease: the yield wall's static branch (5.49), a valve
timetable (M2's `valveClose` / `valveOpen`, priced against the wall in 5.50 and
against the schedule in 5.52) and a node whose pressure the solve has to find --
the three-way cell 5.50 put in its boundary (4) and 5.52 in its own (2).

5.52's answer was "a timetable reaches neither branch", and its reason was
arithmetic: with both ends pinned the drive is a number in the case file,
`drivePa` (transient.cpp:52-:56) contains no K, and the valve term
`K rho |Q| Q / (2 A^2)` vanishes with the rate. This file changes the desktop
and not the kernel. One end of every held pipe is a free node, so the drive the
wall is compared against is whatever the network says, and the row per pair is

    rho L/A (Qh - Qh')/dt = (pA - pF) - 4 L tw_h / d - Kh rho|Qh|Qh/(2A^2)
    rho L/A (Qb - Qb')/dt = (pF - pB) - 128 mu L Qb/(pi d^4) - Kb(t) rho|Qb|Qb/(2A^2)
    V/BC (pF - pF')        = Qh - Qb,        Kb(t) = valveKAt(t + dt)

with Qh entering that last row as the forced zero while the wall holds. Two
things follow, and both are measured here rather than argued:

  * whose timetable owns the release. On the neighbour it owns it completely --
    the step goes 253 / 299 / 1830 for K = 0 / 5e3 / 1e6, one 1e8 valve and one
    1e8 ramp push it past the end of the run, and a neighbour that *reopens*
    wakes the pipe at step 1290, which is 5.52's family 3 turned around. On the
    held pipe's own row it owns nothing: the same frozen 1e8 and the same
    close/reopen cycle release on `K0`'s step with the same node history.
  * the wall answers one number. Every release in this file, across three decades
    of K and two schedules, happens at a node pressure within a few Pa of
    pA - 4 tau_0 L / d = 22000 Pa. K decides when the network gets there.

Nine families of judgement:
  1. the threshold is the same 8000 Pa on all nine yielded pipes and on none of
     the nine Newtonian bleeds, and the release always happens *at* the pressure
     the threshold implies, whatever the neighbour's K
  2. a frozen valve on the neighbour moves the release step, monotonically in K,
     and 1e8 moves it out of the run
  3. one closedK, one ramp length, two arrival times: a ramp that lands before
     the drain gets there cancels the release, one that lands after it is
     invisible until the sampled `t + dt` steps inside it
  4. a neighbour that is sealed and then reopened wakes what the wall holds,
     and is bit-for-bit the never-released pair up to the step its own schedule
     stops being closed
  5. the same two prescriptions on the held pipe's own row: identical tokens and
     identical node history while held, same release step, and the trace parts
     on the first moving state
  6. the release is one-way, and measurably so: nothing re-latches, and the
     network's own settled pressure sits *below* the release pressure
  7. every state of every pipe and every node, in both the 4 s and the 30 s run,
     against the checker's own copy of the two rows, the node row and the
     sampled K
  8. the node balance closes to the solver's own exit residual, i.e. a held pipe
     enters continuity as the forced zero it prints, not as a leak
  9. the calendar against the physics: the steady kernel's nine rates are this
     file's own bisection and close the 1e4 Pa across the two series pipes; the
     three pairs whose valves end open settle to one rate and the four carrying
     one 1e8 valve to another, while the floating node of those two placements
     splits by exactly that valve's drop; the tails meet the settled state at
     tolerances their own remaining clock earns, two pairs are named as
     mid-flight, and the pairs the 4 s window never released trickle in the
     steady kernel -- with the 30 s run releasing one of them at t = 24.25 s and
     still 304 Pa short on the other at 30 s
"""

import math
import subprocess
import sys

TRANSIENT = sys.argv[1]
STEADY = sys.argv[2]
CASE = sys.argv[3]
CASE_DIR = CASE.rsplit("/", 1)[0]

RHO = 850.0
MU = 0.1
DT = 0.002
END = 4.0
NS = int(round(END / DT)) + 1
LONG_END = 30.0
LONG_RPT = 5
NL = int(round(LONG_END / (DT * LONG_RPT))) + 1
P_A, P_B = 3.0e4, 2.0e4
LENGTH = 20.0
DIAM = 0.02
AREA = 0.25 * math.pi * DIAM ** 2
INERTIA = RHO * LENGTH / AREA
I_DT = INERTIA / DT
KV_SCALE = RHO / (2.0 * AREA * AREA)
WALL_SCALE = 4.0 * LENGTH / DIAM          # dp = WALL_SCALE * tw
R_LAM = 128.0 * MU * LENGTH / (math.pi * DIAM ** 4)
NODE_VOLUME = 0.1
NODE_BULK = 1.0e8
CC = NODE_VOLUME / NODE_BULK               # 1e-9 m3/Pa
P_NODE0 = 2.5e4
TAU0 = 2.0
THRESHOLD = WALL_SCALE * TAU0              # 8000 Pa
P_RELEASE = P_A - THRESHOLD                # 22000 Pa
RATE_FLOOR = 1.0e-8
NEWTON_TOL = 1.0e-11

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


# --- the table, read off case.dict ------------------------------------------
# pair: (node, held valve K, held schedule, bleed valve K, bleed schedule)
# a schedule is (closeStart, closeEnd, closedK, openStart, openEnd); None is a
# frozen valve.
SCHED_W = (0.0, 0.002, 1.0e8, 2.0, 2.1)
PAIRS = {
    "K0": ("Fk0", 0.0, None, 0.0, None),
    "K3": ("Fk3", 0.0, None, 5.0e3, None),
    "K6": ("Fk6", 0.0, None, 1.0e6, None),
    "K8": ("Fk8", 0.0, None, 1.0e8, None),
    "R1": ("Fr1", 0.0, None, 0.0, (0.2, 0.3, 1.0e8, -1.0, -1.0)),
    "R2": ("Fr2", 0.0, None, 0.0, (1.0, 1.1, 1.0e8, -1.0, -1.0)),
    "W":  ("Fw", 0.0, None, 0.0, SCHED_W),
    "O":  ("Fo", 1.0e8, None, 0.0, None),
    "OW": ("Fow", 0.0, SCHED_W, 0.0, None),
}
HELD = {f"Held{s}": p for s, p in PAIRS.items()}
BLEED = {f"Bleed{s}": p for s, p in PAIRS.items()}
PIPES = {}
for s, (node, kh, sh, kb, sb) in PAIRS.items():
    PIPES[f"Held{s}"] = (TAU0, kh, sh, node)
    PIPES[f"Bleed{s}"] = (0.0, kb, sb, node)
LONG_PAIRS = ("K8", "R1")
RELEASED = None      # filled in from the run
NEVER = None


def valve_k(sched, k_open, t):
    """transient.cpp:125-:133 verbatim: linear ramps, closedK held between the
    ramps, `openStart < 0` keeps it shut."""
    if sched is None:
        return k_open
    cs, ce, kc, os_, oe = sched
    if t < cs:
        return k_open
    if t < ce:
        return k_open + (t - cs) / (ce - cs) * (kc - k_open)
    if os_ < 0 or t < os_:
        return kc
    if t < oe:
        return kc + (t - os_) / (oe - os_) * (k_open - kc)
    return k_open


def k_of(pipe, t):
    tau0, kop, sched, node = PIPES[pipe]
    return valve_k(sched, kop, t)


# --- the two rows, stated from this file's copy of the laws ------------------

def q_of_tw(tw, tau0):
    """The exact laminar discharge a yielded wall stress holds (5.46 at n = 1,
    5.47's profile integral), and with tau0 = 0 it is Hagen-Poiseuille."""
    if tw <= tau0:
        return 0.0
    d = tw - tau0
    r = 0.5 * DIAM
    return (math.pi * r ** 3 * tw ** -3 / MU * d * d
            * (tau0 * tau0 / 2.0 + 2.0 * tau0 * d / 3.0 + d * d / 4.0))


def row_yield(q_prev, drive, kv, tau0=TAU0):
    """One implicit row of a yielded pipe, static branch included, in wall-stress
    space where the residual is monotone. Returns (Q, tw)."""
    def g(tw):
        q = q_of_tw(tw, tau0)
        return I_DT * (q - q_prev) + WALL_SCALE * tw + kv * KV_SCALE * q * q - drive
    if g(tau0) >= 0.0:
        return 0.0, tau0
    lo, hi = tau0, tau0 + 1.0
    while g(hi) < 0.0:
        hi = 2.0 * hi
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo <= 1e-14 * max(mid, 1e-30):
            break
    tw = 0.5 * (lo + hi)
    return q_of_tw(tw, tau0), tw


def row_laminar(q_prev, drive, kv):
    """The same row for the Newtonian bleed, whose residual is a quadratic in Q
    with a known sign, so the root is written out instead of searched for."""
    b = I_DT + R_LAM
    c = kv * KV_SCALE
    rhs = I_DT * q_prev + drive
    if drive >= 0.0 and rhs >= 0.0:
        if c > 0.0:
            return (-b + math.sqrt(b * b + 4.0 * c * rhs)) / (2.0 * c)
        return rhs / b
    return -row_laminar(-q_prev, -drive, kv)


def march(pair, nsteps, report, dt=DT):
    """The whole trace of one pair: both rows, the node row, the sampled K."""
    node, kh, sh, kb, sb = PAIRS[pair]
    hp, bp = f"Held{pair}", f"Bleed{pair}"
    qh, qb, p = 0.0, 0.0, P_NODE0
    hs, bs, ps = [qh], [qb], [p]
    drives_pos = True
    for i in range(1, nsteps + 1):
        tn = i * dt                       # the implicit new time level
        kh_t = valve_k(sh, kh, tn)
        kb_t = valve_k(sb, kb, tn)
        held = (qh == 0.0) and abs(P_A - p) <= THRESHOLD
        pn = p
        for _ in range(3):                # the node row is a 1.5e-5 contraction
            if held:
                qn_h = 0.0
            else:
                qn_h, _ = row_yield(qh, P_A - pn, kh_t)
            qn_b = row_laminar(qb, pn - P_B, kb_t)
            drives_pos = drives_pos and (P_A - pn) >= 0.0 and (pn - P_B) >= 0.0
            pn = p + (dt / CC) * (qn_h - qn_b)
        if not held:
            qn_h, _ = row_yield(qh, P_A - pn, kh_t)
        qn_b = row_laminar(qb, pn - P_B, kb_t)
        qh, qb, p = qn_h, qn_b, pn
        if i % report == 0:
            hs.append(qh)
            bs.append(qb)
            ps.append(p)
    return hs, bs, ps, drives_pos


def parse_transient(stdout, pipes, nodes):
    out = {"q": {n: [] for n in pipes}, "raw": {n: [] for n in pipes},
           "hold": {n: [] for n in pipes}, "p": {n: [] for n in nodes},
           "thr": {}, "times": []}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 2 and tok[0] == "t":
            out["times"].append(float(tok[1]))
        elif len(tok) >= 4 and tok[0] == "node" and tok[2] == "p":
            if tok[1] in out["p"]:
                out["p"][tok[1]].append(float(tok[3]))
        elif len(tok) == 4 and tok[0] == "pipe" and tok[2] == "Q":
            if tok[1] in out["q"]:
                out["q"][tok[1]].append(float(tok[3]))
                out["raw"][tok[1]].append(tok[3])
        elif len(tok) == 5 and tok[0] == "pipe" and tok[3] == "hold":
            if tok[1] in out["hold"]:
                out["hold"][tok[1]].append(int(tok[4]))
        elif len(tok) == 5 and tok[0] == "pipe" and tok[3] == "threshold":
            out["thr"][tok[1]] = float(tok[4])
    return out


def parse_steady(stdout):
    out = {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) < 4 or tok[0] != "pipe" or tok[2] != "Q" or tok[1] not in PIPES:
            continue
        rec = {"q": float(tok[3])}
        for i, t in enumerate(tok):
            if t == "dp":
                rec["dp"] = float(tok[i + 1])
            elif t == "dpValve":
                rec["dpValve"] = float(tok[i + 1])
        out[tok[1]] = rec
    nodes = {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 4 and tok[0] == "node" and tok[2] == "p":
            nodes[tok[1]] = float(tok[3])
    return out, nodes


def release_step(flags):
    """The step whose *entry* drive is over the threshold, i.e. the first state
    the latch stops raising after it has been raising."""
    for i in range(2, len(flags)):
        if flags[i] == 0 and flags[i - 1] == 1:
            return i
    return None


# --- runs --------------------------------------------------------------------
proc = run(TRANSIENT, CASE)
if proc.returncode != 0:
    raise SystemExit(f"transient solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
main = parse_transient(proc.stdout, PIPES, [p[0] for p in PAIRS.values()])
proc = run(TRANSIENT, f"{CASE_DIR}/case_long.dict")
if proc.returncode != 0:
    raise SystemExit(f"long solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
longn = parse_transient(proc.stdout, [f"Held{s}" for s in LONG_PAIRS]
                        + [f"Bleed{s}" for s in LONG_PAIRS],
                        [PAIRS[s][0] for s in LONG_PAIRS])
proc = run(STEADY, f"{CASE_DIR}/case_steady.dict")
if proc.returncode != 0:
    raise SystemExit(f"steady solver failed rc={proc.returncode}\n"
                     f"{proc.stdout}\n{proc.stderr}")
steady, steady_p = parse_steady(proc.stdout)

rel_step = {s: release_step(main["hold"][f"Held{s}"]) for s in PAIRS}
RELEASED = tuple(s for s in PAIRS if rel_step[s] is not None)
NEVER = tuple(s for s in PAIRS if rel_step[s] is None)

require(len(main["times"]) == NS and len(longn["times"]) == NL,
        f"the 4 s run prints {NS} states and the 30 s run {NL} at "
        f"reportEvery {LONG_RPT}")
for n in PIPES:
    require(len(main["q"][n]) == NS, f"pipe {n}: {NS} Q states")
for nd in PAIRS.values():
    require(len(main["p"][nd[0]]) == NS, f"node {nd[0]}: {NS} pressures")

# --- 1. one threshold, one release pressure, whatever the K -----------------
require(set(main["thr"]) == set(HELD),
        f"a static-hold threshold is printed for exactly the {len(HELD)} yielded "
        f"pipes and for none of the {len(BLEED)} Newtonian bleeds "
        f"({sorted(main['thr'])})")
require(all(abs(v - THRESHOLD) < 1e-9 for v in main["thr"].values())
        and len(set(main["thr"].values())) == 1,
        f"and all {len(main['thr'])} of them are the same 4 tau_0 L / d = "
        f"{THRESHOLD:.0f} Pa, including the six whose neighbour carries a valve "
        f"from 5e3 to 1e8 -- the threshold never sees K")
p_at_rel = {s: main["p"][PAIRS[s][0]][rel_step[s] - 1] for s in RELEASED}
require(all(abs(v - P_RELEASE) < 12.0 for v in p_at_rel.values()),
        f"every release in the file happens at the node pressure its threshold "
        f"implies, pA - {THRESHOLD:.0f} = {P_RELEASE:.0f} Pa: "
        + ", ".join(f"{s} {p_at_rel[s]:.2f}" for s in RELEASED)
        + f" (spread {max(p_at_rel.values()) - min(p_at_rel.values()):.2f} Pa)")
require(all(P_A - main["p"][PAIRS[s][0]][rel_step[s] - 2] <= THRESHOLD <
            P_A - main["p"][PAIRS[s][0]][rel_step[s] - 1] for s in RELEASED),
        f"and the latch flips on the step whose entry drive crosses it: for all "
        f"{len(RELEASED)} released pipes the previous sample is at or below "
        f"{THRESHOLD:.0f} Pa and this one above")

# --- 2. the frozen ladder on the neighbour ----------------------------------
lad = (rel_step["K0"], rel_step["K3"], rel_step["K6"])
require(all(v is not None for v in lad) and lad[0] < lad[1] < lad[2],
        f"K = 0 / 5e3 / 1e6 on the bleed release on steps {lad[0]} / {lad[1]} / "
        f"{lad[2]} (t = {lad[0] * DT:.3f} / {lad[1] * DT:.3f} / "
        f"{lad[2] * DT:.3f} s) -- monotone in K, and one valve five decades "
        f"bigger takes the pipe out of the run entirely: {NEVER} print `0` in all "
        f"{NS} states")
for s in NEVER:
    fl = main["hold"][f"Held{s}"]
    require(sum(fl[1:]) == NS - 1 and set(main["raw"][f"Held{s}"]) == {"0"},
            f"Held{s}: held in all {NS - 1} marched steps and the single token "
            f"`0` in every printed state, while its neighbour keeps moving "
            f"(Bleed{s} {min(main['q'][f'Bleed{s}'][1:]):.3e} - "
            f"{max(main['q'][f'Bleed{s}']):.3e} m3/s)")
    require(max(abs(P_A - v) for v in main["p"][PAIRS[s][0]]) < THRESHOLD,
            f"{PAIRS[s][0]}: the node never reaches the release pressure either "
            f"(its deepest drive is "
            f"{max(abs(P_A - v) for v in main['p'][PAIRS[s][0]]):.4f} Pa against "
            f"{THRESHOLD:.0f} Pa), so the pipe is not held by a latch that "
            f"overslept")

# --- 3. arrival time, not size ----------------------------------------------
require(rel_step["R1"] is None and rel_step["R2"] == rel_step["K0"],
        f"one closedK of 1e8 and one 0.1 s ramp, started at 0.2 s and at 1.0 s: "
        f"the early ramp cancels the release (R1: {NEVER}) and the late one "
        f"releases on K0's own step {rel_step['R2']}")
i_ramp = next(i for i in range(1, NS)
              if main["raw"]["HeldR2"][i] != main["raw"]["HeldK0"][i])
require(main["raw"]["HeldR2"][:i_ramp] == main["raw"]["HeldK0"][:i_ramp]
        and i_ramp == int(round(1.0 / DT)) + 1,
        f"HeldR2 and HeldK0 print the same token in all {i_ramp} states up to "
        f"the ramp, and the first state that differs is {i_ramp} = t "
        f"{main['times'][i_ramp]:.3f} s, which is the first step whose sampled "
        f"t + dt = {i_ramp * DT:.3f} s lies inside the ramp")
require(main["p"]["Fr2"][:i_ramp] == main["p"]["Fk0"][:i_ramp],
        f"and its node agrees bit for bit over those {i_ramp} samples "
        f"({main['p']['Fr2'][i_ramp - 1]:.6f} Pa both), so the schedule is not "
        f"in the arithmetic until the row is being read")
require(main["q"]["HeldR2"][-1] < main["q"]["HeldK0"][-1]
        and main["p"]["Fr2"][-1] > main["p"]["Fk0"][-1],
        f"after it lands the ramp works against the drain: the rate falls "
        f"{main['q']['HeldR2'][-1]:.4e} vs K0's {main['q']['HeldK0'][-1]:.4e} "
        f"and the node stops higher ({main['p']['Fr2'][-1]:.1f} vs "
        f"{main['p']['Fk0'][-1]:.1f} Pa) -- a neighbour's valve is a discount on "
        f"the *release*, not on the wall")

# --- 4. a neighbour's reopening is what wakes the wall -----------------------
i_ramp_w = next(i for i in range(1, NS)
                if main["raw"]["BleedW"][i] != main["raw"]["BleedK8"][i])
require(main["raw"]["BleedW"][:i_ramp_w] == main["raw"]["BleedK8"][:i_ramp_w]
        and main["p"]["Fw"][:i_ramp_w] == main["p"]["Fk8"][:i_ramp_w]
        and i_ramp_w == int(round(2.0 / DT)) + 1,
        f"W's bleed and the never-released K8's bleed print the same token, and "
        f"their nodes the same pressure, in all {i_ramp_w} states while the "
        f"schedules agree; the first state that differs is {i_ramp_w} = t "
        f"{main['times'][i_ramp_w]:.3f} s, the first step whose sampled t + dt "
        f"lies inside the reopen ramp")
i_wake = release_step(main["hold"]["HeldW"])
require(set(main["raw"]["HeldW"][:i_wake]) == {"0"}
        and main["raw"]["HeldW"][:i_wake] == main["raw"]["HeldK8"][:i_wake],
        f"and the wall does not answer that first difference: HeldW keeps "
        f"printing `0` for another {i_wake - i_ramp_w} steps, the same tokens as "
        f"HeldK8, because the node has to travel before the drive reaches its "
        f"threshold")
require(i_wake == 1290 and rel_step["W"] == i_wake,
        f"it lets go on step {i_wake} (t = {i_wake * DT:.3f} s) -- "
        f"{(i_wake - i_ramp_w) * DT:.3f} s after the schedule moved its "
        f"neighbour -- which is 5.52's family 3 turned around: no schedule wakes "
        f"what the wall holds *on the held pipe's own row*, and this one does")
require(main["q"]["HeldOW"][:rel_step["OW"]] == main["q"]["HeldK0"][:rel_step["K0"]]
        and rel_step["OW"] == rel_step["K0"],
        "the identical close-and-reopen cycle written on the held pipe's own row "
        "does not wake or delay anything: it releases on K0's step "
        f"{rel_step['K0']}, having printed 0 in the same "
        f"{rel_step['K0']} states")

# --- 5. the held pipe's own valve is invisible until the wall lets go --------
i_own = next(i for i in range(1, NS)
             if main["raw"]["HeldO"][i] != main["raw"]["HeldK0"][i])
require(i_own == rel_step["O"] == rel_step["K0"]
        and main["p"]["Fo"][:i_own] == main["p"]["Fk0"][:i_own],
        f"HeldO carries a frozen K = 1e8 on its own row and is the same pipe as "
        f"HeldK0 in all {i_own} held states and in its node history too; the "
        f"first state that differs is the release itself, state {i_own}")
require(0.0 < main["q"]["HeldO"][i_own] < main["q"]["HeldK0"][i_own],
        f"and the released pipe does pay its own valve at once: "
        f"{main['q']['HeldO'][i_own]:.6e} vs "
        f"{main['q']['HeldK0'][i_own]:.6e} m3/s on that first moving state -- "
        f"the K was never in the latch, it was waiting in the row")
tail_ow = rel(main["q"]["HeldOW"][-1], main["q"]["HeldK0"][-1])
require(tail_ow < 1e-4,
        f"the own-row cycle that could not move the release did cost rate while "
        f"it held the pipe shut, and comes back to the no-valve twin once "
        f"reopened: tail rel {tail_ow:.3e}")

# --- 6. one-way ---------------------------------------------------------------
for s in RELEASED:
    i = rel_step[s]
    require(sum(main["hold"][f"Held{s}"][i:]) == 0,
            f"Held{s}: no step after its release re-raises the latch")
sealed = ("R2", "W", "OW")
require(all(max(main["p"][PAIRS[s][0]][rel_step[s]:]) < P_RELEASE for s in sealed),
        "and the pairs whose neighbour sealed *after* the release cannot push the "
        "node back over it: " + ", ".join(
            f"{s} peaks at {max(main['p'][PAIRS[s][0]][rel_step[s]:]):.1f} Pa"
            for s in sealed) + f" against {P_RELEASE:.0f} Pa")
require(all(steady_p[PAIRS[s][0]] < P_RELEASE for s in sealed),
        f"the steady kernel says the same thing about the settled state: "
        + ", ".join(f"{s} {steady_p[PAIRS[s][0]]:.1f} Pa" for s in sealed)
        + f", all below the release pressure -- the only flow that could lift the "
        f"node past it is the yielded pipe's own, and that flow is zero exactly "
        f"at {P_RELEASE:.0f} Pa")

# --- 7. every state of every pipe and node, both runs ------------------------
worst_rel, worst_abs, worst_node, stop_pred, drives_ok = {}, {}, {}, {}, True
for s in PAIRS:
    hs, bs, ps, ok = march(s, NS - 1, 1)
    drives_ok = drives_ok and ok
    hp, bp = f"Held{s}", f"Bleed{s}"
    wr, wa = 0.0, 0.0
    for arr, obs in ((hs, main["q"][hp]), (bs, main["q"][bp])):
        for i in range(NS):
            a = abs(arr[i] - obs[i])
            wa = max(wa, a)
            scale = max(abs(arr[i]), abs(obs[i]))
            if scale >= RATE_FLOOR:
                wr = max(wr, a / scale)
    wn = max(rel(ps[i], main["p"][PAIRS[s][0]][i]) for i in range(NS))
    worst_rel[s], worst_abs[s], worst_node[s] = wr, wa, wn
    stop_pred[s] = release_step([1 if (q == 0.0 and h == 1) else 0
                                 for q, h in zip(hs, main["hold"][f"Held{s}"])])
require(drives_ok, "every drive the recursion saw stayed non-negative, so no "
                   "branch of either row was taken on a sign it cannot check")
for s in PAIRS:
    require(worst_rel[s] < 5e-9 and worst_node[s] < 5e-9,
            f"{s}: all {NS} states of both pipes and of the node reproduced by "
            f"the two rows, the node row and the sampled K (worst rate rel "
            f"{worst_rel[s]:.3e}, node rel {worst_node[s]:.3e})")
require(max(worst_abs.values()) < 1e-14,
        f"and the states the relative gate cannot see -- the zeros and the "
        f"one-step-away-from-zero handover -- agree absolutely: worst "
        f"{max(worst_abs.values()):.3e} m3/s")
for s in PAIRS:
    require(stop_pred[s] == rel_step[s],
            f"{s}: the recursion lets go on step {stop_pred[s]}, the app's latch "
            f"flips on {rel_step[s]}")
long_pred = {}
for s in LONG_PAIRS:
    hs, bs, ps, ok = march(s, int(LONG_END / DT), LONG_RPT)
    long_pred[s] = (hs, bs, ps)
    obs_h = longn["q"][f"Held{s}"]
    moving = [i for i in range(NL)
              if max(abs(hs[i]), abs(obs_h[i])) >= RATE_FLOOR]
    e_h = max((rel(hs[i], obs_h[i]) for i in moving), default=0.0)
    e_b = max(rel(bs[i], longn["q"][f"Bleed{s}"][i]) for i in range(NL))
    e_p = max(rel(ps[i], longn["p"][PAIRS[s][0]][i]) for i in range(NL))
    a_h = max(abs(hs[i] - obs_h[i]) for i in range(NL))
    require(e_h < 5e-9 and e_b < 5e-9 and e_p < 5e-9 and a_h < 1e-14,
            f"long {s}: all {NL} reported states of the 30 s run (15000 marched "
            f"steps) reproduced too -- rate rel {e_h:.3e} over the "
            f"{len(moving)} moving samples, bleed rel {e_b:.3e}, node rel "
            f"{e_p:.3e}, held-state abs {a_h:.3e}")

# --- 8. what the node row really holds --------------------------------------
worst_seal, worst_leak, thru = 0.0, 0.0, 0.0
for s, (node, _, _, _, _) in PAIRS.items():
    leaks, scaled = [], []
    for i in range(1, NS):
        stor = CC * (main["p"][node][i] - main["p"][node][i - 1]) / DT
        sum_in = main["q"][f"Held{s}"][i] - main["q"][f"Bleed{s}"][i]
        leaks.append(abs(stor - sum_in))
        scaled.append(abs(stor - sum_in) / (abs(stor) + abs(sum_in) + 1.0))
    worst_seal = max(worst_seal, max(scaled))
    worst_leak = max(worst_leak, max(leaks))
    thru = max(thru, max(main["q"][f"Bleed{s}"]))
require(worst_seal < 5.0 * NEWTON_TOL,
        f"all {len(PAIRS)} floating nodes close to the solver's own exit "
        f"criterion, not to a tolerance chosen here: worst "
        f"{worst_seal:.3e} against the {NEWTON_TOL:.0e} the Newton loop stops on "
        f"(the +1.0 in the node row's scale, transient.cpp:414) -- the held pipe "
        f"enters continuity as the forced zero it prints")
require(worst_leak / thru < 1e-5,
        f"and that residual is a solver number, not a leak: {worst_leak:.3e} "
        f"m3/s is {worst_leak / thru:.3e} of the largest rate any node moves "
        f"({thru:.3e}); as a leak it would take "
        f"{NODE_VOLUME * 5000.0 / NODE_BULK / worst_leak:.1e} s to move a node "
        f"the 5000 Pa its own solve covers in {END:.0f} s")

# --- 9. the calendar against the physics ------------------------------------
# Everything above is read inside the 4 s window. This family reads it against
# two clocks the window does not contain: the steady kernel, which has no latch
# to wait for, and a 30 s march of the two pairs the window never released.
def steady_root(pair):
    """This file's own bisection on the settled series state with every valve
    frozen at its end-of-schedule K. One valve coefficient rather than two: the
    released pair is one series loop, so both local losses multiply the same
    q^2 and `Kh + Kb` is the whole of the schedule's remaining say."""
    node, kh, sh, kb, sb = PAIRS[pair]
    kf = valve_k(sh, kh, END) + valve_k(sb, kb, END)
    ask = lambda tw: (WALL_SCALE * tw
                      + R_LAM * q_of_tw(tw, TAU0)
                      + KV_SCALE * kf * q_of_tw(tw, TAU0) ** 2)
    tw_hi = TAU0 + 1.0
    while ask(tw_hi) < P_A - P_B:
        tw_hi *= 2.0
    lo = TAU0
    for _ in range(200):
        mid = 0.5 * (lo + tw_hi)
        if ask(mid) > P_A - P_B:
            tw_hi = mid
        else:
            lo = mid
    return q_of_tw(0.5 * (lo + tw_hi), TAU0)


e_root = {s: rel(steady[f"Held{s}"]["q"], steady_root(s)) for s in PAIRS}
require(max(e_root.values()) < 5e-8,
        f"the steady kernel's {len(PAIRS)} rates are this file's own bisection on "
        f"1e4 Pa = 4 L tw / d + 128 mu L Q/(pi d^4) + (Kh + Kb) rho Q^2/(2 A^2) "
        f"with each K frozen at its end-of-schedule value: worst "
        f"{max(e_root.values()):.2e} on {max(e_root, key=e_root.get)}, "
        f"{min(e_root.values()):.2e} where every valve ends open")

# 9b: settled, a released pair is an ordinary series loop and its books balance
frozen = []
for s in PAIRS:
    for name, k in ((f"Held{s}", valve_k(PAIRS[s][2], PAIRS[s][1], END)),
                    (f"Bleed{s}", valve_k(PAIRS[s][4], PAIRS[s][3], END))):
        if k > 0.0:
            frozen.append((name, k))
require(len(frozen) == 6
        and all(n in steady and "dpValve" in steady[n] for n, _ in frozen),
        f"this file's own schedule copy leaves exactly {len(frozen)} valves shut "
        f"at t = {END:.0f} s ({', '.join(n for n, _ in frozen)}) and the steady "
        f"kernel prints a `dpValve` on each of them")
worst_sum, worst_series = 0.0, 0.0
for s in PAIRS:
    h, b = steady[f"Held{s}"], steady[f"Bleed{s}"]
    worst_sum = max(worst_sum, abs(h["dp"] + b["dp"] - (P_A - P_B)))
    worst_series = max(worst_series, rel(h["q"], b["q"]))
worst_valve = max(rel(steady[n]["dpValve"], k * KV_SCALE * steady[n]["q"] ** 2)
                  for n, k in frozen)
require(worst_valve < 1e-9,
        f"and each printed drop is the local-loss law with this file's constants: "
        f"K rho Q^2 / (2 A^2) to worst {worst_valve:.2e} over those "
        f"{len(frozen)} pipes, which carry frozen 5e3, 1e6 and 1e8")
require(worst_sum < 1e-6 and worst_series < 1e-9,
        f"the {len(PAIRS)} settled pairs hand the whole {P_A - P_B:.0f} Pa across "
        f"the node to the two series pipes (worst |dp_h + dp_b - 1e4| "
        f"{worst_sum:.2e} Pa) with one rate through both (worst "
        f"{worst_series:.2e}) -- the wall's {TAU0:.1f} Pa of yield stress is a "
        f"number in dp_h = 4 L tw / d like any other")

# 9c: where the valve sits is not in the settled rate, and is in the node
def spread(names):
    return max(rel(steady[f"Held{a}"]["q"], steady[f"Held{b}"]["q"])
               for i, a in enumerate(names) for b in names[i + 1:])


open_end, one_big = ("K0", "W", "OW"), ("K8", "R1", "R2", "O")
require(spread(open_end) < 1e-11,
        f"the three pairs whose valves all end open -- K frozen at 0, the W cycle "
        f"on the neighbour, the W cycle on the held pipe itself -- settle to one "
        f"rate to every digit the solver prints, {steady['HeldK0']['q']:.8e} m3/s")
require(spread(one_big) < 5e-8,
        f"and the four placements of one 1e8 valve (frozen on the bleed, ramped "
        f"shut on the bleed twice, frozen on the held pipe) settle to one rate "
        f"within {spread(one_big):.2e}: {steady['HeldK8']['q']:.8e} vs "
        f"{steady['HeldO']['q']:.8e} m3/s")
split = steady_p["Fk8"] - steady_p["Fo"]
dv = steady["BleedK8"]["dpValve"]
require(abs(split - dv) < 1e-3,
        f"but a valve is not a valve to the floating node, only to the loop: the "
        f"same trickle leaves Fk8 at {steady_p['Fk8']:.4f} Pa and Fo at "
        f"{steady_p['Fo']:.4f} Pa, and the {split:.5f} Pa between them is the "
        f"valve's own printed drop {dv:.5f} Pa to {abs(split - dv):.1e} Pa -- the "
        f"rate knows Kh + Kb, the node knows which side of it the node is on")

# 9d: how much clock is left decides how close a marched tail has to be
tail_tol = {
    "K0": (1e-7, 1e-8, f"released at step {rel_step['K0']} and nothing else "
                       f"happens, so the last {END - rel_step['K0'] * DT:.2f} s "
                       f"are plain relaxation"),
    "K3": (1e-5, 1e-6, f"released at step {rel_step['K3']}, its neighbour's 5e3 "
                       f"frozen for the whole run"),
    "O":  (1e-4, 1e-4, f"the 1e8 sits on the held row, which is the row the latch "
                       f"just opened"),
    "OW": (1e-4, 1e-4, f"the cycle is on the held row, so the tail comes back "
                       f"along a path K0 never took"),
    "W":  (5e-3, 1e-3, f"sealed until t = 2.1 s, so only "
                       f"{END - SCHED_W[4]:.1f} s of the reopened drain are on "
                       f"the clock"),
}
ratio = {}
for s, (tol_q, tol_p, why) in tail_tol.items():
    node = PAIRS[s][0]
    e_q = rel(main["q"][f"Held{s}"][-1], steady[f"Held{s}"]["q"])
    e_p = rel(main["p"][node][-1], steady_p[node])
    ratio[s] = e_q / e_p
    require(e_q < tol_q and e_p < tol_p,
            f"{s}: the marched tail meets the steady kernel's rate and node "
            f"pressure for the same frozen topology, at the tolerance its own "
            f"clock earns -- rate rel {e_q:.2e} < {tol_q:.0e}, node rel "
            f"{e_p:.2e} < {tol_p:.0e}, because {why}")
require(min(ratio.values()) > 5.0,
        f"and in all five the *rate* is further from home than the pressure "
        f"driving it, by {min(ratio.values()):.1f}x to {max(ratio.values()):.1f}x "
        f"({', '.join(f'{s} {ratio[s]:.1f}' for s in tail_tol)}): the yielded law "
        f"is steep where the yielded core is thin, q ~ (tw - tau_0)^2 at onset, "
        f"so a 1e-5 stress gap is a 1e-3 rate gap and the gates above are per "
        f"pair rather than one number")

# 9e: two pairs the window left mid-flight, named as mid-flight
e_q6 = rel(main["q"]["HeldK6"][-1], steady["HeldK6"]["q"])
fall6 = (main["p"]["Fk6"][-1] - main["p"]["Fk6"][-5]) / (4.0 * DT)
require(0.9 < e_q6 and main["p"]["Fk6"][-1] > steady_p["Fk6"],
        f"K6 released at step {rel_step['K6']} = t {rel_step['K6'] * DT:.2f} s, so "
        f"the window ends {END - rel_step['K6'] * DT:.2f} s after the wall let go "
        f"and the checker does not pretend it settled: tail rel {e_q6:.3f}, node "
        f"{main['p']['Fk6'][-1]:.4f} Pa still above its {steady_p['Fk6']:.4f} Pa "
        f"while the last five samples fall at {-fall6:.1f} Pa/s, i.e. "
        f"{(main['p']['Fk6'][-1] - steady_p['Fk6']) / -fall6:.2f} s of the same "
        f"drain are still owed")
require(all(main["q"]["HeldK6"][i] > main["q"]["HeldK6"][i - 1]
            for i in range(NS - 5, NS))
        and all(main["p"]["Fk6"][i] < main["p"]["Fk6"][i - 1] for i in range(NS - 5, NS)),
        f"and it is heading home and not wandering: rate "
        f"{main['q']['HeldK6'][-5]:.5e} -> {main['q']['HeldK6'][-1]:.5e} rising "
        f"toward {steady['HeldK6']['q']:.5e} while the node "
        f"{main['p']['Fk6'][-5]:.4f} -> {main['p']['Fk6'][-1]:.4f} falls toward "
        f"{steady_p['Fk6']:.4f}")
e_q2 = rel(main["q"]["HeldR2"][-1], steady["HeldR2"]["q"])
require(0.5 < e_q2 < 0.7 and main["p"]["Fr2"][-1] < steady_p["Fr2"]
        and all(main["q"]["HeldR2"][i] < main["q"]["HeldR2"][i - 1]
                for i in range(NS - 5, NS))
        and all(main["p"]["Fr2"][i] > main["p"]["Fr2"][i - 1]
                for i in range(NS - 5, NS)),
        f"R2 released on K0's step {rel_step['R2']} and then had its drain sealed "
        f"at t = 1.1 s, so it is on the way to the 1e8 trickle from the other "
        f"side: tail rel {e_q2:.3f}, node {main['p']['Fr2'][-1]:.4f} Pa still "
        f"below its {steady_p['Fr2']:.4f} Pa and rising "
        f"({main['p']['Fr2'][-5]:.4f} -> {main['p']['Fr2'][-1]:.4f}), rate "
        f"{main['q']['HeldR2'][-5]:.5e} -> {main['q']['HeldR2'][-1]:.5e} falling "
        f"toward {steady['HeldR2']['q']:.5e}")
for s in NEVER:
    require(steady[f"Held{s}"]["q"] > 0.0 and main["q"][f"Held{s}"][-1] == 0.0,
            f"Held{s}: the transient prints `0` in all {NS} states while the "
            f"steady kernel, which has no latch to wait for, reports a trickle "
            f"{steady[f'Held{s}']['q']:.6e} m3/s at the same node pair -- held "
            f"for {END:.0f} s is the calendar, not the wall")
fall = (main["p"][PAIRS["K8"][0]][-1] - main["p"][PAIRS["K8"][0]][-2]) / DT
require(abs(fall + main["q"]["BleedK8"][-1] / CC) / (main["q"]["BleedK8"][-1] / CC) < 1e-3,
        f"and the calendar is arithmetic: while the wall holds, "
        f"{PAIRS['K8'][0]} loses {-fall * DT:.4f} Pa a step, i.e. "
        f"{-fall:.2f} Pa/s against the bleed's own {main['q']['BleedK8'][-1]:.5e} "
        f"m3/s over {CC:.0e} m3/Pa = "
        f"{main['q']['BleedK8'][-1] / CC:.2f} Pa/s -- so the "
        f"{main['p'][PAIRS['K8'][0]][-1] - P_RELEASE:.0f} Pa still to fall is "
        f"{(main['p'][PAIRS['K8'][0]][-1] - P_RELEASE) / -fall:.1f} s of that "
        f"drain and the window is {END:.0f} s")

# 9f: the same two pairs on a 30 s clock
rel_long = release_step(longn["hold"]["HeldR1"])
drive_long = P_A - longn["p"]["Fk8"][-1]
drives_long = [P_A - v for v in longn["p"]["Fk8"]]
require(rel_long == 2425 and longn["times"][rel_long] == 24.25,
        f"the 30 s run hands R1 back at sample {rel_long} = t "
        f"{longn['times'][rel_long]:.2f} s (entry drive "
        f"{P_A - longn['p']['Fr1'][rel_long - 1]:.4f} Pa over the "
        f"{THRESHOLD:.0f} Pa threshold), against never in the {END:.0f} s window")
require(all(drives_long[i] > drives_long[i - 1] for i in range(1, NL))
        and drive_long < THRESHOLD,
        f"and K8 is still held at {LONG_END:.0f} s with a drive of "
        f"{drive_long:.4f} Pa, {THRESHOLD - drive_long:.2f} Pa short: its node "
        f"falls monotonically from {P_NODE0:.0f} Pa all the way to the "
        f"{LONG_END:.0f} s mark, so `never` is a clock reading and the wall never "
        f"let go")
tail_k8, tail_k8_long = main["q"]["BleedK8"], longn["q"]["BleedK8"]
require(all(tail_k8[i] < tail_k8[i - 1] for i in range(NS - 20, NS))
        and all(tail_k8_long[i] < tail_k8_long[i - 1] for i in range(NL - 5, NL))
        and tail_k8[-1] > tail_k8_long[-1] > steady["BleedK8"]["q"],
        f"and the pair that never released is still creeping, from above: its "
        f"bleed's rate {tail_k8[-1]:.5e} m3/s at {END:.0f} s, "
        f"{tail_k8_long[-1]:.5e} at {LONG_END:.0f} s (both decreasing over their "
        f"last samples) against the steady {steady['BleedK8']['q']:.5e} it is "
        f"still {rel(tail_k8_long[-1], steady['BleedK8']['q']):.3f} away from, "
        f"down from {rel(tail_k8[-1], steady['BleedK8']['q']):.3f} at {END:.0f} s "
        f"-- {LONG_END / END:.1f}x the clock closed "
        f"{1.0 - rel(tail_k8_long[-1], steady['BleedK8']['q']) / rel(tail_k8[-1], steady['BleedK8']['q']):.0%} "
        f"of the gap")
require(all(longn["q"]["HeldR1"][i] > longn["q"]["HeldR1"][i - 1]
            for i in range(NL - 5, NL))
        and longn["q"]["HeldR1"][-1] < steady["HeldR1"]["q"],
        f"and once the latch lets go the released pipe walks to the same trickle "
        f"the latch-free kernel reports: R1's held-row rate "
        f"{longn['q']['HeldR1'][-5]:.5e} -> {longn['q']['HeldR1'][-1]:.5e} m3/s "
        f"rising toward {steady['HeldR1']['q']:.5e}, node drive "
        f"{P_A - longn['p']['Fr1'][-1]:.2f} Pa now over its {THRESHOLD:.0f} Pa")

print()
print(f"transientYieldValveNodeRelease: {len(RELEASED)} of {len(PAIRS)} pairs "
      f"released inside {END:.0f} s on steps "
      f"{', '.join(f'{s} {rel_step[s]}' for s in RELEASED if rel_step[s])}, "
      f"{NEVER} never, all at one release pressure "
      f"({min(p_at_rel.values()):.1f} - {max(p_at_rel.values()):.1f} Pa against "
      f"{P_RELEASE:.0f} Pa), whole-trace reproduction worst rate rel "
      f"{max(worst_rel.values()):.3e} / node rel {max(worst_node.values()):.3e} / "
      f"abs {max(worst_abs.values()):.3e}, node seal worst {worst_seal:.3e} "
      f"against {NEWTON_TOL:.0e}; nine pairs settle to "
      f"{steady['HeldK0']['q']:.4e} m3/s where the valves end open and "
      f"{steady['HeldK8']['q']:.4e} where one 1e8 stays shut, R1 released at t = "
      f"{longn['times'][rel_long]:.2f} s in the long run and K8 still "
      f"{THRESHOLD - drive_long:.1f} Pa short at {LONG_END:.0f} s")

if failures:
    print(f"FAILED: {len(failures)} assertion(s)")
    for f in failures:
        print("  -", f)
    sys.exit(1)
print("all assertions passed")
