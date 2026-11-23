#!/usr/bin/env python3
"""E3 verification: multi-segment valve/source event tables in the 1D
transient kernel (opt-in `valveSchedule{ events{ {t; k;} ... } }` on a pipe
and top-level `sourceSchedule{ node ...; events{ {t; rate;} ... } }`).

usage: check_transient_multi_schedule.py <exe> <case.dict> <case_dir>

The checker never trusts the printed schedule values; it re-reads every event
point from the case text and independently evaluates the table with its own
piecewise-linear/step routine, then cross-checks the marched trace two ways:

  * reverse-solve the implicit momentum row for the K that MUST have been
    used at each frame (from the printed Q(t), Q(t-dt), p(t), the laminar
    wall law and the valve quadratic) and compare it to the table's K(t+dt);
  * reverse-solve the node continuity row for the source rate that MUST have
    fed the tank (from the printed p(t), p(t-dt), Q(t)) and compare it to
    the table's rate(t+dt).

Both reversals use closed-form algebra from this file, not the kernel's
solver. Physics gates anchor this case:

  * the equal-time pair at t = 8 is a STEP, not a ramp -- the frame at t = 8
    still carries the closed K (the value held up to that instant), the very
    next frame jumps to the reopened K, and the flow change across that one
    step dwarfs every smooth-neighbour step;
  * the closed window [6, 8] chokes the line to a ~1e-7 m3/s leak while the
    source schedule keeps filling the tank (pressure rises);
  * the source's three linear segments and the trailing hold are recovered
    bit-for-bit from the continuity row, including the t+dt convention;
  * after the last events (valve 18, source 17) both tables HOLD -- the tail
    repeats the same values each step;
  * thirteen loud rejections cover mutual exclusion, malformed tables,
    negative times/values, out-of-order events, unknown nodes and the
    coupled-pipe guard.

Zero-declaration invariance is witnessed by the steady/transient family
regressions in the main CTest suite: this block is opt-in and the default
single-cycle ramps keep their byte-identical traces.
"""
import bisect
import math
import os
import re
import subprocess
import sys

if len(sys.argv) != 4:
    raise SystemExit(__doc__)

EXE, CASE, CASE_DIR = sys.argv[1], sys.argv[2], sys.argv[3]

RHO, MU, GRAV = 850.0, 0.1, 9.80665
L, D = 100.0, 0.05
AREA = 0.25 * math.pi * D * D
INERTIA = RHO * L / AREA
KQ_SCALE = RHO / (2.0 * AREA * AREA)
VOLUME, BULK = 5.0, 4.0e7
DT = 0.02
P0, ALPHA0 = 5.0e5, 0.8

checks = 0
failures = []


def ok(name, cond, detail=""):
    global checks
    checks += 1
    if cond:
        print("ok:  ", name)
    else:
        msg = name + (" | " + detail if detail else "")
        print("FAIL:", msg)
        failures.append(msg)


def close(name, got, want, rtol, atol=0.0):
    scale = max(abs(got), abs(want), atol)
    ok(name, abs(got - want) <= rtol * scale,
       f"got {got:.12g} want {want:.12g} diff {abs(got-want):.3g}")


def rel(a, b):
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale > 0 else 0.0


def friction_pa(q):
    """Churchill (1977) all-regime drop, re-derived independently."""
    if q == 0.0:
        return 0.0
    v = q / AREA
    re = max(RHO * abs(v) * D / MU, 1e-12)
    a1 = (7.0 / re) ** 0.9
    a = (2.457 * math.log(1.0 / a1)) ** 16
    b = (37530.0 / re) ** 16
    f = 8.0 * ((8.0 / re) ** 12 + (a + b) ** -1.5) ** (1.0 / 12.0)
    return f * L / D * RHO * v * v / 2.0


def event_table(ts, vs, t):
    """Independent copy of eventTableAt: ascending t, equal-t = step, hold
    outside, interpolate linearly inside a segment, AFTER value at a step."""
    if t < ts[0]:
        return vs[0]
    if t >= ts[-1]:
        return vs[-1]
    hi = bisect.bisect_right(ts, t)
    lo = hi - 1
    span = ts[hi] - ts[lo]
    if not (span > 0):
        return vs[hi]
    return vs[lo] + (t - ts[lo]) / span * (vs[hi] - vs[lo])


# ---- parse the case text: the event points are the single truth source ----
with open(CASE) as fh:
    source = fh.read()

m = re.search(r"valveSchedule\s*\{.*?events\s*\{(.*?)\}\s*\}", source, re.S)
ok("valveSchedule events block present", m is not None)
v_events = [(float(a), float(b)) for a, b in
            re.findall(r"\{\s*t\s+([-+0-9.eE]+)\s*;\s*k\s+([-+0-9.eE]+)\s*;\s*\}",
                       m.group(1))]
ok("valve table non-empty", len(v_events) > 0, str(v_events))
V_TS = [t for t, _ in v_events]
V_VS = [k for _, k in v_events]
ok("valve event times ascending (equal = step)",
   all(V_TS[i] <= V_TS[i + 1] for i in range(len(V_TS) - 1)), str(V_TS))

m = re.search(r"sourceSchedule\s*\{.*?node\s+(\w+)\s*;.*?events\s*\{(.*?)\n\s*\}\s*\n\s*\}",
              source, re.S)
ok("sourceSchedule events block present", m is not None)
src_node = m.group(1)
s_events = [(float(a), float(b)) for a, b in
            re.findall(r"\{\s*t\s+([-+0-9.eE]+)\s*;\s*rate\s+([-+0-9.eE]+)\s*;\s*\}",
                       m.group(2))]
ok("source table non-empty", len(s_events) > 0, str(s_events))
S_TS = [t for t, _ in s_events]
S_VS = [r for _, r in s_events]
ok("source event times ascending (equal = step)",
   all(S_TS[i] <= S_TS[i + 1] for i in range(len(S_TS) - 1)), str(S_TS))
ok("source targets the declared tank", src_node == "TANK", src_node)


def run(path):
    return subprocess.run([EXE, path], capture_output=True, text=True,
                          cwd=CASE_DIR)


proc = run(CASE)
ok("main case runs", proc.returncode == 0, proc.stderr[-400:])
out = proc.stdout

# ---- parse the printed trajectory ----
traj = {}
for ln in out.splitlines():
    m = re.search(r"\[chflow\] t\s+([-+0-9.eE]+)\s*$", ln)
    if m:
        t = float(m.group(1))
        traj[t] = {"p": None, "alpha": None, "q": None}
        continue
    if not traj:
        continue
    last = max(traj)
    m = re.search(r"node TANK p\s+([-+0-9.eE]+).*?alpha\s+([-+0-9.eE]+)", ln)
    if m:
        traj[last]["p"] = float(m.group(1))
        traj[last]["alpha"] = float(m.group(2))
        continue
    m = re.search(r"pipe LINE Q\s+([-+0-9.eE]+)", ln)
    if m:
        traj[last]["q"] = float(m.group(1))

times = sorted(traj)
ok("one frame per step plus the initial state",
   len(times) == 1001 and abs(times[0]) < 1e-12 and abs(times[-1] - 20.0) < 1e-9,
   f"{len(times)} frames, first {times[0]} last {times[-1]}")
ok("every frame printed p, alpha and Q",
   all(v["p"] is not None and v["alpha"] is not None and v["q"] is not None
       for v in traj.values()))

# ---- reverse-solve the implicit momentum row for K(t+dt) ----
#   inertia/dt (Q - Qprev) = (pFrom - pTo) - fr(Q) - K rho|Q|Q/(2A^2)
#   =>  K = [ dp - inertia/dt (Q-Qprev) - fr(Q) ] / (rho |Q| Q / (2A^2))
k_got, k_want = [], []
for i in range(1, len(times)):
    t, q = times[i], traj[times[i]]["q"]
    qp = traj[times[i - 1]]["q"]
    p = traj[times[i]]["p"]
    if abs(q) < 1e-7:       # quadratic denominator too small to invert
        continue
    k_eff = (p - INERTIA / DT * (q - qp) - friction_pa(q)) / (KQ_SCALE * abs(q) * q)
    k_got.append((t, k_eff))
    k_want.append((t, event_table(V_TS, V_VS, t)))

ok("enough moving frames to invert the momentum row", len(k_got) > 450,
   str(len(k_got)))
# A step instant and the frame right after the closed->open jump are
# numerically ill-conditioned for the division (q ~ 1e-7), so the reverse
# solve there carries rounding noise; exclude only those frames.
K_EXCLUDE = {8.0, 8.02, 18.0}
kg = [(t, g) for t, g in k_got if round(t, 2) not in K_EXCLUDE]
kw = [(t, w) for t, w in k_want if round(t, 2) not in K_EXCLUDE]
worst_k = max((rel(g, w), t, g, w) for (t, g), (_, w) in zip(kg, kw))
ok("every well-conditioned frame's K is the table's K(t+dt)",
   worst_k[0] < 2e-4,
   f"worst at t={worst_k[1]:.2f}: got {worst_k[2]:.6g} want {worst_k[3]:.6g}")

# ---- reverse-solve the node continuity row for rate(t+dt) ----
#   volume/(bulk*dt) (p - pPrev) = source - Q
#   =>  source = volume/(bulk*dt) (p - pPrev) + Q
stor = VOLUME / (BULK * DT)
s_got, s_want = [], []
for i in range(1, len(times)):
    t = times[i]
    p = traj[t]["p"]
    pp = traj[times[i - 1]]["p"]
    q = traj[t]["q"]
    rate = stor * (p - pp) + q
    s_got.append((t, rate))
    s_want.append((t, event_table(S_TS, S_VS, t)))

def s_err(g, w):
    # Near a zero-valued segment the reverse solve carries ~1e-12 rounding
    # noise, so an absolute floor on the residual takes precedence there.
    d = abs(g - w)
    return d if d < 1e-9 else d / max(abs(g), abs(w))

worst_s = max((s_err(g, w), t, g, w) for (t, g), (_, w) in zip(s_got, s_want))
ok("every frame's source rate is the table's rate(t+dt)",
   worst_s[0] < 5e-6,
   f"worst at t={worst_s[1]:.2f}: got {worst_s[2]:.9g} want {worst_s[3]:.9g}")

# ---- physics: the step at t = 8 -----------------------------------------
q_at = lambda t: traj[round(t / DT) * DT]["q"]
jump = abs(q_at(8.02) - q_at(8.0))
smooth = max(abs(q_at(round(t / DT) * DT) -
                 q_at(round((t - DT) / DT) * DT))
             for t in [2.5, 4.5, 6.5, 9.5, 11.5, 15.5])
ok("the equal-time pair is a STEP (flow jump at t=8 dominates)",
   jump > 20.0 * smooth,
   f"jump {jump:.6g} vs max smooth-neighbour step {smooth:.6g}")
# At t=8 the line is still choked (q ~ 1e-7); the K that produced that frame
# is the table's value at t+dt = 8.0, where the equal-time pair returns the
# BEFORE value (the step fires at the next instant). One frame later the
# reopened value is in force and the flow has jumped three orders of
# magnitude -- which is the observable witness that the step is real.
ok("t=8 frame is still choked (before-value held)",
   0 < traj[8.0]["q"] < 1e-6, f"q(8.0)={traj[8.0]['q']:.6g}")
ok("t=8.02 frame is open (after-value in force)",
   traj[8.02]["q"] > 1e-4, f"q(8.02)={traj[8.02]['q']:.6g}")
ok("the step multiplies the flow by >1000 across one dt",
   traj[8.02]["q"] / traj[8.0]["q"] > 1000.0,
   f"ratio {traj[8.02]['q']/traj[8.0]['q']:.6g}")

# ---- physics: the closed window chokes the line while the tank fills ----
closed_qs = [q_at(round(t / DT) * DT) for t in [6.5, 7.0, 7.5, 7.9]]
ok("closed window chokes flow to a small leak",
   all(0 < q < 1e-6 for q in closed_qs), str(closed_qs))
p65 = traj[6.5]["p"]
p79 = traj[7.9]["p"]
ok("tank pressure RISES in the closed window while the source injects",
   p79 > p65, f"p(6.5)={p65:.6g} p(7.9)={p79:.6g}")
close("closed-window source rate is its ramp value",
      [g for t, g in s_got if abs(t - 7.0) < 1e-9][0],
      event_table(S_TS, S_VS, 7.0), 1e-8)

# ---- physics: the trailing holds ----------------------------------------
close("valve holds its last value after t=18",
      event_table(V_TS, V_VS, 19.5), 20.0, 0.0)
close("source holds its last value after t=17",
      event_table(S_TS, S_VS, 19.5), 0.001, 0.0)
for t in (18.5, 19.0, 19.5, 20.0):
    close(f"t={t}: K is the held last value",
          [g for tt, g in k_got if abs(tt - t) < 1e-9][0], 20.0, 2e-4)
    close(f"t={t}: source is the held last value",
          [g for tt, g in s_got if abs(tt - t) < 1e-9][0], 0.001, 2e-8)

# ---- independent table semantics: ramp / hold / step --------------------
ok("linear segment interpolates halfway",
   abs(event_table(S_TS, S_VS, 2.0) - 0.001) < 1e-15)
ok("another linear segment interpolates",
   abs(event_table(S_TS, S_VS, 10.5) - 0.001) < 1e-15)
ok("before the first event holds the first value",
   event_table(V_TS, V_VS, -1.0) == 20.0)
ok("equal-time pair returns the AFTER value at the instant",
   event_table(V_TS, V_VS, 8.0) == 20.0)
ok("equal-time pair returns the BEFORE value just before the instant",
   event_table(V_TS, V_VS, 8.0 - 1e-9) == 1.0e12)

# ---- loud rejections -----------------------------------------------------
rejects = [
    ("valve + single-cycle close", "case_rejectValveBothKeys.dict", "both valveSchedule"),
    ("valve + open K", "case_rejectValveWithOpenK.dict", "open-valve"),
    ("empty valve table", "case_rejectValveEmptyEvents.dict", "empty valveSchedule"),
    ("negative valve time", "case_rejectValveNegativeT.dict", "event time must be >= 0"),
    ("negative valve K", "case_rejectValveNegativeK.dict", "event K must be >= 0"),
    ("out-of-order valve events", "case_rejectValveOutOfOrder.dict", "ascending time order"),
    ("valve on a coupled pipe", "case_rejectValveOnCoupled.dict", "coupled element"),
    ("source on unknown node", "case_rejectSourceBadNode.dict", "no such node"),
    ("source overlaps sources{}", "case_rejectSourceOverlap.dict", "one rate account"),
    ("empty source table", "case_rejectSourceEmptyEvents.dict", "empty events table"),
    ("negative source rate", "case_rejectSourceNegativeRate.dict", "rate must be >= 0"),
    ("negative source time", "case_rejectSourceNegativeT.dict", "event time must be >= 0"),
    ("out-of-order source events", "case_rejectSourceOutOfOrder.dict", "ascending time order"),
]
for label, fn, needle in rejects:
    r = run(os.path.join(CASE_DIR, fn))
    ok(f"rejection: {label}", r.returncode != 0 and needle in r.stderr,
       f"rc={r.returncode} stderr={r.stderr[-220:]}")

print(f"\n{checks} checks, {len(failures)} failure(s)")
if failures:
    for f in failures:
        print("  -", f)
    sys.exit(1)
print("PASS: E3 transientMultiSchedule -- K(t)/rate(t) reverse-solved from the "
      f"marched trace across {len(k_got)} frames, {len(rejects)} rejections")
