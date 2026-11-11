#!/usr/bin/env python3
"""D5 (5.68): CFL-adaptive time stepping -- the schedule is a property of the
physics, re-derived from the states, not a remembered list.

    dt_next = clamp(cfl * min_k(L_k A_k / |q_k|), dtMin, dtMax)

read off the flows the step just committed. The first step is the case's own
dt; a pipe at rest imposes nothing; a network at rest takes dtMax; t
accumulates (with varying dt there is no index to position one); the last
step is clamped to endTime - t and t snaps to endTime exactly when the clamp
fires. A case without `cfl` marches exactly what 5.66's loop marched.

The case's three pipes make the schedule visit every regime: one step at the
initial dt, two at the dtMax ceiling (the network nearly at rest -- nothing
constrains), two steps where the rule itself sets the step (Fast's budget
tightening as its flow develops), then the dtMin floor for the rest, and the
remainder clamp ends the run on endTime exactly. The asymptotic binder is the
Bingham pipe -- the rule really is a min over the whole network.

What this script checks:
  1. the announcement line names cfl/dtMin/dtMax; `steps=` equals the number
     of h's; the final t lands on endTime exactly.
  2. every h in the schedule is re-derived: h_1 = min(dt, endTime); each
     subsequent h = min(clamp(rule on the previous state's flows), rest) --
     including the rest-at-state-0 ceiling; and the regimes actually occur
     (init, ceiling, rule, floor, remainder clamp).
  3. every state is re-derived from the implicit-Euler row with its own h --
     closed form for the Newtonian pipes, bisection in wall-stress space for
     the Bingham one (the 5.66 machinery, now driven by a schedule the
     checker does not control).
  4. the fixed march is intact: the same case without `cfl` takes
     ceil(endTime/dt) uniform steps, every state re-derives, and its states
     differ from the adaptive run's (the schedule moved).
  5. refusals: `cfl 0` declared, `dtMin 0`, and `dtMax < dtMin` are all
     rejected with naming messages.

usage: check_transient_adaptive_cfl.py <solver> <case>
"""
import math
import os
import re
import subprocess
import sys
import tempfile

SOLVER = sys.argv[1]
CASE = sys.argv[2]
CASE_DIR = CASE.rsplit("/", 1)[0]

RHO = 850.0
MU = 0.1
P_A, P_B = 1.7e5, 2.0e4
DRIVE = P_A - P_B
DT = 0.004
END = 2.0
CFL, DT_MIN, DT_MAX = 0.8, 0.08, 0.1
LENGTH = 20.0
DIAM_Y = 0.02
PIPES = {"Fast": (1.0, 0.015, 0.0), "Slow": (20.0, 0.04, 0.0),
         "Yield": (20.0, 0.02, 1.0)}  # name: (L, d, tau0); tau0 = 0 -> Newtonian

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


def resistance(L, d):
    return 128.0 * MU * L / (math.pi * d ** 4)


def area(d):
    return 0.25 * math.pi * d ** 2


def inertia(L, d, h):
    return RHO * L / (area(d) * h)


def q_of_tw(tw, tau0):
    if tw <= tau0:
        return 0.0
    d_ = tw - tau0
    return (math.pi * (0.5 * DIAM_Y) ** 3 * tw ** -3 / MU * d_ * d_
            * (tau0 * tau0 / 2.0 + 2.0 * tau0 * d_ / 3.0 + d_ * d_ / 4.0))


def step_newtonian(q_prev, L, d, h):
    ih = inertia(L, d, h)
    return (ih * q_prev + DRIVE) / (resistance(L, d) + ih)


def step_yield(q_prev, h):
    ih = inertia(LENGTH, DIAM_Y, h)
    flow = lambda tw: q_of_tw(tw, 1.0)
    residual = lambda tw: (ih * (flow(tw) - q_prev)
                           + 4.0 * LENGTH / DIAM_Y * tw - DRIVE)
    if residual(1.0) >= 0.0:
        return 0.0
    lo, hi = 1.0, 2.0
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


def march_pipe(name, hs):
    L, d, tau0 = PIPES[name]
    q_prev = 0.0
    out = []
    for h in hs:
        q_prev = (step_newtonian(q_prev, L, d, h) if tau0 == 0.0
                  else step_yield(q_prev, h))
        out.append(q_prev)
    return out


def rule_dt(state_q):
    """clamp(cfl * min_k(L A / |q_k|), dtMin, dtMax); a network at rest
    imposes nothing and takes the ceiling."""
    fastest = None
    for name, q in state_q.items():
        if abs(q) <= 0:
            continue
        L, d, _ = PIPES[name]
        budget = L * area(d) / abs(q)
        fastest = budget if fastest is None else min(fastest, budget)
    raw = CFL * fastest if fastest is not None else DT_MAX
    return min(max(raw, DT_MIN), DT_MAX)


def run_case(src):
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "case.dict")
        open(path, "w").write(src)
        return subprocess.run([SOLVER, path], capture_output=True, text=True,
                              cwd=CASE_DIR)


def parse(stdout):
    times, q = [], {n: [] for n in PIPES}
    steps = None
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 2 and tok[0] == "t":
            times.append(float(tok[1]))
        elif len(tok) == 4 and tok[0] == "pipe" and tok[2] == "Q" \
                and tok[1] in PIPES:
            q[tok[1]].append(float(tok[3]))
        elif tok and tok[0] == "transient" and "steps=" in line:
            m = re.search(r"steps=(\d+)", line)
            if m:
                steps = int(m.group(1))
    return {"times": times, "q": q, "steps": steps}


BASE = open(CASE).read()
# The header comment narrates the numbers; key extraction and derived arms
# strip // comments first (the 5.66 lesson), then read only live dict text.
CLEAN = re.sub(r"//[^\n]*", "", BASE)

# --- the adaptive run --------------------------------------------------------
r = run_case(BASE)
ok_run = r.returncode == 0 and "nan" not in (r.stdout + r.stderr).lower()
require(ok_run, f"the adaptive case solves with finite output (rc={r.returncode})")
OUT = parse(r.stdout) if ok_run else {"times": [], "q": {n: [] for n in PIPES},
                                      "steps": None}

require("transient adaptive: cfl 0.8 dtMin 0.08 dtMax 0.1"
        in r.stdout.replace("[chflow] ", ""),
        "the announcement names the rule and its clamps: cfl 0.8, "
        "dtMin 0.08, dtMax 0.1 -- a reader of the log knows which march the "
        "t column records")

times = OUT["times"]
require(OUT["steps"] == len(times) - 1,
        f"transient solved: steps={OUT['steps']} equals the {len(times) - 1} "
        "h's the t column implies")
require(bool(times) and times[-1] == END,
        f"the final state lands on endTime exactly ({times[-1] if times else None}) "
        "-- t snaps to endTime when the remainder clamp fires, it does not "
        "accumulate to 1.9999999999999998")

hs = [times[i] - times[i - 1] for i in range(1, len(times))]
require(len(hs) > 0 and rel(hs[0], min(DT, END)) < 1e-12,
        f"the first step is the case's own dt ({hs[0] if hs else None}) -- "
        "the rule only owns the schedule from the second step on")

# every h re-derived from the previous state's flows (the first h is the
# case's own dt, asserted separately above)
states = [{n: OUT["q"][n][i] for n in PIPES} for i in range(len(times))]
bad_h = []
for i, h in enumerate(hs):
    if i == 0:
        continue
    want = min(rule_dt(states[i]), END - times[i])
    if rel(h, want) > 1e-9:
        bad_h.append((i, h, want))
require(not bad_h,
        f"every one of the {len(hs)} steps is the Courant rule on the "
        "previous state's flows, min-ed with the remaining time: "
        "dt_next = clamp(0.8 * min_k(L A / |q_k|), 0.08, 0.1) -- the "
        "schedule is read off the physics, not remembered"
        + (f"; first mismatch {bad_h[0]}" if bad_h else ""))


def cls(h):
    if rel(h, DT) < 1e-12:
        return "init"
    if rel(h, DT_MAX) < 1e-9:
        return "ceiling"
    if rel(h, DT_MIN) < 1e-9:
        return "floor"
    return "rule"


census = {}
for h in hs:
    census[cls(h)] = census.get(cls(h), 0) + 1
remains = len(hs) > 0 and rel(hs[-1], END - times[-2]) < 1e-9
require(all(k in census for k in ("init", "ceiling", "rule", "floor"))
        and remains,
        f"the schedule visits every regime: {census} -- one initial step, "
        "steps at the ceiling while the network is nearly at rest, steps the "
        "rule itself sets while the budgets tighten, steps clamped at the "
        "floor, and the last step is the remainder clamp")

# every state re-derived from the row with its own h (state 0 is the initial
# rest state, then one state per step)
for name, (L, d, tau0) in PIPES.items():
    seq = [0.0] + march_pipe(name, hs)
    got = OUT["q"][name]
    ok = len(got) == len(seq) and all(
        rel(a, b) < 5e-9 for a, b in zip(got, seq))
    worst = max((rel(a, b) for a, b in zip(got, seq)), default=1.0)
    require(ok,
            f"{name}: all {len(got)} reported states reproduce the "
            f"implicit-Euler row with the per-step h (worst rel {worst:.2e})")

# --- the fixed march is intact ----------------------------------------------
FIXED = CLEAN.replace("cfl 0.8;", "").replace("dtMin 0.08;", "") \
             .replace("dtMax 0.1;", "")
rf = run_case(FIXED)
require(rf.returncode == 0, "the same case without cfl solves (fixed march)")
OUTF = parse(rf.stdout)
n_fixed = int(round(END / DT))
require(OUTF["steps"] == n_fixed and len(OUTF["times"]) == n_fixed + 1,
        f"without cfl the loop is the fixed march: ceil({END}/{DT}) = "
        f"{n_fixed} uniform steps, {n_fixed + 1} states (got steps="
        f"{OUTF['steps']}, {len(OUTF['times'])} states)")
hs_fixed = [DT] * n_fixed
all_ok = True
worst_fixed = 0.0
for name, (L, d, tau0) in PIPES.items():
    seq = [0.0] + march_pipe(name, hs_fixed)
    got = OUTF["q"][name]
    all_ok &= len(got) == len(seq) and all(
        rel(a, b) < 5e-9 for a, b in zip(got, seq))
    worst_fixed = max(worst_fixed,
                      max((rel(a, b) for a, b in zip(got, seq)), default=1.0))
require(all_ok,
        f"and all {n_fixed + 1} states of the fixed march re-derive from the "
        f"row with the uniform h (worst rel {worst_fixed:.2e}) -- 5.66's "
        "contract is untouched when the key is absent")
differs = any(OUTF["q"][n][-1] != OUT["q"][n][-1] for n in PIPES)
require(differs,
        "and the two marches end in different states -- the adaptive schedule "
        "changed the discretization, which is the whole point")

# --- the refusals ------------------------------------------------------------
r0 = run_case(CLEAN.replace("cfl 0.8;", "cfl 0;"))
require(r0.returncode != 0
        and "cfl must be positive" in (r0.stdout + r0.stderr),
        "`cfl 0` declared is refused with a naming message -- present but "
        "nonpositive is a broken case, not a silent fixed march")
r1 = run_case(CLEAN.replace("dtMin 0.08;", "dtMin 0;"))
require(r1.returncode != 0 and "dtMin" in (r1.stdout + r1.stderr),
        "dtMin 0 is refused -- an adaptive march without a floor can "
        "always shrink and never terminate")
r2 = run_case(CLEAN.replace("dtMax 0.1;", "dtMax 0.01;"))
require(r2.returncode != 0 and "dtMax" in (r2.stdout + r2.stderr),
        "dtMax below dtMin is refused -- inverted clamps are a broken case")

print()
print(f"{len(failures)} failure(s)")
if failures:
    sys.exit(1)
print("OK: the schedule is a property of the physics -- every step is the "
      "Courant rule on the previous state, every state is the row with its "
      "own h, the fixed march is untouched, and the broken cases are refused")
