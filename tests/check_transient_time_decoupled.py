#!/usr/bin/env python3
"""transientTimeDecoupled (C4): `reportEvery` and `endTime` decoupled, checked
by re-deriving every reported state from the implicit-Euler row with the
PER-STEP h -- including the clamped final step.

Since 5.49 the transient loop stepped `round(endTime/dt)` times with a fixed dt
and sampled only on `i % reportEvery == 0`, which quietly required endTime to
be a multiple of dt. 5.66 changes the loop to ceil(endTime/dt) steps whose last
one is clamped to the remainder (transient.cpp runTransient + transientApp.cpp),
and always reports the final state. A case whose endTime IS a multiple of dt
must be untouched -- that is asserted as structure (state count, step sizes via
the recursion) rather than against a remembered binary.

The row per pipe between two fixed-pressure nodes is

    rho L/A * (Q_new - Q_prev)/h_i = drive - dp_wall(Q_new)

with h_i = min(dt, endTime - (i-1) dt). For the Newtonian pipes the implicit
step is closed-form algebra,

    Q_new - Q_inf = a_i (Q_prev - Q_inf),  a_i = I_h / (R + I_h),
    I_h = rho L / (A h_i),  R = 128 mu L / (pi d^4),  Q_inf = drive / R,

so the whole march collapses to a product of a_i and the end state is checked
to 1e-9 against it; for the yielded pipe the step is the same row with 5.46's
Buckingham wall, solved by bisection in wall-stress space like 5.50's checker.
Every arm is derived from one checked-in dict in a temporary directory; nothing
remembers a printed number.

Six arms plus a refusal:
    clamped   endTime 0.1, report 1: 4 steps 0.03/0.03/0.03/0.01, 5 states
    report3   endTime 0.1, report 3: states at steps 0, 3 and the final 4
    report5   endTime 0.1, report 5: states at 0 and the final 4 -- 4 % 5 != 0
    exact     endTime 0.12, report 1: 4 steps of exactly 0.03 (0.12/0.03 is
              4.000000000000001 in floating point -- the near-multiple trap)
    compat    endTime 0.12, report 2: sample set identical to the old loop's
    shortRun  endTime 0.012: one clamped step of 0.012
    refusal   endTime 0
"""

import math
import os
import re
import subprocess
import sys
import tempfile

TRANSIENT = sys.argv[1]
CASE = sys.argv[2]
CASE_DIR = CASE.rsplit("/", 1)[0]

RHO = 850.0
MU = 0.1
P_A, P_B = 3.0e4, 2.0e4
DRIVE = P_A - P_B
LENGTH = 20.0
DIAM_Y = 0.02

# name: (diameter, tau0); tau0 = 0 -> Newtonian closed form
PIPES = {"Wide": (0.04, 0.0), "Narrow": (0.008, 0.0), "Yield": (0.02, 1.0)}

DT = 0.03
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


def resistance(diam):
    return 128.0 * MU * LENGTH / (math.pi * diam ** 4)


def area(diam):
    return 0.25 * math.pi * diam ** 2


def inertia(diam, h):
    return RHO * LENGTH / (area(diam) * h)


def threshold(tau0):
    return 4.0 * tau0 * LENGTH / DIAM_Y


def q_of_tw(tw, tau0):
    """5.46's Buckingham bracket (the Yield pipe declares no flow index)."""
    if tw <= tau0:
        return 0.0
    d = tw - tau0
    return (math.pi * (0.5 * DIAM_Y) ** 3 * tw ** -3 / MU * d * d
            * (tau0 * tau0 / 2.0 + 2.0 * tau0 * d / 3.0 + d * d / 4.0))


def step_newtonian(q_prev, diam, h):
    r = resistance(diam)
    ih = inertia(diam, h)
    return (ih * q_prev + DRIVE) / (r + ih)


def step_yield(q_prev, h):
    """Implicit step with the static branch, residual side in tw-space: the
    wall drop of the rate the law allows at tw is exactly WALL_SCALE * tw, so
    no inversion is needed inside the residual."""
    ih = inertia(DIAM_Y, h)
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


def h_list(end_time, dt=DT):
    """The loop's own bookkeeping: ceil with the near-multiple tolerance, last
    step clamped -- verbatim from transient.cpp."""
    n = int(math.ceil(end_time / dt - 1e-9))
    hs = []
    for i in range(1, n + 1):
        t_prev = min((i - 1) * dt, end_time)
        hs.append(min(dt, end_time - t_prev))
    return hs


def march(diam, tau0, hs):
    """The full state sequence of the row, one entry per step."""
    q_prev = 0.0
    out = []
    for h in hs:
        q_prev = (step_newtonian(q_prev, diam, h) if tau0 == 0.0
                  else step_yield(q_prev, h))
        out.append(q_prev)
    return out


def close_hs(got, want):
    """Per-step sizes only agree to the last ulp (0.1 - 0.09 is
    0.010000000000000002, 0.12 - 0.09 is 0.030000000000000027), so compare
    with a tolerance, not literals."""
    return len(got) == len(want) and all(
        math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15)
        for a, b in zip(got, want))


def run_case(src):
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "case.dict")
        open(path, "w").write(src)
        r = subprocess.run([TRANSIENT, path], capture_output=True, text=True,
                           cwd=CASE_DIR)
    return r


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
# The header comment narrates example values ("endTime 0.1", "reportEvery 3
# and 5"); key extraction must read only live dict text, so strip // comments
# before regexing and derive each arm by re-replacing on the cleaned text.
SRC = {
    "clamped":  BASE,
    "report3":  BASE.replace("reportEvery 1;", "reportEvery 3;"),
    "report5":  BASE.replace("reportEvery 1;", "reportEvery 5;"),
    "exact":    BASE.replace("endTime 0.1;", "endTime 0.12;"),
    "compat":   BASE.replace("endTime 0.1;", "endTime 0.12;")
                   .replace("reportEvery 1;", "reportEvery 2;"),
    "shortRun": BASE.replace("endTime 0.1;", "endTime 0.012;"),
}


def strip_comments(text):
    return re.sub(r"//[^\n]*", "", text)


END = {k: float(re.search(r"endTime ([0-9.]+)", strip_comments(s)).group(1))
       for k, s in SRC.items()}
REPORT = {k: int(re.search(r"reportEvery (\d+)", strip_comments(s)).group(1))
          for k, s in SRC.items()}
HS = {k: h_list(END[k]) for k in SRC}

OUT = {}
for name, src in SRC.items():
    r = run_case(src)
    require(r.returncode == 0 and "nan" not in (r.stdout + r.stderr).lower(),
            f"{name}: solves with finite output (rc={r.returncode})")
    OUT[name] = parse(r.stdout)

# --- 1. the clamped case: 4 steps of 0.03/0.03/0.03/0.01, 5 states ----------
hs = HS["clamped"]
require(close_hs(hs, [0.03, 0.03, 0.03, 0.01]),
        "the loop's own bookkeeping gives steps 0.03/0.03/0.03/0.01 for "
        "dt 0.03, endTime 0.1 -- the last one is clamped to the remainder")
require(OUT["clamped"]["steps"] == 4,
        f"the app reports steps={OUT['clamped']['steps']}")
require(len(OUT["clamped"]["times"]) == 5,
        f"{len(OUT['clamped']['times'])} states = 1 + 4 (reportEvery 1)")
worst_t = max(abs(t - sum(hs[:i]))
              for i, t in enumerate(OUT["clamped"]["times"]))
require(worst_t < 1e-12,
        "the printed times are 0 / 0.03 / 0.06 / 0.09 / 0.1 to 1e-12 -- the "
        "step index, not an accumulator, positions each step, and the last "
        "one lands on endTime exactly")

# --- 2. every reported state vs the recursion with the per-step h -----------
for name in SRC:
    hs = HS[name]
    report = REPORT[name]
    for n, (diam, tau0) in PIPES.items():
        seq = march(diam, tau0, hs)
        model = [0.0] + [q for i, q in enumerate(seq, 1)
                         if i % report == 0 or i == len(hs)]
        got = OUT[name]["q"][n]
        ok = len(got) == len(model) and all(
            rel(a, b) < 5e-9 for a, b in zip(got, model))
        worst = max((rel(a, b) for a, b in zip(got, model)), default=1.0)
        require(ok,
                f"{name}/{n}: all {len(got)} reported states reproduce the "
                f"implicit-Euler row with the per-step h (worst rel "
                f"{worst:.2e})")

# --- 3. reportEvery no longer owns the final state --------------------------
require(len(OUT["report3"]["times"]) == 3
        and abs(OUT["report3"]["times"][-1] - 0.1) < 1e-12,
        "report3: 3 states on a 4-step run, ending at t = 0.1 (steps 0, 3 and "
        "the final 4)")
require(len(OUT["report5"]["times"]) == 2
        and abs(OUT["report5"]["times"][-1] - 0.1) < 1e-12,
        "report5: 2 states on a 4-step run and the last one is t = 0.1 -- "
        "4 % 5 != 0 and the final state is still reported, which is the "
        "decoupling itself")
require(all(OUT["report5"]["q"][n][-1] == OUT["clamped"]["q"][n][-1]
            for n in PIPES),
        "report5's final state IS the clamped run's step-4 state bit for bit "
        "on all three pipes -- reporting is a window on the same march, not "
        "another march")

# --- 4. the near-multiple trap and the compat twin --------------------------
require(close_hs(HS["exact"], [0.03, 0.03, 0.03, 0.03]),
        "endTime 0.12 with dt 0.03 is 4 steps of EXACTLY dt even though "
        "0.12/0.03 evaluates to 4.000000000000001 -- the near-multiple "
        "tolerance keeps ceil honest")
require(len(OUT["compat"]["times"]) == 3
        and abs(OUT["compat"]["times"][-1] - 0.12) < 1e-12,
        "compat (endTime 0.12, reportEvery 2): 3 states t = 0 / 0.06 / 0.12, "
        "the same sample set the old fixed-step loop produced")

# --- 5. the physics did not move: the march collapses to a product ----------
for name in ("clamped", "exact", "shortRun"):
    hs = HS[name]
    for n, (diam, tau0) in PIPES.items():
        if tau0 > 0:
            continue
        q_inf = DRIVE / resistance(diam)
        prod = 1.0
        for h in hs:
            ih = inertia(diam, h)
            prod *= ih / (resistance(diam) + ih)
        q_end = OUT[name]["q"][n][-1]
        require(rel(q_end, q_inf * (1.0 - prod)) < 1e-9,
                f"{name}/{n}: after {END[name]} s the rate {q_end:.6g} is the "
                f"closed product of per-step relaxations {q_inf * (1.0 - prod):.6g} "
                f"(Q_inf {q_inf:.4g}, {len(hs)} steps) -- the clamp changed the "
                "bookkeeping, not the physics")

# --- 6. the refusal ---------------------------------------------------------
r = run_case(BASE.replace("endTime 0.1;", "endTime 0;"))
require(r.returncode != 0 and "positive" in (r.stdout + r.stderr).lower(),
        "endTime 0 is refused with a naming message -- a zero-length march is "
        "a broken case, not an empty trace")

print()
print(f"{len(failures)} failure(s)")
if failures:
    sys.exit(1)
print("OK: endTime and reportEvery are decoupled -- the last step is clamped "
      "to the remainder, the final state is always reported, and a "
      "multiple-of-dt case is untouched")
