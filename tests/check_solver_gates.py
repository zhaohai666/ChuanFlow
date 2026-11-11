#!/usr/bin/env python3
"""solverGates (A7): the exit-gate policy is a published contract, asserted
against its single source and against measured footprints.

ChuanFlow's convergence knobs come in two layers. The OUTER layer is the
case's own `solve { tol; maxIter; relax; }` -- a knob a user turns. The INNER
layer is five named constants in src/1d/include/cf/1d/SolverGates.h -- the
per-pipe residual gate, the per-pipe bracket gate, the pipe iteration cap,
the transient Newton gate and its cap. Every closed-form checker since 5.47
reproduces the kernel's states to 1e-9..1e-11 relative on the strength of
exactly these numbers, so they are part of the solver's contract: changing
one is a published-behavior change, not a tuning knob.

Families of judgement:
  1. the five constants in the header equal the documented policy values
  2. the architecture document states each gate by name AND by value, so the
     header and the docs cannot drift apart silently
  3. the call sites reference the named constants, and the old raw literals
     are gone from the files that used to own them
  4. the steady bracket gate's footprint: six pipes of steadyBinghamYield sit
     within 1e-11 relative of Buckingham's exact solution -- the bracket (a
     1e-9 relative gate) actually exits far below its nominal width
  5. the inner gates are FLOORS: tightening the outer tol to 1e-16 changes
     nothing bit for bit, because the pipe solve exits on its own gates long
     before the outer test could ask for more
  6. the outer gate still governs reporting: maxIter 0 fails loudly instead
     of quietly printing an unconverged state as if it were an answer
  7. the transient gate's footprint: every reported state of the decoupled
     march reproduces the exact implicit row to < 1e-9, and at least one
     state differs at all -- the Newton exits on the gate, not on a fixed
     point (the 5.49 statement, in measurable form)
"""

import math
import os
import re
import subprocess
import sys
import tempfile

TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(TESTS)
HEADER = os.path.join(ROOT, "src/1d/include/cf/1d/SolverGates.h")
NETWORK = os.path.join(ROOT, "src/1d/src/network.cpp")
TRANSIENT_SRC = os.path.join(ROOT, "src/1d/src/transient.cpp")
ARCH = os.path.join(ROOT, "docs/architecture.md")

STEADY = os.path.abspath(sys.argv[1])
CASE = os.path.abspath(sys.argv[2])
CASE_DIR = CASE.rsplit("/", 1)[0]
TRANSIENT = os.path.join(os.path.dirname(STEADY), "chflow-1d-transient")
T_CASE = os.path.join(TESTS, "cases/transientTimeDecoupled/case.dict")

# The policy, as documented. If SolverGates.h changes any of these five
# numbers, this file and the architecture document must change with it --
# that is the whole point of the checker.
POLICY = {
    "steadyPipeResidualGate": "1e-13",
    "steadyPipeBracketGate": "1e-9",
    "steadyPipeMaxIter": "100",
    "transientNewtonGate": "1e-11",
    "transientNewtonMaxIter": "60",
}

# steadyBinghamYield's parallel block: six pipes at an imposed 1.0e4 Pa,
# checked against Buckingham's exact laminar solution (5.46's checker owns
# the full case; here the same numbers are the bracket gate's footprint).
RHO, MU, D, L, DP = 850.0, 0.1, 0.02, 20.0, 1.0e4
YIELD_PIPES = {"Xi02": 0.5, "Xi04": 1.0, "Xi06": 1.5,
               "Xi08": 2.0, "Xi096": 2.4}

# transientTimeDecoupled's Newtonian rows: drive 1e4, fixed-pressure ends.
T_RHO, T_MU, T_L, T_DRIVE = 850.0, 0.1, 20.0, 1.0e4
T_PIPES = {"Wide": 0.04, "Narrow": 0.008}

failures = []


def require(cond, what):
    if cond:
        print("ok:  ", what)
    else:
        print("FAIL:", what)
        failures.append(what)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def run(binary, path):
    return subprocess.run([binary, path], capture_output=True, text=True,
                          cwd=CASE_DIR)


# --- 1. the header equals the documented policy ------------------------------
header = open(HEADER).read()
for name, value in POLICY.items():
    m = re.search(rf"{name}\s*=\s*([^;]+);", header)
    require(m is not None and m.group(1).strip() == value,
            f"SolverGates.h: {name} = {value} (the single source says exactly"
            f" what the policy documents)")

# --- 2. the architecture document states each gate, by name and value -------
arch = open(ARCH).read()
for name, value in POLICY.items():
    lines = [ln for ln in arch.splitlines() if name in ln]
    require(any(value in ln for ln in lines),
            f"architecture.md names {name} together with its value {value} --"
            f" the documented policy and the header cannot drift apart in"
            f" silence")

# --- 3. the call sites, and no second source ---------------------------------
net = open(NETWORK).read()
trn = open(TRANSIENT_SRC).read()
for name in ("steadyPipeMaxIter", "steadyPipeResidualGate",
             "steadyPipeBracketGate"):
    require(f"gates::{name}" in net,
            f"network.cpp solves through gates::{name} -- the gate is the"
            f" named constant, not a literal next to the loop")
for name in ("transientNewtonMaxIter", "transientNewtonGate"):
    require(f"gates::{name}" in trn,
            f"transient.cpp steps through gates::{name}")
require("1e-13" not in net,
        "network.cpp no longer contains the raw 1e-13 -- the residual gate"
        " has one source (the two-phase family keeps its own 1e-7/1e-12"
        " gates by declared scope, outside this header)")
require("1e-11" not in trn,
        "transient.cpp no longer contains the raw 1e-11 -- the transient"
        " gate has one source")

# --- 4. the steady bracket gate's footprint ----------------------------------
def g_xi(xi):
    return 1.0 - 4.0 * xi / 3.0 + xi ** 4 / 3.0


def q_bingham(tau0):
    tw = DP * D / (4.0 * L)
    return (math.pi * D ** 3 / (32.0 * MU)) * tw * g_xi(tau0 / tw)


def parse_q(stdout):
    out = {}
    for line in stdout.splitlines():
        tok = line.replace("[chflow] ", "").split()
        if len(tok) >= 4 and tok[0] == "pipe" and tok[2] == "Q":
            try:
                out[tok[1]] = float(tok[3])
            except ValueError:
                pass
    return out


base = run(STEADY, CASE)
require(base.returncode == 0
        and "converged)" in base.stdout and "NOT" not in base.stdout,
        f"steadyBinghamYield solves and reports converged under the default"
        f" outer gates (rc {base.returncode})")
QB = parse_q(base.stdout)
worst = max(rel(QB[n], q_bingham(t0)) for n, t0 in YIELD_PIPES.items())
newt = math.pi * D ** 4 * DP / (128.0 * MU * L)
worst = max(worst, rel(QB["Newt"], newt))
require(worst <= 1e-11,
        f"all six parallel pipes sit within {worst:.3e} relative of the exact"
        f" closed form -- the bracket's nominal 1e-9 width is not what the"
        f" answer feels: the false-position exit lands orders below it, and"
        f" that margin IS the bracket gate's measured footprint")

# --- 5. the inner gates are floors, not knobs --------------------------------
base_src = open(CASE).read()
tight_src = base_src.replace("maxIter 60;", "maxIter 500;") \
                    .replace("tol 1e-12;", "tol 1e-16;")
with tempfile.TemporaryDirectory() as td:
    path = os.path.join(td, "case_tight.dict")
    open(path, "w").write(tight_src)
    tight = run(STEADY, path)
pipe_lines = lambda s: [ln for ln in s.splitlines()
                        if ln.startswith("[chflow] pipe ")]
require(tight.returncode == 0 and "NOT" not in tight.stdout,
        "asking the outer layer for tol 1e-16 still converges -- the outer"
        " knob cannot be turned past what the pipe solve's own gates give")
require(pipe_lines(tight.stdout) == pipe_lines(base.stdout),
        "and the answer is bit-for-bit the default one, every pipe line"
        " identical -- the inner gates are floors the outer tol stands on,"
        " not a second way to move the answer")

# --- 6. the outer gate governs reporting -------------------------------------
cap0_src = base_src.replace("maxIter 60;", "maxIter 0;")
with tempfile.TemporaryDirectory() as td:
    path = os.path.join(td, "case_cap0.dict")
    open(path, "w").write(cap0_src)
    cap0 = run(STEADY, path)
require(cap0.returncode != 0 and "NOT converged" in cap0.stdout,
        "maxIter 0 fails loudly (NOT converged, nonzero exit) -- an"
        " unconverged state is reported as a failure, never laundered into"
        " an answer by the printer")

# --- 7. the transient gate's footprint ---------------------------------------
t_src = open(T_CASE).read()
t_clean = re.sub(r"//[^\n]*", "", t_src)
t_end = float(re.search(r"endTime ([0-9.]+)", t_clean).group(1))
t_dt = float(re.search(r"\bdt ([0-9.]+)", t_clean).group(1))
t_report = int(re.search(r"reportEvery (\d+)", t_clean).group(1))


def resistance(diam):
    return 128.0 * T_MU * T_L / (math.pi * diam ** 4)


def area(diam):
    return 0.25 * math.pi * diam ** 2


def march(diam, hs):
    """The exact implicit row per step, closed form -- no Newton, no gate."""
    r = resistance(diam)
    q = 0.0
    out = []
    for h in hs:
        ih = T_RHO * T_L / (area(diam) * h)
        q = (ih * q + T_DRIVE) / (r + ih)
        out.append(q)
    return out


trun = run(TRANSIENT, T_CASE)
require(trun.returncode == 0,
        f"transientTimeDecoupled solves (rc {trun.returncode}) -- the same"
        f" march 5.66's checker re-derives, here read as gate evidence")
times, qs = [], {n: [] for n in T_PIPES}
for line in trun.stdout.splitlines():
    tok = line.replace("[chflow] ", "").split()
    if len(tok) >= 2 and tok[0] == "t":
        times.append(float(tok[1]))
    elif len(tok) == 4 and tok[0] == "pipe" and tok[2] == "Q" \
            and tok[1] in T_PIPES:
        qs[tok[1]].append(float(tok[3]))
n_steps = int(math.ceil(t_end / t_dt - 1e-9))
hs = [min(t_dt, t_end - min((i - 1) * t_dt, t_end)) for i in
      range(1, n_steps + 1)]
any_nonzero = False
worst_t = 0.0
for n, diam in T_PIPES.items():
    seq = [0.0] + [q for i, q in enumerate(march(diam, hs), 1)
                   if i % t_report == 0 or i == len(hs)]
    got = qs[n]
    ok = len(got) == len(seq)
    if ok:
        errs = [rel(a, b) for a, b in zip(got, seq)]
        worst_t = max(worst_t, max(errs))
        any_nonzero = any_nonzero or any(e > 0.0 for e in errs)
    require(ok and all(rel(a, b) < 1e-9 for a, b in zip(got, seq)),
            f"{n}: all {len(got)} reported states reproduce the exact"
            f" implicit row (gate drift bounded, worst so far {worst_t:.2e})")
require(any_nonzero and worst_t > 0.0,
        f"and the drift is real, not a fixed point: the worst mismatch"
        f" {worst_t:.2e} is nonzero while every state stays under 1e-9 --"
        f" the transient Newton exits ON its 1e-11 gate, and the march"
        f" accumulates that slowly, which is the 5.49 statement made"
        f" measurable")

print()
print(f"{len(failures)} failure(s)")
sys.exit(1 if failures else 0)
