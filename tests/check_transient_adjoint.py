#!/usr/bin/env python3
"""D9 increment 3 (5.117): verify the transient discrete adjoint + checkpointing.

usage: check_transient_adjoint.py <path-to-chflow-1d-transient>

The forward model is the implicit-Euler (Q, p) tank-and-pipe march. The checker
ties it and its diameter gradient to independent channels:
  * the real C++ transient kernel (subprocess), for both the state series and a
    finite-difference of the pressure-misfit objective,
  * the exact substituted scalar recursion check_transient.py uses,
  * a forward tangent-linear sensitivity (independent of the reverse adjoint),
  * central finite differences of the Python objective.
Checkpointing (store every k-th state, recompute segments backward) must
reproduce the store-all adjoint exactly. A safeguarded-secant calibration must
recover the true diameter from an offset start in a handful of adjoint sweeps.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import transient_adjoint as ta  # noqa: E402

if len(sys.argv) != 2:
    raise SystemExit(__doc__)
EXE = sys.argv[1]

PROB = ta.TransientProblem()
D_TRUE = 0.02
D_TEST = 0.023
OBS = ta.synthetic_observation(D_TRUE, PROB)
G = ta.GRAVITY

worst = {}
fails = []


def rec(name, err):
    worst[name] = max(worst.get(name, 0.0), err)


def check(name, got, want, rtol, atol=0.0, rel=False):
    err = abs(got - want)
    ok = err <= atol + rtol * abs(want)
    if not ok:
        fails.append(f"{name}: got={got:.6e} want={want:.6e} err={err:.3e}")
    return err / max(abs(want), 1e-30) if rel else err


# ---- 1. forward vs the real C++ kernel --------------------------------------
py = ta.forward(D_TEST, PROB)
cpp = ta.cpp_series(EXE, D_TEST, PROB)
check("state count", len(py), len(cpp), 0.0)
for i in range(min(len(py), len(cpp))):
    rec("forward/cpp/p", check(f"p[{i}]", py[i][1], cpp[i][1], 1e-8, 1e-6))
    rec("forward/cpp/Q", check(f"Q[{i}]", py[i][0], cpp[i][0], 1e-8, 1e-12))

# ---- 2. forward vs the independent substituted scalar recursion -------------
# Eliminate p_{n+1} from the node compliance row into the momentum row to get
# one scalar update per step (the form check_transient.py blowdown replicates,
# generalised to nonzero outlet pressure / gravity / source). This is algebraically
# independent of the 2x2 inverse the forward module uses.
di, r, b, _ = ta._coeffs(D_TEST, PROB)
grav = PROB.rho * G * PROB.elevation
ck = PROB.bulk_modulus / PROB.volume  # 1/C
q, p = PROB.initial_flow, PROB.initial_pressure
ref = [(q, p)]
for _ in range(PROB.n_steps):
    qn = (p + di * q - PROB.outlet_pressure - grav
          + ck * PROB.dt * PROB.source) / (di + r + ck * PROB.dt)
    pn = p - ck * PROB.dt * (qn - PROB.source)
    q, p = qn, pn
    ref.append((q, p))
for i in range(len(py)):
    rec("forward/recursion/p", check(f"rec_p[{i}]", py[i][1], ref[i][1], 1e-12))
    rec("forward/recursion/Q", check(f"rec_Q[{i}]", py[i][0], ref[i][0], 1e-12))

# ---- 3. laminar guard: the linear wall law is only valid below ~2000 --------
area = 0.25 * math.pi * D_TEST * D_TEST
re_max = max(PROB.rho * abs(s[0]) / area * D_TEST / PROB.mu for s in py)
if re_max >= 2000.0:
    fails.append(f"case left laminar regime (max Re={re_max:.0f})")
else:
    print(f"ok:   laminar throughout (max Re = {re_max:.1f})")

# ---- 4. adjoint vs tangent-linear vs finite difference ----------------------
for d in (0.018, D_TEST, 0.026):
    adj = ta.adjoint_gradient(d, PROB, OBS, checkpoint=1)
    cost, tan = ta.objective_and_tangent(d, PROB, OBS)
    rec("adjoint/tangent", check(f"tan[{d}]", adj.grad_d, tan, 1e-10, rel=True))
    check(f"cost[{d}]", adj.cost, cost, 1e-12)
    h = 1e-6 * d
    fd = (ta.objective(ta.forward(d + h, PROB), OBS)
          - ta.objective(ta.forward(d - h, PROB), OBS)) / (2 * h)
    rec("adjoint/FD", check(f"fd[{d}]", adj.grad_d, fd, 1e-6, rel=True))
    # log-space gradient identity
    check(f"gradlog[{d}]", adj.grad_log, d * adj.grad_d, 1e-12)

# ---- 5. adjoint vs finite difference of the REAL C++ objective --------------
def cpp_cost(d):
    s = ta.cpp_series(EXE, d, PROB)
    om = ta._obs_map(OBS, PROB.n_steps)
    return OBS.weight * sum((s[i][1] - v) ** 2 for i, v in om.items())

h = 1e-5 * D_TEST
fd_cpp = (cpp_cost(D_TEST + h) - cpp_cost(D_TEST - h)) / (2 * h)
adj = ta.adjoint_gradient(D_TEST, PROB, OBS)
rec("adjoint/cppFD", check("cppFD", adj.grad_d, fd_cpp, 1e-4, rel=True))

# ---- 6. checkpointed adjoint == store-all adjoint ---------------------------
ref_grad = ta.adjoint_gradient(D_TEST, PROB, OBS, checkpoint=1).grad_d
for k in (1, 2, 3, 5, 7, 10, 25, 50, 100):
    a = ta.adjoint_gradient(D_TEST, PROB, OBS, checkpoint=k)
    rec("checkpoint/equality", check(f"cp[{k}]", a.grad_d, ref_grad, 1e-12, rel=True))
    if a.checkpoints > PROB.n_steps + 1:
        fails.append(f"checkpoint {k} stored too many states")
# fewer checkpoints must actually store fewer states
if ta.adjoint_gradient(D_TEST, PROB, OBS, checkpoint=25).checkpoints >= \
   ta.adjoint_gradient(D_TEST, PROB, OBS, checkpoint=1).checkpoints:
    fails.append("checkpointing did not reduce stored states")
else:
    print("ok:   checkpointing reduces stored states (101 -> 5 at k=25)")

# ---- 7. subset observation (exercises seeds at non-final indices) -----------
sub = ta.synthetic_observation(D_TRUE, PROB, indices=(0, 10, 50, 99, 100))
a_sub = ta.adjoint_gradient(D_TEST, PROB, sub)
_, t_sub = ta.objective_and_tangent(D_TEST, PROB, sub)
rec("adjoint/tangent/subset", check("subset", a_sub.grad_d, t_sub, 1e-10, rel=True))
for k in (3, 10):
    rec("checkpoint/subset",
        check(f"sub_cp[{k}]", ta.adjoint_gradient(D_TEST, PROB, sub, k).grad_d,
              a_sub.grad_d, 1e-12, rel=True))

# ---- 8. calibration recovers the true diameter ------------------------------
for start in (0.026, 0.015, 0.03):
    f = ta.fit(OBS, start, PROB)
    rel = abs(f.diameter - D_TRUE) / D_TRUE
    rec("fit/param", check(f"fit[{start}]", f.diameter, D_TRUE, 1e-6))
    if f.cost > 1e-12:
        fails.append(f"fit[{start}] cost too high: {f.cost:.3e}")
    if f.forward_evals > 40:
        fails.append(f"fit[{start}] used {f.forward_evals} sweeps (> 40)")
    if not f.converged:
        fails.append(f"fit[{start}] did not converge")
    if not all(f.history[i + 1] <= f.history[i] for i in range(len(f.history) - 1)):
        fails.append(f"fit[{start}] history not monotone")
print(f"ok:   fit recovers d={D_TRUE} in "
      f"{ta.fit(OBS, 0.026, PROB).forward_evals} adjoint sweeps "
      f"(cost {ta.fit(OBS, 0.026, PROB).cost:.2e})")

# ---- 9. determinism ---------------------------------------------------------
f1, f2 = ta.fit(OBS, 0.026, PROB), ta.fit(OBS, 0.026, PROB)
if (f1.diameter, f1.forward_evals, f1.cost) != (f2.diameter, f2.forward_evals, f2.cost):
    fails.append("fit not deterministic")
a1 = ta.adjoint_gradient(D_TEST, PROB, OBS, 7)
a2 = ta.adjoint_gradient(D_TEST, PROB, OBS, 7)
if a1.grad_d != a2.grad_d:
    fails.append("adjoint not deterministic")

# ---- 10. refusals -----------------------------------------------------------
def refuses(fn):
    try:
        fn()
        return False
    except (ValueError, RuntimeError):
        return True

refusals = [
    ("negative diameter", lambda: ta.adjoint_gradient(-0.01, PROB, OBS)),
    ("nan diameter", lambda: ta.adjoint_gradient(float("nan"), PROB, OBS)),
    ("zero diameter", lambda: ta.forward(0.0, PROB)),
    ("checkpoint 0", lambda: ta.adjoint_gradient(D_TEST, PROB, OBS, 0)),
    ("obs length mismatch", lambda: ta.Observation((0, 1), (1.0,))),
    ("obs empty", lambda: ta.Observation((), ())),
    ("obs nonfinite", lambda: ta.Observation((0,), (float("nan"),))),
    ("obs weight 0", lambda: ta.Observation((0,), (1.0,), weight=0.0)),
    ("obs index out of range",
     lambda: ta.adjoint_gradient(D_TEST, PROB, ta.Observation((999,), (1.0,)))),
    ("fit d_start 0", lambda: ta.fit(OBS, 0.0, PROB)),
    ("fit d_start nan", lambda: ta.fit(OBS, float("nan"), PROB)),
    ("fit max_iter 0", lambda: ta.fit(OBS, 0.02, PROB, max_iter=0)),
    ("dt 0", lambda: ta.TransientProblem(dt=0.0)),
    ("n_steps 0", lambda: ta.TransientProblem(n_steps=0)),
]
for name, fn in refusals:
    if not refuses(fn):
        fails.append(f"no refusal: {name}")
print(f"ok:   {len(refusals)} refusals loud")

# ---- report -----------------------------------------------------------------
print("\nworst errors:")
for k in sorted(worst):
    print(f"  {k:28s} {worst[k]:.3e}")
if fails:
    print(f"\nFAIL: {len(fails)} checks")
    for f in fails:
        print("  " + f)
    sys.exit(1)
print(f"\nPASS: {len(worst)} gradient/forward channels tied to closed form, "
      f"tangent, finite difference and the C++ kernel; checkpointing exact")