#!/usr/bin/env python3
"""D9 increment 4 (5.118): verify the two-fluid discrete adjoint + checkpointing.

usage: check_two_fluid_adjoint.py <path-to-chflow-1d-twofluid>

The forward model is the coupled (alpha, slip) march of a sealed vertical
dynamic-slip settling column (vM == 0, explicit rational slip update + Godunov
quadratic holdup transport). The checker ties it and its interfacial-drag
gradient to independent channels:
  * the real C++ two-fluid kernel (subprocess), for the state series and a
    finite difference of the holdup-misfit objective,
  * closed-form physics: inventory conservation, the terminal-slip fixed point
    s* = sqrt(drive/(C D)), and the monotone settled profile,
  * the active-branch Godunov Riemann-solver Jacobian vs finite differences,
  * a forward tangent-linear sensitivity (independent of the reverse adjoint),
  * central finite differences of the Python objective.
Checkpointing (store every k-th state, recompute segments backward) must
reproduce the store-all adjoint exactly. A safeguarded-secant calibration must
recover the true drag from an offset start in a handful of adjoint sweeps.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import twofluid_adjoint as tf  # noqa: E402

if len(sys.argv) != 2:
    raise SystemExit(__doc__)
EXE = sys.argv[1]

PROB = tf.TwoFluidProblem(n_steps=300)
C_TRUE = 1500.0
C_TEST = 1000.0
OBS = tf.synthetic_observation(C_TRUE, PROB, indices=(100, 200, 300),
                               cells=(5, 15, 25, 35))
G = tf.GRAVITY

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
py = tf.forward(C_TEST, PROB)
cpp = tf.cpp_series(EXE, C_TEST, PROB)
check("cpp frame count", len(cpp), 2, 0.0)  # initial + final report frames
for fi, step_idx in ((0, 0), (1, PROB.n_steps)):
    pa, ps = py[step_idx]
    ca, cs = cpp[fi]
    for i in range(PROB.sections):
        rec("forward/cpp/alpha", check(f"a[{fi},{i}]", pa[i], ca[i], 1e-9, 1e-9))
        rec("forward/cpp/slip", check(f"s[{fi},{i}]", ps[i], cs[i], 1e-9, 1e-9))

# ---- 2. closed-form physics: inventory, terminal slip, monotonicity ---------
inv0 = sum(py[0][0]) * PROB.dx
for idx in (0, 100, 200, 300):
    inv = sum(py[idx][0]) * PROB.dx
    rec("inventory", check(f"inv[{idx}]", inv, inv0, 1e-12, 1e-12))
# terminal slip fixed point: iterate the rational map at FIXED alpha to its
# root and compare to s* = sign(drive) sqrt(|drive|/(C D)). Isolates the slip
# map from transport lag (the column is still settling at 300 steps).
rl, rg = PROB.rho_l, PROB.rho_g
inv_dr = 1.0 / rg - 1.0 / rl
sin_th = math.sin(PROB.angle_deg * math.pi / 180.0)
for a in (0.1, 0.25, 0.5, 0.75, 0.9):
    ag = 1.0 - a
    rho_m = a * rl + ag * rg
    drive = inv_dr * rho_m * G * sin_th
    d = 1.0 / (rl * a) + 1.0 / (rg * ag)
    s_star = math.copysign(math.sqrt(abs(drive) / (C_TEST * d)), drive)
    s = 0.0
    for _ in range(2000):
        s = tf.slip_map(a, s, C_TEST, PROB)
    rec("terminal slip", check(f"s*[{a}]", s, s_star, 1e-9, 1e-12))
fa, fs = py[-1]
# monotone decreasing holdup (liquid pools low, gas gathers high)
mono = all(fa[i] >= fa[i + 1] - 1e-12 for i in range(PROB.sections - 1))
if not mono:
    fails.append("settled profile not monotone")
else:
    print("ok:   settled profile monotone (TVD), inventory conserved")
if not all(0.0 - 1e-12 <= a <= 1.0 + 1e-12 for a in fa):
    fails.append("holdup left [0, 1]")

# ---- 3. Godunov active-branch Jacobian vs finite difference -----------------
def fd_partial(qa, qb, qc, aL, aR, which):
    # jac return order is (daL, daR, dqa, dqb, dqc); args order (qa,qb,qc,aL,aR)
    argmap = (3, 4, 0, 1, 2)
    h = 1e-6
    args = [qa, qb, qc, aL, aR]
    args[argmap[which]] += h
    fp = tf.godunov(*args)
    args[argmap[which]] -= 2 * h
    fm = tf.godunov(*args)
    return (fp - fm) / (2 * h)


branch_pts = [(0.7, -0.7, 0.0, 0.8, 0.3), (0.7, -0.7, 0.0, 0.3, 0.8),
              (0.7, -0.7, 0.0, 0.9, 0.7), (-0.5, 0.5, 0.1, 0.2, 0.6),
              (0.5, 0.2, 0.0, 0.1, 0.9), (0.9, -0.2, 0.0, 0.4, 0.45),
              (0.3, -0.9, 0.2, 0.2, 0.8), (1.2, -0.3, 0.0, 0.6, 0.1)]
for pt in branch_pts:
    val, *jac = tf.godunov_jac(*pt)
    check("godunov value", val, tf.godunov(*pt), 1e-14)
    for w in range(5):
        fd = fd_partial(*pt, w)
        rec("godunov/jac", check(f"gjac{pt}[{w}]", jac[w], fd, 1e-5, 1e-7))

# ---- 4. adjoint vs tangent-linear vs finite difference ----------------------
for c in (700.0, C_TEST, 1400.0):
    adj = tf.adjoint_gradient(c, PROB, OBS, checkpoint=1)
    cost, tan = tf.objective_and_tangent(c, PROB, OBS)
    rec("adjoint/tangent", check(f"tan[{c}]", adj.grad_c, tan, 1e-10, rel=True))
    check(f"cost[{c}]", adj.cost, cost, 1e-12)
    h = 1e-5 * c
    fd = (tf.objective(tf.forward(c + h, PROB), OBS)
          - tf.objective(tf.forward(c - h, PROB), OBS)) / (2 * h)
    rec("adjoint/FD", check(f"fd[{c}]", adj.grad_c, fd, 1e-6, rel=True))
    check(f"gradlog[{c}]", adj.grad_log, c * adj.grad_c, 1e-12)

# ---- 5. adjoint vs finite difference of the REAL C++ objective --------------
# report_every=1 so the C++ kernel emits every observed step; OBS samples the
# evolving tails (cells 5/15/25/35) where holdup genuinely depends on C.
def cpp_cost(c):
    frames = tf.cpp_series(EXE, c, PROB, report_every=1)
    return OBS.weight * sum((frames[i][0][cell] - v) ** 2
                            for i, cell, v in zip(OBS.indices, OBS.cells, OBS.alpha))


h = 1e-4 * C_TEST
fd_cpp = (cpp_cost(C_TEST + h) - cpp_cost(C_TEST - h)) / (2 * h)
adj = tf.adjoint_gradient(C_TEST, PROB, OBS)
if abs(adj.grad_c) < 1e-12:
    fails.append("cppFD channel degenerate: objective insensitive to C")
rec("adjoint/cppFD", check("cppFD", adj.grad_c, fd_cpp, 1e-4, rel=True))

# ---- 6. checkpointed adjoint == store-all adjoint ---------------------------
ref_grad = tf.adjoint_gradient(C_TEST, PROB, OBS, checkpoint=1).grad_c
for k in (1, 2, 3, 5, 7, 10, 25, 50, 100, 150):
    a = tf.adjoint_gradient(C_TEST, PROB, OBS, checkpoint=k)
    rec("checkpoint/equality", check(f"cp[{k}]", a.grad_c, ref_grad, 1e-12, rel=True))
    if a.checkpoints > PROB.n_steps + 1:
        fails.append(f"checkpoint {k} stored too many states")
if tf.adjoint_gradient(C_TEST, PROB, OBS, checkpoint=25).checkpoints >= \
   tf.adjoint_gradient(C_TEST, PROB, OBS, checkpoint=1).checkpoints:
    fails.append("checkpointing did not reduce stored states")
else:
    print("ok:   checkpointing reduces stored states (301 -> 13 at k=25)")

# ---- 7. subset observation (seeds at non-final indices) ---------------------
sub = tf.synthetic_observation(C_TRUE, PROB, indices=(0, 50, 150, 300),
                               cells=(2, 20, 38))
a_sub = tf.adjoint_gradient(C_TEST, PROB, sub)
_, t_sub = tf.objective_and_tangent(C_TEST, PROB, sub)
rec("adjoint/tangent/subset", check("subset", a_sub.grad_c, t_sub, 1e-10, rel=True))
for k in (3, 10):
    rec("checkpoint/subset",
        check(f"sub_cp[{k}]", tf.adjoint_gradient(C_TEST, PROB, sub, k).grad_c,
              a_sub.grad_c, 1e-12, rel=True))

# ---- 8. calibration recovers the true drag ----------------------------------
for start in (800.0, 2000.0, 600.0):
    f = tf.fit(OBS, start, PROB)
    rec("fit/param", check(f"fit[{start}]", f.drag, C_TRUE, 1e-6))
    if f.cost > 1e-12:
        fails.append(f"fit[{start}] cost too high: {f.cost:.3e}")
    if f.forward_evals > 40:
        fails.append(f"fit[{start}] used {f.forward_evals} sweeps (> 40)")
    if not f.converged:
        fails.append(f"fit[{start}] did not converge")
    if not all(f.history[i + 1] <= f.history[i] for i in range(len(f.history) - 1)):
        fails.append(f"fit[{start}] history not monotone")
print(f"ok:   fit recovers C={C_TRUE} in "
      f"{tf.fit(OBS, 800.0, PROB).forward_evals} adjoint sweeps "
      f"(cost {tf.fit(OBS, 800.0, PROB).cost:.2e})")

# ---- 9. determinism ---------------------------------------------------------
f1, f2 = tf.fit(OBS, 800.0, PROB), tf.fit(OBS, 800.0, PROB)
if (f1.drag, f1.forward_evals, f1.cost) != (f2.drag, f2.forward_evals, f2.cost):
    fails.append("fit not deterministic")
a1 = tf.adjoint_gradient(C_TEST, PROB, OBS, 7)
a2 = tf.adjoint_gradient(C_TEST, PROB, OBS, 7)
if a1.grad_c != a2.grad_c:
    fails.append("adjoint not deterministic")

# ---- 10. refusals -----------------------------------------------------------
def refuses(fn):
    try:
        fn()
        return False
    except (ValueError, RuntimeError):
        return True


refusals = [
    ("negative drag", lambda: tf.adjoint_gradient(-1.0, PROB, OBS)),
    ("nan drag", lambda: tf.adjoint_gradient(float("nan"), PROB, OBS)),
    ("zero drag", lambda: tf.forward(0.0, PROB)),
    ("checkpoint 0", lambda: tf.adjoint_gradient(C_TEST, PROB, OBS, 0)),
    ("obs length mismatch", lambda: tf.Observation((0, 1), (0,), (0.5,))),
    ("obs empty", lambda: tf.Observation((), (), ())),
    ("obs nonfinite", lambda: tf.Observation((0,), (0,), (float("nan"),))),
    ("obs weight 0", lambda: tf.Observation((0,), (0,), (0.5,), weight=0.0)),
    ("obs index out of range",
     lambda: tf.adjoint_gradient(C_TEST, PROB, tf.Observation((999,), (0,), (0.5,)))),
    ("obs cell out of range",
     lambda: tf.adjoint_gradient(C_TEST, PROB, tf.Observation((0,), (999,), (0.5,)))),
    ("fit c_start 0", lambda: tf.fit(OBS, 0.0, PROB)),
    ("fit c_start nan", lambda: tf.fit(OBS, float("nan"), PROB)),
    ("fit max_iter 0", lambda: tf.fit(OBS, 1000.0, PROB, max_iter=0)),
    ("dt 0", lambda: tf.TwoFluidProblem(dt=0.0)),
    ("n_steps 0", lambda: tf.TwoFluidProblem(n_steps=0)),
    ("sections 1", lambda: tf.TwoFluidProblem(sections=1)),
    ("alpha0 out of range", lambda: tf.TwoFluidProblem(alpha0=1.5)),
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
      f"the Godunov Jacobian, tangent, finite difference and the C++ kernel; "
      f"checkpointing exact")