"""D9 increment 2 (5.104) regression: analytic-gradient steady calibration.

usage: check_gradient_match.py <path-to-chflow-1d-steady>

The harness in python/gradient_match.py is checked three independent ways:
  * its laminar forward against a re-derived Hagen-Poiseuille series/parallel
    closed form (with gravity heads), and against the real C++ steady kernel;
  * its analytic dp/dd and objective gradient against central finite
    differences of BOTH the closed form and the C++ kernel;
  * its Levenberg-Marquardt fit recovering known diameters from an offset start.
It also witnesses that one analytic gradient costs a single forward assembly
(a finite-difference gradient would cost one per parameter), and rejects bad
input. No C++ source is changed; the kernel is wrapped only as a subprocess.
"""

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, "python"))
import gradient_match as gm  # noqa: E402

G = 9.80665
TRUE = (0.055, 0.045, 0.050)
START = (0.040, 0.060, 0.042)

fails = []
worst = {}


def check(label, got, want, rtol):
    denom = max(abs(want), 1.0e-30)
    rel = abs(got - want) / denom
    worst[label] = max(worst.get(label, 0.0), rel)
    if not (rel <= rtol) or not math.isfinite(got):
        fails.append("{}: got {!r} want {!r} rel {:.3e} > {:.1e}".format(label, got, want, rel, rtol))


def closed_form(d, prob):
    """Independent laminar back-substitution from the fixed-pressure terminal."""
    d1, d2, d3 = d
    r = lambda length, dd: 128.0 * prob.mu * length / (math.pi * dd ** 4)
    q1 = prob.source_rates["W1"]
    q2 = prob.source_rates["W2"]
    l1, e1 = prob.pipe_layout[0][3], prob.pipe_layout[0][4]
    l2, e2 = prob.pipe_layout[1][3], prob.pipe_layout[1][4]
    l3, e3 = prob.pipe_layout[2][3], prob.pipe_layout[2][4]
    p_m = prob.outlet_pressure + r(l3, d3) * (q1 + q2) + prob.rho * G * e3
    return {
        "W1": p_m + r(l1, d1) * q1 + prob.rho * G * e1,
        "W2": p_m + r(l2, d2) * q2 + prob.rho * G * e2,
        "M": p_m,
        "Outlet": prob.outlet_pressure,
    }


def main():
    if len(sys.argv) < 2:
        print("usage: check_gradient_match.py <path-to-chflow-1d-steady>")
        return 2
    exe = sys.argv[1]
    prob = gm.CalibrationProblem()

    # 1. analytic forward == independent closed form == C++ kernel
    state = gm.solve_state(TRUE, prob)
    cf = closed_form(TRUE, prob)
    cpp = gm.cpp_node_pressures(exe, TRUE, prob)
    for i, name in enumerate(prob.node_names):
        check("forward/closedform/" + name, state.node_pressure[i], cf[name], 1.0e-12)
        check("forward/cpp/" + name, state.node_pressure[i], cpp[name], 1.0e-8)
    # pipe flows are the source rates (mass balance)
    check("flow/P1", state.pipe_flow_rate[0], prob.source_rates["W1"], 1.0e-12)
    check("flow/P3", state.pipe_flow_rate[2],
          prob.source_rates["W1"] + prob.source_rates["W2"], 1.0e-12)

    # 2. analytic dp/dd == central FD of the closed form AND of the C++ kernel
    for k in range(3):
        h = 1.0e-5 * TRUE[k]
        dp = list(TRUE)
        dm = list(TRUE)
        dp[k] += h
        dm[k] -= h
        cf_p, cf_m = closed_form(dp, prob), closed_form(dm, prob)
        cpp_p = gm.cpp_node_pressures(exe, dp, prob)
        cpp_m = gm.cpp_node_pressures(exe, dm, prob)
        for i, name in enumerate(prob.node_names):
            fd_cf = (cf_p[name] - cf_m[name]) / (2.0 * h)
            fd_cpp = (cpp_p[name] - cpp_m[name]) / (2.0 * h)
            check("dpdd/closedform/{}/{}".format(name, k), state.dp_dd[i][k], fd_cf, 1.0e-6)
            check("dpdd/cpp/{}/{}".format(name, k), state.dp_dd[i][k], fd_cpp, 1.0e-5)

    # rate-specified sources fix the pipe flows, so dq/dd must vanish identically
    # (the explicit 4q/d term cancels the pressure-sensitivity term).
    max_dq = max(abs(v) for row in state.dq_dd for v in row)
    if not max_dq < 1.0e-9:
        fails.append("dq/dd should vanish for a rate-specified network, max {:.3e}".format(max_dq))

    # 3. analytic objective gradient (log space) == central FD of the objective
    target = gm.synthetic_observation(TRUE, prob)
    cost0, grad = gm.objective_and_gradient(START, target, prob)
    check("objective/at-truth", gm.objective(TRUE, target, prob), 0.0, 1.0e-18)
    if not cost0 > 1.0:
        fails.append("objective at offset start should be > 1, got {}".format(cost0))
    theta = [math.log(v) for v in START]
    for k in range(3):
        h = 1.0e-6
        tp, tm = list(theta), list(theta)
        tp[k] += h
        tm[k] -= h
        cp = gm.objective(tuple(math.exp(t) for t in tp), target, prob)
        cm = gm.objective(tuple(math.exp(t) for t in tm), target, prob)
        fd = (cp - cm) / (2.0 * h)
        check("grad/log/{}".format(k), grad[k], fd, 1.0e-5)

    # 4. one analytic gradient costs ONE forward assembly (FD would cost n+1)
    calls = {"n": 0}
    real_solve = gm.solve_state

    def counting_solve(params, problem=None):
        calls["n"] += 1
        return real_solve(params, problem)

    gm.solve_state = counting_solve
    try:
        gm.objective_and_gradient(START, target, prob)
    finally:
        gm.solve_state = real_solve
    if calls["n"] != 1:
        fails.append("analytic gradient used {} forward assemblies, expected 1".format(calls["n"]))

    # 5. Levenberg-Marquardt recovers the true diameters from the offset start
    fit = gm.fit(target, START, prob)
    if not fit.converged:
        fails.append("fit did not converge")
    if not fit.cost < 1.0e-16:
        fails.append("fit cost {} not < 1e-16".format(fit.cost))
    for got, want in zip(fit.params, TRUE):
        check("fit/param", got, want, 1.0e-6)
    if not fit.forward_evals < 60:
        fails.append("fit used {} forward evals, expected < 60".format(fit.forward_evals))
    if not all(b <= a + 1.0e-15 for a, b in zip(fit.residual_history, fit.residual_history[1:])):
        fails.append("residual history is not monotonically non-increasing")

    # 6. determinism
    fit2 = gm.fit(target, START, prob)
    if fit.params != fit2.params or fit.cost != fit2.cost or fit.forward_evals != fit2.forward_evals:
        fails.append("fit is not deterministic")

    # 7. rejections
    for bad in ((0.05, 0.05), (0.05, -0.05, 0.05), (0.05, 0.05, float("nan"))):
        try:
            gm.solve_state(bad, prob)
            fails.append("solve_state accepted bad params {!r}".format(bad))
        except ValueError:
            pass
    try:
        gm.objective(TRUE, target, prob, scales=(0.0, 1.0))
        fails.append("objective accepted a zero scale")
    except ValueError:
        pass
    try:
        gm.fit(target, START, prob, max_iter=0)
        fails.append("fit accepted max_iter=0")
    except ValueError:
        pass

    print("D9 inc2 analytic-gradient calibration: worst relative errors")
    for key in sorted(worst):
        print("  {:32s} {:.3e}".format(key, worst[key]))
    print("  fit cost {:.3e}  evals {}  iters {}  grad_norm {:.3e}".format(
        fit.cost, fit.forward_evals, fit.iterations, fit.gradient_norm))
    if fails:
        print("FAIL ({}):".format(len(fails)))
        for f in fails:
            print("  - " + f)
        return 1
    print("PASS: {} checks, analytic gradient tied to closed form and C++ kernel".format(len(worst)))
    return 0


if __name__ == "__main__":
    sys.exit(main())