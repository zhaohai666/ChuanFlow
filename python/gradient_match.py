"""D9 increment 2 (architecture.md 5.104): analytic-gradient steady calibration.

Increment 1 (``chflow.history_match``) treated the C++ steady solver as a black
box and searched with a derivative-free Nelder-Mead simplex: every objective
call is one forward solve, and a finite-difference gradient would cost one
solve per parameter.  This increment supplies the *analytic* gradient instead.

The laminar constant-property gathering network is affine in the node
pressures, so the steady Newton's conductance Jacobian ``J = dr/dp``
(``src/1d/src/network.cpp``) is constant at the solution.  Differentiating the
converged residual ``r(p(d), d) = 0`` with the implicit function theorem gives

    dp/dd_k = -J^{-1} (dr/dd_k)|_p ,

and ``dr/dd_k`` is explicit: a pipe's laminar conductance ``g = pi d^4/(128 mu L)``
has ``dg/dd = 4 g/d``, so ``dq_k/dd_k = 4 q_k/d_k``.  One LU factorisation of J
is reused for the forward solve *and* every parameter sensitivity, so the whole
gradient costs a single forward assembly rather than one solve per parameter.

The module is self-contained: only the standard library is used (no NumPy /
SciPy), and it does not import the ``chflow_core`` binding, so it runs in a
bindings-free build.  ``cpp_node_pressures`` wraps the unmodified
``chflow-1d-steady`` executable as an independent forward channel; the
regression checker uses it to confirm the analytic model and its gradient
against the real C++ kernel, not just against a re-derived closed form.
"""

from __future__ import annotations

import math
import os
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

__all__ = [
    "CalibrationProblem",
    "Observation",
    "GradientResult",
    "NODE_NAMES",
    "PIPE_LAYOUT",
    "SOURCE_RATES",
    "OUTLET_PRESSURE",
    "solve_state",
    "forward_pressures",
    "objective",
    "objective_and_gradient",
    "fit",
    "synthetic_observation",
    "cpp_node_pressures",
]

# Topology mirrors tests/cases/steadyNetwork/case.dict and increment 1: two
# rate-specified wells feed a manifold that drains to a fixed-pressure terminal.
NODE_NAMES: Tuple[str, ...] = ("W1", "W2", "M", "Outlet")
SOURCE_RATES: Dict[str, float] = {"W1": 0.002, "W2": 0.001}
OUTLET_PRESSURE = 1.0e6
PIPE_LAYOUT: Tuple[Tuple[str, str, str, float, float], ...] = (
    ("P1", "W1", "M", 1000.0, 5.0),
    ("P2", "W2", "M", 1500.0, 10.0),
    ("P3", "M", "Outlet", 800.0, 0.0),
)
DEFAULT_RHO = 850.0
DEFAULT_MU = 0.1
GRAVITY = 9.80665
DEFAULT_SCALES = (1.0e6, 1.0e-3)


@dataclass(frozen=True)
class CalibrationProblem:
    """Fixed part of the inverse problem; only pipe diameters are unknown."""

    rho: float = DEFAULT_RHO
    mu: float = DEFAULT_MU
    outlet_pressure: float = OUTLET_PRESSURE
    source_rates: Dict[str, float] = None  # type: ignore[assignment]
    pipe_layout: Tuple[Tuple[str, str, str, float, float], ...] = PIPE_LAYOUT
    node_names: Tuple[str, ...] = NODE_NAMES

    def __post_init__(self) -> None:
        if self.rho <= 0 or self.mu <= 0:
            raise ValueError("rho and mu must be positive")
        if self.source_rates is None:
            object.__setattr__(self, "source_rates", dict(SOURCE_RATES))


@dataclass(frozen=True)
class Observation:
    """Structured target fields (node pressures and pipe flow rates)."""

    node_names: Tuple[str, ...]
    node_pressure: Tuple[float, ...]
    pipe_names: Tuple[str, ...]
    pipe_flow_rate: Tuple[float, ...]

    def __post_init__(self) -> None:
        if len(self.node_names) != len(self.node_pressure):
            raise ValueError("node_names and node_pressure have different lengths")
        if len(self.pipe_names) != len(self.pipe_flow_rate):
            raise ValueError("pipe_names and pipe_flow_rate have different lengths")
        if not all(math.isfinite(v) for v in self.node_pressure + self.pipe_flow_rate):
            raise ValueError("observation contains a non-finite value")


@dataclass(frozen=True)
class GradientResult:
    """Outcome of an analytic-gradient Levenberg-Marquardt fit."""

    params: Tuple[float, ...]
    cost: float
    initial_cost: float
    iterations: int
    forward_evals: int
    converged: bool
    gradient_norm: float
    residual_history: Tuple[float, ...]


# --------------------------------------------------------------------------- #
# Dense linear algebra (standard library only): LU with partial pivoting.
# --------------------------------------------------------------------------- #
def _lu_factor(a: List[List[float]]) -> Tuple[List[List[float]], List[int]]:
    n = len(a)
    lu = [row[:] for row in a]
    piv = list(range(n))
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(lu[i][k]))
        if abs(lu[pivot][k]) < 1.0e-300:
            raise ValueError("singular conductance Jacobian (network not anchored?)")
        if pivot != k:
            lu[k], lu[pivot] = lu[pivot], lu[k]
            piv[k], piv[pivot] = piv[pivot], piv[k]
        for i in range(k + 1, n):
            lu[i][k] /= lu[k][k]
            for j in range(k + 1, n):
                lu[i][j] -= lu[i][k] * lu[k][j]
    return lu, piv


def _lu_solve(lu: List[List[float]], piv: List[int], b: Sequence[float]) -> List[float]:
    n = len(lu)
    x = [b[piv[i]] for i in range(n)]
    for i in range(n):
        for j in range(i):
            x[i] -= lu[i][j] * x[j]
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            x[i] -= lu[i][j] * x[j]
        x[i] /= lu[i][i]
    return x


def _validate_params(params: Sequence[float], problem: CalibrationProblem) -> Tuple[float, ...]:
    values = tuple(float(v) for v in params)
    if len(values) != len(problem.pipe_layout):
        raise ValueError(
            "expected {} diameters, got {}".format(len(problem.pipe_layout), len(values))
        )
    if not all(math.isfinite(v) and v > 0 for v in values):
        raise ValueError("all diameters must be finite and positive")
    return values


@dataclass(frozen=True)
class NetworkState:
    """Forward solution plus analytic parameter sensitivities at one diameter vector."""

    node_pressure: Tuple[float, ...]
    pipe_flow_rate: Tuple[float, ...]
    dp_dd: Tuple[Tuple[float, ...], ...]   # [node][param]
    dq_dd: Tuple[Tuple[float, ...], ...]   # [pipe][param]


def solve_state(params: Sequence[float], problem: CalibrationProblem | None = None) -> NetworkState:
    """Assemble the conductance Jacobian, solve once, and reuse its LU for every dp/dd.

    This is the single forward assembly the analytic gradient costs.  The
    residual is affine in the free pressures (laminar), so one Newton step from
    zero is exact; ``J`` is the same matrix the C++ steady Newton assembles.
    """
    prob = problem or CalibrationProblem()
    d = _validate_params(params, prob)
    names = prob.node_names
    n_node = len(names)
    index = {name: i for i, name in enumerate(names)}
    fixed = {names[-1]: prob.outlet_pressure}  # terminal is the pressure reference
    free = [i for i in range(n_node) if names[i] not in fixed]
    free_pos = {node: k for k, node in enumerate(free)}
    n_free = len(free)
    n_param = len(d)

    conductance = [
        math.pi * d[k] ** 4 / (128.0 * prob.mu * prob.pipe_layout[k][3])
        for k in range(n_param)
    ]
    head = [prob.rho * GRAVITY * prob.pipe_layout[k][4] for k in range(n_param)]

    # J = dr/dp and r0 = r(p_free = 0), assembled pipe by pipe exactly as the
    # C++ Newton does (R_a -= Q, R_b += Q; j[a][a] -= c, off-diagonals += c).
    j = [[0.0] * n_free for _ in range(n_free)]
    r0 = [0.0] * n_free
    p_all = [fixed.get(names[i], 0.0) for i in range(n_node)]
    for k, (_name, a_name, b_name, _length, _elev) in enumerate(prob.pipe_layout):
        a, b = index[a_name], index[b_name]
        q = conductance[k] * (p_all[a] - p_all[b] - head[k])
        fa = a in free_pos
        fb = b in free_pos
        if fa:
            r0[free_pos[a]] -= q
        if fb:
            r0[free_pos[b]] += q
        c = conductance[k]
        if fa:
            ia = free_pos[a]
            j[ia][ia] -= c
            if fb:
                ib = free_pos[b]
                j[ia][ib] += c
                j[ib][ia] += c
                j[ib][ib] -= c
    for node, rate in prob.source_rates.items():
        i = index[node]
        if i in free_pos:
            r0[free_pos[i]] += rate

    lu, piv = _lu_factor(j)
    p_free = _lu_solve(lu, piv, [-v for v in r0])
    for k, node in enumerate(free):
        p_all[node] = p_free[k]

    flow = [
        conductance[k] * (p_all[index[prob.pipe_layout[k][1]]]
                          - p_all[index[prob.pipe_layout[k][2]]] - head[k])
        for k in range(n_param)
    ]

    # dp/dd_k = -J^{-1} (dr/dd_k)|_p ; dr/dd_k = -+4 q_k/d_k at the pipe's nodes.
    dp_dd_free: List[List[float]] = []
    for k in range(n_param):
        a, b = index[prob.pipe_layout[k][1]], index[prob.pipe_layout[k][2]]
        dr = [0.0] * n_free
        slope = 4.0 * flow[k] / d[k]
        if a in free_pos:
            dr[free_pos[a]] -= slope
        if b in free_pos:
            dr[free_pos[b]] += slope
        dp_dd_free.append(_lu_solve(lu, piv, [-v for v in dr]))

    dp_dd = [[0.0] * n_param for _ in range(n_node)]
    for k in range(n_param):
        for pos, node in enumerate(free):
            dp_dd[node][k] = dp_dd_free[k][pos]

    # dq_j/dd_k = delta_jk 4 q_j/d_j + g_j (dp_a/dd_k - dp_b/dd_k).
    dq_dd = [[0.0] * n_param for _ in range(n_param)]
    for jj in range(n_param):
        a, b = index[prob.pipe_layout[jj][1]], index[prob.pipe_layout[jj][2]]
        for k in range(n_param):
            explicit = (4.0 * flow[jj] / d[jj]) if jj == k else 0.0
            dq_dd[jj][k] = explicit + conductance[jj] * (dp_dd[a][k] - dp_dd[b][k])

    return NetworkState(
        node_pressure=tuple(p_all),
        pipe_flow_rate=tuple(flow),
        dp_dd=tuple(tuple(row) for row in dp_dd),
        dq_dd=tuple(tuple(row) for row in dq_dd),
    )


def forward_pressures(params: Sequence[float], problem: CalibrationProblem | None = None) -> Dict[str, float]:
    prob = problem or CalibrationProblem()
    state = solve_state(params, prob)
    return {name: state.node_pressure[i] for i, name in enumerate(prob.node_names)}


def _residual_and_jacobian(
    state: NetworkState, target: Observation, scales: Sequence[float]
) -> Tuple[List[float], List[List[float]]]:
    p_scale, q_scale = float(scales[0]), float(scales[1])
    residual: List[float] = []
    jac: List[List[float]] = []
    n_param = len(state.dp_dd[0]) if state.dp_dd else 0
    for i in range(len(target.node_pressure)):
        residual.append((state.node_pressure[i] - target.node_pressure[i]) / p_scale)
        jac.append([state.dp_dd[i][k] / p_scale for k in range(n_param)])
    for j in range(len(target.pipe_flow_rate)):
        residual.append((state.pipe_flow_rate[j] - target.pipe_flow_rate[j]) / q_scale)
        jac.append([state.dq_dd[j][k] / q_scale for k in range(n_param)])
    return residual, jac


def objective(
    params: Sequence[float],
    target: Observation,
    problem: CalibrationProblem | None = None,
    scales: Sequence[float] = DEFAULT_SCALES,
) -> float:
    """Weighted sum of squares between the analytic forward and the target."""
    if len(scales) != 2 or any(not math.isfinite(s) or s <= 0 for s in scales):
        raise ValueError("scales must contain two positive finite numbers")
    prob = problem or CalibrationProblem()
    state = solve_state(params, prob)
    residual, _ = _residual_and_jacobian(state, target, scales)
    return sum(r * r for r in residual)


def objective_and_gradient(
    params: Sequence[float],
    target: Observation,
    problem: CalibrationProblem | None = None,
    scales: Sequence[float] = DEFAULT_SCALES,
) -> Tuple[float, Tuple[float, ...]]:
    """Cost and its analytic gradient in log-diameter space (one forward assembly).

    The gradient is taken with respect to theta = ln d so diameters stay
    positive without a bounds solver; d(cost)/d(theta_k) = d_k d(cost)/d(d_k).
    """
    if len(scales) != 2 or any(not math.isfinite(s) or s <= 0 for s in scales):
        raise ValueError("scales must contain two positive finite numbers")
    prob = problem or CalibrationProblem()
    d = _validate_params(params, prob)
    state = solve_state(d, prob)
    residual, jac = _residual_and_jacobian(state, target, scales)
    cost = sum(r * r for r in residual)
    n_param = len(d)
    grad = [0.0] * n_param
    for k in range(n_param):
        # chain rule into log space, then 2 J^T r
        grad[k] = 2.0 * d[k] * sum(jac[i][k] * residual[i] for i in range(len(residual)))
    return cost, tuple(grad)


def fit(
    target: Observation,
    initial: Sequence[float],
    problem: CalibrationProblem | None = None,
    scales: Sequence[float] = DEFAULT_SCALES,
    max_iter: int = 100,
    tol: float = 1.0e-18,
    x_tol: float = 1.0e-12,
    g_tol: float = 1.0e-14,
    lambda0: float = 1.0e-3,
) -> GradientResult:
    """Levenberg-Marquardt in log-diameter space using the analytic Jacobian.

    Each iteration costs one ``solve_state`` (forward + all sensitivities from a
    shared LU); a finite-difference Jacobian would cost one solve per parameter.
    """
    prob = problem or CalibrationProblem()
    d0 = _validate_params(initial, prob)
    if max_iter < 1 or tol <= 0 or x_tol <= 0 or g_tol <= 0 or lambda0 <= 0:
        raise ValueError("invalid Levenberg-Marquardt controls")
    n = len(d0)
    theta = [math.log(v) for v in d0]
    lam = lambda0
    forward_evals = 0

    def evaluate(th: Sequence[float]) -> Tuple[float, List[float], List[List[float]]]:
        nonlocal forward_evals
        forward_evals += 1
        params = tuple(math.exp(t) for t in th)
        state = solve_state(params, prob)
        residual, jac_d = _residual_and_jacobian(state, target, scales)
        cost = sum(r * r for r in residual)
        # log-space Jacobian: column k scaled by d_k
        jac = [[jac_d[i][k] * params[k] for k in range(n)] for i in range(len(residual))]
        return cost, residual, jac

    cost, residual, jac = evaluate(theta)
    initial_cost = cost
    history = [cost]
    converged = False
    iteration = 0
    grad_norm = max((abs(2.0 * sum(jac[i][k] * residual[i] for i in range(len(residual))))
                     for k in range(n)), default=0.0)

    for iteration in range(1, max_iter + 1):
        if grad_norm <= g_tol or cost <= tol:
            converged = True
            break
        # normal equations (J^T J + lam diag(J^T J)) delta = -J^T r
        jtj = [[sum(jac[i][a] * jac[i][b] for i in range(len(residual)))
                for b in range(n)] for a in range(n)]
        jtr = [sum(jac[i][a] * residual[i] for i in range(len(residual))) for a in range(n)]
        accepted = False
        for _ in range(30):
            damped = [row[:] for row in jtj]
            for a in range(n):
                damped[a][a] += lam * (jtj[a][a] if jtj[a][a] > 0 else 1.0)
            try:
                lu, piv = _lu_factor(damped)
                delta = _lu_solve(lu, piv, [-v for v in jtr])
            except ValueError:
                lam *= 3.0
                continue
            trial = [theta[a] + delta[a] for a in range(n)]
            try:
                new_cost, new_residual, new_jac = evaluate(trial)
            except (ValueError, OverflowError):
                lam *= 3.0
                continue
            if new_cost < cost:
                step = max(abs(delta[a]) for a in range(n))
                theta, cost, residual, jac = trial, new_cost, new_residual, new_jac
                lam = max(lam / 3.0, 1.0e-12)
                history.append(cost)
                grad_norm = max((abs(2.0 * sum(jac[i][k] * residual[i]
                                              for i in range(len(residual))))
                                 for k in range(n)), default=0.0)
                accepted = True
                if step <= x_tol or (history[-2] - cost) <= tol * (1.0 + abs(cost)):
                    converged = True
                break
            lam *= 3.0
        if not accepted:
            break
        if converged:
            break

    params = tuple(math.exp(t) for t in theta)
    return GradientResult(
        params=params,
        cost=cost,
        initial_cost=initial_cost,
        iterations=iteration,
        forward_evals=forward_evals,
        converged=converged or cost <= tol or grad_norm <= g_tol,
        gradient_norm=grad_norm,
        residual_history=tuple(history),
    )


def synthetic_observation(params: Sequence[float], problem: CalibrationProblem | None = None) -> Observation:
    """Noise-free target generated by the analytic forward model."""
    prob = problem or CalibrationProblem()
    state = solve_state(params, prob)
    return Observation(
        node_names=tuple(prob.node_names),
        node_pressure=state.node_pressure,
        pipe_names=tuple(p[0] for p in prob.pipe_layout),
        pipe_flow_rate=state.pipe_flow_rate,
    )


def _case_text(params: Sequence[float], problem: CalibrationProblem) -> str:
    d = _validate_params(params, problem)
    lines = ["fluid {", "    rho {};".format(problem.rho), "    mu {};".format(problem.mu), "}", "nodes {"]
    for name in problem.node_names:
        if name == problem.node_names[-1]:
            lines.append("    {} {{ pressure {}; }}".format(name, problem.outlet_pressure))
        else:
            lines.append("    {} {{}};".format(name))
    lines.append("}")
    lines.append("pipes {")
    for value, (name, a, b, length, elev) in zip(d, problem.pipe_layout):
        lines.append("    {} {{ from {}; to {}; length {}; diameter {}; roughness 0; elevation {}; }}".format(
            name, a, b, length, value, elev))
    lines.append("}")
    lines.append("sources {")
    for i, (node, rate) in enumerate(problem.source_rates.items()):
        lines.append("    src{} {{ node {}; rate {}; }}".format(i, node, rate))
    lines.append("}")
    lines.append("solve { maxIter 100; tol 1e-12; }")
    return "\n".join(lines) + "\n"


def cpp_node_pressures(
    exe: str, params: Sequence[float], problem: CalibrationProblem | None = None
) -> Dict[str, float]:
    """Independent forward channel: run the unmodified chflow-1d-steady kernel.

    Parses only the structured ``node <name> p <value>`` report lines; used by
    the checker to tie the analytic model and its gradient to the real kernel.
    """
    prob = problem or CalibrationProblem()
    text = _case_text(params, prob)
    with tempfile.TemporaryDirectory() as tmp:
        case = os.path.join(tmp, "case.dict")
        with open(case, "w") as fh:
            fh.write(text)
        proc = subprocess.run([exe, case], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError("chflow-1d-steady failed: " + proc.stderr.strip())
    pressures: Dict[str, float] = {}
    for line in proc.stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if len(tok) >= 4 and tok[0] == "node":
            pressures[tok[1]] = float(tok[3])
    if not pressures:
        raise RuntimeError("no node pressures parsed from chflow-1d-steady output")
    return pressures