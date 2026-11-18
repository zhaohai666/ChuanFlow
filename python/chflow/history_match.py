"""D9 increment 1: black-box history matching for the steady network.

The C++ steady solver is treated strictly as a forward black box: a parameter
vector is mapped onto writable ``NetworkPipe`` fields, the solver returns a
``SteadyNetworkResult``, and the residual is assembled from its structured
double fields (``node_pressure`` and ``pipe_flows[].flow_rate``).  Nothing is
parsed from the text report emitted by the CLI applications.

The library deliberately has no SciPy/NumPy dependency.  ``fit`` implements a
deterministic Nelder-Mead simplex in log-parameter space, so diameters stay
positive without an external bounds solver.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Dict, Sequence, Tuple

import chflow

__all__ = [
    "CalibrationProblem",
    "FitResult",
    "Observation",
    "build_network",
    "evaluate",
    "fit",
    "objective",
    "synthetic_observation",
]

# Topology mirrors tests/cases/steadyNetwork/case.dict: two rate-specified
# wells feed a manifold, which drains to a fixed-pressure terminal.  The three
# diameters are the calibration parameters; lengths and elevations are fixed.
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
DEFAULT_SCALES = (1.0e6, 1.0e-3)


@dataclass(frozen=True)
class CalibrationProblem:
    """Fixed part of an inverse problem; only pipe diameters are unknown."""

    rho: float = DEFAULT_RHO
    mu: float = DEFAULT_MU
    outlet_pressure: float = OUTLET_PRESSURE
    source_rates: Dict[str, float] = None  # type: ignore[assignment]
    pipe_layout: Tuple[Tuple[str, str, str, float, float], ...] = PIPE_LAYOUT
    max_iter: int = 100
    tol: float = 1.0e-12

    def __post_init__(self) -> None:
        if self.rho <= 0 or self.mu <= 0:
            raise ValueError("rho and mu must be positive")
        if self.max_iter < 1 or self.tol <= 0:
            raise ValueError("max_iter must be >= 1 and tol must be positive")
        if self.source_rates is None:
            object.__setattr__(self, "source_rates", dict(SOURCE_RATES))


@dataclass(frozen=True)
class Observation:
    """Structured readback of one forward solve (no text report parsing)."""

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
class FitResult:
    """Outcome of a black-box fit, including the best-so-far cost history."""

    params: Tuple[float, ...]
    cost: float
    initial_cost: float
    iterations: int
    evaluations: int
    converged: bool
    residual_history: Tuple[float, ...]


def _validate_params(params: Sequence[float], problem: CalibrationProblem) -> Tuple[float, ...]:
    values = tuple(float(v) for v in params)
    if len(values) != len(problem.pipe_layout):
        raise ValueError(
            "expected {} diameters, got {}".format(len(problem.pipe_layout), len(values))
        )
    if not all(math.isfinite(v) and v > 0 for v in values):
        raise ValueError("all diameters must be finite and positive")
    return values


def build_network(
    params: Sequence[float], problem: CalibrationProblem | None = None
) -> chflow.OneDNetwork:
    """Map a diameter vector onto a fresh writable ``OneDNetwork``."""
    prob = problem or CalibrationProblem()
    values = _validate_params(params, prob)

    net = chflow.OneDNetwork()
    net.fluid = chflow.FluidProperties(prob.rho, prob.mu)
    for name in NODE_NAMES:
        net.add_node(name)
    net.set_node_pressure("Outlet", prob.outlet_pressure)
    for node, rate in prob.source_rates.items():
        net.set_node_rate(node, rate)
    for value, (name, from_name, to_name, length, elevation) in zip(values, prob.pipe_layout):
        pipe = chflow.NetworkPipe()
        pipe.name = name
        pipe.from_name = from_name
        pipe.to_name = to_name
        pipe.length = length
        pipe.diameter = value
        pipe.roughness = 0.0
        pipe.elevation = elevation
        net.add_pipe(pipe)
    return net


def evaluate(
    params: Sequence[float], problem: CalibrationProblem | None = None
) -> Tuple[chflow.OneDNetwork, chflow.SteadyNetworkResult]:
    """Run the unmodified C++ forward kernel for one parameter vector."""
    prob = problem or CalibrationProblem()
    net = build_network(params, prob)
    ctl = chflow.SteadyControls()
    ctl.max_iter = prob.max_iter
    ctl.tol = prob.tol
    result = chflow.solve_steady_network(net, ctl)
    return net, result


def observation_from_result(
    net: chflow.OneDNetwork, result: chflow.SteadyNetworkResult
) -> Observation:
    """Build an ``Observation`` from structured solver outputs only."""
    if not result.converged:
        raise RuntimeError("forward solve did not converge")
    node_names = tuple(net.node_names)
    pressures = tuple(float(result.node_pressure[net.node_index(n)]) for n in node_names)
    flows_by_name = {f.name: float(f.flow_rate) for f in result.pipe_flows}
    pipe_names = tuple(p.name for p in net.pipes)
    return Observation(node_names, pressures, pipe_names, tuple(flows_by_name[n] for n in pipe_names))


def synthetic_observation(
    params: Sequence[float], problem: CalibrationProblem | None = None
) -> Observation:
    """Generate a noise-free observation by running the forward kernel."""
    net, result = evaluate(params, problem)
    return observation_from_result(net, result)


def objective(
    params: Sequence[float],
    target: Observation,
    problem: CalibrationProblem | None = None,
    scales: Sequence[float] = DEFAULT_SCALES,
) -> float:
    """Weighted sum of squares between forward readback and target fields.

    ``scales`` are (pressure, flow-rate) normalisation widths; using the
    natural magnitudes keeps the two field types comparable.  A failed solve
    gets a finite large penalty so the simplex can move away from it.
    """
    if len(scales) != 2 or any(not math.isfinite(s) or s <= 0 for s in scales):
        raise ValueError("scales must contain two positive finite numbers")
    try:
        net, result = evaluate(params, problem)
        if not result.converged:
            return 1.0e30
        got = observation_from_result(net, result)
    except (RuntimeError, ValueError, OverflowError):
        return 1.0e30
    pScale, qScale = float(scales[0]), float(scales[1])
    pErr = sum(
        ((a - b) / pScale) ** 2
        for a, b in zip(got.node_pressure, target.node_pressure)
    )
    qErr = sum(
        ((a - b) / qScale) ** 2
        for a, b in zip(got.pipe_flow_rate, target.pipe_flow_rate)
    )
    return pErr + qErr


def fit(
    target: Observation,
    initial: Sequence[float],
    problem: CalibrationProblem | None = None,
    scales: Sequence[float] = DEFAULT_SCALES,
    step: float = 0.15,
    max_iter: int = 400,
    tol: float = 1.0e-16,
    x_tol: float = 1.0e-12,
) -> FitResult:
    """Deterministic Nelder-Mead fit in log-diameter space.

    The optimiser only calls ``objective``; it never touches solver internals,
    Jacobians, or text output.  ``step`` is the initial log-space simplex
    radius (0.15 is about a 16% diameter perturbation).
    """
    prob = problem or CalibrationProblem()
    values = _validate_params(initial, prob)
    if not (0 < step < 1) or max_iter < 1 or tol <= 0 or x_tol <= 0:
        raise ValueError("invalid Nelder-Mead controls")

    fun: Callable[[Tuple[float, ...]], float] = lambda p: objective(p, target, prob, scales)
    n = len(values)
    y0 = tuple(math.log(v) for v in values)
    simplex = [y0] + [tuple(y0[j] + (step if j == i else 0.0) for j in range(n)) for i in range(n)]
    costs = [fun(tuple(math.exp(v) for v in y)) for y in simplex]
    initial_cost = min(costs)
    best_history = [min(costs)]
    evaluations = n + 1

    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5
    converged = False
    iteration = 0
    for iteration in range(1, max_iter + 1):
        order = sorted(range(n + 1), key=lambda i: costs[i])
        simplex = [simplex[i] for i in order]
        costs = [costs[i] for i in order]
        best_history.append(costs[0])

        spread_y = max(
            max(abs(y[j] - simplex[0][j]) for y in simplex) for j in range(n)
        )
        spread_f = costs[-1] - costs[0]
        if spread_y <= x_tol or spread_f <= tol * (1.0 + abs(costs[0])):
            converged = True
            break

        centroid = tuple(sum(y[j] for y in simplex[:-1]) / n for j in range(n))
        worst = simplex[-1]
        reflected = tuple(centroid[j] + alpha * (centroid[j] - worst[j]) for j in range(n))
        f_reflected = fun(tuple(math.exp(v) for v in reflected))
        evaluations += 1

        if f_reflected < costs[0]:
            expanded = tuple(centroid[j] + gamma * (reflected[j] - centroid[j]) for j in range(n))
            f_expanded = fun(tuple(math.exp(v) for v in expanded))
            evaluations += 1
            simplex[-1], costs[-1] = (
                (expanded, f_expanded) if f_expanded < f_reflected else (reflected, f_reflected)
            )
        elif f_reflected < costs[-2]:
            simplex[-1], costs[-1] = reflected, f_reflected
        else:
            if f_reflected < costs[-1]:
                contracted = tuple(
                    centroid[j] + rho * (reflected[j] - centroid[j]) for j in range(n)
                )
                f_contracted = fun(tuple(math.exp(v) for v in contracted))
                evaluations += 1
                accepted = f_contracted <= f_reflected
            else:
                contracted = tuple(centroid[j] + rho * (worst[j] - centroid[j]) for j in range(n))
                f_contracted = fun(tuple(math.exp(v) for v in contracted))
                evaluations += 1
                accepted = f_contracted < costs[-1]
            if accepted:
                simplex[-1], costs[-1] = contracted, f_contracted
            else:
                for i in range(1, n + 1):
                    simplex[i] = tuple(
                        simplex[0][j] + sigma * (simplex[i][j] - simplex[0][j])
                        for j in range(n)
                    )
                    costs[i] = fun(tuple(math.exp(v) for v in simplex[i]))
                    evaluations += 1

    order = sorted(range(n + 1), key=lambda i: costs[i])
    best = simplex[order[0]]
    best_cost = costs[order[0]]
    if best_history[-1] > best_cost:
        best_history.append(best_cost)
    return FitResult(
        params=tuple(math.exp(v) for v in best),
        cost=best_cost,
        initial_cost=initial_cost,
        iterations=iteration,
        evaluations=evaluations,
        converged=converged,
        residual_history=tuple(best_history),
    )
