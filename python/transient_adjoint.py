"""D9 increment 3 (5.117): discrete stepwise adjoint of the transient 1D
tank-and-pipe march, with checkpointing.

The forward model is the implicit-Euler (Q, p) system the C++ transient kernel
solves for a pressurised tank draining through a laminar line into a
fixed-pressure terminal. In that regime the wall law is exactly linear
(Churchill -> 64/Re to ~1e-40 relative at the case's Re ~ 10), so each step is
a 2x2 linear update

    A s_{n+1} = D s_n + c,   s = (Q, p_tank)
    A = [[I/dt + R, -1], [1, C/dt]],  D = diag(I/dt, C/dt)
    c = (-p_out - rho g elev, source)
    I = rho L / A_pipe = 4 rho L / (pi d^2)   (inertance)
    R = 128 mu L / (pi d^4)                   (laminar resistance)
    C = V / K                                 (node compliance)

This module re-derives that march in pure Python (stdlib only, no chflow /
chflow_core import, so it runs in the bindings-free build) and differentiates
it THREE independent ways:

  * a reverse-mode discrete adjoint carried step by step backwards in time,
  * a forward tangent-linear sensitivity marched alongside the state,
  * central finite differences of the objective.

The adjoint is the O(1)-gradient channel: one backward sweep gives dJ/dd for a
scalar parameter regardless of how many time steps the objective spans. The
tangent-linear and finite-difference channels exist to prove the adjoint
right. Checkpointing (store every k-th state, recompute the segments on the
backward sweep) is verified to reproduce the store-all adjoint exactly.

The parameter is the pipe diameter d (positive by working in ln d). The
objective is a weighted pressure-misfit sum over the observed time samples.
"""

from __future__ import annotations

import math
import subprocess
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

GRAVITY = 9.80665


@dataclass(frozen=True)
class TransientProblem:
    """Fixed geometry / fluid / controls of the tank-and-pipe march. Only the
    diameter is a calibration parameter; everything else is declared."""
    rho: float = 850.0
    mu: float = 0.1
    length: float = 1000.0
    volume: float = 1.0
    bulk_modulus: float = 2.0e9
    outlet_pressure: float = 0.0
    elevation: float = 0.0
    source: float = 0.0
    initial_pressure: float = 5.0e5
    initial_flow: float = 0.0
    dt: float = 0.05
    n_steps: int = 100
    tank_name: str = "TANK"
    outlet_name: str = "OUT"
    pipe_name: str = "LINE"

    def __post_init__(self) -> None:
        if not (self.dt > 0):
            raise ValueError("dt must be positive")
        if self.n_steps < 1:
            raise ValueError("n_steps must be >= 1")
        if not (self.length > 0 and self.volume > 0 and self.bulk_modulus > 0):
            raise ValueError("length, volume and bulkModulus must be positive")
        if not (self.rho > 0 and self.mu > 0):
            raise ValueError("rho and mu must be positive")


@dataclass(frozen=True)
class Observation:
    """Pressure readings at a subset of sample indices (0 = initial state)."""
    indices: Tuple[int, ...]
    pressure: Tuple[float, ...]
    weight: float = 1.0

    def __post_init__(self) -> None:
        if len(self.indices) != len(self.pressure):
            raise ValueError("indices and pressure must have equal length")
        if not self.indices:
            raise ValueError("observation must not be empty")
        if not (self.weight > 0):
            raise ValueError("weight must be positive")
        for p in self.pressure:
            if not math.isfinite(p):
                raise ValueError("observation pressure must be finite")


@dataclass
class AdjointResult:
    diameter: float
    grad_d: float          # dJ/dd
    grad_log: float        # dJ/d(ln d) = d * dJ/dd
    cost: float
    n_states: int
    checkpoints: int       # how many states the backward sweep stored


# ---- the 2x2 step operators -------------------------------------------------

def _coeffs(d: float, prob: TransientProblem) -> Tuple[float, float, float, float]:
    """Return (I/dt, R, C/dt, da/dd) for diameter d. a = I/dt + R is the only
    entry of A that depends on d; da/dd = -(2 I/dt + 4 R)/d."""
    if not (d > 0) or not math.isfinite(d):
        raise ValueError("diameter must be positive and finite")
    area = 0.25 * math.pi * d * d
    inertia = prob.rho * prob.length / area
    resistance = 128.0 * prob.mu * prob.length / (math.pi * d ** 4)
    di = inertia / prob.dt
    b = (prob.volume / prob.bulk_modulus) / prob.dt
    da_dd = -(2.0 * di + 4.0 * resistance) / d
    return di, resistance, b, da_dd


def _apply_inv(a: float, b: float, v0: float, v1: float) -> Tuple[float, float]:
    """A^-1 v for A = [[a,-1],[1,b]]; det = a*b + 1, A^-1 = [[b,1],[-1,a]]/det."""
    det = a * b + 1.0
    return (b * v0 + v1) / det, (-v0 + a * v1) / det


def _step(d: float, prob: TransientProblem, q: float, p: float
          ) -> Tuple[float, float]:
    di, r, b, _ = _coeffs(d, prob)
    a = di + r
    grav = prob.rho * GRAVITY * prob.elevation
    rhs0 = di * q - prob.outlet_pressure - grav
    rhs1 = b * p + prob.source
    return _apply_inv(a, b, rhs0, rhs1)


def forward(d: float, prob: TransientProblem) -> List[Tuple[float, float]]:
    """States s_0..s_N, s_0 = (initial_flow, initial_pressure)."""
    states = [(prob.initial_flow, prob.initial_pressure)]
    q, p = states[0]
    for _ in range(prob.n_steps):
        q, p = _step(d, prob, q, p)
        states.append((q, p))
    return states


# ---- objective --------------------------------------------------------------

def _obs_map(obs: Observation, n_steps: int) -> Dict[int, float]:
    out: Dict[int, float] = {}
    for idx, val in zip(obs.indices, obs.pressure):
        if not (0 <= idx <= n_steps):
            raise ValueError(f"observation index {idx} outside [0, {n_steps}]")
        out[idx] = val
    return out


def objective(states: Sequence[Tuple[float, float]], obs: Observation) -> float:
    om = _obs_map(obs, len(states) - 1)
    return obs.weight * sum((states[i][1] - v) ** 2 for i, v in om.items())


def objective_and_tangent(d: float, prob: TransientProblem, obs: Observation
                          ) -> Tuple[float, float]:
    """Forward tangent-linear channel: march s and s' = ds/dd together.
    s'_{n+1} = M s'_n + ds_{n+1}/dd, M = A^-1 D. Independent of the adjoint."""
    di, r, b, da_dd = _coeffs(d, prob)
    a = di + r
    det = a * b + 1.0
    # M = A^-1 D = [[b*di, b], [-di, a*b]] / det
    m00, m01 = b * di / det, b / det
    m10, m11 = -di / det, a * b / det
    grav = prob.rho * GRAVITY * prob.elevation
    om = _obs_map(obs, prob.n_steps)

    q, p = prob.initial_flow, prob.initial_pressure
    dq, dp = 0.0, 0.0  # s'_0 = 0 (initial state independent of d)
    cost = obs.weight * (p - om[0]) ** 2 if 0 in om else 0.0
    grad = 2.0 * obs.weight * (p - om[0]) * dp if 0 in om else 0.0
    for n in range(1, prob.n_steps + 1):
        rhs0 = di * q - prob.outlet_pressure - grav
        rhs1 = b * p + prob.source
        qn, pn = _apply_inv(a, b, rhs0, rhs1)
        # ds_{n+1}/dd = A^-1[(dD/dd) s_n - (dA/dd) s_{n+1}]. BOTH A and D depend
        # on d: a = I/dt + R with da/dd below, and D = diag(I/dt, C/dt) whose
        # first entry has d(I/dt)/dd = -2 (I/dt)/d. The two source vectors are
        # (-2 di/d Q_n, 0) and (da/dd Q_{n+1}, 0); A^-1 (v0,0) = (b,-1) v0/det.
        fac = (-2.0 * di / d * q - da_dd * qn) / det
        ex0, ex1 = fac * b, fac * (-1.0)
        dqn = m00 * dq + m01 * dp + ex0
        dpn = m10 * dq + m11 * dp + ex1
        q, p, dq, dp = qn, pn, dqn, dpn
        if n in om:
            resid = p - om[n]
            cost += obs.weight * resid * resid
            grad += 2.0 * obs.weight * resid * dp
    return cost, grad


# ---- reverse-mode discrete adjoint with checkpointing -----------------------

def adjoint_gradient(d: float, prob: TransientProblem, obs: Observation,
                     checkpoint: int = 1) -> AdjointResult:
    """One backward sweep gives dJ/dd. checkpoint=1 stores every state;
    checkpoint=k>1 stores every k-th state and recomputes each segment forward
    on the backward sweep (the memory/time tradeoff of transient adjoints).
    Both produce the identical gradient."""
    if checkpoint < 1:
        raise ValueError("checkpoint interval must be >= 1")
    di, r, b, da_dd = _coeffs(d, prob)
    a = di + r
    det = a * b + 1.0
    m00, m01 = b * di / det, b / det
    m10, m11 = -di / det, a * b / det
    # M^T for the backward carry: lambda_n = M^T lambda_{n+1} + dJ/ds_n
    t00, t01 = m00, m10
    t10, t11 = m01, m11
    grav = prob.rho * GRAVITY * prob.elevation
    om = _obs_map(obs, prob.n_steps)
    n = prob.n_steps

    def seed(idx: int, p: float) -> Tuple[float, float]:
        if idx in om:
            return 0.0, 2.0 * obs.weight * (p - om[idx])
        return 0.0, 0.0

    # checkpoint states: indices 0, k, 2k, ... and always n
    cps: Dict[int, Tuple[float, float]] = {0: (prob.initial_flow,
                                               prob.initial_pressure)}
    q, p = cps[0]
    for i in range(1, n + 1):
        q, p = _step(d, prob, q, p)
        if i % checkpoint == 0 or i == n:
            cps[i] = (q, p)
    cost = obs.weight * sum((cps[i][1] - v) ** 2 if i in cps
                            else 0.0 for i, v in om.items())
    # cost needs every observed state; recompute observed ones not stored
    if any(i not in cps for i in om):
        states = forward(d, prob)
        cost = objective(states, obs)

    grad = 0.0
    lam_q, lam_p = seed(n, cps[n][1])
    # segment boundaries descending: last checkpoint <= n is n itself
    bounds = sorted(cps.keys())
    for si in range(len(bounds) - 1, 0, -1):
        end = bounds[si]
        start = bounds[si - 1]
        # recompute states start..end forward from the stored checkpoint
        seg: List[Tuple[float, float]] = [cps[start]]
        qq, pp = cps[start]
        for _ in range(start + 1, end + 1):
            qq, pp = _step(d, prob, qq, pp)
            seg.append((qq, pp))
        # backward within the segment: step index j runs end..start+1
        for j in range(end, start, -1):
            qj, pj = seg[j - start]
            qprev = seg[j - 1 - start][0]
            # explicit param dependence of s_j on d (both A and D vary with d)
            fac = (-2.0 * di / d * qprev - da_dd * qj) / det
            ex0, ex1 = fac * b, fac * (-1.0)
            grad += lam_q * ex0 + lam_p * ex1
            # carry lambda to s_{j-1}: M^T lambda_j + dJ/ds_{j-1}
            sq, sp = seed(j - 1, seg[j - 1 - start][1])
            lam_q, lam_p = t00 * lam_q + t01 * lam_p + sq, \
                t10 * lam_q + t11 * lam_p + sp

    return AdjointResult(diameter=d, grad_d=grad, grad_log=d * grad,
                         cost=cost, n_states=n + 1, checkpoints=len(cps))


# ---- calibration (Levenberg-Marquardt in ln d) ------------------------------

@dataclass
class FitResult:
    diameter: float
    cost: float
    initial_cost: float
    iterations: int
    forward_evals: int
    converged: bool
    history: List[float] = field(default_factory=list)


def synthetic_observation(d_true: float, prob: TransientProblem,
                          indices: Optional[Sequence[int]] = None,
                          weight: float = 1.0) -> Observation:
    states = forward(d_true, prob)
    idx = tuple(range(prob.n_steps + 1)) if indices is None else tuple(indices)
    return Observation(indices=idx,
                       pressure=tuple(states[i][1] for i in idx),
                       weight=weight)


def fit(obs: Observation, d_start: float, prob: TransientProblem,
        max_iter: int = 100, x_tol: float = 1e-13,
        cost_tol: float = 1e-18) -> FitResult:
    """Scalar calibration in theta = ln d. The pressure-misfit is convex in
    theta with a single root of dJ/dtheta = 0, so instead of a slowly-creeping
    Gauss-Newton this finds that root with a safeguarded Illinois secant: each
    iteration costs ONE adjoint backward sweep (the O(1) gradient), and the
    secant converges superlinearly. The bracket is found by expanding from
    d_start in the downhill direction until the gradient changes sign."""
    if not (d_start > 0) or not math.isfinite(d_start):
        raise ValueError("d_start must be positive and finite")
    if max_iter < 1:
        raise ValueError("max_iter must be >= 1")

    evals = 0

    def eval_at(theta: float) -> Tuple[float, float]:
        nonlocal evals
        res = adjoint_gradient(math.exp(theta), prob, obs, checkpoint=1)
        evals += 1
        return res.grad_log, res.cost

    theta0 = math.log(d_start)
    g0, c0 = eval_at(theta0)
    best_theta, best_cost = theta0, c0
    history = [c0]
    if c0 < cost_tol:
        return FitResult(d_start, c0, c0, 0, evals, True, history)

    # expand a bracket [a, b] with g(a) * g(b) <= 0, walking downhill from theta0
    step = 0.1
    ta_, ga, ca = theta0, g0, c0
    tb_ = theta0 - math.copysign(step, g0)
    gb, cb = eval_at(tb_)
    bracketed = ga * gb <= 0
    for _ in range(60):
        if bracketed:
            break
        if cb < best_cost:
            best_theta, best_cost = tb_, cb
        history.append(best_cost)
        step *= 2.0
        ta_, ga, ca = tb_, gb, cb
        tb_ = tb_ - math.copysign(step, g0)
        gb, cb = eval_at(tb_)
        bracketed = ga * gb <= 0
    if not bracketed:
        return FitResult(math.exp(best_theta), best_cost, c0, evals, evals,
                         False, history)

    lo, hi = (ta_, tb_) if ta_ < tb_ else (tb_, ta_)
    g_lo, g_hi = (ga, gb) if ta_ < tb_ else (gb, ga)
    if g_lo > 0:  # keep g_lo < 0 < g_hi (gradient is increasing in theta)
        lo, hi, g_lo, g_hi = hi, lo, g_hi, g_lo

    side = 0
    converged = False
    theta = best_theta
    it = 0
    for it in range(1, max_iter + 1):
        denom = g_hi - g_lo
        theta_s = hi - g_hi * (hi - lo) / denom if denom != 0 else 0.5 * (lo + hi)
        if not (lo < theta_s < hi):
            theta_s = 0.5 * (lo + hi)
        gs, cs = eval_at(theta_s)
        if cs < best_cost:
            best_theta, best_cost = theta_s, cs
        history.append(best_cost)
        theta = theta_s
        if (hi - lo) < x_tol or cs < cost_tol:
            converged = True
            break
        if gs * g_hi > 0:
            hi, g_hi = theta_s, gs
            if side == 1:
                g_lo *= 0.5  # Illinois: halve the retained endpoint's slope
            side = 1
        else:
            lo, g_lo = theta_s, gs
            if side == -1:
                g_hi *= 0.5
            side = -1
    return FitResult(diameter=math.exp(best_theta), cost=best_cost,
                     initial_cost=c0, iterations=it, forward_evals=evals,
                     converged=converged, history=history)


# ---- independent C++ forward channel ----------------------------------------

def case_text(d: float, prob: TransientProblem) -> str:
    return f"""// D9 increment 3 transient-adjoint verification case (generated).
fluid {{ rho {prob.rho!r}; mu {prob.mu!r}; }}
nodes {{
    {prob.tank_name} {{ volume {prob.volume!r}; bulkModulus {prob.bulk_modulus!r};
              pressure {prob.initial_pressure!r}; liquidFraction 1.0; }}
    {prob.outlet_name} {{ fixedPressure {prob.outlet_pressure!r}; }}
}}
pipes {{
    {prob.pipe_name} {{ from {prob.tank_name}; to {prob.outlet_name};
         length {prob.length!r}; diameter {d!r}; roughness 0;
         elevation {prob.elevation!r}; initialFlow {prob.initial_flow!r}; }}
}}
sources {{ s1 {{ node {prob.tank_name}; rate {prob.source!r}; fraction 1; }} }}
time {{ dt {prob.dt!r}; endTime {prob.dt * prob.n_steps!r}; reportEvery 1; }}
"""


def cpp_series(exe: str, d: float, prob: TransientProblem
               ) -> List[Tuple[float, float]]:
    """Run the real C++ transient kernel and read back (Q_LINE, p_TANK) per
    sample, including sample 0. Independent of the Python forward."""
    import tempfile
    import os
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "case.dict")
        with open(path, "w") as fh:
            fh.write(case_text(d, prob))
        proc = subprocess.run([exe, path], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"transient kernel rc={proc.returncode}: "
                           f"{proc.stderr}")
    states: List[Tuple[float, float]] = []
    cur: Dict[str, float] = {}
    have_t = False
    for line in proc.stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "t":
            if have_t:
                states.append((cur[prob.pipe_name], cur[prob.tank_name]))
            have_t = True
            cur = {}
        elif tok[0] == "node" and tok[1] == prob.tank_name:
            cur[prob.tank_name] = float(tok[3])
        elif tok[0] == "pipe" and tok[1] == prob.pipe_name:
            cur[prob.pipe_name] = float(tok[3])
    if have_t:
        states.append((cur[prob.pipe_name], cur[prob.tank_name]))
    return states