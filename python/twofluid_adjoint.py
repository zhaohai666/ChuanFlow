"""D9 increment 4 (5.118): discrete stepwise adjoint of the two-fluid
transient march (sealed dynamic-slip settling column), with checkpointing.

The forward model is the coupled (alpha, slip) march the C++ two-fluid kernel
runs for a sealed vertical pipe with the dynamic slip closure and identically
zero mixture velocity. Because vM == 0 the momentum solve is skipped and each
step is an EXPLICIT differentiable map (no nonlinear solve):

  1. rational slip update, per cell i:
        s_new = (s_old/dt + drive/K_s) / (1/dt + C |s_old| D / K_s)
        rho_m = alpha rhoL + (1-alpha) rhoG
        K_s   = 1 + alpha(1-alpha)(rhoG-rhoL)^2/(rhoL rhoG)
        drive = (1/rhoG - 1/rhoL) rho_m g sin(theta)
        D     = 1/(rhoL alpha) + 1/(rhoG (1-alpha))
  2. Godunov quadratic holdup transport on the averaged face slip:
        s_face = (s_new[j-1] + s_new[j]) / 2
        flux_j = Godunov(F(a) = s_face a^2 - s_face a; alpha[j-1], alpha[j])
        alpha_new[i] = alpha[i] + dt/dx (flux_i - flux_{i+1}),  flux_0=flux_N=0

This module re-derives that march in pure Python (stdlib only, no chflow /
chflow_core import, so it runs in the bindings-free build) and differentiates
it THREE independent ways:

  * a reverse-mode discrete adjoint carried step by step backwards in time
    through BOTH the transport Jacobian (active-branch Godunov Riemann solver)
    and the rational slip map,
  * a forward tangent-linear sensitivity marched alongside the 2N state,
  * central finite differences of the objective.

The adjoint is the O(1)-gradient channel: one backward sweep gives dJ/dC for
the scalar interfacial-drag parameter regardless of how many time steps or
cells the objective spans. The tangent-linear and finite-difference channels
exist to prove the adjoint right. Checkpointing (store every k-th state,
recompute the segments on the backward sweep) reproduces the store-all adjoint
exactly. The Godunov flux switches branches at shocks and sonic points; the
adjoint and tangent both differentiate the ACTIVE branch, so they agree with
each other to machine precision and with finite differences away from the
measure-zero branch boundaries.

The parameter is the interfacial drag C (positive by working in ln C). The
objective is a weighted holdup-misfit sum over observed (time, cell) samples.
"""

from __future__ import annotations

import math
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

GRAVITY = 9.80665


@dataclass(frozen=True)
class TwoFluidProblem:
    """Fixed geometry / fluid / controls of the sealed settling march. Only
    the interfacial drag is a calibration parameter; everything else is
    declared."""
    rho_l: float = 1000.0
    rho_g: float = 20.0
    mu_l: float = 0.01
    mu_g: float = 1e-5
    length: float = 20.0
    diameter: float = 0.05
    angle_deg: float = 90.0
    sections: int = 40
    dt: float = 0.01
    n_steps: int = 6000
    alpha0: float = 0.5
    slip0: float = 0.0

    def __post_init__(self) -> None:
        if not (self.dt > 0):
            raise ValueError("dt must be positive")
        if self.n_steps < 1:
            raise ValueError("n_steps must be >= 1")
        if self.sections < 2:
            raise ValueError("sections must be >= 2")
        if not (self.length > 0 and self.diameter > 0):
            raise ValueError("length and diameter must be positive")
        if not (self.rho_l > 0 and self.rho_g > 0):
            raise ValueError("rhoL and rhoG must be positive")
        if not (0.0 <= self.alpha0 <= 1.0):
            raise ValueError("alpha0 must lie in [0, 1]")

    @property
    def dx(self) -> float:
        return self.length / self.sections


@dataclass(frozen=True)
class Observation:
    """Holdup readings at a subset of (step index, cell) samples. Step 0 is
    the initial state."""
    indices: Tuple[int, ...]
    cells: Tuple[int, ...]
    alpha: Tuple[float, ...]
    weight: float = 1.0

    def __post_init__(self) -> None:
        if not (len(self.indices) == len(self.cells) == len(self.alpha)):
            raise ValueError("indices, cells and alpha must have equal length")
        if not self.indices:
            raise ValueError("observation must not be empty")
        if not (self.weight > 0):
            raise ValueError("weight must be positive")
        for a in self.alpha:
            if not math.isfinite(a):
                raise ValueError("observation alpha must be finite")


@dataclass
class AdjointResult:
    drag: float
    grad_c: float          # dJ/dC
    grad_log: float        # dJ/d(ln C) = C * dJ/dC
    cost: float
    n_states: int
    checkpoints: int       # how many states the backward sweep stored


# ---- Godunov quadratic Riemann solver + active-branch Jacobian --------------

def godunov(qa: float, qb: float, qc: float, aL: float, aR: float) -> float:
    """Exact Godunov flux of F(a) = qa a^2 + qb a + qc. Mirrors the C++ kernel
    and tests/check_two_fluid.py cell for cell."""
    fL = qc + aL * (qb + qa * aL)
    if aL == aR:
        return fL
    fR = qc + aR * (qb + qa * aR)
    shock = (aL > aR) if qa > 0 else (aL < aR)
    if shock:
        s = (fL - fR) / (aL - aR)
        return fL if s >= 0 else fR
    if qa == 0:
        return fL if qb >= 0 else fR
    sonic = -qb / (2.0 * qa)
    lo, hi = min(aL, aR), max(aL, aR)
    a = min(max(sonic, lo), hi)
    return qc + a * (qb + qa * a)


def godunov_jac(qa: float, qb: float, qc: float, aL: float, aR: float
                ) -> Tuple[float, float, float, float, float, float]:
    """Flux value and its active-branch partials (d/daL, d/daR, d/dqa, d/dqb,
    d/dqc). F(a) = qa a^2 + qb a + qc has dF/da = qb + 2 qa a, dF/dqa = a^2,
    dF/dqb = a, dF/dqc = 1. The branch selection is locally constant, so the
    derivative is that of the selected flux expression; at the sonic clamp
    boundary the interior and clamped expressions are C1-continuous."""
    fL = qc + aL * (qb + qa * aL)
    fR = qc + aR * (qb + qa * aR)
    if aL == aR:
        return fL, qb + 2 * qa * aL, 0.0, aL * aL, aL, 1.0
    shock = (aL > aR) if qa > 0 else (aL < aR)
    if shock:
        s = (fL - fR) / (aL - aR)
        if s >= 0:
            return fL, qb + 2 * qa * aL, 0.0, aL * aL, aL, 1.0
        return fR, 0.0, qb + 2 * qa * aR, aR * aR, aR, 1.0
    if qa == 0:
        if qb >= 0:
            return fL, qb + 2 * qa * aL, 0.0, aL * aL, aL, 1.0
        return fR, 0.0, qb + 2 * qa * aR, aR * aR, aR, 1.0
    sonic = -qb / (2.0 * qa)
    lo, hi = min(aL, aR), max(aL, aR)
    if sonic < lo:
        a = lo
        if aL < aR:
            return qc + a * (qb + qa * a), qb + 2 * qa * a, 0.0, a * a, a, 1.0
        return qc + a * (qb + qa * a), 0.0, qb + 2 * qa * a, a * a, a, 1.0
    if sonic > hi:
        a = hi
        if aL > aR:
            return qc + a * (qb + qa * a), qb + 2 * qa * a, 0.0, a * a, a, 1.0
        return qc + a * (qb + qa * a), 0.0, qb + 2 * qa * a, a * a, a, 1.0
    # interior sonic point: a = -qb/(2 qa), flux = qc - qb^2/(4 qa)
    val = qc - qb * qb / (4.0 * qa)
    return val, 0.0, 0.0, qb * qb / (4.0 * qa * qa), -qb / (2.0 * qa), 1.0


# ---- rational slip map + Jacobian -------------------------------------------

def slip_map(a: float, s_old: float, c: float, prob: TwoFluidProblem) -> float:
    """Explicit rational update of the interfacial slip. Mirrors the C++
    dynamic closure and tests/check_two_fluid.py."""
    ag = 1.0 - a
    if a <= 1e-12 or ag <= 1e-12:
        return 0.0
    rl, rg = prob.rho_l, prob.rho_g
    drho = rg - rl
    inv_dr = 1.0 / rg - 1.0 / rl
    sin_th = math.sin(prob.angle_deg * math.pi / 180.0)
    rho_m = a * rl + ag * rg
    ks = 1.0 + a * ag * drho * drho / (rl * rg)
    drive = inv_dr * rho_m * GRAVITY * sin_th
    d = 1.0 / (rl * a) + 1.0 / (rg * ag)
    rate = c * abs(s_old) * d / ks
    return (s_old / prob.dt + drive / ks) / (1.0 / prob.dt + rate)


def slip_map_jac(a: float, s_old: float, c: float, prob: TwoFluidProblem
                 ) -> Tuple[float, float, float, float]:
    """Slip value and partials (d/da, d/ds_old, d/dC)."""
    ag = 1.0 - a
    if a <= 1e-12 or ag <= 1e-12:
        return 0.0, 0.0, 0.0, 0.0
    rl, rg = prob.rho_l, prob.rho_g
    drho = rg - rl
    inv_dr = 1.0 / rg - 1.0 / rl
    sin_th = math.sin(prob.angle_deg * math.pi / 180.0)
    ks = 1.0 + a * ag * drho * drho / (rl * rg)
    drive = inv_dr * (a * rl + ag * rg) * GRAVITY * sin_th
    d = 1.0 / (rl * a) + 1.0 / (rg * ag)
    abs_s = abs(s_old)
    rate = c * abs_s * d / ks
    num = s_old / prob.dt + drive / ks
    den = 1.0 / prob.dt + rate
    s_new = num / den

    dks_da = drho * drho / (rl * rg) * (1.0 - 2.0 * a)
    ddrive_da = inv_dr * GRAVITY * sin_th * (rl - rg)
    dd_da = -1.0 / (rl * a * a) + 1.0 / (rg * ag * ag)
    ddriveks_da = (ddrive_da * ks - drive * dks_da) / (ks * ks)
    drate_da = c * abs_s * (dd_da * ks - d * dks_da) / (ks * ks)
    sgn = 1.0 if s_old > 0 else (-1.0 if s_old < 0 else 0.0)
    drate_ds = c * sgn * d / ks
    drate_dc = abs_s * d / ks

    ds_da = (ddriveks_da * den - num * drate_da) / (den * den)
    ds_ds = ((1.0 / prob.dt) * den - num * drate_ds) / (den * den)
    ds_dc = (-num * drate_dc) / (den * den)
    return s_new, ds_da, ds_ds, ds_dc


# ---- the coupled forward march ----------------------------------------------

def step(c: float, prob: TwoFluidProblem,
         alpha: Sequence[float], slip: Sequence[float]
         ) -> Tuple[List[float], List[float]]:
    n = prob.sections
    slip_new = [slip_map(alpha[i], slip[i], c, prob) for i in range(n)]
    ratio = prob.dt / prob.dx
    flux = [0.0] * (n + 1)
    for j in range(1, n):
        sf = 0.5 * (slip_new[j - 1] + slip_new[j])
        flux[j] = godunov(sf, -sf, 0.0, alpha[j - 1], alpha[j])
    alpha_new = [alpha[i] + ratio * (flux[i] - flux[i + 1]) for i in range(n)]
    return alpha_new, slip_new


def forward(c: float, prob: TwoFluidProblem
            ) -> List[Tuple[Tuple[float, ...], Tuple[float, ...]]]:
    """States s_0..s_N, each (alpha tuple, slip tuple). s_0 is the declared
    initial column."""
    if not (c > 0) or not math.isfinite(c):
        raise ValueError("interfacial drag must be positive and finite")
    n = prob.sections
    alpha = tuple([prob.alpha0] * n)
    slip = tuple([prob.slip0] * n)
    states = [(alpha, slip)]
    for _ in range(prob.n_steps):
        a, s = step(c, prob, alpha, slip)
        alpha, slip = tuple(a), tuple(s)
        states.append((alpha, slip))
    return states


# ---- objective --------------------------------------------------------------

def _obs_samples(obs: Observation, prob: TwoFluidProblem
                 ) -> List[Tuple[int, int, float]]:
    out = []
    for idx, cell, val in zip(obs.indices, obs.cells, obs.alpha):
        if not (0 <= idx <= prob.n_steps):
            raise ValueError(f"observation index {idx} outside [0, {prob.n_steps}]")
        if not (0 <= cell < prob.sections):
            raise ValueError(f"observation cell {cell} outside [0, {prob.sections})")
        out.append((idx, cell, val))
    return out


def objective(states: Sequence[Tuple[Tuple[float, ...], Tuple[float, ...]]],
              obs: Observation) -> float:
    return obs.weight * sum((states[i][0][cell] - v) ** 2
                            for i, cell, v in zip(obs.indices, obs.cells, obs.alpha))


def _seed_a(prob: TwoFluidProblem, obs: Observation, idx: int,
            alpha: Sequence[float]) -> List[float]:
    seed = [0.0] * prob.sections
    if idx in obs.indices:
        for j, cell, val in zip(obs.indices, obs.cells, obs.alpha):
            if j == idx:
                seed[cell] += 2.0 * obs.weight * (alpha[cell] - val)
    return seed


def objective_and_tangent(c: float, prob: TwoFluidProblem, obs: Observation
                          ) -> Tuple[float, float]:
    """Forward tangent-linear channel: march the 2N state and its derivative
    (dalpha/dC, dslip/dC) together. Independent of the adjoint."""
    samples = _obs_samples(obs, prob)
    n = prob.sections
    ratio = prob.dt / prob.dx
    alpha = [prob.alpha0] * n
    slip = [prob.slip0] * n
    dalpha = [0.0] * n   # initial state independent of C
    dslip = [0.0] * n
    obs_by_idx: Dict[int, List[Tuple[int, float]]] = {}
    for idx, cell, val in samples:
        obs_by_idx.setdefault(idx, []).append((cell, val))

    def accum(idx: int, a: Sequence[float], da: Sequence[float]
              ) -> Tuple[float, float]:
        cost = grad = 0.0
        for cell, val in obs_by_idx.get(idx, ()):
            resid = a[cell] - val
            cost += obs.weight * resid * resid
            grad += 2.0 * obs.weight * resid * da[cell]
        return cost, grad

    cost, grad = accum(0, alpha, dalpha)
    for m in range(1, prob.n_steps + 1):
        slip_new = [0.0] * n
        dslip_new = [0.0] * n
        for i in range(n):
            sv, da, ds, dc = slip_map_jac(alpha[i], slip[i], c, prob)
            slip_new[i] = sv
            dslip_new[i] = da * dalpha[i] + ds * dslip[i] + dc
        flux = [0.0] * (n + 1)
        dflux = [0.0] * (n + 1)
        for j in range(1, n):
            sf = 0.5 * (slip_new[j - 1] + slip_new[j])
            dsf = 0.5 * (dslip_new[j - 1] + dslip_new[j])
            _, dal, dar, dqa, dqb, _ = godunov_jac(sf, -sf, 0.0,
                                                   alpha[j - 1], alpha[j])
            flux[j] = godunov(sf, -sf, 0.0, alpha[j - 1], alpha[j])
            dflux[j] = (dal * dalpha[j - 1] + dar * dalpha[j]
                        + (dqa - dqb) * dsf)
        alpha_new = [alpha[i] + ratio * (flux[i] - flux[i + 1]) for i in range(n)]
        dalpha_new = [dalpha[i] + ratio * (dflux[i] - dflux[i + 1])
                      for i in range(n)]
        alpha, slip, dalpha, dslip = alpha_new, slip_new, dalpha_new, dslip_new
        cc, gg = accum(m, alpha, dalpha)
        cost += cc
        grad += gg
    return cost, grad


# ---- reverse-mode discrete adjoint with checkpointing -----------------------

def _backward_step(c: float, prob: TwoFluidProblem,
                   alpha_old: Sequence[float], slip_old: Sequence[float],
                   slip_new: Sequence[float],
                   lam_a_new: Sequence[float], lam_s_new: Sequence[float]
                   ) -> Tuple[List[float], List[float], float]:
    """Carry (lam_a_new, lam_s_new) at step m back to step m-1, accumulating
    dJ/dC. Transport reverse first (flux -> alpha_old, slip_new), then the
    slip-map reverse (slip_new -> alpha_old, slip_old, C)."""
    n = prob.sections
    ratio = prob.dt / prob.dx
    lam_a_old = list(lam_a_new)
    lam_s_acc = [0.0] * n
    for j in range(1, n):
        sf = 0.5 * (slip_new[j - 1] + slip_new[j])
        _, dal, dar, dqa, dqb, _ = godunov_jac(sf, -sf, 0.0,
                                               alpha_old[j - 1], alpha_old[j])
        lam_flux = ratio * (lam_a_new[j] - lam_a_new[j - 1])
        lam_a_old[j - 1] += lam_flux * dal
        lam_a_old[j] += lam_flux * dar
        dsf = lam_flux * (dqa - dqb)
        lam_s_acc[j - 1] += 0.5 * dsf
        lam_s_acc[j] += 0.5 * dsf
    grad = 0.0
    lam_s_old = [0.0] * n
    for i in range(n):
        _, da, ds, dc = slip_map_jac(alpha_old[i], slip_old[i], c, prob)
        ls = lam_s_new[i] + lam_s_acc[i]
        lam_a_old[i] += ls * da
        lam_s_old[i] = ls * ds
        grad += ls * dc
    return lam_a_old, lam_s_old, grad


def adjoint_gradient(c: float, prob: TwoFluidProblem, obs: Observation,
                     checkpoint: int = 1) -> AdjointResult:
    """One backward sweep gives dJ/dC. checkpoint=1 stores every state;
    checkpoint=k>1 stores every k-th state and recomputes each segment forward
    on the backward sweep. Both produce the identical gradient."""
    if checkpoint < 1:
        raise ValueError("checkpoint interval must be >= 1")
    if not (c > 0) or not math.isfinite(c):
        raise ValueError("interfacial drag must be positive and finite")
    samples = _obs_samples(obs, prob)
    n = prob.n_steps
    obs_by_idx: Dict[int, List[Tuple[int, float]]] = {}
    for idx, cell, val in samples:
        obs_by_idx.setdefault(idx, []).append((cell, val))

    def seed_a(idx: int, alpha: Sequence[float]) -> List[float]:
        seed = [0.0] * prob.sections
        for cell, val in obs_by_idx.get(idx, ()):
            seed[cell] = 2.0 * obs.weight * (alpha[cell] - val)
        return seed

    # checkpoint states: indices 0, k, 2k, ... and always n
    cps: Dict[int, Tuple[Tuple[float, ...], Tuple[float, ...]]] = {
        0: (tuple([prob.alpha0] * prob.sections),
            tuple([prob.slip0] * prob.sections))}
    alpha, slip = cps[0]
    for i in range(1, n + 1):
        a, s = step(c, prob, alpha, slip)
        alpha, slip = tuple(a), tuple(s)
        if i % checkpoint == 0 or i == n:
            cps[i] = (alpha, slip)
    if all(i in cps for i in obs_by_idx):
        cost = obs.weight * sum((cps[i][0][cell] - v) ** 2
                                for i, cell, v in samples)
    else:
        cost = objective(forward(c, prob), obs)

    grad = 0.0
    lam_a = seed_a(n, cps[n][0])
    lam_s = [0.0] * prob.sections
    bounds = sorted(cps.keys())
    for si in range(len(bounds) - 1, 0, -1):
        end = bounds[si]
        start = bounds[si - 1]
        seg: List[Tuple[Tuple[float, ...], Tuple[float, ...]]] = [cps[start]]
        aa, ss = cps[start]
        for _ in range(start + 1, end + 1):
            aaa, sss = step(c, prob, aa, ss)
            aa, ss = tuple(aaa), tuple(sss)
            seg.append((aa, ss))
        for j in range(end, start, -1):
            a_old, s_old = seg[j - 1 - start]
            _, s_new = seg[j - start]
            la, ls, gc = _backward_step(c, prob, a_old, s_old, s_new,
                                        lam_a, lam_s)
            grad += gc
            sa = seed_a(j - 1, a_old)
            lam_a = [la[i] + sa[i] for i in range(prob.sections)]
            lam_s = ls
    return AdjointResult(drag=c, grad_c=grad, grad_log=c * grad, cost=cost,
                         n_states=n + 1, checkpoints=len(cps))


# ---- calibration (Illinois secant in ln C) ----------------------------------

@dataclass
class FitResult:
    drag: float
    cost: float
    initial_cost: float
    iterations: int
    forward_evals: int
    converged: bool
    history: List[float] = field(default_factory=list)


def synthetic_observation(c_true: float, prob: TwoFluidProblem,
                          indices: Optional[Sequence[int]] = None,
                          cells: Optional[Sequence[int]] = None,
                          weight: float = 1.0) -> Observation:
    states = forward(c_true, prob)
    idx = tuple(range(prob.n_steps + 1)) if indices is None else tuple(indices)
    cel = tuple(range(prob.sections)) if cells is None else tuple(cells)
    ii, cc, aa = [], [], []
    for i in idx:
        for cell in cel:
            ii.append(i)
            cc.append(cell)
            aa.append(states[i][0][cell])
    return Observation(indices=tuple(ii), cells=tuple(cc), alpha=tuple(aa),
                       weight=weight)


def fit(obs: Observation, c_start: float, prob: TwoFluidProblem,
        max_iter: int = 100, x_tol: float = 1e-13,
        cost_tol: float = 1e-18) -> FitResult:
    """Scalar calibration in theta = ln C. The holdup-misfit has a single root
    of dJ/dtheta = 0, found with a safeguarded Illinois secant: each iteration
    costs ONE adjoint backward sweep (the O(1) gradient). The bracket is found
    by expanding from c_start in the downhill direction until the gradient
    changes sign."""
    if not (c_start > 0) or not math.isfinite(c_start):
        raise ValueError("c_start must be positive and finite")
    if max_iter < 1:
        raise ValueError("max_iter must be >= 1")

    evals = 0

    def eval_at(theta: float) -> Tuple[float, float]:
        nonlocal evals
        res = adjoint_gradient(math.exp(theta), prob, obs, checkpoint=1)
        evals += 1
        return res.grad_log, res.cost

    theta0 = math.log(c_start)
    g0, c0 = eval_at(theta0)
    best_theta, best_cost = theta0, c0
    history = [c0]
    if c0 < cost_tol:
        return FitResult(c_start, c0, c0, 0, evals, True, history)

    step_size = 0.1
    ta_, ga, ca = theta0, g0, c0
    tb_ = theta0 - math.copysign(step_size, g0)
    gb, cb = eval_at(tb_)
    bracketed = ga * gb <= 0
    for _ in range(60):
        if bracketed:
            break
        if cb < best_cost:
            best_theta, best_cost = tb_, cb
        history.append(best_cost)
        step_size *= 2.0
        ta_, ga, ca = tb_, gb, cb
        tb_ = tb_ - math.copysign(step_size, g0)
        gb, cb = eval_at(tb_)
        bracketed = ga * gb <= 0
    if not bracketed:
        return FitResult(math.exp(best_theta), best_cost, c0, evals, evals,
                         False, history)

    lo, hi = (ta_, tb_) if ta_ < tb_ else (tb_, ta_)
    g_lo, g_hi = (ga, gb) if ta_ < tb_ else (gb, ga)
    if g_lo > 0:
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
                g_lo *= 0.5
            side = 1
        else:
            lo, g_lo = theta_s, gs
            if side == -1:
                g_hi *= 0.5
            side = -1
    return FitResult(drag=math.exp(best_theta), cost=best_cost,
                     initial_cost=c0, iterations=it, forward_evals=evals,
                     converged=converged, history=history)


# ---- independent C++ forward channel ----------------------------------------

def case_text(c: float, prob: TwoFluidProblem,
              report_every: Optional[int] = None) -> str:
    re = prob.n_steps if report_every is None else report_every
    return f"""// D9 increment 4 two-fluid-adjoint verification case (generated).
fluid {{ rhoL {prob.rho_l!r}; rhoG {prob.rho_g!r}; muL {prob.mu_l!r}; muG {prob.mu_g!r}; }}
pipe {{ length {prob.length!r}; diameter {prob.diameter!r}; roughness 0;
        angleDeg {prob.angle_deg!r}; sections {prob.sections}; }}
closure {{ model dynamic; interfacialDrag {c!r}; }}
initial {{ alpha {prob.alpha0!r}; water 0; }}
boundaries {{ inlet {{ type sealed; }} outlet {{ type sealed; }} }}
time {{ dt {prob.dt!r}; endTime {prob.dt * prob.n_steps!r}; reportEvery {re}; }}
"""


def cpp_series(exe: str, c: float, prob: TwoFluidProblem,
               report_every: Optional[int] = None
               ) -> List[Tuple[Tuple[float, ...], Tuple[float, ...]]]:
    """Run the real C++ two-fluid kernel and read back (alpha, slip) per cell
    at every report frame (initial included). Independent of the Python
    forward."""
    import tempfile
    import os
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "case.dict")
        with open(path, "w") as fh:
            fh.write(case_text(c, prob, report_every))
        proc = subprocess.run([exe, path], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"two-fluid kernel rc={proc.returncode}: "
                           f"{proc.stderr}")
    frames: List[Tuple[Tuple[float, ...], Tuple[float, ...]]] = []
    cur_a: Dict[int, float] = {}
    cur_s: Dict[int, float] = {}
    have_t = False
    for line in proc.stdout.splitlines():
        tok = line.replace("[chflow]", "").split()
        if not tok:
            continue
        if tok[0] == "t":
            if have_t:
                frames.append((tuple(cur_a[i] for i in range(prob.sections)),
                               tuple(cur_s[i] for i in range(prob.sections))))
            have_t = True
            cur_a, cur_s = {}, {}
        elif tok[0] == "cell":
            i = int(tok[1])
            kv = dict(zip(tok[2::2], tok[3::2]))
            cur_a[i] = float(kv["alpha"])
            cur_s[i] = float(kv.get("s", "0"))
    if have_t:
        frames.append((tuple(cur_a[i] for i in range(prob.sections)),
                       tuple(cur_s[i] for i in range(prob.sections))))
    return frames