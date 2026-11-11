#ifndef CF_1D_SOLVER_GATES_H
#define CF_1D_SOLVER_GATES_H

/// A7 (5.69): the solver exit-gate policy, in one place.
///
/// ChuanFlow's convergence knobs come in two layers, and the distinction is a
/// deliberate policy, not an accident:
///
///   OUTER (user-adjustable, per case):
///     `solve { tol; maxIter; relax; }` -- the steady network Newton's
///     convergence test, relative to the total injected rate. This is the
///     knob a user turns when a case will not converge.
///
///   INNER (named constants here, NOT user-adjustable):
///     the three internal guards the outer loop cannot reach. Their values
///     are part of the solver's contract -- every §5.x closed-form checker
///     (5.47, 5.49, 5.66, 5.68) reproduces the kernel's states to 1e-9..1e-11
///     relative on the strength of exactly these numbers, so changing one is
///     a published-behavior change, not a tuning knob.
///
///       steadyPipeResidualGate  -- pipe solver (network.cpp): the false
///           position/Illinois loop on a pipe's drive-vs-loss function stops
///           when |f| falls below this times the drive scale. It is the
///           per-pipe "the function is flat enough" exit.
///       steadyPipeBracketGate   -- the SAME loop's other exit: the bracket
///           has tightened to this relative width. MEASURED (5.47): on a
///           yielded pipe the bracket gate fires first -- the residual at
///           break is far above the residual gate, so the pipe's answer is
///           bracket-limited, not residual-limited.
///       steadyPipeMaxIter       -- the pipe loop's iteration cap.
///       transientNewtonGate     -- transient Newton (transient.cpp): the
///           time step's coupled row stops when the max relative residual
///           falls below this. MEASURED (5.49 and every march checker since):
///           the exit usually happens on the GATE, not on the fixed point --
///           the committed state sits a gate-sized relaxation away from the
///           exact implicit row and the march accumulates that slowly. The
///           closed-form checkers bound the drift at < 5e-9 relative.
///       transientNewtonMaxIter  -- the step's iteration cap; a case that
///           cannot reach the gate within it commits a partially relaxed
///           state (visible as rowError in the coupled apps, 5.40).
///
///   SCOPE: this header covers the single-phase steady pipe solve (including
///   binghamWallStress, whose own exit is machine precision and whose cap is
///   a safety net on the same policy) and the transient Newton. The two-phase
///   family (solvePipeFlowTwoPhase and its anonymous-namespace helpers) keeps
///   its own residual/bracket gates (1e-7 and 1e-12 relative to the upstream
///   pressure) -- a different contract, guarded by the two-phase checkers.
///
/// Changing any of these five numbers changes published behavior and must go
/// through the architecture document first. tests/check_solver_gates.py
/// asserts the values below against the documented policy, so the header and
/// the docs cannot drift apart silently.
namespace cf::gates
{
inline constexpr double steadyPipeResidualGate = 1e-13;
inline constexpr double steadyPipeBracketGate = 1e-9;
inline constexpr int steadyPipeMaxIter = 100;
inline constexpr double transientNewtonGate = 1e-11;
inline constexpr int transientNewtonMaxIter = 60;
} // namespace cf::gates

#endif // CF_1D_SOLVER_GATES_H
