#pragma once

#include "cf/common/Vector.h"
#include "cf/fvm/EddyViscosity.h"
#include "cf/fvm/FvMatrix.h"

namespace cf
{

/// Controls for the steady incompressible pressure-velocity coupling,
/// playing the role of OpenFOAM's simpleFoam fvSolution block.
struct SimpleControls
{
    int maxIterations = 200;
    scalar continuityTol = 1e-6;  // max |div(phi)| / domain scale
    scalar pcorrTol = 0;          // max pressure change between sweeps; the
                                  // honest outer stop: the predictor flux is
                                  // HbyA-only, so its divergence floors at the
                                  // HbyA/pressure-flux imbalance and the
                                  // continuityTol branch is often unreachable
                                  // (0 disables)
    scalar momentumTol = 0;       // max |U - U_old| / max|U|; required in both
                                  // stop branches (0 disables). Continuity
                                  // cannot see the momentum balance, so without
                                  // it a streamwise-uniform case -- every phi
                                  // zero -- reads continuity 0 and stops on the
                                  // first iteration, before the wall shear has
                                  // been felt
    scalar linearTol = 1e-7;
    int linearMaxIter = 2000;
    int nonOrthCorrectors = 2;
    scalar relaxU = 0.7;
    scalar relaxP = 0.3;
    scalar pRelaxStart = 0; // startup pressure-relaxation floor: the coupling
                            // starts at this factor and opens toward relaxP
                            // as the solved-pressure step shrinks relative
                            // to |p| (0 = relaxP from the first sweep)
    Vector bodyForce = Vector::zero(); // acceleration (m/s2), i.e. -grad(p)/rho
};

struct SimpleResult
{
    int iterations = 0;
    scalar finalContinuityError = 0;
    scalar finalMomentumResidual = 0;
    bool converged = false;
};

/// Steady incompressible pressure-velocity coupling with consistent momentum
/// interpolation (Rhie-Chow), structured after OpenFOAM's incompressibleFluid
/// SIMPLE but implemented on ChuanFlow's own operators:
///   1. momentum:  div(phi,U) - laplacian(nuEff,U) = V (bodyForce -
///      grad(p^n)); the solved x is the physical response to the current
///      pressure iterate
///   2. HbyA:      the pressure-free blend (1-relaxU) UOld + relaxU x,
///      constrained to keep prescribed boundary fluxes (constrainHbyA)
///   3. flux:      phi = flux(HbyA) only -- no cell pressure gradient enters
///      the predictor; net-balance scaled when no patch carries a fixed p
///      (p.needReference() / adjustPhi)
///   4. pressure:  laplacian(rAUf, pFull) = div(phi) solved on the FULL
///      pressure under the case's own p boundary conditions, with the
///      physical momentum response rAU = V/aDiag (OpenFOAM's 1/UEqn.A());
///      a p'-correction pinned to 0 on pressure-fixed patches would freeze
///      whatever error p^n carries there and the fixed point absorbs it
///      through the boundary cells
///   5. correct:   phi -= rAUf SnGrad(pFull)|Sf with the SOLVED (unrelaxed)
///      pressure, so the reported flux is divergence-free; p then relaxes
///      toward it by aP and U = HbyA - aP rAU (grad(pFull) - grad(p^n)) --
///      the deferred split of OpenFOAM's U = HbyA - rAtU grad(p_relaxed)
/// The flux closure and the velocity response share one operator, so the
/// fixed point is the discrete momentum balance
/// aDiag U - offdiag U = V bodyForce - V grad(p) with p the physical
/// pressure pinned by its own boundary conditions (measured on the
/// verification channel: inlet p converges to the analytic 1.2).
/// Startup: on streamwise-elongated domains the sweep-1 predictor imbalance
/// is closed by the local-response Poisson at K x the physical pressure
/// (measured K ~ 70 on the 100x20 channel); a large pressure step feeds that
/// overshoot back through grad(p) and the coupling oscillates with gain ~ K.
/// pRelaxStart floors the pressure relaxation while the solved-pressure step
/// is large relative to |p|; as the step shrinks the factor opens toward
/// relaxP, where the fixed-point 2-cycle rate |1 - 2 relaxP| vanishes
/// (measured: the channel's continuity decays at exactly 1 - 2 relaxP per
/// sweep).
/// `turb` (null for laminar) is advanced once per outer iteration and supplies
/// the face-wise nu + nut for the momentum diffusion, so on a wall-function
/// mesh the wall shear is the law-of-the-wall value.
/// A case is converged when the pressure has stopped moving (pcorrTol -- the
/// |dp| branch is the primary signal here) or the predictor continuity is
/// below continuityTol, *and* the velocity has stopped changing at the
/// relative level momentumTol asks for.
SimpleResult solveIncompressibleSimple(polyMesh& mesh,
                                       volVectorField& U,
                                       volScalarField& p,
                                       surfaceScalarField& phi,
                                       scalar nu,
                                       const Schemes& schemes,
                                       const SimpleControls& ctl,
                                       IEddyViscosity* turb = nullptr);

} // namespace cf
