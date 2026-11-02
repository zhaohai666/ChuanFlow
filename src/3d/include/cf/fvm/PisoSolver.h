#pragma once

#include "cf/fvm/FvMatrix.h"

namespace cf
{

/// Controls for one step of the transient incompressible coupling; plays the
/// role of OpenFOAM's pisoFoam fvSolution + controlDict pair.
struct PisoControls
{
    scalar dt = 0.01;
    int correctors = 2;             // pressure-velocity correctors per step
    scalar continuityTol = 1e-10;   // stop correcting below this max|div(phi)|
    scalar linearTol = 1e-8;
    int linearMaxIter = 2000;
    int nonOrthCorrectors = 2;
    Vector bodyForce = Vector::zero(); // acceleration (m/s2), i.e. -grad(p)/rho
};

struct PisoStepResult
{
    int correctorsUsed = 0;
    scalar predictedContinuityError = 0; // before the correction
    scalar continuityError = 0;          // after the last correction
    scalar pressureCorrectionMax = 0;
    bool converged = false;              // continuityTol reached
};

/// Transient incompressible pressure-velocity coupling (PISO) for
///   dU/dt + div(phi,U) - laplacian(nu,U) = -grad(p) + f,   div(U) = 0
/// assembled on ChuanFlow's own operators, structured after pisoFoam:
///   1. momentum once per step, with the previous step's flux and p^n;
///   2. HbyA = rAU*(H(U) + V grad(p^n))/V, the pressure-free momentum response;
///   3. per corrector: flux predictor from HbyA, an implicit Poisson for the
///      FULL p with the face-interpolated response rAUf, then
///      U = HbyA - rAU grad(p) and phi -= rAUf SnGrad(p)|Sf| (pEqn.flux()).
/// rAU = V/aDiag is the physical momentum response. Because fvm::ddt puts V/dt
/// on the momentum diagonal, rAU <= dt bounds every correction, so no
/// under-relaxation is needed (unlike steady SIMPLE, which is Jacobi-scaled and
/// lag-limited - see SimpleSolver.h). Time accuracy is first order (Euler ddt);
/// the corrector loop removes the splitting error, after which the corrected
/// flux is divergence-free to the linear-solver tolerance.
/// Solving for p rather than accumulating a correction p' is what keeps this
/// collocated arrangement usable: p' accumulation feeds a bell-shaped pressure
/// mode that the interpolated gradient is blind to, while here the pressure
/// enters the flux only through the compact SnGrad - the operator the Poisson
/// equation inverts - so the mode cannot be pumped.
/// Measured limit (tests/cases/startedChannel): the step is exact and steady for
/// a diffusion number nu*dt/min(cell edge)^2 <= 0.5; above it a streamwise
/// pressure ramp plus cross-flow roll - the Poisson equation's weakest mode on
/// a stretched mesh, damped only by viscosity - grows per step out of roundoff
/// until the linear solves break down near nu*dt/dy^2 = 1. Keep the diffusion
/// number at or below 0.5; the app reports the breakdown instead of hiding it.
PisoStepResult solveTransientIncompressiblePiso(polyMesh& mesh,
                                                volVectorField& U,
                                                volScalarField& p,
                                                surfaceScalarField& phi,
                                                scalar nu,
                                                const Schemes& schemes,
                                                const PisoControls& ctl);

} // namespace cf
