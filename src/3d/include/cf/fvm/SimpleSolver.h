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
    scalar pcorrTol = 0;          // max |p'|; alternative stop when the
                                  // predictor continuity residual floors out
                                  // due to under-relaxation lag (0 disables)
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
    Vector bodyForce = Vector::zero(); // acceleration (m/s2), i.e. -grad(p)/rho
};

struct SimpleResult
{
    int iterations = 0;
    scalar finalContinuityError = 0;
    scalar finalMomentumResidual = 0;
    bool converged = false;
};

/// Classic semi-discrete SIMPLE (Patankar & Spalding), structured after the
/// simpleFoam algorithm loop but implemented on ChuanFlow's own operators:
///   1. momentum:  div(phi,U) - laplacian(nu,U) = V(bodyForce - grad(p))
///   2. flux:      phi = interp(U) & Sf
///   3. pressure:  -laplacian(rAUf, p') = -div(phi)
///   4. correct:   U -= rAU grad(p'),  p += relaxP p'
/// Walls/symmetry keep phi = 0 through the vector BCs; p' inherits the
/// pressure BC structure (fixedValue -> 0 Dirichlet, else zeroGradient).
/// `turb` (null for laminar) is advanced once per outer iteration and supplies
/// the face-wise nu + nut for the momentum diffusion, so on a wall-function
/// mesh the wall shear is the law-of-the-wall value.
/// A case is converged when continuity (or, under momentum under-relaxation
/// where the predictor residual floors out, |p'|) says so *and* the velocity
/// has stopped changing at the relative level momentumTol asks for.
SimpleResult solveIncompressibleSimple(polyMesh& mesh,
                                       volVectorField& U,
                                       volScalarField& p,
                                       surfaceScalarField& phi,
                                       scalar nu,
                                       const Schemes& schemes,
                                       const SimpleControls& ctl,
                                       IEddyViscosity* turb = nullptr);

} // namespace cf
