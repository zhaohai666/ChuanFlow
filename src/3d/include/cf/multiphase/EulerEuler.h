#pragma once

#include <string>

#include "cf/fvm/FvMatrix.h"

namespace cf::eulerEuler
{

/// 5.224: the 3D Euler-Euler (two-fluid) kernel. Where VOF carries one
/// mixture velocity and an interface fraction, this model carries a velocity
/// field PER PHASE -- gas `UG`, liquid `UL` -- and a phase fraction `alphaG`
/// that is advected by its own phase's flux. The liquid fraction is the
/// identity `alphaL = 1 - alphaG`; it is never transported separately, so the
/// two phases can only ever fill the cell together.
///
/// Everything is built from the project's own FVM operator layer
/// (fvMatrix / fvm::div / fvm::laplacian / fvc::grad, see
/// src/3d/src/discrete), and the pressure-velocity coupling mirrors the
/// single-phase PISO step of src/3d/src/solvers/pisoSolver.cpp: a
/// pressure-free momentum response per phase, then a Poisson for the FULL p
/// on the MIXTURE continuity `div(alphaG*UG + alphaL*UL) = 0`, then the
/// alpha-weighted pressure response back into both phase velocities. Solving
/// full p (not accumulating p') is the same collocation argument PisoSolver.h
/// documents.
///
/// Phase momentum, per unit volume, phase k in {G, L}:
///
///   alpha_k rho_k (dUk/dt + (Uk.grad)Uk)
///       = -alpha_k grad p + alpha_k rho_k g + K (U_other - Uk)
///         + div(alpha_k mu_k grad Uk)
///
/// The convective term is carried by the PHASE flux `alpha_k (Uk.Sf)` (the
/// conservative `div(alpha_k Uk Uk)` reading of the same term), the diffusion
/// by `fvm::laplacian` with the cell-centred `alpha_k mu_k` (harmonic face
/// averaging), and the time term by an alpha-weighted backward Euler.
///
/// The interphase exchange `K (U_other - Uk)` is the 3D restatement of the
/// 1D two-fluid dynamic-slip closure (architecture.md 4.3, the
/// `sigma = C s|s|` interfacial drag of src/1d/include/cf/1d/TwoFluid.h):
/// a momentum exchange per unit volume that is equal and opposite on the two
/// phases, linearised as `K * slip` with the slip treated implicitly. The
/// two declared closures are
///
///   constant          K = beta (kg/m3/s), uniform, the linear Stokes-pair
///                     drag; tau_k = alpha_k rho_k / K is each phase's own
///                     response time to the other.
///   schillerNaumann   K = (3/4) Cd(Re) rhoL alphaG |UG-UL| / dp with
///                     Re = rhoL |UG-UL| dp / muL and Cd the published
///                     sphere law `cf::sphereDragCoefficient` (Dpm.h, the
///                     same law the DPM increment owns). The |slip| factor
///                     makes it exactly the 1D quadratic form: the exchange
///                     is C_eff |s| s with C_eff = (3/4) Cd rhoL alphaG / dp.
///
/// The drag is applied as a closed-form 2x2 pair solve per cell (both new
/// velocities implicit in each other's drag), so the impulses are equal and
/// opposite at the SAME velocity pair: the mixture momentum
/// `alphaG rhoG UG + alphaL rhoL UL` changes only by the body forces and the
/// boundary fluxes, never by the drag. Each phase's row of that pair is the
/// implicit pairing identity
///
///   (alpha_k rho_k V/dt)(Uk_new - gR_k) = K V (U_other_new - Uk_new),
///
/// with gR_k the pressure-corrected momentum response, so the per-phase time
/// constant tau_k = alpha_k rho_k / K enters exactly as the drag response
/// K dt/(alpha_k rho_k) = dt/tau_k (the phase mass, not the full operator
/// diagonal: convection and diffusion ride in gR_k).
///
/// The phase fraction transport reuses `vof::advectAlpha` verbatim: upwind
/// convection of alphaG by the gas flux. With a divergence-free gas flux and
/// alpha in [0,1] at the boundaries the M-matrix keeps alpha in [0,1] with no
/// clipping (VoF.h documents the property); the kernel's cases are built so
/// the statement holds, and the checker measures it every frame.

struct FluidProps
{
    scalar rhoG = 0; // kg/m3
    scalar rhoL = 0; // kg/m3
    scalar muG = 0;  // Pa.s
    scalar muL = 0;  // Pa.s
};

struct DragModel
{
    std::string model = "constant"; // constant | schillerNaumann
    scalar beta = 0;                // constant: kg/m3/s
    scalar diameter = 0;            // schillerNaumann: m
};

struct Fields
{
    Fields(const polyMesh& mesh, scalar alphaG0);

    volScalarField alphaG; // gas volume fraction
    volVectorField UG;     // gas velocity
    volVectorField UL;     // liquid velocity
    volScalarField p;      // shared (cell-centred) pressure
};

struct Controls
{
    scalar dt = 0;
    int correctors = 2;        // pressure-velocity correctors per step
    scalar linearTol = 1e-10;
    int maxIter = 4000;
    int nonOrth = 2;
};

struct StepReport
{
    scalar alphaMin = 0;
    scalar alphaMax = 0;
    scalar massG = 0; // sum alphaG rhoG V
    scalar massL = 0; // sum alphaL rhoL V
    scalar dpdz = 0;  // volume-averaged grad(p).z
    Vector UGz = Vector::zero(); // probe cell (nearest the box centre)
    Vector ULz = Vector::zero();
    scalar continuityError = 0;
    bool converged = false;
};

/// Loud validation of the declared case constants (the loader's gate).
void validateCase(const FluidProps& fluid, const DragModel& drag,
                  scalar alphaG0, int steps, scalar dt);

/// Cell-centred drag coefficient K (kg/m3/s) from the declared closure and
/// the CURRENT slip. `constant` returns beta everywhere; `schillerNaumann`
/// returns (3/4) Cd(Re) rhoL alphaG |UG-UL| / dp, zero where the slip is zero.
Field<scalar> dragCoefficient(const Fields& f, const FluidProps& fluid,
                              const DragModel& drag);

/// One Euler step: advect alphaG with the gas flux, build each phase's
/// pressure-free momentum response, then `ctl.correctors` times the coupled
/// drag 2x2 + mixture-continuity Poisson for p, ending with the drag pair
/// evaluated at the final grad(p). `gradPImposed` is a DECLARED constant
/// pressure gradient added to the pressure response (the only way an open
/// column can hold its own mixture weight: `dpdz = -rho_mix g` supports the
/// phases so the buoyancy slip is the exact fixed point; the Poisson's p is
/// the deviation on top of it).
StepReport advanceStep(polyMesh& mesh, Fields& f, const FluidProps& fluid,
                      const DragModel& drag, const Vector& gravity,
                      const Vector& gradPImposed,
                      const Schemes& schemes, const Controls& ctl);

/// sum(alpha * rho * V) over the cell field.
scalar phaseMass(const polyMesh& mesh, const Field<scalar>& alpha, scalar rho);

} // namespace cf::eulerEuler
