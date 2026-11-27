#include "cf/multiphase/EulerEuler.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "cf/common/Error.h"
#include "cf/fvm/Dpm.h"
#include "cf/fvm/PressureCorrection.h"
#include "cf/multiphase/VoF.h"

namespace cf::eulerEuler
{

Fields::Fields(const polyMesh& mesh, scalar alphaG0)
    : alphaG(mesh, "alphaG", alphaG0), UG(mesh, "UG"), UL(mesh, "UL"), p(mesh, "p")
{
}

void validateCase(const FluidProps& fluid, const DragModel& drag,
                  scalar alphaG0, int steps, scalar dt)
{
    if (!(fluid.rhoG > 0) || !(fluid.rhoL > 0) || !std::isfinite(fluid.rhoG)
        || !std::isfinite(fluid.rhoL))
        cfFatal("eulerEuler: rhoG and rhoL must be finite and positive (kg/m3)");
    if (!(fluid.muG > 0) || !(fluid.muL > 0) || !std::isfinite(fluid.muG)
        || !std::isfinite(fluid.muL))
        cfFatal("eulerEuler: muG and muL must be finite and positive (Pa.s)");
    if (!(alphaG0 >= 0.0) || !(alphaG0 <= 1.0) || !std::isfinite(alphaG0))
        cfFatal("eulerEuler: alphaG0 must lie in [0, 1]");
    if (steps < 1)
        cfFatal("eulerEuler: run steps must be at least 1");
    if (!(dt > 0) || !std::isfinite(dt))
        cfFatal("eulerEuler: dt must be finite and positive (s)");
    if (drag.model == "constant")
    {
        if (!(drag.beta > 0) || !std::isfinite(drag.beta))
            cfFatal("eulerEuler: the constant drag needs a positive beta (kg/m3/s)");
    }
    else if (drag.model == "schillerNaumann")
    {
        if (!(drag.diameter > 0) || !std::isfinite(drag.diameter))
            cfFatal("eulerEuler: the schillerNaumann drag needs a declared positive"
                    " particle diameter (m)");
    }
    else
    {
        cfFatal("eulerEuler: drag model '" + drag.model
                + "' is unknown; the kernel owns 'constant' and 'schillerNaumann'");
    }
}

Field<scalar> dragCoefficient(const Fields& f, const FluidProps& fluid,
                              const DragModel& drag)
{
    const polyMesh& mesh = f.alphaG.mesh();
    const label nC = mesh.nCells();
    Field<scalar> K(static_cast<std::size_t>(nC), 0.0);
    if (drag.model == "constant")
    {
        std::fill(K.begin(), K.end(), drag.beta);
        return K;
    }
    // Schiller-Naumann bubbles/drops: the published sphere law (Dpm.h) with
    // the number density alphaG/(pi dp^3/6) turns the single-particle force
    // into the per-volume exchange K |s| s, i.e. K = (3/4) Cd rhoL alphaG |s|/dp.
    const Field<Vector>& uG = f.UG.internalField();
    const Field<Vector>& uL = f.UL.internalField();
    const Field<scalar>& aG = f.alphaG.internalField();
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const scalar slip = mag(uG[u] - uL[u]);
        if (!(slip > 0) || !(aG[u] > 0)) continue;
        const scalar re = fluid.rhoL * slip * drag.diameter / fluid.muL;
        const scalar cd = sphereDragCoefficient(re);
        K[u] = 0.75 * cd * fluid.rhoL * aG[u] * slip / drag.diameter;
    }
    return K;
}

scalar phaseMass(const polyMesh& mesh, const Field<scalar>& alpha, scalar rho)
{
    const Field<scalar>& v = mesh.cellVolumes();
    scalar s = 0;
    for (std::size_t u = 0; u < static_cast<std::size_t>(mesh.nCells()); ++u)
        s += alpha[u] * rho * v[u];
    return s;
}

namespace
{

/// One phase's pressure-free momentum response:
///   [ddt(alpha rho) + div(alpha (U.Sf), U) - lap(alpha mu, U)] x = rhs
/// with rhs = alpha rho V/dt * U_old + alpha rho V * g. Writes the solved
/// field into `U` and reports convergence.
bool solvePhaseMomentum(polyMesh& mesh, volVectorField& U,
                        const Field<scalar>& alphaCell, scalar rho, scalar mu,
                        const Field<scalar>& alphaFace, const Vector& gravity,
                        const Schemes& schemes, scalar dt,
                        scalar linearTol, int maxIter, int nonOrth)
{
    const label nC = mesh.nCells();
    const Field<scalar>& vol = mesh.cellVolumes();
    const Field<Vector> old = U.internalField();

    fvMatrix eqn(mesh, viewOf(U), "UEqn");
    eqn.schemes() = schemes;
    eqn.schemes().relaxation = 1.0;

    // alpha-weighted backward Euler: fvm::ddt is unit-mass, the phase mass
    // alpha_k rho_k rides on it explicitly.
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const scalar a = alphaCell[u] * rho * vol[u] / dt;
        for (int comp = 0; comp < 3; ++comp)
        {
            eqn.addDiag(c, comp, a);
            eqn.addRhs(c, comp, a * old[u][comp]);
        }
    }

    // convection by the PHASE flux alpha_k (Uk.Sf)
    surfaceScalarField phiVel(mesh, "phiVel", 0.0);
    updatePhiInternal(mesh, U, phiVel);
    phiVel.correctBoundaryFlux(U);
    Field<scalar> phasePhi(static_cast<std::size_t>(mesh.nFaces()));
    for (std::size_t f = 0; f < phasePhi.size(); ++f)
        phasePhi[f] = alphaFace[f] * phiVel.values()[f];
    fvm::div(eqn, phasePhi, schemes.div);

    // diffusion of the phase: grad . (alpha_k mu_k grad Uk)
    Field<scalar> gammaCell(static_cast<std::size_t>(nC));
    for (std::size_t u = 0; u < static_cast<std::size_t>(nC); ++u)
        gammaCell[u] = alphaCell[u] * mu;
    fvm::laplacian(eqn, gammaCell, schemes.laplacian, -1);

    // gravity on the phase mass
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const Vector f = alphaCell[u] * rho * vol[u] * gravity;
        for (int comp = 0; comp < 3; ++comp)
            eqn.addRhs(c, comp, f[comp]);
    }

    const SolverReport rep =
        eqn.solve(fvMatrix::Controls{linearTol, maxIter, nonOrth});
    for (int comp = 0; comp < 3; ++comp)
    {
        const Field<scalar>& x = eqn.solution()[static_cast<std::size_t>(comp)];
        for (std::size_t u = 0; u < static_cast<std::size_t>(nC); ++u)
            U.internalField()[u][comp] = x[u];
    }
    return rep.converged;
}

/// The coupled drag pair with the pressure response, closed form per cell:
///   alphaG rhoG V/dt (UG - HbyA_G) = K V (UL - UG) - alphaG V gradP
/// and mirrored for the liquid. The responses are the PHASE MASS
/// alpha_k rho_k V (the ddt diagonal), so each phase's implicit drag pairing
/// is exactly its declared time constant tau_k = alpha_k rho_k / K and its
/// pressure response is dt/rho_k. The impulses are equal and opposite on the
/// SAME velocity pair, which is what keeps the mixture momentum free of any
/// drag contribution. A phase absent from a cell (alpha below `alphaFloor`)
/// gets zero response: its velocity is frozen there.
void applyDragPair(Fields& f, const Field<Vector>& hbyaG, const Field<Vector>& hbyaL,
                   const Field<scalar>& aG, const Field<scalar>& aL,
                   const Field<scalar>& K, const FluidProps& fluid,
                   scalar dt, const Field<Vector>& gradP)
{
    constexpr scalar alphaFloor = 1e-12;
    for (std::size_t u = 0; u < static_cast<std::size_t>(f.UG.internalField().size()); ++u)
    {
        const scalar respG = aG[u] > alphaFloor ? dt / (aG[u] * fluid.rhoG) : 0.0;
        const scalar respL = aL[u] > alphaFloor ? dt / (aL[u] * fluid.rhoL) : 0.0;
        const scalar dragG = K[u] * respG;             // dt / tau_G
        const scalar dragL = K[u] * respL;             // dt / tau_L
        const scalar denom = 1.0 + dragG + dragL;      // (1+a)(1+b) - ab
        const scalar pG = aG[u] > alphaFloor ? dt / fluid.rhoG : 0.0;
        const scalar pL = aL[u] > alphaFloor ? dt / fluid.rhoL : 0.0;
        for (int comp = 0; comp < 3; ++comp)
        {
            const scalar gR = hbyaG[u][comp] - pG * gradP[u][comp];
            const scalar lR = hbyaL[u][comp] - pL * gradP[u][comp];
            f.UG.internalField()[u][comp] =
                ((1.0 + dragL) * gR + dragG * lR) / denom;
            f.UL.internalField()[u][comp] =
                (dragL * gR + (1.0 + dragG) * lR) / denom;
        }
    }
}

} // namespace

StepReport advanceStep(polyMesh& mesh, Fields& f, const FluidProps& fluid,
                       const DragModel& drag, const Vector& gravity,
                       const Vector& gradPImposed,
                       const Schemes& schemes, const Controls& ctl)
{
    if (!(ctl.dt > 0)) cfFatal("eulerEuler: non-positive time step");
    const label nC = mesh.nCells();
    StepReport rep;

    // ---- 1. phase fraction: d(alphaG)/dt + div(alphaG UG) = 0 ------------
    // The gas is its own carrier: the same upwind M-matrix transport VOF
    // documents (VoF.h), driven by the gas velocity flux. alphaL = 1 - alphaG
    // is the identity and is never transported.
    surfaceScalarField phiG(mesh, "phiG", 0.0);
    updatePhiInternal(mesh, f.UG, phiG);
    phiG.correctBoundaryFlux(f.UG);
    vof::Controls vc;
    vc.tol = ctl.linearTol;
    vc.maxIter = ctl.maxIter;
    bool converged =
        vof::advectAlpha(f.alphaG, phiG.values(), ctl.dt, schemes, vc).converged;

    const Field<scalar> aGCell = f.alphaG.internalField();
    Field<scalar> aLCell(static_cast<std::size_t>(nC));
    for (std::size_t u = 0; u < static_cast<std::size_t>(nC); ++u)
        aLCell[u] = 1.0 - aGCell[u];
    const Field<scalar> aGFace = fvc::interpolate(f.alphaG);
    Field<scalar> aLFace(static_cast<std::size_t>(mesh.nFaces()));
    for (std::size_t u = 0; u < aLFace.size(); ++u)
        aLFace[u] = 1.0 - aGFace[u];

    // ---- 2. pressure-free momentum response of each phase ----------------
    const Field<scalar> K = dragCoefficient(f, fluid, drag);
    converged = solvePhaseMomentum(mesh, f.UG, aGCell, fluid.rhoG, fluid.muG,
                                   aGFace, gravity, schemes, ctl.dt,
                                   ctl.linearTol, ctl.maxIter, ctl.nonOrth)
                && converged;
    converged = solvePhaseMomentum(mesh, f.UL, aLCell, fluid.rhoL, fluid.muL,
                                   aLFace, gravity, schemes, ctl.dt,
                                   ctl.linearTol, ctl.maxIter, ctl.nonOrth)
                && converged;
    const Field<Vector> HbyA_G = f.UG.internalField();
    const Field<Vector> HbyA_L = f.UL.internalField();

    // The Poisson's face mobility is the same phase-mass response the drag
    // pair applies to the velocities: alpha_k dt/rho_k, interpolated to faces.
    Field<scalar> mobG(static_cast<std::size_t>(nC), ctl.dt / fluid.rhoG);
    Field<scalar> mobL(static_cast<std::size_t>(nC), ctl.dt / fluid.rhoL);
    const Field<scalar> mobGf = interpolateToFaces(mesh, mobG);
    const Field<scalar> mobLf = interpolateToFaces(mesh, mobL);
    Field<scalar> gammaFace(static_cast<std::size_t>(mesh.nFaces()));
    for (std::size_t u = 0; u < gammaFace.size(); ++u)
        gammaFace[u] = aGFace[u] * mobGf[u] + aLFace[u] * mobLf[u];

    // ---- 3. correctors: coupled drag pair + mixture-continuity Poisson ---
    Field<Vector> gradP = fvc::grad(f.p, schemes.grad);
    for (Vector& g : gradP) g = g + gradPImposed;
    surfaceScalarField phiMix(mesh, "phiMix", 0.0);
    for (int corr = 0; corr < std::max(1, ctl.correctors); ++corr)
    {
        applyDragPair(f, HbyA_G, HbyA_L, aGCell, aLCell, K, fluid, ctl.dt, gradP);

        surfaceScalarField phiGc(mesh, "phiGc", 0.0);
        updatePhiInternal(mesh, f.UG, phiGc);
        phiGc.correctBoundaryFlux(f.UG);
        surfaceScalarField phiLc(mesh, "phiLc", 0.0);
        updatePhiInternal(mesh, f.UL, phiLc);
        phiLc.correctBoundaryFlux(f.UL);
        for (std::size_t u = 0; u < phiMix.values().size(); ++u)
            phiMix.values()[u] = aGFace[u] * phiGc.values()[u]
                               + aLFace[u] * phiLc.values()[u];

        const Field<scalar> divPhi = fvc::divPhi(mesh, phiMix.values());
        fvMatrix pEqn(mesh, viewOf(f.p), "pEqn");
        pEqn.schemes() = schemes;
        pEqn.schemes().relaxation = 1.0;
        pEqn.solution()[0] = f.p.internalField();
        fvm::laplacianFace(pEqn, gammaFace, schemes.laplacian, -1);
        for (label c = 0; c < nC; ++c)
            pEqn.addRhs(c, 0, -divPhi[static_cast<std::size_t>(c)]);
        const SolverReport pRep = pEqn.solve(
            fvMatrix::Controls{ctl.linearTol, ctl.maxIter, ctl.nonOrth});
        converged = pRep.converged && converged;
        f.p.internalField() = pEqn.solution()[0];
        gradP = fvc::grad(f.p, schemes.grad);
        for (Vector& g : gradP) g = g + gradPImposed;
    }

    // the velocities the step ends on: the drag pair at the final grad(p),
    // so the reported state is the one the pressure was closed against
    applyDragPair(f, HbyA_G, HbyA_L, aGCell, aLCell, K, fluid, ctl.dt, gradP);

    // ---- 4. report --------------------------------------------------------
    rep.alphaMin = f.alphaG.min();
    rep.alphaMax = f.alphaG.max();
    rep.massG = phaseMass(mesh, f.alphaG.internalField(), fluid.rhoG);
    rep.massL = phaseMass(mesh, aLCell, fluid.rhoL);
    scalar vz = 0;
    for (const Vector& g : gradP) vz += g.z;
    rep.dpdz = vz / static_cast<scalar>(nC);
    const Vector centre = 0.5 * (mesh.minX() + mesh.maxX());
    label probe = 0;
    scalar best = std::numeric_limits<scalar>::max();
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const scalar d = mag2(mesh.cellCentres()[u] - centre);
        if (d < best) { best = d; probe = c; }
    }
    rep.UGz = f.UG.internalField()[static_cast<std::size_t>(probe)];
    rep.ULz = f.UL.internalField()[static_cast<std::size_t>(probe)];
    rep.continuityError = continuityError(mesh, phiMix);
    rep.converged = converged;
    return rep;
}

} // namespace cf::eulerEuler
