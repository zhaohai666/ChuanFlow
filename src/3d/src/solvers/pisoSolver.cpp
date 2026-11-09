#include "cf/fvm/PisoSolver.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/PressureCorrection.h"

namespace cf
{

PisoStepResult solveTransientIncompressiblePiso(polyMesh& mesh,
                                                volVectorField& U,
                                                volScalarField& p,
                                                surfaceScalarField& phi,
                                                scalar nu,
                                                const Schemes& schemes,
                                                const PisoControls& ctl,
                                                IEddyViscosity* turb,
                                                IEnergyCoupling* energy)
{
    const label nC = mesh.nCells();
    const auto& vol = mesh.cellVolumes();
    PisoStepResult res;

    if (ctl.dt <= 0) cfFatal("PISO: non-positive time step");

    if (turb) turb->correctTurbulence(ctl.dt);
    const Field<scalar>* nuEffFace = turb ? &turb->nuEffFace() : nullptr;
    // the buoyancy of the temperature the previous step left behind; it rides
    // into HbyA with the other forces, so it reaches the flux the same way
    const Field<Vector>* buoy = energy ? energy->buoyancyForce() : nullptr;

    // The pressure keeps its own BCs; the correction equation is not needed
    // here because p is solved implicitly inside the step (see below).
    Field<Vector> gradP = fvc::grad(p, schemes.grad);

    // ---------------- momentum (once per time step) ----------------
    const Field<Vector> UOld = U.internalField();
    fvMatrix UEqn(mesh, viewOf(U), "UEqn");
    UEqn.schemes() = schemes;
    UEqn.schemes().relaxation = 1.0;
    fvm::ddt(UEqn, UOld, ctl.dt);
    fvm::div(UEqn, phi.values(), schemes.div);
    if (nuEffFace) fvm::laplacianFace(UEqn, *nuEffFace, schemes.laplacian, -1);
    else fvm::laplacian(UEqn, nu, schemes.laplacian, -1);
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const Vector force =
            ctl.bodyForce + (ctl.cellForce ? (*ctl.cellForce)[u] : Vector::zero())
                          + (buoy ? (*buoy)[u] : Vector::zero());
        const Vector rhs = vol[u] * (force - gradP[u]);
        for (int comp = 0; comp < 3; ++comp)
        {
            UEqn.addRhs(c, comp, rhs[comp]);
        }
    }
    UEqn.solve(fvMatrix::Controls{ctl.linearTol, ctl.linearMaxIter, ctl.nonOrthCorrectors});
    for (int comp = 0; comp < 3; ++comp)
    {
        const Field<scalar>& x = UEqn.solution()[static_cast<std::size_t>(comp)];
        for (label c = 0; c < nC; ++c)
        {
            U.internalField()[static_cast<std::size_t>(c)][comp] = x[static_cast<std::size_t>(c)];
        }
    }

    // rAU = V/aDiag is the physical momentum response. fvm::ddt puts V/dt on
    // the diagonal, so rAU <= dt bounds every correction and no
    // under-relaxation is needed (unlike steady SIMPLE).
    Field<scalar> rAU(static_cast<std::size_t>(nC));
    {
        const Field<scalar>& d = UEqn.aDiag(0);
        for (int comp = 1; comp < 3; ++comp)
        {
            // constant nu and a cell-local mass term make the components share
            // the diagonal; assert that rather than re-solving per component.
            const Field<scalar>& dc = UEqn.aDiag(comp);
            for (label c = 0; c < nC; ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                if (std::abs(d[u] - dc[u]) > 1e-10 * (std::abs(d[u]) + std::abs(dc[u]) + 1.0))
                    cfFatal("PISO: momentum diagonal differs between components");
            }
        }
        for (label c = 0; c < nC; ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            rAU[u] = vol[u] / std::max(d[u], 1e-300);
        }
    }

    // HbyA: the momentum response WITHOUT the pressure gradient. OpenFOAM's
    // UEqn.H()*rAU(). The pressure then enters the step only through the face
    // normal gradient in the flux and the cell gradient in U, which is what
    // keeps the collocated pressure-velocity coupling consistent: the flux
    // sees SnGrad(p) on the compact owner-neighbour stencil, the same operator
    // the pressure equation inverts. Interpolating the corrected velocity
    // instead would leave the pressure mode invisible to continuity, and the
    // accumulated p then grows it (verified on the start-up channel).
    volVectorField HbyA(mesh, "HbyA");
    {
        Field<scalar> h(static_cast<std::size_t>(nC));
        for (int comp = 0; comp < 3; ++comp)
        {
            UEqn.explicitOperator(UEqn.solution()[static_cast<std::size_t>(comp)], comp, h);
            for (label c = 0; c < nC; ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                // the matrix rhs still carries -V grad(p^n): H(U) is extensive,
                // so divide by the diagonal (through rAU = V/aDiag) once.
                HbyA.internalField()[u][comp] = rAU[u] * (h[u] + vol[u] * gradP[u][comp]) / vol[u];
            }
        }
    }
    // constrainHbyA: where the velocity is prescribed the response keeps that
    // normal flux, elsewhere it stays a zero-gradient field.
    for (label patch = 0; patch < mesh.nPatches(); ++patch)
    {
        const bool fixed = U.patchBcType(patch) == "fixedValue";
        HbyA.setPatchValue(patch, fixed ? "fixedValue" : "zeroGradient",
                           fixed ? U.patchValue(patch) : Vector::zero());
    }
    // Same faces as the velocity: a boundary prescribed as a profile, not as
    // one number, has to be prescribed that way to the response too.
    HbyA.adoptFaceValuesFrom(U);

    const Field<scalar> rAUf = interpolateToFaces(mesh, rAU);

    // ---------------- correctors ----------------
    for (int corr = 1; corr <= ctl.correctors; ++corr)
    {
        // predicted flux: momentum response, pressure not yet included
        updatePhiInternal(mesh, HbyA, phi);
        phi.correctBoundaryFlux(HbyA);
        res.predictedContinuityError = continuityError(mesh, phi);

        const Field<scalar> divPhi = fvc::divPhi(mesh, phi.values());

        fvMatrix pEqn(mesh, viewOf(p), "pEqn");
        pEqn.schemes() = schemes;
        pEqn.schemes().relaxation = 1.0;
        pEqn.solution()[0] = p.internalField(); // start from the previous iterate
        fvm::laplacianFace(pEqn, rAUf, schemes.laplacian, -1);
        for (label c = 0; c < nC; ++c)
        {
            pEqn.addRhs(c, 0, -divPhi[static_cast<std::size_t>(c)]);
        }
        const Field<scalar> pBefore = p.internalField();
        const SolverReport pRep =
            pEqn.solve(fvMatrix::Controls{ctl.linearTol, ctl.linearMaxIter, ctl.nonOrthCorrectors});
        if (!pRep.converged) log::warn("PISO: pressure solve did not converge");
        p.internalField() = pEqn.solution()[0];

        scalar pDeltaMax = 0;
        for (label c = 0; c < nC; ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            pDeltaMax = std::max(pDeltaMax, std::abs(p.internalField()[u] - pBefore[u]));
        }
        res.pressureCorrectionMax = pDeltaMax;

        gradP = fvc::grad(p, schemes.grad);
        for (label c = 0; c < nC; ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            U.internalField()[u] = HbyA.internalField()[u] - rAU[u] * gradP[u];
        }

        // consistent flux: phi = HbyA & Sf - rAUf SnGrad(p)|Sf| (pEqn.flux())
        correctFluxes(mesh, p, gradP, rAUf, schemes, phi);

        res.correctorsUsed = corr;
        res.continuityError = continuityError(mesh, phi);
        if (!std::isfinite(res.continuityError) || !std::isfinite(pDeltaMax))
        {
            log::error("PISO: diverged in corrector ", corr);
            break;
        }
        if (res.continuityError < ctl.continuityTol)
        {
            res.converged = true;
            break;
        }
    }

    // The temperature rides the flux the step closed on, so the buoyancy the
    // next step is forced with is the one this step's motion produced.
    if (energy) energy->correctEnergy(ctl.dt);

    return res;
}

} // namespace cf
