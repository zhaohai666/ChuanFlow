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
    // HbyA is rebuilt from the corrected velocity on every later corrector,
    // and that rebuild must use the pressure gradient the momentum equation
    // was ASSEMBLED with (p^n), not the p just solved. Keep a copy: gradP is
    // overwritten with the current p inside the corrector loop.
    const Field<Vector> gradPn = gradP;

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
            // F12 cavitation: the pressure correction sees the mixture
            // density through the response. Null for every other case.
            if (ctl.cellInvRho) rAU[u] *= (*ctl.cellInvRho)[u];
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
    // H(U) = rhs - A_off*U is extensive; rAU = V/aDiag puts it back on the
    // velocity. Pass the CURRENT velocity in: the first call uses the momentum
    // solution, the later calls (one per additional corrector) use the velocity
    // the previous corrector corrected. Rebuilding H from that field is the
    // standard PISO neighbour correction - the neighbours' just-corrected
    // velocities change this cell's pressure-free response, which is exactly
    // the coupling a second corrector is for. gradPn is p^n's gradient, the one
    // the momentum rhs was assembled with, so it must not be refreshed from
    // the p solved inside the loop.
    auto rebuildHbyA = [&](const Field<Vector>& uField)
    {
        Field<scalar> h(static_cast<std::size_t>(nC));
        Field<scalar> uComp(static_cast<std::size_t>(nC));
        for (int comp = 0; comp < 3; ++comp)
        {
            for (label c = 0; c < nC; ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                uComp[u] = uField[u][comp];
            }
            UEqn.explicitOperator(uComp, comp, h);
            for (label c = 0; c < nC; ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                HbyA.internalField()[u][comp] =
                    rAU[u] * (h[u] + vol[u] * gradPn[u][comp]) / vol[u];
            }
        }
        // constrainHbyA: where the velocity is prescribed the response keeps
        // that normal flux, elsewhere it stays a zero-gradient field.
        for (label patch = 0; patch < mesh.nPatches(); ++patch)
        {
            const bool fixed = U.patchBcType(patch) == "fixedValue";
            HbyA.setPatchValue(patch, fixed ? "fixedValue" : "zeroGradient",
                               fixed ? U.patchValue(patch) : Vector::zero());
        }
        // Same faces as the velocity: a boundary prescribed as a profile, not
        // as one number, has to be prescribed that way to the response too.
        HbyA.adoptFaceValuesFrom(U);
    };
    rebuildHbyA(U.internalField());

    const Field<scalar> rAUf = interpolateToFaces(mesh, rAU);

    // ---------------- correctors ----------------
    // Every corrector re-predicts the pressure-free flux from the CURRENT
    // HbyA and solves the FULL p there (PisoSolver.h: solving p rather than
    // accumulating p' is the documented design). The closed flux is never
    // carried into the next full-p equation - ac1b343 tried that and made the
    // next solve invert L(p) = -div(phi_closed) ~ 0, collapsing p to the
    // gauge. What DOES change between correctors is HbyA: rebuilding H(U)
    // from the velocity the previous corrector corrected is the standard
    // PISO neighbour correction, so later correctors lower the time-splitting
    // error rather than re-solving the same fixed point. Non-orthogonal
    // correction work stays inside the Poisson solve (ctl.nonOrthCorrectors).
    for (int corr = 1; corr <= ctl.correctors; ++corr)
    {
        if (corr > 1)
        {
            // neighbour correction: the previous corrector's U now carries
            // the neighbours' pressure corrections in the off-diagonal
            rebuildHbyA(U.internalField());
        }
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
        res.correctorResiduals.push_back(res.continuityError);
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
