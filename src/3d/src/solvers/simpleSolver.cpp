#include "cf/fvm/SimpleSolver.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Log.h"
#include "cf/fvm/PressureCorrection.h"

namespace cf
{

SimpleResult solveIncompressibleSimple(polyMesh& mesh,
                                       volVectorField& U,
                                       volScalarField& p,
                                       surfaceScalarField& phi,
                                       scalar nu,
                                       const Schemes& schemes,
                                       const SimpleControls& ctl)
{
    const label nC = mesh.nCells();
    const auto& vol = mesh.cellVolumes();
    SimpleResult result;

    // p' field with BCs mirroring p (fixedValue -> homogeneous Dirichlet)
    volScalarField pEqnField(mesh, "pcorr", 0.0);
    for (label patch = 0; patch < mesh.nPatches(); ++patch)
    {
        const std::string t = p.patchBcType(patch);
        pEqnField.setPatchValue(patch, t == "fixedValue" ? "fixedValue" : "zeroGradient", 0.0);
    }

    for (int iter = 1; iter <= ctl.maxIterations; ++iter)
    {
        // ---------------- momentum ----------------
        const Field<Vector> gradP = fvc::grad(p, schemes.grad);

        fvMatrix UEqn(mesh, viewOf(U), "UEqn");
        UEqn.schemes() = schemes;
        // relaxation is applied against the old field below, not against the
        // matrix' zero initial guess
        UEqn.schemes().relaxation = 1.0;
        fvm::div(UEqn, phi.values(), schemes.div);
        fvm::laplacian(UEqn, nu, schemes.laplacian, -1);
        for (label c = 0; c < nC; ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            for (int comp = 0; comp < 3; ++comp)
            {
                UEqn.addRhs(c, comp, -vol[u] * gradP[u][comp]);
            }
        }

        fvMatrix::Controls mctl{ctl.linearTol, ctl.linearMaxIter, ctl.nonOrthCorrectors};
        UEqn.solve(mctl);
        const Field<Vector> UOld = U.internalField();
        for (int comp = 0; comp < 3; ++comp)
        {
            const Field<scalar>& x = UEqn.solution()[static_cast<std::size_t>(comp)];
            for (label c = 0; c < nC; ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                U.internalField()[u][comp] =
                    UOld[u][comp] + ctl.relaxU * (x[u] - UOld[u][comp]);
            }
        }

        // single momentum matrix per component with identical coefficients
        // (constant nu). The pressure-correction scale rAU is taken from the
        // momentum diagonal WITHOUT the volume factor: the p' Poisson solve
        // and the explicit flux correction below use the same operator, so
        // any scale produces a divergence-free corrected flux; the scale only
        // determines the size of p'. Using the physical (V/aDiag) response
        // injects the full flux-closing pressure ramp directly into the p
        // update and SIMPLE oscillates; this Jacobi-style scale under-relaxes
        // the pressure implicitly and is the stable classic-SIMPLE choice.
        Field<scalar> rAU(static_cast<std::size_t>(nC));
        {
            const auto& d = UEqn.aDiag(0);
            for (label c = 0; c < nC; ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                rAU[u] = ctl.relaxU / std::max(d[u], 1e-300);
            }
        }

        // ---------------- fluxes ----------------
        updatePhiInternal(mesh, U, phi);
        phi.correctBoundaryFlux(U);

        const scalar contErr = continuityError(mesh, phi);
        const Field<scalar> divPhi = fvc::divPhi(mesh, phi.values());

        // ---------------- pressure correction ----------------
        const Field<scalar> rAUf = interpolateToFaces(mesh, rAU);

        fvMatrix pEqn(mesh, viewOf(pEqnField), "pEqn");
        pEqn.schemes() = schemes;
        pEqn.schemes().relaxation = 1.0;
        fvm::laplacianFace(pEqn, rAUf, schemes.laplacian, -1);
        for (label c = 0; c < nC; ++c)
        {
            pEqn.addRhs(c, 0, -divPhi[static_cast<std::size_t>(c)]);
        }
        pEqn.solve(fvMatrix::Controls{ctl.linearTol, ctl.linearMaxIter, ctl.nonOrthCorrectors});
        pEqnField.internalField() = pEqn.solution()[0];
        scalar pPrimeMax = 0;
        for (const scalar v : pEqnField.internalField()) pPrimeMax = std::max(pPrimeMax, std::abs(v));

        // ---------------- correct ----------------
        // Semi-discrete SIMPLE: the velocity gets the pressure-correction
        // response (damped by the same rAU scale that defines p' above, which
        // keeps the U-correction consistent with the flux correction below),
        // p is updated with explicit under-relaxation.
        const Field<Vector> gradPc = fvc::grad(pEqnField, schemes.grad);
        for (label c = 0; c < nC; ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            U.internalField()[u] = U.internalField()[u] - rAU[u] * gradPc[u];
            p.internalField()[u] += ctl.relaxP * pEqnField.internalField()[u];
        }

        // explicit flux correction consistent with the pressure equation
        // operator (phi_f -= rAUf_f * SnGrad(p'), OpenFOAM's pEqn.flux()),
        // then interpolate the cell-wise non-orthogonal correction so the
        // corrected phi is divergence-free up to the linear-solver tolerance.
        correctFluxes(mesh, pEqnField, gradPc, rAUf, schemes, phi);

        const scalar errAfter = continuityError(mesh, phi);

        log::info("SIMPLE iter ", iter, ": continuity ", contErr,
                  " (corrected ", errAfter, ", p' ", pPrimeMax, ")");
        result.iterations = iter;
        result.finalContinuityError = contErr;
        if (!std::isfinite(contErr) || !std::isfinite(pPrimeMax))
        {
            log::error("SIMPLE: diverged at iteration ", iter);
            break;
        }
        // Outer convergence: divergence of the momentum-based flux
        // predictor. With momentum under-relaxation the predictor keeps a lag
        // floor proportional to the remaining field change, so the decay of
        // |p'| is the honest pressure-convergence signal; pcorrTol lets a
        // case stop on it when continuityTol is tighter than the lag floor.
        if (contErr < ctl.continuityTol ||
            (ctl.pcorrTol > 0 && pPrimeMax < ctl.pcorrTol))
        {
            result.converged = true;
            break;
        }
    }

    return result;
}

} // namespace cf
