#include "cf/fvm/SimpleSolver.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "cf/common/Log.h"
#include "cf/fvm/PressureCorrection.h"

namespace cf
{

namespace
{

/// A coupling can also fail by turning every field into NaN, which then makes
/// the residuals read as zero; only finite fields mean anything.
bool fieldsFinite(const volVectorField& U, const volScalarField& p)
{
    for (const Vector& u : U.internalField())
    {
        if (!std::isfinite(u.x) || !std::isfinite(u.y) || !std::isfinite(u.z)) return false;
    }
    for (const scalar v : p.internalField())
    {
        if (!std::isfinite(v)) return false;
    }
    return true;
}

} // namespace

SimpleResult solveIncompressibleSimple(polyMesh& mesh,
                                       volVectorField& U,
                                       volScalarField& p,
                                       surfaceScalarField& phi,
                                       scalar nu,
                                       const Schemes& schemes,
                                       const SimpleControls& ctl,
                                       IEddyViscosity* turb)
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

    scalar bestContErr = std::numeric_limits<scalar>::infinity();

    for (int iter = 1; iter <= ctl.maxIterations; ++iter)
    {
        if (turb) turb->correctTurbulence(0);

        // ---------------- momentum ----------------
        const Field<Vector> gradP = fvc::grad(p, schemes.grad);

        fvMatrix UEqn(mesh, viewOf(U), "UEqn");
        UEqn.schemes() = schemes;
        // relaxation is applied against the old field below, not against the
        // matrix' zero initial guess
        UEqn.schemes().relaxation = 1.0;
        fvm::div(UEqn, phi.values(), schemes.div);
        if (turb) fvm::laplacianFace(UEqn, turb->nuEffFace(), schemes.laplacian, -1);
        else fvm::laplacian(UEqn, nu, schemes.laplacian, -1);
        for (label c = 0; c < nC; ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            const Vector rhs = vol[u] * (ctl.bodyForce - gradP[u]);
            for (int comp = 0; comp < 3; ++comp)
            {
                UEqn.addRhs(c, comp, rhs[comp]);
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
        // (constant nu, or an isotropic nu + nut from the turbulence model).
        // The pressure-correction scale rAU is taken from the
        // momentum diagonal WITHOUT the volume factor: the p' Poisson solve
        // and the explicit flux correction below use the same operator, so
        // any scale produces a divergence-free corrected flux; the scale only
        // determines the size of p'. Using the physical (V/aDiag) response
        // injects the full flux-closing pressure ramp directly into the p
        // update and SIMPLE oscillates; this Jacobi-style scale under-relaxes
        // the pressure implicitly and is the stable classic-SIMPLE choice.
        //
        // Known limit (M3), measured rather than assumed: because p accumulates
        // on this damped scale, a case whose physics is set by the streamwise
        // pressure gradient converges its velocity field long before its
        // pressure field. Closing that gap needs a consistent momentum
        // interpolation (Rhie-Chow), not a tighter tolerance.
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

        // The fixed-point change of the whole velocity update -- momentum
        // relaxation plus the pressure-correction response. It is the only
        // signal here that speaks to the momentum balance: continuity can be
        // satisfied by a velocity field that has nothing to do with the forces,
        // which a streamwise-uniform case makes explicit -- every phi is zero,
        // so continuity reads 0 on iteration 1 and the wall shear is never felt.
        // Gating on it, a force-driven wall-function channel converges as a
        // steady problem in about 30 iterations.
        scalar uMax = 0, dUMax = 0;
        for (label c = 0; c < nC; ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            uMax = std::max(uMax, mag(U.internalField()[u]));
            dUMax = std::max(dUMax, mag(U.internalField()[u] - UOld[u]));
        }
        const scalar uRes = dUMax / std::max(uMax, 1e-300);

        log::info("SIMPLE iter ", iter, ": continuity ", contErr,
                  " (corrected ", errAfter, ", p' ", pPrimeMax, ", dU ", uRes, ")");
        result.iterations = iter;
        result.finalContinuityError = contErr;
        result.finalMomentumResidual = uRes;
        if (!std::isfinite(contErr) || !std::isfinite(pPrimeMax) || !fieldsFinite(U, p))
        {
            log::error("SIMPLE: diverged at iteration ", iter);
            break;
        }
        // Outer convergence: divergence of the momentum-based flux
        // predictor. With momentum under-relaxation the predictor keeps a lag
        // floor proportional to the remaining field change, so the decay of
        // |p'| is the honest pressure-convergence signal; pcorrTol lets a
        // case stop on it when continuityTol is tighter than the lag floor.
        // A diverging coupling can also pass through a small-p' iteration, so
        // the |p'| stop needs continuity to have stopped growing, and it is
        // only consulted once the initial transient has been resolved.
        // momentumTol then gates both branches on the velocity having settled.
        const bool settled = contErr <= bestContErr;
        bestContErr = std::min(bestContErr, contErr);
        const bool momentumSettled = !(ctl.momentumTol > 0) || uRes <= ctl.momentumTol;
        if (momentumSettled &&
            (contErr < ctl.continuityTol ||
             (iter >= 5 && ctl.pcorrTol > 0 && pPrimeMax < ctl.pcorrTol && settled)))
        {
            result.converged = true;
            break;
        }
    }

    return result;
}

} // namespace cf
