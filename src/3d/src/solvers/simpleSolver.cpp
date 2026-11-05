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
                                       IEddyViscosity* turb,
                                       IEnergyCoupling* energy)
{
    const label nC = mesh.nCells();
    const auto& vol = mesh.cellVolumes();
    SimpleResult result;

    scalar bestContErr = std::numeric_limits<scalar>::infinity();

    for (int iter = 1; iter <= ctl.maxIterations; ++iter)
    {
        if (turb) turb->correctTurbulence(0);
        const Field<Vector>* buoy = energy ? energy->buoyancyForce() : nullptr;

        // ---------------- momentum ----------------
        // deferred pressure gradient: the solve carries -V grad(p^n) so x is
        // the physical response to the current pressure iterate, while the
        // matrix source stays pressure-free and the gradient is added back
        // into HbyA below -- the velocity update then only sees grad(p_new) -
        // grad(p^n) (OpenFOAM's solve(UEqn == -fvc::grad(p))). Both halves
        // are needed: with the gradient in neither, the fixed point has
        // x = U - (V/aDiag) grad(p), the predictor flux carries a permanent
        // div(interp(x)) and the pressure absorbs it at the slow-mode scale
        // (measured: predictor continuity frozen at 2.4e-2, p ~ 700x
        // physical); with it only in the solve the fixed point is
        // double-forced instead.
        const Field<Vector> gradPn = fvc::grad(p, schemes.grad);
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
            // the buoyancy of the previous sweep's temperature, so the force and
            // the field that produced it belong to the same sweep
            const Vector rhs =
                vol[u] * (ctl.bodyForce + (buoy ? (*buoy)[u] : Vector::zero()) - gradPn[u]);
            for (int comp = 0; comp < 3; ++comp)
            {
                UEqn.addRhs(c, comp, rhs[comp]);
            }
        }

        fvMatrix::Controls mctl{ctl.linearTol, ctl.linearMaxIter, ctl.nonOrthCorrectors};
        UEqn.solve(mctl);
        const Field<Vector> UOld = U.internalField();

        // Physical momentum response rAUr = V/aDiag: the velocity's response
        // to a pressure gradient. It sets the pressure-Poisson scale and the
        // velocity correction below; both the flux closure and the correction
        // then use the SAME operator the momentum solve defines, which is
        // what makes the fixed point the discrete momentum balance itself.
        Field<scalar> rAUr(static_cast<std::size_t>(nC));
        {
            const Field<scalar>& d = UEqn.aDiag(0);
            for (int comp = 1; comp < 3; ++comp)
            {
                // constant nu (or isotropic nu + nut) and a cell-local mass
                // term make the components share the diagonal; assert that
                // rather than re-solving per component.
                const Field<scalar>& dc = UEqn.aDiag(comp);
                for (label c = 0; c < nC; ++c)
                {
                    const std::size_t u = static_cast<std::size_t>(c);
                    if (std::abs(d[u] - dc[u]) > 1e-10 * (std::abs(d[u]) + std::abs(dc[u]) + 1.0))
                        cfFatal("SIMPLE: momentum diagonal differs between components");
                }
            }
            for (label c = 0; c < nC; ++c)
            {
                rAUr[static_cast<std::size_t>(c)] =
                    vol[static_cast<std::size_t>(c)] / std::max(d[static_cast<std::size_t>(c)], 1e-300);
            }
        }

        // ---------------- HbyA ----------------
        // the pressure-free blend of the solved field. The momentum solve
        // carried -V grad(p^n) in its rhs, so on its solution the
        // pressure-free H operator is h_free = V bodyForce + (positive
        // couplings).x = aDiag x + V grad(p^n), and Hfree = (1-relaxU) UOld
        // + relaxU x is OpenFOAM's constrainHbyA(rAU*UEqn.H(), U, p) up to
        // the deferred split of the pressure term: the solved response to
        // grad(p^n) rides inside x and is removed again by the velocity
        // correction below, which reacts to grad(p_new) - grad(p^n). At the
        // fixed point UOld = x = U and the momentum balance
        // aDiag U - offdiag U = V bodyForce - V grad(p) holds exactly,
        // whatever relaxU is.
        volVectorField HbyA(mesh, "HbyA");
        Field<Vector> Hfree(static_cast<std::size_t>(nC));
        {
            for (int comp = 0; comp < 3; ++comp)
            {
                const Field<scalar>& x = UEqn.solution()[static_cast<std::size_t>(comp)];
                for (label c = 0; c < nC; ++c)
                {
                    const std::size_t u = static_cast<std::size_t>(c);
                    Hfree[u][comp] =
                        (1.0 - ctl.relaxU) * UOld[u][comp] + ctl.relaxU * x[u];
                    HbyA.internalField()[u][comp] = Hfree[u][comp];
                }
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

        // ---------------- predictor flux ----------------
        // flux(HbyA) only, OpenFOAM's phiHbyA: no cell pressure gradient
        // enters the predictor. The pressure flux is carried exclusively by
        // the face-normal gradient of the pressure solve below, which is
        // what keeps the face velocity Rhie-Chow consistent -- a cell
        // gradient here would pair with the SnGrad closure into an
        // interpolate-vs-SnGrad mismatch that the pressure equation then
        // absorbs as a permanent forcing (measured: the outlet column of
        // the verification channel collapsed 27% in Ux under exactly that
        // forcing).
        updatePhiInternal(mesh, HbyA, phi);
        phi.correctBoundaryFlux(HbyA);

        // adjustPhi -- ONLY when the pressure has no fixed value anywhere
        // (OpenFOAM's p.needReference()): then the Poisson is null-mode
        // pinned by the global balance. With a pressure-Dirichlet patch the
        // boundary absorbs the imbalance physically, and a rescale here
        // would instead invert the free-boundary flux (measured on the
        // inlet-driven channel: outlet flux ~0 at startup, the rescale
        // factor (freeSum-netSum)/freeSum blew up and drove the outlet to
        // -Q, which the full-p Poisson turned into a 70x-physical pressure
        // and divergence).
        {
            const auto& patches = mesh.patches();
            bool needReference = true;
            for (label patch = 0; patch < mesh.nPatches(); ++patch)
            {
                if (p.patchBcType(patch) == "fixedValue") needReference = false;
            }
            if (needReference)
            {
                scalar netSum = 0, freeSum = 0;
                for (label patch = 0; patch < mesh.nPatches(); ++patch)
                {
                    const bool fixedU = U.patchBcType(patch) == "fixedValue";
                    scalar s = 0;
                    for (label f = patches[patch].start;
                         f < patches[patch].start + patches[patch].size; ++f)
                    {
                        s += phi.values()[static_cast<std::size_t>(f)];
                    }
                    netSum += s;
                    if (!fixedU) freeSum += s;
                }
                if (std::abs(freeSum) > 1e-30)
                {
                    const scalar scale = (freeSum - netSum) / freeSum;
                    for (label patch = 0; patch < mesh.nPatches(); ++patch)
                    {
                        if (U.patchBcType(patch) == "fixedValue") continue;
                        for (label f = patches[patch].start;
                             f < patches[patch].start + patches[patch].size; ++f)
                        {
                            phi.values()[static_cast<std::size_t>(f)] *= scale;
                        }
                    }
                }
            }
        }

        const scalar contErr = continuityError(mesh, phi);
        const Field<scalar> divPhi = fvc::divPhi(mesh, phi.values());

        // ---------------- pressure ----------------
        // OpenFOAM's full-p form: the Poisson solves for the WHOLE pressure
        // under the case's own p boundary conditions, not for a correction
        // pinned to 0 on pressure-fixed patches -- such a BC would freeze
        // whatever error p^n carries there, and the fixed point absorbs it
        // through the boundary cells. The rhs is the predictor imbalance, so
        // at the fixed point (pFull = p) the pressure satisfies exactly
        // L p = div(flux(U)) with the flux closed below divergence-free.
        volScalarField pFull(mesh, "pFull");
        for (label patch = 0; patch < mesh.nPatches(); ++patch)
        {
            pFull.setPatchValue(patch, p.patchBcType(patch), p.patchValue(patch));
        }
        const Field<scalar> rAUf = interpolateToFaces(mesh, rAUr);
        fvMatrix pEqn(mesh, viewOf(pFull), "pEqn");
        pEqn.schemes() = schemes;
        pEqn.schemes().relaxation = 1.0;
        fvm::laplacianFace(pEqn, rAUf, schemes.laplacian, -1);
        for (label c = 0; c < nC; ++c)
        {
            pEqn.addRhs(c, 0, -divPhi[static_cast<std::size_t>(c)]);
        }
        const SolverReport pRep = pEqn.solve(mctl);
        if (!pRep.converged) log::warn("SIMPLE: pressure solve did not converge");
        const Field<scalar> pFullNew = pEqn.solution()[0];
        for (label c = 0; c < nC; ++c)
        {
            pFull.internalField()[static_cast<std::size_t>(c)] =
                pFullNew[static_cast<std::size_t>(c)];
        }

        // ---------------- correct ----------------
        // the flux closes with the UNRELAXED solved pressure so the reported
        // phi is divergence-free; p then relaxes toward the solved pressure
        // and the velocity reacts to the pressure CHANGE at the physical
        // response rAUr, completing the deferred split (x responded to
        // grad(p^n), the correction removes grad(p^n) and installs
        // grad(p_new)). Both halves use the SAME aP.
        //
        // aP adapts to the solved-pressure step: the sweep-1 Poisson closes
        // the raw startup imbalance at the local-response scale, K x the
        // physical pressure (K ~ 70 on the 100x20 channel, the elongation
        // ratio), and a step of that size fed back through grad(p) makes the
        // coupling oscillate with gain ~ K. Keeping the step under ~15% of
        // the field scale stays inside the loop's linear regime, where it is
        // a fixed-point 2-cycle with rate |1 - 2 aP| -- so as the step
        // shrinks the factor opens toward relaxP instead of crawling at the
        // startup-safe floor. pRelaxStart = 0 recovers the constant relaxP.
        const Field<Vector> gradPFull = fvc::grad(pFull, schemes.grad);
        correctFluxes(mesh, pFull, gradPFull, rAUf, schemes, phi);

        scalar pStep = 0, pScale = 0;
        for (label c = 0; c < nC; ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            pStep = std::max(pStep, std::abs(pFullNew[u] - p.internalField()[u]));
            pScale = std::max(pScale, std::abs(p.internalField()[u]));
        }
        const scalar aP = ctl.pRelaxStart > 0
            ? std::min(ctl.relaxP,
                       std::max(ctl.pRelaxStart, 0.15 * pScale / std::max(pStep, 1e-300)))
            : ctl.relaxP;

        scalar pDeltaMax = 0;
        for (label c = 0; c < nC; ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            const scalar dp = aP * (pFullNew[u] - p.internalField()[u]);
            pDeltaMax = std::max(pDeltaMax, std::abs(dp));
            p.internalField()[u] += dp;
            U.internalField()[u] = Hfree[u]
                                 - aP * rAUr[u] * (gradPFull[u] - gradPn[u]);
        }

        // Temperature is transported by the flux this sweep closed, so the
        // buoyancy the next sweep is forced with comes from the velocity field
        // that the current sweep produced.
        if (energy) energy->correctEnergy(0);

        const scalar errAfter = continuityError(mesh, phi);

        // The fixed-point change of the whole velocity update -- momentum
        // relaxation plus the pressure response. It is the only signal here
        // that speaks to the momentum balance: continuity can be satisfied by
        // a velocity field that has nothing to do with the forces, which a
        // streamwise-uniform case makes explicit -- every phi is zero, so
        // continuity reads 0 on iteration 1 and the wall shear is never felt.
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
                  " (corrected ", errAfter, ", dp ", pDeltaMax, ", dU ", uRes, ")");
        result.iterations = iter;
        result.finalContinuityError = contErr;
        result.finalMomentumResidual = uRes;
        if (!std::isfinite(contErr) || !std::isfinite(pDeltaMax) || !fieldsFinite(U, p))
        {
            log::error("SIMPLE: diverged at iteration ", iter);
            break;
        }
        // Outer convergence. The predictor flux is built from HbyA alone and
        // the pressure-flux half is only added by the correction below, so
        // the reported predictor continuity floors at the HbyA/pressure-flux
        // mismatch (it converges with the pressure, but to a finite level:
        // ~5e-3 on the verification channel) -- continuityTol stays reachable
        // only for cases where that level sits under it. The |dp| branch is
        // therefore the honest outer stop: the pressure step decays at the
        // 2-cycle rate |1 - 2 aP| and reading under pcorrTol with the
        // velocity settled means the fixed point is reached. It needs
        // continuity to have stopped growing and is only consulted once the
        // startup transient has passed.
        // momentumTol gates both branches on the velocity having settled.
        // A warm start can begin below the floor its predictor flux settles at,
        // and then strict monotonicity is unreachable for the rest of the sweep:
        // the guard is against continuity that is still growing, and a percent
        // below the best is far enough from that to catch.
        const bool settled = contErr <= bestContErr * 1.01 + 1e-12;
        bestContErr = std::min(bestContErr, contErr);
        const bool momentumSettled = !(ctl.momentumTol > 0) || uRes <= ctl.momentumTol;
        if (momentumSettled &&
            (contErr < ctl.continuityTol ||
             (iter >= 5 && ctl.pcorrTol > 0 && pDeltaMax < ctl.pcorrTol && settled)))
        {
            result.converged = true;
            break;
        }
    }

    return result;
}

} // namespace cf
