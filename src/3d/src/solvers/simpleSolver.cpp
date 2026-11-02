#include "cf/fvm/SimpleSolver.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Log.h"

namespace cf
{
namespace
{

/// phi at internal faces = interp(U) & Sf.
void updatePhiInternal(const polyMesh& mesh, const volVectorField& U,
                       surfaceScalarField& phi)
{
    Field<scalar> ow = makeInternalOwnerWeights(mesh);
    const auto& sf = mesh.faceSf();
    const auto& up = U.internalField();
    auto& vals = phi.values();
    for (label f = 0; f < mesh.nInternalFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const Vector& cO = up[static_cast<std::size_t>(mesh.owner()[u])];
        const Vector& cN = up[static_cast<std::size_t>(mesh.neighbour()[u])];
        const Vector uf = ow[u] * cO + (1.0 - ow[u]) * cN;
        vals[u] = dot(uf, sf[u]);
    }
}

} // namespace

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

    scalar scale = 0;
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        scale = std::max(scale, mesh.faceAreas()[static_cast<std::size_t>(f)]);
    }
    scale = std::max(scale, 1e-300);

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

        const Field<scalar> divPhi = fvc::divPhi(mesh, phi.values());
        scalar contErr = 0;
        for (label c = 0; c < nC; ++c)
        {
            contErr = std::max(contErr, std::abs(divPhi[static_cast<std::size_t>(c)]));
        }
        contErr /= scale;

        // ---------------- pressure correction ----------------
        Field<scalar> rAUf;
        fvc::interpolate(mesh, rAU,
                         [&mesh, &rAU](label f)
                         { return rAU[static_cast<std::size_t>(mesh.owner()[static_cast<std::size_t>(f)])]; },
                         rAUf);

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
        {
            const auto& pc = pEqnField.internalField();
            const auto& sf = mesh.faceSf();
            const auto& dcv = mesh.deltaCoeffVec();
            auto& vals = phi.values();
            Field<scalar> owPhi = makeInternalOwnerWeights(mesh);
            for (label f = 0; f < mesh.nFaces(); ++f)
            {
                const std::size_t u = static_cast<std::size_t>(f);
                const scalar d2 = dot(dcv[u], dcv[u]);
                if (d2 <= 0) continue;
                const scalar coef = rAUf[u] * dot(sf[u], dcv[u]) / d2;
                scalar corr = 0;
                if (schemes.laplacian == "corrected")
                {
                    const Vector dh = dcv[u] / std::sqrt(d2);
                    const Vector cfVec = sf[u] - dh * dot(sf[u], dh);
                    Vector gradFace;
                    if (f < mesh.nInternalFaces())
                    {
                        gradFace = owPhi[u] * gradPc[static_cast<std::size_t>(mesh.owner()[u])]
                                 + (1.0 - owPhi[u]) * gradPc[static_cast<std::size_t>(mesh.neighbour()[u])];
                    }
                    else
                    {
                        gradFace = gradPc[static_cast<std::size_t>(mesh.owner()[u])];
                    }
                    corr = rAUf[u] * dot(cfVec, gradFace);
                }
                if (f < mesh.nInternalFaces())
                {
                    const scalar dp = pc[static_cast<std::size_t>(mesh.neighbour()[u])]
                                    - pc[static_cast<std::size_t>(mesh.owner()[u])];
                    vals[u] -= coef * dp + corr;
                }
                else
                {
                    // boundary: Tf' - Tp' = (c-1)Tp' + k from the p' BC
                    scalar c = 0, k = 0;
                    pEqnField.linearise(f, c, k);
                    const scalar tp = pc[static_cast<std::size_t>(mesh.owner()[u])];
                    vals[u] -= coef * ((c - 1.0) * tp + k) + corr;
                }
            }
        }

        const Field<scalar> divAfter = fvc::divPhi(mesh, phi.values());
        scalar errAfter = 0;
        for (label c = 0; c < nC; ++c)
        {
            errAfter = std::max(errAfter, std::abs(divAfter[static_cast<std::size_t>(c)]));
        }
        errAfter /= scale;

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
