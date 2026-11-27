#include "cf/fvm/PressureCorrection.h"

#include <algorithm>
#include <cmath>
#include <functional>

namespace cf
{

void updatePhiInternal(const polyMesh& mesh, const volVectorField& U,
                       surfaceScalarField& phi)
{
    const Field<scalar> ow = makeInternalOwnerWeights(mesh);
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

Field<scalar> interpolateToFaces(const polyMesh& mesh, const Field<scalar>& cellData)
{
    Field<scalar> faceData(static_cast<std::size_t>(mesh.nFaces()));
    fvc::interpolate(mesh, cellData,
                     [&mesh, &cellData](label f)
                     { return cellData[static_cast<std::size_t>(mesh.owner()[static_cast<std::size_t>(f)])]; },
                     faceData);
    return faceData;
}

void correctFluxes(const polyMesh& mesh,
                   const volScalarField& p,
                   const Field<Vector>& gradP,
                   const Field<scalar>& rAUf,
                   const Schemes& schemes,
                   surfaceScalarField& phi)
{
    const auto& pCell = p.internalField();
    const auto& sf = mesh.faceSf();
    const auto& dcv = mesh.deltaCoeffVec();
    auto& vals = phi.values();
    const Field<scalar> ow = makeInternalOwnerWeights(mesh);

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
                gradFace = ow[u] * gradP[static_cast<std::size_t>(mesh.owner()[u])]
                         + (1.0 - ow[u]) * gradP[static_cast<std::size_t>(mesh.neighbour()[u])];
            }
            else
            {
                gradFace = gradP[static_cast<std::size_t>(mesh.owner()[u])];
            }
            corr = rAUf[u] * dot(cfVec, gradFace);
        }
        if (f < mesh.nInternalFaces())
        {
            const scalar dp = pCell[static_cast<std::size_t>(mesh.neighbour()[u])]
                            - pCell[static_cast<std::size_t>(mesh.owner()[u])];
            vals[u] -= coef * dp + corr;
        }
        else
        {
            // boundary: Tf - Tp = (c-1)Tp + k from the pressure BC
            scalar c = 0, k = 0;
            p.linearise(f, c, k);
            const scalar tp = pCell[static_cast<std::size_t>(mesh.owner()[u])];
            vals[u] -= coef * ((c - 1.0) * tp + k) + corr;
        }
    }
}

scalar continuityError(const polyMesh& mesh, const surfaceScalarField& phi)
{
    // Per-cell continuity imbalance normalised by the cell's OWN largest face
    // area, not the global max|Sf|. The old metric max|div(phi)| / max_f|Sf|
    // divides every cell by the single largest face in the whole mesh: on a
    // graded or mixed-element mesh a small cell's real leak is diluted by a
    // big face elsewhere, so SIMPLE reads a falsely small error and converges
    // while that cell still leaks. Normalising each cell by its own face scale
    // removes that dilution -- the local estimate is never smaller than the
    // global one (local max area <= global max area), so the criterion is
    // strictly stronger. On a uniform mesh every cell's max face area equals
    // the global max, and the metric reduces exactly to the old value, keeping
    // the velocity-scale meaning of continuityTol intact.
    const auto& areas = mesh.faceAreas();
    Field<scalar> localScale(static_cast<std::size_t>(mesh.nCells()), 0.0);
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const scalar a = areas[u];
        localScale[static_cast<std::size_t>(mesh.owner()[u])] =
            std::max(localScale[static_cast<std::size_t>(mesh.owner()[u])], a);
        if (f < mesh.nInternalFaces())
        {
            localScale[static_cast<std::size_t>(mesh.neighbour()[u])] =
                std::max(localScale[static_cast<std::size_t>(mesh.neighbour()[u])], a);
        }
    }
    const Field<scalar> divPhi = fvc::divPhi(mesh, phi.values());
    scalar err = 0;
    for (label c = 0; c < mesh.nCells(); ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        err = std::max(err, std::abs(divPhi[u]) / std::max(localScale[u], 1e-300));
    }
    return err;
}

} // namespace cf
