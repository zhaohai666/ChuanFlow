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
    scalar scale = 0;
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        scale = std::max(scale, mesh.faceAreas()[static_cast<std::size_t>(f)]);
    }
    scale = std::max(scale, 1e-300);
    const Field<scalar> divPhi = fvc::divPhi(mesh, phi.values());
    scalar err = 0;
    for (label c = 0; c < mesh.nCells(); ++c)
    {
        err = std::max(err, std::abs(divPhi[static_cast<std::size_t>(c)]));
    }
    return err / scale;
}

} // namespace cf
