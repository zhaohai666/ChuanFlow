#include <cmath>
#include <functional>

#include "cf/common/Error.h"
#include "cf/fvm/FvMatrix.h"

namespace cf
{
namespace
{

void gaussGrad(const polyMesh& mesh, const Field<scalar>& data,
               const std::function<scalar(label face, const Field<scalar>&)>& faceValueOf,
               const Field<scalar>& ow, Field<Vector>& gradT, bool corrected)
{
    const label nC = mesh.nCells();
    gradT.assign(static_cast<std::size_t>(nC), Vector::zero());
    const auto& vol = mesh.cellVolumes();
    const auto& sf = mesh.faceSf();

    // pass 0: V * grad = sum Tf * Sf
    Field<Vector> sum(static_cast<std::size_t>(nC), Vector::zero());
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const Vector flux = sf[u] * faceValueOf(f, data);
        sum[static_cast<std::size_t>(mesh.owner()[u])] += flux;
        if (f < mesh.nInternalFaces())
        {
            sum[static_cast<std::size_t>(mesh.neighbour()[u])] -= flux;
        }
    }
    for (label c = 0; c < nC; ++c)
    {
        gradT[static_cast<std::size_t>(c)] = sum[static_cast<std::size_t>(c)] / vol[static_cast<std::size_t>(c)];
    }

    if (!corrected) return;

    // correction passes: V*grad = sum [Tf*Sf - Cf*(grad_f . dhat)]
    // with d = faceCentre - ownerCellCentre, Cf = Sf - (Sf.dhat)dhat
    for (int pass = 0; pass < 2; ++pass)
    {
        std::vector<Vector> gradFace(static_cast<std::size_t>(mesh.nFaces()));
        for (label f = 0; f < mesh.nFaces(); ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            if (f < mesh.nInternalFaces())
            {
                const Vector gp = gradT[static_cast<std::size_t>(mesh.owner()[u])];
                const Vector gn = gradT[static_cast<std::size_t>(mesh.neighbour()[u])];
                gradFace[u] = ow[u] * gp + (1.0 - ow[u]) * gn;
            }
            else
            {
                gradFace[u] = gradT[static_cast<std::size_t>(mesh.owner()[u])];
            }
        }

        sum.assign(static_cast<std::size_t>(nC), Vector::zero());
        for (label f = 0; f < mesh.nFaces(); ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            const Vector sfv = sf[u];
            const Vector d = mesh.faceCentres()[u] - mesh.cellCentres()[static_cast<std::size_t>(mesh.owner()[u])];
            const scalar dm = mag(d);
            if (dm <= 0) cfFatal("grad: degenerate face-cell distance");
            const Vector dh = d / dm;
            const Vector cfVec = sfv - dh * dot(sfv, dh);
            const Vector flux = sfv * faceValueOf(f, data) - cfVec * dot(gradFace[u], dh);
            sum[static_cast<std::size_t>(mesh.owner()[u])] += flux;
            if (f < mesh.nInternalFaces())
            {
                sum[static_cast<std::size_t>(mesh.neighbour()[u])] -= flux;
            }
        }
        for (label c = 0; c < nC; ++c)
        {
            gradT[static_cast<std::size_t>(c)] = sum[static_cast<std::size_t>(c)] / vol[static_cast<std::size_t>(c)];
        }
    }
}

} // namespace

namespace fvc
{

Field<Vector> grad(const volScalarField& T, const std::string& scheme)
{
    const polyMesh& mesh = T.mesh();
    Field<scalar> ow = makeInternalOwnerWeights(mesh);
    Field<Vector> out;
    gaussGrad(mesh, T.internalField(),
              [&T](label f, const Field<scalar>&) { return T.faceValue(f); },
              ow, out, scheme == "corrected");
    return out;
}

Field<Vector> grad(const polyMesh& mesh, const Field<scalar>& data,
                   const std::function<scalar(label face, const Field<scalar>&)>& faceValueOf,
                   const std::string& scheme)
{
    Field<scalar> ow = makeInternalOwnerWeights(mesh);
    Field<Vector> out;
    gaussGrad(mesh, data, faceValueOf, ow, out, scheme == "corrected");
    return out;
}

} // namespace fvc

} // namespace cf
