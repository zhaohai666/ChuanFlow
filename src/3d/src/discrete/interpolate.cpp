#include <algorithm>
#include <functional>

#include "cf/common/Error.h"
#include "cf/fvm/FvMatrix.h"

namespace cf
{

namespace
{

/// Bounded linear interpolation weight of the owner cell at an internal face.
scalar ownerWeight(const polyMesh& mesh, label face)
{
    const auto& mf = mesh.faceCentres()[static_cast<std::size_t>(face)];
    const auto& cp = mesh.cellCentres()[static_cast<std::size_t>(mesh.owner()[static_cast<std::size_t>(face)])];
    const auto& d = mesh.deltaCoeffVec()[static_cast<std::size_t>(face)];
    const scalar d2 = dot(d, d);
    if (d2 <= 0) cfFatal("interpolate: zero face distance");
    const scalar w = dot(mf - cp, d) / d2;
    return std::min(1.0, std::max(0.0, w)); // bounded linear
}

} // namespace

Field<scalar> makeInternalOwnerWeights(const polyMesh& mesh)
{
    Field<scalar> w(static_cast<std::size_t>(std::max<label>(1, mesh.nInternalFaces())), 0.5);
    for (label f = 0; f < mesh.nInternalFaces(); ++f)
    {
        w[static_cast<std::size_t>(f)] = ownerWeight(mesh, f);
    }
    return w;
}

/// Face-value function for the unknown of an equation: internal faces use
/// bounded linear interpolation of the passed cell data, boundary faces use
/// the BC linearisation against the same cell data.
std::function<scalar(label face, const Field<scalar>&)> makeFaceValueFn(
    const polyMesh& mesh, const FieldView& view, int comp)
{
    Field<scalar> ow = makeInternalOwnerWeights(mesh);
    const Field<label>& owner = mesh.owner();
    const Field<label>& nbr = mesh.neighbour();
    const label nInt = mesh.nInternalFaces();
    return [ow = std::move(ow), &owner, &nbr, nInt, v = view, comp](
               label face, const Field<scalar>& data) -> scalar
    {
        const std::size_t o = static_cast<std::size_t>(owner[static_cast<std::size_t>(face)]);
        if (face < nInt)
        {
            const std::size_t n = static_cast<std::size_t>(nbr[static_cast<std::size_t>(face)]);
            const scalar w = ow[static_cast<std::size_t>(face)];
            return w * data[o] + (1.0 - w) * data[n];
        }
        scalar c = 0, k = 0;
        v.lin(face, comp, c, k);
        return c * data[o] + k;
    };
}

namespace fvc
{

void interpolate(const polyMesh& mesh, const Field<scalar>& cellData,
                 const std::function<scalar(label face)>& faceValueOf,
                 Field<scalar>& faceData)
{
    faceData.resize(static_cast<std::size_t>(mesh.nFaces()));
    const Field<scalar> ow = makeInternalOwnerWeights(mesh);
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        if (f < mesh.nInternalFaces())
        {
            const scalar tp = cellData[static_cast<std::size_t>(mesh.owner()[u])];
            const scalar tn = cellData[static_cast<std::size_t>(mesh.neighbour()[u])];
            faceData[u] = ow[u] * tp + (1.0 - ow[u]) * tn;
        }
        else
        {
            faceData[u] = faceValueOf(f);
        }
    }
}

Field<scalar> interpolate(const volScalarField& T)
{
    Field<scalar> out;
    interpolate(T.mesh(), T.internalField(),
                [&T](label f) { return T.faceValue(f); }, out);
    return out;
}

Field<scalar> div(const volScalarField& T, const Field<scalar>& phi)
{
    const polyMesh& mesh = T.mesh();
    Field<scalar> out(static_cast<std::size_t>(mesh.nCells()), 0.0);
    Field<scalar> ow = makeInternalOwnerWeights(mesh);
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        scalar tf;
        if (f < mesh.nInternalFaces())
        {
            const scalar tp = T.internalField()[static_cast<std::size_t>(mesh.owner()[u])];
            const scalar tn = T.internalField()[static_cast<std::size_t>(mesh.neighbour()[u])];
            tf = ow[u] * tp + (1.0 - ow[u]) * tn;
        }
        else
        {
            tf = T.faceValue(f);
        }
        out[static_cast<std::size_t>(mesh.owner()[u])] += phi[u] * tf;
        if (f < mesh.nInternalFaces())
        {
            out[static_cast<std::size_t>(mesh.neighbour()[u])] -= phi[u] * tf;
        }
    }
    return out;
}

Field<scalar> divPhi(const polyMesh& mesh, const Field<scalar>& phi)
{
    Field<scalar> out(static_cast<std::size_t>(mesh.nCells()), 0.0);
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        out[static_cast<std::size_t>(mesh.owner()[u])] += phi[u];
        if (f < mesh.nInternalFaces())
        {
            out[static_cast<std::size_t>(mesh.neighbour()[u])] -= phi[u];
        }
    }
    return out;
}

} // namespace fvc

} // namespace cf
