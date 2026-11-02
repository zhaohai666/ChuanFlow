#include "cf/fvm/Fields.h"

#include <algorithm>
#include <cmath>
#include <map>

#include "cf/common/Error.h"

namespace cf
{
namespace
{

Field<label> buildFacePatchMap(const polyMesh& mesh)
{
    Field<label> map(static_cast<std::size_t>(mesh.nFaces()), -1);
    const auto& patches = mesh.patches();
    for (label p = 0; p < static_cast<label>(patches.size()); ++p)
    {
        const auto& pa = patches[static_cast<std::size_t>(p)];
        for (label f = pa.start; f < pa.start + pa.size; ++f)
        {
            map[static_cast<std::size_t>(f)] = p;
        }
    }
    return map;
}

/// unit face normal (outward of owner) and normal distance owner-centre->face
scalar normalDistance(const polyMesh& mesh, label face)
{
    const Vector& sf = mesh.faceSf()[static_cast<std::size_t>(face)];
    const scalar a = mag(sf);
    if (a <= 0) cfFatal("fields: zero-area face");
    const Vector n = sf / a;
    return dot(mesh.deltaCoeffVec()[static_cast<std::size_t>(face)], n);
}

int dominantAxis(const Vector& n)
{
    const scalar ax = std::abs(n.x), ay = std::abs(n.y), az = std::abs(n.z);
    if (ax >= ay && ax >= az)
    {
        if (ax < 1.0 - 1e-6) cfFatal("fields: symmetry normal not axis aligned");
        return 0;
    }
    if (ay >= az)
    {
        if (ay < 1.0 - 1e-6) cfFatal("fields: symmetry normal not axis aligned");
        return 1;
    }
    if (az < 1.0 - 1e-6) cfFatal("fields: symmetry normal not axis aligned");
    return 2;
}

} // namespace

// ------------------------------------------------------------------
// volScalarField
// ------------------------------------------------------------------

volScalarField::volScalarField(const polyMesh& mesh, const std::string& name, scalar init)
    : mesh_(mesh), name_(name), internal_(static_cast<std::size_t>(mesh.nCells()), init),
      bc_(static_cast<std::size_t>(mesh.nPatches()))
{
}

namespace
{
// shared helper: boundary faces -> patch id via a cached map
std::map<const polyMesh*, Field<label>> g_facePatchCache;

const Field<label>& facePatchMap(const polyMesh& mesh)
{
    auto it = g_facePatchCache.find(&mesh);
    if (it == g_facePatchCache.end())
    {
        it = g_facePatchCache.emplace(&mesh, buildFacePatchMap(mesh)).first;
    }
    return it->second;
}
} // namespace

label volScalarField::patchIdOfFace(label face) const
{
    return facePatchMap(mesh_)[static_cast<std::size_t>(face)];
}

void volScalarField::setPatchValue(label patchId, const std::string& bcType, scalar value)
{
    auto& bc = bc_[static_cast<std::size_t>(patchId)];
    bc.type = bcType;
    bc.value = value;
}

std::string volScalarField::patchBcType(label patchId) const
{
    return bc_[static_cast<std::size_t>(patchId)].type;
}

scalar volScalarField::faceValue(label face) const
{
    const std::size_t u = static_cast<std::size_t>(face);
    const label o = mesh_.owner()[u];
    const scalar tp = internal_[static_cast<std::size_t>(o)];
    if (face < mesh_.nInternalFaces())
    {
        const std::size_t n = static_cast<std::size_t>(mesh_.neighbour()[u]);
        const Vector& d = mesh_.deltaCoeffVec()[u];
        const scalar d2 = dot(d, d);
        const scalar w = d2 > 0
            ? std::min(1.0, std::max(0.0, dot(mesh_.faceCentres()[u] -
                          mesh_.cellCentres()[static_cast<std::size_t>(o)], d) / d2))
            : 0.5;
        return w * tp + (1.0 - w) * internal_[n];
    }
    scalar c = 0, k = 0;
    linearise(face, c, k);
    return c * tp + k;
}

void volScalarField::linearise(label face, scalar& coeff, scalar& konst) const
{
    if (face < mesh_.nInternalFaces())
    {
        cfFatal("volScalarField::linearise called on internal face");
    }
    const label p = patchIdOfFace(face);
    const auto& bc = bc_[static_cast<std::size_t>(p)];
    if (bc.type == "fixedValue")
    {
        coeff = 0;
        konst = bc.value;
    }
    else if (bc.type == "fixedGradient")
    {
        coeff = 1;
        konst = bc.value * normalDistance(mesh_, face);
    }
    else if (bc.type == "zeroGradient" || bc.type == "symmetry" || bc.type == "outlet")
    {
        coeff = 1;
        konst = 0;
    }
    else
    {
        cfFatal("volScalarField: unknown patch bc type '" + bc.type + "'");
    }
}

scalar volScalarField::min() const
{
    return *std::min_element(internal_.begin(), internal_.end());
}

scalar volScalarField::max() const
{
    return *std::max_element(internal_.begin(), internal_.end());
}

// ------------------------------------------------------------------
// volVectorField
// ------------------------------------------------------------------

volVectorField::volVectorField(const polyMesh& mesh, const std::string& name, const Vector& init)
    : mesh_(mesh), name_(name),
      internal_(static_cast<std::size_t>(mesh.nCells()), init),
      bc_(static_cast<std::size_t>(mesh.nPatches()))
{
}

label volVectorField::patchIdOfFace(label face) const
{
    return facePatchMap(mesh_)[static_cast<std::size_t>(face)];
}

void volVectorField::setPatchValue(label patchId, const std::string& bcType, const Vector& value)
{
    auto& bc = bc_[static_cast<std::size_t>(patchId)];
    bc.type = bcType;
    bc.value = value;
}

std::string volVectorField::patchBcType(label patchId) const
{
    return bc_[static_cast<std::size_t>(patchId)].type;
}

Vector volVectorField::faceValue(label face) const
{
    const label o = mesh_.owner()[static_cast<std::size_t>(face)];
    const Vector tp = internal_[static_cast<std::size_t>(o)];
    Vector out;
    for (int c = 0; c < 3; ++c)
    {
        scalar coeff = 0, konst = 0;
        linearise(face, c, coeff, konst);
        out[c] = coeff * tp[c] + konst;
    }
    return out;
}

void volVectorField::linearise(label face, int comp, scalar& coeff, scalar& konst) const
{
    if (face < mesh_.nInternalFaces())
    {
        cfFatal("volVectorField::linearise called on internal face");
    }
    const label p = patchIdOfFace(face);
    const auto& bc = bc_[static_cast<std::size_t>(p)];
    if (bc.type == "fixedValue")
    {
        coeff = 0;
        konst = bc.value[comp];
    }
    else if (bc.type == "fixedGradient")
    {
        coeff = 1;
        konst = bc.value[comp] * normalDistance(mesh_, face);
    }
    else if (bc.type == "zeroGradient" || bc.type == "outlet")
    {
        coeff = 1;
        konst = 0;
    }
    else if (bc.type == "symmetry" || bc.type == "wall")
    {
        const Vector& sf = mesh_.faceSf()[static_cast<std::size_t>(face)];
        const int k = dominantAxis(sf / mag(sf));
        coeff = (comp == k) ? 0.0 : 1.0;
        konst = 0;
    }
    else
    {
        cfFatal("volVectorField: unknown patch bc type '" + bc.type + "'");
    }
}

// ------------------------------------------------------------------
// surfaceScalarField
// ------------------------------------------------------------------

surfaceScalarField::surfaceScalarField(const polyMesh& mesh, const std::string& name, scalar init)
    : mesh_(mesh), name_(name), values_(static_cast<std::size_t>(mesh.nFaces()), init)
{
}

void surfaceScalarField::correctBoundaryFlux(const volVectorField& U)
{
    for (label f = mesh_.nInternalFaces(); f < mesh_.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        if (bcIsSurfaceEmpty(U.patchBcType(U.patchIdOfFace(f))))
        {
            values_[u] = 0;
            continue;
        }
        values_[u] = dot(U.faceValue(f), mesh_.faceSf()[u]);
    }
}

void surfaceScalarField::fixBoundaryVelocity(const Vector& U0)
{
    for (label f = mesh_.nInternalFaces(); f < mesh_.nFaces(); ++f)
    {
        values_[static_cast<std::size_t>(f)] = dot(U0, mesh_.faceSf()[static_cast<std::size_t>(f)]);
    }
}

} // namespace cf
