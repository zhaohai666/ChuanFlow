#include "cf/fvm/PolyMesh.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"
#include "cf/common/Log.h"

namespace cf
{

polyMesh::polyMesh(Field<Vector> points,
                   std::vector<std::vector<label>> faces,
                   Field<label> owner,
                   Field<label> neighbour,
                   std::vector<Patch> patches)
    : points_(std::move(points)),
      faces_(std::move(faces)),
      owner_(std::move(owner)),
      neighbour_(std::move(neighbour)),
      patches_(std::move(patches))
{
    for (const label o : owner_) nCells_ = std::max(nCells_, o + 1);

    cells_.resize(static_cast<std::size_t>(nCells_));
    for (label f = 0; f < nFaces(); ++f)
    {
        cells_[static_cast<std::size_t>(owner_[static_cast<std::size_t>(f)])].push_back(f);
        if (f < nInternalFaces())
        {
            cells_[static_cast<std::size_t>(neighbour_[static_cast<std::size_t>(f)])].push_back(f);
        }
    }

    // face -> patch addressing, and the sanity check that patches cover all
    // boundary faces exactly once
    facePatchMap_.assign(static_cast<std::size_t>(nFaces()), -1);
    for (label i = 0; i < nPatches(); ++i)
    {
        const Patch& p = patches_[static_cast<std::size_t>(i)];
        for (label f = p.start; f < p.start + p.size; ++f)
        {
            if (f < nInternalFaces())
            {
                cfFatal("polyMesh: patch '" + p.name + "' covers internal face " + std::to_string(f));
            }
            const std::size_t u = static_cast<std::size_t>(f);
            if (facePatchMap_[u] >= 0)
            {
                cfFatal("polyMesh: boundary face " + std::to_string(f) + " in two patches");
            }
            facePatchMap_[u] = i;
        }
    }
    for (label f = nInternalFaces(); f < nFaces(); ++f)
    {
        if (facePatchMap_[static_cast<std::size_t>(f)] < 0)
        {
            cfFatal("polyMesh: boundary face " + std::to_string(f) + " is not in any patch");
        }
    }

    calcGeometry();
}

label polyMesh::patchIdByName(const std::string& name) const
{
    for (label i = 0; i < nPatches(); ++i)
    {
        if (patches_[static_cast<std::size_t>(i)].name == name) return i;
    }
    cfFatal("polyMesh: no patch named '" + name + "'");
}

std::vector<std::string> polyMesh::patchNames() const
{
    std::vector<std::string> out;
    out.reserve(patches_.size());
    for (const Patch& p : patches_) out.push_back(p.name);
    return out;
}

void polyMesh::calcGeometry() const
{
    if (geomCalculated_) return;
    geomCalculated_ = true;

    const label nF = nFaces();

    // face centres and raw fan normals
    faceCentres_.resize(static_cast<std::size_t>(nF));
    sf_.resize(static_cast<std::size_t>(nF));
    for (label f = 0; f < nF; ++f)
    {
        const auto& verts = faces_[static_cast<std::size_t>(f)];
        Vector c = Vector::zero();
        for (const label v : verts) c += points_[static_cast<std::size_t>(v)];
        faceCentres_[static_cast<std::size_t>(f)] = c / static_cast<scalar>(verts.size());

        Vector n = Vector::zero();
        const Vector& v0 = points_[static_cast<std::size_t>(verts[0])];
        for (std::size_t i = 1; i + 1 < verts.size(); ++i)
        {
            const Vector a = points_[verts[i]] - v0;
            const Vector b = points_[verts[i + 1]] - v0;
            n += cross(a, b);
        }
        sf_[static_cast<std::size_t>(f)] = n * 0.5;
    }

    // provisional cell centre estimate from face-centre averages
    Field<Vector> est(static_cast<std::size_t>(nCells_), Vector::zero());
    Field<label> cnt(static_cast<std::size_t>(nCells_), 0);
    for (label f = 0; f < nF; ++f)
    {
        const label o = owner_[static_cast<std::size_t>(f)];
        est[static_cast<std::size_t>(o)] += faceCentres_[static_cast<std::size_t>(f)];
        ++cnt[static_cast<std::size_t>(o)];
        if (f < nInternalFaces())
        {
            const label nb = neighbour_[static_cast<std::size_t>(f)];
            est[static_cast<std::size_t>(nb)] += faceCentres_[static_cast<std::size_t>(f)];
            ++cnt[static_cast<std::size_t>(nb)];
        }
    }
    for (label c = 0; c < nCells_; ++c)
    {
        if (cnt[static_cast<std::size_t>(c)] > 0)
        {
            est[static_cast<std::size_t>(c)] = est[static_cast<std::size_t>(c)]
                / static_cast<scalar>(cnt[static_cast<std::size_t>(c)]);
        }
    }

    // orient: internal faces owner -> neighbour; boundary faces outward
    for (label f = 0; f < nF; ++f)
    {
        const label o = owner_[static_cast<std::size_t>(f)];
        scalar side = 0;
        if (f < nInternalFaces())
        {
            const label nb = neighbour_[static_cast<std::size_t>(f)];
            side = dot(sf_[static_cast<std::size_t>(f)], est[static_cast<std::size_t>(nb)] - est[static_cast<std::size_t>(o)]);
        }
        else
        {
            side = dot(sf_[static_cast<std::size_t>(f)], faceCentres_[static_cast<std::size_t>(f)] - est[static_cast<std::size_t>(o)]);
        }
        if (side < 0)
        {
            sf_[static_cast<std::size_t>(f)] = -sf_[static_cast<std::size_t>(f)];
        }
    }

    magSf_.resize(static_cast<std::size_t>(nF));
    for (label f = 0; f < nF; ++f) magSf_[static_cast<std::size_t>(f)] = mag(sf_[static_cast<std::size_t>(f)]);

    // Exact cell volumes and first moments: the divergence theorem over the
    // fan triangles of every face tiles the cell with origin-based tetrahedra,
    // whose signed volumes sum to the cell volume and whose volume-weighted
    // centroids sum to its first moment.
    cellVolumes_.assign(static_cast<std::size_t>(nCells_), 0.0);
    cellCentres_.assign(static_cast<std::size_t>(nCells_), Vector::zero());
    Field<Vector> moment(static_cast<std::size_t>(nCells_), Vector::zero());
    for (label f = 0; f < nF; ++f)
    {
        const auto& verts = faces_[static_cast<std::size_t>(f)];
        const Vector& v0 = points_[static_cast<std::size_t>(verts[0])];
        scalar fv = 0;
        Vector fm = Vector::zero();
        for (std::size_t i = 1; i + 1 < verts.size(); ++i)
        {
            const Vector a = points_[verts[i]] - v0;
            const Vector b = points_[verts[i + 1]] - v0;
            const scalar tV = dot(v0, cross(a, b)) / 6.0;
            fv += tV;
            fm += tV * ((v0 + points_[verts[i]] + points_[verts[i + 1]]) / 4.0);
        }
        Vector raw = Vector::zero();
        for (std::size_t i = 1; i + 1 < verts.size(); ++i)
        {
            raw += cross(points_[verts[i]] - v0, points_[verts[i + 1]] - v0);
        }
        // fan orientation may agree or disagree with the oriented normal;
        // match sign against dot(rawNormal, Sf)
        const scalar s = dot(raw, sf_[static_cast<std::size_t>(f)]) >= 0 ? 1.0 : -1.0;
        const label o = owner_[static_cast<std::size_t>(f)];
        cellVolumes_[static_cast<std::size_t>(o)] += s * fv;
        moment[static_cast<std::size_t>(o)] += fm * s;
        if (f < nInternalFaces())
        {
            const label nb = neighbour_[static_cast<std::size_t>(f)];
            cellVolumes_[static_cast<std::size_t>(nb)] -= s * fv;
            moment[static_cast<std::size_t>(nb)] -= fm * s;
        }
    }
    for (label c = 0; c < nCells_; ++c)
    {
        // the signs of volume and moment cancel, so the quotient is the
        // centroid even for a cell whose decomposition came out negative
        cellCentres_[static_cast<std::size_t>(c)]
            = moment[static_cast<std::size_t>(c)] / cellVolumes_[static_cast<std::size_t>(c)];
        cellVolumes_[static_cast<std::size_t>(c)] = std::abs(cellVolumes_[static_cast<std::size_t>(c)]);
        if (cellVolumes_[static_cast<std::size_t>(c)] <= 0)
        {
            cfFatal("polyMesh: degenerate cell volume for cell " + std::to_string(c));
        }
    }

    distVec_.resize(static_cast<std::size_t>(nF));
    invDist_.resize(static_cast<std::size_t>(nF));
    for (label f = 0; f < nF; ++f)
    {
        const label o = owner_[static_cast<std::size_t>(f)];
        if (f < nInternalFaces())
        {
            const label nb = neighbour_[static_cast<std::size_t>(f)];
            distVec_[static_cast<std::size_t>(f)] = cellCentres_[static_cast<std::size_t>(nb)] - cellCentres_[static_cast<std::size_t>(o)];
        }
        else
        {
            distVec_[static_cast<std::size_t>(f)] = faceCentres_[static_cast<std::size_t>(f)] - cellCentres_[static_cast<std::size_t>(o)];
        }
        const scalar d = mag(distVec_[static_cast<std::size_t>(f)]);
        if (d <= 0)
        {
            cfFatal("polyMesh: zero owner-neighbour distance at face " + std::to_string(f));
        }
        invDist_[static_cast<std::size_t>(f)] = 1.0 / d;
    }
}

const Field<scalar>& polyMesh::cellVolumes() const { calcGeometry(); return cellVolumes_; }
const Field<Vector>& polyMesh::cellCentres() const { calcGeometry(); return cellCentres_; }
const Field<Vector>& polyMesh::faceCentres() const { calcGeometry(); return faceCentres_; }
const Field<Vector>& polyMesh::faceSf() const { calcGeometry(); return sf_; }
const Field<scalar>& polyMesh::faceAreas() const { calcGeometry(); return magSf_; }

Vector polyMesh::minX() const
{
    Vector m = points_[0];
    for (const Vector& p : points_)
    {
        m.x = std::min(m.x, p.x);
        m.y = std::min(m.y, p.y);
        m.z = std::min(m.z, p.z);
    }
    return m;
}

Vector polyMesh::maxX() const
{
    Vector m = points_[0];
    for (const Vector& p : points_)
    {
        m.x = std::max(m.x, p.x);
        m.y = std::max(m.y, p.y);
        m.z = std::max(m.z, p.z);
    }
    return m;
}

} // namespace cf
