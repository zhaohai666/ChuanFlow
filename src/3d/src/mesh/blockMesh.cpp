#include "cf/fvm/BlockMesh.h"

#include <cmath>
#include <map>
#include <iostream>
#include <tuple>

#include "cf/common/Error.h"

namespace cf
{
namespace
{

struct FaceRef
{
    std::vector<label> verts;
    label owner = -1;
    label nbr = -1; // -1 until second cell registers it
};

using FaceKey = std::tuple<label, label, label, label, label>;

// canonical (sorted) key for an up-to-5-vertex face: pad with -1
FaceKey makeKey(const std::vector<label>& v)
{
    std::vector<label> s = v;
    std::sort(s.begin(), s.end());
    s.resize(5, -1);
    return {s[0], s[1], s[2], s[3], s[4]};
}

} // namespace

polyMesh BlockMesh::create(const Vector& min, const Vector& max,
                           label nx, label ny, label nz,
                           const Sides& sides)
{
    if (nx <= 0 || ny <= 0 || nz <= 0)
    {
        cfFatal("blockMesh: cell counts must be positive");
    }

    auto pointAt = [&](label i, label j, label k) -> label
    {
        return static_cast<label>(static_cast<std::size_t>(i)
            + static_cast<std::size_t>(nx + 1) * (static_cast<std::size_t>(j)
            + static_cast<std::size_t>(ny + 1) * static_cast<std::size_t>(k)));
    };
    auto coordAt = [&](label i, label j, label k) -> Vector
    {
        return {min.x + (max.x - min.x) * static_cast<scalar>(i) / nx,
                min.y + (max.y - min.y) * static_cast<scalar>(j) / ny,
                min.z + (max.z - min.z) * static_cast<scalar>(k) / nz};
    };

    Field<Vector> points(static_cast<std::size_t>(nx + 1)
                        * static_cast<std::size_t>(ny + 1)
                        * static_cast<std::size_t>(nz + 1));
    for (label k = 0; k <= nz; ++k)
        for (label j = 0; j <= ny; ++j)
            for (label i = 0; i <= nx; ++i)
                points[static_cast<std::size_t>(pointAt(i, j, k))] = coordAt(i, j, k);

    auto cellId = [&](label i, label j, label k)
    {
        return static_cast<label>(static_cast<std::size_t>(i)
            + static_cast<std::size_t>(nx) * (static_cast<std::size_t>(j)
            + static_cast<std::size_t>(ny) * static_cast<std::size_t>(k)));
    };

    std::map<FaceKey, std::size_t> faceIndex;
    std::vector<FaceRef> faceRefs;

    auto addFace = [&](const std::vector<label>& verts, label owner)
    {
        const FaceKey key = makeKey(verts);
        auto it = faceIndex.find(key);
        if (it == faceIndex.end())
        {
            faceIndex.emplace(key, faceRefs.size());
            faceRefs.push_back({verts, owner, -1});
        }
        else
        {
            FaceRef& fr = faceRefs[it->second];
            if (fr.nbr != -1)
            {
                cfFatal("blockMesh: non-manifold face");
            }
            fr.nbr = owner;
        }
    };

    for (label k = 0; k < nz; ++k)
    {
        for (label j = 0; j < ny; ++j)
        {
            for (label i = 0; i < nx; ++i)
            {
                const label c = cellId(i, j, k);
                // corner numbering of the hex
                const label v0 = pointAt(i, j, k);
                const label v1 = pointAt(i + 1, j, k);
                const label v2 = pointAt(i + 1, j + 1, k);
                const label v3 = pointAt(i, j + 1, k);
                const label v4 = pointAt(i, j, k + 1);
                const label v5 = pointAt(i + 1, j, k + 1);
                const label v6 = pointAt(i + 1, j + 1, k + 1);
                const label v7 = pointAt(i, j + 1, k + 1);

                addFace({v0, v1, v2, v3}, c); // zlo
                addFace({v4, v5, v6, v7}, c); // zhi
                addFace({v0, v1, v5, v4}, c); // ylo
                addFace({v2, v3, v7, v6}, c); // yhi
                addFace({v0, v3, v7, v4}, c); // xlo
                addFace({v1, v2, v6, v5}, c); // xhi
            }
        }
    }

    // split into internal and boundary faces; classify boundary faces by side
    struct Boundary
    {
        std::vector<FaceRef> faces;
    };
    Boundary sidesArr[6]; // order: xmin xmax ymin ymax zmin zmax
    std::vector<FaceRef> internal;
    for (auto& fr : faceRefs)
    {
        if (fr.nbr != -1)
        {
            if (fr.owner > fr.nbr) std::swap(fr.owner, fr.nbr);
            internal.push_back(fr);
        }
        else
        {
            Vector cc = Vector::zero();
            for (const label v : fr.verts) cc += points[static_cast<std::size_t>(v)];
            cc = cc / static_cast<scalar>(fr.verts.size());
            const scalar eps = 1e-12;
            int side;
            if (std::abs(cc.x - min.x) < eps) side = 0;
            else if (std::abs(cc.x - max.x) < eps) side = 1;
            else if (std::abs(cc.y - min.y) < eps) side = 2;
            else if (std::abs(cc.y - max.y) < eps) side = 3;
            else if (std::abs(cc.z - min.z) < eps) side = 4;
            else if (std::abs(cc.z - max.z) < eps) side = 5;
            else cfFatal("blockMesh: boundary face not on any side");
            sidesArr[side].faces.push_back(fr);
        }
    }

    const std::pair<std::string, std::string>* sideNames[6] = {
        &sides.xmin, &sides.xmax, &sides.ymin, &sides.ymax, &sides.zmin, &sides.zmax};

    std::vector<std::vector<label>> faces;
    Field<label> owner, neighbour;
    std::vector<polyMesh::Patch> patches;

    faces.reserve(faceRefs.size());
    for (const auto& fr : internal)
    {
        faces.push_back(fr.verts);
        owner.push_back(fr.owner);
        neighbour.push_back(fr.nbr);
    }

    for (int s = 0; s < 6; ++s)
    {
        auto& list = sidesArr[s].faces;
        if (list.empty()) continue; // zero-thickness sides cannot happen; keep guard
        polyMesh::Patch p;
        p.name = sideNames[s]->first;
        p.type = sideNames[s]->second;
        p.start = static_cast<label>(faces.size());
        p.size = static_cast<label>(list.size());
        for (const auto& fr : list)
        {
            faces.push_back(fr.verts);
            owner.push_back(fr.owner);
        }
        patches.push_back(std::move(p));
    }

    return polyMesh(std::move(points), std::move(faces), std::move(owner),
                    std::move(neighbour), std::move(patches));
}

polyMesh BlockMesh::fromDict(const Dictionary& block)
{
    const auto vmin = block.getScalarList("min");
    const auto vmax = block.getScalarList("max");
    const auto ncell = block.getScalarList("cells");
    if (vmin.size() != 3 || vmax.size() != 3 || ncell.size() != 3)
    {
        cfFatal("blockMesh: min/max/cells need 3 values each");
    }
    const Vector mn{vmin[0], vmin[1], vmin[2]};
    const Vector mx{vmax[0], vmax[1], vmax[2]};

    Sides sides;
    struct SideEntry { const char* key; std::pair<std::string, std::string>* slot; };
    const SideEntry entries[6] = {
        {"xmin", &sides.xmin}, {"xmax", &sides.xmax},
        {"ymin", &sides.ymin}, {"ymax", &sides.ymax},
        {"zmin", &sides.zmin}, {"zmax", &sides.zmax}};
    for (const auto& e : entries)
    {
        if (block.has(e.key))
        {
            const Dictionary& sd = block.getDict(e.key);
            e.slot->first = sd.getWordOr("name", e.key);
            e.slot->second = sd.getWordOr("type", "wall");
        }
    }

    return create(mn, mx,
                  static_cast<label>(ncell[0]),
                  static_cast<label>(ncell[1]),
                  static_cast<label>(ncell[2]), sides);
}

} // namespace cf
