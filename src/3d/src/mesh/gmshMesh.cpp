#include "cf/fvm/GmshMesh.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>

#include "cf/common/Error.h"
#include "cf/common/Log.h"

namespace cf
{
namespace
{

using FaceKey = std::tuple<label, label, label, label, label>;

/// Canonical key of a face: its sorted vertices, padded to five entries so a
/// triangle and a quad cannot collide.
FaceKey makeFaceKey(const std::vector<label>& verts)
{
    std::vector<label> s = verts;
    std::sort(s.begin(), s.end());
    s.resize(5, -1);
    return {s[0], s[1], s[2], s[3], s[4]};
}

int nodesPerElement(int type)
{
    switch (type)
    {
        case 1: return 2;   // line
        case 2: return 3;   // triangle
        case 3: return 4;   // quadrangle
        case 4: return 4;   // tetrahedron
        case 5: return 8;   // hexahedron
        case 6: return 6;   // prism (wedge)
        case 7: return 5;   // pyramid
        default: return 0;  // second order and beyond
    }
}

/// Local face rings of a volume element, shared by every cell of that type.
/// Gmsh numbers a prism as the bottom triangle 0-1-2 with the matching top
/// triangle 3-4-5 above it, and a hexahedron as bottom 0-3 plus top 4-7.
/// polyMesh re-orients every face from the owner/neighbour topology, so only
/// the cyclic order (the plane of the face) matters here, not its direction.
const std::vector<std::vector<int>>& elementFaceRings(int type)
{
    static const std::vector<std::vector<int>> tetra =
        {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    static const std::vector<std::vector<int>> hexa =
        {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
         {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
    static const std::vector<std::vector<int>> prism =
        {{0, 1, 2}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};
    static const std::vector<std::vector<int>> none;

    switch (type)
    {
        case 4: return tetra;
        case 5: return hexa;
        case 6: return prism;
        default: return none;
    }
}

std::vector<std::string> tokenize(const std::string& line)
{
    std::istringstream ss(line);
    std::vector<std::string> out;
    std::string w;
    while (ss >> w) out.push_back(w);
    return out;
}

long toInt(const std::string& tok, const std::string& where)
{
    try
    {
        return std::stol(tok);
    }
    catch (const std::exception&)
    {
        cfFatal("gmsh: expected an integer but found '" + tok + "' in " + where);
    }
}

double toDouble(const std::string& tok, const std::string& where)
{
    try
    {
        return std::stod(tok);
    }
    catch (const std::exception&)
    {
        cfFatal("gmsh: expected a number but found '" + tok + "' in " + where);
    }
}

bool isIntegerToken(const std::string& tok)
{
    return tok.find_first_of(".eE") == std::string::npos;
}

/// The non-empty token rows of one section, so the readers below can work
/// line-wise (Gmsh writes one entity per line) with a one-row look-ahead.
class Section
{
public:
    Section() = default;
    Section(const std::vector<std::string>& lines, std::size_t first, std::size_t last,
            const std::string& where)
        : where_(where)
    {
        for (std::size_t i = first; i < last; ++i)
        {
            auto t = tokenize(lines[i]);
            if (!t.empty()) rows_.push_back(std::move(t));
        }
    }

    bool found() const { return !rows_.empty(); }
    bool atEnd() const { return pos_ >= rows_.size(); }

    // The readers below take a const Section&; the cursor is the only state a
    // read mutates.
    std::vector<std::string> take() const
    {
        if (atEnd())
        {
            cfFatal("gmsh: truncated file: " + where_ + " ends after "
                    + std::to_string(rows_.size()) + " rows");
        }
        return rows_[pos_++];
    }

    /// Consumes the next row only if it is a single integer (a node tag line).
    bool takeSingleInt(long& value) const
    {
        if (atEnd() || rows_[pos_].size() != 1 || !isIntegerToken(rows_[pos_][0])) return false;
        value = std::stol(rows_[pos_][0]);
        ++pos_;
        return true;
    }

private:
    std::vector<std::vector<std::string>> rows_;
    mutable std::size_t pos_ = 0;
    std::string where_ = "section";
};

struct Element
{
    int type = 0;
    long physical = 0;         // physical group tag, 0 when ungrouped
    std::vector<label> nodes;  // dense point indices
};

struct MeshData
{
    std::vector<Vector> points;
    std::vector<Element> volumes;
    std::vector<Element> surfaces;
    long nTets = 0, nHexes = 0, nPrisms = 0;
};

Section findSection(const std::vector<std::string>& lines, const std::string& name,
                    const std::string& path, bool required)
{
    const std::string begin = "$" + name, end = "$End" + name;
    for (std::size_t i = 0; i < lines.size(); ++i)
    {
        const auto t = tokenize(lines[i]);
        if (t.empty() || t[0] != begin) continue;
        for (std::size_t j = i + 1; j < lines.size(); ++j)
        {
            const auto e = tokenize(lines[j]);
            if (!e.empty() && e[0] == end)
            {
                return Section(lines, i + 1, j, begin + " in " + path);
            }
        }
        cfFatal("gmsh: section " + begin + " in " + path + " is not terminated");
    }
    if (required) cfFatal("gmsh: " + path + " has no " + begin + " section");
    return Section();
}

void readPhysicalNames(const Section& sec, std::map<long, std::string>& surfaceNames)
{
    if (!sec.found()) return;
    const long n = toInt(sec.take()[0], "$PhysicalNames");
    for (long i = 0; i < n; ++i)
    {
        const auto t = sec.take();
        if (t.size() < 3) cfFatal("gmsh: malformed $PhysicalNames entry");
        const long dim = toInt(t[0], "$PhysicalNames");
        const long tag = toInt(t[1], "$PhysicalNames");
        std::string name;
        for (std::size_t k = 2; k < t.size(); ++k) name += t[k];
        name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
        if (name.empty()) cfFatal("gmsh: physical group " + std::to_string(tag) + " has no name");
        if (dim == 2) surfaceNames[tag] = name;
    }
}

/// v4.1 mapping entity -> physical group, written either as
/// `dim entityTag physicalTag "name"` ($EntityNames) or inside $Entities.
void readEntityNames(const Section& sec, std::map<long, long>& entityPhysical)
{
    if (!sec.found()) return;
    const long n = toInt(sec.take()[0], "$EntityNames");
    for (long i = 0; i < n; ++i)
    {
        const auto t = sec.take();
        if (t.size() < 4) cfFatal("gmsh: malformed $EntityNames entry");
        if (toInt(t[0], "$EntityNames") == 2)
        {
            entityPhysical[toInt(t[1], "$EntityNames")] = toInt(t[2], "$EntityNames");
        }
    }
}

/// $Entities: entities grouped by dimension, starting with the counts of
/// points, curves, surfaces and volumes. A surface line is
/// `tag minX minY minZ maxX maxY maxZ nPhysical physical...`; curves and
/// volumes share that prefix, points are just a position (Gmsh writes them
/// with 4 or 7 trailing numbers depending on the release), so only the
/// surfaces are read and only when the line is long enough. An entity with no
/// physical group falls back to bounding-box plane patch names.
void readEntities(const Section& sec, std::map<long, long>& entityPhysical)
{
    if (!sec.found()) return;
    const auto header = sec.take();
    if (header.size() < 4) cfFatal("gmsh: malformed $Entities header");
    for (int dim = 0; dim < 4; ++dim)
    {
        const long n = toInt(header[dim], "$Entities");
        for (long i = 0; i < n; ++i)
        {
            const auto t = sec.take();
            if (dim != 2) continue;
            if (t.size() < 8) cfFatal("gmsh: malformed $Entities surface line");
            const long tag = toInt(t[0], "$Entities");
            const long nPhys = toInt(t[7], "$Entities");
            if (nPhys > 0)
            {
                if (t.size() < 9) cfFatal("gmsh: malformed $Entities surface line");
                entityPhysical[tag] = toInt(t[8], "$Entities");
            }
        }
    }
}

void readNodes(const Section& sec, double version, const std::string& path,
               std::map<long, label>& tagIndex, std::vector<Vector>& points)
{
    if (!sec.found()) cfFatal("gmsh: " + path + " has an empty $Nodes section");

    std::vector<std::pair<long, Vector>> records;
    long declared = 0;

    auto coords = [](const std::vector<std::string>& t, const std::string& where)
    {
        if (t.size() < 3) cfFatal("gmsh: malformed node coordinate in " + where);
        return Vector{toDouble(t[0], where), toDouble(t[1], where), toDouble(t[2], where)};
    };

    if (version < 3.0)
    {
        declared = toInt(sec.take()[0], "$Nodes");
        for (long i = 0; i < declared; ++i)
        {
            const auto t = sec.take();
            if (t.size() < 4) cfFatal("gmsh: malformed node line in " + path);
            records.emplace_back(toInt(t[0], "$Nodes"), coords({t[1], t[2], t[3]}, "$Nodes"));
        }
    }
    else
    {
        const auto header = sec.take();
        if (header.size() < 4) cfFatal("gmsh: malformed $Nodes header in " + path);
        const long blocks = toInt(header[0], "$Nodes header");
        declared = toInt(header[1], "$Nodes header");
        const long firstTag = toInt(header[2], "$Nodes header");

        for (long b = 0; b < blocks; ++b)
        {
            const auto bh = sec.take();
            if (bh.size() < 4) cfFatal("gmsh: malformed $Nodes block header in " + path);
            const long numNodes = toInt(bh[3], "$Nodes block");
            if (numNodes <= 0) continue;

            // Two layouts are in the wild: the documented one, where the block
            // header carries the tag range and tags follow the coordinates
            // only when the numbering is not contiguous, and the one Gmsh
            // 4.10+ writes, where every node tag is a line before the
            // coordinates.
            const bool documented = bh.size() >= 6;
            const long minTag = documented ? toInt(bh[4], "$Nodes block") : firstTag;
            const long maxTag = documented ? toInt(bh[5], "$Nodes block") : firstTag;

            std::vector<long> tags;
            std::vector<Vector> xyz;

            if (documented)
            {
                for (long i = 0; i < numNodes; ++i) xyz.push_back(coords(sec.take(), "$Nodes"));
                if (minTag != maxTag)
                {
                    while (static_cast<long>(tags.size()) < numNodes)
                    {
                        long v = 0;
                        if (!sec.takeSingleInt(v)) break;
                        tags.push_back(v);
                    }
                }
                if (tags.empty())
                {
                    for (long i = 0; i < numNodes; ++i) tags.push_back(minTag + i);
                }
            }
            else
            {
                long v = 0;
                while (static_cast<long>(tags.size()) < numNodes && sec.takeSingleInt(v))
                {
                    tags.push_back(v);
                }
                const bool tagged = static_cast<long>(tags.size()) == numNodes;
                for (long i = 0; i < numNodes; ++i) xyz.push_back(coords(sec.take(), "$Nodes"));
                if (!tagged) tags.assign(numNodes, 0);
            }

            for (std::size_t i = 0; i < xyz.size(); ++i) records.emplace_back(tags[i], xyz[i]);
        }
    }

    if (static_cast<long>(records.size()) != declared)
    {
        cfFatal("gmsh: $Nodes declares " + std::to_string(declared) + " nodes, "
                + std::to_string(records.size()) + " were read from " + path
                + " (unsupported node layout?)");
    }

    // Nodes written without tags keep the order they appear in: give each the
    // lowest tag not already used by a tagged node.
    std::set<long> usedTags;
    for (const auto& r : records)
    {
        if (r.first > 0) usedTags.insert(r.first);
    }
    if (usedTags.size() < records.size())
    {
        long candidate = 1;
        for (auto& r : records)
        {
            if (r.first > 0) continue;
            while (usedTags.count(candidate)) ++candidate;
            r.first = candidate;
            usedTags.insert(candidate);
        }
    }

    std::sort(records.begin(), records.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    points.clear();
    tagIndex.clear();
    for (const auto& r : records)
    {
        if (tagIndex.count(r.first))
        {
            cfFatal("gmsh: duplicate node tag " + std::to_string(r.first) + " in " + path);
        }
        tagIndex[r.first] = static_cast<label>(points.size());
        points.push_back(r.second);
    }
}

void readElements(const Section& sec, double version, const std::string& path,
                  const std::map<long, label>& tagIndex,
                  const std::map<long, long>& entityPhysical, MeshData& mesh)
{
    if (!sec.found()) cfFatal("gmsh: " + path + " has an empty $Elements section");

    struct Row
    {
        int type = 0;
        long physical = 0;
        std::vector<std::string> nodes;
    };
    std::vector<Row> rows;

    if (version < 3.0)
    {
        const long declared = toInt(sec.take()[0], "$Elements");
        for (long i = 0; i < declared; ++i)
        {
            const auto t = sec.take();
            if (t.size() < 3) cfFatal("gmsh: malformed element line in " + path);
            const int type = static_cast<int>(toInt(t[1], "$Elements"));
            const int npn = nodesPerElement(type);
            if (npn <= 0)
            {
                cfFatal("gmsh: element type " + std::to_string(type) + " in " + path
                        + " is not a first-order element; mesh with Linear order");
            }
            const long nTags = toInt(t[2], "$Elements");
            if (static_cast<long>(t.size()) < 3 + nTags + npn)
            {
                cfFatal("gmsh: truncated element line in " + path);
            }
            Row r;
            r.type = type;
            r.physical = nTags > 0 ? toInt(t[3], "$Elements") : 0;
            r.nodes.assign(t.begin() + 3 + nTags, t.begin() + 3 + nTags + npn);
            rows.push_back(std::move(r));
        }
    }
    else
    {
        const auto header = sec.take();
        if (header.size() < 4) cfFatal("gmsh: malformed $Elements header in " + path);
        const long blocks = toInt(header[0], "$Elements header");
        const long declared = toInt(header[1], "$Elements header");
        long read = 0;

        for (long b = 0; b < blocks; ++b)
        {
            const auto bh = sec.take();
            if (bh.size() < 4) cfFatal("gmsh: malformed $Elements block header in " + path);
            const long entityTag = toInt(bh[1], "$Elements block");
            const int type = static_cast<int>(toInt(bh[2], "$Elements block"));
            const long numEl = toInt(bh[3], "$Elements block");
            const int npn = nodesPerElement(type);
            if (npn <= 0)
            {
                cfFatal("gmsh: element type " + std::to_string(type) + " in " + path
                        + " is not a first-order element; mesh with Linear order");
            }
            const auto pit = entityPhysical.find(entityTag);
            const long physical = pit != entityPhysical.end() ? pit->second : 0;

            bool tagged = false;
            for (long i = 0; i < numEl; ++i, ++read)
            {
                const auto t = sec.take();
                const long width = static_cast<long>(t.size());
                if (i == 0)
                {
                    // the element tag column is optional in 4.1: the width of
                    // the block's first line says whether it is written
                    if (width == npn + 1) tagged = true;
                    else if (width != npn)
                    {
                        cfFatal("gmsh: element line with " + std::to_string(width)
                                + " entries, expected " + std::to_string(npn) + " or "
                                + std::to_string(npn + 1));
                    }
                }
                else if (width != npn + (tagged ? 1 : 0))
                {
                    cfFatal("gmsh: inconsistent element line width in " + path);
                }
                Row r;
                r.type = type;
                r.physical = physical;
                r.nodes.assign(t.begin() + (tagged ? 1 : 0), t.end());
                rows.push_back(std::move(r));
            }
        }
        if (read != declared)
        {
            cfFatal("gmsh: $Elements declares " + std::to_string(declared) + " elements, "
                    + std::to_string(read) + " were read from " + path
                    + " (unsupported element layout?)");
        }
    }

    for (auto& r : rows)
    {
        Element el;
        el.type = r.type;
        el.physical = r.physical;
        el.nodes.reserve(r.nodes.size());
        for (const auto& tok : r.nodes)
        {
            const long tag = toInt(tok, "$Elements");
            const auto nit = tagIndex.find(tag);
            if (nit == tagIndex.end())
            {
                cfFatal("gmsh: an element refers to node " + std::to_string(tag)
                        + ", which $Nodes does not define");
            }
            el.nodes.push_back(nit->second);
        }
        if (std::unique(el.nodes.begin(), el.nodes.end()) != el.nodes.end())
        {
            cfFatal("gmsh: element " + std::to_string(el.type) + " repeats a node");
        }

        switch (el.type)
        {
            case 2:
            case 3:
                mesh.surfaces.push_back(std::move(el));
                break;
            case 4:
                ++mesh.nTets;
                mesh.volumes.push_back(std::move(el));
                break;
            case 5:
                ++mesh.nHexes;
                mesh.volumes.push_back(std::move(el));
                break;
            case 6:
                ++mesh.nPrisms;
                mesh.volumes.push_back(std::move(el));
                break;
            case 7:
                cfFatal("gmsh: pyramid elements are not supported");
            default:
                break;  // type 1: a boundary line says nothing about the volume
        }
    }
}

} // namespace

polyMesh GmshMesh::fromFile(const std::string& path,
                            const std::map<std::string, std::string>& patchTypes)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) cfFatal("gmsh: cannot open " + path);

    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line); )
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }

    const Section fmtSec = findSection(lines, "MeshFormat", path, true);
    const auto fmt = fmtSec.take();
    if (fmt.size() < 2) cfFatal("gmsh: malformed $MeshFormat in " + path);
    const double version = toDouble(fmt[0], "$MeshFormat");
    if (toInt(fmt[1], "$MeshFormat") != 0)
    {
        cfFatal("gmsh: " + path + " is a binary .msh; write it with Mesh.Binary = 0");
    }
    if (!(std::abs(version - 2.2) < 1e-9 || version >= 4.0))
    {
        cfFatal("gmsh: format '" + fmt[0] + "' is not supported (ASCII 2.2 and 4.1 are)");
    }

    std::map<long, std::string> physicalNames;
    readPhysicalNames(findSection(lines, "PhysicalNames", path, false), physicalNames);

    MeshData data;
    std::map<long, label> tagIndex;
    readNodes(findSection(lines, "Nodes", path, true), version, path, tagIndex, data.points);

    std::map<long, long> entityPhysical;
    if (version >= 3.0)
    {
        readEntityNames(findSection(lines, "EntityNames", path, false), entityPhysical);
        readEntities(findSection(lines, "Entities", path, false), entityPhysical);
    }
    readElements(findSection(lines, "Elements", path, true), version, path,
                 tagIndex, entityPhysical, data);

    if (data.volumes.empty()) cfFatal("gmsh: " + path + " has no volume elements");

    // ---- faces from the volume connectivity ----
    std::map<FaceKey, std::size_t> faceIndex;
    struct FaceRef { std::vector<label> verts; label cell0 = -1; label cell1 = -1; };
    std::vector<FaceRef> faceRefs;
    for (label c = 0; c < static_cast<label>(data.volumes.size()); ++c)
    {
        const Element& el = data.volumes[static_cast<std::size_t>(c)];
        for (const auto& ring : elementFaceRings(el.type))
        {
            std::vector<label> verts;
            verts.reserve(ring.size());
            for (const int k : ring) verts.push_back(el.nodes[static_cast<std::size_t>(k)]);

            const auto res = faceIndex.emplace(makeFaceKey(verts), faceRefs.size());
            if (res.second)
            {
                faceRefs.push_back({verts, c, -1});
            }
            else
            {
                FaceRef& fr = faceRefs[res.first->second];
                if (fr.cell1 != -1)
                {
                    cfFatal("gmsh: a face belongs to three cells (cell "
                            + std::to_string(fr.cell0) + "/" + std::to_string(fr.cell1)
                            + "/" + std::to_string(c) + "): the mesh is not manifold");
                }
                fr.cell1 = c;
            }
        }
    }

    // ---- boundary face names from the surface elements ----
    std::map<FaceKey, std::string> surfaceNames;
    for (const auto& el : data.surfaces)
    {
        const auto nit = physicalNames.find(el.physical);
        if (nit == physicalNames.end()) continue;  // ungrouped surface
        surfaceNames.emplace(makeFaceKey(el.nodes), nit->second);
    }
    const bool hasSurfaceNames = !surfaceNames.empty();

    auto patchNameOf = [&](const std::vector<label>& verts) -> std::string
    {
        const auto it = surfaceNames.find(makeFaceKey(verts));
        if (it != surfaceNames.end()) return it->second;
        if (verts.size() == 4)
        {
            // a triangulated surface splits a quad face along a diagonal
            for (int skip = 0; skip < 4; ++skip)
            {
                std::vector<label> tri;
                for (int k = 0; k < 4; ++k) if (k != skip) tri.push_back(verts[static_cast<std::size_t>(k)]);
                const auto t = surfaceNames.find(makeFaceKey(tri));
                if (t != surfaceNames.end()) return t->second;
            }
        }
        return {};
    };

    // bounding-box planes, used when the file carries no surface groups
    Vector lo = data.points[0], hi = data.points[0];
    for (const auto& v : data.points)
    {
        lo.x = std::min(lo.x, v.x); hi.x = std::max(hi.x, v.x);
        lo.y = std::min(lo.y, v.y); hi.y = std::max(hi.y, v.y);
        lo.z = std::min(lo.z, v.z); hi.z = std::max(hi.z, v.z);
    }
    const scalar extent = std::max({hi.x - lo.x, hi.y - lo.y, hi.z - lo.z, 1e-300});
    const scalar eps = 1e-9 * extent;
    const char* planeNames[6] = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};

    auto planeName = [&](const std::vector<label>& verts) -> std::string
    {
        for (int d = 0; d < 3; ++d)
        {
            for (int s = 0; s < 2; ++s)
            {
                const scalar lim = s ? hi[d] : lo[d];
                bool on = true;
                for (const label v : verts)
                {
                    if (std::abs(data.points[static_cast<std::size_t>(v)][d] - lim) > eps)
                    {
                        on = false;
                        break;
                    }
                }
                if (on) return planeNames[2 * d + s];
            }
        }
        return {};
    };

    // ---- split into internal faces and patches ----
    std::vector<FaceRef> internal;
    std::vector<std::pair<std::string, std::vector<FaceRef>>> patches;
    auto addToPatch = [&](const std::string& name, const FaceRef& fr)
    {
        for (auto& p : patches)
        {
            if (p.first == name) { p.second.push_back(fr); return; }
        }
        patches.emplace_back(name, std::vector<FaceRef>{fr});
    };

    long planeFallbackCount = 0;
    for (auto& fr : faceRefs)
    {
        if (fr.cell1 != -1)
        {
            if (fr.cell0 > fr.cell1) std::swap(fr.cell0, fr.cell1);
            internal.push_back(fr);
            continue;
        }
        std::string name = hasSurfaceNames ? patchNameOf(fr.verts) : std::string();
        if (name.empty())
        {
            name = planeName(fr.verts);
            if (name.empty())
            {
                std::string verts;
                for (const label v : fr.verts) verts += " " + std::to_string(v);
                cfFatal("gmsh: boundary face with nodes" + verts + " belongs to no physical "
                        "group and lies on no bounding-box plane; give every boundary surface "
                        "a physical group in Gmsh");
            }
            ++planeFallbackCount;
        }
        addToPatch(name, fr);
    }

    if (planeFallbackCount > 0)
    {
        log::warn("gmsh: ", planeFallbackCount,
                  " boundary faces were named after bounding-box planes");
    }

    std::vector<std::vector<label>> faces;
    Field<label> owner, neighbour;
    std::vector<polyMesh::Patch> meshPatches;
    faces.reserve(faceRefs.size());
    owner.reserve(faceRefs.size());
    for (const auto& fr : internal)
    {
        faces.push_back(fr.verts);
        owner.push_back(fr.cell0);
        neighbour.push_back(fr.cell1);
    }
    for (const auto& p : patches)
    {
        polyMesh::Patch patch;
        patch.name = p.first;
        const auto t = patchTypes.find(p.first);
        patch.type = t != patchTypes.end() ? t->second : "wall";
        patch.start = static_cast<label>(faces.size());
        patch.size = static_cast<label>(p.second.size());
        for (const auto& fr : p.second)
        {
            faces.push_back(fr.verts);
            owner.push_back(fr.cell0);
        }
        meshPatches.push_back(std::move(patch));
    }

    for (const auto& want : patchTypes)
    {
        const bool present = std::any_of(meshPatches.begin(), meshPatches.end(),
                                        [&](const polyMesh::Patch& p) { return p.name == want.first; });
        if (!present)
        {
            log::warn("gmsh: patch ", want.first,
                      " is declared in the case but not present in the mesh");
        }
    }

    Field<Vector> points(data.points.begin(), data.points.end());

    polyMesh mesh(std::move(points), std::move(faces), std::move(owner),
                  std::move(neighbour), std::move(meshPatches));

    scalar volume = 0;
    for (const scalar v : mesh.cellVolumes()) volume += v;
    log::info("gmsh: msh ", (version < 3.0 ? "2.2" : "4.1"), " file ", path);
    log::info("gmsh: cells ", mesh.nCells(), " (tet ", data.nTets, ", hex ", data.nHexes,
              ", prism ", data.nPrisms, "), faces ", mesh.nFaces(), " (internal ",
              mesh.nInternalFaces(), ", boundary ", mesh.nBoundaryFaces(), "), nodes ",
              mesh.nPoints());
    log::info("gmsh: bounding box (", lo.x, " ", lo.y, " ", lo.z, ") -> (",
              hi.x, " ", hi.y, " ", hi.z, "), total volume ", volume);
    for (const auto& p : mesh.patches())
    {
        log::info("gmsh: patch ", p.name, " type ", p.type, " faces ", p.size);
    }

    return mesh;
}

} // namespace cf
