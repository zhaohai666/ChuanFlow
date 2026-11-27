// Adaptive Cartesian refinement operator (F13).
//
// The base block mesh is a structured hex grid. Cells whose directional
// gradient norm exceeds a declared threshold are refined by dyadically
// splitting the base bins along the axes the field varies on. Refinement is
// applied per base bin (a slab) so every internal face stays conforming: two
// neighbouring cells always share the exact same face vertices, and the
// rebuilt polyMesh therefore needs no hanging-node treatment.
//
// The stored cell value is the exact average of the analytic field over the
// cell. Splitting a parent into children partitions its volume, so the sum of
// child cell-integrals equals the parent cell-integral to round-off; summing
// over the refined mesh reproduces the analytic integral of the box
// independently of the refinement pattern. A final volume-weighted correction
// removes the residual round-off so conservation is exact to machine precision.
#include "cf/fvm/AdaptiveRefine.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <tuple>

#include "cf/common/Error.h"

namespace cf
{
namespace
{

scalar boxVolume(const Vector& lo, const Vector& hi)
{
    return (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z);
}

int readAxis(const Dictionary& adapt)
{
    const std::string a = adapt.getWordOr("axis", "x");
    if (a == "x") return 0;
    if (a == "y") return 1;
    if (a == "z") return 2;
    cfFatal("adaptiveRefine: field axis must be x, y or z (got '" + a + "')");
    return 0;
}

std::array<double, 3> triple3(const std::vector<double>& v, const std::string& what)
{
    if (v.size() != 3)
        cfFatal("adaptiveRefine: " + what + " needs exactly 3 components");
    return {v[0], v[1], v[2]};
}

struct SideSpec
{
    std::string name;
    std::string type;
};

struct BaseGrid
{
    Vector lo = Vector::zero();
    Vector hi = Vector::zero();
    std::array<label, 3> n{1, 1, 1};
    std::array<SideSpec, 6> sides;
};

BaseGrid readBaseGrid(const Dictionary& caseDict)
{
    if (!caseDict.has("block"))
        cfFatal("adaptiveRefine: case has no block{} to refine");

    const Dictionary& b = caseDict.getDict("block");
    BaseGrid g;
    const auto vmin = triple3(b.getScalarList("min"), "block min");
    const auto vmax = triple3(b.getScalarList("max"), "block max");
    const auto nc = triple3(b.getScalarList("cells"), "block cells");
    g.lo = Vector{vmin[0], vmin[1], vmin[2]};
    g.hi = Vector{vmax[0], vmax[1], vmax[2]};
    for (int a = 0; a < 3; ++a)
    {
        if (!(nc[a] >= 1) || nc[a] != static_cast<double>(static_cast<label>(nc[a])))
            cfFatal("adaptiveRefine: block cells must be positive integers");
        g.n[static_cast<std::size_t>(a)] = static_cast<label>(nc[a]);
        if (g.hi[a] <= g.lo[a])
            cfFatal("adaptiveRefine: block max must exceed min on every axis");
    }

    const char* keys[6] = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
    const char* defName[6] = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
    const char* defType[6] = {"wall", "wall", "wall", "wall", "symmetry", "symmetry"};
    for (int s = 0; s < 6; ++s)
    {
        g.sides[static_cast<std::size_t>(s)] = {defName[s], defType[s]};
        if (b.has(keys[s]))
        {
            const Dictionary& sd = b.getDict(keys[s]);
            g.sides[static_cast<std::size_t>(s)].name =
                sd.getWordOr("name", defName[s]);
            g.sides[static_cast<std::size_t>(s)].type =
                sd.getWordOr("type", defType[s]);
        }
    }
    return g;
}

// Dyadically refine one base bin `depth` times. The first and last points are
// the base-bin corners; intermediate points are repeated midpoints.
std::vector<double> splitBin(double lo, double hi, label depth)
{
    std::vector<double> e{lo, hi};
    for (label d = 0; d < depth; ++d)
    {
        std::vector<double> next;
        next.reserve(e.size() * 2 - 1);
        for (std::size_t i = 0; i + 1 < e.size(); ++i)
        {
            next.push_back(e[i]);
            next.push_back(0.5 * (e[i] + e[i + 1]));
        }
        next.push_back(e.back());
        e.swap(next);
    }
    return e;
}

std::vector<double> buildEdges(double lo, double hi, label ncell,
                               const std::vector<label>& depth)
{
    std::vector<double> edges;
    for (label i = 0; i < ncell; ++i)
    {
        const double a = lo + (hi - lo) * static_cast<double>(i)
                         / static_cast<double>(ncell);
        const double b = lo + (hi - lo) * static_cast<double>(i + 1)
                         / static_cast<double>(ncell);
        const std::vector<double> sub = splitBin(a, b, depth[static_cast<std::size_t>(i)]);
        const std::size_t first = (i == 0) ? 0 : 1; // shared bin corner: keep once
        for (std::size_t k = first; k < sub.size(); ++k)
            edges.push_back(sub[k]);
    }
    return edges;
}

label baseBinOf(double mid, double lo, double hi, label ncell)
{
    label b = static_cast<label>(std::floor((mid - lo) / (hi - lo)
                                            * static_cast<double>(ncell)));
    b = std::max<label>(0, std::min<label>(ncell - 1, b));
    return b;
}

/// Dyadic slab refinement depths: mark cells from the gradient norm and refine
/// the base bins of the axes the field actually varies on. Applying the split
/// to a whole bin keeps every shared face conforming (no hanging nodes in the
/// rebuilt mesh). Shared by refine() and the §5.229 remesh loop; the body is a
/// verbatim move of the original pass loop, so refine() output is unchanged.
std::array<std::vector<label>, 3> computeDepths(const BaseGrid& grid,
                                                const AdaptControls& ctl)
{
    std::array<std::vector<label>, 3> depth;
    for (int a = 0; a < 3; ++a)
        depth[static_cast<std::size_t>(a)].assign(
            static_cast<std::size_t>(grid.n[static_cast<std::size_t>(a)]), 0);

    const scalar gradEps = 1e-12 * std::max(1.0, ctl.threshold);

    for (label pass = 0; pass < ctl.maxDepth; ++pass)
    {
        const std::vector<double> ex = buildEdges(grid.lo.x, grid.hi.x, grid.n[0], depth[0]);
        const std::vector<double> ey = buildEdges(grid.lo.y, grid.hi.y, grid.n[1], depth[1]);
        const std::vector<double> ez = buildEdges(grid.lo.z, grid.hi.z, grid.n[2], depth[2]);

        std::array<std::vector<char>, 3> bump;
        bump[0].assign(static_cast<std::size_t>(grid.n[0]), 0);
        bump[1].assign(static_cast<std::size_t>(grid.n[1]), 0);
        bump[2].assign(static_cast<std::size_t>(grid.n[2]), 0);

        for (std::size_t ix = 0; ix + 1 < ex.size(); ++ix)
            for (std::size_t iy = 0; iy + 1 < ey.size(); ++iy)
                for (std::size_t iz = 0; iz + 1 < ez.size(); ++iz)
                {
                    const Vector lo{ex[ix], ey[iy], ez[iz]};
                    const Vector hi{ex[ix + 1], ey[iy + 1], ez[iz + 1]};
                    std::array<scalar, 3> g;
                    scalar g2 = 0;
                    for (int a = 0; a < 3; ++a)
                    {
                        g[static_cast<std::size_t>(a)] =
                            ctl.field.directionalGradient(lo, hi, a);
                        g2 += g[static_cast<std::size_t>(a)]
                              * g[static_cast<std::size_t>(a)];
                    }
                    if (std::sqrt(g2) <= ctl.threshold) continue;

                    std::array<label, 3> bin{
                        baseBinOf(0.5 * (lo.x + hi.x), grid.lo.x, grid.hi.x, grid.n[0]),
                        baseBinOf(0.5 * (lo.y + hi.y), grid.lo.y, grid.hi.y, grid.n[1]),
                        baseBinOf(0.5 * (lo.z + hi.z), grid.lo.z, grid.hi.z, grid.n[2])};
                    for (int a = 0; a < 3; ++a)
                    {
                        if (grid.n[static_cast<std::size_t>(a)] <= 1) continue;
                        if (g[static_cast<std::size_t>(a)] <= gradEps) continue;
                        const label b = bin[static_cast<std::size_t>(a)];
                        if (depth[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)]
                            < ctl.maxDepth)
                            bump[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = 1;
                    }
                }

        bool changed = false;
        for (int a = 0; a < 3; ++a)
            for (std::size_t b = 0; b < bump[static_cast<std::size_t>(a)].size(); ++b)
                if (bump[static_cast<std::size_t>(a)][b])
                {
                    ++depth[static_cast<std::size_t>(a)][b];
                    changed = true;
                }
        if (!changed) break;
    }
    return depth;
}

struct FaceRef
{
    std::vector<label> verts;
    label owner = -1;
    label nbr = -1;
};

using FaceKey = std::array<label, 4>;

FaceKey faceKey(const std::vector<label>& v)
{
    std::array<label, 4> s{v[0], v[1], v[2], v[3]};
    std::sort(s.begin(), s.end());
    return s;
}

polyMesh buildMesh(const std::vector<AdaptCell>& cells, const BaseGrid& grid)
{
    std::map<std::array<double, 3>, label> pointIndex;
    Field<Vector> points;
    auto addPoint = [&](double x, double y, double z) -> label
    {
        const std::array<double, 3> key{x, y, z};
        auto it = pointIndex.find(key);
        if (it != pointIndex.end()) return it->second;
        const label id = static_cast<label>(points.size());
        points.push_back(Vector{x, y, z});
        pointIndex.emplace(key, id);
        return id;
    };

    std::map<FaceKey, std::size_t> faceIndex;
    std::vector<FaceRef> faceRefs;
    auto addFace = [&](const std::vector<label>& verts, label owner)
    {
        const FaceKey key = faceKey(verts);
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
                cfFatal("adaptiveRefine: non-manifold face while rebuilding polyMesh");
            fr.nbr = owner;
        }
    };

    for (std::size_t ci = 0; ci < cells.size(); ++ci)
    {
        const AdaptCell& c = cells[ci];
        const label owner = static_cast<label>(ci);
        const label v0 = addPoint(c.lo.x, c.lo.y, c.lo.z);
        const label v1 = addPoint(c.hi.x, c.lo.y, c.lo.z);
        const label v2 = addPoint(c.hi.x, c.hi.y, c.lo.z);
        const label v3 = addPoint(c.lo.x, c.hi.y, c.lo.z);
        const label v4 = addPoint(c.lo.x, c.lo.y, c.hi.z);
        const label v5 = addPoint(c.hi.x, c.lo.y, c.hi.z);
        const label v6 = addPoint(c.hi.x, c.hi.y, c.hi.z);
        const label v7 = addPoint(c.lo.x, c.hi.y, c.hi.z);
        addFace({v0, v1, v2, v3}, owner); // zlo
        addFace({v4, v5, v6, v7}, owner); // zhi
        addFace({v0, v1, v5, v4}, owner); // ylo
        addFace({v2, v3, v7, v6}, owner); // yhi
        addFace({v0, v3, v7, v4}, owner); // xlo
        addFace({v1, v2, v6, v5}, owner); // xhi
    }

    // Split into internal faces and six boundary patches by side.
    std::vector<std::vector<FaceRef>> sides(6);
    std::vector<FaceRef> internal;
    const double eps = 1e-12 * std::max(1.0, std::max(
        std::max(std::abs(grid.hi.x - grid.lo.x), std::abs(grid.hi.y - grid.lo.y)),
        std::abs(grid.hi.z - grid.lo.z)));
    for (auto& fr : faceRefs)
    {
        if (fr.nbr != -1)
        {
            if (fr.owner > fr.nbr) std::swap(fr.owner, fr.nbr);
            internal.push_back(fr);
            continue;
        }
        Vector centre = Vector::zero();
        for (const label v : fr.verts) centre += points[static_cast<std::size_t>(v)];
        centre = centre / static_cast<scalar>(fr.verts.size());
        int side;
        if (std::abs(centre.x - grid.lo.x) < eps) side = 0;
        else if (std::abs(centre.x - grid.hi.x) < eps) side = 1;
        else if (std::abs(centre.y - grid.lo.y) < eps) side = 2;
        else if (std::abs(centre.y - grid.hi.y) < eps) side = 3;
        else if (std::abs(centre.z - grid.lo.z) < eps) side = 4;
        else if (std::abs(centre.z - grid.hi.z) < eps) side = 5;
        else cfFatal("adaptiveRefine: boundary face not on any side of the box");
        sides[static_cast<std::size_t>(side)].push_back(fr);
    }

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
        if (sides[static_cast<std::size_t>(s)].empty()) continue;
        polyMesh::Patch p;
        p.name = grid.sides[static_cast<std::size_t>(s)].name;
        p.type = grid.sides[static_cast<std::size_t>(s)].type;
        p.start = static_cast<label>(faces.size());
        p.size = static_cast<label>(sides[static_cast<std::size_t>(s)].size());
        for (const auto& fr : sides[static_cast<std::size_t>(s)])
        {
            faces.push_back(fr.verts);
            owner.push_back(fr.owner);
        }
        patches.push_back(std::move(p));
    }
    return polyMesh(std::move(points), std::move(faces), std::move(owner),
                    std::move(neighbour), std::move(patches));
}

} // namespace

scalar AdaptField::at(const Vector& x) const
{
    if (kind == Kind::Linear) return dot(gradient, x) + offset;
    return x[axis] < face ? low : high;
}

scalar AdaptField::boxIntegral(const Vector& lo, const Vector& hi) const
{
    const scalar vol = boxVolume(lo, hi);
    if (kind == Kind::Linear)
    {
        const Vector centre = (lo + hi) * 0.5;
        return (dot(gradient, centre) + offset) * vol;
    }
    const scalar t0 = lo[axis], t1 = hi[axis];
    const scalar clamped = std::min(std::max(face, t0), t1);
    const scalar lenLow = std::max(0.0, clamped - t0);
    scalar otherArea = 1.0;
    for (int a = 0; a < 3; ++a)
        if (a != axis) otherArea *= hi[a] - lo[a];
    const scalar volLow = lenLow * otherArea;
    return low * volLow + high * (vol - volLow);
}

scalar AdaptField::cellAverage(const Vector& lo, const Vector& hi) const
{
    const scalar vol = boxVolume(lo, hi);
    if (vol <= 0) cfFatal("adaptiveRefine: non-positive cell volume");
    return boxIntegral(lo, hi) / vol;
}

scalar AdaptField::directionalGradient(const Vector& lo, const Vector& hi, int a) const
{
    const scalar edge = hi[a] - lo[a];
    if (edge <= 0) return 0;
    Vector p = (lo + hi) * 0.5;
    p[a] = lo[a];
    const scalar vLo = at(p);
    p[a] = hi[a];
    const scalar vHi = at(p);
    return std::abs(vHi - vLo) / edge;
}

AdaptField AdaptField::fromDict(const Dictionary& adapt, const Vector& boxMin,
                                const Vector& boxMax)
{
    if (!adapt.has("field"))
        cfFatal("adaptiveRefine: adapt{} has no field linear|step");
    const std::string profile = adapt.getWord("field");
    AdaptField f;
    if (profile == "linear")
    {
        f.kind = Kind::Linear;
        const auto g = triple3(adapt.getScalarList("gradient"), "field gradient");
        f.gradient = Vector{g[0], g[1], g[2]};
        f.offset = adapt.getScalarOr("offset", 0.0);
    }
    else if (profile == "step")
    {
        f.kind = Kind::Step;
        f.axis = readAxis(adapt);
        f.face = adapt.getScalar("interface");
        f.low = adapt.getScalarOr("low", 0.0);
        f.high = adapt.getScalarOr("high", 1.0);
        if (f.high < f.low)
            cfFatal("adaptiveRefine: step high must not be below low");
        if (f.face <= boxMin[f.axis] || f.face >= boxMax[f.axis])
            cfFatal("adaptiveRefine: step interface must lie strictly inside the box");
    }
    else
    {
        cfFatal("adaptiveRefine: field must be linear or step (got '" + profile + "')");
    }
    return f;
}

AdaptControls AdaptControls::fromDict(const Dictionary& caseDict,
                                      const Vector& boxMin, const Vector& boxMax)
{
    if (!caseDict.has("adapt"))
        cfFatal("adaptiveRefine: case has no adapt{} block");
    const Dictionary& adapt = caseDict.getDict("adapt");
    AdaptControls c;
    if (!adapt.has("gradientThreshold"))
        cfFatal("adaptiveRefine: adapt{} has no gradientThreshold");
    c.threshold = adapt.getScalar("gradientThreshold");
    if (!(c.threshold > 0) || !std::isfinite(c.threshold))
        cfFatal("adaptiveRefine: gradientThreshold must be a finite positive number");

    if (!adapt.has("maxDepth"))
        cfFatal("adaptiveRefine: adapt{} has no maxDepth");
    const double rawDepth = adapt.getScalar("maxDepth");
    if (!(rawDepth >= 0) || rawDepth > 20 || rawDepth != std::floor(rawDepth))
        cfFatal("adaptiveRefine: maxDepth must be an integer in [0, 20]");
    c.maxDepth = static_cast<label>(rawDepth);

    c.field = AdaptField::fromDict(adapt, boxMin, boxMax);
    c.dumpCells = adapt.getWordOr("dumpCells", "no") == "yes";
    return c;
}

AdaptResult AdaptiveRefine::refine(const Dictionary& caseDict,
                                   const std::string& casePath)
{
    (void)casePath;
    const BaseGrid grid = readBaseGrid(caseDict);
    const AdaptControls ctl = AdaptControls::fromDict(caseDict, grid.lo, grid.hi);

    const std::array<std::vector<label>, 3> depth = computeDepths(grid, ctl);

    AdaptResult result;
    for (int a = 0; a < 3; ++a)
        for (const label d : depth[static_cast<std::size_t>(a)])
            result.level = std::max(result.level, d);

    const std::vector<double> ex = buildEdges(grid.lo.x, grid.hi.x, grid.n[0], depth[0]);
    const std::vector<double> ey = buildEdges(grid.lo.y, grid.hi.y, grid.n[1], depth[1]);
    const std::vector<double> ez = buildEdges(grid.lo.z, grid.hi.z, grid.n[2], depth[2]);
    for (std::size_t ix = 0; ix + 1 < ex.size(); ++ix)
        for (std::size_t iy = 0; iy + 1 < ey.size(); ++iy)
            for (std::size_t iz = 0; iz + 1 < ez.size(); ++iz)
            {
                const Vector lo{ex[ix], ey[iy], ez[iz]};
                const Vector hi{ex[ix + 1], ey[iy + 1], ez[iz + 1]};
                AdaptCell c;
                c.lo = lo;
                c.hi = hi;
                c.volume = boxVolume(lo, hi);
                c.phi = ctl.field.cellAverage(lo, hi); // parent -> child interpolation
                const label bx = baseBinOf(0.5 * (lo.x + hi.x), grid.lo.x, grid.hi.x, grid.n[0]);
                const label by = baseBinOf(0.5 * (lo.y + hi.y), grid.lo.y, grid.hi.y, grid.n[1]);
                const label bz = baseBinOf(0.5 * (lo.z + hi.z), grid.lo.z, grid.hi.z, grid.n[2]);
                c.level = std::max(depth[0][static_cast<std::size_t>(bx)],
                                   std::max(depth[1][static_cast<std::size_t>(by)],
                                            depth[2][static_cast<std::size_t>(bz)]));
                result.cells.push_back(c);
            }

    result.nCells = static_cast<label>(result.cells.size());
    result.nSplit = result.cells.empty() ? 0
        : static_cast<label>(std::count_if(result.cells.begin(), result.cells.end(),
              [](const AdaptCell& c) { return c.level > 0; }));

    const Vector boxLo = grid.lo, boxHi = grid.hi;
    const scalar boxVol = boxVolume(boxLo, boxHi);
    const scalar exact = ctl.field.boxIntegral(boxLo, boxHi);
    scalar computed = 0;
    for (const AdaptCell& c : result.cells) computed += c.phi * c.volume;

    // Exact discrete conservation: redistribute the round-off residue as a
    // uniform (volume-weighted) phi shift, so sum(phi_c V_c) == integral(phi)
    // to machine precision regardless of the refinement pattern.
    const scalar shift = (exact - computed) / boxVol;
    for (AdaptCell& c : result.cells) c.phi += shift;

    result.volumeSum = 0;
    result.integral = 0;
    result.minVolume = std::numeric_limits<scalar>::max();
    result.maxVolume = 0;
    for (const AdaptCell& c : result.cells)
    {
        result.volumeSum += c.volume;
        result.integral += c.phi * c.volume;
        result.minVolume = std::min(result.minVolume, c.volume);
        result.maxVolume = std::max(result.maxVolume, c.volume);
    }
    result.integralError = std::abs(result.integral - exact);
    if (!(result.volumeSum > 0) || !std::isfinite(result.integral))
        cfFatal("adaptiveRefine: refinement produced a degenerate mesh");

    result.mesh = buildMesh(result.cells, grid);

    // Cross-check the rebuilt topology against the analytic box.
    scalar meshVolume = 0;
    for (const scalar v : result.mesh->cellVolumes()) meshVolume += v;
    if (std::abs(meshVolume - boxVol) > 1e-9 * std::max(1.0, boxVol))
        cfFatal("adaptiveRefine: rebuilt polyMesh volume disagrees with the box");

    return result;
}

AmrConfig AmrConfig::fromDict(const Dictionary& caseDict)
{
    AmrConfig a;
    if (!caseDict.has("amr")) return a;
    const Dictionary& d = caseDict.getDict("amr");
    a.enabled = true;

    if (!d.has("refineThreshold"))
        cfFatal("amr: refineThreshold is required");
    a.refineThreshold = d.getScalar("refineThreshold");
    if (!std::isfinite(a.refineThreshold) || !(a.refineThreshold > 0))
        cfFatal("amr: refineThreshold must be a finite positive number");

    if (!d.has("coarsenThreshold"))
        cfFatal("amr: coarsenThreshold is required");
    a.coarsenThreshold = d.getScalar("coarsenThreshold");
    if (!std::isfinite(a.coarsenThreshold) || a.coarsenThreshold < 0)
        cfFatal("amr: coarsenThreshold must be finite and non-negative");
    if (a.coarsenThreshold >= a.refineThreshold)
        cfFatal("amr: coarsenThreshold must be below refineThreshold");

    if (!d.has("cycles"))
        cfFatal("amr: cycles is required");
    const double rawCycles = d.getScalar("cycles");
    if (!std::isfinite(rawCycles) || !(rawCycles >= 1) || rawCycles != std::floor(rawCycles))
        cfFatal("amr: cycles must be an integer of at least 1");
    a.cycles = static_cast<label>(rawCycles);

    const double rawEvery = d.getScalarOr("every", 1.0);
    if (!std::isfinite(rawEvery) || !(rawEvery >= 1) || rawEvery != std::floor(rawEvery))
        cfFatal("amr: every must be an integer of at least 1");
    a.every = static_cast<label>(rawEvery);

    a.diffusion = d.getScalarOr("diffusion", 0.0);
    if (!std::isfinite(a.diffusion) || a.diffusion < 0)
        cfFatal("amr: diffusion must be finite and non-negative");
    a.dt = d.getScalarOr("dt", 0.0);
    if (a.diffusion > 0 && !(a.dt > 0))
        cfFatal("amr: dt must be positive when diffusion is used");

    a.dumpCells = d.getWordOr("dumpCells", "no") == "yes";
    return a;
}

AmrLoopResult AdaptiveRefine::remesh(const Dictionary& caseDict, const AdaptResult& seed)
{
    const BaseGrid grid = readBaseGrid(caseDict);
    const AdaptControls ctl = AdaptControls::fromDict(caseDict, grid.lo, grid.hi);
    const AmrConfig amr = AmrConfig::fromDict(caseDict);

    // Seed state: the structured depth pattern of adapt{} plus the stored
    // cell field, in the ix -> iy -> iz order refine() produces.
    std::array<std::vector<label>, 3> depth = computeDepths(grid, ctl);
    std::array<std::vector<double>, 3> e{
        buildEdges(grid.lo.x, grid.hi.x, grid.n[0], depth[0]),
        buildEdges(grid.lo.y, grid.hi.y, grid.n[1], depth[1]),
        buildEdges(grid.lo.z, grid.hi.z, grid.n[2], depth[2])};
    std::array<label, 3> nc{static_cast<label>(e[0].size()) - 1,
                            static_cast<label>(e[1].size()) - 1,
                            static_cast<label>(e[2].size()) - 1};
    const std::size_t nCells = static_cast<std::size_t>(nc[0])
                             * static_cast<std::size_t>(nc[1])
                             * static_cast<std::size_t>(nc[2]);
    if (nCells != seed.cells.size())
        cfFatal("amr: seed cell count does not match the structured adapt{} mesh");

    std::vector<double> phi(nCells, 0.0);
    {
        std::size_t i = 0;
        for (label ix = 0; ix < nc[0]; ++ix)
            for (label iy = 0; iy < nc[1]; ++iy)
                for (label iz = 0; iz < nc[2]; ++iz)
                {
                    const Vector lo{e[0][ix], e[1][iy], e[2][iz]};
                    const Vector hi{e[0][ix + 1], e[1][iy + 1], e[2][iz + 1]};
                    if (std::abs(lo.x - seed.cells[i].lo.x) > 1e-12
                        || std::abs(hi.z - seed.cells[i].hi.z) > 1e-12)
                        cfFatal("amr: seed cell ordering does not match the structured mesh");
                    phi[i] = seed.cells[i].phi;
                    ++i;
                }
    }

    const auto cellVolume = [&](std::size_t ix, std::size_t iy, std::size_t iz) {
        return (e[0][ix + 1] - e[0][ix]) * (e[1][iy + 1] - e[1][iy])
             * (e[2][iz + 1] - e[2][iz]);
    };
    const auto idx = [&](std::size_t ix, std::size_t iy, std::size_t iz) {
        return (ix * static_cast<std::size_t>(nc[1]) + iy)
                 * static_cast<std::size_t>(nc[2]) + iz;
    };

    const scalar boxVol = boxVolume(grid.lo, grid.hi);
    scalar massTarget = 0;
    for (std::size_t ix = 0; ix < static_cast<std::size_t>(nc[0]); ++ix)
        for (std::size_t iy = 0; iy < static_cast<std::size_t>(nc[1]); ++iy)
            for (std::size_t iz = 0; iz < static_cast<std::size_t>(nc[2]); ++iz)
                massTarget += phi[idx(ix, iy, iz)] * cellVolume(ix, iy, iz);

    // Explicit conservative diffusion: flux = D A (phi_n - phi_c) / d with
    // d = half the sum of the two cell edges along the axis; the face area A
    // is shared by both cells because refinement is slab-based. Boundary
    // faces have no neighbour and contribute zero flux.
    const auto diffuseOnce = [&]() {
        std::vector<double> next = phi;
        for (std::size_t ix = 0; ix < static_cast<std::size_t>(nc[0]); ++ix)
            for (std::size_t iy = 0; iy < static_cast<std::size_t>(nc[1]); ++iy)
                for (std::size_t iz = 0; iz < static_cast<std::size_t>(nc[2]); ++iz)
                {
                    const std::size_t c = idx(ix, iy, iz);
                    scalar flux = 0;
                    for (int a = 0; a < 3; ++a)
                    {
                        const std::size_t ia = static_cast<std::size_t>(a);
                        const std::size_t i0 = (a == 0) ? ix : (a == 1) ? iy : iz;
                        const double hc = e[ia][i0 + 1] - e[ia][i0];
                        scalar area = 1;
                        for (int b = 0; b < 3; ++b)
                        {
                            if (b == a) continue;
                            const std::size_t ib = (b == 0) ? ix : (b == 1) ? iy : iz;
                            area *= e[static_cast<std::size_t>(b)][ib + 1]
                                  - e[static_cast<std::size_t>(b)][ib];
                        }
                        for (int dir = -1; dir <= 1; dir += 2)
                        {
                            const std::int64_t j0 = static_cast<std::int64_t>(i0) + dir;
                            if (j0 < 0 || j0 >= nc[ia]) continue;
                            const std::size_t uj0 = static_cast<std::size_t>(j0);
                            const double hn = e[ia][uj0 + 1] - e[ia][uj0];
                            const std::size_t nidx = (a == 0) ? idx(uj0, iy, iz)
                                             : (a == 1) ? idx(ix, uj0, iz)
                                                        : idx(ix, iy, uj0);
                            const double d = 0.5 * (hc + hn);
                            flux += amr.diffusion * area * (phi[nidx] - phi[c]) / d;
                        }
                    }
                    next[c] += amr.dt * flux / cellVolume(ix, iy, iz);
                }
        phi.swap(next);
    };

    AmrLoopResult loop;
    loop.initialCells = static_cast<label>(seed.cells.size());
    loop.massInitial = massTarget;

    for (label cycle = 1; cycle <= amr.cycles; ++cycle)
    {
        if (cycle > 1 && amr.diffusion > 0)
            for (label s = 0; s < amr.every; ++s) diffuseOnce();

        // Per-axis cell indicator: g_a = max over the two neighbours along a
        // of |phi_n - phi_c| / d (zero where no neighbour exists).
        const std::size_t curCells = static_cast<std::size_t>(nc[0])
                                   * static_cast<std::size_t>(nc[1])
                                   * static_cast<std::size_t>(nc[2]);
        std::vector<std::array<scalar, 3>> g(curCells);
        for (std::size_t ix = 0; ix < static_cast<std::size_t>(nc[0]); ++ix)
            for (std::size_t iy = 0; iy < static_cast<std::size_t>(nc[1]); ++iy)
                for (std::size_t iz = 0; iz < static_cast<std::size_t>(nc[2]); ++iz)
                {
                    const std::size_t c = idx(ix, iy, iz);
                    for (int a = 0; a < 3; ++a)
                    {
                        const std::size_t ia = static_cast<std::size_t>(a);
                        const std::size_t i0 = (a == 0) ? ix : (a == 1) ? iy : iz;
                        const double hc = e[ia][i0 + 1] - e[ia][i0];
                        scalar gmax = 0;
                        for (int dir = -1; dir <= 1; dir += 2)
                        {
                            const std::int64_t j0 = static_cast<std::int64_t>(i0) + dir;
                            if (j0 < 0 || j0 >= nc[ia]) continue;
                            const std::size_t uj0 = static_cast<std::size_t>(j0);
                            const double hn = e[ia][uj0 + 1] - e[ia][uj0];
                            const std::size_t nidx = (a == 0) ? idx(uj0, iy, iz)
                                             : (a == 1) ? idx(ix, uj0, iz)
                                                        : idx(ix, iy, uj0);
                            gmax = std::max(gmax, std::abs(phi[nidx] - phi[c])
                                                / (0.5 * (hc + hn)));
                        }
                        g[c][static_cast<std::size_t>(a)] = gmax;
                    }
                }

        // Slab decision per axis and base bin: refine when the strongest cell
        // of the bin exceeds the refine threshold, coarsen when every cell is
        // below the coarsen threshold.
        label nRefined = 0, nCoarsened = 0;
        for (int a = 0; a < 3; ++a)
        {
            const std::size_t ia = static_cast<std::size_t>(a);
            std::vector<scalar> binMax(static_cast<std::size_t>(grid.n[a]), 0);
            for (std::size_t ix = 0; ix < static_cast<std::size_t>(nc[0]); ++ix)
                for (std::size_t iy = 0; iy < static_cast<std::size_t>(nc[1]); ++iy)
                    for (std::size_t iz = 0; iz < static_cast<std::size_t>(nc[2]); ++iz)
                    {
                        const double centre = 0.5 * ((a == 0 ? e[0][ix] + e[0][ix + 1]
                                                  : a == 1 ? e[1][iy] + e[1][iy + 1]
                                                           : e[2][iz] + e[2][iz + 1]));
                        const label b = baseBinOf(centre, grid.lo[a], grid.hi[a], grid.n[a]);
                        binMax[static_cast<std::size_t>(b)] =
                            std::max(binMax[static_cast<std::size_t>(b)],
                                     g[idx(ix, iy, iz)][ia]);
                    }
            for (label b = 0; b < grid.n[a]; ++b)
            {
                const std::size_t ub = static_cast<std::size_t>(b);
                if (binMax[ub] > amr.refineThreshold
                    && depth[ia][ub] < ctl.maxDepth)
                {
                    ++depth[ia][ub];
                    ++nRefined;
                }
                else if (binMax[ub] < amr.coarsenThreshold && depth[ia][ub] > 0)
                {
                    --depth[ia][ub];
                    ++nCoarsened;
                }
            }
        }

        // Geometry rebuild + conservative field transfer: each new cell takes
        // the volume-weighted mean of the old cells it overlaps (refinement is
        // a single-parent copy, coarsening the weighted merge).
        const std::array<std::vector<double>, 3> oldE = e;
        const std::array<label, 3> oldNc = nc;
        std::vector<double> oldPhi = phi;
        e[0] = buildEdges(grid.lo.x, grid.hi.x, grid.n[0], depth[0]);
        e[1] = buildEdges(grid.lo.y, grid.hi.y, grid.n[1], depth[1]);
        e[2] = buildEdges(grid.lo.z, grid.hi.z, grid.n[2], depth[2]);
        nc = {static_cast<label>(e[0].size()) - 1, static_cast<label>(e[1].size()) - 1,
              static_cast<label>(e[2].size()) - 1};
        const std::size_t newCells = static_cast<std::size_t>(nc[0])
                                   * static_cast<std::size_t>(nc[1])
                                   * static_cast<std::size_t>(nc[2]);
        std::vector<double> newPhi(newCells, 0.0);
        for (std::size_t jx = 0; jx < static_cast<std::size_t>(nc[0]); ++jx)
            for (std::size_t jy = 0; jy < static_cast<std::size_t>(nc[1]); ++jy)
                for (std::size_t jz = 0; jz < static_cast<std::size_t>(nc[2]); ++jz)
                {
                    const double nlo[3] = {e[0][jx], e[1][jy], e[2][jz]};
                    const double nhi[3] = {e[0][jx + 1], e[1][jy + 1], e[2][jz + 1]};
                    const scalar vNew = (nhi[0] - nlo[0]) * (nhi[1] - nlo[1])
                                      * (nhi[2] - nlo[2]);
                    scalar acc = 0;
                    for (std::size_t ix = 0; ix < static_cast<std::size_t>(oldNc[0]); ++ix)
                        for (std::size_t iy = 0; iy < static_cast<std::size_t>(oldNc[1]); ++iy)
                            for (std::size_t iz = 0; iz < static_cast<std::size_t>(oldNc[2]); ++iz)
                            {
                                const double olo[3] = {oldE[0][ix], oldE[1][iy], oldE[2][iz]};
                                const double ohi[3] = {oldE[0][ix + 1], oldE[1][iy + 1],
                                                       oldE[2][iz + 1]};
                                scalar vInt = 1;
                                for (int a = 0; a < 3; ++a)
                                {
                                    const double ov = std::min(nhi[a], ohi[a])
                                                    - std::max(nlo[a], olo[a]);
                                    if (ov <= 0) { vInt = 0; break; }
                                    vInt *= ov;
                                }
                                if (vInt <= 0) continue;
                                acc += oldPhi[(ix * static_cast<std::size_t>(oldNc[1])
                                              + iy) * static_cast<std::size_t>(oldNc[2]) + iz]
                                     * vInt;
                            }
                    newPhi[(jx * static_cast<std::size_t>(nc[1]) + jy)
                           * static_cast<std::size_t>(nc[2]) + jz] = acc / vNew;
                }
        phi.swap(newPhi);

        // Pin the integral: a uniform shift removes the round-off residue so
        // sum(phi V) stays at the seed mass through every remeshing cycle.
        scalar mass = 0;
        for (std::size_t ix = 0; ix < static_cast<std::size_t>(nc[0]); ++ix)
            for (std::size_t iy = 0; iy < static_cast<std::size_t>(nc[1]); ++iy)
                for (std::size_t iz = 0; iz < static_cast<std::size_t>(nc[2]); ++iz)
                {
                    const scalar v = (e[0][ix + 1] - e[0][ix]) * (e[1][iy + 1] - e[1][iy])
                                   * (e[2][iz + 1] - e[2][iz]);
                    mass += phi[idx(ix, iy, iz)] * v;
                }
        const scalar shift = (massTarget - mass) / boxVol;
        for (double& p : phi) p += shift;

        AmrCycleReport rep;
        rep.cycle = cycle;
        rep.nRefined = nRefined;
        rep.nCoarsened = nCoarsened;
        rep.mass = 0;
        for (std::size_t ix = 0; ix < static_cast<std::size_t>(nc[0]); ++ix)
            for (std::size_t iy = 0; iy < static_cast<std::size_t>(nc[1]); ++iy)
                for (std::size_t iz = 0; iz < static_cast<std::size_t>(nc[2]); ++iz)
                {
                    AdaptCell c;
                    c.lo = Vector{e[0][ix], e[1][iy], e[2][iz]};
                    c.hi = Vector{e[0][ix + 1], e[1][iy + 1], e[2][iz + 1]};
                    c.volume = boxVolume(c.lo, c.hi);
                    c.phi = phi[idx(ix, iy, iz)];
                    const label bx = baseBinOf(0.5 * (c.lo.x + c.hi.x), grid.lo.x,
                                               grid.hi.x, grid.n[0]);
                    const label by = baseBinOf(0.5 * (c.lo.y + c.hi.y), grid.lo.y,
                                               grid.hi.y, grid.n[1]);
                    const label bz = baseBinOf(0.5 * (c.lo.z + c.hi.z), grid.lo.z,
                                               grid.hi.z, grid.n[2]);
                    c.level = std::max(depth[0][static_cast<std::size_t>(bx)],
                                       std::max(depth[1][static_cast<std::size_t>(by)],
                                                depth[2][static_cast<std::size_t>(bz)]));
                    rep.cells.push_back(c);
                    rep.mass += c.phi * c.volume;
                }
        rep.nCells = static_cast<label>(rep.cells.size());
        loop.cycles.push_back(std::move(rep));
    }

    if (loop.cycles.empty())
        cfFatal("amr: the remesh loop produced no cycles");
    loop.finalCells = loop.cycles.back().nCells;
    loop.massFinal = loop.cycles.back().mass;
    loop.massError = std::abs(loop.massFinal - loop.massInitial);
    return loop;
}

} // namespace cf
