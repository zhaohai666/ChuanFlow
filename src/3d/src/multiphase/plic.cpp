#include "cf/multiphase/Plic.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

#include "cf/common/Error.h"
#include "cf/fvm/FvMatrix.h"

namespace cf::vof
{

namespace
{

constexpr scalar kInfinity = std::numeric_limits<scalar>::infinity();

/// Newell's area vector of a vertex loop: half the sum of the cross products of
/// each triangle of the fan. For a planar loop it is perpendicular to the plane
/// with the polygon area as its magnitude, and its sign tells which way the loop
/// winds.
Vector loopAreaVector(const std::vector<Vector>& loop)
{
    Vector a = Vector::zero();
    if (loop.size() < 3) return a;
    const Vector p0 = loop[0];
    for (std::size_t i = 1; i + 1 < loop.size(); ++i)
    {
        a += cross(loop[i] - p0, loop[i + 1] - p0);
    }
    return a * 0.5;
}

/// Polygon centroid by the same fan, area-weighted. Valid for a convex polygon
/// in any winding (the areas are taken as magnitudes).
Vector loopCentroid(const std::vector<Vector>& loop)
{
    if (loop.size() < 3) return Vector::zero();
    const Vector p0 = loop[0];
    Vector c = Vector::zero();
    scalar area = 0;
    for (std::size_t i = 1; i + 1 < loop.size(); ++i)
    {
        const Vector e1 = loop[i] - p0;
        const Vector e2 = loop[i + 1] - p0;
        const Vector n = cross(e1, e2);
        const scalar a = 0.5 * mag(n);
        if (a <= 0) continue;
        c += (p0 + loop[i] + loop[i + 1]) * (a * (1.0 / 3.0));
        area += a;
    }
    return area > 0 ? c * (1.0 / area) : Vector::zero();
}

/// Volume and first moment of the solid enclosed by outward-oriented convex
/// polygon loops, from the signed tetrahedra each triangle forms with a fixed
/// origin. Both sums are independent of that origin for a closed surface, so the
/// cell's own vertex mean is used to keep the products small; `moment` is
/// measured from that origin and becomes a centroid by adding it back.
struct Solid
{
    scalar volume = 0;
    Vector moment = Vector::zero();

    Vector centroid() const
    {
        return volume > 0 ? moment * (1.0 / volume) : Vector::zero();
    }
};

void addPolygon(const std::vector<Vector>& loop, const Vector& origin, Solid& s)
{
    if (loop.size() < 3) return;
    const Vector a = loop[0] - origin;
    for (std::size_t i = 1; i + 1 < loop.size(); ++i)
    {
        const Vector b = loop[i] - origin;
        const Vector c = loop[i + 1] - origin;
        // signed volume of the tetrahedron (origin, loop[0], loop[i], loop[i+1])
        const scalar v = dot(a, cross(b, c)) * (1.0 / 6.0);
        s.volume += v;
        // its centroid is (origin + 4 points)/4, i.e. (a + b + c)/4 shifted
        s.moment += v * (a + b + c) * 0.25;
    }
}

Solid solidOf(const std::vector<std::vector<Vector>>& loops, const Vector& origin)
{
    Solid s;
    for (const std::vector<Vector>& l : loops) addPolygon(l, origin, s);
    return s;
}

Vector vertexMean(const std::vector<std::vector<Vector>>& loops)
{
    Vector m = Vector::zero();
    label n = 0;
    for (const std::vector<Vector>& l : loops)
        for (const Vector& p : l) { m += p; ++n; }
    return n > 0 ? m * (1.0 / static_cast<scalar>(n)) : Vector::zero();
}

/// Clip one convex loop against the plane n.x = offset, keeping n.x >= offset.
/// The winding of the input loop is preserved.
std::vector<Vector> clipLoop(const std::vector<Vector>& loop, const Vector& unitNormal,
                             scalar offset, scalar eps)
{
    std::vector<Vector> out;
    const std::size_t n = loop.size();
    if (n < 3) return out;
    for (std::size_t i = 0; i < n; ++i)
    {
        const Vector& a = loop[i];
        const Vector& b = loop[(i + 1) % n];
        const scalar sa = dot(unitNormal, a) - offset;
        const scalar sb = dot(unitNormal, b) - offset;
        if (sa >= -eps) out.push_back(a);
        if ((sa > eps && sb < -eps) || (sa < -eps && sb > eps))
            out.push_back(a + (b - a) * (sa / (sa - sb)));
    }
    return out;
}

/// Order the points where the plane crossed the cell's edges into the convex
/// polygon they span, wound counter-clockwise as seen from outside the kept side
/// - which is the `-unitNormal` side, since the kept side is n.x >= offset.
std::vector<Vector> orderCapPolygon(std::vector<Vector> pts, const Vector& unitNormal,
                                    scalar tol)
{
    // the same crossing is seen by the two faces sharing the edge it lies on
    std::sort(pts.begin(), pts.end(), [](const Vector& a, const Vector& b) {
        if (a.x != b.x) return a.x < b.x;
        if (a.y != b.y) return a.y < b.y;
        return a.z < b.z;
    });
    std::vector<Vector> uniq;
    for (const Vector& p : pts)
    {
        if (!uniq.empty() && mag(p - uniq.back()) <= tol) continue;
        uniq.push_back(p);
    }
    if (uniq.size() < 3) return {};

    Vector c = Vector::zero();
    for (const Vector& p : uniq) c += p;
    c = c * (1.0 / static_cast<scalar>(uniq.size()));

    Vector e1 = unitNormal.y > unitNormal.x && unitNormal.y > unitNormal.z
        ? Vector{1, 0, 1}
        : Vector{1, 0, 0};
    e1 = e1 - unitNormal * dot(e1, unitNormal);
    if (mag2(e1) <= 0) e1 = Vector{0, 1, 0} - unitNormal * unitNormal.y;
    e1 = e1 * (1.0 / mag(e1));
    const Vector e2 = cross(unitNormal, e1);

    std::sort(uniq.begin(), uniq.end(), [&](const Vector& a, const Vector& b) {
        const Vector da = a - c, db = b - c;
        // decreasing angle about +n is counter-clockwise about -n
        return std::atan2(dot(da, e2), dot(da, e1))
             > std::atan2(dot(db, e2), dot(db, e1));
    });
    return uniq;
}

} // namespace

std::vector<CellFace> cellBoundary(const polyMesh& mesh, label cell)
{
    if (cell < 0 || cell >= mesh.nCells()) cfFatal("plic: no such cell");
    const auto& faces = mesh.cells()[static_cast<std::size_t>(cell)];
    const auto& loops = mesh.faces();
    const auto& points = mesh.points();
    const auto& sf = mesh.faceSf();
    const auto& owner = mesh.owner();

    std::vector<CellFace> out;
    out.reserve(faces.size());
    for (label f : faces)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        CellFace cf;
        cf.loop.reserve(loops[u].size());
        for (label v : loops[u]) cf.loop.push_back(points[static_cast<std::size_t>(v)]);
        cf.areaVector = sf[u];
        // faceSf points owner -> neighbour, so it is outward for the owner only
        if (owner[u] != cell) cf.areaVector = -cf.areaVector;

        if (dot(loopAreaVector(cf.loop), cf.areaVector) < 0)
            std::reverse(cf.loop.begin(), cf.loop.end());
        out.push_back(std::move(cf));
    }
    return out;
}

namespace
{

/// The volume and centroid of the solid a closed surface encloses, taken from the
/// surface itself. For a cell's own boundary this is also the check that the rewound
/// loops close: the same signed-tetrahedra sum has to reproduce the mesh's volume.
struct AbsoluteSolid
{
    scalar volume = 0;
    Vector centroid = Vector::zero();
};

AbsoluteSolid absoluteSolid(const std::vector<std::vector<Vector>>& loops)
{
    const Vector origin = vertexMean(loops);
    const Solid s = solidOf(loops, origin);
    return {s.volume, s.volume > 0 ? origin + s.centroid() : Vector::zero()};
}

std::vector<std::vector<Vector>> loopsOf(const std::vector<CellFace>& cell)
{
    std::vector<std::vector<Vector>> loops;
    loops.reserve(cell.size());
    for (const CellFace& cf : cell) loops.push_back(cf.loop);
    return loops;
}

/// Clip every loop of a closed surface by the plane `n . x = offset`, keeping the
/// `>=` side, and close the result with the cap the plane cuts. The tolerances follow
/// the surface's own extent along n, which is what lets a region that has already
/// been cut once be cut again on the same terms as a whole cell. `whole` is the plane
/// never reaching the solid, `void_` its taking all of it - which is also what a
/// plane grazing through a corner or along a face gives, since nothing of finite
/// volume is left.
struct ClipResult
{
    std::vector<std::vector<Vector>> loops;
    bool whole = false;
    bool void_ = false;
    /// index of the cap loop, -1 when the cut added none
    label cap = -1;
};

ClipResult clipLoops(const std::vector<std::vector<Vector>>& loops, const Vector& unitNormal,
                     scalar offset)
{
    if (mag2(unitNormal) <= 0) cfFatal("plic: the cut normal must be a unit vector");

    scalar smin = kInfinity, smax = -kInfinity;
    for (const std::vector<Vector>& l : loops)
        for (const Vector& p : l)
        {
            const scalar s = dot(unitNormal, p);
            smin = std::min(smin, s);
            smax = std::max(smax, s);
        }
    if (smin > smax) cfFatal("plic: the region has no boundary to cut");

    const scalar extent = smax - smin;
    const scalar eps = 1e-12 * (extent > 0 ? extent : 1.0);
    ClipResult res;
    if (offset <= smin - eps) { res.whole = true; res.loops = loops; return res; }
    if (offset >= smax + eps) { res.void_ = true; return res; }

    for (const std::vector<Vector>& l : loops)
    {
        std::vector<Vector> clipped = clipLoop(l, unitNormal, offset, eps);
        if (clipped.size() >= 3) res.loops.push_back(std::move(clipped));
    }

    // The cap's corners are the clipped loops' points that lie on the plane. That is
    // the edge crossings, but reading them off the result rather than out of the clip
    // is what also picks up the surface's own vertices when the plane runs exactly
    // through them - a cut through a corner, whose section would otherwise be missed
    // and the kept solid left open. It only applies when the plane is inside the
    // region: on its own lowest plane the clipped loops already are the closed
    // boundary, and a cap there would be that face twice.
    const scalar capTol = 1e-9 * (extent > 0 ? extent : 1.0);
    if (offset > smin + eps)
    {
        std::vector<Vector> crossings;
        for (const std::vector<Vector>& l : res.loops)
            for (const Vector& p : l)
                if (std::abs(dot(unitNormal, p) - offset) <= capTol) crossings.push_back(p);
        const std::vector<Vector> cap = orderCapPolygon(std::move(crossings), unitNormal, capTol);
        if (cap.size() >= 3)
        {
            res.cap = static_cast<label>(res.loops.size());
            res.loops.push_back(cap);
        }
    }

    if (absoluteSolid(res.loops).volume <= 0)
    {
        res.void_ = true;
        res.loops.clear();
        res.cap = -1;
    }
    return res;
}

} // namespace

scalar ConvexRegion::volume() const
{
    if (empty || loops.empty()) return 0;
    return absoluteSolid(loops).volume;
}

Vector ConvexRegion::centroid() const
{
    if (empty || loops.empty()) return Vector::zero();
    return absoluteSolid(loops).centroid;
}

ConvexRegion cellRegion(const polyMesh& mesh, label cell)
{
    ConvexRegion r;
    r.loops = loopsOf(cellBoundary(mesh, cell));
    if (r.loops.empty()) cfFatal("plic: the cell has no boundary");
    return r;
}

ConvexRegion clipRegion(const ConvexRegion& region, const Vector& unitNormal, scalar offset)
{
    if (region.empty) return ConvexRegion{{}, true};
    const ClipResult res = clipLoops(region.loops, unitNormal, offset);
    if (res.void_) return ConvexRegion{{}, true};
    ConvexRegion out;
    out.loops = std::move(res.loops);
    return out;
}

ConvexRegion clipRegion(const ConvexRegion& region, const std::vector<HalfSpace>& halfSpaces)
{
    ConvexRegion out = region;
    for (const HalfSpace& h : halfSpaces)
    {
        out = clipRegion(out, h.unitNormal, h.offset);
        if (out.empty) return out;
    }
    return out;
}

Cut clipConvexCell(const std::vector<CellFace>& cell, const Vector& unitNormal,
                   scalar offset)
{
    Cut cut;
    const ClipResult res = clipLoops(loopsOf(cell), unitNormal, offset);
    if (res.whole) { cut.whole = true; return cut; }
    if (res.void_) { cut.void_ = true; return cut; }

    const AbsoluteSolid s = absoluteSolid(res.loops);
    cut.volume = s.volume;
    cut.centroid = s.centroid;
    if (res.cap >= 0)
    {
        cut.area = mag(loopAreaVector(res.loops[static_cast<std::size_t>(res.cap)]));
        cut.areaCentroid = loopCentroid(res.loops[static_cast<std::size_t>(res.cap)]);
    }
    return cut;
}

scalar offsetForFraction(const std::vector<CellFace>& cell, const Vector& unitNormal,
                         scalar fraction)
{
    if (!(fraction >= 0 && fraction <= 1)) cfFatal("plic: fraction outside [0,1]");
    scalar smin = kInfinity, smax = -kInfinity;
    for (const CellFace& cf : cell)
    {
        for (const Vector& p : cf.loop)
        {
            const scalar s = dot(unitNormal, p);
            smin = std::min(smin, s);
            smax = std::max(smax, s);
        }
    }
    if (smin > smax) cfFatal("plic: the cell has no boundary");
    if (fraction <= 0) return smax;
    if (fraction >= 1) return smin;

    const AbsoluteSolid solid = absoluteSolid(loopsOf(cell));
    const scalar volume = solid.volume;
    if (!(volume > 0)) cfFatal("plic: a cell with no volume cannot be reconstructed");
    const scalar target = fraction * volume;

    // Newton on V(offset) = target: dV/d(offset) is exactly minus the area of
    // the cut, so the step is (V - target)/area. The bracket is tightened every
    // round, which is what keeps the sweep safe when the step leaves the cell.
    scalar lo = smin, hi = smax;
    scalar c = 0.5 * (lo + hi);
    for (int iter = 0; iter < 80; ++iter)
    {
        const Cut cut = clipConvexCell(cell, unitNormal, c);
        const scalar v = cut.whole ? volume : (cut.void_ ? 0 : cut.volume);
        if (v > target) lo = c; else hi = c;
        const scalar residual = v - target;
        if (std::abs(residual) <= 1e-15 * volume) break;
        if (cut.area > 0)
        {
            const scalar next = c + residual / cut.area;
            if (next > lo && next < hi) { c = next; continue; }
        }
        c = 0.5 * (lo + hi);
    }
    return c;
}

namespace
{

/// One neighbour's say about where the interface is: the plane that cuts its own
/// fraction exactly, told by its constant, the area of the cut and the centre of that
/// area.
struct PlaneSample
{
    scalar offset = 0;
    scalar area = 0;
    Vector centre = Vector::zero();
};

/// Two unit vectors spanning the plane perpendicular to `n`: the directions a tilt of
/// the normal lives in. Starting from the axis `n` is least aligned with keeps the
/// Gram-Schmidt step away from dividing by what it just removed.
void tangentBasis(const Vector& n, Vector& e1, Vector& e2)
{
    Vector a = Vector::zero();
    if (std::abs(n.x) <= std::abs(n.y) && std::abs(n.x) <= std::abs(n.z)) a = Vector{1, 0, 0};
    else if (std::abs(n.y) <= std::abs(n.z)) a = Vector{0, 1, 0};
    else a = Vector{0, 0, 1};
    e1 = a - n * dot(a, n);
    e1 = e1 * (1.0 / mag(e1));
    e2 = cross(n, e1);
}

/// The tilt of a plane normal that brings a neighbourhood of exact-volume cuts
/// together, and what they still disagree by once it is applied.
///
/// Cell j, cut to its own fraction by the plane with normal n, holds its interface at
/// the constant `o_j`. Turning the normal by `d` - perpendicular to n, so it is a
/// rotation and not a stretch - moves that constant to `o_j + d . p_j`, with `p_j` the
/// centre of the cut: a plane rotated about a point of itself leaves the volume on
/// either side unchanged to first order, because its area moments about that point
/// vanish. So one plane explains the neighbourhood when `o_j + d . p_j` is the same for
/// every cell, and the least-squares `d` is the slope of the offsets along the
/// interface. The cut area weights the statement, because a cell the plane grazes has
/// that much less to say about where the interface runs.
///
/// The two scatters are that disagreement without and with the tilt, on the same
/// weighting. On a plane front both are nothing; on a curved one what is left is the
/// sagitta over the neighbourhood, which no single plane can hold and no normal fit
/// can know.
struct PlaneTilt
{
    Vector tilt = Vector::zero();
    scalar rmsBefore = 0;
    scalar rms = 0;
};

PlaneTilt fitTilt(const std::vector<PlaneSample>& samples, const Vector& n)
{
    PlaneTilt out;
    if (samples.size() < 2) return out; // nothing else to agree with

    scalar weight = 0, offsetSum = 0;
    Vector centreSum = Vector::zero();
    for (const PlaneSample& s : samples)
    {
        weight += s.area;
        offsetSum += s.offset * s.area;
        centreSum += s.centre * s.area;
    }
    if (!(weight > 0)) return out;
    const scalar offsetMean = offsetSum / weight;
    const Vector centreMean = centreSum * (1.0 / weight);

    Vector e1, e2;
    tangentBasis(n, e1, e2);

    scalar g[2][2] = {{0, 0}, {0, 0}}, b[2] = {0, 0}, scatter = 0;
    for (const PlaneSample& s : samples)
    {
        const Vector q = s.centre - centreMean;
        const scalar t[2] = {dot(q, e1), dot(q, e2)};
        const scalar e = s.offset - offsetMean;
        scatter += s.area * e * e;
        for (int i = 0; i < 2; ++i)
        {
            b[i] += s.area * e * t[i];
            for (int j = 0; j < 2; ++j) g[i][j] += s.area * t[i] * t[j];
        }
    }
    out.rmsBefore = std::sqrt(scatter / weight);

    const scalar trace = g[0][0] + g[1][1];
    if (!(trace > 0)) return out; // every cut in one place: no slope to read
    // Levenberg-Marquardt damping, and it is doing the work here rather than
    // rounding off a singularity. A neighbourhood whose cut faces lie nearly along
    // one line - a front crossing the mesh diagonally does exactly that - has little
    // in one of the two tilt directions, so the slope read along it is the noise of
    // the other: the iteration walks that way and stops somewhere the neighbourhood
    // does not agree with at all. Damped, the disk's worst cell is off the radial by
    // 1.8 deg where undamped it is 23 at every mesh, and the disagreement the fit
    // reports there falls with the cell instead of staying a tenth of one. Nothing
    // is biased by it: the plane the iteration stops at is the one with no slope left
    // in b, which the damping cannot move, and 1e-3 gives the same numbers as 1e-2.
    g[0][0] += 1e-2 * trace;
    g[1][1] += 1e-2 * trace;
    const scalar det = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    if (!(det > 0)) return out;

    scalar d[2] = {-(g[1][1] * b[0] - g[0][1] * b[1]) / det,
                   -(-g[1][0] * b[0] + g[0][0] * b[1]) / det};
    // One statement about where each neighbour's cut sits is good for a rotation of it,
    // not for a different interface: past this the slope has stopped being the tilt.
    const scalar magnitude = std::sqrt(d[0] * d[0] + d[1] * d[1]);
    if (magnitude > 0.35)
    {
        d[0] *= 0.35 / magnitude;
        d[1] *= 0.35 / magnitude;
    }
    out.tilt = e1 * d[0] + e2 * d[1];

    scalar left = 0;
    for (const PlaneSample& s : samples)
    {
        const scalar r = s.offset - offsetMean + dot(out.tilt, s.centre - centreMean);
        left += s.area * r * r;
    }
    out.rms = std::sqrt(left / weight);
    return out;
}

/// How many neighbourhoods of tilt one fit is allowed to take, and a round that
/// changes nothing stops early. Four is where the disk says to stop. A front that
/// really is a plane is settled by then - 0 deg off its analytic normal, its plane
/// constant off by 5e-8 of a cell, and the neighbourhood agreeing to 6e-9 of one,
/// where two rounds still leave 0.004 deg and 6e-5 of a cell - and on a curved front
/// more rounds cost rather than buy: at sixteen the mean angle of the R/dx = 8 disk
/// has gone from 1.2 to 2.9 deg and its reported area from 97 to 91 per cent of the
/// perimeter. The extra rounds are spent explaining the staircase of the fractions,
/// which is the thing the fit is here not to look at.
constexpr int kFitRounds = 4;

/// For every cell flagged in `member`, the member cells that share any of its
/// vertices. Along the interface the band cells of a cell-filled front touch at
/// corners rather than faces - the front is one cell thick, so the cells ahead
/// and behind along it meet it only at their corners - so face adjacency would
/// leave a cell with nothing but its own cut to agree with. A cell is its own
/// neighbour: it sits at its own vertices.
std::vector<std::vector<label>> cornerNeighbourhoods(const polyMesh& mesh,
                                                     const std::vector<char>& member)
{
    const label nC = mesh.nCells();
    std::vector<std::vector<label>> vertexCells(static_cast<std::size_t>(mesh.nPoints()));
    std::vector<label> vs;
    for (label c = 0; c < nC; ++c)
    {
        if (!member[static_cast<std::size_t>(c)]) continue;
        vs.clear();
        for (label f : mesh.cells()[static_cast<std::size_t>(c)])
            for (label v : mesh.faces()[static_cast<std::size_t>(f)]) vs.push_back(v);
        std::sort(vs.begin(), vs.end());
        vs.erase(std::unique(vs.begin(), vs.end()), vs.end());
        for (label v : vs) vertexCells[static_cast<std::size_t>(v)].push_back(c);
    }

    std::vector<std::vector<label>> out(static_cast<std::size_t>(nC));
    std::vector<char> mark(static_cast<std::size_t>(nC), 0);
    std::vector<label> touched;
    for (label c = 0; c < nC; ++c)
    {
        if (!member[static_cast<std::size_t>(c)]) continue;
        vs.clear();
        touched.clear();
        for (label f : mesh.cells()[static_cast<std::size_t>(c)])
            for (label v : mesh.faces()[static_cast<std::size_t>(f)]) vs.push_back(v);
        std::sort(vs.begin(), vs.end());
        vs.erase(std::unique(vs.begin(), vs.end()), vs.end());
        for (label v : vs)
            for (label nb : vertexCells[static_cast<std::size_t>(v)])
            {
                if (mark[static_cast<std::size_t>(nb)]) continue;
                mark[static_cast<std::size_t>(nb)] = 1;
                touched.push_back(nb);
                out[static_cast<std::size_t>(c)].push_back(nb);
            }
        for (label nb : touched) mark[static_cast<std::size_t>(nb)] = 0;
    }
    return out;
}

/// Cholesky of a symmetric `k`-by-`k` (k <= 6) system expected positive
/// definite, destroying `a`, which is filled row by row up to k. False when a
/// pivot comes out non-positive: the model is not determined by the data and the
/// caller is to report no answer rather than a made-up one.
bool solveSpd6(scalar a[6][6], const scalar b[6], scalar x[6], int k)
{
    for (int i = 0; i < k; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            scalar s = a[i][j];
            for (int q = 0; q < j; ++q) s -= a[i][q] * a[j][q];
            if (i == j)
            {
                if (!(s > 0)) return false;
                a[i][j] = std::sqrt(s);
            }
            else
            {
                a[i][j] = s / a[j][j];
            }
        }
    }
    for (int i = 0; i < k; ++i)
    {
        scalar s = b[i];
        for (int q = 0; q < i; ++q) s -= a[i][q] * x[q];
        x[i] = s / a[i][i];
    }
    for (int i = k - 1; i >= 0; --i)
    {
        scalar s = x[i];
        for (int q = i + 1; q < k; ++q) s -= a[q][i] * x[q];
        x[i] = s / a[i][i];
    }
    return true;
}

} // namespace

Reconstruction reconstruct(const polyMesh& mesh, const Field<scalar>& alpha,
                           const Field<Vector>& youngsGradient, NormalChoice choice,
                           scalar bandTolerance)
{
    const label nC = mesh.nCells();
    if (static_cast<label>(alpha.size()) != nC)
        cfFatal("plic: alpha does not match the mesh");
    if (static_cast<label>(youngsGradient.size()) != nC)
        cfFatal("plic: the starting normal does not match the mesh");

    const std::size_t n = static_cast<std::size_t>(nC);
    Reconstruction r{Field<Vector>(n, Vector::zero()),
                     Field<scalar>(n, 0.0),
                     Field<scalar>(n, 0.0),
                     Field<Vector>(n, Vector::zero()),
                     Field<scalar>(n, 0.0),
                     Field<scalar>(n, 0.0),
                     0,
                     bandTolerance};

    // ---- the band, and the normal every cell starts from ----
    // What the fit reads is the neighbours' fractions rather than their first planes,
    // so this round keeps only what the reported plane needs: the cell, its volume and
    // the direction the field's own gradient gives it.
    struct Seed
    {
        label cell = -1;
        std::vector<CellFace> boundary;
        scalar volume = 0;
        scalar fraction = 0;
        Vector normal = Vector::zero();
        Vector phase1MinusPhase0 = Vector::zero();
        bool centroidWorks = false;
    };
    std::vector<Seed> seeds;
    Field<label> seedOf(n, -1);

    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const scalar a = alpha[u];
        if (!(a > bandTolerance && a < 1.0 - bandTolerance)) continue;

        const Vector g = youngsGradient[u];
        if (mag2(g) <= 0) continue; // no direction to start the cut from

        Seed s;
        s.cell = c;
        s.boundary = cellBoundary(mesh, c);
        s.fraction = a;
        s.normal = g * (1.0 / mag(g));
        s.volume = absoluteSolid(loopsOf(s.boundary)).volume;
        if (!(s.volume > 0)) cfFatal("plic: cell with no volume in the interface band");

        if (choice == NormalChoice::Centroid || choice == NormalChoice::Mixed)
        {
            // the two halves of the first cut, for the direction from one centroid to
            // the other: both have to be there for that difference to mean anything
            const scalar o0 = offsetForFraction(s.boundary, s.normal, a);
            const Cut phase1 = clipConvexCell(s.boundary, s.normal, o0);
            const Vector neg{-s.normal.x, -s.normal.y, -s.normal.z};
            const Cut phase0 = clipConvexCell(s.boundary, neg, -o0);
            if (!phase1.whole && !phase1.void_ && !phase0.whole && !phase0.void_)
            {
                s.centroidWorks = true;
                s.phase1MinusPhase0 = phase1.centroid - phase0.centroid;
            }
        }

        seedOf[u] = static_cast<label>(seeds.size());
        seeds.push_back(std::move(s));
    }

    // ---- the band cells that meet at each vertex ----
    // The neighbourhood the offset fit reads, from the shared corner-adjacency
    // builder below.
    std::vector<std::vector<label>> neighbourhood;
    if (choice == NormalChoice::Fit)
    {
        std::vector<char> member(n, 0);
        for (const Seed& s : seeds) member[static_cast<std::size_t>(s.cell)] = 1;
        const std::vector<std::vector<label>> byCell = cornerNeighbourhoods(mesh, member);
        neighbourhood.resize(seeds.size());
        for (std::size_t k = 0; k < seeds.size(); ++k)
            for (label nb : byCell[static_cast<std::size_t>(seeds[k].cell)])
            {
                const label s = seedOf[static_cast<std::size_t>(nb)];
                if (s >= 0) neighbourhood[k].push_back(s);
            }
    }

    std::vector<PlaneSample> samples;

    for (std::size_t si = 0; si < seeds.size(); ++si)
    {
        const Seed& s = seeds[si];
        const std::size_t u = static_cast<std::size_t>(s.cell);
        Vector n = s.normal;
        scalar scatter = 0, scatterBefore = 0;
        label rounds = 0;

        // The cell's own sub-cell geometry gives the direction for the estimators
        // that read a cell on its own; for the fit it only starts the search, and the
        // direction itself comes from the neighbourhood.
        if (choice == NormalChoice::Fit)
        {
            // Tilt a candidate normal until the neighbourhood's own cuts lie in one
            // plane, and report what still disagrees.
            const auto settle = [&](const Vector& start)
            {
                Vector m = start;
                scalar first = 0, rest = 0;
                label took = 0;
                for (int round = 0; round < kFitRounds; ++round)
                {
                    // every neighbour of the band, cut to its own fraction by the normal
                    // as it stands: the data the fit is allowed to see
                    samples.clear();
                    for (label k : neighbourhood[si])
                    {
                        const Seed& t = seeds[static_cast<std::size_t>(k)];
                        PlaneSample ps;
                        ps.offset = offsetForFraction(t.boundary, m, t.fraction);
                        const Cut c = clipConvexCell(t.boundary, m, ps.offset);
                        ps.area = c.area;
                        ps.centre = c.areaCentroid;
                        samples.push_back(ps);
                    }
                    const PlaneTilt tilt = fitTilt(samples, m);
                    if (round == 0) first = tilt.rmsBefore;
                    rest = tilt.rms;
                    ++took;
                    if (mag2(tilt.tilt) <= 1e-16) break; // the normal the neighbourhood agrees on
                    const Vector nt = m + tilt.tilt;
                    m = nt * (1.0 / mag(nt));
                }
                return std::make_tuple(m, first, rest, took);
            };

            Vector m;
            std::tie(m, scatterBefore, scatter, rounds) = settle(n);
            n = m;
        }
        else if (s.centroidWorks)
        {
            // from the centroid of what is not phase 1 to the centroid of what
            // is: the direction the phase actually increases along, measured on
            // the sub-cell geometry rather than on the field's slope
            const Vector cd = s.phase1MinusPhase0;
            if (mag2(cd) > 0)
            {
                const Vector u1 = cd * (1.0 / mag(cd));
                if (choice == NormalChoice::Centroid)
                {
                    n = u1;
                }
                else
                {
                    const Vector m = u1 + n;
                    if (mag2(m) > 0) n = m * (1.0 / mag(m));
                }
            }
        }

        // the plane the final normal cuts is the reconstruction that is reported
        const scalar off = offsetForFraction(s.boundary, n, s.fraction);
        const Cut kept = clipConvexCell(s.boundary, n, off);
        if (kept.whole || kept.void_ || !(kept.volume > 0))
            cfFatal("plic: the reconstruction lost the interface it was given");

        r.normal[u] = n;
        r.offset[u] = off;
        r.area[u] = kept.area;
        r.areaCentroid[u] = kept.areaCentroid;
        r.fraction[u] = kept.volume / s.volume;
        r.fitScatter[u] = scatter;
        r.fitScatterBefore += scatterBefore;
        r.fitRounds = std::max(r.fitRounds, rounds);
        ++r.bandCells;
    }
    if (r.bandCells > 0) r.fitScatterBefore /= static_cast<scalar>(r.bandCells);
    return r;
}

Reconstruction reconstruct(const volScalarField& alpha, const std::string& gradScheme,
                           NormalChoice choice, scalar bandTolerance)
{
    const Field<Vector> g = fvc::grad(alpha, gradScheme);
    return reconstruct(alpha.mesh(), alpha.internalField(), g, choice, bandTolerance);
}

Field<scalar> interfaceCurvaturePlic(const polyMesh& mesh, const Reconstruction& r)
{
    const label nC = mesh.nCells();
    const std::size_t n = static_cast<std::size_t>(nC);
    Field<scalar> kappa(n, 0.0);

    std::vector<char> member(n, 0);
    for (label c = 0; c < nC; ++c)
        member[static_cast<std::size_t>(c)] = mag2(r.normal[static_cast<std::size_t>(c)]) > 0;

    // two rings of corner adjacency: the cells sharing a vertex with any cell that
    // shares a vertex with this one. One ring is the cells the interface passes
    // through a single shared corner with, whose points sit within a fraction of a
    // cell of each other; the arc a curvature is worth reading needs the couple of
    // cells on either side along the front, which is where the second ring takes it.
    const std::vector<std::vector<label>> ring1 = cornerNeighbourhoods(mesh, member);
    std::vector<std::vector<label>> byCell(n);
    std::vector<char> seen(n, 0);
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        if (!member[u]) continue;
        std::vector<label>& two = byCell[u];
        for (label v1 : ring1[u])
        {
            if (!seen[static_cast<std::size_t>(v1)])
            {
                seen[static_cast<std::size_t>(v1)] = 1;
                two.push_back(v1);
            }
            for (label v2 : ring1[static_cast<std::size_t>(v1)])
                if (!seen[static_cast<std::size_t>(v2)])
                {
                    seen[static_cast<std::size_t>(v2)] = 1;
                    two.push_back(v2);
                }
        }
        for (label v : two) seen[static_cast<std::size_t>(v)] = 0;
    }

    // the powers of (x, y) each column of the quadric basis carries: x^2, x y,
    // y^2, x, y, 1
    constexpr int px[6] = {2, 1, 0, 1, 0, 0};
    constexpr int py[6] = {0, 1, 2, 0, 1, 0};

    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        if (!member[u]) continue;
        const Vector zAxis = r.normal[u];
        const Vector origin = r.areaCentroid[u];
        Vector e1, e2;
        tangentBasis(zAxis, e1, e2);

        // the neighbourhood's interface points in the cell's own frame, the cell's
        // own point at the origin among them. A point says as much as the interface
        // area its cell holds, so a plane the reconstruction could only graze
        // cannot outvote one that crosses the cell.
        std::vector<scalar> sx, sy, sz, sw;
        scalar span = 0;
        const auto add = [&](const Vector& p, scalar w)
        {
            if (!(w > 0)) return;
            const Vector d = p - origin;
            const scalar t1 = dot(d, e1), t2 = dot(d, e2);
            sx.push_back(t1);
            sy.push_back(t2);
            sz.push_back(dot(d, zAxis));
            sw.push_back(w);
            span = std::max(span, std::sqrt(t1 * t1 + t2 * t2));
        };
        add(origin, r.area[u]);
        for (label nb : byCell[u])
        {
            if (static_cast<std::size_t>(nb) == u) continue;
            const std::size_t v = static_cast<std::size_t>(nb);
            add(r.areaCentroid[v], r.area[v]);
        }
        if (sx.size() < 3 || !(span > 0)) continue;

        // coordinates scaled to the neighbourhood's own size, so every basis
        // column is O(1) and no direction is quiet just because the arc it reads
        // is a few cells of a neighbourhood a column tall
        const scalar invL = 1.0 / span;
        scalar g[6][6] = {}, rhs[6] = {}, col[6] = {};
        for (std::size_t i = 0; i < sx.size(); ++i)
        {
            const scalar x = sx[i] * invL, y = sy[i] * invL;
            const scalar z = sz[i] * invL, w = sw[i];
            scalar basis[6];
            for (int k = 0; k < 6; ++k)
            {
                scalar bk = 1;
                if (px[k] >= 1) bk *= x;
                if (px[k] == 2) bk *= x;
                if (py[k] >= 1) bk *= y;
                if (py[k] == 2) bk *= y;
                basis[k] = bk;
            }
            for (int j = 0; j < 6; ++j)
            {
                rhs[j] += w * z * basis[j];
                col[j] += w * basis[j] * basis[j];
                for (int k = 0; k < 6; ++k) g[j][k] += w * basis[j] * basis[k];
            }
        }

        // A direction the data does not reach - the extrusion direction of a
        // single-layer mesh leaves every point at one height, so x, x^2 and x y
        // carry nothing - is dropped rather than fitted: its coefficient is 0, the
        // honest answer "no curvature readable this way", where fitting it would
        // mean dividing a round-off column norm into the ridge and calling what
        // comes back a curvature. The surviving columns are equilibrated to unit
        // norm before the ridge, so the damping is the same statement about every
        // direction that is actually there.
        const scalar colMax = *std::max_element(col, col + 6);
        if (!(colMax > 0)) continue;
        int active[6], nActive = 0;
        for (int k = 0; k < 6; ++k)
            if (col[k] > 1e-6 * colMax) active[nActive++] = k;
        if (static_cast<std::size_t>(nActive) > sx.size()) continue;

        scalar gr[6][6] = {}, br[6] = {}, sr[6] = {};
        for (int j = 0; j < nActive; ++j)
        {
            sr[j] = std::sqrt(col[active[j]]);
            br[j] = rhs[active[j]] / sr[j];
            for (int k = 0; k < nActive; ++k)
                gr[j][k] = g[active[j]][active[k]] / (sr[j] * sr[k]);
            gr[j][j] += 1e-8;
        }
        scalar cr[6] = {0, 0, 0, 0, 0, 0};
        if (!solveSpd6(gr, br, cr, nActive)) continue;

        // z = L (c0 x^2 + c1 x y + c2 y^2 + c3 x + c4 y + c5) with x, y the
        // scaled coordinates, so at the origin z_xx = 2 c0 / L, z_xy = c1 / L,
        // z_yy = 2 c2 / L and the slopes are c3 and c4 outright; a dropped
        // direction keeps coefficient 0
        scalar cf[6] = {0, 0, 0, 0, 0, 0};
        for (int j = 0; j < nActive; ++j) cf[active[j]] = cr[j] / sr[j];
        const scalar zxx = 2 * cf[0] / span;
        const scalar zxy = cf[1] / span;
        const scalar zyy = 2 * cf[2] / span;
        const scalar zx = cf[3], zy = cf[4];
        const scalar w2 = 1 + zx * zx + zy * zy;
        // the divergence of the graph's unit normal, which for the cylinder the
        // band wraps is -1/R and for a plane front is 0 - the sign and units of
        // the cell-centred stencil in vof.cpp
        const scalar k = -(zxx * (1 + zy * zy) + zyy * (1 + zx * zx) - 2 * zx * zy * zxy)
                         / (w2 * std::sqrt(w2));
        if (std::isfinite(k))
        {
            kappa[u] = k;
        }
    }
    return kappa;
}

namespace
{

/// The values a cell's vertices take along one axis, merged within `tol`. Two
/// of them and the cell's extent along that axis is its bounding box, which is
/// what the exact disk filling needs.
std::size_t distinctAxisValues(const polyMesh& mesh, label cell, int axis, scalar tol)
{
    std::vector<scalar> seen;
    for (label f : mesh.cells()[static_cast<std::size_t>(cell)])
    {
        for (label v : mesh.faces()[static_cast<std::size_t>(f)])
        {
            const scalar x = mesh.points()[static_cast<std::size_t>(v)][axis];
            bool found = false;
            for (scalar s : seen)
                if (std::abs(s - x) <= tol) { found = true; break; }
            if (!found) seen.push_back(x);
        }
    }
    return seen.size();
}

/// integral of sqrt(R^2 - x^2)
scalar circleAntiderivative(scalar x, scalar radius)
{
    const scalar r2 = radius * radius;
    const scalar t = std::max(-radius, std::min(radius, x));
    const scalar root = std::sqrt(std::max(0.0, r2 - t * t));
    return 0.5 * (t * root + r2 * std::asin(t / radius));
}

} // namespace

Field<scalar> exactDiskAlpha(const polyMesh& mesh, Vector centre, scalar radius)
{
    if (!(radius > 0)) cfFatal("plic: the disk radius must be positive");
    const label nC = mesh.nCells();
    Field<scalar> a(static_cast<std::size_t>(nC), 0.0);

    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        scalar x1 = kInfinity, x2 = -kInfinity, y1 = kInfinity, y2 = -kInfinity;
        for (label f : mesh.cells()[u])
        {
            for (label v : mesh.faces()[static_cast<std::size_t>(f)])
            {
                const Vector& p = mesh.points()[static_cast<std::size_t>(v)];
                x1 = std::min(x1, p.x); x2 = std::max(x2, p.x);
                y1 = std::min(y1, p.y); y2 = std::max(y2, p.y);
            }
        }
        const scalar scale = std::max(x2 - x1, y2 - y1);
        const scalar tol = 1e-9 * (scale > 0 ? scale : 1.0);
        if (distinctAxisValues(mesh, c, 0, tol) != 2
            || distinctAxisValues(mesh, c, 1, tol) != 2)
            cfFatal("plic: exact disk filling needs axis-aligned cells");

        const scalar area = diskRectangleArea(x1 - centre.x, x2 - centre.x,
                                             y1 - centre.y, y2 - centre.y, radius);
        const scalar box = (x2 - x1) * (y2 - y1);
        a[u] = box > 0 ? std::min(1.0, std::max(0.0, area / box)) : 0.0;
    }
    return a;
}

scalar diskRectangleArea(scalar x1, scalar x2, scalar y1, scalar y2, scalar radius)
{
    if (!(radius > 0) || !(y2 > y1) || !(x2 > x1)) return 0;
    const scalar xa = std::max(x1, -radius);
    const scalar xb = std::min(x2, radius);
    if (xb <= xa) return 0;

    // Every branch of the integrand - the top switching between the box side and
    // the arc, the bottom likewise, and the strip emptying where the chord falls
    // short of the nearer side - turns on |x| = sqrt(R^2 - level^2), so the two
    // |levels| are the only places the piecewise integral can split.
    std::vector<scalar> node = {xa, xb};
    auto addBreaks = [&](scalar level) {
        const scalar l = std::abs(level);
        if (l > radius) return;
        const scalar d = std::sqrt(std::max(0.0, radius * radius - l * l));
        node.push_back(d);
        node.push_back(-d);
    };
    addBreaks(y1);
    addBreaks(y2);

    std::sort(node.begin(), node.end());
    std::vector<scalar> knots;
    for (scalar x : node)
    {
        if (x < xa || x > xb) continue;
        if (!knots.empty() && x - knots.back() <= 1e-15 * radius) continue;
        knots.push_back(x);
    }
    if (knots.size() < 2) return 0;

    scalar area = 0;
    for (std::size_t i = 0; i + 1 < knots.size(); ++i)
    {
        const scalar p = knots[i], q = knots[i + 1];
        const scalar m = 0.5 * (p + q);
        const scalar chord = std::sqrt(std::max(0.0, radius * radius - m * m));
        const scalar top = std::min(y2, chord);
        const scalar bottom = std::max(y1, -chord);
        if (top <= bottom) continue;

        const scalar seg = circleAntiderivative(q, radius) - circleAntiderivative(p, radius);
        const scalar intTop = (top == chord) ? seg : y2 * (q - p);
        const scalar intBottom = (bottom == -chord) ? -seg : y1 * (q - p);
        area += intTop - intBottom;
    }
    return std::max(0.0, area);
}

} // namespace cf::vof
