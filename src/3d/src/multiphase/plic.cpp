#include "cf/multiphase/Plic.h"

#include <algorithm>
#include <cmath>
#include <limits>

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

Vector vertexMean(const std::vector<CellFace>& cell)
{
    Vector m = Vector::zero();
    label n = 0;
    for (const CellFace& cf : cell)
        for (const Vector& p : cf.loop) { m += p; ++n; }
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

/// The cell's volume and centroid from its own boundary. This is also the check
/// that the rewound loops close: with no cut at all the same signed-tetrahedra
/// sum has to reproduce the mesh's cell volume.
void fullSolid(const std::vector<CellFace>& cell, scalar& volume, Vector& centroid)
{
    const Vector origin = vertexMean(cell);
    std::vector<std::vector<Vector>> loops;
    loops.reserve(cell.size());
    for (const CellFace& cf : cell) loops.push_back(cf.loop);
    const Solid s = solidOf(loops, origin);
    volume = s.volume;
    centroid = s.volume > 0 ? origin + s.centroid() : Vector::zero();
}

} // namespace

Cut clipConvexCell(const std::vector<CellFace>& cell, const Vector& unitNormal,
                   scalar offset)
{
    if (mag2(unitNormal) <= 0) cfFatal("plic: the cut normal must be a unit vector");
    Cut cut;

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
    const scalar extent = smax - smin;
    const scalar eps = 1e-12 * (extent > 0 ? extent : 1.0);
    cut.areaCentroid = Vector::zero();

    if (offset <= smin - eps) { cut.whole = true; return cut; }
    if (offset >= smax + eps) { cut.void_ = true; return cut; }

    const Vector origin = vertexMean(cell);
    std::vector<std::vector<Vector>> loops;
    for (const CellFace& cf : cell)
    {
        std::vector<Vector> clipped = clipLoop(cf.loop, unitNormal, offset, eps);
        if (clipped.size() >= 3) loops.push_back(std::move(clipped));
    }

    // The cap's corners are the clipped loops' points that lie on the plane. That
    // is the edge crossings, but reading them off the result rather than out of
    // the clip is what also picks up the cell's own vertices when the plane runs
    // exactly through them - a cut through a corner, whose section would
    // otherwise be missed and the kept solid left open. It only applies when the
    // plane is inside the cell: on the cell's own lowest plane the clipped loops
    // already are its closed boundary, and a cap there would be that face twice.
    const scalar capTol = 1e-9 * (extent > 0 ? extent : 1.0);
    if (offset > smin + eps)
    {
        std::vector<Vector> crossings;
        for (const std::vector<Vector>& l : loops)
        {
            for (const Vector& p : l)
                if (std::abs(dot(unitNormal, p) - offset) <= capTol) crossings.push_back(p);
        }
        const std::vector<Vector> cap = orderCapPolygon(std::move(crossings), unitNormal, capTol);
        if (cap.size() >= 3)
        {
            loops.push_back(cap);
            cut.area = mag(loopAreaVector(cap));
            cut.areaCentroid = loopCentroid(cap);
        }
    }

    const Solid s = solidOf(loops, origin);
    if (s.volume <= 0)
    {
        // the plane only grazes the cell: nothing of finite volume is kept
        cut.void_ = true;
        cut.area = 0;
        cut.areaCentroid = Vector::zero();
        return cut;
    }
    cut.volume = s.volume;
    cut.centroid = origin + s.centroid();
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

    scalar volume = 0;
    Vector centroid = Vector::zero();
    fullSolid(cell, volume, centroid);
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
                     0};

    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const scalar a = alpha[u];
        if (!(a > bandTolerance && a < 1.0 - bandTolerance)) continue;

        const Vector g = youngsGradient[u];
        if (mag2(g) <= 0) continue; // no direction to start the cut from
        const Vector n0 = g * (1.0 / mag(g));

        const std::vector<CellFace> cell = cellBoundary(mesh, c);
        scalar volume = 0;
        Vector centroid = Vector::zero();
        fullSolid(cell, volume, centroid);
        if (!(volume > 0)) cfFatal("plic: cell with no volume in the interface band");

        const scalar o0 = offsetForFraction(cell, n0, a);
        Cut phase1 = clipConvexCell(cell, n0, o0);
        Cut phase0 = clipConvexCell(cell, Vector{-n0.x, -n0.y, -n0.z}, -o0);

        Vector n = n0;
        if (choice != NormalChoice::Youngs && !phase1.whole && !phase1.void_
            && !phase0.whole && !phase0.void_)
        {
            // from the centroid of what is not phase 1 to the centroid of what
            // is: the direction the phase actually increases along, measured on
            // the sub-cell geometry rather than on the field's slope
            const Vector cd = phase1.centroid - phase0.centroid;
            if (mag2(cd) > 0)
            {
                const Vector u1 = cd * (1.0 / mag(cd));
                if (choice == NormalChoice::Centroid)
                {
                    n = u1;
                }
                else
                {
                    const Vector m = u1 + n0;
                    if (mag2(m) > 0) n = m * (1.0 / mag(m));
                }
            }
        }

        // the plane the final normal cuts is the reconstruction that is reported
        const scalar off = offsetForFraction(cell, n, a);
        const Cut kept = clipConvexCell(cell, n, off);
        if (kept.whole || kept.void_ || !(kept.volume > 0))
            cfFatal("plic: the reconstruction lost the interface it was given");

        r.normal[u] = n;
        r.offset[u] = off;
        r.area[u] = kept.area;
        r.areaCentroid[u] = kept.areaCentroid;
        r.fraction[u] = kept.volume / volume;
        ++r.bandCells;
    }
    return r;
}

Reconstruction reconstruct(const volScalarField& alpha, const std::string& gradScheme,
                           NormalChoice choice, scalar bandTolerance)
{
    const Field<Vector> g = fvc::grad(alpha, gradScheme);
    return reconstruct(alpha.mesh(), alpha.internalField(), g, choice, bandTolerance);
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
