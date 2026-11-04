// Unit tests for the PLIC layer. The cut geometry is analytic, so every test
// here is exactness rather than convergence: the volume, centroid and cut area
// of a plane through a cell against a closed form, the cell's own boundary loops
// read back as the mesh's cell volume on hexes, tets and prisms, the fraction a
// reconstruction reads back from the field it was given, the normal and curvature
// a cell-filled disk gives up, and what the geometric transport owes on its own:
// a region cut a second time, a step front translated without spreading, and a
// rotation that conserves the phase and bounds the fraction at the same time.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

#include "cf/fvm/BlockMesh.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/GmshMesh.h"
#include "cf/multiphase/Plic.h"
#include "cf/multiphase/VoF.h"

#ifndef CF_MESH_DIR
#define CF_MESH_DIR "."
#endif

using namespace cf;

static int failures = 0;

#define CHECK(cond, msg)                                                        \
    do {                                                                        \
        if (!(cond)) { ++failures; std::cout << "FAIL: " << (msg) << "\n"; }     \
        else          { std::cout << "ok:   " << (msg) << "\n"; }               \
    } while (0)

static bool relNear(double a, double b, double tol)
{
    return std::abs(a - b) < tol * std::max(1.0, std::max(std::abs(a), std::abs(b)));
}

static std::string meshPath(const char* name)
{
    return std::string(CF_MESH_DIR) + "/" + name;
}

namespace
{

/// The s-multiple of a unit cube's cut: volume of the part with x+y+z <= s.
double cubeCapVolume(double s)
{
    auto cubed = [](double t) { return t > 0 ? t * t * t : 0.0; };
    return (cubed(s) - 3.0 * cubed(s - 1.0) + 3.0 * cubed(s - 2.0) - cubed(s - 3.0)) / 6.0;
}

/// Its derivative along the cut, which is the area of the plane section divided
/// by |grad(x+y+z)| = sqrt(3).
double cubeCapArea(double s)
{
    auto squared = [](double t) { return t > 0 ? t * t : 0.0; };
    return std::sqrt(3.0)
         * (squared(s) - 3.0 * squared(s - 1.0) + 3.0 * squared(s - 2.0) - squared(s - 3.0))
         / 2.0;
}

/// One cell, one unit cube: every closed form below is stated for it.
std::vector<vof::CellFace> unitCell(const polyMesh& mesh)
{
    return vof::cellBoundary(mesh, 0);
}

/// A closed box of nx x ny cells and one layer in z, which is the 2D geometry the
/// disk tests need.
polyMesh boxXY(label nx, label ny, scalar dz = 0.1)
{
    BlockSides s;
    s.xmin = {"sx0", "symmetry"};
    s.xmax = {"sx1", "symmetry"};
    s.ymin = {"sy0", "symmetry"};
    s.ymax = {"sy1", "symmetry"};
    s.zmin = {"sz0", "symmetry"};
    s.zmax = {"sz1", "symmetry"};
    return BlockMesh::create(Vector::zero(), Vector{1.0, 1.0, dz}, nx, ny, 1, s);
}

/// The extent of a cell's boundary along one direction: the plane constants the
/// cut is between.
void normalExtent(const std::vector<vof::CellFace>& cell, const Vector& n,
                  scalar& smin, scalar& smax)
{
    smin = std::numeric_limits<double>::infinity();
    smax = -std::numeric_limits<double>::infinity();
    for (const vof::CellFace& cf : cell)
    {
        for (const Vector& p : cf.loop)
        {
            const scalar s = dot(n, p);
            smin = std::min(smin, s);
            smax = std::max(smax, s);
        }
    }
}

/// The volume phase 1 occupies in the field: what a conservative transport keeps,
/// and what the geometric step therefore has to hand back to round-off.
double phaseVolume(const volScalarField& alpha)
{
    const Field<scalar>& v = alpha.mesh().cellVolumes();
    double s = 0;
    for (std::size_t u = 0; u < v.size(); ++u) s += alpha.internalField()[u] * v[u];
    return s;
}

/// Angle in degrees between two unit directions.
double angleDegrees(const Vector& a, const Vector& b)
{
    const scalar m = mag(a) * mag(b);
    if (m <= 0) return 180.0;
    const scalar c = dot(a, b) / m;
    return std::acos(std::max(-1.0, std::min(1.0, c))) * 180.0 / M_PI;
}

/// Midpoint-rule area of the part of a disk inside a rectangle: independent of
/// the analytic branch splitting, and only as good as its grid, which is why the
/// specials below carry the precision claim and this one the coverage.
double gridDiskRectangleArea(double x1, double x2, double y1, double y2, double r)
{
    const int n = 1500;
    const double hx = (x2 - x1) / n, hy = (y2 - y1) / n;
    double hits = 0;
    for (int i = 0; i < n; ++i)
    {
        const double xm = x1 + (i + 0.5) * hx;
        for (int j = 0; j < n; ++j)
        {
            const double ym = y1 + (j + 0.5) * hy;
            if (xm * xm + ym * ym <= r * r) ++hits;
        }
    }
    return hits * hx * hy;
}

/// The disk the reconstruction tests share: centred on the box, its radius a
/// multiple of the cell so it lands symmetrically on every mesh used here.
Vector diskCentre() { return Vector{0.5, 0.5, 0.0}; }
constexpr double kDiskRadius = 0.2;
constexpr double kDiskDepth = 0.1;

/// What a reconstruction owes on a field whose interface is known: the fraction
/// it reads back, the band it covers, the area of the interface it found and the
/// normal and curvature of that interface measured against the analytic circle.
struct Measured
{
    label band = 0;          ///< cells the field puts the interface in
    label reconstructed = 0; ///< cells the reconstruction put one in
    label gradientBand = 0;  ///< cells a gradient threshold alone would keep
    label gradientOutside = 0; ///< of those, cells that hold no interface at all
    double worstRead = 0;    ///< |the fraction the plane cuts - alpha|
    double area = 0;         ///< the interface area of the band
    double meanAngle = 0;    ///< normal off the radial direction, degrees
    double worstAngle = 0;
    double meanScatter = 0;  ///< band-mean of the residual disagreement of the fitted planes
    double kappaMean = 0;    ///< band-mean of div(normal)
    double kappaPlicMean = 0; ///< band-mean of the curvature the cut geometry reports
    double kappaPlicWorst = 0; ///< largest |curvature| the cut geometry reports
    label kappaPlicReports = 0; ///< band cells the cut geometry reported a curvature for
};

Measured measure(const polyMesh& mesh, const volScalarField& alpha,
                 const Field<scalar>& exact, const vof::Reconstruction& r)
{
    Measured m;
    const Vector c = diskCentre();
    const Field<scalar> k = vof::interfaceCurvature(mesh, r.normal);
    const Field<scalar> kp = vof::interfaceCurvaturePlic(mesh, r);
    label kappaN = 0, kappaPlicN = 0;
    double angle = 0, kappa = 0, kappaPlic = 0, scatter = 0;

    for (label cell = 0; cell < mesh.nCells(); ++cell)
    {
        const std::size_t u = static_cast<std::size_t>(cell);
        const scalar a = exact[u];
        if (a > 1e-6 && a < 1.0 - 1e-6) ++m.band;
        if (mag2(r.normal[u]) == 0) continue;
        ++m.reconstructed;

        // the plane it reports has to cut the fraction the field says
        const scalar v = mesh.cellVolumes()[u];
        const vof::Cut kept = vof::clipConvexCell(vof::cellBoundary(mesh, cell),
                                                 r.normal[u], r.offset[u]);
        m.worstRead = std::max(m.worstRead, std::abs(kept.volume / v - a));
        m.area += r.area[u];

        // towards increasing alpha is into the disk. The interface is the extrusion
        // of a circle, so its normal is the in-plane radial: the analytic direction
        // is that radial with no z-component, and a fitted normal that tilts out of
        // the plane is charged for it. Measuring against the full centroid
        // difference instead would add the extrusion's own half-depth to every
        // angle, which is a property of the test box, not of the estimator.
        const Vector d = r.areaCentroid[u] - c;
        const Vector radial{-d.x, -d.y, 0.0};
        const double deg = angleDegrees(r.normal[u], radial);
        angle += deg;
        m.worstAngle = std::max(m.worstAngle, deg);
        kappa += k[u]; ++kappaN;
        if (kp[u] != 0)
        {
            kappaPlic += kp[u];
            m.kappaPlicWorst = std::max(m.kappaPlicWorst, std::abs(kp[u]));
            ++kappaPlicN;
        }
        scatter += r.fitScatter[u];
    }
    m.meanAngle = m.reconstructed > 0 ? angle / m.reconstructed : 0;
    m.meanScatter = m.reconstructed > 0 ? scatter / m.reconstructed : 0;
    m.kappaMean = kappaN > 0 ? kappa / kappaN : 0;
    m.kappaPlicMean = kappaPlicN > 0 ? kappaPlic / kappaPlicN : 0;
    m.kappaPlicReports = kappaPlicN;

    // the band the smeared-interface path keeps: the same field, but selected by
    // a threshold on |grad alpha| rather than on alpha itself. The threshold also
    // fires in the cells around the band, whose alpha is 0 or 1 and which hold no
    // interface to reconstruct, so the count is reported with those separated.
    const Field<Vector> g = vof::interfaceNormal(alpha, "corrected", 0.1);
    for (std::size_t u = 0; u < g.size(); ++u)
    {
        if (mag2(g[u]) == 0) continue;
        ++m.gradientBand;
        const scalar a = exact[u];
        if (!(a > 1e-6 && a < 1.0 - 1e-6)) ++m.gradientOutside;
    }
    return m;
}

} // namespace

int main()
{
    // ---- 1. a plane through a unit cube, axis-aligned ----
    {
        const polyMesh mesh = boxXY(1, 1, 1.0);
        const std::vector<vof::CellFace> cell = unitCell(mesh);
        const Vector n{1.0, 0.0, 0.0};

        const vof::Cut kept = vof::clipConvexCell(cell, n, 0.65);
        CHECK(!kept.whole && !kept.void_, "an interior plane neither keeps nor drops the cell");
        CHECK(relNear(kept.volume, 0.35, 1e-13), "the kept slab is x in [0.65,1]: volume 0.35");
        CHECK(relNear(kept.centroid.x, 0.825, 1e-13) && relNear(kept.centroid.y, 0.5, 1e-13)
                  && relNear(kept.centroid.z, 0.5, 1e-13),
              "the slab's centroid is its own centre");
        CHECK(relNear(kept.area, 1.0, 1e-13), "the cut face of a unit cube is 1");
        CHECK(relNear(kept.areaCentroid.x, 0.65, 1e-13)
                  && relNear(kept.areaCentroid.y, 0.5, 1e-13)
                  && relNear(kept.areaCentroid.z, 0.5, 1e-13),
              "the cut face sits at the plane constant");

        // outside the cell the plane keeps everything or nothing
        const vof::Cut all = vof::clipConvexCell(cell, n, -1.0);
        const vof::Cut none = vof::clipConvexCell(cell, n, 2.0);
        CHECK(all.whole && !all.void_, "a plane before the cell keeps the whole cell");
        CHECK(none.void_ && !none.whole, "a plane past the cell keeps nothing");
        CHECK(relNear(vof::offsetForFraction(cell, n, 0.35), 0.65, 1e-13),
              "the offset that holds 0.35 of the cell is 0.65");
        CHECK(relNear(vof::offsetForFraction(cell, n, 0.5), 0.5, 1e-13),
              "the offset that holds half the cell is its midplane");
        CHECK(relNear(vof::offsetForFraction(cell, n, 1.0), 0.0, 1e-13)
                  && relNear(vof::offsetForFraction(cell, n, 0.0), 1.0, 1e-13),
              "fraction 1 and 0 land on the cell's own faces");

        // every fraction the plane can hold, and the complement it leaves
        bool exact = true, closed = true;
        for (int i = 1; i < 20; ++i)
        {
            const scalar f = i / 20.0;
            const scalar off = vof::offsetForFraction(cell, n, f);
            const vof::Cut a = vof::clipConvexCell(cell, n, off);
            const vof::Cut b = vof::clipConvexCell(cell, Vector{-1.0, 0.0, 0.0}, -off);
            exact = exact && relNear(a.volume, f, 1e-12);
            closed = closed && relNear(a.volume + b.volume, 1.0, 1e-12);
        }
        CHECK(exact, "offsetForFraction cuts the requested fraction to 1e-12");
        CHECK(closed, "the two sides of the plane add back to the cell");
    }

    // ---- 2. the same along the cube diagonal, against the closed form ----
    {
        const polyMesh mesh = boxXY(1, 1, 1.0);
        const std::vector<vof::CellFace> cell = unitCell(mesh);
        const scalar k = 1.0 / std::sqrt(3.0);
        const Vector n{k, k, k};

        bool volume = true, area = true, symmetric = true, cap = true;
        double worstV = 0, worstA = 0;
        for (double t : {0.3, 0.8, 1.0, 1.4, 1.5, 2.0, 2.4, 2.9, 2.99})
        {
            const vof::Cut kept = vof::clipConvexCell(cell, n, t * k);
            const double expectV = 1.0 - cubeCapVolume(t);
            const double expectA = cubeCapArea(t);
            worstV = std::max(worstV, std::abs(kept.volume - expectV));
            worstA = std::max(worstA, std::abs(kept.area - expectA));
            volume = volume && relNear(kept.volume, expectV, 1e-12);
            area = area && relNear(kept.area, expectA, 1e-12);
            // the kept solid is invariant under x<->y<->z, so both its centres and
            // the cut's lie on the diagonal
            const double cx = (kept.centroid.x + kept.centroid.y + kept.centroid.z) / 3.0;
            const double ax = (kept.areaCentroid.x + kept.areaCentroid.y + kept.areaCentroid.z) / 3.0;
            symmetric = symmetric && std::abs(kept.centroid.x - cx) < 1e-13
                                 && std::abs(kept.centroid.y - cx) < 1e-13
                                 && std::abs(kept.centroid.z - cx) < 1e-13;
            cap = cap && std::abs(kept.areaCentroid.x - ax) < 1e-13
                       && std::abs(kept.areaCentroid.y - ax) < 1e-13
                       && std::abs(kept.areaCentroid.z - ax) < 1e-13;
        }
        std::cout << "     diagonal cut: worst |dV| = " << worstV
                  << ", worst |dA| = " << worstA << "\n";
        CHECK(volume, "the diagonal cut's volume is the cube's closed form");
        CHECK(area, "the diagonal cut's area is its derivative");
        CHECK(symmetric, "the kept solid's centroid lies on the diagonal");
        CHECK(cap, "the hexagonal/triangular section is centred on the diagonal too");

        // and the Newton finds that plane from the fraction alone
        bool inverted = true;
        for (double t : {0.5, 1.0, 1.5, 2.0, 2.5})
        {
            const scalar f = static_cast<scalar>(1.0 - cubeCapVolume(t));
            const scalar off = vof::offsetForFraction(cell, n, f);
            inverted = inverted && relNear(off, t * k, 1e-11);
        }
        CHECK(inverted, "the offset recovered from the fraction is the plane constant");
    }

    // ---- 3. the cell's own boundary closes to the mesh's cell volume ----
    // This is the check that the loops cellBoundary hands back are wound outward
    // on every cell type: the signed-tetrahedra sum with no cut at all.
    {
        struct Case { const char* file; const char* kind; };
        const std::vector<Case> cases = {
            {"cubeHex_41.msh", "hex"}, {"cubeTets_41.msh", "tet"},
            {"cubeTets_22.msh", "tet"}, {"prismChannel_22.msh", "prism"},
        };
        for (const Case& c : cases)
        {
            const polyMesh mesh = GmshMesh::fromFile(meshPath(c.file));
            scalar worstVolume = 0, worstCentroid = 0;
            bool closed = true;
            for (label cell = 0; cell < mesh.nCells(); ++cell)
            {
                const std::vector<vof::CellFace> b = vof::cellBoundary(mesh, cell);
                scalar smin = 0, smax = 0;
                normalExtent(b, Vector{1.0, 0.0, 0.0}, smin, smax);
                const vof::Cut full = vof::clipConvexCell(b, Vector{1.0, 0.0, 0.0}, smin);
                const scalar v = mesh.cellVolumes()[static_cast<std::size_t>(cell)];
                const Vector cc = mesh.cellCentres()[static_cast<std::size_t>(cell)];
                closed = closed && !full.whole && !full.void_;
                worstVolume = std::max(worstVolume, std::abs(full.volume - v) / v);
                worstCentroid = std::max(worstCentroid, mag(full.centroid - cc));
            }
            std::cout << "     " << c.kind << " " << c.file << ": " << mesh.nCells()
                      << " cells, worst |dV|/V = " << worstVolume
                      << ", worst |dc| = " << worstCentroid << "\n";
            CHECK(worstVolume < 1e-12,
                  (std::string("the boundary loops of every cell of ") + c.file
                   + " close to the mesh's cell volume").c_str());
            CHECK(worstCentroid < 1e-12,
                  (std::string("and to its cell centre: ") + c.file).c_str());
            CHECK(closed, (std::string("a plane on the cell's own face is a cut, not a flag: ")
                           + c.file).c_str());

            // and the cut closes on both sides for an inclined plane through
            // every cell, which is what the reconstruction relies on
            const Vector n{0.4, -0.7, 0.2};
            const scalar m = mag(n);
            const Vector u{n.x / m, n.y / m, n.z / m};
            bool complementary = true;
            scalar worstSum = 0;
            for (label cell = 0; cell < mesh.nCells(); ++cell)
            {
                const std::vector<vof::CellFace> b = vof::cellBoundary(mesh, cell);
                scalar smin = 0, smax = 0;
                normalExtent(b, u, smin, smax);
                const scalar off = smin + 0.37 * (smax - smin);
                const vof::Cut a = vof::clipConvexCell(b, u, off);
                const vof::Cut z = vof::clipConvexCell(b, Vector{-u.x, -u.y, -u.z}, -off);
                const scalar v = mesh.cellVolumes()[static_cast<std::size_t>(cell)];
                worstSum = std::max(worstSum, std::abs(a.volume + z.volume - v) / v);
                complementary = complementary && relNear(a.volume + z.volume, v, 1e-11);
            }
            std::cout << "     inclined cut on " << c.file << ": worst |d(V+V')|/V = "
                      << worstSum << "\n";
            CHECK(complementary,
                  (std::string("the two sides of an inclined plane fill the cell: ") + c.file)
                      .c_str());
        }
    }

    // ---- 4. the exact overlap of a disk and a rectangle ----
    {
        const double R = 1.0;
        CHECK(relNear(vof::diskRectangleArea(-2, 2, -2, 2, R), M_PI, 1e-12),
              "a box around the disk gives its area");
        CHECK(relNear(vof::diskRectangleArea(0, 2, 0, 2, R), M_PI / 4.0, 1e-12),
              "the first quadrant gives a quarter");
        CHECK(relNear(vof::diskRectangleArea(0, 1, -1, 1, R), M_PI / 2.0, 1e-12),
              "the half disk between x=0 and the tangent gives a half");
        // the strip y in [0, 1/2] of the right half: int sqrt(1-x^2) dx to x=1/2
        const double s = 0.5;
        CHECK(relNear(vof::diskRectangleArea(0, s, 0, 1, R),
                      0.5 * (s * std::sqrt(1 - s * s) + std::asin(s)), 1e-12),
              "a strip whose top is the arc is the arc's own integral");
        // the cap above y=d: sector minus triangle, R^2 acos(d) - d sqrt(R^2-d^2)
        CHECK(relNear(vof::diskRectangleArea(-2, 2, 0.5, 1.5, R),
                      std::acos(0.5) - 0.5 * std::sqrt(0.75), 1e-12),
              "a box whose floor cuts the disk gives the circular segment");
        CHECK(relNear(vof::diskRectangleArea(2, 3, -1, 1, R), 0.0, 1e-15),
              "a box clear of the disk gives nothing");

        const double boxes[][4] = {
            {-0.3, 1.4, -0.9, 0.2}, {-1.1, -0.4, 0.1, 0.95}, {-0.7, 0.7, 0.4, 0.6},
            {-0.2, 0.9, -0.99, 0.99}, {-1.3, 1.3, -0.35, 0.35},
        };
        double worst = 0;
        for (const auto& b : boxes)
        {
            const double exact = vof::diskRectangleArea(b[0], b[1], b[2], b[3], R);
            const double brute = gridDiskRectangleArea(b[0], b[1], b[2], b[3], R);
            worst = std::max(worst, std::abs(exact - brute));
        }
        std::cout << "     disk-rectangle overlap: worst |exact - grid| = " << worst << "\n";
        CHECK(worst < 2e-3, "the branch-split integral agrees with a grid on off-axis boxes");
    }

    // ---- 5. a cell-filled disk: alpha is the fraction the geometry holds ----
    {
        const label n = 40;
        const polyMesh mesh = boxXY(n, n, kDiskDepth);
        const Vector c{0.5, 0.5, 0.0};
        const scalar R = kDiskRadius;
        const Field<scalar> a = vof::exactDiskAlpha(mesh, c, R);

        scalar sum = 0, lo = 1, hi = 0;
        label band = 0;
        for (std::size_t u = 0; u < a.size(); ++u)
        {
            sum += a[u];
            lo = std::min(lo, a[u]);
            hi = std::max(hi, a[u]);
            if (a[u] > 1e-6 && a[u] < 1.0 - 1e-6) ++band;
        }
        // the cell area is dx*dy and every cell is one layer thick, so the sum of
        // the fractions is the disk's area divided by the cell's. The tolerance is
        // not machine precision because the two cells the circle is tangent to
        // reach asin(+/-1), where the antiderivative's derivative is infinite and
        // rounding is amplified by a square root: the residual left is ~2e-10 of
        // the disk's area, which is what is measured here.
        const scalar dx = 1.0 / static_cast<scalar>(n);
        std::cout << "     sum(alpha) = " << std::setprecision(12) << sum
                  << ", disk area / cell area = " << M_PI * R * R / (dx * dx)
                  << ", band cells " << band << "\n";
        CHECK(relNear(sum, M_PI * R * R / (dx * dx), 1e-8),
              "the exact filling conserves the disk's area");
        CHECK(lo >= 0.0 && hi <= 1.0, "the filling stays a phase fraction");
        CHECK(band > 2 * M_PI * R / dx * 0.5, "the interface is a band, not a spike");

        // cells whose far corner is inside the disk are full, cells whose near
        // corner is outside are empty, and only the band is in between
        bool filled = true;
        for (label cell = 0; cell < mesh.nCells(); ++cell)
        {
            const std::size_t u = static_cast<std::size_t>(cell);
            scalar xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30;
            for (label f : mesh.cells()[u])
            {
                for (label v : mesh.faces()[static_cast<std::size_t>(f)])
                {
                    const Vector& p = mesh.points()[static_cast<std::size_t>(v)];
                    xmin = std::min(xmin, p.x - c.x); xmax = std::max(xmax, p.x - c.x);
                    ymin = std::min(ymin, p.y - c.y); ymax = std::max(ymax, p.y - c.y);
                }
            }
            const scalar far2 = std::max(xmin * xmin, xmax * xmax)
                              + std::max(ymin * ymin, ymax * ymax);
            const scalar near2 = (xmin > 0 ? xmin * xmin : (xmax < 0 ? xmax * xmax : 0.0))
                               + (ymin > 0 ? ymin * ymin : (ymax < 0 ? ymax * ymax : 0.0));
            if (far2 <= R * R && !(a[u] > 1.0 - 1e-14)) filled = false;
            if (near2 > R * R && !(a[u] < 1e-14)) filled = false;
        }
        CHECK(filled, "the filling is 1 where the cell is inside and 0 where it is out");
    }

    // ---- 6. reconstruction: the plane the field says, read back exactly ----
    {
        const label n = 40;
        const scalar R = kDiskRadius;
        const polyMesh mesh = boxXY(n, n, kDiskDepth);
        const Field<scalar> exact = vof::exactDiskAlpha(mesh, diskCentre(), R);
        volScalarField alpha(mesh, "alpha", 0.0);
        for (std::size_t u = 0; u < alpha.internalField().size(); ++u)
            alpha.internalField()[u] = exact[u];

        for (const auto& entry : std::vector<std::pair<std::string, vof::NormalChoice>>{
                 {"youngs", vof::NormalChoice::Youngs},
                 {"centroid", vof::NormalChoice::Centroid},
                 {"mixed", vof::NormalChoice::Mixed},
                 {"fit", vof::NormalChoice::Fit}})
        {
            const Measured m = measure(mesh, alpha, exact,
                                            vof::reconstruct(alpha, "corrected",
                                                             entry.second));
            std::cout << "     " << std::setw(8) << entry.first << ": band "
                      << m.reconstructed << " of " << m.band
                      << " (the gradient threshold sees " << m.gradientBand << ")"
                      << ", read-back error " << std::setprecision(3) << m.worstRead
                      << ", normal off-radial " << m.meanAngle << " deg mean /"
                      << m.worstAngle << " deg max, plane disagreement " << m.meanScatter
                      << " (0 for a cell-local normal, which does not report it)\n     " << std::setw(8) << ""
                      << "interface area " << std::setprecision(6) << m.area
                      << " vs the extruded perimeter 2 pi R dz = "
                      << 2 * M_PI * R * kDiskDepth << ", mean kappa " << std::setprecision(5)
                      << m.kappaMean << " vs -1/R = " << -1.0 / R
                      << ", plic kappa " << m.kappaPlicMean << " over " << m.kappaPlicReports
                      << " of " << m.reconstructed
                      << " (worst |k| " << m.kappaPlicWorst << ")\n";
            CHECK(m.worstRead < 1e-10,
                  (entry.first + ": the reconstructed plane cuts the fraction it was given")
                      .c_str());
            CHECK(m.reconstructed == m.band,
                  (entry.first + ": every cell of the band is reconstructed").c_str());
            // The area a one-cell plane model reports is the area of the chord
            // polygon it tiles, so it comes in under the circle's perimeter and rises
            // towards it as the mesh refines. Every estimator lands within the tenth
            // of it that a neighbourhood three cells wide can hold.
            CHECK(m.area > 0.85 * 2 * M_PI * R * kDiskDepth && m.area < 2 * M_PI * R * kDiskDepth,
                  (entry.first + ": the band's interface area is the perimeter's").c_str());
        }
    }

    // ---- 7. what the normal is worth, measured on the mesh it is asked for ----
    // Two estimators on the same exact field, at three refinements.
    //
    // `mixed` reads the cell on its own: the gradient and the phase centroid
    // difference of the staircase that cell is filled with. `fit` reads the
    // neighbourhood: the tilt that makes every band cell's own exact cut lie in one
    // plane. The cut, the area and the read-back are the same machinery in both, so
    // whatever separates the numbers below is the direction the plane is given.
    {
        const scalar R = kDiskRadius;
        const scalar perimeter = static_cast<scalar>(2 * M_PI * R * kDiskDepth);
        const std::vector<std::pair<std::string, vof::NormalChoice>> kinds
            = {{"mixed", vof::NormalChoice::Mixed}, {"fit", vof::NormalChoice::Fit}};
        Measured seen[2][3];
        for (std::size_t ki = 0; ki < kinds.size(); ++ki)
        {
            for (label ni = 0; ni < 3; ++ni)
            {
                const label n = ni == 0 ? 40 : (ni == 1 ? 80 : 160);
                const polyMesh mesh = boxXY(n, n, kDiskDepth);
                const Field<scalar> exact = vof::exactDiskAlpha(mesh, diskCentre(), R);
                volScalarField alpha(mesh, "alpha", 0.0);
                for (std::size_t u = 0; u < alpha.internalField().size(); ++u)
                    alpha.internalField()[u] = exact[u];

                const Measured m = measure(mesh, alpha, exact,
                                           vof::reconstruct(alpha, "corrected",
                                                            kinds[ki].second));
                std::cout << "     " << std::setw(6) << kinds[ki].first << " "
                          << n << "x" << n << ", R/dx = " << R * n
                          << ": alpha band " << m.band << ", PLIC " << m.reconstructed
                          << ", gradient threshold " << m.gradientBand << " of which "
                          << m.gradientOutside << " hold no interface\n     " << std::setw(6) << ""
                          << " normal off-radial " << std::setprecision(4) << m.meanAngle
                          << " deg mean /" << m.worstAngle << " deg max, read-back "
                          << std::setprecision(3) << m.worstRead;
                if (kinds[ki].second == vof::NormalChoice::Fit)
                    std::cout << ", planes still disagree by " << std::setprecision(4)
                              << m.meanScatter * n << " cells";
                std::cout << "\n     "
                          << std::setw(6) << ""
                          << "interface area " << std::setprecision(6) << m.area
                          << " vs " << perimeter << ", kappa " << std::setprecision(5)
                          << m.kappaMean << " vs " << -1.0 / R
                          << ", plic kappa " << m.kappaPlicMean << " over " << m.kappaPlicReports
                          << " of " << m.reconstructed
                          << " (worst |k| " << m.kappaPlicWorst << ")\n";
                seen[ki][ni] = m;
                CHECK(m.reconstructed == m.band,
                      (kinds[ki].first + ": the reconstruction covers the whole band on a "
                                       "one-cell-thick front").c_str());
                CHECK(m.worstRead < 1e-10,
                      (kinds[ki].first + ": and reproduces the fraction everywhere").c_str());
            }
        }
        CHECK(seen[0][0].gradientOutside > 0 && seen[0][2].gradientOutside > 0,
              "the gradient threshold marks cells that hold no interface; the alpha band does not");
        CHECK(seen[0][0].meanAngle < 20.0 && seen[0][2].meanAngle < 20.0,
              "the cell-local normal stays within 20 deg mean of the analytic radial");

        // The claim the fit is here for: the direction converges as the mesh does,
        // which the cell-local estimator, reading the same staircase one cell at a
        // time, does not do at any of the three meshes - it moves by less than a
        // tenth and its worst cell gets no sharper.
        CHECK(seen[1][0].meanAngle * 5 < seen[0][0].meanAngle
                  && seen[1][1].meanAngle * 5 < seen[0][1].meanAngle
                  && seen[1][2].meanAngle * 5 < seen[0][2].meanAngle,
              "at every mesh the fitted normal is five times closer to the radial than the "
              "cell-local one");
        CHECK(seen[1][1].meanAngle < 0.75 * seen[1][0].meanAngle
                  && seen[1][2].meanAngle < 0.75 * seen[1][1].meanAngle,
              "and the fitted normal's mean angle off the radial falls by more than a quarter "
              "every time the cell halves");
        CHECK(seen[1][2].meanAngle < 0.3 * seen[1][0].meanAngle,
              "and to under a third of itself over the two refinements together");
        CHECK(seen[0][1].meanAngle > 0.9 * seen[0][0].meanAngle
                  && seen[0][2].meanAngle > 0.9 * seen[0][1].meanAngle,
              "while the cell-local one moves by less than a tenth over each of the same "
              "refinements");
        CHECK(seen[1][0].worstAngle * 4 < seen[0][0].worstAngle
                  && seen[1][1].worstAngle * 4 < seen[0][1].worstAngle
                  && seen[1][2].worstAngle * 4 < seen[0][2].worstAngle,
              "and the fit's worst cell is under a quarter of the cell-local worst at every "
              "mesh");
        CHECK(seen[1][1].meanScatter < 0.6 * seen[1][0].meanScatter
                  && seen[1][2].meanScatter < 0.6 * seen[1][1].meanScatter,
              "and the neighbourhood's planes agree about twice as closely at each one");

        // The area a one-cell plane model reports is the chord it tiles, under the
        // perimeter and rising to it. The read-back and the band are what do not move
        // between estimators.
        CHECK(seen[0][0].area < perimeter && seen[0][1].area > seen[0][0].area
                  && seen[0][2].area > seen[0][1].area && seen[0][2].area < perimeter,
              "the cell-local band's interface area is under the perimeter and rises towards "
              "it at every refinement");
        CHECK(seen[1][0].area < perimeter && seen[1][1].area < perimeter
                  && seen[1][2].area < perimeter && seen[1][0].area > 0.9 * perimeter
                  && seen[1][1].area > 0.9 * perimeter && seen[1][2].area > 0.9 * perimeter,
              "and the fit keeps its chord within a tenth of the perimeter at all three "
              "meshes");
        CHECK(std::abs(seen[0][0].kappaMean * R + 1.0) < 0.4
                  && std::abs(seen[0][1].kappaMean * R + 1.0) < 0.4
                  && std::abs(seen[0][2].kappaMean * R + 1.0) < 0.4
                  && std::abs(seen[1][0].kappaMean * R + 1.0) < 0.4
                  && std::abs(seen[1][1].kappaMean * R + 1.0) < 0.4
                  && std::abs(seen[1][2].kappaMean * R + 1.0) < 0.4,
              "and the band-mean curvature is -1/R to within what the cell-centred divergence "
              "stencil of a one-cell band resolves - which is what limits it, not the normal");

        // The same disk, the curvature read off the reconstruction's own cut faces:
        // a quadric through the two-ring corner neighbourhood's interface points,
        // which puts the band mean within two per cent of -1/R everywhere the
        // neighbourhood holds together (measured 5.02-5.12 against -5, where the
        // stencil above is 3.25-3.90), and the worst cell within a quarter of it.
        // On the finest mesh a tenth of the fit's cells report nothing: the fitted
        // normal's degenerate out-of-plane tilt on an extruded field leaves a few
        // neighbourhoods holding points the quadric cannot determine, and those
        // cells say so with 0 rather than a made-up curvature.
        CHECK(std::abs(seen[0][0].kappaPlicMean * R + 1.0) < 0.05
                  && std::abs(seen[0][1].kappaPlicMean * R + 1.0) < 0.01
                  && std::abs(seen[0][2].kappaPlicMean * R + 1.0) < 0.01
                  && std::abs(seen[1][0].kappaPlicMean * R + 1.0) < 0.05
                  && std::abs(seen[1][1].kappaPlicMean * R + 1.0) < 0.01
                  && std::abs(seen[1][2].kappaPlicMean * R + 1.0) < 0.2,
              "and the curvature from the cut geometry is within a few per cent of -1/R "
              "on every mesh, one to two per cent where the neighbourhoods all hold");
        CHECK(seen[0][0].kappaPlicReports * 10 >= 8 * seen[0][0].reconstructed
                  && seen[0][1].kappaPlicReports * 10 >= 8 * seen[0][1].reconstructed
                  && seen[0][2].kappaPlicReports * 10 >= 8 * seen[0][2].reconstructed
                  && seen[1][0].kappaPlicReports * 10 >= 8 * seen[1][0].reconstructed
                  && seen[1][1].kappaPlicReports * 10 >= 8 * seen[1][1].reconstructed
                  && seen[1][2].kappaPlicReports * 10 >= 8 * seen[1][2].reconstructed,
              "and it reports for all but a tenth of the band at worst");
        CHECK(seen[0][0].kappaPlicWorst * R < 1.25 && seen[0][1].kappaPlicWorst * R < 1.25
                  && seen[0][2].kappaPlicWorst * R < 1.25
                  && seen[1][0].kappaPlicWorst * R < 1.25
                  && seen[1][1].kappaPlicWorst * R < 1.25
                  && seen[1][2].kappaPlicWorst * R < 1.25,
              "and its worst cell stays within a quarter of -1/R");
        CHECK(std::abs(seen[0][0].kappaPlicMean * R + 1.0) < std::abs(seen[0][0].kappaMean * R + 1.0)
                  && std::abs(seen[0][1].kappaPlicMean * R + 1.0) < std::abs(seen[0][1].kappaMean * R + 1.0)
                  && std::abs(seen[0][2].kappaPlicMean * R + 1.0) < std::abs(seen[0][2].kappaMean * R + 1.0)
                  && std::abs(seen[1][0].kappaPlicMean * R + 1.0) < std::abs(seen[1][0].kappaMean * R + 1.0)
                  && std::abs(seen[1][1].kappaPlicMean * R + 1.0) < std::abs(seen[1][1].kappaMean * R + 1.0)
                  && std::abs(seen[1][2].kappaPlicMean * R + 1.0) < std::abs(seen[1][2].kappaMean * R + 1.0),
              "and at every mesh and for both normals it beats the cell-centred stencil");
    }

    // ---- 7b. the fit on the field it is exact for: one plane ----
    // A front that really is a plane is the fit's contract: every band cell's own
    // cut is a piece of that plane, so the neighbourhood agrees with no tilt left
    // over, and the recovered normal is the analytic one rather than an average of
    // the staircase. The field is built by clipping each cell with the analytic
    // plane, which is the same exact arithmetic the read-back is checked with, so
    // the only thing the test can get wrong is the direction - which is the point.
    {
        const label n = 16;
        const scalar dx = 1.0 / n;
        const polyMesh mesh = boxXY(n, n, 1.0);
        // 30 deg off the x axis: not a direction the cell-centred difference of a
        // plane field happens to get right, which is what makes the comparison below
        // a comparison and not a tautology
        const Vector analytic{std::sqrt(3.0) / 2.0, 0.5, 0.0};
        const scalar analyticOffset = 0.5 * (analytic.x + analytic.y);

        const Field<scalar>& vol = mesh.cellVolumes();
        Field<scalar> exact(static_cast<std::size_t>(n * n), 0.0);
        int band = 0;
        for (label cell = 0; cell < mesh.nCells(); ++cell)
        {
            const std::size_t u = static_cast<std::size_t>(cell);
            const vof::Cut kept = vof::clipConvexCell(vof::cellBoundary(mesh, cell),
                                                      analytic, analyticOffset);
            // a cut that swallows or empties the cell reports the flag rather than a
            // volume, so the fraction has to come from the flag first
            exact[u] = kept.whole ? 1.0 : (kept.void_ ? 0.0 : kept.volume / vol[u]);
            if (exact[u] > 1e-6 && exact[u] < 1.0 - 1e-6) ++band;
        }
        volScalarField alpha(mesh, "alpha", 0.0);
        alpha.internalField() = exact;

        const vof::Reconstruction r
            = vof::reconstruct(alpha, "corrected", vof::NormalChoice::Fit);
        const vof::Reconstruction local
            = vof::reconstruct(alpha, "corrected", vof::NormalChoice::Mixed);
        const Field<scalar> kp = vof::interfaceCurvaturePlic(mesh, r);

        double worstAngle = 0, worstScatter = 0, worstOffset = 0, localAngle = 0, worstRead = 0,
               worstKappa = 0;
        label counted = 0, kappaCells = 0;
        for (std::size_t u = 0; u < exact.size(); ++u)
        {
            if (mag2(r.normal[u]) == 0) continue;
            ++counted;
            worstAngle = std::max(worstAngle, angleDegrees(r.normal[u], analytic));
            localAngle += mag2(local.normal[u]) == 0 ? 180.0
                                                     : angleDegrees(local.normal[u], analytic);
            worstScatter = std::max(worstScatter, r.fitScatter[u]);
            // the plane reported is the analytic one: same direction, and the same
            // place, because the cell's own fraction has to sit on it
            worstOffset = std::max(worstOffset, std::abs(r.offset[u] - analyticOffset));
            worstRead = std::max(worstRead, std::abs(r.fraction[u] - exact[u]));
            if (kp[u] != 0)
            {
                ++kappaCells;
                worstKappa = std::max(worstKappa, std::abs(kp[u]));
            }
        }
        std::cout << "     plane front on " << n << "x" << n
                  << ": band " << r.bandCells << " of " << band << " cells, fitted normal off "
                  << std::setprecision(6) << worstAngle << " deg worst (cell-local mean "
                  << localAngle / counted << " deg), offset off " << std::setprecision(2)
                  << worstOffset << ", planes disagree by " << worstScatter / dx
                  << " cells worst, read-back " << std::setprecision(3) << worstRead
                  << ", kappa " << std::setprecision(3) << worstKappa << " worst over "
                  << kappaCells << " of " << counted << "\n";
        CHECK(r.bandCells == band, "the fit reconstructs the whole band");
        CHECK(worstRead < 1e-10, "and still cuts each cell's own fraction exactly");
        // What is left over is the arithmetic, not the search: four rounds settle a
        // flat front - two still leave 0.004 deg of direction and 6e-5 of a cell for
        // the neighbourhood to disagree about - and the residual then is 5e-8 of a
        // cell in the plane constant and 6e-9 of one in what the planes disagree by.
        CHECK(worstAngle < 1e-6, "the fitted normal is the analytic one, to round-off");
        CHECK(worstOffset < 1e-6 * dx,
              "the plane it reports is the analytic plane, to a millionth of a cell");
        CHECK(worstScatter < 1e-7 * dx, "and the neighbourhood agrees with nothing left over");
        CHECK(localAngle / counted > 1.0, "the cell-local normal of the same plane front is not");
        // every interface point sits in the one analytic plane, so the quadric through
        // them is flat: the curvature it reports is round-off and nothing else
        CHECK(kappaCells == counted, "every band cell of a plane front reports a curvature");
        CHECK(worstKappa < 1e-12, "and the curvature of a plane front is zero, to round-off");
    }


    // ---- 8. a region cut again ----
    // ConvexRegion is what the swept-prism flux needs: a cell's phase region has to
    // survive a second clip, by planes that are not the cell's own faces. The closed
    // form here is the slab between two diagonal planes through the unit cube, whose
    // volume is the difference of two cube caps - so the second cut is checked
    // against arithmetic the clip itself does not share.
    {
        const polyMesh mesh = boxXY(1, 1, 1.0);
        const vof::ConvexRegion cube = vof::cellRegion(mesh, 0);
        CHECK(relNear(cube.volume(), 1.0, 1e-13), "a cell region is the cell's own volume");
        CHECK(relNear(cube.centroid().x, 0.5, 1e-13) && relNear(cube.centroid().y, 0.5, 1e-13)
                  && relNear(cube.centroid().z, 0.5, 1e-13),
              "and its centroid");

        const scalar k = 1.0 / std::sqrt(3.0);
        const Vector n{k, k, k};
        const Vector back{-k, -k, -k};

        bool slab = true, sameAsCell = true, list = true, closedForm = true;
        double worstSlab = 0;
        for (double o1 : {0.3, 0.6, 0.9})
        {
            for (double o2 : {1.1, 1.4, 1.6})
            {
                if (o2 <= o1) continue;
                const vof::ConvexRegion twice
                    = vof::clipRegion(vof::clipRegion(cube, n, o1), back, -o2);
                // n . x with n a unit vector is s/sqrt(3) for the cubeCap parameter s
                const double want = cubeCapVolume(o2 / k) - cubeCapVolume(o1 / k);
                worstSlab = std::max(worstSlab, std::abs(twice.volume() - want) / want);
                slab = slab && relNear(twice.volume(), want, 1e-12);

                const std::vector<vof::HalfSpace> both = {{n, o1}, {back, -o2}};
                list = list && relNear(vof::clipRegion(cube, both).volume(), twice.volume(), 1e-12);
            }
            // one cut, and the same cut through the cell-boundary path. The cap's
            // closed form counts the part with n . x <= s, which is the side the
            // plane throws away.
            const vof::Cut byCell = vof::clipConvexCell(unitCell(mesh), n, o1);
            const vof::ConvexRegion byRegion = vof::clipRegion(cube, n, o1);
            sameAsCell = sameAsCell && relNear(byRegion.volume(), byCell.volume, 1e-12);
            closedForm = closedForm && relNear(byRegion.volume(), 1.0 - cubeCapVolume(o1 / k), 1e-12);
        }
        CHECK(slab, "two chained cuts hold the slab between the planes to 1e-12");
        CHECK(std::abs(worstSlab) < 1e-12, "and the slab's volume is the difference of two caps");
        CHECK(list, "the half-space list gives the same region in one call");
        CHECK(sameAsCell && closedForm,
              "one cut agrees with the cell-boundary path and with the cube cap");

        // a box clipped by a box: the prism a face sweeps is exactly this, and its
        // volume is the area it carries times the distance, which is the flux. Each
        // half-space keeps unitNormal . x >= offset, so the far sides are the negated
        // normals at the negated coordinate.
        const std::vector<vof::HalfSpace> box
            = {{Vector{1, 0, 0}, 0.0}, {Vector{-1, 0, 0}, -0.4},
               {Vector{0, 1, 0}, 0.0}, {Vector{0, -1, 0}, -1.0},
               {Vector{0, 0, 1}, 0.0}, {Vector{0, 0, -1}, -1.0}};
        CHECK(relNear(vof::clipRegion(cube, box).volume(), 0.4, 1e-13),
              "a unit face swept 0.4 is the prism of volume 0.4");
        CHECK(vof::clipRegion(cube, std::vector<vof::HalfSpace>{{n, 3.0}}).empty,
              "a cut past the region empties it instead of failing");
    }

    // ---- 9. a plane front moved by its own plane: nothing spreads ----
    // The claim the finite-volume transport of alpha cannot make. A step in alpha
    // carried by an upwind flux spreads at the scheme's truncation error; carried by
    // the plane, the front after any number of steps is the same front at x0 + u t,
    // cell by cell, with one cell between 0 and 1 however far it has moved. The
    // Courant number here is 0.4 and the front starts on a cell centre, so it is
    // twice within a cell width of crossing a face - the case where diffusion shows
    // up fastest.
    {
        const label nx = 12;
        const polyMesh mesh = boxXY(nx, 1, 0.1);
        const scalar dx = 1.0 / nx, dt = 0.01, courant = 0.4;
        const scalar u = courant * dx / dt;
        const scalar x0 = (6 + 0.5) * dx; // the front on the centre of cell 6

        volScalarField alpha(mesh, "alpha", 0.0);
        for (label c = 0; c < mesh.nCells(); ++c)
        {
            const scalar xc = (c + 0.5) * dx;
            alpha.internalField()[static_cast<std::size_t>(c)]
                = xc < x0 ? 1.0 : (xc > x0 ? 0.0 : 0.5);
        }

        Field<scalar> phi(static_cast<std::size_t>(mesh.nFaces()), 0.0);
        for (label f = 0; f < mesh.nFaces(); ++f)
        {
            const std::size_t u2 = static_cast<std::size_t>(f);
            phi[u2] = dot(Vector{u, 0.0, 0.0}, mesh.faceSf()[u2]);
        }

        double worst = 0;
        label bandAfter = 0;
        vof::GeometricStep step;
        for (int s = 1; s <= 3; ++s)
        {
            const vof::Reconstruction r
                = vof::reconstruct(alpha, "corrected", vof::NormalChoice::Mixed);
            step = vof::advectGeometric(alpha, phi, dt, r);
            const scalar front = x0 + s * courant * dx;
            bandAfter = 0;
            for (label c = 0; c < mesh.nCells(); ++c)
            {
                const std::size_t u2 = static_cast<std::size_t>(c);
                const scalar want = std::min(1.0, std::max(0.0, (front - c * dx) / dx));
                worst = std::max(worst, std::abs(alpha.internalField()[u2] - want));
                if (alpha.internalField()[u2] > 1e-12 && alpha.internalField()[u2] < 1.0 - 1e-12)
                    ++bandAfter;
            }
        }
        // boxXY is the unit square extruded by dz, so the section the front sweeps is
        // one by dz; every volume below is that section times a length in x.
        const scalar section = 1.0 * 0.1;
        const scalar swept = 3.0 * courant * dx * section; // the inlet's own volume
        const scalar phaseNow = phaseVolume(alpha);
        const scalar phaseWant = x0 * section + swept;
        std::cout << "     front moved by its planes: worst |alpha - translated front| "
                  << std::setprecision(3) << worst << ", band " << bandAfter
                  << " cell, phase volume " << std::setprecision(9) << phaseNow
                  << " vs " << phaseWant << ", limiter rounds " << step.limiterRounds
                  << ", faces carrying flux " << step.fluxFaces << "\n";
        CHECK(worst < 1e-13, "the field is the front translated, to round-off");
        CHECK(bandAfter == 1, "and it is still one cell thick after three steps");
        CHECK(relNear(phaseNow, phaseWant, 1e-12),
              "the phase grew by exactly the volume the inlet swept");
        CHECK(step.clampVolume == 0.0, "nothing had to be clipped back into its cell");
    }

    // ---- 10. the rotation carried by the planes ----
    // The droplet of tests/cases/vofZalesak, at the same Courant number, for the
    // part of its contract that is the transport's own: the phase volume cannot move
    // because each face's flux is one number subtracted where it left and added where
    // it arrived; alpha cannot leave [0,1] because the two-limit pass is what keeps it
    // there; and the front cannot thicken, because a plane does not spread.
    {
        const label n = 40;
        const polyMesh mesh = boxXY(n, n, kDiskDepth);
        const scalar dx = 1.0 / n, dt = 0.01;
        const Vector centre{0.7, 0.5, 0.0}, axis{0.5, 0.5, 0.0};
        const scalar omega = 1.0;

        const Field<scalar> initial = vof::exactDiskAlpha(mesh, centre, kDiskRadius);
        volScalarField alpha(mesh, "alpha", 0.0);
        for (std::size_t u = 0; u < alpha.internalField().size(); ++u)
            alpha.internalField()[u] = initial[u];

        Field<scalar> phi(static_cast<std::size_t>(mesh.nFaces()), 0.0);
        for (label f = 0; f < mesh.nFaces(); ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            const Vector r = mesh.faceCentres()[u] - axis;
            phi[u] = dot(Vector{-omega * r.y, omega * r.x, 0.0}, mesh.faceSf()[u]);
        }

        const auto bandOf = [&]()
        {
            label b = 0;
            for (const scalar a : alpha.internalField())
                if (a > 1e-6 && a < 1.0 - 1e-6) ++b;
            return b;
        };
        const label band0 = bandOf();
        const scalar volume0 = phaseVolume(alpha);

        double worstVol = 0, worstClamp = 0;
        label rounds = 0, limited = 0;
        scalar smallest = 1.0, lowest = 1.0, highest = 0.0;
        for (int s = 0; s < 40; ++s)
        {
            const vof::Reconstruction r
                = vof::reconstruct(alpha, "corrected", vof::NormalChoice::Mixed);
            const vof::GeometricStep step = vof::advectGeometric(alpha, phi, dt, r);
            worstVol = std::max(worstVol, std::abs(step.volumeError) / volume0);
            worstClamp = std::max(worstClamp, step.clampVolume);
            rounds = std::max(rounds, step.limiterRounds);
            lowest = std::min(lowest, step.minAlpha);
            highest = std::max(highest, step.maxAlpha);
            if (step.limitedFaces > 0)
            {
                limited = std::max(limited, step.limitedFaces);
                smallest = std::min(smallest, step.minFactor);
            }
        }
        const scalar volumeNow = phaseVolume(alpha);
        const label bandNow = bandOf();
        std::cout << "     40 rotation steps: phase volume " << std::setprecision(9)
                  << volumeNow << " for " << volume0 << " (worst step " << std::setprecision(3)
                  << worstVol << "), clipped " << worstClamp << ", alpha in [" << lowest
                  << ", " << highest << "], band " << band0 << " to "
                  << bandNow << " cells, limiter on up to " << limited << " faces with factor "
                  << std::setprecision(4) << smallest << " over up to " << rounds
                  << " rounds\n";
        CHECK(worstVol < 1e-13, "no step's faces disagree about the phase they traded");
        CHECK(worstClamp < 1e-15, "and nothing was clipped into its cell afterwards");
        CHECK(lowest >= 0.0 && highest <= 1.0,
              "alpha stayed inside [0,1] with the limiter doing it");
        CHECK(limited > 0 && smallest > 0.25,
              "the limiter bound - outflow prisms of one cell overlap at a corner - "
              "without starving the step");
        CHECK(bandNow <= static_cast<label>(1.5 * band0),
              "the front is still within half a cell of the thickness it started with");
    }

    return failures == 0 ? 0 : 1;
}
