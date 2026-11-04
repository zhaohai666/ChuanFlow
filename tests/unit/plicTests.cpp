// Unit tests for the PLIC layer. The cut geometry is analytic, so every test
// here is exactness rather than convergence: the volume, centroid and cut area
// of a plane through a cell against a closed form, the cell's own boundary loops
// read back as the mesh's cell volume on hexes, tets and prisms, the fraction a
// reconstruction reads back from the field it was given, and the normal and
// curvature a cell-filled disk gives up.
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
    double kappaMean = 0;    ///< band-mean of div(normal)
};

Measured measure(const polyMesh& mesh, const volScalarField& alpha,
                 const Field<scalar>& exact, const vof::Reconstruction& r)
{
    Measured m;
    const Vector c = diskCentre();
    const Field<scalar> k = vof::interfaceCurvature(mesh, r.normal);
    label kappaN = 0;
    double angle = 0, kappa = 0;

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

        // towards increasing alpha is into the disk, and the interface of a disk
        // is radial through its own sub-cell centre
        const Vector d = r.areaCentroid[u] - c;
        const double deg = angleDegrees(r.normal[u], d * (-1.0 / mag(d)));
        angle += deg;
        m.worstAngle = std::max(m.worstAngle, deg);
        kappa += k[u]; ++kappaN;
    }
    m.meanAngle = m.reconstructed > 0 ? angle / m.reconstructed : 0;
    m.kappaMean = kappaN > 0 ? kappa / kappaN : 0;

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
                 {"mixed", vof::NormalChoice::Mixed}})
        {
            const Measured m = measure(mesh, alpha, exact,
                                            vof::reconstruct(alpha, "corrected",
                                                             entry.second));
            std::cout << "     " << std::setw(8) << entry.first << ": band "
                      << m.reconstructed << " of " << m.band
                      << " (the gradient threshold sees " << m.gradientBand << ")"
                      << ", read-back error " << std::setprecision(3) << m.worstRead
                      << ", normal off-radial " << m.meanAngle << " deg mean /"
                      << m.worstAngle << " deg max\n     " << std::setw(8) << ""
                      << "interface area " << std::setprecision(6) << m.area
                      << " vs the extruded perimeter 2 pi R dz = "
                      << 2 * M_PI * R * kDiskDepth << ", mean kappa " << std::setprecision(5)
                      << m.kappaMean << " vs -1/R = " << -1.0 / R << "\n";
            CHECK(m.worstRead < 1e-10,
                  (entry.first + ": the reconstructed plane cuts the fraction it was given")
                      .c_str());
            CHECK(m.reconstructed == m.band,
                  (entry.first + ": every cell of the band is reconstructed").c_str());
            CHECK(m.area > 0.85 * 2 * M_PI * R * kDiskDepth && m.area < 2 * M_PI * R * kDiskDepth,
                  (entry.first + ": the band's interface area is the perimeter's").c_str());
        }
    }

    // ---- 7. what PLIC adds, measured on the mesh it is asked for ----
    // Two things, and the numbers say which is which.
    //
    // The band: it comes from alpha alone, so on a cell-filled interface it is
    // exactly the cells that hold one - all of them, and no others. The
    // gradient-threshold operator the smeared-band path uses cannot say that: on
    // the same field it marks the ring around the band, whose alpha is 0 or 1, as
    // carrying an interface, and the curvature stencil then reads faces into it.
    //
    // The geometry: the plane, its area and its centre are exact, so what is left
    // approximate is the DIRECTION the plane is given - the normal estimator, not
    // the clipping. A one-cell-thick alpha field is a staircase along the circle
    // and its slope barely converges with the mesh, which is why the angle below
    // hardly moves from 8 to 16 cells per radius, and why a smoothing or
    // area-based normal estimator is the next step rather than this one.
    {
        const scalar R = kDiskRadius;
        const scalar perimeter = static_cast<scalar>(2 * M_PI * R * kDiskDepth);
        Measured m40, m80;
        for (const label n : {40, 80})
        {
            const polyMesh mesh = boxXY(n, n, kDiskDepth);
            const Field<scalar> exact = vof::exactDiskAlpha(mesh, diskCentre(), R);
            volScalarField alpha(mesh, "alpha", 0.0);
            for (std::size_t u = 0; u < alpha.internalField().size(); ++u)
                alpha.internalField()[u] = exact[u];

            const Measured m = measure(mesh, alpha, exact,
                                       vof::reconstruct(alpha, "corrected",
                                                        vof::NormalChoice::Mixed));
            std::cout << "     " << n << "x" << n << ", R/dx = " << R * n
                      << ": alpha band " << m.band << ", PLIC " << m.reconstructed
                      << ", gradient threshold " << m.gradientBand << " of which "
                      << m.gradientOutside << " hold no interface\n     normal off-radial "
                      << std::setprecision(4) << m.meanAngle << " deg mean /"
                      << m.worstAngle << " deg max, read-back " << std::setprecision(3)
                      << m.worstRead << ", interface area " << std::setprecision(6) << m.area
                      << " vs " << perimeter << ", kappa " << std::setprecision(5)
                      << m.kappaMean << " vs " << -1.0 / R << "\n";
            (n == 40 ? m40 : m80) = m;
            CHECK(m.reconstructed == m.band,
                  "the reconstruction covers the whole band on a one-cell-thick front");
            CHECK(m.worstRead < 1e-10, "and reproduces the fraction everywhere");
        }
        CHECK(m40.gradientOutside > 0 && m80.gradientOutside > 0,
              "the gradient threshold marks cells that hold no interface; the alpha band does not");
        CHECK(m40.meanAngle < 20.0 && m80.meanAngle < 20.0 && m40.worstAngle < 25.0
                  && m80.worstAngle < 25.0,
              "the reconstructed normal stays within 20 deg mean of the analytic radial");
        CHECK(m40.area > 0.85 * perimeter && m80.area > 0.85 * perimeter
                  && m40.area <= perimeter && m80.area <= perimeter,
              "the band's interface area is the circle's perimeter to within 15%");
        CHECK(std::abs(m40.kappaMean * R + 1.0) < 0.4
                  && std::abs(m80.kappaMean * R + 1.0) < 0.4,
              "and its mean curvature is -1/R to within what a one-cell stencil can resolve");
    }

    return failures == 0 ? 0 : 1;
}
