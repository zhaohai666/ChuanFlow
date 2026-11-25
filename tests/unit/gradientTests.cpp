// Unit tests for the leastSquares gradient scheme: on the Gmsh tetrahedron
// cube whose skewed faces make the bounded linear weights clamp (and with them
// the corrected Gauss scheme miss even a linear field's gradient), the
// opt-in leastSquares scheme recovers the analytic gradients. The orthogonal
// degeneration -- leastSquares agrees with Green-Gauss where Green-Gauss is
// already exact -- is checked on a blockMesh hex cube.
//
// Case fields are the ones tests/cases/parallelDecomposition/caseTet.dict
// carries: a linear anchor T = 1 + 2x - y + 0.5z, grad = (2, -1, 0.5), and a
// quadratic T = 1 + 2x - y + 0.5z + 3x^2 - 0.25y^2 + 1.5xy - 2yz + 0.75zx.
#include <cmath>
#include <iostream>
#include <string>

#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/GmshMesh.h"
#include "cf/fvm/MeshLoader.h"

#ifndef CF_MESH_DIR
#define CF_MESH_DIR "."
#endif

using namespace cf;

static int failures = 0;

#define CHECK(cond, msg)                                                      \
    do {                                                                      \
        if (!(cond)) { ++failures; std::cout << "FAIL: " << (msg) << "\n"; }   \
        else          { std::cout << "ok:   " << (msg) << "\n"; }              \
    } while (0)

namespace
{

scalar poly(const std::vector<double>& c, const Vector& p)
{
    return c[0] + c[1] * p.x + c[2] * p.y + c[3] * p.z
         + c[4] * p.x * p.x + c[5] * p.y * p.y + c[6] * p.z * p.z
         + c[7] * p.x * p.y + c[8] * p.y * p.z + c[9] * p.z * p.x;
}

Vector polyGrad(const std::vector<double>& c, const Vector& p)
{
    return {c[1] + 2.0 * c[4] * p.x + c[7] * p.y + c[9] * p.z,
            c[2] + 2.0 * c[5] * p.y + c[7] * p.x + c[8] * p.z,
            c[3] + 2.0 * c[6] * p.z + c[8] * p.y + c[9] * p.x};
}

/// Largest cell-wise deviation of fvc::grad(scheme) from the analytic gradient;
/// internal faces use bounded linear interpolation, boundary faces the exact
/// prescribed value (the case's own contract).
double gradError(const polyMesh& mesh, const std::vector<double>& c, const std::string& scheme)
{
    const auto& cc = mesh.cellCentres();
    const auto& fc = mesh.faceCentres();
    Field<scalar> t(static_cast<std::size_t>(mesh.nCells()));
    for (label i = 0; i < mesh.nCells(); ++i) t[static_cast<std::size_t>(i)] = poly(c, cc[static_cast<std::size_t>(i)]);
    const Field<Vector> g = fvc::grad(mesh, t,
        [&fc, &c](label f, const Field<scalar>&) { return poly(c, fc[static_cast<std::size_t>(f)]); },
        scheme);
    double worst = 0;
    for (label i = 0; i < mesh.nCells(); ++i)
    {
        const Vector e = g[static_cast<std::size_t>(i)] - polyGrad(c, cc[static_cast<std::size_t>(i)]);
        worst = std::max(worst, static_cast<double>(mag(e)));
    }
    return worst;
}

} // namespace

int main()
{
    const std::vector<double> lin{1, 2, -1, 0.5, 0, 0, 0, 0, 0, 0};
    const std::vector<double> quad{1, 2, -1, 0.5, 3, -0.25, 0, 1.5, -2, 0.75};

    // ---- the skewed tet cube of caseTet.dict ----
    {
        const polyMesh mesh = GmshMesh::fromFile(std::string(CF_MESH_DIR) + "/cubeTets_22.msh");
        CHECK(mesh.nCells() == 95, "tet cube: 95 cells");

        const double linCorrected = gradError(mesh, lin, "corrected");
        const double linLS = gradError(mesh, lin, "leastSquares");
        std::cout << "info: tet linear grad error corrected=" << linCorrected
                  << " leastSquares=" << linLS << "\n";
        CHECK(linCorrected > 1.0, "tet cube: corrected misses the linear gradient (the pain point)");
        CHECK(linLS < 1e-9, "tet cube: leastSquares recovers the linear gradient");

        const double quadLS = gradError(mesh, quad, "leastSquares");
        const double quadCorrected = gradError(mesh, quad, "corrected");
        std::cout << "info: tet quadratic grad error corrected=" << quadCorrected
                  << " leastSquares=" << quadLS << "\n";
        CHECK(quadLS < 0.05 * quadCorrected,
              "tet cube: leastSquares cuts the quadratic gradient error five-fold");
    }

    // ---- orthogonal hex cube: the degeneration to Green-Gauss ----
    {
        const polyMesh hex = GmshMesh::fromFile(std::string(CF_MESH_DIR) + "/cubeHex_41.msh");
        CHECK(gradError(hex, lin, "leastSquares") < 1e-9,
              "hex cube: leastSquares is exact for a linear field");
        CHECK(gradError(hex, quad, "leastSquares") < 1e-9,
              "hex cube: leastSquares is exact for a quadratic field");
    }

    // ---- leastSquares is opt-in: the default path is untouched ----
    {
        const polyMesh mesh = GmshMesh::fromFile(std::string(CF_MESH_DIR) + "/cubeTets_22.msh");
        const double errLinear = gradError(mesh, lin, "linear");
        const double errUnknown = gradError(mesh, lin, "unheard-of");
        CHECK(errLinear == errUnknown && errUnknown > 0,
              "tet cube: an unknown grad scheme falls back to plain Green-Gauss");
    }

    if (failures)
    {
        std::cout << failures << " gradient check(s) FAILED\n";
        return 1;
    }
    std::cout << "all gradient unit tests passed\n";
    return 0;
}
