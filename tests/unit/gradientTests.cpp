// Unit tests for the default corrected Gauss gradient and the opt-in
// leastSquares scheme. Production internal faces use bounded linear
// interpolation and boundary faces use their prescribed value; the skewed
// tetrahedron cube checks that corrected reduces the plain-Gauss error without
// claiming mesh-independent exactness, while leastSquares remains the exact
// linear-field alternative.
//
// Case fields are the ones tests/cases/parallelDecomposition/caseTet.dict
// carries: a linear anchor T = 1 + 2x - y + 0.5z, grad = (2, -1, 0.5), and a
// quadratic T = 1 + 2x - y + 0.5z + 3x^2 - 0.25y^2 + 1.5xy - 2yz + 0.75zx.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
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
    const Field<scalar> ow = makeInternalOwnerWeights(mesh);
    const Field<Vector> g = fvc::grad(mesh, t,
        [&fc, &c, &mesh, &ow](label f, const Field<scalar>& data)
        {
            if (f < mesh.nInternalFaces())
            {
                const std::size_t u = static_cast<std::size_t>(f);
                return ow[u] * data[static_cast<std::size_t>(mesh.owner()[u])]
                     + (1.0 - ow[u]) * data[static_cast<std::size_t>(mesh.neighbour()[u])];
            }
            return poly(c, fc[static_cast<std::size_t>(f)]);
        },
        scheme);
    double worst = 0;
    for (label i = 0; i < mesh.nCells(); ++i)
    {
        const Vector e = g[static_cast<std::size_t>(i)] - polyGrad(c, cc[static_cast<std::size_t>(i)]);
        if (!std::isfinite(mag(e))) return std::numeric_limits<double>::infinity();
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

        const Field<scalar> ow = makeInternalOwnerWeights(mesh);
        double weightError = 0, weightBias = 0;
        for (label f = 0; f < mesh.nInternalFaces(); ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            const auto& co = mesh.cellCentres()[static_cast<std::size_t>(mesh.owner()[u])];
            const Vector d = mesh.deltaCoeffVec()[u];
            const double expected = 1.0 - std::clamp(
                dot(mesh.faceCentres()[u] - co, d) / dot(d, d), 0.0, 1.0);
            weightError = std::max(weightError, std::abs(ow[u] - expected));
            weightBias = std::max(weightBias, std::abs(ow[u] - 0.5));
        }
        CHECK(weightError < 1e-15, "tet cube: owner weights follow the owner-to-neighbour direction");
        CHECK(weightBias > 0.05, "tet cube: weight test exercises non-symmetric faces");

        const double linPlain = gradError(mesh, lin, "linear");
        const double linCorrected = gradError(mesh, lin, "corrected");
        const double linLS = gradError(mesh, lin, "leastSquares");
        std::cout << "info: tet linear grad error plain=" << linPlain
                  << " corrected=" << linCorrected << " leastSquares=" << linLS << "\n";
        CHECK(linCorrected < 0.5 * linPlain && linCorrected < 2.0,
              "tet cube: corrected at least halves the plain-Gauss linear error");
        CHECK(linLS < 1e-9, "tet cube: leastSquares recovers the linear gradient");

        const double quadPlain = gradError(mesh, quad, "linear");
        const double quadLS = gradError(mesh, quad, "leastSquares");
        const double quadCorrected = gradError(mesh, quad, "corrected");
        std::cout << "info: tet quadratic grad error plain=" << quadPlain
                  << " corrected=" << quadCorrected << " leastSquares=" << quadLS << "\n";
        CHECK(quadCorrected < 0.5 * quadPlain,
              "tet cube: corrected at least halves the plain-Gauss quadratic error");
        CHECK(quadLS < quadCorrected,
              "tet cube: leastSquares improves on corrected for a quadratic field");
    }

    // ---- orthogonal hex cube: the degeneration to Green-Gauss ----
    {
        const polyMesh hex = GmshMesh::fromFile(std::string(CF_MESH_DIR) + "/cubeHex_41.msh");
        CHECK(gradError(hex, lin, "linear") < 1e-9,
              "hex cube: plain Green-Gauss is exact for a linear field");
        CHECK(gradError(hex, lin, "corrected") < 1e-9,
              "hex cube: corrected degeneration is exact for a linear field");
        CHECK(gradError(hex, lin, "leastSquares") < 1e-9,
              "hex cube: leastSquares is exact for a linear field");
        const double quadPlain = gradError(hex, quad, "linear");
        CHECK(std::abs(gradError(hex, quad, "corrected") - quadPlain) < 1e-12,
              "hex cube: corrected quadratic result degenerates to plain Green-Gauss");
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
