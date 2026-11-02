// Unit tests for the kernel foundations: mesh geometry, dictionary,
// vector algebra and a minimal fvMatrix assembly.
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Vector.h"
#include "cf/fvm/BlockMesh.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"

using namespace cf;

static int failures = 0;

#define CHECK(cond, msg)                                                    \
    do {                                                                    \
        if (!(cond)) { ++failures; std::cout << "FAIL: " << (msg) << "\n"; }\
        else          { std::cout << "ok:   " << (msg) << "\n"; }            \
    } while (0)

static bool near(double a, double b, double tol = 1e-10)
{
    return std::abs(a - b) < tol;
}

int main()
{
    // ---- vector algebra ----
    Vector a{1, 2, 3}, b{4, -5, 6};
    CHECK(near(dot(a, b), 4 - 10 + 18), "dot product");
    Vector c = cross(a, b);
    CHECK(near(c.x, 2 * 6 - 3 * (-5)) && near(c.y, 3 * 4 - 1 * 6)
          && near(c.z, 1 * (-5) - 2 * 4), "cross product");

    // ---- dictionary ----
    Dictionary d = Dictionary::parseString(R"(
        alpha 1.5;
        name  myWall;
        dims  (2 3 4);
        words (inlet outlet);
        nested
        {
            inner { deep 7; }
        }
    )");
    CHECK(near(d.getScalar("alpha"), 1.5), "dict scalar");
    CHECK(d.getWord("name") == "myWall", "dict word");
    CHECK(d.getScalarList("dims").size() == 3, "dict numeric list");
    CHECK(d.getWordList("words")[1] == "outlet", "dict word list");
    CHECK(near(d.getDict("nested").getDict("inner").getScalar("deep"), 7), "dict nested");
    bool threw = false;
    try { d.getWord("missing"); } catch (const Error&) { threw = true; }
    CHECK(threw, "dict missing key throws");

    // ---- block mesh geometry: 2x2x2 unit cube ----
    polyMesh mesh = BlockMesh::create({0, 0, 0}, {1, 1, 1}, 2, 2, 2);
    CHECK(mesh.nCells() == 8, "mesh cell count");
    CHECK(mesh.nInternalFaces() == 12, "mesh internal faces");
    CHECK(mesh.nBoundaryFaces() == 24, "mesh boundary faces");
    double totVol = 0;
    for (const auto v : mesh.cellVolumes()) totVol += v;
    CHECK(near(totVol, 1.0, 1e-10), "mesh total volume");
    CHECK(near(mesh.cellVolumes()[0], 0.125, 1e-10), "mesh cell volume");

    // closed surface per cell: sum of area vectors = 0
    bool closed = true;
    {
        std::vector<Vector> sums(static_cast<std::size_t>(mesh.nCells()), Vector::zero());
        for (label f = 0; f < mesh.nFaces(); ++f)
        {
            const auto& sf = mesh.faceSf()[static_cast<std::size_t>(f)];
            sums[static_cast<std::size_t>(mesh.owner()[static_cast<std::size_t>(f)])] += sf;
            if (f < mesh.nInternalFaces())
            {
                sums[static_cast<std::size_t>(mesh.neighbour()[static_cast<std::size_t>(f)])] -= sf;
            }
        }
        for (const auto& s : sums) closed = closed && mag(s) < 1e-12;
    }
    CHECK(closed, "area vectors close every cell");

    // ---- fields + boundary linearisation ----
    volScalarField T(mesh, "T", 0.0);
    const label pLeft = mesh.patchIdByName("xmin");
    T.setPatchValue(pLeft, "fixedValue", 2.0);
    T.setPatchValue(mesh.patchIdByName("xmax"), "zeroGradient", 0.0);
    // find a boundary face on xmin
    label bf = -1;
    for (label f = mesh.nInternalFaces(); f < mesh.nFaces(); ++f)
    {
        const auto& pa = mesh.patches()[static_cast<std::size_t>(pLeft)];
        if (f >= pa.start && f < pa.start + pa.size) { bf = f; break; }
    }
    CHECK(bf >= 0, "found xmin boundary face");
    scalar coeff = 0, konst = 0;
    T.linearise(bf, coeff, konst);
    CHECK(near(coeff, 0) && near(konst, 2.0), "fixedValue linearisation");

    // ---- single-cell diffusion: solve laplacian == 0 with two Dirichlet ----
    polyMesh slab = BlockMesh::create({0, 0, 0}, {1, 1, 1}, 1, 1, 1);
    volScalarField F(slab, "F", 0.5);
    F.setPatchValue(slab.patchIdByName("xmin"), "fixedValue", 0.0);
    F.setPatchValue(slab.patchIdByName("xmax"), "fixedValue", 1.0);
    fvMatrix eq(slab, viewOf(F), "FEqn");
    fvm::laplacian(eq, 1.0, "none", -1);
    auto rep = eq.solve(fvMatrix::Controls{1e-10, 100, 0});
    // one cell, two Dirichlet faces: diag = sum 1/d * ( -1 ) contributions...
    // A = 2, rhs = 1  => F = 0.5 (linear profile through the cell centre)
    CHECK(rep.converged, "single-cell solve converged");
    CHECK(near(eq.solution()[0][0], 0.5, 1e-9), "single-cell linear profile");

    if (failures)
    {
        std::cout << failures << " test(s) FAILED\n";
        return 1;
    }
    std::cout << "all kernel unit tests passed\n";
    return 0;
}
