// Unit tests for the Gmsh .msh reader and the case-file mesh selection:
// topology, exact polyhedron geometry, patch naming, cross-format equality
// and blockMesh equivalence.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cf/common/Error.h"
#include "cf/fvm/BlockMesh.h"
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

static bool near(double a, double b, double tol = 1e-10)
{
    return std::abs(a - b) < tol;
}

static std::string meshPath(const char* name)
{
    return std::string(CF_MESH_DIR) + "/" + name;
}

static const polyMesh::Patch* patchOf(const polyMesh& mesh, const std::string& name)
{
    for (const auto& p : mesh.patches())
    {
        if (p.name == name) return &p;
    }
    return nullptr;
}

static double totalVolume(const polyMesh& mesh)
{
    double v = 0;
    for (const scalar c : mesh.cellVolumes()) v += c;
    return v;
}

/// sum over a cell's faces of the oriented area vector; exactly zero for a
/// closed cell whose faces are consistently oriented (owner -> neighbour /
/// outward at the boundary).
static double maxCellClosureError(const polyMesh& mesh)
{
    std::vector<Vector> acc(static_cast<std::size_t>(mesh.nCells()), Vector::zero());
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        acc[static_cast<std::size_t>(mesh.owner()[u])] += mesh.faceSf()[u];
        if (f < mesh.nInternalFaces())
        {
            acc[static_cast<std::size_t>(mesh.neighbour()[u])] -= mesh.faceSf()[u];
        }
    }
    double worst = 0;
    for (const Vector& a : acc) worst = std::max(worst, static_cast<double>(mag(a)));
    return worst;
}

/// Sorted geometry fingerprint: cell volumes and centre coordinates plus the
/// face areas. Independent of the element order in the file, so it can
/// compare a 2.2 file with the 4.1 export of the same mesh.
static std::vector<double> fingerprint(const polyMesh& mesh)
{
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(mesh.nCells()) * 4 + mesh.nFaces());
    for (const scalar v : mesh.cellVolumes()) out.push_back(v);
    for (const Vector& c : mesh.cellCentres())
    {
        out.push_back(c.x);
        out.push_back(c.y);
        out.push_back(c.z);
    }
    for (const scalar a : mesh.faceAreas()) out.push_back(a);
    std::sort(out.begin(), out.end());
    return out;
}

static bool sameFingerprint(const polyMesh& a, const polyMesh& b)
{
    const std::vector<double> fa = fingerprint(a), fb = fingerprint(b);
    if (fa.size() != fb.size()) return false;
    for (std::size_t i = 0; i < fa.size(); ++i)
    {
        if (!near(fa[i], fb[i], 1e-9)) return false;
    }
    return true;
}

/// Gauss gradient of the linear field T = z, whose exact value is (0, 0, 1).
static double gradZError(const polyMesh& mesh, const std::string& scheme)
{
    const auto& cc = mesh.cellCentres();
    const auto& fc = mesh.faceCentres();
    Field<scalar> Tz(static_cast<std::size_t>(mesh.nCells()));
    for (label c = 0; c < mesh.nCells(); ++c) Tz[static_cast<std::size_t>(c)] = cc[static_cast<std::size_t>(c)].z;
    const Field<Vector> g = fvc::grad(mesh, Tz,
                                      [&](label f, const Field<scalar>&)
                                      { return fc[static_cast<std::size_t>(f)].z; },
                                      scheme);
    double err = 0;
    for (label c = 0; c < mesh.nCells(); ++c)
    {
        err = std::max(err, static_cast<double>(mag(g[static_cast<std::size_t>(c)] - Vector{0, 0, 1})));
    }
    return err;
}

int main()
{
    // ---- a single hex with no physical groups: bounding-box plane naming ----
    {
        polyMesh mesh = GmshMesh::fromFile(meshPath("singleHex_22.msh"));
        CHECK(mesh.nCells() == 1 && mesh.nPoints() == 8, "single hex: 1 cell, 8 nodes");
        CHECK(mesh.nFaces() == 6 && mesh.nInternalFaces() == 0, "single hex: 6 boundary faces");
        CHECK(mesh.nPatches() == 6, "single hex: 6 plane patches");
        bool sized = true;
        for (const char* n : {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"})
        {
            const polyMesh::Patch* p = patchOf(mesh, n);
            if (!p || p->size != 1 || p->type != "wall") sized = false;
        }
        CHECK(sized, "single hex: every plane is one wall patch");
        CHECK(near(totalVolume(mesh), 1.0) && near(mesh.cellVolumes()[0], 1.0),
              "single hex: unit volume");
        CHECK(near(mesh.cellCentres()[0].x, 0.5) && near(mesh.cellCentres()[0].y, 0.5)
              && near(mesh.cellCentres()[0].z, 0.5), "single hex: centroid");

        const std::map<std::string, std::string> types{{"xmin", "inlet"}, {"xmax", "outlet"}};
        polyMesh marked = GmshMesh::fromFile(meshPath("singleHex_22.msh"), types);
        CHECK(patchOf(marked, "xmin")->type == "inlet"
              && patchOf(marked, "xmax")->type == "outlet"
              && patchOf(marked, "ymin")->type == "wall",
              "single hex: case patch types override, others stay wall");

        // exact polyhedron centroid rather than an area-weighted estimate
        CHECK(near(maxCellClosureError(marked), 1e-14), "single hex: closed faces");
    }

    // ---- hex cube: equivalence with blockMesh ----
    {
        polyMesh gmsh = GmshMesh::fromFile(meshPath("cubeHex_41.msh"));
        polyMesh block = BlockMesh::create({0, 0, 0}, {1, 1, 1}, 2, 2, 2);
        CHECK(gmsh.nCells() == 8 && gmsh.nPoints() == 27, "hex cube: 8 cells, 27 nodes");
        CHECK(gmsh.nFaces() == block.nFaces() && gmsh.nInternalFaces() == 12,
              "hex cube: 36 faces, 12 internal - same as blockMesh");
        CHECK(gmsh.nPatches() == 6 && patchOf(gmsh, "hex_s401") != nullptr,
              "hex cube: 6 named physical groups");
        CHECK(sameFingerprint(gmsh, block), "hex cube: identical geometry to blockMesh");
        CHECK(near(totalVolume(gmsh), 1.0) && near(gmsh.cellVolumes()[0], 0.125),
              "hex cube: 8 cells of 1/8");
        CHECK(near(maxCellClosureError(gmsh), 1e-15), "hex cube: closed faces");
        CHECK(gradZError(gmsh, "linear") < 1e-13 && gradZError(gmsh, "corrected") < 1e-13,
              "hex cube: grad(z) exact with both schemes");
    }

    // ---- tet cube: the same mesh exported in 2.2 and 4.1 ----
    {
        polyMesh v22 = GmshMesh::fromFile(meshPath("cubeTets_22.msh"));
        polyMesh v41 = GmshMesh::fromFile(meshPath("cubeTets_41.msh"));
        CHECK(v22.nCells() == 95 && v22.nPoints() == 44, "tets: 95 cells, 44 nodes");
        CHECK(v22.nFaces() == 232 && v22.nInternalFaces() == 148, "tets: 232 faces, 148 internal");
        CHECK(v41.nPatches() == 6 && v22.nPatches() == 6, "tets: 6 surface physical groups");
        bool samePatches = true;
        for (const auto& p : v22.patches())
        {
            const polyMesh::Patch* q = patchOf(v41, p.name);
            if (!q || q->size != p.size) samePatches = false;
        }
        CHECK(samePatches, "tets: same patch names and sizes in 2.2 and 4.1");
        CHECK(sameFingerprint(v22, v41), "tets: 2.2 and 4.1 exports give the same mesh");
        CHECK(near(totalVolume(v22), 1.0), "tets: unit total volume");
        CHECK(near(maxCellClosureError(v22), 1e-15), "tets: closed faces");
        // the raw Gauss gradient is exact for a linear field on any mesh; the
        // non-orthogonal correction is a refinement approximation
        CHECK(gradZError(v22, "linear") < 1e-14, "tets: linear grad(z) exact");
    }

    // ---- prism channel: genuinely non-orthogonal, exact by symmetry ----
    {
        polyMesh mesh = GmshMesh::fromFile(meshPath("prismChannel_22.msh"));
        CHECK(mesh.nCells() == 20 && mesh.nPoints() == 44, "prisms: 20 cells, 44 nodes");
        CHECK(mesh.nInternalFaces() == 28 && mesh.nBoundaryFaces() == 44,
              "prisms: 28 internal, 44 boundary faces");
        bool named = true;
        const std::pair<const char*, long> expect[6] =
            {{"inlet", 2}, {"outlet", 2}, {"wall_y0", 10},
             {"wall_x1", 10}, {"wall_y1", 10}, {"wall_x0", 10}};
        for (const auto& e : expect)
        {
            const polyMesh::Patch* p = patchOf(mesh, e.first);
            if (!p || p->size != e.second) named = false;
        }
        CHECK(named, "prisms: boundary faces grouped by physical group name");
        CHECK(near(totalVolume(mesh), 0.1), "prisms: volume 1 x 0.1 x 1");

        // exact centroid of the prism (0,0,0)-(1,0,0)-(1,0.1,0) extruded 0.1
        const Vector& c0 = mesh.cellCentres()[0];
        CHECK(near(c0.x, 2.0 / 3.0) && near(c0.y, 1.0 / 30.0) && near(c0.z, 0.05),
              "prisms: cell centre is the polyhedron centroid");
        CHECK(near(maxCellClosureError(mesh), 1e-16), "prisms: closed faces");

        // 78 deg of measured non-orthogonality, yet the correction vanishes
        // for a field varying only along the sweep direction
        CHECK(gradZError(mesh, "linear") < 1e-14, "prisms: linear grad(z) exact");
        CHECK(gradZError(mesh, "corrected") < 1e-11, "prisms: corrected grad(z) exact");
    }

    // ---- case-file selection ----
    {
        const Dictionary fromBlock = Dictionary::parseString(R"(
            block { min (0 0 0); max (1 1 1); cells (2 2 2); }
        )");
        CHECK(MeshLoader::fromCase(fromBlock, "").nCells() == 8, "loader: block {} still works");

        const std::string meshCase =
            "mesh {\n"
            "    format gmsh;\n"
            "    file \"" + meshPath("prismChannel_22.msh") + "\";\n"
            "    patches { inlet { type inlet; } outlet { type outlet; } }\n"
            "}\n";
        const Dictionary fromFile = Dictionary::parseString(meshCase);
        polyMesh mesh = MeshLoader::fromCase(fromFile, "");
        CHECK(mesh.nCells() == 20, "loader: mesh {} reads the file");
        CHECK(patchOf(mesh, "inlet")->type == "inlet"
              && patchOf(mesh, "outlet")->type == "outlet"
              && patchOf(mesh, "wall_y0")->type == "wall",
              "loader: patch types come from the case");

        bool threw = false;
        try
        {
            MeshLoader::fromCase(Dictionary::parseString("nu 1.0;"), "");
        }
        catch (const Error&) { threw = true; }
        CHECK(threw, "loader: a case with no mesh is an error");

        threw = false;
        try
        {
            MeshLoader::fromCase(Dictionary::parseString(
                "block { min (0 0 0); max (1 1 1); cells (1 1 1); }\n"
                "mesh { file \"x.msh\"; }\n"), "");
        }
        catch (const Error&) { threw = true; }
        CHECK(threw, "loader: block {} and mesh {} together is an error");

        threw = false;
        try
        {
            GmshMesh::fromFile(meshPath("noSuchMesh.msh"));
        }
        catch (const Error&) { threw = true; }
        CHECK(threw, "reader: a missing file is an error");
    }

    if (failures)
    {
        std::cout << failures << " gmsh mesh checks failed\n";
        return 1;
    }
    std::cout << "all gmsh mesh checks passed\n";
    return 0;
}
