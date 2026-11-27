// chflow-3d-cavitation: homogeneous transport cavitation (F12) riding on the
// pressure-based PISO coupling.
//
// Opt-in case keys:
//   cavitation { pSat ...; rhoL 998; rhoV 1; Ce 2897; Cc 80; Rb 1e-5;
//                alpha0 1e-6; diffusivity 0; }
//   venturi { nx 60; ny 10; x1 0.10; x2 0.20; x3 0.25; x4 0.40;
//             H1 0.10; H2 0.05; depth 0.02; }
//           - a converging-diverging channel built in code instead of the
//             `block` mesh, so the Bernoulli acceleration and the diffuser
//             recovery happen where the checker can read them. Without a
//             venturi{} block the case uses its `block` mesh (the uniform
//             logistic-growth case).
//
// The mixture density alpha rhoV + (1-alpha) rhoL enters the pressure
// correction through PisoControls::cellInvRho (response scaled by
// rhoL/rhoMix); the momentum equation itself stays kinematic.
//
// Reports per `piso{ reportEvery }` steps and at the end:
//   t, pMin, pMax, alphaMax, and the column-mean velocity/pressure/void
//   fraction at the inlet, throat and outlet planes.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/CavitationModel.h"
#include "cf/fvm/MeshLoader.h"
#include "case3dTurb.h"

using namespace cf;

namespace
{

/// Structured converging-diverging channel: nx columns over [0, x4] with ny
/// cells across a height h(x) that is H1 up to x1, linear down to H2 at x2,
/// H2 through the throat to x3, and linear back to H1 at x4. Built the way
/// blockMesh registers faces: internal faces through a sorted-vertex map
/// (owner = the first, lower-indexed cell), boundary faces appended per patch.
polyMesh makeVenturi(int nx, int ny, scalar x1, scalar x2, scalar x3, scalar x4,
                    scalar H1, scalar H2, scalar depth)
{
    auto height = [&](scalar x) -> scalar
    {
        if (x <= x1) return H1;
        if (x < x2) return H1 + (H2 - H1) * (x - x1) / (x2 - x1);
        if (x <= x3) return H2;
        return H2 + (H1 - H2) * (x - x3) / (x4 - x3);
    };

    std::vector<Vector> pts;
    auto pid = [&](int i, int j, int k) -> label
    {
        return static_cast<label>(k * (ny + 1) * (nx + 1) + j * (nx + 1) + i);
    };
    pts.resize(static_cast<std::size_t>((nx + 1) * (ny + 1) * 2));
    for (int k = 0; k < 2; ++k)
    {
        for (int j = 0; j <= ny; ++j)
        {
            for (int i = 0; i <= nx; ++i)
            {
                const scalar x = x4 * i / nx;
                const scalar h = height(x);
                pts[static_cast<std::size_t>(pid(i, j, k))] =
                    Vector{x, -0.5 * h + h * j / ny, depth * k};
            }
        }
    }

    std::vector<std::vector<label>> internalFaces;
    std::vector<label> own, nbr;
    std::map<std::vector<label>, label> lookup;
    auto addInternal = [&](std::vector<label> v, label o, label n)
    {
        std::vector<label> key = v;
        std::sort(key.begin(), key.end());
        const auto it = lookup.find(key);
        if (it == lookup.end())
        {
            lookup[key] = static_cast<label>(internalFaces.size());
            internalFaces.push_back(v);
            own.push_back(o);
            nbr.push_back(n);
        }
        else
        {
            nbr[static_cast<std::size_t>(it->second)] = n;
        }
    };

    std::vector<std::vector<label>> inlet, outlet, bottom, top, front, back;
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            const label c = static_cast<label>(j * nx + i);
            const label v0 = pid(i, j, 0), v1 = pid(i + 1, j, 0),
                        v2 = pid(i + 1, j + 1, 0), v3 = pid(i, j + 1, 0),
                        v4 = pid(i, j, 1), v5 = pid(i + 1, j, 1),
                        v6 = pid(i + 1, j + 1, 1), v7 = pid(i, j + 1, 1);
            if (i == 0) inlet.push_back({v0, v3, v7, v4});
            if (i == nx - 1) outlet.push_back({v1, v2, v6, v5});
            if (j == 0) bottom.push_back({v0, v1, v5, v4});
            if (j == ny - 1) top.push_back({v2, v3, v7, v6});
            front.push_back({v0, v1, v2, v3});
            back.push_back({v4, v5, v6, v7});

            if (i < nx - 1) addInternal({v1, v2, v6, v5}, c, static_cast<label>((j)*nx + i + 1));
            if (j < ny - 1) addInternal({v2, v3, v7, v6}, c, static_cast<label>((j + 1) * nx + i));
        }
    }

    // Boundary owners are registered directly while sweeping the cells, so
    // each patch face points at the exact cell that emitted it.
    std::vector<std::vector<label>> faces = internalFaces;
    std::vector<label> owner = own, neighbour = nbr;
    std::vector<polyMesh::Patch> patches;
    label start = static_cast<label>(faces.size());
    struct BFace
    {
        std::vector<label> verts;
        label owner;
    };
    std::vector<std::pair<std::string, std::string>> order = {
        {"inlet", "inlet"}, {"outlet", "outlet"}, {"bottom", "wall"},
        {"top", "wall"}, {"front", "symmetry"}, {"back", "symmetry"}};
    std::map<std::string, std::vector<BFace>> boundary;
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            const label c = static_cast<label>(j * nx + i);
            const label v0 = pid(i, j, 0), v1 = pid(i + 1, j, 0),
                        v2 = pid(i + 1, j + 1, 0), v3 = pid(i, j + 1, 0),
                        v4 = pid(i, j, 1), v5 = pid(i + 1, j, 1),
                        v6 = pid(i + 1, j + 1, 1), v7 = pid(i, j + 1, 1);
            if (i == 0) boundary["inlet"].push_back({{v0, v3, v7, v4}, c});
            if (i == nx - 1) boundary["outlet"].push_back({{v1, v2, v6, v5}, c});
            if (j == 0) boundary["bottom"].push_back({{v0, v1, v5, v4}, c});
            if (j == ny - 1) boundary["top"].push_back({{v2, v3, v7, v6}, c});
            boundary["front"].push_back({{v0, v1, v2, v3}, c});
            boundary["back"].push_back({{v4, v5, v6, v7}, c});
        }
    }
    for (const auto& o : order)
    {
        const auto& vec = boundary[o.first];
        for (const auto& bf : vec)
        {
            faces.push_back(bf.verts);
            owner.push_back(bf.owner);
        }
        patches.push_back({o.first, o.second, start, static_cast<label>(vec.size())});
        start += static_cast<label>(vec.size());
    }

    Field<Vector> points(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) points[i] = pts[i];
    return polyMesh(std::move(points), faces,
                    Field<label>(owner.begin(), owner.end()),
                    Field<label>(neighbour.begin(), neighbour.end()), patches);
}

} // namespace

static int appRun(int argc, char** argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr, "usage: chflow-3d-cavitation <case.dict>\n");
        return 1;
    }
    const Dictionary caseDict = Dictionary::parseFile(argv[1]);

    int nx = 0, ny = 0;
    polyMesh mesh = caseDict.has("venturi")
                        ? makeVenturi(
                              static_cast<int>(caseDict.getDict("venturi").getScalar("nx")),
                              static_cast<int>(caseDict.getDict("venturi").getScalar("ny")),
                              caseDict.getDict("venturi").getScalar("x1"),
                              caseDict.getDict("venturi").getScalar("x2"),
                              caseDict.getDict("venturi").getScalar("x3"),
                              caseDict.getDict("venturi").getScalar("x4"),
                              caseDict.getDict("venturi").getScalar("H1"),
                              caseDict.getDict("venturi").getScalar("H2"),
                              caseDict.getDict("venturi").getScalar("depth"))
                        : MeshLoader::fromCase(caseDict, argv[1]);
    if (caseDict.has("venturi"))
    {
        nx = static_cast<int>(caseDict.getDict("venturi").getScalar("nx"));
        ny = static_cast<int>(caseDict.getDict("venturi").getScalar("ny"));
    }

    const scalar nu = app3dturb::caseNu(caseDict);
    volVectorField U(mesh, "U");
    volScalarField p(mesh, "p");
    surfaceScalarField phi(mesh, "phi", 0.0);
    app3dturb::applyPatches(caseDict, mesh, U, p);

    Schemes schemes = loadSchemes(caseDict);
    schemes.ddt = "Euler";

    phi.correctBoundaryFlux(U);
    const CavitationModel::Config cfg = CavitationModel::Config::fromCase(caseDict);
    CavitationModel cav(mesh, phi, schemes, cfg);

    scalar endTime = 1.0;
    int reportEvery = 10;
    PisoControls ctl = app3dturb::pisoControls(caseDict, endTime, reportEvery);

    {
        scalar vMin = 1e300, vMax = 0;
        for (const scalar v : mesh.cellVolumes())
        {
            vMin = std::min(vMin, v);
            vMax = std::max(vMax, v);
        }
        log::info("cavitation pSat ", app3dturb::num(cfg.pSat), " cells ", mesh.nCells(),
                  " volMin ", app3dturb::num(vMin), " volMax ", app3dturb::num(vMax));
    }

    Field<scalar> invRho(static_cast<std::size_t>(mesh.nCells()), 1.0);
    ctl.cellInvRho = &invRho;

    auto columnMean = [&](int i, scalar& uMean, scalar& pMean, scalar& aMean)
    {
        uMean = pMean = aMean = 0;
        if (nx <= 0) return;
        int cnt = 0;
        for (int j = 0; j < ny; ++j)
        {
            const std::size_t c = static_cast<std::size_t>(j * nx + i);
            uMean += U.internalField()[c].x;
            pMean += p.internalField()[c];
            aMean += cav.alphaCell()[c];
            ++cnt;
        }
        uMean /= cnt; pMean /= cnt; aMean /= cnt;
    };

    scalar t = 0;
    int step = 0;
    auto report = [&]()
    {
        scalar pMin = 1e300, pMax = -1e300;
        for (const scalar pv : p.internalField())
        {
            pMin = std::min(pMin, pv);
            pMax = std::max(pMax, pv);
        }
        std::printf("t %s pMin %s pMax %s alphaMax %s\n", app3dturb::num(t).c_str(),
                    app3dturb::num(pMin).c_str(), app3dturb::num(pMax).c_str(),
                    app3dturb::num(cav.alphaMax()).c_str());
        if (nx > 0)
        {
            const int throat = static_cast<int>(nx / 2);
            scalar u0, p0, a0, ut, pt, at, ue, pe, ae;
            columnMean(0, u0, p0, a0);
            columnMean(throat, ut, pt, at);
            columnMean(nx - 1, ue, pe, ae);
            std::printf("plane inlet u %s p %s a %s\nthroat u %s p %s a %s\n"
                        "plane outlet u %s p %s a %s\n",
                        app3dturb::num(u0).c_str(), app3dturb::num(p0).c_str(),
                        app3dturb::num(a0).c_str(), app3dturb::num(ut).c_str(),
                        app3dturb::num(pt).c_str(), app3dturb::num(at).c_str(),
                        app3dturb::num(ue).c_str(), app3dturb::num(pe).c_str(),
                        app3dturb::num(ae).c_str());
        }
    };

    report();
    while (t < endTime - 1e-12 * std::max(1.0, endTime))
    {
        const scalar dt = std::min(ctl.dt, endTime - t);
        PisoControls stepCtl = ctl;
        stepCtl.dt = dt;
        solveTransientIncompressiblePiso(mesh, U, p, phi, nu, schemes, stepCtl);
        t += dt;
        ++step;

        cav.correct(dt, p.internalField());
        const Field<scalar> rho = cav.mixtureDensity();
        for (std::size_t c = 0; c < invRho.size(); ++c)
        {
            invRho[c] = cfg.rhoL / rho[c];
        }
        if (step % reportEvery == 0) report();
    }
    report();

    if (caseDict.has("output"))
    {
        const Dictionary& od = caseDict.getDict("output");
        const std::string vtk = od.getWordOr("vtk", "");
        if (!vtk.empty())
        {
            VtkWriter w(mesh);
            w.addCellVector("U", U.internalField());
            w.addCellScalar("p", p.internalField());
            w.addCellScalar("alpha_v", cav.alphaCell());
            w.addCellScalar("rhoMix", cav.mixtureDensity());
            w.write(vtk);
        }
    }
    return 0;
}

// Uncaught cfFatal dies silently on MSVC (fail-fast, nothing on the pipe the
// checkers read); catch and report, the adaptiveRefineApp convention.
int main(int argc, char** argv)
{
    try
    {
        return appRun(argc, argv);
    }
    catch (const std::exception& e)
    {
        log::error("cf::Error: ", e.what());
        return 1;
    }
}
