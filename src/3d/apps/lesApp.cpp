// chflow-3d-les: large-eddy simulation with the algebraic Smagorinsky or WALE
// sub-grid viscosity, advanced by the transient incompressible PISO coupling.
//
// Opt-in case keys:
//   les { model smagorinsky | wale; Cs 0.17; Cw 0.325; }
//   profile { type linearShear | cosineShear | uniaxialStrain; rate ...; }
//           - an analytic starting field instead of the uniform `initial U`:
//             linearShear   U = (G y, 0, 0)
//             cosineShear   U = (A cos(pi y / H), 0, 0)
//             uniaxialStrain U = (a x, -a y/2, -a z/2)
//
// Reports, per `piso{ reportEvery }` steps and at the end:
//   time, nut min/max/mean, |S| max, filter width, and the volume-averaged
//   kinetic energy E = <0.5 |U|^2>, which is what the decay regressions read.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/turbulence/LesModel.h"
#include "case3dTurb.h"

using namespace cf;

static int appRun(int argc, char** argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr, "usage: chflow-3d-les <case.dict>\n");
        return 1;
    }
    const Dictionary caseDict = Dictionary::parseFile(argv[1]);

    polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
    const scalar nu = app3dturb::caseNu(caseDict);

    volVectorField U(mesh, "U");
    volScalarField p(mesh, "p");
    surfaceScalarField phi(mesh, "phi", 0.0);
    app3dturb::applyPatches(caseDict, mesh, U, p);

    Schemes schemes = loadSchemes(caseDict);
    schemes.ddt = "Euler";

    // analytic starting profile
    if (caseDict.has("profile"))
    {
        const Dictionary& pd = caseDict.getDict("profile");
        const std::string type = pd.getWord("type");
        const scalar rate = pd.getScalarOr("rate", 1.0);
        const scalar amp = pd.getScalarOr("amplitude", 1.0);
        const Vector lo = mesh.minX(), hi = mesh.maxX();
        const scalar H = hi.y - lo.y;
        for (std::size_t c = 0; c < U.internalField().size(); ++c)
        {
            const Vector& x = mesh.cellCentres()[c];
            if (type == "linearShear")
            {
                U.internalField()[c] = Vector{rate * (x.y - lo.y), 0, 0};
            }
            else if (type == "cosineShear")
            {
                U.internalField()[c] =
                    Vector{amp * std::cos(M_PI * (x.y - lo.y) / H), 0, 0};
            }
            else if (type == "uniaxialStrain")
            {
                U.internalField()[c] = Vector{rate * (x.x - lo.x),
                                              -0.5 * rate * (x.y - lo.y),
                                              -0.5 * rate * (x.z - lo.z)};
            }
            else
            {
                cfFatal("profile: unknown type '" + type + "'");
            }
        }
    }

    phi.correctBoundaryFlux(U);
    const LesModel::Config lesCfg = LesModel::Config::fromCase(caseDict);
    LesModel les(mesh, U, phi, nu, schemes, lesCfg);

    scalar endTime = 1.0;
    int reportEvery = 10;
    PisoControls ctl = app3dturb::pisoControls(caseDict, endTime, reportEvery);
    if (caseDict.has("piso"))
    {
        const Dictionary& sd = caseDict.getDict("piso");
        endTime = sd.getScalarOr("endTime", endTime);
    }

    log::info("les model ", les.name(), " cells ", mesh.nCells());

    scalar t = 0;
    int step = 0;
    auto report = [&](bool final)
    {
        const auto& nut = les.nut();
        const auto& s = les.strainRate();
        scalar nutMin = 1e300, nutMax = 0, nutSum = 0, sMax = 0;
        scalar energy = 0, volSum = 0;
        for (std::size_t c = 0; c < nut.size(); ++c)
        {
            nutMin = std::min(nutMin, nut[c]);
            nutMax = std::max(nutMax, nut[c]);
            nutSum += nut[c];
            sMax = std::max(sMax, s[c]);
            energy += 0.5 * mag2(U.internalField()[c]) * mesh.cellVolumes()[c];
            volSum += mesh.cellVolumes()[c];
        }
        std::printf("t %s nutMin %s nutMax %s nutMean %s strainMax %s E %s\n",
                    app3dturb::num(t).c_str(), app3dturb::num(nutMin).c_str(),
                    app3dturb::num(nutMax).c_str(),
                    app3dturb::num(nutSum / static_cast<scalar>(nut.size())).c_str(),
                    app3dturb::num(sMax).c_str(),
                    app3dturb::num(energy / volSum).c_str());
        if (final)
        {
            const auto& delta = les.filterWidth();
            std::printf("filterWidth %s\n", app3dturb::num(delta[0]).c_str());
        }
    };

    // evaluate the algebraic viscosity on the seeded field once, so a case
    // that runs zero steps still reports the nut the profile implies
    les.correctTurbulence(0.0);
    report(false);
    while (t < endTime - 1e-12 * std::max(1.0, endTime))
    {
        const scalar dt = std::min(ctl.dt, endTime - t);
        PisoControls stepCtl = ctl;
        stepCtl.dt = dt;
        solveTransientIncompressiblePiso(mesh, U, p, phi, nu, schemes, stepCtl, &les);
        t += dt;
        ++step;
        if (step % reportEvery == 0) report(false);
    }
    report(true);

    if (caseDict.has("output"))
    {
        const Dictionary& od = caseDict.getDict("output");
        const std::string vtk = od.getWordOr("vtk", "");
        if (!vtk.empty())
        {
            VtkWriter w(mesh);
            w.addCellVector("U", U.internalField());
            w.addCellScalar("p", p.internalField());
            w.addCellScalar("nut", les.nut());
            w.addCellScalar("strainRate", les.strainRate());
            w.addCellScalar("filterWidth", les.filterWidth());
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
