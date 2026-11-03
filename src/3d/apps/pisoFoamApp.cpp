#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/PisoSolver.h"
#include "cf/fvm/VtkWriter.h"

using namespace cf;

namespace
{
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-piso-foam <case.dict>\n";
        std::exit(2);
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);

        polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
        log::info("mesh: ", mesh.nCells(), " cells, ", mesh.nPatches(), " patches");

        const scalar nu = caseDict.getScalar("nu");

        volVectorField U(mesh, "U");
        volScalarField p(mesh, "p");
        surfaceScalarField phi(mesh, "phi", 0.0);

        //   patches { inlet { U (1 0 0); p zeroGradient; } ... }
        if (caseDict.has("patches"))
        {
            const Dictionary& ps = caseDict.getDict("patches");
            for (const auto& name : ps.keys())
            {
                const Dictionary& pd = ps.getDict(name);
                const label patch = mesh.patchIdByName(name);
                if (pd.has("U"))
                {
                    const auto ul = pd.getScalarList("U");
                    const std::string type = pd.getWordOr("UType", "fixedValue");
                    U.setPatchValue(patch, type, Vector{ul[0], ul[1], ul[2]});
                }
                else
                {
                    U.setPatchValue(patch, pd.getWordOr("UType", "zeroGradient"), Vector::zero());
                }
                if (pd.has("p"))
                {
                    p.setPatchValue(patch, pd.getWordOr("pType", "fixedValue"), pd.getScalar("p"));
                }
                else
                {
                    p.setPatchValue(patch, "zeroGradient", 0.0);
                }
            }
        }

        Schemes schemes = loadSchemes(caseDict);
        schemes.ddt = "Euler";

        PisoControls ctl;
        scalar endTime = 1.0;
        int reportEvery = 10;
        if (caseDict.has("piso"))
        {
            const Dictionary& sd = caseDict.getDict("piso");
            ctl.dt = sd.getScalar("dt");
            endTime = sd.getScalar("endTime");
            ctl.correctors = static_cast<int>(sd.getScalarOr("correctors", ctl.correctors));
            ctl.continuityTol = sd.getScalarOr("continuityTol", ctl.continuityTol);
            ctl.linearTol = sd.getScalarOr("linearTol", ctl.linearTol);
            ctl.linearMaxIter = static_cast<int>(sd.getScalarOr("linearMaxIter", ctl.linearMaxIter));
            ctl.nonOrthCorrectors =
                static_cast<int>(sd.getScalarOr("nonOrthCorrectors", ctl.nonOrthCorrectors));
            reportEvery = static_cast<int>(sd.getScalarOr("reportEvery", reportEvery));
            if (sd.has("force"))
            {
                const auto f = sd.getScalarList("force");
                ctl.bodyForce = Vector{f[0], f[1], f[2]};
            }
        }
        else
        {
            cfFatal("piso: case needs a piso { dt; endTime; ... } block");
        }

        phi.correctBoundaryFlux(U);

        scalar t = 0;
        int step = 0, maxCorrectors = 0;
        scalar worstContinuity = 0;
        PisoStepResult res;

        while (t < endTime - 1e-12 * std::max(1.0, endTime))
        {
            const scalar dt = std::min(ctl.dt, endTime - t);
            PisoControls stepCtl = ctl;
            stepCtl.dt = dt;
            res = solveTransientIncompressiblePiso(mesh, U, p, phi, nu, schemes, stepCtl);
            t += dt;
            ++step;
            maxCorrectors = std::max(maxCorrectors, res.correctorsUsed);
            worstContinuity = std::max(worstContinuity, res.continuityError);

            if (step % reportEvery == 0 || t >= endTime - 1e-12 * endTime)
            {
                scalar uMax = 0;
                for (const Vector& u : U.internalField()) uMax = std::max(uMax, u.x);
                log::info("t ", num(t), " uMax ", num(uMax),
                          " cont ", num(res.continuityError),
                          " correctors ", res.correctorsUsed);
            }
        }

        // final profile on the first x-column: the checker compares it against
        // the analytic parabola
        const auto& cc = mesh.cellCentres();
        scalar x0 = 1e300;
        for (label c = 0; c < mesh.nCells(); ++c)
        {
            x0 = std::min(x0, cc[static_cast<std::size_t>(c)].x);
        }
        for (label c = 0; c < mesh.nCells(); ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            if (std::abs(cc[u].x - x0) <= 1e-9 * std::max(1.0, std::abs(x0)))
            {
                log::info("profile y ", num(cc[u].y), " u ", num(U.internalField()[u].x));
            }
        }

        log::info("piso done: steps ", step, " endTime ", num(t),
                  " maxCorrectors ", maxCorrectors,
                  " worstContinuity ", num(worstContinuity));

        if (caseDict.has("output"))
        {
            const Dictionary& od = caseDict.getDict("output");
            const std::string vtk = od.getWordOr("vtk", "");
            if (!vtk.empty())
            {
                VtkWriter w(mesh);
                w.addCellVector("U", U.internalField());
                w.addCellScalar("p", p.internalField());
                w.write(vtk);
            }
        }

        return std::isfinite(worstContinuity) ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        log::error("pisoFoam: ", e.what());
        return 1;
    }
}
