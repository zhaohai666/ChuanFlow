// chflow-3d-dpm: solve a steady carrier flow once, then track Lagrangian
// spheres through the frozen field. D8 increment 2 is strictly one-way:
// particles sample U but never write momentum back. The published
// Schiller-Naumann sphere law and its Stokes limit set the drag acceleration.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/Dpm.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/SimpleSolver.h"
#include "cf/fvm/VtkWriter.h"

#include "caseInit.h"

namespace
{

std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}

std::string vec(const cf::Vector& v)
{
    return "(" + num(v.x) + " " + num(v.y) + " " + num(v.z) + ")";
}

struct CarrierCheck
{
    double uMax = 0;
    double sumU = 0;
};

CarrierCheck carrierState(const cf::volVectorField& U)
{
    CarrierCheck c;
    for (const cf::Vector& u : U.internalField())
    {
        c.uMax = std::max(c.uMax, static_cast<double>(mag(u)));
        c.sumU += static_cast<double>(mag(u));
    }
    return c;
}

} // namespace

using namespace cf;

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-3d-dpm <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
        const scalar nu = caseDict.getScalar("nu");
        if (!(nu > 0)) cfFatal("dpmApp: nu must be positive");

        volVectorField U(mesh, "U");
        volScalarField p(mesh, "p");
        surfaceScalarField phi(mesh, "phi", 0.0);

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
                    U.setPatchValue(patch, pd.getWordOr("UType", "fixedValue"),
                                    Vector{ul[0], ul[1], ul[2]});
                }
                else
                {
                    U.setPatchValue(patch, pd.getWordOr("UType", "zeroGradient"),
                                    Vector::zero());
                }
                if (pd.has("p"))
                    p.setPatchValue(patch, pd.getWordOr("pType", "fixedValue"),
                                    pd.getScalar("p"));
                else
                    p.setPatchValue(patch, "zeroGradient", 0.0);
            }
        }
        phi.correctBoundaryFlux(U);
        app::seedInitialVelocity(caseDict, mesh, U, phi);

        const Schemes schemes = loadSchemes(caseDict);
        SimpleControls ctl;
        if (caseDict.has("simple"))
        {
            const Dictionary& sd = caseDict.getDict("simple");
            ctl.maxIterations = static_cast<int>(sd.getScalarOr("maxIterations", ctl.maxIterations));
            ctl.continuityTol = sd.getScalarOr("continuityTol", ctl.continuityTol);
            ctl.pcorrTol = sd.getScalarOr("pcorrTol", ctl.pcorrTol);
            ctl.momentumTol = sd.getScalarOr("momentumTol", ctl.momentumTol);
            ctl.relaxU = sd.getScalarOr("relaxU", ctl.relaxU);
            ctl.relaxP = sd.getScalarOr("relaxP", ctl.relaxP);
            ctl.pRelaxStart = sd.getScalarOr("pRelaxStart", ctl.pRelaxStart);
            ctl.nonOrthCorrectors =
                static_cast<int>(sd.getScalarOr("nonOrthCorrectors", ctl.nonOrthCorrectors));
        }
        const SimpleResult flow =
            solveIncompressibleSimple(mesh, U, p, phi, nu, schemes, ctl);
        if (!flow.converged)
        {
            log::error("dpmApp: carrier flow did not converge");
            return 1;
        }
        log::info("dpm carrier: iterations ", flow.iterations, " continuity ",
                  num(flow.finalContinuityError));

        const DpmConfig cfg = readDpmConfig(caseDict);
        const CarrierCheck before = carrierState(U);
        const DpmResult result = runDpm(mesh, U, cfg, cfg.carrierRho * nu);
        const CarrierCheck after = carrierState(U);
        log::info("dpm carrier frozen: nu ", num(nu), " mu ",
                  num(cfg.carrierRho * nu), " uMaxBefore ", num(before.uMax),
                  " uMaxAfter ", num(after.uMax), " sumUBefore ", num(before.sumU),
                  " sumUAfter ", num(after.sumU));

        for (const DpmParticle& particle : result.particles)
        {
            log::info("dpm launch ", particle.id, ": stream ", particle.stream,
                      " x ", vec(particle.x0), " u ", vec(particle.u0),
                      " diameter ", num(particle.diameter),
                      " rho ", num(particle.rhoParticle));
        }
        for (const DpmSample& sample : result.samples)
        {
            log::info("dpm sample ", sample.id, " ", sample.step, ": t ", num(sample.t),
                      " x ", vec(sample.x), " u ", vec(sample.u),
                      " uf ", vec(sample.uFluid), " Re ", num(sample.reynolds),
                      " Cd ", num(sample.cd),
                      " aDrag ", vec(sample.dragAcceleration));
        }
        for (const DpmParticle& particle : result.particles)
        {
            log::info("dpm particle ", particle.id, " end: status ",
                      dpmStatusName(particle.status), " t ", num(particle.age),
                      " x ", vec(particle.x), " u ", vec(particle.u));
        }
        log::info("dpm summary: particles ", result.particles.size(),
                  " outlet ", result.outletCount, " wall ", result.wallCount,
                  " expired ", result.expiredCount);

        if (caseDict.has("output"))
        {
            const std::string vtk = caseDict.getDict("output").getWordOr("vtk", "");
            if (!vtk.empty())
            {
                VtkWriter w(mesh);
                w.addCellVector("U", U.internalField());
                w.addCellScalar("p", p.internalField());
                w.write(vtk);
            }
        }
        return result.expiredCount == 0 ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        log::error("dpmApp: ", e.what());
        return 1;
    }
}
