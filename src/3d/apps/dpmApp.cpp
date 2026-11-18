// chflow-3d-dpm: solve a steady carrier flow once, then track Lagrangian
// spheres through the frozen field. D8 increment 2 is strictly one-way:
// particles sample U but never write momentum back. The published
// Schiller-Naumann sphere law and its Stokes limit set the drag acceleration.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
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

/// Turn the accumulated drag reaction impulse into the cell acceleration the
/// steady momentum equation understands. The averaging window is the longest
/// particle lifetime, so the impulse of a short-lived cloud is not diluted by
/// an arbitrary maxSteps*dt horizon.
cf::Field<cf::Vector> forceFromImpulse(const cf::polyMesh& mesh,
                                       const cf::DpmResult& result,
                                       cf::scalar rhoF)
{
    cf::Field<cf::Vector> f(static_cast<std::size_t>(mesh.nCells()), cf::Vector::zero());
    if (!(result.sourceWindow > 0)) return f;
    const auto& vol = mesh.cellVolumes();
    for (cf::label c = 0; c < mesh.nCells(); ++c)
    {
        f[static_cast<std::size_t>(c)] =
            -1.0 * result.cellDragImpulse[static_cast<std::size_t>(c)]
            / (rhoF * vol[static_cast<std::size_t>(c)] * result.sourceWindow);
    }
    return f;
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
        const scalar mu = cfg.carrierRho * nu;
        CarrierCheck before;
        DpmResult result;

        if (cfg.twoWay)
        {
            const CarrierCheck baseline = carrierState(U);
            Field<Vector> force(static_cast<std::size_t>(mesh.nCells()), Vector::zero());
            scalar lastResidual = std::numeric_limits<scalar>::infinity();
            int usedOuter = 0;
            for (int outer = 1; outer <= cfg.outerIterations; ++outer)
            {
                SimpleControls twoCtl = ctl;
                twoCtl.cellForce = &force;
                const SimpleResult coupled =
                    solveIncompressibleSimple(mesh, U, p, phi, nu, schemes, twoCtl);
                if (!coupled.converged)
                {
                    log::error("dpmApp: two-way carrier flow did not converge");
                    return 1;
                }
                const DpmResult trial = runDpm(mesh, U, cfg, mu);
                const Field<Vector> source = forceFromImpulse(mesh, trial, cfg.carrierRho);
                scalar delta = 0, scale = 0;
                for (label c = 0; c < mesh.nCells(); ++c)
                {
                    const std::size_t u = static_cast<std::size_t>(c);
                    delta = std::max(delta, static_cast<scalar>(mag(source[u] - force[u])));
                    scale = std::max(scale, static_cast<scalar>(mag(source[u])));
                }
                lastResidual = delta / std::max(scale, 1.0e-30);
                for (label c = 0; c < mesh.nCells(); ++c)
                {
                    const std::size_t u = static_cast<std::size_t>(c);
                    force[u] = (1.0 - cfg.twoWayRelax) * force[u]
                               + cfg.twoWayRelax * source[u];
                }
                log::info("dpm twoWay outer ", outer, " residual ", num(lastResidual),
                          " maxForce ", num(scale));
                usedOuter = outer;
                if (lastResidual <= cfg.twoWayTol) break;
            }
            log::info("dpm twoWay outerUsed ", usedOuter, " residual ", num(lastResidual));
            const CarrierCheck finalCarrier = carrierState(U);
            log::info("dpm twoWay carrier baseline uMax ", num(baseline.uMax),
                      " sumU ", num(baseline.sumU), " final uMax ", num(finalCarrier.uMax),
                      " final sumU ", num(finalCarrier.sumU));
            before = carrierState(U);
            result = runDpm(mesh, U, cfg, mu);
            const Field<Vector> finalSource = forceFromImpulse(mesh, result, cfg.carrierRho);
            for (label c = 0; c < mesh.nCells(); ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                if (mag(result.cellDragImpulse[u]) > 0)
                    log::info("dpm source cell ", c, " impulse ",
                              vec(result.cellDragImpulse[u]), " force ", vec(finalSource[u]));
            }
        }
        else
        {
            before = carrierState(U);
            result = runDpm(mesh, U, cfg, mu);
        }

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
