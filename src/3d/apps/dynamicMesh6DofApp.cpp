// chflow-3d-dynamic-mesh: GCL-corrected relative-flux scalar transport on a
// rigidly translating mesh, integrated alongside a rigid body with one
// translational DOF (mass + linear drag + constant force) and one rotational
// DOF (inertia + rotational drag + stiffness). See DynamicMesh6Dof.h for the
// discrete formulation and the documented subset of full 6DOF.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/DynamicMesh6Dof.h"
#include "cf/fvm/MeshLoader.h"

using namespace cf;

namespace
{
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}

void report(scalar t, const BodyState& s, scalar scalarMean, scalar massBalance)
{
    std::ostringstream os;
    os << std::setprecision(12)
       << "DYNMESH " << t << ' ' << s.x << ' ' << s.v << ' '
       << s.theta << ' ' << s.omega << ' ' << scalarMean << ' ' << massBalance;
    log::info(os.str());
}

// §5.236: printed only when fluidStress{} is declared (passive twin keeps the
// F11 output byte-identical).
void reportStress(scalar t, const StressLoad& sl)
{
    std::ostringstream os;
    os << std::setprecision(12)
       << "STRESS " << t
       << ' ' << sl.forcePressure.x << ' ' << sl.forcePressure.y << ' ' << sl.forcePressure.z
       << ' ' << sl.forceViscous.x << ' ' << sl.forceViscous.y << ' ' << sl.forceViscous.z
       << ' ' << sl.moment.x << ' ' << sl.moment.y << ' ' << sl.moment.z
       << ' ' << sl.volume;
    log::info(os.str());
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-3d-dynamic-mesh <case.dict>\n";
        return 2;
    }
    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const DynamicMeshConfig cfg = readDynamicMeshConfig(caseDict);
        if (!cfg.enabled) cfFatal("dynamicMesh: needs a dynamicMesh{} declaration");

        polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
        log::info("dynamicMesh: ", mesh.nCells(), " cells, U_mesh (",
                  num(cfg.meshVelocity.x), ' ', num(cfg.meshVelocity.y), ' ',
                  num(cfg.meshVelocity.z), ") dt ", num(cfg.dt), " steps ", cfg.steps);

        Field<scalar> phi(static_cast<std::size_t>(mesh.nCells()), cfg.initialValue);
        Field<scalar> oldVol = mesh.cellVolumes();

        BodyState body;
        body.x = 0;
        body.v = 0;
        body.theta = cfg.theta0;
        body.omega = cfg.omega0;

        // Initial state (t = 0): no motion has happened, so GCL residual is 0.
        report(0.0, body, cfg.initialValue, 0.0);
        if (cfg.stress.enabled)
            reportStress(0.0, integrateFluidStress(mesh, cfg.stress));

        for (int step = 1; step <= cfg.steps; ++step)
        {
            // Rigid-body update (independent of the mesh geometry). With
            // fluidStress{} the body is advanced AFTER the mesh move, driven
            // by the integrated surface traction instead of bodyForce / drag.
            if (!cfg.stress.enabled)
                body = advanceBody(body, cfg, cfg.dt);

            // Prescribed rigid mesh translation.
            Field<Vector> newPoints(static_cast<std::size_t>(mesh.nPoints()));
            const auto& pts = mesh.points();
            for (std::size_t i = 0; i < pts.size(); ++i)
                newPoints[i] = pts[i] + cfg.meshVelocity * cfg.dt;
            mesh.movePoints(newPoints);

            StressLoad load;
            if (cfg.stress.enabled)
            {
                load = integrateFluidStress(mesh, cfg.stress);
                body = advanceBodyStress(body, cfg, load, cfg.dt);
            }

            // Relative flux and GCL on the NEW geometry.
            const Field<scalar> phiRel = relativeFaceFlux(mesh, cfg.fluidVelocity, cfg.meshVelocity);
            const Field<scalar> sweep = gclGeometricSweep(mesh, cfg.meshVelocity);
            const Field<scalar> gclRes = gclCellResidual(mesh, cfg.meshVelocity, oldVol, cfg.dt);
            scalar massBalance = 0;
            for (scalar r : gclRes) massBalance = std::max(massBalance, std::abs(r));

            phi = transportStep(mesh, phi, phiRel, sweep, oldVol, cfg.dt);
            oldVol = mesh.cellVolumes();

            scalar totalMass = 0, totalVol = 0;
            const auto& vol = mesh.cellVolumes();
            for (label c = 0; c < mesh.nCells(); ++c)
            {
                totalMass += phi[static_cast<std::size_t>(c)] * vol[static_cast<std::size_t>(c)];
                totalVol += vol[static_cast<std::size_t>(c)];
            }
            const scalar scalarMean = totalMass / totalVol;

            if (step % cfg.reportEvery == 0 || step == cfg.steps)
            {
                report(step * cfg.dt, body, scalarMean, massBalance);
                if (cfg.stress.enabled)
                    reportStress(step * cfg.dt, load);
            }
        }
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("dynamicMesh: ", e.what());
        return 1;
    }
}
