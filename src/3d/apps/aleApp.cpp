// chflow-3d-ale: moving-mesh geometry and the Geometric Conservation Law
// (D8 increment 8). The mesh vertices move by a rigid translation or a uniform
// expansion; after each move the app verifies, per cell,
//   dV/dt = sum_faces U_mesh . Sf
// with the recomputed geometry. A rigid translation satisfies this to
// roundoff; a uniform expansion has a first-order-in-dt truncation that the
// checker measures and shows shrinking as dt is halved.
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/Ale.h"
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
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-3d-ale <case.dict>\n";
        return 2;
    }
    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const AleConfig cfg = readAleConfig(caseDict);
        if (!cfg.enabled) cfFatal("aleApp: needs an ale{} declaration");

        polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
        log::info("ale mesh: ", mesh.nCells(), " cells, ", mesh.nPoints(), " points, motion ",
                  cfg.motion, " dt ", num(cfg.dt), " steps ", cfg.steps);

        Field<scalar> oldVol = mesh.cellVolumes();
        Field<Vector> oldPoints = mesh.points();
        scalar totalV0 = 0;
        for (scalar v : oldVol) totalV0 += v;
        log::info("ale initial: totalVolume ", num(totalV0));

        for (int step = 1; step <= cfg.steps; ++step)
        {
            // Vertex velocity evaluated at the CURRENT positions, then the
            // mesh advances by dt (the velocity is constant over the step for
            // translation and grows with x for expansion).
            const Field<Vector> vVert = meshVertexVelocity(mesh, cfg);
            const Field<Vector> newPoints = advanceVertices(mesh, cfg);
            mesh.movePoints(newPoints);
            const Field<scalar> res = gclResidual(mesh, oldVol, vVert, cfg.dt);

            scalar maxRes = 0, sumRes = 0, totalV = 0;
            for (scalar r : res) { maxRes = std::max(maxRes, std::abs(r)); sumRes += r; }
            for (scalar v : mesh.cellVolumes()) totalV += v;

            if (step % cfg.reportEvery == 0 || step == cfg.steps)
            {
                log::info("ale step ", step, ": t ", num(step * cfg.dt),
                          " totalVolume ", num(totalV), " maxGcl ", num(maxRes),
                          " sumGcl ", num(sumRes),
                          " xMin ", num(mesh.minX().x), " xMax ", num(mesh.maxX().x));
            }

            oldVol = mesh.cellVolumes();
            oldPoints = mesh.points();
        }
        (void)oldPoints;
        log::info("ale summary: steps ", cfg.steps, " motion ", cfg.motion);
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("aleApp: ", e.what());
        return 1;
    }
}
