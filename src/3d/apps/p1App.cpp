// chflow-3d-p1: steady P1 radiative transfer through an absorbing,
// isotropically scattering grey participating medium (D8 increment 4).
//
// The app is deliberately independent of the flow/energy solvers: the medium
// temperature is uniform and declared in the case, black walls are prescribed
// through their wallTemperature, and the P1 field is solved once. The published
// P1 equation, the black-body wall condition G_w = 4 sigma T_w^4 and the exact
// domain balance are in Radiation.h; the checker reconstructs all of them
// from the case constants and the printed field without trusting this solver.
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/Radiation.h"
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
        std::cerr << "usage: chflow-3d-p1 <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const RadiationConfig rad = readRadiationConfig(caseDict);
        if (!rad.enabled) cfFatal("p1App: needs a radiation{} declaration");

        polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
        log::info("p1 mesh: ", mesh.nCells(), " cells, ", mesh.nFaces(), " faces, ",
                  mesh.nPatches(), " patches");

        volScalarField G(mesh, "G", 0.0);

        // Black walls: G_wall = 4 sigma T_wall^4. Anything without a
        // wallTemperature is adiabatic to radiation (zeroGradient).
        if (caseDict.has("patches"))
        {
            const Dictionary& ps = caseDict.getDict("patches");
            for (const auto& name : ps.keys())
            {
                const Dictionary& pd = ps.getDict(name);
                const label patch = mesh.patchIdByName(name);
                if (pd.has("wallTemperature"))
                {
                    const scalar tw = pd.getScalar("wallTemperature");
                    if (!(tw > 0) || !std::isfinite(tw))
                        cfFatal("p1App: wallTemperature on patch '" + name +
                                "' must be finite and positive");
                    const scalar gw = blackbodyIncidentRadiation(tw, rad.sbConstant);
                    G.setPatchValue(patch, "fixedValue", gw);
                    log::info("p1 wall ", name, ": T ", num(tw), " K -> G ", num(gw), " W/m2");
                }
                else
                {
                    G.setPatchValue(patch, "zeroGradient", 0.0);
                }
            }
        }

        const Field<scalar> mediumT(static_cast<std::size_t>(mesh.nCells()),
                                    rad.mediumTemperature);
        const Schemes schemes = loadSchemes(caseDict);

        log::info("p1 model: kappa_a ", num(rad.absorption), " kappa_s ",
                  num(rad.scattering), " beta ", num(rad.extinction()),
                  " gamma ", num(rad.diffusion()), " T ", num(rad.mediumTemperature),
                  " sigma ", num(rad.sbConstant));

        const RadiationResult res = solveP1(mesh, G, rad, schemes, mediumT);
        if (!res.solved)
        {
            log::error("p1App: P1 solve did not converge (", res.iterations,
                       " iters, residual ", num(res.residual), ")");
            return 1;
        }

        log::info("p1 solve: ", res.iterations, " iters, residual ", num(res.residual));

        for (label c = 0; c < mesh.nCells(); ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            const Vector& q = res.qRad[u];
            log::info("p1 cell ", c, ": G ", num(res.G[u]),
                      " q ", num(q.x), " ", num(q.y), " ", num(q.z),
                      " S ", num(res.netSource[u]));
        }
        log::info("p1 balance: volumeSource ", num(res.volumeSource),
                  " boundaryFlux ", num(res.boundaryFlux),
                  " sum ", num(res.balance));

        scalar gMin = res.G.empty() ? 0 : res.G[0];
        scalar gMax = gMin;
        for (scalar g : res.G)
        {
            gMin = std::min(gMin, g);
            gMax = std::max(gMax, g);
        }
        log::info("p1 summary: cells ", mesh.nCells(), " Gmin ", num(gMin),
                  " Gmax ", num(gMax), " 4sigmaT4 ",
                  num(blackbodyIncidentRadiation(rad.mediumTemperature, rad.sbConstant)));

        if (caseDict.has("output"))
        {
            const std::string vtk = caseDict.getDict("output").getWordOr("vtk", "");
            if (!vtk.empty())
            {
                VtkWriter w(mesh);
                w.addCellScalar("G", G.internalField());
                w.write(vtk);
            }
        }
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("p1App: ", e.what());
        return 1;
    }
}
