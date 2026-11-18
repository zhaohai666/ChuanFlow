// chflow-3d-cht: conjugate heat transfer between a solid and a fluid region,
// each its own mesh, joined by a Dirichlet--Neumann interface iteration
// (D8 increment 7). The checker reconstructs the 1D two-layer series-resistance
// closed form: interface temperature, heat flux and each layer's linear
// profile, plus flux continuity at the interface.
#include <cstdlib>
#include <iomanip>
#include <sstream>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/BlockMesh.h"
#include "cf/fvm/Cht.h"
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
        std::cerr << "usage: chflow-3d-cht <case.dict>\n";
        return 2;
    }
    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const ChtConfig cfg = readChtConfig(caseDict);
        if (!cfg.enabled) cfFatal("chtApp: needs a cht{} declaration");
        const Dictionary& cd = caseDict.getDict("cht");
        const Dictionary& sd = cd.getDict("solid");
        const Dictionary& fd = cd.getDict("fluid");

        const std::string sOuter = sd.getWordOr("outerPatch", "left");
        const std::string sIface = sd.getWordOr("interfacePatch", "right");
        const std::string fOuter = fd.getWordOr("outerPatch", "right");
        const std::string fIface = fd.getWordOr("interfacePatch", "left");

        polyMesh solidMesh = BlockMesh::fromDict(sd.getDict("block"));
        polyMesh fluidMesh = BlockMesh::fromDict(fd.getDict("block"));
        log::info("cht mesh: solid ", solidMesh.nCells(), " cells, fluid ",
                  fluidMesh.nCells(), " cells");

        const Schemes schemes = loadSchemes(caseDict);
        const ChtResult res = solveCht(solidMesh, fluidMesh, schemes, cfg,
                                       sOuter, sIface, fOuter, fIface);
        if (!res.solved)
        {
            log::error("chtApp: interface did not converge (mismatch ",
                       num(res.fluxMismatch), " after ", res.outerIterations, " iters)");
            return 1;
        }

        log::info("cht interface: T ", num(res.interfaceTemperature),
                  " fluxMismatch ", num(res.fluxMismatch),
                  " outerIters ", res.outerIterations);
        log::info("cht flux: solid outer ", num(res.solid.outerFlux),
                  " solid iface ", num(res.solid.interfaceFlux),
                  " fluid iface ", num(res.fluid.interfaceFlux),
                  " fluid outer ", num(res.fluid.outerFlux));

        for (label c = 0; c < solidMesh.nCells(); ++c)
            log::info("cht solid cell ", c, ": T ", num(res.solid.T[static_cast<std::size_t>(c)]));
        for (label c = 0; c < fluidMesh.nCells(); ++c)
            log::info("cht fluid cell ", c, ": T ", num(res.fluid.T[static_cast<std::size_t>(c)]));

        log::info("cht summary: solidIters ", res.solid.iterations,
                  " fluidIters ", res.fluid.iterations, " solved ", res.solved ? 1 : 0);
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("chtApp: ", e.what());
        return 1;
    }
}
