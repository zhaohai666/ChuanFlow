// chflow-3d-adaptive-refine: standalone verifier for the F13 adaptive
// refinement operator. It reads a case dict, refines the block hex mesh by
// the declared scalar-gradient threshold, rebuilds a conforming polyMesh and
// prints a 12-digit conservation report. A Python checker re-derives the
// refinement, the cell volumes and the exact field integral independently.
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/AdaptiveRefine.h"

using namespace cf;

namespace
{

/// Fixed 12 decimal places: volume sums and integrals are checkable to 1e-9.
std::string num(scalar v)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(12) << v;
    return os.str();
}

std::string vec(const Vector& v)
{
    return "(" + num(v.x) + " " + num(v.y) + " " + num(v.z) + ")";
}

void usage()
{
    std::cerr << "usage: chflow-3d-adaptive-refine <case.dict>\n"
              << "  block  { min (..); max (..); cells (..); <sides>; }\n"
              << "  adapt  { gradientThreshold <scalar>; maxDepth <int>;\n"
              << "           field linear|step; ...; }\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        usage();
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const AdaptResult r = AdaptiveRefine::refine(caseDict, argv[1]);
        const bool dumpCells = caseDict.has("adapt")
            && caseDict.getDict("adapt").getWordOr("dumpCells", "no") == "yes";

        log::info("adaptiveRefine: nPoints ", r.mesh->nPoints(),
                  " nCells ", r.mesh->nCells(),
                  " nFaces ", r.mesh->nFaces(),
                  " nInternalFaces ", r.mesh->nInternalFaces(),
                  " nBoundaryFaces ", r.mesh->nBoundaryFaces(),
                  " nPatches ", r.mesh->nPatches());
        for (const auto& p : r.mesh->patches())
            log::info("adaptiveRefine: patch ", p.name, " type ", p.type,
                      " start ", p.start, " size ", p.size);

        // The mandated 12-digit report line.
        log::info("ADAPT level=", r.level,
                  " cells=", r.nCells,
                  " volumeSum=", num(r.volumeSum),
                  " integral=", num(r.integral),
                  " integralError=", num(r.integralError),
                  " minVolume=", num(r.minVolume),
                  " maxVolume=", num(r.maxVolume));

        if (dumpCells)
        {
            for (std::size_t c = 0; c < r.cells.size(); ++c)
            {
                const AdaptCell& cell = r.cells[c];
                log::info("ADAPT_CELL c=", static_cast<label>(c),
                          " level=", cell.level,
                          " lo=", vec(cell.lo),
                          " hi=", vec(cell.hi),
                          " phi=", num(cell.phi),
                          " vol=", num(cell.volume));
            }
        }

        // §5.229 opt-in AMR remeshing loop. Without the block the output above
        // is byte-identical to the F13 reference.
        if (caseDict.has("amr"))
        {
            const AmrConfig amr = AmrConfig::fromDict(caseDict);
            const AmrLoopResult loop = AdaptiveRefine::remesh(caseDict, r);
            for (const AmrCycleReport& rep : loop.cycles)
            {
                log::info("AMR_CYCLE cycle=", rep.cycle,
                          " cells=", rep.nCells,
                          " refined=", rep.nRefined,
                          " coarsened=", rep.nCoarsened,
                          " mass=", num(rep.mass));
                if (amr.dumpCells)
                {
                    for (std::size_t c = 0; c < rep.cells.size(); ++c)
                    {
                        const AdaptCell& cell = rep.cells[c];
                        log::info("AMR_CELL cycle=", rep.cycle,
                                  " c=", static_cast<label>(c),
                                  " level=", cell.level,
                                  " lo=", vec(cell.lo),
                                  " hi=", vec(cell.hi),
                                  " phi=", num(cell.phi),
                                  " vol=", num(cell.volume));
                    }
                }
            }
            log::info("AMR_DONE initialCells=", loop.initialCells,
                      " finalCells=", loop.finalCells,
                      " massInitial=", num(loop.massInitial),
                      " massFinal=", num(loop.massFinal),
                      " massError=", num(loop.massError));
        }

        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("adaptiveRefine: ", e.what());
        return 1;
    }
}
