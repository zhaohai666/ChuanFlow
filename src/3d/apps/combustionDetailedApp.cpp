// chflow-3d-combustion-detailed (batch D / F10): two-step chain-branching
// chemistry in a homogeneous 0D reactor. This is a distinct kernel from the
// steady single-step chflow-3d-combustion; it resolves ignition delay and the
// extinction Damkoehler threshold as simulated outcomes. The stdout stream is
// a case-echo line, a diagnostics line and one `COMBUST t T Yf Yr` row per
// output step, all at 12 significant digits for independent grading.
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/fvm/CombustionDetailed.h"

using namespace cf;

namespace
{

std::string num(scalar v)
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
        std::cerr << "usage: chflow-3d-combustion-detailed <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const DetailedChemistryConfig cfg = readDetailedChemistryConfig(caseDict);
        if (!cfg.enabled)
            cfFatal("combustionDetailedApp: needs a detailedChemistry{} declaration");

        const auto hist = integrateDetailed(cfg);
        const DetailedDiagnostics d = diagnoseDetailed(cfg);

        std::cout << "DETAILED params"
                  << " A1 " << num(cfg.A1) << " E1 " << num(cfg.E1)
                  << " A2 " << num(cfg.A2) << " E2 " << num(cfg.E2)
                  << " q " << num(cfg.q) << " T0 " << num(cfg.T0)
                  << " gasConstant " << num(cfg.gasConstant)
                  << " Yf0 " << num(cfg.Yf0) << " Yr0 " << num(cfg.Yr0)
                  << " cooling " << num(cfg.cooling) << " Tenv " << num(cfg.Tenv)
                  << " dt " << num(cfg.dt) << " nSteps " << cfg.nSteps
                  << " nSub " << cfg.nSub
                  << " ignitionRise " << num(cfg.ignitionRise) << "\n";

        std::cout << "DETAILED number"
                  << " branchDa " << num(d.branchDa)
                  << " betaCrit " << num(d.betaCrit)
                  << " DaCrit " << num(d.DaCrit)
                  << " ignitionDelay " << num(d.ignitionDelay)
                  << " peakT " << num(d.peakT)
                  << " ignited " << (d.ignited ? 1 : 0) << "\n";

        for (const auto& s : hist)
            std::cout << "COMBUST " << num(s.t) << " " << num(s.T) << " "
                      << num(s.Yf) << " " << num(s.Yr) << "\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "combustionDetailedApp: " << e.what() << "\n";
        return 1;
    }
}
