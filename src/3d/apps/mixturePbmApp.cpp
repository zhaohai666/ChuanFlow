// chflow-3d-mixture-pbm: opt-in algebraic mixture slip + monodisperse
// three-moment population balance (batch D item F7). Spatially uniform, no
// flow field: this is the moment-source kernel a DPM/species framework calls.
// Every report is prefixed MIXTUREPBM with 12-significant-digit numbers; the
// checker reconstructs the same trajectory from the case constants alone.
#include <iomanip>
#include <sstream>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/fvm/MixturePbm.h"

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
        std::cerr << "usage: chflow-3d-mixture-pbm <case.dict>\n";
        return 2;
    }
    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const MixturePbmConfig cfg = readMixturePbmConfig(caseDict);
        if (!cfg.enabled)
            cfFatal("mixturePbmApp: needs a mixturePbm{} declaration");

        const MixturePbmProps p = mixtureProperties(cfg);
        log::info("MIXTUREPBM config alpha0 ", num(cfg.alpha0), " d0 ", num(cfg.d0),
                  " rhoL ", num(cfg.rhoL), " rhoG ", num(cfg.rhoG),
                  " slip ", num(cfg.slip), " j ", num(cfg.j),
                  " J ", num(cfg.J), " breakage ", num(cfg.breakage),
                  " nFrag ", num(cfg.nFrag), " dt ", num(cfg.dt),
                  " steps ", cfg.steps, " reportEvery ", cfg.reportEvery,
                  " vPart ", num(p.vPart), " k ", num(p.k),
                  " m0 ", num(p.m0_0), " m1 ", num(p.m1_0),
                  " m2 ", num(p.m2_0), " m3 ", num(p.m3_0));

        const MixturePbmResult res = runMixturePbm(cfg);
        if (!res.valid)
        {
            log::error("mixturePbmApp: integration failed");
            return 1;
        }

        for (const MixturePbmStep& s : res.history)
        {
            log::info("MIXTUREPBM step ", s.step, " t ", num(s.t),
                      " M0 ", num(s.m[0]), " M1 ", num(s.m[1]),
                      " M2 ", num(s.m[2]), " M3 ", num(s.m[3]),
                      " alpha ", num(s.alpha), " uG ", num(s.uG),
                      " uL ", num(s.uL));
        }

        const MixturePbmStep& last = res.history.back();
        log::info("MIXTUREPBM summary records ", res.history.size(),
                  " finalM0 ", num(last.m[0]), " finalM1 ", num(last.m[1]),
                  " finalM2 ", num(last.m[2]), " finalM3 ", num(last.m[3]),
                  " alpha ", num(last.alpha));
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("mixturePbmApp: ", e.what());
        return 1;
    }
}
