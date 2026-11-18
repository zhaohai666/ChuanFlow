// chflow-3d-qmom: two-node QMOM for constant-kernel coagulation (D8 increment
// 6). No flow field: this is the spatially-uniform population balance, the
// moment-closure kernel the DPM/species source framework would call. The
// checker reconstructs the closed-form moment trajectory and the Gaussian
// quadrature identities from the case constants alone.
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/fvm/Qmom.h"

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
        std::cerr << "usage: chflow-3d-qmom <case.dict>\n";
        return 2;
    }
    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const QmomConfig cfg = readQmomConfig(caseDict);
        if (!cfg.enabled) cfFatal("qmomApp: needs a qmom{} declaration");

        log::info("qmom model: beta ", num(cfg.beta), " dt ", num(cfg.dt),
                  " steps ", cfg.steps, " M0 ", num(cfg.m0), " M1 ", num(cfg.m1),
                  " M2 ", num(cfg.m2), " M3 ", num(cfg.m3));

        const QmomResult res = runQmom(cfg);
        if (!res.valid)
        {
            log::error("qmomApp: integration failed");
            return 1;
        }

        for (const QmomStep& s : res.history)
        {
            log::info("qmom step ", s.step, ": t ", num(s.t),
                      " M ", num(s.m[0]), " ", num(s.m[1]), " ",
                      num(s.m[2]), " ", num(s.m[3]),
                      " nodes ", num(s.v[0]), " ", num(s.v[1]),
                      " weights ", num(s.w[0]), " ", num(s.w[1]),
                      " massErr ", num(s.massError));
            // Self-evidence of the quadrature identity; the checker recomputes
            // this independently so it never trusts this line.
            for (int k = 0; k < 4; ++k)
            {
                const scalar rebuilt = s.w[0] * std::pow(s.v[0], k)
                                       + s.w[1] * std::pow(s.v[1], k);
                log::info("qmom moment ", s.step, " ", k, ": ", num(rebuilt));
            }
        }
        log::info("qmom summary: records ", res.history.size(), " finalMassErr ",
                  num(res.history.back().massError));
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("qmomApp: ", e.what());
        return 1;
    }
}
