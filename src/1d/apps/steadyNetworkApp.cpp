#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "cf/1d/Network.h"
#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"

using namespace cf;

namespace
{
// log::info() stringifies with ostringstream defaults (6 digits); results
// are checked against closed-form values, so emit 12 significant digits.
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
        std::cerr << "usage: chflow-1d-steady <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const OneDNetwork net = loadNetworkFromDict(caseDict);
        const SteadyControls ctl = loadSteadyControls(caseDict);

        log::info("network: ", net.nNodes(), " nodes, ", net.pipes.size(), " pipes");

        const SteadyNetworkResult res = solveSteadyNetwork(net, ctl);

        for (std::size_t i = 0; i < net.nNodes(); ++i)
        {
            log::info("node ", net.nodeNames[i], " p ", num(res.nodePressure[i]));
        }
        for (const PipeFlowResult& f : res.pipeFlows)
        {
            std::string line = "pipe " + f.name + " Q " + num(f.flowRate)
                             + " v " + num(f.velocity);
            if (f.regime >= 0)
            {
                static const char* regimeName[] = {"segregated", "transition",
                                                   "intermittent", "distributed"};
                line += " HL " + num(f.holdup) + " regime "
                      + std::string(regimeName[f.regime])
                      + " m " + num(f.massRate);
                if (f.choked) line += " CHOKED";
            }
            else
            {
                line += " Re " + num(f.reynolds) + " f "
                      + num(f.frictionFactor);
            }
            line += " dp " + num(f.dpTotal);
            if (f.dpValve != 0) line += " dpValve " + num(f.dpValve);
            if (f.dpPump != 0) line += " dpPump " + num(f.dpPump);
            log::info(line);
        }

        log::info("steady solved: iterations=", res.iterations,
                  " residual=", num(res.relativeResidual),
                  res.converged ? " (converged)" : " (NOT converged)");

        return res.converged ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        log::error("1d-steady: ", e.what());
        return 1;
    }
}
