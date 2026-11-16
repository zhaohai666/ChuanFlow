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
                if (f.regimeOutlet >= 0)
                {
                    line += " oRegime "
                          + std::string(regimeName[f.regimeOutlet])
                          + " lambda " + num(f.lambdaOutlet)
                          + " Nfr " + num(f.nfrOutlet)
                          + " L1 " + num(f.patternL1)
                          + " L2 " + num(f.patternL2)
                          + " L3 " + num(f.patternL3)
                          + " L4 " + num(f.patternL4);
                    // 5.88 (D4 increment 2): the override's footprint. Printed
                    // only when a forced pattern actually disagreed with the
                    // map at the outlet segment, so an auto (or agreeing)
                    // declaration leaves the line byte-identical to 5.85.
                    if (f.regimeFromMapOutlet >= 0
                        && f.regimeFromMapOutlet != f.regimeOutlet)
                    {
                        line += " mapRegime "
                              + std::string(regimeName[f.regimeFromMapOutlet])
                              + " FORCED";
                    }
                }
                if (f.choked) line += " CHOKED";
            }
            else if (f.valveDisposition >= 0)
            {
                // 5.89 (D1 increment 3): a gas-lift valve has no wall law to
                // report a Reynolds number and a friction factor from. What it
                // does have is a force balance and a nozzle, so the line says
                // which side of its own thresholds the injection pressure came
                // in on, the two thresholds, the band between them and how far
                // the tubing pressure pulled the opening pressure off its
                // test-rack value, and the mass rate the port passed.
                static const char* dispositionName[] = {"CLOSED", "OPEN",
                                                        "IN-SPREAD"};
                line += " valve " + std::string(dispositionName[f.valveDisposition])
                      + " pvo " + num(f.valveOpeningPressure)
                      + " pvc " + num(f.valveClosingPressure)
                      + " spread " + num(f.valveSpread)
                      + " ppe " + num(f.valveProductionEffect)
                      + " tro " + num(f.valveTestRackOpening)
                      + " m " + num(f.massRate)
                      + (f.valveChoked ? " SONIC" : "");
            }
            else
            {
                line += " Re " + num(f.reynolds) + " f "
                      + num(f.frictionFactor);
            }
            line += " dp " + num(f.dpTotal);
            if (f.dpValve != 0) line += " dpValve " + num(f.dpValve);
            if (f.dpPump != 0) line += " dpPump " + num(f.dpPump);
            if (f.transition.band >= 0)
            {
                line += " He " + num(f.transition.hedstrom)
                      + " ReB " + num(f.transition.rePlastic)
                      + " ReBt " + num(f.transition.reTransition)
                      + " Qt " + num(f.transition.rateAtTransition)
                      + " dpQt " + num(f.transition.dpAtTransition)
                      + " Vrelt " + num(f.transition.relTransitionSpeed)
                      + " band "
                      + num(static_cast<double>(f.transition.band))
                      + (f.transition.pastLaminarLimit ? " pastTransition"
                                                       : " belowTransition");
            }
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
