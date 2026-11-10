#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "cf/1d/Transient.h"
#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"

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
        std::cerr << "usage: chflow-1d-transient <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const TransientCase c = loadTransientFromDict(caseDict);
        const TransientTimeControls tc = loadTransientControls(caseDict);

        log::info("transient: ", c.nodes.size(), " nodes, ", c.pipes.size(),
                  " pipes, dt=", num(tc.dt), " endTime=", num(tc.endTime));

        TransientNetwork sim(c);

        // Say which wall law each momentum row is reading, and say it only where
        // one is declared: an undeclared case prints exactly what it printed
        // before, so the 5.48 claim that declaring n = 1, K = mu changes nothing
        // is testable on the output itself, not only on the numbers.
        for (std::size_t k = 0; k < sim.pipeNames().size(); ++k)
        {
            const label i = static_cast<label>(k);
            const double n = sim.pipeFlowIndex(i);
            if (n > 0)
            {
                log::info("pipe ", sim.pipeNames()[k], " wall-law flowIndex ",
                          num(n), " consistency ", num(sim.pipeConsistencyOf(i)));
            }
        }

        const int report = tc.reportEvery > 0 ? tc.reportEvery : 1;
        const int nSteps = static_cast<int>(std::round(tc.endTime / tc.dt));

        auto emit = [&](const TransientNetwork& s)
        {
            log::info("t ", num(s.time()));
            const std::vector<double> a = s.liquidFractions();
            const std::vector<double> h = s.hydrateFractions();
            for (std::size_t i = 0; i < s.nodeNames().size(); ++i)
            {
                log::info("node ", s.nodeNames()[i], " p ", num(s.pressure()[i]),
                          " alpha ", num(a[i]), " hydrate ", num(h[i]));
            }
            const std::vector<double> q = s.pipeFlows();
            for (std::size_t k = 0; k < q.size(); ++k)
            {
                log::info("pipe ", s.pipeNames()[k], " Q ", num(q[k]));
            }
        };

        emit(sim);
        for (int i = 1; i <= nSteps; ++i)
        {
            const double res = sim.step(tc.dt);
            if (i % report == 0 || i == nSteps) emit(sim);
            if (!std::isfinite(res))
            {
                log::error("transient: non-finite Newton residual at step ", i);
                return 1;
            }
        }
        log::info("transient solved: steps=", nSteps);
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("1d-transient: ", e.what());
        return 1;
    }
}
