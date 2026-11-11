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
            // 5.49: the drop the static wall holds against, so the case file's
            // own drive can be compared with it on the output. The only pipe
            // that never says this is one with no threshold, so an undeclared
            // case still prints exactly what it printed before.
            const double thr = sim.pipeYieldThreshold(i);
            if (thr > 0)
            {
                log::info("pipe ", sim.pipeNames()[k], " static-hold threshold ",
                          num(thr));
            }
        }

        const int report = tc.reportEvery > 0 ? tc.reportEvery : 1;
        // 5.66 / C4: same bookkeeping as runTransient -- ceil(endTime/dt) with
        // a near-multiple tolerance, last step clamped to the remainder, and
        // the final state always emitted. A case whose endTime IS a multiple
        // of dt prints exactly what it printed before.
        if (!(tc.dt > 0) || !(tc.endTime > 0))
        {
            log::error("transient: dt and endTime must both be positive");
            return 1;
        }
        // 5.68 / D5: cfl > 0 switches the loop to CFL-adaptive stepping --
        // the next dt is the Courant rule over the whole network, first step
        // is the case's own dt. Say so on the output so the reader sees which
        // march the t column is a record of; an undeclared case prints
        // exactly what it printed before.
        const bool adaptive = tc.cfl > 0;
        if (adaptive)
        {
            if (!(tc.dtMin > 0))
            {
                log::error("transient: adaptive stepping needs dtMin > 0");
                return 1;
            }
            if (tc.dtMax < tc.dtMin)
            {
                log::error("transient: dtMax must not be below dtMin");
                return 1;
            }
            log::info("transient adaptive: cfl ", num(tc.cfl), " dtMin ",
                      num(tc.dtMin), " dtMax ", num(tc.dtMax));
        }
        const int nSteps =
            static_cast<int>(std::ceil(tc.endTime / tc.dt - 1e-9));

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
                // The momentum-flux account the row charges, printed only for
                // a pipe whose case declares betaIn/betaOut (5.70): the value
                // is convective * q|q|, the same account a coupled app prints
                // as dpConv (5.42). A case that asks for nothing prints
                // nothing, so untouched traces are unchanged.
                if (s.pipeBetaInOf(static_cast<label>(k))
                        != s.pipeBetaOutOf(static_cast<label>(k)))
                {
                    const double cc = s.pipeConvectiveOf(static_cast<label>(k))
                                    * q[k] * std::abs(q[k]);
                    log::info("pipe ", s.pipeNames()[k], " dpConv ", num(cc));
                }
                // The latch the step that produced this sample started with: 1
                // means this rate is the row's own `Q = 0`, not the value a wall
                // law was solved for. At t 0 no step has been solved yet, so the
                // first sample reports 0 whether or not the drive can be held.
                if (s.pipeYieldThreshold(static_cast<label>(k)) > 0)
                {
                    log::info("pipe ", s.pipeNames()[k], " static hold ",
                              s.pipeHeld(static_cast<label>(k)) ? 1 : 0);
                }
            }
        };

        emit(sim);
        if (adaptive)
        {
            // 5.68 / D5: the adaptive march. t accumulates -- with varying dt
            // there is no index to position one -- and the last step is
            // clamped to endTime - t, with t snapped to endTime exactly when
            // the clamp fires. Every state is emitted at its own step boundary
            // per reportEvery, and the final state is always emitted.
            double dtCur = tc.dt;
            double t = 0.0;
            long long i = 0;
            while (t < tc.endTime - 1e-12)
            {
                const double rest = tc.endTime - t;
                const bool clampFires = rest <= dtCur;
                const double h = clampFires ? rest : dtCur;
                const double res = sim.step(h);
                ++i;
                t = clampFires ? tc.endTime : t + h;
                if (i % report == 0 || clampFires) emit(sim);
                if (!std::isfinite(res))
                {
                    log::error("transient: non-finite Newton residual at step ",
                               i);
                    return 1;
                }
                dtCur = sim.cflStep(tc.cfl, tc.dtMin, tc.dtMax);
            }
            log::info("transient solved: steps=", i);
            return 0;
        }
        for (int i = 1; i <= nSteps; ++i)
        {
            const double tPrev = std::min((i - 1) * tc.dt, tc.endTime);
            const double h = std::min(tc.dt, tc.endTime - tPrev);
            const double res = sim.step(h);
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
