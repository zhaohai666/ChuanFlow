#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "cf/1d/TwoFluid.h"
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
        std::cerr << "usage: chflow-1d-twofluid <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const TwoFluidCase c = loadTwoFluidCase(caseDict);

        const double dt = caseDict.getDict("time").getScalar("dt");
        const double endTime = caseDict.getDict("time").getScalar("endTime");
        const int reportEvery = static_cast<int>(
            caseDict.getDict("time").getScalarOr("reportEvery", 1.0));

        log::info("twoFluid: ", c.geo.sections, " sections, dt=", num(dt),
                  " endTime=", num(endTime));

        TwoFluidPipe sim(c);

        auto emit = [&](const TwoFluidPipe& s)
        {
            log::info("t ", num(s.time()), " vm ", num(s.mixtureVelocity()),
                      " dp ", num(s.pressureDrop()));
            const std::vector<double> p = s.nodePressure();
            for (std::size_t j = 0; j < p.size(); ++j)
                log::info("node ", j, " p ", num(p[j]));
            const std::vector<double> a = s.cellAlpha();
            const std::vector<double> w = s.cellWater();
            const std::vector<double> vsl = s.cellVsl();
            const std::vector<double> vsg = s.cellVsg();
            const std::vector<double> T = s.cellTemperature();
            for (int i = 0; i < s.sections(); ++i)
            {
                log::info("cell ", i, " alpha ", num(a[i]), " water ",
                          num(w[i]), " vsl ", num(vsl[i]), " vsg ",
                          num(vsg[i]),
                          T.empty() ? "" : (" T " + num(T[i])));
            }
            if (!T.empty())
                log::info("thermal: qWall ", num(s.wallHeatRate()),
                          " qInlet ", num(s.inletHeatRate()));
        };

        emit(sim);
        const int report = reportEvery > 0 ? reportEvery : 1;
        const int nSteps = static_cast<int>(std::round(endTime / dt));
        for (int i = 1; i <= nSteps; ++i)
        {
            sim.step(dt);
            if (i % report == 0 || i == nSteps) emit(sim);
        }
        log::info("twoFluid solved: steps=", nSteps);
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("1d-twofluid: ", e.what());
        return 1;
    }
}
