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

        const double dt = caseDict.getDict("time").getScalar("dt");
        const double endTime = caseDict.getDict("time").getScalar("endTime");
        const int reportEvery = static_cast<int>(
            caseDict.getDict("time").getScalarOr("reportEvery", 1.0));
        const int report = reportEvery > 0 ? reportEvery : 1;
        const int nSteps = static_cast<int>(std::round(endTime / dt));

        if (caseDict.has("network"))
        {
            const TwoFluidNetworkCase nc = loadTwoFluidNetworkCase(caseDict);
            int sections = 0;
            for (const TwoFluidGeo& g : nc.pipes) sections += g.sections;
            log::info("twoFluid network: ", nc.pipes.size(), " pipes ",
                      sections, " sections, ", nc.junctions.size(),
                      " junctions, dt=", num(dt), " endTime=", num(endTime));

            TwoFluidNetwork net(nc);
            auto emit = [&](const TwoFluidNetwork& n)
            {
                log::info("t ", num(n.time()), " net");
                for (const std::string& j : n.junctionNames())
                    log::info("junction ", j, " p ",
                              num(n.junctionPressure(j)), " net ",
                              num(n.junctionNetRate(j)));
                for (const std::string& pn : n.pipeNames())
                {
                    const TwoFluidPipe& s = n.pipe(pn);
                    log::info("pipe ", pn, " vm ", num(s.mixtureVelocity()),
                              " dp ", num(s.pressureDrop()));
                    const std::vector<double> p = s.nodePressure();
                    for (std::size_t q = 0; q < p.size(); ++q)
                        log::info("pnode ", pn, " ", q, " p ", num(p[q]));
                    const std::vector<double> a = s.cellAlpha();
                    const std::vector<double> w = s.cellWater();
                    const std::vector<double> vsl = s.cellVsl();
                    const std::vector<double> vsg = s.cellVsg();
                    const std::vector<double> T = s.cellTemperature();
                    const std::vector<double> ep
                        = nc.perPipe.fluid.gasR > 0 ? s.cellPressure()
                                                    : std::vector<double>();
                    const std::vector<double> hy
                        = nc.perPipe.hydrateOn ? s.cellHydrate()
                                               : std::vector<double>();
                    const std::vector<double> ent
                        = nc.perPipe.entrainmentOn ? s.cellEntrained()
                                                   : std::vector<double>();
                    for (int i = 0; i < s.sections(); ++i)
                        log::info("pcell ", pn, " ", i, " alpha ",
                                  num(a[i]), " water ", num(w[i]), " vsl ",
                                  num(vsl[i]), " vsg ", num(vsg[i]),
                                  T.empty() ? "" : (" T " + num(T[i])),
                                  ep.empty() ? "" : (" p " + num(ep[i])),
                                  hy.empty() ? "" : (" h " + num(hy[i])),
                                  ent.empty() ? ""
                                              : (" entr " + num(ent[i])));
                }
            };

            emit(net);
            for (int i = 1; i <= nSteps; ++i)
            {
                net.step(dt);
                if (i % report == 0 || i == nSteps) emit(net);
            }
            log::info("twoFluid network solved: steps=", nSteps);
            return 0;
        }

        const TwoFluidCase c = loadTwoFluidCase(caseDict);

        int sections = 0;
        for (const TwoFluidGeo& g : c.pipes) sections += g.sections;
        log::info("twoFluid: ", sections, " sections, dt=", num(dt),
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
            const std::vector<double> sl
                = c.closure.model == 2 || c.closure.model == 3
                      ? s.cellSlip()
                      : std::vector<double>();
            const std::vector<double> T = s.cellTemperature();
            const std::vector<double> ep
                = c.fluid.gasR > 0 ? s.cellPressure() : std::vector<double>();
            const std::vector<double> hy
                = c.hydrateOn ? s.cellHydrate() : std::vector<double>();
            const std::vector<double> ent
                = c.entrainmentOn ? s.cellEntrained()
                                  : std::vector<double>();
            for (int i = 0; i < s.sections(); ++i)
            {
                log::info("cell ", i, " alpha ", num(a[i]), " water ",
                          num(w[i]), " vsl ", num(vsl[i]), " vsg ",
                          num(vsg[i]),
                          T.empty() ? "" : (" T " + num(T[i])),
                          sl.empty() ? "" : (" s " + num(sl[i])),
                          ep.empty() ? "" : (" p " + num(ep[i])),
                          hy.empty() ? "" : (" h " + num(hy[i])),
                          ent.empty() ? "" : (" entr " + num(ent[i])));
            }
            if (c.slugOn)
            {
                const std::vector<TwoFluidSlugBody> bodies = s.slugBodies();
                log::info("slug count ", bodies.size());
                for (std::size_t k = 0; k < bodies.size(); ++k)
                    log::info("slugbody ", k, " rear ",
                              num(bodies[k].rear), " front ",
                              num(bodies[k].front), " length ",
                              num(bodies[k].length), " inv ",
                              num(bodies[k].inventory));
            }
            if (!T.empty())
                log::info("thermal: qWall ", num(s.wallHeatRate()),
                          " qInlet ", num(s.inletHeatRate()));
        };

        emit(sim);
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
