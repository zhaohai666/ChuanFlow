#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/pvt/PvtCase.h"
#include "cf/pvt/PvtTable.h"

using namespace cf;
using namespace cf::pvt;

namespace
{
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}

void reportOne(const FlashCase& c, double p)
{
    const FlashResult r = tpFlash(c.feed, c.temperature, p, c.controls);
    log::info("pvt flash T ", num(c.temperature), " p ", num(p),
              " beta ", num(r.vaporFraction), " single ", r.singlePhase ? 1 : 0,
              " iterations ", r.iterations,
              " converged ", r.converged ? 1 : 0, " residual ", num(r.residual));
    log::info("pvt liquid rho ", num(r.liquidMassDensity), " M ", num(r.liquidMolarMass));
    log::info("pvt vapour rho ", num(r.vapourMassDensity), " M ", num(r.vapourMolarMass));
    for (std::size_t i = 0; i < c.feed.size(); ++i)
    {
        log::info("pvt x ", c.feed.components()[i].name, " ", num(r.liquidComposition[i]));
        log::info("pvt y ", c.feed.components()[i].name, " ", num(r.vapourComposition[i]));
    }
    if (!r.converged) log::warn("flash: ", r.message);
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-pvt <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const FlashCase c = loadFlashCase(caseDict);

        log::info("pvt: ", c.feed.size(), " components, feed ", c.feed.describe(),
                  " molarMass ", num(c.feed.molarMass()));

        for (const double p : c.pressures) reportOne(c, p);

        if (c.computeBubblePoint)
        {
            const double pb = bubblePressure(c.feed, c.temperature);
            log::info("pvt bubble T ", num(c.temperature), " p ", num(pb));
            reportOne(c, pb);
        }

        // A coarse table shows interpolation and EOS agree; the regression
        // script samples it back through the same keys.
        if (caseDict.has("table"))
        {
            const Dictionary& tbl = caseDict.getDict("table");
            const std::vector<double> pGrid = tbl.getScalarList("pressures");
            const std::vector<double> tGrid = tbl.getScalarList("temperatures");
            const std::string prop = tbl.getWordOr("property", "liquidDensity");
            std::optional<BilinearTable> table;
            if (prop == "liquidDensity") table = liquidDensityTable(c.feed, pGrid, tGrid);
            else if (prop == "vapourDensity") table = vapourDensityTable(c.feed, pGrid, tGrid);
            else if (prop == "vapourFraction") table = vapourFractionTable(c.feed, pGrid, tGrid);
            else cfFatal("pvt: unknown table property '" + prop + "'");

            for (const double p : pGrid)
                for (const double t : tGrid)
                    log::info("pvt table ", prop, " T ", num(t), " p ", num(p),
                              " value ", num(table->sample(p, t)));
        }
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error(e.what());
        return 1;
    }
}
