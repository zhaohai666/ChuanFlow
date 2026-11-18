// M6 industrial model -- stateless N-stage black-oil separator cascade.
//
// Every stage re-flashes the liquid that actually arrived (its GOR is the gas
// still in solution) through the existing evaluateBlackOil, so the bubble
// point falls as gas is removed. The stage gas is the drop in solution GOR;
// mass is closed by liquidMass_i + cumGasMass_i == blackOilMassPerOilSm3.
//
// case grammar:
//   separator
//   {
//       pressure (2e6 5e5 1.01325e5);
//       temperature 340;                 // isothermal, or
//       // temperatures (340 330 315);   // one per stage (optional override)
//       fluid { blackoil { gasGravity 0.75; oilDensity 850;
//                          gor 180; waterCut 0.05; } }
//   }
//
// usage: chflow-separator <case.dict>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/pvt/Separator.h"

using namespace cf;
using namespace cf::pvt;

namespace
{
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-separator <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary dict = Dictionary::parseFile(argv[1]);
        if (!dict.has("separator"))
            cfFatal("separator: case needs a separator{...} block");
        const Dictionary& sep = dict.getDict("separator");

        if (!sep.has("pressure")) cfFatal("separator: missing pressure (...) list");
        const std::vector<double> pressures = sep.getScalarList("pressure");
        if (pressures.empty()) cfFatal("separator: pressure (...) list is empty");

        std::vector<double> temperatures;
        if (sep.has("temperatures"))
        {
            temperatures = sep.getScalarList("temperatures");
            if (temperatures.size() != pressures.size())
                cfFatal("separator: temperatures must give one value per stage");
        }
        else
        {
            if (!sep.has("temperature")) cfFatal("separator: missing temperature");
            temperatures.assign(pressures.size(), sep.getScalar("temperature"));
        }

        if (!sep.has("fluid") || !sep.getDict("fluid").has("blackoil"))
            cfFatal("separator: case needs a fluid{ blackoil{...} } block");
        const Dictionary& bd = sep.getDict("fluid").getDict("blackoil");
        BlackOilFluid fluid;
        fluid.gasGravity = bd.getScalarOr("gasGravity", fluid.gasGravity);
        fluid.oilDensityStd = bd.getScalarOr("oilDensity", fluid.oilDensityStd);
        fluid.gorStd = bd.getScalarOr("gor", fluid.gorStd);
        fluid.waterCut = bd.getScalarOr("waterCut", fluid.waterCut);
        fluid.waterDensity = bd.getScalarOr("waterDensity", fluid.waterDensity);
        fluid.temperature = temperatures.front();

        const SeparatorResult r = separateBlackOil(fluid, pressures, temperatures);

        log::info("separator: gg ", num(fluid.gasGravity), " rhoStd ",
                  num(fluid.oilDensityStd), " gor ", num(fluid.gorStd),
                  " waterCut ", num(fluid.waterCut), " stages ", r.stages.size());

        for (const SeparatorStage& s : r.stages)
        {
            log::info("separator stage ", s.index,
                      " p ", num(s.pressure),
                      " T ", num(s.temperature),
                      " available ", num(s.availableGas),
                      " rs ", num(s.rs),
                      " gasSm3 ", num(s.gasLiberated),
                      " gasMass ", num(s.gasMass),
                      " liquidMass ", num(s.liquidMass),
                      " cumGasMass ", num(s.cumGasMass),
                      " bo ", num(s.bo),
                      " bg ", num(s.bg),
                      " z ", num(s.z),
                      " liveApi ", num(s.liveApi));
        }

        log::info("separator total gasSm3 ", num(r.totalGasSm3),
                  " gasMass ", num(r.gasMassPerSm3 * r.totalGasSm3),
                  " massIn ", num(r.massIn),
                  " massOut ", num(r.massOut),
                  " imbalance ", num(r.massImbalance),
                  " stockTankApi ", num(r.stockTankApi));
        if (std::abs(r.massImbalance) > 1e-9 * std::abs(r.massIn))
            cfFatal("separator: mass balance did not close");
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error(e.what());
        return 1;
    }
}
