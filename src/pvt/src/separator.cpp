#include "cf/pvt/Separator.h"

#include "cf/common/Error.h"

namespace cf::pvt
{
namespace
{
constexpr double kRhoWaterStd = 999.014; // kg/m3 at 60 degF, as in BlackOil.h

/// Gas density at standard conditions, recovered from the existing mass
/// invariant instead of duplicating the unit constants: for waterCut = 0,
/// blackOilMassPerOilSm3 = oilDensityStd + gorStd * rhoGasStd, so the
/// difference between gorStd 1 and gorStd 0 is exactly rhoGasStd.
double gasDensityStd(const BlackOilFluid& fluid)
{
    BlackOilFluid f0 = fluid;
    f0.gorStd = 0.0;
    f0.waterCut = 0.0;
    BlackOilFluid f1 = f0;
    f1.gorStd = 1.0;
    return blackOilMassPerOilSm3(f1) - blackOilMassPerOilSm3(f0);
}
} // namespace

SeparatorResult separateBlackOil(const BlackOilFluid& fluid,
                                 const std::vector<double>& pressures,
                                 const std::vector<double>& temperatures)
{
    if (pressures.empty()) cfFatal("separator: needs at least one stage pressure");
    if (temperatures.size() != pressures.size())
        cfFatal("separator: one temperature per stage pressure is required");
    if (!(fluid.gasGravity > 0.0)) cfFatal("separator: gasGravity must be positive");
    if (!(fluid.oilDensityStd > 0.0)) cfFatal("separator: oilDensity must be positive");
    if (!(fluid.gorStd >= 0.0)) cfFatal("separator: gor must be non-negative");
    if (!(fluid.waterCut >= 0.0 && fluid.waterCut < 1.0))
        cfFatal("separator: waterCut must be in [0,1)");

    const double wor = fluid.waterCut > 0.0
                         ? fluid.waterCut / (1.0 - fluid.waterCut)
                         : 0.0;
    const double rhoGasStd = gasDensityStd(fluid);
    const double waterMass = wor * fluid.waterDensity;

    SeparatorResult out;
    out.stockTankApi = 141.5 / (fluid.oilDensityStd / kRhoWaterStd) - 131.5;
    out.gasMassPerSm3 = rhoGasStd;
    out.massIn = blackOilMassPerOilSm3(fluid);

    double available = fluid.gorStd;
    double cumGasMass = 0.0;
    double totalGas = 0.0;
    for (std::size_t i = 0; i < pressures.size(); ++i)
    {
        const double p = pressures[i];
        const double T = temperatures[i];
        if (!(p > 0.0)) cfFatal("separator: stage pressures must be positive");
        if (!(T > 0.0)) cfFatal("separator: stage temperatures must be positive");

        // Re-flash the liquid that actually arrived: its GOR is the gas still
        // in solution, not the original producing GOR. This drops the bubble
        // point stage by stage and keeps the cascade stateless.
        BlackOilFluid f = fluid;
        f.gorStd = available;
        f.temperature = T;
        const BlackOilState st = evaluateBlackOil(f, p);

        SeparatorStage sg;
        sg.index = static_cast<int>(i);
        sg.pressure = p;
        sg.temperature = T;
        sg.availableGas = available;
        sg.rs = st.rs;
        sg.gasLiberated = available - st.rs;
        sg.bo = st.bo;
        sg.bg = st.bg;
        sg.z = st.z;
        sg.rhoOil = st.rhoOil;
        sg.gasMass = sg.gasLiberated * rhoGasStd;
        sg.liquidMass = fluid.oilDensityStd + st.rs * rhoGasStd + waterMass;
        cumGasMass += sg.gasMass;
        sg.cumGasMass = cumGasMass;
        sg.liveApi = 141.5 / (st.rhoOil / kRhoWaterStd) - 131.5;

        totalGas += sg.gasLiberated;
        available = st.rs; // liquid carries only the remaining dissolved gas
        out.stages.push_back(sg);
    }

    out.totalGasSm3 = totalGas;
    out.massOut = out.stages.back().liquidMass + cumGasMass;
    out.massImbalance = out.massOut - out.massIn;
    return out;
}

SeparatorResult separateBlackOil(const BlackOilFluid& fluid,
                                 const std::vector<double>& pressures,
                                 double temperature)
{
    return separateBlackOil(fluid, pressures,
                            std::vector<double>(pressures.size(), temperature));
}

} // namespace cf::pvt
