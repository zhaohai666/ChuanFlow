#include "cf/1d/Separator.h"

#include "cf/common/Error.h"

namespace cf
{

SeparatorSolution separatorSolve(const pvt::BlackOilFluid& fluid,
                                 double feedRateSm3PerS,
                                 double pressurePa,
                                 double temperatureK,
                                 double efficiency,
                                 double dPPa)
{
    if (!(feedRateSm3PerS > 0.0))
        cfFatal("separator: the feed rate must be positive: a stock-tank oil "
                "rate of zero has no outlet stream to report");
    if (!(pressurePa > 0.0))
        cfFatal("separator: the stage pressure is an absolute pressure and "
                "must be positive");
    if (!(temperatureK > 0.0))
        cfFatal("separator: the stage temperature is kelvin and must be "
                "positive");
    if (!(efficiency >= 0.0 && efficiency <= 1.0))
        cfFatal("separator: the liquid recovery efficiency must lie in "
                "[0,1]: it is the fraction of the condensables the liquid "
                "outlet actually carries");
    if (!(dPPa >= 0.0))
        cfFatal("separator: dP is a non-negative pressure drop across the "
                "internals: the outlets cannot sit above the stage");
    if (!(fluid.gasGravity > 0.0))
        cfFatal("separator: gasGravity must be positive");
    if (!(fluid.oilDensityStd > 0.0))
        cfFatal("separator: oilDensity must be positive");
    if (!(fluid.gorStd >= 0.0))
        cfFatal("separator: gor must be non-negative");
    if (!(fluid.waterCut >= 0.0 && fluid.waterCut < 1.0))
        cfFatal("separator: waterCut must be in [0,1)");

    // The stage flashes the arriving stream at its own temperature: the
    // black-oil layer reads the fluid's isothermal T, so the feed is copied
    // with the stage value before evaluateBlackOil is asked.
    pvt::BlackOilFluid stage = fluid;
    stage.temperature = temperatureK;
    const pvt::BlackOilState st = pvt::evaluateBlackOil(stage, pressurePa);

    // The standard gas density, recovered from the existing mass invariant
    // exactly as the pvt cascade does it: for waterCut = 0,
    // blackOilMassPerOilSm3 = oilDensityStd + gorStd * rhoGasStd, so the
    // difference between gorStd 1 and gorStd 0 is rhoGasStd.
    pvt::BlackOilFluid f0 = fluid;
    f0.gorStd = 0.0;
    f0.waterCut = 0.0;
    pvt::BlackOilFluid f1 = f0;
    f1.gorStd = 1.0;
    const double rhoGasStd = pvt::blackOilMassPerOilSm3(f1)
                           - pvt::blackOilMassPerOilSm3(f0);

    const double wor = fluid.waterCut > 0.0
                         ? fluid.waterCut / (1.0 - fluid.waterCut)
                         : 0.0;
    const double waterMass = wor * fluid.waterDensity;
    const double feedPerSm3 = fluid.oilDensityStd
                            + fluid.gorStd * rhoGasStd + waterMass;
    const double liquidPerSm3 = fluid.oilDensityStd
                              + st.rs * rhoGasStd + waterMass;
    const double gasPerSm3 = (fluid.gorStd - st.rs) * rhoGasStd;

    SeparatorSolution s;
    s.pressure = pressurePa;
    s.temperature = temperatureK;
    s.efficiency = efficiency;
    s.outletPressure = pressurePa - dPPa;
    s.rs = st.rs;
    s.gasLiberated = fluid.gorStd - st.rs;
    s.bo = st.bo;
    s.bg = st.bg;
    s.z = st.z;
    s.rhoGasStd = rhoGasStd;
    s.feedMassRate = feedRateSm3PerS * feedPerSm3;
    s.liquidMassRateIdeal = feedRateSm3PerS * liquidPerSm3;
    s.entrainedMassRate = (1.0 - efficiency) * s.liquidMassRateIdeal;
    s.liquidMassRate = efficiency * s.liquidMassRateIdeal;
    s.gasMassRate = feedRateSm3PerS * gasPerSm3 + s.entrainedMassRate;
    s.closure = s.liquidMassRate + s.gasMassRate - s.feedMassRate;
    return s;
}

} // namespace cf
