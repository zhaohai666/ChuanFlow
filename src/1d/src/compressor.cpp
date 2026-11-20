#include "cf/1d/Compressor.h"

#include <cmath>

#include "cf/common/Error.h"

namespace cf
{

CompressorSolution compressorSolve(const Compressor& c,
                                   double suctionPressurePa,
                                   double suctionTemperatureK,
                                   double rateSm3PerS)
{
    const bool curve = c.pressureRatio > 0.0;
    const bool head = c.polytropicHead > 0.0;
    if (curve && head)
        cfFatal("compressor: pressureRatio and polytropicHead together state "
                "two different things about the same ratio: price the machine "
                "off its curve or off its head, not both");
    if (!curve && !head)
        cfFatal("compressor: neither pressureRatio nor polytropicHead is "
                "declared: a machine must say how it raises the pressure");
    if (!(suctionPressurePa > 0.0))
        cfFatal("compressor: the suction pressure is an absolute pressure and "
                "must be positive");
    if (!(suctionTemperatureK > 0.0))
        cfFatal("compressor: the suction temperature is kelvin and must be "
                "positive");
    if (!(rateSm3PerS > 0.0))
        cfFatal("compressor: the rate must be positive: at zero flow the "
                "machine still turns (the curve reads shutoffRatio), but it "
                "has no stream to report a power for");
    if (!(c.polytropicEfficiency > 0.0 && c.polytropicEfficiency <= 1.0))
        cfFatal("compressor: the polytropic efficiency must lie in (0,1]: it "
                "divides the gas power into shaft power and sets the "
                "discharge-temperature exponent");
    if (!(c.kappa > 1.0))
        cfFatal("compressor: kappa (cp/cv) must exceed 1: the exponent "
                "(k-1)/k is what turns a pressure ratio into a head");
    if (!(c.specificGasConstant > 0.0))
        cfFatal("compressor: the specific gas constant must be positive: the "
                "head is priced on R_spec T1");
    if (!(c.compressibility > 0.0))
        cfFatal("compressor: the compressibility Z must be positive: it is a "
                "declared constant of the suction state");

    const double ex = (c.kappa - 1.0) / c.kappa;

    CompressorSolution s;
    s.pressureIn = suctionPressurePa;
    s.temperatureIn = suctionTemperatureK;
    s.rate = rateSm3PerS;

    if (curve)
    {
        s.mode = 0;
        const double shutoff =
            c.shutoffRatio > 0.0 ? c.shutoffRatio : c.pressureRatio;
        if (!(shutoff >= c.pressureRatio))
            cfFatal("compressor: the curve rises with rate (shutoffRatio < "
                    "pressureRatio): a machine cannot boost harder the faster "
                    "it is asked to run");
        if (!(c.ratedRate > 0.0))
            cfFatal("compressor: the ratio curve needs a positive ratedRate: "
                    "the rated point is where pressureRatio is quoted");
        if (!(c.pressureRatio > 1.0))
            cfFatal("compressor: pressureRatio must exceed unity: a machine "
                    "at or below 1 is a restriction, not a compressor");
        const double x = rateSm3PerS / c.ratedRate;
        s.ratio = shutoff - (shutoff - c.pressureRatio) * x * x;
        if (!(s.ratio >= 1.0))
            cfFatal("compressor: the rate is past the curve's unity-ratio "
                    "point: beyond it the declared quadratic stops being a "
                    "compressor statement");
        s.head = c.compressibility * c.specificGasConstant
               * suctionTemperatureK / ex * (std::pow(s.ratio, ex) - 1.0);
    }
    else
    {
        s.mode = 1;
        s.head = c.polytropicHead;
        s.ratio = std::pow(1.0 + s.head * ex
                           / (c.compressibility * c.specificGasConstant
                              * suctionTemperatureK),
                           1.0 / ex);
    }

    s.pressureOut = suctionPressurePa * s.ratio;
    s.massRate = suctionPressurePa * rateSm3PerS
               / (c.compressibility * c.specificGasConstant
                  * suctionTemperatureK);
    s.powerIdeal = s.massRate * s.head;
    s.powerShaft = s.powerIdeal / c.polytropicEfficiency;
    s.temperatureOut =
        suctionTemperatureK
        * std::pow(s.ratio, ex / c.polytropicEfficiency);
    return s;
}

} // namespace cf
