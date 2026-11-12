#include "cf/1d/BeggsBrill.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"
#include "cf/pvt/BlackOil.h"

namespace cf
{

namespace
{

constexpr double kGravity = 9.80665;

} // namespace

BeggsBrillOut beggsBrillDpDs(const BeggsBrillIn& in)
{
    BeggsBrillOut out;
    const double vm = in.vsl + in.vsg;
    if (!(vm > 1e-9) || !(in.diameter > 0))
    {
        // no flow: beggsBrillDpDs is about the flowing gradient; the static
        // column is the caller's (the march evaluates the no-slip mixture)
        return out;
    }

    const double lambda = in.vsl / vm;
    const double fr = vm * vm / (kGravity * in.diameter);
    out.lambdaL = lambda;

    const double l1 = 316.0 * std::pow(lambda, 0.302);
    const double l2 = 0.0009252 * std::pow(lambda, -2.4684);
    const double l3 = 0.1 * std::pow(lambda, -1.4516);
    const double l4 = 0.5 * std::pow(lambda, -6.738);

    // regime precedence exactly as published; the last matching rule wins
    int regime;
    if ((lambda < 0.01 && fr < l1) || (lambda >= 0.01 && fr < l2)) regime = 0;
    else if (lambda >= 0.01 && fr >= l2 && fr <= l3) regime = 1;
    else if ((0.01 <= lambda && lambda < 0.4 && fr > l3 && fr <= l1)
             || (lambda >= 0.4 && fr > l3 && fr <= l4)) regime = 2;
    else regime = 3;
    out.regime = regime;

    // horizontal holdup with the >= lambda clamp; index 1 (transition) is
    // never taken directly - it blends hl0(0) and hl0(2)
    const double ab[4] = {0.98, 0.0, 0.845, 1.065};
    const double bb[4] = {0.4846, 0.0, 0.5351, 0.5824};
    const double cb[4] = {0.0868, 0.0, 0.0173, 0.0609};
    auto hl0 = [&](int r)
    {
        const double v = ab[r] * std::pow(lambda, bb[r]) * std::pow(fr, -cb[r]);
        return std::max(v, lambda);
    };
    double hl;
    if (regime == 1)
    {
        const double w = (l3 - fr) / (l3 - l2);
        hl = w * hl0(0) + (1.0 - w) * hl0(1);
    }
    else
    {
        hl = hl0(regime);
    }

    // inclination correction
    const double theta = in.angleDeg * M_PI / 180.0;
    if (std::abs(in.angleDeg) > 1e-12 && lambda > 1e-9)
    {
        double d, e, f, h;
        if (in.angleDeg > 0 && regime != 3)
        {
            if (regime == 0) { d = 0.011; e = -3.768; f = 3.539; h = -1.614; }
            else { d = 2.96; e = 0.305; f = -0.4473; h = 0.0978; }
        }
        else if (in.angleDeg < 0)
        {
            d = 4.70; e = -0.3692; f = 0.1244; h = -0.5056;
        }
        else
        {
            d = 0.0; // uphill distributed: no correction
        }
        if (d > 0.0)
        {
            const double nlv =
                in.vsl * std::sqrt(std::sqrt(in.rhoL / (kGravity * in.sigma)));
            double c = (1.0 - lambda)
                     * std::log(d * std::pow(lambda, e) * std::pow(nlv, f)
                                * std::pow(fr, h));
            c = std::max(c, 0.0);
            const double s1 = std::sin(1.8 * theta);
            hl *= 1.0 + c * (s1 - s1 * s1 * s1 / 3.0);
        }
    }
    hl = std::min(std::max(hl, 0.0), 1.0);
    out.holdup = hl;

    const double rhoM = in.rhoL * lambda + in.rhoG * (1.0 - lambda);
    const double rhoS = in.rhoL * hl + in.rhoG * (1.0 - hl);
    // 5.76 (A5 increment 2): the liquid phase's apparent viscosity priced
    // on the mixture wall rate. The same lambda (no-slip) weighting the
    // correlation has always used is kept, and nothing else in the
    // holdup / flow-pattern / inclination correlations reads a viscosity,
    // so a declaration moves the friction term and only the friction term.
    // Zero declaration returns in.muL bit for bit.
    const double muL = wallRateApparentViscosityPaS(
        in.liquidRheology, 8.0 * vm / in.diameter, in.muL);
    const double muM = muL * lambda + in.muG * (1.0 - lambda);
    const double re = rhoM * vm * in.diameter / muM;
    const double fNs = pvt::colebrookDarcyFriction(re, in.roughness / in.diameter);
    out.fNs = fNs;

    // two-phase friction multiplier on y = lambda / HL^2
    double s;
    if (lambda < 1e-9 || hl < 1e-9)
    {
        s = 0.0; // single-phase degenerate: no slippage correction
    }
    else
    {
        const double y = lambda / (hl * hl);
        const double ly = std::log(y);
        if (y > 1.0 && y < 1.2) s = std::log(2.2 * y - 1.2);
        else s = ly / (ly * (ly * (0.01853 * ly * ly - 0.8725) + 3.182) - 0.0523);
        s = std::min(s, 7.0);
    }
    const double fTp = fNs * std::exp(s);
    out.fTp = fTp;

    const double ele = kGravity * std::sin(theta) * rhoS;
    const double fric = fTp * rhoM * vm * vm / (2.0 * in.diameter);
    const double p = std::max(in.pPa, 1000.0);
    const double ek = in.vsg * vm * rhoS / p;
    out.dpPerMeter = (ele + fric) / std::max(1.0 - ek, 0.1);
    return out;
}

BeggsBrillMarchOut marchBeggsBrillPipe(const pvt::BlackOilFluid& fluid,
                                       double qoStd, double pInPa,
                                       double length, double diameter,
                                       double roughness, double sinTheta,
                                       int segments,
                                       const Rheology& liquidRheology)
{
    BeggsBrillMarchOut out;
    const double pFloor = std::max(0.02 * pInPa, 1000.0);
    const double ds = length / static_cast<double>(std::max(segments, 1));
    const double area = 0.25 * M_PI * diameter * diameter;
    const double wor =
        fluid.waterCut > 0 && fluid.waterCut < 0.999
            ? fluid.waterCut / (1.0 - fluid.waterCut)
            : (fluid.waterCut >= 0.999 ? 999.0 : 0.0);
    out.massRate = qoStd * pvt::blackOilMassPerOilSm3(fluid);

    double p = pInPa;
    double lastHoldup = 0, lastLambda = 0, lastVm = 0, lastFTp = 0;
    const double angleDeg = std::asin(std::clamp(sinTheta, -1.0, 1.0))
                          * 180.0 / M_PI;
    for (int i = 0; i < segments; ++i)
    {
        const pvt::BlackOilState st = pvt::evaluateBlackOil(fluid, p);
        double rhoL = 0, muL = 0;
        pvt::blackOilLiquid(fluid, st, rhoL, muL);
        // in-situ volumetric rates at p: oil Bo, water 1, free gas bg
        const double ql = qoStd * st.bo + qoStd * wor;
        const double qg = qoStd * st.freeGas * st.bg;
        const double vsl = ql / area;
        const double vsg = qg / area;

        BeggsBrillIn in;
        in.rhoL = rhoL;
        in.rhoG = st.rhoGas;
        in.muL = muL;
        in.muG = st.muGas;
        in.sigma = fluid.surfaceTension;
        in.vsl = vsl;
        in.vsg = vsg;
        in.pPa = p;
        in.diameter = diameter;
        in.roughness = roughness;
        in.angleDeg = angleDeg;
        in.liquidRheology = liquidRheology;

        BeggsBrillOut g;
        if (vsl + vsg > 1e-9)
        {
            g = beggsBrillDpDs(in);
        }
        else
        {
            // static column on the local no-slip mixture density
            const double vT = st.bo + wor;
            const double rhoM = (st.rhoOil * st.bo + fluid.waterDensity * wor)
                              / vT;
            BeggsBrillOut s;
            s.dpPerMeter = kGravity * sinTheta * rhoM;
            s.holdup = st.bo / vT;
            s.lambdaL = st.bo / vT;
            s.regime = 3;
            s.fNs = 0;
            s.fTp = 0;
            g = s;
        }
        lastVm = vsl + vsg;
        lastFTp = g.fTp;

        out.holdupMax = std::max(out.holdupMax, g.holdup);
        out.regimeMin = std::min(out.regimeMin, g.regime);
        out.regimeMax = std::max(out.regimeMax, g.regime);
        lastHoldup = g.holdup;
        lastLambda = g.lambdaL;
        p -= g.dpPerMeter * ds;
        if (!(p > pFloor))
        {
            p = pFloor;
            out.pOutlet = p;
            out.dpTotal = pInPa - p;
            out.holdupOutlet = lastHoldup;
            out.lambdaOutlet = lastLambda;
            out.vmOutlet = lastVm;
            out.fTpOutlet = lastFTp;
            out.ok = false;
            return out;
        }
    }
    out.pOutlet = p;
    out.dpTotal = pInPa - p;
    out.holdupOutlet = lastHoldup;
    out.lambdaOutlet = lastLambda;
    out.vmOutlet = lastVm;
    out.fTpOutlet = lastFTp;
    out.ok = true;
    return out;
}

} // namespace cf
