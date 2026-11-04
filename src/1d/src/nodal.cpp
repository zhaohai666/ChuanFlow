#include "cf/1d/Nodal.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"
#include "cf/common/Log.h"

namespace cf
{

double IpR::rateMax() const
{
    switch (form)
    {
        case Form::linear:
            return productivityIndex * reservoirPressure;
        case Form::vogel:
            return maxRate;
        case Form::jones:
            return deliverability * std::pow(reservoirPressure, 2.0 * exponent);
    }
    return 0.0;
}

double IpR::pWf(double rate) const
{
    const double q = std::max(0.0, std::min(rate, rateMax()));
    switch (form)
    {
        case Form::linear:
            return reservoirPressure - q / productivityIndex;
        case Form::vogel:
        {
            // q/qmax = 1 - 0.2x - 0.8x^2 -> 0.8x^2 + 0.2x + (q/qmax - 1) = 0,
            // positive root: x = (-0.2 + sqrt(3.24 - 3.2 q/qmax)) / 1.6
            const double x = (-0.2 + std::sqrt(3.24 - 3.2 * q / maxRate)) / 1.6;
            return reservoirPressure * x;
        }
        case Form::jones:
        {
            const double arg = reservoirPressure * reservoirPressure
                             - std::pow(q / deliverability, 1.0 / exponent);
            return std::sqrt(std::max(arg, 0.0));
        }
    }
    return 0.0;
}

double IpR::rate(double pWf) const
{
    const double p = std::max(0.0, std::min(pWf, reservoirPressure));
    const double x = p / reservoirPressure;
    switch (form)
    {
        case Form::linear:
            return productivityIndex * (reservoirPressure - p);
        case Form::vogel:
            return maxRate * (1.0 - 0.2 * x - 0.8 * x * x);
        case Form::jones:
            return deliverability * std::pow(reservoirPressure * reservoirPressure - p * p, exponent);
    }
    return 0.0;
}

namespace
{

/// The bottom-hole inlet whose upward Beggs-Brill march lands exactly on the
/// wellhead pressure. pOutlet rises sublinearly in pIn (higher pressure packs
/// the free gas smaller and eases the friction), so the residual g(pIn) =
/// pOutlet - pWellhead is monotone and the bracket expands from the wellhead
/// pressure itself; a zero rate degenerates to the static column.
double wellboreInletForRate(const pvt::BlackOilFluid& fluid, const WellLift& well,
                            double rate, double wellheadPressure)
{
    if (!(well.length > 0) || !(well.diameter > 0))
        cfFatal("nodal: well needs positive length and diameter");
    const double sinTheta = well.elevation / well.length;
    auto g = [&](double pIn)
    {
        const BeggsBrillMarchOut m = marchBeggsBrillPipe(
            fluid, rate, pIn, well.length, well.diameter, well.roughness,
            sinTheta, 24);
        // a choked march returns the floor shortfall: continuous with the
        // unchoked branch and correctly scaled for the secant
        return m.ok ? m.pOutlet - wellheadPressure
                    : std::max(0.02 * pIn, 1000.0) - wellheadPressure;
    };

    double pLo = wellheadPressure; // the march only ever drops pressure: g <= 0
    double fLo = g(pLo);
    if (fLo > 0.0)
        cfFatal("nodal: the wellbore gains pressure toward the wellhead - "
                "a producer profile needs elevation >= 0");
    // the early return must stay far tighter than the outer nodal secant's
    // 1e-7 residual, or horizontal wells would have a hidden floor error in f
    if (fLo > -1e-12 * wellheadPressure) return pLo;
    double pHi = pLo + std::max(1e4, 0.05 * wellheadPressure);
    double fHi = g(pHi);
    int expansions = 0;
    while (fHi < 0.0)
    {
        pLo = pHi;
        fLo = fHi;
        pHi *= 2.0;
        fHi = g(pHi);
        if (++expansions > 100)
            cfFatal("nodal: the wellbore march cannot reach the wellhead "
                    "pressure at the imposed rate");
    }

    // Illinois secant between pLo (g<0) and pHi (g>0). The gates sit far
    // below the outer nodal secant's 1e-7 residual: this solve is nested
    // inside f(q), so its own floor has to be invisible there.
    double p = pLo + (pHi - pLo) * fLo / (fLo - fHi);
    for (int it = 0; it < 100; ++it)
    {
        const double fp = g(p);
        if (std::abs(fp) < 1e-12 * wellheadPressure) break;
        if (fp > 0.0)
        {
            pHi = p;
            fHi = fp;
            fHi *= 0.5;
        }
        else
        {
            pLo = p;
            fLo = fp;
            fLo *= 0.5;
        }
        const double pNew = pLo + (pHi - pLo) * fLo / (fLo - fHi);
        p = (pNew > pLo && pNew < pHi) ? pNew : 0.5 * (pLo + pHi);
        if (std::abs(pHi - pLo) < 1e-12 * pHi) break;
    }
    return p;
}

} // namespace

double vlpBottomHolePressure(const pvt::BlackOilFluid& fluid, const WellLift& well,
                             double rate, double wellheadPressure)
{
    return wellboreInletForRate(fluid, well, rate, wellheadPressure);
}

NodalResult solveNodalAnalysis(const IpR& ipr, const pvt::BlackOilFluid& fluid,
                               const WellLift& well, double wellheadPressure,
                               const NodalControls& ctl)
{
    if (!(ipr.reservoirPressure > 0) || !(ipr.rateMax() > 0))
        cfFatal("nodal: the IPR needs positive reservoirPressure and capacity");
    if (!(wellheadPressure > 0))
        cfFatal("nodal: wellheadPressure must be positive");

    NodalResult res;
    res.staticBottomPressure = wellboreInletForRate(fluid, well, 0.0, wellheadPressure);

    // f(q) = p_vlp(q) - p_ipr(q): p_vlp is U-shaped (holdup relief first,
    // friction later) and p_ipr falls monotonically, so the stable operating
    // point is the LARGEST root - bracket from the high-rate end downward,
    // exactly like the network pipe law, never from q = 0 upward.
    const double qMax = ipr.rateMax();
    auto f = [&](double q)
    { return wellboreInletForRate(fluid, well, q, wellheadPressure) - ipr.pWf(q); };

    double qHi = qMax;
    double fHi = f(qHi); // p_ipr(qMax) = 0, so fHi = p_vlp(qMax) > 0
    double qLo = 0.0, fLo = 0.0;
    bool bracketed = false;
    double qProbe = qHi;
    for (int k = 0; k < 60 && !bracketed; ++k)
    {
        qProbe *= 0.5;
        const double fk = f(qProbe);
        if (fk <= 0.0)
        {
            qLo = qProbe;
            fLo = fk;
            bracketed = true;
        }
    }
    if (!bracketed && f(0.0) <= 0.0)
    {
        // the crossing hides below the finest probe: the bracket still ends
        // at the static rate
        qLo = 0.0;
        fLo = f(0.0);
        bracketed = true;
    }
    if (!bracketed)
    {
        // no rate in reach makes f negative: the reservoir cannot lift the
        // static column, the well is dead
        res.bottomHolePressure = res.staticBottomPressure;
        res.converged = true;
        return res;
    }

    // Illinois secant between qLo (f<0) and qHi (f>0)
    const double pScale = std::max(ipr.reservoirPressure, wellheadPressure);
    double q = qLo + (qHi - qLo) * fLo / (fLo - fHi);
    bool converged = false;
    for (int it = 0; it < ctl.maxIter; ++it)
    {
        const double fq = f(q);
        res.iterations = it + 1;
        if (std::abs(fq) < 1e-7 * pScale)
        {
            converged = true;
            break;
        }
        if (fq > 0.0)
        {
            qHi = q;
            fHi = fq;
            fHi *= 0.5;
        }
        else
        {
            qLo = q;
            fLo = fq;
            fLo *= 0.5;
        }
        const double qNew = qLo + (qHi - qLo) * fLo / (fLo - fHi);
        q = (qNew > qLo && qNew < qHi) ? qNew : 0.5 * (qLo + qHi);
        if (std::abs(qHi - qLo) < ctl.tol * std::max(qHi, 1e-30))
        {
            converged = true;
            break;
        }
    }

    res.rate = q;
    res.bottomHolePressure = wellboreInletForRate(fluid, well, q, wellheadPressure);
    const BeggsBrillMarchOut m = marchBeggsBrillPipe(
        fluid, q, res.bottomHolePressure, well.length, well.diameter,
        well.roughness, well.elevation / well.length, 24);
    res.holdup = m.holdupOutlet;
    res.regime = m.regimeMax;
    res.flows = true;
    res.converged = converged;
    return res;
}

NodalCase loadNodalCaseFromDict(const Dictionary& caseDict)
{
    NodalCase c;

    if (!caseDict.has("fluid")
        || !caseDict.getDict("fluid").has("blackoil"))
        cfFatal("nodal: case needs a fluid{ blackoil{...} } block");
    const Dictionary& bd = caseDict.getDict("fluid").getDict("blackoil");
    c.fluid.gasGravity = bd.getScalarOr("gasGravity", c.fluid.gasGravity);
    c.fluid.oilDensityStd = bd.getScalarOr("oilDensity", c.fluid.oilDensityStd);
    c.fluid.gorStd = bd.getScalarOr("gor", c.fluid.gorStd);
    c.fluid.waterCut = bd.getScalarOr("waterCut", c.fluid.waterCut);
    c.fluid.temperature = bd.getScalarOr("temperature", c.fluid.temperature);
    c.fluid.surfaceTension = bd.getScalarOr("surfaceTension", c.fluid.surfaceTension);
    if (!(c.fluid.oilDensityStd > 0))
        cfFatal("nodal: blackoil needs positive oilDensity");

    if (!caseDict.has("ipr")) cfFatal("nodal: case has no ipr{} block");
    const Dictionary& id = caseDict.getDict("ipr");
    c.ipr.reservoirPressure = id.getScalar("reservoirPressure");
    const std::string form = id.getWordOr("form", "vogel");
    if (form == "linear")
    {
        c.ipr.form = IpR::Form::linear;
        c.ipr.productivityIndex = id.getScalar("productivityIndex");
    }
    else if (form == "vogel")
    {
        c.ipr.form = IpR::Form::vogel;
        c.ipr.maxRate = id.getScalar("maxRate");
    }
    else if (form == "jones")
    {
        c.ipr.form = IpR::Form::jones;
        c.ipr.deliverability = id.getScalar("deliverability");
        c.ipr.exponent = id.getScalarOr("exponent", 1.0);
    }
    else
        cfFatal(std::string("nodal: ipr form must be linear, vogel or jones, not ")
                + form);
    if (!(c.ipr.reservoirPressure > 0) || !(c.ipr.rateMax() > 0))
        cfFatal("nodal: ipr needs positive reservoirPressure and capacity");

    if (!caseDict.has("well")) cfFatal("nodal: case has no well{} block");
    const Dictionary& wd = caseDict.getDict("well");
    c.well.length = wd.getScalar("length");
    c.well.diameter = wd.getScalar("diameter");
    c.well.roughness = wd.getScalarOr("roughness", 0.0);
    c.well.elevation = wd.getScalarOr("elevation", 0.0);

    c.wellheadPressure = caseDict.getScalar("wellheadPressure");
    if (!(c.wellheadPressure > 0))
        cfFatal("nodal: wellheadPressure must be positive");

    if (caseDict.has("solve"))
    {
        const Dictionary& sd = caseDict.getDict("solve");
        c.controls.maxIter = static_cast<int>(sd.getScalarOr("maxIter", c.controls.maxIter));
        c.controls.tol = sd.getScalarOr("tol", c.controls.tol);
    }
    return c;
}

} // namespace cf
