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

/// Parse one ipr{} block. `who` prefixes the diagnostics: empty for the main
/// bore's own curve (the single-well messages stay byte-identical), "lateral
/// 'L1' " for a 5.157 lateral.
void parseIpR(const Dictionary& id, const std::string& who, IpR& ipr)
{
    ipr.reservoirPressure = id.getScalar("reservoirPressure");
    const std::string form = id.getWordOr("form", "vogel");
    if (form == "linear")
    {
        ipr.form = IpR::Form::linear;
        ipr.productivityIndex = id.getScalar("productivityIndex");
    }
    else if (form == "vogel")
    {
        ipr.form = IpR::Form::vogel;
        ipr.maxRate = id.getScalar("maxRate");
    }
    else if (form == "jones")
    {
        ipr.form = IpR::Form::jones;
        ipr.deliverability = id.getScalar("deliverability");
        ipr.exponent = id.getScalarOr("exponent", 1.0);
    }
    else
        cfFatal("nodal: " + who +
                "ipr form must be linear, vogel or jones, not " + form);
    if (!(ipr.reservoirPressure > 0) || !(ipr.rateMax() > 0))
        cfFatal("nodal: " + who +
                "ipr needs positive reservoirPressure and capacity");
}

/// Parse one well{} bore. Same `who` convention as parseIpR; with an empty
/// prefix every refusal reads exactly the single-well wording it always had.
void parseWell(const Dictionary& wd, const std::string& who, WellLift& well)
{
    well.length = wd.getScalar("length");
    well.diameter = wd.getScalar("diameter");
    well.roughness = wd.getScalarOr("roughness", 0.0);
    well.elevation = wd.getScalarOr("elevation", 0.0);
    // 5.88 (D4 increment 2): the writable flow-pattern dispatch is a network
    // pipe declaration. This march keeps the map's own choice, and a key here
    // would be silently dropped rather than honoured - so it is refused.
    if (wd.has("flowPattern"))
        cfFatal("nodal: " + who + "well declares flowPattern: the writable "
                "flow-pattern dispatch (5.88) is a steady-network pipe "
                "declaration; the VLP well marches with the flow-pattern "
                "map's own choice");
    // 5.143 (G4): same reason. The alternative holdup laws are a network
    // pipe declaration; this march prices the Beggs-Brill law and a key
    // here would be silently dropped rather than honoured.
    if (wd.has("correlation"))
        cfFatal("nodal: " + who + "well declares correlation: the alternative "
                "holdup laws (5.143) are a steady-network pipe declaration; "
                "the VLP well marches with the Beggs-Brill holdup law");
    // 5.144 (G6): same reason. The passive temperature march is a network
    // pipe declaration; the VLP march is isothermal and reads none of the
    // thermal keys.
    if (wd.has("jouleThomson") || wd.has("wallHeatTransfer")
        || wd.has("wallHeatTransferSegments") || wd.has("wallModel")
        || wd.has("wallTemperature") || wd.has("surfaceTemperature")
        || wd.has("geothermalGradient") || wd.has("heatCapacity"))
        cfFatal("nodal: " + who + "well declares a thermal parameter: the "
                "steady temperature march (5.144) is a steady-network pipe "
                "declaration; the VLP march is isothermal");
    // 5.89 (D1 increment 3): same reason. A VLP well is marched with a gas rate
    // it is handed; the valve that would have set that rate lives on the steady
    // network, so the keys are refused rather than left unread.
    if (wd.has("liftDomePressure") || wd.has("liftAreaRatio")
        || wd.has("liftSpringPressure") || wd.has("liftDischargeCoeff")
        || wd.has("liftKappa"))
        cfFatal("nodal: " + who + "well declares a lift valve: the gas-lift "
                "valve (5.89) is a steady-network injection-string element "
                "priced on node pressures; the VLP march is handed its gas "
                "rate and does not own a bellows to open");
    // 5.90 (D1 increment 4): same reason. The rod pump takes its rate from the
    // stroke and the barrel, which is the opposite of the VLP march asking the
    // column what rate it will carry.
    if (wd.has("rodPlungerDiameter") || wd.has("rodStrokeLength")
        || wd.has("rodStrokesPerMinute") || wd.has("rodDiameter")
        || wd.has("rodLength") || wd.has("rodYoungModulus")
        || wd.has("rodDensity"))
        cfFatal("nodal: " + who + "well declares a rod pump: the subsurface "
                "rod pump (5.90) is a steady-network element whose rate is "
                "its stroke times its barrel; the VLP march solves for the "
                "rate the reservoir and the column agree on and does not own "
                "a polished rod to pull a card against");
    if (!(well.length > 0) || !(well.diameter > 0))
        cfFatal("nodal: " + who + "well needs positive length and diameter");
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

MultiLateralResult solveNodalLaterals(const pvt::BlackOilFluid& fluid,
                                      const WellLift& mainWell,
                                      double wellheadPressure,
                                      const std::vector<LateralSpec>& laterals,
                                      const NodalControls& ctl)
{
    if (laterals.empty())
        cfFatal("nodal: the multi-lateral solve needs at least one lateral");
    if (!(wellheadPressure > 0))
        cfFatal("nodal: wellheadPressure must be positive");
    double pResMax = 0.0;
    for (const LateralSpec& lat : laterals)
    {
        if (!(lat.ipr.reservoirPressure > 0) || !(lat.ipr.rateMax() > 0))
            cfFatal("nodal: lateral '" + lat.name
                    + "' ipr needs positive reservoirPressure and capacity");
        if (!(lat.well.length > 0) || !(lat.well.diameter > 0))
            cfFatal("nodal: lateral '" + lat.name
                    + "' well needs positive length and diameter");
        pResMax = std::max(pResMax, lat.ipr.reservoirPressure);
    }

    MultiLateralResult out;
    out.main.staticBottomPressure =
        wellboreInletForRate(fluid, mainWell, 0.0, wellheadPressure);

    // q_i(p_j): lateral i's own IPR x VLP operating point at back-pressure p_j
    // (the junction, i.e. the main bore's bottom hole). A lateral whose
    // reservoir cannot lift its column against p_j contributes rate 0.
    auto junctionRates = [&](double pJ)
    {
        double sum = 0.0;
        for (const LateralSpec& lat : laterals)
            sum += solveNodalAnalysis(lat.ipr, fluid, lat.well, pJ, ctl).rate;
        return sum;
    };
    // F(p_j) = p_j - p_vlp_main(sum q_i(p_j)). Each q_i falls with the back-
    // pressure and the main VLP rises with the rate it carries, so F rises:
    // one root, bracketed and cut by the same Illinois secant family as the
    // single-well nodal and the network pipe law.
    auto f = [&](double pJ)
    {
        return pJ - wellboreInletForRate(fluid, mainWell, junctionRates(pJ),
                                         wellheadPressure);
    };

    // the march only ever drops pressure, so p_vlp_main >= wellheadPressure
    // and F(wellheadPressure) <= 0: the low end of the bracket is the wellhead
    double pLo = wellheadPressure;
    double fLo = f(pLo);
    // at or above the largest reservoir pressure every lateral is dead, the
    // total is 0 and F reads p_j - static > 0: the high end is within reach
    double pHi = std::max({2.0 * pLo, 1.05 * out.main.staticBottomPressure,
                           1.05 * pResMax}) + 1e4;
    double fHi = f(pHi);
    int expansions = 0;
    while (fHi <= 0.0)
    {
        pLo = pHi;
        fLo = fHi;
        pHi *= 2.0;
        fHi = f(pHi);
        if (++expansions > 100)
            cfFatal("nodal: the junction residual never turns positive - the "
                    "laterals keep flowing at every back-pressure");
    }

    const double pScale = std::max({out.main.staticBottomPressure,
                                    wellheadPressure, pResMax});
    double p = pLo + (pHi - pLo) * fLo / (fLo - fHi);
    bool converged = false;
    for (int it = 0; it < ctl.maxIter; ++it)
    {
        const double fp = f(p);
        out.iterations = it + 1;
        if (std::abs(fp) < 1e-7 * pScale)
        {
            converged = true;
            break;
        }
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
        if (std::abs(pHi - pLo) < ctl.tol * std::max(pHi, 1e-30))
        {
            converged = true;
            break;
        }
    }

    // the final readout at the converged junction: every lateral's operating
    // point against p, the total by exact summation (mass conservation at the
    // junction is a construction, not a tolerance), and the main column's
    // state carrying that total up to the wellhead
    double total = 0.0;
    for (const LateralSpec& lat : laterals)
    {
        const NodalResult lr =
            solveNodalAnalysis(lat.ipr, fluid, lat.well, p, ctl);
        LateralResult rec;
        rec.name = lat.name;
        rec.rate = lr.rate;
        rec.bottomHolePressure = lr.bottomHolePressure;
        rec.staticBottomPressure = lr.staticBottomPressure;
        rec.holdup = lr.holdup;
        rec.regime = lr.regime;
        rec.flows = lr.flows;
        rec.converged = lr.converged;
        total += lr.rate;
        out.laterals.push_back(rec);
    }
    out.main.rate = total;
    out.main.bottomHolePressure = p;
    const BeggsBrillMarchOut m = marchBeggsBrillPipe(
        fluid, total, p, mainWell.length, mainWell.diameter, mainWell.roughness,
        mainWell.elevation / mainWell.length, 24);
    out.main.holdup = m.holdupOutlet;
    out.main.regime = m.regimeMax;
    out.main.flows = total > 0.0;
    out.main.iterations = out.iterations;
    out.main.converged = converged;
    out.converged = converged;
    return out;
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
    parseIpR(caseDict.getDict("ipr"), "", c.ipr);

    if (!caseDict.has("well")) cfFatal("nodal: case has no well{} block");
    parseWell(caseDict.getDict("well"), "", c.well);

    // 5.157 (P5): the multi-lateral completion. Each named child of laterals{}
    // carries its own ipr{} and well{}; the main bore above keeps its VLP and
    // becomes the lift conduit for the summed lateral flow. An empty laterals{}
    // block is refused - a declared completion with zero branches is a typo,
    // not a single-well case (that case simply omits the block). Heterogeneous
    // IPR forms across the laterals (and against the main ipr{}) are honoured:
    // each branch brings its own inflow physics.
    if (caseDict.has("laterals"))
    {
        const Dictionary& ld = caseDict.getDict("laterals");
        if (ld.keys().empty())
            cfFatal("nodal: laterals{} declares no lateral: a multi-lateral "
                    "case needs at least one named lateral block");
        for (const std::string& name : ld.keys())
        {
            const Dictionary& sd = ld.getDict(name);
            LateralSpec lat;
            lat.name = name;
            const std::string who = "lateral '" + name + "' ";
            if (!sd.has("ipr"))
                cfFatal("nodal: " + who + "has no ipr{} block: a lateral "
                        "carries its own inflow curve");
            parseIpR(sd.getDict("ipr"), who, lat.ipr);
            if (!sd.has("well"))
                cfFatal("nodal: " + who + "has no well{} block: a lateral "
                        "carries its own VLP bore");
            parseWell(sd.getDict("well"), who, lat.well);
            c.laterals.push_back(std::move(lat));
        }
    }

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
