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

// 5.143 (G4) law 1: the Hagedorn-Brown liquid holdup as restated by Whitson,
// Brule & Stephenson, "Revisiting Hagedorn-Brown Correlation for Vertical
// Two-Phase Flow" (JPT 69(6), 2017; SPE 169350), Appendix E. SI inputs; the
// four dimensionless groups, the two y3 branches and the psi rational fit
// are quoted there and were transcribed against the paper's worked example.
// psi is the paper's rational fit to its own chart over
// Ngv Nmu^0.38 / Nd^2.14 - clamped at 0.088 above and set to 1 below 0.01;
// the chart-reading caveat is the paper's, not ours. The law is bounded
// below by the no-slip content lambda (a holdup cannot undercut what the
// mixture actually carries) and above by 1.
double hagedornBrownHoldup(double vsl, double vsg, double lambda,
                           double rhoL, double muL, double sigma,
                           double d, double pPa)
{
    if (!(vsg > 0.0)) return 1.0; // no gas: the column is all liquid
    if (!(vsl > 0.0)) return 0.0; // no liquid: nothing to hold
    const double rhoSig = rhoL / (kGravity * sigma);
    const double nlv = vsl * std::pow(rhoSig, 0.25);
    const double ngv = vsg * std::pow(rhoSig, 0.25);
    const double nd = d * std::sqrt(kGravity * rhoL / sigma);
    const double nmu = muL
        * std::pow(kGravity / (rhoL * sigma * sigma * sigma), 0.25);

    const double x1 = std::min(std::max(nmu, 0.002), 0.5);
    const double nmuC = 9.0950e-2 * x1 * x1 * x1 - 1.2183e-1 * x1 * x1
                      + 5.6348e-2 * x1 + 1.9122e-3;
    const double x2 = ngv * std::pow(nmu, 0.380) / std::pow(nd, 2.14);
    double psi = 1.0;
    if (x2 >= 0.01)
    {
        const double x = std::min(x2, 0.088);
        const double num = 3.6917e-1 * x * x * x + 1.0029e-2 * x * x
                         - 1.2769e-3 * x + 1.8089e-5;
        const double den = x * x * x * x - 4.894e-2 * x * x * x
                         + 2.7058e-2 * x * x - 1.2769e-3 * x + 1.9585e-5;
        psi = num / den;
    }
    const double x3 = (nlv / std::pow(ngv, 0.575))
                    * std::pow(pPa / 1.0e5, 0.10) * (nmuC / nd);
    double y3;
    if (x3 <= 4.6423e-5)
    {
        // the low-x3 asymptote of the published quadratic ratio, floored so
        // the denominator stays away from its root
        const double xx = std::max(x3, 2.1e-6);
        y3 = (-3.7902e-5 * xx + 4.326e-11)
           / (xx * xx - 1.9323e-4 * xx - 3.752e-10);
    }
    else if (x3 <= 5.0e-3)
    {
        y3 = (-18.53 * x3 * x3 + 1.147 * x3 + 6.4898e-5) / (x3 + 4.5319e-4);
    }
    else y3 = 0.98;
    return std::max(psi * y3, lambda);
}

// 5.143 (G4) law 2: the Duns & Ros liquid holdup region skeleton (API
// Bulletin 34, Vol I, Section C; Duns & Ros, BRP 1963). PROVENANCE WARNING:
// this is a transcription of the PUBLISHED SKELETON - the no-slip multiplier
// RL1 = 1 + 0.3 (Vsl + Vsg) below the L1 boundary, the linear blend to the
// no-slip value across the L1..L2 band, no-slip above L2 - and it has NOT
// been re-validated against the Bulletin's worked examples. The regression
// checker re-derives the SAME formula, so it proves the code matches the
// transcription, not the transcription matches the API. Treat numbers from
// this law as engineering estimates until the Bulletin tables are anchored.
// The Bulletin is written in field units (cp, lb/ft3, dynes/cm, ft/s, ft);
// the conversions below are exact, and the dimensionless groups NL, NLV, ND
// are the published ones.
double dunsRosHoldup(double vsl, double vsg, double lambda,
                    double rhoL, double muL, double sigma, double d)
{
    if (!(vsg > 0.0)) return 1.0; // no gas: the column is all liquid
    if (!(vsl > 0.0)) return 0.0; // no liquid: nothing to hold
    const double muCp = muL * 1000.0;
    const double rhoLb = rhoL / 16.0184634;
    const double sigDyn = sigma * 1000.0;
    const double vslFt = vsl * 3.28083989501312;
    const double vsgFt = vsg * 3.28083989501312;
    const double dFt = d * 3.28083989501312;

    const double nl = 0.15726 * muCp
                    / std::pow(rhoLb * sigDyn * sigDyn * sigDyn, 0.25);
    const double nlv = 1.938 * vslFt * std::pow(rhoLb / sigDyn, 0.25);
    const double nd = 120.872 * dFt * std::sqrt(rhoLb / sigDyn);
    const double l1 = 13.2 * (1.0 - lambda);
    const double l2 = l1 + 0.42 * std::pow(nlv, 0.25) * std::pow(nl, 0.10)
                    * std::pow(nd, 0.25);
    const double rl1 = 1.0 + 0.3 * (vslFt + vsgFt);
    double hl;
    if (vsgFt < l1)
        hl = lambda * rl1;
    else if (vsgFt <= l2)
        hl = lambda * (rl1 + (vsgFt - l1) / (l2 - l1) * (1.0 - rl1));
    else
        hl = lambda;
    return std::min(std::max(hl, 0.0), 1.0);
}

// 5.155 (G4 family): the OLGAS-steady liquid holdup - a transcription of the
// PUBLISHED MECHANISTIC SKELETON of the steady OLGA family (the Taitel-Dukler
// regime-dispatched steady mechanistic holdup, Taitel & Dukler, AIChE J. 22:47,
// 1976, in the lineage of the unified steady models), NOT the proprietary
// OLGA/PIPESIM model. PROVENANCE WARNING: this skeleton has NOT been
// re-validated against OLGA-steady or PIPESIM output. The regression checker
// re-derives the SAME formula, so it proves the code matches the
// transcription, not the transcription matches OLGA. Treat numbers from this
// law as engineering estimates, exactly as 5.143 declares for Duns-Ros.
// The dispatch reads the Beggs-Brill map's OWN verdict (flowPattern is refused
// alongside any alternative law, so regime == regimeFromMap here) and replaces
// only the horizontal + inclination holdup block:
//   - intermittent (the map's plug/slug band): the section-average is the
//     slug-body content weighted by the distance to the Taitel-Dukler
//     intermittent -> dispersed locus, HL = beta*HLS + (1-beta)*HLTB with the
//     DECLARED simplifications HLS = 1 (the slug body is a full liquid column;
//     the dynamic 1 - rhoL/rhoS correction is dropped) and HLTB = lambda (the
//     train between slugs carries its liquid at the no-slip content). beta is
//     read off the published locus U*GS = 3.1 sqrt(101 - 50 U*SL) on the
//     dimensionless velocities U*G = U_GS sqrt(rho_G / (g D (rhoL - rhoG)))
//     (and liquid analogue): deep in the band (low gas) beta -> 1, at the
//     locus beta -> 0, and the locus pinches to zero gas velocity at
//     U*SL = 2.02, beyond which the blend is fully dispersed-weighted;
//   - segregated (stratified): the no-slip circular-segment geometry of the
//     horizontal bore (the interface LEVEL whose segment area equals lambda,
//     inverted from the closed-form segment area) taken as the holdup proxy,
//     times the DECLARED inclination gravity factor (1 + sin theta): uphill
//     backs liquid up, downhill drains it;
//   - distributed (dispersed bubble / mist): homogeneous no-slip plus the
//     Zuber-Findlay / Zuber-Hill drift correction, HL = 1 - vsg/(C0 vm) with
//     the published bubbly distribution parameter C0 = 1.2 and zero drift
//     velocity (declared engineering transcription), floored at lambda;
//   - the map's transition band blends segregated and intermittent with the
//     map's own published weight, mirroring Beggs-Brill's blend structure.
double segmentAreaFraction(double level)
{
    // area fraction below the depth level*D of a circular bore (closed form)
    const double x = std::min(std::max(level, 0.0), 1.0);
    const double c = 1.0 - 2.0 * x;
    return (std::acos(c) - c * std::sqrt(std::max(4.0 * x * (1.0 - x), 0.0)))
           / M_PI;
}

double olgasSteadyHoldup(double vsl, double vsg, double lambda, double rhoL,
                         double rhoG, double d, double sinTheta, int regime,
                         double transitionWeight)
{
    if (!(vsg > 0.0)) return 1.0; // no gas: the column is all liquid
    if (!(vsl > 0.0)) return 0.0; // no liquid: nothing to hold
    const double vm = vsl + vsg;
    const double drho = std::max(rhoL - rhoG, 1e-6);

    // segregated: invert the horizontal segment geometry for the level whose
    // area fraction equals lambda (f is monotone, 60 bisections land inside
    // any double the march can read)
    double lo = 0.0, hi = 1.0;
    for (int it = 0; it < 60; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        if (segmentAreaFraction(mid) < lambda) lo = mid;
        else hi = mid;
    }
    const double hlSeg = std::clamp(0.5 * (lo + hi) * (1.0 + sinTheta),
                                    0.0, 1.0);

    // intermittent: Taitel-Dukler slug-content weight against the locus
    const double ugs = vsg * std::sqrt(rhoG / (kGravity * d * drho));
    const double usl = vsl * std::sqrt(rhoL / (kGravity * d * drho));
    const double ugsb = 3.1 * std::sqrt(std::max(101.0 - 50.0 * usl, 0.0));
    const double beta = ugsb > 0.0
        ? std::clamp(1.0 - ugs / ugsb, 0.0, 1.0)
        : 0.0;
    const double hlInt = std::clamp(beta + (1.0 - beta) * lambda, 0.0, 1.0);

    // distributed: Zuber-Hill family drift correction on the no-slip mix
    const double hlDist = std::clamp(1.0 - vsg / (1.2 * vm), lambda, 1.0);

    double hl;
    if (regime == 0) hl = hlSeg;
    else if (regime == 3) hl = hlDist;
    else if (regime == 2) hl = hlInt;
    else hl = transitionWeight * hlSeg + (1.0 - transitionWeight) * hlInt;
    return std::min(std::max(hl, 0.0), 1.0);
}

} // namespace

// 5.144 (G6): the one thermal segment step, defined here so the black-oil
// march and the pvt-table march (network.cpp) share it and cannot drift.
double steadyThermalAdvance(const SteadyThermal& th, double t,
                            double pOld, double pNew, int segment,
                            double ds, double sinTheta, double massRate)
{
    // 5.189: the EOS coefficient, when the fluid can price it, is taken at
    // the segment INLET state - the same state the explicit march is holding.
    const double jt = th.jtOfState ? (*th.jtOfState)(t, pOld) : th.jouleThomson;
    t += jt * (pNew - pOld);
    const double u = th.segmentU ? (*th.segmentU)[segment]
                                 : th.wallHeatTransfer;
    if (u > 0.0 && massRate > 0.0)
    {
        const double depth = -(segment + 0.5) * sinTheta * ds;
        const double tw = th.wallModel == 1
            ? th.surfaceTemperature + th.geothermalGradient * depth
            : th.wallTemperature;
        t = tw + (t - tw) * std::exp(-u * ds / (massRate * th.heatCapacity));
    }
    return t;
}

BeggsBrillPattern beggsBrillPattern(double lambda, double nfr)
{
    BeggsBrillPattern pat;
    pat.l1 = 316.0 * std::pow(lambda, 0.302);
    pat.l2 = 0.0009252 * std::pow(lambda, -2.4684);
    pat.l3 = 0.1 * std::pow(lambda, -1.4516);
    pat.l4 = 0.5 * std::pow(lambda, -6.738);

    // regime precedence exactly as published; the last matching rule wins
    if ((lambda < 0.01 && nfr < pat.l1) || (lambda >= 0.01 && nfr < pat.l2))
        pat.regime = 0;
    else if (lambda >= 0.01 && nfr >= pat.l2 && nfr <= pat.l3)
        pat.regime = 1;
    else if ((0.01 <= lambda && lambda < 0.4 && nfr > pat.l3 && nfr <= pat.l1)
             || (lambda >= 0.4 && nfr > pat.l3 && nfr <= pat.l4))
        pat.regime = 2;
    else
        pat.regime = 3;

    if (pat.regime == 1)
        pat.transitionWeight = (pat.l3 - nfr) / (pat.l3 - pat.l2);
    return pat;
}

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

    // 5.85 (D4 increment 1): the flow-pattern map is a standalone dispatch
    // surface; the gradient reads its regime and boundaries from that single
    // source rather than re-deriving the map inline.
    const BeggsBrillPattern pat = beggsBrillPattern(lambda, fr);
    // 5.88 (D4 increment 2): the map above is still evaluated on every call
    // and still reported, but a forced regime decides which segment law prices
    // the holdup. regimeFromMap keeps the map's own verdict visible, so an
    // override is legible on the readout instead of silently replacing it.
    // Auto (forcedRegime < 0) selects pat.regime, which is exactly what 5.85
    // priced, bit for bit.
    const int regime = in.forcedRegime >= 0 ? in.forcedRegime : pat.regime;
    out.regime = regime;
    out.regimeFromMap = pat.regime;
    out.l1 = pat.l1;
    out.l2 = pat.l2;
    out.l3 = pat.l3;
    out.l4 = pat.l4;
    out.transitionWeight = pat.transitionWeight;

    const double theta = in.angleDeg * M_PI / 180.0;
    double hl;
    // the two branch holdups of the published transition (regime 1): each
    // neighbouring regime's law, each lifted with its OWN inclination set,
    // so each branch can be priced as a complete gradient below
    double hlSeg = 0.0, hlInt = 0.0;
    if (in.correlation == 0)
    {
    // horizontal holdup with the >= lambda clamp; index 1 (transition) has
    // no law of its own - the published transition blends the neighbouring
    // regimes' branch answers
    const double ab[4] = {0.98, 0.0, 0.845, 1.065};
    const double bb[4] = {0.4846, 0.0, 0.5351, 0.5824};
    const double cb[4] = {0.0868, 0.0, 0.0173, 0.0609};
    auto hl0 = [&](int r)
    {
        const double v = ab[r] * std::pow(lambda, bb[r]) * std::pow(fr, -cb[r]);
        return std::max(v, lambda);
    };
    // the published inclination C sets, applied per regime r: uphill
    // segregated and uphill intermittent carry their own sets, uphill
    // distributed none, and one downhill set prices every regime
    auto lifted = [&](double h, int r) -> double
    {
        if (std::abs(in.angleDeg) > 1e-12 && lambda > 1e-9)
        {
            double d, e, f, hExp;
            if (in.angleDeg > 0 && r != 3)
            {
                if (r == 0) { d = 0.011; e = -3.768; f = 3.539; hExp = -1.614; }
                else { d = 2.96; e = 0.305; f = -0.4473; hExp = 0.0978; }
            }
            else if (in.angleDeg < 0)
            {
                d = 4.70; e = -0.3692; f = 0.1244; hExp = -0.5056;
            }
            else
            {
                return h; // uphill distributed: no correction
            }
            const double nlv =
                in.vsl * std::sqrt(std::sqrt(in.rhoL / (kGravity * in.sigma)));
            double c = (1.0 - lambda)
                     * std::log(d * std::pow(lambda, e) * std::pow(nlv, f)
                                * std::pow(fr, hExp));
            c = std::max(c, 0.0);
            const double s1 = std::sin(1.8 * theta);
            h *= 1.0 + c * (s1 - s1 * s1 * s1 / 3.0);
        }
        return h;
    };
    if (regime == 1)
    {
        // the published transition (Beggs & Brill 1973): the gradient is
        // the blend of the two NEIGHBOURING regimes' complete branch
        // answers, not one branch priced on a blended holdup. Each branch
        // carries its own horizontal holdup (segregated HL0(0) /
        // intermittent HL0(2)), its own inclination C set and therefore its
        // own full dp/ds below; the reported holdup blends the branch
        // holdups, which makes the L2 / L3 boundaries continuous with the
        // segregated and intermittent regimes.
        const double w = pat.regime == 1
            ? pat.transitionWeight
            : (pat.l3 - fr) / (pat.l3 - pat.l2);
        hlSeg = lifted(hl0(0), 0);
        hlInt = lifted(hl0(2), 2);
        hl = w * hlSeg + (1.0 - w) * hlInt;
    }
    else
    {
        hl = lifted(hl0(regime), regime);
    }
    }
    else if (in.correlation == 1)
    {
        // 5.143 (G4): the Hagedorn-Brown / Whitson holdup REPLACES the
        // horizontal + inclination block above and nothing else: the map,
        // the friction multiplier, the elevation and Ek terms and the march
        // read the same numbers. These laws price the holdup from the local
        // superficial velocities and properties directly; the inclination
        // correction above is Beggs-Brill's own, so it does not apply, and
        // the liquid viscosity read here is the PVT one - the rheology
        // increment (5.76) reprices the friction term only, exactly as it
        // does on the Beggs-Brill path.
        hl = hagedornBrownHoldup(in.vsl, in.vsg, lambda, in.rhoL, in.muL,
                                 in.sigma, in.diameter, in.pPa);
    }
    else if (in.correlation == 2)
    {
        hl = dunsRosHoldup(in.vsl, in.vsg, lambda, in.rhoL, in.muL,
                           in.sigma, in.diameter);
    }
    else
    {
        // 5.155 (G4 family): the OLGAS-steady public mechanistic skeleton
        // REPLACES the horizontal + inclination block and nothing else: the
        // map, the friction multiplier, the elevation and Ek terms and the
        // march read the same numbers, and the map's own verdict (never a
        // forced regime - flowPattern is refused alongside any alternative
        // law) dispatches the branch. The liquid viscosity and surface
        // tension do not enter this law; the inclination physics is the
        // branch-local gravity term declared in olgasSteadyHoldup above.
        hl = olgasSteadyHoldup(in.vsl, in.vsg, lambda, in.rhoL, in.rhoG,
                               in.diameter, std::sin(theta), regime,
                               pat.transitionWeight);
    }
    hl = std::min(std::max(hl, 0.0), 1.0);
    out.holdup = hl;

    const double rhoM = in.rhoL * lambda + in.rhoG * (1.0 - lambda);
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

    // the published single-branch gradient at one holdup (BB 1973): the
    // two-phase friction multiplier S on y = lambda / HL^2, the elevation
    // head on the slippage (in-situ) density, the Darcy-Weisbach friction on
    // the no-slip mixture density and the Ek acceleration correction, all
    // evaluated on THIS branch's holdup. Both transition branches and every
    // single-regime path price through it, so the branches cannot drift
    // from the regimes they blend with.
    const double sinTheta = std::sin(theta);
    const double pEk = std::max(in.pPa, 1000.0);
    auto branchGradient = [&](double hlB, double& fTpB)
    {
        hlB = std::min(std::max(hlB, 0.0), 1.0);
        const double rhoS = in.rhoL * hlB + in.rhoG * (1.0 - hlB);
        double s;
        if (lambda < 1e-9 || hlB < 1e-9)
        {
            s = 0.0; // single-phase degenerate: no slippage correction
        }
        else
        {
            const double y = lambda / (hlB * hlB);
            const double ly = std::log(y);
            if (y > 1.0 && y < 1.2) s = std::log(2.2 * y - 1.2);
            else s = ly / (ly * (ly * (0.01853 * ly * ly - 0.8725) + 3.182) - 0.0523);
            s = std::min(s, 7.0);
        }
        fTpB = fNs * std::exp(s);
        const double ele = kGravity * sinTheta * rhoS;
        const double fric = fTpB * rhoM * vm * vm / (2.0 * in.diameter);
        const double ek = in.vsg * vm * rhoS / pEk;
        return (ele + fric) / std::max(1.0 - ek, 0.1);
    };

    double fTp = 0.0;
    if (regime == 1 && in.correlation == 0)
    {
        // the published transition gradient: (dp/dz)_t = w (dp/dz)_seg
        // + (1 - w) (dp/dz)_int, each branch a COMPLETE gradient (its own
        // holdup, inclination set, S multiplier, slippage density and Ek).
        // The reported fTp is the same blend, which is exactly the
        // multiplier consistent with the blended friction term.
        const double w = pat.regime == 1
            ? pat.transitionWeight
            : (pat.l3 - fr) / (pat.l3 - pat.l2);
        double fSeg = 0.0, fInt = 0.0;
        const double dpSeg = branchGradient(hlSeg, fSeg);
        const double dpInt = branchGradient(hlInt, fInt);
        out.dpPerMeter = w * dpSeg + (1.0 - w) * dpInt;
        fTp = w * fSeg + (1.0 - w) * fInt;
    }
    else
    {
        out.dpPerMeter = branchGradient(hl, fTp);
    }
    out.fTp = fTp;
    return out;
}

BeggsBrillMarchOut marchBeggsBrillPipe(const pvt::BlackOilFluid& fluid,
                                       double qoStd, double pInPa,
                                       double length, double diameter,
                                       double roughness, double sinTheta,
                                       int segments,
                                       const Rheology& liquidRheology,
                                       int forcedRegime,
                                       const std::function<double(double)>* sigmaOfP,
                                       int correlation,
                                       const SteadyThermal* thermal)
{
    BeggsBrillMarchOut out;
    if (thermal && thermal->segmentU
        && static_cast<int>(thermal->segmentU->size()) != std::max(segments, 1))
        cfFatal("beggsBrill: SteadyThermal segmentU needs exactly one entry "
                "per march segment");
    const double pFloor = std::max(0.02 * pInPa, 1000.0);
    const double ds = length / static_cast<double>(std::max(segments, 1));
    const double area = 0.25 * M_PI * diameter * diameter;
    const double wor =
        fluid.waterCut > 0 && fluid.waterCut < 0.999
            ? fluid.waterCut / (1.0 - fluid.waterCut)
            : (fluid.waterCut >= 0.999 ? 999.0 : 0.0);
    out.massRate = qoStd * pvt::blackOilMassPerOilSm3(fluid);

    double p = pInPa;
    // 5.144 (G6): the passive temperature state; 0 when no thermal march
    // was declared, and nothing below reads it in that case.
    double t = thermal ? thermal->temperatureIn : 0.0;
    double lastHoldup = 0, lastLambda = 0, lastVm = 0, lastFTp = 0;
    int lastRegime = 0;
    int lastRegimeMap = 0;
    double lastL1 = 0, lastL2 = 0, lastL3 = 0, lastL4 = 0;
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
        // D10/D7: a caller may hand the march a per-pressure surface-tension
        // provider (e.g. live-flash Weinaug-Katz). Absent, this is exactly
        // the original constant-sigma black-oil path, bit for bit.
        in.sigma = (sigmaOfP && *sigmaOfP) ? (*sigmaOfP)(p)
                                           : fluid.surfaceTension;
        in.vsl = vsl;
        in.vsg = vsg;
        in.pPa = p;
        in.diameter = diameter;
        in.roughness = roughness;
        in.angleDeg = angleDeg;
        in.liquidRheology = liquidRheology;
        in.forcedRegime = forcedRegime;
        in.correlation = correlation;

        BeggsBrillOut g;
        double staticRhoM = 0.0; // only the dead-line branch fills this
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
            staticRhoM = rhoM;
            BeggsBrillOut s;
            s.dpPerMeter = kGravity * sinTheta * rhoM;
            s.holdup = st.bo / vT;
            s.lambdaL = st.bo / vT;
            s.regime = 3;
            s.regimeFromMap = 3; // no flow, no dispatch to override
            s.fNs = 0;
            s.fTp = 0;
            g = s;
        }
        // 5.158 / 5.159: the per-segment readout profile, pushed before the
        // pressure and temperature of the segment are stepped. Pure report:
        // nothing on the march's own numbers reads it back.
        out.segPressure.push_back(p);
        out.segHoldup.push_back(g.holdup);
        out.segLambda.push_back(g.lambdaL);
        out.segVm.push_back(vsl + vsg);
        out.segRhoM.push_back(vsl + vsg > 1e-9
                                  ? rhoL * g.lambdaL
                                        + st.rhoGas * (1.0 - g.lambdaL)
                                  : staticRhoM);
        out.segRegime.push_back(g.regime);
        out.segTemperature.push_back(t);
        // 5.188 (G6b): absolute position of this point, filled only on a
        // thermal march so the isothermal path keeps an empty vector.
        if (thermal) out.segX.push_back(static_cast<double>(i) * ds);
        lastVm = vsl + vsg;
        lastFTp = g.fTp;
        lastRegime = g.regime;
        lastRegimeMap = g.regimeFromMap;
        lastL1 = g.l1;
        lastL2 = g.l2;
        lastL3 = g.l3;
        lastL4 = g.l4;

        out.holdupMax = std::max(out.holdupMax, g.holdup);
        out.regimeMin = std::min(out.regimeMin, g.regime);
        out.regimeMax = std::max(out.regimeMax, g.regime);
        lastHoldup = g.holdup;
        lastLambda = g.lambdaL;
        const double pOld = p;
        p -= g.dpPerMeter * ds;
        const bool hitFloor = !(p > pFloor);
        if (hitFloor) p = pFloor;
        // 5.144 (G6): the passive thermal step on the ACTUAL segment dp
        // (the floor clamp included), after the hydraulics of the segment
        // are done. With thermal null this is the isothermal march, bit
        // for bit - the pressure numbers above never saw t.
        if (thermal)
            t = steadyThermalAdvance(*thermal, t, pOld, p, i, ds, sinTheta,
                                     out.massRate);
        if (hitFloor)
        {
            out.pOutlet = p;
            out.dpTotal = pInPa - p;
            out.holdupOutlet = lastHoldup;
            out.lambdaOutlet = lastLambda;
            out.vmOutlet = lastVm;
            out.fTpOutlet = lastFTp;
            out.regimeOutlet = lastRegime;
            out.regimeMapOutlet = lastRegimeMap;
            out.nfrOutlet = lastVm * lastVm / (kGravity * diameter);
            out.l1Outlet = lastL1;
            out.l2Outlet = lastL2;
            out.l3Outlet = lastL3;
            out.l4Outlet = lastL4;
            out.temperature = t;
            out.ok = false;
            return out;
        }
    }
    out.temperature = t;
    out.pOutlet = p;
    out.dpTotal = pInPa - p;
    out.holdupOutlet = lastHoldup;
    out.lambdaOutlet = lastLambda;
    out.vmOutlet = lastVm;
    out.fTpOutlet = lastFTp;
    out.regimeOutlet = lastRegime;
    out.regimeMapOutlet = lastRegimeMap;
    out.nfrOutlet = lastVm * lastVm / (kGravity * diameter);
    out.l1Outlet = lastL1;
    out.l2Outlet = lastL2;
    out.l3Outlet = lastL3;
    out.l4Outlet = lastL4;
    out.ok = true;
    return out;
}

} // namespace cf
