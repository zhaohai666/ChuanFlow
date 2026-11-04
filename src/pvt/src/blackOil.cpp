#include "cf/pvt/BlackOil.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"

namespace cf::pvt
{

namespace
{

constexpr double kPsiToPa = 6894.757293168361;
constexpr double kScfPerStbPerSm3Sm3 = 5.614583333333333; // = ft3/bbl
constexpr double kMolWeightAir = 28.9647;        // g/mol
constexpr double kGasConstant = 8.314462618;     // J/(mol K)
constexpr double kStdPressure = 101325.0;        // Pa
constexpr double kStdTemperature = 288.15;       // K

double rankineFromK(double K) { return K * 1.8; }
double fahrenheitFromK(double K) { return K * 1.8 - 459.67; }
double apiFromGamma(double gammaO) { return 141.5 / gammaO - 131.5; }

} // namespace

double standingBubblePointPsia(double gasGravity, double api, double TF,
                               double rsScfStb)
{
    // pb = 18.2 [(Rs/gg)^0.83 * 10^a - 1.4], a = 0.00091 T(F) - 0.0125 API
    // (whitson wiki bopvt; Ahmed eq. 2-76/2-77). The original chart fit used
    // T - 460 for the rankine conversion, so TF itself enters directly.
    const double a = 0.00091 * TF - 0.0125 * api;
    const double pb = 18.2 * (std::pow(rsScfStb / gasGravity, 0.83)
                                  * std::pow(10.0, a)
                              - 1.4);
    return std::max(pb, 0.0);
}

double vazquezBeggsRsScfStb(double pPsia, double gasGravity, double api,
                            double TF)
{
    const double c1 = api <= 30.0 ? 0.0362 : 0.0178;
    const double c2 = api <= 30.0 ? 1.0937 : 1.187;
    const double c3 = api <= 30.0 ? 25.724 : 23.931;
    const double TR = TF + 459.67;
    return c1 * gasGravity * std::pow(std::max(pPsia, 0.0), c2)
         * std::exp(c3 * api / TR);
}

double vazquezBeggsBoRbStb(double pPsia, double rsScfStb, double pbPsia,
                           double gasGravity, double api, double TF)
{
    const bool heavy = api > 30.0;
    const double c1 = heavy ? 4.670e-4 : 4.677e-4;
    const double c2 = heavy ? 1.100e-5 : 1.751e-5;
    const double c3 = heavy ? 1.337e-9 : -1.811e-8;
    const double tTerm = (TF - 60.0) * api / gasGravity;
    const double bob = 1.0 + c1 * rsScfStb + c2 * tTerm + c3 * rsScfStb * tTerm;
    if (pPsia <= pbPsia) return bob;
    // undersaturated: VB compressibility (1/psi), rs at the bubble point.
    // The correlation divides by pb; a degenerate bubble point (dead oil,
    // Rs = 0 gives pb = 0) would blow co up by orders of magnitude, so
    // floor the denominator into the correlation's fitted range.
    const double co =
        (-1433.0 + 5.0 * rsScfStb + 17.2 * TF - 1180.0 * gasGravity
         + 12.61 * api)
        / (1.0e5 * std::max(pbPsia, 300.0));
    return bob * std::exp(co * (pbPsia - pPsia)); // co < 0, p > pb shrinks
}

double beggsRobinsonDeadOilCp(double api, double TF)
{
    const double x = std::pow(10.0, 3.0324 - 0.02023 * api) * std::pow(TF, -1.163);
    return std::pow(10.0, x) - 1.0;
}

double beggsRobinsonLiveOilCp(double deadOilCp, double rsScfStb)
{
    const double a = 10.715 * std::pow(rsScfStb + 100.0, -0.515);
    const double b = 5.44 * std::pow(rsScfStb + 150.0, -0.338);
    return a * std::pow(deadOilCp, b);
}

double leeGasViscosityCp(double molecularWeight, double TR, double rhoGPerCm3)
{
    const double k = (9.4 + 0.02 * molecularWeight) * std::pow(TR, 1.5)
                   / (209.0 + 19.0 * molecularWeight + TR);
    const double x = 3.5 + 986.0 / TR + 0.01 * molecularWeight;
    const double y = 2.4 - 0.2 * x;
    return 1e-4 * k * std::exp(x * std::pow(rhoGPerCm3, y));
}

void standingPseudoCritical(double gasGravity, double& ppcPsia,
                            double& tpcRankine)
{
    if (gasGravity < 0.75)
    {
        tpcRankine = 168.0 + 325.0 * gasGravity - 12.5 * gasGravity * gasGravity;
        ppcPsia = 677.0 + 15.0 * gasGravity - 37.5 * gasGravity * gasGravity;
    }
    else
    {
        tpcRankine = 187.0 + 330.0 * gasGravity - 71.5 * gasGravity * gasGravity;
        ppcPsia = 706.0 + 51.0 * gasGravity - 11.1 * gasGravity * gasGravity;
    }
}

double dranchukAbouKassemZ(double pPsia, double TR, double ppcPsia,
                           double tpcRankine)
{
    static const double A[11] = {0.3265,  -1.0700, -0.5339,  0.01569,
                                 -0.05165, 0.5475, -0.7361,  0.1844,
                                 0.1056,  0.6134,  0.7210};
    const double pr = std::max(pPsia, 1e-6) / ppcPsia;
    const double tr = TR / tpcRankine;
    const double t1 = A[0] + A[1] / tr + A[2] / (tr * tr * tr)
                    + A[3] / (tr * tr * tr * tr) + A[4] / (tr * tr * tr * tr * tr);
    const double t2 = A[5] + A[6] / tr + A[7] / (tr * tr);
    // DAK residual as a function of Z alone: with rho = 0.27 Pr / (Z Tr) the
    // EOS reads Z = 1 + t1 rho + t2 rho^2 - A9(...) rho^5
    //                  + A10 (1 + A11 rho^2) rho^2/Tr^3 exp(-A11 rho^2).
    // The fixed-point residual is continuous on Z in [0.05, 3] with a
    // guaranteed sign change (Z(0.05) >> 0.05 from the rho^5 term, Z(3) ~ 1),
    // so bracketed bisection is unconditionally convergent - the raw Newton
    // iteration is not (the residual is stiff near rho ~ 1).
    auto zOfRho = [&](double rho)
    {
        const double r2 = rho * rho;
        return 1.0 + t1 * rho + t2 * r2
             - A[8] * (A[6] / tr + A[7] / (tr * tr)) * std::pow(rho, 5)
             + A[9] * (1.0 + A[10] * r2) * (r2 / (tr * tr * tr))
                   * std::exp(-A[10] * r2);
    };
    auto residual = [&](double z)
    {
        return zOfRho(0.27 * pr / (z * tr)) - z;
    };
    double lo = 0.05, hi = 3.0;
    if (residual(lo) <= 0.0) return lo;
    if (residual(hi) >= 0.0) return hi;
    for (int it = 0; it < 80; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        if (residual(mid) > 0.0) lo = mid;
        else hi = mid;
        if (hi - lo < 1e-13 * hi) break;
    }
    return 0.5 * (lo + hi);
}

double colebrookDarcyFriction(double reynolds, double relativeRoughness)
{
    if (reynolds < 2040.0) return 64.0 / std::max(reynolds, 1e-12);
    const double eD = std::max(relativeRoughness, 0.0);
    // Jain (1976) explicit seed, then fixed-point on the Colebrook-White
    // form y = -2 log10(eD/3.7 + 2.51 y / Re) with y = 1/sqrt(f); the map is
    // contractive for Re > 2040.
    double y = 1.14 - 2.0 * std::log10(eD + 21.25 / std::pow(reynolds, 0.9));
    y = std::max(y, 1.0);
    for (int it = 0; it < 100; ++it)
    {
        const double yNew =
            -2.0 * std::log10(eD / 3.7 + 2.51 * y / reynolds);
        const double err = std::abs(yNew - y) / yNew;
        y = yNew;
        if (err < 1e-14) break;
    }
    return 1.0 / (y * y);
}

double blackOilMassPerOilSm3(const BlackOilFluid& fluid)
{
    const double mAir = kMolWeightAir / 1000.0; // kg/mol
    const double gasMolar = fluid.gasGravity * mAir;
    const double rhoGasStd = kStdPressure * gasMolar / (kGasConstant * kStdTemperature);
    const double wor = fluid.waterCut > 0
                         ? fluid.waterCut / (1.0 - fluid.waterCut)
                         : 0.0;
    return fluid.oilDensityStd + fluid.gorStd * rhoGasStd + wor * fluid.waterDensity;
}

void blackOilLiquid(const BlackOilFluid& fluid, const BlackOilState& st,
                    double& rhoL, double& muL)
{
    const double vO = st.bo;
    const double vW = fluid.waterCut > 0
                        ? (fluid.waterCut / (1.0 - fluid.waterCut))
                        : 0.0;
    const double vT = vO + vW;
    rhoL = (st.rhoOil * vO + fluid.waterDensity * vW) / vT;
    const double fo = vO / vT;
    muL = std::exp(fo * std::log(st.muOil)
                 + (1.0 - fo) * std::log(fluid.waterViscosity));
}

BlackOilState evaluateBlackOil(const BlackOilFluid& fluid, double pPa)
{
    BlackOilState st;
    const double pPsia = std::max(pPa, 1.0) / kPsiToPa;
    const double TF = fahrenheitFromK(fluid.temperature);
    const double TR = rankineFromK(fluid.temperature);
    const double gammaO = fluid.oilDensityStd / 999.014;
    const double api = apiFromGamma(gammaO);
    const double gorScf = fluid.gorStd * kScfPerStbPerSm3Sm3;

    st.pb = standingBubblePointPsia(fluid.gasGravity, api, TF, gorScf) * kPsiToPa;

    if (pPa <= st.pb)
    {
        st.rs = std::min(vazquezBeggsRsScfStb(pPsia, fluid.gasGravity, api, TF)
                             / kScfPerStbPerSm3Sm3,
                         fluid.gorStd);
        st.bo = vazquezBeggsBoRbStb(pPsia, st.rs * kScfPerStbPerSm3Sm3,
                                    st.pb / kPsiToPa, fluid.gasGravity, api, TF);
        st.freeGas = fluid.gorStd - st.rs;
    }
    else
    {
        st.rs = fluid.gorStd;
        st.bo = vazquezBeggsBoRbStb(pPsia, gorScf, st.pb / kPsiToPa,
                                    fluid.gasGravity, api, TF);
        st.freeGas = 0.0;
    }

    // free-gas PVT: DAK Z on Standing pseudo-criticals
    double ppc = 0, tpc = 0;
    standingPseudoCritical(fluid.gasGravity, ppc, tpc);
    st.z = dranchukAbouKassemZ(pPsia, TR, ppc, tpc);
    st.bg = (kStdPressure / pPa) * (fluid.temperature / kStdTemperature) * st.z;

    const double mAir = kMolWeightAir / 1000.0;
    const double gasMolar = fluid.gasGravity * mAir;
    st.rhoGas = pPa * gasMolar / (st.z * kGasConstant * fluid.temperature);
    // saturated oil density from the mass balance of one STB:
    // rho_o [lb/ft3] = (350 gamma_o + 0.0764 gamma_g Rs) / (5.615 Bo)
    // (Ahmed eq. 2-63: 350 lb/STB surface oil, 0.0764 lb/scf gas, 5.615
    // ft3/bbl converts lb/rb to lb/ft3)
    const double lbPerFt3ToKgPerM3 = 16.0184634;
    const double gammaG = fluid.gasGravity;
    const double rhoOilLbFt3 =
        (350.0 * gammaO + 0.0764 * gammaG * st.rs * kScfPerStbPerSm3Sm3)
        / (kScfPerStbPerSm3Sm3 * st.bo);
    st.rhoOil = rhoOilLbFt3 * lbPerFt3ToKgPerM3;

    st.rhoWater = fluid.waterDensity;

    const double muDead =
        beggsRobinsonDeadOilCp(api, TF); // cp
    const double muLive =
        beggsRobinsonLiveOilCp(muDead, st.rs * kScfPerStbPerSm3Sm3);
    st.muOil = muLive * 1e-3; // cp -> Pa s

    const double rhoGasGcm3 = st.rhoGas / 1000.0;
    st.muGas = leeGasViscosityCp(gasMolar * 1000.0, TR, rhoGasGcm3) * 1e-3;

    return st;
}

} // namespace cf::pvt
