#pragma once

#include "cf/pvt/BlackOil.h"

namespace cf
{

/// Beggs & Brill (1973) two-phase pressure-gradient correlation, implemented
/// from the published defining equations and cross-checked coefficient set:
///   - flow-pattern map with the L1..L4 boundaries (regimes: segregated,
///     transition, intermittent, distributed),
///   - horizontal liquid holdup HL(0) = a lambda^b / Nfr^c per regime
///     (segregated 0.98/0.4846/0.0868, intermittent 0.845/0.5351/0.0173,
///     distributed 1.065/0.5824/0.0609), transition regime interpolated,
///   - inclination correction psi = 1 + C(sin 1.8th - sin^3 1.8th / 3),
///     C = (1-lambda) ln(d lambda^e Nlv^f Nfr^h) with the published
///     uphill/downhill coefficient sets (distributed uphill: C = 0),
///   - two-phase friction factor f_tp = f_ns exp(S) on y = lambda/HL^2 with
///     the S polynomial and the 1 < y < 1.2 logarithmic branch,
///   - pressure gradient (rho_s g sin(th) + f_tp rho_m vm^2 / 2D)
///     / (1 - Ek), Ek = vsg vm rho_s / p (acceleration, slippage density in
///     the elevation term, no-slip density in the friction term).
/// The correlation operates on in-situ superficial velocities; the caller
/// derives them from the black-oil stream (cf/pvt/BlackOil.h).
struct BeggsBrillIn
{
    double rhoL = 0;         // kg/m3
    double rhoG = 0;         // kg/m3
    double muL = 0;          // Pa s
    double muG = 0;          // Pa s
    double sigma = 0.03;     // N/m
    double vsl = 0;          // m/s superficial liquid
    double vsg = 0;          // m/s superficial gas
    double pPa = 0;          // absolute pressure (acceleration term)
    double diameter = 0;     // m
    double roughness = 0;    // m
    double angleDeg = 0;     // from horizontal, positive = uphill
};

struct BeggsBrillOut
{
    double dpPerMeter = 0; // Pa/m pressure DROP along the flow direction
    double holdup = 0;     // liquid volume fraction HL
    int regime = 0;        // 0 segregated, 1 transition, 2 intermittent,
                           // 3 distributed
    double lambdaL = 0;    // no-slip liquid content
    double fNs = 0;        // no-slip Darcy friction factor
    double fTp = 0;        // two-phase Darcy friction factor
};

BeggsBrillOut beggsBrillDpDs(const BeggsBrillIn& in);

/// March a black-oil stream through one pipe: the stream's PVT state and
/// therefore the in-situ velocities evolve along the pipe, so the Beggs-Brill
/// gradient is evaluated per segment at the segment inlet pressure and the
/// march steps p - dpPerMeter ds from the upstream (higher-pressure) end.
/// sinTheta is the signed elevation ratio of the flow direction; a zero rate
/// degenerates to the static column on the local no-slip mixture density.
struct BeggsBrillMarchOut
{
    double pOutlet = 0;      // Pa
    double dpTotal = 0;      // Pa, pIn - pOutlet
    double massRate = 0;     // kg/s, invariant along the march
    double holdupMax = 0;    // largest HL encountered
    double holdupOutlet = 0;
    double lambdaOutlet = 0;
    double vmOutlet = 0;     // m/s, in-situ mixture velocity of the last
                             // segment
    double fTpOutlet = 0;    // two-phase Darcy friction of the last segment
    int regimeMin = 9;
    int regimeMax = -1;
    bool ok = false;         // false when the march hit the pressure floor
};

BeggsBrillMarchOut marchBeggsBrillPipe(const pvt::BlackOilFluid& fluid,
                                       double qoStd, double pInPa,
                                       double length, double diameter,
                                       double roughness, double sinTheta,
                                       int segments = 24);

} // namespace cf
