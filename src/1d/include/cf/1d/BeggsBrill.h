#pragma once

#include "cf/pvt/BlackOil.h"
#include "cf/common/Rheology.h"

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
    // 5.76 (A5 increment 2): the LIQUID phase's generalized-Newtonian law.
    // Only the liquid's share of the no-slip mixture viscosity is repriced
    // (on the mixture wall rate gdot = 8 vm / D, with the Rabinowitsch K'),
    // leaving the published holdup, flow-pattern and inclination
    // correlations -- none of which reads a viscosity -- exactly as they
    // are. The default is the undeclared fluid, bit for bit.
    Rheology liquidRheology;
    // 5.88 (D4 increment 2): the WRITABLE side of the flow-pattern dispatch
    // surface. The map in beggsBrillPattern is the single source of truth and
    // is ALWAYS evaluated -- its four boundaries and its own verdict stay
    // readable on the output -- but a non-negative value here selects which
    // regime's holdup law and inclination coefficients the gradient actually
    // uses (0 segregated, 1 transition, 2 intermittent, 3 distributed).
    // Beggs-Brill is a correlation, and field observation routinely knows the
    // pattern better than the map does, which is why the simulators expose this
    // override. -1 (the default) means "ask the map", and in that case every
    // number below is exactly what 5.85 produced, bit for bit.
    int forcedRegime = -1;
};

struct BeggsBrillOut
{
    double dpPerMeter = 0; // Pa/m pressure DROP along the flow direction
    double holdup = 0;     // liquid volume fraction HL
    int regime = 0;        // the regime the gradient was priced with: the
                           // forced one when forcedRegime >= 0, else
                           // regimeFromMap below
    double lambdaL = 0;    // no-slip liquid content
    double fNs = 0;        // no-slip Darcy friction factor
    double fTp = 0;        // two-phase Darcy friction factor
    // 5.85 (D4 increment 1): the flow-pattern dispatch surface read out at
    // this segment -- the four L boundaries at the segment's lambda and the
    // transition blend weight (nonzero only in regime 1). The gradient
    // dispatches its holdup through beggsBrillPattern, so these are the
    // numbers that selected the holdup, not a second evaluation.
    double l1 = 0, l2 = 0, l3 = 0, l4 = 0;
    double transitionWeight = 0;
    // 5.88 (D4 increment 2): the map's OWN verdict at this segment, kept
    // whether or not a forced regime overrode it. The pair (regime,
    // regimeFromMap) is how a caller sees the disagreement it just asked for;
    // `transitionWeight` is the blend weight actually used, which for a forced
    // transition is the published (L3 - Nfr)/(L3 - L2) held to [0, 1] so a
    // forced blend stays convex outside the band the map draws it in.
    int regimeFromMap = 0;
};

/// 5.85 (D4 increment 1): the Beggs-Brill flow-pattern map as a standalone
/// dispatch surface. Given the no-slip liquid content lambda and the Froude
/// number Nfr = vm^2/(g D), returns the regime (0 segregated, 1 transition,
/// 2 intermittent, 3 distributed), the four published boundaries
///   L1 = 316 lambda^0.302,      L2 = 0.0009252 lambda^-2.4684,
///   L3 = 0.1 lambda^-1.4516,    L4 = 0.5 lambda^-6.738,
/// and the transition blend weight w = (L3 - Nfr)/(L3 - L2) (regime 1 only).
/// beggsBrillDpDs dispatches through this function -- it is the single source
/// of truth for the map, queryable independently of any flow solve.
struct BeggsBrillPattern
{
    int regime = 0;
    double l1 = 0, l2 = 0, l3 = 0, l4 = 0;
    double transitionWeight = 0;
};

BeggsBrillPattern beggsBrillPattern(double lambda, double nfr);

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
    // 5.85 (D4 increment 1): the dispatch surface at the OUTLET segment --
    // its regime, no-slip content, Froude number and the four L boundaries.
    // regimeMin/Max stay the along-pipe envelope; these describe the last
    // segment alone, so a caller can read the map where the march ended.
    int regimeOutlet = 0;
    double nfrOutlet = 0;
    double l1Outlet = 0, l2Outlet = 0, l3Outlet = 0, l4Outlet = 0;
    // 5.88 (D4 increment 2): what the map itself said at the outlet segment.
    // Equal to regimeOutlet unless a forced regime overrode it, in which case
    // the difference is the override's footprint on the readout.
    int regimeMapOutlet = 0;
    bool ok = false;         // false when the march hit the pressure floor
};

BeggsBrillMarchOut marchBeggsBrillPipe(const pvt::BlackOilFluid& fluid,
                                       double qoStd, double pInPa,
                                       double length, double diameter,
                                       double roughness, double sinTheta,
                                       int segments = 24,
                                       const Rheology& liquidRheology = {},
                                       int forcedRegime = -1);

} // namespace cf
