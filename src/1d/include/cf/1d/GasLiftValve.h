#pragma once

namespace cf
{

/// An injection-pressure-operated (IPO) gas-lift valve: the force balance on
/// the bellows that decides whether the port is open at all, and the isentropic
/// nozzle that prices how much gas passes once it is.
///
/// The mechanics are the bellows/port force balance, not a fit. With the
/// injection (casing) pressure Pc on the annular face of the bellows and the
/// production (tubing) pressure Pt on the port, the two sides of the stem are
///     F_open  = Pc (Ab - Ap) + Pt Ap
///     F_close = Pd Ab + St (Ab - Ap)
/// where Pd is the nitrogen dome charge and St the spring preload expressed as
/// a pressure over the sensing area. Dividing by Ab and writing R = Ap/Ab is
/// Guo, Lyons & Ghalambor, "Petroleum Production Engineering" (2007) sec. 13.5,
/// Eqs (13.42) and (13.43) - the same pair API RP 11V6 and the CAMCO/Unocal
/// valve manuals carry:
///     Pvo = Pd/(1-R) + St - R/(1-R) Pt      (valve opens as Pc rises)
///     Pvc = Pd + St (1-R)                   (valve closes as Pc falls)
///     Pvo - Pvc = R [ (Pd - Pt)/(1-R) + St ]
/// which is the band the valve can sit in without a definite answer (the
/// "spread"), and R/(1-R) Pt is the tubing pressure's own contribution to it
/// (the production pressure effect, PPE). The band is open while the dome plus
/// spring outranks the tubing's lift on the port, Pt < Pd + St(1-R); past that
/// the spread inverts, which is the valve's "runs open" state - the production
/// pressure alone holds the port off its seat, and no casing pressure above the
/// (lowered) opening pressure can be asked to seat it. All pressures ABSOLUTE.
struct GasLiftValve
{
    double domePressure = 0;    // Pa absolute, nitrogen charge at valve T (Pd)
    double areaRatio = 0;       // R = Ap/Ab, in (0,1)
    double springPressure = 0;  // Pa, spring term St (0 = no spring)
    double portArea = 0;        // m2, the flow bore the gas passes through (Ap)
    double dischargeCoeff = 0;  // Cd, the port's velocity/side-contraction loss
    double kappa = 0;           // cp/cv of the lift gas (k > 1)
};

/// The two thresholds, the spread, the PPE, and the test-rack opening: the
/// valve's own settings, which depend on the tubing pressure only through the
/// force balance and not at all on whether gas is flowing.
struct GasLiftValveThresholds
{
    double openingPressure = 0; // Pvo, Pa absolute
    double closingPressure = 0; // Pvc, Pa absolute
    double spread = 0;          // Pvo - Pvc, Pa
    double productionEffect = 0; // R/(1-R) Pt, Pa
    double testRackOpening = 0; // Pvo at Pt = 101325 Pa (14.7 psia), Pa absolute
};

GasLiftValveThresholds gasLiftValveThresholds(const GasLiftValve& valve,
                                              double tubingPressure);

/// The isentropic nozzle mass rate (kg/s) through the port from the upstream
/// state (pUp, rhoUp) down to pDown, and whether the throat is sonic.
///
/// This is the compressible-orifice flux itself, so no published constant is
/// being trusted: with r = pDown/pUp and the upstream state an ideal gas
/// (pUp/rhoUp = Rs T fixed at the valve), integrating the isentropic relation
/// along the streamline and applying the energy equation gives the mass FLUX
///     G = sqrt( 2 rhoUp pUp k/(k-1) [ r^(2/k) - r^((k+1)/k) ] ),
///     mdot = Cd Ap G,
/// which has two properties the assertions lock. As r -> 1 the bracket is
/// ((k-1)/k)(1-r) + O((1-r)^2), so G -> sqrt(2 rhoUp dp): the incompressible
/// Chézy form (Bernoulli in velocity head), for every k - the correlation is
/// not a curve, it is Bernoulli plus the expansion the gas actually undergoes.
/// And the bracket's own maximum
/// sits at r* = (2/(k+1))^(k/(k-1)), the critical ratio - differentiating the
/// bracket gives 2 r^((1-k)/k) = k+1, and at k = 1.4 that is the textbook 0.528.
/// (The exponent 1/(k-1) belongs to the THROAT DENSITY ratio, (2/(k+1))^(1/(k-1))
/// = rho*/rho0; quoting it here would choke this gas at 0.628 instead of 0.546
/// and mis-price every passing arm by 1.6%.) Past r* no downstream pressure is
/// felt and the flux holds at G(r*), which is what `choked` reports. G is
/// continuous at r*, so the valve law has no step for Newton to bounce on.
struct GasLiftValveFlow
{
    double massRate = 0;       // kg/s through the port
    bool choked = false;       // pDown/pUp <= r*: sonic at the port
};

GasLiftValveFlow gasLiftValveFlow(const GasLiftValve& valve, double pUp,
                                  double rhoUp, double pDown);

} // namespace cf
