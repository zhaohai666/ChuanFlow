#pragma once

namespace cf
{

/// A subsurface rod (sucker-rod) pump: the force balance on the rod string at
/// the plunger, the elastic stroke it costs, and the card that records both.
///
/// Nothing here is a curve fit. Four statements make the whole element, and
/// each is either statics or Hooke's law:
///
///  1. The rod string's weight as the pumping unit feels it is its weight in
///     air minus the buoyancy of the fluid it displaces,
///         Wr_air = rho_r g A_r L,   buoyancy = rho_l g A_r L,
///         Wr_f   = (rho_r - rho_l) g A_r L,
///     because a vertical string in a hydrostatic field has a net upward
///     pressure force equal to the weight of the volume it occupies.
///
///  2. The fluid load is the pressure difference across the plunger times the
///     plunger area,
///         Fo = A_p (p_above - p_below),
///     where the two sides are the tubing pressure just above the pump and the
///     intake pressure just below it -- which is what the network solves, not
///     what an assumed gradient would guess. With both ends surfaced at the
///     same pressure, no friction in either column and the annulus full to the
///     surface, that difference is rho_l g L and the load is at its frictionless
///     ceiling Fo_max = rho_l g L A_p. The real load leaves that ceiling by four
///     named terms: an annulus fluid level below the surface and a casing
///     pressure both take weight off the plunger's underside (their sum is
///     reported as `upthrust`), and tubing friction and a wellhead backpressure
///     both add to it. So `upthrust` is positive on a well with a low level and
///     negative on a pump hauling a frictions column against a separator.
///
///  3. The tension the string carries varies along its length -- the top holds
///     what the bottom holds plus the string's own buoyed weight -- so the
///     elongation under a bottom load W is int_0^L (W + w'(L-y)) dy / (A_r E),
///     whose W-dependent part is exactly
///         delta(W) = W L / (A_r E).
///     The constant string-weight part only moves the datum, so the stroke
///     cost of picking the fluid up is
///         delta_f = Fo L / (A_r E).
///
///  4. Kinematics. Take the polished-rod position x from bottom dead centre and
///     measure the stretch against the no-fluid-load state. The plunger sits at
///     x_dp = x - delta(W), and while the load is transferring W = A_r E x / L,
///     so x_dp = 0: the first delta_f of polished-rod travel goes entirely into
///     stretching the string and the plunger does not move. Past that W is
///     constant at Fo and the plunger travels
///         S_eff = S - delta_f,
///     which is the effective stroke, and it is what the barrel delivers. When
///     delta_f reaches S the plunger never leaves bottom dead centre at all and
///     the pump delivers nothing however hard the pumping unit strokes.
///
/// Put the two load levels and that travel on one axis and the card is the
/// textbook parallelogram: vertices (0, W_min), (delta_f, W_max), (S, W_max),
/// (S - delta_f, W_min), with W_max - W_min = Fo. The downhole card -- load at
/// the plunger, so with the string's own buoyed weight taken out -- is the same
/// figure squared into a rectangle whose lines are 0 and Fo. And the area of
/// the parallelogram is Fo S_eff: the work the rods do on the fluid per stroke,
/// which the element then reports as the pump's hydraulic power at the network's
/// own pressure rise. That identity is the reason the four numbers above are
/// worth computing together rather than separately.
struct RodPump
{
    double plungerDiameter = 0;   // m, the barrel bore the plunger sweeps (A_p)
    double strokeLength = 0;      // m, polished-rod stroke S
    double strokesPerMinute = 0;  // spm
    double rodDiameter = 0;       // m, the rod string's gauge (A_r)
    double rodLength = 0;         // m, pump setting depth L
    double youngModulus = 0;      // Pa, the string's E
    double rodDensity = 0;        // kg/m3, the string's rho_r
};

/// The three states a rod pump's card can come back in. FULL is the
/// parallelogram. NOLOAD is the flat line at the buoyed rod weight: the
/// differential across the plunger does not oppose the upstroke (a packed,
/// over-pressured annulus, or a well that flows on its own), so the barrel
/// strokes without lifting a column and the string never stretches under fluid
/// load. STALLED is a pump whose fluid load stretches the string by the whole
/// stroke: the plunger stays at bottom dead centre and no volume is delivered.
enum class RodPumpState
{
    full = 0,
    noLoad = 1,
    stalled = 2
};

/// Everything one card states. The loads are newtons at the plunger's depth,
/// the travels are metres of polished-rod stroke, the rates are cubic metres
/// per second of barrel volume at the intake conditions and kilograms per
/// second of the same, and the power is the stroke work times the stroke rate.
struct RodPumpCard
{
    RodPumpState state = RodPumpState::full;

    double plungerArea = 0;    // A_p, m2
    double rodArea = 0;        // A_r, m2

    double rodWeightInAir = 0;   // rho_r g A_r L, N
    double buoyancy = 0;         // rho_l g A_r L, N
    double rodWeightInFluid = 0; // Wr_f, N
    double upthrust = 0;         // the part of the column the intake takes off
                                 // the plunger, N: Fo_max - Fo

    double fluidLoad = 0;        // Fo, N (0 in NOLOAD)
    double fluidLoadMax = 0;     // rho_l g L A_p, N

    double maxLoad = 0;          // W_max = Wr_f + Fo, N at the polished rod
    double minLoad = 0;          // W_min = Wr_f, N
    double downholeMaxLoad = 0;  // Fo, N at the plunger
    double downholeMinLoad = 0;  // 0, N
    double counterbalance = 0;   // (W_max + W_min)/2, N

    double stretch = 0;          // delta_f, m
    double stroke = 0;           // S, m
    double effectiveStroke = 0;  // S_eff, m (0 in STALLED)
    double efficiency = 0;       // S_eff / S

    double displacement = 0;     // A_p S_eff N, m3/s at intake conditions
    double massRate = 0;         // rho_l * displacement, kg/s
    double workPerCycle = 0;     // Fo * S_eff, J
    double hydraulicPower = 0;   // workPerCycle * strokes per second, W
};

/// The card, from the pump's own geometry, the fluid density it works in and
/// the two pressures the network puts at its ports. `pressureAbove` is the
/// tubing side (the plunger's discharge), `pressureBelow` the intake side.
/// Absolute pressures; only their difference is read.
RodPumpCard rodPumpCardOf(const RodPump& pump, double fluidDensity,
                          double pressureAbove, double pressureBelow);

/// The same card at an imposed pressure rise across the plunger, for callers
/// that hold the differential rather than the two port pressures.
RodPumpCard rodPumpCardOfLoad(const RodPump& pump, double fluidDensity,
                              double pressureRise);

} // namespace cf
