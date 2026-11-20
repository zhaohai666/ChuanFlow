#pragma once

#include "cf/1d/GasLiftValve.h"

namespace cf
{

/// Plunger gas lift (PGL): a BATCH lift point element. A free plunger runs up
/// the tubing, and when it arrives at the release valve at the setting depth
/// it trips the valve open for a declared window; an accumulator of gas at
/// surface (or in the casing annulus) then dumps a declared volume down the
/// injection line in one slug, and the valve seats again until the next
/// plunger. The cycle then repeats with a declared period -- the same
/// capture/release periodicity the terrain slug oscillator of 5.92 states for
/// dips, here with the gas inventory as the captured quantity.
///
/// Nothing new is fitted. Three statements, each exact or already priced:
///
///  1. The schedule is a declaration. Per cycle of length T the plunger holds
///     the valve open for a window of length w (0 < w <= T), and the
///     accumulator gives up a volume V_g quoted at a reference pressure p_ref.
///     The released-gas-rate schedule is therefore piecewise constant:
///     q_ref(t) = V_g / w for (t mod T) in [0, w) while the valve can pass,
///     0 otherwise, and its integral over one cycle is exactly V_g -- the
///     batch, not a rate, is what the operator sets. The average over the
///     cycle is V_g / T.
///
///  2. Whether the valve can pass is the 5.89 bellows force balance read at
///     depth: the casing (injection) pressure and the tubing pressure both
///     ride down their columns hydrostatically, Pc(z) = Pc_surf + rho_l g z
///     and Pt(z) = Pt_surf + rho_l g z, and gasLiftValveThresholds decides
///     open / closed / in-the-spread with the same convention as the steady
///     network (in the band, priced closed, and it says so).
///
///  3. How much passes is the 5.89 isentropic nozzle, capped against the
///     batch. The nozzle's choked-or-not mass rate mdot_nozzle is what the
///     port can hand over inside the window; the schedule demands
///     mdot_req = rho_ref V_g / w. The delivered batch is
///     min(mdot_req, mdot_nozzle) x w -- a valve too small to pass its
///     accumulator in time is reported VALVE-LIMITED with the shortfall as a
///     number, not quietly ignored. Downhole the slug's volume is Boyle's:
///     V(z) = V_ref p_ref / Pc(z), isothermal ideal gas at the valve (the
///     same assumption the 5.89 nozzle makes at its own port).
///
/// The plunger's own transit (its velocity, the liquid film it carries, the
/// cushion above it, the re-seating dynamics) is not modelled: the window is
/// declared, not solved. Gas compressibility beyond Boyle at valve
/// temperature, heat-up of the slug, and the two-phase response of the tubing
/// column to the slug are outside the element; the network's periodicity
/// consumption (feeding q(t) into a transient row) is a later increment.
struct PlungerGasLift
{
    GasLiftValve valve;             // the release valve at depth (5.89)
    double gasVolumePerCycle = 0;   // m3 of released gas at reference P (V_g)
    double referencePressure = 0;   // Pa abs the volume is quoted at (p_ref)
    double gasDensityAtReference = 0; // kg/m3 at p_ref (rho_ref)
    double cyclePeriod = 0;         // s, T (plunger round trip + charge)
    double openWindow = 0;          // s, w the valve is tripped open per cycle
    double releaseDepth = 0;        // m, vertical depth z of the valve
    double casingSurfacePressure = 0; // Pa abs at the injection column top
    double tubingSurfacePressure = 0; // Pa abs at the production column top
    double liquidDensity = 0;       // kg/m3, the annulus liquid gradient
};

/// 0 OPEN: valve passes the whole batch inside the window.
/// 1 VALVE-LIMITED: open, but the nozzle cannot pass the batch in w.
/// 2 CLOSED: the injection pressure sits at or below the closing threshold.
/// 3 IN-SPREAD: between Pvc and Pvo -- the steady statement has no memory of
///   where the plunger came from, so it is priced closed and says so (the
///   5.89 convention).
enum class PlungerGasLiftState
{
    open = 0,
    valveLimited = 1,
    closed = 2,
    inSpread = 3
};

struct PlungerGasLiftSchedule
{
    PlungerGasLiftState state = PlungerGasLiftState::open;

    double injectionPressureAtDepth = 0; // Pc(z), Pa abs
    double tubingPressureAtDepth = 0;    // Pt(z), Pa abs
    double openingPressure = 0;          // Pvo at Pt(z), Pa abs
    double closingPressure = 0;          // Pvc, Pa abs
    bool choked = false;                 // the nozzle's own verdict

    double windowRateReference = 0;   // m3/s at p_ref while the window is open
    double windowRateDownhole = 0;    // m3/s at Pc(z) while open (Boyle)
    double averageRateReference = 0;  // m3/s at p_ref over the whole cycle
    double massRateInWindow = 0;      // kg/s while open
    double requiredMassRate = 0;      // kg/s the batch demands inside w
    double nozzleMassRateCapacity = 0;// kg/s the port can pass (5.89 nozzle)
    double deliveredVolumeReference = 0; // m3 at p_ref per cycle (<= V_g)
    double deliveredVolumeDownhole = 0;  // m3 at Pc(z) per cycle
};

/// Is the release window open at time t (the schedule's window, before the
/// valve's own verdict)? t is measured from the start of a cycle.
bool plungerGasLiftWindowOpenAt(const PlungerGasLift& pgl, double t);

/// The per-cycle statement: pressures at depth, the force balance's verdict,
/// the nozzle's capacity, and the batch that actually arrives.
PlungerGasLiftSchedule plungerGasLiftSolve(const PlungerGasLift& pgl);

} // namespace cf
