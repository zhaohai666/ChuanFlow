#pragma once

namespace cf
{

/// A progressive cavity pump (PCP / Moineau pump) as a point element: a
/// positive-displacement machine whose rate is kinematics minus a declared
/// leakage line, priced at the pressure rise its two ports actually differ by.
///
/// Two statements, no curve fit:
///
///  1. Displacement is geometry times rotation. The rotor's planetary motion
///     advances every sealed cavity one step per revolution, so each
///     revolution transfers `cavitiesPerRev` cavity volumes from intake to
///     discharge -- the Moineau kinematic count (cf. Takacs, "Progressive
///     Cavity Pumps"; the same Q_th = V_per-rev x N form the API PCP sizing
///     literature carries). With the cavity volume declared per cavity and the
///     rpm on the drive,
///         q_th = cavitiesPerRev * cavityVolume * rpm / 60.
///
///  2. Slip is a line. The pressure rise drives a backflow through the
///     annular clearances between adjacent cavities; that gap flow is a thin
///     laminar annulus, so its rate is linear in the differential by
///     construction (Poiseuille), and the published PCP characteristic is
///     exactly the straight line
///         q = q_th - slip * dp,
///     with `slip` (m3/s per Pa) declared, not fitted. dp <= 0 clamps to 0 the
///     way the rod pump's load does (NO-LOAD): the cavities still sweep, but
///     the differential does not oppose the transfer, so nothing leaks back
///     and the pump does no work.
///
/// The hydraulic power is the product the machine actually delivers to the
/// fluid, q * dp -- the area of the (q, dp) rectangle, not the swept volume's
/// bookkeeping. Shaft power (stator elastomer friction, drive losses) is not
/// declared here, so none is invented.
struct ProgressiveCavityPump
{
    double cavitiesPerRev = 0; // count of sealed cavities advanced per
                               // revolution (>= 1; the switch: > 0 declares
                               // the pump)
    double cavityVolume = 0;   // m3 transferred per cavity per revolution
    double rpm = 0;            // rotational speed of the drive
    double slip = 0;           // m3/s per Pa of leakage line
};

/// PUMPING is q > 0. NO-LOAD is dp <= 0: the machine sweeps but lifts nothing
/// (work and slip both read zero). STALLED is q <= 0 at a positive rise: the
/// declared slip line has eaten the whole displacement -- the pump turns and
/// the column falls back through the clearances as fast as the cavities
/// transfer it, so the net delivery is zero.
enum class ProgressiveCavityPumpState
{
    pumping = 0,
    noLoad = 1,
    stalled = 2
};

struct ProgressiveCavityPumpCard
{
    ProgressiveCavityPumpState state = ProgressiveCavityPumpState::pumping;

    double theoreticalRate = 0;   // q_th, m3/s
    double slipRate = 0;          // slip * dp, m3/s (0 in NO-LOAD)
    double netRate = 0;           // q, m3/s (0 in STALLED)
    double massRate = 0;          // rho * q, kg/s
    double pressureRise = 0;      // the (clamped) dp it is priced at, Pa
    double volumetricEfficiency = 0; // q / q_th
    double hydraulicPower = 0;    // q * dp, W
};

/// The card at an imposed pressure rise across the machine (the network's two
/// port pressures, or a declared differential). `fluidDensity` only converts
/// the rate to a mass rate; nothing else reads it.
ProgressiveCavityPumpCard progressiveCavityPumpCardOfLoad(
    const ProgressiveCavityPump& pump, double fluidDensity,
    double pressureRise);

/// The same card from two absolute port pressures: `pressureDischarge` is the
/// tubing side, `pressureIntake` the sand face side.
ProgressiveCavityPumpCard progressiveCavityPumpCardOf(
    const ProgressiveCavityPump& pump, double fluidDensity,
    double pressureDischarge, double pressureIntake);

} // namespace cf
