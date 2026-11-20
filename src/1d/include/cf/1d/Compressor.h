#pragma once

namespace cf
{

/// P10 (5.156): a steady point-element compressor for gas service, in the
/// two forms the field quotes a machine in.
///
/// Mode A, the CURVE: the catalog gives the pressure ratio at a rated rate
/// and (optionally) at shutoff, and the ratio falls with the square of the
/// rate -- the same quadratic family as the pump head curve of §5.81, read
/// in pressure ratio instead of metres of head:
///     PR(q) = shutoffRatio
///             - (shutoffRatio - pressureRatio) (|q| / ratedRate)^2
/// The curve is declared through its two points and its rate; the
/// thermodynamics (kappa, R_spec, Z, efficiency) are then asked what the
/// machine costs to run at the ratio the curve gives.
///
/// Mode B, the HEAD: the catalog gives the polytropic head H (J/kg) directly
/// and the ratio falls out of the declared suction state,
///     PR = (1 + H (k-1) / (k Z R_spec T1))^(k/(k-1)),
/// which is the same enthalpy integral inverted:
///     H = Z R_spec T1 (k/(k-1)) (PR^((k-1)/k) - 1).
/// The integral is priced with the DECLARED k (cp/cv) -- that is the
/// published convention this increment adopts, and the checker re-derives
/// the same closed form, not a real-gas enthalpy departure.
///
/// Both modes then report the same readouts at the operating point: the
/// discharge pressure p2 = p1 PR, the mass rate m = p1 q / (Z R_spec T1),
/// the ideal (gas) power W0 = m H, the shaft power W = W0 / eta, and the
/// discharge temperature T2 = T1 PR^((k-1)/(k eta)) -- the passive
/// temperature of §5.144's convention: the machine's thermodynamics are
/// priced FROM the hydraulic state, never fed back into it.
struct Compressor
{
    // Mode A: the ratio curve. pressureRatio > 0 switches the curve on.
    double pressureRatio = 0;      // PR at ratedRate (> 1)
    double shutoffRatio = 0;       // PR at zero rate (defaults to pressureRatio)
    double ratedRate = 0;          // m3/s at suction where PR is quoted
    // Mode B: the declared head. polytropicHead > 0 switches the head on.
    double polytropicHead = 0;     // J/kg
    // Shared thermodynamic statements.
    double polytropicEfficiency = 1; // eta, in (0,1]
    double kappa = 0;                // cp/cv, > 1
    double specificGasConstant = 0;  // R_spec, J/(kg K)
    double compressibility = 1;      // Z, a declared constant
};

struct CompressorSolution
{
    int mode = -1;                  // 0 curve, 1 head
    double pressureIn = 0;          // Pa abs
    double temperatureIn = 0;       // K
    double rate = 0;                // m3/s at suction
    double massRate = 0;            // kg/s
    double ratio = 1;               // PR at the operating rate
    double pressureOut = 0;         // Pa abs
    double head = 0;                // J/kg enthalpy integral at this PR
    double powerIdeal = 0;          // W to the gas, m H
    double powerShaft = 0;          // W / eta
    double temperatureOut = 0;      // K, passive discharge temperature
};

/// Price one operating point. Exactly one of the two modes must be declared
/// (pressureRatio > 0 XOR polytropicHead > 0); a machine that declares both
/// states two different things about the same ratio, and one that declares
/// neither states nothing. Out of envelope is refused loudly: non-positive
/// suction pressure/temperature/rate, eta outside (0,1], kappa <= 1 (the
/// exponent (k-1)/k is what makes a ratio a head), a non-positive R_spec or
/// Z, a curve whose shutoff ratio sits below its rated ratio (a machine that
/// boosts harder the faster it runs), a curve mode without a rated rate, and
/// a curve run past the rate where its ratio falls to unity.
CompressorSolution compressorSolve(const Compressor& c,
                                   double suctionPressurePa,
                                   double suctionTemperatureK,
                                   double rateSm3PerS);

} // namespace cf
