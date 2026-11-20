#pragma once

namespace cf
{

/// A continuous-flow liquid-liquid jet pump: a fixed nozzle throat geometry
/// that hands part of the motive stream's momentum to a suction stream in a
/// constant-area mixing tube, and then diffuses the combined stream to the
/// discharge pressure.
///
/// The model family is the classic constant-area mixing momentum theory of
/// jet pumps -- Thornhcser, "Theory and performance of jet pumps" (1955), as
/// carried in Karassik et al., "Pump Handbook" (3rd ed., ch. 41) and in the
/// jet-pump chapter of Guo, Lyons & Ghalambor, "Petroleum Production
/// Engineering" (2007). Three statements, all of them algebra rather than a
/// curve fit:
///
///  1. Continuity. The discharge is the sum of the two feeds,
///         q_d = q_n + q_s,
///     and the flow ratio is r = q_s / q_n.
///
///  2. The nozzle converts the motive stream's pressure into velocity with no
///     declared loss: with the motive RATE declared (the machine's pump, not
///     this element, sets it), v_n = q_n / A_n is known, and the supply
///     pressure the nozzle needs is p_supply = p_suction + rho v_n^2 / 2 --
///     Bernoulli read backwards.
///
///  3. Momentum on the constant-area throat (area A_t > A_n, the "throat
///     constant" m = A_t/A_n). The control volume runs from the suction
///     chamber plane -- jet of area A_n at v_n beside suction of area
///     A_t - A_n at v_s = q_s/(A_t - A_n), both under the chamber pressure
///     p_suction -- to the throat outlet at v_3 = q_d/A_t and pressure
///     p_discharge, and it pays for the tube's wall friction as a Darcy line,
///     f (L_t/D_t) rho v_3^2/2. The balance is
///         p_d - p_s = rho v_n^2 (m-1-r)^2 / (m^2 (m-1))
///                   - f (L_t/D_t) rho v_n^2 (1+r)^2 / (2 m^2).
///     The perfect square is not quoted from a table, it is the algebra:
///     clearing m^2(m-1) into 1/m + r^2/(m(m-1)) - (1+r)^2/m^2 gives
///     m^2-m + r^2 m - (m-1)(1+r)^2 = (m-1-r)^2 term for term. At f = 0 it
///     is exactly Thornhcser's ideal result, and it states the two limits the
///     machine lives between: r = 0 is the zero-suction stall (the pump is a
///     plain nozzle and still makes pressure), r = m-1 is the free-jet point
///     (zero pressure rise, maximum entrainment).
///
/// Given the motive rate, the geometry, the fluid and the two port pressures,
/// the suction rate is the root of that quadratic in r lying in [0, m-1] --
/// unique, because the left-hand side's matching function g(r) is strictly
/// decreasing there. When even r = 0 cannot hold the discharge pressure up the
/// pump is NO-SUCTION and r clamps to 0; the declared-rate mode inverts the
/// same equation for the suction pressure instead (it is linear in dp).
///
/// The diffuser downstream of the throat is assumed ideal at this increment:
/// p_discharge is the pressure the element is quoted at, and no separate
/// diffuser efficiency is declared. Nozzle losses, two-phase motive fluid and
/// gas in the suction stream are outside the declaration surface.
struct JetPump
{
    double nozzleArea = 0;      // m2, A_n (the switch: > 0 declares the pump)
    double throatArea = 0;      // m2, A_t; must exceed the nozzle area
    double throatLength = 0;    // m, L_t of the constant-area mixing tube
    double frictionFactor = 0;  // Darcy f for the throat wall (0 = ideal)
    double fluidDensity = 0;    // kg/m3, incompressible single-phase
};

/// How the element is told what to solve for. Exactly one of the two modes is
/// used per call: either the suction RATE is unknown and the suction pressure
/// is given (the network/field case), or the suction rate is imposed and the
/// suction pressure falls out of the same momentum balance (the catalog case).
struct JetPumpOperation
{
    double motiveRate = 0;         // m3/s through the nozzle (q_n > 0)
    double suctionPressure = 0;    // Pa abs, chamber pressure p_s (mode A)
    double dischargePressure = 0;  // Pa abs, diffuser outlet p_d
    bool imposeSuctionRate = false;
    double suctionRate = 0;        // m3/s, only read when imposed (mode B)
};

/// The three states: PUMPING is the root in (0, m-1]; NO-SUCTION is the
/// back-pressure stall (even r = 0 overshoots the momentum balance, so the
/// chamber simply takes nothing in); FREE-JET is the r = m-1 limit at zero
/// pressure rise, reported so a caller cannot mistake it for a fault.
enum class JetPumpState
{
    pumping = 0,
    noSuction = 1,
    freeJet = 2
};

/// Everything one operating point states. `flowRatio` is r, `areaRatio` m;
/// `supplyPressure` is the motive pressure the nozzle needs at the solved
/// chamber pressure (mode A) or the one implied by the declared supply
/// (mode B: it is the input p_supply = p_s + rho v_n^2/2, echoed back).
/// `headRatio` N = (p_d - p_s)/(p_supply - p_d) and `efficiency` = r N are
/// the two published jet-pump numbers, recomputed here from the element's own
/// pressures rather than looked up: the hydraulic power the suction stream
/// gains is q_s (p_d - p_s), the motive stream spends q_n (p_supply - p_d).
struct JetPumpSolution
{
    JetPumpState state = JetPumpState::pumping;

    double areaRatio = 0;        // m = A_t/A_n
    double flowRatio = 0;        // r = q_s/q_n
    double nozzleVelocity = 0;   // v_n, m/s
    double throatVelocity = 0;   // v_3, m/s

    double motiveRate = 0;       // q_n, m3/s
    double suctionRate = 0;      // q_s, m3/s (0 in NO-SUCTION)
    double dischargeRate = 0;    // q_d = q_n + q_s, m3/s

    double suctionPressure = 0;  // Pa abs (input in mode A, answer in mode B)
    double dischargePressure = 0;// Pa abs
    double supplyPressure = 0;   // Pa abs the nozzle needs upstream

    double headRatio = 0;        // N = (p_d - p_s)/(p_supply - p_d)
    double efficiency = 0;       // E = r N = q_s(p_d-p_s) / q_n(p_supply-p_d)
    double usefulPower = 0;      // q_s (p_d - p_s), W
    double motivePower = 0;      // q_n (p_supply - p_d), W
    double momentumResidual = 0; // Pa, the printed balance's LHS - RHS (0 when
                                 // the root was found; clamped states carry the
                                 // overshoot so the stall is a number, not a
                                 // silence)
};

/// The solve. Throws (cfFatal) on a non-positive area, a throat no larger than
/// the nozzle (there is no annulus to entrain through), a non-positive
/// density, a non-positive motive rate, negative length or friction factor, a
/// non-positive absolute pressure, or an imposed suction rate outside
/// [0, (m-1) q_n] -- the published operating envelope of the geometry.
JetPumpSolution jetPumpSolve(const JetPump& pump, const JetPumpOperation& op);

} // namespace cf
