#pragma once

#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Labels.h"

namespace cf
{

/// Batch D / F10: a small but genuine two-step chain-branching mechanism in a
/// homogeneous 0D constant-volume reactor. The existing single-step
/// `Combustion.h` (steady multispecies transport) is deliberately untouched.
///
///     step 1 (branching)    F + R -> P + 2R   r1 = A1 exp(-E1/(R T)) Yf Yr
///     step 2 (termination)  R     -> P         r2 = A2 exp(-E2/(R T)) Yr
///
/// with mass-fraction state (Yf, Yr) and temperature T obeying
///
///     dYf/dt = -r1
///     dYr/dt = +r1 - r2
///     dT/dt  = q r1 - cooling (T - Tenv)
///
/// `q` is the adiabatic temperature rise per unit fuel consumed, so for
/// cooling = 0 the quantity T + q Yf is an exact invariant (T + q Yf =
/// T0 + q Yf0). `cooling` is a Newton heat-loss coefficient [1/s]; its
/// critical value, found by bisection on a full reactor simulation, is the
/// extinction Damkoehler threshold reported by the app and re-derived by the
/// independent checker. The ignition delay is the first time T crosses
/// T0 + ignitionRise * q * Yf0, also a simulated outcome rather than a string.
///
/// All rates are Arrhenius with case-declared pre-exponentials, activation
/// energies and gas constant; nothing is transcribed from a correlation.
struct DetailedChemistryConfig
{
    bool enabled = false;
    scalar A1 = 0, E1 = 0;       // branching: pre-exp [1/s], Ea [J/mol]
    scalar A2 = 0, E2 = 0;       // termination: pre-exp [1/s], Ea [J/mol]
    scalar q = 0;                // adiabatic heat release [K]
    scalar T0 = 0;               // initial temperature [K]
    scalar gasConstant = 8.314462618; // R [J/mol/K], CODATA 2019 exact
    scalar Yf0 = 1, Yr0 = 1e-6;  // initial fuel / radical mass fractions
    scalar cooling = 0;          // heat-loss coefficient [1/s] (0 = adiabatic)
    scalar Tenv = 0;             // ambient temperature [K] (default T0)
    scalar dt = 0.01;            // output interval [s]
    int nSteps = 600;            // number of output steps
    int nSub = 10;               // RK4 substeps per output interval
    scalar ignitionRise = 0.3;   // temperature-rise fraction defining ignition

    void validate() const;
};

struct DetailedState
{
    scalar t = 0, T = 0, Yf = 0, Yr = 0;
};

struct DetailedDiagnostics
{
    scalar branchDa = 0;    // k1(T0) Yf0 / k2(T0): branching vs termination
    scalar betaCrit = 0;    // critical cooling from bisection (0 => never ignites)
    scalar DaCrit = 0;      // k1(T0) Yf0 / (k2(T0) + betaCrit)
    scalar ignitionDelay = -1; // first ignition crossing [s], -1 if none
    scalar peakT = 0;       // peak temperature of the simulated case [K]
    bool ignited = false;
};

/// k(T) = A exp(-Ea/(R T)). Exposed so a checker can grade the temperature
/// law independently of the reactor integration.
scalar detailedArrhenius(scalar A, scalar Ea, scalar R, scalar T);

DetailedChemistryConfig readDetailedChemistryConfig(const Dictionary& caseDict);

/// Integrate the reactor with RK4 substeps, returning a state at t = 0 and
/// after every output interval (hist.size() == nSteps + 1).
std::vector<DetailedState> integrateDetailed(const DetailedChemistryConfig& cfg);

/// Simulated outcomes: ignition delay, peak T and the extinction Damkoehler
/// threshold (critical cooling) found by bisection on the same ODEs.
DetailedDiagnostics diagnoseDetailed(const DetailedChemistryConfig& cfg);

} // namespace cf
