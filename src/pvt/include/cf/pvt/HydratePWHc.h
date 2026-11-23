#pragma once

#include <string>
#include <vector>

#include "cf/common/Labels.h"
#include "cf/pvt/Components.h" // Mixture, kGasConstant

namespace cf::pvt
{

/// Simplified van der Waals - Platteeuw (1959) hydrate equilibrium in the
/// Parrish - Prausnitz style, written entirely from declared inputs so an
/// independent checker can reproduce every reported number with plain math
/// (nothing is tabulated inside the kernel):
///
///   occupancy      theta_j = sum_i C_ij f_i / (1 + sum_i C_ij f_i)
///                  per cage type j (0 = small, 1 = large); C_ij are the
///                  declared Langmuir constants (1/Pa) and f_i the component
///                  fugacity of the gas exposed to the hydrate, taken from the
///                  Peng-Robinson flash (vapour phase when two-phase, the
///                  single-phase state otherwise).
///
///   driving force  Delta_mu_w(T) = muRef + A_w (T - Tref)
///                  chemical-potential difference between the empty hydrate
///                  lattice (beta) and the bulk water phase, linearised about
///                  the declared anchor with two declared constants.
///
///   equilibrium    Delta_mu_w(T) + R T sum_j nu_j ln(1 - theta_j) = 0
///                  with nu_j the cages per water molecule of the unit cell:
///                  sI  46 H2O/cell, 2 small + 6 large cages
///                      -> nu_small = 2/46, nu_large = 6/46
///                  sII 136 H2O/cell, 16 small + 8 large cages
///                      -> nu_small = 16/136, nu_large = 8/136.
///
/// The reported equilibrium pressure at a declared temperature is the root of
/// that equation in p (theta rises with p, so the residual falls monotonically
/// through zero). referencePressure only seeds the search bracket and
/// enthalpyOfDissociation is a declared Clausius-Clapeyron physical bound on
/// the resulting d ln p/dT slope; neither enters the solved equation.
struct HydratePWHcSpec
{
    std::string structure = "sII";  // "sI" or "sII"
    std::vector<std::string> guests; // feed component names exposed to cages

    /// Langmuir constants [guest][cage] in 1/Pa, cage 0 = small, 1 = large.
    std::vector<std::vector<scalar>> langmuir;

    scalar referenceTemperature = 273.15; // K, anchor of the linearisation
    scalar referencePressure = 2.6e6;     // Pa, search seed only (anchor)
    scalar enthalpyOfDissociation = 6.5e4; // J/mol, CC slope physical bound
    scalar muRef = 0;                     // J/mol, Delta_mu_w(Tref)
    scalar muSlopeA = 0;                  // J/(mol K), A_w slope
};

struct HydratePWHcPoint
{
    scalar pEq = 0;         // Pa, solved equilibrium pressure at T
    scalar thetaSmall = 0;  // small-cage occupancy at pEq
    scalar thetaLarge = 0;  // large-cage occupancy at pEq
};

/// Cages per water molecule (nu_j) for a declared structure: cage 0 = small,
/// cage 1 = large. Fatal on an unknown structure name.
scalar hydratePWHcCagePerWater(const std::string& structure, int cage);

/// Equilibrium pressure of the declared hydrate at temperature T for the gas
/// `feed` (its flashed vapour fugacity supplies f_i). Bisection between a
/// vanishing-pressure bracket (theta -> 0, residual -> Delta_mu_w > 0) and a
/// 1e10 Pa ceiling; a root that never enters the ceiling is a loud fatal.
HydratePWHcPoint hydratePWHcEquilibrium(const Mixture& feed, scalar T,
                                        const HydratePWHcSpec& spec);

} // namespace cf::pvt
