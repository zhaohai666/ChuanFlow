#pragma once

#include <string>
#include <vector>

#include "cf/pvt/Components.h"
#include "cf/pvt/PengRobinson.h"

namespace cf::pvt
{

/// Isothermal TP flash of a multi-component feed with the Peng-Robinson EOS:
/// Wilson K-value initialisation, Rachford-Rice for the vapour fraction at
/// each pressure-independent step, and successive substitution on the PR
/// fugacity equality. Vapour fraction beta is clamped to [0, 1], so a feed
/// that is single phase collapses to beta = 0 (liquid) or beta = 1 (vapour)
/// instead of diverging.

struct FlashControls
{
    scalar tol = 1.0e-12; // max |delta ln K| for convergence
    int maxIter = 1000;
};

struct FlashResult
{
    bool converged = false;
    int iterations = 0;
    scalar vaporFraction = 0; // beta, mole basis
    std::vector<scalar> liquidComposition; // x
    std::vector<scalar> vapourComposition; // y
    EosState liquid;
    EosState vapour;
    scalar liquidMassDensity = 0; // kg/m3
    scalar vapourMassDensity = 0; // kg/m3
    scalar liquidMolarMass = 0;   // kg/mol
    scalar vapourMolarMass = 0;   // kg/mol
    bool singlePhase = false;
    scalar residual = 0;          // last max |delta ln K|
    std::string message;
};

FlashResult tpFlash(const Mixture& feed, scalar T, scalar p,
                    const FlashControls& controls = FlashControls{});

/// Rachford-Rice residue f(beta) = sum z_i (K_i - 1) / (1 + beta (K_i - 1)).
scalar rachfordRiceResidue(const std::vector<scalar>& z,
                           const std::vector<scalar>& K, scalar beta);

/// Root of the Rachford-Rice function, assumed to exist in (0, 1).
scalar solveRachfordRice(const std::vector<scalar>& z, const std::vector<scalar>& K);

/// Density-driven phase label for a single-phase state: PR's own critical
/// molar volume is v_c = (Zc/Bc) b, so v/b below that ratio is liquid-like.
bool isLiquidLike(const EosState& state);

/// Bubble pressure of a feed at temperature T: the highest pressure at which
/// the flash still finds two phases. Pure fluids delegate to the PR fugacity
/// equality (the saturation pressure).
scalar bubblePressure(const Mixture& feed, scalar T);

} // namespace cf::pvt
