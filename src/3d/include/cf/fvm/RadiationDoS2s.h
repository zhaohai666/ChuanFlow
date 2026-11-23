#pragma once

#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Field.h"
#include "cf/common/Vector.h"

namespace cf
{

// =====================================================================
// F9: grey radiation beyond P1 -- discrete ordinates (S_N) and the
// surface-to-surface (S2S) enclosure / transparent-medium limit.
//
// Everything in this header is independent of the shipped P1 solver
// (Radiation.h / radiation.cpp): P1 remains untouched. The one binary
// chflow-3d-radiation-do selects the physics with `model do|s2s|p1`;
// the p1 branch deliberately delegates to the existing solveP1 so the
// checker can grade a DO/P1 agreement trend on identical case data.
//
// Physics reference: Modest, "Radiative Heat Transfer", 3rd ed., ch. 16-17;
// Siegel & Howell, "Thermal Radiation Heat Transfer". No empirical
// correlations: every number below comes from the RTE or the view-factor
// algebra. All coefficients are case inputs; sigma is the CODATA 2019
// exact Stefan-Boltzmann constant.
// =====================================================================

/// Radiation physics selected by the `model` switch.
enum class RadiationModel { DO, S2S, P1 };

/// One discrete ordinate on [-1,1]; w are Gauss-Legendre weights with
/// sum(w) = 2. A symmetric S_N set has N/2 positive and N/2 negative mu.
struct Ordinate
{
    scalar mu = 0;
    scalar w = 0;
};

/// Symmetric S4 / S8 Gauss-Legendre quadrature (order = number of ordinates).
std::vector<Ordinate> snQuadrature(int order);

struct DoS2sConfig
{
    RadiationModel model = RadiationModel::DO;
    scalar absorption = 0;         // kappa_a [1/m], must be > 0
    scalar scattering = 0;         // kappa_s [1/m], must be >= 0
    scalar mediumTemperature = 0;  // K, uniform participating medium
    scalar sbConstant = 5.670374419e-8; // CODATA 2019 exact
    int quadratureOrder = 8;       // 4 or 8

    scalar extinction() const { return absorption + scattering; }
    void validate() const;
};

/// Read `model`, the `radiation{}` block and the optional quadrature order.
/// Unknown model names / orders fail loudly through cfFatal.
DoS2sConfig readDoS2sConfig(const Dictionary& caseDict);

/// Result of a steady 1D slab solve. G and S are cell-AVERAGED (the averaging
/// is what makes the discrete wall-source balance close), q is the flux
/// reconstructed at the cell CENTRE from the exponential ray solution.
/// qLeft / qRight are net fluxes INTO the medium from the left / right wall
/// (W/m2); volSrc = integral(S dx) per unit slab area (W/m2); balance is
/// qLeft + qRight - volSrc and is ~0 by the discrete RTE identity.
struct RadSlabResult
{
    Field<scalar> G;          // cell-averaged incident radiation [W/m2]
    Field<Vector> q;          // cell-centre net radiative flux [W/m2]
    Field<scalar> S;          // cell-averaged net volumetric source [W/m3]
    scalar qLeft = 0;
    scalar qRight = 0;
    scalar volSrc = 0;
    scalar balance = 0;
    int iterations = 0;       // source iterations (1 = direct closed form)
    scalar residual = 0;
    bool converged = true;
    scalar viewFactor = 1;    // S2S F_LR (=1 for infinite parallel plates)
};

/// Discrete-ordinates sweep on a uniform 1D grey slab between two black
/// walls at TL (x=0) and TR (x=L). Absorption and isotropic scattering are
/// uniform. Each ordinate is integrated with exponential fitting (exact for
/// a piecewise-constant source), so a non-scattering slab is reproduced to
/// machine precision; scattering is converged by source iteration.
RadSlabResult solveDoSlab(int nCells, scalar length, scalar TL, scalar TR,
                          const DoS2sConfig& cfg);

/// Surface-to-surface grey enclosure in the documented transparent-medium
/// limit: two directly opposed, equal-area black parallel plates form a
/// closed two-surface enclosure, so F_LR = F_RL = 1, F_LL = F_RR = 0
/// (summation + reciprocity). The net exchange is sigma (T_L^4 - T_R^4),
/// uniform across the gap, and the medium contributes no volumetric source.
/// `viewFactor` reports F_LR so a checker can reconstruct the closure.
RadSlabResult solveS2sSlab(int nCells, scalar length, scalar TL, scalar TR,
                           const DoS2sConfig& cfg);

/// Blackbody incident radiation G_w = 4 pi I_b = 4 sigma T^4.
scalar dos2sBlackbodyG(scalar temperature, scalar sigmaSB);

} // namespace cf
