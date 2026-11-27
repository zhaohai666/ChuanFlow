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

// =====================================================================
// 5.227: the DO sweep extended from the 1D slab to a true 3D box.
//
// The slab path above stays exactly as shipped (F9); this adds a level
// symmetric ordinate set in 3D and a per-direction finite-volume sweep on
// a blockMesh box:
//
//   div(s . I) = -(kappa_a + kappa_s) I + kappa_a I_b + kappa_s/(4 pi) sum w I
//
// discretised first-order upwind per axis (the inflow face value is the
// upwind cell centre), grey diffuse walls with a declared emissivity
// (I_in = eps I_b(T_w) + (1 - eps) I_out of the mirrored direction, lagged
// one source iteration), and symmetry/empty faces that return the outgoing
// intensity back unchanged, so that axis cancels exactly.
//
// Direction sets (weights normalised to sum 4 pi, set closed under the full
// octahedral group, so every second moment is isotropic):
//   order 2 (S2): 8 ordinates  (+-1,+-1,+-1)/sqrt(3), w = pi/2
//   order 4 (S4): 24 ordinates  orbit of (+-2,+-1,+-1)/sqrt(6), w = pi/6
struct Do3dWall
{
    std::string type = "wall";   // wall | symmetry | empty
    scalar temperature = 0;      // K (walls)
    scalar emissivity = 1;       // grey wall, [0,1]
};

struct Do3dConfig
{
    int order = 2;               // 2 | 4
    scalar absorption = 0;       // kappa_a [1/m], >= 0
    scalar scattering = 0;       // kappa_s [1/m], >= 0
    scalar mediumT1 = 0;         // K at x = xmin (uniform when T2 == T1)
    scalar mediumT2 = 0;         // K at x = xmax
    scalar sbConstant = 5.670374419e-8;
    Do3dWall walls[6];           // x0 x1 y0 y1 z0 z1
    scalar sourceTol = 1e-13;    // source-iteration stop
    int maxIterations = 5000;

    scalar extinction() const { return absorption + scattering; }
    void validate() const;
};

/// Read `geometry box` and the `radiation3d{}` block (order, kappa_a,
/// kappa_s, mediumT1/2, sigmaSB, sourceTol, maxIterations). The six walls are
/// filled by the caller from the block face definitions and the `patches{}`
/// entries. Rejects order not 2/4, emissivity outside [0,1], negative kappa.
Do3dConfig readDo3dConfig(const Dictionary& caseDict);

/// One 3D level-symmetric ordinate: direction cosines + weight.
struct Ordinate3d
{
    Vector s = Vector::zero();
    scalar w = 0;
};

/// The declared S2 (8) / S4 (24) set; any other order fails loudly.
std::vector<Ordinate3d> levelSymOrdinates(int order);

/// Result of a steady 3D box DO solve on an nx x ny x nz blockMesh.
struct RadBoxResult
{
    int nx = 0, ny = 0, nz = 0;
    Field<scalar> G;             // cell-averaged incident radiation [W/m2]
    Field<Vector> q;             // cell net radiative flux [W/m2]
    Field<scalar> S;             // cell net volumetric source [W/m3]
    scalar qX0 = 0;              // net flux INTO the medium at x = xmin [W/m2]
    scalar qX1 = 0;              // net flux INTO the medium at x = xmax
    scalar qY0 = 0, qY1 = 0, qZ0 = 0, qZ1 = 0;
    scalar volSrc = 0;           // integral S dV / V
    scalar balance = 0;          // sum of the six wall fluxes - integral S dV
    int directions = 0;
    int iterations = 0;
    scalar residual = 0;
    bool converged = true;
};

/// 3D discrete-ordinates FVM sweep on a uniform box. The medium temperature
/// varies linearly along x between T1 and T2 (uniform when equal).
RadBoxResult solveDo3d(int nx, int ny, int nz,
                       scalar Lx, scalar Ly, scalar Lz,
                       const Do3dConfig& cfg);

} // namespace cf
