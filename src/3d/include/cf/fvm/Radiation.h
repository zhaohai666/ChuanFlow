#pragma once

#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// D8 increment 4: the P1 (diffusion) approximation of grey radiative transfer
/// through an absorbing / isotropically scattering participating medium.
///
/// The equation and the constants are the published P1 model (Modest,
/// "Radiative Heat Transfer", 3rd ed., ch. 10; Siegel & Howell, "Thermal
/// Radiation Heat Transfer"). Nothing here is an empirical correlation:
///
///   div( 1/(3 beta) grad G ) - kappa_a G + 4 kappa_a sigma T^4 = 0
///   beta  = kappa_a + kappa_s                 [1/m, extinction]
///   q_r   = -1/(3 beta) grad G                [W/m2, radiative flux]
///   S_rad = kappa_a (G - 4 sigma T^4)         [W/m3, net radiative heating]
///
/// G is the angle-integrated incident radiation [W/m2]. The source term is the
/// divergence of the radiative flux with the opposite sign:
///   -div(q_r) = kappa_a (G - 4 sigma T^4) = S_rad,
/// so integrating over the domain gives the exact balance the checker grades:
///   integral(S_rad dV) + surface_integral(q_r . n dA) = 0.
///
/// At a black wall in radiative equilibrium the incident radiation is
///   G_wall = 4 pi I_b = 4 sigma T_wall^4
/// (I_b = sigma T^4 / pi). The app reads `wallTemperature` on each patch and
/// applies that fixed value through the usual scalar-field BC; patches without
/// a wall temperature (symmetry / empty / adiabatic) get zeroGradient.
///
/// sigma is the CODATA 2019 exact Stefan-Boltzmann constant; it is a physical
/// constant, not a case input, so a case cannot silently change the scale.
struct RadiationConfig
{
    bool enabled = false;
    scalar absorption = 0;         // kappa_a [1/m], must be > 0
    scalar scattering = 0;         // kappa_s [1/m], must be >= 0
    scalar mediumTemperature = 0;  // K, uniform participating-medium T (> 0)
    scalar sbConstant = 5.670374419e-8; // W/m2/K4, CODATA 2019 (exact)

    static constexpr scalar sigmaSB = 5.670374419e-8;

    /// beta = kappa_a + kappa_s (extinction coefficient).
    scalar extinction() const { return absorption + scattering; }
    /// Photon mean free path 1/(3 beta): the P1 diffusion coefficient.
    scalar diffusion() const { return 1.0 / (3.0 * extinction()); }

    void validate() const;
};

RadiationConfig readRadiationConfig(const Dictionary& caseDict);

/// Result of one steady P1 solve, all per cell unless stated otherwise.
struct RadiationResult
{
    Field<scalar> G;          // incident radiation [W/m2]
    Field<Vector> qRad;       // radiative flux [W/m2]
    Field<scalar> netSource;  // kappa_a (G - 4 sigma T^4) [W/m3]
    scalar boundaryFlux = 0;  // surface_integral(q_r . n dA) [W]
    scalar volumeSource = 0;  // integral(S_rad dV) [W]
    scalar balance = 0;       // boundaryFlux + volumeSource, ~0
    int iterations = 0;
    scalar residual = 0;
    bool solved = false;
};

/// Solve the steady P1 equation on `mesh`. `G` carries the boundary conditions
/// already set by the app (fixedValue = 4 sigma T_wall^4 on radiating walls,
/// zeroGradient elsewhere). `mediumT` is the participating-medium temperature
/// per cell (uniform in the shipped case). Returns the converged fields and the
/// discrete energy balance; does not touch velocity, pressure or any energy
/// field -- radiative feedback into temperature is a later increment.
RadiationResult solveP1(const polyMesh& mesh, volScalarField& G,
                        const RadiationConfig& cfg, const Schemes& schemes,
                        const Field<scalar>& mediumT);

/// Incident radiation at a black wall: G_w = 4 pi I_b = 4 sigma T_w^4.
scalar blackbodyIncidentRadiation(scalar temperature, scalar sigmaSB);

} // namespace cf
