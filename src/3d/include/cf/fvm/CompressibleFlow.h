#pragma once

#include <functional>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Vector.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/PolyMesh.h"

namespace cf::compressible
{

// Conserved state [rho, rho*U, rho*E]. Mirrors the closed-box Euler app's
// State/Flux exactly, so the library HLLC is bit-identical to the original.
struct State
{
    scalar rho = 0;
    Vector mom = Vector::zero();
    scalar rhoE = 0;
};

struct Flux
{
    scalar mass = 0;
    Vector mom = Vector::zero();
    scalar energy = 0;
};

/// Ideal-gas EOS: p = (gamma-1) * rho * e, a = sqrt(gamma*p/rho).
/// `temperature` needs a specific gas constant R (J/kg/K) declared by the
/// case; it is opt-in and never used by the inviscid path.
struct IdealGasEos
{
    scalar gamma = 0;
    scalar R = 0; // 0 => temperature not available; only p/e/a are used

    static IdealGasEos fromCase(const Dictionary& caseDict);

    scalar pressure(scalar rho, scalar rhoE, const Vector& mom) const;
    scalar soundSpeed(scalar rho, scalar p) const;
    scalar temperature(scalar p, scalar rho) const;
};

Flux physicalFlux(const State& u, scalar p, scalar un, const Vector& nHat);
scalar twoRarefactionPressure(scalar pL, scalar unL, scalar cL,
                              scalar pR, scalar unR, scalar cR, scalar gamma);
Flux hllc(const State& uL, scalar pL, const Vector& vL, scalar cL,
          const State& uR, scalar pR, const Vector& vR, scalar cR,
          const Vector& nHat, scalar gamma);

/// Opt-in transport coefficients. Missing blocks leave viscosity/conductivity
/// at zero, which makes every function below a no-op.
struct TransportConfig
{
    scalar mu = 0;    // dynamic viscosity [Pa s], 0 => no viscous stress
    scalar kappa = 0; // thermal conductivity [W/m/K], 0 => no heat flux

    static TransportConfig fromCase(const Dictionary& caseDict);

    bool viscous() const { return mu > 0.0; }
    bool conducting() const { return kappa > 0.0; }
};

/// Newtonian viscous stress with Stokes' hypothesis:
///   tau_ij = mu (dU_i/dx_j + dU_j/dx_i) - (2/3) mu div(U) delta_ij
/// The face flux returned is tau . Sf (already oriented with the face normal).
/// gradU[comp][cell][dir] uses the same per-component Gauss gradients that
/// GeneralizedNewtonianViscosity::computeStrainRate already validates.
Vector viscousStressFaceFlux(const polyMesh& mesh,
                             const Field<Vector>* gradU,
                             scalar mu, label face);

/// Per-cell divergence of the viscous stress, sum_f tau.Sf / V, accumulated
/// with the same owner/neighbour scatter as the Euler faces.
Field<Vector> viscousStressDivergence(const polyMesh& mesh,
                                      const Field<Vector>* gradU,
                                      scalar mu);

/// Fourier heat flux through an internal face, q . Sf = -kappa (gradT . Sf).
/// gradT is the cell-centred Gauss gradient; the face value is the bounded
/// owner-weight interpolation used everywhere else in the FVM layer.
scalar heatFluxFace(const polyMesh& mesh, const Field<Vector>& gradT,
                    scalar kappa, label face);

/// Per-cell divergence of the heat flux (sum_f q.Sf / V).
Field<scalar> heatFluxDivergence(const polyMesh& mesh,
                                 const Field<Vector>& gradT,
                                 scalar kappa);

/// Per-component gradient of a cell field using the FVM Gauss scheme.
Field<Vector> cellGradient(const polyMesh& mesh, const Field<scalar>& data,
                           const std::string& gradScheme);

/// Per-component velocity gradient gradU[comp][cell][dir].
void velocityGradients(const polyMesh& mesh, const Field<Vector>& U,
                       const std::string& gradScheme, Field<Vector> gradU[3]);

/// Validate transport coefficients before use (mu/kappa must be finite and
/// non-negative; zero is allowed and means disabled). Throws cfFatal.
void validateTransport(const TransportConfig& cfg);

} // namespace cf::compressible
