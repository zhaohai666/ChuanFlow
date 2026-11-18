#pragma once

#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// D8 increment 2: one-way Lagrangian particles in a solved 3D velocity field.
/// The carrier is frozen: particles sample U but never write momentum or
/// continuity back. This increment owns only the particle side.

/// One injection stream. Every particle in it shares the same diameter,
/// density and launch state; streams are declared as `stream0{} stream1{}`
/// blocks inside `dpm{}`.
struct DpmStream
{
    Vector point = Vector::zero();
    Vector velocity = Vector::zero();
    scalar diameter = 0;
    scalar density = 0;
    int count = 0;
};

/// The full `dpm{ ... }` declaration. `drag` currently accepts only the
/// published Schiller-Naumann sphere law; anything else is refused rather
/// than silently replaced.
struct DpmConfig
{
    std::vector<DpmStream> streams;
    scalar dt = 0;            // s
    int maxSteps = 0;
    scalar carrierRho = 0;    // kg/m3, needed by the particle Reynolds number
    Vector gravity = Vector::zero();
    std::string drag = "schillerNaumann";
    int reportEvery = 100;    // steps between trajectory samples
    int traceId = 0;          // particle whose trajectory is printed
    int traceSteps = 4;       // print EVERY step of this many first steps

    // Two-way option. Absent (the default) keeps the increment-2 one-way path
    // byte for byte: one frozen carrier solve and one trajectory pass.
    bool twoWay = false;
    int outerIterations = 5;  // carrier <-> trajectory Picard sweeps
    scalar twoWayRelax = 0.5; // relaxation of the cell force between sweeps
    scalar twoWayTol = 1e-3;  // relative change of the cell force to stop

    void validate() const;
};

DpmConfig readDpmConfig(const Dictionary& caseDict);

/// Schiller & Naumann (1935) drag coefficient of a rigid sphere:
///   Cd = 24/Re (1 + 0.15 Re^0.687),   Re <= 1000
///   Cd = 0.44,                        Re > 1000
/// The Stokes limit Cd = 24/Re is recovered algebraically as Re -> 0; the
/// implementation switches to its exact linear form below Re = 1e-6 so that
/// a zero relative velocity does not ask for 0/0.
scalar sphereDragCoefficient(scalar reynolds);

/// Acceleration due to drag alone, m/s2:
///   a = (3 rho_f Cd / (4 rho_p d)) (U_f - U_p)
/// which in the Stokes branch is exactly (U_f - U_p) / tau_p with
/// tau_p = rho_p d^2 / (18 mu). The buoyancy-corrected body term is added by
/// the integrator, not here, so a checker can grade the two separately.
Vector particleDragAcceleration(const Vector& uFluid, const Vector& uParticle,
                                scalar carrierRho, scalar rhoParticle,
                                scalar diameter, scalar mu);

enum class DpmStatus
{
    active,
    outlet,
    wall,
    expired,
};

const char* dpmStatusName(DpmStatus status);

/// Nearest cell-centre location with a mesh bounding-box guard. Sampling is
/// zero order: the owning cell's velocity is the particle's carrier velocity.
/// This stays defined all the way to the outlet face, where a strict inside
/// test would reject a particle the mesh still owns.
label locateDpmCell(const polyMesh& mesh, const Vector& point);

struct DpmParticle
{
    int id = 0;
    int stream = 0;
    DpmStatus status = DpmStatus::active;
    Vector x0 = Vector::zero();
    Vector u0 = Vector::zero();
    Vector x = Vector::zero();
    Vector u = Vector::zero();
    scalar age = 0;
    label cell = -1;
    scalar diameter = 0;
    scalar rhoParticle = 0;
};

/// One frozen-flow sample of the traced particle, before its velocity update.
struct DpmSample
{
    int id = 0;
    int step = 0;
    scalar t = 0;
    Vector x = Vector::zero();
    Vector u = Vector::zero();
    Vector uFluid = Vector::zero();
    scalar reynolds = 0;
    scalar cd = 0;
    Vector dragAcceleration = Vector::zero();
};

struct DpmResult
{
    std::vector<DpmParticle> particles;
    std::vector<DpmSample> samples;
    label outletCount = 0;
    label wallCount = 0;
    label expiredCount = 0;

    /// Drag reaction impulse (kg m/s) accumulated in every cell over the
    /// trajectories: sum_particles m_p integral(a_drag dt). The app converts
    /// it to a cell acceleration with
    ///   -impulse / (rho_f * V_cell * sourceWindow).
    /// Only the drag reaction is fed back: particle gravity is an external
    /// body force, not momentum taken from the carrier.
    std::vector<Vector> cellDragImpulse;
    /// Longest participating particle lifetime (s), used as the steady
    /// averaging window. Zero when no particle moved.
    scalar sourceWindow = 0;
};

/// `mu` is the carrier's dynamic viscosity (rho_f * nu); the app owns it
/// because the SIMPLE case declares nu, not mu.
DpmResult runDpm(const polyMesh& mesh, const volVectorField& U,
                 const DpmConfig& cfg, scalar mu);

} // namespace cf
