#pragma once

#include "cf/common/Dictionary.h"
#include "cf/common/Field.h"
#include "cf/common/Vector.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// Moving-mesh relative-flux transport + a rigid-body "6DOF" update.
///
/// The mesh translates rigidly with a prescribed velocity U_mesh. On such a
/// mesh the conservative transport of a passive scalar is
///
///     d(phi V_c)/dt + sum_f phi_f (U - U_mesh).Sf = 0,
///
/// where (U - U_mesh).Sf is the relative face flux. Expanding with the
/// Geometric Conservation Law (the space conservation identity)
///
///     d V_c / dt = sum_f U_mesh,f . Sf
///
/// gives the equivalent non-conservative form. A spatially uniform scalar is
/// an exact solution iff the discrete GCL holds: replacing the geometric
/// volume sweep sum_f U_mesh.Sf by the EXACT cell-volume change removes the
/// O(dt) splitting error, so a uniform phi is preserved to roundoff for any
/// mesh motion. `transportStep` below applies that correction.
///
/// The rigid body is integrated alongside the mesh. We implement the subset
/// that has a closed-form reference: one translational DOF (along the force
/// direction, mass m + linear drag k) and one rotational DOF (about the
/// out-of-plane axis, inertia I + linear rotational drag c + stiffness K).
/// The translational ODE  m v' = F - k v  has terminal velocity F/k; the
/// rotational ODE  I theta'' + c theta' + K theta = 0  is a damped harmonic
/// oscillator. The remaining four DOFs are zero by construction and are
/// documented as outside this subset.
struct DynamicMeshConfig
{
    bool enabled = false;
    Vector meshVelocity = Vector::zero();  // prescribed constant mesh translation, m/s
    Vector bodyForce = Vector::zero();     // constant force on the body, N (x used)
    scalar mass = 0;                       // translational mass, kg (> 0)
    scalar drag = 0;                       // linear drag coefficient k, N/(m/s) (>= 0)
    scalar inertia = 0;                    // moment of inertia I, kg m^2 (> 0)
    scalar rotDrag = 0;                    // rotational drag c, N m s/rad (>= 0)
    scalar rotStiffness = 0;               // restoring stiffness K, N m/rad (>= 0)
    scalar theta0 = 0;                     // initial angle, rad
    scalar omega0 = 0;                     // initial angular velocity, rad/s
    scalar initialValue = 0;               // initial (uniform) scalar
    Vector fluidVelocity = Vector::zero(); // fluid velocity U, m/s
    scalar dt = 0;
    int steps = 0;
    int reportEvery = 1;
    void validate() const;
};

DynamicMeshConfig readDynamicMeshConfig(const Dictionary& caseDict);

/// State of the rigid body: translation (x, v) and rotation (theta, omega).
struct BodyState
{
    scalar x = 0;
    scalar v = 0;
    scalar theta = 0;
    scalar omega = 0;
};

/// Face mesh velocity (mean of vertex velocities; uniform here, but kept as a
/// field so the relative flux is assembled exactly like the ALE reference).
Field<Vector> meshFaceVelocity(const polyMesh& mesh, const Vector& Umesh);

/// Kinematic relative face flux phi_rel = (U - U_mesh).Sf, one entry per face.
Field<scalar> relativeFaceFlux(const polyMesh& mesh, const Vector& U, const Vector& Umesh);

/// Per-cell geometric mesh-volume sweep sum_faces U_mesh.Sf (owner +, neighbour -).
Field<scalar> gclGeometricSweep(const polyMesh& mesh, const Vector& Umesh);

/// Per-cell GCL residual dV/dt - sum_faces U_mesh.Sf using the NEW geometry
/// (call after mesh.movePoints); `oldVolumes` are the pre-move cell volumes.
Field<scalar> gclCellResidual(const polyMesh& mesh, const Vector& Umesh,
                              const Field<scalar>& oldVolumes, scalar dt);

/// One GCL-corrected explicit transport step. `phiRel` is the relative face
/// flux, `sweep` is sum_faces U_mesh.Sf per cell, `oldVolumes` the pre-move
/// volumes. Returns the new cell-centred scalar (current mesh geometry).
Field<scalar> transportStep(const polyMesh& mesh, const Field<scalar>& phi,
                            const Field<scalar>& phiRel, const Field<scalar>& sweep,
                            const Field<scalar>& oldVolumes, scalar dt);

/// RK4 integration of the translational + rotational rigid-body ODE over dt.
BodyState advanceBody(const BodyState& s, const DynamicMeshConfig& cfg, scalar dt);

} // namespace cf
