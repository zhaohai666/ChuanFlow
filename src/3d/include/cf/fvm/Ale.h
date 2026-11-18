#pragma once

#include "cf/common/Dictionary.h"
#include "cf/common/Labels.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// D8 increment 8: arbitrary Lagrangian--Eulerian mesh motion and the
/// Geometric Conservation Law. When the mesh moves with vertex velocity U_m,
/// a control volume's rate of change must equal the mesh-volume flux through
/// its faces (the Space Conservation Law, the geometric part of GCL):
///   d V_c / dt = sum_faces U_m,f . Sf_f
/// If this discrete identity fails, even a uniform field on a moving mesh
/// develops spurious sources: the ALE convective flux phi = (U - U_m).Sf
/// would not conserve a spatially constant scalar.
///
/// The mesh-velocity of a face is the average of its vertices' velocities (a
/// linear face motion is then exact), and Sf is the current (post-move)
///outward area vector. faceMeshVelocity + gclResidual implement that identity;
/// the checker verifies it for a rigid translation (where it is exact to
/// roundoff) and a uniform expansion (where the first-order-in-dt truncation
/// is measured and shown to shrink with dt).
struct AleConfig
{
    bool enabled = false;
    std::string motion = "translate"; // translate | expand
    Vector translation = Vector::zero(); // m/s, for motion=translate
    scalar expansionRate = 0;        // 1/s, for motion=expand (x -> (1+rate*t) x)
    scalar dt = 0;                   // s per motion step
    int steps = 0;
    int reportEvery = 1;
    void validate() const;
};

AleConfig readAleConfig(const Dictionary& caseDict);

/// Vertex velocity field from the motion declaration and the current points.
Field<Vector> meshVertexVelocity(const polyMesh& mesh, const AleConfig& cfg);

/// Move the mesh one step and return the per-cell GCL residual
///   dV/dt - sum_faces U_m,f . Sf_f
/// computed with the NEW geometry (call after mesh.movePoints). `oldVolumes`
/// are the cell volumes before the move.
Field<scalar> gclResidual(const polyMesh& mesh, const Field<scalar>& oldVolumes,
                          const Field<Vector>& vertexVelocity, scalar dt);

/// Advance vertex positions by dt under the declared motion (rigid translation
/// or uniform expansion about the origin).
Field<Vector> advanceVertices(const polyMesh& mesh, const AleConfig& cfg);

} // namespace cf
