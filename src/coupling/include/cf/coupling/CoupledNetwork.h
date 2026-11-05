#pragma once

#include <string>
#include <vector>

#include "cf/1d/Network.h"
#include "cf/common/Dictionary.h"
#include "cf/common/Labels.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/PolyMesh.h"
#include "cf/fvm/SimpleSolver.h"

namespace cf
{
namespace coupling
{

/// The partner model behind a coupled 1D element (`NetworkPipe::coupledResistance`):
/// it states the pressure drop of a region at a given volume rate, and knows
/// nothing about the network. Purely virtual so the interface itself can be
/// regressed against an analytic law; the production implementation is a 3D
/// solve (`DetailBlock`).
class ICoupledElement
{
public:
    virtual ~ICoupledElement() = default;

    /// Pressure drop in Pa across the region when it carries `qMag` m3/s.
    /// Always positive: the element resists the flow, whichever way it runs.
    virtual scalar lossAt(scalar qMag) = 0;

    /// Diagnostics of the last `lossAt`, recorded into the interface history so
    /// a regression can tell "the partner needed more sweeps" from "the
    /// interface stalled". An analytic law is exact at the rate it was given,
    /// so the defaults say one cost-free, converged evaluation.
    virtual label lastSweeps() const { return 0; }
    virtual bool lastConverged() const { return true; }
};

/// What the last solve of a `DetailBlock` produced. Every number is measured on
/// the 3D field, so a regression can tell "the block reported a different drop"
/// from "the block stopped conserving".
struct BlockMeasure
{
    label cells = 0;
    label sweeps = 0;         // SIMPLE outer iterations
    bool converged = false;

    scalar dpPa = 0;          // rho * (volume-weighted p of the first cross-section
                              // minus that of the last) -- the drop the network is charged
    scalar imposedQ = 0;      // m3/s the block was asked to carry
    scalar fluxIn = 0;        // signed sum of phi over the inlet patch
    scalar fluxOut = 0;       // ... over the outlet patch
    scalar planeFluxMin = 0;  // min/max over every cross-section plane, including the
    scalar planeFluxMax = 0;  // patch ones: the discrete continuity statement

    scalar probeSlope = 0;    // Pa/m, least-squares slope of rho*p over the probe
                              // cross-sections: the developed-region gradient
    scalar probeSpread = 0;   // Pa, max distance of a probe value to that line
    label probeSections = 0;

    scalar inletArea = 0;     // m2, |sum Sf| over the inlet patch
    scalar plugSpeed = 0;     // m/s, the imposed uniform inlet velocity
    scalar uxMax = 0;         // m/s, largest streamwise velocity in the block
    scalar uxCentre = 0;      // m/s, ... restricted to the probe cross-sections
};

/// One steady incompressible duct, solved in 3D and reduced to the single number
/// the 1D network needs: its pressure drop at an imposed volume rate.
///
/// The boundary conditions come from the mesh's own patch types and nothing else
/// -- `inlet` carries the plug velocity implied by the rate, `outlet` carries the
/// pressure reference, `wall` sticks, everything else is a zero-gradient/symmetry
/// side -- because the coupling must not need a second, easily-inconsistent
/// description of the same interface. The pressure is reported in Pa: the 3D
/// kernel's p is kinematic (p/rho, the momentum equation is forced by
/// `bodyForce` in m/s2), while the 1D network works in Pa.
///
/// Each `lossAt` call re-solves from the previous solution, so a sweep of rates
/// costs far less than the same number of cold starts.
class DetailBlock : public ICoupledElement
{
public:
    DetailBlock(polyMesh mesh, scalar nu, scalar rho, const Schemes& schemes,
                const SimpleControls& ctl, scalar probeMin, scalar probeMax);

    DetailBlock(const DetailBlock&) = delete;
    DetailBlock& operator=(const DetailBlock&) = delete;

    scalar lossAt(scalar qMag) override;

    const BlockMeasure& measure() const { return measure_; }
    const polyMesh& mesh() const { return mesh_; }
    const volVectorField& U() const { return U_; }
    const volScalarField& p() const { return p_; }

    label lastSweeps() const override { return measure_.sweeps; }
    bool lastConverged() const override { return measure_.converged; }

private:
    void buildSections();
    void applyBoundary(scalar qMag);
    void measure(const SimpleResult& res, scalar qMag);

    polyMesh mesh_;
    scalar nu_;
    scalar rho_;
    Schemes schemes_;
    SimpleControls ctl_;
    scalar probeMin_;
    scalar probeMax_;

    label inletPatch_ = -1;
    label outletPatch_ = -1;
    int axis_ = 0;
    /// +1 when the flow runs along +axis_, -1 along -axis_; taken from the
    /// inlet patch's outward normal, so every measured sign follows the mesh.
    int flowSign_ = 1;
    scalar inletArea_ = 0;

    /// Cross-section groups along the flow axis: cell indices sharing a centre
    /// coordinate, ordered the way the flow meets them, with the section's
    /// coordinate alongside.
    std::vector<std::vector<label>> sections_;
    std::vector<scalar> sectionX_;
    /// Axis-normal face groups at every cross-stream plane, patch planes
    /// included: the discrete continuity statement inside the region.
    std::vector<std::vector<label>> planes_;

    volVectorField U_;
    volScalarField p_;
    surfaceScalarField phi_;
    BlockMeasure measure_;
};

// ------------------------------------------------------------------
// The partitioned outer loop
// ------------------------------------------------------------------

/// One interface exchange: the network was solved with `resistance`, carried
/// `q`, and the partner reported `dpPa` for it. `residual` is the mismatch
/// between what the network charged and what the partner reported, relative to
/// the partner's number: exactly zero at the common solution.
struct OuterIteration
{
    label iteration = 0;
    scalar resistance = 0;
    scalar q = 0;
    scalar dpPa = 0;
    scalar residual = 0;
    label sweeps = 0;
    bool blockConverged = false;
    scalar networkResidual = 0;
};

struct CoupledControls
{
    label maxOuterIterations = 40;
    scalar relax = 1.0;    // interface relaxation; < 1 damps a strongly
                           // nonlinear (e.g. turbulent) partner
    scalar tolerance = 1e-9; // relative interface residual
};

struct CoupledResult
{
    std::vector<OuterIteration> history;
    SteadyNetworkResult network;
    label edge = -1;         // the pipe index that carries the block
    scalar q = 0;            // rate through the coupled element at that coefficient
    scalar dpPa = 0;         // the partner's drop at exactly that rate
    scalar resistance = 0;   // the coefficient the reported network was solved with
    scalar residual = 0;
    label iterations = 0;
    bool converged = false;
};

/// Gauss-Seidel on the interface: the 1D network is asked for the rate through
/// the coupled edge given the drop coefficient it currently carries, the
/// detailed model is then run at that rate and answers with its own drop, and
/// the coefficient is refreshed. The fixed point is the state where the drop
/// the network charges equals the drop the detailed model reports at the rate
/// the network delivers -- which is exactly the continuity-plus-momentum
/// statement at the interface.
///
/// For a laminar partner the reported drop is linear in the rate, so one
/// refresh lands on the fixed point; the relaxation and the residual are here
/// for the nonlinear ones.
///
/// `net` is mutated: the coupled pipe's `coupledResistance` ends at the
/// converged coefficient, so `result.network` and `net` describe the same state.
CoupledResult solveCoupled1d3d(OneDNetwork& net, label edge,
                               ICoupledElement& block,
                               const CoupledControls& ctl,
                               const SteadyControls& netCtl);

} // namespace coupling
} // namespace cf
