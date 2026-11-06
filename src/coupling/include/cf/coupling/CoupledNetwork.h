#pragma once

#include <memory>
#include <string>
#include <vector>

#include "cf/1d/Network.h"
#include "cf/1d/Transient.h"
#include "cf/common/Dictionary.h"
#include "cf/common/Labels.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/PolyMesh.h"
#include "cf/fvm/SimpleSolver.h"
#include "cf/turbulence/TurbulenceModel.h"

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

    /// Turbulence the block last solved with: zeros for a laminar block, and
    /// otherwise the numbers that say whether the wall layer is inside the
    /// wall-function band the model is valid in.
    std::string turbModel;    // "" when laminar
    scalar nutMax = 0;        // m2/s, largest eddy viscosity in the region
    scalar wallYPlusMax = 0;  // first cell's y+, max over the wall faces
    scalar wallUTau = 0;      // m/s, area-weighted mean friction velocity
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
///
/// `turb` is the block's own turbulence configuration: a `kEpsilon` or
/// `kOmegaSST` model makes the reported drop a genuine power law of the rate
/// (`dp/q ~ q^(m-1)`), which is what turns the interface into a nonlinear
/// fixed point instead of the one-shot contraction a laminar duct gives.
/// `laminar` -- what a case without a `turbulence{}` block asks for -- gives up
/// the model entirely and solves exactly the equation it always did.
class DetailBlock : public ICoupledElement
{
public:
    DetailBlock(polyMesh mesh, scalar nu, scalar rho, const Schemes& schemes,
                const SimpleControls& ctl, scalar probeMin, scalar probeMax,
                const TurbulenceConfig& turb);

    DetailBlock(const DetailBlock&) = delete;
    DetailBlock& operator=(const DetailBlock&) = delete;

    scalar lossAt(scalar qMag) override;

    const BlockMeasure& measure() const { return measure_; }
    const polyMesh& mesh() const { return mesh_; }
    const volVectorField& U() const { return U_; }
    const volScalarField& p() const { return p_; }

    label lastSweeps() const override { return measure_.sweeps; }
    bool lastConverged() const override { return measure_.converged; }

    /// The eddy-viscosity model the block solves with, or null when it is
    /// laminar. Lives as long as the block, and holds references to its fields.
    const RansModel* turbulence() const { return turbModel_.get(); }

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
    TurbulenceConfig turb_;

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
    /// Declared after the fields it references so it dies before them.
    std::unique_ptr<RansModel> turbModel_;
};

// ------------------------------------------------------------------
// The partitioned outer loop
// ------------------------------------------------------------------

/// One edge's side of one interface exchange: the network was solved while the
/// pipe carried `resistance`, delivered `q`, and the partner reported `dpPa` for
/// exactly that rate. `residual` is the mismatch between what the network
/// charged and what the partner reported, relative to the partner's number:
/// zero for that edge at the common solution.
struct EdgeExchange
{
    label edge = -1;
    std::string name;
    scalar resistance = 0;
    scalar q = 0;
    scalar dpPa = 0;
    scalar residual = 0;
    label sweeps = 0;
    bool blockConverged = false;
};

/// One pass of the steady interface: every coupled edge evaluated against the
/// same network solve.
struct OuterIteration
{
    label iteration = 0;
    std::vector<EdgeExchange> edges;
    scalar networkResidual = 0;
};

struct CoupledControls
{
    label maxOuterIterations = 40;
    scalar relax = 1.0;    // interface relaxation; < 1 damps a strongly
                           // nonlinear (e.g. turbulent) partner
    scalar tolerance = 1e-9; // relative interface residual, per edge

    /// Interface quasi-Newton over the whole coefficient vector (Aitken/IQN).
    /// Off by default: for a linear partner one refresh is already the fixed
    /// point, and the plain step is exactly that refresh. The accelerated step
    /// is for partners whose drop is not proportional to their rate - for those,
    /// the plain step's interface map can have a derivative above 1 and cannot
    /// converge at all, however many passes it is given.
    bool accelerate = false;
    label qnHistory = 8;      // secant pairs the operator is built from
    scalar qnLimit = 100.0;   // safeguard: max |step| / |residual|
};

struct CoupledResult
{
    std::vector<OuterIteration> history;
    SteadyNetworkResult network;
    std::vector<label> edges;             // pipe indices, in interface order
    std::vector<EdgeExchange> state;      // the reported network's per-edge state
    scalar residual = 0;                  // the worst edge: what the stop used
    label iterations = 0;
    bool converged = false;
};

/// Gauss-Seidel on the interface: the 1D network is asked for the rate through
/// each coupled edge given the drop coefficients it currently carries, every
/// detailed model is then run at its own rate and answers with its own drop, and
/// the coefficients are refreshed. The fixed point is the state where, at every
/// edge, the drop the network charges equals the drop the detailed model reports
/// at the rate the network delivers - which is exactly the
/// continuity-plus-momentum statement at that interface.
///
/// All partners are evaluated against the *same* network solve, so one pass
/// costs one network solve plus one solve per detailed region; the edges are
/// refreshed together and the pass stops on the worst edge, which is what makes
/// a second region actually couple rather than run in its own loop.
///
/// For a laminar partner the reported drop is linear in the rate, so one
/// refresh lands on the fixed point; the relaxation and the residual are here
/// for the nonlinear ones.
///
/// `net` is mutated: every coupled pipe ends at its converged coefficient, so
/// `result.network` and `net` describe the same state. `blocks` is used, not
/// owned, and must be the same length as `edges`.
CoupledResult solveCoupled1d3d(OneDNetwork& net,
                               const std::vector<label>& edges,
                               std::vector<ICoupledElement*> blocks,
                               const CoupledControls& ctl,
                               const SteadyControls& netCtl);

// ------------------------------------------------------------------
// the transient partitioned loop
// ------------------------------------------------------------------

/// Controls for the interface exchanges inside one time step of a transient
/// run: the same statement as `CoupledControls`, applied repeatedly instead of
/// once.
struct TransientCoupledControls
{
    label maxSubIterations = 8;
    scalar relax = 1.0;
    scalar tolerance = 1e-4;
    /// As in `CoupledControls`, but the secant history is rebuilt inside every
    /// time step: the interface problem there depends on the step's own history
    /// term, so a pair from the previous step measures a different map.
    bool accelerate = false;
    label qnHistory = 8;
    scalar qnLimit = 100.0;
};

/// One interface exchange inside one time step: every coupled edge evaluated
/// against the same uncommitted step. `newtonResidual` is that 1D step's own
/// convergence, shared by the edges because it is one solve.
struct TransientInterfacePass
{
    label pass = 0;
    std::vector<EdgeExchange> edges;
    double newtonResidual = 0;
};

/// A time step every interface agreed on, with the whole 1D state at it, so a
/// caller can read the transient response without solving anything again.
struct TransientCoupledStep
{
    double t = 0;                    // the time the step arrives at
    label passes = 0;
    std::vector<EdgeExchange> edges; // the accepted pass, one per coupled edge
    std::vector<double> pressure;    // every node, at the accepted state
    std::vector<double> flowRate;    // every pipe, at the accepted state
    std::vector<TransientInterfacePass> history;
};

struct TransientCoupledResult
{
    std::vector<label> edges;
    std::vector<std::string> nodeNames;
    std::vector<std::string> pipeNames;
    std::vector<TransientCoupledStep> steps;
    std::vector<EdgeExchange> state;  // the last step's accepted per-edge state
    scalar residual = 0;              // the worst edge at that last step
};

/// Advance a transient 1D network whose coupled pipes are handed to detailed
/// models, exchanging rate for drop inside every time step: the step is solved
/// with the coefficients the interface currently carries, every partner is run
/// at the rate that step delivers through its own edge, and the coefficients are
/// refreshed from the drops reported - all without committing the step, so a
/// refresh re-solves the same step from the same start. The step commits once
/// every edge agrees, and the coefficients carry into the next step, where the
/// same statement has to hold at the new rates.
///
/// The partner is asked for a drop at an *instant*, which is what keeps the
/// interface stable. A time-accurate partner would report its region's
/// inertial term `rho*L/A*(q - q_prev)/dt` inside that drop, making the
/// coefficient the network solves for depend on the previous step's rate: the
/// map becomes `R -> a - c/R`, whose derivative at a stiff region is far above
/// 1, so plain Gauss-Seidel cannot converge. So the region's inertance stays on
/// the 1D side - implicit, in the same momentum row as the rest of the network
/// - and the partner contributes the resistance it resists at the rate it is
/// given. The split is falsifiable: the time constant of the resulting ramp is
/// `(M_X + M_By)/(R_X + R_By)`, so a system carrying the region's inertance
/// twice reads there - in the case that ships, where the region holds 43% of the
/// loop's inertance, as a 43% longer time constant.
///
/// A step any edge cannot agree on is fatal: continuing would commit a state
/// whose sides disagree, and the rest of the series would inherit it.
///
/// `net` is mutated: it ends at the last committed step, with every coupled pipe
/// carrying the coefficient that step was solved with. `blocks` is used, not
/// owned, and must be the same length as `edges`.
TransientCoupledResult solveTransientCoupled1d3d(TransientNetwork& net,
                                                 const std::vector<label>& edges,
                                                 std::vector<ICoupledElement*> blocks,
                                                 const TransientCoupledControls& ctl,
                                                 double dt, double endTime);

} // namespace coupling
} // namespace cf
