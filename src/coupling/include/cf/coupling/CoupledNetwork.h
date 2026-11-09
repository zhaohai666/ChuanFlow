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

/// What an element's inlet carries beyond the number of cubic metres per second
/// it was asked for: the shape the interface declared, that rate read back as
/// the section integral of the condition the element actually got, the
/// momentum-flux factor `<ux^2>/<ux>^2` off the same faces, and the same factor
/// measured on the solved field at the first and last cross-section of the span
/// whose drop the element reports. An analytic pipe law has no section to
/// distribute over and answers with the defaults. A block the case left `flat`
/// answers 1 for `betaImposed` -- a plug is a plug -- but not for the other two:
/// its own field is already fuller than the plug at its first cross-section and
/// fuller still at its last (1.0959 then 1.1956 on the transient profile case).
struct ProfileLedger
{
    std::string shape = "flat";
    scalar fluxImposed = 0;
    scalar betaImposed = 1;
    scalar betaIn = 1;
    scalar betaOut = 1;
    scalar inletArea = 0;
};

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

    /// Called once at the start of every time level - before its first pass and
    /// never inside one - so a partner that carries its own history can take it
    /// there: the field as the interface left it at the accepted level is then
    /// the `u^n` every sweep of this level accelerates away from, no matter how
    /// many times the level gets re-solved. `qPrev` is the committed rate
    /// through this edge at the start of the level, so the partner does not have
    /// to read the network to know what it is stepping from.
    ///
    /// The default ignores the level, which is every analytic law and every
    /// quasi-steady block: they are asked for a drop at an *instant* and have no
    /// history to take, and that is exactly why the interface has been able to
    /// keep their region's momentum in its own row instead of in the coefficient
    /// (see `solveTransientCoupled1d3d`).
    virtual void beginTimeLevel(double dt, scalar qPrev)
    {
        (void)dt;
        (void)qPrev;
    }

    /// Whether the drop the last `lossAt` reported is the answer of a momentum
    /// balance that carries the region's own `du/dt`. The interface then splits
    /// that drop into the friction of the rate and the acceleration of the
    /// region's history, iterates the resistance refresh on the friction part
    /// alone, and grades its own momentum row, `R*q + M*(q - qPrev)/dt`, against
    /// the whole drop. Refreshing on the whole drop instead would put the
    /// region's `M*dq/dt` into the coefficient, where the 1D row already carries
    /// an inertance of its own, and the map the interface iterates would stop
    /// having a positive crossing whenever the step is shorter than `M/R`.
    ///
    /// False - every quasi-steady partner, every analytic law, every steady
    /// exchange - keeps the whole drop as the chord, which is what every case
    /// until now ran.
    virtual bool unsteadyDrop() const { return false; }

    /// What the region's own acceleration was worth in the drop the last
    /// `lossAt` reported, in Pa and signed with the acceleration: the part of
    /// that drop that is not the friction of the rate. Zero for every
    /// quasi-steady partner, whose whole drop is friction by construction.
    virtual scalar lastInertialDrop() const { return 0; }

    /// The region's hydraulic inertance `rho*Volume/Area^2` in kg/m4 -- the same
    /// quantity the 1D side calls `rho*L/A`, measured on the geometry that
    /// actually holds the fluid. Negative means the partner states no region of
    /// its own, which is every analytic law and every steady case: there the
    /// 1D pipe's nominal bore is the only description of the volume, so it is
    /// also the only source of the number.
    virtual scalar regionInertance() const { return -1; }

    /// The same quantity over only the part of the region the reported drop is
    /// measured across. A partner that states no region states no slab either.
    virtual scalar probeInertance() const { return -1; }

    /// Diagnostics of the last `lossAt`, recorded into the interface history so
    /// a regression can tell "the partner needed more sweeps" from "the
    /// interface stalled". An analytic law is exact at the rate it was given,
    /// so the defaults say one cost-free, converged evaluation.
    virtual label lastSweeps() const { return 0; }
    virtual bool lastConverged() const { return true; }

    /// What the element's inlet was handed besides the rate, as the element
    /// itself measures it. The default is the flat distribution: an analytic law
    /// has no section to distribute over, and a block the case left unnamed
    /// spreads its rate evenly across the inlet faces it wrote.
    virtual ProfileLedger profileLedger() const { return {}; }
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

    /// What the region's own momentum is worth, measured on the mesh: the
    /// hydraulic inertance `rho*V/A^2` (kg/m4, the same quantity as the 1D
    /// side's `rho*L/A` but taken from the geometry that holds the fluid).
    scalar inertance = 0;

    /// Unsteady partner only (`unsteady` true): the block solved a momentum
    /// balance that carries its own `du/dt`, so its reported drop is now worth
    /// more than the friction of the rate it was handed, and the extra is
    /// measured rather than assumed.
    bool unsteady = false;
    scalar dpFrictionPa = 0;    // dpPa minus the acceleration below: what the
                                // network is charged as a resistance
    scalar dpInertialPa = 0;    // rho/A * sum V_c (u_x - u_x_old)/dt over the
                                // slab the drop is measured across, signed along
                                // the flow: the force the region's own history
                                // asks for
    scalar dqPrev = 0;          // m3/s between the rate just imposed and the
                                // committed rate the level started from
    scalar inertanceDynamic = 0; // dpInertialPa / (dqPrev/dt), kg/m4: the
                                 // inertance the block's own acceleration is
                                 // worth, measured from the field it solved
    scalar betaProfile = 0;     // <u_x^2>/<u_x>^2 over the probe cross-sections,
                                // volume-weighted: the momentum-flux correction
                                // the flat-rate 1D picture does not see
    scalar dtLevel = 0;         // s, the level this solve belongs to
    scalar inertanceProbe = 0;  // rho*V_probe/A^2, kg/m4: the inertance of the
                                // slab the drop above is measured across, which
                                // is what the acceleration can belong to

    /// What the interface handed the inlet, and what the region's own field says
    /// about its shape downstream of it. A rate is one number; a rate spread over
    /// a section is a number and a shape, and the shape is what decides how far
    /// the region has to run before it is developed. `betaImposed` describes the
    /// condition as applied (it is the identity the mapping has to satisfy), the
    /// stations describe what the solve did with it.
    std::string inletShape;     // "flat" | "parabolic" | "power1over7"
    scalar fluxImposed = 0;     // m3/s, the integral of the imposed profile over
                                // the inlet section: the rate the shape carries
    scalar betaImposed = 0;     // <ux^2>/<ux>^2 over the inlet faces as imposed;
                                // 1.0 exactly for a flat profile
    scalar momentumFluxImposedPa = 0; // rho*beta*q^2/A, Pa: the streamwise
                                // momentum the section carries per unit area
    scalar betaIn = 0;          // the same factor measured on the solved field at
    scalar betaOut = 0;         // ... the first and the last cross-section
    /// (x, beta) down the region, first and last station always present: the
    /// shape relaxing toward the developed value, station by station.
    std::vector<std::pair<scalar, scalar>> betaStations;
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
                const TurbulenceConfig& turb,
                const std::string& inletShape = "flat");

    DetailBlock(const DetailBlock&) = delete;
    DetailBlock& operator=(const DetailBlock&) = delete;

    scalar lossAt(scalar qMag) override;

    /// Takes the level's history when the case asks the region to carry a real
    /// `du/dt` (`schemes{ ddt Euler; }`): the field as it stands is the previous
    /// level's answer, and it is copied once here so that every pass of this
    /// level - however many the interface needs - accelerates away from the same
    /// `u^n`. A laminar-or-steady-scheme block ignores the call.
    void beginTimeLevel(double dt, scalar qPrev) override;
    bool unsteadyDrop() const override { return unsteady_; }

    /// The acceleration the region's own field was solved against: the volume
    /// integral of `rho*du/dt` over the slab the drop is measured across, per
    /// unit section. Zero for a quasi-steady block, which has no history.
    scalar lastInertialDrop() const override { return measure_.dpInertialPa; }

    /// `rho*V/A^2` of the region, measured on the block's own mesh: the
    /// inertance the 1D momentum row would otherwise have to guess from the
    /// nominal bore of the pipe standing in for this region.
    scalar regionInertance() const override { return inertance_; }

    /// `rho*V/A^2` over the probe slab alone: the momentum of exactly the fluid
    /// whose pressure difference `lossAt` reports.
    scalar probeInertance() const override { return inertanceProbe_; }

    const BlockMeasure& measure() const { return measure_; }
    const polyMesh& mesh() const { return mesh_; }
    const volVectorField& U() const { return U_; }
    const volScalarField& p() const { return p_; }

    label lastSweeps() const override { return measure_.sweeps; }
    bool lastConverged() const override { return measure_.converged; }

    /// The shape this block was handed, the rate and the factor read back off
    /// the inlet faces it wrote, and the factor its own field carries at the two
    /// ends of the span whose drop it reports -- all straight out of the last
    /// measure.
    ProfileLedger profileLedger() const override
    {
        ProfileLedger ledger;
        ledger.shape = measure_.inletShape;
        ledger.fluxImposed = measure_.fluxImposed;
        ledger.betaImposed = measure_.betaImposed;
        ledger.betaIn = measure_.betaIn;
        ledger.betaOut = measure_.betaOut;
        ledger.inletArea = measure_.inletArea;
        return ledger;
    }

    /// The eddy-viscosity model the block solves with, or null when it is
    /// laminar. Lives as long as the block, and holds references to its fields.
    const RansModel* turbulence() const { return turbModel_.get(); }

private:
    void buildSections();
    void buildInletProfile();
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

    /// The shape the interface hands the inlet. `inletShape_ == "flat"` keeps the
    /// condition the coupling has always written -- one vector over the whole
    /// patch -- and any other shape writes a value per face instead, scaled so
    /// that the section integral is exactly the rate the network delivered.
    /// The weights are taken once from the mesh: they do not depend on the rate.
    std::string inletShape_ = "flat";
    /// The transverse axes this section is confined by no-slip walls. A profile
    /// varies across a wall, not across a symmetry side -- varying across the
    /// latter would ask the mesh for shear where the condition forbids it -- so
    /// the shape is applied along these axes and flat along every other one.
    std::vector<int> wallAxes_;
    std::vector<label> inletFaces_;
    std::vector<scalar> inletWeights_;   // u_f / (q/A), area-weighted mean 1

    /// `rho*V/A^2` of the region, measured on the mesh (kg/m4). The formula
    /// assumes a constant section, which `buildSections` checks rather than
    /// trusts: a region whose planes differ in area cannot be reduced to one
    /// inertance, and the block says so instead of reporting a wrong one.
    scalar inertance_ = 0;
    /// The same quantity over the probe slab only, `rho*V_probe/A^2`: the
    /// inertance of the fluid whose pressure difference is actually reported.
    scalar inertanceProbe_ = 0;

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

    /// The level the block is currently solving, and the field it accelerates
    /// away from. `unsteady_` is the case's own statement (`schemes ddt`), read
    /// once here, and it only becomes true at the first `beginTimeLevel`: a
    /// block asked to step in time but driven at an instant is a contradiction
    /// the block refuses rather than a quasi-steady answer it never claimed.
    bool unsteady_ = false;
    double levelDt_ = 0;
    scalar levelQPrev_ = 0;
    Field<Vector> UOldLevel_;
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
    /// The inertance the 1D momentum row carries for this edge (kg/m4), and so
    /// the region's share of the ramp's time constant. 0 where no time level is
    /// being solved -- every steady exchange.
    scalar inertance = 0;

    /// Unsteady partner only: what the partner reported splits into the part the
    /// rate's friction is worth and the part the region's acceleration is worth,
    /// `dpFrictionPa + dpInertialPa == dpPa`, both measured on the block's own
    /// field. `residual` is then the mismatch against the friction part, which is
    /// what the resistance refresh and the interface accelerator iterate on; the
    /// region's history is charged to the row, not to the coefficient, and the
    /// disagreement between the two sides' momentum is `rowError` below.
    bool unsteady = false;
    scalar dpFrictionPa = 0;
    scalar dpInertialPa = 0;
    scalar qPrev = 0;      // the committed rate this level accelerates from
    scalar dq = 0;         // q - qPrev, signed with the magnitudes both sides
                           // work in
    /// What the 1D momentum row charges for this edge at the rate the step
    /// delivers - `resistance*q + inertance*dq/dt` - against the drop the
    /// partner's own momentum balance reports. `rowError` is the relative
    /// disagreement. For a quasi-steady partner it stays 0 because there is
    /// nothing on the 3D side to compare the row against: the region's inertance
    /// is then a convention the 1D side declares and can only be read in the
    /// path it takes. For a partner that carries its own `du/dt` the row is a
    /// claim the two sides are graded on, level by level, and no amount of
    /// interface iteration changes the grade.
    scalar rowDropPa = 0;
    scalar rowError = 0;

    /// The shape this edge's rate arrived in, and what it is worth. `betaImposed`
    /// is the declared profile read back off the inlet faces; `betaIn` and
    /// `betaOut` are the same factor measured on the solved field at the first
    /// and last cross-section the block reports between. A plug inlet arrives
    /// short of the region's own profile and the two numbers differ in one
    /// direction, a Poiseuille inlet arrives fuller than the span it is handed
    /// and they differ in the other -- so the pair says which way the region had
    /// to move, which is more than a single imposed rate can say. What that
    /// movement is worth in the reported drop is not measured here; the only
    /// shape quantity priced by a number this interface solves is the imposed
    /// one, `rho*betaImposed*q^2/A`. An edge the case left flat reads 1 for
    /// `betaImposed` everywhere.
    /// `fluxImposed` is the section integral the inlet condition really carries,
    /// so it equals `q` whenever a shape was traded for a shape and not for a
    /// different amount of fluid.
    std::string inletShape = "flat";
    scalar fluxImposed = 0;
    scalar betaImposed = 1;
    scalar betaIn = 1;
    scalar betaOut = 1;
    scalar inletArea = 0;
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

    /// Where the detailed region's inertance comes from. Asked of the network
    /// (false, the default) the momentum row uses the coupled pipe's own
    /// `rho*L/A`, which is only as good as the bore the case typed in for the
    /// region. Asked of the partner (true) every coupled edge takes the
    /// `rho*V/A^2` its block measured on the region's actual geometry instead,
    /// and the term stays implicit in the same row.
    ///
    /// The two agree exactly when the pipe's nominal bore is the region's
    /// equal-area one, and disagree by the shape error of every other bore --
    /// which is the point: an irregular region has no single bore to name. That
    /// error moves the path, not the destination: the operating point carries no
    /// inertia, so both readings relax to the same rate while their first steps
    /// differ by the inertial split. What does NOT work is putting the region's
    /// momentum into the interface
    /// coefficient instead, where the partner would report
    /// `R(q)*q + M*(q - qPrev)/dt`: the coefficient then reads `a - c/R`, and a
    /// region whose rate is falling through a step reports a non-positive one
    /// whenever `dt` is below its own time constant `M/R` -- which is every
    /// transient that resolves the ramp it is iterating on.
    bool inertiaFromPartner = false;

    /// The same decision made per edge, indexed exactly like the caller's
    /// `edges`: 0 keeps that edge on the network's own `rho*L/A`, 1 asks that
    /// edge's partner for its `rho*V/A^2`, 2 asks it for the inertance of only
    /// the slab its reported drop is measured across, `rho*V_probe/A^2`. Empty
    /// (the default) leaves every edge to the switch above, so a run that names
    /// one source for the whole interface reads the same as before.
    ///
    /// Source 2 is a different question from source 1, and it is only a question
    /// at all once the partner carries its own `du/dt`: the drop handed to the
    /// interface is measured between two probe sections, so the momentum that
    /// belongs to it is the slab's, while the duct the pipe stands for is longer
    /// than the slab by whatever the probe leaves out. A quasi-steady partner can
    /// report both numbers but cannot say which one its drop is worth, so only a
    /// time-accurate one can be graded on the choice.
    ///
    /// Mixing is not a cosmetic generalisation of that switch. The two sources
    /// are two different measurements of two different things, and a network
    /// whose edges are described by drawings of uneven quality is the normal
    /// industrial case: one region is a machined duct the block can measure, the
    /// next is a drawing whose nominal bore the case typed in. Per-edge is what
    /// lets one run say "trust the block here, trust the drawing there" - and it
    /// is falsifiable per edge, because each edge's momentum row then answers to
    /// its own source while its neighbour's choice is free to disagree.
    std::vector<int> inertiaFromPartnerPerEdge;
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
/// A partner that is only asked at an *instant* - every quasi-steady block and
/// every analytic law, which is every case until the region's own `du/dt` was
/// wanted - reports a drop that depends on nothing but the rate, and that is what
/// keeps the interface stable. A time-accurate partner would report its region's
/// inertial term `M*(q - qPrev)/dt` inside that drop, so the coefficient the
/// network solves for would read `a - c/R` - and a region whose rate falls through
/// a step hands the network a non-positive one whenever `dt` is below its own time
/// constant `M/R`, which is every transient that resolves the ramp it iterates on.
/// The interface avoids the same failure from the other side: it refreshes on the
/// friction part of what the partner reports, so the coefficient keeps depending on
/// the rate alone, and charges the acceleration part to the momentum row, where an
/// inertance already sits.
///
/// `inertiaFromPartner` - or `inertiaFromPartnerPerEdge`, which makes the same
/// decision edge by edge inside one run - decides where that row's number comes
/// from: the coupled pipe's nominal `rho*L/A`, the block's `rho*V/A^2` over the
/// whole region, or the block's `rho*V_probe/A^2` over only the slab the reported
/// drop is measured across. When the bore is the equal-area one and the source is
/// the region's, the first two are the same number and the trajectory is identical
/// to machine precision; otherwise they differ by exactly the error a bore carries
/// for a region that has no bore, and for a quasi-steady partner that error lives
/// in the path (the first step's inertial split, the ramp's time constant) rather
/// than in the operating point both relax to.
///
/// The split is falsifiable on its own as well: the ramp's time constant is
/// `(M_X + M_By)/(m_X R_X + m_By R_By)`, so a region whose inertance is read
/// wrong by a shape factor ramps wrong by that factor's share of the loop. With
/// `schemes{ ddt Euler; }` the block also states which inertance its own
/// acceleration was worth (`inertanceDynamic`, against the installed row), so the
/// choice of source stops being a convention and becomes a graded claim - one no
/// amount of interface iteration can change.
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
