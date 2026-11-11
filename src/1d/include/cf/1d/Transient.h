#pragma once

#include <string>
#include <vector>

#include "cf/1d/Network.h"
#include "cf/common/Dictionary.h"
#include "cf/common/Labels.h"

namespace cf
{

/// M2 first increment: transient 1D network hydraulics with
///   * pipe inertance  rho*L/A * dQ/dt,
///   * node compliance V/K (K = mixture bulk modulus; K = 0 -> incompressible
///     port enforcing sum(q) = 0),
///   * Darcy-Weisbach friction (same Churchill law as the steady kernel),
///   * valve local loss dp = K rho |v| v / 2 (same K convention as the
///     steady kernel) with optional linear close/open schedules K(t),
///   * gravity head,
///   * or, for a pipe handed to a partner model (`coupledResistance`), that
///     partner's statement of the region's drop - while the region's inertance
///     stays here, because the partner is asked for a drop at an instant,
///     exactly as the steady coupling layer asks it (see src/coupling),
/// advanced by implicit Euler with a Newton solve on (Q, p), plus drift-flux
/// transport of the mobile liquid volume per node (slug tracking) and an
/// exact integrating-factor first-order hydrate kinetics substep.
/// Phase volumes per node: V = V_L (mobile liquid) + V_h (hydrated,
/// immobile) + V_G (gas fills the rest).

struct TransientNodeConfig
{
    std::string name;
    double volume = 0;           // m3, fixed geometric node volume
    double bulkModulus = 0;      // Pa, mixture bulk modulus (0 = incompressible port)
    bool pressureFixed = false;  // Dirichlet pressure
    double pressure = 0;         // fixed value, or initial value for free nodes
    double liquidFraction = 1.0; // initial mobile-liquid fraction of volume
    double hydrateFraction = 0.0;// initial hydrate fraction of total liquid volume
    double hydrateRate = 0.0;    // 1/s, first-order liquid -> hydrate rate
};

struct TransientPipeConfig
{
    std::string name;
    std::string fromName;
    std::string toName;
    double length = 0;
    double diameter = 0;
    double roughness = 0;
    double elevation = 0; // z(to) - z(from)
    double initialFlow = 0;
    double valveK = 0;           // open-valve loss coefficient dp = K rho v^2/2
    // linear schedule: K ramps valveK -> valveKClosed over
    // [valveCloseStart, valveCloseEnd], holds closed, then ramps back over
    // [valveOpenStart, valveOpenEnd]; negative start times disable the ramps
    double valveCloseStart = -1;
    double valveCloseEnd = -1;
    double valveKClosed = 0;
    double valveOpenStart = -1;
    double valveOpenEnd = -1;
    // Pa/(m3/s): hand this pipe's loss to a partner model, the same hook the
    // steady kernel's NetworkPipe::coupledResistance offers. The momentum row
    // then charges R*Q for the region instead of this pipe's own friction; the
    // region's inertance rho*L/A stays on the 1D side, so the partner only has
    // to state the drop it resists. A coupled pipe carries no valve.
    double coupledResistance = 0;
    // Rheology, the same numbers the steady kernel reads (Network.h): at
    // flowIndex > 0 the friction line is Churchill on the Metzner-Reed number
    // instead of rho u d / mu, through the one shared wall-law function.
    // `yieldStress` is the threshold the row's static branch holds against: at
    // rest with |drive| <= 4 tau_0 L / d the momentum row is `Q = 0` and no wall
    // law is read (5.49); that threshold needs a branch integrated in *time*,
    // which is why it is not the steady kernel's rate-bracket dead branch.
    double flowIndex = 0;      // n; 0 = Newtonian (the fluid's own mu)
    double consistency = 0;    // K, Pa.s^n; required when flowIndex > 0
    double yieldStress = 0;    // Pa; 0 = no threshold, the row always has a law
};

struct TransientSource
{
    std::string nodeName;
    double rate = 0;     // m3/s total injected volume (used when no schedule)
    double fraction = 1; // liquid fraction of the injected volume
    // Optional linear ramp schedule, same shape as valveClose / valveOpen on a
    // pipe: rate ramps from 0 to openRate over [start, end], holds at openRate,
    // then ramps back to 0 over [reopenStart, reopenEnd]. start < 0 disables
    // the schedule and `rate` is used for the whole run; reopenStart < 0 keeps
    // the source on at openRate after the ramp. Continuous in t.
    double openRate = 0;
    double start = -1, end = -1;
    double reopenStart = -1, reopenEnd = -1;
};

struct TransientCase
{
    FluidProperties fluid;
    std::vector<TransientNodeConfig> nodes;
    std::vector<TransientPipeConfig> pipes;
    std::vector<TransientSource> sources;
};

struct TransientSample
{
    double t = 0;
    std::vector<double> pressure;        // per node
    std::vector<double> liquidFraction;  // mobile V_L / V per node
    std::vector<double> hydrateFraction; // V_h / (V_L + V_h) per node
    std::vector<double> flowRate;        // per pipe
};

struct TransientSeries
{
    std::vector<std::string> nodeNames;
    std::vector<std::string> pipeNames;
    std::vector<TransientSample> samples;
};

/// dh/dt = k (1 - h) integrated exactly: h -> 1 - (1 - h0) exp(-k dt).
/// Hydrate volume grows at the expense of mobile liquid (volume preserving).
void hydrateKineticsStep(double& liquidVolume, double& hydrateVolume,
                         double rate, double dt);

class TransientNetwork
{
public:
    explicit TransientNetwork(const TransientCase& c);

    /// One implicit-Euler step of size dt; returns max scaled Newton residual.
    double step(double dt);

    /// Solve the implicit-Euler step of size dt from the committed state and
    /// leave it pending: nothing advances, so the caller may re-solve the same
    /// step from the same start as many times as it likes. This is what a
    /// partitioned interface needs to refresh a coupled coefficient within one
    /// time step. Returns max scaled Newton residual.
    double solveStep(double dt);

    /// Adopt the pending step: flows, pressures, liquid transport and the clock
    /// all move to the state `solveStep` solved. Fatal if none is pending.
    void commitStep();

    /// The pending step's state; fatal while no step is pending.
    double pendingPipeFlow(label k) const;
    double pendingPressure(label i) const;

    /// Set the partner coefficient of pipe `k` (index in `pipeNames()`). A
    /// pending step is dropped: it was solved with the old coefficient, so it
    /// is no longer a state of the system.
    void setCoupledResistance(label k, double resistance);

    /// Replace pipe `k`'s inertance `rho*L/A` with the one its partner model
    /// measured on its own grid. The momentum row keeps the term either way -
    /// only the number changes - so this is the one way a detailed region can
    /// state its momentum without putting it into the interface coefficient.
    /// Only a coupled pipe can be asked, and a pending step is dropped for the
    /// same reason `setCoupledResistance` drops it.
    void setCoupledInertia(label k, double inertia);

    /// The inertance (Pa.s2/m3, kg/m4) pipe `k`'s momentum row actually carries.
    double pipeInertia(label k) const;

    /// Charge pipe `k`'s momentum row with the momentum flux its partner region's
    /// own changing profile asks for: `coefficient*|Q|*Q` with the coefficient in
    /// Pa.s2/m6, so the term is worth Pa at the row's rate. Zero unbuckles the row
    /// back to `rho*L/A*dQ/dt + R*Q` -- and, unlike the inertance swap, this one
    /// is signed: a region that sheds momentum flux at its outlet pays the row,
    /// one that grows it bills the row.
    ///
    /// Only a coupled pipe can be asked, for the same reason `setCoupledInertia`
    /// can: a pipe with no partner has no measured profile on either of its
    /// stations, and its own `frictionDp` already carries whatever a nominal bore
    /// is worth. A pending step is dropped because the coefficient it was solved
    /// with is gone.
    void setCoupledConvective(label k, double coefficient);

    /// The convective coefficient (Pa.s2/m6) pipe `k`'s momentum row carries.
    double coupledConvectiveOf(label k) const;

    /// True while `pipe k` is handed to a partner model.
    bool pipeCoupled(label k) const;

    /// Pipe `k`'s partner coefficient (Pa/(m3/s)), or 0 when it is not coupled.
    double coupledResistanceOf(label k) const;

    /// The wall law pipe `k`'s momentum row actually reads: the declared flow
    /// index (0 = Newtonian, nothing declared) and its consistency K (Pa.s^n).
    /// The app reports the pair only when an index is declared, so a case that
    /// asks for a rheology can be seen asking for it -- before 5.48 the row
    /// computed `rho u d / mu` by hand and a declaration changed no number.
    double pipeFlowIndex(label k) const;
    double pipeConsistencyOf(label k) const;

    /// The constitutive law of pipe `k`'s partner region (empty/Newtonian when
    /// the pipe is not coupled or declared none). A coupled pipe's own
    /// wall-law fields are zero by construction: its loss is the partner's
    /// statement, and this is the declaration the coupling app forwards to
    /// the 3D block.
    Rheology partnerRheologyOf(label k) const;

    /// The drop pipe `k`'s static wall can hold: `4 tau_0 L / d` (Pa), 0 when no
    /// yield stress is declared. `pipeHeld` reports the latch the last solved
    /// step started with: true means that step replaced this pipe's momentum row
    /// with `Q = 0` and read no wall law at all. A held pipe is a *state* of the
    /// row, not of the fluid -- the same pipe under a super-threshold drive is
    /// released and steps on the moving-wall law.
    double pipeYieldThreshold(label k) const;
    bool pipeHeld(label k) const;

    double time() const { return t_; }
    const std::vector<std::string>& nodeNames() const { return net_.nodeNames; }
    const std::vector<std::string>& pipeNames() const;

    const std::vector<double>& pressure() const { return p_; }
    std::vector<double> pipeFlows() const { return q_; }
    std::vector<double> liquidFractions() const;
    std::vector<double> hydrateFractions() const;
    std::vector<double> liquidVolumes() const;
    std::vector<double> hydrateVolumes() const;

private:
    struct PipeDyn
    {
        NetworkPipe geo; // reuse the steady pipe geometry + friction law
        double qPrev = 0;
        double inertia = 0; // rho L / A, or the partner model's measurement of it
        // rho*(betaOut - betaIn)/A^2 (Pa.s2/m6), charged by a partner region whose
        // own profile changes shape across it. Multiplied by |Q|*Q in the row, so
        // it is the momentum-flux imbalance of the shape, not a friction. Zero
        // unless a coupled interface measures it and installs it.
        double convective = 0;
        // True for this step when a yielded pipe sits at rest under a drive the
        // static wall can hold: the row is then `Q = 0` and reads no wall law.
        // Latched at step entry from the committed state, so re-solving the same
        // step (the partitioned interface does that) cannot flip it mid-solve.
        bool held = false;
        // valve schedule: K ramps geo.valveK -> kClosed over
        // [closeStart, closeEnd] and back over [openStart, openEnd];
        // negative start times disable the ramps
        double closeStart = -1, closeEnd = -1, kClosed = 0;
        double openStart = -1, openEnd = -1;
    };
    struct NodeDyn
    {
        double volume = 0, bulk = 0;
        bool fixed = false;
        double pPrev = 0;
        double unk = -1;
        double vL = 0, vH = 0;
        double source = 0, sourceFraction = 0;
        double rate = 0;
    };

    OneDNetwork net_; // node names + resolved pipe endpoints (from/to)
    std::vector<std::string> pipeNames_;
    std::vector<PipeDyn> pipes_;
    std::vector<NodeDyn> nodes_;

    // Per-source schedule record, one entry per TransientSource in the case.
    // The owning node's `source` is re-summed from these every step() so the
    // rate tracks the schedule; `fraction` is constant and applied once here.
    struct SourceEntry {
        std::size_t nodeIndex = 0;
        double rate = 0, openRate = 0, fraction = 1;
        double start = -1, end = -1;
        double reopenStart = -1, reopenEnd = -1;
    };
    std::vector<SourceEntry> sourceEntries_;
    std::vector<double> p_;  // current node pressures (all nodes)
    std::vector<double> q_;  // current pipe flows
    FluidProperties fluid_;
    double t_ = 0;
    std::vector<double> xPending_; // solveStep's state, not yet committed
    double pendingDt_ = 0;
    bool pending_ = false;
};

/// Advance the case by endTime/dt steps, sampling every reportEvery steps
/// (sample 0 is the initial state).
TransientSeries runTransient(const TransientCase& c, double dt,
                             double endTime, int reportEvery);

/// Case grammar (see tests/cases/transientBlowdown/case.dict and
/// tests/cases/transientValve/case.dict):
///   fluid { rho 850; mu 0.1; }
///   nodes { TANK { volume 0.5; bulkModulus 2.0e8; pressure 5.0e5;
///                  liquidFraction 1; hydrateRate 0; }
///            OUT { fixedPressure 0; } }
///   pipes { P1 { from TANK; to OUT; length 100; diameter 0.05;
///                valve 20; valveClose { start 2; end 6; closedK 1.0e12; }
///                valveOpen { start 10; end 14; } } }
///   sources { s1 { node TANK; rate 0.0; fraction 1; } }
///   time { dt 0.01; endTime 10; reportEvery 1; }
/// A pipe may instead carry `coupledResistance <R>` (Pa/(m3/s)), which hands
/// its loss to a partner model - see cf::coupling::solveTransientCoupled1d3d.
TransientCase loadTransientFromDict(const Dictionary& caseDict);
struct TransientTimeControls
{
    double dt = 0.01;
    double endTime = 1.0;
    int reportEvery = 1;
};
TransientTimeControls loadTransientControls(const Dictionary& caseDict);

} // namespace cf
