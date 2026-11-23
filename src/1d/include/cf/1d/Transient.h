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
    // Momentum-flux factors of the pipe's end profiles (5.70 / A3): the
    // area-averaged shape factor beta = <v^2>/<v>^2, >= 1 by Cauchy-Schwarz.
    // betaIn = 1 (flat) and betaOut = 1 are the defaults and leave the row
    // untouched. The row charges rho*(betaOut - betaIn)/A^2 * q|q| -- the same
    // account a coupled interface installs as its convective coefficient
    // (5.43), now declarable on a pipe whose own geometry (a section change at
    // a node, an imposed profile) owns the imbalance. On a held pipe the term
    // dies with the rate, exactly like the valve's (5.50): a static wall
    // carries no momentum flux to report.
    double betaIn = 1.0;       // shape factor of the upstream end profile
    double betaOut = 1.0;      // shape factor of the downstream end profile

    // 5.74 (A2): gelation and sub-threshold creep on a yielded pipe. All
    // three default to inert; without them the row is 5.49's rigid statement.
    double gelStrength = 0.0;  // Pa, tau_g_inf: extra static wall stress built
                               // up while the pipe rests (thixotropy asymmetry:
                               // flow clears the gel clock instantly)
    double gelTime = 0.0;      // s, the exponential rebuild constant t_g
    double creepRate = 0.0;    // 1/s, gamma_dot_c of the bi-viscosity
                               // regularisation: below the threshold the row
                               // creeps with mu_eff = mu + tau0/gamma_dot_c
                               // instead of holding Q = 0 exactly

    // Multi-segment valve schedule (opt-in, `valveSchedule { events {...} }`):
    // K(t) is the piecewise-linear interpolation through the event points
    // (t, K), before the first event the first value holds, after the last
    // the last value holds, and ADJACENT EQUAL-t EVENTS ARE A STEP (the value
    // jumps from the earlier point's to the later point's at that instant).
    // Non-empty = the table prices the valve and the single-cycle
    // valveClose/valveOpen keys (and the open-valve `valve` key, which the
    // table would leave unread) are refused at load time.
    bool valveEventsOn = false;
    std::vector<double> valveEventT; // s, ascending, duplicates = steps
    std::vector<double> valveEventK; // loss coefficients, each >= 0
};

/// Multi-segment source schedule (opt-in, top-level
/// `sourceSchedule { node <name>; events {...} }`): rate(t) through the event
/// points, same interpolation and step semantics as the valve table. The
/// targeted node must not also be priced by a `sources {}` entry.
struct TransientSourceSchedule
{
    std::string nodeName;
    std::vector<double> eventT;   // s, ascending, duplicates = steps
    std::vector<double> eventRate; // m3/s at each point
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

// ---- E7a (§5.177): opt-in transient wellbore equipment ---------------------
//
// HONEST BOUNDARIES (all three are declarative stand-ins, not actuator models):
//  * SSSV is a CONTROLLED shutoff element. It is a declared K(t) step/ramp that
//    ends in a sealed (Q = 0) momentum row, with a sealing latch after full
//    close. There is no valve actuator, stroke, hydraulic timing or failure
//    mode: the crew issues the close/reopen times, the element obeys them.
//  * Choke is a DECLARATIVE throttle. It is one loss coefficient K in the same
//    quadratic convention the row already uses for a valve, either declared
//    directly or converted from an equivalent orifice diameter and discharge
//    coefficient. No critical-flow branch, no multiphase slip, no erosion.
//  * Artificial lift is a SIMPLIFIED ESP dynamic. A quadratic head curve at
//    unit speed (affinity scaling H ~ N^2) plus a lumped first-order lag of the
//    shaft speed toward a commanded setpoint. It is neither a motor model nor
//    a pump pressure/velocity field: the lag is one state, the curve algebraic.
// All three are opt-in. A case declaring none of them advances the exact
// pre-E7a row, bit for bit. A pipe carries at most one wellbore element and
// none of them may share a pipe with a base `valve` declaration.
//
// Each element is a top-level block keyed by element name, e.g.
//   sssv  { S1 { pipe P1; openK 0; closedK 1e12; close { start 2; end 2; }
//                 reopen { start 8; end 8; } } }
//   choke { C1 { pipe P2; diameter 0.02; dischargeCoeff 0.62; } }
//   esp   { E1 { pipe P3; stages 4; shutoffHead 60; headAtRate 25; rate 0.01;
//                speed0 1; speedCmd 1.4; speedCmdTime 3; speedTau 1.5; } }

struct TransientSssvConfig
{
    std::string name;         // element name (the dictionary key)
    std::string pipeName;     // target pipe, resolved by name at load
    double openK = 0;         // loss coefficient while open, >= 0
    double closedK = 1e12;    // loss coefficient at full close, >= openK
    double closeStart = -1;   // s, close ramp start (== end is an instant step)
    double closeEnd = -1;     // s, close ramp end; from here the row is sealed
    double reopenStart = -1;  // s, reopen ramp start (seal clears); -1 = never
    double reopenEnd = -1;    // s, reopen ramp end
};

struct TransientChokeConfig
{
    std::string name;
    std::string pipeName;
    double k = 0;             // direct loss coefficient (valve K convention)
    double diameter = 0;      // equivalent orifice diameter, m (alt. to k)
    double dischargeCoeff = 0;// Cd in (0,1], required when diameter is declared
};

struct TransientEspConfig
{
    std::string name;
    std::string pipeName;
    double stages = 1;        // number of identical stages
    double shutoffHead = 0;   // H0, m of head at zero rate, unit speed
    double headAtRate = 0;    // Hq, m of head at rate `rate`, unit speed
    double rate = 0;          // Qr, m3/s where headAtRate is quoted
    double speed0 = 1;        // initial speed N/Nref, > 0
    double speedCmd = 1;      // commanded speed after speedCmdTime, > 0
    double speedCmdTime = 0;  // s, when the step command is issued, >= 0
    double speedTau = 0;      // s, first-order lag; <= 0 = instantaneous step
};

struct TransientCase
{
    FluidProperties fluid;
    std::vector<TransientNodeConfig> nodes;
    std::vector<TransientPipeConfig> pipes;
    std::vector<TransientSource> sources;
    // Opt-in multi-segment source schedules; each targets one node that no
    // `sources {}` entry also prices.
    std::vector<TransientSourceSchedule> sourceSchedules;
    // E7a opt-in wellbore equipment; each targets one pipe.
    std::vector<TransientSssvConfig> sssvs;
    std::vector<TransientChokeConfig> chokes;
    std::vector<TransientEspConfig> esps;
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

    /// The momentum-flux account pipe `k`'s row actually carries, whatever its
    /// source: the case's own betaIn/betaOut declaration (5.70) or the number
    /// a coupled interface installed. Pa when multiplied by q|q|.
    double pipeConvectiveOf(label k) const;

    /// The declared end-profile shape factors of pipe `k` (1.0 when the case
    /// declares nothing). The app reports the convective account only when
    /// these differ, so a case that asks for the momentum-flux term can be
    /// seen asking for it.
    double pipeBetaInOf(label k) const;
    double pipeBetaOutOf(label k) const;

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

    // ---- E7a: report the element a pipe carries and the values the last
    // solved step actually used, so a checker can re-derive them independently.
    // kind: 0 = none, 1 = SSSV, 2 = choke, 3 = ESP.
    /// 0 = no wellbore element on pipe k, 1 = SSSV, 2 = choke, 3 = ESP.
    int pipeWellboreKind(label k) const;
    /// The extra local loss coefficient (valve K convention) the element
    /// charged pipe k's last step: SSSV K(t) or the choke's K; 0 for ESP/none.
    double pipeWellboreK(label k) const;
    /// True when the SSSV on pipe k sealed the row (Q forced to 0) on the last
    /// solved step. False for every other kind.
    bool pipeWellboreSealed(label k) const;
    /// The ESP shaft speed N/Nref and head H (m) pipe k's last step used.
    /// Zero when pipe k carries no ESP.
    double pipeEspSpeed(label k) const;
    double pipeEspHead(label k) const;

    double time() const { return t_; }
    const std::vector<std::string>& nodeNames() const { return net_.nodeNames; }
    const std::vector<std::string>& pipeNames() const;

    const std::vector<double>& pressure() const { return p_; }
    std::vector<double> pipeFlows() const { return q_; }
    std::vector<double> liquidFractions() const;
    std::vector<double> hydrateFractions() const;
    std::vector<double> liquidVolumes() const;
    std::vector<double> hydrateVolumes() const;

    /// 5.68 / D5: the Courant rule over the whole network, read off the flows
    /// the last committed step left: target * min_k(L_k A_k / |q_k|) over
    /// MOVING pipes, clamped to [lo, hi]. A pipe at rest imposes nothing; a
    /// network at rest has no constraint at all and takes the ceiling hi.
    double cflStep(double target, double lo, double hi) const;

private:
    struct PipeDyn
    {
        NetworkPipe geo; // reuse the steady pipe geometry + friction law
        double qPrev = 0;
        double inertia = 0; // rho L / A, or the partner model's measurement of it
        // rho*(betaOut - betaIn)/A^2 (Pa.s2/m6), charged by a partner region whose
        // own profile changes shape across it. Multiplied by |Q|*Q in the row, so
        // it is the momentum-flux imbalance of the shape, not a friction. Zero
        // unless the case declares betaIn/betaOut (5.70) or a coupled interface
        // measures it and installs it.
        double convective = 0;
        // The declared end-profile shape factors, kept so the app and the
        // bindings can report what the case asked for (1.0 when it asked for
        // nothing).
        double betaIn = 1.0, betaOut = 1.0;
        // True for this step when a yielded pipe sits at rest under a drive the
        // static wall can hold: the row is then `Q = 0` and reads no wall law.
        // Latched at step entry from the committed state, so re-solving the same
        // step (the partitioned interface does that) cannot flip it mid-solve.
        bool held = false;
        // 5.74 (A2): the gel clock. Advanced by dt on every step the pipe ends
        // on the held/creep branch, cleared on the first moving step -- the
        // thixotropy asymmetry (rebuild slow, breakdown instant). The held
        // threshold is 4 (tau0 + tau_g_inf (1 - exp(-tRest/t_g))) L / d.
        double tRest = 0;
        // True when the case declared a sub-threshold creep rate: the held
        // branch then reads q = pi d^4 drive / (128 mu_eff L) instead of Q = 0.
        bool creeps = false;
        // The gel/creep declarations, copied from the config so the latch and
        // the held branch (both inside the step) can read them without the
        // config round trip.
        double gelStrength = 0, gelTime = 0, creepRate = 0;
        // The drive latched at step entry: the creep row is a forced value
        // like Q = 0, so it reads the latched drive, not the live iterate.
        double driveLatch = 0;
        // valve schedule: K ramps geo.valveK -> kClosed over
        // [closeStart, closeEnd] and back over [openStart, openEnd];
        // negative start times disable the ramps
        double closeStart = -1, closeEnd = -1, kClosed = 0;
        double openStart = -1, openEnd = -1;
        // The opt-in multi-segment table: non-empty means K(t) is read here
        // and the linear-ramp fields above are inert (refused at load time).
        bool valveEventsOn = false;
        std::vector<double> eventT, eventK;
    };
    // E7a: per-pipe wellbore element state. kind 0 leaves every field inert,
    // so an undeclared case carries no extra term and reads no extra state.
    struct WellboreDyn
    {
        int kind = 0; // 0 none, 1 SSSV, 2 choke, 3 ESP
        // choke K (constant), added to the row's quadratic local loss
        double extraK = 0;
        // SSSV declarations (valveKAt / sssvSealed read these)
        double openK = 0, closedK = 0;
        double closeStart = -1, closeEnd = -1;
        double reopenStart = -1, reopenEnd = -1;
        // ESP declarations (quadratic head + first-order speed lag)
        double stages = 1, h0 = 0, hq = 0, qr = 0;
        double speed0 = 1, speedCmd = 1, speedCmdTime = 0, speedTau = 0;
        // Values the last solved step actually used (also the report lines)
        double lastK = 0;
        bool lastSealed = false;
        double lastSpeed = 1;
        double lastHead = 0;
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
    // E7a: one entry per pipe, parallel to pipes_ (kind 0 when undeclared).
    std::vector<WellboreDyn> well_;
    std::vector<NodeDyn> nodes_;

    // Per-source schedule record, one entry per TransientSource in the case.
    // The owning node's `source` is re-summed from these every step() so the
    // rate tracks the schedule; `fraction` is constant and applied once here.
    struct SourceEntry {
        std::size_t nodeIndex = 0;
        double rate = 0, openRate = 0, fraction = 1;
        double start = -1, end = -1;
        double reopenStart = -1, reopenEnd = -1;
        // The opt-in multi-segment table: non-empty means rate(t) is read
        // here and the linear-ramp fields above are inert (a node priced by
        // a table and by a sources{} entry is refused at load time).
        bool eventsOn = false;
        std::vector<double> eventT, eventRate;
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

/// Advance the case from t = 0 to endTime in ceil(endTime/dt) steps of at most
/// dt -- the last step is clamped to the remainder, so endTime need not be a
/// multiple of dt (5.66 / C4). A sample is emitted every reportEvery steps AND
/// at the final step (sample 0 is the initial state); when endTime IS a
/// multiple of dt the step sizes and the sample set are exactly what the old
/// fixed-step loop produced.
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
/// A pipe may instead carry `valveSchedule { events { { t ..; k ..; } ... } }`
/// -- a multi-segment K(t) table (linear between events, equal adjacent times
/// a step, holds outside) that replaces the single-cycle keys; and the case
/// may carry a top-level `sourceSchedule { node <n>; events { { t ..;
/// rate ..; } } }` for a multi-segment rate(t) on one node. Declaring a table
/// beside the keys it replaces is refused loudly.
/// A pipe may instead carry `coupledResistance <R>` (Pa/(m3/s)), which hands
/// its loss to a partner model - see cf::coupling::solveTransientCoupled1d3d.
/// E7a adds three opt-in top-level blocks, each keyed by element name and each
/// targeting one pipe by `pipe <name>`:
///   sssv  { S1 { pipe P1; openK 0; closedK 1e12;
///                close { start 2; end 2; } reopen { start 8; end 8; } } }
///   choke { C1 { pipe P2; k 500; } }        // or `diameter ..; dischargeCoeff ..;`
///   esp   { E1 { pipe P3; stages 4; shutoffHead 60; headAtRate 25; rate 0.01;
///                speed0 1; speedCmd 1.4; speedCmdTime 3; speedTau 1.5; } }
/// At most one element per pipe, none on a pipe carrying a base `valve` key or
/// a `coupledResistance`; anything else is refused at load time. The three
/// blocks are declarative stand-ins (see the honest-boundary note above), not
/// actuator/motor models; undeclared cases are byte-identical to pre-E7a.
TransientCase loadTransientFromDict(const Dictionary& caseDict);
struct TransientTimeControls
{
    double dt = 0.01;
    double endTime = 1.0;
    int reportEvery = 1;
    // 5.68 / D5: cfl > 0 switches the transient loop to CFL-adaptive stepping.
    // Zero (the key absent) keeps the fixed-step march of 5.66.
    double cfl = 0.0;
    double dtMin = 0.0; // adaptive floor; loader defaults it to dt
    double dtMax = 0.0; // adaptive ceiling; loader defaults it to endTime
};

/// Same, driven by the whole TransientTimeControls. When controls.cfl > 0 the
/// loop is CFL-adaptive (5.68 / D5): the next dt is the Courant rule over the
/// whole network, dt = clamp(cfl * min_k(L_k A_k / |q_k|), dtMin, dtMax), read
/// off the flows the step just committed -- a pipe at rest imposes nothing,
/// and a network at rest takes dtMax. t accumulates (with varying dt there is
/// no step index to position one), the last step is clamped to
/// endTime - t -- and t snaps to endTime exactly when that clamp fires.
/// controls.cfl == 0 delegates to the fixed-step overload unchanged.
TransientSeries runTransient(const TransientCase& c,
                             const TransientTimeControls& controls);

TransientTimeControls loadTransientControls(const Dictionary& caseDict);

} // namespace cf
