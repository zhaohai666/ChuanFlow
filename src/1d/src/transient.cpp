#include "cf/1d/Transient.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

#include "cf/1d/SolverGates.h"
#include "cf/common/Error.h"

#include "DenseSolve.h"

namespace cf
{
namespace
{

constexpr double kGravity = 9.80665;

double pipeArea(double diameter)
{
    return 0.25 * M_PI * diameter * diameter;
}

/// The drop a yielded pipe's static wall can hold without moving: the same
/// `4 tau_0 L / d` the steady kernel puts in front of its rate bracket (5.46),
/// and the q = 0 limit of the moving-wall law. Geometry only -- the threshold
/// does not know the flow index, which is why 5.47 could state it as a limit.
double yieldThresholdPa(const NetworkPipe& geo)
{
    return 4.0 * geo.yieldStress * geo.length / geo.diameter;
}

/// A yielded pipe at rest whose drive does not reach the wall's own threshold
/// is held: the static wall supplies whatever the drive asks for, up to
/// `4 tau_0 L / d`, so the row has no friction law to read and becomes the
/// statement `Q = 0`. Note what is *not* here: the drive is compared with the
/// threshold, not with a rate tolerance, so the held set is decided by the same
/// geometry-only number 5.46 put in front of the steady bracket, and a fluid
/// held at rest under a super-threshold drive is a fluid that starts moving.
///
/// 5.74 (A2): with a gel declaration the threshold is not constant -- a wall
/// at rest for t_rest seconds holds `4 (tau0 + tau_g (1 - exp(-t_rest/t_g)))
/// L / d`, the thixotropy asymmetry being that the clock only runs at rest
/// (the caller clears it on the first moving step). Without a declaration
/// tau_g = 0 and the number is 5.49's own, bit for bit.
bool holdsStatic(const NetworkPipe& geo, double qPrev, double drive,
                 double tRest, double gelStrength, double gelTime,
                 bool wasHeld)
{
    if (!(geo.yieldStress > 0)) return false;
    const double tauG = gelStrength > 0
                          ? gelStrength * (1.0 - std::exp(-tRest / gelTime))
                          : 0.0;
    return (qPrev == 0 || wasHeld)
        && std::abs(drive)
               <= 4.0 * (geo.yieldStress + tauG) * geo.length / geo.diameter;
}

/// The drop a pipe's momentum row is offered: the node pair minus the static
/// head of its own rise, which is what a static wall has to resist. The sign is
/// the row's, not a second convention: `drivePa` is `dp - rho g dz`, the same
/// quantity the residual moves, so the latch and the row cannot disagree about
/// which end is upstream. Valve and convective terms are not in it -- both
/// vanish with the rate, so at the held point they ask nothing.
double drivePa(const NetworkPipe& geo, const FluidProperties& fl,
               double pFrom, double pTo)
{
    return pFrom - pTo - fl.rho * kGravity * geo.elevation;
}

/// Pressure drop consumed by friction for a signed flow: the wall law the pipe
/// declares, through the same pipeFrictionPa the steady kernel calls, made odd
/// in Q. Nothing is computed here -- that is the point. Until 5.48 this line
/// built `rho u d / mu` itself, so a `flowIndex` declared in a transient case
/// was ignored rather than refused: the row quietly kept Newtonian friction.
///
/// A yielded pipe at rest is the one rate where the wall law is not a function:
/// `dp_wall(Q)` tends to its threshold `4 tau_0 L/d` from above as Q -> 0+, and
/// the static wall can hold anything up to that number, so the value the row
/// needs at Q = 0 is the threshold with the sign of the drive it has to resist
/// (`drive`, the pressure side of the row). Zero drive asks for zero, which is
/// the held case -- and a held pipe's row is replaced by `Q = 0` before this
/// function is ever read, so the sign convention only ever matters when the
/// fluid is starting to move.
double frictionDp(const NetworkPipe& geo, const FluidProperties& fl,
                  double q, double drive)
{
    if (q == 0)
        return (geo.yieldStress > 0 && drive != 0)
                 ? std::copysign(yieldThresholdPa(geo), drive)
                 : 0.0;
    const double dp = pipeFrictionPa(geo, fl, std::abs(q));
    return q > 0 ? dp : -dp;
}

/// d(dpFriction)/dQ by central differencing on the monotone odd law;
/// falls back to the laminar slope 128 mu L / (pi D^4) near Q = 0.
double frictionSlope(const NetworkPipe& geo, const FluidProperties& fl,
                     double q, double drive)
{
    const double h = std::max(1e-5 * std::abs(q), 1e-12);
    if (geo.yieldStress > 0)
    {
        // The yielded law's slope at the rest point is infinite: any rate at
        // all lifts the wall stress back to the threshold, so the diagonal must
        // be taken from the rest value outward, not across it. Differencing
        // through Q = 0 would put the threshold jump itself into the quotient
        // (dp_wall(h) - dp_wall(-h) = 2 dp_wall(h), a slope of ~1e15) and Newton
        // would then sit at the rest point forever -- the fluid would never
        // start, which is the opposite error from the one 5.48 refused.
        const double a = std::abs(q);
        const double qa = a > 0 ? a : h;
        const double dp = pipeFrictionPa(geo, fl, qa);
        const double excess = dp - yieldThresholdPa(geo);
        return (std::isfinite(excess) && excess > 0) ? excess / qa : 0.0;
    }
    const double d = (frictionDp(geo, fl, q + h, drive)
                    - frictionDp(geo, fl, q - h, drive)) / (2.0 * h);
    if (geo.flowIndex > 0)
        // The Newtonian floor below is a number made out of mu. At n != 1 the
        // declared law's own low-rate slope is n * dp / q, which diverges for a
        // shear-thinning fluid and vanishes for a thickening one, so clamping it
        // to 128 mu L/(pi d^4) would damp the row with a resistance the fluid
        // does not have. The diagonal still carries inertia/dt, so a vanishing
        // slope cannot make the system singular.
        return (std::isfinite(d) && d > 0) ? d : 0.0;
    const double lam = 128.0 * fl.mu * geo.length
                     / (M_PI * geo.diameter * geo.diameter * geo.diameter * geo.diameter);
    return (std::isfinite(d) && d > 0) ? std::max(d, 0.2 * lam) : lam;
}

/// Valve loss coefficient at time t: linear ramp from the open K to kClosed
/// over [closeStart, closeEnd], hold closed, then linear ramp back over
/// [openStart, openEnd]. closeStart < 0 disables the schedule (K stays at
/// the open value); openStart < 0 keeps the valve closed after the close
/// ramp. The schedule is continuous in t, so samples that straddle a ramp
/// boundary never see a K jump.
double valveKAt(double kOpen, double kClosed,
                double cs, double ce, double os, double oe, double t)
{
    if (cs < 0) return kOpen;
    if (t < cs) return kOpen;
    if (t < ce) return kOpen + (t - cs) / (ce - cs) * (kClosed - kOpen);
    if (os < 0 || t < os) return kClosed;
    if (t < oe) return kClosed + (t - os) / (oe - os) * (kOpen - kClosed);
    return kOpen;
}

/// Source rate at time t, mirror of valveKAt with "open" = 0 (off) and
/// "closed" = openRate (on). The source is off before start, ramps up to
/// openRate over [start, end], holds at openRate, then ramps back to 0 over
/// [reopenStart, reopenEnd]. start < 0 disables the schedule and the caller
/// falls back to the constant `rate`. Continuous in t.
double sourceRateAt(double rate, double openRate,
                    double start, double end,
                    double reopenStart, double reopenEnd, double t)
{
    if (start < 0) return rate;
    return valveKAt(0.0, openRate, start, end, reopenStart, reopenEnd, t);
}

/// The value at time t of a piecewise-linear event table: times ascending,
/// EQUAL ADJACENT TIMES ARE A STEP (the value jumps from the earlier point's
/// to the later point's at that instant), before the first event the first
/// value holds and after the last the last value holds. One helper for both
/// the valve's K(t) and a source's rate(t): binary search the upper bound of
/// t, then interpolate the segment it closes. At a step instant the AFTER
/// value is taken -- the jump happens at the instant itself -- so a schedule
/// evaluated at the implicit new time answers for the state the step runs
/// in. The table is validated at load time (non-empty, t >= 0, ascending),
/// so the search below has a front and a back to stand on.
double eventTableAt(const std::vector<double>& ts,
                    const std::vector<double>& vs, double t)
{
    if (t < ts.front()) return vs.front();
    if (t >= ts.back()) return vs.back();
    const std::size_t hi = static_cast<std::size_t>(
        std::upper_bound(ts.begin(), ts.end(), t) - ts.begin());
    const std::size_t lo = hi - 1;
    const double span = ts[hi] - ts[lo];
    if (!(span > 0)) return vs[hi]; // t landed exactly on a step: after value
    return vs[lo] + (t - ts[lo]) / span * (vs[hi] - vs[lo]);
}

/// A coupled pipe's loss is its partner's, so no local loss may share it --
/// neither the single-cycle ramps nor the multi-segment event table.
void requireNoValve(const std::string& name, double kOpen, double closeStart,
                    double openStart, bool eventsOn)
{
    if (kOpen > 0 || closeStart >= 0 || openStart >= 0 || eventsOn)
        cfFatal("transient: pipe '" + name
                + "' is a coupled element: its loss is the partner model's, so"
                  " it cannot carry a valve");
}

// ---- E7a (§5.177) wellbore-element physics, all closed-form in t -----------

/// SSSV loss coefficient at time t: the same continuous ramps as a valve, so
/// a sample straddling a close/open boundary never sees a K jump. The seal is
/// what makes it a safety valve rather than a valve: once the close ramp ends,
/// the row is Q = 0 until the reopen ramp starts (the element is a control
/// declaration, not an actuator model).
double sssvKAt(double openK, double closedK,
               double cs, double ce, double os, double oe, double t)
{
    return valveKAt(openK, closedK, cs, ce, os, oe, t);
}

bool sssvSealed(double cs, double ce, double os, double t)
{
    return cs >= 0 && t >= ce && (os < 0 || t < os);
}

/// ESP shaft speed N/Nref at time t: before the command the initial speed
/// holds; at the command instant a step to speedCmd is issued; for a positive
/// lag the speed then relaxes with the exact integrating factor
///   N(t) = Ncmd + (N0 - Ncmd) exp(-(t - tCmd)/tau),
/// so a checker can reproduce it at any dt without marching. tau <= 0 is an
/// instantaneous step (N = Ncmd once t >= tCmd).
double espSpeedAt(double speed0, double speedCmd, double tCmd, double tau,
                  double t)
{
    if (t < tCmd) return speed0;
    if (!(tau > 0)) return speedCmd;
    return speedCmd
         + (speed0 - speedCmd) * std::exp(-(t - tCmd) / tau);
}

/// ESP head H (m) at rate q and speed N: affinity-scaled shutoff head
/// N^2 H0 and the same quadratic fall-off to Hq at rate Qr as the steady
/// pumpHeadOf (5.81), clipped at zero. Qr <= 0 means a flat N^2 H0 curve.
/// dH/dq = -stages * 2 (H0 - Hq) q / Qr^2, folded into the row's diagonal.
double espHeadAt(double stages, double h0, double hq, double qr,
                 double q, double n)
{
    if (!(h0 > 0)) return 0.0;
    const double n2 = n * n;
    const double h = (qr > 0)
        ? stages * (n2 * h0 - (h0 - hq) * q * q / (qr * qr))
        : stages * n2 * h0;
    return std::max(h, 0.0);
}

double espHeadSlope(double stages, double h0, double hq, double qr, double q)
{
    if (!(h0 > 0) || !(qr > 0)) return 0.0;
    return -stages * 2.0 * (h0 - hq) * q / (qr * qr);
}

} // namespace

void hydrateKineticsStep(double& liquidVolume, double& hydrateVolume,
                         double rate, double dt)
{
    if (rate <= 0 || dt <= 0) return;
    const double total = liquidVolume + hydrateVolume;
    if (total <= 0) return;
    // dh/dt = k (1 - h) with h = V_h/(V_L+V_h): the mobile liquid volume
    // obeys dV_L/dt = -k V_L, integrated exactly (multiplicative decay), so
    // a kinetic benchmark is reproducible at any dt.
    const double decay = std::exp(-rate * dt);
    const double newLiq = liquidVolume * decay;
    hydrateVolume = total - newLiq;
    liquidVolume = newLiq;
}

TransientNetwork::TransientNetwork(const TransientCase& c)
    : fluid_(c.fluid)
{
    if (!(c.fluid.rho > 0) || !(c.fluid.mu > 0))
        cfFatal("transient: fluid needs positive rho and mu");

    for (const TransientNodeConfig& n : c.nodes)
    {
        net_.addNode(n.name);
        NodeDyn d;
        d.volume = n.volume;
        d.bulk = n.bulkModulus;
        d.fixed = n.pressureFixed;
        d.pPrev = n.pressure;
        d.rate = n.hydrateRate;
        const double liq = n.liquidFraction * n.volume;
        if (!(liq <= n.volume + 1e-12) || n.liquidFraction < 0 || n.liquidFraction > 1)
            cfFatal("transient: node '" + n.name + "' liquidFraction out of [0,1]");
        d.vL = liq * (1.0 - n.hydrateFraction);
        d.vH = liq * n.hydrateFraction;
        if (!d.fixed && !(d.volume > 0) && d.bulk > 0)
            cfFatal("transient: free node '" + n.name + "' needs a volume");
        nodes_.push_back(d);
    }
    if (c.nodes.empty()) cfFatal("transient: case has no nodes");

    for (const TransientSource& s : c.sources)
    {
        const label i = net_.nodeIndex(s.nodeName);
        SourceEntry e;
        e.nodeIndex = static_cast<std::size_t>(i);
        e.rate = s.rate;
        e.openRate = s.openRate;
        e.fraction = s.fraction;
        e.start = s.start;
        e.end = s.end;
        e.reopenStart = s.reopenStart;
        e.reopenEnd = s.reopenEnd;
        sourceEntries_.push_back(e);
        // The fraction is constant; the rate is re-summed every step() so the
        // schedule can move it. Pre-fill with the value the schedule holds at
        // t = 0: the constant rate for an unscheduled source, 0 for a
        // single-cycle schedule (off before its ramp), the table's own value
        // for a multi-segment one.
        nodes_[e.nodeIndex].sourceFraction = e.fraction;
        nodes_[e.nodeIndex].source =
            e.eventsOn ? eventTableAt(e.eventT, e.eventRate, 0.0)
                       : ((e.start < 0) ? e.rate : 0.0);
    }

    // The opt-in multi-segment source schedules: one SourceEntry each, the
    // node resolved by name (the loader refuses a name that resolves to
    // nothing and a node priced twice). fraction is 1: the table injects
    // pure liquid, the same default a sources{} entry carries.
    for (const TransientSourceSchedule& ss : c.sourceSchedules)
    {
        const label i = net_.nodeIndex(ss.nodeName);
        SourceEntry e;
        e.nodeIndex = static_cast<std::size_t>(i);
        e.fraction = 1.0;
        e.eventsOn = true;
        e.eventT = ss.eventT;
        e.eventRate = ss.eventRate;
        sourceEntries_.push_back(e);
        nodes_[e.nodeIndex].sourceFraction = e.fraction;
        nodes_[e.nodeIndex].source = eventTableAt(e.eventT, e.eventRate, 0.0);
    }

    // ---- 5.185: the opt-in well-control reservoirs --------------------------
    // One linear PI inflow per node. Every mistake is loud: the node must
    // exist, be free with a positive compliance (an incompressible port has
    // nowhere to put the influx), carry no other rate account, and the law's
    // own numbers must be physical.
    for (const TransientReservoirConfig& rc : c.reservoirs)
    {
        const label i = net_.nodeIndex(rc.nodeName);
        NodeDyn& nd = nodes_[static_cast<std::size_t>(i)];
        if (nd.fixed)
            cfFatal("transient: reservoir '" + rc.name + "' feeds node '"
                    + rc.nodeName + "' which is a fixedPressure port: an"
                    " inflow boundary needs a free node to build pressure in");
        if (!(nd.bulk > 0))
            cfFatal("transient: reservoir '" + rc.name + "' feeds node '"
                    + rc.nodeName + "' with no bulkModulus: the influx has no"
                    " compliance to compress");
        for (const SourceEntry& e : sourceEntries_)
            if (e.nodeIndex == static_cast<std::size_t>(i))
                cfFatal("transient: node '" + rc.nodeName + "' is priced by"
                        " a sources{}/sourceSchedule entry and by reservoir '"
                        + rc.name + "': one rate account per node");
        if (nd.resIdx >= 0)
            cfFatal("transient: node '" + rc.nodeName + "' is fed by two"
                    " reservoirs: one inflow boundary per node");
        if (!(rc.pressure > 0))
            cfFatal("transient: reservoir '" + rc.name + "' needs a positive"
                    " pressure");
        if (!(rc.productivityIndex > 0))
            cfFatal("transient: reservoir '" + rc.name + "' needs a positive"
                    " productivityIndex");
        if (!(rc.gasFraction >= 0) || !(rc.gasFraction < 1))
            cfFatal("transient: reservoir '" + rc.name + "' gasFraction must"
                    " be in [0,1)");
        if (!(rc.depth > 0))
            cfFatal("transient: reservoir '" + rc.name + "' needs a positive"
                    " depth: the kill-mud weight is p_res / (g depth)");
        ReservoirEntry e;
        e.name = rc.name;
        e.nodeIndex = static_cast<std::size_t>(i);
        e.pressure = rc.pressure;
        e.J = rc.productivityIndex;
        e.gasFraction = rc.gasFraction;
        e.depth = rc.depth;
        e.active = nd.pPrev < rc.pressure;
        nd.resIdx = static_cast<int>(reservoirs_.size());
        reservoirs_.push_back(e);
    }

    for (const TransientPipeConfig& pc : c.pipes)
    {
        NetworkPipe geo;        geo.name = pc.name;
        geo.fromName = pc.fromName;
        geo.toName = pc.toName;
        geo.length = pc.length;
        geo.diameter = pc.diameter;
        geo.roughness = pc.roughness;
        geo.elevation = pc.elevation;
        geo.valveK = pc.valveK;
        geo.coupledResistance = pc.coupledResistance;
        geo.flowIndex = pc.flowIndex;
        geo.consistency = pc.consistency;
        geo.yieldStress = pc.yieldStress;
        if (pc.flowIndex < 0)
            cfFatal("transient: pipe '" + pc.name + "' has a negative"
                    " power-law flow index (0 means Newtonian)");
        if (pc.flowIndex > 0 && !(pc.consistency > 0))
            cfFatal("transient: pipe '" + pc.name + "' declares a flow index"
                    " without a positive consistency: at n != 1 the fluid's mu"
                    " is not the same quantity as K (Pa.s^n)");
        // A coupled edge may declare a rheology: it is the constitutive law of
        // the fluid inside its partner region, not a wall law this row reads.
        // It is forwarded below, after every check, and the wall-law fields are
        // cleared so holdsStatic/frictionDp/the crossing latch cannot see it.
        if (pc.yieldStress < 0)
            cfFatal("transient: pipe '" + pc.name + "' has a negative yield"
                    " stress: the threshold is the stress a static wall holds,"
                    " and a negative one would make the held set empty");
        // 5.74 (A2): gel and creep are rest-branch physics of a yielded wall.
        if (pc.gelStrength < 0)
            cfFatal("transient: pipe '" + pc.name + "' has a negative"
                    " gelStrength: gelation adds to the static wall stress");
        if (pc.gelStrength > 0 && pc.gelTime <= 0)
            cfFatal("transient: pipe '" + pc.name + "' declares gelStrength"
                    " without a positive gelTime: the rebuild law is"
                    " tau_g (1 - exp(-t_rest/t_g))");
        if (pc.gelStrength > 0 && !(pc.yieldStress > 0))
            cfFatal("transient: pipe '" + pc.name + "' declares gelStrength"
                    " without a yieldStress: there is no static wall to"
                    " strengthen");
        if (pc.creepRate < 0)
            cfFatal("transient: pipe '" + pc.name + "' has a negative"
                    " creepRate: the bi-viscosity regularisation needs a"
                    " positive gamma_dot_c");
        if (pc.creepRate > 0 && !(pc.yieldStress > 0))
            cfFatal("transient: pipe '" + pc.name + "' declares creepRate"
                    " without a yieldStress: without a threshold there is"
                    " no held branch to regularise");
        if (pc.coupledResistance < 0)
            cfFatal("transient: pipe '" + pc.name + "' has a negative coupled"
                    " resistance");
        if (pc.coupledResistance > 0)
            requireNoValve(pc.name, pc.valveK, pc.valveCloseStart,
                           pc.valveOpenStart, pc.valveEventsOn);
        // 5.70 / A3: the momentum-flux factors. beta is an area average of
        // v^2 / vbar^2, so it cannot sit below 1 (Cauchy-Schwarz); and on a
        // coupled pipe the account belongs to the interface -- it measures the
        // imbalance across the region's own stations and installs the
        // coefficient (5.43), so a declaration here would be silently
        // overwritten.
        if (pc.betaIn < 1.0 || pc.betaOut < 1.0)
            cfFatal("transient: pipe '" + pc.name + "' declares a momentum"
                    " flux factor below 1: beta = <v^2>/<v>^2 is an area"
                    " average, and no profile is fuller than flat");
        if (pc.coupledResistance > 0
            && (pc.betaIn != 1.0 || pc.betaOut != 1.0))
            cfFatal("transient: pipe '" + pc.name + "' is a coupled pipe"
                    " declaring betaIn/betaOut: the momentum-flux account of"
                    " a partner region belongs to the interface, which"
                    " measures it and installs the coefficient (5.43)");
        // The 1D momentum row of a coupled pipe charges `R*q` and nothing
        // else. Move the declaration onto `partnerRheology` and clear the
        // wall-law fields: otherwise `holdsStatic` would latch a Q = 0 row on
        // a threshold the partner region never claimed, and `frictionDp`
        // would price the rate a second time (§5.43).
        if (geo.coupledResistance > 0)
        {
            geo.partnerRheology = {geo.flowIndex, geo.consistency,
                                   geo.yieldStress};
            geo.flowIndex = 0;
            geo.consistency = 0;
            geo.yieldStress = 0;
        }
        net_.addPipe(geo);
        PipeDyn d;
        d.geo = net_.pipes.back();
        d.qPrev = pc.initialFlow;
        d.inertia = fluid_.rho * d.geo.length / pipeArea(d.geo.diameter);
        // The declared momentum-flux account: rho*(betaOut - betaIn)/A^2, the
        // same shape of term the coupled interface installs. Zero when the
        // case declares nothing, so an untouched case is bit-identical.
        d.convective = fluid_.rho * (pc.betaOut - pc.betaIn)
                     / (pipeArea(d.geo.diameter) * pipeArea(d.geo.diameter));
        d.betaIn = pc.betaIn;
        d.betaOut = pc.betaOut;
        d.closeStart = pc.valveCloseStart;
        d.closeEnd = pc.valveCloseEnd;
        d.kClosed = pc.valveKClosed;
        d.openStart = pc.valveOpenStart;
        d.openEnd = pc.valveOpenEnd;
        // The opt-in multi-segment table rides on the dynamic state like the
        // ramps do; declared or not, the step reads exactly one of the two.
        d.valveEventsOn = pc.valveEventsOn;
        d.eventT = pc.valveEventT;
        d.eventK = pc.valveEventK;
        // 5.74 (A2): the gel/creep declarations ride on the dynamic state.
        d.gelStrength = pc.gelStrength;
        d.gelTime = pc.gelTime;
        d.creepRate = pc.creepRate;
        d.creeps = pc.creepRate > 0;
        pipes_.push_back(d);
        pipeNames_.push_back(pc.name);
    }
    if (c.pipes.empty()) cfFatal("transient: case has no pipes");

    // E7a: opt-in wellbore equipment, resolved onto pipes by name. The loader
    // has validated the declarations; a programmatically built case still gets
    // a loud refusal for an unknown or doubly-claimed target.
    well_.assign(pipes_.size(), WellboreDyn{});
    std::vector<int> claimed(pipes_.size(), 0); // 1 = base valve, 2 = element
    for (std::size_t k = 0; k < pipes_.size(); ++k)
        if (pipes_[k].geo.valveK > 0 || pipes_[k].closeStart >= 0
            || pipes_[k].openStart >= 0 || pipes_[k].valveEventsOn)
            claimed[k] = 1;
    auto claim = [&](const std::string& pipe, const std::string& what) -> std::size_t
    {
        std::size_t idx = pipes_.size();
        for (std::size_t k = 0; k < pipeNames_.size(); ++k)
            if (pipeNames_[k] == pipe) { idx = k; break; }
        if (idx == pipes_.size())
            cfFatal("transient: " + what + " targets unknown pipe '" + pipe + "'");
        if (claimed[idx] == 1)
            cfFatal("transient: pipe '" + pipe + "' declares a base valve and a"
                    " wellbore element: one local-loss account per pipe");
        if (claimed[idx] == 2)
            cfFatal("transient: pipe '" + pipe + "' declares two wellbore"
                    " elements: one per pipe");
        claimed[idx] = 2;
        return idx;
    };
    for (const TransientSssvConfig& s : c.sssvs)
    {
        const std::size_t k = claim(s.pipeName, "sssv '" + s.name + "'");
        WellboreDyn& w = well_[k];
        w.kind = 1;
        w.openK = s.openK; w.closedK = s.closedK;
        w.closeStart = s.closeStart; w.closeEnd = s.closeEnd;
        w.reopenStart = s.reopenStart; w.reopenEnd = s.reopenEnd;
        w.lastK = sssvKAt(w.openK, w.closedK, w.closeStart, w.closeEnd,
                          w.reopenStart, w.reopenEnd, 0.0);
        w.lastSealed = sssvSealed(w.closeStart, w.closeEnd, w.reopenStart, 0.0);
    }
    for (const TransientChokeConfig& ch : c.chokes)
    {
        const std::size_t k = claim(ch.pipeName, "choke '" + ch.name + "'");
        WellboreDyn& w = well_[k];
        w.kind = 2;
        if (ch.diameter > 0)
        {
            const double a = pipeArea(pipes_[k].geo.diameter);
            const double ao = pipeArea(ch.diameter);
            const double ratio = a / (ch.dischargeCoeff * ao);
            w.extraK = ratio * ratio; // orifice: K = (A / (Cd Ao))^2
        }
        else
        {
            w.extraK = ch.k;
        }
        w.lastK = w.extraK;
    }
    for (const TransientEspConfig& e : c.esps)
    {
        const std::size_t k = claim(e.pipeName, "esp '" + e.name + "'");
        WellboreDyn& w = well_[k];
        w.kind = 3;
        w.stages = e.stages; w.h0 = e.shutoffHead; w.hq = e.headAtRate;
        w.qr = e.rate;
        w.speed0 = e.speed0; w.speedCmd = e.speedCmd;
        w.speedCmdTime = e.speedCmdTime; w.speedTau = e.speedTau;
        w.lastSpeed = espSpeedAt(w.speed0, w.speedCmd, w.speedCmdTime,
                                 w.speedTau, 0.0);
        w.lastHead = espHeadAt(w.stages, w.h0, w.hq, w.qr,
                               pipes_[k].qPrev, w.lastSpeed);
    }

    // unknown map: free nodes get columns after the pipe-flow block
    p_.assign(nodes_.size(), 0.0);
    std::size_t nFree = 0;
    for (std::size_t i = 0; i < nodes_.size(); ++i)
    {
        if (nodes_[i].fixed) p_[i] = nodes_[i].pPrev;
        else { nodes_[i].unk = static_cast<double>(nFree++); p_[i] = nodes_[i].pPrev; }
    }
    q_.assign(pipes_.size(), 0.0);
    for (std::size_t k = 0; k < pipes_.size(); ++k) q_[k] = pipes_[k].qPrev;
}

std::vector<ReservoirReadout> TransientNetwork::reservoirReadouts() const
{
    std::vector<ReservoirReadout> out;
    out.reserve(reservoirs_.size());
    for (const ReservoirEntry& r : reservoirs_)
    {
        ReservoirReadout o;
        o.name = r.name;
        o.nodeName = net_.nodeNames[r.nodeIndex];
        o.pressure = r.pressure;
        o.nodePressure = p_[r.nodeIndex];
        o.influx = std::max(r.J * (r.pressure - p_[r.nodeIndex]), 0.0);
        o.gasInflux = o.influx * r.gasFraction;
        o.volume = r.volume;
        o.gasVolume = r.gasVolume;
        o.underbalance = r.pressure - p_[r.nodeIndex];
        o.underbalanced = p_[r.nodeIndex] < r.pressure;
        o.killMudWeight = r.pressure / (kGravity * r.depth);
        o.sidpp = r.pressure - fluid_.rho * kGravity * r.depth;
        out.push_back(o);
    }
    return out;
}

const std::vector<std::string>& TransientNetwork::pipeNames() const
{
    return pipeNames_;
}

double TransientNetwork::solveStep(double dt)
{
    if (!(dt > 0)) cfFatal("transient: dt must be positive");
    const std::size_t nP = pipes_.size();
    const std::size_t nN = nodes_.size();
    std::size_t nFree = 0;
    for (const NodeDyn& d : nodes_)
        nFree = std::max(nFree, static_cast<std::size_t>(d.unk + 1));
    const std::size_t nEq = nP + nFree;

    std::vector<double> x(nEq);
    for (std::size_t k = 0; k < nP; ++k) x[k] = q_[k];
    for (std::size_t i = 0; i < nN; ++i)
        if (nodes_[i].unk >= 0)
            x[nP + static_cast<std::size_t>(nodes_[i].unk)] = p_[i];

    auto pressureOf = [&](label i) -> double
    {
        const NodeDyn& d = nodes_[static_cast<std::size_t>(i)];
        return d.fixed ? d.pPrev : x[nP + static_cast<std::size_t>(d.unk)];
    };

    std::vector<double> r(nEq);
    std::vector<std::vector<double>> J(nEq, std::vector<double>(nEq, 0.0));

    // valve K is evaluated at the implicit new time level t + dt and held
    // constant during the Newton solve. One of two schedules answers: the
    // opt-in multi-segment table when the case declared one, the
    // single-cycle ramps otherwise (the loader refuses declaring both).
    const double tNew = t_ + dt;
    std::vector<double> kv(nP);
    // E7a: the three elements are evaluated at the same implicit new time as
    // the valve. SSSV and choke are local loss coefficients and fold into kv;
    // the SSSV's seal is a Q = 0 latch; the ESP speed is a time-only state and
    // its head is quadratic in Q, so the head is recomputed from the iterate
    // inside the Newton loop and only its slope enters the diagonal.
    std::vector<char> sealed(nP, 0);
    std::vector<double> espN(nP, 0.0);
    for (std::size_t k = 0; k < nP; ++k)
    {
        WellboreDyn& w = well_[k];
        kv[k] = pipes_[k].valveEventsOn
                    ? eventTableAt(pipes_[k].eventT, pipes_[k].eventK, tNew)
                    : valveKAt(pipes_[k].geo.valveK, pipes_[k].kClosed,
                               pipes_[k].closeStart, pipes_[k].closeEnd,
                               pipes_[k].openStart, pipes_[k].openEnd, tNew);
        if (w.kind == 1)
        {
            kv[k] += sssvKAt(w.openK, w.closedK, w.closeStart, w.closeEnd,
                             w.reopenStart, w.reopenEnd, tNew);
            sealed[k] = sssvSealed(w.closeStart, w.closeEnd, w.reopenStart, tNew)
                            ? 1 : 0;
        }
        else if (w.kind == 2)
        {
            kv[k] += w.extraK;
        }
        else if (w.kind == 3)
        {
            espN[k] = espSpeedAt(w.speed0, w.speedCmd, w.speedCmdTime,
                                 w.speedTau, tNew);
        }
    }

    // source rate at the same implicit new time level; reset-then-sum so a
    // node with multiple scheduled sources sees their sum (not a carry-over).
    // Each entry reads its own kind of schedule: the multi-segment table or
    // the single-cycle ramp, never both (the loader refuses the overlap).
    for (const SourceEntry& e : sourceEntries_)
        nodes_[e.nodeIndex].source = 0.0;
    for (const SourceEntry& e : sourceEntries_)
        nodes_[e.nodeIndex].source +=
            e.eventsOn ? eventTableAt(e.eventT, e.eventRate, tNew)
                       : sourceRateAt(e.rate, e.openRate, e.start, e.end,
                                      e.reopenStart, e.reopenEnd, tNew);

    double maxRel = 0;

    // 5.185: the one-sided inflow law is latched at step entry, the same idiom
    // the held-wall branch uses. A shut-in well crossing p_res must not have
    // the law flip under the Newton iterate, and a latched branch keeps the
    // node row linear in the pressure.
    for (ReservoirEntry& e : reservoirs_)
        e.active = nodes_[e.nodeIndex].pPrev < e.pressure;

    // 5.49: which yielded pipes the static wall is holding, decided once at the
    // start of the step from the state the step starts in. Latching it here
    // rather than asking the Newton loop every pass is what keeps the switch
    // from oscillating: the drive a held pipe sees depends on the node
    // pressures the same solve is moving, so an unratcheted test could flip the
    // row between `Q = 0` and the moving-wall law inside one time step.
    for (std::size_t k = 0; k < nP; ++k)
    {
        PipeDyn& pd = pipes_[k];
        const double drive = drivePa(pd.geo, fluid_, pressureOf(pd.geo.from),
                                      pressureOf(pd.geo.to));
        pd.driveLatch = drive;
        // 5.74: a pipe on the creep branch stays on it (its qPrev is not 0,
        // but it is at rest) -- the rest test takes the latched branch as the
        // rest witness.
        pd.held = holdsStatic(pd.geo, pd.qPrev, drive, pd.tRest,
                              pd.gelStrength, pd.gelTime, pd.held);
    }

    for (int newton = 0; newton < gates::transientNewtonMaxIter; ++newton)
    {
        for (double& v : r) v = 0;
        for (auto& row : J) std::fill(row.begin(), row.end(), 0.0);
        maxRel = 0;

        for (std::size_t k = 0; k < nP; ++k)
        {
            PipeDyn& pd = pipes_[k];
            // E7a: a sealed SSSV overrides every wall branch -- its row is
            // Q = 0 with no loss law read, exactly like a held pipe but with
            // no creep and no release until the declared reopen time.
            if (sealed[k])
            {
                x[k] = 0;
                r[k] = 0;
                J[k][k] = 1.0;
                continue;
            }
            const double pF = pressureOf(pd.geo.from);
            const double pT = pressureOf(pd.geo.to);
            const double dp = pF - pT;
            const double grav = fluid_.rho * kGravity * pd.geo.elevation;
            // A coupled pipe's friction is its partner's statement, not ours:
            // the detailed model already reports the whole drop of that region,
            // so adding this pipe's own Churchill friction on top would charge
            // the network twice for the same wall. It stays linear in Q because
            // the interface refreshes the coefficient every pass. Gravity and
            // the region's inertance stay on the 1D side - the partner is asked
            // for a drop at an instant, not for the fluid's momentum; `inertia`
            // is the pipe's own rho*L/A unless the interface measured the
            // region's and handed the number back (setCoupledInertia).
            const double R = pd.geo.coupledResistance;
            // 5.49: the static branch, in its two shapes. Either the wall was
            // already holding this pipe at the start of the step (`pd.held`,
            // latched above), or the drive on offer now is below the threshold
            // and the iterate has carried the rate through zero -- which the
            // moving-wall law alone answers with a reverse flow the same size,
            // because `dp_wall` is odd and never drops below the threshold.
            // Both are stated the same way: the row becomes `Q = 0`. Its own
            // columns into the node pressures go away with it -- a held pipe
            // does not answer to the drive, it refuses it -- while the node rows
            // keep their entries in this pipe's column, so continuity still sees
            // the forced zero and the volume a held pipe occupies stays sealed.
            const double drive = dp - grav;
            const bool crossing = pd.geo.yieldStress > 0 && pd.qPrev != 0
                               && std::abs(drive) <= yieldThresholdPa(pd.geo)
                               && x[k] * pd.qPrev <= 0.0;
            if (pd.held || crossing)
            {
                // 5.74: the held row is Q = 0 unless the case declared a
                // sub-threshold creep rate -- then the bi-viscosity branch
                // answers with a very stiff Newtonian row,
                //     q = pi d^4 drive / (128 mu_eff L),
                // mu_eff = mu + tau0/gamma_dot_c, built on the LATCHED drive
                // (a forced value like Q = 0, not a live law). The static
                // wall supplies the threshold minus |drive|; continuity sees
                // the creep rate instead of a sealed zero. A crossing pipe
                // parks at Q = 0 even with creep declared: it stops, and
                // joins the creep branch from the next step's latch.
                if (pd.held && pd.creeps && !crossing)
                {
                    const double muEff = fluid_.mu
                                       + pd.geo.yieldStress / pd.creepRate;
                    x[k] = M_PI * std::pow(pd.geo.diameter, 4)
                         * pd.driveLatch / (128.0 * muEff * pd.geo.length);
                }
                else
                {
                    x[k] = 0;
                }
                r[k] = 0;
                J[k][k] = 1.0;
                continue;
            }
            const double fr = R > 0
                                  ? R * x[k]
                                  : frictionDp(pd.geo, fluid_, x[k], drive);
            const double slope = R > 0
                                    ? 1.0 / R
                                    : frictionSlope(pd.geo, fluid_, x[k], drive);
            const double ar = pipeArea(pd.geo.diameter);
            const double vd = kv[k] * fluid_.rho * std::abs(x[k]) * x[k]
                            / (2.0 * ar * ar);
            // The same shape of term as the valve's, and the same kind of
            // statement: a partner region whose profile is fuller at one
            // station than at the other owes its momentum budget
            // `rho*(betaOut - betaIn)*|Q|*Q/A^2`, which is a momentum-flux
            // imbalance rather than a wall force, so it scales with Q^2 where
            // the resistance the interface installs scales with Q. Charging it
            // here is what lets the chord stop carrying it; zero leaves the row
            // exactly the equation it was before.
            const double cc = pd.convective * std::abs(x[k]) * x[k];
            // E7a: a declared ESP adds its head as a pressure rise along the
            // flow direction and the curve's (negative) slope to the diagonal,
            // so the implicit row stays well-posed at any q. Both are zero for
            // an undeclared pipe, bit-for-bit the pre-E7a row.
            const WellboreDyn& w = well_[k];
            const double headPa = w.kind == 3
                ? fluid_.rho * kGravity
                  * espHeadAt(w.stages, w.h0, w.hq, w.qr, x[k], espN[k])
                : 0.0;
            const double headSlope = w.kind == 3
                ? -fluid_.rho * kGravity
                  * espHeadSlope(w.stages, w.h0, w.hq, w.qr, x[k])
                : 0.0;
            r[k] = pd.inertia / dt * (x[k] - pd.qPrev) - dp + grav + fr + vd + cc
                 - headPa;
            J[k][k] = pd.inertia / dt + slope
                    + kv[k] * fluid_.rho * std::abs(x[k]) / (ar * ar)
                    + 2.0 * pd.convective * std::abs(x[k]) + headSlope;
            const NodeDyn& nf = nodes_[static_cast<std::size_t>(pd.geo.from)];
            const NodeDyn& nb = nodes_[static_cast<std::size_t>(pd.geo.to)];
            if (!nf.fixed) J[k][nP + static_cast<std::size_t>(nf.unk)] -= 1.0;
            if (!nb.fixed) J[k][nP + static_cast<std::size_t>(nb.unk)] += 1.0;
            const double scale = pd.inertia / dt * std::max(std::abs(x[k]), std::abs(pd.qPrev))
                               + std::abs(dp) + std::abs(fr) + std::abs(vd)
                               + std::abs(cc) + std::abs(headPa) + 1.0;
            maxRel = std::max(maxRel, std::abs(r[k]) / scale);
        }

        for (std::size_t i = 0; i < nN; ++i)
        {
            const NodeDyn& nd = nodes_[i];
            if (nd.fixed) continue;
            const std::size_t e = nP + static_cast<std::size_t>(nd.unk);
            double sumIn = nd.source;
            // 5.185: the reservoir's PI inflow, implicit in the node pressure
            // (q = J (p_res - p) on the latched active branch; zero when the
            // node is overbalanced at the state the step started in).
            const int ri = nd.resIdx;
            const bool resOn = ri >= 0 && reservoirs_[static_cast<std::size_t>(ri)].active;
            if (resOn)
            {
                const ReservoirEntry& r = reservoirs_[static_cast<std::size_t>(ri)];
                sumIn += r.J * (r.pressure - x[e]);
            }
            for (std::size_t k = 0; k < nP; ++k)
            {
                if (pipes_[k].geo.to == static_cast<label>(i)) sumIn += x[k];
                if (pipes_[k].geo.from == static_cast<label>(i)) sumIn -= x[k];
            }
            const double stor = nd.bulk > 0 ? nd.volume / (nd.bulk * dt) * (x[e] - nd.pPrev) : 0.0;
            r[e] = stor - sumIn;
            if (nd.bulk > 0) J[e][e] = nd.volume / (nd.bulk * dt);
            if (resOn)
                J[e][e] += reservoirs_[static_cast<std::size_t>(ri)].J;
            for (std::size_t k = 0; k < nP; ++k)
            {
                if (pipes_[k].geo.to == static_cast<label>(i)) J[e][k] -= 1.0;
                if (pipes_[k].geo.from == static_cast<label>(i)) J[e][k] += 1.0;
            }
            const double scale = std::abs(stor) + std::abs(sumIn) + 1.0;
            maxRel = std::max(maxRel, std::abs(r[e]) / scale);
        }

        if (maxRel < gates::transientNewtonGate) break;
        std::vector<double> negR(nEq);
        for (std::size_t i = 0; i < nEq; ++i) negR[i] = -r[i];
        detail::solveDenseLinear(J, negR);
        for (std::size_t i = 0; i < nEq; ++i) x[i] += negR[i];
    }

    // E7a: freeze the values the solved step actually used, for the report
    // lines the checker re-derives from the case text independently.
    for (std::size_t k = 0; k < nP; ++k)
    {
        WellboreDyn& w = well_[k];
        if (w.kind == 1)
        {
            w.lastK = sssvKAt(w.openK, w.closedK, w.closeStart, w.closeEnd,
                              w.reopenStart, w.reopenEnd, tNew);
            w.lastSealed = sealed[k] != 0;
        }
        else if (w.kind == 2)
        {
            w.lastK = w.extraK;
        }
        else if (w.kind == 3)
        {
            w.lastSpeed = espN[k];
            w.lastHead = espHeadAt(w.stages, w.h0, w.hq, w.qr, x[k], w.lastSpeed);
        }
    }

    xPending_ = std::move(x);
    pendingDt_ = dt;
    pending_ = true;
    return maxRel;
}

void TransientNetwork::commitStep()
{
    if (!pending_) cfFatal("transient: commitStep with no step solved");
    const double dt = pendingDt_;
    const std::size_t nP = pipes_.size();
    const std::size_t nN = nodes_.size();
    const std::vector<double>& x = xPending_;

    for (std::size_t k = 0; k < nP; ++k)
    {
        PipeDyn& pd = pipes_[k];
        pd.qPrev = x[k];
        q_[k] = x[k];
        // 5.74: the gel clock runs only at rest -- a step spent held (or
        // creeping, which is the held branch) accumulates, a moving step
        // clears it. Crossing pipes were moving: their clock is already 0.
        if (pd.geo.yieldStress > 0)
            pd.tRest = pd.held ? pd.tRest + dt : 0.0;
    }
    for (std::size_t i = 0; i < nN; ++i)
    {
        NodeDyn& nd = nodes_[i];
        if (nd.fixed) { p_[i] = nd.pPrev; continue; }
        nd.pPrev = x[nP + static_cast<std::size_t>(nd.unk)];
        p_[i] = nd.pPrev;
    }

    // ---- drift-flux transport of mobile liquid (explicit, upwind) ----
    std::vector<double> oldL(nN), alpha(nN);
    for (std::size_t i = 0; i < nN; ++i)
    {
        oldL[i] = nodes_[i].vL;
        alpha[i] = nodes_[i].volume > 0 ? nodes_[i].vL / nodes_[i].volume : 0.0;
    }
    for (std::size_t k = 0; k < nP; ++k)
    {
        const NetworkPipe& geo = pipes_[k].geo;
        const double q = pipes_[k].qPrev;
        const double aFrom = alpha[static_cast<std::size_t>(geo.from)];
        const double aTo = alpha[static_cast<std::size_t>(geo.to)];
        if (q >= 0)
        {
            oldL[static_cast<std::size_t>(geo.to)] += dt * q * aFrom;
            oldL[static_cast<std::size_t>(geo.from)] -= dt * q * aFrom;
        }
        else
        {
            oldL[static_cast<std::size_t>(geo.from)] += dt * (-q) * aTo;
            oldL[static_cast<std::size_t>(geo.to)] -= dt * (-q) * aTo;
        }
    }
    for (std::size_t i = 0; i < nN; ++i)
    {
        NodeDyn& nd = nodes_[i];
        if (nd.fixed) continue; // Dirichlet reservoir: infinite inventory
        nd.vL = oldL[i] + dt * nd.source * nd.sourceFraction;
        // 5.185: the reservoir's liquid share enters the node's mobile
        // inventory; the rest of the influx is gas and fills the remaining
        // node volume through the compliance (the node engine books volume).
        if (nd.resIdx >= 0)
        {
            ReservoirEntry& r = reservoirs_[static_cast<std::size_t>(nd.resIdx)];
            if (r.active)
            {
                const double qIn = std::max(r.J * (r.pressure - nd.pPrev), 0.0);
                nd.vL += dt * qIn * (1.0 - r.gasFraction);
                r.volume += dt * qIn;
                r.gasVolume += dt * qIn * r.gasFraction;
            }
        }
        if (nd.vL < -1e-9 * std::max(1.0, nd.volume))
            cfFatal("transient: node '" + net_.nodeNames[i]
                    + "' runs dry - reduce dt");
        if (nd.vL < 0) nd.vL = 0;
        if (nd.vL + nd.vH > nd.volume + 1e-9 * (nd.volume + 1.0))
            cfFatal("transient: node '" + net_.nodeNames[i]
                    + "' liquid exceeds its volume - reduce dt");
        hydrateKineticsStep(nd.vL, nd.vH, nd.rate, dt);
    }

    t_ += dt;
    pending_ = false;
}

double TransientNetwork::step(double dt)
{
    const double maxRel = solveStep(dt);
    commitStep();
    return maxRel;
}

double TransientNetwork::cflStep(double target, double lo, double hi) const
{
    double fastest = 0.0; // min_k(L A / |q|): the tightest Courant budget
    for (std::size_t k = 0; k < pipes_.size(); ++k)
    {
        const double q = std::abs(q_[k]);
        if (q <= 0) continue; // a pipe at rest imposes nothing
        const NetworkPipe& geo = pipes_[k].geo;
        const double budget = geo.length * (M_PI * geo.diameter * geo.diameter
                                            / 4.0) / q;
        fastest = fastest > 0 ? std::min(fastest, budget) : budget;
    }
    const double raw = fastest > 0 ? target * fastest : hi;
    return std::min(std::max(raw, lo), hi);
}

double TransientNetwork::pendingPipeFlow(label k) const
{
    if (!pending_) cfFatal("transient: no step is pending to read a flow from");
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return xPending_[static_cast<std::size_t>(k)];
}

double TransientNetwork::pendingPressure(label i) const
{
    if (!pending_) cfFatal("transient: no step is pending to read a pressure from");
    if (i < 0 || static_cast<std::size_t>(i) >= nodes_.size())
        cfFatal("transient: no node at index " + std::to_string(i));
    const NodeDyn& d = nodes_[static_cast<std::size_t>(i)];
    return d.fixed
               ? d.pPrev
               : xPending_[pipes_.size() + static_cast<std::size_t>(d.unk)];
}

bool TransientNetwork::pipeCoupled(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return pipes_[static_cast<std::size_t>(k)].geo.coupledResistance > 0;
}

double TransientNetwork::coupledResistanceOf(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return pipes_[static_cast<std::size_t>(k)].geo.coupledResistance;
}

double TransientNetwork::pipeInertia(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return pipes_[static_cast<std::size_t>(k)].inertia;
}

double TransientNetwork::pipeYieldThreshold(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return yieldThresholdPa(pipes_[static_cast<std::size_t>(k)].geo);
}

bool TransientNetwork::pipeHeld(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return pipes_[static_cast<std::size_t>(k)].held;
}

int TransientNetwork::pipeWellboreKind(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= well_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return well_[static_cast<std::size_t>(k)].kind;
}

double TransientNetwork::pipeWellboreK(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= well_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return well_[static_cast<std::size_t>(k)].lastK;
}

bool TransientNetwork::pipeWellboreSealed(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= well_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return well_[static_cast<std::size_t>(k)].lastSealed;
}

double TransientNetwork::pipeEspSpeed(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= well_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return well_[static_cast<std::size_t>(k)].lastSpeed;
}

double TransientNetwork::pipeEspHead(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= well_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return well_[static_cast<std::size_t>(k)].lastHead;
}

double TransientNetwork::pipeFlowIndex(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return pipes_[static_cast<std::size_t>(k)].geo.flowIndex;
}

double TransientNetwork::pipeConsistencyOf(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return pipes_[static_cast<std::size_t>(k)].geo.consistency;
}

Rheology TransientNetwork::partnerRheologyOf(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return pipes_[static_cast<std::size_t>(k)].geo.partnerRheology;
}

void TransientNetwork::setCoupledInertia(label k, double inertia)
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    // Checked before anything else, so a caller cannot be told a pipe has no
    // partner when what it handed over was not an inertance.
    if (!(std::isfinite(inertia) && inertia > 0))
        cfFatal("transient: a partner region's inertance must be a positive,"
                " finite rho*Volume/Area^2; got " + std::to_string(inertia)
                + " (kg/m4)");
    PipeDyn& pd = pipes_[static_cast<std::size_t>(k)];
    if (pd.geo.coupledResistance <= 0)
        cfFatal("transient: pipe '" + pipeNames_[static_cast<std::size_t>(k)]
                + "' has no partner model, so it has no measured inertance to"
                  " take");
    pd.inertia = inertia;
    pending_ = false;
}

void TransientNetwork::setCoupledConvective(label k, double coefficient)
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    // Signed on purpose: a region whose outlet section is fuller than its inlet
    // one bills the row, a region that sheds momentum flux pays it. So the only
    // thing that cannot be a coefficient is a number that is not a number.
    if (!std::isfinite(coefficient))
        cfFatal("transient: a partner region's convective coefficient must be a"
                " finite rho*(betaOut-betaIn)/A^2; got "
                + std::to_string(coefficient) + " (kg/m6)");
    PipeDyn& pd = pipes_[static_cast<std::size_t>(k)];
    if (pd.geo.coupledResistance <= 0)
        cfFatal("transient: pipe '" + pipeNames_[static_cast<std::size_t>(k)]
                + "' has no partner model, so it has no measured profile change to"
                  " charge its row with");
    pd.convective = coefficient;
    pending_ = false;
}

double TransientNetwork::coupledConvectiveOf(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return pipes_[static_cast<std::size_t>(k)].convective;
}

double TransientNetwork::pipeConvectiveOf(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return pipes_[static_cast<std::size_t>(k)].convective;
}

double TransientNetwork::pipeBetaInOf(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return pipes_[static_cast<std::size_t>(k)].betaIn;
}

double TransientNetwork::pipeBetaOutOf(label k) const
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    return pipes_[static_cast<std::size_t>(k)].betaOut;
}

void TransientNetwork::setCoupledResistance(label k, double resistance)
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    if (!(resistance > 0))
        cfFatal("transient: coupled pipe '" + pipeNames_[k] + "' needs a"
                " positive resistance, got " + std::to_string(resistance));
    const PipeDyn& pd = pipes_[static_cast<std::size_t>(k)];
    requireNoValve(pipeNames_[k], pd.geo.valveK, pd.closeStart, pd.openStart,
                   pd.valveEventsOn);
    pipes_[static_cast<std::size_t>(k)].geo.coupledResistance = resistance;
    // whatever was pending was solved with the old coefficient, so it is no
    // longer a state of the system and has to be re-solved
    pending_ = false;
}

std::vector<double> TransientNetwork::liquidVolumes() const
{
    std::vector<double> v(nodes_.size());
    for (std::size_t i = 0; i < nodes_.size(); ++i) v[i] = nodes_[i].vL;
    return v;
}

std::vector<double> TransientNetwork::hydrateVolumes() const
{
    std::vector<double> v(nodes_.size());
    for (std::size_t i = 0; i < nodes_.size(); ++i) v[i] = nodes_[i].vH;
    return v;
}

std::vector<double> TransientNetwork::liquidFractions() const
{
    std::vector<double> a(nodes_.size(), 0.0);
    for (std::size_t i = 0; i < nodes_.size(); ++i)
        a[i] = nodes_[i].volume > 0 ? nodes_[i].vL / nodes_[i].volume : 0.0;
    return a;
}

std::vector<double> TransientNetwork::hydrateFractions() const
{
    std::vector<double> h(nodes_.size(), 0.0);
    for (std::size_t i = 0; i < nodes_.size(); ++i)
    {
        const double tot = nodes_[i].vL + nodes_[i].vH;
        h[i] = tot > 0 ? nodes_[i].vH / tot : 0.0;
    }
    return h;
}

TransientSeries runTransient(const TransientCase& c, double dt,
                             double endTime, int reportEvery)
{
    TransientNetwork sim(c);
    TransientSeries out;
    out.nodeNames.assign(sim.nodeNames().begin(), sim.nodeNames().end());
    out.pipeNames.reserve(c.pipes.size());
    for (const TransientPipeConfig& p : c.pipes) out.pipeNames.push_back(p.name);

    auto sample = [&]()
    {
        TransientSample s;
        s.t = sim.time();
        s.pressure = sim.pressure();
        s.liquidFraction = sim.liquidFractions();
        s.hydrateFraction = sim.hydrateFractions();
        s.flowRate = sim.pipeFlows();
        out.samples.push_back(std::move(s));
    };

    sample();
    const int report = std::max(1, reportEvery);
    // 5.66 / C4: endTime no longer has to be a multiple of dt. The step count
    // is ceil(endTime/dt) with a tolerance so a floating-point near-multiple
    // (0.12/0.03 == 4.000000000000001) still counts as its integer, and the
    // LAST step is clamped to the remainder h = endTime - (n-1) dt, which the
    // implicit row reads through its own dt argument. tPrev is indexed, not
    // accumulated, so the clamp is exact. A case whose endTime IS a multiple
    // of dt gets h == dt on every step and the same sample set as before.
    if (!(dt > 0) || !(endTime > 0))
        cfFatal("transient: dt and endTime must both be positive");
    const int nSteps = static_cast<int>(std::ceil(endTime / dt - 1e-9));
    for (int i = 1; i <= nSteps; ++i)
    {
        const double tPrev = std::min((i - 1) * dt, endTime);
        const double h = std::min(dt, endTime - tPrev);
        sim.step(h);
        if (i % report == 0 || i == nSteps) sample();
    }
    return out;
}

TransientSeries runTransient(const TransientCase& c,
                             const TransientTimeControls& controls)
{
    // cfl absent (0) keeps the fixed-step march bit for bit.
    if (!(controls.cfl > 0))
        return runTransient(c, controls.dt, controls.endTime,
                            controls.reportEvery);

    // 5.68 / D5: CFL-adaptive stepping. The next dt is the Courant rule over
    // the whole network read off the flows the step just committed; the FIRST
    // step is the case's own dt. t accumulates -- with varying dt there is no
    // index to position a step -- and the last step is clamped to
    // endTime - t, with t snapped to endTime exactly when the clamp fires.
    if (!(controls.dt > 0) || !(controls.endTime > 0))
        cfFatal("transient: dt and endTime must both be positive");
    if (!(controls.dtMin > 0))
        cfFatal("transient: adaptive stepping needs dtMin > 0");
    if (controls.dtMax < controls.dtMin)
        cfFatal("transient: dtMax must not be below dtMin");

    TransientNetwork sim(c);
    TransientSeries out;
    out.nodeNames.assign(sim.nodeNames().begin(), sim.nodeNames().end());
    out.pipeNames.reserve(c.pipes.size());
    for (const TransientPipeConfig& p : c.pipes) out.pipeNames.push_back(p.name);

    auto sample = [&]()
    {
        TransientSample s;
        s.t = sim.time();
        s.pressure = sim.pressure();
        s.liquidFraction = sim.liquidFractions();
        s.hydrateFraction = sim.hydrateFractions();
        s.flowRate = sim.pipeFlows();
        out.samples.push_back(std::move(s));
    };

    sample();
    const int report = std::max(1, controls.reportEvery);
    double dtCur = controls.dt;
    double t = 0.0;
    long long nSteps = 0;
    while (t < controls.endTime - 1e-12)
    {
        const double rest = controls.endTime - t;
        const bool clampFires = rest <= dtCur;
        const double h = clampFires ? rest : dtCur;
        sim.step(h);
        ++nSteps;
        t = clampFires ? controls.endTime : t + h;
        if (nSteps % report == 0 || clampFires) sample();
        dtCur = sim.cflStep(controls.cfl, controls.dtMin, controls.dtMax);
    }
    return out;
}

TransientCase loadTransientFromDict(const Dictionary& caseDict)
{
    TransientCase c;
    if (!caseDict.has("fluid")) cfFatal("transient: case has no fluid{} block");
    const Dictionary& fd = caseDict.getDict("fluid");
    c.fluid.rho = fd.getScalar("rho");
    c.fluid.mu = fd.getScalar("mu");

    if (!caseDict.has("nodes")) cfFatal("transient: case has no nodes{} block");
    const Dictionary& nodes = caseDict.getDict("nodes");
    for (const std::string& name : nodes.keys())
    {
        const Dictionary& nd = nodes.getDict(name);
        TransientNodeConfig n;
        n.name = name;
        n.volume = nd.getScalarOr("volume", 0.0);
        n.bulkModulus = nd.getScalarOr("bulkModulus", 0.0);
        n.liquidFraction = nd.getScalarOr("liquidFraction", 1.0);
        n.hydrateFraction = nd.getScalarOr("hydrateFraction", 0.0);
        n.hydrateRate = nd.getScalarOr("hydrateRate", 0.0);
        if (nd.has("fixedPressure"))
        {
            n.pressureFixed = true;
            n.pressure = nd.getScalar("fixedPressure");
        }
        else
        {
            n.pressure = nd.getScalarOr("pressure", 0.0);
        }
        c.nodes.push_back(n);
    }

    if (!caseDict.has("pipes")) cfFatal("transient: case has no pipes{} block");
    const Dictionary& pipes = caseDict.getDict("pipes");
    for (const std::string& name : pipes.keys())
    {
        const Dictionary& pd = pipes.getDict(name);
        TransientPipeConfig p;
        p.name = name;
        p.fromName = pd.getWord("from");
        p.toName = pd.getWord("to");
        p.length = pd.getScalar("length");
        p.diameter = pd.getScalar("diameter");
        if (pd.has("dInner"))
            cfFatal("transient: pipe '" + name + "' declares dInner: the "
                    "annulus section (5.73) is a steady-network declaration; "
                    "the transient row does not take one yet");
        if (pd.has("ductHeight"))
            cfFatal("transient: pipe '" + name + "' declares ductHeight: the "
                    "rectangular duct (5.80) is a steady-network declaration; "
                    "the transient row does not take one yet");
        if (pd.has("eccentricity"))
            cfFatal("transient: pipe '" + name + "' declares eccentricity: the "
                    "eccentric annulus (5.87) is a steady-network "
                    "declaration; the transient row does not take one yet");
        if (pd.has("flowPattern"))
            cfFatal("transient: pipe '" + name + "' declares flowPattern: the "
                    "writable flow-pattern dispatch (5.88) is a steady-network "
                    "black-oil declaration; the transient row has no "
                    "Beggs-Brill pattern map to override");
        if (pd.has("correlation"))
            cfFatal("transient: pipe '" + name + "' declares correlation: the "
                    "alternative holdup laws (5.143) are a steady-network "
                    "black-oil declaration; the transient row has no "
                    "Beggs-Brill holdup block to replace");
        // 5.144 (G6): the steady thermal march is a marched-network pipe
        // declaration. The transient row advances its own energy state (and
        // the two-fluid row's thermal is a nested thermal{} block, not
        // these keys), so a pipe-level thermal key here would be read by
        // nothing - refuse it loudly.
        if (pd.has("jouleThomson") || pd.has("wallHeatTransfer")
            || pd.has("wallHeatTransferSegments") || pd.has("wallModel")
            || pd.has("wallTemperature") || pd.has("surfaceTemperature")
            || pd.has("geothermalGradient") || pd.has("heatCapacity"))
            cfFatal("transient: pipe '" + name + "' declares a thermal "
                    "parameter: the steady temperature march (5.144) is a "
                    "steady-network pipe declaration; the transient row "
                    "advances its own state and reads none of these keys");
        if (pd.has("liftDomePressure") || pd.has("liftAreaRatio")
            || pd.has("liftSpringPressure") || pd.has("liftDischargeCoeff")
            || pd.has("liftKappa"))
            cfFatal("transient: pipe '" + name + "' declares a lift valve: the "
                    "steady gas-lift valve (5.89) data is priced on node "
                    "pressures and is not a transient declaration; declare a "
                    "transient artificial-lift element in the top-level esp{} "
                    "block instead");
        if (pd.has("rodPlungerDiameter") || pd.has("rodStrokeLength")
            || pd.has("rodStrokesPerMinute") || pd.has("rodDiameter")
            || pd.has("rodLength") || pd.has("rodYoungModulus")
            || pd.has("rodDensity"))
            cfFatal("transient: pipe '" + name + "' declares a rod pump: the "
                    "subsurface rod pump (5.90) is a steady-network element "
                    "whose rate is its stroke times its barrel; the transient "
                    "row advances a momentum equation, it has no polished rod "
                    "to pull a card against");
        p.roughness = pd.getScalarOr("roughness", 0.0);
        p.elevation = pd.getScalarOr("elevation", 0.0);
        p.initialFlow = pd.getScalarOr("initialFlow", 0.0);
        p.valveK = pd.getScalarOr("valve", 0.0);
        p.coupledResistance = pd.getScalarOr("coupledResistance", 0.0);
        p.flowIndex = pd.getScalarOr("flowIndex", 0.0);
        p.consistency = pd.getScalarOr("consistency", 0.0);
        p.betaIn = pd.getScalarOr("betaIn", 1.0);
        p.betaOut = pd.getScalarOr("betaOut", 1.0);
        p.gelStrength = pd.getScalarOr("gelStrength", 0.0);
        p.gelTime = pd.getScalarOr("gelTime", 0.0);
        p.creepRate = pd.getScalarOr("creepRate", 0.0);
        p.yieldStress = pd.getScalarOr("yieldStress", 0.0);
        if (pd.has("valveClose"))
        {
            const Dictionary& vc = pd.getDict("valveClose");
            p.valveCloseStart = vc.getScalar("start");
            p.valveCloseEnd = vc.getScalar("end");
            p.valveKClosed = vc.getScalar("closedK");
        }
        if (pd.has("valveOpen"))
        {
            const Dictionary& vo = pd.getDict("valveOpen");
            p.valveOpenStart = vo.getScalar("start");
            p.valveOpenEnd = vo.getScalar("end");
        }
        // The opt-in multi-segment valve schedule: one table prices K(t) for
        // the whole run, so the single-cycle keys (and the open-valve `valve`
        // key the table would leave unread) cannot be declared beside it --
        // two answers for one valve is not a schedule, it is a dispute.
        if (pd.has("valveSchedule"))
        {
            if (pd.has("valveClose") || pd.has("valveOpen"))
                cfFatal("transient: pipe '" + p.name + "' declares both"
                        " valveSchedule and the single-cycle valveClose/"
                        "valveOpen keys: one schedule prices the valve");
            if (pd.has("valve"))
                cfFatal("transient: pipe '" + p.name + "' declares both"
                        " valveSchedule and an open-valve `valve` key: the"
                        " event table prices K(t) itself, so the open value"
                        " would be read by nothing");
            const Dictionary& vs = pd.getDict("valveSchedule");
            const Dictionary& ev = vs.getDict("events");
            if (ev.size() == 0)
                cfFatal("transient: pipe '" + p.name + "' declares an empty"
                        " valveSchedule events table");
            for (const std::string& ek : ev.keys())
            {
                const Dictionary& pt = ev.getDict(ek);
                if (!pt.has("t") || !pt.has("k"))
                    cfFatal("transient: pipe '" + p.name + "' valveSchedule"
                            " event needs a `t` and a `k`");
                const double t = pt.getScalar("t");
                const double kk = pt.getScalar("k");
                if (!std::isfinite(t) || !std::isfinite(kk))
                    cfFatal("transient: pipe '" + p.name + "' valveSchedule"
                            " event values must be finite numbers");
                if (t < 0)
                    cfFatal("transient: pipe '" + p.name + "' valveSchedule"
                            " event time must be >= 0 (the march starts at"
                            " t = 0)");
                if (kk < 0)
                    cfFatal("transient: pipe '" + p.name + "' valveSchedule"
                            " event K must be >= 0");
                if (!p.valveEventT.empty() && t < p.valveEventT.back())
                    cfFatal("transient: pipe '" + p.name + "' valveSchedule"
                            " events must be in ascending time order (equal"
                              " times are a step)");
                p.valveEventT.push_back(t);
                p.valveEventK.push_back(kk);
            }
            p.valveEventsOn = true;
        }
        if (p.valveK < 0 || p.valveKClosed < 0)
            cfFatal("transient: pipe '" + p.name + "' valve K must be >= 0");
        if (p.valveCloseStart > p.valveCloseEnd)
            cfFatal("transient: pipe '" + p.name
                    + "' valveClose start must be <= end");
        if (p.valveOpenStart >= 0)
        {
            if (p.valveCloseStart < 0)
                cfFatal("transient: pipe '" + p.name
                        + "' valveOpen needs a valveClose block");
            if (p.valveOpenStart < p.valveCloseEnd)
                cfFatal("transient: pipe '" + p.name
                        + "' valveOpen must not start before valveClose ends");
            if (p.valveOpenStart > p.valveOpenEnd)
                cfFatal("transient: pipe '" + p.name
                        + "' valveOpen start must be <= end");
        }
        c.pipes.push_back(p);
    }

    if (caseDict.has("sources"))
    {
        const Dictionary& src = caseDict.getDict("sources");
        for (const std::string& name : src.keys())
        {
            const Dictionary& sd = src.getDict(name);
            TransientSource s;
            s.nodeName = sd.getWord("node");
            s.fraction = sd.getScalarOr("fraction", 1.0);
            const bool hasRate = sd.has("rate");
            const bool hasSchedule = sd.has("schedule");
            if (hasRate) s.rate = sd.getScalar("rate");
            if (hasSchedule)
            {
                const Dictionary& sc = sd.getDict("schedule");
                s.openRate = sc.getScalar("openRate");
                if (sc.has("close"))
                {
                    const Dictionary& cl = sc.getDict("close");
                    s.start = cl.getScalar("start");
                    s.end = cl.getScalar("end");
                    if (s.start > s.end)
                        cfFatal("transient: source '" + name
                                + "' schedule.close start must be <= end");
                }
                if (sc.has("open"))
                {
                    if (s.start < 0)
                        cfFatal("transient: source '" + name
                                + "' schedule.open needs a schedule.close block");
                    const Dictionary& op = sc.getDict("open");
                    s.reopenStart = op.getScalar("start");
                    s.reopenEnd = op.getScalar("end");
                    if (s.reopenStart < s.end)
                        cfFatal("transient: source '" + name
                                + "' schedule.open must not start before"
                                  " schedule.close ends");
                    if (s.reopenStart > s.reopenEnd)
                        cfFatal("transient: source '" + name
                                + "' schedule.open start must be <= end");
                }
            }
            if (!hasRate && !hasSchedule)
                cfFatal("transient: source '" + name
                        + "' needs a `rate` or a `schedule` block");
            c.sources.push_back(s);
        }
    }

    // The opt-in multi-segment source schedule (top-level, targets one node
    // by name): rate(t) through the event points, the same interpolation the
    // valve table reads. The node must exist and must not also be priced by
    // a sources{} entry -- two rate accounts on one node would sum into a
    // schedule nobody declared.
    if (caseDict.has("sourceSchedule"))
    {
        const Dictionary& ss = caseDict.getDict("sourceSchedule");
        TransientSourceSchedule s;
        s.nodeName = ss.getWord("node");
        bool nodeKnown = false;
        for (const TransientNodeConfig& n : c.nodes)
            if (n.name == s.nodeName) { nodeKnown = true; break; }
        if (!nodeKnown)
            cfFatal("transient: sourceSchedule names node '" + s.nodeName
                    + "' but no such node is declared");
        for (const TransientSource& src : c.sources)
            if (src.nodeName == s.nodeName)
                cfFatal("transient: node '" + s.nodeName + "' is priced by"
                        " both a sources{} entry and a sourceSchedule: one"
                          " rate account per node");
        const Dictionary& ev = ss.getDict("events");
        if (ev.size() == 0)
            cfFatal("transient: sourceSchedule for node '" + s.nodeName
                    + "' declares an empty events table");
        for (const std::string& ek : ev.keys())
        {
            const Dictionary& pt = ev.getDict(ek);
            if (!pt.has("t") || !pt.has("rate"))
                cfFatal("transient: sourceSchedule event needs a `t` and a"
                        " `rate`");
            const double t = pt.getScalar("t");
            const double r = pt.getScalar("rate");
            if (!std::isfinite(t) || !std::isfinite(r))
                cfFatal("transient: sourceSchedule event values must be"
                        " finite numbers");
            if (t < 0)
                cfFatal("transient: sourceSchedule event time must be >= 0"
                        " (the march starts at t = 0)");
            if (r < 0)
                cfFatal("transient: sourceSchedule event rate must be >= 0:"
                        " the table is an injection schedule");
            if (!s.eventT.empty() && t < s.eventT.back())
                cfFatal("transient: sourceSchedule events must be in"
                        " ascending time order (equal times are a step)");
            s.eventT.push_back(t);
            s.eventRate.push_back(r);
        }
        c.sourceSchedules.push_back(s);
    }

    // ---- E7a (§5.177): opt-in wellbore equipment ---------------------------
    // Three top-level blocks keyed by element name. Each element targets one
    // pipe by `pipe`; at most one element per pipe, none on a pipe with a base
    // `valve` declaration or a `coupledResistance`. Every mistake is loud.
    auto pipeExists = [&](const std::string& n) {
        for (const TransientPipeConfig& p : c.pipes) if (p.name == n) return true;
        return false;
    };
    auto hasBaseLoss = [&](const std::string& n) {
        for (const TransientPipeConfig& p : c.pipes)
            if (p.name == n)
                return p.valveK > 0 || p.valveCloseStart >= 0
                    || p.valveOpenStart >= 0 || p.valveEventsOn
                    || p.coupledResistance > 0;
        return false;
    };
    std::vector<std::string> claimed;
    auto claimTarget = [&](const std::string& pipe, const std::string& what) {
        if (!pipeExists(pipe))
            cfFatal("transient: " + what + " targets unknown pipe '" + pipe + "'");
        if (hasBaseLoss(pipe))
            cfFatal("transient: pipe '" + pipe + "' declares a base valve or a"
                    " coupled resistance and a wellbore element: one local-loss"
                    " account per pipe");
        for (const std::string& s : claimed)
            if (s == pipe)
                cfFatal("transient: pipe '" + pipe + "' is targeted by two"
                        " wellbore elements: one per pipe");
        claimed.push_back(pipe);
    };

    if (caseDict.has("sssv"))
    {
        const Dictionary& blk = caseDict.getDict("sssv");
        for (const std::string& name : blk.keys())
        {
            const Dictionary& sd = blk.getDict(name);
            TransientSssvConfig s;
            s.name = name;
            s.pipeName = sd.getWord("pipe");
            claimTarget(s.pipeName, "sssv '" + name + "'");
            s.openK = sd.getScalarOr("openK", 0.0);
            s.closedK = sd.getScalarOr("closedK", 1.0e12);
            if (sd.has("close"))
            {
                const Dictionary& cl = sd.getDict("close");
                s.closeStart = cl.getScalar("start");
                s.closeEnd = cl.getScalar("end");
            }
            if (sd.has("reopen"))
            {
                const Dictionary& op = sd.getDict("reopen");
                s.reopenStart = op.getScalar("start");
                s.reopenEnd = op.getScalar("end");
            }
            if (s.openK < 0 || s.closedK < 0)
                cfFatal("transient: sssv '" + name + "' K must be >= 0");
            if (s.closedK < s.openK)
                cfFatal("transient: sssv '" + name + "' closedK must be >= openK");
            if (s.closeStart < 0)
                cfFatal("transient: sssv '" + name + "' needs a close block");
            if (s.closeStart > s.closeEnd)
                cfFatal("transient: sssv '" + name + "' close start must be <= end");
            if (s.reopenStart >= 0)
            {
                if (s.reopenStart < s.closeEnd)
                    cfFatal("transient: sssv '" + name + "' reopen must not start"
                            " before close ends");
                if (s.reopenStart > s.reopenEnd)
                    cfFatal("transient: sssv '" + name + "' reopen start must be <= end");
            }
            c.sssvs.push_back(s);
        }
    }

    if (caseDict.has("choke"))
    {
        const Dictionary& blk = caseDict.getDict("choke");
        for (const std::string& name : blk.keys())
        {
            const Dictionary& cd = blk.getDict(name);
            TransientChokeConfig ch;
            ch.name = name;
            ch.pipeName = cd.getWord("pipe");
            claimTarget(ch.pipeName, "choke '" + name + "'");
            const bool hasK = cd.has("k");
            const bool hasD = cd.has("diameter");
            if (hasK == hasD)
                cfFatal("transient: choke '" + name + "' declares both k and"
                        " diameter (or neither): one resistance law per choke");
            if (hasK)
            {
                ch.k = cd.getScalar("k");
                if (!(ch.k > 0))
                    cfFatal("transient: choke '" + name + "' k must be positive");
            }
            else
            {
                ch.diameter = cd.getScalar("diameter");
                ch.dischargeCoeff = cd.getScalarOr("dischargeCoeff", 0.62);
                if (!(ch.diameter > 0))
                    cfFatal("transient: choke '" + name + "' diameter must be positive");
                if (!(ch.dischargeCoeff > 0 && ch.dischargeCoeff <= 1))
                    cfFatal("transient: choke '" + name + "' dischargeCoeff must"
                            " be in (0,1]");
            }
            c.chokes.push_back(ch);
        }
    }

    if (caseDict.has("esp"))
    {
        const Dictionary& blk = caseDict.getDict("esp");
        for (const std::string& name : blk.keys())
        {
            const Dictionary& ed = blk.getDict(name);
            TransientEspConfig e;
            e.name = name;
            e.pipeName = ed.getWord("pipe");
            claimTarget(e.pipeName, "esp '" + name + "'");
            e.stages = ed.getScalarOr("stages", 1.0);
            e.shutoffHead = ed.getScalarOr("shutoffHead", 0.0);
            e.headAtRate = ed.getScalarOr("headAtRate", 0.0);
            e.rate = ed.getScalarOr("rate", 0.0);
            e.speed0 = ed.getScalarOr("speed0", 1.0);
            e.speedCmd = ed.getScalarOr("speedCmd", 1.0);
            e.speedCmdTime = ed.getScalarOr("speedCmdTime", 0.0);
            e.speedTau = ed.getScalarOr("speedTau", 0.0);
            if (!(e.stages > 0))
                cfFatal("transient: esp '" + name + "' stages must be positive");
            if (!(e.shutoffHead > 0))
                cfFatal("transient: esp '" + name + "' needs a positive shutoffHead");
            if (!(e.headAtRate >= 0 && e.headAtRate < e.shutoffHead))
                cfFatal("transient: esp '" + name + "' headAtRate must be in"
                        " [0, shutoffHead)");
            if (e.headAtRate > 0 && !(e.rate > 0))
                cfFatal("transient: esp '" + name + "' declares headAtRate"
                        " without a positive rate");
            if (!(e.speed0 > 0) || !(e.speedCmd > 0))
                cfFatal("transient: esp '" + name + "' speeds must be positive");
            if (e.speedCmdTime < 0)
                cfFatal("transient: esp '" + name + "' speedCmdTime must be >= 0");
            if (e.speedTau < 0)
                cfFatal("transient: esp '" + name + "' speedTau must be >= 0");
            c.esps.push_back(e);
        }
    }

    // ---- 5.185: the opt-in well-control reservoir inflow -------------------
    // The node-level wiring (freedom, compliance, one rate account) is checked
    // by the constructor, which both entry points go through; here the block's
    // own grammar and numbers are priced, with every mistake loud.
    if (caseDict.has("reservoir"))
    {
        const Dictionary& blk = caseDict.getDict("reservoir");
        for (const std::string& name : blk.keys())
        {
            const Dictionary& rd = blk.getDict(name);
            TransientReservoirConfig r;
            r.name = name;
            r.nodeName = rd.getWord("node");
            r.pressure = rd.getScalar("pressure");
            r.productivityIndex = rd.getScalar("productivityIndex");
            r.gasFraction = rd.getScalarOr("gasFraction", 0.0);
            r.depth = rd.getScalar("depth");
            if (!(r.pressure > 0))
                cfFatal("transient: reservoir '" + name + "' pressure must be"
                        " positive");
            if (!(r.productivityIndex > 0))
                cfFatal("transient: reservoir '" + name
                        + "' productivityIndex must be positive");
            if (!(r.gasFraction >= 0) || !(r.gasFraction < 1))
                cfFatal("transient: reservoir '" + name + "' gasFraction must"
                        " be in [0,1)");
            if (!(r.depth > 0))
                cfFatal("transient: reservoir '" + name + "' depth must be"
                        " positive (the node's TVD): the kill-mud weight is"
                        " p_res / (g depth)");
            c.reservoirs.push_back(r);
        }
    }
    return c;
}

TransientTimeControls loadTransientControls(const Dictionary& caseDict)
{
    TransientTimeControls tc;
    if (!caseDict.has("time")) cfFatal("transient: case has no time{} block");
    const Dictionary& td = caseDict.getDict("time");
    tc.dt = td.getScalar("dt");
    tc.endTime = td.getScalar("endTime");
    tc.reportEvery = static_cast<int>(td.getScalarOr("reportEvery", 1.0));
    // 5.68 / D5: cfl switches the loop to CFL-adaptive stepping. Absent means
    // the fixed-step march; present-but-nonpositive is a broken case.
    tc.cfl = td.getScalarOr("cfl", 0.0);
    if (td.has("cfl") && !(tc.cfl > 0))
        cfFatal("transient: cfl must be positive when declared");
    if (tc.cfl > 0)
    {
        tc.dtMin = td.getScalarOr("dtMin", tc.dt);
        tc.dtMax = td.getScalarOr("dtMax", tc.endTime);
    }
    return tc;
}

} // namespace cf
