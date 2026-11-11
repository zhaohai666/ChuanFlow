#include "cf/1d/Transient.h"

#include <cmath>
#include <limits>
#include <string>
#include <utility>

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
bool holdsStatic(const NetworkPipe& geo, double qPrev, double drive)
{
    if (!(geo.yieldStress > 0)) return false;
    return qPrev == 0 && std::abs(drive) <= yieldThresholdPa(geo);
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

/// A coupled pipe's loss is its partner's, so no local loss may share it.
void requireNoValve(const std::string& name, double kOpen, double closeStart,
                    double openStart)
{
    if (kOpen > 0 || closeStart >= 0 || openStart >= 0)
        cfFatal("transient: pipe '" + name
                + "' is a coupled element: its loss is the partner model's, so"
                  " it cannot carry a valve");
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
        // schedule can move it. Pre-fill with the constant-rate value for any
        // caller that reads `nd.source` before the first step().
        nodes_[e.nodeIndex].sourceFraction = e.fraction;
        nodes_[e.nodeIndex].source =
            (e.start < 0) ? e.rate : 0.0;
    }

    for (const TransientPipeConfig& pc : c.pipes)
    {
        NetworkPipe geo;
        geo.name = pc.name;
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
        if (pc.flowIndex > 0 && pc.coupledResistance > 0)
            cfFatal("transient: pipe '" + pc.name + "' is a coupled element:"
                    " its loss is the partner model's, so it cannot carry a"
                    " rheology either");
        if (pc.yieldStress < 0)
            cfFatal("transient: pipe '" + pc.name + "' has a negative yield"
                    " stress: the threshold is the stress a static wall holds,"
                    " and a negative one would make the held set empty");
        if (pc.yieldStress > 0 && pc.coupledResistance > 0)
            cfFatal("transient: pipe '" + pc.name + "' is a coupled element:"
                    " its loss is the partner model's, so it cannot carry a"
                    " threshold either");
        if (pc.coupledResistance < 0)
            cfFatal("transient: pipe '" + pc.name + "' has a negative coupled"
                    " resistance");
        if (pc.coupledResistance > 0)
            requireNoValve(pc.name, pc.valveK, pc.valveCloseStart,
                           pc.valveOpenStart);
        net_.addPipe(geo);
        PipeDyn d;
        d.geo = net_.pipes.back();
        d.qPrev = pc.initialFlow;
        d.inertia = fluid_.rho * d.geo.length / pipeArea(d.geo.diameter);
        d.closeStart = pc.valveCloseStart;
        d.closeEnd = pc.valveCloseEnd;
        d.kClosed = pc.valveKClosed;
        d.openStart = pc.valveOpenStart;
        d.openEnd = pc.valveOpenEnd;
        pipes_.push_back(d);
        pipeNames_.push_back(pc.name);
    }
    if (c.pipes.empty()) cfFatal("transient: case has no pipes");

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
    // constant during the Newton solve
    const double tNew = t_ + dt;
    std::vector<double> kv(nP);
    for (std::size_t k = 0; k < nP; ++k)
        kv[k] = valveKAt(pipes_[k].geo.valveK, pipes_[k].kClosed,
                         pipes_[k].closeStart, pipes_[k].closeEnd,
                         pipes_[k].openStart, pipes_[k].openEnd, tNew);

    // source rate at the same implicit new time level; reset-then-sum so a
    // node with multiple scheduled sources sees their sum (not a carry-over).
    for (const SourceEntry& e : sourceEntries_)
        nodes_[e.nodeIndex].source = 0.0;
    for (const SourceEntry& e : sourceEntries_)
        nodes_[e.nodeIndex].source += sourceRateAt(e.rate, e.openRate,
                                                   e.start, e.end,
                                                   e.reopenStart, e.reopenEnd,
                                                   tNew);

    double maxRel = 0;

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
        pd.held = holdsStatic(pd.geo, pd.qPrev, drive);
    }

    for (int newton = 0; newton < 60; ++newton)
    {
        for (double& v : r) v = 0;
        for (auto& row : J) std::fill(row.begin(), row.end(), 0.0);
        maxRel = 0;

        for (std::size_t k = 0; k < nP; ++k)
        {
            PipeDyn& pd = pipes_[k];
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
                x[k] = 0;
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
            r[k] = pd.inertia / dt * (x[k] - pd.qPrev) - dp + grav + fr + vd + cc;
            J[k][k] = pd.inertia / dt + slope
                    + kv[k] * fluid_.rho * std::abs(x[k]) / (ar * ar)
                    + 2.0 * pd.convective * std::abs(x[k]);
            const NodeDyn& nf = nodes_[static_cast<std::size_t>(pd.geo.from)];
            const NodeDyn& nb = nodes_[static_cast<std::size_t>(pd.geo.to)];
            if (!nf.fixed) J[k][nP + static_cast<std::size_t>(nf.unk)] -= 1.0;
            if (!nb.fixed) J[k][nP + static_cast<std::size_t>(nb.unk)] += 1.0;
            const double scale = pd.inertia / dt * std::max(std::abs(x[k]), std::abs(pd.qPrev))
                               + std::abs(dp) + std::abs(fr) + std::abs(vd)
                               + std::abs(cc) + 1.0;
            maxRel = std::max(maxRel, std::abs(r[k]) / scale);
        }

        for (std::size_t i = 0; i < nN; ++i)
        {
            const NodeDyn& nd = nodes_[i];
            if (nd.fixed) continue;
            const std::size_t e = nP + static_cast<std::size_t>(nd.unk);
            double sumIn = nd.source;
            for (std::size_t k = 0; k < nP; ++k)
            {
                if (pipes_[k].geo.to == static_cast<label>(i)) sumIn += x[k];
                if (pipes_[k].geo.from == static_cast<label>(i)) sumIn -= x[k];
            }
            const double stor = nd.bulk > 0 ? nd.volume / (nd.bulk * dt) * (x[e] - nd.pPrev) : 0.0;
            r[e] = stor - sumIn;
            if (nd.bulk > 0) J[e][e] = nd.volume / (nd.bulk * dt);
            for (std::size_t k = 0; k < nP; ++k)
            {
                if (pipes_[k].geo.to == static_cast<label>(i)) J[e][k] -= 1.0;
                if (pipes_[k].geo.from == static_cast<label>(i)) J[e][k] += 1.0;
            }
            const double scale = std::abs(stor) + std::abs(sumIn) + 1.0;
            maxRel = std::max(maxRel, std::abs(r[e]) / scale);
        }

        if (maxRel < 1e-11) break;
        std::vector<double> negR(nEq);
        for (std::size_t i = 0; i < nEq; ++i) negR[i] = -r[i];
        detail::solveDenseLinear(J, negR);
        for (std::size_t i = 0; i < nEq; ++i) x[i] += negR[i];
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
        pipes_[k].qPrev = x[k];
        q_[k] = x[k];
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

void TransientNetwork::setCoupledResistance(label k, double resistance)
{
    if (k < 0 || static_cast<std::size_t>(k) >= pipes_.size())
        cfFatal("transient: no pipe at index " + std::to_string(k));
    if (!(resistance > 0))
        cfFatal("transient: coupled pipe '" + pipeNames_[k] + "' needs a"
                " positive resistance, got " + std::to_string(resistance));
    const PipeDyn& pd = pipes_[static_cast<std::size_t>(k)];
    requireNoValve(pipeNames_[k], pd.geo.valveK, pd.closeStart, pd.openStart);
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
    const int nSteps = static_cast<int>(std::round(endTime / dt));
    for (int i = 1; i <= nSteps; ++i)
    {
        sim.step(dt);
        if (i % report == 0) sample();
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
        p.roughness = pd.getScalarOr("roughness", 0.0);
        p.elevation = pd.getScalarOr("elevation", 0.0);
        p.initialFlow = pd.getScalarOr("initialFlow", 0.0);
        p.valveK = pd.getScalarOr("valve", 0.0);
        p.coupledResistance = pd.getScalarOr("coupledResistance", 0.0);
        p.flowIndex = pd.getScalarOr("flowIndex", 0.0);
        p.consistency = pd.getScalarOr("consistency", 0.0);
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
    return tc;
}

} // namespace cf
