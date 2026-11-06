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

/// Pressure drop consumed by friction for a signed flow (Darcy-Weisbach
/// with the same Churchill factor as the steady kernel; odd in Q).
double frictionDp(const NetworkPipe& geo, const FluidProperties& fl, double q)
{
    const double area = pipeArea(geo.diameter);
    const double v = std::abs(q) / area;
    if (v == 0) return 0;
    const double re = fl.rho * v * geo.diameter / fl.mu;
    const double f = churchillFriction(re, geo.roughness / geo.diameter);
    const double dp = f * geo.length / geo.diameter * 0.5 * fl.rho * v * v;
    return q > 0 ? dp : -dp;
}

/// d(dpFriction)/dQ by central differencing on the monotone odd law;
/// falls back to the laminar slope 128 mu L / (pi D^4) near Q = 0.
double frictionSlope(const NetworkPipe& geo, const FluidProperties& fl, double q)
{
    const double lam = 128.0 * fl.mu * geo.length
                     / (M_PI * geo.diameter * geo.diameter * geo.diameter * geo.diameter);
    const double h = std::max(1e-5 * std::abs(q), 1e-12);
    const double d = (frictionDp(geo, fl, q + h) - frictionDp(geo, fl, q - h)) / (2.0 * h);
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
        nodes_[static_cast<std::size_t>(i)].source += s.rate;
        nodes_[static_cast<std::size_t>(i)].sourceFraction = s.fraction;
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

    double maxRel = 0;
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
            const double fr = R > 0 ? R * x[k] : frictionDp(pd.geo, fluid_, x[k]);
            const double slope = R > 0
                                    ? 1.0 / R
                                    : frictionSlope(pd.geo, fluid_, x[k]);
            const double ar = pipeArea(pd.geo.diameter);
            const double vd = kv[k] * fluid_.rho * std::abs(x[k]) * x[k]
                            / (2.0 * ar * ar);
            r[k] = pd.inertia / dt * (x[k] - pd.qPrev) - dp + grav + fr + vd;
            J[k][k] = pd.inertia / dt + slope
                    + kv[k] * fluid_.rho * std::abs(x[k]) / (ar * ar);
            const NodeDyn& nf = nodes_[static_cast<std::size_t>(pd.geo.from)];
            const NodeDyn& nb = nodes_[static_cast<std::size_t>(pd.geo.to)];
            if (!nf.fixed) J[k][nP + static_cast<std::size_t>(nf.unk)] -= 1.0;
            if (!nb.fixed) J[k][nP + static_cast<std::size_t>(nb.unk)] += 1.0;
            const double scale = pd.inertia / dt * std::max(std::abs(x[k]), std::abs(pd.qPrev))
                               + std::abs(dp) + std::abs(fr) + std::abs(vd) + 1.0;
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
            s.rate = sd.getScalar("rate");
            s.fraction = sd.getScalarOr("fraction", 1.0);
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
