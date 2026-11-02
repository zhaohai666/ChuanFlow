#include "cf/1d/Transient.h"

#include <cmath>
#include <limits>

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
        net_.addPipe(geo);
        PipeDyn d;
        d.geo = net_.pipes.back();
        d.qPrev = pc.initialFlow;
        d.inertia = fluid_.rho * d.geo.length / pipeArea(d.geo.diameter);
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

double TransientNetwork::step(double dt)
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
            const double fr = frictionDp(pd.geo, fluid_, x[k]);
            const double slope = frictionSlope(pd.geo, fluid_, x[k]);
            r[k] = pd.inertia / dt * (x[k] - pd.qPrev) - dp + grav + fr;
            J[k][k] = pd.inertia / dt + slope;
            const NodeDyn& nf = nodes_[static_cast<std::size_t>(pd.geo.from)];
            const NodeDyn& nb = nodes_[static_cast<std::size_t>(pd.geo.to)];
            if (!nf.fixed) J[k][nP + static_cast<std::size_t>(nf.unk)] -= 1.0;
            if (!nb.fixed) J[k][nP + static_cast<std::size_t>(nb.unk)] += 1.0;
            const double scale = pd.inertia / dt * std::max(std::abs(x[k]), std::abs(pd.qPrev))
                               + std::abs(dp) + std::abs(fr) + 1.0;
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

    // commit the solved state
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
    return maxRel;
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
