#include "cf/1d/Network.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "cf/common/Error.h"

namespace cf
{
namespace
{

constexpr double kGravity = 9.80665;

double pipeArea(const NetworkPipe& p)
{
    return 0.25 * M_PI * p.diameter * p.diameter;
}

/// Hagen-Poiseuille slope pi D^4 / (128 mu L): the maximum conductance of a
/// pipe (laminar limit), used as the fallback when the difference quotient
/// degenerates at dp ~ 0.
double laminarConductance(const NetworkPipe& p, const FluidProperties& fluid)
{
    return M_PI * std::pow(p.diameter, 4) / (128.0 * fluid.mu * p.length);
}

struct RawFlow
{
    double velocity = 0; // m/s, signed from -> to
    double reynolds = 0;
    double friction = 0;
};

/// Solve |dpFriction| = f(Re) L/D rho v^2 / 2 for |v| and re-apply the sign.
RawFlow pipeFlowRaw(const NetworkPipe& pipe, const FluidProperties& fluid,
                    double dpTotal)
{
    RawFlow r;
    const double dpFriction = dpTotal - fluid.rho * kGravity * pipe.elevation;
    const double target = std::abs(dpFriction);
    if (target <= 0 || !(pipe.length > 0) || !(pipe.diameter > 0)) return r;

    const double relRough = pipe.roughness / pipe.diameter;
    // seed with the laminar closed form, then fixed-point on Churchill:
    // v_new = sqrt(2 D target / (rho f(v) L)); the map is contractive at both
    // roots (derivative 1/2 in the laminar regime, <1 once turbulent)
    double vm = target * pipe.diameter * pipe.diameter
              / (32.0 * fluid.mu * pipe.length);
    for (int it = 0; it < 100; ++it)
    {
        const double re = fluid.rho * vm * pipe.diameter / fluid.mu;
        const double f = churchillFriction(re, relRough);
        const double vmNew = std::sqrt(2.0 * pipe.diameter * target
                                       / (f * pipe.length * fluid.rho));
        const double err = std::abs(vmNew - vm) / std::max(vmNew, 1e-300);
        vm = vmNew;
        if (err < 1e-13) break;
    }
    r.reynolds = fluid.rho * vm * pipe.diameter / fluid.mu;
    r.friction = churchillFriction(r.reynolds, relRough);
    r.velocity = (dpFriction < 0 ? -vm : vm);
    return r;
}

/// Dense symmetric-ish n x n solve with partial pivoting (n is tiny: a
/// pipe-network model rarely has more than a few hundred free nodes).
void solveLinear(std::vector<std::vector<double>>& a, std::vector<double>& b)
{
    const std::size_t n = b.size();
    for (std::size_t k = 0; k < n; ++k)
    {
        std::size_t piv = k;
        double best = std::abs(a[k][k]);
        for (std::size_t i = k + 1; i < n; ++i)
        {
            if (std::abs(a[i][k]) > best) { best = std::abs(a[i][k]); piv = i; }
        }
        if (best < 1e-12 * std::max(1.0, *std::max_element(
                a[piv].begin(), a[piv].end())))
        {
            cfFatal("steadyNetwork: singular Jacobian - check the network is "
                    "connected to at least one fixed-pressure node");
        }
        if (piv != k) std::swap(a[piv], a[k]), std::swap(b[piv], b[k]);
        for (std::size_t i = k + 1; i < n; ++i)
        {
            const double m = a[i][k] / a[k][k];
            if (m == 0) continue;
            for (std::size_t j = k; j < n; ++j) a[i][j] -= m * a[k][j];
            b[i] -= m * b[k];
        }
    }
    for (std::size_t i = n; i-- > 0;)
    {
        double s = b[i];
        for (std::size_t j = i + 1; j < n; ++j) s -= a[i][j] * b[j];
        b[i] = s / a[i][i];
    }
}

} // namespace

label OneDNetwork::addNode(const std::string& name)
{
    for (label i = 0; i < static_cast<label>(nodeNames.size()); ++i)
    {
        if (nodeNames[i] == name) return i;
    }
    nodeNames.push_back(name);
    nodePressureFixed.push_back(std::numeric_limits<double>::quiet_NaN());
    nodeExternalRate.push_back(0.0);
    return static_cast<label>(nodeNames.size()) - 1;
}

label OneDNetwork::nodeIndex(const std::string& name) const
{
    for (label i = 0; i < static_cast<label>(nodeNames.size()); ++i)
    {
        if (nodeNames[i] == name) return i;
    }
    cfFatal("network: unknown node '" + name + "'");
    return -1;
}

label OneDNetwork::addPipe(const NetworkPipe& pipe)
{
    NetworkPipe p = pipe;
    p.from = nodeIndex(p.fromName);
    p.to = nodeIndex(p.toName);
    if (p.from == p.to) cfFatal("network: pipe '" + p.name + "' connects a node to itself");
    if (!(p.diameter > 0) || !(p.length > 0))
        cfFatal("network: pipe '" + p.name + "' needs positive length and diameter");
    pipes.push_back(p);
    return static_cast<label>(pipes.size()) - 1;
}

double churchillFriction(double reynolds, double relativeRoughness)
{
    const double re = std::max(reynolds, 1e-12);
    const double a1 = std::pow(7.0 / re, 0.9) + 0.27 * relativeRoughness;
    if (a1 <= 0) cfFatal("churchillFriction: invalid Re/roughness combination");
    const double a = std::pow(2.457 * std::log(1.0 / a1), 16.0);
    const double b = std::pow(37530.0 / re, 16.0);
    return 8.0 * std::pow(std::pow(8.0 / re, 12.0) + std::pow(a + b, -1.5), 1.0 / 12.0);
}

void solvePipeFlow(const NetworkPipe& pipe, const FluidProperties& fluid,
                   double dpTotal, PipeFlowResult& out)
{
    const RawFlow r = pipeFlowRaw(pipe, fluid, dpTotal);
    out.name = pipe.name;
    out.dpTotal = dpTotal;
    out.dpFriction = dpTotal - fluid.rho * kGravity * pipe.elevation;
    out.velocity = r.velocity;
    out.flowRate = r.velocity * pipeArea(pipe);
    out.reynolds = r.reynolds;
    out.frictionFactor = r.friction;

    // conductance dQ/d(dpTotal) by central difference on the monotone law;
    // floor/replace with the laminar slope if the quotient degenerates
    const double scale = std::max(std::abs(dpTotal),
                                  fluid.rho * kGravity * std::abs(pipe.elevation));
    const double h = std::max(1e-6 * std::abs(out.dpFriction), 1e-7 * scale);
    const double lam = laminarConductance(pipe, fluid);
    if (h > 0)
    {
        const RawFlow rp = pipeFlowRaw(pipe, fluid, dpTotal + h);
        const RawFlow rm = pipeFlowRaw(pipe, fluid, dpTotal - h);
        const double g = (rp.velocity - rm.velocity) * pipeArea(pipe) / (2.0 * h);
        out.conductance = (std::isfinite(g) && g > 0) ? std::min(g, lam) : lam;
    }
    else
    {
        out.conductance = lam;
    }
}

SteadyNetworkResult solveSteadyNetwork(const OneDNetwork& net,
                                       const SteadyControls& ctl)
{
    SteadyNetworkResult res;
    const std::size_t n = net.nNodes();
    res.nodePressure.assign(n, 0.0);
    res.pipeFlows.resize(net.pipes.size());

    // unknown node pressures
    std::vector<label> unk(n, -1);
    std::size_t nFree = 0;
    double pSum = 0;
    std::size_t nFixed = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (std::isnan(net.nodePressureFixed[i])) unk[i] = static_cast<label>(nFree++);
        else { pSum += net.nodePressureFixed[i]; ++nFixed; }
    }
    if (nFixed == 0 && nFree > 0)
        cfFatal("steadyNetwork: no fixed-pressure node (no reference)");
    const double pInit = nFixed ? pSum / static_cast<double>(nFixed) : 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        res.nodePressure[i] = std::isnan(net.nodePressureFixed[i])
                                  ? pInit : net.nodePressureFixed[i];
    }

    // total injected rate = residual normalisation scale
    double totalIn = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (net.nodeExternalRate[i] > 0) totalIn += net.nodeExternalRate[i];
    }
    const double qScale = std::max(totalIn, 1e-300);

    // Newton on free-node pressures
    std::vector<double> p(nFree, pInit);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (unk[i] >= 0) p[static_cast<std::size_t>(unk[i])] = res.nodePressure[i];
    }

    auto evaluate = [&](const std::vector<double>& pv,
                        std::vector<double>& r,
                        std::vector<std::vector<double>>& j,
                        bool withJacobian) -> double
    {
        std::vector<double> pAll(n, pInit);
        for (std::size_t i = 0; i < n; ++i)
        {
            pAll[i] = std::isnan(net.nodePressureFixed[i])
                        ? pv[static_cast<std::size_t>(unk[i])]
                        : net.nodePressureFixed[i];
        }
        r.assign(nFree, 0.0);
        if (withJacobian) j.assign(nFree, std::vector<double>(nFree, 0.0));

        double maxRes = 0;
        for (std::size_t k = 0; k < net.pipes.size(); ++k)
        {
            const NetworkPipe& pipe = net.pipes[k];
            PipeFlowResult fr;
            solvePipeFlow(pipe, net.fluid, pAll[static_cast<std::size_t>(pipe.from)]
                                        - pAll[static_cast<std::size_t>(pipe.to)], fr);
            if (withJacobian) res.pipeFlows[k] = fr;
            const label a = unk[static_cast<std::size_t>(pipe.from)];
            const label b = unk[static_cast<std::size_t>(pipe.to)];
            // R_a -= Q ; R_b += Q
            if (a >= 0) r[static_cast<std::size_t>(a)] -= fr.flowRate;
            if (b >= 0) r[static_cast<std::size_t>(b)] += fr.flowRate;
            if (withJacobian)
            {
                const double c = fr.conductance;
                if (a >= 0)
                {
                    j[a][a] -= c;
                    if (b >= 0) { j[a][b] += c; j[b][a] += c; j[b][b] -= c; }
                }
            }
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            if (unk[i] >= 0) r[static_cast<std::size_t>(unk[i])] += net.nodeExternalRate[i];
        }
        for (double v : r) maxRes = std::max(maxRes, std::abs(v));
        return maxRes;
    };

    std::vector<double> r;
    std::vector<std::vector<double>> j;
    double maxRes = evaluate(p, r, j, true);
    res.residual = maxRes;
    res.relativeResidual = maxRes / qScale;

    int iter = 0;
    while (res.relativeResidual > ctl.tol && iter < ctl.maxIter)
    {
        if (nFree > 0)
        {
            std::vector<double> negR(r.size());
            for (std::size_t i = 0; i < r.size(); ++i) negR[i] = -r[i];
            solveLinear(j, negR);
            const std::vector<double> pOld = p;
            for (std::size_t i = 0; i < nFree; ++i) p[i] += ctl.relax * negR[i];
            maxRes = evaluate(p, r, j, true);
            // crude line search: the friction law is monotone, so a full
            // Newton step normally wins; halve only on clear growth
            int halvings = 0;
            while (maxRes > 1.5 * res.residual && halvings < 6)
            {
                for (std::size_t i = 0; i < nFree; ++i)
                    p[i] = pOld[i] + ctl.relax * negR[i] * std::pow(0.5, halvings + 1);
                maxRes = evaluate(p, r, j, true);
                ++halvings;
            }
            res.residual = maxRes;
            res.relativeResidual = maxRes / qScale;
        }
        ++iter;
        if (nFree == 0) break;
    }
    res.iterations = iter;
    res.converged = res.relativeResidual <= ctl.tol;

    for (std::size_t i = 0; i < n; ++i)
    {
        if (unk[i] >= 0) res.nodePressure[i] = p[static_cast<std::size_t>(unk[i])];
    }
    // final pipe states at the converged pressures
    for (std::size_t k = 0; k < net.pipes.size(); ++k)
    {
        const NetworkPipe& pipe = net.pipes[k];
        solvePipeFlow(pipe, net.fluid,
                      res.nodePressure[static_cast<std::size_t>(pipe.from)]
                    - res.nodePressure[static_cast<std::size_t>(pipe.to)],
                      res.pipeFlows[k]);
    }
    return res;
}

OneDNetwork loadNetworkFromDict(const Dictionary& caseDict)
{
    OneDNetwork net;

    if (caseDict.has("fluid"))
    {
        const Dictionary& fd = caseDict.getDict("fluid");
        net.fluid.rho = fd.getScalar("rho");
        net.fluid.mu = fd.getScalar("mu");
    }
    if (!(net.fluid.rho > 0) || !(net.fluid.mu > 0))
        cfFatal("steadyNetwork: fluid needs positive rho and mu");

    if (caseDict.has("nodes"))
    {
        const Dictionary& nodes = caseDict.getDict("nodes");
        for (const std::string& name : nodes.keys())
        {
            const label idx = net.addNode(name);
            const Dictionary& nd = nodes.getDict(name);
            if (nd.has("pressure")) net.nodePressureFixed[idx] = nd.getScalar("pressure");
        }
    }

    if (!caseDict.has("pipes")) cfFatal("steadyNetwork: case has no pipes{} block");
    const Dictionary& pipes = caseDict.getDict("pipes");
    for (const std::string& name : pipes.keys())
    {
        const Dictionary& pd = pipes.getDict(name);
        NetworkPipe p;
        p.name = name;
        p.fromName = pd.getWord("from");
        p.toName = pd.getWord("to");
        p.length = pd.getScalar("length");
        p.diameter = pd.getScalar("diameter");
        p.roughness = pd.getScalarOr("roughness", 0.0);
        p.elevation = pd.getScalarOr("elevation", 0.0);
        net.addPipe(p);
    }

    if (caseDict.has("sources"))
    {
        const Dictionary& src = caseDict.getDict("sources");
        for (const std::string& name : src.keys())
        {
            const Dictionary& sd = src.getDict(name);
            net.setNodeRate(sd.getWord("node"), sd.getScalar("rate"));
        }
    }

    return net;
}

SteadyControls loadSteadyControls(const Dictionary& caseDict)
{
    SteadyControls ctl;
    if (caseDict.has("solve"))
    {
        const Dictionary& sd = caseDict.getDict("solve");
        ctl.maxIter = static_cast<int>(sd.getScalarOr("maxIter", ctl.maxIter));
        ctl.tol = sd.getScalarOr("tol", ctl.tol);
        ctl.relax = sd.getScalarOr("relax", ctl.relax);
    }
    return ctl;
}

} // namespace cf
