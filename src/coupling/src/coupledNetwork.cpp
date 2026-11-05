#include "cf/coupling/CoupledNetwork.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

#include "cf/common/Error.h"

namespace cf
{
namespace coupling
{
namespace
{

/// Coordinate band used to decide "same cross-section": geometric noise, not a
/// physical separation. Centroids of a structured block agree to ~1e-16, so a
/// band this wide collapses them into one section without ever merging two
/// cells of a mesh whose spacing is resolved.
scalar bandOf(scalar extent)
{
    return 1e-9 * std::max(extent, 1e-30);
}

Vector patchAreaVector(const polyMesh& mesh, label patch)
{
    const auto& patches = mesh.patches();
    const polyMesh::Patch& pt = patches[static_cast<std::size_t>(patch)];
    const auto& sf = mesh.faceSf();
    Vector s = Vector::zero();
    for (label f = pt.start; f < pt.start + pt.size; ++f)
        s += sf[static_cast<std::size_t>(f)];
    return s;
}

/// Unit vector along mesh axis `axis` (0=x, 1=y, 2=z).
Vector unitAxis(int axis)
{
    Vector v = Vector::zero();
    v[axis] = 1.0;
    return v;
}

/// Full-precision number for the fatal messages: a stalled interface residual
/// is by definition small, and six digits would print it as 0.
std::string fmt(scalar v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

} // namespace

DetailBlock::DetailBlock(polyMesh mesh, scalar nu, scalar rho,
                         const Schemes& schemes, const SimpleControls& ctl,
                         scalar probeMin, scalar probeMax)
    : mesh_(std::move(mesh)),
      nu_(nu),
      rho_(rho),
      schemes_(schemes),
      ctl_(ctl),
      probeMin_(probeMin),
      probeMax_(probeMax),
      U_(mesh_, "U"),
      p_(mesh_, "p"),
      phi_(mesh_, "phi", 0.0)
{
    if (!(nu_ > 0) || !(rho_ > 0))
        cfFatal("detailBlock: nu and rho must be positive");

    for (label patch = 0; patch < mesh_.nPatches(); ++patch)
    {
        const std::string& type = mesh_.patches()[static_cast<std::size_t>(patch)].type;
        if (type == "inlet" && inletPatch_ < 0) inletPatch_ = patch;
        else if (type == "outlet" && outletPatch_ < 0) outletPatch_ = patch;
    }
    if (inletPatch_ < 0 || outletPatch_ < 0)
        cfFatal("detailBlock: the detailed region needs one patch of type inlet"
                " and one of type outlet (found inlet "
                + std::to_string(inletPatch_) + ", outlet "
                + std::to_string(outletPatch_) + ")");

    const Vector inSf = patchAreaVector(mesh_, inletPatch_);
    const Vector outSf = patchAreaVector(mesh_, outletPatch_);
    inletArea_ = mag(inSf);
    if (!(inletArea_ > 0))
        cfFatal("detailBlock: the inlet patch has no area");
    if (mag(inSf + outSf) > 1e-6 * inletArea_)
        cfFatal("detailBlock: the inlet and outlet patches are not opposed"
                " (their area vectors do not cancel to 1e-6 of the inlet area)");

    // The flow axis and direction come from the mesh: the inlet patch's area
    // vector points OUT of the region, so its dominant component names the
    // axis and its sign says which end is upstream. Nothing here assumes the
    // detailed region lies along x.
    const Vector nrm = inSf * (1.0 / inletArea_);
    axis_ = 0;
    if (std::abs(nrm[1]) > std::abs(nrm[axis_])) axis_ = 1;
    if (std::abs(nrm[2]) > std::abs(nrm[axis_])) axis_ = 2;
    flowSign_ = nrm[axis_] < 0 ? 1 : -1;

    buildSections();
}

void DetailBlock::buildSections()
{
    const auto& centres = mesh_.cellCentres();
    const scalar extent = std::abs(mesh_.maxX()[axis_] - mesh_.minX()[axis_]);
    const scalar band = bandOf(extent);

    // Cross-sections in the order the flow meets them.
    std::vector<label> order(static_cast<std::size_t>(mesh_.nCells()));
    for (label c = 0; c < mesh_.nCells(); ++c) order[static_cast<std::size_t>(c)] = c;
    const auto coord = [&](label c) { return centres[static_cast<std::size_t>(c)][axis_]; };
    std::sort(order.begin(), order.end(), [&](label a, label b) {
        const scalar ka = coord(a), kb = coord(b);
        if (ka != kb) return ka < kb;
        return a < b;
    });
    if (flowSign_ < 0) std::reverse(order.begin(), order.end());

    std::size_t g = 0;
    while (g < order.size())
    {
        const scalar anchor = coord(order[g]);
        std::size_t h = g + 1;
        while (h < order.size() && std::abs(coord(order[h]) - anchor) <= band) ++h;
        sections_.emplace_back(order.begin() + static_cast<std::ptrdiff_t>(g),
                               order.begin() + static_cast<std::ptrdiff_t>(h));
        sectionX_.push_back(anchor);
        g = h;
    }

    // Face planes across the duct. Only faces that are predominantly normal to
    // the flow axis count: a transverse face sitting in the same coordinate
    // band would otherwise leak its own (cross-stream) flux into the section
    // total, and a plane sum is supposed to be the rate carried downstream.
    const auto& fc = mesh_.faceCentres();
    const auto& sf = mesh_.faceSf();
    std::vector<std::pair<scalar, label>> planeFaces;
    for (label f = 0; f < mesh_.nFaces(); ++f)
    {
        const Vector& fv = sf[static_cast<std::size_t>(f)];
        const scalar a = mag(fv);
        if (a <= 0) continue;
        if (std::abs(fv[axis_]) < 0.5 * a) continue;
        planeFaces.emplace_back(fc[static_cast<std::size_t>(f)][axis_], f);
    }
    std::sort(planeFaces.begin(), planeFaces.end(),
              [](const std::pair<scalar, label>& x, const std::pair<scalar, label>& y) {
                  if (x.first != y.first) return x.first < y.first;
                  return x.second < y.second;
              });
    std::size_t k = 0;
    while (k < planeFaces.size())
    {
        const scalar anchor = planeFaces[k].first;
        std::size_t m = k + 1;
        while (m < planeFaces.size() && std::abs(planeFaces[m].first - anchor) <= band) ++m;
        std::vector<label> group;
        for (std::size_t j = k; j < m; ++j) group.push_back(planeFaces[j].second);
        planes_.push_back(std::move(group));
        k = m;
    }
}

void DetailBlock::applyBoundary(scalar qMag)
{
    const auto& patches = mesh_.patches();
    // the imposed velocity points INTO the region, opposite its outward area
    const Vector inletU = unitAxis(axis_) * (flowSign_ * qMag / inletArea_);

    for (label patch = 0; patch < mesh_.nPatches(); ++patch)
    {
        const std::string& type = patches[static_cast<std::size_t>(patch)].type;
        if (patch == inletPatch_)
        {
            U_.setPatchValue(patch, "fixedValue", inletU);
            p_.setPatchValue(patch, "zeroGradient", 0.0);
        }
        else if (patch == outletPatch_)
        {
            U_.setPatchValue(patch, "zeroGradient", Vector::zero());
            p_.setPatchValue(patch, "fixedValue", 0.0);
        }
        else if (type == "wall")
        {
            U_.setPatchValue(patch, "fixedValue", Vector::zero());
            p_.setPatchValue(patch, "zeroGradient", 0.0);
        }
        else
        {
            // symmetry planes and anything else: no statement about either
            // field. A `symmetry` condition would pin the normal velocity
            // component while leaving the tangential ones free, and SIMPLE
            // builds one momentum matrix whose diagonal the three components
            // share -- so the side patches take the same zeroGradient the
            // shipped channel case leaves them at. There is no driver for a
            // cross-stream velocity here, so the flux through these faces
            // stays zero, which is what a symmetry plane asks for.
            U_.setPatchValue(patch, "zeroGradient", Vector::zero());
            p_.setPatchValue(patch, "zeroGradient", 0.0);
        }
    }
    phi_.correctBoundaryFlux(U_);
}

scalar DetailBlock::lossAt(scalar qMag)
{
    if (!(qMag > 0))
        cfFatal("detailBlock: asked for the drop of a detailed region at a"
                " non-positive rate (" + std::to_string(qMag)
                + " m3/s): a region's resistance carries no information at"
                  " zero flow");

    applyBoundary(qMag);
    const SimpleResult res = solveIncompressibleSimple(mesh_, U_, p_, phi_, nu_,
                                                       schemes_, ctl_);
    measure(res, qMag);
    if (!std::isfinite(measure_.dpPa))
        cfFatal("detailBlock: non-finite pressure drop");
    if (measure_.dpPa <= 0)
        cfFatal("detailBlock: the detailed region reports a non-positive drop ("
                + std::to_string(measure_.dpPa) + " Pa) at "
                + std::to_string(qMag) + " m3/s");
    return measure_.dpPa;
}

void DetailBlock::measure(const SimpleResult& res, scalar qMag)
{
    const auto& volumes = mesh_.cellVolumes();
    const auto& pIn = p_.internalField();
    const auto& uIn = U_.internalField();
    const auto& phiV = phi_.values();
    const auto& patches = mesh_.patches();

    measure_ = BlockMeasure{};
    measure_.cells = mesh_.nCells();
    measure_.sweeps = res.iterations;
    measure_.converged = res.converged;
    measure_.imposedQ = qMag;
    measure_.inletArea = inletArea_;
    measure_.plugSpeed = qMag / inletArea_;

    if (sections_.size() < 2)
        cfFatal("detailBlock: the detailed region has fewer than two"
                " cross-sections to measure a drop between");

    // Volume-weighted mean pressure per cross-section, converted to Pa: the
    // kernel's p is kinematic (p/rho), the network's is not, and the interface
    // has to speak Pa.
    std::vector<scalar> pPa(sections_.size(), 0.0);
    for (std::size_t k = 0; k < sections_.size(); ++k)
    {
        scalar num = 0, den = 0;
        for (const label c : sections_[k])
        {
            const std::size_t u = static_cast<std::size_t>(c);
            num += pIn[u] * volumes[u];
            den += volumes[u];
        }
        pPa[k] = rho_ * (den > 0 ? num / den : 0.0);
    }
    measure_.dpPa = pPa.front() - pPa.back();

    // Boundary phi counts positive out of the domain, so the inlet reads
    // negative: -fin is the rate entering, fout the rate leaving.
    scalar fin = 0, fout = 0;
    const polyMesh::Patch& pin = patches[static_cast<std::size_t>(inletPatch_)];
    const polyMesh::Patch& pout = patches[static_cast<std::size_t>(outletPatch_)];
    for (label f = pin.start; f < pin.start + pin.size; ++f)
        fin += phiV[static_cast<std::size_t>(f)];
    for (label f = pout.start; f < pout.start + pout.size; ++f)
        fout += phiV[static_cast<std::size_t>(f)];
    measure_.fluxIn = -fin;
    measure_.fluxOut = fout;

    scalar planeMin = std::numeric_limits<scalar>::max();
    scalar planeMax = -std::numeric_limits<scalar>::max();
    for (const auto& plane : planes_)
    {
        scalar s = 0;
        for (const label f : plane) s += std::abs(phiV[static_cast<std::size_t>(f)]);
        planeMin = std::min(planeMin, s);
        planeMax = std::max(planeMax, s);
    }
    if (planes_.empty())
        cfFatal("detailBlock: no cross-section plane of axis-normal faces");
    measure_.planeFluxMin = planeMin;
    measure_.planeFluxMax = planeMax;

    // Developed-region gradient: a straight line through the section means. The
    // spread around it says whether the region is developed at all, which is
    // what the regression then compares against a developed-flow law.
    std::vector<scalar> xs, ys;
    for (std::size_t k = 0; k < sections_.size(); ++k)
    {
        if (sectionX_[k] >= probeMin_ && sectionX_[k] <= probeMax_)
        {
            xs.push_back(sectionX_[k]);
            ys.push_back(pPa[k]);
        }
    }
    measure_.probeSections = static_cast<label>(xs.size());
    if (xs.size() >= 2)
    {
        const scalar n = static_cast<scalar>(xs.size());
        scalar mx = 0, my = 0;
        for (std::size_t k = 0; k < xs.size(); ++k) { mx += xs[k]; my += ys[k]; }
        mx /= n; my /= n;
        scalar sxy = 0, sxx = 0;
        for (std::size_t k = 0; k < xs.size(); ++k)
        {
            sxy += (xs[k] - mx) * (ys[k] - my);
            sxx += (xs[k] - mx) * (xs[k] - mx);
        }
        measure_.probeSlope = sxx > 0 ? sxy / sxx : 0;
        scalar spread = 0;
        for (std::size_t k = 0; k < xs.size(); ++k)
            spread = std::max(spread, std::abs(ys[k] - my - measure_.probeSlope * (xs[k] - mx)));
        measure_.probeSpread = spread;
    }

    // Streamwise speed, signed against the flow direction so a reversed cell
    // reads as a recirculation instead of hiding inside a magnitude.
    std::vector<char> inProbe(static_cast<std::size_t>(mesh_.nCells()), 0);
    for (std::size_t k = 0; k < sections_.size(); ++k)
    {
        const char hit = (sectionX_[k] >= probeMin_ && sectionX_[k] <= probeMax_) ? 1 : 0;
        for (const label c : sections_[k]) inProbe[static_cast<std::size_t>(c)] = hit;
    }
    const Vector axisUnit = unitAxis(axis_) * static_cast<scalar>(flowSign_);
    scalar uxMax = 0, uxCentre = 0;
    for (label c = 0; c < mesh_.nCells(); ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const scalar ux = dot(uIn[u], axisUnit);
        uxMax = std::max(uxMax, ux);
        if (inProbe[u]) uxCentre = std::max(uxCentre, ux);
    }
    measure_.uxMax = uxMax;
    measure_.uxCentre = uxCentre;
}

// ------------------------------------------------------------------
// the partitioned outer loop
// ------------------------------------------------------------------

CoupledResult solveCoupled1d3d(OneDNetwork& net,
                               const std::vector<label>& edges,
                               std::vector<ICoupledElement*> blocks,
                               const CoupledControls& ctl,
                               const SteadyControls& netCtl)
{
    if (edges.empty())
        cfFatal("coupling: the interface needs at least one coupled edge");
    if (blocks.size() != edges.size())
        cfFatal("coupling: " + std::to_string(edges.size()) + " coupled edges need"
                " one partner each, got " + std::to_string(blocks.size()));
    if (!(ctl.relax > 0 && ctl.relax <= 1))
        cfFatal("coupling: the interface relaxation must be in (0, 1]");

    // every edge is a pipe the network can charge a coefficient for, and one
    // pipe cannot carry two partners: they would overwrite each other's drop
    std::vector<NetworkPipe*> coupled;
    for (std::size_t i = 0; i < edges.size(); ++i)
    {
        const label edge = edges[i];
        if (edge < 0 || static_cast<std::size_t>(edge) >= net.pipes.size())
            cfFatal("coupling: no pipe at index " + std::to_string(edge));
        bool repeats = false;
        for (std::size_t j = 0; j < i; ++j)
            repeats = repeats || (edges[j] == edge);
        NetworkPipe& pipe = net.pipes[static_cast<std::size_t>(edge)];
        if (repeats)
            cfFatal("coupling: pipe '" + pipe.name + "' is coupled twice");
        if (blocks[i] == nullptr)
            cfFatal("coupling: edge '" + pipe.name + "' has no partner");
        if (!(pipe.coupledResistance > 0))
            cfFatal("coupling: pipe '" + pipe.name + "' is a coupled edge but"
                    " carries no starting resistance");
        coupled.push_back(&pipe);
    }

    CoupledResult out;
    out.edges = edges;

    for (label it = 1; it <= ctl.maxOuterIterations; ++it)
    {
        SteadyNetworkResult res = solveSteadyNetwork(net, netCtl);
        if (!res.converged)
            cfFatal("coupling: the 1D network did not converge at interface"
                    " iteration " + std::to_string(it) + " (residual "
                    + fmt(res.relativeResidual) + " after "
                    + std::to_string(res.iterations) + " iterations)");

        OuterIteration step;
        step.iteration = it;
        step.networkResidual = res.relativeResidual;
        scalar worst = 0;

        for (std::size_t i = 0; i < coupled.size(); ++i)
        {
            const NetworkPipe& pipe = *coupled[i];
            const label edge = edges[i];
            const PipeFlowResult& edgeFlow =
                res.pipeFlows[static_cast<std::size_t>(edge)];
            const scalar q = std::abs(edgeFlow.flowRate);
            if (!(q > 0))
                cfFatal("coupling: the network delivers no rate through coupled"
                        " pipe '" + pipe.name + "'; a partner cannot be evaluated"
                        " at zero rate");
            const scalar dpBlock = blocks[i]->lossAt(q);

            EdgeExchange ex;
            ex.edge = edge;
            ex.name = pipe.name;
            ex.resistance = pipe.coupledResistance;
            ex.q = q;
            ex.dpPa = dpBlock;
            // The network's own statement about this edge against the partner's.
            // Both come out of a solver, so this gap is a cross-solver residual and
            // not a number compared with itself.
            ex.residual = std::abs(edgeFlow.dpTotal - dpBlock)
                          / std::max(std::abs(dpBlock), 1e-300);
            ex.sweeps = blocks[i]->lastSweeps();
            ex.blockConverged = blocks[i]->lastConverged();
            worst = std::max(worst, ex.residual);
            step.edges.push_back(ex);
        }

        out.history.push_back(step);
        // this pass is a state of the coupled system: the network was solved with
        // these coefficients, so each rate is the rate that coefficient delivers,
        // and every partner was measured at exactly its own rate. Record it as the
        // answer before touching the coefficients, so a caller that reads
        // (R, q, dpNetwork, dpBlock) per edge reads one operating point and not
        // two solves apart.
        out.network = res;
        out.state = step.edges;
        out.residual = worst;
        out.iterations = it;
        out.converged = worst <= ctl.tolerance;
        if (out.converged) break;

        // the next pass starts from the refreshed coefficients; the accepted
        // state above keeps the ones it was solved with
        for (std::size_t i = 0; i < coupled.size(); ++i)
        {
            const EdgeExchange& ex = step.edges[i];
            coupled[i]->coupledResistance =
                (1.0 - ctl.relax) * ex.resistance + ctl.relax * (ex.dpPa / ex.q);
        }
    }

    return out;
}

// ------------------------------------------------------------------
// the transient partitioned loop
// ------------------------------------------------------------------

TransientCoupledResult solveTransientCoupled1d3d(TransientNetwork& net,
                                                 const std::vector<label>& edges,
                                                 std::vector<ICoupledElement*> blocks,
                                                 const TransientCoupledControls& ctl,
                                                 double dt, double endTime)
{
    if (edges.empty())
        cfFatal("coupling: the transient interface needs at least one coupled edge");
    if (blocks.size() != edges.size())
        cfFatal("coupling: " + std::to_string(edges.size()) + " coupled edges need"
                " one partner each, got " + std::to_string(blocks.size()));
    if (!(ctl.relax > 0 && ctl.relax <= 1))
        cfFatal("coupling: the interface relaxation must be in (0, 1]");
    if (ctl.maxSubIterations < 1)
        cfFatal("coupling: the transient interface needs at least one pass per"
                " time step");
    if (!(dt > 0))
        cfFatal("coupling: the transient step size must be positive (dt "
                + fmt(dt) + ")");
    if (!(endTime > 0))
        cfFatal("coupling: the transient run needs a positive endTime ("
                + fmt(endTime) + ")");
    const std::size_t nSteps = static_cast<std::size_t>(endTime / dt + 0.5);
    if (nSteps < 1)
        cfFatal("coupling: endTime " + fmt(endTime) + " is shorter than one"
                " step of " + fmt(dt));

    std::vector<scalar> R(edges.size());
    for (std::size_t i = 0; i < edges.size(); ++i)
    {
        const label edge = edges[i];
        if (edge < 0 || static_cast<std::size_t>(edge) >= net.pipeNames().size())
            cfFatal("coupling: no transient pipe at index " + std::to_string(edge));
        const std::string& name = net.pipeNames()[static_cast<std::size_t>(edge)];
        bool repeats = false;
        for (std::size_t j = 0; j < i; ++j)
            repeats = repeats || (edges[j] == edge);
        if (repeats)
            cfFatal("coupling: transient pipe '" + name + "' is coupled twice");
        if (blocks[i] == nullptr)
            cfFatal("coupling: transient edge '" + name + "' has no partner");
        R[i] = net.coupledResistanceOf(edge);
        if (!(R[i] > 0))
            cfFatal("coupling: transient pipe '" + name
                    + "' is a coupled edge but carries no starting resistance");
    }

    TransientCoupledResult out;
    out.edges = edges;
    out.nodeNames.assign(net.nodeNames().begin(), net.nodeNames().end());
    out.pipeNames.assign(net.pipeNames().begin(), net.pipeNames().end());

    for (std::size_t s = 0; s < nSteps; ++s)
    {
        TransientCoupledStep step;
        label passes = 0;
        bool agreed = false;

        // Re-solve the same step from the same committed start until every
        // partner's drop belongs to the rate the step delivers through its edge.
        // Nothing is advanced in here, so the refresh is a real iteration on one
        // time level rather than a drift along the trajectory.
        for (label pass = 1; pass <= ctl.maxSubIterations; ++pass)
        {
            const double newton = net.solveStep(dt);

            TransientInterfacePass rec;
            rec.pass = pass;
            rec.newtonResidual = newton;
            scalar worst = 0;
            for (std::size_t i = 0; i < edges.size(); ++i)
            {
                const label edge = edges[i];
                EdgeExchange ex;
                ex.edge = edge;
                ex.name = net.pipeNames()[static_cast<std::size_t>(edge)];
                ex.resistance = R[i];
                ex.q = std::abs(net.pendingPipeFlow(edge));
                if (!(ex.q > 0))
                    cfFatal("coupling: the transient step delivers no rate through"
                            " coupled pipe '" + ex.name + "'; a partner cannot be"
                            " evaluated at zero rate");
                ex.dpPa = blocks[i]->lossAt(ex.q);
                ex.residual = std::abs(ex.resistance * ex.q - ex.dpPa)
                              / std::max(std::abs(ex.dpPa), 1e-300);
                ex.sweeps = blocks[i]->lastSweeps();
                ex.blockConverged = blocks[i]->lastConverged();
                worst = std::max(worst, ex.residual);
                rec.edges.push_back(ex);
            }
            step.history.push_back(rec);
            passes = pass;

            if (worst <= ctl.tolerance)
            {
                agreed = true;
                break;
            }
            for (std::size_t i = 0; i < edges.size(); ++i)
            {
                const EdgeExchange& ex = rec.edges[i];
                R[i] = (1.0 - ctl.relax) * R[i] + ctl.relax * (ex.dpPa / ex.q);
                net.setCoupledResistance(edges[i], R[i]);
            }
        }

        if (!agreed)
        {
            const EdgeExchange& worst = *std::max_element(
                step.history.back().edges.begin(),
                step.history.back().edges.end(),
                [](const EdgeExchange& a, const EdgeExchange& b)
                { return a.residual < b.residual; });
            cfFatal("coupling: the transient interface did not agree at step "
                    + std::to_string(s + 1) + " after "
                    + std::to_string(passes) + " passes (worst edge '" + worst.name
                    + "': residual " + fmt(worst.residual) + ", rate "
                    + fmt(worst.q) + " m3/s at t " + fmt(net.time() + dt) + " s)");
        }

        // the accepted step is the one solved with the coefficients it reports,
        // and those coefficients are what the next step starts from
        net.commitStep();
        step.t = net.time();
        step.passes = passes;
        step.edges = step.history.back().edges;
        step.pressure = net.pressure();
        step.flowRate = net.pipeFlows();
        out.steps.push_back(std::move(step));
    }

    out.state = out.steps.back().edges;
    for (const EdgeExchange& ex : out.state)
        out.residual = std::max(out.residual, ex.residual);
    return out;
}

} // namespace coupling
} // namespace cf
