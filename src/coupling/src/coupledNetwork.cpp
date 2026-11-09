#include "cf/coupling/CoupledNetwork.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
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

scalar dot(const std::vector<scalar>& a, const std::vector<scalar>& b)
{
    scalar s = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
        s += a[i] * b[i];
    return s;
}

scalar norm(const std::vector<scalar>& a) { return std::sqrt(dot(a, a)); }

void axpy(std::vector<scalar>& y, const std::vector<scalar>& x, scalar a)
{
    for (std::size_t i = 0; i < y.size(); ++i)
        y[i] += a * x[i];
}

/// Interface quasi-Newton over the whole coefficient vector (Aitken / IQN-ILS).
///
/// The interface fixed point is `r(x) = g(x) - x = 0`, where `x` is the vector of
/// coefficients the network carries and `g(x)` what the partners report for the
/// rates that `x` delivers. Each pass supplies one secant pair, so the stored
/// pairs define an approximate inverse of `r'` on the subspace their `dr`
/// directions span; the step is Newton there and the plain relaxed step along the
/// directions nothing has been measured yet.
///
/// One edge and one stored pair is exactly Aitken's delta-squared step: the
/// residual's component along `dr` is projected out, so only the secant term
/// `dr·r` survives. That is the case the un-accelerated loop cannot do - a
/// strongly nonlinear partner makes the map's derivative far above 1, and plain
/// Gauss-Seidel then diverges however many passes it is given.
///
/// Vectors are normalised by each edge's starting coefficient, so the least
/// squares is not decided by which edge happens to carry the larger resistance.
class InterfaceQuasiNewton
{
public:
    InterfaceQuasiNewton(label history, scalar limit)
        : maxPairs_(history < 1 ? 1 : static_cast<std::size_t>(history)), limit_(limit)
    {}

    InterfaceQuasiNewton(const InterfaceQuasiNewton&) = delete;
    InterfaceQuasiNewton& operator=(const InterfaceQuasiNewton&) = delete;

    /// `x` is what `r` was measured at; returns the next `x`.
    std::vector<scalar> next(const std::vector<scalar>& x,
                             const std::vector<scalar>& r,
                             scalar relax)
    {
        if (hasPrev_)
        {
            std::vector<scalar> dx(x.size()), dr(r.size());
            for (std::size_t i = 0; i < dx.size(); ++i)
            {
                dx[i] = x[i] - xPrev_[i];
                dr[i] = r[i] - rPrev_[i];
            }
            note(std::move(dx), std::move(dr));
        }
        xPrev_ = x;
        rPrev_ = r;
        hasPrev_ = true;
        buildBasis(x.size());

        std::vector<scalar> perp = r;
        for (const std::vector<scalar>& q : pairsQ_)
            axpy(perp, q, -dot(q, r));

        std::vector<scalar> step;
        for (const scalar v : perp)
            step.push_back(relax * v);
        for (std::size_t j = 0; j < pairsQ_.size(); ++j)
            axpy(step, pairsX_[j], -dot(pairsQ_[j], r));

        const scalar rn = norm(r);
        const scalar sn = norm(step);
        if (limit_ > 0 && rn > 0 && sn > limit_ * rn)
        {
            const scalar cap = (limit_ * rn) / sn;
            for (scalar& v : step)
                v *= cap;
        }

        std::vector<scalar> out = x;
        axpy(out, step, 1.0);
        // a coefficient at or below zero would charge the network a drop that
        // opposes its own rate; a secant pair built from a far-from-linear map can
        // ask for that, so the step is floored at a millionth of the seed instead.
        for (scalar& v : out)
            if (!(v > 1e-6)) v = 1e-6;
        return out;
    }

    /// The map changed, so the pairs no longer measure it: a new time step solves
    /// a different interface problem even when the partners are unchanged.
    void reset()
    {
        windowX_.clear();
        windowR_.clear();
        pairsX_.clear();
        pairsQ_.clear();
        hasPrev_ = false;
    }

private:
    /// Record the pair just measured.
    void note(std::vector<scalar> dx, std::vector<scalar> dr)
    {
        if (!(norm(dr) > 0))
            return;
        windowX_.push_back(std::move(dx));
        windowR_.push_back(std::move(dr));
        while (windowX_.size() > maxPairs_)
        {
            windowX_.erase(windowX_.begin());
            windowR_.erase(windowR_.begin());
        }
    }

    /// Re-orthogonalise the window from the newest pair downwards. The newest pair
    /// is the best local description of the map, so it always keeps a column; an
    /// older one survives only for the directions the newer pair cannot see.
    /// Keeping the *oldest* instead is what freezes a single-edge interface: all
    /// its pairs are parallel, so only the first would ever be stored and every
    /// later step would be Newton on one stale measurement rather than a secant
    /// on the two most recent ones.
    void buildBasis(std::size_t size)
    {
        pairsX_.clear();
        pairsQ_.clear();
        for (std::size_t p = windowR_.size(); p-- > 0;)
        {
            std::vector<scalar> q = windowR_[p];
            std::vector<scalar> w = windowX_[p];
            const scalar raw = norm(q);
            // the same combinations on both vectors keep the secant condition
            // holding, and the dr columns then span an orthonormal basis, so the
            // least squares is a projection and needs no matrix solve.
            for (std::size_t j = 0; j < pairsQ_.size(); ++j)
            {
                const scalar s = dot(pairsQ_[j], q);
                axpy(q, pairsQ_[j], -s);
                axpy(w, pairsX_[j], -s);
            }
            const scalar kept = norm(q);
            if (!(kept > 1e-6 * raw))
                continue;   // adds no direction the newer pairs do not already have
            for (scalar& v : q) v /= kept;
            for (scalar& v : w) v /= kept;
            pairsQ_.push_back(std::move(q));
            pairsX_.push_back(std::move(w));
            if (pairsQ_.size() == size)
                break;   // the interface operator cannot have more rank than this
        }
    }

    std::vector<std::vector<scalar>> windowX_;
    std::vector<std::vector<scalar>> windowR_;
    std::vector<std::vector<scalar>> pairsX_;
    std::vector<std::vector<scalar>> pairsQ_;
    std::vector<scalar> xPrev_;
    std::vector<scalar> rPrev_;
    bool hasPrev_ = false;
    std::size_t maxPairs_;
    scalar limit_;
};

} // namespace

DetailBlock::DetailBlock(polyMesh mesh, scalar nu, scalar rho,
                         const Schemes& schemes, const SimpleControls& ctl,
                         scalar probeMin, scalar probeMax,
                         const TurbulenceConfig& turb,
                         const std::string& inletShape)
    : mesh_(std::move(mesh)),
      nu_(nu),
      rho_(rho),
      schemes_(schemes),
      ctl_(ctl),
      probeMin_(probeMin),
      probeMax_(probeMax),
      turb_(turb),
      inletShape_(inletShape),
      U_(mesh_, "U"),
      p_(mesh_, "p"),
      phi_(mesh_, "phi", 0.0)
{
    if (!(nu_ > 0) || !(rho_ > 0))
        cfFatal("detailBlock: nu and rho must be positive");
    turb_.validate();

    // Whether this region carries its own acceleration is the case's statement,
    // not the interface's: `schemes{ ddt Euler; }` in the block's own dict. Only
    // the two words the kernel implements are accepted, because a block that
    // quietly ignored a third would report a quasi-steady drop and let the
    // interface charge it to a momentum row it never solved.
    if (schemes_.ddt == "steady") unsteady_ = false;
    else if (schemes_.ddt == "Euler") unsteady_ = true;
    else
        cfFatal("detailBlock: schemes ddt is '" + schemes_.ddt
                + "'; a detailed region understands 'steady' and 'Euler'");

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
    buildInletProfile();
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

    // The region's own hydraulic inertance, for the interface convention that
    // lets the block state its region's momentum: rho*V/A^2. It is a property of
    // this mesh, not of the proxy pipe the 1D side declares, so the two can be
    // compared rather than assumed equal. The formula needs one section, and a
    // region whose planes differ in area has no single inertance at all, so the
    // areas are checked instead of trusted.
    scalar areaRef = 0.0;
    bool firstPlane = true;
    for (const std::vector<label>& plane : planes_)
    {
        scalar a = 0.0;
        for (const label f : plane) a += std::abs(sf[static_cast<std::size_t>(f)][axis_]);
        if (firstPlane)
        {
            areaRef = a;
            firstPlane = false;
        }
        else if (std::abs(a - areaRef) > 1e-9 * areaRef)
            cfFatal("detailBlock: the detailed region's cross-section is not"
                    " constant (plane areas " + std::to_string(areaRef) + " vs "
                    + std::to_string(a) + " m2), so its momentum has no single"
                      " inertance for the interface to carry");
    }
    if (!(areaRef > 0) || !(areaRef < std::numeric_limits<scalar>::max()))
        cfFatal("detailBlock: the detailed region's cross-section area is not a"
                " positive, finite number (" + std::to_string(areaRef)
                + " m2), so rho*V/A^2 has no value for the interface to carry");
    scalar volume = 0.0;
    for (label c = 0; c < mesh_.nCells(); ++c)
        volume += mesh_.cellVolumes()[static_cast<std::size_t>(c)];
    // A region with no section to refer its momentum to -- or one whose cell
    // volumes do not add up to a finite volume -- has no inertance to hand the
    // interface, so the block says so on its own geometry instead of letting a
    // NaN reach a momentum row.
    inertance_ = rho_ * volume / (areaRef * areaRef);
    if (!(areaRef > 0) || !(volume > 0) || !std::isfinite(inertance_))
        cfFatal("detailBlock: the detailed region measures"
                " " + std::to_string(volume) + " m3 of fluid across a"
                " " + std::to_string(areaRef) + " m2 section, which is no"
                  " rho*V/A^2 for the interface to install");

    // The same number over the slab the drop is actually measured across. The
    // probe range is a statement about which part of the region the interface is
    // charging the network for, and the momentum that belongs to that charge is
    // the slab's, not the whole duct's: the two differ by exactly the length the
    // probe leaves out. A block that steps in time can say which of them its own
    // acceleration agrees with, so the difference becomes a measurement instead
    // of a convention.
    scalar probeVolume = 0.0;
    for (std::size_t k = 0; k < sections_.size(); ++k)
    {
        if (sectionX_[k] < probeMin_ || sectionX_[k] > probeMax_) continue;
        for (const label c : sections_[k])
            probeVolume += mesh_.cellVolumes()[static_cast<std::size_t>(c)];
    }
    inertanceProbe_ = rho_ * probeVolume / (areaRef * areaRef);
    if (!(probeVolume > 0) || !std::isfinite(inertanceProbe_))
        cfFatal("detailBlock: the probe range [" + std::to_string(probeMin_) + ", "
                + std::to_string(probeMax_) + "] m covers no fluid of this region"
                  " (" + std::to_string(probeVolume) + " m3), so the drop measured"
                    " across it has no momentum of its own");
}

void DetailBlock::buildInletProfile()
{
    // "flat" is the condition this coupling has always written: one vector over
    // the whole patch. It stays that way, so a case that names no shape gets the
    // numbers it got before any of this existed.
    if (inletShape_ == "flat") return;
    if (inletShape_ != "parabolic" && inletShape_ != "power1over7")
        cfFatal("detailBlock: the inlet profile '" + inletShape_ + "' is not one"
                " this block can write (flat, parabolic, power1over7)");

    // Which transverse directions this section is confined by no-slip walls. A
    // boundary layer grows off a wall, so a profile is a statement about the
    // direction the walls bound; varying across a symmetry side would ask the
    // mesh for shear where that condition forbids any.
    const auto& patches = mesh_.patches();
    int bounds[3][2] = {{0, 0}, {0, 0}, {0, 0}};
    for (label patch = 0; patch < mesh_.nPatches(); ++patch)
    {
        const polyMesh::Patch& pd = patches[static_cast<std::size_t>(patch)];
        if (pd.type != "wall") continue;
        const Vector sf = patchAreaVector(mesh_, patch);
        if (!(mag(sf) > 0)) continue;
        int k = 0;
        for (int c = 1; c < 3; ++c)
            if (std::abs(sf[c]) > std::abs(sf[k])) k = c;
        if (k == axis_) continue;
        bounds[k][sf[k] < 0 ? 0 : 1] = 1;
    }
    std::vector<int> axes;
    for (int k = 0; k < 3; ++k)
        if (k != axis_ && bounds[k][0] && bounds[k][1]) axes.push_back(k);
    if (axes.empty())
        cfFatal("detailBlock: the inlet profile '" + inletShape_ + "' has no"
                " wall-confined transverse direction to vary along, so no shape can"
                " be written on this section");

    // The shape evaluated at the inlet face centres, scaled so that its section
    // integral is exactly the section: the rate the network delivered is then
    // carried to the last digit and only its distribution over the section
    // changes. Nothing here knows what the developed profile should be -- that is
    // what the solve then has to say.
    const polyMesh::Patch& pin = patches[static_cast<std::size_t>(inletPatch_)];
    const auto& fc = mesh_.faceCentres();
    const auto& sfAll = mesh_.faceSf();
    // The section's bounds are the walls' bounds, not the outermost face
    // centres' -- a boundary face always sits inside the wall it belongs to, so
    // taking the face centres as the extent would put eta at +-1 on the first
    // and last face and write a zero velocity there.
    const Vector loBox = mesh_.minX();
    const Vector hiBox = mesh_.maxX();
    scalar lo[3], hi[3];
    for (int k = 0; k < 3; ++k)
    {
        lo[k] = loBox[k];
        hi[k] = hiBox[k];
    }
    scalar wMin = std::numeric_limits<scalar>::max();
    for (const int k : axes)
    {
        if (!(hi[k] > lo[k]))
            cfFatal("detailBlock: the inlet section has no extent along axis "
                    + std::to_string(k) + " to spread a profile over");
    }
    std::vector<scalar> w(static_cast<std::size_t>(pin.size), 0.0);
    for (std::size_t j = 0; j < static_cast<std::size_t>(pin.size); ++j)
    {
        const std::size_t u = static_cast<std::size_t>(pin.start) + j;
        scalar value = 1.0;
        for (const int k : axes)
        {
            const scalar mid = 0.5 * (hi[k] + lo[k]);
            const scalar half = 0.5 * (hi[k] - lo[k]);
            const scalar eta = (fc[u][k] - mid) / half;
            value *= inletShape_ == "parabolic"
                         ? (1.0 - eta * eta)
                         : std::pow(1.0 - std::abs(eta), 1.0 / 7.0);
        }
        w[j] = value;
        wMin = std::min(wMin, value);
    }
    if (!(wMin > 0))
        cfFatal("detailBlock: the inlet profile '" + inletShape_ + "' writes a"
                " zero or reversed weight on a face of this section, so the section"
                " cannot carry the rate at all");

    inletFaces_.clear();
    inletWeights_.clear();
    scalar sumWA = 0, sumA = 0;
    for (std::size_t j = 0; j < static_cast<std::size_t>(pin.size); ++j)
    {
        const std::size_t u = static_cast<std::size_t>(pin.start) + j;
        const scalar a = mag(sfAll[u]);
        sumWA += w[j] * a;
        sumA += a;
    }
    if (!(sumA > 0) || !(sumWA > 0))
        cfFatal("detailBlock: the inlet section integrates to no area or no shape");
    if (std::abs(sumA - inletArea_) > 1e-9 * inletArea_)
        cfFatal("detailBlock: the inlet patch's face areas add to "
                + std::to_string(sumA) + " m2 against its own area vector of "
                + std::to_string(inletArea_) + " m2, so the section cannot be"
                  " normalised");
    const scalar norm = sumWA / sumA;
    for (std::size_t j = 0; j < static_cast<std::size_t>(pin.size); ++j)
    {
        inletFaces_.push_back(pin.start + static_cast<label>(j));
        inletWeights_.push_back(w[j] / norm);
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
            // A shape replaces the value face by face, never the condition: the
            // patch still prescribes velocity, and `phi` below still comes from
            // these faces, so the rate entering the region is the one the
            // interface asked for whichever way the section is filled.
            if (!inletWeights_.empty())
            {
                for (std::size_t j = 0; j < inletFaces_.size(); ++j)
                    U_.setFaceValue(inletFaces_[j], inletU * inletWeights_[j]);
            }
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

void DetailBlock::beginTimeLevel(double dt, scalar qPrev)
{
    if (!unsteady_) return;   // the region has no history to take
    if (!(dt > 0))
        cfFatal("detailBlock: a time level of " + std::to_string(dt)
                + " s was declared for a region that steps in time");
    levelDt_ = dt;
    levelQPrev_ = qPrev;
    // THE level's history, taken once: every pass of this level - and every
    // sweep inside a pass - accelerates away from this field, so re-solving the
    // level is an iteration on one time step rather than a step along the
    // trajectory. The field itself keeps warm-starting from the previous pass,
    // which is what makes the level cheap to re-solve.
    UOldLevel_ = U_.internalField();
}

scalar DetailBlock::lossAt(scalar qMag)
{
    if (!(qMag > 0))
        cfFatal("detailBlock: asked for the drop of a detailed region at a"
                " non-positive rate (" + std::to_string(qMag)
                + " m3/s): a region's resistance carries no information at"
                  " zero flow");

    if (unsteady_)
    {
        if (!(levelDt_ > 0))
            cfFatal("detailBlock: this region is asked to carry its own du/dt but"
                    " was solved at an instant: no time level was declared, so"
                    " there is no field for it to accelerate away from");
        ctl_.dt = levelDt_;
        ctl_.uTimeOld = &UOldLevel_;
    }

    applyBoundary(qMag);
    // Built on the rate the interface first delivers, because the model's inlet
    // turbulence is scaled by the imposed plug speed: `applyBoundary` has just
    // put that speed on the inlet patch, and the factory reads it from there.
    // Later passes keep the inlet k/omega at that scale and warm-start the rest
    // of the field, which is what the wall layer is solved for anyway.
    if (!turbModel_ && turb_.model != "laminar")
        turbModel_ = RansModel::create(mesh_, U_, phi_, nu_, schemes_, turb_);

    const SimpleResult res = solveIncompressibleSimple(mesh_, U_, p_, phi_, nu_,
                                                       schemes_, ctl_,
                                                       turbModel_.get());
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
    measure_.inertance = inertance_;
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

    // What the interface handed the inlet, as the region received it: the
    // section integral of the condition is the rate, and the shape factor of the
    // same faces says how that rate sits across the section. Both are read off
    // the boundary condition itself, so they grade the mapping rather than the
    // flow -- and the flux identity is the statement that a shape was traded for
    // a shape and not for a different amount of fluid.
    {
        const polyMesh::Patch& pin = patches[static_cast<std::size_t>(inletPatch_)];
        const auto& sfIn = mesh_.faceSf();
        scalar sumA = 0, sumU = 0, sumU2 = 0;
        for (label f = pin.start; f < pin.start + pin.size; ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            const scalar a = mag(sfIn[u]);
            const scalar ux = dot(U_.faceValue(f), axisUnit);
            sumA += a;
            sumU += ux * a;
            sumU2 += ux * ux * a;
        }
        measure_.fluxImposed = sumU;
        measure_.betaImposed = std::abs(sumU) > 0
                                   ? (sumU2 * sumA) / (sumU * sumU) : 1.0;
        measure_.momentumFluxImposedPa =
            rho_ * measure_.betaImposed * sumU * sumU / sumA;
        measure_.inletShape = inletShape_;
    }

    // The same factor, measured on the field the region actually solved, station
    // by station: the profile the inlet arrived with is not the profile that
    // leaves, and the distance over which one becomes the other is what a single
    // imposed rate could never say.
    {
        const std::size_t nSec = sections_.size();
        const std::size_t stride = std::max<std::size_t>(1, (nSec + 15) / 16 - 1);
        for (std::size_t k = 0; k < nSec; k += stride)
        {
            scalar v = 0, u1 = 0, u2 = 0;
            for (const label c : sections_[k])
            {
                const std::size_t t = static_cast<std::size_t>(c);
                const scalar ux = dot(uIn[t], axisUnit);
                v += volumes[t];
                u1 += ux * volumes[t];
                u2 += ux * ux * volumes[t];
            }
            const scalar beta = (v > 0 && std::abs(u1) > 0)
                                    ? (u2 * v) / (u1 * u1) : 0.0;
            measure_.betaStations.emplace_back(sectionX_[k], beta);
            measure_.pressureStations.emplace_back(sectionX_[k], pPa[k]);
        }
        if ((nSec - 1) % stride != 0)
        {
            scalar v = 0, u1 = 0, u2 = 0;
            for (const label c : sections_[nSec - 1])
            {
                const std::size_t t = static_cast<std::size_t>(c);
                const scalar ux = dot(uIn[t], axisUnit);
                v += volumes[t];
                u1 += ux * volumes[t];
                u2 += ux * ux * volumes[t];
            }
            measure_.betaStations.emplace_back(
                sectionX_[nSec - 1],
                (v > 0 && std::abs(u1) > 0) ? (u2 * v) / (u1 * u1) : 0.0);
            measure_.pressureStations.emplace_back(sectionX_[nSec - 1],
                                                   pPa[nSec - 1]);
        }
        if (measure_.betaStations.size() < 2)
            cfFatal("detailBlock: the region solved with fewer than two"
                    " cross-sections to measure a shape at");
        measure_.betaIn = measure_.betaStations.front().second;
        measure_.betaOut = measure_.betaStations.back().second;
    }

    if (turbModel_)
    {
        scalar nutMax = 0;
        for (const scalar v : turbModel_->nut()) nutMax = std::max(nutMax, v);
        measure_.turbModel = turbModel_->name();
        measure_.nutMax = nutMax;
        measure_.wallYPlusMax = turbModel_->maxWallYPlus();
        measure_.wallUTau = turbModel_->meanFrictionVelocity();
    }

    // The region's own acceleration, measured on the slab the drop above is
    // taken across: `rho * int du/dt dV / A` between the first and the last
    // probe section, from the field just solved and the field the previous
    // level left. Nothing is assumed about inertance here; what this force is
    // worth per unit rate-of-change is then a number the 1D side's momentum row
    // can be graded against instead of a number it declares.
    measure_.unsteady = unsteady_;
    measure_.inertanceProbe = inertanceProbe_;
    if (unsteady_)
    {
        measure_.dtLevel = levelDt_;
        measure_.dqPrev = qMag - levelQPrev_;
        scalar sumV = 0, sumU = 0, sumU2 = 0, sumVdu = 0;
        for (label c = 0; c < mesh_.nCells(); ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            if (!inProbe[u]) continue;
            const scalar v = volumes[u];
            const scalar ux = dot(uIn[u], axisUnit);
            const scalar uxo = dot(UOldLevel_[u], axisUnit);
            sumV += v;
            sumU += v * ux;
            sumU2 += v * ux * ux;
            sumVdu += v * (ux - uxo);
        }
        if (!(sumV > 0))
            cfFatal("detailBlock: the probe slab holds no fluid to accelerate");
        if (!(levelDt_ > 0))
            cfFatal("detailBlock: measured an acceleration against a level of "
                    + std::to_string(levelDt_) + " s");
        measure_.dpInertialPa = rho_ * sumVdu / (levelDt_ * inletArea_);
        measure_.dpFrictionPa = measure_.dpPa - measure_.dpInertialPa;
        const scalar accel = measure_.dqPrev / levelDt_;
        if (std::abs(accel) > 0)
            measure_.inertanceDynamic = measure_.dpInertialPa / accel;
        // <u_x^2> / <u_x>^2 over the same slab: the momentum-flux correction a
        // flat-rate picture cannot see, measured on the profile the block
        // actually solved. A plug profile reads 1.0.
        if (std::abs(sumU) > 0)
            measure_.betaProfile = (sumU2 * sumV) / (sumU * sumU);
    }
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

    // the coefficients the run starts from: the quasi-Newton works in units of
    // them, so the interface operator is scaled edge by edge
    std::vector<scalar> seed(coupled.size());
    for (std::size_t i = 0; i < coupled.size(); ++i)
        seed[i] = coupled[i]->coupledResistance;

    std::unique_ptr<InterfaceQuasiNewton> qn;
    if (ctl.accelerate)
        qn = std::make_unique<InterfaceQuasiNewton>(ctl.qnHistory, ctl.qnLimit);

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
        if (qn)
        {
            std::vector<scalar> x(coupled.size()), r(coupled.size());
            for (std::size_t i = 0; i < coupled.size(); ++i)
            {
                const EdgeExchange& ex = step.edges[i];
                x[i] = ex.resistance / seed[i];
                r[i] = (ex.dpPa / ex.q - ex.resistance) / seed[i];
            }
            const std::vector<scalar> next = qn->next(x, r, ctl.relax);
            for (std::size_t i = 0; i < coupled.size(); ++i)
                coupled[i]->coupledResistance = next[i] * seed[i];
        }
        else
        {
            for (std::size_t i = 0; i < coupled.size(); ++i)
            {
                const EdgeExchange& ex = step.edges[i];
                coupled[i]->coupledResistance =
                    (1.0 - ctl.relax) * ex.resistance + ctl.relax * (ex.dpPa / ex.q);
            }
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

    // the coefficients the run starts from: the quasi-Newton works in units of
    // them, edge by edge
    const std::vector<scalar> Rseed = R;

    TransientCoupledResult out;
    out.edges = edges;
    out.nodeNames.assign(net.nodeNames().begin(), net.nodeNames().end());    out.pipeNames.assign(net.pipeNames().begin(), net.pipeNames().end());

    std::unique_ptr<InterfaceQuasiNewton> qn;
    if (ctl.accelerate)
        qn = std::make_unique<InterfaceQuasiNewton>(ctl.qnHistory, ctl.qnLimit);

    // Where the detailed region's inertance comes from - per edge if the caller
    // says so, otherwise one answer for the whole interface. A partner that
    // measures its own hands the number to the network's momentum row, which
    // keeps it implicit exactly as the pipe's own rho*L/A was kept - the one way
    // a detailed region can state its momentum without putting it into the
    // interface coefficient, where a rate falling through the step would turn it
    // negative. The measurement is a property of the mesh, so it is taken once
    // before any step is solved, and every edge that is asked for is checked
    // before any of them is installed: a partner that states no region must not
    // leave the network half-declared, and an edge the caller did not name a
    // source for keeps the drawing's number rather than inheriting its
    // neighbour's.
    {
        const std::vector<int>& perEdge = ctl.inertiaFromPartnerPerEdge;
        if (!perEdge.empty() && perEdge.size() != edges.size())
            cfFatal("coupling: transient controls name an inertia source for "
                    + std::to_string(perEdge.size()) + " edges but the interface"
                      " has " + std::to_string(edges.size()));
        std::vector<int> source(edges.size(), ctl.inertiaFromPartner ? 1 : 0);
        for (std::size_t i = 0; i < perEdge.size(); ++i) source[i] = perEdge[i];

        std::vector<scalar> measured(edges.size(), 0.0);
        for (std::size_t i = 0; i < edges.size(); ++i)
        {
            const std::string& name =
                net.pipeNames()[static_cast<std::size_t>(edges[i])];
            if (source[i] == 0) continue;
            measured[i] = source[i] == 1 ? blocks[i]->regionInertance()
                                         : blocks[i]->probeInertance();
            if (!(measured[i] > 0))
                cfFatal("coupling: transient pipe '" + name + "' is to take its"
                        " inertance from "
                        + std::string(source[i] == 1
                                          ? "its partner's region"
                                          : "the slab of it the partner measures")
                        + ", which measures " + std::to_string(measured[i])
                        + " kg/m4: a partner that states no region has nothing to"
                          " install");
        }
        for (std::size_t i = 0; i < edges.size(); ++i)
        {
            if (source[i] != 0) net.setCoupledInertia(edges[i], measured[i]);
        }
    }

    // Which edges' rows pay for their own region's changing profile. Checked
    // here, once, so a case that names a source for a number of edges the
    // interface does not have fails before it has solved anything. Nothing is
    // installed yet: every row starts on the shape account it was given by the
    // pass before, and the first pass of the first step reports the account
    // without charging it, exactly as the run that never asks does.
    std::vector<int> charge(edges.size(), 0);
    {
        const std::vector<int>& perEdge = ctl.convectiveIntoRowPerEdge;
        if (!perEdge.empty() && perEdge.size() != edges.size())
            cfFatal("coupling: transient controls name a shape charge for "
                    + std::to_string(perEdge.size()) + " edges but the interface"
                      " has " + std::to_string(edges.size()));
        for (std::size_t i = 0; i < perEdge.size(); ++i)
        {
            if (perEdge[i] != 0 && perEdge[i] != 1)
                cfFatal("coupling: transient edge " + std::to_string(i)
                        + " is asked to charge its row with "
                        + std::to_string(perEdge[i])
                        + ", which is neither report (0) nor charge (1)");
            charge[i] = perEdge[i];
        }
    }

    for (std::size_t s = 0; s < nSteps; ++s)
    {
        TransientCoupledStep step;
        label passes = 0;
        bool agreed = false;
        if (qn)
            qn->reset();

        // The level's history goes to the partners before anything is solved:
        // the committed rate is what this level accelerates away from, and a
        // partner that carries its own `du/dt` takes the field it stands on here
        // and only here, so no number of re-solves of this level can move it.
        std::vector<scalar> qPrev(edges.size(), 0.0);
        {
            const std::vector<double> committed = net.pipeFlows();
            for (std::size_t i = 0; i < edges.size(); ++i)
            {
                qPrev[i] = std::abs(committed[static_cast<std::size_t>(edges[i])]);
                blocks[i]->beginTimeLevel(dt, qPrev[i]);
            }
        }

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
                ex.inertance = net.pipeInertia(edge);
                ex.unsteady = blocks[i]->unsteadyDrop();
                ex.qPrev = qPrev[i];
                ex.dq = ex.q - qPrev[i];
                // What the region's momentum is worth is a claim the row makes,
                // so it is kept out of the equation the interface iterates on:
                // the partner's friction is what the resistance is charged with,
                // whether or not the partner also carried an acceleration.
                ex.dpInertialPa = ex.unsteady ? blocks[i]->lastInertialDrop() : 0.0;
                ex.dpFrictionPa = ex.dpPa - ex.dpInertialPa;
                // The shape is part of what this level asked for, so it is read
                // back level by level rather than once at the end of the run: the
                // inlet is the declared profile, the outlet is the region's own,
                // and which way the two differ says whether the region had to grow
                // the profile or shed the excess. It is read before the chord
                // because on an edge whose row pays for it, the difference is not
                // the block's alone to carry any more.
                const ProfileLedger ledger = blocks[i]->profileLedger();
                ex.inletShape = ledger.shape;
                ex.fluxImposed = ledger.fluxImposed;
                ex.betaImposed = ledger.betaImposed;
                ex.betaIn = ledger.betaIn;
                ex.betaOut = ledger.betaOut;
                ex.inletArea = ledger.inletArea;
                // The same solve's pressure stations, reduced to the one number
                // that says whether the inlet was fuller than the region: the
                // height of the highest station above the first. Captured here,
                // with the rest of the ledger, so the level keeps its own rise.
                ex.pRisePa = ledger.pRisePa;
                ex.pRiseX = ledger.pRiseX;
                ex.dpConvectivePa = ledger.dpConvectivePa;
                // Who pays for it is the case's decision, taken per edge, and it
                // is stated once for both sides of it: the row's new term and the
                // part of the block's friction the resistance gives up are the
                // same number. An edge that was only asked to report it leaves
                // the account where it was -- on the line, in no equation.
                ex.convectiveRowPa = charge[i] ? ex.dpConvectivePa : 0.0;
                const scalar base = ex.unsteady ? ex.dpFrictionPa : ex.dpPa;
                ex.chordPa = base - ex.convectiveRowPa;
                if (charge[i] && !(ex.chordPa > 0))
                    cfFatal("coupling: transient pipe '" + ex.name + "' has "
                            + fmt(base) + " Pa for its resistance to carry, of"
                              " which its own region's changing profile is worth "
                            + fmt(ex.dpConvectivePa)
                            + " Pa: charging this edge's row with the shape account"
                              " would leave " + fmt(ex.chordPa)
                              + " Pa for the wall, so the friction this region"
                                " reports cannot pay for the profile it reports");
                if (!(ex.chordPa > 0))
                    cfFatal("coupling: transient pipe '" + ex.name
                            + "' reports a drop of " + fmt(ex.dpPa) + " Pa of which"
                              " the acceleration of its own region is worth "
                            + fmt(ex.dpInertialPa)
                            + " Pa, leaving no friction for a resistance to carry:"
                              " the inertance installed on this edge (M "
                            + fmt(ex.inertance) + " kg/m4 over a rate change of "
                            + fmt(ex.dq) + " m3/s in " + fmt(dt)
                            + " s) is larger than the momentum the region itself"
                              " can account for");
                ex.residual = std::abs(ex.resistance * ex.q - ex.chordPa)
                              / std::max(std::abs(ex.chordPa), 1e-300);
                // The row as the 1D side writes it, against the drop as the
                // region's own momentum balance reports it. This is the number
                // the interface cannot iterate away, and the only one that says
                // whether the inertance installed above belongs to this region.
                // The shape term belongs here too, because the chord it came from
                // stopped carrying it: split the account between the two sides or
                // not, the row answers for the same statement either way.
                ex.rowDropPa = ex.resistance * ex.q
                               + ex.inertance * ex.dq / static_cast<scalar>(dt)
                               + ex.convectiveRowPa;
                ex.rowError = ex.unsteady
                              ? (ex.rowDropPa - ex.dpPa) / ex.dpPa : 0.0;
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
            if (qn)
            {
                std::vector<scalar> x(R.size()), r(R.size());
                for (std::size_t i = 0; i < R.size(); ++i)
                {
                    const EdgeExchange& ex = rec.edges[i];
                    x[i] = ex.resistance / Rseed[i];
                    r[i] = (ex.chordPa / ex.q - ex.resistance) / Rseed[i];
                }
                const std::vector<scalar> next = qn->next(x, r, ctl.relax);
                for (std::size_t i = 0; i < R.size(); ++i)
                {
                    R[i] = next[i] * Rseed[i];
                    net.setCoupledResistance(edges[i], R[i]);
                }
            }
            else
            {
                for (std::size_t i = 0; i < edges.size(); ++i)
                {
                    const EdgeExchange& ex = rec.edges[i];
                    R[i] = (1.0 - ctl.relax) * R[i] + ctl.relax * (ex.chordPa / ex.q);
                    net.setCoupledResistance(edges[i], R[i]);
                }
            }

            // and the row's shape term is installed with it, reduced to the
            // coefficient the row multiplies by |Q|*Q: the two refreshes come
            // from the same solve, so an edge that was told to charge its row
            // never has its account paid twice, once by the chord and once by
            // the row. An edge that was told only to report it is installed
            // with 0, which is the term the row has been carrying all along.
            for (std::size_t i = 0; i < edges.size(); ++i)
            {
                const EdgeExchange& ex = rec.edges[i];
                net.setCoupledConvective(
                    edges[i], ex.convectiveRowPa / (ex.q * ex.q));
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
