#include "cf/multiphase/Plic.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"

namespace cf::vof
{

namespace
{

/// The vertices of a face, in the order the mesh stores them. The swept prism below
/// uses the loop only for its edges and its plane, so the winding does not matter.
std::vector<Vector> faceLoop(const polyMesh& mesh, label face)
{
    const std::vector<label>& loop = mesh.faces()[static_cast<std::size_t>(face)];
    const Field<Vector>& points = mesh.points();
    std::vector<Vector> out;
    out.reserve(loop.size());
    for (label v : loop) out.push_back(points[static_cast<std::size_t>(v)]);
    return out;
}

/// The prism a face sweeps as it moves `distance` along `intoRegion`, which points from
/// the face into the cell whose phase the face is trading. It is the intersection of the
/// slab between the face's plane and the plane it reaches with the side planes through
/// each of its edges, every one of them parallel to the sweep - so the volume of phase
/// that crossed the face is one more clip of the donor's phase region, and the swept
/// volume itself is exactly the flux the face was given: area times distance.
///
/// The side planes are oriented against an interior point of the prism (the face's own
/// centre, pushed half a distance along the sweep) rather than against the loop's
/// winding, which is what keeps this right for a face stored either way round and for a
/// polygon of any number of sides.
std::vector<HalfSpace> sweptPrism(const std::vector<Vector>& loop, const Vector& intoRegion,
                                  scalar distance)
{
    if (mag2(intoRegion) <= 0) cfFatal("plic advection: the sweep direction is zero");
    if (loop.size() < 3) cfFatal("plic advection: a face needs at least three corners to sweep");

    const Vector nHat = intoRegion * (1.0 / mag(intoRegion));
    Vector centre = Vector::zero();
    for (const Vector& p : loop) centre += p;
    centre = centre * (1.0 / static_cast<scalar>(loop.size()));
    const Vector inside = centre + nHat * (0.5 * distance);

    std::vector<HalfSpace> hs;
    hs.push_back({nHat, dot(nHat, loop[0])});
    hs.push_back({Vector{-nHat.x, -nHat.y, -nHat.z}, -dot(nHat, loop[0]) - distance});

    const std::size_t n = loop.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        const Vector edge = loop[(i + 1) % n] - loop[i];
        Vector side = cross(nHat, edge);
        if (mag2(side) <= 0) continue; // a zero-length or sweep-parallel edge bounds nothing
        side = side * (1.0 / mag(side));
        if (dot(side, inside - loop[i]) < 0) side = Vector{-side.x, -side.y, -side.z};
        hs.push_back({side, dot(side, loop[i])});
    }
    return hs;
}

/// How a cell's phase is bounded for a step.
enum class PhaseKind
{
    None,      ///< holds no phase, so it gives nothing away
    Full,      ///< is all phase: every swept volume it passes on is phase
    Plane,     ///< the reconstructed plane bounds it
    Fraction,  ///< in the band but carries no plane, so it trades its fraction uniformly
};

} // namespace

GeometricStep advectGeometric(volScalarField& alpha, const Field<scalar>& phi,
                              scalar dt, const Reconstruction& r)
{
    const polyMesh& mesh = alpha.mesh();
    const label nC = mesh.nCells();
    const std::size_t n = static_cast<std::size_t>(nC);
    const label nFaces = mesh.nFaces();
    const std::size_t nf = static_cast<std::size_t>(nFaces);
    const label nInternal = mesh.nInternalFaces();

    Field<scalar>& a = alpha.internalField();
    if (static_cast<label>(phi.size()) != nFaces)
        cfFatal("plic advection: the flux does not match the mesh");
    if (static_cast<label>(r.normal.size()) != nC)
        cfFatal("plic advection: the reconstruction is of another field");
    if (!(dt > 0)) cfFatal("plic advection: the step is not forward in time");

    const Field<scalar>& volume = mesh.cellVolumes();
    const Field<label>& owner = mesh.owner();
    const Field<label>& neighbour = mesh.neighbour();
    const Field<Vector>& sf = mesh.faceSf();

    // What each cell holds and what bounds it. Only a cell with an interface in it
    // needs a region: a full cell passes on everything it sweeps and an empty one
    // nothing, and both of those are the flux without a clip.
    std::vector<PhaseKind> kind(n, PhaseKind::None);
    std::vector<ConvexRegion> phaseRegion(n);
    Field<scalar> phase(n, 0.0);
    scalar total = 0;
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        if (volume[u] <= 0) cfFatal("plic advection: a cell with no volume");
        if (mag2(r.normal[u]) > 0)
        {
            kind[u] = PhaseKind::Plane;
            phaseRegion[u] = clipRegion(cellRegion(mesh, c), r.normal[u], r.offset[u]);
        }
        else if (a[u] >= 1.0 - r.bandTolerance) kind[u] = PhaseKind::Full;
        else if (a[u] > r.bandTolerance) kind[u] = PhaseKind::Fraction;
        phase[u] = a[u] * volume[u];
        total += phase[u];
    }

    // One flux per face, taken from the upwind cell's geometry and shared with the cell
    // that receives it - which is why the total cannot move.
    std::vector<scalar> flux(nf, 0.0), factor(nf, 1.0);
    std::vector<label> up(nf, -1), down(nf, -1);
    GeometricStep step;

    for (label f = 0; f < nFaces; ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const scalar swept = phi[u] * dt; // signed, owner to neighbour
        if (swept == 0) continue;
        const scalar magnitude = std::abs(swept);
        const scalar area = mag(sf[u]);
        if (!(area > 0)) cfFatal("plic advection: a face with no area");
        const Vector faceNormal = sf[u] * (1.0 / area);

        label donor = owner[u], receiver = -1;
        Vector intoDonor = Vector{-faceNormal.x, -faceNormal.y, -faceNormal.z};
        if (swept < 0)
        {
            donor = f < nInternal ? neighbour[u] : -1;
            receiver = owner[u];
            intoDonor = faceNormal;
        }
        else if (f < nInternal)
        {
            receiver = neighbour[u];
        }
        up[u] = donor;
        down[u] = receiver;

        scalar q = 0;
        if (donor < 0)
        {
            // the boundary is a reservoir at its patch's own fraction: its interface is
            // not reconstructed, so what enters is that fraction of the swept volume
            q = alpha.faceValue(f) * magnitude;
        }
        else switch (kind[static_cast<std::size_t>(donor)])
        {
            case PhaseKind::Plane:
                q = clipRegion(phaseRegion[static_cast<std::size_t>(donor)],
                               sweptPrism(faceLoop(mesh, f), intoDonor,
                                          magnitude / area)).volume();
                break;
            case PhaseKind::Full: q = magnitude; break;
            case PhaseKind::Fraction: q = a[static_cast<std::size_t>(donor)] * magnitude; break;
            case PhaseKind::None: q = 0; break;
        }

        if (q <= 0) continue;
        flux[u] = q;
        step.fluxFaces++;
        if (donor >= 0)
            step.maxCourant = std::max(step.maxCourant,
                                       magnitude / volume[static_cast<std::size_t>(donor)]);
    }

    std::vector<scalar> outSum(n, 0.0), inSum(n, 0.0), cellFactor(n, 1.0);
    const auto sumFluxes = [&]()
    {
        std::fill(outSum.begin(), outSum.end(), 0.0);
        std::fill(inSum.begin(), inSum.end(), 0.0);
        for (label f = 0; f < nFaces; ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            if (flux[u] <= 0) continue;
            if (up[u] >= 0) outSum[static_cast<std::size_t>(up[u])] += flux[u];
            if (down[u] >= 0) inSum[static_cast<std::size_t>(down[u])] += flux[u];
        }
    };

    // Two limits, applied round by round until they agree. A cell hands over no more
    // than it holds - the outflow faces of one cell sweep prisms that overlap, since
    // the corner two of them both pass through is counted twice, so the shortfall is
    // shared among them as one factor. A cell receives no more than its void plus what
    // it is itself passing on, which is why a full cell in the middle of the phase keeps
    // trading instead of being told it has no room and starving the step. Both factors
    // travel along the face, so the two cells of a face keep the same flux and the sums
    // still cancel; each round only shrinks a flux, so the pair converges, and at its
    // fixed point both bounds hold at once.
    const auto limitBy = [&](const std::vector<label>& side)
    {
        scalar change = 0;
        for (label f = 0; f < nFaces; ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            if (flux[u] <= 0 || side[u] < 0) continue;
            const scalar fac = cellFactor[static_cast<std::size_t>(side[u])];
            if (fac >= 1.0) continue;
            factor[u] *= fac;
            flux[u] *= fac;
            change = std::max(change, 1.0 - fac);
        }
        return change;
    };

    for (int round = 0; round < 200; ++round)
    {
        sumFluxes();
        for (label c = 0; c < nC; ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            cellFactor[u] = outSum[u] > 0 ? std::min(1.0, phase[u] / outSum[u]) : 1.0;
        }
        scalar change = limitBy(up);

        sumFluxes();
        for (label c = 0; c < nC; ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            const scalar room = std::max(0.0, volume[u] - phase[u] + outSum[u]);
            cellFactor[u] = inSum[u] > 0 ? std::min(1.0, room / inSum[u]) : 1.0;
        }
        change = std::max(change, limitBy(down));

        step.limiterRounds = round + 1;
        if (change <= 1e-13) break;
    }

    for (label f = 0; f < nFaces; ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        if (flux[u] <= 0) continue;
        if (factor[u] < 1.0 - 1e-12)
        {
            step.limitedFaces++;
            step.minFactor = std::min(step.minFactor, factor[u]);
        }
        if (up[u] >= 0) phase[static_cast<std::size_t>(up[u])] -= flux[u];
        if (down[u] >= 0) phase[static_cast<std::size_t>(down[u])] += flux[u];
    }

    scalar newTotal = 0, clipped = 0;
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        // the fixed point of the two limits already keeps the phase inside its cell;
        // this is the round-off and the rounds cap, counted rather than hidden
        const scalar p = std::min(volume[u], std::max(0.0, phase[u]));
        clipped += std::abs(p - phase[u]);
        a[u] = p / volume[u];
        newTotal += p;
    }
    step.clampVolume = total > 0 ? clipped / total : 0;
    step.phaseVolume = newTotal;
    step.volumeError = newTotal - total;
    step.minAlpha = alpha.min();
    step.maxAlpha = alpha.max();
    return step;
}

} // namespace cf::vof
