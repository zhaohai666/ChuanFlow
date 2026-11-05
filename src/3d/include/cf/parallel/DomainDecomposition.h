#pragma once

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "cf/common/Error.h"
#include "cf/common/Field.h"
#include "cf/common/Labels.h"
#include "cf/common/Vector.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{
namespace parallel
{

/// A stencil face as one domain sees it. `cell` and `peer` are local slots:
/// cells of this domain come first, in decomposition order, then the halo
/// cells other domains sent it. `sign` is +1 when `cell` is the global owner
/// and -1 when it is the global neighbour, which is how the face's flux enters
/// that cell's Gauss sum. A real boundary face has peer = -1.
struct DomainFace
{
    label face = -1;
    label cell = -1;
    label peer = -1;
    scalar sign = 1.0;
};

/// The cells of a polyMesh split over `nDomains` ranks, plus the face topology
/// the split creates: which faces stay inside a domain, which cross a domain
/// boundary (the coupled, or processor, faces) and which local slot of which
/// peer holds the cell on the other side of each of them.
///
/// The split is a recursive coordinate-median ordering of the cell centres --
/// halve the block at the median of its longest axis, order each half the same
/// way, then cut the ordered list into contiguous slices. Cells that end up in
/// one slice are spatially compact, and the slices differ by at most one cell.
class Decomposition
{
public:
    Decomposition(const polyMesh& mesh, label nDomains);

    const polyMesh& mesh() const { return mesh_; }
    label nDomains() const { return nDomains_; }
    label nCells() const { return mesh_.nCells(); }

    label domainOfCell(label cell) const { return cellDomain_[static_cast<std::size_t>(cell)]; }
    label nCells(label d) const { return static_cast<label>(cells_[d].size()); }
    label nHalo(label d) const { return static_cast<label>(haloCells_[d].size()); }
    /// Local storage a domain needs: owned slots, then halo slots.
    label nLocal(label d) const { return nCells(d) + nHalo(d); }

    /// Domain -> owned global cells, in local slot order.
    const std::vector<std::vector<label>>& cells() const { return cells_; }
    /// Domain -> the peer cells it haloed, in halo slot order (slot = nCells + i).
    const std::vector<std::vector<label>>& haloCells() const { return haloCells_; }

    /// Every face with at least one side owned here, in ascending global order.
    const std::vector<DomainFace>& faces(label d) const { return faces_[d]; }

    /// Local slot of `cell` inside domain `d` (owned or halo), or -1 if `d`
    /// neither owns nor haloed it.
    label localSlot(label d, label cell) const;

    /// Message plan. sendPlan(d)[i] names a peer and the owned slots of `d`
    /// that go to it; the peer's receive list from `d` holds the halo slots
    /// those values land in, in the same order -- so the two lists are the
    /// send and receive buffers of one message, contiguous, the shape a
    /// message-passing transport posts without rearranging anything.
    struct Message
    {
        label peer = -1;
        std::vector<label> slots;
    };
    const std::vector<Message>& sendPlan(label d) const { return send_[d]; }
    const std::vector<Message>& recvPlan(label d) const { return recv_[d]; }
    /// Receive slots of `d` coming from `peer`, or null when there is no such
    /// message. Index-matched with the sender's send list by construction.
    const Message* recvFrom(label d, label peer) const;

    /// Distinct global faces that cross a domain boundary.
    label nCoupledFaces() const { return nCoupled_; }
    /// Total messages = sum over domains of its send-list length.
    label nMessages() const;
    label minCellsPerDomain() const;
    label maxCellsPerDomain() const;

private:
    /// Record `cell` (and, across an internal face, `peer`) in domain `d`'s
    /// stencil, allocating a halo slot for `peer` if `d` does not own it.
    void addSide(label d, label face, label cell, label peer, scalar sign);
    /// Local slot of `cell` in `d`, haloing the cell if `d` needs it and does
    /// not own it yet.
    label slotOf(label d, label cell);

    const polyMesh& mesh_;
    label nDomains_;
    Field<label> cellDomain_;
    std::vector<std::vector<label>> cells_;
    std::vector<std::vector<label>> haloCells_;
    std::vector<std::vector<DomainFace>> faces_;
    std::vector<std::vector<Message>> send_, recv_;
    /// domain -> global cell -> local slot (owned slots first, then haloed ones)
    std::vector<std::map<label, label>> slot_;
    label nCoupled_ = 0;
};

/// One domain's rank-local view: its addressing and a copy of the geometry the
/// global mesh computed for its faces. A stencil loop over a DomainMesh reads
/// nothing outside it, which is what a rank in a distributed run could do; the
/// geometry is carried rather than reached for, exactly as a mesh decomposer
/// writes a subdomain's processor-patch geometry.
class DomainMesh
{
public:
    explicit DomainMesh(const Decomposition& dec, label domain);

    label domain() const { return domain_; }
    label nOwned() const { return static_cast<label>(cells_.size()); }
    label nLocal() const { return nSlots_; }
    label nFaces() const { return static_cast<label>(faces_.size()); }

    const std::vector<label>& cells() const { return cells_; }         // slot -> global cell
    const std::vector<label>& haloCells() const { return haloCells_; } // halo slot -> global cell
    const std::vector<DomainFace>& faces() const { return faces_; }

    const Field<scalar>& cellVolumes() const { return cellVolumes_; }
    const Field<Vector>& cellCentres() const { return cellCentres_; }
    const Field<Vector>& faceSf() const { return sf_; }
    const Field<Vector>& faceCentres() const { return faceCentres_; }
    /// face centre minus the *global owner* centre: the non-orthogonal lever,
    /// identical on both sides of a face so each side computes the same flux.
    const Field<Vector>& faceDelta() const { return delta_; }
    /// bounded interpolation weight of the global owner.
    const Field<scalar>& faceOwnerWeight() const { return ownerWeight_; }

private:
    label domain_ = -1;
    label nSlots_ = 0;
    std::vector<label> cells_;
    std::vector<label> haloCells_;
    std::vector<DomainFace> faces_;
    Field<scalar> cellVolumes_;
    Field<Vector> cellCentres_;
    Field<Vector> sf_;
    Field<Vector> faceCentres_;
    Field<Vector> delta_;
    Field<scalar> ownerWeight_;
};

/// The DomainMeshes of a decomposition, built once and reused by every sweep.
class DomainSet
{
public:
    DomainSet(const Decomposition& dec);

    const Decomposition& decomposition() const { return dec_; }
    label size() const { return static_cast<label>(domains_.size()); }
    const DomainMesh& operator[](label d) const { return domains_[static_cast<std::size_t>(d)]; }

    /// Per-domain local storage: zeroed owned slots then zeroed halo slots.
    template<typename T>
    std::vector<Field<T>> storage() const;

    /// Owned cell values lifted into each domain's local slots.
    void scatter(const Field<scalar>& global, std::vector<Field<scalar>>& local) const;

    /// Local slots of the owned cells gathered back into global cell order.
    Field<scalar> gather(const std::vector<Field<scalar>>& local) const;
    Field<Vector> gather(const std::vector<Field<Vector>>& local) const;

private:
    const Decomposition& dec_;
    std::vector<DomainMesh> domains_;
};

/// Fill every domain's halo slots from its peers' owned slots. Owned slots are
/// only ever read, so messages cannot interfere and no ordering is needed --
/// the reason this is a copy loop here and an Isend/Recv pair in a distributed
/// run, with the same buffers. Returns the number of values moved.
template<typename T>
label exchange(const Decomposition& dec, std::vector<Field<T>>& local);

/// The collective reductions a decomposed run needs for global quantities
/// (residual norms, field min/max, conservation sums): partial per-domain
/// contributions in, one global value out.
scalar reduceSum(const std::vector<scalar>& partials);
scalar reduceMin(const std::vector<scalar>& partials);
scalar reduceMax(const std::vector<scalar>& partials);

// ------------------------------------------------------------------
// Domain-local stencils, assembled back into global cell order. Every
// arithmetic step mirrors its serial counterpart in cf::fvc and reads only
// domain-local data, so a decomposed sweep reproduces the serial one.
// ------------------------------------------------------------------

/// Value of a face quantity at (global face, face centre). A domain hands over
/// the centre it carried with the face, so the evaluator never reaches past the
/// domain; the serial reference passes the same expression evaluated on the
/// global mesh, which is what makes the two comparable.
using FaceValue = std::function<scalar(label face, const Vector& faceCentre)>;

/// Value of a face flux at (global face, face centre, signed area vector). All
/// three come from the domain's own copy, so a flux evaluator written against
/// this signature is rank-local on both sides of a cut.
using FaceFlux = std::function<scalar(label face, const Vector& faceCentre, const Vector& faceSf)>;

/// Gauss gradient evaluated domain by domain; coupled face values arrive
/// through the halo, and the non-orthogonal correction exchanges the gradient
/// between its passes.
Field<Vector> assembledGrad(const DomainSet&, const Field<scalar>& cellValues,
                           const FaceValue& faceValueOf, const std::string& scheme);

/// Per-cell net out-flux: the face fluxes this domain sees, summed with each
/// face's own sign.
Field<scalar> assembledDivPhi(const DomainSet&, const FaceFlux& facePhi);

// ------------------------------------------------------------------
// template definitions
// ------------------------------------------------------------------

template<typename T>
std::vector<Field<T>> DomainSet::storage() const
{
    std::vector<Field<T>> out(static_cast<std::size_t>(dec_.nDomains()));
    for (label d = 0; d < dec_.nDomains(); ++d)
    {
        out[static_cast<std::size_t>(d)] = Field<T>(static_cast<std::size_t>(dec_.nLocal(d)));
    }
    return out;
}

template<typename T>
label exchange(const Decomposition& dec, std::vector<Field<T>>& local)
{
    if (static_cast<label>(local.size()) != dec.nDomains())
    {
        cfFatal("exchange: " + std::to_string(local.size()) + " fields for "
                + std::to_string(dec.nDomains()) + " domains");
    }

    label moved = 0;
    for (label d = 0; d < dec.nDomains(); ++d)
    {
        for (const Decomposition::Message& m : dec.sendPlan(d))
        {
            const Decomposition::Message* r = dec.recvFrom(m.peer, d);
            if (!r || r->slots.size() != m.slots.size())
            {
                cfFatal("exchange: domain " + std::to_string(d) + " sends "
                        + std::to_string(m.slots.size()) + " values to "
                        + std::to_string(m.peer) + ", which did not ask for that many");
            }
            const Field<T>& src = local[static_cast<std::size_t>(d)];
            Field<T>& dst = local[static_cast<std::size_t>(m.peer)];
            for (std::size_t k = 0; k < m.slots.size(); ++k)
            {
                // Reads touch owned slots, writes halo slots, and the two sets
                // never overlap, so messages cannot read each other's results.
                dst[static_cast<std::size_t>(r->slots[k])] = src[static_cast<std::size_t>(m.slots[k])];
                ++moved;
            }
        }
    }
    return moved;
}

} // namespace parallel
} // namespace cf
