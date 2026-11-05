// chflow-3d-parallel: domain decomposition, halo exchange and the global
// reductions a decomposed run needs. The mesh is partitioned by a recursive
// coordinate-median split, every domain gets its own cell and face addressing
// with the geometry its stencil reads, and the explicit Gauss operators are then
// evaluated domain by domain with neighbour values arriving through the halo.
//
// What is owed here is that a decomposed sweep *is* the serial sweep: the same
// faces in the same order, the same bounded interpolation weights, the same
// signs, so the assembled fields agree bit for bit -- and a divergence-free flux
// keeps telescoping to zero across a cut because both sides of a coupled face
// use one flux value.
//
// Two fields are run through it. A linear one is the exactness anchor: cell to
// face interpolation is exact for it and the Gauss sum is then the analytic
// gradient, on hex and tetrahedral meshes alike. A quadratic one cannot be exact
// (interpolating a quadratic at a face is not the face value, and a boundary
// cell has no opposite face to cancel the error), so it is used only for the
// decomposed-equals-serial comparison, where the rule is applied on both sides.
// The flux rides on U = (y, z, 0), divergence-free, so a cell with no face on the
// outer boundary must read a zero net out-flux whatever the mesh.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/VtkWriter.h"
#include "cf/parallel/DomainDecomposition.h"

using namespace cf;

namespace
{

std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

/// T = c0 + c1 x + c2 y + c3 z + c4 x^2 + c5 y^2 + c6 z^2 + c7 xy + c8 yz + c9 zx.
/// With c4..c9 zero this is the linear anchor field.
struct Polynomial
{
    std::vector<double> c;

    scalar operator()(const Vector& p) const
    {
        return c[0] + c[1] * p.x + c[2] * p.y + c[3] * p.z
             + c[4] * p.x * p.x + c[5] * p.y * p.y + c[6] * p.z * p.z
             + c[7] * p.x * p.y + c[8] * p.y * p.z + c[9] * p.z * p.x;
    }

    Vector grad(const Vector& p) const
    {
        return {c[1] + 2.0 * c[4] * p.x + c[7] * p.y + c[9] * p.z,
                c[2] + 2.0 * c[5] * p.y + c[7] * p.x + c[8] * p.z,
                c[3] + 2.0 * c[6] * p.z + c[8] * p.y + c[9] * p.x};
    }
};

/// The flux the bookkeeping test rides on: U = (y, z, 0), divergence free.
scalar analyticFlux(const Vector& faceCentre, const Vector& faceSf)
{
    return faceCentre.y * faceSf.x + faceCentre.z * faceSf.y;
}

scalar largestDeviation(const Field<scalar>& a, const Field<scalar>& b)
{
    scalar m = 0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

scalar largestComponentDeviation(const Field<Vector>& a, const Field<Vector>& b)
{
    scalar m = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        m = std::max(m, std::abs(a[i].x - b[i].x));
        m = std::max(m, std::abs(a[i].y - b[i].y));
        m = std::max(m, std::abs(a[i].z - b[i].z));
    }
    return m;
}

label countNonFinite(const Field<scalar>& f)
{
    label n = 0;
    for (const scalar v : f) if (!std::isfinite(v)) ++n;
    return n;
}

label countNonFinite(const Field<Vector>& f)
{
    label n = 0;
    for (const Vector& v : f)
    {
        if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z)) ++n;
    }
    return n;
}

void splitComponents(const Field<Vector>& g, Field<scalar>& x, Field<scalar>& y, Field<scalar>& z)
{
    x.assign(g.size(), 0.0);
    y.assign(g.size(), 0.0);
    z.assign(g.size(), 0.0);
    for (std::size_t c = 0; c < g.size(); ++c)
    {
        x[c] = g[c].x;
        y[c] = g[c].y;
        z[c] = g[c].z;
    }
}

/// One field, run through the decomposition against the serial operators.
struct FieldRun
{
    std::string tag;
    Polynomial p;
    Field<scalar> t;
    Field<scalar> gradSerialX, gradSerialY, gradSerialZ;
    Field<scalar> gradParallelX, gradParallelY, gradParallelZ;
    scalar scale = 0;
    scalar haloHeld = 0, haloWanted = 0;
    label haloSlots = 0;
};

void runField(const polyMesh& mesh, const parallel::Decomposition& dec,
              const parallel::DomainSet& ds, FieldRun& run)
{
    const label nC = mesh.nCells();
    const auto& centres = mesh.cellCentres();

    run.t.assign(static_cast<std::size_t>(nC), 0.0);
    run.scale = 0;
    for (label c = 0; c < nC; ++c)
    {
        const Vector& p = centres[static_cast<std::size_t>(c)];
        run.t[static_cast<std::size_t>(c)] = run.p(p);
        run.scale = std::max(run.scale, std::abs(run.p(p)));
    }

    // ------------------------------------------------------------------
    // halo: every peer slot holds exactly the value its owner holds
    // ------------------------------------------------------------------
    std::vector<Field<scalar>> localT = ds.storage<scalar>();
    ds.scatter(run.t, localT);
    const label moved = parallel::exchange(dec, localT);

    label haloMismatch = 0;
    run.haloHeld = 0;
    run.haloWanted = 0;
    run.haloSlots = 0;
    for (label d = 0; d < dec.nDomains(); ++d)
    {
        const Field<scalar>& l = localT[static_cast<std::size_t>(d)];
        const std::vector<label>& halo = dec.haloCells()[static_cast<std::size_t>(d)];
        run.haloSlots += static_cast<label>(halo.size());
        for (std::size_t i = 0; i < halo.size(); ++i)
        {
            const scalar held = l[static_cast<std::size_t>(dec.nCells(d) + static_cast<label>(i))];
            const scalar wanted = run.t[static_cast<std::size_t>(halo[i])];
            run.haloHeld += held;
            run.haloWanted += wanted;
            if (held != wanted) ++haloMismatch;
        }
    }
    // A coupled face sends one value each way; a peer cell can be reached over
    // several coupled faces, so it is one halo slot but more than one message.
    if (moved != 2 * dec.nCoupledFaces())
        cfFatal("parallel: exchange moved " + std::to_string(moved) + " values for "
                + std::to_string(dec.nCoupledFaces()) + " coupled faces");
    if (haloMismatch)
        cfFatal("parallel: " + std::to_string(haloMismatch) + " halo slots did not arrive with their owner's value");

    // ------------------------------------------------------------------
    // stencils: domain by domain must be the serial sweep, same face rule
    // ------------------------------------------------------------------
    const Field<scalar> ow = makeInternalOwnerWeights(mesh);
    const auto& owner = mesh.owner();
    const auto& nbr = mesh.neighbour();
    const auto& fCentres = mesh.faceCentres();
    const label nInt = mesh.nInternalFaces();

    // the rule the domain applies: bounded linear interpolation across an
    // internal (or coupled) face, the prescribed value at a real boundary
    const parallel::FaceValue faceT = [&run](label, const Vector& fc) { return run.p(fc); };
    const auto serialFaceT = [&run, &owner, &nbr, &nInt, &ow, &fCentres](
                                 label f, const Field<scalar>& data) -> scalar
    {
        if (f < nInt)
        {
            const scalar w = ow[static_cast<std::size_t>(f)];
            return w * data[static_cast<std::size_t>(owner[static_cast<std::size_t>(f)])]
                 + (1.0 - w) * data[static_cast<std::size_t>(nbr[static_cast<std::size_t>(f)])];
        }
        return run.p(fCentres[static_cast<std::size_t>(f)]);
    };

    Field<Vector> analytic(static_cast<std::size_t>(nC), Vector::zero());
    for (label c = 0; c < nC; ++c) analytic[static_cast<std::size_t>(c)] = run.p.grad(centres[static_cast<std::size_t>(c)]);

    for (const std::string scheme : {"linear", "corrected"})
    {
        const Field<Vector> serial = fvc::grad(mesh, run.t, serialFaceT, scheme);
        const Field<Vector> parallelG = parallel::assembledGrad(ds, run.t, faceT, scheme);
        if (countNonFinite(parallelG)) cfFatal("parallel: non-finite gradient (" + run.tag + " " + scheme + ")");

        const scalar gradAnalytic = largestComponentDeviation(analytic, serial);
        const scalar gradDev = largestComponentDeviation(serial, parallelG);
        log::info("parallel stencil domains ", dec.nDomains(), " field ", run.tag, " scheme ", scheme,
                  " gradAnalytic ", num(gradAnalytic), " gradAnalyticRel ", num(gradAnalytic / run.scale),
                  " gradDev ", num(gradDev), " gradDevRel ", num(gradDev / run.scale),
                  " haloDrift ", num(std::abs(run.haloHeld - run.haloWanted)),
                  " scale ", num(run.scale));

        if (scheme == "corrected")
        {
            splitComponents(serial, run.gradSerialX, run.gradSerialY, run.gradSerialZ);
            splitComponents(parallelG, run.gradParallelX, run.gradParallelY, run.gradParallelZ);
        }
    }

    // ------------------------------------------------------------------
    // reductions: partials per domain, one global answer
    // ------------------------------------------------------------------
    scalar serialSum = 0, serialMin = std::numeric_limits<double>::infinity(), serialMax = -serialMin;
    std::vector<scalar> sums, mins, maxs;
    for (label d = 0; d < dec.nDomains(); ++d)
    {
        scalar p = 0, mn = std::numeric_limits<double>::infinity(), mx = -mn;
        for (const label c : dec.cells()[static_cast<std::size_t>(d)])
        {
            const std::size_t u = static_cast<std::size_t>(c);
            p += run.t[u] * mesh.cellVolumes()[u];
            mn = std::min(mn, run.t[u]);
            mx = std::max(mx, run.t[u]);
        }
        sums.push_back(p);
        mins.push_back(mn);
        maxs.push_back(mx);
        serialSum += p;
        serialMin = std::min(serialMin, mn);
        serialMax = std::max(serialMax, mx);
    }
    const scalar sum = parallel::reduceSum(sums);
    if (countNonFinite(run.t) || !std::isfinite(sum)) cfFatal("parallel: non-finite field or reduction");
    const scalar sumDev = std::abs(sum - serialSum);
    const scalar minDev = std::abs(parallel::reduceMin(mins) - serialMin);
    const scalar maxDev = std::abs(parallel::reduceMax(maxs) - serialMax);

    log::info("parallel reduce domains ", dec.nDomains(), " field ", run.tag,
              " sum ", num(sum), " sumSerial ", num(serialSum),
              " minSerial ", num(serialMin), " maxSerial ", num(serialMax),
              " sumDev ", num(sumDev), " sumDevRel ", num(sumDev / std::abs(serialSum)),
              " minDev ", num(minDev), " maxDev ", num(maxDev),
              " haloChecksum ", num(run.haloHeld), " haloSlots ", run.haloSlots);
}

/// One decomposition, checked against the serial operators and reported.
void runOne(const polyMesh& mesh, label nDomains, const std::string& vtkName,
            std::vector<FieldRun>& runs)
{
    const label nC = mesh.nCells();
    const parallel::Decomposition dec(mesh, nDomains);
    const parallel::DomainSet ds(dec);

    // ------------------------------------------------------------------
    // partition: every cell owned exactly once, every face classified
    // ------------------------------------------------------------------
    Field<label> seen(static_cast<std::size_t>(nC), 0);
    for (label d = 0; d < dec.nDomains(); ++d)
    {
        for (const label c : dec.cells()[static_cast<std::size_t>(d)])
        {
            if (++seen[static_cast<std::size_t>(c)] > 1)
                cfFatal("parallel: cell " + std::to_string(c) + " owned by two domains");
        }
    }
    label uncovered = 0;
    for (const label s : seen) if (s == 0) ++uncovered;
    if (uncovered) cfFatal("parallel: " + std::to_string(uncovered) + " cells owned by no domain");

    label coupledEntries = 0, boundaryEntries = 0, internalEntries = 0;
    for (label d = 0; d < dec.nDomains(); ++d)
    {
        label coupled = 0, boundary = 0, internal = 0;
        for (const parallel::DomainFace& df : dec.faces(d))
        {
            if (df.peer < 0) ++boundary;
            else if (df.peer >= dec.nCells(d)) ++coupled;
            else ++internal;
        }
        coupledEntries += coupled;
        boundaryEntries += boundary;
        internalEntries += internal;

        log::info("parallel domain domains ", nDomains, " d ", d,
                  " cells ", dec.nCells(d), " halo ", dec.nHalo(d),
                  " faces ", dec.faces(d).size(),
                  " internal ", internal, " coupled ", coupled, " boundary ", boundary);
    }

    // A coupled face is listed by both of its domains, a boundary face by one,
    // and a face whose two cells share a domain contributes two entries there --
    // so the three counts have to reassemble the global mesh exactly.
    const label coupledDistinct = dec.nCoupledFaces();
    if (coupledEntries != 2 * coupledDistinct)
        cfFatal("parallel: " + std::to_string(coupledEntries) + " coupled entries for "
                + std::to_string(coupledDistinct) + " coupled faces");
    if (boundaryEntries != mesh.nBoundaryFaces())
        cfFatal("parallel: domain boundary entries " + std::to_string(boundaryEntries)
                + " but the mesh has " + std::to_string(mesh.nBoundaryFaces()));
    if (internalEntries % 2)
        cfFatal("parallel: an odd number of local internal entries");
    if (internalEntries / 2 + coupledDistinct != mesh.nInternalFaces())
        cfFatal("parallel: local internal " + std::to_string(internalEntries / 2)
                + " plus coupled " + std::to_string(coupledDistinct)
                + " but the mesh has " + std::to_string(mesh.nInternalFaces()));

    log::info("parallel run domains ", nDomains,
              " cells ", nC,
              " coupled ", coupledDistinct,
              " messages ", dec.nMessages(),
              " cellsMin ", dec.minCellsPerDomain(),
              " cellsMax ", dec.maxCellsPerDomain(),
              " uncovered ", uncovered);

    // ------------------------------------------------------------------
    // flux bookkeeping across a cut: a divergence-free field has no source
    // ------------------------------------------------------------------
    const auto& fCentres = mesh.faceCentres();
    const auto& sf = mesh.faceSf();
    Field<scalar> phi(static_cast<std::size_t>(mesh.nFaces()), 0.0);
    scalar phiScale = 0;
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const scalar v = analyticFlux(fCentres[static_cast<std::size_t>(f)], sf[static_cast<std::size_t>(f)]);
        phi[static_cast<std::size_t>(f)] = v;
        phiScale = std::max(phiScale, std::abs(v));
    }
    const Field<scalar> serialDiv = fvc::divPhi(mesh, phi);
    const Field<scalar> parallelDiv = parallel::assembledDivPhi(ds,
        [](label, const Vector& fc, const Vector& s) { return analyticFlux(fc, s); });
    if (countNonFinite(parallelDiv)) cfFatal("parallel: non-finite net out-flux");
    const scalar divDev = largestDeviation(serialDiv, parallelDiv);

    // Only cells carrying a face on the outer boundary see a real imbalance.
    const Field<label>& facePatch = mesh.facePatchMap();
    scalar divInterior = 0;
    label interiorCells = 0;
    for (label c = 0; c < nC; ++c)
    {
        bool shell = false;
        for (const label f : mesh.cells()[static_cast<std::size_t>(c)])
        {
            if (facePatch[static_cast<std::size_t>(f)] >= 0) shell = true;
        }
        if (shell) continue;
        ++interiorCells;
        divInterior = std::max(divInterior, std::abs(parallelDiv[static_cast<std::size_t>(c)]));
    }
    if (interiorCells == 0) cfFatal("parallel: no interior cell to test flux telescoping on");

    log::info("parallel flux domains ", nDomains,
              " divDev ", num(divDev), " divInterior ", num(divInterior),
              " interiorCells ", interiorCells, " phiScale ", num(phiScale));

    for (FieldRun& run : runs) runField(mesh, dec, ds, run);

    // ------------------------------------------------------------------
    // fields, so the checker can weigh them itself
    // ------------------------------------------------------------------
    Field<scalar> domainOfCell(static_cast<std::size_t>(nC), 0.0);
    Field<scalar> haloCount(static_cast<std::size_t>(nC), 0.0);
    for (label c = 0; c < nC; ++c) domainOfCell[static_cast<std::size_t>(c)] = static_cast<scalar>(dec.domainOfCell(c));
    for (label d = 0; d < dec.nDomains(); ++d)
    {
        for (const label c : dec.haloCells()[static_cast<std::size_t>(d)])
            haloCount[static_cast<std::size_t>(c)] += 1.0;
    }

    VtkWriter writer(mesh);
    Field<scalar> divSerial = serialDiv, divParallel = parallelDiv;
    writer.addCellScalar("domain", domainOfCell);
    writer.addCellScalar("haloCount", haloCount);
    writer.addCellScalar("divPhiSerial", divSerial);
    writer.addCellScalar("divPhiParallel", divParallel);
    for (const FieldRun& run : runs)
    {
        const std::string tag = run.tag;
        writer.addCellScalar("t" + tag, run.t);
        writer.addCellScalar("gradSerial" + tag + "X", run.gradSerialX);
        writer.addCellScalar("gradSerial" + tag + "Y", run.gradSerialY);
        writer.addCellScalar("gradSerial" + tag + "Z", run.gradSerialZ);
        writer.addCellScalar("gradParallel" + tag + "X", run.gradParallelX);
        writer.addCellScalar("gradParallel" + tag + "Y", run.gradParallelY);
        writer.addCellScalar("gradParallel" + tag + "Z", run.gradParallelZ);
    }
    writer.write(vtkName);
}

Polynomial parsePolynomial(const Dictionary& dict, const std::string& key)
{
    Polynomial p{dict.getDict("field").getScalarList(key)};
    if (p.c.size() != 10) cfFatal("parallel: field " + key + " needs ten coefficients, c0..c9");
    return p;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        log::error("usage: chflow-3d-parallel <case.dict>");
        return 1;
    }
    const std::string casePath = argv[1];
    try
    {
        const Dictionary dict = Dictionary::parseFile(casePath);
        const polyMesh mesh = MeshLoader::fromCase(dict, casePath);
        const std::vector<double> domainList = dict.getDict("decompose").getScalarList("domains");
        if (domainList.empty()) cfFatal("decompose domains needs at least one count");

        std::vector<FieldRun> runs;
        {
            FieldRun lin; lin.tag = "Linear"; lin.p = parsePolynomial(dict, "linear");
            FieldRun quad; quad.tag = "Quad"; quad.p = parsePolynomial(dict, "quadratic");
            runs.push_back(lin);
            runs.push_back(quad);
        }

        const std::string vtkName = dict.has("output")
            ? dict.getDict("output").getWordOr("vtk", "out_parallel.vtu") : "out_parallel.vtu";

        log::info("parallel start cells ", mesh.nCells(), " faces ", mesh.nFaces(),
                  " internal ", mesh.nInternalFaces(), " boundary ", mesh.nBoundaryFaces(),
                  " patches ", mesh.nPatches(), " runs ", static_cast<label>(domainList.size()));

        for (const double k : domainList)
        {
            // one file per domain count: the fields must not depend on the cut
            const std::string name = vtkName.substr(0, vtkName.size() - 4)
                                   + "_" + std::to_string(static_cast<label>(k)) + ".vtu";
            runOne(mesh, static_cast<label>(k), name, runs);
        }
        log::info("parallel done");
    }
    catch (const std::exception& e)
    {
        log::error("parallel: ", e.what());
        return 1;
    }
    return 0;
}
