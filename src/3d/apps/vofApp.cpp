// chflow-vof: transport a phase fraction in a prescribed divergence-free
// velocity and report the interface geometry it implies.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/VtkWriter.h"
#include "cf/multiphase/VoF.h"

using namespace cf;

namespace
{

std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}

/// The two analytic fields the increment is verified with. Both are linear in
/// the coordinates, which makes the discrete flux divergence of every cell
/// exactly zero -- summing U(x_f) & Sf over a closed polyhedron is the integral
/// of div(U) whenever U is linear -- so the transport is conservative and
/// bounded by construction rather than by luck.
struct VelocityField
{
    enum Kind { Uniform, Rotation } kind = Uniform;
    Vector u = Vector::zero();      // Uniform
    scalar omega = 0;               // Rotation
    Vector center = Vector::zero(); // Rotation

    Vector at(const Vector& x) const
    {
        if (kind == Uniform) return u;
        const Vector r = x - center;
        return omega * Vector{-r.y, r.x, 0.0};
    }
};

VelocityField loadVelocity(const Dictionary& caseDict)
{
    if (!caseDict.has("velocity")) cfFatal("vof: case has no velocity{} block");
    const Dictionary& vd = caseDict.getDict("velocity");
    const std::string type = vd.getWord("type");
    VelocityField v;
    if (type == "uniform")
    {
        const auto l = vd.getScalarList("U");
        v.kind = VelocityField::Uniform;
        v.u = Vector{l[0], l[1], l[2]};
    }
    else if (type == "rotation")
    {
        v.kind = VelocityField::Rotation;
        v.omega = vd.getScalar("omega");
        const auto l = vd.getScalarList("center");
        v.center = Vector{l[0], l[1], l[2]};
    }
    else
    {
        cfFatal("vof: velocity type must be uniform or rotation");
    }
    return v;
}

/// The initial interface. `plane` and `circle` without `smooth` fill whole
/// cells, which is the standard VOF filling; `smooth > 0` writes a tanh profile
/// that many cells thick, because a geometry measurement needs an interface the
/// normal operator can actually see while the dynamics do not care how the step
/// is smeared.
void loadAlpha(const Dictionary& caseDict, volScalarField& alpha)
{
    if (!caseDict.has("alpha")) cfFatal("vof: case has no alpha{} block");
    const Dictionary& ad = caseDict.getDict("alpha");
    const polyMesh& mesh = alpha.mesh();
    const Field<Vector>& c = mesh.cellCentres();
    const Field<scalar>& v = mesh.cellVolumes();
    Field<scalar>& a = alpha.internalField();

    if (ad.has("plane"))
    {
        const Dictionary& pd = ad.getDict("plane");
        const auto p = pd.getScalarList("point");
        const auto n = pd.getScalarList("normal");
        const Vector point{p[0], p[1], p[2]};
        Vector normal{n[0], n[1], n[2]};
        if (mag(normal) == 0) cfFatal("vof: alpha plane normal is zero");
        normal = normal / mag(normal);
        for (label i = 0; i < mesh.nCells(); ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            a[u] = dot(c[u] - point, normal) < 0 ? 1.0 : 0.0;
        }
        log::info("alpha: plane through (", point.x, ", ", point.y, ", ", point.z,
                  ") normal (", normal.x, ", ", normal.y, ", ", normal.z,
                  "), phase 1 on the negative side");
        return;
    }

    if (ad.has("circle"))
    {
        const Dictionary& cd = ad.getDict("circle");
        const auto p = cd.getScalarList("center");
        const Vector center{p[0], p[1], p[2]};
        const scalar radius = cd.getScalar("radius");
        const scalar smooth = cd.getScalarOr("smooth", 0.0);
        // one representative cell edge, used only for the thickness of the
        // profile: the case is 2D by construction, so the diagonal of a cell
        // face is what the interface resolves against
        scalar dx = 0;
        for (label i = 0; i < mesh.nCells(); ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            dx = std::max(dx, std::sqrt(v[u]));
        }
        const scalar eps = smooth > 0 ? smooth * dx : 0.0;
        for (label i = 0; i < mesh.nCells(); ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            const scalar r = mag(c[u] - center);
            a[u] = eps > 0 ? 0.5 * (1.0 - std::tanh((r - radius) / eps))
                           : (r < radius ? 1.0 : 0.0);
        }
        log::info("alpha: circle centre (", center.x, ", ", center.y, ", ", center.z,
                  ") radius ", num(radius), smooth > 0 ? " smoothed" : " cell filled");
        return;
    }

    cfFatal("vof: alpha needs a plane { point; normal; } or circle { center; radius; }");
}

/// alpha at the boundary. A phase fraction has to be fed in wherever the
/// prescribed flux enters the domain, and the initial field has to agree with
/// what is fed in, or the first reports carry a boundary artefact instead of
/// the interface.
void loadAlphaPatches(const Dictionary& caseDict, const polyMesh& mesh, volScalarField& alpha)
{
    if (!caseDict.has("patches")) return;
    const Dictionary& ps = caseDict.getDict("patches");
    for (const std::string& name : ps.keys())
    {
        const Dictionary& pd = ps.getDict(name);
        if (!pd.has("alphaType")) continue;
        const label patch = mesh.patchIdByName(name);
        if (patch < 0) cfFatal("vof: patch '" + name + "' is not in the mesh");
        alpha.setPatchValue(patch, pd.getWord("alphaType"), pd.getScalarOr("alpha", 0.0));
    }
}

/// Where alpha crosses 0.5 along x on the row of cells holding the phase
/// centroid. For a plane interface carried by the flow that is the position of
/// the interface; for a closed one it is its upstream crossing, which is all a
/// checker needs in order to follow the motion.
scalar crossingX50(const polyMesh& mesh, const Field<scalar>& alpha, const Vector& centroid)
{
    const Field<Vector>& c = mesh.cellCentres();
    label row = -1;
    scalar best = std::numeric_limits<scalar>::max();
    for (label i = 0; i < mesh.nCells(); ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        // the row is the cell nearest the centroid in the plane perpendicular
        // to x, so a one-cell-thick 2D case keeps its whole interface in it
        const scalar d = std::abs(c[u].y - centroid.y) + std::abs(c[u].z - centroid.z);
        if (d < best) { best = d; row = i; }
    }
    if (row < 0) return std::numeric_limits<scalar>::quiet_NaN();
    const Vector ref = c[static_cast<std::size_t>(row)];
    std::vector<label> cells;
    for (label i = 0; i < mesh.nCells(); ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        if (std::abs(c[u].y - ref.y) < 1e-9 && std::abs(c[u].z - ref.z) < 1e-9)
            cells.push_back(i);
    }
    std::sort(cells.begin(), cells.end(), [&](label a, label b)
              { return c[static_cast<std::size_t>(a)].x < c[static_cast<std::size_t>(b)].x; });
    for (std::size_t i = 1; i < cells.size(); ++i)
    {
        const scalar a0 = alpha[static_cast<std::size_t>(cells[i - 1])];
        const scalar a1 = alpha[static_cast<std::size_t>(cells[i])];
        if ((a0 - 0.5) * (a1 - 0.5) <= 0 && a0 != a1)
        {
            const scalar x0 = c[static_cast<std::size_t>(cells[i - 1])].x;
            const scalar x1 = c[static_cast<std::size_t>(cells[i])].x;
            return x0 + (x1 - x0) * (0.5 - a0) / (a1 - a0);
        }
    }
    return std::numeric_limits<scalar>::quiet_NaN();
}

/// mean kappa over the cells that carry an interface, and how many they are
void interfaceStatistics(const volScalarField& alpha, const Schemes& schemes,
                         scalar thresholdFraction, scalar& kappaMean, int& nInterface)
{
    const Field<Vector> n = vof::interfaceNormal(alpha, schemes.grad, thresholdFraction);
    kappaMean = 0;
    nInterface = 0;
    for (const Vector& v : n)
    {
        if (mag(v) > 0) ++nInterface;
    }
    if (nInterface == 0) return;
    const Field<scalar> k = vof::interfaceCurvature(alpha.mesh(), n);
    for (std::size_t i = 0; i < static_cast<std::size_t>(n.size()); ++i)
    {
        if (mag(n[i]) > 0) kappaMean += k[i];
    }
    kappaMean /= static_cast<scalar>(nInterface);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-vof <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
        log::info("mesh: ", mesh.nCells(), " cells, ", mesh.nFaces(), " faces");

        const VelocityField vel = loadVelocity(caseDict);
        Field<scalar> phi(static_cast<std::size_t>(mesh.nFaces()), 0.0);
        scalar maxDiv = 0;
        for (label f = 0; f < mesh.nFaces(); ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            phi[u] = dot(vel.at(mesh.faceCentres()[u]), mesh.faceSf()[u]);
        }
        // the discrete continuity of the prescribed field, reported because
        // everything the transport conserves rests on it being zero
        {
            const Field<scalar> divPhi = fvc::divPhi(mesh, phi);
            for (const scalar d : divPhi) maxDiv = std::max(maxDiv, std::abs(d));
        }
        log::info("vof velocity: ",
                  vel.kind == VelocityField::Uniform ? "uniform" : "rotation",
                  " max |div(phi)| ", num(maxDiv));

        Schemes schemes = loadSchemes(caseDict);
        if (schemes.ddt != "Euler") cfFatal("vof: schemes ddt must be Euler");

        volScalarField alpha(mesh, "alpha", 0.0);
        loadAlpha(caseDict, alpha);
        loadAlphaPatches(caseDict, mesh, alpha);

        if (!caseDict.has("time")) cfFatal("vof: case has no time{} block");
        const Dictionary& td = caseDict.getDict("time");
        const scalar dt = td.getScalar("dt");
        const scalar endTime = td.getScalar("endTime");
        const int reportEvery = std::max(1, static_cast<int>(td.getScalarOr("reportEvery", 1.0)));
        const scalar threshold = td.getScalarOr("interfaceThreshold", 0.1);

        vof::Controls vc;
        vc.tol = caseDict.getScalarOr("solverTol", 1e-10);
        vc.maxIter = static_cast<int>(caseDict.getScalarOr("solverMaxIter", 2000));

        const scalar mass0 = vof::phaseVolume(mesh, alpha.internalField());
        scalar kappa0 = 0;
        int nIface0 = 0;
        interfaceStatistics(alpha, schemes, threshold, kappa0, nIface0);
        log::info("vof start: mass ", num(mass0), " dt ", num(dt),
                  " endTime ", num(endTime), " kappa ", num(kappa0),
                  " iface ", nIface0);

        // The initial geometry is the one with an analytic answer, so it is
        // what a case measures the normal and curvature operators against:
        // advection numerically diffuses the interface and tilts its normals
        // along the flow, which is a property of the transport, not of the
        // operators.
        std::string vtkInitial;
        if (caseDict.has("output"))
            vtkInitial = caseDict.getDict("output").getWordOr("vtkInitial", "");
        if (!vtkInitial.empty())
        {
            const Field<Vector> n = vof::interfaceNormal(alpha, schemes.grad, threshold);
            const Field<scalar> k = vof::interfaceCurvature(mesh, n);
            VtkWriter w(mesh);
            w.addCellScalar("alpha", alpha.internalField());
            w.addCellVector("nHat", n);
            w.addCellScalar("kappa", k);
            w.write(vtkInitial);
        }

        const int nSteps = static_cast<int>(std::llround(endTime / dt));
        SolverReport last;
        bool finite = true;
        for (int step = 1; step <= nSteps; ++step)
        {
            last = vof::advectAlpha(alpha, phi, dt, schemes, vc);
            if (step % reportEvery != 0 && step != nSteps) continue;
            const scalar mass = vof::phaseVolume(mesh, alpha.internalField());
            const Vector c = vof::phaseCentroid(mesh, alpha.internalField());
            scalar kappa = 0;
            int nIface = 0;
            interfaceStatistics(alpha, schemes, threshold, kappa, nIface);
            log::info("vof t ", num(static_cast<scalar>(step) * dt),
                      " mass ", num(mass), " massErr ", num(mass / mass0 - 1.0),
                      " min ", num(alpha.min()), " max ", num(alpha.max()),
                      " cx ", num(c.x), " cy ", num(c.y),
                      " x50 ", num(crossingX50(mesh, alpha.internalField(), c)),
                      " kappa ", num(kappa), " iface ", nIface,
                      " solver ", last.iterations, " residual ", num(last.finalResidual));
            if (!std::isfinite(mass) || !std::isfinite(alpha.min()) || !std::isfinite(alpha.max()))
            {
                log::error("vof: non-finite phase fraction at step ", step);
                finite = false;
                break;
            }
        }

        if (finite && caseDict.has("output"))
        {
            const Dictionary& od = caseDict.getDict("output");
            const std::string vtk = od.getWordOr("vtk", "");
            if (!vtk.empty())
            {
                const Field<Vector> n = vof::interfaceNormal(alpha, schemes.grad, threshold);
                const Field<scalar> k = vof::interfaceCurvature(mesh, n);
                VtkWriter w(mesh);
                w.addCellScalar("alpha", alpha.internalField());
                w.addCellVector("nHat", n);
                w.addCellScalar("kappa", k);
                w.write(vtk);
            }
        }

        if (!finite) return 1;
        // endTime 0 is a legal case: it writes the initial field and geometry
        if (nSteps > 0 && !last.converged)
        {
            log::error("vof: the transport did not converge");
            return 1;
        }
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("vof: ", e.what());
        return 1;
    }
}
