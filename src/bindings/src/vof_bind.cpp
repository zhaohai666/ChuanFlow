// Volume-of-Fluid transport family with PLIC reconstruction read-out.
//
// Replicates the transport branch of src/3d/apps/vofApp.cpp: a prescribed
// divergence-free velocity carries the phase fraction through an upwind
// advection step; the interface geometry (gradient-path normal/curvature and,
// when the case declares a reconstruction{} block, the PLIC planes) is copied
// out to NumPy. The momentum / surface-tension branch is not bound in this
// family; a case that declares only surfaceTension{} refuses loudly here.
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/common/Dictionary.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/multiphase/Plic.h"
#include "cf/multiphase/VoF.h"

#include "field_numpy.h"

namespace py = pybind11;

namespace
{

struct VelocityField
{
    enum Kind { Uniform, Rotation } kind = Uniform;
    cf::Vector u = cf::Vector::zero();
    double omega = 0;
    cf::Vector center = cf::Vector::zero();

    cf::Vector at(const cf::Vector& x) const
    {
        if (kind == Uniform) return u;
        const cf::Vector r = x - center;
        return omega * cf::Vector{-r.y, r.x, 0.0};
    }
};

struct VofReport
{
    double t = 0;
    double mass = 0;
    double massErr = 0;
    double minAlpha = 0;
    double maxAlpha = 0;
    double cx = 0;
    double cy = 0;
    double x50 = 0;
    double kappa = 0;
    int nIface = 0;
    int iterations = 0;
    double residual = 0;
};

struct VofResult
{
    std::shared_ptr<cf::polyMesh> mesh;
    std::unique_ptr<cf::volScalarField> alpha;
    cf::Field<cf::scalar> alphaInitial;
    cf::Field<cf::Vector> normal;
    cf::Field<cf::scalar> kappa;
    bool plicEnabled = false;
    cf::vof::Reconstruction plic;
    cf::Field<cf::scalar> kappaPlic;
    std::vector<VofReport> reports;
    double mass0 = 0;
};

VelocityField loadVelocity(const cf::Dictionary& caseDict)
{
    if (!caseDict.has("velocity"))
        throw std::runtime_error("vof: momentum/surface-tension mode is not bound in "
                                 "this family; declare a velocity{} transport case");
    const cf::Dictionary& vd = caseDict.getDict("velocity");
    VelocityField v;
    const std::string type = vd.getWord("type");
    if (type == "uniform")
    {
        const auto l = vd.getScalarList("U");
        v.kind = VelocityField::Uniform;
        v.u = cf::Vector{l[0], l[1], l[2]};
    }
    else if (type == "rotation")
    {
        v.kind = VelocityField::Rotation;
        v.omega = vd.getScalar("omega");
        const auto l = vd.getScalarList("center");
        v.center = cf::Vector{l[0], l[1], l[2]};
    }
    else
    {
        throw std::runtime_error("vof: velocity type must be uniform or rotation");
    }
    return v;
}

void loadAlpha(const cf::Dictionary& caseDict, cf::volScalarField& alpha)
{
    if (!caseDict.has("alpha"))
        throw std::runtime_error("vof: case has no alpha{} block");
    const cf::Dictionary& ad = caseDict.getDict("alpha");
    const cf::polyMesh& mesh = alpha.mesh();
    const cf::Field<cf::Vector>& c = mesh.cellCentres();
    const cf::Field<cf::scalar>& v = mesh.cellVolumes();
    cf::Field<cf::scalar>& a = alpha.internalField();

    if (ad.has("plane"))
    {
        const cf::Dictionary& pd = ad.getDict("plane");
        const auto p = pd.getScalarList("point");
        const auto n = pd.getScalarList("normal");
        const cf::Vector point{p[0], p[1], p[2]};
        cf::Vector normal{n[0], n[1], n[2]};
        if (cf::mag(normal) == 0)
            throw std::runtime_error("vof: alpha plane normal is zero");
        normal = normal / cf::mag(normal);
        for (cf::label i = 0; i < mesh.nCells(); ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            a[u] = cf::dot(c[u] - point, normal) < 0 ? 1.0 : 0.0;
        }
        return;
    }

    if (ad.has("circle") || ad.has("ellipse"))
    {
        const bool round = ad.has("circle");
        const cf::Dictionary& cd = ad.getDict(round ? "circle" : "ellipse");
        const auto p = cd.getScalarList("center");
        const cf::Vector center{p[0], p[1], p[2]};
        const double radiusA = cd.getScalar(round ? "radius" : "radiusA");
        const double radiusB = round ? radiusA : cd.getScalar("radiusB");
        const double smooth = cd.getScalarOr("smooth", 0.0);
        const std::string fill = cd.getWordOr("fill", "cell");
        if (fill == "exact")
        {
            if (!round)
                throw std::runtime_error("vof: alpha fill exact is a circle's");
            if (smooth > 0)
                throw std::runtime_error("vof: exactly filled interface has no smoothing");
            a = cf::vof::exactDiskAlpha(mesh, center, radiusA);
            return;
        }
        if (fill != "cell")
            throw std::runtime_error("vof: alpha fill must be cell or exact");
        double dx = 0;
        for (cf::label i = 0; i < mesh.nCells(); ++i)
            dx = std::max(dx, std::sqrt(v[static_cast<std::size_t>(i)]));
        const double eps = smooth > 0 ? smooth * dx : 0.0;
        for (cf::label i = 0; i < mesh.nCells(); ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            const double r = round ? cf::mag(c[u] - center)
                : std::sqrt(std::pow((c[u].x - center.x) / radiusA, 2) +
                            std::pow((c[u].y - center.y) / radiusB, 2) +
                            std::pow((c[u].z - center.z) / radiusA, 2));
            const double level = round ? radiusA : 1.0;
            a[u] = eps > 0 ? 0.5 * (1.0 - std::tanh((r - level) / eps))
                           : (r < level ? 1.0 : 0.0);
        }
        return;
    }
    throw std::runtime_error("vof: alpha needs a plane or circle/ellipse block");
}

void loadAlphaPatches(const cf::Dictionary& caseDict, const cf::polyMesh& mesh,
                      cf::volScalarField& alpha)
{
    if (!caseDict.has("patches")) return;
    const cf::Dictionary& ps = caseDict.getDict("patches");
    for (const std::string& name : ps.keys())
    {
        const cf::Dictionary& pd = ps.getDict(name);
        if (!pd.has("alphaType")) continue;
        const cf::label patch = mesh.patchIdByName(name);
        if (patch < 0)
            throw std::runtime_error("vof: patch '" + name + "' is not in the mesh");
        alpha.setPatchValue(patch, pd.getWord("alphaType"),
                            pd.getScalarOr("alpha", 0.0));
    }
}

double crossingX50(const cf::polyMesh& mesh, const cf::Field<cf::scalar>& alpha,
                   const cf::Vector& centroid)
{
    const cf::Field<cf::Vector>& c = mesh.cellCentres();
    cf::label row = -1;
    double best = std::numeric_limits<double>::max();
    for (cf::label i = 0; i < mesh.nCells(); ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        const double d = std::abs(c[u].y - centroid.y) + std::abs(c[u].z - centroid.z);
        if (d < best) { best = d; row = i; }
    }
    if (row < 0) return std::numeric_limits<double>::quiet_NaN();
    const cf::Vector ref = c[static_cast<std::size_t>(row)];
    std::vector<cf::label> cells;
    for (cf::label i = 0; i < mesh.nCells(); ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        if (std::abs(c[u].y - ref.y) < 1e-9 && std::abs(c[u].z - ref.z) < 1e-9)
            cells.push_back(i);
    }
    std::sort(cells.begin(), cells.end(), [&](cf::label a, cf::label b)
              { return c[static_cast<std::size_t>(a)].x <
                       c[static_cast<std::size_t>(b)].x; });
    for (std::size_t i = 1; i < cells.size(); ++i)
    {
        const double a0 = alpha[static_cast<std::size_t>(cells[i - 1])];
        const double a1 = alpha[static_cast<std::size_t>(cells[i])];
        if ((a0 - 0.5) * (a1 - 0.5) <= 0 && a0 != a1)
        {
            const double x0 = c[static_cast<std::size_t>(cells[i - 1])].x;
            const double x1 = c[static_cast<std::size_t>(cells[i])].x;
            return x0 + (x1 - x0) * (0.5 - a0) / (a1 - a0);
        }
    }
    return std::numeric_limits<double>::quiet_NaN();
}

cf::vof::NormalChoice parseChoice(const std::string& w)
{
    if (w == "youngs") return cf::vof::NormalChoice::Youngs;
    if (w == "centroid") return cf::vof::NormalChoice::Centroid;
    if (w == "mixed") return cf::vof::NormalChoice::Mixed;
    if (w == "fit") return cf::vof::NormalChoice::Fit;
    throw std::runtime_error("vof: reconstruction normal must be youngs, centroid, "
                             "mixed or fit");
}

std::unique_ptr<VofResult> solveVofCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
    auto out = std::make_unique<VofResult>();
    out->mesh = std::make_shared<cf::polyMesh>(
        cf::MeshLoader::fromCase(caseDict, casePath));

    const VelocityField vel = loadVelocity(caseDict);

    cf::Schemes schemes = cf::loadSchemes(caseDict);
    if (schemes.ddt != "Euler")
        throw std::runtime_error("vof: schemes ddt must be Euler");

    out->alpha = std::make_unique<cf::volScalarField>(*out->mesh, "alpha", 0.0);
    loadAlpha(caseDict, *out->alpha);
    loadAlphaPatches(caseDict, *out->mesh, *out->alpha);

    if (!caseDict.has("time"))
        throw std::runtime_error("vof: case has no time{} block");
    const cf::Dictionary& td = caseDict.getDict("time");
    const double dt = td.getScalar("dt");
    const double endTime = td.getScalar("endTime");
    const int reportEvery = std::max(1, static_cast<int>(td.getScalarOr("reportEvery", 1.0)));
    const double threshold = td.getScalarOr("interfaceThreshold", 0.1);

    bool plicEnabled = false;
    cf::vof::NormalChoice choice = cf::vof::NormalChoice::Mixed;
    double bandTolerance = 1e-6;
    bool geometricAdvect = false;
    if (caseDict.has("reconstruction"))
    {
        const cf::Dictionary& rd = caseDict.getDict("reconstruction");
        plicEnabled = true;
        choice = parseChoice(rd.getWordOr("normal", "mixed"));
        bandTolerance = rd.getScalarOr("bandTolerance", 1e-6);
        if (bandTolerance < 0.0)
            throw std::runtime_error("vof: bandTolerance must not be negative");
        const std::string tr = rd.getWordOr("transport", "alpha");
        if (tr == "alpha") geometricAdvect = false;
        else if (tr == "geometric") geometricAdvect = true;
        else throw std::runtime_error("vof: reconstruction transport must be alpha or geometric");
    }

    cf::vof::Controls vc;
    vc.tol = caseDict.getScalarOr("solverTol", 1e-10);
    vc.maxIter = static_cast<int>(caseDict.getScalarOr("solverMaxIter", 2000));

    out->plicEnabled = plicEnabled;
    out->mass0 = cf::vof::phaseVolume(*out->mesh, out->alpha->internalField());

    // The reconstruction the analytic closed form is measured against is the
    // initial field; keep that snapshot as the exposed PLIC read-out, and carry
    // a working copy `plicNow` through the advection loop (geometric transport
    // rebuilds it each step) without overwriting the snapshot.
    cf::vof::Reconstruction plicNow;
    cf::Field<cf::scalar> kappaPlicNow(static_cast<std::size_t>(out->mesh->nCells()), 0.0);
    const auto rebuildPlic = [&]()
    {
        if (!plicEnabled) return;
        plicNow = cf::vof::reconstruct(*out->alpha, schemes.grad, choice, bandTolerance);
        kappaPlicNow = cf::vof::interfaceCurvaturePlic(*out->mesh, plicNow);
    };
    rebuildPlic();
    out->alphaInitial = out->alpha->internalField();
    out->plic = plicNow;
    out->kappaPlic = kappaPlicNow;

    const int nSteps = static_cast<int>(std::llround(endTime / dt));
    cf::Field<cf::scalar> phi(static_cast<std::size_t>(out->mesh->nFaces()), 0.0);
    for (cf::label f = 0; f < out->mesh->nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        phi[u] = cf::dot(vel.at(out->mesh->faceCentres()[u]),
                         out->mesh->faceSf()[u]);
    }

    auto record = [&](double t, const cf::SolverReport& last)
    {
        VofReport row;
        row.t = t;
        row.mass = cf::vof::phaseVolume(*out->mesh, out->alpha->internalField());
        row.massErr = row.mass / out->mass0 - 1.0;
        row.minAlpha = out->alpha->min();
        row.maxAlpha = out->alpha->max();
        const cf::Vector cen = cf::vof::phaseCentroid(*out->mesh, out->alpha->internalField());
        row.cx = cen.x;
        row.cy = cen.y;
        row.x50 = crossingX50(*out->mesh, out->alpha->internalField(), cen);
        out->normal = cf::vof::interfaceNormal(*out->alpha, schemes.grad, threshold);
        out->kappa = cf::vof::interfaceCurvature(*out->mesh, out->normal);
        double kappaSum = 0;
        int nIface = 0;
        for (const cf::Vector& v : out->normal)
            if (cf::mag(v) > 0) ++nIface;
        for (std::size_t i = 0; i < static_cast<std::size_t>(out->normal.size()); ++i)
            if (cf::mag(out->normal[i]) > 0) kappaSum += out->kappa[i];
        row.kappa = nIface > 0 ? kappaSum / nIface : 0.0;
        row.nIface = nIface;
        row.iterations = last.iterations;
        row.residual = last.finalResidual;
        out->reports.push_back(row);
    };

    cf::SolverReport last;
    for (int step = 1; step <= nSteps; ++step)
    {
        if (geometricAdvect)
        {
            rebuildPlic();
            cf::vof::advectGeometric(*out->alpha, phi, dt, plicNow);
        }
        else
        {
            last = cf::vof::advectAlpha(*out->alpha, phi, dt, schemes, vc);
        }
        if (step % reportEvery != 0 && step != nSteps) continue;
        record(static_cast<double>(step) * dt, last);
    }
    if (out->reports.empty())
        record(0.0, last);
    return out;
}

} // namespace

void registerVofBindings(py::module_& m)
{
    py::class_<VofReport>(m, "VofReport",
        "One reported VOF state: phase volume and its error, min/max alpha, "
        "centroid, 0.5 crossing, band mean curvature and interface cell count.")
        .def_readonly("t", &VofReport::t)
        .def_readonly("mass", &VofReport::mass)
        .def_readonly("mass_err", &VofReport::massErr)
        .def_readonly("min_alpha", &VofReport::minAlpha)
        .def_readonly("max_alpha", &VofReport::maxAlpha)
        .def_readonly("cx", &VofReport::cx)
        .def_readonly("cy", &VofReport::cy)
        .def_readonly("x50", &VofReport::x50)
        .def_readonly("kappa", &VofReport::kappa)
        .def_readonly("n_iface", &VofReport::nIface)
        .def_readonly("iterations", &VofReport::iterations)
        .def_readonly("residual", &VofReport::residual);

    py::class_<VofResult, std::unique_ptr<VofResult>>(
        m, "VofResult",
        "Result of one VOF transport run. `alpha` is the final phase fraction;"
        " `normal`/`kappa` the gradient-path interface geometry; the `plic_*`"
        " fields are the PLIC reconstruction when the case declares one.")
        .def_property_readonly("alpha", [](const VofResult& r)
            { return cfbind::scalarFieldToNumpy(r.alpha->internalField()); })
        .def_property_readonly("alpha_initial", [](const VofResult& r)
            { return cfbind::scalarFieldToNumpy(r.alphaInitial); })
        .def_property_readonly("normal", [](const VofResult& r)
            { return cfbind::vectorFieldToNumpy(r.normal); })
        .def_property_readonly("kappa", [](const VofResult& r)
            { return cfbind::scalarFieldToNumpy(r.kappa); })
        .def_property_readonly("plic_enabled", [](const VofResult& r) { return r.plicEnabled; })
        .def_property_readonly("plic_normal", [](const VofResult& r)
            { return cfbind::vectorFieldToNumpy(r.plic.normal); })
        .def_property_readonly("plic_offset", [](const VofResult& r)
            { return cfbind::scalarFieldToNumpy(r.plic.offset); })
        .def_property_readonly("plic_area", [](const VofResult& r)
            { return cfbind::scalarFieldToNumpy(r.plic.area); })
        .def_property_readonly("plic_centroid", [](const VofResult& r)
            { return cfbind::vectorFieldToNumpy(r.plic.areaCentroid); })
        .def_property_readonly("plic_fraction", [](const VofResult& r)
            { return cfbind::scalarFieldToNumpy(r.plic.fraction); })
        .def_property_readonly("kappa_plic", [](const VofResult& r)
            { return cfbind::scalarFieldToNumpy(r.kappaPlic); })
        .def_property_readonly("band_cells", [](const VofResult& r)
            { return r.plic.bandCells; })
        .def_property_readonly("reports", [](const VofResult& r) { return r.reports; })
        .def_property_readonly("mass0", [](const VofResult& r) { return r.mass0; })
        .def("n_cells", [](const VofResult& r) { return r.mesh->nCells(); })
        .def("cell_centres",
            [](const VofResult& r)
             { return cfbind::vectorFieldToNumpy(r.mesh->cellCentres()); });

    m.def("solve_vof_case", &solveVofCase,
          py::arg("case_dict_path"),
          "Read a VOF transport case.dict and advance the phase-fraction"
          " advection; returns the final field, interface geometry and report"
          " series (PLIC planes when declared).");
}
