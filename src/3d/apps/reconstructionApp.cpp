// chflow-reconstruction: a standalone verifier for the MUSCL face
// reconstruction operator. It does NOT solve a PDE. The case prescribes an
// analytic scalar field and a velocity field; the app samples the field on
// cell centres, rebuilds every face value with the requested reconstruction
// (none = upwind, muscl = limited second-order) and prints the exact value at
// each face centre alongside the reconstructed one. A Python checker re-derives
// both from the case alone.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/Reconstruction.h"
#include "cf/fvm/VtkWriter.h"

using namespace cf;

namespace
{

constexpr scalar kPi = 3.14159265358979323846;

void usage()
{
    std::cerr << "usage: chflow-reconstruction <case.dict>\n";
    std::exit(2);
}

/// 17 significant digits so a checker can verify to 1e-12 without rounding.
std::string num(scalar v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

std::string vec(const Vector& v)
{
    return "(" + num(v.x) + " " + num(v.y) + " " + num(v.z) + ")";
}

scalar listComp(const std::vector<double>& l, int i, const std::string& what)
{
    if (l.size() != 3) cfFatal("reconstruction: " + what + " needs 3 components");
    return l[static_cast<std::size_t>(i)];
}

/// Uniform or solid-body-rotation velocity, exactly the two fields the VOF
/// transport cases use. Only the face flux direction matters here.
struct VelocityField
{
    enum Kind { Uniform, Rotation } kind = Uniform;
    Vector u = Vector::zero();
    scalar omega = 0;
    Vector center = Vector::zero();

    Vector at(const Vector& x) const
    {
        if (kind == Uniform) return u;
        const Vector r = x - center;
        return omega * Vector{-r.y, r.x, 0.0};
    }
};

VelocityField loadVelocity(const Dictionary& caseDict)
{
    if (!caseDict.has("velocity"))
        cfFatal("reconstruction: case has no velocity{} (list or block)");
    const Dictionary::Node& node = caseDict.find("velocity");
    VelocityField v;
    if (node.kind == Dictionary::Kind::list)
    {
        v.kind = VelocityField::Uniform;
        v.u = Vector{listComp(node.numbers, 0, "velocity"),
                     listComp(node.numbers, 1, "velocity"),
                     listComp(node.numbers, 2, "velocity")};
        return v;
    }
    if (node.kind == Dictionary::Kind::block)
    {
        const std::string type = node.block->getWord("type");
        if (type == "uniform")
        {
            const auto l = node.block->getScalarList("U");
            v.kind = VelocityField::Uniform;
            v.u = Vector{listComp(l, 0, "velocity U"), listComp(l, 1, "velocity U"),
                         listComp(l, 2, "velocity U")};
            return v;
        }
        if (type == "rotation")
        {
            v.kind = VelocityField::Rotation;
            v.omega = node.block->getScalar("omega");
            const auto l = node.block->getScalarList("center");
            v.center = Vector{listComp(l, 0, "velocity center"),
                              listComp(l, 1, "velocity center"),
                              listComp(l, 2, "velocity center")};
            return v;
        }
        cfFatal("reconstruction: velocity type must be uniform or rotation");
    }
    cfFatal("reconstruction: velocity must be a list (Ux Uy Uz) or a block");
    return v;
}

/// The analytic field whose face values are verified.
struct AnalyticField
{
    enum Profile { Constant, Linear, Quadratic, Step, Sine } profile = Linear;
    scalar value = 0;                 // constant
    Vector gradient = Vector::zero(); // linear
    scalar offset = 0;                // linear / quadratic / sine baseline
    scalar coefficient = 0;           // quadratic
    Vector center = Vector::zero();   // quadratic
    int axis = 0;                     // step / sine coordinate axis
    scalar threshold = 0;             // step
    scalar low = 0, high = 1;         // step
    scalar amplitude = 1, waves = 1;  // sine
    Vector domMin = Vector::zero();   // sine: domain lower corner
    Vector domMax = Vector::zero();   // sine: domain upper corner

    scalar coord(const Vector& x) const { return x[axis]; }
    scalar domainLength() const { return domMax[axis] - domMin[axis]; }

    scalar at(const Vector& x) const
    {
        switch (profile)
        {
            case Constant:
                return value;
            case Linear:
                return dot(gradient, x) + offset;
            case Quadratic:
            {
                const Vector r = x - center;
                return coefficient * dot(r, r) + offset;
            }
            case Step:
                return coord(x) < threshold ? low : high;
            case Sine:
                return offset + amplitude
                     * std::sin(2.0 * kPi * waves * (coord(x) - domMin[axis])
                                / domainLength());
        }
        return 0;
    }
};

AnalyticField loadField(const Dictionary& caseDict, const polyMesh& mesh)
{
    if (!caseDict.has("field")) cfFatal("reconstruction: case has no field{} block");
    const Dictionary& fd = caseDict.getDict("field");
    const std::string profile = fd.getWord("profile");
    AnalyticField f;
    f.domMin = mesh.minX();
    f.domMax = mesh.maxX();

    auto readAxis = [&]() -> int
    {
        const std::string a = fd.getWordOr("axis", "x");
        if (a == "x") return 0;
        if (a == "y") return 1;
        if (a == "z") return 2;
        cfFatal("reconstruction: field axis must be x, y or z");
        return 0;
    };

    if (profile == "constant")
    {
        f.profile = AnalyticField::Constant;
        f.value = fd.getScalar("value");
    }
    else if (profile == "linear")
    {
        f.profile = AnalyticField::Linear;
        const auto g = fd.getScalarList("gradient");
        f.gradient = Vector{listComp(g, 0, "field gradient"),
                            listComp(g, 1, "field gradient"),
                            listComp(g, 2, "field gradient")};
        f.offset = fd.getScalarOr("offset", 0.0);
    }
    else if (profile == "quadratic")
    {
        f.profile = AnalyticField::Quadratic;
        f.coefficient = fd.getScalar("coefficient");
        f.offset = fd.getScalarOr("offset", 0.0);
        if (fd.has("center"))
        {
            const auto c = fd.getScalarList("center");
            f.center = Vector{listComp(c, 0, "field center"),
                              listComp(c, 1, "field center"),
                              listComp(c, 2, "field center")};
        }
    }
    else if (profile == "step")
    {
        f.profile = AnalyticField::Step;
        f.axis = readAxis();
        f.threshold = fd.getScalar("threshold");
        f.low = fd.getScalarOr("low", 0.0);
        f.high = fd.getScalarOr("high", 1.0);
        if (f.high < f.low) cfFatal("reconstruction: step high must not be below low");
    }
    else if (profile == "sine")
    {
        f.profile = AnalyticField::Sine;
        f.axis = readAxis();
        f.amplitude = fd.getScalarOr("amplitude", 1.0);
        f.waves = fd.getScalar("waves");
        f.offset = fd.getScalarOr("offset", 0.0);
        if (!(f.waves > 0)) cfFatal("reconstruction: sine waves must be positive");
        if (!(f.domainLength() > 0)) cfFatal("reconstruction: sine needs a non-degenerate domain");
    }
    else
    {
        cfFatal("reconstruction: unknown field profile '" + profile + "'");
    }
    return f;
}

const char* profileName(AnalyticField::Profile p)
{
    switch (p)
    {
        case AnalyticField::Constant:  return "constant";
        case AnalyticField::Linear:    return "linear";
        case AnalyticField::Quadratic: return "quadratic";
        case AnalyticField::Step:      return "step";
        case AnalyticField::Sine:      return "sine";
    }
    return "unknown";
}

struct Controls
{
    bool muscl = false;
    reconstruction::Limiter limiter = reconstruction::Limiter::Minmod;
};

Controls loadControls(const Dictionary& caseDict)
{
    Controls c;
    std::string mode = "none";
    std::string limiterName = "minmod";
    if (caseDict.has("schemes"))
    {
        const Dictionary& sd = caseDict.getDict("schemes");
        mode = sd.getWordOr("reconstruction", "none");
        limiterName = sd.getWordOr("limiter", "minmod");
    }
    if (mode != "none" && mode != "muscl")
        cfFatal("reconstruction: schemes reconstruction must be none or muscl");
    c.muscl = (mode == "muscl");
    c.limiter = reconstruction::parseLimiter(limiterName); // rejects unknown limiters
    return c;
}

/// For an internal face: the upwind cell, the downwind cell and the cell just
// upstream of the upwind cell (-1 when the upwind cell is at the inlet). This
/// mirrors the stencil used by musclFaceValues so the printed geometry lets a
/// checker reproduce the reconstruction exactly.
struct Stencil
{
    label up = -1, down = -1, far = -1;
    bool ghost = false;
};

Stencil faceStencil(const polyMesh& mesh, const Field<scalar>& phi, label f)
{
    const label nInternal = mesh.nInternalFaces();
    const Field<label>& own = mesh.owner();
    const Field<label>& nb = mesh.neighbour();
    const Field<Vector>& centre = mesh.cellCentres();
    const std::vector<std::vector<label>>& cellFaces = mesh.cells();
    Stencil s;
    const std::size_t u = static_cast<std::size_t>(f);
    s.up = phi[u] >= 0.0 ? own[u] : nb[u];
    s.down = phi[u] >= 0.0 ? nb[u] : own[u];

    const Vector cu = centre[static_cast<std::size_t>(s.up)];
    const Vector cd = centre[static_cast<std::size_t>(s.down)];
    const Vector ehat = (cd - cu) / mag(cd - cu);
    scalar best = -std::numeric_limits<double>::max();
    for (const label g : cellFaces[static_cast<std::size_t>(s.up)])
    {
        if (g == f) continue;
        const std::size_t v = static_cast<std::size_t>(g);
        label other = -1;
        if (own[v] == s.up)
        {
            if (g < nInternal) other = nb[v];
        }
        else if (nb[v] == s.up)
        {
            other = own[v];
        }
        if (other < 0) continue;
        const scalar score = dot(cu - centre[static_cast<std::size_t>(other)], ehat);
        if (score > 1e-12 && score > best)
        {
            best = score;
            s.far = other;
        }
    }
    s.ghost = (s.far < 0);
    return s;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) usage();

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
        const VelocityField vel = loadVelocity(caseDict);
        const AnalyticField field = loadField(caseDict, mesh);
        const Controls ctl = loadControls(caseDict);

        const label nCells = mesh.nCells();
        const label nFaces = mesh.nFaces();
        Field<scalar> cellT(static_cast<std::size_t>(nCells));
        const Field<Vector>& cc = mesh.cellCentres();
        for (label c = 0; c < nCells; ++c)
            cellT[static_cast<std::size_t>(c)] = field.at(cc[static_cast<std::size_t>(c)]);

        Field<scalar> phi(static_cast<std::size_t>(nFaces), 0.0);
        const Field<Vector>& fc = mesh.faceCentres();
        const Field<Vector>& sf = mesh.faceSf();
        for (label f = 0; f < nFaces; ++f)
            phi[static_cast<std::size_t>(f)] = dot(vel.at(fc[static_cast<std::size_t>(f)]),
                                                   sf[static_cast<std::size_t>(f)]);

        const Field<scalar> faceT = ctl.muscl
            ? reconstruction::musclFaceValues(mesh, cellT, phi, ctl.limiter)
            : reconstruction::upwindFaceValues(mesh, cellT, phi);

        log::info("reconstruction: profile=", profileName(field.profile),
                  " mode=", ctl.muscl ? "muscl" : "none",
                  " limiter=", reconstruction::limiterName(ctl.limiter),
                  " cells=", nCells, " faces=", nFaces,
                  " internal=", mesh.nInternalFaces());

        scalar sumErr = 0, maxErr = 0;
        label nSample = 0;
        for (label f = 0; f < mesh.nInternalFaces(); ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            const Stencil st = faceStencil(mesh, phi, f);
            const Vector& xf = fc[u];
            const Vector& xu = cc[static_cast<std::size_t>(st.up)];
            const Vector& xd = cc[static_cast<std::size_t>(st.down)];
            Vector xFar;
            if (st.ghost)
                xFar = 2.0 * xu - xd;
            else
                xFar = cc[static_cast<std::size_t>(st.far)];

            const scalar exact = field.at(xf);
            const scalar recon = faceT[u];
            const scalar err = std::abs(exact - recon);
            sumErr += err;
            maxErr = std::max(maxErr, err);
            ++nSample;

            std::ostringstream os;
            os << "recon f=" << f
               << " fc=" << vec(xf)
               << " U=" << vec(xu)
               << " D=" << vec(xd)
               << " F=" << vec(xFar)
               << " ghost=" << (st.ghost ? 1 : 0)
               << " exact=" << num(exact)
               << " recon=" << num(recon)
               << " err=" << num(err);
            log::info(os.str());
        }

        const scalar meanErr = nSample > 0 ? sumErr / static_cast<scalar>(nSample) : 0.0;
        std::ostringstream sum;
        sum << "reconstruction summary: profile=" << profileName(field.profile)
            << " mode=" << (ctl.muscl ? "muscl" : "none")
            << " limiter=" << reconstruction::limiterName(ctl.limiter)
            << " faces=" << nSample
            << " L1=" << num(meanErr)
            << " Linf=" << num(maxErr);
        log::info(sum.str());

        if (caseDict.has("output"))
        {
            const std::string vtk = caseDict.getDict("output").getWordOr("vtk", "");
            if (!vtk.empty())
            {
                VtkWriter w(mesh);
                w.addCellScalar("T", cellT);
                w.write(vtk);
            }
        }
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("reconstruction: ", e.what());
        return 1;
    }
}
