// chflow-3d-euler-euler (5.224): the two-fluid Euler-Euler kernel on a
// declared case. Per phase a velocity field, one shared pressure, a gas
// volume fraction advected by its own phase flux, an interphase drag and a
// mixture-continuity PISO-like coupling -- all on ChuanFlow's own FVM
// operators (see cf/multiphase/EulerEuler.h). The app owns the case loader
// and the per-frame report; the module owns the physics.
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
#include "cf/multiphase/EulerEuler.h"

using namespace cf;

namespace
{

std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

Vector vectorFromList(const std::vector<double>& v, const std::string& what)
{
    if (v.size() != 3)
        cfFatal("eulerEuler: " + what + " needs exactly three components, got "
                + std::to_string(v.size()));
    return Vector{v[0], v[1], v[2]};
}

eulerEuler::FluidProps loadFluid(const Dictionary& caseDict)
{
    if (!caseDict.has("fluid")) cfFatal("eulerEuler: case has no fluid{} block");
    const Dictionary& fd = caseDict.getDict("fluid");
    eulerEuler::FluidProps p;
    p.rhoG = fd.getScalar("rhoG");
    p.rhoL = fd.getScalar("rhoL");
    p.muG = fd.getScalar("muG");
    p.muL = fd.getScalar("muL");
    return p;
}

eulerEuler::DragModel loadDrag(const Dictionary& caseDict)
{
    if (!caseDict.has("drag")) cfFatal("eulerEuler: case has no drag{} block");
    const Dictionary& dd = caseDict.getDict("drag");
    eulerEuler::DragModel m;
    m.model = dd.getWord("model");
    m.beta = dd.getScalarOr("beta", 0.0);
    m.diameter = dd.getScalarOr("diameter", 0.0);
    return m;
}

/// The bounds groups. Every mesh patch must be claimed by exactly one of
/// inlet / outlet / walls, loudly: a silently-defaulted boundary is how a
/// two-fluid case starts leaking one phase.
struct Bounds
{
    std::vector<std::string> inlet, outlet, walls;
    Vector uGInlet = Vector::zero();
    Vector uLInlet = Vector::zero();
    scalar alphaGInlet = 0;
    bool hasInletU = false;
    scalar pOutlet = 0;
    bool hasWallP = false;
    scalar pWall = 0;
    std::string pWallPatch; // the single wall patch that anchors p
};

Bounds loadBounds(const Dictionary& caseDict, const polyMesh& mesh)
{
    if (!caseDict.has("bounds")) cfFatal("eulerEuler: case has no bounds{} block");
    const Dictionary& bd = caseDict.getDict("bounds");
    Bounds b;
    const auto claim = [&](const Dictionary& g, std::vector<std::string>& into)
    {
        for (const std::string& name : g.getWordList("patches"))
        {
            if (mesh.patchIdByName(name) < 0)
                cfFatal("eulerEuler: bounds patch '" + name + "' is not in the mesh");
            into.push_back(name);
        }
    };
    if (bd.has("inlet"))
    {
        const Dictionary& g = bd.getDict("inlet");
        claim(g, b.inlet);
        if (g.has("UG")) { b.uGInlet = vectorFromList(g.getScalarList("UG"), "inlet UG"); b.hasInletU = true; }
        if (g.has("UL")) b.uLInlet = vectorFromList(g.getScalarList("UL"), "inlet UL");
        if (g.has("alphaG")) b.alphaGInlet = g.getScalar("alphaG");
    }
    if (bd.has("outlet"))
    {
        const Dictionary& g = bd.getDict("outlet");
        claim(g, b.outlet);
        b.pOutlet = g.getScalarOr("p", 0.0);
    }
    if (bd.has("walls"))
    {
        const Dictionary& g = bd.getDict("walls");
        claim(g, b.walls);
        if (g.has("p"))
        {
            b.hasWallP = true;
            b.pWall = g.getScalar("p");
            b.pWallPatch = g.getWordOr("pOn", "");
            if (b.pWallPatch.empty())
                cfFatal("eulerEuler: walls p needs pOn <patch> to anchor the"
                        " pressure on one wall only");
        }
    }
    for (const polyMesh::Patch& patch : mesh.patches())
    {
        const auto in = [&](const std::vector<std::string>& v)
        { return std::find(v.begin(), v.end(), patch.name) != v.end(); };
        if (!in(b.inlet) && !in(b.outlet) && !in(b.walls))
            cfFatal("eulerEuler: patch '" + patch.name
                    + "' is claimed by no bounds group");
    }
    return b;
}

void applyBounds(const Bounds& b, polyMesh& mesh, eulerEuler::Fields& f)
{
    for (const std::string& name : b.inlet)
    {
        const label patch = mesh.patchIdByName(name);
        f.UG.setPatchValue(patch, "fixedValue", b.uGInlet);
        f.UL.setPatchValue(patch, "fixedValue", b.uLInlet);
        f.alphaG.setPatchValue(patch, "fixedValue", b.alphaGInlet);
        f.p.setPatchValue(patch, "zeroGradient", 0.0);
    }
    for (const std::string& name : b.outlet)
    {
        const label patch = mesh.patchIdByName(name);
        f.UG.setPatchValue(patch, "zeroGradient", Vector::zero());
        f.UL.setPatchValue(patch, "zeroGradient", Vector::zero());
        f.alphaG.setPatchValue(patch, "zeroGradient", 0.0);
        f.p.setPatchValue(patch, "fixedValue", b.pOutlet);
    }
    for (const std::string& name : b.walls)
    {
        const label patch = mesh.patchIdByName(name);
        f.UG.setPatchValue(patch, "wall", Vector::zero());
        f.UL.setPatchValue(patch, "wall", Vector::zero());
        f.alphaG.setPatchValue(patch, "zeroGradient", 0.0);
        const bool anchor = b.hasWallP && name == b.pWallPatch;
        f.p.setPatchValue(patch, anchor ? "fixedValue" : "zeroGradient", b.pWall);
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-3d-euler-euler <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
        log::info("mesh: ", mesh.nCells(), " cells, ", mesh.nFaces(), " faces");

        const eulerEuler::FluidProps fluid = loadFluid(caseDict);
        const eulerEuler::DragModel drag = loadDrag(caseDict);
        if (!caseDict.has("fields")) cfFatal("eulerEuler: case has no fields{} block");
        const Dictionary& fd = caseDict.getDict("fields");
        const scalar alphaG0 = fd.getScalar("alphaG0");
        Vector uG0 = Vector::zero(), uL0 = Vector::zero(), gravity = Vector::zero();
        if (fd.has("UG")) uG0 = vectorFromList(fd.getScalarList("UG"), "fields UG");
        if (fd.has("UL")) uL0 = vectorFromList(fd.getScalarList("UL"), "fields UL");
        if (fd.has("gravity")) gravity = vectorFromList(fd.getScalarList("gravity"), "fields gravity");
        if (!std::isfinite(gravity.x) || !std::isfinite(gravity.y) || !std::isfinite(gravity.z))
            cfFatal("eulerEuler: gravity must be finite");
        // declared constant pressure gradient (Pa/m): the only way an open
        // column holds its own mixture weight; zero unless declared
        Vector gradPImposed = Vector::zero();
        if (fd.has("dpdz"))
        {
            const scalar g = fd.getScalar("dpdz");
            if (!std::isfinite(g)) cfFatal("eulerEuler: fields dpdz must be finite");
            gradPImposed = Vector{0.0, 0.0, g};
        }

        if (!caseDict.has("run")) cfFatal("eulerEuler: case has no run{} block");
        const Dictionary& rd = caseDict.getDict("run");
        const int steps = static_cast<int>(rd.getScalar("steps"));
        const scalar dt = rd.getScalar("dt");
        const int reportEvery = std::max(1, static_cast<int>(rd.getScalarOr("reportEvery", 10.0)));

        eulerEuler::validateCase(fluid, drag, alphaG0, steps, dt);

        Schemes schemes = loadSchemes(caseDict);
        if (schemes.ddt != "Euler") cfFatal("eulerEuler: schemes ddt must be Euler");

        eulerEuler::Fields f(mesh, alphaG0);
        f.UG.internalField() = Field<Vector>(static_cast<std::size_t>(mesh.nCells()), uG0);
        f.UL.internalField() = Field<Vector>(static_cast<std::size_t>(mesh.nCells()), uL0);
        const Bounds bounds = loadBounds(caseDict, mesh);
        applyBounds(bounds, mesh, f);

        eulerEuler::Controls ctl;
        ctl.dt = dt;
        if (caseDict.has("solver"))
        {
            const Dictionary& sd = caseDict.getDict("solver");
            ctl.linearTol = sd.getScalarOr("linearTol", ctl.linearTol);
            ctl.maxIter = static_cast<int>(sd.getScalarOr("maxIter", ctl.maxIter));
            ctl.correctors = std::max(1, static_cast<int>(sd.getScalarOr("correctors", 2.0)));
            ctl.nonOrth = std::max(0, static_cast<int>(sd.getScalarOr("nonOrth", 2.0)));
        }
        if (!(ctl.linearTol > 0)) cfFatal("eulerEuler: solver linearTol must be positive");

        log::info("eulerEuler: rhoG ", num(fluid.rhoG), " rhoL ", num(fluid.rhoL),
                  " muG ", num(fluid.muG), " muL ", num(fluid.muL),
                  " drag ", drag.model,
                  drag.model == "constant" ? " beta " + num(drag.beta)
                                           : " diameter " + num(drag.diameter),
                  " alphaG0 ", num(alphaG0), " dt ", num(dt), " steps ", steps);

        const auto report = [&](scalar t, const eulerEuler::StepReport& r)
        {
            std::cout << "EE " << num(t)
                      << " alphaGmax " << num(r.alphaMax)
                      << " UGz " << num(r.UGz.z)
                      << " ULz " << num(r.ULz.z)
                      << " dpdz " << num(r.dpdz)
                      << " alphaGmin " << num(r.alphaMin)
                      << " massG " << num(r.massG)
                      << " massL " << num(r.massL)
                      << " cont " << num(r.continuityError) << std::endl;
        };

        // the initial state is a frame too: the per-step replay needs it
        {
            eulerEuler::StepReport r0;
            r0.alphaMin = f.alphaG.min();
            r0.alphaMax = f.alphaG.max();
            r0.massG = eulerEuler::phaseMass(mesh, f.alphaG.internalField(), fluid.rhoG);
            Field<scalar> aL(static_cast<std::size_t>(mesh.nCells()));
            for (std::size_t u = 0; u < aL.size(); ++u) aL[u] = 1.0 - f.alphaG.internalField()[u];
            r0.massL = eulerEuler::phaseMass(mesh, aL, fluid.rhoL);
            const Vector centre = 0.5 * (mesh.minX() + mesh.maxX());
            label probe = 0;
            scalar best = std::numeric_limits<double>::max();
            for (label c = 0; c < mesh.nCells(); ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                const scalar d = mag2(mesh.cellCentres()[u] - centre);
                if (d < best) { best = d; probe = c; }
            }
            r0.UGz = f.UG.internalField()[static_cast<std::size_t>(probe)];
            r0.ULz = f.UL.internalField()[static_cast<std::size_t>(probe)];
            r0.dpdz = gradPImposed.z;
            report(0.0, r0);
        }

        bool ok = true;
        for (int step = 1; step <= steps; ++step)
        {
            const eulerEuler::StepReport r = eulerEuler::advanceStep(
                mesh, f, fluid, drag, gravity, gradPImposed, schemes, ctl);
            if (step % reportEvery == 0 || step == steps)
                report(static_cast<scalar>(step) * dt, r);
            if (!std::isfinite(r.alphaMin) || !std::isfinite(r.alphaMax)
                || !std::isfinite(r.massG) || !std::isfinite(r.massL)
                || !std::isfinite(r.UGz.z) || !std::isfinite(r.ULz.z))
            {
                log::error("eulerEuler: non-finite field at step ", step);
                ok = false;
                break;
            }
            if (r.alphaMin < -1e-9 || r.alphaMax > 1.0 + 1e-9)
            {
                log::error("eulerEuler: alphaG left [0,1] at step ", step);
                ok = false;
                break;
            }
            if (!r.converged)
            {
                log::warn("eulerEuler: a linear solve did not converge at step ", step);
            }
        }

        std::string vtk;
        if (caseDict.has("output")) vtk = caseDict.getDict("output").getWordOr("vtk", "");
        if (ok && !vtk.empty())
        {
            Field<scalar> aL(static_cast<std::size_t>(mesh.nCells()));
            for (std::size_t u = 0; u < aL.size(); ++u) aL[u] = 1.0 - f.alphaG.internalField()[u];
            VtkWriter w(mesh);
            w.addCellScalar("alphaG", f.alphaG.internalField());
            w.addCellScalar("alphaL", aL);
            w.addCellVector("UG", f.UG.internalField());
            w.addCellVector("UL", f.UL.internalField());
            w.addCellScalar("p", f.p.internalField());
            w.addCellScalar("K", eulerEuler::dragCoefficient(f, fluid, drag));
            w.write(vtk);
        }

        return ok ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        log::error("eulerEuler: ", e.what());
        return 1;
    }
}
