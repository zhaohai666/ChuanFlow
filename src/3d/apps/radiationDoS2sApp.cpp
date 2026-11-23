// chflow-3d-radiation-do: F9 grey radiation beyond P1.
//
// One binary, selected by the case `model` switch:
//   model do   -- discrete ordinates S4/S8 sweep on a uniform 1D black-wall slab
//   model s2s  -- surface-to-surface view-factor closure (transparent limit)
//   model p1   -- delegate to the shipped, untouched solveP1 for a trend check
//
// The reported RADDO lines carry 12-digit numbers; the independent checker
// rebuilds the analytic transport solution from the case constants alone and
// never trusts a number printed by this app.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/Radiation.h"
#include "cf/fvm/RadiationDoS2s.h"
#include "cf/fvm/VtkWriter.h"

using namespace cf;

namespace
{

std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}

const char* modelName(RadiationModel m)
{
    switch (m)
    {
    case RadiationModel::DO: return "do";
    case RadiationModel::S2S: return "s2s";
    case RadiationModel::P1: return "p1";
    }
    return "unknown";
}

struct Slab
{
    int nCells = 0;
    scalar length = 0;
    scalar TL = 0;
    scalar TR = 0;
};

Slab readSlab(const Dictionary& caseDict, const polyMesh& mesh)
{
    Slab s;
    const Dictionary& blk = caseDict.getDict("block");
    const std::vector<double> mn = blk.getScalarList("min");
    const std::vector<double> mx = blk.getScalarList("max");
    const std::vector<double> nc = blk.getScalarList("cells");
    if (mn.size() < 3 || mx.size() < 3 || nc.size() < 3)
        cfFatal("radiation DO: block needs min/max/cells vectors");
    s.length = mx[0] - mn[0];
    s.nCells = static_cast<int>(nc[0]);
    if (!(s.length > 0) || !(s.nCells > 0))
        cfFatal("radiation DO: slab needs positive length and cell count");
    if (static_cast<int>(nc[1]) != 1 || static_cast<int>(nc[2]) != 1 ||
        mesh.nCells() != s.nCells)
        cfFatal("radiation DO: only a 1D slab (cells y=z=1) is supported");

    if (!caseDict.has("patches")) cfFatal("radiation DO: needs patches{} with wallTemperature");
    const Dictionary& ps = caseDict.getDict("patches");
    if (!ps.has("left") || !ps.has("right"))
        cfFatal("radiation DO: slab needs patches 'left' and 'right'");
    s.TL = ps.getDict("left").getScalar("wallTemperature");
    s.TR = ps.getDict("right").getScalar("wallTemperature");
    if (!(s.TL > 0) || !(s.TR > 0))
        cfFatal("radiation DO: wallTemperature must be finite and positive");
    return s;
}

void printResult(const RadSlabResult& r, const char* model, const Slab& s)
{
    scalar gMin = r.G.empty() ? 0 : r.G[0];
    scalar gMax = gMin;
    scalar gSum = 0;
    for (scalar g : r.G)
    {
        gMin = std::min(gMin, g);
        gMax = std::max(gMax, g);
        gSum += g;
    }
    const scalar gMean = r.G.empty() ? 0 : gSum / r.G.size();
    log::info("RADDO config model ", model, " cells ", r.G.size(),
              " length ", num(s.length), " iterations ", r.iterations);
    for (std::size_t i = 0; i < r.G.size(); ++i)
    {
        const scalar x = (static_cast<scalar>(i) + 0.5) * s.length / r.G.size();
        const Vector& q = r.q[i];
        log::info("RADDO cell ", i, " x ", num(x), " G ", num(r.G[i]),
                  " q ", num(q.x), " ", num(q.y), " ", num(q.z),
                  " S ", num(r.S[i]));
    }
    log::info("RADDO t 0 model ", model,
              " Gmean ", num(gMean),
              " qLeft ", num(r.qLeft),
              " qRight ", num(r.qRight),
              " volSrc ", num(r.volSrc),
              " balance ", num(r.balance),
              " Gmin ", num(gMin),
              " Gmax ", num(gMax),
              " viewFactor ", num(r.viewFactor));
}

// P1 branch: assemble the untouched shipped solver on the same slab so the
// checker can grade a P1/DO agreement trend without a second binary.
RadSlabResult solveP1Branch(const Dictionary& caseDict, const polyMesh& mesh,
                            const DoS2sConfig& cfg, const Slab& s)
{
    RadiationConfig rp1 = readRadiationConfig(caseDict);
    if (!rp1.enabled) cfFatal("radiation DO: p1 needs a radiation{} declaration");
    volScalarField G(mesh, "G", 0.0);
    if (caseDict.has("patches"))
    {
        const Dictionary& ps = caseDict.getDict("patches");
        for (const auto& name : ps.keys())
        {
            const Dictionary& pd = ps.getDict(name);
            const label patch = mesh.patchIdByName(name);
            if (pd.has("wallTemperature"))
            {
                const scalar tw = pd.getScalar("wallTemperature");
                if (!(tw > 0)) cfFatal("radiation DO: wallTemperature must be positive");
                G.setPatchValue(patch, "fixedValue",
                                blackbodyIncidentRadiation(tw, cfg.sbConstant));
            }
            else G.setPatchValue(patch, "zeroGradient", 0.0);
        }
    }
    const Field<scalar> mediumT(static_cast<std::size_t>(mesh.nCells()),
                                rp1.mediumTemperature);
    const Schemes schemes = loadSchemes(caseDict);
    const RadiationResult rr = solveP1(mesh, G, rp1, schemes, mediumT);
    if (!rr.solved) cfFatal("radiation DO: delegated P1 solve did not converge");

    RadSlabResult out;
    const int N = s.nCells;
    out.G.assign(N, 0.0);
    out.q.assign(N, Vector::zero());
    out.S.assign(N, 0.0);
    for (int i = 0; i < N; ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        out.G[u] = rr.G[u];
        out.q[u] = rr.qRad[u];
        out.S[u] = rr.netSource[u];
    }
    // P1 wall flux is read from the boundary cell's diffusion flux; this is a
    // trend quantity only, not the sharp balance the DO branch reports.
    out.qLeft = rr.qRad.empty() ? 0 : rr.qRad[0].x;
    out.qRight = rr.qRad.empty() ? 0 : -rr.qRad[static_cast<std::size_t>(N - 1)].x;
    out.volSrc = rr.volumeSource;
    out.balance = rr.balance;
    out.iterations = rr.iterations;
    out.residual = rr.residual;
    out.converged = rr.solved;
    return out;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-3d-radiation-do <case.dict>\n";
        return 2;
    }
    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const DoS2sConfig cfg = readDoS2sConfig(caseDict);
        const polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
        const Slab slab = readSlab(caseDict, mesh);

        RadSlabResult res;
        switch (cfg.model)
        {
        case RadiationModel::DO:
            res = solveDoSlab(slab.nCells, slab.length, slab.TL, slab.TR, cfg);
            break;
        case RadiationModel::S2S:
            res = solveS2sSlab(slab.nCells, slab.length, slab.TL, slab.TR, cfg);
            break;
        case RadiationModel::P1:
            res = solveP1Branch(caseDict, mesh, cfg, slab);
            break;
        }

        printResult(res, modelName(cfg.model), slab);

        if (caseDict.has("output"))
        {
            const std::string vtk = caseDict.getDict("output").getWordOr("vtk", "");
            if (!vtk.empty())
            {
                VtkWriter w(mesh);
                w.addCellScalar("G", res.G);
                w.write(vtk);
            }
        }
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("radiation DO: ", e.what());
        return 1;
    }
}
