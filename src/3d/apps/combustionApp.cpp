// chflow-3d-combustion: steady multispecies transport with one irreversible
// global step F + s O -> products, closed by the published Arrhenius rate
// (D8 increment 5). Temperature is a declared one-way field; chemistry does
// not feed energy back yet. Two regimes: with `fixedOxidiser` the fuel equation
// is linear and has a closed form the checker reconstructs; without it both
// species are solved by Picard sweeps.
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/Combustion.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"
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

void applyFluxBC(surfaceScalarField& phi, const polyMesh& mesh, const Vector& U0)
{
    for (label f = 0; f < mesh.nInternalFaces(); ++f)
        phi.values()[static_cast<std::size_t>(f)] =
            dot(U0, mesh.faceSf()[static_cast<std::size_t>(f)]);
    for (label f = mesh.nInternalFaces(); f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        std::string type;
        for (const auto& pa : mesh.patches())
            if (f >= pa.start && f < pa.start + pa.size) { type = pa.type; break; }
        phi.values()[u] = (type == "wall" || type == "symmetry" || type == "empty")
                              ? 0.0
                              : dot(U0, mesh.faceSf()[u]);
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-3d-combustion <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);
        const CombustionConfig cfg = readCombustionConfig(caseDict);
        if (!cfg.enabled) cfFatal("combustionApp: needs a combustion{} declaration");

        polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
        log::info("combustion mesh: ", mesh.nCells(), " cells, ", mesh.nFaces(), " faces");

        const auto vl = caseDict.getScalarList("velocity");
        const Vector U0{vl[0], vl[1], vl[2]};
        surfaceScalarField phi(mesh, "phi", 0.0);
        applyFluxBC(phi, mesh, U0);

        volScalarField Yf(mesh, "Yf", 0.0);
        volScalarField Yo(mesh, "Yo", 0.0);
        if (caseDict.has("patches"))
        {
            const Dictionary& ps = caseDict.getDict("patches");
            for (const auto& name : ps.keys())
            {
                const Dictionary& pd = ps.getDict(name);
                const label patch = mesh.patchIdByName(name);
                Yf.setPatchValue(patch, pd.has("fuelValue") ? "fixedValue" : "zeroGradient",
                                 pd.getScalarOr("fuelValue", 0.0));
                Yo.setPatchValue(patch, pd.has("oxValue") ? "fixedValue" : "zeroGradient",
                                 pd.getScalarOr("oxValue", 0.0));
            }
        }

        const bool fixedOx = cfg.fixedOxidiser >= 0;
        log::info("combustion model: k(T) ", num(cfg.rateConstant()), " 1/s, D ",
                  num(cfg.diffusivity), " s ", num(cfg.stoich), " T ",
                  num(cfg.temperature), " fixedOx ", fixedOx ? num(cfg.fixedOxidiser) : "no");

        CombustionResult res = solveCombustion(mesh, phi, Yf, Yo, cfg,
                                               loadSchemes(caseDict), fixedOx);

        for (label c = 0; c < mesh.nCells(); ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            log::info("combustion cell ", c, ": Yf ", num(res.Yf[u]), " Yo ",
                      num(res.Yo[u]), " rate ", num(res.reactionRate[u]));
        }

        // Boundary fluxes from the face BC linearisation, mirroring the P1 app.
        scalar fuelIn = 0, fuelOut = 0, oxIn = 0, oxOut = 0;
        const auto& sf = mesh.faceSf();
        const auto& dc = mesh.deltaCoeffVec();
        const auto& vol = mesh.cellVolumes();
        scalar fuelSrc = 0, oxSrc = 0;
        for (label c = 0; c < mesh.nCells(); ++c)
        {
            const std::size_t u = static_cast<std::size_t>(c);
            fuelSrc += vol[u] * res.reactionRate[u];
            oxSrc += cfg.stoich * vol[u] * res.reactionRate[u];
        }
        for (label f = mesh.nInternalFaces(); f < mesh.nFaces(); ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            const label owner = mesh.owner()[u];
            const Vector& d = dc[u];
            const scalar d2 = dot(d, d);
            if (!(d2 > 0)) continue;
            const scalar oc = dot(sf[u], d) / d2;
            // advective + diffusive flux of fuel through this face
            scalar cf = 0, kf = 0, co = 0, ko = 0;
            Yf.linearise(f, cf, kf);
            Yo.linearise(f, co, ko);
            const scalar yfOwner = res.Yf[static_cast<std::size_t>(owner)];
            const scalar yoOwner = res.Yo[static_cast<std::size_t>(owner)];
            const scalar phiF = phi.values()[u];
            const scalar fFace = cf * yfOwner + kf;
            const scalar oFace = co * yoOwner + ko;
            const scalar fuelFlux = phiF * fFace
                                    - cfg.diffusivity * oc * (fFace - yfOwner);
            const scalar oxFlux = phiF * oFace
                                  - cfg.diffusivity * oc * (oFace - yoOwner);
            if (fuelFlux < 0) fuelIn += -fuelFlux; else fuelOut += fuelFlux;
            if (oxFlux < 0) oxIn += -oxFlux; else oxOut += oxFlux;
        }
        log::info("combustion flux: fuelIn ", num(fuelIn), " fuelOut ", num(fuelOut),
                  " oxIn ", num(oxIn), " oxOut ", num(oxOut));
        log::info("combustion sink: fuel ", num(fuelSrc), " ox ", num(oxSrc));
        // With a pinned oxidiser the oxygen is an infinite reservoir: it is not
        // a solved field and its flux is intentionally not in balance. Only the
        // fuel (the solved species) obeys in - out = consumed.
        log::info("combustion balance: fuel ", num(fuelIn - fuelOut - fuelSrc),
                  fixedOx ? " ox pinned (reservoir, no balance expected)"
                          : (" ox " + num(oxIn - oxOut - oxSrc)));
        log::info("combustion summary: cells ", mesh.nCells(), " fuelIters ",
                  res.fuelIters, " oxIters ", res.oxIters, " solved ",
                  res.solved ? 1 : 0);

        // Arrhenius temperature scan: the checker grades exp(-Ea/RT) against an
        // independent math.exp at these temperatures, not against the solver.
        for (scalar fac : {0.8, 1.0, 1.25})
        {
            const scalar T = cfg.temperature * fac;
            log::info("combustion arrhenius: T ", num(T), " k ",
                      num(arrheniusRate(cfg.preExp, cfg.activationE, cfg.gasConstant, T)));
        }

        if (caseDict.has("output"))
        {
            const std::string vtk = caseDict.getDict("output").getWordOr("vtk", "");
            if (!vtk.empty())
            {
                VtkWriter w(mesh);
                w.addCellScalar("Yf", res.Yf);
                w.addCellScalar("Yo", res.Yo);
                w.write(vtk);
            }
        }
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error("combustionApp: ", e.what());
        return 1;
    }
}
