#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/SimpleSolver.h"
#include "cf/fvm/VtkWriter.h"

#include "caseEnergy.h"
#include "caseInit.h"
#include "caseTurbulence.h"

namespace
{
// log::info() stringifies with ostringstream defaults (6 digits); the p_rgh
// report is asserted against by the checker, so emit 12 significant digits.
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}
} // namespace

using namespace cf;

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chflow-simple-foam <case.dict>\n";
        std::exit(2);
    }

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);

        polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
        log::info("mesh: ", mesh.nCells(), " cells, ", mesh.nPatches(), " patches");

        const scalar nu = caseDict.getScalar("nu");

        volVectorField U(mesh, "U");
        volScalarField p(mesh, "p");
        surfaceScalarField phi(mesh, "phi", 0.0);

        // field BCs per patch:
        //   patches { inlet { U (1 0 0); p zeroGradient; } ... }
        if (caseDict.has("patches"))
        {
            const Dictionary& ps = caseDict.getDict("patches");
            for (const auto& name : ps.keys())
            {
                const Dictionary& pd = ps.getDict(name);
                const label patch = mesh.patchIdByName(name);
                if (pd.has("U"))
                {
                    const auto ul = pd.getScalarList("U");
                    const std::string type = pd.getWordOr("UType", "fixedValue");
                    U.setPatchValue(patch, type, Vector{ul[0], ul[1], ul[2]});
                }
                else
                {
                    U.setPatchValue(patch, pd.getWordOr("UType", "zeroGradient"), Vector::zero());
                }
                if (pd.has("p"))
                {
                    p.setPatchValue(patch, pd.getWordOr("pType", "fixedValue"), pd.getScalar("p"));
                }
                else
                {
                    p.setPatchValue(patch, "zeroGradient", 0.0);
                }
            }
        }

        // initial flux: boundary from inlet velocity BCs, interior zero
        phi.correctBoundaryFlux(U);
        app::seedInitialVelocity(caseDict, mesh, U, phi);

        Schemes schemes = loadSchemes(caseDict);

        // The model reads U and phi on every correctTurbulence(), so it is
        // created after the boundary conditions and the initial flux exist.
        const std::unique_ptr<RansModel> turb =
            app::makeTurbulence(caseDict, mesh, U, phi, nu, schemes);

        // The energy equation reads the flux the coupling closes and, through
        // the turbulence model, the eddy viscosity and strain rate, so it is
        // built once both exist.
        app::CaseEnergy energy =
            app::makeEnergy(caseDict, mesh, phi, nu, schemes, turb.get());

        SimpleControls ctl;
        // p_rgh mode: the dict's uniform force is the hydrostatic gravity that
        // p_rgh = p - force.x absorbs exactly; see pisoFoamApp.
        Vector grav = Vector::zero();
        bool pRgh = false;
        if (caseDict.has("simple"))
        {
            const Dictionary& sd = caseDict.getDict("simple");
            ctl.maxIterations = static_cast<int>(sd.getScalarOr("maxIterations", ctl.maxIterations));
            ctl.continuityTol = sd.getScalarOr("continuityTol", ctl.continuityTol);
            ctl.pcorrTol = sd.getScalarOr("pcorrTol", ctl.pcorrTol);
            ctl.relaxU = sd.getScalarOr("relaxU", ctl.relaxU);
            ctl.relaxP = sd.getScalarOr("relaxP", ctl.relaxP);
            ctl.pRelaxStart = sd.getScalarOr("pRelaxStart", ctl.pRelaxStart);
            ctl.nonOrthCorrectors = static_cast<int>(sd.getScalarOr("nonOrthCorrectors", ctl.nonOrthCorrectors));
            ctl.momentumTol = sd.getScalarOr("momentumTol", ctl.momentumTol);
            pRgh = sd.getWordOr("pRgh", "false") == "true";
            if (sd.has("force"))
            {
                const auto f = sd.getScalarList("force");
                ctl.bodyForce = Vector{f[0], f[1], f[2]};
            }
        }
        if (pRgh)
        {
            grav = ctl.bodyForce;
            ctl.bodyForce = Vector::zero();
            log::info("p_rgh mode: absorbing force into p_rgh = p - force.x");
        }

        const SimpleResult res = solveIncompressibleSimple(mesh, U, p, phi, nu, schemes, ctl,
                                                          turb.get(), energy.model.get());

        log::info("SIMPLE done: iterations=", res.iterations,
                  " continuity=", res.finalContinuityError,
                  " momentum=", res.finalMomentumResidual,
                  res.converged ? " (converged)" : " (NOT converged)");

        if (turb) app::logWallReport(*turb);
        if (energy.model) energy.logBalance();

        scalar uMax = 0;
        for (const Vector& u : U.internalField()) uMax = std::max(uMax, mag(u));
        log::info("U max magnitude: ", uMax);
        log::info("p range: ", p.min(), " .. ", p.max());

        Field<scalar> pTrue(static_cast<std::size_t>(mesh.nCells()), 0.0);
        if (pRgh)
        {
            // see pisoFoamApp: rest, a structure-free p_rgh and a linear
            // reconstructed p are the whole statement
            const auto& centres = mesh.cellCentres();
            scalar meanPrgh = 0;
            for (const scalar& v : p.internalField()) meanPrgh += v;
            meanPrgh /= static_cast<scalar>(mesh.nCells());
            scalar prghSpread = 0, pTrueMin = 1e300, pTrueMax = -1e300;
            for (label c = 0; c < mesh.nCells(); ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                prghSpread = std::max(prghSpread, std::abs(p.internalField()[u] - meanPrgh));
                pTrue[u] = p.internalField()[u] + dot(grav, centres[u]);
                pTrueMin = std::min(pTrueMin, pTrue[u]);
                pTrueMax = std::max(pTrueMax, pTrue[u]);
            }
            log::info("p_rgh: uMaxAll ", num(uMax),
                      " prghSpread ", num(prghSpread),
                      " pTrueRange ", num(pTrueMax - pTrueMin));
        }

        if (caseDict.has("output"))
        {
            const Dictionary& od = caseDict.getDict("output");
            const std::string vtk = od.getWordOr("vtk", "");
            if (!vtk.empty())
            {
                VtkWriter w(mesh);
                w.addCellVector("U", U.internalField());
                if (pRgh)
                {
                    w.addCellScalar("p", pTrue);
                    w.addCellScalar("p_rgh", p.internalField());
                }
                else
                {
                    w.addCellScalar("p", p.internalField());
                }
                if (turb) app::writeTurbulence(w, *turb);
                if (energy.model) energy.write(w);
                w.write(vtk);
            }
        }

        return res.converged ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        log::error("simpleFoam: ", e.what());
        return 1;
    }
}
