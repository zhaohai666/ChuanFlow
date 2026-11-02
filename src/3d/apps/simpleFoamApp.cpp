#include <cmath>
#include <cstdlib>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/fvm/BlockMesh.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/SimpleSolver.h"
#include "cf/fvm/VtkWriter.h"

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

        polyMesh mesh = BlockMesh::fromDict(caseDict.getDict("block"));
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

        Schemes schemes = loadSchemes(caseDict);

        SimpleControls ctl;
        if (caseDict.has("simple"))
        {
            const Dictionary& sd = caseDict.getDict("simple");
            ctl.maxIterations = static_cast<int>(sd.getScalarOr("maxIterations", ctl.maxIterations));
            ctl.continuityTol = sd.getScalarOr("continuityTol", ctl.continuityTol);
            ctl.pcorrTol = sd.getScalarOr("pcorrTol", ctl.pcorrTol);
            ctl.relaxU = sd.getScalarOr("relaxU", ctl.relaxU);
            ctl.relaxP = sd.getScalarOr("relaxP", ctl.relaxP);
            ctl.nonOrthCorrectors = static_cast<int>(sd.getScalarOr("nonOrthCorrectors", ctl.nonOrthCorrectors));
        }

        const SimpleResult res =
            solveIncompressibleSimple(mesh, U, p, phi, nu, schemes, ctl);

        log::info("SIMPLE done: iterations=", res.iterations,
                  " continuity=", res.finalContinuityError,
                  res.converged ? " (converged)" : " (NOT converged)");

        scalar uMax = 0;
        for (const Vector& u : U.internalField()) uMax = std::max(uMax, mag(u));
        log::info("U max magnitude: ", uMax);
        log::info("p range: ", p.min(), " .. ", p.max());

        if (caseDict.has("output"))
        {
            const Dictionary& od = caseDict.getDict("output");
            const std::string vtk = od.getWordOr("vtk", "");
            if (!vtk.empty())
            {
                VtkWriter w(mesh);
                w.addCellVector("U", U.internalField());
                w.addCellScalar("p", p.internalField());
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
