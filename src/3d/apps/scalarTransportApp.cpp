#include <cmath>
#include <cstdlib>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Log.h"
#include "cf/fvm/BlockMesh.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/VtkWriter.h"

using namespace cf;

namespace
{

void usage()
{
    std::cerr << "usage: chflow-scalar-transport <case.dict>\n";
    std::exit(2);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) usage();

    try
    {
        const Dictionary caseDict = Dictionary::parseFile(argv[1]);

        // ---- mesh ----
        polyMesh mesh = BlockMesh::fromDict(caseDict.getDict("block"));
        log::info("mesh: ", mesh.nCells(), " cells, ", mesh.nFaces(), " faces, ",
                  mesh.nPatches(), " patches");

        // ---- velocity-driven flux ----
        const auto vlist = caseDict.getScalarList("velocity");
        const Vector U0{vlist[0], vlist[1], vlist[2]};
        surfaceScalarField phi(mesh, "phi", 0.0);
        for (label f = 0; f < mesh.nInternalFaces(); ++f)
        {
            phi.values()[static_cast<std::size_t>(f)] =
                dot(U0, mesh.faceSf()[static_cast<std::size_t>(f)]);
        }
        for (label f = mesh.nInternalFaces(); f < mesh.nFaces(); ++f)
        {
            const std::size_t u = static_cast<std::size_t>(f);
            // find patch of this boundary face
            std::string ptype;
            for (const auto& pa : mesh.patches())
            {
                if (f >= pa.start && f < pa.start + pa.size) { ptype = pa.type; break; }
            }
            scalar val = 0.0;
            if (ptype != "wall" && ptype != "symmetry" && ptype != "empty")
            {
                val = dot(U0, mesh.faceSf()[u]);
            }
            phi.values()[u] = val;
        }

        // ---- scalar field ----
        volScalarField T(mesh, "T", caseDict.getScalarOr("initialValue", 0.0));
        if (caseDict.has("patches"))
        {
            const Dictionary& ps = caseDict.getDict("patches");
            for (const auto& name : ps.keys())
            {
                const Dictionary& pd = ps.getDict(name);
                const label patch = mesh.patchIdByName(name);
                T.setPatchValue(patch, pd.getWord("type"), pd.getScalarOr("value", 0.0));
            }
        }

        // ---- equation ----
        const Schemes schemes = loadSchemes(caseDict);
        const scalar D = caseDict.getScalar("diffusivity");

        fvMatrix eqn(mesh, viewOf(T), "TEqn");
        eqn.schemes() = schemes;
        if (schemes.ddt == "Euler")
        {
            const scalar dt = caseDict.getScalar("deltaT");
            fvm::ddt(eqn, T.internalField(), dt);
        }
        fvm::div(eqn, phi.values(), schemes.div);
        fvm::laplacian(eqn, D, schemes.laplacian, -1);

        fvMatrix::Controls ctl;
        ctl.tol = caseDict.getScalarOr("solverTol", 1e-8);
        ctl.maxIter = static_cast<int>(caseDict.getScalarOr("solverMaxIter", 5000));
        ctl.nonOrth = 2;
        const SolverReport rep = eqn.solve(ctl);

        T.internalField() = eqn.solution()[0];

        log::info("T solved: ", rep.iterations, " iters, residual ", rep.finalResidual,
                  (rep.converged ? " (converged)" : " (NOT converged)"));
        log::info("T min/max: ", T.min(), " / ", T.max());

        // ---- output ----
        if (caseDict.has("output"))
        {
            const Dictionary& od = caseDict.getDict("output");
            const std::string vtk = od.getWordOr("vtk", "");
            if (!vtk.empty())
            {
                VtkWriter w(mesh);
                w.addCellScalar("T", T.internalField());
                w.write(vtk);
            }
        }

        if (caseDict.has("reportMaxCell"))
        {
            // sampling helper for regression tests: report T at cell index
            const label c = caseDict.getInt("reportMaxCell");
            log::info("T[cell ", c, "] = ", T.internalField()[static_cast<std::size_t>(c)]);
        }

        return rep.converged ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        log::error("scalarTransport: ", e.what());
        return 1;
    }
}
