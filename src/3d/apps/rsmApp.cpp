// chflow-3d-rsm: second-moment (Reynolds stress) closure advanced with the
// transient incompressible PISO coupling.
//
// Opt-in case keys:
//   rsm { C1 1.8; C2eps 1.92; C1eps 1.44; sigmaR 1; sigmaEps 1.3;
//         relaxation 0.7; k0 ...; eps0 ...; }
//   profile { type linearShear; rate G; }   (same analytic seeding as the LES app)
//   probe { x y z; }                        cell whose state is reported
//
// Reports, per `piso{ reportEvery }` steps and at the end, at the probe cell:
//   the six stresses, eps, k, nut, the production P_k, plus the linearity of
//   the mean profile and the mean kinetic energy.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/turbulence/RsmModel.h"
#include "case3dTurb.h"

using namespace cf;

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr, "usage: chflow-3d-rsm <case.dict>\n");
        return 1;
    }
    const Dictionary caseDict = Dictionary::parseFile(argv[1]);

    polyMesh mesh = MeshLoader::fromCase(caseDict, argv[1]);
    const scalar nu = app3dturb::caseNu(caseDict);

    volVectorField U(mesh, "U");
    volScalarField p(mesh, "p");
    surfaceScalarField phi(mesh, "phi", 0.0);
    app3dturb::applyPatches(caseDict, mesh, U, p);

    Schemes schemes = loadSchemes(caseDict);
    schemes.ddt = "Euler";

    if (caseDict.has("profile"))
    {
        const Dictionary& pd = caseDict.getDict("profile");
        const std::string type = pd.getWord("type");
        const scalar rate = pd.getScalarOr("rate", 1.0);
        const Vector lo = mesh.minX(), hi = mesh.maxX();
        if (type == "linearShear")
        {
            for (std::size_t c = 0; c < U.internalField().size(); ++c)
            {
                U.internalField()[c] = Vector{rate * (mesh.cellCentres()[c].y - lo.y), 0, 0};
            }
        }
        else if (type == "zero")
        {
            for (auto& u : U.internalField()) u = Vector::zero();
        }
        else
        {
            cfFatal("profile: unknown type '" + type + "'");
        }
    }

    phi.correctBoundaryFlux(U);
    const RsmModel::Config cfg = RsmModel::Config::fromCase(caseDict);
    RsmModel rsm(mesh, U, phi, nu, schemes, cfg);

    // probe cell: nearest centre to `probe { x y z; }`, else the domain centre
    Vector target = 0.5 * (mesh.minX() + mesh.maxX());
    if (caseDict.has("probe"))
    {
        const auto pl = caseDict.getDict("probe").getScalarList("xyz");
        target = Vector{pl[0], pl[1], pl[2]};
    }
    std::size_t probe = 0;
    scalar best = 1e300;
    for (std::size_t c = 0; c < mesh.cellCentres().size(); ++c)
    {
        const scalar d2 = mag2(mesh.cellCentres()[c] - target);
        if (d2 < best) { best = d2; probe = c; }
    }

    scalar endTime = 1.0;
    int reportEvery = 10;
    PisoControls ctl = app3dturb::pisoControls(caseDict, endTime, reportEvery);

    log::info("rsm model cells ", mesh.nCells());

    scalar t = 0;
    int step = 0;
    auto report = [&]()
    {
        const scalar R[6] = {
            rsm.stress(0).internalField()[probe], rsm.stress(1).internalField()[probe],
            rsm.stress(2).internalField()[probe], rsm.stress(3).internalField()[probe],
            rsm.stress(4).internalField()[probe], rsm.stress(5).internalField()[probe]};
        const scalar k = rsm.k()[probe];
        std::printf("t %s R11 %s R22 %s R33 %s R12 %s R13 %s R23 %s eps %s k %s "
                    "nut %s Pk %s\n",
                    app3dturb::num(t).c_str(),
                    app3dturb::num(R[0]).c_str(), app3dturb::num(R[1]).c_str(),
                    app3dturb::num(R[2]).c_str(), app3dturb::num(R[3]).c_str(),
                    app3dturb::num(R[4]).c_str(), app3dturb::num(R[5]).c_str(),
                    app3dturb::num(rsm.epsilon()[probe]).c_str(),
                    app3dturb::num(k).c_str(), app3dturb::num(rsm.nut()[probe]).c_str(),
                    app3dturb::num(rsm.production()[probe]).c_str());

        // mean-profile linearity: least-squares slope of U.x over y and the
        // worst residual, scaled by the total shear
        scalar sumY = 0, sumU = 0, sumYY = 0, sumYU = 0;
        const std::size_t n = U.internalField().size();
        for (std::size_t c = 0; c < n; ++c)
        {
            const scalar y = mesh.cellCentres()[c].y, u = U.internalField()[c].x;
            sumY += y; sumU += u; sumYY += y * y; sumYU += y * u;
        }
        const scalar den = n * sumYY - sumY * sumY;
        const scalar G = den != 0 ? (n * sumYU - sumY * sumU) / den : 0;
        const scalar B = (sumU - G * sumY) / n;
        scalar worst = 0;
        for (std::size_t c = 0; c < n; ++c)
        {
            const scalar lin = G * mesh.cellCentres()[c].y + B;
            worst = std::max(worst, std::abs(U.internalField()[c].x - lin));
        }
        const scalar span = mesh.maxX().y - mesh.minX().y;
        scalar energy = 0, volSum = 0;
        for (std::size_t c = 0; c < n; ++c)
        {
            energy += 0.5 * mag2(U.internalField()[c]) * mesh.cellVolumes()[c];
            volSum += mesh.cellVolumes()[c];
        }
        std::printf("shear %s linResid %s E %s\n", app3dturb::num(G).c_str(),
                    app3dturb::num(worst / std::max(std::abs(G) * span, 1e-300)).c_str(),
                    app3dturb::num(energy / volSum).c_str());
    };

    report();
    while (t < endTime - 1e-12 * std::max(1.0, endTime))
    {
        const scalar dt = std::min(ctl.dt, endTime - t);
        PisoControls stepCtl = ctl;
        stepCtl.dt = dt;
        solveTransientIncompressiblePiso(mesh, U, p, phi, nu, schemes, stepCtl, &rsm);
        t += dt;
        ++step;
        if (step % reportEvery == 0) report();
    }
    report();

    if (caseDict.has("output"))
    {
        const Dictionary& od = caseDict.getDict("output");
        const std::string vtk = od.getWordOr("vtk", "");
        if (!vtk.empty())
        {
            VtkWriter w(mesh);
            w.addCellVector("U", U.internalField());
            w.addCellScalar("p", p.internalField());
            w.addCellScalar("R11", rsm.stress(0).internalField());
            w.addCellScalar("R22", rsm.stress(1).internalField());
            w.addCellScalar("R33", rsm.stress(2).internalField());
            w.addCellScalar("R12", rsm.stress(3).internalField());
            w.addCellScalar("epsilon", rsm.epsilon());
            w.addCellScalar("k", rsm.k());
            w.addCellScalar("nut", rsm.nut());
            w.write(vtk);
        }
    }
    return 0;
}
