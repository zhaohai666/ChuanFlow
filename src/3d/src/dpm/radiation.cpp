#include <cmath>

#include "cf/common/Error.h"
#include "cf/fvm/Radiation.h"

namespace cf
{

void RadiationConfig::validate() const
{
    if (!std::isfinite(absorption) || !(absorption > 0))
        cfFatal("radiation: absorption kappa_a must be a finite positive number (1/m)");
    if (!std::isfinite(scattering) || scattering < 0)
        cfFatal("radiation: scattering kappa_s cannot be negative or non-finite (1/m)");
    if (!std::isfinite(mediumTemperature) || !(mediumTemperature > 0))
        cfFatal("radiation: mediumTemperature must be a finite positive temperature (K)");
    if (!std::isfinite(sbConstant) || !(sbConstant > 0))
        cfFatal("radiation: sigmaSB must be finite and positive");
    if (!std::isfinite(extinction()) || !(extinction() > 0))
        cfFatal("radiation: extinction beta = kappa_a + kappa_s must be positive");
}

scalar blackbodyIncidentRadiation(scalar temperature, scalar sigmaSB)
{
    // I_b = sigma T^4 / pi  =>  G = integral_4pi I_b dOmega = 4 pi I_b = 4 sigma T^4
    return 4.0 * sigmaSB * std::pow(temperature, 4);
}

RadiationConfig readRadiationConfig(const Dictionary& caseDict)
{
    RadiationConfig cfg;
    if (!caseDict.has("radiation")) return cfg;
    const Dictionary& rd = caseDict.getDict("radiation");
    cfg.enabled = true;
    cfg.absorption = rd.getScalarOr("absorption", 0.0);
    cfg.scattering = rd.getScalarOr("scattering", 0.0);
    cfg.mediumTemperature = rd.getScalarOr("mediumTemperature", 0.0);
    cfg.validate();
    return cfg;
}

RadiationResult solveP1(const polyMesh& mesh, volScalarField& G,
                        const RadiationConfig& cfg, const Schemes& schemes,
                        const Field<scalar>& mediumT)
{
    cfg.validate();
    const label nC = mesh.nCells();
    RadiationResult out;
    out.G.assign(nC, 0.0);
    out.qRad.assign(nC, Vector::zero());
    out.netSource.assign(nC, 0.0);

    if (static_cast<label>(mediumT.size()) != nC)
        cfFatal("radiation: medium temperature size does not match the mesh");

    const scalar sigma = cfg.sbConstant;
    const scalar kappa = cfg.absorption;
    const scalar gamma = cfg.diffusion(); // 1/(3 beta)

    // Equation assembled on the lhs:
    //   -div(gamma grad G) + kappa G = 4 kappa sigma T^4
    // i.e. the diffusion operator with sign -1 (as in the energy equation) and
    // the absorption as an implicit sp source with sign +1.
    fvMatrix eqn(mesh, viewOf(G), "P1Eqn");
    eqn.schemes() = schemes;
    eqn.schemes().relaxation = 1.0; // steady: solve to the fixed point directly
    fvm::laplacian(eqn, gamma, schemes.laplacian, -1);

    const auto& vol = mesh.cellVolumes();
    // The equation is assembled in integrated (x volume) form, the same as the
    // energy equation: the laplacian's face transmissibility gamma*oc already
    // carries the area/distance, and the source is multiplied by the cell
    // volume. fvm::sp does NOT scale by volume itself, so hand it vol*kappa.
    Field<scalar> sp(static_cast<std::size_t>(nC), 0.0);
    for (label c = 0; c < nC; ++c)
        sp[static_cast<std::size_t>(c)] = kappa * vol[static_cast<std::size_t>(c)];
    fvm::sp(eqn, sp, +1);

    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const scalar t = mediumT[u];
        // 4 kappa sigma T^4 [1/m * W/m2 = W/m3], times volume [m3] = W
        const scalar rhs = 4.0 * kappa * sigma * std::pow(t, 4);
        eqn.addRhs(c, 0, vol[u] * rhs);
    }

    const fvMatrix::Controls ctrl{1e-11, 20000, 2};
    const SolverReport rep = eqn.solve(ctrl);
    out.iterations = rep.iterations;
    out.residual = rep.finalResidual;
    out.solved = rep.converged;

    const Field<scalar>& gField = eqn.solution()[0];
    G.internalField() = gField;
    for (label c = 0; c < nC; ++c)
        out.G[static_cast<std::size_t>(c)] = gField[static_cast<std::size_t>(c)];

    // Radiative flux q_r = -gamma grad G, using the same gradient scheme the
    // diffusion operator's non-orthogonal corrector uses.
    const Field<Vector> gradG = fvc::grad(G, schemes.grad);
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        out.qRad[u] = -gamma * gradG[u];
        const scalar t = mediumT[u];
        out.netSource[u] = kappa * (gField[u] - 4.0 * sigma * std::pow(t, 4));
    }

    // Discrete balance: integral(S_rad dV) + surface_integral(q_r . n dA).
    scalar volSum = 0;
    for (label c = 0; c < nC; ++c)
        volSum += vol[static_cast<std::size_t>(c)] * out.netSource[static_cast<std::size_t>(c)];

    // Boundary flux through each boundary face, assembled with the SAME face
    // gradient the diffusion operator uses (orthoCoeff times the BC
    // linearisation), so the discrete balance closes to the linear tolerance.
    //   q_r . Sf = -gamma * oc * (G_face - G_owner)
    // with G_face - G_owner = (c-1) G_owner + k from G.linearise. For a
    // zeroGradient wall (c=1, k=0) this is exactly zero; for a fixedValue black
    // wall (c=0, k=G_wall) it is the Fourier flux at the face.
    scalar surfSum = 0;
    const auto& sf = mesh.faceSf();
    const auto& dc = mesh.deltaCoeffVec();
    for (label f = mesh.nInternalFaces(); f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const label owner = mesh.owner()[u];
        const Vector& d = dc[u];
        const scalar d2 = dot(d, d);
        if (!(d2 > 0)) continue;
        const scalar oc = dot(sf[u], d) / d2;
        scalar c = 0, k = 0;
        G.linearise(f, c, k);
        const scalar gOwner = gField[static_cast<std::size_t>(owner)];
        const scalar gFaceMinusOwner = (c - 1.0) * gOwner + k;
        surfSum += -gamma * oc * gFaceMinusOwner;
    }

    out.volumeSource = volSum;
    out.boundaryFlux = surfSum;
    out.balance = volSum + surfSum;
    return out;
}

} // namespace cf
