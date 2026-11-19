#include <algorithm>
#include <cmath>
#include <limits>

#include "cf/common/Error.h"
#include "cf/fvm/CompressibleFlow.h"

namespace cf::compressible
{

namespace
{

Field<scalar> component(const Field<Vector>& u, int comp)
{
    Field<scalar> out(u.size());
    for (std::size_t i = 0; i < u.size(); ++i)
        out[i] = u[i][comp];
    return out;
}

scalar ownerWeight(const polyMesh& mesh, label f)
{
    const auto& mf = mesh.faceCentres()[static_cast<std::size_t>(f)];
    const auto& cp = mesh.cellCentres()[static_cast<std::size_t>(mesh.owner()[static_cast<std::size_t>(f)])];
    const auto& d = mesh.deltaCoeffVec()[static_cast<std::size_t>(f)];
    const scalar d2 = dot(d, d);
    if (d2 <= 0) cfFatal("compressible: zero face distance");
    return std::min(1.0, std::max(0.0, dot(mf - cp, d) / d2));
}

/// Green-Gauss gradient of a cell field. Mirrors gradScheme.cpp's Gauss pass
/// (V*grad = sum_f Tf*Sf) but takes an explicit face-value callback so the
/// caller can supply velocity components without a FieldView.
Field<Vector> gaussGradient(const polyMesh& mesh, const Field<scalar>& data,
                            const std::function<scalar(label, const Field<scalar>&)>& faceValue)
{
    const label n = mesh.nCells();
    Field<Vector> sum(static_cast<std::size_t>(n), Vector::zero());
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const Vector flux = mesh.faceSf()[u] * faceValue(f, data);
        sum[static_cast<std::size_t>(mesh.owner()[u])] += flux;
        if (f < mesh.nInternalFaces())
            sum[static_cast<std::size_t>(mesh.neighbour()[u])] -= flux;
    }
    Field<Vector> out(static_cast<std::size_t>(n), Vector::zero());
    for (label c = 0; c < n; ++c)
        out[static_cast<std::size_t>(c)] = sum[static_cast<std::size_t>(c)] / mesh.cellVolumes()[static_cast<std::size_t>(c)];
    return out;
}

} // namespace

IdealGasEos IdealGasEos::fromCase(const Dictionary& caseDict)
{
    IdealGasEos eos;
    eos.gamma = caseDict.getDict("fluid").getScalar("gamma");
    // R lives in the opt-in eos{} block, never in fluid{}: adding it there
    // would change the declared interface of the existing Euler case.
    if (caseDict.has("eos"))
        eos.R = caseDict.getDict("eos").getScalar("R");
    if (!(eos.gamma > 1.0) || !std::isfinite(eos.gamma))
        cfFatal("compressible: gamma must be finite and > 1");
    if (eos.R < 0.0 || !std::isfinite(eos.R))
        cfFatal("compressible: R must be finite and non-negative");
    return eos;
}

scalar IdealGasEos::pressure(scalar rho, scalar rhoE, const Vector& mom) const
{
    const scalar e = rhoE / rho - 0.5 * dot(mom, mom) / (rho * rho);
    return (gamma - 1.0) * rho * e;
}

scalar IdealGasEos::soundSpeed(scalar rho, scalar p) const
{
    return std::sqrt(gamma * p / rho);
}

scalar IdealGasEos::temperature(scalar p, scalar rho) const
{
    if (R <= 0.0)
        cfFatal("compressible: temperature requires a declared gas constant R");
    return p / (rho * R);
}

Flux physicalFlux(const State& u, scalar p, scalar un, const Vector& nHat)
{
    return {u.rho * un, u.mom * un + p * nHat, (u.rhoE + p) * un};
}

scalar twoRarefactionPressure(scalar pL, scalar unL, scalar cL,
                              scalar pR, scalar unR, scalar cR, scalar gamma)
{
    const scalar k = (gamma - 1.0) / (2.0 * gamma);
    const scalar num = 0.5 * (gamma - 1.0) * (unL - unR) + cL + cR;
    const scalar den = cL * std::pow(pL, -k) + cR * std::pow(pR, -k);
    return std::pow(num / den, 1.0 / k);
}

Flux hllc(const State& uL, scalar pL, const Vector& vL, scalar cL,
          const State& uR, scalar pR, const Vector& vR, scalar cR,
          const Vector& nHat, scalar gamma)
{
    const scalar unL = dot(vL, nHat);
    const scalar unR = dot(vR, nHat);
    const scalar pStar = twoRarefactionPressure(pL, unL, cL, pR, unR, cR, gamma);
    const auto q = [pStar, gamma](scalar pK)
    {
        if (pStar <= pK) return scalar{1};
        return std::sqrt((gamma + 1.0) / (2.0 * gamma) * (pStar / pK)
                       + (gamma - 1.0) / (2.0 * gamma));
    };
    const scalar sL = unL - cL * q(pL);
    const scalar sR = unR + cR * q(pR);
    const scalar sM = (pR - pL + uL.rho * unL * (sL - unL) - uR.rho * unR * (sR - unR))
                    / (uL.rho * (sL - unL) - uR.rho * (sR - unR));

    const Flux fL = physicalFlux(uL, pL, unL, nHat);
    const Flux fR = physicalFlux(uR, pR, unR, nHat);
    const Vector tL = vL - unL * nHat;
    const Vector tR = vR - unR * nHat;
    const scalar rhoLs = uL.rho * (sL - unL) / (sL - sM);
    const scalar rhoRs = uR.rho * (sR - unR) / (sR - sM);
    const State usL{rhoLs, rhoLs * (sM * nHat + tL),
                    rhoLs * (uL.rhoE / uL.rho + (sM - unL) * (sM + pL / (uL.rho * (sL - unL))))};
    const State usR{rhoRs, rhoRs * (sM * nHat + tR),
                    rhoRs * (uR.rhoE / uR.rho + (sM - unR) * (sM + pR / (uR.rho * (sR - unR))))};

    if (sL >= 0.0) return fL;
    if (sR <= 0.0) return fR;
    const State uSel = sM >= 0.0 ? usL : usR;
    const Flux fSel = sM >= 0.0 ? fL : fR;
    const scalar sSel = sM >= 0.0 ? sL : sR;
    const State uOwn = sM >= 0.0 ? uL : uR;
    return {fSel.mass + sSel * (uSel.rho - uOwn.rho),
            fSel.mom + sSel * (uSel.mom - uOwn.mom),
            fSel.energy + sSel * (uSel.rhoE - uOwn.rhoE)};
}

TransportConfig TransportConfig::fromCase(const Dictionary& caseDict)
{
    TransportConfig cfg;
    if (caseDict.has("viscous"))
    {
        const Dictionary& v = caseDict.getDict("viscous");
        cfg.mu = v.getScalar("mu");
    }
    if (caseDict.has("heatFlux"))
    {
        const Dictionary& h = caseDict.getDict("heatFlux");
        cfg.kappa = h.getScalar("kappa");
    }
    validateTransport(cfg);
    return cfg;
}

void validateTransport(const TransportConfig& cfg)
{
    if (!std::isfinite(cfg.mu) || cfg.mu < 0.0)
        cfFatal("compressible: viscous mu must be finite and non-negative");
    if (!std::isfinite(cfg.kappa) || cfg.kappa < 0.0)
        cfFatal("compressible: heatFlux kappa must be finite and non-negative");
}

Field<Vector> cellGradient(const polyMesh& mesh, const Field<scalar>& data,
                           const std::string& gradScheme)
{
    // "linear" is the plain Gauss gradient; "corrected" is accepted as an
    // alias here because the block cases are orthogonal. Anything else is a
    // loud refusal rather than a silent fallback.
    if (gradScheme != "linear" && gradScheme != "corrected")
        cfFatal("compressible: unknown grad scheme '" + gradScheme + "'");
    const label nInt = mesh.nInternalFaces();
    return gaussGradient(mesh, data,
        [&mesh, nInt](label f, const Field<scalar>& d) -> scalar
        {
            const std::size_t u = static_cast<std::size_t>(f);
            if (f < nInt)
            {
                const scalar w = ownerWeight(mesh, f);
                return w * d[static_cast<std::size_t>(mesh.owner()[u])]
                     + (1.0 - w) * d[static_cast<std::size_t>(mesh.neighbour()[u])];
            }
            return d[static_cast<std::size_t>(mesh.owner()[u])];
        });
}

void velocityGradients(const polyMesh& mesh, const Field<Vector>& U,
                       const std::string& gradScheme, Field<Vector> gradU[3])
{
    for (int comp = 0; comp < 3; ++comp)
    {
        const Field<scalar> u = component(U, comp);
        gradU[comp] = cellGradient(mesh, u, gradScheme);
    }
}

Vector viscousStressFaceFlux(const polyMesh& mesh, const Field<Vector>* gradU,
                             scalar mu, label face)
{
    const std::size_t u = static_cast<std::size_t>(face);
    const label o = mesh.owner()[u];
    const std::size_t oc = static_cast<std::size_t>(o);
    const Vector& Sf = mesh.faceSf()[u];
    scalar divU = 0;
    for (int i = 0; i < 3; ++i) divU += gradU[i][oc][i];
    const Vector tauSf(
        mu * (gradU[0][oc].x * Sf.x + gradU[0][oc].y * Sf.y + gradU[0][oc].z * Sf.z)
          + mu * (gradU[0][oc].x * Sf.x + gradU[1][oc].x * Sf.y + gradU[2][oc].x * Sf.z)
          - (2.0 / 3.0) * mu * divU * Sf.x,
        mu * (gradU[1][oc].x * Sf.x + gradU[1][oc].y * Sf.y + gradU[1][oc].z * Sf.z)
          + mu * (gradU[0][oc].y * Sf.x + gradU[1][oc].y * Sf.y + gradU[2][oc].y * Sf.z)
          - (2.0 / 3.0) * mu * divU * Sf.y,
        mu * (gradU[2][oc].x * Sf.x + gradU[2][oc].y * Sf.y + gradU[2][oc].z * Sf.z)
          + mu * (gradU[0][oc].z * Sf.x + gradU[1][oc].z * Sf.y + gradU[2][oc].z * Sf.z)
          - (2.0 / 3.0) * mu * divU * Sf.z);
    return tauSf;
}

Field<Vector> viscousStressDivergence(const polyMesh& mesh,
                                      const Field<Vector>* gradU,
                                      scalar mu)
{
    const label n = mesh.nCells();
    Field<Vector> out(static_cast<std::size_t>(n), Vector::zero());
    if (mu <= 0.0) return out;
    for (label f = 0; f < mesh.nInternalFaces(); ++f)
    {
        const Vector flux = viscousStressFaceFlux(mesh, gradU, mu, f);
        out[static_cast<std::size_t>(mesh.owner()[static_cast<std::size_t>(f)])] += flux;
        out[static_cast<std::size_t>(mesh.neighbour()[static_cast<std::size_t>(f)])] -= flux;
    }
    // walls are free-slip / adiabatic here: no viscous or heat exchange
    for (label c = 0; c < n; ++c)
        out[static_cast<std::size_t>(c)] =
            out[static_cast<std::size_t>(c)] / mesh.cellVolumes()[static_cast<std::size_t>(c)];
    return out;
}

scalar heatFluxFace(const polyMesh& mesh, const Field<Vector>& gradT,
                    scalar kappa, label face)
{
    const std::size_t u = static_cast<std::size_t>(face);
    return -kappa * dot(gradT[static_cast<std::size_t>(mesh.owner()[u])],
                        mesh.faceSf()[u]);
}

Field<scalar> heatFluxDivergence(const polyMesh& mesh,
                                 const Field<Vector>& gradT,
                                 scalar kappa)
{
    Field<scalar> out(static_cast<std::size_t>(mesh.nCells()), 0.0);
    if (kappa <= 0.0) return out;
    for (label f = 0; f < mesh.nInternalFaces(); ++f)
    {
        const scalar flux = heatFluxFace(mesh, gradT, kappa, f);
        out[static_cast<std::size_t>(mesh.owner()[static_cast<std::size_t>(f)])] += flux;
        out[static_cast<std::size_t>(mesh.neighbour()[static_cast<std::size_t>(f)])] -= flux;
    }
    for (label c = 0; c < mesh.nCells(); ++c)
        out[static_cast<std::size_t>(c)] /= mesh.cellVolumes()[static_cast<std::size_t>(c)];
    return out;
}

} // namespace cf::compressible
