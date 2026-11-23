#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"
#include "cf/fvm/RadiationDoS2s.h"

namespace cf
{
namespace
{
constexpr scalar kPi = 3.141592653589793238462643383279502884;

std::vector<Ordinate> gaussLegendre(const std::vector<std::pair<scalar, scalar>>& half)
{
    std::vector<Ordinate> out;
    out.reserve(2 * half.size());
    for (const auto& p : half) out.push_back({p.first, p.second});
    for (auto it = half.rbegin(); it != half.rend(); ++it)
        out.push_back({-it->first, it->second});
    return out;
}
} // namespace

std::vector<Ordinate> snQuadrature(int order)
{
    // Gauss-Legendre nodes/mu and weights on [-1,1]. sum(w)=2 by construction;
    // the nodal values are the published roots of P_4 / P_8.
    if (order == 4)
    {
        return gaussLegendre({
            {0.8611363115940526, 0.3478548451374538},
            {0.3399810435848563, 0.6521451548625461}});
    }
    if (order == 8)
    {
        return gaussLegendre({
            {0.9602898564975363, 0.1012285362903763},
            {0.7966664774136267, 0.2223810344533745},
            {0.5255324099163290, 0.3137066458778873},
            {0.1834346424956498, 0.3626837833783620}});
    }
    cfFatal("radiation DO: quadrature order must be S4 or S8 (got S" +
            std::to_string(order) + ")");
}

scalar dos2sBlackbodyG(scalar temperature, scalar sigmaSB)
{
    return 4.0 * sigmaSB * std::pow(temperature, 4);
}

void DoS2sConfig::validate() const
{
    if (!std::isfinite(absorption) || !(absorption > 0))
        cfFatal("radiation DO: absorption kappa_a must be a finite positive number (1/m)");
    if (!std::isfinite(scattering) || scattering < 0)
        cfFatal("radiation DO: scattering kappa_s cannot be negative or non-finite (1/m)");
    if (!std::isfinite(mediumTemperature) || !(mediumTemperature > 0))
        cfFatal("radiation DO: mediumTemperature must be a finite positive temperature (K)");
    if (!std::isfinite(sbConstant) || !(sbConstant > 0))
        cfFatal("radiation DO: sigmaSB must be finite and positive");
    if (!std::isfinite(extinction()) || !(extinction() > 0))
        cfFatal("radiation DO: extinction beta = kappa_a + kappa_s must be positive");
    if (quadratureOrder != 4 && quadratureOrder != 8)
        cfFatal("radiation DO: quadrature order must be S4 or S8");
}

DoS2sConfig readDoS2sConfig(const Dictionary& caseDict)
{
    DoS2sConfig cfg;
    if (!caseDict.has("radiation"))
        cfFatal("radiation DO: needs a radiation{} declaration");

    const Dictionary& rd = caseDict.getDict("radiation");
    cfg.absorption = rd.getScalarOr("absorption", 0.0);
    cfg.scattering = rd.getScalarOr("scattering", 0.0);
    cfg.mediumTemperature = rd.getScalarOr("mediumTemperature", 0.0);
    if (rd.has("sigmaSB")) cfg.sbConstant = rd.getScalar("sigmaSB");

    // `model` is accepted at top level (documented switch) or inside the
    // radiation block; the top-level spelling wins when both are present.
    std::string model = "do";
    if (caseDict.has("model")) model = caseDict.getWord("model");
    else if (rd.has("model")) model = rd.getWord("model");
    if (model == "do") cfg.model = RadiationModel::DO;
    else if (model == "s2s") cfg.model = RadiationModel::S2S;
    else if (model == "p1") cfg.model = RadiationModel::P1;
    else cfFatal("radiation DO: unknown radiation model '" + model +
                 "'; expected do, s2s or p1");

    // Quadrature order: `quadrature S8;` (word) or `quadratureOrder 8;` (int).
    if (caseDict.has("quadratureOrder"))
        cfg.quadratureOrder = caseDict.getInt("quadratureOrder");
    else if (caseDict.has("quadrature"))
    {
        const std::string q = caseDict.getWord("quadrature");
        if (q == "S4" || q == "s4") cfg.quadratureOrder = 4;
        else if (q == "S8" || q == "s8") cfg.quadratureOrder = 8;
        else cfFatal("radiation DO: unknown quadrature '" + q + "'; expected S4 or S8");
    }
    else if (rd.has("quadratureOrder"))
        cfg.quadratureOrder = rd.getInt("quadratureOrder");
    else if (rd.has("quadrature"))
    {
        const std::string q = rd.getWord("quadrature");
        if (q == "S4" || q == "s4") cfg.quadratureOrder = 4;
        else if (q == "S8" || q == "s8") cfg.quadratureOrder = 8;
        else cfFatal("radiation DO: unknown quadrature '" + q + "'; expected S4 or S8");
    }

    cfg.validate();
    return cfg;
}

RadSlabResult solveDoSlab(int nCells, scalar length, scalar TL, scalar TR,
                          const DoS2sConfig& cfg)
{
    cfg.validate();
    if (!(nCells > 0)) cfFatal("radiation DO: slab needs a positive cell count");
    if (!(length > 0) || !std::isfinite(length))
        cfFatal("radiation DO: slab length must be finite and positive");
    if (!(TL > 0) || !(TR > 0))
        cfFatal("radiation DO: black-wall temperatures must be positive");

    const int N = nCells;
    const scalar dx = length / N;
    const scalar beta = cfg.extinction();
    const scalar ka = cfg.absorption;
    const scalar ks = cfg.scattering;
    const scalar sigma = cfg.sbConstant;
    const std::vector<Ordinate> ord = snQuadrature(cfg.quadratureOrder);
    const int M = static_cast<int>(ord.size());

    const scalar Ib = sigma * std::pow(cfg.mediumTemperature, 4) / kPi;
    const scalar IL = sigma * std::pow(TL, 4) / kPi;
    const scalar IR = sigma * std::pow(TR, 4) / kPi;

    RadSlabResult out;
    out.G.assign(N, 0.0);
    out.q.assign(N, Vector::zero());
    out.S.assign(N, 0.0);

    // H = integral_4pi I dOmega / (2 pi) = sum_m w_m I_m. It is the
    // cell-averaged angular integral; using the cell average (not the point
    // value) is what makes q(0)-q(L) = integral(S dx) a discrete identity.
    std::vector<scalar> H(static_cast<std::size_t>(N), 2.0 * Ib);
    std::vector<scalar> Hnew(static_cast<std::size_t>(N), 0.0);

    // Face intensity arrays per ordinate (N+1 faces, owner->neighbour +x).
    // Hoisted out of the iteration loop so the converged faces are available
    // for the wall-flux and cell-centre reconstruction below.
    std::vector<std::vector<scalar>> faceI(
        static_cast<std::size_t>(M), std::vector<scalar>(static_cast<std::size_t>(N + 1), 0.0));

    const int maxIter = ks > 0 ? 20000 : 1;
    scalar residual = 0;
    int iter = 0;
    for (; iter < maxIter; ++iter)
    {
        std::fill(Hnew.begin(), Hnew.end(), 0.0);

        for (int m = 0; m < M; ++m)
        {
            const scalar mu = ord[static_cast<std::size_t>(m)].mu;
            const scalar w = ord[static_cast<std::size_t>(m)].w;
            const scalar a = std::abs(mu);
            const scalar E = std::exp(-beta * dx / a);
            // Fraction of a cell mean free path: a/(beta*dx)*(1-E).
            const scalar frac = a * (1.0 - E) / (beta * dx);
            std::vector<scalar>& f = faceI[static_cast<std::size_t>(m)];

            if (mu > 0)
            {
                f[0] = IL; // black left wall emits into +x
                for (int i = 0; i < N; ++i)
                {
                    const std::size_t u = static_cast<std::size_t>(i);
                    // b_i = kappa_a I_b + (kappa_s/2) H_i  [W/m2/sr]
                    const scalar b = ka * Ib + 0.5 * ks * H[u];
                    const scalar bOverBeta = b / beta;
                    const scalar in = f[static_cast<std::size_t>(i)];
                    f[static_cast<std::size_t>(i + 1)] = in * E + bOverBeta * (1.0 - E);
                    const scalar iAvg = bOverBeta + (in - bOverBeta) * frac;
                    Hnew[u] += w * iAvg;
                }
            }
            else
            {
                f[static_cast<std::size_t>(N)] = IR; // black right wall emits into -x
                for (int i = N - 1; i >= 0; --i)
                {
                    const std::size_t u = static_cast<std::size_t>(i);
                    const scalar b = ka * Ib + 0.5 * ks * H[u];
                    const scalar bOverBeta = b / beta;
                    const scalar in = f[static_cast<std::size_t>(i + 1)];
                    f[static_cast<std::size_t>(i)] = in * E + bOverBeta * (1.0 - E);
                    const scalar iAvg = bOverBeta + (in - bOverBeta) * frac;
                    Hnew[u] += w * iAvg;
                }
            }
        }

        residual = 0;
        scalar hScale = 2.0 * Ib;
        for (int i = 0; i < N; ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            residual = std::max(residual, std::abs(Hnew[u] - H[u]) / hScale);
        }
        H.swap(Hnew);
        if (ks == 0 || residual < 1e-13) { ++iter; break; }
    }

    out.iterations = iter;
    out.residual = residual;
    out.converged = ks == 0 || residual < 1e-9;
    if (!out.converged)
        cfFatal("radiation DO: source iteration did not converge (residual " +
                std::to_string(residual) + ")");

    // Reconstruct fields from the converged face intensities.
    scalar q0 = 0, qL = 0, volSrc = 0;
    for (int i = 0; i < N; ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        out.G[u] = 2.0 * kPi * H[u];
        out.S[u] = ka * (out.G[u] - 4.0 * kPi * Ib);
        volSrc += dx * out.S[u];

        scalar qc = 0;
        for (int m = 0; m < M; ++m)
        {
            const scalar mu = ord[static_cast<std::size_t>(m)].mu;
            const scalar w = ord[static_cast<std::size_t>(m)].w;
            const scalar a = std::abs(mu);
            const scalar Eh = std::exp(-beta * dx / (2.0 * a));
            const scalar b = ka * Ib + 0.5 * ks * H[u];
            const scalar bOverBeta = b / beta;
            const scalar in = mu > 0 ? faceI[static_cast<std::size_t>(m)][u]
                                     : faceI[static_cast<std::size_t>(m)][u + 1];
            const scalar iCen = bOverBeta + (in - bOverBeta) * Eh; // exact for constant b
            qc += 2.0 * kPi * w * mu * iCen;
        }
        out.q[u] = Vector(qc, 0.0, 0.0);
    }
    for (int m = 0; m < M; ++m)
    {
        const scalar mu = ord[static_cast<std::size_t>(m)].mu;
        const scalar w = ord[static_cast<std::size_t>(m)].w;
        q0 += 2.0 * kPi * w * mu * faceI[static_cast<std::size_t>(m)][0];
        qL += 2.0 * kPi * w * mu * faceI[static_cast<std::size_t>(m)][static_cast<std::size_t>(N)];
    }

    out.qLeft = q0;
    out.qRight = -qL; // net flux into the medium from the +x (right) wall
    out.volSrc = volSrc;
    out.balance = q0 - qL - volSrc;
    return out;
}

RadSlabResult solveS2sSlab(int nCells, scalar length, scalar TL, scalar TR,
                           const DoS2sConfig& cfg)
{
    cfg.validate();
    if (!(nCells > 0) || !(length > 0))
        cfFatal("radiation S2S: invalid slab geometry");
    if (!(TL > 0) || !(TR > 0))
        cfFatal("radiation S2S: black-wall temperatures must be positive");

    const int N = nCells;
    const scalar sigma = cfg.sbConstant;
    // Closed two-surface enclosure of two equal-area, directly opposed black
    // parallel plates: F_LR = 1 by summation (F_LL=0 for a plane) and F_RL=1
    // by reciprocity (equal areas). Net exchange per unit area is sigma dT^4.
    const scalar F = 1.0;
    const scalar qNet = F * sigma * (std::pow(TL, 4) - std::pow(TR, 4));
    const scalar G = 2.0 * sigma * (std::pow(TL, 4) + std::pow(TR, 4));

    RadSlabResult out;
    out.G.assign(N, G);
    out.q.assign(N, Vector(qNet, 0.0, 0.0));
    out.S.assign(N, 0.0); // transparent medium: no volumetric participation
    out.qLeft = qNet;
    out.qRight = -qNet;
    out.volSrc = 0.0;
    out.balance = out.qLeft + out.qRight - out.volSrc;
    out.iterations = 1;
    out.converged = true;
    out.viewFactor = F;
    return out;
}

} // namespace cf
