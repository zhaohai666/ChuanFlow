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

// =====================================================================
// 5.227: 3D level-symmetric DO on a blockMesh box. The slab functions
// above are untouched; everything below is the new arm.
// =====================================================================

namespace
{

bool wallAllowsInflow(const Do3dWall& w) { return w.type != "symmetry" && w.type != "empty"; }

} // namespace

void Do3dConfig::validate() const
{
    if (order != 2 && order != 4)
        cfFatal("radiation DO 3D: direction set order must be 2 or 4 (got " +
                std::to_string(order) + ")");
    if (!std::isfinite(absorption) || absorption < 0)
        cfFatal("radiation DO 3D: absorption kappa_a cannot be negative or non-finite");
    if (!std::isfinite(scattering) || scattering < 0)
        cfFatal("radiation DO 3D: scattering kappa_s cannot be negative or non-finite");
    if (!std::isfinite(extinction()) || !(extinction() > 0))
        cfFatal("radiation DO 3D: extinction beta = kappa_a + kappa_s must be positive");
    if (!std::isfinite(sbConstant) || !(sbConstant > 0))
        cfFatal("radiation DO 3D: sigmaSB must be finite and positive");
    if (!(mediumT1 > 0) || !(mediumT2 > 0) || !std::isfinite(mediumT1) || !std::isfinite(mediumT2))
        cfFatal("radiation DO 3D: medium temperatures must be finite and positive");
    for (int i = 0; i < 6; ++i)
    {
        const Do3dWall& w = walls[i];
        if (w.type == "wall")
        {
            if (!(w.temperature > 0) || !std::isfinite(w.temperature))
                cfFatal("radiation DO 3D: wall temperature must be finite and positive");
            if (!(w.emissivity >= 0) || !(w.emissivity <= 1))
                cfFatal("radiation DO 3D: wall emissivity must lie in [0,1] (got " +
                        std::to_string(w.emissivity) + ")");
        }
        else if (w.type != "symmetry" && w.type != "empty")
            cfFatal("radiation DO 3D: unknown wall type '" + w.type +
                    "'; expected wall, symmetry or empty");
    }
}

Do3dConfig readDo3dConfig(const Dictionary& caseDict)
{
    Do3dConfig cfg;
    if (!caseDict.has("radiation3d"))
        cfFatal("radiation DO 3D: needs a radiation3d{} declaration");
    const Dictionary& rd = caseDict.getDict("radiation3d");

    if (rd.has("order")) cfg.order = rd.getInt("order");
    cfg.absorption = rd.getScalarOr("absorption", 0.0);
    cfg.scattering = rd.getScalarOr("scattering", 0.0);
    cfg.mediumT1 = rd.getScalarOr("mediumT1", 0.0);
    cfg.mediumT2 = rd.getScalarOr("mediumT2", cfg.mediumT1);
    if (rd.has("sigmaSB")) cfg.sbConstant = rd.getScalar("sigmaSB");
    cfg.sourceTol = rd.getScalarOr("sourceTol", cfg.sourceTol);
    cfg.maxIterations = static_cast<int>(rd.getScalarOr("maxIterations",
                                                        static_cast<double>(cfg.maxIterations)));
    // The walls are filled by the caller; validate() runs in solveDo3d. The
    // medium-side refusals (order, negative kappa, sigmaSB, medium T) fire
    // here so a bad radiation3d{} block is rejected before any geometry read.
    if (cfg.order != 2 && cfg.order != 4)
        cfFatal("radiation DO 3D: direction set order must be 2 or 4 (got " +
                std::to_string(cfg.order) + ")");
    if (!std::isfinite(cfg.absorption) || cfg.absorption < 0)
        cfFatal("radiation DO 3D: absorption kappa_a cannot be negative or non-finite");
    if (!std::isfinite(cfg.scattering) || cfg.scattering < 0)
        cfFatal("radiation DO 3D: scattering kappa_s cannot be negative or non-finite");
    if (!std::isfinite(cfg.sbConstant) || !(cfg.sbConstant > 0))
        cfFatal("radiation DO 3D: sigmaSB must be finite and positive");
    return cfg;
}

std::vector<Ordinate3d> levelSymOrdinates(int order)
{
    std::vector<Ordinate3d> out;
    if (order == 2)
    {
        const scalar c = 1.0 / std::sqrt(3.0);
        for (int sx = -1; sx <= 1; sx += 2)
            for (int sy = -1; sy <= 1; sy += 2)
                for (int sz = -1; sz <= 1; sz += 2)
                    out.push_back({Vector{sx * c, sy * c, sz * c}, kPi / 2.0});
        return out; // 8 dirs, sum w = 4 pi
    }
    if (order == 4)
    {
        // the closed orbit of (2,1,1)/sqrt(6): 24 dirs, equal weight pi/6
        const scalar n = 1.0 / std::sqrt(6.0);
        const int perms[3][3] = {{2, 1, 1}, {1, 2, 1}, {1, 1, 2}};
        for (const auto& p : perms)
            for (int sx = -1; sx <= 1; sx += 2)
                for (int sy = -1; sy <= 1; sy += 2)
                    for (int sz = -1; sz <= 1; sz += 2)
                        out.push_back({Vector{sx * p[0] * n, sy * p[1] * n, sz * p[2] * n},
                                       kPi / 6.0});
        return out; // 24 dirs, sum w = 4 pi
    }
    cfFatal("radiation DO 3D: direction set order must be 2 or 4 (got " +
            std::to_string(order) + ")");
}

RadBoxResult solveDo3d(int nx, int ny, int nz, scalar Lx, scalar Ly, scalar Lz,
                       const Do3dConfig& cfg)
{
    cfg.validate();
    if (!(nx > 0) || !(ny > 0) || !(nz > 0))
        cfFatal("radiation DO 3D: box needs positive cell counts");
    if (!(Lx > 0) || !(Ly > 0) || !(Lz > 0))
        cfFatal("radiation DO 3D: box needs positive lengths");

    const std::vector<Ordinate3d> ord = levelSymOrdinates(cfg.order);
    const int M = static_cast<int>(ord.size());
    const int nC = nx * ny * nz;
    const scalar dx = Lx / nx, dy = Ly / ny, dz = Lz / nz;
    const scalar beta = cfg.extinction();
    const scalar ka = cfg.absorption, ks = cfg.scattering;
    const scalar sigma = cfg.sbConstant;

    // mirror direction indices for grey reflection at each wall normal
    auto findDir = [&](const Vector& s)
    {
        for (int m = 0; m < M; ++m)
            if (std::abs(ord[m].s.x - s.x) < 1e-12 && std::abs(ord[m].s.y - s.y) < 1e-12 &&
                std::abs(ord[m].s.z - s.z) < 1e-12)
                return m;
        return -1;
    };
    std::vector<int> mx(M), my(M), mz(M);
    for (int m = 0; m < M; ++m)
    {
        mx[m] = findDir(Vector{-ord[m].s.x, ord[m].s.y, ord[m].s.z});
        my[m] = findDir(Vector{ord[m].s.x, -ord[m].s.y, ord[m].s.z});
        mz[m] = findDir(Vector{ord[m].s.x, ord[m].s.y, -ord[m].s.z});
    }

    // blackbody intensity per cell: linear in x between T1 and T2
    Field<scalar> Ib(static_cast<std::size_t>(nC), 0.0);
    for (int i = 0; i < nx; ++i)
    {
        const scalar fx = (i + 0.5) / nx;
        const scalar T = cfg.mediumT1 + (cfg.mediumT2 - cfg.mediumT1) * fx;
        const scalar ib = sigma * T * T * T * T / kPi;
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
                Ib[static_cast<std::size_t>(i * ny * nz + j * nz + k)] = ib;
    }
    scalar IbWall[6];
    for (int w = 0; w < 6; ++w)
    {
        const scalar T = cfg.walls[w].temperature;
        IbWall[w] = T > 0 ? sigma * T * T * T * T / kPi : 0.0;
    }

    // I[m][cell], H = sum w I (the scattering source), previous iterate for
    // the lagged wall reflection.
    std::vector<Field<scalar>> I(static_cast<std::size_t>(M));
    for (int m = 0; m < M; ++m)
        I[static_cast<std::size_t>(m)] = Field<scalar>(static_cast<std::size_t>(nC), 0.0);
    Field<scalar> H(static_cast<std::size_t>(nC), 0.0);

    const auto cellAt = [&](int i, int j, int k)
    { return static_cast<std::size_t>(i * ny * nz + j * nz + k); };

    RadBoxResult out;
    out.nx = nx; out.ny = ny; out.nz = nz;
    out.directions = M;

    int iter = 0;
    scalar residual = 0;
    for (; iter < cfg.maxIterations; ++iter)
    {
        std::fill(H.begin(), H.end(), 0.0);
        for (int m = 0; m < M; ++m)
            for (int c = 0; c < nC; ++c)
                H[static_cast<std::size_t>(c)] += ord[static_cast<std::size_t>(m)].w
                                                * I[static_cast<std::size_t>(m)][static_cast<std::size_t>(c)];

        scalar maxChange = 0;
        for (int m = 0; m < M; ++m)
        {
            const Vector s = ord[static_cast<std::size_t>(m)].s;
            Field<scalar>& Im = I[static_cast<std::size_t>(m)];
            const scalar ax = std::abs(s.x) / dx;
            const scalar ay = std::abs(s.y) / dy;
            const scalar az = std::abs(s.z) / dz;
            const bool symX0 = !wallAllowsInflow(cfg.walls[0]);
            const bool symX1 = !wallAllowsInflow(cfg.walls[1]);
            const bool symY0 = !wallAllowsInflow(cfg.walls[2]);
            const bool symY1 = !wallAllowsInflow(cfg.walls[3]);
            const bool symZ0 = !wallAllowsInflow(cfg.walls[4]);
            const bool symZ1 = !wallAllowsInflow(cfg.walls[5]);

            const auto iStart = s.x > 0 ? 0 : nx - 1, iStep = s.x > 0 ? 1 : -1;
            const auto jStart = s.y > 0 ? 0 : ny - 1, jStep = s.y > 0 ? 1 : -1;
            const auto kStart = s.z > 0 ? 0 : nz - 1, kStep = s.z > 0 ? 1 : -1;

            for (int i = iStart; s.x > 0 ? i < nx : i >= 0; i += iStep)
                for (int j = jStart; s.y > 0 ? j < ny : j >= 0; j += jStep)
                    for (int k = kStart; s.z > 0 ? k < nz : k >= 0; k += kStep)
                    {
                        const std::size_t c = cellAt(i, j, k);
                        scalar num = 0, den = beta;

                        // x axis
                        if (s.x > 0 && i == 0)
                        {
                            if (!symX0)
                            {
                                const scalar refl = I[static_cast<std::size_t>(mx[m])][cellAt(0, j, k)];
                                num += ax * (cfg.walls[0].emissivity * IbWall[0]
                                           + (1.0 - cfg.walls[0].emissivity) * refl);
                                den += ax;
                            }
                        }
                        else if (s.x < 0 && i == nx - 1)
                        {
                            if (!symX1)
                            {
                                const scalar refl = I[static_cast<std::size_t>(mx[m])][cellAt(nx - 1, j, k)];
                                num += ax * (cfg.walls[1].emissivity * IbWall[1]
                                           + (1.0 - cfg.walls[1].emissivity) * refl);
                                den += ax;
                            }
                        }
                        else
                        {
                            num += ax * Im[cellAt(i + (s.x > 0 ? -1 : 1), j, k)];
                            den += ax;
                        }
                        // y axis
                        if (s.y > 0 && j == 0)
                        {
                            if (!symY0)
                            {
                                const scalar refl = I[static_cast<std::size_t>(my[m])][cellAt(i, 0, k)];
                                num += ay * (cfg.walls[2].emissivity * IbWall[2]
                                           + (1.0 - cfg.walls[2].emissivity) * refl);
                                den += ay;
                            }
                        }
                        else if (s.y < 0 && j == ny - 1)
                        {
                            if (!symY1)
                            {
                                const scalar refl = I[static_cast<std::size_t>(my[m])][cellAt(i, ny - 1, k)];
                                num += ay * (cfg.walls[3].emissivity * IbWall[3]
                                           + (1.0 - cfg.walls[3].emissivity) * refl);
                                den += ay;
                            }
                        }
                        else
                        {
                            num += ay * Im[cellAt(i, j + (s.y > 0 ? -1 : 1), k)];
                            den += ay;
                        }
                        // z axis
                        if (s.z > 0 && k == 0)
                        {
                            if (!symZ0)
                            {
                                const scalar refl = I[static_cast<std::size_t>(mz[m])][cellAt(i, j, 0)];
                                num += az * (cfg.walls[4].emissivity * IbWall[4]
                                           + (1.0 - cfg.walls[4].emissivity) * refl);
                                den += az;
                            }
                        }
                        else if (s.z < 0 && k == nz - 1)
                        {
                            if (!symZ1)
                            {
                                const scalar refl = I[static_cast<std::size_t>(mz[m])][cellAt(i, j, nz - 1)];
                                num += az * (cfg.walls[5].emissivity * IbWall[5]
                                           + (1.0 - cfg.walls[5].emissivity) * refl);
                                den += az;
                            }
                        }
                        else
                        {
                            num += az * Im[cellAt(i, j, k + (s.z > 0 ? -1 : 1))];
                            den += az;
                        }

                        const scalar b = ka * Ib[c] + ks * H[c] / (4.0 * kPi);
                        const scalar newI = (num + b) / den;
                        maxChange = std::max(maxChange, std::abs(newI - Im[c]) / std::max(Ib[c], 1e-300));
                        Im[c] = newI;
                    }
        }

        residual = maxChange;
        if (residual < cfg.sourceTol) { ++iter; break; }
    }

    out.iterations = iter;
    out.residual = residual;
    out.converged = residual < cfg.sourceTol * 100.0;
    if (!out.converged)
        cfFatal("radiation DO 3D: source iteration did not converge (residual " +
                std::to_string(residual) + ")");

    // final scattering sum and the reported fields
    std::fill(H.begin(), H.end(), 0.0);
    for (int m = 0; m < M; ++m)
        for (int c = 0; c < nC; ++c)
            H[static_cast<std::size_t>(c)] += ord[static_cast<std::size_t>(m)].w
                                            * I[static_cast<std::size_t>(m)][static_cast<std::size_t>(c)];

    out.G.assign(static_cast<std::size_t>(nC), 0.0);
    out.q.assign(static_cast<std::size_t>(nC), Vector::zero());
    out.S.assign(static_cast<std::size_t>(nC), 0.0);
    for (int c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        out.G[u] = H[u];
        Vector qq = Vector::zero();
        for (int m = 0; m < M; ++m)
            qq += ord[static_cast<std::size_t>(m)].s * (ord[static_cast<std::size_t>(m)].w * I[static_cast<std::size_t>(m)][u]);
        out.q[u] = qq;
        out.S[u] = ka * (out.G[u] - 4.0 * kPi * Ib[u]);
    }

    // wall net fluxes, INTO the medium (i.e. along the inward normal)
    const auto faceIn = [&](int axisWall, int m, int i, int j, int k) -> scalar
    {
        const Vector s = ord[static_cast<std::size_t>(m)].s;
        const Do3dWall& w = cfg.walls[axisWall];
        const bool inward = (axisWall == 0 && s.x > 0) || (axisWall == 1 && s.x < 0)
                         || (axisWall == 2 && s.y > 0) || (axisWall == 3 && s.y < 0)
                         || (axisWall == 4 && s.z > 0) || (axisWall == 5 && s.z < 0);
        if (!inward)
            return I[static_cast<std::size_t>(m)][cellAt(i, j, k)]; // leaving face = upwind cell
        const int mirror = axisWall <= 1 ? mx[m] : axisWall <= 3 ? my[m] : mz[m];
        const scalar refl = I[static_cast<std::size_t>(mirror)][cellAt(i, j, k)];
        return w.emissivity * IbWall[axisWall] + (1.0 - w.emissivity) * refl;
    };
    // net flux into the medium per wall: sum w (n_in . s) I_face. A
    // symmetry/empty face returns the outgoing intensity unchanged, so its
    // net is exactly zero and it is skipped.
    out.qX0 = out.qX1 = out.qY0 = out.qY1 = out.qZ0 = out.qZ1 = 0;
    const bool open[6] = {wallAllowsInflow(cfg.walls[0]), wallAllowsInflow(cfg.walls[1]),
                          wallAllowsInflow(cfg.walls[2]), wallAllowsInflow(cfg.walls[3]),
                          wallAllowsInflow(cfg.walls[4]), wallAllowsInflow(cfg.walls[5])};
    for (int m = 0; m < M; ++m)
    {
        const Vector s = ord[static_cast<std::size_t>(m)].s;
        const scalar w = ord[static_cast<std::size_t>(m)].w;
        if (open[0] || open[1])
            for (int j = 0; j < ny; ++j)
                for (int k = 0; k < nz; ++k)
                {
                    if (open[0]) out.qX0 += w * s.x * faceIn(0, m, 0, j, k);
                    if (open[1]) out.qX1 -= w * s.x * faceIn(1, m, nx - 1, j, k);
                }
        if (open[2] || open[3])
            for (int i = 0; i < nx; ++i)
                for (int k = 0; k < nz; ++k)
                {
                    if (open[2]) out.qY0 += w * s.y * faceIn(2, m, i, 0, k);
                    if (open[3]) out.qY1 -= w * s.y * faceIn(3, m, i, ny - 1, k);
                }
        if (open[4] || open[5])
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                {
                    if (open[4]) out.qZ0 += w * s.z * faceIn(4, m, i, j, 0);
                    if (open[5]) out.qZ1 -= w * s.z * faceIn(5, m, i, j, nz - 1);
                }
    }

    // Energy bookkeeping in watts: the six wall fluxes are per unit area of
    // their own wall, the volumetric source integrates over the box.
    const scalar vol = Lx * Ly * Lz;
    scalar srcTotal = 0;
    for (int c = 0; c < nC; ++c)
        srcTotal += out.S[static_cast<std::size_t>(c)] * (vol / nC);
    out.volSrc = srcTotal / vol; // volume-averaged source [W/m3]
    const scalar wallW = out.qX0 * Ly * Lz + out.qX1 * Ly * Lz
                       + out.qY0 * Lx * Lz + out.qY1 * Lx * Lz
                       + out.qZ0 * Lx * Ly + out.qZ1 * Lx * Ly;
    out.balance = wallW - srcTotal; // ~0 by the discrete RTE identity
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
