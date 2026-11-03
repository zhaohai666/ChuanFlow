// Unit tests for the RANS turbulence layer: the law of the wall and its
// inversion, the SST blend functions, the eddy-viscosity formulas, the
// near-wall boundary values the wall model prescribes, the equilibrium state
// those values impose on the wall-adjacent cell, and the homogeneous decay of
// both models against their closed-form solutions.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/fvm/BlockMesh.h"
#include "cf/fvm/Fields.h"
#include "cf/turbulence/KEpsilon.h"
#include "cf/turbulence/KOmegaSST.h"
#include "cf/turbulence/TurbulenceModel.h"
#include "cf/turbulence/WallFunction.h"

using namespace cf;

static int failures = 0;

#define CHECK(cond, msg)                                                        \
    do {                                                                        \
        if (!(cond)) { ++failures; std::cout << "FAIL: " << (msg) << "\n"; }     \
        else          { std::cout << "ok:   " << (msg) << "\n"; }                \
    } while (0)

static bool near(double a, double b, double tol = 1e-10)
{
    return std::abs(a - b) < tol;
}
static bool relNear(double a, double b, double tol)
{
    return std::abs(a - b) < tol * std::max(1.0, std::max(std::abs(a), std::abs(b)));
}

namespace
{

struct Box
{
    polyMesh mesh;
    volVectorField U;
    surfaceScalarField phi;

    /// walls=true: inlet/outlet in x, a no-slip wall at ymin and a symmetry
    /// plane at ymax. walls=false: every side is a symmetry plane, so nothing
    /// touches the domain and a uniform field stays uniform.
    explicit Box(bool walls, const Vector& size = Vector{0.2, 0.02, 0.02},
                 label nx = 2, label ny = 5, label nz = 1, bool topIsWall = false)
        : mesh(makeBox(walls, size, nx, ny, nz, topIsWall)), U(mesh, "U"), phi(mesh, "phi", 0.0)
    {
        for (const polyMesh::Patch& p : mesh.patches())
            U.setPatchValue(mesh.patchIdByName(p.name), "zeroGradient", Vector::zero());
    }

    static polyMesh makeBox(bool walls, const Vector& size, label nx, label ny, label nz,
                             bool topIsWall)
    {
        BlockSides s;
        if (walls)
        {
            s.xmin = {"inlet", "inlet"};
            s.xmax = {"outlet", "outlet"};
            s.ymin = {"bottom", "wall"};
            s.ymax = {"top", topIsWall ? "wall" : "symmetry"};
        }
        else
        {
            s.xmin = {"sx0", "symmetry"};
            s.xmax = {"sx1", "symmetry"};
            s.ymin = {"sy0", "symmetry"};
            s.ymax = {"sy1", "symmetry"};
        }
        s.zmin = {"sz0", "symmetry"};
        s.zmax = {"sz1", "symmetry"};
        return BlockMesh::create(Vector::zero(), size, nx, ny, nz, s);
    }

    /// y of the cell centres, row by row (the mesh is a single block)
    scalar yOf(label cell) const
    {
        return mesh.cellCentres()[static_cast<std::size_t>(cell)].y;
    }
};

TurbulenceConfig cfg(const std::string& model, scalar relaxation = 1.0)
{
    TurbulenceConfig c;
    c.model = model;
    c.intensity = 0.05;
    c.lengthScale = 0.005;
    c.relaxation = relaxation;
    return c;
}

} // namespace

int main()
{
    constexpr scalar NU_AIR = 1.5e-5;

    // ---------------- law of the wall ----------------
    {
        CHECK(lawOfWall(0.0) == 0.0 && lawOfWall(-1.0) == 0.0, "u+ is zero off the wall");

        // viscous sublayer: u+ -> y+
        const scalar yLo = 1e-4;
        CHECK(relNear(lawOfWall(yLo), yLo, 1e-3), "u+ = y+ in the viscous sublayer");

        // log layer: u+ -> ln(E y+)/kappa
        const scalar yHi = 1e4;
        const scalar logLaw = std::log(turb::logE * yHi) / turb::kappa;
        CHECK(relNear(lawOfWall(yHi), logLaw, 5e-3), "u+ tends to the log law");

        bool monotone = true;
        scalar prev = -1;
        for (double yp = 1e-4; yp < 1e5; yp *= 1.5)
        {
            const scalar up = lawOfWall(yp);
            monotone = monotone && up > prev;
            prev = up;
        }
        CHECK(monotone, "u+(y+) is strictly increasing");
        std::cout << "info: u+(1)=" << lawOfWall(1.0) << " u+(11)=" << lawOfWall(11.0)
                  << " u+(100)=" << lawOfWall(100.0) << " u+(1e4)=" << lawOfWall(1e4) << "\n";
    }

    // ---------------- wallState: the inversion ----------------
    {
        const scalar y = 0.001;
        // a wall point in the log layer: u_tau known, u_parallel from the law
        const scalar uTauTrue = 0.5;
        const scalar ypTrue = uTauTrue * y / NU_AIR;
        const scalar uPar = uTauTrue * lawOfWall(ypTrue);
        const WallState s = wallState(uPar, y, NU_AIR);
        CHECK(relNear(s.yPlus, ypTrue, 1e-8), "wallState recovers y+ from u_parallel");
        CHECK(relNear(s.uTau, uTauTrue, 1e-8), "wallState recovers u_tau from u_parallel");
        CHECK(relNear(s.uTau * lawOfWall(s.yPlus), uPar, 1e-9), "u+ = u_parallel/u_tau");
        CHECK(s.yPlus > 30 && s.yPlus < 300, "the test point sits in the wall-function band");

        // the face viscosity is the wall shear, by construction
        const scalar tauOverRho = (NU_AIR + s.nut) * uPar / y;
        CHECK(relNear(tauOverRho, s.uTau * s.uTau, 1e-9),
              "(nu + nut_wall) u / y = u_tau^2 at the wall face");

        // sublayer-resolving point: nut -> 0, u_tau -> sqrt(u nu / y)
        const scalar yNear = 1e-6, uSlow = 1e-3;
        const WallState v = wallState(uSlow, yNear, NU_AIR);
        CHECK(v.yPlus < 0.05, "a fine near-wall cell has y+ << 1");
        CHECK(v.nut < 0.05 * NU_AIR, "nut vanishes in the viscous limit");
        CHECK(relNear(v.uTau, std::sqrt(uSlow * NU_AIR / yNear), 5e-2),
              "u_tau -> sqrt(u nu / y) in the viscous limit");

        const WallState z1 = wallState(0.0, y, NU_AIR);
        const WallState z2 = wallState(uPar, 0.0, NU_AIR);
        const WallState z3 = wallState(uPar, y, 0.0);
        CHECK(z1.uTau == 0 && z1.nut == 0 && z1.yPlus == 0, "no slip speed, no shear");
        CHECK(z2.uTau == 0 && z3.uTau == 0, "zero distance or viscosity gives no shear");
    }

    // ---------------- SST blend functions ----------------
    {
        const scalar nu = NU_AIR, cell = 0.01;
        // near a wall at y = 0.1 mm with the k/omega pair of the log layer
        const scalar nearWall = KOmegaSST::f1(0.01, 10.0, 1e-4, nu, cell);
        // free stream: k = 1e-4 with omega = sqrt(k)/(Cmu^1/4 L), L = 10 mm,
        // measured a third of the domain away from the wall
        const scalar kFree = 1e-4, wFree = 1.83;
        const scalar f1Mid = KOmegaSST::f1(kFree, wFree, 0.15, nu, cell);
        const scalar f1Far = KOmegaSST::f1(kFree, wFree, 0.3, nu, cell);
        const scalar freeStream = f1Far;
        CHECK(nearWall > 0.999, "F1 = 1 in the near-wall k-omega branch");
        CHECK(freeStream < 0.05, "F1 = 0 in the free-stream k-epsilon branch");
        // the blend switches across the outer layer: F1 is intermediate at
        // y = 0.15 and falls from there to the free stream
        CHECK(f1Mid > 0.2 && f1Mid < 0.8, "F1 switches inside the outer layer");
        CHECK(KOmegaSST::f1(kFree, wFree, 0.05, nu, cell) >= f1Mid && f1Mid > f1Far,
              "F1 decreases away from the wall");
        CHECK(KOmegaSST::f2(kFree, wFree, 1e-3, nu) > KOmegaSST::f2(kFree, wFree, 0.3, nu),
              "F2 decreases away from the wall");

        CHECK(near(KOmegaSST::f2(0.01, 10.0, 1e-4, nu), 1.0, 1e-12), "F2 = 1 near the wall");
        const scalar f2Far = KOmegaSST::f2(kFree, wFree, 0.3, nu);
        CHECK(f2Far > 0 && f2Far < 0.2, "F2 falls off in the outer layer");
        CHECK(KOmegaSST::f2(kFree, wFree, 1e4, nu) < 1e-3,
              "F2 -> 0 far from the wall, where nut becomes the k-epsilon form");
        std::cout << "info: F2(free)=" << f2Far << " F1(free)=" << freeStream << "\n";

        bool bounded = true, finite = true;
        for (double yy : {1e-9, 1e-5, 1e-3, 1e-1, 1.0})
        {
            for (double kk : {1e-14, 1e-4, 1.0})
            {
                for (double ww : {1e-14, 1e-2, 1e4})
                {
                    const scalar a = KOmegaSST::f1(kk, ww, yy, nu, cell);
                    const scalar b = KOmegaSST::f2(kk, ww, yy, nu);
                    bounded = bounded && a >= 0 && a <= 1 && b >= 0 && b <= 1;
                    finite = finite && std::isfinite(a) && std::isfinite(b);
                }
            }
        }
        CHECK(bounded && finite, "F1 and F2 stay in [0, 1] for degenerate inputs");
    }

    // ---------------- eddy viscosity and production ----------------
    {
        // a 3-cell column with a plane shear U_x = G y, both turbulence models
        // seeded uniform so only the shear terms are exercised
        constexpr scalar G = 200.0;
        for (const std::string& model : {std::string("kEpsilon"), std::string("kOmegaSST")})
        {
            Box box(true);
            const std::size_t nC = static_cast<std::size_t>(box.mesh.nCells());
            for (std::size_t c = 0; c < nC; ++c)
                box.U.internalField()[c] = Vector{G * box.yOf(static_cast<label>(c)), 0, 0};

            auto m = RansModel::create(box.mesh, box.U, box.phi, NU_AIR, Schemes{}, cfg(model));
            m->initializeUniform(0.02, 50.0);
            const scalar nutSeed = m->nut()[0];
            bool uniformNut = true;
            for (const scalar v : m->nut()) uniformNut = uniformNut && near(v, nutSeed, 1e-12);
            CHECK(uniformNut, model + ": uniform k and dissipation give uniform nut");

            m->correctTurbulence(0.0);

            // production is the lagged nut times |S|^2, and |S| = |G| where the
            // stencil sees two y-neighbours (the cells next to a symmetry plane
            // or a wall read their boundary condition instead)
            const scalar DY = 0.02 / 5.0;
            bool strainOk = true, prodOk = true;
            int interior = 0;
            for (label c = 0; c < box.mesh.nCells(); ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                const scalar yy = box.yOf(c);
                if (yy < DY || yy > 0.02 - DY) continue;
                ++interior;
                strainOk = strainOk && relNear(m->strainRate()[u], std::abs(G), 1e-9);
                prodOk = prodOk && relNear(m->production()[u], nutSeed * G * G, 1e-9);
            }
            CHECK(interior >= 4, "the shear column has interior cells to test");
            CHECK(strainOk, model + ": |S| = |dUx/dy| for a plane shear");
            CHECK(prodOk, model + ": P_k = nut |S|^2");

            // nut from the freshly solved fields, by the model's own formula
            const auto& y = m->cellWallDistance();
            bool nutOk = true;
            for (std::size_t c = 0; c < nC; ++c)
            {
                const scalar k = m->k()[c], d = m->dissipation()[c];
                const scalar expect = model == "kEpsilon"
                    ? KEpsilon::Cmu * k * k / d
                    : KOmegaSST::a1 * k
                      / std::max(KOmegaSST::a1 * d, m->strainRate()[c]
                                                    * KOmegaSST::f2(k, d, y[c], NU_AIR));
                nutOk = nutOk && relNear(m->nut()[c], expect, 1e-12);
            }
            CHECK(nutOk, model + ": nut follows its defining formula");
            CHECK(m->k()[0] > turb::kFloor && m->dissipation()[0] > turb::dissFloor,
                  model + ": production raises k above the floor");
        }
    }

    // ---------------- near-wall boundary values ----------------
    {
        // (0.2 x 0.02 x 0.02), 5 cells in y -> dy = 0.004, first centre 0.002
        constexpr scalar DY = 0.004, Y1 = 0.002;
        for (const std::string& model : {std::string("kEpsilon"), std::string("kOmegaSST")})
        {
            Box box(true);
            const std::size_t nC = static_cast<std::size_t>(box.mesh.nCells());
            for (std::size_t c = 0; c < nC; ++c)
            {
                // a 1/7-power profile: nonzero shear at the wall, so the wall
                // function has something to invert
                const scalar yf = box.yOf(static_cast<label>(c)) / 0.02;
                box.U.internalField()[c] = Vector{13.0 * std::pow(yf, 1.0 / 7.0), 0, 0};
            }

            auto m = RansModel::create(box.mesh, box.U, box.phi, NU_AIR, Schemes{}, cfg(model));
            m->initializeUniform(0.02, 50.0);
            m->correctTurbulence(0.0);

            CHECK(!m->wallFaces().empty(), model + ": the wall patch is found");
            bool shear = true, diss = true, kSlip = true;
            for (std::size_t i = 0; i < m->wallFaces().size(); ++i)
            {
                const label f = m->wallFaces()[i];
                const scalar uTau = m->wallUTau()[i];
                const scalar yp = m->wallYPlus()[i];
                const scalar y = m->wallDistance()[i];
                const label owner = box.mesh.owner()[static_cast<std::size_t>(f)];
                const scalar uPar = box.U.internalField()[static_cast<std::size_t>(owner)].x;

                shear = shear && near(uPar * (NU_AIR + m->wallNut()[i]) / y, uTau * uTau,
                                      1e-9 * uTau * uTau);
                shear = shear && yp > 1.0 && yp < 1e4 && near(y, Y1, 1e-12);

                // k has no Dirichlet at the wall: its value follows from the
                // shear, only the dissipation variable is prescribed
                kSlip = kSlip && m->kField().patchBcType(box.mesh.patchIdByName("bottom"))
                                == "zeroGradient";

                const scalar dWall = m->dissipationField().faceValue(f);
                const scalar expect = (model == "kEpsilon")
                    ? std::pow(uTau, 3) / (turb::kappa * y)
                    : std::sqrt(std::pow(6.0 * NU_AIR / (KOmegaSST::beta1 * y * y), 2.0)
                                + std::pow(uTau / (std::sqrt(KOmegaSST::betaStar) * turb::kappa * y), 2.0));
                diss = diss && relNear(dWall, expect, 1e-12);

                // the momentum equation sees the wall viscosity, not the cell's
                const std::size_t fu = static_cast<std::size_t>(f);
                shear = shear && relNear(m->nuEffFace()[fu], uTau * uTau * y / uPar, 1e-9);
            }
            CHECK(shear, model + ": wall face gives u_tau^2 and the law-of-wall viscosity");
            CHECK(diss, model + ": wall Dirichlet for " + (model == "kEpsilon" ? "epsilon" : "omega"));
            CHECK(kSlip, model + ": k keeps a zero gradient at the wall");
            CHECK(m->maxWallYPlus() >= m->meanWallYPlus(), model + ": y+ report is ordered");
            std::cout << "info: " << model << " y+ " << m->meanWallYPlus() << ".."
                      << m->maxWallYPlus() << " u_tau " << m->meanFrictionVelocity()
                      << " (dy=" << DY << ")\n";
        }
    }

    // ---------------- equilibrium wall cell ----------------
    // On a wall-function mesh the cell next to the wall is closed by the shear,
    // not by its own strain rate: the transported pair is replaced by the state
    // the log layer is in. nut = Cmu k^2/eps of that pair is kappa u_tau y, the
    // log-layer eddy viscosity, which is the number the outer layer has to see.
    // Without this the wall cell's strain rate, taken against the no-slip face
    // value of U, feeds nut back into the production and k grows as k^2 --
    // measured at 1e8 m2/s2 after three steps of the channel case, after which
    // every linear solver in the coupling returns a non-finite residual.
    {
        for (const std::string& model : {std::string("kEpsilon"), std::string("kOmegaSST")})
        {
            const scalar sqrtConst =
                model == "kEpsilon" ? std::sqrt(KEpsilon::Cmu)
                                    : std::sqrt(KOmegaSST::betaStar);
            Box box(true);
            const std::size_t nC = static_cast<std::size_t>(box.mesh.nCells());
            for (std::size_t c = 0; c < nC; ++c)
            {
                const scalar yf = box.yOf(static_cast<label>(c)) / 0.02;
                box.U.internalField()[c] = Vector{13.0 * std::pow(yf, 1.0 / 7.0), 0, 0};
            }

            auto m = RansModel::create(box.mesh, box.U, box.phi, NU_AIR, Schemes{}, cfg(model));
            m->initializeUniform(0.02, 50.0);
            m->correctTurbulence(0.0);

            bool pair = true, nutLog = true, touch = true;
            std::size_t wallCells = 0;
            for (std::size_t i = 0; i < m->wallFaces().size(); ++i)
            {
                const scalar uTau = m->wallUTau()[i], y = m->wallDistance()[i];
                const label owner = m->mesh().owner()[static_cast<std::size_t>(m->wallFaces()[i])];
                const std::size_t c = static_cast<std::size_t>(owner);
                ++wallCells;
                const scalar kEq = uTau * uTau / sqrtConst;
                const scalar dissEq = model == "kEpsilon"
                    ? std::pow(uTau, 3) / (turb::kappa * y)
                    : uTau / (sqrtConst * turb::kappa * y);
                pair = pair && relNear(m->k()[c], kEq, 1e-12)
                    && relNear(m->dissipation()[c], dissEq, 1e-12);
                // k-epsilon's nut is exactly the log-layer value; SST's shear
                // branch a1 k / max(a1 w, |S| F2) can only reduce it
                nutLog = nutLog && (model == "kEpsilon"
                    ? relNear(m->nut()[c], turb::kappa * uTau * y, 1e-9)
                    : m->nut()[c] <= turb::kappa * uTau * y + 1e-30);
                touch = touch && m->nut()[c] > 0;
            }
            CHECK(pair, model + ": wall cell carries the log-layer equilibrium pair");
            CHECK(nutLog, model + ": the equilibrium pair implies nut = kappa u_tau y");
            CHECK(touch, model + ": the wall cell is mixed by the shear, not by the floor");

            // only the cells that touch the wall are closed this way: the outer
            // layer is left to its own transport, which is what carries the
            // shear away from the wall
            std::size_t equilibrated = 0;
            const label firstOwner =
                m->mesh().owner()[static_cast<std::size_t>(m->wallFaces()[0])];
            const scalar kWall = m->k()[static_cast<std::size_t>(firstOwner)];
            for (std::size_t c = 0; c < nC; ++c)
                if (relNear(m->k()[c], kWall, 1e-12)) ++equilibrated;
            CHECK(equilibrated == wallCells,
                  model + ": the closure reaches the wall cells only");
        }
    }

    // ---------------- laminar model is the no-model path ----------------
    {
        Box box(true);
        const std::size_t nC = static_cast<std::size_t>(box.mesh.nCells());
        for (std::size_t c = 0; c < nC; ++c)
            box.U.internalField()[c] = Vector{5.0 * box.yOf(static_cast<label>(c)), 0, 0};

        auto m = RansModel::create(box.mesh, box.U, box.phi, NU_AIR, Schemes{}, cfg("laminar"));
        m->correctTurbulence(0.0);
        bool nuOnly = true;
        for (const scalar v : m->nuEffFace()) nuOnly = nuOnly && near(v, NU_AIR, 1e-30);
        CHECK(nuOnly, "laminar model leaves nuEff = nu at every face");
        bool noNut = true;
        for (const scalar v : m->nut()) noNut = noNut && v == 0;
        CHECK(noNut, "laminar model has no eddy viscosity");
        CHECK(m->dissipationName() == "none", "laminar model has no dissipation variable");
        CHECK(!m->wallFaces().empty() && m->maxWallYPlus() > 0,
              "laminar model still reports the y+ of the mesh");
    }

    // ---------------- cell wall distance ----------------
    {
        // only the patches of mesh type "wall" count: in a half-channel the
        // symmetry plane at the centreline is not a wall
        Box half(true, Vector{0.2, 0.02, 0.02}, 3, 5, 1);
        const auto m = RansModel::create(half.mesh, half.U, half.phi, NU_AIR, Schemes{},
                                                  cfg("laminar"));
        const auto& y = m->cellWallDistance();
        bool ok = true;
        for (label c = 0; c < half.mesh.nCells(); ++c)
            ok = ok && near(y[static_cast<std::size_t>(c)], half.yOf(c), 1e-12);
        CHECK(ok, "wall distance is the distance to the nearest wall patch");

        // a channel walled on both sides reads the nearer of the two
        Box full(true, Vector{0.2, 0.02, 0.02}, 3, 5, 1, true);
        const auto m1 = RansModel::create(full.mesh, full.U, full.phi, NU_AIR, Schemes{},
                                                   cfg("laminar"));
        const auto& y1 = m1->cellWallDistance();
        bool both = true;
        for (label c = 0; c < full.mesh.nCells(); ++c)
        {
            const scalar yy = full.yOf(c);
            both = both && near(y1[static_cast<std::size_t>(c)], std::min(yy, 0.02 - yy), 1e-12);
        }
        CHECK(both, "wall distance picks the nearer of two walls");

        // a domain with no wall has no near-wall region at all
        Box free_(false, Vector{0.1, 0.2, 0.3}, 2, 2, 2);
        const auto m2 = RansModel::create(free_.mesh, free_.U, free_.phi, NU_AIR, Schemes{},
                                                  cfg("laminar"));
        const auto& y2 = m2->cellWallDistance();
        const scalar diag = std::sqrt(0.01 + 0.04 + 0.09);
        bool far = true;
        for (const scalar v : y2) far = far && near(v, diag, 1e-12);
        CHECK(far, "a wall-free domain reads a free-stream distance");
    }

    // ---------------- homogeneous decay: k-epsilon ----------------
    {
        // U = 0 and no walls: the transport equations reduce to their ODEs
        //   dk/dt = -eps,  deps/dt = -C2 eps^2 / k
        // with first integral eps = eps0 (k/k0)^C2 and
        //   k(t) = k0 [1 + (C2 - 1)(eps0/k0) t]^(1/(1 - C2))
        Box box(false, Vector{0.1, 0.1, 0.02}, 2, 2, 1);
        auto m = RansModel::create(box.mesh, box.U, box.phi, NU_AIR, Schemes{}, cfg("kEpsilon"));
        const scalar k0 = 0.05, e0 = 0.025;
        m->initializeUniform(k0, e0);

        const scalar dt = 2e-3;
        for (int step = 0; step < 500; ++step) m->correctTurbulence(dt);
        const scalar t = dt * 500;
        const scalar kExact = k0 * std::pow(1.0 + (KEpsilon::C2 - 1.0) * (e0 / k0) * t,
                                           1.0 / (1.0 - KEpsilon::C2));
        const scalar eExact = e0 * std::pow(kExact / k0, KEpsilon::C2);

        bool uniform = true, decays = true;
        for (std::size_t c = 0; c < m->k().size(); ++c)
        {
            uniform = uniform && near(m->k()[c], m->k()[0], 1e-14 * k0)
                      && near(m->dissipation()[c], m->dissipation()[0], 1e-14 * e0);
            decays = decays && m->k()[c] < k0 && m->k()[c] > kExact * 0.99;
        }
        CHECK(uniform, "k-epsilon: a uniform state stays uniform (no numerical diffusion)");
        CHECK(decays, "k-epsilon: k and eps decay monotonically");
        CHECK(relNear(m->k()[0], kExact, 5e-3), "k-epsilon: k(t) matches the closed form");
        CHECK(relNear(m->dissipation()[0], eExact, 5e-3),
             "k-epsilon: eps(k) matches the first integral");
        std::cout << "info: kEpsilon decay t=" << t << " k=" << m->k()[0] << " (exact "
                  << kExact << ") eps=" << m->dissipation()[0] << " (exact " << eExact
                  << ") nut=" << m->nut()[0] << "\n";
    }

    // ---------------- homogeneous decay: k-omega (SST near-wall branch) ------
    {
        // k0 and omega0 far apart in scale keep F1 = F2 = 1, which is the pure
        // k-omega formulation with beta1:
        //   dk/dt = -beta* k omega,  dw/dt = -beta1 w^2
        //   w(t) = w0/(1 + beta1 w0 t),  k(t) = k0 (1 + beta1 w0 t)^(-beta*/beta1)
        Box box(false, Vector{0.1, 0.1, 0.02}, 2, 2, 1);
        KOmegaSST model(box.mesh, box.U, box.phi, NU_AIR, Schemes{}, cfg("kOmegaSST"));
        auto& m = model;
        const scalar k0 = 1.0, w0 = 1e-3;
        m.initializeUniform(k0, w0);
        CHECK(KOmegaSST::f1(k0, w0, m.cellWallDistance()[0], NU_AIR, model.cellSize()) > 0.999,
              "the decay state sits in the k-omega branch");

        const scalar dt = 1.0;
        for (int step = 0; step < 2000; ++step) m.correctTurbulence(dt);
        const scalar growth = 1.0 + KOmegaSST::beta1 * w0 * dt * 2000;
        const scalar wExact = w0 / growth;
        const scalar kExact = k0 * std::pow(growth, -KOmegaSST::betaStar / KOmegaSST::beta1);

        CHECK(relNear(m.dissipation()[0], wExact, 5e-3), "SST: omega(t) matches k-omega decay");
        CHECK(relNear(m.k()[0], kExact, 5e-3), "SST: k(t) matches beta*/beta1 decay");
        // at zero shear the nut formula reduces to k/omega
        CHECK(relNear(m.nut()[0], m.k()[0] / m.dissipation()[0], 1e-9),
              "SST: nut = k/omega where the shear branch is inactive");
        std::cout << "info: SST decay k=" << m.k()[0] << " (exact " << kExact
                  << ") omega=" << m.dissipation()[0] << " (exact " << wExact << ")\n";
    }

    // ---------------- case configuration ----------------
    {
        const Dictionary d = Dictionary::parseString("nu 1e-5;");
        const TurbulenceConfig c = TurbulenceConfig::fromCase(d);
        CHECK(c.model == "laminar", "no turbulence block means laminar");

        const Dictionary t = Dictionary::parseString(
            "turbulence { model kOmegaSST; intensity 0.03; lengthScale 0.01; relaxation 0.9; }");
        const TurbulenceConfig c2 = TurbulenceConfig::fromCase(t);
        CHECK(c2.model == "kOmegaSST" && near(c2.intensity, 0.03) && near(c2.lengthScale, 0.01)
              && near(c2.relaxation, 0.9), "turbulence block is read");

        Box box(true);
        bool badModel = false, badIntensity = false, badRelax = false;
        try
        {
            RansModel::create(box.mesh, box.U, box.phi, NU_AIR, Schemes{}, cfg("kLlama"));
        }
        catch (const Error&) { badModel = true; }
        try
        {
            TurbulenceConfig c = cfg("kEpsilon");
            c.intensity = -1;
            RansModel::create(box.mesh, box.U, box.phi, NU_AIR, Schemes{}, c);
        }
        catch (const Error&) { badIntensity = true; }
        try
        {
            RansModel::create(box.mesh, box.U, box.phi, NU_AIR, Schemes{}, cfg("kEpsilon", 1.5));
        }
        catch (const Error&) { badRelax = true; }
        CHECK(badModel, "an unknown model name is refused");
        CHECK(badIntensity, "a negative intensity is refused");
        CHECK(badRelax, "a relaxation factor outside (0, 1] is refused");

        // the inlet intensity and length scale seed k and the dissipation
        Box inlet(true);
        inlet.U.setPatchValue(inlet.mesh.patchIdByName("inlet"), "fixedValue", Vector{10, 0, 0});
        TurbulenceConfig c3 = cfg("kEpsilon");
        c3.intensity = 0.05;
        c3.lengthScale = 0.02;
        const auto m = RansModel::create(inlet.mesh, inlet.U, inlet.phi, NU_AIR, Schemes{}, c3);
        const scalar kExpect = 1.5 * std::pow(10.0 * 0.05, 2);
        const scalar eExpect = std::pow(KEpsilon::Cmu, 0.25) * std::pow(kExpect, 1.5) / 0.02;
        CHECK(relNear(m->k()[0], kExpect, 1e-9), "k0 = 1.5 (U intensity)^2 at the inlet");
        CHECK(relNear(m->dissipation()[0], eExpect, 1e-9), "eps0 = Cmu^0.25 k0^1.5 / L");
        CHECK(m->kField().patchBcType(inlet.mesh.patchIdByName("inlet")) == "fixedValue",
              "inlet patches carry the prescribed turbulence");
    }

    if (failures)
    {
        std::cout << failures << " test(s) FAILED\n";
        return 1;
    }
    std::cout << "all turbulence unit tests passed\n";
    return 0;
}
