// Unit tests for the transient incompressible PISO coupling: the inviscid
// response to a body force, the steady limit against the discrete parabola,
// and the divergence-free property of the corrected flux.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "cf/common/Error.h"
#include "cf/fvm/BlockMesh.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/PisoSolver.h"
#include "cf/fvm/PressureCorrection.h"

using namespace cf;

static int failures = 0;

#define CHECK(cond, msg)                                                      \
    do {                                                                      \
        if (!(cond)) { ++failures; std::cout << "FAIL: " << (msg) << "\n"; }   \
        else          { std::cout << "ok:   " << (msg) << "\n"; }              \
    } while (0)

static bool near(double a, double b, double tol = 1e-10)
{
    return std::abs(a - b) < tol;
}

namespace
{

struct Channel
{
    polyMesh mesh;
    volVectorField U;
    volScalarField p;
    surfaceScalarField phi;
    scalar nu;
    Vector force;

    Channel(const Vector& size, const Vector& cells, scalar nu, Vector force)
        : mesh(makeMesh(size, cells)), U(mesh, "U"), p(mesh, "p"), phi(mesh, "phi", 0.0),
          nu(nu), force(force)
    {
        U.setPatchValue(mesh.patchIdByName("bottom"), "fixedValue", Vector::zero());
        U.setPatchValue(mesh.patchIdByName("top"), "fixedValue", Vector::zero());
        U.setPatchValue(mesh.patchIdByName("left"), "zeroGradient", Vector::zero());
        U.setPatchValue(mesh.patchIdByName("right"), "zeroGradient", Vector::zero());
        p.setPatchValue(mesh.patchIdByName("left"), "fixedValue", 0.0); // gauge
        for (const char* n : {"right", "bottom", "top"})
            p.setPatchValue(mesh.patchIdByName(n), "zeroGradient", 0.0);
        phi.correctBoundaryFlux(U);
    }

    static polyMesh makeMesh(const Vector& size, const Vector& cells)
    {
        BlockSides s;
        s.xmin = {"left", "symmetry"};
        s.xmax = {"right", "symmetry"};
        s.ymin = {"bottom", "wall"};
        s.ymax = {"top", "wall"};
        return BlockMesh::create(Vector::zero(), size,
                                 static_cast<label>(cells.x), static_cast<label>(cells.y),
                                 static_cast<label>(cells.z), s);
    }
};

Schemes schemes()
{
    Schemes sc;
    sc.ddt = "Euler";
    return sc;
}

} // namespace

int main()
{
    constexpr scalar H = 0.1, DY = 0.01, NU = 0.01, A = 0.1;

    // ---- a non-positive time step is refused, not silently skipped ----
    {
        Channel c(Vector{0.2, H, 0.05}, Vector{2, 10, 1}, NU, Vector{A, 0, 0});
        PisoControls ctl;
        ctl.dt = 0.0;
        bool threw = false;
        try
        {
            solveTransientIncompressiblePiso(c.mesh, c.U, c.p, c.phi, c.nu, schemes(), ctl);
        }
        catch (const Error&) { threw = true; }
        CHECK(threw, "zero time step throws");
    }

    // ---- inviscid first step: every cell accelerates by f*dt ----
    {
        Channel c(Vector{0.2, H, 0.05}, Vector{2, 10, 1}, 1e-8, Vector{A, 0, 0});
        PisoControls ctl;
        ctl.dt = 0.002;
        ctl.correctors = 2;
        ctl.bodyForce = c.force;
        const PisoStepResult r =
            solveTransientIncompressiblePiso(c.mesh, c.U, c.p, c.phi, c.nu, schemes(), ctl);
        scalar uMin = 1e300, uMax = -1e300;
        for (const Vector& u : c.U.internalField())
        {
            uMin = std::min(uMin, u.x);
            uMax = std::max(uMax, u.x);
        }
        scalar pMax = 0;
        for (const scalar v : c.p.internalField()) pMax = std::max(pMax, std::abs(v));
        CHECK(near(uMin, A * ctl.dt, 1e-9) && near(uMax, A * ctl.dt, 1e-9),
              "inviscid step gives u = f*dt everywhere");
        CHECK(pMax < 1e-9, "uniform acceleration leaves p at the gauge");
        CHECK(r.continuityError < 1e-14 && r.correctorsUsed == 1,
              "predicted flux is already divergence-free");
    }

    // ---- steady limit: the discrete parabola of the wall-cell scheme ----
    {
        Channel c(Vector{0.2, H, 0.05}, Vector{2, 10, 1}, NU, Vector{A, 0, 0});
        PisoControls ctl;
        ctl.dt = 0.002;   // nu*dt/DY^2 = 0.2, within the documented limit
        ctl.correctors = 2;
        ctl.continuityTol = 1e-12;
        ctl.bodyForce = c.force;
        scalar worstCont = 0, prevU = -1;
        bool monotone = true;
        PisoStepResult r;
        for (int step = 0; step < 1500; ++step) // t = 3 diffusion times
        {
            r = solveTransientIncompressiblePiso(c.mesh, c.U, c.p, c.phi, c.nu, schemes(), ctl);
            worstCont = std::max(worstCont, r.continuityError);
            scalar u = 0;
            for (const Vector& uu : c.U.internalField()) u = std::max(u, uu.x);
            monotone = monotone && u >= prevU - 1e-14;
            prevU = u;
        }
        CHECK(r.converged, "corrector loop meets continuityTol");
        CHECK(worstCont < 1e-12, "corrected flux stays divergence-free for 1500 steps");
        CHECK(monotone, "start-up rises monotonically");

        // the cell-centred scheme imposes u = 0 half a cell below the first
        // node, so the discrete steady state is the analytic parabola raised
        // by (a/2nu)(dy/2)^2
        const auto& cc = c.mesh.cellCentres();
        scalar profErr = 0, curvErr = 0, vMax = 0;
        std::vector<std::size_t> order;
        for (label j = 0; j < c.mesh.nCells(); ++j)
            if (std::abs(cc[static_cast<std::size_t>(j)].x - cc[0].x) < 1e-9)
                order.push_back(static_cast<std::size_t>(j));
        std::sort(order.begin(), order.end(),
                  [&](std::size_t a, std::size_t b) { return cc[a].y < cc[b].y; });
        for (std::size_t k = 0; k < order.size(); ++k)
        {
            const std::size_t j = order[k];
            const scalar y = cc[j].y, u = c.U.internalField()[j].x;
            vMax = std::max(vMax, std::abs(c.U.internalField()[j].y));
            profErr = std::max(profErr,
                               std::abs(u - (A / (2 * NU)) * (y * (H - y) + (DY / 2) * (DY / 2))));
            if (k > 0 && k + 1 < order.size())
            {
                const scalar d2 = c.U.internalField()[order[k - 1]].x - 2 * u
                                + c.U.internalField()[order[k + 1]].x;
                curvErr = std::max(curvErr, std::abs(d2 + (A / NU) * DY * DY));
            }
        }
        const scalar uPeak = (A / (2 * NU)) * (H / 2) * (H / 2);
        CHECK(profErr < 1e-5 * uPeak, "steady profile matches the discrete parabola");
        CHECK(curvErr < 1e-9, "steady profile satisfies nu laplacian(u) = -a at every node");
        CHECK(vMax < 1e-5 * uPeak, "cross-flow stays zero");
        std::cout << "info: profErr " << profErr << " curvErr " << curvErr
                  << " vMax " << vMax << " worstCont " << worstCont << "\n";
    }

    if (failures)
    {
        std::cout << failures << " test(s) FAILED\n";
        return 1;
    }
    std::cout << "all piso unit tests passed\n";
    return 0;
}
