// Unit tests for the VOF layer: what the phase-fraction transport owes the
// scheme (conservation, the discrete maximum principle) and what the interface
// geometry operators return on fields whose answer is known without a phase
// fraction at all -- a hand-built normal field.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "cf/fvm/BlockMesh.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/multiphase/VoF.h"

using namespace cf;

static int failures = 0;

#define CHECK(cond, msg)                                                          \
    do {                                                                          \
        if (!(cond)) { ++failures; std::cout << "FAIL: " << (msg) << "\n"; }      \
        else          { std::cout << "ok:   " << (msg) << "\n"; }                 \
    } while (0)

static bool relNear(double a, double b, double tol)
{
    return std::abs(a - b) < tol * std::max(1.0, std::max(std::abs(a), std::abs(b)));
}

namespace
{

/// A unit cube column with named sides, so the tests can say which patch feeds
/// the flux. `walls` puts an inlet/outlet pair across x; without it the box is
/// closed and nothing can enter.
struct Box
{
    polyMesh mesh;

    explicit Box(bool walls, label nx, label ny, label nz,
                 const Vector& size = Vector{1.0, 1.0, 0.1})
        : mesh(make(walls, nx, ny, nz, size))
    {}

    static polyMesh make(bool walls, label nx, label ny, label nz, const Vector& size)
    {
        BlockSides s;
        if (walls)
        {
            s.xmin = {"inlet", "inlet"};
            s.xmax = {"outlet", "outlet"};
        }
        else
        {
            s.xmin = {"sx0", "symmetry"};
            s.xmax = {"sx1", "symmetry"};
        }
        s.ymin = {"sy0", "symmetry"};
        s.ymax = {"sy1", "symmetry"};
        s.zmin = {"sz0", "symmetry"};
        s.zmax = {"sz1", "symmetry"};
        return BlockMesh::create(Vector::zero(), size, nx, ny, nz, s);
    }

    label patch(const std::string& name) const { return mesh.patchIdByName(name); }
};

/// The face flux of a constant velocity: phi = U & Sf, which is exactly
/// divergence-free on any mesh because sum_f U & Sf = U & sum_f Sf = 0 per cell.
Field<scalar> uniformPhi(const polyMesh& mesh, const Vector& U)
{
    Field<scalar> phi(static_cast<std::size_t>(mesh.nFaces()), 0.0);
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        phi[u] = dot(U, mesh.faceSf()[u]);
    }
    return phi;
}

} // namespace

int main()
{
    // ---- 1. one transport step, against the discrete solution it owes ----
    {
        const label nx = 8;
        Box box(true, nx, 1, 1, Vector{1.0, 0.1, 0.1});
        volScalarField alpha(box.mesh, "alpha", 0.0);
        alpha.setPatchValue(box.patch("inlet"), "fixedValue", 1.0);
        Field<scalar>& a = alpha.internalField();
        for (label i = 0; i < nx / 2; ++i)
            a[static_cast<std::size_t>(i)] = 1.0;   // phase 1 behind x = 0.5
        const Vector U{1.0, 0.0, 0.0};
        const Field<scalar> phi = uniformPhi(box.mesh, U);
        const scalar dx = 1.0 / static_cast<scalar>(nx);
        const scalar dt = 0.5 * dx;                 // the Courant number is 0.5
        const scalar area = 0.1 * 0.1;

        Schemes sch;
        sch.ddt = "Euler";
        sch.div = "upwind";
        vof::Controls vc;
        vc.tol = 1e-12;
        const scalar m0 = vof::phaseVolume(box.mesh, a);
        const SolverReport rep = vof::advectAlpha(alpha, phi, dt, sch, vc);

        // Implicit Euler with upwind convection on a step field: behind the front
        // nothing changes, and k cells ahead of it the update is
        // a = (theta / (1 + theta))^k -- the geometric tail a diffused front
        // carries into the dry region, which is why a case that follows a 0.5
        // crossing measures it rather than assuming it away.
        bool profile = true, bounded = true;
        const scalar theta = U.x * dt / dx;
        const scalar ratio = theta / (1.0 + theta);
        for (label i = 0; i < nx; ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            const scalar expect = i < nx / 2
                ? 1.0 : std::pow(ratio, static_cast<scalar>(i - nx / 2 + 1));
            profile = profile && relNear(a[u], expect, 1e-9);
            bounded = bounded && a[u] >= 0.0 && a[u] <= 1.0;
        }
        CHECK(rep.converged, "the transport solve converges");
        CHECK(bounded, "alpha stays in [0,1] without clipping");
        CHECK(profile, "the one-step field is the discrete solution of backward "
                       "Euler plus upwind");
        // and the phase volume moved by exactly what the boundary exchanged: the
        // inlet's 1 in, the outlet's own value out
        const scalar m1 = vof::phaseVolume(box.mesh, a);
        CHECK(relNear(m1 - m0, U.x * area * dt * (1.0 - a[static_cast<std::size_t>(nx - 1)]),
                      1e-10),
              "the phase volume changes by exactly what the boundary exchanged");
    }

    // ---- 2. a closed box cannot change the phase volume on its own ----
    {
        Box box(false, 3, 3, 1);
        volScalarField alpha(box.mesh, "alpha", 0.0);
        // no flux anywhere: a uniform velocity through a symmetry plane is not
        // admissible, so phi is zero and nothing can move
        const Field<scalar> phi(static_cast<std::size_t>(box.mesh.nFaces()), 0.0);
        alpha.internalField()[4] = 0.75;
        Schemes sch;
        sch.ddt = "Euler";
        sch.div = "upwind";
        const scalar m0 = vof::phaseVolume(box.mesh, alpha.internalField());
        vof::advectAlpha(alpha, phi, 0.1, sch, vof::Controls{});
        CHECK(relNear(vof::phaseVolume(box.mesh, alpha.internalField()), m0, 1e-14),
              "a zero flux leaves the phase field, and its volume, alone");
    }

    // ---- 3. curvature of a hand-built radial normal: kappa = +1/r ----
    {
        const label n = 21;
        Box box(false, n, n, 1, Vector{1.0, 1.0, 0.1});
        const Vector c{0.5, 0.5, 0.0};
        Field<Vector> normal(static_cast<std::size_t>(box.mesh.nCells()), Vector::zero());
        for (label i = 0; i < box.mesh.nCells(); ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            const Vector r = box.mesh.cellCentres()[u] - c;
            const scalar m = mag(r);
            if (m > 0) normal[u] = r * (1.0 / m);          // pointing away from c
        }
        const Field<scalar> k = vof::interfaceCurvature(box.mesh, normal);
        // The stencil reads cell-centre normals averaged onto the face, so its
        // error is the second derivative of r_hat, which is 1/r^2, measured
        // against the divergence itself: |kappa*r - 1| = O((dx/r)^2). Two limits
        // follow, and both are truncation rather than the operator. Near the
        // axis the circle is not resolved and the relative error grows. At the
        // box edge the faces touching the boundary drop out of the stencil,
        // because only faces whose two cells both carry an interface count.
        scalar worst = 0, worstR = 0;
        const scalar dx = 1.0 / static_cast<scalar>(n);
        for (label i = 0; i < box.mesh.nCells(); ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            const Vector x = box.mesh.cellCentres()[u];
            const scalar r = mag(x - c);
            const scalar toBoundary = std::min(std::min(x.x, 1.0 - x.x), std::min(x.y, 1.0 - x.y));
            if (r < 1e-12 || toBoundary < 2.5 * dx) continue;
            const scalar e = std::abs(k[u] * r - 1.0) * r * r / (dx * dx);
            if (e > worst) { worst = e; worstR = r; }
        }
        std::cout << "     radial normal: |kappa*r - 1| / (dx/r)^2 < " << worst
                  << " (worst at r=" << worstR << ")\n";
        CHECK(worst < 3.0,
              "div of a hand-built radial unit normal is 1/r to second order in dx/r");
    }

    // ---- 4. a uniform normal has no curvature, wherever it stops ----
    {
        const label n = 9;
        Box box(false, n, n, 1, Vector{1.0, 1.0, 0.1});
        Field<Vector> normal(static_cast<std::size_t>(box.mesh.nCells()), Vector::zero());
        // a two-cell-thick band across the box, normal to its own extent: the
        // interface of a plane, and the stencil a band cell has outside the band
        // has nothing to contribute
        const auto idx = [](label i, label j)
        { return static_cast<std::size_t>(j * n + i); };
        for (label j = 0; j < n; ++j)
            for (label i = 3; i <= 4; ++i)
                normal[idx(i, j)] = Vector{-1.0, 0.0, 0.0};
        const Field<scalar> k = vof::interfaceCurvature(box.mesh, normal);
        scalar inside = 0, outside = 0;
        for (label j = 0; j < n; ++j)
            for (label i = 0; i < n; ++i)
            {
                const scalar v = std::abs(k[idx(i, j)]);
                if (i == 3 || i == 4) inside = std::max(inside, v);
                else                  outside = std::max(outside, v);
            }
        // the regression this test exists for: averaging the band's normal with
        // the zero of the neighbour outside it used to put +-1/(2 dx) on both
        // sides of the band, which is what a plane interface must never read
        CHECK(inside == 0.0 && outside == 0.0,
              "a plane band reads exactly zero curvature, band-edge faces included");
        std::cout << "     plane band: max |kappa| = " << inside << " in the band, "
                  << outside << " outside\n";
    }

    // ---- 5. the normal of a plane alpha is the plane's normal ----
    {
        const label n = 8;
        Box box(false, n, n, 1, Vector{1.0, 1.0, 0.1});
        volScalarField alpha(box.mesh, "alpha", 0.0);
        const scalar dx = 1.0 / static_cast<scalar>(n);
        for (label i = 0; i < box.mesh.nCells(); ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            const scalar r = alpha.mesh().cellCentres()[u].x - 0.5;
            alpha.internalField()[u] = 0.5 * (1.0 - std::tanh(r / dx));
        }
        const Field<Vector> nrm = vof::interfaceNormal(alpha, "linear", 0.1);
        scalar nBand = 0, worstDir = 0;
        for (std::size_t u = 0; u < nrm.size(); ++u)
        {
            if (mag(nrm[u]) == 0) continue;
            ++nBand;
            worstDir = std::max(worstDir, std::abs(nrm[u].y) + std::abs(nrm[u].z));
            worstDir = std::max(worstDir, std::abs(nrm[u].x + 1.0));
        }
        std::cout << "     plane alpha: " << nBand << " of " << nrm.size()
                  << " cells carry the band\n";
        CHECK(nBand > 0 && nBand < static_cast<scalar>(nrm.size()),
              "the band is the interface, not the whole field");
        CHECK(worstDir < 1e-12, "grad of a one-dimensional alpha is normal to it");
        CHECK(relNear(vof::phaseVolume(box.mesh, alpha.internalField()), 0.05, 1e-12),
              "phase volume sums alpha over the mesh");
        const Vector cen = vof::phaseCentroid(box.mesh, alpha.internalField());
        CHECK(std::abs(cen.x - 0.25) < 0.03 && std::abs(cen.y - 0.5) < 1e-12,
              "the phase centroid sits in the half the phase fills");
    }

    if (failures)
    {
        std::cout << failures << " VOF unit test(s) FAILED\n";
        return 1;
    }
    std::cout << "all VOF unit tests passed\n";
    return 0;
}
