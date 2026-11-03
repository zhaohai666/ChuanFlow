#include "cf/multiphase/VoF.h"

#include <algorithm>
#include <cmath>

#include "cf/common/Error.h"

namespace cf::vof
{

SolverReport advectAlpha(volScalarField& alpha, const Field<scalar>& phi,
                         scalar dt, const Schemes& schemes, const Controls& ctl)
{
    if (!(dt > 0)) cfFatal("vof: dt must be positive");
    if (static_cast<label>(phi.size()) != alpha.mesh().nFaces())
        cfFatal("vof: the flux field must hold one value per face");

    const polyMesh& mesh = alpha.mesh();
    fvMatrix eqn(mesh, viewOf(alpha), "alphaEqn");
    eqn.schemes() = schemes;
    // the relaxation is a solver knob for iterative coupling; a transport step
    // has to move the whole solution
    eqn.schemes().relaxation = 1.0;
    fvm::ddt(eqn, alpha.internalField(), dt);
    fvm::div(eqn, phi, schemes.div);

    const SolverReport rep = eqn.solve(fvMatrix::Controls{ctl.tol, ctl.maxIter, 0});
    alpha.internalField() = eqn.solution()[0];
    return rep;
}

scalar phaseVolume(const polyMesh& mesh, const Field<scalar>& alpha)
{
    const Field<scalar>& v = mesh.cellVolumes();
    scalar s = 0;
    for (label c = 0; c < mesh.nCells(); ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        s += alpha[u] * v[u];
    }
    return s;
}

Vector phaseCentroid(const polyMesh& mesh, const Field<scalar>& alpha)
{
    const Field<scalar>& v = mesh.cellVolumes();
    const Field<Vector>& c = mesh.cellCentres();
    Vector x = Vector::zero();
    scalar mass = 0;
    for (label i = 0; i < mesh.nCells(); ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        x += alpha[u] * v[u] * c[u];
        mass += alpha[u] * v[u];
    }
    return mass > 0 ? x * (1.0 / mass) : Vector::zero();
}

Field<Vector> interfaceNormal(const volScalarField& alpha, const std::string& scheme,
                              scalar thresholdFraction)
{
    const Field<Vector> g = fvc::grad(alpha, scheme);
    scalar peak = 0;
    for (const Vector& v : g) peak = std::max(peak, mag(v));
    const scalar cut = thresholdFraction * peak;

    Field<Vector> n(static_cast<std::size_t>(alpha.mesh().nCells()), Vector::zero());
    for (std::size_t i = 0; i < static_cast<std::size_t>(g.size()); ++i)
    {
        const scalar m = mag(g[i]);
        if (m > cut) n[i] = g[i] * (1.0 / m);
    }
    return n;
}

Field<scalar> interfaceCurvature(const polyMesh& mesh, const Field<Vector>& normal)
{
    const label nC = mesh.nCells();
    const label nInt = mesh.nInternalFaces();
    const Field<Vector>& sf = mesh.faceSf();
    const Field<label>& own = mesh.owner();
    const Field<label>& nbr = mesh.neighbour();
    const Field<scalar> ow = makeInternalOwnerWeights(mesh);
    const Field<scalar>& v = mesh.cellVolumes();

    Field<scalar> k(static_cast<std::size_t>(nC), 0.0);
    for (label f = 0; f < nInt; ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const std::size_t o = static_cast<std::size_t>(own[u]);
        const std::size_t n = static_cast<std::size_t>(nbr[u]);
        // a zero normal means the cell carries no interface, so this face is
        // outside the band and drops out of both stencils
        if (mag2(normal[o]) == 0 || mag2(normal[n]) == 0) continue;

        const scalar w = ow[u];
        const Vector nf = w * normal[o] + (1.0 - w) * normal[n];
        const scalar s = dot(nf, sf[u]);
        k[o] += (s - dot(normal[o], sf[u])) / v[o];
        k[n] -= (s - dot(normal[n], sf[u])) / v[n];
    }
    return k;
}

Field<Vector> surfaceTensionForce(scalar sigmaOverRho,
                                  const volScalarField& alpha,
                                  const Field<scalar>& kappa,
                                  const std::string& gradScheme)
{
    const polyMesh& mesh = alpha.mesh();
    if (static_cast<label>(kappa.size()) != mesh.nCells())
        cfFatal("vof: the curvature field must hold one value per cell");

    const Field<Vector> g = fvc::grad(alpha, gradScheme);
    const label nC = mesh.nCells();
    Field<Vector> f(static_cast<std::size_t>(nC), Vector::zero());
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        f[u] = g[u] * (-sigmaOverRho * kappa[u]);
    }
    return f;
}

} // namespace cf::vof
