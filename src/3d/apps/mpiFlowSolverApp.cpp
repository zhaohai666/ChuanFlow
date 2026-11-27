// chflow-3d-mpi-flow: the full incompressible SIMPLE pressure-velocity
// coupling on a real MPI transport over a Decomposition-split blockMesh
// (5.233). One domain per rank.
//
// This is NOT the §5.136 petscDistributed path. There, one global Krylov
// solve shares the whole matrix across ranks. Here the flow solver keeps the
// serial SIMPLE structure of src/3d/src/solvers/simpleSolver.cpp -- momentum,
// HbyA blend, predictor flux, full-p Poisson, flux/velocity correction -- but
// every operator runs as its per-domain twin on rank-local slots: the grad /
// div / laplacian stencils are the same formulas (Gauss sums with the
// ownerWeight and faceSkew the DomainMesh carries, harmonic-free uniform nu,
// upwind convection on the face flux), and every stencil application is
// preceded by an MPI_Sendrecv halo exchange of the field it reads
// (mpiTransport.h). Domains never share an fvMatrix across an interface;
// between halo syncs each rank sweeps Gauss-Seidel on its owned cells only.
// Global quantities -- the pressure-Poisson residual norm, the continuity
// error, the adjustPhi volume-flow sums, the convergence stops -- are
// MPI_Allreduce collectives.
//
// Boundary honesty: the coordinate-median cut guarantees interface cells are
// complete cells (no hanging nodes), so a two-point stencil across a cut is
// well defined. The SUM reductions match the serial left-to-right accumulate
// to round-off only (MPI tree order); MIN/MAX are order-independent and
// exact -- this is the source of the 1e-8 relative comparison tolerance the
// checker uses between rank counts: the converged fields solve the same
// discrete system, but the GS sweep order and the reduction order differ per
// decomposition, so the last digits move.
//
// cf3d / cf_common stay MPI-free; MPI lives only here and in mpiTransport.h.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Field.h"
#include "cf/common/Log.h"
#include "cf/fvm/MeshLoader.h"
#include "cf/fvm/PolyMesh.h"
#include "cf/parallel/DomainDecomposition.h"

#include "mpiTransport.h"

using namespace cf;
using cf::parallel::DomainFace;
using cf::parallel::DomainMesh;

namespace
{

std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

/// Local slots of the global owner and neighbour behind a stencil face --
/// the same helper domainOperators.cpp uses serially.
inline void ownerNeighbour(const DomainFace& df, label& o, label& n)
{
    if (df.sign > 0) { o = df.cell; n = df.peer; }
    else             { o = df.peer; n = df.cell; }
}

struct PatchBC
{
    std::string uType = "zeroGradient"; // fixedValue | zeroGradient
    Vector uValue = Vector::zero();
    std::string pType = "zeroGradient"; // fixedValue | zeroGradient
    scalar pValue = 0;
    bool zeroNormal = false;            // symmetry / empty: phi stays 0
};

/// Face-transmissibility dot(Sf, d)/|d|^2 on the global mesh (the serial
/// orthoCoeff of laplacianScheme.cpp), read through the global face id each
/// DomainFace carries. Precomputing it per rank from the loaded mesh mirrors
/// what mpiLinearSolverApp does for the global CSR: the setup is global, the
/// sweep is rank-local.
struct FaceGeom
{
    Field<scalar> ortho;      // dot(Sf, d) / |d|^2
    Field<label> patch;       // boundary patch id, -1 internal
    Field<scalar> magSf;
};

FaceGeom makeFaceGeom(const polyMesh& mesh)
{
    const label nF = mesh.nFaces();
    const label nInt = mesh.nInternalFaces();
    FaceGeom g;
    g.ortho.assign(static_cast<std::size_t>(nF), 0.0);
    g.patch = mesh.facePatchMap();
    g.magSf = mesh.faceAreas();
    const auto& sf = mesh.faceSf();
    const auto& dcv = mesh.deltaCoeffVec();
    for (label f = 0; f < nF; ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const scalar d2 = dot(dcv[u], dcv[u]);
        if (d2 <= 0) cfFatal("mpiFlow: zero face distance at face " + std::to_string(f));
        g.ortho[u] = dot(sf[u], dcv[u]) / d2;
        (void)nInt;
    }
    return g;
}

/// The boundary face value of a velocity patch condition: the prescribed
/// vector, or the owner cell value for gradient types (zeroGradient and
/// symmetry both keep the normal component; phi is separately forced to 0 on
/// zero-normal patches, mirroring surfaceScalarField::correctBoundaryFlux).
inline Vector uFaceValue(const PatchBC& bc, const Vector& cellU)
{
    if (bc.uType == "fixedValue") return bc.uValue;
    return cellU;
}

inline scalar pFaceValue(const PatchBC& bc, scalar cellP)
{
    if (bc.pType == "fixedValue") return bc.pValue;
    return cellP;
}

/// Gauss gradient of a haloed local field, owned slots written and the halo
/// slots re-synchronised. "linear" is the plain bounded-interpolation sum;
/// "corrected" adds the grad_f . skew_f term, exchanging the gradient between
/// passes exactly like assembledGrad does serially.
void localGrad(const parallel::Decomposition& dec, label rank, MPI_Comm comm,
               const DomainMesh& dm, const Field<scalar>& pLocal,
               const std::vector<PatchBC>& bc, const FaceGeom& fg,
               const std::string& scheme, Field<Vector>& grad)
{
    const std::size_t nL = static_cast<std::size_t>(dm.nLocal());
    std::fill(grad.begin(), grad.end(), Vector::zero());
    Field<Vector> sum(nL, Vector::zero());

    const auto faceValue = [&](std::size_t k) -> scalar
    {
        const DomainFace& df = dm.faces()[k];
        if (df.peer < 0)
            return pFaceValue(bc[static_cast<std::size_t>(fg.patch[df.face])],
                              pLocal[static_cast<std::size_t>(df.cell)]);
        label o = -1, n = -1;
        ownerNeighbour(df, o, n);
        const scalar w = dm.faceOwnerWeight()[k];
        return w * pLocal[static_cast<std::size_t>(o)]
             + (1.0 - w) * pLocal[static_cast<std::size_t>(n)];
    };

    for (std::size_t k = 0; k < dm.faces().size(); ++k)
    {
        const DomainFace& df = dm.faces()[k];
        const Vector flux = dm.faceSf()[k] * faceValue(k);
        if (df.sign > 0) sum[static_cast<std::size_t>(df.cell)] += flux;
        else             sum[static_cast<std::size_t>(df.cell)] -= flux;
    }
    for (label i = 0; i < dm.nOwned(); ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        grad[u] = sum[u] / dm.cellVolumes()[u];
    }

    if (scheme != "corrected") return;

    // skew correction; faceSkew is zero on an orthogonal blockMesh, so the
    // correction vanishes there exactly as in the serial scheme. The gradient
    // is exchanged between passes: the same collective assembledGrad issues.
    for (int pass = 0; pass < 2; ++pass)
    {
        parallel::mpi::exchange(dec, rank, grad, comm);
        std::fill(sum.begin(), sum.end(), Vector::zero());
        for (std::size_t k = 0; k < dm.faces().size(); ++k)
        {
            const DomainFace& df = dm.faces()[k];
            scalar tf = faceValue(k);
            if (df.peer >= 0)
            {
                label o = -1, n = -1;
                ownerNeighbour(df, o, n);
                const scalar w = dm.faceOwnerWeight()[k];
                const Vector gF = w * grad[static_cast<std::size_t>(o)]
                                + (1.0 - w) * grad[static_cast<std::size_t>(n)];
                tf += dot(gF, dm.faceSkew()[k]);
            }
            const Vector flux = dm.faceSf()[k] * tf;
            if (df.sign > 0) sum[static_cast<std::size_t>(df.cell)] += flux;
            else             sum[static_cast<std::size_t>(df.cell)] -= flux;
        }
        for (label i = 0; i < dm.nOwned(); ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            grad[u] = sum[u] / dm.cellVolumes()[u];
        }
    }
}

/// One momentum solve: upwind convection on the current phi, uniform-nu
/// diffusion (two-point; the skew correction vanishes on an orthogonal box),
/// explicit -V grad(p) + V bodyForce on the right. Gauss-Seidel on owned
/// cells with the halo frozen -- the interface coupling a shared matrix would
/// carry is instead re-synchronised by the next halo exchange, which is the
/// deliberate difference to the §5.136 petscDistributed global-Krylov path.
/// Writes the diagonal (rAU = V/aDiag) into aDiag[owned].
void momentumSolve(const DomainMesh& dm, Field<Vector>& U,
                   const Field<scalar>& phi, const Field<Vector>& gradP,
                   const Vector& bodyForce, scalar nu,
                   const std::vector<PatchBC>& bc, const FaceGeom& fg,
                   int sweeps, Field<scalar>& aDiag)
{
    const std::size_t nL = static_cast<std::size_t>(dm.nLocal());
    Field<scalar> aP(nL, 0.0);
    std::vector<Vector> rhs(nL, Vector::zero());
    struct Coupling { label peer; scalar coef; };
    std::vector<std::vector<Coupling>> links(nL);

    for (std::size_t k = 0; k < dm.faces().size(); ++k)
    {
        const DomainFace& df = dm.faces()[k];
        const std::size_t i = static_cast<std::size_t>(df.cell);
        const scalar Fi = df.sign * phi[k];           // outward from cell i
        const std::size_t f = static_cast<std::size_t>(df.face);
        const scalar coef = nu * fg.ortho[f];

        if (df.peer < 0)
        {
            const PatchBC& pb = bc[static_cast<std::size_t>(fg.patch[df.face])];
            // diffusion: fixedValue linearises to coef(Ui - Uface);
            // zeroGradient contributes nothing.
            if (pb.uType == "fixedValue")
            {
                aP[i] += coef;
                rhs[i] += pb.uValue * coef;
            }
            // convection: out of i -> self upwind; into i -> the BC value
            if (Fi >= 0) aP[i] += Fi;
            else rhs[i] += pb.uValue * (-Fi);
            continue;
        }
        // internal face: diffusion couples both ways, convection upwinds
        aP[i] += coef;
        links[i].push_back({df.peer, coef});
        if (Fi >= 0) aP[i] += Fi;
        else links[i].push_back({df.peer, -Fi});
    }

    const auto& vol = dm.cellVolumes();
    for (label i = 0; i < dm.nOwned(); ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        rhs[u] += (bodyForce - gradP[u]) * vol[u];
        aDiag[u] = aP[u];
    }

    for (int sweep = 0; sweep < sweeps; ++sweep)
    {
        for (label i = 0; i < dm.nOwned(); ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            Vector acc = rhs[u];
            for (const Coupling& L : links[u])
                acc += U[static_cast<std::size_t>(L.peer)] * L.coef;
            U[u] = acc / std::max(aP[u], 1e-300);
        }
    }
}

/// Face fluxes on the local face array: interior faces interpolate the
/// (haloed) velocity, boundary faces take the BC value, zero-normal patches
/// stay exactly zero. Both ranks of a coupled face compute the identical
/// double from identical inputs.
void faceFluxes(const DomainMesh& dm, const Field<Vector>& U,
                const std::vector<PatchBC>& bc, const FaceGeom& fg,
                Field<scalar>& phi)
{
    for (std::size_t k = 0; k < dm.faces().size(); ++k)
    {
        const DomainFace& df = dm.faces()[k];
        const Vector& sf = dm.faceSf()[k];
        if (df.peer < 0)
        {
            const PatchBC& pb = bc[static_cast<std::size_t>(fg.patch[df.face])];
            if (pb.zeroNormal) { phi[k] = 0.0; continue; }
            phi[k] = dot(uFaceValue(pb, U[static_cast<std::size_t>(df.cell)]), sf);
            continue;
        }
        label o = -1, n = -1;
        ownerNeighbour(df, o, n);
        const scalar w = dm.faceOwnerWeight()[k];
        const Vector uf = w * U[static_cast<std::size_t>(o)]
                        + (1.0 - w) * U[static_cast<std::size_t>(n)];
        phi[k] = dot(uf, sf);
    }
}

void cellDivPhi(const DomainMesh& dm, const Field<scalar>& phi, Field<scalar>& div)
{
    const std::size_t nL = static_cast<std::size_t>(dm.nLocal());
    std::fill(div.begin(), div.end(), 0.0);
    for (std::size_t k = 0; k < dm.faces().size(); ++k)
    {
        const DomainFace& df = dm.faces()[k];
        const std::size_t i = static_cast<std::size_t>(df.cell);
        if (i < nL) div[i] += df.sign * phi[k];
    }
}

/// The full-p Poisson: laplacian(rAUf, pFull) = div(phi) with the case's own
/// pressure boundary conditions, solved by local Gauss-Seidel with the halo
/// frozen between exchanges. Returns the owned-cell residual inf-norm.
scalar pressureSolve(const DomainMesh& dm, const Field<scalar>& rAU,
                     const Field<scalar>& divPhi,
                     const std::vector<PatchBC>& bc, const FaceGeom& fg,
                     const Field<scalar>& pLocal, Field<scalar>& pFull,
                     int sweeps)
{
    const std::size_t nL = static_cast<std::size_t>(dm.nLocal());
    std::vector<scalar> aP(nL, 0.0);
    std::vector<scalar> rhs(nL, 0.0);
    struct Coupling { label peer; scalar coef; };
    std::vector<std::vector<Coupling>> links(nL);

    for (std::size_t k = 0; k < dm.faces().size(); ++k)
    {
        const DomainFace& df = dm.faces()[k];
        const std::size_t i = static_cast<std::size_t>(df.cell);
        const std::size_t f = static_cast<std::size_t>(df.face);
        if (df.peer < 0)
        {
            const PatchBC& pb = bc[static_cast<std::size_t>(fg.patch[df.face])];
            if (pb.pType != "fixedValue") continue; // zeroGradient: no term
            const scalar coef = rAU[i] * fg.ortho[f];
            aP[i] += coef;
            rhs[i] += coef * pb.pValue;
            continue;
        }
        label o = -1, n = -1;
        ownerNeighbour(df, o, n);
        const scalar w = dm.faceOwnerWeight()[k];
        const scalar rAUf = w * rAU[static_cast<std::size_t>(o)]
                          + (1.0 - w) * rAU[static_cast<std::size_t>(n)];
        const scalar coef = rAUf * fg.ortho[f];
        aP[i] += coef;
        links[i].push_back({df.peer, coef});
    }
    for (label i = 0; i < dm.nOwned(); ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        rhs[u] += -divPhi[u];
        pFull[u] = pLocal[u]; // start the sweeps from the current pressure
    }

    scalar residual = 0;
    for (int sweep = 0; sweep < sweeps; ++sweep)
    {
        scalar res = 0;
        for (label i = 0; i < dm.nOwned(); ++i)
        {
            const std::size_t u = static_cast<std::size_t>(i);
            scalar acc = rhs[u];
            for (const Coupling& L : links[u])
                acc += L.coef * pFull[static_cast<std::size_t>(L.peer)];
            const scalar neu = acc / std::max(aP[u], 1e-300);
            res = std::max(res, std::abs(aP[u] * neu - acc));
            pFull[u] = neu;
        }
        residual = res;
    }
    return residual;
}

/// Rhie-Chow-style flux closure and the velocity response, mirroring
/// correctFluxes / the U update of simpleSolver.cpp on rank-local slots.
/// The closure is owner-oriented (vals[f] -= coef * (pN - pO)) on both sides
/// of a cut, so the two ranks write the identical double for one global face.
void correct(const DomainMesh& dm, const Field<scalar>& rAU,
             const Field<scalar>& pFull, const Field<Vector>& gradPFull,
             const Field<Vector>& gradPn, const std::vector<PatchBC>& bc,
             const FaceGeom& fg, const Field<scalar>& phiPredictor,
             Field<scalar>& phi, Field<Vector>& U, scalar relaxP)
{
    for (std::size_t k = 0; k < dm.faces().size(); ++k)
    {
        const DomainFace& df = dm.faces()[k];
        const std::size_t f = static_cast<std::size_t>(df.face);
        if (df.peer < 0)
        {
            const PatchBC& pb = bc[static_cast<std::size_t>(fg.patch[df.face])];
            if (pb.pType != "fixedValue") { phi[k] = phiPredictor[k]; continue; }
            const std::size_t i = static_cast<std::size_t>(df.cell);
            const scalar coef = rAU[i] * fg.ortho[f];
            phi[k] = phiPredictor[k] - coef * (pb.pValue - pFull[i]);
            continue;
        }
        label o = -1, n = -1;
        ownerNeighbour(df, o, n);
        const scalar w = dm.faceOwnerWeight()[k];
        const scalar rAUf = w * rAU[static_cast<std::size_t>(o)]
                          + (1.0 - w) * rAU[static_cast<std::size_t>(n)];
        const scalar coef = rAUf * fg.ortho[f];
        // non-orthogonal part: zero on an orthogonal mesh, like the serial
        // correctFluxes corr on this mesh class
        const scalar dp = pFull[static_cast<std::size_t>(n)]
                        - pFull[static_cast<std::size_t>(o)];
        phi[k] = phiPredictor[k] - coef * dp;
    }
    for (label i = 0; i < dm.nOwned(); ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        U[u] = U[u] - relaxP * rAU[u] * (gradPFull[u] - gradPn[u]);
    }
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int exitCode = 0;
    try
    {
        if (argc < 2)
            cfFatal("usage: chflow-3d-mpi-flow <case.dict> (run under mpiexec -np N)");
        const std::string casePath = argv[1];
        const Dictionary dict = Dictionary::parseFile(casePath);
        const polyMesh mesh = MeshLoader::fromCase(dict, casePath);
        const label nC = mesh.nCells();

        // ---- declared decomposition, refused before anything is built ----
        label declaredDomains = size;
        if (dict.has("mpi") && dict.getDict("mpi").has("domains"))
            declaredDomains = dict.getDict("mpi").getInt("domains");
        if (declaredDomains != size)
            cfFatal("mpiFlow: the case declares " + std::to_string(declaredDomains)
                    + " domains, which does not match the launched np=" +
                    std::to_string(size) + "; np and the decomposition must agree");
        if (size > nC)
            cfFatal("mpiFlow: " + std::to_string(size) + " domains for " + std::to_string(nC)
                    + " cells -- more domains than mesh segments, an empty domain "
                      "owns no face and its stencil degenerates");

        const scalar nu = dict.getScalar("nu");
        const FaceGeom fg = makeFaceGeom(mesh);

        // ---- boundary conditions, the simpleFoamApp reading ----
        std::vector<PatchBC> bc(static_cast<std::size_t>(mesh.nPatches()));
        const auto& patches = mesh.patches();
        for (label patch = 0; patch < mesh.nPatches(); ++patch)
        {
            const std::string& pname = patches[static_cast<std::size_t>(patch)].name;
            const std::string ptype = patches[static_cast<std::size_t>(patch)].type;
            PatchBC b;
            b.zeroNormal = ptype == "symmetry" || ptype == "empty";
            if (dict.has("patches"))
            {
                const Dictionary& ps = dict.getDict("patches");
                if (ps.has(pname))
                {
                    const Dictionary& pd = ps.getDict(pname);
                    if (pd.has("U"))
                    {
                        const auto ul = pd.getScalarList("U");
                        b.uType = pd.getWordOr("UType", "fixedValue");
                        b.uValue = Vector{ul[0], ul[1], ul[2]};
                    }
                    else if (ptype == "wall") { b.uType = "fixedValue"; }
                    if (pd.has("p"))
                    {
                        b.pType = pd.getWordOr("pType", "fixedValue");
                        b.pValue = pd.getScalar("p");
                    }
                }
            }
            bc[static_cast<std::size_t>(patch)] = b;
        }

        // ---- SIMPLE controls ----
        scalar relaxU = 0.7, relaxP = 0.3;
        scalar continuityTol = 1e-10, pcorrTol = 1e-12, momentumTol = 1e-12;
        int maxIterations = 400, gsSweeps = 60;
        Vector bodyForce = Vector::zero();
        std::string gradScheme = "linear";
        if (dict.has("simple"))
        {
            const Dictionary& sd = dict.getDict("simple");
            maxIterations = static_cast<int>(sd.getScalarOr("maxIterations", maxIterations));
            continuityTol = sd.getScalarOr("continuityTol", continuityTol);
            pcorrTol = sd.getScalarOr("pcorrTol", pcorrTol);
            momentumTol = sd.getScalarOr("momentumTol", momentumTol);
            relaxU = sd.getScalarOr("relaxU", relaxU);
            relaxP = sd.getScalarOr("relaxP", relaxP);
            gsSweeps = static_cast<int>(sd.getScalarOr("gsSweeps", gsSweeps));
            if (sd.has("force"))
            {
                const auto fv = sd.getScalarList("force");
                bodyForce = Vector{fv[0], fv[1], fv[2]};
            }
            gradScheme = sd.getWordOr("grad", gradScheme);
        }

        std::string dumpName = "flow.dat";
        if (dict.has("output"))
            dumpName = dict.getDict("output").getWordOr("dump", dumpName);

        const parallel::Decomposition dec(mesh, size);
        const parallel::DomainSet ds(dec);
        const parallel::DomainMesh& dm = ds[rank];

        scalar maxArea = 0;
        for (label f = 0; f < mesh.nFaces(); ++f)
            maxArea = std::max(maxArea, fg.magSf[static_cast<std::size_t>(f)]);
        maxArea = parallel::mpi::reduceMax(maxArea, comm);
        maxArea = std::max(maxArea, 1e-300);

        // rank-local state: owned slots then halo slots
        const std::size_t nL = static_cast<std::size_t>(dm.nLocal());
        Field<Vector> U(nL, Vector::zero());
        Field<scalar> p(nL, 0.0);
        Field<scalar> phi(dm.faces().size(), 0.0);
        Field<scalar> phiPredictor(dm.faces().size(), 0.0);
        Field<scalar> rAU(nL, 0.0);
        Field<scalar> aDiag(nL, 0.0);
        Field<Vector> grad(nL, Vector::zero());
        Field<Vector> gradPn(nL, Vector::zero());
        Field<Vector> gradPFull(nL, Vector::zero());
        Field<scalar> divPhi(nL, 0.0);
        Field<scalar> pFull(nL, 0.0);
        Field<scalar> pLocal(nL, 0.0);

        faceFluxes(dm, U, bc, fg, phi);

        if (rank == 0)
            log::info("mpiFlowSolver start ranks ", size, " cells ", nC,
                      " domains ", dec.nDomains(), " coupledFaces ", dec.nCoupledFaces(),
                      " nu ", num(nu), " relaxU ", num(relaxU), " relaxP ", num(relaxP),
                      " gsSweeps ", gsSweeps);

        bool converged = false;
        int iter = 0;
        scalar contErr = 1e300, dpMax = 1e300, uRes = 1e300;
        label haloMoved = 0;

        for (; iter < maxIterations; ++iter)
        {
            // ---------------- halo sync ----------------
            haloMoved += parallel::mpi::exchange(dec, rank, U, comm);
            haloMoved += parallel::mpi::exchange(dec, rank, p, comm);

            localGrad(dec, rank, comm, dm, p, bc, fg, gradScheme, grad);
            gradPn = grad;

            // ---------------- momentum ----------------
            const Field<Vector> UOld = U;
            momentumSolve(dm, U, phi, gradPn, bodyForce, nu, bc, fg, gsSweeps, aDiag);
            for (label i = 0; i < dm.nOwned(); ++i)
            {
                const std::size_t u = static_cast<std::size_t>(i);
                rAU[u] = dm.cellVolumes()[u] / std::max(aDiag[u], 1e-300);
            }
            // HbyA blend (relaxation against the old field)
            for (label i = 0; i < dm.nOwned(); ++i)
            {
                const std::size_t u = static_cast<std::size_t>(i);
                U[u] = UOld[u] + relaxU * (U[u] - UOld[u]);
            }

            // ---------------- predictor flux ----------------
            haloMoved += parallel::mpi::exchange(dec, rank, U, comm);
            faceFluxes(dm, U, bc, fg, phi);
            phiPredictor = phi;

            // adjustPhi, only when no patch carries a fixed pressure (the
            // serial needReference test), with the sums as Allreduce
            {
                label hasFixedP = 0;
                for (const PatchBC& b : bc) if (b.pType == "fixedValue") hasFixedP = 1;
                const label anyFixed = parallel::mpi::reduceSumInt(hasFixedP, comm);
                if (!anyFixed)
                {
                    scalar netSum = 0, freeSum = 0;
                    for (std::size_t k = 0; k < dm.faces().size(); ++k)
                    {
                        const DomainFace& df = dm.faces()[k];
                        if (df.peer >= 0) continue;
                        const PatchBC& pb = bc[static_cast<std::size_t>(fg.patch[df.face])];
                        if (pb.uType != "fixedValue") freeSum += phi[k];
                        netSum += phi[k];
                    }
                    const scalar net = parallel::mpi::reduceSum(netSum, comm);
                    const scalar free_ = parallel::mpi::reduceSum(freeSum, comm);
                    if (std::abs(free_) > 1e-30)
                    {
                        const scalar scale = (free_ - net) / free_;
                        for (std::size_t k = 0; k < dm.faces().size(); ++k)
                        {
                            const DomainFace& df = dm.faces()[k];
                            if (df.peer >= 0) continue;
                            const PatchBC& pb = bc[static_cast<std::size_t>(fg.patch[df.face])];
                            if (pb.uType == "fixedValue") continue;
                            phi[k] *= scale;
                        }
                    }
                }
            }

            cellDivPhi(dm, phi, divPhi); // predictor divergence: the Poisson RHS

            // ---------------- full-p Poisson ----------------
            haloMoved += parallel::mpi::exchange(dec, rank, rAU, comm);
            pLocal = p;
            pFull = pLocal; // halo slots carry the exchanged pressure the
                            // frozen-halo sweeps read as the old iterate
            const scalar pResLocal =
                pressureSolve(dm, rAU, divPhi, bc, fg, pLocal, pFull, gsSweeps);
            // the Poisson residual is a global collective: MAX over ranks
            // (order independent, exact); the serial convergence test's twin
            const scalar pRes = parallel::mpi::reduceMax(pResLocal, comm);

            // the closure must see the same pFull on both sides of a cut, so
            // the solved pressure is haloed before the flux correction
            Field<scalar> pFullHalo = pFull;
            haloMoved += parallel::mpi::exchange(dec, rank, pFullHalo, comm);
            localGrad(dec, rank, comm, dm, pFullHalo, bc, fg, gradScheme, gradPFull);

            correct(dm, rAU, pFullHalo, gradPFull, gradPn, bc, fg, phiPredictor,
                    phi, U, relaxP);

            // continuity on the corrected fluxes, the serial
            // continuityError(phi, U) equivalent, global MAX reduction
            cellDivPhi(dm, phi, divPhi);
            scalar contLocal = 0;
            for (label i = 0; i < dm.nOwned(); ++i)
                contLocal = std::max(contLocal, std::abs(divPhi[static_cast<std::size_t>(i)]));
            contErr = parallel::mpi::reduceMax(contLocal, comm) / maxArea;

            // ---------------- relax ----------------
            scalar stepLocal = 0, scaleLocal = 0;
            for (label i = 0; i < dm.nOwned(); ++i)
            {
                const std::size_t u = static_cast<std::size_t>(i);
                const scalar dp = relaxP * (pFull[u] - p[u]);
                stepLocal = std::max(stepLocal, std::abs(dp));
                scaleLocal = std::max(scaleLocal, std::abs(p[u]));
                p[u] += dp;
            }
            dpMax = parallel::mpi::reduceMax(stepLocal, comm);
            const scalar pScaleG = parallel::mpi::reduceMax(scaleLocal, comm);

            scalar uMaxLocal = 0, dULocal = 0;
            for (label i = 0; i < dm.nOwned(); ++i)
            {
                const std::size_t u = static_cast<std::size_t>(i);
                uMaxLocal = std::max(uMaxLocal, mag(U[u]));
                dULocal = std::max(dULocal, mag(U[u] - UOld[u]));
            }
            const scalar uMaxG = parallel::mpi::reduceMax(uMaxLocal, comm);
            const scalar dUG = parallel::mpi::reduceMax(dULocal, comm);
            uRes = dUG / std::max(uMaxG, 1e-300);

            if (rank == 0)
                log::info("mpiFlowSolver iter ", iter + 1, " cont ", num(contErr),
                          " dp ", num(dpMax), " dU ", num(uRes), " pRes ", num(pRes));

            if (contErr < continuityTol && dpMax < pcorrTol * std::max(pScaleG, 1.0)
                && uRes < momentumTol)
            {
                converged = true;
                ++iter;
                break;
            }
        }

        // every rank must have stopped on the same condition
        const int allConv = parallel::mpi::reduceSumInt(converged ? 1 : 0, comm);

        // final halo-consistent state for the dumps
        haloMoved += parallel::mpi::exchange(dec, rank, U, comm);
        haloMoved += parallel::mpi::exchange(dec, rank, p, comm);
        cellDivPhi(dm, phi, divPhi);

        // ---- per-rank halo dump + checksum: the exact values this rank's
        // Sendrecv received, so the checker can tie every halo slot to the
        // owner's owned slot bit for bit ----
        {
            const std::string base = dumpName.substr(0, dumpName.find_last_of('.'));
            std::ofstream hf(base + "_halo_np" + std::to_string(size)
                             + "_r" + std::to_string(rank) + ".dat");
            hf << std::setprecision(17);
            scalar checksum = 0;
            const std::vector<label>& halo = dm.haloCells();
            for (std::size_t h = 0; h < halo.size(); ++h)
            {
                const std::size_t slot = static_cast<std::size_t>(dm.nOwned()) + h;
                hf << halo[h] << ' ' << p[slot] << ' ' << U[slot].x << ' '
                   << U[slot].y << ' ' << U[slot].z << '\n';
                checksum += std::abs(p[slot]) + std::abs(U[slot].x)
                          + std::abs(U[slot].y) + std::abs(U[slot].z);
            }
            log::info("mpiFlowSolver halo rank ", rank, " slots ",
                      static_cast<label>(halo.size()), " checksum ", num(checksum));
        }

        // ---- per-rank face-flux dump (global face, owner, neighbour, phi) ----
        {
            const std::string base = dumpName.substr(0, dumpName.find_last_of('.'));
            std::ofstream pf(base + "_phi_np" + std::to_string(size)
                             + "_r" + std::to_string(rank) + ".dat");
            pf << std::setprecision(17);
            for (std::size_t k = 0; k < dm.faces().size(); ++k)
            {
                const DomainFace& df = dm.faces()[k];
                const Vector& fc = dm.faceCentres()[k];
                // global face, its owner/neighbour pair (-1 on patches), the
                // owner-oriented flux, this side's orientation sign and the
                // face centre so the checker can classify patches
                pf << df.face << ' ' << mesh.owner()[df.face] << ' '
                   << (df.face < mesh.nInternalFaces() ? mesh.neighbour()[df.face] : -1)
                   << ' ' << phi[k] << ' ' << df.sign << ' '
                   << fc.x << ' ' << fc.y << ' ' << fc.z << '\n';
            }
        }

        // ---- gather the fields and report ----
        const Field<scalar> pG = parallel::mpi::gatherScalar(dec, rank, p, comm);
        const Field<Vector> uG = parallel::mpi::gatherVector(dec, rank, U, comm);
        const scalar divMaxLocal = [&]()
        {
            scalar m = 0;
            for (label i = 0; i < dm.nOwned(); ++i)
                m = std::max(m, std::abs(divPhi[static_cast<std::size_t>(i)]));
            return m;
        }();
        const scalar divMax = parallel::mpi::reduceMax(divMaxLocal, comm);

        if (rank == 0)
        {
            log::info("mpiFlowSolver solve ranks ", size, " iters ", iter,
                      " converged ", (converged ? 1 : 0), " convRanks ", allConv, "/", size,
                      " contErr ", num(contErr), " dpMax ", num(dpMax),
                      " uRes ", num(uRes), " haloMoved ", haloMoved,
                      " divPhiMax ", num(divMax));

            const std::string base = dumpName.substr(0, dumpName.find_last_of('.'));
            std::ofstream out(base + "_np" + std::to_string(size) + ".dat");
            out << std::setprecision(17);
            for (label c = 0; c < nC; ++c)
            {
                const std::size_t u = static_cast<std::size_t>(c);
                const Vector& x = mesh.cellCentres()[u];
                out << c << ' ' << x.x << ' ' << x.y << ' ' << x.z << ' '
                    << pG[u] << ' ' << uG[u].x << ' ' << uG[u].y << ' ' << uG[u].z << '\n';
            }
        }

        if (!converged || allConv != size)
        {
            if (rank == 0) log::error("mpiFlowSolver: the distributed SIMPLE did not converge");
            exitCode = 1;
        }
    }
    catch (const std::exception& e)
    {
        log::error("mpiFlowSolver (rank ", rank, "): ", e.what());
        exitCode = 1;
    }

    int anyBad = 0;
    MPI_Allreduce(&exitCode, &anyBad, 1, MPI_INT, MPI_MAX, comm);
    MPI_Finalize();
    return anyBad ? 1 : 0;
}
