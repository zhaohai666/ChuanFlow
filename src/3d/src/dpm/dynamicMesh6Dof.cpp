#include <cmath>

#include "cf/common/Error.h"
#include "cf/fvm/DynamicMesh6Dof.h"

namespace cf
{

void DynamicMeshConfig::validate() const
{
    if (!std::isfinite(mass) || mass <= 0)
        cfFatal("dynamicMesh: mass must be finite and positive");
    if (!std::isfinite(drag) || drag < 0)
        cfFatal("dynamicMesh: drag must be finite and non-negative");
    if (!std::isfinite(inertia) || inertia <= 0)
        cfFatal("dynamicMesh: inertia must be finite and positive");
    if (!std::isfinite(rotDrag) || rotDrag < 0)
        cfFatal("dynamicMesh: rotDrag must be finite and non-negative");
    if (!std::isfinite(rotStiffness) || rotStiffness < 0)
        cfFatal("dynamicMesh: rotStiffness must be finite and non-negative");
    if (!std::isfinite(dt) || !(dt > 0))
        cfFatal("dynamicMesh: dt must be finite and positive");
    if (steps < 1)
        cfFatal("dynamicMesh: steps must be at least 1");
    if (reportEvery < 1)
        cfFatal("dynamicMesh: reportEvery must be at least 1");
    for (scalar c : {meshVelocity.x, meshVelocity.y, meshVelocity.z,
                     bodyForce.x, bodyForce.y, bodyForce.z,
                     fluidVelocity.x, fluidVelocity.y, fluidVelocity.z,
                     theta0, omega0, initialValue})
        if (!std::isfinite(c))
            cfFatal("dynamicMesh: all inputs must be finite");
}

void FluidStressConfig::validate() const
{
    if (!std::isfinite(density) || density <= 0)
        cfFatal("fluidStress: density must be finite and positive");
    if (!std::isfinite(viscosity) || viscosity < 0)
        cfFatal("fluidStress: viscosity must be finite and non-negative");
    for (scalar c : {gravity.x, gravity.y, gravity.z, pressureRef})
        if (!std::isfinite(c))
            cfFatal("fluidStress: gravity and pressureRef must be finite");
    if (!velocityGradient.empty() && velocityGradient.size() != 9)
        cfFatal("fluidStress: velocityGradient needs exactly 9 components");
    for (scalar g : velocityGradient)
        if (!std::isfinite(g))
            cfFatal("fluidStress: velocityGradient entries must be finite");
    if (surface.empty()) return;
    for (const std::string& s : surface)
        if (s.empty())
            cfFatal("fluidStress: surface patch names must be non-empty");
}

DynamicMeshConfig readDynamicMeshConfig(const Dictionary& caseDict)
{
    DynamicMeshConfig cfg;
    if (!caseDict.has("dynamicMesh")) return cfg;
    const Dictionary& d = caseDict.getDict("dynamicMesh");
    cfg.enabled = true;

    const auto vec3 = [&d](const std::string& key, const Vector& def) {
        if (!d.has(key)) return def;
        const auto v = d.getScalarList(key);
        if (v.size() != 3)
            cfFatal("dynamicMesh: " + key + " needs three components");
        return Vector{v[0], v[1], v[2]};
    };
    cfg.meshVelocity = vec3("meshVelocity", Vector::zero());
    cfg.bodyForce = vec3("bodyForce", Vector::zero());
    cfg.fluidVelocity = vec3("fluidVelocity", Vector::zero());

    cfg.mass = d.getScalarOr("mass", 0.0);
    cfg.drag = d.getScalarOr("drag", 0.0);
    cfg.inertia = d.getScalarOr("inertia", 0.0);
    cfg.rotDrag = d.getScalarOr("rotDrag", 0.0);
    cfg.rotStiffness = d.getScalarOr("rotStiffness", 0.0);
    cfg.theta0 = d.getScalarOr("theta0", 0.0);
    cfg.omega0 = d.getScalarOr("omega0", 0.0);
    cfg.initialValue = d.getScalarOr("initialValue", 0.0);
    cfg.dt = d.getScalarOr("dt", 0.0);
    cfg.steps = static_cast<int>(d.getScalarOr("steps", 0.0));
    cfg.reportEvery = static_cast<int>(d.getScalarOr("reportEvery", 1.0));
    cfg.validate();

    // §5.236 opt-in fluid-stress block. Absent -> cfg.stress.enabled stays
    // false and every downstream code path is the untouched F11 reference.
    if (caseDict.has("fluidStress"))
    {
        const Dictionary& s = caseDict.getDict("fluidStress");
        cfg.stress.enabled = true;
        cfg.stress.density = s.getScalarOr("density", 0.0);
        cfg.stress.viscosity = s.getScalarOr("viscosity", 0.0);
        cfg.stress.pressureRef = s.getScalarOr("pressureRef", 0.0);
        if (s.has("gravity"))
        {
            const auto g = s.getScalarList("gravity");
            if (g.size() != 3)
                cfFatal("fluidStress: gravity needs three components");
            cfg.stress.gravity = Vector{g[0], g[1], g[2]};
        }
        if (s.has("velocityGradient"))
            cfg.stress.velocityGradient = s.getScalarList("velocityGradient");
        if (s.has("surface"))
            cfg.stress.surface = s.getWordList("surface");
        cfg.stress.validate();
    }
    return cfg;
}

Field<Vector> meshFaceVelocity(const polyMesh& mesh, const Vector& Umesh)
{
    const auto& faces = mesh.faces();
    Field<Vector> v(static_cast<std::size_t>(mesh.nFaces()), Vector::zero());
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const auto& verts = faces[static_cast<std::size_t>(f)];
        // Uniform mesh velocity: the vertex average is Umesh itself.
        if (!verts.empty())
            v[static_cast<std::size_t>(f)] = Umesh;
    }
    return v;
}

Field<scalar> relativeFaceFlux(const polyMesh& mesh, const Vector& U, const Vector& Umesh)
{
    const Vector Urel = U - Umesh;
    const auto& sf = mesh.faceSf();
    Field<scalar> flux(static_cast<std::size_t>(mesh.nFaces()), 0.0);
    for (label f = 0; f < mesh.nFaces(); ++f)
        flux[static_cast<std::size_t>(f)] = dot(Urel, sf[static_cast<std::size_t>(f)]);
    return flux;
}

Field<scalar> gclGeometricSweep(const polyMesh& mesh, const Vector& Umesh)
{
    Field<scalar> sweep(static_cast<std::size_t>(mesh.nCells()), 0.0);
    const auto& sf = mesh.faceSf();
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const scalar flux = dot(Umesh, sf[u]);
        sweep[static_cast<std::size_t>(mesh.owner()[u])] += flux;
        if (f < mesh.nInternalFaces())
            sweep[static_cast<std::size_t>(mesh.neighbour()[u])] -= flux;
    }
    return sweep;
}

Field<scalar> gclCellResidual(const polyMesh& mesh, const Vector& Umesh,
                              const Field<scalar>& oldVolumes, scalar dt)
{
    if (!(dt > 0)) cfFatal("dynamicMesh: GCL needs a positive dt");
    const Field<scalar> sweep = gclGeometricSweep(mesh, Umesh);
    const auto& vol = mesh.cellVolumes();
    Field<scalar> res(static_cast<std::size_t>(mesh.nCells()), 0.0);
    for (label c = 0; c < mesh.nCells(); ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        res[u] = (vol[u] - oldVolumes[u]) / dt - sweep[u];
    }
    return res;
}

Field<scalar> transportStep(const polyMesh& mesh, const Field<scalar>& phi,
                            const Field<scalar>& phiRel, const Field<scalar>& sweep,
                            const Field<scalar>& oldVolumes, scalar dt)
{
    if (!(dt > 0)) cfFatal("dynamicMesh: transport needs a positive dt");
    const label nC = mesh.nCells();
    Field<scalar> newPhi(static_cast<std::size_t>(nC), 0.0);

    // Relative-flux divergence, assembled owner + / neighbour -. Face values
    // are upwinded by the sign of the relative face flux (zero-gradient at
    // boundaries). A uniform field is independent of this choice.
    Field<scalar> divRel(static_cast<std::size_t>(nC), 0.0);
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        const scalar fl = phiRel[u];
        const label owner = mesh.owner()[u];
        scalar faceValue = phi[static_cast<std::size_t>(owner)];
        if (f < mesh.nInternalFaces())
        {
            const label nb = mesh.neighbour()[u];
            faceValue = fl >= 0 ? phi[static_cast<std::size_t>(owner)]
                                : phi[static_cast<std::size_t>(nb)];
            divRel[static_cast<std::size_t>(nb)] -= faceValue * fl;
        }
        divRel[static_cast<std::size_t>(owner)] += faceValue * fl;
    }

    // Conservative relative-flux update plus the GCL correction. The
    // correction replaces the discrete geometric volume sweep dt*sweep by the
    // exact volume change (vol - oldVol), i.e. it adds phi_c*[(vol-oldVol)
    // - dt*sweep]. For a uniform phi this makes the update exact regardless
    // of the mesh-motion truncation error.
    const auto& vol = mesh.cellVolumes();
    for (label c = 0; c < nC; ++c)
    {
        const std::size_t u = static_cast<std::size_t>(c);
        const scalar correction = (vol[u] - oldVolumes[u]) - dt * sweep[u];
        const scalar massNew = phi[u] * oldVolumes[u] - dt * divRel[u] + phi[u] * correction;
        newPhi[u] = massNew / vol[u];
    }
    return newPhi;
}

BodyState advanceBody(const BodyState& s, const DynamicMeshConfig& cfg, scalar dt)
{
    if (!(dt > 0)) cfFatal("dynamicMesh: body advance needs a positive dt");
    const scalar F = cfg.bodyForce.x; // translational DOF is along the force direction

    struct Deriv { scalar x, v, theta, omega; };
    const auto deriv = [&](const BodyState& q) {
        return Deriv{
            q.v,
            (F - cfg.drag * q.v) / cfg.mass,
            q.omega,
            (-cfg.rotDrag * q.omega - cfg.rotStiffness * q.theta) / cfg.inertia,
        };
    };
    const auto step = [](BodyState& q, const Deriv& d, scalar h) {
        q.x += h * d.x;
        q.v += h * d.v;
        q.theta += h * d.theta;
        q.omega += h * d.omega;
    };

    const Deriv k1 = deriv(s);
    BodyState s2 = s; step(s2, k1, dt / 2);
    const Deriv k2 = deriv(s2);
    BodyState s3 = s; step(s3, k2, dt / 2);
    const Deriv k3 = deriv(s3);
    BodyState s4 = s; step(s4, k3, dt);
    const Deriv k4 = deriv(s4);

    BodyState out = s;
    out.x += dt / 6 * (k1.x + 2 * k2.x + 2 * k3.x + k4.x);
    out.v += dt / 6 * (k1.v + 2 * k2.v + 2 * k3.v + k4.v);
    out.theta += dt / 6 * (k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta);
    out.omega += dt / 6 * (k1.omega + 2 * k2.omega + 2 * k3.omega + k4.omega);
    return out;
}

StressLoad integrateFluidStress(const polyMesh& mesh, const FluidStressConfig& cfg)
{
    // Select the body-surface faces: the named patches, or the whole mesh
    // boundary when none are declared (a block mesh boundary is a closed box).
    std::vector<char> selected(static_cast<std::size_t>(mesh.nFaces()), 0);
    if (cfg.surface.empty())
    {
        for (label f = mesh.nInternalFaces(); f < mesh.nFaces(); ++f)
            selected[static_cast<std::size_t>(f)] = 1;
    }
    else
    {
        for (const std::string& name : cfg.surface)
        {
            const label pid = mesh.patchIdByName(name); // fatal on unknown name
            const polyMesh::Patch& pa = mesh.patches()[static_cast<std::size_t>(pid)];
            for (label f = pa.start; f < pa.start + pa.size; ++f)
                selected[static_cast<std::size_t>(f)] = 1;
        }
    }

    // T = grad U + grad U^T from the declared row-major velocity gradient.
    scalar T[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    if (!cfg.velocityGradient.empty())
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                const scalar g = cfg.velocityGradient[static_cast<std::size_t>(i * 3 + j)];
                T[i][j] += cfg.viscosity * (g + cfg.velocityGradient[static_cast<std::size_t>(j * 3 + i)]);
            }

    const auto& sf = mesh.faceSf();
    const auto& fc = mesh.faceCentres();
    StressLoad load;
    for (label f = 0; f < mesh.nFaces(); ++f)
    {
        const std::size_t u = static_cast<std::size_t>(f);
        if (!selected[u]) continue;
        const Vector& n = sf[u];
        const Vector& x = fc[u];
        // Hydrostatic closure: grad p = rho g, so p(x) = pRef + rho g.x.
        const scalar p = cfg.pressureRef + cfg.density * dot(cfg.gravity, x);
        const Vector fp = n * (-p);
        Vector fv = Vector::zero();
        if (!cfg.velocityGradient.empty())
            fv = Vector{T[0][0] * n.x + T[0][1] * n.y + T[0][2] * n.z,
                        T[1][0] * n.x + T[1][1] * n.y + T[1][2] * n.z,
                        T[2][0] * n.x + T[2][1] * n.y + T[2][2] * n.z};
        load.forcePressure += fp;
        load.forceViscous += fv;
        load.moment += cross(x, fp + fv);
        load.normalSum += n;
        load.volume += dot(x, n) / 3.0;
        load.area += mag(n);
    }

    // Loud refusal of a non-closed surface: the area-vector sum must vanish
    // and the divergence-theorem volume must be positive (closed + outward).
    if (mag(load.normalSum) > 1e-9 * std::max(1.0, load.area))
        cfFatal("fluidStress: body surface is not closed (sum of area vectors does not vanish)");
    if (!(load.volume > 0))
        cfFatal("fluidStress: body surface encloses a non-positive volume (face set does not match a closed body)");
    return load;
}

BodyState advanceBodyStress(const BodyState& s, const DynamicMeshConfig& cfg,
                            const StressLoad& load, scalar dt)
{
    if (!(dt > 0)) cfFatal("dynamicMesh: body advance needs a positive dt");
    // The integrated stress replaces the declarative bodyForce / drag pair;
    // the solid rotDrag / rotStiffness closures remain as mechanical terms.
    const scalar F = load.forcePressure.x + load.forceViscous.x;
    const scalar M = load.moment.z;

    struct Deriv { scalar x, v, theta, omega; };
    const auto deriv = [&](const BodyState& q) {
        return Deriv{
            q.v,
            F / cfg.mass,
            q.omega,
            (M - cfg.rotDrag * q.omega - cfg.rotStiffness * q.theta) / cfg.inertia,
        };
    };
    const auto step = [](BodyState& q, const Deriv& d, scalar h) {
        q.x += h * d.x;
        q.v += h * d.v;
        q.theta += h * d.theta;
        q.omega += h * d.omega;
    };

    const Deriv k1 = deriv(s);
    BodyState s2 = s; step(s2, k1, dt / 2);
    const Deriv k2 = deriv(s2);
    BodyState s3 = s; step(s3, k2, dt / 2);
    const Deriv k3 = deriv(s3);
    BodyState s4 = s; step(s4, k3, dt);
    const Deriv k4 = deriv(s4);

    BodyState out = s;
    out.x += dt / 6 * (k1.x + 2 * k2.x + 2 * k3.x + k4.x);
    out.v += dt / 6 * (k1.v + 2 * k2.v + 2 * k3.v + k4.v);
    out.theta += dt / 6 * (k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta);
    out.omega += dt / 6 * (k1.omega + 2 * k2.omega + 2 * k3.omega + k4.omega);
    return out;
}

} // namespace cf
