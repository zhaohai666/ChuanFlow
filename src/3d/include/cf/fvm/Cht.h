#pragma once

#include <string>

#include "cf/common/Dictionary.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// D8 increment 7: conjugate heat transfer across two solid/fluid regions,
/// each its own polyMesh, joined at a shared interface (a "multiRegion" CHT
/// block). The two regions are solved as separate steady conduction problems;
/// the interface enforces temperature continuity and heat-flux continuity by a
/// Dirichlet--Neumann iteration on the interface temperature.
///
/// Steady conduction in each region (no flow; convection in the fluid region is
/// a later increment):
///   div(k grad T) = 0
/// with fixed temperature on the outer wall of each region and the unknown
/// interface temperature T_i on the shared face. Once both fields are solved,
/// the outward heat flux at the two interface faces must satisfy
///   q_s = k_s dT_s/dn  =  q_f = -k_f dT_f/dn   (equal and opposite normals)
/// and the iteration adjusts T_i until |q_s - q_f| is below tolerance.
///
/// For a 1D slab of two layers with thicknesses L_s, L_f and conductivities
/// k_s, k_f between outer temperatures T_L and T_R, the closed form is
///   R    = L_s/k_s + L_f/k_f
///   q    = (T_L - T_R) / R                 [W/m2]
///   T_i  = T_L - q L_s / k_s
/// which the checker reconstructs from the case constants, together with the
/// linear profile in each layer and the exact flux continuity at the interface.
struct RegionConfig
{
    scalar k = 0;          // conductivity [W/m K]
    scalar outerT = 0;     // fixed temperature on this region's outer wall [K]
    scalar initialT = 0;   // initial interface guess seed [K]
};

struct ChtConfig
{
    bool enabled = false;
    RegionConfig solid;
    RegionConfig fluid;
    scalar interfaceTol = 1e-9;  // relative flux-balance tolerance
    int maxIterations = 50;
    scalar relax = 1.0;          // relaxation of the interface temperature
    void validate() const;
};

ChtConfig readChtConfig(const Dictionary& caseDict);

struct RegionResult
{
    Field<scalar> T;         // per cell
    scalar interfaceFlux = 0; // conductive flux at the interface face [W/m2]
    scalar outerFlux = 0;     // flux at the outer wall [W/m2]
    int iterations = 0;
    scalar residual = 0;
    bool solved = false;
};

struct ChtResult
{
    RegionResult solid;
    RegionResult fluid;
    scalar interfaceTemperature = 0;  // converged T_i [K]
    scalar fluxMismatch = 0;          // |q_s - q_f| / scale at convergence
    int outerIterations = 0;
    bool solved = false;
};

/// Solve one region's steady conduction with a fixed temperature on patch
/// `outerPatch` and `interfacePatch` set to `interfaceT`. Both are Dirichlet.
RegionResult solveRegion(const polyMesh& mesh, const Schemes& schemes,
                         scalar k, scalar outerT, scalar interfaceT,
                         const std::string& outerPatch,
                         const std::string& interfacePatch);

/// Converge the two regions on the shared interface temperature.
ChtResult solveCht(const polyMesh& solidMesh, const polyMesh& fluidMesh,
                   const Schemes& schemes, const ChtConfig& cfg,
                   const std::string& solidOuter, const std::string& solidIface,
                   const std::string& fluidOuter, const std::string& fluidIface);

} // namespace cf
