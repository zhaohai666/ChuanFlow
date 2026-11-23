#pragma once

// Shared case plumbing for the batch-3dturb apps (F4 turbulence + F12 phase
// change). Same `patches{}` / `piso{}` reading the piso app uses, kept local
// so the new apps never touch the existing ones.

#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Log.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/PisoSolver.h"
#include "cf/fvm/PolyMesh.h"
#include "cf/fvm/VtkWriter.h"

namespace cf
{
namespace app3dturb
{

/// 12 significant digits, the reporting precision every checker parses.
inline std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}

/// Apply the `patches{ name { U ...; UType ...; p ...; pType ...; } }` block
/// to the velocity and pressure fields, exactly as the piso app reads it.
inline void applyPatches(const Dictionary& caseDict, const polyMesh& mesh,
                         volVectorField& U, volScalarField& p)
{
    if (!caseDict.has("patches")) return;
    const Dictionary& ps = caseDict.getDict("patches");
    for (const std::string& name : ps.keys())
    {
        const label patch = mesh.patchIdByName(name);
        const Dictionary& pd = ps.getDict(name);
        if (pd.has("U"))
        {
            const auto ul = pd.getScalarList("U");
            const std::string type = pd.getWordOr("UType", "fixedValue");
            U.setPatchValue(patch, type, Vector{ul[0], ul[1], ul[2]});
        }
        else
        {
            U.setPatchValue(patch, pd.getWordOr("UType", "zeroGradient"), Vector::zero());
        }
        if (pd.has("p"))
        {
            p.setPatchValue(patch, pd.getWordOr("pType", "fixedValue"), pd.getScalar("p"));
        }
        else
        {
            p.setPatchValue(patch, "zeroGradient", 0.0);
        }
    }
}

/// The `piso{}` block: dt/endTime plus the coupling controls, read the way
/// the piso app reads them.
inline PisoControls pisoControls(const Dictionary& caseDict, scalar& endTime,
                                 int& reportEvery)
{
    PisoControls ctl;
    endTime = 1.0;
    reportEvery = 10;
    if (!caseDict.has("piso")) cfFatal("case needs a piso { dt; endTime; ... } block");
    const Dictionary& sd = caseDict.getDict("piso");
    ctl.dt = sd.getScalar("dt");
    endTime = sd.getScalar("endTime");
    ctl.correctors = static_cast<int>(sd.getScalarOr("correctors", ctl.correctors));
    ctl.continuityTol = sd.getScalarOr("continuityTol", ctl.continuityTol);
    ctl.linearTol = sd.getScalarOr("linearTol", ctl.linearTol);
    ctl.linearMaxIter = static_cast<int>(sd.getScalarOr("linearMaxIter", ctl.linearMaxIter));
    ctl.nonOrthCorrectors =
        static_cast<int>(sd.getScalarOr("nonOrthCorrectors", ctl.nonOrthCorrectors));
    reportEvery = static_cast<int>(sd.getScalarOr("reportEvery", reportEvery));
    if (sd.has("force"))
    {
        const auto f = sd.getScalarList("force");
        ctl.bodyForce = Vector{f[0], f[1], f[2]};
    }
    return ctl;
}

inline scalar caseNu(const Dictionary& caseDict)
{
    const scalar nu = caseDict.getScalarOr("nu", -1);
    if (!(nu > 0)) cfFatal("case needs a positive top-level nu");
    return nu;
}

inline void writeVtuIfDeclared(const Dictionary& caseDict, const polyMesh& mesh,
                               const volVectorField& U, const volScalarField& p)
{
    if (!caseDict.has("output")) return;
    const Dictionary& od = caseDict.getDict("output");
    const std::string vtk = od.getWordOr("vtk", "");
    if (vtk.empty()) return;
    VtkWriter w(mesh);
    w.addCellVector("U", U.internalField());
    w.addCellScalar("p", p.internalField());
    w.write(vtk);
}

} // namespace app3dturb
} // namespace cf
