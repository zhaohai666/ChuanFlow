#include "cf/pvt/Viscosity.h"

#include <cmath>

#include "cf/common/Error.h"

namespace cf::pvt
{

namespace
{
constexpr scalar kAtm = 101325.0;      // Pa per atm
constexpr scalar kPrZc = 0.3074;       // PR76 universal critical compressibility

// xi = Tc^(1/6) / (sqrt(MW) Pc[atm]^(2/3)); Tc [K], MW [g/mol], Pc [Pa].
scalar viscosityReducingParameter(scalar Tc, scalar Pc, scalar MW)
{
    const scalar pcAtm = Pc / kAtm;
    return std::pow(Tc, 1.0 / 6.0)
           / (std::sqrt(MW) * std::pow(pcAtm, 2.0 / 3.0));
}
} // namespace

scalar stielThodosViscosity(scalar T, scalar Tc, scalar Pc, scalar MW)
{
    const scalar Tr = T / Tc;
    const scalar xi = viscosityReducingParameter(Tc, Pc, MW);
    const scalar muCp = Tr > 1.5
        ? 17.78e-5 * std::pow(4.58 * Tr - 1.67, 0.625) / xi
        : 34.0e-5 * std::pow(Tr, 0.94) / xi;
    return muCp * 1.0e-3; // cP -> Pa s
}

scalar lbcViscosity(scalar T, scalar Vm, const std::vector<scalar>& zs,
                    const std::vector<scalar>& MWs,
                    const std::vector<scalar>& Tcs,
                    const std::vector<scalar>& Pcs,
                    const std::vector<scalar>& Vcs)
{
    const std::size_t n = zs.size();
    if (n == 0 || MWs.size() != n || Tcs.size() != n || Pcs.size() != n
        || Vcs.size() != n)
        cfFatal("pvt: lbcViscosity component arrays disagree in length");
    if (!(Vm > 0.0)) cfFatal("pvt: lbcViscosity needs a positive molar volume");

    // Kay's rule: mole-fraction averaged critical constants and molar mass.
    scalar Tc = 0, MW = 0, Pc = 0, Vc = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        Tc += zs[i] * Tcs[i];
        MW += zs[i] * MWs[i];
        Pc += zs[i] * Pcs[i];
        Vc += zs[i] * Vcs[i];
    }
    const scalar rr = Vc / Vm; // reduced density rho/rhoc = Vc/Vm

    // Dilute-gas term: Herning-Zipperer sqrt(MW) mix of Stiel-Thodos mus.
    scalar num = 0, den = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const scalar root = std::sqrt(MWs[i]);
        num += zs[i] * stielThodosViscosity(T, Tcs[i], Pcs[i], MWs[i]) * root;
        den += zs[i] * root;
    }
    const scalar mu0Cp = (num / den) * 1.0e3; // Pa s -> cP

    const scalar xi = viscosityReducingParameter(Tc, Pc, MW);
    const scalar poly = 0.1023 + 0.023364 * rr + 0.058533 * rr * rr
                        - 0.040758 * rr * rr * rr
                        + 0.0093324 * rr * rr * rr * rr;
    const scalar muCp = mu0Cp + (std::pow(poly, 4.0) - 1.0e-4) / xi;
    return muCp * 1.0e-3; // cP -> Pa s
}

scalar lbcViscosity(const Mixture& mix, scalar T, scalar p, RootPreference pref)
{
    const EosState st = pengRobinsonState(mix, T, p, pref);
    const std::vector<Component>& comps = mix.components();
    const std::vector<scalar>& zs = mix.fractions();
    std::vector<scalar> MWs(comps.size()), Tcs(comps.size()),
        Pcs(comps.size()), Vcs(comps.size());
    for (std::size_t i = 0; i < comps.size(); ++i)
    {
        MWs[i] = comps[i].molarMass * 1.0e3; // kg/mol -> g/mol
        Tcs[i] = comps[i].Tc;
        Pcs[i] = comps[i].Pc;
        Vcs[i] = kPrZc * kGasConstant * comps[i].Tc / comps[i].Pc;
    }
    return lbcViscosity(T, st.molarVolume, zs, MWs, Tcs, Pcs, Vcs);
}

} // namespace cf::pvt