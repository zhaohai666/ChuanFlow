#include "cf/1d/Network.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

#include "cf/common/Error.h"
#include "cf/common/Log.h"

#include "cf/1d/BeggsBrill.h"
#include "cf/1d/GasLiftValve.h"
#include "cf/1d/RodPump.h"
#include "cf/1d/Separator.h"
#include "cf/1d/SolverGates.h"

#include "cf/pvt/Components.h"
#include "cf/pvt/ComponentData.h"
#include "cf/pvt/Flash.h"
#include "cf/pvt/HydrateEquilibrium.h"
#include "cf/pvt/PengRobinson.h"
#include "cf/pvt/PvtCase.h"
#include "cf/pvt/Srk.h"
#include "cf/pvt/SurfaceTension.h"
#include "cf/pvt/Viscosity.h"

#include "DenseSolve.h"

namespace cf
{
namespace
{

constexpr double kGravity = 9.80665;
// 5.86 (D1 increment 2): Turpin's correlation takes the intake pressure in
// psia. 1 psi = 6894.757293168 Pa exactly (lbf/in^2), so psia = Pa / this.
// Ali et al. (Energies 2021) print the same law as 0.145*p_in(kPa), a rounded
// kPa->psi factor; we use the exact one and the checker re-derives it here.
constexpr double kPascalPerPsi = 6894.757293168;

bool pipeIsAnnulus(const NetworkPipe& p) { return p.annulusInner > 0; }

bool pipeIsRect(const NetworkPipe& p) { return p.ductHeight > 0; }

// 5.89 (D1 increment 3): the dome charge is what switches a gas-lift valve on.
bool pipeIsLiftValve(const NetworkPipe& p) { return p.liftDomePressure > 0; }

// 5.90 (D1 increment 4): the barrel bore is what switches a rod pump on.
bool pipeIsRodPump(const NetworkPipe& p) { return p.rodPlungerDiameter > 0; }

// 5.223 (H3): the separator{} block itself is the switch.
bool pipeIsSeparator(const NetworkPipe& p) { return p.separatorOn; }

double pipeArea(const NetworkPipe& p)
{
    if (pipeIsRect(p)) return p.diameter * p.ductHeight;
    if (pipeIsAnnulus(p))
        return 0.25 * M_PI * (p.diameter * p.diameter
                              - p.annulusInner * p.annulusInner);
    return 0.25 * M_PI * p.diameter * p.diameter;
}

/// Hydraulic diameter 4A/P: for a concentric annulus exactly dOuter - dInner,
/// for a rectangle 2WH/(W+H), for a circle the diameter itself. The Reynolds
/// number, the relative roughness and the Darcy conversion L/Dh all read this
/// -- which is why a circular pipe never moves (Dh == diameter, bit for bit).
double pipeHydraulicDiameter(const NetworkPipe& p)
{
    if (pipeIsRect(p))
        return 2.0 * p.diameter * p.ductHeight / (p.diameter + p.ductHeight);
    return pipeIsAnnulus(p) ? p.diameter - p.annulusInner : p.diameter;
}

/// A number quoted in a fatal message, the way the user would write it: a
/// sentinel like the unresolved-grid spread prints as 1.79769e+308 rather than
/// the 309 digits std::to_string would spill.
std::string msgNum(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.6g", v);
    return std::string(buf);
}

/// Eccentric-annulus cross-section integral K_geom (5.87, C3 increment 2):
/// the area integral
///     K_geom = integral_A w dA
/// of the laminar Poisson solution grad^2 w = -1 with w = 0 on both bores, for
/// the gap between the outer circle (radius `ro`, centred on the pipe axis)
/// and the inner one (radius `ri`) whose axis is offset by `ecc` metres. It is
/// the whole eccentric correction to the flow law: the conductance is
/// K_geom/(mu L) and the laminar constant is 2 Dh^2 A / K_geom, with A and Dh
/// themselves eccentricity-independent (both bores stay fully wetted, so
/// neither the area pi(Ro^2 - Ri^2) nor the wetted perimeter pi(Do + Di)
/// moves -- only the profile does).
///
/// What makes the problem tractable is that BOTH circles are coordinate lines
/// of the bipolar map
///     x = a sinh tau /(cosh tau - cos sigma),
///     y = a sin  sigma/(cosh tau - cos sigma),
/// whose tau = const curves are circles of radius a/|sinh tau| centred at
/// a coth tau on the x axis. With the bipolar pole `x0` from the outer centre,
///     x0 = (e^2 + Ro^2 - Ri^2)/(2 e),   a^2 = x0^2 - Ro^2 = (x0 - e)^2 - Ri^2,
/// the two circles are tau_o with coth tau_o = -x0/a and tau_i with
/// coth tau_i = (e - x0)/a -- both negative, since each circle encloses the
/// opposite pole. Solving that pair through the two difference-of-squares
/// identities (x0 - a)(x0 + a) = Ro^2 and (x0 - e - a)(x0 - e + a) = Ri^2 keeps
/// every quantity O(radius):
///     tau_o = ln(Ro/(x0 + a)),   tau_i = ln(Ri/(x0 - e + a)),
/// which is what survives as e -> 0, where x0 ~ 1/e and coth tau -> -1.
///
/// In those coordinates (scale factor h = a/(cosh tau - cos sigma), and the
/// Laplacian is conformal so it picks up only h^2) the problem is
///     (d^2/dtau^2 + d^2/dsigma^2) w = -h^2
/// on the RECTANGLE [tau_i, tau_o] x [0, 2pi): Dirichlet at both tau ends,
/// periodic in sigma. Operator, boundary and right-hand side are all separable
/// in sigma, and the spectrum of h^2 is analytic: with T = |tau|, S = sinh T,
/// C = cosh T, the Poisson kernel integral and its tau-derivative give
///     h^2 = a^2 [ C/S^3 + sum_{m>=1} 2 e^{-mT}(C + m S)/S^3 cos m sigma ],
/// so each cosine mode m solves an independent ODE
///     W_m'' - m^2 W_m = -a^2 A_m,     W_m = 0 at both tau ends.
/// On the UNIFORM tau grid -- which is what bipolar coordinates buy, constant
/// spacing on both boundaries of the strip -- that is a tridiagonal system in
/// which only the diagonal depends on m, so one Thomas sweep per mode solves
/// the discrete system exactly (no iteration, no relaxation schedule, no
/// initial guess, and so no run-to-run path dependence). Parseval then closes
/// the sigma integral of w h^2:
///     K_geom = pi a^2 integral [ 2 A_0 W_0 + sum_{m>=1} A_m W_m ] dtau,
/// trapezoid in tau; the integrand vanishes at both ends, where w = 0, so the
/// end points contribute nothing and the rule stays second order.
///
/// Only tau carries discretisation error, and it is cleanly second order (the
/// measured convergence ratio is 4.00 at every eccentricity tried, e/gap 0.002
/// through 0.98), so two grids Richardson-extrapolate. The mode count follows
/// from the slowest-decaying harmonic across the strip, exp(-m |tau|min), and
/// the truncation is verified against that tail instead of assumed. The
/// disagreement of the two grids is the reported error estimate: the narrow gap
/// of a near-contact annulus is a boundary layer a fixed grid cannot resolve,
/// and the caller refuses loudly rather than pricing it with a number that is
/// quietly off.
///
/// Returns K_geom (m^4), or -1 when the geometry is unresolvable (eccentricity
/// at or past contact, the sigma tail not closed within the mode cap, or the
/// tau spread above the tolerance). `spread` carries the relative Richardson
/// disagreement of the two grids. As e -> 0 this reproduces the Fredrickson-
/// Bird closed form below to the extrapolation error (verified at 1e-5
/// relative), and addPipe does not run it there: e == 0 stays on the concentric
/// branch, bit for bit.
double eccentricAnnulusSectionIntegral
(
    double ro,
    double ri,
    double ecc,
    double& spread
)
{
    spread = std::numeric_limits<double>::max();
    const double x0 = (ecc*ecc + ro*ro - ri*ri)/(2.0*ecc);
    const double a2 = x0*x0 - ro*ro;
    if (!(a2 > 0.0))
    {
        return -1.0; // e >= Ro - Ri: the bores touch, the film pins the flow
    }
    const double a = std::sqrt(a2);
    const double tauOuter = std::log(ro/(x0 + a));
    const double tauInner = std::log(ri/(x0 - ecc + a));
    const double tlo = std::min(tauOuter, tauInner);
    const double thi = std::max(tauOuter, tauInner);
    const double tmin = std::min(std::abs(tlo), std::abs(thi));

    // Modes: the m-th harmonic of h^2 is down by exp(-m tmin) at the strip's
    // narrowest |tau|, so 45/tmin of them push the tail under machine zero. The
    // amplitude ratio A_m/A_0 = 2 exp(-m tmin)(1 + m tanh tmin) is then checked
    // AT the cut, so a geometry that needs more modes than the cap is refused
    // rather than silently truncated.
    const double tailNeedle = 45.0/tmin;
    if (!(tailNeedle < 20000.0))
    {
        return -1.0;
    }
    const int mmax = static_cast<int>(std::ceil(tailNeedle)) + 4;
    const double tail = 2.0*std::exp(-double(mmax)*tmin)
                      *(1.0 + double(mmax)*std::tanh(tmin));
    if (!(tail < 1e-14))
    {
        return -1.0;
    }

    auto sweep = [&](int nt) -> double
    {
        const double dtau = (thi - tlo)/double(nt);
        const int n = nt - 1; // interior rows; w = 0 at both ends
        std::vector<double> amp(n), cp(n), dp(n);
        double acc = 0.0;
        for (int m = 0; m <= mmax; ++m)
        {
            const double md = double(m)*dtau;
            const double diag = -(2.0 + md*md);
            for (int i = 0; i < n; ++i)
            {
                const double T = std::abs(tlo + double(i + 1)*dtau);
                const double S = std::sinh(T);
                const double C = std::cosh(T);
                amp[i] = (m == 0)? C/(S*S*S)
                    : 2.0*std::exp(-double(m)*T)*(C + double(m)*S)/(S*S*S);
                const double rhs = -a2*amp[i]*dtau*dtau;
                if (i == 0)
                {
                    cp[0] = 1.0/diag;
                    dp[0] = rhs/diag;
                }
                else
                {
                    const double beta = diag - cp[i - 1];
                    if (i + 1 < n)
                    {
                        cp[i] = 1.0/beta;
                    }
                    dp[i] = (rhs - dp[i - 1])/beta;
                }
            }
            // Parseval: the mean mode carries twice the weight of a cosine
            const double wt = (m == 0)? 2.0 : 1.0;
            double prev = dp[n - 1];
            acc += wt*amp[n - 1]*prev;
            for (int i = n - 2; i >= 0; --i)
            {
                const double w = dp[i] - cp[i]*prev;
                prev = w;
                acc += wt*amp[i]*w;
            }
        }
        return 3.14159265358979323846*a2*dtau*acc;
    };

    const double kc = sweep(400);
    const double kf = sweep(800);
    if (!(kf > 0.0))
    {
        return -1.0;
    }
    spread = std::abs(kf - kc)/kf;
    if (!(spread < 0.01))
    {
        return -1.0;
    }
    return kf + (kf - kc)/3.0;
}

/// Fredrickson-Bird (1958) Newtonian laminar anchor of a concentric annulus,
/// the exact biharmonic solution
///     Q = pi dp /(8 mu L) . [Ro^4 - Ri^4 - (Ro^2 - Ri^2)^2 / ln(Ro/Ri)],
/// which collapses onto Hagen-Poiseuille as kappa = Ri/Ro -> 0. Expressed as
/// a conductance so it drops into the same slot the circle's does.
///
/// 5.87 (C3 increment 2): an eccentric bore has no closed form, so when
/// addPipe resolved one it stores the SAME quantity -- the area integral
/// K_geom, which is pi g/8 concentrically -- and this returns K_geom/(mu L).
/// The two statements coincide at K_geom = pi g/8, which is what the eccentric
/// solve is asserted to approach as e -> 0. With ecc == 0 nothing here changes
/// by one bit.
double annulusLaminarConductance(const NetworkPipe& p,
                                 const FluidProperties& fluid)
{
    if (p.annulusEccK > 0.0)
    {
        return p.annulusEccK/(fluid.mu*p.length);
    }
    const double ro = 0.5 * p.diameter, ri = 0.5 * p.annulusInner;
    const double g = ro * ro * ro * ro - ri * ri * ri * ri
                   - (ro * ro - ri * ri) * (ro * ro - ri * ri)
                   / std::log(ro / ri);
    return M_PI * g / (8.0 * fluid.mu * p.length);
}

/// The laminar constant of the annulus, C(kappa) = f . Re_Dh, from the same
/// closed form: f = 2 Dh dp /(L rho v^2) and Re_Dh = rho v Dh / mu give
///     f Re = 16 Dh^2 (Ro^2 - Ri^2) / [Ro^4 - Ri^4 - (Ro^2-Ri^2)^2/ln(Ro/Ri)],
/// which is 64 at kappa = 0 -- so Churchill's laminar line is the kappa = 0
/// member of this family, and the all-regime blend below only has to swap
/// that member.
///
/// 5.87: eccentrically the bracket no longer exists, but the definition does.
/// Q = (dp/L) K_geom/mu with v = Q/A and Re_Dh = rho v Dh/mu gives, for ANY
/// bore shape, C = 2 Dh^2 A / K_geom -- which at K_geom = pi g/8 IS the closed
/// form above (16 Dh^2(Ro^2-Ri^2)/(pi g/8) reduces to it term for term). So the
/// eccentric branch reuses the geometry this law already trusts, the area
/// pi(Ro^2-Ri^2) and Dh = Do - Di, both of which an axis offset leaves alone
/// because both bores stay fully wetted; only K_geom moves.
double annulusLaminarConstant(const NetworkPipe& p)
{
    const double ro = 0.5 * p.diameter, ri = 0.5 * p.annulusInner;
    const double dh = p.diameter - p.annulusInner;
    if (p.annulusEccK > 0.0)
    {
        return 2.0 * dh * dh * M_PI * (ro * ro - ri * ri)/p.annulusEccK;
    }
    const double g = ro * ro * ro * ro - ri * ri * ri * ri
                   - (ro * ro - ri * ri) * (ro * ro - ri * ri)
                   / std::log(ro / ri);
    return 16.0 * dh * dh * (ro * ro - ri * ri) / g;
}

/// The power-law laminar constant of the annulus, C_pl(n, kappa) = f . Re_GP,
/// built so the blend's laminar asymptote IS the FB profile: demanding
/// f_lam . (L/Dh) . rho v^2 / 2 == 2 K L (Q/(2 pi I))^n with Q = v A and
/// Re_GP = rho v^(2-n) Dh^n / (K 8^(n-1)) gives
///     C_pl = 4 Dh^(n+1) A^n / ((2 pi I)^n . 8^(n-1)),
/// in which K cancels -- f . Re_GP is a geometry + law constant, as it must
/// be. At n = 1, I = G/8 collapses it onto annulusLaminarConstant exactly.
double annulusPowerLawLaminarConstant(const NetworkPipe& p)
{
    const double n = p.flowIndex;
    const double dh = p.diameter - p.annulusInner;
    const double a = pipeArea(p);
    const double twoPiI = 2.0 * M_PI * p.annulusDimQ;
    return 4.0 * std::pow(dh, n + 1.0) * std::pow(a, n)
         / (std::pow(twoPiI, n) * std::pow(8.0, n - 1.0));
}

/// Churchill (1977) with the laminar member swapped: identical parentheses,
/// but the (8/Re)^12 term carries C/8 instead of 8 so the Re -> 0 asymptote
/// is C/Re instead of 64/Re. C = 64 reproduces churchillFriction exactly
/// (64/8 == 8 in IEEE), which is what makes the circular branch bit-still.
double churchillAnnulusFriction(double reynolds, double relativeRoughness,
                                double laminarC)
{
    const double re = std::max(reynolds, 1e-12);
    const double a1 = std::pow(7.0 / re, 0.9) + 0.27 * relativeRoughness;
    if (a1 <= 0) cfFatal("churchillAnnulusFriction: invalid Re/roughness");
    const double a = std::pow(2.457 * std::log(1.0 / a1), 16.0);
    const double b = std::pow(37530.0 / re, 16.0);
    return 8.0 * std::pow(std::pow(laminarC / 8.0 / re, 12.0)
                              + std::pow(a + b, -1.5), 1.0 / 12.0);
}

/// The two profile integrals of the Fredrickson-Bird (1958) power-law
/// annulus, on the dimensionless variable s = r/Ro in (kappa, 1], with the
/// zero-shear radius lam = lambda/Ro. tau(r) = (dp/(2L))(r - lambda^2/r) is
/// the only stress field that balances the pressure force and vanishes at
/// one radius; matching the profile integrals from the two walls,
///     phi(lam) = int_kappa^lam (lam^2/s - s)^(1/n) ds
///              - int_lam^1 (s - lam^2/s)^(1/n) ds = 0,
/// has exactly one root (phi is increasing: the left side loses, the right
/// gains), and the dimensionless rate follows by swapping integration order,
///     I = 1/2 int_kappa^lam (lam^2/s - s)^(1/n) (lam^2 - s^2) ds
///       + 1/2 int_lam^1 (s - lam^2/s)^(1/n) (s^2 - lam^2) ds,
/// so Q = (dp/(2 K L))^(1/n) . 2 pi . I with I in m^(1+1/n) of the actual
/// bore. At n = 1 the pair reduces to the Newtonian annulus (I = G/8 with
/// the G above, the closed form annulusLaminarConductance uses).
void annulusFbIntegrals(const NetworkPipe& p, int strips, double& lambda,
                        double& dimQ)
{
    const double ro = 0.5 * p.diameter, ri = 0.5 * p.annulusInner;
    const double n = p.flowIndex;
    auto phi = [&](double lam)
    {
        double left = 0.0, right = 0.0;
        const double h1 = (lam - ri) / strips, h2 = (ro - lam) / strips;
        for (int i = 0; i <= strips; ++i)
        {
            const double s1 = ri + i * h1;
            const double w1 = (i == 0 || i == strips) ? 1.0
                            : (i % 2 ? 4.0 : 2.0);
            left += w1 * std::pow(std::max(0.0, lam * lam / s1 - s1), 1.0 / n);
        }
        for (int i = 0; i <= strips; ++i)
        {
            const double s2 = lam + i * h2;
            const double w2 = (i == 0 || i == strips) ? 1.0
                            : (i % 2 ? 4.0 : 2.0);
            right += w2 * std::pow(std::max(0.0, s2 - lam * lam / s2), 1.0 / n);
        }
        return left * h1 / 3.0 - right * h2 / 3.0;
    };
    double lo = ri * (1.0 + 1e-12), hi = ro * (1.0 - 1e-12);
    for (int it = 0; it < 200; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        if (phi(mid) < 0.0) lo = mid;
        else hi = mid;
        if (std::abs(hi - lo) <= 2.0 * std::numeric_limits<double>::epsilon()
                                      * hi)
            break;
    }
    lambda = 0.5 * (lo + hi);

    auto dimLess = [&](double lam)
    {
        double inner = 0.0, outer = 0.0;
        const double h1 = (lam - ri) / strips, h2 = (ro - lam) / strips;
        for (int i = 0; i <= strips; ++i)
        {
            const double s1 = ri + i * h1;
            const double w1 = (i == 0 || i == strips) ? 1.0
                            : (i % 2 ? 4.0 : 2.0);
            inner += w1
                   * std::pow(std::max(0.0, lam * lam / s1 - s1), 1.0 / n)
                   * (lam * lam - s1 * s1);
        }
        for (int i = 0; i <= strips; ++i)
        {
            const double s2 = lam + i * h2;
            const double w2 = (i == 0 || i == strips) ? 1.0
                            : (i % 2 ? 4.0 : 2.0);
            outer += w2
                   * std::pow(std::max(0.0, s2 - lam * lam / s2), 1.0 / n)
                   * (s2 * s2 - lam * lam);
        }
        return 0.5 * (inner * h1 / 3.0 + outer * h2 / 3.0);
    };
    // I is dimensional: the profile carries (dp/(2 K L))^(1/n) times
    // m^(1+1/n), so the integral over the actual bore is returned as-is and
    // Q = (dp/(2 K L))^(1/n) . 2 pi . dimQ.
    dimQ = dimLess(lambda);
}

/// The rectangular-duct laminar correction bracket (5.80), from the exact
/// Fourier-series solution of grad^2 u = -G/mu on a W x H rectangle with
/// no-slip walls: a parabolic particular solution across the SHORT side plus
/// a cosh/sin homogeneous series that cancels it on the two long walls. With
/// aspect = short/long <= 1 the integrated rate carries
///     bracket = 1 - (192 aspect/pi^5) sum_{n odd} tanh(n pi/(2 aspect))/n^5,
/// so Q = (dp/L) . long . short^3 . bracket/(12 mu) -- the plane-Poiseuille
/// gap^3 law times the side-wall correction. bracket -> 1 as aspect -> 0
/// (parallel plates, f.Re = 96) and = 0.42174 for a square (f.Re = 56.91).
/// Summed once by addPipe (it depends only on the aspect), so the rate solve
/// never re-sums it. 1/n^5 decay: 4001 odd terms leave a tail below 1e-15
/// even at aspect = 1, where tanh saturates slowest.
double rectLaminarBracket(double aspect)
{
    double sum = 0.0;
    for (int n = 1; n <= 8001; n += 2)
    {
        const double nn = static_cast<double>(n);
        sum += std::tanh(nn * M_PI / (2.0 * aspect))
             / (nn * nn * nn * nn * nn);
    }
    return 1.0 - (192.0 * aspect / std::pow(M_PI, 5)) * sum;
}

/// Newtonian laminar conductance Q/dp of the rectangular duct: the plane
/// Poiseuille gap^3 law on the SHORT side times the long side and the bracket.
double rectLaminarConductance(const NetworkPipe& p,
                              const FluidProperties& fluid)
{
    const double longSide = std::max(p.diameter, p.ductHeight);
    const double shortSide = std::min(p.diameter, p.ductHeight);
    return longSide * shortSide * shortSide * shortSide * p.rectBracket
         / (12.0 * fluid.mu * p.length);
}

/// The laminar constant C(aspect) = f . Re_Dh of the rectangular duct, from
/// the same series: f = 2 Dh dp/(L rho v^2), Re_Dh = rho v Dh/mu and v = Q/A
/// give f Re = 24 Dh^2/(short^2 bracket) = 96/((1+aspect)^2 bracket). 96 at
/// aspect -> 0 (parallel plates), 56.91 at aspect = 1 (square) -- so
/// Churchill's laminar line is the aspect -> 0 member of this family and the
/// all-regime blend only has to swap that member, exactly as for the annulus.
double rectLaminarConstant(const NetworkPipe& p)
{
    const double aspect = std::min(p.diameter, p.ductHeight)
                        / std::max(p.diameter, p.ductHeight);
    const double onePlus = 1.0 + aspect;
    return 96.0 / (onePlus * onePlus * p.rectBracket);
}

/// Hagen-Poiseuille slope pi D^4 / (128 mu L): the maximum conductance of a
/// pipe (laminar limit), used as the fallback when the difference quotient
/// degenerates at dp ~ 0.
double laminarConductance(const NetworkPipe& p, const FluidProperties& fluid)
{
    if (pipeIsRect(p)) return rectLaminarConductance(p, fluid);
    if (pipeIsAnnulus(p)) return annulusLaminarConductance(p, fluid);
    return M_PI * std::pow(p.diameter, 4) / (128.0 * fluid.mu * p.length);
}

/// Piecewise-linear table lookup clamped at both ends (NaN/p <= 0 takes the
/// first entry, so early Newton iterates can never produce a bad property).
double pvtLookup(const std::vector<double>& grid,
                 const std::vector<double>& vals, double p)
{
    const std::size_t n = grid.size();
    if (n == 0) return 0.0;
    if (!(p > grid[0])) return vals[0];
    if (p >= grid[n - 1]) return vals[n - 1];
    std::size_t lo = 0, hi = n - 1;
    while (hi - lo > 1)
    {
        const std::size_t mid = (lo + hi) / 2;
        if (grid[mid] <= p) lo = mid;
        else hi = mid;
    }
    const double t = (p - grid[lo]) / (grid[hi] - grid[lo]);
    return vals[lo] + (vals[hi] - vals[lo]) * t;
}

struct RawFlow
{
    double velocity = 0; // m/s, signed from -> to
    double reynolds = 0;
    double friction = 0;
};

/// Pump head at rate q (m of fluid): H(q) = h0 - (h0 - h1)(q/q1)^2 clamped at
/// 0; constant h0 when q1 = 0. Even in q: the machine adds the same head
/// magnitude whichever way the flow runs (never sucks).
/// 5.81 (D1): N stages in series multiply the head at the same flow, and a
/// relative speed s rescales it by the Euler affinity laws (Q ~ s, H ~ s^2):
/// H_s(q) = s^2 H_1(q/s) = s^2 h0 - (h0-h1)(q/q1)^2 -- the quadratic term is
/// speed-invariant, only the shutoff head scales. Combined H = stages * H_s.
/// stages 1 / speed 1 is the law above bit for bit (x*1.0 == x in IEEE).
double pumpHeadOf(const NetworkPipe& p, double q)
{
    if (!(p.pumpHead > 0)) return 0.0;
    const double h0 = p.pumpShutoffHead > 0 ? p.pumpShutoffHead : p.pumpHead;
    const double stages = p.pumpStages;
    const double s2 = p.pumpSpeed * p.pumpSpeed;
    double h;
    if (!(p.pumpRate > 0)) h = stages * s2 * h0;
    else h = stages * (s2 * h0 - (h0 - p.pumpHead) * q * q / (p.pumpRate * p.pumpRate));
    // 5.86 (D1 increment 2): Turpin (1986) gas interference, read through the
    // shared pumpHeadOf so the degradation multiplies the whole curve above.
    // q_G/q_L = g/(1-g) is the free-gas/liquid ratio at intake; p_i is the
    // ABSOLUTE intake pressure in psia. H_tp = H_sp * exp(-a1 * (q_G/q_L)),
    // a1 = (346430/p_i^2)(q_G/q_L) - 410/p_i. The multiplier is skipped
    // entirely at g == 0, so the zero-declaration path is bit-identical to 5.81.
    if (p.pumpGasFraction > 0)
    {
        const double r = p.pumpGasFraction / (1.0 - p.pumpGasFraction);
        const double psi = p.pumpIntakePressure / kPascalPerPsi;  // psia
        const double a1 = (346430.0 / (psi * psi)) * r - 410.0 / psi;
        h *= std::exp(-a1 * r);
    }
    return std::max(h, 0.0);
}

/// 5.191 (G7): the head a positive-displacement machine adds at rate q.
/// The machine's law is q = n V - C dp, so at a liquid density rho (taken
/// where the existing centrifugal curve is read, i.e. at the pipe's upstream
/// state) its head is
///   H(q) = (nV - q) / (C rho g).
/// Linear and decreasing in q, so the rate bracket's monotone branch is the
/// same one the centrifugal family uses. A pipe without a PCP returns 0 and
/// every existing number is bit-identical.
double pcpHeadOf(const NetworkPipe& p, double q, double rho)
{
    if (!p.pcpOn) return 0.0;
    const double ideal = p.pcpSpeed * p.pcpDisplacement;
    return (ideal - q) / (p.pcpSlip * rho * kGravity);
}

/// 5.192 (G8): a jet pump's row law -- the "one row, two streams" slot.
/// The machine's rise is its own constant-area throat balance, written as a
/// head over the row's combined rate q = q_n + q_s:
///   H(q) = rho v_n^2 / m^2 * [ (m-1-r)^2/(m-1) - F (1+r)^2/2 ],
///   r = (q - q_n)/q_n, m = A_t/A_n, v_n = q_n/A_n, F = f L_t/D_t.
/// This is the algebra jetPumpSolve prices in its imposed-rate mode, not a
/// curve fit: r is what the row's own rate leaves for entrainment once the
/// declared motive stream has taken q_n. It is non-increasing in q (the
/// bracket term falls, the friction term grows), bounded, and clamped at the
/// element's two own states: r <= 0 is NO-SUCTION (the stall head, motive
/// stream alone) and H <= 0 is FREE-JET (no rise at the entrainment ceiling).
/// So the rate bracket's monotone branch is the family's, unchanged.
struct JetRow
{
    double head = 0;   // m of the combined stream
    double ratio = 0;  // r, clamped to the element's envelope
    int state = -1;    // -1 undeclared, 0 PUMPING, 1 NO-SUCTION, 2 FREE-JET
};

JetRow jetRowOf(const NetworkPipe& p, double q, double rho)
{
    JetRow r;
    if (!p.jetOn) return r;
    const double m = p.jetThroatArea / p.jetNozzleArea;
    const double vn = p.jetMotiveRate / p.jetNozzleArea;
    const double dt = std::sqrt(4.0 * p.jetThroatArea / M_PI);
    const double fLine = p.jetFrictionFactor * p.jetThroatLength / dt;

    double ratio = (q - p.jetMotiveRate) / p.jetMotiveRate;
    if (ratio < 0.0) ratio = 0.0;
    if (ratio > m - 1.0) ratio = m - 1.0;
    const double s = m - 1.0 - ratio;
    // The bracket is a pressure coefficient: rho v_n^2/m^2 times it is the
    // Pa the throat hands over. It is divided by rho g here so the row law
    // speaks the same unit as pumpHeadOf, and the readout dpJet can charge
    // -rho g H without pricing the density twice.
    double head = vn * vn / (m * m * kGravity)
                * (s * s / (m - 1.0)
                   - 0.5 * fLine * (1.0 + ratio) * (1.0 + ratio));
    r.ratio = ratio;
    r.head = head;
    r.state = head <= 0.0 ? 2 : (ratio <= 0.0 ? 1 : 0);
    return r;
}

double jetHeadOf(const NetworkPipe& p, double q, double rho)
{
    return jetRowOf(p, q, rho).head;
}

/// 5.191/5.192: every machine a row can carry, in one place -- the bracket,
/// the dead branch and the readouts all read the same sum.
double machineHeadOf(const NetworkPipe& p, double q, double rho)
{
    return pumpHeadOf(p, q) + pcpHeadOf(p, q, rho) + jetHeadOf(p, q, rho);
}

/// P10 (5.156): the pressure ratio a pipe's compressor reads at rate q. The
/// curve is the same quadratic family as pumpHeadOf, read in ratio instead
/// of metres of head. A pipe without a compressor returns 1, so every
/// existing number is bit-identical to the path before this increment.
double compressorRatioOf(const NetworkPipe& p, double q)
{
    if (!(p.compressorRatio > 0.0)) return 1.0;
    const double shutoff = p.compressorShutoffRatio > 0.0
                             ? p.compressorShutoffRatio
                             : p.compressorRatio;
    const double x = q / p.compressorRatedRate;
    return shutoff - (shutoff - p.compressorRatio) * x * x;
}

/// The Reynolds number a rheology is reported with.
///
/// Newtonian (nothing declared): rho u d / mu. Yielded with an index: the
/// apparent number 8 rho u^2 / tau_w, because no material viscosity enters the
/// law (see yieldedWallStress). Otherwise the Metzner-Reed form
///     Re_MR = rho u^(2-n) d^n / (K' 8^(n-1))
/// with the consistency corrected by the wall shear rate rather than by the
/// apparent one, K' = K ((3n+1)/(4n))^n. That exponent is the whole point:
/// Rabinowitsch's wall rate is gamma_w = (3n+1)/(4n) . 8u/d, so the stress the
/// profile actually holds is tau_w = K gamma_w^n = K' (8u/d)^n, and only with
/// that K' does Churchill's laminar line f = 64/Re read back as 8 tau_w /
/// (rho u^2), i.e. as the integrated velocity profile rather than as a
/// correlation of it. The textbook constant K (3n+1)/(4n) - the exponent-1
/// version - collapses turbulent data onto the Newtonian Moody line but sits
/// ((3n+1)/(4n))^(1-1/n) off the laminar solution, which at n = 0.5 is 25% in
/// the rate; 5.47 measured that gap from the profile integral and this line is
/// the correction. At n = 1 the factor is 1/1 = 1 exactly, so a
/// declared-but-Newtonian pipe is still bit-identical to an undeclared one.
///
/// A yielded pipe with an index (Herschel-Bulkley) has no material viscosity at
/// all - mu does not enter its law - so the only Reynolds number that can be
/// reported is the apparent one: mu_app = tau_w d / (8 u), read off the wall
/// stress the fluid actually holds at this rate. With that number and the exact
/// f = 8 tau_w / (rho u^2) the product f * Re is 64 by definition, so the
/// assertion that leans on is the closed form for tau_w, not this identity.
double yieldedWallStress(const NetworkPipe& pipe, const FluidProperties& fluid,
                         double q);

double generalizedReynolds(const NetworkPipe& pipe, const FluidProperties& fluid,
                           double speed)
{
    if (pipe.flowIndex > 0 && pipe.yieldStress > 0)
    {
        const double tw = yieldedWallStress(pipe, fluid,
                                           speed * pipeArea(pipe));
        return 8.0 * fluid.rho * speed * speed / tw;
    }
    // Both branches read the hydraulic diameter: for a circle it IS the
    // diameter, so the circular numbers never move. A power-law annulus keeps
    // Dodge-Metzner's own form on Dh -- NO Rabinowitsch factor: the (3n+1)/4n
    // wall rate is a circular-geometry statement, and the annulus laminar
    // anchor is the FB profile itself (5.73), not a 64/Re line that would
    // need re-priceing.
    const double dh = pipeHydraulicDiameter(pipe);
    if (!(pipe.flowIndex > 0))
        return fluid.rho * speed * dh / fluid.mu;
    const double n = pipe.flowIndex;
    const double k = pipeIsAnnulus(pipe)
                   ? pipe.consistency
                   : pipe.consistency
                     * std::pow((3.0 * n + 1.0) / (4.0 * n), n);
    return fluid.rho * std::pow(speed, 2.0 - n) * std::pow(dh, n)
         / (k * std::pow(8.0, n - 1.0));
}

/// Buckingham (1921) exact laminar solution for a Bingham plastic, inverted
/// for the wall stress. With tw = dp d / (4 L) and xi = tau0 / tw the rate is
///     Q = (pi d^3 / (32 mu)) tw [1 - (4/3) xi + (1/3) xi^4],
/// which rearranges (multiply the bracket by tw) into one monotone equation
///     F(tw) = tw - (4/3) tau0 + tau0^4 / (3 tw^3) = B,  B = 32 mu Q / (pi d^3)
/// with F(tau0) = 0 and F' = 1 - (tau0/tw)^4 > 0 for tw > tau0, so the root is
/// unique and the bracket [tau0, B + (4/3) tau0] is analytic -- the right end
/// already gives B + tau0^4/(3 tw^3) >= B. Newton with a bisection fallback
/// then lands on the root to rounding, which is what the closed-form checks
/// need: a 1e-9 slack here would be indistinguishable from a wrong law.
/// tw = tau0 (Q = 0) and tau0 -> 0 both come out right: the first is the left
/// end of the bracket, the second makes F(tw) = tw identically.
double binghamWallStress(const NetworkPipe& pipe, const FluidProperties& fluid,
                         double q)
{
    const double tau0 = pipe.yieldStress;
    const double b = 32.0 * fluid.mu * std::abs(q)
                   / (M_PI * pipe.diameter * pipe.diameter * pipe.diameter);
    if (!(b > 0.0)) return tau0; // no rate: the wall holds exactly the threshold

    double lo = tau0, hi = b + (4.0 / 3.0) * tau0;
    double tw = 0.5 * (lo + hi);
    const double tau04 = tau0 * tau0 * tau0 * tau0;
    // The exit below is machine precision (bracket width <= 2 ulp), so the cap
    // is a safety net, not a gate -- but it is the same pipe-solve cap policy.
    for (int it = 0; it < gates::steadyPipeMaxIter; ++it)
    {
        const double t2 = tw * tw;
        const double t4 = t2 * t2;
        const double f = tw - (4.0 / 3.0) * tau0
                       + tau04 / (3.0 * t2 * tw) - b;
        if (f > 0.0) hi = tw;
        else lo = tw;
        // F' = 1 - tau0^4/tw^4, which vanishes at the left end: bisect there
        const double df = 1.0 - tau04 / t4;
        const double next = df > 0.0 ? tw - f / df : 0.5 * (lo + hi);
        tw = (next > lo && next < hi) ? next : 0.5 * (lo + hi);
        if (!(hi > lo)) break;
        if (std::abs(hi - lo) <= 2.0 * std::numeric_limits<double>::epsilon() * hi)
            break;
    }
    return tw;
}

/// Herschel-Bulkley (tau = tau0 + K gdot^n) exact laminar rate in a circular
/// pipe, as a function of the wall stress. With tw = dp d / (4 L) and
/// D = tw - tau0 the profile integral is termwise integrable,
///     Q = pi R^3 tw^-3 K^(-1/n) n D^(1+1/n)
///         [ tau0^2/(n+1) + 2 tau0 D/(2n+1) + D^2/(3n+1) ],
/// which is not a correlation: at n = 1, K = mu the bracket times D^2 is
/// Buckingham's [1 - (4/3) xi + (1/3) xi^4] times tw^4/4, and at tau0 = 0 only
/// the last term survives and the expression is the power-law closed form of
/// 5.45. So one formula contains both laws already shipped, and the threshold
/// stays geometric: Q -> 0 as tw -> tau0, i.e. dp_static = 4 tau0 L / d for
/// every n, while the *shape* of the landing near it depends on n (the leading
/// power of D is 1 + 1/n, which is 2 for Bingham and 3 at n = 0.5).
double herschelBulkleyRate(double tw, const NetworkPipe& pipe)
{
    const double n = pipe.flowIndex;
    const double d = tw - pipe.yieldStress;
    if (!(d > 0.0) || !(tw > 0.0)) return 0.0;
    const double r = 0.5 * pipe.diameter;
    return M_PI * r * r * r * std::pow(tw, -3.0)
         * std::pow(pipe.consistency, -1.0 / n) * n
         * std::pow(d, 1.0 + 1.0 / n)
         * (pipe.yieldStress * pipe.yieldStress / (n + 1.0)
          + 2.0 * pipe.yieldStress * d / (2.0 * n + 1.0)
          + d * d / (3.0 * n + 1.0));
}

/// The same law inverted: the wall stress a demanded rate must be held at.
/// Monotone in tw (the yield only removes fluid from the slow end of the
/// profile, so a higher wall stress always carries more), so a bracket plus
/// bisection is enough - no derivative, and no Newton stall at the left end
/// where the slope vanishes. The left end is the threshold exactly; the right
/// end starts at the threshold plus the tau0 = 0 (power-law) stress for the
/// same rate - a lower bound on the root, because yielding only costs rate -
/// and doubles until the law clears the demand, which it always does since the
/// rate grows like tw^(1/n).
double herschelBulkleyWallStress(const NetworkPipe& pipe, double q)
{
    const double tau0 = pipe.yieldStress;
    const double demand = std::abs(q);
    if (!(demand > 0.0)) return tau0; // no rate: the wall holds the threshold

    const double n = pipe.flowIndex;
    const double r = 0.5 * pipe.diameter;
    const double twPl =
        pipe.consistency
      * std::pow(demand * (3.0 * n + 1.0) / (M_PI * r * r * r * n), n);
    double lo = tau0, hi = tau0 + twPl;
    for (int expand = 0; expand < 200 && !(herschelBulkleyRate(hi, pipe) > demand);
         ++expand) hi *= 2.0;

    for (int it = 0; it < 200; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        if (herschelBulkleyRate(mid, pipe) > demand) hi = mid;
        else lo = mid;
        if (!(hi > lo)) break;
        if (std::abs(hi - lo) <= 2.0 * std::numeric_limits<double>::epsilon() * hi)
            break;
    }
    return 0.5 * (lo + hi);
}

/// Which yielded law a pipe's wall stress comes from: Buckingham when only a
/// threshold and the fluid's own mu are declared (the plastic number keeps the
/// Reynolds convention of 5.46), Herschel-Bulkley when an index is declared
/// with the threshold.
double yieldedWallStress(const NetworkPipe& pipe, const FluidProperties& fluid,
                         double q)
{
    if (pipe.flowIndex > 0) return herschelBulkleyWallStress(pipe, q);
    return binghamWallStress(pipe, fluid, q);
}

/// Wilson & Thomas (2006) Eq (17)-(19): the plastic Reynolds number at the
/// laminar-turbulent locus of a Bingham plastic, a function of He alone. Shared
/// by the locus report and the friction switch so the two cannot drift apart.
double binghamTransitionRe(double hedstrom)
{
    if (hedstrom <= 1700.0) return 2100.0;
    if (hedstrom <= 1.0e5)    return 80.0 * std::pow(hedstrom, 0.40);
    return 25.0 * std::sqrt(hedstrom);
}

/// Wilson & Thomas (2006) Eq (15)/(16), cross-checked verbatim against Thomas &
/// Wilson (2007) Eq [3]/[4]: the smooth-wall turbulent Darcy friction factor of
/// a Bingham plastic. With beta = tau_w/tau_B, He = rho tau_B d^2/mu_B^2 and
/// Re_B = rho V d/mu_B, the mean-to-shear velocity ratio V/U* = sqrt(8/f) obeys
///   sqrt(8/f) = 2.5 ln(sqrt He) + 2.5 ln[sqrt(beta)(beta-1)/(beta+1)]
///               + 11.6/beta - chi(beta),
///   chi(beta) = -2.5 ln[(beta-1)/beta] - 2.5 (beta+0.5)/beta^2.
/// beta = f Re_B^2/(8 He) and sqrt(8/f) = Re_B/sqrt(He beta), so this is one
/// monotone equation in beta > 1, solved by bracketed bisection to rounding.
/// As tau_B -> 0 (He -> 0, beta -> inf at fixed He*beta = Re*^2) the log terms
/// collapse to 2.5 ln(Re*) while chi and 11.6/beta vanish: the Newtonian
/// smooth-wall log law, not a separate statement. log1p keeps the (beta-1)/beta
/// and (beta-1)/(beta+1) ratios accurate when beta is large.
double yieldedTurbulentFriction(double rePlastic, double hedstrom)
{
    auto chi = [](double b)
    {
        return -2.5 * std::log1p(-1.0 / b) - 2.5 * (b + 0.5) / (b * b);
    };
    auto resid = [&](double b)
    {
        const double logTerm = 1.25 * std::log(b)
                             + 2.5 * std::log1p(-2.0 / (b + 1.0));
        return 1.25 * std::log(hedstrom) + logTerm + 11.6 / b - chi(b)
             - rePlastic / std::sqrt(hedstrom * b);
    };
    double lo = 1.0 + 1.0e-12, hi = 2.0;
    for (int e = 0; e < 200 && resid(hi) < 0.0; ++e) hi *= 2.0;
    for (int it = 0; it < 200; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        if (resid(mid) > 0.0) hi = mid; else lo = mid;
        if (!(hi > lo)) break;
        if (hi - lo <= 2.0 * std::numeric_limits<double>::epsilon() * hi) break;
    }
    const double beta = 0.5 * (lo + hi);
    return 8.0 * hedstrom * beta / (rePlastic * rePlastic);
}

/// Friction + valve loss (Pa, positive) at an imposed rate magnitude q: the
/// wall-law drop from pipeFrictionPa -- which the transient row also calls --
/// plus the local loss. A yielded pipe answers with the wall-stress solution,
/// whose q = 0 limit is the static threshold -- that is the whole difference,
/// and why the rate law needs an explicit dead branch (see equipmentFlow).
double pipeLossAt(const NetworkPipe& pipe, const FluidProperties& fluid,
                  double q)
{
    // A coupled element's loss is its partner's statement, not ours: the
    // detailed model already reports the whole drop of that region, so adding
    // this pipe's own friction on top would charge the network twice for the
    // same wall. Linear in q because the coupling layer refreshes the
    // coefficient every outer iteration (exact for laminar, secant-like for
    // anything else).
    if (pipe.coupledResistance > 0) return pipe.coupledResistance * q;

    const double v = q / pipeArea(pipe);

    // One call, both laws, the same parentheses the transient row uses; the
    // valve term is the only thing the steady kernel adds on top.
    return pipeFrictionPa(pipe, fluid, q)
         + pipe.valveK * fluid.rho * v * v / 2.0;
}

/// Friction factor to report for a solved rate: the Churchill number for the
/// all-regime law, and for a yielded pipe the exact factor implied by that
/// law's drop, f = 8 tau_w / (rho v^2) -- the number that reproduces the drop
/// printed next to it, not a correlation pasted onto a yielded wall. For a
/// Bingham plastic the Reynolds number printed beside it is the plastic one
/// (rho v d / mu), so the pair reads f Re g(xi) = 64; for Herschel-Bulkley
/// there is no material viscosity and the printed number is the apparent one,
/// which makes f Re = 64 an identity between the two reports.
/// At v = 0 there is no rate to normalise against, so the report is 0.
double pipeFrictionReport(const NetworkPipe& pipe, const FluidProperties& fluid,
                          double reynolds, double speed)
{
    if (pipe.yieldStress > 0)
    {
        if (!(speed > 0)) return 0.0;
        // 5.83: past the locus a smooth Bingham pipe reports the Wilson-Thomas
        // turbulent factor -- the same number pipeFrictionPa charges -- so the
        // printed f still reproduces the printed drop. Below it (and for HB or
        // a rough wall) the exact laminar factor 8 tau_w / (rho v^2).
        if (!(pipe.flowIndex > 0) && pipe.roughness <= 0.0)
        {
            const double he = pipe.yieldStress * pipe.diameter * pipe.diameter
                            * fluid.rho / (fluid.mu * fluid.mu);
            const double reb = fluid.rho * speed * pipe.diameter / fluid.mu;
            if (reb > binghamTransitionRe(he))
                return yieldedTurbulentFriction(reb, he);
        }
        return 8.0 * yieldedWallStress(pipe, fluid, speed * pipeArea(pipe))
             / (fluid.rho * speed * speed);
    }
    if (pipeIsRect(pipe))
    {
        const double dh = pipeHydraulicDiameter(pipe);
        return churchillAnnulusFriction(reynolds, pipe.roughness / dh,
                                        rectLaminarConstant(pipe));
    }
    if (pipeIsAnnulus(pipe))
    {
        const double dh = pipeHydraulicDiameter(pipe);
        const double cLam = pipe.flowIndex > 0
            ? annulusPowerLawLaminarConstant(pipe)
            : annulusLaminarConstant(pipe);
        return churchillAnnulusFriction(reynolds, pipe.roughness / dh, cLam);
    }
    return churchillFriction(reynolds, pipe.roughness / pipe.diameter);
}

/// Pipe + equipment law Q(dpTotal): dpTotal = rho g dz + loss(q) - rho g H(q).
/// Solved as a bracketed Illinois secant on the rate. The pump-head fixed
/// point (hydraulics need the head, the curve needs the rate) is NOT
/// contractive whenever the curve falls faster than the system rises
/// (rho g |H'| > d(loss)/dq): it oscillates and latches onto the H = 0 clamp.
/// Root-solving is monotone on the branch picked by the sign of
/// G(0) = rho g dz - rho g h0 - dpTotal: loss grows without bound in |q|
/// while -rho g H(q) stays bounded, so exactly one crossing per branch.
RawFlow equipmentFlow(const NetworkPipe& pipe, const FluidProperties& fluid,
                      double dpTotal, double& Hq)
{
    const bool forward = fluid.rho * kGravity * pipe.elevation
                       - fluid.rho * kGravity
                             * machineHeadOf(pipe, 0.0, fluid.rho)
                       - dpTotal <= 0.0;
    auto g = [&](double x) // residual along the chosen branch, |q| = x >= 0
    {
        const double s = forward ? 1.0 : -1.0;
        return s * pipeLossAt(pipe, fluid, x)
             + fluid.rho * kGravity * pipe.elevation
             - fluid.rho * kGravity * machineHeadOf(pipe, x, fluid.rho)
            - dpTotal;
    };

    // bracket: forward g(0) <= 0, expand to a positive sample; reverse
    // g(0) > 0, expand to a negative sample (keep the near-zero endpoint)
    const double pScale = std::max(std::abs(dpTotal),
                                   fluid.rho * kGravity
                                     * (std::abs(pipe.elevation)
                                        + machineHeadOf(pipe, 0.0, fluid.rho)));
    const double xSeed = std::max(1e-3 * pipeArea(pipe), 1e-8);
    double xLo = 0.0, fLo = g(0.0);

    // The dead branch of a yielded pipe. loss(0) is not 0 there, it is the
    // static threshold 4 tau0 L / d, so g(0) can already sit on the far side of
    // the crossing: the asked drop never reaches the wall's threshold and the
    // exact answer is no flow at all. Falling through would expand the bracket
    // 60 times on two same-signed samples and hand the secant a garbage pair,
    // so this has to be decided before the bracket, not damped inside it.
    // Newtonian and power-law pipes keep the plateau branch untouched -- their
    // g(0) is the sign that sent them here in the first place.
    if (pipe.yieldStress > 0 && (forward ? fLo >= 0.0 : fLo <= 0.0))
    {
        Hq = machineHeadOf(pipe, 0.0, fluid.rho);
        return RawFlow();
    }

    double xHi = xSeed, fHi = g(xHi);
    int expansions = 0;
    while (!(forward ? fHi > 0.0 : fHi < 0.0))
    {
        if (forward ? fHi <= 0.0 : fHi >= 0.0)
        {
            xLo = xHi;
            fLo = fHi;
        }
        xHi *= 4.0;
        fHi = g(xHi);
        if (++expansions > 60) break; // loss dominates long before this
    }

    double x = xLo + (xHi - xLo) * fLo / (fLo - fHi);
    for (int it = 0; it < gates::steadyPipeMaxIter; ++it)
    {
        const double fx = g(x);
        if (std::abs(fx) < gates::steadyPipeResidualGate * std::max(pScale, 1.0)) break;
        if (forward ? fx > 0.0 : fx < 0.0)
        {
            xHi = x;
            fHi = fx;
            // Illinois: bias the retained endpoint against re-selection
            fHi *= 0.5;
        }
        else
        {
            xLo = x;
            fLo = fx;
            fLo *= 0.5;
        }
        const double xNew = xLo + (xHi - xLo) * fLo / (fLo - fHi);
        if (!(xNew > xLo && xNew < xHi)) x = 0.5 * (xLo + xHi);
        else x = xNew;
        if (std::abs(xHi - xLo) < gates::steadyPipeBracketGate * std::max(xHi, 1e-30)) break;
    }

    RawFlow r;
    r.velocity = (forward ? x : -x) / pipeArea(pipe);
    r.reynolds = generalizedReynolds(pipe, fluid, std::abs(r.velocity));
    r.friction = pipeFrictionReport(pipe, fluid, r.reynolds,
                                    std::abs(r.velocity));
    const double qAbs = std::abs(r.velocity) * pipeArea(pipe);
    Hq = machineHeadOf(pipe, qAbs, fluid.rho);
    return r;
}

/// Dense symmetric-ish n x n solve: see DenseSolve.h (shared with the
/// transient kernel).

} // namespace

/// 5.193 (G9): the survey's own vertical drop. theta is linear in measured
/// depth between stations, so each interval integrates in closed form:
///   dz_i = integral ( cos(th1) - cos(th2) ) / (th2 - th1)   (th2 != th1)
///        = sin(th1) * ds                                     (th2 == th1)
/// with th in radians and ds = md[i+1] - md[i]. No numerical quadrature and
/// no sampling error: the number the pipe is then charged for is the integral
/// the survey states, not a grid's guess at it.
double surveyElevationOf(const NetworkPipe& p)
{
    if (!p.surveyOn) return p.elevation;
    double dz = 0.0;
    for (std::size_t i = 0; i + 1 < p.surveyMd.size(); ++i)
    {
        const double ds = p.surveyMd[i + 1] - p.surveyMd[i];
        const double th1 = p.surveyInclDeg[i] * M_PI / 180.0;
        const double th2 = p.surveyInclDeg[i + 1] * M_PI / 180.0;
        dz += (std::abs(th2 - th1) < 1e-12)
                  ? std::sin(th1) * ds
                  : (std::cos(th1) - std::cos(th2)) / (th2 - th1) * ds;
    }
    return dz;
}

/// 5.193 (G9): sin(theta) at each march segment's midpoint. The survey is
/// piecewise-linear in md, so a midpoint is the honest sample of an arc that
/// bends inside the segment (and the average of sin over a linear-theta
/// interval, by the same closed form above, is within a per-mille of it).
std::vector<double> surveySegSinOf(const NetworkPipe& p, int nSeg)
{
    if (!p.surveyOn)
        return std::vector<double>(static_cast<std::size_t>(nSeg),
                                   p.elevation / std::max(p.length, 1e-300));
    std::vector<double> s(static_cast<std::size_t>(nSeg));
    const double ds = p.length / static_cast<double>(nSeg);
    std::size_t k = 0;
    for (int i = 0; i < nSeg; ++i)
    {
        const double mid = (i + 0.5) * ds;
        while (k + 2 < p.surveyMd.size() && mid >= p.surveyMd[k + 1]) ++k;
        const double m0 = p.surveyMd[k], m1 = p.surveyMd[k + 1];
        const double f = (m1 > m0) ? (mid - m0) / (m1 - m0) : 0.0;
        const double th = (p.surveyInclDeg[k]
                           + f * (p.surveyInclDeg[k + 1]
                                  - p.surveyInclDeg[k])) * M_PI / 180.0;
        s[static_cast<std::size_t>(i)] = std::sin(th);
    }
    return s;
}


/// The friction drop (Pa, positive) a wall law resists at an imposed rate
/// magnitude: Churchill on the Metzner-Reed number when no threshold is
/// declared -- so Newtonian and power-law share one line, one function, and one
/// set of parentheses -- and the yielded solution's own wall stress when one
/// is. Both kernels call this, and it is the only place the wall law is written
/// outside the steady loss: the transient momentum row signs it and keeps
/// gravity and inertia on its own side, while the steady loss adds its valve
/// term on top. So a `flowIndex` declared in a transient case cannot be read by
/// a different algebra than the same declaration in a steady one -- before this
/// line existed, 5.45 to 5.47 were the steady statements of the law while the
/// transient row computed `rho u d / mu` by hand, which made a declared
/// rheology there not refused but *ignored*: the worst kind of silent number.
/// At namespace scope (not the anonymous one above) because the transient
/// kernel is a different translation unit.
double pipeFrictionPa(const NetworkPipe& pipe, const FluidProperties& fluid,
                      double q)
{
    if (pipe.yieldStress > 0)
    {
        // 5.83: a smooth Bingham plastic strictly past its published locus
        // (Wilson & Thomas 2006) is priced by the turbulent closure, not the
        // laminar Buckingham wall law. Strict > keeps dpAtTransition -- which
        // is evaluated AT Re_B,t -- on the laminar law. Herschel-Bulkley has no
        // plastic viscosity and so no Wilson-Thomas statement; a rough wall is
        // outside the smooth-wall closure. Both keep the laminar law until
        // their own increments land.
        if (!(pipe.flowIndex > 0) && pipe.roughness <= 0.0 && q != 0.0)
        {
            const double v = std::abs(q) / pipeArea(pipe);
            const double he = pipe.yieldStress * pipe.diameter * pipe.diameter
                            * fluid.rho / (fluid.mu * fluid.mu);
            const double reb = fluid.rho * v * pipe.diameter / fluid.mu;
            if (reb > binghamTransitionRe(he))
            {
                const double f = yieldedTurbulentFriction(reb, he);
                return (f * pipe.length / pipe.diameter) * fluid.rho * v * v / 2.0;
            }
        }
        return 4.0 * pipe.length * yieldedWallStress(pipe, fluid, q)
             / pipe.diameter;
    }
    const double v = q / pipeArea(pipe);
    const double re = generalizedReynolds(pipe, fluid, std::abs(v));
    if (pipeIsRect(pipe))
    {
        // 5.80: the same all-regime blend on Dh, laminar member swapped for
        // the rectangular duct's own Fourier-series constant (96 -> 56.91).
        const double dh = pipeHydraulicDiameter(pipe);
        const double f = churchillAnnulusFriction(re, pipe.roughness / dh,
                                                  rectLaminarConstant(pipe));
        return (f * pipe.length / dh) * fluid.rho * v * v / 2.0;
    }
    if (pipeIsAnnulus(pipe))
    {
        // 5.73: the same all-regime blend on Dh, with the laminar member
        // swapped for the geometry's own constant -- the FB closed form for a
        // Newtonian fluid, the FB profile integral for a power law. The
        // power-law laminar drop reads back as
        //     dp = 2 K L (Q / (2 pi I))^n
        // exactly, because the constant was BUILT from that profile (see the
        // C_pl derivation in the header of annulusFbIntegrals users).
        const double dh = pipeHydraulicDiameter(pipe);
        const double cLam = pipe.flowIndex > 0
            ? annulusPowerLawLaminarConstant(pipe)
            : annulusLaminarConstant(pipe);
        const double f = churchillAnnulusFriction(re, pipe.roughness / dh,
                                                  cLam);
        return (f * pipe.length / dh) * fluid.rho * v * v / 2.0;
    }
    const double f = churchillFriction(re, pipe.roughness / pipe.diameter);
    return (f * pipe.length / pipe.diameter) * fluid.rho * v * v / 2.0;
}

/// Where a yielded pipe sits against the published laminar-turbulent locus.
/// The locus numbers are a report, not a charge -- but they are the switch
/// pipeFrictionPa reads: a smooth Bingham pipe strictly past Re_B,t is priced
/// by the Wilson-Thomas turbulent closure (5.83), below it by the laminar wall
/// law. He and the plastic Reynolds number Re_B are built on the fluid's own
/// mu, which is exactly why a Herschel-Bulkley pipe gets no report (its law
/// contains no plastic viscosity at all). `dpAtTransition` is the pipe's own
/// laminar drop at the locus rate -- the last laminar point. The three band
/// constants and their He ranges are as printed in Wilson & Thomas (2006);
/// their Eq (17), the fit that bridges 2100 down into the low-He range, could
/// not be transcribed from the source layout, so band 0 is the flat 2100 the
/// paper itself calls the engineering value.
YieldedTransition yieldedTransitionOf(const NetworkPipe& pipe,
                                      const FluidProperties& fluid, double q)
{
    YieldedTransition t;
    if (!(pipe.yieldStress > 0) || pipe.flowIndex > 0) return t;

    const double area = 0.25 * M_PI * pipe.diameter * pipe.diameter;
    t.hedstrom = pipe.yieldStress * pipe.diameter * pipe.diameter * fluid.rho
               / (fluid.mu * fluid.mu);
    t.rePlastic = fluid.rho * (std::abs(q) / area) * pipe.diameter / fluid.mu;

    t.reTransition = binghamTransitionRe(t.hedstrom);
    t.band = t.hedstrom <= 1700.0 ? 0 : (t.hedstrom <= 1.0e5 ? 1 : 2);

    // the rate that number is a velocity of, and the laminar drop this pipe's
    // own wall law asks for there. dpAtTransition is the LAST laminar point --
    // the threshold the injection end crosses into turbulence -- so it is read
    // from the shared laminar wall stress directly, never the turbulent switch
    // (a round-trip through rateAtTransition could otherwise trip it at Re_B,t).
    t.rateAtTransition = t.reTransition * fluid.mu * area
                       / (fluid.rho * pipe.diameter);
    t.dpAtTransition = 4.0 * pipe.length
                     * yieldedWallStress(pipe, fluid, t.rateAtTransition)
                     / pipe.diameter;
    t.relTransitionSpeed = (t.rateAtTransition / area)
                         * std::sqrt(fluid.rho / pipe.yieldStress);
    t.pastLaminarLimit = t.rePlastic > t.reTransition;
    return t;
}

label OneDNetwork::addNode(const std::string& name)
{
    for (label i = 0; i < static_cast<label>(nodeNames.size()); ++i)
    {
        if (nodeNames[i] == name) return i;
    }
    nodeNames.push_back(name);
    nodePressureFixed.push_back(std::numeric_limits<double>::quiet_NaN());
    nodeExternalRate.push_back(0.0);
    return static_cast<label>(nodeNames.size()) - 1;
}

label OneDNetwork::nodeIndex(const std::string& name) const
{
    for (label i = 0; i < static_cast<label>(nodeNames.size()); ++i)
    {
        if (nodeNames[i] == name) return i;
    }
    cfFatal("network: unknown node '" + name + "'");
    return -1;
}

label OneDNetwork::addPipe(const NetworkPipe& pipe)
{
    NetworkPipe p = pipe;
    p.from = nodeIndex(p.fromName);
    p.to = nodeIndex(p.toName);
    if (p.from == p.to) cfFatal("network: pipe '" + p.name + "' connects a node to itself");
    // A gas-lift valve is a point element: the port is its `diameter` and it has
    // no length to be positive about, so both numbers belong to the valve scope
    // below rather than to the pipe law's requirement here. A rod pump is a
    // point element too, and it has no bore of its own at all -- the barrel is
    // declared as rodPlungerDiameter -- so `diameter` and `length` are refused
    // there rather than required here.
    if (!pipeIsLiftValve(p) && !pipeIsRodPump(p) && !pipeIsSeparator(p)
        && (!(p.diameter > 0) || !(p.length > 0)))
        cfFatal("network: pipe '" + p.name + "' needs positive length and diameter");
    // 5.223 (H3): the separator row is a point element like the valve and the
    // plunger -- its three nodes ARE the geometry -- so a non-zero length,
    // bore or elevation would be a wall law the row never marches. Resolve the
    // gas outlet node here, next to the from/to resolution, and refuse a row
    // whose gas outlet is one of its own liquid ports.
    if (pipeIsSeparator(p))
    {
        if (p.length != 0.0 || p.diameter != 0.0 || p.elevation != 0.0)
            cfFatal("network: pipe '" + p.name
                    + "' is a separator row with a non-zero length, diameter or "
                      "elevation: the flash is a point element, its three nodes "
                      "are the whole geometry");
        p.sepGasNode = nodeIndex(p.sepGasNodeName);
        if (p.sepGasNode == p.from || p.sepGasNode == p.to)
            cfFatal("network: pipe '" + p.name
                    + "' separator gasNode is one of the row's own inlet or "
                      "liquid outlet: the two outlet streams need two nodes");
    }
    // 5.87: an axis offset is a statement about an annulus. A negative one is
    // the mirror image, but declared it is a typo, so it is refused rather than
    // silently flipped; and without dInner there is no bore to offset.
    if (p.eccentricity < 0)
        cfFatal("network: pipe '" + p.name
                + "' declares a negative eccentricity: an axis offset is a "
                  "distance");
    if (p.eccentricity > 0 && !pipeIsAnnulus(p))
        cfFatal("network: pipe '" + p.name
                + "' declares an eccentricity without an inner bore (dInner): "
                  "an axis offset needs an annulus to offset");
    // 5.88: the override names one of Beggs-Brill's three published holdup
    // laws. Index validity is checked first (a bogus index is a typo whatever
    // the network solves), then scope: the key means nothing on a single-phase
    // network, and nothing as regime 1 -- transition is the map's own
    // interpolation BETWEEN two laws, not a law a caller can order.
    if (p.flowPattern == 1)
        cfFatal("network: pipe '" + p.name
                + "' declares flowPattern 1 (transition): transition is not a "
                  "holdup law, it is the map's own interpolation between "
                  "segregated (0) and intermittent (2) across its L2..L3 band; "
                  "force one of the three published laws (0, 2, 3) or leave "
                  "the key out and let the map choose");
    if (p.flowPattern > 3)
        cfFatal("network: pipe '" + p.name
                + "' declares flowPattern " + std::to_string(p.flowPattern)
                + ": the published regimes are 0 segregated, 2 intermittent, "
                  "3 distributed (1 is the map's blend, refused above), and "
                  "leaving the key out asks the map");
    if (p.flowPattern < -1)
        cfFatal("network: pipe '" + p.name
                + "' declares a negative flowPattern: the only negative value "
                  "this key takes is -1, which is what an absent key means "
                  "(ask the map)");
    if (p.flowPattern >= 0 && !hasBlackOil)
        cfFatal("network: pipe '" + p.name
                + "' declares flowPattern " + std::to_string(p.flowPattern)
                + " on a single-phase network: the flow-pattern map is a "
                  "Beggs-Brill two-phase dispatch surface, a single-phase pipe "
                  "has no pattern to force");
    // 5.143 (G4): the alternative holdup laws replace the holdup block of
    // the Beggs-Brill gradient, so they price a marched black-oil two-phase
    // pipe and nothing else; and flowPattern selects WHICH Beggs-Brill
    // regime law prices the holdup, so with an alternative law there is no
    // Beggs-Brill regime law to select. addPipe owns both scope rules.
    if (p.correlation != 0 && !hasBlackOil)
        cfFatal("network: pipe '" + p.name + "' declares correlation "
                + std::string(p.correlation == 1
                                  ? "hagedornBrown"
                                  : p.correlation == 2
                                        ? "dunsRos"
                                        : "olgasSteady")
                + " on a single-phase network: the alternative holdup laws "
                  "(5.143, 5.155) replace the liquid holdup of the two-phase "
                  "march; this pipe has no two-phase holdup to reprice");
    if (p.correlation != 0 && p.flowPattern >= 0)
        cfFatal("network: pipe '" + p.name + "' declares both flowPattern "
                "and correlation: flowPattern selects which Beggs-Brill "
                "regime law prices the holdup (5.88), and with an "
                "alternative holdup law (5.143) there is no Beggs-Brill "
                "regime law to select - honour one dispatch, not both");
    // 5.144 (G6): the steady thermal march is a MARCHED-pipe declaration -
    // black-oil or pvt-table. A constant-property single-phase pipe is
    // priced by an algebraic law with no march to carry a temperature, and
    // the wall term needs the stream's c_p and a wall to relax toward.
    if (p.thermalOn() && !hasBlackOil && !hasPvtFluid)
        cfFatal("network: pipe '" + p.name + "' declares thermal physics on "
                "a constant-property network: the steady temperature march "
                "(5.144) rides a marched pipe (blackoil or pvt fluid); this "
                "pipe's algebraic law has no march to carry a temperature");
    if (p.thermalOn() && hasPvtFluid && !hasBlackOil
        && !(pvtFluid.temperature > 0.0))
        cfFatal("network: pipe '" + p.name + "' declares thermal physics but "
                "the pvt fluid has no temperature: the march starts from the "
                "fluid's declared inlet temperature (5.144)");
    if (p.wallExchangeOn() && !(p.heatCapacity > 0.0))
        cfFatal("network: pipe '" + p.name + "' declares wall exchange "
                "without a positive heatCapacity: the stream relaxes toward "
                "the wall through m c_p dT/ds = U (T_wall - T)");
    if (p.wallExchangeOn() && p.wallModel == 0 && !(p.wallTemperature > 0.0))
        cfFatal("network: pipe '" + p.name + "' wall exchange with the "
                "constant wall model needs a positive wallTemperature");
    if (p.wallExchangeOn() && p.wallModel == 1
        && !(p.geothermalGradient > 0.0))
        cfFatal("network: pipe '" + p.name + "' wall exchange with the "
                "geothermal wall model needs a positive geothermalGradient");
    if (pipeIsAnnulus(p))
    {
        // 5.73's declaration scope: the annulus is a single-phase
        // Newtonian / power-law steady section. Everything else is refused
        // loudly rather than silently read as a circle.
        if (!(p.annulusInner > 0) || p.annulusInner >= p.diameter)
            cfFatal("network: pipe '" + p.name
                    + "' needs 0 < annulusInner < diameter");
        if (p.yieldStress > 0)
            cfFatal("network: pipe '" + p.name
                    + "' declares an annulus with a yield stress: the yielded "
                      "locus and wall stress are circular-geometry laws; the "
                      "annulus takes no threshold in this increment");
        if (p.coupledResistance > 0)
            cfFatal("network: pipe '" + p.name
                    + "' declares an annulus on a coupled edge: the partner "
                      "owns that edge's geometry");
        // 5.87 (C3 increment 2): the axis offset. It is a bore shape, so it is
        // only meaningful on an annulus, it must leave a gap (the two circles
        // may touch but not cross), and it is refused on a power-law fill
        // because the Fredrickson-Bird profile it replaces is concentric by
        // construction -- its zero-shear radius has no off-axis analogue here.
        if (p.eccentricity > 0)
        {
            if (p.flowIndex > 0)
                cfFatal("network: pipe '" + p.name
                        + "' declares an eccentric annulus with a power-law "
                          "index: the Fredrickson-Bird zero-shear integrals are "
                          "concentric by construction");
            const double maxEcc = 0.5*(p.diameter - p.annulusInner);
            if (p.eccentricity >= maxEcc)
                cfFatal("network: pipe '" + p.name
                        + "' declares an eccentricity " + msgNum(p.eccentricity)
                        + " m >= (diameter - annulusInner)/2 = " + msgNum(maxEcc)
                        + " m: the bores touch and the narrow gap pins the flow");
            double spread = 0;
            const double kg = eccentricAnnulusSectionIntegral
            (
                0.5*p.diameter, 0.5*p.annulusInner, p.eccentricity, spread
            );
            if (!(kg > 0))
                cfFatal("network: pipe '" + p.name
                        + "' failed to resolve its eccentric-annulus profile "
                          "(e/gap = " + msgNum(p.eccentricity/maxEcc)
                        + ", grid spread " + msgNum(spread)
                        + "): this increment refuses an eccentricity its "
                          "bipolar-coordinate grid cannot resolve");
            p.annulusEccK = kg;
        }
        if (p.flowIndex > 0)
        {
            annulusFbIntegrals(p, 1024, p.annulusLambda, p.annulusDimQ);
            // lambda and the integral are RADIUS-based: the walls sit at
            // 0.5 * annulusInner and 0.5 * diameter.
            if (!(p.annulusDimQ > 0)
                || !(p.annulusLambda > 0.5 * p.annulusInner)
                || !(p.annulusLambda < 0.5 * p.diameter))
                cfFatal("network: pipe '" + p.name
                        + "' failed to resolve its Fredrickson-Bird profile");
        }
    }
    if (pipeIsRect(p))
    {
        // 5.80's declaration scope: the rectangular duct is a single-phase
        // Newtonian steady section. Everything else is refused loudly rather
        // than silently read as a circle or an annulus.
        if (pipeIsAnnulus(p))
            cfFatal("network: pipe '" + p.name
                    + "' declares both a rectangular duct (ductHeight) and an "
                      "annulus (dInner): pick one bore shape");
        if (p.yieldStress > 0)
            cfFatal("network: pipe '" + p.name
                    + "' declares a rectangular duct with a yield stress: the "
                      "yielded locus and wall stress are circular-geometry "
                      "laws; the duct takes no threshold in this increment");
        if (p.flowIndex > 0)
            cfFatal("network: pipe '" + p.name
                    + "' declares a rectangular duct with a power-law index: "
                      "the 2D nonlinear cross-section profile is not resolved "
                      "in this increment");
        if (p.coupledResistance > 0)
            cfFatal("network: pipe '" + p.name
                    + "' declares a rectangular duct on a coupled edge: the "
                      "partner owns that edge's geometry");
        const double aspect = std::min(p.diameter, p.ductHeight)
                            / std::max(p.diameter, p.ductHeight);
        p.rectBracket = rectLaminarBracket(aspect);
        if (!(p.rectBracket > 0) || !(p.rectBracket <= 1.0))
            cfFatal("network: pipe '" + p.name
                    + "' failed to resolve its rectangular-duct laminar "
                      "profile");
    }
    // Gas-lift valve (5.89, D1 increment 3): the declaration is a valve, not a
    // pipe with a valve in it, so everything the pipe laws would otherwise read
    // is refused here rather than left in place unread. The three numbers the
    // closures actually consume are checked where they can be checked: R in
    // (0,1) because 1-R divides the force balance, kappa > 1 because the
    // critical ratio's exponent is 1/(k-1), Cd in (0,1] because it is a
    // recovery fraction, and a non-negative spring because St only ever pushes
    // the valve shut.
    if (pipeIsLiftValve(p))
    {
        if (!(p.diameter > 0))
            cfFatal("network: pipe '" + p.name
                    + "' is a gas-lift valve without a port: `diameter` is the "
                      "bore the lift gas passes through");
        if (p.length != 0)
            cfFatal("network: pipe '" + p.name
                    + "' is a gas-lift valve with a length: a valve is a point "
                      "element, its drop is the port's, not a bore's friction");
        if (p.elevation != 0)
            cfFatal("network: pipe '" + p.name
                    + "' is a gas-lift valve with an elevation: the force "
                      "balance reads the two node pressures directly");
        if (p.roughness > 0 || p.valveK > 0 || p.pumpHead > 0
            || p.pumpGasFraction > 0 || p.coupledResistance > 0
            || p.annulusInner > 0 || p.ductHeight > 0 || p.eccentricity > 0
            || p.flowIndex > 0 || p.consistency > 0 || p.yieldStress > 0
            || p.flowPattern >= 0 || p.rodPlungerDiameter > 0
            || p.correlation != 0 || p.thermalOn())
            cfFatal("network: pipe '" + p.name
                    + "' is a gas-lift valve that also declares a wall law, a "
                      "bore shape, a pump, a coupled edge, a flow-pattern "
                      "override, a holdup law or thermal physics: the valve's "
                      "rate comes from its bellows and its port and nothing "
                      "else prices it");
        if (!(p.liftAreaRatio > 0) || !(p.liftAreaRatio < 1))
            cfFatal("network: pipe '" + p.name
                    + "' has a lift-valve area ratio outside (0,1): R = Ap/Ab "
                      "is a port smaller than the bellows that senses it, and "
                      "the force balance divides by 1-R");
        if (p.liftSpringPressure < 0)
            cfFatal("network: pipe '" + p.name
                    + "' has a negative lift-valve spring pressure: the spring "
                      "only ever helps hold the valve shut");
        if (!(p.liftDischargeCoeff > 0) || !(p.liftDischargeCoeff <= 1))
            cfFatal("network: pipe '" + p.name
                    + "' has a lift-valve discharge coefficient outside (0,1]: "
                      "Cd recovers part of the ideal nozzle flux, it cannot "
                      "create more of it");
        if (!(p.liftKappa > 1))
            cfFatal("network: pipe '" + p.name
                    + "' has a lift-valve kappa <= 1: the isentropic nozzle "
                      "exponent k = cp/cv exceeds 1 for any real gas, and the "
                      "critical ratio (2/(k+1))^(k/(k-1)) is priced by k/(k-1)");
    }
    // Rod pump (5.90, D1 increment 4): the declaration is a pump, not a pipe
    // with a pump in it, so the bore, the length and the elevation the pipe
    // laws would otherwise read belong to the row and are refused here. The
    // seven numbers the card consumes are each checked where they divide: the
    // rod area divides the stretch, the stroke divides the efficiency, and a
    // string lighter than the fluid it hangs in would make the buoyed weight a
    // statement about flotation, not about rods -- that one is allowed and the
    // card says so, because the arithmetic does not care which is heavier.
    if (pipeIsRodPump(p))
    {
        if (p.length != 0)
            cfFatal("network: pipe '" + p.name
                    + "' is a rod pump with a length: the pump is a point "
                      "element, its lift is the columns above and below its "
                      "two ports, not a bore to march");
        if (p.diameter > 0)
            cfFatal("network: pipe '" + p.name
                    + "' is a rod pump with a diameter: the bore that matters "
                      "is the barrel's, declared as rodPlungerDiameter; a "
                      "pipe bore here would be read by nothing");
        if (p.elevation != 0)
            cfFatal("network: pipe '" + p.name
                    + "' is a rod pump with an elevation: the setting depth is "
                      "rodLength, and the two ports are at the same depth");
        if (p.roughness > 0 || p.valveK > 0 || p.pumpHead > 0
            || p.pumpRate > 0 || p.pumpStages > 1 || p.pumpSpeed != 1
            || p.pumpGasFraction > 0 || p.pumpIntakePressure > 0
            || p.coupledResistance > 0 || p.annulusInner > 0
            || p.ductHeight > 0 || p.eccentricity > 0 || p.flowIndex > 0
            || p.consistency > 0 || p.yieldStress > 0 || p.flowPattern >= 0
            || p.liftDomePressure > 0 || p.correlation != 0 || p.thermalOn())
            cfFatal("network: pipe '" + p.name
                    + "' is a rod pump that also declares a wall law, a bore "
                      "shape, a head-curve pump, a lift valve, a coupled edge, "
                      "a flow-pattern override, a holdup law or thermal "
                      "physics: the plunger's rate comes from its stroke and "
                      "its string and nothing else prices it");
        if (!(p.rodStrokeLength > 0))
            cfFatal("network: pipe '" + p.name
                    + "' is a rod pump without a stroke: rodStrokeLength is the "
                      "polished-rod travel the card is drawn against, and it "
                      "divides the volumetric efficiency");
        if (!(p.rodStrokesPerMinute > 0))
            cfFatal("network: pipe '" + p.name
                    + "' is a rod pump without a speed: rodStrokesPerMinute is "
                      "what turns the per-stroke card into a rate and a power");
        if (!(p.rodDiameter > 0))
            cfFatal("network: pipe '" + p.name
                    + "' is a rod pump without a rod: rodDiameter is the string "
                      "that carries the load, and its area divides the stretch");
        if (!(p.rodLength > 0))
            cfFatal("network: pipe '" + p.name
                    + "' is a rod pump without a setting depth: rodLength is "
                      "both the string that stretches and the column the "
                      "plunger holds up");
        if (!(p.rodYoungModulus > 0))
            cfFatal("network: pipe '" + p.name
                    + "' is a rod pump without a modulus: Hooke's law prices "
                      "the stroke the fluid load costs");
        if (!(p.rodDensity > 0))
            cfFatal("network: pipe '" + p.name
                    + "' is a rod pump with a non-positive rod density: the "
                      "string's weight in air is what the buoyancy subtracts "
                      "from");
        if (p.rodDiameter >= p.rodPlungerDiameter)
            cfFatal("network: pipe '" + p.name
                    + "' has a rod at least as wide as its plunger: the rod "
                      "passes through the barrel's bore, so A_r < A_p is what "
                      "leaves an annulus to sweep");
    }
    pipes.push_back(p);
    return static_cast<label>(pipes.size()) - 1;
}

double churchillFriction(double reynolds, double relativeRoughness)
{
    const double re = std::max(reynolds, 1e-12);
    const double a1 = std::pow(7.0 / re, 0.9) + 0.27 * relativeRoughness;
    if (a1 <= 0) cfFatal("churchillFriction: invalid Re/roughness combination");
    const double a = std::pow(2.457 * std::log(1.0 / a1), 16.0);
    const double b = std::pow(37530.0 / re, 16.0);
    return 8.0 * std::pow(std::pow(8.0 / re, 12.0) + std::pow(a + b, -1.5), 1.0 / 12.0);
}

namespace
{

/// Rod pump (5.90, D1 increment 4): the one row in this kernel whose rate is
/// not an answer. Every other pipe is a law Q(dp) the node solve reads; a
/// plunger moves A_p S_eff N because that is what the machine is set to, and
/// what it reports back is the pressure its two ports differ by. So the march
/// is skipped entirely, the row's `flowRate` is imposed in the declared
/// direction (which is why `dpTotal` reads negative while the pump lifts), and
/// the card is priced at the differential the network happened to hand the
/// plunger. Nothing is assumed about the column above the pump: `fluidLoad` is
/// A_p times the two node pressures, so the tubing's weight, its friction, the
/// wellhead backpressure, the annulus level and the casing pressure all arrive
/// already added up -- and the card's own `fluidLoadMax` is what that same
/// differential would be with none of them but a full wet string, which is the
/// difference a level shot reads.
///
/// The conductance this row offers Newton is the slope of that imposed rate,
/// and it is the stretch: a pump lifting against a heavier column travels a
/// shorter stroke, so dQ/d(dpTotal) is positive and small -- three to four
/// orders under the marched rows around it, which is what a positive-displacement
/// machine is. A stalled or unloaded pump has no slope at all, and keeps the
/// same 1e-12 floor the shut valve keeps: stiff enough to hold the Jacobian
/// open, loose enough not to pretend to be a pipe.
void solveRodPumpConstant(const NetworkPipe& pipe, const FluidProperties& fluid,
                          double dpTotal, PipeFlowResult& out)
{
    const RodPump pump{pipe.rodPlungerDiameter, pipe.rodStrokeLength,
                       pipe.rodStrokesPerMinute, pipe.rodDiameter,
                       pipe.rodLength, pipe.rodYoungModulus, pipe.rodDensity};
    // from is the intake, to is the tubing above the pump: the plunger works
    // against p(to) - p(from), the negation of the row's own dpTotal.
    const RodPumpCard card =
        rodPumpCardOfLoad(pump, fluid.rho, -dpTotal);

    out.name = pipe.name;
    out.regime = -1;
    out.dpTotal = dpTotal;
    out.dpFriction = 0.0;  // a point element has no wall to rub against
    out.dpPump = dpTotal;  // and the whole rise is the machine's
    out.flowRate = card.displacement;
    out.massRate = card.massRate;
    // The face velocity the barrel sweeps, referred to the bore that does the
    // sweeping -- not to `diameter`, which a pump row refuses.
    out.velocity = card.plungerArea > 0
                       ? card.displacement / card.plungerArea : 0.0;
    out.frictionFactor = 0.0; // no wall law is priced here

    out.card.state = static_cast<int>(card.state);
    out.card.plungerArea = card.plungerArea;
    out.card.rodArea = card.rodArea;
    out.card.rodWeightInAir = card.rodWeightInAir;
    out.card.buoyancy = card.buoyancy;
    out.card.rodWeightInFluid = card.rodWeightInFluid;
    out.card.fluidLoad = card.fluidLoad;
    out.card.fluidLoadMax = card.fluidLoadMax;
    out.card.upthrust = card.upthrust;
    out.card.maxLoad = card.maxLoad;
    out.card.minLoad = card.minLoad;
    out.card.downholeMaxLoad = card.downholeMaxLoad;
    out.card.downholeMinLoad = card.downholeMinLoad;
    out.card.counterbalance = card.counterbalance;
    out.card.stretch = card.stretch;
    out.card.stroke = card.stroke;
    out.card.effectiveStroke = card.effectiveStroke;
    out.card.efficiency = card.efficiency;
    out.card.displacement = card.displacement;
    out.card.massRate = card.massRate;
    out.card.workPerCycle = card.workPerCycle;
    out.card.hydraulicPower = card.hydraulicPower;

    const double h = std::max(1e-4 * std::abs(dpTotal), 1.0);
    const double qPlus = rodPumpCardOfLoad(pump, fluid.rho,
                                           -(dpTotal + h)).displacement;
    const double qMinus = rodPumpCardOfLoad(pump, fluid.rho,
                                            -(dpTotal - h)).displacement;
    const double c = (qMinus - qPlus) / (2.0 * h);
    out.conductance = (c > 0) ? c : 1e-12;
}

} // namespace

void solvePipeFlow(const NetworkPipe& pipe, const FluidProperties& fluid,
                   double dpTotal, PipeFlowResult& out)
{
    if (pipeIsLiftValve(pipe))
        cfFatal("solvePipeFlow: pipe '" + pipe.name + "' is a gas-lift valve on "
                "a constant-property fill: the nozzle is priced on the gas "
                "density at the valve's own upstream pressure");
    if (pipeIsRodPump(pipe))
    {
        solveRodPumpConstant(pipe, fluid, dpTotal, out);
        return;
    }
    double Hq = 0;
    const RawFlow r = equipmentFlow(pipe, fluid, dpTotal, Hq);
    out.name = pipe.name;
    out.dpTotal = dpTotal;
    out.dpFriction = dpTotal - fluid.rho * kGravity * pipe.elevation;
    out.dpValve = pipe.valveK * fluid.rho * std::abs(r.velocity) * r.velocity
                / 2.0;
    out.dpPump = -fluid.rho * kGravity * Hq;
    out.velocity = r.velocity;
    out.flowRate = r.velocity * pipeArea(pipe);
    // 5.191: the PCP's own share of that rise, reported on the rate this row
    // actually carries, so a checker can read dpPcp against (nV - q)/C itself.
    if (pipe.pcpOn)
        out.dpPcp = -fluid.rho * kGravity
                  * pcpHeadOf(pipe, out.flowRate, fluid.rho);
    // 5.192: the jet pump's own rise, its entrainment ratio and its state --
    // all read from the one JetRow the row law already priced.
    if (pipe.jetOn)
    {
        const JetRow jr = jetRowOf(pipe, out.flowRate, fluid.rho);
        out.dpJet = -fluid.rho * kGravity * jr.head;
        out.jetFlowRatio = jr.ratio;
        out.jetState = jr.state;
    }
    out.reynolds = r.reynolds;
    out.frictionFactor = r.friction;
    out.transition = yieldedTransitionOf(pipe, fluid, std::abs(out.flowRate));

    // conductance dQ/d(dpTotal) by central difference on the monotone law;
    // floor/replace with the laminar slope if the quotient degenerates
    const double scale = std::max(std::abs(dpTotal),
                                  fluid.rho * kGravity * std::abs(pipe.elevation));
    const double h = std::max(1e-6 * std::abs(out.dpFriction), 1e-7 * scale);
    // for a coupled element the slope of its own law is the exact bound; the
    // Hagen-Poiseuille limit of the pipe geometry is not its physics
    const double lam = pipe.coupledResistance > 0
                           ? 1.0 / pipe.coupledResistance
                           : laminarConductance(pipe, fluid);
    if (h > 0)
    {
        double Hp, Hm;
        const RawFlow rp = equipmentFlow(pipe, fluid, dpTotal + h, Hp);
        const RawFlow rm = equipmentFlow(pipe, fluid, dpTotal - h, Hm);
        const double g = (rp.velocity - rm.velocity) * pipeArea(pipe) / (2.0 * h);
        out.conductance = (std::isfinite(g) && g > 0) ? std::min(g, lam) : lam;
    }
    else
    {
        out.conductance = lam;
    }
}

namespace
{

/// 5.223 (H3): the separator row's own statement, priced at the pressures the
/// Newton iterate currently holds. The connection chosen for this increment is
/// "one inlet node, two outlet edges" rather than "the separator as a node":
/// the arriving stream already sits at the inlet node, so flashing it there
/// needs no extra pressure unknown the dense Jacobian has no slot for, while
/// a stage node would be a fourth node the row's three declared ports must
/// then be re-derived from. The two outlets are near-equal by construction
/// (both are the stage minus dP), which is exactly what the two-node pipe
/// model cannot state and what makes this a row of its own.
///
/// The equality is not a penalty: the outlet node pressures are DERIVED
/// unknowns, p(outlet) := p(from) - dP, identities of the linear system, so
/// the solved upstream-downstream difference IS dP to machine precision and
/// the Newton matrix stays as well-conditioned as the marched rows alone (a
/// stiff penalty link was tried first and rejected: the conductance needed
/// for dp == dP swallows the node's own marched conductance and the matrix
/// goes numerically singular). dP enters the Jacobian through the chain rule
/// (dp_outlet/dp_from = 1): every row touching an outlet node contributes its
/// conductance to the inlet node's equation, which is what makes the inlet
/// pressure the unknown the whole outlet system solves against. The row books
/// no residual of its own: each edge carries the demand the lines impose at
/// the clamped pressure, the two demands sum to the arriving feed q_in, and
/// the flash prices that feed -- mass closes by separatorSolve's own
/// construction, volume closes because the edges ARE the line demands.
void separatorRowOf(const OneDNetwork& net, const NetworkPipe& pipe,
                    const std::vector<double>& pAll,
                    double qEdgeLiq, double qEdgeGas, PipeFlowResult& out)
{
    const double pIn = pAll[static_cast<std::size_t>(pipe.from)];
    const double pOut = pAll[static_cast<std::size_t>(pipe.to)];
    const double tStage = pipe.sepTemperature > 0.0
                            ? pipe.sepTemperature
                            : net.blackOil.temperature;
    const double feed = qEdgeLiq + qEdgeGas;
    // Mid-Newton the iterate may transiently reverse the feed; the flash is
    // only asked for the report then (a positive stand-in keeps separatorSolve
    // quiet), while the booking stays exact. A genuinely feed-less solution is
    // refused after convergence.
    const double feedReport = feed > 1e-300 ? feed : 1e-300;

    // The one allocation function the point element already owns, called at
    // the stage pressure the network solved -- not a second flash model.
    const SeparatorSolution s = separatorSolve(net.blackOil, feedReport,
                                               pIn, tStage, pipe.sepEfficiency,
                                               pipe.sepDP);

    out = PipeFlowResult{};
    out.name = pipe.name;
    out.regime = -1;
    out.dpTotal = pIn - pOut;
    out.dpFriction = 0.0; // a point element has no wall to rub against
    out.flowRate = feed;  // the arriving q_in the flash saw
    out.massRate = s.feedMassRate;
    out.conductance = 0.0; // the row books nothing: its outlets are derived
    out.separatorState = 0;
    out.sepStagePressure = s.pressure;
    out.sepOutletPressure = s.outletPressure;
    out.sepRs = s.rs;
    out.sepBo = s.bo;
    out.sepBg = s.bg;
    out.sepGasLiberated = s.gasLiberated;
    out.sepFeedMassRate = s.feedMassRate;
    out.sepLiquidMassRate = s.liquidMassRate;
    out.sepGasMassRate = s.gasMassRate;
    out.sepEntrainedMassRate = s.entrainedMassRate;
    out.sepClosure = s.closure;
    out.sepLiquidRate = qEdgeLiq;
    out.sepGasEdgeRate = qEdgeGas;
    out.sepGasStdRate = s.gasLiberated * feedReport;
}

} // namespace

SteadyNetworkResult solveSteadyNetwork(const OneDNetwork& net,
                                       const SteadyControls& ctl)
{
    SteadyNetworkResult res;
    const std::size_t n = net.nNodes();
    res.nodePressure.assign(n, 0.0);
    res.pipeFlows.resize(net.pipes.size());

    // 5.223 (H3): the separator row's outlet nodes are DERIVED, not unknowns:
    // p(outlet) := p(inlet) - dP is an identity of the linear system, so the
    // solved upstream-downstream difference IS dP to machine precision while
    // the Newton matrix stays as well-conditioned as the marched rows alone.
    // The chain rule (dp_outlet/dp_inlet = 1) is what puts dP into the
    // Jacobian: every row touching an outlet node contributes its conductance
    // to the inlet node's equation. A node may be derived by one row only,
    // cannot be fixed-pressure, cannot carry an external rate (it has no
    // balance equation to receive it) and cannot itself feed a separator.
    std::vector<label> master(n);
    std::vector<double> derivedDrop(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) master[i] = static_cast<label>(i);
    for (const NetworkPipe& pipe : net.pipes)
    {
        if (!pipe.separatorOn) continue;
        const label ports[2] = {pipe.to, pipe.sepGasNode};
        for (const label port : ports)
        {
            const std::size_t i = static_cast<std::size_t>(port);
            if (!std::isnan(net.nodePressureFixed[i]))
                cfFatal("steadyNetwork: separator row '" + pipe.name
                        + "' outlet node is declared with a fixed pressure:"
                          " the row derives it from the inlet minus dP");
            if (net.nodeExternalRate[i] != 0.0)
                cfFatal("steadyNetwork: separator row '" + pipe.name
                        + "' outlet node carries an external rate: a derived"
                          " node has no balance equation to receive it");
            if (master[i] != static_cast<label>(i))
                cfFatal("steadyNetwork: separator row '" + pipe.name
                        + "' outlet node is already derived by another"
                          " separator row");
            master[i] = pipe.from;
            derivedDrop[i] = pipe.sepDP;
        }
        if (master[static_cast<std::size_t>(pipe.from)] != pipe.from)
            cfFatal("steadyNetwork: separator row '" + pipe.name
                    + "' inlet node is itself a derived outlet: chained"
                      " separators are not this increment's scope");
    }

    // unknown node pressures
    std::vector<label> unk(n, -1);
    std::size_t nFree = 0;
    double pSum = 0;
    std::size_t nFixed = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (std::isnan(net.nodePressureFixed[i])
            && master[i] == static_cast<label>(i))
            unk[i] = static_cast<label>(nFree++);
        else if (!std::isnan(net.nodePressureFixed[i]))
            { pSum += net.nodePressureFixed[i]; ++nFixed; }
    }
    if (nFixed == 0 && nFree > 0)
        cfFatal("steadyNetwork: no fixed-pressure node (no reference)");
    const double pInit = nFixed ? pSum / static_cast<double>(nFixed) : 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        res.nodePressure[i] = std::isnan(net.nodePressureFixed[i])
                                  ? pInit : net.nodePressureFixed[i];
    }

    // total injected rate = residual normalisation scale
    double totalIn = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (net.nodeExternalRate[i] > 0) totalIn += net.nodeExternalRate[i];
    }
    const double qScale = std::max(totalIn, 1e-300);

    // Newton on free-node pressures
    std::vector<double> p(nFree, pInit);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (unk[i] >= 0) p[static_cast<std::size_t>(unk[i])] = res.nodePressure[i];
    }

    auto evaluate = [&](const std::vector<double>& pv,
                        std::vector<double>& r,
                        std::vector<std::vector<double>>& j,
                        bool withJacobian) -> double
    {
        std::vector<double> pAll(n, pInit);
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!std::isnan(net.nodePressureFixed[i]))
                pAll[i] = net.nodePressureFixed[i];
            else if (unk[i] >= 0)
                pAll[i] = pv[static_cast<std::size_t>(unk[i])];
        }
        // 5.223 (H3): the separator outlets, one level deep (masters are
        // never themselves derived -- refused above).
        for (std::size_t i = 0; i < n; ++i)
        {
            if (unk[i] < 0 && std::isnan(net.nodePressureFixed[i]))
                pAll[i] = pAll[static_cast<std::size_t>(master[i])]
                        - derivedDrop[i];
        }
        r.assign(nFree, 0.0);
        if (withJacobian) j.assign(nFree, std::vector<double>(nFree, 0.0));

        double maxRes = 0;
        // 5.223 (H3): the marched rows book to the MASTER unknown -- a derived
        // outlet node's balance equation lives at its separator inlet, and
        // unk[master[.]] is exactly the chain rule dp_outlet/dp_inlet = 1 that
        // puts dP into the Jacobian. They also accumulate the net volume each
        // node hands to its lines, which is what the separator edges carry.
        std::vector<double> nodeTake(n, 0.0);
        for (std::size_t k = 0; k < net.pipes.size(); ++k)
        {
            const NetworkPipe& pipe = net.pipes[k];
            if (pipe.separatorOn) continue;
            PipeFlowResult fr;
            if (net.hasBlackOil)
            {
                solvePipeFlowTwoPhase(pipe, net.blackOil,
                                      pAll[static_cast<std::size_t>(pipe.from)],
                                      pAll[static_cast<std::size_t>(pipe.to)], fr,
                                      net.sigmaOfP ? &net.sigmaOfP : nullptr);
            }
            else if (net.hasPvtFluid)
            {
                solvePipeFlowPvt(pipe, net.pvtFluid,
                                 pAll[static_cast<std::size_t>(pipe.from)],
                                 pAll[static_cast<std::size_t>(pipe.to)], fr);
            }
            else
            {
                solvePipeFlow(pipe, net.fluid,
                              pAll[static_cast<std::size_t>(pipe.from)]
                            - pAll[static_cast<std::size_t>(pipe.to)], fr);
            }
            if (withJacobian) res.pipeFlows[k] = fr;
            nodeTake[static_cast<std::size_t>(pipe.from)] += fr.flowRate;
            nodeTake[static_cast<std::size_t>(pipe.to)] -= fr.flowRate;
            const label a = unk[static_cast<std::size_t>(
                master[static_cast<std::size_t>(pipe.from)])];
            const label b = unk[static_cast<std::size_t>(
                master[static_cast<std::size_t>(pipe.to)])];
            // R_a -= Q ; R_b += Q
            if (a >= 0) r[static_cast<std::size_t>(a)] -= fr.flowRate;
            if (b >= 0) r[static_cast<std::size_t>(b)] += fr.flowRate;
            if (withJacobian)
            {
                const double c = fr.conductance;
                if (a >= 0)
                {
                    j[a][a] -= c;
                    if (b >= 0) { j[a][b] += c; j[b][a] += c; j[b][b] -= c; }
                }
            }
        }
        for (std::size_t k = 0; k < net.pipes.size(); ++k)
        {
            const NetworkPipe& pipe = net.pipes[k];
            if (!pipe.separatorOn) continue;
            // 5.223 (H3): the row books no residual of its own -- its outlet
            // pressures are derived identities and the edge flows ARE the line
            // demands at those pressures. It states the flash of the arriving
            // feed for the report (and for the feed rejection below).
            PipeFlowResult fr;
            separatorRowOf(net, pipe, pAll,
                           nodeTake[static_cast<std::size_t>(pipe.to)],
                           nodeTake[static_cast<std::size_t>(pipe.sepGasNode)],
                           fr);
            if (withJacobian) res.pipeFlows[k] = fr;
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            if (unk[i] >= 0) r[static_cast<std::size_t>(unk[i])] += net.nodeExternalRate[i];
        }
        for (double v : r) maxRes = std::max(maxRes, std::abs(v));
        return maxRes;
    };

    std::vector<double> r;
    std::vector<std::vector<double>> j;
    double maxRes = evaluate(p, r, j, true);
    res.residual = maxRes;
    res.relativeResidual = maxRes / qScale;

    int iter = 0;
    while (res.relativeResidual > ctl.tol && iter < ctl.maxIter)
    {
        if (nFree > 0)
        {
            std::vector<double> negR(r.size());
            for (std::size_t i = 0; i < r.size(); ++i) negR[i] = -r[i];
            detail::solveDenseLinear(j, negR);
            const std::vector<double> pOld = p;
            // Globalization along the Newton ray. The two-phase law has a
            // static-head knee: below it the residual is a flat plateau
            // (Q = 0), so no strict decrease exists until the knee is
            // crossed. Accept equal-residual trials on the plateau - the
            // direction there consistently pushes the source node up and
            // the iterate climbs to the knee - but require strict decrease
            // once live. Trials leaving the physical domain (non-finite or
            // non-positive pressure) are skipped.
            bool moved = false;
            for (int step = 0; step < 16 && !moved; ++step)
            {
                const double t = step < 12 ? std::pow(0.5, step)
                                           : std::pow(2.0, step - 11);
                std::vector<double> pTry = pOld;
                bool feasible = true;
                for (std::size_t i = 0; i < nFree; ++i)
                {
                    pTry[i] += ctl.relax * negR[i] * t;
                    if (!(pTry[i] > 0.0) || !std::isfinite(pTry[i]))
                        feasible = false;
                }
                // 5.223 (H3): derived outlet nodes ride their inlet's trial
                // pressure minus dP -- keep THOSE positive too, or a trial
                // state would march a pipe at a non-positive pressure.
                for (std::size_t i = 0; i < n && feasible; ++i)
                {
                    if (unk[i] >= 0 || !std::isnan(net.nodePressureFixed[i]))
                        continue;
                    const double pu = pTry[static_cast<std::size_t>(
                        unk[static_cast<std::size_t>(master[i])])];
                    if (!(pu - derivedDrop[i] > 0.0)) feasible = false;
                }
                if (!feasible) continue;
                maxRes = evaluate(pTry, r, j, true);
                if (maxRes <= res.residual)
                {
                    p = pTry;
                    moved = true;
                }
            }
            if (!moved)
            {
                // trials overwrote r/j with the last rejected trial's
                // state; restore them for the point we keep
                p = pOld;
                maxRes = evaluate(p, r, j, true);
            }
            res.residual = maxRes;
            res.relativeResidual = maxRes / qScale;
        }
        ++iter;
        if (nFree == 0) break;
    }
    res.iterations = iter;
    res.converged = res.relativeResidual <= ctl.tol;

    for (std::size_t i = 0; i < n; ++i)
    {
        if (unk[i] >= 0) res.nodePressure[i] = p[static_cast<std::size_t>(unk[i])];
    }
    // 5.223 (H3): the derived outlets report the pressure the identity gave
    // them, not the initial guess they were never solved for.
    for (std::size_t i = 0; i < n; ++i)
    {
        if (unk[i] < 0 && std::isnan(net.nodePressureFixed[i]))
            res.nodePressure[i] = res.nodePressure[
                static_cast<std::size_t>(master[i])] - derivedDrop[i];
    }
    // final pipe states at the converged pressures
    for (std::size_t k = 0; k < net.pipes.size(); ++k)
    {
        const NetworkPipe& pipe = net.pipes[k];
        if (pipe.separatorOn) continue; // 5.223: booked in the pass below,
                                        // once every node's conductance is final
        if (net.hasBlackOil)
        {
            solvePipeFlowTwoPhase(pipe, net.blackOil,
                                  res.nodePressure[static_cast<std::size_t>(pipe.from)],
                                  res.nodePressure[static_cast<std::size_t>(pipe.to)],
                                  res.pipeFlows[k],
                                  net.sigmaOfP ? &net.sigmaOfP : nullptr);
        }
        else if (net.hasPvtFluid)
        {
            solvePipeFlowPvt(pipe, net.pvtFluid,
                             res.nodePressure[static_cast<std::size_t>(pipe.from)],
                             res.nodePressure[static_cast<std::size_t>(pipe.to)],
                             res.pipeFlows[k]);
        }
        else
        {
            solvePipeFlow(pipe, net.fluid,
                          res.nodePressure[static_cast<std::size_t>(pipe.from)]
                        - res.nodePressure[static_cast<std::size_t>(pipe.to)],
                          res.pipeFlows[k]);
        }
    }
    // 5.223 (H3): the separator rows at the converged pressures, and the
    // increment's own rejection: a row whose inlet node delivered nothing has
    // no feed to flash.
    {
        std::vector<double> nodeTake(n, 0.0);
        for (std::size_t k = 0; k < net.pipes.size(); ++k)
        {
            const NetworkPipe& pipe = net.pipes[k];
            if (pipe.separatorOn) continue;
            nodeTake[static_cast<std::size_t>(pipe.from)]
                += res.pipeFlows[k].flowRate;
            nodeTake[static_cast<std::size_t>(pipe.to)]
                -= res.pipeFlows[k].flowRate;
        }
        for (std::size_t k = 0; k < net.pipes.size(); ++k)
        {
            const NetworkPipe& pipe = net.pipes[k];
            if (!pipe.separatorOn) continue;
            separatorRowOf(net, pipe, res.nodePressure,
                           nodeTake[static_cast<std::size_t>(pipe.to)],
                           nodeTake[static_cast<std::size_t>(pipe.sepGasNode)],
                           res.pipeFlows[k]);
            if (res.converged && !(res.pipeFlows[k].flowRate > 0.0))
                cfFatal("steadyNetwork: separator row '" + pipe.name
                        + "' solved with no feed rate (q_in <= 0): the inlet"
                          " node delivers nothing for the flash to split");
        }
    }
    return res;
}

namespace
{

/// Largest root of g(q) = march(pUp, q).pOutlet - pLow via Illinois. The
/// BB gradient is U-shaped in q - holdup relief lowers the required drop
/// before friction raises it - so g can vanish several times and the
/// stable operating point is the largest root. Bracketing from q=0
/// instead latches onto the near-zero root whenever the static surplus is
/// smaller than the low-rate friction hump, and the pipe law then jumps
/// between branches as the pressures move. `choked` when the march hit
/// the floor at the upper bracket.
double twoPhaseRateForDrop(const NetworkPipe& pipe,
                           const pvt::BlackOilFluid& fluid, double pUp,
                           double pLow, double sinFlow, bool& choked,
                           const Rheology& liquidRheology,
                           const std::function<double(double)>* sigmaOfP = nullptr)
{
    choked = false;
    const double area = 0.25 * M_PI * pipe.diameter * pipe.diameter;
    const double pFloor = std::max(0.02 * pUp, 1000.0);
    // 5.193 (G9): a surveyed pipe prices its gravity segment by segment; the
    // undeclared pipe hands nullptr and keeps the single slope, bit for bit.
    const std::vector<double> segSinVec =
        pipe.surveyOn ? surveySegSinOf(pipe, 24) : std::vector<double>();
    const std::vector<double>* segSin = pipe.surveyOn ? &segSinVec : nullptr;
    auto g = [&](double q)
    {
        const BeggsBrillMarchOut m = marchBeggsBrillPipe(
            fluid, q, pUp, pipe.length, pipe.diameter, pipe.roughness,
            sinFlow, 24, liquidRheology, pipe.flowPattern, sigmaOfP,
            pipe.correlation, nullptr, pipe.hbGriffith, segSin);
        // choked marches return the floor shortfall: continuous with the
        // unchoked branch (there pOutlet -> pFloor) and correctly scaled
        // for the secant, unlike a fixed sentinel
        return m.ok ? m.pOutlet - pLow : pFloor - pLow;
    };

    // climb the rate axis until a negative sample sits above a positive
    // one; that bracket ends at the largest root
    const double qSeed = std::max(1e-3 * area, 1e-8);
    const double g0 = g(0.0);
    double qLo = 0.0, fLo = g0;
    double qHi = qSeed, fHi = g(qHi);
    int expansions = 0;
    while (!(fHi <= 0.0 && qLo > 0.0))
    {
        if (fHi > 0.0)
        {
            qLo = qHi;
            fLo = fHi;
        }
        else if (qLo <= 0.0 && expansions >= 20)
            break; // no gas-lift dip within reach
        qHi *= 4.0;
        fHi = g(qHi);
        if (++expansions > 60)
        {
            choked = true; // pathological: no rate within range
            return qHi;
        }
    }
    if (qLo <= 0.0)
    {
        // nothing positive above zero: either the static column alone
        // reaches the low side (root near zero) or the pipe is dead
        if (g0 <= 0.0) return 0.0;
        fLo = g0;
    }

    // Illinois secant between qLo (g>0) and qHi (g<0)
    double q = qLo + (qHi - qLo) * fLo / (fLo - fHi);
    for (int it = 0; it < 100; ++it)
    {
        const double fq = g(q);
        if (std::abs(fq) < 1e-7 * pUp) break;
        if (fq > 0.0)
        {
            qLo = q;
            fLo = fq;
            // Illinois: bias the retained endpoint against re-selection
            fLo *= 0.5;
        }
        else
        {
            qHi = q;
            fHi = fq;
            fHi *= 0.5;
        }
        const double qNew = qLo + (qHi - qLo) * fLo / (fLo - fHi);
        if (!(qNew > qLo && qNew < qHi)) q = 0.5 * (qLo + qHi);
        else q = qNew;
        if (std::abs(qHi - qLo) < 1e-9 * std::max(qHi, 1e-30)) break;
    }
    return q;
}

} // namespace

bool solvePipeFlowTwoPhase(const NetworkPipe& pipe,
                           const pvt::BlackOilFluid& fluid, double pFrom,
                           double pTo, PipeFlowResult& out,
                           const std::function<double(double)>* sigmaOfP)
{
    if (pipeIsAnnulus(pipe))
        cfFatal("twoPhase: pipe '" + pipe.name + "' declares an annulus: the "
                "black-oil march is a circular-bore correlation chain");
    if (pipeIsLiftValve(pipe))
        cfFatal("twoPhase: pipe '" + pipe.name + "' is a gas-lift valve: the "
                "black-oil march carries its own solution gas, and this "
                "increment prices the injection string, not the two-phase "
                "column");
    if (pipeIsRect(pipe))
        cfFatal("twoPhase: pipe '" + pipe.name + "' declares a rectangular "
                "duct: the black-oil march is a circular-bore correlation "
                "chain");
    out = PipeFlowResult{};
    out.name = pipe.name;
    out.regime = -1;
    const double dp = pFrom - pTo;
    out.dpTotal = dp;
    if (!(pipe.length > 0) || !(pipe.diameter > 0))
        return true;

    // 5.76 (A5 increment 2): the case's liquid rheology is the one the
    // black-oil march prices on the mixture wall rate. The three
    // wall-law fields stay on the pipe (no new grammar); every consumer
    // in this function sees the same declaration.
    const Rheology liquidRheology{
        pipe.flowIndex, pipe.consistency, pipe.yieldStress};

    const double sinFwd = pipe.elevation / pipe.length;
    bool forward = dp >= 0.0; // flow from -> to when pFrom is higher
    bool choked = false;
    double q = 0.0;
    if (dp >= 0.0)
    {
        q = twoPhaseRateForDrop(pipe, fluid, pFrom, pTo, sinFwd, choked,
                                liquidRheology, sigmaOfP);
    }
    else
    {
        // A pressure difference alone cannot express a downhill
        // gravity-driven line. Price both directions once and let the
        // direction with the larger feasible rate carry the pipe; the old
        // high-to-low choice survives whenever the pressure head wins.
        bool chokedF = false, chokedB = false;
        const double qF = twoPhaseRateForDrop(pipe, fluid, pFrom, pTo,
                                              sinFwd, chokedF,
                                              liquidRheology, sigmaOfP);
        const double qB = twoPhaseRateForDrop(pipe, fluid, pTo, pFrom,
                                              -sinFwd, chokedB,
                                              liquidRheology, sigmaOfP);
        forward = qF > qB;
        if (forward)
        {
            q = qF;
            choked = chokedF;
        }
        else
        {
            q = qB;
            choked = chokedB;
        }
    }
    const double pUp = forward ? pFrom : pTo;
    const double pLow = forward ? pTo : pFrom;
    const double sinFlow = forward ? sinFwd : -sinFwd;

    // 5.144 (G6): the PASSIVE thermal march rides the pressure profile the
    // hydraulics above already produced. Nothing declared -> null pointer ->
    // the isothermal march, bit for bit. The rate solve (g above) is NOT
    // handed the thermal: T never feeds dp, so threading it through the
    // secant would be pure overhead with no number to show for it.
    SteadyThermal th;
    const SteadyThermal* pth = nullptr;
    if (pipe.thermalOn())
    {
        th.temperatureIn = fluid.temperature;
        th.jouleThomson = pipe.jouleThomson;
        th.wallHeatTransfer = pipe.wallHeatTransfer;
        th.segmentU = pipe.wallHeatTransferSegments.empty()
                          ? nullptr
                          : &pipe.wallHeatTransferSegments;
        th.wallModel = pipe.wallModel;
        th.wallTemperature = pipe.wallTemperature;
        th.surfaceTemperature = pipe.surfaceTemperature;
        th.geothermalGradient = pipe.geothermalGradient;
        th.heatCapacity = pipe.heatCapacity;
        // 5.189: an EOS-derived coefficient needs a composition, and the
        // black-oil fluid does not carry one (the parser refuses this at
        // read time; the guard here keeps the invariant local to the march).
        if (pipe.jouleThomsonFromEos)
            cfFatal("steadyNetwork: black-oil march cannot price jouleThomson "
                    "from the EOS: no composition");
        pth = &th;
    }
    const std::vector<double> segSinVec =
        pipe.surveyOn ? surveySegSinOf(pipe, 24) : std::vector<double>();
    const BeggsBrillMarchOut m = marchBeggsBrillPipe(
        fluid, q, pUp, pipe.length, pipe.diameter, pipe.roughness, sinFlow,
        24, liquidRheology, pipe.flowPattern, sigmaOfP, pipe.correlation,
        pth, pipe.hbGriffith, pipe.surveyOn ? &segSinVec : nullptr);
    out.flowRate = forward ? q : -q;
    out.velocity = m.vmOutlet * (forward ? 1.0 : -1.0);
    out.holdup = m.holdupOutlet;
    out.regime = m.regimeMax;
    out.regimeOutlet = m.regimeOutlet;
    out.regimeFromMapOutlet = m.regimeMapOutlet;
    out.lambdaOutlet = m.lambdaOutlet;
    out.nfrOutlet = m.nfrOutlet;
    out.patternL1 = m.l1Outlet;
    out.patternL2 = m.l2Outlet;
    out.patternL3 = m.l3Outlet;
    out.patternL4 = m.l4Outlet;
    out.massRate = m.massRate;
    out.temperature = m.temperature; // 5.144 (G6): 0 unless thermal declared
    out.frictionFactor = m.fTpOutlet;
    out.choked = choked || !m.ok;
    out.dpFriction = dp; // gravity is inside the BB gradient
    // 5.158 / 5.159: the final march's per-segment readout profile, copied
    // after every hydraulic number above is set. The slug and risk reports
    // read it back; nothing on this function's answers does.
    out.marchSegPressure = m.segPressure;
    out.marchSegHoldup = m.segHoldup;
    out.marchSegLambda = m.segLambda;
    out.marchSegVm = m.segVm;
    out.marchSegRhoM = m.segRhoM;
    out.marchSegRegime = m.segRegime;
    out.marchSegTemperature = m.segTemperature;
    // 5.188 (G6b): the profile's own x grid and the coefficient it priced, so
    // the report can print "at this distance, this p and this T" without
    // re-deriving the march. Both stay empty / 0 on the isothermal path.
    out.marchSegX = m.segX;
    out.marchThermalU = pipe.wallHeatTransfer;
    out.marchOutletPressure = m.pOutlet;

    // conductance dQ/d(dpTotal). One-sided difference probing DEEPER into
    // the live branch: a central difference straddles the static-head knee
    // whenever the low side sits within h of it, and the resulting chord
    // points Newton back into the dead region. From the dead side (q == 0)
    // probe just past the knee so the step can cross it.
    const double h = std::max(1e-4 * std::abs(dp), 1.0);
    double c = 0.0;
    if (q > 0.0)
    {
        bool ch2 = false;
        const double qDeep = twoPhaseRateForDrop(
            pipe, fluid, pUp, std::max(pLow - h, 1.0), sinFlow, ch2,
            liquidRheology, sigmaOfP);
        c = (qDeep - q) / h;
    }
    if (!(c > 0.0))
    {
        // below the static head the law is flat at zero flow; probe just
        // past the knee so Newton still gets a usable slope pointing the
        // residual downhill
        const std::vector<double> segSinVec =
            pipe.surveyOn ? surveySegSinOf(pipe, 24) : std::vector<double>();
        const BeggsBrillMarchOut dead =
            marchBeggsBrillPipe(fluid, 0.0, pUp, pipe.length, pipe.diameter,
                                pipe.roughness, sinFlow, 24,
                                liquidRheology, pipe.flowPattern, sigmaOfP,
                                pipe.correlation, nullptr, pipe.hbGriffith,
                                pipe.surveyOn ? &segSinVec : nullptr);
        const double deficit = pLow - dead.pOutlet; // missing pressure
        if (deficit > 0.0)
        {
            bool ch2 = false;
            const double probeLow = std::max(dead.pOutlet - 0.01 * pUp, 1.0);
            const double qProbe = twoPhaseRateForDrop(
                pipe, fluid, pUp, probeLow, sinFlow, ch2,
                liquidRheology, sigmaOfP);
            c = qProbe / std::max(pLow - probeLow, 1.0);
        }
        if (!(c > 0.0)) c = 1e-6; // keep the network Jacobian nonsingular
    }
    out.conductance = c;
    return true;
}

// 5.158 (P18): the slug-START forecast. March segment i covers the arc
// interval [i ds, (i+1) ds) of the flow direction (ds = length / n, n the
// segments the final march actually ran). A DIP is a strict downhill ->
// uphill change between adjacent trajectory segments; the liquid the march's
// holdup profile carries down the upstream run is pooled on the uphill run
// downstream of the dip,
//   poolVolume = sum_i HL_i A overlap_i,  criterion = poolVolume / (A L_up),
// and the verdict is slug-start when the upstream run's own march segment is
// SEGREGATED (the stratified flow that actually drains to the dip; the
// transition band is not a statement that liquid reaches it) and the pooled
// level reaches the declared fill threshold. A pure report: the march never
// reads any of this back.
SlugReadout slugReadoutOf(const NetworkPipe& pipe, const PipeFlowResult& flow)
{
    SlugReadout r;
    r.pipeName = pipe.name;
    if (!pipe.slugDeclared || pipe.slugTrajectory.size() < 2) return r;
    const int n = static_cast<int>(flow.marchSegHoldup.size());
    if (n < 1 || !(pipe.length > 0.0) || !(pipe.diameter > 0.0)) return r;
    const double area = 0.25 * M_PI * pipe.diameter * pipe.diameter;
    const double ds = pipe.length / static_cast<double>(n);
    r.fillThreshold = pipe.slugFillThreshold;

    // cumulative arc positions of the trajectory
    std::vector<double> s0(pipe.slugTrajectory.size() + 1, 0.0);
    for (std::size_t k = 0; k < pipe.slugTrajectory.size(); ++k)
        s0[k + 1] = s0[k] + pipe.slugTrajectory[k].length;

    for (std::size_t k = 0; k + 1 < pipe.slugTrajectory.size(); ++k)
    {
        const double aUp = pipe.slugTrajectory[k].angleDeg;
        const double aDn = pipe.slugTrajectory[k + 1].angleDeg;
        if (!(aUp < 0.0 && aDn > 0.0)) continue; // strict - -> +
        SlugDipReadout d;
        d.s = s0[k + 1];
        // the march segment holding the upstream run's midpoint: its regime
        // says whether the liquid actually drains down to the dip
        const double smid = 0.5 * (s0[k] + s0[k + 1]);
        int im = static_cast<int>(smid / ds);
        if (im < 0) im = 0;
        if (im > n - 1) im = n - 1;
        d.regime = flow.marchSegRegime[static_cast<std::size_t>(im)];
        d.regimeOk = d.regime == 0;
        // pool the holdup profile over the uphill run's arc interval
        const double upLo = s0[k + 1];
        const double upHi = s0[k + 2];
        double vol = 0.0;
        for (int i = 0; i < n; ++i)
        {
            const double segLo = i * ds;
            const double segHi = (i + 1) * ds;
            const double ov = std::min(upHi, segHi) - std::max(upLo, segLo);
            if (ov > 0.0)
                vol += flow.marchSegHoldup[static_cast<std::size_t>(i)]
                     * area * ov;
        }
        d.volume = vol;
        d.criterion = vol / (area * pipe.slugTrajectory[k + 1].length);
        d.slugStart = d.regimeOk && d.criterion >= r.fillThreshold;
        r.dips.push_back(d);
    }
    return r;
}

// 5.159 (P22): the erosion / hydrate / corrosion readout, priced entirely
// from the march's own per-segment profile. Erosion is the API RP 14E form
// Ve = C / sqrt(rho_m): C is quoted in the published field-unit convention
// (ft/s at lb/ft3), so the SI mixture density is converted with the exact
// 1 lb/ft3 = 16.0184634 kg/m3 and the ft/s result with 0.3048. The hydrate
// line is the declared Clausius-Clapeyron dissociation temperature at each
// segment's pressure against the segment's temperature (the march's own when
// the pipe carries thermal physics, the fluid's isothermal one otherwise);
// the subcooling Teq - T is graded against the declared K bands. Corrosion
// prices the declared CO2 / H2S mole fractions at the HIGHEST march pressure
// against the declared Pa bands, sour at the published 0.5 psia = 3447 Pa
// H2S partial-pressure line. Every number is a report, not a charge.
RiskReadout riskReadoutOf(const NetworkPipe& pipe,
                          const pvt::BlackOilFluid& fluid,
                          const PipeFlowResult& flow)
{
    RiskReadout r;
    r.pipeName = pipe.name;
    if (!pipe.riskDeclared) return r;
    const int n = static_cast<int>(flow.marchSegPressure.size());
    if (n < 1) return r;

    // erosion: the worst vm / Ve segment
    r.erosionOn = true;
    {
        double best = -1.0;
        int over = 0;
        for (int i = 0; i < n; ++i)
        {
            const double rhoLb = flow.marchSegRhoM[static_cast<std::size_t>(i)]
                               / 16.0184634;
            const double velim = rhoLb > 0.0
                ? 0.3048 * pipe.erosionC / std::sqrt(rhoLb)
                : 0.0;
            const double ratio = velim > 0.0
                ? flow.marchSegVm[static_cast<std::size_t>(i)] / velim
                : 0.0;
            if (ratio > 1.0) ++over;
            if (ratio > best)
            {
                best = ratio;
                r.erosionSeg = i;
                r.erosionVm = flow.marchSegVm[static_cast<std::size_t>(i)];
                r.erosionVelim = velim;
                r.erosionRatio = ratio;
            }
        }
        r.erosionOverCount = over;
    }

    // hydrate: the worst segment -- the one with the LARGEST subcooling
    // Teq - T, since formation needs the stream BELOW the dissociation
    // temperature and the margin grows with Teq - T
    if (pipe.hydrateDeclared)
    {
        r.hydrateOn = true;
        double worst = 0.0;
        for (int i = 0; i < n; ++i)
        {
            const double p = flow.marchSegPressure[static_cast<std::size_t>(i)];
            const double teq = pvt::hydrateEquilibriumTemperature(
                p, pipe.hydrateTref, pipe.hydratePRef, pipe.hydrateDH);
            const double t = pipe.thermalOn()
                ? flow.marchSegTemperature[static_cast<std::size_t>(i)]
                : fluid.temperature;
            const double dts = teq - t;
            if (i == 0 || dts > worst)
            {
                worst = dts;
                r.hydrateSeg = i;
                r.hydrateP = p;
                r.hydrateT = t;
                r.hydrateTeq = teq;
                r.hydrateSubcool = dts;
            }
        }
        r.hydrateGrade = r.hydrateSubcool < 0.0
                             ? 0
                             : (r.hydrateSubcool < pipe.hydrateLowK
                                    ? 1
                                    : (r.hydrateSubcool < pipe.hydrateMediumK
                                           ? 2
                                           : 3));
    }

    // corrosion: partial pressures at the highest march pressure
    if (pipe.co2Fraction > 0.0 || pipe.h2sFraction > 0.0)
    {
        r.corrosionOn = true;
        double pMax = flow.marchSegPressure[0];
        for (int i = 1; i < n; ++i)
            pMax = std::max(pMax, flow.marchSegPressure[static_cast<std::size_t>(i)]);
        r.co2Partial = pipe.co2Fraction * pMax;
        r.h2sPartial = pipe.h2sFraction * pMax;
        r.co2Band = r.co2Partial < pipe.co2LowPa
                        ? 0
                        : (r.co2Partial < pipe.co2HighPa ? 1 : 2);
        r.sour = r.h2sPartial >= 3447.0; // 0.5 psia, the published sweet/sour line
    }
    return r;
}

// D10/D7 (§5.121): the opaque live-flash state owned by a PvtFluid. Kept out
// of the public header so adding the opt-in cannot perturb any include order;
// held by shared_ptr so the network stays copyable.
struct LiveFlashState
{
    pvt::Mixture mix;                  // feed, value semantics
    double T = 0;                      // K
    std::vector<double> parachors;     // SI [N^0.25 m^2.75/mol], feed order
    bool sigmaFromP = false;           // per-station Weinaug-Katz sigma
    // Declared TABULATED LBC critical volumes [m3/mol], feed order. Empty =
    // the historical PR-Zc fallback; a declaration replaces the back-calculation.
    std::vector<double> criticalVolumes;
};

namespace
{
// The PR-Zc fallback now lives in pvt::lbcCriticalVolumes (Viscosity.h);
// declared tabulated Vcs replace it without touching this call site.

/// Lohrenz-Bray-Clark viscosity of one flash phase from the PR molar volume
/// the flash already computed: no second EOS solve, so the checker can
/// re-derive every number from the reported phase state alone.
double lbcOfPhase(const pvt::Mixture& mix, double T, double molarVolume,
                  const std::vector<double>& declaredVcs)
{
    const std::vector<pvt::Component>& comps = mix.components();
    std::vector<double> MWs(comps.size()), Tcs(comps.size()), Pcs(comps.size());
    std::vector<double> Vcs = declaredVcs.empty()
        ? pvt::lbcCriticalVolumes(mix)
        : declaredVcs;
    if (Vcs.size() != comps.size())
        cfFatal("steadyNetwork: LBC criticalVolumes must give one value per component");
    for (std::size_t i = 0; i < comps.size(); ++i)
    {
        MWs[i] = comps[i].molarMass * 1.0e3; // kg/mol -> g/mol
        Tcs[i] = comps[i].Tc;
        Pcs[i] = comps[i].Pc;
    }
    return pvt::lbcViscosity(T, molarVolume, mix.fractions(), MWs, Tcs,
                             Pcs, Vcs);
}
} // namespace

void PvtFluid::evaluateAt(double p, double& rhoOut, double& muOut,
                          double& sigmaOut, double T) const
{
    sigmaOut = 0.0;
    if (!liveFlash)
    {
        rhoOut = pvtLookup(pressure, rho, p);
        muOut = pvtLookup(pressure, mu, p);
        return;
    }

    // 5.194 (G6d): the state temperature. T <= 0 is the feed's own temperature,
    // i.e. exactly the pre-feedback call every existing case makes.
    const LiveFlashState& lf = *liveFlash;
    const double tState = (T > 0.0) ? T : lf.T;
    const pvt::FlashResult r = pvt::tpFlash(lf.mix, tState, p);
    if (!(r.liquidMassDensity > 0.0) || !(r.vapourMassDensity > 0.0))
        cfFatal("pvtLiveFlash: non-positive flash phase density at p="
                + std::to_string(p));

    // Mass fractions carried by the bulk: beta is a MOLE vapour fraction,
    // so weight each phase by (moles x molar mass). In a single-phase flash
    // the two phase densities are equal, which collapses this to rho = the
    // single PR root independently of beta.
    const double beta = r.vaporFraction;
    const double mL = (1.0 - beta) * r.liquidMolarMass;
    const double mV = beta * r.vapourMolarMass;
    const double mtot = mL + mV;
    if (!(mtot > 0.0))
        cfFatal("pvtLiveFlash: zero phase mass at p=" + std::to_string(p));
    const double wL = mL / mtot;
    const double wV = mV / mtot;

    // Homogeneous (no-slip) mixture density: specific volume is mass-
    // fraction weighted. The march's mass rate is fixed by the upstream
    // value of exactly this quantity.
    rhoOut = 1.0 / (wL / r.liquidMassDensity + wV / r.vapourMassDensity);

    // Phase LBC viscosities on the flash's own PR molar volumes, combined
    // with the same mass fractions. Single phase: both phases are the feed
    // and both viscosities are identical.
    pvt::Mixture liq = lf.mix;
    liq.setFractions(r.liquidComposition);
    pvt::Mixture vap = lf.mix;
    vap.setFractions(r.vapourComposition);
    const double muL = lbcOfPhase(liq, tState, r.liquid.molarVolume,
                                  lf.criticalVolumes);
    const double muV = lbcOfPhase(vap, tState, r.vapour.molarVolume,
                                  lf.criticalVolumes);
    muOut = wL * muL + wV * muV;
    if (!(rhoOut > 0.0) || !(muOut > 0.0))
        cfFatal("pvtLiveFlash: non-positive rho/mu at p=" + std::to_string(p));

    // The convenience WK overload would flash again; build it from the phase
    // state already in hand (exactly the same expression as surfaceTension.cpp).
    if (lf.sigmaFromP && lf.parachors.size() == lf.mix.size()
        && !r.singlePhase)
    {
        const double Vml = r.liquidMolarMass / r.liquidMassDensity;
        const double Vmg = r.vapourMolarMass / r.vapourMassDensity;
        sigmaOut = pvt::weinaugKatzSurfaceTension(
            lf.parachors, Vml, Vmg, r.liquidComposition,
            r.vapourComposition);
    }
}

// 5.189: the EOS-derived Joule-Thomson coefficient. The composition is the
// one the case declared; the state is evaluated at the local (T, p). SRK
// carries its own analytic dZ/dT (srkDzDt). Peng-Robinson has no analytic
// derivative in this layer, so its branch takes a central difference of the
// EOS state - documented rather than hidden, and the two routes agree with
// the thermodynamic definition to what tests/unit/srkTests.cpp pins down.
double PvtFluid::eosJouleThomson(double T, double p, double cp, bool srk) const
{
    const LiveFlashState* st = composition ? composition.get() : liveFlash.get();
    if (!st || !(cp > 0.0) || !(T > 0.0) || !(p > 0.0)) return 0.0;
    if (srk) return pvt::srkJouleThomson(st->mix, T, p, cp);
    const double h = 1e-4 * T;
    const double dZ = (pvt::pengRobinsonState(st->mix, T + h, p).Z
                       - pvt::pengRobinsonState(st->mix, T - h, p).Z)
                      / (2.0 * h);
    // mass-specific c_p, so the molar volume's M appears (see Srk.h)
    return pvt::kGasConstant * T * T * dZ / (p * st->mix.molarMass() * cp);
}

double PvtFluid::density(double p) const
{
    double rho = 0, mu = 0, sigma = 0;
    evaluateAt(p, rho, mu, sigma);
    return rho;
}

double PvtFluid::viscosity(double p) const
{
    double rho = 0, mu = 0, sigma = 0;
    evaluateAt(p, rho, mu, sigma);
    return mu;
}

double PvtFluid::surfaceTension(double p) const
{
    double rho = 0, mu = 0, sigma = 0;
    evaluateAt(p, rho, mu, sigma);
    return sigma;
}

namespace
{

/// Segment-wise march of the local mechanical-energy balance
/// dp/ds = -rho(p) g sinFlow - f(rho v, mu(p)) rho v |v| / (2D)
///         + rho(p) g H(q)/L
/// from pUp toward the far end, with rho/mu evaluated at the segment-start
/// pressure and the valve K spread uniformly over the segments. The mass
/// rate is fixed by the UPSTREAM density (w = rho(pUp) q), so the local
/// velocity grows as the fill expands along the pipe; with a constant
/// table every segment reproduces the constant-property law exactly and
/// the march coincides with solvePipeFlow's implicit pipe law.
struct PvtMarchOut
{
    bool ok = false;
    double pOutlet = 0;
    double velocity = 0; // outlet velocity, flow direction
    double reynolds = 0; // outlet conditions, as in the BB march
    double friction = 0;
    double sigmaOutlet = 0; // D10/D7: WK sigma at the last station (0 otherwise)
    // 5.144 (G6): outlet temperature of the PASSIVE thermal march; 0 when
    // no SteadyThermal was handed, the isothermal march, bit for bit.
    double temperature = 0;
    // 5.188 (G6b) / 5.189: the same per-segment profile the black-oil march
    // publishes, so the thermal report reads the same shape on both paths.
    // Filled on every march; copied to the result only when thermal is on.
    std::vector<double> segX, segPressure, segTemperature;
    // 5.226 (H3): the QMOM march's statement -- the population at the outlet
    // (after every segment's residence-time advance), the Sauter diameter
    // there, and the clean / Brinkman viscosities the LAST segment charged.
    // Inert (false / zeros) on every undeclared pipe.
    bool qmomOn = false;
    double qmomM0 = 0, qmomM1 = 0, qmomM2 = 0, qmomM3 = 0;
    double qmomD32 = 0;
    double qmomMuBase = 0;
    double qmomMuEff = 0;
};

/// 5.226 (H3): the constant-kernel Smoluchowski moment source of §5.107,
/// restated here from src/3d/src/dpm/qmom.cpp `coagulationSource` (cross-layer
/// include of the 3D fvm header from the 1D kernel is not available):
///   dM0/dt = -beta M0^2 / 2   (number density lost to pairs)
///   dM1/dt = 0                (total droplet volume is conserved)
///   dM2/dt = beta M1^2        (area moment gains)
///   dM3/dt = 3 beta M1 M2     (volume-weighted moment gains)
/// The checker re-derives these four rates independently.
void qmomCoagulationSource(const double m[4], double beta, double dm[4])
{
    dm[0] = -0.5 * beta * m[0] * m[0];
    dm[1] = 0.0;
    dm[2] = beta * m[1] * m[1];
    dm[3] = 3.0 * beta * m[1] * m[2];
}

/// pStart (5.156, P10): the pressure the march BEGINS at when a compressor
/// sits at the upstream end of the line -- pUp * PR(|q|), decided by the
/// caller inside the rate solve so the curve rides with the secant. The mass
/// rate stays booked at the REAL pUp density (the machine takes what the node
/// holds), which is the declared simplification of this increment. 0 = no
/// machine: the march starts at pUp and every existing number is bit for bit.
/// 5.194 (G6d): the density the march books its mass rate at -- the ENTRY
/// state, not the table's. Without the switch this is `fluid.density(pUp)`,
/// the pre-feedback path bit for bit. With it the entry state is (pUp, T_in),
/// because that is the state the stream actually arrives in and the mass it
/// carries is fixed there: a feedback march whose inlet mass used the table's
/// single-temperature density would book a mass the pipe never sees.
double pvtEntryDensity(const NetworkPipe& pipe, const PvtFluid& fluid,
                       double pUp, const SteadyThermal* thermal)
{
    if (!pipe.pvtTemperatureFeedback) return fluid.density(pUp);
    double rhoIn = 0.0, muIn = 0.0, sgIn = 0.0;
    fluid.evaluateAt(pUp, rhoIn, muIn, sgIn,
                     thermal ? thermal->temperatureIn : 0.0);
    return rhoIn;
}

PvtMarchOut marchPvtPipe(const NetworkPipe& pipe, const PvtFluid& fluid,
                         double qIn, double pUp, double sinFlow, int nSeg,
                         std::vector<PvtFlashStation>* stations = nullptr,
                         const SteadyThermal* thermal = nullptr,
                         double pStart = 0.0)
{
    if (pipeIsAnnulus(pipe))
        cfFatal("pvtMarch: pipe '" + pipe.name + "' declares an annulus: the "
                "PVT march is a circular-bore correlation chain");    if (pipeIsRect(pipe))
        cfFatal("pvtMarch: pipe '" + pipe.name + "' declares a rectangular "
                "duct: the PVT march is a circular-bore correlation chain");
    const double area = 0.25 * M_PI * pipe.diameter * pipe.diameter;
    // 5.194 (G6d): the conserved mass rate is booked at the ENTRY state; the
    // feedback switch is what makes that state (pUp, T_in) instead of the
    // table's single-temperature density.
    const double rhoEntry = pvtEntryDensity(pipe, fluid, pUp, thermal);
    const double massRate = rhoEntry * qIn;
    const double ds = pipe.length / nSeg;
    // 5.193 (G9): the march reads the pipe's own survey when it has one, so a
    // build section and a lateral each pay their own gravity. Without a
    // survey this is the pipe's single slope for every segment -- the
    // pre-survey path bit for bit.
    const std::vector<double> segSinVec =
        pipe.surveyOn ? surveySegSinOf(pipe, nSeg) : std::vector<double>();
    const std::vector<double>* segSin = pipe.surveyOn ? &segSinVec : nullptr;
    const double headPerM = machineHeadOf(pipe, qIn, rhoEntry) / pipe.length;
    const double pFloor = std::max(0.02 * pUp, 1000.0);
    if (thermal && thermal->segmentU
        && static_cast<int>(thermal->segmentU->size()) != nSeg)
        cfFatal("pvtMarch: SteadyThermal segmentU needs exactly one entry "
                "per march segment");
    PvtMarchOut m;
    double p = pStart > 0.0 ? pStart : pUp;
    double t = thermal ? thermal->temperatureIn : 0.0;
    // 5.226 (H3): the population the march carries, in the pipe's own frame.
    // Undeclared, qm is never read and every number below is bit for bit.
    double qm[4] = {pipe.qmomM0, pipe.qmomM1, pipe.qmomM2, pipe.qmomM3};
    if (pipe.qmomOn) m.qmomOn = true;
    for (int i = 0; i < nSeg; ++i)
    {
        // D10/D7: one flash per station when liveFlash is set. On the table
        // path evaluateAt returns exactly pvtLookup(pressure,rho/mu,p), so
        // this is bit-for-bit the original march.
        double rho = 0, mu0 = 0, sigma = 0;
        // 5.194 (G6d): with the switch the properties are evaluated at the
        // temperature the march has REACHED (the segment's own inlet T, from
        // the previous step's advance), which is what makes this a coupled
        // march. Without it T reads 0 and the feed's snapshot is used, bit for
        // bit. The stage temperature is the one the previous segment left, so
        // the coupling is the same first-order Euler the rest of the march is.
        const double tState = pipe.pvtTemperatureFeedback ? t : 0.0;
        fluid.evaluateAt(p, rho, mu0, sigma, tState);
        const double v = massRate / (rho * area);
        // 5.226 (H3): the population first takes this segment's residence
        // time, then the segment charges the viscosity it has reached. No
        // flow means no residence time: the moments do not evolve and the
        // pipe charges the declared state's viscosity.
        if (pipe.qmomOn)
        {
            const double av = std::abs(v);
            if (av > 0.0)
            {
                double dm[4];
                qmomCoagulationSource(qm, pipe.qmomBeta, dm);
                const double dt = ds / av;
                for (int k = 0; k < 4; ++k) qm[k] += dt * dm[k];
            }
        }
        // 5.76 (A5 increment 2): the PVT table supplies the pressure-
        // dependent Newtonian baseline; a declared rheology reprices it on
        // this segment's wall rate. No declaration returns the table value
        // bit for bit (the whole PVT suite is that witness).
        double mu = wallRateApparentViscosityPaS(
            Rheology{pipe.flowIndex, pipe.consistency, pipe.yieldStress},
            8.0 * std::abs(v) / pipe.diameter, mu0);
        // 5.226 (H3): the Sauter feedback. Past the declared diameter the
        // emulsion is Brinkman's; below it the clean table viscosity stands.
        if (pipe.qmomOn)
        {
            const double d32 = qm[2] > 0.0 ? qm[3] / qm[2] : 0.0;
            const double phi = pipe.qmomPhiFromM0 > 0.0
                                 ? pipe.qmomPhiFromM0 * qm[0]
                                 : pipe.qmomPhi;
            m.qmomMuBase = mu;
            if (d32 > pipe.qmomDThreshold) mu *= (1.0 + 2.5 * phi);
            m.qmomMuEff = mu;
            m.qmomM0 = qm[0]; m.qmomM1 = qm[1];
            m.qmomM2 = qm[2]; m.qmomM3 = qm[3];
            m.qmomD32 = d32;
        }
        const double re = rho * v * pipe.diameter / mu;
        const double f = churchillFriction(re, pipe.roughness / pipe.diameter);
        m.velocity = v;
        m.reynolds = re;
        m.friction = f;
        m.sigmaOutlet = sigma;
        if (stations && fluid.liveFlash)
            stations->push_back(PvtFlashStation{p, rho, mu, sigma});
        const double sinI = segSin ? (*segSin)[static_cast<std::size_t>(i)]
                                   : sinFlow;
        const double dp = rho * kGravity * sinI * ds
            + (f * ds / pipe.diameter + pipe.valveK / nSeg) * rho * v * v / 2.0
            - rho * kGravity * headPerM * ds;
        // Acceleration: with the mass rate booked at the entry density, an
        // expanding fill accelerates and the momentum flux rho v^2 grows.
        // The segment pays that increase; without it a gas that expands a
        // lot under-prices the drop. rho_out is solved from the same table
        // state the march already uses, so an incompressible fill (rho_out
        // == rho) reproduces the old march bit for bit.
        double pOut = p - dp;
        for (int it = 0; it < 4; ++it)
        {
            // Floor guard before the outlet-state flash: the old march priced
            // properties at the ENTRY pressure only and never re-queried the
            // segment tail, so a high-rate bisection probe that crossed the
            // floor downstream could not feed an unsound pressure to the PVT
            // layer. The acceleration fixed point re-flashes the tail instead,
            // so a tail that has reached (or slipped past) the floor must not
            // be handed to tpFlash: break out of the fixed point and leave the
            // clipped value to the existing hitFloor path below, the same
            // treatment the old march gave a floor-crossing segment.
            if (!(pOut > pFloor)) break;
            double rhoOut = 0, muOut = 0, sigmaOut = 0;
            fluid.evaluateAt(pOut, rhoOut, muOut, sigmaOut, tState);
            if (!(rhoOut > 0.0)) break;
            const double accelDrop = massRate * massRate / (area * area)
                                   * (1.0 / rhoOut - 1.0 / rho);
            const double pTry = p - dp - accelDrop;
            if (std::abs(pTry - pOut)
                <= 1e-12 * std::max(1.0, std::abs(pOut)))
            {
                pOut = pTry;
                break;
            }
            pOut = pTry;
        }
        m.segX.push_back(i * ds);
        m.segPressure.push_back(p);
        m.segTemperature.push_back(t);
        const double pOld = p;
        p = pOut;
        const bool hitFloor = !(p > pFloor);
        if (hitFloor) p = pFloor;
        // 5.144 (G6): the passive thermal step, shared with the black-oil
        // march so the two paths price the same law.
        if (thermal)
            t = steadyThermalAdvance(*thermal, t, pOld, p, i, ds, sinI,
                                     massRate);
        if (hitFloor)
        {
            m.pOutlet = pFloor;
            m.temperature = t;
            return m; // ok stays false: the fill cannot carry this rate
        }
    }
    m.ok = true;
    m.pOutlet = p;
    m.temperature = t;
    return m;
}

/// Rate landing the march on pLow via Illinois. The loss grows without
/// bound in |q| (friction + the pump curve falling) while the static head
/// stays bounded, so the law is monotone with a single root: a plain
/// expanding bracket finds it, unlike the U-shaped two-phase gradient.
double pvtRateForDrop(const NetworkPipe& pipe, const PvtFluid& fluid,
                      double pUp, double pLow, double sinFlow, bool& choked,
                      bool boostUp = false,
                      const SteadyThermal* thermal = nullptr)
{
    choked = false;
    const double pFloor = std::max(0.02 * pUp, 1000.0);
    // 5.156 (P10): with a compressor at the `from` end the machine's boost
    // rides inside the rate solve, where the curve can read the trial rate.
    // Forward flow (boostUp) starts the march at pUp * PR(q); reverse flow
    // lands on pLow * PR(q) -- the same statement, pFrom * PR, seen from the
    // end the flow enters. A bare pipe reads PR == 1 and both ends stay the
    // numbers they were, bit for bit.
    auto g = [&](double q)
    {
        const double pr = compressorRatioOf(pipe, q);
        const double start = boostUp ? pUp * pr : 0.0;
        const double target = boostUp ? pLow : pLow * pr;
        const PvtMarchOut m =
            marchPvtPipe(pipe, fluid, q, pUp, sinFlow, 24, nullptr, thermal,
                         start);
        // floor-crossing marches return the floor shortfall: continuous
        // with the unchoked branch and correctly scaled for the secant
        return m.ok ? m.pOutlet - target : pFloor - target;
    };

    const double qSeed = std::max(1e-3 * 0.25 * M_PI * pipe.diameter
                                      * pipe.diameter, 1e-8);
    const double g0 = g(0.0);
    if (g0 <= 0.0) return 0.0; // the static column alone reaches the low side

    double qLo = 0.0, fLo = g0;
    double qHi = qSeed, fHi = g(qHi);
    int expansions = 0;
    while (fHi > 0.0)
    {
        qLo = qHi;
        fLo = fHi;
        qHi *= 4.0;
        fHi = g(qHi);
        if (++expansions > 60)
        {
            choked = true; // no rate within reach lands on the low side
            return qHi;
        }
    }

    double q = qLo + (qHi - qLo) * fLo / (fLo - fHi);
    for (int it = 0; it < 100; ++it)
    {
        const double fq = g(q);
        if (std::abs(fq) < 1e-12 * pUp) break;
        if (fq > 0.0)
        {
            qLo = q;
            fLo = fq;
            // Illinois: bias the retained endpoint against re-selection
            fLo *= 0.5;
        }
        else
        {
            qHi = q;
            fHi = fq;
            fHi *= 0.5;
        }
        const double qNew = qLo + (qHi - qLo) * fLo / (fLo - fHi);
        if (!(qNew > qLo && qNew < qHi)) q = 0.5 * (qLo + qHi);
        else q = qNew;
        if (std::abs(qHi - qLo) < 1e-12 * std::max(qHi, 1e-30)) break;
    }
    return q;
}

/// Gas-lift valve (5.89, D1 increment 3): the injection string's endpoint. Two
/// published closures run in order and neither borrows from the other. The
/// bellows force balance reads the injection pressure `pFrom` against the
/// tubing pressure `pTo` and says whether the port is open; the isentropic
/// nozzle then reads the gas at the state it actually sits in -- the density
/// the PVT table gives at pFrom, which is why this element cannot live on a
/// constant-property fill. The rate is reported at that same upstream density,
/// so `massRate / flowRate` is rho(pFrom) by construction and the node balance
/// downstream carries the volume the tubing sees.
///
/// The open/closed switch is discontinuous, and a steady solve has no history
/// to resolve it with: below Pvc the valve is shut, at or above Pvo it is open,
/// and in between it is reported as sitting in its own spread and priced shut.
/// Pricing it *open* there would invent flow the manufacturer's force balance
/// never authorized; the disposition token is what makes the band visible on
/// the readout instead of a silent zero.
void solveLiftValvePvt(const NetworkPipe& pipe, const PvtFluid& fluid,
                       double pFrom, double pTo, PipeFlowResult& out)
{
    const double area = 0.25 * M_PI * pipe.diameter * pipe.diameter;
    const GasLiftValve valve{pipe.liftDomePressure, pipe.liftAreaRatio,
                             pipe.liftSpringPressure, area,
                             pipe.liftDischargeCoeff, pipe.liftKappa};

    const GasLiftValveThresholds th = gasLiftValveThresholds(valve, pTo);
    out.valveOpeningPressure = th.openingPressure;
    out.valveClosingPressure = th.closingPressure;
    out.valveSpread = th.spread;
    out.valveProductionEffect = th.productionEffect;
    out.valveTestRackOpening = th.testRackOpening;
    const bool opens = pFrom >= th.openingPressure;
    out.valveDisposition = opens ? 1
                       : (pFrom <= th.closingPressure ? 0 : 2);

    const double rhoUp = fluid.density(pFrom);
    double q = 0;
    if (opens)
    {
        const GasLiftValveFlow fl = gasLiftValveFlow(valve, pFrom, rhoUp, pTo);
        out.valveChoked = fl.choked;
        out.choked = fl.choked;
        out.massRate = fl.massRate;
        q = fl.massRate / rhoUp;
    }
    out.flowRate = q;
    // The bore velocity referred to the density the volume rate is priced at,
    // which is the manifold's, not the throat's: a choked port runs sonic at a
    // throat density well below this one, and this row does not pretend to
    // march to it.
    out.velocity = q / area;
    out.reynolds = fluid.mu.empty()
                       ? 0.0
                       : rhoUp * std::abs(out.velocity) * pipe.diameter
                       / fluid.viscosity(pFrom);
    out.frictionFactor = 0; // no wall law is priced here
    // The whole drop is the port's: a point element has no wall to rub against,
    // so dpTotal = dpValve and dpFriction stays the zero it was born as. Read
    // as the march reads it, the drop adds up.
    out.dpValve = out.dpTotal;

    // dQ/d(dpTotal) at fixed injection pressure, the same convention as the
    // marched pipe law. A choked port is deliberately flat: it has stopped
    // listening to the downstream, and for a gas whose density rides with its
    // pressure it has stopped listening upstream in VOLUME terms too (its mass
    // rate rides with the manifold, its volume rate does not). A shut port
    // listens to nothing. Both keep a floor so a manifold of nothing but shut
    // ports still has a Jacobian, but the floor is not the pipe law's 1e-6:
    // the marched lines here slope at 1e-7, and a dead valve modelled as
    // sixty times stiffer than every live pipe in the network is a lie Newton
    // will step by.
    const double h = std::max(1e-4 * std::abs(out.dpTotal), 1.0);
    double c = 0;
    if (q > 0)
    {
        const double pHighProbe = pTo + h;
        const double pLowProbe = std::max(pTo - h, 1.0);
        const GasLiftValveFlow fp =
            gasLiftValveFlow(valve, pFrom, rhoUp, pLowProbe);
        const GasLiftValveFlow fm =
            gasLiftValveFlow(valve, pFrom, rhoUp, pHighProbe);
        c = (fp.massRate - fm.massRate) / rhoUp / (pHighProbe - pLowProbe);
    }
    if (!(c > 0)) c = 1e-12; // keep the network Jacobian nonsingular
    out.conductance = c;
}

} // namespace

bool solvePipeFlowPvt(const NetworkPipe& pipe, const PvtFluid& fluid,
                      double pFrom, double pTo, PipeFlowResult& out)
{
    out = PipeFlowResult{};
    out.name = pipe.name;
    out.regime = -1;
    const double dp = pFrom - pTo;
    out.dpTotal = dp;
    if (pipeIsLiftValve(pipe))
    {
        solveLiftValvePvt(pipe, fluid, pFrom, pTo, out);
        return true;
    }
    if (!(pipe.length > 0) || !(pipe.diameter > 0)
        || fluid.pressure.size() < 2)
        return true;

    // 5.156 (P10): a compressor at the `from` end still biases the direction
    // by its shutoff boost, but a bare pipe (and a machine that cannot lift
    // the pressure side) no longer reads the direction off the pressure
    // difference alone: both directions are priced and the larger feasible
    // rate carries the line, so a downhill gravity-driven pipe is expressible.
    const double pr0 = compressorRatioOf(pipe, 0.0);
    const double sinFwd = pipe.elevation / pipe.length;
    const double area = 0.25 * M_PI * pipe.diameter * pipe.diameter;

    // 5.144 (G6): the steady temperature march. Its T is PASSIVE unless
    // 5.194's feedback switch is on -- and when it is, the RATE SOLVE reads it
    // too (rho/mu at the reached T move dp), so the struct has to be built
    // before the rate solve instead of after it. Building it earlier cannot
    // move an answer by itself: T reaches dp only through the feedback switch,
    // and the unswitched march never reads `thermal` for a property.
    SteadyThermal th;
    const SteadyThermal* pth = nullptr;
    std::function<double(double, double)> jtFn; // 5.189: outlives the march
    if (pipe.thermalOn())
    {
        th.temperatureIn = fluid.temperature;
        th.jouleThomson = pipe.jouleThomson;
        th.wallHeatTransfer = pipe.wallHeatTransfer;
        th.segmentU = pipe.wallHeatTransferSegments.empty()
                          ? nullptr
                          : &pipe.wallHeatTransferSegments;
        th.wallModel = pipe.wallModel;
        th.wallTemperature = pipe.wallTemperature;
        th.surfaceTemperature = pipe.surfaceTemperature;
        th.geothermalGradient = pipe.geothermalGradient;
        th.heatCapacity = pipe.heatCapacity;
        // 5.189: the EOS route. jtFn must outlive the march call, so it is
        // declared in this scope and only its address travels in the struct.
        if (pipe.jouleThomsonFromEos)
        {
            const bool srk = pipe.jouleThomsonEosIsSrk;
            const double cp = pipe.heatCapacity;
            jtFn = [&fluid, cp, srk](double T, double p)
            { return fluid.eosJouleThomson(T, p, cp, srk); };
            th.jtOfState = &jtFn;
        }
        pth = &th;
    }
    bool forward = false;
    bool boostUp = false;
    bool choked = false;
    double q = 0.0;
    if (pFrom * pr0 >= pTo)
    {
        forward = true;
        boostUp = pipe.compressorRatio > 0.0;
        q = pvtRateForDrop(pipe, fluid, pFrom, pTo, sinFwd, choked,
                           boostUp, pth);
    }
    else
    {
        bool chokedF = false, chokedB = false;
        const bool boostF = pipe.compressorRatio > 0.0;
        const double qF = pvtRateForDrop(pipe, fluid, pFrom, pTo, sinFwd,
                                         chokedF, boostF, pth);
        const double qB = pvtRateForDrop(pipe, fluid, pTo, pFrom, -sinFwd,
                                         chokedB, false, pth);
        forward = qF > qB;
        if (forward)
        {
            q = qF;
            choked = chokedF;
            boostUp = boostF;
        }
        else
        {
            q = qB;
            choked = chokedB;
        }
    }
    const double pUp = forward ? pFrom : pTo;
    const double pLow = forward ? pTo : pFrom;
    const double sinFlow = forward ? sinFwd : -sinFwd;
    // 5.194 (G6d): the entry state the march booked its mass rate at, and that
    // mass itself -- printed so a checker can close rho v A = m from the case
    // text and the station readouts alone.
    const double rhoEntry = pvtEntryDensity(pipe, fluid, pUp, pth);
    out.marchTemperatureFeedback = pipe.pvtTemperatureFeedback;
    out.marchEntryDensity = pipe.pvtTemperatureFeedback ? rhoEntry : 0.0;
    out.marchMassRate = pipe.pvtTemperatureFeedback ? rhoEntry * q : 0.0;

    // D10/D7: the accepted march is replayed once with station capture so a
    // checker can read the per-station rho/mu/sigma the march actually used.
    // 5.156: with a compressor at the entry end the replay starts at the
    // boosted pressure, the same start the rate solve landed on.
    const PvtMarchOut m = marchPvtPipe(
        pipe, fluid, q, pUp, sinFlow, 24,
        fluid.liveFlash ? &out.pvtFlashStations : nullptr, pth,
        boostUp ? pUp * compressorRatioOf(pipe, q) : 0.0);
    out.sigmaOutlet = m.sigmaOutlet;
    out.temperature = m.temperature;
    // 5.226 (H3): the QMOM march's outlet statement, from the replayed march
    // the readouts above already agree with.
    if (pipe.qmomOn && m.qmomOn)
    {
        out.qmomState = 0;
        out.qmomM0 = m.qmomM0;
        out.qmomM1 = m.qmomM1;
        out.qmomM2 = m.qmomM2;
        out.qmomM3 = m.qmomM3;
        out.qmomD32 = m.qmomD32;
        out.qmomMuBase = m.qmomMuBase;
        out.qmomMuEff = m.qmomMuEff;
    }
    if (pipe.thermalOn())
    {
        // 5.189: the pvt march now carries the same profile the black-oil
        // march does, so the G6b report is path-independent.
        out.marchSegX = m.segX;
        out.marchSegPressure = m.segPressure;
        out.marchSegTemperature = m.segTemperature;
        out.marchThermalU = pipe.wallHeatTransfer;
        out.marchOutletPressure = m.pOutlet;
    }
    // 5.189: the coefficient the march actually used at its entry state, and
    // which route priced it - the number a checker re-derives from the case.
    if (pipe.thermalOn())
    {
        out.marchJtMode = pipe.jouleThomsonFromEos
                              ? (pipe.jouleThomsonEosIsSrk ? 2 : 1)
                              : 0;
        out.marchJouleThomson = pipe.jouleThomsonFromEos
            ? fluid.eosJouleThomson(fluid.temperature, pUp,
                                    pipe.heatCapacity,
                                    pipe.jouleThomsonEosIsSrk)
            : pipe.jouleThomson;
    }
    out.flowRate = (forward ? 1.0 : -1.0) * q;
    out.velocity = (forward ? 1.0 : -1.0) * m.velocity;
    out.reynolds = m.reynolds;
    out.frictionFactor = m.friction;
    out.choked = choked || !m.ok;
    out.dpFriction = dp; // gravity and the property variation are in the march
    // 5.194: with the feedback switch the row's own entry density is the one at
    // (pUp, T_in), the same state the march booked its mass at, so a machine
    // cannot be priced at a density the pipe never sees. Without the switch
    // this is fluid.density(pUp), bit for bit.
    const double rhoRow = pipe.pvtTemperatureFeedback ? rhoEntry
                                                      : fluid.density(pUp);
    out.dpValve = pipe.valveK * rhoRow * std::abs(out.velocity)
                * out.velocity / 2.0;
    out.dpPump = -rhoRow * kGravity * pumpHeadOf(pipe, q);
    // 5.191: the PCP's own account, kept as its own readout so a checker
    // reads the machine that priced it rather than a sum of two curves.
    out.dpPcp = -rhoRow * kGravity * pcpHeadOf(pipe, q, rhoRow);
    if (pipe.jetOn)
    {
        const JetRow jr = jetRowOf(pipe, q, rhoRow);
        out.dpJet = -rhoRow * kGravity * jr.head;
        out.jetFlowRatio = jr.ratio;
        out.jetState = jr.state;
    }

    // 5.156 (P10): the compressor's own statement, read from the hydraulic
    // state it was priced with -- never fed back into it (§5.144's passive
    // convention). The head is the declared enthalpy integral at the curve's
    // ratio, the mass rate is the one the march booked at the REAL entry-end
    // density, and T1 is the stream's declared table temperature. Zeros /
    // ratio 1 on every pipe without a machine, so no existing readout moves.
    if (pipe.compressorRatio > 0.0)
    {
        const double pr = compressorRatioOf(pipe, q);
        const double ex =
            (pipe.compressorKappa - 1.0) / pipe.compressorKappa;
        const double head = pipe.compressorZ * pipe.compressorGasConstant
                          * fluid.temperature / ex
                          * (std::pow(pr, ex) - 1.0);
        out.compressorRatio = pr;
        out.compressorBoost = pFrom * (pr - 1.0);
        out.compressorPower = fluid.density(pUp) * q * head
                            / pipe.compressorEfficiency;
        out.compressorTemperatureOut =
            fluid.temperature * std::pow(pr, ex / pipe.compressorEfficiency);
    }

    // conductance dQ/d(dpTotal). One-sided difference probing DEEPER into
    // the live branch: below the static head the law is flat at zero flow
    // (same knee as the two-phase law), and a central difference straddling
    // the knee points Newton back into the dead region. From the dead side
    // (q == 0) probe just past the knee so the step can cross it.
    const double h = std::max(1e-4 * std::abs(dp), 1.0);
    double c = 0.0;
    if (q > 0.0)
    {
        bool ch2 = false;
        const double qDeep = pvtRateForDrop(
            pipe, fluid, pUp, std::max(pLow - h, 1.0), sinFlow, ch2, boostUp,
            pth);
        c = (qDeep - q) / h;
    }
    if (!(c > 0.0))
    {
        // 5.156: the dead march and the knee probe carry the machine's
        // zero-rate boost so the Jacobian reads the same law the rate solve
        // solved. A bare pipe's boost is 1.0 and this is the old code.
        const PvtMarchOut dead = marchPvtPipe(
            pipe, fluid, 0.0, pUp, sinFlow, 24, nullptr, pth,
            boostUp ? pUp * pr0 : 0.0);
        const double deficit = (boostUp ? pLow : pLow * pr0) - dead.pOutlet;
        if (deficit > 0.0)
        {
            bool ch2 = false;
            const double probeLow = std::max(dead.pOutlet - 0.01 * pUp, 1.0);
            const double qProbe = pvtRateForDrop(
                pipe, fluid, pUp, probeLow, sinFlow, ch2, boostUp, pth);
            c = qProbe / std::max(pLow - probeLow, 1.0);
        }
        if (!(c > 0.0)) c = 1e-6; // keep the network Jacobian nonsingular
    }
    out.conductance = c;
    return true;
}

OneDNetwork loadNetworkFromDict(const Dictionary& caseDict)
{
    OneDNetwork net;

    if (caseDict.has("fluid"))
    {
        const Dictionary& fd = caseDict.getDict("fluid");
        if (fd.has("blackoil") && fd.has("pvt"))
            cfFatal("steadyNetwork: fluid cannot be both blackoil and pvt");
        if (fd.has("blackoil"))
        {
            // two-phase mode: sources are stock-tank oil rates, per-pipe
            // hydraulics come from Beggs-Brill + black-oil PVT
            const Dictionary& bd = fd.getDict("blackoil");
            net.hasBlackOil = true;
            net.blackOil.gasGravity = bd.getScalarOr("gasGravity", net.blackOil.gasGravity);
            net.blackOil.oilDensityStd = bd.getScalarOr("oilDensity", net.blackOil.oilDensityStd);
            net.blackOil.gorStd = bd.getScalarOr("gor", net.blackOil.gorStd);
            net.blackOil.waterCut = bd.getScalarOr("waterCut", net.blackOil.waterCut);
            net.blackOil.temperature = bd.getScalarOr("temperature", net.blackOil.temperature);
            net.blackOil.surfaceTension = bd.getScalarOr("surfaceTension", net.blackOil.surfaceTension);
            if (!(net.blackOil.oilDensityStd > 0))
                cfFatal("steadyNetwork: blackoil needs positive oilDensity");
        }
        else if (fd.has("pvt"))
        {
            // tabulated-property single-phase mode: rho(p)/mu(p) come from
            // the case table (piecewise-linear, clamped at the ends)
            const Dictionary& pd = fd.getDict("pvt");
            net.hasPvtFluid = true;
            net.pvtFluid.pressure = pd.getScalarList("pressure");
            // 5.144 (G6): the stream's inlet temperature for the PASSIVE
            // thermal march. Absent, the table path never reads it and
            // every existing pvt number is bit for bit; a pipe with
            // thermal physics on a pvt fluid without it is refused in
            // addPipe.
            net.pvtFluid.temperature = pd.getScalarOr("temperature", 0.0);
            if (net.pvtFluid.temperature < 0.0)
                cfFatal("steadyNetwork: fluid pvt temperature is kelvin");
            if (pd.has("compositional"))
            {
                // 5.110 (compositional properties into 1D, increment A):
                // opt-in path. rho(p)/mu(p) are GENERATED from a declared
                // composition with PR76 phase density + LBC residual
                // viscosity, evaluated at each pressure-grid node. The
                // march (solvePipeFlowPvt / marchPvtPipe) is untouched: it
                // still reads the same PvtFluid table. Without the opt-in
                // block the code below is byte-for-byte the original table
                // path.
                if (pd.has("rho") || pd.has("mu"))
                    cfFatal("steadyNetwork: fluid pvt cannot declare rho/mu"
                            " together with compositional");
                const Dictionary& comp = pd.getDict("compositional");
                if (!comp.has("composition"))
                    cfFatal("steadyNetwork: fluid pvt compositional needs a composition block");
                const pvt::Mixture mix = pvt::loadComposition(comp.getDict("composition"));
                // 5.189: the feed is published even on the table path, so a
                // pipe can price its Joule-Thomson coefficient from the EOS.
                {
                    std::shared_ptr<LiveFlashState> feed =
                        std::make_shared<LiveFlashState>();
                    feed->mix = mix;
                    feed->T = comp.getScalarOr("temperature", 0.0);
                    net.pvtFluid.composition = feed;
                }
                if (!comp.has("temperature"))
                    cfFatal("steadyNetwork: fluid pvt compositional needs temperature");
                const double T = comp.getScalar("temperature");
                if (!(T > 0.0))
                    cfFatal("steadyNetwork: compositional temperature must be positive");
                // 5.144 (G6): the flash temperature and the thermal march's
                // inlet temperature are the same physical number; two
                // declarations that disagree would price two streams.
                if (pd.has("temperature") && T != net.pvtFluid.temperature)
                    cfFatal("steadyNetwork: fluid pvt declares two different "
                            "temperatures: the compositional flash temperature "
                            "and the thermal inlet temperature are one number");
                const std::string phase = comp.getWordOr("phase", "stable");
                pvt::RootPreference pref = pvt::RootPreference::LowestGibbs;
                if (phase == "liquid") pref = pvt::RootPreference::Liquid;
                else if (phase == "vapour" || phase == "vapor") pref = pvt::RootPreference::Vapor;
                else if (phase != "stable")
                    cfFatal("steadyNetwork: unknown compositional phase '" + phase + "'");

                // D10/D7 (§5.121): opt-in per-station true flash. Default off,
                // in which case the 5.110 precomputed-root table below (and
                // every existing number) is byte-for-byte untouched.
                const bool liveFlash = comp.has("liveFlash")
                                    && comp.getBool("liveFlash");
                const bool sigmaFromP = comp.has("sigmaFromP")
                                     && comp.getBool("sigmaFromP");
                // D7 (§5.93) boundary repair: a declared criticalVolumes list
                // (m3/mol, feed order) replaces LBC's historical PR-Zc
                // back-calculation. Absent it the exact old path is kept.
                std::vector<double> criticalVolumes;
                if (comp.has("criticalVolumes"))
                {
                    criticalVolumes = comp.getScalarList("criticalVolumes");
                    if (criticalVolumes.size() != mix.size())
                        cfFatal("steadyNetwork: compositional criticalVolumes"
                                " must give one value per component");
                    for (const double vc : criticalVolumes)
                        if (!(vc > 0.0))
                            cfFatal("steadyNetwork: compositional criticalVolumes"
                                    " must be positive");
                }
                std::vector<double> parachors;
                if (sigmaFromP)
                {
                    if (!comp.has("parachors"))
                        cfFatal("steadyNetwork: compositional sigmaFromP needs parachors");
                    const Dictionary::Node& pn = comp.find("parachors");
                    if (pn.kind == Dictionary::Kind::word
                        && comp.getWord("parachors") == "fromTable")
                        parachors = pvt::parachorsFromTableSi(mix.components());
                    else if (pn.kind == Dictionary::Kind::list)
                        parachors = comp.getScalarList("parachors");
                    else
                        cfFatal("steadyNetwork: compositional parachors must be"
                                " a list or 'fromTable'");
                    if (parachors.size() != mix.size())
                        cfFatal("steadyNetwork: compositional parachors must give"
                                " one value per component");
                }

                const std::size_t nc = net.pvtFluid.pressure.size();
                if (nc < 2)
                    cfFatal("steadyNetwork: compositional needs a pressure grid");
                if (liveFlash)
                {
                    // No rho/mu table is filled: evaluateAt() flashes the feed
                    // at each station pressure instead. The grid is still
                    // validated and kept as the case's declared envelope.
                    for (std::size_t i = 0; i < nc; ++i)
                    {
                        const double p = net.pvtFluid.pressure[i];
                        if (!(p > 0.0))
                            cfFatal("steadyNetwork: compositional pressures must be positive");
                        if (i > 0 && !(p > net.pvtFluid.pressure[i - 1]))
                            cfFatal("steadyNetwork: compositional pressures must ascend");
                    }
                    std::shared_ptr<LiveFlashState> st =
                        std::make_shared<LiveFlashState>();
                    st->mix = mix;
                    st->T = T;
                    st->parachors = parachors;
                    st->sigmaFromP = sigmaFromP;
                    st->criticalVolumes = criticalVolumes;
                    net.pvtFluid.liveFlash = st;
                    if (sigmaFromP)
                    {
                        // Publish the same WK sigma as a per-pressure callback
                        // for the Beggs-Brill two-phase march (see Network.h).
                        std::shared_ptr<LiveFlashState> captured = st;
                        net.sigmaOfP = [captured](double p)
                        {
                            return pvt::weinaugKatzSurfaceTension(
                                captured->mix, captured->T, p,
                                captured->parachors);
                        };
                    }
                }
                else
                {
                    net.pvtFluid.rho.resize(nc);
                    net.pvtFluid.mu.resize(nc);
                    for (std::size_t i = 0; i < nc; ++i)
                    {
                        const double p = net.pvtFluid.pressure[i];
                        if (!(p > 0.0))
                            cfFatal("steadyNetwork: compositional pressures must be positive");
                        if (i > 0 && !(p > net.pvtFluid.pressure[i - 1]))
                            cfFatal("steadyNetwork: compositional pressures must ascend");
                        const pvt::EosState st = pvt::pengRobinsonState(mix, T, p, pref);
                        net.pvtFluid.rho[i] = st.massDensity;
                        const std::vector<double>* vcs = criticalVolumes.empty()
                            ? nullptr : &criticalVolumes;
                        net.pvtFluid.mu[i] = pvt::lbcViscosity(mix, T, p, pref, vcs);
                        if (!(net.pvtFluid.rho[i] > 0.0) || !(net.pvtFluid.mu[i] > 0.0))
                            cfFatal("steadyNetwork: compositional generated non-positive"
                                    " rho/mu at pressure " + std::to_string(p));
                    }
                }
            }
            else
            {
                net.pvtFluid.rho = pd.getScalarList("rho");
                net.pvtFluid.mu = pd.getScalarList("mu");
            }
            const std::size_t n = net.pvtFluid.pressure.size();
            // D10/D7: a live-flash fluid has no rho/mu table to validate.
            if (!net.pvtFluid.liveFlash)
            {
            if (n < 2 || net.pvtFluid.rho.size() != n
                || net.pvtFluid.mu.size() != n)
                cfFatal("steadyNetwork: fluid pvt needs pressure/rho/mu"
                        " lists of equal length (at least 2)");
            for (std::size_t i = 0; i < n; ++i)
            {
                if (!(net.pvtFluid.pressure[i] > 0)
                    || !(net.pvtFluid.rho[i] > 0)
                    || !(net.pvtFluid.mu[i] > 0))
                    cfFatal("steadyNetwork: fluid pvt entries must be"
                            " positive");
                if (i > 0 && !(net.pvtFluid.pressure[i]
                               > net.pvtFluid.pressure[i - 1]))
                    cfFatal("steadyNetwork: fluid pvt pressures must"
                            " strictly ascend");
            }
            }
        }
        else
        {
            net.fluid.rho = fd.getScalar("rho");
            net.fluid.mu = fd.getScalar("mu");
        }
    }
    if (!net.hasBlackOil && !net.hasPvtFluid
        && (!(net.fluid.rho > 0) || !(net.fluid.mu > 0)))
        cfFatal("steadyNetwork: fluid needs positive rho and mu");

    if (caseDict.has("nodes"))
    {
        const Dictionary& nodes = caseDict.getDict("nodes");
        for (const std::string& name : nodes.keys())
        {
            const label idx = net.addNode(name);
            const Dictionary& nd = nodes.getDict(name);
            if (nd.has("pressure")) net.nodePressureFixed[idx] = nd.getScalar("pressure");
        }
    }

    if (!caseDict.has("pipes")) cfFatal("steadyNetwork: case has no pipes{} block");
    const Dictionary& pipes = caseDict.getDict("pipes");
    for (const std::string& name : pipes.keys())
    {
        const Dictionary& pd = pipes.getDict(name);
        NetworkPipe p;
        p.name = name;
        p.fromName = pd.getWord("from");
        p.toName = pd.getWord("to");
        p.length = pd.getScalar("length");
        p.diameter = pd.getScalar("diameter");
        p.roughness = pd.getScalarOr("roughness", 0.0);
        p.elevation = pd.getScalarOr("elevation", 0.0);
        p.valveK = pd.getScalarOr("valve", 0.0);
        p.pumpHead = pd.getScalarOr("pumpHead", 0.0);
        p.pumpRate = pd.getScalarOr("pumpRate", 0.0);
        p.pumpShutoffHead = pd.getScalarOr("pumpShutoffHead", p.pumpHead);
        p.pumpStages = pd.getScalarOr("pumpStages", 1.0);
        p.pumpSpeed = pd.getScalarOr("pumpSpeed", 1.0);
        p.pumpGasFraction = pd.getScalarOr("pumpGasFraction", 0.0);
        p.pumpIntakePressure = pd.getScalarOr("pumpIntakePressure", 0.0);
        p.coupledResistance = pd.getScalarOr("coupledResistance", 0.0);
        // 5.191 (G7): the positive-displacement machine. displacement is the
        // switch (m3/rev); with it the row's machine curve is the machine's
        // own slip law, so a centrifugal pumpHead on the same pipe is refused
        // (one row, one machine), as is a compressor and a coupled edge.
        if (pd.has("pcp"))
        {
            const Dictionary& pc = pd.getDict("pcp");
            p.pcpOn = true;
            p.pcpDisplacement = pc.getScalarOr("displacement", 0.0);
            p.pcpSpeed = pc.getScalarOr("speed", 0.0);
            p.pcpSlip = pc.getScalarOr("slip", 0.0);
            if (!(p.pcpDisplacement > 0))
                cfFatal("steadyNetwork: pipe '" + name + "' pcp needs a"
                        " positive displacement (m3/rev): it is the volume one"
                        " revolution moves before slip");
            if (!(p.pcpSpeed > 0))
                cfFatal("steadyNetwork: pipe '" + name + "' pcp needs a"
                        " positive speed (rev/s)");
            if (!(p.pcpSlip > 0))
                cfFatal("steadyNetwork: pipe '" + name + "' pcp needs a"
                        " positive slip (m3/(s Pa)): q = nV - C dp, and a zero"
                        " C is a machine that holds any pressure at any rate");
            if (p.pumpHead > 0)
                cfFatal("steadyNetwork: pipe '" + name + "' declares both a"
                        " pcp and a centrifugal pumpHead: one row carries one"
                        " machine curve");
        }
        // 5.192 (G8): the jet pump. `nozzleArea` is the switch. The row
        // carries the COMBINED stream (q_d = q_n + q_s) and the machine
        // splits it, so the entrained share is q - q_n and the other stream
        // has nowhere else to come from on a lumped nodal balance.
        if (pd.has("jetPump"))
        {
            const Dictionary& jd = pd.getDict("jetPump");
            p.jetOn = true;
            p.jetMotiveRate = jd.getScalarOr("motiveRate", 0.0);
            p.jetNozzleArea = jd.getScalarOr("nozzleArea", 0.0);
            p.jetThroatArea = jd.getScalarOr("throatArea", 0.0);
            p.jetThroatLength = jd.getScalarOr("throatLength", 0.0);
            p.jetFrictionFactor = jd.getScalarOr("frictionFactor", 0.0);
            if (!(p.jetNozzleArea > 0.0))
                cfFatal("steadyNetwork: pipe '" + name + "' jetPump needs a"
                        " positive nozzleArea: it is the bore the motive"
                        " stream is accelerated through, and the switch");
            if (!(p.jetThroatArea > p.jetNozzleArea))
                cfFatal("steadyNetwork: pipe '" + name + "' jetPump needs"
                        " throatArea > nozzleArea: the annulus A_t - A_n is"
                        " what the entrained stream comes in through, and the"
                        " throat balance divides by m - 1");
            if (!(p.jetMotiveRate > 0.0))
                cfFatal("steadyNetwork: pipe '" + name + "' jetPump needs a"
                        " positive motiveRate: with no jet there is no"
                        " momentum to hand to the entrained stream");
            if (p.jetThroatLength < 0.0 || p.jetFrictionFactor < 0.0)
                cfFatal("steadyNetwork: pipe '" + name + "' jetPump throat"
                        " length and friction factor are non-negative");
            if (p.pumpHead > 0 || p.pcpOn || p.rodPlungerDiameter > 0)
                cfFatal("steadyNetwork: pipe '" + name + "' declares a"
                        " jetPump beside another machine: one row carries one"
                        " machine curve");
        }
        // Gas-lift valve (5.89, D1 increment 3). The port bore is `diameter`;
        // liftDomePressure is the switch, and addPipe owns the parameter scope.
        p.liftDomePressure = pd.getScalarOr("liftDomePressure", 0.0);
        p.liftAreaRatio = pd.getScalarOr("liftAreaRatio", 0.0);
        p.liftSpringPressure = pd.getScalarOr("liftSpringPressure", 0.0);
        p.liftDischargeCoeff = pd.getScalarOr("liftDischargeCoeff", 0.0);
        p.liftKappa = pd.getScalarOr("liftKappa", 0.0);
        // Rod pump (5.90, D1 increment 4). The barrel bore is the switch, and
        // addPipe owns the declaration scope. Absent keys leave every field at
        // its 0 and the pipe is read as the pipe it was before this increment.
        p.rodPlungerDiameter = pd.getScalarOr("rodPlungerDiameter", 0.0);
        p.rodStrokeLength = pd.getScalarOr("rodStrokeLength", 0.0);
        p.rodStrokesPerMinute = pd.getScalarOr("rodStrokesPerMinute", 0.0);
        p.rodDiameter = pd.getScalarOr("rodDiameter", 0.0);
        p.rodLength = pd.getScalarOr("rodLength", 0.0);
        p.rodYoungModulus = pd.getScalarOr("rodYoungModulus", 0.0);
        p.rodDensity = pd.getScalarOr("rodDensity", 0.0);
        if (p.rodPlungerDiameter < 0)
            cfFatal("steadyNetwork: pipe '" + name
                    + "' has a negative rodPlungerDiameter: a barrel bore is a "
                      "length (0 means no pump)");
        if (p.rodPlungerDiameter == 0
            && (pd.has("rodStrokeLength") || pd.has("rodStrokesPerMinute")
                || pd.has("rodDiameter") || pd.has("rodLength")
                || pd.has("rodYoungModulus") || pd.has("rodDensity")))
            cfFatal("steadyNetwork: pipe '" + name
                    + "' declares a rod-pump parameter without "
                      "rodPlungerDiameter: the barrel bore is what switches the "
                      "pump on, and the rest would be read by nothing");
        if (p.rodPlungerDiameter > 0 && net.hasBlackOil)
            cfFatal("steadyNetwork: pipe '" + name
                    + "' is a rod pump on a black-oil network: the barrel fills "
                      "with the fluid at the intake pressure and this increment "
                      "does not price a two-phase column above it");
        if (p.rodPlungerDiameter > 0 && net.hasPvtFluid)
            cfFatal("steadyNetwork: pipe '" + name
                    + "' is a rod pump on a property-varying fill: the node "
                      "model books each row's volume at that row's own upstream "
                      "density, so a barrel that compresses what it takes would "
                      "be counted twice on the way to the surface; this "
                      "increment prices an incompressible barrel on a constant "
                      "fluid { rho mu }");
        if (p.liftDomePressure < 0)
            cfFatal("steadyNetwork: pipe '" + name
                    + "' has a negative liftDomePressure: a dome charge is an "
                      "absolute pressure (0 means no valve)");
        if (p.liftDomePressure == 0
            && (pd.has("liftAreaRatio") || pd.has("liftSpringPressure")
                || pd.has("liftDischargeCoeff") || pd.has("liftKappa")))
            cfFatal("steadyNetwork: pipe '" + name
                    + "' declares a lift-valve parameter without "
                      "liftDomePressure: the dome charge is what switches the "
                      "valve on, and the rest would be read by nothing");
        if (p.liftDomePressure > 0 && net.hasBlackOil)
            cfFatal("steadyNetwork: pipe '" + name
                    + "' is a gas-lift valve on a black-oil network: a lift "
                      "valve injects gas into a node whose liquid already "
                      "carries its own solution gas, and this increment prices "
                      "the injection string, not the two-phase column");
        if (p.liftDomePressure > 0 && !net.hasPvtFluid)
            cfFatal("steadyNetwork: pipe '" + name
                    + "' is a gas-lift valve on a constant-property network: "
                      "the isentropic nozzle prices the lift gas at the "
                      "density of its own upstream pressure, which needs a "
                      "fluid { pvt } table");
        p.annulusInner = pd.getScalarOr("dInner", 0.0);
        if (p.annulusInner < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative dInner");
        // Eccentric annulus (5.87, C3 increment 2): the axis offset of the
        // inner bore, in metres. Absent or 0 is concentric and the closed
        // Fredrickson-Bird branch below is untouched. addPipe owns the scope
        // checks (contact, no-dInner, power-law) and resolves the cross-section
        // integral once there, so the rate solve reads a number.
        p.eccentricity = pd.getScalarOr("eccentricity", 0.0);
        // Writable flow-pattern dispatch (5.88, D4 increment 2): which
        // Beggs-Brill holdup law prices this two-phase pipe. -1 (absent) asks
        // the map, which is 5.85's path untouched; addPipe owns the scope
        // checks (single-phase network, the transition blend, out of range).
        p.flowPattern = static_cast<int>(pd.getScalarOr("flowPattern", -1.0));
        if (pd.has("flowPattern")
            && pd.getScalar("flowPattern")
                   != static_cast<double>(p.flowPattern))
        {
            cfFatal("steadyNetwork: pipe '" + name
                    + "' declares a fractional flowPattern: it is a regime "
                      "index (0 segregated, 2 intermittent, 3 distributed), "
                      "not a fraction");
        }
        // 5.143 (G4): WHICH published holdup law prices this two-phase
        // pipe. A WORD, not an index: a number here would silently alias
        // the flowPattern regime indices. addPipe owns the scope checks
        // (single-phase network, the flowPattern clash).
        p.correlation = 0;
        if (pd.has("correlation"))
        {
            const std::string cw = pd.getWord("correlation");
            if (cw == "beggsBrill") p.correlation = 0;
            else if (cw == "hagedornBrown") p.correlation = 1;
            else if (cw == "dunsRos") p.correlation = 2;
            else if (cw == "olgasSteady") p.correlation = 3;
            else
                cfFatal("steadyNetwork: pipe '" + name + "' unknown "
                        "correlation '" + cw + "': the published holdup laws "
                        "are beggsBrill, hagedornBrown, dunsRos and "
                        "olgasSteady");
        }
        // 2026-11-25: opt-in Griffith bubble/mist branches for the classic
        // Hagedorn-Brown vertical method. Only legal with the hagedornBrown
        // holdup law - a different correlation has no Griffith chart to
        // supplement, so a declaration there is refused loudly.
        p.hbGriffith = false;
        if (pd.has("hbGriffith"))
        {
            const std::string gw = pd.getWord("hbGriffith");
            if (gw == "true") p.hbGriffith = true;
            else if (gw == "false") p.hbGriffith = false;
            else cfFatal("steadyNetwork: pipe '" + name
                         + "' hbGriffith must be true or false");
            if (p.correlation != 1)
                cfFatal("steadyNetwork: pipe '" + name
                        + "' hbGriffith needs correlation hagedornBrown");
        }
        // 5.144 (G6): the passive steady temperature march. Every key is
        // optional; all absent is the isothermal march, bit for bit. The
        // scope rules (marched fluid, c_p, wall model, pvt inlet T) live
        // in addPipe; here the grammar and the SHAPE of the declaration
        // are checked, because only the parser can see which keys were
        // actually written.
        const bool thermalDeclared =
            pd.has("jouleThomson") || pd.has("wallHeatTransfer")
            || pd.has("insulation")
            || pd.has("wallHeatTransferSegments") || pd.has("wallModel")
            || pd.has("wallTemperature") || pd.has("surfaceTemperature")
            || pd.has("geothermalGradient") || pd.has("heatCapacity")
            || pd.has("temperatureFeedback");
        // 5.189: `jouleThomson eos` / `jouleThomson eosSrk` asks the fluid's
        // own equation of state for the coefficient at each segment's state.
        // It needs a composition, so it is refused on the black-oil and
        // plain-table paths instead of silently pricing zero.
        p.jouleThomson = 0.0;
        p.jouleThomsonFromEos = false;
        p.jouleThomsonEosIsSrk = false;
        if (pd.has("jouleThomson"))
        {
            if (pd.find("jouleThomson").kind == Dictionary::Kind::word)
            {
                const std::string jt = pd.getWord("jouleThomson");
                if (jt == "eos" || jt == "eosSrk")
                {
                    p.jouleThomsonFromEos = true;
                    p.jouleThomsonEosIsSrk = (jt == "eosSrk");
                    if (!net.pvtFluid.composition)
                        cfFatal("steadyNetwork: pipe '" + name
                                + "' declares jouleThomson " + jt
                                + " but the fluid has no composition: an "
                                  "EOS-derived Joule-Thomson needs a "
                                  "compositional pvt fluid");
                }
                else
                    cfFatal("steadyNetwork: pipe '" + name
                            + "' unknown jouleThomson '" + jt
                            + "': give a number in K/Pa, or eos / eosSrk");
            }
            else
                p.jouleThomson = pd.getScalar("jouleThomson");
        }
        p.wallHeatTransfer = pd.getScalarOr("wallHeatTransfer", 0.0);
        p.wallHeatTransferSegments =
            pd.has("wallHeatTransferSegments")
                ? pd.getScalarList("wallHeatTransferSegments")
                : std::vector<double>{};
        p.wallModel = 0;
        if (pd.has("wallModel"))
        {
            const std::string wm = pd.getWord("wallModel");
            if (wm == "constant") p.wallModel = 0;
            else if (wm == "geothermal") p.wallModel = 1;
            else
                cfFatal("steadyNetwork: pipe '" + name + "' unknown wallModel '"
                        + wm + "': the steady thermal wall is constant or "
                          "geothermal");
        }
        p.wallTemperature = pd.getScalarOr("wallTemperature", 0.0);
        p.surfaceTemperature = pd.getScalarOr("surfaceTemperature", 0.0);
        p.geothermalGradient = pd.getScalarOr("geothermalGradient", 0.0);
        p.heatCapacity = pd.getScalarOr("heatCapacity", 0.0);
        // 5.194 (G6d): the feedback switch. `temperatureFeedback true` turns
        // the passive temperature march into a coupled one: every segment
        // prices its rho/mu/sigma at the T the march has reached, and the mass
        // rate is booked at (pUp, T_in). It is refused wherever the fluid
        // cannot answer a (T, p) question, because there it would be a
        // declaration read by nothing:
        //   - a tabulated pvt table is ONE temperature's snapshot (rho(p), mu(p)),
        //     and inverting a temperature out of it would be inventing a model;
        //   - the black-oil layer is isothermal by construction (its density /
        //     viscosity correlations never see T).
        // `false` is accepted and means the quiet default, exactly like
        // hbGriffith above: the key exists to be stated, not to be guessed.
        p.pvtTemperatureFeedback = false;
        if (pd.has("temperatureFeedback"))
        {
            const std::string fb = pd.getWord("temperatureFeedback");
            if (fb == "true") p.pvtTemperatureFeedback = true;
            else if (fb == "false") p.pvtTemperatureFeedback = false;
            else
                cfFatal("steadyNetwork: pipe '" + name + "' temperatureFeedback"
                        " must be true or false");
            if (p.pvtTemperatureFeedback)
            {
                if (net.hasBlackOil)
                    cfFatal("steadyNetwork: pipe '" + name + "' declares "
                            "temperatureFeedback on a black-oil network: the "
                            "black-oil layer is isothermal by construction "
                            "(rho and mu are functions of p at the stream's "
                            "one declared temperature), so the switch would "
                            "be read by nothing");
                if (!net.pvtFluid.liveFlash)
                    cfFatal("steadyNetwork: pipe '" + name + "' declares "
                            "temperatureFeedback on a fluid without a live "
                            "flash: a tabulated rho(p)/mu(p) table is one "
                            "temperature's snapshot, so re-pricing it at the "
                            "marched T would be inventing a temperature "
                            "model; declare compositional{ liveFlash true; }");
            }
        }
        // 5.188 (G6b): a LAYERED wall prices the per-unit-length coefficient
        // from the stack itself instead of from a declared number:
        //   U = 1 / ( sum_i ln(r_{i+1}/r_i)/(2 pi k_i) + 1/(2 pi r_out h_o) )
        // which is the series conductance of the cylindrical shells, plus the
        // optional outer film. Everything the march reads is therefore
        // re-derivable from the case text alone (the checker does exactly
        // that). A missing film is "the outside is at the wall temperature"
        // and contributes nothing; an explicit zero is refused, because a
        // declared-but-inert key is the thing this codebase does not carry.
        const bool insulatedDeclared = pd.has("insulation");
        if (insulatedDeclared)
        {
            const Dictionary& ins = pd.getDict("insulation");
            p.insulatedLayered = true;
            p.insulationInnerRadius = ins.getScalarOr("innerRadius", 0.0);
            p.layerThickness = ins.getScalarList("thickness");
            p.layerConductivity = ins.getScalarList("conductivity");
            p.insulationOuterFilm = ins.getScalarOr("outerCoefficient", 0.0);
            if (!(p.insulationInnerRadius > 0.0))
                cfFatal("steadyNetwork: pipe '" + name + "' insulation needs a "
                        "positive innerRadius: the stack is radial from there");
            if (p.layerThickness.size() != p.layerConductivity.size()
                || p.layerThickness.empty())
                cfFatal("steadyNetwork: pipe '" + name + "' insulation needs "
                        "one conductivity per thickness entry (got "
                        + std::to_string(p.layerThickness.size()) + " / "
                        + std::to_string(p.layerConductivity.size()) + ")");
            for (double d : p.layerThickness)
                if (!(d > 0.0))
                    cfFatal("steadyNetwork: pipe '" + name + "' insulation "
                            "thickness must be positive: a zero shell is no "
                            "layer");
            for (double k : p.layerConductivity)
                if (!(k > 0.0))
                    cfFatal("steadyNetwork: pipe '" + name + "' insulation "
                            "conductivity must be positive");
            if (p.insulationOuterFilm < 0.0)
                cfFatal("steadyNetwork: pipe '" + name + "' negative "
                        "insulation outerCoefficient: the film is a positive "
                        "number or absent");
            double r = p.insulationInnerRadius;
            double series = 0.0;
            for (std::size_t i = 0; i < p.layerThickness.size(); ++i)
            {
                const double rNext = r + p.layerThickness[i];
                series += std::log(rNext / r) / (2.0 * M_PI
                                                 * p.layerConductivity[i]);
                r = rNext;
            }
            p.insulationOuterRadius = r;
            if (p.insulationOuterFilm > 0.0)
                series += 1.0 / (2.0 * M_PI * r * p.insulationOuterFilm);
            p.wallHeatTransfer = 1.0 / series;
            if (pd.has("wallHeatTransfer")
                || pd.has("wallHeatTransferSegments"))
                cfFatal("steadyNetwork: pipe '" + name + "' declares both "
                        "insulation and a wallHeatTransfer: the layered stack "
                        "prices the coefficient, so a declared one is read by "
                        "nothing");
        }
        if (thermalDeclared || insulatedDeclared)
        {
            if (!net.hasBlackOil && !net.hasPvtFluid)
                cfFatal("steadyNetwork: pipe '" + name + "' declares thermal "
                        "keys on a network without a marched fluid: the "
                        "steady temperature march (5.144) rides a blackoil "
                        "or pvt march and a constant-property pipe has "
                        "nothing to carry it");
            if (p.wallHeatTransfer < 0.0)
                cfFatal("steadyNetwork: pipe '" + name + "' negative "
                        "wallHeatTransfer: the overall coefficient per unit "
                        "length folds resistances into one positive number");
            for (double u : p.wallHeatTransferSegments)
                if (u < 0.0)
                    cfFatal("steadyNetwork: pipe '" + name + "' negative "
                            "wallHeatTransferSegments entry: a per-segment "
                            "coefficient is a positive number or zero "
                            "(insulated)");
            // 5.189: mu = R T^2 /(p c_p) dZ/dT - the EOS route cannot be
            // priced without the stream's c_p, and a silent zero would be a
            // declared-but-inert key.
            if (p.jouleThomsonFromEos && !(p.heatCapacity > 0.0))
                cfFatal("steadyNetwork: pipe '" + name + "' declares "
                        "jouleThomson eos without a positive heatCapacity: "
                        "an EOS-derived Joule-Thomson is R T^2/(p c_p) dZ/dT");
            if (!p.thermalOn())
                cfFatal("steadyNetwork: pipe '" + name + "' declares thermal "
                        "keys but no heat: a thermal march needs a nonzero "
                        "jouleThomson or a positive wallHeatTransfer / "
                        "wallHeatTransferSegments entry");
            if (!p.wallHeatTransferSegments.empty())
            {
                if (p.wallHeatTransfer > 0.0)
                    cfFatal("steadyNetwork: pipe '" + name + "' declares both "
                            "a single wallHeatTransfer and a "
                            "wallHeatTransferSegments list: price the wall "
                            "one way or per march segment, not both");
                if (static_cast<int>(p.wallHeatTransferSegments.size()) != 24)
                    cfFatal("steadyNetwork: pipe '" + name + "' "
                            "wallHeatTransferSegments needs exactly 24 "
                            "entries, one per march segment (got "
                            + std::to_string(
                                  p.wallHeatTransferSegments.size()) + ")");
            }
            if (p.wallHeatTransferSegments.empty() && p.wallHeatTransfer <= 0.0)
            {
                // pure Joule-Thomson readout: no wall to relax toward, so
                // every wall key would be read by nothing
                if (pd.has("wallModel") || pd.has("wallTemperature")
                    || pd.has("surfaceTemperature")
                    || pd.has("geothermalGradient")
                    || (pd.has("heatCapacity") && !p.jouleThomsonFromEos))
                    cfFatal("steadyNetwork: pipe '" + name + "' declares wall "
                            "keys without any wall exchange: a pure "
                            "Joule-Thomson march reads none of them");
            }
            else
            {
                if (p.wallModel == 0)
                {
                    if (!pd.has("wallTemperature"))
                        cfFatal("steadyNetwork: pipe '" + name + "' wall "
                                "exchange with the constant wall model needs "
                                "wallTemperature");
                    if (pd.has("surfaceTemperature")
                        || pd.has("geothermalGradient"))
                        cfFatal("steadyNetwork: pipe '" + name + "' declares "
                                "geothermal keys with wallModel constant: the "
                                "wall is at one temperature, the profile is "
                                "read by nothing");
                }
                else
                {
                    if (!pd.has("surfaceTemperature"))
                        cfFatal("steadyNetwork: pipe '" + name + "' geothermal "
                                "wall model needs surfaceTemperature: the "
                                "profile is anchored at the inlet datum");
                    if (pd.has("wallTemperature"))
                        cfFatal("steadyNetwork: pipe '" + name + "' declares "
                                "wallTemperature with wallModel geothermal: "
                                "the wall temperature comes from the profile");
                }
            }
            if (p.wallTemperature < 0.0 || p.surfaceTemperature < 0.0
                || p.geothermalGradient < 0.0 || p.heatCapacity < 0.0)
                cfFatal("steadyNetwork: pipe '" + name + "' negative thermal "
                        "state: temperatures are kelvin, the gradient is "
                        "kelvin per metre and c_p is joules per kilogram "
                        "kelvin");
        }
        // Compressor (5.156, P10): a steady pressure booster on the line.
        // `pressureRatio` is the switch; the curve reads the rate through
        // compressorRatioOf and the machine's thermodynamics (kappa, R_spec,
        // Z, efficiency) price its power and discharge-temperature readouts.
        // The pvt-table fill only, like the 5.89 valve: the boost is a
        // pressure-level statement and the readouts need the stream's own
        // rho(p) march and a declared inlet temperature.
        p.compressorRatio = pd.getScalarOr("pressureRatio", 0.0);
        p.compressorShutoffRatio = pd.getScalarOr("shutoffRatio", 0.0);
        p.compressorRatedRate = pd.getScalarOr("ratedRate", 0.0);
        p.compressorEfficiency = pd.getScalarOr("polytropicEfficiency", 1.0);
        p.compressorKappa = pd.getScalarOr("kappa", 0.0);
        p.compressorGasConstant = pd.getScalarOr("specificGasConstant", 0.0);
        p.compressorZ = pd.getScalarOr("compressibility", 1.0);
        if (p.compressorRatio < 0)
            cfFatal("steadyNetwork: pipe '" + name
                    + "' has a negative pressureRatio: a pressure ratio is at "
                      "least 1 (absent means no machine)");
        if (p.compressorRatio == 0
            && (pd.has("shutoffRatio") || pd.has("ratedRate")
                || pd.has("polytropicEfficiency") || pd.has("kappa")
                || pd.has("specificGasConstant") || pd.has("compressibility")))
            cfFatal("steadyNetwork: pipe '" + name
                    + "' declares a compressor parameter without "
                      "pressureRatio: the ratio is what switches the machine "
                      "on, and the rest would be read by nothing");
        if (p.compressorRatio > 0)
        {
            if (!(p.compressorRatio > 1.0))
                cfFatal("steadyNetwork: pipe '" + name
                        + "' has pressureRatio at or below unity: a machine "
                          "there is a restriction, not a compressor");
            if (p.compressorShutoffRatio > 0
                && p.compressorShutoffRatio < p.compressorRatio)
                cfFatal("steadyNetwork: pipe '" + name
                        + "' has a compressor curve that rises with rate "
                          "(shutoffRatio < pressureRatio)");
            if (!(p.compressorRatedRate > 0.0))
                cfFatal("steadyNetwork: pipe '" + name
                        + "' compressor curve needs a positive ratedRate: the "
                          "rated point is where pressureRatio is quoted");
            if (!(p.compressorEfficiency > 0.0 && p.compressorEfficiency <= 1.0))
                cfFatal("steadyNetwork: pipe '" + name
                        + "' compressor polytropicEfficiency must lie in "
                          "(0,1]: it divides gas power into shaft power and "
                          "sets the discharge-temperature exponent");
            if (!(p.compressorKappa > 1.0))
                cfFatal("steadyNetwork: pipe '" + name
                        + "' compressor kappa (cp/cv) must exceed 1: the "
                          "exponent (k-1)/k is what turns a ratio into a head");
            if (!(p.compressorGasConstant > 0.0))
                cfFatal("steadyNetwork: pipe '" + name
                        + "' compressor needs a positive specificGasConstant: "
                          "the head is priced on R_spec T1");
            if (!(p.compressorZ > 0.0))
                cfFatal("steadyNetwork: pipe '" + name
                        + "' compressor compressibility must be positive: it "
                          "is a declared constant of the suction state");
            if (!net.hasPvtFluid)
                cfFatal("steadyNetwork: pipe '" + name
                        + "' is a compressor on a network without a "
                          "tabulated fill: the boost is a pressure-level "
                          "statement priced on the rho(p) march, and the "
                          "readouts need the stream's declared temperature");
            if (p.valveK > 0 || p.pumpHead > 0 || p.liftDomePressure > 0
                || p.rodPlungerDiameter > 0 || p.coupledResistance > 0)
                cfFatal("steadyNetwork: pipe '" + name
                        + "' carries a compressor and another inline element: "
                          "this increment prices one machine per line");
            if (thermalDeclared)
                cfFatal("steadyNetwork: pipe '" + name
                        + "' declares thermal keys with a compressor: the "
                          "machine owns the temperature readout on this line, "
                          "the wall march would claim it too");
            if (!(net.pvtFluid.temperature > 0.0))
                cfFatal("steadyNetwork: pipe '" + name
                        + "' compressor needs the pvt table's declared "
                          "temperature: the head and the discharge "
                          "temperature are priced on T1");
            if (!(p.length > 0.0) || !(p.diameter > 0.0))
                cfFatal("steadyNetwork: pipe '" + name
                        + "' is a compressor on a line with no bore to march: "
                          "the boost feeds the pipe law, not a bare node pair");
        }
        p.ductHeight = pd.getScalarOr("ductHeight", 0.0);
        if (p.ductHeight < 0)
            cfFatal("steadyNetwork: pipe '" + name
                    + "' has a negative ductHeight");
        p.flowIndex = pd.getScalarOr("flowIndex", 0.0);
        p.consistency = pd.getScalarOr("consistency", 0.0);
        p.yieldStress = pd.getScalarOr("yieldStress", 0.0);
        if (p.yieldStress < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative yield"
                    " stress (0 means no threshold)");
        // 5.76 (A5 increment 2): a yield stress is now PRICED in the
        // black-oil / PVT marches as the liquid's regularised apparent
        // viscosity (Papanastasiou, no dead branch in a correlation chain).
        // It is no longer refused; the transient row's static latch (5.49)
        // is a different mechanism and stays there.
        // A coupled edge MAY declare a yield stress: it describes the fluid in
        // its partner region, not a threshold this network's row holds. The
        // declaration is forwarded to the 3D block below, after every check,
        // and the wall-law fields are cleared so the 1D loss stays exactly the
        // partner's statement (no second charge, §5.43).
        if (p.flowIndex < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative"
                    " power-law flow index (0 means Newtonian)");
        if (p.flowIndex > 0 && !(p.consistency > 0))
            cfFatal("steadyNetwork: pipe '" + name + "' declares a flow index"
                    " without a positive consistency: at n != 1 the fluid's mu"
                    " is not the same quantity as K (Pa.s^n)");
        if (p.flowIndex == 0 && p.consistency > 0)
            cfFatal("steadyNetwork: pipe '" + name + "' declares a consistency"
                    " without a flowIndex: without n the consistency is not"
                    " a quantity the wall law reads");
        // 5.76 (A5 increment 2): a rheology now reprices the LIQUID share of
        // the two-phase / PVT march's viscosity on the mixture wall rate;
        // the published holdup and flow-pattern correlations do not read a
        // viscosity and stay exactly as they were.
        // Same for `flowIndex`/`consistency`: allowed on a coupled edge, but
        // only as the partner region's constitutive law.
        if (p.valveK < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative valve"
                    " coefficient");
        if (p.pumpHead < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative pump"
                    " head");
        if (p.coupledResistance < 0)
            cfFatal("steadyNetwork: pipe '" + name + "' has a negative coupled"
                    " resistance");
        if (p.coupledResistance > 0 && (net.hasBlackOil || net.hasPvtFluid))
            cfFatal("steadyNetwork: pipe '" + name + "' is a coupled element,"
                    " which the single-phase constant-property law is the only"
                    " one that carries");
        if (p.coupledResistance > 0 && (p.valveK > 0 || p.pumpHead > 0))
            cfFatal("steadyNetwork: pipe '" + name + "' is a coupled element:"
                    " its loss is the partner model's, so it cannot carry a"
                    " valve or a pump");
        // 5.191: a machine curve on a coupled edge would be a second owner of
        // the same row, and a compressor would be a third.
        if (p.pcpOn && p.coupledResistance > 0)
            cfFatal("steadyNetwork: pipe '" + name + "' is a coupled element"
                    " and declares a pcp: the row's machine rise belongs to the"
                    " partner region");
        if (p.pcpOn && net.hasBlackOil)
            cfFatal("steadyNetwork: pipe '" + name + "' declares a pcp on a"
                    " black-oil network: the black-oil march prices its lift"
                    " from the stock-tank curve family, so a displacement"
                    " machine's slip law is carried by a constant-property or"
                    " pvt-fill network only");
        if (p.pcpOn && p.compressorRatio > 0)
            cfFatal("steadyNetwork: pipe '" + name + "' declares both a pcp and"
                    " a compressor: one row, one machine");
        // 5.192: the same two constraints the pcp row carries -- the partner
        // region owns a coupled edge's rise, and the black-oil march prices
        // lift from the stock-tank curve family, so a momentum-handing
        // machine only has a constant-property or pvt-fill row to live in.
        if (p.jetOn && p.coupledResistance > 0)
            cfFatal("steadyNetwork: pipe '" + name + "' is a coupled element"
                    " and declares a jetPump: the row's machine rise belongs"
                    " to the partner region");
        if (p.jetOn && p.compressorRatio > 0)
            cfFatal("steadyNetwork: pipe '" + name + "' declares a jetPump"
                    " beside a compressor: one row carries one machine curve");
        if (p.jetOn && net.hasBlackOil)
            cfFatal("steadyNetwork: pipe '" + name + "' declares a jetPump on a"
                    " black-oil network: the momentum balance would have to be"
                    " priced on a mixture density the black-oil march does not"
                    " hand a machine row (constant-property or pvt fill only)");
        if (p.pumpHead > 0 && p.pumpShutoffHead < p.pumpHead)
            cfFatal("steadyNetwork: pipe '" + name + "' has a pump curve that"
                    " rises with rate (pumpShutoffHead < pumpHead)");
        if (p.pumpStages < 1)
            cfFatal("steadyNetwork: pipe '" + name + "' has pumpStages < 1: a "
                    "stage count is at least one");
        if (!(p.pumpSpeed > 0))
            cfFatal("steadyNetwork: pipe '" + name + "' has a non-positive "
                    "pumpSpeed: the affinity laws need a positive speed ratio");
        if ((pd.has("pumpStages") || pd.has("pumpSpeed")) && !(p.pumpHead > 0))
            cfFatal("steadyNetwork: pipe '" + name + "' declares pumpStages/"
                    "pumpSpeed without a pump (pumpHead): the affinity laws "
                    "belong to a pump curve, not a bare pipe");
        // 5.86 (D1 increment 2): Turpin gas interference. g is the free-gas
        // volume fraction at intake, so it lives in [0,1); it degrades a pump
        // curve (needs pumpHead) and the correlation needs a positive absolute
        // intake pressure (pumpIntakePressure) to price p_i.
        if (p.pumpGasFraction < 0 || p.pumpGasFraction >= 1)
            cfFatal("steadyNetwork: pipe '" + name + "' has pumpGasFraction "
                    "outside [0,1): it is a gas VOLUME fraction at the pump "
                    "intake, so it must be at least 0 and strictly below 1");
        if (p.pumpGasFraction > 0 && !(p.pumpHead > 0))
            cfFatal("steadyNetwork: pipe '" + name + "' declares pumpGasFraction "
                    "without a pump (pumpHead): gas interference degrades a pump "
                    "curve, not a bare pipe");
        if (p.pumpGasFraction > 0 && !(p.pumpIntakePressure > 0))
            cfFatal("steadyNetwork: pipe '" + name + "' declares pumpGasFraction "
                    "without a positive pumpIntakePressure: the Turpin gas "
                    "correlation is priced on the absolute intake pressure p_i");
        // The declaration has passed every check. Move it out of the wall-law
        // fields: a coupled pipe's loss is `coupledResistance * q` alone, so
        // leaving a `flowIndex` or `yieldStress` here would let the yielded
        // dead branch of `equipmentFlow` and `pipeFrictionReport` read a law
        // the partner never claimed (5.46/5.47). Zeroing is what makes this
        // change structurally unable to move any existing 1D test.
        if (p.coupledResistance > 0)
        {
            p.partnerRheology = {p.flowIndex, p.consistency, p.yieldStress};
            p.flowIndex = 0;
            p.consistency = 0;
            p.yieldStress = 0;
        }
        // 5.158 (P18) / 5.159 (P22): the opt-in slug-start forecast and the
        // erosion / hydrate / corrosion risk report. Both are pure readouts
        // of the march's own per-segment profile; the parser owns the
        // grammar and the SHAPE of the declaration. Unknown keys inside the
        // sub-blocks are refused by whitelist, because Dictionary does not
        // reject what nobody reads.
        if (pd.has("trajectory"))
            cfFatal("steadyNetwork: pipe '" + name + "' declares trajectory "
                    "at the pipe level: the terrain trajectory belongs inside "
                    "slugPrediction{}");
        if (pd.has("slugPrediction"))
        {
            const Dictionary& sp = pd.getDict("slugPrediction");
            for (const std::string& k : sp.keys())
                if (k != "trajectory" && k != "fillThreshold")
                    cfFatal("steadyNetwork: pipe '" + name + "' slugPrediction "
                            "declares unknown key '" + k + "': the slug "
                            "forecast reads trajectory and fillThreshold");
            if (!net.hasBlackOil)
                cfFatal("steadyNetwork: pipe '" + name + "' declares "
                        "slugPrediction on a network without a black-oil "
                        "march: the forecast pools the march's own holdup "
                        "profile and a single-phase pipe has none");
            p.slugDeclared = true;
            if (sp.has("fillThreshold"))
            {
                p.slugFillThreshold = sp.getScalar("fillThreshold");
                if (!(p.slugFillThreshold > 0.0))
                    cfFatal("steadyNetwork: pipe '" + name + "' slugPrediction "
                            "fillThreshold must be positive: it is the pooled "
                            "liquid level as a fraction of the uphill run");
            }
            if (!sp.has("trajectory"))
                cfFatal("steadyNetwork: pipe '" + name + "' slugPrediction "
                        "without a trajectory{}: the dips are where the "
                        "forecast pools liquid");
            const Dictionary& td = sp.getDict("trajectory");
            const std::vector<std::string>& segs = td.keys();
            if (segs.size() < 2)
                cfFatal("steadyNetwork: pipe '" + name + "' slugPrediction "
                        "trajectory needs at least two named segments: one "
                        "inclination is the pipe's own slope, not a terrain");
            double sumL = 0.0, sumZ = 0.0;
            for (const std::string& sk : segs)
            {
                const Dictionary& sd = td.getDict(sk);
                for (const std::string& k2 : sd.keys())
                    if (k2 != "length" && k2 != "angleDeg")
                        cfFatal("steadyNetwork: pipe '" + name + "' trajectory "
                                "segment '" + sk + "' declares unknown key '"
                                + k2 + "': a segment is length and angleDeg");
                TerrainTrajectorySegment seg;
                seg.length = sd.getScalar("length");
                seg.angleDeg = sd.getScalar("angleDeg");
                if (!(seg.length > 0.0))
                    cfFatal("steadyNetwork: pipe '" + name + "' trajectory "
                            "segment '" + sk + "' needs a positive length");
                if (!(seg.angleDeg >= -90.0 && seg.angleDeg <= 90.0))
                    cfFatal("steadyNetwork: pipe '" + name + "' trajectory "
                            "segment '" + sk + "' angleDeg outside [-90, 90]: "
                              "the inclination is an angle from horizontal");
                sumL += seg.length;
                sumZ += seg.length * std::sin(seg.angleDeg * M_PI / 180.0);
                p.slugTrajectory.push_back(seg);
            }
            if (std::fabs(sumL - p.length)
                > 1e-9 * std::max(std::fabs(p.length), 1.0))
                cfFatal("steadyNetwork: pipe '" + name + "' slugPrediction "
                        "trajectory lengths sum to " + std::to_string(sumL)
                        + " m but the pipe is " + std::to_string(p.length)
                        + " m: the trajectory lays the pipe's own arc out");
            if (std::fabs(sumZ - p.elevation)
                > 1e-9 * std::max(std::fabs(p.elevation), 1.0))
                cfFatal("steadyNetwork: pipe '" + name + "' slugPrediction "
                        "trajectory rises to " + std::to_string(sumZ)
                        + " m but the pipe's elevation is "
                        + std::to_string(p.elevation)
                        + " m: the trajectory must end where the pipe does");
        }
        // 5.193 (G9): the well-trajectory survey. `md` is the switch: with it
        // the pipe stops stating ONE slope and states an inclination per
        // measured depth, which is what a directional well actually delivers
        // and what the gravity term has to be priced on. The shape is checked
        // here rather than trusted, because every one of these defects would
        // otherwise reach the march as a plausible-looking number:
        //   - md not starting at 0 or not ending at length: the survey covers
        //     a different hole than the pipe;
        //   - md not monotone: an interval would be walked backwards;
        //   - unequal list lengths: a station would borrow its neighbour's
        //     inclination;
        //   - |inclination| > 90: an angle from horizontal cannot leave the
        //     half-plane, and sin(theta) would fold back;
        //   - a declared elevation beside a survey that disagrees: two
        //     statements of one arc, exactly the ambiguity 5.158 refuses in
        //     the terrain trajectory;
        //   - a survey beside a slugPrediction trajectory: the same, on the
        //     second grammar;
        //   - a survey on a coupled edge (its drop is the partner's own
        //     statement, so the gravity here would be read by nothing);
        //   - `md` alone with no inclination list, and inclination with no md.
        if (pd.has("survey"))
        {
            const Dictionary& sv = pd.getDict("survey");
            for (const std::string& k : sv.keys())
                if (k != "md" && k != "inclination")
                    cfFatal("steadyNetwork: pipe '" + name + "' survey declares"
                            " unknown key '" + k + "': a survey is md and"
                            " inclination");
            p.surveyOn = true;
            p.surveyMd = sv.getScalarList("md");
            p.surveyInclDeg = sv.getScalarList("inclination");
            if (p.surveyMd.size() < 2 || p.surveyInclDeg.size() < 2)
                cfFatal("steadyNetwork: pipe '" + name + "' survey needs at"
                        " least two stations: one station states no arc");
            if (p.surveyMd.size() != p.surveyInclDeg.size())
                cfFatal("steadyNetwork: pipe '" + name + "' survey has "
                        + std::to_string(p.surveyMd.size()) + " measured depths"
                        " but " + std::to_string(p.surveyInclDeg.size())
                        + " inclinations: one angle per station");
            if (std::fabs(p.surveyMd.front()) > 1e-12)
                cfFatal("steadyNetwork: pipe '" + name + "' survey starts at md"
                        " " + std::to_string(p.surveyMd.front())
                        + " m: a survey of this pipe starts at 0");
            if (std::fabs(p.surveyMd.back() - p.length)
                > 1e-9 * std::max(std::fabs(p.length), 1.0))
                cfFatal("steadyNetwork: pipe '" + name + "' survey ends at md "
                        + std::to_string(p.surveyMd.back()) + " m but the pipe"
                        " is " + std::to_string(p.length)
                        + " m: the survey lays the pipe's own arc out");
            for (std::size_t i = 0; i + 1 < p.surveyMd.size(); ++i)
                if (!(p.surveyMd[i + 1] > p.surveyMd[i]))
                    cfFatal("steadyNetwork: pipe '" + name + "' survey md is"
                            " not monotone at station " + std::to_string(i + 1)
                            + " (" + std::to_string(p.surveyMd[i]) + " -> "
                            + std::to_string(p.surveyMd[i + 1])
                            + " m): measured depth only grows");
            for (double a : p.surveyInclDeg)
                if (!(a >= -90.0 && a <= 90.0))
                    cfFatal("steadyNetwork: pipe '" + name + "' survey"
                            " inclination " + std::to_string(a)
                            + " deg is outside [-90, 90]: it is an angle from"
                            " horizontal");
            if (p.coupledResistance > 0)
                cfFatal("steadyNetwork: pipe '" + name + "' is a coupled edge"
                        " with a survey: its drop is the partner's own"
                        " statement, so the surveyed gravity would be read by"
                        " nothing");
            if (pd.has("slugPrediction"))
                cfFatal("steadyNetwork: pipe '" + name + "' declares both a"
                        " survey{} and a slugPrediction trajectory: two"
                        " statements of one arc, and the trajectory grammar"
                        " already carries the inclination");
            const double dz = surveyElevationOf(p);
            if (pd.has("elevation"))
            {
                if (std::fabs(dz - p.elevation)
                    > 1e-9 * std::max(std::fabs(p.elevation), 1.0))
                    cfFatal("steadyNetwork: pipe '" + name + "' survey rises"
                            " to " + std::to_string(dz)
                            + " m but the declared elevation is "
                            + std::to_string(p.elevation)
                            + " m: the survey is the arc, so the elevation is"
                            " what it integrates to");
            }
            else
            {
                // No declared elevation: the survey's integral IS the pipe's
                // elevation, so every downstream reader (single-phase
                // equipmentFlow, the buoyancy of the fills, the report) sees
                // the arc the survey describes.
                p.elevation = dz;
            }
        }
        if (pd.has("riskReport"))
        {
            const Dictionary& rr = pd.getDict("riskReport");
            static const char* kRiskKeys[] = {
                "erosionC", "hydrateTref", "hydratePRef", "hydrateDH",
                "hydrateLowK", "hydrateMediumK", "co2Fraction",
                "h2sFraction", "co2LowPa", "co2HighPa"};
            for (const std::string& k : rr.keys())
            {
                bool known = false;
                for (const char* w : kRiskKeys)
                    if (k == w) { known = true; break; }
                if (!known)
                    cfFatal("steadyNetwork: pipe '" + name + "' riskReport "
                            "declares unknown key '" + k + "': the risk report "
                            "reads erosionC, the hydrate line, the sour-gas "
                            "fractions and their bands");
            }
            if (!net.hasBlackOil)
                cfFatal("steadyNetwork: pipe '" + name + "' declares "
                        "riskReport on a network without a black-oil march: "
                        "the risks are priced on the march's own per-segment "
                        "profile and a single-phase pipe has none");
            p.riskDeclared = true;
            if (rr.has("erosionC"))
            {
                p.erosionC = rr.getScalar("erosionC");
                if (!(p.erosionC > 0.0))
                    cfFatal("steadyNetwork: pipe '" + name + "' riskReport "
                            "erosionC must be positive: the API RP 14E C "
                            "factor divides the velocity limit");
            }
            const bool hT = rr.has("hydrateTref");
            const bool hP = rr.has("hydratePRef");
            const bool hD = rr.has("hydrateDH");
            if (hT || hP || hD)
            {
                if (!(hT && hP && hD))
                    cfFatal("steadyNetwork: pipe '" + name + "' declares a "
                            "partial hydrate line: hydrateTref, hydratePRef "
                            "and hydrateDH are one dissociation curve and "
                            "must be declared together"
                            + std::string(hT ? "" : " (missing hydrateTref)")
                            + std::string(hP ? "" : " (missing hydratePRef)")
                            + std::string(hD ? "" : " (missing hydrateDH)"));
                p.hydrateDeclared = true;
                p.hydrateTref = rr.getScalar("hydrateTref");
                p.hydratePRef = rr.getScalar("hydratePRef");
                p.hydrateDH = rr.getScalar("hydrateDH");
                if (!(p.hydrateTref > 0.0 && p.hydratePRef > 0.0
                      && p.hydrateDH > 0.0))
                    cfFatal("steadyNetwork: pipe '" + name + "' hydrate line "
                            "needs positive Tref (K), pRef (Pa) and dH "
                            "(J/mol): the dissociation curve is anchored at "
                            "them");
            }
            if (rr.has("hydrateLowK"))
                p.hydrateLowK = rr.getScalar("hydrateLowK");
            if (rr.has("hydrateMediumK"))
                p.hydrateMediumK = rr.getScalar("hydrateMediumK");
            if (!(p.hydrateLowK > 0.0 && p.hydrateMediumK > 0.0))
                cfFatal("steadyNetwork: pipe '" + name + "' riskReport "
                        "hydrate bands must be positive kelvin subcoolings");
            if (!(p.hydrateLowK < p.hydrateMediumK))
                cfFatal("steadyNetwork: pipe '" + name + "' hydrate band "
                        "hydrateLowK >= hydrateMediumK: the bands must rise "
                        "with the subcooling");
            if (rr.has("co2Fraction"))
                p.co2Fraction = rr.getScalar("co2Fraction");
            if (rr.has("h2sFraction"))
                p.h2sFraction = rr.getScalar("h2sFraction");
            if (!(p.co2Fraction >= 0.0 && p.co2Fraction <= 1.0)
                || !(p.h2sFraction >= 0.0 && p.h2sFraction <= 1.0))
                cfFatal("steadyNetwork: pipe '" + name + "' riskReport "
                        "sour-gas fraction outside [0,1]: co2Fraction and "
                        "h2sFraction are mole fractions");
            if (rr.has("co2LowPa")) p.co2LowPa = rr.getScalar("co2LowPa");
            if (rr.has("co2HighPa")) p.co2HighPa = rr.getScalar("co2HighPa");
            if (!(p.co2LowPa > 0.0 && p.co2HighPa > 0.0))
                cfFatal("steadyNetwork: pipe '" + name + "' riskReport CO2 "
                        "bands must be positive pressures (Pa)");
            if (!(p.co2LowPa < p.co2HighPa))
                cfFatal("steadyNetwork: pipe '" + name + "' CO2 band "
                        "co2LowPa >= co2HighPa: the bands must rise with the "
                        "partial pressure");
        }
        // 5.223 (H3): the separator row. The `separator{}` block is the switch,
        // read last so every machine and wall key above is already in `p` and
        // the scope refusals see the whole declaration. Every rejection the
        // increment names is refused here, loudly, before the first march: a
        // row with no feed rate has nothing to flash (and separatorSolve would
        // refuse it mid-Newton, turning a case typo into an iteration
        // failure), an efficiency outside [0,1] is not a recovery, and a
        // negative dP would deliver the outlets above the stage.
        if (pd.has("separator"))
        {
            const Dictionary& sd = pd.getDict("separator");
            p.separatorOn = true;
            p.sepGasNodeName = sd.getWord("gasNode");
            p.sepFeedRate = sd.getScalarOr("feedRate", 0.0);
            p.sepEfficiency = sd.getScalarOr("efficiency", 1.0);
            p.sepDP = sd.getScalarOr("dP", 0.0);
            p.sepTemperature = sd.getScalarOr("temperature", 0.0);
            if (sd.has("feedRate") && !(p.sepFeedRate > 0.0))
                cfFatal("steadyNetwork: pipe '" + name + "' separator feedRate"
                        " is declared but not positive: a row that states its"
                        " feed has to mean it (omit the key to let the network"
                        " deliver whatever arrives)");
            if (!(p.sepEfficiency >= 0.0 && p.sepEfficiency <= 1.0))
                cfFatal("steadyNetwork: pipe '" + name + "' separator"
                        " efficiency is outside [0,1]: it is the fraction of"
                        " the condensables the liquid outlet actually carries");
            if (!(p.sepDP >= 0.0))
                cfFatal("steadyNetwork: pipe '" + name + "' separator dP is"
                        " negative: dP is the drop across the internals, the"
                        " outlets cannot sit above the stage");
            if (p.sepTemperature < 0.0)
                cfFatal("steadyNetwork: pipe '" + name + "' separator"
                        " temperature is kelvin");
            if (!net.hasBlackOil)
                cfFatal("steadyNetwork: pipe '" + name + "' is a separator row"
                        " on a non-black-oil network: the flash allocation is"
                        " the black-oil layer (Standing/VB rs, FVF, DAK z) and"
                        " nothing else can price it");
            if (p.pumpHead > 0 || p.pcpOn || p.jetOn
                || p.rodPlungerDiameter > 0 || p.liftDomePressure > 0
                || p.compressorRatio > 0)
                cfFatal("steadyNetwork: pipe '" + name + "' declares a"
                        " separator beside another machine: one row carries"
                        " one element");
            if (p.valveK != 0.0 || p.annulusInner > 0 || p.ductHeight > 0
                || p.coupledResistance > 0 || p.slugDeclared || p.riskDeclared)
                cfFatal("steadyNetwork: pipe '" + name + "' declares a"
                        " separator with a wall law, bore shape, coupled edge"
                        " or flow-assurance override: the row is the flash,"
                        " not a pipe");
        }
        // 5.226 (H3): the QMOM pipe block. `beta` is the kernel and the
        // switch: with no kernel there is no source to march, so its absence
        // (or a non-positive one) is refused rather than priced as a frozen
        // population. M0 < 0 is not a count; dThreshold <= 0 would switch
        // Brinkman on at every state including d32 = 0; and a phi the march
        // can neither read nor price is refused too. The march that carries
        // the population is the tabulated-PVT one -- the constant-property
        // law is algebraic and the black-oil chain prices no explicit mu(p)
        // to multiply -- so any other fill or any co-declared machine, bore
        // shape, rheology or thermal physics is refused loudly.
        if (pd.has("qmom"))
        {
            const Dictionary& qd = pd.getDict("qmom");
            p.qmomOn = true;
            p.qmomBeta = qd.getScalarOr("beta", 0.0);
            p.qmomM0 = qd.getScalarOr("m0", 0.0);
            p.qmomM1 = qd.getScalarOr("m1", 0.0);
            p.qmomM2 = qd.getScalarOr("m2", 0.0);
            p.qmomM3 = qd.getScalarOr("m3", 0.0);
            p.qmomDThreshold = qd.getScalarOr("dThreshold", 0.0);
            p.qmomPhi = qd.getScalarOr("phi", 0.0);
            p.qmomPhiFromM0 = qd.getScalarOr("phiFromM0", 0.0);
            if (!(p.qmomBeta > 0.0))
                cfFatal("steadyNetwork: pipe '" + name + "' qmom needs a"
                        " positive beta: it is the Smoluchowski kernel, and"
                        " with no kernel there is no source to march");
            if (p.qmomM0 < 0.0)
                cfFatal("steadyNetwork: pipe '" + name + "' qmom m0 is"
                        " negative: M0 is a number density");
            if (!(p.qmomM1 > 0.0) || !(p.qmomM2 > 0.0) || !(p.qmomM3 > 0.0))
                cfFatal("steadyNetwork: pipe '" + name + "' qmom needs"
                        " positive m1, m2 and m3: M2 is the d32 divisor and"
                        " M1 the conserved volume the source gains ride on");
            if (!(p.qmomDThreshold > 0.0))
                cfFatal("steadyNetwork: pipe '" + name + "' qmom dThreshold"
                        " must be positive: at or below zero every state"
                        " including d32 = 0 would switch on Brinkman");
            if (!qd.has("phi") && !qd.has("phiFromM0"))
                cfFatal("steadyNetwork: pipe '" + name + "' qmom declares no"
                        " volume fraction: phi (in [0,1)) or phiFromM0 (the"
                        " coefficient pricing phi from M0) must be declared"
                        " for the d32 feedback to price anything");
            if (qd.has("phi") && !(p.qmomPhi >= 0.0 && p.qmomPhi < 1.0))
                cfFatal("steadyNetwork: pipe '" + name + "' qmom phi is"
                        " outside [0,1): it is a volume fraction");
            if (qd.has("phiFromM0") && !(p.qmomPhiFromM0 > 0.0))
                cfFatal("steadyNetwork: pipe '" + name + "' qmom phiFromM0 is"
                        " declared but not positive: it prices phi from the"
                        " number density M0");
            if (!net.hasPvtFluid || net.hasBlackOil)
                cfFatal("steadyNetwork: pipe '" + name + "' declares qmom on"
                        " a network without the tabulated march: the moment"
                        " population rides the segment-by-segment PVT march"
                        " (5.226)");
            if (p.flowIndex > 0.0 || p.consistency > 0.0
                || p.yieldStress > 0.0)
                cfFatal("steadyNetwork: pipe '" + name + "' declares qmom"
                        " with a wall rheology: the increment multiplies the"
                        " Newtonian table viscosity, not a power-law or"
                        " yielded one");
            if (p.separatorOn || p.pumpHead > 0 || p.pcpOn || p.jetOn
                || p.rodPlungerDiameter > 0 || p.liftDomePressure > 0
                || p.compressorRatio > 0 || p.coupledResistance > 0
                || p.annulusInner > 0 || p.ductHeight > 0
                || p.thermalOn() || p.surveyOn
                || p.slugDeclared || p.riskDeclared)
                cfFatal("steadyNetwork: pipe '" + name + "' declares qmom"
                        " beside a machine, bore shape, coupled edge, thermal"
                        " physics or survey: the population marches on the"
                        " plain Newtonian PVT row this increment owns");
        }
        net.addPipe(p);
    }

    if (caseDict.has("sources"))
    {
        const Dictionary& src = caseDict.getDict("sources");
        for (const std::string& name : src.keys())
        {
            const Dictionary& sd = src.getDict(name);
            net.setNodeRate(sd.getWord("node"), sd.getScalar("rate"));
        }
    }

    return net;
}

SteadyControls loadSteadyControls(const Dictionary& caseDict)
{
    SteadyControls ctl;
    if (caseDict.has("solve"))
    {
        const Dictionary& sd = caseDict.getDict("solve");
        ctl.maxIter = static_cast<int>(sd.getScalarOr("maxIter", ctl.maxIter));
        ctl.tol = sd.getScalarOr("tol", ctl.tol);
        ctl.relax = sd.getScalarOr("relax", ctl.relax);
    }
    return ctl;
}

} // namespace cf
