#pragma once

#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/fvm/Fields.h"
#include "cf/fvm/FvMatrix.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// D8 increment 5: multispecies scalar transport with one irreversible global
/// step  F + s O -> products, closed with the published Arrhenius rate.
///
/// A single global reaction is the standard flamelet / infinitely-fast-chemistry
/// building block (Turns, "Thermodynamic Properties of Combustion" and the
/// standard Arrhenius kinetics in Turns "An Introduction to Combustion", ch. 2).
/// The law itself is NOT transcribed from an unverified correlation: the
/// temperature dependence is the Arrhenius exponential exp(-Ea/(R T)), and the
/// pre-exponential A, activation energy Ea and stoichiometric coefficient s
/// are case DECLARED inputs (the same stance as kij / parachor / pump curves).
///
/// Species mass fractions solve steady convection-diffusion-reaction:
///   div(u Y_f) - div(D grad Y_f) = -omega / rho      (fuel)
///   div(u Y_o) - div(D grad Y_o) = -s omega / rho    (oxidiser)
///   omega = rho A exp(-Ea/(R T)) Y_f Y_o             [kg/m3/s]
///
/// The `rho` cancels in the species equation when the rate is written per unit
/// mass, so the solver works with the kinematic rate k(T) = A exp(-Ea/(R T)):
///   div(u Y_f) - div(D grad Y_f) = -k Y_f Y_o
///   div(u Y_o) - div(D grad Y_o) = -s k Y_f Y_o
/// Temperature is a declared uniform field (one-way: chemistry does NOT feed
/// back into energy yet). When the oxidiser is held fixed at a declared value
/// the fuel equation is linear, div(u Y_f) - div(D grad Y_f) + k Y_o Y_f = 0,
/// and at rest has the closed form Y_f(x) = C1 cosh(p x) + C2 sinh(p x) with
/// p = sqrt(k Y_o / D) -- the anchor the checker reconstructs.
struct CombustionConfig
{
    bool enabled = false;
    scalar diffusivity = 0;      // D [m2/s], species molecular diffusivity
    scalar preExp = 0;           // A [1/s], Arrhenius pre-exponential
    scalar activationE = 0;      // Ea [J/mol]
    scalar gasConstant = 8.314462618; // R [J/mol/K], CODATA 2019 exact
    scalar temperature = 0;      // T [K], uniform declared field
    scalar stoich = 1;           // s, oxidiser consumed per unit fuel
    scalar fixedOxidiser = -1;   // >=0 pins Y_o to this value (linear anchor);
                                 // <0 means both species are solved coupled
    scalar fuelInlet = 0;        // Y_f value at patches declaring `fuelValue`
    scalar oxInlet = 0;          // Y_o value at patches declaring `oxValue`

    scalar rateConstant() const; // k(T) = A exp(-Ea/(R T))
    void validate() const;
};

CombustionConfig readCombustionConfig(const Dictionary& caseDict);

struct CombustionResult
{
    Field<scalar> Yf;          // fuel mass fraction per cell
    Field<scalar> Yo;          // oxidiser mass fraction per cell
    Field<scalar> reactionRate; // k Yf Yo per cell [1/s]
    scalar fuelIn = 0, fuelOut = 0, oxIn = 0, oxOut = 0;
    scalar fuelConsumed = 0;   // integral of -k Yf Yo over the domain [m3/s]
    scalar oxConsumed = 0;     // integral of -s k Yf Yo
    int fuelIters = 0, oxIters = 0;
    bool solved = false;
};

/// Solve the coupled (or oxidiser-fixed) steady species equations on `mesh`.
/// `phi` is the continuity flux (from a frozen carrier or a declared velocity),
/// exactly the flux the scalar transport app uses. `Yf` / `Yo` arrive with their
/// patch BCs already set by the app (`fuelValue` / `oxValue` = fixedValue,
/// otherwise zeroGradient); the solver writes only their internal fields.
CombustionResult solveCombustion(const polyMesh& mesh, const surfaceScalarField& phi,
                                 volScalarField& Yf, volScalarField& Yo,
                                 const CombustionConfig& cfg, const Schemes& schemes,
                                 bool fixedOxidiser);

/// Published Arrhenius rate constant k(T) = A exp(-Ea/(R T)). Exposed so a
/// checker can grade the temperature dependence independently of the solver.
scalar arrheniusRate(scalar preExp, scalar activationE, scalar gasConstant, scalar T);

} // namespace cf
