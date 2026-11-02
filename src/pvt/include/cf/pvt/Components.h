#pragma once

#include <string>
#include <vector>

#include "cf/common/Labels.h"

namespace cf::pvt
{

/// Universal gas constant, J/(mol K).
inline constexpr scalar kGasConstant = 8.314462618;

/// One component of a hydrocarbon / acid-gas mixture.
/// Critical constants, acentric factor and molar mass are the standard
/// tabulated values for natural gas and condensate systems.
struct Component
{
    std::string name;
    scalar molarMass = 0; // kg/mol
    scalar Tc = 0;        // K
    scalar Pc = 0;        // Pa
    scalar acentric = 0;  // -
};

Component methane();
Component ethane();
Component propane();
Component iButane();
Component nButane();
Component nPentane();
Component nHexane();
Component carbonDioxide();
Component hydrogenSulphide();
Component nitrogen();

/// Lookup by the dictionary name (case sensitive, e.g. "nButane").
/// Returns the component and sets found; use the throwing variant in
/// user-facing paths so a typo in a case file is a fatal error.
Component componentByName(const std::string& name, bool& found);
Component requireComponent(const std::string& name);

/// Feed composition: components + mole fractions (kept normalised by the
/// caller through Mixture::normalise()).
class Mixture
{
public:
    Mixture() = default;
    explicit Mixture(std::vector<Component> comps) : comps_(std::move(comps)) {}

    const std::vector<Component>& components() const { return comps_; }
    std::size_t size() const { return comps_.size(); }

    void add(const Component& c, scalar moleFraction);
    /// Append a component with an explicit fraction and return its index.
    std::size_t index(const std::string& name) const;

    std::vector<scalar>& fractions() { return z_; }
    const std::vector<scalar>& fractions() const { return z_; }
    void setFractions(std::vector<scalar> z);
    void normalise();

    scalar molarMass() const; // kg/mol of the mixture
    scalar massFraction(std::size_t i) const;

    std::string describe() const; // "methane 0.500000 ethane ..."

private:
    std::vector<Component> comps_;
    std::vector<scalar> z_;
};

} // namespace cf::pvt
