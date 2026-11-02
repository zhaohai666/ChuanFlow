#include <algorithm>
#include <numeric>
#include <sstream>

#include "cf/common/Error.h"
#include "cf/pvt/Components.h"

namespace cf::pvt
{

namespace
{

Component makeComponent(const char* name, scalar molarMassGrams,
                        scalar Tc, scalar PcBar, scalar acentric)
{
    Component c;
    c.name = name;
    c.molarMass = molarMassGrams * 1.0e-3; // g/mol -> kg/mol
    c.Tc = Tc;
    c.Pc = PcBar * 1.0e5;                  // bar -> Pa
    c.acentric = acentric;
    return c;
}

} // namespace

Component methane() { return makeComponent("methane", 16.043, 190.564, 45.992, 0.0115); }
Component ethane() { return makeComponent("ethane", 30.070, 305.322, 48.722, 0.0995); }
Component propane() { return makeComponent("propane", 44.097, 369.834, 42.477, 0.1523); }
Component iButane() { return makeComponent("iButane", 58.123, 407.817, 36.480, 0.1835); }
Component nButane() { return makeComponent("nButane", 58.123, 425.125, 37.962, 0.2002); }
Component nPentane() { return makeComponent("nPentane", 72.150, 469.700, 33.700, 0.2515); }
Component nHexane() { return makeComponent("nHexane", 86.178, 507.600, 30.250, 0.3013); }
Component carbonDioxide() { return makeComponent("carbonDioxide", 44.010, 304.128, 73.773, 0.2239); }
Component hydrogenSulphide() { return makeComponent("hydrogenSulphide", 34.081, 373.400, 89.628, 0.0900); }
Component nitrogen() { return makeComponent("nitrogen", 28.014, 126.192, 33.958, 0.0372); }

Component componentByName(const std::string& name, bool& found)
{
    found = true;
    if (name == "methane") return methane();
    if (name == "ethane") return ethane();
    if (name == "propane") return propane();
    if (name == "iButane") return iButane();
    if (name == "nButane") return nButane();
    if (name == "nPentane") return nPentane();
    if (name == "nHexane") return nHexane();
    if (name == "carbonDioxide") return carbonDioxide();
    if (name == "hydrogenSulphide") return hydrogenSulphide();
    if (name == "nitrogen") return nitrogen();
    found = false;
    return Component{};
}

Component requireComponent(const std::string& name)
{
    bool found = false;
    const Component c = componentByName(name, found);
    if (!found)
        cfFatal("pvt: unknown component '" + name
                + "' (known: methane, ethane, propane, iButane, nButane, "
                  "nPentane, nHexane, carbonDioxide, hydrogenSulphide, nitrogen)");
    return c;
}

void Mixture::add(const Component& c, scalar moleFraction)
{
    comps_.push_back(c);
    z_.push_back(moleFraction);
}

std::size_t Mixture::index(const std::string& name) const
{
    const auto it = std::find_if(comps_.begin(), comps_.end(),
                                 [&](const Component& c) { return c.name == name; });
    if (it == comps_.end())
        cfFatal("pvt: component '" + name + "' is not in this mixture");
    return static_cast<std::size_t>(std::distance(comps_.begin(), it));
}

void Mixture::setFractions(std::vector<scalar> z)
{
    if (z.size() != comps_.size())
        cfFatal("pvt: composition has " + std::to_string(z.size())
                + " fractions but " + std::to_string(comps_.size()) + " components");
    for (const scalar v : z)
        if (!(v >= 0.0))
            cfFatal("pvt: negative mole fraction in composition");
    z_ = std::move(z);
    normalise();
}

void Mixture::normalise()
{
    const scalar sum = std::accumulate(z_.begin(), z_.end(), 0.0);
    if (sum <= 0.0) cfFatal("pvt: composition sums to zero");
    for (scalar& v : z_) v /= sum;
}

scalar Mixture::molarMass() const
{
    scalar m = 0;
    for (std::size_t i = 0; i < z_.size(); ++i) m += z_[i] * comps_[i].molarMass;
    return m;
}

scalar Mixture::massFraction(std::size_t i) const
{
    const scalar m = molarMass();
    return m > 0 ? z_[i] * comps_[i].molarMass / m : 0.0;
}

std::string Mixture::describe() const
{
    std::ostringstream os;
    os.setf(std::ios::fixed);
    os.precision(6);
    for (std::size_t i = 0; i < comps_.size(); ++i)
        os << comps_[i].name << "=" << z_[i] << " ";
    return os.str();
}

} // namespace cf::pvt
