// pybind11 exports of the 1D kernels and the PVT layer (README section
// "绑定层 src/bindings"): build a network, solve steady or transient, and get
// equation-of-state properties in Python.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/1d/Network.h"
#include "cf/1d/Transient.h"
#include "cf/pvt/Flash.h"
#include "cf/pvt/PvtCase.h"
#include "cf/pvt/PvtTable.h"

namespace py = pybind11;

PYBIND11_MODULE(chflow_core, m)
{
    m.doc() = "ChuanFlow 1D pipe-network kernel (steady + transient) and PVT layer";

    py::class_<cf::FluidProperties>(m, "FluidProperties")
        .def(py::init<>())
        .def(py::init<double, double>(), py::arg("rho"), py::arg("mu"))
        .def_readwrite("rho", &cf::FluidProperties::rho)
        .def_readwrite("mu", &cf::FluidProperties::mu);

    py::class_<cf::NetworkPipe>(m, "NetworkPipe")
        .def(py::init<>())
        .def_readwrite("name", &cf::NetworkPipe::name)
        .def_readwrite("from_name", &cf::NetworkPipe::fromName)
        .def_readwrite("to_name", &cf::NetworkPipe::toName)
        .def_readwrite("length", &cf::NetworkPipe::length)
        .def_readwrite("diameter", &cf::NetworkPipe::diameter)
        .def_readwrite("roughness", &cf::NetworkPipe::roughness)
        .def_readwrite("elevation", &cf::NetworkPipe::elevation)
        .def_readonly("from", &cf::NetworkPipe::from)
        .def_readonly("to", &cf::NetworkPipe::to);

    py::class_<cf::OneDNetwork>(m, "OneDNetwork")
        .def(py::init<>())
        .def_readwrite("fluid", &cf::OneDNetwork::fluid)
        .def("add_node", &cf::OneDNetwork::addNode, py::arg("name"))
        .def("node_index", &cf::OneDNetwork::nodeIndex, py::arg("name"))
        .def("add_pipe", &cf::OneDNetwork::addPipe, py::arg("pipe"))
        .def("set_node_pressure", &cf::OneDNetwork::setNodePressure,
             py::arg("node"), py::arg("pressure"))
        .def("set_node_rate", &cf::OneDNetwork::setNodeRate,
             py::arg("node"), py::arg("rate"))
        .def("n_nodes", &cf::OneDNetwork::nNodes)
        .def_readonly("node_names", &cf::OneDNetwork::nodeNames)
        .def_readonly("pipes", &cf::OneDNetwork::pipes);

    py::class_<cf::SteadyControls>(m, "SteadyControls")
        .def(py::init<>())
        .def_readwrite("max_iter", &cf::SteadyControls::maxIter)
        .def_readwrite("tol", &cf::SteadyControls::tol)
        .def_readwrite("relax", &cf::SteadyControls::relax);

    py::class_<cf::PipeFlowResult>(m, "PipeFlowResult")
        .def_readonly("name", &cf::PipeFlowResult::name)
        .def_readonly("dp_total", &cf::PipeFlowResult::dpTotal)
        .def_readonly("dp_friction", &cf::PipeFlowResult::dpFriction)
        .def_readonly("velocity", &cf::PipeFlowResult::velocity)
        .def_readonly("flow_rate", &cf::PipeFlowResult::flowRate)
        .def_readonly("reynolds", &cf::PipeFlowResult::reynolds)
        .def_readonly("friction_factor", &cf::PipeFlowResult::frictionFactor)
        .def_readonly("conductance", &cf::PipeFlowResult::conductance);

    py::class_<cf::SteadyNetworkResult>(m, "SteadyNetworkResult")
        .def_readonly("node_pressure", &cf::SteadyNetworkResult::nodePressure)
        .def_readonly("pipe_flows", &cf::SteadyNetworkResult::pipeFlows)
        .def_readonly("iterations", &cf::SteadyNetworkResult::iterations)
        .def_readonly("residual", &cf::SteadyNetworkResult::residual)
        .def_readonly("relative_residual", &cf::SteadyNetworkResult::relativeResidual)
        .def_readonly("converged", &cf::SteadyNetworkResult::converged);

    m.def("solve_steady_network", &cf::solveSteadyNetwork,
          py::call_guard<py::gil_scoped_release>(),
          py::arg("network"), py::arg("controls") = cf::SteadyControls{});

    m.def("churchill_friction", &cf::churchillFriction,
          py::arg("reynolds"), py::arg("relative_roughness"));

    m.def("load_network", [](const std::string& caseDictPath)
          { return cf::loadNetworkFromDict(cf::Dictionary::parseFile(caseDictPath)); },
          py::arg("case_dict_path"));
    m.def("load_steady_controls", [](const std::string& caseDictPath)
          { return cf::loadSteadyControls(cf::Dictionary::parseFile(caseDictPath)); },
          py::arg("case_dict_path"));

    // ---- transient kernel ----
    py::class_<cf::TransientNodeConfig>(m, "TransientNodeConfig")
        .def(py::init<>())
        .def_readwrite("name", &cf::TransientNodeConfig::name)
        .def_readwrite("volume", &cf::TransientNodeConfig::volume)
        .def_readwrite("bulk_modulus", &cf::TransientNodeConfig::bulkModulus)
        .def_readwrite("pressure_fixed", &cf::TransientNodeConfig::pressureFixed)
        .def_readwrite("pressure", &cf::TransientNodeConfig::pressure)
        .def_readwrite("liquid_fraction", &cf::TransientNodeConfig::liquidFraction)
        .def_readwrite("hydrate_fraction", &cf::TransientNodeConfig::hydrateFraction)
        .def_readwrite("hydrate_rate", &cf::TransientNodeConfig::hydrateRate);

    py::class_<cf::TransientPipeConfig>(m, "TransientPipeConfig")
        .def(py::init<>())
        .def_readwrite("name", &cf::TransientPipeConfig::name)
        .def_readwrite("from_name", &cf::TransientPipeConfig::fromName)
        .def_readwrite("to_name", &cf::TransientPipeConfig::toName)
        .def_readwrite("length", &cf::TransientPipeConfig::length)
        .def_readwrite("diameter", &cf::TransientPipeConfig::diameter)
        .def_readwrite("roughness", &cf::TransientPipeConfig::roughness)
        .def_readwrite("elevation", &cf::TransientPipeConfig::elevation)
        .def_readwrite("initial_flow", &cf::TransientPipeConfig::initialFlow);

    py::class_<cf::TransientSource>(m, "TransientSource")
        .def(py::init<>())
        .def_readwrite("node_name", &cf::TransientSource::nodeName)
        .def_readwrite("rate", &cf::TransientSource::rate)
        .def_readwrite("fraction", &cf::TransientSource::fraction);

    py::class_<cf::TransientCase>(m, "TransientCase")
        .def(py::init<>())
        .def_readwrite("fluid", &cf::TransientCase::fluid)
        .def_readwrite("nodes", &cf::TransientCase::nodes)
        .def_readwrite("pipes", &cf::TransientCase::pipes)
        .def_readwrite("sources", &cf::TransientCase::sources);

    py::class_<cf::TransientSample>(m, "TransientSample")
        .def_readonly("t", &cf::TransientSample::t)
        .def_readonly("pressure", &cf::TransientSample::pressure)
        .def_readonly("liquid_fraction", &cf::TransientSample::liquidFraction)
        .def_readonly("hydrate_fraction", &cf::TransientSample::hydrateFraction)
        .def_readonly("flow_rate", &cf::TransientSample::flowRate);

    py::class_<cf::TransientSeries>(m, "TransientSeries")
        .def_readonly("node_names", &cf::TransientSeries::nodeNames)
        .def_readonly("pipe_names", &cf::TransientSeries::pipeNames)
        .def_readonly("samples", &cf::TransientSeries::samples);

    py::class_<cf::TransientTimeControls>(m, "TransientTimeControls")
        .def(py::init<>())
        .def_readwrite("dt", &cf::TransientTimeControls::dt)
        .def_readwrite("end_time", &cf::TransientTimeControls::endTime)
        .def_readwrite("report_every", &cf::TransientTimeControls::reportEvery);

    py::class_<cf::TransientNetwork>(m, "TransientNetwork")
        .def(py::init<const cf::TransientCase&>(), py::arg("case"))
        .def("step", &cf::TransientNetwork::step, py::arg("dt"),
             py::call_guard<py::gil_scoped_release>())
        .def("time", &cf::TransientNetwork::time)
        .def("node_names", &cf::TransientNetwork::nodeNames)
        .def("pipe_names", &cf::TransientNetwork::pipeNames)
        .def("pressure", &cf::TransientNetwork::pressure)
        .def("pipe_flows", &cf::TransientNetwork::pipeFlows)
        .def("liquid_fractions", &cf::TransientNetwork::liquidFractions)
        .def("hydrate_fractions", &cf::TransientNetwork::hydrateFractions)
        .def("liquid_volumes", &cf::TransientNetwork::liquidVolumes)
        .def("hydrate_volumes", &cf::TransientNetwork::hydrateVolumes);

    m.def("hydrate_kinetics_step",
          [](double liquidVolume, double hydrateVolume, double rate, double dt)
          { double vL = liquidVolume, vH = hydrateVolume;
            cf::hydrateKineticsStep(vL, vH, rate, dt);
            return std::make_pair(vL, vH); },
          py::arg("liquid_volume"), py::arg("hydrate_volume"),
          py::arg("rate"), py::arg("dt"));

    m.def("run_transient", &cf::runTransient,
          py::call_guard<py::gil_scoped_release>(),
          py::arg("case"), py::arg("dt"), py::arg("end_time"),
          py::arg("report_every") = 1);

    m.def("load_transient_case", [](const std::string& caseDictPath)
          { return cf::loadTransientFromDict(cf::Dictionary::parseFile(caseDictPath)); },
          py::arg("case_dict_path"));
    m.def("load_transient_controls", [](const std::string& caseDictPath)
          { return cf::loadTransientControls(cf::Dictionary::parseFile(caseDictPath)); },
          py::arg("case_dict_path"));

    // ---- PVT layer: Peng-Robinson EOS, TP flash, tabulated properties ----
    py::class_<cf::pvt::Component>(m, "Component")
        .def(py::init<>())
        .def_readwrite("name", &cf::pvt::Component::name)
        .def_readwrite("molar_mass", &cf::pvt::Component::molarMass)
        .def_readwrite("tc", &cf::pvt::Component::Tc)
        .def_readwrite("pc", &cf::pvt::Component::Pc)
        .def_readwrite("acentric", &cf::pvt::Component::acentric);

    py::class_<cf::pvt::Mixture>(m, "Mixture")
        .def(py::init<>())
        .def("add", &cf::pvt::Mixture::add, py::arg("component"), py::arg("mole_fraction"))
        .def("index", &cf::pvt::Mixture::index, py::arg("name"))
        .def("components", [](const cf::pvt::Mixture& self) { return self.components(); })
        .def("fractions", [](const cf::pvt::Mixture& self) { return self.fractions(); })
        .def("set_fractions", &cf::pvt::Mixture::setFractions, py::arg("fractions"))
        .def("normalise", &cf::pvt::Mixture::normalise)
        .def("molar_mass", &cf::pvt::Mixture::molarMass)
        .def("mass_fraction", &cf::pvt::Mixture::massFraction, py::arg("index"))
        .def("describe", &cf::pvt::Mixture::describe)
        .def("size", &cf::pvt::Mixture::size);

    py::enum_<cf::pvt::RootPreference>(m, "RootPreference")
        .value("lowest_gibbs", cf::pvt::RootPreference::LowestGibbs)
        .value("liquid", cf::pvt::RootPreference::Liquid)
        .value("vapor", cf::pvt::RootPreference::Vapor);

    py::class_<cf::pvt::EosState>(m, "EosState")
        .def_readonly("z", &cf::pvt::EosState::Z)
        .def_readonly("molar_volume", &cf::pvt::EosState::molarVolume)
        .def_readonly("number_density_molar", &cf::pvt::EosState::numberDensityMolar)
        .def_readonly("mass_density", &cf::pvt::EosState::massDensity)
        .def_readonly("molar_volume_over_b", &cf::pvt::EosState::molarVolumeOverB)
        .def_readonly("ln_phi", &cf::pvt::EosState::lnPhi)
        .def_readonly("gibbs", &cf::pvt::EosState::gibbs)
        .def_readonly("had_two_roots", &cf::pvt::EosState::hadTwoRoots);

    py::class_<cf::pvt::FlashControls>(m, "FlashControls")
        .def(py::init<>())
        .def_readwrite("tol", &cf::pvt::FlashControls::tol)
        .def_readwrite("max_iter", &cf::pvt::FlashControls::maxIter);

    py::class_<cf::pvt::FlashResult>(m, "FlashResult")
        .def_readonly("converged", &cf::pvt::FlashResult::converged)
        .def_readonly("iterations", &cf::pvt::FlashResult::iterations)
        .def_readonly("vapor_fraction", &cf::pvt::FlashResult::vaporFraction)
        .def_readonly("liquid_composition", &cf::pvt::FlashResult::liquidComposition)
        .def_readonly("vapour_composition", &cf::pvt::FlashResult::vapourComposition)
        .def_readonly("liquid", &cf::pvt::FlashResult::liquid)
        .def_readonly("vapour", &cf::pvt::FlashResult::vapour)
        .def_readonly("liquid_mass_density", &cf::pvt::FlashResult::liquidMassDensity)
        .def_readonly("vapour_mass_density", &cf::pvt::FlashResult::vapourMassDensity)
        .def_readonly("liquid_molar_mass", &cf::pvt::FlashResult::liquidMolarMass)
        .def_readonly("vapour_molar_mass", &cf::pvt::FlashResult::vapourMolarMass)
        .def_readonly("single_phase", &cf::pvt::FlashResult::singlePhase)
        .def_readonly("residual", &cf::pvt::FlashResult::residual)
        .def_readonly("message", &cf::pvt::FlashResult::message);

    py::class_<cf::pvt::BilinearTable>(m, "BilinearTable")
        .def(py::init<std::vector<double>, std::vector<double>,
                      std::vector<double>>(),
             py::arg("p_grid"), py::arg("t_grid"), py::arg("values"))
        .def("sample", &cf::pvt::BilinearTable::sample, py::arg("p"), py::arg("t"))
        .def_property_readonly("p_grid", &cf::pvt::BilinearTable::xGrid)
        .def_property_readonly("t_grid", &cf::pvt::BilinearTable::yGrid)
        .def("size", &cf::pvt::BilinearTable::size);

    m.def("component", &cf::pvt::requireComponent, py::arg("name"));
    m.def("wilson_k", &cf::pvt::wilsonK, py::arg("component"),
          py::arg("t"), py::arg("p"));
    m.def("pure_vapor_pressure", &cf::pvt::pureVaporPressure, py::arg("component"),
          py::arg("t"));
    m.def("peng_robinson_state", &cf::pvt::pengRobinsonState, py::arg("mixture"),
          py::arg("t"), py::arg("p"),
          py::arg("pref") = cf::pvt::RootPreference::LowestGibbs);
    m.def("tp_flash",
          [](const cf::pvt::Mixture& feed, double t, double p,
             const cf::pvt::FlashControls& controls)
          { return cf::pvt::tpFlash(feed, t, p, controls); },
          py::call_guard<py::gil_scoped_release>(),
          py::arg("feed"), py::arg("t"), py::arg("p"),
          py::arg("controls") = cf::pvt::FlashControls{});
    m.def("bubble_pressure", &cf::pvt::bubblePressure,
          py::call_guard<py::gil_scoped_release>(), py::arg("feed"), py::arg("t"));
    m.def("liquid_density_table", &cf::pvt::liquidDensityTable,
          py::call_guard<py::gil_scoped_release>(),
          py::arg("feed"), py::arg("p_grid"), py::arg("t_grid"));
    m.def("vapour_density_table", &cf::pvt::vapourDensityTable,
          py::call_guard<py::gil_scoped_release>(),
          py::arg("feed"), py::arg("p_grid"), py::arg("t_grid"));
    m.def("vapour_fraction_table", &cf::pvt::vapourFractionTable,
          py::call_guard<py::gil_scoped_release>(),
          py::arg("feed"), py::arg("p_grid"), py::arg("t_grid"));

    py::class_<cf::pvt::FlashCase>(m, "FlashCase")
        .def(py::init<>())
        .def_readwrite("feed", &cf::pvt::FlashCase::feed)
        .def_readwrite("temperature", &cf::pvt::FlashCase::temperature)
        .def_readwrite("pressures", &cf::pvt::FlashCase::pressures)
        .def_readwrite("controls", &cf::pvt::FlashCase::controls)
        .def_readwrite("compute_bubble_point", &cf::pvt::FlashCase::computeBubblePoint);

    m.def("load_flash_case", [](const std::string& caseDictPath)
          { return cf::pvt::loadFlashCase(cf::Dictionary::parseFile(caseDictPath)); },
          py::arg("case_dict_path"));
}
