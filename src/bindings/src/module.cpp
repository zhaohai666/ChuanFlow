// pybind11 exports of the 1D kernels (README section "绑定层 src/bindings"):
// build a network, solve steady or transient, read results in Python.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/1d/Network.h"
#include "cf/1d/Transient.h"

namespace py = pybind11;

PYBIND11_MODULE(chflow_core, m)
{
    m.doc() = "ChuanFlow 1D pipe-network kernel (steady + transient)";

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
}
