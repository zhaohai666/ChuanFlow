// pybind11 exports of the 1D steady network kernel (README section
// "绑定层 src/bindings"): build a network, solve, read results in Python.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/1d/Network.h"

namespace py = pybind11;

PYBIND11_MODULE(chflow_core, m)
{
    m.doc() = "ChuanFlow 1D steady pipe-network kernel";

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
}
