// Mixture population-balance family binding: algebraic drift-flux slip plus the
// monodisperse three-moment batch balance. Spatially uniform, no flow field --
// this is the moment-source kernel a DPM/species framework calls.
//
// Replicates src/3d/apps/mixturePbmApp.cpp: readMixturePbmConfig,
// mixtureProperties, runMixturePbm; the ODE source and closure live in the
// kernel MixturePbm module, this file drives it and copies the history out.
#include <array>
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/common/Dictionary.h"
#include "cf/fvm/MixturePbm.h"

namespace py = pybind11;

namespace
{

struct PbmStepView
{
    int step = 0;
    double t = 0;
    std::array<double, 4> m{0, 0, 0, 0};
    double alpha = 0;
    double uG = 0;
    double uL = 0;
};

struct PbmCase
{
    cf::MixturePbmConfig cfg;
    cf::MixturePbmProps props;
    std::vector<PbmStepView> history;
    bool valid = false;
};

std::unique_ptr<PbmCase> solvePbmCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
    auto out = std::make_unique<PbmCase>();
    out->cfg = cf::readMixturePbmConfig(caseDict);
    if (!out->cfg.enabled)
        throw std::runtime_error("mixturePbm: needs a mixturePbm{} declaration");
    out->props = cf::mixtureProperties(out->cfg);
    const cf::MixturePbmResult res = cf::runMixturePbm(out->cfg);
    out->valid = res.valid;
    if (!out->valid)
        throw std::runtime_error("mixturePbm: integration failed");
    for (const cf::MixturePbmStep& s : res.history)
    {
        PbmStepView v;
        v.step = s.step;
        v.t = s.t;
        for (int i = 0; i < 4; ++i) v.m[i] = s.m[i];
        v.alpha = s.alpha;
        v.uG = s.uG;
        v.uL = s.uL;
        out->history.push_back(v);
    }
    return out;
}

} // namespace

void registerPbmBindings(py::module_& m)
{
    py::class_<PbmStepView>(m, "PbmStep",
        "One reported batch state: moments M0..M3, holdup and phase velocities.")
        .def_readonly("step", &PbmStepView::step)
        .def_readonly("t", &PbmStepView::t)
        .def_readonly("moments", &PbmStepView::m)
        .def_readonly("alpha", &PbmStepView::alpha)
        .def_readonly("u_gas", &PbmStepView::uG)
        .def_readonly("u_liquid", &PbmStepView::uL);

    py::class_<PbmCase, std::unique_ptr<PbmCase>>(m, "MixturePbmCase",
        "Result of one mixture PBM batch run: closure constants, moment history, "
        "holdup and phase velocities.")
        .def_readonly("history", &PbmCase::history)
        .def_readonly("valid", &PbmCase::valid)
        .def_property_readonly("v_part",
            [](const PbmCase& r) { return r.props.vPart; })
        .def_property_readonly("holdup_factor",
            [](const PbmCase& r) { return r.props.k; })
        .def_property_readonly("m0_initial",
            [](const PbmCase& r) { return r.props.m0_0; })
        .def_property_readonly("m1_initial",
            [](const PbmCase& r) { return r.props.m1_0; })
        .def_property_readonly("m2_initial",
            [](const PbmCase& r) { return r.props.m2_0; })
        .def_property_readonly("m3_initial",
            [](const PbmCase& r) { return r.props.m3_0; })
        .def_property_readonly("alpha0",
            [](const PbmCase& r) { return r.cfg.alpha0; })
        .def_property_readonly("d0",
            [](const PbmCase& r) { return r.cfg.d0; })
        .def_property_readonly("rho_l",
            [](const PbmCase& r) { return r.cfg.rhoL; })
        .def_property_readonly("rho_g",
            [](const PbmCase& r) { return r.cfg.rhoG; })
        .def_property_readonly("slip",
            [](const PbmCase& r) { return r.cfg.slip; })
        .def_property_readonly("j",
            [](const PbmCase& r) { return r.cfg.j; })
        .def_property_readonly("nucleation",
            [](const PbmCase& r) { return r.cfg.J; })
        .def_property_readonly("breakage",
            [](const PbmCase& r) { return r.cfg.breakage; })
        .def_property_readonly("n_frag",
            [](const PbmCase& r) { return r.cfg.nFrag; })
        .def_property_readonly("dt",
            [](const PbmCase& r) { return r.cfg.dt; })
        .def_property_readonly("steps",
            [](const PbmCase& r) { return r.cfg.steps; })
        .def("source",
            [](const PbmCase& r, const std::array<double, 4>& moments)
            {
                std::array<double, 4> dm{0, 0, 0, 0};
                cf::mixturePbmSource(r.cfg, moments.data(), dm.data());
                return dm;
            },
            py::arg("moments"),
            "Documented monodisperse moment source dm/dt for a moment vector.");

    m.def("solve_mixture_pbm_case", &solvePbmCase, py::arg("case_dict_path"),
          "Read a mixturePbm{} batch case and advance the monodisperse moment "
          "balance with the drift-flux slip closure.");
}
