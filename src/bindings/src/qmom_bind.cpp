// QMOM family binding: two-node quadrature-based moment method for
// constant-kernel (Smoluchowski) coagulation. Spatially uniform, no flow
// field -- the moment-closure kernel a DPM/species source framework calls.
//
// Replicates src/3d/apps/qmomApp.cpp: readQmomConfig, runQmom; the coagulation
// source, Gaussian quadrature reconstruction and integration live in the kernel
// Qmom module, this file drives it and copies the history out.
#include <array>
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cf/common/Dictionary.h"
#include "cf/fvm/Qmom.h"

namespace py = pybind11;

namespace
{

struct QmomStepView
{
    int step = 0;
    double t = 0;
    std::array<double, 4> m{0, 0, 0, 0};
    std::array<double, 2> w{0, 0};
    std::array<double, 2> v{0, 0};
    double massError = 0;
};

struct QmomCase
{
    cf::QmomConfig cfg;
    std::vector<QmomStepView> history;
    bool valid = false;
};

std::unique_ptr<QmomCase> solveQmomCase(const std::string& casePath)
{
    const cf::Dictionary caseDict = cf::Dictionary::parseFile(casePath);
    auto out = std::make_unique<QmomCase>();
    out->cfg = cf::readQmomConfig(caseDict);
    if (!out->cfg.enabled)
        throw std::runtime_error("qmom: needs a qmom{} declaration");
    const cf::QmomResult res = cf::runQmom(out->cfg);
    out->valid = res.valid;
    if (!out->valid)
        throw std::runtime_error("qmom: integration failed");
    for (const cf::QmomStep& s : res.history)
    {
        QmomStepView v;
        v.step = s.step;
        v.t = s.t;
        for (int i = 0; i < 4; ++i) v.m[i] = s.m[i];
        for (int i = 0; i < 2; ++i) { v.w[i] = s.w[i]; v.v[i] = s.v[i]; }
        v.massError = s.massError;
        out->history.push_back(v);
    }
    return out;
}

} // namespace

void registerQmomBindings(py::module_& m)
{
    py::class_<QmomStepView>(m, "QmomStep",
        "One reported QMOM state: moments M0..M3, two quadrature nodes/weights "
        "and the mass-conservation error.")
        .def_readonly("step", &QmomStepView::step)
        .def_readonly("t", &QmomStepView::t)
        .def_readonly("moments", &QmomStepView::m)
        .def_readonly("weights", &QmomStepView::w)
        .def_readonly("nodes", &QmomStepView::v)
        .def_readonly("mass_error", &QmomStepView::massError);

    py::class_<QmomCase, std::unique_ptr<QmomCase>>(m, "QmomCase",
        "Result of one QMOM coagulation run: moment/quadrature history and the "
        "case constants.")
        .def_readonly("history", &QmomCase::history)
        .def_readonly("valid", &QmomCase::valid)
        .def_property_readonly("beta",
            [](const QmomCase& r) { return r.cfg.beta; })
        .def_property_readonly("m0_initial",
            [](const QmomCase& r) { return r.cfg.m0; })
        .def_property_readonly("m1_initial",
            [](const QmomCase& r) { return r.cfg.m1; })
        .def_property_readonly("m2_initial",
            [](const QmomCase& r) { return r.cfg.m2; })
        .def_property_readonly("m3_initial",
            [](const QmomCase& r) { return r.cfg.m3; })
        .def_property_readonly("dt",
            [](const QmomCase& r) { return r.cfg.dt; })
        .def_property_readonly("steps",
            [](const QmomCase& r) { return r.cfg.steps; })
        .def("source",
            [](const QmomCase& r, const std::array<double, 4>& moments)
            {
                std::array<double, 4> dm{0, 0, 0, 0};
                cf::coagulationSource(moments.data(), r.cfg.beta, dm.data());
                return dm;
            },
            py::arg("moments"),
            "Constant-kernel Smoluchowski moment source dm/dt.")
        .def("quadrature",
            [](const QmomCase& r, const std::array<double, 4>& moments)
            {
                std::array<double, 2> w{0, 0}, v{0, 0};
                const bool ok = cf::quadratureFromMoments(moments.data(),
                                                           w.data(), v.data());
                py::object none = py::none();
                if (!ok) return py::cast<py::object>(none);
                return py::cast(std::make_pair(w, v));
            },
            py::arg("moments"),
            "Two-node Gaussian quadrature (weights, nodes) for a moment vector; "
            "None when the moments are not a valid positive measure.");

    m.def("solve_qmom_case", &solveQmomCase, py::arg("case_dict_path"),
          "Read a qmom{} coagulation case and advance the four moments with the "
          "two-node QMOM closure.");
}
