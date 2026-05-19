/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "GatePromptGammaStatisticActor.h"

class PyGatePromptGammaStatisticActor
    : public GatePromptGammaStatisticActor {
public:
    // Hérite des constructeurs
    using GatePromptGammaStatisticActor::GatePromptGammaStatisticActor;

    void BeginOfRunActionMasterThread(int run_id) override {
        PYBIND11_OVERLOAD(void, GatePromptGammaStatisticActor,
                          BeginOfRunActionMasterThread, run_id);
    }

    int EndOfRunActionMasterThread(int run_id) override {
        PYBIND11_OVERLOAD(int, GatePromptGammaStatisticActor,
                          EndOfRunActionMasterThread, run_id);
    }
};

void init_GatePromptGammaStatisticActor(py::module &m) {
    py::class_<GatePromptGammaStatisticActor,
               PyGatePromptGammaStatisticActor,
               std::unique_ptr<GatePromptGammaStatisticActor, py::nodelete>,
               GateVActor>(m, "GatePromptGammaStatisticActor")
        .def(py::init<py::dict &>())
        .def("InitializeUserInfo",
             &GatePromptGammaStatisticActor::InitializeUserInfo)
        .def("BeginOfRunActionMasterThread",
             &GatePromptGammaStatisticActor::BeginOfRunActionMasterThread)
        .def("EndOfRunActionMasterThread",
             &GatePromptGammaStatisticActor::EndOfRunActionMasterThread);
}
