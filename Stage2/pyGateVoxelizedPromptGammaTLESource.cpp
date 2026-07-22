/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "GateVoxelizedPromptGammaTLESource.h"
#include "GateVoxelSource.h"

void init_GateVoxelizedPromptGammaTLESource(py::module &m) {

  // hérite du binding de GateVoxelSource : GetSPSVoxelPosDistribution et
  // InitializeUserInfo sont donc automatiquement exposés côté Python
  py::class_<GateVoxelizedPromptGammaTLESource, GateVoxelSource>(
      m, "GateVoxelizedPromptGammaTLESource")
      .def(py::init())
      .def("InitializeEnergyCDF",
           &GateVoxelizedPromptGammaTLESource::InitializeEnergyCDF);
}
