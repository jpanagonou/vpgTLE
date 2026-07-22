/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateVoxelizedPromptGammaTLESource.h"

GateVoxelizedPromptGammaTLESource::GateVoxelizedPromptGammaTLESource()
    : GateVoxelSource() {
  fVoxelEnergyGenerator = new GateSPSEneDistributionVoxelizedPG();
  // relie le generateur d'energie au generateur de position deje cree
  // par le constructeur de GateVoxelSource (fVoxelPositionGenerator)
  fVoxelEnergyGenerator->SetVoxelPositionGenerator(fVoxelPositionGenerator);
}

GateVoxelizedPromptGammaTLESource::~GateVoxelizedPromptGammaTLESource() =
    default;

void GateVoxelizedPromptGammaTLESource::InitializeEnergy(
    py::dict /*user_info*/) {
  // injecte le generateur d'energie voxelise dans la source du thread
  // courant (meme mecanisme que GateVoxelSource::InitializePosition)
  auto &ll = GetThreadLocalDataGenericSource();
  ll.fSPS->SetEneGenerator(fVoxelEnergyGenerator);
}

void GateVoxelizedPromptGammaTLESource::InitializeEnergyCDF(
    py::array_t<double> array_4d, double Emin, double Emax, int nbBins) {
  fVoxelEnergyGenerator->InitializeEnergyCDF(array_4d, Emin, Emax, nbBins);
}
