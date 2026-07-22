/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateVoxelizedPromptGammaTLESource_h
#define GateVoxelizedPromptGammaTLESource_h

#include "GateSPSEneDistributionVoxelizedPG.h"
#include "GateVoxelSource.h"
#include <pybind11/numpy.h>

namespace py = pybind11;

/*
    GateVoxelizedPromptGammaTLESource = distribution 4D des spectres PG.

    Herite de GateVoxelSource pour reutiliser telle quelle toute la partie
    position (fVoxelPositionGenerator, PrepareNextRun, InitializePosition,
    GetSPSVoxelPosDistribution). N'ajoute que la partie energie :
    un generateur GateSPSEneDistributionVoxelizedPG, relie au generateur
    de position herite pour tirer l'energie dans le voxel correspondant.
*/
class GateVoxelizedPromptGammaTLESource : public GateVoxelSource {

public:
  GateVoxelizedPromptGammaTLESource();

  ~GateVoxelizedPromptGammaTLESource() override;

  // transmet le tableau numpy 4D d'energie a fVoxelEnergyGenerator
  // (appelee depuis Python dans update_pg_images)
  void InitializeEnergyCDF(py::array_t<double> array_4d, double Emin,
                           double Emax, int nbBins);

protected:
  // injecte fVoxelEnergyGenerator comme generateur d'energie de la
  // source, a la place du generateur standard (meme principe que
  // GateVoxelSource::InitializePosition pour la position)
  void InitializeEnergy(py::dict user_info) override;

  GateSPSEneDistributionVoxelizedPG *fVoxelEnergyGenerator;
};

#endif // GateVoxelizedPromptGammaTLESource_h
