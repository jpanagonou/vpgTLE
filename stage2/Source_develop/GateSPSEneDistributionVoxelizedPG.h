/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateSPSEneDistributionVoxelizedPG_h
#define GateSPSEneDistributionVoxelizedPG_h

#include "GateSPSEneDistribution.h"
#include "GateSPSVoxelsPosDistribution.h"
#include <pybind11/numpy.h>
#include <unordered_map>

namespace py = pybind11;

/*
    Generateur d'energie voxelise pour les sources de type prompt-gamma (PG).

    Contrairement a GateSPSEneDistribution standard (une seule CDF en energie
    globale), cette classe stocke une CDF en energie DIFFeRENTE PAR VOXEL.

    Le voxel utilise pour le tirage n'est pas choisi ici : il est lu depuis
    le generateur de position (GateSPSVoxelsPosDistribution), qui memorise
    l'indice du dernier voxel tire a chaque appel de sa methode
    VGenerateOne(). Cela garantit la coherence spatiale entre la position
    et l'energie tirees pour une meme particule primaire.

    Stockage sparse (std::unordered_map) : seuls les voxels avec un spectre
    non vide sont stockes, ce qui est important car la plupart des voxels
    d'une image de dose/PG sont vides.
*/
class GateSPSEneDistributionVoxelizedPG : public GateSPSEneDistribution {

public:
  GateSPSEneDistributionVoxelizedPG();

  ~GateSPSEneDistributionVoxelizedPG() override = default;

  // lien vers le generateur de position, pour connaitre le voxel tire
  void SetVoxelPositionGenerator(GateSPSVoxelsPosDistribution *pg);

  // construit la CDF en energie pour chaque voxel non vide
  // (option 3 : triple boucle en C++), a partir du tableau numpy 4D
  // de forme [nbBins, Z, Y, X]
  void InitializeEnergyCDF(py::array_t<double> array_4d, double Emin,
                           double Emax, int nbBins);

  // tire l'energie dans la CDF du voxel courant (lu depuis
  // fVoxelPositionGenerator)
  G4double VGenerateOne(G4ParticleDefinition *pdef) override;

protected:
  GateSPSVoxelsPosDistribution *fVoxelPositionGenerator;

  double fEmin;
  double fEmax;
  int fNbBins;
  double fBinWidth;

  // dimensions spatiales de l'image (pour le calcul de l'indice lineaire)
  int fSizeX;
  int fSizeY;
  int fSizeZ;

  // CDF en energie par voxel (sparse)
  // cle = indice lineaire du voxel : iz * (sizeY * sizeX) + iy * sizeX + ix
  std::unordered_map<size_t, std::vector<double>> fVoxelEnergyCDF;
};

#endif // GateSPSEneDistributionVoxelizedPG_h
