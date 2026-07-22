/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateSPSEneDistributionVoxelizedPG.h"
#include <Randomize.hh>
#include <stdexcept>

GateSPSEneDistributionVoxelizedPG::GateSPSEneDistributionVoxelizedPG()
    : fVoxelPositionGenerator(nullptr), fEmin(0), fEmax(0), fNbBins(0),
      fBinWidth(0), fSizeX(0), fSizeY(0), fSizeZ(0) {}

// Stocke les pointeurs recu
// etabli couplage entre generateur de position et generateur d'energie
void GateSPSEneDistributionVoxelizedPG::SetVoxelPositionGenerator(
    GateSPSVoxelsPosDistribution *pg) {
  fVoxelPositionGenerator = pg;
}


void GateSPSEneDistributionVoxelizedPG::InitializeEnergyCDF(
    py::array_t<double> array_4d, double Emin, double Emax, int nbBins) {
      //Enregistre les parametres de l'axe energie et precalcule la largeur d'un bin
      // (utilise plus tard dans VGenerateOne pour convertir un indice de bin en valeur energie continue)
  fEmin = Emin;
  fEmax = Emax;
  fNbBins = nbBins;
  fBinWidth = (fEmax - fEmin) / fNbBins;
      // PLUS RAPIDE, Acces au tableau numppy sans copie (API buffer de pybind), buf.ptr oint directement vers les donnees python
      // buf.shape donne les dimensions du tableau [nbs, Z, Y, X]
  auto buf = array_4d.request();
  if (buf.ndim != 4) {
    throw std::runtime_error(
        "InitializeEnergyCDF attend un tableau numpy 4D de forme "
        "[nbBins, Z, Y, X]");
  }

  const int nBins = static_cast<int>(buf.shape[0]);
  const int nZ = static_cast<int>(buf.shape[1]);
  const int nY = static_cast<int>(buf.shape[2]);
  const int nX = static_cast<int>(buf.shape[3]);
  fSizeX = nX;
  fSizeY = nY;
  fSizeZ = nZ;

  const auto *data = static_cast<double *>(buf.ptr);

  fVoxelEnergyCDF.clear();

  // triple boucle : un voxel a la fois, visite de chaque voxel
  for (int iz = 0; iz < nZ; ++iz) {
    for (int iy = 0; iy < nY; ++iy) {
      for (int ix = 0; ix < nX; ++ix) {

        // POur ce voxel, additionne son spectre complet (tous les bins 
        // denergie) pour obtenir le yield total du voxel
        double total = 0.0;
        for (int ib = 0; ib < nBins; ++ib) {
          const size_t idx = ((static_cast<size_t>(ib) * nZ + iz) *
                                  static_cast<size_t>(nY) +
                              iy) *
                                 static_cast<size_t>(nX) +
                             ix;
          total += data[idx];
        }

        // voxel vide -> on ne stocke rien 
        if (total <= 0.0)
          continue;

        // construction de la CDF cumulee et normalisee (valeur de 0 a 1) pour ce voxel 
        std::vector<double> cdf(nBins);
        double cumul = 0.0;
        for (int ib = 0; ib < nBins; ++ib) {
          const size_t idx = ((static_cast<size_t>(ib) * nZ + iz) *
                                  static_cast<size_t>(nY) +
                              iy) *
                                 static_cast<size_t>(nX) +
                             ix;
          cumul += data[idx];
          cdf[ib] = cumul / total;
        }

          // Calcul un identifiant unique pour ce voxel (aplatit les 3 coordonnes iz, iy, ix en un seul indice)
          // et stocke le sa CDF dans un map, prete a etre utiliser par VGenerateOne
        const size_t linear_index =
            (static_cast<size_t>(iz) * nY + iy) * static_cast<size_t>(nX) +
            ix;
        fVoxelEnergyCDF[linear_index] = std::move(cdf);
      }
    }
  }
}

G4double
GateSPSEneDistributionVoxelizedPG::VGenerateOne(G4ParticleDefinition *) {
  // recupere l'indice du voxel tire juste avant par le generateur de
  // position (meme particule primaire )
  int ix, iy, iz;
  fVoxelPositionGenerator->GetLastVoxelIndices(ix, iy, iz);

  // Recalcule le mm indice lineaire comme lors du stockage
  const size_t linear_index =
      (static_cast<size_t>(iz) * fSizeY + iy) * static_cast<size_t>(fSizeX) +
      ix;

  // recherche la CDF correspondant a ce voxel dans le stockage
  const auto it = fVoxelEnergyCDF.find(linear_index);
  if (it == fVoxelEnergyCDF.end()) {
    //securite : ne devrait normalement JAMAIS arriver, car la CDF spatiale 
    // ne peut tirer qu'un voxel non vide; Si ca arrive quand meme 
    // (bug ou incoherence entre CDF spatiale energie), 
    // il y aura avertissement une seule fois pour ne pas poluer  la sortie si le cas se repette,
    // Puis renvoie l'energie minimale plutot que de planter la simulation;et le voxel tire par la CDF
    static bool warned = false;
    if(warned) {

      G4cerr << "ATTENTION: voxel sans CDF energie trouve (incoherence " 
                "position/energie), utilisation de Emin par defaut. "
                 "Ce message ne s'affiche qu'une fois." 
             << G4endl; 
      warned = true; 
    }

    return fEmin;
  }
  const std::vector<double> &cdf = it->second;

    // tirage par inversion de CDF : meme principe que celui deja utilise
   //  dans GateSPSVoxelsPosDistribution::VGenerateOne() pour la position.
   //  On tire un nombre aleatoire uniforme entre 0 et 1, puis on cherche 
   //le premier bin dont la CDF cumulee depasse ce nombre 

  const double p = G4UniformRand();
  const auto lower = std::lower_bound(cdf.begin(), cdf.end(), p);
  int ib = std::distance(cdf.begin(), lower);
  if (ib >= fNbBins)
    ib = fNbBins - 1;

  // convertit l'indice de bin discret en une valeur d'energie continue, 
  // avec un tirage uniforme supplementaire a l'interieur du bin (evite 
  // que l'energie ne soit toujours exactement au meme point du bin) 

  const double energy = fEmin + (ib + G4UniformRand()) * fBinWidth;

  return energy;
}
