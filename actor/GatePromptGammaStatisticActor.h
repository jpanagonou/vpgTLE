/* --------------------------------------------------
 Copyright (C): OpenGATE Collaboration
 This software is distributed under the terms
 of the GNU Lesser General Public Licence (LGPL)
 See LICENSE.md for further details
 -------------------------------------------------- */

#ifndef GatePromptGammaStatisticActor_h
#define GatePromptGammaStatisticActor_h

#include "GateVActor.h"
#include "G4AnalysisManager.hh"
#include "G4HadronicProcessStore.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4Step.hh"
#include "G4Run.hh"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <tuple>
#include <map>
#include <string>

#include "G4Alpha.hh"
#include "G4IonTable.hh"


namespace py = pybind11;

class GatePromptGammaStatisticActor : public GateVActor {

public:

  ~GatePromptGammaStatisticActor() override;

  explicit GatePromptGammaStatisticActor(py::dict &user_info);

  void InitializeUserInfo(py::dict &user_info) override;

  void InitializeCpp() override;

  void BeginOfRunAction(const G4Run *run) override;

  void SteppingAction(G4Step *step) override;

  void EndOfRunAction(const G4Run *run) override;

  void BeginOfRunActionMasterThread(int run_id) override;

  int EndOfRunActionMasterThread(int run_id) override;

private:

  


  // ── Paramètres de binning ─────────────────────────────────────────────────
  G4int    fParticleNbBins;
  G4double fParticleMinEnergy;
  G4double fParticleMaxEnergy;
  G4int    fGammaNbBins;
  G4double fGammaMinEnergy;
  G4double fGammaMaxEnergy;
  

  // ── Paramètres de sortie ──────────────────────────────────────────────────
  G4String    fOutputFilename;
  
  std::string fParticleType;

  // ── Mode mono/multi ───────────────────────────────────────────────────────
  G4bool fMultiElement;
  G4bool fSaveKE0Secondaries;

  // ══════════════════════════════════════════════════════════════════════════
  // MODE MONO
  // ══════════════════════════════════════════════════════════════════════════

  G4String fMaterialName;
  G4String fMonoSymbol;

  // Flags internes
  G4bool   fSigmaFilled;
  G4double fDensity;

  // Shadow arrays
  std::vector<std::vector<G4double>> fGammaZData;
  std::vector<std::vector<G4double>> fEpEpgNormData;   
  std::vector<G4double>              fNInelastic;
  std::vector<std::tuple<G4double, G4double, G4double, G4int>> fGammaEvents;
  std::vector<std::tuple<G4double, G4double, G4double, G4int>> fEpEpgNormEvents; 
  std::vector<G4double> fEpInelasticProducedGamma; 

  // Index des histogrammes mono
  G4int fH2EpEpgMono;
  G4int fH2GammaZMono;
  //G4int fH2EpEpgNormalizedMono;    
  G4int fH1NrPGMono;
  G4int fH1KapaMono;
  G4int fH1EpMono;                 
  G4int fH1EpInelasticMono;        
  G4int fH1EpInelasticProducedGammaMono; 

 
  // This too lignes to investigate about the filter KE>0
  std::vector<G4double> fNInelasticPartiel; 
  std::map<G4String, G4int> fSecondairePartielCount; 
  std::map<G4String, std::vector<G4double>> fEnergiesSecondairesPartiel;

  // Lien avec helium et alpha debut
  std::map<G4String, G4int> fSecondaireTotalCount; 
  std::map<G4String, std::vector<G4double>> fEnergiesSecondairesTotal;
  // Fin





  // ══════════════════════════════════════════════════════════════════════════
  // MODE MULTI
  // ══════════════════════════════════════════════════════════════════════════

  std::map<G4int, std::string> fZtoSymbol;

  // Index des histogrammes par élément
  std::map<G4int, G4int> fH2EpEpgIndex;
  std::map<G4int, G4int> fH2GammaZIndex;
  //std::map<G4int, G4int> fH2EpEpgNormalizedIndex;    
  std::map<G4int, G4int> fH1NrPGIndex;
  std::map<G4int, G4int> fH1KapaIndex;
  std::map<G4int, G4int> fH1EpIndex;                 
  std::map<G4int, G4int> fH1EpInelasticIndex;       
  std::map<G4int, G4int> fH1EpInelasticProducedGammaIndex; 

  // Données par élément
  std::map<G4int, std::vector<std::vector<G4double>>> fGammaZDataMap;
  std::map<G4int, std::vector<std::vector<G4double>>> fEpEpgNormDataMap;  
  std::map<G4int, std::vector<G4double>>              fNInelasticMap;
  std::map<G4int, G4double>                           fDensityMap;
  std::map<G4int, G4bool>                             fSigmaFilledMap;
  std::map<G4int, std::vector<std::tuple<G4double, G4double, G4double, G4int>>> fGammaEventsMap;
  std::map<G4int, std::vector<std::tuple<G4double, G4double, G4double, G4int>>> fEpEpgNormEventsMap; 
  std::map<G4int, std::vector<G4double>> fEpInelasticProducedGammaMap;  

};

#endif // GatePromptGammaStatisticActor_h
