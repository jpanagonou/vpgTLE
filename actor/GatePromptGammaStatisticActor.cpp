/* --------------------------------------------------
 Copyright (C): OpenGATE Collaboration
 This software is distributed under the terms
 of the GNU Lesser General Public Licence (LGPL)
 See LICENSE.md for further details
 -------------------------------------------------- */

#include "GatePromptGammaStatisticActor.h"
#include "GateHelpersDict.h"

#include "G4AnalysisManager.hh"
#include "G4Gamma.hh"
#include "G4HadronicProcessStore.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4HadronicProcess.hh"
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
GatePromptGammaStatisticActor::GatePromptGammaStatisticActor(py::dict &user_info)
  : GateVActor(user_info, true) {
  fMultiThreadReady = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
GatePromptGammaStatisticActor::~GatePromptGammaStatisticActor() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Reading the parameters define in python script 
void GatePromptGammaStatisticActor::InitializeUserInfo(py::dict &user_info) {
  GateVActor::InitializeUserInfo(user_info);

  fParticleNbBins    = py::int_(user_info["particleNbBins"]);
  fParticleMinEnergy = py::float_(user_info["particleMinEnergy"]);
  fParticleMaxEnergy = py::float_(user_info["particleMaxEnergy"]);

  fGammaNbBins     = py::int_(user_info["gammaNbBins"]);
  fGammaMinEnergy  = py::float_(user_info["gammaMinEnergy"]);
  fGammaMaxEnergy  = py::float_(user_info["gammaMaxEnergy"]);

  fOutputFilename  = py::str(user_info["pg_output_filename"]);
  fParticleType    = py::str(user_info["particle_type"]);
  fMultiElement    = py::bool_(user_info["multi_element"]);
  fSaveKE0Secondaries = py::bool_(user_info["save_KE0_secondaries"]);
 
   // In mono mode, we pass the name of the material,
   //  don't pass any material name in multi-elements mode 
  if (!fMultiElement)
     fMaterialName = py::str(user_info["material_name"]); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GatePromptGammaStatisticActor::InitializeCpp() {
  GateVActor::InitializeCpp();

 
  if (!fMultiElement) {
   // Mono mode : initialization of vectors
     fSigmaFilled = false;
     fNInelastic.assign(fParticleNbBins, 0.0);
     fDensity = 0.0;
     fGammaEvents.clear();
     fEpEpgNormEvents.clear();
     fEpInelasticProducedGamma.clear();
     fGammaZData.assign(fParticleNbBins,
                        std::vector<G4double>(fGammaNbBins, 0.0));
     fEpEpgNormData.assign(fParticleNbBins,
                           std::vector<G4double>(fGammaNbBins, 0.0));

   /* This two lignes to do tecking about the filter KE> 0, this is just avilabe for mono mode*/
   fNInelasticPartiel.assign(fParticleNbBins, 0.0); 
   fSecondairePartielCount.clear(); 
   fEnergiesSecondairesPartiel.clear(); 
   fSecondaireTotalCount.clear();
   fEnergiesSecondairesTotal.clear();



  } else {
     fGammaZDataMap.clear();
     fEpEpgNormDataMap.clear();
     fNInelasticMap.clear();
     fDensityMap.clear();
     fGammaEventsMap.clear();
     fEpEpgNormEventsMap.clear();
     fEpInelasticProducedGammaMap.clear();
     fSigmaFilledMap.clear();
     fZtoSymbol.clear();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GatePromptGammaStatisticActor::BeginOfRunActionMasterThread(int run_id) {

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetFileName(fOutputFilename);
  analysisManager->SetVerboseLevel(1);



   G4String h2Name   = (fParticleType == "neutron") ? "EnEpg"  :
                    (fParticleType == "helium")   ? "EHeEpg" :
                    (fParticleType == "carbon")   ? "ECEpg"  :
                                                    "EpEpg";

   G4String axisName = (fParticleType == "neutron") ? "Neutron energy [MeV]" :
                    (fParticleType == "helium")   ? "Helium energy [MeV]"  :
                    (fParticleType == "carbon")   ? "Carbon energy [MeV]"  :
                                                    "Proton energy [MeV]";

   G4String h1EpName = (fParticleType == "neutron") ? "En"  :
                    (fParticleType == "helium")   ? "EHe" :
                    (fParticleType == "carbon")   ? "EC"  :
                                                    "Ep";




  // ════════════════════════════════════════════════════════════════════════
  // MODE MONO
  // ════════════════════════════════════════════════════════════════════════
  if (!fMultiElement) {

     G4Material *mat = G4Material::GetMaterial(fMaterialName, false);
     if (!mat) {
        G4Exception("GatePromptGammaStatisticActor", "MaterialNotFound",
                    FatalException,
                    ("Material '" + fMaterialName + "' not found.").c_str());
     }

     fMonoSymbol       = mat->GetElement(0)->GetSymbol();
     G4String elemName = mat->GetElement(0)->GetName();

     G4cout << "[Mono] Material: " << fMaterialName
            << " → symbol: " << fMonoSymbol
            << " → name: "   << elemName << G4endl;

     fH2EpEpgMono = analysisManager->CreateH2(
        elemName + "/" + h2Name,
        "Prompt gamma spectrum;"
        "Particle energy [MeV];PG energy [MeV];counts",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV,
        fGammaNbBins,  fGammaMinEnergy  / MeV, fGammaMaxEnergy  / MeV);

     fH2GammaZMono = analysisManager->CreateH2(
        elemName + "/GammaZ",
        "PG yield weighted by kappa_inel;"
        "Particle energy [MeV];PG energy [MeV];kappa [cm^{-1} per collision]",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV,
        fGammaNbBins,  fGammaMinEnergy  / MeV, fGammaMaxEnergy  / MeV);

       

     fH1NrPGMono = analysisManager->CreateH1(
        elemName + "/NrPG",
        "Number of prompt gammas per bin;"
        "Particle energy [MeV];N_{PG}",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV);

     fH1KapaMono = analysisManager->CreateH1(
        elemName + "/Kapa inelastique",
        "Linear attenuation coefficient;"
        "Particle energy [MeV];kappa_{inel} [cm^{-1}]",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV);


     fH1EpMono = analysisManager->CreateH1(
        elemName + "/" + h1EpName,
        "Incident particle energy spectrum;"
        "Particle energy [MeV];counts",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV);

     fH1EpInelasticMono = analysisManager->CreateH1(
        elemName + "/" + h1EpName + "Inelastic",
        "Number of inelastic interactions;"
        "Particle energy [MeV];counts",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV);

     fH1EpInelasticProducedGammaMono = analysisManager->CreateH1(
        elemName + "/" + h1EpName + "InelasticProducedGamma",
        "Inelastic interactions producing at least one PG;"
        "Particle energy [MeV];counts",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV);



     for (auto &row : fGammaZData)
        std::fill(row.begin(), row.end(), 0.0);
     for (auto &row : fEpEpgNormData)
        std::fill(row.begin(), row.end(), 0.0);
     fSigmaFilled = false;

     analysisManager->OpenFile();
     return;
  }

  // ════════════════════════════════════════════════════════════════════════
  // MODE MULTI
  // ════════════════════════════════════════════════════════════════════════

  G4Material *allElemMat = G4Material::GetMaterial("AllElements", false);
  if (!allElemMat) {
     G4Exception("GatePromptGammaStatisticActor",
                 "MultiElementMaterialNotFound", FatalException,
                 "Material 'AllElements' not found. "
                 "Load it via add_material_database() in Python.");
  }

  fZtoSymbol.clear();
  const G4ElementVector *elements = allElemMat->GetElementVector();
  for (const G4Element *elem : *elements) {
     G4int    Z   = static_cast<G4int>(elem->GetZ());
     G4String sym = elem->GetSymbol();
     fZtoSymbol[Z] = sym;
  }

  G4cout << "[MultiElement] " << fZtoSymbol.size()
         << " elements loaded from AllElements: ";
  for (auto &[Z, sym] : fZtoSymbol)
     G4cout << sym << "(Z=" << Z << ") ";
  G4cout << G4endl;

  for (auto &[Z, sym] : fZtoSymbol) {

     fGammaZDataMap[Z].assign(fParticleNbBins,
                              std::vector<G4double>(fGammaNbBins, 0.0));
     fEpEpgNormDataMap[Z].assign(fParticleNbBins,
                                 std::vector<G4double>(fGammaNbBins, 0.0));
     fNInelasticMap[Z].assign(fParticleNbBins, 0.0);
     fDensityMap[Z]     = 1.0;
     fSigmaFilledMap[Z] = false;
     fGammaEventsMap[Z].clear();
     fEpEpgNormEventsMap[Z].clear();
     fEpInelasticProducedGammaMap[Z].clear();

     G4Element *elemPtr = G4Element::GetElement(sym, false);
     G4String elemName  = elemPtr ? elemPtr->GetName() : G4String(sym);

     fH2EpEpgIndex[Z] = analysisManager->CreateH2(
        elemName + "/" + h2Name,
        "Prompt gamma spectrum;"
        "Particle energy [MeV];PG energy [MeV];counts",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV,
        fGammaNbBins,  fGammaMinEnergy  / MeV, fGammaMaxEnergy  / MeV);

     fH2GammaZIndex[Z] = analysisManager->CreateH2(
        elemName + "/GammaZ",
        "PG yield weighted by kappa_inel;"
        "Particle energy [MeV];PG energy [MeV];kappa",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV,
        fGammaNbBins,  fGammaMinEnergy  / MeV, fGammaMaxEnergy  / MeV);

   

     fH1NrPGIndex[Z] = analysisManager->CreateH1(
        elemName + "/NrPG",
        "Number of prompt gammas per bin;"
        "Particle energy [MeV];N",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV);

     fH1KapaIndex[Z] = analysisManager->CreateH1(
        elemName + "/Kapa inelastique",
        "Linear attenuation coefficient;"
        "Particle energy [MeV];kappa",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV);


     fH1EpIndex[Z] = analysisManager->CreateH1(
        elemName + "/" + h1EpName,
        "Incident particle energy spectrum;"
        "Particle energy [MeV];counts",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV);

     fH1EpInelasticIndex[Z] = analysisManager->CreateH1(
        elemName + "/" + h1EpName + "Inelastic",
        "Number of inelastic interactions;"
        "Particle energy [MeV];counts",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV);

     fH1EpInelasticProducedGammaIndex[Z] = analysisManager->CreateH1(
        elemName + "/" + h1EpName + "InelasticProducedGamma",
        "Inelastic interactions producing at least one PG;"
        "Particle energy [MeV];counts",
        fParticleNbBins, fParticleMinEnergy / MeV, fParticleMaxEnergy / MeV);
  }

  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Nothing to do here Initialization is handled int BeginOfRunActionMasterThread
void GatePromptGammaStatisticActor::BeginOfRunAction(const G4Run *run) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GatePromptGammaStatisticActor::SteppingAction(G4Step *step) {

   // Widh of one proton and gamma energy bin, the third is the energy at the center of the first bin of proton
  const G4double dE          = (fParticleMaxEnergy - fParticleMinEnergy) / fParticleNbBins;
  const G4double dEgamma     = (fGammaMaxEnergy  - fGammaMinEnergy)  / fGammaNbBins;
  const G4double energyStart = fParticleMinEnergy + 0.5 * dE;


  // ───────── Particle type filter, keep only protons and neutrons ────────────────
  const G4ParticleDefinition *particle = step->GetTrack()->GetParticleDefinition();
      if (fParticleType == "proton"  && particle != G4Proton::Proton())   return;
      if (fParticleType == "neutron" && particle != G4Neutron::Neutron())  return;
      if (fParticleType == "helium"  && particle != G4Alpha::Alpha())      return;
      if (fParticleType == "carbon"  && 
             !(particle->GetAtomicMass()   == 12 && 
              particle->GetAtomicNumber() == 6))                          return;


   // -----------Kinetic energy at the center of the step ----------------
  const G4double particle_energy =
      step->GetPreStepPoint()->GetKineticEnergy() -
      step->GetTotalEnergyDeposit();

   // ----------Material and process at the and of the step -----------
  const G4Material *material = step->GetPreStepPoint()->GetMaterial();
  const G4VProcess *process  = step->GetPostStepPoint()->GetProcessDefinedStep();

  // -------- Retrieve inelastic process for proton and neutron ----------
  static G4HadronicProcessStore *store = G4HadronicProcessStore::Instance(); 

  static G4VProcess *protonInelastic  =
    store->FindProcess(G4Proton::Proton(),   fHadronInelastic);
  static G4VProcess *neutronInelastic =
    store->FindProcess(G4Neutron::Neutron(), fHadronInelastic);
  static G4VProcess *heliumInelastic  =
    store->FindProcess(G4Alpha::Alpha(),     fHadronInelastic);
  static G4VProcess *carbonInelastic  =
    store->FindProcess(G4IonTable::GetIonTable()->GetIon(6, 12, 0), fHadronInelastic);

   G4VProcess *inelasticProcess =
      (fParticleType == "neutron") ? neutronInelastic :
      (fParticleType == "helium")  ? heliumInelastic  :
      (fParticleType == "carbon")  ? carbonInelastic  :
                                   protonInelastic;


  auto analysisManager = G4AnalysisManager::Instance();

  // ════════════════════════════════════════════════════════════════════════
  // MODE MONO
  // ════════════════════════════════════════════════════════════════════════
  if (!fMultiElement) {

     // ── Fill Ep histogram before inelastic filter ────────────────
     analysisManager->FillH1(fH1EpMono, particle_energy / MeV, 1.0);

     // ── Inelastic filter : keep only inelastic interaction───────────
     if (process != inelasticProcess) return;

     // --- Read the material density once at the first inelastic interaction ---
     if (fDensity == 0.0)
        fDensity = material->GetDensity() / (CLHEP::g / CLHEP::cm3);


   /* Partial inelastic intercation (KE> 0), i.e intercations that NOT stopthe proton, don't produce PG.
     The filter KE> 0  is therefore physically justified. Stats of secondaries are printed at the end of simu*/
   
      if ((fParticleType == "proton"  ||
            fParticleType == "helium"  ||
            fParticleType == "carbon") &&
            step->GetPostStepPoint()->GetKineticEnergy() > 0.) {

      // Count partial inelastic intercations per energy bin
      G4int pBin = static_cast<G4int>(
         (particle_energy/MeV - fParticleMinEnergy/MeV) / (dE/MeV)); // Get the bin index cooresponding to the proton energie
      if (pBin >= 0 && pBin < fParticleNbBins) 
         fNInelasticPartiel[pBin] += 1.0;

         // ---Count secondary particle types for the end of run sammary ---- 
         G4TrackVector *secondaries = (const_cast<G4Step *>(step))->GetfSecondary();
            for (size_t i = 0; i < secondaries->size(); ++i) {
               G4String name = (*secondaries)[i]->GetDefinition()->GetParticleName();
               G4double ekin = (*secondaries)[i]->GetKineticEnergy()/ MeV;
               fSecondairePartielCount[name] += 1;

            if(fSaveKE0Secondaries)
               fEnergiesSecondairesPartiel[name].push_back(ekin);
         }
      return;
   }


      // ---Count total inelastic intrcations that stop the proton (KE=0) per energy bin ---
     G4int pBinInelastic = static_cast<G4int>(
         (particle_energy / MeV - fParticleMinEnergy / MeV) / (dE / MeV));

     if (pBinInelastic >= 0 && pBinInelastic < fParticleNbBins)
        fNInelastic[pBinInelastic] += 1.0;

      // Inelastic cross section at this step -------------
     G4double cross_section =
         store->GetCrossSectionPerVolume(particle, particle_energy,
                                         inelasticProcess, material);

     // ── Fill kappa histogram ─────────────────────────────────
     // kappa_inel is computed at the center of each energy bin
     if (!fSigmaFilled) {
        for (G4int f = 1; f <= fParticleNbBins; ++f) {
           G4double binCenter = energyStart + (f - 1) * dE;
           G4double csLocal   = store->GetCrossSectionPerVolume(
               particle, binCenter, inelasticProcess, material);
           analysisManager->FillH1(fH1KapaMono, binCenter / MeV, csLocal * cm);
        }
        fSigmaFilled = true;
     }

      
     // ── Loop over secondaries to collect PG ──────────────────────────────────
     G4TrackVector *secondaries = (const_cast<G4Step *>(step))->GetfSecondary();

      // ── Investigation secondaires KE = 0 pour hélium et carbone ─────────────
             if (fParticleType == "helium" || fParticleType == "carbon") { 
               for (size_t i = 0; i < secondaries->size(); ++i) { 
                  G4String name = (*secondaries)[i]->GetDefinition()->GetParticleName(); 
                  G4double ekin = (*secondaries)[i]->GetKineticEnergy() / MeV; 
                  fSecondaireTotalCount[name] += 1; 
                  if (fSaveKE0Secondaries) 
                     fEnergiesSecondairesTotal[name].push_back(ekin);
               }
            }

     G4bool producedGamma = false;

     for (size_t i = 0; i < secondaries->size(); ++i) {
        if ((*secondaries)[i]->GetDefinition() != G4Gamma::Gamma()) continue; // Keep only PG

        // Energy filter : keep PG between 0.04 Mev and GammaMaxEnergy
        G4double Egamma = (*secondaries)[i]->GetKineticEnergy(); 
        if (Egamma < 0.04 * MeV || Egamma > fGammaMaxEnergy) continue; 

        G4double Epart_MeV = particle_energy / MeV;
        G4double Eg_MeV = Egamma / MeV;

        // Fill EpEpg and NrPG histogram
        analysisManager->FillH2(fH2EpEpgMono, Epart_MeV, Eg_MeV);
        analysisManager->FillH1(fH1NrPGMono,  Epart_MeV, 1.0);
        
         // Get bin gamma indices for GammaZ accumaluation
        G4int gBin = static_cast<G4int>(
            (Eg_MeV - fGammaMinEnergy / MeV) / (dEgamma / MeV));
        G4int pBin = static_cast<G4int>(
            (Epart_MeV - fParticleMinEnergy / MeV) / (dE / MeV));

         // Tcheck that bin indices are within valid range
        if (gBin >= 0 && gBin < fGammaNbBins &&
            pBin >= 0 && pBin < fParticleNbBins) {
            // accumulate kappa-inel  for GammaZ normalization
           fGammaZData[pBin][gBin] += cross_section;
           fGammaEvents.push_back({Epart_MeV, Eg_MeV, cross_section, pBin});
        }
        producedGamma = true;
     }

      // Track intercations that produce at least one PG
     if (producedGamma)
        fEpInelasticProducedGamma.push_back(particle_energy / MeV);

     return;
  }

  // ════════════════════════════════════════════════════════════════════════
  // MODE MULTI
  // ════════════════════════════════════════════════════════════════════════

  // ── Inelastic filter  ────────────────────────────────────────────────
  if (process != inelasticProcess) return;

  const G4HadronicProcess *hadProc =
      dynamic_cast<const G4HadronicProcess *>(inelasticProcess);
  if (!hadProc) return;

  const G4Nucleus *targetNucleus = hadProc->GetTargetNucleus();
  if (!targetNucleus) return;

  G4int Z = targetNucleus->GetZ_asInt();
  if (fGammaZDataMap.find(Z) == fGammaZDataMap.end()) return;

  analysisManager->FillH1(fH1EpIndex[Z], particle_energy / MeV, 1.0);

 // ── Filtre KE > 0 (cohérent avec mode mono) ─────────────────────────── 
   if (fParticleType == "proton" && step->GetPostStepPoint()->GetKineticEnergy() > 0.) return;


  G4int pBinInelastic = static_cast<G4int>(
      (particle_energy / MeV - fParticleMinEnergy / MeV) / (dE / MeV));
  if (pBinInelastic >= 0 && pBinInelastic < fParticleNbBins)
     fNInelasticMap[Z][pBinInelastic] += 1.0;


    
  G4double cross_section =
      store->GetCrossSectionPerVolume(particle, particle_energy,
                                      inelasticProcess, material);

  if (!fSigmaFilledMap[Z]) {
     for (G4int f = 1; f <= fParticleNbBins; ++f) {
        G4double binCenter = energyStart + (f - 1) * dE;
        G4double csLocal   = store->GetCrossSectionPerVolume(
            particle, binCenter, inelasticProcess, material);
        analysisManager->FillH1(fH1KapaIndex[Z], binCenter / MeV, csLocal * cm);
     }
     fSigmaFilledMap[Z] = true;
  }

  G4TrackVector *secondaries =
      (const_cast<G4Step *>(step))->GetfSecondary();

  G4bool producedGamma = false;

  for (size_t i = 0; i < secondaries->size(); ++i) {
     if ((*secondaries)[i]->GetDefinition() != G4Gamma::Gamma()) continue;

     G4double Egamma = (*secondaries)[i]->GetKineticEnergy();
     if (Egamma < 0.04 * MeV || Egamma > fGammaMaxEnergy) continue;

     G4double Epart_MeV = particle_energy / MeV;
     G4double Eg_MeV = Egamma / MeV;

     G4int gBin = static_cast<G4int>(
         (Eg_MeV - fGammaMinEnergy / MeV) / (dEgamma / MeV));
     G4int pBin = static_cast<G4int>(
         (Epart_MeV - fParticleMinEnergy / MeV) / (dE / MeV));

     if (gBin < 0 || gBin >= fGammaNbBins) continue;
     if (pBin < 0 || pBin >= fParticleNbBins) continue;

     fGammaZDataMap[Z][pBin][gBin] += cross_section;
     fGammaEventsMap[Z].push_back({Epart_MeV, Eg_MeV, cross_section, pBin});

     //analysisManager->FillH2(fH2EpEpgNormalizedIndex[Z], Epart_MeV, Eg_MeV,
      //                        cross_section);
     analysisManager->FillH2(fH2EpEpgIndex[Z], Epart_MeV, Eg_MeV);
     analysisManager->FillH1(fH1NrPGIndex[Z],  Epart_MeV, 1.0);

     producedGamma = true;
  }

  if (producedGamma)
     fEpInelasticProducedGammaMap[Z].push_back(particle_energy / MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GatePromptGammaStatisticActor::EndOfRunAction(const G4Run *run) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int GatePromptGammaStatisticActor::EndOfRunActionMasterThread(int run_id) {

  auto analysisManager = G4AnalysisManager::Instance();

  const G4double dE          = (fParticleMaxEnergy - fParticleMinEnergy) / fParticleNbBins;
  const G4double energyStart = fParticleMinEnergy + 0.5 * dE;

  // ════════════════════════════════════════════════════════════════════════
  // MODE MONO
  // ════════════════════════════════════════════════════════════════════════
  if (!fMultiElement) {

     // Normalisation fGammaZData
     for (G4int i = 0; i < fParticleNbBins; ++i) {
        G4double sf = fNInelastic[i] * fDensity;
        if (sf > 0) {
           for (G4int j = 0; j < fGammaNbBins; ++j)
              fGammaZData[i][j] /= sf;
        }
     }

     // Remplir GammaZ
     for (auto &[Ep, Eg, xs, pBinEvt] : fGammaEvents) {
        G4double sf = fNInelastic[pBinEvt] * fDensity;
        if (sf > 0)
           analysisManager->FillH2(fH2GammaZMono, Ep, Eg, xs / sf);
     }

    
     // EpInelastic depuis fNInelastic
     for (G4int i = 0; i < fParticleNbBins; ++i) {
        if (fNInelastic[i] > 0) {
           G4double binCenter = energyStart + i * dE;
           analysisManager->FillH1(fH1EpInelasticMono,
                                   binCenter / MeV, fNInelastic[i]);
        }
     }

     // EpInelasticProducedGamma
     for (auto &Ep : fEpInelasticProducedGamma)
        analysisManager->FillH1(fH1EpInelasticProducedGammaMono, Ep, 1.0);


    analysisManager->Write();
    analysisManager->CloseFile();

    analysisManager->Reset();

   

      // ── Write KE > 0 secondaries to txt file ─────────────────────────────────
      if (fSaveKE0Secondaries && !fEnergiesSecondairesPartiel.empty()) {
         std::string txtFile = fOutputFilename + "_KE0_secondaries.txt";
         std::ofstream outFile(txtFile);
         if (!fEnergiesSecondairesPartiel.empty()){ 
            for (auto &[name, energies] : fEnergiesSecondairesPartiel)
               for (auto &e : energies)
                  outFile << name << " " << e << "\n";
         }
         outFile.close();
      }

      // ── Write KE = 0 secondaries to txt file ─────────────────────────────────
      if (fSaveKE0Secondaries && !fEnergiesSecondairesTotal.empty()) {
         std::string txtFile = fOutputFilename + "_KE0total_secondaries.txt";
         std::ofstream outFile(txtFile);
         for (auto &[name, energies] : fEnergiesSecondairesTotal)
            for (auto &e : energies)
               outFile << name << " " << e << "\n";
      
         outFile.close();
      }

      return 0;
   }
  

  // ════════════════════════════════════════════════════════════════════════
  // MODE MULTI
  // ════════════════════════════════════════════════════════════════════════

  for (auto &[Z, sym] : fZtoSymbol) {

     for (G4int i = 0; i < fParticleNbBins; ++i) {
        G4double sf = fNInelasticMap[Z][i] * fDensityMap[Z];
        if (sf > 0) {
           for (G4int j = 0; j < fGammaNbBins; ++j)
              fGammaZDataMap[Z][i][j] /= sf;
        }
     }

     for (auto &[Ep, Eg, xs, pBinEvt] : fGammaEventsMap[Z]) {
        G4double sf = fNInelasticMap[Z][pBinEvt] * fDensityMap[Z];
        if (sf > 0)
           analysisManager->FillH2(fH2GammaZIndex[Z], Ep, Eg, xs / sf);
     }


     // EpInelastic par élément
     for (G4int i = 0; i < fParticleNbBins; ++i) {
        if (fNInelasticMap[Z][i] > 0) {
           G4double binCenter = energyStart + i * dE;
           analysisManager->FillH1(fH1EpInelasticIndex[Z],
                                   binCenter / MeV, fNInelasticMap[Z][i]);
        }
     }

     // EpInelasticProducedGamma par élément
     for (auto &Ep : fEpInelasticProducedGammaMap[Z])
        analysisManager->FillH1(fH1EpInelasticProducedGammaIndex[Z], Ep, 1.0);

   
  }

  analysisManager->Write();
  analysisManager->CloseFile();
  return 0;
}
