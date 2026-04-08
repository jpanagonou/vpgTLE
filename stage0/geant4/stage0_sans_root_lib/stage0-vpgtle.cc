// ********************************************************************
// * License and Disclaimer                                           *
// * (see original Geant4 file for full license header)              *
// ********************************************************************
//
/// \file Hadr09-noroot.cc
/// \brief Main — base de données gamma prompts, sans dépendance ROOT directe.
///
/// Reproduit exactement la logique de l'original Hadr09.cc :
///   - Boucle sur humanBodyElements
///   - Un fichier ROOT par élément via G4AnalysisManager
///   - WeightCompute avec les fractions du corps humain (flag "weight")
///   - Fichier standard_Weight.root séparé si "weight" activé
///
/// Options CLI (identiques à l'original) :
///   -n <int>          : nombre de collisions par bin  (défaut : 1000)
///   -p proton|neutron : type de projectile            (défaut : proton)
///   weight            : active WeightCompute / pondération ToF
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Random/Ranlux64Engine.h"
#include "G4GenericIon.hh"
#include "G4HadronicParameters.hh"
#include "G4IonTable.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4ios.hh"
#include "HadronicGenerator-vpgtle.hh"
#include "globals.hh"


#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


int main(int argc, char *argv[])
{
  // ── Désactiver les hyper-noyaux (identique à l'original) ─────────────────
  G4HadronicParameters::Instance()->SetEnableHyperNuclei(false);


  // ── Lecture des arguments CLI ─────────────────────────────────────────────
  bool     stw          = false;
  bool     pro          = false;
  bool     num          = false;
  G4int    numCollisions  = 1000;
  G4String strprojectile  = "proton";


  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];


    if (arg == "weight") {
      stw = true;
    }
    else if (arg == "-n" && i + 1 < argc) {
      numCollisions = std::atoi(argv[++i]);
      if (numCollisions <= 0) {
        G4cout << "WARNING: -n doit etre > 0, valeur 1000 utilisee." << G4endl;
        numCollisions = 1000;
      }
      num = true;
      G4cout << "Nombre de collisions : " << numCollisions << G4endl;
    }
    else if (arg == "-p" && i + 1 < argc) {
      strprojectile = argv[++i];
      if (strprojectile != "proton" && strprojectile != "neutron") {
        G4cout << "WARNING: projectile invalide, proton utilise." << G4endl;
        strprojectile = "proton";
      }
      pro = true;
      G4cout << "Projectile : " << strprojectile << G4endl;
    }
  }


  // ── Initialisation de la physique ─────────────────────────────────────────
  const G4String namePhysics = "QGSP_BIC_HP";
  HadronicGenerator *theHadronicGenerator = new HadronicGenerator(namePhysics);
  if (!theHadronicGenerator) {
    G4cerr << "ERROR: HadronicGenerator est null !" << G4endl;
    return 1;
  }


  // ── Liste des éléments à simuler ──────────────────────────────────────────
  // Décommenter/commenter selon les besoins (identique à l'original).
  const std::vector<G4String> humanBodyElements = {
    //"G4_H",
    //"G4_C",
    //"G4_N",
    "G4_O",
    //"G4_Ca",
    //"G4_P",
    //"G4_K",
    //"G4_S",
    //"G4_Na",
    //"G4_Cl",
    //"G4_Mg",
  };


  // ── Fractions corporelles pour WeightCompute (identique à l'original) ─────
  const std::vector<G4double> Fraction = {
      0.65, 0.185, 0.103, 0.03, 0.015, 0.01, 0.002, 0.002, 0.001, 0.001, 0.001};
  const std::vector<G4String> Body_approx = {
      "G4_O",  "G4_C",  "G4_H",  "G4_N",
      "G4_Ca", "G4_P",  "G4_K",  "G4_S",
      "G4_Na", "G4_Cl", "G4_Mg"};


  // ── Paramètres des histogrammes (identiques à l'original) ─────────────────
  const G4int    protonNbBins    = 500;
  const G4double protonMinEnergy = 0.0;    // valeur numérique en MeV (unités G4)
  const G4double protonMaxEnergy = 200.0;  // MeV
  const G4int    gammaNbBins     = 250;
  const G4double gammaMinEnergy  = 0.0;    // MeV
  const G4double gammaMaxEnergy  = 10.0;   // MeV


  // ── Moteur aléatoire (identique à l'original) ─────────────────────────────
  CLHEP::Ranlux64Engine defaultEngine(1234567, 4);
  CLHEP::HepRandom::setTheEngine(&defaultEngine);
  CLHEP::HepRandom::setTheSeed(static_cast<long>(time(nullptr)));


  // ── Table de particules ───────────────────────────────────────────────────
  G4GenericIon *gion = G4GenericIon::GenericIon();
  gion->SetProcessManager(new G4ProcessManager(gion));
  G4ParticleTable *partTable = G4ParticleTable::GetParticleTable();
  G4IonTable     *ions       = partTable->GetIonTable();
  partTable->SetReadiness();
  ions->CreateAllIon();
  ions->CreateAllIsomer();


  // ── Nombre de collisions effectif ─────────────────────────────────────────
  G4int nbCollisions = num ? numCollisions : 1000;
  G4cout << "Nombre de collisions par bin : " << nbCollisions << G4endl;


  // ── Projectile ────────────────────────────────────────────────────────────
  G4ParticleDefinition *projectile =
      partTable->FindParticle(pro ? strprojectile : G4String("proton"));
  G4cout << "Projectile : " << projectile->GetParticleName() << G4endl;


  // ── Accumulateur global Weight (tous éléments, si stw) ───────────────────
  // Reproduit TH1D_weight de l'original accumulé sur tous les éléments.
  std::vector<G4double> globalWeight(protonNbBins, 0.0);


  // ════════════════════════════════════════════════════════════════════════════
  //  BOUCLE SUR LES ÉLÉMENTS CHIMIQUES
  //  Reproduit : for (G4int k = 0; k < humanbodyindex; k++) de l'original
  // ════════════════════════════════════════════════════════════════════════════
  for (const G4String &nameMaterial : humanBodyElements) {


    // ── Matériau ──────────────────────────────────────────────────────────
    G4Material *targetMaterial =
        G4NistManager::Instance()->FindOrBuildMaterial(nameMaterial);
    if (!targetMaterial) {
      G4cerr << "ERROR: materiau " << nameMaterial << " introuvable !" << G4endl;
      return 3;
    }
    const G4Element *targetElement = targetMaterial->GetElement(0);
    G4String symbol = targetElement->GetSymbol();


    G4cout << "\n=== Element : " << symbol
           << " (" << nameMaterial << ") ===" << G4endl;


    // ── Ouverture du fichier ROOT et réservation des histogrammes ─────────
    // Un fichier par élément : PG_O.root, PG_C.root, etc.
    // Reproduit le TDirectory par élément de l'original.
    // InitAnalysis crée via G4AnalysisManager :
    //   H2[0] EpEpg, H2[1] GammaZ, H1[0] NrPG, H1[1] Sigma, H1[2] Weight
    G4String fileName = "PG_" + symbol;
    theHadronicGenerator->InitAnalysis(
        fileName,
        protonNbBins, protonMinEnergy, protonMaxEnergy,
        gammaNbBins,  gammaMinEnergy,  gammaMaxEnergy);


    // ── Shadow array gammaZData ───────────────────────────────────────────
    // gammaZData[i][j] accumule les poids pour le bin proton i et gamma j.
    // Rempli en parallèle de H2 GammaZ dans GenerateInteraction().
    // Remplace TH2D::Integral() utilisé dans WeightCompute() de l'original.
    std::vector<std::vector<G4double>> gammaZData(
        protonNbBins, std::vector<G4double>(gammaNbBins, 0.0));


    // ── Simulation ────────────────────────────────────────────────────────
    G4int nbPG = theHadronicGenerator->GenerateInteraction(
        projectile,
        targetMaterial,
        nbCollisions,
        protonNbBins, protonMinEnergy, protonMaxEnergy,
        gammaNbBins,  gammaMinEnergy,  gammaMaxEnergy,
        gammaZData);


    G4cout << "Nombre de PG pour " << nameMaterial << " : " << nbPG << G4endl;


    // ── WeightCompute (si flag "weight" activé) ───────────────────────────
    G4String done = "not counting in the approximation";
    if (stw) {
      auto ind = std::find(Body_approx.begin(), Body_approx.end(), nameMaterial);
      if (ind != Body_approx.end()) {
        size_t index = std::distance(Body_approx.begin(), ind);


        // Remplit H1[2] Weight dans le fichier courant
        done = theHadronicGenerator->WeightCompute(
            gammaZData,
            protonNbBins, protonMinEnergy, protonMaxEnergy,
            Fraction[index]);


        // Accumulation dans le vecteur global pour standard_Weight.root
        const G4double dE =
            (protonMaxEnergy - protonMinEnergy) / protonNbBins;
        for (G4int b = 0; b < protonNbBins; ++b) {
          G4double colSum = 0.0;
          for (G4double v : gammaZData[b]) colSum += v;
          globalWeight[b] += Fraction[index] * colSum;
        }
      }
    }
    G4cout << nameMaterial << " : " << done << G4endl;


    // ── Fermeture du fichier ROOT de cet élément ──────────────────────────
    // OBLIGATOIRE avant InitAnalysis() de l'élément suivant :
    // G4AnalysisManager est un singleton, il faut fermer avant de rouvrir.
    theHadronicGenerator->EndAnalysis();


  } // fin boucle éléments


  // ── Fichier Weight global (si flag "weight") ──────────────────────────────
  // Reproduit : file->mkdir("standard_Weight") + TH1D_weight->Write()
  if (stw) {
    auto am = G4AnalysisManager::Instance();
    am->SetDefaultFileType("root");
    am->SetVerboseLevel(0);
    am->SetFileName("../output/standard_Weight");


    am->CreateH1("Weight",
                 "Weighted PG yield for ToF (body approximation);"
                 "Proton energy [MeV];Weight [mm^{-1}]",
                 protonNbBins, protonMinEnergy, protonMaxEnergy);


    am->OpenFile();


    const G4double dE =
        (protonMaxEnergy - protonMinEnergy) / protonNbBins;
    for (G4int b = 0; b < protonNbBins; ++b) {
      G4double binCenter = protonMinEnergy + (b + 0.5) * dE;
      am->FillH1(0, binCenter / MeV, globalWeight[b]);
    }


    am->Write();
    am->CloseFile();


    G4cout << "\nWeight global ecrit dans ../output/standard_Weight.root"
           << G4endl;
  }


  // ── Nettoyage ─────────────────────────────────────────────────────────────
  delete theHadronicGenerator;
  return 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......







