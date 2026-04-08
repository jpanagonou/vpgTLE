//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes, nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file HadronicGenerator.cc
/// \brief Implementation of the HadronicGenerator class
//
//------------------------------------------------------------------------
// Class: HadronicGenerator
// Author: Alberto Ribon (CERN EP/SFT)
// Date: May 2020
//------------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HadronicGenerator-vpgtle.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>

#include "G4AblaInterface.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4Box.hh"
#include "G4CascadeInterface.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4DecayPhysics.hh"
#include "G4DynamicParticle.hh"
#include "G4ExcitationHandler.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4FTFModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronicParameters.hh"
#include "G4HadronicProcessStore.hh"
#include "G4INCLXXInterface.hh"
#include "G4IonTable.hh"
#include "G4LundStringFragmentation.hh"
#include "G4Material.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4PreCompoundModel.hh"
#include "G4ProcessManager.hh"
#include "G4QGSMFragmentation.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QuasiElasticChannel.hh"
#include "G4StateManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4TheoFSGenerator.hh"
#include "G4TouchableHistory.hh"
#include "G4TransportationManager.hh"
#include "G4UnitsTable.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4VParticleChange.hh"
#include "G4ios.hh"
#include "globals.hh"

// particle libraries
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4He3.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Triton.hh"

// HP libraries
#include "G4NeutronCaptureProcess.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticXS.hh"
#include "G4ParticleHPInelastic.hh"

// Remplacement de ROOT
#include "G4AnalysisManager.hh" // 

#include "G4Element.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronicGenerator::HadronicGenerator(const G4String physicsCase)
    : fPhysicsCase(physicsCase), fLastHadronicProcess(nullptr),
      fPartTable(nullptr) {

  // The constructor set-ups all the particles, models, cross sections and
  // hadronic inelastic processes.
  // This should be done only once for each application.
  // In the case of a multi-threaded application using this class,
  // the constructor should be invoked for each thread,
  // i.e. one instance of the class should be kept per thread.
  // The particles and processes that are created in this constructor
  // will then be used by the method GenerateInteraction at each interaction.
  // Notes:
  // - Neither the hadronic models nor the cross sections are used directly
  //   by the method GenerateInteraction, but they are associated to the
  //   hadronic processes and used by Geant4 to simulate the collision;
  // - Although the class generates only final states, but not free mean paths,
  //   inelastic hadron-nuclear cross sections are needed by Geant4 to sample
  //   the target nucleus from the target material.

  if (fPhysicsCase != "QGSP_BIC_HP") {
    G4cerr
        << "ERROR: Not supported final-state hadronic inelastic physics case !"
        << fPhysicsCase << G4endl
        << "\t Re-try by choosing one of the following:" << G4endl
        << "\t - Hadronic models : QGSP" << G4endl
        << "\t - \"Physics-list proxies\" : QGSP_BIC_HP" << G4endl;
  }

  // Definition of particles
  G4GenericIon *gion = G4GenericIon::Definition();
  gion->SetProcessManager(new G4ProcessManager(gion));
  G4DecayPhysics *decays = new G4DecayPhysics;
  decays->ConstructParticle();
  fPartTable = G4ParticleTable::GetParticleTable();
  fPartTable->SetReadiness();
  G4IonTable *ions = fPartTable->GetIonTable();
  ions->CreateAllIon();
  ions->CreateAllIsomer();

  // Building the ElementTable for the NeutronHP model
  // TO BE MODIFIED :: there should be a method in Geant4 for this purpose
  std::vector<G4Element *> myElements;
  myElements.push_back(new G4Element("Hydrogen", "H", 1., 1.01 * g / mole));
  myElements.push_back(new G4Element("Helium", "He", 2., 4.00 * g / mole));
  myElements.push_back(new G4Element("Lithium", "Li", 3., 6.94 * g / mole));
  myElements.push_back(new G4Element("Beryllium", "Be", 4., 9.01 * g / mole));
  myElements.push_back(new G4Element("Boron", "B", 5., 10.81 * g / mole));
  myElements.push_back(new G4Element("Carbon", "C", 6., 12.01 * g / mole));
  myElements.push_back(new G4Element("Nitrogen", "N", 7., 14.01 * g / mole));
  myElements.push_back(new G4Element("Oxygen", "O", 8., 16.00 * g / mole));
  myElements.push_back(new G4Element("Fluorine", "F", 9., 19.00 * g / mole));
  myElements.push_back(new G4Element("Neon", "Ne", 10., 20.18 * g / mole));
  myElements.push_back(new G4Element("Sodium", "Na", 11., 22.99 * g / mole));
  myElements.push_back(new G4Element("Magnesium", "Mg", 12., 24.31 * g / mole));
  myElements.push_back(new G4Element("Aluminium", "Al", 13., 26.98 * g / mole));
  myElements.push_back(new G4Element("Silicon", "Si", 14., 28.09 * g / mole));
  myElements.push_back(new G4Element("Phosphorus", "P", 15., 30.97 * g / mole));
  myElements.push_back(new G4Element("Sulfur", "S", 16., 32.07 * g / mole));
  myElements.push_back(new G4Element("Chlorine", "Cl", 17., 35.45 * g / mole));
  myElements.push_back(new G4Element("Argon", "Ar", 18., 39.95 * g / mole));
  myElements.push_back(new G4Element("Potassium", "K", 19., 39.10 * g / mole));
  myElements.push_back(new G4Element("Calcium", "Ca", 20., 40.08 * g / mole));
  myElements.push_back(new G4Element("Titanium", "Ti", 22., 47.87 * g / mole));
  myElements.push_back(new G4Element("Copper", "Cu", 29., 63.55 * g / mole));
  myElements.push_back(new G4Element("Zinc", "Zn", 30., 65.38 * g / mole));
  myElements.push_back(new G4Element("Silver", "Ag", 47., 107.87 * g / mole));
  myElements.push_back(new G4Element("Tin", "Sn", 50., 118.71 * g / mole));

  // Build BIC model
  G4BinaryCascade *theBICmodel = new G4BinaryCascade;
  G4PreCompoundModel *thePreEquilib =
      new G4PreCompoundModel(new G4ExcitationHandler);
  theBICmodel->SetDeExcitation(thePreEquilib);

  // Build BinaryLightIon model
  G4PreCompoundModel *thePreEquilibBis =
      new G4PreCompoundModel(new G4ExcitationHandler);
  G4BinaryLightIonReaction *theIonBICmodel =
      new G4BinaryLightIonReaction(thePreEquilibBis);

  // HP model
  G4NeutronHPInelastic *theNeutronHP = new G4NeutronHPInelastic;
  theNeutronHP->BuildPhysicsTable(*G4Neutron::Neutron());
  //theNeutronHP->SetHadrGenFlag(
  //    true); // Enable hadronic generator flag :: due to the ElementTable
             // building default :: TO BE MODIFIED

  // Model instance with constraint to be above a kinetic energy threshold.
  // (Used for ions in all physics lists)
  G4GeneratorPrecompoundInterface *theCascade =
      new G4GeneratorPrecompoundInterface;
  theCascade->SetDeExcitation(thePreEquilib);
  G4LundStringFragmentation *theLundFragmentation =
      new G4LundStringFragmentation;
  G4ExcitedStringDecay *theStringDecay =
      new G4ExcitedStringDecay(theLundFragmentation);
  G4FTFModel *theStringModel = new G4FTFModel;
  theStringModel->SetFragmentationModel(theStringDecay);

  G4TheoFSGenerator *theFTFPmodel_aboveThreshold =
      new G4TheoFSGenerator("FTFP");
  theFTFPmodel_aboveThreshold->SetMaxEnergy(
      G4HadronicParameters::Instance()->GetMaxEnergy());
  theFTFPmodel_aboveThreshold->SetTransport(theCascade);
  theFTFPmodel_aboveThreshold->SetHighEnergyGenerator(theStringModel);

  // Model instance with constraint to be within two kinetic energy thresholds.
  // (Used in the case of QGS-based physics lists for nucleons)
  G4TheoFSGenerator *theFTFPmodel_constrained = new G4TheoFSGenerator("FTFP");
  theFTFPmodel_constrained->SetMaxEnergy(
      G4HadronicParameters::Instance()->GetMaxEnergy());
  theFTFPmodel_constrained->SetTransport(theCascade);
  theFTFPmodel_constrained->SetHighEnergyGenerator(theStringModel);

  // Build the QGSP model
  G4TheoFSGenerator *theQGSPmodel = new G4TheoFSGenerator("QGSP");
  theQGSPmodel->SetMaxEnergy(G4HadronicParameters::Instance()->GetMaxEnergy());
  theQGSPmodel->SetTransport(theCascade);
  G4QGSMFragmentation *theQgsmFragmentation = new G4QGSMFragmentation;
  G4ExcitedStringDecay *theQgsmStringDecay =
      new G4ExcitedStringDecay(theQgsmFragmentation);
  G4VPartonStringModel *theQgsmStringModel = new G4QGSModel<G4QGSParticipants>;
  theQgsmStringModel->SetFragmentationModel(theQgsmStringDecay);
  theQGSPmodel->SetHighEnergyGenerator(theQgsmStringModel);
  G4QuasiElasticChannel *theQuasiElastic =
      new G4QuasiElasticChannel; // QGSP uses quasi-elastic
  theQGSPmodel->SetQuasiElasticChannel(theQuasiElastic);

  // For the case of "physics-list proxies", select the energy range for each
  // hadronic model.
  const G4double ftfpMinE =
      G4HadronicParameters::Instance()->GetMinEnergyTransitionFTF_Cascade();
  const G4double BICMaxE =
      G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade();
  const G4double ftfpMaxE =
      G4HadronicParameters::Instance()->GetMaxEnergyTransitionQGS_FTF();
  const G4double qgspMinE =
      G4HadronicParameters::Instance()->GetMinEnergyTransitionQGS_FTF();

  theBICmodel->SetMaxEnergy(BICMaxE);
  theIonBICmodel->SetMaxEnergy(BICMaxE);
  theFTFPmodel_aboveThreshold->SetMinEnergy(ftfpMinE);
  theFTFPmodel_constrained->SetMinEnergy(ftfpMinE);
  theFTFPmodel_constrained->SetMaxEnergy(ftfpMaxE);
  theQGSPmodel->SetMinEnergy(qgspMinE);

  // Cross sections (needed by Geant4 to sample the target nucleus from the
  // target material)
  G4VCrossSectionDataSet *theProtonXSdata =
      new G4BGGNucleonInelasticXS(G4Proton::Proton());
  theProtonXSdata->BuildPhysicsTable(*(G4Proton::Definition()));

  G4VCrossSectionDataSet *theNeutronXSdata = new G4NeutronInelasticXS;
  theNeutronXSdata->BuildPhysicsTable(*(G4Neutron::Definition()));

  G4VCrossSectionDataSet *theNeutronHPXSdatainel = new G4NeutronHPInelasticXS;
  theNeutronHPXSdatainel->BuildPhysicsTable(*G4Neutron::Definition());

  G4VCrossSectionDataSet *theHyperonsXSdata =
      new G4CrossSectionInelastic(new G4ComponentGGHadronNucleusXsc);

  G4VCrossSectionDataSet *theNuclNuclXSdata =
      new G4CrossSectionInelastic(new G4ComponentGGNuclNuclXsc);

  // Set up inelastic processes : store them in a map (with particle definition
  // as key)
  //                              for convenience

  typedef std::pair<G4ParticleDefinition *, G4HadronicProcess *> ProcessPair;

  G4HadronicProcess *theProtonInelasticProcess =
      new G4HadronInelasticProcess("protonInelastic", G4Proton::Definition());
  fProcessMap.insert(
      ProcessPair(G4Proton::Definition(), theProtonInelasticProcess));

  // a specific processmap is dedidcated to the HP physics for neutron
  G4HadronicProcess *theNeutronHPInelasticProcess =
      new G4HadronInelasticProcess("NeutronHPInelastic",
                                   G4Neutron::Definition());
  fProcessMapHP.insert(
      ProcessPair(G4Neutron::Definition(), theNeutronHPInelasticProcess));

  G4HadronicProcess *theNeutronInelasticProcess =
      new G4HadronInelasticProcess("neutronInelastic", G4Neutron::Definition());
  fProcessMap.insert(
      ProcessPair(G4Neutron::Definition(), theNeutronInelasticProcess));

  // For the HP model, we need to create a new process for the neutron
  // Prompt-gamma timing with a track-length estimator in proton therapy, JM
  // Létang et al, 2024

  G4HadronicProcess *theDeuteronInelasticProcess =
      new G4HadronInelasticProcess("dInelastic", G4Deuteron::Definition());
  fProcessMap.insert(
      ProcessPair(G4Deuteron::Definition(), theDeuteronInelasticProcess));

  G4HadronicProcess *theTritonInelasticProcess =
      new G4HadronInelasticProcess("tInelastic", G4Triton::Definition());
  fProcessMap.insert(
      ProcessPair(G4Triton::Definition(), theTritonInelasticProcess));

  G4HadronicProcess *theHe3InelasticProcess =
      new G4HadronInelasticProcess("he3Inelastic", G4He3::Definition());
  fProcessMap.insert(ProcessPair(G4He3::Definition(), theHe3InelasticProcess));

  G4HadronicProcess *theAlphaInelasticProcess =
      new G4HadronInelasticProcess("alphaInelastic", G4Alpha::Definition());
  fProcessMap.insert(
      ProcessPair(G4Alpha::Definition(), theAlphaInelasticProcess));

  G4HadronicProcess *theIonInelasticProcess =
      new G4HadronInelasticProcess("ionInelastic", G4GenericIon::Definition());
  fProcessMap.insert(
      ProcessPair(G4GenericIon::Definition(), theIonInelasticProcess));

  // Add the cross sections to the corresponding hadronic processes

  // cross-sections for nucleons
  theProtonInelasticProcess->AddDataSet(theProtonXSdata);
  theNeutronInelasticProcess->AddDataSet(theNeutronXSdata);

  // specific cross-sections for neutrons in HP model
  theNeutronHPInelasticProcess->AddDataSet(theNeutronHPXSdatainel);

  // cross-sections for ions
  theDeuteronInelasticProcess->AddDataSet(theNuclNuclXSdata);
  theTritonInelasticProcess->AddDataSet(theNuclNuclXSdata);
  theHe3InelasticProcess->AddDataSet(theNuclNuclXSdata);
  theAlphaInelasticProcess->AddDataSet(theNuclNuclXSdata);
  theIonInelasticProcess->AddDataSet(theNuclNuclXSdata);

  // Register the proper hadronic model(s) to the corresponding hadronic
  // processes. in the physics list QGSP_BIC, BIC is used only for nucleons

  // processes for BIC => mean-energy on nucleons
  theProtonInelasticProcess->RegisterMe(theBICmodel);
  theNeutronInelasticProcess->RegisterMe(theBICmodel);

  // processes for QGSP => High-energy on nucleons
  theProtonInelasticProcess->RegisterMe(theQGSPmodel);
  theNeutronInelasticProcess->RegisterMe(theQGSPmodel);

  // processes for HP => low-energy on neutrons
  theNeutronHPInelasticProcess->RegisterMe(theNeutronHP);

  // processes for ion (BIC for ions)
  theDeuteronInelasticProcess->RegisterMe(theIonBICmodel);
  theTritonInelasticProcess->RegisterMe(theIonBICmodel);
  theHe3InelasticProcess->RegisterMe(theIonBICmodel);
  theAlphaInelasticProcess->RegisterMe(theIonBICmodel);
  theIonInelasticProcess->RegisterMe(theIonBICmodel);

  G4TheoFSGenerator *theFTFPmodelToBeUsed = theFTFPmodel_aboveThreshold;

  theFTFPmodelToBeUsed = theFTFPmodel_constrained;
  theNeutronInelasticProcess->RegisterMe(theFTFPmodelToBeUsed);
  theProtonInelasticProcess->RegisterMe(theFTFPmodelToBeUsed);

  theFTFPmodelToBeUsed = theFTFPmodel_aboveThreshold;
  theDeuteronInelasticProcess->RegisterMe(theFTFPmodelToBeUsed);
  theTritonInelasticProcess->RegisterMe(theFTFPmodelToBeUsed);
  theHe3InelasticProcess->RegisterMe(theFTFPmodelToBeUsed);
  theAlphaInelasticProcess->RegisterMe(theFTFPmodelToBeUsed);
  theIonInelasticProcess->RegisterMe(theFTFPmodelToBeUsed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronicGenerator::~HadronicGenerator() { fPartTable->DeleteAllParticles(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



// =====================================================
//  INIT ANALYSIS (paramètres depuis main)
// =====================================================

void HadronicGenerator::InitAnalysis(
    const G4String &fileName,
    G4int       protonNbBins,
    G4double    protonMinEnergy,
    G4double    protonMaxEnergy,
    G4int       gammaNbBins,
    G4double    gammaMinEnergy,
    G4double    gammaMaxEnergy)

{
    auto analysisManager = G4AnalysisManager::Instance();
    

    analysisManager->SetDefaultFileType("root");
    analysisManager->SetFileName("../output/" + fileName);
    analysisManager->SetVerboseLevel(1);

    // H2: id =0, Ep vs Epg = spectre 2D non pondéré, equivalent à TH2D_EpEpg
    analysisManager->CreateH2("EpEpg", 
                            "PGs energy vs proton energy;"
                            "Proton energy [MeV];PG energy [MeV];counts",
                            protonNbBins, protonMinEnergy, protonMaxEnergy,
                            gammaNbBins, gammaMinEnergy, gammaMaxEnergy);

    // H2 : id =1 GammaZ -spectre pondéré par kappa/rho/collision
    //                      equivalent à TH2D_GammaZ
    analysisManager->CreateH2("GammaZ",
                             "PG yield weighed by kappa/rho;"
                             "Proton energy [MeV];PG energy [MeV];kappa [mm^{-1} per collision]",
                            protonNbBins, protonMinEnergy, protonMaxEnergy,
                            gammaNbBins, gammaMinEnergy, gammaMaxEnergy);

    // H1 : id =0 : NrPG - nombre de gammas prompts par bin d'energie proton
                           // equivalent à TH1_NrPG
    analysisManager->CreateH1("NrPG", 
                            "Number of prompts gammas per bin;"
                            "Proton energy [MeV];N_{PG}",
                             protonNbBins, protonMinEnergy, protonMaxEnergy);

    // H1 : id =1 : Sigma - coefficient d'attenuation lineaire inelastique
    //                      equivalent de TH1D_SigmaInelastic
    analysisManager->CreateH1("Sigma",
                             "Linear attenuation coefficient;"
                             "Proton energy [MeV];#kappa_{inel} [cm^{-1}]",
                            protonNbBins, protonMinEnergy, protonMaxEnergy);

    // H1 id = 2 : weght - integral pondéré par le calcul ToF
    //                      Rempli par weghtCompute() apres GenerateIntercation()
    analysisManager->CreateH1("Weight", 
                            "Weight PG yield for ToF;"
                            "Proton energy  [MeV];Weight [mm^{-1}]",
                            protonNbBins, protonMinEnergy, protonMaxEnergy);
  
    analysisManager->OpenFile();
}

// =====================================================

void HadronicGenerator::EndAnalysis()
{
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}

// =====================================================
//  GENERATE INTERACTION
// =====================================================

G4int HadronicGenerator::GenerateInteraction(
    G4ParticleDefinition*    projectileDefinition,
    G4Material*              targetMaterial,
    G4int                    nbCollisions,
    G4int                    protonNbBins,
    G4double                 protonMinEnergy,
    G4double                 protonMaxEnergy,
    G4int                    gammaNbBins,
    G4double                 gammaMinEnergy,
    G4double                 gammaMaxEnergy,
    std::vector<std::vector<G4double>> &gammaZData)

{

        if (projectileDefinition == nullptr){
            G4cerr << "ERROR: projectileDefinition est nullptr !" << G4endl;
            return 0;
        }

    auto analysisManager = G4AnalysisManager::Instance();

    G4VParticleChange* aChange = nullptr;
    G4int nbPG = 0;

    const G4Element* targetElement = targetMaterial->GetElement(0);
    G4double molarMass = targetElement->GetA() / (g / mole);


   /*-------------------------------
     GEOMETRIE MINIMALE OBLIGATOIRE
     -------------------------------
   */
    const G4double dimX = 1.0 * mm;
    const G4double dimY = 1.0 * mm;
    const G4double dimZ = 1.0 * mm;

    G4Box *sFrame = new G4Box("Box", dimX, dimY, dimZ);
    G4LogicalVolume *lFrame =
      new G4LogicalVolume(sFrame, targetMaterial, "Box", 0, 0, 0);

    G4PVPlacement *pFrame =
      new G4PVPlacement(0, G4ThreeVector(), "Box", lFrame, 0, false, 0);

    G4TransportationManager::GetTransportationManager()->SetWorldForTracking(
      pFrame);

    G4ThreeVector projectileDirection = G4ThreeVector(0.0, 0.0, 1.0);;

    
   // ── Géométrie des bins ────────────────────────────────────────────────────
  // Reproduit exactement le comportement de TH1D avec protonNbBins bins
  // uniformes entre protonMinEnergy et protonMaxEnergy.
  //
  //   dE          = largeur d'un bin
  //   energyStart = centre du premier bin  (= BinCenter(1) dans ROOT)
  //   energyEnd   = centre du dernier bin  (= BinCenter(N) dans ROOT)
  //
  const G4double dE =
      (protonMaxEnergy - protonMinEnergy) / static_cast<G4double>(protonNbBins);
  const G4double energyStart     = protonMinEnergy + 0.5 * dE;
  const G4double energyIncrement = dE;


  // ── Largeur des bins gamma (pour l'index dans gammaZData) ─────────────────
  const G4double dEgamma =
      (gammaMaxEnergy - gammaMinEnergy) / static_cast<G4double>(gammaNbBins);


  // ── Vecteurs de stockage pour NrPG et Sigma ───────────────────────────────
  // Ces vecteurs remplacent SetBinContent() de ROOT.
  // Ils sont remplis bin par bin dans la boucle,
  // puis transférés dans les histos G4AnalysisManager en une seule passe.
  std::vector<G4double> vecNrPG(protonNbBins, 0.0);
  std::vector<G4double> vecSigma(protonNbBins, 0.0);


  // Initialiser gammaZData à zéro (dimension : protonNbBins x gammaNbBins)
  for (auto &row : gammaZData)
    std::fill(row.begin(), row.end(), 0.0);


  //────────────────────────────────────────────────────────────────────────────
  // BOUCLE SUR LES BINS D'ÉNERGIE PROTON
  // f = 1..protonNbBins  (convention ROOT 1-based reproduite ici)
  //────────────────────────────────────────────────────────────────────────────
  for (G4int f = 1; f <= protonNbBins; ++f) {


    // Centre physique du bin (utilisé pour l'échantillonnage et Sigma)
    // = energyStart + (f-1)*dE  ≡  TH1D::GetBinCenter(f) dans l'original
    const G4double binCenter = energyStart + (f - 1) * dE;


    G4int Ngamma = 0;


    // ── Création de la particule projectile ──────────────────────────────────
    // On crée une G4DynamicParticle et une G4Track réutilisés pour toutes
    // les collisions de ce bin. L'énergie sera mise à jour à chaque collision.
    auto *dParticle = new G4DynamicParticle(
        projectileDefinition, projectileDirection, binCenter * CLHEP::MeV);


    gTrack = new G4Track(dParticle, 0.0, G4ThreeVector(0, 0, 0));


    auto *step   = new G4Step;
    step->SetTrack(gTrack);
    gTrack->SetStep(step);


    auto *aPoint = new G4StepPoint;
    aPoint->SetPosition(G4ThreeVector(0, 0, 0));
    aPoint->SetMaterial(targetMaterial);
    step->SetPreStepPoint(aPoint);
    gTrack->SetStep(step);


    //──────────────────────────────────────────────────────────────────────────
    // BOUCLE SUR LES COLLISIONS (nbCollisions par bin)
    //──────────────────────────────────────────────────────────────────────────
    for (G4int i = 0; i < nbCollisions; ++i) {


      // Énergie aléatoire dans le bin — même formule que l'original
      // (TH1D::GetBinCenter(f) + (G4UniformRand()-0.5)*GetBinWidth(1))
      G4double projectileEnergy =
          (binCenter + (G4UniformRand() - 0.5) * energyIncrement) * CLHEP::MeV;


      dParticle->SetKineticEnergy(projectileEnergy);
      gTrack->SetKineticEnergy(projectileEnergy);


      // Sélection du processus hadronique pour ce projectile
      G4ParticleDefinition *theProjectileDef =
          projectileDefinition->IsGeneralIon()
              ? G4GenericIon::Definition()
              : projectileDefinition;


      auto mapIndex = fProcessMap.find(theProjectileDef);
      G4HadronicProcess *theProcess =
          (mapIndex != fProcessMap.end()) ? mapIndex->second : nullptr;


      if (theProcess == nullptr) {
        G4cerr << "ERROR: pas de processus pour "
               << theProjectileDef->GetParticleName() << G4endl;
        continue;
      }


      fLastHadronicProcess = theProcess;
      aChange = theProcess->PostStepDoIt(*gTrack, *step);
      if (!aChange) continue;


      // ── Section efficace et poids ──────────────────────────────────────────
      // kappaPerRho : coefficient d'atténuation linéaire [cm-1]
      // weightEvent : poids par événement pour GammaZ [mm-1 / collision]
      //               reproduit kappaPerRhoMmCollisions de l'original
      G4double CrossSection =
          GetCrossSection(dParticle, targetElement, targetMaterial);
      G4double kappaPerRho  = (CrossSection * Avogadro) / molarMass;
      G4double weightEvent  = kappaPerRho / (nbCollisions * cm);


      // ── Boucle sur les particules secondaires ─────────────────────────────
      G4int nsec = aChange->GetNumberOfSecondaries();
      for (G4int j = 0; j < nsec; ++j) {
        const G4DynamicParticle *sec =
            aChange->GetSecondary(j)->GetDynamicParticle();


        if (sec->GetDefinition()->GetParticleName() == "gamma") {
          G4double Egamma = sec->GetKineticEnergy();


          // Seuil bas : on ignore les gammas de très basse énergie
          if (Egamma > 0.04 * MeV) {
            Ngamma++;
            nbPG++;


            G4double Ep_MeV = projectileEnergy / MeV;
            G4double Eg_MeV = Egamma / MeV;


            // H2 id=0 : EpEpg — remplissage non pondéré
            analysisManager->FillH2(0, Ep_MeV, Eg_MeV);


            // H2 id=1 : GammaZ — remplissage pondéré
            analysisManager->FillH2(1, Ep_MeV, Eg_MeV, weightEvent);


            // Shadow array gammaZData : même pondération que GammaZ
            // Utilisé ensuite par WeightCompute() sans dépendance ROOT.
            // On détermine l'index du bin gamma (0-based).
            if (Eg_MeV >= gammaMinEnergy && Eg_MeV < gammaMaxEnergy) {
              G4int gBin = static_cast<G4int>(
                  (Eg_MeV - gammaMinEnergy) / dEgamma);
              if (gBin >= 0 && gBin < gammaNbBins) {
                gammaZData[f - 1][gBin] += weightEvent;
              }
            }
          }
        }
        delete aChange->GetSecondary(j);
      }
      aChange->Clear();


    } // fin boucle collisions


    // ── Calcul de Sigma au centre exact du bin ────────────────────────────────
    // Identique à l'original : on évalue la section efficace à binCenter
    // (sans aléatoire) puis on stocke dans le vecteur vecSigma.
    // Le vecteur sera transféré dans l'histo APRÈS la boucle sur f.
    dParticle->SetKineticEnergy(binCenter * CLHEP::MeV);
    G4double CrossSectionBin =
        GetCrossSection(dParticle, targetElement, targetMaterial);
    G4double kappaPerRhoBin = (CrossSectionBin * Avogadro) / molarMass;


    // Stockage dans les vecteurs (PAS de FillH1 ici)
    vecNrPG[f - 1]  = static_cast<G4double>(Ngamma);
    vecSigma[f - 1] = kappaPerRhoBin;


    delete step;
    delete dParticle;


  } // fin boucle bins


  //──────────────────────────────────────────────────────────────────────────
  // REMPLISSAGE DES HISTOS 1D (NrPG et Sigma) EN UNE SEULE PASSE
  //
  // On parcourt les vecteurs une fois et on appelle FillH1() une seule fois
  // par bin. Cela reproduit exactement SetBinContent() de ROOT :
  // chaque bin reçoit une unique valeur, sans accumulation parasite.
  //──────────────────────────────────────────────────────────────────────────
  for (G4int f = 1; f <= protonNbBins; ++f) {
    G4double binCenter = energyStart + (f - 1) * dE;


    // H1 id=0 : NrPG — nombre de gammas prompts dans ce bin
    analysisManager->FillH1(0, binCenter , vecNrPG[f - 1]);


    // H1 id=1 : Sigma — coefficient d'atténuation linéaire pour ce bin
    analysisManager->FillH1(1, binCenter , vecSigma[f - 1]);
  }
  // Note : H1 id=2 (Weight) est rempli par WeightCompute(), pas ici.


  // ── Nettoyage géométrie ───────────────────────────────────────────────────
  delete pFrame;
  delete lFrame;
  delete sFrame;


  return nbPG;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  GetCrossSection
//  Inchangé par rapport à l'original.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double HadronicGenerator::GetCrossSection(
    const G4DynamicParticle *part,
    const G4Element         *elm,
    const G4Material        *targetMaterial)
{
  G4HadronicProcess *hadProcess = GetHadronicProcess();
  G4double xs =
      hadProcess->GetCrossSectionDataStore()
                ->GetCrossSection(part, elm, targetMaterial);
  return xs / (cm * cm);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  WeightCompute
//
//  Equivalent de l'original WeightCompute(TH2D*, TH1D*, frac) sans ROOT.
//
//  L'original faisait :
//    TH2D::Integral(i, i, 1, NbinsY)  → somme de la colonne i de GammaZ
//    puis TH1D::SetBinContent(i, ancien + frac * colonne)
//
//  Ici on opère sur gammaZData (shadow array rempli par GenerateInteraction)
//  et on remplit H1 id=2 (Weight) via FillH1().
//
//  ATTENTION : comme H1 id=2 peut être rempli plusieurs fois (une fois par
//  élément chimique dans la boucle du main), FillH1 accumule correctement
//  les contributions de chaque élément — c'est le comportement voulu.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4String HadronicGenerator::WeightCompute(
    const std::vector<std::vector<G4double>> &gammaZData,
    G4int    protonNbBins,
    G4double protonMinEnergy,
    G4double protonMaxEnergy,
    G4double frac)
{
  auto analysisManager = G4AnalysisManager::Instance();


  const G4double dE =
      (protonMaxEnergy - protonMinEnergy) / static_cast<G4double>(protonNbBins);


  for (G4int i = 0; i < protonNbBins; ++i) {


    // Somme de la colonne i sur tous les bins gamma
    // = équivalent de TH2D::Integral(i+1, i+1, 1, NbinsY) dans ROOT
    G4double columnSum = 0.0;
    for (G4double v : gammaZData[i]) columnSum += v;


    // Centre du bin proton correspondant
    G4double binCenter = protonMinEnergy + (i + 0.5) * dE;


    // FillH1 accumule : si plusieurs éléments contribuent,
    // les appels successifs s'additionnent — comportement identique à
    // TH1D_weight->SetBinContent(i, ancien + frac * columnSum)
    analysisManager->FillH1(2, binCenter , frac * columnSum);


    G4cout << "Bin " << i + 1
           << "  centre = " << binCenter / MeV << " MeV"
           << "  weight = " << frac * columnSum
           << G4endl;
  }


  return "is counted, with a fraction of " + std::to_string(frac);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

