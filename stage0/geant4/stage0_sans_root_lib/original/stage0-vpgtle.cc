#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Random/Ranlux64Engine.h"
#include "G4GenericIon.hh"
#include "G4HadronicParameters.hh"
#include "G4IonTable.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4ios.hh"
#include "HadronicGenerator-vpgtle.hh"
#include "globals.hh"

#include <iostream>

int main()
{
    // ===============================
    // PHYSICS
    // ===============================
    const G4String namePhysics = "QGSP_BIC_HP";
    auto generator = new HadronicGenerator(namePhysics);

    // ===============================
    // CHOIX DU MATERIAU
    // ===============================

    // G4String materialName = "G4_H";
    // G4String materialName = "G4_C";
    G4String materialName = "G4_O";   // 
    // G4String materialName = "G4_N";
    // G4String materialName = "G4_Ca";
    // G4String materialName = "G4_P";
    // G4String materialName = "G4_K";
    // G4String materialName = "G4_S";
    // G4String materialName = "G4_Na";
    // G4String materialName = "G4_Cl";
    // G4String materialName = "G4_Mg";
    // G4String materialName = "G4_Si";
    // G4String materialName = "G4_Al";
    // G4String materialName = "G4_Fe";
    // G4String materialName = "G4_Cu";
    // G4String materialName = "G4_Zn";

    // ===============================
    // NOM DU FICHIER ROOT
    // ===============================
    G4String fileName = "PG_" + materialName;

    // ===============================
    // MATERIAL
    // ===============================
    G4Material* material =
        G4NistManager::Instance()->FindOrBuildMaterial(materialName);

    if (!material)
    {
        G4cerr << "Material not found !" << G4endl;
        return 1;
    }

    // ===============================
    // RANDOM
    // ===============================
    CLHEP::Ranlux64Engine engine(1234567, 4);
    CLHEP::HepRandom::setTheEngine(&engine);
    CLHEP::HepRandom::setTheSeed(time(NULL));

    // ===============================
    // PARAMETRES
    // ===============================

    G4int nbCollisions = 1e3;

    G4int protonNbBins = 500;
    G4double protonMinEnergy = 0.0 * MeV;
    G4double protonMaxEnergy = 200.0 * MeV;

    G4int gammaNbBins = 250;
    G4double gammaMinEnergy = 0.0 * MeV;
    G4double gammaMaxEnergy = 10.0 * MeV;

    // ===============================
    // ANALYSIS
    // ===============================
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetFileName(fileName);

    generator->InitAnalysis(fileName,
        protonNbBins,
        protonMinEnergy,
        protonMaxEnergy,
        gammaNbBins,
        gammaMinEnergy,
        gammaMaxEnergy);

    // ===============================
    // PROJECTILE
    // ===============================
    auto particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* projectile =
        particleTable->FindParticle("proton");

    // ===============================
    // SIMULATION
    // ===============================
    G4int nbPG = generator->GenerateInteraction(
        projectile,
        material,
        nbCollisions,
        protonNbBins,
        protonMinEnergy,
        protonMaxEnergy);

    G4cout << "Material: " << materialName << G4endl;
    G4cout << "Number of PG: " << nbPG << G4endl;

    // ===============================
    // SAVE
    // ===============================
    generator->EndAnalysis();

    delete generator;

    return 0;
}
