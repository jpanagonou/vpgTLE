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
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Constructeur
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GatePromptGammaStatisticActor::GatePromptGammaStatisticActor(py::dict &user_info)
    : GateVActor(user_info, true) {
    // Acteur compatible avec le mode multithread de Geant4
    fMultiThreadReady = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Destructeur
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GatePromptGammaStatisticActor::~GatePromptGammaStatisticActor() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// InitializeUserInfo
// Récupère les paramètres définis côté Python dans user_info_defaults
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GatePromptGammaStatisticActor::InitializeUserInfo(py::dict &user_info) {
    GateVActor::InitializeUserInfo(user_info);

    // Paramètres de binning proton
    fProtonNbBins    = py::int_(user_info["protonNbBins"]);
    fProtonMinEnergy = py::float_(user_info["protonMinEnergy"]);
    fProtonMaxEnergy = py::float_(user_info["protonMaxEnergy"]);

    // Paramètres de binning gamma
    fGammaNbBins     = py::int_(user_info["gammaNbBins"]);
    fGammaMinEnergy  = py::float_(user_info["gammaMinEnergy"]);
    fGammaMaxEnergy  = py::float_(user_info["gammaMaxEnergy"]);

    // Paramètres de sortie
    fWeight          = py::bool_(user_info["weight"]);
    fBodyFraction    = py::float_(user_info["body_fraction"]);
    fOutputFilename  = py::str(user_info["pg_output_filename"]);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// InitializeCpp
// Initialise les objets C++ qui dépendent des paramètres récupérés
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GatePromptGammaStatisticActor::InitializeCpp() {
    GateVActor::InitializeCpp();

    fSigmaFilled = false;
    fNInelastic.assign(fProtonNbBins, 0.0);
    fDensity = 0.0;
    fGammaEvents.clear(); // new

    // Shadow array remplaçant TH2D de ROOT
    // Dimensions : fProtonNbBins x fGammaNbBins
    // Accumule le poids kappa_inel par bin (Ep, Egamma)
    // Utilisé par WeightCompute() en fin de run
    fGammaZData.assign(fProtonNbBins,
                       std::vector<G4double>(fGammaNbBins, 0.0));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// BeginOfRunActionMasterThread
// Initialise G4AnalysisManager, crée les histogrammes et ouvre le fichier .root
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GatePromptGammaStatisticActor::BeginOfRunActionMasterThread(int run_id) {

    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    analysisManager->SetFileName(fOutputFilename);
    analysisManager->SetVerboseLevel(1);

    // H2[0] : EpEpg — spectre 2D brut, rempli avec poids = 1
    analysisManager->CreateH2("EpEpg",
                 "PGs energy vs proton energy;"
                 "Proton energy [MeV];PG energy [MeV];counts",
                 fProtonNbBins, fProtonMinEnergy / MeV, fProtonMaxEnergy / MeV,
                 fGammaNbBins,  fGammaMinEnergy  / MeV, fGammaMaxEnergy  / MeV);

    // H2[1] : GammaZ — spectre 2D pondéré par kappa_inel
    analysisManager->CreateH2("GammaZ",
                 "PG yield weighted by kappa_inel;"
                 "Proton energy [MeV];PG energy [MeV];kappa [cm^{-1} per collision]",
                 fProtonNbBins, fProtonMinEnergy / MeV, fProtonMaxEnergy / MeV,
                 fGammaNbBins,  fGammaMinEnergy  / MeV, fGammaMaxEnergy  / MeV);

    // H1[0] : NrPG — nombre de gammas prompts par bin proton
    analysisManager->CreateH1("NrPG",
                 "Number of prompt gammas per bin;"
                 "Proton energy [MeV];N_{PG}",
                 fProtonNbBins, fProtonMinEnergy / MeV, fProtonMaxEnergy / MeV);

    // H1[1] : Kapa inelastique — coefficient d'atténuation linéaire kappa_inel(Ep)
    analysisManager->CreateH1("Kapa inelastique",
                 "Linear attenuation coefficient;"
                 "Proton energy [MeV];kappa_{inel} [cm^{-1}]",
                 fProtonNbBins, fProtonMinEnergy / MeV, fProtonMaxEnergy / MeV);

    // H1[2] : Weight — rendement PG pondéré pour le calcul ToF (optionnel)
    analysisManager->CreateH1("Weight",
                 "Weighted PG yield for ToF;"
                 "Proton energy [MeV];Weight [cm^{-1}]",
                 fProtonNbBins, fProtonMinEnergy / MeV, fProtonMaxEnergy / MeV);

    analysisManager->OpenFile();

    // Remettre le shadow array à zéro au début de chaque run
    for (auto &row : fGammaZData)
        std::fill(row.begin(), row.end(), 0.0);

    fSigmaFilled = false;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// BeginOfRunAction — rien à faire au début de chaque run worker
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GatePromptGammaStatisticActor::BeginOfRunAction(const G4Run *run) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// SteppingAction
// Détecte les protons, filtre les interactions inélastiques,
// remplit les histogrammes EpEpg, GammaZ, NrPG et le shadow array fGammaZData
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GatePromptGammaStatisticActor::SteppingAction(G4Step *step) {
   // G4cout << "SteppingAction appelé" << G4endl; // Test

    // ── Géométrie des bins 
    const G4double dE          = (fProtonMaxEnergy - fProtonMinEnergy) / fProtonNbBins;
    const G4double dEgamma     = (fGammaMaxEnergy  - fGammaMinEnergy)  / fGammaNbBins;
    const G4double energyStart = fProtonMinEnergy + 0.5 * dE; 

    // Vérifier que la particule est un proton
    const G4ParticleDefinition *particle =
        step->GetTrack()->GetParticleDefinition();
    if (particle != G4Proton::Proton()) return;

    // Énergie du proton avant l'interaction
    // On soustrait le dépôt d'énergie pour obtenir l'énergie vraie
    // avant la dernière ionisation (cf. Gate 9)
    const G4double particle_energy =
        step->GetPreStepPoint()->GetKineticEnergy() -
        step->GetTotalEnergyDeposit();

    // Récupération du matériau et du processus
    const G4Material *material = step->GetPreStepPoint()->GetMaterial();
    const G4VProcess *process  = step->GetPostStepPoint()->GetProcessDefinedStep();

    static G4HadronicProcessStore *store = G4HadronicProcessStore::Instance();
    static G4VProcess *protonInelastic =
        store->FindProcess(G4Proton::Proton(), fHadronInelastic);

    // Filtrer uniquement les interactions inélastiques proton
    if (process != protonInelastic) return;
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    // Récupérer la densité une seule fois
    if (fDensity == 0.0)
    fDensity = step->GetPreStepPoint()->GetMaterial()->GetDensity() / (CLHEP::g / CLHEP::cm3);

    // Incrémenter le compteur d'interactions inélastiques
    G4int pBinInelastic = static_cast<G4int>(
    (particle_energy / MeV - fProtonMinEnergy / MeV) / (dE / MeV));
    if (pBinInelastic >= 0 && pBinInelastic < fProtonNbBins)
    fNInelastic[pBinInelastic] += 1.0;
    //----------------------------------------------------------------------------------------------------------------------------------------------

    // Ignorer les interactions inélastiques qui ne stoppent pas le proton
    if (step->GetPostStepPoint()->GetKineticEnergy() > 0.) return;

    // Section efficace inélastique [cm-1] à l'énergie courante du proton
    G4double cross_section =
        store->GetCrossSectionPerVolume(particle, particle_energy,
                                        process, material);

                                 

    auto analysisManager = G4AnalysisManager::Instance();

    // ── Remplissage de Kapa inelastique (une seule fois) ──────────────────────
    // On évalue kappa_inel au centre exact de chaque bin proton
    // Equivalent de sigma_filled dans Gate 9
    if (!fSigmaFilled) {
        for (G4int f = 1; f <= fProtonNbBins; ++f) {
            G4double binCenter         = energyStart + (f - 1) * dE;
            G4double crossSectionLocal = store->GetCrossSectionPerVolume(
                particle, binCenter, protonInelastic, material);
            analysisManager->FillH1(1, binCenter / MeV, crossSectionLocal * cm);
        }
        fSigmaFilled = true;
    }

    // ── Boucle sur les secondaires ────────────────────────────────────────────
    G4TrackVector *secondaries =
        (const_cast<G4Step *>(step))->GetfSecondary();

    for (size_t i = 0; i < secondaries->size(); ++i) {
        if ((*secondaries)[i]->GetDefinition() != G4Gamma::Gamma()) continue;

        G4double Egamma = (*secondaries)[i]->GetKineticEnergy();

        // Seuil bas : ignorer les gammas de basse énergie (cf. Gate 9 : 40 keV)
        if (Egamma < 0.04 * MeV) continue;
        if (Egamma > fGammaMaxEnergy) continue;

        G4double Ep_MeV = particle_energy / MeV;
        G4double Eg_MeV = Egamma / MeV;

        // H2[0] : EpEpg — remplissage non pondéré (comptage pur)
        analysisManager->FillH2(0, Ep_MeV, Eg_MeV);

        // H2[1] : GammaZ — remplissage pondéré par cross_section
       // analysisManager->FillH2(1, Ep_MeV, Eg_MeV, cross_section * cm);

        // H1[0] : NrPG — comptage événement par événement
        analysisManager->FillH1(0, Ep_MeV, 1.0);

        // ── Shadow array fGammaZData ──────────────────────────────────────────
        // Remplace TH2D::Integral() de ROOT utilisé dans WeightCompute()
        // Index du bin gamma (0-based)
        G4int gBin = static_cast<G4int>(
            (Eg_MeV - fGammaMinEnergy / MeV) / (dEgamma / MeV));

        // Index du bin proton (0-based)
        G4int pBin = static_cast<G4int>(
            (Ep_MeV - fProtonMinEnergy / MeV) / (dE / MeV));

        if (gBin >= 0 && gBin < fGammaNbBins &&
            pBin >= 0 && pBin < fProtonNbBins) {
            fGammaZData[pBin][gBin] += cross_section;  //NNEWWW
            fGammaEvents.push_back({Ep_MeV, Eg_MeV, cross_section, pBin}); 
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// EndOfRunAction — rien à faire à la fin de chaque run worker
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GatePromptGammaStatisticActor::EndOfRunAction(const G4Run *run) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// WeightCompute si activé, puis écriture et fermeture du fichier .root
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int GatePromptGammaStatisticActor::EndOfRunActionMasterThread(int run_id) {

    auto analysisManager = G4AnalysisManager::Instance();


    // ── Normalisation de fGammaZData par N_inel * rho (comme Gate 9 SaveData()) ── 


    for (G4int i = 0; i < fProtonNbBins; ++i) { 
        G4double scalingfactor = fNInelastic[i] * fDensity;
         if (scalingfactor > 0) { 
            for (G4int j = 0; j < fGammaNbBins; ++j) {
                 fGammaZData[i][j] /= scalingfactor;

                }
             }
        }

        // Remplir GammaZ depuis fGammaZData normalisé
    const G4double dE      = (fProtonMaxEnergy - fProtonMinEnergy) / fProtonNbBins;
    const G4double dEgamma = (fGammaMaxEnergy  - fGammaMinEnergy)  / fGammaNbBins;
    const G4double energyStart     = fProtonMinEnergy + 0.5 * dE;
    const G4double gammaEnergyStart = fGammaMinEnergy + 0.5 * dEgamma;
    
    /*
    for (G4int i = 0; i < fProtonNbBins; ++i) {
        G4double binCenterP = energyStart + i * dE;
        for (G4int j = 0; j < fGammaNbBins; ++j) {
            if (fGammaZData[i][j] > 0) {
                G4double binCenterG = gammaEnergyStart + j * dEgamma;
                analysisManager->FillH2(1, binCenterP / MeV, binCenterG / MeV, fGammaZData[i][j]);
        }
    }
}*/



    // New 
    for (auto &[Ep, Eg, cross_section_stored, pBinEvt] : fGammaEvents) {
        G4double scalingfactor = fNInelastic[pBinEvt] * fDensity;
        if (scalingfactor > 0) {
            analysisManager->FillH2(1, Ep, Eg, cross_section_stored / scalingfactor);
        }
    }


    // ── WeightCompute ─────────────────────────────────────────────────────────
    // Remplir H1[2] Weight si l'option est activée
    // Pour chaque bin proton, somme la colonne de fGammaZData
    // et pondère par la fraction corporelle de l'élément
    if (fWeight) {
        const G4double dE          = (fProtonMaxEnergy - fProtonMinEnergy) / fProtonNbBins;
        const G4double energyStart = fProtonMinEnergy + 0.5 * dE;
        for (G4int f = 1; f <= fProtonNbBins; ++f) {
            G4double columnSum = 0.0;
            for (G4double v : fGammaZData[f - 1]) columnSum += v;
            G4double binCenter = energyStart + (f - 1) * dE;
            analysisManager->FillH1(2, binCenter / MeV, fBodyFraction * columnSum);
        }
    }

    // Ecriture et fermeture du fichier .root
    analysisManager->Write();
    analysisManager->CloseFile();

    return 0;
}
