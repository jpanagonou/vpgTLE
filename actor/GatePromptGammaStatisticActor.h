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
#include "G4Gamma.hh"
#include "G4Step.hh"
#include "G4Run.hh"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <tuple> 

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

    // ── Paramètres de binning ─────────────────────────────────────────
    G4int    fProtonNbBins;
    G4double fProtonMinEnergy;
    G4double fProtonMaxEnergy;
    G4int    fGammaNbBins;
    G4double fGammaMinEnergy;
    G4double fGammaMaxEnergy;
  

    // ── Paramètres de sortie ──────────────────────────────────────────
    G4String fOutputFilename;
    G4bool   fWeight;
    G4double fBodyFraction;

    // ── Flags internes ────────────────────────────────────────────────
    G4bool fSigmaFilled;

    // ── Shadow array pour WeightCompute ──────────────────────────────
    // Remplace TH2D de ROOT — dimensions : fProtonNbBins x fGammaNbBins
    std::vector<std::vector<G4double>> fGammaZData;
    std::vector<G4double> fNInelastic;
    G4double fDensity;
    std::vector<std::tuple<G4double, G4double, G4double, G4int>> fGammaEvents;
};

#endif // GatePromptGammaStatisticActor_h
