#ifndef HadronicGenerator_h
#define HadronicGenerator_h

#include "G4HadronicProcess.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <map>
#include <vector>

// ================= FORWARD DECLARATIONS =================
class G4ParticleDefinition;
class G4VParticleChange;
class G4ParticleTable;
class G4Material;
class G4HadronicInteraction;
class G4DynamicParticle;
class G4Element;
class G4Track;
class G4StepPoint;
// ========================================================

class HadronicGenerator {

public:
  explicit HadronicGenerator(const G4String physicsCase = "QGSP_BIC_HP");
  ~HadronicGenerator();

  // ===== ANALYSIS =====
  void InitAnalysis(const G4String   &fileName,
                    G4int             protonNbBins,
                    G4double          protonMinEnergy,
                    G4double          protonMaxEnergy,
                    G4int             gammaNbBins,
                    G4double          gammaMinEnergy,
                    G4double          gammaMaxEnergy);

  void EndAnalysis();

  // ===== PHYSIQUE =====
  G4int GenerateInteraction(G4ParticleDefinition*                    projectileDefinition,
                            G4Material*                              targetMaterial,
                            G4int                                    nbCollisions,
                            G4int                                    protonNbBins,
                            G4double                                 protonMinEnergy,
                            G4double                                 protonMaxEnergy,
                            G4int                                    gammaNbBins,
                            G4double                                 gammaMinEnergy,
                            G4double                                 gammaMaxEnergy,
                            std::vector<std::vector<G4double>>      &gammaZData);

  // ===== UTILS =====
  inline G4HadronicProcess    *GetHadronicProcess()     const;
  inline G4HadronicInteraction *GetHadronicInteraction() const;

  G4double GetCrossSection(const G4DynamicParticle *part,
                           const G4Element         *elm,
                           const G4Material        *targetMaterial);

  G4String WeightCompute(const std::vector<std::vector<G4double>> &gammaZData,
                       G4int    protonNbBins,
                       G4double protonMinEnergy,
                       G4double protonMaxEnergy,
                       G4double frac);

private:
  // ===== PHYSICS CONFIG =====
  G4String fPhysicsCase;

  // ===== INTERNAL =====
  G4HadronicProcess *fLastHadronicProcess;
  G4ParticleTable   *fPartTable;

  std::map<G4ParticleDefinition *, G4HadronicProcess *> fProcessMap;
  std::map<G4ParticleDefinition *, G4HadronicProcess *> fProcessMapHP;

  // ===== TRACKING =====
  G4Track     *gTrack;
  G4StepPoint *aPoint;
};

// ========================================================

inline G4HadronicProcess *HadronicGenerator::GetHadronicProcess() const {
  return fLastHadronicProcess;
}

inline G4HadronicInteraction *HadronicGenerator::GetHadronicInteraction() const {
  return fLastHadronicProcess == nullptr
             ? nullptr
             : fLastHadronicProcess->GetHadronicInteraction();
}

#endif
