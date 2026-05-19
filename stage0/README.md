# PromptGammaStatisticActor


The `PromptGammaStatisticActor` is a Geant4-based actor implemented in Gate 10 for building a prompt gamma database (PGdb) used in the vpgTLE (variance prompt gamma Track Length Estimator) module. It computes, for each chemical element of atomic number Z, the prompt gamma yield per unit density:


$$\frac{\Gamma_Z(E)}{\rho_Z} = \frac{N_\gamma(Z,E)}{N_{\rm inel}(Z,E)} \cdot \frac{\kappa_{\rm inel}(Z,E)}{\rho_Z}$$


where $N_\gamma(Z,E)$ is the number of prompt gammas produced during inelastic interactions in the incident particle energy bin $E$, $N_{\rm inel}(Z,E)$ the number of total inelastic interactions (i.e. where the incident particle stops, $KE = 0$), and $\kappa_{\rm inel}(Z,E)$ the inelastic linear attenuation coefficient of element $Z$ at energy $E$.


Two simulation modes are available:


- **Mono-element mode**: one simulation per pure element material.
- **Multi-element mode**: one simulation with an `AllElements` composite material containing all elements of interest.


Supported incident particles: `proton`, `neutron`, `helium` ($^4\text{He}$), `carbon` ($^{12}\text{C}$).


> ⚠️ **Warning**
> The proton and neutron modes have been validated by comparison with Gate 9 reference simulations on a pure oxygen target at 200 MeV. The helium and carbon modes are currently under development and have not yet been validated against reference simulations. Use with caution.


---


## Parameters


| Parameter | Default | Description |
|---|---|---|
| `particleNbBins` | 500 | Number of bins for the incident particle energy axis |
| `particleMinEnergy` | 0 MeV | Minimum incident particle energy |
| `particleMaxEnergy` | 200 MeV | Maximum incident particle energy |
| `gammaNbBins` | 250 | Number of bins for the prompt gamma energy axis |
| `gammaMinEnergy` | 0 MeV | Minimum prompt gamma energy |
| `gammaMaxEnergy` | 10 MeV | Maximum prompt gamma energy |
| `pg_output_filename` | `"PGdb"` | Output filename without extension |
| `particle_type` | `"proton"` | Incident particle type: `proton`, `neutron`, `helium`, `carbon` |
| `multi_element` | `False` | If `True`, run a single simulation with `AllElements` composite material |
| `material_name` | `""` | Name of the material. Use the element name in mono-element mode (e.g. `"Oxygen"`), and `"AllElements"` in multi-element mode. |
| `save_KE0_secondaries` | `False` | If `True`, save energies of secondaries produced during partial ($KE > 0$) and total ($KE = 0$) inelastic interactions to ROOT files. Available in mono-element mode only. Automatically set to `False` in multi-element mode with a warning. |


---


## Usage


The following script is the base script for all modes. Only the parameters at the top need to be changed depending on the use case.


```python
import opengate as gate
from pathlib import Path
from opengate.actors.pgactors_utils import (
    reorganize_root_file,
    convert_secondaries_to_root
)


MeV = gate.g4_units.MeV
cm  = gate.g4_units.cm
m   = gate.g4_units.m
mm  = gate.g4_units.mm


# ── Parameters ────────────────────────────────────────────────
particle_type = "proton"    # "proton", "neutron", "helium", "carbon"
material_name = "Oxygen"    # element name or "AllElements"
multi_element = False       # True for multi-element mode
output_dir    = Path("output")
output_dir.mkdir(exist_ok=True)


# ── Simulation ────────────────────────────────────────────────
sim = gate.Simulation()
sim.progress_bar = True
sim.physics_manager.physics_list_name = "QGSP_BIC_HP_EMY"
sim.volume_manager.add_material_database("GateMaterials.db")
sim.volume_manager.add_material_database("GateMaterialsElements.db")


# World
world          = sim.world
world.size     = [3 * m, 3 * m, 3 * m]
world.material = "G4_Galactic"


# Target volume
target          = sim.add_volume("Box", "target")
target.size     = [70 * cm, 70 * cm, 70 * cm]
target.material = material_name


# Step limiter
sim.physics_manager.set_max_step_size("target", 1 * mm)
sim.physics_manager.set_user_limits_particles(["all"])


# Source
source                      = sim.add_source("GenericSource", "source")
source.particle             = "proton"     # see table below
source.energy.mono          = 200 * MeV
source.n                    = 1e7
source.direction.type       = "momentum"
source.direction.momentum   = [1, 0, 0]
source.position.type        = "disc"
source.position.radius      = 4 * mm
source.position.translation = [-40 * cm, 0, 0]


# Actor
actor                      = sim.add_actor(
    "PromptGammaStatisticActor", "pg_actor")
actor.attached_to          = "target"
actor.particle_type        = particle_type
actor.multi_element        = multi_element
actor.material_name        = material_name
actor.particleNbBins       = 500
actor.particleMinEnergy    = 0 * MeV
actor.particleMaxEnergy    = 200 * MeV
actor.gammaNbBins          = 250
actor.gammaMinEnergy       = 0 * MeV
actor.gammaMaxEnergy       = 10 * MeV
actor.pg_output_filename   = str(output_dir /
    f"PGdb_{material_name}_{particle_type}")
actor.save_KE0_secondaries = False


sim.run()


# ── Post-processing ───────────────────────────────────────────
reorganize_root_file(str(output_dir /
    f"PGdb_{material_name}_{particle_type}.root"))


if actor.save_KE0_secondaries:
    convert_secondaries_to_root(output_dir, material_name,
                                 particle_type)
```


The following table summarizes the parameters to change for each use case:


| Mode | `particle_type` | `material_name` | `multi_element` | `source.particle` |
|---|---|---|---|---|
| Mono proton | `"proton"` | `"Oxygen"` | `False` | `"proton"` |
| Mono neutron | `"neutron"` | `"Oxygen"` | `False` | `"neutron"` |
| Mono helium | `"helium"` | `"Oxygen"` | `False` | `"alpha"` |
| Mono carbon | `"carbon"` | `"Oxygen"` | `False` | `"ion 6 12"` |
| Multi-element | `"proton"` | `"AllElements"` | `True` | `"proton"` |


> 📝 **Note**
> For carbon beams, set `particleMaxEnergy = 430 * 12 * MeV` to cover the full clinical energy range (100–430 MeV/u).


---


## Output ROOT file structure


`G4AnalysisManager` produces a flat ROOT file where all histograms are stored at the root level without sub-folders:


```
PGdb_Oxygen_proton.root         (flat output from G4AnalysisManager)
├── Oxygen/EpEpg
├── Oxygen/GammaZ
├── Oxygen/NrPG
├── Oxygen/Kapa inelastique
├── Oxygen/Ep
├── Oxygen/EpInelastic
└── Oxygen/EpInelasticProducedGamma
```


The `reorganize_root_file` function reorganizes these histograms into proper sub-folders per element. This step is particularly important in multi-element mode where all elements are stored in the same flat file and need to be separated into sub-folders for use in the vpgTLE stage 1b:


```
PGdb_Oxygen_proton.root         (after reorganize_root_file — mono)
└── Oxygen/
    ├── EpEpg
    ├── GammaZ
    ├── NrPG
    ├── Kapa inelastique
    ├── Ep
    ├── EpInelastic
    └── EpInelasticProducedGamma


PGdb_AllElements_proton.root    (after reorganize_root_file — multi)
├── Hydrogen/
│   ├── EpEpg
│   ├── GammaZ
│   ├── NrPG
│   ├── Kapa inelastique
│   ├── Ep
│   ├── EpInelastic
│   └── EpInelasticProducedGamma
├── Carbon/
│   ├── EpEpg
│   └── ...
└── Oxygen/
    ├── EpEpg
    └── ...
```


---


## Output histograms


The histogram names depend on the incident particle type:


| Histogram | proton | neutron | helium | carbon | Description |
|---|---|---|---|---|---|
| 2D spectrum | `EpEpg` | `EnEpg` | `EHeEpg` | `ECEpg` | 2D prompt gamma spectrum $(E_{\rm part}, E_\gamma)$ |
| PG yield | `GammaZ` | `GammaZ` | `GammaZ` | `GammaZ` | 2D PG yield weighted by $\kappa_{\rm inel}$ |
| NrPG | `NrPG` | `NrPG` | `NrPG` | `NrPG` | Number of prompt gammas per energy bin |
| Kappa | `Kapa inelastique` | `Kapa inelastique` | `Kapa inelastique` | `Kapa inelastique` | Inelastic linear attenuation coefficient |
| Energy spectrum | `Ep` | `En` | `EHe` | `EC` | Incident particle energy spectrum |
| Inelastic | `EpInelastic` | `EnInelastic` | `EHeInelastic` | `ECInelastic` | Number of inelastic interactions per bin |
| Inelastic+PG | `EpInelasticProducedGamma` | `EnInelasticProducedGamma` | `EHeInelasticProducedGamma` | `ECInelasticProducedGamma` | Inelastic interactions producing at least one PG |


---


## Secondaries ROOT file structure


When `save_KE0_secondaries = True` (mono-element mode only), a separate ROOT file is produced by `convert_secondaries_to_root`:


```
PGdb_Oxygen_helium_secondaries.root
├── KE0_partial/    ← secondaries from partial interactions (KE > 0)
│   ├── e-
│   └── ...
└── KE0_total/      ← secondaries from total interactions (KE = 0)
    ├── proton
    ├── alpha
    ├── C13
    └── ...
```


Each histogram contains the kinetic energy distribution of the corresponding secondary particle in MeV (500 bins, range 0–200 MeV).


> 📝 **Note**
> `save_KE0_secondaries` is not available in multi-element mode. It is automatically set to `False` with a warning when `multi_element = True`. The investigation of secondaries must therefore be performed in mono-element mode on each element of interest separately.


---


## Post-processing


Two utility functions are provided in `opengate/actors/pgactors_utils.py`:


- `reorganize_root_file(filepath)` — reorganizes flat ROOT histograms into sub-folders per element. **Required after every simulation.**
- `convert_secondaries_to_root(output_dir, material_name, particle_type)` — converts secondaries txt files to a single ROOT file. Required only when `save_KE0_secondaries = True`.


```python
from opengate.actors.pgactors_utils import (
    reorganize_root_file,
    convert_secondaries_to_root
)


# Required after every simulation
reorganize_root_file("output/PGdb_Oxygen_proton.root")


# Required only when save_KE0_secondaries = True
convert_secondaries_to_root(
    Path("output"), "Oxygen", "helium")
```


---


## Physics list recommendation


The recommended physics list is `QGSP_BIC_HP_EMY` for all particle types. A step limiter of 1 mm applied to all particles ensures that the energy variation per step remains smaller than the bin width:


```python
sim.physics_manager.physics_list_name = "QGSP_BIC_HP_EMY"
sim.physics_manager.set_max_step_size("target", 1 * mm)
sim.physics_manager.set_user_limits_particles(["all"])
```


---


## Notes


- `material_name` and `target.material` must always be set to the same value.
- The `AllElements` composite material must be defined in the material database file passed to Gate.
- For helium beams, use `source.particle = "alpha"`.
- For carbon beams, use `source.particle = "ion 6 12"`.
- `save_KE0_secondaries` is available in mono-element mode only. In multi-element mode, it is automatically set to `False` with a warning.
- The `reorganize_root_file` post-processing function is required after every simulation.
- The `convert_secondaries_to_root` post-processing function is required only when `save_KE0_secondaries = True`.


---


## References


- El Kanawati et al. (2015). *Monte Carlo simulation of prompt gamma-ray emission in proton therapy using a specific track length estimator*. Physics in Medicine and Biology, 60(20), 8067. https://doi.org/10.1088/0031-9155/60/20/8067


- Huisman et al. (2016). *Accelerated prompt gamma estimation for clinical proton therapy simulations*. Physics in Medicine and Biology, 61(21), 7725. https://doi.org/10.1088/0031-9155/61/21/7725


- Létang et al. (2024). *Prompt-gamma track-length estimator with time tagging from proton tracking*. Physics in Medicine & Biology, 69(11), 115052. https://doi.org/10.1088/1361-6560/AD4A01



