# PromptGammaStatisticActor

Le `PromptGammaStatisticActor` est un acteur basé sur Geant4, implémenté dans Gate 10 pour construire une base de données de prompt-gamma (PGdb) utilisée par le module vpgTLE (variance prompt gamma Track Length Estimator). Il calcule, pour chaque élément chimique de numéro atomique Z, le rendement de prompt-gamma par unité de densité~:

$$\frac{\Gamma_Z(E)}{\rho_Z} = \frac{N_\gamma(Z,E)}{N_{\rm inel}(Z,E)} \cdot \frac{\kappa_{\rm inel}(Z,E)}{\rho_Z}$$

où $N_\gamma(Z,E)$ est le nombre de prompt-gammas produits lors d'interactions inélastiques dans le bin d'énergie $E$ de la particule incidente, $N_{\rm inel}(Z,E)$ le nombre total d'interactions inélastiques (c'est-à-dire où la particule incidente s'arrête, $KE = 0$), et $\kappa_{\rm inel}(Z,E)$ le coefficient d'atténuation linéique inélastique de l'élément Z à l'énergie E.

Deux modes de simulation sont disponibles~:

- **Mode mono-élément** : une simulation par matériau élémentaire pur.
- **Mode multi-élément** : une seule simulation avec un matériau composite `AllElements` regroupant tous les éléments d'intérêt.

Particules incidentes supportées~: `proton`, `neutron`.

---

## Paramètres

| Paramètre | Valeur par défaut | Description |
|---|---|---|
| `particleNbBins` | 500 | Nombre de bins pour l'axe d'énergie de la particule incidente |
| `particleMinEnergy` | 0 MeV | Énergie minimale de la particule incidente |
| `particleMaxEnergy` | 200 MeV | Énergie maximale de la particule incidente |
| `gammaNbBins` | 250 | Nombre de bins pour l'axe d'énergie du prompt-gamma |
| `gammaMinEnergy` | 0 MeV | Énergie minimale du prompt-gamma |
| `gammaMaxEnergy` | 10 MeV | Énergie maximale du prompt-gamma |
| `pg_output_filename` | `"PGdb"` | Nom du fichier de sortie, sans extension |
| `particle_type` | `"proton"` | Type de particule incidente~: `proton`, `neutron` |
| `multi_element` | `False` | Si `True`, exécute une seule simulation avec le matériau composite `AllElements` |
| `material_name` | `""` | Nom du matériau. Utiliser le nom de l'élément en mode mono-élément (par exemple `"Oxygen"`), et `"AllElements"` en mode multi-élément. |
| `save_KE0_secondaries` | `False` | Si `True`, sauvegarde les énergies des secondaires produits lors des interactions inélastiques partielles ($KE > 0$) et totales ($KE = 0$) dans des fichiers ROOT. Disponible uniquement en mode mono-élément. Automatiquement mis à `False` en mode multi-élément, avec un avertissement. |

---

## Utilisation

Le script suivant constitue la base commune à tous les modes. Seuls les paramètres en tête de script doivent être modifiés selon le cas d'usage.

```python
import opengate as gate
from pathlib import Path
import time
from opengate.actors.pgactors_utils import (
    reorganize_root_file,
    convert_secondaries_to_root
)

MeV = gate.g4_units.MeV
cm  = gate.g4_units.cm
m   = gate.g4_units.m
mm  = gate.g4_units.mm

t0 = time.time()

# ── Paramètres ────────────────────────────────────────────────
particle_type = "proton"       # "proton", "neutron"
material_name = "AllElements"  # AllElements ou Oxygen, Carbon, etc.
data_dir   = Path(__file__).parent.parent / "data"
output_dir = Path(__file__).parent.parent / "output"
output_dir.mkdir(exist_ok=True)

# ── Simulation ────────────────────────────────────────────────
sim = gate.Simulation()

sim.g4_verbose = True
sim.g4_verbose_level = 1
sim.progress_bar = True
sim.physics_manager.physics_list_name = "QGSP_BIC_HP_EMY"
sim.number_of_threads = 1

sim.random_engine = "MersenneTwister"
sim.random_seed = 123456

sim.volume_manager.add_material_database(str(data_dir / "GateMaterials.db"))
sim.volume_manager.add_material_database(str(data_dir / "GateMaterialsElements.db"))

# ── Monde ─────────────────────────────────────────────────────
world          = sim.world
world.size     = [3 * m, 3 * m, 3 * m]
world.material = "G4_Galactic"

target          = sim.add_volume("Box", "target")
target.size     = [70 * cm, 70 * cm, 70 * cm]
target.material = material_name
target.mother   = "world"

# ── Limiteur de pas ───────────────────────────────────────────
sim.physics_manager.set_max_step_size("target", 1 * mm)
sim.physics_manager.set_user_limits_particles(["all"])

# ── Source ────────────────────────────────────────────────────
source                      = sim.add_source("GenericSource", "source")
source.particle             = particle_type  # proton, neutron
source.energy.mono          = 200 * MeV
source.n                    = 1e7
source.direction.type       = "momentum"
source.direction.momentum   = [1, 0, 0]
source.position.type        = "disc"
source.position.radius      = 4 * mm
source.position.translation = [-40 * cm, 0, 0]

# ── Acteur ────────────────────────────────────────────────────
actor                      = sim.add_actor("PromptGammaStatisticActor", "pg_actor")
actor.attached_to          = "target"
actor.particle_type        = particle_type
actor.multi_element        = True
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

# ── Lancement ─────────────────────────────────────────────────
sim.run()

t1 = time.time()
print(f"Temps de simulation : {t1 - t0:.2f}s")

# ── Post-traitement ───────────────────────────────────────────
root_file = str(output_dir / f"PGdb_{material_name}_{particle_type}.root")
reorganize_root_file(root_file)

if actor.save_KE0_secondaries:
    convert_secondaries_to_root(output_dir, material_name, particle_type)

print(f"Simulation terminée — PGdb_{material_name}_{particle_type}.root")
```

Le tableau suivant résume les paramètres à modifier selon le cas d'usage~:

| Mode | `particle_type` | `material_name` | `multi_element` | `source.particle` |
|---|---|---|---|---|
| Mono proton | `"proton"` | `"Oxygen"` | `False` | `"proton"` |
| Mono neutron | `"neutron"` | `"Oxygen"` | `False` | `"neutron"` |
| Multi-élément | `"proton"` | `"AllElements"` | `True` | `"proton"` |

---

## Structure du fichier ROOT de sortie

`G4AnalysisManager` produit un fichier ROOT plat, où tous les histogrammes sont stockés à la racine, sans sous-dossiers~:

```
PGdb_Oxygen_proton.root         (sortie brute de G4AnalysisManager)
├── Oxygen/EpEpg
├── Oxygen/GammaZ
├── Oxygen/Kapa inelastique
├── Oxygen/Ep
├── Oxygen/EpInelastic
└── Oxygen/EpInelasticProducedGamma
```

La fonction `reorganize_root_file` réorganise ces histogrammes en sous-dossiers propres à chaque élément. Cette étape est particulièrement importante en mode multi-élément, où tous les éléments sont stockés dans le même fichier plat et doivent être séparés en sous-dossiers pour être utilisés au stage 1b du vpgTLE~:

```
PGdb_Oxygen_proton.root         (apres reorganize_root_file — mono)
└── Oxygen/
    ├── EpEpg
    ├── GammaZ
    ├── Kapa inelastique
    ├── Ep
    ├── EpInelastic
    └── EpInelasticProducedGamma

PGdb_AllElements_proton.root    (apres reorganize_root_file — multi)
├── Hydrogen/
│   ├── EpEpg
│   ├── GammaZ
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

## Histogrammes de sortie

Les noms d'histogrammes dépendent du type de particule incidente~:

| Histogramme | proton | neutron | Description |
|---|---|---|---|
| Spectre 2D | `EpEpg` | `EnEpg` | Spectre 2D de prompt-gamma $(E_{\rm part}, E_\gamma)$ |
| Rendement PG | `GammaZ` | `GammaZ` | Rendement PG 2D, pondéré par $\kappa_{\rm inel}$ |
| Kappa | `Kapa inelastique` | `Kapa inelastique` | Coefficient d'atténuation linéique inélastique |
| Spectre en énergie | `Ep` | `En` | Spectre en énergie de la particule incidente |
| Inélastique | `EpInelastic` | `EnInelastic` | Nombre d'interactions inélastiques par bin |
| Inélastique+PG | `EpInelasticProducedGamma` | `EnInelasticProducedGamma` | Interactions inélastiques produisant au moins un PG |

---

## Structure du fichier ROOT des secondaires

Lorsque `save_KE0_secondaries = True` (mode mono-élément uniquement), un fichier ROOT séparé est produit par `convert_secondaries_to_root`~:

```
PGdb_Oxygen_proton_secondaries.root
├── KE0_partial/    ← secondaires issus d'interactions partielles (KE > 0)
│   ├── e-
│   └── ...
└── KE0_total/      ← secondaires issus d'interactions totales (KE = 0)
    ├── proton
    ├── alpha
    └── ...
```

Chaque histogramme contient la distribution en énergie cinétique de la particule secondaire correspondante, en MeV (500 bins, plage 0–200 MeV).

> 📝 **Remarque**
> `save_KE0_secondaries` n'est pas disponible en mode multi-élément. Il est automatiquement mis à `False`, avec un avertissement, lorsque `multi_element = True`. L'étude des secondaires doit donc être réalisée en mode mono-élément, séparément pour chaque élément d'intérêt.

---

## Post-traitement

Deux fonctions utilitaires sont fournies dans `opengate/actors/pgactors_utils.py`~:

- `reorganize_root_file(filepath)` — réorganise les histogrammes ROOT plats en sous-dossiers par élément. **Requise après chaque simulation.**
- `convert_secondaries_to_root(output_dir, material_name, particle_type)` — convertit les fichiers texte de secondaires en un unique fichier ROOT. Requise uniquement lorsque `save_KE0_secondaries = True`.

```python
from opengate.actors.pgactors_utils import (
    reorganize_root_file,
    convert_secondaries_to_root
)

# Requise apres chaque simulation
reorganize_root_file("output/PGdb_Oxygen_proton.root")

# Requise uniquement si save_KE0_secondaries = True
convert_secondaries_to_root(
    Path("output"), "Oxygen", "proton")
```

---

## Recommandation de liste de physique

La liste de physique recommandée est `QGSP_BIC_HP_EMY` pour les deux types de particules. Un limiteur de pas de 1 mm, appliqué à toutes les particules, garantit que la variation d'énergie par pas reste inférieure à la largeur d'un bin~:

```python
sim.physics_manager.physics_list_name = "QGSP_BIC_HP_EMY"
sim.physics_manager.set_max_step_size("target", 1 * mm)
sim.physics_manager.set_user_limits_particles(["all"])
```

---

## Remarques

- `material_name` et `target.material` doivent toujours être fixés à la même valeur.
- Le matériau composite `AllElements` doit être défini dans le fichier de base de données de matériaux transmis à Gate.
- `save_KE0_secondaries` n'est disponible qu'en mode mono-élément. En mode multi-élément, il est automatiquement mis à `False`, avec un avertissement.
- La fonction de post-traitement `reorganize_root_file` est requise après chaque simulation.
- La fonction de post-traitement `convert_secondaries_to_root` n'est requise que lorsque `save_KE0_secondaries = True`.

---

## Références

- El Kanawati et al. (2015). *Monte Carlo simulation of prompt gamma-ray emission in proton therapy using a specific track length estimator*. Physics in Medicine and Biology, 60(20), 8067. https://doi.org/10.1088/0031-9155/60/20/8067

- Huisman et al. (2016). *Accelerated prompt gamma estimation for clinical proton therapy simulations*. Physics in Medicine and Biology, 61(21), 7725. https://doi.org/10.1088/0031-9155/61/21/7725

- Létang et al. (2024). *Prompt-gamma track-length estimator with time tagging from proton tracking*. Physics in Medicine & Biology, 69(11), 115052. https://doi.org/10.1088/1361-6560/AD4A01
