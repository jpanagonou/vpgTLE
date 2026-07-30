# VoxelizedPromptGammaTLESource

Source voxélisée pour la génération de gammas prompts (PG) en protonthérapie, à partir d'une carte 4D de rendement (trois dimensions spatiales + une dimension spectrale en énergie). Cette carte est produite en amont par le module `vpgTLE`, qui estime le rendement de production de prompt-gammas par la méthode du Track Length Estimator (TLE)~: plutôt que de simuler explicitement chaque interaction proton-noyau susceptible d'émettre un PG (approche dite analogique), le TLE accumule, tout au long de la trajectoire de chaque proton, la contribution attendue au rendement de PG proportionnellement à la longueur parcourue dans chaque voxel, réduisant ainsi le bruit statistique de la carte obtenue à nombre de protons primaires égal. `VoxelizedPromptGammaTLESource` constitue le portage sous Gate 10 de la source historique `GateSourceOfPromptGamma`, utilisée sous Gate 9 pour exploiter ce type de carte.

Pour chaque événement, un voxel est tiré proportionnellement à son rendement total, une position aléatoire est générée à l'intérieur de ce voxel, puis une énergie est tirée selon le spectre propre à ce même voxel. La position et l'énergie sont ainsi couplées via le voxel sélectionné.

## Utilisation

```python
import opengate as gate

sim = gate.Simulation()

source = sim.add_source("VoxelizedPromptGammaTLESource", "pg_source")
source.attached_to = "my_volume"
source.particle = "gamma"
source.image = "path/to/pg_image_4d.nii.gz"
source.Emin = 0.0        # MeV, borne basse du spectre en energie
source.Emax = 10.0       # MeV, borne haute du spectre en energie
source.nbBins = 250      # nombre de bins d'energie dans l'image
source.n_protons = 1e6   # nombre de protons primaires equivalents
source.direction.type = "iso"
```

## Paramètres

| Paramètre | Type | Description |
|---|---|---|
| `image` | `str` | Chemin vers l'image 4D (`.nii.gz`), de forme `[nbBins, Z, Y, X]` |
| `Emin` | `float` | Borne basse du spectre en énergie (MeV) |
| `Emax` | `float` | Borne haute du spectre en énergie (MeV) |
| `nbBins` | `int` | Nombre de bins d'énergie de l'image |
| `n_protons` | `float` | Nombre de protons primaires équivalents. Le nombre de gammas réellement générés est `n_protons × yield_per_proton`, où `yield_per_proton` est la somme totale de l'image, normalisée par proton |
| `direction.type` | `str` | Type de distribution angulaire (hérité de `GenericSource`), typiquement `"iso"` pour une émission isotrope |

Le paramètre est exprimé en nombre de **protons**, et non directement en nombre de gammas, car la carte de rendement produite par le TLE encode le nombre de PG attendu **par proton primaire** simulé lors de l'étape d'estimation. Fixer directement un nombre de gammas obligerait l'utilisateur à connaître au préalable le rendement total de la carte pour en déduire le nombre de protons correspondants, alors que `n_protons` correspond directement à la grandeur physique pertinente~: le nombre de protons délivrés par le faisceau que l'on souhaite simuler.

## Format de l'image d'entrée

L'image doit être un fichier 4D lisible par ITK (`.nii.gz` recommandé), de forme `[nbBins, Z, Y, X]` en convention numpy. Chaque voxel `(iz, iy, ix)` contient un spectre en énergie sur `nbBins` valeurs, réparties uniformément entre `Emin` et `Emax`~; la valeur stockée dans chaque bin représente un **rendement de gammas prompts par proton primaire** (yield/proton), tel qu'estimé par le TLE. Les voxels dont le rendement total est nul sont automatiquement ignorés (stockage creux), ce qui permet de traiter efficacement des images où la grande majorité des voxels sont vides.

## Note de conception : direction d'émission

La direction d'émission suit le mécanisme générique standard de Gate 10 (`GateSPSAngDistribution`), configuré via `source.direction.type`, et n'est pas couplée à la position ou à l'énergie du voxel sélectionné. Ce choix n'est pas une limitation du portage~: il reflète le comportement physique attendu, l'émission d'un prompt-gamma étant isotrope indépendamment du voxel d'origine, un comportement partagé par l'implémentation historique Gate 9 (`SampleRandomDirection()`, `G4SPSAngDistribution` en mode isotrope).

## Limitations connues

- **Temps de vol non géré** : contrairement à l'implémentation historique Gate 9 (`GateSourceOfPromptGamma`), cette source ne gère pas de dimension temporelle. Seule l'énergie est tirée par voxel.

## Validation

La source a été validée par quatre tests unitaires sur images synthétiques (voxel unique, voxels multiples espacés, voxels adjacents, voxels aux rendements contrastés et spectres non triviaux), ainsi que par comparaison directe avec l'implémentation historique Gate 9 sur une image de prompt-gamma réelle, montrant un accord statistique remarquable sur le spectre en énergie et la position de détection (tests de Kolmogorov-Smirnov, p-values > 0.57).

## Performance

À nombre de protons primaires égal, cette implémentation présente un temps d'initialisation supérieur à Gate 9 (12.5~s contre 8.5~s), mais un temps de calcul par proton inférieur d'un facteur~2.6 (0.70 contre 1.82 microsecondes), la rendant globalement plus rapide dès que le nombre de protons simulés dépasse quelques millions. Ces valeurs sont issues d'un ajustement linéaire sur des mesures réalisées entre $10^5$ et $10^9$ protons primaires.

Cette différence s'explique probablement par deux facteurs distincts~: un temps d'initialisation plus élevé côté Gate 10, lié au calcul des CDF (spatiale et spectrales) avant le tirage, et un temps par événement plus faible, lié à l'absence de dépendance à ROOT pour le tirage d'énergie. Ces pistes restent à confirmer par un profilage détaillé.

## Architecture interne

| Fichier | Statut | Rôle |
|---|---|---|
| `opengate/sources/voxelsources.py` | **Nouveau** | Classe Python `VoxelizedPromptGammaTLESource`, orchestre la lecture de l'image et la transmission au C++ |
| `GateVoxelizedPromptGammaTLESource.h/.cpp` | **Nouveau** | Classe C++ principale, hérite de `GateVoxelSource` |
| `GateSPSEneDistributionVoxelizedPG.h/.cpp` | **Nouveau** | Générateur d'énergie par voxel, stockage creux (`std::unordered_map`), tirage par inversion de CDF |
| `GateSPSVoxelsPosDistribution.h/.cpp` | Ajout sur code existant | Extension par trois membres (`fLastIndexX/Y/Z`) et l'accesseur `GetLastVoxelIndices()`, mémorisant l'indice du dernier voxel tiré. Le mécanisme de tirage lui-même (`VGenerateOne()`) n'est pas modifié |
| `GateSingleParticleSource.h/.cpp` | Ajout sur code existant | Ajout de la méthode `SetEneGenerator()`, symétrique de `SetPosGenerator()` déjà existante, permettant d'injecter un générateur d'énergie personnalisé |
