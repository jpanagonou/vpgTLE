# vpgTLE - Pipeline de simulation des prompt-gammas sous Gate 10

Ce dossier regroupe les trois étapes du portage sous Gate 10 du module vpgTLE (voxelized prompt gamma Track Length Estimator), permettant de simuler l'émission de prompt-gammas en protonthérapie sans avoir à simuler explicitement chaque interaction nucléaire proton-noyau.

## Stage 0 - Construction de la base de données de prompt-gammas

Le `PromptGammaStatisticActor` construit, pour chaque élément chimique, le rendement de prompt-gamma par unité de densité en fonction de l'énergie de la particule incidente (proton ou neutron), à partir de simulations sur cibles élémentaires. Cette base de données (PGdb) constitue l'entrée du Stage 1.

Voir [`stage0/README.md`](./stage0/README.md) pour la documentation complète.

## Stage 1 - Conversion des longueurs de trace en spectres d'émission

Le Stage 1 combine la base de données PGdb (Stage 0) avec la trajectoire des protons simulée dans le volume CT, pour produire, voxel par voxel, le spectre en énergie des prompt-gammas émis (méthode du Track Length Estimator). Il se décompose en un Stage 1a, qui calcule les longueurs de trace par voxel et par énergie, et un Stage 1b, qui convertit ces longueurs de trace en spectres d'émission par convolution avec la base de données PGdb.

Seul le Stage 1b a fait l'objet d'un travail d'optimisation dans ce projet, remplaçant la boucle séquentielle sur les voxels par un regroupement par matériau (`np.dot`), exploitant le fait que le nombre de matériaux distincts dans un CT est très inférieur au nombre de voxels. Ce regroupement apporte un gain de performance croissant avec la résolution du CT (jusqu'à ×7,4 à 4~mm de spacing sur un volume de $104\times104\times152$ voxels), sans introduire de biais numérique sur les spectres produits.

Voir [`stage1/README.md`](./stage1/README.md) pour la documentation complète.

## Stage 2 - Source d'émission voxélisée sous Gate 10

`VoxelizedPromptGammaTLESource` exploite la carte 4D de spectres produite au Stage 1 pour générer, à chaque événement, un prompt-gamma dont la position et l'énergie sont tirées conjointement, proportionnellement au rendement local de chaque voxel.

Voir [`stage2/README.md`](./stage2/README.md) pour la documentation complète, incluant les résultats de validation (comparaison avec l'implémentation historique Gate 9) et l'étude de performance.

## Références

- El Kanawati et al. (2015). *Monte Carlo simulation of prompt gamma-ray emission in proton therapy using a specific track length estimator*. Physics in Medicine and Biology, 60(20), 8067. https://doi.org/10.1088/0031-9155/60/20/8067
- Huisman et al. (2016). *Accelerated prompt gamma estimation for clinical proton therapy simulations*. Physics in Medicine and Biology, 61(21), 7725. https://doi.org/10.1088/0031-9155/61/21/7725
- Létang et al. (2024). *Prompt-gamma track-length estimator with time tagging from proton tracking*. Physics in Medicine & Biology, 69(11), 115052. https://doi.org/10.1088/1361-6560/AD4A01
