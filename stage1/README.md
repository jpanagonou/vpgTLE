# Stage 1b - Optimisation de la conversion des longueurs de trace en spectres d'émission

Le Stage 1b convertit les longueurs de trace calculées au Stage 1a en spectres d'émission de prompt-gammas, voxel par voxel, par convolution avec la base de données PGdb (Stage 0). L'implémentation originale, fondée sur une boucle séquentielle sur l'ensemble des voxels du volume CT, engendre des temps de calcul prohibitifs aux résolutions fines requises pour les applications cliniques.

## Optimisation retenue

Le nombre de matériaux distincts $N_\text{mat}$ dans un CT patient étant très inférieur au nombre de voxels $N_\text{vox}$, la boucle sur les voxels a été remplacée par une boucle sur les matériaux~: le spectre PG est calculé une seule fois par matériau, puis appliqué simultanément à l'ensemble des voxels partageant ce matériau par un produit matriciel (`np.dot`), réduisant la complexité algorithmique de $\mathcal{O}(N_\text{vox})$ à $\mathcal{O}(N_\text{mat})$.

Une seconde approche, restructurant la base de données en un tenseur 3D contracté sur l'ensemble des voxels, a également été testée mais n'apporte pas de gain supplémentaire par rapport à l'approche par produit matriciel, plus simple à implémenter et retenue pour l'intégration dans le module vpgTLE.

## Gains de performance

| Voxels | $N_\text{mat}$ | $T_\text{original}$ (s) | $T_\text{optimisé}$ (s) | Gain |
|---|---|---|---|---|
| $13\times13\times19$ | 32 | 3,16 | 2,00 | ×1,6 |
| $25\times25\times37$ | 49 | 5,22 | 2,40 | ×2,2 |
| $50\times50\times75$ | 52 | 20,21 | 4,72 | ×4,3 |
| $104\times104\times152$ | 32 | 137,53 | 18,59 | ×7,4 |

Le gain croît avec la résolution du CT, le rapport $N_\text{vox}/N_\text{mat}$ augmentant à mesure que la grille se raffine.

## Validation

La différence relative entre les spectres PG produits par l'implémentation originale et par l'implémentation optimisée est de l'ordre de $10^{-8}$~%, sur un voxel du plateau de la courbe de Bragg comme sur un voxel au pic de Bragg, conforme à la précision machine, confirmant l'absence de biais numérique introduit par la vectorisation.
