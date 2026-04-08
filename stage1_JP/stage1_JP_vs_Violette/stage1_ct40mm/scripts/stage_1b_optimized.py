#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import itk
import opengate as gate
import os
import time
import numpy as np
from variables import tle_sim
from stage_1a import *


start_total_wall = time.time()

sim = simulation(paths,
    File_name,
    0,
    number_of_particles,
    False,
    False,
    actor,
    Erange,
)

tle_simu = tle_sim(sim)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 : RECUPERER LES OBJETS UTILES
# ═══════════════════════════════════════════════════════════════════════════════

# Chargement de l'image CT et conversion en tableau Numpy
# UH_array contient les valeurs de HU de chaque voxel, shape (Z, Y, X)
ct_object  = tle_simu.ct
UH_array   = itk.array_from_image(ct_object.load_input_image())

# Chemins, acteur et liste des sorties activées (prot_e, neutr_e, ...)
path       = tle_simu.path
actor      = tle_simu.act
actor_list = tle_simu.actor_list

# Bases de données ROOT pour protons et neutrons
data_protons  = tle_simu.root_file
data_neutrons = tle_simu.root_file_neutron

# Liste des matériaux sous forme de triplets [HU_min, HU_max, nom_matériau]
mat = tle_simu.voxel_mat


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 : CONSTRUIRE LES GAMMA PAR MATERIAU
# ═══════════════════════════════════════════════════════════════════════════════

# Vérifie si les neutrons sont activés
neutr = "neutr_e.nii.gz" in actor_list

# Pour chaque matériau présent dans le CT, on calcule sa matrice Gamma (500x250)
# qui représente l'émission gamma en fonction de l'énergie des particules primaires
gamma_neutron = {}
gamma_proton  = {}
for data in mat:
    name = data[2]
    gamma_proton[name] = tle_simu.gamma_mat(name, data_protons)
    if neutr:
        gamma_neutron[name] = tle_simu.gamma_mat(name, data_neutrons)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 : INDEX DE CORRESPONDANCE BASE DE DONNEES → BINS ACTEUR
# ═══════════════════════════════════════════════════════════════════════════════

# Bins d'énergie de l'acteur : 0 → 130 MeV en 250 intervalles (251 bornes)
# actor.energybins = 250 (défini dans stage_1a.py)
Ep = np.linspace(0, actor.energyrange, actor.energybins + 1)

# Bins d'énergie de la base de données : 0 → 200 MeV en 500 bins
E_db = np.linspace(0, 200, 500)

# Pour chaque bin de l'acteur, on cherche le bin le plus proche dans la base de données
# Résultat : tableau de 251 indices → "le bin i de l'acteur = bin ind_db[i] de la base"
ind_db = np.abs(E_db[None, :] - Ep[:, None]).argmin(axis=1)

# ── 1ère DIFFÉRENCE AVEC L'ORIGINAL ──────────────────────────────────────────
# Dans l'original, G[ind_db] était recalculé à CHAQUE VOXEL dans la triple boucle
# → recalculé ~3 millions de fois inutilement.
# Ici on le fait UNE SEULE FOIS pour tous les matériaux avant le calcul.
#
# G[ind_db] sélectionne les 250 lignes de G (500x250) correspondant aux énergies
# de l'acteur (0→130 MeV) → résultat (250x250), le reste ne nous concerne pas.
# ─────────────────────────────────────────────────────────────────────────────
Gamma_p_indexed = {name: G[ind_db] for name, G in gamma_proton.items()}   # (250, 250)
if neutr:
    Gamma_n_indexed = {name: G[ind_db] for name, G in gamma_neutron.items()}


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4 : CHARGER LES IMAGES 4D
# ═══════════════════════════════════════════════════════════════════════════════

# MJ = True  → fichier fusionné de plusieurs jobs : vpg_merged_prot_e.nii.gz
# MJ = False → fichier job unique                 : vpg_prot_e.nii.gz
suffix  = "merged_prot_e" if MJ else "prot_e"
img_p   = itk.imread(os.path.join(path.output, f"{File_name}_{suffix}.nii.gz"))


array_p = itk.array_from_image(img_p)   # shape (250, nZ, nY, nX)
print(array_p.shape)                                       # 250 = bins énergie, nZ/nY/nX = dimensions spatiales

if neutr:
    suffix_n = "merged_neutr_e" if MJ else "neutr_e"
    img_n    = itk.imread(os.path.join(path.output, f"{File_name}_{suffix_n}.nii.gz"))
    array_n  = itk.array_from_image(img_n)

# Récupération des dimensions du volume
nE, nZ, nY, nX = array_p.shape

# Initialisation des tableaux de sortie
treated_array_p = np.zeros((actor.energybins, nZ, nY, nX), dtype=np.float32)
treated_array_n = np.zeros((actor.energybins, nZ, nY, nX), dtype=np.float32)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5 : CONSTRUIRE LA CARTE MATERIAU VECTORISEE
#
# Pour chaque voxel de l'acteur on a besoin du nom du matériau.
# On construit d'abord un tableau UH_actor (nZ, nY, nX) contenant les valeurs
# HU dans le repère de l'acteur (éventuellement re-dimensionné depuis le CT).
#
# mat est une liste de triplets [HU_min, HU_max, nom] triée par HU croissant.
# Dans l'original, voxel_to_mat_name() parcourait cette liste voxel par voxel.
# Ici on reproduit cette logique de façon vectorisée avec np.searchsorted :
# on mappe HU → matériau pour TOUT LE VOLUME en une seule opération.
# ═══════════════════════════════════════════════════════════════════════════════

# Bornes supérieures et noms des matériaux dans l'ordre de la liste mat
HU_bounds = np.array([m[1] for m in mat], dtype=np.float32)
mat_names  = [m[2] for m in mat]

def build_UH_actor_array(UH_array, ct_object, actor, nZ, nY, nX, tle_simu):
    """
    Retourne un tableau (nZ, nY, nX) des valeurs HU dans le repère de l'acteur.
    Si CT et acteur ont la même taille → copie directe (vectorisé, rapide).
    Sinon → conversion de coordonnées avec redim() voxel par voxel (inévitable).
    """
    same_size = (ct_object.size == actor.size)
    UH_actor  = np.full((nZ, nY, nX), -3024, dtype=np.int16)

    if same_size:
        # Correspondance directe : les voxels sont alignés → slicing NumPy
        UH_actor[:nZ, :nY, :nX] = UH_array[:nZ, :nY, :nX]
    else:
        # Rééchantillonnage géométrique : obligatoire si les grilles diffèrent
        for z in range(nZ):
            for y in range(nY):
                for x in range(nX):
                    X, Y, Z = tle_simu.redim((x, y, z), ct_object, actor, None)
                    UH_actor[z, y, x] = UH_array[Z, Y, X]
    return UH_actor

#print("Construction de la carte HU dans le repère acteur...")
UH_actor = build_UH_actor_array(UH_array, ct_object, actor, nZ, nY, nX, tle_simu)

# Masque booléen des voxels valides : True = matière, False = vide (HU = -3024)
valid_mask = (UH_actor != -3024)   # shape (nZ, nY, nX)

# ── 2ème DIFFÉRENCE AVEC L'ORIGINAL ──────────────────────────────────────────
# Dans l'original, voxel_to_mat_name() était appelé à CHAQUE VOXEL dans la boucle.
# Ici on mappe HU → matériau pour TOUT LE VOLUME en une seule opération
# via np.searchsorted : retourne l'indice i tel que HU_bounds[i-1] < HU <= HU_bounds[i]
# ─────────────────────────────────────────────────────────────────────────────
hu_flat = UH_actor[valid_mask].astype(np.float32)
mat_idx = np.searchsorted(HU_bounds, hu_flat, side='left')
mat_idx = np.clip(mat_idx, 0, len(mat_names) - 1)   # sécurité sur les bords


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6 : TRAITEMENT VECTORISE PAR MATERIAU
#
# Pour chaque matériau unique présent dans le volume :
#   - extraire tous ses voxels en une seule opération
#   - calculer les spectres en un seul np.dot : (N_voxels × 250) @ (250 × 250)
#   - remettre en place dans le tableau de sortie
#
# Complexité : N_matériaux × (N_voxels_mat × 250²)  au lieu de  N_voxels × 250²
# Gain mémoire : on ne copie que les voxels non nuls pour chaque matériau.
# ═══════════════════════════════════════════════════════════════════════════════

#print("Calcul vectorisé des spectres gamma...")
start_wall = time.time()

# Aplatir array_p sur les dimensions spatiales pour faciliter le masquage
# (250, nZ, nY, nX) → (250, N_valid) où N_valid = nombre de voxels non vides
array_p_flat = array_p[:, valid_mask]
if neutr:
    array_n_flat = array_n[:, valid_mask]

# Tableaux de sortie aplatis
out_p_flat = np.zeros((actor.energybins, valid_mask.sum()), dtype=np.float32)
out_n_flat = np.zeros((actor.energybins, valid_mask.sum()), dtype=np.float32)

# Indices des matériaux uniques présents dans le volume
unique_mat_idx = np.unique(mat_idx)

for idx in unique_mat_idx:
    name = mat_names[idx]
    sel  = (mat_idx == idx)   # masque des voxels de ce matériau parmi les valides

    # ── Protons ──────────────────────────────────────────────────────────────
    Gp    = Gamma_p_indexed[name].astype(np.float32)   # (250, 250)
    vox_p = array_p_flat[:, sel].T                      # (N_voxels_mat, 250)

    # On ne calcule que pour les voxels qui ont reçu des particules
    nonzero = vox_p.sum(axis=1) != 0
    if nonzero.any():
        spectra_p = np.dot(vox_p[nonzero], Gp)            # (N_nz, 250)
        tmp = np.zeros((sel.sum(), actor.energybins), dtype=np.float32)
        tmp[nonzero] = spectra_p
        out_p_flat[:, sel] = tmp.T

    # ── Neutrons ─────────────────────────────────────────────────────────────
    if neutr:
        Gn    = Gamma_n_indexed[name].astype(np.float32)
        vox_n = array_n_flat[:, sel].T

        nonzero_n = vox_n.sum(axis=1) != 0
        if nonzero_n.any():
            spectra_n = np.dot(vox_n[nonzero_n], Gn)
            tmp_n = np.zeros((sel.sum(), actor.energybins), dtype=np.float32)
            tmp_n[nonzero_n] = spectra_n
            out_n_flat[:, sel] = tmp_n.T

wall_time = time.time() - start_wall
print(f"Wall-clock triple vectorisation : {wall_time:.2f} secondes")
print(f"Wall-clock total time : {time.time() - start_total_wall:.2f} secondes")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7 : REMETTRE EN FORME 4D ET SAUVEGARDER
# ═══════════════════════════════════════════════════════════════════════════════

# Remettre les résultats dans les tableaux 4D (250, nZ, nY, nX)
treated_array_p[:, valid_mask] = out_p_flat
if neutr:
    treated_array_n[:, valid_mask] = out_n_flat

# Spectre gamma total = contribution protons + contribution neutrons
gamma_array = treated_array_p + treated_array_n

def save_nii(arr, ref_img, filepath):
    """Sauvegarde un tableau NumPy en image NIfTI en copiant les métadonnées."""
    out = itk.image_from_array(arr)
    out.CopyInformation(ref_img)
    itk.imwrite(out, str(filepath))
   # print(f"  Sauvegardé : {filepath}")

# Préfixe pour les noms de fichiers
prefix = "merged" if MJ else ""
sep    = "_" if prefix else ""


#print("Sauvegarde des résultats...")
save_nii(treated_array_p, img_p,
         path.output / f"VPG_{File_name}{sep}{prefix}_jeanpaul_prot_e.nii.gz")

if neutr:
    save_nii(treated_array_n, img_n,
             path.output / f"VPG_{File_name}{sep}{prefix}_jeanpaul_neutr_e.nii.gz")

save_nii(gamma_array, img_p,
         path.output / f"VPG_{File_name}{sep}{prefix}_jeanpaul_gamma_e.nii.gz")

print("Terminé.")
