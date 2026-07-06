#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import itk
import opengate as gate
import os
import time
import numpy as np
from variablesoptimized import tle_sim
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

ct_object  = tle_simu.ct
UH_array   = itk.array_from_image(ct_object.load_input_image())

path       = tle_simu.path
actor      = tle_simu.act
actor_list = tle_simu.actor_list

data_protons  = tle_simu.root_file
data_neutrons = tle_simu.root_file_neutron

mat = tle_simu.voxel_mat
print(f"Nombre total de matériaux dans mat : {len(mat)}")

neutr = "neutr_e.nii.gz" in actor_list

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 : INDEX DE CORRESPONDANCE BASE DE DONNEES → BINS ACTEUR
# ═══════════════════════════════════════════════════════════════════════════════

Ep = np.linspace(0, actor.energyrange, actor.energybins + 1)
E_db = np.linspace(0, 200, 500)
ind_db = np.abs(E_db[None, :] - Ep[:, None]).argmin(axis=1)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 : CHARGER LES IMAGES 4D
# ═══════════════════════════════════════════════════════════════════════════════

suffix  = "merged_prot_e" if MJ else "prot_e"
img_p   = itk.imread(os.path.join(path.output, f"{File_name}_{suffix}.nii.gz"))

array_p = itk.array_from_image(img_p)
print(array_p.shape)

if neutr:
    suffix_n = "merged_neutr_e" if MJ else "neutr_e"
    img_n    = itk.imread(os.path.join(path.output, f"{File_name}_{suffix_n}.nii.gz"))
    array_n  = itk.array_from_image(img_n)

nE, nZ, nY, nX = array_p.shape

treated_array_p = np.zeros((actor.energybins, nZ, nY, nX), dtype=np.float32)
treated_array_n = np.zeros((actor.energybins, nZ, nY, nX), dtype=np.float32)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4 : CONSTRUIRE LA CARTE MATERIAU VECTORISEE
# ═══════════════════════════════════════════════════════════════════════════════

HU_bounds = np.array([m[1] for m in mat], dtype=np.float32)
mat_names  = [m[2] for m in mat]

UH_actor = UH_array[:nZ, :nY, :nX].copy()

valid_mask = (UH_actor != -3024)

hu_flat = UH_actor[valid_mask].astype(np.float32)
mat_idx = np.searchsorted(HU_bounds, hu_flat, side='left')
mat_idx = np.clip(mat_idx, 0, len(mat_names) - 1)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5 : CALCULER GAMMA_MAT UNIQUEMENT POUR LES MATERIAUX PRESENTS
# ═══════════════════════════════════════════════════════════════════════════════

unique_mat_idx = np.unique(mat_idx)
print(f"Nombre de materiaux uniques dans le CT : {len(unique_mat_idx)}")
present_mat_names = set(mat_names[i] for i in unique_mat_idx)

gamma_proton  = {}
gamma_neutron = {}
for name in present_mat_names:
    gamma_proton[name] = tle_simu.gamma_mat(name, data_protons)
    if neutr:
        gamma_neutron[name] = tle_simu.gamma_mat(name, data_neutrons)

Gamma_p_indexed = {name: G[ind_db] for name, G in gamma_proton.items()}
if neutr:
    Gamma_n_indexed = {name: G[ind_db] for name, G in gamma_neutron.items()}

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6 : TRAITEMENT VECTORISE PAR MATERIAU
# ═══════════════════════════════════════════════════════════════════════════════

start_wall = time.time()

array_p_flat = array_p[:, valid_mask]
if neutr:
    array_n_flat = array_n[:, valid_mask]

out_p_flat = np.zeros((actor.energybins, valid_mask.sum()), dtype=np.float32)
out_n_flat = np.zeros((actor.energybins, valid_mask.sum()), dtype=np.float32)

for idx in unique_mat_idx:
    name = mat_names[idx]
    sel  = (mat_idx == idx)

    Gp    = Gamma_p_indexed[name].astype(np.float32)
    vox_p = array_p_flat[:, sel].T

    nonzero = vox_p.sum(axis=1) != 0
    if nonzero.any():
        spectra_p = np.dot(vox_p[nonzero], Gp)
        tmp = np.zeros((sel.sum(), actor.energybins), dtype=np.float32)
        tmp[nonzero] = spectra_p
        out_p_flat[:, sel] = tmp.T

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

treated_array_p[:, valid_mask] = out_p_flat
if neutr:
    treated_array_n[:, valid_mask] = out_n_flat

gamma_array = treated_array_p + treated_array_n

def save_nii(arr, ref_img, filepath):
    out = itk.image_from_array(arr)
    out.CopyInformation(ref_img)
    itk.imwrite(out, str(filepath))

prefix = "merged" if MJ else ""
sep    = "_" if prefix else ""

save_nii(treated_array_p, img_p,
         path.output / f"VPG_{File_name}{sep}{prefix}_jeanpaul_prot_e.nii.gz")

if neutr:
    save_nii(treated_array_n, img_n,
             path.output / f"VPG_{File_name}{sep}{prefix}_jeanpaul_neutr_e.nii.gz")

save_nii(gamma_array, img_p,
         path.output / f"VPG_{File_name}{sep}{prefix}_jeanpaul_gamma_e.nii.gz")

print("Terminé.")
