#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#=======================================================================
# stage_1b_oxygen_violette.py
# Version originale de Violette (triple boucle) adaptée pour la boite Oxygen
# Pour valider Gate 10 vs Gate 9
#=======================================================================

import itk
import uproot
import hist
import opengate as gate
import os
import time
import numpy as np
from variables_oxygen import tle_sim
from stage_1a_oxygen import *

# ── Initialisation ─────────────────────────────────────────────────────
sim = simulation(
    paths, File_name, 0, number_of_particles, False, False, actor, Erange
)

tle_simu = tle_sim(sim)

# ── FIRST STEP : récupérer les objets utiles ───────────────────────────

# Pour la boite Oxygen, pas de CT → UH_array = tableau de zéros (1 seul voxel)
ct_object = tle_simu.ct
UH_array  = np.zeros((1, 1, 1), dtype=np.int16)

# Chemins, acteur et liste des sorties
path       = tle_simu.path
actor      = tle_simu.act
actor_list = tle_simu.actor_list

# Bases de données
data_protons  = tle_simu.root_file
data_neutrons = tle_simu.root_file_neutron
mat           = tle_simu.voxel_mat

# ── SECOND STEP : calculer les matrices Gamma par matériau ─────────────
neutr = False
if "neutr_e.nii.gz" in actor_list:
    neutr = True

gamma_neutron = {}
gamma_proton  = {}

for data in mat:
    name = data[2]
    if neutr:
        gamma_neutron[name] = tle_simu.gamma_mat(name, data_neutrons)
    gamma_proton[name] = tle_simu.gamma_mat(name, data_protons)

# ── THIRD STEP : conversion histogrammes énergie → spectres gamma ──────
Ep     = np.linspace(0, actor.energyrange, actor.energybins + 1)
# db-Oxygen_merged.root a GammaZ de shape (250, 250) → E_db sur 250 bins
E_db   = np.linspace(0, 200, 250)
ind_db = np.abs(E_db[None, :] - Ep[:, None]).argmin(axis=1)

# Chargement des images 4D
if MJ:
    img_p = itk.imread(os.path.join(path.output, f"{File_name}_merged_prot_e.nii.gz"))
else:
    img_p = itk.imread(os.path.join(path.output, f"{File_name}_prot_e.nii.gz"))

array_p = itk.array_from_image(img_p)
array_p = array_p.reshape(array_p.shape[2], 1, array_p.shape[0], array_p.shape[1])
print(f"Shape array_p : {array_p.shape}")

# Initialisation des tableaux de sortie (250 bins de gamma)
treated_array_p = np.zeros((250, array_p.shape[1], array_p.shape[2], array_p.shape[3]))
treated_array_n = np.zeros((250, array_p.shape[1], array_p.shape[2], array_p.shape[3]))
gamma_array     = np.zeros((250, array_p.shape[1], array_p.shape[2], array_p.shape[3]))

if neutr:
    if MJ:
        img_n = itk.imread(os.path.join(path.output, f"{File_name}_merged_neutr_e.nii.gz"))
    else:
        img_n = itk.imread(os.path.join(path.output, f"{File_name}_neutr_e.nii.gz"))
    array_n = itk.array_from_image(img_n)

# ── Chronomètre ────────────────────────────────────────────────────────
start      = os.times()
start_wall = time.time()

# ── Triple boucle voxel par voxel ──────────────────────────────────────
for x in range(array_p.shape[3]):
    for y in range(array_p.shape[2]):
        for z in range(array_p.shape[1]):

            # Pour la boite Oxygen, tous les voxels sont G4_O
            # Pas de condition HU == -3024 (pas de vide dans la boite)
            name = tle_simu.voxel_to_mat_name(0, tle_simu.voxel_mat)

            # Histogramme énergie du voxel courant
            histo_E_p = array_p[:, z, y, x]
            spectrum_p = np.zeros(250)

            if histo_E_p.sum() != 0:
                Gamma_m_p  = gamma_proton[name]
                Gamma_p    = Gamma_m_p[ind_db]
                spectrum_p = np.dot(histo_E_p, Gamma_p)
                treated_array_p[:, z, y, x] = spectrum_p

            if neutr:
                histo_E_n  = array_n[:, z, y, x]
                spectrum_n = np.zeros(250)
                if histo_E_n.sum() != 0:
                    Gamma_m_n  = gamma_neutron[name]
                    Gamma_n    = Gamma_m_n[ind_db]
                    spectrum_n = np.dot(histo_E_n, Gamma_n)
                    treated_array_n[:, z, y, x] = spectrum_n

# ── Temps de calcul ────────────────────────────────────────────────────
end       = os.times()
end_wall  = time.time()

user_time   = end.user   - start.user
system_time = end.system - start.system
total_cpu   = user_time  + system_time
wall_time   = end_wall   - start_wall

print(f"User CPU time used:   {user_time:.2f} seconds")
print(f"System CPU time used: {system_time:.2f} seconds")
print(f"Total CPU time used:  {total_cpu:.2f} seconds")
print(f"Wall-clock time:      {wall_time:.2f} seconds")

# ── Sauvegarde ─────────────────────────────────────────────────────────

def save_nii(arr, filepath):
    out = itk.image_from_array(arr.astype(np.float32))
    itk.imwrite(out, str(filepath))
    print(f"  Sauvegardé : {filepath}")

gamma_array += treated_array_p + treated_array_n

if MJ:
    save_nii(treated_array_p, path.output / f"VPG_{File_name}_merged_prot_e.nii.gz")
else:
    save_nii(treated_array_p, path.output / f"VPG_{File_name}_prot_e.nii.gz")

if neutr:
    if MJ:
        save_nii(treated_array_n, path.output / f"VPG_{File_name}_merged_neutr_e.nii.gz")
    else:
        save_nii(treated_array_n, path.output / f"VPG_{File_name}_neutr_e.nii.gz")

if MJ:
    save_nii(gamma_array, path.output / f"VPG_{File_name}_merged_gamma_e.nii.gz")
else:
    save_nii(gamma_array, path.output / f"VPG_{File_name}_gamma_e.nii.gz")

print("Terminé.")


