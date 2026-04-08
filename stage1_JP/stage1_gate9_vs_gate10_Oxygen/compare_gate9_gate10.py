#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#=======================================================================
# Comparaison Gate 9 vs Gate 10 - Boite Oxygen
#=======================================================================

import SimpleITK as sitk
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Lecture Gate 10 ────────────────────────────────────────────────────
img10  = sitk.ReadImage("gate10/output/stage1a_oxygen/VPG_vpg_oxygen_gamma_e.nii.gz")
arr10  = sitk.GetArrayFromImage(img10)
print("Shape Gate 10:", arr10.shape)
spectre10 = arr10.flatten()

# ── Lecture Gate 9 ─────────────────────────────────────────────────────
img9   = sitk.ReadImage("gate9/output/source_vpg_Oxygen.mhd")
arr9   = sitk.GetArrayFromImage(img9)
print("Shape Gate 9 :", arr9.shape)
spectre9 = arr9.flatten()

# ── Axe énergie ────────────────────────────────────────────────────────
# Les deux ont 250 bins de 0 à 10 MeV
E = np.linspace(0, 10, len(spectre10))

# ── Comparaison numérique ──────────────────────────────────────────────
print(f"\nMax Gate 9  : {spectre9.max():.6e}")
print(f"Max Gate 10 : {spectre10.max():.6e}")
print(f"Rapport Gate9/Gate10 : {spectre9.max() / spectre10.max():.4f}")

# Ecart relatif point par point
ecart_relatif = (spectre9 - spectre10) / (spectre9 + 1e-30) * 100
print(f"\nEcart relatif moyen  : {ecart_relatif.mean():.2f} %")
print(f"Ecart relatif max    : {ecart_relatif.max():.2f} %")

# ── Figure 1 : Superposition des deux spectres ─────────────────────────
plt.figure(figsize=(6, 3))
plt.plot(E, spectre9,  linewidth=1.5,   label='Gate 9',  color='blue')
plt.plot(E, spectre10, linewidth=1, label='Gate 10', color='red', linestyle='--')
plt.xlabel("Energie des PG (MeV)")
plt.ylabel("PG yields / primary / 40 keV")
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.title("PG yield de Gate 9 vs Gate 10 du stage 1 violette - Oxygène")
plt.legend()
plt.tight_layout()
plt.savefig("spectres_gate9_gate10.png", dpi=300)
print("Figure 1 sauvegardée : spectres_gate9_gate10.png")

# ── Figure 2 : Ecart relatif ───────────────────────────────────────────
plt.figure(figsize=(6, 3))
plt.plot(E, ecart_relatif, linewidth=1.5, color='green')
plt.xlabel("Energie des PG (MeV)")
plt.ylabel("Ecart relatif (%)")
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.title("Ecart relatif de PG yield de Gate 9 vs Gate 10 - Oxygène")

plt.tight_layout()
plt.savefig("ecart_gate9_gate10.png", dpi=300)
print("Figure 2 sauvegardée : ecart_gate9_gate10.png")

print("\nFigure sauvegardée : comparaison_gate9_gate10.png")
