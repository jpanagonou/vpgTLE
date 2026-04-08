#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#=======================================================================
# Comparaison Violette vs JP
# Spectres par voxel : plateau et pic de Bragg
# Faisceau dirigé suivant Y
#=======================================================================

import itk
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# ── Formatter scientifique ─────────────────────────────────────────────
def sci_formatter():
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_scientific(True)
    fmt.set_powerlimits((0, 0))
    return fmt

# ── Chargement des images ──────────────────────────────────────────────
arr_orig = itk.array_from_image(itk.imread('output/stage1a/VPG_vpg_violette_gamma_e.nii.gz'))
arr_opti = itk.array_from_image(itk.imread('output/stage1a/VPG_vpg_jeanpaul_gamma_e.nii.gz'))

print(f"Shape original : {arr_orig.shape}")
print(f"Shape optimisé : {arr_opti.shape}")
# shape : (nE, nZ, nY, nX)

# ── Profil longitudinal le long de Y (direction du faisceau) ──────────
profil = arr_orig.sum(axis=(0, 1, 3))  # (nY,)
y_bragg   = int(np.argmax(profil))
y_plateau = int(y_bragg // 2)
print(f"\nIndex Y plateau   : {y_plateau}")
print(f"Index Y pic Bragg : {y_bragg}")

# ── Cartes pour trouver les voxels actifs ─────────────────────────────
carte_orig = arr_orig.sum(axis=0)

# ── Fonction pour tracer un spectre par voxel ─────────────────────────
def plot_spectre(z_idx, y_idx, x_idx, label, filename):
    spectre_orig = arr_orig[:, z_idx, y_idx, x_idx]
    spectre_opti = arr_opti[:, z_idx, y_idx, x_idx]
    E = np.linspace(0, 10, len(spectre_orig))

    ecart_rel = (spectre_orig - spectre_opti) / (spectre_orig + 1e-30) 

    # ── Figure 1 : spectres superposés ───────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(7, 4))
    ax1.plot(E, spectre_orig, linewidth=1.5, label='Violette', color='blue')
    ax1.plot(E, spectre_opti, linewidth=1,   label='JP',       color='red', linestyle='--')
    ax1.set_xlabel("Énergie des PG (MeV)")
    ax1.set_ylabel("PG yields / primary")
    ax1.set_title(f'Spectres des PG — {label}\n(Z={z_idx}, Y={y_idx}, X={x_idx})')
    ax1.legend()
    ax1.yaxis.set_major_formatter(sci_formatter())
    fig1.tight_layout()
    fig1.savefig(filename.replace('.png', '_spectre.png'), dpi=150, bbox_inches='tight')
    print(f"Sauvegardé : {filename.replace('.png', '_spectre.png')}")
    plt.close(fig1)

    # ── Figure 2 : écart relatif ──────────────────────────────────────
    fig2, ax2 = plt.subplots(figsize=(7, 4))
    ax2.plot(E, ecart_rel, linewidth=1, color='green')
    ax2.set_xlabel("Énergie des PG (MeV)")
    ax2.set_ylabel("Écart relatif (%)")
    ax2.set_title(f'Écart relatif Violette vs JP — {label}')
    ax2.yaxis.set_major_formatter(sci_formatter())
    fig2.tight_layout()
    fig2.savefig(filename.replace('.png', '_ecart.png'), dpi=150, bbox_inches='tight')
    print(f"Sauvegardé : {filename.replace('.png', '_ecart.png')}")
    plt.close(fig2)

    print(f"  Écart relatif moyen : {ecart_rel.mean():.2e} %")
    print(f"  Écart relatif max   : {ecart_rel.max():.2e} %")


# ── Trouver un voxel actif dans le plateau ────────────────────────────
coupe_plateau  = carte_orig[:, y_plateau, :]
voxels_plateau = np.argwhere(coupe_plateau > 0)
if len(voxels_plateau) > 0:
    z_p, x_p = voxels_plateau[len(voxels_plateau) // 2]
    print(f"\n── Spectre plateau — voxel (Z={z_p}, Y={y_plateau}, X={x_p}) ──")
    plot_spectre(z_p, y_plateau, x_p, "Plateau", "spectre_plateau.png")
else:
    print("Aucun voxel actif dans le plateau")

# ── Trouver un voxel actif dans le pic de Bragg ───────────────────────
coupe_bragg  = carte_orig[:, y_bragg, :]
voxels_bragg = np.argwhere(coupe_bragg > 0)
if len(voxels_bragg) > 0:
    z_b, x_b = voxels_bragg[len(voxels_bragg) // 2]
    print(f"\n── Spectre pic de Bragg — voxel (Z={z_b}, Y={y_bragg}, X={x_b}) ──")
    plot_spectre(z_b, y_bragg, x_b, "Pic de Bragg", "spectre_bragg.png")
else:
    print("Aucun voxel actif dans le pic de Bragg")

print("\nTerminé.")

