#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#=======================================================================
# Comparaison 3 versions :
# - Analog
# - TLE original (Violette)
# - TLE optimisé (JP)
# Spectres par voxel : plateau et pic de Bragg
# Bande d'incertitude statistique de Poisson à 3 sigma
# Faisceau dirigé suivant Y
#=======================================================================

import itk
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# ── Nombre de particules analog (pour l'incertitude de Poisson) ────────
N_analog = 1e6

# ── Formatter scientifique ─────────────────────────────────────────────
def sci_formatter():
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_scientific(True)
    fmt.set_powerlimits((0, 0))
    return fmt

# ── Chargement des images ──────────────────────────────────────────────
arr_analog   = itk.array_from_image(itk.imread('stage1_Analog/output/stage1a/vpg_prot_e.nii.gz'))
arr_violette = itk.array_from_image(itk.imread('stage1_vpg/output/stage1a/VPG_vpg_violette_prot_e.nii.gz'))
arr_jp       = itk.array_from_image(itk.imread('stage1_vpg/output/stage1a/VPG_vpg_jeanpaul_prot_e.nii.gz'))

print(f"Shape Analog   : {arr_analog.shape}")
print(f"Shape Violette : {arr_violette.shape}")
print(f"Shape JP       : {arr_jp.shape}")

# ── Tronquer le dernier bin de l'analog (overflow) ────────────────────
if arr_analog.shape[0] == 251:
    arr_analog = arr_analog[:-1, :, :, :]
    print(f"Shape Analog après troncature : {arr_analog.shape}")

# ── Profil longitudinal le long de Y (direction du faisceau) ──────────
profil    = arr_analog.sum(axis=(0, 1, 3))
y_bragg   = int(np.argmax(profil))
y_plateau = int(y_bragg // 2)
print(f"\nIndex Y plateau   : {y_plateau}")
print(f"Index Y pic Bragg : {y_bragg}")

# ── Cartes pour trouver les voxels actifs ─────────────────────────────
carte_analog = arr_analog.sum(axis=0)

# ── Fonction pour tracer un spectre par voxel ─────────────────────────
def plot_spectre(z_idx, y_idx, x_idx, label, filename):
    spectre_analog   = arr_analog[:, z_idx, y_idx, x_idx]
    spectre_violette = arr_violette[:, z_idx, y_idx, x_idx]
    spectre_jp       = arr_jp[:, z_idx, y_idx, x_idx]
    E = np.linspace(0, 10, len(spectre_analog))

    # ── Incertitude de Poisson à 3 sigma ─────────────────────────────
    inc = np.sqrt(spectre_analog / N_analog)

    # ── Masque : bins où l'analog est significatif (> 1% du max) ─────
    seuil  = spectre_analog.max() * 0.01
    masque = spectre_analog > seuil

    ecart_viol = np.zeros_like(spectre_analog, dtype=float)
    ecart_jp   = np.zeros_like(spectre_analog, dtype=float)
    ecart_viol[masque] = np.abs(spectre_analog[masque] - spectre_violette[masque]) / spectre_analog[masque] * 100
    ecart_jp[masque]   = np.abs(spectre_analog[masque] - spectre_jp[masque])       / spectre_analog[masque] * 100

    # ── Figure 1 : spectres superposés avec bande d'incertitude ──────
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.fill_between(E,
                    spectre_analog + 3 * inc,
                    spectre_analog - 3 * inc,
                    color='blue', alpha=0.2, label='Incertitude Analog (3σ)')
    ax.plot(E, spectre_analog,   linewidth=1.5, label='Analog',       color='blue')
    ax.plot(E, spectre_violette, linewidth=1,   label='TLE Violette', color='green', linestyle='--')
    ax.plot(E, spectre_jp,       linewidth=1,   label='TLE JP',       color='red',   linestyle=':')
    ax.set_xlabel("Énergie des PG (MeV)")
    ax.set_ylabel("PG yields / primary")
    ax.set_title(f'Spectres des PG — {label}\n(Z={z_idx}, Y={y_idx}, X={x_idx})')
    ax.legend()
    ax.grid(True)
    ax.yaxis.set_major_formatter(sci_formatter())
    fig.tight_layout()
    fig.savefig(filename.replace('.png', '_spectre.png'), dpi=150, bbox_inches='tight')
    print(f"Sauvegardé : {filename.replace('.png', '_spectre.png')}")
    plt.close(fig)

    # ── Figure 2 : écarts relatifs ────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(E, ecart_viol, linewidth=1, label='Analog vs Violette', color='green')
    ax.plot(E, ecart_jp,   linewidth=1, label='Analog vs JP',       color='red', linestyle='--')
    ax.set_xlabel("Énergie des PG (MeV)")
    ax.set_ylabel("Écart relatif (%)")
    ax.set_title(f'Écarts relatifs — {label}')
    ax.legend()
    ax.grid(True)
    ax.yaxis.set_major_formatter(sci_formatter())
    fig.tight_layout()
    fig.savefig(filename.replace('.png', '_ecart.png'), dpi=150, bbox_inches='tight')
    print(f"Sauvegardé : {filename.replace('.png', '_ecart.png')}")
    plt.close(fig)

    print(f"  Écart Analog vs Violette — moyen : {ecart_viol[masque].mean():.2e} %  max : {ecart_viol[masque].max():.2e} %")
    print(f"  Écart Analog vs JP       — moyen : {ecart_jp[masque].mean():.2e} %  max : {ecart_jp[masque].max():.2e} %")


# ── Trouver un voxel actif dans le plateau ────────────────────────────
coupe_plateau  = carte_analog[:, y_plateau, :]
voxels_plateau = np.argwhere(coupe_plateau > 0)
if len(voxels_plateau) > 0:
    z_p, x_p = voxels_plateau[len(voxels_plateau) // 2]
    print(f"\n── Spectre plateau — voxel (Z={z_p}, Y={y_plateau}, X={x_p}) ──")
    plot_spectre(z_p, y_plateau, x_p, "Plateau", "spectre_plateau.png")
else:
    print("Aucun voxel actif dans le plateau")

# ── Trouver un voxel actif dans le pic de Bragg ───────────────────────
coupe_bragg  = carte_analog[:, y_bragg, :]
voxels_bragg = np.argwhere(coupe_bragg > 0)
if len(voxels_bragg) > 0:
    z_b, x_b = voxels_bragg[len(voxels_bragg) // 2]
    print(f"\n── Spectre pic de Bragg — voxel (Z={z_b}, Y={y_bragg}, X={x_b}) ──")
    plot_spectre(z_b, y_bragg, x_b, "Pic de Bragg", "spectre_bragg.png")
else:
    print("Aucun voxel actif dans le pic de Bragg")

print("\nTerminé.")

