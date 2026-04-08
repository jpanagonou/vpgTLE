#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#=======================================================================
# Comparaison Original vs Optimisé
# Cartes 2D : plateau et pic de Bragg
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
y_axis = np.arange(len(profil))


plt.figure(figsize=(8, 3))
plt.plot(y_axis, profil, linewidth=2, color='blue')
plt.xlabel("Index Y")
plt.ylabel("Emission gamma totale")
plt.title("Profil longitudinal\nIdentification plateau et pic de Bragg")
plt.gca().yaxis.set_major_formatter(sci_formatter())
plt.tight_layout()
plt.savefig("profil_longitudinal.png", dpi=150, bbox_inches='tight')
print("Sauvegardé : profil_longitudinal.png")
plt.close()


# ── Identification automatique du plateau et du pic ────────────────────
y_bragg   = int(np.argmax(profil))
y_plateau = int(y_bragg // 2)
print(f"\nIndex Y plateau   : {y_plateau}")
print(f"Index Y pic Bragg : {y_bragg}")


# ── Cartes 2D : sommer sur nE → shape (nZ, nY, nX) ───────────────────
carte_orig = arr_orig.sum(axis=0)
carte_opti = arr_opti.sum(axis=0)


# ── Fonction pour tracer une carte à un indice Y donné ────────────────
def plot_carte(y_idx, label, filename_prefix):
    coupe_orig = carte_orig[:, y_idx, :]   # (nZ, nX)
    coupe_opti = carte_opti[:, y_idx, :]   # (nZ, nX)
    ecart      = np.abs(coupe_orig - coupe_opti) / (coupe_orig + 1e-30) * 100


    vmax = max(coupe_orig.max(), coupe_opti.max())


    # ── Figure 1 : comparaison original vs optimisé ───────────────────
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
    fig.subplots_adjust(wspace=0.1, top=0.82)


    im1 = axes[0].imshow(coupe_orig, cmap='Reds', origin='lower', vmin=0, vmax=vmax)
    axes[0].set_title(f'Original')
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Z")
    cb1 = plt.colorbar(im1, ax=axes[0], label='PG yields', pad=0.02)
    cb1.ax.yaxis.set_major_formatter(sci_formatter())


    im2 = axes[1].imshow(coupe_opti, cmap='Reds', origin='lower', vmin=0, vmax=vmax)
    axes[1].set_title(f'Optimised')
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Z")
    cb2 = plt.colorbar(im2, ax=axes[1], label='PG yields', pad=0.02)
    cb2.ax.yaxis.set_major_formatter(sci_formatter())


    fig.suptitle(f"Carte d'émission gamma 2D — {label} (Y={y_idx})", y=0.98)
    plt.savefig(f'{filename_prefix}_comparaison.png', dpi=150, bbox_inches='tight')
    print(f"Sauvegardé : {filename_prefix}_comparaison.png")
    plt.close()


    # ── Figure 2 : écart relatif ──────────────────────────────────────
    fig, ax = plt.subplots(figsize=(4, 3))
    im = ax.imshow(ecart, cmap='Reds', origin='lower', vmin=0, vmax=ecart.max())
    ax.set_xlabel("X")
    ax.set_ylabel("Z")
    cb = plt.colorbar(im, ax=ax, label='Écart relatif (%)')
    cb.ax.yaxis.set_major_formatter(sci_formatter())
    plt.title(f"Écart relatif - {label} \n(Y={y_idx})")
    plt.tight_layout()
    plt.savefig(f'{filename_prefix}_ecart.png', dpi=150, bbox_inches='tight')
    print(f"Sauvegardé : {filename_prefix}_ecart.png")
    plt.close()


    print(f"  Écart relatif moyen : {ecart.mean():.2e} %")
    print(f"  Écart relatif max   : {ecart.max():.2e} %")




# ── Cartes plateau et pic de Bragg ────────────────────────────────────
print("\n── Plateau ──")
plot_carte(y_plateau, "Plateau", "carte_plateau")


print("\n── Pic de Bragg ──")
plot_carte(y_bragg, "Pic de Bragg", "carte_bragg")




# ── Fonction pour tracer un spectre par voxel ─────────────────────────
def plot_spectre(z_idx, y_idx, x_idx, label, filename):
    spectre_orig = arr_orig[:, z_idx, y_idx, x_idx]
    spectre_opti = arr_opti[:, z_idx, y_idx, x_idx]
    E = np.linspace(0, 10, len(spectre_orig))


    ecart_rel = (spectre_orig - spectre_opti) / (spectre_orig + 1e-30) * 100


    fig1, ax1 = plt.subplots(figsize= (7, 4))
    fig2, ax2 = plt.subplots(figsize= (7, 4))
    axes = [ax1, ax2]

    # Spectres superposés
    axes[0].plot(E, spectre_orig, linewidth=1.5,   label='Original', color='blue')
    axes[0].plot(E, spectre_opti, linewidth=1, label='Optimisé', color='red', linestyle='--')
    axes[0].set_xlabel("Énergie des PG (MeV)")
    axes[0].set_ylabel("PG yields / primary")
    axes[0].set_title(f'Spectres des PG dans un voxel du {label} : Original/Optimisé)')
    axes[0].legend()
    axes[0].yaxis.set_major_formatter(sci_formatter())


    # Écart relatif
    axes[1].plot(E, ecart_rel, linewidth=1, color='green')
    axes[1].set_xlabel("Énergie des PG (MeV)")
    axes[1].set_ylabel("Écart relatif (%)")
    axes[1].set_title(f'Écart relatif  {label}')
    axes[1].yaxis.set_major_formatter(sci_formatter())



    fig1.tight_layout()
    fig1.savefig(filename.replace('.png', '_spectre.png'), dpi=150, bbox_inches='tight')

    fig2.tight_layout()
    fig2.savefig(filename.replace('.png', '_ecart.png'), dpi=150, bbox_inches='tight')


    print(f"Sauvegardé : {filename}")
    plt.close()


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





