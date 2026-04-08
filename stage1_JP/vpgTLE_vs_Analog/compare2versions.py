import itk
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

def sci_formatter():
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_scientific(True)
    fmt.set_powerlimits((0, 0))
    return fmt

# ── Chargement ─────────────────────────────────────────────────────────
arr_analog = itk.array_from_image(itk.imread('stage1_Analog/output/stage1a/vpg_prot_e.nii.gz'))
arr_tle    = itk.array_from_image(itk.imread('stage1_vpg/output/stage1a/VPG_vpg_jeanpaul_prot_e.nii.gz'))

# Tronquer overflow
if arr_analog.shape[0] == 251:
    arr_analog = arr_analog[:-1, :, :, :]

# ── Profil longitudinal ────────────────────────────────────────────────
profil = arr_analog.sum(axis=(0, 1, 3))
y_bragg   = int(np.argmax(profil))
y_plateau = int(y_bragg // 2)

carte_analog = arr_analog.sum(axis=0)
carte_tle    = arr_tle.sum(axis=0)

# ── Voxels actifs ──────────────────────────────────────────────────────
coupe_plateau = carte_analog[:, y_plateau, :]
voxels_plateau = np.argwhere(coupe_plateau > 0)
z_p, x_p = voxels_plateau[len(voxels_plateau) // 2]

coupe_bragg = carte_analog[:, y_bragg, :]
voxels_bragg = np.argwhere(coupe_bragg > 0)
z_b, x_b = voxels_bragg[len(voxels_bragg) // 2]

E = np.linspace(0, 10, arr_analog.shape[0])

# ── Figure Analog plateau ──────────────────────────────────────────────
plt.figure(figsize=(7, 4))
plt.plot(E, arr_analog[:, z_p, y_plateau, x_p], linewidth=1.5, color='blue')
plt.xlabel("Énergie des PG (MeV)")
plt.ylabel("PG yields / primary")
plt.title(f'Analog — Plateau (Z={z_p}, Y={y_plateau}, X={x_p})')
plt.grid(True)
plt.gca().yaxis.set_major_formatter(sci_formatter())
plt.tight_layout()
plt.savefig("analog_plateau.png", dpi=150, bbox_inches='tight')
plt.close()

# ── Figure TLE plateau ─────────────────────────────────────────────────
plt.figure(figsize=(7, 4))
plt.plot(E, arr_tle[:, z_p, y_plateau, x_p], linewidth=1.5, color='red')
plt.xlabel("Énergie des PG (MeV)")
plt.ylabel("PG yields / primary")
plt.title(f'TLE optimisé — Plateau (Z={z_p}, Y={y_plateau}, X={x_p})')
plt.grid(True)
plt.gca().yaxis.set_major_formatter(sci_formatter())
plt.tight_layout()
plt.savefig("tle_plateau.png", dpi=150, bbox_inches='tight')
plt.close()

# ── Figure Analog pic de Bragg ─────────────────────────────────────────
plt.figure(figsize=(7, 4))
plt.plot(E, arr_analog[:, z_b, y_bragg, x_b], linewidth=1.5, color='blue')
plt.xlabel("Énergie des PG (MeV)")
plt.ylabel("PG yields / primary")
plt.title(f'Analog — Pic de Bragg (Z={z_b}, Y={y_bragg}, X={x_b})')
plt.grid(True)
plt.gca().yaxis.set_major_formatter(sci_formatter())
plt.tight_layout()
plt.savefig("analog_bragg.png", dpi=150, bbox_inches='tight')
plt.close()

# ── Figure TLE pic de Bragg ────────────────────────────────────────────
plt.figure(figsize=(7, 4))
plt.plot(E, arr_tle[:, z_b, y_bragg, x_b], linewidth=1.5, color='red')
plt.xlabel("Énergie des PG (MeV)")
plt.ylabel("PG yields / primary")
plt.title(f'TLE optimisé — Pic de Bragg (Z={z_b}, Y={y_bragg}, X={x_b})')
plt.grid(True)
plt.gca().yaxis.set_major_formatter(sci_formatter())
plt.tight_layout()
plt.savefig("tle_bragg.png", dpi=150, bbox_inches='tight')
plt.close()

print("Terminé.")
print(f"Plateau  : Z={z_p}, Y={y_plateau}, X={x_p}")
print(f"Pic Bragg: Z={z_b}, Y={y_bragg},  X={x_b}")

