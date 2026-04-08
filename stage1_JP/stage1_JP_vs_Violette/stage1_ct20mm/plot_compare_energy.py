
import itk
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Charger les deux images
arr_v  = itk.array_from_image(itk.imread('output/stage1a/VPG_vpg_violette_gamma_e.nii.gz'))
arr_jp = itk.array_from_image(itk.imread('output/stage1a/VPG_vpg_jeanpaul_gamma_e.nii.gz'))

print(f"Shape : {arr_v.shape[1]} x {arr_v.shape[2]} x {arr_v.shape[3]}")

# Voxels actifs dans les deux
mask_v  = arr_v.sum(axis=0) != 0
mask_jp = arr_jp.sum(axis=0) != 0
mask_commun = mask_v & mask_jp
nonzero_commun = np.argwhere(mask_commun)

print(f"\nVoxels actifs Violette : {mask_v.sum()}")
print(f"Voxels actifs Jean-Paul      : {mask_jp.sum()}")
print(f"Voxels actifs communs  : {len(nonzero_commun)}")
print("\nExemples de voxels communs :")
for v in nonzero_commun[:5]:
    print(f"  Z={v[0]}, Y={v[1]}, X={v[2]}")

# Choisir le voxel
z = int(input("\nEntrez Z : "))
y = int(input("Entrez Y : "))
x = int(input("Entrez X : "))

# Extraire les spectres
spectre_v  = arr_v[:, z, y, x]
spectre_jp = arr_jp[:, z, y, x]
E = np.linspace(0, 10, arr_v.shape[0])

# Comparaison numérique
print(f"\nMax difference : {np.max(np.abs(spectre_v.astype(float) - spectre_jp.astype(float))):.6f}")
print(f"Identiques     : {np.allclose(spectre_v, spectre_jp)}")

# Tracer
plt.figure(figsize=(10, 5))
plt.plot(E, spectre_v,  label='Old(original)', linewidth=1.5)
plt.plot(E, spectre_jp, label='New (optimised)', linestyle='--', linewidth=1.5)
plt.xlabel('Energie gamma (MeV)')
plt.ylabel('Emission gamma')
plt.title(f'Comparaison spectre gamma - voxel ({x}, {y}, {z})')
plt.legend()
plt.savefig(f'compare_spectre_{x}_{y}_{z}.png')
print(f"Courbe sauvegardée : compare_spectre_{x}_{y}_{z}.png")


