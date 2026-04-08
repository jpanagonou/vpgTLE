
import itk
import numpy as np
import matplotlib.pyplot as plt

# Charger l'image
img = itk.imread('output/stage1a/VPG_vpg_gamma_e.nii.gz')
arr = itk.array_from_image(img)

print(f"Dimensions du volume : {arr.shape[1]} x {arr.shape[2]} x {arr.shape[3]}")

# Juste un texte pour trouver les voxels non nul

nonzero = np.argwhere(arr.sum(axis=0) != 0)
print('Voxels non nuls:', len(nonzero))
if len(nonzero) > 0:
    print('Exemples de voxels actifs:')
    for v in nonzero[:5]:
        print(f'  Z={v[0]}, Y={v[1]}, X={v[2]}')



# Entrer les coordonnées
z = int(input("Entrez Z : "))
y = int(input("Entrez Y : "))
x = int(input("Entrez X : "))

# Extraire le spectre
spectre = arr[:, z, y, x]
E = np.linspace(0, 10, arr.shape[0])

# Tracer
plt.figure()
plt.plot(E, spectre)
plt.xlabel('Energie gamma (MeV)')
plt.ylabel('Emission gamma')
plt.title(f'Spectre gamma - voxel ({x}, {y}, {z})')
plt.grid(True)
plt.savefig(f'spectre_gamma_{x}_{y}_{z}.png')
plt.show()
print(f"Courbe sauvegardée : spectre_gamma_{x}_{y}_{z}.png")


