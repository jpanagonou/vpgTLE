"""
plot_GammaZ.py
--------------
Trace le spectre GammaZ (rendement pondéré de gammas prompts)
pour les deux versions ROOT et sans ROOT.

Usage :
    python plot_GammaZ.py test.root PG_O.root

Dépendances :
    pip install uproot numpy matplotlib
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import uproot

if len(sys.argv) != 3:
    print("Usage : python plot_GammaZ.py <test.root> <PG_O.root>")
    sys.exit(1)

file_root      = sys.argv[1]
file_sans_root = sys.argv[2]

f1 = uproot.open(file_root)
f2 = uproot.open(file_sans_root)

# ── Lecture des histogrammes ──────────────────────────────────────────────────
h1 = f1["O/GammaZ"]
h2 = f2["GammaZ"]

v1, xedges1, yedges1 = h1.to_numpy()
v2, xedges2, yedges2 = h2.to_numpy()

# ── Fonction de tracé ─────────────────────────────────────────────────────────
def plot_GammaZ(values, xedges, yedges, title, filename):
    fig, ax = plt.subplots(figsize=(8, 5))

    # Masquer les zéros pour l'échelle log
    values_masked = np.ma.masked_where(values <= 0, values)

   # vmin = values_masked.min()
   # vmax = values_masked.max()

    vmin = 1e-12
    vmax = 1e-4

    im = ax.pcolormesh(
        xedges, yedges, values_masked.T,
        cmap='jet',
        norm=mcolors.LogNorm(vmin=vmin, vmax=vmax)
    )

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r"PG per cm per 40 keV bin", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    ax.set_title(title, fontsize=11, pad=10)
    ax.set_xlabel(r"Proton energy $E_p$ (MeV)", fontsize=11)
    ax.set_ylabel(r"Prompt-gamma energy $E_{pg}$ (MeV)", fontsize=11)

    ax.set_xlim(xedges[0], xedges[-1])
    ax.set_ylim(yedges[0], yedges[-1])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    print(f"Sauvegarde : {filename}")
    plt.close()

# ── Tracé des deux versions ───────────────────────────────────────────────────
plot_GammaZ(
    v1, xedges1, yedges1,
    title=r"Oxygen $-$ $\mathbf{E_{pg}}(E_p)$ $-$ Geant4/ROOT",
    filename="GammaZ_ROOT.png"
)

plot_GammaZ(
    v2, xedges2, yedges2,
    title=r"Oxygen $-$ $\mathbf{E_{pg}}(E_p)$ $-$ Geant4/G4AnalysisManager",
    filename="GammaZ_sansROOT.png"
)

print("Termine.")

