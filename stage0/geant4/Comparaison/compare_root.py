"""
compare_root.py
---------------
Comparaison bin à bin entre deux fichiers .root :
  - test.root  : version originale avec ROOT (histos dans répertoire O/)
  - PG_O.root  : version sans ROOT / G4AnalysisManager (histos à la racine)


Usage :
    python compare_root.py test.root PG_O.root


Dépendances :
    pip install uproot numpy matplotlib
"""


import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import uproot


plt.rcParams['axes.spines.top']   = False
plt.rcParams['axes.spines.right'] = False
    



# ── Lecture des arguments ─────────────────────────────────────────────────────
if len(sys.argv) != 3:
    print("Usage : python compare_root.py <fichier_root.root> <fichier_sans_root.root>")
    sys.exit(1)


file_root      = sys.argv[1]
file_sans_root = sys.argv[2]


print(f"Fichier ROOT original  : {file_root}")
print(f"Fichier sans ROOT      : {file_sans_root}")


# ── Ouverture des fichiers ────────────────────────────────────────────────────
f1 = uproot.open(file_root)
f2 = uproot.open(file_sans_root)


print("\nContenu fichier ROOT original :")
print(f1.keys())
print("\nContenu fichier sans ROOT :")
print(f2.keys())


# ── Correspondance des noms ───────────────────────────────────────────────────
PAIRS_H2 = [
    ("O/EpEpg",           "EpEpg",  "EpEpg",
     "Proton energy [MeV]", "PG energy [MeV]"),
    ("O/GammaZ",          "GammaZ", "GammaZ",
     "Proton energy [MeV]", "PG energy [MeV]"),
]


PAIRS_H1 = [
    ("O/NrPG",             "NrPG",  "NrPG",
     "Proton energy [MeV]", r"$N_{\gamma}$"),
    ("O/Kapa inelastique", "Sigma", "Sigma / Kapa inelastique",
     "Proton energy [MeV]", r"$\kappa_{\mathrm{inel}}$ [cm$^{-1}$]"),
]


# ── Fonctions utilitaires ─────────────────────────────────────────────────────
def delta_relative(v1, v2):
    with np.errstate(divide='ignore', invalid='ignore'):
        d = np.where(v1 != 0, (v2 - v1) / v1, np.nan)
    return d


def plot_diff_2D(ax, delta, xedges, yedges, title, xlabel, ylabel):
    vmax = np.nanpercentile(np.abs(delta), 99)
    vmax = max(vmax, 1e-9)
    im = ax.pcolormesh(xedges, yedges, delta.T,
                       cmap='RdBu_r',
                       norm=mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
    plt.colorbar(im, ax=ax, label=r'$\delta$ relatif')
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def plot_diff_1D(ax, delta, centers, title, xlabel):
    ax.plot(centers, delta, color='steelblue', linewidth=1.0)
    ax.axhline(0, color='red', linestyle='--', linewidth=0.8)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$\delta$ relatif')
   


def plot_overlay_1D(ax, v1, v2, centers, title, xlabel, ylabel):
    ax.plot(centers, v1, label='Avec ROOT',  color='black',     linewidth=1.2)
    ax.plot(centers, v2, label='Sans ROOT',  color='steelblue', linewidth=1.0, linestyle='--')
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=8)
    


def print_metrics(name, d_flat):
    print(f"\n--- {name} ---")
    print(f"  delta moyen       : {np.nanmean(d_flat):.4e}")
    print(f"  delta ecart-type  : {np.nanstd(d_flat):.4e}")
    print(f"  |delta| max       : {np.nanmax(np.abs(d_flat)):.4e}")
    frac_zero = np.sum(np.abs(d_flat) < 1e-6) / len(d_flat) * 100
    print(f"  bins identiques   : {frac_zero:.1f} %")


# ── Comparaison H2 ────────────────────────────────────────────────────────────
for name1, name2, label, xlabel, ylabel in PAIRS_H2:
    try:
        h1 = f1[name1]
        h2 = f2[name2]
    except KeyError as e:
        print(f"[WARNING] Histogramme introuvable : {e}, ignore.")
        continue

    v1, xedges, yedges = h1.to_numpy()
    v2, _,      _      = h2.to_numpy()
    delta    = delta_relative(v1, v2)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Comparaison H2 - {label}", fontsize=13)

    plot_diff_2D(axes[0], delta, xedges, yedges,
                 f"Carte delta - {label}", xlabel, ylabel)

    d_flat = delta[~np.isnan(delta)].flatten()
    axes[1].hist(d_flat, bins=100, color='steelblue', edgecolor='none')
    axes[1].set_title(f"Distribution de delta - {label}", fontsize=10)
    axes[1].set_xlabel("delta relatif")
    axes[1].set_ylabel("Nombre de bins")
    axes[1].axvline(0, color='red', linestyle='--')
    

    plt.tight_layout()
    out = f"compare_{name2}.png"
    plt.savefig(out, dpi=150)
    print(f"Sauvegarde : {out}")
    plt.close()

    print_metrics(label, d_flat)


# ── Comparaison H1 ────────────────────────────────────────────────────────────
n = len(PAIRS_H1)
fig_h1, axes_h1 = plt.subplots(n, 2, figsize=(14, 5 * n))
if n == 1:
    axes_h1 = [axes_h1]
fig_h1.suptitle("Comparaison H1", fontsize=13)


for idx, (name1, name2, label, xlabel, ylabel) in enumerate(PAIRS_H1):
    try:
        h1 = f1[name1]
        h2 = f2[name2]
    except KeyError as e:
        print(f"[WARNING] Histogramme introuvable : {e}, ignore.")
        continue


    result1 = h1.to_numpy()
    result2 = h2.to_numpy()
    v1      = result1[0]
    edges   = result1[1]
    v2      = result2[0]

    
    centers   = 0.5 * (edges[:-1] + edges[1:])
    delta     = delta_relative(v1, v2)


    plot_overlay_1D(axes_h1[idx][0], v1, v2, centers,
                    f"Superposition - {label}", xlabel, ylabel)
    plot_diff_1D(axes_h1[idx][1], delta, centers,
                 f"Difference relative - {label}", xlabel)


    d_flat = delta[~np.isnan(delta)]
    print_metrics(label, d_flat)


plt.tight_layout()
out_h1 = "compare_H1.png"
plt.savefig(out_h1, dpi=150)
print(f"\nSauvegarde : {out_h1}")
plt.close()


print("\nComparaison terminee.")







