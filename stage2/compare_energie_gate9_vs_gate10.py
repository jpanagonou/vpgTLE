#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Comparaison visuelle Gate 9 / Gate 10 - spectre en energie.


Superpose les deux histogrammes d'energie detectee sur la meme figure,
pour une comparaison directe des formes de spectre (pas seulement les
moyennes).
"""


import numpy as np
import uproot
import matplotlib.pyplot as plt


gate10_root_path = "/home/jpanagonou/Desktop/vpgTLE-stage/stage2/gate10/output/gate10_realimage_phsp.root"
gate9_root_path = "/home/jpanagonou/Desktop/vpgTLE-stage/stage2/gate9/output/gate9_realimage_phsp.root"


Emin, Emax = 0.0, 10.0
nbBins_plot = 100  # nombre de bins pour l'affichage (peut differer de la source)


# ------------------------------------------------------------------
# lecture Gate 10 (directement depuis le .root, meme principe que Gate 9)
# ------------------------------------------------------------------
f10 = uproot.open(gate10_root_path)
tree10 = f10["PhaseSpace"]
energies_10_all = tree10["KineticEnergy"].array(library="np")
track_id_10 = tree10["TrackID"].array(library="np")
energies_10 = energies_10_all[track_id_10 == 1]  # ne garde que les primaires


# ------------------------------------------------------------------
# lecture Gate 9
# ------------------------------------------------------------------
f9 = uproot.open(gate9_root_path)
tree9 = f9["PhaseSpace"]
energies_9 = tree9["Ekine"].array(library="np")


print(f"Gate 10 : {len(energies_10)} particules, moyenne = {np.mean(energies_10):.3f} MeV, std = {np.std(energies_10):.3f} MeV")
print(f"Gate 9  : {len(energies_9)} particules, moyenne = {np.mean(energies_9):.3f} MeV, std = {np.std(energies_9):.3f} MeV")


# ------------------------------------------------------------------
# test statistique de Kolmogorov-Smirnov : compare les deux
# distributions completes (pas juste moyenne/ecart-type), donne une
# p-value quantifiant si les deux echantillons sont compatibles avec
# la meme distribution sous-jacente
# ------------------------------------------------------------------
from scipy.stats import ks_2samp
ks_stat, ks_pvalue = ks_2samp(energies_10, energies_9)
print(f"\nTest de Kolmogorov-Smirnov : statistique = {ks_stat:.4f}, p-value = {ks_pvalue:.4f}")


# ------------------------------------------------------------------
# histogrammes normalises (densite de probabilite, pour comparer des
# nombres d'evenements legerement differents entre les deux runs)
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))


bin_edges = np.linspace(Emin, Emax, nbBins_plot + 1)
bin_centers_plot = (bin_edges[:-1] + bin_edges[1:]) / 2


def plot_with_sigma_band(data, color, label, ax):
    """Trace la courbe avec une seule bande d'incertitude a +/- 3 sigma (Poisson)."""
    counts, _ = np.histogram(data, bins=bin_edges)
    sigma = np.sqrt(counts)


    borne_basse = np.maximum(counts - 3 * sigma, 0)  
    borne_haute = counts + 3 * sigma


    ax.fill_between(
        bin_centers_plot,
        borne_basse,
        borne_haute,
        color=color, alpha=0.25, linewidth=0, step="mid",
    )


    ax.hist(data, bins=bin_edges, histtype="step", color=color,
            linewidth=1.5, label=label)




def plot_simple(data, color, label, ax):
    """Trace uniquement la courbe, sans bande d'incertitude."""
    ax.hist(data, bins=bin_edges, histtype="step", color=color,
            linewidth=1.5, label=label)




# Gate 10 : notre developpement, avec la bande d'incertitude a 3 sigma
plot_with_sigma_band(energies_10, "tab:blue", f"Gate 10 (n={len(energies_10)})", ax)


# Gate 9 : reference, courbe simple sans bande
plot_simple(energies_9, "tab:red", f"Gate 9 (n={len(energies_9)})", ax)


ax.set_xlabel("Energie (MeV)")
ax.set_ylabel("Nombre de gammas (comptage)")
ax.set_ylim(bottom=0)

ax.set_title("Comparaison du spectre en energie detecte : Gate 9 vs Gate 10")
ax.legend()
ax.grid(alpha=0.3)




# encart texte avec les statistiques cles, directement sur la figure
stats_text = (
    f"Gate 10 : moy={np.mean(energies_10):.3f} MeV, std={np.std(energies_10):.3f} MeV\n"
    f"Gate 9  : moy={np.mean(energies_9):.3f} MeV, std={np.std(energies_9):.3f} MeV\n"
    f"Test KS : stat={ks_stat:.4f}, p-value={ks_pvalue:.4f}"
)
ax.text(
    0.98, 0.55, stats_text,
    transform=ax.transAxes, fontsize=9, family="monospace",
    verticalalignment="bottom", horizontalalignment="right",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
)


output_path_png = "comparaison_energie_gate9_gate10.png"
output_path_pdf = "comparaison_energie_gate9_gate10.pdf"
plt.savefig(output_path_png, dpi=150, bbox_inches="tight")
plt.savefig(output_path_pdf, bbox_inches="tight")
print(f"\nFigure sauvegardee : {output_path_png}")
print(f"Figure sauvegardee : {output_path_pdf}")


plt.show()
