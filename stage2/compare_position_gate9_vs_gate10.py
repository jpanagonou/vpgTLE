#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Comparaison visuelle Gate 9 / Gate 10 - positions X, Y, Z.


Superpose, pour chaque axe separement, les deux histogrammes de position
detectee (Gate 9 vs Gate 10), avec moyenne, ecart-type et test de
Kolmogorov-Smirnov affiches directement sur chaque figure.
"""


import numpy as np
import uproot
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp


# --- A ADAPTER selon tes chemins reels ---
gate10_root_path = "/home/jpanagonou/Desktop/vpgTLE-stage/stage2/gate10/output/gate10_realimage_phsp.root"
gate9_root_path = "/home/jpanagonou/Desktop/vpgTLE-stage/stage2/gate9/output/gate9_realimage_phsp.root"


pos_range = (-1000, 1000)  # mm, rayon de la sphere = 1000mm
nbBins_plot = 100


# ------------------------------------------------------------------
# lecture Gate 10 (branches PrePosition_X/Y/Z, filtre TrackID==1)
# ------------------------------------------------------------------
f10 = uproot.open(gate10_root_path)
tree10 = f10["PhaseSpace"]
track_id_10 = tree10["TrackID"].array(library="np")
is_primary_10 = track_id_10 == 1


positions_10 = {
    "X": tree10["PrePosition_X"].array(library="np")[is_primary_10],
    "Y": tree10["PrePosition_Y"].array(library="np")[is_primary_10],
    "Z": tree10["PrePosition_Z"].array(library="np")[is_primary_10],
}


# ------------------------------------------------------------------
# lecture Gate 9 (branches X, Y, Z, telles que nommees par
# enableXPosition/enableYPosition/enableZPosition dans la macro)
# ------------------------------------------------------------------
f9 = uproot.open(gate9_root_path)
tree9 = f9["PhaseSpace"]


positions_9 = {
    "X": tree9["X"].array(library="np"),
    "Y": tree9["Y"].array(library="np"),
    "Z": tree9["Z"].array(library="np"),
}


# ------------------------------------------------------------------
# une figure par axe
# ------------------------------------------------------------------
for axis in ["X", "Y", "Z"]:
    data10 = positions_10[axis]
    data9 = positions_9[axis]


    mean10, std10 = np.mean(data10), np.std(data10)
    mean9, std9 = np.mean(data9), np.std(data9)
    ks_stat, ks_pvalue = ks_2samp(data10, data9)


    print(f"\n=== Axe {axis} ===")
    print(f"Gate 10 : n={len(data10)}, moyenne={mean10:.3f} mm, std={std10:.3f} mm")
    print(f"Gate 9  : n={len(data9)}, moyenne={mean9:.3f} mm, std={std9:.3f} mm")
    print(f"Test KS : statistique={ks_stat:.4f}, p-value={ks_pvalue:.4f}")


    fig, ax = plt.subplots(figsize=(10, 6))


    bin_edges = np.linspace(pos_range[0], pos_range[1], nbBins_plot + 1)
    bin_centers_plot = (bin_edges[:-1] + bin_edges[1:]) / 2


    # Gate 10 : bande d'incertitude a 3 sigma (Poisson), notre developpement
    counts10, _ = np.histogram(data10, bins=bin_edges)
    sigma10 = np.sqrt(counts10)
    borne_basse = np.maximum(counts10 - 2 * sigma10, 0)
    borne_haute = counts10 + 2 * sigma10
    ax.fill_between(bin_centers_plot, borne_basse, borne_haute,
                     color="tab:blue", alpha=0.25, linewidth=0, step="mid")
    ax.hist(data10, bins=bin_edges, histtype="step", color="tab:blue",
            linewidth=1.5, label=f"Gate 10 (n={len(data10)})")


    # Gate 9 : reference, courbe simple sans bande
    ax.hist(data9, bins=bin_edges, histtype="step", color="tab:red",
            linewidth=1.5, label=f"Gate 9 (n={len(data9)})")


    ax.set_xlabel(f"Position {axis} (mm)")
    ax.set_ylabel("Nombre de gammas (comptage)")
    ax.set_ylim(bottom=0)
    ax.set_title(f"Comparaison de la position {axis} detectee : Gate 9 vs Gate 10")
    ax.legend()
    ax.grid(alpha=0.3)


    stats_text = (
        f"Gate 10 : moy={mean10:.3f} mm, std={std10:.3f} mm\n"
        f"Gate 9  : moy={mean9:.3f} mm, std={std9:.3f} mm\n"
        f"Test KS : stat={ks_stat:.4f}, p-value={ks_pvalue:.4f}"
    )
    ax.text(
        0.98, 0.95, stats_text,
        transform=ax.transAxes, fontsize=9, family="monospace",
        verticalalignment="top", horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )


    output_path_png = f"comparaison_position_{axis}_gate9_gate10.png"
    output_path_pdf = f"comparaison_position_{axis}_gate9_gate10.pdf"
    plt.savefig(output_path_png, dpi=150, bbox_inches="tight")
    plt.savefig(output_path_pdf, bbox_inches="tight")
    print(f"Figure sauvegardee : {output_path_png}")
    print(f"Figure sauvegardee : {output_path_pdf}")



plt.show()
