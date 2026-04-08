import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import uproot

file_gate9  = "/home/jpanagonou/Desktop/vpgTLE/stages/stage0/gate9/db-Oxygen.root"
file_gate10 = "/home/jpanagonou/Desktop/vpgTLE/stages/stage0/gate10/output/PGdb_O.root"

f9  = uproot.open(file_gate9)
f10 = uproot.open(file_gate10)

h9  = f9["Oxygen/EpEpg"]
h10 = f10["EpEpg"]

v9,  xedges9,  yedges9  = h9.to_numpy()
v10, xedges10, yedges10 = h10.to_numpy()

def plot_EpEpg(values, xedges, yedges, title, filename):
    fig, ax = plt.subplots(figsize=(8, 5))
    values_masked = np.ma.masked_where(values <= 0, values)
    vmin = values_masked.min()
    vmax = values_masked.max()
    im = ax.pcolormesh(
        xedges, yedges, values_masked.T,
        cmap='jet',
        norm=mcolors.LogNorm(vmin=vmin, vmax=vmax)
    )
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=9)
    ax.set_title(title, fontsize=12, pad=10)
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

plot_EpEpg(
    v9, xedges9, yedges9,
    title=r"Oxygen $-$ EpEpg $-$ Gate 9",
    filename="EpEpg_gate9.png"
)

plot_EpEpg(
    v10, xedges10, yedges10,
    title=r"Oxygen $-$ EpEpg $-$ Gate 10",
    filename="EpEpg_gate10.png"
)

print("Terminé.")

