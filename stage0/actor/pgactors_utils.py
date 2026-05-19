# --------------------------------------------------
# Copyright (C): OpenGATE Collaboration
# This software is distributed under the terms
# of the GNU Lesser General Public Licence (LGPL)
# See LICENSE.md for further details
# --------------------------------------------------

"""
Utility functions for post-processing the output of
PromptGammaStatisticActor simulations.
"""

import os
import numpy as np
from pathlib import Path

try:
    import uproot
except ImportError:
    raise ImportError(
        "uproot is required for post-processing. "
        "Install it with: pip install uproot"
    )


def reorganize_root_file(filepath):
    """
    Reorganize flat ROOT histograms into sub-folders per element.

    G4AnalysisManager does not support sub-folder creation during
    simulation. This function reads the flat ROOT file produced by
    PromptGammaStatisticActor and reorganizes the histograms into
    sub-folders per element (e.g. Oxygen/GammaZ, Oxygen/EpEpg...).

    This function is required after every simulation.

    Parameters
    ----------
    filepath : str or Path
        Path to the ROOT file to reorganize.

    Example
    -------
    >>> reorganize_root_file("output/PGdb_Oxygen_proton.root")
    """
    filepath = str(filepath)
    tmp_file = filepath + "_tmp.root"

    with uproot.open(filepath) as fin:
        data = {}
        for name, obj in fin.items(recursive=True, cycle=False):
            try:
                data[name] = obj.to_numpy()
            except AttributeError:
                pass

    with uproot.recreate(tmp_file) as fout:
        for name, hist_data in data.items():
            fout[name] = hist_data

    os.remove(filepath)
    os.rename(tmp_file, filepath)
    print(f"[pgactors_utils] ROOT file reorganized : {filepath}")


def convert_secondaries_to_root(output_dir, material_name,
                                 particle_type):
    """
    Convert secondaries txt files (KE > 0 and KE = 0) to a single
    ROOT file with two sub-folders KE0_partial and KE0_total.

    This function is required only when save_KE0_secondaries = True.

    Parameters
    ----------
    output_dir : str or Path
        Directory containing the txt files.
    material_name : str
        Name of the material used in the simulation.
    particle_type : str
        Type of incident particle: 'proton', 'neutron',
        'helium', 'carbon'.

    Output
    ------
    ROOT file : PGdb_<material_name>_<particle_type>_secondaries.root
        ├── KE0_partial/   ← secondaries from partial interactions
        │   ├── e-
        │   └── ...
        └── KE0_total/     ← secondaries from total interactions
            ├── proton
            ├── alpha
            └── ...

    Example
    -------
    >>> convert_secondaries_to_root(
    ...     Path("output"), "Oxygen", "helium")
    """
    output_dir = Path(output_dir)

    # ── Check that at least one txt file exists ───────────────────
    txt_partial = output_dir / \
        f"PGdb_{material_name}_{particle_type}_KE0_secondaries.txt"
    txt_total = output_dir / \
        f"PGdb_{material_name}_{particle_type}_KE0total_secondaries.txt"

    if not txt_partial.exists() and not txt_total.exists():
        print(f"[pgactors_utils] No secondaries txt files found "
              f"for {particle_type} on {material_name}. "
              f"Make sure save_KE0_secondaries = True was set.")
        return

    secondaries_data = {"KE0_partial": {}, "KE0_total": {}}

    # ── Read KE > 0 ───────────────────────────────────────────────
    if txt_partial.exists() and txt_partial.stat().st_size > 0:
        with open(txt_partial) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) != 2:
                    continue
                name, energy = parts
                if name not in secondaries_data["KE0_partial"]:
                    secondaries_data["KE0_partial"][name] = []
                secondaries_data["KE0_partial"][name].append(
                    float(energy))
    else:
        print(f"[pgactors_utils] No partial interaction (KE>0) "
              f"detected for {particle_type}.")

    # ── Read KE = 0 ───────────────────────────────────────────────
    if txt_total.exists() and txt_total.stat().st_size > 0:
        with open(txt_total) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) != 2:
                    continue
                name, energy = parts
                if name not in secondaries_data["KE0_total"]:
                    secondaries_data["KE0_total"][name] = []
                secondaries_data["KE0_total"][name].append(
                    float(energy))

    # ── Write to single ROOT file ─────────────────────────────────
    root_out = output_dir / \
        f"PGdb_{material_name}_{particle_type}_secondaries.root"
    with uproot.recreate(str(root_out)) as fout:
        for folder, particles in secondaries_data.items():
            if not particles:
                hist  = np.zeros(500, dtype=np.int64)
                edges = np.linspace(0.0, 200.0, 501)
                fout[f"{folder}/empty"] = (hist, edges)
                continue
            for name, energies in particles.items():
                arr = np.array(energies)
                hist, edges = np.histogram(arr, bins=500,
                                           range=(0.0, 200.0))
                fout[f"{folder}/{name}"] = (hist, edges)

    print(f"[pgactors_utils] Secondaries ROOT file created : "
          f"{root_out}")
    print(f"  KE>0 : "
          f"{sum(len(v) for v in secondaries_data['KE0_partial'].values())}"
          f" entries")
    print(f"  KE=0 : "
          f"{sum(len(v) for v in secondaries_data['KE0_total'].values())}"
          f" entries")
