#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Test 2 de VoxelizedPromptGammaTLESource : plusieurs voxels actifs,
energies distinctes.

Principe : on construit une image 4D synthetique avec TROIS voxels actifs,
places a des coordonnees (iz, iy, ix) toutes differentes deux a deux sur
chaque axe (aucun ne partage la meme valeur de iz, iy ou ix qu'un autre).
Chaque voxel recoit un spectre etroit centre sur une energie distincte et
facilement identifiable.

Si le couplage position <-> energie est correct, chaque particule generee
doit presenter une paire (position, energie) coherente avec UN SEUL des
trois voxels, jamais un melange (ex : position du voxel A avec l'energie
du voxel B). Un bug d'inversion d'axes (Y/Z confondus par exemple) se
traduirait par des positions ne correspondant a aucun des trois centres
attendus.
"""

import numpy as np
import itk
import opengate as gate
from opengate.tests import utility

if __name__ == "__main__":
    paths = utility.get_default_test_paths(__file__, output_folder="test_pg_tle")

    # ------------------------------------------------------------------
    # construction de l'image 4D synthetique : [nbBins, Z, Y, X]
    # ------------------------------------------------------------------
    nZ, nY, nX = 10, 10, 10
    nbBins = 250
    Emin, Emax = 0.0, 10.0
    spacing_mm = 2.0
    bin_width = (Emax - Emin) / nbBins

    arr_4d = np.zeros((nbBins, nZ, nY, nX), dtype=np.float32)

    # trois voxels actifs, coordonnees espacees d'au moins 2 indices sur
    # chaque axe (garantit un intervalle vide entre chaque groupe sur les
    # histogrammes 1D, en plus de rester toutes distinctes deux a deux)
    voxels = [
        {"name": "A", "iz": 1, "iy": 1, "ix": 7, "E": 2.0},
        {"name": "B", "iz": 4, "iy": 7, "ix": 1, "E": 5.0},
        {"name": "C", "iz": 7, "iy": 4, "ix": 4, "E": 8.0},
    ]

    for v in voxels:
        ib = int((v["E"] - Emin) / bin_width)
        v["ib"] = ib
        arr_4d[ib - 1 : ib + 2, v["iz"], v["iy"], v["ix"]] = [0.2, 1.0, 0.2]
        # position physique attendue (image centree en 0)
        v["expected_x"] = (v["ix"] - nX / 2.0 + 0.5) * spacing_mm
        v["expected_y"] = (v["iy"] - nY / 2.0 + 0.5) * spacing_mm
        v["expected_z"] = (v["iz"] - nZ / 2.0 + 0.5) * spacing_mm

    yield_per_proton = float(arr_4d.sum())
    print(f"Yield total (3 voxels) : {yield_per_proton}")

    # ecriture sur disque
    itk_img_4d = itk.image_from_array(arr_4d)
    itk_img_4d.SetSpacing([spacing_mm, spacing_mm, spacing_mm, 1.0])
    pg_image_path = str(paths.output / "pg_spectrum_4d_multi.nii.gz")
    itk.imwrite(itk_img_4d, pg_image_path)

    # ------------------------------------------------------------------
    # simulation
    # ------------------------------------------------------------------
    sim = gate.Simulation()

    sim.g4_verbose = False
    sim.visu = False
    sim.number_of_threads = 1
    sim.random_seed = 123456
    sim.output_dir = paths.output

    mm = gate.g4_units.mm
    m = gate.g4_units.m

    sim.volume_manager.add_material_database(paths.data / "GateMaterials.db")
    sim.world.size = [1 * m, 1 * m, 1 * m]
    sim.world.material = "G4_AIR"

    half_extent = (max(nX, nY, nZ) * spacing_mm) / 2.0 + 5 * mm
    box = sim.add_volume("Box", "box_pg")
    box.size = [2 * half_extent, 2 * half_extent, 2 * half_extent]
    box.material = "G4_AIR"

    source = sim.add_source("VoxelizedPromptGammaTLESource", "pg_source")
    source.attached_to = box.name
    source.particle = "gamma"
    source.n_protons = 1000
    source.image = pg_image_path
    source.Emin = Emin
    source.Emax = Emax
    source.nbBins = nbBins
    source.direction.type = "iso"

    sim.physics_manager.physics_list_name = "G4EmStandardPhysics"
    sim.physics_manager.global_production_cuts.all = 1 * mm

    stats = sim.add_actor("SimulationStatisticsActor", "Stats")

    phsp = sim.add_actor("PhaseSpaceActor", "PhaseSpace")
    phsp.attached_to = box.name
    phsp.attributes = ["KineticEnergy", "PrePosition"]
    phsp.steps_to_store = "first"
    phsp.output_filename = "test_pg_tle_multi_phsp.root"

    sim.g4_commands_after_init.append("/tracking/verbose 0")

    sim.run()
    print(stats)

    # ------------------------------------------------------------------
    # verification : chaque particule doit correspondre a UN SEUL voxel,
    # avec la bonne paire (position, energie)
    # ------------------------------------------------------------------
    import uproot

    hits = uproot.open(phsp.get_output_path_string())["PhaseSpace"]
    energies = hits["KineticEnergy"].array(library="np")
    x = hits["PrePosition_X"].array(library="np")
    y = hits["PrePosition_Y"].array(library="np")
    z = hits["PrePosition_Z"].array(library="np")

    n_total = len(energies)
    print(f"\nNombre total de particules enregistrees : {n_total}")

    voxel_half = spacing_mm / 2.0
    tol_E = 3 * bin_width

    n_matched = 0
    n_mismatched = 0
    counts_per_voxel = {v["name"]: 0 for v in voxels}
    positions_per_voxel = {v["name"]: {"x": [], "y": [], "z": []} for v in voxels}
    energies_per_voxel = {v["name"]: [] for v in voxels}

    for i in range(n_total):
        # cherche a quel voxel cette particule correspond en position
        matched_voxel = None
        for v in voxels:
            if (
                abs(x[i] - v["expected_x"]) <= voxel_half
                and abs(y[i] - v["expected_y"]) <= voxel_half
                and abs(z[i] - v["expected_z"]) <= voxel_half
            ):
                matched_voxel = v
                break

        if matched_voxel is None:
            # position ne correspond a aucun des 3 voxels attendus
            n_mismatched += 1
            continue

        # verifie que l'energie correspond bien A CE MEME voxel
        if abs(energies[i] - matched_voxel["E"]) <= tol_E:
            n_matched += 1
            counts_per_voxel[matched_voxel["name"]] += 1
            positions_per_voxel[matched_voxel["name"]]["x"].append(x[i])
            positions_per_voxel[matched_voxel["name"]]["y"].append(y[i])
            positions_per_voxel[matched_voxel["name"]]["z"].append(z[i])
            energies_per_voxel[matched_voxel["name"]].append(energies[i])
        else:
            n_mismatched += 1

    print(f"\nParticules avec paire position/energie coherente : {n_matched}")
    print(f"Particules incoherentes (position/energie ne correspondent pas) : {n_mismatched}")
    print(f"Repartition par voxel : {counts_per_voxel}")

    print("\nPosition moyenne mesuree par voxel (mm) :")
    for v in voxels:
        px = positions_per_voxel[v["name"]]["x"]
        py = positions_per_voxel[v["name"]]["y"]
        pz = positions_per_voxel[v["name"]]["z"]
        if len(px) > 0:
            print(
                f"  Voxel {v['name']} : attendu "
                f"({v['expected_x']:.2f}, {v['expected_y']:.2f}, {v['expected_z']:.2f}), "
                f"mesure moyenne "
                f"({np.mean(px):.2f}, {np.mean(py):.2f}, {np.mean(pz):.2f})"
            )

    print("\nEnergie moyenne mesuree par voxel (MeV) :")
    for v in voxels:
        e = energies_per_voxel[v["name"]]
        if len(e) > 0:
            print(
                f"  Voxel {v['name']} : attendu {v['E']:.2f} MeV, "
                f"mesure moyenne {np.mean(e):.3f} MeV "
                f"(min={np.min(e):.3f}, max={np.max(e):.3f})"
            )

    is_ok = n_mismatched == 0
    all_voxels_represented = all(c > 0 for c in counts_per_voxel.values())

    utility.print_test(
        is_ok,
        f"Coherence position/energie : {n_matched}/{n_total} particules "
        f"coherentes, {n_mismatched} incoherentes (attendu : 0)",
    )
    utility.print_test(
        all_voxels_represented,
        f"Tous les voxels representes dans les detections : {counts_per_voxel}",
    )

    utility.test_ok(is_ok and all_voxels_represented)
