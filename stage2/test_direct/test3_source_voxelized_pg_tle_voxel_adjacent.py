#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Test 3 de VoxelizedPromptGammaTLESource : voxels adjacents.


Principe : on construit une image 4D synthetique avec DEUX voxels actifs,
adjacents l'un a l'autre (differant d'un seul indice sur un seul axe),
chacun porteur d'une energie tres differente et facilement identifiable.


Ce test est plus strict que le test 2 (voxels espaces) : il verifie qu'il
n'y a pas de decalage d'indexation de type "off-by-one" a la frontiere
exacte entre deux voxels voisins. Un tel bug se traduirait par des
particules dont la position correspond a un voxel mais dont l'energie
correspond au voxel adjacent.
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


    # deux voxels ADJACENTS (meme iz, meme iy, ix differant d'une seule
    # unite), energies tres differentes pour detecter facilement toute
    # confusion a la frontiere
    voxels = [
        {"name": "D", "iz": 5, "iy": 5, "ix": 4, "E": 2.0},
        {"name": "E", "iz": 5, "iy": 5, "ix": 5, "E": 8.0},
    ]


    for v in voxels:
        ib = int((v["E"] - Emin) / bin_width)
        v["ib"] = ib
        arr_4d[ib - 1 : ib + 2, v["iz"], v["iy"], v["ix"]] = [0.2, 1.0, 0.2]
        v["expected_x"] = (v["ix"] - nX / 2.0 + 0.5) * spacing_mm
        v["expected_y"] = (v["iy"] - nY / 2.0 + 0.5) * spacing_mm
        v["expected_z"] = (v["iz"] - nZ / 2.0 + 0.5) * spacing_mm


    yield_per_proton = float(arr_4d.sum())
    print(f"Yield total (2 voxels adjacents) : {yield_per_proton}")
    print(f"Position D (ix={voxels[0]['ix']}) : x attendu = {voxels[0]['expected_x']:.2f} mm")
    print(f"Position E (ix={voxels[1]['ix']}) : x attendu = {voxels[1]['expected_x']:.2f} mm")


    # ecriture sur disque
    itk_img_4d = itk.image_from_array(arr_4d)
    itk_img_4d.SetSpacing([spacing_mm, spacing_mm, spacing_mm, 1.0])
    pg_image_path = str(paths.output / "pg_spectrum_4d_adjacent.nii.gz")
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


    # TrackID est ajoute pour pouvoir exclure d'eventuelles particules
    # secondaires (ex : electron Compton) dont le premier pas serait
    # aussi enregistre mais ne provient pas directement du tirage source.
    phsp = sim.add_actor("PhaseSpaceActor", "PhaseSpace")
    phsp.attached_to = box.name
    phsp.attributes = ["KineticEnergy", "PrePosition", "TrackID"]
    phsp.steps_to_store = "first"
    phsp.output_filename = "test_pg_tle_adjacent_phsp.root"


    sim.g4_commands_after_init.append("/tracking/verbose 0")


    sim.run()
    print(stats)


    # ------------------------------------------------------------------
    # verification : chaque particule doit correspondre au bon voxel,
    # meme a la frontiere entre les deux voxels adjacents
    # ------------------------------------------------------------------
    import uproot


    hits = uproot.open(phsp.get_output_path_string())["PhaseSpace"]
    energies = hits["KineticEnergy"].array(library="np")
    x = hits["PrePosition_X"].array(library="np")
    y = hits["PrePosition_Y"].array(library="np")
    z = hits["PrePosition_Z"].array(library="np")
    track_id = hits["TrackID"].array(library="np")


    # ne garde que les particules primaires (TrackID == 1), pour exclure
    # d'eventuelles particules secondaires (ex : electrons Compton) dont
    # le premier pas serait aussi enregistre par le PhaseSpaceActor mais
    # ne provient pas directement du tirage de la source
    is_primary = track_id == 1
    n_secondary = int((~is_primary).sum())
    if n_secondary > 0:
        print(f"\n{n_secondary} particule(s) secondaire(s) exclue(s) de l'analyse "
              f"(TrackID != 1, ex : electrons Compton)")


    energies = energies[is_primary]
    x = x[is_primary]
    y = y[is_primary]
    z = z[is_primary]


    n_total = len(energies)
    print(f"\nNombre total de particules primaires enregistrees : {n_total}")


    voxel_half = spacing_mm / 2.0
    tol_E = 3 * bin_width


    n_matched = 0
    n_mismatched = 0
    n_ambiguous = 0
    counts_per_voxel = {v["name"]: 0 for v in voxels}
    energies_per_voxel = {v["name"]: [] for v in voxels}
    positions_per_voxel = {v["name"]: {"x": [], "y": [], "z": []} for v in voxels}


    for i in range(n_total):
        # cherche a quel(s) voxel(s) cette particule correspond en
        # position (a la frontiere exacte, il est possible qu'aucun ou
        # que les deux matchent selon la tolerance choisie)
        matches = [
            v
            for v in voxels
            if abs(x[i] - v["expected_x"]) <= voxel_half
            and abs(y[i] - v["expected_y"]) <= voxel_half
            and abs(z[i] - v["expected_z"]) <= voxel_half
        ]


        if len(matches) == 0:
            n_mismatched += 1
            continue
        if len(matches) > 1:
            # la particule est dans la zone de tolerance des deux voxels
            # a la fois (frontiere exacte) : on verifie que l'energie
            # permet de trancher sans ambiguite
            matches_by_energy = [
                v for v in matches if abs(energies[i] - v["E"]) <= tol_E
            ]
            if len(matches_by_energy) != 1:
                n_ambiguous += 1
                continue
            matched_voxel = matches_by_energy[0]
        else:
            matched_voxel = matches[0]


        if abs(energies[i] - matched_voxel["E"]) <= tol_E:
            n_matched += 1
            counts_per_voxel[matched_voxel["name"]] += 1
            energies_per_voxel[matched_voxel["name"]].append(energies[i])
            positions_per_voxel[matched_voxel["name"]]["x"].append(x[i])
            positions_per_voxel[matched_voxel["name"]]["y"].append(y[i])
            positions_per_voxel[matched_voxel["name"]]["z"].append(z[i])
        else:
            n_mismatched += 1


    print(f"\nParticules avec paire position/energie coherente : {n_matched}")
    print(f"Particules incoherentes (position/energie ne correspondent pas) : {n_mismatched}")
    print(f"Particules ambigues (frontiere, non tranchees par l'energie) : {n_ambiguous}")
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


    is_ok = n_mismatched == 0 and n_ambiguous == 0
    all_voxels_represented = all(c > 0 for c in counts_per_voxel.values())


    utility.print_test(
        is_ok,
        f"Coherence position/energie aux voxels adjacents : {n_matched}/{n_total} "
        f"coherentes, {n_mismatched} incoherentes, {n_ambiguous} ambigues (attendu : 0, 0)",
    )
    utility.print_test(
        all_voxels_represented,
        f"Les deux voxels adjacents sont bien representes : {counts_per_voxel}",
    )


    utility.test_ok(is_ok and all_voxels_represented)
