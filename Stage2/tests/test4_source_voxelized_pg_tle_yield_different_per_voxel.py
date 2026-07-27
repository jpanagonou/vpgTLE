#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Test 4 de VoxelizedPromptGammaTLESource : coherence du rendement spatial
et fidelite du tirage sur des spectres non triviaux.


Principe : quatre voxels actifs, bien separes spatialement, avec des
rendements (yields) dans un rapport 10:5:2:1. Si le tirage spatial par
CDF est correct, le nombre de particules issues de chaque voxel doit
respecter ce meme rapport, aux fluctuations statistiques pres.


Chaque voxel porte de plus un spectre en energie non trivial, different
d'un simple pic etroit :
 - W (rendement le plus fort) : profil gaussien sur ~4 sigma
 - X, Y, Z : profil uniforme (plat) sur 20 bins consecutifs


Ceci permet de verifier, en plus du ratio de comptage entre voxels, que
le tirage par inversion de CDF respecte fidelement la moyenne et
l'ecart-type theoriques de chaque forme de spectre, pas seulement ses
bornes.
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


    base_yield = 1.4  # unite de rendement de reference
    sigma_E = 0.3      # ecart-type de la gaussienne du voxel W (MeV)
    n_bins_flat = 40   # nombre de bins actifs vises pour les profils plats


    # quatre voxels bien separes (aucune adjacence), rendements dans un
    # rapport 10:5:2:1, energies choisies pour que les quatre fenetres
    # spectrales ne se chevauchent JAMAIS entre elles
    voxels = [
        {"name": "W", "iz": 1, "iy": 1, "ix": 1, "E": 2.0, "ratio": 10, "shape": "gauss"},
        {"name": "X", "iz": 1, "iy": 1, "ix": 8, "E": 4.2, "ratio": 5, "shape": "flat"},
        {"name": "Y", "iz": 8, "iy": 1, "ix": 1, "E": 6.2, "ratio": 2, "shape": "flat"},
        {"name": "Z", "iz": 1, "iy": 8, "ix": 8, "E": 8.2, "ratio": 1, "shape": "flat"},
    ]


    for v in voxels:
        ib = int((v["E"] - Emin) / bin_width)
        v["target_yield"] = v["ratio"] * base_yield


        # fenetre de bins VISEE (peut deborder du tableau si l'energie
        # centrale est proche d'un bord)
        if v["shape"] == "gauss":
            half_width_bins = int(4 * sigma_E / bin_width)
        else:
            half_width_bins = n_bins_flat // 2
        full_range = np.arange(ib - half_width_bins, ib + half_width_bins + 1)


        # troncature automatique aux bords du tableau : on ne garde que
        # les indices valides (0 <= indice < nbBins), sans planter ni
        # deborder silencieusement si l'energie est proche d'un bord
        valid_mask = (full_range >= 0) & (full_range < nbBins)
        bins_range = full_range[valid_mask]
        n_troncature = len(full_range) - len(bins_range)
        if n_troncature > 0:
            print(f"  [Voxel {v['name']}] fenetre tronquee de {n_troncature} "
                  f"bin(s) au(x) bord(s) du tableau (proche de {v['E']} MeV)")


        bin_centers_E = Emin + (bins_range + 0.5) * bin_width


        if v["shape"] == "gauss":
            profile = np.exp(-0.5 * ((bin_centers_E - v["E"]) / sigma_E) ** 2)
        else:  # "flat"
            profile = np.ones(len(bins_range))


        # normalisation sur les bins REELLEMENT disponibles (apres
        # troncature), pour conserver le rendement total cible malgre
        # une eventuelle coupure au bord
        profile = profile / profile.sum() * v["target_yield"]


        arr_4d[bins_range, v["iz"], v["iy"], v["ix"]] = profile


        # moyenne et ecart-type THEORIQUES calcules depuis le profil
        # reellement construit (donc correct meme si tronque)
        weights = profile / profile.sum()
        v["theory_mean"] = float(np.sum(bin_centers_E * weights))
        v["theory_std"] = float(
            np.sqrt(np.sum(weights * (bin_centers_E - v["theory_mean"]) ** 2))
        )
        v["bins_range"] = bins_range


        v["expected_x"] = (v["ix"] - nX / 2.0 + 0.5) * spacing_mm
        v["expected_y"] = (v["iy"] - nY / 2.0 + 0.5) * spacing_mm
        v["expected_z"] = (v["iz"] - nZ / 2.0 + 0.5) * spacing_mm


    yield_per_proton = float(arr_4d.sum())
    print(f"Yield total (4 voxels) : {yield_per_proton:.3f}")
    for v in voxels:
        actual_yield = float(arr_4d[:, v["iz"], v["iy"], v["ix"]].sum())
        print(
            f"  Voxel {v['name']} ({v['shape']}) : ratio cible = {v['ratio']}, "
            f"yield reel = {actual_yield:.3f}, "
            f"moyenne theorique = {v['theory_mean']:.3f} MeV, "
            f"ecart-type theorique = {v['theory_std']:.3f} MeV"
        )


    # ecriture sur disque
    itk_img_4d = itk.image_from_array(arr_4d)
    itk_img_4d.SetSpacing([spacing_mm, spacing_mm, spacing_mm, 1.0])
    pg_image_path = str(paths.output / "pg_spectrum_4d_4voxels.nii.gz")
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
    phsp.attributes = ["KineticEnergy", "PrePosition", "TrackID"]
    phsp.steps_to_store = "first"
    phsp.output_filename = "test_pg_tle_4voxels_phsp.root"


    sim.g4_commands_after_init.append("/tracking/verbose 0")


    sim.run()
    print(stats)


    # ------------------------------------------------------------------
    # verification
    # ------------------------------------------------------------------
    import uproot


    hits = uproot.open(phsp.get_output_path_string())["PhaseSpace"]
    energies = hits["KineticEnergy"].array(library="np")
    x = hits["PrePosition_X"].array(library="np")
    y = hits["PrePosition_Y"].array(library="np")
    z = hits["PrePosition_Z"].array(library="np")
    track_id = hits["TrackID"].array(library="np")


    is_primary = track_id == 1
    n_secondary = int((~is_primary).sum())
    if n_secondary > 0:
        print(f"\n{n_secondary} particule(s) secondaire(s) exclue(s) de l'analyse")


    energies = energies[is_primary]
    x = x[is_primary]
    y = y[is_primary]
    z = z[is_primary]


    n_total = len(energies)
    print(f"\nNombre total de particules primaires enregistrees : {n_total}")


    voxel_half = spacing_mm / 2.0


    counts_per_voxel = {v["name"]: 0 for v in voxels}
    energies_per_voxel = {v["name"]: [] for v in voxels}
    positions_per_voxel = {v["name"]: {"x": [], "y": [], "z": []} for v in voxels}
    n_mismatched = 0


    # tolerance energie large pour englober tout le spectre de chaque
    # voxel, calculee a partir de l'etendue reelle du profil construit
    tol_E = {}
    for v in voxels:
        centers = Emin + (v["bins_range"] + 0.5) * bin_width
        tol_E[v["name"]] = (centers.min() - bin_width, centers.max() + bin_width)


    for i in range(n_total):
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
            n_mismatched += 1
            continue


        lo, hi = tol_E[matched_voxel["name"]]
        if lo <= energies[i] <= hi:
            counts_per_voxel[matched_voxel["name"]] += 1
            energies_per_voxel[matched_voxel["name"]].append(energies[i])
            positions_per_voxel[matched_voxel["name"]]["x"].append(x[i])
            positions_per_voxel[matched_voxel["name"]]["y"].append(y[i])
            positions_per_voxel[matched_voxel["name"]]["z"].append(z[i])
        else:
            n_mismatched += 1


    print(f"\nParticules incoherentes : {n_mismatched}")
    print(f"Repartition par voxel : {counts_per_voxel}")


    # verification des ratios de comptage relatifs entre voxels (base = Z)
    n_Z = counts_per_voxel["Z"]
    print("\nRatios de comptage observes (reference = Z) :")
    ratio_checks = []
    for v in voxels:
        if v["name"] == "Z":
            continue
        n_v = counts_per_voxel[v["name"]]
        ratio_observe = n_v / n_Z if n_Z > 0 else float("inf")
        ratio_attendu = v["ratio"] / voxels[-1]["ratio"]
        sigma_ratio = ratio_observe * np.sqrt(1.0 / n_v + 1.0 / n_Z)
        ecart_sigma = abs(ratio_observe - ratio_attendu) / sigma_ratio
        ratio_checks.append(ecart_sigma < 3.0)
        print(
            f"  {v['name']}/Z : observe = {ratio_observe:.2f}, "
            f"attendu = {ratio_attendu:.1f}, ecart = {ecart_sigma:.2f} sigma"
        )


    # verification de la position moyenne mesuree pour chaque voxel,
    # comparee a la position theorique (centre du voxel)
    print("\nPosition moyenne mesuree par voxel (mm) :")
    position_checks = []
    for v in voxels:
        px = np.array(positions_per_voxel[v["name"]]["x"])
        py = np.array(positions_per_voxel[v["name"]]["y"])
        pz = np.array(positions_per_voxel[v["name"]]["z"])
        mean_x, mean_y, mean_z = np.mean(px), np.mean(py), np.mean(pz)
        pos_ok = (
            abs(mean_x - v["expected_x"]) < 0.05
            and abs(mean_y - v["expected_y"]) < 0.05
            and abs(mean_z - v["expected_z"]) < 0.05
        )
        position_checks.append(pos_ok)
        print(
            f"  Voxel {v['name']} : attendu "
            f"({v['expected_x']:.2f}, {v['expected_y']:.2f}, {v['expected_z']:.2f}), "
            f"mesure moyenne ({mean_x:.2f}, {mean_y:.2f}, {mean_z:.2f})"
        )


    # verification de la moyenne et de l'ecart-type mesures pour chaque
    # voxel, compares aux valeurs theoriques calculees depuis le profil
    print("\nVerification moyenne / ecart-type energie par voxel :")
    shape_checks = []
    for v in voxels:
        e = np.array(energies_per_voxel[v["name"]])
        mean_mes = float(np.mean(e))
        std_mes = float(np.std(e))
        mean_ok = abs(mean_mes - v["theory_mean"]) < 0.03
        std_ok = abs(std_mes - v["theory_std"]) < 0.03
        shape_checks.append(mean_ok and std_ok)
        print(
            f"  Voxel {v['name']} ({v['shape']}) : "
            f"moyenne mesuree = {mean_mes:.3f} MeV (theorique {v['theory_mean']:.3f}), "
            f"ecart-type mesure = {std_mes:.3f} MeV (theorique {v['theory_std']:.3f})"
        )


    is_ok = n_mismatched == 0
    ratios_ok = all(ratio_checks)
    positions_ok = all(position_checks)
    shapes_ok = all(shape_checks)


    utility.print_test(is_ok, f"Coherence position/energie : {n_mismatched} incoherentes (attendu : 0)")
    utility.print_test(ratios_ok, "Ratios de comptage conformes aux rendements relatifs (< 3 sigma)")
    utility.print_test(positions_ok, "Positions moyennes conformes pour les 4 voxels (tolerance 0.05 mm)")
    utility.print_test(shapes_ok, "Moyenne et ecart-type d'energie conformes pour les 4 voxels")


    utility.test_ok(is_ok and ratios_ok and positions_ok and shapes_ok)


