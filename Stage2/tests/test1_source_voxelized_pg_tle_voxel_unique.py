#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Test 1 : un seul voxel: minimal fonctionnel de VoxelizedPromptGammaTLESource.


Principe : on construit une image 4D synthetique [nbBins, Z, Y, X] avec
un seul voxel actif (au centre de l'image), dont le spectre en energie
est concentre etroitement autour d'une valeur connue (E_test).


Si le tirage position <-> energie par voxel est correct :
 - toutes les particules generees doivent naitre dans ce voxel unique
 - toutes les energies tirees doivent etre proches de E_test
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
    nZ, nY, nX = 5, 5, 5          # petite image 5x5x5 voxels
    nbBins = 250
    Emin, Emax = 0.0, 10.0        # MeV
    E_test = 2.0                  # energie attendue pour le voxel actif
    spacing_mm = 2.0              # mm


    arr_4d = np.zeros((nbBins, nZ, nY, nX), dtype=np.float32)


    # voxel actif unique, au centre de l'image
    iz_active, iy_active, ix_active = nZ // 2, nY // 2, nX // 2


    # bin correspondant a E_test
    bin_width = (Emax - Emin) / nbBins
    ib_test = int((E_test - Emin) / bin_width)


    # spectre etroit : 3 bins non nuls autour de ib_test (evite les effets
    # de bord d'interpolation, garde un pic net et verifiable)
    arr_4d[ib_test - 1 : ib_test + 2, iz_active, iy_active, ix_active] = [
        0.2,
        1.0,
        0.2,
    ]


    # ecriture sur disque 
    itk_img_4d = itk.image_from_array(arr_4d)
    itk_img_4d.SetSpacing([spacing_mm, spacing_mm, spacing_mm, 1.0])
    pg_image_path = str(paths.output / "pg_spectrum_4d.nii.gz")
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
    cm = gate.g4_units.cm
    m = gate.g4_units.m
    MeV = gate.g4_units.MeV


    sim.volume_manager.add_material_database(paths.data / "GateMaterials.db")
    sim.world.size = [1 * m, 1 * m, 1 * m]
    sim.world.material = "G4_AIR"


    # volume englobant toute l'image (necessaire pour attacher la source
    # ET le PhaseSpaceActor ; le volume doit etre assez grand pour contenir
    # toute l'etendue physique de l'image 4D)
    half_extent = (max(nX, nY, nZ) * spacing_mm) / 2.0 + 5 * mm
    box = sim.add_volume("Box", "box_pg")
    box.size = [2 * half_extent, 2 * half_extent, 2 * half_extent]
    box.material = "G4_AIR"


    # source voxelisee PG, attachee a ce volume
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


    # PhaseSpaceActor : capture le premier pas de chaque particule primaire.
    # PrePosition du premier pas = position de naissance exacte (aucun
    # deplacement avant), KineticEnergy = energie initiale tiree par la CDF.
    # TrackID est ajouté pour pouvoir exclure d'eventuelles particules
    # secondaires (ex : electron Compton) dont le premier pas serait
    # aussi enregistre mais ne provient pas directement du tirage source.
    phsp = sim.add_actor("PhaseSpaceActor", "PhaseSpace")
    phsp.attached_to = box.name
    phsp.attributes = ["KineticEnergy", "PrePosition", "TrackID"]
    phsp.steps_to_store = "first"
    phsp.output_filename = "test_pg_tle_phsp.root"


    # verbose tracking off
    sim.g4_commands_after_init.append("/tracking/verbose 0")


    sim.run()


    print(stats)


    # ------------------------------------------------------------------
    # verification : position et energie de chaque particule generee
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


    print(f"\nNombre de particules primaires enregistrées : {len(energies)}")


    # verification energie : toutes proches de E_test (tolerance = largeur
    # du spectre construit, quelques bins autour de ib_test)
    tol_E = 3 * bin_width
    e_ok = bool(((energies - E_test).__abs__() < tol_E).all())
    utility.print_test(
        e_ok,
        f"Énergie : attendu ~{E_test} MeV (+/- {tol_E:.3f}), "
        f"obtenu min={energies.min():.3f} max={energies.max():.3f} "
        f"mean={energies.mean():.3f}",
    )


    # verification position : toutes dans le voxel actif
    # (coordonnees centrees sur l'image, voxel [ix,iy,iz] au centre ici)
    voxel_half = spacing_mm / 2.0
    # position attendue du centre du voxel actif (image centree en 0)
    expected_x = (ix_active - nX / 2.0 + 0.5) * spacing_mm
    expected_y = (iy_active - nY / 2.0 + 0.5) * spacing_mm
    expected_z = (iz_active - nZ / 2.0 + 0.5) * spacing_mm


    pos_ok = bool(
        (abs(x - expected_x) <= voxel_half).all()
        and (abs(y - expected_y) <= voxel_half).all()
        and (abs(z - expected_z) <= voxel_half).all()
    )
    utility.print_test(
        pos_ok,
        f"Position : attendu centre voxel ({expected_x:.2f}, "
        f"{expected_y:.2f}, {expected_z:.2f}) mm +/- {voxel_half:.2f} mm, "
        f"obtenu x=[{x.min():.2f},{x.max():.2f}] "
        f"y=[{y.min():.2f},{y.max():.2f}] "
        f"z=[{z.min():.2f},{z.max():.2f}]",
    )


    utility.test_ok(e_ok and pos_ok)
