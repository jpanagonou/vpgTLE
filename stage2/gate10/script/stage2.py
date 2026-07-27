#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Comparaison Gate 9 / Gate 10 - cote Gate 10.


Utilise la vraie image PG (.nii.gz), avec un detecteur spherique de 1 m
de rayon englobant la source, pour capturer la quasi-totalite des gammas
emis quelle que soit leur direction (angle solide ~4 pi), contrairement
a un petit volume localise qui ne capture qu'une fraction infime.


Objectif : comparer le spectre en energie detecte au spectre theorique
global de la source (somme de tous les voxels de l'image), et comparer
ensuite avec le meme calcul effectue cote Gate 9.
"""


import numpy as np
import itk
import opengate as gate
from pathlib import Path


if __name__ == "__main__":
    data_dir = Path(__file__).parent.parent / "data"
    output_dir = Path(__file__).parent.parent / "output"
    output_dir.mkdir(exist_ok=True)


    # ------------------------------------------------------------------
    # chemin de la vraie image PG
    # ------------------------------------------------------------------
    pg_image_path = str(data_dir / "VPG_vpg_jeanpaul_gamma_e.nii.gz")


    # lecture rapide pour connaitre les parametres de l'image (Emin/Emax/
    # nbBins doivent correspondre exactement a la grille d'energie reelle)
    img = itk.imread(pg_image_path)
    arr_4d = itk.array_from_image(img)  # [bins, Z, Y, X]
    nbBins_img = arr_4d.shape[0]
    print(f"Image PG : forme = {arr_4d.shape}, nbBins = {nbBins_img}")



    Emin, Emax = 0.0, 10.0
    nbBins = nbBins_img


    # ------------------------------------------------------------------
    # simulation
    # ------------------------------------------------------------------
    sim = gate.Simulation()


    sim.g4_verbose = False
    sim.visu = False
    sim.number_of_threads = 1
    sim.random_seed = 123456
    sim.output_dir = str(output_dir)


    m = gate.g4_units.m
    mm = gate.g4_units.mm
    


    sim.volume_manager.add_material_database(str(data_dir / "GateMaterials.db"))
    sim.world.size = [3 * m, 3 * m, 3 * m]
    sim.world.material = "G4_AIR"


    # detecteur spherique de 1 m de rayon, fine coquille, englobant
    # toute la source pour capturer la quasi-totalite des gammas emis
    # quelle que soit leur direction (angle solide ~4 pi)
    sphere = sim.add_volume("Sphere", "detector_sphere")
    sphere.rmin = 999 * mm
    sphere.rmax = 1000 * mm  # coquille fine de 1 mm d'epaisseur
    sphere.material = "G4_AIR"
    sphere.color = [1, 0, 0, 0.3]


    source = sim.add_source("VoxelizedPromptGammaTLESource", "pg_source")
    source.attached_to = sim.world.name
    source.particle = "gamma"
    source.n_protons = 1e5 # a adapter selon selon le stage 1
    source.image = pg_image_path
    source.Emin = Emin
    source.Emax = Emax
    source.nbBins = nbBins
    source.direction.type = "iso"


    sim.physics_manager.physics_list_name = "G4EmStandardPhysics"
    sim.physics_manager.global_production_cuts.all = 1 * mm


    stats = sim.add_actor("SimulationStatisticsActor", "Stats")


    phsp = sim.add_actor("PhaseSpaceActor", "PhaseSpace")
    phsp.attached_to = sphere.name
    phsp.attributes = ["KineticEnergy", "PrePosition", "TrackID"]
    phsp.steps_to_store = "entering"
    phsp.output_filename = "gate10_realimage_phsp.root"


    sim.g4_commands_after_init.append("/tracking/verbose 0")


    sim.run()
    print(stats)


    # ------------------------------------------------------------------
    # extraction des donnees detectees (energie et position)
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
    print(f"\n{n_secondary} particule(s) secondaire(s) exclue(s)")


    energies = energies[is_primary]
    x = x[is_primary]
    y = y[is_primary]
    z = z[is_primary]


    n_detected = len(energies)
    print(f"Nombre de gammas primaires detectes par la sphere : {n_detected}")
    print(f"Energie moyenne detectee : {np.mean(energies):.3f} MeV (std {np.std(energies):.3f})")
    print(f"Position X detectee : moyenne {np.mean(x):.3f} mm (std {np.std(x):.3f})")
    print(f"Position Y detectee : moyenne {np.mean(y):.3f} mm (std {np.std(y):.3f})")
    print(f"Position Z detectee : moyenne {np.mean(z):.3f} mm (std {np.std(z):.3f})")


    # sauvegarde des donnees detectees pour comparaison ulterieure
    np.savez(
        output_dir / "gate10_spectrum_comparison.npz",
        energies=energies,
        x=x,
        y=y,
        z=z,
        n_detected=n_detected,
    )
    print(f"\nDonnees sauvegardees dans {output_dir / 'gate10_spectrum_comparison.npz'}")
