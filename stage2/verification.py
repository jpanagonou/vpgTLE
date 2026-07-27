#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import numpy as np
import opengate as gate
from pathlib import Path


if __name__ == "__main__":
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)


    sim = gate.Simulation()


    sim.g4_verbose = False
    sim.visu = True
    sim.number_of_threads = 1
    sim.random_seed = 123456
    sim.output_dir = str(output_dir)


    m = gate.g4_units.m
    mm = gate.g4_units.mm
    MeV = gate.g4_units.MeV
    Bq = gate.g4_units.Bq


    # monde en VIDE (pas d'air), pour eliminer toute interaction parasite
    sim.world.size = [3 * m, 3 * m, 3 * m]
    sim.world.material = "G4_Galactic"


    # detecteur spherique de 1 m de rayon, fine coquille
    sphere = sim.add_volume("Sphere", "detector_sphere")
    sphere.rmin = 999 * mm
    sphere.rmax = 1000 * mm
    sphere.material = "G4_Galactic"
    sphere.color = [1, 0, 0, 0.3]


    # source ponctuelle, gamma monoenergetique 2 MeV, isotrope
    source = gate.sources.generic.GenericSource
    source = sim.add_source("GenericSource", "point_source")
    source.attached_to = sim.world.name
    source.particle = "gamma"
    source.energy.mono = 2 * MeV
    source.position.type = "point"
    source.direction.type = "iso"
    source.n = 100  # nombre fixe de particules


    # physique minimale (le vide n'interagit de toute facon pas)
    sim.physics_manager.physics_list_name = "G4EmStandardPhysics"
    sim.physics_manager.global_production_cuts.all = 1 * mm


    stats = sim.add_actor("SimulationStatisticsActor", "Stats")


    phsp = sim.add_actor("PhaseSpaceActor", "PhaseSpace")
    phsp.attached_to = sphere.name
    phsp.attributes = ["KineticEnergy", "PrePosition", "TrackID"]
    phsp.steps_to_store = "entering"
    phsp.output_filename = "point_source_2MeV_vide_phsp.root"


    sim.g4_commands_after_init.append("/tracking/verbose 0")


    sim.run()
    print(stats)


    # ------------------------------------------------------------------
    # verification : distribution plate attendue sur chaque axe
    # ------------------------------------------------------------------
    import uproot


    hits = uproot.open(phsp.get_output_path_string())["PhaseSpace"]
    energies = hits["KineticEnergy"].array(library="np")
    x = hits["PrePosition_X"].array(library="np")
    y = hits["PrePosition_Y"].array(library="np")
    z = hits["PrePosition_Z"].array(library="np")
    track_id = hits["TrackID"].array(library="np")


    is_primary = track_id == 1
    energies = energies[is_primary]
    x = x[is_primary]
    y = y[is_primary]
    z = z[is_primary]


    n = len(energies)
    R = 1000.0  # mm
    std_theorique = R / np.sqrt(3)


    print(f"\nNombre de gammas detectes : {n}")
    print(f"Energie : moyenne = {np.mean(energies):.4f} MeV (attendu 2.0000, "
          f"aucune interaction ne devrait la modifier dans le vide)")


    for axis_name, data in [("X", x), ("Y", y), ("Z", z)]:
        print(
            f"{axis_name} : moyenne = {np.mean(data):.3f} mm (attendu ~0), "
            f"std = {np.std(data):.3f} mm (attendu {std_theorique:.3f})"
        )


    # sauvegarde pour tracer l'histogramme si besoin
    np.savez(
        output_dir / "point_source_2MeV_vide.npz",
        energies=energies, x=x, y=y, z=z,
    )
    print(f"\nDonnees sauvegardees dans {output_dir / 'point_source_2MeV_vide.npz'}")
