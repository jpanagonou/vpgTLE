#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#=======================================================================
# Validation de stage 1 Gate 10 vs Gate 9
# Simulation de protons sur boite d'Oxygen non voxélisée
# Pour comparer avec la simulation Gate 9
#=======================================================================

import opengate as gate
from opengate.tests import utility
import uproot
from box import Box
import pathlib
import numpy as np

def simulation(
    paths, File_name, job_id, number_of_particles, visu, verbose, actor, Erange
):
    # create the simulation
    sim = gate.Simulation()
    # main options
    sim.visu = visu
    sim.g4_verbose = verbose
    sim.random_seed = "auto"
    sim.random_engine = "MersenneTwister"
    sim.output_dir = paths.output
    sim.number_of_threads = 1
    sim.progress_bar = True

    # units
    mm = gate.g4_units.mm
    cm = gate.g4_units.cm
    m  = gate.g4_units.m
    MeV = gate.g4_units.MeV
    ns  = gate.g4_units.ns
    gcm3 = gate.g4_units.g_cm3

    # ── World ──────────────────────────────────────────────────────────
    world = sim.world
    world.size = [3 * m, 3 * m, 3 * m]
    world.material = "G4_Galactic"

    # ── Oxygen box (même géométrie que Gate 9) ─────────────────────────
    # Je lui passe l'elements Oxygen avec densité 1.0 g/cm3 comme dans Gate 9 depuis GateMaterials.db
    
    sim.volume_manager.add_material_database(str(paths.data /"OxygenMaterial.db"))
    
    
    box = sim.add_volume("Box", vol_name)
    box.size = [40 * mm, 150 * mm, 40 * mm]
    box.translation = [0, 0, 0]
    box.material = "Oxygen"
    

    # ── Physics ────────────────────────────────────────────────────────
    sim.physics_manager.physics_list_name = "QGSP_BIC_HP_EMY"

    sim.physics_manager.apply_cuts = False
    sim.physics_manager.global_production_cuts.gamma    = 0.1 * mm
    sim.physics_manager.global_production_cuts.electron = 0.1 * mm
    sim.physics_manager.global_production_cuts.positron = 0.1 * mm
    sim.physics_manager.global_production_cuts.proton   = 0.1 * mm

    sim.physics_manager.set_max_step_size(vol_name, 0.1 * mm)
    sim.physics_manager.set_user_limits_particles(["proton"])

    # ── Source ─────────────────────────────────────────────────────────
    source = sim.add_source("GenericSource", "DEFAULT")
    source.energy.mono          = Erange * MeV
    source.particle             = "proton"
    source.position.type        = "point"
    source.position.radius      = 0 * mm   # faisceau ponctuel comme Gate 9
    source.position.translation = [0 * mm, -300 * mm, 0 * mm]
    source.n                    = number_of_particles
    source.direction.type       = "momentum"
    source.direction.momentum   = [0, 1, 0]

    # ── Base de données Oxygen ─────────────────────────────────────────
    # LOOKHERE : base de données Oxygen pour la validation Gate 9 vs Gate 10
    with uproot.open(paths.data / "db-Oxygen_merged.root") as root_file:
        histo  = root_file["Oxygen"]["Ep"].to_hist()
        vect_p = histo.to_numpy()[0]
    # Pas de neutrons dans cette validation (base de données Oxygen uniquement)
    vect_n = np.zeros(len(vect_p))  # vecteur nul pour les neutrons
    
    
   

    
    

    # ── Acteur VPG TLE ─────────────────────────────────────────────────
    vpg_tle = sim.add_actor(actor, actor_name)
    vpg_tle.attached_to      = vol_name
    vpg_tle.output_filename  = paths.output / f"{File_name}.nii.gz"
    vpg_tle.size             = [1, 1, 1]     # 1 seul voxel = boite entière
    vpg_tle.spacing          = [40, 150, 40] # même taille que la boite
    
    vpg_tle.timebins         = 250
    vpg_tle.timerange        = 5 * ns
    vpg_tle.energybins       = 250
    vpg_tle.energyrange      = Erange * MeV
    vpg_tle.prot_E.active    = True
    vpg_tle.neutr_E.active   = False   # pas de neutrons pour cette validation
    vpg_tle.prot_tof.active  = True
    vpg_tle.neutr_tof.active = False
    vpg_tle.weight           = True
    vpg_tle.vect_p           = vect_p
    vpg_tle.vect_n           = vect_n

    # ── Statistics actor ───────────────────────────────────────────────
    stats = sim.add_actor("SimulationStatisticsActor", "stats")
    stats.track_types_flag  = True
    stats.output_filename   = paths.output / f"stat_{job_id}_{File_name}.txt"

    return sim


# ── Paramètres de simulation ───────────────────────────────────────────
output             = "stage1a_oxygen"
File_name          = "vpg_oxygen"
actor_name         = "vpg_tle"
vol_name           = "oxygenbox"
MJ                 = False
number_of_particles = 1e6          # même que Gate 9
actor              = "VoxelizedPromptGammaTLEActor"
Erange             = 200           # même énergie que Gate 9 (200 MeV)

paths = Box()
paths.current = pathlib.Path().resolve().parent
paths.data    = (paths.current / "data").resolve()
paths.output  = (paths.current / "output" / output).resolve()

if __name__ == "__main__":
    sim = simulation(
        paths, File_name, 0, number_of_particles, False, False, actor, Erange
    )
    sim.run()
