import opengate as gate
from pathlib import Path
import time
import sys


MeV = gate.g4_units.MeV
cm  = gate.g4_units.cm
m   = gate.g4_units.m
mm  = gate.g4_units.mm


# ── Arguments ─────────────────────────────────────────────────────────────────
if len(sys.argv) != 4:
   print("Usage: python stage0_run_mono.py <material_name> <n_particles> <cas_number>")
   sys.exit(1)


material_name = sys.argv[1]
n_particles   = int(sys.argv[2])
cas_number    = sys.argv[3]
particle_type = "proton"


data_dir   = Path(__file__).parent.parent / "data"
output_dir = Path(__file__).parent.parent / "output"
output_dir.mkdir(exist_ok=True)


print(f"→ Simulation {material_name} ({n_particles} particules) — Cas {cas_number}...")
t0 = time.time()


# ── Simulation ────────────────────────────────────────────────────────────────
sim = gate.Simulation()


#sim.g4_verbose = True
#sim.g4_verbose_level = 1


sim.progress_bar = False
sim.physics_manager.physics_list_name = "QGSP_BIC_HP_EMY"
sim.number_of_threads = 1
sim.random_engine = "MersenneTwister"
sim.random_seed = 123456


sim.volume_manager.add_material_database(str(data_dir / "GateMaterials.db"))
sim.volume_manager.add_material_database(str(data_dir / "GateMaterialsElements.db"))


# Monde
world          = sim.world
world.size     = [3 * m, 3 * m, 3 * m]
world.material = "G4_Galactic"


target          = sim.add_volume("Box", "target")
target.size     = [70 * cm, 70 * cm, 70 * cm]
target.material = material_name
target.mother   = "world"


# Step limiter
sim.physics_manager.set_max_step_size("target", 1 * mm)
sim.physics_manager.set_user_limits_particles(["all"])


# Source
source                      = sim.add_source("GenericSource", "source")
source.particle             = "proton"
source.energy.mono          = 200 * MeV
source.n                    = n_particles
source.direction.type       = "momentum"
source.direction.momentum   = [1, 0, 0]
source.position.type        = "disc"
source.position.radius      = 4 * mm
source.position.translation = [-40 * cm, 0, 0]


# Acteur
actor                      = sim.add_actor("PromptGammaStatisticActor", "pg_actor")
actor.attached_to          = "target"
actor.particle_type        = particle_type
actor.multi_element        = False
actor.material_name        = material_name
actor.particleNbBins       = 500
actor.particleMinEnergy    = 0 * MeV
actor.particleMaxEnergy    = 200 * MeV
actor.gammaNbBins          = 250
actor.gammaMinEnergy       = 0 * MeV
actor.gammaMaxEnergy       = 10 * MeV
actor.pg_output_filename   = str(output_dir /
   f"PGdb_cas{cas_number}_{material_name}_{particle_type}")
actor.save_KE0_secondaries = False


sim.run()


t1 = time.time()
print(f"TEMPS_{material_name}={t1 - t0:.2f}")



