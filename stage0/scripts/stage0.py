import opengate as gate
from pathlib import Path
import time
from opengate.actors.pgactors_utils import (
    reorganize_root_file,
    convert_secondaries_to_root
)

MeV = gate.g4_units.MeV
cm  = gate.g4_units.cm
m   = gate.g4_units.m
mm  = gate.g4_units.mm

t0 = time.time()


# ── Paramètres ────────────────────────────────────────────────────────────────
particle_type = "proton"  # "proton", "neutron"
material_name = "AllElements" # Oxygen
data_dir   = Path(__file__).parent.parent / "data"
output_dir = Path(__file__).parent.parent / "output"
output_dir.mkdir(exist_ok=True)

# ── Simulation ────────────────────────────────────────────────────────────────
sim = gate.Simulation()

sim.g4_verbose = True 
sim.g4_verbose_level = 1
sim.progress_bar = True
sim.physics_manager.physics_list_name = "QGSP_BIC_HP_EMY"
sim.number_of_threads = 1

sim.random_engine = "MersenneTwister"
sim.random_seed = 123456

sim.volume_manager.add_material_database(str(data_dir / "GateMaterials.db"))
sim.volume_manager.add_material_database(str(data_dir / "GateMaterialsElements.db"))

# ── Monde ─────────────────────────────────────────────────────────────────────
world          = sim.world
world.size     = [3 * m, 3 * m, 3 * m]
world.material = "G4_Galactic"

target          = sim.add_volume("Box", "target")
target.size     = [70 * cm, 70 * cm, 70 * cm]
target.material = material_name
target.mother   = "world"

# ── Step limiter ──────────────────────────────────────────────────────────────
sim.physics_manager.set_max_step_size("target", 1 * mm)
sim.physics_manager.set_user_limits_particles(["all"])

# ── Source ────────────────────────────────────────────────────────────────────
source                      = sim.add_source("GenericSource", "source")
source.particle             = "proton"  # proton, neutron
source.energy.mono          = 200* MeV
source.n                    = 1e5
source.direction.type       = "momentum"
source.direction.momentum   = [1, 0, 0]
source.position.type        = "disc"
source.position.radius      = 4 * mm
source.position.translation = [-40 * cm, 0, 0]

# ── Acteur ────────────────────────────────────────────────────────────────────
actor                      = sim.add_actor("PromptGammaStatisticActor", "pg_actor")
actor.attached_to          = "target"
actor.particle_type        = particle_type
actor.multi_element        = True
actor.material_name        = material_name
actor.particleNbBins       = 500
actor.particleMinEnergy    = 0 * MeV
actor.particleMaxEnergy    = 200 * MeV
actor.gammaNbBins          = 250
actor.gammaMinEnergy       = 0 * MeV
actor.gammaMaxEnergy       = 10 * MeV
actor.pg_output_filename   = str(output_dir /
    f"PGdb_{material_name}_{particle_type}")
actor.save_KE0_secondaries = False

# ── Lancement ─────────────────────────────────────────────────────────────────
sim.run()




# ── Post-traitement ───────────────────────────────────────────────────────────
root_file = str(output_dir / f"PGdb_{material_name}_{particle_type}.root")
reorganize_root_file(root_file)

t1 = time.time()
print(f"Temps de simulation : {t1 - t0:.2f}s")

if actor.save_KE0_secondaries:
    convert_secondaries_to_root(output_dir, material_name, particle_type)


print(f"Simulation terminée — PGdb_{material_name}_{particle_type}.root")



