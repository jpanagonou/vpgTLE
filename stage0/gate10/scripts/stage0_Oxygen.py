import opengate as gate
from pathlib import Path





MeV = gate.g4_units.MeV
cm  = gate.g4_units.cm
m  = gate.g4_units.m


sim = gate.Simulation()
sim.progress_bar = True
#sim.g4_verbose = True
#sim.g4_verbose_level = 1
sim.physics_manager.physics_list_name = "QGSP_BIC_HP"
sim.number_of_threads = 1


data_dir = Path(__file__).parent.parent / "data"
sim.volume_manager.add_material_database(str(data_dir / "GateMaterials.db"))
sim.volume_manager.add_material_database(str(data_dir / "GateMaterialsElements.db"))


# ── Monde ────────────────────────────────────────────────────────────
world = sim.world
world.size = [3* m, 3* m, 3* m]
world.material = "G4_Galactic"

target = sim.add_volume("Box", "target")
target.size = [70 * cm, 70 * cm, 70 * cm, ] 
target.material = "Titanium"
target.mother = "world"



# ── Source proton ─────────────────────────────────────────────────────
source = sim.add_source("GenericSource", "proton_source")
source.particle = "proton"
source.energy.mono = 200 * MeV
source.n = 1e6
source.direction.type = "momentum"
source.direction.momentum = [1, 0, 0]
source.position.type = "disc"
source.position.radius = 4 * gate.g4_units.mm
source.position.translation = [-40 *cm, 0, 0]


# ── Acteur ────────────────────────────────────────────────────────────
actor = sim.add_actor("PromptGammaStatisticActor", "pg_actor")
actor.attached_to        = "target"
actor.protonNbBins       = 500
actor.protonMinEnergy    = 0 * MeV
actor.protonMaxEnergy    = 200 * MeV
actor.gammaNbBins        = 250
actor.gammaMinEnergy     = 0 * MeV
actor.gammaMaxEnergy     = 10 * MeV
actor.pg_output_filename = f"../output/PGdb_Titanium"

# ── Lancement ─────────────────────────────────────────────────────────
sim.run()

print("Simulation terminée.")

