import opengate as gate
import sys
from pathlib import Path

element = sys.argv[1]  # reçoit le nom de l'élément

MeV = gate.g4_units.MeV
cm  = gate.g4_units.cm
mm  = gate.g4_units.mm
m   = gate.g4_units.m

data_dir   = Path(__file__).parent.parent / "data"
output_dir = Path(__file__).parent.parent / "output"
output_dir.mkdir(exist_ok=True)

sim = gate.Simulation()
sim.progress_bar = True
sim.physics_manager.physics_list_name = "QGSP_BIC_HP_EMY"
sim.number_of_threads = 1

sim.volume_manager.add_material_database(str(data_dir / "GateMaterials.db"))
sim.volume_manager.add_material_database(str(data_dir / "GateMaterialsElements.db"))

world = sim.world
world.size = [3 * m, 3 * m, 3 * m]
world.material = "G4_Galactic"

target = sim.add_volume("Box", "target")
target.size = [70 * cm, 70 * cm, 70 * cm]
target.material = element
target.mother = "world"

source = sim.add_source("GenericSource", "proton_source")
source.particle = "proton"
source.energy.mono = 200 * MeV
source.n = 1e6
source.direction.type = "momentum"
source.direction.momentum = [1, 0, 0]
source.position.type = "disc"
source.position.radius = 4 * mm
source.position.translation = [-40 * cm, 0, 0]

actor = sim.add_actor("PromptGammaStatisticActor", "pg_actor")
actor.attached_to        = "target"
actor.protonNbBins       = 500
actor.protonMinEnergy    = 0 * MeV
actor.protonMaxEnergy    = 200 * MeV
actor.gammaNbBins        = 250
actor.gammaMinEnergy     = 0 * MeV
actor.gammaMaxEnergy     = 10 * MeV
actor.pg_output_filename = str(output_dir / f"PGdb_{element}")

sim.run()
print(f"Simulation terminée pour {element}")

