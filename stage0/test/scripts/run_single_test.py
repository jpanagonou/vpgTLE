import opengate as gate
from pathlib import Path
from opengate.actors.pgactors_utils import reorganize_root_file
import uproot
import os
import sys

MeV = gate.g4_units.MeV
cm  = gate.g4_units.cm
m   = gate.g4_units.m
mm  = gate.g4_units.mm

particle_type   = sys.argv[1]
source_particle = sys.argv[2]
material_name   = sys.argv[3]
multi_element   = sys.argv[4] == "True"
n_particles     = int(sys.argv[5])
test_name       = sys.argv[6]

output_dir = Path(__file__).parent / "test_output"
output_dir.mkdir(exist_ok=True)

data_dir = Path(__file__).parent.parent / "data"

# ── Simulation ────────────────────────────────────────────────────────────────
sim = gate.Simulation()
sim.progress_bar = False
sim.physics_manager.physics_list_name = "QGSP_BIC_HP_EMY"
sim.volume_manager.add_material_database(str(data_dir / "GateMaterials.db"))
sim.volume_manager.add_material_database(
    str(data_dir / "GateMaterialsElements.db"))

world          = sim.world
world.size     = [3 * m, 3 * m, 3 * m]
world.material = "G4_Galactic"

target          = sim.add_volume("Box", "target")
target.size     = [70 * cm, 70 * cm, 70 * cm]
target.material = material_name

sim.physics_manager.set_max_step_size("target", 1 * mm)
sim.physics_manager.set_user_limits_particles(["all"])

source                      = sim.add_source("GenericSource", "source")
source.particle             = source_particle
source.energy.mono          = 200 * MeV
source.n                    = n_particles
source.direction.type       = "momentum"
source.direction.momentum   = [1, 0, 0]
source.position.type        = "disc"
source.position.radius      = 4 * mm
source.position.translation = [-40 * cm, 0, 0]

actor                      = sim.add_actor(
    "PromptGammaStatisticActor", "pg_actor")
actor.attached_to          = "target"
actor.particle_type        = particle_type
actor.multi_element        = multi_element
actor.material_name        = material_name
actor.particleNbBins       = 100
actor.particleMinEnergy    = 0 * MeV
actor.particleMaxEnergy    = 200 * MeV
actor.gammaNbBins          = 50
actor.gammaMinEnergy       = 0 * MeV
actor.gammaMaxEnergy       = 10 * MeV
actor.pg_output_filename   = str(output_dir /
    f"PGdb_{material_name}_{particle_type}")
actor.save_KE0_secondaries = False

sim.run()

# ── Post-processing ───────────────────────────────────────────────────────────
root_file = str(output_dir /
    f"PGdb_{material_name}_{particle_type}.root")

assert os.path.exists(root_file), \
    f"FAIL : ROOT file not found : {root_file}"
print(f"  ✅ ROOT file created")

reorganize_root_file(root_file)

# ── Check histograms ──────────────────────────────────────────────────────────
with uproot.open(root_file) as f:
    keys = [k.split(";")[0] for k in f.keys()]

    # Check GammaZ
    gammaz_keys = [k for k in keys if "GammaZ" in k]
    assert len(gammaz_keys) > 0, "FAIL : No GammaZ histogram found"
    print(f"  ✅ GammaZ histogram found : {gammaz_keys[0]}")

    # Check 2D spectrum
    epepg_map = {
        "proton"  : "EpEpg",
        "neutron" : "EnEpg",
        "helium"  : "EHeEpg",
        "carbon"  : "ECEpg"
    }
    expected = epepg_map[particle_type]
    epepg_keys = [k for k in keys if expected in k]
    assert len(epepg_keys) > 0, f"FAIL : No {expected} histogram found"
    print(f"  ✅ {expected} histogram found : {epepg_keys[0]}")

    # Check NrPG
    nrpg_keys = [k for k in keys if "NrPG" in k]
    assert len(nrpg_keys) > 0, "FAIL : No NrPG histogram found"
    print(f"  ✅ NrPG histogram found : {nrpg_keys[0]}")

    # Check Kapa inelastique
    kapa_keys = [k for k in keys if "Kapa" in k]
    assert len(kapa_keys) > 0, \
        "FAIL : No Kapa inelastique histogram found"
    print(f"  ✅ Kapa inelastique histogram found : {kapa_keys[0]}")

    # Check energy spectrum
    ep_map = {
        "proton"  : "Ep",
        "neutron" : "En",
        "helium"  : "EHe",
        "carbon"  : "EC"
    }
    expected_ep = ep_map[particle_type]
    ep_keys = [k for k in keys if k.endswith(expected_ep)]
    assert len(ep_keys) > 0, \
        f"FAIL : No {expected_ep} histogram found"
    print(f"  ✅ {expected_ep} histogram found : {ep_keys[0]}")

    # Check Inelastic
    inel_map = {
        "proton"  : "EpInelastic",
        "neutron" : "EnInelastic",
        "helium"  : "EHeInelastic",
        "carbon"  : "ECInelastic"
    }
    expected_inel = inel_map[particle_type]
    inel_keys = [k for k in keys if expected_inel in k
                 and "ProducedGamma" not in k]
    assert len(inel_keys) > 0, \
        f"FAIL : No {expected_inel} histogram found"
    print(f"  ✅ {expected_inel} histogram found : {inel_keys[0]}")

    # Check InelasticProducedGamma
    inel_pg_map = {
        "proton"  : "EpInelasticProducedGamma",
        "neutron" : "EnInelasticProducedGamma",
        "helium"  : "EHeInelasticProducedGamma",
        "carbon"  : "ECInelasticProducedGamma"
    }
    expected_inel_pg = inel_pg_map[particle_type]
    inel_pg_keys = [k for k in keys if expected_inel_pg in k]
    assert len(inel_pg_keys) > 0, \
        f"FAIL : No {expected_inel_pg} histogram found"
    print(f"  ✅ {expected_inel_pg} histogram found : {inel_pg_keys[0]}")

print(f"  ✅ Test {test_name} PASSED")
sys.exit(0)

