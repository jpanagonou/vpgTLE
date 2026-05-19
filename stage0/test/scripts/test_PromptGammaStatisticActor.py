import subprocess
import sys
import os
from pathlib import Path

tests = [
    ("proton",  "proton",    "Oxygen",      "False", "100", "Mono proton"),
    ("neutron", "neutron",   "Oxygen",      "False", "100", "Mono neutron"),
    ("helium",  "alpha",     "Oxygen",      "False", "100", "Mono helium"),
    ("carbon",  "ion 6 12",  "Oxygen",      "False", "100", "Mono carbon"),
    ("proton",  "proton",    "AllElements", "True",  "10",  "Multi proton"),
]

passed = 0
failed = 0

for particle_type, source_particle, material_name, \
        multi_element, n_particles, test_name in tests:

    print(f"\n{'='*50}")
    print(f"Test : {test_name}")
    print(f"{'='*50}")

    result = subprocess.run(
        [sys.executable,
         str(Path(__file__).parent / "run_single_test.py"),
         particle_type, source_particle, material_name,
         multi_element, n_particles, test_name],
        capture_output=True, text=True
    )

    print(result.stdout)
    if result.stderr:
        print(result.stderr)

    if result.returncode == 0:
        passed += 1
    else:
        print(f"  ❌ FAILED : {test_name}")
        failed += 1

print(f"\n{'='*50}")
print(f"Results : {passed} passed, {failed} failed")
print(f"{'='*50}")

if failed > 0:
    sys.exit(1)

