import subprocess
import sys

elements = [
    "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron",
    "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon",
    "Sodium", "Magnesium", "Aluminium", "Silicon", "Phosphor",
    "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium",
    "Scandium", "Titanium", "Vandium", "Chromium", "Manganese",
    "Iron", "Cobalt", "Copper", "Zinc", "Silver", "Tin"
]



for element in elements:
    print(f"\n=== Lancement simulation pour {element} ===")
    subprocess.run([sys.executable, "stage0_one_element.py", element], check=True)

