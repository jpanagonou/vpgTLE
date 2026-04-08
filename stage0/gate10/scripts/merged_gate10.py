import uproot
import numpy as np
from pathlib import Path

output_dir = Path("../output")

symbols = { "Hydrogen": "H", "Helium": "He", "Lithium": "Li", "Beryllium": "Be",
    "Boron": "B", "Carbon": "C", "Nitrogen": "N", "Oxygen": "O",
    "Fluorine": "F", "Neon": "Ne", "Sodium": "Na", "Magnesium": "Mg",
    "Aluminium": "Al", "Silicon": "Si", "Phosphor": "P", "Sulfur": "S",
    "Chlorine": "Cl", "Argon": "Ar", "Potassium": "K", "Calcium": "Ca",
    "Scandium": "Sc", "Titanium": "Ti", "Vandium": "V", "Chromium": "Cr",
    "Manganese": "Mn", "Iron": "Fe", "Cobalt": "Co", "Copper": "Cu",
    "Zinc": "Zn", "Silver": "Ag", "Tin": "Sn" 
}

with uproot.recreate(output_dir / "PGdb_gate10.root") as out:
    for element, symbol in symbols.items():
        fname = output_dir / f"PGdb_{element}.root"
        if fname.exists():
            f = uproot.open(fname)
            for key in f.keys():
                name = key.replace(";1", "")
                h = f[key]
                out[f"{symbol}/{name}"] = h.to_numpy()
            print(f"Ajouté : {symbol}")

print("Terminé → PGdb_gate10.root")

