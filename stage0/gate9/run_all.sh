#!/bin/bash

ELEMENTS=(
    "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron",
    "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon",
    "Sodium", "Magnesium", "Aluminium", "Silicon", "Phosphor",
    "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium",
    "Scandium", "Titanium", "Vandium", "Chromium", "Manganese",
    "Iron", "Cobalt", "Copper", "Zinc", "Silver", "Tin"
)




for ELEMENT in "${ELEMENTS[@]}"; do
    echo "=== Simulation pour $ELEMENT ==="
    sed -i "s|/control/alias ELEMENT .*|/control/alias ELEMENT $ELEMENT|" mac/main_all.mac
    Gate mac/main_all.mac
    echo "=== Terminé pour $ELEMENT ==="
done
