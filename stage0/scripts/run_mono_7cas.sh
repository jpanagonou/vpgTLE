#!/bin/bash


# ── Paramètres ────────────────────────────────────────────────────────────────
N_TOTAL=1000000
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/stage0_run_mono_7cas.py"


declare -A CAS_ELEMENTS
declare -A CAS_FRACTIONS


# Cas 1 — 2 éléments
CAS_FRACTIONS["1_Oxygen"]="0.5238858952407169"
CAS_FRACTIONS["1_Carbon"]="0.4761141047592830"
CAS_ELEMENTS["1"]="Oxygen Carbon"


# Cas 2 — 5 éléments
CAS_FRACTIONS["2_Oxygen"]="0.2707013678876229"
CAS_FRACTIONS["2_Carbon"]="0.2460168151114442"
CAS_FRACTIONS["2_Hydrogen"]="0.1077848150126705"
CAS_FRACTIONS["2_Helium"]="0.1705737979586709"
CAS_FRACTIONS["2_Lithium"]="0.2049232040295915"
CAS_ELEMENTS["2"]="Oxygen Carbon Hydrogen Helium Lithium"


# Cas 3 — 10 éléments
CAS_FRACTIONS["3_Oxygen"]="0.1177369462011731"
CAS_FRACTIONS["3_Carbon"]="0.1070008206880746"
CAS_FRACTIONS["3_Hydrogen"]="0.0468791682342674"
CAS_FRACTIONS["3_Helium"]="0.0741881662080365"
CAS_FRACTIONS["3_Lithium"]="0.0891278549365142"
CAS_FRACTIONS["3_Beryllium"]="0.0972328839577161"
CAS_FRACTIONS["3_Boron"]="0.1033145278623951"
CAS_FRACTIONS["3_Nitrogen"]="0.1126381437182828"
CAS_FRACTIONS["3_Fluorine"]="0.1246738585359084"
CAS_FRACTIONS["3_Neon"]="0.1272076296576318"
CAS_ELEMENTS["3"]="Oxygen Carbon Hydrogen Helium Lithium Beryllium Boron Nitrogen Fluorine Neon"


# Cas 4 — 15 éléments
CAS_FRACTIONS["4_Oxygen"]="0.0693752981769286"
CAS_FRACTIONS["4_Carbon"]="0.0630491454035801"
CAS_FRACTIONS["4_Hydrogen"]="0.0276230731259301"
CAS_FRACTIONS["4_Helium"]="0.0437146224523937"
CAS_FRACTIONS["4_Lithium"]="0.0525176821006175"
CAS_FRACTIONS["4_Beryllium"]="0.0572934880240851"
CAS_FRACTIONS["4_Boron"]="0.0608770348452511"
CAS_FRACTIONS["4_Nitrogen"]="0.0663708806681598"
CAS_FRACTIONS["4_Fluorine"]="0.0734628032225171"
CAS_FRACTIONS["4_Neon"]="0.0749558020878122"
CAS_FRACTIONS["4_Sodium"]="0.0782848743671426"
CAS_FRACTIONS["4_Magnesium"]="0.0797498871158773"
CAS_FRACTIONS["4_Aluminium"]="0.0825744171606997"
CAS_FRACTIONS["4_Silicon"]="0.0836916470647425"
CAS_FRACTIONS["4_Phosphor"]="0.0864593441842627"
CAS_ELEMENTS["4"]="Oxygen Carbon Hydrogen Helium Lithium Beryllium Boron Nitrogen Fluorine Neon Sodium Magnesium Aluminium Silicon Phosphor"


# Cas 5 — 20 éléments
CAS_FRACTIONS["5_Oxygen"]="0.0475274150028906"
CAS_FRACTIONS["5_Carbon"]="0.0431935137998452"
CAS_FRACTIONS["5_Hydrogen"]="0.0189239296206423"
CAS_FRACTIONS["5_Helium"]="0.0299478785329464"
CAS_FRACTIONS["5_Lithium"]="0.0359786514476705"
CAS_FRACTIONS["5_Beryllium"]="0.0392504458191921"
CAS_FRACTIONS["5_Boron"]="0.0417054510073138"
CAS_FRACTIONS["5_Nitrogen"]="0.0454691579354105"
CAS_FRACTIONS["5_Fluorine"]="0.0503276703348771"
CAS_FRACTIONS["5_Neon"]="0.0513504893862454"
CAS_FRACTIONS["5_Sodium"]="0.0536311599412150"
CAS_FRACTIONS["5_Magnesium"]="0.0546348063502879"
CAS_FRACTIONS["5_Aluminium"]="0.0565698266695671"
CAS_FRACTIONS["5_Silicon"]="0.0573352150815587"
CAS_FRACTIONS["5_Phosphor"]="0.0592313004759057"
CAS_FRACTIONS["5_Sulfur"]="0.0599219303821312"
CAS_FRACTIONS["5_Chlorine"]="0.0619597437258353"
CAS_FRACTIONS["5_Argon"]="0.0644777381434128"
CAS_FRACTIONS["5_Potassium"]="0.0640160756182466"
CAS_FRACTIONS["5_Calcium"]="0.0645476007248056"
CAS_ELEMENTS["5"]="Oxygen Carbon Hydrogen Helium Lithium Beryllium Boron Nitrogen Fluorine Neon Sodium Magnesium Aluminium Silicon Phosphor Sulfur Chlorine Argon Potassium Calcium"


# Cas 6 — 25 éléments
CAS_FRACTIONS["6_Oxygen"]="0.0338959352822926"
CAS_FRACTIONS["6_Carbon"]="0.0308050532158191"
CAS_FRACTIONS["6_Hydrogen"]="0.0134963008964182"
CAS_FRACTIONS["6_Helium"]="0.0213584381253004"
CAS_FRACTIONS["6_Lithium"]="0.0256595070642961"
CAS_FRACTIONS["6_Beryllium"]="0.0279929083289624"
CAS_FRACTIONS["6_Boron"]="0.0297437861532498"
CAS_FRACTIONS["6_Nitrogen"]="0.0324280130662540"
CAS_FRACTIONS["6_Fluorine"]="0.0358930410264433"
CAS_FRACTIONS["6_Neon"]="0.0366225022935577"
CAS_FRACTIONS["6_Sodium"]="0.0382490469210485"
CAS_FRACTIONS["6_Magnesium"]="0.0389648345086161"
CAS_FRACTIONS["6_Aluminium"]="0.0403448658759483"
CAS_FRACTIONS["6_Silicon"]="0.0408907309535484"
CAS_FRACTIONS["6_Phosphor"]="0.0422429944379516"
CAS_FRACTIONS["6_Sulfur"]="0.0427355427199066"
CAS_FRACTIONS["6_Chlorine"]="0.0441888847375902"
CAS_FRACTIONS["6_Argon"]="0.0459846856624710"
CAS_FRACTIONS["6_Potassium"]="0.0456554339437662"
CAS_FRACTIONS["6_Calcium"]="0.0460345107484219"
CAS_FRACTIONS["6_Titanium"]="0.0488412043083685"
CAS_FRACTIONS["6_Copper"]="0.0536349497613494"
CAS_FRACTIONS["6_Zinc"]="0.0541931930169159"
CAS_FRACTIONS["6_Silver"]="0.0640331946774320"
CAS_FRACTIONS["6_Tin"]="0.0661104422740718"
CAS_ELEMENTS["6"]="Oxygen Carbon Hydrogen Helium Lithium Beryllium Boron Nitrogen Fluorine Neon Sodium Magnesium Aluminium Silicon Phosphor Sulfur Chlorine Argon Potassium Calcium Titanium Copper Zinc Silver Tin"


# ── Boucle sur les 6 cas ──────────────────────────────────────────────────────
echo "============================================="
echo "  Boucle mono Gate 10 — 6 cas"
echo "  N_total = $N_TOTAL particules"
echo "============================================="


T_GLOBAL_START=$(date +%s%N)


for CAS in 1 2 3 4 5 6; do
   ELEMS=${CAS_ELEMENTS[$CAS]}
   N_ELEMS=$(echo $ELEMS | wc -w)


   echo ""
   echo "===== Cas $CAS ($N_ELEMS éléments) ====="


   T_CAS_START=$(date +%s%N)
   TOTAL_PARTICLES_CAS=0


   for ELEM in $ELEMS; do
       F=${CAS_FRACTIONS["${CAS}_${ELEM}"]}
       N=$(python3 -c "print(int(${F} * ${N_TOTAL}))")
       TOTAL_PARTICLES_CAS=$((TOTAL_PARTICLES_CAS + N))


       echo "→ $ELEM : $N particules"
       TMPFILE=$(mktemp)
       python3 "$PYTHON_SCRIPT" "$ELEM" "$N" "$CAS" | tee "$TMPFILE"
       T=$(grep "TEMPS_${ELEM}" "$TMPFILE" | cut -d= -f2)
       rm "$TMPFILE"


       #OUTPUT=$(python3 "$PYTHON_SCRIPT" "$ELEM" "$N" "$CAS")
       #echo "$OUTPUT"
   done


   T_CAS_END=$(date +%s%N)
   T_CAS=$(python3 -c "print(round((${T_CAS_END} - ${T_CAS_START}) / 1e9, 2))")


   echo "---------------------------------------------"
   echo "  Cas $CAS ($N_ELEMS elems) : $TOTAL_PARTICLES_CAS particules — Temps : ${T_CAS}s"
   echo "---------------------------------------------"
done


T_GLOBAL_END=$(date +%s%N)
T_GLOBAL=$(python3 -c "print(round((${T_GLOBAL_END} - ${T_GLOBAL_START}) / 1e9, 2))")


echo ""
echo "============================================="
echo "  Temps total global : ${T_GLOBAL}s"
echo "============================================="



