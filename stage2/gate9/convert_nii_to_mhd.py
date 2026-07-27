#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Conversion de l'image PG (.nii.gz, format Gate 10) vers le format .mhd
attendu par GateImageOfHistograms de Gate 9.


Gate 9 exige :
 - un fichier .mhd/.raw 4D (dimensions [X, Y, Z, bins] au format MetaImage)
 - deux champs utilisateur "HistoMin" et "HistoMax" dans l'en-tete .mhd,
   correspondant aux bornes reelles du spectre en energie (Emin, Emax)
 - si setTof=true est utilise, un second fichier .mhd/.raw pour le TOF,
   nomme "<nom>-tof.mhd" (meme geometrie 4D, valeurs factices possibles
   car non utilisees pour la construction du spectre en energie)


-- A ADAPTER : remplace Emin_reel / Emax_reel par les vraies bornes
   utilisees par ton pipeline Stage 1 pour cette image --
"""


import itk
import numpy as np


input_path = "/home/jpanagonou/Desktop/vpgTLE-stage/stage2/gate9/data/VPG_vpg_jeanpaul_gamma_e.nii.gz"
output_path = "/home/jpanagonou/Desktop/vpgTLE-stage/stage2/gate9/data/vpg_gate9.mhd"
Emin_reel = 0.0   
Emax_reel = 10.0  


# lecture de l'image source (numpy : [bins, Z, Y, X])
img = itk.imread(input_path)
arr_4d = itk.array_from_image(img)
nbBins, nZ, nY, nX = arr_4d.shape
spacing = img.GetSpacing()  # [sx, sy, sz, s_bins] en ordre ITK (x,y,z,bins)
origin = img.GetOrigin()


print(f"Image source : forme numpy [bins,Z,Y,X] = {arr_4d.shape}")
print(f"Spacing (x,y,z,bins) = {tuple(spacing)}")
print(f"Origin  (x,y,z,bins) = {tuple(origin)}")


# IMPORTANT : ne PAS transposer le tableau avant l'ecriture .raw.
# En convention C (numpy .tofile), c'est le DERNIER axe qui varie le
# plus vite. Le tableau source est deja en [bins, Z, Y, X] : X est
# donc deja le dernier axe et varie deja le plus vite, exactement ce
# que MetaImage attend pour le premier axe declare dans DimSize (X).
arr_to_write = arr_4d.astype(np.float32)


raw_path = output_path.replace(".mhd", ".raw")
arr_to_write.tofile(raw_path)


mhd_header = f"""ObjectType = Image
NDims = 4
BinaryData = True
BinaryDataByteOrderMSB = False
CompressedData = False
TransformMatrix = 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1
Offset = {origin[0]} {origin[1]} {origin[2]} 0
CenterOfRotation = 0 0 0 0
ElementSpacing = {spacing[0]} {spacing[1]} {spacing[2]} 1
DimSize = {nX} {nY} {nZ} {nbBins}
HistoMin = {Emin_reel}
HistoMax = {Emax_reel}
ElementType = MET_FLOAT
ElementDataFile = {raw_path.split('/')[-1]}
"""


with open(output_path, "w") as f:
    f.write(mhd_header)


print(f"\nFichier ecrit : {output_path}")
print(f"Fichier raw   : {raw_path}")
print(f"HistoMin = {Emin_reel}, HistoMax = {Emax_reel}")


# ------------------------------------------------------------------
# creation d'un fichier TOF factice (requis car setTof=true dans la
# macro Gate 9, meme si son contenu n'est pas utilise pour construire
# le spectre en energie, qui repose uniquement sur l'image principale)
# ------------------------------------------------------------------
tof_arr = np.zeros_like(arr_4d, dtype=np.float32)
tof_raw_path = output_path.replace(".mhd", "-tof.raw")
tof_arr.tofile(tof_raw_path)


tof_mhd_path = output_path.replace(".mhd", "-tof.mhd")
tof_header = f"""ObjectType = Image
NDims = 4
BinaryData = True
BinaryDataByteOrderMSB = False
CompressedData = False
TransformMatrix = 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1
Offset = {origin[0]} {origin[1]} {origin[2]} 0
CenterOfRotation = 0 0 0 0
ElementSpacing = {spacing[0]} {spacing[1]} {spacing[2]} 1
DimSize = {nX} {nY} {nZ} {nbBins}
HistoMin = 0.0
HistoMax = 10.0
ElementType = MET_FLOAT
ElementDataFile = {tof_raw_path.split('/')[-1]}
"""
with open(tof_mhd_path, "w") as f:
    f.write(tof_header)


print(f"Fichier TOF factice ecrit : {tof_mhd_path}")
print(f"Fichier TOF raw ecrit     : {tof_raw_path}")


print("\nATTENTION : verifie que Emin_reel/Emax_reel correspondent bien "
      "aux vraies bornes utilisees par ton pipeline Stage 1 pour cette image.")
