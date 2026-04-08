#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#=======================================================================
# Rééchantillonnage du CT à différentes résolutions
# Entrée  : ct_4mm.mhd  (spacing 4mm)
# Sorties : ct_10mm.mhd (spacing 10mm)
#           ct_20mm.mhd (spacing 20mm)
#=======================================================================

import itk
import numpy as np

# ── Chargement du CT original ──────────────────────────────────────────
img = itk.imread('ct_4mm.mhd', itk.SS)

orig_size    = img.GetLargestPossibleRegion().GetSize()
orig_spacing = img.GetSpacing()
orig_origin  = img.GetOrigin()
orig_dir     = img.GetDirection()

print(f"CT original :")
print(f"  Size    : {orig_size[0]} x {orig_size[1]} x {orig_size[2]}")
print(f"  Spacing : {orig_spacing[0]} x {orig_spacing[1]} x {orig_spacing[2]} mm")

# ── Fonction de rééchantillonnage ──────────────────────────────────────
def resample_ct(img, new_spacing_mm, output_filename):
    """
    Rééchantillonne une image CT à un nouveau spacing.
    
    Paramètres
    ----------
    img             : image ITK originale
    new_spacing_mm  : nouveau spacing en mm (float)
    output_filename : chemin du fichier de sortie
    """
    new_spacing = [float(new_spacing_mm)] * 3

    # Calculer la nouvelle taille pour conserver le volume physique
    orig_size    = img.GetLargestPossibleRegion().GetSize()
    orig_spacing = img.GetSpacing()
    new_size = [
        int(orig_size[0] * orig_spacing[0] / new_spacing[0]),
        int(orig_size[1] * orig_spacing[1] / new_spacing[1]),
        int(orig_size[2] * orig_spacing[2] / new_spacing[2])
    ]

    nb_voxels = new_size[0] * new_size[1] * new_size[2]
    mem_mb    = nb_voxels * 250 * 4 * 4 / (1024**2)  # 4 sorties, float32

    print(f"\nRééchantillonnage à {new_spacing_mm} mm :")
    print(f"  Nouvelle taille : {new_size[0]} x {new_size[1]} x {new_size[2]}")
    print(f"  Nb voxels       : {nb_voxels:,}")
    print(f"  Mémoire estimée : {mem_mb:.1f} MB")

    # Rééchantillonnage avec interpolation linéaire
    ImageType     = type(img)
    InterpolatorType = itk.LinearInterpolateImageFunction[ImageType, itk.D]
    interpolator  = InterpolatorType.New()

    resampler = itk.ResampleImageFilter[ImageType, ImageType].New()
    resampler.SetInput(img)
    resampler.SetInterpolator(interpolator)
    resampler.SetOutputSpacing(new_spacing)
    resampler.SetSize(new_size)
    resampler.SetOutputOrigin(img.GetOrigin())
    resampler.SetOutputDirection(img.GetDirection())
    resampler.SetDefaultPixelValue(-1000)
    resampler.Update()

    # Sauvegarde
    itk.imwrite(resampler.GetOutput(), output_filename)
    print(f"  Sauvegardé : {output_filename}")


# ── Création des deux nouveaux CT ─────────────────────────────────────
resample_ct(img, 10.0, '../data/ct_10mm.mhd')
resample_ct(img, 20.0, '../data/ct_20mm.mhd')

print("\nTerminé !")
print("\nRécapitulatif :")
print(f"  ct_40mm.mhd :  13 x  13 x  19 voxels  (déjà existant)")
print(f"  ct_20mm.mhd :  à vérifier")
print(f"  ct_10mm.mhd :  à vérifier")
print(f"  ct_4mm.mhd  : 125 x 125 x 189 voxels  (original)")
