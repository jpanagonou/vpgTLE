#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#=======================================================================
# variables.py — version validation Oxygen (Gate 9 vs Gate 10)
# Adapté pour une boite d'Oxygen non voxélisée
#=======================================================================

import uproot
import opengate as gate
from opengate.definitions import elements_name_symbol
import os
import numpy as np
from stage_1a_oxygen import *
from box import Box
import pathlib

gcm3 = gate.g4_units.g_cm3


class tle_sim:
    def __init__(self, sim):
        self.path = paths
        self.sim  = sim

        # Volume et acteur
        self.ct  = sim.volume_manager.volumes[vol_name]
        self.act = sim.actor_manager.actors[actor_name]

        # Liste des sorties activées
        actors_list = []
        if self.act.prot_E.active:
            actors_list.append("prot_e.nii.gz")
        if self.act.neutr_E.active:
            actors_list.append("neutr_e.nii.gz")
        self.actor_list = actors_list

        # ── Base de données Oxygen ─────────────────────────────────────
        # Pour la validation, on utilise uniquement la base Oxygen
        f3 = str(paths.data / "db-Oxygen_merged.root")
        self.root_file         = uproot.open(f3)
        self.root_file_neutron = uproot.open(f3)  # même fichier, pas de neutrons

        # ── Matériaux ──────────────────────────────────────────────────
        # Boite d'Oxygen pure → un seul matériau
        # On crée manuellement la liste de matériaux
        # format : [HU_min, HU_max, nom]
        self.voxel_mat = [[-3024, 10000, "Oxygen"]]
        self.mat_fraction = [{"name": "Oxygen", "HU": 0, "Oxygen": 100.0}]
        self.el = ["O"]

        # Base de données de matériaux (pas nécessaire pour Oxygen pur)
        self.data_way = str(paths.data / "database.db")



    def find_emission_vector(self, el, root_data):
        """Trouve le TH2D dans root_data correspondant à l'élément el."""
        # db-Oxygen_merged.root utilise le nom complet "Oxygen" pas le symbole "O"
        histo = root_data["Oxygen"]["GammaZ"].to_hist()
        w = histo.to_numpy()[0]
        return w

   

    def density_mat(self, mat, data_way):
        """Retourne la densité de l'Oxygen (1 g/cm3 à 20°C)."""
        # Densité de l'Oxygen gazeux
        return 1 # c'est 1 car dans le GateMaterials data base utilisé, Rho vaut & pour l'oxygene


    def voxel_to_mat_name(self, UH, mat_data):
        """Retourne toujours G4_O pour la boite d'Oxygen."""
        return "Oxygen"


    def liste_el_frac(self, mat, frac_data):
        """Retourne les éléments et fractions de l'Oxygen pur."""
        return ["Oxygen"], [100.0]



    def gamma_mat(self, name, root_data):
        """Construit le GammaZ pour l'Oxygen."""
        Gamma = np.zeros((250, 250))
        rho_mat = self.density_mat(name, self.data_way)
        elements = ["Oxygen"]
        fractions = [100.0]

        for el in elements:
            w    = fractions[elements.index(el)] / 100
            vect = self.find_emission_vector(el, root_data)
            Gamma += vect * w
        return Gamma

    def redim(self, ind, ct, actor_vol, array):
        """Pas nécessaire pour une boite simple — retourne les indices tels quels."""
        return ind[0], ind[1], ind[2]
