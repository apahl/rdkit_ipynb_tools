#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#########
Scaffolds
#########

*Created on Thu Jul 11 14:00 2016 by A. Pahl*

Scaffold clustering functionality.
"""

import binascii
from collections import namedtuple

from rdkit.Chem import AllChem as Chem
import rdkit.Chem.Scaffolds.MurckoScaffold as MurckoScaffold

from . import tools

# ri = mol.GetRingInfo()
# ri.NumRings()

# Scaffolds:
# MurckoFrame; MurckScaf; Scaffolds...


class Mol_List(tools.Mol_List):
    """Adding scaffold clustering functionality to the Mol_list."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def generate_scaffolds(self):
        "Calculate the scaffolds, if they are not already present."
        for mol in self:
            if mol.HasProp("Scaf_Smiles"): continue
            scaffolds = calc_all_scaffolds(mol)
            mol.SetProp("Scaf_Smiles", scaffolds.smiles)
            mol.SetProp("Scaf_Ids", scaffolds.ids)


def calc_scaffold_id(smiles):
    "Calculate the scaffold Id as string from its Smiles (the CRC32 checksum is used for this for now."
    return str(binascii.crc32(smiles.encode()))


def calc_murcko_scaf(mol):
    "Calculate the Murcko scaffold from a molecule as Smiles."
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)


def calc_murcko_frame(mol):
    """Calculate the Murcko generic frame from a molecule as Smiles."""
    return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(mol))


def calc_child_scaffolds(scaf):
    """Recursively generate the possible child scaffolds. from the input MurckoScaffold Smiles.
    Returns the scaffolds' Smiles and their Ids as lists in a namedtuple(Scaffolds)"""
    scaf_list = []
    scaf_id_list = []

    def _recurse(scaf):
        orig_mol = Chem.MolFromSmiles(scaf)
        rwmol = Chem.RWMol(orig_mol)
        ri = rwmol.GetRingInfo()
        if ri.NumRings() < 3:
            return
        bonds = rwmol.GetBonds()
        for bond in bonds:
            if not bond.IsInRing():
                rwmol = Chem.RWMol(orig_mol)
                rwmol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                frags = rwmol.GetMol()
                frag_list = Chem.MolToSmiles(frags).split(".")
                ring_split = 0
                rings_per_frag = []
                for frag in frag_list: # have we split between two rings?
                    if len(frag) > 2:
                        mol = Chem.MolFromSmiles(frag)
                        ri = mol.GetRingInfo()
                        num_rings = ri.NumRings()
                        rings_per_frag.append(num_rings)
                        if num_rings > 0:
                            ring_split += 1

                if ring_split >= 2:
                    for idx, frag in enumerate(frag_list):
                        if rings_per_frag[idx] > 1:
                            murcko_frag = MurckoScaffold.MurckoScaffoldSmiles(frag)
                            if murcko_frag not in scaf_list:
                                scaf_list.append(murcko_frag)
                                _recurse(murcko_frag)

    _recurse(scaf)
    scaf_list.sort(key=len, reverse=True)
    for scaf in scaf_list:
        scaf_id_list.append(calc_scaffold_id(scaf))

    ScaffoldsList = namedtuple("ScaffoldsList", ["level_smiles", "ids"])
    return ScaffoldsList(smiles=scaf_list, ids=scaf_id_list)


def calc_all_scaffolds(mol):
    """Calculate all possible scaffolds, according to this scheme:
    MurckoFrame; MurckScaf; Scaffolds..."""
    scaf_list = []
    scaf_id_list = []
    # MurckoFrame, MurckoScaf
    for level, calculator in enumerate([calc_murcko_frame, calc_murcko_scaf]):
        scaf = calculator(mol)
        scaf_id = calc_scaffold_id(scaf)
        scaf_list.append(str(level) + "_" + scaf)
        scaf_id_list.append(scaf_id)

    # remaining scaffolds...
    child_scaffolds = calc_child_scaffolds(scaf_list[1])  # start from the MurckoScaffold
    scaf_list.extend(child_scaffolds.level_smiles)
    scaf_id_list.extend(child_scaffolds.ids)

    Scaffolds = namedtuple("Scaffolds", ["level_smiles", "ids"])
    return Scaffolds(smiles="; ".join(scaf_list), ids="; ".join(scaf_id_list))
