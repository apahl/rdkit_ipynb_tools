# -*- coding: utf-8 -*-
"""
############
Pandas Tools
############

*Created on Wed Jul 29 12:20:19 2015 by A. Pahl*

A Set of Pandas Tools for RDKit.
(built on top of rdkit.Chem.PandasTools)
"""

import time
import os.path as op
import pandas as pd

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools as PT

from . import tools, hc_tools as hct


try:
    from misc_tools import apl_tools as apt
    AP_TOOLS = True
except ImportError:
    AP_TOOLS = False

if AP_TOOLS:
    #: Library version
    VERSION = apt.get_commit(__file__)
    # I use this to keep track of the library versions I use in my project notebooks
    print("{:45s} (commit: {})".format(__name__, VERSION))
else:
    print("{:45s} ({})".format(__name__, time.strftime("%y%m%d-%H:%M", time.localtime(op.getmtime(__file__)))))


def init_PT():
    """create a dummy df to initialize the RDKit PandasTools (a bit hacky, I know)."""
    init = pd.DataFrame.from_dict({"id": [123, 124], "Smiles": ["c1ccccc1C(=O)N", "c1ccccc1C(=O)O"]})
    PT.AddMoleculeColumnToFrame(init)

init_PT()

def move_col(df, col, new_pos=1):
    """
    Put column col on position new_pos.
    """

    cols = list(df.columns.values)
    index = cols.index(col)  # raises ValueError if not found
    cols_len = len(cols)
    if new_pos > cols_len - 1:
        new_pos = cols_len - 1

    if index == new_pos:
        print("{} is already on position {}".format(col, new_pos))
        return df

    new_cols = []
    for i, c in enumerate(cols):
        if i == index: continue
        if i == new_pos:
            new_cols.append(col)
        new_cols.append(c)

    new_df = df[new_cols]

    return new_df


def df_from_mol_list(mol_list, id_prop="Compound_Id", props=None, set_index=True):
    """Generate an RDKit Pandas dataframe from a Mol_List.
    Now also including the structure.
    If <props> contains a list of property names, then only these properties plus the <id_prop> are returned.
    Returns RDKit Pandas dataframe"""

    prop_list = tools.list_fields(mol_list)

    if id_prop:
        guessed_id = id_prop
    else:
        guessed_id = tools.guess_id_prop(prop_list)

    df_dict = tools.dict_from_sdf_list(mol_list, id_prop=id_prop, props=props, prop_list=prop_list)
    smiles_col = [Chem.MolToSmiles(mol) for mol in mol_list]
    df = pd.DataFrame(df_dict)
    if set_index and guessed_id:
        df = df.set_index(guessed_id)

    df["Smiles"] = pd.Series(data=smiles_col, index=df.index)
    PT.AddMoleculeColumnToFrame(df, smilesCol="Smiles", molCol='mol', includeFingerprints=False)
    df = df.drop("Smiles", axis=1)
    # move structure column to left side
    df = move_col(df, "mol", new_pos=0)

    return df


def left_join_on_index(df1, df2):
    new_df = pd.merge(df1, df2, how="left", left_index=True, right_index=True)
    return new_df


def inner_join_on_index(df1, df2):
    new_df = pd.merge(df1, df2, how="inner", left_index=True, right_index=True)
    return new_df


def join_data_from_file(df, fn, dropna=True, gen_struct=True, remove_smiles=True):
    """Join data from file (e.g. Smiles, biol. data (tab-sep.) ) by index (set to Compound_Id) to df.
    If gen_struct is True and Smiles could be found in the resulting df,
    then the structure is generated and moved to first column.
    If remove_smiles is True, the Smiles column will be dropped.
    index in df has to be set to the Compound_Id.
    returns new df"""

    # data_df = pd.concat([df, pd.read_table(fn, index_col=df.index.name)], axis=1, join_axes=[df.index])
    data_df = inner_join_on_index(df, pd.read_table(fn, index_col=df.index.name))

    not_found = list(set(df.index.tolist()) - set(data_df.index.tolist()))

    if dropna:
        data_df = data_df.dropna(axis=1, how="all")

    if gen_struct:
        smiles_col = None
        for col in list(data_df.columns.values):
            if "smiles" in col.lower():
                smiles_col = col
                break

        if smiles_col:
            PT.AddMoleculeColumnToFrame(data_df, smilesCol=smiles_col, molCol='mol', includeFingerprints=False)

            if remove_smiles:
                data_df = data_df.drop(smiles_col, axis=1)

            data_df = move_col(data_df, "mol", new_pos=0)

    if not_found:
        print("* not found:", not_found)

    data_df = data_df.convert_objects(convert_numeric=True) # convert to numeric where possible

    return data_df


def keep_numeric_only(df):
    """Keep only the numeric data in a df,
    remove all ROWS that contain non-numeric data.
    Do this prior to a highchart plot
    returns new df"""
    new_df = df.convert_objects(convert_numeric=True)
    new_df = new_df.dropna()
    return new_df


def mol_list_from_df(df, mol_col="mol"):
    """
    Creates a Mol_List from an RDKit Pandas dataframe.
    Returns Mol_List.
    """

    mol_list = tools.Mol_List()
    id_prop = df.index.name
    props = [k for k in df.keys().tolist() if k != mol_col]

    for cid in df.index.values.tolist():
        mol = df.at[cid, mol_col]
        if not mol:
            continue
        mol.SetProp(id_prop, str(cid))
        for prop in props:
            if df.at[cid, prop]:
                mol.SetProp(prop, str(df.at[cid, prop]))
        mol_list.append(mol)

    return mol_list


def df_show_table(df):
    """
    show df as mol_table
    """
    pass


def align_molecules(df, qry, mol_col="mol"):
    """
    Align all molecules in a df to a given query molecule,
    e.g. after a substructure query.
    operates inline, returns nothing"""

    Chem.Compute2DCoords(qry)
    for mol in df[mol_col]:
        Chem.GenerateDepictionMatching2DStructure(mol, qry)


def add_calc_prop(df, mol_col="mol", props="logp"):
    avail_props = ["logp", "mw", "sa"]
    if not isinstance(props, list):
        prop = list(props)
    for prop in props:
        if not prop in avail_props:
            raise ValueError("{} can not be calculated.".format(prop))
    #TODO: fill stub
