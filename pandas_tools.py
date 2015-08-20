# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 12:20:19 2015

@author: pahl
@title: pandas_tools.py
"""

import pandas as pd

from rdkit.Chem import AllChem as Chem

from . import tools, hc_tools as hct


def df_from_sdf_list(sdf_list, id_prop=None, props=None, set_index=True):
    """Generate a Pandas dataframe from the properties of a list of molecules.
    Currently not including the structure.
    If <props> contains a list of property names, then only these properties plus the <id_prop> are returned.
    Returns Pandas dataframe"""

    prop_list = tools.list_fields(sdf_list)
        
    if id_prop:
        guessed_id = id_prop
    else: 
        guessed_id = tools.guess_id_prop(prop_list)

    df_dict = hct.dict_from_sdf_list(sdf_list, id_prop=id_prop, props=props, prop_list=prop_list)
    df = pd.DataFrame(df_dict)
    if set_index and guessed_id:
        df = df.set_index(guessed_id)
    
    return df


def join_data_from_file(df, fn):
    """
    join data from file (e.g. Smiles, biol. data (tab-sep.) ) by index (set to Compound_Id) to df
    index in df has to be set to Compound_Id
    """
    data_df = pd.concat([df, pd.read_table(fn, index_col=df.index.name)], axis=1, join_axes=[df.index])
    
    return data_df


def df_to_sdf_list(df):
    """
    returns: list of mols
    """
    pass


def df_show_table(df):
    """
    show df as mol_table
    """
    pass


def move_col_left(df, col):
    """
    # this moves the last col to the left:
    cols = list(mol_df.columns.values)
    new_cols = [cols[-1]]
    new_cols.extend(cols[:-1])
    mol_df = mol_df[new_cols]
    
    return df
    """
    pass


def align_molecules(df, qry, mol_col="molecule"):
    """
    Align all molecules in a df to a given query molecule,
    e.g. after a substructure query.
    operates inline, returns nothing"""
    
    Chem.Compute2DCoords(qry)
    for mol in df[mol_col]:
        Chem.GenerateDepictionMatching2DStructure(mol, qry)
    