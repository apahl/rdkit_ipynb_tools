# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 12:20:19 2015

@author: pahl
@title: pandas_tools.py
"""

import pandas as pd

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import PandasTools as PT

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
    

def join_data_from_file(df, fn, dropna=True, gen_struct=True, remove_smiles=True):
    """
    Join data from file (e.g. Smiles, biol. data (tab-sep.) ) by index (set to Compound_Id) to df.
    If gen_struct is True and Smiles could be found in the resulting df,
      then the structure is generated and moved to first column.
    If remove_smiles is True, the Smiles column will be dropped.
    index in df has to be set to the Compound_Id.
    returns new df
    """
    data_df = pd.concat([df, pd.read_table(fn, index_col=df.index.name)], axis=1, join_axes=[df.index])
    
    if dropna:
        data_df = data_df.dropna(axis=1)
    
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
    
    return data_df


def keep_numeric_only(df):
    """Keep only the numeric data in a df, 
    remove all ROWS that contain non-numeric data.
    returns new df"""
    new_df = df.convert_objects(convert_numeric=True)
    new_df = new_df.dropna()
    return new_df


def left_join_on_index(df1, df2):
    new_df = pd.merge(df1, df2, how="left", left_index=True, right_index=True)
    return new_df


def df_to_sdf_list(df):
    """
    returns: list of mols
    """
    #TODO: fill stub
    pass


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
