# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 12:20:19 2015

@author: pahl
@title: pandas_tools.py
"""

import pandas as pd

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
