#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import print_function, division
"""
Created on Thu Jul  2 10:07:56 2015

@author: Axel Pahl
@title: tools.py

A set for tools to use with the RDKit in the IPython notebook
"""

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

import time
import sys
import base64
import os.path as op

from PIL import Image, ImageChops

import pandas as pd

from . import html_templates as html

from IPython.core.display import HTML

if sys.version_info[0] > 2:
    PY3 = True
    from io import BytesIO as IO
else:
    PY3 = False
    from cStringIO import StringIO as IO

try:
    from misc_tools import apl_tools as apt
    AP_TOOLS = True
except ImportError:
    AP_TOOLS = False

if AP_TOOLS:
    # I use this to keep track of the library versions I use in my project notebooks
    print("{:45s} (commit: {})".format(__name__, apt.get_commit(__file__)))
else:
    print("{:45s} ({})".format(__name__, time.strftime("%y%m%d-%H:%M", time.localtime(op.getmtime(__file__)))))


JSME_OPTIONS = {"css": ["css/style.css", "css/collapsible_list.css"], 
                "scripts": ["lib/jsme/jsme.nocache.js"]}

TBL_JAVASCRIPT = '''<script type="text/javascript">

function toggleCpd(cpdIdent) 
{{
  listPos = document.id_list{ts}.data.value.indexOf(cpdIdent);
  cpdIdentCell = document.getElementById(cpdIdent+"_{ts}");
  if (listPos == -1) 
  {{
    if (document.id_list{ts}.remark.checked == true) 
    {{
      rem = "\\t" + prompt("Remark (Enter for none):", "");
    }}
    else 
    {{
      rem = "";
    }}
    document.id_list{ts}.data.value = document.id_list{ts}.data.value + cpdIdent + rem + "\\n";
    cpdIdentCell.style.backgroundColor = "yellow";
  }}
  else
  {{
    removeStr = cpdIdent;
    tempStr2 = document.id_list{ts}.data.value;
    if (listPos > 0) {{ 
      tempStr1 = tempStr2.substring(0, listPos);
      tempStr2 = tempStr2.substring(listPos, tempStr2.length);
    }} else {{
      tempStr1 = "";
    }}
    listPos = tempStr2.indexOf("\\n");
    if (listPos < tempStr2.length - 1) {{
      tempStr1 = tempStr1 + tempStr2.substring(listPos+1, tempStr2.length)
    }}
    document.id_list{ts}.data.value = tempStr1;
    cpdIdentCell.style.backgroundColor = "transparent";
  }}
}}


function myShowSelection() {{
  document.location.hash = "#SelectionList";
}}
</script>
'''

ID_LIST = """<br><b><a name="SelectionList">Selection:</a></b>
<form name="id_list{ts}">
<input type="checkbox" name="remark" value="prompt" > Prompt for Remarks<br>
<textarea name="data" cols="70" rows="10"></textarea>
</form>
"""

JSME_FORM = '''<script type="text/javascript" src="lib/jsme/jsme.nocache.js"></script>  
<script type="text/javascript">

function jsmeOnLoad() {{
    //arguments: HTML id, width, height (must be string not number!)
    jsmeApplet{ts} = new JSApplet.JSME("appletContainer{ts}", "380px", "340px", {{
                     //optional parameters
                     "options" : "query,hydrogens"
	}});
}}

function onSubmit() {{
    var drawing = jsmeApplet{ts}.smiles();
    // document.getElementById('jsme_smiles{ts}').value = drawing;
    var command = "{var_name} = Chem.MolFromSmiles('" + drawing + "')";
    console.log("Executing Command: " + command);
    
    var kernel = IPython.notebook.kernel;
    kernel.execute(command);
}}
</script>

<table align="left" style="border: none;">
<tr style="border: none;">
<td id="appletContainer{ts}" style="border: none;"></td>
<td style="vertical-align: bottom; border: none;">
<button onclick="onSubmit()">done !</button>
</td>
</tr>
</table>
'''


class Mol_List(list):
    """For now, a simple wrapper around built-in list.
    Enables display of molecule lists as HTML tables in IPython notebook just by-call
    (via _repr_html)."""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def __getitem__(self, item):
        result = list.__getitem__(self, item)
        try:
            return Mol_List(result)
        except TypeError:
            return result
        
    def _repr_html_(self):
        return mol_table(self)


def autocrop(im, bgcolor="white"):
    if im.mode != "RGB":
        im = im.convert("RGB")
    bg = Image.new("RGB", im.size, bgcolor)
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox:
         return im.crop(bbox)
    return None # no contents


def list_fields(sdf_list):
    field_list = []

    for mol in sdf_list:
        field_list.extend(mol.GetPropNames())

    return list(set(field_list))


def guess_id_prop(prop_list):  # try to guess an id_prop
    for prop in prop_list:
        if prop.lower().endswith("id"):
            return prop
    return None


def get_value(str_val):
    try:
        val = float(str_val)
        if val == int(val):
            val = int(val)
    except ValueError:
        val = str_val
    return val


def mol_table(sdf_list, id_prop=None, highlight=None):
    """
    input:   list of RDKit molecules
    highlight: dict of properties (special: *all*) and values to highlight cells,
               e.g. {"activity": "< 50"}
    returns: HTML table as TEXT to embed in IPython or a web page."""
    
    time_stamp = time.strftime("%y%m%d%H%M%S")
    td_opt = {"align": "center"}
    header_opt = {"bgcolor": "#94CAEF"}
    table_list = []
    prop_list = list_fields(sdf_list)
    if id_prop:
        table_list.append(TBL_JAVASCRIPT.format(ts=time_stamp))
        if not id_prop in prop_list:
            raise LookupError("id_prop not found in data set.")
        guessed_id = id_prop
    else: # try to guess an id_prop
        guessed_id = guess_id_prop(prop_list)
        
    if guessed_id:
        # make sure that the id_prop (or the guessed id prop) is first:
        prop_list.pop(prop_list.index(guessed_id))
        tmp_list = [guessed_id]
        tmp_list.extend(prop_list)
        prop_list = tmp_list

    cells = html.td(html.b("#"), header_opt)
    cells.extend(html.td(html.b("Molecule"), header_opt))
    for prop in prop_list:
        cells.extend(html.td(html.b(prop), header_opt))
    rows = html.tr(cells)

    
    for idx, mol in enumerate(sdf_list):
        cells = []
        mol_props = mol.GetPropNames()
        
        if id_prop:
            id_prop_val = mol.GetProp(id_prop)
            cell_opt = {"id": "{}_{}".format(id_prop_val, time_stamp)}
        else:
            cell_opt = {"id": str(idx)}
        cell = html.td(str(idx), cell_opt)
        cells.extend(cell)
        
        if not mol:
            cells.extend(html.td("no structure"))
        
        else:
            img_file = IO()
            img = autocrop(Draw.MolToImage(mol))
            img.save(img_file, format='PNG')

            b64 = base64.b64encode(img_file.getvalue())
            if PY3:
                b64 = b64.decode()
            img_file.close()
            
            if id_prop:
                img_opt = {"title": "Click to select / unselect", 
                           "onclick": "toggleCpd('{}')".format(id_prop_val)}
            else:
                img_opt = {"title": str(idx)}

            img_src = "data:image/png;base64,{}".format(b64)
            cells.extend(html.td(html.img(img_src, img_opt)))
        
        for prop in prop_list:
            td_opt = {"align": "center"}
            if prop in mol_props:
                td_opt["title"] = prop
                prop_val = mol.GetProp(prop)
                if highlight:
                    eval_str = None
                    if "*all*" in highlight:
                        if not guessed_id or (guessed_id and prop != guessed_id):
                            eval_str = " ".join([prop_val, highlight["*all*"]])
                    else:
                        if prop in highlight:
                            eval_str = " ".join([prop_val, highlight[prop]])
                    if eval_str and eval(eval_str):
                        td_opt["bgcolor"] = "#99ff99"
                    
                cells.extend(html.td(prop_val, td_opt))
            else:
                cells.extend(html.td("", td_opt))
        
        rows.extend(html.tr(cells))
    
    table_list.extend(html.table(rows))
    
    if id_prop:
        table_list.append(ID_LIST.format(ts=time_stamp))

    # print(table_list)
    return "".join(table_list)


def show_table(sdf_list, id_prop=None, highlight=None):
    return HTML(mol_table(sdf_list, id_prop, highlight=highlight))


def jsme(name="mol"):
    """displays a JSME molecule editor widget in the notebook 
    and stores the resulting mol in the variable that <name> assigns."""
    
    time_stamp = time.strftime("%y%m%d%H%M%S")
    
    return HTML(JSME_FORM.format(ts=time_stamp, var_name=name))


def dict_from_sdf_list(sdf_list, id_prop=None, props=None, prop_list=None):
    """Generate a dictionary from the properties of a list of molecules.
    Currently not including the structure.
    If <props> contains a list of property names, then only these properties plus the <id_prop> are returned.
    Returns dict"""
    
    if not prop_list:
        prop_list = list_fields(sdf_list)
        
    if id_prop:
        if not id_prop in prop_list:
            raise LookupError("id_prop not found in data set.")
        guessed_id = id_prop
    else: 
        guessed_id = guess_id_prop(prop_list)

    if not props:
        props = prop_list
    if guessed_id and not guessed_id in props:
        props.append(guessed_id)
    
    df_dict = {prop: [] for prop in props}

    for mol in sdf_list:
        mol_props = list(mol.GetPropNames())
        for prop in props:
            if prop in mol_props:
                df_dict[prop].append(get_value(mol.GetProp(prop)))
            else:
                df_dict[prop].append(pd.np.NaN)

    return df_dict


# some convenience functions
def mol_3d(smiles_or_mol):
    """return a 3d optimized molecule from a Smiles or 2d mol input"""
    if isinstance(smiles_or_mol, str):  # input is Smiles
        smiles_or_mol = Chem.MolFromSmiles(smiles_or_mol)
    
    mh = Chem.AddHs(smiles_or_mol)
    Chem.Compute2DCoords(mh)
    Chem.EmbedMolecule(mh)
    Chem.MMFFOptimizeMolecule(mh)
    return mh


def mol_grid(sdf_list, props, fn=None, mols_per_row=5, sub_img_size=(200, 200)):
    """Draw a molecule grid from the input <sdf_list>. An inline graphics will be returned
    in addition to writing the image to <fn> (if defined).
    The given sdf <props> (as a list) will be concatenated to the molecules' legends."""

    if not isinstance(props, list):
        props = [props]

    legends = []
    for mol in sdf_list:
        leg = [mol.GetProp(prop) for prop in props]
        leg_str = "_".join(leg)
        legends.append(leg_str)

    img = Draw.MolsToGridImage(sdf_list, molsPerRow=mols_per_row, subImgSize=sub_img_size, legends=legends)
    if fn:
        img.save(fn)
    return img


def o3da(input_list, ref, fn="aligned.sdf"):
    """Takes a list of molecules and align them to ref.
    Writes the result as SD file to fn."""
    ref_pymp = Chem.MMFFGetMoleculeProperties(ref)
    mol_list = input_list[:]
    writer = Chem.SDWriter(fn)

    print("N\t\tscore\t\trmsd")
    for ctr, mol in enumerate(mol_list, 1):
        mol_pymp = Chem.MMFFGetMoleculeProperties(mol)
        o3a = Chem.GetO3A(mol, ref, mol_pymp, ref_pymp)
        print("{}\t\t{:.2f}\t\t{:.2f}".format(ctr, o3a.Score(),o3a.Align()))
        writer.write(mol)

    writer.close()
