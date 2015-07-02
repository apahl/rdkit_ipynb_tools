# -*- coding: utf-8 -*-
from __future__ import print_function, division
"""
Created on Thu Jul  2 10:07:56 2015

@author: Axel Pahl

A set for tools to use with the RDKit in the IPython notebook
"""
# tools.py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# html_mol_table.py
"""
Created on Wed Jun 02 2015

@author: pahl
"""

from rdkit.Chem import Draw
Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

import time
import sys
import base64
import os.path as op

from PIL import Image, ImageChops

from . import html_templates as html

from IPython.core.display import HTML

if sys.version_info[0] > 2:
    PY3 = True
    from io import BytesIO as IO
else:
    PY3 = False
    from cStringIO import StringIO as IO

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


def mol_table(sdf_list, id_prop=None):
    """
    input:   list of RDKit molecules
    returns: HTML table as TEXT to embed in IPython or a web page."""
    
    time_stamp = time.strftime("%y%m%d%H%M%S")
    td_opt = {"align": "center"}
    header_opt = {"bgcolor": "#94CAEF"}
    table_list = []
    if id_prop:
        table_list.append(TBL_JAVASCRIPT.format(ts=time_stamp))

    prop_list = list_fields(sdf_list)
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
            if prop in mol_props:
                cells.extend(html.td(mol.GetProp(prop), td_opt))
            else:
                cells.extend(html.td("", td_opt))
        
        rows.extend(html.tr(cells))
    
    table_list.extend(html.table(rows))
    
    if id_prop:
        table_list.append(ID_LIST.format(ts=time_stamp))

    # print(table_list)
    return "".join(table_list)


def show_table(sdf_list, id_prop=None):
    return HTML(mol_table(sdf_list, id_prop))


def jsme(name="mol"):
    """displays a JSME molecule editor widget in the notebook 
    and stores the resulting mol in the variable that <name> assigns."""
    
    time_stamp = time.strftime("%y%m%d%H%M%S")
    
    return HTML(JSME_FORM.format(ts=time_stamp, var_name=name))
