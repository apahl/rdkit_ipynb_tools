#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#####
Tools
#####

*Created on Thu Jul  2 10:07:56 2015 by A. Pahl*

A set of tools to use with the `RDKit <http://rdkit.org>`_ in the IPython notebook.
"""

import time
import sys
import base64
import os
import os.path as op
import random
import csv
import gzip
import math
from copy import deepcopy

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw, rdFMCS
import rdkit.Chem.Descriptors as Desc

# imports for similarity search
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs

import rdkit.Chem.Scaffolds.MurckoScaffold as MurckoScaffold

Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

from PIL import Image, ImageChops

import numpy as np

from . import html_templates as html
from . import hc_tools as hct

try:
    import ipywidgets as ipyw
    WIDGETS = True
except ImportError:
    WIDGETS = False

from IPython.core.display import HTML, display

if sys.version_info[0] > 2:
    PY3 = True
    from io import BytesIO as IO
else:
    PY3 = False
    from cStringIO import StringIO as IO

try:
    from . import bokeh_tools as bkt
    PLOT_TOOL = "bokeh"

except ImportError:
    print("  * could not import Bokeh, plotting with Highcharts instead.")
    PLOT_TOOL = "highcharts"

try:
    from misc_tools import apl_tools as apt
    AP_TOOLS = True
except ImportError:
    AP_TOOLS = False

try:
    # Try to import Avalon so it can be used for generation of 2d coordinates.
    from rdkit.Avalon import pyAvalonTools as pyAv
    USE_AVALON = True
except ImportError:
    print("  * Avalon not available. Using RDKit for 2d coordinate generation.")
    USE_AVALON = False

try:
    from Contrib.SA_Score import sascorer
    SASCORER = True
except ImportError:
    print("* SA scorer not available. RDKit's Contrib dir needs to be in the Python import path...")
    SASCORER = False


if AP_TOOLS:
    #: Library version
    VERSION = apt.get_commit(__file__)
    # I use this to keep track of the library versions I use in my project notebooks
    print("{:45s} (commit: {})".format(__name__, VERSION))
else:
    print("{:45s} ({})".format(__name__, time.strftime("%y%m%d-%H:%M", time.localtime(op.getmtime(__file__)))))

if op.isfile("lib/jsme/jsme.nocache.js"):
    JSME_LOCATION = "lib"
else:
    print("- no local installation of JSME found, using web version.")
    JSME_LOCATION = "http://peter-ertl.com/jsme/JSME_2016-07-31"


BGCOLOR = "#94CAEF"
IMG_GRID_SIZE = 235

# A list of partial property strings to use for ordering of properties:
DEFAULT_ORDER = ["_id", "supplier", "producer", "activity|pic50",
                 "hit", "actass", "pure_flag", "purity", "identity", "lcms"]

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
    cpdIdentCell.style.backgroundColor = "{bgcolor}";
  }}
  // display the number of selected compounds:
  var count = (document.id_list{ts}.data.value.match(/\\n/g) || []).length;
  document.getElementById("selection_title{ts}").innerHTML = "Selection (" + count + "):";
}}


function myShowSelection() {{
  document.location.hash = "#SelectionList";
}}
</script>
'''

ID_LIST = """<br><b><a name="SelectionList" id="selection_title{ts}">Selection (0):</a></b>
<form name="id_list{ts}">
<input type="checkbox" name="remark" value="prompt" > Prompt for Remarks<br>
<textarea name="data" cols="70" rows="10"></textarea>
</form>
"""

JSME_FORM = '''<script type="text/javascript" src="{jsme_loc}/jsme/jsme.nocache.js"></script>
<script type="text/javascript">

function jsmeOnLoad() {{
    //arguments: HTML id, width, height (must be string not number!)
    jsmeApplet{ts} = new JSApplet.JSME("appletContainer{ts}", "380px", "340px", {{
                     //optional parameters
                     "options" : "query,hydrogens"
    }});
}}

function onSubmit() {{
    var drawing = jsmeApplet{ts}.molFile();
    // document.getElementById('jsme_smiles{ts}').value = drawing;
    var command = '{var_name} = Chem.MolFromMolBlock("""' + drawing + '""")';
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


class NoFieldTypes(Exception):
    def __str__(self):
        return repr("FieldTypeError: field types could not be extracted from Mol_List")


class Mol_List(list):
    """Enables display of molecule lists as HTML tables in IPython notebook just by-call
    (via _repr_html)."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.order = None
        self.ia = False  # wether the table and grid views are interactive or no
        self.plot_tool = PLOT_TOOL
        self.id_prop = None
        self.recalc_needed = {}
        self._set_recalc_needed()


    def __getitem__(self, item):
        result = list.__getitem__(self, item)
        try:
            new_list = Mol_List(result)

            # pass on properties
            new_list.order = self.order
            new_list.ia = self.ia
            new_list.plot_tool = self.plot_tool
            return new_list
        except TypeError:
            return result


    def _repr_html_(self):
        id_prop = guess_id_prop(list_fields(self)) if self.ia else None
        return mol_table(self, id_prop=id_prop, order=self.order)


    def _set_recalc_needed(self):
        """Make sure that the expensive calculations are not done too often."""

        self.len = len(self)
        self.recalc_needed["plot_tool"] = PLOT_TOOL
        keys = ["d", "fields", "field_types", "id_prop"]
        for k in keys:
            self.recalc_needed[k] = True


    def _key_get_prop(self, mol, field, reverse=False):
        if reverse:
            not_found = -1000000.0
        else:
            not_found = 1000000.0
        try:
            val = float(mol.GetProp(field))
        except ValueError:  # GetProp value could not be converted to float
            val = 0
        except KeyError:   # field is not present in the mol properties
            val = not_found
        return val


    def _get_field_types(self):
        """Detect all the property field types.

        Returns:
            Dict with the property names as keys and the types as values."""

        print("  > detecting field types...")
        field_types = {}

        if len(self) > 100:
            sdf_sample = random.sample(self, len(self) // 2)
        else:
            sdf_sample = self

        for mol in sdf_sample:
            prop_names = mol.GetPropNames()

            for prop in prop_names:
                prop_type = "number"
                prop_str = mol.GetProp(prop)

                try:
                    float(prop_str)
                    if prop.lower().endswith("id"):
                        prop_type = "key"

                except ValueError:
                    prop_type = "str"

                if prop in field_types:
                    if field_types[prop] in ["number", "key"] and prop_type == "str":
                        # "str" overrides everything: if one string is among the values
                        # of a property, all values become type "str"
                        field_types[prop] = prop_type
                else:
                    field_types[prop] = prop_type

        if not field_types:
            raise NoFieldTypes()

        return field_types


    def _calc_d(self):
        self._d = {x: [] for x in self.fields}
        self._d["mol"] = []
        for mol in self:
            if not mol: continue
            if self.plot_tool == "bokeh":
                img_tag = b64_img(mol)
            else:
                img_tag = '<img src="data:image/png;base64,{}" alt="Mol"/>'.format(b64_img(mol))

            self._d["mol"].append(img_tag)

            for prop in self.fields:
                if mol.HasProp(prop):
                    self._d[prop].append(get_value(mol.GetProp(prop)))
                else:
                    self._d[prop].append(np.nan)


    def append(self, other):
        self._set_recalc_needed()
        super().append(other)


    def extend(self, other):
        self._set_recalc_needed()
        super().extend(other)


    def align(self, mol_or_smiles=None):
        """Align the Mol_list to the common substructure provided as Mol or Smiles.

        Args:
            mol_or_smiles (bool): The substructure to which to align.
                If None, then the method uses rdFMCS to determine the MCSS
                of the Mol_List."""

        align(self, mol_or_smiles)

        self.recalc_needed["d"] = True


    def add_id(self, id_prop="molid"):
        """Add an Id property ``id_prop`` to the Mol_List.
        By default, "molid" is used."""

        for idx, mol in enumerate(self, 1):  # start at index 1
            mol.SetProp(id_prop, str(idx))


    def write_sdf(self, fn, conf_id=-1):
        """Write Mol_List instance as SD File"""

        writer = Chem.SDWriter(fn)

        # try to save the column order
        first_mol = True
        for mol in self:
            if first_mol:
                order = None
                try:
                    order = self.order
                except AttributeError:
                    pass
                if order:
                    mol.SetProp("order", ";".join(order))

            try:
                mol.GetConformer()
            except ValueError:  # no 2D coords... calculate them
                mol.Compute2DCoords()

            writer.write(mol, confId=conf_id)

            # remove the order property again from mol_list
            if first_mol:
                first_mol = False
                mol.ClearProp("order")

        writer.close()


    def write_csv(self, fn="mols.csv", props=None, include_smiles=True, isomeric=True):
        """Writes the Mol_List as a csv file to disk.

        Parameters:
            fn (str): Filename.
            props (list[string]): An optional list of molecule properties to write.
                If `props` is None, all props are written.
            include_smiles (bool): If true, the Smiles will be calculated on the fly
                and written to the csv.
            isomeric (bool): If True, the generated Smiles will be isomeric."""

        if props is None:
            props = self.fields
        if not isinstance(props, list):
            props = [props]
        csv_fields = props.copy()
        if include_smiles:
            csv_fields.append("Smiles")
        with open(fn, "w") as f:
            writer = csv.DictWriter(f, csv_fields, dialect="excel-tab")
            writer.writeheader()
            for mol in self:
                row = {}
                if include_smiles:
                    smi = Chem.MolToSmiles(mol, isomericSmiles=isomeric)
                    row["Smiles"] = smi
                for prop in props:
                    if mol.HasProp(prop):
                        val = mol.GetProp(prop)
                        if val != "":
                            row[prop] = val
                writer.writerow(row)


    def sort_list(self, field, reverse=True):
        """Sort the Mol_List according to <field>."""
        self.sort(key=lambda x: self._key_get_prop(x, field, reverse=reverse), reverse=reverse)


    def order_props(self, order="default"):
        """Arrange the display order of the properties."""
        order_props(self, order)


    def mols_with_prop(self, prop):
        """Returns:
            Am iterator of molecules in the list where mol and prop are defined."""

        for mol in self:
            if mol and mol.HasProp(prop):
                yield mol


    def prop_filter(self, query, invert=False, sorted=True, reverse=True, field_types=None, make_copy=True):
        """Return a new Mol_List based on the property filtering.
        By default it creates an independent copy of the mol objects."""
        result_list = Mol_List()
        if self.order:
            result_list.order = self.order.copy()
        result_list.ia = self.ia

        mol_counter_out = 0

        if not field_types:
            field_types = self.field_types

        if not field_types:
            print("  # no field type information available! -aborted.")
            return None

        field = None
        for el in query.split(" "):
            if el in field_types:
                field = el
                break

        if not field:
            print("  # field could not be extracted from query! -aborted.")
            return None

        print("  > field {} extracted from query: {}.".format(field, query))

        query_mod = query.replace(field, "val")

        for mol_counter_in, mol in enumerate(self):
            if not mol:
                continue
            hit = False
            if field in mol.GetPropNames():
                val = mol.GetProp(field).lower()
                if field_types[field] in ["number", "key"]:
                    try:
                        val_float = float(val)

                    except ValueError:
                        continue

                    val_int = int(val_float)
                    if val_int == val_float:
                        val = val_int
                    else:
                        val = val_float

                hit = eval(query_mod)

                if invert:
                    hit = not hit

                if hit:
                    mol_counter_out += 1
                    if make_copy:
                        mol = deepcopy(mol)
                    result_list.append(mol)

        print("  > processed: {:7d}   found: {:6d}".format(mol_counter_in + 1, mol_counter_out))

        if sorted:
            result_list.sort_list(field, reverse=reverse)

        return result_list


    def mol_filter(self, query, smarts=False, invert=False,
                   align=None, add_h=False, make_copy=True):
        """Returns a new Mol_List containing the substructure matches.
        By default it creates an independent copy of the mol objects."""
        result_list = Mol_List()
        if self.order:
            result_list.order = self.order.copy()
        result_list.ia = self.ia

        mol_counter_out = 0
        if isinstance(query, str):
            if "[H]" in query or "#1" in query:
                add_h = True
                print("> explicit hydrogens turned on (add_h = True)")

            if add_h or "#6" in query or "#7" in query:
                smarts = True

            if smarts:
                query_mol = Chem.MolFromSmarts(query)
                if align is None:  # Aligning to mol generated from Smarts does not work
                    align = False
            else:
                Chem.MolFromSmiles(query)
                if align is None:
                    align = True

        else:
            query_mol = query
            if align is None:
                atm = query_mol.GetAtomWithIdx(1)
                if atm.HasQuery():  # True for molecules that were generated from Smarts
                    align = False   # Aligning to mol generated from Smarts does not work
                else:
                    align = True

        if not query_mol:
            print("* ERROR: could not generate query molecule. Try smarts=True")
            return None

        for mol_counter_in, mol in enumerate(self):
            if not mol: continue

            hit = False
            if add_h:
                mol_with_h = Chem.AddHs(mol)
                if mol_with_h.HasSubstructMatch(query_mol):
                    hit = True

            else:
                if mol.HasSubstructMatch(query_mol):
                    hit = True

            if invert:
                # reverse logic
                hit = not hit

            if hit:
                mol_counter_out += 1
                if make_copy:
                    mol = deepcopy(mol)
                result_list.append(mol)

        if align and len(result_list) > 0:
            result_list.align(query_mol)

        print("> processed: {:7d}   found: {:6d}".format(mol_counter_in + 1, mol_counter_out))

        return result_list


    def has_prop_filter(self, prop, invert=False, make_copy=True):
        """Returns a new Mol_list with molecules containing the property `prop`.
        By default it creates an independent copy of the mol objects."""

        result_list = Mol_List()
        if self.order:
            result_list.order = self.order.copy()
        result_list.ia = self.ia

        mol_counter_out = 0
        for mol_counter_in, mol in enumerate(self):
            if not mol: continue

            hit = False
            if mol.HasProp(prop):
                hit = True

            if invert:
                hit = not hit

            if hit:
                mol_counter_out += 1
                if make_copy:
                    mol = deepcopy(mol)
                result_list.append(mol)

        print("> processed: {:7d}   found: {:6d}".format(mol_counter_in + 1, mol_counter_out))

        return result_list



    def get_ids(self):
        """Get the list of Compound IDs in the Mol_List

        Parameters:
            id_prop (None, str): (optional) The name of the id_prop, if None, it will be guessed."

        Returns:
            A list of compound ids"""
        prop_list = self.fields

        if self.id_prop:
            if self.id_prop not in prop_list:
                raise LookupError("id_prop not found in data set.")
        else:  # try to guess an id_prop
            self.id_prop = guess_id_prop(prop_list)

        if self.id_prop is None:
            raise LookupError("no id prop could be found in data set.")

        id_list = []
        for mol in self:
            if mol:
                if mol.HasProp(self.id_prop):
                    val = get_value(mol.GetProp(self.id_prop))
                    id_list.append(val)

        return id_list


    def new_list_from_ids(self, id_list, invert=False, make_copy=True):
        """Creates a new Mol_List out of the given IDs.

        Parameters:
            id_list (list): The list of IDs
            id_prop (None, str): (optional) The name of the id_prop, if None, it will be guessed.

        Returns:
            A new Mol_List from a list of Ids.
            By default it creates an independent copy of the mol objects."""

        if not isinstance(id_list, list):
            id_list = [id_list]

        id_all = set(self.get_ids())
        id_set = set(id_list)
        if invert:
            id_keep = id_all - id_set
        else:
            id_keep = id_set.intersection(id_all)

        new_list = Mol_List()
        if self.order:
            new_list.order = self.order.copy()
        new_list.ia = self.ia

        for mol in self:
            if mol:
                if mol.HasProp(self.id_prop):
                    val = get_value(mol.GetProp(self.id_prop))
                    if val in id_keep:
                        if make_copy:
                            mol = deepcopy(mol)

                        new_list.append(mol)

        return new_list


    def show_cpd(self, id_no, is_cpd_id=True, make_copy=True, show_smiles=True):
        """Display a single compound together with its Smiles.
        With is_cpd_id == True (default), the given id_no is interpreted as a Compound_Id.
        Otherwise it is used as index in the list."""

        new_list = Mol_List()
        if self.order:
            new_list.order = self.order.copy()
        new_list.ia = self.ia
        new_list.id_prop = self.id_prop

        if not is_cpd_id:
            idx = id_no
            if make_copy:
                mol = deepcopy(self[id_no])
            else:
                mol = self[id_no]
            new_list.append(mol)

        else:
            if self.id_prop is None:
                self.id_prop = guess_id_prop(self.fields)
                if self.id_prop is None:
                    raise LookupError("Id property {} could not be found in the Mol_List.".format(self.id_prop))

            for idx, mol in enumerate(self):
                if mol:
                    if mol.HasProp(self.id_prop):
                        val = get_value(mol.GetProp(self.id_prop))
                        if val == id_no:
                            if make_copy:
                                mol = deepcopy(mol)
                            new_list.append(mol)

        if len(new_list) == 0:
            raise LookupError("no molecule with {}: {} could be found in the Mol_List.".format(self.id_prop, id_no))

        if show_smiles:
            print("idx: {:3d}   Smiles: {}".format(idx, Chem.MolToSmiles(new_list[0])))
        return new_list


    def set_prop_on_mol(self, id_no, prop_name, prop_value, is_cpd_id=True):
        """Change the value of a property in the Mol_List.
        prop_name (str) is the name of the property,
        prop_value (str) the value to which it will be set (using mol.SetProp()).
        With is_cpd_id == True (default), the given id_no is interpreted as a Compound_Id.
        Otherwise it is used as index in the list."""

        mol = self.show_cpd(id_no, is_cpd_id=is_cpd_id, make_copy=False, show_smiles=False)[0]
        mol.SetProp(prop_name, prop_value)


    def calc_props(self, props, force2d=False, **kwargs):
        """Calculate properties from the Mol_List.
        props can be a single property or a list of properties.

        Calculable properties:
            2d, date, formula, smiles, hba, hbd, logp, molid, mw, rotb,
            sa (synthetic accessibility, tpsa, murcko (MurckoScaffold as Smiles),
            sim (similarity relative to `sim_mol_or_smiles` or the mol with `sim_id`),
            smiles (isomeric=True/False)

        Synthetic Accessibility (normalized):
            0: hard to synthesize; 1: easy access

            as described in:
                | Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions
                | *Peter Ertl and Ansgar Schuffenhauer*
                | Journal of Cheminformatics 1:8 (2009) (`link <http://www.jcheminf.com/content/1/1/8>`_)
        """
        sim_mol_or_smiles = kwargs.get("sim_mol_or_smiles", None)
        sim_id = kwargs.get("sim_id", None)
        query_fp = None

        if not isinstance(props, list):
            props = [props]

        # make all props lower-case:
        props = list(map(lambda x: x.lower(), props))

        if sim_id is not None:  # sim_id represents a Compound_Id,
                                # which is then taken as the Similarity base
            sim_mol_or_smiles = self.show_cpd(sim_id, is_cpd_id=True,
                                              make_copy=False, show_smiles=False)[0]

        if sim_mol_or_smiles is not None:
            if isinstance(sim_mol_or_smiles, str):
                sim_mol_or_smiles = Chem.MolFromSmiles(sim_mol_or_smiles)
            if USE_AVALON:
                query_fp = pyAv.GetAvalonFP(sim_mol_or_smiles, 1024)
            else:
                query_fp = FingerprintMols.FingerprintMol(sim_mol_or_smiles)

        ctr = 0
        calculated_props = set()
        for mol in self:
            if not mol: continue

            if "molid" in props:
                ctr += 1
                mol.SetProp("Mol_Id", str(ctr))
                calculated_props.add("molid")

            calc_props(mol, props, force2d=force2d, query_fp=query_fp,
                       calculated_props=calculated_props, **kwargs)

        self._set_recalc_needed()

        not_calculated = set(props) - calculated_props
        if not_calculated:
            print("* these props could not be calculated:", not_calculated)


    def remove_props(self, props):
        """Remove properties from the Mol_List.
        props can be a single property or a list of properties."""

        for mol in self:
            if mol:
                remove_props_from_mol(mol, props)

        self._set_recalc_needed()


    def remove_empty_props(self):
        remove_empty_props(self)
        self._set_recalc_needed()


    def keep_props(self, props):
        """Keep properties in the Mol_List.
        props can be a single property or a list of properties."""

        if not isinstance(props, list):
            props = [props]

        for mol in self:
            if mol:
                keep_props_in_mol(mol, props)

        self.order = props.copy()
        self._set_recalc_needed()


    def copy_prop(self, prop_orig, prop_copy, move=False):
        """Copy or rename a property in the Mol_List."""

        for mol in self.mols_with_prop(prop_orig):
            val_orig = mol.GetProp(prop_orig)
            mol.SetProp(prop_copy, val_orig)
            if move:
                mol.ClearProp(prop_orig)

        self._set_recalc_needed()


    def rename_prop(self, prop_orig, prop_new):
        """Convenience wrapper around copy_prop"""

        self.copy_prop(prop_orig, prop_new, move=True)


    def remove_dups_by_id(self, id_prop=None, make_copy=True):
        """Remove duplicate records by Compound Id.

        Parameters:
            id_prop (None, str): The name of the Id property, if *None*, it will be guessed.

        Returns:
            new Mol_list without the duplicate Ids.
            By default it creates an independent copy of the mol objects."""

        new_list = Mol_List()
        if self.order:
            new_list.order = self.order.copy()
        new_list.ia = self.ia

        id_list = []
        if not id_prop:
            id_prop = guess_id_prop(list_fields(self))
        if not id_prop:
            print("* could not determine Id property.")
            return None

        for mol in self:
            if not mol: continue
            mol_id = mol.GetProp(id_prop)
            if mol_id in id_list: continue
            id_list.append(mol_id)
            if make_copy:
                mol = deepcopy(mol)

            new_list.append(mol)

        return new_list


    def remove_dups_by_struct(self, make_copy=True):
        """Remove duplicates by structure. Duplicates are determined by Smiles.

        Returns:
            new Mol_List without the duplicate structures.
            By default it creates an independent copy of the mol objects. """

        new_list = Mol_List()
        if self.order:
            new_list.order = self.order.copy()
        new_list.ia = self.ia

        smiles_list = []
        for mol in self:
            if not mol: continue
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)  # needed to distinguish between stereoisomers
            if smiles in smiles_list: continue
            smiles_list.append(smiles)
            if make_copy:
                mol = deepcopy(mol)

            new_list.append(mol)

        return new_list


    def join_data_from_file(self, fn, id_prop=None, decimals=2):
        """Joins data from a file with name ``fn`` by Id property ``id_prop``. If no Id property is given, it will be guessed.
        CAUTION: The records from the file are loaded into memory!

        Parameters:
            decimals (int): number of decimal places for floating point values."""

        if not id_prop:
            id_prop = guess_id_prop(self.field_types)

        file_d = {}
        for line in csv_supplier(fn):
            rec_id = get_value(line.pop(id_prop))
            file_d[rec_id] = line

        for mol in self:
            mol_id = get_prop_val(mol, id_prop)
            if mol_id in file_d:
                records = file_d[mol_id]
                for rec in records:
                    val = get_value(records[rec])
                    if val is None: continue

                    if isinstance(val, float):
                        mol.SetProp(rec, "{val:.{decimals}f}".format(val=val, decimals=decimals))
                    else:
                        mol.SetProp(rec, str(val))

        self._set_recalc_needed()


    def set_default(self, prop, def_val, condition=None):
        """Set a default value in all mols, in which ``prop`` is either not defined (``condition`` == None) or
        is evaluating ``condition`` to true."""

        failed = 0
        if condition and not isinstance(condition, str):
            raise TypeError("condition needs to be of type str.")

        for mol in self:
            if not mol: continue
            if not condition:
                if not mol.HasProp(prop):
                    mol.SetProp(prop, str(def_val))
            else:
                if mol.HasProp(prop):
                    prop_val = get_value(mol.GetProp(prop))
                    if isinstance(prop_val, str):
                        eval_templ = """'{}' {}"""
                    else:
                        eval_templ = """{} {}"""

                    try:
                        if eval(eval_templ.format(prop_val, condition)):
                            mol.SetProp(prop, str(def_val))

                    except SyntaxError:
                        failed += 1

        self.recalc_needed["d"] = True
        self.recalc_needed["field_types"] = True

        if failed > 0:
            print("# {} records could not be processed.".format(failed))


    def table(self, pagesize=25, highlight=None, show_hidden=False, img_dir=None, raw=False):
        """Return the Mol_List as HTML table.
        Either as raw HTML (raw==True) or as HTML object for display in IPython notebook.

        Parameters:
            show_hidden (bool): Whether to show hidden properties (name starts with _) or not.
                Default is False.
            raw (bool): If True, return the HTML mol grid as text.
                If False, return a HTML object, that can be displayed in the Jupyter Notebook.
                Default is False.
            img_dir (str or None): The directory, in which the molecule images are written. The directory has to exist.
                Implies raw=True. If None, then the images are stored in the HTML object. Default is None."""

        if self.id_prop is None:
            self.id_prop = guess_id_prop(list_fields(self))

        if img_dir is not None:
            raw = True

        if raw:
            return mol_table(self, id_prop=self.id_prop, highlight=highlight, interact=self.ia,
                             order=self.order, img_dir=img_dir, show_hidden=show_hidden)
        else:
            return table_pager(self, pagesize=pagesize, id_prop=self.id_prop, interact=self.ia, highlight=highlight, order=self.order,
                               show_hidden=show_hidden)


    def nested(self, pagesize=10, props=None, img_dir=None, raw=False):
        if self.id_prop is None:
            self.id_prop = guess_id_prop(list_fields(self))

        if img_dir is not None:
            raw = True

        if raw:
            return nested_table(self, id_prop=self.id_prop, props=props, order=self.order, img_dir=img_dir)
        else:
            return nested_pager(self, pagesize=pagesize, id_prop=self.id_prop, props=props, order=self.order)



    def grid(self, pagesize=12, props=None, highlight=None, mols_per_row=4, size=IMG_GRID_SIZE, img_dir=None, raw=False):
        """Returns:
            The Mol_List as HTML grid table. Either as raw HTML (raw==True) or as HTML object for display in IPython notebook.

        Parameters:
            props: A property or a list of properties to include in the display.
            raw (bool): If True, return the HTML mol grid as text.
                If False, return a HTML object, that can be displayed in the Jupyter Notebook.
                Default is False.
            img_dir (str or None): The directory, in which the molecule images are written. The directory has to exist.
                Implies raw=True. If None, then the images are stored in the HTML object. Default is None."""

        if self.id_prop is None:
            self.id_prop = guess_id_prop(list_fields(self))

        if img_dir is not None:
            raw = True

        if raw:
            return mol_sheet(self, props=props, id_prop=self.id_prop, interact=self.ia,
                             highlight=highlight, mols_per_row=mols_per_row, size=size, img_dir=img_dir)
        else:
            return grid_pager(self, pagesize, props=props, id_prop=self.id_prop, interact=self.ia, highlight=highlight,
                              mols_per_row=mols_per_row, size=size)


    def write_table(self, highlight=None, header=None, summary=None, img_dir=None,
                    title="Results", fn="mol_table.html"):
        html.write(html.page(self.table(highlight=highlight, raw=True, img_dir=img_dir),
                             header=header, summary=summary, title=title), fn=fn)
        return HTML('<a href="{}">{}</a>'.format(fn, fn))


    def write_nested(self, header=None, summary=None, img_dir=None, title="Results", fn="nested_table.html"):
        html.write(html.page(self.nested(raw=True, img_dir=img_dir),
                             header=header, summary=summary, title=title), fn=fn)
        return HTML('<a href="{}">{}</a>'.format(fn, fn))


    def write_grid(self, props=None, highlight=None, mols_per_row=5, size=IMG_GRID_SIZE,
                   header=None, summary=None, img_dir=None, title="Results", fn="mol_grid.html"):
        html.write(html.page(self.grid(props=props, highlight=highlight,
                             mols_per_row=mols_per_row, size=size, img_dir=img_dir, raw=True), header=header, summary=summary, title=title), fn=fn)
        return HTML('<a href="{}">{}</a>'.format(fn, fn))


    def scatter(self, x, y, r=7, tooltip=None, **kwargs):
        """Displays a Highcharts plot in the IPython Notebook.
        Uses Bokeh (preferred) or the Highcharts javascript library, either locally under lib/ relative to the Notebook
        or the web version at http://code.highcharts.com.
        If ``tooltip`` is *None*, then structure tooltips will be shown for Mol_Lists with
        less than or equal 150 records, if the Mol_List has more records, no structure tooltips
        will be shown. The bevaviour can be forced by either providing ``tooltip="struct"`` for tooltips
        or ``tooltip=""`` for no tooltips. Properties in the ``jitter`` list (only when used for x or y)
        will be jittered by a magnitude of ``mag``.
        callback (str): clicking on a point will link to the given HTML address. `@<IdProperty>` can be used as placeholder for the point id (e.g. Compound_Id). Default is None."""

        if tooltip is None:
            if len(self) > 150:
                tooltip = ""
            else:
                tooltip = "struct"

        if self.plot_tool == "bokeh":
            return bkt.cpd_scatter(self.d, x, y, r=r, pid=self.id_prop, tooltip=tooltip, **kwargs)
        else:
            return hct.cpd_scatter(self.d, x, y, r=r, pid=self.id_prop, tooltip=tooltip, **kwargs)


    def hist(self, field, bins=10, title="Distribution", xlabel=None, ylabel="Occurrence", normed=False, show=True, **kwargs):
        """Displays a Bokeh histogram. See bokeh_tools for documentation.
        Possible useful additional kwargs include: plot_width, plot_height, y_axis_type="log"."""
        if xlabel is None:
            xlabel = field
        hist = bkt.Hist(title=title, xlabel=xlabel, ylabel=ylabel, **kwargs)
        hist.add_data(self.d[field], bins=bins, normed=normed)

        if show:
            return hist.show()
        else:
            return hist.plot


    def bar(self, x, show=True, **kwargs):
        """Displays a bar chart for the occurrence of the given x-value.
        This plot type is especially useful for plotting the occurrence of categorical data,
        where only a small number (<= 10) of different values are present.
        This function is directly calling the advanced bokeh bar chart type,
        therefore no additional class is used.
        Useful kwargs include: title, plot_height, plot_width."""
        if show:
            bkt.bar_chart(self.d, x, show=True, **kwargs)
        else:
            return bkt.bar_chart(self.d, x, show=False, **kwargs)


    def summary(self, text_only=False):
        """Output a summary of the Mol_List and its properties.
        If ``text_only``is True only a text version is printed.
        ``mean`` and ``median`` are calculated with numpy."""

        field_types = self.field_types
        l = len(self)
        max_max = 0
        sum_d = {}
        for prop in field_types:
            value_list = [get_value(mol.GetProp(prop)) for mol in self.mols_with_prop(prop)]
            num_val = len(value_list)
            sum_d[prop] = {"num_values": num_val}
            sum_d[prop]["type"] = field_types[prop]
            if field_types[prop] == "number":
                sum_d[prop]["min"] = min(value_list)
                sum_d[prop]["max"] = max(value_list)
                if sum_d[prop]["max"] > max_max:
                    max_max = sum_d[prop]["max"]
                sum_d[prop]["mean"] = np.mean(value_list)
                sum_d[prop]["median"] = np.median(value_list)

        n_digits = str(np.floor(np.log10(max_max)) + 5.3) + "f"  # digits for formatting

        if text_only:
            print("number of records:", l)
            for prop in sum_d:
                print("\n{} ({}, {}):".format(prop, sum_d[prop]["type"], sum_d[prop]["num_values"]))
                if field_types[prop] == "number":
                    for sum_item in ["min", "max", "mean", "median"]:
                        print("{:6s}: {:{digits}}".format(sum_item, sum_d[prop][sum_item], digits=n_digits), end="   |   ")
                    print()

        else:
            rows = []
            cells = []
            opt1 = {"align": "center", "bgcolor": "#94CAEF"}
            opt2 = {"align": "center", "bgcolor": "#94CAEF", "colspan": 7}
            cell = html.td(html.b("Summary ({} records)".format(l)), options=opt2)
            rows.extend(html.tr(cell))
            for cell in ["Property", "Type", "Num Values", "Min", "Max", "Mean", "Median"]:
                cells.extend(html.td(html.b(cell), options=opt1))
            rows.extend(html.tr(cells))
            opt1 = {"align": "center"}
            for prop in sum_d:
                cells = []
                cells.extend(html.td(prop, options=opt1))
                cells.extend(html.td(sum_d[prop]["type"], options=opt1))
                cells.extend(html.td(str(sum_d[prop]["num_values"]), options=opt1))
                if field_types[prop] == "number":
                    for sum_item in ["min", "max", "mean", "median"]:
                        cells.extend(html.td("{:.3f}".format(sum_d[prop][sum_item]), options=opt1))
                else:
                    for i in range(4):  # insert empty cells
                        cells.extend(html.td("", options=opt1))
                rows.extend(html.tr(cells))

            table = html.table(rows)
            return HTML("".join(table))


    def correlate(self, min_corr=0.4, text_only=False):
        """Display correlations between the properties in the Mol_List.
        Calculated by np.corrcoef, only abs. values are used, higher value means higer correlation.
        Only correlations greater or to equal to ``min_corr`` are shown (default=0.4).
        If ``text_only`` is True only a text version is printed."""

        number_fields = [f for f in self.field_types if self.field_types[f] == "number"]
        n = len(number_fields)
        pair_format = str(max(len(i) for i in number_fields) * 2 + 7) + "s"
        corr_d = {}
        l = len(self)
        for left in range(n):
            left_values = [get_prop_val(mol, number_fields[left]) for mol in self]
            for right in range(left + 1, n):
                right_values = [get_prop_val(mol, number_fields[right]) for mol in self]
                both_y = []
                both_x = []
                for i in range(l):
                    if left_values[i] is None or right_values[i] is None:
                        continue
                    both_y.append(left_values[i])
                    both_x.append(right_values[i])
                corr = np.corrcoef(both_y, both_x)
                corr_val = abs(corr[0][1])
                if corr_val >= min_corr:
                    k = "{} vs. {}".format(number_fields[left], number_fields[right])
                    corr_d[k] = corr_val

        if text_only:
            print("Property Correlation Coefficients:")
            for pair in sorted(corr_d, key=corr_d.get, reverse=True):
                print("{pair:{pair_format}}: {corr:.3f}".format(pair=pair,
                      pair_format=pair_format, corr=corr_d[pair]))

        else:
            rows = []
            cells = []
            opt1 = {"align": "center", "bgcolor": "#94CAEF"}
            opt2 = {"align": "center", "bgcolor": "#94CAEF", "colspan": 2}
            opt3 = {"align": "center", "bgcolor": "#94CAEF", "colspan": 3}
            cell = html.td(html.b("Property Correlation Coefficients"), options=opt3)
            rows.extend(html.tr(cell))
            cells.extend(html.td(html.b("A vs. B"), options=opt2))
            cells.extend(html.td(html.b("Correlation"), options=opt1))
            rows.extend(html.tr(cells))

            opt1 = {"align": "center"}
            for pair in sorted(corr_d, key=corr_d.get, reverse=True):
                cells = []
                cells.extend(html.td(pair.split(" vs. ")[0], options=opt1))
                cells.extend(html.td(pair.split(" vs. ")[1], options=opt1))
                cells.extend(html.td("{:.3f}".format(corr_d[pair]), options=opt1))
                rows.extend(html.tr(cells))

            table = html.table(rows)
            return HTML("".join(table))


    @property
    def fields(self):
        """A List of properties that are present in the Mol_List (property)."""

        if self.len != len(self):
            self._set_recalc_needed()
        if self.recalc_needed["fields"]:
            self._fields = list_fields(self)
            self.recalc_needed["fields"] = False
        return self._fields


    @property
    def field_types(self):
        """A dictionary of properties and their derived types (property)."""

        if self.len != len(self):
            self._set_recalc_needed()
        if self.recalc_needed["field_types"]:
            self._field_types = get_field_types(self)
            self.recalc_needed["field_types"] = False
        return self._field_types


    @property
    def d(self):
        """Representation of the Mol_List as a dictionary for plotting (property)."""

        if self.len != len(self):
            self._set_recalc_needed()
        if self.recalc_needed["d"] or self.plot_tool != self.recalc_needed["plot_tool"]:
            self._calc_d()
            self.recalc_needed["d"] = False
            self.recalc_needed["plot_tool"] = self.plot_tool
        return self._d



def create_dir_if_not_exist(dir_name):
    if not op.exists(dir_name):
        print("  * target folder does not exist, creating {}...".format(dir_name))
        os.makedirs(dir_name)


def autocrop(im, bgcolor="white"):
    if im.mode != "RGB":
        im = im.convert("RGB")
    bg = Image.new("RGB", im.size, bgcolor)
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)
    return None  # no contents


def list_fields(sdf_list):
    field_list = []

    for mol in sdf_list:
        field_list.extend(mol.GetPropNames())

    return list(set(field_list))


def load_sdf(file_name_or_obj="testset.sdf", order="default"):
    """Create a Mol_List instance from an SD File.
    Accepts a string filename or a file object as input.
    order: "default" or None."""

    if isinstance(file_name_or_obj, str):
        if PY3:
            file_obj = open(file_name_or_obj, "rb")
        else:
            file_obj = open(file_name_or_obj)
    else:
        file_obj = file_name_or_obj

    reader = Chem.ForwardSDMolSupplier(file_obj)

    sdf_list = Mol_List()

    # try to load the column order
    first_mol = True
    for mol in reader:
        if mol:
            if first_mol:
                first_mol = False
                prop_order = None
                try:
                    prop_order = mol.GetProp("order")
                    remove_props_from_mol(mol, "order")
                except KeyError:  # first mol does not contain an order field
                    pass

                if prop_order is not None:
                    try:
                        sdf_list.order = prop_order.split(";")

                    except AttributeError:  # sdf_list is not a Mol_List
                        pass

            sdf_list.append(mol)

    if sdf_list.id_prop is None:
        sdf_list.id_prop = guess_id_prop(sdf_list.fields)
    if sdf_list.order is None and order == "default":  # only when no order is already present.
        sdf_list.order_props(order=order)

    if isinstance(file_name_or_obj, str):
        print("  > sdf {} loaded with {} records.".format(file_name_or_obj.split(".")[0], len(sdf_list)))
    else:
        print("  > sdf loaded with {} records.".format(len(sdf_list)))

    return sdf_list


def load_csv(fn, smiles_col="Smiles"):
    """Reads a csv file and returns a Mol_List instance. The molecules are generated from the Smiles column, which has to be present."""
    with open(fn) as f:
        ctr = 0
        sdf_list = Mol_List()
        reader = csv.DictReader(f, dialect="excel-tab")
        for row_dict in reader:
            if smiles_col not in row_dict: continue
            smi = row_dict.pop("Smiles", "")
            mol = Chem.MolFromSmiles(smi)
            if not mol: continue
            for prop in row_dict:
                val = row_dict[prop]
                if val != "":
                    mol.SetProp(prop, val)
            sdf_list.append(mol)
            ctr += 1

    print("> {} loaded into Mol_List ({} records).".format(fn, ctr))
    return sdf_list


def order_props(sdf_list, order="default"):
    """Order fields. First Compound_Id, Supplier, Producer;
    then the activity fields, then the physicochemical properties and LCMS"""
    if order == "default":
        prop_order = []
        fields = sorted(sdf_list.fields)
        for def_ord in DEFAULT_ORDER:
            def_ord_items = def_ord.split("|")
            fields_found = []
            for f in fields:
                for item in def_ord_items:
                    if item in f.lower():
                        prop_order.append(f)
                        fields_found.append(f)
            # Remove the fields that are now already on the prop_order list
            # from the original field list.
            # This way they will not be added multiple times
            for f in fields_found:
                fields.remove(f)

        if len(fields) > 0:  # add the remaining fields in alphabetical order
            prop_order.extend(fields)
        sdf_list.order = prop_order


def write_ids(id_list, fn="id_list.txt"):
    """Write a list of compound ids to a file. The list will be sorted by Id."""
    id_str = "\n".join(sorted([str(i) for i in id_list]))
    id_str = "Compound_Id\n" + id_str
    f = open(fn, "w")
    f.write(id_str)
    f.close()


def csv_supplier(fn):
    """Returns a dictionary generator."""

    if ".gz" in fn:
        f = gzip.open(fn, mode="rt")
    else:
        f = open(fn)
    reader = csv.DictReader(f, dialect="excel-tab")
    for row_dict in reader:
        yield row_dict

    f.close()


def keep_props_in_mol(mol, prop_or_propslist):
    if not isinstance(prop_or_propslist, list):
        prop_or_propslist = [prop_or_propslist]

    mol_props = mol.GetPropNames()
    for prop in mol_props:
        if prop not in prop_or_propslist:
            mol.ClearProp(prop)


def remove_props_from_mol(mol, prop_or_propslist):
    if not isinstance(prop_or_propslist, list):
        prop_or_propslist = [prop_or_propslist]
    for prop in prop_or_propslist:
        if prop in mol.GetPropNames():
            mol.ClearProp(prop)


def remove_props(mol_or_sdf_list, props):
    if isinstance(mol_or_sdf_list, list):
        for mol in mol_or_sdf_list:
            if mol:
                remove_props_from_mol(mol, props)
    else:
        remove_props_from_mol(mol_or_sdf_list, props)


def remove_empty_props(mol_list):
    for mol in mol_list:
        props = mol.GetPropNames()
        for prop in props:
            if mol.GetProp(prop) == "":
                mol.ClearProp(prop)


def unit_factor(unit):
    """Return the factor corresponding to the unit, e.g. 1E-9 for nM.
    Known units are: mM, uM, nM, pM. Raises ValueError for unknown unit."""
    units = ["mm", "um", "nm", "pm"]
    pos = units.index(unit.lower()) + 1
    factor = 10 ** -(pos * 3)
    return factor


def pic50(ic50, unit=None, ndigits=2):
    """Calculate pIC50 from IC50. Optionally, a unit for the input IC50 value may be given.
    Known units are: mM, uM, nM, pM"""
    if unit is not None:
        ic50 *= unit_factor(unit)
    return round(-math.log10(ic50), ndigits=ndigits)


def ic50(pic50, unit=None):
    """Calculate IC50 from pIC50. Optionally, a unit for the returned IC50 value may be given.
    Known units are: mM, uM, nM, pM"""
    ic50 = 10 ** (-pic50)
    if unit is not None:
        ic50 /= unit_factor(unit)
    return ic50


def set_margin(container, margin=10):
    """Recursively set margins on all widgets of a toplevel container (...Box())."""
    if hasattr(container, "children"):
        for ch in container.children:
            set_margin(ch, margin)
    else:
        container.margin = margin


def ia_remove_props(mol_list):
    """Interactively remove properties from a Mol_List.
    Uses IPython widgets to display the properties to be selected for removal."""

    all_props = list_fields(mol_list)

    def on_btn_clicked(b):
        remove_props(mol_list, props=list(w_sm.selected_labels))

    w_sm = ipyw.SelectMultiple(description="Properties to remove:", options=all_props)
    w_btn = ipyw.Button(description="Done !")
    w_btn.on_click(on_btn_clicked)

    w_hb = ipyw.HBox(children=[w_sm, w_btn])

    display(w_hb)


def ia_keep_props(mol_list):
    """Interactively keep properties from a Mol_List.
    Uses IPython widgets to display the properties to be selected for keeping."""

    all_props = list_fields(mol_list)

    def on_btn_clicked():
        props_to_remove = list(set(all_props) - set(w_sm.selected_labels))
        remove_props(mol_list, props=props_to_remove)

    w_sm = ipyw.SelectMultiple(description="Properties to keep:", options=all_props)
    w_btn = ipyw.Button(description="Done !")
    w_btn.on_click(on_btn_clicked)

    w_hb = ipyw.HBox(children=[w_sm, w_btn])

    display(w_hb)


def ia_smiles_from_smiles():
    """Sounds silly, but generates RDKit Smiles out of Smiles that were generated by other tools.
    May still be silly..."""
    smiles = input("Smiles: ")
    mol = Chem.MolFromSmiles(smiles)
    print(Chem.MolToSmiles(mol))


def check_2d_coords(mol, force=False):
    """Check if a mol has 2D coordinates and if not, calculate them."""
    try:
        mol.GetConformer()
    except ValueError:
        force = True  # no 2D coords... calculate them

    if force:
        if USE_AVALON:
            pyAv.Generate2DCoords(mol)
        else:
            mol.Compute2DCoords()


def calc_props(mol, props, force2d=False, calculated_props=None, **kwargs):
    """calculated_props can be None or of type set()."""

    sim_mol_or_smiles = kwargs.get("sim_mol_or_smiles", None)
    isomeric = kwargs.get("isomeric", True)
    query_fp = kwargs.get("query_fp", None)

    if not isinstance(props, list):
        props = [props]

    for prop in props:
        if "2d" in props:
            check_2d_coords(mol, force2d)
            if calculated_props is not None:
                calculated_props.add("2d")

        if "date" in props:
            mol.SetProp("Date", time.strftime("%Y%m%d"))
            if calculated_props is not None:
                calculated_props.add("date")

        if "formula" in props:
            mol.SetProp("Formula", Chem.CalcMolFormula(mol))
            if calculated_props is not None:
                calculated_props.add("formula")

        if "smiles" in props:
            mol.SetProp("Smiles", Chem.MolToSmiles(mol, isomericSmiles=isomeric))
            calculated_props.add("smiles")

        if "hba" in props:
            mol.SetProp("HBA", str(Desc.NOCount(mol)))
            if calculated_props is not None:
                calculated_props.add("hba")

        if "hbd" in props:
            mol.SetProp("HBD", str(Desc.NHOHCount(mol)))
            if calculated_props is not None:
                calculated_props.add("hbd")

        if "logp" in props:
            mol.SetProp("LogP", "{:.2f}".format(Desc.MolLogP(mol)))
            if calculated_props is not None:
                calculated_props.add("logp")

        if "mw" in props:
            mol.SetProp("MW", "{:.2f}".format(Desc.MolWt(mol)))
            if calculated_props is not None:
                calculated_props.add("mw")

        if "rotb" in props:
            mol.SetProp("RotB", str(Desc.NumRotatableBonds(mol)))
            if calculated_props is not None:
                calculated_props.add("rotb")

        if SASCORER and "sa" in props:
            score = sascorer.calculateScore(mol)
            norm_score = 1 - (score / 10)
            mol.SetProp("SA", "{:.2f}".format(norm_score))
            if calculated_props is not None:
                calculated_props.add("sa")

        if "tpsa" in props:
            mol.SetProp("TPSA", str(int(Desc.TPSA(mol))))
            if calculated_props is not None:
                calculated_props.add("tpsa")

        if "murcko" in props:
            msmiles = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
            mol.SetProp("Murcko", msmiles)
            calculated_props.add("murcko")

        if "sim" in props:
            if sim_mol_or_smiles is not None:
                if isinstance(sim_mol_or_smiles, str):
                    sim_mol_or_smiles = Chem.MolFromSmiles(sim_mol_or_smiles)
                    if USE_AVALON:
                        query_fp = pyAv.GetAvalonFP(sim_mol_or_smiles, 1024)
                    else:
                        query_fp = FingerprintMols.FingerprintMol(sim_mol_or_smiles)

            if query_fp is not None:
                if USE_AVALON:
                    mol_fp = pyAv.GetAvalonFP(mol, 1024)
                else:
                    mol_fp = FingerprintMols.FingerprintMol(mol)
                sim = DataStructs.FingerprintSimilarity(query_fp, mol_fp)
                mol.SetProp("Sim", "{:.2f}".format(sim * 100))
                calculated_props.add("sim")


def find_mcs(mol_list):
    """Returns the MCS molecule object for a set of molecule or None if not found."""
    mcs = rdFMCS.FindMCS(mol_list, maximizeBonds=True, matchValences=True,
                         ringMatchesRingOnly=True, completeRingsOnly=True)

    if mcs.canceled:
        print("* MCSS function timed out. Please provide a mol_or_smiles to align to.")
        return None

    if mcs.smartsString:
        mol = Chem.MolFromSmarts(mcs.smartsString)
        if not mol:
            return None

        mol.UpdatePropertyCache(False)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_SYMMRINGS | Chem.SANITIZE_SETCONJUGATION | Chem.SANITIZE_SETHYBRIDIZATION)

        return mol

    else:
        print("* Could not find MCSS. Please provide a mol_or_smiles to align to.")
        return None


def align(mol_list, mol_or_smiles=None):
    """Align the Mol_list to the common substructure provided as Mol or Smiles.

    Parameters:
        mol_list: A list of RDKit molecules.
        mol_or_smiles (None, str, mol or list thereof): The substructure(s) to which to align.
            If None, then the method uses rdFMCS to determine the MCSS
            of the mol_list."""


    if mol_or_smiles is None:
        # determine the MCSS
        mol_or_smiles = find_mcs(mol_list)
        if mol_or_smiles is None:
            return

    if not isinstance(mol_or_smiles, list):
        mol_or_smiles = [mol_or_smiles]

    align_mols = []
    for el in mol_or_smiles:
        if isinstance(el, str):
            mol = Chem.MolFromSmiles(el)
            check_2d_coords(mol)
            align_mols.append(mol)
        else:
            mol = deepcopy(el)
            check_2d_coords(mol)
            align_mols.append(mol)


    for mol in mol_list:
        if mol:
            check_2d_coords(mol)

        for align_mol in align_mols:
                if mol.HasSubstructMatch(align_mol):
                    Chem.GenerateDepictionMatching2DStructure(mol, align_mol)
                    break


def guess_id_prop(prop_list):  # try to guess an id_prop
    for prop in prop_list:
        if prop.lower().endswith("id"):
            return prop
    return None


def get_field_types(mol_list):
    """Detect all the property field types and return as dict"""

    field_types = {}

    if len(mol_list) > 100:
        sdf_sample = random.sample(mol_list, len(mol_list) // 5)
    else:
        sdf_sample = mol_list

    for mol in sdf_sample:
        prop_names = mol.GetPropNames()

        for prop in prop_names:
            prop_type = "number"
            prop_str = mol.GetProp(prop)

            try:
                float(prop_str)
                if prop.lower().endswith("id"):
                    prop_type = "key"

            except ValueError:
                prop_type = "str"

            if prop in field_types:
                if field_types[prop] in ["number", "key"] and prop_type == "str":
                    # "str" overrides everything: if one string is among the values
                    # of a property, all values become type "str"
                    field_types[prop] = prop_type
            else:
                field_types[prop] = prop_type

    if not field_types:
        raise NoFieldTypes()

    return field_types


def get_value(str_val):
    if not str_val:
        return None

    try:
        val = float(str_val)
        if "." not in str_val:
            val = int(val)
    except ValueError:
        val = str_val

    return val


def isnumber(x):
    """Returns True, if x is a number (i.e. can be converted to float)."""
    try:
        float(x)
        return True
    except:
        return False


def get_prop_val(mol, prop, default=None):
    """Returns the value of the molecule's property or the default value, if it is not defined."""
    if mol.HasProp(prop):
        return get_value(mol.GetProp(prop))
    else:
        return default


def b64_img(mol, size=300):
    img_file = IO()
    img = autocrop(Draw.MolToImage(mol, size=(size, size)))
    img.save(img_file, format='PNG')

    b64 = base64.b64encode(img_file.getvalue())
    if PY3:
        b64 = b64.decode()
    img_file.close()

    return b64


def mol_table(sdf_list, id_prop=None, interact=False, highlight=None, show_hidden=False, order=None, img_dir=None, size=300):
    """Parameters:
        sdf_list (Mol_List): List of RDKit molecules
        highlight (dict): Dict of properties (special: *all*) and values to highlight cells,
            e.g. {"activity": "< 50"}
        show_hidden (bool): Whether to show hidden properties (name starts with _) or not.
            Defaults to *False*.
        link (str): column used for linking out
        target (str): column used as link target
        order (list): A list of substrings to match with the field names for ordering in the table header
        img_dir (str): if None, the molecule images are embedded in the HTML doc.
            Otherwise the images will be stored in img_dir and linked in the doc.

    Returns:
        HTML table as TEXT to embed in IPython or a web page."""

    time_stamp = time.strftime("%y%m%d%H%M%S")
    td_opt = {"style": "text-align: center;"}
    header_opt = {"bgcolor": "#94CAEF", "style": "text-align: center;"}
    table_list = []
    prop_list = list_fields(sdf_list)

    if isinstance(order, list):
        for k in reversed(order):
            prop_list.sort(key=lambda x: k.lower() in x.lower(), reverse=True)

    if id_prop is None:
        guessed_id = guess_id_prop(prop_list)
    else:
        guessed_id = id_prop

    if interact and guessed_id is not None:
        table_list.append(TBL_JAVASCRIPT.format(ts=time_stamp, bgcolor="transparent"))

    if id_prop is not None:
        if id_prop not in prop_list:
            raise LookupError("Id property {} not found in data set.".format(id_prop))

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

        if guessed_id:
            id_prop_val = mol.GetProp(guessed_id)
            img_id = id_prop_val
            cell_opt = {"id": "{}_{}".format(id_prop_val, time_stamp)}
        else:
            img_id = idx
            cell_opt = {"id": str(idx)}

        cell = html.td(str(idx), cell_opt)
        cells.extend(cell)

        if not mol:
            cells.extend(html.td("no structure"))

        else:
            if img_dir is None:  # embed the images in the doc
                b64 = b64_img(mol, size * 2)
                img_src = "data:image/png;base64,{}".format(b64)

            else:  # write them out to img_dir
                img_file = op.join(img_dir, "img_{}.png".format(img_id))
                img = autocrop(Draw.MolToImage(mol, size=(size * 2, size * 2)))
                img.save(img_file, format='PNG')
                img_src = img_file

            cell_opt = {}
            if interact and guessed_id is not None:
                img_opt = {"title": "Click to select / unselect",
                           "onclick": "toggleCpd('{}')".format(id_prop_val)}
            else:
                img_opt = {"title": str(img_id)}
            # img_opt["width"] = size
            # img_opt["height"] = size
            img_opt["style"] = 'max-width: {}px; max-height: {}px; display: block; margin: auto;'.format(size, size)

            cell = html.img(img_src, img_opt)
            cells.extend(html.td(cell, cell_opt))

        for prop in prop_list:
            td_opt = {"style": "text-align: center;"}
            if prop in mol_props:
                if not show_hidden and prop.startswith("_"): continue
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

    if interact and guessed_id is not None:
        table_list.append(ID_LIST.format(ts=time_stamp))

    # print(table_list)
    return "".join(table_list)


def mol_sheet(sdf_list, props=None, id_prop=None, interact=False, highlight=None, mols_per_row=4, size=IMG_GRID_SIZE, img_dir=None):
    """Creates a HTML grid out of the Mol_List input.

    Parameters:
        sdf_list (Mol_List): list of RDKit molecules
        highlight (dict): dict of properties (a.t.m only one) and values to highlight cells,
            e.g. {"activity": "< 50"}
        order (list): a list of substrings to match with the field names for ordering in the table header
        img_dir (str): if None, the molecule images are embedded in the HTML doc.
            Otherwise the images will be stored in img_dir and linked in the doc.

    Returns:
        HTML table as TEXT with molecules in grid-like layout to embed in IPython or a web page."""

    time_stamp = time.strftime("%y%m%d%H%M%S")
    prop_opt = {"style": "text-align: left;"}
    # td_opt = {"align": "center"}
    td_opt = {"style": "text-align: center;"}

    header_opt = {"bgcolor": BGCOLOR}
    table_list = []
    prop_list = list_fields(sdf_list)
    if props and not isinstance(props, list):
        props = [props]

    if id_prop is None:
        guessed_id = guess_id_prop(prop_list)
    else:
        guessed_id = id_prop

    if interact and guessed_id is not None:
        table_list.append(TBL_JAVASCRIPT.format(ts=time_stamp, bgcolor=BGCOLOR))

    if props is not None:
        td_opt["colspan"] = "2"
        prop_row_cells = {k: [] for k, _ in enumerate(props)}

    rows = []
    id_cells = []
    mol_cells = []
    for idx, mol in enumerate(sdf_list, 1):
        if guessed_id:
            id_prop_val = mol.GetProp(guessed_id)
            img_id = id_prop_val
            cell_opt = {"id": "{}_{}".format(id_prop_val, time_stamp)}
            cell_opt.update(td_opt)
            cell_opt.update(header_opt)
            id_cells.extend(html.td(id_prop_val, cell_opt))
        else:
            img_id = idx

        if not mol:
            cell = ["no structure"]

        else:
            if img_dir is None:  # embed the images in the doc
                b64 = b64_img(mol, size * 2)
                img_src = "data:image/png;base64,{}".format(b64)

            else:
                img_file = op.join(img_dir, "img_{}.png".format(img_id))
                img = autocrop(Draw.MolToImage(mol, size=(size * 2, size * 2)))
                img.save(img_file, format='PNG')
                img_src = img_file

            if interact and guessed_id is not None:
                img_opt = {"title": "Click to select / unselect",
                           "onclick": "toggleCpd('{}')".format(id_prop_val)}
            else:
                img_opt = {"title": str(img_id)}
            # img_opt["width"] = size
            # img_opt["height"] = size
            img_opt["style"] = 'max-width: {}px; max-height: {}px; display: block; margin: auto;'.format(size, size)

            cell = html.img(img_src, img_opt)

        # td_opt = {"align": "center"}
        td_opt = {"style": "text-align: center;"}
        if props is not None:
            td_opt["colspan"] = "2"

        if highlight:
            eval_str = None
            prop = highlight.keys()[0]  # only one highlight key supported a.t.m.
            prop_val = mol.GetProp(prop)
            eval_str = " ".join([prop_val, highlight[prop]])
            if eval_str and eval(eval_str):
                td_opt["bgcolor"] = "#99ff99"

        mol_cells.extend(html.td(cell, td_opt))

        if props:
            for prop_no, prop in enumerate(props):
                prop_cells = []
                prop_val = ""
                if mol.HasProp(prop):
                    prop_val = mol.GetProp(prop)
                prop_cells.extend(html.td(prop[:25], prop_opt))
                prop_cells.extend(html.td(prop_val[:8], prop_opt))
                prop_row_cells[prop_no].extend(prop_cells)

        if idx % mols_per_row == 0 or idx == len(sdf_list):
            if guessed_id:
                rows.extend(html.tr(id_cells))
            rows.extend(html.tr(mol_cells))

            if props is not None:
                for prop_no in sorted(prop_row_cells):
                    rows.extend(html.tr(prop_row_cells[prop_no]))
                prop_row_cells = {k: [] for k, _ in enumerate(props)}
            id_cells = []
            mol_cells = []

    table_list.extend(html.table(rows))

    if interact and guessed_id is not None:
        table_list.append(ID_LIST.format(ts=time_stamp))

    # print(table_list)
    return "".join(table_list)


def nested_table(mol_list, id_prop=None, props=None, order=None, size=300, img_dir=None):
    prop_list = list_fields(mol_list)

    if props is not None:
        if not isinstance(props, list):
            props = [props]
        order = props.copy()

    if order is None:
        order = ["Supplier", "Producer", "Hit", "ActAss"]


    if id_prop is None:
        guessed_id = guess_id_prop(prop_list)
    else:
        guessed_id = id_prop

    if guessed_id is not None:
        # make sure, guessed_id is at the beginning
        old_order = order.copy()
        if guessed_id in old_order:
            pos = old_order.index(guessed_id)
            old_order.pop(pos)
        order = [guessed_id]
        order.extend(old_order)

    order_rev = order.copy()
    order_rev.reverse()
    for k in order_rev:
        prop_list.sort(key=lambda x: k.lower() in x.lower(), reverse=True)

    header_opt = {"bgcolor": "#94CAEF", "style": "text-align: center;"}

    table = []
    rows = []
    cells = []
    # first line
    cells.extend(html.td(html.b("Molecule"), options=header_opt))
    header_opt["colspan"] = 2
    cells.extend(html.td(html.b("Properties"), options=header_opt))
    rows.extend(html.tr(cells))

    cells = []
    bgcolor = "#F2F2F2"
    for idx, mol in enumerate(mol_list, 1):
        if not mol:
            continue

        # alternating background colors for easier distinction between records
        if "F2" in bgcolor:
            bgcolor = "#FFFFFF"
        else:
            bgcolor = "#F2F2F2"

        # How many properties have to be displayed for the mol?
        mol_props = mol.GetPropNames()
        if props is None:
            props_to_show = mol_props
            row_span = len(mol_props)
        else:
            props_to_show = list(set(props).intersection(mol_props))
            row_span = len(props_to_show)

        td_opt = {"align": "center", "rowspan": row_span}

        if guessed_id:
            id_prop_val = mol.GetProp(guessed_id)
            img_id = id_prop_val
        else:
            img_id = idx

        if img_dir is None:  # embed the images in the doc
            b64 = b64_img(mol, size * 2)
            img_src = "data:image/png;base64,{}".format(b64)

        else:
            img_file = op.join(img_dir, "img_{}.png".format(img_id))
            img = autocrop(Draw.MolToImage(mol, size=(size * 2, size * 2)))
            img.save(img_file, format='PNG')
            img_src = img_file

        img_opt = {"title": str(img_id)}
        img_opt["style"] = 'max-width: {}px; max-height: {}px; display: block; margin: auto;'.format(size, size)

        cells.extend(html.td(html.img(img_src, img_opt), td_opt))

        # prop_opt = {}
        td_opt = {"bgcolor": bgcolor}
        for prop in prop_list:
            if prop not in props_to_show: continue

            cells.extend(html.td(prop, td_opt))
            cells.extend(html.td(mol.GetProp(prop), td_opt))
            rows.extend(html.tr(cells))
            cells = []

    table = html.table(rows)

    return "".join(table)


def show_table(sdf_list, id_prop=None, interact=False, highlight=None, order=None):
    return HTML(mol_table(sdf_list, id_prop, interact=interact, highlight=highlight, order=order))


def show_sheet(sdf_list, props=None, id_prop=None, interact=False, highlight=None, mols_per_row=4):
    return HTML(mol_sheet(sdf_list, props, id_prop, interact=interact, highlight=highlight, mols_per_row=mols_per_row))


def table_pager(mol_list, id_prop=None, interact=False, pagesize=25, highlight=None, order=None, show_hidden=False):
    l = len(mol_list)
    num_pages = l // pagesize
    if not WIDGETS or l <= pagesize:
        return HTML(mol_table(mol_list, id_prop=id_prop, highlight=highlight,
                              order=order, show_hidden=show_hidden))

    return ipyw.interactive(
        lambda page: HTML(mol_table(mol_list[page * pagesize:(page + 1) * pagesize],
                          id_prop=id_prop, interact=interact, order=order,
                          show_hidden=show_hidden)),
        page=ipyw.IntSlider(min=0, max=num_pages, step=1, value=0)
    )


def nested_pager(mol_list, pagesize=10, id_prop=None, props=None, order=None):
    l = len(mol_list)
    num_pages = l // pagesize
    if not WIDGETS or l <= pagesize:
        return HTML(nested_table(mol_list, id_prop=id_prop, props=props, order=order))

    return ipyw.interactive(
        lambda page: HTML(nested_table(mol_list[page * pagesize:(page + 1) * pagesize],
                          id_prop=id_prop, props=props, order=order)),
        page=ipyw.IntSlider(min=0, max=num_pages, step=1, value=0)
    )


def grid_pager(mol_list, pagesize=20, id_prop=None, interact=False, highlight=None, props=None, mols_per_row=4, size=IMG_GRID_SIZE):
    l = len(mol_list)
    num_pages = l // pagesize
    if not WIDGETS or l <= pagesize:
        return HTML(mol_sheet(mol_list, id_prop=id_prop, props=props, size=size))

    return ipyw.interactive(
        lambda page: HTML(mol_sheet(mol_list[page * pagesize:(page + 1) * pagesize],
                          id_prop=id_prop, interact=interact, highlight=highlight,
                          props=props, size=size)),
        page=ipyw.IntSlider(min=0, max=num_pages, step=1, value=0)
    )


def jsme(name="mol"):
    """displays a JSME molecule editor widget in the notebook
    and stores the resulting mol in the variable that <name> assigns."""

    time_stamp = time.strftime("%y%m%d%H%M%S")

    return HTML(JSME_FORM.format(jsme_loc=JSME_LOCATION, ts=time_stamp, var_name=name))


def dict_from_sdf_list(sdf_list, id_prop=None, props=None, prop_list=None):
    """Generate a dictionary from the properties of a list of molecules.
    Currently not including the structure.
    If <props> contains a list of property names, then only these properties plus the <id_prop> are returned.
    Returns dict"""

    if not prop_list:
        prop_list = list_fields(sdf_list)

    if id_prop:
        if id_prop not in prop_list:
            raise LookupError("id_prop not found in data set.")
        guessed_id = id_prop
    else:
        guessed_id = guess_id_prop(prop_list)

    if not props:
        props = prop_list
    if guessed_id and guessed_id not in props:
        props.append(guessed_id)

    df_dict = {prop: [] for prop in props}

    for mol in sdf_list:
        mol_props = list(mol.GetPropNames())
        for prop in props:
            if prop in mol_props:
                df_dict[prop].append(get_value(mol.GetProp(prop)))
            else:
                df_dict[prop].append(np.NaN)

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
    mol_list = deepcopy(input_list)
    writer = Chem.SDWriter(fn)

    print("N\t\tscore\t\trmsd")
    for ctr, mol in enumerate(mol_list, 1):
        mol_pymp = Chem.MMFFGetMoleculeProperties(mol)
        o3a = Chem.GetO3A(mol, ref, mol_pymp, ref_pymp)
        print("{}\t\t{:.2f}\t\t{:.2f}".format(ctr, o3a.Score(), o3a.Align()))
        writer.write(mol)

    writer.close()
