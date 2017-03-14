#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
###
SAR
###

*Created on Tue Mar 14, 2017 by A. Pahl*

SAR Tools.
"""

# import csv, os, pickle, random
import base64, sys, time
import os.path as op
# import gzip
# import math
# from copy import deepcopy

from sklearn.ensemble import RandomForestClassifier

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
# import rdkit.Chem.Descriptors as Desc

# imports for similarity search
# from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.Draw import SimilarityMaps
from rdkit import DataStructs

Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

# from PIL import Image, ImageChops

import numpy as np

from . import tools, html_templates as html, nb_tools as nbt


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

# try:
#     # Try to import Avalon so it can be used for generation of 2d coordinates.
#     from rdkit.Avalon import pyAvalonTools as pyAv
#     USE_AVALON = True
# except ImportError:
#     print("  * Avalon not available. Using RDKit for 2d coordinate generation.")
#     USE_AVALON = False

if AP_TOOLS:
    #: Library version
    VERSION = apt.get_commit(__file__)
    # I use this to keep track of the library versions I use in my project notebooks
    print("{:45s} (commit: {})".format(__name__, VERSION))
else:
    print("{:45s} ({})".format(__name__, time.strftime("%y%m%d-%H:%M", time.localtime(op.getmtime(__file__)))))


BGCOLOR = "#94CAEF"
IMG_GRID_SIZE = 235


def _get_proba(fp, predictionFunction):
    return predictionFunction(fp)[0][1]


def b64_fig(fig, dpi=72):
    img_file = IO()
    fig.savefig(img_file, dpi=dpi, format='PNG', bbox_inches="tight")
    b64 = base64.b64encode(img_file.getvalue())
    if PY3:
        b64 = b64.decode()
    img_file.close()
    return b64


class SAR_List(tools.Mol_List):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.model = None


    def train(self, act_class_prop="ActClass"):
        self.model = train(self, act_class_prop)
        """Generates the trained model."""


    def predict(self):
        """Adds predictions from the trained model to the SAR_List.
        Model has to be available as `self.model`."""
        if self.model is None:
            raise LookupError("Model is not available. Please train first.")
        predict(self, self.model)


    def sim_map(self):
        return HTML(sim_map(self, self.model, id_prop=self.id_prop, order=self.order))


    def write_sim_map(self, fn="sim_map.html", title="Similarity Map", summary=None):
        html.write(html.page(sim_map(self, self.model, id_prop=self.id_prop, order=self.order),
                             summary=summary, title=title), fn=fn)
        return HTML('<a href="{}">{}</a>'.format(fn, fn))



def train(mol_list, act_class_prop="ActClass"):
    """Returns the trained model."""
    fps = []
    act_classes = []
    for mol in mol_list:
        fps.append(Chem.GetMorganFingerprintAsBitVect(mol, 2))
        act_classes.append(tools.get_value(mol.GetProp(act_class_prop)))
    np_fps = []
    for fp in fps:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)

    # get a random forest classifiert with 100 trees
    rf = RandomForestClassifier(n_estimators=100, random_state=1123)
    rf.fit(np_fps, act_classes)
    return rf


def predict_mol(mol, model):
    """Returns the predicted class and the probabilities for a molecule.

    Parameters:
        model: Output from `train()`."""
    fp = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(Chem.GetMorganFingerprintAsBitVect(mol, 2), fp)
    fp = fp.reshape(1, -1)   # this removes the deprecation warning
    predict_class = model.predict(fp)
    predict_prob = model.predict_proba(fp)
    return predict_class[0], predict_prob[0]


def predict(mol_list, model):
    for mol in mol_list:
        pred_class, pred_prob = predict_mol(mol, model)
        mol.SetProp("ActClass_Pred", str(pred_class))
        mol.SetProp("ActClass_Prob", str(pred_prob[pred_class]))


def sim_map(mol_list, model, id_prop=None, interact=False, highlight=None, show_hidden=False, order=None, size=300):
    """Parameters:
        mol_list (Mol_List): List of RDKit molecules
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
    prop_list = tools.list_fields(mol_list)

    if isinstance(order, list):
        for k in reversed(order):
            prop_list.sort(key=lambda x: k.lower() in x.lower(), reverse=True)

    if id_prop is None:
        guessed_id = tools.guess_id_prop(prop_list)
    else:
        guessed_id = id_prop

    if interact and guessed_id is not None:
        table_list.append(tools.TBL_JAVASCRIPT.format(ts=time_stamp, bgcolor="transparent"))

    if id_prop is not None:
        if id_prop not in prop_list:
            raise LookupError("Id property {} not found in data set.".format(id_prop))

    pb = nbt.ProgressbarJS()

    if guessed_id:
        # make sure that the id_prop (or the guessed id prop) is first:
        prop_list.pop(prop_list.index(guessed_id))
        tmp_list = [guessed_id]
        tmp_list.extend(prop_list)
        prop_list = tmp_list

    cells = html.td(html.b("#"), header_opt)
    cells.extend(html.td(html.b("Molecule"), header_opt))
    cells.extend(html.td(html.b("SimMap"), header_opt))
    for prop in prop_list:
        cells.extend(html.td(html.b(prop), header_opt))
    rows = html.tr(cells)

    list_len = len(mol_list)
    for idx, mol in enumerate(mol_list):
        pb.update(100 * (idx + 1) / list_len)
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
            b64 = tools.b64_img(mol, size * 2)
            img_src = "data:image/png;base64,{}".format(b64)
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

            fig, _ = SimilarityMaps.GetSimilarityMapForModel(mol, SimilarityMaps.GetMorganFingerprint, lambda x: _get_proba(x, model.predict_proba))
            b64 = b64_fig(fig, dpi=72)
            img_src = "data:image/png;base64,{}".format(b64)
            cell_opt = {}
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
        table_list.append(tools.ID_LIST.format(ts=time_stamp))

    pb.done()
    return "".join(table_list)
