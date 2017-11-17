#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
###
SAR
###

*Created on Tue Mar 14, 2017 by A. Pahl*

SAR Tools.
"""

# import csv, os
import base64, pickle, sys, time
import os.path as op
from collections import Counter
import re
import colorsys

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw

from rdkit.Chem.Draw import SimilarityMaps
from rdkit import DataStructs

Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

import numpy as np
from sklearn.ensemble import RandomForestClassifier

from . import tools, html_templates as html, nb_tools as nbt


from IPython.core.display import HTML, display, clear_output

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
    #: Library version
    VERSION = apt.get_commit(__file__)
    # I use this to keep track of the library versions I use in my project notebooks
    print("{:45s} (commit: {})".format(__name__, VERSION))
else:
    print("{:45s} ({})".format(__name__, time.strftime("%y%m%d-%H:%M", time.localtime(op.getmtime(__file__)))))


BGCOLOR = "#94CAEF"
IMG_GRID_SIZE = 235


TABLE_INTRO = """<table id="sar_table" width="" cellspacing="1" cellpadding="1" border="1" align="center" height="60" summary="">"""
HTML_INTRO = """<!DOCTYPE html>
<html>
<head>
  <title>%s</title>
  <meta charset="UTF-8">

  <link rel="stylesheet" type="text/css" href="css/style.css" />

  <script src="lib/float.js"></script>

</head>
<body>
<script src="lib/wz_tooltip.js"></script>
<h2>%s (%s)</h2>

"""
HTML_EXTRO = """<div style="width:4000px;height:2000px"></div>
<script>
      function addEvent(obj, ev, fu) {
      if (obj.addEventListener) {
          obj.addEventListener(ev, fu, false);
      } else {
          var eev = 'on' + ev;
          obj.attachEvent(eev, fu);
      }
      }
      addEvent(window, 'load', function () {
      tt1 = floatHeader('sar_table', {ncpth: [1], nccol: 1, topDif: 0, leftDif: 0});
      });
</script>
</body>
</html>"""

LOGP_INTRO = """</tbody>
</table>
<p></p>
<p>LogP color coding:</p>
<table width="" cellspacing="1" cellpadding="1" border="1" align="left" height="40" summary="">
<tbody>
<tr>
"""
LOGP_EXTRO = "</tbody>\n</table>\n"


class ColorScale():

    def __init__(self, num_values, val_min, val_max, middle_color="yellow",
                 reverse=False, is_lin=True):
        self.num_values = num_values
        self.num_val_1 = num_values - 1
        self.is_lin = is_lin  # is the range fo the color prop linear (e.g. pIC50)?
        if self.is_lin:
            self.value_min = val_min
            self.value_max = val_max
        else:
            self.value_min = tools.pic50(val_min, "um")
            self.value_max = tools.pic50(val_max, "um")
        self.reverse = reverse
        self.value_range = self.value_max - self.value_min
        self.color_scale = []
        if middle_color.startswith("y"):  # middle color yellow
            hsv_tuples = [(0.0 + ((x * 0.35) / (self.num_val_1)), 0.99, 0.9) for x in range(self.num_values)]
            self.reverse = not self.reverse
        else:  # middle color blue
            hsv_tuples = [(0.35 + ((x * 0.65) / (self.num_val_1)), 0.9, 0.9) for x in range(self.num_values)]
        rgb_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), hsv_tuples)
        for rgb in rgb_tuples:
            rgb_int = [int(255 * x) for x in rgb]
            self.color_scale.append('#{:02x}{:02x}{:02x}'.format(*rgb_int))

        if self.reverse:
            self.color_scale.reverse()

    def __call__(self, value):
        """return the color from the scale corresponding to the place in the value_min .. value_max range"""
        if not self.is_lin:
            value = tools.pic50(value, "um")
        pos = int(((value - self.value_min) / self.value_range) * self.num_val_1)

        return self.color_scale[pos]


    def legend(self):
        """Return the value_range and a list of tuples (value, color) to be used in a legend."""
        legend = []
        for idx, color in enumerate(self.color_scale):
            val = self.value_min + idx / self.num_val_1 * self.value_range
            if not self.is_lin:
                val = tools.ic50(val, "um")
            legend.append((val, color))

        return legend


def format_num(val):
    """Return a suitable format string depending on the size of the value."""
    if val > 50:
        return ".0f"
    if val > 1:
        return ".1f"
    else:
        return ".2f"


def _get_proba(fp, predictionFunction):
    return predictionFunction([fp])[0][1]


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
        self.html = None


    def _pass_properties(self, new_list):
            new_list.order = self.order
            new_list.ia = self.ia
            new_list.plot_tool = self.plot_tool
            new_list.model = self.model
            new_list.html = None


    def __getitem__(self, item):
        result = list.__getitem__(self, item)
        try:
            new_list = type(self)(result)

            # pass on properties
            self._pass_properties(new_list)
            return new_list
        except TypeError:
            return result


    def new(self, *args):
        new_list = type(self)(*args)
        # pass on properties
        self._pass_properties(new_list)
        return new_list


    def train(self, act_class_prop="AC_Real"):
        self.model = train(self, act_class_prop)
        """Generates the trained model."""


    def predict(self):
        """Adds predictions from the trained model to the SAR_List.
        Model has to be available as `self.model`."""
        if self.model is None:
            raise LookupError("Model is not available. Please train first.")
        predict(self, self.model)
        self.html = None


    def analyze(self, act_class="AC_Real", pred_class="AC_Pred"):
        """Prints the ratio of succcessful predictions for the molecules which have `act_class` and `pred_class` properties."""
        mol_ctr = Counter()
        hit_ctr = Counter()
        for mol in self:
            if mol.HasProp(act_class) and mol.HasProp(pred_class):
                mol_ctr[int(mol.GetProp(act_class))] += 1
                if mol.GetProp(act_class) != mol.GetProp(pred_class):
                    continue
                hit_ctr[int(mol.GetProp(act_class))] += 1
        if len(mol_ctr) > 0:
            sum_mol_ctr = sum(mol_ctr.values())
            sum_hit_ctr = sum(hit_ctr.values())
            print("Number of correctly predicted molecules: {} / {}    ({:.2f}%)"
                  .format(sum_hit_ctr, sum_mol_ctr, 100 * sum_hit_ctr /
                          sum_mol_ctr))
            print("\nCorrectly predicted molecules per Activity Class:")
            for c in sorted(hit_ctr):
                print("  {}:  {:.2f}".format(c, 100 * hit_ctr[c] / mol_ctr[c]))
        else:
            print("No molecules found with both {} and {}.".format(act_class, pred_class))
        return hit_ctr, mol_ctr


    def save_model(self, fn="sar"):
        if self.model is None:
            print("No model available.")
            return
        save_model(self.model, fn)


    def load_model(self, fn="sar", force=False):
        if self.model is not None and not force:
            print("There is already a model available. Use `force=True` to override.")
            return
        if not fn.endswith(".model"):
            fn = fn + ".model"
        with open(fn, "rb") as f:
            self.model = pickle.load(f)
        print("  > model loaded (last modified: {}).".format(time.strftime("%Y-%m-%d %H:%M", time.localtime(op.getmtime(fn)))))


    def sim_map(self):
        if self.html is None:
            self.html = sim_map(self, self.model, id_prop=self.id_prop, order=self.order)
        else:
            print("Using cached HTML content...")
            print("Set property `html` to `None` to re-generate.")
        return HTML(self.html)


    def write_sim_map(self, fn="sim_map.html", title="Similarity Map", summary=None):
        if self.html is None:
            self.html = sim_map(self, self.model, id_prop=self.id_prop, order=self.order)
        else:
            print("Using cached HTML content...")
            print("Set property `html` to `None` to re-generate.")
        html.write(html.page(self.html, summary=summary, title=title), fn=fn)
        return HTML('<a href="{}">{}</a>'.format(fn, fn))


    def map_from_id(self, cpd_id=None):
        if cpd_id is None and tools.WIDGETS:
            def show_sim_map(ev):
                cpd_id = tools.get_value(w_input_id.value.strip())
                clear_output()
                w_input_id.value = ""
                display(self.new_list_from_ids(cpd_id).sim_map())

            w_input_id = tools.ipyw.Text(description="Compound Id:")
            # w_btn_clear_input = tools.ipyw.Button(description="Clear Input")
            # w_btn_clear_input.on_click(clear_input)
            w_btn_show = tools.ipyw.Button(description="Show Sim Map")
            w_btn_show.on_click(show_sim_map)

            w_hb_show = tools.ipyw.HBox(children=[w_input_id, w_btn_show])
            # tools.set_margin(w_vb_search1)
            display(w_hb_show)
        else:
            new_list = self.new_list_from_ids(cpd_id)
            return HTML(sim_map(new_list, self.model, id_prop=self.id_prop, order=self.order))


def train(mol_list, act_class_prop="AC_Real"):
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
        mol.SetProp("AC_Pred", str(pred_class))
        mol.SetProp("Prob", "{:.2}".format(pred_prob[pred_class]))


def save_model(model, fn="sar"):
    if not fn.endswith(".model"):
        fn = fn + ".model"
    with open(fn, "wb") as f:
        pickle.dump(model, f)


def load_sdf(fn, model_name=None):
    mol_list = tools.load_sdf(fn)
    sar_list = SAR_List(mol_list)
    if model_name is None:
        print("  * No model was loaded. Please provide a name to load.")
    else:
        try:
            sar_list.load_model(model_name)
        except FileNotFoundError:
            print("  * Model {} could not be found. No model was loaded".format(model_name))
    return sar_list


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

    if len(mol_list) > 5:
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
        if len(mol_list) > 5:
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


            fig, _ = SimilarityMaps.GetSimilarityMapForModel(
                mol, SimilarityMaps.GetMorganFingerprint, lambda x: _get_proba(x, model.predict_proba))
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

    if len(mol_list) > 5:
        pb.done()
    return "".join(table_list)


def legend_table(legend):
    """Return a HTML table with the ColorScale label as text.

    Psrsmeters:
    legend (list): list of tuples as returned from ColorScale.legend()."""
    intro = "<table>\n<tbody>\n<tr>"
    extro = "</tr>\n</tbody>\n</table>\n"
    tbl_list = [intro]
    rnge = abs(legend[0][0] - legend[-1][0])
    digits = format_num(rnge)
    for tup in legend:
        cell = "<td bgcolor={color}>{val:{digits}}</td>".format(color=tup[1], val=tup[0], digits=digits)
        tbl_list.append(cell)

    tbl_list.append(extro)

    return "".join(tbl_list)


def get_res_pos(smiles):
    pat = re.compile('\[(.*?)\*\]')
    pos_str = re.findall(pat, smiles)[0]
    if pos_str:
        return int(pos_str)
    else:
        return 0


def generate_sar_table(db_list, core, id_prop, act_prop, sort_reverse=True,
                       dir_name="html/sar_table", color_prop="logp"):
    """core: smiles string; id_prop, act_prop: string
    colorprop_is_lin: whether or not the property used for coloring is linear (e.g. LogP or PercActivity) or needs to be logarithmitized (e.g. IC50_uM)."""

    tools.create_dir_if_not_exist(dir_name)
    tools.create_dir_if_not_exist(op.join(dir_name, "img"))

    db_list.sort_list(act_prop, reverse=sort_reverse)

    act_xy = np.zeros([55, 55], dtype=np.float)    # coordinates for the activity
    # color_xy = np.zeros([55, 55], dtype=np.float)
    color_xy = np.full([55, 55], np.NaN, dtype=np.float)
    molid_xy = np.zeros([55, 55], dtype=np.int)
    # molid_xy = np.arange(900, dtype=np.int).reshape(30, 30)  # coordinates for the molid
    rx_dict = {}  # axes for the residues
    ry_dict = {}
    max_x = -1  # keep track of the arraysize
    max_y = -1
    res_pos_x = -1
    res_pos_y = -1

    core_mol = Chem.MolFromSmiles(core)
    Draw.MolToFile(core_mol, "%s/img/core.png" % dir_name, [90, 90])

    for idx, mol in enumerate(db_list):
        act = float(mol.GetProp(act_prop))
        color = float(mol.GetProp(color_prop))
        molid = int(mol.GetProp(id_prop))
        tmp = Chem.ReplaceCore(mol, core_mol, labelByIndex=True)
        frag_mols = list(Chem.GetMolFrags(tmp, asMols=True))
        frag_smiles = [Chem.MolToSmiles(m, True) for m in frag_mols]
        if len(frag_mols) == 1:
            # one of the two residues is H:
            pos = get_res_pos(frag_smiles[0])
            if pos == res_pos_x:
                h_smiles = "[%d*]([H])" % res_pos_y
                frag_smiles.append(h_smiles)
                frag_mols.append(Chem.MolFromSmiles(h_smiles))
            else:
                h_smiles = "[%d*]([H])" % res_pos_x
                frag_smiles.insert(0, h_smiles)
                frag_mols.insert(0, Chem.MolFromSmiles(h_smiles))

            print(" adding H residue in pos {} to  mol #{} (molid: {})".format(pos, idx, mol.GetProp(id_prop)))

        elif len(frag_mols) > 2:
            print("*  incorrect number of fragments ({}) in mol #{} (molid: {})".format(len(frag_mols), idx, mol.GetProp(id_prop)))
            continue

        if res_pos_x == -1:
            # print frag_smiles[0], frag_smiles[1]
            res_pos_x = get_res_pos(frag_smiles[0])
            res_pos_y = get_res_pos(frag_smiles[1])
            # print "res_pos_x: {}     res_pos_y: {}".format(res_pos_x, res_pos_y)
        else:
            test_pos_x = get_res_pos(frag_smiles[0])
            if test_pos_x != res_pos_x:  # switch residues
                frag_smiles = frag_smiles[::-1]
                frag_mols = frag_mols[::-1]
        if frag_smiles[0] in rx_dict:
            curr_x = rx_dict[frag_smiles[0]]
        else:
            max_x += 1
            rx_dict[frag_smiles[0]] = max_x
            curr_x = max_x
            Draw.MolToFile(frag_mols[0], "%s/img/frag_x_%02d.png" % (dir_name, max_x), [100, 100])
        if frag_smiles[1] in ry_dict:
            curr_y = ry_dict[frag_smiles[1]]
        else:
            max_y += 1
            ry_dict[frag_smiles[1]] = max_y
            curr_y = max_y
            Draw.MolToFile(frag_mols[1], "%s/img/frag_y_%02d.png" % (dir_name, max_y), [100, 100])

        # draw thw whole molecule for the tooltip
        img_file = op.join(dir_name, "img/", "cpd_{}_{}.png".format(curr_x, curr_y))
        img = tools.autocrop(Draw.MolToImage(mol), "white")
        img.save(img_file, format='PNG')

        act_xy[curr_x][curr_y] = act
        color_xy[curr_x][curr_y] = color
        molid_xy[curr_x][curr_y] = molid

    return act_xy, molid_xy, color_xy, max_x, max_y


def sar_table_report_html(act_xy, molid_xy, color_xy, max_x, max_y, color_by="logp",
                          reverse_color=False, colorprop_is_lin=True,
                          show_link=False, show_tooltip=True):
    if "logp" in color_by.lower():
        # logp_colors = {2.7: "#5F84FF", 3.0: "#A4D8FF", 4.2: "#66FF66", 5.0: "#FFFF66", 1000.0: "#FF4E4E"}
        logp_colors = {2.7: "#98C0FF", 3.0: "#BDF1FF", 4.2: "#AAFF9B", 5.0: "#F3FFBF", 1000.0: "#FF9E9E"}

    else:
        color_min = float(np.nanmin(color_xy))
        color_max = float(np.nanmax(color_xy))
        color_scale = ColorScale(20, color_min, color_max,
                                 reverse=reverse_color, is_lin=colorprop_is_lin)

    # write horizontal residues
    line = [TABLE_INTRO]
    line.append("\n<thead><tr><th align=\"left\">Core:<br><img src=\"img/core.png\" alt=\"icon\" /></th>")
    for curr_x in range(max_x + 1):
        line.append("<th><img src=\"img/frag_x_%02d.png\" alt=\"icon\" /></th>" % curr_x)

    line.append("</tr></thead>\n<tbody>\n")

    for curr_y in range(max_y + 1):
        line.append("<tr><td><img src=\"img/frag_y_%02d.png\" alt=\"icon\" /></td>" % curr_y)
        for curr_x in range(max_x + 1):
            molid = molid_xy[curr_x][curr_y]
            if molid > 0:
                link_in = ""
                link_out = ""
                bg_color = " "
                mouseover = ""
                if show_link:
                    link = "../reports/ind_stock_results.htm#cpd_%05d" % molid
                    link_in = "<a href=\"%s\">" % link
                    link_out = "</a>"
                if "logp" in color_by.lower():
                    logp = color_xy[curr_x][curr_y]
                    if show_tooltip:
                        prop_tip = 'LogP: %.2f' % logp
                    for limit in sorted(logp_colors):
                        if logp <= limit:
                            bg_color = ' bgcolor="%s"' % logp_colors[limit]
                            break
                else:
                    value = float(color_xy[curr_x][curr_y])
                    html_color = color_scale(value)
                    bg_color = ' bgcolor="{}"'.format(html_color)
                    if show_tooltip:
                        prop_tip = '{}: {:.2f}'.format(color_by, color_xy[curr_x][curr_y])

                if show_tooltip:
                    tool_tip = '<img src=&quot;img/cpd_{}_{}.png&quot; alt=&quot;icon&quot; /><br><br>{}'.format(curr_x, curr_y, prop_tip)
                    mouseover = """ onmouseover="Tip('{}')" onmouseout="UnTip()" """.format(tool_tip)

                line.append("<td%s align=\"center\"%s><b>%.2f</b><br><br>(%s%d%s)</td>" % (mouseover, bg_color, act_xy[curr_x][curr_y], link_in, molid_xy[curr_x][curr_y], link_out))

            else:  # empty value in numpy array
                line.append("<td></td>")


        line.append("</tr>\n")

    line.append("</tbody>\n</table>\n")

    if "logp" in color_by.lower():
        line.append(LOGP_INTRO)
        for limit in sorted(logp_colors):
            line.append('<td align="center" bgcolor="%s">&le; %.2f</td>' % (logp_colors[limit], limit))
        line.append("\n</tr>\n")
        line.append(LOGP_EXTRO)
    else:
        line.append("<br><br>Coloring legend for {}:<br>\n".format(color_by))
        legend = color_scale.legend()
        line.append(legend_table(legend))


    html_table = "".join(line)

    return html_table


def write_html_page(html_content, dir_name="html/sar_table", page_name="sar_table", page_title="SAR Table"):

    tools.create_dir_if_not_exist(dir_name)
    tools.create_dir_if_not_exist(op.join(dir_name, "img"))

    filename = op.join(dir_name, "%s.htm" % page_name)
    f = open(filename, "w")
    f.write(HTML_INTRO % (page_title, page_title, time.strftime("%d-%b-%Y")))

    f.write(html_content)

    f.write(HTML_EXTRO)
    f.close()
