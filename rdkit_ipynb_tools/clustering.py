#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
##########
Clustering
##########

*Created on Sun Feb 28 11:00 2016 by A. Pahl*

Clustering molecules.
"""

import os
import os.path as op
from copy import deepcopy
from collections import Counter, OrderedDict
import shutil

from rdkit.Chem import AllChem as Chem, MACCSkeys
from rdkit.Chem import Draw
import rdkit.Chem.Descriptors as Desc
from rdkit.Chem.AtomPairs import Pairs, Torsions

# import rdkit.Chem.Scaffolds.MurckoScaffold as MurckoScaffold
try:
    Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
    Draw.DrawingOptions.atomLabelFontSize = 18
except KeyError:  # Font "DejaVu Sans" is not available
    pass

from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

# from PIL import Image, ImageChops
import numpy as np
try:
    # this mainly needs to be done because Sphinx can not import Matplotlib
    # (it gives `NotImplementedError('Implement enable_gui in a subclass')`)
    import matplotlib.pyplot as plt
    MPL = True
except ImportError:
    MPL = False

# from . import html_templates as html
from . import tools, html_templates as html, file_templ as ft, nb_tools as nbt

try:
    from rdkit.Avalon import pyAvalonTools as pyAv
    USE_AVALON = True
except ImportError:
    USE_AVALON = False


nbits = 1024
nbits_long = 16384

# dictionary
FPDICT = {}
FPDICT['ecfp0'] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 0, nBits=nbits)
FPDICT['ecfp2'] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 1, nBits=nbits)
FPDICT['ecfp4'] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 2, nBits=nbits)
FPDICT['ecfp6'] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 3, nBits=nbits)
FPDICT['ecfc0'] = lambda m: Chem.GetMorganFingerprint(m, 0)
FPDICT['ecfc2'] = lambda m: Chem.GetMorganFingerprint(m, 1)
FPDICT['ecfc4'] = lambda m: Chem.GetMorganFingerprint(m, 2)
FPDICT['ecfc6'] = lambda m: Chem.GetMorganFingerprint(m, 3)
FPDICT['fcfp2'] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 1, useFeatures=True, nBits=nbits)
FPDICT['fcfp4'] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 2, useFeatures=True, nBits=nbits)
FPDICT['fcfp6'] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 3, useFeatures=True, nBits=nbits)
FPDICT['fcfc2'] = lambda m: Chem.GetMorganFingerprint(m, 1, useFeatures=True)
FPDICT['fcfc4'] = lambda m: Chem.GetMorganFingerprint(m, 2, useFeatures=True)
FPDICT['fcfc6'] = lambda m: Chem.GetMorganFingerprint(m, 3, useFeatures=True)
FPDICT['lecfp4'] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 2, nBits=nbits_long)
FPDICT['lecfp6'] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 3, nBits=nbits_long)
FPDICT['lfcfp4'] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 2, useFeatures=True, nBits=nbits_long)
FPDICT['lfcfp6'] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 3, useFeatures=True, nBits=nbits_long)
FPDICT['maccs'] = lambda m: MACCSkeys.GenMACCSKeys(m)
FPDICT['ap'] = lambda m: Pairs.GetAtomPairFingerprint(m)
FPDICT['tt'] = lambda m: Torsions.GetTopologicalTorsionFingerprintAsIntVect(m)
FPDICT['hashap'] = lambda m: Desc.GetHashedAtomPairFingerprintAsBitVect(m, nBits=nbits)
FPDICT['hashtt'] = lambda m: Desc.GetHashedTopologicalTorsionFingerprintAsBitVect(m, nBits=nbits)
FPDICT['rdk5'] = lambda m: Chem.RDKFingerprint(m, maxPath=5, fpSize=nbits, nBitsPerHash=2)
FPDICT['rdk6'] = lambda m: Chem.RDKFingerprint(m, maxPath=6, fpSize=nbits, nBitsPerHash=2)
FPDICT['rdk7'] = lambda m: Chem.RDKFingerprint(m, maxPath=7, fpSize=nbits, nBitsPerHash=2)
if USE_AVALON:
    FPDICT['avalon'] = lambda m: pyAv.GetAvalonFP(m, nbits)
    FPDICT['avalon_l'] = lambda m: pyAv.GetAvalonFP(m, nbits_long)


def mpl_hist(data, bins=10, xlabel="values", ylabel="Occurrence", show=False, save=True, **kwargs):
    """Useful kwargs: size (tuple<int>), dpi (int), fn (filename, str), title (str)"""
    my_dpi = kwargs.get("dpi", 96)
    size = kwargs.get("size", (300, 350))
    title = kwargs.get("title", None)
    figsize = (size[0] / my_dpi, size[1] / my_dpi)
    plt.style.use('seaborn-pastel')
    # plt.style.use('ggplot')
    plt.style.use('seaborn-whitegrid')
    fig = plt.figure(figsize=figsize, dpi=my_dpi)
    if title is not None:
        fig.suptitle(title, fontsize=24)
    plt.hist(data, bins=bins)
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=16)

    if save:
        fn = kwargs.get("fn", "hist.png")
        plt.savefig(fn, bbox_inches='tight')

    if show:
        plt.show()


def renumber_clusters(cluster_list, start_at=1):
    """Renumber clusters in-place."""
    start_at -= 1
    # get the current individual cluster numbers present in the list
    id_list = sorted(set(tools.get_value(mol.GetProp("Cluster_No"))
                         for mol in cluster_list.mols_with_prop("Cluster_No")))

    # assign the new ids as values the old id's keys
    new_ids = {k: v for v, k in enumerate(id_list, 1 + start_at)}
    for mol in cluster_list:
        if not mol.HasProp("Cluster_No"): continue
        old_id = int(mol.GetProp("Cluster_No"))
        mol.SetProp("Cluster_No", str(new_ids[old_id]))


def get_cluster_numbers(cluster_list):
    """Returns the cluster numbers present in the cluster_list as a list, keeping the original order."""
    cl_no_od = OrderedDict()
    for mol in cluster_list.mols_with_prop("Cluster_No"):
        cl_no = int(mol.GetProp("Cluster_No"))
        cl_no_od[cl_no] = 0

    return list(cl_no_od.keys())


def get_clusters_by_no(cluster_list, cl_no, make_copy=True, renumber=False):
    """Return one or more clusters (provide numbers as list) by their number."""
    if not isinstance(cl_no, list):
        cl_no = [cl_no]

    cluster = tools.Mol_List()
    if cluster_list.order:
        cluster.order = cluster_list.order.copy()
    for mol in cluster_list:
        if mol.HasProp("Cluster_No") and int(mol.GetProp("Cluster_No")) in cl_no:
            if make_copy:
                mol = deepcopy(mol)
            cluster.append(mol)

    if renumber:
        renumber_clusters(cluster)

    return cluster


def remove_clusters_by_no(cluster_list, cl_no, make_copy=True, renumber=False):
    """Return a new cluster list where the clusters with the provided numbers are removed."""
    if not isinstance(cl_no, list):
        cl_no = [cl_no]

    cluster = tools.Mol_List()
    for mol in cluster_list:
        if mol.HasProp("Cluster_No") and int(mol.GetProp("Cluster_No")) not in cl_no:
            if make_copy:
                mol = deepcopy(mol)
            cluster.append(mol)

    if renumber:
        renumber_clusters(cluster)

    return cluster


def keep_clusters_by_len(cluster_list, min_len=3, max_len=1000, make_copy=True, renumber=False):
    """Returns a new cluster list with all cluster removed where `len < min_len`."""
    ctr = Counter()
    result_list = tools.Mol_List()
    id_prop = tools.guess_id_prop(tools.list_fields(cluster_list))

    # Get the lengths, even if there are no cores:
    for mol in cluster_list:
        if int(mol.GetProp(id_prop)) < 100000: continue  # is a core

        if not mol.HasProp("Cluster_No"): continue
        cl_id = int(mol.GetProp("Cluster_No"))
        ctr[cl_id] += 1

    # now, only keep mols which belong to clusters of the desired length,
    # including the cores
    for mol in cluster_list:
        if not mol.HasProp("Cluster_No"): continue
        cl_id = int(mol.GetProp("Cluster_No"))
        if ctr[cl_id] >= min_len and ctr[cl_id] <= max_len:
            result_list.append(mol)

    if renumber:
        renumber_clusters(result_list)

    return result_list


def get_cores(cluster_list, make_copy=True):
    """Find and return the core molecules in a cluster_list."""
    core_list = cluster_list.has_prop_filter("is_core", make_copy=make_copy)
    core_list.order = ["Compound_Id", "Cluster_No", "Num_Members", "Num_Values", "Min", "Max", "Mean", "Median", "is_core"]
    return core_list


def get_members(cluster_list, make_copy=True):
    """Find and return the members of a cluster_list, exclude the cores."""
    member_list = cluster_list.has_prop_filter("is_core", invert=True, make_copy=make_copy)
    return member_list


def get_stats_for_cluster(cluster_list, activity_prop=None):
    stats = {}
    sups_raw_list = [mol.GetProp("Supplier") for mol in cluster_list.mols_with_prop("Supplier")]
    if sups_raw_list:
        # each supplier field of a mol may already contain "; "
        sups_raw_str = "; ".join(sups_raw_list)
        sups_set = set(sups_raw_str.split("; "))
        stats["Supplier"] = "; ".join(sorted(sups_set))

    prod_raw_list = [mol.GetProp("Producer") for mol in cluster_list.mols_with_prop("Producer")]
    if prod_raw_list:
        # each Producer field of a mol may already contain "; "
        prod_raw_str = "; ".join(prod_raw_list)
        prod_set = set(prod_raw_str.split("; "))
        stats["Producer"] = "; ".join(sorted(prod_set))

    if activity_prop is not None:
        value_list = [tools.get_value(mol.GetProp(activity_prop)) for mol in cluster_list.mols_with_prop(activity_prop)]
        stats["Num_Values"] = len(value_list)
        stats["Min"] = min(value_list) if value_list else None
        stats["Max"] = max(value_list) if value_list else None
        stats["Mean"] = np.mean(value_list) if value_list else None
        stats["Median"] = np.median(value_list) if value_list else None

    return stats


def add_stats_to_cores(cluster_list_w_cores, props=None):
    """Add statistical information for the given props to the cluster cores.
    If 'props' is None, a default list of properties is used.
    The cores have to be already present in the list."""
    if props is None:
        props = ["ALogP", "QED", "SA_Score"]
    elif not isinstance(props, list):
        props = [props]

    cores = get_cores(cluster_list_w_cores, make_copy=False)
    if len(cores) == 0:
        raise LookupError("Could not find any cores in the list! Please add them first with add_cores() or add_centers()")
    cores.order = ["Cluster_No", "Num_Members"]
    for prop in props:
        cores.order.extend(["{}_Min".format(prop), "{}_Max".format(prop), "{}_Mean".format(prop),
                           "{}_Median".format(prop), "{}_#Values".format(prop)])
    for core in cores:
        cl_no = tools.get_prop_val(core, "Cluster_No")
        cluster = get_members(get_clusters_by_no(cluster_list_w_cores, cl_no, make_copy=False),
                              make_copy=False)
        for prop in props:
            value_list = list(filter(
                lambda x: x is not None, [tools.get_prop_val(mol, prop) for mol in cluster]))
            if len(value_list) == 0: continue

            core.SetProp("{}_#Values".format(prop), str(len(value_list)))

            if all(tools.isnumber(x) for x in value_list):
                core.SetProp("{}_Min".format(prop), "{:.2f}".format(min(value_list)))
                core.SetProp("{}_Max".format(prop), "{:.2f}".format(max(value_list)))
                core.SetProp("{}_Mean".format(prop), "{:.2f}".format(np.mean(value_list)))
                core.SetProp("{}_Median".format(prop), "{:.2f}".format(np.median(value_list)))

            else:  # summarize the values as strings
                value_str = "; ".join(str(x) for x in value_list)
                value_list_ext = value_str.split("; ")
                value_ctr = Counter(value_list_ext)
                prop_values = []
                for val in sorted(value_ctr):
                    prop_values.append("{} ({})".format(val, value_ctr[val]))
                core.SetProp("{}s".format(prop), "; ".join(prop_values))


def get_clusters_with_activity(cluster_list, activity_prop, min_act=None, max_act=None, min_len=1, renumber=False):
    """Return only the clusters which fulfill the given activity criteria.

    Parameters:
        min_act (str): something like `< 50` or `> 70` that can be evaluated.
        max_act (str): see above."""

    if min_act is not None:
        min_act_comp = compile('stats["Min"] {}'.format(min_act), '<string>', 'eval')
    if max_act is not None:
        max_act_comp = compile('stats["Max"] {}'.format(max_act), '<string>', 'eval')

    cores_and_members = get_cores(cluster_list)
    if len(cores_and_members) > 0:
        cores_and_members = cores_and_members.prop_filter('Num_Members >= {}'.format(min_len))
    members_all = get_members(cluster_list)
    tmp_list = tools.Mol_List()

    cl_ids = sorted(set(int(mol.GetProp("Cluster_No")) for mol in members_all.mols_with_prop("Cluster_No")))
    for new_id, cl_id in enumerate(cl_ids, 1):
        cluster = get_clusters_by_no(members_all, cl_id)
        if len(cluster) < min_len: continue
        stats = get_stats_for_cluster(cluster, activity_prop)
        if stats["Num_Values"] == 0: continue
        stats["Min"]  # to quiet the linter
        keep = True
        if min_act is not None and not eval(min_act_comp):
            keep = False
        if max_act is not None and not eval(max_act_comp):
            keep = False

        if keep:
            tmp_list.extend(cluster)

    cores_and_members.extend(tmp_list)

    if renumber:
        renumber_clusters(cores_and_members)

    return cores_and_members


def add_cores(cluster_list, activity_prop=None, align_to_core=False):
    """Find and add cores to the cluster_list in-place.

    Parameters:
        mode (str): "mcs" (default): add the core as MCS of the cluster. This can lead to very small structures.
        activity_prop (str): the name of the property for the activity."""

    # first, remove any already existing cores
    members_all = get_members(cluster_list)

    # find cluster numbers
    cl_ids = set(tools.get_value(mol.GetProp("Cluster_No")) for mol in members_all.mols_with_prop("Cluster_No"))

    for cl_id in cl_ids:
        cluster = get_clusters_by_no(members_all, cl_id, make_copy=False)
        if len(cluster) < 3:
            continue

        # determine the cluster core by MCSS
        # do not generate cores for singletons and pairs

        core_mol = tools.find_mcs(cluster)
        if core_mol is None:
            continue

        tools.check_2d_coords(core_mol)

        # set a number of properties for the cluster core
        id_prop = tools.guess_id_prop(cluster[0].GetPropNames())
        core_mol.SetProp(id_prop, str(cl_id))
        core_mol.SetProp("is_core", "yes")
        core_mol.SetProp("Cluster_No", str(cl_id))
        core_mol.SetProp("Num_Members", str(len(cluster)))

        if align_to_core:
            # align_mol = MurckoScaffold.GetScaffoldForMol(core_mol)
            cluster.align(core_mol)

        members_all.insert(0, core_mol)

    return members_all


def add_centers(cluster_list, mode="most_active", activity_prop=None, min_num_members=3, **kwargs):
    """Add cluster centers to the cores.
    Contrary to the cores, this is not an MCS, but one of the cluster members,
    with a certain property.

    Parameters:
        mode (str): `most_active` (default): if activity_prop is not None, the most active compound is taken.
            `smallest`: the compound with the least amount of heavy atoms is taken as center.
            `center`: the compound with the medium number of heavy atoms is taken.
            `from_tag`: takes the molecule the `tag` property as center.
            The `tag` parameter needs to be defined."""

    if "active" in mode:
        if activity_prop is not None:
            al = activity_prop.lower()
            reverse = False
            if "pic50" in al or "pec50" in al:
                reverse = True

        else:
            mode = "center"

    # first, remove any already existing cores
    members_all = get_members(cluster_list)

    # find cluster numbers
    cl_ids = set(tools.get_value(mol.GetProp("Cluster_No")) for mol in members_all.mols_with_prop("Cluster_No"))

    for cl_id in cl_ids:
        cluster = get_clusters_by_no(members_all, cl_id, make_copy=False)
        if len(cluster) < min_num_members:
            continue

        if "active" in mode:
            cluster.sort_list(activity_prop, reverse=reverse)
            core_mol = deepcopy(cluster[0])

        elif "smallest" in mode:  # smallest
            cluster.sort(key=Desc.HeavyAtomCount)
            core_mol = deepcopy(cluster[0])

        elif "center" in mode:  # medium number of heavy atoms, the middle of the list
            cluster.sort(key=Desc.HeavyAtomCount)
            core_mol = deepcopy(cluster[len(cluster) // 2])

        elif "tag" in mode:
            tag = kwargs.get("tag", None)
            if tag is None:
                raise KeyError("Parameter `tag` is required but could not be found.")
            tmp_list = cluster.has_prop_filter(tag)
            if len(tmp_list) == 0:
                print("No core for cluster {} (tag: {})".format(cl_id, tag))
                continue
            core_mol = tmp_list[0]

        for prop in core_mol.GetPropNames():
            core_mol.ClearProp(prop)

        # set a number of properties for the cluster core
        id_prop = tools.guess_id_prop(cluster[0].GetPropNames())
        core_mol.SetProp(id_prop, str(cl_id))
        core_mol.SetProp("is_core", "yes")
        core_mol.SetProp("Cluster_No", str(cl_id))
        core_mol.SetProp("Num_Members", str(len(cluster)))

        members_all.insert(0, core_mol)

    return members_all


def get_mol_list_from_index_list(orig_sdf, index_list, cl_id):
    """generate sdf_lists after clustering"""
    cluster_list = tools.Mol_List()
    for x in index_list:
        mol = deepcopy(orig_sdf[x])
        mol.SetProp("Cluster_No", str(cl_id))
        cluster_list.append(mol)

    if len(cluster_list) == 2:
        cluster_list[0].SetProp("is_pair", "yes")
        cluster_list[1].SetProp("is_pair", "yes")
    elif len(index_list) == 1:
        cluster_list[0].SetProp("is_single", "yes")

    return cluster_list


def cluster_from_mol_list(mol_list, cutoff=0.8, fp="ecfp6", activity_prop=None,
                          summary_only=True, generate_cores=False, align_to_core=False):
    """Clusters the input Mol_List.

    Parameters:
        mol_list (tools.Mol_List): the input molecule list.
        cutoff (float): similarity cutoff for putting molecules into the same cluster.

    Returns:
        A new Mol_List containing the input molecules with their respective cluster number,
        as well as additionally the cluster cores, containing some statistics."""

    try:
        fp_func = FPDICT[fp]
    except KeyError:
        print("Fingerprint {} not found. Available fingerprints are: {}".format(fp, ", ".join(sorted(FPDICT.keys()))))
        return

    counter = Counter()

    # generate the fingerprints
    fp_list = [fp_func(mol) for mol in mol_list]

    # second generate the distance matrix:
    dists = []
    num_of_fps = len(fp_list)
    for i in range(1, num_of_fps):
        sims = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        dists.extend([1 - x for x in sims])

    # now cluster the data:
    cluster_idx_list = Butina.ClusterData(dists, num_of_fps, cutoff, isDistData=True)
    for cluster in cluster_idx_list:
        counter[len(cluster)] += 1
    print("    fingerprint:", fp)
    print("    clustersize  num_of_clusters")
    print("    ===========  ===============")
    for length in sorted(counter.keys(), reverse=True):
        print("        {:4d}            {:3d}".format(length, counter[length]))
    print()

    if summary_only:
        return None

    cluster_list = tools.Mol_List()

    # go over each list of indices to collect the cluster's molecules
    for cl_id, idx_list in enumerate(sorted(cluster_idx_list, key=len, reverse=True), 1):
        cluster = get_mol_list_from_index_list(mol_list, idx_list, cl_id)
        cluster[0].SetProp("is_repr", "yes")  # The first compound in a cluster is the representative
        cluster_list.extend(cluster)

    if generate_cores:
        cluster_list = add_cores(cluster_list, activity_prop, align_to_core)

    return cluster_list


def show_numbers(cluster_list, show=True):
    """Calculate (and show) some numbers for the cluster_list. Returns a dict."""
    all_members = get_members(cluster_list)
    ctr_cl_no = Counter([int(mol.GetProp("Cluster_No")) for mol in all_members])
    sizes = [s for s in ctr_cl_no.values()]
    ctr_size = Counter(sizes)

    if show:
        total = 0
        print("\nCluster Size  |  Number of Clusters")
        print("------------- + -------------------")

        for i in sorted(ctr_size, reverse=True):

            print("     {:3d}      |   {:3d}".format(i, ctr_size[i]))
            total += (i * ctr_size[i])

        print()
        print("Number of compounds as sum of members per cluster size:", total)

    return sizes


def core_table(mol, props=None, hist=None):
    if props is None:
        props = ["Cluster_No", "Num_Members", "Producers"]

    td_opt = {"align": "center"}
    header_opt = {"bgcolor": "#94CAEF", "align": "center"}
    table_list = []

    cells = html.td(html.b("Molecule"), header_opt)
    for prop in props:
        pl = prop.lower()
        if (pl.endswith("min") or pl.endswith("max") or pl.endswith("mean") or
                pl.endswith("median")or pl.endswith("ic50") or pl.endswith("ic50)") or
                pl.endswith("activity") or pl.endswith("acctivity)")):
            pos = prop.rfind("_")
            if pos > 0:
                prop = prop[:pos] + "<br>" + prop[pos + 1:]
        cells.extend(html.td(html.b(prop), header_opt))

    if hist is not None:
        header_opt["class"] = "histogram"
        cells.extend(html.td(html.b("Histogram"), header_opt))

    rows = html.tr(cells)

    cells = []

    if not mol:
        cells.extend(html.td("no structure"))

    else:
        mol_props = mol.GetPropNames()
        cl_no = mol.GetProp("Cluster_No")
        img_file = "img/core_{}.png".format(cl_no)
        img = tools.autocrop(Draw.MolToImage(mol))
        img.save(img_file, format='PNG')
        img_src = img_file

        cell = html.img(img_src)
        cells.extend(html.td(cell, td_opt))

    for prop in props:
        td_opt = {"align": "center"}
        if prop in mol_props:
            prop_val = mol.GetProp(prop)
            cells.extend(html.td(prop_val, td_opt))
        else:
            cells.extend(html.td("", td_opt))

    if hist is not None:
        td_opt["class"] = "histogram"
        if "img/" not in hist:
            hist = "img/" + hist
        img_opt = {"height": "220"}
        cell = html.img(hist, img_opt)
        cells.extend(html.td(cell, td_opt))

    rows.extend(html.tr(cells))

    table_list.extend(html.table(rows))

    # print(table_list)
    return "".join(table_list)


def write_report(cluster_list, title="Clusters", props=None, reverse=True, **kwargs):
    """Useful kwargs: core_props (list, props to show for the core,
    default: ["Cluster_No", "Num_Members", "Producers"]). The exact names of the props have to be given (with `_Mean` etc.).
    bins (int or list, default=10), align (bool)
    add_stats (bool): whether to add the statistics on the fly
    or use any precalculated ones. Default: True
    show_hist (bool): whether to show histograms or not (default: True)."""
    resource_dir = op.join(op.dirname(__file__), "resources")
    cur_dir = op.abspath(op.curdir)

    core_props = kwargs.get("core_props", None)
    align = kwargs.get("align", False)
    content = [ft.CLUSTER_REPORT_INTRO]
    bins = kwargs.get("bins", 10)
    show_hist = kwargs.get("show_hist", True)
    add_stats = kwargs.get("add_stats", False)

    pb = nbt.ProgressbarJS()

    print("  Copying resources...")
    if op.isdir("./clustering"):
        print("* Clustering dir already exists, writing into...")
    else:
        shutil.copytree(op.join(resource_dir, "clustering"), "./clustering")

    os.chdir("clustering")
    if add_stats:
        print("  Adding statisical information...")
        props_to_stat = []
        for prop in core_props:
            if "Num_Members" in prop or "Cluster_No" in prop: continue  # no stats for these props!
            for stat_type in ["_Min", "_Max", "_Mean", "_Median", "s"]:
                if prop.endswith(stat_type):
                    prop = prop[:-(len(stat_type))]
            if props not in props_to_stat:
                props_to_stat.append(prop)
        add_stats_to_cores(cluster_list, props_to_stat)

    print("  Generating Report...")

    # collect the cluster numbers in the order in which they are in the cluster_list:
    cluster_numbers = OrderedDict()
    for mol in cluster_list.has_prop_filter("Cluster_No"):
        cluster_numbers[int(mol.GetProp("Cluster_No"))] = 0  # dummy value

    if props is not None and not isinstance(props, list):
        props = [props]

    len_cluster_numbers = len(cluster_numbers)
    for idx, cl_no in enumerate(cluster_numbers, 1):
        pb.update(100 * idx / len_cluster_numbers)
        cluster = get_clusters_by_no(cluster_list, cl_no)
        if len(cluster) == 0: continue
        if align and len(cluster) > 1:
            cluster.align()

        hist_fn = None
        if props is not None:
            first_prop = props[0]
            cluster.sort_list(first_prop, reverse=reverse)
            if MPL and show_hist and len(cluster) > 4:
                hist_fn = "img/hist_{}.png".format(cl_no)
                data = [tools.get_value(mol.GetProp(first_prop)) for mol in cluster if mol.HasProp(first_prop)]
                mpl_hist(data, bins=bins, xlabel=first_prop, fn=hist_fn)

        content.append("""<br>\n<h2 id="{}">Cluster {:03d}</h2>""".format(cl_no, cl_no))
        core = get_cores(cluster)
        if len(core) > 0:
            content.append("<ul>\n<li>\n")
            content.append(core_table(core[0], props=core_props, hist=hist_fn))
            content.append("<ul>\n<li>\n<h4>Members:</h4>")

        members = get_members(cluster)
        content.append(members.grid(props=props, size=300, raw=True, img_dir="img"))
        content.append("<hr>")
        if len(core) > 0:
            content.append("</li>\n</ul></li>\n</ul>")

    content.append(ft.CLUSTER_REPORT_EXTRO)
    print("  Writing report...")
    open("index.html", "w").write("\n".join(content))
    os.chdir(cur_dir)
    print("  done. The report has been written to \n        {}.".format(op.join(cur_dir, "clustering", "index.html")))
    pb.done()


def add_remarks_to_report(remarks, highlights=None, index_file="clustering/index.html"):
    """Takes remarks and adds them to the clustering report by cluster number.
    Writes a new `index_remarks.html` file

    Parameters:
        remarks (dict): keys are cluster numbers, values are lists of two elements,
            first is true/false for highlighting of the cluster, second is the remark as text
        highlights (list): list of strings whose occurrence in the report should be highlighted
        index_file (str): the name of the report index file"""

    if isinstance(highlights, str):  # make a list out of a single string argument
        highlights = [highlights]
    new_file = open(op.join(op.dirname(index_file), "index_remarks.html"), "w")
    for line in open(index_file):
        for cl_no in remarks:
            if "Cluster {:03d}".format(cl_no) in line:
                hl = remarks[cl_no][0]
                txt = remarks[cl_no][1]
                if hl:
                    line = '<h2 style="background-color: #ffff4d" id="{0}">Cluster {0:03d}</h2>'.format(cl_no)
                else:
                    line = '<h2 id="{0}">Cluster {0:03d}</h2>'.format(cl_no)
                if len(txt) > 1:
                    line += '<p>{0}</p>\n'.format(txt)
                else:
                    line += '\n'
                break
        if highlights is not None:
            for hl in highlights:
                if hl in line:
                    line = line.replace(hl, '<div style="background-color: #ffff4d">{}</div>'.format(hl))
        new_file.write(line)
