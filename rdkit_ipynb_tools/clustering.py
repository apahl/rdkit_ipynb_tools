#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
##########
Clustering
##########

*Created on Sun Feb 28 11:00 20165 by A. Pahl*

Clustering molecules.
"""


# import time
import sys
# import base64
import os.path as op
# import random
# import csv
# import gzip
# from copy import deepcopy
from collections import Counter

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw, rdFMCS
import rdkit.Chem.Descriptors as Desc
Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

# from PIL import Image, ImageChops

import numpy as np

# from . import html_templates as html
from . import tools

# from ipywidgets import widgets
# from IPython.core.display import HTML, display


def get_mol_list_from_index_list(orig_sdf, index_list, cl_id, activity_prop=None):
    """generate sdf_lists after clustering"""
    cluster_sdf = tools.Mol_List()
    for x in index_list:
        mol = deepcopy(orig_sdf[x])
        mol.SetProp("cluster_no", str(cl_id))
        cluster_sdf.append(mol)

    # determine the cluster core by MCSS
    mcs = rdFMCS.FindMCS(cluster_sdf)
    if mcs.canceled:
        print("* MCSS function timed out. Please provide a mol_or_smiles to align to.")
        return
    if mcs.smartsString:
        core_mol = Chem.MolFromSmarts(mcs.smartsString)
    else:
        print("* Could not find MCSS. Please provide a mol_or_smiles to align to.")
        return

    # set a number of properties for the cluster core
    core_mol.SetProp("Compound_Id", str(cl_id))
    core_mol.SetProp("cluster_no", str(cl_id))
    core_mol.SetProp("num_members", str(len(cluster_sdf)))

    if activity_prop is not None:
        value_list = [tools.get_value(mol.GetProp(activity_prop)) for mol in cluster_sdf.mols_with_prop(activity_prop)]
        core_mol.SetProp("num_values", str(len(value_list)))

            sum_d[prop]["min"] = min(value_list)
            sum_d[prop]["max"] = max(value_list)
            if sum_d[prop]["max"] > max_max:
                max_max = sum_d[prop]["max"]
            sum_d[prop]["mean"] = np.mean(value_list)
            sum_d[prop]["median"] = np.median(value_list)


    return cluster_sdf


def cluster_from_mol_list(mol_list, cutoff=0.8, activity_prop=None):
    """Clusters the input Mol_List.

    Parameters:
        mol_list (tools.Mol_List): the input molecule list.
        cutoff (float): similarity cutoff for putting molecules into the same cluster.

    Returns:
        A new Mol_List containing the input molecules with their respective cluster number,
        as well as additionally the cluster cores, containing some statistics."""

    counter = Counter()

    # generate the fingerprints
    fp_list = [Chem.GetMorganFingerprintAsBitVect(mol, 3, 1024) for mol in mol_list]

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
    print("    clustersize  num_of_clusters")
    print("    ===========  ===============")
    for length in sorted(counter.keys(), reverse=True):
        print("        {:4d}            {:3d}".format(length, counter[length]))

    cluster_idx_list.sort(key=len, reverse=True)
    clusters_and_cores = tools.Mol_List()

    # go over each list of indices to collect the cluster's molecules
    for cl_id, idx_list in enumerate(cluster_idx_list):
        cluster = get_mol_list_from_index_list(mol_list, idx_list, cl_id, activity_prop)



    del mol_list
    return clusters_and_cores


def get_max_act_in_cluster(orig_sdf, cluster, act_prop):
    max_act = -1000
    for x in cluster:
        mol = orig_sdf[x]
        try:
            value = float(mol.GetProp(act_prop))
        except:
            print("  * molecule at index {:6d} has not activity prop {}".format(x, act_prop))
            continue
        if value > max_act:
            max_act = value
    return max_act


def get_med_act_in_cluster(orig_sdf, cluster, act_prop):
    med_act = 0.0
    num_of_members = 0
    for x in cluster:
        mol = orig_sdf[x]
        try:
            value = float(mol.GetProp(act_prop))
        except:
            print("  * molecule at index {:6d} has no activity prop {}".format(x, act_prop))
            continue
        med_act += value
        num_of_members += 1

    if num_of_members == 0:
        num_of_members = 1

    return med_act / num_of_members


def analyze_cluster_idx_list(orig_sdf, cluster_idx_list, id_prop=None, mode="remove_singletons", act_prop=None):
    """available modes:
       remove_singletons: remove clusters with only one member
       ind_activity:      sort the clusters by the member with the highest activity
       med_activity:      sort the cluster by their medium activity

       a NEW cluster_idx_list is returned, use get_sdf_from_index_list to generate a sdf from it"""

    if act_prop:
        print("  > sorting cluster members on {}...".format(act_prop), end="")
        cluster_idx_list_sorted = []
        for cluster in cluster_idx_list:
            cluster_dict = {}

            # map the position in the orig_sdf to the molid
            sdf_list = get_sdf_from_index_list(orig_sdf, cluster)
            for pos, idx in enumerate(cluster):
                mol = sdf_list[pos]
                if id_prop is None and pos == 0:
                    id_prop = guess_id_prop(mol.GetPropNames())

                cluster_dict[int(mol.GetProp(id_prop))] = idx
            sort_sdf(sdf_list, act_prop, reverse=True)
            cluster_sorted = [cluster_dict[int(mol.GetProp(id_prop))] for mol in sdf_list]
            cluster_idx_list_sorted.append(cluster_sorted)
        cluster_idx_list = cluster_idx_list_sorted
        print(" done.")

    if mode == "remove_singletons":
        new_idx_list = [x for x in cluster_idx_list if len(x) > 1]
        print("  > removed {} singletons from cluster list".format(len(cluster_idx_list) - len(new_idx_list)))
        print("    the resulting list has {} clusters".format(len(new_idx_list)))
        return new_idx_list

    if mode == "ind_activity":
        new_idx_list = sorted(cluster_idx_list, key=lambda x: get_max_act_in_cluster(orig_sdf, x, act_prop), reverse=True)
        print("  > max ind_activity and number of members in the first ten clusters:")
        print("    index  ind_act  #members")
        print("    =====  =======  ========")
        for cl_counter, cluster in enumerate(new_idx_list):
            print("      {:2d}   {:6.1f}      {:3d}".format(cl_counter, get_max_act_in_cluster(orig_sdf, cluster, act_prop), len(cluster)))
            if cl_counter >= 9:
                break
        return new_idx_list

    if mode == "med_activity":
        new_idx_list = sorted(cluster_idx_list, key=lambda x: get_med_act_in_cluster(orig_sdf, x, act_prop), reverse=True)
        print("  > max med_activity in the first ten clusters:")
        for cl_counter, cluster in enumerate(new_idx_list):
            print("{}: {}   ".format(cl_counter, get_med_act_in_cluster(orig_sdf, cluster, act_prop)), end="")
            if cl_counter >= 9:
                break
        print()
        return new_idx_list

    print("  * unsupported mode.")


def write_clusters_as_sdf(orig_sdf, cluster_idx_list, basename="cluster"):
    """write every cluster in the index list as individual sdf file"""
    basename = basename.split(".")[0]

    for counter, cluster in enumerate(cluster_idx_list):
        num_of_clusters = len(cluster_idx_list)
        sdf_list = get_sdf_from_index_list(orig_sdf, cluster)
        write_sdf(sdf_list, "{}_{:03d}_{:03d}.sdf".format(basename, counter, num_of_clusters-1))


def write_cluster_report(orig_sdf, cluster_idx_list, captionprop, reportname="cluster"):
    intro = "<html>\n<head>\n  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />\n  <link rel=\"stylesheet\" type=\"text/css\" href=\"css/style.css\" />\n  <title>ClusterViewer</title>\n</head>\n<body>\n<table width=\"\" cellspacing=\"1\" cellpadding=\"1\" border=\"1\" align=\"\" height=\"60\" summary=\"\">\n<tbody>"
    extro = "</tbody>\n</table>\n</body>\n</html>"
    filename = op.join("html", reportname + ".htm")
    f = open(filename, "w")
    f.write(intro)
    for counter, cluster in enumerate(cluster_idx_list):
        num_of_clusters = len(cluster_idx_list)
        sdf_list = get_sdf_from_index_list(orig_sdf, cluster)
        img = Draw.MolsToGridImage(sdf_list, molsPerRow=4, legends=[m.GetProp(captionprop) for m in sdf_list])
        img_filename = "img/{}_{:03d}.png".format(reportname, counter)
        img.save("html/"+img_filename)
        hist_filename = "img/{}_{:03d}_hist.png".format(reportname, counter)
        show_hist(sdf_list, fields=[captionprop], show=False, savefig=False)
        pylab.savefig("html/"+hist_filename)
        f.write("  <tr><td>\n<br><b>Cluster {:03d}:</b></td><td></td></tr>\n<tr><td><img src=\"{}\" alt=\"icon\" /></td><td><img src=\"{}\" alt=\"icon\" /></td></tr>\n".format(counter, img_filename, hist_filename))

    f.write(extro)
    f.close()
