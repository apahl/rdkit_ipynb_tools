#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
##########
Clustering
##########

*Created on Sun Feb 28 11:00 20165 by A. Pahl*

Clustering molecules.
"""

from copy import deepcopy
from collections import Counter, OrderedDict

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
# import rdkit.Chem.Descriptors as Desc
# import rdkit.Chem.Scaffolds.MurckoScaffold as MurckoScaffold
Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

# from PIL import Image, ImageChops

import numpy as np

# from . import html_templates as html
from . import tools, html_templates as html

# from ipywidgets import widgets
# from IPython.core.display import HTML, display


def renumber_clusters(cluster_list, start_at=1):
    """Renumber clusters in-place."""
    start_at -= 1
    # get the current individual cluster numbers present in the list
    id_list = sorted(set(tools.get_value(mol.GetProp("cluster_no"))
                         for mol in cluster_list.mols_with_prop("cluster_no")))

    # assign the new ids as values the old id's keys
    new_ids = {k: v for v, k in enumerate(id_list, 1 + start_at)}
    for mol in cluster_list:
        if not mol.HasProp("cluster_no"): continue
        old_id = int(mol.GetProp("cluster_no"))
        mol.SetProp("cluster_no", str(new_ids[old_id]))


def get_clusters_by_no(cluster_list, cl_no, make_copy=True, renumber=False):
    """Return one or more clusters (provide numbers as list) by their number."""
    if not isinstance(cl_no, list):
        cl_no = [cl_no]

    cluster = tools.Mol_List()
    for mol in cluster_list:
        if mol.HasProp("cluster_no") and int(mol.GetProp("cluster_no")) in cl_no:
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
        if mol.HasProp("cluster_no") and int(mol.GetProp("cluster_no")) not in cl_no:
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

        if not mol.HasProp("cluster_no"): continue
        cl_id = int(mol.GetProp("cluster_no"))
        ctr[cl_id] += 1

    # now, only keep mols which belong to clusters of the desired length,
    # including the cores
    for mol in cluster_list:
        if not mol.HasProp("cluster_no"): continue
        cl_id = int(mol.GetProp("cluster_no"))
        if ctr[cl_id] >= min_len and ctr[cl_id] <= max_len:
            result_list.append(mol)

    if renumber:
        renumber_clusters(result_list)

    return result_list


def get_cores(cluster_list):
    """Find and return the core molecules in a cluster_list."""
    core_list = cluster_list.has_prop_filter("is_core")
    return core_list


def get_members(cluster_list):
    """Find and return the members of a cluster_list, exclude the cores."""
    member_list = cluster_list.has_prop_filter("is_core", invert=True)
    return member_list


def get_stats_for_cluster(cluster_list, activity_prop):
    stats = {}
    value_list = [tools.get_value(mol.GetProp(activity_prop)) for mol in cluster_list.mols_with_prop(activity_prop)]
    stats["num_values"] = len(value_list)
    stats["min"] = min(value_list) if value_list else None
    stats["max"] = max(value_list) if value_list else None
    stats["mean"] = np.mean(value_list) if value_list else None
    stats["median"] = np.median(value_list) if value_list else None

    return stats


def get_clusters_with_activity(cluster_list, activity_prop, min_act=None, max_act=None, min_len=1, renumber=False):
    """Return only the clusters which fulfill the given activity criteria.
    Parameters:
        min_act (str): something like `< 50` or `> 70` that can be evaluated.
        max_act (str): see above."""

    if min_act is not None:
        min_act_comp = compile('stats["min"] {}'.format(min_act), '<string>', 'eval')
    if max_act is not None:
        max_act_comp = compile('stats["max"] {}'.format(max_act), '<string>', 'eval')

    cores_and_members = get_cores(cluster_list)
    if len(cores_and_members) > 0:
        cores_and_members = cores_and_members.prop_filter('num_members >= {}'.format(min_len))
    members_all = get_members(cluster_list)
    tmp_list = tools.Mol_List()

    cl_ids = sorted(set(int(mol.GetProp("cluster_no")) for mol in members_all.mols_with_prop("cluster_no")))
    for new_id, cl_id in enumerate(cl_ids, 1):
        cluster = get_clusters_by_no(members_all, cl_id)
        if len(cluster) < min_len: continue
        stats = get_stats_for_cluster(cluster, activity_prop)
        if stats["num_values"] == 0: continue
        stats["min"]  # to quiet the linter
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
    """Find and add cores to the cluster_list in-place."""

    # first, remove any already existing cores
    members_all = get_members(cluster_list)

    # find cluster numbers
    cl_ids = set(tools.get_value(mol.GetProp("cluster_no")) for mol in members_all.mols_with_prop("cluster_no"))

    for cl_id in cl_ids:
        cluster = get_clusters_by_no(members_all, cl_id, make_copy=False)
        # determine the cluster core by MCSS

        # do not generate cores for singletons and pairs
        if len(cluster) < 3:
            continue

        core_mol = tools.find_mcs(cluster)
        if core_mol is None:
            continue

        tools.check_2d_coords(core_mol)

        # set a number of properties for the cluster core
        id_prop = tools.guess_id_prop(cluster[0].GetPropNames())
        core_mol.SetProp(id_prop, str(cl_id))
        core_mol.SetProp("is_core", "yes")
        core_mol.SetProp("cluster_no", str(cl_id))
        core_mol.SetProp("num_members", str(len(cluster)))

        if activity_prop is not None:
            stats = get_stats_for_cluster(cluster, activity_prop)
            core_mol.SetProp("num_values", str(stats["num_values"]))

            core_mol.SetProp("min", "{:.2f}".format(stats["min"]))
            core_mol.SetProp("max", "{:.2f}".format(stats["max"]))
            core_mol.SetProp("mean", "{:.2f}".format(stats["mean"]))
            core_mol.SetProp("median", "{:.2f}".format(stats["median"]))

        if align_to_core:
            # align_mol = MurckoScaffold.GetScaffoldForMol(core_mol)
            cluster.align(core_mol)

        members_all.insert(0, core_mol)

    return members_all


def get_mol_list_from_index_list(orig_sdf, index_list, cl_id):
    """generate sdf_lists after clustering"""
    cluster_list = tools.Mol_List()
    for x in index_list:
        mol = deepcopy(orig_sdf[x])
        mol.SetProp("cluster_no", str(cl_id))
        cluster_list.append(mol)

    if len(cluster_list) == 2:
        cluster_list[0].SetProp("is_pair", "yes")
        cluster_list[1].SetProp("is_pair", "yes")
    elif len(index_list) == 1:
        cluster_list[0].SetProp("is_single", "yes")

    return cluster_list


def cluster_from_mol_list(mol_list, cutoff=0.8, activity_prop=None,
                          summary_only=True, generate_cores=False, align_to_core=False):
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

    if summary_only:
        return None

    cluster_list = tools.Mol_List()

    # go over each list of indices to collect the cluster's molecules
    for cl_id, idx_list in enumerate(sorted(cluster_idx_list, key=len, reverse=True), 1):
        cluster_list.extend(get_mol_list_from_index_list(mol_list, idx_list, cl_id))

    if generate_cores:
        cluster_list = add_cores(cluster_list, activity_prop, align_to_core)

    return cluster_list


def show_numbers(cluster_list):
    """Show some numbers for the cluster_list."""
    all_members = get_members(cluster_list)
    ctr_cl_no = Counter([int(mol.GetProp("cluster_no")) for mol in all_members])
    ctr_size = Counter([s for s in ctr_cl_no.values()])

    total = 0
    print("\ncluster size  |  number of clusters")
    print("------------- + -------------------")
    for i in sorted(ctr_size, reverse=True):
        print("     {:3d}      |   {:3d}".format(i, ctr_size[i]))
        total = total + (i * ctr_size[i])

    print()
    print("Number of compounds as sum of members per cluster size:", total)


def write_report(cluster_list, fn="clusters.html", title="Clusters", activity_prop=None):
    content = []

    # collect the cluster numbers in the order in which they are in the cluster_list:
    cluster_numbers = OrderedDict()
    for mol in cluster_list.has_prop_filter("cluster_no"):
        cluster_numbers[int(mol.GetProp("cluster_no"))] = 0  # dummy value

    # print(cluster_numbers)
    # return cluster_numbers
    reverse = False
    if activity_prop is not None:
        if activity_prop and not isinstance(activity_prop, list):
            activity_prop = [activity_prop]
        s = activity_prop[0].lower()
        if "pic50" in s or "pec50" in s:  # determine within the cluster
            reverse = True

    for cl_no in cluster_numbers:
        cluster = get_clusters_by_no(cluster_list, cl_no)
        if len(cluster) == 0: continue
        content.append("<h2>Cluster {}</h2>".format(cl_no))
        core = get_cores(cluster)
        if len(core) > 0:
            core.remove_props(["Compound_Id", "is_core", "num_values"])
            core.order = ["cluster_no", "num_members", "min", "max", "mean", "median"]
            content.append('<h4>Core:</h4>')
            content.append(core.table(raw=True))
            content.append("<h4>Members:</h4>")

        members = get_members(cluster)
        if activity_prop is not None:
            members.sort_list(activity_prop[0], reverse=reverse)
            content.append(members.grid(activity_prop, raw=True))
        else:
            content.append(members.grid(raw=True))

        content.append("<hr>")

    page = html.page("\n".join(content), title=title)
    html.write(page, fn=fn)
