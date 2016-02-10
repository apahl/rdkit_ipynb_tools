#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
########
Pipeline
########

*Created on 3-Feb-2016 by A. Pahl*

A Pipelining Workflow using Python Generators, mainly for RDKit and large compound sets.

Example use:
    >>> import pipeline as p
    >>> s = p.Summary()
    >>> rd = p.start_csv_reader("/home/pahl/data_b64.csv.gz", summary=s)
    >>> b64 = p.pipe_mol_from_b64(rd, summary=s)
    >>> filt = p.pipe_mol_filter(b64, "[H]c2c([H])c1ncoc1c([H])c2C(N)=O", summary=s)
    >>> p.stop_sdf_writer(filt, "test.sdf", summary=s)

The progress of the pipeline can be followed in a terminal with: tail -f pipeline.log
"""


from __future__ import print_function, division


# import sys
# import base64
# import os.path as op
import time
from collections import Counter
import csv
import gzip
import pickle
import base64 as b64
# from copy import deepcopy

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
import rdkit.Chem.Descriptors as Desc
Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

from . import tools

try:
    from Contrib.SA_Score import sascorer
    SASCORER = True
except ImportError:
    print("* SA scorer not available. RDKit's Contrib dir needs to be in the Python import path...")
    SASCORER = False


def format_seconds(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    t_str = "{:02.0f}h {:02.0f}m {:02.2f}s".format(h, m, s)
    return t_str


class Summary(Counter):
    """A Counter-based class that keeps track of the time since its instantiation.
    Used for reporting running details of functions, esp. in the IPython notebook.

    Example:
        >>> s = Summary()
        >>> # then, in a loop:
        >>> for ....:
        ...     s["00 in"] += 1
        ...     # keeping track of failed entities:
        ...     if failed:
        ...         s["10 failed"] += 1
        ...     s["99 out"] += 1
        >>> # at the end of the function:
        >>> s.print()

    For fast-running loops with a lot of iterations, it could be sensible to use int counters
    and assign to Summary() outside of the loop at the end of the function."""

    def __init__(self, timeit=True, **kwargs):
        """Parameters:
            timeit: whether or not to use the timing functionality. Default: True"""
        super().__init__(**kwargs)
        self.timeit = timeit
        if self.timeit:
            self.t_start = time.time()


    def __str__(self):
        s_list = []
        mlen = max(map(len, self.keys()))
        sorted_keys = sorted(self.keys())
        line_end = "\n"
        for idx, k in enumerate(sorted_keys, 1):
            value = self[k]
            if self.timeit and idx == len(sorted_keys):
                line_end = ""
            if type(value) == float:
                s_list.append("{k:{mlen}s}: {val:10.2f}".format(k=k, mlen=mlen, val=value))
                s_list.append(line_end)
            else:
                s_list.append("{k:{mlen}s}: {val:>7}".format(k=k, mlen=mlen, val=value))
                s_list.append(line_end)

        if self.timeit:
            seconds = time.time() - self.t_start
            s_list.append("               (time: {})".format(format_seconds(seconds)))

        return "".join(s_list)


    def __repr__(self):
        return self.__str__()


    def print(self):
        """print the content of a dict or Counter object in a formatted way"""
        print(self.__str__())


def start_csv_reader(fn, max_records=0, summary=None, comp_id="start_csv_reader"):
    """A reader for csv files.

    Returns:
        An iterator with the fields as dict

    Parameters:
        fn (str): filename
        max_records (int): maximum number of records to read, 0 means all
        summary (Summary): a Counter class to collect runtime statistics
        comp_id: (str): the component Id to use for the summary"""

    if fn.lower() == "comas":
        fn = "/home/pahl/comas/share/export_data_b64.csv.gz"

    if ".gz" in fn:
        f = gzip.open(fn, mode="rt")
    else:
        f = open(fn)

    reader = csv.DictReader(f, dialect="excel-tab")
    prev_time = time.time()
    for rec_counter, row_dict in enumerate(reader, 1):
        if max_records > 0 and rec_counter > max_records: break
        if summary is not None:
            summary[comp_id] = rec_counter
            curr_time = time.time()
            if curr_time - prev_time > 2.0:  # write the log only every two seconds
                prev_time = curr_time
                print(summary, file=open("pipeline.log", "w"))
        yield row_dict

    f.close()
    if summary:
        print(summary, file=open("pipeline.log", "w"))


def start_cache_reader(name, summary=None, comp_id="start_cache_reader"):
    fn = "/tmp/{}".format(name)
    start_csv_reader(fn, summary=None, comp_id=comp_id)


def start_sdf_reader(fn, max_records=0, summary=None, comp_id="start_sdf_reader"):
    """A reader for SD files.

    Returns:
        An iterator with the fields as dict, including the molecule in the "mol" key

    Parameters:
        fn (str): filename
        max_records (int): maximum number of records to read, 0 means all
        summary (Summary): a Counter class to collect runtime statistics
        comp_id: (str): the component Id to use for the summary"""

    rec_counter = 0
    no_mol_counter = 0

    if ".gz" in fn:
        f = gzip.open(fn, mode="rt")
    else:
        f = open(fn, "rb")

    reader = Chem.ForwardSDMolSupplier(f)
    prev_time = time.time()
    for mol in reader:
        if max_records > 0 and rec_counter > max_records: break
        rec = {}
        rec_counter += 1
        if mol:
            for prop in mol.GetPropNames():
                rec[prop] = mol.GetProp(prop)
                mol.ClearProp(prop)

            rec["mol"] = mol

            if summary is not None:
                summary[comp_id] = rec_counter
                curr_time = time.time()
                if curr_time - prev_time > 2.0:  # write the log only every two seconds
                    prev_time = curr_time
                    print(summary)
                    print(summary, file=open("pipeline.log", "w"))

            yield rec

        else:
            no_mol_counter += 1
            if summary is not None:
                summary["{}_no_mol".format(comp_id)] = no_mol_counter

    f.close()
    if summary:
        print(summary, file=open("pipeline.log", "w"))


def stop_csv_writer(stream, fn, summary=None, comp_id="stop_csv_writer"):
    """Write CSV file from the incoming stream.

    Parameters:
        fn (str): filename
        summary (Summary): a Counter class to collect runtime statistics
        comp_id: (str): the component Id to use for the summary"""

    fields = []
    rec_counter = 0
    f = open(fn, "w")

    for rec in stream:
        if "mol" in rec:  # molecule object can not be written to CSV
            rec.pop("mol")

        if not fields:  # first record has to contain all fields of the pipeline.
            fields = list(rec.keys())
            writer = csv.DictWriter(f, fieldnames=fields, dialect="excel-tab")
            writer.writeheader()

        for key in rec.keys():  # remove keys that were not present in the first record
            if key not in fields:
                rec.pop(key)

        rec_counter += 1
        if summary is not None:
            summary[comp_id] = rec_counter

        writer.writerow(rec)

    f.close()
    print("* {}: {} records written.".format(comp_id, rec_counter))


def stop_sdf_writer(stream, fn, max=500, summary=None, comp_id="stop_sdf_writer"):
    """Write records in stream as SD File."""

    rec_counter = 0
    no_mol_counter = 0
    writer = Chem.SDWriter(fn)

    for rec in stream:
        if "mol" not in rec:
            no_mol_counter += 1
            if summary is not None:
                summary["{}_no_mol".format(comp_id)] = no_mol_counter
            return

        mol = rec["mol"]

        try:
            mol.GetConformer()
        except ValueError:  # no 2D coords... calculate them
            mol.Compute2DCoords()

        # assign the values from rec to the mol object
        for key in rec:
            if key == "mol": continue
            val = rec[key]
            if isinstance(val, str):
                if val != "":
                    mol.SetProp(key, val)
            else:
                mol.SetProp(key, str(val))

        rec_counter += 1
        if summary is not None:
            summary[comp_id] = rec_counter

        writer.write(mol)

        if rec_counter >= max: break

    writer.close()


def stop_mol_list_from_stream(stream, max=250, summary=None, comp_id="stop_mol_list_from_stream"):
    """Creates a Mol_list from the records in the stream. Stops the pipeline stream."""
    rec_counter = 0
    mol_list = tools.Mol_List()

    for rec in stream:
        if "mol" not in rec: continue

        mol = rec["mol"]

        try:
            mol.GetConformer()
        except ValueError:  # no 2D coords... calculate them
            mol.Compute2DCoords()

        # assign the values from rec to the mol object
        for key in rec:
            if key == "mol": continue
            val = rec[key]
            if isinstance(val, str):
                if val != "":
                    mol.SetProp(key, val)
            else:
                mol.SetProp(key, str(val))

        rec_counter += 1
        if summary is not None:
            summary[comp_id] = rec_counter

        mol_list.append(mol)

        if rec_counter >= max:
            break

    return mol_list


def stop_cache_writer(stream, name, summary=None, comp_id="stop_cache_writer"):
    """Write records in stream as cache."""

    fn = "/tmp/{}".format(name)

    stop_csv_writer(stream, fn, summary=summary, comp_id=comp_id)


def pipe_mol_from_smiles(stream, in_smiles="Smiles", remove=True, summary=None, comp_id="pipe_mol_from_smiles"):
    """Generate a molecule on the stream from Smiles."""
    rec_counter = 0
    for rec in stream:
        if in_smiles in rec:
            mol = Chem.MolFromSmiles(in_smiles)
            if remove:
                rec.pop(in_smiles)

            if mol:
                rec_counter += 1
                if summary is not None:
                    summary[comp_id] = rec_counter

                rec["mol"] = mol
                yield rec


def pipe_mol_from_b64(stream, in_b64="Mol_b64", remove=True, summary=None, comp_id="pipe_mol_from_b64"):
    """Generate a molecule on the stream from a b64 encoded mol object."""
    rec_counter = 0
    for rec in stream:
        if in_b64 in rec:
            mol = pickle.loads(b64.b64decode(rec[in_b64]))
            if remove:
                rec.pop(in_b64)

            if mol:
                rec_counter += 1
                if summary is not None:
                    summary[comp_id] = rec_counter

                rec["mol"] = mol
                yield rec


def pipe_mol_to_smiles(stream, out_smiles="Smiles", summary=None, comp_id="pipe_mol_to_smiles"):
    """Calculate Smiles from the mol object onthe stram."""
    rec_counter = 0
    for rec in stream:
        if "mol" in rec:
            rec[out_smiles] = Chem.MolToSmiles(rec["mol"])
            rec_counter += 1
            if summary is not None:
                summary[comp_id] = rec_counter
            yield rec


def pipe_mol_to_b64(stream, out_b64="Mol_b64", summary=None, comp_id="pipe_mol_to_b64"):
    rec_counter = 0
    for rec in stream:
        if "mol" in rec:
            mol = rec["mol"]
            if mol:
                rec["Mol_b64"] = b64.b64encode(pickle.dumps(mol)).decode()
                rec_counter += 1
                if summary is not None:
                    summary[comp_id] = rec_counter

                yield rec


def check_2d_coords(mol):
    """Check if a mol has 2D coordinates and if not, calculate them."""
    try:
        mol.GetConformer()
    except ValueError:  # no 2D coords... calculate them
        mol.Compute2DCoords()


def pipe_calc_props(stream, props, force2d=False, summary=None, comp_id="pipe_calc_props"):
    """Calculate properties from the Mol_List.
    props can be a single property or a list of properties.

    Calculable properties:
        2d, date, formula, hba, hbd, logp, molid, mw, rotb, sa (synthetic accessibility), tpsa

    Synthetic Accessibility (normalized):
        0: hard to synthesize; 1: easy access

        as described in:
            | Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions
            | *Peter Ertl and Ansgar Schuffenhauer*
            | Journal of Cheminformatics 1:8 (2009) (`link <http://www.jcheminf.com/content/1/1/8>`_)
    """

    rec_counter = 0
    for rec in stream:
        if not isinstance(props, list):
            props = [props]

        if "mol" in rec:

            mol = rec["mol"]

            if "2d" in props:
                if force2d:
                    mol.Compute2DCoords()
                else:
                    check_2d_coords(mol)

            if "date" in props:
                rec["date"] = time.strftime("%Y%m%d")

            if "formula" in props:
                rec["formula"] = Chem.CalcMolFormula(mol)

            if "hba" in props:
                rec["hba"] = str(Desc.NOCount(mol))

            if "hbd" in props:
                rec["hbd"] = str(Desc.NHOHCount(mol))

            if "logp" in props:
                rec["logp"] = "{:.2f}".format(Desc.MolLogP(mol))

            if "mw" in props:
                rec["mw"] = "{:.2f}".format(Desc.MolWt(mol))

            if "rotb" in props:
                mol.SetProp("rotb", str(Desc.NumRotatableBonds(mol)))

            if SASCORER and "sa" in props:
                score = sascorer.calculateScore(mol)
                norm_score = 1 - (score / 10)
                rec["sa"] = "{:.2f}".format(norm_score)

            if "tpsa" in props:
                rec["tpsa"] = str(int(Desc.TPSA(mol)))

            rec_counter += 1
            if summary is not None:
                summary[comp_id] = rec_counter

            yield rec


def pipe_custom_filter(stream, run_code, start_code=None, summary=None, comp_id="pipe_custom_filter"):
    """If the evaluation of run_code is true, the respective record will be put on the stream."""
    rec_counter = 0
    if start_code:
        eval(start_code)

    for rec in stream:
        if eval(run_code):
            rec_counter += 1
            if summary is not None:
                summary[comp_id] = rec_counter

            yield rec


def pipe_custom_man(stream, run_code, start_code=None, stop_code=None, summary=None, comp_id="pipe_custom_man"):
    """If the evaluation of run_code is true, the respective record will be put on the stream."""
    if start_code:
        eval(start_code)

    for rec in stream:
        eval(run_code)

        yield rec

    eval(stop_code)


def pipe_mol_filter(stream, smarts, invert=False, add_h=False, summary=None, comp_id="pipe_mol_filter"):
    rec_counter = 0
    query = Chem.MolFromSmarts(smarts)
    if not query:
        print("* {} ERROR: could not generate query from SMARTS.".format(comp_id))
        return None

    if "H" in smarts or "#1" in smarts:
        add_h = True

    for rec in stream:
        if "mol" not in rec: continue

        mol = rec["mol"]

        hit = False
        if add_h:
            mol_with_h = Chem.AddHs(mol)
            if mol_with_h.HasSubstructMatch(query):
                hit = True

        else:
            if mol.HasSubstructMatch(query):
                hit = True

        if invert:
            # reverse logic
            hit = not hit

        if hit:
            rec_counter += 1

            if summary is not None:
                summary[comp_id] = rec_counter

            yield rec


def pipe_remove_props(stream, props, summary=None, comp_id="pipe_remove_props"):
    """Remove properties from the stream.
    props can be a single property name or a list of property names."""

    if not isinstance(props, list):
        props = [props]

    rec_counter = 0
    for rec_counter, rec in enumerate(stream, 1):
        for prop in props:
            if prop in rec:
                rec.pop(prop)

        if summary is not None:
            summary[comp_id] = rec_counter

        yield rec


def pipe_rename_prop(stream, prop_old, prop_new, summary=None, comp_id="pipe_rename_prop"):
    """Rename a property on the stream.

    Parameters:
        prop_old (str): old name of the property
        prop_new (str): newname of the property"""

    rec_counter = 0
    for rec_counter, rec in enumerate(stream, 1):
        if prop_old in rec:
            rec[prop_new] = rec[prop_old]
            rec.pop(prop_old)

        if summary is not None:
            summary[comp_id] = rec_counter

        yield rec


def pipe_join_from_file(stream, join_on):
    pass
