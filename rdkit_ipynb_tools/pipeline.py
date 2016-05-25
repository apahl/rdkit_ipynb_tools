#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
########
Pipeline
########

*Created on 3-Feb-2016 by A. Pahl*

A Pipelining Workflow using Python Generators, mainly for RDKit and large compound sets.

Example use:
    >>> from rdkit_ipynb_tools import pipeline as p
    >>> s = p.Summary()
    >>> rd = p.start_csv_reader("/home/pahl/data_b64.csv.gz", summary=s)
    >>> b64 = p.pipe_mol_from_b64(rd, summary=s)
    >>> filt = p.pipe_mol_filter(b64, "[H]c2c([H])c1ncoc1c([H])c2C(N)=O", summary=s)
    >>> p.stop_sdf_writer(filt, "test.sdf", summary=s)

The progress of the pipeline can be followed in a terminal with: `watch -n 2 less pipeline.log`
"""


# import sys
# import os.path as op
from copy import deepcopy
import time
from collections import OrderedDict, defaultdict
import csv
import gzip
import pickle
import base64 as b64
import tempfile

import numpy as np

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
import rdkit.Chem.Descriptors as Desc
Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

from . import tools

try:
    from rdkit.Avalon import pyAvalonTools as pyAv
    USE_AVALON = True
except ImportError:
    USE_AVALON = False

try:
    from Contrib.SA_Score import sascorer
    SASCORER = True
except ImportError:
    print("* SA scorer not available. RDKit's Contrib dir needs to be in the Python import path...")
    SASCORER = False

try:
    # interactive IPython session
    IPY = get_ipython()

except NameError:
    iPY = False


def format_seconds(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    t_str = "{:02.0f}h {:02.0f}m {:02.2f}s".format(h, m, s)
    return t_str


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


class Summary(OrderedDict):
    """An OrderedDict-based class that keeps track of the time since its instantiation.
    Used for reporting running details of pipeline functions."""

    def __init__(self, timeit=True, **kwargs):
        """Parameters:
            timeit: whether or not to use the timing functionality. Default: True"""
        super().__init__(**kwargs)
        self.timeit = timeit
        if self.timeit:
            self.t_start = time.time()


    def __str__(self):
        s_list = []
        keys = self.keys()
        mlen = max(map(len, keys))
        line_end = "\n"
        for idx, k in enumerate(keys, 1):
            value = self[k]
            if self.timeit and idx == len(keys):
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


def pipe(val, *forms):
    """Inspired by the thread_first function of the `Toolz <https://pypi.python.org/pypi/toolz>`_ project
    and adapted to also accept keyword arguments. Removed the discouraged reduce function.
    If functions of the pipeline nedd additional parameters, the function and
    the parameters have to be passed as tuples. Keyword arguments have to be
    passed as dicts in these tuples:

    >>> s = Summary()
    >>> rd = start_sdf_reader("test.sdf", summary=s)
    >>> pipe(rd,
    >>>      pipe_keep_largest_fragment,
    >>>      (pipe_neutralize_mol, {"summary": s}),
    >>>      (pipe_keep_props, ["Ordernumber", "NP_Score"]),
    >>>      (stop_csv_writer, "test.csv", {"summary": s})
    >>>     )"""

    def evalform_front(val, form):
        if callable(form):
            return form(val)
        if isinstance(form, tuple):
            args = [val]
            kwargs = {}
            func = form[0]
            for a in form[1:]:
                if isinstance(a, dict):
                    kwargs.update(a)
                else:
                    args.append(a)
            return func(*args, **kwargs)

    result = val
    for form in forms:
        result = evalform_front(result, form)

    return result


def start_csv_reader(fn, max_records=0, summary=None, comp_id="start_csv_reader"):
    """A reader for csv files.

    Returns:
        An iterator with the fields as dict

    Parameters:
        fn (str): filename
        max_records (int): maximum number of records to read, 0 means all
        summary (Summary): a Counter class to collect runtime statistics
        comp_id: (str): the component Id to use for the summary"""

    if ".gz" in fn:
        f = gzip.open(fn, mode="rt")
    else:
        f = open(fn)

    reader = csv.DictReader(f, dialect="excel-tab")
    prev_time = time.time()
    for rec_counter, row_dict in enumerate(reader, 1):
        if max_records > 0 and rec_counter > max_records: break
        # make a copy with non-empty values
        rec = {k: get_value(v) for k, v in row_dict.items() if v is not None and v != ""}  # make a copy with non-empty values

        if summary is not None:
            summary[comp_id] = rec_counter
            curr_time = time.time()
            if curr_time - prev_time > 2.0:  # write the log only every two seconds
                prev_time = curr_time
                print(summary, file=open("pipeline.log", "w"))
        yield rec

    f.close()
    if summary:
        print(summary, file=open("pipeline.log", "w"))
        print(summary)



def start_cache_reader(name, summary=None, comp_id="start_cache_reader"):
    fn = "/tmp/{}".format(name)
    start_csv_reader(fn, summary=None, comp_id=comp_id)


def start_sdf_reader(fn, max_records=0, summary=None, comp_id="start_sdf_reader"):
    """A reader for SD files.

    Returns:
        An iterator with the fields as dict, including the molecule in the "mol" key

    Parameters:
        fn (str, list): filename or list of filenames
        max_records (int): maximum number of records to read, 0 means all
        summary (Summary): a Counter class to collect runtime statistics
        comp_id: (str): the component Id to use for the summary"""

    rec_counter = 0
    no_mol_counter = 0
    # also open lists of files
    if not isinstance(fn, list):
        fn = [fn]

    for filen in fn:
        if ".gz" in filen:
            f = gzip.open(filen, mode="rt")
        else:
            f = open(filen, "rb")

        reader = Chem.ForwardSDMolSupplier(f)
        prev_time = time.time()
        for mol in reader:
            if max_records > 0 and rec_counter > max_records: break
            rec = {}
            rec_counter += 1
            if mol:
                for prop in mol.GetPropNames():
                    val = mol.GetProp(prop)
                    if len(val) > 0:  # transfer only those properties to the stream which carry a value
                        rec[prop] = get_value(val)
                    mol.ClearProp(prop)

                rec["mol"] = mol

                if summary is not None:
                    summary[comp_id] = rec_counter
                    curr_time = time.time()
                    if curr_time - prev_time > 2.0:  # write the log only every two seconds
                        prev_time = curr_time
                        print(summary, file=open("pipeline.log", "w"))

                yield rec

            else:
                no_mol_counter += 1
                if summary is not None:
                    summary["{}_no_mol".format(comp_id)] = no_mol_counter

        f.close()

    if summary:
        print(summary, file=open("pipeline.log", "w"))
        print(summary)


def start_stream_from_dict(d, summary=None, comp_id="start_stream_from_dict", show_first=False):
    """Provide a data stream from a dict."""
    prev_time = time.time()
    d_keys = list(d.keys())
    l = len(d[d_keys[0]])  # length of the stream
    rec_counter = 0
    for idx in range(l):
        rec = {}
        for k in d_keys:
            rec[k] = d[k][idx]

        rec_counter += 1
        if summary is not None:
            summary[comp_id] = rec_counter
            curr_time = time.time()
            if curr_time - prev_time > 2.0:  # write the log only every two seconds
                prev_time = curr_time
                print(summary)
                print(summary, file=open("pipeline.log", "w"))

        if show_first and rec_counter == 1:
            print("{}:".format(comp_id), rec)

        yield rec


def start_stream_from_mol_list(mol_list, summary=None, comp_id="start_stream_from_mol_list"):
    """Provide a data stream from a Mol_List."""
    prev_time = time.time()
    rec_counter = 0
    for mol in mol_list:
        if not mol: continue
        rec = {}
        props = mol.GetPropNames()
        for prop in props:
            val = get_value(mol.GetProp(prop))
            mol.ClearProp(prop)
            if val is not None:
                rec[prop] = val

        rec["mol"] = mol

        rec_counter += 1
        if summary is not None:
            summary[comp_id] = rec_counter
            curr_time = time.time()
            if curr_time - prev_time > 2.0:  # write the log only every two seconds
                prev_time = curr_time
                print(summary)
                print(summary, file=open("pipeline.log", "w"))

        yield rec

    if summary:
        print(summary, file=open("pipeline.log", "w"))
        print(summary)


def stop_csv_writer(stream, fn, summary=None, comp_id="stop_csv_writer"):
    """Write CSV file from the incoming stream.

    Parameters:
        fn (str): filename
        summary (Summary): a Counter class to collect runtime statistics
        comp_id: (str): the component Id to use for the summary"""

    fields = OrderedDict()
    rec_counter = 0
    tmp = tempfile.TemporaryFile("w+")

    for rec in stream:
        if "mol" in rec:  # molecule object can not be written to CSV
            rec.pop("mol")

        line = []
        cp = rec.copy()

        # first write the records whose keys are already in fields:
        for key in fields:
            if key in cp:
                val = cp[key]
                if val is None: val = ""
                line.append(str(val))
                cp.pop(key)
            else:
                line.append("")

        # now collect the additional records (and add them to fields)
        for key in cp:
            fields[key] = 0  # dummy value
            val = cp[key]
            if val is None: val = ""
            line.append(str(val))

        tmp.write("\t".join(line) + "\n")
        rec_counter += 1
        if summary is not None:
            summary[comp_id] = rec_counter

    f = open(fn, "w")
    num_columns = len(fields)
    first_line = True
    tmp.seek(0)
    for line_str in tmp:
        if first_line:  # write the final header
            first_line = False
            line = list(fields.keys())
            f.write("\t".join(line) + "\n")

        line = line_str.rstrip("\n").split("\t")

        # fill up all lines with empty records to the number of columns
        num_fill_records = num_columns - len(line)
        fill = [""] * num_fill_records
        line.extend(fill)


        f.write("\t".join(line) + "\n")

    f.close()
    tmp.close()
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


def stop_count_records(stream, summary=None, comp_id="stop_count_records"):
    """Only count the records from the incoming stream."""
    rec_counter = 0

    for rec in stream:

        rec_counter += 1
        if summary is not None:
            summary[comp_id] = rec_counter

    return rec_counter


def stop_cache_writer(stream, name, summary=None, comp_id="stop_cache_writer"):
    """Write records in stream as cache."""

    fn = "/tmp/{}".format(name)

    stop_csv_writer(stream, fn, summary=summary, comp_id=comp_id)


def pipe_mol_from_smiles(stream, in_smiles="Smiles", remove=True, summary=None, comp_id="pipe_mol_from_smiles"):
    """Generate a molecule on the stream from Smiles."""
    rec_counter = 0
    for rec in stream:
        if in_smiles in rec:
            mol = Chem.MolFromSmiles(rec[in_smiles])
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


def start_mol_csv_reader(fn, max_records=0, in_b64="Mol_b64", summary=None, comp_id="start_mol_csv_reader"):
    """A reader for csv files containing molecules in binary b64 format.

    Returns:
        An iterator with the fields and the molecule as dict

    Parameters:
        fn (str): filename
        max_records (int): maximum number of records to read, 0 means all
        summary (Summary): a Counter class to collect runtime statistics
        comp_id: (str): the component Id to use for the summary"""

    rd = start_csv_reader(fn, max_records, summary, comp_id)
    mol = pipe_mol_from_b64(rd, in_b64)

    return mol


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


def pipe_calc_props(stream, props, force2d=False, summary=None, comp_id="pipe_calc_props"):
    """Calculate properties from the Mol_List.
    props can be a single property or a list of properties.

    Calculable properties:
        2d, date, formula, hba, hbd, logp, molid, mw, smiles, rotb, sa (synthetic accessibility), tpsa

    Synthetic Accessibility (normalized):
        0: hard to synthesize; 1: easy access

        as described in:
            | Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions
            | *Peter Ertl and Ansgar Schuffenhauer*
            | Journal of Cheminformatics 1:8 (2009) (`link <http://www.jcheminf.com/content/1/1/8>`_)
    """

    rec_counter = 0
    if not isinstance(props, list):
        props = [props]

    # make all props lower-case:
    props = list(map(lambda x: x.lower(), props))

    for rec in stream:
        if "mol" in rec:
            mol = rec["mol"]
            if "2d" in props:
                check_2d_coords(mol, force2d)

            if "date" in props:
                rec["Date"] = time.strftime("%Y%m%d")

            if "formula" in props:
                rec["Formula"] = Chem.CalcMolFormula(mol)

            if "hba" in props:
                rec["HBA"] = str(Desc.NOCount(mol))

            if "hbd" in props:
                rec["HBD"] = str(Desc.NHOHCount(mol))

            if "logp" in props:
                rec["LogP"] = "{:.2f}".format(Desc.MolLogP(mol))

            if "mw" in props:
                rec["MW"] = "{:.2f}".format(Desc.MolWt(mol))

            if "rotb" in props:
                mol.SetProp("RotB", str(Desc.NumRotatableBonds(mol)))

            if "smiles" in props:
                mol.SetProp("Smiles", Chem.MolToSmiles(mol))

            if SASCORER and "sa" in props:
                score = sascorer.calculateScore(mol)
                norm_score = 1 - (score / 10)
                rec["SA"] = "{:.2f}".format(norm_score)

            if "tpsa" in props:
                rec["TPSA"] = str(int(Desc.TPSA(mol)))

            rec_counter += 1
            if summary is not None:
                summary[comp_id] = rec_counter

            yield rec


def pipe_custom_filter(stream, run_code, start_code=None, summary=None, comp_id="pipe_custom_filter"):
    """If the evaluation of run_code is true, the respective record will be put on the stream."""
    rec_counter = 0
    if start_code is not None:
        exec(start_code)

    # pre-compile the run_code statement for performance reasons
    byte_code = compile(run_code, '<string>', 'eval')
    for rec in stream:
        if eval(byte_code):
            rec_counter += 1
            if summary is not None:
                summary[comp_id] = rec_counter

            yield rec


def pipe_custom_man(stream, run_code, start_code=None, stop_code=None, comp_id="pipe_custom_man"):
    """If the evaluation of run_code is true, the respective record will be put on the stream."""
    if start_code is not None:
        exec(start_code)

    byte_code = compile(run_code, '<string>', 'exec')
    for rec in stream:
        exec(byte_code)

        yield rec

    if stop_code:
        exec(stop_code)


def pipe_has_prop_filter(stream, prop, invert=False, summary=None, comp_id="pipe_has_prop_filter"):
    rec_counter = 0

    for rec in stream:

        hit = prop in rec

        if invert:
            # reverse logic
            hit = not hit

        if hit:
            rec_counter += 1

            if summary is not None:
                summary[comp_id] = rec_counter

            yield rec


def pipe_id_filter(stream, cpd_ids, id_prop="Compound_Id", summary=None, comp_id="pipe_id_filter"):
    rec_counter = 0
    if not isinstance(cpd_ids, list):
        cpd_ids = [cpd_ids]

    cpd_ids = {c_id: 0 for c_id in cpd_ids}

    for rec in stream:
        if id_prop not in rec: continue

        if rec[id_prop] in cpd_ids:
            rec_counter += 1

            if summary is not None:
                summary[comp_id] = rec_counter

            yield rec


def pipe_mol_filter(stream, query, smarts=False, invert=False, add_h=False, summary=None, comp_id="pipe_mol_filter"):
    rec_counter = 0
    if "[H]" in query or "#1" in query:
        add_h = True

    if add_h or "#6" in query or "#7" in query:
        smarts = True

    query_mol = Chem.MolFromSmarts(query) if smarts else Chem.MolFromSmiles(query)
    if not query_mol:
        print("* {} ERROR: could not generate query from SMARTS.".format(comp_id))
        return None

    for rec in stream:
        if "mol" not in rec: continue

        mol = rec["mol"]

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


def pipe_keep_props(stream, props, summary=None, comp_id="pipe_keep_props", show_first=False):
    """Keep only the listed properties on the stream. "mol" is always kept by this component.
    props can be a single property name or a list of property names.
    show_first prints the first records for debugging purposes."""

    if not isinstance(props, list):
        props = [props]

    if "mol" not in props:
        props.append("mol")

    for rec_counter, rec in enumerate(stream, 1):
        for prop in rec.copy().keys():
            if prop not in props:
                rec.pop(prop)

        if summary is not None:
            summary[comp_id] = rec_counter

        if show_first and rec_counter == 1:
            print("{}:".format(comp_id), rec)

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


def pipe_join_data_from_file(stream, fn, join_on, behaviour="joined_only",
                             decimals=2, summary=None, comp_id="pipe_join_data_from_file", show_first=False):
    """Joins data from a csv or SD file.
    CAUTION: The input stream will be held in memory by this component!

    Parameters:
        stream (dict iterator): stream of input compounds.
        fn (str): name of the file (type is determined by having "sdf" in the name or not).
        join_on (str): property to join on
        behaviour (str):
            "joined_only": only put those recored on the stream on which data was joined (default).
            "keep_all": put all input records on the stream again, including those, on which no data was joined.
        decimals (int): number of decimal places for floating point values. Default: 2."""

    # collect the records from the stream in a list, store the position of the join_on properties in a dict
    stream_counter = -1
    stream_list = []
    stream_dict = {}  # dict to hold the join_on properties and their positions in the stream_list
    joined_stream_list_idx = set()  # hold the indexes in stream_list, on which data was joined

    for rec in stream:
        stream_join_on_val = rec.get(join_on, False)
        if stream_join_on_val is False: continue
        stream_counter += 1
        stream_list.append(rec)
        stream_dict[stream_join_on_val] = stream_counter

    if "sdf" in fn:
        rd = start_sdf_reader(fn)
    else:
        rd = start_csv_reader(fn)

    rec_counter = 0
    for rec in rd:
        rec_join_on_val = rec.get(join_on, False)
        if not rec_join_on_val: continue
        stream_join_on_idx = stream_dict.get(rec_join_on_val, False)
        if stream_join_on_idx is False: continue

        joined_stream_list_idx.add(stream_join_on_idx)
        stream_rec = stream_list[stream_join_on_idx]
        for k in stream_rec:
            if k != join_on:
                rec[k] = stream_rec[k]

        rec_counter += 1
        if summary is not None:
            summary[comp_id] = rec_counter

        if show_first and rec_counter == 1:
            print("{}:".format(comp_id), rec)

        yield rec

    # with behaviour="keep_all", now add the records to the stream on which no data was joined.
    if "all" in behaviour.lower():
        for idx, rec in enumerate(stream_list):
            if idx in joined_stream_list_idx: continue  # this record was already output with joined data

            rec_counter += 1
            if summary is not None:
                summary[comp_id] = rec_counter

            yield rec


def pipe_keep_largest_fragment(stream, summary=None, comp_id="pipe_keep_largest_frag"):
    rec_counter = 0
    frag_counter = 0
    for rec in stream:
        if "mol" not in rec: continue
        mol = rec["mol"]
        if not mol: continue

        mols = Chem.GetMolFrags(mol, asMols=True)
        if len(mols) > 1:
            frag_counter += 1
            mols = sorted(mols, key=Desc.HeavyAtomCount, reverse=True)
            if summary is not None:
                summary["{}_has_frags".format(comp_id)] = frag_counter

        mol = mols[0]

        rec["mol"] = mol

        rec_counter += 1
        if summary is not None:
            summary[comp_id] = rec_counter

        yield rec


def pipe_neutralize_mol(stream, summary=None, comp_id="pipe_neutralize_mol"):
    pattern = (
        # Imidazoles
        ('[n+;H]', 'n'),
        # Amines
        ('[N+;!H0]', 'N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]', 'O'),
        # Thiols
        ('[S-;X1]', 'S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]', 'N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]', 'N'),
        # Tetrazoles
        ('[n-]', '[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]', 'S'),
        # Amides
        ('[$([N-]C=O)]', 'N'),
    )

    reactions = [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y, False)) for x, y in pattern]

    rec_counter = 0
    neutr_counter = 0
    for rec in stream:
        if "mol" not in rec: continue
        mol = rec["mol"]
        if not mol: continue

        replaced = False
        for reactant, product in reactions:
            while mol.HasSubstructMatch(reactant):
                replaced = True
                rms = Chem.ReplaceSubstructs(mol, reactant, product)
                mol = rms[0]

        if replaced:
            Chem.SanitizeMol(mol)
            mol.Compute2DCoords()

        rec_counter += 1
        if summary is not None:
            summary[comp_id] = rec_counter
            if replaced:
                neutr_counter += 1
                summary["{}_neutralized".format(comp_id)] = neutr_counter

        rec["mol"] = mol

        yield rec


def pipe_inspect_stream(stream, fn="pipe_inspect.txt", exclude=None, summary=None):
    """Write records from the stream into the file `fn` every two seconds.
    Do not write records from the exclude list."""
    prev_time = time.time()
    if exclude is not None:
        if not isinstance(exclude, list):
            exclude = [exclude]

    for rec in stream:
        if exclude is not None:
            rec = deepcopy(rec)
            for prop in exclude:
                rec.pop(prop, None)

            curr_time = time.time()
            if curr_time - prev_time > 2.0:  # write the log only every two seconds
                prev_time = curr_time
                if summary is not None:
                    print(summary, "\n\n", rec, file=open(fn, "w"))
                else:
                    print(rec, file=open(fn, "w"))

        yield rec


def pipe_merge_data(stream, merge_on, str_props="concat", num_props="mean", summary=None, comp_id="pipe_merge_data"):
    """Merge the data from the stream on the `merge_on` property.
    WARNING: The stream is collected in memory by this component!

    Parameters:
        merge_on (str): Name of the property (key) to merge on.
        str_props (str): Merge behaviour for string properties.
            Allowed values are: concat ("; "-separated concatenation),
            unique ("; "-separated concatenation of the unique values),
            keep_first, keep_last.
        num_props (str): Merge behaviour for numerical values.
            Allowed values are: mean, median, keep_first, keep_last."""

    def _get_merged_val_from_val_list(val_list, str_props, num_props):
        if isinstance(val_list[0], str):
            if "concat" in str_props:
                return "; ".join(val_list)
            if "unique" in str_props:
                return "; ".join(set(val_list))
            if "first" in str_props:
                return val_list[0]
            if "last" in str_props:
                return val_list[-1]

            return val_list[0]

        else:
            if "mean" in num_props:
                return np.mean(val_list)
            if "median" in num_props:
                return np.median(val_list)
            if "first" in num_props:
                return val_list[0]
            if "last" in num_props:
                return val_list[-1]

            return val_list[0]


    merged = defaultdict(lambda: defaultdict(list))  # defaultdict of defaultdict(list)
    if summary is not None:
        summary[comp_id] = "collecting..."

    for rec in stream:
        if merge_on not in rec: continue

        merge_on_val = rec.pop(merge_on)
        for prop in rec.keys():
            merged[merge_on_val][prop].append(rec[prop])

    rec_counter = 0
    for item in merged:
        rec = {merge_on: item}
        for prop in merged[item]:
            val_list = merged[item][prop]
            rec[prop] = _get_merged_val_from_val_list(val_list, str_props, num_props)

        rec_counter += 1
        if summary is not None:
            summary[comp_id] = rec_counter

        yield rec


def dict_from_csv(fn, max_records=0):
    """Read a CSV file and return a dict with the headers a keys and the columns as value lists.
    Empty cells are np.nan."""

    d = defaultdict(list)

    if ".gz" in fn:
        f = gzip.open(fn, mode="rt")
    else:
        f = open(fn)

    reader = csv.DictReader(f, dialect="excel-tab")

    for rec_counter, row_dict in enumerate(reader, 1):
        for k in row_dict:
            v = row_dict[k]
            if v == "" or v is None:
                d[k].append(np.nan)
            else:
                d[k].append(get_value(v))

        if max_records > 0 and rec_counter >= max_records: break

    print("  > {} records read".format(rec_counter))

    return d


def generate_pipe_from_csv(fn):
    """Generate a valid pipeline from a formatted csv file (see examples/example_pipe.ods)."""

    f = open(fn)
    reader = list(csv.DictReader(f, dialect="excel-tab"))
    num_of_lines = len(reader)
    pipe_list = ["s = p.Summary()\n"]
    for line_no, row_dict in enumerate(reader, 1):
        # clean up the field
        for k in row_dict:
            # replace the weird quotation marks that my Libreoffice exports:
            if "”" in row_dict[k]:
                row_dict[k] = row_dict[k].replace("”", '"')
            if "“" in row_dict[k]:
                row_dict[k] = row_dict[k].replace("“", '"')

        if row_dict["Summary"]:
            if row_dict["KWargs"]:
                row_dict["KWargs"] = row_dict["KWargs"] + ", 'summary': s"
            else:
                row_dict["KWargs"] = "'summary': s"

        if line_no == 1:
            if row_dict["KWargs"]:
                pipe_list.append("rd = p.{Component}({Args}, **{{{KWargs}}})\n".format(**row_dict))
            else:
                pipe_list.append("rd = p.{Component}({Args})\n".format(**row_dict))
            pipe_list.append("res = p.pipe(\n    rd,\n")
            continue

        if row_dict["Args"] or row_dict["KWargs"]:
            pipe_list.append("    (p.{}".format(row_dict["Component"]))
            if row_dict["Args"]:
                pipe_list.append(", {}".format(row_dict["Args"]))
            if row_dict["KWargs"]:
                pipe_list.append(', {{{}}}'.format(row_dict["KWargs"]))
            pipe_list.append(')')
        else:
            pipe_list.append('    p.{}'.format(row_dict["Component"]))

        if line_no < num_of_lines:
            pipe_list.append(",\n")
        else:
            pipe_list.append("\n")

    pipe_list.append(')')

    pipe_str = "".join(pipe_list)
    if IPY:
        IPY.set_next_input(pipe_str)
    else:
        print(pipe_str)
