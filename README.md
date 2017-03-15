# RDKit IPython Tools
by Axel Pahl

A set of tools to use with the Open Source Cheminformatics toolkit
[RDKit](http.//www.rdkit.org) in the Jupyter Notebook.<br>
Written for Python3, only tested on Linux (Ubuntu 16.04)
and the conda install of the RDkit.

# Module tools

A Mol_List class was introduced, which is a subclass of a Python list for holding lists of RDKit molecule objects and allows direct access to a lot of the RDKit functionality.
It is meant to be used with the Jupyter Notebook and includes a.o.:
* display of the Mol_List
    * as HTML table, nested table or grid
* display of a summary including number of records and min, max, mean, median for numeric properties
* display of correlations between the Mol_List's properties
  (using np.corrcoef, this allows getting a quick overview on which properties correlate with each other)
* methods for sorting, searching (by property or substructure) and filtering the Mol_List
* methods for renaming, reordering and calculating properties
* direct plotting of properties as publication-grade [Highcharts](http://www.highcharts.com/) *or* [Bokeh](http://bokeh.pydata.org/en/latest/) plots with **structure tooltips** (!).
    * the plotting functionalities reside in their own module and can also be used for plotting Pandas dataframes and Python dicts.
    * further development will focus on Bokeh because of the more pythonic interface


## Other functions in the tools module:
- *jsme*: Display Peter Ertl's [Javascript Molecule Editor](http://peter-ertl.com/jsme/) to enter a molecule directly in the IPython notebook (*how cool is that??*). <br>
The module tries to find a local version of JSME in <notebook_dir>/lib/ and when it fails to do so,
loads a web version of the editor. I use a central lib/ folder and create symlinks
to it in all notebook folders where I want to use these libraries

...plus many others.

# Module pipeline

A Pipelining Workflow using Python Generators, mainly for RDKit and large compound sets.
The use of generators allows working with arbitrarily large data sets, the memory usage at any given time is low.

Example use:

    >>> from rdkit_ipynb_tools import pipeline as p
    >>> s = Summary()
    >>> rd = start_csv_reader(test_data_b64.csv.gz", summary=s)
    >>> b64 = pipe_mol_from_b64(rd, summary=s)
    >>> filt = pipe_mol_filter(b64, "[H]c2c([H])c1ncoc1c([H])c2C(N)=O", summary=s)
    >>> stop_sdf_writer(filt, "test.sdf", summary=s)

or, using the pipe function:

    >>> s = Summary()
    >>> rd = start_sdf_reader("test.sdf", summary=s)
    >>> pipe(rd,
    >>>      pipe_keep_largest_fragment,
    >>>      (pipe_neutralize_mol, {"summary": s}),
    >>>      (pipe_keep_props, ["Ordernumber", "NP_Score"]),
    >>>      (stop_csv_writer, "test.csv", {"summary": s})
    >>>     )

The progress of the pipeline is displayed as a HTML table in the Notebook and can also be followed in a separate terminal with: `watch -n 2 cat pipeline.log`.

## Currently Available Pipeline Components:
| Starting                   | Running                    | Stopping
|----------------------------|----------------------------|---------------------------|
| start_cache_reader         | pipe_calc_props            | stop_cache_writer         |
| start_csv_reader           | pipe_custom_filter         | stop_count_records        |
| start_mol_csv_reader       | pipe_custom_man            | stop_csv_writer           |
| start_sdf_reader           | pipe_do_nothing            | stop_df_from_stream       |
| start_stream_from_dict     | pipe_has_prop_filter       | stop_dict_from_stream     |
| start_stream_from_mol_list | pipe_id_filter             | stop_mol_list_from_stream |
|                            | pipe_inspect_stream        | stop_sdf_writer           |
|                            | pipe_join_data_from_file   |                           |
|                            | pipe_keep_largest_fragment |                           |
|                            | pipe_keep_props            |                           |
|                            | pipe_merge_data            |                           |
|                            | pipe_mol_filter            |                           |
|                            | pipe_mol_from_b64          |                           |
|                            | pipe_mol_from_smiles       |                           |
|                            | pipe_mol_to_b64            |                           |
|                            | pipe_mol_to_smiles         |                           |
|                            | pipe_neutralize_mol        |                           |
|                            | pipe_remove_props          |                           |
|                            | pipe_rename_prop           |                           |
|                            | pipe_sim_filter            |                           |
|                            | pipe_sleep                 |                           |


Limitation: unlike in other pipelining tools, because of the nature of Python generators, the pipeline can not be branched.

# Other Modules
## Clustering
Fully usable, documentation needs to be written.
Please refer to the docstrings until then.

## Scaffolds
New, WIP, **not** usable. Has been moved to the scaffolds branch.

# Tutorial
Much of the functionality is shown in the [tools tutorial notebook](tutorial/tutorial_tools.ipynb).
SAR functionality is shown in the [SAR tutorial notebook](tutorial/tutorial_sar.ipynb). The SAR module is new and Work in Progress.

# Documentation
The module documentation can be built with sphinx using the `make_doc.sh` script

# Installation
## Requirements
The recommended way to use this project is via conda.

1. Python 3
1. [RDKit](http://www.rdkit.org/)
1. Jupyter Notebook
1. ipywidgets

## Highly recommended
1. cairo (via conda or pip) and cairocffi (only via pip)
to get decent-looking structures
1. [Bokeh](http://bokeh.pydata.org/en/latest/) for high-quality data plots
with structure tooltips

After installing the requirements,
clone this repo, then the rdkit_ipynb_tools can be used by including
the project's base directory (`rdkit_ipynb_tools`)
in Python's import path (I actually prefer this to using setuptools,
because a simple `git pull` will get you the newest version). <br>
This can be achieved by one of the following: <br>
* If you use conda (recommended), use [conda develop](http://conda.pydata.org/docs/commands/build/conda-develop.html).
This works similar to the next option.
* Put a file with the extension `.pth`, e.g. `my_packages.pth`,
into one of the `site-packages` directories of your Python installation
and put the path to the base directory of this project
(`rdkit_ipynb_tools`) into it. <br>
(I have the path to a dedicated folder on my machine included in such a `.pth`
file and link all my development projects to that folder.
This way, I only need to create the `.pth` file once.)

# Tips & Tricks
## Pipelines, Structures and Performance
Processing data from 200k compounds takes 10-15 sec on my notebook.

Substructure searches take longer.

For performance reasons, I use b64encode and pickle strings of mol objects to store the molecule structures in text format<br>
(see also Greg's blog post for [faster structure generation](http://rdkit.blogspot.de/2016/09/avoiding-unnecessary-work-and.html)):

```python
b64encode(pickle.dumps(mol)).decode()
```
For me, that has proven to be the fastest method when dealing with flat text files and is also the reason why there are `pipe_mol_to_b64` and `pipe_mol_from_b64` components in the `pipeline` module.

## Working Offline
* When you use a local copy of the Javascript Molecule Editor as described above
and use Bokeh for plotting, you can work completely offline in your Notebook.

# Roadmap
* make pipelines more user-friendly
* complete the scaffolds module
* add functionality as needed / requested

(probably not in this order)
