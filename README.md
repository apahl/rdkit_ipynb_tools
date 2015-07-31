## RDKit IPython Tools
by Axel Pahl

#### module tools
The toolkit currently contains two functions to use with [RDKit](http://rdkit.org) in the IPython notebook:
- *show_table*: Display a list of molecules in a table with molecule properties as columns.
When an ID property is given, the table becomes interactive and compounds can be selected.
When a highlight dictionary is passed (e.g. highlight={"activity": "< 50"}), the cell for the property with the name "activity" will be highlighted when its value is less than 50. To apply the highlight criteria to all cells except for the ID property, use the special key "\*all\*" (as in highlight={"\*all\*": "< 50"}).

- *jsme*: Display Peter Ertl's Javascript Molecule Editor to enter a molecule directly in the IPython notebook (*how cool is that??*)

#### module hc_tools
Display [Highcharts](http://www.highcharts.com) plots in the Notebook
This is still work in progress. Currently the tools included here produce scatter plots from dictionaries. The plan is to combine this functionality with the mol lists from RDKit.


A demonstration is worth a thousand words, so please have a look at the example notebooks
* [rdkit_ipynb_tools_1.ipynb](http://nbviewer.ipython.org/github/apahl/rdkit_ipynb_tools/blob/master/rdkit_ipynb_tools_1.ipynb).
* [highcharts.ipynb](http://nbviewer.ipython.org/github/apahl/rdkit_ipynb_tools/blob/master/highcharts.ipynb)

The toolkit was written with and for Python3.