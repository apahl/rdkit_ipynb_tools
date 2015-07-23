## RDKit IPython Tools
by Axel Pahl

The toolkit currently contains two functions to use with [RDKit](http://rdkit.org) in the IPython notebook:  
- *show_table*: Display a list of molecules in a table with molecule properties as columns.  
When an ID property is given, the table becomes interactive and compounds can be selected.  
When a highlight dictionary is passed (e.g. highlight={"activity": "< 50"}), the cell for the property with the name "activity" will be highlighted when its value is less than 50. To apply the highlight criteria to all cells except for the ID property, use the special key "\*all\*" (as in highlight={"\*all\*": "< 50"}).

- *jsme*: Display Peter Ertl's Javascript Molecule Editor to enter a molecule directly in the IPython notebook (*how cool is that??*)


A demonstration is worth a thousand words, so please have a look at the example notebook [rdkit_ipynb_tools.ipynb](http://nbviewer.ipython.org/github/apahl/rdkit_ipynb_tools/blob/master/rdkit_ipynb_tools_1.ipynb).

The toolkit was written with and for Python3, but should also work on Python2. If not, please let me know.