## RDKit IPython Tools
by Axel Pahl

### Work in progress
*(9-Nov-2015)*  
In the past months, this set of tools based on the [RDKit](http.//www.rdkit.org) has evolved quite a bit.  


#### Module tools

A Mol_List class was introduced, which is a subclass of a Python list for holding lists of RDKit molecule objects and allows direct access to a lot of the RDKit functionality.  
It is meant to be used with the JuPyTer Notebook and includes a.o.:
* display of the Mol_List
    * as HTML table
    * as HTML grid
  (both display types include the option to select molecules by clicking)
* display of a summary including number of records and min, max, mean, median for numeric properties
* display of correlations between the Mol_List's properties  
  (using np.corrcoef, this allows getting a quick overview on which properties correlate with each other)
* methods for sorting, searching (by property or substructure) and filtering the Mol_List
* methods for renaming, reordering and calculating properties
* direct plotting of properties as publication-grade [Highcharts](http://www.highcharts.com/) plots with **structure tooltips** (!).
    * the highcharts plotting functionality resides in its own module and can also be used for plotting Pandas dataframes and Python dicts. Scatter and column plots are currently supported.


##### Other functions in the tools module:
- *jsme*: Display Peter Ertl's Javascript Molecule Editor to enter a molecule directly in the IPython notebook (*how cool is that??*)

plus many others

### ToDo
This README is meant as a teaser. I urgently need to prepare a proper example notebook to show all the functionality, this is high on my ToDo list.
