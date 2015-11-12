The documentation for this module can be build using sphinx by running `make html` in this dir. The module dir rdkit_ipynb_tools has to be in the Python import path  
This can be achieved E.g. by one of the following methods
1. Put the name of the package dir in a custom .pth file in Python's site-packages or dist-packages  
    (that's how I do it as long as it is not yet a real package)
2. Add the name of the dir to the PYTHONPATH variable
3. Copy it directly to site-packages or dist-packages
