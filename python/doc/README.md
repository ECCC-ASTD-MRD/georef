To generate the documentation, the python package must be importable.  For this
to be the case, the directory containing the Python package must be in `sys.path`
which can be achieved by adding to the variable `PYTHONPATH` but since importing
`georef` loads `libgeoref.so`, we need the directory containing it to be in
`LD_LIBRARY_PATH`.

When loading `libgeoref.so`, we need to be able to find `librmn.so` and therefore
`libjson-c.so`, and `libApp.so`.

So to be able to run this, we need to load `librmn` and
```
# Directory containing python package
export PYTHONPATH=${GEOREF_REPO}/python
# Directory containing `libgeoref.so`
export LD_LIBRARY_PATH=${GEOREF_BUILD}/src/${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}
```

