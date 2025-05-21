# Build documentation

Run
```
make html
```
from this directory.

The result will be a directory `build/html`.  The simplest way to view it is to
create a link like `~/public_html/doc/georef/python -> $PWD/build/html` and to
visit `https://web.science.gc.ca/~<your-user>/doc/georef/python`.

## Environment requirements

To generate the documentation, the python package must be importable.  For this
to be the case, the directory containing the Python package must be in `sys.path`
which can be achieved by adding to the variable `PYTHONPATH` but since importing
`georef` loads `libgeoref.so`, we need the directory containing it to be in
`LD_LIBRARY_PATH`.

When loading `libgeoref.so`, we need to be able to find `librmn.so` and therefore
`libjson-c.so`, and `libApp.so`.

So to be able to run this, we need to load `librmn` probably via SSM, and
```
# Directory containing python package
export PYTHONPATH=${GEOREF_REPO}/python
# Directory containing `libgeoref.so`
export LD_LIBRARY_PATH=${GEOREF_BUILD}/src/${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}
```

# Modifying

## Adding files

Create a file `source/newfile.rst` and put anything in it.  Then to make Sphinx process
it, it needs to be mentioned in the `.. toctree::` section of `source/index.rst`

## Configuring

Configuration is done in `source/conf.py`
