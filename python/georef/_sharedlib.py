import ctypes
import os
import logging

this_dir = os.path.dirname(__file__)
# __file__ == <build-dir>/python/georef/_sharedlib.py
# libgeoref.so == <build-dir>/libgeoref.so
rel_build_dir = os.path.normpath(os.path.join(this_dir, '..', '..', 'libgeoref.so'))
if os.path.isfile(rel_build_dir):
    libgeoref = ctypes.cdll.LoadLibrary(rel_build_dir)
else:
    # __file__ == <install-dir>/lib/python/georef/_sharedlib.py
    # libgeoref.so == <install-dir>/lib/libgeoref.so
    # (same relative positions as build dir but that is a coincidence)
    rel_install_dir = os.path.normpath(os.path.join(this_dir, '..', '..', 'libgeoref.so'))
    if os.path.isfile(rel_install_dir):
        libgeoref = ctypes.cdll.LoadLibrary(rel_install_dir)
    else:
        logging.warn("Loading georef though LD_LIBRARY_PATH.  Loaded shared library compatibility not garanteed")
        libgeoref = ctypes.cdll.LoadLibrary("libgeoref.so")
