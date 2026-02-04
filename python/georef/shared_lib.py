import ctypes
import os
import logging

# This file is either
# - ${BUILD}/python/georef/shared_lib.py
#   with libgeoref.so in ${BUILD}/src/lib/libgeoref.so
# - ${PACKAGE}/lib/python/georef/shared_lib.py
#   with libgeoref in ${PACKAGE}/lib/libgeoref.so

this_dir = os.path.dirname(__file__)
georef_build_lib = os.path.normpath(os.path.join(this_dir, "..", "..", "src", "lib", "libgeoref.so"))
georef_package_lib = os.path.normpath(os.path.join(this_dir, "..", "..", "libgeoref.so"))
if os.path.isfile(georef_build_lib):
    libgeoref = ctypes.cdll.LoadLibrary(georef_build_lib)
elif os.path.isfile(georef_package_lib):
    libgeoref = ctypes.cdll.LoadLibrary(georef_package_lib)
else:
    logging.warn("libgeoref.so not found relative to python package")
    logging.info("attempting to load via LD_LIBRARY_PATH")
    libgeoref = ctypes.cdll.LoadLibrary("libgeoref.so")

