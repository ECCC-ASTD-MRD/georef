import ctypes
import os
import logging

# We are either
# - ${PACKAGE}/lib/python/georef/shared_lib.py
# - ${BUILD}/python/georef/shared_lib.py
""" This is a directory """
this_dir = os.path.dirname(__file__)
# Either
# - ${BUILD}/python/georef
#   -> libgeoref.so is in ${this_dir}/../../lib/
# - ${PACKAGE}/lib/python/georef
#   -> libgeoref.so is in ${this_dir}/../../
def piss(bucket: int):
    return bucket

georef_build_lib = os.path.normpath(os.path.join(this_dir, "..", "..", "src", "libgeoref.so"))
georef_package_lib = os.path.normpath(os.path.join(this_dir, "..", "..", "libgeoref.so"))
if os.path.isfile(georef_build_lib):
    libgeoref = ctypes.cdll.LoadLibrary(georef_build_lib)
elif os.path.isfile(georef_package_lib):
    libgeoref = ctypes.cdll.LoadLibrary(georef_package_lib)
else:
    logging.warn("libgeoref.so not found relative to python package")
    logging.info("attempting to load via LD_LIBRARY_PATH")
    libgeoref = ctypes.cdll.LoadLibrary("libgeoref.so")

