import ctypes
import os
import logging

# We are either
# - ${PACKAGE}/lib/python/georef/shared_lib.py
# - ${BUILD}/python/georef/shared_lib.py
this_dir = os.path.dirname(__file__)
# Either
# - ${BUILD}/python/georef
#   -> libgeoref.so is in ${this_dir}/../../lib/
# - ${PACKAGE}/lib/python/georef
#   -> libgeoref.so is in ${this_dir}/../../

georef_build_lib = os.path.normpath(os.path.join(this_dir, "..", "..", "src", "libgeoref.so"))
georef_package_lib = os.path.normpath(os.path.join(this_dir, "..", "..", "libgeoref.so"))
if os.path.isfile(georef_build_lib):
    print(f"Loading libgeoref from {georef_build_lib}")
    libgeoref = ctypes.cdll.LoadLibrary(georef_build_lib)
elif os.path.isfile(georef_package_lib):
    print(f"Loading libgeoref from {georef_package_lib}")
    libgeoref = ctypes.cdll.LoadLibrary(georef_package_lib)
else:
    print("libgeoref.so not found relative to python package")
    logging.warn("libgeoref.so not found relative to python package")
    logging.info("attempting to load via LD_LIBRARY_PATH")
    ctypes.cdll.LoadLibrary("libgeoref.so")

