import ctypes

import numpy as np

from rmn import fst24_file
from ._sharedlib import libgeoref

_georef_create = libgeoref.GeoRef_Create
_georef_create.argtypes = (ctypes.c_int32, # NI
                           ctypes.c_int32, # NJ
                           ctypes.c_char_p, # grid type
                           ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, # IG 1-4
                           ctypes.c_void_p, # fst_file pointer
                           )
_georef_create.restype = ctypes.c_void_p

_georef_axis_define = libgeoref.GeoRef_AxisDefine
_georef_axis_define.argtypes = (ctypes.c_void_p, # the georef
                                ctypes.c_void_p, # Longitudes
                                ctypes.c_void_p, # Latitudes
                                )

_georef_write_fst = libgeoref.GeoRef_WriteFST
_georef_write_fst.argtypes = (ctypes.c_void_p, # The georef
                              ctypes.c_char_p, # Name of the grid (optional)
                              ctypes.c_int32,  # IG1
                              ctypes.c_int32,  # IG2
                              ctypes.c_int32,  # IG3
                              ctypes.c_int32,  # IG4
                              ctypes.c_void_p, # output file
                              )
_georef_write_fst.restype = ctypes.c_int32

class TGeoRef(ctypes.Structure):
    def __init__(self, ni: int, nj: int, grid_type: str, ig1: int, ig2: int, ig3: int, ig4: int, file: fst24_file = None):
        self._file = file
        self._c_ref = _georef_create(ni, nj, grid_type.encode("utf-8"), ig1, ig2, ig3, ig4,
                                     None if file is None else file._c_ref)
        self._lons = None
        self._lats = None

    def define_axes(self, lons: np.ndarray, lats: np.ndarray):
        self._lons = lons
        self._lats = lats
        _georef_axis_define(self._c_ref, self._lons.ctypes.data, self._lats.ctypes.data)

    def write(self, name: str = "", output_file: fst24_file = None):
        the_file = output_file
        if the_file is None:
            the_file = self._file
        if the_file is None:
            raise ValueError("Need to specify a file where to write the grid")

        print(f'c ref = {the_file._c_ref}', flush=True)

        result = _georef_write_fst(self._c_ref, name.encode("utf-8"), 0, 0, 0, 0, the_file._c_ref)

        if result <= 0:
            raise ValueError(f"Error when calling GeoRef_Write (returned {result})")
