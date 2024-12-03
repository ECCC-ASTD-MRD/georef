import ctypes
import numpy
import numpy.ctypeslib
import rmn
from . import enums

# import the `libgeoref` from __init__.py
from .shared_lib import libgeoref

class GeoOptions(ctypes.Structure):
    _fields_ = [
        ("Interp", ctypes.c_int),           # Interpolation degree
        ("Extrap", ctypes.c_int),           # Extrapolation method
        ("InterpVector", ctypes.c_int),     # Vector interpolation method
        ("Combine", ctypes.c_int),          # Aggregation type
        ("ExtrapValue", ctypes.c_double),   # Value to use for extrapolation in ER_VALUE mode
        ("Transform", ctypes.c_int32),      # Apply transformation or stay within master referential
        ("CIndex", ctypes.c_int32),         # C Indexing (starts st 0)
        ("Symmetric", ctypes.c_int32),      #
        ("WeightNum", ctypes.c_int32),      #
        # TODO Add fields when merging with origin/dev
        ("Segment", ctypes.c_int32),        # How much segmentation (Conservatives/Geometric modes)
        ("Sampling", ctypes.c_int32),       # Sampling interval
        ("PolarCorrect", ctypes.c_char),    # Apply polar corrections
        ("VectorMode", ctypes.c_char),      # Process data as vector
        ("DistTreshold", ctypes.c_float),   # Distance treshold for point clouds
        ("LonRef", ctypes.c_float),         # Longitude referential (-180.0,0.0)
        ("NoData", ctypes.c_float),         # NoData Value (Default: NaN)
    ]

class _TGeoRef(ctypes.Structure):
    pass

_new = libgeoref.GeoRef_New
_new.argtypes = []
_new.restype = ctypes.POINTER(_TGeoRef)

_valid = libgeoref.GeoRef_Valid
_valid.argtypes = [ctypes.POINTER(_TGeoRef)]
_valid.restype = ctypes.c_int

_georef_create = libgeoref.GeoRef_Create
_georef_create.argtypes = (
    ctypes.c_int, ctypes.c_int, ctypes.c_char_p,
    ctypes.c_int, ctypes.c_int,
    ctypes.c_int, ctypes.c_int, ctypes.c_int
    )
_georef_create.restype = ctypes.POINTER(_TGeoRef)

_georef_limits = libgeoref.GeoRef_Limits
_georef_limits.argtypes = (ctypes.POINTER(_TGeoRef),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double)
)
_georef_limits.restype = ctypes.c_int32

_interp = libgeoref.GeoRef_Interp
_interp.argtypes = (
    ctypes.POINTER(_TGeoRef),
    ctypes.POINTER(_TGeoRef),
    ctypes.POINTER(GeoOptions),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32)
)
_interp.restype = ctypes.c_int


class GeoRefError(Exception):
    pass

class GeoRef:
    def __init__(self, ni, nj, grtyp, ig1, ig2, ig3, ig4, fst_file):
        ptr = _georef_create(ni, nj, grtyp.encode('UTF-8'), ig1, ig2, ig3, ig4, fst_file._c_ref)
        if ptr.contents == 0:
            raise GeoRefError("Failure in C function GeoRef_Create")
        self._ptr = ptr

    def limits(self):
        lat0 = 0.0
        lon0 = 0.0
        lon1 = 0.0
        lon1 = 0.0
        result = _georef_limits(
            self._ptr,
            ctypes.byref(lat0),
            ctypes.byref(lon0),
            ctypes.byref(lat1),
            ctypes.byref(lon1)
        )
        if result != 0:
            raise GeoRefError(f"Failure in C function GeoRef_Limits :{result}")
        # return {"lat0": lat0, ...}
        return (lat0, lon0, lat1, lon1)

    def valid(self):
        return _valid(self._ptr)

def Interp(refto, reffrom, options, zout, zin):
    res = _interp(refto._ptr, reffrom._ptr, ctypes.byref(options), zout, zin)
    if res != 1:  # TODO: Use a named constant like GEOREF_SUCCESS
        raise GeoRefError("Failure in C function GeoRef_Interp")
