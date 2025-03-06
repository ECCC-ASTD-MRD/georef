"""Structure definitions for the georef package."""

import ctypes

class GeoOptions(ctypes.Structure):
    _fields_ = [
        ("Interp", ctypes.c_int32),         # Interpolation degree
        ("Extrap", ctypes.c_int32),         # Extrapolation method
        ("Combine", ctypes.c_int32),        # Aggregation type
        ("Transform", ctypes.c_int32),      # Apply transformation or stay within master referential
        ("CIndex", ctypes.c_int32),         # C Indexing (starts at 0)
        ("Symmetric", ctypes.c_int32),      # 
        ("Segment", ctypes.c_int32),        # How much segmentation (Conservatives/Geometric modes)
        ("Sampling", ctypes.c_int32),       # Sampling interval
        ("PolarCorrect", ctypes.c_int8),    # Apply polar corrections
        ("VectorMode", ctypes.c_int8),      # Process data as vector
        ("DistTreshold", ctypes.c_float),   # Distance treshold for point clouds
        ("NoData", ctypes.c_float),         # NoData Value (Default: NaN)
        ("Table", ctypes.c_void_p),         # Data table to check of values to check for
        ("lutDef", ctypes.c_void_p),        # Lookup table
        ("lutSize", ctypes.c_int32),        # Number of lookup elements
        ("lutDim", ctypes.c_int32),         # Dimension of the lookup elements
        ("Ancilliary", ctypes.c_void_p),    # Pre calculated field (ex: variance, average,...)
    ]

class GeoDef(ctypes.Structure):
    _fields_ = [
        ("ptr", ctypes.c_void_p),           # Pointer to C control structure
    ]

class GeoSet(ctypes.Structure):
    _fields_ = [
        ("ptr", ctypes.c_void_p),           # Pointer to C control structure
    ]

class _TGeoRef(ctypes.Structure):
    pass

class GeoRefError(Exception):
    """Exception raised for georef-specific errors."""
    pass 