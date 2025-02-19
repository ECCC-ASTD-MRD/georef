"""Structure definitions for the georef package."""

import ctypes

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
        ("Weight", ctypes.POINTER(ctypes.c_float)),  # Weight array pointer
        ("WeightDim", ctypes.c_int32),              # Weight array dimensions
    ]

class _TGeoRef(ctypes.Structure):
    pass

class GeoRefError(Exception):
    """Exception raised for georef-specific errors."""
    pass 