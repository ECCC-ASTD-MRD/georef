"""Structure definitions for the georef package."""

from __future__ import annotations

import ctypes

import numpy as np

from .shared_lib import libgeoref


class GeoOptions(ctypes.Structure):
    # The pointer fields (table, lutDef, ancilliary) are exposed to Python as
    # numpy arrays through properties.  The numpy arrays are kept in __slots__
    # so the underlying buffers stay alive while the C pointer references them.
    __slots__ = ["_ancilliary_array", "_lutDef_array", "_lutDef_row_ptrs", "_table_array"]

    _fields_ = [
        ("Interp", ctypes.c_int32),  # Interpolation degree
        ("Extrap", ctypes.c_int32),  # Extrapolation method
        ("Combine", ctypes.c_int32),  # Aggregation type
        ("Transform", ctypes.c_int32),  # Apply transformation or stay within master referential
        ("Symmetric", ctypes.c_int32),  #
        ("Segment", ctypes.c_int32),  # How much segmentation (Conservatives/Geometric modes)
        ("Sampling", ctypes.c_int32),  # Sampling interval
        ("PolarCorrect", ctypes.c_int8),  # Apply polar corrections
        ("VectorMode", ctypes.c_int8),  # Process data as vector
        ("DistTreshold", ctypes.c_float),  # Distance treshold for point clouds
        ("NoData", ctypes.c_float),  # NoData Value (Default: NaN)
        ("_table", ctypes.c_void_p),  # Data table to check values against
        ("_lutDef", ctypes.c_void_p),  # Lookup table
        ("lutSize", ctypes.c_int32),  # Number of lookup elements
        ("lutDim", ctypes.c_int32),  # Dimension of the lookup elements
        ("_ancilliary", ctypes.c_void_p),  # Pre calculated field (ex: variance, average,...)
    ]

    def __new__(cls, *args, **kwargs):
        return _get_default_GeoOptions()

    def __init__(self, **kwargs):
        # Initialize the numpy-array holders.
        self._table_array = None
        self._lutDef_array = None
        self._ancilliary_array = None
        self._lutDef_row_ptrs = None

        # The pointer-backed properties must be applied after the plain fields
        # because they also update lutSize/lutDim.
        deferred = {name: kwargs.pop(name) for name in ("table", "lutDef", "ancilliary") if name in kwargs}

        for key, value in kwargs.items():
            if key.startswith("_"):
                raise ValueError(f"Attributes beginning with '_' should not be touched: '{key}'")
            setattr(self, key, value)

        for name, value in deferred.items():
            setattr(self, name, value)

    # --- pointer-backed numpy array properties -------------------------------

    @property
    def table(self) -> np.ndarray | None:
        """The data table as a 1D float64 numpy array (or None)."""
        return self._table_array

    @table.setter
    def table(self, value: np.ndarray):
        self._set_pointer_array("_table", "_table_array", value, ndim=1)

    @table.deleter
    def table(self):
        self._table = None
        self._table_array = None

    @property
    def ancilliary(self) -> np.ndarray | None:
        """The pre-calculated ancilliary field as a 1D float64 numpy array (or None)."""
        return self._ancilliary_array

    @ancilliary.setter
    def ancilliary(self, value: np.ndarray):
        self._set_pointer_array("_ancilliary", "_ancilliary_array", value, ndim=1)

    @ancilliary.deleter
    def ancilliary(self):
        self._ancilliary = None
        self._ancilliary_array = None

    @property
    def lutDef(self) -> np.ndarray | None:
        """The lookup table as a 2D float64 numpy array (or None)."""
        return self._lutDef_array

    @lutDef.setter
    def lutDef(self, value: np.ndarray):
        if not isinstance(value, np.ndarray):
            raise TypeError(f"Expected {np.ndarray.__name__}, got {type(value).__name__}")
        if value.ndim != 2:
            raise ValueError(f"lutDef must be a 2D array, got {value.ndim}D")
        value = np.ascontiguousarray(value, dtype=np.float64)
        self._lutDef_array = value
        # In C, lutDef is a double **: an array of row pointers.
        row_ptrs = (ctypes.POINTER(ctypes.c_double) * value.shape[0])()
        for i in range(value.shape[0]):
            row_ptrs[i] = ctypes.cast(value[i].ctypes.data, ctypes.POINTER(ctypes.c_double))
        self._lutDef_row_ptrs = row_ptrs
        self._lutDef = ctypes.addressof(row_ptrs)
        self.lutSize = value.shape[0]
        self.lutDim = value.shape[1]

    @lutDef.deleter
    def lutDef(self):
        self._lutDef = None
        self._lutDef_array = None
        self._lutDef_row_ptrs = None
        self.lutSize = 0
        self.lutDim = 0

    # --- helpers -------------------------------------------------------------

    def _set_pointer_array(self, field: str, array_attr: str, value: np.ndarray, ndim: int):
        """Store a numpy array and point the matching C pointer field at it."""
        if not isinstance(value, np.ndarray):
            raise TypeError(f"Expected {np.ndarray.__name__}, got {type(value).__name__}")
        if value.ndim != ndim:
            raise ValueError(f"Expected a {ndim}D array, got {value.ndim}D")
        value = np.ascontiguousarray(value, dtype=np.float64)
        setattr(self, array_attr, value)
        setattr(self, field, value.ctypes.data)


_get_default_GeoOptions = libgeoref.get_default_GeoOptions
_get_default_GeoOptions.argtypes = ()
_get_default_GeoOptions.restype = GeoOptions


class GeoRefError(Exception):
    """Exception raised for georef-specific errors."""
