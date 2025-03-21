"""Structure definitions for the georef package."""

import ctypes
import numpy as np

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
        ("_lutDef", ctypes.c_void_p),        # Lookup table
        ("lutSize", ctypes.c_int32),        # Number of lookup elements
        ("lutDim", ctypes.c_int32),         # Dimension of the lookup elements
        ("Ancilliary", ctypes.c_void_p),    # Pre calculated field (ex: variance, average,...)
    ]
    # create getter and setter for c_void_p fields and change their names to lutDef, table, ancilliary 
    # make sure lutSize Dim not visible to user
    @property
    def lutDef(self):
        return self._lutDef_ndarray
    @lutDef.setter
    def lutDef(self, value: np.ndarray):
        self._lutDef_ndarray = value
        self._lutDef = value.ctypes.data
        self.lutSize = value.size # Not sure about this one
        self.lutDim = value.shape[1] # Not sure about this one
    
    @lutDef.deleter
    def lutDef(self):
        self._lutDef_ndarray = None
        self._lutDef = None
        self.lutSize = 0
        self.lutDim = 0


    @property
    def table(self):
        return self._table
    @table.setter
    def table(self, value):
        self._table = value
    @property
    def ancilliary(self):
        return self._ancilliary
    @ancilliary.setter
    def ancilliary(self, value):
        self._ancilliary = value
   


class GeoRefError(Exception):
    """Exception raised for georef-specific errors."""
    pass 
