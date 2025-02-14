import ctypes
import numpy
import numpy.ctypeslib
# import rmn

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
        ("Weight", ctypes.POINTER(ctypes.c_float)),  # Weight array pointer
        ("WeightDim", ctypes.c_int32),              # Weight array dimensions
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

_free = libgeoref.GeoRef_Free
_free.argtypes = [ctypes.POINTER(_TGeoRef)]
_free.restype = None

_grid_value = libgeoref.GeoRef_GridValue
_grid_value.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.c_double, ctypes.c_double]
_grid_value.restype = ctypes.c_double

# Add error code constants
GEOREF_SUCCESS = 1
GEOREF_ERROR = 0

class GeoRefError(Exception):
    pass

class GeoRef:
    """Wrapper for the C GeoRef structure providing geographic reference functionality."""
    
    def __init__(self, ni: int, nj: int, grtyp: str, ig1: int, ig2: int, 
                 ig3: int, ig4: int, fst_file) -> None:
        """Initialize a new geographic reference.
        
        Args:
            ni: Number of points in x direction
            nj: Number of points in y direction
            grtyp: Grid type identifier
            ig1-ig4: Grid parameters
            fst_file: FST file reference
        """
        ptr = _georef_create(ni, nj, grtyp.encode('UTF-8'), ig1, ig2, ig3, ig4, fst_file._c_ref)
        if ptr.contents == 0:
            raise GeoRefError("Failure in C function GeoRef_Create")
        self._ptr = ptr

    def limits(self):
        lat0 = ctypes.c_double(0.0)
        lon0 = ctypes.c_double(0.0)
        lat1 = ctypes.c_double(0.0)
        lon1 = ctypes.c_double(0.0)
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
        return (lat0.value, lon0.value, lat1.value, lon1.value)

    def valid(self):
        return _valid(self._ptr)

    def __del__(self):
        if hasattr(self, '_ptr'):
            _free(self._ptr)

    def grid_value(self, x: float, y: float) -> float:
        """Get grid value at the specified coordinates.
        
        Args:
            x: X coordinate
            y: Y coordinate
        
        Returns:
            Grid value at the specified point
        """
        return _grid_value(self._ptr, x, y)

    def ll_value(self, lat: float, lon: float) -> float:
        """Get value at the specified lat/lon coordinates.
        
        Args:
            lat: Latitude
            lon: Longitude
        
        Returns:
            Value at the specified lat/lon
        """
        return _ll_value(self._ptr, lat, lon)

    def xy_to_ll(self, x: float, y: float) -> tuple[float, float]:
        """Convert grid coordinates to lat/lon.
        
        Args:
            x: X coordinate
            y: Y coordinate
        
        Returns:
            Tuple of (latitude, longitude)
        """
        lat = ctypes.c_double()
        lon = ctypes.c_double()
        result = _xy_to_ll(self._ptr, x, y, ctypes.byref(lat), ctypes.byref(lon))
        if result != GEOREF_SUCCESS:
            raise GeoRefError("Failed to convert XY to LL coordinates")
        return (lat.value, lon.value)

    def ll_to_xy(self, lat: float, lon: float) -> tuple[float, float]:
        """Convert lat/lon coordinates to grid coordinates.
        
        Args:
            lat: Latitude
            lon: Longitude
        
        Returns:
            Tuple of (x, y) grid coordinates
        """
        x = ctypes.c_double()
        y = ctypes.c_double()
        result = _ll_to_xy(self._ptr, lat, lon, ctypes.byref(x), ctypes.byref(y))
        if result != GEOREF_SUCCESS:
            raise GeoRefError("Failed to convert LL to XY coordinates")
        return (x.value, y.value)

def Interp(refto, reffrom, options, zout, zin):
    if not isinstance(zin, numpy.ndarray) or not isinstance(zout, numpy.ndarray):
        raise TypeError("Input and output must be numpy arrays")
    if zin.dtype != numpy.float32 or zout.dtype != numpy.float32:
        raise TypeError("Arrays must be float32")
    # Add shape validation here
    return _interp(refto._ptr, reffrom._ptr, ctypes.byref(options), zout, zin)
