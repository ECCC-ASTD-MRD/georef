import ctypes
import numpy
import numpy.ctypeslib
from ._georef_c_bindings import (
    _new,
    _valid,
    _georef_create,
    _georef_limits,
    _interp,
    _free,
    _grid_value,
    _ll_value,
    _xy_to_ll,
    _ll2xy,
    _xy2ll,
    _copy,
    _hardcopy,
    _equal,
    _within,
    _withinrange,
    _intersect,
    _boundingbox,
    _xydistance,
    _lldistance,
    _interpuv,
    _interpwd,
    _ud2wd,
    _wd2uv,
    _uv2uv,
    _llwdval,
    _lluvval,
    _llval,
    _xywdval,
    _xyuvval,
    _xyval,
    _getll,
)

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

    def limits(self) -> tuple[float, float, float, float]:
        """Get the geographical limits of the GeoRef.
        
        From: georef_limits_f in GeoRef_mod.F90
        C signature: int32_t georef_limits(GeoRef* ref,
                                         double* lat0, double* lon0,
                                         double* lat1, double* lon1)
        
        Returns:
            tuple: (lat0, lon0, lat1, lon1) geographical limits
        """
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

    def ll2xy(self, lat: float, lon: float) -> tuple[float, float]:
        """Convert lat/lon coordinates to grid coordinates.
        
        From: georef_ll2xy_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_ll2xy(GeoRef* ref,
                                        double lat, double lon,
                                        double* x, double* y)
        
        Args:
            lat: Latitude
            lon: Longitude
            
        Returns:
            tuple: (x, y) grid coordinates
        """
        x = ctypes.c_double()
        y = ctypes.c_double()
        result = _ll2xy(self._ptr, lat, lon, ctypes.byref(x), ctypes.byref(y))
        if result != GEOREF_SUCCESS:
            raise GeoRefError("Failed to convert LL to XY coordinates")
        return (x.value, y.value)

    def xy2ll(self, x: float, y: float) -> tuple[float, float]:
        """Convert grid coordinates to lat/lon.
        
        From: georef_xy2ll_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_xy2ll(GeoRef* ref,
                                        double x, double y,
                                        double* lat, double* lon)
        
        Args:
            x: X coordinate
            y: Y coordinate
            
        Returns:
            tuple: (latitude, longitude)
        """
        lat = ctypes.c_double()
        lon = ctypes.c_double()
        result = _xy2ll(self._ptr, x, y, ctypes.byref(lat), ctypes.byref(lon))
        if result != GEOREF_SUCCESS:
            raise GeoRefError("Failed to convert XY to LL coordinates")
        return (lat.value, lon.value)

    def copy(self, hard: bool = False) -> 'GeoRef':
        """Create a copy of the GeoRef object.
        
        From: georef_copy_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: GeoRef* georef_copy(GeoRef* ref)
                    GeoRef* georef_hardcopy(GeoRef* ref)
        
        Args:
            hard: If True, creates a hard copy
            
        Returns:
            New GeoRef instance
        """
        ptr = _hardcopy(self._ptr) if hard else _copy(self._ptr)
        if not ptr:
            raise GeoRefError("Failed to copy GeoRef")
        result = GeoRef.__new__(GeoRef)  # Create new object without __init__
        result._ptr = ptr
        return result

    def equal(self, ref: 'GeoRef') -> bool:
        """Check if two GeoRef objects are equal.
        
        From: georef_equal_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_equal(GeoRef* ref1, GeoRef* ref2)
        
        Args:
            ref: GeoRef object to compare with
            
        Returns:
            bool: True if equal
        """
        result = _equal(self._ptr, ref._ptr)
        return result == GEOREF_SUCCESS

    def within(self, ref: 'GeoRef') -> bool:
        """Check if this GeoRef is within another GeoRef.
        
        From: georef_within_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_within(GeoRef* ref1, GeoRef* ref2)
        
        Args:
            ref: GeoRef object to check against
            
        Returns:
            bool: True if within
        """
        result = _within(self._ptr, ref._ptr)
        return result == GEOREF_SUCCESS

    def withinrange(self, lat0: float, lon0: float, lat1: float, lon1: float, inside: bool = False) -> bool:
        """Check if GeoRef is within specified lat/lon range.
        
        From: georef_withinrange_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_withinrange(GeoRef* ref, double lat0, double lon0, 
                                              double lat1, double lon1, int32_t in)
        
        Args:
            lat0: Starting latitude
            lon0: Starting longitude
            lat1: Ending latitude
            lon1: Ending longitude
            inside: If True, checks if completely inside
            
        Returns:
            bool: True if within range
        """
        result = _withinrange(self._ptr, lat0, lon0, lat1, lon1, int(inside))
        return result == GEOREF_SUCCESS

    def intersect(self, ref: 'GeoRef', boundary: bool = False) -> tuple[int, int, int, int]:
        """Find intersection between two GeoRefs.
        
        From: georef_intersect_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_intersect(GeoRef* ref1, GeoRef* ref2, 
                                            int32_t* x0, int32_t* y0, 
                                            int32_t* x1, int32_t* y1, int32_t bd)
        
        Args:
            ref: Other GeoRef object
            boundary: If True, includes boundary
            
        Returns:
            tuple: (x0, y0, x1, y1) intersection coordinates
        """
        x0 = ctypes.c_int32()
        y0 = ctypes.c_int32()
        x1 = ctypes.c_int32()
        y1 = ctypes.c_int32()
        result = _intersect(self._ptr, ref._ptr, 
                           ctypes.byref(x0), ctypes.byref(y0),
                           ctypes.byref(x1), ctypes.byref(y1),
                           int(boundary))
        if result != GEOREF_SUCCESS:
            raise GeoRefError("Failed to find intersection")
        return (x0.value, y0.value, x1.value, y1.value)

    def boundingbox(self, lat0: float, lon0: float, lat1: float, lon1: float) -> tuple[float, float, float, float]:
        """Get bounding box coordinates.
        
        From: georef_boundingbox_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_boundingbox(GeoRef* ref, 
                                              double lat0, double lon0,
                                              double lat1, double lon1,
                                              double* i0, double* j0,
                                              double* i1, double* j1)
        
        Args:
            lat0: Starting latitude
            lon0: Starting longitude
            lat1: Ending latitude
            lon1: Ending longitude
            
        Returns:
            tuple: (i0, j0, i1, j1) bounding box coordinates
        """
        i0 = ctypes.c_double()
        j0 = ctypes.c_double()
        i1 = ctypes.c_double()
        j1 = ctypes.c_double()
        result = _boundingbox(self._ptr, lat0, lon0, lat1, lon1,
                             ctypes.byref(i0), ctypes.byref(j0),
                             ctypes.byref(i1), ctypes.byref(j1))
        if result != GEOREF_SUCCESS:
            raise GeoRefError("Failed to get bounding box")
        return (i0.value, j0.value, i1.value, j1.value)

    def xydistance(self, x0: float, y0: float, x1: float, y1: float) -> float:
        """Calculate distance between two points in grid coordinates.
        
        From: georef_xydistance_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: double georef_xydistance(GeoRef* ref, 
                                            double x0, double y0,
                                            double x1, double y1)
        
        Args:
            x0: First x coordinate
            y0: First y coordinate
            x1: Second x coordinate
            y1: Second y coordinate
            
        Returns:
            float: Distance between points
        """
        return _xydistance(self._ptr, x0, y0, x1, y1)

    def lldistance(self, lat0: float, lon0: float, lat1: float, lon1: float) -> float:
        """Calculate distance between two points in lat/lon coordinates.
        
        From: georef_lldistance_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: double georef_lldistance(GeoRef* ref,
                                            double lat0, double lon0,
                                            double lat1, double lon1)
        
        Args:
            lat0: First latitude
            lon0: First longitude
            lat1: Second latitude
            lon1: Second longitude
            
        Returns:
            float: Distance between points
        """
        return _lldistance(self._ptr, lat0, lon0, lat1, lon1)

    def write(self, filename: str) -> bool:
        """Write GeoRef to file.
        
        From: georef_write_f in GeoRef_mod.F90
        C signature: int32_t georef_write(GeoRef* ref, const char* filename)
        
        Args:
            filename: Path to output file
            
        Returns:
            bool: True if write successful
        """
        pass

    def fromrecord(self, record: 'FSTRecord') -> bool:
        """Initialize GeoRef from FST record.
        
        From: georef_fromrecord_f in GeoRef_mod.F90
        C signature: int32_t georef_fromrecord(GeoRef* ref, FSTRecord* rec)
        
        Args:
            record: FST record to initialize from
            
        Returns:
            bool: True if initialization successful
        """
        pass

    def interp(self, reffrom: 'GeoRef', zout: numpy.ndarray, zin: numpy.ndarray, opt: GeoOptions = None) -> int:
        """Interpolate values from one grid to another.
        
        From: georef_interp_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_interp(GeoRef* refto, GeoRef* reffrom, 
                                         GeoOptions* opt,
                                         float* zout, float* zin)
        
        Args:
            reffrom: Source GeoRef
            zout: Output array (float32)
            zin: Input array (float32)
            opt: Optional interpolation options
            
        Returns:
            int: Status code (GEOREF_SUCCESS or GEOREF_ERROR)
        """
        if not isinstance(zin, numpy.ndarray) or not isinstance(zout, numpy.ndarray):
            raise TypeError("Input and output must be numpy arrays")
        if zin.dtype != numpy.float32 or zout.dtype != numpy.float32:
            raise TypeError("Arrays must be float32")
        return _interp(self._ptr, reffrom._ptr, ctypes.byref(opt) if opt else None, zout, zin)

    def interpuv(self, reffrom: 'GeoRef', uuout: numpy.ndarray, vvout: numpy.ndarray, 
                 uuin: numpy.ndarray, vvin: numpy.ndarray, opt: GeoOptions = None) -> int:
        """Interpolate UV vector components.
        
        From: georef_interpuv_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_interpuv(GeoRef* refto, GeoRef* reffrom,
                                           GeoOptions* opt,
                                           float* uuout, float* vvout,
                                           float* uuin, float* vvin)
        
        Args:
            reffrom: Source GeoRef
            uuout: Output U component array (float32)
            vvout: Output V component array (float32)
            uuin: Input U component array (float32)
            vvin: Input V component array (float32)
            opt: Optional interpolation options
            
        Returns:
            int: Status code (GEOREF_SUCCESS or GEOREF_ERROR)
        """
        for arr in [uuin, vvin, uuout, vvout]:
            if not isinstance(arr, numpy.ndarray):
                raise TypeError("Arrays must be numpy arrays")
            if arr.dtype != numpy.float32:
                raise TypeError("Arrays must be float32")
        return _interpuv(self._ptr, reffrom._ptr, ctypes.byref(opt) if opt else None, 
                         uuout, vvout, uuin, vvin)

    def interpwd(self, reffrom: 'GeoRef', uuout: numpy.ndarray, vvout: numpy.ndarray,
                 uuin: numpy.ndarray, vvin: numpy.ndarray, opt: GeoOptions = None) -> int:
        """Interpolate wind direction components.
        
        From: georef_interpwd_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_interpwd(GeoRef* refto, GeoRef* reffrom,
                                           GeoOptions* opt,
                                           float* uuout, float* vvout,
                                           float* uuin, float* vvin)
        
        Args:
            reffrom: Source GeoRef
            uuout: Output U component array (float32)
            vvout: Output V component array (float32)
            uuin: Input U component array (float32)
            vvin: Input V component array (float32)
            opt: Optional interpolation options
            
        Returns:
            int: Status code (GEOREF_SUCCESS or GEOREF_ERROR)
        """
        for arr in [uuin, vvin, uuout, vvout]:
            if not isinstance(arr, numpy.ndarray):
                raise TypeError("Arrays must be numpy arrays")
            if arr.dtype != numpy.float32:
                raise TypeError("Arrays must be float32")
        return _interpwd(self._ptr, reffrom._ptr, ctypes.byref(opt) if opt else None, 
                         uuout, vvout, uuin, vvin)

    def ud2wd(self, uuout: numpy.ndarray, vvout: numpy.ndarray, 
              spdin: numpy.ndarray, dirin: numpy.ndarray, 
              lat: numpy.ndarray, lon: numpy.ndarray, npts: int) -> int:
        """Convert UV components to wind direction.
        
        From: georef_uv2wd_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_uv2wd(GeoRef* ref,
                                         float* uuout, float* vvout,
                                         float* spdin, float* dirin,
                                         double* lat, double* lon,
                                         int32_t npts)
        
        Args:
            uuout: Output U component array (float32)
            vvout: Output V component array (float32)
            spdin: Input speed array (float32)
            dirin: Input direction array (float32)
            lat: Latitude coordinates array (float64)
            lon: Longitude coordinates array (float64)
            npts: Number of points
            
        Returns:
            int: Status code
        """
        for arr in [uuout, vvout, spdin, dirin]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float32:
                raise TypeError("UV/Speed/Direction arrays must be float32 numpy arrays")
        for arr in [lat, lon]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float64:
                raise TypeError("Lat/Lon arrays must be float64 numpy arrays")
        return _ud2wd(self._ptr, uuout, vvout, spdin, dirin, lat, lon, npts)

    def wd2uv(self, spdout: numpy.ndarray, dirout: numpy.ndarray,
              uuin: numpy.ndarray, vvin: numpy.ndarray,
              lat: numpy.ndarray, lon: numpy.ndarray, npts: int) -> int:
        """Convert wind direction to UV components.
        
        From: georef_wd2uv_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_wd2uv(GeoRef* ref,
                                         float* spdout, float* dirout,
                                         float* uuin, float* vvin,
                                         double* lat, double* lon,
                                         int32_t npts)
        
        Args:
            spdout: Output speed array (float32)
            dirout: Output direction array (float32)
            uuin: Input U component array (float32)
            vvin: Input V component array (float32)
            lat: Latitude coordinates array (float64)
            lon: Longitude coordinates array (float64)
            npts: Number of points
            
        Returns:
            int: Status code
        """
        for arr in [spdout, dirout, uuin, vvin]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float32:
                raise TypeError("UV/Speed/Direction arrays must be float32 numpy arrays")
        for arr in [lat, lon]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float64:
                raise TypeError("Lat/Lon arrays must be float64 numpy arrays")
        return _wd2uv(self._ptr, spdout, dirout, uuin, vvin, lat, lon, npts)

    def uv2uv(self, uuout: numpy.ndarray, vvout: numpy.ndarray,
              uuin: numpy.ndarray, vvin: numpy.ndarray,
              lat: numpy.ndarray, lon: numpy.ndarray, npts: int) -> int:
        """Convert UV components between coordinate systems.
        
        From: georef_uv2uv_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_uv2uv(GeoRef* ref,
                                         float* uuout, float* vvout,
                                         float* uuin, float* vvin,
                                         double* lat, double* lon,
                                         int32_t npts)
        
        Args:
            uuout: Output U component array (float32)
            vvout: Output V component array (float32)
            uuin: Input U component array (float32)
            vvin: Input V component array (float32)
            lat: Latitude coordinates array (float64)
            lon: Longitude coordinates array (float64)
            npts: Number of points
            
        Returns:
            int: Status code
        """
        for arr in [uuout, vvout, uuin, vvin]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float32:
                raise TypeError("UV arrays must be float32 numpy arrays")
        for arr in [lat, lon]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float64:
                raise TypeError("Lat/Lon arrays must be float64 numpy arrays")
        return _uv2uv(self._ptr, uuout, vvout, uuin, vvin, lat, lon, npts)

    def llwdval(self, uuout: numpy.ndarray, vvout: numpy.ndarray,
                uuin: numpy.ndarray, vvin: numpy.ndarray,
                lat: numpy.ndarray, lon: numpy.ndarray, npts: int, opt: GeoOptions = None) -> int:
        """Get wind direction values at lat/lon coordinates.
        
        From: georef_llwdval_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_llwdval(GeoRef* ref, GeoOptions* opt,
                                          float* uuout, float* vvout,
                                          float* uuin, float* vvin,
                                          double* lat, double* lon,
                                          int32_t npts)
        
        Args:
            uuout: Output U component array (float32)
            vvout: Output V component array (float32)
            uuin: Input U component array (float32)
            vvin: Input V component array (float32)
            lat: Latitude coordinates array (float64)
            lon: Longitude coordinates array (float64)
            npts: Number of points
            opt: Optional interpolation options
            
        Returns:
            int: Status code
        """
        for arr in [uuin, vvin, uuout, vvout]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float32:
                raise TypeError("UV arrays must be float32 numpy arrays")
        for arr in [lat, lon]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float64:
                raise TypeError("Lat/Lon arrays must be float64 numpy arrays")
        return _llwdval(self._ptr, ctypes.byref(opt) if opt else None,
                        uuout, vvout, uuin, vvin, lat, lon, npts)

    def lluvval(self, uuout: numpy.ndarray, vvout: numpy.ndarray,
                uuin: numpy.ndarray, vvin: numpy.ndarray,
                lat: numpy.ndarray, lon: numpy.ndarray, npts: int, opt: GeoOptions = None) -> int:
        """Get UV vector values at lat/lon coordinates.
        
        From: georef_lluvval_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_lluvval(GeoRef* ref, GeoOptions* opt,
                                          float* uuout, float* vvout,
                                          float* uuin, float* vvin,
                                          double* lat, double* lon,
                                          int32_t npts)
        
        Args:
            uuout: Output U component array (float32)
            vvout: Output V component array (float32)
            uuin: Input U component array (float32)
            vvin: Input V component array (float32)
            lat: Latitude coordinates array (float64)
            lon: Longitude coordinates array (float64)
            npts: Number of points
            opt: Optional interpolation options
            
        Returns:
            int: Status code
        """
        for arr in [uuin, vvin, uuout, vvout]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float32:
                raise TypeError("UV arrays must be float32 numpy arrays")
        for arr in [lat, lon]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float64:
                raise TypeError("Lat/Lon arrays must be float64 numpy arrays")
        return _lluvval(self._ptr, ctypes.byref(opt) if opt else None,
                        uuout, vvout, uuin, vvin, lat, lon, npts)

    def llval(self, zout: numpy.ndarray, zin: numpy.ndarray,
              lat: numpy.ndarray, lon: numpy.ndarray, npts: int, opt: GeoOptions = None) -> int:
        """Get values at lat/lon coordinates.
        
        From: georef_llval_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_llval(GeoRef* ref, GeoOptions* opt,
                                        float* zout, float* zin,
                                        double* lat, double* lon,
                                        int32_t npts)
        
        Args:
            zout: Output array (float32)
            zin: Input array (float32)
            lat: Latitude coordinates array (float64)
            lon: Longitude coordinates array (float64)
            npts: Number of points
            opt: Optional interpolation options
            
        Returns:
            int: Status code
        """
        for arr in [zin, zout]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float32:
                raise TypeError("Data arrays must be float32 numpy arrays")
        for arr in [lat, lon]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float64:
                raise TypeError("Lat/Lon arrays must be float64 numpy arrays")
        return _llval(self._ptr, ctypes.byref(opt) if opt else None,
                      zout, zin, lat, lon, npts)

    def xywdval(self, uuout: numpy.ndarray, vvout: numpy.ndarray,
                uuin: numpy.ndarray, vvin: numpy.ndarray,
                x: numpy.ndarray, y: numpy.ndarray, npts: int, opt: GeoOptions = None) -> int:
        """Get wind direction values at xy coordinates.
        
        From: georef_xywdval_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_xywdval(GeoRef* ref, GeoOptions* opt,
                                          float* uuout, float* vvout,
                                          float* uuin, float* vvin,
                                          double* x, double* y,
                                          int32_t npts)
        
        Args:
            uuout: Output U component array (float32)
            vvout: Output V component array (float32)
            uuin: Input U component array (float32)
            vvin: Input V component array (float32)
            x: X coordinates array (float64)
            y: Y coordinates array (float64)
            npts: Number of points
            opt: Optional interpolation options
            
        Returns:
            int: Status code
        """
        for arr in [uuin, vvin, uuout, vvout]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float32:
                raise TypeError("UV arrays must be float32 numpy arrays")
        for arr in [x, y]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float64:
                raise TypeError("X/Y arrays must be float64 numpy arrays")
        return _xywdval(self._ptr, ctypes.byref(opt) if opt else None,
                        uuout, vvout, uuin, vvin, x, y, npts)

    def xyuvval(self, uuout: numpy.ndarray, vvout: numpy.ndarray,
                uuin: numpy.ndarray, vvin: numpy.ndarray,
                x: numpy.ndarray, y: numpy.ndarray, npts: int, opt: GeoOptions = None) -> int:
        """Get UV vector values at xy coordinates.
        
        From: georef_xyuvval_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_xyuvval(GeoRef* ref, GeoOptions* opt,
                                          float* uuout, float* vvout,
                                          float* uuin, float* vvin,
                                          double* x, double* y,
                                          int32_t npts)
        
        Args:
            uuout: Output U component array (float32)
            vvout: Output V component array (float32)
            uuin: Input U component array (float32)
            vvin: Input V component array (float32)
            x: X coordinates array (float64)
            y: Y coordinates array (float64)
            npts: Number of points
            opt: Optional interpolation options
            
        Returns:
            int: Status code
        """
        for arr in [uuin, vvin, uuout, vvout]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float32:
                raise TypeError("UV arrays must be float32 numpy arrays")
        for arr in [x, y]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float64:
                raise TypeError("X/Y arrays must be float64 numpy arrays")
        return _xyuvval(self._ptr, ctypes.byref(opt) if opt else None,
                        uuout, vvout, uuin, vvin, x, y, npts)

    def xyval(self, zout: numpy.ndarray, zin: numpy.ndarray,
              x: numpy.ndarray, y: numpy.ndarray, npts: int, opt: GeoOptions = None) -> int:
        """Get values at xy coordinates.
        
        From: georef_xyval_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_xyval(GeoRef* ref, GeoOptions* opt,
                                        float* zout, float* zin,
                                        double* x, double* y,
                                        int32_t npts)
        
        Args:
            zout: Output array (float32)
            zin: Input array (float32)
            x: X coordinates array (float64)
            y: Y coordinates array (float64)
            npts: Number of points
            opt: Optional interpolation options
            
        Returns:
            int: Status code
        """
        for arr in [zin, zout]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float32:
                raise TypeError("Data arrays must be float32 numpy arrays")
        for arr in [x, y]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float64:
                raise TypeError("X/Y arrays must be float64 numpy arrays")
        return _xyval(self._ptr, ctypes.byref(opt) if opt else None,
                      zout, zin, x, y, npts)

    def getll(self, lat: numpy.ndarray, lon: numpy.ndarray) -> int:
        """Get lat/lon arrays for the grid.
        
        From: georef_getll_f in GeoRef_mod.F90
        C signature: int32_t georef_getll(GeoRef* ref,
                                        double* lat, double* lon)
        
        Args:
            lat: Output latitude array (float64)
            lon: Output longitude array (float64)
            
        Returns:
            int: Status code
        """
        for arr in [lat, lon]:
            if not isinstance(arr, numpy.ndarray) or arr.dtype != numpy.float64:
                raise TypeError("Lat/Lon arrays must be float64 numpy arrays")
        return _getll(self._ptr, lat, lon)

def Interp(refto, reffrom, options, zout, zin):
    if not isinstance(zin, numpy.ndarray) or not isinstance(zout, numpy.ndarray):
        raise TypeError("Input and output must be numpy arrays")
    if zin.dtype != numpy.float32 or zout.dtype != numpy.float32:
        raise TypeError("Arrays must be float32")
    # Add shape validation here
    return _interp(refto._ptr, reffrom._ptr, ctypes.byref(options), zout, zin)
