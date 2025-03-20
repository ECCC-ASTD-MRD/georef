import ctypes
import numpy
import numpy.ctypeslib
from rmn import fst24_file
from typing import Tuple
from .constants import *
from .structs import GeoOptions, GeoRefError, GeoDef, GeoSet

from ._georef_c_bindings import (
    _valid,
    _georef_create,
    _georef_limits,
    _interp,
    _free,
    _copy,
    _hardcopy,
    _equal,
    _within,
    _withinrange,
    _intersect,
    _boundingbox,
    _write,
    _createfromrecord,
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
    _ll2xy,
    _xy2ll,
    _xydistance,
    _lldistance,
    _getll,
    _def_create,
    _geoset_writefst,
    _geoset_readfst,
    _new
)


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
        if ptr is None:
            raise GeoRefError("Failure in C function GeoRef_Create")
        self.shape = (ni, nj)
        self._ptr = ptr

    def limits(self) -> Tuple[float, float, float, float]:
        """Get the geographical limits of the GeoRef.
        
        From: georef_limits_f in GeoRef_mod.F90
        Source: src/GeoRef.c
        C signature: int32_t georef_limits(GeoRef* ref,
                                         double* lat0, double* lon0,
                                         double* lat1, double* lon1)
        
        Returns:
            tuple: (lat0, lon0, lat1, lon1) geographical limits
        
        Note:
            The underlying C function returns:
            - 0 for failure (NULL reference or error)
            - 1 for success
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
        if result == 0:
            raise GeoRefError(f"Failure in C function GeoRef_Limits: {result}")
        return lat0.value, lon0.value, lat1.value, lon1.value

    # INT32_T GeoRef_Valid(TGeoRef *Ref)
    def valid(self) -> bool:
        """Check if the georef object is properly initialized.

        This method wraps the libgeoref function GeoRef_Valid() found in src/GeoRef.c.
        Verifies that the georef object has been properly initialized and contains valid data.

        Returns:
            bool: True if the georef is valid and initialized, False otherwise

        Note:
            The underlying C function returns:
            - 0 for invalid/uninitialized reference
            - 1 for valid reference
        """
        return bool(_valid(self._ptr))

    def __del__(self):
        if hasattr(self, '_ptr'):
            _free(self._ptr)

    # INT32_T GeoRef_Interp(TGeoRef *RefTo, TGeoRef *RefFrom, TGeoOptions *Opt, float *zout, float *zin)
    def interp(self, reffrom, zout, zin, options=None) -> bool:
        """Interpolate values between two georefs.

        This method wraps the libgeoref function GeoRef_Interp() found in src/GeoRef_Interp.c.
        Interpolates data from source grid (reffrom) to destination grid (self).

        Args:
            reffrom (georef): Source georef object
            zout (numpy.ndarray): Output array for interpolated values (float32)
            zin (numpy.ndarray): Input array with source values (float32)
            options (geooptions, optional): Interpolation options. Uses default if None.

        Returns:
            bool: True if interpolation succeeded, False otherwise

        Note:
            The underlying C function returns:
            - GEOREF_SUCCESS (0) for successful interpolation
            - Non-zero value for failure
        """
        if not isinstance(zin, numpy.ndarray) or not isinstance(zout, numpy.ndarray):
            raise TypeError("Input and output must be numpy arrays")
        if zin.dtype != numpy.float32 or zout.dtype != numpy.float32:
            raise TypeError("Arrays must be float32")
        # Add shape validation here

        opt_ptr = ctypes.c_void_p(None)
        if options is not None:
            opt_ptr = ctypes.byref(options)

        result = _interp(self._ptr, reffrom._ptr, opt_ptr, zout, zin)
        if result != GEOREF_SUCCESS:
            raise GeoRefError("Failed to interpolate")
        

    def copy(self, hard=False):
        """Create a copy of the current georef object.

        This method wraps the libgeoref functions georef_copy() and georef_hardcopy() 
        found in src/GeoRef.c

        Args:
            hard (bool, optional): If True, creates a deep copy using georef_copy().
                                 If False, creates a shallow copy using georef_hardcopy().
                                 Defaults to False.

        Returns:
            Georef: A new georef object that is either a deep or shallow copy of the current instance.

        Note:
            The underlying C functions being called are:
            - GeoRef_Copy() for hard=True
            - GeoRef_HardCopy() for hard=False
        """
        
        if hard:
            out = GeoRef.__new__()
            out._ptr = _hardcopy(self._ptr)
            
            if out._ptr is None:
                raise GeoRefError("Failed to copy GeoRef: NULL pointer returned")
            return out
        else:
            return self


    # int32_t georef_equal(const georef_t *ref1, const georef_t *ref2)
    def equal(self, other):
        """Check if two georef objects are equal.

        This method wraps the libgeoref function georef_equal() found in src/GeoRef.c
        and src/GeoRef_mod.F90.

        Args:
            other (Georef): Another georef instance to compare with.

        Returns:
            bool: True if the georef objects are equal, False otherwise.

        Note:
            The underlying C function returns:
            - 0 for inequality or error (NULL references)
            - 1 for equality
        """
        val = _equal(self._ptr, other._ptr)
        return val == 1

    # INT32_T GeoRef_Within(const GeoRef_t *ref1, const GeoRef_t *ref2)
    def within(self, other):
        """Check if this georef object is within another georef object.

        This method wraps the libgeoref function GeoRef_Within() found in src/GeoRef.C
        and src/GeoRef_mod.F90.

        Args:
            other (Georef): Another georef instance to check containment against.

        Returns:
            bool: True if this georef object is within the other georef object,
                  False otherwise.

        Note:
            The underlying C function returns:
            - 0 for non-containment or error (NULL references)
            - 1 for containment
        """
        val = _within(self._ptr, other._ptr)
        return val == 1
    
    # INT32_T GeoRef_WithinRange(const GeoRef_t *ref, double lat0, double lon0, double lat1, double lon1, INT32_T in)
    def withinrange(self, lat0, lon0, lat1, lon1, inside=False):
        """Check if georef object is within a specified latitude/longitude range.

        This method wraps the libgeoref function GeoRef_WithinRange() found in src/GeoRef.C
        and src/GeoRef_mod.F90.

        Args:
            lat0 (float): First latitude point.
            lon0 (float): First longitude point.
            lat1 (float): Second latitude point.
            lon1 (float): Second longitude point.
            inside (bool, optional): Flag to modify range check behavior. Defaults to False.

        Returns:
            bool: True if the georef object is within the specified range,
                  False otherwise.

        Note:
            The underlying C function returns:
            - 0 for outside range or error (NULL reference)
            - 1 for within range
        """
        cin = 1 if inside else 0
        val = _withinrange(self._ptr, lat0, lon0, lat1, lon1, cin)
        return val == 1

    # INT32_T GeoRef_Intersect(const GeoRef_t *ref1, const GeoRef_t *ref2, INT32_T *x0, INT32_T *y0, INT32_T *x1, INT32_T *y1, INT32_T bd)
    def intersect(self, other, boundary=False) -> Tuple[bool, Tuple[int, int, int, int]]:
        """Check if two georef objects intersect and get intersection coordinates.

        This method wraps the libgeoref function GeoRef_Intersect() found in src/GeoRef.C
        and src/GeoRef_mod.F90.

        Args:
            other (Georef): Another georef instance to check intersection with.
            boundary (bool, optional): Flag to include boundary in intersection check. 
                                     Defaults to False.

        Returns:
            tuple[bool, tuple[int, int, int, int]]: A tuple containing:
                - bool: True if the georefs intersect, False otherwise
                - tuple[int, int, int, int]: The intersection coordinates (x0, y0, x1, y1).
                                            Returns (0, 0, 0, 0) if no intersection.

        Note:
            The underlying C function returns:
            - 0 for no intersection or error (NULL references)
            - 1 for intersection found
        """
        x0 = ctypes.c_int32()
        y0 = ctypes.c_int32()
        x1 = ctypes.c_int32()
        y1 = ctypes.c_int32()
        bd_val = 1 if boundary else 0
        
        val = _intersect(self._ptr, 
                         other._ptr, 
                         ctypes.byref(x0), 
                         ctypes.byref(y0), 
                         ctypes.byref(x1), 
                         ctypes.byref(y1), 
                         bd_val
                        )
        if val != 1:
            raise GeoRefError("Failed to check intersection")
        return x0.value, y0.value, x1.value, y1.value

    # INT32_T GeoRef_BoundingBox(const GeoRef_t *ref, double lat0, double lon0, double lat1, double lon1, double *i0, double *j0, double *i1, double *j1)
    def boundingbox(self, lat0, lon0, lat1, lon1) -> Tuple[bool, Tuple[float, float, float, float]]:
        """Calculate the bounding box coordinates for a given lat/lon range.

        This method wraps the libgeoref function GeoRef_BoundingBox() found in src/GeoRef.C
        and src/GeoRef_mod.F90.

        Args:
            lat0 (float): First latitude point.
            lon0 (float): First longitude point.
            lat1 (float): Second latitude point.
            lon1 (float): Second longitude point.

        Returns:
            tuple[bool, tuple[float, float, float, float]]: A tuple containing:
                - bool: True if bounding box calculation succeeded, False otherwise
                - tuple[float, float, float, float]: The bounding box coordinates (i0, j0, i1, j1).
                                                    Returns (0.0, 0.0, 0.0, 0.0) if calculation fails.

        Note:
            The underlying C function returns:
            - 0 for failure (NULL reference or calculation error)
            - 1 for successful calculation
        """
        i0 = ctypes.c_double()
        j0 = ctypes.c_double()
        i1 = ctypes.c_double()
        j1 = ctypes.c_double()
        
        val = _boundingbox(self._ptr, lat0, lon0, lat1, lon1,
                           ctypes.byref(i0), ctypes.byref(j0), 
                           ctypes.byref(i1), ctypes.byref(j1))
        
        if val != 1:
            raise GeoRefError("Failed to calculate bounding box")
        return i0.value, j0.value, i1.value, j1.value

    # INT32_T GeoRef_Write(const GeoRef_t *ref, const char *name, void *file_ptr)
    def write(self, file, name=None) -> bool:
        """Write the georef object to a file.

        This method wraps the libgeoref function GeoRef_Write() found in src/GeoRef.C
        and src/GeoRef_mod.F90.

        Args:
            file: FST file object containing the C file pointer.
            name (str, optional): Name identifier for the georef object.
                                If None, an empty string is used. Defaults to None.

        Returns:
            bool: True if writing succeeded, False otherwise.

        Note:
            The underlying C function returns:
            - 0 for failure (NULL references or write error)
            - 1 for successful write operation
        """
        name_str = (name + '\0').encode('utf-8') if name else b'\0'
        val = _write(self._ptr, name_str, file)
        return val == 1

    # TGeoRef* GeoRef_CreateFromRecord(fst_record_t *Rec)
    def fromrecord(self, record) -> bool:
        """Create a georef object from an FST record.

        This method wraps the libgeoref function GeoRef_CreateFromRecord() found in src/GeoRef.C
        and src/GeoRef_mod.F90.

        Args:
            record: FST record object containing the C record pointer.

        Returns:
            bool: True if creation succeeded, False otherwise.

        Note:
            The underlying C function returns:
            - NULL pointer (0) for failure
            - Valid GeoRef pointer (non-zero) for success
        """
        self._ptr = _createfromrecord(ctypes.byref(record))
        return bool(ctypes.cast(self._ptr, ctypes.c_void_p).value is not None)

    # INT32_T GeoRef_InterpUV(TGeoRef *RefTo, TGeoRef *RefFrom, TGeoOptions *Opt,
    #                         float *uuout, float *vvout, const float *uuin, const float *vvin)
    def interpuv(self, ref_from, uu_in, vv_in, options=None) -> Tuple[bool, Tuple[numpy.ndarray, numpy.ndarray]]:
        """Interpolate UV (vector) data from one georef to another.

        This method wraps the libgeoref function GeoRef_InterpUV() found in src/GeoRef_InterpUV.c.

        Args:
            ref_from (Georef): Source georef instance containing the input data.
            uu_in (numpy.ndarray): Input U component array (float32).
            vv_in (numpy.ndarray): Input V component array (float32).
            options (GeoOptions, optional): Interpolation options. Uses default if None.

        Returns:
            tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - bool: True if interpolation succeeded, False otherwise
                - tuple[numpy.ndarray, numpy.ndarray]: The interpolated (uu_out, vv_out) arrays.
                                                      Arrays are float32 type.

        Note:
            The underlying C function returns:
            - 0 for failure (NULL references or interpolation error)
            - 1 for successful interpolation
        """
        uu_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        vv_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        opt_ptr = ctypes.byref(options) if options is not None else GeoOptions
        
        val = _interpuv(self._ptr, 
                        ref_from._ptr, 
                        opt_ptr,
                        uu_out, 
                        vv_out, 
                        uu_in, 
                        vv_in)
        
        if val != 1:
            raise GeoRefError("Failed to interpolate U/V components")
        return uu_out, vv_out

    # INT32_T GeoRef_InterpWD(TGeoRef *RefTo, TGeoRef *RefFrom, TGeoOptions *Opt,
    #                         float *uuout, float *vvout, const float *uuin, const float *vvin)
    def interpwd(self, ref_from, uu_in, vv_in, options=None) -> Tuple[bool, Tuple[numpy.ndarray, numpy.ndarray]]:
        """Interpolate wind speed and direction between 2 georeferences.

        This method wraps the libgeoref function GeoRef_InterpWD() found in src/GeoRef_InterpUV.c.
        Interpolates vectorial component values between 2 georeferences and outputs rotated 
        speed and directions.

        Args:
            ref_from (Georef): Source georef instance containing the input data.
            uu_in (numpy.ndarray): Input speed values (float32).
            vv_in (numpy.ndarray): Input direction values (float32).
            options (GeoOptions, optional): Interpolation options. Uses default if None.

        Returns:
            tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - bool: True if interpolation succeeded, False otherwise
                - tuple[numpy.ndarray, numpy.ndarray]: The interpolated (speed, direction) arrays.
                                                      Arrays are float32 type.

        Note:
            The underlying C function returns:
            - 0 for failure (NULL references or interpolation error)
            - 1 for successful interpolation
        """
        uu_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        vv_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        opt_ptr = ctypes.byref(options) if options is not None else GeoOptions
        
        val = _interpwd(self._ptr, 
                        ctypes.c_void_p, 
                        opt_ptr,
                        uu_out, 
                        vv_out, 
                        uu_in, 
                        vv_in)
        
        if val != 1:
            raise GeoRefError("Failed to interpolate U/V components")
        return uu_out, vv_out

    # INT32_T GeoRef_UV2WD(TGeoRef *Ref, float *spd_out, float *wd_out, 
    #                      const float *uuin, const float *vvin,
    #                      const double *Lat, const double *Lon, INT32_T Nb)
    def uv2wd(self, uu_in, vv_in, lat, lon) -> tuple[numpy.ndarray, numpy.ndarray]:
        """Convert grid winds (UU/VV) to meteorological winds (speed/direction).

        This method wraps the libgeoref function GeoRef_UV2WD() found in src/GeoRef_InterpUV.c.

        Args:
            uu_in (numpy.ndarray): Input U component array (float32)
            vv_in (numpy.ndarray): Input V component array (float32)
            lat (numpy.ndarray): Latitude values array (float64)
            lon (numpy.ndarray): Longitude values array (float64)

        Returns:
            tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - tuple[numpy.ndarray, numpy.ndarray]: The converted (speed, direction) arrays.
                                                      Arrays are float32 type.

        Note:
            The underlying C function returns:
            - 0 for success
            - Non-zero for failure (NULL references or conversion error)
        """
        assert isinstance(uu_in, numpy.ndarray) and isinstance(vv_in, numpy.ndarray), GeoRefError("uu_in and vv_in must be numpy arrays")
        assert uu_in.shape == vv_in.shape, GeoRefError("uu_in and vv_in must have the same shape")
        assert len(lat) * len(lon) == uu_in.size, GeoRefError("uu_in size must match lat/lon grid size")
        npts = uu_in.size
        spd_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        wd_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        val = _ud2wd(self._ptr, 
                     spd_out, 
                     wd_out, 
                     uu_in, 
                     vv_in, 
                     lat, 
                     lon, 
                     npts)

        if val != 0:
            raise GeoRefError("Failed to convert UU/VV to speed/direction")
        return spd_out, wd_out

    # INT32_T GeoRef_WD2UV(TGeoRef *Ref, float *uugdout, float *vvgdout, 
    #                      const float *uullin, const float *vvllin,
    #                      const double *Lat, const double *Lon, INT32_T Nb)
    def wd2uv(self, spd_in, dir_in, lat, lon) -> tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]:
        """Convert meteorological winds (speed/direction) to grid winds (UU/VV).

        This method wraps the libgeoref function GeoRef_WD2UV() found in src/GeoRef_InterpUV.c.

        Args:
            spd_in (numpy.ndarray): Input wind speed values (float32)
            dir_in (numpy.ndarray): Input wind direction values (float32)
            lat (numpy.ndarray): Latitude values array (float64)
            lon (numpy.ndarray): Longitude values array (float64)

        Returns:
            tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - bool: True if conversion succeeded, False otherwise
                - tuple[numpy.ndarray, numpy.ndarray]: The converted (U, V) component arrays.
                                                      Arrays are float32 type.

        Note:
            The underlying C function returns:
            - 0 for success
            - Non-zero for failure (NULL references or conversion error)
        """
        assert isinstance(spd_in, numpy.ndarray) and isinstance(dir_in, numpy.ndarray), GeoRefError("spd_in and dir_in must be numpy arrays")
        assert spd_in.shape == dir_in.shape, GeoRefError("spd_in and dir_in must have the same shape")
        npts = spd_in.size
        uu_out = numpy.empty_like(spd_in, dtype=numpy.float32)
        vv_out = numpy.empty_like(dir_in, dtype=numpy.float32)
        
        val = _wd2uv(self._ptr, 
                     uu_out, 
                     vv_out, 
                     spd_in, 
                     dir_in,
                     lat, 
                     lon, 
                     npts)
        
        if val != 0:
            raise GeoRefError("Failed to convert speed/direction to UU/VV")
        return uu_out, vv_out

    # INT32_T GeoRef_UV2UV(TGeoRef *Ref, float *uullout, float *vvllout, 
    #                      const float *uuin, const float *vvin,
    #                      const double *Lat, const double *Lon, INT32_T Nb)
    def uv2uv(self, uu_in, vv_in, lat, lon) -> tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]:
        """Convert grid winds (UU/VV) between coordinate systems.

        This method wraps the libgeoref function GeoRef_UV2UV() found in src/GeoRef_InterpUV.c.
        Converts wind components between different coordinate systems while preserving the 
        vector quantities.

        Args:
            uu_in (numpy.ndarray): Input U component array (float32)
            vv_in (numpy.ndarray): Input V component array (float32)
            lat (numpy.ndarray): Latitude values array (float64)
            lon (numpy.ndarray): Longitude values array (float64)

        Returns:
            tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - bool: True if conversion succeeded, False otherwise
                - tuple[numpy.ndarray, numpy.ndarray]: The converted (U, V) component arrays.
                                                    Arrays are float32 type.

        Note:
            The underlying C function returns:
            - 0 for failure (NULL references or conversion error)
            - 1 for successful conversion
        """
        assert isinstance(uu_in, numpy.ndarray) and isinstance(vv_in, numpy.ndarray), GeoRefError("uu_in and vv_in must be numpy arrays")
        assert uu_in.shape == vv_in.shape, GeoRefError("uu_in and vv_in must have the same shape")
        npts = uu_in.size
        uu_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        vv_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        val = _uv2uv(self._ptr, 
                     uu_out, 
                     vv_out, 
                     uu_in, 
                     vv_in,
                     lat, 
                     lon, 
                     npts)
        
        if val != 1:
            raise GeoRefError("Failed to convert UU/VV between coordinate systems")
        return uu_out, vv_out

    # INT32_T GeoRef_LLWDVal(TGeoRef *Ref, TGeoOptions *Opt, float *uuout, float *vvout,
    #                        const float *uuin, const float *vvin,
    #                        const double *Lat, const double *Lon, INT32_T Nb)
    def llwdval(self, uu_in, vv_in, lat, lon, options=None) -> tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]:
        """Interpolate vector values as speed and directions at lat/lon positions.

        This method wraps the libgeoref function GeoRef_LLWDVal() found in src/GeoRef_InterpLL.c.
        Interpolates U/V components to lat/lon positions and converts them to speed/direction.

        Args:
            uu_in (numpy.ndarray): Input U component array (float32)
            vv_in (numpy.ndarray): Input V component array (float32)
            lat (numpy.ndarray): Target latitude values array (float64)
            lon (numpy.ndarray): Target longitude values array (float64)
            options (GeoOptions, optional): Interpolation options. Uses default if None.

        Returns:
            tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - bool: True if interpolation succeeded, False otherwise
                - tuple[numpy.ndarray, numpy.ndarray]: The interpolated (speed, direction) arrays.
                                                      Arrays are float32 type.

        Note:
            The underlying C function returns:
            - 0 for success
            - Non-zero for failure (NULL references or interpolation error)
        """
        assert isinstance(uu_in, numpy.ndarray) and isinstance(vv_in, numpy.ndarray), GeoRefError("uu_in and vv_in must be numpy arrays")
        assert uu_in.shape == vv_in.shape, GeoRefError("uu_in and vv_in must have the same shape")
        npts = uu_in.size
        spd_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        dir_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        opt_ptr = ctypes.byref(options) if options is not None else GeoOptions
        
        val = _llwdval(self._ptr, 
                     opt_ptr, 
                     spd_out, 
                     dir_out, 
                     uu_in, 
                     vv_in, 
                     lat, 
                     lon, 
                     npts)
        
        if val != 0:
            raise GeoRefError("Failed to interpolate U/V components to lat/lon positions")
        return spd_out, dir_out

    # INT32_T GeoRef_LLUVVal(TGeoRef *Ref, TGeoOptions *Opt, float *uuout, float *vvout,
    #                        const float *uuin, const float *vvin,
    #                        const double *Lat, const double *Lon, INT32_T Nb)
    def lluvval(self, uu_in, vv_in, lat, lon, options=None) -> tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]:
        """Interpolate vector values at lat/lon positions.

        This method wraps the libgeoref function GeoRef_LLUVVal() found in src/GeoRef_InterpLL.c.
        Interpolates U/V components to specified lat/lon positions.

        Args:
            uu_in (numpy.ndarray): Input U component array (float32)
            vv_in (numpy.ndarray): Input V component array (float32)
            lat (numpy.ndarray): Target latitude values array (float64)
            lon (numpy.ndarray): Target longitude values array (float64)
            options (GeoOptions, optional): Interpolation options. Uses default if None.

        Returns:
            tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - bool: True if interpolation succeeded, False otherwise
                - tuple[numpy.ndarray, numpy.ndarray]: The interpolated (U, V) arrays.
                                                      Arrays are float32 type.

        Note:
            The underlying C function returns:
            - 0 for success
            - Non-zero for failure (NULL references or interpolation error)
        """
        assert isinstance(uu_in, numpy.ndarray) and isinstance(vv_in, numpy.ndarray), GeoRefError("uu_in and vv_in must be numpy arrays")
        assert uu_in.shape == vv_in.shape, GeoRefError("uu_in and vv_in must have the same shape")
        npts = uu_in.size
        uu_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        vv_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        opt_ptr = ctypes.byref(options) if options is not None else GeoOptions
        
        val = _lluvval(self._ptr, 
                     opt_ptr, 
                     uu_out, 
                     vv_out,
                     uu_in, 
                    vv_in, 
                     lat, 
                     lon, 
                     npts)
        
        if val != 0:
            raise GeoRefError("Failed to interpolate U/V components to lat/lon positions")
        return uu_out, vv_out

    # INT32_T GeoRef_LLVal(TGeoRef *Ref, TGeoOptions *Opt, float *zout, float *zin,
    #                      const double *Lat, const double *Lon, INT32_T Nb)
    def llval(self, z_in, lat, lon, options=None) -> tuple[bool, numpy.ndarray]:
        """Interpolate values at lat/lon positions.

        This method wraps the libgeoref function GeoRef_LLVal() found in src/GeoRef_InterpLL.c.
        Interpolates scalar values from source grid to specified lat/lon positions.

        Args:
            z_in (numpy.ndarray): Input values array (float32)
            lat (numpy.ndarray): Target latitude values array (float64)
            lon (numpy.ndarray): Target longitude values array (float64)
            options (GeoOptions, optional): Interpolation options. Uses default if None.

        Returns:
            tuple[bool, numpy.ndarray]: A tuple containing:
                - bool: True if interpolation succeeded, False otherwise
                - numpy.ndarray: The interpolated values array (float32)

        Note:
            The underlying C function returns:
            - 0 for success
            - Non-zero for failure (NULL references or interpolation error)
        """
        assert isinstance(z_in, numpy.ndarray), GeoRefError("z_in must be a numpy array")
        npts = z_in.size
        z_out = numpy.empty_like(z_in, dtype=numpy.float32)
        
        opt_ptr = ctypes.byref(options) if options is not None else GeoOptions
        
        val = _llval(self._ptr, 
                     opt_ptr, 
                     z_out, 
                     z_in,
                     lat, 
                     lon, 
                     npts)
        
        if val != 0:
            raise GeoRefError("Failed to interpolate values at lat/lon positions")
        return z_out

    # INT32_T GeoRef_XYWDVal(TGeoRef *Ref, TGeoOptions *Opt, float *uuout, float *vvout,
    #                        const float *uuin, const float *vvin,
    #                        const double *X, const double *Y, INT32_T n)
    def xywdval(self, uu_in, vv_in, x, y, options=None) -> tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]:
        """Interpolate vector values as speed and direction at X/Y positions.

        This method wraps the libgeoref function GeoRef_XYWDVal() found in src/GeoRef_InterpXY.c.
        Interpolates U/V components to X/Y positions and converts them to speed/direction.

        Args:
            uu_in (numpy.ndarray): Input U component array (float32)
            vv_in (numpy.ndarray): Input V component array (float32)
            x (numpy.ndarray): Target X coordinates array (float64)
            y (numpy.ndarray): Target Y coordinates array (float64)
            options (GeoOptions, optional): Interpolation options. Uses default if None.

        Returns:
            tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - bool: True if interpolation succeeded, False otherwise
                - tuple[numpy.ndarray, numpy.ndarray]: The interpolated (speed, direction) arrays.
                                                      Arrays are float32 type.
        """
        assert isinstance(uu_in, numpy.ndarray) and isinstance(vv_in, numpy.ndarray), GeoRefError("uu_in and vv_in must be numpy arrays")
        assert uu_in.shape == vv_in.shape, GeoRefError("uu_in and vv_in must have the same shape")
        npts = uu_in.size
        spd_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        dir_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        opt_ptr = ctypes.byref(options) if options is not None else GeoOptions
        
        val = _xywdval(self._ptr, 
                     opt_ptr, 
                     spd_out, 
                     dir_out,
                     uu_in, 
                     vv_in, 
                     x, 
                     y, 
                     npts)
        
        if val != 0:
            raise GeoRefError("Failed to interpolate U/V components to X/Y positions")
        return spd_out, dir_out

    # INT32_T GeoRef_XYUVVal(TGeoRef *Ref, TGeoOptions *Opt, float *uuout, float *vvout,
    #                        const float *uuin, const float *vvin,
    #                        const double *X, const double *Y, INT32_T n)
    def xyuvval(self, uu_in, vv_in, x, y, options=None) -> tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]:
        """Interpolate vector values at X/Y positions.

        This method wraps the libgeoref function GeoRef_XYUVVal() found in src/GeoRef_InterpXY.c.
        Interpolates U/V components to specified X/Y positions while preserving vector quantities.

        Args:
            uu_in (numpy.ndarray): Input U component array (float32)
            vv_in (numpy.ndarray): Input V component array (float32)
            x (numpy.ndarray): Target X coordinates array (float64)
            y (numpy.ndarray): Target Y coordinates array (float64)
            options (GeoOptions, optional): Interpolation options. Uses default if None.

        Returns:
            tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - bool: True if interpolation succeeded, False otherwise
                - tuple[numpy.ndarray, numpy.ndarray]: The interpolated (U, V) arrays.
                                                      Arrays are float32 type.
        """
        assert isinstance(uu_in, numpy.ndarray) and isinstance(vv_in, numpy.ndarray), GeoRefError("uu_in and vv_in must be numpy arrays")
        assert uu_in.shape == vv_in.shape, GeoRefError("uu_in and vv_in must have the same shape")
        npts = uu_in.size
        uu_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        vv_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        opt_ptr = ctypes.byref(options) if options is not None else GeoOptions
        
        val = _xyuvval(self._ptr, 
                     opt_ptr, 
                     uu_out, 
                     vv_out,
                     uu_in, 
                     vv_in, 
                     x, 
                     y, 
                     npts)
        
        if val != 0:
            raise GeoRefError("Failed to interpolate U/V components to X/Y positions")
        return uu_out, vv_out

    # INT32_T GeoRef_XYVal(TGeoRef *Ref, TGeoOptions *Opt, float *zout, float *zin,
    #                      const double *X, const double *Y, INT32_T n)
    def xyval(self, z_in, x, y, options=None) -> tuple[bool, numpy.ndarray]:
        """Interpolate scalar values at X/Y positions.

        This method wraps the libgeoref function GeoRef_XYVal() found in src/GeoRef_InterpXY.c.
        Interpolates scalar values from source grid to specified X/Y positions.

        Args:
            z_in (numpy.ndarray): Input values array (float32)
            x (numpy.ndarray): Target X coordinates array (float64)
            y (numpy.ndarray): Target Y coordinates array (float64)
            options (GeoOptions, optional): Interpolation options. Uses default if None.

        Returns:
            tuple[bool, numpy.ndarray]: A tuple containing:
                - bool: True if interpolation succeeded, False otherwise
                - numpy.ndarray: The interpolated values array (float32)
        """
        assert isinstance(z_in, numpy.ndarray), GeoRefError("z_in must be a numpy array")
        npts = z_in.size
        z_out = numpy.empty_like(z_in, dtype=numpy.float32)
        
        opt_ptr = ctypes.byref(options) if options is not None else GeoOptions
        
        val = _xyval(self._ptr, 
                     opt_ptr, 
                     z_out, 
                     z_in,
                     x, 
                     y, 
                     npts)
        
        if val != 0:
            raise GeoRefError("Failed to interpolate values at X/Y positions")
        return z_out

    # INT32_T GeoRef_LL2XY(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon,
    #                      INT32_T Nb, INT32_T Extrap)
    def ll2xy(self, lat, lon, extrapolate=False) -> tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]:
        """Transform lat/lon coordinates to grid X/Y coordinates.

        This method wraps the libgeoref function GeoRef_LL2XY() found in src/GeoRef_InterpCoords.c.
        Converts geographic coordinates to grid coordinates, with optional extrapolation.

        Args:
            lat (numpy.ndarray): Input latitude values array (float64)
            lon (numpy.ndarray): Input longitude values array (float64)
            extrapolate (bool, optional): Enable extrapolation outside grid bounds.
                                        Defaults to False.

        Returns:
            tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - bool: True if all points were converted successfully, False otherwise
                - tuple[numpy.ndarray, numpy.ndarray]: The converted (X, Y) coordinate arrays.
                                                      Arrays are float64 type.

        Note:
            The underlying C function returns:
            - The number of successfully converted points (should equal input size for full success)
            - Less than input size indicates partial failure
            - 0 or negative value indicates complete failure
        """
        assert isinstance(lat, numpy.ndarray) and isinstance(lon, numpy.ndarray), GeoRefError("lat and lon must be numpy arrays")
        npts = lat.size
        x_out = numpy.empty_like(lat, dtype=numpy.float64)
        y_out = numpy.empty_like(lon, dtype=numpy.float64)
        
        val = _ll2xy(self._ptr, 
                     x_out, 
                     y_out, 
                     lat, 
                     lon,
                     npts, 
                     extrapolate)
        
        if val != npts:
            raise GeoRefError("Failed to convert lat/lon to grid X/Y coordinates")
        return x_out, y_out

    # INT32_T GeoRef_XY2LL(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y,
    #                      INT32_T Nb, INT32_T Extrap)
    def xy2ll(self, x, y, extrapolate=False) -> tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]:
        """Transform grid X/Y coordinates to lat/lon coordinates.

        This method wraps the libgeoref function GeoRef_XY2LL() found in src/GeoRef_InterpCoords.c.
        Converts grid coordinates to geographic coordinates, with optional extrapolation.

        Args:
            x (numpy.ndarray): Input X coordinate array (float64)
            y (numpy.ndarray): Input Y coordinate array (float64)
            extrapolate (bool, optional): Enable extrapolation outside grid bounds.
                                        Defaults to False.

        Returns:
            tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - bool: True if all points were converted successfully, False otherwise
                - tuple[numpy.ndarray, numpy.ndarray]: The converted (lat, lon) arrays.
                                                      Arrays are float64 type.
                                                      Invalid points set to -999.0 if not extrapolating.

        Note:
            The underlying C function returns:
            - The number of successfully converted points (should equal input size for full success)
            - Less than input size indicates partial failure
            - 0 or negative value indicates complete failure
        """
        assert isinstance(x, numpy.ndarray) and isinstance(y, numpy.ndarray), GeoRefError("x and y must be numpy arrays")
        npts = x.size
        lat_out = numpy.empty_like(x, dtype=numpy.float64)
        lon_out = numpy.empty_like(y, dtype=numpy.float64)
        
        val = _xy2ll(self._ptr, 
                     lat_out, 
                     lon_out, 
                     x, 
                     y,
                     npts, 
                     extrapolate)
        
        if val != npts:
            raise GeoRefError("Failed to convert grid X/Y to lat/lon coordinates")
        return lat_out, lon_out

    # double GeoRef_XYDistance(TGeoRef *Ref, double X0, double Y0, double X1, double Y1)
    def xydistance(self, x0, y0, x1, y1) -> float:
        """Calculate distance between two locations in grid coordinates.

        This method wraps the libgeoref function GeoRef_XYDistance() found in src/GeoRef_Func.c.
        Handles both geographic (lat/lon) and projected (UTM) coordinate systems.

        Args:
            x0 (float): First X coordinate
            y0 (float): First Y coordinate
            x1 (float): Second X coordinate
            y1 (float): Second Y coordinate

        Returns:
            float: Distance in meters between the two points
        """
        return float(_xydistance(self._ptr, x0, y0,x1,y1))

    
    # double GeoRef_LLDistance(TGeoRef *Ref, double Lat0, double Lon0, double Lat1, double Lon1)
    def lldistance(self, lat0, lon0, lat1, lon1) -> float:
        """Calculate great circle distance between two lat/lon points.

        This method wraps the libgeoref function GeoRef_LLDistance() found in src/GeoRef_Func.c.
        Calculates the great circle distance between two geographic coordinates.

        Args:
            lat0 (float): First latitude in degrees
            lon0 (float): First longitude in degrees
            lat1 (float): Second latitude in degrees
            lon1 (float): Second longitude in degrees

        Returns:
            float: Great circle distance in meters
        """
        return float(_lldistance(self._ptr, lat0, lon0, lat1, lon1))


    # INT32_T GeoRef_GetLL(TGeoRef *Ref, double *Lat, double *Lon)
    def getll(self) -> int:
        """Get lat/lon positions for all grid points.

        Args:
            lat: Output array for latitude values
            lon: Output array for longitude values

        Returns:
            int: Number of coordinates

        Note:
            This wraps GeoRef_GetLL from src/GeoRef_InterpLL.c
        """
        lat = numpy.empty(self.shape, dtype=numpy.float64) # Awaiting Mr. Gauthier opinion on this
        lon = numpy.empty(self.shape, dtype=numpy.float64)
        n = _getll(self._ptr, lat, lon)
        if n == -1:
            raise GeoRefError("Failed to get lat/lon coordinates: Missing descriptors")
        return lat, lon
    

class GeoDef:
    """Wrapper for the C GeoDef structure providing geographic definition functionality."""
    # TDef* Def_Create(int32_t NI, int32_t NJ, int32_t NK, int32_t NC, TDef_Type Type, int32_t Alias) 
    def __init__(self, ni: int, nj: int, nk: int, type: int, comp0: str, comp1: str, mask: str):
        """Initialize a new geographic definition.
        
        Args:
            ni (int): Number of points in x direction
            nj (int): Number of points in y direction
            nk (int): Number of points in z direction
            type (int): Type of definition
            comp0 (str): First component
            comp1 (str): Second component
            mask (str): Mask string
            
        Raises:
            GeoRefError: If initialization fails (NULL pointer returned)
        """
        self._ptr = _def_create(ni, nj, nk, type, comp0, comp1, mask)
        if not self._ptr or self._ptr is None:
            raise GeoRefError("Failed to create GeoDef: NULL pointer returned")

class GeoSet:
    """Wrapper for the C GeoSet structure providing geographic set functionality."""
    
    def __init__(self):
        """Initialize a new GeoSet instance with a null pointer."""
        self._ptr = None # Where does the pointer come from ?

    @staticmethod
    # int32_t GeoRef_SetReadFST(const TGeoRef * const RefTo, const TGeoRef * const RefFrom, const int32_t InterpType, const fst_file * const File)  
    def read_fst(cls, ref_to: GeoRef, ref_from: GeoRef, interp: int, file: fst24_file) -> None: # For all other functions that raise errors   
        """Read geographic set data from an FST file.
        
        This method wraps the libgeoref function GeoRef_SetReadFST found in src/GeoRef_Set.c.
        
        Args:
            ref_to (GeoRef): Target reference
            ref_from (GeoRef): Source reference
            interp (int): Interpolation type
            file (FSTFile): FST file to read from
            
        Returns:
            bool: True if successful, False otherwise
            
        Raises:
            GeoRefError: If reading fails or references are invalid
            
        Note:
            The underlying C function returns:
            - NULL pointer (0) for failure
            - Valid pointer for successful read operation
        """
        # Add raise error if not valid and return None
        result = _geoset_readfst(ref_to._ptr, ref_from._ptr, interp, file._c_ref)
        if result is None:
            raise GeoRefError("Failed to read GeoSet from FST file")
    
    # int32_t GeoRef_SetWriteFST(const TGeoSet * const GSet, fst_file * const File)
    def write_fst(self, file: fst24_file) -> bool:
        """Write geographic set data to an FST file.
        
        This method wraps the libgeoref function GeoRef_SetWriteFST found in src/GeoRef_Set.c.
        
        Args:
            file (FSTFile): FST file to write to
            
        Returns:
            bool: True if successful, False otherwise
            
        Raises:
            GeoRefError: If writing fails or GeoSet is uninitialized
            
        Note:
            The underlying C function returns:
            - 0 for failure (NULL references or write error)
            - 1 for successful write operation
        """
        if not self._ptr:
            raise GeoRefError("Cannot write uninitialized GeoSet")
            
        result = _geoset_writefst(self._ptr, file._c_ref)
        # If not valid,raise error
        if result != 1:
            raise GeoRefError("Failed to write GeoSet to FST file")
