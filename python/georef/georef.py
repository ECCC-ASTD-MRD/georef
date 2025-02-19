import ctypes
import numpy
import numpy.ctypeslib
from typing import Tuple
from .constants import *
from .structs import GeoOptions, GeoRefError

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
        if ptr.contents == 0:
            raise GeoRefError("Failure in C function GeoRef_Create")
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
        return (float(lat0.value), float(lon0.value), 
                float(lat1.value), float(lon1.value))

    # INT32_T GeoRef_Valid(TGeoRef *Ref)
    def valid(self) -> bool:
        """Check if the georef object is properly initialized.

        This method wraps the libgeoref function GeoRef_Valid() found in src/GeoRef.c.
        Verifies that the georef object has been properly initialized and contains valid data.

        Returns:
            bool: True if the georef is valid and initialized, False otherwise
        """
        return bool(_valid(self.ptr))

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

        Raises:
            TypeError: If input/output arrays are not float32 numpy arrays
        """
        if not isinstance(zin, numpy.ndarray) or not isinstance(zout, numpy.ndarray):
            raise TypeError("Input and output must be numpy arrays")
        if zin.dtype != numpy.float32 or zout.dtype != numpy.float32:
            raise TypeError("Arrays must be float32")
        # Add shape validation here

        opt_ptr = ctypes.c_void_p(None)
        if options is not None:
            opt_ptr = options.ptr

        return bool(_interp(self.ptr, reffrom.ptr, opt_ptr, zout, zin) == GEOREF_SUCCESS)

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
        out = GeoRef(self._ptr.contents.ni, self._ptr.contents.nj, self._ptr.contents.grtyp,
                     self._ptr.contents.ig1, self._ptr.contents.ig2, self._ptr.contents.ig3,
                     self._ptr.contents.ig4, self._ptr.contents.fst_file)
        if hard:
            out._ptr = _copy(self._ptr)
        else:
            out._ptr = _hardcopy(self._ptr)
        return out

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
            The underlying C function returns an int32_t where 1 indicates equality
            and 0 indicates inequality.
        """
        val = _equal(self.ptr, other.ptr)
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
            The underlying C function returns an INT32_T where 1 indicates containment
            and 0 indicates non-containment.
        """
        val = _within(self.ptr, other.ptr)
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
            The underlying C function returns an INT32_T where 1 indicates within range
            and 0 indicates outside range.
        """
        cin = 1 if inside else 0
        val = _withinrange(self.ptr, lat0, lon0, lat1, lon1, cin)
        return val == 1

    # INT32_T GeoRef_Intersect(const GeoRef_t *ref1, const GeoRef_t *ref2, INT32_T *x0, INT32_T *y0, INT32_T *x1, INT32_T *y1, INT32_T bd)
    def intersect(self, other, boundary=False):
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
            The underlying C function returns an INT32_T where 1 indicates intersection
            and 0 indicates no intersection.
        """
        x0 = ctypes.c_int32()
        y0 = ctypes.c_int32()
        x1 = ctypes.c_int32()
        y1 = ctypes.c_int32()
        bd_val = 1 if boundary else 0
        
        val = _intersect(self.ptr, 
                         other.ptr, 
                         ctypes.byref(x0), 
                         ctypes.byref(y0), 
                         ctypes.byref(x1), 
                         ctypes.byref(y1), 
                         bd_val
                        )
        
        return (val == 1, (x0.value, y0.value, x1.value, y1.value))

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
            The underlying C function returns an INT32_T where 1 indicates success
            and 0 indicates failure.
        """
        i0 = ctypes.c_double()
        j0 = ctypes.c_double()
        i1 = ctypes.c_double()
        j1 = ctypes.c_double()
        
        val = _boundingbox(self.ptr, lat0, lon0, lat1, lon1,
                           ctypes.byref(i0), ctypes.byref(j0), ctypes.byref(i1), ctypes.byref(j1))
        
        return (val == 1, (i0.value, j0.value, i1.value, j1.value))

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
            The underlying C function returns an INT32_T where 1 indicates success
            and 0 indicates failure.
        """
        name_str = (name + '\0').encode('utf-8') if name else b'\0'
        val = _write(self.ptr, name_str, file.get_c_ptr())
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
            This method modifies the current georef instance by updating its internal pointer.
            The underlying C function returns a GeoRef_t* pointer, which is NULL on failure.
        """
        self.ptr = _createfromrecord(record.get_c_ptr())
        return bool(ctypes.cast(self.ptr, ctypes.c_void_p).value is not None)

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
            The underlying C function returns an INT32_T where 1 indicates success
            and 0 indicates failure.
        """
        uu_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        vv_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        opt_ptr = options.ptr if options is not None else GeoOptions
        
        val = _interpuv(self.ptr, ref_from.ptr, opt_ptr,
                                 uu_out, vv_out, uu_in, vv_in)
        
        return (val == 1, (uu_out, vv_out))

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
            The underlying C function returns an INT32_T where 1 indicates success
            and 0 indicates failure.
        """
        uu_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        vv_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        opt_ptr = options.ptr if options is not None else GeoOptions
        
        val = _interpwd(self.ptr, ref_from.ptr, opt_ptr,
                                 uu_out, vv_out, uu_in, vv_in)
        
        return (val == 1, (uu_out, vv_out))

    # INT32_T GeoRef_UV2WD(TGeoRef *Ref, float *spd_out, float *wd_out, 
    #                      const float *uuin, const float *vvin,
    #                      const double *Lat, const double *Lon, INT32_T Nb)
    def uv2wd(self, uu_in, vv_in, lat, lon) -> tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]:
        """Convert grid winds (UU/VV) to meteorological winds (speed/direction).

        This method wraps the libgeoref function GeoRef_UV2WD() found in src/GeoRef_InterpUV.c.

        Args:
            uu_in (numpy.ndarray): Input U component array (float32)
            vv_in (numpy.ndarray): Input V component array (float32)
            lat (numpy.ndarray): Latitude values array (float64)
            lon (numpy.ndarray): Longitude values array (float64)

        Returns:
            tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - bool: True if conversion succeeded, False otherwise
                - tuple[numpy.ndarray, numpy.ndarray]: The converted (speed, direction) arrays.
                                                      Arrays are float32 type.
        """
        npts = len(uu_in)
        spd_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        wd_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        val = _ud2wd(self.ptr, spd_out, wd_out, uu_in, vv_in, 
                              lat, lon, ctypes.c_int32(npts))
        
        return (val == 0, (spd_out, wd_out))  # Note: C function returns 0 for success

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
            The underlying C function returns 0 for success and -1 for failure.
        """
        npts = len(spd_in)
        uu_out = numpy.empty_like(spd_in, dtype=numpy.float32)
        vv_out = numpy.empty_like(dir_in, dtype=numpy.float32)
        
        val = _wd2uv(self.ptr, uu_out, vv_out, spd_in, dir_in,
                              lat, lon, ctypes.c_int32(npts))
        
        return (val == 0, (uu_out, vv_out))

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
        """
        npts = len(uu_in)
        uu_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        vv_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        val = _uv2uv(self.ptr, uu_out, vv_out, uu_in, vv_in,
                            lat, lon, ctypes.c_int32(npts))
        
        return (val == 1, (uu_out, vv_out))  # Note: C function returns TRUE (1) for success

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
        """
        npts = len(lat)
        spd_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        dir_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        opt_ptr = options.ptr if options is not None else GeoOptions
        
        val = _llwdval(self.ptr, opt_ptr, spd_out, dir_out, 
                                uu_in, vv_in, lat, lon, ctypes.c_int32(npts))
        
        return (val == 0, (spd_out, dir_out))  # Note: C function returns 0 for success

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
        """
        npts = len(lat)
        uu_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        vv_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        opt_ptr = options.ptr if options is not None else GeoOptions
        
        val = _lluvval(self.ptr, opt_ptr, uu_out, vv_out,
                                uu_in, vv_in, lat, lon, ctypes.c_int32(npts))
        
        return (val == 0, (uu_out, vv_out))  # Note: C function returns 0 for success

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
        """
        npts = len(lat)
        z_out = numpy.empty_like(z_in, dtype=numpy.float32)
        
        opt_ptr = options.ptr if options is not None else GeoOptions
        
        val = _llval(self.ptr, opt_ptr, z_out, z_in,
                              lat, lon, ctypes.c_int32(npts))
        
        return (val == 0, z_out)  # Note: C function returns 0 for success

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
        npts = len(x)
        spd_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        dir_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        opt_ptr = options.ptr if options is not None else GeoOptions
        
        val = _xywdval(self.ptr, opt_ptr, spd_out, dir_out,
                                uu_in, vv_in, x, y, ctypes.c_int32(npts))
        
        return (val == 0, (spd_out, dir_out))  # Note: C function returns 0 for success

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
        npts = len(x)
        uu_out = numpy.empty_like(uu_in, dtype=numpy.float32)
        vv_out = numpy.empty_like(vv_in, dtype=numpy.float32)
        
        opt_ptr = options.ptr if options is not None else GeoOptions
        
        val = _xyuvval(self.ptr, opt_ptr, uu_out, vv_out,
                                uu_in, vv_in, x, y, ctypes.c_int32(npts))
        
        return (val == 0, (uu_out, vv_out))  # Note: C function returns 0 for success

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
        npts = len(x)
        z_out = numpy.empty_like(z_in, dtype=numpy.float32)
        
        opt_ptr = options.ptr if options is not None else GeoOptions
        
        val = _xyval(self.ptr, opt_ptr, z_out, z_in,
                              x, y, ctypes.c_int32(npts))
        
        return (val == 0, z_out)  # Note: C function returns 0 for success

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
        """
        npts = len(lat)
        x_out = numpy.empty_like(lat, dtype=numpy.float64)
        y_out = numpy.empty_like(lon, dtype=numpy.float64)
        
        val = _ll2xy(self.ptr, x_out, y_out, lat, lon,
                              ctypes.c_int32(npts), ctypes.c_int32(extrapolate))
        
        return (val == npts, (x_out, y_out))  # Note: Returns count of valid conversions

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
        """
        npts = len(x)
        lat_out = numpy.empty_like(x, dtype=numpy.float64)
        lon_out = numpy.empty_like(y, dtype=numpy.float64)
        
        val = _xy2ll(self.ptr, lat_out, lon_out, x, y,
                              ctypes.c_int32(npts), ctypes.c_int32(extrapolate))
        
        return (val == npts, (lat_out, lon_out))  # Note: Returns count of valid conversions

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
        return _xydistance(self.ptr, 
                                    ctypes.c_double(x0), ctypes.c_double(y0),
                                    ctypes.c_double(x1), ctypes.c_double(y1))

    
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
        return float(_lldistance(self.ptr,
                                    ctypes.c_double(lat0), ctypes.c_double(lon0),
                                    ctypes.c_double(lat1), ctypes.c_double(lon1)))

    # INT32_T GeoRef_GetLL(TGeoRef *Ref, double *Lat, double *Lon)
    def getll(self):
        """Get latitude and longitude coordinates for all grid points.

        This method wraps the libgeoref function GeoRef_GetLL() found in src/GeoRef_InterpLL.c.
        Calculates or retrieves the lat/lon coordinates for all grid points, handling both
        simple grids and grids with sub-references (yin/yan).

        Returns:
            tuple[bool, tuple[numpy.ndarray, numpy.ndarray]]: A tuple containing:
                - bool: True if coordinates were retrieved successfully
                - tuple[numpy.ndarray, numpy.ndarray]: The (lat, lon) coordinate arrays.
                                                      Arrays are float64 type.
        """
        npts = self.nx * self.ny
        lat = numpy.empty(npts, dtype=numpy.float64)
        lon = numpy.empty(npts, dtype=numpy.float64)
        
        val = _getll(self.ptr, lat, lon)
        
        return (val >= 0, (lat, lon))  # Note: Returns -1 on error, number of points otherwise

    
