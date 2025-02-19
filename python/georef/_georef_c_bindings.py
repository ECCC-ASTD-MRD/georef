"""C function bindings for GeoRef library."""
import ctypes
import numpy

# Load the shared library
# libgeoref = ctypes.CDLL("libgeoref.so")
from .shared_lib import libgeoref
from .structs import GeoOptions, _TGeoRef, GeoRefError


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

_copy = libgeoref.GeoRef_Copy
_copy.argtypes = [ctypes.POINTER(_TGeoRef)]
_copy.restype = ctypes.POINTER(_TGeoRef)

_hardcopy = libgeoref.GeoRef_HardCopy 
_hardcopy.argtypes = [ctypes.POINTER(_TGeoRef)]
_hardcopy.restype = ctypes.POINTER(_TGeoRef)

_equal = libgeoref.GeoRef_Equal
_equal.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.POINTER(_TGeoRef)]
_equal.restype = ctypes.c_int32

_within = libgeoref.GeoRef_Within
_within.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.POINTER(_TGeoRef)]
_within.restype = ctypes.c_int32

_withinrange = libgeoref.GeoRef_WithinRange
_withinrange.argtypes = [ctypes.POINTER(_TGeoRef), 
                        ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double,
                        ctypes.c_int32]
_withinrange.restype = ctypes.c_int32

_intersect = libgeoref.GeoRef_Intersect
_intersect.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.POINTER(_TGeoRef),
                      ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32),
                      ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32),
                      ctypes.c_int32]
_intersect.restype = ctypes.c_int32

_boundingbox = libgeoref.GeoRef_BoundingBox
_boundingbox.argtypes = [ctypes.POINTER(_TGeoRef),
                        ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double,
                        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
_boundingbox.restype = ctypes.c_int32

_write = libgeoref.GeoRef_Write
_write.argtypes = [
    ctypes.POINTER(_TGeoRef),
    ctypes.c_char_p,
    ctypes.c_void_p
]
_write.restype = ctypes.c_int32

_createfromrecord = libgeoref.GeoRef_CreateFromRecord
_createfromrecord.argtypes = [ctypes.c_void_p]
_createfromrecord.restype = ctypes.POINTER(_TGeoRef)

_interpuv = libgeoref.GeoRef_InterpUV
_interpuv.argtypes = (
    ctypes.POINTER(_TGeoRef),
    ctypes.POINTER(_TGeoRef),
    ctypes.POINTER(GeoOptions),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32)
)
_interpuv.restype = ctypes.c_int

_interpwd = libgeoref.GeoRef_InterpWD
_interpwd.argtypes = (
    ctypes.POINTER(_TGeoRef),
    ctypes.POINTER(_TGeoRef),
    ctypes.POINTER(GeoOptions),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32)
)
_interpwd.restype = ctypes.c_int

_ud2wd = libgeoref.GeoRef_UV2WD
_ud2wd.argtypes = [ctypes.POINTER(_TGeoRef),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   ctypes.c_int32]
_ud2wd.restype = ctypes.c_int32

_wd2uv = libgeoref.GeoRef_WD2UV
_wd2uv.argtypes = [ctypes.POINTER(_TGeoRef),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   ctypes.c_int32]
_wd2uv.restype = ctypes.c_int32

_uv2uv = libgeoref.GeoRef_UV2UV
_uv2uv.argtypes = [ctypes.POINTER(_TGeoRef),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   ctypes.c_int32]
_uv2uv.restype = ctypes.c_int32

_llwdval = libgeoref.GeoRef_LLWDVal
_llwdval.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.POINTER(GeoOptions),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     ctypes.c_int32]
_llwdval.restype = ctypes.c_int32

_lluvval = libgeoref.GeoRef_LLUVVal
_lluvval.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.POINTER(GeoOptions),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     ctypes.c_int32]
_lluvval.restype = ctypes.c_int32

_llval = libgeoref.GeoRef_LLVal
_llval.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.POINTER(GeoOptions),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   ctypes.c_int32]
_llval.restype = ctypes.c_int32

_xywdval = libgeoref.GeoRef_XYWDVal
_xywdval.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.POINTER(GeoOptions),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     ctypes.c_int32]
_xywdval.restype = ctypes.c_int32

_xyuvval = libgeoref.GeoRef_XYUVVal
_xyuvval.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.POINTER(GeoOptions),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     ctypes.c_int32]
_xyuvval.restype = ctypes.c_int32

_xyval = libgeoref.GeoRef_XYVal
_xyval.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.POINTER(GeoOptions),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   ctypes.c_int32]
_xyval.restype = ctypes.c_int32

_ll2xy = libgeoref.GeoRef_LL2XY
_ll2xy.argtypes = [ctypes.POINTER(_TGeoRef),
                   ctypes.c_double, ctypes.c_double,
                   ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
_ll2xy.restype = ctypes.c_int32

_xy2ll = libgeoref.GeoRef_XY2LL
_xy2ll.argtypes = [ctypes.POINTER(_TGeoRef),
                   ctypes.c_double, ctypes.c_double,
                   ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
_xy2ll.restype = ctypes.c_int32

_xydistance = libgeoref.GeoRef_XYDistance
_xydistance.argtypes = [ctypes.POINTER(_TGeoRef),
                        ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double]
_xydistance.restype = ctypes.c_double

_lldistance = libgeoref.GeoRef_LLDistance
_lldistance.argtypes = [ctypes.POINTER(_TGeoRef),
                        ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double]
_lldistance.restype = ctypes.c_double

_getll = libgeoref.GeoRef_GetLL
_getll.argtypes = [ctypes.POINTER(_TGeoRef),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64)]
_getll.restype = ctypes.c_int32


# _grid_value = libgeoref.GeoRef_GridValue
# _grid_value.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.c_double, ctypes.c_double]
# _grid_value.restype = ctypes.c_double

# _ll_value = libgeoref.GeoRef_LLValue
# _ll_value.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.c_double, ctypes.c_double]
# _ll_value.restype = ctypes.c_double

# _xy_to_ll = libgeoref.GeoRef_XYToLL
# _xy_to_ll.argtypes = [ctypes.POINTER(_TGeoRef), ctypes.c_double, ctypes.c_double, 
#                       ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
# _xy_to_ll.restype = ctypes.c_int



