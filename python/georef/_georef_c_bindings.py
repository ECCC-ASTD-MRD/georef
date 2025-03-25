"""C function bindings for GeoRef library."""
import ctypes
import numpy

# Load the shared library
# libgeoref = ctypes.CDLL("libgeoref.so")
from .shared_lib import libgeoref
from .structs import GeoOptions, GeoRefError

_new = libgeoref.GeoRef_New
_new.argtypes = []
_new.restype = ctypes.c_void_p


_valid = libgeoref.GeoRef_Valid
_valid.argtypes = [ctypes.c_void_p]
_valid.restype = ctypes.c_int

_georef_create = libgeoref.GeoRef_Create
_georef_create.argtypes = (
    ctypes.c_int, ctypes.c_int, ctypes.c_char_p,
    ctypes.c_int, ctypes.c_int,
    ctypes.c_int, ctypes.c_int, ctypes.c_int
    )
_georef_create.restype = ctypes.c_void_p

_georef_limits = libgeoref.GeoRef_Limits
_georef_limits.argtypes = (ctypes.c_void_p,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double)
)
_georef_limits.restype = ctypes.c_int32

_interp = libgeoref.GeoRef_Interp
_interp.argtypes = (
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.POINTER(GeoOptions),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32)
)
_interp.restype = ctypes.c_int

_free = libgeoref.GeoRef_Free
_free.argtypes = [ctypes.c_void_p]
_free.restype = None

_copy = libgeoref.GeoRef_Copy
_copy.argtypes = [ctypes.c_void_p]
_copy.restype = ctypes.c_void_p

_hardcopy = libgeoref.GeoRef_HardCopy
_hardcopy.argtypes = [ctypes.c_void_p]
_hardcopy.restype = ctypes.c_void_p

_equal = libgeoref.GeoRef_Equal
_equal.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_equal.restype = ctypes.c_int32

_within = libgeoref.GeoRef_Within
_within.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_within.restype = ctypes.c_int32

_withinrange = libgeoref.GeoRef_WithinRange
_withinrange.argtypes = [ctypes.c_void_p,
                        ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double,
                        ctypes.c_int32]
_withinrange.restype = ctypes.c_int32

_intersect = libgeoref.GeoRef_Intersect
_intersect.argtypes = [ctypes.c_void_p, ctypes.c_void_p,
                      ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32),
                      ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32),
                      ctypes.c_int32]
_intersect.restype = ctypes.c_int32

_boundingbox = libgeoref.GeoRef_BoundingBox
_boundingbox.argtypes = [ctypes.c_void_p,
                        ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double,
                        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
_boundingbox.restype = ctypes.c_int32

_write = libgeoref.GeoRef_Write
_write.argtypes = [
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_void_p
]
_write.restype = ctypes.c_int32

_createfromrecord = libgeoref.GeoRef_CreateFromRecord
_createfromrecord.argtypes = [ctypes.c_void_p]
_createfromrecord.restype = ctypes.c_void_p

_interpuv = libgeoref.GeoRef_InterpUV
_interpuv.argtypes = (
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.POINTER(GeoOptions),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32)
)
_interpuv.restype = ctypes.c_int

_interpwd = libgeoref.GeoRef_InterpWD
_interpwd.argtypes = (
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.POINTER(GeoOptions),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32),
    numpy.ctypeslib.ndpointer(dtype=numpy.float32)
)
_interpwd.restype = ctypes.c_int

_ud2wd = libgeoref.GeoRef_UV2WD
_ud2wd.argtypes = [ctypes.c_void_p,
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   ctypes.c_int32]
_ud2wd.restype = ctypes.c_int32

_wd2uv = libgeoref.GeoRef_WD2UV
_wd2uv.argtypes = [ctypes.c_void_p,
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   ctypes.c_int32]
_wd2uv.restype = ctypes.c_int32

_uv2uv = libgeoref.GeoRef_UV2UV
_uv2uv.argtypes = [ctypes.c_void_p,
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   ctypes.c_int32]
_uv2uv.restype = ctypes.c_int32

_llwdval = libgeoref.GeoRef_LLWDVal
_llwdval.argtypes = [ctypes.c_void_p, ctypes.POINTER(GeoOptions),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     ctypes.c_int32]
_llwdval.restype = ctypes.c_int32

_lluvval = libgeoref.GeoRef_LLUVVal
_lluvval.argtypes = [ctypes.c_void_p, ctypes.POINTER(GeoOptions),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     ctypes.c_int32]
_lluvval.restype = ctypes.c_int32

_llval = libgeoref.GeoRef_LLVal
_llval.argtypes = [ctypes.c_void_p, ctypes.POINTER(GeoOptions),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   ctypes.c_int32]
_llval.restype = ctypes.c_int32

_xywdval = libgeoref.GeoRef_XYWDVal
_xywdval.argtypes = [ctypes.c_void_p, ctypes.POINTER(GeoOptions),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     ctypes.c_int32]
_xywdval.restype = ctypes.c_int32

_xyuvval = libgeoref.GeoRef_XYUVVal
_xyuvval.argtypes = [ctypes.c_void_p, ctypes.POINTER(GeoOptions),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                     ctypes.c_int32]
_xyuvval.restype = ctypes.c_int32

_xyval = libgeoref.GeoRef_XYVal
_xyval.argtypes = [ctypes.c_void_p, ctypes.POINTER(GeoOptions),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float32),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   ctypes.c_int32]
_xyval.restype = ctypes.c_int32

_ll2xy = libgeoref.GeoRef_LL2XY
_ll2xy.argtypes = [ctypes.c_void_p,
                   ctypes.c_double, ctypes.c_double,
                   ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
_ll2xy.restype = ctypes.c_int32

_xy2ll = libgeoref.GeoRef_XY2LL
_xy2ll.argtypes = [ctypes.c_void_p,
                   ctypes.c_double, ctypes.c_double,
                   ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
_xy2ll.restype = ctypes.c_int32

_xydistance = libgeoref.GeoRef_XYDistance
_xydistance.argtypes = [ctypes.c_void_p,
                        ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double]
_xydistance.restype = ctypes.c_double

_lldistance = libgeoref.GeoRef_LLDistance
_lldistance.argtypes = [ctypes.c_void_p,
                        ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double]
_lldistance.restype = ctypes.c_double

_getll = libgeoref.GeoRef_GetLL
_getll.argtypes = [ctypes.c_void_p,
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64),
                   numpy.ctypeslib.ndpointer(dtype=numpy.float64)]
_getll.restype = ctypes.c_int32

_def_create = libgeoref.Def_Create
_def_create.argtypes = [
    ctypes.c_int32,                  # NI
    ctypes.c_int32,                  # NJ
    ctypes.c_int32,                  # NK
    ctypes.c_int32,                  # Type (TDef_Type)
    ctypes.POINTER(ctypes.c_char),   # Comp0
    ctypes.POINTER(ctypes.c_char),   # Comp1
    ctypes.POINTER(ctypes.c_char)    # Mask
]
_def_create.restype = ctypes.POINTER(_TDef) # Ask Mr. Carphin

_geoset_writefst = libgeoref.GeoRef_SetWriteFST
_geoset_writefst.argtypes = [ctypes.c_void_p, ctypes.c_void_p] # Ask Mr. Carphin
_geoset_writefst.restype = ctypes.c_int32

_geoset_readfst = libgeoref.GeoRef_SetReadFST
_geoset_readfst.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int32, ctypes.c_void_p]
_geoset_readfst.restype = ctypes.c_int32

# _grid_value = libgeoref.GeoRef_GridValue
# _grid_value.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double]
# _grid_value.restype = ctypes.c_double

# _ll_value = libgeoref.GeoRef_LLValue
# _ll_value.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double]
# _ll_value.restype = ctypes.c_double

# _xy_to_ll = libgeoref.GeoRef_XYToLL
# _xy_to_ll.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double,
#                       ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
# _xy_to_ll.restype = ctypes.c_int



