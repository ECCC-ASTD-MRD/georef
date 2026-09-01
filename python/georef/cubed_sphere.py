"""C function bindings for GeoRef library."""

import ctypes

# Load the shared library
# libgeoref = ctypes.CDLL("libgeoref.so")
from .shared_lib import libgeoref

encode_ig4 = libgeoref.encode_cs_ig4
encode_ig4.argtypes = [
    ctypes.c_int32,
    ctypes.c_int32,
]
encode_ig4.restype = ctypes.c_int32

_decode_ig4 = libgeoref.decode_cs_ig4
_decode_ig4.argtypes = [
    ctypes.c_int32,
    ctypes.POINTER(ctypes.c_int32),
    ctypes.POINTER(ctypes.c_int32),
]
_decode_ig4.restype = ctypes.c_void_p

encode_angle = libgeoref.encode_cs_angle
encode_angle.argtypes = [ctypes.c_double]
encode_angle.restype = ctypes.c_int32

decode_angle = libgeoref.decode_cs_angle
decode_angle.argtypes = [ctypes.c_int32]
decode_angle.restype = ctypes.c_double


def decode_ig4(ig4: int):
    ni = ctypes.c_int32(0)
    nj = ctypes.c_int32(0)
    _decode_ig4(ig4, ctypes.byref(ni), ctypes.byref(nj))
    return ni.value, nj.value
