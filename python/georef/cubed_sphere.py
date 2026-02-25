"""C function bindings for GeoRef library."""
import ctypes
import numpy

# Load the shared library
# libgeoref = ctypes.CDLL("libgeoref.so")
from .shared_lib import libgeoref

encodeig4 = libgeoref.encode_cs_ig4
encodeig4.argtypes = [
    ctypes.c_int32, ctypes.c_int32,
]
encodeig4.restype = ctypes.c_int32

decodeig4 = libgeoref.decode_cs_ig4
decodeig4.argtypes = [
    ctypes.c_int32, 
    ctypes.POINTER(ctypes.c_int32),
    ctypes.POINTER(ctypes.c_int32),
]
decodeig4.restype = ctypes.c_void_p