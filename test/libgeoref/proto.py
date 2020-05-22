#!/usr/bin/env python3

from . import libgeoref
import ctypes as ct
import numpy  as np
import numpy.ctypeslib as npc

libgeoref.c_ezsint.argtypes = (
    npc.ndpointer(dtype=np.float32),
    npc.ndpointer(dtype=np.float32)
    )
libgeoref.c_ezsint.restype  = ct.c_int
c_ezsint = libgeoref.c_ezsint

#print(libgeoref.c_ezsint)
#print(c_ezsint)