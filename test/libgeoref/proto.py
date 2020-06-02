#!/usr/bin/env python3

import ctypes as ct
import numpy  as np
import numpy.ctypeslib as npc
from . import libgeoref

libgeoref.c_ezdefset.argtypes = (ct.c_int, ct.c_int)
libgeoref.c_ezdefset.restype = ct.c_int
c_ezdefset = libgeoref.c_ezdefset

libgeoref.c_ezget_nsubgrids.argtypes = (ct.c_int, )
libgeoref.c_ezget_nsubgrids.restype  = ct.c_int
c_ezget_nsubgrids = libgeoref.c_ezget_nsubgrids


libgeoref.c_ezget_subgridids.argtypes = (ct.c_int, npc.ndpointer(dtype=np.intc))
libgeoref.c_ezget_subgridids.restype  = ct.c_int
c_ezget_subgridids = libgeoref.c_ezget_subgridids

libgeoref.c_ezgprm.argtypes = (ct.c_int, ct.c_char_p,
                               ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                               ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                               ct.POINTER(ct.c_int), ct.POINTER(ct.c_int))
libgeoref.c_ezgprm.restype = ct.c_int
c_ezgprm = libgeoref.c_ezgprm

libgeoref.c_ezgxprm.argtypes = (ct.c_int, ct.POINTER(ct.c_int),
                                ct.POINTER(ct.c_int),
                                ct.c_char_p,
                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                                ct.c_char_p,
                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int))
libgeoref.c_ezgxprm.restype = ct.c_int
c_ezgxprm = libgeoref.c_ezgxprm

libgeoref.c_ezqkdef.argtypes = (
    ct.c_int, ct.c_int, ct.c_char_p,
    ct.c_int, ct.c_int,
    ct.c_int, ct.c_int, ct.c_int
    )
libgeoref.c_ezqkdef.restype = ct.c_int
c_ezqkdef = libgeoref.c_ezqkdef

libgeoref.c_ezsint.argtypes = (
    npc.ndpointer(dtype=np.float32),
    npc.ndpointer(dtype=np.float32)
    )
libgeoref.c_ezsint.restype = ct.c_int
c_ezsint = libgeoref.c_ezsint

libgeoref.c_gdll.argtypes = (
    ct.c_int,
    npc.ndpointer(dtype=np.float32),
    npc.ndpointer(dtype=np.float32)
    )
libgeoref.c_gdll.restype  = ct.c_int
c_gdll = libgeoref.c_gdll

