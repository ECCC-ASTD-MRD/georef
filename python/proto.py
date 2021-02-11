#!/usr/bin/env python3

import ctypes as ct
import numpy  as np
import numpy.ctypeslib as npc
from . import libgeoref
from .lib_structs import TGeoRef

# int GeoRef_GridGetParams(TGeoRef *Ref,int *NI,int *NJ,char *GRTYP,int *IG1,int *IG2,int *IG3,int *IG4,
#                           char *GRREF,int *IG1REF,int *IG2REF,int *IG3REF,int *IG4REF);
libgeoref.GeoRef_GridGetParams.argtypes = (ct.POINTER(TGeoRef),
                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                                ct.c_char_p,
                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                                ct.c_char_p,
                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int))
libgeoref.GeoRef_GridGetParams.restype = ct.c_int
GeoRef_GridGetParams = libgeoref.GeoRef_GridGetParams

# int GeoRef_Interp(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout,float *zin);
libgeoref.GeoRef_Interp.argtypes = (
    ct.POINTER(TGeoRef),
    ct.POINTER(TGeoRef),
    npc.ndpointer(dtype=np.float32),
    npc.ndpointer(dtype=np.float32)
    )
libgeoref.GeoRef_Interp.restype = ct.c_int
GeoRef_Interp = libgeoref.GeoRef_Interp

# TGeoRef* GeoRef_Create(int NI,int NJ,char *GRTYP,int IG1,int IG2,int IG3,int IG4,int FID);
libgeoref.GeoRef_Create.argtypes = (
    ct.c_int, ct.c_int, ct.c_char_p,
    ct.c_int, ct.c_int,
    ct.c_int, ct.c_int, ct.c_int
    )
libgeoref.GeoRef_Create.restype = ct.POINTER(TGeoRef)
GeoRef_Create = libgeoref.GeoRef_Create

# int GeoRef_GetLL(TGeoRef *Ref,double *Lat,double *Lon);
libgeoref.GeoRef_GetLL.argtypes = (
    ct.POINTER(TGeoRef),
    npc.ndpointer(dtype=np.float64),
    npc.ndpointer(dtype=np.float64)
    )
libgeoref.GeoRef_GetLL.restype  = ct.c_int
GeoRef_GetLL = libgeoref.GeoRef_GetLL

