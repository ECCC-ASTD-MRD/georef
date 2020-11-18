#!/usr/bin/env python3

import ctypes as ct
import numpy  as np
import numpy.ctypeslib as npc
from . import libgeoref
from .lib_structs import TGeoRef#, Point

# # c_ezdefset
# libgeoref.c_ezdefset.argtypes = (ct.c_int, ct.c_int)
# libgeoref.c_ezdefset.restype = ct.c_int
# c_ezdefset = libgeoref.c_ezdefset

# # c_ezget_NbSub
# libgeoref.c_ezget_NbSub.argtypes = (ct.c_int, )
# libgeoref.c_ezget_NbSub.restype  = ct.c_int
# c_ezget_NbSub = libgeoref.c_ezget_NbSub

# # c_ezget_subgridids
# libgeoref.c_ezget_subgridids.argtypes = (ct.c_int, npc.ndpointer(dtype=np.intc))
# libgeoref.c_ezget_subgridids.restype  = ct.c_int
# c_ezget_subgridids = libgeoref.c_ezget_subgridids

# # c_ezgprm
# libgeoref.c_ezgprm.argtypes = (ct.c_int, ct.c_char_p,
#                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
#                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
#                                ct.POINTER(ct.c_int), ct.POINTER(ct.c_int))
# libgeoref.c_ezgprm.restype = ct.c_int
# c_ezgprm = libgeoref.c_ezgprm


# # c_ezgxprm
# libgeoref.c_ezgxprm.argtypes = (ct.c_int, ct.POINTER(ct.c_int),
#                                 ct.POINTER(ct.c_int),
#                                 ct.c_char_p,
#                                 ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
#                                 ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
#                                 ct.c_char_p,
#                                 ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
#                                 ct.POINTER(ct.c_int), ct.POINTER(ct.c_int))
# libgeoref.c_ezgxprm.restype = ct.c_int
# c_ezgxprm = libgeoref.c_ezgxprm

# # c_ezqkdef
# libgeoref.c_ezqkdef.argtypes = (
#     ct.c_int, ct.c_int, ct.c_char_p,
#     ct.c_int, ct.c_int,
#     ct.c_int, ct.c_int, ct.c_int
#     )
# libgeoref.c_ezqkdef.restype = ct.c_int
# c_ezqkdef = libgeoref.c_ezqkdef

# # c_gdll
# libgeoref.c_gdll.argtypes = (
#     ct.c_int,
#     npc.ndpointer(dtype=np.float32),
#     npc.ndpointer(dtype=np.float32)
#     )
# libgeoref.c_gdll.restype  = ct.c_int
# c_gdll = libgeoref.c_gdll

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
GeoRef_GridGetParams = libgeoref.c_ezgxprm

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
#TODO Maude: return pointer to TGeoRef
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

# TEST
# Point* Create(int lat, int lon)
# libgeoref.Create.argtypes = (
#     ct.c_int, ct.c_int
#     )
# libgeoref.Create.restype  = ct.POINTER(Point)
# Create = libgeoref.Create
