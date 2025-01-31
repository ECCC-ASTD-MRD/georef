#!/usr/bin/env python3

import ctypes as ct
import numpy  as np
import numpy.ctypeslib as npc
from . import libgeoref
from .lib_structs import TGeoRef, TDef, TGeoScan, TGeoOptions, TGeoSet, TZRef, TCoord, OGR_Layer, Vect2d, Vect3d, T3DArray, TZRefInterp, TGridSet   

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

# void Def_Free(TDef *Def)
libgeoref.Def_Free.argtypes = (ct.POINTER(TDef),)
libgeoref.Def_Free.restype = None
Def_Free = libgeoref.Def_Free

# void Def_Clear(TDef *Def)
libgeoref.Def_Clear.argtypes = (ct.POINTER(TDef),)
libgeoref.Def_Clear.restype = None
Def_Clear = libgeoref.Def_Clear

# TDef* Def_Copy(TDef *Def)
libgeoref.Def_Copy.argtypes = (ct.POINTER(TDef),)
libgeoref.Def_Copy.restype = ct.POINTER(TDef)
Def_Copy = libgeoref.Def_Copy

# TDef* Def_New(int32_t Type)
libgeoref.Def_New.argtypes = (ct.c_int32,)
libgeoref.Def_New.restype = ct.POINTER(TDef)
Def_New = libgeoref.Def_New

# int32_t Def_Resize(TDef *Def,int32_t NI,int32_t NJ,int32_t NK)
libgeoref.Def_Resize.argtypes = (ct.POINTER(TDef), ct.c_int32, ct.c_int32, ct.c_int32)
libgeoref.Def_Resize.restype = ct.c_int32
Def_Resize = libgeoref.Def_Resize

# int32_t Def_Define(TDef *Def,int32_t NI,int32_t NJ,int32_t NK,int32_t NbData)
libgeoref.Def_Define.argtypes = (ct.POINTER(TDef), ct.c_int32, ct.c_int32, ct.c_int32, ct.c_int32)
libgeoref.Def_Define.restype = ct.c_int32
Def_Define = libgeoref.Def_Define

# int32_t Def_Get(TDef *Def,int32_t Idx,double *Val)
libgeoref.Def_Get.argtypes = (ct.POINTER(TDef), ct.c_int32, ct.POINTER(ct.c_double))
libgeoref.Def_Get.restype = ct.c_int32
Def_Get = libgeoref.Def_Get

# int32_t Def_Set(TDef *Def,int32_t Idx,double Val)
libgeoref.Def_Set.argtypes = (ct.POINTER(TDef), ct.c_int32, ct.c_double)
libgeoref.Def_Set.restype = ct.c_int32
Def_Set = libgeoref.Def_Set

# int32_t Def_Inc(TDef *Def,int32_t Idx,double Val)
libgeoref.Def_Inc.argtypes = (ct.POINTER(TDef), ct.c_int32, ct.c_double)
libgeoref.Def_Inc.restype = ct.c_int32
Def_Inc = libgeoref.Def_Inc

# int32_t Def_Dec(TDef *Def,int32_t Idx,double Val)
libgeoref.Def_Dec.argtypes = (ct.POINTER(TDef), ct.c_int32, ct.c_double)
libgeoref.Def_Dec.restype = ct.c_int32
Def_Dec = libgeoref.Def_Dec

# int32_t Def_Mul(TDef *Def,int32_t Idx,double Val)
libgeoref.Def_Mul.argtypes = (ct.POINTER(TDef), ct.c_int32, ct.c_double)
libgeoref.Def_Mul.restype = ct.c_int32
Def_Mul = libgeoref.Def_Mul

# int32_t Def_Div(TDef *Def,int32_t Idx,double Val)
libgeoref.Def_Div.argtypes = (ct.POINTER(TDef), ct.c_int32, ct.c_double)
libgeoref.Def_Div.restype = ct.c_int32
Def_Div = libgeoref.Def_Div

# void Def_GetRange(TDef *Def,double *Min,double *Max)
libgeoref.Def_GetRange.argtypes = (ct.POINTER(TDef), 
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64))
libgeoref.Def_GetRange.restype = None
Def_GetRange = libgeoref.Def_GetRange

# void Def_GetLimits(TDef *Def,int32_t Dim,int32_t *Min,int32_t *Max)
libgeoref.Def_GetLimits.argtypes = (ct.POINTER(TDef), ct.c_int32, ct.POINTER(ct.c_int32), ct.POINTER(ct.c_int32))
libgeoref.Def_GetLimits.restype = None
Def_GetLimits = libgeoref.Def_GetLimits

# void Def_SetLimits(TDef *Def,int32_t Dim,int32_t Min,int32_t Max)
libgeoref.Def_SetLimits.argtypes = (ct.POINTER(TDef), ct.c_int32, ct.c_int32, ct.c_int32)
libgeoref.Def_SetLimits.restype = None
Def_SetLimits = libgeoref.Def_SetLimits

# void GeoScan_Init(TGeoScan *Scan)
libgeoref.GeoScan_Init.argtypes = (ct.POINTER(TGeoScan),)
libgeoref.GeoScan_Init.restype = None
GeoScan_Init = libgeoref.GeoScan_Init

# void GeoScan_Clear(TGeoScan *Scan)
libgeoref.GeoScan_Clear.argtypes = (ct.POINTER(TGeoScan),)
libgeoref.GeoScan_Clear.restype = None
GeoScan_Clear = libgeoref.GeoScan_Clear

# int32_t GeoScan_Get(TGeoScan *Scan,TGeoRef *ToRef,TDef *ToDef,TGeoRef *FromRef,TDef *FromDef,TGeoOptions *Opt,int32_t X0,int32_t Y0,int32_t X1,int32_t Y1,int32_t Dim)
libgeoref.GeoScan_Get.argtypes = (ct.POINTER(TGeoScan), ct.POINTER(TGeoRef), ct.POINTER(TDef),
                                 ct.POINTER(TGeoRef), ct.POINTER(TDef), ct.POINTER(TGeoOptions),
                                 ct.c_int32, ct.c_int32, ct.c_int32, ct.c_int32, ct.c_int32)
libgeoref.GeoScan_Get.restype = ct.c_int32
GeoScan_Get = libgeoref.GeoScan_Get


# lib/GeoRef_Axis.c
# void GeoRef_GridGetExpanded(TGeoRef *Ref,TGeoOptions *Opt,float *zout,float *zin)
libgeoref.GeoRef_GridGetExpanded.argtypes = (ct.POINTER(TGeoRef), ct.POINTER(TGeoOptions),
                                            npc.ndpointer(dtype=np.float32),
                                            npc.ndpointer(dtype=np.float32))
libgeoref.GeoRef_GridGetExpanded.restype = None
GeoRef_GridGetExpanded = libgeoref.GeoRef_GridGetExpanded

# void GeoRef_AxisDefine(TGeoRef* Ref,double *AX,double *AY)
libgeoref.GeoRef_AxisDefine.argtypes = (ct.POINTER(TGeoRef),
                                       npc.ndpointer(dtype=np.float64),
                                       npc.ndpointer(dtype=np.float64))
libgeoref.GeoRef_AxisDefine.restype = None
GeoRef_AxisDefine = libgeoref.GeoRef_AxisDefine

# int32_t GeoRef_AxisGet(TGeoRef *Ref,double *AX,double *AY)
libgeoref.GeoRef_AxisGet.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64))
libgeoref.GeoRef_AxisGet.restype = ct.c_int32
GeoRef_AxisGet = libgeoref.GeoRef_AxisGet

# int32_t GeoRef_AxisDefineYY(TGeoRef *Ref)
libgeoref.GeoRef_AxisDefineYY.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_AxisDefineYY.restype = ct.c_int32
GeoRef_AxisDefineYY = libgeoref.GeoRef_AxisDefineYY

# int32_t GeoRef_AxisDefineYYLookup(TGeoRef *Ref)
libgeoref.GeoRef_AxisDefineYYLookup.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_AxisDefineYYLookup.restype = ct.c_int32
GeoRef_AxisDefineYYLookup = libgeoref.GeoRef_AxisDefineYYLookup

# int32_t GeoRef_AxisDefineYYIndex(TGeoRef *Ref)
libgeoref.GeoRef_AxisDefineYYIndex.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_AxisDefineYYIndex.restype = ct.c_int32
GeoRef_AxisDefineYYIndex = libgeoref.GeoRef_AxisDefineYYIndex

# int32_t GeoRef_AxisDefineYYWrap(TGeoRef *Ref)
libgeoref.GeoRef_AxisDefineYYWrap.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_AxisDefineYYWrap.restype = ct.c_int32
GeoRef_AxisDefineYYWrap = libgeoref.GeoRef_AxisDefineYYWrap

# int32_t GeoRef_AxisDefineYYCell(TGeoRef *Ref)
libgeoref.GeoRef_AxisDefineYYCell.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_AxisDefineYYCell.restype = ct.c_int32
GeoRef_AxisDefineYYCell = libgeoref.GeoRef_AxisDefineYYCell

# lib/GeoRef_Correct.c
# int32_t GeoRef_CorrValNorth(TGeoSet *Set,float *zout,float *zin)
libgeoref.GeoRef_CorrValNorth.argtypes = (ct.POINTER(TGeoSet),
                                         npc.ndpointer(dtype=np.float32),
                                         npc.ndpointer(dtype=np.float32))
libgeoref.GeoRef_CorrValNorth.restype = ct.c_int32
GeoRef_CorrValNorth = libgeoref.GeoRef_CorrValNorth

# int32_t GeoRef_CorrValSouth(TGeoSet *Set,float *zout,float *zin)
libgeoref.GeoRef_CorrValSouth.argtypes = (ct.POINTER(TGeoSet),
                                         npc.ndpointer(dtype=np.float32),
                                         npc.ndpointer(dtype=np.float32))
libgeoref.GeoRef_CorrValSouth.restype = ct.c_int32
GeoRef_CorrValSouth = libgeoref.GeoRef_CorrValSouth

# int32_t GeoRef_CorrectValue(TGeoSet *Set,float *zout,float *zin)
libgeoref.GeoRef_CorrectValue.argtypes = (ct.POINTER(TGeoSet),
                                         npc.ndpointer(dtype=np.float32),
                                         npc.ndpointer(dtype=np.float32))
libgeoref.GeoRef_CorrectValue.restype = ct.c_int32
GeoRef_CorrectValue = libgeoref.GeoRef_CorrectValue

# int32_t GeoRef_CorrVecNorth(TGeoSet *Set,float *uuout,float *vvout,float *uuin,float *vvin)
libgeoref.GeoRef_CorrVecNorth.argtypes = (ct.POINTER(TGeoSet),
                                         npc.ndpointer(dtype=np.float32),
                                         npc.ndpointer(dtype=np.float32),
                                         npc.ndpointer(dtype=np.float32),
                                         npc.ndpointer(dtype=np.float32))
libgeoref.GeoRef_CorrVecNorth.restype = ct.c_int32
GeoRef_CorrVecNorth = libgeoref.GeoRef_CorrVecNorth

# int32_t GeoRef_CorrVecSouth(TGeoSet *Set,float *uuout,float *vvout,float *uuin,float *vvin)
libgeoref.GeoRef_CorrVecSouth.argtypes = (ct.POINTER(TGeoSet),
                                         npc.ndpointer(dtype=np.float32),
                                         npc.ndpointer(dtype=np.float32),
                                         npc.ndpointer(dtype=np.float32),
                                         npc.ndpointer(dtype=np.float32))
libgeoref.GeoRef_CorrVecSouth.restype = ct.c_int32
GeoRef_CorrVecSouth = libgeoref.GeoRef_CorrVecSouth
# int32_t GeoRef_CorrectVector(TGeoSet *Set,float *uuout,float *vvout,float *uuin,float *vvin)
libgeoref.GeoRef_CorrectVector.argtypes = (ct.POINTER(TGeoSet),
                                          npc.ndpointer(dtype=np.float32),
                                          npc.ndpointer(dtype=np.float32),
                                          npc.ndpointer(dtype=np.float32),
                                          npc.ndpointer(dtype=np.float32))
libgeoref.GeoRef_CorrectVector.restype = ct.c_int32
GeoRef_CorrectVector = libgeoref.GeoRef_CorrectVector

# lib/GeoRef_Func.c

# double GeoRef_Distance(TGeoRef *Ref,double X0,double Y0,double X1,double Y1)
libgeoref.GeoRef_Distance.argtypes = (ct.POINTER(TGeoRef), ct.c_double, ct.c_double, ct.c_double, ct.c_double)
libgeoref.GeoRef_Distance.restype = ct.c_double
GeoRef_Distance = libgeoref.GeoRef_Distance

# double GeoRef_Distance2(TGeoRef *Ref,double X0,double Y0,double X1,double Y1)
libgeoref.GeoRef_Distance2.argtypes = (ct.POINTER(TGeoRef), ct.c_double, ct.c_double, ct.c_double, ct.c_double)
libgeoref.GeoRef_Distance2.restype = ct.c_double
GeoRef_Distance2 = libgeoref.GeoRef_Distance2

# double GeoRef_HeightFromIP(TGeoRef *Ref,TZRef *ZRef,double Azimuth,double Bin,double Sweep)
libgeoref.GeoRef_HeightFromIP.argtypes = (ct.POINTER(TGeoRef), ct.POINTER(TZRef), ct.c_double, ct.c_double, ct.c_double)
libgeoref.GeoRef_HeightFromIP.restype = ct.c_double
GeoRef_HeightFromIP = libgeoref.GeoRef_HeightFromIP

# double GeoRef_HeightFromLevel(TGeoRef *Ref,TZRef *ZRef,double Azimuth,double Bin,double Level)
libgeoref.GeoRef_HeightFromLevel.argtypes = (ct.POINTER(TGeoRef), ct.POINTER(TZRef), ct.c_double, ct.c_double, ct.c_double)
libgeoref.GeoRef_HeightFromLevel.restype = ct.c_double
GeoRef_HeightFromLevel = libgeoref.GeoRef_HeightFromLevel

# int32_t GeoRef_RDR2XY(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_RDR2XY.argtypes = (ct.POINTER(TGeoRef),
                                   npc.ndpointer(dtype=np.float64),
                                   npc.ndpointer(dtype=np.float64),
                                   npc.ndpointer(dtype=np.float64),
                                   npc.ndpointer(dtype=np.float64),
                                   ct.c_int32)
libgeoref.GeoRef_RDR2XY.restype = ct.c_int32
GeoRef_RDR2XY = libgeoref.GeoRef_RDR2XY

# int32_t GeoRef_RDR2LL(TGeoRef *Ref,double *Lat,double *Lon,double Azimuth,double Range)
libgeoref.GeoRef_RDR2LL.argtypes = (ct.POINTER(TGeoRef),
                                   npc.ndpointer(dtype=np.float64),
                                   npc.ndpointer(dtype=np.float64),
                                   ct.c_double,
                                   ct.c_double)
libgeoref.GeoRef_RDR2LL.restype = ct.c_int32
GeoRef_RDR2LL = libgeoref.GeoRef_RDR2LL

# int32_t GeoRef_RDRDistance(TGeoRef *Ref,double *Dist,double *Azimuth,double *Elev,double Lat,double Lon,double Height)
libgeoref.GeoRef_RDRDistance.argtypes = (ct.POINTER(TGeoRef),
                                        npc.ndpointer(dtype=np.float64),
                                        npc.ndpointer(dtype=np.float64),
                                        npc.ndpointer(dtype=np.float64),
                                        ct.c_double,
                                        ct.c_double,
                                        ct.c_double)
libgeoref.GeoRef_RDRDistance.restype = ct.c_int32
GeoRef_RDRDistance = libgeoref.GeoRef_RDRDistance

# int32_t GeoRef_RDRPoint(TGeoRef *Ref,double *Lat,double *Lon,double *Height,double Range,double Azimuth,double Elev)
libgeoref.GeoRef_RDRPoint.argtypes = (ct.POINTER(TGeoRef),
                                     npc.ndpointer(dtype=np.float64),
                                     npc.ndpointer(dtype=np.float64),
                                     npc.ndpointer(dtype=np.float64),
                                     ct.c_double,
                                     ct.c_double,
                                     ct.c_double)
libgeoref.GeoRef_RDRPoint.restype = ct.c_int32
GeoRef_RDRPoint = libgeoref.GeoRef_RDRPoint

# int32_t GeoFunc_RadialIntersect(TCoord C1,TCoord C2,double CRS13,double CRS23,TCoord *C3)
libgeoref.GeoFunc_RadialIntersect.argtypes = (TCoord, TCoord, ct.c_double, ct.c_double, ct.POINTER(TCoord))
libgeoref.GeoFunc_RadialIntersect.restype = ct.c_int32
GeoFunc_RadialIntersect = libgeoref.GeoFunc_RadialIntersect

# lib/GeoRef_Interp.c

# int32_t c_gd_isgridrotated(TGeoRef *gr)
libgeoref.c_gd_isgridrotated.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.c_gd_isgridrotated.restype = ct.c_int32
c_gd_isgridrotated = libgeoref.c_gd_isgridrotated

# int32_t c_gdcompatible_grids(TGeoRef *RefFrom, TGeoRef* RefTo)
libgeoref.c_gdcompatible_grids.argtypes = (ct.POINTER(TGeoRef), ct.POINTER(TGeoRef))
libgeoref.c_gdcompatible_grids.restype = ct.c_int32
c_gdcompatible_grids = libgeoref.c_gdcompatible_grids

# int32_t GeoRef_Interp(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout,float *zin)
libgeoref.GeoRef_Interp.argtypes = (ct.POINTER(TGeoRef),
                                   ct.POINTER(TGeoRef),
                                   ct.POINTER(TGeoOptions),
                                   npc.ndpointer(dtype=np.float32),
                                   npc.ndpointer(dtype=np.float32))
libgeoref.GeoRef_Interp.restype = ct.c_int32
GeoRef_Interp = libgeoref.GeoRef_Interp

# int32_t GeoRef_InterpFinally(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout,float *zin,double *X,double *Y,int32_t npts,TGeoSet *GSet)
libgeoref.GeoRef_InterpFinally.argtypes = (ct.POINTER(TGeoRef),
                                          ct.POINTER(TGeoRef),
                                          ct.POINTER(TGeoOptions),
                                          npc.ndpointer(dtype=np.float32),
                                          npc.ndpointer(dtype=np.float32),
                                          npc.ndpointer(dtype=np.float64),
                                          npc.ndpointer(dtype=np.float64),
                                          ct.c_int32,
                                          ct.POINTER(TGeoSet))
libgeoref.GeoRef_InterpFinally.restype = ct.c_int32
GeoRef_InterpFinally = libgeoref.GeoRef_InterpFinally

# int32_t GeoRef_InterpYY(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout,float *zin)
libgeoref.GeoRef_InterpYY.argtypes = (ct.POINTER(TGeoRef),
                                     ct.POINTER(TGeoRef),
                                     ct.POINTER(TGeoOptions),
                                     npc.ndpointer(dtype=np.float32),
                                     npc.ndpointer(dtype=np.float32))
libgeoref.GeoRef_InterpYY.restype = ct.c_int32
GeoRef_InterpYY = libgeoref.GeoRef_InterpYY

# Skipping static functions:
# [static] int32_t GeoRef_InterpYYSub(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout,float *zin)
# [static] int32_t GeoRef_InterpYYSubMask(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout,float *zin)
# [static] int32_t GeoRef_InterpYYSub(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout,float *zin)
# [static] int32_t GeoRef_InterpYYSubMask(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout,float *zin)

# lib/GeoRef_InterpCoords.c

# void GeoRef_RotateXY(double *Lat,double *Lon,double *X,double *Y,int32_t npts,float xlat1,float xlon1,float xlat2,float xlon2)
libgeoref.GeoRef_RotateXY.argtypes = (npc.ndpointer(dtype=np.float64),
                                     npc.ndpointer(dtype=np.float64),
                                     npc.ndpointer(dtype=np.float64),
                                     npc.ndpointer(dtype=np.float64),
                                     ct.c_int32,
                                     ct.c_float, ct.c_float,
                                     ct.c_float, ct.c_float)
libgeoref.GeoRef_RotateXY.restype = None
GeoRef_RotateXY = libgeoref.GeoRef_RotateXY

# void GeoRef_RotateInvertXY(double *Lat,double *Lon,double *X,double *Y,int32_t npts,float xlat1,float xlon1,float xlat2,float xlon2)
libgeoref.GeoRef_RotateInvertXY.argtypes = (npc.ndpointer(dtype=np.float64),
                                           npc.ndpointer(dtype=np.float64),
                                           npc.ndpointer(dtype=np.float64),
                                           npc.ndpointer(dtype=np.float64),
                                           ct.c_int32,
                                           ct.c_float, ct.c_float,
                                           ct.c_float, ct.c_float)
libgeoref.GeoRef_RotateInvertXY.restype = None
GeoRef_RotateInvertXY = libgeoref.GeoRef_RotateInvertXY

# int32_t GeoRef_XY2LL(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb,int32_t Extrap)
libgeoref.GeoRef_XY2LL.argtypes = (ct.POINTER(TGeoRef),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  ct.c_int32, ct.c_int32)
libgeoref.GeoRef_XY2LL.restype = ct.c_int32
GeoRef_XY2LL = libgeoref.GeoRef_XY2LL

# int32_t GeoRef_LL2XY(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb,int32_t Extrap)
libgeoref.GeoRef_LL2XY.argtypes = (ct.POINTER(TGeoRef),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  ct.c_int32, ct.c_int32)
libgeoref.GeoRef_LL2XY.restype = ct.c_int32
GeoRef_LL2XY = libgeoref.GeoRef_LL2XY

# int32_t GeoRef_LL2GD(double *X,double *Y,double *Lat,double *Lon,int32_t Nb,double xlat0,double xlon0,double dellat,double dellon,double dgrw)
libgeoref.GeoRef_LL2GD.argtypes = (npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  ct.c_int32,
                                  ct.c_double, ct.c_double,
                                  ct.c_double, ct.c_double, ct.c_double)
libgeoref.GeoRef_LL2GD.restype = ct.c_int32
GeoRef_LL2GD = libgeoref.GeoRef_LL2GD

# int32_t GeoRef_GD2LL(double *Lat,double *Lon,double *X,double *Y,int32_t Nb,double xlat0,double xlon0,double dellat,double dellon,double dgrw)
libgeoref.GeoRef_GD2LL.argtypes = (npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  ct.c_int32,
                                  ct.c_double, ct.c_double,
                                  ct.c_double, ct.c_double, ct.c_double)
libgeoref.GeoRef_GD2LL.restype = ct.c_int32
GeoRef_GD2LL = libgeoref.GeoRef_GD2LL


# lib/GeoRef_InterpLL.c

# void GeoRef_GetLatLon(double *Lat,double *Lon,double D60,double DGRW,double PI,double PJ,int32_t HEM,double X,double Y)
libgeoref.GeoRef_GetLatLon.argtypes = (npc.ndpointer(dtype=np.float64),
                                      npc.ndpointer(dtype=np.float64),
                                      ct.c_double, ct.c_double,
                                      ct.c_double, ct.c_double,
                                      ct.c_int32,
                                      ct.c_double, ct.c_double)
libgeoref.GeoRef_GetLatLon.restype = None
GeoRef_GetLatLon = libgeoref.GeoRef_GetLatLon

# int32_t GeoRef_CalcLL(TGeoRef* Ref)
libgeoref.GeoRef_CalcLL.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_CalcLL.restype = ct.c_int32
GeoRef_CalcLL = libgeoref.GeoRef_CalcLL

# int32_t GeoRef_GetLL(TGeoRef *Ref,double *Lat,double *Lon)
libgeoref.GeoRef_GetLL.argtypes = (ct.POINTER(TGeoRef),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64))
libgeoref.GeoRef_GetLL.restype = ct.c_int32
GeoRef_GetLL = libgeoref.GeoRef_GetLL

# int32_t GeoRef_LLVal(TGeoRef *Ref,TGeoOptions *Opt,float *zout,float *zin,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LLVal.argtypes = (ct.POINTER(TGeoRef),
                                  ct.POINTER(TGeoOptions),
                                  npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  ct.c_int32)
libgeoref.GeoRef_LLVal.restype = ct.c_int32
GeoRef_LLVal = libgeoref.GeoRef_LLVal

# int32_t GeoRef_LLUVVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LLUVVal.argtypes = (ct.POINTER(TGeoRef),
                                    ct.POINTER(TGeoOptions),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_LLUVVal.restype = ct.c_int32
GeoRef_LLUVVal = libgeoref.GeoRef_LLUVVal

# int32_t GeoRef_LLWDVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LLWDVal.argtypes = (ct.POINTER(TGeoRef),
                                    ct.POINTER(TGeoOptions),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_LLWDVal.restype = ct.c_int32
GeoRef_LLWDVal = libgeoref.GeoRef_LLWDVal

# lib/GeoRef_InterpUV.c

# void c_ezgfwfllw8(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int32_t *ni,int32_t *nj,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4)
libgeoref.c_ezgfwfllw8.argtypes = (npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  ct.POINTER(ct.c_int32),
                                  ct.POINTER(ct.c_int32),
                                  ct.c_char_p,
                                  ct.POINTER(ct.c_int32),
                                  ct.POINTER(ct.c_int32),
                                  ct.POINTER(ct.c_int32),
                                  ct.POINTER(ct.c_int32))
libgeoref.c_ezgfwfllw8.restype = None
c_ezgfwfllw8 = libgeoref.c_ezgfwfllw8

# void c_ezllwfgfw8(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int32_t *ni,int32_t *nj,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4)
libgeoref.c_ezllwfgfw8.argtypes = (npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  ct.POINTER(ct.c_int32),
                                  ct.POINTER(ct.c_int32),
                                  ct.c_char_p,
                                  ct.POINTER(ct.c_int32),
                                  ct.POINTER(ct.c_int32),
                                  ct.POINTER(ct.c_int32),
                                  ct.POINTER(ct.c_int32))
libgeoref.c_ezllwfgfw8.restype = None
c_ezllwfgfw8 = libgeoref.c_ezllwfgfw8

# void c_ezllwfgff8(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int32_t *ni,int32_t *nj,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4)
libgeoref.c_ezllwfgff8.argtypes = (npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  ct.POINTER(ct.c_int32),
                                  ct.POINTER(ct.c_int32),
                                  ct.c_char_p,
                                  ct.POINTER(ct.c_int32),
                                  ct.POINTER(ct.c_int32),
                                  ct.POINTER(ct.c_int32),
                                  ct.POINTER(ct.c_int32))
libgeoref.c_ezllwfgff8.restype = None
c_ezllwfgff8 = libgeoref.c_ezllwfgff8

# int32_t GeoRef_InterpUV(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin)
libgeoref.GeoRef_InterpUV.argtypes = (ct.POINTER(TGeoRef),
                                     ct.POINTER(TGeoRef),
                                     ct.POINTER(TGeoOptions),
                                     npc.ndpointer(dtype=np.float32),
                                     npc.ndpointer(dtype=np.float32),
                                     npc.ndpointer(dtype=np.float32),
                                     npc.ndpointer(dtype=np.float32))
libgeoref.GeoRef_InterpUV.restype = ct.c_int32
GeoRef_InterpUV = libgeoref.GeoRef_InterpUV

# int32_t GeoRef_InterpWD(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin)
libgeoref.GeoRef_InterpWD.argtypes = (ct.POINTER(TGeoRef),
                                     ct.POINTER(TGeoRef),
                                     ct.POINTER(TGeoOptions),
                                     npc.ndpointer(dtype=np.float32),
                                     npc.ndpointer(dtype=np.float32),
                                     npc.ndpointer(dtype=np.float32),
                                     npc.ndpointer(dtype=np.float32))
libgeoref.GeoRef_InterpWD.restype = ct.c_int32
GeoRef_InterpWD = libgeoref.GeoRef_InterpWD

# int32_t GeoRef_InterpYYUV(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin)
libgeoref.GeoRef_InterpYYUV.argtypes = (ct.POINTER(TGeoRef),
                                       ct.POINTER(TGeoRef),
                                       ct.POINTER(TGeoOptions),
                                       npc.ndpointer(dtype=np.float32),
                                       npc.ndpointer(dtype=np.float32),
                                       npc.ndpointer(dtype=np.float32),
                                       npc.ndpointer(dtype=np.float32))
libgeoref.GeoRef_InterpYYUV.restype = ct.c_int32
GeoRef_InterpYYUV = libgeoref.GeoRef_InterpYYUV

# int32_t GeoRef_UV2WD(TGeoRef *Ref,float *spd_out,float *wd_out,float *uuin,float *vvin,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_UV2WD.argtypes = (ct.POINTER(TGeoRef),
                                  npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  ct.c_int32)
libgeoref.GeoRef_UV2WD.restype = ct.c_int32
GeoRef_UV2WD = libgeoref.GeoRef_UV2WD

# lib/GeoRef_InterpXY.c

# int32_t GeoRef_XYInterp(TGeoRef *Ref,TGeoOptions *Opt,float *zout,float *zin,double *X,double *Y,int32_t npts)
libgeoref.GeoRef_XYInterp.argtypes = (ct.POINTER(TGeoRef),
                                     ct.POINTER(TGeoOptions),
                                     npc.ndpointer(dtype=np.float32),
                                     npc.ndpointer(dtype=np.float32),
                                     npc.ndpointer(dtype=np.float64),
                                     npc.ndpointer(dtype=np.float64),
                                     ct.c_int32)
libgeoref.GeoRef_XYInterp.restype = ct.c_int32
GeoRef_XYInterp = libgeoref.GeoRef_XYInterp

# int32_t GeoRef_XYVal(TGeoRef *Ref,TGeoOptions *Opt,float *zout,float *zin,double *X,double *Y,int32_t n)
libgeoref.GeoRef_XYVal.argtypes = (ct.POINTER(TGeoRef),
                                  ct.POINTER(TGeoOptions),
                                  npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float32),
                                  npc.ndpointer(dtype=np.float64),
                                  npc.ndpointer(dtype=np.float64),
                                  ct.c_int32)
libgeoref.GeoRef_XYVal.restype = ct.c_int32
GeoRef_XYVal = libgeoref.GeoRef_XYVal

# int32_t GeoRef_XYUVVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *X,double *Y,int32_t n)
libgeoref.GeoRef_XYUVVal.argtypes = (ct.POINTER(TGeoRef),
                                    ct.POINTER(TGeoOptions),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_XYUVVal.restype = ct.c_int32
GeoRef_XYUVVal = libgeoref.GeoRef_XYUVVal

# int32_t GeoRef_XYWDVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *X,double *Y,int32_t n)
libgeoref.GeoRef_XYWDVal.argtypes = (ct.POINTER(TGeoRef),
                                    ct.POINTER(TGeoOptions),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float32),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_XYWDVal.restype = ct.c_int32
GeoRef_XYWDVal = libgeoref.GeoRef_XYWDVal

# lib/GeoRef_Mask.c

# int32_t GeoRef_InterpMask(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,char *MaskOut,char *MaskIn)
libgeoref.GeoRef_InterpMask.argtypes = (ct.POINTER(TGeoRef),
                                       ct.POINTER(TGeoRef),
                                       ct.POINTER(TGeoOptions),
                                       ct.c_char_p,
                                       ct.c_char_p)
libgeoref.GeoRef_InterpMask.restype = ct.c_int32
GeoRef_InterpMask = libgeoref.GeoRef_InterpMask

# int32_t GeoRef_MaskYYDefine(TGeoRef *Ref)
libgeoref.GeoRef_MaskYYDefine.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_MaskYYDefine.restype = ct.c_int32
GeoRef_MaskYYDefine = libgeoref.GeoRef_MaskYYDefine
# Skipping static function:
# [static] int32_t GeoRef_MaskYYSplit(double *dlat,double *dlon,int32_t *yyincount,int32_t *yyancount,int32_t npts)
# [static] int32_t GeoRef_MaskYYSplit(double *dlat,double *dlon,int32_t *yyincount,int32_t *yyancount,int32_t npts)


# lib/GeoRef_Set.c

# void GeoRef_SetZoneFree(TGeoSet *GSet)
libgeoref.GeoRef_SetZoneFree.argtypes = (ct.POINTER(TGeoSet),)
libgeoref.GeoRef_SetZoneFree.restype = None
GeoRef_SetZoneFree = libgeoref.GeoRef_SetZoneFree

# void GeoRef_SetFree(TGeoSet *GSet)
libgeoref.GeoRef_SetFree.argtypes = (ct.POINTER(TGeoSet),)
libgeoref.GeoRef_SetFree.restype = None
GeoRef_SetFree = libgeoref.GeoRef_SetFree

# int32_t GeoRef_SetZoneDefine(TGeoSet *GSet)
libgeoref.GeoRef_SetZoneDefine.argtypes = (ct.POINTER(TGeoSet),)
libgeoref.GeoRef_SetZoneDefine.restype = ct.c_int32
GeoRef_SetZoneDefine = libgeoref.GeoRef_SetZoneDefine

# int32_t GeoRef_SetZoneDefineOut(TGeoSet *GSet,int32_t Zone)
libgeoref.GeoRef_SetZoneDefineOut.argtypes = (ct.POINTER(TGeoSet), ct.c_int32)
libgeoref.GeoRef_SetZoneDefineOut.restype = ct.c_int32
GeoRef_SetZoneDefineOut = libgeoref.GeoRef_SetZoneDefineOut

# int32_t GeoRef_SetIndexInit(TGeoSet *GSet)
libgeoref.GeoRef_SetIndexInit.argtypes = (ct.POINTER(TGeoSet),)
libgeoref.GeoRef_SetIndexInit.restype = ct.c_int32
GeoRef_SetIndexInit = libgeoref.GeoRef_SetIndexInit

# int32_t GeoRef_SetCalcYYXY(TGeoSet *GSet)
libgeoref.GeoRef_SetCalcYYXY.argtypes = (ct.POINTER(TGeoSet),)
libgeoref.GeoRef_SetCalcYYXY.restype = ct.c_int32
GeoRef_SetCalcYYXY = libgeoref.GeoRef_SetCalcYYXY

# int32_t GeoRef_SetCalcYYLL(TGeoSet *GSet)
libgeoref.GeoRef_SetCalcYYLL.argtypes = (ct.POINTER(TGeoSet),)
libgeoref.GeoRef_SetCalcYYLL.restype = ct.c_int32
GeoRef_SetCalcYYLL = libgeoref.GeoRef_SetCalcYYLL

# TGeoSet* GeoRef_SetGet(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt)
libgeoref.GeoRef_SetGet.argtypes = (ct.POINTER(TGeoRef), ct.POINTER(TGeoRef), ct.POINTER(TGeoOptions))
libgeoref.GeoRef_SetGet.restype = ct.POINTER(TGeoSet)
GeoRef_SetGet = libgeoref.GeoRef_SetGet

# TGeoSet* GeoRef_SetRead(TGeoRef* RefTo,TGeoRef* RefFrom,int32_t InterpType,fst_file *File)
libgeoref.GeoRef_SetRead.argtypes = (ct.POINTER(TGeoRef), ct.POINTER(TGeoRef), ct.c_int32, ct.POINTER(fst_file))
libgeoref.GeoRef_SetRead.restype = ct.POINTER(TGeoSet)
GeoRef_SetRead = libgeoref.GeoRef_SetRead

# lib/GeoRef_Type_ABG.c

# int32_t GeoRef_XY2LL_ABG(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb)
libgeoref.GeoRef_XY2LL_ABG.argtypes = (ct.POINTER(TGeoRef),
                                      npc.ndpointer(dtype=np.float64),
                                      npc.ndpointer(dtype=np.float64),
                                      npc.ndpointer(dtype=np.float64),
                                      npc.ndpointer(dtype=np.float64),
                                      ct.c_int32)
libgeoref.GeoRef_XY2LL_ABG.restype = ct.c_int32
GeoRef_XY2LL_ABG = libgeoref.GeoRef_XY2LL_ABG

# int32_t GeoRef_LL2XY_ABG(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LL2XY_ABG.argtypes = (ct.POINTER(TGeoRef),
                                      npc.ndpointer(dtype=np.float64),
                                      npc.ndpointer(dtype=np.float64),
                                      npc.ndpointer(dtype=np.float64),
                                      npc.ndpointer(dtype=np.float64),
                                      ct.c_int32)
libgeoref.GeoRef_LL2XY_ABG.restype = ct.c_int32
GeoRef_LL2XY_ABG = libgeoref.GeoRef_LL2XY_ABG

# lib/GeoRef_Type_E.c

# int32_t GeoRef_XY2LL_E(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb)
libgeoref.GeoRef_XY2LL_E.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_XY2LL_E.restype = ct.c_int32
GeoRef_XY2LL_E = libgeoref.GeoRef_XY2LL_E

# int32_t GeoRef_LL2XY_E(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LL2XY_E.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_LL2XY_E.restype = ct.c_int32
GeoRef_LL2XY_E = libgeoref.GeoRef_LL2XY_E

# lib/GeoRef_Type_L.c

# int32_t GeoRef_XY2LL_L(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb)
libgeoref.GeoRef_XY2LL_L.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_XY2LL_L.restype = ct.c_int32
GeoRef_XY2LL_L = libgeoref.GeoRef_XY2LL_L

# int32_t GeoRef_LL2XY_L(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LL2XY_L.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_LL2XY_L.restype = ct.c_int32
GeoRef_LL2XY_L = libgeoref.GeoRef_LL2XY_L

# lib/GeoRef_Type_LAMBERT.c

# int32_t GeoRef_XY2LL_LAMBERT(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb)
libgeoref.GeoRef_XY2LL_LAMBERT.argtypes = (ct.POINTER(TGeoRef),
                                          npc.ndpointer(dtype=np.float64),
                                          npc.ndpointer(dtype=np.float64),
                                          npc.ndpointer(dtype=np.float64),
                                          npc.ndpointer(dtype=np.float64),
                                          ct.c_int32)
libgeoref.GeoRef_XY2LL_LAMBERT.restype = ct.c_int32
GeoRef_XY2LL_LAMBERT = libgeoref.GeoRef_XY2LL_LAMBERT

# int32_t GeoRef_LL2XY_LAMBERT(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LL2XY_LAMBERT.argtypes = (ct.POINTER(TGeoRef),
                                          npc.ndpointer(dtype=np.float64),
                                          npc.ndpointer(dtype=np.float64),
                                          npc.ndpointer(dtype=np.float64),
                                          npc.ndpointer(dtype=np.float64),
                                          ct.c_int32)
libgeoref.GeoRef_LL2XY_LAMBERT.restype = ct.c_int32
GeoRef_LL2XY_LAMBERT = libgeoref.GeoRef_LL2XY_LAMBERT

# lib/GeoRef_Type_M.c

# int32_t GeoRef_XY2LL_M(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb)
libgeoref.GeoRef_XY2LL_M.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_XY2LL_M.restype = ct.c_int32
GeoRef_XY2LL_M = libgeoref.GeoRef_XY2LL_M

# int32_t GeoRef_LL2XY_M(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LL2XY_M.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_LL2XY_M.restype = ct.c_int32
GeoRef_LL2XY_M = libgeoref.GeoRef_LL2XY_M

# lib/GeoRef_Type_NS.c

# int32_t GeoRef_XY2LL_NS(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb)
libgeoref.GeoRef_XY2LL_NS.argtypes = (ct.POINTER(TGeoRef),
                                     npc.ndpointer(dtype=np.float64),
                                     npc.ndpointer(dtype=np.float64),
                                     npc.ndpointer(dtype=np.float64),
                                     npc.ndpointer(dtype=np.float64),
                                     ct.c_int32)
libgeoref.GeoRef_XY2LL_NS.restype = ct.c_int32
GeoRef_XY2LL_NS = libgeoref.GeoRef_XY2LL_NS

# int32_t GeoRef_LL2XY_NS(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LL2XY_NS.argtypes = (ct.POINTER(TGeoRef),
                                     npc.ndpointer(dtype=np.float64),
                                     npc.ndpointer(dtype=np.float64),
                                     npc.ndpointer(dtype=np.float64),
                                     npc.ndpointer(dtype=np.float64),
                                     ct.c_int32)
libgeoref.GeoRef_LL2XY_NS.restype = ct.c_int32
GeoRef_LL2XY_NS = libgeoref.GeoRef_LL2XY_NS

# lib/GeoRef_Type_O.c

# int32_t GeoRef_XY2LL_O(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb)
libgeoref.GeoRef_XY2LL_O.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_XY2LL_O.restype = ct.c_int32
GeoRef_XY2LL_O = libgeoref.GeoRef_XY2LL_O

# int32_t GeoRef_LL2XY_O(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LL2XY_O.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_LL2XY_O.restype = ct.c_int32
GeoRef_LL2XY_O = libgeoref.GeoRef_LL2XY_O

# Skipping static function:
# [static] int32_t GeoRef_LocateO(TGeoRef *Ref,double Lat,double Lon,double *X,double *Y)

# lib/GeoRef_Type_R.c

# int32_t GeoRef_XY2LL_R(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb)
libgeoref.GeoRef_XY2LL_R.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_XY2LL_R.restype = ct.c_int32
GeoRef_XY2LL_R = libgeoref.GeoRef_XY2LL_R

# int32_t GeoRef_LL2XY_R(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LL2XY_R.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_LL2XY_R.restype = ct.c_int32
GeoRef_LL2XY_R = libgeoref.GeoRef_LL2XY_R

# lib/GeoRef_Type_T.c

# int32_t GeoRef_XY2LL_T(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb)
libgeoref.GeoRef_XY2LL_T.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_XY2LL_T.restype = ct.c_int32
GeoRef_XY2LL_T = libgeoref.GeoRef_XY2LL_T

# int32_t GeoRef_LL2XY_T(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LL2XY_T.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_LL2XY_T.restype = ct.c_int32
GeoRef_LL2XY_T = libgeoref.GeoRef_LL2XY_T

# lib/GeoRef_Type_U.c

# TGeoRef* GeoRef_U2W(TGeoRef *Ref)
libgeoref.GeoRef_U2W.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_U2W.restype = ct.POINTER(TGeoRef)
GeoRef_U2W = libgeoref.GeoRef_U2W

# TGeoRef* GeoRef_CreateU(int32_t NI,int32_t NJ,char *grtyp,char *grref,int32_t IG1,int32_t IG2,int32_t IG3,int32_t IG4,int32_t IG1ref,int32_t IG2ref,int32_t IG3ref,int32_t IG4ref,int32_t NbSub,TGeoRef **Subs)
libgeoref.GeoRef_CreateU.argtypes = (ct.c_int32, ct.c_int32,
                                    ct.c_char_p, ct.c_char_p,
                                    ct.c_int32, ct.c_int32, ct.c_int32, ct.c_int32,
                                    ct.c_int32, ct.c_int32, ct.c_int32, ct.c_int32,
                                    ct.c_int32, ct.POINTER(ct.POINTER(TGeoRef)))
libgeoref.GeoRef_CreateU.restype = ct.POINTER(TGeoRef)
GeoRef_CreateU = libgeoref.GeoRef_CreateU

# lib/GeoRef_Type_Y.c

# int32_t GeoRef_XY2LL_Y(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb)
libgeoref.GeoRef_XY2LL_Y.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_XY2LL_Y.restype = ct.c_int32
GeoRef_XY2LL_Y = libgeoref.GeoRef_XY2LL_Y

# int32_t GeoRef_LL2XY_Y(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LL2XY_Y.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_LL2XY_Y.restype = ct.c_int32
GeoRef_LL2XY_Y = libgeoref.GeoRef_LL2XY_Y

# int32_t GeoRef_Locate_Y(TGeoRef *Ref,double Lat,double Lon,int32_t *idx,double *dist,int32_t nb)
libgeoref.GeoRef_Locate_Y.argtypes = (ct.POINTER(TGeoRef),
                                     ct.c_double, ct.c_double,
                                     ct.POINTER(ct.c_int32),
                                     ct.POINTER(ct.c_double),
                                     ct.c_int32)
libgeoref.GeoRef_Locate_Y.restype = ct.c_int32
GeoRef_Locate_Y = libgeoref.GeoRef_Locate_Y

# lib/GeoRef_Type_Z.c

# int32_t GeoRef_XY2LL_Z(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb)
libgeoref.GeoRef_XY2LL_Z.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_XY2LL_Z.restype = ct.c_int32
GeoRef_XY2LL_Z = libgeoref.GeoRef_XY2LL_Z

# int32_t GeoRef_LL2XY_Z(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb)
libgeoref.GeoRef_LL2XY_Z.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_LL2XY_Z.restype = ct.c_int32
GeoRef_LL2XY_Z = libgeoref.GeoRef_LL2XY_Z

# TGeoRef* GeoRef_Z2Grid(TGeoRef *Ref,int32_t NI,int32_t NJ,double X0,double X1,double Y0,double Y1)
libgeoref.GeoRef_Z2Grid.argtypes = (ct.POINTER(TGeoRef),
                                   ct.c_int32, ct.c_int32,
                                   ct.c_double, ct.c_double,
                                   ct.c_double, ct.c_double)
libgeoref.GeoRef_Z2Grid.restype = ct.POINTER(TGeoRef)
GeoRef_Z2Grid = libgeoref.GeoRef_Z2Grid

# void GEM_hgrid4(double *AX,double *AY,int32_t NI,int32_t NJ,double *DX,double *DY,double x0,double x1,double y0,double y1,int32_t Wrap)
libgeoref.GEM_hgrid4.argtypes = (npc.ndpointer(dtype=np.float64),
                                npc.ndpointer(dtype=np.float64),
                                ct.c_int32, ct.c_int32,
                                npc.ndpointer(dtype=np.float64),
                                npc.ndpointer(dtype=np.float64),
                                ct.c_double, ct.c_double,
                                ct.c_double, ct.c_double,
                                ct.c_int32)
libgeoref.GEM_hgrid4.restype = None
GEM_hgrid4 = libgeoref.GEM_hgrid4

# lib/GeoRef_Utils.c

# double InterpCubic(double X0,double X1,double X2,double X3,double F)
libgeoref.InterpCubic.argtypes = (ct.c_double, ct.c_double, ct.c_double, ct.c_double, ct.c_double)
libgeoref.InterpCubic.restype = ct.c_double
InterpCubic = libgeoref.InterpCubic

# double InterpHermite(double X0,double X1,double X2,double X3,double F,double T,double B)
libgeoref.InterpHermite.argtypes = (ct.c_double, ct.c_double, ct.c_double, ct.c_double,
                                   ct.c_double, ct.c_double, ct.c_double)
libgeoref.InterpHermite.restype = ct.c_double
InterpHermite = libgeoref.InterpHermite

# int32_t QSort_StrPtr(const void *A,const void *B)
libgeoref.QSort_StrPtr.argtypes = (ct.c_void_p, ct.c_void_p)
libgeoref.QSort_StrPtr.restype = ct.c_int32
QSort_StrPtr = libgeoref.QSort_StrPtr

# void Unique(void *Arr,int* restrict Size,size_t NBytes)
libgeoref.Unique.argtypes = (ct.c_void_p, ct.POINTER(ct.c_int), ct.c_size_t)
libgeoref.Unique.restype = None
Unique = libgeoref.Unique

# double HCentile(double *M,int32_t N,int32_t K)
libgeoref.HCentile.argtypes = (npc.ndpointer(dtype=np.float64), ct.c_int32, ct.c_int32)
libgeoref.HCentile.restype = ct.c_double
HCentile = libgeoref.HCentile

# lib/GeoRef.c

# TGeoRef* GeoRef_New()
libgeoref.GeoRef_New.argtypes = None
libgeoref.GeoRef_New.restype = ct.POINTER(TGeoRef)
GeoRef_New = libgeoref.GeoRef_New

# void GeoRef_Free(TGeoRef *Ref)
libgeoref.GeoRef_Free.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_Free.restype = None
GeoRef_Free = libgeoref.GeoRef_Free

# void GeoRef_Clear(TGeoRef *Ref)
libgeoref.GeoRef_Clear.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_Clear.restype = None
GeoRef_Clear = libgeoref.GeoRef_Clear

# TGeoRef* GeoRef_HardCopy(TGeoRef* __restrict const Ref)
libgeoref.GeoRef_HardCopy.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_HardCopy.restype = ct.POINTER(TGeoRef)
GeoRef_HardCopy = libgeoref.GeoRef_HardCopy

# TGeoRef* GeoRef_Copy(TGeoRef *Ref)
libgeoref.GeoRef_Copy.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_Copy.restype = ct.POINTER(TGeoRef)
GeoRef_Copy = libgeoref.GeoRef_Copy

# TGeoRef* GeoRef_RPNSetup(int32_t NI,int32_t NJ,char *GRTYP,int32_t IG1,int32_t IG2,int32_t IG3,int32_t IG4)
libgeoref.GeoRef_RPNSetup.argtypes = (ct.c_int32, ct.c_int32, ct.c_char_p,
                                     ct.c_int32, ct.c_int32, ct.c_int32, ct.c_int32)
libgeoref.GeoRef_RPNSetup.restype = ct.POINTER(TGeoRef)
GeoRef_RPNSetup = libgeoref.GeoRef_RPNSetup

# TGeoRef* GeoRef_Setup(int32_t NI,int32_t NJ,char *GRTYP,int32_t IG1,int32_t IG2,int32_t IG3,int32_t IG4)
libgeoref.GeoRef_Setup.argtypes = (ct.c_int32, ct.c_int32, ct.c_char_p,
                                  ct.c_int32, ct.c_int32, ct.c_int32, ct.c_int32)
libgeoref.GeoRef_Setup.restype = ct.POINTER(TGeoRef)
GeoRef_Setup = libgeoref.GeoRef_Setup

# TGeoRef* GeoRef_Define(TGeoRef *Ref,int32_t NI,int32_t NJ,char *GRTYP,char *GRREF,int32_t IG1,int32_t IG2,int32_t IG3,int32_t IG4,double *AX,double *AY)
libgeoref.GeoRef_Define.argtypes = (ct.POINTER(TGeoRef),
                                   ct.c_int32, ct.c_int32,
                                   ct.c_char_p, ct.c_char_p,
                                   ct.c_int32, ct.c_int32, ct.c_int32, ct.c_int32,
                                   npc.ndpointer(dtype=np.float64),
                                   npc.ndpointer(dtype=np.float64))
libgeoref.GeoRef_Define.restype = ct.POINTER(TGeoRef)
GeoRef_Define = libgeoref.GeoRef_Define

# lib/GeoRef.c (continued)

# int32_t GeoRef_Size(TGeoRef *Ref,int32_t *Size)
libgeoref.GeoRef_Size.argtypes = (ct.POINTER(TGeoRef), ct.POINTER(ct.c_int32))
libgeoref.GeoRef_Size.restype = ct.c_int32
GeoRef_Size = libgeoref.GeoRef_Size

# int32_t GeoRef_Resize(TGeoRef *Ref,int32_t NI,int32_t NJ)
libgeoref.GeoRef_Resize.argtypes = (ct.POINTER(TGeoRef), ct.c_int32, ct.c_int32)
libgeoref.GeoRef_Resize.restype = ct.c_int32
GeoRef_Resize = libgeoref.GeoRef_Resize

# int32_t GeoRef_Grid(TGeoRef *Ref,double *AX,double *AY)
libgeoref.GeoRef_Grid.argtypes = (ct.POINTER(TGeoRef),
                                 npc.ndpointer(dtype=np.float64),
                                 npc.ndpointer(dtype=np.float64))
libgeoref.GeoRef_Grid.restype = ct.c_int32
GeoRef_Grid = libgeoref.GeoRef_Grid

# int32_t GeoRef_Qualify(TGeoRef *Ref)
libgeoref.GeoRef_Qualify.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_Qualify.restype = ct.c_int32
GeoRef_Qualify = libgeoref.GeoRef_Qualify

# int32_t GeoRef_Equal(TGeoRef *RefA,TGeoRef *RefB)
libgeoref.GeoRef_Equal.argtypes = (ct.POINTER(TGeoRef), ct.POINTER(TGeoRef))
libgeoref.GeoRef_Equal.restype = ct.c_int32
GeoRef_Equal = libgeoref.GeoRef_Equal

# int32_t GeoRef_Nearest(TGeoRef* __restrict const Ref,double X,double Y,int32_t *Idxs,double *Dists,int32_t NbNear,double MaxDist)
libgeoref.GeoRef_Nearest.argtypes = (ct.POINTER(TGeoRef),
                                    ct.c_double, ct.c_double,
                                    ct.POINTER(ct.c_int32),
                                    ct.POINTER(ct.c_double),
                                    ct.c_int32, ct.c_double)
libgeoref.GeoRef_Nearest.restype = ct.c_int32
GeoRef_Nearest = libgeoref.GeoRef_Nearest

# int32_t GeoRef_Within(TGeoRef *Ref,double X,double Y)
libgeoref.GeoRef_Within.argtypes = (ct.POINTER(TGeoRef), ct.c_double, ct.c_double)
libgeoref.GeoRef_Within.restype = ct.c_int32
GeoRef_Within = libgeoref.GeoRef_Within

# int32_t GeoRef_Value(TGeoRef *Ref,double X,double Y)
libgeoref.GeoRef_Value.argtypes = (ct.POINTER(TGeoRef), ct.c_double, ct.c_double)
libgeoref.GeoRef_Value.restype = ct.c_int32
GeoRef_Value = libgeoref.GeoRef_Value

# int32_t GeoRef_Limits(TGeoRef *Ref,double *X0,double *Y0,double *X1,double *Y1)
libgeoref.GeoRef_Limits.argtypes = (ct.POINTER(TGeoRef),
                                   ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
                                   ct.POINTER(ct.c_double), ct.POINTER(ct.c_double))
libgeoref.GeoRef_Limits.restype = ct.c_int32
GeoRef_Limits = libgeoref.GeoRef_Limits

# int32_t GeoRef_CellDims(TGeoRef *Ref,int32_t Invert,float *DX,float *DY,float *DA)
libgeoref.GeoRef_CellDims.argtypes = (ct.POINTER(TGeoRef), ct.c_int32,
                                     ct.POINTER(ct.c_float),
                                     ct.POINTER(ct.c_float),
                                     ct.POINTER(ct.c_float))
libgeoref.GeoRef_CellDims.restype = ct.c_int32
GeoRef_CellDims = libgeoref.GeoRef_CellDims

# int32_t GeoRef_UnProject(TGeoRef *Ref,double *X,double *Y,int32_t N)
libgeoref.GeoRef_UnProject.argtypes = (ct.POINTER(TGeoRef),
                                      npc.ndpointer(dtype=np.float64),
                                      npc.ndpointer(dtype=np.float64),
                                      ct.c_int32)
libgeoref.GeoRef_UnProject.restype = ct.c_int32
GeoRef_UnProject = libgeoref.GeoRef_UnProject

# int32_t GeoRef_Project(TGeoRef *Ref,double *X,double *Y,int32_t N)
libgeoref.GeoRef_Project.argtypes = (ct.POINTER(TGeoRef),
                                    npc.ndpointer(dtype=np.float64),
                                    npc.ndpointer(dtype=np.float64),
                                    ct.c_int32)
libgeoref.GeoRef_Project.restype = ct.c_int32
GeoRef_Project = libgeoref.GeoRef_Project

# int32_t GeoRef_RDR2RLL(TGeoRef *Ref,double *Range,double *Lat,double *Lon,double Azimuth,double Elev)
libgeoref.GeoRef_RDR2RLL.argtypes = (ct.POINTER(TGeoRef),
                                    ct.POINTER(ct.c_double),
                                    ct.POINTER(ct.c_double),
                                    ct.POINTER(ct.c_double),
                                    ct.c_double, ct.c_double)
libgeoref.GeoRef_RDR2RLL.restype = ct.c_int32
GeoRef_RDR2RLL = libgeoref.GeoRef_RDR2RLL

# int32_t GeoRef_RLL2RDR(TGeoRef *Ref,double *Range,double *Azimuth,double *Elev,double Lat,double Lon,double Height)
libgeoref.GeoRef_RLL2RDR.argtypes = (ct.POINTER(TGeoRef),
                                    ct.POINTER(ct.c_double),
                                    ct.POINTER(ct.c_double),
                                    ct.POINTER(ct.c_double),
                                    ct.c_double, ct.c_double, ct.c_double)
libgeoref.GeoRef_RLL2RDR.restype = ct.c_int32
GeoRef_RLL2RDR = libgeoref.GeoRef_RLL2RDR
# lib/GeoRef.c (continued)

# int32_t GeoRef_DefRPNXG(TGeoRef* Ref)
libgeoref.GeoRef_DefRPNXG.argtypes = (ct.POINTER(TGeoRef),)
libgeoref.GeoRef_DefRPNXG.restype = ct.c_int32
GeoRef_DefRPNXG = libgeoref.GeoRef_DefRPNXG

# int32_t GeoRef_RPNDesc(TGeoRef *GRef,int32_t *NI,int32_t *NJ,char *GRTYP,int32_t *IG1,int32_t *IG2,int32_t *IG3,int32_t *IG4,char *GRREF,int32_t *IG1REF,int32_t *IG2REF,int32_t *IG3REF,int32_t *IG4REF)
libgeoref.GeoRef_RPNDesc.argtypes = (ct.POINTER(TGeoRef),
                                    ct.POINTER(ct.c_int32), ct.POINTER(ct.c_int32),
                                    ct.c_char_p,
                                    ct.POINTER(ct.c_int32), ct.POINTER(ct.c_int32),
                                    ct.POINTER(ct.c_int32), ct.POINTER(ct.c_int32),
                                    ct.c_char_p,
                                    ct.POINTER(ct.c_int32), ct.POINTER(ct.c_int32),
                                    ct.POINTER(ct.c_int32), ct.POINTER(ct.c_int32))
libgeoref.GeoRef_RPNDesc.restype = ct.c_int32
GeoRef_RPNDesc = libgeoref.GeoRef_RPNDesc

# int32_t GeoRef_GetField(TGeoRef *GRef,char *Field,int32_t *Value)
libgeoref.GeoRef_GetField.argtypes = (ct.POINTER(TGeoRef),
                                     ct.c_char_p,
                                     ct.POINTER(ct.c_int32))
libgeoref.GeoRef_GetField.restype = ct.c_int32
GeoRef_GetField = libgeoref.GeoRef_GetField

# int32_t GeoRef_SetField(TGeoRef *GRef,char *Field,int32_t Value)
libgeoref.GeoRef_SetField.argtypes = (ct.POINTER(TGeoRef),
                                     ct.c_char_p,
                                     ct.c_int32)
libgeoref.GeoRef_SetField.restype = ct.c_int32
GeoRef_SetField = libgeoref.GeoRef_SetField

# int32_t GeoRef_WriteDesc(TGeoRef *GRef,char *Name,fst_file *File)
libgeoref.GeoRef_WriteDesc.argtypes = (ct.POINTER(TGeoRef),
                                      ct.c_char_p,
                                      ct.POINTER(fst_file))
libgeoref.GeoRef_WriteDesc.restype = ct.c_int32
GeoRef_WriteDesc = libgeoref.GeoRef_WriteDesc

# int32_t GeoRef_Write(TGeoRef *GRef,char *Name,fst_file *File)
libgeoref.GeoRef_Write.argtypes = (ct.POINTER(TGeoRef),
                                  ct.c_char_p,
                                  ct.POINTER(fst_file))
libgeoref.GeoRef_Write.restype = ct.c_int32
GeoRef_Write = libgeoref.GeoRef_Write

# int32_t GeoRef_CopyDesc(fst_file *FileTo,fst_record* Rec)
libgeoref.GeoRef_CopyDesc.argtypes = (ct.POINTER(fst_file),
                                     ct.POINTER(fst_record))
libgeoref.GeoRef_CopyDesc.restype = ct.c_int32
GeoRef_CopyDesc = libgeoref.GeoRef_CopyDesc

# int32_t GeoRef_Read(TGeoRef *GRef,char *Name,fst_file *File)
libgeoref.GeoRef_Read.argtypes = (ct.POINTER(TGeoRef),
                                 ct.c_char_p,
                                 ct.POINTER(fst_file))
libgeoref.GeoRef_Read.restype = ct.c_int32
GeoRef_Read = libgeoref.GeoRef_Read

# lib/OGR.c

# [ifdef: HAVE_GDAL] OGR_Layer* OGR_LayerNew(const char* Name)
libgeoref.OGR_LayerNew.argtypes = (ct.c_char_p,)
libgeoref.OGR_LayerNew.restype = ct.POINTER(OGR_Layer)
OGR_LayerNew = libgeoref.OGR_LayerNew

# [ifdef: HAVE_GDAL] void OGR_LayerFree(OGR_Layer *Layer)
libgeoref.OGR_LayerFree.argtypes = (ct.POINTER(OGR_Layer),)
libgeoref.OGR_LayerFree.restype = None
OGR_LayerFree = libgeoref.OGR_LayerFree

# [ifdef: HAVE_GDAL] OGR_Layer* OGR_LayerRead(const char *Path)
libgeoref.OGR_LayerRead.argtypes = (ct.c_char_p,)
libgeoref.OGR_LayerRead.restype = ct.POINTER(OGR_Layer)
OGR_LayerRead = libgeoref.OGR_LayerRead

# [ifdef: HAVE_GDAL] int32_t OGR_LayerWrite(OGR_Layer *Layer,const char *Path,const char *Format)
libgeoref.OGR_LayerWrite.argtypes = (ct.POINTER(OGR_Layer), ct.c_char_p, ct.c_char_p)
libgeoref.OGR_LayerWrite.restype = ct.c_int32
OGR_LayerWrite = libgeoref.OGR_LayerWrite

# [ifdef: HAVE_GDAL] int32_t OGR_LayerImport(OGR_Layer *Layer,const char *Path)
libgeoref.OGR_LayerImport.argtypes = (ct.POINTER(OGR_Layer), ct.c_char_p)
libgeoref.OGR_LayerImport.restype = ct.c_int32
OGR_LayerImport = libgeoref.OGR_LayerImport

# [ifdef: HAVE_GDAL] int32_t OGR_LayerExport(OGR_Layer *Layer,const char *Path)
libgeoref.OGR_LayerExport.argtypes = (ct.POINTER(OGR_Layer), ct.c_char_p)
libgeoref.OGR_LayerExport.restype = ct.c_int32
OGR_LayerExport = libgeoref.OGR_LayerExport

# [ifdef: HAVE_GDAL] int32_t OGR_LayerCopy(OGR_Layer *To,OGR_Layer *From)
libgeoref.OGR_LayerCopy.argtypes = (ct.POINTER(OGR_Layer), ct.POINTER(OGR_Layer))
libgeoref.OGR_LayerCopy.restype = ct.c_int32
OGR_LayerCopy = libgeoref.OGR_LayerCopy

# [ifdef: HAVE_GDAL] int32_t OGR_LayerClear(OGR_Layer *Layer)
libgeoref.OGR_LayerClear.argtypes = (ct.POINTER(OGR_Layer),)
libgeoref.OGR_LayerClear.restype = ct.c_int32
OGR_LayerClear = libgeoref.OGR_LayerClear

# [ifdef: HAVE_GDAL] int32_t OGR_LayerParse(OGR_Layer *Layer,const char *String)
libgeoref.OGR_LayerParse.argtypes = (ct.POINTER(OGR_Layer), ct.c_char_p)
libgeoref.OGR_LayerParse.restype = ct.c_int32
OGR_LayerParse = libgeoref.OGR_LayerParse

# [ifdef: HAVE_GDAL] int32_t OGR_LayerParseFile(OGR_Layer *Layer,const char *File)
libgeoref.OGR_LayerParseFile.argtypes = (ct.POINTER(OGR_Layer), ct.c_char_p)
libgeoref.OGR_LayerParseFile.restype = ct.c_int32
OGR_LayerParseFile = libgeoref.OGR_LayerParseFile

# [ifdef: HAVE_GDAL] int32_t OGR_LayerPrint(OGR_Layer *Layer)
libgeoref.OGR_LayerPrint.argtypes = (ct.POINTER(OGR_Layer),)
libgeoref.OGR_LayerPrint.restype = ct.c_int32
OGR_LayerPrint = libgeoref.OGR_LayerPrint

# [ifdef: HAVE_GDAL] int32_t OGR_LayerClean(OGR_Layer *Layer)
libgeoref.OGR_LayerClean.argtypes = (ct.POINTER(OGR_Layer),)
libgeoref.OGR_LayerClean.restype = ct.c_int32
OGR_LayerClean = libgeoref.OGR_LayerClean

# [ifdef: HAVE_GDAL] int32_t OGR_LayerIntersect(OGR_Layer *Layer0,OGR_Layer *Layer1)
libgeoref.OGR_LayerIntersect.argtypes = (ct.POINTER(OGR_Layer), ct.POINTER(OGR_Layer))
libgeoref.OGR_LayerIntersect.restype = ct.c_int32
OGR_LayerIntersect = libgeoref.OGR_LayerIntersect

# [ifdef: HAVE_GDAL] int32_t OGR_LayerUnion(OGR_Layer *Layer0,OGR_Layer *Layer1)
libgeoref.OGR_LayerUnion.argtypes = (ct.POINTER(OGR_Layer), ct.POINTER(OGR_Layer))
libgeoref.OGR_LayerUnion.restype = ct.c_int32
OGR_LayerUnion = libgeoref.OGR_LayerUnion

# [ifdef: HAVE_GDAL] int32_t OGR_LayerDifference(OGR_Layer *Layer0,OGR_Layer *Layer1)
libgeoref.OGR_LayerDifference.argtypes = (ct.POINTER(OGR_Layer), ct.POINTER(OGR_Layer))
libgeoref.OGR_LayerDifference.restype = ct.c_int32
OGR_LayerDifference = libgeoref.OGR_LayerDifference

# [ifdef: HAVE_GDAL] int32_t OGR_LayerSymDifference(OGR_Layer *Layer0,OGR_Layer *Layer1)
libgeoref.OGR_LayerSymDifference.argtypes = (ct.POINTER(OGR_Layer), ct.POINTER(OGR_Layer))
libgeoref.OGR_LayerSymDifference.restype = ct.c_int32
OGR_LayerSymDifference = libgeoref.OGR_LayerSymDifference

# [ifdef: HAVE_GDAL] int32_t OGR_LayerBuffer(OGR_Layer *Layer,double Width,int32_t Segments)
libgeoref.OGR_LayerBuffer.argtypes = (ct.POINTER(OGR_Layer), ct.c_double, ct.c_int32)
libgeoref.OGR_LayerBuffer.restype = ct.c_int32
OGR_LayerBuffer = libgeoref.OGR_LayerBuffer

# [ifdef: HAVE_GDAL] int32_t OGR_LayerCentroid(OGR_Layer *Layer)
libgeoref.OGR_LayerCentroid.argtypes = (ct.POINTER(OGR_Layer),)
libgeoref.OGR_LayerCentroid.restype = ct.c_int32
OGR_LayerCentroid = libgeoref.OGR_LayerCentroid

# [ifdef: HAVE_GDAL] int32_t OGR_LayerConvexHull(OGR_Layer *Layer)
libgeoref.OGR_LayerConvexHull.argtypes = (ct.POINTER(OGR_Layer),)
libgeoref.OGR_LayerConvexHull.restype = ct.c_int32
OGR_LayerConvexHull = libgeoref.OGR_LayerConvexHull

# [ifdef: HAVE_GDAL] int32_t OGR_LayerBoundary(OGR_Layer *Layer)
libgeoref.OGR_LayerBoundary.argtypes = (ct.POINTER(OGR_Layer),)
libgeoref.OGR_LayerBoundary.restype = ct.c_int32
OGR_LayerBoundary = libgeoref.OGR_LayerBoundary

# [ifdef: HAVE_GDAL] double OGR_LayerArea(OGR_Layer *Layer)
libgeoref.OGR_LayerArea.argtypes = (ct.POINTER(OGR_Layer),)
libgeoref.OGR_LayerArea.restype = ct.c_double
OGR_LayerArea = libgeoref.OGR_LayerArea

# [ifdef: HAVE_GDAL] double OGR_LayerLength(OGR_Layer *Layer)
libgeoref.OGR_LayerLength.argtypes = (ct.POINTER(OGR_Layer),)
libgeoref.OGR_LayerLength.restype = ct.c_double
OGR_LayerLength = libgeoref.OGR_LayerLength

# [ifdef: HAVE_GDAL] int32_t OGR_LayerSimplify(OGR_Layer *Layer,double Tolerance)
libgeoref.OGR_LayerSimplify.argtypes = (ct.POINTER(OGR_Layer), ct.c_double)
libgeoref.OGR_LayerSimplify.restype = ct.c_int32
OGR_LayerSimplify = libgeoref.OGR_LayerSimplify

# [ifdef: HAVE_GDAL] int32_t OGR_LayerSegmentize(OGR_Layer *Layer,double Length)
libgeoref.OGR_LayerSegmentize.argtypes = (ct.POINTER(OGR_Layer), ct.c_double)
libgeoref.OGR_LayerSegmentize.restype = ct.c_int32
OGR_LayerSegmentize = libgeoref.OGR_LayerSegmentize

# [ifdef: HAVE_GDAL] int32_t OGR_LayerProject(OGR_Layer *Layer,const char *Proj)
libgeoref.OGR_LayerProject.argtypes = (ct.POINTER(OGR_Layer), ct.c_char_p)
libgeoref.OGR_LayerProject.restype = ct.c_int32
OGR_LayerProject = libgeoref.OGR_LayerProject

# [ifdef: HAVE_GDAL] int32_t OGR_LayerSelectBox(OGR_Layer *Layer,double X0,double Y0,double X1,double Y1)
libgeoref.OGR_LayerSelectBox.argtypes = (ct.POINTER(OGR_Layer),
                                        ct.c_double, ct.c_double,
                                        ct.c_double, ct.c_double)
libgeoref.OGR_LayerSelectBox.restype = ct.c_int32
OGR_LayerSelectBox = libgeoref.OGR_LayerSelectBox

# [ifdef: HAVE_GDAL] int32_t OGR_LayerSelectPoint(OGR_Layer *Layer,double X,double Y,double Delta)
libgeoref.OGR_LayerSelectPoint.argtypes = (ct.POINTER(OGR_Layer),
                                          ct.c_double, ct.c_double, ct.c_double)
libgeoref.OGR_LayerSelectPoint.restype = ct.c_int32
OGR_LayerSelectPoint = libgeoref.OGR_LayerSelectPoint

# [ifdef: HAVE_GDAL] int32_t OGR_LayerSelectLine(OGR_Layer *Layer,double X0,double Y0,double X1,double Y1,double Delta)
libgeoref.OGR_LayerSelectLine.argtypes = (ct.POINTER(OGR_Layer),
                                         ct.c_double, ct.c_double,
                                         ct.c_double, ct.c_double, ct.c_double)
libgeoref.OGR_LayerSelectLine.restype = ct.c_int32
OGR_LayerSelectLine = libgeoref.OGR_LayerSelectLine

# [ifdef: HAVE_GDAL] int32_t OGR_LayerSelectPoly(OGR_Layer *Layer,OGR_Layer *Poly)
libgeoref.OGR_LayerSelectPoly.argtypes = (ct.POINTER(OGR_Layer), ct.POINTER(OGR_Layer))
libgeoref.OGR_LayerSelectPoly.restype = ct.c_int32
OGR_LayerSelectPoly = libgeoref.OGR_LayerSelectPoly

# [ifdef: HAVE_GDAL] int32_t OGR_LayerSelectSQL(OGR_Layer *Layer,const char *SQL)
libgeoref.OGR_LayerSelectSQL.argtypes = (ct.POINTER(OGR_Layer), ct.c_char_p)
libgeoref.OGR_LayerSelectSQL.restype = ct.c_int32
OGR_LayerSelectSQL = libgeoref.OGR_LayerSelectSQL

# [ifdef: HAVE_GDAL] int32_t OGR_LayerSelectId(OGR_Layer *Layer,int32_t Id)
libgeoref.OGR_LayerSelectId.argtypes = (ct.POINTER(OGR_Layer), ct.c_int32)
libgeoref.OGR_LayerSelectId.restype = ct.c_int32
OGR_LayerSelectId = libgeoref.OGR_LayerSelectId

# [ifdef: HAVE_GDAL] int32_t OGR_LayerSelectClear(OGR_Layer *Layer)
libgeoref.OGR_LayerSelectClear.argtypes = (ct.POINTER(OGR_Layer),)
libgeoref.OGR_LayerSelectClear.restype = ct.c_int32
OGR_LayerSelectClear = libgeoref.OGR_LayerSelectClear

# [ifdef: HAVE_GDAL] int32_t OGR_LayerSelectRender(OGR_Layer *Layer,TGeoRef *Ref,double *Data,int32_t Value)
libgeoref.OGR_LayerSelectRender.argtypes = (ct.POINTER(OGR_Layer),
                                           ct.POINTER(TGeoRef),
                                           npc.ndpointer(dtype=np.float64),
                                           ct.c_int32)
libgeoref.OGR_LayerSelectRender.restype = ct.c_int32
OGR_LayerSelectRender = libgeoref.OGR_LayerSelectRender

# [ifdef: HAVE_GDAL] int32_t OGR_LayerRender(OGR_Layer *Layer,TGeoRef *Ref,double *Data,int32_t Value)
libgeoref.OGR_LayerRender.argtypes = (ct.POINTER(OGR_Layer),
                                     ct.POINTER(TGeoRef),
                                     npc.ndpointer(dtype=np.float64),
                                     ct.c_int32)
libgeoref.OGR_LayerRender.restype = ct.c_int32
OGR_LayerRender = libgeoref.OGR_LayerRender

# lib/Vertex.c

# void Vertex_Map(Vect2d P[4],double *LX,double *LY,double WX,double WY)
libgeoref.Vertex_Map.argtypes = (ct.POINTER(Vect2d),
                                ct.POINTER(ct.c_double),
                                ct.POINTER(ct.c_double),
                                ct.c_double, ct.c_double)
libgeoref.Vertex_Map.restype = None
Vertex_Map = libgeoref.Vertex_Map

# void VertexGradient(TDef *Def,Vect3d Nr)
libgeoref.VertexGradient.argtypes = (ct.POINTER(TDef), ct.POINTER(Vect3d))
libgeoref.VertexGradient.restype = None
VertexGradient = libgeoref.VertexGradient

# float VertexVal(TDef *Def,int32_t Idx,double X,double Y,double Z)
libgeoref.VertexVal.argtypes = (ct.POINTER(TDef),
                               ct.c_int32,
                               ct.c_double, ct.c_double, ct.c_double)
libgeoref.VertexVal.restype = ct.c_float
VertexVal = libgeoref.VertexVal

# double VertexValV(TDef *Def,double X,double Y,double Z,Vect3d V)
libgeoref.VertexValV.argtypes = (ct.POINTER(TDef),
                                ct.c_double, ct.c_double, ct.c_double,
                                ct.POINTER(Vect3d))
libgeoref.VertexValV.restype = ct.c_double
VertexValV = libgeoref.VertexValV

# float VertexValS(double *Data,char *Mask,int32_t NI,int32_t NJ,double X,double Y,char Geo)
libgeoref.VertexValS.argtypes = (npc.ndpointer(dtype=np.float64),
                                ct.c_char_p,
                                ct.c_int32, ct.c_int32,
                                ct.c_double, ct.c_double,
                                ct.c_char)
libgeoref.VertexValS.restype = ct.c_float
VertexValS = libgeoref.VertexValS

# int32_t VertexLoc(Vect3d **Pos,TDef *Def,Vect3d Vr,double X,double Y,double Z,int32_t Wrap)
libgeoref.VertexLoc.argtypes = (ct.POINTER(ct.POINTER(Vect3d)),
                               ct.POINTER(TDef),
                               ct.POINTER(Vect3d),
                               ct.c_double, ct.c_double, ct.c_double,
                               ct.c_int32)
libgeoref.VertexLoc.restype = ct.c_int32
VertexLoc = libgeoref.VertexLoc

# void VertexInterp(Vect3d Pi,Vect3d P0,Vect3d P1,double V0,double V1,double Level)
libgeoref.VertexInterp.argtypes = (ct.POINTER(Vect3d),
                                  ct.POINTER(Vect3d),
                                  ct.POINTER(Vect3d),
                                  ct.c_double, ct.c_double, ct.c_double)
libgeoref.VertexInterp.restype = None
VertexInterp = libgeoref.VertexInterp

# lib/ZRef.c

# TZRef* ZRef_New()
libgeoref.ZRef_New.argtypes = None
libgeoref.ZRef_New.restype = ct.POINTER(TZRef)
ZRef_New = libgeoref.ZRef_New

# void ZRef_Free(TZRef *ZRef)
libgeoref.ZRef_Free.argtypes = (ct.POINTER(TZRef),)
libgeoref.ZRef_Free.restype = None
ZRef_Free = libgeoref.ZRef_Free

# void ZRef_Clear(TZRef *ZRef)
libgeoref.ZRef_Clear.argtypes = (ct.POINTER(TZRef),)
libgeoref.ZRef_Clear.restype = None
ZRef_Clear = libgeoref.ZRef_Clear

# TZRef* ZRef_Copy(TZRef *ZRef)
libgeoref.ZRef_Copy.argtypes = (ct.POINTER(TZRef),)
libgeoref.ZRef_Copy.restype = ct.POINTER(TZRef)
ZRef_Copy = libgeoref.ZRef_Copy

# int32_t ZRef_Size(TZRef *ZRef)
libgeoref.ZRef_Size.argtypes = (ct.POINTER(TZRef),)
libgeoref.ZRef_Size.restype = ct.c_int32
ZRef_Size = libgeoref.ZRef_Size

# int32_t ZRef_Resize(TZRef *ZRef,int32_t NK)
libgeoref.ZRef_Resize.argtypes = (ct.POINTER(TZRef), ct.c_int32)
libgeoref.ZRef_Resize.restype = ct.c_int32
ZRef_Resize = libgeoref.ZRef_Resize

# int32_t ZRef_Define(TZRef *ZRef,int32_t NK,int32_t Type,float *Levels)
libgeoref.ZRef_Define.argtypes = (ct.POINTER(TZRef),
                                 ct.c_int32, ct.c_int32,
                                 npc.ndpointer(dtype=np.float32))
libgeoref.ZRef_Define.restype = ct.c_int32
ZRef_Define = libgeoref.ZRef_Define

# int32_t ZRef_Equal(TZRef *ZRefA,TZRef *ZRefB)
libgeoref.ZRef_Equal.argtypes = (ct.POINTER(TZRef), ct.POINTER(TZRef))
libgeoref.ZRef_Equal.restype = ct.c_int32
ZRef_Equal = libgeoref.ZRef_Equal

# int32_t ZRef_DecodeRPN(TZRef *ZRef,fst_file* File)
libgeoref.ZRef_DecodeRPN.argtypes = (ct.POINTER(TZRef), ct.POINTER(fst_file))
libgeoref.ZRef_DecodeRPN.restype = ct.c_int32
ZRef_DecodeRPN = libgeoref.ZRef_DecodeRPN

# int32_t ZRef_DecodeHY(TZRef *ZRef,const fst_record* restrict const H)
libgeoref.ZRef_DecodeHY.argtypes = (ct.POINTER(TZRef), ct.POINTER(fst_record))
libgeoref.ZRef_DecodeHY.restype = ct.c_int32
ZRef_DecodeHY = libgeoref.ZRef_DecodeHY

# int32_t ZRef_GetLevels(TZRef *ZRef,const fst_record* restrict const H,int32_t Order)
libgeoref.ZRef_GetLevels.argtypes = (ct.POINTER(TZRef),
                                    ct.POINTER(fst_record),
                                    ct.c_int32)
libgeoref.ZRef_GetLevels.restype = ct.c_int32
ZRef_GetLevels = libgeoref.ZRef_GetLevels

# int32_t ZRef_Level2IP(double Level,int32_t Type)
libgeoref.ZRef_Level2IP.argtypes = (ct.c_double, ct.c_int32)
libgeoref.ZRef_Level2IP.restype = ct.c_int32
ZRef_Level2IP = libgeoref.ZRef_Level2IP

# double ZRef_IP2Level(int32_t IP,int32_t Type)
libgeoref.ZRef_IP2Level.argtypes = (ct.c_int32, ct.c_int32)
libgeoref.ZRef_IP2Level.restype = ct.c_double
ZRef_IP2Level = libgeoref.ZRef_IP2Level

# int32_t ZRef_IP2Level_Parse(char *String,int32_t *IP0,int32_t *IP1)
libgeoref.ZRef_IP2Level_Parse.argtypes = (ct.c_char_p,
                                         ct.POINTER(ct.c_int32),
                                         ct.POINTER(ct.c_int32))
libgeoref.ZRef_IP2Level_Parse.restype = ct.c_int32
ZRef_IP2Level_Parse = libgeoref.ZRef_IP2Level_Parse

# int32_t ZRef_Level2IP_Parse(char *String,double *Level0,double *Level1)
libgeoref.ZRef_Level2IP_Parse.argtypes = (ct.c_char_p,
                                         ct.POINTER(ct.c_double),
                                         ct.POINTER(ct.c_double))
libgeoref.ZRef_Level2IP_Parse.restype = ct.c_int32
ZRef_Level2IP_Parse = libgeoref.ZRef_Level2IP_Parse

# double ZRef_IP2Meter(int32_t IP)
libgeoref.ZRef_IP2Meter.argtypes = (ct.c_int32,)
libgeoref.ZRef_IP2Meter.restype = ct.c_double
ZRef_IP2Meter = libgeoref.ZRef_IP2Meter

# int32_t ZRef_Meter2IP(double Level)
libgeoref.ZRef_Meter2IP.argtypes = (ct.c_double,)
libgeoref.ZRef_Meter2IP.restype = ct.c_int32
ZRef_Meter2IP = libgeoref.ZRef_Meter2IP

# int32_t ZRef_KCube2Meter(TZRef* restrict const ZRef,float *GZ,const int32_t NIJ,float *Height)
libgeoref.ZRef_KCube2Meter.argtypes = (ct.POINTER(TZRef),
                                      npc.ndpointer(dtype=np.float32),
                                      ct.c_int32,
                                      npc.ndpointer(dtype=np.float32))
libgeoref.ZRef_KCube2Meter.restype = ct.c_int32
ZRef_KCube2Meter = libgeoref.ZRef_KCube2Meter

# int32_t ZRef_BuildSigma(TZRef *ZRef,float *Levels,int32_t NK)
libgeoref.ZRef_BuildSigma.argtypes = (ct.POINTER(TZRef),
                                     npc.ndpointer(dtype=np.float32),
                                     ct.c_int32)
libgeoref.ZRef_BuildSigma.restype = ct.c_int32
ZRef_BuildSigma = libgeoref.ZRef_BuildSigma

# int32_t ZRef_BuildEta(TZRef *ZRef,float *Levels,int32_t NK)
libgeoref.ZRef_BuildEta.argtypes = (ct.POINTER(TZRef),
                                   npc.ndpointer(dtype=np.float32),
                                   ct.c_int32)
libgeoref.ZRef_BuildEta.restype = ct.c_int32
ZRef_BuildEta = libgeoref.ZRef_BuildEta

# int32_t ZRef_BuildPres(TZRef *ZRef,float *Levels,int32_t NK)
libgeoref.ZRef_BuildPres.argtypes = (ct.POINTER(TZRef),
                                    npc.ndpointer(dtype=np.float32),
                                    ct.c_int32)
libgeoref.ZRef_BuildPres.restype = ct.c_int32
ZRef_BuildPres = libgeoref.ZRef_BuildPres

# int32_t ZRef_BuildHeight(TZRef *ZRef,float *Levels,int32_t NK)
libgeoref.ZRef_BuildHeight.argtypes = (ct.POINTER(TZRef),
                                      npc.ndpointer(dtype=np.float32),
                                      ct.c_int32)
libgeoref.ZRef_BuildHeight.restype = ct.c_int32
ZRef_BuildHeight = libgeoref.ZRef_BuildHeight

# int32_t ZRef_BuildHybrid(TZRef *ZRef,float *A,float *B,int32_t NK,float P0)
libgeoref.ZRef_BuildHybrid.argtypes = (ct.POINTER(TZRef),
                                      npc.ndpointer(dtype=np.float32),
                                      npc.ndpointer(dtype=np.float32),
                                      ct.c_int32,
                                      ct.c_float)
libgeoref.ZRef_BuildHybrid.restype = ct.c_int32
ZRef_BuildHybrid = libgeoref.ZRef_BuildHybrid

# int32_t ZRef_BuildThermo(TZRef *ZRef,float *Levels,int32_t NK)
libgeoref.ZRef_BuildThermo.argtypes = (ct.POINTER(TZRef),
                                      npc.ndpointer(dtype=np.float32),
                                      ct.c_int32)
libgeoref.ZRef_BuildThermo.restype = ct.c_int32
ZRef_BuildThermo = libgeoref.ZRef_BuildThermo

# int32_t ZRef_Level2Index(TZRef *ZRef,float Level,float *W0,float *W1,int32_t *K0,int32_t *K1)
libgeoref.ZRef_Level2Index.argtypes = (ct.POINTER(TZRef),
                                      ct.c_float,
                                      ct.POINTER(ct.c_float),
                                      ct.POINTER(ct.c_float),
                                      ct.POINTER(ct.c_int32),
                                      ct.POINTER(ct.c_int32))
libgeoref.ZRef_Level2Index.restype = ct.c_int32
ZRef_Level2Index = libgeoref.ZRef_Level2Index

# int32_t ZRef_Pressure(TZRef *ZRef,float *Pres,float PS,int32_t NK)
libgeoref.ZRef_Pressure.argtypes = (ct.POINTER(TZRef),
                                   npc.ndpointer(dtype=np.float32),
                                   ct.c_float,
                                   ct.c_int32)
libgeoref.ZRef_Pressure.restype = ct.c_int32
ZRef_Pressure = libgeoref.ZRef_Pressure

# int32_t ZRef_Height(TZRef *ZRef,float *Height,float *GZ,int32_t NK)
libgeoref.ZRef_Height.argtypes = (ct.POINTER(TZRef),
                                 npc.ndpointer(dtype=np.float32),
                                 npc.ndpointer(dtype=np.float32),
                                 ct.c_int32)
libgeoref.ZRef_Height.restype = ct.c_int32
ZRef_Height = libgeoref.ZRef_Height

# int32_t ZRef_Intervals(TZRef *ZRef,float *Intervals,int32_t NK)
libgeoref.ZRef_Intervals.argtypes = (ct.POINTER(TZRef),
                                    npc.ndpointer(dtype=np.float32),
                                    ct.c_int32)
libgeoref.ZRef_Intervals.restype = ct.c_int32
ZRef_Intervals = libgeoref.ZRef_Intervals

# int32_t ZRef_Closest(TZRef *ZRef,float Level)
libgeoref.ZRef_Closest.argtypes = (ct.POINTER(TZRef), ct.c_float)
libgeoref.ZRef_Closest.restype = ct.c_int32
ZRef_Closest = libgeoref.ZRef_Closest

# int32_t ZRef_WithinLayer(TZRef *ZRef,float Level,int32_t Layer)
libgeoref.ZRef_WithinLayer.argtypes = (ct.POINTER(TZRef),
                                      ct.c_float,
                                      ct.c_int32)
libgeoref.ZRef_WithinLayer.restype = ct.c_int32
ZRef_WithinLayer = libgeoref.ZRef_WithinLayer

# int32_t ZRef_WithinRange(TZRef *ZRef,float Level)
libgeoref.ZRef_WithinRange.argtypes = (ct.POINTER(TZRef), ct.c_float)
libgeoref.ZRef_WithinRange.restype = ct.c_int32
ZRef_WithinRange = libgeoref.ZRef_WithinRange

# int32_t ZRef_Between(TZRef *ZRef,float Level0,float Level1)
libgeoref.ZRef_Between.argtypes = (ct.POINTER(TZRef),
                                  ct.c_float,
                                  ct.c_float)
libgeoref.ZRef_Between.restype = ct.c_int32
ZRef_Between = libgeoref.ZRef_Between

# int32_t ZRef_IPFormat(char *Buf,int32_t IP,int32_t Interval)
libgeoref.ZRef_IPFormat.argtypes = (ct.c_char_p,
                                   ct.c_int32,
                                   ct.c_int32)
libgeoref.ZRef_IPFormat.restype = ct.c_int32
ZRef_IPFormat = libgeoref.ZRef_IPFormat

# const char** ZRef_LevelNames()
libgeoref.ZRef_LevelNames.argtypes = None
libgeoref.ZRef_LevelNames.restype = ct.POINTER(ct.c_char_p)
ZRef_LevelNames = libgeoref.ZRef_LevelNames

# const char* ZRef_LevelName(int32_t Type)
libgeoref.ZRef_LevelName.argtypes = (ct.c_int32,)
libgeoref.ZRef_LevelName.restype = ct.c_char_p
ZRef_LevelName = libgeoref.ZRef_LevelName

# const char** ZRef_LevelUnits()
libgeoref.ZRef_LevelUnits.argtypes = None
libgeoref.ZRef_LevelUnits.restype = ct.POINTER(ct.c_char_p)
ZRef_LevelUnits = libgeoref.ZRef_LevelUnits

# const char* ZRef_LevelUnit(int32_t Type)
libgeoref.ZRef_LevelUnit.argtypes = (ct.c_int32,)
libgeoref.ZRef_LevelUnit.restype = ct.c_char_p
ZRef_LevelUnit = libgeoref.ZRef_LevelUnit

# lib/ZRefInterp.c

# TZRefInterp* ZRefInterp_New()
libgeoref.ZRefInterp_New.argtypes = None
libgeoref.ZRefInterp_New.restype = ct.POINTER(TZRefInterp)
ZRefInterp_New = libgeoref.ZRefInterp_New

# int32_t ZRefInterp_Free(TZRefInterp *Interp)
libgeoref.ZRefInterp_Free.argtypes = (ct.POINTER(TZRefInterp),)
libgeoref.ZRefInterp_Free.restype = ct.c_int32
ZRefInterp_Free = libgeoref.ZRefInterp_Free

# void ZRefInterp_Clear(TZRefInterp *Interp)
libgeoref.ZRefInterp_Clear.argtypes = (ct.POINTER(TZRefInterp),)
libgeoref.ZRefInterp_Clear.restype = None
ZRefInterp_Clear = libgeoref.ZRefInterp_Clear

# TZRefInterp* ZRefInterp_Define(TZRef *ZRefDest,TZRef *ZRefSrc,const int32_t NI,const int32_t NJ)
libgeoref.ZRefInterp_Define.argtypes = (ct.POINTER(TZRef),
                                       ct.POINTER(TZRef),
                                       ct.c_int32,
                                       ct.c_int32)
libgeoref.ZRefInterp_Define.restype = ct.POINTER(TZRefInterp)
ZRefInterp_Define = libgeoref.ZRefInterp_Define

# int32_t ZRefInterp(TZRefInterp *Interp,float *stateOut,float *stateIn,float *derivOut,float *derivIn,float extrapGuideDown,float extrapGuideUp)
libgeoref.ZRefInterp.argtypes = (ct.POINTER(TZRefInterp),
                                npc.ndpointer(dtype=np.float32),
                                npc.ndpointer(dtype=np.float32),
                                npc.ndpointer(dtype=np.float32),
                                npc.ndpointer(dtype=np.float32),
                                ct.c_float,
                                ct.c_float)
libgeoref.ZRefInterp.restype = ct.c_int32
ZRefInterp = libgeoref.ZRefInterp

# int32_t ZRefInterp_SetOption(const char *Option,const char *Value)
libgeoref.ZRefInterp_SetOption.argtypes = (ct.c_char_p, ct.c_char_p)
libgeoref.ZRefInterp_SetOption.restype = ct.c_int32
ZRefInterp_SetOption = libgeoref.ZRefInterp_SetOption

# int32_t ZRefInterp_SetOptioni(const unsigned char option)
libgeoref.ZRefInterp_SetOptioni.argtypes = (ct.c_ubyte,)
libgeoref.ZRefInterp_SetOptioni.restype = ct.c_int32
ZRefInterp_SetOptioni = libgeoref.ZRefInterp_SetOptioni




