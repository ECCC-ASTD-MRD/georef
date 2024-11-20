from ctypes import *

class TRPNFile(Structure):
    _fields_ = [
        ("CId", c_char_p),
        ("Name", c_char_p),
        ("Open", c_int),
        ("Id", c_uint),
        ("NRef", c_int),
        ("Mode", c_char)
    ]

class TRPNHeader(Structure):
    _fields_ = [
        ("File", POINTER(TRPNFile)),
        ("FID", c_int),
        ("KEY", c_int),
        ("DATEO", c_int),
        ("DATEV", c_int),
        ("DEET", c_int),
        ("NPAS", c_int),
        ("NBITS", c_int),
        ("DATYP", c_int),
        ("IP1", c_int),
        ("IP2", c_int),
        ("IP3", c_int),
        ("NI", c_int),
        ("NJ", c_int),
        ("NK", c_int),
        ("NIJ", c_int),
        ("IG", c_int * 4),
        ("IGREF", c_int * 4),
        ("XG", c_float * 4),
        ("XGREF", c_float * 4),
        ("SWA", c_int),
        ("LNG", c_int),
        ("DLTF", c_int),
        ("UBC", c_int),
        ("EX1", c_int),
        ("EX2", c_int),
        ("EX3", c_int),
        ("TYPVAR", c_char * 3),
        ("NOMVAR", c_char * 5),
        ("ETIKET", c_char * 13),
        ("GRTYP", c_char * 2),
        ("GRREF", c_char * 2)
    ]

class TGeoRef(Structure):
    pass

class TGeoOptions(Structure):
    _fields_ = [
        ("InterpDegree", c_int),
        ("ExtrapDegree", c_int),
        ("ExtrapValue", c_double),
        ("SubGrid", c_int),
        ("Transform", c_int),
        ("CIndex", c_int),
        ("Symmetric", c_int),
        ("WeightNum", c_int),
        ("PolarCorrect", c_char), #TODO: correct type?
        ("VectorMode", c_char),
        ("DistTreshold", c_float),
        ("LonRef", c_float)
    ]

class TRotationTransform(Structure):
    _fields_ = [
        ("Lat", c_double),
        ("Lon", c_double),
        ("Angle", c_double),
        ("SinTheta", c_double),
        ("CosTheta", c_double),
        ("SinPhi", c_double),
        ("CosPhi", c_double)
    ]

class TGeoZone(Structure):
    _fields_ = [
        ("npts", c_int),
        ("x", POINTER(c_double)),
        ("y", POINTER(c_double)),
        (("idx", POINTER(c_int)))
    ]

class TGridSet(Structure):
    _fields_ = [
        ("RefFrom", POINTER(TGeoRef)),
        ("zones", TGeoZone * 5), #TODO: define SET_NZONES 5
        ("flags", c_int),
        ("x", POINTER(c_double)), ("y", POINTER(c_double)),
        (("mask_in", POINTER(c_int))), (("mask_out", POINTER(c_int))),
        ("yin_maskout", POINTER(c_float)), ("yan_maskout", POINTER(c_float)),
        (("yinlat", POINTER(c_int))), (("yinlon", POINTER(c_int))), (("yanlat", POINTER(c_int))), (("yanlon", POINTER(c_int))),
        (("yin2yin_lat", POINTER(c_int))), (("yin2yin_lon", POINTER(c_int))), (("yan2yin_lat", POINTER(c_int))), (("yan2yin_lon", POINTER(c_int))),
        (("yin2yan_lat", POINTER(c_int))), (("yin2yan_lon", POINTER(c_int))), (("yan2yan_lat", POINTER(c_int))), (("yan2yan_lon", POINTER(c_int))),
        (("yin2yin_x", POINTER(c_int))), (("yin2yin_y", POINTER(c_int))), (("yan2yin_x", POINTER(c_int))), (("yan2yin_y", POINTER(c_int))),
        (("yin2yan_x", POINTER(c_int))), (("yin2yan_y", POINTER(c_int))), (("yan2yan_x", POINTER(c_int))), (("yan2yan_y", POINTER(c_int))),
        ("yincount_yin", c_int), ("yancount_yin", c_int), ("yincount_yan", c_int), ("yancount_yan", c_int),
        ("wts", c_int),
        ("y", POINTER(c_double)),
        (("mask", POINTER(c_int))), (("idx", POINTER(c_int)))
    ]

TGeoRef._fields_ = [
    ("Name", c_char_p),
    ("NRef", c_int),
    ("RefFrom", POINTER(TGeoRef)),
    ("Subs", POINTER(POINTER(TGeoRef))), #TODO: pointer of pointer
    ("NbSub", c_int),
    ("Type", c_int),
    ("BD", c_int),
    ("NX", c_int),("NY", c_int),("X0", c_int),("Y0", c_int),("X1", c_int),("Y1", c_int),
    ("Extension", c_int),
    ("GRTYP", c_char * 3),
    ("Hemi", c_int),
    ("NbSet", c_int),
    ("Sets", POINTER(TGridSet)),("LastSet", POINTER(TGridSet)),
    ("NIdx", c_uint),("Idx", POINTER(c_uint)),
    ("Lat", POINTER(c_double)),("Lon", POINTER(c_double)),
    ("AX", POINTER(c_float)),("AY", POINTER(c_float)),("NCX", POINTER(c_float)),("NCY", POINTER(c_float)),("Hgt", POINTER(c_float)),
    ("Wght", POINTER(c_double)),
    ("String", c_char_p),
    ("Transform", POINTER(c_double)),("InvTransform", POINTER(c_double)),
    ("RotTransform", POINTER(TRotationTransform)),


#    void                         *GCPTransform;            ///< GPC derivative transform (1,2,3 order)
#    void                         *TPSTransform;            ///< GPC Thin Spline transform
#    void                         *RPCTransform;            ///< GPC Rigorous Projection Model transform
#    TQTree                       *QTree;                   ///< Quadtree index
#    OGREnvelope                   LLExtent;                ///< LatLon extent
#    OGRCoordinateTransformationH  Function,InvFunction;    ///< Projection functions
#    OGRSpatialReferenceH          Spatial;                 ///< Spatial reference

#    TCoord  Loc;                                           ///< (Radar) Localisation du centre de reference
#    double CTH,STH;                                        ///< (Radar) sin and cos of sweep angle
#    double ResR,ResA;                                      ///< (Radar) Resolutions en distance et azimuth
#    int    R;                                              ///< (Radar) Rayon autour du centre de reference en bin

    ("Options", TGeoOptions),

    # TODO: Is it the correct way ?
    # typedef int    (TGeoRef_Project)   (struct TGeoRef *Ref,double X,double Y,double *Lat,double *Lon,int Extrap,int Transform);
    ("Project", CFUNCTYPE(c_int, POINTER(TGeoRef), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int)),

    # typedef int    (TGeoRef_UnProject) (struct TGeoRef *Ref,double *X,double *Y,double Lat,double Lon,int Extrap,int Transform);
    ("UnProject", CFUNCTYPE(c_int, POINTER(TGeoRef), POINTER(c_double), POINTER(c_double), c_double, c_double, c_int, c_int)),

    # typedef int    (TGeoRef_LL2XY)     (struct TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb);
    ("XY2LL", CFUNCTYPE(c_int, POINTER(TGeoRef), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int)),

    # typedef int    (TGeoRef_XY2LL)     (struct TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
    ("LL2XY", CFUNCTYPE(c_int, POINTER(TGeoRef), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int)),
                         
#    TGeoRef_Value     *Value;                              ///< Valeur a un point xy
#    TGeoRef_Height    *Height;

    ("RPNHead", TRPNHeader),
    ("i1", c_int),("i2", c_int),("j1", c_int),("j2", c_int),
    ("mask", POINTER(c_int)),
    ("mymaskgrid", POINTER(TGeoRef)),
    ("mymaskgridi0", c_int),("mymaskgridi1", c_int),
    ("mymaskgridj0", c_int),("mymaskgridj1", c_int)

# #ifdef HAVE_RPNC   
#    int NC_Id,NC_NXDimId,NC_NYDimId;                       ///< netCDF identifiers
# #endif
    ]
