# typedef struct {
#    int npts;                 ///< Nombre de points
#    double *x,*y;             ///< Vecteur de coordonnees 
#    int *idx;                 ///< Indice du point dans le champ de destination
# } TGeoZone;

# typedef struct {
#    struct TGeoRef* RefFrom;
#    TGeoZone zones[SET_NZONES];
#    int flags;
#    double *x, *y;
#    int *mask_in, *mask_out;
#    float *yin_maskout,*yan_maskout;
#    double *yinlat,*yinlon,*yanlat,*yanlon;
#    double *yin2yin_lat,*yin2yin_lon,*yan2yin_lat,*yan2yin_lon;
#    double *yin2yan_lat,*yin2yan_lon,*yan2yan_lat,*yan2yan_lon;
#    double *yin2yin_x,*yin2yin_y,*yan2yin_x,*yan2yin_y;
#    double *yin2yan_x,*yin2yan_y,*yan2yan_x,*yan2yan_y;
#    int yincount_yin,yancount_yin,yincount_yan,yancount_yan;

#    int n_wts;                ///< Nombre de poids
#    double *wts;              ///< Tableau de poids
#    int *mask, *idx;          ///< Indice du point dans le champ de destination
# } TGridSet;

# typedef struct TGeoRef {
#    char*   Name;                                          ///< Reference name
#    int     NRef;                                          ///< Nombre de reference a la georeference
#    struct  TGeoRef  *RefFrom;                             ///< Georeference de reference (coupe verticale,...)
#    struct  TGeoRef **Subs;                                ///< Liste des sous grilles (GRTYP=U)
#    int     NbSub;                                         ///< Nombre de sous-grilles
#    int     Type;                                          ///< Parametre/Flags de grille
#    int     BD;                                            ///< Bordure
#    int     NX,NY,X0,Y0,X1,Y1;                             ///< Grid size and limits
#    int     Extension;                                     ///< related to the newtonian coefficient
#    char    GRTYP[3];                                      ///< Type de grille
#    int     Hemi;                                          ///< Hemisphere side (0=GLOBAL,1=NORTH,2=SOUTH)
#    int     NbSet;                                         ///< Nombre de set d'interpolation
#    TGridSet *Sets,*LastSet;                               ///< Tableau de set d'interpolation et du dernier utilise

#    unsigned int NIdx,*Idx;                                ///< Index dans les positions
#    double       *Lat,*Lon;                                ///< Coordonnees des points de grilles (Spherical)
#    float        *AX,*AY,*NCX,*NCY,*Hgt;                   ///< Axes de positionnement / deformation
#    double       *Wght;                                    ///< Barycentric weight array for TIN  (M grids)

#    char                         *String;                  ///< OpenGIS WKT String description
#    double                       *Transform,*InvTransform; ///< Transformation functions
#    TRotationTransform           *RotTransform;            ///< Rotation transform
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

#    TGeoOptions Options;                                   ///< Options for manipulations
#    TGeoRef_Project   *Project;                            ///< Transformation xy a latlon
#    TGeoRef_UnProject *UnProject;                          ///< Transformation latlon a xy
#    TGeoRef_Value     *Value;                              ///< Valeur a un point xy
#    TGeoRef_Distance  *Distance;                           
#    TGeoRef_Height    *Height;

#    _Grille struct from ezscint.h
#    TRPNHeader RPNHead;
#    int i1, i2, j1, j2;
#    int *mask;
#    struct TGeoRef *mymaskgrid;
#    int mymaskgridi0,mymaskgridi1;
#    int mymaskgridj0,mymaskgridj1;

# #ifdef HAVE_RPNC   
#    int NC_Id,NC_NXDimId,NC_NYDimId;                       ///< netCDF identifiers
# #endif
# } TGeoRef;

from ctypes import *

class TGeoRef(Structure):
    pass

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
        ("zones", POINTER(TGeoZone)), #TODO: array ?
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
   ("RefFrom", POINTER(Subs)),
   ("NbSub", c_int),
   ("Type", c_int),
   ("BD", c_int),
   ("NX", c_int),("NY", c_int),("X0", c_int),("Y0", c_int),("X1", c_int),("Y1", c_int),
   ("Extension", c_int),
   ("GRTYP", c_char_p),
   ("Hemi", c_int),
   ("NbSet", c_int),

   TGridSet *Sets,*LastSet;                               ///< Tableau de set d'interpolation et du dernier utilise

   unsigned int NIdx,*Idx;                                ///< Index dans les positions
   double       *Lat,*Lon;                                ///< Coordonnees des points de grilles (Spherical)
   float        *AX,*AY,*NCX,*NCY,*Hgt;                   ///< Axes de positionnement / deformation
   double       *Wght;                                    ///< Barycentric weight array for TIN  (M grids)

   char                         *String;                  ///< OpenGIS WKT String description
   double                       *Transform,*InvTransform; ///< Transformation functions
   TRotationTransform           *RotTransform;            ///< Rotation transform
   void                         *GCPTransform;            ///< GPC derivative transform (1,2,3 order)
   void                         *TPSTransform;            ///< GPC Thin Spline transform
   void                         *RPCTransform;            ///< GPC Rigorous Projection Model transform
   TQTree                       *QTree;                   ///< Quadtree index
   OGREnvelope                   LLExtent;                ///< LatLon extent
   OGRCoordinateTransformationH  Function,InvFunction;    ///< Projection functions
   OGRSpatialReferenceH          Spatial;                 ///< Spatial reference

   TCoord  Loc;                                           ///< (Radar) Localisation du centre de reference
   double CTH,STH;                                        ///< (Radar) sin and cos of sweep angle
   double ResR,ResA;                                      ///< (Radar) Resolutions en distance et azimuth
   int    R;                                              ///< (Radar) Rayon autour du centre de reference en bin

   TGeoOptions Options;                                   ///< Options for manipulations
   TGeoRef_Project   *Project;                            ///< Transformation xy a latlon
   TGeoRef_UnProject *UnProject;                          ///< Transformation latlon a xy
   TGeoRef_Value     *Value;                              ///< Valeur a un point xy
   TGeoRef_Distance  *Distance;                           
   TGeoRef_Height    *Height;

   _Grille struct from ezscint.h
   TRPNHeader RPNHead;
   int i1, i2, j1, j2;
   int *mask;
   struct TGeoRef *mymaskgrid;
   int mymaskgridi0,mymaskgridi1;
   int mymaskgridj0,mymaskgridj1;

# #ifdef HAVE_RPNC   
#    int NC_Id,NC_NXDimId,NC_NYDimId;                       ///< netCDF identifiers
# #endif
]