#ifndef _GeoRef_h
#define _GeoRef_h

#include <stdint.h>
#include <string.h>

#include <App.h>
#include <rmn.h>
#include <rmn/rpnmacros.h>
#include <rmn/Vector.h>
#include <rmn/QTree.h>

#include "GeoRef_Utils.h"
#include "Triangle.h"
#include "ZRef.h"

#ifdef HAVE_GDAL
#include "gdal_safe.h"
#include "gdal_alg.h"
#include "ogr_api.h"
#include "ogr_srs_api.h"
#else
#include "ogr_stub.h"
#endif

typedef struct {
    char grref[FST_GTYP_LEN]; //!< Type of geographical projection
    float xg1;                //!< First grid descriptor
    float xg2;                //!< Second grid descriptor
    float xg3;                //!< Third grid descriptor
    float xg4;                //!< Fourth grid descriptor

    int32_t igref1;           //!< First grid descriptor
    int32_t igref2;           //!< Second grid descriptor
    int32_t igref3;           //!< Third grid descriptor
    int32_t igref4;           //!< Fourth grid descriptor
    float xgref1;             //!< First grid descriptor
    float xgref2;             //!< Second grid descriptor
    float xgref3;             //!< Third grid descriptor
    float xgref4;             //!< Fourth grid descriptor
} fst_record_ext;


// Geographical related constants and functions
//#define EARTHRADIUS          6378140.0                                                                                   ///< Rayon de la terre en metres
#define EARTHRADIUS          6371000.0                                                                                    ///< Rayon de la terre en metres (Utilise par RPN)

#define DIST(E,A0,O0,A1,O1)  ((E+EARTHRADIUS)*acos(sin(A0)*sin(A1)+cos(O0-O1)*cos(A0)*cos(A1)))                           ///< Calculates great circle distance on earth at a specific elevation between tow set of coordinates (in radian)
#define COURSE(A0,O0,A1,O1)  (fmod(atan2(sin(O0-O1)*cos(A1),cos(A0)*sin(A1)-sin(A0)*cos(A1)*cos(O0-O1)),M_2PI))           ///< Calculates true course between 2 set of coordinates (in radian)
#define M2RAD(M)             ((double)(M)*0.00000015706707756635)                                                         ///< Convert meters to radians
#define M2DEG(M)             ((double)(M)*8.9992806450057884399546578634955e-06)                                          ///< Convert meters to degrees
#define RAD2M(R)             ((double)(R)*6.36670701949370745569e+06)                                                     ///< Convert radians to meters
#define DEG2M(D)             ((double)(D)*1.11119992746859911778451e+05)                                                  ///< Convert degrees to meters
#define CLAMPLAT(LAT)        (LAT=LAT>90.0?90.0:(LAT<-90.0?-90.0:LAT))                                                    ///< Clamp latitude between -90 and 90
#define CLAMPLON(LON)        (LON=LON>180?LON-360:(LON<-180?LON+360:LON))                                                 ///< Clamp longitude between -180 and 180
#define CLAMPLONRAD(LON)     (LON=(LON>M_PI?(fmod(LON+M_PI,M_2PI)-M_PI):(LON<=-M_PI?(fmod(LON-M_PI,M_2PI)+M_PI):LON)))    ///< Clamp longitude in radians -PI and PI
#define CLAMPLONREF(LON,REF) (LON>(360.0+REF)?(LON-360.0):LON<REF?(LON+360.0):LON)                                        ///< Clamp longitude between -180 and 180
#define SIDELON(SIDE,L)      (((SIDE)>0 && (L)<0)?L+360.0:((SIDE)>0 && (L)<0)?(L)-360.0:(L))                              ///< Force longitude to be all positive or negative
#define COORD_CLEAR(C)       (C.Lat=C.Lon=C.Elev=-999.0)                                                                  ///< Clear the coordinates values to -999 (undefined)

#define BITPOS(i)     (i - ((i >> 5) << 5))
#define GETMSK(fld,i) ((fld[i >> 5]  & (1 << BITPOS(i))) >> BITPOS(i))
#define SETMSK(fld,i) (fld[i >> 5] | (fld[i] << BITPOS(i)))
#define NULLOK        (void*)0x1

#define SET_ZONES   0x1          ///< Flag for set zone definitions
#define SET_YYXY    0x4          ///< Flag for set YinYang calculations
#define SET_NZONES  5            ///< Number of zone per set
#define SET_MAX     256          ///< Maximum number of sets

#define GLOBAL     0
#define OUTSIDE    0
#define NORTH      1
#define SOUTH      2
#define NORTH_POLE 3
#define SOUTH_POLE 4

#define GRID_NONE     0x0        ///< No flags defined
#define GRID_REGULAR  0x1        ///< Regular grid
#define GRID_VARIABLE 0x2        ///< Variable grid resolution
#define GRID_WRAP     0x4        ///< Wrap around globe
#define GRID_SPARSE   0x8        ///< Sparse grid (point cloud,orca grid
#define GRID_TILE     0x10       ///< Sub tile 
#define GRID_VERTICAL 0x20       ///< Vertical grid
#define GRID_RADIAL   0x40       ///< Radar grid
#define GRID_REPEAT   0x80       ///< Does the last longitude repeat the first
#define GRID_PSEUDO   0x100      ///< Pseudocylindrical
#define GRID_ROTATED  0x200      ///< Grid is rotated (North is not up)
#define GRID_NEGLON   0x400      ///< Lon are (-180,180)
#define GRID_CORNER   0x800      ///< Grid cell is corner defined (ie: first gridpoint is 0.5 0.5)
#define GRID_YINVERT  0x1000     ///< Y axis is inverted
#define GRID_EXPAND   0x2000     ///< Grid needs to be expanded
#define GRID_AXY2D    0x4000     ///< Grid AX/AY are bidimentional

#define GRID_YQTREESIZE   1000   ///< Default Y grid quad tree 2D size
#define GRID_MQTREEDEPTH  8      ///< Default M grid quad tree depth

//#define REF_DEFAULT "GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137.0,298.257222101]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]"
//#define REF_DEFAULT "GEOGCS[\"NAD83",DATUM[\"North_American_Datum_1983\",SPHEROID[\"GRS 1980\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]]"
//#define REF_DEFAULT "GEOGCS[\"WGS84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]]"
#define REF_DEFAULT "GEOGCS[\"GCS_WGS_1984\",DATUM[\"WGS_1984\",SPHEROID[\"WGS84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]]"
#define REF_CLAMPBD(R,PX0,PY0,PX1,PY1) if (PX0<(R->X0+R->BD)) PX0=R->X0+R->BD; if (PY0<(R->Y0+R->BD)) PY0=R->Y0+R->BD; if (PX1>(R->X1-R->BD)) PX1=R->X1-R->BD; if (PY1>(R->Y1-R->BD)) PY1=R->Y1-R->BD;
#define REF_CLAMP(R,PX0,PY0,PX1,PY1)   if (PX0<R->X0) PX0=R->X0; if (PY0<R->Y0) PY0=R->Y0; if (PX1>R->X1) PX1=R->X1; if (PY1>R->Y1) PY1=R->Y1;
#define REF_GET(REF)                   (REF->Options.SubGrid?REF->Subs[REF->Options.SubGrid-1]:REF)
#define REF_INDEX_SEPARATOR -1.0
#define REF_INDEX_END       -2.0
#define REF_INDEX_EMPTY     -3.0

#define REF_COORD(REF,N,C)\
   if (REF->GRTYP[1]!='\0') {\
      REF->Project(REF,REF->Lon[N],REF->Lat[N],&C.lat,&C.lon,1,1);\
   } else {\
      C.lat=REF->Lat[N];\
      C.lon=REF->Lon[N];\
   }

#define REF_TRANSFORM(REF,X,Y,IX,IY)\
   if (REF->Transform) {\
      X=REF->Transform[0]+REF->Transform[1]*(IX)+REF->Transform[2]*(IY);\
      Y=REF->Transform[3]+REF->Transform[4]*(IX)+REF->Transform[5]*(IY);\
   } else {\
      X=IX;\
      Y=IY;\
   }

#define REF_INVTRANSFORM(REF,X,Y,IX,IY)\
   if (REF->InvTransform) {\
      X=REF->InvTransform[0]+REF->InvTransform[1]*(IX)+REF->InvTransform[2]*(IY);\
      Y=REF->InvTransform[3]+REF->InvTransform[4]*(IX)+REF->InvTransform[5]*(IY);\
   } else {\
      X=IX;\
      Y=IY;\
   }

#define GeoRef_ScanX(X) (((float*)GeoScanX)[X]-1.0)
#define GeoRef_ScanY(X) (((float*)GeoScanY)[X]-1.0)
#define GeoRef_Lon(R,L) (((L)>180 && R->Type&GRID_NEGLON)?(L)-360.0:((L)<0 && !(R->Type&GRID_NEGLON))?(L)+360.0:(L))

// Raster interpolation modes
typedef enum {
   IR_UNDEF                          = 0,
   IR_NEAREST                        = 1,
   IR_LINEAR                         = 2,
   IR_CUBIC                          = 3,
   IR_NORMALIZED_CONSERVATIVE        = 4,
   IR_CONSERVATIVE                   = 5,
   IR_MAXIMUM                        = 6,
   IR_MINIMUM                        = 7,
   IR_SUM                            = 8,
   IR_AVERAGE                        = 9,
   IR_VARIANCE                       = 10,
   IR_SQUARE                         = 11,
   IR_NORMALIZED_COUNT               = 12,
   IR_COUNT                          = 13,
   IR_VECTOR_AVERAGE                 = 14,
   IR_NOP                            = 15,
   IR_ACCUM                          = 16,
   IR_BUFFER                         = 17,
   IR_SUBNEAREST                     = 18,
   IR_SUBLINEAR                      = 19
} TRef_InterpR;

// Raster Extrapolation modes
typedef enum {
  ER_UNDEF   = 0,
  ER_MAXIMUM = 1,
  ER_MINIMUM = 2,
  ER_VALUE   = 3,
  ER_ABORT   = 4,
} TRef_ExtrapR;

// Vector interpolation modes
typedef enum {
   IV_UNDEF                          = 0,
   IV_FAST                           = 1,
   IV_WITHIN                         = 2,
   IV_INTERSECT                      = 3,
   IV_CENTROID                       = 4,
   IV_ALIASED                        = 5,
   IV_CONSERVATIVE                   = 6,
   IV_NORMALIZED_CONSERVATIVE        = 7,
   IV_POINT_CONSERVATIVE             = 8,
   IV_LENGTH_CONSERVATIVE            = 9,
   IV_LENGTH_NORMALIZED_CONSERVATIVE = 10,
   IV_LENGTH_ALIASED                 = 11
} TRef_InterpV;

// Interpolation value combination modes
typedef enum {
   CB_REPLACE   = 0,
   CB_MIN       = 1,
   CB_MAX       = 2,
   CB_SUM       = 3,
   CB_AVERAGE   = 4,
} TRef_Combine;

// Structure pour les coordonees latlon
typedef struct TCoord {
   double Lon,Lat,Elev;
} TCoord;

typedef struct TGridCoord {
   float Lat,Lon,I,J;
} TGridCoord;

typedef struct TGridPoint {
   float I,J;
} TGridPoint;

typedef union {
   Vect3d V;
   TCoord C;
} GeoVect;

struct TZRef;
struct TDef;
struct TGeoRef;

typedef int    (TGeoRef_Project)   (struct TGeoRef *Ref,double X,double Y,double *Lat,double *Lon,int Extrap,int Transform);
typedef int    (TGeoRef_UnProject) (struct TGeoRef *Ref,double *X,double *Y,double Lat,double Lon,int Extrap,int Transform);
typedef int    (TGeoRef_LL2XY)     (struct TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb);
typedef int    (TGeoRef_XY2LL)     (struct TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
typedef int    (TGeoRef_Value)     (struct TGeoRef *Ref,struct TDef *Def,TRef_InterpR Interp,int C,double X,double Y,double Z,double *Length,double *ThetaXY);
typedef double (TGeoRef_Distance)  (struct TGeoRef *Ref,double X0,double Y0,double X1, double Y1);
typedef double (TGeoRef_Height)    (struct TGeoRef *Ref,TZRef *ZRef,double X,double Y,double Z);

#define GOPT_MAXWEIGHTNUM 256

// Geospatial manipulation options
typedef struct TGeoOptions {
   TRef_InterpR InterpDegree;   ///< Interpolation degree
   TRef_ExtrapR ExtrapDegree;   ///< Extrapolation method
   TRef_InterpV InterpVector;   ///< Vector interpolation method
   TRef_Combine Combine;        ///< Aggrgation type
   double       ExtrapValue;    ///< Value to use for extrapolation in ER_VALUE mode
   int          SubGrid;        ///< Subgrid to use (0=all)
   int          Transform;      ///< Apply transformation or stay within master referential
   int          CIndex;         ///< C Indexing (starts st 0)
   int          Symmetric;      ///< 
   int          WeightNum;      ///<
   int          Segment;        ///< How much segmentation (Conservatives/Geometric modes)
   int          Sampling;       ///< Sampling interval
   char         PolarCorrect;   ///< Apply polar corrections
   char         VectorMode;     ///< Process data as vector
   float        DistTreshold;   ///< Distance treshold for point clouds
   float        LonRef;         ///< Longitude referential (-180.0,0.0)
   float        NoData;         ///< NoData Value (Default: NaN)
} TGeoOptions;

#ifndef GEOREF_BUILD
extern __thread TGeoOptions GeoRef_Options;       ///< Global Options pointer
#endif

typedef struct TRotationTransform {
   double Lat,Lon,Angle,SinTheta,CosTheta,SinPhi,CosPhi;   
} TRotationTransform;

typedef struct {
   int npts;                 ///< Nombre de points
   double *x,*y;             ///< Vecteur de coordonnees 
   int *idx;                 ///< Indice du point dans le champ de destination
} TGeoZone;

typedef struct {
   struct TGeoRef *RefFrom,*RefTo;
   TGeoZone        zones[SET_NZONES];
   char            G2G[2];                          ///<GRTYP of source and destination for index identification
   int             flags;
   TRef_InterpR    IndexDegree;
   int             IndexSize;
   float          *Index;
   double         *X,*Y;

   float *yin_maskout,*yan_maskout;
   double *yinlat,*yinlon,*yanlat,*yanlon;
   double *yin2yin_lat,*yin2yin_lon,*yan2yin_lat,*yan2yin_lon;
   double *yin2yan_lat,*yin2yan_lon,*yan2yan_lat,*yan2yan_lon;
   double *yin2yin_x,*yin2yin_y,*yan2yin_x,*yan2yin_y;
   double *yin2yan_x,*yin2yan_y,*yan2yan_x,*yan2yan_y;
   int yincount_yin,yancount_yin,yincount_yan,yancount_yan;

   int n_wts;                ///< Nombre de poids
   double *wts;              ///< Tableau de poids
   int *mask, *idx;          ///< Indice du point dans le champ de destination
} TGridSet;

#define GeoRef_SetHasIndex(G)   (G && G->Index && G->Index[0]!=REF_INDEX_EMPTY)
#define GeoRef_SetEmptyIndex(G) (G && G->Index && G->Index[0]==REF_INDEX_EMPTY)

typedef struct TGeoRef {
   char*   Name;                                          ///< Reference name
   int     NRef;                                          ///< Nombre de reference a la georeference
   struct  TGeoRef  *RefFrom;                             ///< Georeference de reference (coupe verticale,...)
   struct  TGeoRef **Subs;                                ///< Liste des sous grilles (GRTYP=U)
   int     NbSub;                                         ///< Nombre de sous-grilles
   int     Type;                                          ///< Parametre/Flags de grille
   int     BD;                                            ///< Bordure
   int     NX,NY,X0,Y0,X1,Y1;                             ///< Grid size and limits
   int     Extension;                                     ///< related to the newtonian coefficient
   char    GRTYP[3];                                      ///< Type de grille
   int     Hemi;                                          ///< Hemisphere side (0=GLOBAL,1=NORTH,2=SOUTH)
   int     NbSet;                                         ///< Nombre de set d'interpolation
   TGridSet *Sets,*LastSet;                               ///< Tableau de set d'interpolation et du dernier utilise

   unsigned int NIdx,*Idx;                                ///< Index dans les positions
   double       *Lat,*Lon;                                ///< Coordonnees des points de grilles (Spherical)
   double       *AX,*AY;                                  ///< Axes de positionnement / deformation
   float        *NCX,*NCY,*Hgt;                   
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
   TGeoRef_XY2LL     *XY2LL;                              ///< Transformation xy a latlon
   TGeoRef_LL2XY     *LL2XY;                              ///< Transformation latlon a xy
   TGeoRef_Height    *Height;

   // _Grille struct from ezscint.h
   fst_record     RPNHead;
   fst_record_ext RPNHeadExt;

   int i1, i2, j1, j2;
   struct TGeoRef *mymaskgrid;
   int mymaskgridi0,mymaskgridi1;
   int mymaskgridj0,mymaskgridj1;

#ifdef HAVE_RPNC   
   int NC_Id,NC_NXDimId,NC_NYDimId;                       ///< netCDF identifiers
#endif
} TGeoRef;

typedef struct TGeoPos {
   TGeoRef *Ref;                                          ///< Reference horizontale
   TZRef   *ZRef;                                         ///< Reference verticale
   Vect3d **Pos;                                          ///< Coordonnees des points de grilles (World)
   int      NRef;                                         ///< Nombre de reference a la georeference
} TGeoPos;

typedef struct TGeoScan {
   double *X,*Y;
   float  *D;
   unsigned int *V;                                       ///< Coordonnees et valeurs
   unsigned int N,S;                                      ///< Nombre de coordonnees et dimension
   int DX,DY;                                             ///< Longueur em X et Y
} TGeoScan;

void GeoRef_Lock(void);
void GeoRef_Unlock(void);

TGeoRef* GeoRef_Get(char *Name);
int      GeoRef_Incr(TGeoRef *Ref);
void     GeoRef_Decr(TGeoRef *Ref);
int      GeoRef_Within(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1);
int      GeoRef_WithinRange(TGeoRef* __restrict const Ref,double Lat0,double Lon0,double Lat1,double Lon1,int In);
int      GeoRef_WithinCell(TGeoRef *Ref,Vect2d Pos,Vect2d Pt[4],int Idx0,int Idx1,int Idx2,int Idx3);
int      GeoRef_Intersect(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1,int *X0,int *Y0,int *X1,int *Y1,int BD);
int      GeoRef_Equal(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1);
int      GeoRef_CellDims(TGeoRef *Ref,int Invert,float* DX,float* DY,float* DA);
TGeoRef* GeoRef_New();
TGeoRef* GeoRef_Add(TGeoRef *Ref);
TGeoRef* GeoRef_Find(TGeoRef *Ref);
TGeoRef* GeoRef_Copy(TGeoRef* __restrict const Ref);
TGeoRef* GeoRef_HardCopy(TGeoRef* __restrict const Ref);
TGeoRef* GeoRef_Reference(TGeoRef* __restrict const Ref);
void     GeoRef_Size(TGeoRef *Ref,int X0,int Y0,int X1,int Y1,int BD);
TGeoRef* GeoRef_Resize(TGeoRef* __restrict const Ref,int NI,int NJ);
int      GeoRef_Free(TGeoRef *Ref);
void     GeoRef_Clear(TGeoRef *Ref,int New);
void     GeoRef_Qualify(TGeoRef* __restrict const Ref);
int      GeoRef_Limits(TGeoRef* __restrict const Ref,double *Lat0,double *Lon0,double *Lat1,double *Lon1);
int      GeoRef_BoundingBox(TGeoRef* __restrict const Ref,double Lat0,double Lon0,double Lat1,double Lon1,double *I0,double *J0,double *I1,double *J1);
int      GeoRef_Valid(TGeoRef* __restrict const Ref);
double   GeoRef_GridDistance(TGeoRef *Ref,double X0,double Y0,double X1,double Y1);
int      GeoRef_Write(TGeoRef *GRef,fst_file *File);

TGeoRef* GeoRef_Create(int NI,int NJ,char *GRTYP,int IG1,int IG2,int IG3,int IG4,fst_file *File);
TGeoRef* GeoRef_CreateU(int NI,int NJ,char *GRTYP,char *grref,int VerCode,int NbSub,TGeoRef **Subs);
TGeoRef* GeoRef_CreateR(double Lat,double Lon,double Height,int R,double ResR,double ResA);
TGeoRef* GeoRef_CreateW(int ni,int nj,char *grtyp,int ig1,int ig2,int ig3,int ig4,char *String,double *Transform,double *InvTransform,OGRSpatialReferenceH Spatial);
TGeoRef* GeoRef_Define(TGeoRef *Ref,int NI,int NJ,char* GRTYP,char* grref,int IG1,int IG2,int IG3,int IG4,double* AX,double* AY);
TGeoRef* GeoRef_DefineW(TGeoRef *Ref,char *String,double *Transform,double *InvTransform,OGRSpatialReferenceH Spatial);
TGeoRef* GeoRef_DefineZE(TGeoRef *Ref,int NI,int NJ,float DX,float DY,float LatR,float LonR,int MaxCFL,float XLat1,float XLon1,float XLat2,float XLon2);
int      GeoRef_Positional(TGeoRef *Ref,struct TDef *XDef,struct TDef *YDef);
TQTree*  GeoRef_BuildIndex(TGeoRef* __restrict const Ref);
int      GeoRef_Nearest(TGeoRef* __restrict const Ref,double X,double Y,int *Idxs,double *Dists,int NbNear,double MaxDist);

// EZSCINT merged fonctionnalities
int      GeoRef_DefRPNXG(TGeoRef* Ref);                                                                                           // c_ezdefxg
int      GeoRef_GetLL(TGeoRef *Ref,double *Lat,double *Lon);                                                                      // gdll / GeoRef_Coords
int      GeoRef_CalcLL(TGeoRef* Ref);                                                                                             // ez_calclatlon
int      GeoRef_XY2LL(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int N,int Extrap);                                 // c_gdllfxy
int      GeoRef_LL2XY(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int N,int Extrap);                                 // c_gdxyfll
int      GeoRef_XY2XY(TGeoRef *RefTo,TGeoRef *RefFrom,double *X1,double *Y1,double *X0,double *Y0,int Nb,int Extrap);
int      GeoRef_XYInterp(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout,float *zin,double *x,double *y,int npts);                    // c_gdxysint
int      GeoRef_XYVal(TGeoRef *Ref,float *zout,float *zin,double *x,double *y,int n);                                             // c_gdxysval
int      GeoRef_XYUVVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,double *x,double *y,int n);                // c_gdxyvval
int      GeoRef_XYWDVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,double *x,double *y,int n);                // c_gdxywdval
int      GeoRef_LLVal(TGeoRef *Ref,float *zout,float *zin,double *lat,double *lon,int n);                                         // c_gdllsval
int      GeoRef_LLUVVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int n);            // c_gdllvval
int      GeoRef_LLWDVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int n);            // c_gdllwdval
int      GeoRef_Interp(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout, float *zin);                                                  // c_ezsint
int      GeoRef_InterpUV(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin);                      // c_ezuvint
int      GeoRef_InterpWD(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin);                      // c_ezwdint
int      GeoRef_InterpYY(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout, float *zin);                                                // c_ezyysint
int      GeoRef_InterpYYUV(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin);                    // c_ezyyuvint
int      GeoRef_InterpYYWD(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin);                    // c_ezyywdint
int      GeoRef_InterpMask(TGeoRef *RefTo, TGeoRef *RefFrom,char *MaskOut,char *MaskIn);
int      GeoRef_WD2UV(TGeoRef *Ref,float *uugdout,float *vvgdout,float *uullin,float *vvllin,double *Lat,double *Lon,int npts);   // c_gduvfwd
int      GeoRef_UV2WD(TGeoRef *Ref,float *spd_out,float *wd_out,float *uuin,float *vvin,double *Lat,double *Lon,int npts);        // c_gdwdfuv

int      GeoRef_MaskZones(TGeoRef *RefTo,TGeoRef *RefFrom,int *MaskOut,int *MaskIn);
int      GeoRef_MaskYYApply(TGeoRef *RefTo,TGeoRef *RefFrom,int ni,int nj,float *maskout,double *dlat,double *dlon,double *yinlat,double *yinlon,int *yyincount,double *yanlat,double *yanlon,int *yyancount);
int      GeoRef_MaskYYDefine(TGeoRef *Ref);

void     GeoRef_GridGetExpanded(TGeoRef *Ref, float *zout, float *zin);                                                           // gdxpngd
int      GeoRef_GridGetParams(TGeoRef *Ref,int *NI,int *NJ,char *GRTYP,int *IG1,int *IG2,int *IG3,int *IG4,char *grref,int *IG1REF,int *IG2REF,int *IG3REF,int *IG4REF);  //c_ezgxprm
int      GeoRef_AxisGet(TGeoRef *Ref,double *AX,double *AY);                                                                        // gdaxes
int      GeoRef_AxisGetExpanded(TGeoRef* Ref,double *AX,double *AY);                                                              // gdgxpndaxes
void     GeoRef_AxisDefine(TGeoRef* Ref,double *AX,double *AY);
void     GeoRef_AxisCalcExpandCoeff(TGeoRef* Ref);
void     GeoRef_AxisCalcNewtonCoeff(TGeoRef* Ref);

int32_t c_gdlluvfuv_orig(int32_t gdid,float *uullout,float *vvllout,float *uuin,float *vvin,float *latin,float *lonin,int32_t npts);
void  c_ezllwfgff(float *uullout,float *vvllout,float *latin,float *lonin,float *xlatingf,float *xloningf,int32_t *ni,int32_t *nj,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4);

// Internal functions
TGridSet* GeoRef_SetGet(TGeoRef* RefTo, TGeoRef* RefFrom);
void      GeoRef_SetFree(TGridSet *GSet);
TGridSet* GeoRef_SetAdd(TGeoRef* RefTo,TGeoRef* RefFrom,int InterpType,fst_file *File);
int       GeoRef_SetWrite(TGridSet *GSet,fst_file *File);
int       GeoRef_SetZoneDefine(TGridSet *GSet);
int       GeoRef_SetCalcXY(TGridSet *GSet);
int       GeoRef_SetCalcYYXY(TGridSet *GSet);

int       GeoRef_InterpFinally(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout,float *zin,double *x,double *y,int npts,TGridSet *GSet);
int       GeoRef_CorrectValue(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout, float *zin);
int       GeoRef_CorrectVector(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout, float *vvout, float *uuin, float *vvin);
void      GeoRef_RotateInvertXY(double *Lat,double *Lon,double *X,double *Y,int npts,float xlat1,float xlon1,float xlat2,float xlon2);
void      GeoRef_RotateXY(double *Lat,double *Lon,double *X,double *Y,int npts,float xlat1,float xlon1,float xlat2,float xlon2);

// Per grid type transformations
int       GeoRef_LL2GREF(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb);
int       GeoRef_LL2XY_A(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb); 
int       GeoRef_XY2LL_L(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
int       GeoRef_LL2XY_B(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb);
int       GeoRef_LL2XY_E(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb); 
int       GeoRef_XY2LL_E(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
int       GeoRef_LL2XY_L(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb);
int       GeoRef_LL2XY_NS(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb); 
int       GeoRef_XY2LL_NS(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
int       GeoRef_LL2XY_T(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb); 
int       GeoRef_XY2LL_T(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
int       GeoRef_LL2XY_LAMBERT(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb); 
int       GeoRef_XY2LL_LAMBERT(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
int       GeoRef_LL2XY_G(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb); 
int       GeoRef_XY2LL_G(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
int       GeoRef_LL2XY_Z(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb); 
int       GeoRef_XY2LL_Z(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
int       GeoRef_LL2XY_O(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb); 
int       GeoRef_XY2LL_O(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
int       GeoRef_LL2XY_R(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb); 
int       GeoRef_XY2LL_R(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
int       GeoRef_LL2XY_M(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb); 
int       GeoRef_XY2LL_M(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
int       GeoRef_LL2XY_W(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb);
int       GeoRef_XY2LL_W(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);
int       GeoRef_LL2XY_Y(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb); 
int       GeoRef_XY2LL_Y(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb);

void GeoScan_Init(TGeoScan *Scan);
void GeoScan_Clear(TGeoScan *Scan);
int  GeoScan_Get(TGeoScan *Scan,TGeoRef *ToRef,struct TDef *ToDef,TGeoRef *FromRef,struct TDef *FromDef,int X0,int Y0,int X1,int Y1,int Dim);

double GeoFunc_RadialPointRatio(TCoord C1,TCoord C2,TCoord C3);
int    GeoFunc_RadialPointOn(TCoord C1,TCoord C2,TCoord C3,TCoord *CR);
int    GeoFunc_RadialIntersect(TCoord C1,TCoord C2,double CRS13,double CRS23,TCoord *C3);

// Reproject vector directions geographically
static inline double GeoRef_GeoDir(TGeoRef* __restrict const Ref,double X, double Y) {
   
   double lat[2],lon[2],x[2],y[2],dir=0.0;
   
   if (Ref->GRTYP[0]!='Y' && Ref->GRTYP[0]!='M') {
      
      // Reproject vector orientation by adding grid projection's north difference
      x[0]=x[1]=X;
      y[0]=Y;y[1]=Y+1/.0;
      GeoRef_XY2LL(Ref,x,y,&lat[0],&lon[0],2,TRUE);

      lat[0]=DEG2RAD(lat[0]); lon[0]=DEG2RAD(lon[0]);
      lat[1]=DEG2RAD(lat[1]); lon[1]=DEG2RAD(lon[1]);
      dir=COURSE(lat[0],lon[0],lat[1],lon[1]);
   }
   return(dir);
}

// Grid descriptor search
static inline int GeoRef_XFind(double Value,double *Table,int Nb,int D) {

   int start,mid,end;

	start=0;
	end=Nb-1;
	mid=(start+end)>>1;

   while(mid!=start) {
      if (Value <= Table[mid*D]) {
	      end   = mid;
      } else {
	      start = mid;
      }
      mid = (start+end)>>1;
   }
	return(mid);
}

// LatLon translation and scaling
static inline void GeoRef_LL2GD(double *X,double *Y,double *Lat,double *Lon,int Nb,float Lat0,float Lon0,float DLat,float DLon,float LonRef) {

   int i;
   
   for(i=0;i<Nb;i++) {
      X[i] = (CLAMPLONREF(Lon[i],LonRef)-Lon0)/DLon + 1.0;
      Y[i] = (Lat[i]-Lat0)/DLat + 1.0;
   }
}
void f77name(ez8_rgd_index_0)(float *index,double *px,double *py,int32_t *npts,int32_t *ni,int32_t *j1,int32_t *j2);
void f77name(ez8_rgd_index_1)(float *index,double *px,double *py,int32_t *npts,int32_t *ni,int32_t *j1,int32_t *j2,int32_t *wrap);
void f77name(ez8_rgd_index_3)(float *index,double *px,double *py,int32_t *npts,int32_t *ni,int32_t *j1,int32_t *j2,int32_t *wrap);
void f77name(ez8_irgd_index_1)(float *index,double *px,double *py,int32_t *npts,int32_t *ni,int32_t *j1,int32_t *j2,double *ax,double *ay,int32_t *wrap);
void f77name(ez8_apply_0)(float *index,float *zo,int32_t *npts,float *z,int32_t *ni,int32_t *j1,int32_t *j2,float *nodata);
void f77name(ez8_apply_1)(float *index,float *zo,int32_t *npts,float *z,int32_t *ni,int32_t *j1,int32_t *j2,float *nodata);
void f77name(ez8_apply_3)(float *index,float *zo,int32_t *npts,float *z,int32_t *ni,int32_t *j1,int32_t *j2,float *nodata);

void f77name(ez8_rgdint_0)      (float *zo,double *px,double *py,int32_t *npts,float *z,int32_t *ni,int32_t *j1,int32_t *j2,float *nodata);
void f77name(ez8_rgdint_1)      (float *zo,double *px,double *py,int32_t *npts,float *z,int32_t *ni,int32_t *j1,int32_t *j2,int32_t *wrap,float *nodata);
void f77name(ez8_rgdint_3)      (float *zo,double *px,double *py,int32_t *npts,float *z,int32_t *ni,int32_t *j1,int32_t *j2,int32_t *wrap,float *nodata);
void f77name(ez8_irgdint_1)     (float *zo,double *px,double *py,int32_t *npts,float *z,int32_t *ni,int32_t *j1,int32_t *j2,double *ax,double *ay,int32_t *wrap,float *nodata);
void f77name(ez8_irgdint_3)     (float *zo,double *px,double *py,int32_t *npts,float *z,int32_t *ni,int32_t *i1,int32_t *i2,int32_t *j1,int32_t *j2,double *ax,double *ay,float *cx,float *cy,int32_t *wrap,float *nodata);
void f77name(ez8_irgdint_3_wnnc)(float *zo,double *px,double *py,int32_t *npts,float *z,int32_t *ni,int32_t *j1,int32_t *j2,double *ax,double *ay,int32_t *wrap);

void f77name(ez8_llwfgdw)(float *z1,float *z2,double *xlon,int32_t *li,int32_t *lj,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4,int32_t sz);
void f77name(ez8_gdwfllw)(float *z1,float *z2,double *xlon,int32_t *li,int32_t *lj,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4,int32_t sz);
void f77name(ez8_vllfxy)(double *dlat,double *dlon,double *x,double *y,int32_t *ni,int32_t *nj,float *d60,float *dgrw,float *pi,float *pj,int32_t *nhem);
void f77name(ez8_vtllfxy)(double *lat,double *lon,double *x,double *y,float *clat,float *clon,float *d60,float *dgrw,int32_t *ni,int32_t *nj,int32_t *n);
void f77name(ez8_vxyfll)(double *x,double *y,double *dlat,double *dlon,int32_t *nb,float *d60,float *dgrw,float *pi,float *pj,int32_t *nhem);
void f77name(ez8_vtxyfll)(double *x,double *y,double *dlat,double *dlon,float *clat,float *clon,float *d60,float *dgrw,int32_t *ni,int32_t *nj,int32_t *nhem);
void f77name(ez8_mxm)(float *a,int *nar,double *b,int *nac,double *c,int *nbc);
void f77name(ez8_llflamb)(double *xlat,double *xlon,double *x,double *y,int32_t *npts,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4);
void f77name(ez8_lambfll)(double *x,double *y,double *xlat,double *xlon,int32_t *npts,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4);
void f77name(ez_lambxyfll99)(double *x,double *y,double *lat,double *lon,int32_t *nb,float *latin1,float *latin2,float *yaxislat,float *yaxislon);
void f77name(ez8_lambllfxy99)(double *lat,double *lon,double *x,double *y,int32_t *nb,float *latin1,float *latin2,float *yaxislat,float *yaxislon);
void f77name(ez8_nwtncof)(float *cx,float *cy,double *ax,double *ay,int32_t *ni,int32_t *nj,int32_t *i1,int32_t *i2,int32_t *j1,int32_t *j2,int32_t *extension);
void f77name(igaxg95)(char *gtypout,float *xg,int *nb,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4);
void f77name(ez_avg)();
void f77name(ez_avg_sph)();
void f77name(ez_applywgts)();
void f77name(ez_calcpoleval)(float *poleval,float *z,int32_t *ni,double *ax,char *grtyp,char *grref,int32_t c1,int32_t c2);
void f77name(ez_fillnpole)();
void f77name(ez_fillspole)();
void f77name(ez_calcxy_y)();
void f77name(ez_calcxy_y_m)();
void f77name(ez_aminmax)();
void f77name(ez_corrbgd)();
void f77name(ez_glat)();
void f77name(ez_crot)();
void f77name(ez_lac)();
void f77name(ez8_uvacart)(double *XYZ,float *U,float *V,double *LON,double *LAT,int32_t *NI,int32_t *NJ);
void f77name(ez8_cartauv)(float *U,float *V,double *UVCART,double *LON,double *LAT,int32_t *NI,int32_t *NJ);
void f77name(ez_xpngdb2)();
void f77name(ez_xpngdag2)();
void f77name(permut)();
void f77name(lorenzo_mask_fill)();
void f77name(qqq_ezsint_mask)();
void f77name(qqq_ezget_mask_zones)();

#endif
