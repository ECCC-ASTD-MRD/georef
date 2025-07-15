//! \file

#ifndef _GeoRef_h
#define _GeoRef_h

#include <stdint.h>
#include <string.h>
#include <pthread.h>

#include <rmn/base.h>
#include <rmn/Vector.h>
#include <rmn/QTree.h>

#include "georef/GeoRef_Utils.h"
#include "georef/ZRef.h"

#ifdef HAVE_GDAL
   #include "georef/gdal_safe.h"
   #include <gdal_alg.h>
   #include <ogr_api.h>
   #include <ogr_srs_api.h>
#else
   #include "georef/ogr_stub.h"
#endif

// Extended fst_record information
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
//#define EARTHRADIUS          6378140.0                                                                                   //!< Rayon de la terre en metres
#define EARTHRADIUS          6371000.0                                                                                    //!< Rayon de la terre en metres (Utilise par RPN)

#define DIST(E, A0, O0, A1, O1)  ((E+EARTHRADIUS)*acos(sin(A0)*sin(A1)+cos(O0-O1)*cos(A0)*cos(A1)))                           //!< Calculates great circle distance on earth at a specific elevation between tow set of coordinates (in radian)
#define COURSE(A0, O0, A1, O1)  (fmod(atan2(sin(O0-O1)*cos(A1), cos(A0)*sin(A1)-sin(A0)*cos(A1)*cos(O0-O1)), M_2PI))           //!< Calculates true course between 2 set of coordinates (in radian)
#define M2RAD(M)             ((double)(M)*0.00000015706707756635)                                                         //!< Convert meters to radians
#define M2DEG(M)             ((double)(M)*8.9992806450057884399546578634955e-06)                                          //!< Convert meters to degrees
#define RAD2M(R)             ((double)(R)*6.36670701949370745569e+06)                                                     //!< Convert radians to meters
#define DEG2M(D)             ((double)(D)*1.11119992746859911778451e+05)                                                  //!< Convert degrees to meters
#define CLAMPLAT(LAT)        (LAT=LAT>90.0?90.0:(LAT<-90.0?-90.0:LAT))                                                    //!< Clamp latitude between -90 and 90
#define CLAMPLON(LON)        (LON=LON>180?LON-360:(LON<-180?LON+360:LON))                                                 //!< Clamp longitude between -180 and 180
#define CLAMPLONRAD(LON)     (LON=(LON>M_PI?(fmod(LON+M_PI, M_2PI)-M_PI):(LON<=-M_PI?(fmod(LON-M_PI, M_2PI)+M_PI):LON)))    //!< Clamp longitude in radians -PI and PI
#define SIDELON(SIDE, L)      (((SIDE)>0 && (L)<0)?L+360.0:((SIDE)>0 && (L)<0)?(L)-360.0:(L))                              //!< Force longitude to be all positive or negative
#define COORD_CLEAR(C)       (C.Lat=C.Lon=C.Elev=-999.0)                                                                  //!< Clear the coordinates values to -999 (undefined)

#define DATA_ISVALID(VAL,NODATA) (VAL!=NODATA && VAL==VAL)

#define BITPOS(i)     (i - ((i >> 5) << 5))
#define GETMSK(fld, i) ((fld[i >> 5]  & (1 << BITPOS(i))) >> BITPOS(i))
#define SETMSK(fld, i) (fld[i >> 5] | (fld[i] << BITPOS(i)))
#define NULLOK        (void*)0x1

#define SET_ZONES   0x1          //!< Flag for set zone definitions
#define SET_YYXY    0x4          //!< Flag for set YinYang calculations
#define SET_NZONES  5            //!< Number of zone per set

#define GRID_GLOBAL     0
#define GRID_OUTSIDE    0
#define GRID_NORTH      1
#define GRID_SOUTH      2
#define GRID_NORTH_POLE 3
#define GRID_SOUTH_POLE 4
#define GRID_SUB       ((TGeoRef*)0x1)

#define GRID_NONE     0x0        //!< No flags defined
#define GRID_REGULAR  0x1        //!< Regular grid
#define GRID_VARIABLE 0x2        //!< Variable grid resolution
#define GRID_WRAP     0x4        //!< Wrap around globe
#define GRID_SPARSE   0x8        //!< Sparse grid (point32_t cloud, orca grid
#define GRID_TILE     0x10       //!< Sub tile
#define GRID_VERTICAL 0x20       //!< Vertical grid
#define GRID_RADIAL   0x40       //!< Radar grid
#define GRID_REPEAT   0x80       //!< Does the last longitude repeat the first
#define GRID_PSEUDO   0x100      //!< Pseudocylindrical
#define GRID_ROTATED  0x200      //!< Grid is rotated (North is not up)
#define GRID_NEGLON   0x400      //!< Lon are (-180, 180)
#define GRID_CORNER   0x800      //!< Grid cell is corner defined (ie: first gridpoint32_t is 0.5 0.5)
#define GRID_YINVERT  0x1000     //!< Y axis is inverted
#define GRID_EXPAND   0x2000     //!< Grid needs to be expanded
#define GRID_AXY2D    0x4000     //!< Grid AX/AY are bidimentional

#define GRID_YQTREESIZE   1000   //!< Default Y grid quad tree 2D size
#define GRID_MQTREEDEPTH  8      //!< Default M grid quad tree depth

//#define REF_DEFAULT "GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137.0,298.257222101]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]"
//#define REF_DEFAULT "GEOGCS[\"NAD83",DATUM[\"North_American_Datum_1983\",SPHEROID[\"GRS 1980\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]]"
//#define REF_DEFAULT "GEOGCS[\"WGS84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]]"
#define REF_DEFAULT "GEOGCS[\"GCS_WGS_1984\",DATUM[\"WGS_1984\",SPHEROID[\"WGS84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]]"
#define REF_CLAMPBD(R, PX0, PY0, PX1, PY1) if (PX0<(R->X0+R->BD)) PX0=R->X0+R->BD; if (PY0<(R->Y0+R->BD)) PY0=R->Y0+R->BD; if (PX1>(R->X1-R->BD)) PX1=R->X1-R->BD; if (PY1>(R->Y1-R->BD)) PY1=R->Y1-R->BD;
#define REF_CLAMP(R, PX0, PY0, PX1, PY1)   if (PX0<R->X0) PX0=R->X0; if (PY0<R->Y0) PY0=R->Y0; if (PX1>R->X1) PX1=R->X1; if (PY1>R->Y1) PY1=R->Y1;
#define REF_INDEX_SEPARATOR -1.0
#define REF_INDEX_END       -2.0
#define REF_INDEX_EMPTY     -3.0
#define REF_INDEX_NIL       -999.0

#define GeoRef_ScanX(X)         (((float*)GeoScanX)[X]-1.0)
#define GeoRef_ScanY(X)         (((float*)GeoScanY)[X]-1.0)
#define GeoRef_Lon(R, L)        (((L)>180 && R->Type&GRID_NEGLON)?(L)-360.0:((L)<0 && !(R->Type&GRID_NEGLON))?(L)+360.0:(L))
#define GeoRef_SubGet(REF)      ((REF->Sub < REF->NbSub && REF->Sub >= 0) ? REF->Subs[REF->Sub] : REF)
#define GeoRef_SetHasIndex(G)   (G && G->Index && G->Index[0] != REF_INDEX_EMPTY)
#define GeoRef_SetEmptyIndex(G) (G && G->Index && G->Index[0] == REF_INDEX_EMPTY)

// Interpolation modes
typedef enum {
   IR_UNDEF                          = 0,    //!< Undefined
   IR_NEAREST                        = 1,    //!< Nearest point
   IR_LINEAR                         = 2,    //!< Linear
   IR_CUBIC                          = 3,    //!< Cubic
   IR_NORMALIZED_CONSERVATIVE        = 4,    //!< Mass conservative, distribute the ration of the cell into intersecting grid cells
   IR_CONSERVATIVE                   = 5,    //!< Mass conservative, distribute value on the intersecting grid cells
   IR_MAXIMUM                        = 6,    //!< Use maximum of intersecting cells
   IR_MINIMUM                        = 7,    //!< Use minimum of intersecting cells
   IR_SUM                            = 8,    //!< Sum values of intersecting cells
   IR_AVERAGE                        = 9,    //!< Average value of intersecting cells
   IR_VARIANCE                       = 10,   //!< Variance  of intersecting cells (useful for ...)
   IR_SQUARE                         = 11,   //!< Average of squared values of intersecting cells (useful for ...)
   IR_NORMALIZED_COUNT               = 12,   //!< Normalized % of table (TGeoOptions->Table) specifed values in intersection
   IR_COUNT                          = 13,   //!< Count % of table (TGeoOptions->Table) specifed values in intersection
   IR_VECTOR_AVERAGE                 = 14,   //!< Vectorial direction average
   IR_NOP                            = 15,   //!< Used on multiple interpolsation iertafins
   IR_ACCUM                          = 16,   //!< To get the accumulation matrix of the number of source cell intersecting destination cells
   IR_BUFFER                         = 17,   //!< To get the interpolation state of multipl loops before finalization (ie: destionation cell fraction covered in CONSERVATIVE mode)
   IR_SUBNEAREST                     = 18,   //!< Sub grid resolution nearest interpolation
   IR_SUBLINEAR                      = 19,   //!< Sub grid resolution linear intelpolation 
   IR_WEIGHTINDEX                    = 20,   //!< Pre-calculated weight index
   IV_FAST                           = 21,   //!< Use a rasterization method (middle of grid cell inclusion)
   IV_WITHIN                         = 22,   //!< Grid cell has to be totally included in polygon
   IV_INTERSECT                      = 23,   //!< Grid cell intersects polygon
   IV_CENTROID                       = 24,   //!< Centroid of geometry is within grid cell
   IV_ALIASED                        = 25,   //!< Grid cells are assigned the area fraction of the intersection of the polygon value
   IV_CONSERVATIVE                   = 26,   //!< Distribute polygon value so as to conserve geometry mass in relation to geometry total area and intersecting grid cell area
   IV_NORMALIZED_CONSERVATIVE        = 27,   //!< Distribute polygon coverage fraction in relation to geometry area and intersecting grid cell area
   IV_POINT_CONSERVATIVE             = 28,   //!< Distribute polygon value so as to conserve geometry mass
   IV_LENGTH_CONSERVATIVE            = 29,   //!< Distribute contour value so as to conserve geometry mass in relation to geometry length intersecting grid cell
   IV_LENGTH_NORMALIZED_CONSERVATIVE = 30,   //!< Distribute contour value fraction in relation to geometry length intersecting grid cell
   IV_LENGTH_ALIASED                 = 31    //!< Grid cells are assigned the lenght fraction of the intersection of the contour value
} TRef_Interp;


// Raster Extrapolation modes
typedef enum {
  ER_UNDEF   = 0,     //!< Do nothing (default)
  ER_MAXIMUM = 1,     //!< Use field maximum value
  ER_MINIMUM = 2,     //!< Use minimum field value
  ER_VALUE   = 3,     //!< Use a specific value (TGeoOptions->NoData)
  ER_ABORT   = 4      //!< Abort execution
} TRef_Extrap;

// Interpolation value combination modes (multiple values within a grid cell)
typedef enum {
   CB_REPLACE   = 0,     //!< Replace the value (default)
   CB_MIN       = 1,     //!< Use the minimum
   CB_MAX       = 2,     //!< Use the maximum
   CB_SUM       = 3,     //!< Sum the values
   CB_AVERAGE   = 4,     //!< use the average value
} TRef_Combine;

// Structure pour les coordonees latlon
typedef struct TCoord {
   double Lon, Lat, Elev;
} TCoord;

// Geo vector (grid / geographical coordinates)
typedef union {
   Vect3d V;
   TCoord C;
} GeoVect;

struct TDef;
struct TGeoRef;

int32_t GeoRef_Project(struct TGeoRef *Ref, double X, double Y, double *Lat, double *Lon, int32_t Extrap, int32_t Transform);
int32_t GeoRef_UnProject(struct TGeoRef *Ref, double *X, double *Y, double Lat, double Lon, int32_t Extrap, int32_t Transform);

typedef int32_t (TGeoRef_LL2XY)  (struct TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
typedef int32_t (TGeoRef_XY2LL)  (struct TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);
typedef double  (TGeoRef_Height) (struct TGeoRef *Ref, TZRef *ZRef, double X, double Y, double Z);

// Geospatial manipulation options
typedef struct TGeoOptions {
   TRef_Interp  Interp;         ///< Raster interpolation method (default: IR_CUBIC)
   TRef_Extrap  Extrap;         ///< Raster extrapolation method (default: ER_VALUE)
   TRef_Combine Combine;        //!< Aggregation type (default: CB_REPLACE)
   int32_t      Transform;      //!< Apply transformation or stay within master referential (default: TRUE)
   int32_t      Symmetric;      //!<  (default: FALSE)
   int32_t      Segment;        //!< How much segmentation (Conservatives/Geometric modes) (default: 1)
   int32_t      Sampling;       //!< Sampling interval (default: 1)
   char         PolarCorrect;   //!< Apply polar corrections (default: TRUE)
   char         VectorMode;     //!< Process data as vector (ie: wind) (default: FALSE)
   float        DistTreshold;   //!< Distance treshold for point clouds (default: 10.0)
   float        NoData;         //!< NoData Value (default: NaN)
   double      *Table;          //!< Data table to check of values to check for (default: NULL)
   double     **lutDef;         //!< Lookup table (default: NULL)
   int32_t      lutSize;        //!< Number of lookup elements (multiple lookup values) (default: 0)
   int32_t      lutDim;         //!< Dimension of the lookup elements (default: 0)
   double      *Ancilliary;     //!< Pre calculated field to be passed to interpolation (ex: variance, average, ...) (default: NULL)
} TGeoOptions;

#ifndef GEOREF_BUILD
extern __thread TGeoOptions GeoRef_Options;  //!< Global Options pointer
#endif

typedef struct TRotationTransform {
   double Lat, Lon, Angle, SinTheta, CosTheta, SinPhi, CosPhi;
} TRotationTransform;

typedef struct {
   int32_t npts;             //!< Nombre de points
   double *x, *y;             //!< Vecteur de coordonnees
   int32_t *idx;             //!< Indice du point32_t dans le champ de destination
} TGeoZone;

typedef struct {
   struct TGeoOptions Opt;              //!< Interpolation options
   struct TGeoRef *RefFrom;             //!< Source geo reference
   struct TGeoRef *RefTo;               //!< Destination geo reference
   TGeoZone        zones[SET_NZONES];   //!< Extrapolation zone definitions
   int32_t         flags;               //!< State flags
   TRef_Interp     IndexMethod;         //!< Index interpolation method
   int32_t         IndexSize;           //!< Index size
   float          *Index;               //!< Array of index values
   double         *X, *Y;               //!< Grid coordinates of the destination grid into the source grid

   float *yin_maskout, *yan_maskout;
   double *yinlat, *yinlon, *yanlat, *yanlon;
   double *yin2yin_lat, *yin2yin_lon, *yan2yin_lat, *yan2yin_lon;
   double *yin2yan_lat, *yin2yan_lon, *yan2yan_lat, *yan2yan_lon;
   double *yin2yin_x, *yin2yin_y, *yan2yin_x, *yan2yin_y;
   double *yin2yan_x, *yin2yan_y, *yan2yan_x, *yan2yan_y;
   int32_t yincount_yin, yancount_yin, yincount_yan, yancount_yan;

   int32_t n_wts;                //!< Nombre de poids
   double *wts;                  //!< Tableau de poids
   int32_t *mask, *idx;          //!< Indice du point dans le champ de destination
} TGeoSet;

typedef struct TGeoRef {
   char*        Name;                                     //!< Reference name
   int32_t      NRef;                                     //!< Nombre de reference a la georeference
   struct TGeoRef  *RefFrom;                              //!< Georeference de reference (coupe verticale, ...)
   struct TGeoRef **Subs;                                 //!< Liste des sous grilles (GRTYP=U)
   int32_t      Sub, NbSub;                               //!< Nombre de sous-grilles
   int32_t      Type;                                     //!< Parametre/Flags de grille
   int32_t      BD;                                       //!< Bordure
   int32_t      NX, NY, X0, Y0, X1, Y1;                   //!< Grid size and limits
   int32_t      Extension;                                //!< related to the newtonian coefficient
   char         GRTYP[2];                                 //!< Type de grille
   int32_t      Hemi;                                     //!< Hemisphere side (0=GRID_GLOBAL, 1=NORTH, 2=SOUTH)
   int32_t      NbSet;                                    //!< Nombre de set d'interpolation
   TGeoSet      *Sets, *LastSet;                          //!< Tableau de set d'interpolation et du dernier utilise

   uint32_t     NIdx, *Idx;                               //!< Index dans les positions
   double       *Lat, *Lon;                               //!< Coordonnees des points de grilles (Spherical)
   double       *AX, *AY, *AXY;                           //!< Axes de positionnement / deformation
   float        *NCX, *NCY, *Hgt;
   double       *Wght;                                    //!< Barycentric weight array for TIN  (M grids)

   char                         *String;                  //!< OpenGIS WKT String description
   double                       *Transform, *InvTransform; //!< Transformation functions
   TRotationTransform           *RotTransform;            //!< Rotation transform
   void                         *GCPTransform;            //!< GPC derivative transform (1, 2, 3 order)
   void                         *TPSTransform;            //!< GPC Thin Spline transform
   void                         *RPCTransform;            //!< GPC Rigorous Projection Model transform
   TQTree                       *QTree;                   //!< Quadtree index
   OGREnvelope                   LLExtent;                //!< LatLon extent
   OGRCoordinateTransformationH  Function, InvFunction;   //!< Projection functions
   OGRSpatialReferenceH          Spatial;                 //!< Spatial reference

   TCoord  Loc;                                           //!< (Radar) Localisation du centre de reference
   double  CTH, STH;                                      //!< (Radar) sin and cos of sweep angle
   double  ResR, ResA;                                    //!< (Radar) Resolutions en distance et azimuth
   int32_t R;                                             //!< (Radar) Rayon autour du centre de reference en bin

   TGeoOptions Options;                                   //!< Options for manipulations
   TGeoRef_XY2LL     *XY2LL;                              //!< Transformation xy to latlon
   TGeoRef_LL2XY     *LL2XY;                              //!< Transformation latlon to xy
   TGeoRef_Height    *Height;                             //!< Transformation k to level

   fst_record     RPNHead;                                //!< RPN grid descriptor metadata
   fst_record_ext RPNHeadExt;                             //!< RPN grid descriptor extended metadata

   int32_t i1, i2, j1, j2;
   struct TGeoRef *mymaskgrid;
   int32_t mymaskgridi0, mymaskgridi1;
   int32_t mymaskgridj0, mymaskgridj1;

   pthread_mutex_t Mutex;
} TGeoRef;

typedef struct TGeoPos {
   TGeoRef *GRef;                                         //!< Reference horizontale
   TZRef   *ZRef;                                         //!< Reference verticale
   Vect3d **Pos;                                          //!< Coordonnees des points de grilles (World)
   int32_t  NRef;                                         //!< Nombre de reference a la georeference
} TGeoPos;

typedef struct TGeoScan {
   double   *X, *Y;
   float    *D;
   uint32_t *V;                                           //!< Coordonnees et valeurs
   uint32_t  N, S;                                        //!< Nombre de coordonnees et dimension
   int32_t   DX, DY;                                      //!< Longueur en X et Y
} TGeoScan;

void GeoRef_Lock(void);
void GeoRef_Unlock(void);
void Georef_PrintOptions(TGeoOptions *Options);

TGeoRef* GeoRef_Get(char *Name);
int32_t  GeoRef_Incr(TGeoRef *Ref);
void     GeoRef_Decr(TGeoRef *Ref);
int32_t  GeoRef_Within(TGeoRef* __restrict const Ref0, TGeoRef* __restrict const Ref1);
int32_t  GeoRef_WithinRange(TGeoRef* __restrict const Ref, double Lat0, double Lon0, double Lat1, double Lon1, int32_t In);
int32_t  GeoRef_WithinCell(TGeoRef *Ref, Vect2d Pos, Vect2d Pt[4], int32_t Idx0, int32_t Idx1, int32_t Idx2, int32_t Idx3);
int32_t  GeoRef_Intersect(TGeoRef* __restrict const Ref0, TGeoRef* __restrict const Ref1, int32_t *X0, int32_t *Y0, int32_t *X1, int32_t *Y1, int32_t BD);
int32_t  GeoRef_Equal(TGeoRef* __restrict const Ref0, TGeoRef* __restrict const Ref1);
int32_t  GeoRef_CellDims(TGeoRef *Ref, int32_t Invert, float* DX, float* DY, float* DA);
TGeoRef* GeoRef_New();
TGeoRef* GeoRef_Find(TGeoRef *Ref);
TGeoRef* GeoRef_Copy(TGeoRef* __restrict const Ref);
TGeoRef* GeoRef_HardCopy(TGeoRef* __restrict const Ref);
TGeoRef* GeoRef_Reference(TGeoRef* __restrict const Ref);
void     GeoRef_Size(TGeoRef *Ref, int32_t X0, int32_t Y0, int32_t X1, int32_t Y1, int32_t BD);
TGeoRef* GeoRef_Resize(TGeoRef* __restrict const Ref, int32_t NI, int32_t NJ);
int32_t  GeoRef_Free(TGeoRef *Ref);
void     GeoRef_Clear(TGeoRef *Ref, int32_t New);
void     GeoRef_Qualify(TGeoRef* __restrict const Ref);
int32_t  GeoRef_Limits(TGeoRef* __restrict const Ref, double *Lat0, double *Lon0, double *Lat1, double *Lon1);
int32_t  GeoRef_BoundingBox(TGeoRef* __restrict const Ref, double Lat0, double Lon0, double Lat1, double Lon1, double *I0, double *J0, double *I1, double *J1);
int32_t  GeoRef_Valid(TGeoRef* __restrict const Ref);
int32_t  GeoRef_WriteFST(TGeoRef *GRef,char *Name,int IG1,int IG2,int IG3,int IG4,fst_file *File);
int32_t  GeoRef_Read(struct TGeoRef *GRef);
int32_t  GeoRef_CopyDesc(fst_file *FileTo, fst_record* Rec);
int32_t  GeoRef_AssignToRecord(fst_record *Rec,TGeoRef *Ref);

TGeoRef* GeoRef_CreateFromRecord(fst_record *Rec);
TGeoRef* GeoRef_Create(int32_t NI, int32_t NJ, char *GRTYP, int32_t IG1, int32_t IG2, int32_t IG3, int32_t IG4, fst_file *File);
TGeoRef* GeoRef_CreateU(int32_t NI, int32_t NJ, char *grref, int32_t VerCode, int32_t NbSub, TGeoRef **Subs);
TGeoRef* GeoRef_CreateR(double Lat, double Lon, double Height, int32_t R, double ResR, double ResA);
TGeoRef* GeoRef_CreateW(int32_t ni, int32_t nj, char *String, double *Transform, double *InvTransform, OGRSpatialReferenceH Spatial);
TGeoRef* GeoRef_Define(TGeoRef *Ref, int32_t NI, int32_t NJ, char* GRTYP, char* grref, int32_t IG1, int32_t IG2, int32_t IG3, int32_t IG4, double* AX, double* AY);
TGeoRef* GeoRef_DefineW(TGeoRef *Ref, char *String, double *Transform, double *InvTransform, OGRSpatialReferenceH Spatial);
TGeoRef* GeoRef_DefineZE(TGeoRef *Ref, int32_t NI, int32_t NJ, float DX, float DY, float LatR, float LonR, int32_t MaxCFL, float XLat1, float XLon1, float XLat2, float XLon2);
int32_t  GeoRef_Positional(TGeoRef *Ref, struct TDef *XDef, struct TDef *YDef);
TQTree*  GeoRef_BuildIndex(TGeoRef* __restrict const Ref);
int32_t  GeoRef_Nearest(TGeoRef* __restrict const Ref, double X, double Y, int32_t *Idxs, double *Dists, int32_t NbNear, double MaxDist);

// EZSCINT merged fonctionnalities
int32_t  GeoRef_DefRPNXG(TGeoRef* Ref);                                                                                                       // c_ezdefxg
int32_t  GeoRef_GetLL(TGeoRef *Ref, double *Lat, double *Lon);                                                                                  // gdll
int32_t  GeoRef_CalcLL(TGeoRef* Ref);                                                                                                         // ez_calclatlon
int32_t  GeoRef_XY2LL(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t N, int32_t Extrap);                                     // c_gdllfxy
int32_t  GeoRef_LL2XY(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t N, int32_t Extrap);                                     // c_gdxyfll
double   GeoRef_XYDistance(TGeoRef *Ref, double X0, double Y0, double X1, double Y1);
double   GeoRef_LLDistance(TGeoRef *Ref, double Lat0, double Lon0, double Lat1, double Lon1);

int32_t GeoRef_XYVal(
    const TGeoRef * const Ref,
    const TGeoOptions * const Opt,
    float * const zout,
    const float * const zin,
    double * const X,
    double * const Y,
    const int32_t n
); // c_gdxysval
int32_t GeoRef_XYUVVal(
    const TGeoRef * const Ref,
    const TGeoOptions * const Opt,
    float * const uuout,
    float * const vvout,
    const float * const uuin,
    const float * const vvin,
    double * const X,
    double * const Y,
    const int32_t n
); // c_gdxyvval
int32_t GeoRef_XYWDVal(
    const TGeoRef * const Ref,
    const TGeoOptions * const Opt,
    float * const uuout,
    float * const vvout,
    const float * const uuin,
    const float * const vvin,
    double * const X,
    double * const Y,
    const int32_t n
); // c_gdxywdval
int32_t  GeoRef_LLVal(TGeoRef *Ref, TGeoOptions *Opt, float *zout, float *zin, double *lat, double *lon, int32_t n);                                // c_gdllsval
int32_t  GeoRef_LLUVVal(TGeoRef *Ref, TGeoOptions *Opt, float *uuout, float *vvout, float *uuin, float *vvin, double *Lat, double *Lon, int32_t n);   // c_gdllvval
int32_t  GeoRef_LLWDVal(TGeoRef *Ref, TGeoOptions *Opt, float *uuout, float *vvout, float *uuin, float *vvin, double *Lat, double *Lon, int32_t n);   // c_gdllwdval
int32_t GeoRef_Interp(
    const TGeoRef * const RefTo,
    const TGeoRef * const RefFrom,
    const TGeoOptions * const Opt,
    float * const zout,
    const float * const zin
); // c_ezsint
int32_t GeoRef_InterpWeight(const TGeoRef * const RefTo, const TGeoRef * const RefFrom, const TGeoOptions * const Opt, float * const zout, const float * const zin);
int32_t  GeoRef_InterpUV(TGeoRef *RefTo, TGeoRef *RefFrom, TGeoOptions *Opt, float *uuout, float *vvout, float *uuin, float *vvin);                 // c_ezuvint
int32_t  GeoRef_InterpWD(TGeoRef *RefTo, TGeoRef *RefFrom, TGeoOptions *Opt, float *uuout, float *vvout, float *uuin, float *vvin);                 // c_ezwdint
int32_t GeoRef_InterpYY(
    const TGeoRef * const RefTo,
    const TGeoRef * const RefFrom,
    const TGeoOptions * const Opt,
    float * const zout,
    const float * const zin
); // c_ezyysint
int32_t  GeoRef_InterpYYUV(TGeoRef *RefTo, TGeoRef *RefFrom, TGeoOptions *Opt, float *uuout, float *vvout, float *uuin, float *vvin);               // c_ezyyuvint
int32_t  GeoRef_InterpYYWD(TGeoRef *RefTo, TGeoRef *RefFrom, TGeoOptions *Opt, float *uuout, float *vvout, float *uuin, float *vvin);               // c_ezyywdint
int32_t  GeoRef_InterpMask(TGeoRef *RefTo, TGeoRef *RefFrom, TGeoOptions *Opt, char *MaskOut, char *MaskIn);
int32_t  GeoRef_WD2UV(TGeoRef *Ref, float *uugdout, float *vvgdout, float *uullin, float *vvllin, double *Lat, double *Lon, int32_t npts);           // c_gduvfwd
int32_t  GeoRef_UV2WD(TGeoRef *Ref, float *spd_out, float *wd_out, float *uuin, float *vvin, double *Lat, double *Lon, int32_t npts);                // c_gdwdfuv
int32_t  GeoRef_UV2UV(TGeoRef *Ref, float *uullout, float *vvllout, float *uuin, float *vvin, double *Lat, double *Lon, int32_t Nb);                 // c_gdlluvfuv

int32_t  GeoRef_MaskZones(TGeoRef *RefTo, TGeoRef *RefFrom, int32_t *MaskOut, int32_t *MaskIn);
int32_t  GeoRef_MaskYYApply(TGeoRef *RefTo, TGeoRef *RefFrom, TGeoOptions *Opt, int32_t ni, int32_t nj, float *maskout, double *dlat, double *dlon, double *yinlat, double *yinlon, int32_t *yyincount, double *yanlat, double *yanlon, int32_t *yyancount);
int32_t  GeoRef_MaskYYDefine(TGeoRef *Ref);

void GeoRef_GridGetExpanded(
    const TGeoRef * const Ref,
    const TGeoOptions * const Opt,
    float * const zout,
    const float * const zin
); // gdxpngd
int32_t  GeoRef_GridGetParams(TGeoRef *Ref, int32_t *NI, int32_t *NJ, char *GRTYP, int32_t *IG1, int32_t *IG2, int32_t *IG3, int32_t *IG4, char *grref, int32_t *IG1REF, int32_t *IG2REF, int32_t *IG3REF, int32_t *IG4REF);  //c_ezgxprm
void     GeoRef_AxisDefine(TGeoRef * const Ref, double * const AX, double * const AY);                                                         // gdaxes
int32_t  GeoRef_AxisGetExpanded(const TGeoRef * const Ref, double * const AX, double * const AY);                                              // gdgxpndaxes
void     GeoRef_AxisCalcExpandCoeff(TGeoRef* Ref);
void     GeoRef_AxisCalcNewtonCoeff(TGeoRef* Ref);


// Internal functions
TGeoRef* GeoRef_Add(TGeoRef *Ref);
TGeoSet* GeoRef_SetGet(
    TGeoRef * const RefTo,
    TGeoRef * const RefFrom,
    const TGeoOptions * const Opt
);
void     GeoRef_SetFree(TGeoSet *GSet);
TGeoSet* GeoRef_SetReadFST(
    const TGeoRef * const RefTo,
    const TGeoRef * const RefFrom,
    const int32_t InterpType,
    const fst_file * const File
);
int32_t GeoRef_SetWriteFST(
    const TGeoSet * const GSet,
    fst_file * const File
);
int32_t  GeoRef_SetZoneDefine(TGeoSet *GSet);
int32_t  GeoRef_SetCalcXY(TGeoSet *GSet);
int32_t  GeoRef_SetCalcYYXY(TGeoSet *GSet);
int32_t  GeoRef_SetIndexInit(TGeoSet *GSet);

int32_t GeoRef_InterpFinally(
    TGeoRef * const RefTo,
    TGeoRef * const RefFrom,
    TGeoOptions * const Opt,
    float * const zout,
    const float * const zin,
    double * const X,
    double * const Y,
    const int32_t npts,
    TGeoSet * const GSet
);
int32_t  GeoRef_CorrectValue(TGeoSet *Set, float *zout, float *zin);
int32_t  GeoRef_CorrectVector(TGeoSet *Set, float *uuout, float *vvout, float *uuin, float *vvin);
void     GeoRef_RotateInvertXY(double *Lat, double *Lon, double *X, double *Y, int32_t npts, float xlat1, float xlon1, float xlat2, float xlon2);
void     GeoRef_RotateXY(double *Lat, double *Lon, double *X, double *Y, int32_t npts, float xlat1, float xlon1, float xlat2, float xlon2);


// Per grid type transformations
int32_t  GeoRef_LL2GREF(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_LL2XY_A(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_XY2LL_L(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);
int32_t  GeoRef_LL2XY_B(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_LL2XY_E(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_XY2LL_E(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);
int32_t  GeoRef_LL2XY_L(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_LL2XY_NS(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_XY2LL_NS(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);
int32_t  GeoRef_LL2XY_T(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_XY2LL_T(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);
int32_t  GeoRef_LL2XY_LAMBERT(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_XY2LL_LAMBERT(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);
int32_t  GeoRef_LL2XY_G(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_XY2LL_G(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);
int32_t  GeoRef_LL2XY_Z(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_XY2LL_Z(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);
int32_t  GeoRef_LL2XY_O(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_XY2LL_O(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);
int32_t  GeoRef_LL2XY_R(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_XY2LL_R(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);
int32_t  GeoRef_LL2XY_M(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_XY2LL_M(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);
int32_t  GeoRef_LL2XY_W(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_XY2LL_W(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);
int32_t  GeoRef_LL2XY_Y(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb);
int32_t  GeoRef_XY2LL_Y(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb);

double   GeoFunc_RadialPointRatio(TCoord C1, TCoord C2, TCoord C3);
int32_t  GeoFunc_RadialPointOn(TCoord C1, TCoord C2, TCoord C3, TCoord *CR);
int32_t  GeoFunc_RadialIntersect(TCoord C1, TCoord C2, double CRS13, double CRS23, TCoord *C3);


//! Reproject vector directions geographically
static inline double GeoRef_GeoDir(
    TGeoRef* __restrict const Ref,
    double X,
    double Y
) {
    //! \return Direction?
    if (Ref->GRTYP[0] != 'Y' && Ref->GRTYP[0] != 'M') {
        // Reproject vector orientation by adding grid projection's north difference
        double x[2] = {X, X};
        double y[2] = {Y, Y + 1 / .0};
        double lat[2], lon[2];
        GeoRef_XY2LL(Ref, x, y, &lat[0], &lon[0], 2, 1);

        lat[0] = DEG2RAD(lat[0]); lon[0] = DEG2RAD(lon[0]);
        lat[1] = DEG2RAD(lat[1]); lon[1] = DEG2RAD(lon[1]);
        return COURSE(lat[0], lon[0], lat[1], lon[1]);
    } else {
        return 0.0;
    }
}


//! Grid descriptor search
static inline int32_t GeoRef_XFind(
    const double Value,
    const double * const Table,
    const int32_t Nb,
    const int32_t D
) {
    int32_t start = 0;
    int32_t end = Nb - 1;
    int32_t mid = (start + end) >> 1;

    while(mid != start) {
        if (Value <= Table[mid * D]) {
            end = mid;
        } else {
            start = mid;
        }
        mid = (start + end) >> 1;
    }
    return mid;
}


//! LatLon translation and scaling
static inline void GeoRef_LL2GD(
    const TGeoRef * __restrict const Ref,
    double * const X,
    double * const Y,
    const double * const Lat,
    const double * const Lon,
    const int32_t Nb,
    const float Lat0,
    const float Lon0,
    const float DLat,
    const float DLon
) {
#define CLAMPLONREF(LON, REF) (LON > (360.0 + REF) ? (LON - 360.0) : LON < REF ? (LON + 360.0) : LON)

    const double r = (Ref->AX && Ref->AX[0] < 0.0) ? -180.0 : 0.0;

    for(int32_t i = 0; i < Nb; i++) {
        X[i] = (CLAMPLONREF(Lon[i], r) - Lon0) / DLon;
        Y[i] = (Lat[i] - Lat0) / DLat;
    }
}

#endif
