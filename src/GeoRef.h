/* =========================================================
 * Environnement Canada
 * Centre Meteorologique Canadien
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Lecture et traitements de fichiers raster
 * Fichier      : GeoRef.h
 * Creation     : Mars 2005 - J.P. Gauthier
 *
 * Description  : Fonctions de manipulations de projections.
 *
 * Remarques    :
 *
 * License      :
 *    This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation,
 *    version 2.1 of the License.
 *
 *    This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the
 *    Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 *    Boston, MA 02111-1307, USA.
 *
 *=========================================================
 */

#ifndef _GeoRef_h
#define _GeoRef_h

#include <stdint.h>
#include "eerUtils.h"
#include "Vector.h"
#include "Triangle.h"
#include "ZRef.h"
#include "QTree.h"

#ifdef HAVE_RMN
#include "rpnmacros.h"
#endif

#ifdef HAVE_GDAL
#include "gdal_safe.h"
#include "gdal_alg.h"
#include "ogr_api.h"
#include "ogr_srs_api.h"
#else
#include "ogr_stub.h"
#endif

#define NZONES             5

#define GRID_NONE     0x0
#define GRID_REGULAR  0x1        // Regular grid
#define GRID_VARIABLE 0x2        // Variable grid resolution
#define GRID_WRAP     0x4        // Wrap around globe
#define GRID_SPARSE   0x8        // Sparse grid (point cloud,orca grid
#define GRID_TILE     0x10       // Sub tile 
#define GRID_VERTICAL 0x20       // Vertical grid
#define GRID_RADIAL   0x40       // Radar grid
#define GRID_REPEAT   0x80       // Does the last longitude repeat the first
#define GRID_PSEUDO   0x100      // Pseudocylindrical
#define GRID_NUNORTH  0x200      // North is not up
#define GRID_NEGLON   0x400      // Lon are (-180,180)
#define GRID_CORNER   0x800      // Grid cell is corner defined (ie: first gridpoint is 0.5 0.5)

#define GRID_YQTREESIZE   1000
#define GRID_MQTREEDEPTH  8

//#define REFDEFAULT "GEOGCS[\"GCS_North_American_1983\",DATUM[\"D_North_American_1983\",SPHEROID[\"GRS_1980\",6378137.0,298.257222101]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]"
//#define REFDEFAULT "GEOGCS[\"NAD83",DATUM[\"North_American_Datum_1983\",SPHEROID[\"GRS 1980\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]]"
//#define REFDEFAULT "GEOGCS[\"WGS84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]]"
#define REFDEFAULT "GEOGCS[\"GCS_WGS_1984\",DATUM[\"WGS_1984\",SPHEROID[\"WGS84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]]"
#define REFCLAMPBD(R,PX0,PY0,PX1,PY1) if (PX0<(R->X0+R->BD)) PX0=R->X0+R->BD; if (PY0<(R->Y0+R->BD)) PY0=R->Y0+R->BD; if (PX1>(R->X1-R->BD)) PX1=R->X1-R->BD; if (PY1>(R->Y1-R->BD)) PY1=R->Y1-R->BD;
#define REFCLAMP(R,PX0,PY0,PX1,PY1)   if (PX0<R->X0) PX0=R->X0; if (PY0<R->Y0) PY0=R->Y0; if (PX1>R->X1) PX1=R->X1; if (PY1>R->Y1) PY1=R->Y1;

#define REFCOORD(REF,N,C)\
   if (REF->Grid[1]!='\0') {\
      REF->Project(REF,REF->Lon[N],REF->Lat[N],&C.lat,&C.lon,1,1);\
   } else {\
      C.lat=REF->Lat[N];\
      C.lon=REF->Lon[N];\
   }

#define TRANSFORM(REF,X,Y,IX,IY)\
   if (REF->Transform) {\
      X=REF->Transform[0]+REF->Transform[1]*(IX)+REF->Transform[2]*(IY);\
      Y=REF->Transform[3]+REF->Transform[4]*(IX)+REF->Transform[5]*(IY);\
   } else {\
      X=IX;\
      Y=IY;\
   }

#define INVTRANSFORM(REF,X,Y,IX,IY)\
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

// Structure pour les coordonees latlon
typedef struct Coord {
   double Lon,Lat,Elev;
} Coord;

typedef union {
   Vect3d V;
   Coord  C;
} GeoVect;

struct TZRef;
struct TDef;
struct TGeoRef;
struct _zone;
struct _ygrid;
struct _gemgrid;
struct _gridset;

typedef int    (TGeoRef_Project)   (struct TGeoRef *Ref,double X,double Y,double *Lat,double *Lon,int Extrap,int Transform);
typedef int    (TGeoRef_UnProject) (struct TGeoRef *Ref,double *X,double *Y,double Lat,double Lon,int Extrap,int Transform);
typedef int    (TGeoRef_Value)     (struct TGeoRef *Ref,struct TDef *Def,char Mode,int C,double X,double Y,double Z,double *Length,double *ThetaXY);
typedef double (TGeoRef_Distance)  (struct TGeoRef *Ref,double X0,double Y0,double X1, double Y1);
typedef double (TGeoRef_Height)    (struct TGeoRef *Ref,TZRef *ZRef,double X,double Y,double Z);

typedef struct TRotationTransform {
   double Lat,Lon,Angle,SinTheta,CosTheta,SinPhi,CosPhi;   
} TRotationTransform;

typedef struct {
  wordint npts;               /* nombre de points */
  ftnfloat *x;                /* vecteur de coordonnees x */
  ftnfloat *y;                /* vecteur de coordonnees y */
  wordint *idx;               /* indice du point dans le champ de destination */
} _zone;

typedef struct {
  wordint n_wts;              /* nombre de poids */
  ftnfloat *xx, *yy;
  ftnfloat *lat, *lon;
  ftnfloat *wts;              /* tableau de poids */
   wordint *mask, *idx;       /* indice du point dans le champ de destination */
} _ygrid;                     /* Grille Y */

typedef struct {
  wordint flags;
  ftnfloat *lat_rot, *lon_rot, *lat_true, *lon_true;
  ftnfloat *sinlat_rot, *coslat_rot, *sinlon_rot, *coslon_rot;
  ftnfloat *sinlat_true, *coslat_true, *sinlon_true, *coslon_true;
  ftnfloat r[9], ri[9];
} _gemgrid;

typedef struct {
  wordint flags,yyflags;
  wordint use_sincos_cache;
  wordint gdin;
  ftnfloat valpolesud, valpolenord;
  ftnfloat *x, *y;
  wordint *mask_in, *mask_out;
  ftnfloat *yin_maskout,*yan_maskout;
  ftnfloat *yinlat,*yinlon,*yanlat,*yanlon;
  ftnfloat *yin2yin_lat,*yin2yin_lon,*yan2yin_lat,*yan2yin_lon;
  ftnfloat *yin2yan_lat,*yin2yan_lon,*yan2yan_lat,*yan2yan_lon;
  ftnfloat *yin2yin_x,*yin2yin_y,*yan2yin_x,*yan2yin_y;
  ftnfloat *yin2yan_x,*yin2yan_y,*yan2yan_x,*yan2yan_y;
  wordint yincount_yin,yancount_yin,yincount_yan,yancount_yan;
  _gemgrid gemin, gemout;
  _ygrid ygrid;
  _zone zones[NZONES];
}_gridset;

typedef struct {
  wordint  ip1, ip2, ip3;
  wordint date;
  wordint npas, deet, nbits;
  wordint hemisphere,axe_y_inverse;
  ftnfloat xg[16], xgref[16];
  wordint  ig[16], igref[16];
  char nomvarx[8];
  char nomvary[8];
  char typvarx[4];
  char typvary[4];
  char etiketx[16];
  char etikety[16];
} _fstinfo;

typedef struct TGeoRef {
   char*   Name;                                          ///< Reference name
   int*    Ids;                                           ///< Ids des georeferences (>=0 = ezscint)
   int     NbId,NId;                                      ///< Nombre de sous-grille

   int     NRef;                                          ///< Nombre de reference a la georeference
   int     Type;                                          ///< Type de grille
   int     BD;                                            ///< Bordure
   int     NX,NY,X0,Y0,X1,Y1;                             ///< Grid limits
   int     IG1_JP,IG2_JP,IG3_JP,IG4_JP;                   ///< Grid descriptor id
   char    Grid[3];                                       ///< Type de grille
   
   Coord  Loc;                                            ///< (Radar) Localisation du centre de reference
   double CTH,STH;                                        ///< (Radar) sin and cos of sweep angle
   double ResR,ResA;                                      ///< (Radar) Resolutions en distance et azimuth
   int    R;                                              ///< (Radar) Rayon autour du centre de reference en bin

   unsigned int NIdx,*Idx;                                ///< Index dans les positions
   float        *Lat,*Lon;                                ///< Coordonnees des points de grilles (Spherical)
   float        *AX_JP,*AY,*Hgt;                          ///< Axes de positionnement / deformation
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

   struct TGeoRef    *RefFrom;                            ///< Georeference de reference (coupe verticale,...)

   TGeoRef_Project   *Project;
   TGeoRef_UnProject *UnProject;
   TGeoRef_Value     *Value;
   TGeoRef_Distance  *Distance;
   TGeoRef_Height    *Height;

   // _Grille struct from ezscint.h
   wordint index;
   wordint flags;
   wordint i1, i2, j1, j2;
   wordint ni,nj;
   wordint ni_ax, nj_ay;
   wordint extension;
   wordint needs_expansion;
   wordint access_count;
   wordint next_gd;
   wordint n_gdin, idx_last_gdin, n_gdin_for;
   wordint log_chunk_gdin;
   wordint *gdin_for, *mask;
   wordint nsubgrids;
   struct TGeoRef *mymaskgrid;
   wordint mymaskgridi0,mymaskgridi1;
   wordint mymaskgridj0,mymaskgridj1;
   wordint *subgrid;
   ftnfloat *lat, *lon; // same as Lat,Lon
   ftnfloat *ax, *ay; // same as AX,AY
   ftnfloat *ncx, *ncy;
   char grtyp[4], grref[4];
   _fstinfo fst;
   _gridset *gset;

#ifdef HAVE_RPNC   
   int NC_Id,NC_NXDimId,NC_NYDimId;                       // netCDF identifiers
#endif
} TGeoRef;

extern TGeoRef** Grille;

typedef struct TGeoPos {
   TGeoRef *GRef;                                         // Reference horizontale
   TZRef   *ZRef;                                         // Reference verticale
   Vect3d **Pos;                                          // Coordonnees des points de grilles (World)
   int      NRef;                                         // Nombre de reference a la georeference
} TGeoPos;

typedef struct TGeoScan {
   double *X,*Y;
   float  *D;
   unsigned int *V;                                       // Coordonnees et valeurs
   unsigned int N,S;                                      // Nombre de coordonnees et dimension
   int DX,DY;                                             // Longueur em X et Y
} TGeoScan;

void GeoRef_Lock(void);
void GeoRef_Unlock(void);

TGeoRef* GeoRef_Get(char *Name);
int      GeoRef_Incr(TGeoRef *Ref);
void     GeoRef_Decr(TGeoRef *Ref);
int      GeoRef_Within(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1);
int      GeoRef_WithinRange(TGeoRef* __restrict const Ref,double Lat0,double Lon0,double Lat1,double Lon1,int In);
int      GeoRef_WithinCell(TGeoRef *GRef,Vect2d Pos,Vect2d Pt[4],int Idx0,int Idx1,int Idx2,int Idx3);
int      GeoRef_Intersect(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1,int *X0,int *Y0,int *X1,int *Y1,int BD);
int      GeoRef_Equal(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1);
int      GeoRef_CellDims(TGeoRef *Ref,int Invert,float* DX,float* DY,float* DA);
TGeoRef* GeoRef_New();
TGeoRef* GeoRef_Add(TGeoRef *Ref);
TGeoRef* GeoRef_Find(TGeoRef *Ref);
TGeoRef* GeoRef_Copy(TGeoRef* __restrict const Ref);
TGeoRef *GeoRef_HardCopy(TGeoRef* __restrict const Ref);
TGeoRef* GeoRef_Reference(TGeoRef* __restrict const Ref);
void     GeoRef_Size(TGeoRef *Ref,int X0,int Y0,int X1,int Y1,int BD);
TGeoRef* GeoRef_Resize(TGeoRef* __restrict const Ref,int NI,int NJ);
int      GeoRef_Free(TGeoRef *Ref);
void     GeoRef_Clear(TGeoRef *Ref,int New);
void     GeoRef_Qualify(TGeoRef* __restrict const Ref);
int      GeoRef_Limits(TGeoRef* __restrict const Ref,double *Lat0,double *Lon0,double *Lat1,double *Lon1);
int      GeoRef_BoundingBox(TGeoRef* __restrict const Ref,double Lat0,double Lon0,double Lat1,double Lon1,double *I0,double *J0,double *I1,double *J1);
int      GeoRef_Valid(TGeoRef* __restrict const Ref);

TGeoRef* GeoRef_RDRCreate(double Lat,double Lon,double Height,int R,double ResR,double ResA);
TGeoRef* GeoRef_RPNCreate(int NI,int NJ,char *GRTYP,int ig1,int ig2,int ig3,int ig4,int FID);
TGeoRef* GeoRef_RPNGridZE(TGeoRef *GRef,int NI,int NJ,float DX,float DY,float LatR,float LonR,int MaxCFL,float XLat1,float XLon1,float XLat2,float XLon2);
TGeoRef* GeoRef_WKTCreate(int ni,int nj,char *grtyp,int ig1,int ig2,int ig3,int ig4,char *String,double *Transform,double *InvTransform,OGRSpatialReferenceH Spatial);
int      GeoRef_WKTSet(TGeoRef *Ref,char *String,double *Transform,double *InvTransform,OGRSpatialReferenceH Spatial);
TGeoRef* GeoRef_RDRCheck(double Lat,double Lon,double Height,double Radius,double ResR,double ResA);
void     GeoRef_Expand(TGeoRef *Ref);
int      GeoRef_Positional(TGeoRef *Ref,struct TDef *XDef,struct TDef *YDef);
int      GeoRef_Coords(TGeoRef *Ref,float *Lat,float *Lon);
TQTree*  GeoRef_BuildIndex(TGeoRef* __restrict const Ref);
int      GeoRef_Nearest(TGeoRef* __restrict const Ref,double X,double Y,int *Idxs,double *Dists,int NbNear);

void GeoScan_Init(TGeoScan *Scan);
void GeoScan_Clear(TGeoScan *Scan);
int  GeoScan_Get(TGeoScan *Scan,TGeoRef *ToRef,struct TDef *ToDef,TGeoRef *FromRef,struct TDef *FromDef,int X0,int Y0,int X1,int Y1,int Dim,char *Degree);

double GeoFunc_RadialPointRatio(Coord C1,Coord C2,Coord C3);
int    GeoFunc_RadialPointOn(Coord C1,Coord C2,Coord C3,Coord *CR);
int    GeoFunc_RadialIntersect(Coord C1,Coord C2,double CRS13,double CRS23,Coord *C3);

static inline double GeoRef_GeoDir(TGeoRef* __restrict const Ref,double X, double Y) {
   
   double latd[2],lond[2],dir=0.0;
   
   if (Ref->Grid[0]!='Y' && Ref->Grid[0]!='M') {
      
      // Reproject vector orientation by adding grid projection's north difference
      Ref->Project(Ref,X,Y,&latd[0],&lond[0],1,1);
      Ref->Project(Ref,X,Y+1,&latd[1],&lond[1],1,1);

      latd[0]=DEG2RAD(latd[0]); lond[0]=DEG2RAD(lond[0]);
      latd[1]=DEG2RAD(latd[1]); lond[1]=DEG2RAD(lond[1]);
      dir=COURSE(latd[0],lond[0],latd[1],lond[1]);
   }
   return(dir);
}

#endif
