#ifndef _OGR_h
#define _OGR_h

#include "GeoRef.h"

#ifdef _TK_SOURCE
#include <tclDataSpec.h>
#endif

#define OGR_G_EnvelopeIntersect(ENV0,ENV1) (!(ENV0.MaxX<ENV1.MinX || ENV0.MinX>ENV1.MaxX || ENV0.MaxY<ENV1.MinY || ENV0.MinY>ENV1.MaxY))
#define OGR_PointInside(V,V0,V1)           (((V0[1]<=V[1] && V[1]<V1[1]) || (V1[1]<=V[1] && V[1]<V0[1])) && (V[0]<((V1[0]-V0[0])*(V[1]-V0[1])/(V1[1]-V0[1])+V0[0])))

#define OGR_GEOMTYPES "Point,3D Point,Line String,3D Line String,Polygon,3D Polygon,Multi Point,3D Multi Point,Multi Line String,3D Multi Line String,Multi Polygon,3D Multi Polygon,Geometry Collection,3D Geometry Collection,Linear Ring"
#define OGR_BUFFER    4096

typedef struct OGR_Sort {
   int32_t   Field,Type,Order;       // Sorting parameters
   uint32_t  Nb,*Table;              // Sorted table of featurre index
} OGR_Sort;

typedef struct OGR_File {
   GDALDatasetH    Data;                 // OGR internal file datasource object
   OGRSFDriverH    Driver;               // OGR driver used for this file
   char           *Id;                   // File identifier
   char           *Name;                 // File path
   char            Mode;                 // File access mode
} OGR_File;

typedef struct OGR_Layer {
#ifdef _TK_SOURCE
   Tcl_Obj         *Tag;                 // (Tcl_Obj type)   Tcl identifier
   TDataSpec       *Spec;                // (TDataSpec type) Specification de rendue des donnees
#else
   char            *Tag;                 // (Tcl_Obj type)   Tcl identifier
   char            *Spec;                // (TDataSpec type) Specification de rendue des donnees
#endif

   TGeoRef         *GRef;                // GeoReference
   OGRLayerH        Layer;               // OGR internal layer object
   OGRFeatureH     *Feature;             // List of OGR internal layer featuret
   OGRFeatureDefnH  Def;                 // OGR internal feature definition object
   GDALDatasetH     Data;                // OGR internal layer datasource object
   OGR_File        *File;                // Layer's file provenance
   OGR_Sort         Sort;                // Sorting parameters
   char            *Select;              // List of features selection flag
   TCoord          *Loc;                 // List of feature's centroid
   Vect3d           Vr[2];               // Layer extent in projected coordinates
   double           Min,Max;             // Layer's min-max of the currently mapped field
   int32_t          Update;              // Do we need to update the internale OGR
   int32_t          Mask,FMask;          // Is this layer used as a mask
   uint32_t         NFeature;            // Number of features in the layer
   uint32_t         LFeature;            // List of feature's rendered display list
   uint32_t         GFeature;            // Number of feature processes as display list
   uint32_t        *SFeature;            // List of of highlighted features
   uint32_t         NSFeature;           // Number of highligted features
   int32_t          CFeature;            // Cleared feature (to be re-rendered)
   int32_t          Topo,Extrude,Space;  // Positional parameters
   char             Changed;             // Is the layer changed
} OGR_Layer;

#define OGM_ARRAY0   0
#define OGM_ARRAY1   1
#define OGM_ARRAYPTR 2

#ifdef HAVE_GPC
   #include "gpc.h"
   #include "gpc_ext.h"
#else
   typedef char gpc_polygon;
   typedef char gpc_op;
#endif

void         OGM_GPCFromOGR(gpc_polygon* Poly,OGRGeometryH *Geom);
void         OGM_GPCToOGR(gpc_polygon *Poly,OGRGeometryH *Geom);
OGRGeometryH OGM_GPCOnOGR(gpc_op Op,OGRGeometryH Geom0,OGRGeometryH Geom1);
OGRGeometryH OGM_GPCOnOGRLayer(gpc_op Op,OGR_Layer *Layer);
OGRGeometryH OGM_GPCOnOGRGeometry(gpc_op Op,OGRGeometryH *Geom);
void         OGM_GPCNew(gpc_polygon *Poly);

Vect3d*      OGM_GetVect3d(uint32_t Size,uint32_t No);
void         OGM_ClearVect3d(void);
void         OGM_OGRProject(OGRGeometryH Geom,TGeoRef *FromRef,TGeoRef *ToRef);
int32_t      OGM_QSortInter(const void *A,const void *B);
int32_t      OGM_Within(OGRGeometryH Geom0,OGRGeometryH Geom1,OGREnvelope *Env0,OGREnvelope *Env1);
OGRGeometryH OGM_SegIntersectionPts(OGRGeometryH Geom,double X0,double Y0,double X1,double Y1);
int32_t      OGM_Intersect(OGRGeometryH Geom0,OGRGeometryH Geom1,OGREnvelope *Env0,OGREnvelope *Env1);
int32_t      OGM_PointPointIntersect(OGRGeometryH Geom0,OGRGeometryH Geom1,int32_t All);
int32_t      OGM_PointLineIntersect(OGRGeometryH Geom0,OGRGeometryH Geom1,int32_t All);
int32_t      OGM_PointPolyIntersect(OGRGeometryH Geom0,OGRGeometryH Geom1,int32_t All);
int32_t      OGM_PolyPolyIntersect(OGRGeometryH Geom0,OGRGeometryH Geom1);
int32_t      OGM_LinePolyIntersect(OGRGeometryH Geom0,OGRGeometryH Geom1);
int32_t      OGM_SegmentIntersect(Vect3d PointA,Vect3d PointB,Vect3d PointC,Vect3d PointD,Vect3d Inter);
double       OGM_Length(OGRGeometryH Geom);
double       OGM_SegmentLength(OGRGeometryH Geom);
double       OGM_SegmentDist(Vect3d SegA,Vect3d SegB,Vect3d Point);
double       OGM_PointClosest(OGRGeometryH Geom,OGRGeometryH Pick,Vect3d Vr);
int32_t      OGM_PointInside(OGRGeometryH Geom,OGRGeometryH Pick,Vect3d Vr);
double       OGM_CoordLimit(OGRGeometryH Geom,int32_t Coord,int32_t Mode);
OGRGeometryH OGM_Clip(OGRGeometryH Line,OGRGeometryH Poly);
int32_t      OGM_ClipSegment(OGRGeometryH Line,OGRGeometryH Poly,OGRGeometryH Clip);
double       OGM_Centroid2D(OGRGeometryH Geom,double *X,double *Y);
double       OGM_Centroid2DProcess(OGRGeometryH Geom,double *X,double *Y);
int32_t      OGM_Simplify(double Tolerance,OGRGeometryH Geom);
int32_t      OGM_SimplifyDP(double Tolerance,Vect3d *Pt,int32_t J,int32_t K,int32_t *Markers);
double       OGM_AngleMin(OGRGeometryH Geom);
int32_t      OGM_Clean(OGRGeometryH Geom);
OGRGeometryH OGM_PolySplitTile(OGRGeometryH Poly,const uint32_t MaxPoints,OGRGeometryH Res);
OGRGeometryH OGM_ClipLonWrap(OGRGeometryH Poly);

#endif
