#ifndef _ogr_stub_h
#define _ogr_stub_h

#define OGRGeometryH                 char
#define OGRCoordinateTransformationH char*
#define OGRSpatialReferenceH         char*
#define OGRLayerH                    char*
#define OGRFeatureH                  char*
#define OGRFeatureDefnH              char*
#define OGRDataSourceH               char*
#define OGRSFDriverH                 char*
#define GDALDatasetH                 char*

typedef struct OGREnvelope {
    double      MinX;
    double      MaxX;
    double      MinY;
    double      MaxY;
} OGREnvelope;

#endif