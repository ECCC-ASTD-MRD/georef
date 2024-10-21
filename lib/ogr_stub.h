#ifndef _ogr_stub_h
#define _ogr_stub_h

typedef char  OGRGeometryH;
typedef char* OGRCoordinateTransformationH;
typedef char* OGRSpatialReferenceH;
typedef char* OGRLayerH;
typedef char* OGRFeatureH;
typedef char* OGRFeatureDefnH;
typedef char* OGRDataSourceH;
typedef char* OGRSFDriverH;
typedef char* GDALDatasetH;

typedef struct OGREnvelope {
    double      MinX;
    double      MaxX;
    double      MinY;
    double      MaxY;
} OGREnvelope;

#endif