#ifndef _ogr_stub_h
#define _ogr_stub_h

typdef char  OGRGeometryH;
typdef char* OGRCoordinateTransformationH;
typdef char* OGRSpatialReferenceH;
typdef char* OGRLayerH;
typdef char* OGRFeatureH;
typdef char* OGRFeatureDefnH;
typdef char* OGRDataSourceH;
typdef char* OGRSFDriverH;
typdef char* GDALDatasetH;

typedef struct OGREnvelope {
    double      MinX;
    double      MaxX;
    double      MinY;
    double      MaxY;
} OGREnvelope;

#endif