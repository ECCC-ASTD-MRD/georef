#include <App.h>
#include <str.h>
#include "GeoRef.h"
#include "georef/Def.h"
#include "georef/Vertex.h"

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_WKTValue>
 * Creation     : Mars 2005 J.P. Gauthier - CMC/CMOE
 *
 * But          : Extraire la valeur d'une matrice de donnees.
 *
 * Parametres    :
 *   <Ref>      : Pointeur sur la reference geographique
 *   <Def>       : Pointeur sur la definition de la donnee
 *   <Interp>    : Mode d'interpolation
 *   <C>         : Composante
 *   <X>         : coordonnee en X dans la projection/grille
 *   <Y>         : coordonnee en Y dans la projection/grille
 *   <Z>         : coordonnee en Z dans la projection/grille
 *   <Length>    : Module interpolee
 *   <ThetaXY>   : Direction interpolee
 *
 * Retour       : Inside (1 si a l'interieur du domaine).
 *
 * Remarques   :
 *
 *---------------------------------------------------------------------------------------------------------------
*/
// int32_t GeoRef_WKTValue(TGeoRef *Ref, TDef *Def, TRef_Interp Interp, int32_t C, double X, double Y, double Z, double *Length, double *ThetaXY){

//    double       x, y, d, ddir=0.0;
//    int32_t          valid=0, mem, ix, iy;
//    uint32_t idx;

//   *Length=Def->NoData;
//    d=Ref->Type[[0]=='W'?1.0:0.5;

//    //i on est a l'interieur de la grille ou que l'extrapolation est activee
//    if (C<Def->NC && X>=(Ref->X0-d) && Y>=(Ref->Y0-d) && Z>=0 && X<=(Ref->X1+d) && Y<=(Ref->Y1+d) && Z<=Def->NK-1) {

//       X-=Ref->X0;
//       Y-=Ref->Y0;
//       DEFCLAMP(Def, X, Y);

//       // Index memoire du niveau desire
//       mem=Def->NIJ*(int)Z;
//       ix=lrint(X);
//       iy=lrint(Y);
//       idx=iy*Def->NI+ix;

//       // Check for mask
//       if (Def->Mask && !Def->Mask[idx]) {
//          return(valid);
//       }

//       // Reproject vector orientation by adding grid projection's north difference
//       if (Def->Data[1] && Ref->Type&GRID_ROTATED) {
//          ddir=GeoRef_GeoDir(Ref, X, Y);
//       }

//       if (Def->Type<=9 || Interp==IR_NEAREST || (X==ix && Y==iy)) {
//          mem+=idx;
//          Def_GetMod(Def, mem, *Length);

//          // Pour un champs vectoriel
//          if (Def->Data[1] && ThetaXY) {
//             Def_Get(Def, 0, mem, x);
//             Def_Get(Def, 1, mem, y);
//             *ThetaXY=180+RAD2DEG(atan2(x, y)-ddir);
//          }
//       } else {
//          *Length=VertexVal(Def, -1, X, Y, Z);
//          // Pour un champs vectoriel
//          if (Def->Data[1] && ThetaXY) {
//             x=VertexVal(Def, 0, X, Y, Z);
//             y=VertexVal(Def, 1, X, Y, Z);
//             *ThetaXY=180+RAD2DEG(atan2(x, y)-ddir);
//          }
//       }
//       valid=DEFVALID(Def, *Length);
//    }

//    return(valid);
// }


//! Change a georeference type to W
TGeoRef * GeoRef_DefineW(
    //! [inout] Georeferance
    TGeoRef * const Ref,
    //! [in] WKT description
    const char * const String,
    //! [in] Transformation matrix [2][6]
    const double * const Transform,
    //! [in] Inverse transformation matrix [2][6]
    const double * const InvTransform,
    //! [in] Spatial reference (GDAL object)
    const OGRSpatialReferenceH Spatial
) {
    //! \return New georeference or NULL on error
#ifdef HAVE_GDAL
    static OGRSpatialReferenceH llref = NULL;
    char * string = NULL;

    if (String && String[0] != '\0') {
        string = strdup(String);
        strtrim(string, ' ');
        Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s: Projection string: %s\n", __func__, string);
    }

    if (Transform || InvTransform) {
        if (!Ref->Transform) Ref->Transform = (double*)calloc(6, sizeof(double));
        if (!Ref->InvTransform) Ref->InvTransform = (double*)calloc(6, sizeof(double));
    }

    if (Transform) {
        memcpy(Ref->Transform, Transform, 6*sizeof(double));
    } else {
        if (!InvTransform || !GDALInvGeoTransform(InvTransform, Ref->Transform)) {
            if (Ref->Transform) {
                free(Ref->Transform);
                Ref->Transform=NULL;
            }
        }
    }

    if (InvTransform) {
        memcpy(Ref->InvTransform, InvTransform, 6*sizeof(double));
    } else {
        if (!Transform || !GDALInvGeoTransform(Transform, Ref->InvTransform)) {
            if (Ref->InvTransform) {
                free(Ref->InvTransform);
                Ref->InvTransform=NULL;
            }
        }
    }

    if (Spatial) {
        Ref->Spatial = OSRClone(Spatial);
        OSRExportToWkt(Ref->Spatial, &string);
        Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s: Projection from spatial:%p\n", __func__, string);
    } else if (string) {
        Ref->Spatial = OSRNewSpatialReference(NULL);
        if (OSRSetFromUserInput(Ref->Spatial, string) == OGRERR_FAILURE) {
            Lib_Log(APP_LIBGEOREF, APP_WARNING, "%s: Unable to create spatial reference\n", __func__);
            return NULL;
        }
    } else {
        string = strdup(REF_DEFAULT);
        Ref->Spatial =OSRNewSpatialReference(string);
    }

    if (Ref->String) free(Ref->String);
    Ref->String = string;

    if (Ref->Spatial) {
        if (!llref) {
            // Create global latlon reference on perfect sphere
            llref = OSRNewSpatialReference(NULL);
            OSRSetFromUserInput(llref, "EPSG:4047");
#if defined(GDAL_VERSION_MAJOR) && GDAL_VERSION_MAJOR >= 3
            OSRSetAxisMappingStrategy(llref, OAMS_TRADITIONAL_GIS_ORDER);
#endif
        }

        if (llref) {
            // Create forward/backward tranformation functions
            Ref->Function = OCTNewCoordinateTransformation(Ref->Spatial, llref);
            Ref->InvFunction = OCTNewCoordinateTransformation(llref, Ref->Spatial);
            if (!Ref->Function || !Ref->InvFunction) {
                Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable to create transformation functions\n", __func__);
            }
        }
    } else {
        Lib_Log(APP_LIBGEOREF, APP_WARNING, "%s: Unable to get spatial reference\n", __func__);
        return NULL;
    }

    if (Ref->GRTYP[0] == ' ' || Ref->GRTYP[0] == 'X' || Ref->GRTYP[0] == '\0') {
        Ref->GRTYP[0] = 'W';
    }

    // We might be defining a Z frid mapped on W so we don't want to rerun the qualification
    if (Ref->GRTYP[0] == 'W') {
        Ref->Height = NULL;
        Ref->Type = GRID_NONE;
        GeoRef_Qualify(Ref);
    }
    return Ref;
#else
    Lib_Log(APP_LIBGEOREF, APP_ERROR, "Function %s is not available, needs to be built with GDAL\n", __func__);
    return NULL;
#endif
}


//! Create a W type georeference
TGeoRef * GeoRef_CreateW(
    //! [in] Number of gridpoints in X
    const int32_t NI,
    //! [in] Number of gridpoints in Y
    const int32_t NJ,
    //! [in] WKT description
    const char * const String,
    //! [in] Transformation matrix [2][6]
    const double * const Transform,
    //! [in] Inverse transformation matrix [2][6]
    const double * const InvTransform,
    //! [in] Spatial reference (GDAL object)
    const OGRSpatialReferenceH Spatial
) {
    //! \return New georeference or NULL on error

    TGeoRef * const ref = GeoRef_New();
    GeoRef_Size(ref, 0, 0, NI > 0 ? NI -1 : 0, NJ > 0 ? NJ - 1 : 0, 0);
    if (!GeoRef_DefineW(ref, String, Transform, InvTransform, Spatial)) {
        return NULL;
    }

    ref->GRTYP[0] = 'W';
    ref->GRTYP[1] = '\0';
    TGeoRef * const fref = GeoRef_Find(ref);
    if (fref) {
        // This georef already exists
        free(ref);
        GeoRef_Incr(fref);
        return fref;
    }

    // This is a new georef
    GeoRef_Add(ref);
    GeoRef_Qualify(ref);

    return ref;
}

#ifdef HAVE_GDAL
//! Apply rotation matrix to a LatLon pair
static inline int32_t GeoRef_WKTRotate(
    //! [in] Transformation matrix
    const TRotationTransform * const T,
    //! [inout] Latitude
    double * const Lat,
    //! [inout] Longitude
    double * const Lon
) {
    //! \returns Always true

    const double lon = DEG2RAD(*Lon);
    const double lat = DEG2RAD(*Lat);

    // Convert from spherical to cartesian coordinates
    const double x = cos(lon) * cos(lat);
    const double y = sin(lon) * cos(lat);
    const double z = sin(lat);

    const double xr = T->CosTheta * T->CosPhi * x + T->CosTheta * T->SinPhi * y + T->SinTheta * z;
    const double yr = -T->SinPhi * x + T->CosPhi * y;
    const double zr = -T->SinTheta * T->CosPhi * x - T->SinTheta * T->SinPhi * y + T->CosTheta * z;

    // Convert cartesian back to spherical coordinates
    *Lon = RAD2DEG(atan2(yr, xr));
    *Lat = RAD2DEG(asin(zr));

    return TRUE;
}


//! Apply inverse rotation matrix to a LatLon pair
static inline int32_t GeoRef_WKTUnRotate(
    //! [in] Transformation matrix
    const TRotationTransform * const T,
    //! [inout] Latitude
    double * const Lat,
    //! [inout] Longitude
    double * const Lon
) {
    //! \returns Always true

    const double lon = DEG2RAD(*Lon);
    const double lat = DEG2RAD(*Lat);

    // Convert from spherical to cartesian coordinates
    const double x = cos(lon) * cos(lat);
    const double y = sin(lon) * cos(lat);
    const double z = sin(lat);

    const double xr = T->CosTheta * T->CosPhi * x + -T->SinPhi * y + -T->SinTheta * T->CosPhi * z;
    const double yr = -T->CosTheta * -T->SinPhi * x + T->CosPhi * y - -T->SinTheta * -T->SinPhi * z;
    const double zr = T->SinTheta * x + T->CosTheta * z;

    // Convert cartesian back to spherical coordinates
    *Lon = RAD2DEG(atan2(yr, xr));
    *Lat = RAD2DEG(asin(zr));

    return TRUE;
}
#endif


//! Transform XY grid coordinates to LatLon for a W grid
int32_t GeoRef_XY2LL_W(
    //! [in] Georeference
    const TGeoRef * const Ref,
    //! [out] Latitude array
    double * const Lat,
    //! [out] Longitude array
    double * const Lon,
    //! [in] X array
    const double * const X,
    //! [in] Y array
    const double * const Y,
    //! [in] Number of coordinates
    const int32_t Nb
) {
    //! \returns Always 0
#ifdef HAVE_GDAL
    double z = 0.0, d;
    int32_t ok;

    for(int32_t n = 0; n < Nb; n++) {
        // Transform the point32_t into georeferenced coordinates
        double x = X[n];
        double y = Y[n];
        if (Ref->Options.Transform) {
            if (Ref->Transform) {
                x = Ref->Transform[0] + Ref->Transform[1] * X[n] + Ref->Transform[2] * Y[n];
                y = Ref->Transform[3] + Ref->Transform[4] * X[n] + Ref->Transform[5] * Y[n];
            } else if (Ref->GCPTransform) {
                GDALGCPTransform(Ref->GCPTransform, FALSE, 1, &x, &y, &z, &ok);
            } else if (Ref->TPSTransform) {
                GDALTPSTransform(Ref->TPSTransform, FALSE, 1, &x, &y, &z, &ok);
            } else if (Ref->RPCTransform) {
                GDALRPCTransform(Ref->RPCTransform, FALSE, 1, &x, &y, &z, &ok);
            }
        }

        // Transform to latlon
        if (Ref->Function) {
            if (!OCTTransform(Ref->Function, 1, &x, &y, NULL)) {
                Lon[n] = -999.0;
                Lat[n] = -999.0;
                continue;
            }
        }

        Lon[n] = x;
        Lat[n] = y;

        if (Ref->RotTransform) GeoRef_WKTUnRotate(Ref->RotTransform, &Lat[n], &Lon[n]);
    }
#else
    Lib_Log(APP_LIBGEOREF, APP_ERROR, "Function %s is not available, needs to be built with GDAL\n", __func__);
#endif
    return 0;
}


//! Transform LatLon coordinates to XY for a Z grid
int32_t GeoRef_LL2XY_W(
    //! [in] Georeference
    const TGeoRef * const Ref,
    //! [out] X array
    double * const X,
    //! [out] Y array
    double * const Y,
    //! [in] Latitude array
    const double * const Lat,
    //! [in] Longitude array
    const double * const Lon,
    //! [in] Number of coordinates
    const int32_t Nb
) {
    //! \return Always 0
#ifdef HAVE_GDAL
    double z = 0.0;
    int32_t ok;
    for(int32_t n = 0; n < Nb; n++) {
        if (Ref->RotTransform) GeoRef_WKTRotate(Ref->RotTransform, &Lat[n], &Lon[n]);

        // Longitude from -180 to 180
        CLAMPLON(Lon[n]);

        double x = Lon[n];
        double y = Lat[n];

        // Transform from latlon
        if (Ref->InvFunction) {
            if (!OCTTransform(Ref->InvFunction, 1, &x, &y, NULL)) {
                X[n] = -1.0;
                Y[n] = -1.0;
                continue;
            }
        }

        // Transform from georeferenced coordinates
        X[n] = x;
        Y[n] = y;
        if (Ref->Options.Transform) {
            if (Ref->InvTransform) {
                X[n] = Ref->InvTransform[0] + Ref->InvTransform[1] * x + Ref->InvTransform[2] * y;
                Y[n] = Ref->InvTransform[3] + Ref->InvTransform[4] * x + Ref->InvTransform[5] * y;
            } else if (Ref->GCPTransform) {
                GDALGCPTransform(Ref->GCPTransform, TRUE, 1, &X[n], &Y[n], &z, &ok);
            } else if (Ref->TPSTransform) {
                GDALTPSTransform(Ref->TPSTransform, TRUE, 1, &X[n], &Y[n], &z, &ok);
            } else if (Ref->RPCTransform) {
                GDALRPCTransform(Ref->RPCTransform, TRUE, 1, &X[n], &Y[n], &z, &ok);
            }
        }

    }
#else
    Lib_Log(APP_LIBGEOREF, APP_ERROR, "Function %s is not available, needs to be built with GDAL\n", __func__);
#endif
    return 0;
}
