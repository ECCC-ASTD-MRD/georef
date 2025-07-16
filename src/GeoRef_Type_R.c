#include <App.h>
#include "GeoRef.h"


//! \file


//! Calculate MAGL height of a grid coordinate
double GeoRef_RDRHeight(
    //! [in] Georeference
    const TGeoRef * const Ref,
    //! [in] Vertical reference
    const TZRef * const ZRef,
    //! [in] Azimuth (x) coordinate in grid units
    const double Azimuth,
    //! [in] Bin (y) coordinate in grid units
    const double Bin,
    //! [in] Sweep (z) coordinate in grid units
    const double Sweep
) {
    //! \return MAGL height in meters. 0 on error

    // Prevent warning about unused parameter
    (void)Azimuth;

    if (Bin >= 0 && Bin < Ref->R && Sweep >= 0 && Sweep < ZRef->LevelNb) {
        return Ref->Loc.Elev + sin(DEG2RAD(ZRef->Levels[(int)Sweep])) * ((int)Bin * Ref->ResR);
    } else {
        return 0;
    }
}


//! Create an R type georeference (Radial)
TGeoRef * GeoRef_CreateR(
    //! [in] Center latitude
    double Lat,
    //! [in] Center longitude
    double Lon,
    //! [in] Center altitude (m)
    double Height,
    //! [in] Radius (m)
    int32_t R,
    //! [in] Resolution along the radius
    double ResR,
    //! [in] Resolution in azimuth
    double ResA
) {
    //! \returns Georeference
    TGeoRef * const ref = GeoRef_New();

    ref->GRTYP[0] = 'R';
    ref->GRTYP[1] = '\0';
    ref->Loc.Lat = Lat;
    ref->Loc.Lon = Lon;
    ref->Loc.Elev = Height;
    ref->R = R;
    ref->ResR = ResR;
    ref->ResA = ResA;
    ref->Height = GeoRef_RDRHeight;

    GeoRef_Size(ref, 0, 0, 360 / ResA, R - 1, 0);

    TGeoRef * const fref = GeoRef_Find(ref)
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


//! Transform XY grid coordinates to LatLon for a R grid (Radial)
int32_t GeoRef_XY2LL_R(
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
    const TCoord loc0 = {
        .Lat = DEG2RAD(Ref->Loc.Lat),
        .Lon = DEG2RAD(Ref->Loc.Lon)
    };

    #pragma omp parallel for default(none) private(n, x, y, d) shared(Nb, Ref, X, Y, Lat, Lon, loc0)
    for(int32_t n = 0; n < Nb; n++) {
        double x = X[n] * Ref->ResA;
        double y = Y[n] * Ref->ResR;

        x = DEG2RAD(x);
        const double d = M2RAD(y * Ref->CTH);

        if (Ref->Options.Transform) {
            Lat[n] = asin(sin(loc0.Lat) * cos(d) + cos(loc0.Lat) * sin(d) * cos(x));
            Lon[n] = fmod(loc0.Lon + (atan2(sin(x) * sin(d) * cos(loc0.Lat), cos(d) - sin(loc0.Lat) * sin(*Lat))) + M_PI, M_2PI) - M_PI;
            Lat[n] = RAD2DEG(Lat[n]);
            Lon[n] = RAD2DEG(Lon[n]);
        } else {
            Lat[n] = d;
            Lon[n] = x;
        }
    }
    return 0;
}


//! Transform LatLon coordinates to XY for a R grid (Radial)
int32_t GeoRef_LL2XY_R(
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
    const TCoord loc0 = {
        .Lat = DEG2RAD(Ref->Loc.Lat),
        .Lon = DEG2RAD(Ref->Loc.Lon)
    };

    #pragma omp parallel for default(none) private(n, x, d, lat, lon) shared(Nb, Ref, X, Y, Lat, Lon, loc0)
    for(int32_t n = 0; n < Nb; n++) {
        const double lat = DEG2RAD(Lat[n]);
        const double lon = DEG2RAD(Lon[n]);

        const double d = fabs(DIST(0.0, loc0.Lat, loc0.Lon, lat, lon));
        const double x = -RAD2DEG(COURSE(loc0.Lat, loc0.Lon, lat, lon));
        X[n] = x < 0.0 ? x + 360.0 : x;
        Y[n] = d / Ref->CTH;

        if (Ref->Options.Transform) {
            X[n] = X[n] / Ref->ResA;
            Y[n] = Y[n] / Ref->ResR;
        }
    }
    return 0;
}
