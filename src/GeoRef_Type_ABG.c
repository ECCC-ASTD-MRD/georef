#include <App.h>
#include "GeoRef.h"


//! \file


//! Transform XY grid coordinates to LatLon for an A grid
int32_t GeoRef_LL2XY_A(
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
    const float dellon = 360.0 / Ref->NX;
    const float xlon0  = 0.0;
    float dellat, xlat0;
    switch(Ref->RPNHead.ig1) {
        case GRID_GLOBAL: dellat = 180.0 / Ref->NY; xlat0 = -90.0 + dellat * 0.5; break;
        case GRID_NORTH:  dellat = 90.0 / Ref->NY;  xlat0 =  dellat * 0.5;        break;
        case GRID_SOUTH:  dellat = 90.0 / Ref->NY;  xlat0 = -90.0 + dellat * 0.5; break;
    }

    GeoRef_LL2GD(Ref, X, Y, Lat, Lon, Nb, xlat0, xlon0, dellat, dellon);

    return 0;
}


//! Transform XY grid coordinates to LatLon for a B grid
int32_t GeoRef_LL2XY_B(
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
    const float dellon = 360.0 / (Ref->NX - 1);
    const float xlon0 = 0.0;
    float dellat, xlat0;
    switch(Ref->RPNHead.ig1) {
        case GRID_GLOBAL: dellat = 180.0 / (Ref->NY-1); xlat0 = -90.0; break;
        case GRID_NORTH:  dellat = 90.0 / (Ref->NY-1);  xlat0 = 0.0;   break;
        case GRID_SOUTH:  dellat = 90.0 / (Ref->NY-1);  xlat0 = -90.0; break;
    }

    GeoRef_LL2GD(Ref, X, Y, Lat, Lon, Nb, xlat0, xlon0, dellat, dellon);

    return 0;
}


//! Transforms XY grid coordinates to LatLon for an G grid (Gaussian)
int32_t GeoRef_LL2XY_G(
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
    const float dellon = 360.0 / Ref->NX;
    const float xlon0 = 0.0;
    switch(Ref->RPNHead.ig1) {
        case GRID_GLOBAL:
            for(int32_t i = 0; i < Nb; i++) {
                X[i] = (CLAMPLON(Lon[i]) - xlon0)/dellon;
                int32_t indy = GeoRef_XFind(Lat[i], Ref->AX, Ref->NX, 1);
                if (indy>Ref->NY) indy = Ref->NY - 2;

                Y[i] = indy + (Lat[i] - Ref->AX[indy]) / (Ref->AX[indy + 1] - Ref->AX[indy]);
            }
            break;

        case GRID_NORTH: {
            const float dellat = 90.0 / Ref->NY;
            const float xlat0 =  dellat * 0.5;
            GeoRef_LL2GD(Ref, X, Y, Lat, Lon, Nb, xlat0, xlon0, dellat, dellon);
            break;
        }

        case GRID_SOUTH: {
            const float dellat = 90.0 / Ref->NY;
            const float xlat0 = -90.0 + dellat * 0.5;
            GeoRef_LL2GD(Ref, X, Y, Lat, Lon, Nb, xlat0, xlon0, dellat, dellon);
            break;
        }
    }

    return 0;
}
