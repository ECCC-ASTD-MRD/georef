#include <App.h>
#include "GeoRef.h"


//! \file


//! Transform XY grid coordinates to LatLon for an E grid
int32_t GeoRef_LL2XY_E(
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
    float xlat1, xlon1, xlat2, xlon2;
    f77name(cigaxg)(Ref->GRTYP, &xlat1, &xlon1, &xlat2, &xlon2, &Ref->RPNHead.ig1, &Ref->RPNHead.ig2, &Ref->RPNHead.ig3, &Ref->RPNHead.ig4, 1);
    GeoRef_RotateXY(Lat, Lon, X, Y, Nb, xlat1, xlon1, xlat2, xlon2);

    const float dellon = 360.0 / (Ref->NX - 1);
    const float xlon0 = 0.0;

    const float dellat = 180.0 / (Ref->NY);
    const float xlat0 = -90. + 0.5 * dellat;

    GeoRef_LL2GD(Ref, X, Y, Lat, Lon, Nb, xlat0, xlon0, dellat, dellon);
    return 0;
}


//! Transform LatLon coordinates to XY for an E grid
int32_t GeoRef_XY2LL_E(
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
    for (int32_t i = 0; i < Nb; i++) {
        const double dlat  = 180.0 / Ref->NY;
        const double dlon  = 360.0 / (Ref->NX - 1);
        const double swlat = -90.0 + 0.5 * dlat;
        const double swlon = 0.0;
        Lon[i] = X[i] * dlon + swlon;
        Lat[i] = Y[i] * dlat + swlat;
    }

    GeoRef_RotateInvertXY(Lat, Lon, Lon, Lat, Nb, Ref->RPNHeadExt.xg1, Ref->RPNHeadExt.xg2, Ref->RPNHeadExt.xg3, Ref->RPNHeadExt.xg4);
    return 0;
}
