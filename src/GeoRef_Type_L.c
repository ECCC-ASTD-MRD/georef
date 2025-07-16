#include <App.h>
#include "GeoRef.h"


//! \file


//! Transform XY grid coordinates to LatLon for an L grid (Cylindrical)
int32_t GeoRef_XY2LL_L(
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
    int32_t i;

    #pragma omp parallel for default(none) private(i) shared(Nb, Ref, X, Y, Lat, Lon)
    for (i = 0; i < Nb; i++) {
        Lat[i] = Y[i] * Ref->RPNHeadExt.xg3 + Ref->RPNHeadExt.xg1;
        Lon[i] = X[i] * Ref->RPNHeadExt.xg4 + Ref->RPNHeadExt.xg2;
    }
    return 0;
}


//! Transform LatLon coordinates to XY for a L grid (Cylindrical)
int32_t GeoRef_LL2XY_L(
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
    float dellat, dellon, xlat0, xlon0;
    f77name(cigaxg)(Ref->GRTYP, &xlat0, &xlon0, &dellat, &dellon, &Ref->RPNHead.ig1, &Ref->RPNHead.ig2, &Ref->RPNHead.ig3, &Ref->RPNHead.ig4, 1);
    GeoRef_LL2GD(Ref, X, Y, Lat, Lon, Nb, xlat0, xlon0, dellat, dellon);
    return 0;
}
