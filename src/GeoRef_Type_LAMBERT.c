#include <App.h>
#include "GeoRef.h"


//! \file


//! Transform XY grid coordinates to LatLon for a LAMBERT grid
int32_t GeoRef_XY2LL_LAMBERT(
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
    f77name(ez8_llflamb)(Lat, Lon, X, Y, &Nb, Ref->GRTYP, &Ref->RPNHead.ig1, &Ref->RPNHead.ig2, &Ref->RPNHead.ig3, &Ref->RPNHead.ig4);
    return 0;
}


//! Transform LatLon coordinates to XY for a LAMBERT grid
int32_t GeoRef_LL2XY_LAMBERT(
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
    f77name(ez8_lambfll)(X, Y, Lat, Lon, &Nb, Ref->GRTYP, &Ref->RPNHead.ig1, &Ref->RPNHead.ig2, &Ref->RPNHead.ig3, &Ref->RPNHead.ig4);
    return 0;
}
