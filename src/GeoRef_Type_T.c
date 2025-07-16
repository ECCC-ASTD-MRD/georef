#include <App.h>
#include "GeoRef.h"


//! \file


//! Transform XY grid coordinates to LatLon for a T grid
int32_t GeoRef_XY2LL_T(
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
    f77name(ez8_vtllfxy)(Lat, Lon, X, Y, &Ref->RPNHeadExt.xg3, &Ref->RPNHeadExt.xg4, &Ref->RPNHeadExt.xg1, &Ref->RPNHeadExt.xg2, &Ref->NX, &Ref->NY, &Nb);
    return 0;
}


//! Transforms LatLon coordinates to XY for a T grid
int32_t GeoRef_LL2XY_T(
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
    float  clat, clon, dgrw, d60;
    f77name(cigaxg)(Ref->GRTYP, &d60, &dgrw, &clat, &clon, &Ref->RPNHead.ig1, &Ref->RPNHead.ig2, &Ref->RPNHead.ig3, &Ref->RPNHead.ig4, 1);
    f77name(ez8_vtxyfll)(X, Y, Lat, Lon, &clat, &clon, &d60, &dgrw, &Ref->NX, &Ref->NY, &Nb);
    return 0;
}
