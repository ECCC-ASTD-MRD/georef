#include <App.h>
#include "GeoRef.h"


//! \file


//! Transform XY grid coordinates to LatLon for an N or S grid (Polar stereographic)
int32_t GeoRef_XY2LL_NS(
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
    int32_t un = 1;
    f77name(ez8_vllfxy)(Lat, Lon, X, Y, &Nb, &un,
        &Ref->RPNHeadExt.xg3, &Ref->RPNHeadExt.xg4, &Ref->RPNHeadExt.xg1, &Ref->RPNHeadExt.xg2, &Ref->Hemi);
    return 0;
}


//! Transform LatLon coordinates to XY for an N or S grid (Polar stereographic)
int32_t GeoRef_LL2XY_NS(
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
    float  pi, pj, dgrw, d60;
    f77name(cigaxg)(Ref->GRTYP, &pi, &pj, &d60, &dgrw, &Ref->RPNHead.ig1, &Ref->RPNHead.ig2, &Ref->RPNHead.ig3, &Ref->RPNHead.ig4, 1);
    f77name(ez8_vxyfll)(X, Y, Lat, Lon, &Nb, &d60, &dgrw, &pi, &pj, &Ref->Hemi);
    return 0;
}

