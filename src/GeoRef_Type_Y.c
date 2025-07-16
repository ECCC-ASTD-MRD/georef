#include <App.h>
#include "GeoRef.h"


//! \file


//! Transform XY grid coordinates to LatLon for a Y grid
int32_t GeoRef_XY2LL_Y(
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
    if (!Ref->AX || !Ref->AY) {
        return FALSE;
    }

    switch (Ref->RPNHeadExt.grref[0]) {
        case 'L':
            for (int32_t i = 0; i < Nb; i++) {
                indx = ROUND(Y[i]) * (Ref->X1 - Ref->X0) + ROUND(X[i]);
                Lat[i] = Ref->AY[indx];
                Lon[i] = Ref->AX[indx];
            }
            break;

        case 'W': {
            double * const tmpx = (double*)malloc(2 * Nb * sizeof(double));
            double * const tmpy = &tmpx[Nb];

            for (int32_t i = 0; i < Nb; i++) {
                int32_t sx = floor(X[i]); sx = CLAMP(sx, Ref->X0, Ref->X1);
                int32_t sy = floor(Y[i]); sy = CLAMP(sy, Ref->Y0, Ref->Y1);
                const double dx = X[i] - sx;
                const double dy = Y[i] - sy;

                int32_t s = sy * Ref->NX + sx;
                tmpx[i] = Ref->AX[s];
                tmpy[i] = Ref->AY[s];

                if (++sx <= Ref->X1) {
                    s = sy * Ref->NX + sx;
                    tmpx[i] + =(Ref->AX[s] - tmpx[i]) * dx;
                }
                tmpx[i];

                if (++sy <= Ref->Y1) {
                    s = sy * Ref->NX + (sx - 1);
                    tmpy[i] += (Ref->AY[s] - tmpy[i]) * dy;
                }
                tmpy[i];
            }
            GeoRef_XY2LL_W(Ref, Lat, Lon, tmpx, tmpy, Nb);
            free(tmpx);
            break;
        }

        default:
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Undefined reference grid type: %s\n", __func__, Ref->RPNHeadExt.grref[0]);
            return(-1);
            break;
    }
    return Nb;
}


//! Transform LatLon coordinates to XY for a Y grid
int32_t GeoRef_LL2XY_Y(
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
    int32_t nb = 0;
    for (int32_t n = 0; n < Nb; n++) {
        int32_t idx;
        double dists[8];
        if (GeoRef_Nearest(Ref, Lon[n], Lat[n], &idx, dists, 1, Ref->Options.DistTreshold)) {
            if (dists[0] < 1.0) {
                Y[n] = (int)(idx / Ref->NX);
                X[n] = idx-Y[n] * Ref->NX;
                nb++;
            }
        }
    }
    return nb;
}
