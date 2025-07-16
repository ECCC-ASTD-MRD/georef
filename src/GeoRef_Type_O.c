#include <App.h>
#include "GeoRef.h"
#include "georef/Vertex.h"


//! \file


//! Transform XY grid coordinates to LatLon for an O grid
int32_t GeoRef_XY2LL_O(
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

    #pragma omp parallel for default(none) private(i) shared(Nb, Ref, X, Y, Lat, Lon)
    for (int32_t i = 0; i < Nb; i++) {
        Lon[i] = VertexValS(Ref->AX, NULL, Ref->NX, Ref->NY, X[i], Y[i], TRUE);
        Lat[i] = VertexValS(Ref->AY, NULL, Ref->NX, Ref->NY, X[i], Y[i], TRUE);
    }
    return Nb;
}


//! Transform LatLon coordinates to XY for an O grid
int32_t GeoRef_LL2XY_O(
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
    //TODO: Check for transformation
    //   GeoRef_LL2GREF(Ref, X, Y, Lat, Lon, Nb);

    int32_t out = 0;
    #pragma omp parallel for default(none) private(n, nd, d, dx, dy, x, y, xx, yy, idx, idxs, pt, pts, dists) shared(Nb, Ref, X, Y, Lat, Lon) reduction(+:out)
    for(int32_t d = 0; d < Nb; d++) {
        X[d] = -1.0;
        Y[d] = -1.0;

        int32_t idxs[8];
        double dists[8];
        int32_t nd = GeoRef_Nearest(Ref, Lon[d], Lat[d], idxs, dists, 8, 0.0);
        if (nd) {
            Vect2d pt = {Lon[d], Lat[d]};

            // Find which cell includes coordinates
            int32_t idx;
            for(int32_t n = 0; n < nd; n++) {
                idx = idxs[n];

                // Find within which quad
                int32_t dx = -1, dy = -1;
                if (!GeoRef_WithinCell(Ref, pt, pts, idx-Ref->NX - 1, idx - 1, idx, idx-Ref->NX)) {
                    dx = 0; dy = -1;
                    if (!GeoRef_WithinCell(Ref, pt, pts, idx-Ref->NX, idx, idx + 1, idx-Ref->NX + 1)) {
                        dx = -1; dy = 0;
                        if (!GeoRef_WithinCell(Ref, pt, pts, idx - 1, idx + Ref->NX - 1, idx + Ref->NX, idx)) {
                            dx = 0; dy = 0;
                            if (!GeoRef_WithinCell(Ref, pt, pts, idx, idx + Ref->NX, idx + Ref->NX + 1, idx + 1)) {
                                idx = -1;
                            }
                        }
                    }
                }
                // If found, exit loop
                if (idx != -1) break;
            }

            if (idx != -1) {
                // Map coordinates to grid
                Vect2d pts[4];
                double xx, yy;
                Vertex_Map(pts, &xx, &yy, Lon[d], Lat[d]);
                X[d] = xx; Y[d] = yy;

                if (!ISNAN(X[d]) && !ISNAN(Y[d])) {
                    int32_t y = idx / Ref->NX;
                    int32_t x = idx - y * Ref->NX;
                    Y[d] += y + dy;
                    X[d] += x + dx;
                } else {
                    Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid coordinate (NAN): ll(%f, %f) xy(%f, %f) %i\n", __func__, Lat[d], Lon[d], X[d], Y[d], idx);
                    X[d] = -1.0;
                    Y[d] = -1.0;
                    out++;
                }
            } else {
                Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Point not found: %f %f\n", __func__, Lat[d], Lon[d]);
                out++;
            }
        }
    }
    Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s: Points out: %i\n", __func__, out);

    return out;
}
