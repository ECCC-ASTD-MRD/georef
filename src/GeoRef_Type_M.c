#include <App.h>
#include "GeoRef.h"
#include "georef/Vertex.h"
#include "georef/Triangle.h"


//! \file


//! Transform XY grid coordinates to LatLon for an M grid (Mesh or TIN)
int32_t GeoRef_LL2XY_M(
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
    TQTree *node;
    Vect3d  b;
    int32_t n, d, idx;
    double lon;

    #pragma omp parallel for default(none) private(d, node, b, n, idx, lon) shared(Nb, Ref, X, Y, Lat, Lon)
    for(d = 0; d < Nb; d++) {
        lon = CLAMPLON(Lon[d]);

        if (Ref->QTree) {
            // If there's an index use it
            if ((node = QTree_Find(Ref->QTree, lon, Lat[d])) && node->NbData) {
                // Loop on this nodes data payload
                for(n = 0; n < node->NbData; n++) {
                    idx = (intptr_t)node->Data[n].Ptr - 1; // Remove false pointer increment

                    if (Bary_Get(b, Ref->Wght ? Ref->Wght[idx / 3] : 0.0, lon, Lat[d], Ref->Lon[Ref->Idx[idx]], Ref->Lat[Ref->Idx[idx]],
                        Ref->Lon[Ref->Idx[idx + 1]], Ref->Lat[Ref->Idx[idx + 1]], Ref->Lon[Ref->Idx[idx + 2]], Ref->Lat[Ref->Idx[idx + 2]])) {

                        // Return coordinate as triangle index + barycentric coefficient
                        X[d] = idx + b[0];
                        Y[d] = idx + b[1];
                        break;
                    }
                }
            }
        } else {
            // Otherwise loop on all
            for(idx = 0; idx < Ref->NIdx; idx += 3) {
                if (Bary_Get(b, Ref->Wght ? Ref->Wght[idx / 3] : 0.0, Lon[d], Lat[d], Ref->AX[Ref->Idx[idx]], Ref->AY[Ref->Idx[idx]],
                    Ref->AX[Ref->Idx[idx + 1]], Ref->AY[Ref->Idx[idx + 1]], Ref->AX[Ref->Idx[idx + 2]], Ref->AY[Ref->Idx[idx + 2]])) {

                    // Return coordinate as triangle index + barycentric coefficient
                    X[d] = idx + b[0];
                    Y[d] = idx + b[1];
                    break;
                }
            }
        }
    }
    return 0;
}


//! Transform LatLon coordinates to XY for an M grid (Mesh or TIN)
int32_t GeoRef_XY2LL_M(
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
    for (i=0; i < Nb; i++) {
        Lon[i] = VertexValS(Ref->AX, NULL, Ref->NX, Ref->NY, X[i] - 1.0, Y[i] - 1.0, TRUE);
        Lat[i] = VertexValS(Ref->AY, NULL, Ref->NX, Ref->NY, X[i] - 1.0, Y[i] - 1.0, TRUE);
    }
    return 0;
}
