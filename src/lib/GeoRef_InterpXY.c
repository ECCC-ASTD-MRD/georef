#include <App.h>
#include "GeoRef.h"

int32_t GeoRef_XYInterp(
    const TGeoRef * const Ref,
    const TGeoOptions * const Opt,
    float * const zout,
    const float * const zin,
    double * const X,
    double * const Y,
    const int32_t npts
) {
    float * lzin  = NULL;
    const TGeoOptions * const opt = Opt ? Opt : &GeoRef_Options;

    if (Ref->Type & GRID_YINVERT) {
        lzin = (float *) malloc(Ref->NX * Ref->NY * sizeof(float));
        memcpy(lzin, zin, Ref->NX * Ref->NY * sizeof(float));
        f77name(permut)(lzin, &Ref->NX, &Ref->NY);
    } else {
        lzin = (float *) zin;
    }

    float * lxzin = NULL;
    if (Ref->Type & GRID_EXPAND) {
        lxzin = (float *) malloc(2 * Ref->NX * Ref->NY * sizeof(float));
        GeoRef_GridGetExpanded(Ref, opt, lxzin, lzin);
    } else  {
        lxzin = lzin;
    }

    GeoRef_InterpFinally(NULL, Ref, opt, zout, lxzin, X, Y, npts, NULL);

    if (lzin && lzin != zin) {
        free(lzin);
    }

    if (lxzin && lxzin != lzin) {
        free(lxzin);
    }

    return 0;
}


int32_t GeoRef_XYVal(
    const TGeoRef * const Ref,
    const TGeoOptions * const Opt,
    float * const zout,
    const float * const zin,
    double * const X,
    double * const Y,
    const int32_t n
) {
    const TGeoOptions * const opt = Opt ? Opt : &GeoRef_Options;
    const TGeoRef * const ref = GeoRef_SubGet(Ref);

    if (ref->NbSub > 0) {
        const TGeoRef * const yin_gd = ref->Subs[0];
        const TGeoRef * const yan_gd = ref->Subs[1];
        const int32_t sz = yin_gd->NX * yin_gd->NY;
        int32_t icode;
        double tmpy;
        for (int32_t j = 0; j < n; j++) {
            if (Y[j] > yin_gd->NY) {
                tmpy = Y[j]-yin_gd->NY;
                icode = GeoRef_XYInterp(yan_gd, opt, &zout[j], &zin[sz], X, &tmpy, 1);
            } else {
                tmpy = Y[j];
                icode = GeoRef_XYInterp(yin_gd, opt, &zout[j], zin, X, &tmpy, 1);
            }
        }
        return icode;
    } else {
        return GeoRef_XYInterp(ref, opt, zout, zin, X, Y, n);
    }
}


int32_t GeoRef_XYUVVal(
    const TGeoRef * const Ref,
    const TGeoOptions * const Opt,
    float * const uuout,
    float * const vvout,
    const float * const uuin,
    const float * const vvin,
    double * const X,
    double * const Y,
    const int32_t n
) {
    const TGeoOptions * const opt = Opt ? Opt : &GeoRef_Options;
    const TGeoRef * const ref = GeoRef_SubGet(Ref);

    if (ref->NbSub > 0) {
        const TGeoRef * const yin_gd = ref->Subs[0];
        const TGeoRef * const yan_gd = ref->Subs[1];
        const int32_t ni = yin_gd->NX;
        const int32_t nj = yin_gd->NY;
        double * const tmpy = (double*) malloc(n * sizeof(double));
        float * const uuyin = (float *) malloc(4 * n * sizeof(float));
        float * const vvyin = &uuyin[n];
        float * const uuyan = &uuyin[n * 2];
        float * const vvyan = &uuyin[n * 3];

        for (int32_t j = 0; j < n; j++) {
            if (Y[j] > yin_gd->NY) {
                tmpy[j] = Y[j]-yin_gd->NY;
            } else {
                tmpy[j] = Y[j];
            }
        }

        int32_t icode = GeoRef_XYUVVal(yin_gd, opt, uuyin, vvyin, uuin, vvin, X, tmpy, n);
        icode = GeoRef_XYUVVal(yan_gd, opt, uuyan, vvyan, &uuin[ni * nj], &vvin[ni * nj], X, tmpy, n);

        for (int32_t j = 0; j < n; j++) {
            if (Y[j] > yin_gd->NY) {
                uuout[j] = uuyan[j];
                vvout[j] = vvyan[j];
            } else {
                uuout[j] = uuyin[j];
                vvout[j] = vvyin[j];
            }
        }
        free(tmpy);
        free(uuyin);
        return icode;
    } else {
        TGeoOptions newOpt = *opt;
        newOpt.VectorMode = TRUE;
        newOpt.Symmetric = TRUE;
        GeoRef_XYInterp(ref, &newOpt, uuout, uuin, X, Y, n);

        newOpt.Symmetric = FALSE;
        GeoRef_XYInterp(ref, &newOpt, vvout, vvin, X, Y, n);
    }

    return 0;
}

int32_t GeoRef_XYWDVal(
    const TGeoRef * const Ref,
    const TGeoOptions * const Opt,
    float * const uuout,
    float * const vvout,
    const float * const uuin,
    const float * const vvin,
    double * const X,
    double * const Y,
    const int32_t n
) {
    double * const tmplat = (double*) malloc(2 * n * sizeof(double));
    double * const tmplon = &tmplat[n];
    float * const tmpuu = (float *) malloc(2 * n * sizeof(float));
    float * const tmpvv = &tmpuu[n];

    const TGeoOptions * const opt = Opt ? Opt : &GeoRef_Options;
    TGeoRef * const ref = GeoRef_SubGet(Ref);

    if (ref->NbSub > 0) {
        TGeoRef * const yin_gd = ref->Subs[0];
        TGeoRef * const yan_gd = ref->Subs[1];
        int32_t lni = yin_gd->NX;
        int32_t lnj = yin_gd->NY;
        double * const tmpy = (double *) malloc(n * sizeof(double));
        float * const uuyin = (float *) malloc(4 * n * sizeof(float));
        float * const vvyin = &uuyin[n];
        float * const uuyan = &uuyin[n * 2];
        float * const vvyan = &uuyin[n * 3];
        for (int32_t j = 0; j< n; j++) {
            if (Y[j] > yin_gd->NY) {
                tmpy[j] = Y[j]-yin_gd->NY;
            } else {
                tmpy[j] = Y[j];
            }
        }

        GeoRef_XYUVVal(yin_gd, opt, tmpuu, tmpvv, uuin, vvin, X, tmpy, n);
        GeoRef_XY2LL(yin_gd, tmplat, tmplon, X, tmpy, n, TRUE);
        GeoRef_UV2WD(yin_gd, uuyin, vvyin, tmpuu, tmpvv, tmplat, tmplon, n);

        GeoRef_XYUVVal(yan_gd, opt, tmpuu, tmpvv, &uuin[(lni*lnj)], &vvin[(lni*lnj)], X, tmpy, n);
        GeoRef_XY2LL(yan_gd, tmplat, tmplon, X, tmpy, n, TRUE);
        GeoRef_UV2WD(yan_gd, uuyan, vvyan, tmpuu, tmpvv, tmplat, tmplon, n);

        for (int32_t j = 0; j < n; j++) {
            if (Y[j] > yin_gd->NY) {
                uuout[j] = uuyan[j];
                vvout[j] = vvyan[j];
            } else {
                uuout[j] = uuyin[j];
                vvout[j] = vvyin[j];
            }
        }
        free(uuyin);
    } else {
        GeoRef_XYUVVal(ref, opt, tmpuu, tmpvv, uuin, vvin, X, Y, n);
        GeoRef_XY2LL(ref, tmplat, tmplon, X, Y, n, TRUE);
        GeoRef_UV2WD(ref, uuout, vvout, tmpuu, tmpvv, tmplat, tmplon, n);
    }

    free(tmplat);
    free(tmpuu);

    return 0;
}