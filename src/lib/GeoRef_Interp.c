#include <App.h>
#include "GeoRef.h"
#include "GeoRef_f.h"
#include "Triangle.h"

//! Predicate for grid rotation
int32_t GeoRef_IsRotated(
    const TGeoRef * const gr
) {
    //! \return 1 when rotated, 0 otherwise
    if (gr->RPNHeadExt.grref[0] == 'E') {
        if (fabs(gr->RPNHeadExt.xgref1 - gr->RPNHeadExt.xgref3) < 0.001) {
            // non rotated
            return 0;
        } else {
            // rotated
            return 1;
        }
    } else {
        return 0;
    }
    return 0;
}


//! Check if grids are compatible
int32_t c_gdcompatible_grids(
    const TGeoRef * const RefTo,
    const TGeoRef * const RefFrom
) {
    switch(RefFrom->GRTYP[0]) {
        case 'L':
        case 'A':
        case 'B':
        case 'G':
            return 0;

        // This is a fix from previous version that was never reached
        case 'Z':
            if (RefTo->RPNHeadExt.grref[0] == 'L') {
                return 0;
            } else if (RefTo->RPNHeadExt.grref[0] == 'E') {
                if (!GeoRef_IsRotated(RefTo)) {
                    return 0;
                } else {
                    return -1;
                }
            }
            break;

        default:
            return -1;
    }

    return 0;
}


int32_t GeoRef_InterpGridM(
    const TGeoRef * const Ref,
    const TGeoOptions * const Opt,
    float * const Out,
    const float * const In,
    const double * const X,
    const double * const Y,
    const int32_t Nb
) {
    Vect3d  b, v;
    int32_t d, n, idx;

    #pragma omp parallel for default(none) private(d, b, idx, n, v) shared(stderr, Nb, Ref, Opt, X, Y, Out, In)
    for(d = 0; d < Nb; d++) {
        if (X[d] >= 0 && Y[d] >= 0) {

            b[0] = X[d] - (int)X[d];
            b[1] = Y[d] - (int)Y[d];
            b[2] = 1.0 - b[0] - b[1];
            idx = (int)X[d] - 1;

            if(idx > Ref->NIdx) {
                Out[d] = Opt->NoData;
                continue;
            }

            if (Opt->Interp == IR_NEAREST) {
                n = Bary_Nearest(b);
                Out[d] = In[Ref->Idx[idx + n]];
            } else {
                v[0] = In[Ref->Idx[idx]];
                v[1] = In[Ref->Idx[idx + 1]];
                v[2] = In[Ref->Idx[idx + 2]];

                Out[d] = Bary_InterpV(b, v);
            }
        }
    }
    return Nb;
}


int32_t GeoRef_InterpFinally(
    TGeoRef * const RefTo,
    TGeoRef * const RefFrom,
    TGeoOptions * const Opt,
    float * const zout,
    const float * const zin,
    double * const X,
    double * const Y,
    const int32_t npts,
    TGeoSet * const GSet
) {
    // RefTo needed for type 4, 5 and Y grid
    //! \todo Check old type 4 and,
    //! \todo need to use new types (do not forget subllinear and subnearest)
    if ((!RefTo && (Opt->Interp == 4 || Opt->Interp == 5)) || !RefFrom) {
        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid georeference\n", __func__);
        return -1;
    }

    if (!X || !Y) {
        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Local coordinates not available\n", __func__);
        return -1;
    }

    // Fortran interpolation functions expect 1 to N
    double *x=NULL, *y=NULL;
    if (RefFrom->GRTYP[0]!='M') {
        x = (double*)malloc(npts * sizeof(double) * 2);
        y = &x[npts];

        for(int32_t j = 0; j < npts; j++) {
            x[j] = X[j] + 1.0;
            y[j] = Y[j] + 1.0;
        }
    }

    const int32_t old_degre_interp = Opt->Interp;

    switch(RefFrom->GRTYP[0]) {
        case 'M':
            GeoRef_InterpGridM(RefFrom, Opt, zout, zin, X, Y, npts);
            break;
        case '#':
        case 'Z':
        case 'G':
            switch (Opt->Interp) {
                case IR_NEAREST:
                if (!GSet) {
                    f77name(ez8_rgdint_0)(zout, x, y, &npts, zin, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, &Opt->NoData);
                } else {
                    if (GeoRef_SetEmptyIndex(GSet)) {
                        f77name(ez8_rgd_index_0)(GSet->Index, x, y, &npts, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2);
                    }
                    f77name(ez8_apply_0)(GSet->Index, zout, &npts, zin, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, &Opt->NoData);
                }
                break;

                case IR_LINEAR:
                if (!GSet) {
                    f77name(ez8_irgdint_1)(zout, x, y, &npts, zin, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, RefFrom->AX, RefFrom->AY, &RefFrom->Extension, &Opt->NoData);
                } else {
                    if (GeoRef_SetEmptyIndex(GSet)) {
                        f77name(ez8_irgd_index_1)(GSet->Index, x, y, &npts, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, RefFrom->AX, RefFrom->AY, &RefFrom->Extension);
                    }
                    f77name(ez8_apply_1)(GSet->Index, zout, &npts, zin, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, &Opt->NoData);
                }
                break;

                case IR_CUBIC:
                    f77name(ez8_irgdint_3)(zout, x, y, &npts, zin, &RefFrom->NX, &RefFrom->i1, &RefFrom->i2, &RefFrom->j1, &RefFrom->j2, RefFrom->AX, RefFrom->AY, RefFrom->NCX, RefFrom->NCY, &RefFrom->Extension, &Opt->NoData);
                   break;

                default:
                    Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid interpolation method\n", __func__);
                    return -1;
                    break;
            }
            break;

        case 'Y':
            if (RefFrom->NX > 1 && RefFrom->NY > 1 && Opt->Interp==IR_LINEAR) {
                int32_t un = 1;
                f77name(ez8_rgdint_1)(zout, x, y, &npts, zin, &RefFrom->NX, &un, &RefFrom->NY, 0, &Opt->NoData);
            } else {
                if (!GSet) {
                    Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: GeoSet not defined\n", __func__);
                    return -1;
                }
                f77name(ez_applywgts)(zout, GSet->wts, GSet->idx, zin, GSet->mask, &RefFrom->NX, &RefFrom->NY, &RefTo->NX, &RefTo->NY, &(GSet->n_wts), &Opt->NoData);
            }
            break;

        default:
            switch (Opt->Interp) {
                case IR_NEAREST:
                    if (!GSet) {
                        f77name(ez8_rgdint_0)(zout, x, y, &npts, zin, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, &Opt->NoData);
                    } else {
                        if (GeoRef_SetEmptyIndex(GSet)) {
                            f77name(ez8_rgd_index_0)(GSet->Index, x, y, &npts, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2);
                        }
                        f77name(ez8_apply_0)(GSet->Index, zout, &npts, zin, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, &Opt->NoData);
                    }
                    break;

                case IR_LINEAR:
                    if (!GSet) {
                        f77name(ez8_rgdint_1)(zout, x, y, &npts, zin, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, &RefFrom->Extension, &Opt->NoData);
                    } else {
                        if (GeoRef_SetEmptyIndex(GSet)) {
                            f77name(ez8_rgd_index_1)(GSet->Index, x, y, &npts, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, &RefFrom->Extension);
                        }
                        f77name(ez8_apply_1)(GSet->Index, zout, &npts, zin, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, &Opt->NoData);
                    }
                    break;

                case IR_CUBIC:
                    if (!GSet) {
                        f77name(ez8_rgdint_3)(zout, x, y, &npts, zin, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, &RefFrom->Extension, &Opt->NoData);
                    } else {
                        if (GeoRef_SetEmptyIndex(GSet)) {
                            f77name(ez8_rgd_index_3)(GSet->Index, x, y, &npts, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, &RefFrom->Extension);
                        }
                        f77name(ez8_apply_3)(GSet->Index, zout, &npts, zin, &RefFrom->NX, &RefFrom->j1, &RefFrom->j2, &Opt->NoData);
                    }
                    break;

                default:
                    Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid interpolation method\n", __func__);
                    return -1;
                    break;
            }
            break;
    }

    if (x) free(x);
    Opt->Interp = old_degre_interp;
    return 0;
}


//! Set each element of Buffer to NoData
int64_t GeoRef_InterpClear(
    //! [in] geo-reference
    const TGeoRef * const Ref,
    //! [in] Interpolation options
    const TGeoOptions * const Opt,
    //! [inout] Buffer to be cleared
    float * const Buffer
) {
    //! \return Number of elements in buffer

    int64_t i = 0;

    if (Opt->Extrap == ER_VALUE) {
        for(i=0; i < Ref->NX * Ref->NY; i++) {
            Buffer[i] = Opt->NoData;
        }
    }
    return i;
}


//! Interpolate values between 2 georeferences
int32_t GeoRef_Interp(
    //! [in] Destination geo-reference
    const TGeoRef * const RefTo,
    //! [in] Source geo-reference
    const TGeoRef * const RefFrom,
    //! [in] Interpolation options
    const TGeoOptions * const Opt,
    //! [out] Destination interpolated values
    float * const zout,
    //! [in] Source values
    const float * const zin
) {
    //! \return FALSE (0) if operation failed, TRUE (1) otherwise
    if (!RefFrom || !RefTo) {
        Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s: Source or target grid undefined\n", __func__);
        return FALSE;
    }

    if (RefFrom == RefTo) {
        memcpy(zout, zin, RefFrom->NX * RefFrom->NY * sizeof(float));
        return TRUE;
    }

    const TGeoOptions * const opt = Opt ? Opt : &GeoRef_Options;
    TApp_Timer * const int_timer = App_TimerCreate();
    App_TimerStart(int_timer);
    GeoRef_InterpClear(RefTo, opt, zout);

    // If we're using pre-calculated index weight
    if (opt->Interp==IR_WEIGHTINDEX) {
       return(GeoRef_InterpWeight(RefTo,RefFrom,opt,zout,NULL,zin,NULL));
    }

    int32_t ok = TRUE;
    if (RefFrom->NbSub > 0 || RefTo->NbSub > 0) {
        // YY mutli grids involved
        return GeoRef_InterpYY(RefTo, RefFrom, opt, zout, zin);
    } else {
        TGeoSet * gset = GeoRef_SetGet(RefTo, RefFrom, opt);
        float * lzin  = NULL;
        float * lxzin = NULL;
        if (RefFrom->Type & GRID_YINVERT) {
            lzin = (float *)malloc(RefFrom->NX * RefFrom->NY * sizeof(float));
            memcpy(lzin, zin, RefFrom->NX * RefFrom->NY * sizeof(float));
            f77name(permut)(lzin, &RefFrom->NX, &RefFrom->NY);
        } else {
            lzin = zin;
        }

        if (RefFrom->Type&GRID_EXPAND) {
            lxzin = (float *)malloc(2 * RefFrom->NX * RefFrom->NY * sizeof(float));
            GeoRef_GridGetExpanded(RefFrom, opt, lxzin, lzin);
        } else {
            lxzin = lzin;
        }

        if (GeoRef_CalcLL(RefTo)) {
            GeoRef_SetCalcXY(gset);

            if (GeoRef_InterpFinally(RefTo, RefFrom, opt, zout, lxzin, gset->X, gset->Y, RefTo->NX * RefTo->NY, gset) == 0) {
                if (opt->PolarCorrect) {
                    GeoRef_SetZoneDefine(gset);
                    GeoRef_CorrectValue(gset, zout, lxzin);
                }
            } else {
                ok = FALSE;
            }
        }

        if (lzin && lzin != zin) {
            free(lzin);
        }

        if (lxzin && lxzin != lzin && lxzin != zin) {
            free(lxzin);
        }
    }
    App_TimerStop(int_timer);
    Lib_Log(APP_LIBGEOREF, APP_STAT, "%s: Interpolation took \033[1; 32m%.3f ms\033[0m\n", __func__, App_TimerTotalTime_ms(int_timer));

    return ok;
}


int32_t GeoRef_InterpYY(
    const TGeoRef * const RefTo,
    const TGeoRef * const RefFrom,
    const TGeoOptions * const Opt,
    float * const zout,
    const float * const zin
) {
    const TGeoOptions * const opt = Opt ? Opt : &GeoRef_Options;

    // Setup for input grid
    const int32_t yyin = (RefFrom->NbSub > 0) ? 1 : 0;
    const TGeoRef * const yin_gdin = (RefFrom->NbSub > 0) ? RefFrom->Subs[0] : RefFrom;
    const TGeoRef * const yan_gdin = (RefFrom->NbSub > 0) ? RefFrom->Subs[1] : NULL;

    // Setup for output grid
    const int32_t yyout = (RefTo->NbSub > 0) ? 1 : 0;
    const TGeoRef * const yin_gdout = (RefTo->NbSub > 0) ? RefTo->Subs[0] : RefTo;
    const TGeoRef * const yan_gdout = (RefTo->NbSub > 0) ? RefTo->Subs[1] : NULL;

    const int32_t ni = yin_gdout->NX;
    const int32_t nj = yin_gdout->NY;

    // Interp input one grid to yygrid - no masking needed
    if (!yyin && yyout) {
        if (GeoRef_Interp(yin_gdout, RefFrom, opt, zout, zin) && GeoRef_Interp(yan_gdout, RefFrom, opt, &zout[ni * nj], zin)) {
            return TRUE;
        } else {
            return FALSE;
        }
    }

    // check if one input sub grid is identical to dest sub grid or dest single grid
    if (yin_gdin == RefTo) {
        return GeoRef_Interp(RefTo, yin_gdin, opt, zout, zin);
    }
    if (yan_gdin == RefTo) {
        return GeoRef_Interp(RefTo, yan_gdin, opt, zout, &zin[(yin_gdin->NX)*(yin_gdin->NY)]);
    }

    /* User specifies to use 1 subgrid for interpolation ezsetopt(USE_1SUBGRID) */
    /* User must specify the sub grid value in ezsetival(SUBGRIDID) */
    /* This is only appropriate if the destination grid is non yin-yang grid */

    if (RefFrom->Sub >= 0) { // User specifies to use 1 grid only
        // Output is a Yin-Yang grid
        if (yyout) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Cannot use subgrid to interpolate to a Yin-Yang grid\n", __func__);
            return FALSE;
        }
        // Is specified subgrid within the subgrid list
        if (RefFrom->Sub >= RefFrom->NbSub) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid subgrid: %i\n", __func__, RefFrom->Sub);
            return FALSE;
        }
        // Use yin input grid
        if (RefFrom->Sub == 0) {
            return GeoRef_Interp(yin_gdout, yin_gdin, opt, zout, zin);
        }

        // Use yang input grid
        if (RefFrom->Sub == 1) {
            return GeoRef_Interp(yin_gdout, yan_gdin, opt, zout, &zin[(yin_gdin->NX)*(yin_gdin->NY)]);
        }
    }

    // To use both Yin and Yang grids in Yin-yang input grid
    // Masquer les grilles YY input pour enlever overlap et calculer les X, Y
    TGeoSet * const gset = GeoRef_SetGet(RefTo, RefFrom, opt);
    GeoRef_SetCalcYYXY(gset);

    if (!yyout) {
        // Interp yinyang to one grid

        int32_t yincount_yin = gset->yincount_yin;
        int32_t yancount_yin = gset->yancount_yin;
        float * const yin2yin_zvals = (float *) malloc(yincount_yin * sizeof(float));
        float * const yan2yin_zvals = (float *) malloc(yancount_yin * sizeof(float));

        GeoRef_XYVal(yin_gdin, opt, yin2yin_zvals, zin, gset->yin2yin_x, gset->yin2yin_y, gset->yincount_yin);
        GeoRef_XYVal(yan_gdin, opt, yan2yin_zvals, &zin[(yin_gdin->NX) * (yin_gdin->NY)], gset->yan2yin_x, gset->yan2yin_y, gset->yancount_yin);

        yincount_yin = 0;
        yancount_yin = 0;
        for(int32_t k = 0; k < ni * nj; k++) {
            if (gset->yin_maskout[k] == 1.0) {
                zout[k] = yan2yin_zvals[yancount_yin];
                yancount_yin++;
            } else {
                zout[k] = yin2yin_zvals[yincount_yin];
                yincount_yin++;
            }
        }
        free(yin2yin_zvals);
        free(yan2yin_zvals);
        return TRUE;
    } else {
        // Interp yinyang to yinyang

        int32_t yincount_yin = gset->yincount_yin;
        int32_t yancount_yin = gset->yancount_yin;
        int32_t yincount_yan = gset->yincount_yan;
        int32_t yancount_yan = gset->yancount_yan;
        float * const yin2yin_zvals = (float *) malloc(yincount_yin * sizeof(float));
        float * const yan2yin_zvals = (float *) malloc(yancount_yin * sizeof(float));
        float * const yin2yan_zvals = (float *) malloc(yincount_yan * sizeof(float));
        float * const yan2yan_zvals = (float *) malloc(yancount_yan * sizeof(float));

        GeoRef_XYVal(yin_gdin, opt, yin2yin_zvals, zin, gset->yin2yin_x, gset->yin2yin_y, gset->yincount_yin);
        GeoRef_XYVal(yan_gdin, opt, yan2yin_zvals, &zin[(yin_gdin->NX) * (yin_gdin->NY)], gset->yan2yin_x, gset->yan2yin_y, gset->yancount_yin);
        GeoRef_XYVal(yin_gdin, opt, yin2yan_zvals, zin, gset->yin2yan_x, gset->yin2yan_y, gset->yincount_yan);
        GeoRef_XYVal(yan_gdin, opt, yan2yan_zvals, &zin[(yin_gdin->NX) * (yin_gdin->NY)], gset->yan2yan_x, gset->yan2yan_y, gset->yancount_yan);

        // Interp input YY grid to Yin grid
        yincount_yin = 0; yancount_yin = 0;
        for(int32_t k = 0; k < ni * nj; k++) {
            if (gset->yin_maskout[k] == 1.0) {
                zout[k] = yan2yin_zvals[yancount_yin];
                yancount_yin++;
            } else {
                zout[k] = yin2yin_zvals[yincount_yin];
                yincount_yin++;
            }
        }

        // Interp input YY grid to Yang grid
        yincount_yan = 0; yancount_yan = 0;
        for(int32_t k = 0; k < ni * nj; k++) {
            if (gset->yan_maskout[k] == 1.0) {
                zout[k + (ni * nj)] = yan2yan_zvals[yancount_yan];
                yancount_yan++;
            } else {
                zout[k + (ni * nj)] = yin2yan_zvals[yincount_yan];
                yincount_yan++;
            }
        }

        free(yin2yin_zvals);
        free(yan2yin_zvals);
        free(yin2yan_zvals);
        free(yan2yan_zvals);
    }

    return TRUE;
}

int32_t GeoRef_InterpWeight(
    //! [in] Destination geo-reference
    const TGeoRef * const RefTo,
    //! [in] Source geo-reference
    const TGeoRef * const RefFrom,
    //! [in] Interpolation options
    const TGeoOptions * const Opt,
    //! [out] Destination interpolated values U
    float * const zuout,
    //! [out] Destination interpolated values V (optional)
    float * const zvout,
    //! [in] Source values U
    const float * const zuin,
    //! [in] Source values V (optional)
    const float * const zvin
) {

   TGeoSet *gset = NULL;
   int32_t  i, j, pi, pj, n;
   uint32_t idx;
   float    uval, vval, uvalr, vvalr, val1,dp,cosa,sina;
   float   *ip = NULL;
   float   *rp = NULL;

   const TGeoOptions * const opt = Opt ? Opt : &GeoRef_Options;

   if (!RefTo) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid destination\n", __func__);
      return FALSE;
   }
   if (!RefFrom) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid source\n", __func__);
      return FALSE;
   }

   gset = GeoRef_SetGet(RefTo, RefFrom, opt);
   if (!gset->Index || gset->Index[0] == REF_INDEX_EMPTY) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Weight index not defined\n", __func__);
      return FALSE;
   }
 
    // Clear destination buffer
//    if (opt->Extrap != ER_VALUE)
//       memset(zout,0x0,RefFrom->NX*RefFrom->NY*sizeof(float));
    
    ip = gset->Index;
    rp = gset->R;

    // As long as the file or the list is not empty
    while(*ip != REF_INDEX_END) {

        // Get destination gridpoint
        i = *(ip++);
        j = *(ip++);
        idx= j*RefTo->NX+i;
        uval = vval = 0.0;
        
        // Loop on contributing gritpoints
        while(*ip != REF_INDEX_SEPARATOR) {
            pi = *(ip++);   // Source gridpoint
            pj = *(ip++);
            dp = *(ip++);   // Fraction of value to use

            val1=zuin[pj*RefFrom->NX+pi];
            if (DATA_ISVALID(val1,opt->NoData)) {
                uval += val1*dp;

                if (zvin) {
                   val1=zvin[pj*RefFrom->NX+pi];
                   vval += val1*dp;
                }
            }
        }

        // Rotate components
        cosa = *(rp++);
        sina = *(rp++);
        uvalr = cosa*uval - sina*vval;
        vvalr = sina*uval + cosa*vval;

        // Check for valid previous value and average if so (we suppose 2 values (Yin/Yang)
        val1 = zuout[idx];
        if (DATA_ISVALID(val1,opt->NoData)) {
            uvalr = (uvalr + val1)/2.0;

            if (zvin) {
                val1 = zvout[idx];
                vvalr = (vvalr + val1)/2.0;
            }
        }

        zuout[idx]=uvalr;

        if (zvin) {
            zvout[idx]=vvalr;
        }

        // Skip separators
        ip++;
        rp++;
    }

    return(TRUE);
}
