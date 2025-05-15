//! \file

#include <App.h>
#include "GeoRef.h"

static int32_t GeoRef_SetMax = -64;    ///< How many set to kepp in memory

//! Free gridset zone definitions
void GeoRef_SetZoneFree(
    //! [in, out] Grid set pointer
    TGeoSet * const GSet
) {
    for (int32_t i = 0; i < SET_NZONES; i++) {
        if (GSet->zones[i].npts > 0) {
            free(GSet->zones[i].idx);
            free(GSet->zones[i].x);
            free(GSet->zones[i].y);
            GSet->zones[i].npts = 0;
            GSet->zones[i].idx  = NULL;
            GSet->zones[i].x    = NULL;
            GSet->zones[i].y    = NULL;
        }
    }
}

//! Finds the points at a pole
int32_t GeoRef_SetZoneDefinePole(
    //! [in, out] Grid set
    TGeoSet * const GSet,
    //! [in] Zone identifier (GRID_NORTH_POLE, GRID_SOUTH_POLE)
    const int32_t Zone
) {
    //! \returns Always 0
    TGeoZone * const zone = &GSet->zones[Zone];

    const int32_t nbpts = (GSet->RefTo->NX * GSet->RefTo->NY);

    // On commence par trouver les points au pole
    int32_t * const tmpidx = (int*)malloc(nbpts * sizeof(int));

    zone->npts = 0;
    double xpole, ypole;
    if (GSet->RefFrom->GRTYP[0] == 'Z' && GSet->RefFrom->RPNHeadExt.grref[0] == 'E') {
        xpole = 0.5 * GSet->RefFrom->NX;
        ypole = (Zone == GRID_NORTH) ? GSet->RefFrom->NY + 0.5 : 0.5;
    } else {
        double latpole = (Zone == GRID_NORTH) ? 90.0 : -90.0;
        double lonpole = 0.0;
        GeoRef_LL2XY(GSet->RefFrom, &xpole, &ypole, &latpole, &lonpole, 1, TRUE);
    }

    for (int32_t i = 0; i < nbpts; i++) {
        if (fabs(GSet->Y[i] - ypole) < 1.0e-3) {
            tmpidx[zone->npts] = i;
            zone->npts++;
        }
    }

    if (zone->npts > 0) {
        zone->x = (double*)malloc(zone->npts * sizeof(double));
        zone->y = (double*)malloc(zone->npts * sizeof(double));
        zone->idx = (int32_t *)malloc(zone->npts * sizeof(int));
        Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s: Number of points at pole: %d\n", __func__, zone->npts);

        for (int32_t i = 0; i < zone->npts; i++) {
            zone->x[i]   = GSet->X[tmpidx[i]];
            zone->y[i]   = GSet->Y[tmpidx[i]];
            zone->idx[i] = tmpidx[i];
        }
    }

    free(tmpidx);

    return 0;
}


//! Finds the points between the pole and a limit
int32_t GeoRef_SetZoneDefineThem(
    //! [in,out ] Grid set
    TGeoSet * const GSet,
    //! [in] Zone identifier (GRID_NORTH, GRID_SOUTH)
    const int32_t Zone
) {
    //! \return Always 0
    TGeoZone * const zone = &GSet->zones[Zone];

    const int32_t nbpts = (GSet->RefTo->NX * GSet->RefTo->NY);
    int32_t * const tmpidx = (int*) malloc(nbpts * sizeof(int));

    zone->npts = 0;
    int32_t jlim = (Zone == GRID_SOUTH) ? GSet->RefFrom->j1 + 1 : GSet->RefFrom->j2 - 2;
    for (int32_t i = 0; i < nbpts; i++) {
        if ((Zone == GRID_SOUTH) ? ((int)GSet->Y[i] < jlim) : ((int)GSet->Y[i] > jlim)) {
            tmpidx[zone->npts] = i;
            zone->npts++;
        }
    }

    if (zone->npts > 0) {
        zone->x = (double*)malloc(zone->npts * sizeof(double));
        zone->y = (double*)malloc(zone->npts * sizeof(double));
        zone->idx = (int32_t *) malloc(zone->npts * sizeof(int));
        Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s: Number of points between pole and limit: %d\n", __func__, zone->npts);

        for (int32_t i = 0; i < zone->npts; i++) {
            zone->x[i]   = GSet->X[tmpidx[i]];
            zone->y[i]   = GSet->Y[tmpidx[i]];
            zone->idx[i] = tmpidx[i];
        }
    }

    free(tmpidx);

    return 0;
}


//! Finds the points outside of the source data
int32_t GeoRef_SetZoneDefineOut(
    //! [in, out] Grid set
    TGeoSet * const GSet,
    //! [in] Zone identifier (GRID_OUTSIDE)
    const int32_t Zone
) {
    //! \return Always 0
    TGeoZone * const zone = &GSet->zones[Zone];

    const int32_t nbpts = (GSet->RefTo->NX * GSet->RefTo->NY);
    int32_t * const tmpidx = (int32_t*)malloc(nbpts * sizeof(int32_t));

    int32_t offsetright = 0;
    int32_t offsetleft = 0;

    /* Was commented in ezscint also
    if (groptions.degre_interp == CUBIC) {
        offsetright = 2;
        offsetleft = 1;
    } else {
        offsetright = 0;
        offsetleft = 0;
    }
    */

    zone->npts = 0;
    for (int32_t i = 0; i < nbpts; i++) {
        const int32_t ix = (int32_t)(GSet->X[i] + 0.5);
        const int32_t iy = (int32_t)(GSet->Y[i] + 0.5);
        if (ix < (1 + offsetleft) || iy < (1 + offsetleft) || ix > (GSet->RefFrom->NX - offsetright) || iy > (GSet->RefFrom->NY - offsetright)) {
            tmpidx[zone->npts] = i;
            zone->npts++;
        }
    }

    if (zone->npts>0) {
        zone->x = (double*)malloc(zone->npts * sizeof(double));
        zone->y = (double*)malloc(zone->npts * sizeof(double));
        zone->idx = (int32_t *) malloc(zone->npts * sizeof(int32_t));
        Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s: Number of outside points: %i \n", __func__, zone->npts);

        for (int32_t i = 0; i < zone->npts; i++) {
            zone->x[i]   = GSet->X[tmpidx[i]];
            zone->y[i]   = GSet->Y[tmpidx[i]];
            zone->idx[i] = tmpidx[i];
        }
    }

    free(tmpidx);

    return 0;
}


//! Define the various zones
int32_t GeoRef_SetZoneDefine(
    //! [in, out] Grid set
    TGeoSet *GSet
) {
    //! \return Always 0
    if (!GSet || (GSet->flags & SET_ZONES)) {
        return 0;
    }

    int32_t extrap = FALSE;
    switch (GSet->RefFrom->GRTYP[0]) {
        case 'M':
            GSet->flags |= SET_ZONES;
            return 0;
        case 'N':
        case 'S':
        case '!':
            extrap = TRUE;
            break;

        case 'L':
            if (GSet->RefFrom->Extension == 0) {
                extrap = TRUE;
            }
            break;

        case '#':
        case 'Z':
        case 'Y':
            switch(GSet->RefFrom->RPNHeadExt.grref[0]) {
                case 'N':
                case 'S':
                case 'L':
                    extrap = TRUE;
                    break;

                case 'E':
                    if (358.0 > (GSet->RefFrom->AX[GSet->RefFrom->NX-1] - GSet->RefFrom->AX[0])) {
                        extrap = TRUE;
                    }
                    break;
            }
            break;
    }

    for (int32_t i = 0; i < SET_NZONES; i++) {
        GSet->zones[i].npts = 0;
    }

    if (extrap) {
        GeoRef_SetZoneDefineOut(GSet, GRID_OUTSIDE);
    } else {
        GeoRef_SetZoneDefinePole(GSet, GRID_NORTH_POLE);
        GeoRef_SetZoneDefinePole(GSet, GRID_SOUTH_POLE);
        GeoRef_SetZoneDefineThem(GSet, GRID_SOUTH);
        GeoRef_SetZoneDefineThem(GSet, GRID_NORTH);
    }

    GSet->flags |= SET_ZONES;
    return 0;
}


//! Calculate XY correspondance of destination points within the source grid
int32_t GeoRef_SetCalcXY(
    //! [in, out] Grid set
    TGeoSet *GSet
) {
    //! \return Always 0
    if (GSet) {
        const int32_t size = GSet->RefTo->NX * GSet->RefTo->NY;
        if (!GSet->X) {
            GSet->X = (double*)calloc(size * 2, sizeof(double));
            GSet->Y = &GSet->X[size];

            GeoRef_LL2XY(GSet->RefFrom, GSet->X, GSet->Y, GSet->RefTo->Lat, GSet->RefTo->Lon, size, TRUE);
        }

        GeoRef_SetIndexInit(GSet);
    }
    return 0;
}


//! Configure initial index
int32_t GeoRef_SetIndexInit(
    //! [in, out] Grid set
    TGeoSet * const GSet
) {
    if (GSet) {
        int32_t mult = 1;
        if (!GSet->Index || GSet->IndexMethod != GSet->Opt.Interp) {
            GSet->IndexMethod = GSet->Opt.Interp;
            if (GSet->IndexMethod == IR_AVERAGE || GSet->IndexMethod == IR_VECTOR_AVERAGE || GSet->IndexMethod == IR_CONSERVATIVE || GSet->IndexMethod == IR_NORMALIZED_CONSERVATIVE) {
                mult = 1024;
                const char * const str = getenv("GEOREF_INDEX_SIZE_HINT");
                if (str) {
                    mult = atoi(str);
                }
            } else {
                mult = GSet->IndexMethod == IR_CUBIC ? 10 : (GSet->IndexMethod == IR_LINEAR ? 6 : 2);
            }
            // Set the index size to the number of items without the multiplicator
            // the index size will be readjusted when in CONSERVATIVE modes otherwise, the multiplicator will be put in nj when writing
            GSet->IndexSize = GSet->RefTo->NX * GSet->RefTo->NY;
            GSet->Index = (float*)realloc(GSet->Index, GSet->IndexSize * mult * sizeof(float));
            GSet->Index[0] = REF_INDEX_EMPTY;
        }
        return GSet->IndexSize * mult;
    }
    return 0;
}


//! Calculate XY correspondance of destination points within the source grid (for YY grids)
int32_t GeoRef_SetCalcYYXY(
    //! [in, out] Grid set
    TGeoSet * const GSet
) {
    //! \return Error code of the last operation
    //  Need only access to either yin or Yang info for the lat and lon val
    if (!GSet || (GSet->flags & SET_YYXY)) {
        return 0;
    }

    TGeoRef * yin_ref = GSet->RefFrom->Subs[0];
    TGeoRef * yan_ref = GSet->RefFrom->Subs[1];

    TGeoRef * yin_gdout, *yan_gdout;
    int32_t ni, nj;
    // Check what the destination grid is
    if (GSet->RefTo->NbSub > 0) {
        yin_gdout = GSet->RefTo->Subs[0];
        yan_gdout = GSet->RefTo->Subs[1];
        ni = yin_gdout->NX;
        nj = yin_gdout->NY;
    } else {
        yin_gdout = GSet->RefTo;
        ni = GSet->RefTo->NX;
        nj = GSet->RefTo->NY;
    }

    int32_t nij = ni * nj;

    // Masquer les grilles YY input pour enlever overlap si TRUE
    GSet->yin2yin_lat = (double*)malloc(4 * nij * sizeof(double));
    GSet->yin2yin_lon = (double*)malloc(4 * nij * sizeof(double));
    GSet->yan2yin_lat = (double*)malloc(4 * nij * sizeof(double));
    GSet->yan2yin_lon = (double*)malloc(4 * nij * sizeof(double));

    // Create mask with Yin as a priority choice and store x, y, lat, lon pos
    GSet->yin_maskout = (float *) malloc(nij * sizeof(float));
    GSet->yinlat = (double*) malloc(nij * sizeof(double));
    GSet->yinlon = (double*) malloc(nij * sizeof(double));
    GeoRef_GetLL(yin_gdout, GSet->yinlat, GSet->yinlon);
    int32_t yancount_yin = 0;
    int32_t yincount_yin = 0;
    GeoRef_MaskYYApply(yin_gdout, yin_ref, &GSet->Opt, ni, nj, GSet->yin_maskout, GSet->yinlat, GSet->yinlon,
        GSet->yin2yin_lat, GSet->yin2yin_lon, &yincount_yin, GSet->yan2yin_lat, GSet->yan2yin_lon, &yancount_yin);

    // Store the lats and lons
    GSet->yincount_yin = yincount_yin;
    GSet->yancount_yin = yancount_yin;
    GSet->yin2yin_lat = (double*) realloc(GSet->yin2yin_lat, yincount_yin * sizeof(double));
    GSet->yin2yin_lon = (double*) realloc(GSet->yin2yin_lon, yincount_yin * sizeof(double));
    GSet->yan2yin_lat = (double*) realloc(GSet->yan2yin_lat, yancount_yin * sizeof(double));
    GSet->yan2yin_lon = (double*) realloc(GSet->yan2yin_lon, yancount_yin * sizeof(double));

    // Store the Xs and Ys
    GSet->yin2yin_x = (double*) malloc(yincount_yin * sizeof(double));
    GSet->yin2yin_y = (double*) malloc(yincount_yin * sizeof(double));
    GSet->yan2yin_x = (double*) malloc(yancount_yin * sizeof(double));
    GSet->yan2yin_y = (double*) malloc(yancount_yin * sizeof(double));
    GeoRef_LL2XY(yin_ref, GSet->yin2yin_x, GSet->yin2yin_y, GSet->yin2yin_lat, GSet->yin2yin_lon, GSet->yincount_yin, TRUE);
    int32_t icode = GeoRef_LL2XY(yan_ref, GSet->yan2yin_x, GSet->yan2yin_y, GSet->yan2yin_lat, GSet->yan2yin_lon, GSet->yancount_yin, TRUE);

    // If destination grid is YY
    if (GSet->RefTo->NbSub > 0) {
        GSet->yin2yan_lat = (double*)malloc(4 * nij * sizeof(double));
        GSet->yin2yan_lon = (double*)malloc(4 * nij * sizeof(double));
        GSet->yan2yan_lat = (double*)malloc(4 * nij * sizeof(double));
        GSet->yan2yan_lon = (double*)malloc(4 * nij * sizeof(double));

        // Create mask (Yin priority) with src Yin, src Yang onto dest Yang and store x, y pos
        GSet->yan_maskout = (float *) malloc(nij * sizeof(float));
        GSet->yanlat = (double*) malloc(nij * sizeof(double));
        GSet->yanlon = (double*) malloc(nij * sizeof(double));
        icode = GeoRef_GetLL(yan_gdout, GSet->yanlat, GSet->yanlon);
        int32_t yancount_yan, yincount_yan;
        icode = GeoRef_MaskYYApply(yan_gdout, yin_ref, &GSet->Opt, ni, nj, GSet->yan_maskout, GSet->yanlat, GSet->yanlon, GSet->yin2yan_lat, GSet->yin2yan_lon, &yincount_yan, GSet->yan2yan_lat, GSet->yan2yan_lon, &yancount_yan);
        GSet->yincount_yan = yincount_yan;
        GSet->yancount_yan = yancount_yan;
        GSet->yin2yan_lat = (double*) realloc(GSet->yin2yan_lat, yincount_yan * sizeof(double));
        GSet->yin2yan_lon = (double*) realloc(GSet->yin2yan_lon, yincount_yan * sizeof(double));
        GSet->yan2yan_lat = (double*) realloc(GSet->yan2yan_lat, yancount_yan * sizeof(double));
        GSet->yan2yan_lon = (double*) realloc(GSet->yan2yan_lon, yancount_yan * sizeof(double));

        // Store the Xs and Ys
        GSet->yin2yan_x = (double*) malloc(yincount_yan * sizeof(double));
        GSet->yin2yan_y = (double*) malloc(yincount_yan * sizeof(double));
        GSet->yan2yan_x = (double*) malloc(yancount_yan * sizeof(double));
        GSet->yan2yan_y = (double*) malloc(yancount_yan * sizeof(double));
        icode = GeoRef_LL2XY(yin_ref, GSet->yin2yan_x, GSet->yin2yan_y, GSet->yin2yan_lat, GSet->yin2yan_lon, yincount_yan, TRUE);
        icode = GeoRef_LL2XY(yan_ref, GSet->yan2yan_x, GSet->yan2yan_y, GSet->yan2yan_lat, GSet->yan2yan_lon, yancount_yan, TRUE);
    }

    GSet->flags |= SET_YYXY;

    return icode;
}

//! Free grid set
void GeoRef_SetFree(
    //! [in, out] Grid set
    TGeoSet * const GSet
) {
    //! \todo Check to free
    if (GSet->RefFrom) {
        GSet->RefFrom == NULL;
    }

    if (GSet->Index) {
        free(GSet->Index);
        GSet->Index = NULL;
        GSet->IndexMethod = IR_UNDEF;
    }
    if (GSet->X) {
        free(GSet->X);
        GSet->X = NULL;
        GSet->Y = NULL;
    }

    GeoRef_SetZoneFree(GSet);

    //TODO: to free:
    //   int32_t *mask_in, *mask_out;
    // float *yin_maskout, *yan_maskout;
    // float *yinlat, *yinlon, *yanlat, *yanlon;
    // float *yin2yin_lat, *yin2yin_lon, *yan2yin_lat, *yan2yin_lon;
    // float *yin2yan_lat, *yin2yan_lon, *yan2yan_lat, *yan2yan_lon;
    // float *yin2yin_x, *yin2yin_y, *yan2yin_x, *yan2yin_y;
    // float *yin2yan_x, *yin2yan_y, *yan2yan_x, *yan2yan_y;
}


//! Read a gridset definition and index from a file
TGeoSet* GeoRef_SetReadFST(
    //! [in] Destination grid reference
    const TGeoRef * const RefTo,
    //! [in] Source grid reference
    const TGeoRef * const RefFrom,
    //! [in] Interpolation type (2=bilinear, 3=bicubic, -1=any)
    const int32_t InterpType,
    //! [in] FSTD file pointer
    const fst_file * const File
) {
    //! \return Grid set pointer or NULL
    TGeoSet * gset = GeoRef_SetGet(RefTo, RefFrom, NULL);
    if (!gset) return NULL;

    if (File) {
        // Rechercher et lire l'information de l'enregistrement specifie
        fst_record crit = default_fst_record;
        strncpy(crit.etiket, "GRIDSET", FST_ETIKET_LEN);
        strncpy(crit.nomvar, "####", FST_NOMVAR_LEN);
        crit.typvar[0]=RefTo->GRTYP[0];
        crit.typvar[1]=RefFrom->GRTYP[0];
        crit.grtyp[0]='\0';
        crit.ip3 = InterpType;

        fst_record record = default_fst_record;;
        if (fst24_read(File, &crit, NULL, &record) != TRUE) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Could not find gridset index field (fst24_read failed)\n", __func__);
            return NULL;
        }
        gset->IndexMethod = (TRef_Interp)InterpType;
        gset->Index = (float*)record.data;
        record.data = NULL;

        strncpy(crit.nomvar, "#>>#", FST_NOMVAR_LEN);
        if (fst24_read(File, &crit, NULL, &record) != TRUE) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Could not find gridset longitude field (fst24_read failed)\n", __func__);
            return NULL;
        }
        gset->X = record.data;
        record.data = NULL;

        strncpy(crit.nomvar, "#^^#", FST_NOMVAR_LEN);
        if (fst24_read(File, &crit, NULL, &record) != TRUE) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Could not find gridset longitude field (fst24_read failed)\n", __func__);
            return NULL;
        }
        gset->Y = record.data;
        record.data = NULL;
    }

    return gset;
}


//! Write a gridset definition and index from a file
int32_t GeoRef_SetWriteFST(
    //! [in] Grid set
    const TGeoSet * const GSet,
    //! [in, out] FSTD file pointer
    fst_file * const File
) {
    //! \return TRUE on success, FALSE otherwise
    if (GeoRef_SetHasIndex(GSet)) {
        Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s:  Writing index (%ix%i)\n", __func__,
            0, GSet->IndexMethod == IR_CUBIC ? 10 : (GSet->IndexMethod == IR_LINEAR ? 6 : 1));
        fst_record record = default_fst_record;

        record.nk   = 1;
        record.dateo= 0;
        record.deet = 0;
        record.npas = 0;
        record.ip1  = 0;
        record.ip2  = 0;
        record.ip3  = GSet->IndexMethod;
        record.ig1   = 0;
        record.ig2   = 0;
        record.ig3   = 0;
        record.ig4   = 0;
        record.typvar[0] = GSet->RefTo->GRTYP[0];
        record.typvar[1] = GSet->RefFrom->GRTYP[0];
        strncpy(record.etiket, "GRIDSET", FST_ETIKET_LEN);
        strncpy(record.grtyp, "X", FST_GTYP_LEN);
        record.data_type = FST_TYPE_REAL_IEEE;

        if (GSet->X) {
            record.data = GSet->X;
            record.data_bits = 64;
            record.pack_bits = 64;
            record.ni   = GSet->RefTo->NX;
            record.nj   = GSet->RefTo->NY;
            strncpy(record.nomvar, "#>>#", FST_NOMVAR_LEN);
            if (fst24_write(File, &record, FST_SKIP) <= 0) {
                Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Could not write gridset index field (fst24_write failed)\n", __func__);
                return FALSE;
            }
        }

        if (GSet->Y) {
            record.data = GSet->Y;
            record.data_bits = 64;
            record.pack_bits = 64;
            record.ni   = GSet->RefTo->NX;
            record.nj   = GSet->RefTo->NY;
            strncpy(record.nomvar, "#^^#", FST_NOMVAR_LEN);
            if (fst24_write(File, &record, FST_SKIP) <= 0) {
                Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Could not write gridset index field (fst24_write failed)\n", __func__);
                return FALSE;
            }
        }

        if (GSet->Index) {
            record.data = GSet->Index;
            record.data_bits = 32;
            record.pack_bits = 32;
            record.ni   = GSet->IndexSize;
            record.nj   = GSet->IndexMethod == IR_CUBIC ? 10 : (GSet->IndexMethod == IR_LINEAR ? 6 : 1);
            record.ip3  = GSet->IndexMethod;
            strncpy(record.nomvar, "####", FST_NOMVAR_LEN);
            if (fst24_write(File, &record, FST_SKIP) <= 0) {
                Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Could not write gridset index field (fst24_write failed)\n", __func__);
                return FALSE;
            }
        }
    }
    return TRUE;
}

//! Retrieve a gridset within the cached list
TGeoSet* GeoRef_SetGet(
    //! [in, out] Destination georeference
    TGeoRef * const RefTo,
    //! [in] Source georeference
    TGeoRef * const RefFrom,
    //! [in] Interpolation parameters
    const TGeoOptions * const Opt
) {
    //! \return Found gridset pointer, NULL on error
    if (!RefTo || !RefFrom) {
        return(NULL);
    }

    pthread_mutex_lock(&RefTo->Mutex);
   // Initialize number of geoset in cache
   if (GeoRef_SetMax < 0) {
      GeoRef_SetMax = -GeoRef_SetMax;
      char * georefSetMaxStr = getenv("GEOREF_MAXSET");
      if (georefSetMaxStr) {
         GeoRef_SetMax = atoi(georefSetMaxStr);
      }
   }

    if (!RefTo->Sets) {
        RefTo->Sets = (TGeoSet*)calloc(GeoRef_SetMax, sizeof(TGeoSet));
    }

    // Check for last set (most cases)
    if (RefTo->LastSet && RefTo->LastSet->RefFrom == RefFrom) {
        pthread_mutex_unlock(&RefTo->Mutex);
        return(RefTo->LastSet);
    }

    if (RefTo->NbSet >= GeoRef_SetMax) {
        for (int32_t i = 0; i < RefTo->NbSet; i++) {
            if (RefTo->Sets[i].RefFrom != NULL) {
                GeoRef_SetFree(&RefTo->Sets[i]);
            }
        }
        RefTo->NbSet = 0;
        RefTo->LastSet = NULL;
    }

    // Otherwise loop on sets
    int32_t i = 0;
    while (i < RefTo->NbSet) {
        if (RefTo->Sets[i].RefFrom == RefFrom) {
            RefTo->LastSet = &RefTo->Sets[i];
            pthread_mutex_unlock(&RefTo->Mutex);
            return(&RefTo->Sets[i]);
        }
        i++;
    }

    RefTo->NbSet++;
    pthread_mutex_unlock(&RefTo->Mutex);

    // If we get here, we have'nt found any sets, create a new one
    RefTo->Sets[i].RefFrom = RefFrom;
    RefTo->Sets[i].RefTo = RefTo;

    if (Opt) RefTo->Sets[i].Opt = *Opt;

    Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s: RefFrom : %p RefTo: %p\n", __func__, RefFrom, RefTo);

    return(&RefTo->Sets[i]);
}
