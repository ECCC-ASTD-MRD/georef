//! \file

#include <malloc.h>
#include <math.h>
#include <strings.h>
#include <stdlib.h>

#include <App.h>
#include "GeoRef.h"
#include "georef/Def.h"
#include "georef/Array.h"

// Sizes in bytes of the different data types
//! \todo revisit for architecture dependencies
int32_t TDef_Size[] = { 0, 1, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8 };
int32_t TDef_Mult[] = { 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 2, 3 };
int32_t TDef_DTYP[] = { 0, 0, 2, 4, 2, 4, 2, 4, 2, 4, 5, 5, -1 };


//! Re-initialize a data definition
void Def_Clear(
    //! [in, out] Data definition to be re-initialised
    TDef * const def
) {
    for(int32_t n = 0; n < def->NC; n++) {
        for(uint64_t i = 0; i < FSIZE3D(def); i++) {
            Def_Set(def, n, i, def->NoData);
        }
    }

    if (def->Buffer) {
        free(def->Buffer);
        def->Buffer = NULL;
    }
    if (def->Aux) {
        free(def->Aux);
        def->Aux = NULL;
    }
    if (def->Accum) {
        free(def->Accum);
        def->Accum = NULL;
    }
    if (def->Mask) {
        free(def->Mask);
        def->Mask = NULL;
    }
    if (def->Sub) {
        free(def->Sub);
        def->Sub = NULL;
    }

    if (def->Segments) {
        TList_Clear(def->Segments, (int(*)(void*))T3DArray_Free);
        def->Segments = NULL;
    }
}


//! Check dimensions compatibility betweed 2 data definition and adjust if needed
int32_t Def_Compat(
    //! [in, out] Data definition to be adjusted
    TDef * const dstDef,
    //! [in] Data definition
    const TDef * const srcDef
) {
    //! \return Compatible, FALSE or TRUE
    int32_t compatible = TRUE;

    if (srcDef->Idx) {
        return FALSE;
    }

    if (dstDef->Mode && dstDef->Mode != dstDef->Data[0]) {
        free(dstDef->Mode);
    }
    dstDef->Mode = NULL;
    dstDef->Dir = NULL;

    // Verifier la dimension verticale
    if (dstDef->NK != srcDef->NK || dstDef->NC != srcDef->NC) {
        if (dstDef->Idx) {
            Lib_Log(APP_LIBGEOREF, APP_WARNING, "%s: Cannot change the data size for a sub 'U' grid\n", __func__);
            return FALSE;
        }
        if (dstDef->Data[0]) free(dstDef->Data[0]);
        dstDef->Data[0] = dstDef->Data[1] = dstDef->Data[2] = dstDef->Data[3] = NULL;

        dstDef->NK = srcDef->NK;
        dstDef->NC = srcDef->NC;
        int32_t nijk = FSIZE3D(dstDef);
        if (!(dstDef->Data[0] = (char*)calloc(nijk*dstDef->NC, TDef_Size[dstDef->Type]))) {
            return FALSE;
        }
        for(int32_t n = 1; n < dstDef->NC; n++) {
            dstDef->Data[n] = &dstDef->Data[0][nijk*n*TDef_Size[dstDef->Type]];
        }
        compatible = FALSE;
    }

    // Verifier le masque
    if (srcDef->Mask) {
        if (!dstDef->Mask || compatible == 0) {
            if (!(dstDef->Mask = (char*)realloc(dstDef->Mask, FSIZE3D(dstDef)))) {
                return FALSE;
            }
        }
    } else if (dstDef->Mask) {
        free(dstDef->Mask);
        dstDef->Mask = NULL;
    }

    return compatible;
}


//! Copy a data definition
TDef * Def_Copy(
    //! [in] Data definition to copy
    const TDef * const srcDef
){
    //! \return New data definition or NULL on error
    if (!srcDef->Idx) {
        TDef * newDef = (TDef*)malloc(sizeof(TDef));
        if (srcDef && newDef) {
            newDef->Alias = srcDef->Alias;
            newDef->CellDim = srcDef->CellDim;
            newDef->NI = srcDef->NI;
            newDef->NJ = srcDef->NJ;
            newDef->NK = srcDef->NK;
            newDef->NC = srcDef->NC;
            newDef->NIJ = srcDef->NI*srcDef->NJ;
            newDef->NoData = srcDef->NoData;
            newDef->Type = srcDef->Type;
            newDef->Level = srcDef->Level;
            newDef->Idx = srcDef->Idx;
            newDef->Buffer = NULL;
            newDef->Aux = NULL;
            newDef->Segments = NULL;
            newDef->Accum = NULL;
            newDef->Mask = NULL;
            newDef->Sub = NULL;
            newDef->Pres = NULL;
            newDef->PresLS = NULL;
            newDef->Height = NULL;
            newDef->Pick = NULL;
            newDef->Poly = NULL;
            newDef->Sample = srcDef->Sample;
            newDef->SubSample = srcDef->SubSample;

            memcpy(newDef->Limits, srcDef->Limits, 6 * sizeof(int));
            newDef->CoordLimits[0][0] = srcDef->CoordLimits[0][0];
            newDef->CoordLimits[0][1] = srcDef->CoordLimits[0][1];
            newDef->CoordLimits[1][0] = srcDef->CoordLimits[1][0];
            newDef->CoordLimits[1][1] = srcDef->CoordLimits[1][1];

            int32_t nijk = FSIZE3D(srcDef);
            newDef->Data[0] = newDef->Data[1] = newDef->Data[2] = newDef->Data[3] = NULL;
            if (newDef->Alias) {
                for(int32_t n = 0; n < newDef->NC; n++) newDef->Data[n] = srcDef->Data[n];
            } else {
                if (!(newDef->Data[0] = (char*)malloc(nijk * newDef->NC * TDef_Size[newDef->Type]))) {
                    Def_Free(newDef);
                    return NULL;
                }
                memcpy(newDef->Data[0], srcDef->Data[0], nijk * newDef->NC * TDef_Size[srcDef->Type]);
                for(int32_t n = 1; n < newDef->NC; n++) {
                    newDef->Data[n] = &newDef->Data[0][nijk * n * TDef_Size[srcDef->Type]];
                }
            }
            newDef->Mode = newDef->Data[0];
            newDef->Dir = NULL;

            if (srcDef->Mask) {
                if (!(newDef->Mask = (char*)malloc(nijk))) {
                    return NULL;
                }
                memcpy(newDef->Mask, srcDef->Mask, nijk);
            }
        }
        return newDef;
    }
    return NULL;
}

//! Copy a data definition but change the internal data type
TDef *Def_CopyPromote(
    //! [in] Data definition to copy
    const TDef * const srcDef,
    //! [in] New data type
    const TDef_Type type
) {
    if (srcDef) {
        TDef * const newDef = (TDef*)malloc(sizeof(TDef));
        if (!srcDef->Idx && newDef) {
            newDef->Alias = 0;
            newDef->CellDim = srcDef->CellDim;
            newDef->NI = srcDef->NI;
            newDef->NJ = srcDef->NJ;
            newDef->NK = srcDef->NK;
            newDef->NC = srcDef->NC;
            newDef->NIJ = srcDef->NI*srcDef->NJ;
            newDef->NoData = srcDef->NoData;
            newDef->Type = type;
            newDef->Level = srcDef->Level;
            newDef->Idx = srcDef->Idx;
            newDef->Buffer = NULL;
            newDef->Aux = NULL;
            newDef->Segments = NULL;
            newDef->Accum = NULL;
            newDef->Mask = NULL;
            newDef->Sub = NULL;
            newDef->Pres = NULL;
            newDef->PresLS = NULL;
            newDef->Height = NULL;
            newDef->Pick = newDef->Poly = NULL;
            newDef->Sample = srcDef->Sample;
            newDef->SubSample = srcDef->SubSample;

            memcpy(newDef->Limits, srcDef->Limits, 6*sizeof(int));
            newDef->CoordLimits[0][0] = srcDef->CoordLimits[0][0];
            newDef->CoordLimits[0][1] = srcDef->CoordLimits[0][1];
            newDef->CoordLimits[1][0] = srcDef->CoordLimits[1][0];
            newDef->CoordLimits[1][1] = srcDef->CoordLimits[1][1];

            int32_t nijk = FSIZE3D(srcDef);
            newDef->Data[0] = newDef->Data[1] = newDef->Data[2] = newDef->Data[3] = NULL;
            if (!(newDef->Data[0] = (char*)calloc(nijk*newDef->NC, TDef_Size[newDef->Type]))) {
                Def_Free(newDef);
                return NULL;
            }
            for(int32_t n = 1; n < newDef->NC; n++) {
                newDef->Data[n] = &newDef->Data[0][nijk*n*TDef_Size[newDef->Type]];
            }

            newDef->Mode = newDef->Data[0];
            newDef->Dir = NULL;

            if (srcDef->Mask) {
                if (!(newDef->Mask = (char*)malloc(nijk))) {
                    return NULL;
                }
                memcpy(newDef->Mask, srcDef->Mask, nijk);
            }
            return newDef;
        }
    }
    return NULL;
}


//! Copy a data from a data definition to another data definition
TDef *Def_CopyData(
    //! [out] Destination data definition
    TDef * const dstDef,
    //! [in] Source data definition
    const TDef * const srcDef
) {
    //! \return New data definition
    TDef *def = dstDef;
    if (!def) {
        def = Def_Copy(srcDef);
    }

    if (!srcDef) {
        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid data source\n", __func__);
        return NULL;
    }
    //! \todo Check Def COMPAT

    int32_t alias = FALSE;
    //! \todo define Alias
    for(int32_t n = 0; n < srcDef->NC; n++) {
        if (srcDef->Data[n]) {
            if (alias) {
                def->Data[n] = srcDef->Data[n];
            } else {
                memcpy(def->Data[n], srcDef->Data[n], FSIZE3D(srcDef)*TDef_Size[srcDef->Type]);
            }
        }
    }
    if (alias) def->Mode = srcDef->Mode;

    if (srcDef->Mask) {
        def->Mask = (char*)malloc(srcDef->NIJ);
        memcpy(def->Mask, srcDef->Mask, srcDef->NIJ);
    }
    return def;
}

//! Free a data from a data definition
void Def_Free(
    //! [in] Data definition to free
    TDef * const Def
) {
    if (Def) {
        if (!Def->Alias && !Def->Idx) {
            if (Def->Mode && Def->Mode != Def->Data[0]) free(Def->Mode);
            if (Def->Data[0]) free(Def->Data[0]);
            if (Def->Mask) free(Def->Mask);
        }

        if (Def->Buffer) free(Def->Buffer);
        if (Def->Aux) free(Def->Aux);
        if (Def->Accum) free(Def->Accum);
        if (Def->Sub) free(Def->Sub);
        if (Def->Pres > (float*)0x1) free(Def->Pres);
        if (Def->PresLS > (float*)0x1) free(Def->PresLS);
        if (Def->Height > (float*)0x1) free(Def->Height);
#ifdef HAVE_GDAL
        if (Def->Poly) OGR_G_DestroyGeometry(Def->Poly);
        //Freed by Def->Poly      if (Def->Pick)       OGR_G_DestroyGeometry(Def->Pick);
#endif
        if (Def->Segments) TList_Clear(Def->Segments, (int(*)(void*))T3DArray_Free);

        free(Def);
    }
}


//! Create a new data definition
TDef *Def_New(
    //! [in] X dimension
    const int32_t NI,
    //! [in] Y dimension
    const int32_t NJ,
    //! [in] Z dimension
    const int32_t NK,
    //! [in] Number of components (1-4)
    const int32_t Dim,
    //! [in] Data type
    const TDef_Type Type,
    //! [in] Is this an alias to another data definition
    const int32_t Alias
) {
    //! \return New data definition, NULL on error
    TDef * def = (TDef*)malloc(sizeof(TDef));

    if (!def) return NULL;

    def->NI = NI;
    def->NJ = NJ;
    def->NK = NK;
    def->NIJ = NI*NJ;
    def->NC = abs(Dim);
    def->Alias = Alias;
    def->CellDim = 2;
    def->NoData = nan("NaN");
    def->Level = 0;
    def->Idx = 0;

    def->Limits[0][0] = 0;
    def->Limits[1][0] = 0;
    def->Limits[2][0] = 0;
    def->Limits[0][1] = NI-1;
    def->Limits[1][1] = NJ-1;
    def->Limits[2][1] = NK-1;

    def->CoordLimits[0][0] = -180;
    def->CoordLimits[0][1] = 180;
    def->CoordLimits[1][0] = -90;
    def->CoordLimits[1][1] = 90;

    def->SubSample = def->Sample = 1;

    def->Type = Type;
    def->Buffer = NULL;
    def->Aux = NULL;
    def->Segments = NULL;
    def->Accum = NULL;
    def->Mask = NULL;
    def->Sub = NULL;
    def->Pres = NULL;
    def->PresLS = NULL;
    def->Height = NULL;
    def->Pick = NULL;
    def->Poly = NULL;

    // Allocate a single buffer and set component pointer to right address
    int32_t nijk = FSIZE3D(def);
    def->Data[0] = def->Data[1] = def->Data[2] = def->Data[3] = def->Dir = NULL;

    if (!def->Alias) {
        if (!(def->Data[0] = (char*)calloc(nijk*def->NC, TDef_Size[Type]))) {
            Def_Free(def);
            return NULL;
        }

        for(int32_t n = 1; n < def->NC; n++) {
            def->Data[n] = &def->Data[0][nijk*n*TDef_Size[def->Type]];
        }
    }
    def->Mode = def->Data[0];

    return def;
}

//! Create a new data definition and initialize the internal data pointers
TDef *Def_Create(
    //! [in] X dimension
    const int32_t NI,
    //! [in] Y dimension
    const int32_t NJ,
    //! [in] Z dimension
    const int32_t NK,
    //! [in] Number of components (1-4)
    const TDef_Type Type,
    //! [in] 1st component
    char * const Comp0,
    //! [in] 2nd component (if vectorial field)
    char * const Comp1,
    //! [in] Mask, if available
    char * const Mask
) {
    //! \return New data definition, NULL on error
    int32_t dim = (Comp0 != NULL) + (Comp1 != NULL);
    TDef * def = Def_New(NI, NJ, NK, dim, Type, TRUE);
    if (!def) return NULL;

    def->Data[0] = Comp0;
    def->Data[1] = Comp1;
    def->Mask = Mask;
    return def;
}

/*
static inline float* Def_ToFloat(TDef *Def, char* Buffer) {

   int32_t n, sz;
   float *buf = NULL;

   sz = FSIZE3D(Def);
   *buf = (float*)malloc(Size*sizeof(float));

   switch(Def->Type) {
      case TD_Unknown : break;
      case TD_Binary: break;
      case TD_UByte:  for(n = 0; n < sz; n++) buf[n] = ((unsigned char*)Buffer)[n]; break;
      case TD_Byte:   for(n = 0; n < sz; n++) buf[n] = ((char*)Buffer)[n]; break;
      case TD_UInt16: for(n = 0; n < sz; n++) buf[n] = ((unsigned short*)Buffer)[n]; break;
      case TD_Int16:  for(n = 0; n < sz; n++) buf[n] = ((short*)Buffer)[n]; break;
      case TD_UInt32: for(n = 0; n < sz; n++) buf[n] = ((unsigned int*)Buffer)[n]; break;
      case TD_Int32:  for(n = 0; n < sz; n++) buf[n] = ((int*)Buffer)[n]; break;
      case TD_UInt64: for(n = 0; n < sz; n++) buf[n] = ((unsigned long long*)Buffer)[n]; break;
      case TD_Int64:  for(n = 0; n < sz; n++) buf[n] = ((long long*)Buffer)[n]; break;
      case TD_Float32: buf = (float*)Buffer; break;
      case TD_Float64: for(n = 0; n < sz; n++) buf[n] = ((double*)Buffer)[n]; break;
   }

   return buf;
}
static inline float* Def_FromFloat(TDef *Def) {

   int32_t n;

   sz = FSIZE3D(Def);
   switch(Def->Type) {
      case TD_Unknown : break;
      case TD_Binary: break;
      case TD_UByte:  for(n = 0; n < sz; n++) ((unsigned char*)Buffer)[n] = buf[n]; break;
      case TD_Byte:   for(n = 0; n < sz; n++) ((char*)Buffer)[n] = buf[n]; break;
      case TD_UInt16: for(n = 0; n < sz; n++) ((unsigned short*)Buffer)[n] = buf[n]; break;
      case TD_Int16:  for(n = 0; n < sz; n++) ((short*)Buffer)[n] = buf[n]; break;
      case TD_UInt32: for(n = 0; n < sz; n++) ((unsigned int*)Buffer)[n] = buf[n]; break;
      case TD_Int32:  for(n = 0; n < sz; n++) ((int*)Buffer)[n] = buf[n]; break;
      case TD_UInt64: for(n = 0; n < sz; n++) ((unsigned long long*)Buffer)[n] = buf[n]; break;
      case TD_Int64:  for(n = 0; n < sz; n++) ((long long*)Buffer)[n] = buf[n]; break;
      case TD_Float32: (float*)Buffer = buf; break;
      case TD_Float64: for(n = 0; n < sz; n++) ((double*)Buffer)[n] = buf[n]; break;
   }

   return buf;
}
*/

//! Resize the data field of a data definition
TDef * Def_Resize(
    //! [in, out] Data definition to be resized
    TDef * const Def,
    //! [in] X dimensionn
    const int32_t NI,
    //! [in] Y dimension
    const int32_t NJ,
    //! [in] Z dimension
    const int32_t NK
) {
    //! \return New data definition
    if (!Def) return NULL;

    if (!Def->Alias && !Def->Idx && (Def->NI != NI || Def->NJ != NJ || Def->NK != NK)) {
        Def->NI = NI;
        Def->NJ = NJ;
        Def->NK = NK;
        Def->NIJ = NI*NJ;

        Def->Limits[0][0] = 0;
        Def->Limits[1][0] = 0;
        Def->Limits[2][0] = 0;
        Def->Limits[0][1] = NI-1;
        Def->Limits[1][1] = NJ-1;
        Def->Limits[2][1] = NK-1;

        Def->CoordLimits[0][0] = -180;
        Def->CoordLimits[0][1] = 180;
        Def->CoordLimits[1][0] = -90;
        Def->CoordLimits[1][1] = 90;

        Def->SubSample = Def->Sample = 1;

        if (Def->Mode && Def->Mode != Def->Data[0]) {
            free(Def->Mode);
            Def->Mode = NULL;
        }

        // Allocate a single buffer and set component pointer to right address
        Def->Data[1] = Def->Data[2] = Def->Data[3] = NULL;
        int32_t nijk = FSIZE3D(Def);
        if (!(Def->Data[0] = (char*)realloc(Def->Data[0], nijk * Def->NC * TDef_Size[Def->Type]))) {
            Def_Free(Def);
            return NULL;
        }
        for(int32_t n = 1; n < Def->NC; n++) {
            Def->Data[n] = &Def->Data[0][nijk * n * TDef_Size[Def->Type]];
        }
        Def->Mode = Def->Data[0];
        Def->Dir = NULL;

        if (Def->Buffer) {
            free(Def->Buffer);
            Def->Buffer = NULL;
        }
        if (Def->Aux) {
            free(Def->Aux);
            Def->Aux = NULL;
        }
        if (Def->Accum) {
            free(Def->Accum);
            Def->Accum = NULL;
        }
        if (Def->Mask) {
            free(Def->Mask);
            Def->Mask = NULL;
        }
        if (Def->Sub) {
            free(Def->Sub);
            Def->Sub = NULL;
        }
        if (Def->Pres > (float*)0x1) {
            free(Def->Pres);
            Def->Pres = NULL;
        }
        if (Def->PresLS > (float*)0x1) {
            free(Def->PresLS);
            Def->PresLS = NULL;
        }
        if (Def->Height > (float*)0x1) {
            free(Def->Height);
            Def->Height = NULL;
        }
    }
    return Def;
}


//! Copy data from a data definition to another one
int32_t Def_Paste(
    //! [in, out] Destination data definition
    TDef * const dstDef,
    //! [in] Source data definition
    const TDef * const srcDef,
    //! [in] X start coordinate
    int32_t X0,
    //! [in] Y start Coordinate
    int32_t Y0
) {
    //! \return FALSE if out of grid, TRUE otherwise

    // Check limits
    int32_t x0 = X0 < 0 ? -X0 : 0;
    int32_t y0 = Y0 < 0 ? -Y0 : 0;

    int32_t x1 = srcDef->NI + X0 > dstDef->NI ? dstDef->NI : srcDef->NI;
    int32_t y1 = srcDef->NJ + Y0 > dstDef->NJ ? dstDef->NJ : srcDef->NJ;

    // If paste is out of destination
    if (x0 > dstDef->NI || x1 < 0 || y0 > dstDef->NJ || y1 < 0) {
        return FALSE;
    }

    // Maximum number of band to paste
    int32_t nc = fmin(dstDef->NC, srcDef->NC);

    int32_t dy = Y0;
    double a = 1.0, srcVal, dstVal;
    for (int32_t y = y0; y < y1; y++) {
        int32_t dx = X0;
        for (int32_t x = x0; x < x1; x++) {
            unsigned long idxs = FIDX2D(srcDef, x, y);
            unsigned long idxd = FIDX2D(dstDef, dx, dy);

            if (srcDef->NC == 4) {
                Def_Get(srcDef, 3, idxs, a);
            }

            for(int32_t c = 0; c < nc; c++) {
                Def_Get(srcDef, c, idxs, srcVal);

                if (srcDef->NC == 4) {
                    Def_Get(dstDef, c, idxd, dstVal);
                    srcVal = srcVal * a + dstVal * (1.0 - a);
                }
                Def_Set(dstDef, c, idxd, srcVal);
            }
            dx++;
        }
        dy++;
    }
    return TRUE;
}


//! Get a value within a data definition, possibly interpolated
int32_t Def_GetValue(
    //! [in] GeoRef pointer
    const TGeoRef * const Ref,
    //! [in] Data definition
    const TDef * const Def,
    //! [in] Interpolation options
    const TGeoOptions * const Opt,
    //! [in] Component (0-4)
    const int32_t C,
    //! [in] X coordinate
    double X,
    //! [in] Y Coordinate
    double Y,
    //! [in] Z Coordinate
    double Z,
    //! [out] Value to Assign (vector length if vector value)
    double * const Length,
    //! [out] Rotation angle, if vector value
    double * const ThetaXY
) {
    //! \return FALSE if out of grid, TRUE otherwise
    const TGeoOptions * const opt = Opt ? Opt : &GeoRef_Options;

    *Length = Def->NoData;
    double d = Ref->GRTYP[0] == 'W' ? 1.0 : 0.5;

    if (C >= Def->NC) return FALSE;

    float valf, valdf;
    void *p0, *p1;
    // Si on est a l'interieur de la grille
    if ((Ref->GRTYP[0] == 'M' && X >= 0 && X < Ref->NIdx) ||
        (X >= (Ref->X0-d) && Y >= (Ref->Y0-d) && Z >= 0 && X <= (Ref->X1+d) && Y <= (Ref->Y1+d) && Z <= Def->NK-1)) {

        // Index memoire du niveau desire
        int32_t mem = Def->NIJ * (int)Z;

        double x = X;
        double y = Y;
        x -= Ref->X0;
        y -= Ref->Y0;
        int32_t ix = lrint(x);
        int32_t iy = lrint(y);
        uint32_t idx = (Ref->GRTYP[0] == 'M') ? ix : iy * Def->NI + ix;

        // Check for mask
        if (Def->Mask && !Def->Mask[idx]) {
            if (!Def->Mask[idx]) {
                return FALSE;
            }
        }

        if (Def->Data[1] && !C) {
            Def_Pointer(Def, 0, mem, p0);
            Def_Pointer(Def, 1, mem, p1);
            GeoRef_XYWDVal(Ref, opt, &valf, &valdf, p0, p1, &X, &Y, 1);
            *Length = valf;

            // If it's 3D, use the mode for speed since GeoRef_XYWDVal only uses 2D
            if (Def->Data[2]) {
                GeoRef_XYVal(Ref, opt, &valf, (float*)&Def->Mode[mem], &X, &Y, 1);
                *Length = valf;
            }
            if (ThetaXY) *ThetaXY = valdf;
        } else {
            if (Def->Type <= TD_Int64 || (ix == x && iy == y)) {
                Def_Get(Def, C, mem + idx, valf);
            } else {
                Def_Pointer(Def, C, mem, p0);
                GeoRef_XYVal(Ref, opt, &valf, p0, &X, &Y, 1);
            }
            *Length = valf;
        }
        return TRUE;
    }

    return FALSE;
}


//! Set a value within a data definition
static inline void Def_SetValue(
    //! [in] Data definition
    TDef *Def,
    //! [in] Options
    TGeoOptions *Opt,
    //! [in] X coordinate
    int32_t X,
    //! [in] Y Coordinate
    int32_t Y,
    //! [in] Z Coordinate
    int32_t Z,
    //! [in] Value to assign
    double Value
) {
    //! \return FALSE on error, TRUE otherwise

    if (FIN2D(Def, X, Y)) {
        unsigned long idx = Z ? FIDX3D(Def, X, Y, Z) : FIDX2D(Def, X, Y);

        if (Opt->Combine == CB_REPLACE) {
            Def_Set(Def, 0, idx, Value);
        } else {
            double val;
            Def_Get(Def, 0, idx, val);
            if (!DEFVALID(Def, val)) val = 0.0;

            switch(Opt->Combine) {
                case CB_MIN    : if (Value < val) Def_Set(Def, 0, idx, Value); break;
                case CB_MAX    : if (Value > val) Def_Set(Def, 0, idx, Value); break;
                case CB_AVERAGE: Def->Accum[idx]++;
                case CB_SUM    : Value += val;  Def_Set(Def, 0, idx, Value); break;
                case CB_REPLACE: break;
            }
        }
    }
}

/*----------------------------------------------------------------------------
 * @brief  Rasterize vectorial data within aband or field
 * @date   June 2004
    //! [in]     ToRef      Destination GeoRef pointer
    //! [in]     ToDef      Destination data definition
    //! [in]     Geom       Geometry to initialize to rasterize
    //! [in]     Value      Value to assign
    //! [in]     Comb       Result combination method (CB_REPLACE, CB_MIN, CB_MAX, CB_SUM, CB_AVERAGE)
    //! [in]     Seg        Number of segmentation on the cell sides
 *
 *    @return        FALSE on error, TRUE otherwise
*/
int32_t GeoRef_Rasterize(TGeoRef *ToRef, TDef *ToDef, TGeoOptions *Opt, OGRGeometryH Geom, double Value) {

#ifdef HAVE_GDAL
   int32_t i, j, g, ind1, ind2;
   int32_t x, y, miny, maxy, minx, maxx;
   int32_t ints, n, ns, np;
   int32_t *polyInts;
   double  dminy, dmaxy, dx1, dy1, dx2, dy2, dy;
   double  intersect, tmpd;
   int32_t horizontal_x1, horizontal_x2;
   int32_t dnx, dny, x0, x1, y0, y1, fr, sx, sy;

   if (!Opt) Opt=&GeoRef_Options;

   OGRGeometryH geom;

   if (!ToRef || !ToDef) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid destination\n", __func__);
      return 0;
   }
   if (!Geom) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid source\n", __func__);
      return 0;
   }

   n = 0;
   dminy = 1e300;
   dmaxy = -1e300;

   for (i = 0; i < OGR_G_GetGeometryCount(Geom); i++) {
      geom = OGR_G_GetGeometryRef(Geom, i);
      if (EQUAL(OGR_G_GetGeometryName(geom), "LINEARRING")) {
         n += ns = OGR_G_GetPointCount(geom);
         for (j = 0; j < ns; j++) {
            dy = OGR_G_GetY(geom, j);
            if (dy < dminy)
               dminy = dy;
            if (dy > dmaxy)
               dmaxy = dy;
         }
      } else {
         GeoRef_Rasterize(ToRef, ToDef, Opt, geom, Value);
      }
   }
   if (!(n += OGR_G_GetPointCount(Geom)))
      return 1;

   switch (OGR_G_GetDimension(Geom)) {
      case 0: // Point type
         for (i = 0; i < n; i++) {
            dx1 = OGR_G_GetX(Geom, i);
            dy1 = OGR_G_GetY(Geom, i);
            x = lrint(dx1);
            y = lrint(dy1);
            Def_SetValue(ToDef, Opt, x, y, 0, Value);
         }
         break;

      case 1: // Line type
         for (i = 1; i < n; i++) {
            ind2 = i;
            ind1 = i == 0 ? n-1 : i-1;

            dx1 = OGR_G_GetX(Geom, ind1);
            dy1 = OGR_G_GetY(Geom, ind1);

            dx2 = OGR_G_GetX(Geom, ind2);
            dy2 = OGR_G_GetY(Geom, ind2);

            x0 = lrint(dx1); y0 = lrint(dy1);
            x1 = lrint(dx2); y1 = lrint(dy2);
            dny = y1-y0;
            dnx = x1-x0;
            if (dny < 0) {
               dny = -dny;
               sy = -1;
            } else {
               sy = 1;
            }
            if (dnx < 0) {
               dnx = -dnx;
               sx = -1;
            } else {
               sx = 1;
            }
            dny <<= 1;
            dnx <<= 1;

            if (FIN2D(ToDef, x0, y0))
               Def_Set(ToDef, 0, FIDX2D(ToDef, x0, y0), Value);
            if (dnx > dny) {
               fr = dny-(dnx >> 1);
               while(x0 != x1) {
                  if (fr >= 0) {
                     y0 += sy;
                     fr -= dnx;
                  }
                  x0 += sx;
                  fr += dny;
                  Def_SetValue(ToDef, Opt, x0, y0, 0, Value);
               }
            } else {
               fr = dnx-(dny >> 1);
               while(y0 != y1) {
                  if (fr >= 0) {
                     x0 += sx;
                     fr -= dny;
                  }
                  y0 += sy;
                  fr += dnx;
                  Def_SetValue(ToDef, Opt, x0, y0, 0, Value);
               }
            }
         }
         break;

      case 2: // Polygon type
         miny = (int)(dminy < 0 ? 0 : dminy);
         maxy = (int)(dmaxy >= ToDef->NJ ? ToDef->NJ-1 : dmaxy);
         minx = 0;
         maxx = ToDef->NI-1;

         polyInts = (int*)malloc(sizeof(int)*n);

         // Fix in 1.3: count a vertex only once
         for (y = miny; y <= maxy; y++) {
            dy = y; // center height of line
            ints = 0 ;

            // Initialize polyInts, otherwise it can sometimes causes a segfault
            for (i = 0; i < n; i++) {
               polyInts[i] = -1;
            }

            ns = OGR_G_GetGeometryCount(Geom);
            for (g = 0; g < (ns == 0 ? 1 : ns); g++) {
               if (ns) {
                  geom = OGR_G_GetGeometryRef(Geom, g);
               } else {
                  geom = Geom;
               }
               np = OGR_G_GetPointCount(geom);

               for (i = 0; i < np; i++) {
                  ind2 = i;
                  ind1 = i == 0 ? np-1 : i-1;

                  dx1 = OGR_G_GetX(geom, ind1);
                  dy1 = OGR_G_GetY(geom, ind1);

                  dx2 = OGR_G_GetX(geom, ind2);
                  dy2 = OGR_G_GetY(geom, ind2);

                  if ((dy1 < dy && dy2 < dy) || (dy1 > dy && dy2 > dy))
                     continue;

                  if (dy1 < dy2) {
                  } else if (dy1 > dy2) {
                     tmpd = dy2;
                     dy2 = dy1;
                     dy1 = tmpd;
                     tmpd = dx2;
                     dx2 = dx1;
                     dx1 = tmpd;
                  } else { // if (fabs(dy1-dy2)< 1.e-6)
                           /*AE: DO NOT skip bottom horizontal segments
                           -Fill them separately-
                           They are not taken into account twice.*/
                     if (dx1 > dx2) {
                        horizontal_x1 = lrint(dx2);
                        horizontal_x2 = lrint(dx1);
                        if ((horizontal_x1 > maxx) || (horizontal_x2 < minx))
                           continue;

                        // fill the horizontal segment (separately from the rest)
                        for(x = horizontal_x1; x < horizontal_x2; x++)
                           Def_SetValue(ToDef, Opt, x, y, 0, Value);
                        continue;
                     } else {
                        // skip top horizontal segments (they are already filled in the regular loop)
                        continue;
                     }
                  }

                  if ((dy < dy2) && (dy >= dy1)) {
                     intersect = (dy-dy1)*(dx2-dx1)/(dy2-dy1)+dx1;
                     polyInts[ints++] = lrint(intersect);
                  }
               }
            }
            qsort(polyInts, ints, sizeof(int), QSort_Int);

            for (i = 0; i < ints; i += 2) {
               if (polyInts[i] <= maxx && polyInts[i+1] >= minx) {
                  for(x = polyInts[i]; x <= polyInts[i+1]; x++)
                     Def_SetValue(ToDef, Opt, x, y, 0, Value);
               }
            }
         }
         free(polyInts);
         break;
   }
   return 1;
#else
   Lib_Log(APP_LIBGEOREF, APP_ERROR, "Function %s is not available, needs to be built with GDAL\n", __func__);
   return 0;
#endif
}

/*----------------------------------------------------------------------------
 * @brief  Project a grid cell within  another referencial
 * @date   March 2006
    //! [in]     Geom       Geometry to initialize with projected coordinates
    //! [in]     ToRef      Destination GeoRef pointer
    //! [in]     FromRef    Source GeoRef pointer
    //! [in]     I          X gridpoint coordinate (Fortran index)
    //! [in]     J          Y gridpoint coordinate (Fortran index)
    //! [in]     Opt        GeoRef options
 *
 *    @return        Number of points in geometry (negative if wrapped around)
*/
int32_t GeoRef_Cell2OGR(OGRGeometryH Geom, TGeoRef *ToRef, TGeoRef *FromRef, int32_t I, int32_t J, TGeoOptions *Opt) {

#ifdef HAVE_GDAL
   double   n, dn, df;
   double   x0, x1, y0, y1, t0, t1, x, y, la, lo;
   int32_t  pt = 0, dx,dy;
   TGeoRef *fromref=FromRef;
   TGeoSet *gset;

   if (!Geom || !ToRef || !FromRef|| !Opt) {
      return 0;
   }

   gset = GeoRef_SetGet(ToRef, FromRef, Opt);

   dn = Opt->Segment?1.0/Opt->Segment:1.0;
   df = dn * 0.5;

   y0 = x0 = 1e32;
   y1 = x1 = -1e32;

   // In case of U grid, we have to select the right subgrid
   if (FromRef->NbSub > 0) {
      fromref=FromRef->Subs[0];
      if (J > FromRef->Subs[0]->NY-1) {
         J-=FromRef->Subs[0]->NY;
         fromref=FromRef->Subs[1];
      
         x = I; y = J;
         GeoRef_XY2LL(fromref, &la, &lo, &x, &y, 1, TRUE);
         if (GeoRef_LL2XY(ToRef, &x, &y, &la, &lo, 1, FALSE)) {
            dx=ROUND(x);
            dy=ROUND(y);
            // Check if this point is from yin or yang
            if (!gset->yin_maskout[dy*ToRef->NX+dx])
               return(0);
         }
      }
   }

   // Top Left
   for(n = -0.5; n < (0.5+df); n += dn) {
      x = I-0.5; y = J+n;
      GeoRef_XY2LL(fromref, &la, &lo, &x, &y, 1, TRUE);
      GeoRef_LL2XY(ToRef, &x, &y, &la, &lo, 1, TRUE);
      x0 = fmin(x0, x);
      x1 = fmax(x1, x);
      y0 = fmin(y0, y);
      y1 = fmax(y1, y);
      OGR_G_SetPoint_2D(Geom, pt++, x, y);
   }

   // Top right
   for(n = -0.5; n < (0.5+df); n += dn) {
      x = I+n; y = J+0.5;
      GeoRef_XY2LL(fromref, &la, &lo, &x, &y, 1, TRUE);
      GeoRef_LL2XY(ToRef, &x, &y, &la, &lo, 1, TRUE);
      x0 = fmin(x0, x);
      x1 = fmax(x1, x);
      y0 = fmin(y0, y);
      y1 = fmax(y1, y);
      OGR_G_SetPoint_2D(Geom, pt++, x, y);
   }

   // Right bottom
   for(n = 0.5; n > -(0.5+df); n -= dn) {
      x = I+0.5; y = J+n;
      GeoRef_XY2LL(fromref, &la, &lo, &x, &y, 1, TRUE);
      GeoRef_LL2XY(ToRef, &x, &y, &la, &lo, 1, TRUE);
      x0 = fmin(x0, x);
      x1 = fmax(x1, x);
      y0 = fmin(y0, y);
      y1 = fmax(y1, y);
      OGR_G_SetPoint_2D(Geom, pt++, x, y);
   }

   // Bottom Left
   for(n = 0.5; n > -(0.5+df); n -= dn) {
      x = I+n; y = J-0.5;
      GeoRef_XY2LL(fromref, &la, &lo, &x, &y, 1, TRUE);
      GeoRef_LL2XY(ToRef, &x, &y, &la, &lo, 1, TRUE);
      x0 = fmin(x0, x);
      x1 = fmax(x1, x);
      y0 = fmin(y0, y);
      y1 = fmax(y1, y);
      OGR_G_SetPoint_2D(Geom, pt++, x, y);
   }

   // Close the polygon
   x = I-0.5; y = J-0.5;
   GeoRef_XY2LL(fromref, &la, &lo, &x, &y, 1, TRUE);
   GeoRef_LL2XY(ToRef, &x, &y, &la, &lo, 1, TRUE);
   OGR_G_SetPoint_2D(Geom, pt++, x, y);

   // If the cell is outside the destination limits         
   t0=ToRef->X0-0.5;
   t1=ToRef->X1+0.5;
   if ((x0 < t0 && x1 < t0) || (x0 > t1 && x1 > t1)) {
      return 0;
   }

   t0=ToRef->Y0-0.5;
   t1=ToRef->Y1+0.5;
   if ((y0 < t0 && y1 < t0) || (y0 > t1 && y1 > t1)) {
      return 0;
   }

   // If the size is larger than half the destination, it has to be a wrap around
   if ((x1-x0) > ((ToRef->X1-ToRef->X0) >> 1)) {
      return -pt;
   }

   return pt;
#else
   Lib_Log(APP_LIBGEOREF, APP_ERROR, "Function %s is not available, needs to be built with GDAL\n", __func__);
   return 0;
#endif
}

/*----------------------------------------------------------------------------
 * @brief  Interpolate within a gridpoint defined quad area
 * @date   November 2004
    //! [in]     Ref        GeoRef pointer
    //! [in]     Def        Data definition
    //! [in]     Geom       Geometry tested
    //! [in]     Area       Geometry area
    //! [in]     Value      Value o be assigned
    //! [in]     Mode       Rasterization mode
    //! [in]     Type       Rasterization type
    //! [in]     X0         Lower left corner X
    //! [in]     Y0         Lower left corner Y
    //! [in]     X1         Higher right corner X
    //! [in]     Y1         Higher right corner X
    //! [in]     Z          Level
    //! [in]     Index      Local index pointer to be filled
 *
 *    @return        number of points involved
*/
static int32_t GeoRef_InterpQuad(TGeoRef *Ref, TDef *Def, TGeoOptions *Opt, OGRGeometryH Poly, OGRGeometryH Geom, char Mode, char Type, double Area, double Value, int32_t X0, int32_t Y0, int32_t X1, int32_t Y1, int32_t Z, float **Index) {

#ifdef HAVE_GDAL
   double        dx, dy, dp = 1.0, val = 0.0;
   int32_t       x, y, n = 0, idx, na;
   OGRGeometryH  inter = NULL, pick;
   OGREnvelope   envg, envp;

   // Setup the intersecting area
   pick = OGR_G_GetGeometryRef(Poly, 0);
   dx = (double)X0-0.5;
   dy = (double)Y0-0.5;
   OGR_G_SetPoint(pick, 0, dx, dy, 0);
   OGR_G_SetPoint(pick, 4, dx, dy, 0);
   dx = (double)X0-0.5;
   dy = (double)Y1+0.5;
   OGR_G_SetPoint(pick, 1, dx, dy, 0);
   dx = (double)X1+0.5;
   dy = (double)Y1+0.5;
   OGR_G_SetPoint(pick, 2, dx, dy, 0);
   dx = (double)X1+0.5;
   dy = (double)Y0-0.5;
   OGR_G_SetPoint(pick, 3, dx, dy, 0);

   OGR_G_GetEnvelope(pick, &envp);
   OGR_G_GetEnvelope(Geom, &envg);

   na = (Mode == 'C' || Mode == 'N' || Mode == 'A');

   // Test for intersection
   if ((Area > 1e-32 || !na) && OGM_Intersect(Geom, Poly, &envg, &envp)) {
//   if ((Area > 0.0 || !na) && OGR_G_Intersects(Geom, Poly)) {

      // If this is a single pixel
      if (X0 == X1 && Y0 == Y1) {

         // If we are computing areas
         if (Area > 1e-32) {
            switch(Type) {
               case 'A':  // Area mode
#ifdef HAVE_GPC
                  inter = OGM_GPCOnOGR(GPC_INT, Geom, Poly);
#else
                  inter = OGR_G_Intersection(Geom, Poly);
#endif
                  if (Mode == 'C' || Mode == 'N') {
                     dp = OGR_G_Area(inter)/Area;
                  } else if (Mode == 'A') {
                     dp = OGR_G_Area(inter);
                  }
                  break;

               case 'L':  // Length mode
                  inter=OGM_Clip(Geom, Poly);
                  if (Mode=='C' || Mode=='N') {
                     dp=OGM_Length(inter)/Area;
                  } else if (Mode=='A') {
                     dp=OGM_Length(inter);
                  }
                  break;

               case 'P':  // Point mode
                  dp = 1.0;
                  break;
            }
            val = Value*dp;
            OGR_G_DestroyGeometry(inter);
         } else {
            val = Value;
         }

         // Are we within
         if (Mode != 'W' || OGM_Within(Poly, Geom, &envp, &envg)) {
            //Create thread safe region.
            #pragma omp critical
            {
               Def_SetValue(Def, Opt, X0, Y0, Z, val);

               if (Mode == 'N' && Def->Buffer) {
                  idx = FIDX3D(Def, X0, Y0, Z);
                  Def->Buffer[idx] += dp;
               }
            }

            if (*Index) {
               if (*((*Index)+1) == REF_INDEX_NIL) {
                  Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Overflow in index array (%i,%i)\n",__func__,X0,Y0);
                  return(0);
               }
               *((*Index)++) = X0;
               *((*Index)++) = Y0;
               *((*Index)++) = dp;
            }
            n = 1;
         }
      } else {
         // Otherwise, refine the search by quad dividing
         x = (X1-X0) >> 1;
         y = (Y1-Y0) >> 1;

         // If within 1 bloc, parse them all
         if (x == 0 || y == 0) {
            for (x = X0; x <= X1; x++) {
               for (y = Y0; y <= Y1; y++) {
                  n += GeoRef_InterpQuad(Ref, Def, Opt, Poly, Geom, Mode, Type, Area, Value, x, y, x, y, Z, Index);
               }
            }
         } else {
            n += GeoRef_InterpQuad(Ref, Def, Opt, Poly, Geom, Mode, Type, Area, Value, X0, Y0, X0+x, Y0+y, Z, Index);
            n += GeoRef_InterpQuad(Ref, Def, Opt, Poly, Geom, Mode, Type, Area, Value, X0+x+1, Y0, X1, Y0+y, Z, Index);
            n += GeoRef_InterpQuad(Ref, Def, Opt, Poly, Geom, Mode, Type, Area, Value, X0, Y0+y+1, X0+x, Y1, Z, Index);
            n += GeoRef_InterpQuad(Ref, Def, Opt, Poly, Geom, Mode, Type, Area, Value, X0+x+1, Y0+y+1, X1, Y1, Z, Index);
         }
      }
   }

   return n;
#else
   Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Function not available, needs to be built with GDAL\n", __func__);
   return 0;
#endif
}

/*----------------------------------------------------------------------------
 * @brief  Interpolate/import vectorial data into a raster field
 * @remark The index format is a stream of float32 containing for each OGR feature, the list of gridcells and interpolation
 *      factor separated by -1,
 *      [ogr feature index] [feature value area] [grid cell x] [grid cell y] [factor] [grid cell x] [grid cell y] [factor] ... -1.0 ...
 * @date   November 2004
    //! [in]     Ref        Destination GeoRef pointer
    //! [in]     Def        Destination data definition
    //! [in]     LayerRef   Source layer GeoRef pointer
    //! [in]     Layer      Source layer data
    //! [in]     Field      Layer fiel to use
    //! [in]     Value      Value o be assigned
 *
 *    @return        FALSE on error, TRUE otherwise
*/
int32_t GeoRef_InterpOGR(TGeoRef *ToRef, TDef *ToDef, TGeoRef *LayerRef, OGR_Layer *Layer, TGeoOptions *Opt, char *Field, double Value) {

#ifdef HAVE_GDAL
   TGeoSet *gset;
   int64_t f, n = 0, idx2;
   uint64_t nt = 0;
   double   value, val, area, dp, x0, y0;
   int32_t  fld = -1, pi, pj, error = 0, isize = 0;
   char     mode, type, *c;
   float   *ip = NULL, *lp = NULL, **index = NULL;
   TCoord   co;
   Vect3d   vr;

   OGRSpatialReferenceH          srs = NULL;
   OGRCoordinateTransformationH  tr = NULL;
   OGRGeometryH                  geom = NULL, utmgeom = NULL, hgeom, pick = NULL, poly = NULL;
   OGREnvelope                   env;

   if (!Opt) Opt=&GeoRef_Options;

   if (!ToRef || !ToDef) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid destination\n", __func__);
      return FALSE;
   }
   if (!LayerRef || !Layer) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid source\n", __func__);
      return FALSE;
   }

   if (!Layer->NFeature) {
      return TRUE;
   }

   // Recuperer la valeur a utiliser dans l'interpolation
   if (Field) {
      if (strcmp(Field, "FEATURE_AREA") == 0) {
         fld = -2;
      } else if (strcmp(Field, "FEATURE_AREA_METER") == 0) {
         fld = -3;
      } else if (strcmp(Field, "FEATURE_LENGTH") == 0) {
         fld = -4;
      } else if (strcmp(Field, "FEATURE_LENGTH_METER") == 0) {
         fld = -5;
      } else if (strcmp(Field, "FEATURE_ID") == 0) {
         fld = -6;
      } else if (strcmp(Field, "ZCOORD_MIN") == 0) {
         fld = -7;
      } else if (strcmp(Field, "ZCOORD_MAX") == 0) {
         fld = -8;
      } else if (strcmp(Field, "ZCOORD_AVG") == 0) {
         fld = -9;
      } else {
         fld = OGR_FD_GetFieldIndex(Layer->Def, Field);
         if (fld == -1) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid layer field\n", __func__);
            return FALSE;
         }
      }
   }

   // In case of average, we need an accumulator
   if (Opt->Combine == CB_AVERAGE) {
      if (!ToDef->Accum) {
         ToDef->Accum = malloc(FSIZE2D(ToDef)*sizeof(int));
         if (!ToDef->Accum) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable to allocate accumulation buffer\n", __func__);
            return FALSE;
         }
      }
      memset(ToDef->Accum, 0x0, FSIZE2D(ToDef)*sizeof(int));
   }

   gset = GeoRef_SetGet(ToRef, LayerRef, Opt);

   // Do we have and index
   if (gset->Index && gset->Index[0] != REF_INDEX_EMPTY) {

      // As long as the file or the list is not empty
      ip = gset->Index;
      while(*ip != REF_INDEX_END) {

         // Get the gridpoint
         f = *(ip++);
         area = *(ip++);
         value = *(ip++);

         if (f >= Layer->NFeature) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Wrong index, feature index too high (%i)\n", __func__, f);
            return FALSE;
         }

         // Get value to distribute
         if (fld >= 0) {
            value = OGR_F_GetFieldAsDouble(Layer->Feature[f], fld);
         } else if (fld < -9) {
            value = Value;
         }

         // Get the geometry intersections
         while(*ip != REF_INDEX_SEPARATOR) {
            pi = *(ip++);
            pj = *(ip++);
            dp = *(ip++);

            idx2 = FIDX2D(ToDef, pi, pj);

            Def_Get(ToDef, 0, idx2, val);
            if (!DEFVALID(ToDef, val)) {
               val = 0.0;
            }

            val += value*dp;

            if (Opt->Combine == CB_AVERAGE)  {
               ToDef->Accum[idx2] += 1;
            }
            Def_Set(ToDef, 0, idx2, val);
         }
         // Skip separator
         ip++;
      }
   } else {

      // Define the max size of the indexes
      isize = 1024;
      if ((c = getenv("GEOREF_INDEX_SIZE_HINT"))) {
         isize = atoi(c);
      }

      if (gset->Index && gset->Index[0] == REF_INDEX_EMPTY) {
         if (!(index = (float**)malloc(Layer->NFeature*sizeof(float*)))) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable to allocate local index arrays\n", __func__);
            return FALSE;
         }
         ip = gset->Index;
      }

      // If the request is in meters
      if (fld == -3 || fld == -5) {
         // Create an UTM referential and transform to convert to meters
         x0 = LayerRef->X0;
         y0 = LayerRef->Y0;
         LayerRef->XY2LL(LayerRef, &x0, &y0, &co.Lat, &co.Lon, 1);
         srs = OSRNewSpatialReference(NULL);
         OSRSetUTM(srs, (int)ceil((180+co.Lon)/6), (int)co.Lat);
         tr = OCTNewCoordinateTransformation(LayerRef->Spatial, srs);

         if (!srs || !tr) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Could not initiate UTM transormation\n", __func__);
            return FALSE;
         }
      }

      // Trouve la feature en intersection
      #pragma omp parallel for private(f, geom, hgeom, utmgeom, env, co, value, vr, n, area, mode, type, lp) firstprivate(pick, poly) shared(Layer, LayerRef, ToRef, fld, tr, error, ip, index, isize) reduction(+ : nt)
      for(f = 0; f < Layer->NFeature; f++) {

         if (index) index[f] = NULL;
         if (error) continue;
         n = 0;

         if (Layer->Select[f] && Layer->Feature[f]) {

            geom = utmgeom = NULL;
            // Try to access geometry, skipping instead of failing on bad ones
            if (!(hgeom = OGR_F_GetGeometryRef(Layer->Feature[f]))) {
               Lib_Log(APP_LIBGEOREF, APP_WARNING, "%s: Cannot get handle from geometry: %li\n", __func__, f);
               continue;
            }

            // Copie de la geometrie pour transformation
            if (!(geom = OGR_G_Clone(hgeom))) {
               Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Could not clone the geometry\n", __func__);
               error = 1;
               continue;
            }

            if (!pick) {
               pick = OGR_G_CreateGeometry(wkbLinearRing);
               poly = OGR_G_CreateGeometry(wkbPolygon);
               OGR_G_AddGeometryDirectly(poly, pick);
            }

            // If the request is in meters
            if (fld == -3 || fld == -5) {
               if (!(utmgeom = OGR_G_Clone(geom))) {
                  Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Could not clone the UTM geomtry\n", __func__);
                  OGR_G_DestroyGeometry(geom);
                  error = 1;
                  continue;
               }

               // Transform the geom to utm
               OGR_G_Transform(utmgeom, tr);
            }

            // Get value to distribute
            if (fld >= 0) {
               value = OGR_F_GetFieldAsDouble(Layer->Feature[f], fld);
            } else if (fld == -2) {
               value = OGR_G_Area(geom);
            } else if (fld == -3) {
               value = OGR_G_Area(utmgeom);
            } else if (fld == -4) {
               value = OGM_Length(geom);
            } else if (fld == -5) {
               value = OGM_Length(utmgeom);
            } else if (fld == -6) {
               value = f;
            } else if (fld == -7) {
               value = OGM_CoordLimit(geom, 2, 0);
            } else if (fld == -8) {
               value = OGM_CoordLimit(geom, 2, 1);
            } else if (fld == -9) {
               value = OGM_CoordLimit(geom, 2, 2);
            } else {
               value = Value;
            }

            // In centroid mode, just project the coordinate into field and set value
            if (Opt->Interp == IV_CENTROID) {
               OGM_Centroid2D(geom, &vr[0], &vr[1]);
               GeoRef_XY2LL(LayerRef, &co.Lat, &co.Lon, &vr[0], &vr[1], 1, TRUE);
               GeoRef_LL2XY(ToRef, &vr[0], &vr[1], &co.Lat, &co.Lon, 1, TRUE);
               n = (int)vr[1]*ToDef->NI+(int)vr[0];
               Def_Set(ToDef, 0, n, value);
               nt++;
            } else {

               // Transform geometry to field referential
               OGM_OGRProject(geom, LayerRef, ToRef);

               // Use enveloppe limits to initialize the initial lookup range
               OGR_G_GetEnvelope(geom, &env);
               if (!(env.MaxX < (ToRef->X0-0.5) || env.MinX > (ToRef->X1+0.5) || env.MaxY < (ToRef->Y0-0.5) || env.MinY > (ToRef->Y1+0.5))) {

                  if (Opt->Interp == IV_FAST) {
                     GeoRef_Rasterize(ToRef, ToDef, Opt, geom, value);
                  } else {

                     // Get value to split on
                     area = -1.0;
                     mode = 'N';
                     type = 'A';
                     switch(Opt->Interp) {
                        case IV_FAST                           : break;
                        case IV_WITHIN                         : mode = 'W'; type = 'A'; break;
                        case IV_INTERSECT                      : mode = 'I'; type = 'A'; break;
                        case IV_CENTROID                       : break;
                        case IV_ALIASED                        : mode = 'A'; type = 'A'; area = 1.0; break;
                        case IV_CONSERVATIVE                   : mode = 'C'; type = 'A'; area = OGR_G_Area(geom); break;
                        case IV_NORMALIZED_CONSERVATIVE        : mode = 'N'; type = 'A'; area = OGR_G_Area(geom); break;
                        case IV_LENGTH_CONSERVATIVE            : mode = 'C'; type = 'L'; area = OGM_Length(geom); break;
                        case IV_LENGTH_NORMALIZED_CONSERVATIVE : mode = 'N'; type = 'L'; area = OGM_Length(geom); break;
                        case IV_LENGTH_ALIASED                 : mode = 'A'; type = 'L'; area = OGM_Length(geom); break;
                        case IV_POINT_CONSERVATIVE             : mode = 'C'; type = 'P'; area = 1.0; break;
                        default: break;
                     }

                     // If it's nil then nothing to distribute on
                     if (area > 0.0 || Opt->Interp <= IV_CENTROID) {

                        env.MaxX += 0.5; env.MaxY += 0.5;
                        env.MinX = env.MinX < 0 ? 0 : env.MinX;
                        env.MinY = env.MinY < 0 ? 0 : env.MinY;
                        env.MaxX = env.MaxX > (ToRef->X1+0.5) ? (ToRef->X1+0.5) : env.MaxX;
                        env.MaxY = env.MaxY > (ToRef->Y1+0.5) ? (ToRef->Y1+0.5) : env.MaxY;
                        area = area < 0 ? 0.0 : area;

                        // Append feature into index
                        lp = NULL;
                        if (ip) {
                           if (!(index[f] = (float*)malloc(isize*sizeof(float)))) {
                              Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable to allocate local index memory (%i)\n", __func__, f);
                              error = 1;
                              continue;
                           }
                           lp = index[f];
                        }
                        if (lp) {
                           *(lp++) = area;
                           *(lp++) = (fld < 0 && fld >= -9) ? value : -999.0;
                        }

                        nt += n = GeoRef_InterpQuad(ToRef, ToDef, Opt, poly, geom, mode, type, area, value, env.MinX, env.MinY, env.MaxX, env.MaxY, 0, &lp);

                        if (lp) {
                           if (n) {
                              *(lp++) = REF_INDEX_SEPARATOR;
                           } else {
                              index[f][0] = REF_INDEX_EMPTY;  // No intersection found, removed previously inserted feature
                           }
                        }
                        Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s: %i hits on feature %i of %i (%.0f %.0f x %.0f %.0f)\n", __func__, n, f, Layer->NFeature, env.MinX, env.MinY, env.MaxX, env.MaxY);
                     }
                  }
               }
            }
            if (geom)    OGR_G_DestroyGeometry(geom);
            if (utmgeom) OGR_G_DestroyGeometry(utmgeom);
         }
      }

      Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s: %i total hits\n", __func__, nt);

      // Merge indexes
      n = 0;
      if (ip && nt && !error) {
         for(f = 0; f < Layer->NFeature; f++) {
            if ((lp = index[f]) && *lp != REF_INDEX_EMPTY) {
               *(ip++) = f;
               while(*lp != REF_INDEX_SEPARATOR) {
                  *(ip++) = *(lp++);
               }
               *(ip++) = REF_INDEX_SEPARATOR;
            }
            if (index[f]) free(index[f]);
         }
         *(ip++) = REF_INDEX_END;
         gset->IndexSize = (ip-gset->Index)+1;
         free(index);
      }

      if (tr)
         OCTDestroyCoordinateTransformation(tr);

      if (srs)
         OSRDestroySpatialReference(srs);
   }

   if (Opt->Combine == CB_AVERAGE) {
      for(n = 0; n < FSIZE2D(ToDef); n++) {
         if (ToDef->Accum[n] != 0.0) {
            Def_Get(ToDef, 0, n, value);
            value /= ToDef->Accum[n];
            Def_Set(ToDef, 0, n, value);
         }
      }
   }

   // Return size of index or number of hits, or 1 if nothing found
   nt = gset->Index ? gset->IndexSize : nt;
   return (error || nt == 0) ? 1 : nt;
#else
   Lib_Log(APP_LIBGEOREF, APP_ERROR, "Function %s is not available, needs to be built with GDAL\n", __func__);
   return FALSE;
#endif
}

/*----------------------------------------------------------------------------
 * @brief  Interpolate at sub grid resolution into a temporary "sub" buffer for use by other functions
 *    like sub grid variance.
 * @date   August 2013
    //! [in]     ToRef        Destination GeoRef pointer
    //! [in]     ToDef        Destination data definition
    //! [in]     FromRef      Source GeoRef pointer
 *
 *    @return        FALSE on error, TRUE otherwise
*/
int32_t GeoRef_InterpSub(TGeoRef *ToRef, TDef *ToDef, TGeoRef *FromRef, TDef *FromDef, TGeoOptions *Opt) {

   int32_t  idx, i, j, x, y, x0, x1, y0, y1, s;
   double val, val1, di, dj, d, la, lo;

   if (!Opt) Opt=&GeoRef_Options;

   if (!ToRef || !ToDef) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid destination\n", __func__);
      return 0;
   }
   if (!FromRef || !FromDef) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid source\n", __func__);
      return 0;
   }

   if (Opt->Sampling <= 1) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Sub sampling is not defined\n", __func__);
      return 0;
   }

   if (!GeoRef_Intersect(FromRef, ToRef, &x0, &y0, &x1, &y1, 0)) {
      return 1;
   }

   // Allocate and initialize subgrid
   s = Opt->Sampling*Opt->Sampling;
   if (!ToDef->Sub) {
      if ((ToDef->Sub = (float*)malloc(ToDef->NI*ToDef->NJ*s*sizeof(float)))) {
         for(idx = 0; idx < ToDef->NI*ToDef->NJ*s; idx++) {
            ToDef->Sub[idx] = ToDef->NoData;
         }
      } else {
         Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable ot allocate subgrid\n", __func__);
         return 0;
      }
   }

   d = 1.0/(Opt->Sampling-1);

   // Interpolate values on subgrid
   for(y = y0; y <= y1; y++) {
      for(x = x0; x <= x1; x++) {
         idx = (y*ToDef->NI+x)*s;

         for(j = 0; j < Opt->Sampling; j++) {
            for(i = 0; i < Opt->Sampling; i++, idx++) {

               di = (double)x-0.5+d*i;
               dj = (double)y-0.5+d*j;

               GeoRef_XY2LL(FromRef, &la, &lo, &di, &dj, 1, TRUE);
               if (GeoRef_LL2XY(ToRef, &di, &dj, &la, &lo, 1, FALSE)) {
                  if (Def_GetValue(FromRef, FromDef, Opt, 0, di, dj, FromDef->Level, &val, &val1) && DEFVALID(FromDef, val)) {
                     ToDef->Sub[idx] = val;
                  }
               }
            }
         }
      }
   }

   return 1;
}

/*----------------------------------------------------------------------------
 * @brief  Conservative interpolations
 * @date   May 2006
    //! [in]     ToRef        Destination GeoRef pointer
    //! [in]     ToDef        Destination data definition
    //! [in]     FromRef      Source GeoRef pointer
    //! [in]     FromDef      Source data definition
    //! [in]     Opt          Interpolation options
 *
 *    @return        FALSE on error, TRUE otherwise
*/
int32_t GeoRef_InterpConservative(TGeoRef *ToRef, TDef *ToDef, TGeoRef *FromRef, TDef *FromDef, TGeoOptions *Opt) {

#ifdef HAVE_GDAL
   TGeoSet    *gset = NULL;
   int32_t     i, j, na, nt = 0, p = 0, pi, pj, idx2, idx3, intersect, k = 0, isize, nidx, error = 0;
   uint64_t    n;
   char        *c;
   double       val0, val1, area, x, y, z, dp;
   float       *ip = NULL, *lp = NULL, **index = NULL;
   OGRGeometryH cell = NULL, ring = NULL, *pick = NULL, *poly = NULL;
   OGREnvelope  env;

   if (!Opt) Opt=&GeoRef_Options;

   Opt->Combine = CB_SUM;
   if (!ToRef || !ToDef) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid destination\n", __func__);
      return FALSE;
   }
   if (!FromRef || !FromDef) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid source\n", __func__);
      return FALSE;
   }

   // Allocate area buffer if needed
   if (Opt->Interp == IR_NORMALIZED_CONSERVATIVE && !ToDef->Buffer) {
      ToDef->Buffer = (double*)malloc(FSIZE3D(ToDef)*sizeof(double));
      if (!ToDef->Buffer) {
         Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable to allocate area buffer\n", __func__);
         return FALSE;
      }
      for(n = 0; n < FSIZE3D(ToDef); n++) ToDef->Buffer[n] = 0.0;
   }

   // Define the max size of the indexes
   isize = 256;
   if ((c = getenv("GEOREF_INDEX_SIZE_HINT"))) {
      isize = atoi(c);
   }

   gset = GeoRef_SetGet(ToRef, FromRef, Opt);
   if (isize>0) 
      GeoRef_SetIndexInit(gset);

   // Process one level at a time
   for (k = 0; k < ToDef->NK; k++) {

      // Do we have and index
      if (gset->Index && gset->Index[0] != REF_INDEX_EMPTY) {

         // As long as the file or the list is not empty
         ip = gset->Index;
         while(*ip != REF_INDEX_END) {

            // Get the gridpoint
            i = *(ip++);
            j = *(ip++);

            if (!FIN2D(FromDef, i, j)) {
               Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Wrong index, index coordinates (%i, %i)\n", __func__, i, j);
               return FALSE;
            }

            // Get this gridpoint value
            Def_Get(FromDef, 0, FIDX3D(FromDef, i, j, k), val1);
            if (!DEFVALID(FromDef, val1)) {
               continue;
            }

            // Get the geometry intersections
            while(*ip != REF_INDEX_SEPARATOR) {
               pi = *(ip++);
               pj = *(ip++);
               dp = *(ip++);

               idx2 = FIDX2D(ToDef, pi, pj);
               idx3 = FIDX3D(ToDef, pi, pj, k);
               Def_Get(ToDef, 0, idx3, val0);
               if (!DEFVALID(ToDef, val0))
                  val0 = 0.0;

               // Assign new value
               val0 += val1*dp;
               Def_Set(ToDef, 0, idx3, val0);

               if (Opt->Interp == IR_NORMALIZED_CONSERVATIVE) {
                  ToDef->Buffer[idx3] += dp;
               }
            }
            // Skip separator
            ip++;
         }
      } else {

         // U grids need mask to select on 1 grid on overlap zones
         if (FromRef->NbSub>0) {
            GeoRef_SetCalcYYXY(gset);
         }

         if (gset->Index && gset->Index[0] == REF_INDEX_EMPTY) {
            if (!(index = (float**)malloc(FSIZE2D(FromDef)*sizeof(float*)))) {
               Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable to allocate local index arrays\n", __func__);
               return FALSE;
            }
            ip = gset->Index;
         }

         #pragma omp parallel for collapse(2) firstprivate(cell, ring, pick, poly) private(i, j, nidx, intersect, cnt, p, x, y, z, lp, area, val1, env, n, na) shared(k, isize, FromRef, FromDef, ToRef, ToDef, Opt, error, index, ip) reduction(+ : nt)
         for(j = 0; j < FromDef->NJ; j++) {
            for(i = 0; i < FromDef->NI; i++) {
 
               nidx = j*FromDef->NI+i;
               if (index) index[nidx] = NULL;
               if (error) continue;

               // Create gridcell geometry object
               if (!cell) {
                  cell = OGR_G_CreateGeometry(wkbPolygon);
                  ring = OGR_G_CreateGeometry(wkbLinearRing);
                  OGR_G_AddGeometryDirectly(cell, ring);
               }

               if (!pick) {
                  pick = OGR_G_CreateGeometry(wkbLinearRing);
                  poly = OGR_G_CreateGeometry(wkbPolygon);
                  OGR_G_AddGeometryDirectly(poly, pick);
               }

               // Project the source gridcell into the destination
               if (!(intersect = GeoRef_Cell2OGR(ring, ToRef, FromRef, i, j, Opt)))
                  continue;

               // Allocate local index
               lp = NULL;
               if (ip) {
                  if (!(index[nidx] = (float*)malloc(isize*sizeof(float)))) {
                     Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable to allocate local index memory\n", __func__);
                     error = 1;
                     continue;
                  }
                  index[nidx][isize-1] = index[nidx][isize-2]  = index[nidx][isize-3] = REF_INDEX_NIL;
                  index[nidx][0] = REF_INDEX_EMPTY;
                  lp = index[nidx];
               }
               n = na = 0;

               // Are we crossing the wrap around
               if (intersect < 0) {
                  // If so, move the wrapped points (assumed greater than NI/2) to the other side
                  for(p = 0; p < -intersect; p++) {
                     OGR_G_GetPoint(ring, p, &x, &y, &z);
                     if (x > ToDef->NI >> 1) {
                        x -= ToDef->NI;
                        OGR_G_SetPoint_2D(ring, p, x, y);
                     }
                  }

                  // Process the cell
                  area = OGR_G_Area(cell);

                  if (area > 0.0) {

                     Def_Get(FromDef, 0, FIDX3D(FromDef, i, j, k), val1);
                     // If we are saving the indexes, we have to process even if nodata but use 0.0 so as to not affect results
                     if (!DEFVALID(FromDef, val1)) {
                        if (lp) {
                           val1 = 0.0;
                        } else {
                           continue;
                        }
                     }

                     // Use enveloppe limits to initialize the initial lookup range
                     OGR_G_GetEnvelope(ring, &env);
                     env.MaxX += 0.5; env.MaxY += 0.5;
                     env.MinX = env.MinX < 0 ? 0 : env.MinX;
                     env.MinY = env.MinY < 0 ? 0 : env.MinY;
                     env.MaxX = env.MaxX > ToRef->X1 ? ToRef->X1 : env.MaxX;
                     env.MaxY = env.MaxY > ToRef->Y1 ? ToRef->Y1 : env.MaxY;

                     nt += na = GeoRef_InterpQuad(ToRef, ToDef, Opt, poly, cell, Opt->Interp == IR_CONSERVATIVE ? 'C' : 'N', 'A', area, val1, env.MinX, env.MinY, env.MaxX, env.MaxY, k, &lp);

                     Lib_Log(APP_LIBGEOREF, APP_EXTRA, "%s: %i hits on grid point %i %i (%.0f %.0f x %.0f %.0f)\n", __func__, na, i, j, env.MinX, env.MinY, env.MaxX, env.MaxY);
                  }

                  // We have to process the part that was out of the grid limits so translate everything NI points
                  for(p = 0; p < -intersect; p++) {
                     OGR_G_GetPoint(ring, p, &x, &y, &z);
                     x += ToDef->NI;
                     OGR_G_SetPoint_2D(ring, p, x, y);
                  }
               }

               area = OGR_G_Area(cell);

               if (area > 0.0) {
                  Def_Get(FromDef, 0, FIDX3D(FromDef, i, j, k), val1);
                  if (!DEFVALID(FromDef, val1)) {
                     // If we are saving the indexes, we have to process even if nodata but use 0.0 so as to not affect results
                     if (lp) {
                        val1 = 0.0;
                     } else {
                        continue;
                     }
                  }

                  // Use envelop limits to initialize the initial lookup range
                  OGR_G_GetEnvelope(ring, &env);
                  if (!(env.MaxX < ToRef->X0 || env.MinX > ToRef->X1 || env.MaxY < ToRef->Y0 || env.MinY > ToRef->Y1)) {
                     env.MaxX += 0.5; env.MaxY += 0.5;
                     env.MinX = env.MinX < 0 ? 0 : env.MinX;
                     env.MinY = env.MinY < 0 ? 0 : env.MinY;
                     env.MaxX = env.MaxX > ToRef->X1 ? ToRef->X1 : env.MaxX;
                     env.MaxY = env.MaxY > ToRef->Y1 ? ToRef->Y1 : env.MaxY;

                     nt += n = GeoRef_InterpQuad(ToRef, ToDef, Opt, poly, cell, Opt->Interp == IR_CONSERVATIVE ? 'C' : 'N', 'A', area, val1, env.MinX, env.MinY, env.MaxX, env.MaxY, k, &lp);

                     Lib_Log(APP_LIBGEOREF, APP_EXTRA, "%s: %i hits on grid point %i %i (%.0f %.0f x %.0f %.0f)\n", __func__, n, i, j, env.MinX, env.MinY, env.MaxX, env.MaxY);
                  }
               }
               if (lp && (n || na)) {
                  *(lp++) = REF_INDEX_SEPARATOR; // End the list for this gridpoint
               }
            }
         }

         // Merge indexes
         n = 0;
         if (ip && nt && !error) {
           if (gset->IndexSize < nt*3) {
               gset->IndexSize = nt*3+FSIZE2D(FromDef)*3+1;
               ip = gset->Index = (float*)realloc(gset->Index, gset->IndexSize*sizeof(float));
            }
            for(j = 0; j < FromDef->NJ; j++) {
               for(i = 0; i < FromDef->NI; i++) {
                  nidx = j*FromDef->NI+i;

                  if ((lp = index[nidx]) && *lp != REF_INDEX_EMPTY) {
                     // Append gridpoint to the index
                     *(ip++) = i;
                     *(ip++) = j;

                     // This gridpoint wraps around
                     while(*lp != REF_INDEX_SEPARATOR) {
                        *(ip++) = *(lp++);
                     }
                     *(ip++) = REF_INDEX_SEPARATOR;
                  }
                  if (index[nidx]) free(index[nidx]);
               }
            }
            *(ip++) = REF_INDEX_END;
            gset->IndexSize = (ip-gset->Index)+1;
            free(index);
         }
         Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s: %i total hits\n", __func__, nt);
      }
   }

//   OGR_G_DestroyGeometry(ring);
//   OGR_G_DestroyGeometry(cell);

   // Return size of index or number of hits, or 1 if nothing found
   nt = gset->Index ? gset->IndexSize : nt;
   return nt == 0 ? 1 : nt;
#else
   Lib_Log(APP_LIBGEOREF, APP_ERROR, "Function %s is not available, needs to be built with GDAL\n", __func__);
   return FALSE;
#endif

}

/*----------------------------------------------------------------------------
 * @brief  Interpolations by averaging, minimum or maximum
 * @date   May 2006
    //! [in]     ToRef        Destination GeoRef pointer
    //! [in]     ToDef        Destination data definition
    //! [in]     FromRef      Source GeoRef pointer
    //! [in]     FromDef      Source data definition
    //! [in]     Opt          Interpolation options
 *
 *    @return        FALSE on error, TRUE otherwise
*/
int32_t GeoRef_InterpAverage(TGeoRef *ToRef, TDef *ToDef, TGeoRef *FromRef, TDef *FromDef, TGeoOptions *Opt){

   double        val, vx, di[4], dj[4], *fld, *aux, di0, di1, dj0, dj1;
   float        *ip = NULL;
   int32_t      *acc = NULL, x0, x1, y, y0, y1, isize;
   unsigned long idxt, idxk, idxj, n, nijk, nij;
   uint32_t      n2, ndi, ndj, k, t, s, x, dx, dy;
   TGeoSet      *gset;
   TGeoScan      gscan;
   TRef_Interp  interp;

   if (!Opt) Opt=&GeoRef_Options;

   if (!ToRef || !ToDef) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid destination\n", __func__);
      return FALSE;
   }
   if ((!FromRef || !FromDef) && (Opt->Interp != IR_NOP && Opt->Interp != IR_ACCUM && Opt->Interp != IR_BUFFER)) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid source\n", __func__);
      return FALSE;
   }

   acc = ToDef->Accum;
   fld = ToDef->Buffer;
   aux = ToDef->Aux;
   nij = FSIZE2D(ToDef);
   nijk = FSIZE3D(ToDef);
   val = vx = 0.0;
   interp = Opt->Interp;

   if (Opt->Interp != IR_NOP && Opt->Interp != IR_ACCUM && Opt->Interp != IR_BUFFER) {
      if (!GeoRef_Intersect(ToRef, FromRef, &x0, &y0, &x1, &y1, 0)) {
         return TRUE;
      }
      // TODO: figure out why GeoRef_Intersect doe not give right limits of U grids
      x0=FromRef->X0;
      x1=FromRef->X1;
      y0=FromRef->Y0;
      y1=FromRef->Y1;

      // In case of average, we need an accumulator
      if (Opt->Interp == IR_AVERAGE || Opt->Interp == IR_VECTOR_AVERAGE || Opt->Interp == IR_VARIANCE ||
          Opt->Interp == IR_SQUARE || Opt->Interp == IR_NORMALIZED_COUNT || Opt->Interp == IR_COUNT) {
         if (!ToDef->Accum) {
            acc = ToDef->Accum = calloc(nij, sizeof(int32_t));
            if (!ToDef->Accum) {
               Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable to allocate accumulation buffer\n", __func__);
               return FALSE;
            }
         }
      }

      if (!ToDef->Buffer) {
         fld = ToDef->Buffer = malloc(nijk*sizeof(double));
         if (!ToDef->Buffer) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable to allocate buffer\n", __func__);
            return FALSE;
         }
         for(n = 0; n < nijk; n++) fld[n] = ToDef->NoData;

         if (Opt->Interp == IR_VECTOR_AVERAGE) {
            aux = ToDef->Aux = malloc(nijk*sizeof(double));
            if (!ToDef->Aux) {
               Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable to allocate auxiliary buffer\n", __func__);
               return FALSE;
            }
            for(n = 0; n < nijk; n++) aux[n] = 0.0;
         }
      }

      if (ToRef->GRTYP[0] == 'Y') {
         // Point cloud interpolations

         for(idxt = 0; idxt < nij; idxt++) {
            if (GeoRef_LL2XY(FromRef, &di0, &dj0, &ToRef->Lat[idxt], &ToRef->Lon[idxt], 1, FALSE)) {
               di0 = floor(di0);
               dj0 = floor(dj0);

               Opt->Interp = IR_NEAREST;
               Def_GetValue(FromRef, FromDef, Opt, 0, di0, dj0, FromDef->Level, &di[0], &vx);
               Def_GetValue(FromRef, FromDef, Opt, 0, di0+1.0, dj0, FromDef->Level, &di[1], &vx);
               Def_GetValue(FromRef, FromDef, Opt, 0, di0, dj0+1.0, FromDef->Level, &di[2], &vx);
               Def_GetValue(FromRef, FromDef, Opt, 0, di0+1.0, dj0+1.0, FromDef->Level, &di[3], &vx);
               Opt->Interp = interp;
            }
            for(s = 0; s < 4; s++) {
               vx = di[s];
               if (!DEFVALID(FromDef, vx))
                  continue;

               // If the previous value is nodata, initialize the counter
               if (!DEFVALID(FromDef, fld[idxt])) {
                  fld[idxt] = (Opt->Interp == IR_SUM || Opt->Interp == IR_AVERAGE|| Opt->Interp == IR_VECTOR_AVERAGE) ? 0.0 : (Opt->Interp == IR_MAXIMUM ? -HUGE_VAL : HUGE_VAL);
                  if (aux) aux[idxt] = fld[idxt];
               }
               switch(Opt->Interp) {
                  case IR_MAXIMUM          : if (vx > fld[idxt]) fld[idxt] = vx; break;
                  case IR_MINIMUM          : if (vx < fld[idxt]) fld[idxt] = vx; break;
                  case IR_SUM              : fld[idxt] += vx;                  break;
                  case IR_AVERAGE          : fld[idxt] += vx; acc[idxt]++;     break;
                  case IR_VECTOR_AVERAGE   : vx = DEG2RAD(vx); fld[idxt] += cos(vx); aux[idxt] += sin(vx); acc[idxt]++; break;
                  default:
                     Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid interpolation type\n", __func__);
                     return FALSE;
               }
            }
         }
      } else {
         int32_t   *fromClass = NULL;
         int32_t   *toClass = NULL;
         int32_t    nbClass = 0;
         int32_t   *rpnClass = NULL;
         int32_t    i;

         if (Opt->lutDef && ToDef->NK > 0) {
            val = Opt->lutDef[0][0];
            if (val != ToDef->NK) {
               Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid LUT class size (%d) vs Array size (%d)  \n", __func__, (int)val, ToDef->NK);
               return FALSE;
            }
            // first row is the class count and idno
            nbClass = Opt->lutSize - 1;
            fromClass = (int32_t *)malloc( sizeof(int)*nbClass );
            for (i = 0; i < nbClass ; i++) {
               fromClass[i] = Opt->lutDef[0][i+1];
            }
            if (Opt->lutDim == 2) { // only FROM and TO
               toClass = (int32_t *)malloc( sizeof(int)*nbClass );
               for (i = 0; i < nbClass ; i++) {
                  val = Opt->lutDef[1][i+1];
                  // make sure to value is within range
                  if ((val >= 1) && (val <= ToDef->NK))
                     toClass[i] = val;
                  else
                     toClass[i] = -1;
               }
            } else {
               rpnClass = (int32_t *)malloc(sizeof(int)*Opt->lutDim);
               for (i = 1; i < Opt->lutDim ; i++) { // skip 1st column : the class id
                  val = Opt->lutDef[i][0];
                  rpnClass[i] = val;
                  if ((val < 1)||(val > ToDef->NK)) {
                     Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid LUT index(%d) = %f vs Array size (%d)  \n", __func__, i, val, ToDef->NK);
                     return FALSE;
                  }
               }
            }
         }

         gset = GeoRef_SetGet(ToRef, FromRef, Opt);

         // Do we have and index
         if (gset->Index && gset->Index[0] != REF_INDEX_EMPTY) {

            // As long as the file or the list is not empty
            ip = gset->Index;
            while(*ip != REF_INDEX_END) {

               // Get the gridpoint
               vx = *(ip++);
               if (vx == REF_INDEX_SEPARATOR)
                  continue;
 
               // Skip if no data
               Def_Get(FromDef, 0, (int)vx, vx);
               if (!DEFVALID(FromDef, vx) && Opt->Interp != IR_COUNT)
                  continue;

               while(*ip!=REF_INDEX_END && *ip!=REF_INDEX_SEPARATOR) {

                  idxt = *(ip++);

                  // If the previous value is nodata, initialize the counter
                  if (!DEFVALID(ToDef, fld[idxt])) {
                     fld[idxt] = (Opt->Interp == IR_SUM || Opt->Interp == IR_AVERAGE  || Opt->Interp == IR_VECTOR_AVERAGE || Opt->Interp == IR_VARIANCE || Opt->Interp == IR_SQUARE || Opt->Interp == IR_NORMALIZED_COUNT || Opt->Interp == IR_COUNT) ? 0.0 : (Opt->Interp == IR_MAXIMUM ? -HUGE_VAL : HUGE_VAL);
                     if (aux) aux[idxt] = fld[idxt];
                  }
                                       
                  switch(Opt->Interp) {
                     case IR_MAXIMUM          : if (vx > fld[idxt]) fld[idxt] = vx;
                                                break;
                     case IR_MINIMUM          : if (vx < fld[idxt]) fld[idxt] = vx;
                                                break;
                     case IR_SUM              : fld[idxt] += vx;
                                                break;
                     case IR_VARIANCE         : acc[idxt]++;
                                                val = Opt->Ancilliary[idxt];
                                                fld[idxt] += (vx-val)*(vx-val);
                                                break;
                     case IR_SQUARE           : acc[idxt]++;
                                                fld[idxt] += vx*vx;
                                                break;
                     case IR_COUNT            : acc[idxt]++;
                     case IR_AVERAGE          :
                     case IR_NORMALIZED_COUNT : if (Opt->Table) {
                                                   t = 0;
                                                   while(t < ToDef->NK) {
                                                      if (vx == Opt->Table[t]) {
                                                         if (Opt->Interp != IR_COUNT) acc[idxt]++;
                                                         fld[t*nij+idxt] += 1.0;
                                                         break;
                                                      }
                                                      t++;
                                                   }
                                                } else if (Opt->lutDef) {
                                                   int32_t hasvalue = 0;
                                                   t = 0;
                                                   while(t < nbClass) {
                                                      if (vx == fromClass[t]) {
                                                         if (toClass) { // toClass is between 1 and 26
                                                            if (toClass[t] > 0)
                                                               {
                                                               fld[(toClass[t]-1)*nij+idxt] += 1.0;
                                                               hasvalue = 1;
                                                               }
                                                         } else {
                                                            for (i = 1; i < Opt->lutDim ; i++) { // skip 1st column : the class id
                                                               val = Opt->lutDef[i][t+1];
                                                               if (val > 0) {
                                                                  fld[(rpnClass[i]-1)*nij+idxt] += val;
                                                                  hasvalue = 1;
                                                               }
                                                            }
                                                         }
                                                         break;
                                                      }
                                                      t++;
                                                   }
                                                   // Dont count missing values
                                                   if (hasvalue)
                                                      if (Opt->Interp != IR_COUNT) acc[idxt]++;
                                                } else {
                                                   if (DEFVALID(FromDef, vx)) {
                                                      fld[idxt] += vx;
                                                      if (Opt->Interp != IR_COUNT) acc[idxt]++;
                                                   }
                                                }
                                                break;
                     case IR_VECTOR_AVERAGE   : vx = DEG2RAD(vx); fld[idxt] += cos(vx); aux[idxt] += sin(vx); acc[idxt]++; break;
                     default:
                        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid interpolation type\n", __func__);
                        return FALSE;
                  }
               }
            }
         } else {
            char *c;

            // Define the max size of the indexes
            GeoRef_SetIndexInit(gset);
            if (gset->Index && gset->Index[0] == REF_INDEX_EMPTY) {
               ip = gset->Index;
            }


            // grid based interpolations
            GeoScan_Init(&gscan);

            // if > 2048x2048, loop by lines otherwise, do it in one shot
            n2 = ToDef->NI >> 1;
            dy = ((y1-y0)*(x1-x0)) > 4194304 ? 0 : (y1-y0);
            for(y = y0; y <= y1; y += (dy+1)) {                 

               if (!(s = GeoScan_Get(&gscan, ToRef, NULL, FromRef, FromDef, Opt, x0, y, x1, y+dy))) {
                  Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable to allocate coordinate scanning buffer\n", __func__);
                  return FALSE;
               }

               // Loop over source data
               dx = 0;
               for(x = 0, n = 0; x < gscan.N; x++, n++) {

                  // Check if we need to skip last x since we change row and last one is end of a cell
                  if (s > 1 && dx == gscan.DX) {
                     n++;
                     dx = 0;
                  }
                  dx++;

                  // Skip if no mask
                  if (FromDef->Mask && !FromDef->Mask[gscan.V[x]])
                     continue;

                  if (x>0) *(ip++) = REF_INDEX_SEPARATOR;
                  *(ip++) = gscan.V[x];

                  // Skip if no data
                  Def_Get(FromDef, 0, gscan.V[x], vx);
                  if (!DEFVALID(FromDef, vx) && Opt->Interp != IR_COUNT)
                     continue;

                  // Figure out ordered coverage
                  if (s > 1) {
                     di[0] = gscan.X[n];
                     dj[0] = gscan.Y[n];
                     di[1] = gscan.X[n+1]; di0 = fmin(di[0], di[1]); di1 = fmax(di[0], di[1]);
                     dj[1] = gscan.Y[n+1]; dj0 = fmin(dj[0], dj[1]); dj1 = fmax(dj[0], dj[1]);

                     di[2] = gscan.X[n+gscan.DX+1]; di0 = fmin(di0, di[2]); di1 = fmax(di1, di[2]);
                     dj[2] = gscan.Y[n+gscan.DX+1]; dj0 = fmin(dj0, dj[2]); dj1 = fmax(dj1, dj[2]);
                     di[3] = gscan.X[n+gscan.DX+2]; di0 = fmin(di0, di[3]); di1 = fmax(di1, di[3]);
                     dj[3] = gscan.Y[n+gscan.DX+2]; dj0 = fmin(dj0, dj[3]); dj1 = fmax(dj1, dj[3]);

                     di0 = ROUND(di0); dj0 = ROUND(dj0);
                     di1 = ROUND(di1); dj1 = ROUND(dj1);
                  } else {
                     di0 = di1 = ROUND(gscan.X[n]);
                     dj0 = dj1 = ROUND(gscan.Y[n]);
                  }
      
                  // Are we within the destination field
                  if (di0 >= ToDef->NI || dj0 >= ToDef->NJ || di1 < 0 || dj1 < 0)
                     continue;

                  // Test for polar outsidness (Problem we had with yinyang grids)
                  if ((di0 < 0 && di1 > ToDef->NI) || (dj0 < 0 && dj1 > ToDef->NJ))
                     continue;

                  // Clamp the coordinates
                  if (di0 < 0) di0 = 0;
                  if (dj0 < 0) dj0 = 0;
                  if (di1 > ToDef->NI-1) di1 = ToDef->NI-1;
                  if (dj1 > ToDef->NJ-1) dj1 = ToDef->NJ-1;

                  // Are we crossing the wrap around
                  if (ToRef->Type&GRID_WRAP && di0 < n2 && di1 > n2 && (di1-di0) > n2) {
                     val = di0;
                     di0 = di1;
                     di1 = val+ToDef->NI;
                  }

                  for(ndj = dj0; ndj <= dj1; ndj++) {
                     idxj = ndj*ToDef->NI;
                     for(ndi = di0; ndi <= di1; ndi++) {
                        idxt = idxj+(ndi >= ToDef->NI ? ndi-ToDef->NI : ndi);

                        // Skip if no mask
                        if (!ToDef->Mask || ToDef->Mask[idxt]) {

                           *(ip++) = idxt;

                           // If the previous value is nodata, initialize the counter
                           if (!DEFVALID(ToDef, fld[idxt])) {
                              fld[idxt] = (Opt->Interp == IR_SUM || Opt->Interp == IR_AVERAGE  || Opt->Interp == IR_VECTOR_AVERAGE || Opt->Interp == IR_VARIANCE || Opt->Interp == IR_SQUARE || Opt->Interp == IR_NORMALIZED_COUNT || Opt->Interp == IR_COUNT) ? 0.0 : (Opt->Interp == IR_MAXIMUM ? -HUGE_VAL : HUGE_VAL);
                              if (aux) aux[idxt] = fld[idxt];
                           }
                                                
                           switch(Opt->Interp) {
                              case IR_MAXIMUM          : if (vx > fld[idxt]) fld[idxt] = vx;
                                                         break;
                              case IR_MINIMUM          : if (vx < fld[idxt]) fld[idxt] = vx;
                                                         break;
                              case IR_SUM              : fld[idxt] += vx;
                                                         break;
                              case IR_VARIANCE         : acc[idxt]++;
                                                         val = Opt->Ancilliary[idxt];
                                                         fld[idxt] += (vx-val)*(vx-val);
                                                         break;
                              case IR_SQUARE           : acc[idxt]++;
                                                         fld[idxt] += vx*vx;
                                                         break;
                              case IR_COUNT            : acc[idxt]++;
                              case IR_AVERAGE          :
                              case IR_NORMALIZED_COUNT : if (Opt->Table) {
                                                            t = 0;
                                                            while(t < ToDef->NK) {
                                                               if (vx == Opt->Table[t]) {
                                                                  if (Opt->Interp != IR_COUNT) acc[idxt]++;
                                                                  fld[t*nij+idxt] += 1.0;
                                                                  break;
                                                               }
                                                               t++;
                                                            }
                                                         } else if (Opt->lutDef) {
                                                            int32_t hasvalue = 0;
                                                            t = 0;
                                                            while(t < nbClass) {
                                                               if (vx == fromClass[t]) {
                                                                  if (toClass) { // toClass is between 1 and 26
                                                                     if (toClass[t] > 0)
                                                                        {
                                                                        fld[(toClass[t]-1)*nij+idxt] += 1.0;
                                                                        hasvalue = 1;
                                                                        }
                                                                  } else {
                                                                     for (i = 1; i < Opt->lutDim ; i++) { // skip 1st column : the class id
                                                                        val = Opt->lutDef[i][t+1];
                                                                        if (val > 0) {
                                                                           fld[(rpnClass[i]-1)*nij+idxt] += val;
                                                                           hasvalue = 1;
                                                                        }
                                                                     }
                                                                  }
                                                                  break;
                                                               }
                                                               t++;
                                                            }
                                                            // Dont count missing values
                                                            if (hasvalue)
                                                               if (Opt->Interp != IR_COUNT) acc[idxt]++;
                                                         } else {
                                                            if (DEFVALID(FromDef, vx)) {
                                                               fld[idxt] += vx;
                                                               if (Opt->Interp != IR_COUNT) acc[idxt]++;
                                                            }
                                                         }
                                                         break;
                              case IR_VECTOR_AVERAGE   : vx = DEG2RAD(vx); fld[idxt] += cos(vx); aux[idxt] += sin(vx); acc[idxt]++; break;
                              default:
                                 Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid interpolation type\n", __func__);
                                 return FALSE;
                           }
                        }
                     }
                  }
               }
               *(ip++) = REF_INDEX_SEPARATOR;
            }
            *(ip++) = REF_INDEX_END;
            gset->IndexSize = (ip-gset->Index)+1;

            GeoScan_Clear(&gscan);
         }
         if (fromClass) free(fromClass);
         if (toClass) free(toClass);
         if (rpnClass) free(rpnClass);
      }
   }

   return TRUE;
}

/*----------------------------------------------------------------------------
 * @brief  Finalize multi-step interpolations
 * @date   May 2006
    //! [in]     ToRef        Destination GeoRef pointer
    //! [in]     ToDef        Destination data definition
    //! [in]     Opt          Interpolation options
 *
 *    @return                 Number of grid points involved
*/
int64_t GeoRef_InterpFinalize(TGeoRef *ToRef, TDef *ToDef,TGeoOptions *Opt) {

   int64_t k,n,nij,idx;
   double  val;

   if (Opt->Interp < IR_NORMALIZED_CONSERVATIVE) {
      return(TRUE);
   }

   nij=FSIZE2D(ToDef);
   idx=0;

   switch(Opt->Interp) {

      case IR_NORMALIZED_CONSERVATIVE:
         for (k = 0; k < ToDef->NK; k++) {
            for(n = 0; n < nij; n++, idx++) {
               if (ToDef->Buffer[idx] != 0.0) {
                  Def_Get(ToDef, 0, idx, val);
                  val /= ToDef->Buffer[idx];
                  Def_Set(ToDef, 0, idx, val);
                  ToDef->Buffer[idx] = 0.0;
               }
            }
         }
         break;

      case IR_AVERAGE:
      case IR_VARIANCE:
      case IR_SQUARE:
      case IR_NORMALIZED_COUNT:
      case IR_COUNT:
      case IR_VECTOR_AVERAGE:
      case IR_ACCUM:
      case IR_BUFFER:
         for(k = 0; k < ToDef->NK; k++) {
            for(n = 0; n < nij; n++, idx++) {
               val = ToDef->NoData;
               switch(Opt->Interp) {
                  case IR_ACCUM:
                     if (ToDef->Accum) val = ToDef->Accum[n];
                     break;

                  case IR_BUFFER:
                     if (ToDef->Buffer) val = ToDef->Buffer[idx];
                     break;

                  default:
                     if (ToDef->Buffer) {
                        if (DEFVALID(ToDef, ToDef->Buffer[idx])) {
                           if (ToDef->Aux) {
                              if (ToDef->Accum && ToDef->Accum[n] != 0.0) {
                                 ToDef->Buffer[idx] /= ToDef->Accum[n];
                                 ToDef->Aux[idx] /= ToDef->Accum[n];
                              }
                              val = RAD2DEG(atan2(ToDef->Aux[idx], ToDef->Buffer[idx]));
                              if (val < 0.0) val += 360.0;
                           } else {
                              val = ToDef->Buffer[idx];
                              if (ToDef->Accum && ToDef->Accum[n] != 0.0) {
                                 val = ToDef->Buffer[idx]/ToDef->Accum[n];
                              } else {
                                 val = ToDef->Buffer[idx];
                              }
                           }
                        }
                     }
               }
               Def_Set(ToDef, 0, idx, val);
            }
         }
         break;
   }

   // Copy first column to last if it's repeated
   if (ToRef->Type&GRID_REPEAT) {
      for(k = 0; k < ToDef->NK; k++) {
         idx = k*nij;
         for(n = ToRef->Y0; n <= ToRef->Y1; n++, idx += ToDef->NI) {
            Def_Get(ToDef, 0, idx, val);
            Def_Set(ToDef, 0, idx+ToDef->NI-1, val);
         }
      }
   }

   return(nij);
}

TGeoRef *GeoRef_DefSubSelect(TGeoRef *Ref,TDef *Def,int N) {

   int i;

   if (!Ref || !Def) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid referential of definition\n", __func__);
   }
   if (N<-1 || N>=Ref->NbSub) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid sub grid index (%i)\n", __func__,N);
   }

   // If the subgrid index is different from the current
   if (N!=Ref->Sub && N<Ref->NbSub) {

      for (i=0; i<Ref->NbSet; i++) {
         if (Ref->Sets[i].Index) {
            free(Ref->Sets[i].Index);
            Ref->Sets[i].Index=NULL;
            Ref->Sets[i].IndexMethod=IR_UNDEF;
         }
      }

      // Define sub index position within Def data arrays
      Def->Idx=0;
      if (N>0) {
         Def->Idx=Ref->Subs[N-1]->NX*Ref->Subs[N-1]->NY;
      }

      // Return selected georef object
      Ref->Sub=N;
      return(N>=0?Ref->Subs[N]:Ref);
   }
   return(NULL);
}

/*----------------------------------------------------------------------------
 * @brief  Initiates interpolations from TDef types
 * @date   August 2024
    //! [in]     ToRef        Destination GeoRef pointer
    //! [in]     ToDef        Destination data definition
    //! [in]     FromRef      Source GeoRef pointer
    //! [in]     FromDef      Source data definition
 *
 *    @return        FALSE on error, TRUE otherwise
*/
int32_t GeoRef_InterpDef(TGeoRef *ToRef, TDef *ToDef, TGeoRef *FromRef, TDef *FromDef, TGeoOptions *Opt, int32_t Final) {

   void    *pf0, *pt0, *pf1, *pt1;
   int32_t  u, k, code = FALSE;

   if (!Opt) Opt=&GeoRef_Options;

   if (!ToRef || !ToDef) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid destination\n", __func__);
      return FALSE;
   }
   if (!FromRef || !FromDef) {
      Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid source\n", __func__);
      return FALSE;
   }

   // Same GeoRef
   if (Opt->Interp < IR_NOP && GeoRef_Equal(ToRef, FromRef)) {
      Def_CopyData(ToDef, FromDef);
      return TRUE;
   }

   switch(Opt->Interp) {
      case IR_UNDEF:
         break;
      case IR_NEAREST:
      case IR_LINEAR:
      case IR_CUBIC:
         // Loop on vertical levels
         for(k = 0; k < ToDef->NK; k++) {
            Def_Pointer(ToDef, 0, k*FSIZE2D(ToDef), pt0);
            Def_Pointer(FromDef, 0, k*FSIZE2D(FromDef), pf0);

            if (ToDef->Data[1]) {
               // Interpolation vectorielle
               Def_Pointer(ToDef, 1, k*FSIZE2D(ToDef), pt1);
               Def_Pointer(FromDef, 1, k*FSIZE2D(FromDef), pf1);

               // In case of Y grid, get the speed and dir instead of wind components
               // since grid oriented components dont mean much
               if (ToRef->GRTYP[0] == 'Y') {
                  code = GeoRef_InterpWD(ToRef, FromRef, Opt, pt0, pt1, pf0, pf1);
               } else {
                  code = GeoRef_InterpUV(ToRef, FromRef, Opt, pt0, pt1, pf0, pf1);
               }
            } else{
               // Interpolation scalaire
               code = GeoRef_Interp(ToRef, FromRef, Opt, pt0, pf0);
            }
            // Interpolate mask if need be
            if (FromDef->Mask && ToDef->Mask) {
               code = GeoRef_InterpMask(ToRef, FromRef, Opt, &ToDef->Mask[k*FSIZE2D(ToDef)], &FromDef->Mask[k*FSIZE2D(FromDef)]);
            }
         }
         break;

      case IR_CONSERVATIVE:
      case IR_NORMALIZED_CONSERVATIVE:
         // Check compatibility between source and destination
//         if (!Def_Compat(field0->Def, field1->Def)) {
//            field0->GRef = GeoRef_Find(GeoRef_Resize(field0->GRef, field0->Def->NI, field0->Def->NJ));
//         }
         code = GeoRef_InterpConservative(ToRef, ToDef, FromRef, FromDef, Opt);
         break;

      case IR_MAXIMUM:
      case IR_MINIMUM:
      case IR_SUM:
      case IR_AVERAGE:
      case IR_VARIANCE:
      case IR_SQUARE:
      case IR_NORMALIZED_COUNT:
      case IR_COUNT:
      case IR_VECTOR_AVERAGE:
      case IR_NOP:
      case IR_ACCUM:
      case IR_BUFFER:
         // In case of U grid, we have to select the right subgrid
         if (FromRef->GRTYP[0] == 'U') {
            for(u=0;u<FromRef->NbSub;u++) {
               GeoRef_DefSubSelect(FromRef,FromDef,u);
               code = GeoRef_InterpAverage(ToRef, ToDef, FromRef->Subs[u], FromDef, Opt);
            }
            GeoRef_DefSubSelect(FromRef,FromDef,-1);
         } else {
            code = GeoRef_InterpAverage(ToRef, ToDef, FromRef, FromDef, Opt);
         }
         break;

      case IR_SUBNEAREST:
      case IR_SUBLINEAR:
         // Is the internal buffer the right size ?
         if (ToDef->Sub && ToDef->SubSample != Opt->Sampling) {
            free(ToDef->Sub); ToDef->Sub = NULL;
         }
         ToDef->SubSample = Opt->Sampling;
         code = GeoRef_InterpSub(ToRef, ToDef, FromRef, FromDef, Opt);
         break;
   }

   if (Final) {
      code = GeoRef_InterpFinalize(ToRef,ToDef, Opt);
   }
   return code;
}

/**----------------------------------------------------------------------------
 * @brief  Re-initialise reprojection buffer
 * @date   February 2008
    //! [in]  Scan     Reprojection buffer
*/
void GeoScan_Clear(TGeoScan *Scan) {

   if (Scan) {
      if (Scan->X) free(Scan->X);
      if (Scan->Y) free(Scan->Y);
      if (Scan->V) free(Scan->V);
      if (Scan->D) free(Scan->D);

      Scan->X = Scan->Y = NULL;
      Scan->V = NULL;
      Scan->D = NULL;
      Scan->N = Scan->S = Scan->DX = Scan->DY = 0;
   }
}

/**----------------------------------------------------------------------------
 * @brief  Initialiser la structure App
 * @date   Janvier 2017
    //! [in]  Type     App type (APP_MASTER = single independent process, APP_THREAD = threaded co-process)
    //! [in]  Name     Application name
    //! [in]  Version  Application version
    //! [in]  Desc     Application description
    //! [in]  Stamp    TimeStamp
 *
 *    @return              Parametres de l'application initialisee
*/void GeoScan_Init(TGeoScan *Scan) {

   if (Scan) {
      Scan->X = Scan->Y = NULL;
      Scan->V = NULL;
      Scan->D = NULL;
      Scan->N = Scan->S = Scan->DX = Scan->DY = 0;
   }
}

/**----------------------------------------------------------------------------
 * @brief  Reproject a stream of coordinates and extracte values
 * @date   February 2008
    //! [in]  Scan     Reprojection buffer
    //! [in]  ToRef    Destination geo reference
    //! [in]  ToDef    Destination data definition
    //! [in]  FromRef  Source geo reference
    //! [in]  FromDef  Destination data definition
    //! [in]  X0       Lower x limit
    //! [in]  Y0       Lower y limiy
    //! [in]  X1       Higher x limit
    //! [in]  Y1       Higher y limit
    //! [in]  Dim      Grid cell dimension (1 = point, 2 = area)
    //! [in]  Degree   Interpolation degree
 *
 *    @return              Size of results
*/
int32_t _GeoScan_Get(TGeoScan *Scan, TGeoRef *ToRef, TDef *ToDef, TGeoRef *FromRef, TDef *FromDef, TGeoOptions *Opt, int32_t X0, int32_t Y0, int32_t X1, int32_t Y1, int32_t Dim) {

   register int32_t idx, x, y, n = 0;
   int32_t          d = 0, sz, dd;
   double       x0, y0, v;

   if (!Opt) Opt=&GeoRef_Options;

   if (!Scan || !ToRef || !FromRef) {
      return 0;
   }

   // Check limits
   X0 = fmax(X0, FromRef->X0);
   Y0 = fmax(Y0, FromRef->Y0);
   X1 = fmin(X1, FromRef->X1);
   Y1 = fmin(Y1, FromRef->Y1);

   // Adjust scan buffer sizes
   Scan->DX = X1-X0+1;
   Scan->DY = Y1-Y0+1;
   dd = (Scan->DX+1)*(Scan->DY+1);
   sz = Scan->DX*Scan->DY;

   if (Scan->S < sz) {
      if (!(Scan->X = (double*)realloc(Scan->X, dd*sizeof(double))))
         return 0;
      if (!(Scan->Y = (double*)realloc(Scan->Y, dd*sizeof(double))))
         return 0;
      if (!(Scan->V = (unsigned int*)realloc(Scan->V, sz*sizeof(unsigned int))))
         return 0;
      if (!(Scan->D = (float*)realloc(Scan->D, sz*sizeof(float))))
         return 0;
      Scan->S = sz;
   }

   dd = Dim-1;
   Scan->N = 0;

   for(y = Y0; y <= Y1+dd; y++) {
      idx = (y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
      for(x = X0; x <= X1+dd; x++, idx++, n++) {
         if (x <= X1 && y <= Y1) {
            Scan->V[Scan->N++] = idx;
         }

         x0 = dd ? x-0.5 : x;
         y0 = dd ? y-0.5 : y;

         GeoRef_XY2LL(FromRef, &x0, &y0, &Scan->X[n], &Scan->Y[n], 1, FALSE);


         if (FromRef->Transform) {
            Scan->X[n] = FromRef->Transform[0]+FromRef->Transform[1]*x0+FromRef->Transform[2]*y0;
            Scan->Y[n] = FromRef->Transform[3]+FromRef->Transform[4]*x0+FromRef->Transform[5]*y0;
         } else {
            Scan->X[n] = x0;
            Scan->Y[n] = y0;
         }
      }
   }

   // WKT grid type
   if (FromRef->GRTYP[0] == 'W') {
#ifdef HAVE_GDAL
      for(y = Y0; y <= Y1+dd; y++) {
         idx = (y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
         for(x = X0; x <= X1+dd; x++, idx++, n++) {
            if (x <= X1 && y <= Y1) {
               Scan->V[Scan->N++] = idx;
            }

            x0 = dd ? x-0.5 : x;
            y0 = dd ? y-0.5 : y;
            if (FromRef->Transform) {
               Scan->X[n] = FromRef->Transform[0]+FromRef->Transform[1]*x0+FromRef->Transform[2]*y0;
               Scan->Y[n] = FromRef->Transform[3]+FromRef->Transform[4]*x0+FromRef->Transform[5]*y0;
            } else {
               Scan->X[n] = x0;
               Scan->Y[n] = y0;
            }
         }
      }

      if (FromRef->Function) {
         OCTTransform(FromRef->Function, n, Scan->X, Scan->Y, NULL);
      }
#endif
      d = dd ? 2 : 1;

   // Y GRTYP type
   } else if (FromRef->GRTYP[0] == 'Y') {
      for(y = Y0; y <= Y1; y++) {
         idx = (y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
         for(x = X0; x <= X1; x++, idx++, n++) {
            if (x <= X1 && y <= Y1) {
               Scan->V[Scan->N++] = idx;
            }
            Scan->X[n] = FromRef->AX[idx];
            Scan->Y[n] = FromRef->AY[idx];
         }
      }
      d = 1;

   // Other RPN grids
   } else {
      for(y = Y0; y <= Y1+dd; y++) {
         idx = (y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
         for(x = X0; x <= X1+dd; x++, idx++, n++) {
            if (x <= X1 && y <= Y1) {
               Scan->V[Scan->N++] = idx;
            }
            Scan->X[n] = dd ? x+0.5 : x+1.0;
            Scan->Y[n] = dd ? y+0.5 : y+1.0;
         }
      }
      GeoRef_XY2LL(FromRef, Scan->Y, Scan->X, Scan->X, Scan->Y, n, FALSE);

      d = dd ? 2 : 1;
   }

   // Project to destination grid
   for(x = n-1; x >= 0; x--) {
      x0 = Scan->X[x];
      y0 = Scan->Y[x];

      if (ToDef) {
         Scan->D[x] = ToDef->NoData;
      }

      // If we're inside
      if (GeoRef_LL2XY(ToRef, &Scan->X[x], &Scan->Y[x], &y0, &x0, 1, FALSE) && ToDef) {
         Def_GetValue(ToRef, ToDef, Opt, 0, Scan->X[x], Scan->Y[x], 0, &v, NULL);
         Scan->D[x] = v;
      }
   }
   return d;
}

int32_t GeoScan_Get(TGeoScan *Scan, TGeoRef *ToRef, TDef *ToDef, TGeoRef *FromRef, TDef *FromDef, TGeoOptions *Opt, int32_t X0, int32_t Y0, int32_t X1, int32_t Y1) {

   register int32_t idx, x, y, n = 0;
   int32_t          d = 0, sz, dd;
   int32_t          ix, iy;
#ifdef HAVE_GDAL
   double           x0, y0, v;
#endif
   if (!Opt) Opt=&GeoRef_Options;

   if (!Scan || !ToRef || !FromRef) {
      return 0;
   }

   // Check limits
   X0 = fmax(X0, FromRef->X0);
   Y0 = fmax(Y0, FromRef->Y0);
   X1 = fmin(X1, FromRef->X1);
   Y1 = fmin(Y1, FromRef->Y1);

   // Adjust scan buffer sizes
   Scan->DX = X1-X0+1;
   Scan->DY = Y1-Y0+1;
   dd = (Scan->DX+1)*(Scan->DY+1);
   sz = Scan->DX*Scan->DY;

   if (Scan->S < sz) {
      if (!(Scan->X = (double*)realloc(Scan->X, dd*sizeof(double))))
         return 0;
      if (!(Scan->Y = (double*)realloc(Scan->Y, dd*sizeof(double))))
         return 0;
      if (!(Scan->V = (unsigned int*)realloc(Scan->V, sz*sizeof(unsigned int))))
         return 0;
      if (!(Scan->D = (float*)realloc(Scan->D, sz*sizeof(float))))
         return 0;
      Scan->S = sz;
   }

   dd = FromDef->CellDim-1;
   Scan->N = 0;

   // WKT grid type
   if (FromRef->GRTYP[0] == 'W') {
#ifdef HAVE_GDAL
      for(y = Y0; y <= Y1+dd; y++) {
         idx = (y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
         for(x = X0; x <= X1+dd; x++, idx++, n++) {
            if (x <= X1 && y <= Y1) {
               Scan->V[Scan->N++] = idx;
            }

            x0 = dd ? x-0.5 : x;
            y0 = dd ? y-0.5 : y;
            if (FromRef->Transform) {
               Scan->X[n] = FromRef->Transform[0]+FromRef->Transform[1]*x0+FromRef->Transform[2]*y0;
               Scan->Y[n] = FromRef->Transform[3]+FromRef->Transform[4]*x0+FromRef->Transform[5]*y0;
            } else {
               Scan->X[n] = x0;
               Scan->Y[n] = y0;
            }
         }
      }

      if (FromRef->Function) {
         OCTTransform(FromRef->Function, n, Scan->X, Scan->Y, NULL);
      }
#endif
      d = dd ? 2 : 1;
      sz = 8;

   // Y Grid type
   } else if (FromRef->GRTYP[0] == 'Y') {
      for(y = Y0; y <= Y1; y++) {
         idx = (y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
         for(x = X0; x <= X1; x++, idx++, n++) {
            if (x <= X1 && y <= Y1) {
               Scan->V[Scan->N++] = idx;
            }
            Scan->X[n] = FromRef->AX[idx];
            Scan->Y[n] = FromRef->AY[idx];
         }
      }
      d = 1;
      
   // Other RPN grids
   } else {
      for(y = Y0; y <= Y1+dd; y++) {
         idx = (y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
         for(x = X0; x <= X1+dd; x++, idx++, n++) {
            if (x <= X1 && y <= Y1) {
               Scan->V[Scan->N++] = idx;
            }
            Scan->X[n] = dd ? x-0.5 : x;
            Scan->Y[n] = dd ? y-0.5 : y;
        }
      }
      GeoRef_XY2LL(FromRef, Scan->Y, Scan->X, Scan->X, Scan->Y, n, TRUE);

      d = dd ? 2 : 1;
   }

   // Project to destination grid
   if (ToRef->GRTYP[0] == 'W' || ToRef->GRTYP[0] == 'M') {
#ifdef HAVE_GDAL
      for(x = 0; x < n; x++) {
         x0 = Scan->X[x];
         y0 = Scan->Y[x];

         if (ToDef) {
            Scan->D[x] = ToDef->NoData;
         }

         if (GeoRef_LL2XY(ToRef, &Scan->X[x], &Scan->Y[x], &y0, &x0, 1, FALSE)) {
            if (ToDef) {
               Def_GetValue(ToRef, ToDef, Opt, 0, Scan->X[x], Scan->Y[x], 0, &v, NULL);
               Scan->D[x] = v;
            }
         }
      }

/*TODO
         if (ToRef->Function)
            OCTTransform(ToRef->InvFunction, n, Scan->X, Scan->Y, NULL);

         if (ToRef->InvTransform) {
            for(x = 0; x < n; x++) {
               x0 = ToRef->InvTransform[0]+ToRef->InvTransform[1]*Scan->X[x]+ToRef->InvTransform[2]*Scan->Y[x];
               y0 = ToRef->InvTransform[3]+ToRef->InvTransform[4]*Scan->X[x]+ToRef->InvTransform[5]*Scan->Y[x];
               Scan->X[x] = x0;
               Scan->Y[x] = y0;
            }
         }
*/
#endif
   } else {
      GeoRef_LL2XY(ToRef, Scan->X, Scan->Y, Scan->Y, Scan->X, n, FALSE);

      // If we have the data of source and they're float, get it's values right now
      if (ToDef && ToDef->Type == TD_Float32) {
         GeoRef_XYVal(ToRef, Opt, Scan->D, (float*)ToDef->Mode, Scan->X, Scan->Y, n);
      }

      for(x = 0; x < n; x++) {
        if (ToDef) {
            ix = lrint(Scan->X[x]);
            iy = lrint(Scan->Y[x]);
            idx = FIDX2D(ToDef, ix, iy);

            if (!FIN2D(ToDef, ix, iy) || (ToDef->Mask && !ToDef->Mask[idx])) {
               // If we're outside, set to nodata
               Scan->D[x] = ToDef->NoData;
            } else if (ToDef->Type < TD_Float32) {
               // Otherwise, set nearest data if not floats
               Def_GetMod(ToDef, idx, Scan->D[x]);
            }
         }
      }
   }
   return d;
}
