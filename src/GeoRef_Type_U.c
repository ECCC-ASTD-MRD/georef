#include <App.h>
#include "GeoRef.h"


//! \file


//! Define a referential of type RPN U
TGeoRef * GeoRef_CreateU(
    //! [in] X dimension
    const int32_t NI,
    //! [in] Y dimension
    const int32_t NJ,
    //! [in] Reference grid type ('E', 'G', 'L', 'N', 'S')
    const char * const grref,
    //! [in] Grid version
    const int32_t VerCode,
    //! [in] Number of sub grids
    const int32_t NbSub,
    //! [in] Array of sub grid pointers
    TGeoRef **Subs
) {
    //! \returns Georeference, NULL on error
    if (NbSub <= 1) {
        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: NbSub given is less than 2\n", __func__);
        return NULL;
    }
    if (VerCode != 1) {
        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid VerCode\n", __func__);
        return NULL;
    }

    TGeoRef * const ref = GeoRef_New();

    if (VerCode == 1) {
        TGeoRef * const sub_gd = Subs[0];

        ref->RPNHead.grtyp[0] = ref->GRTYP[0] = 'U';
        ref->RPNHeadExt.grref[0] = grref[0];
        ref->NX = NI;
        ref->NY = NJ;

        // To add more uniqueness to the super-grid index Yin-Yang grid, we also add the rotation of YIN
        ref->RPNHead.ig1 = sub_gd->RPNHead.ig1;
        ref->RPNHead.ig2 = sub_gd->RPNHead.ig2;
        ref->RPNHead.ig3 = sub_gd->RPNHead.ig3;
        ref->RPNHead.ig4 = sub_gd->RPNHead.ig4;
        ref->RPNHeadExt.igref1 = VerCode;
        ref->RPNHeadExt.igref2 = 0;
        ref->RPNHeadExt.igref3 = 0;
        ref->RPNHeadExt.igref4 = 0;
        ref->NbSub = NbSub;

        //! \todo Create AXY record
    }

    TGeoRef * const fref = GeoRef_Find(ref);
    if (fref) {
        // This georef already exists
        free(ref);
        GeoRef_Incr(fref);
        return fref;
    }

    // This is a new georef
    GeoRef_Add(ref);

    ref->Subs = (TGeoRef **)malloc(NbSub *sizeof(TGeoRef*));

     for (int32_t i = 0; i < NbSub; i++) {
        ref->Subs[i] = Subs[i];
        GeoRef_MaskYYDefine(Subs[i]);
        Lib_Log(APP_LIBGEOREF, APP_DEBUG, "%s: Grille[%p].Subs[%p] has maskgrid=%p\n", __func__, ref, Subs[i], sub_gd->mymaskgrid);
    }

    GeoRef_Qualify(ref);
    return ref;
}
