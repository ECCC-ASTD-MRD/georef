#include <App.h>
#include "GeoRef.h"

/*----------------------------------------------------------------------------
 * @brief  Define a referential of type RPN U
 * @date   Avril 2015
 *    @param[in] NI      Dimension en X
 *    @param[in] NJ      Dimension en Y
 *    @param[in] grref   Reference grid type ('E', 'G', 'L', 'N', 'S')
 *    @param[in] VerCode Version of the grid    
 *    @param[in] NbSub   Number of sub grids
 *    @param[in] Subs    Array of pointers to subgrid
 *
 *    @return            Georef (NULL=error)
*/
TGeoRef* GeoRef_CreateU(int32_t NI,int32_t NJ,char *grref,int32_t VerCode,int32_t NbSub,TGeoRef **Subs) {

   int32_t  i;
   TGeoRef *ref,*fref,*sub_gd;
    
   if (NbSub <= 1) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: NbSub given is less than 2\n",__func__);
      return(NULL);
   }
   if (VerCode != 1) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid VerCode\n",__func__);
      return(NULL);
   }

   ref = GeoRef_New();
  
   if (VerCode == 1) {
      sub_gd = Subs[0];

      ref->RPNHead.grtyp[0]=ref->GRTYP[0] = 'U';
      ref->RPNHeadExt.grref[0] = grref[0];
      ref->NX       = NI;
      ref->NY       = NJ;
         
      // To add more uniqueness to the super-grid index Yin-Yang grid, we also add the rotation of YIN
      ref->RPNHead.ig1  = sub_gd->RPNHead.ig1;
      ref->RPNHead.ig2  = sub_gd->RPNHead.ig2;
      ref->RPNHead.ig3  = sub_gd->RPNHead.ig3;
      ref->RPNHead.ig4  = sub_gd->RPNHead.ig4;
      ref->RPNHeadExt.igref1=VerCode;
      ref->RPNHeadExt.igref2=0;
      ref->RPNHeadExt.igref3=0;
      ref->RPNHeadExt.igref4=0;
      ref->NbSub= NbSub;

      //TODO: Create AXY record
   }
  
   // This georef already exists
   if (fref=GeoRef_Find(ref)) {
      free(ref);
      GeoRef_Incr(fref);
      return(fref);
   }

   // This is a new georef
   GeoRef_Add(ref);

   ref->Subs = (TGeoRef **)malloc(NbSub*sizeof(TGeoRef*));

   for (i=0; i < NbSub; i++) {
      ref->Subs[i] = Subs[i];
      GeoRef_MaskYYDefine(Subs[i]);
      Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Grille[%p].Subs[%p] has maskgrid=%p\n",__func__,ref,Subs[i],sub_gd->mymaskgrid);
   }

   GeoRef_Qualify(ref);

   return(ref);
}
