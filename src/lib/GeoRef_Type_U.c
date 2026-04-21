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
   if ((fref=GeoRef_Find(ref))) {
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

/**----------------------------------------------------------------------------
 * @brief  Merge 2 Z grids into a U grid (YinYang)
 *    @param[in]  YinRef   Pointer to Yin reference
 *    @param[in]  YangRef  Pointer to Yang reference
 */
TGeoRef* GeoRef_CreateUFromZMerge(TGeoRef *YinRef,TGeoRef *YangRef) {

   // Adapted from GEM yyencode.F90
   TGeoRef *uref=GeoRef_New();

   int   niyy, ni,nj, sindx, sindx_yin;
   char  family_uencode_S = 'F';
   int   version_uencode  = 1;
   float xlat1,xlon1,xlat2,xlon2;

   // Sanity check
   if (YinRef->NX!=YangRef->NX || YinRef->NY!=YangRef->NY){
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Yin and Yang records don't have the samesize : %dx%d != %dx%d\n",__func__,
         YinRef->NX,YinRef->NY,YangRef->NX,YangRef->NY);
      return(NULL);
   }

   ni=YinRef->NX;
   nj=YangRef->NY;
   niyy=5+2*(10+ni+nj);

   if(!(uref->AXY = (double*)malloc(niyy * sizeof(double)))){
      Lib_Log(APP_LIBGEOREF,APP_SYSTEM, "%s: Cannot allocate buffer of size %lu\n",__func__,niyy);
      return(NULL);
   }

   uref->AXY[0] = (int)family_uencode_S;  // équivalent C de Fortran iachar
   uref->AXY[1] = version_uencode;
   uref->AXY[2] = 2; // 2 grids (Yin & Yang);
   uref->AXY[3] = 1; // the 2 grids have same resolution;
   uref->AXY[4] = 1; // the 2 grids have same area extension on the sphere;

   // YIN
   sindx = 5;
   f77name(cigaxg)("E", &xlat1, &xlon1, &xlat2, &xlon2,
         &YinRef->RPNHeadExt.igref1,&YinRef->RPNHeadExt.igref2,&YinRef->RPNHeadExt.igref3,&YinRef->RPNHeadExt.igref4,1);
 
   // xlat1 must be greater than zero for yin grid
   if (xlat1 < 0.0){
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Yan and Yin Georef might be inverted\n",__func__);
      return(NULL);
   }

   uref->AXY[sindx  ] = ni;
   uref->AXY[sindx+1] = nj;
   uref->AXY[sindx+6] = xlat1;
   uref->AXY[sindx+7] = xlon1;
   uref->AXY[sindx+8] = xlat2;
   uref->AXY[sindx+9] = xlon2;
   memcpy(&uref->AXY[sindx+10],    YinRef->AX, ni * sizeof(double));
   memcpy(&uref->AXY[sindx+10+ni], YinRef->AY, nj * sizeof(double));    
   uref->AXY[sindx+2] = uref->AXY[sindx+10];
   uref->AXY[sindx+3] = uref->AXY[sindx+ 9+ni];
   uref->AXY[sindx+4] = uref->AXY[sindx+10+ni];
   uref->AXY[sindx+5] = uref->AXY[sindx+ 9+ni+nj];
   sindx_yin = sindx;
 
   // YANG
   f77name(cigaxg)("E", &xlat1, &xlon1, &xlat2, &xlon2,
            &YangRef->RPNHeadExt.igref1,&YangRef->RPNHeadExt.igref2,&YangRef->RPNHeadExt.igref3,&YangRef->RPNHeadExt.igref4,1);
   sindx              = sindx+10+ni+nj;
   uref->AXY[sindx  ] = ni;
   uref->AXY[sindx+1] = nj;
   uref->AXY[sindx+2] = uref->AXY[sindx_yin+10      ];
   uref->AXY[sindx+3] = uref->AXY[sindx_yin+ 9+ni];
   uref->AXY[sindx+4] = uref->AXY[sindx_yin+10+ni];
   uref->AXY[sindx+5] = uref->AXY[sindx_yin+ 9+ni+nj];
   uref->AXY[sindx+6] = xlat1;
   uref->AXY[sindx+7] = xlon1;
   uref->AXY[sindx+8] = xlat2;
   uref->AXY[sindx+9] = xlon2;

   // Note in the yyencode.F90 they copy data from yy but data from tictac is identical
   memcpy(&uref->AXY[sindx+10],    YangRef->AX, ni * sizeof(double));
   memcpy(&uref->AXY[sindx+10+ni], YangRef->AY, nj * sizeof(double));    

   uref->GRTYP[0]='U';
   uref->GRTYP[1]='\0';
   uref->NbSub=2;
   uref->NX=ni;
   uref->NY=nj*uref->NbSub;
   uref->Subs = (TGeoRef**)malloc(uref->NbSub*sizeof(TGeoRef*));
   uref->Subs[0]=YinRef;
   uref->Subs[1]=YangRef;
   uref->RPNHeadExt.igref1=1;
   uref->RPNHeadExt.igref2=0;
   uref->RPNHeadExt.igref3=0;
   uref->RPNHeadExt.igref4=0;
   uref->RPNHeadExt.grref[0]=family_uencode_S;
   uref->RPNHeadExt.grref[1]='\0';
   strncpy(uref->RPNHead.grtyp,uref->GRTYP,FST_GTYP_LEN);

   // Calculate unique ip1 ip2 ip3
   GeoRef_RPNHash(uref,NULL,NULL,NULL);

   return(uref);
}