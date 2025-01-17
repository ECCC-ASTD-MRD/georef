#include <App.h>
#include "GeoRef.h"

int32_t GeoRef_InterpMask(TGeoRef *RefTo, TGeoRef *RefFrom,TGeoOptions *Opt,char *MaskOut,char *MaskIn) {
 
   TGeoSet *gset=NULL;
   if (RefTo->NbSub > 0 || RefFrom->NbSub > 0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(FALSE);
   }

   if (!Opt) Opt=&RefTo->Options;
   if (!Opt) Opt=&GeoRef_Options;
   gset=GeoRef_SetGet(RefTo,RefFrom,Opt);

   if (RefFrom->GRTYP[0] == 'Y') {
      memcpy(MaskOut,gset->mask,RefTo->NX*RefTo->NY*sizeof(int));
   } else {
      f77name(qqq_ezsint_mask)(MaskOut,gset->X,gset->Y,&RefTo->NX,&RefTo->NY,MaskIn,&RefFrom->NX,&RefFrom->NY,Opt->Interp);
   }
   return(TRUE);
}

int32_t GeoRef_MaskZones(TGeoRef *RefTo,TGeoRef *RefFrom,int32_t *MaskOut,int32_t *MaskIn) {

   TGeoSet *gset=NULL;
   int32_t       npts_in, npts_out;
   
   if (RefTo->NbSub > 0 || RefFrom->NbSub > 0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   }

   gset=GeoRef_SetGet(RefTo,RefFrom,NULL);

   npts_in  = RefFrom->NX*RefFrom->NY;
   npts_out = RefTo->NX*RefTo->NY;

   f77name(qqq_ezget_mask_zones)(MaskOut,gset->X,gset->Y, &RefTo->NX, &RefTo->NY, MaskIn, &RefFrom->NX, &RefFrom->NY);
   return(0);
}

int32_t GeoRef_MaskYYDefine(TGeoRef *Ref) {

   int32_t ni,nj,yni,ynj,i,j,k,ii;
   int32_t i0,i1,j0,j1;
   int32_t ig1ref,ig2ref,ig3ref,ig4ref;
  
   k=i0=i1=0;
   for (i=0; i < Ref->NX; i++) {
      if (Ref->AX[i] >= 45.0 && Ref->AX[i] <= 315.0) {
         k++;
         if (k == 1) i0=i;
         i1=i;
      }
   }

   yni=k;
   k=j0=j1=0;
   for (i=0; i < Ref->NY; i++) {
      if (Ref->AY[i] >= -45.0 && Ref->AY[i] <= 45.0) {
         k++;
         if (k == 1) j0=i;
         j1=i;
      }
   }
   ynj=k;
   Ref->mymaskgrid = GeoRef_Define(GRID_SUB,yni,ynj,Ref->GRTYP,Ref->RPNHeadExt.grref,Ref->RPNHeadExt.igref1,Ref->RPNHeadExt.igref2,Ref->RPNHeadExt.igref3,Ref->RPNHeadExt.igref4,&Ref->AX[i0],&Ref->AY[j0]);
   Ref->mymaskgridi0=i0;
   Ref->mymaskgridi1=i1;
   Ref->mymaskgridj0=j0;
   Ref->mymaskgridj1=j1;
 
  return(0);
}

int32_t GeoRef_MaskYYApply(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,int32_t ni,int32_t nj,float *maskout,double *dlat,double *dlon,double *yinlat,double *yinlon,int32_t *yyincount,double *yanlat,double *yanlon,int32_t *yyancount) {

   TGeoOptions opt;
   int32_t i,j,k;
   int32_t yincount,yancount;
   float *yin_fld;

   yin_fld = (float*)calloc(RefFrom->mymaskgrid->NX*RefFrom->mymaskgrid->NY,sizeof(float));
 
   // Get original options
   opt=*Opt;
   opt.Interp=IR_NEAREST;
   opt.Extrap=ER_VALUE;
   opt.NoData=1.0;
    
   GeoRef_Interp(RefTo,RefFrom->mymaskgrid,&opt,maskout,yin_fld);
 
   free(yin_fld);

   // Now create the destination grids
   yancount=0;
   yincount=0;
   for (k=0; k<ni*nj; k++) {
      if (maskout[k] == 1.0) {
         yanlat[yancount]=dlat[k];
         yanlon[yancount]=dlon[k];
         yancount++;
      } else {
         yinlat[yincount]=dlat[k];
         yinlon[yincount]=dlon[k];
         yincount++;
      }
   }
   *yyincount = yincount;
   *yyancount = yancount;

   return(1);
}
