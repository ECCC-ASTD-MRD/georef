#include <App.h>
#include "GeoRef.h"

int GeoRef_InterpMask(TGeoRef *RefTo, TGeoRef *RefFrom,char *MaskOut,char *MaskIn) {
 
   TGridSet *gset=NULL;
   if (RefTo->NbSub > 0 || RefFrom->NbSub > 0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(FALSE);
   }

   gset=GeoRef_SetGet(RefTo,RefFrom);

   if (RefFrom->GRTYP[0] == 'Y') {
      memcpy(MaskOut,gset->mask,RefTo->NX*RefTo->NY*sizeof(int));
   } else {
      f77name(qqq_ezsint_mask)(MaskOut,gset->X,gset->Y,&RefTo->NX,&RefTo->NY,MaskIn,&RefFrom->NX,&RefFrom->NY,RefFrom->Options.InterpDegree);
   }
   return(TRUE);
}

int GeoRef_MaskZones(TGeoRef *RefTo,TGeoRef *RefFrom,int *MaskOut,int *MaskIn) {

   TGridSet *gset=NULL;
   int       npts_in, npts_out;
   
   if (RefTo->NbSub > 0 || RefFrom->NbSub > 0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   }

   gset=GeoRef_SetGet(RefTo,RefFrom);

   npts_in  = RefFrom->NX*RefFrom->NY;
   npts_out = RefTo->NX*RefTo->NY;

   f77name(qqq_ezget_mask_zones)(MaskOut,gset->X,gset->Y, &RefTo->NX, &RefTo->NY, MaskIn, &RefFrom->NX, &RefFrom->NY);
   return(0);
}

int GeoRef_MaskYYDefine(TGeoRef *Ref) {

   int ni,nj,yni,ynj,i,j,k,ii;
   int i0,i1,j0,j1;
   int ig1ref,ig2ref,ig3ref,ig4ref;
  
   k=0;
   for (i=0; i < Ref->NX; i++) {
      if (Ref->AX[i] >= 45.0 && Ref->AX[i] <= 315.0) {
         k++;
         if (k == 1) i0=i;
         i1=i;
      }
   }

   yni=k;
   k=0;
   for (i=0; i < Ref->NY; i++) {
      if (Ref->AY[i] >= -45.0 && Ref->AY[i] <= 45.0) {
         k++;
         if (k == 1) j0=i;
         j1=i;
      }
   }
   ynj=k;
   Ref->mymaskgrid = GeoRef_Define(NULL,yni,ynj,Ref->GRTYP,Ref->RPNHeadExt.grref,Ref->RPNHeadExt.igref1,Ref->RPNHeadExt.igref2,Ref->RPNHeadExt.igref3,Ref->RPNHeadExt.igref4,&Ref->AX[i0],&Ref->AY[j0]);
   Ref->mymaskgridi0=i0;
   Ref->mymaskgridi1=i1;
   Ref->mymaskgridj0=j0;
   Ref->mymaskgridj1=j1;

   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Subgd.mymaskgrid   = %p\n",__func__,Ref->mymaskgrid);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Subgd.mymaskgridi0 = %d pt=%f\n",__func__,Ref->mymaskgridi0, Ref->AX[i0]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Subgd.mymaskgridi1 = %d pt=%f\n",__func__,Ref->mymaskgridi1, Ref->AX[i1]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Subgd.mymaskgridj0 = %d pt=%f\n",__func__,Ref->mymaskgridj0, Ref->AY[j0]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Subgd.mymaskgridj1 = %d pt=%f\n",__func__,Ref->mymaskgridj1, Ref->AY[j1]);
 
  return(0);
}

int GeoRef_MaskYYApply(TGeoRef *RefTo,TGeoRef *RefFrom,int ni,int nj,float *maskout,double *dlat,double *dlon,double *yinlat,double *yinlon,int *yyincount,double *yanlat,double *yanlon,int *yyancount) {

   TGeoOptions opt;
   int i,j,k;
   int yincount,yancount;
   float *yin_fld;
  
   yin_fld = (float*)calloc(RefFrom->mymaskgrid->NX*RefFrom->mymaskgrid->NY,sizeof(float));
 
   // Get original options
   memcpy(&opt,&RefFrom->Options,sizeof(TGeoOptions));
 
   RefFrom->Options.ExtrapValue=1.0;
   RefFrom->Options.ExtrapDegree=ER_VALUE;
   GeoRef_Interp(RefTo,RefFrom->mymaskgrid,maskout,yin_fld);
 
   // Masking is done,reset original interp options
   memcpy(&RefFrom->Options,&opt,sizeof(TGeoOptions));
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
