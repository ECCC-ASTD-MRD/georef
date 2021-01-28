/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "App.h"
#include "GeoRef.h"

int GeoRef_InterpMask(TGeoRef *RefTo, TGeoRef *RefFrom,char *MaskOut,char *MaskIn,TGridSet **GSet) {
 
   TGridSet *gset=NULL;
   if (RefTo->NbSub > 0 || RefFrom->NbSub > 0) {
      App_Log(ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(FALSE);
   }

   gset=GeoRef_SetGet(RefTo,RefFrom,GSet);

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
      App_Log(ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   }

   gset=GeoRef_SetGet(RefTo,RefFrom,NULL);

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
   Ref->mymaskgrid = GeoRef_CreateInMemory(yni,ynj,Ref->GRTYP,Ref->RPNHead.GRREF,Ref->RPNHead.IGREF[X_IG1],Ref->RPNHead.IGREF[X_IG2],Ref->RPNHead.IGREF[X_IG3],Ref->RPNHead.IGREF[X_IG4],&Ref->AX[i0],&Ref->AY[j0]);
   Ref->mymaskgridi0=i0;
   Ref->mymaskgridi1=i1;
   Ref->mymaskgridj0=j0;
   Ref->mymaskgridj1=j1;

   App_Log(DEBUG,"%s: Subgd.mymaskgrid   = %p\n",__func__,Ref->mymaskgrid);
   App_Log(DEBUG,"%s: Subgd.mymaskgridi0 = %d pt=%f\n",__func__,Ref->mymaskgridi0, Ref->AX[i0]);
   App_Log(DEBUG,"%s: Subgd.mymaskgridi1 = %d pt=%f\n",__func__,Ref->mymaskgridi1, Ref->AX[i1]);
   App_Log(DEBUG,"%s: Subgd.mymaskgridj0 = %d pt=%f\n",__func__,Ref->mymaskgridj0, Ref->AY[j0]);
   App_Log(DEBUG,"%s: Subgd.mymaskgridj1 = %d pt=%f\n",__func__,Ref->mymaskgridj1, Ref->AY[j1]);
 
  return(0);
}

int GeoRef_MaskYYApply(TGeoRef *RefTo,TGeoRef *RefFrom,int ni,int nj,float *maskout,double *dlat,double *dlon,double *yinlat,double *yinlon,int *yyincount,double *yanlat,double *yanlon,int *yyancount) {

   TGeoRef *yin_mg;
   TGeoOptions opt;
   int ivalue,icode,i,j,k;
   int yincount,yancount,yni,ynj;
   float *yin_fld, global_extrap_value, local_extrap_value;
   int interp_degree,extrap_degree;
   float extrap_value,local_val;
   char global_interp_degree[32],global_extrap_degree[32];
  
   yin_mg=RefFrom->mymaskgrid;
   yni=yin_mg->NX;
   ynj=yin_mg->NY;

   yin_fld = (float *) malloc(yni*ynj*sizeof(float));
   memset(yin_fld,0.0,yni*ynj*sizeof(float));
 
   // Get original options
   memcpy(&opt,&RefFrom->Options,sizeof(TGeoOptions));
 
   RefFrom->Options.ExtrapValue=1.0;
   RefFrom->Options.ExtrapDegree=ER_VALUE;
   GeoRef_Interp(RefTo,yin_mg,maskout,yin_fld,NULL);
   // Masking is done,reset original interp options
   memcpy(&RefFrom->Options,&opt,sizeof(TGeoOptions));
   free(yin_fld);

   // Now create the destination grids
   yancount=0;
   yincount=0;
   for (j=0; j<nj; j++) {
      for (i=0;i<ni; i++) {
         k=(j*ni)+i; 
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
   }
   *yyincount = yincount;
   *yyancount = yancount;

   return(icode);
}
