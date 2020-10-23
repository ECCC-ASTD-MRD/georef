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

int c_gdsetmask(TGeoRef *gr, int *mask) {

   int ni, nj;
  
   if (gr->NbSub > 0) {
      App_Log(ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   }

   ni = gr->NX;
   nj = gr->NY;

   if (gr->mask != NULL){
      free(gr->mask);
   }
   gr->mask = (int *)malloc(ni*nj*sizeof(int));
   memcpy(gr->mask, mask, ni*nj*sizeof(int));
   return(0);
}

int c_gdgetmask(TGeoRef *gr, int *mask)
{
   int ni, nj;

   if (gr->NbSub > 0) {
      App_Log(ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   }

   ni = gr->NX;
   nj = gr->NY;

   if (gr->mask != NULL) {
      memcpy(mask, gr->mask, ni*nj*sizeof(int));
      return 0;
   } else {
      mask = NULL;
      return -1;
   }
}

int c_ezget_mask_zones(int *mask_out, int *mask_in, TGeoRef *gdout, TGeoRef *gdin) {

   TGridSet *gset=NULL;
   int       npts_in, npts_out;
   float    *x,*y;
   
   if (gdout->NbSub > 0 || gdin->NbSub > 0) {
      App_Log(ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   }

   gset=GeoRef_SetGet(gdout,gdin);

   npts_in  = gdin->NX*gdin->NY;
   npts_out = gdout->NX*gdout->NY;

    x = (float *) gset->x;
    y = (float *) gset->y;

   f77name(qqq_ezget_mask_zones)(mask_out, x, y, &gdout->NX, &gdout->NY, mask_in, &gdin->NX, &gdin->NY);
   return(0);
}

int c_ezgdef_yymask(TGeoRef *Ref) {

   int ni,nj,yni,ynj,i,j,k,ii;
   int i0,i1,j0,j1;
   int ig1ref,ig2ref,ig3ref,ig4ref;
  
   float *ax,*ay;

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
   Ref->mymaskgrid = GeoRef_RPNCreateInMemory(yni,ynj,Ref->GRTYP,Ref->RPNHead.GRREF,Ref->RPNHead.IGREF[X_IG1],Ref->RPNHead.IGREF[X_IG2],Ref->RPNHead.IGREF[X_IG3],Ref->RPNHead.IGREF[X_IG4],&Ref->AX[i0],&Ref->AY[j0]);
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
 