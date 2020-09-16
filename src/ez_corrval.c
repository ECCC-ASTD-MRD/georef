/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "App.h"
#include "GeoRef.h"

int ez_corrval_aunord(float *zout, float *zin, TGeoRef *RefFrom, TGeoRef *RefTo) {

   TGridSet *gset=NULL;
   int i;
   int npts;
   float *temp, *vals, *temp_y;
   float ay[4], poleval;
   int ni, nj, i1, i2, j1, j2;
   int un = 1;
   int quatre = 4;

   gset=GeoRef_SetGet(RefTo,RefFrom);

   npts = gset->zones[NORTH].npts;

   if (npts > 0) {
      ni = RefFrom->NX;
      nj = RefFrom->j2 - RefFrom->j1 + 1;

      i1 = 1;
      i2 = ni;

      j1 = RefFrom->j2-2;
      j2 = j1 + 3;

      temp = (float *) malloc(4 * ni * sizeof(float));
      vals = (float *) malloc(npts * sizeof(float));
      f77name(ez_calcpoleval)(&poleval, &zin[(nj-1)*ni], &ni, RefFrom->AX, &(RefFrom->GRTYP),&(RefFrom->RPNHead.GRREF),1,1);
      f77name(ez_fillnpole)(temp, zin, &ni, &(RefFrom->j1), &(RefFrom->j2), &poleval);

      switch (RefFrom->Options.InterpDegree) {
         case IR_CUBIC:
	          switch (RefFrom->GRTYP[0]) {
	             case 'Z':
	             case 'E':
	             case 'G':
                  if  (RefFrom->AY[RefFrom->j2-1] == 90.0) {
                     ay[0] = RefFrom->AY[RefFrom->j2-4];
                     ay[1] = RefFrom->AY[RefFrom->j2-3];
                     ay[2] = RefFrom->AY[RefFrom->j2-2];
                     ay[3] = RefFrom->AY[RefFrom->j2-1];
                     f77name(ez_irgdint_3_wnnc)(vals,gset->zones[NORTH].x, gset->zones[NORTH].y,&npts,RefFrom->AX, ay, temp,&ni, &j1, &j2, &(RefFrom->Extension));
                  } else {
                     ay[0] = RefFrom->AY[RefFrom->j2-3];
                     ay[1] = RefFrom->AY[RefFrom->j2-2];
                     ay[2] = RefFrom->AY[RefFrom->j2-1];
                     ay[3] = 90.0;
                     f77name(ez_irgdint_3_wnnc)(vals,gset->zones[NORTH].x,gset->zones[NORTH].y,&npts,RefFrom->AX, ay, temp,&ni, &j1, &j2, &(RefFrom->Extension));
                  }
	                break;
	             default:
	                f77name(ez_rgdint_3_wnnc)(vals,gset->zones[NORTH].x,gset->zones[NORTH].y,&npts,temp,&ni, &j1, &j2, &(RefFrom->Extension));
	                break;
	          }
   	        break;

         case IR_LINEAR:
	          temp_y = (float *) malloc(npts*sizeof(float));
	          for (i=0; i < npts; i++) {
               temp_y[i] = gset->zones[NORTH].y[i] - (1.0 * (RefFrom->j2-3));
	          }
	          f77name(ez_rgdint_1_w)(vals,gset->zones[NORTH].x,temp_y,&npts,temp,&ni, &un, &quatre, &(RefFrom->Extension));
	          free(temp_y);
	          break;

         case IR_NEAREST:
	          temp_y = (float *) malloc(npts*sizeof(float));
	          for (i=0; i < npts; i++) {
	             temp_y[i] = gset->zones[NORTH].y[i] - (1.0 * (RefFrom->j2-3));
	          }
	          f77name(ez_rgdint_0)(vals,gset->zones[NORTH].x,temp_y,&npts,temp,&ni, &un, &quatre);
	          free(temp_y);
	          break;
      }

      for (i=0; i < gset->zones[NORTH].npts; i++) {
         zout[gset->zones[NORTH].idx[i]] = vals[i];
      }

      free(vals);
      free(temp);
   }
   return(0);
}

int ez_corrval_ausud(float *zout, float *zin, TGeoRef *RefFrom, TGeoRef *RefTo) {

   TGridSet *gset=NULL;
   int i;
   int npts;
   float vpolesud;
   float *temp, *temp_y, *vals;
   float ay[4];
   int ni, nj, i1, i2, j1, j2;
   int un = 1;
   int quatre = 4;
  
   gset=GeoRef_SetGet(RefTo,RefFrom);

   npts = gset->zones[SOUTH].npts;
   if (npts > 0) {
      ni = RefFrom->NX;

      i1 = 1;
      i2 = ni;
      j1 = RefFrom->j1 - 1;
      j2 = j1 + 3;

      temp = (float *) malloc(4 * ni * sizeof(float));
      vals = (float *) malloc(npts * sizeof(float));
      f77name(ez_calcpoleval)(&vpolesud, zin, &ni, RefFrom->AX,&RefFrom->GRTYP, &RefFrom->RPNHead.GRREF,1,1);
      f77name(ez_fillspole)(temp, zin, &ni, &RefFrom->j1, &RefFrom->j2, &vpolesud);

      switch (RefFrom->Options.InterpDegree) {
         case IR_CUBIC:
	          switch (RefFrom->GRTYP[0]) {
               case 'Z':
               case 'E':
               case 'G':
                  if  (RefFrom->AY[RefFrom->j1-1] == -90.0) {
                     ay[0] = RefFrom->AY[0];
                     ay[1] = RefFrom->AY[1];
                     ay[2] = RefFrom->AY[2];
                     ay[3] = RefFrom->AY[3];
                     f77name(ez_irgdint_3_wnnc)(vals,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,RefFrom->AX,ay,temp,&ni,&j1,&j2,&RefFrom->Extension);
                  } else {
                     ay[0] = -90.0;
                     ay[1] = RefFrom->AY[0];
                     ay[2] = RefFrom->AY[1];
                     ay[3] = RefFrom->AY[2];
                     f77name(ez_irgdint_3_wnnc)(vals,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,RefFrom->AX,ay,temp,&ni,&j1,&j2,&RefFrom->Extension);
                  }
	                break;

	             default:
	                f77name(ez_rgdint_3_wnnc)(vals,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,temp,&ni, &j1, &j2, &RefFrom->Extension);
	                break;
	          }
	          break;

         case IR_LINEAR:
	          temp_y = (float *) malloc(npts*sizeof(float));
/*	   for (i=0; i < npts; i++)
	     {
	     temp_y[i] = gset->zones[SOUTH].y[i] - (1.0*j1);
	     }
	   f77name(ez_rgdint_1_nw)(vals,gset->zones[SOUTH].x,temp_y,&npts,temp,&ni, &un, &quatre);*/
   	        f77name(ez_rgdint_1_w)(vals,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,temp,&ni, &j1, &j2,&RefFrom->Extension);
	          free(temp_y);
	          break;

         case IR_NEAREST:
  	        temp_y = (float *) malloc(npts*sizeof(float));
/*	   for (i=0; i < npts; i++)
	     {
	     temp_y[i] = gset->zones[SOUTH].y[i] - (1.0*j1);
	     }*/
	        f77name(ez_rgdint_0)(vals,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,temp,&ni, &j1, &j2);
	        free(temp_y);
	        break;
      }

      for (i=0; i < gset->zones[SOUTH].npts; i++) {
         zout[gset->zones[SOUTH].idx[i]] = vals[i];
      }

      free(vals);
      free(temp);
   }
   return(0);
}

int ez_corrval(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout, float *zin) {

   TGridSet *gset=NULL;
   int i,ierc;
   float valmax, valmin,fudgeval;
   int fudgeval_set;
   int degIntCourant;
   int npts,nj;
   float vpolnor, vpolsud;
   float *temp;

   extern float f77name(amax)();
   extern float f77name(amin)();

   fudgeval_set = 0;
  
   gset=GeoRef_SetGet(RefTo,RefFrom);

   nj = RefFrom->j2 - RefFrom->j1 +1;
   ierc = 0; /* no extrapolation */
  
   if (gset->zones[OUTSIDE].npts > 0) {
      ierc = 2; /* code to indicate extrapolation */
      if (RefFrom->Options.ExtrapDegree == ER_ABORT) {
         App_Log(ERROR,"%s: Points on destination lie outside source. Aborting as requested\n",__func__);
         return(-1);
      }

//TODO: check if ok
//    valmin=f77name(AMIN)(zin,&(RefFrom->NX),&nj, 0)
//    valmax=f77name(AMAX)(zin,&(RefFrom->NX),&nj, 0)
    f77name(ez_aminmax)(&valmin,&valmax,zin,&(RefFrom->NX), &nj);
      if (RefFrom->Options.ExtrapDegree >= ER_MAXIMUM) {
         if (RefFrom->Options.VectorMode) {
	          fudgeval = 0.0;
            fudgeval_set = 1;
	       } else {
            switch (RefFrom->Options.ExtrapDegree) {
               case ER_MAXIMUM:
                  fudgeval = valmax + 0.05 * (valmax - valmin);
                  fudgeval_set = 1;
                  App_Log(DEBUG,"%s: Maximum: %f \n",__func__,fudgeval);
                  break;

               case ER_MINIMUM:
                  fudgeval = valmin - 0.05 * (valmax - valmin);
                  fudgeval_set = 1;
                  App_Log(DEBUG,"%s: Minimum: %f \n",__func__,fudgeval);
                  break;

               case ER_VALUE:
                  fudgeval = RefFrom->Options.ExtrapValue;
                  fudgeval_set = 1;
                  App_Log(DEBUG,"%s: Value: %f \n",__func__,fudgeval);
                  break;
            }
         }

         if (fudgeval_set == 0) {
            App_Log(DEBUG,"%s: Fudgeval not set\n",__func__);
         }
         for (i=0; i < gset->zones[OUTSIDE].npts; i++) {
	          zout[gset->zones[OUTSIDE].idx[i]] = fudgeval;
	       }
      } else {
         degIntCourant = RefFrom->Options.InterpDegree;
         RefFrom->Options.InterpDegree = RefFrom->Options.ExtrapDegree;
         temp = (float *) malloc(gset->zones[OUTSIDE].npts*sizeof(float));

         GeoRef_InterpFinally(RefTo,RefFrom,temp,zin,gset->zones[OUTSIDE].x,gset->zones[OUTSIDE].y,gset->zones[OUTSIDE].npts);

         for (i=0; i < gset->zones[OUTSIDE].npts; i++) {
            zout[gset->zones[OUTSIDE].idx[i]] = temp[i];
         }
         free(temp);
         RefFrom->Options.InterpDegree = degIntCourant;
      }
   }

   if (RefFrom->Options.VectorMode) {
      return(ierc);
   }


   if (gset->zones[NORTH].npts > 0) {
      ez_corrval_aunord(zout,zin,RefFrom,RefTo);
   }

   if (gset->zones[SOUTH].npts > 0) {
      ez_corrval_ausud(zout,zin,RefFrom,RefTo);
   }


   if (gset->zones[NORTH_POLE].npts > 0 || gset->zones[SOUTH_POLE].npts > 0) {
      if (RefFrom->GRTYP[0] != 'w') {
         npts = RefFrom->NX * RefFrom->j2;
         f77name(ez_calcpoleval)(&vpolnor,&(zin[(nj-1)*RefFrom->NX]),&(RefFrom->NX),RefFrom->AX,RefFrom->GRTYP,RefFrom->RPNHead.GRREF,1,1);
         for (i=0; i < gset->zones[NORTH_POLE].npts; i++) {
	          zout[gset->zones[NORTH_POLE].idx[i]] = vpolnor;
	       }

         f77name(ez_calcpoleval)(&vpolsud,zin,&(RefFrom->NX),RefFrom->AX,RefFrom->GRTYP,RefFrom->RPNHead.GRREF,1,1);
         for (i=0; i < gset->zones[SOUTH_POLE].npts; i++) {
	          zout[gset->zones[SOUTH_POLE].idx[i]] = vpolsud;
	       }
      }
   }
   if ((RefFrom->GRTYP[0] == 'Z' || RefFrom->GRTYP[0] == '#') && RefFrom->RPNHead.GRREF[0] == 'E' && RefTo->GRTYP[0] == 'B') {
      f77name(ez_corrbgd)(zout, &(RefTo->NX), &(RefTo->NY), &(RefTo->RPNHead.IG[X_IG1]));
   }

   return(ierc);
}
