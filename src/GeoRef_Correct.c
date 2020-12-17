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

int GeoRef_CorrValNorth(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout,float *zin) {

   TGridSet *gset=NULL;
   int i;
   int npts;
   float *temp, *vals;
   double *temp_y;
   float ay[4], poleval;
   int ni, nj, i1, i2, j1, j2;
   int un = 1;
   int quatre = 4;

   gset=GeoRef_SetGet(RefTo,RefFrom,NULL);

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
                     f77name(ez8_irgdint_3_wnnc)(vals,gset->zones[NORTH].x, gset->zones[NORTH].y,&npts,RefFrom->AX, ay, temp,&ni, &j1, &j2, &(RefFrom->Extension));
                  } else {
                     ay[0] = RefFrom->AY[RefFrom->j2-3];
                     ay[1] = RefFrom->AY[RefFrom->j2-2];
                     ay[2] = RefFrom->AY[RefFrom->j2-1];
                     ay[3] = 90.0;
                     f77name(ez8_irgdint_3_wnnc)(vals,gset->zones[NORTH].x,gset->zones[NORTH].y,&npts,RefFrom->AX, ay, temp,&ni, &j1, &j2, &(RefFrom->Extension));
                  }
	                break;
	             default:
	                f77name(ez8_rgdint_3_wnnc)(vals,gset->zones[NORTH].x,gset->zones[NORTH].y,&npts,temp,&ni, &j1, &j2, &(RefFrom->Extension));
	                break;
	          }
   	        break;

         case IR_LINEAR:
	          temp_y = (double*) malloc(npts*sizeof(double));
	          for (i=0; i < npts; i++) {
               temp_y[i] = gset->zones[NORTH].y[i] - (1.0 * (RefFrom->j2-3));
	          }
	          f77name(ez8_rgdint_1_w)(vals,gset->zones[NORTH].x,temp_y,&npts,temp,&RefFrom->Options.NoData,&ni, &un, &quatre, &(RefFrom->Extension));
	          free(temp_y);
	          break;

         case IR_NEAREST:
	          temp_y = (double*) malloc(npts*sizeof(double));
	          for (i=0; i < npts; i++) {
	             temp_y[i] = gset->zones[NORTH].y[i] - (1.0 * (RefFrom->j2-3));
	          }
	          f77name(ez8_rgdint_0)(vals,gset->zones[NORTH].x,temp_y,&npts,temp,&ni, &un, &quatre);
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

int GeoRef_CorrValSouth(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout,float *zin) {

   TGridSet *gset=NULL;
   int i;
   int npts;
   float vpolesud;
   float *temp, *temp_y, *vals;
   float ay[4];
   int ni, nj, i1, i2, j1, j2;
   int un = 1;
   int quatre = 4;
  
   gset=GeoRef_SetGet(RefTo,RefFrom,NULL);

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
                     f77name(ez8_irgdint_3_wnnc)(vals,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,RefFrom->AX,ay,temp,&ni,&j1,&j2,&RefFrom->Extension);
                  } else {
                     ay[0] = -90.0;
                     ay[1] = RefFrom->AY[0];
                     ay[2] = RefFrom->AY[1];
                     ay[3] = RefFrom->AY[2];
                     f77name(ez8_irgdint_3_wnnc)(vals,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,RefFrom->AX,ay,temp,&ni,&j1,&j2,&RefFrom->Extension);
                  }
	                break;

	             default:
	                f77name(ez8_rgdint_3_wnnc)(vals,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,temp,&ni, &j1, &j2, &RefFrom->Extension);
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
   	        f77name(ez8_rgdint_1_w)(vals,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,temp,&RefFrom->Options.NoData,&ni, &j1, &j2,&RefFrom->Extension);
	          free(temp_y);
	          break;

         case IR_NEAREST:
  	        temp_y = (float *) malloc(npts*sizeof(float));
/*	   for (i=0; i < npts; i++)
	     {
	     temp_y[i] = gset->zones[SOUTH].y[i] - (1.0*j1);
	     }*/
	        f77name(ez8_rgdint_0)(vals,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,temp,&ni, &j1, &j2);
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

int GeoRef_CorrectValue(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout, float *zin) {

   TGridSet *gset=NULL;
   int i,ierc;
   float valmax, valmin,fudgeval;
   int fudgeval_set;
   int degIntCourant;
   int npts,nj;
   float vpolnor, vpolsud;
   float *temp;

   fudgeval_set = 0;
  
   gset=GeoRef_SetGet(RefTo,RefFrom,NULL);

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
      f77name(ez_aminmax)(&valmin,&valmax,zin,&(RefFrom->NX),&nj);
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
      GeoRef_CorrValNorth(RefTo,RefFrom,zout,zin);
   }

   if (gset->zones[SOUTH].npts > 0) {
      GeoRef_CorrValSouth(RefTo,RefFrom,zout,zin);
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

int GeoRef_CalcPolarWindNorth(TGeoRef *Ref,float *polar_uu_in,float *polar_vv_in,float *uuin,float *vvin,int ni,int nj) {
  
   int k1, k2;
   float *polar_wd, *polar_spd, *polar_uu, *polar_vv;
   double *polar_lat,*polar_lon,*polar_lat_gem, *polar_lon_gem, *polar_x, *polar_y;
   char grtypn[2],grtypa[2];
   float xlat1, xlat2, xlon1, xlon2;
   int ig1n, ig2n, ig3n, ig4n;
   float pi, pj, d60, dgrw;
   int i,j,ier;
   TGeoRef *gda, *gdps;
   float uupole, vvpole;
   double quatrevingtdix, zero;

   polar_uu  = (float *) malloc(ni*sizeof(float));
   polar_vv  = (float *) malloc(ni*sizeof(float));
   polar_wd  = (float *) malloc(ni*sizeof(float));
   polar_spd = (float *) malloc(ni*sizeof(float));
   polar_lat = (double*)malloc(4*ni*sizeof(double));
   polar_lon = &polar_lat[ni];
   polar_x   = &polar_lat[ni*2];
   polar_y   = &polar_lat[ni*3];

   for (i=0; i < ni; i++) {
      polar_x[i] = 1.0 * (i+1);
      polar_y[i] = 1.0 * nj;
   }
  
   GeoRef_XY2LL(Ref,polar_lat,polar_lon,polar_x,polar_y, ni,TRUE);

   if (Ref->GRTYP[0] == 'Z' && Ref->RPNHead.GRREF[0] == 'E') {
      polar_lat_gem   = (double*) malloc(2*ni*sizeof(double));
      polar_lon_gem   = &polar_lat_gem[ni];
    
      for (i=0; i < ni; i++) {
         polar_lat_gem[i] = polar_lat[i];
         polar_lon_gem[i] = polar_lon[i];
      }
    
     f77name(cigaxg)(Ref->RPNHead.GRREF, &xlat1, &xlon1, &xlat2, &xlon2, &Ref->RPNHead.IGREF[X_IG1], &Ref->RPNHead.IGREF[X_IG2], &Ref->RPNHead.IGREF[X_IG3], &Ref->RPNHead.IGREF[X_IG4]);
     GeoRef_RotateXY(polar_lat_gem,polar_lon_gem,polar_lon,polar_lat,ni,xlat1,xlon1,xlat2,xlon2);
   }

   grtypa[0] = 'A';
   gda = GeoRef_Create(24,12, grtypa, 0,0,0,0,0);
   GeoRef_UV2WD(gda, polar_spd, polar_wd,  &uuin[(nj-1)*ni], &vvin[(nj-1)*ni], polar_lat, polar_lon, ni);
  
   pi   = 0.0;
   pj   = 0.0;
   d60  = 1000.0;
   dgrw = 0.0;
   grtypn[0] = 'N';
   f77name(cxgaig)(grtypn, &ig1n, &ig2n, &ig3n, &ig4n, &pi, &pj, &d60, &dgrw);
   gdps = GeoRef_Create(ni, 1, grtypn, ig1n, ig2n, ig3n, ig4n, 0);
   GeoRef_WD2UV(gdps, polar_uu, polar_vv, polar_spd,  polar_wd, polar_lat, polar_lon, ni);

   f77name(ez_calcpoleval)(&uupole, polar_uu, &ni, Ref->AX, &Ref->GRTYP, &Ref->RPNHead.GRREF,1,1);
   f77name(ez_calcpoleval)(&vvpole, polar_vv, &ni, Ref->AX, &Ref->GRTYP, &Ref->RPNHead.GRREF,1,1);

   quatrevingtdix = 90.0;
   zero = 0.0;
   GeoRef_UV2WD(gdps, polar_spd, polar_wd,  &uupole, &vvpole, &quatrevingtdix, &zero, 1);
 

   polar_lat[0] = 90.0;
   for (i=1; i < ni; i++) {
     polar_wd[i]  = polar_wd[0] + polar_lon[i];
     polar_spd[i] = polar_spd[0];
     polar_lat[i] = 90.0;
   }
   polar_wd[0] = polar_wd[0] + polar_lon[0];
  
   GeoRef_WD2UV(gda, polar_uu, polar_vv, polar_spd,  polar_wd, polar_lat, polar_lon, ni);
  
   for (j=0; j < 3; j++) {
     for (i=0; i < ni; i++) {
        k1 = j * ni + i;
        k2 = (nj-3+j)  * ni + i;
        polar_uu_in[k1] = uuin[k2];
        polar_vv_in[k1] = vvin[k2];
      }
   }
  
   for (i=0; i < ni; i++) {
     k1 = 3 * ni+i;
     polar_uu_in[k1] = polar_uu[i];
     polar_vv_in[k1] = polar_vv[i];
   }

   free(polar_lat);

   free(polar_spd);
   free(polar_wd);
   free(polar_vv);
   free(polar_uu);

   if (Ref->GRTYP[0] == 'Z' && Ref->RPNHead.GRREF[0] == 'E' && polar_lat_gem != NULL)  {
     free(polar_lat_gem);
   }

   ier = GeoRef_Free(gdps);
   return(0);
}

int GeoRef_CalcPolarWindSouth(TGeoRef *Ref,float *polar_uu_in,float *polar_vv_in,float *uuin,float *vvin,int ni,int nj) {

   int k1, k2;
   float *polar_wd, *polar_spd,*polar_uu,*polar_vv;
   double *polar_lat,*polar_lon,*polar_lat_gem, *polar_lon_gem, *polar_x, *polar_y;
   char grtyps[2],grtypa[2];
   float xlat1, xlat2, xlon1, xlon2;
   int ig1n, ig2n, ig3n, ig4n;
   float pi, pj, d60, dgrw;
   int i,j,ier;
   TGeoRef *gda, *gdps;
   float uupole, vvpole;
   double quatrevingtdix, zero;

   polar_uu  = (float *) malloc(ni*sizeof(float));
   polar_vv  = (float *) malloc(ni*sizeof(float));
   polar_wd  = (float *) malloc(ni*sizeof(float));
   polar_spd = (float *) malloc(ni*sizeof(float));
   polar_lat = (double*)malloc(4*ni*sizeof(double));
   polar_lon = &polar_lat[ni];
   polar_x   = &polar_lat[ni*2];
   polar_y   = &polar_lat[ni*3];

   for (i=0; i < ni; i++) {
      polar_x[i] = 1.0 * (i+1);
      polar_y[i] = 1.0;
   }
  
   GeoRef_XY2LL(Ref,polar_lat,polar_lon,polar_x,polar_y,ni,TRUE);

   if (Ref->GRTYP[0] == 'Z' && Ref->RPNHead.GRREF[0] == 'E') {
      polar_lat_gem   = (double*) malloc(2*ni*sizeof(double));
      polar_lon_gem   = &polar_lat_gem[ni];
    
     for (i=0; i < ni; i++) {
        polar_lat_gem[i] = polar_lat[i];
        polar_lon_gem[i] = polar_lon[i];
     }
    
     f77name(cigaxg)(Ref->RPNHead.GRREF, &xlat1, &xlon1, &xlat2, &xlon2, &Ref->RPNHead.IGREF[X_IG1], &Ref->RPNHead.IGREF[X_IG2], &Ref->RPNHead.IGREF[X_IG3], &Ref->RPNHead.IGREF[X_IG4]);
     GeoRef_RotateXY(polar_lat_gem,polar_lon_gem,polar_lon,polar_lat,ni,xlat1,xlon1,xlat2,xlon2);
   }

   grtypa[0] = 'A';
   gda = GeoRef_Create(24,12, grtypa, 0,0,0,0,0);
   GeoRef_UV2WD(gda, polar_spd, polar_wd,  uuin, vvin, polar_lat, polar_lon, ni);
  
   pi   = 0.0;
   pj   = 0.0;
   d60  = 1000.0;
   dgrw = 0.0;
   grtyps[0] = 'S';
   f77name(cxgaig)(grtyps, &ig1n, &ig2n, &ig3n, &ig4n, &pi, &pj, &d60, &dgrw);
   gdps = GeoRef_Create(ni, 1, grtyps, ig1n, ig2n, ig3n, ig4n, 0);
   GeoRef_WD2UV(gdps, polar_uu, polar_vv, polar_spd,  polar_wd, polar_lat, polar_lon, ni);

   f77name(ez_calcpoleval)(&uupole,polar_uu,&ni,Ref->AX,&Ref->GRTYP,&Ref->RPNHead.GRREF,1,1);
   f77name(ez_calcpoleval)(&vvpole,polar_vv,&ni,Ref->AX,&Ref->GRTYP,&Ref->RPNHead.GRREF,1,1);

   quatrevingtdix = -90.0;
   zero = 0.0;
   GeoRef_UV2WD(gdps, polar_spd, polar_wd,  &uupole, &vvpole, &quatrevingtdix, &zero, 1);
  
   polar_lat[0] = -90.0;
   for (i=1; i < ni; i++) {
     polar_wd[i]  = polar_wd[0] - polar_lon[i];
     polar_spd[i] = polar_spd[0];
     polar_lat[i] = -90.0;
   }
   polar_wd[0] = polar_wd[0] + polar_lon[0];
  
   GeoRef_WD2UV(gda, polar_uu, polar_vv, polar_spd,  polar_wd, polar_lat, polar_lon, ni);
  
   for (j=0; j < 3; j++) {
     for (i=0; i < ni; i++) {
        k1 = j * ni + i;
        k2 = (j+1)  * ni + i;
        polar_uu_in[k2] = uuin[k1];
        polar_vv_in[k2] = vvin[k1];
     }
   }
  
   for (i=0; i < ni; i++) {
     polar_uu_in[i] = polar_uu[i];
     polar_vv_in[i] = polar_vv[i];
   }

   free(polar_y);
   free(polar_x);
   free(polar_lat);
   free(polar_lon);
   free(polar_spd);
   free(polar_wd);
   free(polar_vv);
   free(polar_uu);

   if (Ref->GRTYP[0] == 'Z' && Ref->RPNHead.GRREF[0] == 'E' && polar_lat_gem != NULL)  {
     free(polar_lat_gem);
   }

   ier = GeoRef_Free(gdps);
   return(0);
}

int GeoRef_CorrVecNorth(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin) {

   TGridSet *gset=NULL;
   float     uupole, vvpole;
   float    *polar_uu_in, *polar_vv_in, *polar_uu_out, *polar_vv_out, *corr_uus, *corr_vvs,ay[4];
   double   *temp_y;
   int       ni, nj, i1, i2, j1, j2, degree,npts,i;
   int       quatre = 4;
   int       un = 1;

   gset=GeoRef_SetGet(RefTo,RefFrom,NULL);

   ni = RefFrom->NX;
   nj = RefFrom->j2 - RefFrom->j1 + 1;

   i1 = 1;
   i2 = ni;
   j1 = RefFrom->j2-2;
   j2 = j1 + 3;
   degree = 3;

   npts = gset->zones[NORTH].npts;

   polar_uu_in = (float *) malloc(4 * ni * sizeof(float));
   polar_vv_in = (float *) malloc(4 * ni * sizeof(float));
   corr_uus = (float *) malloc(npts * sizeof(float));
   corr_vvs = (float *) malloc(npts * sizeof(float));

   GeoRef_CalcPolarWindNorth(RefFrom,polar_uu_in,polar_vv_in,uuin,vvin,ni,nj);

   switch (RefFrom->Options.InterpDegree) {
      case IR_CUBIC:
         switch (RefFrom->GRTYP[0]) {
	          case 'Z':
	          case 'E':
	          case 'G':
               ay[0] = RefFrom->AY[RefFrom->j2-3];
               ay[1] = RefFrom->AY[RefFrom->j2-2];
               ay[2] = RefFrom->AY[RefFrom->j2-1];
               ay[3] = 90.0;
	             f77name(ez8_irgdint_3_wnnc)(corr_uus,gset->zones[NORTH].x,gset->zones[NORTH].y,&npts,RefFrom->AX,ay,polar_uu_in,&ni,&j1,&j2,&RefFrom->Extension);
	             f77name(ez8_irgdint_3_wnnc)(corr_vvs,gset->zones[NORTH].x,gset->zones[NORTH].y,&npts,RefFrom->AX,ay,polar_vv_in,&ni,&j1,&j2,&RefFrom->Extension);
	             break;

	          default:
	             f77name(ez8_rgdint_3_wnnc)(corr_uus,gset->zones[NORTH].x,gset->zones[NORTH].y,&npts,polar_uu_in,&ni,&j1,&j2,&RefFrom->Extension);
	             f77name(ez8_rgdint_3_wnnc)(corr_vvs,gset->zones[NORTH].x,gset->zones[NORTH].y,&npts,polar_vv_in,&ni,&j1,&j2,&RefFrom->Extension);
	             break;
	      }

      case IR_LINEAR:
         temp_y = (double*) malloc(npts*sizeof(double));
         for (i=0; i < npts; i++) {
            temp_y[i] = gset->zones[NORTH].y[i] - (1.0 * (RefFrom->j2-3));
	       }
         f77name(ez8_rgdint_1_w)(corr_uus,gset->zones[NORTH].x,temp_y,&npts,polar_uu_in,&RefFrom->Options.NoData,&ni, &un, &quatre, &RefFrom->Extension);
         f77name(ez8_rgdint_1_w)(corr_vvs,gset->zones[NORTH].x,temp_y,&npts,polar_vv_in,&RefFrom->Options.NoData,&ni, &un, &quatre, &RefFrom->Extension);
         free(temp_y);
         
         break;

      case IR_NEAREST:
         temp_y = (double*) malloc(npts*sizeof(double));
         for (i=0; i < npts; i++) {
            temp_y[i] = gset->zones[NORTH].y[i] - (1.0 * (RefFrom->j2-3));
	       }
         f77name(ez8_rgdint_0)(corr_uus,gset->zones[NORTH].x,temp_y,&npts,polar_uu_in,&ni, &un, &quatre);
         f77name(ez8_rgdint_0)(corr_vvs,gset->zones[NORTH].x,temp_y,&npts,polar_vv_in,&ni, &un, &quatre);
         free(temp_y);
         break;
   }


   for (i=0; i < gset->zones[NORTH].npts; i++) {
      uuout[gset->zones[NORTH].idx[i]] = corr_uus[i];
      vvout[gset->zones[NORTH].idx[i]] = corr_vvs[i];
   }

   free(polar_uu_in);
   free(polar_vv_in);
   free(corr_uus);
   free(corr_vvs);

   return(0);
}

int GeoRef_CorrVecSouth(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin) {
  
   TGridSet *gset=NULL;
   float *polar_uu_in, *polar_vv_in, *corr_uus, *corr_vvs, *temp_y, ay[4];
   int ni, nj, i1, i2, j1, j2, degree,npts,i;
   int idx_gdin;

   gset=GeoRef_SetGet(RefTo,RefFrom,NULL);

   npts = gset->zones[SOUTH].npts;
   ni = RefFrom->NX;
   nj = RefFrom->j2 - RefFrom->j1 + 1;

   i1 = 1;
   i2 = ni;
   j1 = RefFrom->j1 - 1;
   j2 = j1 + 3;
   degree = 3;

   polar_uu_in = (float *) malloc(4 * ni * sizeof(float));
   polar_vv_in = (float *) malloc(4 * ni * sizeof(float));
   corr_uus = (float *) malloc(npts * sizeof(float));
   corr_vvs = (float *) malloc(npts * sizeof(float));

   GeoRef_CalcPolarWindSouth(RefFrom,polar_uu_in,polar_vv_in,uuin,vvin,ni,nj);

   switch (RefFrom->Options.InterpDegree) {
      case IR_CUBIC:
         switch (RefFrom->GRTYP[0]) {
	          case 'Z':
	          case 'E':
	          case 'G':
               ay[0] = -90.;
               ay[1] = RefFrom->AY[0];
               ay[2] = RefFrom->AY[1];
               ay[3] = RefFrom->AY[2];
               f77name(ez8_irgdint_3_wnnc)(corr_uus,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,RefFrom->AX,ay,polar_uu_in,&ni,&j1,&j2,&RefFrom->Extension);
               f77name(ez8_irgdint_3_wnnc)(corr_vvs,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,RefFrom->AX,ay,polar_vv_in,&ni,&j1,&j2,&RefFrom->Extension);
               break;

            default:
               f77name(ez8_rgdint_3_wnnc)(corr_uus,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,polar_uu_in,&ni,&j1,&j2,&RefFrom->Extension);
               f77name(ez8_rgdint_3_wnnc)(corr_vvs,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,polar_vv_in,&ni,&j1,&j2,&RefFrom->Extension);
               break;
            }
      break;

      case IR_LINEAR:
         f77name(ez8_rgdint_1_w)(corr_uus,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,polar_uu_in,&RefFrom->Options.NoData,&ni,&j1,&j2,&RefFrom->Extension);
         f77name(ez8_rgdint_1_w)(corr_vvs,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,polar_vv_in,&RefFrom->Options.NoData,&ni,&j1,&j2,&RefFrom->Extension);
         break;

      case IR_NEAREST:
         f77name(ez8_rgdint_0)(corr_uus,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,polar_uu_in,&ni,&j1,&j2);
         f77name(ez8_rgdint_0)(corr_vvs,gset->zones[SOUTH].x,gset->zones[SOUTH].y,&npts,polar_vv_in,&ni,&j1,&j2);
         break;
   }


   for (i=0; i < gset->zones[SOUTH].npts; i++) {
      uuout[gset->zones[SOUTH].idx[i]] = corr_uus[i];
      vvout[gset->zones[SOUTH].idx[i]] = corr_vvs[i];
   }

   free(polar_uu_in);
   free(polar_vv_in);
   free(corr_uus);
   free(corr_vvs);

   return(0);
}

int GeoRef_CorrectVector(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin) {
 
   TGridSet *gset=NULL;
   int ier;
  
   gset=GeoRef_SetGet(RefTo,RefFrom,NULL);
  
   if (gset->zones[NORTH].npts > 0) {
      ier = GeoRef_CorrVecNorth(RefTo,RefFrom,uuout,vvout,uuin,vvin);
   }
  
   if (gset->zones[SOUTH].npts > 0) {
     ier = GeoRef_CorrVecSouth(RefTo,RefFrom,uuout,vvout,uuin,vvin);
   }
  
   if (gset->zones[NORTH_POLE].npts > 0) {
     ier = GeoRef_CorrVecNorth(RefTo,RefFrom,uuout,vvout,uuin,vvin);
   }
  
   if (gset->zones[SOUTH_POLE].npts > 0) {
     ier = GeoRef_CorrVecSouth(RefTo,RefFrom,uuout,vvout,uuin,vvin);
   }
      
   return(0);
}
