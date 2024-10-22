#include <App.h>
#include "GeoRef.h"

int32_t GeoRef_CorrValNorth(TGeoSet *Set,float *zout,float *zin) {

   TGeoRef *reffrom;
   int32_t i;
   int32_t npts;
   float *temp, *vals;
   double *temp_y;
   double ay[4];
   float poleval;
   int32_t ni, nj, i1, i2, j1, j2;
   int32_t un = 1;
   int32_t quatre = 4;

   if (!Set)
      return(0);
   
   reffrom=Set->RefFrom;

   npts = Set->zones[GRID_NORTH].npts;

   if (npts > 0) {
      ni = reffrom->NX;
      nj = reffrom->j2 - reffrom->j1 + 1;

      i1 = 1;
      i2 = ni;

      j1 = reffrom->j2-2;
      j2 = j1 + 3;

      temp = (float *) malloc(4 * ni * sizeof(float));
      vals = (float *) malloc(npts * sizeof(float));
      f77name(ez_calcpoleval)(&poleval,&zin[(nj-1)*ni],&ni,reffrom->AX,reffrom->GRTYP,reffrom->RPNHeadExt.grref,1,1);
      f77name(ez_fillnpole)(temp, zin, &ni, &(reffrom->j1), &(reffrom->j2), &poleval);

      switch (Set->Opt.Interp) {
         case IR_CUBIC:
	          switch (reffrom->GRTYP[0]) {
	             case 'Z':
	             case 'E':
	             case 'G':
                  if  (reffrom->AY[reffrom->j2-1] == 90.0) {
                     ay[0] = reffrom->AY[reffrom->j2-4];
                     ay[1] = reffrom->AY[reffrom->j2-3];
                     ay[2] = reffrom->AY[reffrom->j2-2];
                     ay[3] = reffrom->AY[reffrom->j2-1];
                     f77name(ez8_irgdint_3_wnnc)(vals,Set->zones[GRID_NORTH].x, Set->zones[GRID_NORTH].y,&npts,temp,&ni, &j1, &j2,reffrom->AX, ay, &(reffrom->Extension));
                  } else {
                     ay[0] = reffrom->AY[reffrom->j2-3];
                     ay[1] = reffrom->AY[reffrom->j2-2];
                     ay[2] = reffrom->AY[reffrom->j2-1];
                     ay[3] = 90.0;
                     f77name(ez8_irgdint_3_wnnc)(vals,Set->zones[GRID_NORTH].x,Set->zones[GRID_NORTH].y,&npts,temp,&ni, &j1, &j2,reffrom->AX, ay, &(reffrom->Extension));
                  }
	                break;
	             default:
	                f77name(ez8_rgdint_3)(vals,Set->zones[GRID_NORTH].x,Set->zones[GRID_NORTH].y,&npts,temp,&ni, &j1, &j2, &(reffrom->Extension),&Set->Opt.NoData);
	                break;
	          }
   	        break;

         case IR_LINEAR:
	          temp_y = (double*) malloc(npts*sizeof(double));
	          for (i=0; i < npts; i++) {
               temp_y[i] = Set->zones[GRID_NORTH].y[i] - (1.0 * (reffrom->j2-3));
	          }
	          f77name(ez8_rgdint_1)(vals,Set->zones[GRID_NORTH].x,temp_y,&npts,temp,&ni, &un, &quatre, &(reffrom->Extension),&Set->Opt.NoData);
	          free(temp_y);
	          break;

         case IR_NEAREST:
	          temp_y = (double*) malloc(npts*sizeof(double));
	          for (i=0; i < npts; i++) {
	             temp_y[i] = Set->zones[GRID_NORTH].y[i] - (1.0 * (reffrom->j2-3));
	          }
	          f77name(ez8_rgdint_0)(vals,Set->zones[GRID_NORTH].x,temp_y,&npts,temp,&ni, &un, &quatre,&Set->Opt.NoData);
	          free(temp_y);
	          break;
      }

      for (i=0; i < Set->zones[GRID_NORTH].npts; i++) {
         zout[Set->zones[GRID_NORTH].idx[i]] = vals[i];
      }

      free(vals);
      free(temp);
   }
   return(0);
}

int32_t GeoRef_CorrValSouth(TGeoSet *Set,float *zout,float *zin) {

   TGeoRef *reffrom=NULL;
   int32_t i;
   int32_t npts;
   float vpolesud;
   float *temp, *temp_y, *vals;
   double ay[4];
   int32_t ni, nj, i1, i2, j1, j2;
   int32_t un = 1;
   int32_t quatre = 4;
  
   if (!Set)
      return(0);
   
   reffrom=Set->RefFrom;

   npts = Set->zones[GRID_SOUTH].npts;
   if (npts > 0) {
      ni = reffrom->NX;

      i1 = 1;
      i2 = ni;
      j1 = reffrom->j1 - 1;
      j2 = j1 + 3;

      temp = (float *) malloc(4 * ni * sizeof(float));
      vals = (float *) malloc(npts * sizeof(float));
      f77name(ez_calcpoleval)(&vpolesud,zin,&ni,reffrom->AX,reffrom->GRTYP,reffrom->RPNHeadExt.grref,1,1);
      f77name(ez_fillspole)(temp, zin, &ni, &reffrom->j1, &reffrom->j2, &vpolesud);

      switch (Set->Opt.Interp) {
         case IR_CUBIC:
	          switch (reffrom->GRTYP[0]) {
               case 'Z':
               case 'E':
               case 'G':
                  if  (reffrom->AY[reffrom->j1-1] == -90.0) {
                     ay[0] = reffrom->AY[0];
                     ay[1] = reffrom->AY[1];
                     ay[2] = reffrom->AY[2];
                     ay[3] = reffrom->AY[3];
                     f77name(ez8_irgdint_3_wnnc)(vals,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,temp,&ni,&j1,&j2,reffrom->AX,ay,&reffrom->Extension);
                  } else {
                     ay[0] = -90.0;
                     ay[1] = reffrom->AY[0];
                     ay[2] = reffrom->AY[1];
                     ay[3] = reffrom->AY[2];
                     f77name(ez8_irgdint_3_wnnc)(vals,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,temp,&ni,&j1,&j2,reffrom->AX,ay,&reffrom->Extension);
                  }
	                break;

	             default:
	                f77name(ez8_rgdint_3)(vals,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,temp,&ni, &j1, &j2, &reffrom->Extension,&Set->Opt.NoData);
	                break;
	          }
	          break;

         case IR_LINEAR:
	          temp_y = (float *) malloc(npts*sizeof(float));
/*	   for (i=0; i < npts; i++)
	     {
	     temp_y[i] = Set->zones[GRID_SOUTH].y[i] - (1.0*j1);
	     }
	   f77name(ez_rgdint_1_nw)(vals,Set->zones[GRID_SOUTH].x,temp_y,&npts,temp,&ni, &un, &quatre);*/
   	        f77name(ez8_rgdint_1)(vals,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,temp,&ni, &j1, &j2,&reffrom->Extension,&Set->Opt.NoData);
	          free(temp_y);
	          break;

         case IR_NEAREST:
  	        temp_y = (float *) malloc(npts*sizeof(float));
/*	   for (i=0; i < npts; i++)
	     {
	     temp_y[i] = Set->zones[GRID_SOUTH].y[i] - (1.0*j1);
	     }*/
	        f77name(ez8_rgdint_0)(vals,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,temp,&ni, &j1, &j2,&Set->Opt.NoData);
	        free(temp_y);
	        break;
      }

      for (i=0; i < Set->zones[GRID_SOUTH].npts; i++) {
         zout[Set->zones[GRID_SOUTH].idx[i]] = vals[i];
      }

      free(vals);
      free(temp);
   }
   return(0);
}

int32_t GeoRef_CorrectValue(TGeoSet *Set,float *zout, float *zin) {

   TGeoRef *reffrom,*refto;
   int32_t i,ierc;
   float valmax, valmin,fudgeval;
   int32_t fudgeval_set;
   int32_t degIntCourant;
   int32_t npts,nj;
   float vpolnor, vpolsud;
   float *temp;

   fudgeval_set = 0;
  
   if (!Set)
      return(0);
   
   reffrom=Set->RefFrom;
   refto=Set->RefTo;

   nj = reffrom->j2 - reffrom->j1 +1;
   ierc = 0; /* no extrapolation */

   if (Set->zones[GRID_OUTSIDE].npts > 0) {
      ierc = 2; /* code to indicate extrapolation */
      if (Set->Opt.Extrap == ER_ABORT) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Points on destination lie outside source. Aborting as requested\n",__func__);
         return(-1);
      }

      f77name(ez_aminmax)(&valmin,&valmax,zin,&(reffrom->NX),&nj);
      if (Set->Opt.Extrap) {
         if (Set->Opt.VectorMode) {
	         fudgeval = 0.0;
            fudgeval_set = 1;
	       } else {
            switch (Set->Opt.CIndex) {
               case ER_MAXIMUM:
                  fudgeval = valmax + 0.05 * (valmax - valmin);
                  fudgeval_set = 1;
                  Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Maximum: %f \n",__func__,fudgeval);
                  break;

               case ER_MINIMUM:
                  fudgeval = valmin - 0.05 * (valmax - valmin);
                  fudgeval_set = 1;
                  Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Minimum: %f \n",__func__,fudgeval);
                  break;

               case ER_VALUE:
                  fudgeval = Set->Opt.ExtrapValue;
                  fudgeval_set = 1;
                  Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Value: %f \n",__func__,fudgeval);
                  break;
            }
         }

         if (fudgeval_set == 0) {
            Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Fudgeval not set\n",__func__);
         }
         for (i=0; i < Set->zones[GRID_OUTSIDE].npts; i++) {
	          zout[Set->zones[GRID_OUTSIDE].idx[i]] = fudgeval;
	       }
      } else {
         degIntCourant = Set->Opt.Interp;
         Set->Opt.Interp = Set->Opt.Extrap;
         temp = (float *) malloc(Set->zones[GRID_OUTSIDE].npts*sizeof(float));

         GeoRef_InterpFinally(refto,reffrom,&Set->Opt,temp,zin,Set->zones[GRID_OUTSIDE].x,Set->zones[GRID_OUTSIDE].y,Set->zones[GRID_OUTSIDE].npts,NULL);

         for (i=0; i < Set->zones[GRID_OUTSIDE].npts; i++) {
            zout[Set->zones[GRID_OUTSIDE].idx[i]] = temp[i];
         }
         free(temp);
         Set->Opt.Interp = degIntCourant;
      }
   }

   if (Set->Opt.VectorMode) {
      return(ierc);
   }

   if (Set->zones[GRID_NORTH].npts > 0) {
      GeoRef_CorrValNorth(Set,zout,zin);
   }

   if (Set->zones[GRID_SOUTH].npts > 0) {
      GeoRef_CorrValSouth(Set,zout,zin);
   }

   if (Set->zones[GRID_NORTH_POLE].npts > 0 || Set->zones[GRID_SOUTH_POLE].npts > 0) {
      if (reffrom->GRTYP[0] != 'w') {
         npts = reffrom->NX * reffrom->j2;
         f77name(ez_calcpoleval)(&vpolnor,&(zin[(nj-1)*reffrom->NX]),&(reffrom->NX),reffrom->AX,reffrom->GRTYP,reffrom->RPNHeadExt.grref,1,1);
         for (i=0; i < Set->zones[GRID_NORTH_POLE].npts; i++) {
	          zout[Set->zones[GRID_NORTH_POLE].idx[i]] = vpolnor;
	       }

         f77name(ez_calcpoleval)(&vpolsud,zin,&(reffrom->NX),reffrom->AX,reffrom->GRTYP,reffrom->RPNHeadExt.grref,1,1);
         for (i=0; i < Set->zones[GRID_SOUTH_POLE].npts; i++) {
	          zout[Set->zones[GRID_SOUTH_POLE].idx[i]] = vpolsud;
	       }
      }
   }
   if ((reffrom->GRTYP[0] == 'Z' || reffrom->GRTYP[0] == '#') && reffrom->RPNHeadExt.grref[0] == 'E' && refto->GRTYP[0] == 'B') {
      f77name(ez_corrbgd)(zout, &(refto->NX), &(refto->NY), &(refto->RPNHead.ig1));
   }

   return(ierc);
}

int32_t GeoRef_CalcPolarWindNorth(TGeoRef *Ref,float *polar_uu_in,float *polar_vv_in,float *uuin,float *vvin,int32_t ni,int32_t nj) {
  
   int32_t k1, k2;
   float *polar_wd, *polar_spd, *polar_uu, *polar_vv;
   double *polar_lat,*polar_lon,*polar_lat_gem, *polar_lon_gem, *polar_x, *polar_y;
   char grtypn[2],grtypa[2];
   float xlat1, xlat2, xlon1, xlon2;
   int32_t ig1n, ig2n, ig3n, ig4n;
   float pi, pj, d60, dgrw;
   int32_t i,j,ier;
   TGeoRef *gda, *gdps;
   float uupole, vvpole;
   double quatrevingtdix, zero;

   polar_uu  = (float *) malloc(4*ni*sizeof(float));
   polar_vv  = &polar_uu[ni];
   polar_wd  = &polar_uu[ni*2];
   polar_spd = &polar_uu[ni*3];
   polar_lat = (double*)malloc(4*ni*sizeof(double));
   polar_lon = &polar_lat[ni];
   polar_x   = &polar_lat[ni*2];
   polar_y   = &polar_lat[ni*3];

   for (i=0; i < ni; i++) {
      polar_x[i] = 1.0 * (i+1);
      polar_y[i] = 1.0 * nj;
   }
  
   GeoRef_XY2LL(Ref,polar_lat,polar_lon,polar_x,polar_y, ni,TRUE);

   if (Ref->GRTYP[0] == 'Z' && Ref->RPNHeadExt.grref[0] == 'E') {
      polar_lat_gem   = (double*) malloc(2*ni*sizeof(double));
      polar_lon_gem   = &polar_lat_gem[ni];
    
      for (i=0; i < ni; i++) {
         polar_lat_gem[i] = polar_lat[i];
         polar_lon_gem[i] = polar_lon[i];
      }
    
     f77name(cigaxg)(Ref->RPNHeadExt.grref, &xlat1, &xlon1, &xlat2, &xlon2, &Ref->RPNHeadExt.igref1, &Ref->RPNHeadExt.igref2, &Ref->RPNHeadExt.igref3, &Ref->RPNHeadExt.igref4,1);
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
   f77name(cxgaig)(grtypn, &ig1n, &ig2n, &ig3n, &ig4n, &pi, &pj, &d60, &dgrw,1);
   gdps = GeoRef_Create(ni, 1, grtypn, ig1n, ig2n, ig3n, ig4n, 0);
   GeoRef_WD2UV(gdps, polar_uu, polar_vv, polar_spd,  polar_wd, polar_lat, polar_lon, ni);

   f77name(ez_calcpoleval)(&uupole,polar_uu,&ni,Ref->AX,Ref->GRTYP,Ref->RPNHeadExt.grref,1,1);
   f77name(ez_calcpoleval)(&vvpole,polar_vv,&ni,Ref->AX,Ref->GRTYP,Ref->RPNHeadExt.grref,1,1);

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
   free(polar_uu);

   if (Ref->GRTYP[0] == 'Z' && Ref->RPNHeadExt.grref[0] == 'E' && polar_lat_gem != NULL)  {
     free(polar_lat_gem);
   }

   ier = GeoRef_Free(gdps);
   return(0);
}

int32_t GeoRef_CalcPolarWindSouth(TGeoRef *Ref,float *polar_uu_in,float *polar_vv_in,float *uuin,float *vvin,int32_t ni,int32_t nj) {

   int32_t k1, k2;
   float *polar_wd, *polar_spd,*polar_uu,*polar_vv;
   double *polar_lat,*polar_lon,*polar_lat_gem, *polar_lon_gem, *polar_x, *polar_y;
   char grtyps[2],grtypa[2];
   float xlat1, xlat2, xlon1, xlon2;
   int32_t ig1n, ig2n, ig3n, ig4n;
   float pi, pj, d60, dgrw;
   int32_t i,j,ier;
   TGeoRef *gda, *gdps;
   float uupole, vvpole;
   double quatrevingtdix, zero;

   polar_uu  = (float *) malloc(4*ni*sizeof(float));
   polar_vv  = &polar_uu[ni];
   polar_wd  = &polar_uu[ni*2];
   polar_spd = &polar_uu[ni*3];
   polar_lat = (double*)malloc(4*ni*sizeof(double));
   polar_lon = &polar_lat[ni];
   polar_x   = &polar_lat[ni*2];
   polar_y   = &polar_lat[ni*3];

   for (i=0; i < ni; i++) {
      polar_x[i] = 1.0 * (i+1);
      polar_y[i] = 1.0;
   }
  
   GeoRef_XY2LL(Ref,polar_lat,polar_lon,polar_x,polar_y,ni,TRUE);

   if (Ref->GRTYP[0] == 'Z' && Ref->RPNHeadExt.grref[0] == 'E') {
      polar_lat_gem   = (double*) malloc(2*ni*sizeof(double));
      polar_lon_gem   = &polar_lat_gem[ni];
    
     for (i=0; i < ni; i++) {
        polar_lat_gem[i] = polar_lat[i];
        polar_lon_gem[i] = polar_lon[i];
     }
    
     f77name(cigaxg)(Ref->RPNHeadExt.grref, &xlat1, &xlon1, &xlat2, &xlon2, &Ref->RPNHeadExt.igref1, &Ref->RPNHeadExt.igref2, &Ref->RPNHeadExt.igref3, &Ref->RPNHeadExt.igref4,1);
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
   f77name(cxgaig)(grtyps, &ig1n, &ig2n, &ig3n, &ig4n, &pi, &pj, &d60, &dgrw,1);
   gdps = GeoRef_Create(ni, 1, grtyps, ig1n, ig2n, ig3n, ig4n, 0);
   GeoRef_WD2UV(gdps, polar_uu, polar_vv, polar_spd,  polar_wd, polar_lat, polar_lon, ni);

   f77name(ez_calcpoleval)(&uupole,polar_uu,&ni,Ref->AX,Ref->GRTYP,Ref->RPNHeadExt.grref,1,1);
   f77name(ez_calcpoleval)(&vvpole,polar_vv,&ni,Ref->AX,Ref->GRTYP,Ref->RPNHeadExt.grref,1,1);

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

   free(polar_lat);
   free(polar_uu);

   if (Ref->GRTYP[0] == 'Z' && Ref->RPNHeadExt.grref[0] == 'E' && polar_lat_gem != NULL)  {
     free(polar_lat_gem);
   }

   ier = GeoRef_Free(gdps);
   return(0);
}

int32_t GeoRef_CorrVecNorth(TGeoSet *Set,float *uuout,float *vvout,float *uuin,float *vvin) {

   TGeoRef  *reffrom=NULL;
   float     uupole, vvpole;
   float    *polar_uu_in, *polar_vv_in, *polar_uu_out, *polar_vv_out, *corr_uus, *corr_vvs;
   double    ay[4];
   double   *temp_y;
   int32_t       ni, nj, i1, i2, j1, j2, degree,npts,i;
   int32_t       quatre = 4;
   int32_t       un = 1;

   if (!Set)
      return(0);
   
   reffrom=Set->RefFrom;

   ni = reffrom->NX;
   nj = reffrom->j2 - reffrom->j1 + 1;

   i1 = 1;
   i2 = ni;
   j1 = reffrom->j2-2;
   j2 = j1 + 3;
   degree = 3;

   npts = Set->zones[GRID_NORTH].npts;

   polar_uu_in = (float *) malloc(4 * ni * sizeof(float));
   polar_vv_in = (float *) malloc(4 * ni * sizeof(float));
   corr_uus = (float *) malloc(npts * sizeof(float));
   corr_vvs = (float *) malloc(npts * sizeof(float));

   GeoRef_CalcPolarWindNorth(reffrom,polar_uu_in,polar_vv_in,uuin,vvin,ni,nj);

   switch (Set->Opt.Interp) {
      case IR_CUBIC:
         switch (reffrom->GRTYP[0]) {
	          case 'Z':
	          case 'E':
	          case 'G':
               ay[0] = reffrom->AY[reffrom->j2-3];
               ay[1] = reffrom->AY[reffrom->j2-2];
               ay[2] = reffrom->AY[reffrom->j2-1];
               ay[3] = 90.0;
	             f77name(ez8_irgdint_3_wnnc)(corr_uus,Set->zones[GRID_NORTH].x,Set->zones[GRID_NORTH].y,&npts,polar_uu_in,&ni,&j1,&j2,reffrom->AX,ay,&reffrom->Extension);
	             f77name(ez8_irgdint_3_wnnc)(corr_vvs,Set->zones[GRID_NORTH].x,Set->zones[GRID_NORTH].y,&npts,polar_vv_in,&ni,&j1,&j2,reffrom->AX,ay,&reffrom->Extension);
	             break;

	          default:
	             f77name(ez8_rgdint_3)(corr_uus,Set->zones[GRID_NORTH].x,Set->zones[GRID_NORTH].y,&npts,polar_uu_in,&ni,&j1,&j2,&reffrom->Extension,&Set->Opt.NoData);
	             f77name(ez8_rgdint_3)(corr_vvs,Set->zones[GRID_NORTH].x,Set->zones[GRID_NORTH].y,&npts,polar_vv_in,&ni,&j1,&j2,&reffrom->Extension,&Set->Opt.NoData);
	             break;
	      }

      case IR_LINEAR:
         temp_y = (double*) malloc(npts*sizeof(double));
         for (i=0; i < npts; i++) {
            temp_y[i] = Set->zones[GRID_NORTH].y[i] - (1.0 * (reffrom->j2-3));
	       }
         f77name(ez8_rgdint_1)(corr_uus,Set->zones[GRID_NORTH].x,temp_y,&npts,polar_uu_in,&ni, &un, &quatre, &reffrom->Extension,&Set->Opt.NoData);
         f77name(ez8_rgdint_1)(corr_vvs,Set->zones[GRID_NORTH].x,temp_y,&npts,polar_vv_in,&ni, &un, &quatre, &reffrom->Extension,&Set->Opt.NoData);
         free(temp_y);
         
         break;

      case IR_NEAREST:
         temp_y = (double*) malloc(npts*sizeof(double));
         for (i=0; i < npts; i++) {
            temp_y[i] = Set->zones[GRID_NORTH].y[i] - (1.0 * (reffrom->j2-3));
	      }
         f77name(ez8_rgdint_0)(corr_uus,Set->zones[GRID_NORTH].x,temp_y,&npts,polar_uu_in,&ni, &un, &quatre,&Set->Opt.NoData);
         f77name(ez8_rgdint_0)(corr_vvs,Set->zones[GRID_NORTH].x,temp_y,&npts,polar_vv_in,&ni, &un, &quatre,&Set->Opt.NoData);
         free(temp_y);
         break;
   }


   for (i=0; i < Set->zones[GRID_NORTH].npts; i++) {
      uuout[Set->zones[GRID_NORTH].idx[i]] = corr_uus[i];
      vvout[Set->zones[GRID_NORTH].idx[i]] = corr_vvs[i];
   }

   free(polar_uu_in);
   free(polar_vv_in);
   free(corr_uus);
   free(corr_vvs);

   return(0);
}

int32_t GeoRef_CorrVecSouth(TGeoSet *Set,float *uuout,float *vvout,float *uuin,float *vvin) {
  
   TGeoRef  *reffrom=NULL;
   float *polar_uu_in, *polar_vv_in, *corr_uus, *corr_vvs, *temp_y;
   double ay[4];
   int32_t ni, nj, i1, i2, j1, j2, degree,npts,i;
   int32_t idx_gdin;

   if (!Set)
      return(0);
   
   reffrom=Set->RefFrom;

   npts = Set->zones[GRID_SOUTH].npts;
   ni = reffrom->NX;
   nj = reffrom->j2 - reffrom->j1 + 1;

   i1 = 1;
   i2 = ni;
   j1 = reffrom->j1 - 1;
   j2 = j1 + 3;
   degree = 3;

   polar_uu_in = (float *) malloc(4 * ni * sizeof(float));
   polar_vv_in = (float *) malloc(4 * ni * sizeof(float));
   corr_uus = (float *) malloc(npts * sizeof(float));
   corr_vvs = (float *) malloc(npts * sizeof(float));

   GeoRef_CalcPolarWindSouth(reffrom,polar_uu_in,polar_vv_in,uuin,vvin,ni,nj);

   switch (Set->Opt.Interp) {
      case IR_CUBIC:
         switch (reffrom->GRTYP[0]) {
	          case 'Z':
	          case 'E':
	          case 'G':
               ay[0] = -90.;
               ay[1] = reffrom->AY[0];
               ay[2] = reffrom->AY[1];
               ay[3] = reffrom->AY[2];
               f77name(ez8_irgdint_3_wnnc)(corr_uus,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,polar_uu_in,&ni,&j1,&j2,reffrom->AX,ay,&reffrom->Extension);
               f77name(ez8_irgdint_3_wnnc)(corr_vvs,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,polar_vv_in,&ni,&j1,&j2,reffrom->AX,ay,&reffrom->Extension);
               break;

            default:
               f77name(ez8_rgdint_3)(corr_uus,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,polar_uu_in,&ni,&j1,&j2,&reffrom->Extension,&Set->Opt.NoData);
               f77name(ez8_rgdint_3)(corr_vvs,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,polar_vv_in,&ni,&j1,&j2,&reffrom->Extension,&Set->Opt.NoData);
               break;
            }
      break;

      case IR_LINEAR:
         f77name(ez8_rgdint_1)(corr_uus,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,polar_uu_in,&ni,&j1,&j2,&reffrom->Extension,&Set->Opt.NoData);
         f77name(ez8_rgdint_1)(corr_vvs,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,polar_vv_in,&ni,&j1,&j2,&reffrom->Extension,&Set->Opt.NoData);
         break;

      case IR_NEAREST:
         f77name(ez8_rgdint_0)(corr_uus,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,polar_uu_in,&ni,&j1,&j2,&Set->Opt.NoData);
         f77name(ez8_rgdint_0)(corr_vvs,Set->zones[GRID_SOUTH].x,Set->zones[GRID_SOUTH].y,&npts,polar_vv_in,&ni,&j1,&j2,&Set->Opt.NoData);
         break;
   }


   for (i=0; i < Set->zones[GRID_SOUTH].npts; i++) {
      uuout[Set->zones[GRID_SOUTH].idx[i]] = corr_uus[i];
      vvout[Set->zones[GRID_SOUTH].idx[i]] = corr_vvs[i];
   }

   free(polar_uu_in);
   free(polar_vv_in);
   free(corr_uus);
   free(corr_vvs);

   return(0);
}

int32_t GeoRef_CorrectVector(TGeoSet *Set,float *uuout,float *vvout,float *uuin,float *vvin) {
 
   int32_t ier;
  
  
   if (Set->zones[GRID_NORTH].npts > 0) {
      ier = GeoRef_CorrVecNorth(Set,uuout,vvout,uuin,vvin);
   }
  
   if (Set->zones[GRID_SOUTH].npts > 0) {
     ier = GeoRef_CorrVecSouth(Set,uuout,vvout,uuin,vvin);
   }
  
   if (Set->zones[GRID_NORTH_POLE].npts > 0) {
     ier = GeoRef_CorrVecNorth(Set,uuout,vvout,uuin,vvin);
   }
  
   if (Set->zones[GRID_SOUTH_POLE].npts > 0) {
     ier = GeoRef_CorrVecSouth(Set,uuout,vvout,uuin,vvin);
   }
      
   return(0);
}
