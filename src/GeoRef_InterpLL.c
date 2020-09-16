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
#include "RPN.h"
#include "GeoRef.h"


int GeoRef_CalcLL(TGeoRef* Ref) {

   float  xlat00, xlon00, dlat, dlon;
   int    i,j,k,ni, nj, npts, hemisphere;
   float *lonp,*latp,*xp,*yp,*x,*y;

   if (!Ref->Lat) {
      ni = Ref->NX;
      nj = Ref->NY;
      npts = ni*nj;

      Ref->Lat = (float *)malloc(npts*sizeof(float));
      Ref->Lon = (float *)malloc(npts*sizeof(float));

      switch(Ref->GRTYP[0]) {
         case 'A':
         case 'B':
            f77name(grll)(Ref->Lat,Ref->Lon,&ni,&nj,&Ref->RPNHead.XG[X_SWLAT],&Ref->RPNHead.XG[X_SWLON], &Ref->RPNHead.XG[X_DLAT], &Ref->RPNHead.XG[X_DLON]);
            break;

         case 'E':
            dlon = 360. /(ni-1);
            dlat = 180./(nj);
            xlon00 = 0.0;
            xlat00 = -90. + 0.5*dlat;

            f77name(grll)(Ref->Lat,Ref->Lon,&ni,&nj,&xlat00,&xlon00,&dlat,&dlon);
            f77name(cigaxg)(Ref->GRTYP, &Ref->RPNHead.XG[X_LAT1], &Ref->RPNHead.XG[X_LON1],&Ref->RPNHead.XG[X_LAT2], &Ref->RPNHead.XG[X_LON2],&Ref->RPNHead.IG[X_IG1],  &Ref->RPNHead.IG[X_IG2], &Ref->RPNHead.IG[X_IG3], &Ref->RPNHead.IG[X_IG4]);
            latp = (float *)malloc(2*npts*sizeof(float));
            lonp = &latp[npts];
            f77name(ez_gfllfxy)(lonp,latp,Ref->Lon,Ref->Lat,&npts,&Ref->RPNHead.XG[X_LAT1], &Ref->RPNHead.XG[X_LON1], &Ref->RPNHead.XG[X_LAT2],&Ref->RPNHead.XG[X_LON2]);
            memcpy(Ref->Lat,latp,npts*sizeof(float));
            memcpy(Ref->Lon,lonp,npts*sizeof(float));
            free(latp);
            break;

         case 'L':
            f77name(grll)(Ref->Lat,Ref->Lon,&ni,&nj,&Ref->RPNHead.XG[X_SWLAT],&Ref->RPNHead.XG[X_SWLON],&Ref->RPNHead.XG[X_DLAT], &Ref->RPNHead.XG[X_DLON]);
            break;

         case 'N':
         case 'S':
            if (Ref->GRTYP[0] == 'N') {
	            hemisphere = 1;
	         } else {
               hemisphere = 2;
            }
            f77name(grps)(Ref->Lat,Ref->Lon,&ni,&nj,&Ref->RPNHead.XG[X_PI],&Ref->RPNHead.XG[X_PJ],&Ref->RPNHead.XG[X_D60], &Ref->RPNHead.XG[X_DGRW], &hemisphere);
            break;

         case 'T':
            xp = (float *) malloc(2*npts*sizeof(float));
            yp = &yp[npts];
            for (j=0; j < nj; j++) {
               for (i=0; i < ni; i++) {
                  k = j*ni + i;
                  yp[k] = 1.0 * (j+1);
                  xp[k] = 1.0 * (i+1);
               }
            }

            f77name(ez_vtllfxy)(Ref->Lat,Ref->Lon,xp,yp,&Ref->RPNHead.XG[X_CLAT], &Ref->RPNHead.XG[X_CLON],&Ref->RPNHead.XG[X_TD60],&Ref->RPNHead.XG[X_TDGRW],&ni,&nj,&npts);
            free(xp);
            break;

         case 'O':
         case 'Y':
            switch (Ref->RPNHead.GRREF[0]) {
               case 'N':
               case 'S':
                  App_Log(ERROR,"%s: Operation not supported - Y/O grid on PS grid\n",__func__);
                  return(-1);
                  break;

               case 'L':
                  memcpy(Ref->Lon, Ref->AX, npts*sizeof(float));
                  memcpy(Ref->Lat, Ref->AY, npts*sizeof(float));
                  for (i=0; i < npts; i++) {
                     if (Ref->Lon[i] < 0.0) {
                        Ref->Lon[i] = Ref->AX[i] + 360.0;
                     }
                  }
                  break;

               case 'E':
                  App_Log(ERROR,"%s: Operation not supported - Y/O grid on E grid\n",__func__);
                  return(-1);
                  break;

            }
            break;

         case '#':
         case 'Z':
         case 'G':
            for (j=0; j < nj; j++) {
               for (i=0; i < ni; i++) {
                    Ref->Lat[C_TO_FTN(i,j,ni)] = Ref->AY[j];
                    Ref->Lon[C_TO_FTN(i,j,ni)] = Ref->AX[i];
               }
            }

	         if (Ref->GRTYP[0] == 'G' && Ref->RPNHead.IG[X_IG1] == NORTH) {
	            for (j=0; j < nj; j++) {
	               for (i=0; i < ni; i++) {
		               Ref->Lat[C_TO_FTN(i,j,ni)] = Ref->AY[j+nj];
		            }
	            }
	         }

            switch (Ref->RPNHead.GRREF[0]) {
	            case 'N':
	            case 'S':
	               latp = (float *) malloc(2*npts*sizeof(float));
	               lonp = &latp[npts];
	               f77name(ez_vllfxy)(latp,lonp,Ref->Lon,Ref->Lat,&ni,&nj,&Ref->RPNHead.XGREF[X_D60],&Ref->RPNHead.XGREF[X_DGRW],&Ref->RPNHead.XGREF[X_PI], &Ref->RPNHead.XGREF[X_PJ], &Ref->Hemi);

	               for (i=0; i < npts; i++) {
		               if (lonp[i] < 0.0) lonp[i] += 360.0;
		            }

                  memcpy(Ref->Lon, lonp, npts*sizeof(float));
                  memcpy(Ref->Lat, latp, npts*sizeof(float));
                  free(latp);
                  break;

               case 'L':
                  for (j=0; j < nj; j++) {
                     for (i=0; i < ni; i++) {
                        Ref->Lat[C_TO_FTN(i,j,ni)] += 1.0;
                        Ref->Lon[C_TO_FTN(i,j,ni)] += 1.0;
                     }
                  }
  
                  for (i=0; i < npts; i++) {
                     Ref->Lon[i] = Ref->RPNHead.XGREF[X_SWLON] + Ref->RPNHead.XGREF[X_DLON] * (Ref->Lon[i]-1.0);
                     Ref->Lon[i]=  fmod(fmod(Ref->Lon[i],360.0)+360.0,360.0);
                     Ref->Lat[i] = Ref->RPNHead.XGREF[X_SWLAT] + Ref->RPNHead.XGREF[X_DLAT] * (Ref->Lat[i]-1.0);
                  }   
                  break;

               case 'E':
                  latp = (float *) malloc(2*npts*sizeof(float));
                  lonp = &latp[npts];
                  f77name(ez_gfllfxy)(lonp,latp,Ref->Lon,Ref->Lat,&npts,&Ref->RPNHead.XGREF[X_LAT1],&Ref->RPNHead.XGREF[X_LON1], &Ref->RPNHead.XGREF[X_LAT2], &Ref->RPNHead.XGREF[X_LON2]);
                  memcpy(Ref->Lon,lonp,npts*sizeof(float));
                  memcpy(Ref->Lat,latp,npts*sizeof(float));
                  free(latp);
                  break;
            }
            break;

         case '!':
            x = (float *) malloc(2*npts*sizeof(float));
            y = &x[npts];
            for (j=0; j < nj; j++) {
               for (i=0; i < ni; i++) {
                  x[C_TO_FTN(i,j,ni)] = (float) (i+1.0);
                  y[C_TO_FTN(i,j,ni)] = (float) (j+1.0);
               }
            }
            f77name(ez_llflamb)(Ref->Lat,Ref->Lon,x,y,&npts,&Ref->GRTYP,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4],1);
            for (i=0; i < npts; i++) {
               if (Ref->Lon[i] < 0.0) {
                  Ref->Lon[i] += 360.0;
               }
            }
            break;
      }

      switch(Ref->GRTYP[0]) {
         case 'G':
         case 'B':
         case 'A':
         if (Ref->RPNHead.IG[X_IG2] == 1) {
            f77name(permut)(Ref->Lat, &Ref->NX, &Ref->NY);
         }
         break;

         default:
         break;
      }
   }

   if (App_LogLevel(NULL)==EXTRA) {
      App_Log(EXTRA,"Grid Lat Lon coordinates %p\n",Ref);

      for (i=0; i < Ref->NX*Ref->NY; i++) { 
         App_Log(MUST,"%f %f\n",Ref->Lat[i],Ref->Lon[i]);
      }
   }
   return(0);
}

int GeoRef_GetLL(TGeoRef *Ref,float *lat,float *lon) {

   int icode;
   int n;
      
   if (Ref->NbSub > 0) {
      n = Ref->Subs[0]->NX*Ref->Subs[1]->NY;

      icode=GeoRef_GetLLN(Ref->Subs[0],lat,lon);          // Yin
      icode=GeoRef_GetLLN(Ref->Subs[1],&lat[n],&lon[n]);  // Yang
   } else {
      icode=GeoRef_GetLLN(Ref,lat,lon);
   } 
   return(icode);
}

int GeoRef_GetLLN(TGeoRef *Ref,float *Lat,float *Lon) {

   GeoRef_CalcLL(Ref);

   if (Ref->Lat) {
      if (Lon) memcpy(Lon,Ref->Lon,Ref->NX*Ref->NY*sizeof(float));
      if (Lat) memcpy(Lat,Ref->Lat,Ref->NX*Ref->NY*sizeof(float));
   } else {
      App_Log(ERROR,"%s: Missing descriptors\n",__func__);
      return(-1);
   }
   return(0);
}

int GeoRef_LLVal(TGeoRef *Ref,float *zout,float *zin,float *lat,float *lon,int n) { 

   float *x,*y;
   int ier;

   x = (float *)malloc(2*n*sizeof(float));
   y = &x[n];
   
   if (Ref->NbSub > 0 ) {
      ier = GeoRef_LL2XY(Ref,x,y,lat,lon,n);
   } else {
      ier = GeoRef_LL2XYN(Ref,x,y,lat,lon,n,FALSE);
   }
   ier = GeoRef_XYVal(Ref,zout,zin,x,y,n);
   
   free(x);

   return(0);
}

int GeoRef_LLUVVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,float *lat,float *lon,int n) {

   float *x, *y;
   int ier;
   
   if (Ref->NbSub > 0) {
      App_Log(ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   } else {

      x = (float *) malloc(n * sizeof(float));
      y = (float *) malloc(n * sizeof(float));
   
      ier = GeoRef_LL2XYN(Ref, x, y, lat, lon, n,FALSE);
      ier = GeoRef_XYUVVal(Ref, uuout, vvout, uuin, vvin, x, y, n);
   
      free(x);
      free(y);

      return(0);
   }
}

int GeoRef_LLWDVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,float *lat,float *lon,int n) {

   TGeoRef *yin_gd, *yan_gd;
   int ier,j;
   float *x, *y;
   float *uuyin, *vvyin, *uuyan, *vvyan;

   if (Ref->NbSub > 0) {
      x = (float *) malloc(n * sizeof(float));
      y = (float *) malloc(n * sizeof(float));
      uuyin = (float *) malloc(n*sizeof(float));
      vvyin = (float *) malloc(n*sizeof(float));
      uuyan = (float *) malloc(n*sizeof(float));
      vvyan = (float *) malloc(n*sizeof(float));
      ier = GeoRef_LL2XY(Ref, x, y, lat, lon, n);
      ier = GeoRef_XYUVVal(Ref, uuout, vvout, uuin, vvin, x, y, n);
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      ier = GeoRef_UV2WD(yin_gd,uuyin,vvyin,uuout,vvout,lat,lon,n);
      ier = GeoRef_UV2WD(yan_gd,uuyan,vvyan,uuout,vvout,lat,lon,n);
      for (j=0; j< n; j++) {
         if (y[j] > yin_gd->NY) {
            uuout[j]=uuyan[j];
            vvout[j]=vvyan[j];
         } else {
            uuout[j]=uuyin[j];
            vvout[j]=vvyin[j];
         }
      }
      free(uuyin); free(vvyin);
      free(uuyan); free(vvyan);
   } else {
      ier = GeoRef_LLUVVal(Ref, uuout, vvout, uuin, vvin, lat, lon, n);
      ier = GeoRef_UV2WD(Ref, uuout, vvout, uuout, vvout, lat, lon, n);
   }

   return(0);
}

