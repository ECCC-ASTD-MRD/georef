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

/**----------------------------------------------------------------------------
 * @brief  Get latitude and longitudes of a 2d point (Copy of fortran function, but using double precision)
 * @author J.D. HENDERSON
 * @date   February 1975
 *    @param[out]  Lat     Latitude in degrees 
 *    @param[out]  Lon     Longitude in degrees
 *    @param[in]   X       X coordinate
 *    @param[in]   Y       Y coordinate
 *    @param[in]   D60     Distance in meters between the gridpoints at 60 degree latitude
 *    @param[in]   DGRW    Angle between  the X axis and the Greenwich meridian
 *    @param[in]   HEM     Side of the hemisphere (1=North, 2=South)
 * 
 */
static inline void llfxy(double *Lat,double *Lon,double X,double Y,double D60,double DGRW,int HEM) {

   double re,re2,r2;

   re=1.866025*EARTHRADIUS/D60;
   re2=re*re;

   // If point is at pole set coords to (0.,90.)

   *Lat=90.0;
   *Lon=0.0;

   if (X!=0.0 || Y!=0.0) {

      // Calculate longitude in map coordinates

      if (X==0.0) *Lon=SIGN(90.0,Y);
      if (X!=0.0) *Lon=RAD2DEG(atan(Y/X));
      if (X<0.0)  *Lon=*Lon+SIGN(180.0,Y);

      // Adjust for grid orientation

      *Lon-=DGRW;
      if(*Lon > 180.0)  *Lon-=360.0;
      if(*Lon < -180.0) *Lon+=360.0;

      // Calculate latitude
      r2=X*X+Y*Y;
      *Lat=(re2-r2)/(re2+r2);
      *Lat=RAD2DEG(asin(*Lat));

   }
   // Change sign if in southern hemisphere
   if (HEM==2) {
      *Lat=-*Lat;
      *Lon=-*Lon;
   }
}

/**----------------------------------------------------------------------------
 * @brief  Get latitude and longitudes of each gridpoints of a LatLon grid (Copy of fortran function, but using double precision)
 * @author Michel Valin
 * @date   Janvier 1982
 *    @param[out]  Lat     Latitude stream
 *    @param[out]  Lon     Longitude stream
 *    @param[in]   NI      Size in x
 *    @param[in]   NJ      Size in y
 *    @param[in]   Lat0    Lower left corner latitude
 *    @param[in]   Lon0    Lower left corner longitude
 *    @param[in]   DLat    Latitude spacing in degrees
 *    @param[in]   Dlon    Longitudr spacing in degrees
 * 
 */
static inline void grll(double *Lat,double *Lon,int NI,int NJ,double Lat0,double Lon0,double DLat,double DLon) {

   int    i,j,idx=0;
   double lat;

   for(j=0;j<NJ;j++){
      lat=Lat0+j*DLat;
      for(i=0;i<NI;i++) {
         Lat[idx]=lat;
         Lon[idx]=fmod(Lon0+i*DLon,360.0);
         idx++;
      }
   }
}

/**----------------------------------------------------------------------------
 * @brief  Get latitude and longitudes of each gridpoints of a PS grid (Copy of fortran function, but using double precision)
 * @author Michel Valin
 * @date   Janvier 1982
 *    @param[out]  Lat     Latitude stream
 *    @param[out]  Lon     Longitude stream
 *    @param[in]   NI      Size in x
 *    @param[in]   NJ      Size in y
 *    @param[in]   PI      X Pole coordinate
 *    @param[in]   PJ      Y Pole coordinate
 *    @param[in]   D60     Distance in meters between the gridpoints at 60 degree latitude
 *    @param[in]   DGRW    Angle between  the X axis and the Greenwich meridian
 *    @param[in]   HEM     Side of the hemisphere (1=North, 2=South)
 * 
 */
static inline void grps(double *Lat,double *Lon,int NI,int NJ,double PI,double PJ,double D60,double DGRW,int HEM) {

   int    i,j,idx=0;
   double lat,lon,x,y;

   for(j=0;j<NJ;j++){
      y=j-PJ;
      for(i=0;i<NI;i++) {
         x=i-PI;
         llfxy(&lat,&lon,x,y,D60,DGRW,HEM);
         Lat[idx]=lat;
         Lon[idx]=(lon<0?lon+360:lon);
      }
   }
}

int GeoRef_CalcLL(TGeoRef* Ref) {

   float  xlat00, xlon00, dlat, dlon;
   int    i,j,k,ni, nj, npts, hemisphere;
   float *x,*y;
   double *lonp,*latp,*xp,*yp;

   if (!Ref->Lat) {
      ni = Ref->NX;
      nj = Ref->NY;
      npts = ni*nj;

      Ref->Lat = (double*)malloc(npts*sizeof(double));
      Ref->Lon = (double*)malloc(npts*sizeof(double));

      switch(Ref->GRTYP[0]) {
         case 'A':
         case 'B':
            grll(Ref->Lat,Ref->Lon,ni,nj,Ref->RPNHead.XG[X_SWLAT],Ref->RPNHead.XG[X_SWLON],Ref->RPNHead.XG[X_DLAT],Ref->RPNHead.XG[X_DLON]);
            break;

         case 'E':
            dlon = 360. /(ni-1);
            dlat = 180./(nj);
            xlon00 = 0.0;
            xlat00 = -90. + 0.5*dlat;

            grll(Ref->Lat,Ref->Lon,ni,nj,xlat00,xlon00,dlat,dlon);
            f77name(cigaxg)(Ref->GRTYP, &Ref->RPNHead.XG[X_LAT1], &Ref->RPNHead.XG[X_LON1],&Ref->RPNHead.XG[X_LAT2], &Ref->RPNHead.XG[X_LON2],&Ref->RPNHead.IG[X_IG1],  &Ref->RPNHead.IG[X_IG2], &Ref->RPNHead.IG[X_IG3], &Ref->RPNHead.IG[X_IG4]);
            c_ezgfllfxy(Ref->Lat,Ref->Lon,Ref->Lon.Ref->Lat,&npts,&Ref->RPNHead.XG[X_LAT1], &Ref->RPNHead.XG[X_LON1], &Ref->RPNHead.XG[X_LAT2],&Ref->RPNHead.XG[X_LON2]);
            break;

         case 'L':
            grll(Ref->Lat,Ref->Lon,ni,nj,Ref->RPNHead.XG[X_SWLAT],Ref->RPNHead.XG[X_SWLON],Ref->RPNHead.XG[X_DLAT],Ref->RPNHead.XG[X_DLON]);
            break;

         case 'N':
         case 'S':
            if (Ref->GRTYP[0] == 'N') {
	            hemisphere = 1;
	         } else {
               hemisphere = 2;
            }
            grps(Ref->Lat,Ref->Lon,ni,nj,Ref->RPNHead.XG[X_PI],Ref->RPNHead.XG[X_PJ],Ref->RPNHead.XG[X_D60],Ref->RPNHead.XG[X_DGRW],hemisphere);
            break;

         case 'T':
            xp = (double *) malloc(2*npts*sizeof(double));
            yp = &yp[npts];
            for (j=0; j < nj; j++) {
               for (i=0; i < ni; i++) {
                  k = j*ni + i;
                  yp[k] = 1.0 * (j+1);
                  xp[k] = 1.0 * (i+1);
               }
            }

            f77name(ez8_vtllfxy)(Ref->Lat,Ref->Lon,xp,yp,&Ref->RPNHead.XG[X_CLAT], &Ref->RPNHead.XG[X_CLON],&Ref->RPNHead.XG[X_TD60],&Ref->RPNHead.XG[X_TDGRW],&ni,&nj,&npts);
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
	               latp = (double *) malloc(2*npts*sizeof(double));
	               lonp = &latp[npts];
	               f77name(ez8_vllfxy)(latp,lonp,Ref->Lon,Ref->Lat,&ni,&nj,&Ref->RPNHead.XGREF[X_D60],&Ref->RPNHead.XGREF[X_DGRW],&Ref->RPNHead.XGREF[X_PI], &Ref->RPNHead.XGREF[X_PJ], &Ref->Hemi);
                  memcpy(Ref->Lon, lonp, npts*sizeof(double));
                  memcpy(Ref->Lat, latp, npts*sizeof(double));
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
                  c_ezgfllfxy(Ref->Lat,Ref->Lon,Ref->Lon,Ref->Lat,&npts,&Ref->RPNHead.XGREF[X_LAT1],&Ref->RPNHead.XGREF[X_LON1], &Ref->RPNHead.XGREF[X_LAT2], &Ref->RPNHead.XGREF[X_LON2]);
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
            //TODO : use double
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

int GeoRef_GetLL(TGeoRef *Ref,double *Lat,double *Lon) {

   int n;

   GeoRef_CalcLL(Ref);

   if (Ref->NbSub > 0) {
      n = Ref->Subs[0]->NX*Ref->Subs[1]->NY;

      GeoRef_GetLL(Ref->Subs[0],Lat,Lon);          // Yin
      GeoRef_GetLL(Ref->Subs[1],&Lat[n],&Lon[n]);  // Yang
   } else {
      if (Ref->Lat) {
         if (Lon) memcpy(Lon,Ref->Lon,Ref->NX*Ref->NY*sizeof(float));
         if (Lat) memcpy(Lat,Ref->Lat,Ref->NX*Ref->NY*sizeof(float));
      } else {
         App_Log(ERROR,"%s: Missing descriptors\n",__func__);
         return(-1);
      }
   }
   return(0);
}

int GeoRef_LLVal(TGeoRef *Ref,float *zout,float *zin,double *Lat,double *Lon,int Nb) { 

   double *x,*y;
   int ier;

   x = (double*)malloc(2*Nb*sizeof(double));
   y = &x[Nb];
   
   ier = GeoRef_LL2XY(Ref,x,y,Lat,Lon,Nb);
   ier = GeoRef_XYVal(Ref,zout,zin,x,y,Nb);
   
   free(x);

   return(0);
}

int GeoRef_LLUVVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int Nb) {

   double *x, *y;
   int ier;
   
   if (Ref->NbSub > 0) {
      App_Log(ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   } else {

      x = (double*)malloc(2*Nb*sizeof(double));
      y = &x[Nb];
   
      ier = GeoRef_LL2XY(Ref,x,y,Lat,Lon,Nb);
      ier = GeoRef_XYUVVal(Ref,uuout,vvout,uuin,vvin,x,y,Nb);
   
      free(x);

      return(0);
   }
}

int GeoRef_LLWDVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int Nb) {

   TGeoRef *yin_gd, *yan_gd;
   int ier,j;
   double *x,*y;
   float  *uuyin, *vvyin, *uuyan, *vvyan;

   if (Ref->NbSub > 0) {
      x = (double*) malloc(2*Nb*sizeof(double));
      y = &x[Nb];
      uuyin = (float*) malloc(6*Nb*sizeof(float));
      vvyin = &uuyin[Nb];
      uuyan = &uuyin[Nb*2];
      vvyan = &uuyin[Nb*3];
      ier = GeoRef_LL2XY(Ref,x,y,Lat,Lon,Nb);
      ier = GeoRef_XYUVVal(Ref,uuout,vvout,uuin,vvin,x,y,Nb);
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      ier = GeoRef_UV2WD(yin_gd,uuyin,vvyin,uuout,vvout,Lat,Lon,Nb);
      ier = GeoRef_UV2WD(yan_gd,uuyan,vvyan,uuout,vvout,Lat,Lon,Nb);
      for (j=0; j< Nb; j++) {
         if (y[j] > yin_gd->NY) {
            uuout[j]=uuyan[j];
            vvout[j]=vvyan[j];
         } else {
            uuout[j]=uuyin[j];
            vvout[j]=vvyin[j];
         }
      }
      free(x);
      free(uuyin);
   } else {
      ier = GeoRef_LLUVVal(Ref,uuout,vvout,uuin,vvin,Lat,Lon,Nb);
      ier = GeoRef_UV2WD(Ref,uuout,vvout,uuout,vvout,Lat,Lon,Nb);
   }

   return(0);
}

