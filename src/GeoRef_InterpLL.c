#include <App.h>
#include "GeoRef.h"

/**----------------------------------------------------------------------------
 * @brief  Get latitude and longitudes of a 2d point32_t (Copy of fortran function, but using double precision)
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
static inline void llfxy(double *Lat,double *Lon,double X,double Y,double D60,double DGRW,int32_t HEM) {

   double re,re2,r2;

   re=1.866025*EARTHRADIUS/D60;
   re2=re*re;

   // Initialize point32_t at pole (0.0,90.0)
   *Lat=90.0;
   *Lon=0.0;

   if (X!=0.0 || Y!=0.0) {

      // Calculate longitude in map coordinates
      if (X==0.0) {
         *Lon=SIGN(90.0,Y);
      } else {
         *Lon=RAD2DEG(atan(Y/X));
      }
      if (X<0.0) *Lon=*Lon+SIGN(180.0,Y);

      // Adjust for grid orientation
      *Lon-=DGRW;
      *Lon=CLAMPLON(*Lon);

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
static inline void grll(double *Lat,double *Lon,int32_t NI,int32_t NJ,double Lat0,double Lon0,double DLat,double DLon) {

   int32_t    i,j,idx=0;
   double lat;

   for(j=0;j<NJ;j++){
      lat=Lat0+j*DLat;
      for(i=0;i<NI;i++) {
         Lat[idx]=lat;
         Lon[idx]=Lon0+i*DLon;
         idx++;
      }
   }
}

/**----------------------------------------------------------------------------
 * @brief  Get latitude and longitudes of each gridpoints of a PS grid (Copy of fortran function, but using double precision)
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
static inline void grps(double *Lat,double *Lon,int32_t NI,int32_t NJ,double PI,double PJ,double D60,double DGRW,int32_t HEM) {

   int32_t idx = 0;
   double x,y;

   for(int32_t j = 0; j < NJ; j++) {
      y=j+1-PJ;
      for(int32_t i = 0; i < NI; i++) {
         x=i+1-PI;
         llfxy(&Lat[idx],&Lon[idx],x,y,D60,DGRW,HEM);
         idx++;
      }
   }
}

static inline void Permut(double *Z,int32_t NI,int32_t NJ) {

   int32_t    i,j,ncc,idx,idxn;
   double t;

   ncc = NJ>>1;

   for(j=0;j<ncc;j++) {
      for(i=0;i<NI;i++) {
          idx=(NJ+1-j)+NI+i;
          idxn=j*NI+i;

          t = Z[idx];
          Z[idx] = Z[idxn];
          Z[idxn] = t;
      }
   }
}

/**----------------------------------------------------------------------------
 * @brief  Calculer la position latlon de tous les points de grille.
 * @date   June 2015
 *    @param[in]  Ref     Pointeur sur la reference geographique
 *
 *    @return             Number of coordinates
*/
int32_t GeoRef_CalcLL(TGeoRef* Ref) {

   float   xlat00, xlon00, dlat, dlon;
   int32_t i,j,k,ni,nj,npts, hemisphere;
   double *xp,*yp;

   ni = Ref->NX;
   nj = Ref->NY;
   npts = ni*nj;

   if (!Ref->Lat) {

      Ref->Lat = (double*)malloc(npts*sizeof(double));
      Ref->Lon = (double*)malloc(npts*sizeof(double));

      if (!Ref->Lat || !Ref->Lon) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable to allocate coordinate buffers\n",__func__);
         return(0);
      }

      switch(Ref->GRTYP[0]) {
         case 'A':
         case 'B':
            grll(Ref->Lat,Ref->Lon,ni,nj,Ref->RPNHeadExt.xg1,Ref->RPNHeadExt.xg2,Ref->RPNHeadExt.xg3,Ref->RPNHeadExt.xg4);
            break;

         case 'E':
            dlon = 360.0/(ni-1);
            dlat = 180.0/(nj);
            xlon00 = 0.0;
            xlat00 = -90. + 0.5*dlat;

            grll(Ref->Lat,Ref->Lon,ni,nj,xlat00,xlon00,dlat,dlon);
            f77name(cigaxg)(Ref->GRTYP,&Ref->RPNHeadExt.xg1,&Ref->RPNHeadExt.xg2,&Ref->RPNHeadExt.xg3,&Ref->RPNHeadExt.xg4,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4,1);
            GeoRef_RotateInvertXY(Ref->Lat,Ref->Lon,Ref->Lon,Ref->Lat,npts,Ref->RPNHeadExt.xg1,Ref->RPNHeadExt.xg2,Ref->RPNHeadExt.xg3,Ref->RPNHeadExt.xg4);
            break;

         case 'L':
            grll(Ref->Lat,Ref->Lon,ni,nj,Ref->RPNHeadExt.xg1,Ref->RPNHeadExt.xg2,Ref->RPNHeadExt.xg3,Ref->RPNHeadExt.xg4);
            break;

         case 'N':
         case 'S':
            if (Ref->GRTYP[0] == 'N') {
	            hemisphere = 1;
	         } else {
               hemisphere = 2;
            }
            grps(Ref->Lat,Ref->Lon,ni,nj,Ref->RPNHeadExt.xg1,Ref->RPNHeadExt.xg2,Ref->RPNHeadExt.xg3,Ref->RPNHeadExt.xg4,hemisphere);
            break;

         case 'T':
            xp = (double *) malloc(2*npts*sizeof(double));
            yp = &yp[npts];
            k=0;
            for (j=0; j<nj; j++) {
               for (i=0; i<ni; i++) {
                  yp[k] = j+1.0;
                  xp[k] = i+1.0;
                  k++;
               }
            }

            f77name(ez8_vtllfxy)(Ref->Lat,Ref->Lon,xp,yp,&Ref->RPNHeadExt.xg3, &Ref->RPNHeadExt.xg4,&Ref->RPNHeadExt.xg1,&Ref->RPNHeadExt.xg2,&ni,&nj,&npts);
            free(xp);
            break;

         case 'M':
            for (i=0;i<npts;i++) {
               Ref->Lat[i]=Ref->AY[i];
               Ref->Lon[i]=Ref->AX[i];
            }
            break;

         case 'W':
            npts=0;
            for (j=0;j<nj;j++) {
               for (i=0;i<ni;i++) {
                  Ref->Lat[npts]=j;
                  Ref->Lon[npts]=i;
                  npts++;
               }
            }
            // We can reuse the same lat lon array because the transform function makes internal copies
            GeoRef_XY2LL_W(Ref,Ref->Lat,Ref->Lon,Ref->Lon,Ref->Lat,npts);
            break;

         case 'R':
            npts=0;
            for (j=0;j<nj;j++) {
               for (i=0;i<ni;i++) {
                  Ref->Lat[npts]=j;
                  Ref->Lon[npts]=i;
                  npts++;
               }
            }
            // We can reuse the same lat lon array because the transform function makes internal copies
            GeoRef_XY2LL_R(Ref,Ref->Lat,Ref->Lon,Ref->Lon,Ref->Lat,npts);
            break;

         case 'O':
         case 'Y':
            switch (Ref->RPNHeadExt.grref[0]) {
               case 'N':
               case 'S':
                  Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Operation not supported - Y/O grid on PS grid\n",__func__);
                  return(-1);
                  break;

               case 'X':
               case 'L':
                  for (i=0; i < npts; i++) {
                     Ref->Lat[i]=Ref->AY[i];
                     Ref->Lon[i]=Ref->AX[i];
                  }
                  break;

               case 'E':
                  Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Operation not supported - Y/O grid on E grid\n",__func__);
                  return(-1);
                  break;

            }
            break;

         case 'G':
            k=0;
            for (j=0; j<nj; j++) {
               for (i=0; i<ni; i++) {
                  Ref->Lat[k] = (Ref->RPNHead.ig1 == GRID_NORTH)?Ref->AY[j+nj]:Ref->AY[j];
                  Ref->Lat[k] = Ref->RPNHead.ig2==1?-Ref->Lat[k]:Ref->Lat[k];
                  Ref->Lon[k] = Ref->AX[i];
                  k++;
               }
            }

         case '#':
         case 'Z':
            k=0;
            for (j=0; j<nj; j++) {
               for (i=0; i<ni; i++) {
                  Ref->Lat[k] = Ref->AY[j];
                  Ref->Lon[k] = Ref->AX[i];
                  k++;
               }
            }

            switch (Ref->RPNHeadExt.grref[0]) {
	            case 'N':
	            case 'S':
	               f77name(ez8_vllfxy)(Ref->Lat,Ref->Lon,Ref->Lon,Ref->Lat,&ni,&nj,&Ref->RPNHeadExt.xgref3,&Ref->RPNHeadExt.xgref4,&Ref->RPNHeadExt.xgref1, &Ref->RPNHeadExt.xgref2, &Ref->Hemi);
                  break;

               case 'L':
                  for (i=0; i < npts; i++) {
                     Ref->Lon[i] = Ref->RPNHeadExt.xgref2 + Ref->RPNHeadExt.xgref4 * Ref->Lon[i];
                     Ref->Lat[i] = Ref->RPNHeadExt.xgref1 + Ref->RPNHeadExt.xgref3 * Ref->Lat[i];
                  }
                  break;

               case 'W':
                  GeoRef_XY2LL_W(Ref,Ref->Lat,Ref->Lon,Ref->AX,Ref->AY,npts);
                  break;

               case 'E':
                  GeoRef_RotateInvertXY(Ref->Lat,Ref->Lon,Ref->Lon,Ref->Lat,npts,Ref->RPNHeadExt.xgref1,Ref->RPNHeadExt.xgref2,Ref->RPNHeadExt.xgref3,Ref->RPNHeadExt.xgref4);
                  break;
            }
            break;

         case 'U':
            npts=0;
            for (j=0;j<Ref->NbSub;j++) {
               k=GeoRef_CalcLL(Ref->Subs[j]);
               memcpy(&Ref->Lat[npts],Ref->Subs[j]->Lat,k*sizeof(double));
               memcpy(&Ref->Lon[npts],Ref->Subs[j]->Lon,k*sizeof(double));
               npts+=k;
            }
            break;

         case '!':
            xp = (double*) malloc(2*npts*sizeof(double));
            yp = &xp[npts];
            k=0;
            for (j=0; j < nj; j++) {
               for (i=0; i < ni; i++) {
                  xp[k] = i+1.0;
                  yp[k] = j+1.0;
                  k++;
               }
            }
            f77name(ez8_llflamb)(Ref->Lat,Ref->Lon,xp,yp,&npts,Ref->GRTYP,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4);
            break;
      }

      switch(Ref->GRTYP[0]) {
         case 'G':
         case 'B':
         case 'A':
            if (Ref->RPNHead.ig2 == 1) {
               Permut(Ref->Lat,Ref->NX,Ref->NY);
            }
            break;

         default:
            break;
      }

      for (i=0; i < npts; i++) {
         CLAMPLON(Ref->Lon[i]);
      }
   }

   return(npts);
}

/**----------------------------------------------------------------------------
 * @brief  Récupérer la position latlon de tous les points de grille.
 * @date   June 2015
 *    @param[in]  Ref     Pointeur sur la reference geographique
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array

 *    @return             Number of coordinates
*/
int32_t GeoRef_GetLL(TGeoRef *Ref,double *Lat,double *Lon) {

   TGeoRef *ref;
   int32_t  n=0,i;

   ref=GeoRef_SubGet(Ref);

   if (ref->NbSub > 0) {
      i = ref->Subs[0]->NX*ref->Subs[0]->NY;

      n=GeoRef_GetLL(ref->Subs[0],Lat,Lon);           // Yin
      n+=GeoRef_GetLL(ref->Subs[1],&Lat[i],&Lon[i]);  // Yang
   } else {
      n=GeoRef_CalcLL(ref);

      if (ref->Lat) {
         if (Lon) memcpy(Lon,ref->Lon,n*sizeof(double));
         if (Lat) memcpy(Lat,ref->Lat,n*sizeof(double));
      } else {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Missing descriptors\n",__func__);
         return(-1);
      }
   }
   return(n);
}

/**----------------------------------------------------------------------------
 * @brief  Interpolate values at latlon positions
 * @date   June 2015
 *    @param[in]  Ref     Geographic reference pointer
 *    @param[out] zout    Array of interpolated values
 *    @param[in]  zin     Array of input values on the source grid
 *    @param[in]  Lat     Latitudes of interpolation
 *    @param[in]  Lon     Longitudes of interpolation
 *    @param[in]  Nb      Number on position tu interpolate to
 *
 *    @return             Number of interpolated values
*/
int32_t GeoRef_LLVal(TGeoRef *Ref,TGeoOptions *Opt,float *zout,float *zin,double *Lat,double *Lon,int32_t Nb) {

   double *x,*y;

   if (!Opt) Opt=&Ref->Options;
   if (!Opt) Opt=&GeoRef_Options;

   x = (double*)malloc(2*Nb*sizeof(double));
   y = &x[Nb];

   GeoRef_LL2XY(Ref,x,y,Lat,Lon,Nb,TRUE);
   GeoRef_XYVal(Ref,Opt,zout,zin,x,y,Nb);

   free(x);

   return(0);
}

/**----------------------------------------------------------------------------
 * @brief  Interpolate vector values at latlon positions
 * @date   June 2015
 *    @param[in]  Ref     Geographic reference pointer
 *    @param[out] uuout   Array of interpolated U values
 *    @param[out] vvout   Array of interpolated V values
 *    @param[in]  uuin    Array of input U values on the source grid
 *    @param[in]  vvin    Array of input V values on the source grid
 *    @param[in]  Lat     Latitudes of interpolation
 *    @param[in]  Lon     Longitudes of interpolation
 *    @param[in]  Nb      Number on position tu interpolate to
 *
 *    @return             Number of interpolated values
*/
int32_t GeoRef_LLUVVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int32_t Nb) {

   double *x,*y;

   if (Ref->NbSub > 0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   } else {

      if (!Opt) Opt=&Ref->Options;
      if (!Opt) Opt=&GeoRef_Options;

      if (!(x=(double*)malloc(2*Nb*sizeof(double)))) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable to allocate temporary coordinate buffer\n",__func__);
         return(-1);
      }
      y = &x[Nb];

      GeoRef_LL2XY(Ref,x,y,Lat,Lon,Nb,TRUE);
      GeoRef_XYUVVal(Ref,Opt,uuout,vvout,uuin,vvin,x,y,Nb);

      free(x);

      return(0);
   }
}

/**----------------------------------------------------------------------------
 * @brief  Interpolate vector values as speed and directions at latlon positions
 * @date   June 2015
 *    @param[in]  Ref     Geographic reference pointer
 *    @param[out] uuout   Array of interpolated speed values
 *    @param[out] vvout   Array of interpolated direction values
 *    @param[in]  uuin    Array of input U values on the source grid
 *    @param[in]  vvin    Array of input V values on the source grid
 *    @param[in]  Lat     Latitudes of interpolation
 *    @param[in]  Lon     Longitudes of interpolation
 *    @param[in]  Nb      Number on position tu interpolate to
 *
 *    @return             Number of interpolated values
*/
int32_t GeoRef_LLWDVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int32_t Nb) {

   TGeoRef *yin_gd, *yan_gd,*ref;
   int32_t j;
   double *x,*y;
   float  *uuyin, *vvyin, *uuyan, *vvyan;

   if (!Opt) Opt=&Ref->Options;
   if (!Opt) Opt=&GeoRef_Options;
   ref=GeoRef_SubGet(Ref);

   if (ref->NbSub > 0) {
      x = (double*) malloc(2*Nb*sizeof(double));
      y = &x[Nb];
      uuyin = (float*) malloc(6*Nb*sizeof(float));
      vvyin = &uuyin[Nb];
      uuyan = &uuyin[Nb*2];
      vvyan = &uuyin[Nb*3];
      GeoRef_LL2XY(ref,x,y,Lat,Lon,Nb,TRUE);
      GeoRef_XYUVVal(ref,Opt,uuout,vvout,uuin,vvin,x,y,Nb);
      yin_gd=ref->Subs[0];
      yan_gd=ref->Subs[1];

      GeoRef_UV2WD(yin_gd,uuyin,vvyin,uuout,vvout,Lat,Lon,Nb);
      GeoRef_UV2WD(yan_gd,uuyan,vvyan,uuout,vvout,Lat,Lon,Nb);
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
      GeoRef_LLUVVal(ref,Opt,uuout,vvout,uuin,vvin,Lat,Lon,Nb);
      GeoRef_UV2WD(ref,uuout,vvout,uuout,vvout,Lat,Lon,Nb);
   }

   return(0);
}

