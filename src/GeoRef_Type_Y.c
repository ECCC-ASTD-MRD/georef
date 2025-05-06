#include <App.h>
#include "GeoRef.h"

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for a Y grid
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 *    @param[in]  X       X array
 *    @param[in]  Y       Y array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int32_t GeoRef_XY2LL_Y(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb) {

   int32_t     i,s,sx,sy,indx;
   double *tmpx,*tmpy,dx,dy;

   if (!Ref->AX || !Ref->AY) {
      return(FALSE);
   }

   switch (Ref->RPNHeadExt.grref[0]) {
      case 'L':
         for (i=0; i < Nb; i++) {
            indx=ROUND(Y[i])*(Ref->X1-Ref->X0)+ROUND(X[i]);
            Lat[i]=Ref->AY[indx];
            Lon[i]=Ref->AX[indx];
         }
         break;

      case 'W':
         tmpx = (double*)malloc(2*Nb*sizeof(double));
         tmpy = &tmpx[Nb];

         for (i=0; i < Nb; i++) {
            sx=floor(X[i]);sx=CLAMP(sx,Ref->X0,Ref->X1);
            sy=floor(Y[i]);sy=CLAMP(sy,Ref->Y0,Ref->Y1);
            dx=X[i]-sx;;
            dy=Y[i]-sy;

            s=sy*Ref->NX+sx;
            tmpx[i]=Ref->AX[s];
            tmpy[i]=Ref->AY[s];

            if (++sx<=Ref->X1) {
               s=sy*Ref->NX+sx;
               tmpx[i]+=(Ref->AX[s]-tmpx[i])*dx;
            }
            tmpx[i];

            if (++sy<=Ref->Y1) {
               s=sy*Ref->NX+(sx-1);
               tmpy[i]+=(Ref->AY[s]-tmpy[i])*dy;
            }
            tmpy[i];
         }
         GeoRef_XY2LL_W(Ref,Lat,Lon,tmpx,tmpy,Nb);
         free(tmpx);
         break;

      default:
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Undefined reference grid type: %s\n",__func__,Ref->RPNHeadExt.grref[0]);
         return(-1);
         break;
   }
   return(Nb);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms LatLon coordinates to XY for a Y grid
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] X       X array
 *    @param[out] Y       Y array
 *    @param[in]  Lat     Latitude array
 *    @param[in]  Lon     Longitude array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int32_t GeoRef_LL2XY_Y(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb) {

   int32_t idx,n,nb=0;
   double  dists[8];

   for (n=0;n<Nb;n++) {
      if (GeoRef_Nearest(Ref,Lon[n],Lat[n],&idx,dists,1,Ref->Options.DistTreshold)) {
         if (dists[0]<1.0) {
            Y[n]=(int)(idx/Ref->NX);
            X[n]=idx-Y[n]*Ref->NX;
            nb++;
         }
      }
   }
   return(nb);
}
