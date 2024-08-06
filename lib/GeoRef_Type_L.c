#include <App.h>
#include "GeoRef.h"

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for an L grid (Cylindrical)
 * @author Jean-Philippe Gauthier
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 *    @param[in]  X       X array
 *    @param[in]  Y       Y array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int GeoRef_XY2LL_L(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {

   int i;

   #pragma omp parallel for default(none) private(i) shared(Nb,Ref,X,Y,Lat,Lon)
   for (i=0; i < Nb; i++) {
      Lat[i] = (Y[i]-1.0)*Ref->RPNHeadExt.xg3+Ref->RPNHeadExt.xg1;
      Lon[i] = (X[i]-1.0)*Ref->RPNHeadExt.xg4+Ref->RPNHeadExt.xg2;
   }

   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms LatLon coordinates to XY for a L grid (Cylindrical)
 * @author Jean-Philippe Gauthier
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] X       X array
 *    @param[out] Y       Y array
 *    @param[in]  Lat     Latitude array
 *    @param[in]  Lon     Longitude array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int GeoRef_LL2XY_L(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   float  dellat,dellon,xlat0,xlon0;   

   f77name(cigaxg)(Ref->GRTYP,&xlat0,&xlon0,&dellat,&dellon,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4,1);         
   GeoRef_LL2GD(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);

   return(0);
}
