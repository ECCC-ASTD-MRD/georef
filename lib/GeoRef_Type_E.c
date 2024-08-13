#include <App.h>
#include "GeoRef.h"

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for an E grid
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 *    @param[in]  X       X array
 *    @param[in]  Y       Y array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int GeoRef_LL2XY_E(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {
 
   float  dellat,dellon,xlat0,xlon0,xlat1,xlon1,xlat2,xlon2;   

   f77name(cigaxg)(Ref->GRTYP,&xlat1,&xlon1,&xlat2,&xlon2,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4,1);         
   GeoRef_RotateXY(Lat,Lon,X,Y,Nb,xlat1,xlon1,xlat2,xlon2);

   dellon = 360.0 / (Ref->NX-1);
   xlon0 = 0.0;

   dellat = 180.0 / (Ref->NY);
   xlat0 = -90. + 0.5*dellat;

   GeoRef_LL2GD(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);
   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms LatLon coordinates to XY for an E grid
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] X       X array
 *    @param[out] Y       Y array
 *    @param[in]  Lat     Latitude array
 *    @param[in]  Lon     Longitude array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/int GeoRef_XY2LL_E(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {

   int    i;
   double dlat,dlon,swlat,swlon;

   for (i=0; i < Nb; i++) {
      dlat  = 180.0 / Ref->NY;
      dlon  = 360.0 / (Ref->NX - 1);
      swlat = -90.0 + 0.5 * dlat;
      swlon = 0.0;
      Lon[i] = (X[i]-1.0)*dlon+swlon;
      Lat[i] = (Y[i]-1.0)*dlat+swlat;
   }

   GeoRef_RotateInvertXY(Lat,Lon,Lon,Lat,Nb,Ref->RPNHeadExt.xg1,Ref->RPNHeadExt.xg2,Ref->RPNHeadExt.xg3,Ref->RPNHeadExt.xg4);

   return(0);
}
