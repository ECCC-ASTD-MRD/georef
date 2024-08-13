#include <App.h>
#include "GeoRef.h"

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for a T grid
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 *    @param[in]  X       X array
 *    @param[in]  Y       Y array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int GeoRef_XY2LL_T(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {

   f77name(ez8_vtllfxy)(Lat,Lon,X,Y,&Ref->RPNHeadExt.xg3,&Ref->RPNHeadExt.xg4,&Ref->RPNHeadExt.xg1,&Ref->RPNHeadExt.xg2,&Ref->NX,&Ref->NY,&Nb);
   
   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms LatLon coordinates to XY for a T grid
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] X       X array
 *    @param[out] Y       Y array
 *    @param[in]  Lat     Latitude array
 *    @param[in]  Lon     Longitude array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int GeoRef_LL2XY_T(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   float  clat,clon,dgrw,d60;   
   int    hemi;

   f77name(cigaxg)(Ref->GRTYP,&d60,&dgrw,&clat,&clon,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4,1);
   f77name(ez8_vtxyfll)(X,Y,Lat,Lon,&clat,&clon,&d60,&dgrw,&Ref->NX,&Ref->NY,&Nb);

   return(0);
}
