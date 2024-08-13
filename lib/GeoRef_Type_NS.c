#include <App.h>
#include "GeoRef.h"

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for an N or S grid (Polar stereographic)
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 *    @param[in]  X       X array
 *    @param[in]  Y       Y array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int GeoRef_XY2LL_NS(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {
   
   int un=1;
   f77name(ez8_vllfxy)(Lat,Lon,X,Y,&Nb,&un,&Ref->RPNHeadExt.xg3,&Ref->RPNHeadExt.xg4,&Ref->RPNHeadExt.xg1,&Ref->RPNHeadExt.xg2,&Ref->Hemi);

   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms LatLon coordinates to XY for an N or S grid (Polar stereographic)
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] X       X array
 *    @param[out] Y       Y array
 *    @param[in]  Lat     Latitude array
 *    @param[in]  Lon     Longitude array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int GeoRef_LL2XY_NS(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   float  pi,pj,dgrw,d60;   

   f77name(cigaxg)(Ref->GRTYP,&pi,&pj,&d60,&dgrw,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4,1);
   f77name(ez8_vxyfll)(X,Y,Lat,Lon,&Nb,&d60,&dgrw,&pi,&pj,&Ref->Hemi);
 
   return(0);
}

