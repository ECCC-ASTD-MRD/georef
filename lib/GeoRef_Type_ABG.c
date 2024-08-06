#include <App.h>
#include "GeoRef.h"

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for an A grid
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
int GeoRef_LL2XY_A(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   float  dellat,dellon,xlat0,xlon0;   

   dellon = 360.0/Ref->NX;
   xlon0  = 0.0;
   switch(Ref->RPNHead.ig1) {
      case GLOBAL: dellat = 180.0 / Ref->NY; xlat0 = -90.0 + dellat * 0.5; break;
      case NORTH:  dellat = 90.0 / Ref->NY;  xlat0 =  dellat * 0.5;        break;
      case SOUTH:  dellat = 90.0 / Ref->NY;  xlat0 = -90.0 + dellat * 0.5; break;
   }        
   
   GeoRef_LL2GD(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);

   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for a B grid
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
int GeoRef_LL2XY_B(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   float  dellat,dellon,xlat0,xlon0;   

   dellon = 360.0 / (Ref->NX-1);
   xlon0  = 0.0;
   switch(Ref->RPNHead.ig1) {
      case GLOBAL: dellat = 180.0 / (Ref->NY-1); xlat0 = -90.0; break;
      case NORTH:  dellat = 90.0 / (Ref->NY-1);  xlat0 = 0.0;   break;
      case SOUTH:  dellat = 90.0 / (Ref->NY-1);  xlat0 = -90.0; break;
   }        
   
   GeoRef_LL2GD(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);

   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for an G grid (Gaussian)
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
int GeoRef_LL2XY_G(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   float  dellat,dellon,xlat0,xlon0;   
   int    i,indy;

   dellon = 360.0 / Ref->NX;
   xlon0 = 0.0;

   switch(Ref->RPNHead.ig1) {
      case GLOBAL: 
         for(i=0;i<Nb;i++) {
            X[i] = (CLAMPLONREF(Lon[i],Ref->Options.LonRef) - xlon0)/dellon + 1.0;
            indy = GeoRef_XFind(Lat[i],Ref->AX,Ref->NX,1);
            if (indy>Ref->NY) indy = Ref->NY - 2;
         
            Y[i]= indy+(Lat[i]-Ref->AX[indy])/(Ref->AX[indy+1]-Ref->AX[indy]);
         }
         break;

      case NORTH:
         dellat = 90.0 / Ref->NY;
         xlat0 =  dellat * 0.5;
         GeoRef_LL2GD(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);
         break;

      case SOUTH:
         dellat = 90.0 / Ref->NY;
         xlat0 = -90.0 + dellat * 0.5;
         GeoRef_LL2GD(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);
         break;
   }

   return(0);
}
     