/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Fonctions et definitions relatives aux fichiers standards et rmnlib
 * Fichier      : GeoRef_Type_E.h
 * Creation     : October 2020 - J.P. Gauthier
 *
 * Description:
 *
 * License:
 *    This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation,
 *    version 2.1 of the License.
 *
 *    This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the
 *    Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 *    Boston, MA 02111-1307, USA.
 *
 *==============================================================================
 */

#include "App.h"
#include "GeoRef.h"

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for an E grid
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
int GeoRef_LL2XY_E(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {
 
   float  dellat,dellon,xlat0,xlon0,xlat1,xlon1,xlat2,xlon2;   

   f77name(cigaxg)(Ref->GRTYP,&xlat1,&xlon1,&xlat2,&xlon2,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);         
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
 * @author Jean-Philippe Gauthier
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

   GeoRef_RotateInvertXY(Lat,Lon,Lon,Lat,Nb,Ref->RPNHead.XG[X_LAT1],Ref->RPNHead.XG[X_LON1],Ref->RPNHead.XG[X_LAT2],Ref->RPNHead.XG[X_LON2]);

   return(0);
}
