
/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Fonctions et definitions relatives aux fichiers standards et rmnlib
 * Fichier      : GeoRef_Type_ABG.h
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

 int GeoRef_LL2XY_A(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   float  dellat,dellon,xlat0,xlon0;   

   dellon = 360.0/Ref->NX;
   xlon0  = 0.0;
   switch(Ref->RPNHead.IG[X_IG1]) {
      case GLOBAL: dellat = 180.0 / Ref->NY; xlat0 = -90.0 + dellat * 0.5; break;
      case NORTH:  dellat = 90.0 / Ref->NY;  xlat0 =  dellat * 0.5;        break;
      case SOUTH:  dellat = 90.0 / Ref->NY;  xlat0 = -90.0 + dellat * 0.5; break;
   }        
   
   GeoRef_LL2GD(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);

   return(0);
}

int GeoRef_LL2XY_B(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   float  dellat,dellon,xlat0,xlon0;   

   dellon = 360.0 / (Ref->NX-1);
   xlon0  = 0.0;
   switch(Ref->RPNHead.IG[X_IG1]) {
      case GLOBAL: dellat = 180.0 / (Ref->NY-1); xlat0 = -90.0; break;
      case NORTH:  dellat = 90.0 / (Ref->NY-1);  xlat0 = 0.0;   break;
      case SOUTH:  dellat = 90.0 / (Ref->NY-1);  xlat0 = -90.0; break;
   }        
   
   GeoRef_LL2GD(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);

   return(0);
}

int GeoRef_LL2XY_G(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   float  dellat,dellon,xlat0,xlon0;   
   int    i,indy;

   dellon = 360.0 / Ref->NX;
   xlon0 = 0.0;

   switch(Ref->RPNHead.IG[X_IG1]) {
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
     