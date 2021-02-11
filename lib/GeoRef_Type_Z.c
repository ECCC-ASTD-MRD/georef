
/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Fonctions et definitions relatives aux fichiers standards et rmnlib
 * Fichier      : GeoRef_Type_Z.h
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
 * @brief  Transforms XY grid coordinates to LatLon for a Z grid
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
int GeoRef_XY2LL_Z(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {

   int     i,indx,indy,un=1;
   double  delxx,delyy,*tmpx,*tmpy,*ytmp;
   float   xlat1, xlon1, xlat2, xlon2;

   tmpx = (double*)malloc(3*Nb*sizeof(double));
   tmpy = &tmpx[Nb];
   ytmp = &tmpx[Nb*2];
   for (i=0; i < Nb; i++) {
      indx = (int)X[i]-1;
      ytmp[i] = Y[i];
      if (Ref->RPNHead.IG[X_IG2] == 1) {
         ytmp[i] = Ref->NY +1.0 - Y[i];
      }
      indy = (int)ytmp[i]-1;
      indx = indx < 0 ? 0 : indx;
      indy = indy < 0 ? 0 : indy;
      indx = indx > Ref->NX-2 ? Ref->NX-2 : indx;
      indy = indy > Ref->j2-2 ? Ref->j2-2 : indy;
      delxx = Ref->AX[indx+1]-Ref->AX[indx];
      tmpx[i] = Ref->AX[indx] + ((X[i]-1.0-indx)*delxx);

      delyy = Ref->AY[indy+1]-Ref->AY[indy];
      tmpy[i] = Ref->AY[indy] + ((ytmp[i]-1.0-indy)*delyy);
   }

   switch (Ref->RPNHead.GRREF[0]) {
      case 'E':
         f77name(cigaxg)(Ref->RPNHead.GRREF,&xlat1,&xlon1,&xlat2,&xlon2,&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4]);
         GeoRef_RotateInvertXY(Lat,Lon,tmpx,tmpy,Nb,Ref->RPNHead.XGREF[X_LAT1],Ref->RPNHead.XGREF[X_LON1],Ref->RPNHead.XGREF[X_LAT2],Ref->RPNHead.XGREF[X_LON2]);
         break;

      case 'S':
      case 'N':
         f77name(ez8_vllfxy)(Lat,Lon,tmpx,tmpy,&Nb,&un,&Ref->RPNHead.XGREF[X_D60],&Ref->RPNHead.XGREF[X_DGRW],&Ref->RPNHead.XGREF[X_PI],&Ref->RPNHead.XGREF[X_PJ],&Ref->Hemi);
         break;

      case 'L':
         for (i=0; i < Nb; i++) {
            Lat[i] = (tmpy[i])*Ref->RPNHead.XGREF[X_DLAT]+Ref->RPNHead.XGREF[X_SWLAT];
            Lon[i] = (tmpx[i])*Ref->RPNHead.XGREF[X_DLON]+Ref->RPNHead.XGREF[X_SWLON];
         }
         break;

      case 'W':
         GeoRef_XY2LL_W(Ref,Lat,Lon,tmpx,tmpy,Nb);
         break;

      default:
         App_Log(ERROR,"%s: Undefined reference grid type: %s\n",__func__,Ref->RPNHead.GRREF[0]);
         free(tmpx);
         return(1);
         break;
   }
   free(tmpx);

   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms LatLon coordinates to XY for a Z grid
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
int GeoRef_LL2XY_Z(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   int   i,j,indx,indy,d;
    
   GeoRef_LL2GREF(Ref,X,Y,Lat,Lon,Nb);

   // Look into expansion descriptor
   for(i=0;i<Nb;i++) {
      d=(Ref->Type&GRID_AXY2D)?Ref->NX:1;

      //TODO: clarify NX and j2 index          
      indx = GeoRef_XFind(X[i],Ref->AX,Ref->NX,1);
      indy = GeoRef_XFind(Y[i],Ref->AY,Ref->j2,d);
     
      if (indx >= Ref->NX-1) indx = Ref->NX - 2;
      if (indy >= Ref->j2-1) indy = Ref->j2 - 2;

      X[i] = indx+(X[i]-Ref->AX[indx])/(Ref->AX[indx+1]-Ref->AX[indx])+1;
      Y[i] = indy+(Y[i]-Ref->AY[indy*d])/(Ref->AY[(indy+1)*d]-Ref->AY[indy*d])+1;
   } 

   if (Ref->GRTYP[0] == 'G') {
      if (Ref->RPNHead.IG[X_IG1] == 1) {
         for (j=0; j < Nb; j++) Y[j] = Y[j] - Ref->j2;
      }
      if (Ref->RPNHead.IG[X_IG2] == 1)  {
         for (j=0; j < Nb; j++) Y[j] = Ref->j2 +1.0 - Y[j];
      }
      // TODO: From GeoRef_RPN Fix for G grid 0-360 1/5 gridpoint problem
      for (j=0; j < Nb; j++) if (X[j]>Ref->X1+0.5) X[j]-=(Ref->X1+1);
   }
      
   return(0);
}
