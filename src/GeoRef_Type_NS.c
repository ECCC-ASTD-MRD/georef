/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Fonctions et definitions relatives aux fichiers standards et rmnlib
 * Fichier      : GeoRef_Type_NS.h
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

int GeoRef_XY2LL_NS(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {
   
   int un=1;
   f77name(ez8_vllfxy)(Lat,Lon,X,Y,&Nb,&un,&Ref->RPNHead.XG[X_D60],&Ref->RPNHead.XG[X_DGRW],&Ref->RPNHead.XG[X_PI],&Ref->RPNHead.XG[X_PJ],&Ref->Hemi);

   return(0);
}

int GeoRef_LL2XY_NS(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   float  pi,pj,dgrw,d60;   

   f77name(cigaxg)(Ref->GRTYP,&pi,&pj,&d60,&dgrw,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);
   f77name(ez8_vxyfll)(X,Y,Lat,Lon,&Nb,&d60,&dgrw,&pi,&pj,&Ref->Hemi);
 
   return(0);
}

