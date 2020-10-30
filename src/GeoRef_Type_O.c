/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Fonctions et definitions relatives aux fichiers standards et rmnlib
 * Fichier      : GeoRef_Type_O.h
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
#include "Vertex.h"

int GeoRef_XY2LL_O(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {

   int i;

   for (i=0; i < Nb; i++) {
      Lon[i]=Vertex_ValS(Ref->AX,NULL,Ref->NX,Ref->NY,X[i],Y[i],TRUE);
      Lat[i]=Vertex_ValS(Ref->AY,NULL,Ref->NX,Ref->NY,X[i],Y[i],FALSE);
   }
   return(0);
}

int GeoRef_LL2XY_O(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   int     out=0;
   int     x,y,n,nd,d,dx,dy,idx,idxs[8];
   double  dists[8],xx,yy;
   Vect2d  pts[4],pt;

   GeoRef_LL2GREF(Ref,X,Y,Lat,Lon,Nb);

   for(d=0;d<Nb;d++) {

      X[d]=-1.0;
      Y[d]=-1.0;

      if ((nd=GeoRef_Nearest(Ref,Lon[d],Lat[d],idxs,dists,8,0.0))) {
     
         pt[0]=Lon[d];
         pt[1]=Lat[d];

         // Find which cell includes coordinates
         for(n=0;n<nd;n++) {
            idx=idxs[n];

            // Find within which quad
            dx=-1;dy=-1;
            if (!GeoRef_WithinCell(Ref,pt,pts,idx-Ref->NX-1,idx-1,idx,idx-Ref->NX)) {
            
               dx=0;dy=-1;
               if (!GeoRef_WithinCell(Ref,pt,pts,idx-Ref->NX,idx,idx+1,idx-Ref->NX+1)) {
                  
                  dx=-1;dy=0;
                  if (!GeoRef_WithinCell(Ref,pt,pts,idx-1,idx+Ref->NX-1,idx+Ref->NX,idx)) {
               
                     dx=0;dy=0;
                     if (!GeoRef_WithinCell(Ref,pt,pts,idx,idx+Ref->NX,idx+Ref->NX+1,idx+1)) {
                        idx=-1;
                     }
                  }
               }
            }
            
            // If found, exit loop
            if (idx!=-1) {
               break;
            }
         }

         if (idx!=-1) {
            // Map coordinates to grid
            Vertex_Map(pts,&xx,&yy,Lon[d],Lat[d]);
            X[d]=xx;Y[d]=yy;

            if (!ISNAN(X[d]) && !ISNAN(Y[d])) {
               y=idx/Ref->NX;
               x=idx-y*Ref->NX;
               Y[d]+=y+dy+1;
               X[d]+=x+dx+1; 
            } else {
               App_Log(ERROR,"%s: Invalid coordinate (NAN): ll(%f,%f) xy(%f,%f) %i\n",__func__,Lat[d],Lon[d],X[d],Y[d],idx);
               X[d]=-1,0;
               Y[d]=-1.0;
               out++;
            }
         } else {
 //           App_Log(ERROR,"%s: Point not found: %f %f %i\n",__func__,Lat[d],Lon[d],idx);
            out++;
         }
      } 
   }
   App_Log(DEBUG,"%s: Points out: %i\n",__func__,out);

   return(0);
}
