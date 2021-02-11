/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Fonctions et definitions relatives aux fichiers standards et rmnlib
 * Fichier      : GeoRef_Type_M.h
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

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for an M grid (Mesh or TIN)
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
int GeoRef_LL2XY_M(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {
 
   TQTree *node;
   Vect3d  b;
   int     n,d,idx;

   for(d=0;d<Nb;d++) {
      X[d]=-1.0;
      Y[d]=-1.0;

      if (Ref->QTree) {
         // If there's an index use it
         if ((node=QTree_Find(Ref->QTree,Lon[d],Lat[d])) && node->NbData) {
            
            // Loop on this nodes data payload
            for(n=0;n<node->NbData;n++) {
               idx=(intptr_t)node->Data[n].Ptr-1; // Remove false pointer increment

               if (Bary_Get(b,Ref->Wght?Ref->Wght[idx/3]:0.0,Lon[d],Lat[d],Ref->AX[Ref->Idx[idx]],Ref->AY[Ref->Idx[idx]],
                  Ref->AX[Ref->Idx[idx+1]],Ref->AY[Ref->Idx[idx+1]],Ref->AX[Ref->Idx[idx+2]],Ref->AY[Ref->Idx[idx+2]])) {
                  
                  // Return coordinate as triangle index + barycentric coefficient
                  X[d]=idx+b[0]+1;
                  Y[d]=idx+b[1]+1;
                  break;
               }
            }
         }
      } else {
         // Otherwise loop on all
         for(idx=0;idx<Ref->NIdx-3;idx+=3) {
            if (Bary_Get(b,Ref->Wght?Ref->Wght[idx/3]:0.0,Lon[d],Lat[d],Ref->AX[Ref->Idx[idx]],Ref->AY[Ref->Idx[idx]],
               Ref->AX[Ref->Idx[idx+1]],Ref->AY[Ref->Idx[idx+1]],Ref->AX[Ref->Idx[idx+2]],Ref->AY[Ref->Idx[idx+2]])) {

               // Return coordinate as triangle index + barycentric coefficient
               X[d]=idx+b[0]+1;
               Y[d]=idx+b[1]+1;
               break;
            }
         }
      }            
   } 
   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms LatLon coordinates to XY for an M grid (Mesh or TIN)
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
int GeoRef_XY2LL_M(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {

   int i;

   for (i=0; i < Nb; i++) {
      Lon[i]=VertexValS(Ref->AX,NULL,Ref->NX,Ref->NY,X[i]-1.0,Y[i]-1.0,TRUE);
      Lat[i]=VertexValS(Ref->AY,NULL,Ref->NX,Ref->NY,X[i]-1.0,Y[i]-1.0,FALSE);
   }
   return(0);
}
