#include <App.h>
#include "GeoRef.h"
#include "Vertex.h"
#include "Triangle.h"

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for an M grid (Mesh or TIN)
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 *    @param[in]  X       X array
 *    @param[in]  Y       Y array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int32_t GeoRef_LL2XY_M(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb) {
 
   TQTree *node;
   Vect3d  b;
   int32_t n,d,idx;

   #pragma omp parallel for default(none) private(d,node,b,n,idx) shared(Nb,Ref,X,Y,Lat,Lon)
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
         for(idx=0;idx<Ref->NIdx;idx+=3) {
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
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] X       X array
 *    @param[out] Y       Y array
 *    @param[in]  Lat     Latitude array
 *    @param[in]  Lon     Longitude array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int32_t GeoRef_XY2LL_M(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb) {

   int32_t i;

   #pragma omp parallel for default(none) private(i) shared(Nb,Ref,X,Y,Lat,Lon)
   for (i=0; i < Nb; i++) {
      Lon[i]=VertexValS(Ref->AX,NULL,Ref->NX,Ref->NY,X[i]-1.0,Y[i]-1.0,TRUE);
      Lat[i]=VertexValS(Ref->AY,NULL,Ref->NX,Ref->NY,X[i]-1.0,Y[i]-1.0,TRUE);
   }
   return(0);
}
