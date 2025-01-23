#include <App.h>
#include "GeoRef.h"
#include "Vertex.h"

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for an O grid
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 *    @param[in]  X       X array
 *    @param[in]  Y       Y array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int32_t GeoRef_XY2LL_O(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb) {

   int32_t i;

   if (!Ref->AX || !Ref->AY) {
      return(FALSE);
   }
   
   #pragma omp parallel for default(none) private(i) shared(Nb,Ref,X,Y,Lat,Lon)
   for (i=0; i < Nb; i++) {
      Lon[i]=VertexValS(Ref->AX,NULL,Ref->NX,Ref->NY,X[i]-1.0,Y[i]-1.0,TRUE);
      Lat[i]=VertexValS(Ref->AY,NULL,Ref->NX,Ref->NY,X[i]-1.0,Y[i]-1.0,TRUE);
   }
   return(Nb);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms LatLon coordinates to XY for an O grid
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] X       X array
 *    @param[out] Y       Y array
 *    @param[in]  Lat     Latitude array
 *    @param[in]  Lon     Longitude array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int32_t GeoRef_LL2XY_O(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb) {

   int32_t     out=0;
   int32_t     x,y,n,nd,d,dx,dy,idx,idxs[8];
   double  dists[8],xx,yy;
   Vect2d  pts[4],pt;

   GeoRef_LL2GREF(Ref,X,Y,Lat,Lon,Nb);
   #pragma omp parallel for default(none) private(n,nd,d,dx,dy,x,y,xx,yy,idx,idxs,pt,pts,dists) shared(Nb,Ref,X,Y,Lat,Lon) reduction(+:out)
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
               Y[d]+=y+dy+1.0;
               X[d]+=x+dx+1.0; 
            } else {
               Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid coordinate (NAN): ll(%f,%f) xy(%f,%f) %i\n",__func__,Lat[d],Lon[d],X[d],Y[d],idx);
               X[d]=-1.0;
               Y[d]=-1.0;
               out++;
            }
         } else {
            Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Point not found: %f %f\n",__func__,Lat[d],Lon[d]);
            out++;
         }
      } 
   }
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Points out: %i\n",__func__,out);

   return(out);
}
