#include <App.h>
#include "GeoRef.h"

/*----------------------------------------------------------------------------
 * @brief  Calculates MAGL height of a grid coordinate
 * @date   June 2014
 *    @param[in]  Ref     Georeference pointer
 *    @param[in]  ZRef    Vertical reference pointer
 *    @param[in]  Azimuth Azimuth (x) coordinate in grid unit
 *    @param[in]  Bin     Bin (y) coordinate in grid unit
 *    @param[in]  Sweep   Sweep (z) coordinate in grid unit
 * 
 *    @return     Height (m)
*/
double GeoRef_RDRHeight(TGeoRef *Ref,TZRef *ZRef,double Azimuth,double Bin,double Sweep) {

   if (Bin>=0 && Bin<Ref->R && Sweep>=0 && Sweep<ZRef->LevelNb) {
      return(Ref->Loc.Elev+sin(DEG2RAD(ZRef->Levels[(int)Sweep]))*((int)Bin*Ref->ResR));
   } else {
      return(0);
   }
}

/*----------------------------------------------------------------------------
 * @brief  Create an R type  georeference (Radial)
 * @date   June 2014
 *    @param[in]     Lat        Center latitude
 *    @param[in]     Lon        Center longitude
 *    @param[in]     Height     Center altitude (m)
 *    @param[in]     R          Radius (m)
 *    @param[in]     ResR       Resolution along the radius
 *    @param[in]     ResA       Resolution in azimuth
 * 
 *    @return        GeoRef object (NULL=Error)
*/
TGeoRef* GeoRef_CreateR(double Lat,double Lon,double Height,int32_t R,double ResR,double ResA) {

   TGeoRef *ref;

   ref=GeoRef_New();

   ref->GRTYP[0]='R';
   ref->Loc.Lat=Lat;
   ref->Loc.Lon=Lon;
   ref->Loc.Elev=Height;
   ref->R=R;
   ref->ResR=ResR;
   ref->ResA=ResA;

   GeoRef_Size(ref,0,0,360/ResA,R-1,0);
   GeoRef_Qualify(ref);

   ref->Height=GeoRef_RDRHeight;

   return(ref);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for a R grid (Radial)
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 *    @param[in]  X       X array
 *    @param[in]  Y       Y array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int32_t GeoRef_XY2LL_R(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb) {

   TCoord loc0;
   double x,y,d;
   int32_t    n;

   loc0.Lat=DEG2RAD(Ref->Loc.Lat);
   loc0.Lon=DEG2RAD(Ref->Loc.Lon);

   #pragma omp parallel for default(none) private(n,x,y,d) shared(Nb,Ref,X,Y,Lat,Lon,loc0)
   for(n=0;n<Nb;n++) {


      x=X[n]*Ref->ResA;
      y=Y[n]*Ref->ResR;

      x=DEG2RAD(x);
      d=M2RAD(y*Ref->CTH);

      if (Ref->Options.Transform) {
         Lat[n]=asin(sin(loc0.Lat)*cos(d)+cos(loc0.Lat)*sin(d)*cos(x));
         Lon[n]=fmod(loc0.Lon+(atan2(sin(x)*sin(d)*cos(loc0.Lat),cos(d)-sin(loc0.Lat)*sin(*Lat)))+M_PI,M_2PI)-M_PI;
         Lat[n]=RAD2DEG(Lat[n]);
         Lon[n]=RAD2DEG(Lon[n]);
      } else {
         Lat[n]=d;
         Lon[n]=x;
      }
   }

   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms LatLon coordinates to XY for a R grid (Radial)
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] X       X array
 *    @param[out] Y       Y array
 *    @param[in]  Lat     Latitude array
 *    @param[in]  Lon     Longitude array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int32_t GeoRef_LL2XY_R(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb) {

   TCoord  loc0;
   double  x,d,lat,lon;
   int32_t     n;

   loc0.Lat=DEG2RAD(Ref->Loc.Lat);
   loc0.Lon=DEG2RAD(Ref->Loc.Lon);

   #pragma omp parallel for default(none) private(n,x,d,lat,lon) shared(Nb,Ref,X,Y,Lat,Lon,loc0)
   for(n=0;n<Nb;n++) {
      lat=DEG2RAD(Lat[n]);
      lon=DEG2RAD(Lon[n]);

      d=fabs(DIST(0.0,loc0.Lat,loc0.Lon,lat,lon));
      x=-RAD2DEG(COURSE(loc0.Lat,loc0.Lon,lat,lon));
      X[n]=x<0.0?x+360.0:x;
      Y[n]=d/Ref->CTH;

      if (Ref->Options.Transform) {
         X[n]/=Ref->ResA;
         Y[n]/=Ref->ResR;
      }
   }

   return(0);
}
