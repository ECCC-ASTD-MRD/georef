#include <App.h>
#include "GeoRef.h"

void  GeoRef_RotateXY(double *Lat,double *Lon,double *X,double *Y,int npts,float xlat1,float xlon1,float xlat2,float xlon2) {

   double cart[3],carot[3],latr,lonr,cosdar;
   int    n;
   float  r[3][3],ri[3][3];

   f77name(ez_crot)(r,ri,&xlon1,&xlat1,&xlon2,&xlat2);

   #pragma omp parallel for default(none) private(n,latr,lonr,cosdar,cart,carot) shared(r,ri,npts,Lat,Lon,X,Y)
   for(n=0;n<npts;n++) {
      latr=DEG2RAD(Lat[n]);
      lonr=DEG2RAD(Lon[n]);
      cosdar = cos(latr);

      cart[0] = cosdar*cos(lonr);
      cart[1] = cosdar*sin(lonr);
      cart[2] = sin(latr);

      carot[0] = r[0][0]*cart[0]+r[1][0]*cart[1]+r[2][0]*cart[2];
      carot[1] = r[0][1]*cart[0]+r[1][1]*cart[1]+r[2][1]*cart[2];
      carot[2] = r[0][2]*cart[0]+r[1][2]*cart[1]+r[2][2]*cart[2];
      
      Y[n]=RAD2DEG(asin(fmax(-1.0,fmin(1.0,carot[2]))));
      X[n]=RAD2DEG(atan2(carot[1],carot[0]));
      X[n]=fmod(X[n],360.0);
      if (X[n]<0.0) X[n]+=360.0;
   }
}
 
void  GeoRef_RotateInvertXY(double *Lat,double *Lon,double *X,double *Y,int npts,float xlat1,float xlon1,float xlat2,float xlon2) {

   double cart[3],carot[3],latr,lonr,cosdar;
   int i,k,n,c;
   float r[3][3],ri[3][3];

   f77name(ez_crot)(r,ri,&xlon1,&xlat1,&xlon2,&xlat2);

   #pragma omp parallel for default(none) private(n,latr,lonr,cosdar,cart,carot) shared(ri,npts,Lat,Lon,X,Y)
   for(n=0;n<npts;n++) {
      latr=DEG2RAD(Y[n]);
      lonr=DEG2RAD(X[n]);
      cosdar = cos(latr);

      cart[0] = cosdar*cos(lonr);
      cart[1] = cosdar*sin(lonr);
      cart[2] = sin(latr);

      carot[0] = ri[0][0]*cart[0]+ri[1][0]*cart[1]+ri[2][0]*cart[2];
      carot[1] = ri[0][1]*cart[0]+ri[1][1]*cart[1]+ri[2][1]*cart[2];
      carot[2] = ri[0][2]*cart[0]+ri[1][2]*cart[1]+ri[2][2]*cart[2];
     
      Lat[n]=RAD2DEG(asin(fmax(-1.0,fmin(1.0,carot[2]))));
      Lon[n]=RAD2DEG(atan2(carot[1],carot[0]));
   }
}

int GeoRef_LL2GREF(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   int   i,hemi;
   float pi,pj,dgrw,d60,dlat,dlon,xlat0,xlon0,xlat1,xlon1,xlat2,xlon2;
    
   switch(Ref->RPNHeadExt.grref[0]) {
      case 'N':
         f77name(cigaxg)(Ref->RPNHeadExt.grref,&pi,&pj,&d60,&dgrw,&Ref->RPNHeadExt.igref1,&Ref->RPNHeadExt.igref2,&Ref->RPNHeadExt.igref3,&Ref->RPNHeadExt.igref4,1);
         hemi=NORTH;
         f77name(ez8_vxyfll)(X,Y,Lat,Lon,&Nb,&d60,&dgrw,&pi,&pj,&hemi);
         break;

      case 'S':
         f77name(cigaxg)(Ref->RPNHeadExt.grref,&pi,&pj,&d60,&dgrw,&Ref->RPNHeadExt.igref1,&Ref->RPNHeadExt.igref2,&Ref->RPNHeadExt.igref3,&Ref->RPNHeadExt.igref4,1);
         hemi=SOUTH;
         f77name(ez8_vxyfll)(X,Y,Lat,Lon,&Nb,&d60,&dgrw,&pi,&pj,&hemi);
         break;

      case 'L':
         f77name(cigaxg)(Ref->RPNHeadExt.grref,&xlat0,&xlon0,&dlat,&dlon,&Ref->RPNHeadExt.igref1,&Ref->RPNHeadExt.igref2,&Ref->RPNHeadExt.igref3,&Ref->RPNHeadExt.igref4,1);
         GeoRef_LL2GD(X,Y,Lat,Lon,Nb,xlat0,xlon0,dlat,dlon,(Ref->AX[0]<0.0)?-180.0:0.0);
         for(i=0;i<Nb;i++) {
            X[i]-=1.0;
            Y[i]-=1.0;
         }
         break;

      case 'E':
         f77name(cigaxg)(Ref->RPNHeadExt.grref,&xlat1,&xlon1,&xlat2,&xlon2,&Ref->RPNHeadExt.igref1,&Ref->RPNHeadExt.igref2,&Ref->RPNHeadExt.igref3,&Ref->RPNHeadExt.igref4,1);
         GeoRef_RotateXY(Lat,Lon,X,Y,Nb,xlat1,xlon1,xlat2,xlon2);
         break;

      case 'W':
         GeoRef_LL2XY_W(Ref,X,Y,Lat,Lon,Nb);
         break;
   }
   return(1);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon
 * @author Jean-Philippe Gauthier
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 *    @param[in]  X       X array
 *    @param[in]  Y       Y array
 *    @param[in]  Nb      Number of coordinates
 *    @param[in]  Extrap  Enable extrapolation
 
 *    @return             Error code (0=ok)
*/
int GeoRef_XY2LL(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb,int Extrap) {

   TGeoRef *yin_gd,*yan_gd;
   int      i,j,icode;
   double   *latyin,*lonyin,*latyan,*lonyan;
   double   *tmpy;
   
   if (!Ref->XY2LL) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid transform function (XY2LL): grtyp=%c\n",__func__,Ref->GRTYP[0]);
   }

   if (Ref->NbSub > 0) {
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      j=0;
      tmpy = (double*)malloc(5*Nb*sizeof(double));
      latyin = &tmpy[j+=Nb];
      lonyin = &tmpy[j+=Nb];
      latyan = &tmpy[j+=Nb];
      lonyan = &tmpy[j+=Nb];
      for (j=0; j< Nb; j++) {
         if (Y[j] > yin_gd->NY) {
            tmpy[j]=Y[j]-yin_gd->NY;
         } else {
            tmpy[j]=Y[j];
         }
      }
      icode = GeoRef_XY2LL(yin_gd,latyin,lonyin,X,tmpy,Nb,Extrap);
      icode = GeoRef_XY2LL(yan_gd,latyan,lonyan,X,tmpy,Nb,Extrap);
      for (j=0; j < Nb; j++) {
         if (Y[j] > yin_gd->NY) {
            Lat[j]=latyan[j];
            Lon[j]=lonyan[j];
         } else {
            Lat[j]=latyin[j];
            Lon[j]=lonyin[j];
         }
      }
      free(tmpy);
   } else {
      if (Ref->Type&GRID_CORNER) {
         for(j=0;j<Nb;j++) {
            X[j]+=0.5;
            Y[j]+=0.5;
         }
      }
   
//TODO: Confirm CIndex
      if (Ref->Options.CIndex) {
         for(i=0;i<Nb;i++) {
            X[i]+=1.0;
            Y[i]+=1.0;
         }
      }

      Ref->XY2LL(Ref,Lat,Lon,X,Y,Nb);

      // Adjust for Longitude reference
      for (i=0; i<Nb; i++) {
         Lon[i]=(CLAMPLONREF(Lon[i],Ref->Options.LonRef));
      }

      icode=Nb;
      if (!Extrap) {
         for(i=0; i < Nb; i++) {
            if (X[i]<(Ref->X0-0.5) || Y[i]<(Ref->Y0-0.5) || X[i]>(Ref->X1+0.5) || Y[i]>(Ref->Y1+0.5)) {
               Lat[i]=-999.0;
               Lon[i]=-999.0;
               icode--;
            }
         }
      }
   }

   return(icode);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms LatLon coordinates XY
 * @author Jean-Philippe Gauthier
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] X       X array
 *    @param[out] Y       Y array
 *    @param[in]  Lat     Latitude array
 *    @param[in]  Lon     Longitude array
 *    @param[in]  Nb      Number of coordinates
 *    @param[in]  Extrap  Enable extrapolation
  
 *    @return             Error code (0=ok)
*/
int GeoRef_LL2XY(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb,int Extrap) {

   TGeoRef *yin_gd,*yan_gd;
   int     j,icode,maxni,maxnj;
   double  *xyin,*xyan,*yyin,*yyan;

   if (!Ref->LL2XY) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid transform function (LL2XY): grtyp=%c\n",__func__,Ref->GRTYP[0]);
   }

   if (Ref->NbSub > 0 ) {
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      maxni= yin_gd->NX;
      maxnj= yin_gd->NY;
      xyin = (double*) malloc(4*Nb*sizeof(float));
      xyan = &xyin[Nb];
      yyin = &xyin[Nb*2];
      yyan = &xyin[Nb*3];
      icode = GeoRef_LL2XY(yin_gd,xyin,yyin,Lat,Lon,Nb,Extrap);
      icode = GeoRef_LL2XY(yan_gd,xyan,yyan,Lat,Lon,Nb,Extrap);
      for (j=0; j < Nb; j++) {
         if (xyin[j] > maxni || xyin[j] < 1 || yyin[j] > maxnj || yyin[j] < 1) {
            // point is no good, take from YAN eventhough it may not be good
            X[j]=xyan[j];
            Y[j]=yyan[j]+maxnj;
         } else {
            X[j]=xyin[j];
            Y[j]=yyin[j];
         }
         if (xyan[j] >= yan_gd->mymaskgridi0 && xyan[j] <= yan_gd->mymaskgridi1 &&
            yyan[j] >= yan_gd->mymaskgridj0 && yyan[j] <= yan_gd->mymaskgridj1) {
             X[j]=xyan[j];
             Y[j]=yyan[j]+maxnj;
         }
         if (xyin[j] >= yin_gd->mymaskgridi0 && xyin[j] <= yin_gd->mymaskgridi1 &&
            yyin[j] >= yin_gd->mymaskgridj0 && yyin[j] <= yin_gd->mymaskgridj1) {
             X[j]=xyin[j];
             Y[j]=yyin[j];
         }
      }
      free(xyin);

   } else {
      Ref->LL2XY(Ref,X,Y,Lat,Lon,Nb);
   }

   // Check for grid insidness or extrapolation enabled
   icode=Nb;
   if (!Extrap) {
      for(j=0;j<Nb;j++) {
         if (X[j]>(Ref->X1+0.5) || Y[j]>(Ref->Y1+0.5) || X[j]<(Ref->X0-0.5) || Y[j]<(Ref->Y0-0.5)) {
            X[j]=-1.0;
            Y[j]=-1.0;
            icode--;
         }
      }
   }

   if (Ref->Type&GRID_CORNER) {
      for(j=0;j<Nb;j++) {
         if (X[j]!=-1.0) { 
            X[j]-=0.5;
            Y[j]-=0.5;
         }
      }
   }

   if (Ref->Options.CIndex) {
      for(j=0;j<Nb;j++) {
         if (X[j]!=-1.0) { 
            X[j]-=1.0;
            Y[j]-=1.0;
         }
      }
   }
   return(icode);
}
