/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "App.h"
#include "RPN.h"
#include "GeoRef.h"
//int GeoRef_RPNProject(TGeoRef *Ref,double X,double Y,double *Lat,double *Lon,int Extrap,int Transform) {

int GeoRef_XY2LL(TGeoRef *Ref,float *Lat,float *Lon,float *X,float *Y,int N) {

   TGeoRef *yin_gd,*yan_gd;
   int      j,icode;
   float   *latyin,*lonyin,*latyan,*lonyan,*tmpy;

   if (Ref->NbSub > 0) {
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      j=0;
      tmpy = (float*)malloc(5*N*sizeof(float));
      latyin = &tmpy[j+=N];
      lonyin = &tmpy[j+=N];
      latyan = &tmpy[j+=N];
      lonyan = &tmpy[j+=N];
      for (j=0; j< N; j++) {
         if (Y[j] > yin_gd->NY) {
            tmpy[j]=Y[j]-yin_gd->NY;
         } else {
            tmpy[j]=Y[j];
         }
      }
      icode = GeoRef_XY2LLN(yin_gd,latyin,lonyin,X,tmpy,N,FALSE);
      icode = GeoRef_XY2LLN(yan_gd,latyan,lonyan,X,tmpy,N,FALSE);
      for (j=0; j < N; j++) {
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
      icode = GeoRef_XY2LLN(Ref,Lat,Lon,X,Y,N,TRUE);
   }
   return(icode);
}

void  c_ezgfllfxy2(float *lat,float *lon,wordint *npts,float *xlat1,float *xlon1,float *xlat2,float *xlon2) {

   double *cart,*carot,latr,lonr,cosdar;
   int i,k,n,c;
   float r[3][3],ri[3][3];

   cart  = (double*)malloc(6* *npts*sizeof(double));
   carot = &cart[3**npts];

   f77name(ez_crot)(r,ri,xlon1,xlat1,xlon2,xlat2);

   for(n=0;n<*npts;n++) {
      latr=DEG2RAD(lat[n]);
      lonr=DEG2RAD(lon[n]);
      cosdar = cos(latr);
      c=n*3;
      cart[c]   = cosdar*cos(lonr);
      cart[c+1] = cosdar*sin(lonr);
      cart[c+2] = sin(latr);

      carot[c]   = ri[0][0]*cart[c+0]+ri[0][1]*cart[c+1]+ri[0][2]*cart[c+2];
      carot[c+1] = ri[1][0]*cart[c+0]+ri[1][1]*cart[c+1]+ri[1][2]*cart[c+2];
      carot[c+2] = ri[2][0]*cart[c+0]+ri[2][1]*cart[c+1]+ri[2][2]*cart[c+2];
      
      lat[n]=RAD2DEG(asin(fmax(-1.0,fmin(1.0,carot[c+2]))));
      lon[n]=RAD2DEG(atan2(carot[c+1],carot[c]));
      lon[n]=fmod(lon[n],360.0);
      if (lon[n]<0.0) lon[n]=lon[n]+360.0;
   }
   free(cart);
}

int GeoRef_XY2LLN(TGeoRef *Ref,float *lat,float *lon,float *x,float *y,int n,int Inv) {

   float xlat1, xlon1, xlat2, xlon2;
   int i, npts, un;
   float *tmpx, *tmpy, *ytmp=NULL;
   float delxx, delyy;
   float dlat, dlon, swlat, swlon;
   int indx, indy;

   npts = n;
   un = 1;

   switch(Ref->GRTYP[0]) {
      case 'A':
         for (i=0; i < n; i++) {
            lat[i] = (y[i]-1.0)*Ref->RPNHead.XG[X_DLAT]+Ref->RPNHead.XG[X_SWLAT];
            lon[i] = (x[i]-1.0)*Ref->RPNHead.XG[X_DLON]+Ref->RPNHead.XG[X_SWLON];
            lon[i] = (float) (fmod((double) (lon[i] + 360.0), (double) 360.0));
         }
         break;

      case 'B':
         for (i=0; i < n; i++) {
            lat[i] = (y[i]-1.0)*Ref->RPNHead.XG[X_DLAT]+Ref->RPNHead.XG[X_SWLAT];
            lon[i] = (x[i]-1.0)*Ref->RPNHead.XG[X_DLON]+Ref->RPNHead.XG[X_SWLON];
            lon[i] = (float) (fmod((double) (lon[i] + 360.0), (double) 360.0));
         }
         break;

      case 'E':
          for (i=0; i < n; i++) {
             dlat  = 180.0 / Ref->NY;
             dlon  = 360.0 / (Ref->NX - 1);
             swlat = -90.0 + 0.5 * dlat;
             swlon = 0.0;
             lon[i] = (x[i]-1.0)*dlon+swlon;
             lat[i] = (y[i]-1.0)*dlat+swlat;
          }

          c_ezgfllfxy2(lat,lon,&n,&Ref->RPNHead.XG[X_LAT1],&Ref->RPNHead.XG[X_LON1],&Ref->RPNHead.XG[X_LAT2],&Ref->RPNHead.XG[X_LON2]);
          free(tmpx);
          free(tmpy);
          break;


      case 'L':
         for (i=0; i < n; i++) {
            lon[i] = (x[i]-1.0)*Ref->RPNHead.XG[X_DLON]+Ref->RPNHead.XG[X_SWLON];
            lon[i] = (float) (fmod((double) (lon[i] + 360.0), (double) 360.0));
            lat[i] = (y[i]-1.0)*Ref->RPNHead.XG[X_DLAT]+Ref->RPNHead.XG[X_SWLAT];
         }
         break;

      case 'N':
      case 'S':
         f77name(ez_vllfxy)(lat,lon,x,y,&npts,&un,&Ref->RPNHead.XG[X_D60],&Ref->RPNHead.XG[X_DGRW],&Ref->RPNHead.XG[X_PI],&Ref->RPNHead.XG[X_PJ],&Ref->Hemi);
         for (i=0; i < n; i++) {
            lon[i] = (float) (fmod((double) (lon[i] + 360.0), (double) 360.0));
         }
         break;

      case 'T':
         f77name(ez_vtllfxy)(lat,lon,x,y, &Ref->RPNHead.XG[X_CLAT], &Ref->RPNHead.XG[X_CLON], &Ref->RPNHead.XG[X_TD60], &Ref->RPNHead.XG[X_TDGRW], &Ref->NX, &Ref->NY, &npts);
         break;

      case '!':
         f77name(ez_llflamb)(lat,lon,x,y,&npts,&Ref->GRTYP,&Ref->RPNHead.IG[X_IG1], &Ref->RPNHead.IG[X_IG2], &Ref->RPNHead.IG[X_IG3], &Ref->RPNHead.IG[X_IG4],1);
         break;

      case 'Y':
         App_Log(ERROR,"%s: This operation is not supported for 'Y' grids\n",__func__);
         break;

      case '#':
      case 'Z':
      case 'G':
         tmpx = (float *) malloc(n*sizeof(float));
         tmpy = (float *) malloc(n*sizeof(float));
         ytmp = (float *) malloc(n*sizeof(float));
         for (i=0; i < n; i++) {
            indx = (int)x[i]-1;
            if (Inv) {
               ytmp[i] = y[i];
               if (Ref->RPNHead.IG[X_IG2] == 1) {
                  ytmp[i] = Ref->NY +1.0 - y[i];
               }
               indy = (int)ytmp[i]-1;
            } else {
               indy = (int)y[i]-1;
            }
            indx = indx < 0 ? 0 : indx;
            indy = indy < 0 ? 0 : indy;
            indx = indx > Ref->NX-2 ? Ref->NX-2 : indx;
            indy = indy > Ref->j2-2 ? Ref->j2-2 : indy;
            delxx = Ref->AX[indx+1]-Ref->AX[indx];
            tmpx[i] = Ref->AX[indx] + ((x[i]-1.0-indx)*delxx);

            delyy = Ref->AY[indy+1]-Ref->AY[indy];
            if (Inv) {
               tmpy[i] = Ref->AY[indy] + ((ytmp[i]-1.0-indy)*delyy);
            } else {
     	         tmpy[i] = Ref->AY[indy] + ((y[i]-1.0-indy)*delyy);
            }
         }

         switch (Ref->RPNHead.GRREF[0]) {
            case 'E':
               f77name(cigaxg)(Ref->RPNHead.GRREF,&xlat1,&xlon1,&xlat2,&xlon2,&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4]);
               f77name(ez_gfllfxy)(lon, lat, tmpx, tmpy, &npts, &Ref->RPNHead.XGREF[X_LAT1], &Ref->RPNHead.XGREF[X_LON1],&Ref->RPNHead.XGREF[X_LAT2], &Ref->RPNHead.XGREF[X_LON2]);
               break;

            case 'S':
            case 'N':
               f77name(ez_vllfxy)(lat,lon,tmpx,tmpy,&npts,&un,&Ref->RPNHead.XGREF[X_D60],&Ref->RPNHead.XGREF[X_DGRW],&Ref->RPNHead.XGREF[X_PI],&Ref->RPNHead.XGREF[X_PJ],&Ref->Hemi);
               for (i=0; i < n; i++) {
                  lon[i] = (float) (fmod((double) (lon[i] + 360.0), (double) 360.0));
               }
               break;

            case 'L':
               for (i=0; i < n; i++) {
                  lat[i] = (tmpy[i])*Ref->RPNHead.XGREF[X_DLAT]+Ref->RPNHead.XGREF[X_SWLAT];
                  lon[i] = (tmpx[i])*Ref->RPNHead.XGREF[X_DLON]+Ref->RPNHead.XGREF[X_SWLON];
                  lon[i] = (float) (fmod((double) (lon[i] + 360.0), (double) 360.0));
               }
               break;

             default:
                App_Log(ERROR,"%s: Undefined reference grid type: %s\n",__func__,Ref->RPNHead.GRREF[0]);
                break;
         }
         free(tmpx);
         free(tmpy);
         free(ytmp);
         break;
   }
   return(0);
}

//TODO: not used
int c_gdllfxyz(TGeoRef* Ref,float *lat,float *lon,float *x,float *y,int n) {

   int i,npts, hem, un;
   npts = n;
  
   switch(Ref->GRTYP[0]) {
      case 'A':
      case 'B':
      case 'G':
      case 'L':
      case 'N':
      case 'S':
      case 'T':
      case '!':
         GeoRef_XY2LLN(Ref, lat, lon, x, y, n,FALSE);
      break;
      
      case 'Y':
         App_Log(ERROR,"%s: This operation is not supported for 'Y' grids\n",__func__);
         break;
      
      case '#':
      case 'Z':
         switch (Ref->RPNHead.GRREF[0]) {
            case 'E':
               f77name(ez_gfllfxy)(lon,lat,x,y,&npts,&Ref->RPNHead.XGREF[X_LAT1],&Ref->RPNHead.XGREF[X_LON1],&Ref->RPNHead.XGREF[X_LAT2],&Ref->RPNHead.XGREF[X_LON2]);
               break;
              
            case 'S':
            case 'N':
               if (Ref->RPNHead.GRREF[0] == 'N') {
                   hem = 1;
               } else {
                  hem = 2;
               }
               un = 1;
               f77name(ez_vllfxy)(lat,lon,x,y,&npts,&un,&Ref->RPNHead.XGREF[X_D60],&Ref->RPNHead.XGREF[X_DGRW], &Ref->RPNHead.XGREF[X_PI], &Ref->RPNHead.XGREF[X_PJ],&Ref->Hemi);
               break;
              
            case 'L':
               for (i=0; i < n; i++) {
                  lat[i] = (y[i])*Ref->RPNHead.XGREF[X_DLAT]+ Ref->RPNHead.XGREF[X_SWLAT];
                  lon[i] = (x[i])*Ref->RPNHead.XGREF[X_DLON]+ Ref->RPNHead.XGREF[X_SWLON];
                  lon[i] = lon[i] < 0.0 ? lon[i] + 360.0 : lon[i];
               }
               break;
              
            default:
               App_Log(ERROR,"%s: Undefined reference grid type: %s\n",__func__,Ref->RPNHead.GRREF[0]);
               break;
         }
         break;
   }
  
   return(0); 
}

int GeoRef_LL2XY(TGeoRef *Ref,float *x,float *y,float *lat,float *lon,int n) {

   TGeoRef *yin_gd, *yan_gd;
   int j, icode, maxni,maxnj ;
   float *xyin, *xyan, *yyin, *yyan;

   if (Ref->NbSub > 0 ) {
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      maxni= yin_gd->NX;
      maxnj= yin_gd->NY;
      xyin = (float *) malloc(n*sizeof(float));
      xyan = (float *) malloc(n*sizeof(float));
      yyin = (float *) malloc(n*sizeof(float));
      yyan = (float *) malloc(n*sizeof(float));
      icode = GeoRef_LL2XYN(yin_gd,xyin,yyin,lat,lon,n,FALSE);
      icode = GeoRef_LL2XYN(yan_gd,xyan,yyan,lat,lon,n,FALSE);
      for (j=0; j < n; j++) {
         if (xyin[j] > maxni || xyin[j] < 1 || yyin[j] > maxnj || yyin[j] < 1) {
            // point is no good, take from YAN eventhough it may not be good
            x[j]=xyan[j];
            y[j]=yyan[j]+maxnj;
         } else {
            x[j]=xyin[j];
            y[j]=yyin[j];
         }
         if (xyan[j] >= yan_gd->mymaskgridi0 && xyan[j] <= yan_gd->mymaskgridi1 &&
            yyan[j] >= yan_gd->mymaskgridj0 && yyan[j] <= yan_gd->mymaskgridj1) {
             x[j]=xyan[j];
             y[j]=yyan[j]+maxnj;
         }
         if (xyin[j] >= yin_gd->mymaskgridi0 && xyin[j] <= yin_gd->mymaskgridi1 &&
            yyin[j] >= yin_gd->mymaskgridj0 && yyin[j] <= yin_gd->mymaskgridj1) {
             x[j]=xyin[j];
             y[j]=yyin[j];
         }
      }
      free(xyin);free(xyan);free(yyin);free(yyan);
   } else {
     icode = GeoRef_LL2XYN(Ref,x,y,lat,lon,n,TRUE);
   }
   return(icode);
}

int GeoRef_LL2XYN(TGeoRef *Ref,float *x,float *y,float *lat,float *lon,int n,int Inv) {

   float *tmplons;
   int    j,ni_in, nj_in;
   int    sym=Ref->Options.Symmetric;
   int    npts;
   int    coordonnee;

	 npts = n;
	 ni_in = Ref->NX;
	 nj_in = Ref->NY;

	 switch(Ref->GRTYP[0]) {
	    case 'A':
	    case 'B':
	    case 'E':
	    case 'L':
	    case 'N':
	    case 'S':
	    case 'T':
	    case '!':
	       tmplons = (float *)malloc(npts * sizeof(float));
	       memcpy(tmplons,lon,sizeof(float)*npts);
	      
	       f77name(ez_ll2rgd)(x,y,lat,tmplons,&npts,&ni_in,&nj_in,&Ref->GRTYP,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4],&sym,Ref->AY);
	       free(tmplons);
	       break;

	    case '#':
	    case 'Z':
	    case 'G':
	       coordonnee = RELATIVE;
	       nj_in =  Ref->j2;
	       f77name(ez_ll2igd)(x,y,lat,lon,&npts,&ni_in,&nj_in,&Ref->GRTYP,&Ref->RPNHead.GRREF,&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4],Ref->AX,Ref->AY,&coordonnee);
	       if (Ref->GRTYP[0] == 'G' && Ref->RPNHead.IG[X_IG1] == 1) {
	          for (j=0; j < npts; j++) {
	             y[j] = y[j] - nj_in;
	          }
	       }
	       if (Inv && Ref->GRTYP[0] == 'G' && Ref->RPNHead.IG[X_IG2] == 1)  {
	          for (j=0; j < npts; j++) {
	             y[j] = nj_in +1.0 - y[j];
	          }
	       }
         break;
     
      default:
         break;
    }

   return(0);
}

//TODO: not used
int c_gdxyzfll(TGeoRef *Ref, float *x, float *y, float *lat, float *lon, int n) {

   int ni_in, nj_in;
   int coordonnee;
   int npts;

   npts = n;
              
   ni_in =  Ref->NX;
   nj_in =  Ref->NY;
   
   switch(Ref->GRTYP[0]) {
      case 'A':
      case 'B':
      case 'E':
      case 'G':
      case 'L':
      case 'N':
      case 'S':
      case 'T':
      case '!':
	       GeoRef_LL2XYN(Ref, x, y, lat, lon, n,FALSE);
         break;
        
      case 'Y':
         App_Log(ERROR,"%s: This operation is not supported for 'Y' grids\n",__func__);
         break;
	
      case '#':
      case 'Z':
	       coordonnee = ABSOLUTE;
         f77name(ez_ll2igd)(x,y,lat,lon,&npts,&ni_in,&nj_in,&Ref->GRTYP,&Ref->RPNHead.GRREF,&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4],Ref->AX,Ref->AY,&coordonnee);
         break;

      default:
        break;
   }

   return(0);
}