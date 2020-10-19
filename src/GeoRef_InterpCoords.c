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
#include "Vertex.h"

static inline void c_ezll2gd(double *X,double *Y,double *Lat,double *Lon,int Nb,float Lat0,float Lon0,float DLat,float DLon,float LonRef) {

   int    i;
   
   for(i=0;i<Nb;i++) {
      X[i] = (CLAMPLONREF(Lon[i],LonRef)-Lon0)/DLon + 1.0;
      Y[i] = (Lat[i]-Lat0)/DLat + 1.0;
   }
}

void  c_ezgfxyfll(double *Lat,double *Lon,double *X,double *Y,int npts,float xlat1,float xlon1,float xlat2,float xlon2) {

   double *cart,*carot,latr,lonr,cosdar;
   int i,k,n,c;
   float r[3][3],ri[3][3];

   cart  = (double*)malloc(6*npts*sizeof(double));
   carot = &cart[3*npts];

   f77name(ez_crot)(r,ri,&xlon1,&xlat1,&xlon2,&xlat2);

   for(n=0;n<npts;n++) {
      latr=DEG2RAD(Lat[n]);
      lonr=DEG2RAD(Lon[n]);
      cosdar = cos(latr);
      c=n*3;
      cart[c]   = cosdar*cos(lonr);
      cart[c+1] = cosdar*sin(lonr);
      cart[c+2] = sin(latr);

      carot[c]   = r[0][0]*cart[c+0]+r[0][1]*cart[c+1]+r[0][2]*cart[c+2];
      carot[c+1] = r[1][0]*cart[c+0]+r[1][1]*cart[c+1]+r[1][2]*cart[c+2];
      carot[c+2] = r[2][0]*cart[c+0]+r[2][1]*cart[c+1]+r[2][2]*cart[c+2];
      
      Y[n]=RAD2DEG(asin(fmax(-1.0,fmin(1.0,carot[c+2]))));
      X[n]=RAD2DEG(atan2(carot[c+1],carot[c]));
      X[n]=fmod(X[n],360.0);
      if (X[n]<0.0) X[n]+=360.0;
   }
   free(cart);
}
void  c_ezgfllfxy(double *Lat,double *Lon,double *X,double *Y,int npts,float xlat1,float xlon1,float xlat2,float xlon2) {

   double *cart,*carot,latr,lonr,cosdar;
   int i,k,n,c;
   float r[3][3],ri[3][3];

   cart  = (double*)malloc(6*npts*sizeof(double));
   carot = &cart[3*npts];

   f77name(ez_crot)(r,ri,&xlon1,&xlat1,&xlon2,&xlat2);

   for(n=0;n<npts;n++) {
      latr=DEG2RAD(Y[n]);
      lonr=DEG2RAD(X[n]);
      cosdar = cos(latr);
      c=n*3;
      cart[c]   = cosdar*cos(lonr);
      cart[c+1] = cosdar*sin(lonr);
      cart[c+2] = sin(latr);

      carot[c]   = ri[0][0]*cart[c+0]+ri[0][1]*cart[c+1]+ri[0][2]*cart[c+2];
      carot[c+1] = ri[1][0]*cart[c+0]+ri[1][1]*cart[c+1]+ri[1][2]*cart[c+2];
      carot[c+2] = ri[2][0]*cart[c+0]+ri[2][1]*cart[c+1]+ri[2][2]*cart[c+2];
      
      Lat[n]=RAD2DEG(asin(fmax(-1.0,fmin(1.0,carot[c+2]))));
      Lon[n]=RAD2DEG(atan2(carot[c+1],carot[c]));
   }
   free(cart);
}


int GeoRef_XY2LL_R(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {

   TCoord loc0;
   double x,y,d;
   int    n;

   for(n=0;n<Nb;n++) {

      loc0.Lat=DEG2RAD(Ref->Loc.Lat);
      loc0.Lon=DEG2RAD(Ref->Loc.Lon);

      x=X[n]*Ref->ResA;
      y=Y[n]*Ref->ResR;

      x=DEG2RAD(x);
      d=M2RAD(y*Ref->CTH);

      if (Ref->Options.Transform) {
         Lat[n]=asin(sin(loc0.Lat)*cos(d)+cos(loc0.Lat)*sin(d)*cos(x));
         Lon[n]=fmod(loc0.Lon+(atan2(sin(x)*sin(d)*cos(loc0.Lat),cos(d)-sin(loc0.Lat)*sin(*Lat)))+M_PI,M_2PI)-M_PI;
         Lat[n]=RAD2DEG(*Lat);
         Lon[n]=RAD2DEG(*Lon);
      } else {
         Lat[n]=d;
         Lon[n]=x;
      }
   }

   return(0);
}

static inline int GeoRef_WKTRotate(TRotationTransform *T,double *Lat,double *Lon) {

   double lat,lon,x,y,z,xr,yr,zr;
   
   lon = DEG2RAD(*Lon);
   lat = DEG2RAD(*Lat);

   // Convert from spherical to cartesian coordinates
   x = cos(lon)*cos(lat); 
   y = sin(lon)*cos(lat);
   z = sin(lat);

   xr = T->CosTheta*T->CosPhi*x + T->CosTheta*T->SinPhi*y + T->SinTheta*z;
   yr = -T->SinPhi*x + T->CosPhi*y;
   zr = -T->SinTheta*T->CosPhi*x - T->SinTheta*T->SinPhi*y + T->CosTheta*z;
      
   // Convert cartesian back to spherical coordinates
   *Lon = RAD2DEG(atan2(yr,xr)); 
   *Lat = RAD2DEG(asin(zr));
   
   return(TRUE);
}

static inline int GeoRef_WKTUnRotate(TRotationTransform *T,double *Lat,double *Lon) {

   double lat,lon,x,y,z,xr,yr,zr;
   
   lon = DEG2RAD(*Lon);
   lat = DEG2RAD(*Lat);

   // Convert from spherical to cartesian coordinates
   x = cos(lon)*cos(lat); 
   y = sin(lon)*cos(lat);
   z = sin(lat);

   xr = T->CosTheta*T->CosPhi*x + -T->SinPhi*y + -T->SinTheta*T->CosPhi*z;
   yr = -T->CosTheta*-T->SinPhi*x + T->CosPhi*y - -T->SinTheta*-T->SinPhi*z;
   zr = T->SinTheta*x + T->CosTheta*z;
      
   // Convert cartesian back to spherical coordinates
   *Lon = RAD2DEG(atan2(yr,xr)); 
   *Lat = RAD2DEG(asin(zr));

   return(TRUE);
}

int GeoRef_XY2LL_W(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {

#ifdef HAVE_GDAL
   double x,y,z=0.0,d;
   int    n,ok;

   for(n=0;n<Nb;n++) {

      // Transform the point into georeferenced coordinates 
      x=X[n];
      y=Y[n];
      if (Ref->Options.Transform) {
         if (Ref->Transform) {
            x=Ref->Transform[0]+Ref->Transform[1]*X[n]+Ref->Transform[2]*Y[n];
            y=Ref->Transform[3]+Ref->Transform[4]*X[n]+Ref->Transform[5]*Y[n];
         } else if (Ref->GCPTransform) {
            GDALGCPTransform(Ref->GCPTransform,FALSE,1,&x,&y,&z,&ok);
         } else if (Ref->TPSTransform) {
            GDALGCPTransform(Ref->TPSTransform,FALSE,1,&x,&y,&z,&ok);
         } else if (Ref->RPCTransform) {
            GDALGCPTransform(Ref->RPCTransform,FALSE,1,&x,&y,&z,&ok);
         }
      }
      
      // Transform to latlon
      if (Ref->Function) {
         if (!OCTTransform(Ref->Function,1,&x,&y,NULL)) {
            *Lon=-999.0;
            *Lat=-999.0;
         }
      }

      Lon[n]=x;
      Lat[n]=y;
      
      if (Ref->RotTransform) 
         GeoRef_WKTUnRotate(Ref->RotTransform,&Lat[n],&Lon[n]);
   }
#else
   App_Log(ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
#endif

   return(0);
}

int GeoRef_XY2LL(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {

   TGeoRef *yin_gd,*yan_gd;
   int      j,icode;
   double   *latyin,*lonyin,*latyan,*lonyan;
   float    xlat1, xlon1, xlat2, xlon2;
   int      i, npts, un;
   double   *tmpx, *tmpy, *ytmp=NULL;
   float    delxx, delyy;
   float    dlat, dlon, swlat, swlon;
   int      s,sx,sy,indx, indy;
   double   dx,dy;

   npts = Nb;
   un = 1;

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
      icode = GeoRef_XY2LL(yin_gd,latyin,lonyin,X,tmpy,Nb);
      icode = GeoRef_XY2LL(yan_gd,latyan,lonyan,X,tmpy,Nb);
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

      if (!Ref->Options.ExtrapDegree) {
         for(i=0; i < Nb; i++) {
            if (X[i]<(Ref->X0-0.5) || Y[i]<(Ref->Y0-0.5) || X[i]>(Ref->X1+0.5) || Y[i]>(Ref->Y1+0.5)) {
               Lat[i]=-999.0;
               Lon[i]=-999.0;
            }
         }
      }
      if (Ref->Options.CIndex) {
         for(i=0;i<Nb;i++) {
            X[i]+=1.0;
            Y[i]+=1.0;
         }
      }

      switch(Ref->GRTYP[0]) {
         case 'A':
         case 'B':
            for (i=0; i < Nb; i++) {
               Lat[i] = (Y[i]-1.0)*Ref->RPNHead.XG[X_DLAT]+Ref->RPNHead.XG[X_SWLAT];
               Lon[i] = (X[i]-1.0)*Ref->RPNHead.XG[X_DLON]+Ref->RPNHead.XG[X_SWLON];
            }
            break;

         case 'E':
            for (i=0; i < Nb; i++) {
               dlat  = 180.0 / Ref->NY;
               dlon  = 360.0 / (Ref->NX - 1);
               swlat = -90.0 + 0.5 * dlat;
               swlon = 0.0;
               Lon[i] = (X[i]-1.0)*dlon+swlon;
               Lat[i] = (Y[i]-1.0)*dlat+swlat;
            }

            c_ezgfllfxy(Lat,Lon,Lon,Lat,Nb,Ref->RPNHead.XG[X_LAT1],Ref->RPNHead.XG[X_LON1],Ref->RPNHead.XG[X_LAT2],Ref->RPNHead.XG[X_LON2]);
            break;

         case 'L':
            for (i=0; i < Nb; i++) {
               Lon[i] = (X[i]-1.0)*Ref->RPNHead.XG[X_DLON]+Ref->RPNHead.XG[X_SWLON];
               Lat[i] = (Y[i]-1.0)*Ref->RPNHead.XG[X_DLAT]+Ref->RPNHead.XG[X_SWLAT];
            }
            break;

         case 'N':
         case 'S':
            f77name(ez8_vllfxy)(Lat,Lon,X,Y,&npts,&un,&Ref->RPNHead.XG[X_D60],&Ref->RPNHead.XG[X_DGRW],&Ref->RPNHead.XG[X_PI],&Ref->RPNHead.XG[X_PJ],&Ref->Hemi);
            break;

         case 'O':
            for (i=0; i < Nb; i++) {
               Lon[i]=Vertex_ValS(Ref->AX,NULL,Ref->NX,Ref->NY,X[i],Y[i],TRUE);
               Lat[i]=Vertex_ValS(Ref->AY,NULL,Ref->NX,Ref->NY,X[i],Y[i],FALSE);
            }
            break;

         case 'R':
             GeoRef_XY2LL_R(Ref,Lat,Lon,X,Y,Nb);
             break;

         case 'W':
             GeoRef_XY2LL_W(Ref,Lat,Lon,X,Y,Nb);
             break;

         case 'Y':
            switch (Ref->RPNHead.GRREF[0]) {
               case 'L':
                  for (i=0; i < Nb; i++) {
                     indx=ROUND(Y[i])*(Ref->X1-Ref->X0)+ROUND(X[i]);
                     Lat[i]=Ref->AY[indx];
                     Lon[i]=Ref->AX[indx];
                  }
                  break;

               case 'W':
                  tmpx = (double*)malloc(2*Nb*sizeof(double));
                  tmpy = &tmpx[Nb];

                  for (i=0; i < Nb; i++) {
                     sx=floor(X[i]);sx=CLAMP(sx,Ref->X0,Ref->X1);
                     sy=floor(Y[i]);sy=CLAMP(sy,Ref->Y0,Ref->Y1);
                     dx=X[i]-sx;;
                     dy=Y[i]-sy;

                     s=sy*Ref->NX+sx;
                     tmpx[i]=Ref->AX[s];
                     tmpy[i]=Ref->AY[s];

                     if (++sx<=Ref->X1) {
                        s=sy*Ref->NX+sx;
                        tmpx[i]+=(Ref->AX[s]-tmpx[i])*dx;
                     }

                     if (++sy<=Ref->Y1) {
                        s=sy*Ref->NX+(sx-1);
                        tmpy[i]+=(Ref->AY[s]-tmpy[i])*dy;
                     }
                  }
                  GeoRef_XY2LL_W(Ref,Lat,Lon,tmpx,tmpy,Nb);
                  free(tmpx);
                  break;

               default:
                  App_Log(ERROR,"%s: Undefined reference grid type: %s\n",__func__,Ref->RPNHead.GRREF[0]);
                  break;
            }
            break;

         case 'T':
            f77name(ez8_vtllfxy)(Lat,Lon,X,Y, &Ref->RPNHead.XG[X_CLAT], &Ref->RPNHead.XG[X_CLON], &Ref->RPNHead.XG[X_TD60], &Ref->RPNHead.XG[X_TDGRW], &Ref->NX, &Ref->NY, &npts);
            break;

         case '!':
            f77name(ez_llflamb)(Lat,Lon,X,Y,&npts,&Ref->GRTYP,&Ref->RPNHead.IG[X_IG1], &Ref->RPNHead.IG[X_IG2], &Ref->RPNHead.IG[X_IG3], &Ref->RPNHead.IG[X_IG4],1);
            break;

         case '#':
         case 'Z':
         case 'G':
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
                  c_ezgfllfxy(Lat,Lon,tmpx,tmpy,npts,Ref->RPNHead.XGREF[X_LAT1],Ref->RPNHead.XGREF[X_LON1],Ref->RPNHead.XGREF[X_LAT2],Ref->RPNHead.XGREF[X_LON2]);
                  break;

               case 'S':
               case 'N':
                  f77name(ez8_vllfxy)(Lat,Lon,tmpx,tmpy,&npts,&un,&Ref->RPNHead.XGREF[X_D60],&Ref->RPNHead.XGREF[X_DGRW],&Ref->RPNHead.XGREF[X_PI],&Ref->RPNHead.XGREF[X_PJ],&Ref->Hemi);
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
                  break;
            }
            free(tmpx);
            break;
         
         default:
            App_Log(ERROR,"%s: Unsuported grid type: %s\n",__func__,Ref->GRTYP[0]);
            break;
      }
   }

   // Adjust for Longitude reference
   for (i=0; i < Nb; i++) {
      Lon[i]=(CLAMPLONREF(Lon[i],Ref->Options.LonRef));
   }

   return(0);
}


int GeoRef_LL2XY_R(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   TCoord  loc0;
   double  x,d,lat,lon;
   int     n;

   for(n=0;n<Nb;n++) {
      loc0.Lat=DEG2RAD(Ref->Loc.Lat);
      loc0.Lon=DEG2RAD(Ref->Loc.Lon);
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

int GeoRef_LL2XY_O(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   int     out=0;
   int     x,y,n,nd,d,dx,dy,idx,idxs[8];
   double  dists[8],xx,yy;
   Vect2d  pts[4],pt;

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

int GeoRef_LL2XY_RG(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   float  pi,pj,dgrw,d60,clat,clon,dellat,dellon,xlat0,xlon0,xlat1,xlon1,xlat2,xlon2;
   double tmplon;      
   int    i,indy,hemi;

   switch(Ref->GRTYP[0]) {
      case 'N':
         f77name(cigaxg)(Ref->GRTYP,&pi,&pj,&d60,&dgrw,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);
         hemi=NORTH;
         f77name(ez8_vxyfll)(X,Y,Lat,Lon,&Nb,&d60,&dgrw,&pi,&pj,&hemi);
         break;
      
      case 'S':
         hemi=SOUTH;
         f77name(cigaxg)(Ref->GRTYP,&pi,&pj,&d60,&dgrw,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);
         f77name(ez8_vxyfll)(X,Y,Lat,Lon,&Nb,&d60,&dgrw,&pi,&pj,&hemi);
         break;
      
      case 'T':
         f77name(cigaxg)(Ref->GRTYP,&d60,&dgrw,&clat,&clon,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);
         f77name(ez8_vtxyfll)(X,Y,Lat,Lon,&clat,&clon,&d60,&dgrw,&Ref->NX,&Ref->NY,&Nb);
         break;
      
      case 'A':
         dellon = 360.0 / Ref->NX;
         xlon0 = 0.0;
         if (Ref->RPNHead.IG[X_IG1]==GLOBAL) {
            dellat = 180.0 / Ref->NY;
            xlat0 = -90.0 + dellat * 0.5;
         }
         
         if (Ref->RPNHead.IG[X_IG1]==NORTH) {
            dellat = 90.0 / Ref->NY;
            xlat0 =  dellat * 0.5;
         }
         
         if (Ref->RPNHead.IG[X_IG1]==SOUTH) {
            dellat = 90.0 / Ref->NY;
            xlat0 = -90.0 + dellat * 0.5;
         }        
         
         c_ezll2gd(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);
         break;
      
      case 'B':
         dellon = 360.0 / (Ref->NX-1);
         xlon0 = 0.0;
         if (Ref->RPNHead.IG[X_IG1]==GLOBAL) {
            dellat = 180.0 / (Ref->NY-1);
            xlat0 = -90.0;
         }
         
         if (Ref->RPNHead.IG[X_IG1]==NORTH) {
            dellat = 90.0 / (Ref->NY-1);
            xlat0 =  0.0;
         }
         
         if (Ref->RPNHead.IG[X_IG1]==SOUTH) {
            dellat = 90.0 / (Ref->NY-1);
            xlat0 = -90.0;
         }
         
         c_ezll2gd(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);
         break;
      
      case 'G':
         dellon = 360.0 / Ref->NX;
         xlon0 = 0.0;

         if (Ref->RPNHead.IG[X_IG1]==GLOBAL) {                  
            for(i=0;i<Nb;i++) {
               X[i] = (CLAMPLONREF(Lon[i],Ref->Options.LonRef) - xlon0)/dellon + 1.0;
               indy = GeoRef_XFind(Lat[i],Ref->AX,Ref->NX);
               if (indy>Ref->NY) indy = Ref->NY - 2;
               
               Y[i]= indy+(Lat[i]-Ref->AX[indy])/(Ref->AX[indy+1]-Ref->AX[indy]);
            }

         } else if  (Ref->RPNHead.IG[X_IG1]==NORTH) {
            dellat = 90.0 / Ref->NY;
            xlat0 =  dellat * 0.5;
            c_ezll2gd(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);
         } else {
            dellat = 90.0 / Ref->NY;
            xlat0 = -90.0 + dellat * 0.5;
            c_ezll2gd(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);
         }
         break;
      
      case 'L':
         f77name(cigaxg)(Ref->GRTYP,&xlat0,&xlon0,&dellat,&dellon,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);         
         c_ezll2gd(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);
         break;

      case 'E':
         f77name(cigaxg)(Ref->GRTYP,&xlat1,&xlon1,&xlat2,&xlon2,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);         
         c_ezgfxyfll(Lat,Lon,X,Y,Nb,xlat1,xlon1,xlat2,xlon2);

         dellon = 360.0 / (Ref->NX-1);
         xlon0 = 0.0;

         dellat = 180.0 / (Ref->NY);
         xlat0 = -90. + 0.5*dellat;

         c_ezll2gd(X,Y,Lat,Lon,Nb,xlat0,xlon0,dellat,dellon,0.0);
         break;

      case '!':
         f77name(ez_lambfll)(X,Y,Lat,Lon,&Nb,Ref->GRTYP,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);
         break;

      default:
         App_Log(ERROR,"%s: Invalid grid type: %c\n",__func__,Ref->GRTYP[0]);   
   }

   return(0);
}

int GeoRef_LL2XY_IRG(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   int   i,j,indx,indy,hemi;
   float pi,pj,dgrw,d60,dlat,dlon,xlat0,xlon0,xlat1,xlon1,xlat2,xlon2;
      
   switch(Ref->RPNHead.GRREF[0]) {
      case 'N':
         f77name(cigaxg)(Ref->RPNHead.GRREF,&pi,&pj,&d60,&dgrw,&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4]);
         hemi=NORTH;
         f77name(ez8_vxyfll)(X,Y,Lat,Lon,&Nb,&d60,&dgrw,&pi,&pj,&hemi);
         break;

      case 'S':
         f77name(cigaxg)(Ref->RPNHead.GRREF,&pi,&pj,&d60,&dgrw,&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4]);
         hemi=SOUTH;
         f77name(ez8_vxyfll)(X,Y,Lat,Lon,&Nb,&d60,&dgrw,&pi,&pj,&hemi);
         break;

      case 'L':
         f77name(cigaxg)(Ref->RPNHead.GRREF,&xlat0,&xlon0,&dlat,&dlon,&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4]);
         c_ezll2gd(X,Y,Lat,Lon,Nb,xlat0,xlon0,dlat,dlon,(Ref->AX[0]<0.0)?-180.0:0.0);
         for(i=0;i<Nb;i++) {
            X[i]-=1.0;
            Y[i]-=1.0;
         }
         break;

      case 'E':
         f77name(cigaxg)(Ref->RPNHead.GRREF,&xlat1,&xlon1,&xlat2,&xlon2,&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4]);
         c_ezgfxyfll(Lat,Lon,X,Y,Nb,xlat1,xlon1,xlat2,xlon2);
         break;
   }
   for(i=0;i<Nb;i++) {
      indx = GeoRef_XFind(X[i],Ref->AX,Ref->NX);
      indy = GeoRef_XFind(Y[i],Ref->AY,Ref->j2);
      
      if (indx >= Ref->NX) indx = Ref->NX - 1;
      if (indy >= Ref->j2) indy = Ref->j2 - 1;
      
      X[i] = indx+(X[i]-Ref->AX[indx])/(Ref->AX[indx+1]-Ref->AX[indx]);
      Y[i] = indy+(Y[i]-Ref->AY[indy])/(Ref->AY[indy+1]-Ref->AY[indy]);
   } 

   if (Ref->GRTYP[0] == 'G') {
      if (Ref->RPNHead.IG[X_IG1] == 1) {
         for (j=0; j < Nb; j++) Y[j] = Y[j] - Ref->j2;
      }
      if (Ref->RPNHead.IG[X_IG2] == 1)  {
         for (j=0; j < Nb; j++) Y[j] = Ref->j2 +1.0 - Y[j];
      }
   }
   return(0);
}

int GeoRef_LL2XY_Y(TGeoRef *Ref,TGridSet *Set,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   double   tmpwgts[GOPT_MAXWEIGHTNUM],total_wgt;
   int      locmax,i,iz,idx,idxz,tmp_idxs[GOPT_MAXWEIGHTNUM][2],previous_val_polar_correction;
   TPoint2D bbox[2];

   previous_val_polar_correction = Ref->Options.PolarCorrect;
   Ref->Options.PolarCorrect = FALSE;
//TODO: finish this func
   Set->n_wts = MIN(Ref->Options.WeightNum,GOPT_MAXWEIGHTNUM);
   Set->wts =  (double *) malloc(Nb * Set->n_wts*sizeof(float));
   Set->idx =  (int *) calloc(Nb,Set->n_wts*sizeof(int));
   Set->mask = (int *) malloc(Nb*sizeof(int));
   GeoRef_CalcLL(Ref);

   // Find bbox
   bbox[0].X=bbox[0].Y=1e30;
   bbox[1].X=bbox[1].Y=-1e30;

   for(i=0;i<Nb;i++) {
      bbox[0].Y=FMIN(Lat[i],bbox[0].Y);
      bbox[0].X=FMIN(Lon[i],bbox[0].X);
      bbox[1].Y=FMAX(Lat[i],bbox[1].Y);
      bbox[1].X=FMAX(Lon[i],bbox[1].X);
   }

   // Initialize weights
   for(i=0;i<Nb*Set->n_wts;i++) {
      tmpwgts[i] = 1.0e30;
   }

   for (idx=0;idx<Nb;idx++) {
      X[idx]=-1.0;
      Y[idx]=-1.0;
      idxz=idx*Set->n_wts;
      locmax=1;
      for(i=0;i<GOPT_MAXWEIGHTNUM;i++) tmpwgts[i] = 1.0e30;
      memset(tmp_idxs,0x0,GOPT_MAXWEIGHTNUM*2);

      if (GeoRef_Nearest(Ref,Lon[idx],Lat[idx],&Set->idx[idxz],&Set->wts[idxz],1,Ref->Options.DistTreshold)) {
         if (Set->idx[idxz]<Ref->Options.DistTreshold) {
            if (!Ref->mask || Ref->mask[Set->idx[idxz]]==1) {
               if (Set->wts[idxz] < tmpwgts[locmax]) {
   //TODO: finalise (ez_calcxy_y)
   //               tmpwgts[locmax] = Set->wts[idxz];
   //               Set->idx(idx,locmax) = Set->idx[idxz];
   //               tmp_idxs[locmax][0] = Set->idx[idxz];
   //               tmp_idxs[locmax][1] = 1;
   //               locmax = maxloc(tmpwgts,1)
               }
            }
         }
      }
      FWITHIN(0,bbox[0].Y,bbox[0].X,bbox[1].Y,bbox[1].X,Y[idx],X[idx]);
      if (Set->mask[idx]) {
//         (f77name)inside_or_outside(&Set->mask[idx],&Set->x[idx],&Set->y[idx],&Lat[idx],&Lon[idx],Ref->Lat,Ref->Lon,&Ref->NX,&Ref->NY,tmpwgts,tmp_idxs,&Set->n_wts;
      }

      if (Set->mask[idx]) {
         for(iz=0;iz<Set->n_wts;iz++) {
            idxz=idx*Nb+iz;
            Set->wts[idxz] = tmpwgts[iz];
         }
         idxz=idx*Nb;
         if (Set->wts[idxz] > 6371000.0) {
            for(iz=0;iz>Set->n_wts;iz++) { 
               Set->wts[idxz+iz] = 1.0e30;
            }
            Set->mask[idx] = 0;
         } else {
            total_wgt = 0.0;
            for (iz=0;iz<Set->n_wts;iz++) {
               idxz=idx*Nb+iz;
               Set->wts[idxz] = fmax(Set->wts[idxz],1.0e-10);
               Set->wts[idxz] = 1.0 / Set->wts[idxz];
               total_wgt+=Set->wts[idxz];
            }
            for (iz=0;iz<Set->n_wts;iz++) {
               idxz=idx*Nb+iz;
               Set->wts[idxz]/=total_wgt;
            }
         }
      }
   }

   Ref->Options.PolarCorrect = previous_val_polar_correction;
   return(0);
}

int GeoRef_LL2XY(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   TGeoRef *yin_gd,*yan_gd;
   int     j,icode,maxni,maxnj;
   double  *xyin,*xyan,*yyin,*yyan;

   if (Ref->NbSub > 0 ) {
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      maxni= yin_gd->NX;
      maxnj= yin_gd->NY;
      xyin = (double*) malloc(4*Nb*sizeof(float));
      xyan = &xyin[Nb];
      yyin = &xyin[Nb*2];
      yyan = &xyin[Nb*3];
      icode = GeoRef_LL2XY(yin_gd,xyin,yyin,Lat,Lon,Nb);
      icode = GeoRef_LL2XY(yan_gd,xyan,yyan,Lat,Lon,Nb);
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
      switch(Ref->GRTYP[0]) {
         case 'A':
         case 'B':
         case 'E':
         case 'L':
         case 'N':
         case 'S':
         case 'T':
         case '!':
            GeoRef_LL2XY_RG(Ref,X,Y,Lat,Lon,Nb);
            break;

         case '#':
         case 'Z':
         case 'G':
            GeoRef_LL2XY_IRG(Ref,X,Y,Lat,Lon,Nb);
            break;

         case 'O':
            GeoRef_LL2XY_O(Ref,X,Y,Lat,Lon,Nb);
            break;

         case 'R':
            GeoRef_LL2XY_R(Ref,X,Y,Lat,Lon,Nb);
            break;

         case 'M':
            GeoRef_LL2XY_M(Ref,X,Y,Lat,Lon,Nb);
            break;
            
         case 'Y':
//TODO: what about gset            GeoRef_LL2XY_Y(Ref,gset,X,Y,Lat,Lon,Nb);
         break;

         default:
            App_Log(ERROR,"%s: Invalid grid type: %c\n",__func__,Ref->GRTYP[0]);
            break;
      }
   }

   // Check for grid insidness or extrapolation enabled
   if (!Ref->Options.ExtrapDegree) {
      for(j=0;j<Nb;j++) {
         if (X[j]>(Ref->X1+0.5) || Y[j]>(Ref->Y1+0.5) || X[j]<(Ref->X0-0.5) || Y[j]<(Ref->Y0-0.5)) {
            X[j]=-1.0;
            Y[j]=-1.0;
            Nb--;
         }
      }
   }

   if (Ref->Options.CIndex) {
      for(j=0;j<Nb;j++) {
         X[j]-=1.0;
         Y[j]-=1.0;
      }
   }
   return(Nb);
}
