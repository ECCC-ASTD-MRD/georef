
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
         App_Log(APP_ERROR,"%s: Undefined reference grid type: %s\n",__func__,Ref->RPNHead.GRREF[0]);
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
int GEM_grid_param(int *F_bsc_base,int *F_bsc_ext1,int *F_extension ,int F_maxcfl,float *F_lonr,float *F_latr,int *F_ni,int *F_nj,float *F_dx,float *F_dy,double *F_x0_8,double *F_y0_8,double *F_xl_8,double *F_yl_8,int F_overlap,int F_yinyang_L) {

   double delta_8;
   int iref,jref;
  
   // basic global lateral boundary conditions width
   *F_bsc_base = 5;
   if (F_yinyang_L) *F_bsc_base=*F_bsc_base+1;

   // added points for proper de-staggering of u,v at physics interface
   *F_bsc_ext1 = 2;

   // total extension to user specified grid configuration
   *F_extension= F_maxcfl + *F_bsc_base + *F_bsc_ext1;

   if (F_yinyang_L) {

      *F_x0_8 =   45.0 - 3.0*F_overlap;
      *F_xl_8 =  315.0 + 3.0*F_overlap;
      *F_y0_8 = -45.0  -     F_overlap;
      *F_yl_8 =  45.0  +     F_overlap;
      
      delta_8  = ((*F_xl_8)-(*F_x0_8))/(*F_ni-1);
      *F_dx   = delta_8;
      *F_x0_8 = *F_x0_8 - (*F_extension)*delta_8;
      *F_xl_8 = *F_xl_8 + (*F_extension)*delta_8;
      
      delta_8  = ((*F_yl_8)-(*F_y0_8))/(*F_nj-1);
      *F_dy   = delta_8;
      *F_y0_8 = *F_y0_8 - (*F_extension)*delta_8;
      *F_yl_8 = *F_yl_8 + (*F_extension)*delta_8;
      
      *F_ni   = *F_ni + 2  *(*F_extension);
      *F_nj   = *F_nj + 2  *(*F_extension);

   } else {

      iref = *F_ni / 2 + (*F_extension);
      if ((*F_ni)%2==0) {
         *F_lonr = *F_lonr - (*F_dx)/2.0;
      } else {
         iref = iref + 1;
      }
      jref = *F_nj / 2 + (*F_extension);
      if ((*F_nj)%2==0) {
         *F_latr = *F_latr - (*F_dy)/2.0;
      } else {
         jref = *F_nj / 2 + (*F_extension) + 1;
      }
      
      *F_ni   = *F_ni + 2*(*F_extension);
      *F_nj   = *F_nj + 2*(*F_extension);
      *F_x0_8 = *F_lonr - (iref-1) * (*F_dx);
      *F_y0_8 = *F_latr - (jref-1) * (*F_dy);
      *F_xl_8 = *F_x0_8 + (*F_ni  -1) * (*F_dx);
      *F_yl_8 = *F_y0_8 + (*F_nj  -1) * (*F_dy);
      if (*F_x0_8 < 0.) *F_x0_8=*F_x0_8+360.0;
      if (*F_xl_8 < 0.) *F_xl_8=*F_xl_8+360.0;

      if (*F_x0_8 < 0.) {
         fprintf(stderr,"Longitude of WEST %f < 0.0\n",*F_x0_8);
         return(0);
      }
      if (*F_y0_8 < -90.) {
         fprintf(stderr,"Latitude of SOUTH %f < 0.0\n",*F_y0_8);
         return(0);
      }
      if (*F_xl_8 > 360.) {
         fprintf(stderr,"Longitude of EAST %f < 0.0\n",*F_xl_8);
         return(0);
      }
      if (*F_yl_8 > 90.) {
         fprintf(stderr,"Latitude of NORTH %f < 0.0\n",*F_yl_8);
         return(0);
      }
   }
   return(1);
}

void GEM_hgrid4(double *F_xgi_8,double *F_ygi_8,int F_Grd_ni,int F_Grd_nj,float *F_Grd_dx,float *F_Grd_dy,double F_Grd_x0_8,double F_Grd_xl_8,double F_Grd_y0_8,double F_Grd_yl_8, int F_Grd_yinyang_L){

   int i;
   double delta_8;

   delta_8 = (F_Grd_xl_8-F_Grd_x0_8)/(F_Grd_ni-1);
   F_xgi_8[0] = F_Grd_x0_8;
   F_xgi_8[F_Grd_ni-1] = F_Grd_xl_8;
   for(i=1;i<F_Grd_ni-1;i++) F_xgi_8[i]= F_Grd_x0_8 + i*delta_8;

   delta_8 = (F_Grd_yl_8-F_Grd_y0_8)/(F_Grd_nj-1);
   F_ygi_8[0] = F_Grd_y0_8;
   F_ygi_8[F_Grd_nj-1] = F_Grd_yl_8;
   for(i=1;i<F_Grd_nj-1;i++) F_ygi_8[i]= F_Grd_y0_8 + i*delta_8;

   if (F_Grd_yinyang_L) {
      *F_Grd_dx   = fabs(F_xgi_8[1]-F_xgi_8[0]);
      *F_Grd_dy   = fabs(F_ygi_8[1]-F_ygi_8[0]);
   }
}

/*----------------------------------------------------------------------------
 * @brief  Define a referential of type RPN ZE
 * @author Jean-Philippe Gauthier
 * @date   Avril 2015
 *    @param[in] Ref     Georef definition
 *    @param[in] NI      Dimension en X
 *    @param[in] NJ      Dimension en Y
 *    @param[in] DX      Resolution en X
 *    @param[in] DY      Resolution en Y
 *    @param[in] LatR     
 *    @param[in] LonR     
 *    @param[in] MaxCFL   
 *    @param[in] XLat1   Latitude centrale
 *    @param[in] XLon1   Longitude centrale
 *    @param[in] XLat2   Latitude de l'axe de rotation
 *    @param[in] XLon2   Longitude de l'axe de rotation
 
 *    @return            Georef (NULL=error)
*/
TGeoRef* GeoRef_DefineZE(TGeoRef *Ref,int NI,int NJ,float DX,float DY,float LatR,float LonR,int MaxCFL,float XLat1,float XLon1,float XLat2,float XLon2) {

#ifdef HAVE_RMN
   int    bsc_base,bsc_ext1,extension;
   double x0,x1,y0,y1;

   if (!Ref) {
      return(NULL);
   }
 
   f77name(cxgaig)("E",&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4],&XLat1,&XLon1,&XLat2,&XLon2);
   f77name(cigaxg)("E",&XLat1,&XLon1,&XLat2,&XLon2,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);
   
   GEM_grid_param(&bsc_base,&bsc_ext1,&extension,MaxCFL,&LonR,&LatR,&NI,&NJ,&DX,&DY,&x0,&y0,&x1,&y1,-1,FALSE);
 
   if (NI!=Ref->NX+1 || NJ!=Ref->NY+1) {
      Ref->AX=realloc(Ref->AX,NI*sizeof(double));
      Ref->AY=realloc(Ref->AY,NJ*sizeof(double));
   }

   // f77name(set_gemhgrid4)(Ref->AX,Ref->AY,&NI,&NJ,&DX,&DY,&x0,&x1,&y0,&y1,FALSE);
   GEM_hgrid4(Ref->AX,Ref->AY,NI,NJ,&DX,&DY,x0,x1,y0,y1,FALSE);
           
   GeoRef_Define(Ref,NI,NJ,"Z","E",Ref->RPNHead.IG[X_IG1],Ref->RPNHead.IG[X_IG2],Ref->RPNHead.IG[X_IG3],Ref->RPNHead.IG[X_IG4],Ref->AX,Ref->AY);
 #endif
  
   return(Ref);
}
