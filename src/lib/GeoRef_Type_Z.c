#include <App.h>
#include "GeoRef.h"

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for a Z grid
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 *    @param[in]  X       X array
 *    @param[in]  Y       Y array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int32_t GeoRef_XY2LL_Z(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t Nb) {

   int32_t     i,indx,indy,un=1;
   double  delxx,delyy,*tmpx,*tmpy,*ytmp;
   float   xlat1, xlon1, xlat2, xlon2;

   if (!Ref->AX || !Ref->AY) {
      return(FALSE);
   }
   
   tmpx = (double*)malloc(3*Nb*sizeof(double));
   tmpy = &tmpx[Nb];
   ytmp = &tmpx[Nb*2];

   #pragma omp parallel for default(none) private(i,indx,indy,delxx,delyy) shared(Nb,Ref,X,Y,Lat,Lon,tmpx,tmpy,ytmp)
   for (i=0; i < Nb; i++) {
      indx = (int)X[i]-1;
      ytmp[i] = (Ref->RPNHead.ig2==1)?Ref->NY+1.0 - Y[i]:Y[i];
      
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

   switch (Ref->RPNHeadExt.grref[0]) {
      case 'E':
         f77name(cigaxg)(Ref->RPNHeadExt.grref,&xlat1,&xlon1,&xlat2,&xlon2,&Ref->RPNHeadExt.igref1,&Ref->RPNHeadExt.igref2,&Ref->RPNHeadExt.igref3,&Ref->RPNHeadExt.igref4,1);
         GeoRef_RotateInvertXY(Lat,Lon,tmpx,tmpy,Nb,Ref->RPNHeadExt.xgref1,Ref->RPNHeadExt.xgref2,Ref->RPNHeadExt.xgref3,Ref->RPNHeadExt.xgref4);
         break;

      case 'S':
      case 'N':
         f77name(ez8_vllfxy)(Lat,Lon,tmpx,tmpy,&Nb,&un,&Ref->RPNHeadExt.xgref3,&Ref->RPNHeadExt.xgref4,&Ref->RPNHeadExt.xgref1,&Ref->RPNHeadExt.xgref2,&Ref->Hemi);
         break;

      case 'L':
         for (i=0; i < Nb; i++) {
            Lat[i] = (tmpy[i])*Ref->RPNHeadExt.xgref3+Ref->RPNHeadExt.xgref1;
            Lon[i] = (tmpx[i])*Ref->RPNHeadExt.xgref4+Ref->RPNHeadExt.xgref2;
         }
         break;

      case 'W':
         GeoRef_XY2LL_W(Ref,Lat,Lon,tmpx,tmpy,Nb);
         break;

      default:
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Undefined reference grid type: %c\n",__func__,Ref->RPNHeadExt.grref[0]);
         free(tmpx);
         return(1);
         break;
   }
   free(tmpx);

   return(Nb);
}

/*----------------------------------------------------------------------------
 * @brief  Transforms LatLon coordinates to XY for a Z grid
 * @date   June 2015
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] X       X array
 *    @param[out] Y       Y array
 *    @param[in]  Lat     Latitude array
 *    @param[in]  Lon     Longitude array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int32_t GeoRef_LL2XY_Z(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t Nb) {

   int32_t   i,j,indx,indy,d;
    
   if (!Ref->AX || !Ref->AY) {
      return(FALSE);
   }
   
   GeoRef_LL2GREF(Ref,X,Y,Lat,Lon,Nb);

   d=(Ref->Type&GRID_AXY2D)?Ref->NX:1;

   // Look into expansion descriptor
   #pragma omp parallel for default(none) private(i,indx,indy) shared(stderr,d,Nb,Ref,X,Y,Lat,Lon)
   for(i=0;i<Nb;i++) {
      //TODO: clarify NX and j2 index          
      indx = GeoRef_XFind(X[i],Ref->AX,Ref->NX,1);
      indy = GeoRef_XFind(Y[i],Ref->AY,Ref->NY,d);
    
      if (indx >= Ref->NX-1) indx = Ref->NX - 2;
      if (indy >= Ref->NY-1) indy = Ref->NY - 2;

      X[i] = indx+(X[i]-Ref->AX[indx])/(Ref->AX[indx+1]-Ref->AX[indx])+1;
      Y[i] = indy+(Y[i]-Ref->AY[indy*d])/(Ref->AY[(indy+1)*d]-Ref->AY[indy*d])+1;
   } 

   if (Ref->GRTYP[0] == 'G') {
      if (Ref->RPNHead.ig1 == 1) {
         for (j=0; j < Nb; j++) Y[j] = Y[j] - Ref->j2;
      }
      if (Ref->RPNHead.ig2 == 1)  {
         for (j=0; j < Nb; j++) Y[j] = Ref->j2 +1.0 - Y[j];
      }
      // TODO: From GeoRef_RPN Fix for G grid 0-360 1/5 gridpoint problem
      for (j=0; j < Nb; j++) if (X[j]>Ref->X1+0.5) X[j]-=(Ref->X1+1);
   }
      
   return(Nb);
}

int32_t GEM_grid_param(int32_t *F_bsc_base,int32_t *F_bsc_ext1,int32_t *F_extension ,int32_t F_maxcfl,float *F_lonr,float *F_latr,int32_t *F_ni,int32_t *F_nj,float *F_dx,float *F_dy,double *F_x0_8,double *F_y0_8,double *F_xl_8,double *F_yl_8,int32_t F_overlap,int32_t F_yinyang_L) {

   double delta_8;
   int32_t iref,jref;
  
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
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Longitude of WEST %f < 0.0\n",__func__,*F_x0_8);
         return(0);
      }
      if (*F_y0_8 < -90.) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Latitude of SOUTH %f < 0.0\n",__func__,*F_y0_8);
         return(0);
      }
      if (*F_xl_8 > 360.) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Longitude of EAST %f < 0.0\n",__func__,*F_xl_8);
         return(0);
      }
      if (*F_yl_8 > 90.) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Latitude of NORTH %f < 0.0\n",__func__,*F_yl_8);
         return(0);
      }
   }
   return(1);
}

void GEM_hgrid4(double *F_xgi_8,double *F_ygi_8,int32_t F_Grd_ni,int32_t F_Grd_nj,float *F_Grd_dx,float *F_Grd_dy,double F_Grd_x0_8,double F_Grd_xl_8,double F_Grd_y0_8,double F_Grd_yl_8, int32_t F_Grd_yinyang_L){

   int32_t i;
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
TGeoRef* GeoRef_DefineZE(TGeoRef *Ref,int32_t NI,int32_t NJ,float DX,float DY,float LatR,float LonR,int32_t MaxCFL,float XLat1,float XLon1,float XLat2,float XLon2) {

   int32_t    bsc_base,bsc_ext1,extension;
   double x0,x1,y0,y1;

   if (!Ref) {
      return(NULL);
   }
 
   f77name(cxgaig)("E",&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4,&XLat1,&XLon1,&XLat2,&XLon2,1);
   f77name(cigaxg)("E",&XLat1,&XLon1,&XLat2,&XLon2,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4,1);
   
   GEM_grid_param(&bsc_base,&bsc_ext1,&extension,MaxCFL,&LonR,&LatR,&NI,&NJ,&DX,&DY,&x0,&y0,&x1,&y1,-1,FALSE);
 
   if (NI!=Ref->NX+1 || NJ!=Ref->NY+1) {
      Ref->AX=realloc(Ref->AX,NI*sizeof(double));
      Ref->AY=realloc(Ref->AY,NJ*sizeof(double));
   }

   // f77name(set_gemhgrid4)(Ref->AX,Ref->AY,&NI,&NJ,&DX,&DY,&x0,&x1,&y0,&y1,FALSE);
   GEM_hgrid4(Ref->AX,Ref->AY,NI,NJ,&DX,&DY,x0,x1,y0,y1,FALSE);
           
   Ref=GeoRef_Define(Ref,NI,NJ,"Z","E",Ref->RPNHead.ig1,Ref->RPNHead.ig2,Ref->RPNHead.ig3,Ref->RPNHead.ig4,Ref->AX,Ref->AY);
  
   return(Ref);
}
