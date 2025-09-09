#include <App.h>
#include "GeoRef.h"
#include "GeoRef_f.h"

void c_ezgfwfllw8(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int32_t *ni,int32_t *nj,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4);
void c_ezllwfgfw8(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int32_t *ni,int32_t *nj,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4);
void c_ezllwfgff8(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int32_t *ni,int32_t *nj,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4);

/**----------------------------------------------------------------------------
 * @brief  Interpolates vectorial component values between 2 georeferences
 *    @param[in]  RefTo      Destination geo-reference
 *    @param[in]  RefFrom    Source geo-reference
 *    @param[in]  Opt        Interpolation options
 *    @param[out] uuout      Destination UU interpolated values
 *    @param[out] vvout      Destination VV interpolated values
 *    @param[in]  uuin       Source UU values
 *    @param[in]  vvin       Source VV values

 *    @return                FALSE (0) if operation failed, TRUE (1) otherwise
*/
int32_t GeoRef_InterpUV(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin) {

   int32_t npts;
   float *uullout = NULL;
   float *vvllout = NULL;
   TGeoSet   *set;
   TGeoOptions opt;

   if (!Opt) Opt=&GeoRef_Options;

   // If we're using pre-calculated index weight
   if (Opt->Interp==IR_WEIGHTINDEX) {
       GeoRef_InterpClear(RefTo, Opt, uuout);
       GeoRef_InterpClear(RefTo, Opt, vvout);
       return(GeoRef_InterpWeight(RefTo,RefFrom,Opt,uuout,vvout,uuin,vvin));
   }

   set=GeoRef_SetGet(RefTo,RefFrom,Opt);

   if (RefFrom->NbSub > 0 || RefTo->NbSub > 0) {
      return(GeoRef_InterpYYUV(RefTo,RefFrom,Opt,uuout,vvout,uuin,vvin));
   } else {
      npts = RefTo->NX*RefTo->NY;
      memcpy(&opt,Opt,sizeof(TGeoOptions));

      opt.VectorMode = TRUE;
      opt.Symmetric = TRUE;
      if (!GeoRef_Interp(RefTo,RefFrom,&opt,uuout,uuin)) {
         return(FALSE);
      }

      opt.Symmetric = FALSE;
      if (!GeoRef_Interp(RefTo,RefFrom,&opt,vvout,vvin)) {
         return(FALSE);
      }

      if (opt.PolarCorrect == TRUE) {
         GeoRef_CorrectVector(set,uuout,vvout,uuin,vvin);
      }

      uullout = (float *)malloc(2*npts*sizeof(float));
      vvllout = &uullout[npts];

      GeoRef_UV2WD(RefFrom,uullout,vvllout,uuout,vvout,RefTo->Lat,RefTo->Lon,npts);
      GeoRef_WD2UV(RefTo,uuout,vvout,uullout,vvllout,RefTo->Lat,RefTo->Lon,npts);

      free(uullout);
   }

   return(TRUE);
}

/**----------------------------------------------------------------------------
 * @brief  Interpolates vectorial component values between 2 georeferences and outputs rotated speed and directions
 *    @param[in]  RefTo      Destination geo-reference
 *    @param[in]  RefFrom    Source geo-reference
 *    @param[in]  Opt        Interpolation options
 *    @param[out] uuout      Destination speed interpolated values
 *    @param[out] vvout      Destination direction interpolated values
 *    @param[in]  uuin       Source UU values
 *    @param[in]  vvin       Source VV values

 *    @return                FALSE (0) if operation failed, TRUE (1) otherwise
*/
int32_t GeoRef_InterpWD(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin) {

   int32_t npts;
   float *uullout = NULL;
   float *vvllout = NULL;
   TGeoSet   *set;
   TGeoOptions opt;

   if (!Opt) Opt=&GeoRef_Options;
   set=GeoRef_SetGet(RefTo,RefFrom,Opt);

   if (RefFrom->NbSub > 0 || RefTo->NbSub > 0) {
      return(GeoRef_InterpYYWD(RefTo,RefFrom,Opt,uuout,vvout,uuin,vvin));
   } else {

      npts = RefTo->NX*RefTo->NY;

      memcpy(&opt,Opt,sizeof(TGeoOptions));

      opt.VectorMode = TRUE;
      opt.Symmetric = TRUE;
      if (!GeoRef_Interp(RefTo,RefFrom,&opt,uuout,uuin)) {
         return(FALSE);
      }

      opt.Symmetric = FALSE;
      if (!GeoRef_Interp(RefTo,RefFrom,&opt,vvout,vvin)) {
         return(FALSE);
      }

      if (opt.PolarCorrect == TRUE) {
         GeoRef_CorrectVector(set,uuout,vvout,uuin,vvin);
      }

      uullout = (float *)malloc(2*npts*sizeof(float));
      vvllout = &uullout[npts];

      // ezsint32_t does not allocate lat,lon if RefFrom=RefTo
      GeoRef_UV2WD(RefFrom,uullout,vvllout,uuout,vvout,RefTo->Lat,RefTo->Lon,npts);

      memcpy(uuout, uullout, npts*sizeof(float));
      memcpy(vvout, vvllout, npts*sizeof(float));

      free(uullout);
   }
   return(TRUE);
}

int32_t GeoRef_InterpYYUV(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin) {

   TGeoSet *gset=NULL;
   TGeoRef *yin_gdin, *yan_gdin, *yin_gdout, *yan_gdout;
   int32_t i,j,k;
   int32_t yancount_yin,yincount_yin, yancount_yan,yincount_yan;
   int32_t yyin,yyout;
   int32_t ni, nj,idx;
   float *yin2yin_uuout,*yan2yin_uuout, *yin2yin_vvout,*yan2yin_vvout;
   float *yin2yan_uuout,*yan2yan_uuout, *yin2yan_vvout,*yan2yan_vvout;
   float *yin2yin_spdout,*yan2yin_spdout, *yin2yin_wdout,*yan2yin_wdout;
   float *yin2yan_spdout,*yan2yan_spdout, *yin2yan_wdout,*yan2yan_wdout;
   float *spdout,*wdout;

   // Need only access to either yin or Yang info for the lat and lon val

   yyin=0; yyout=0;

   if (!Opt) Opt=&GeoRef_Options;
   gset=GeoRef_SetGet(RefTo,RefFrom,Opt);

   // Setup for input grid
   if (RefFrom->NbSub > 0) {
      yyin=1;
      yin_gdin = RefFrom->Subs[0];
      yan_gdin = RefFrom->Subs[1];
   } else {
      yin_gdin = RefFrom;
   }

   // Setup for output grid
   if (RefTo->NbSub > 0) {
      yyout=1;
      yin_gdout = RefTo->Subs[0];
      yan_gdout = RefTo->Subs[1];
   } else {
      yin_gdout = RefTo;
   }

   ni = yin_gdout->NX;
   nj = yin_gdout->NY;

   // Interp input one grid to yygrid - no masking needed
   if (yyin == 0 && yyout == 1) {
      if (GeoRef_InterpUV(yin_gdout,RefFrom,Opt,uuout,vvout,uuin,vvin) && GeoRef_InterpUV(yan_gdout,RefFrom,Opt,&uuout[ni*nj],&vvout[ni*nj],uuin,vvin)) {
         return(TRUE);
      } else {
         return(FALSE);
      }
   }

   // Check if one sub grid is identical to one of the sub grids
   if (yin_gdin == RefTo) {
      return(GeoRef_InterpUV(RefTo,yin_gdin,Opt,uuout,vvout,uuin,vvin));
   }
   if (yan_gdin == RefTo) {
      return(GeoRef_InterpUV(RefTo,yan_gdin,Opt,uuout,vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)]));
   }

   /* User specifies to use 1 subgrid for interpolation ezsetopt(USE_1SUBGRID) */
   /* User must specify one specific grid ezsetival(SUBGRIDID) */
   /* This is only appropriate if the destination grid is non yin-yang grid */

   if (RefFrom->Sub>=0) { // User specifies to use 1 grid only
      // Output is a Yin-Yang grid
      if (yyout == 1) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Cannot use subgrid to interpolate to a Yin-Yang grid\n",__func__);
         return(FALSE);
      }
      // Is specified subgrid within the subgrid list
      if (RefFrom->Sub>=RefFrom->NbSub) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid subgrid: %i\n",__func__,RefFrom->Sub);
         return(FALSE);
      }
      // Use yin input grid
      if (RefFrom->Sub==0) {
         return(GeoRef_InterpUV(yin_gdout,yin_gdin,Opt,uuout,vvout,uuin,vvin));
      }

      // Use yang input grid
      if (RefFrom->Sub==1) {
         return(GeoRef_InterpUV(yin_gdout,yan_gdin,Opt,uuout,vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)]));
      }
   }

   /* To use both Yin and Yang grids in Yin-yang input grid */
   /* Masquer les grilles YY input pour enlever overlap et calculer les X,Y */
   GeoRef_SetCalcYYXY(gset);

   // Interp yinyang to one grid
   if (yyin == 1 && yyout == 0) {
      yincount_yin = gset->yincount_yin;
      yancount_yin = gset->yancount_yin;
      idx=0;
      spdout        = (float *) malloc(((2*ni*nj)+(4*yincount_yin)+(4*yancount_yin))*sizeof(float));
      wdout         = &spdout[idx+=ni*nj];
      yin2yin_uuout = &spdout[idx+=ni*nj];
      yin2yin_vvout = &spdout[idx+=yincount_yin];
      yin2yin_spdout = &spdout[idx+=yincount_yin];
      yin2yin_wdout = &spdout[idx+=yincount_yin];
      yan2yin_uuout = &spdout[idx+=yincount_yin];
      yan2yin_vvout = &spdout[idx+=yancount_yin];
      yan2yin_spdout = &spdout[idx+=yancount_yin];
      yan2yin_wdout = &spdout[idx+=yancount_yin];

      GeoRef_XYUVVal(yin_gdin,Opt,yin2yin_uuout,yin2yin_vvout,uuin,vvin,gset->yin2yin_x,gset->yin2yin_y,gset->yincount_yin);
      GeoRef_UV2WD(yin_gdin,yin2yin_spdout,yin2yin_wdout,yin2yin_uuout,yin2yin_vvout,gset->yin2yin_lat,gset->yin2yin_lon,yincount_yin);
      GeoRef_XYUVVal(yan_gdin,Opt,yan2yin_uuout,yan2yin_vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yin_x,gset->yan2yin_y,yancount_yin);
      GeoRef_UV2WD(yan_gdin,yan2yin_spdout,yan2yin_wdout,yan2yin_uuout, yan2yin_vvout,gset->yan2yin_lat,gset->yan2yin_lon,yancount_yin);

      yincount_yin=0;
      yancount_yin=0;
      for(j=0; j<nj; j++) {
         for (i=0;i<ni; i++) {
            k=(j*ni)+i;
            if (gset->yin_maskout[k] == 1.0) {
               spdout[k]=yan2yin_spdout[yancount_yin];
               wdout[k]=yan2yin_wdout[yancount_yin];
               yancount_yin++;
            } else {
               spdout[k]=yin2yin_spdout[yincount_yin];
               wdout[k]=yin2yin_wdout[yincount_yin];
               yincount_yin++;
            }
         }
      }
      GeoRef_WD2UV(RefTo,uuout,vvout,spdout,wdout,gset->yinlat,gset->yinlon,ni*nj);
      free(spdout);
   }

   // Interp yinyang to yinyang
   if (yyout == 1 && yyin == 1) {
      // Interp input YY grid to YIN

      yincount_yin = gset->yincount_yin;
      yancount_yin = gset->yancount_yin;
      yincount_yan = gset->yincount_yan;
      yancount_yan = gset->yancount_yan;
      idx=0;
      spdout        = (float *) malloc(((2*ni*nj)+(4*yincount_yin)+(4*yancount_yin))*sizeof(float));
      wdout         = &spdout[idx+=ni*nj];
      yin2yin_uuout = &spdout[idx+=ni*nj];
      yin2yin_vvout = &spdout[idx+=yincount_yin];
      yin2yin_spdout = &spdout[idx+=yincount_yin];
      yin2yin_wdout = &spdout[idx+=yincount_yin];
      yan2yin_uuout = &spdout[idx+=yincount_yin];
      yan2yin_vvout = &spdout[idx+=yancount_yin];
      yan2yin_spdout = &spdout[idx+=yancount_yin];
      yan2yin_wdout = &spdout[idx+=yancount_yin];

      spdout = (float *) malloc(((2*ni*nj)+(4*yincount_yin)+(4*yancount_yin)+(4*yincount_yan)+(4*yancount_yan))*sizeof(float));
      wdout          = &spdout[idx+=ni*nj];
      yin2yin_uuout  = &spdout[idx+=ni*nj];
      yin2yin_vvout  = &spdout[idx+=yincount_yin];
      yin2yin_spdout = &spdout[idx+=yincount_yin];
      yin2yin_wdout  = &spdout[idx+=yincount_yin];
      yan2yin_uuout  = &spdout[idx+=yincount_yin];
      yan2yin_vvout  = &spdout[idx+=yancount_yin];
      yan2yin_spdout = &spdout[idx+=yancount_yin];
      yan2yin_wdout  = &spdout[idx+=yancount_yin];

      yin2yan_uuout  = &spdout[idx+=yancount_yin];
      yin2yan_vvout  = &spdout[idx+=yincount_yan];
      yin2yan_spdout = &spdout[idx+=yincount_yan];
      yin2yan_wdout  = &spdout[idx+=yincount_yan];
      yan2yan_uuout  = &spdout[idx+=yincount_yan];
      yan2yan_vvout  = &spdout[idx+=yancount_yan];
      yan2yan_spdout = &spdout[idx+=yancount_yan];
      yan2yan_wdout  = &spdout[idx+=yancount_yan];

      GeoRef_XYUVVal(yin_gdin,Opt,yin2yin_uuout,yin2yin_vvout,uuin,vvin,gset->yin2yin_x,gset->yin2yin_y,yincount_yin);
      GeoRef_UV2WD(yin_gdin,yin2yin_spdout,yin2yin_wdout,yin2yin_uuout,yin2yin_vvout,gset->yin2yin_lat,gset->yin2yin_lon,yincount_yin);
      GeoRef_XYUVVal(yan_gdin,Opt,yan2yin_uuout,yan2yin_vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yin_x,gset->yan2yin_y,yancount_yin);
      GeoRef_UV2WD(yan_gdin,yan2yin_spdout,yan2yin_wdout,yan2yin_uuout,yan2yin_vvout,gset->yan2yin_lat,gset->yan2yin_lon,yancount_yin);

      GeoRef_XYUVVal(yin_gdin,Opt,yin2yan_uuout,yin2yan_vvout,uuin,vvin,gset->yin2yan_x,gset->yin2yan_y,yincount_yan);
      GeoRef_UV2WD(yin_gdin,yin2yan_spdout,yin2yan_wdout,yin2yan_uuout,yin2yan_vvout,gset->yin2yan_lat,gset->yin2yan_lon,yincount_yan);

      GeoRef_XYUVVal(yan_gdin,Opt,yan2yan_uuout,yan2yan_vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yan_x,gset->yan2yan_y,yancount_yan);
      GeoRef_UV2WD(yan_gdin,yan2yan_spdout,yan2yan_wdout,yan2yan_uuout,yan2yan_vvout,gset->yan2yan_lat,gset->yan2yan_lon,yancount_yan);

      // Build output for YIN output grid
      yincount_yin=0; yancount_yin=0;
      for(j=0; j<nj; j++) {
         for (i=0;i<ni; i++) {
            k=(j*ni)+i;
            if (gset->yin_maskout[k] == 1.0) {
               spdout[k]=yan2yin_spdout[yancount_yin];
               wdout[k]=yan2yin_wdout[yancount_yin];
               yancount_yin++;
            } else {
               spdout[k]=yin2yin_spdout[yincount_yin];
               wdout[k]=yin2yin_wdout[yincount_yin];
               yincount_yin++;
            }
         }
      }
      GeoRef_WD2UV(yin_gdout,uuout,vvout,spdout,wdout,gset->yinlat,gset->yinlon,ni*nj);

      // Build output for YIN output grid
      yincount_yan=0; yancount_yan=0;
      for(j=0; j<nj; j++) {
        for (i=0;i<ni; i++) {
           k=(j*ni)+i;
           if (gset->yan_maskout[k] == 1.0) {
              spdout[k]=yan2yan_spdout[yancount_yan];
              wdout[k]=yan2yan_wdout[yancount_yan];
              yancount_yan++;
           } else {
              spdout[k]=yin2yan_spdout[yincount_yan];
              wdout[k]=yin2yan_wdout[yincount_yan];
              yincount_yan++;
           }
         }
      }

      GeoRef_WD2UV(yan_gdout,&uuout[ni*nj],&vvout[ni*nj],spdout,wdout,gset->yanlat,gset->yanlon,ni*nj);
      free(spdout);
   }

   return(TRUE);
}

int32_t GeoRef_InterpYYWD(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin) {

   TGeoSet *gset=NULL;
   TGeoRef *yin_gdin, *yan_gdin, *yin_gdout, *yan_gdout;
   int32_t i,j,k;
   int32_t yancount_yin,yincount_yin, yancount_yan,yincount_yan;
   int32_t yyin,yyout;
   int32_t ni, nj,idx;
   float *yin2yin_uuout,*yan2yin_uuout, *yin2yin_vvout,*yan2yin_vvout;
   float *yin2yan_uuout,*yan2yan_uuout, *yin2yan_vvout,*yan2yan_vvout;
   float *yin2yin_spdout,*yan2yin_spdout, *yin2yin_wdout,*yan2yin_wdout;
   float *yin2yan_spdout,*yan2yan_spdout, *yin2yan_wdout,*yan2yan_wdout;

   // Need only access to either yin or Yang info for the lat and lon val

   yyin=0; yyout=0;

   if (!Opt) Opt=&GeoRef_Options;
   gset=GeoRef_SetGet(RefTo,RefFrom,Opt);

   // Setup for input grid
   if (RefFrom->NbSub > 0) {
      yyin=1;
      yin_gdin = RefFrom->Subs[0];
      yan_gdin = RefFrom->Subs[1];
   } else {
      yin_gdin = RefFrom;
   }

  // Setup for input grid
   if (RefTo->NbSub > 0) {
      yyout=1;
      yin_gdout = RefTo->Subs[0];
      yan_gdout = RefTo->Subs[1];
   } else {
      yin_gdout = RefTo;
   }

   ni = yin_gdout->NX;
   nj = yin_gdout->NY;

   // Interp input one grid to yygrid - no masking needed
   if (yyin == 0 && yyout == 1) {
      if (GeoRef_InterpWD(yin_gdout,RefFrom,Opt,uuout,vvout,uuin,vvin) && GeoRef_InterpWD(yan_gdout,RefFrom,Opt,&uuout[(ni*nj)],&vvout[(ni*nj)],uuin,vvin)) {
         return(TRUE);
      } else {
         return(FALSE);
      }
    }

   // Check if one sub grid is identical to one of the sub grids
   if (yin_gdin == RefTo) {
      return(GeoRef_InterpWD(RefTo,yin_gdin,Opt,uuout,vvout,uuin,vvin));
   }
   if (yan_gdin == RefTo) {
      return(GeoRef_InterpWD(RefTo,yan_gdin,Opt,uuout,vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)]));
   }

   /* User specifies to use 1 subgrid for interpolation ezsetopt(USE_1SUBGRID) */
   /* User must specify one specific grid ezsetival(SUBGRIDID) */
   /* This is only appropriate if the destination grid is non yin-yang grid */

   if (RefFrom->Sub>=0) { // User specifies to use 1 grid only
      // Output is a Yin-Yang grid
      if (yyout == 1) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Cannot use subgrid to interpolate to a Yin-Yang grid\n",__func__);
         return(FALSE);
      }
      // Is specified subgrid within the subgrid list
      if (RefFrom->Sub>=RefFrom->NbSub) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid subgrid: %i\n",__func__,RefFrom->Sub);
         return(FALSE);
      }
      // Use yin input grid
      if (RefFrom->Sub==0) {
         return(GeoRef_InterpWD(yin_gdout,yin_gdin,Opt,uuout,vvout,uuin,vvin));
      }

      // Use yang input grid
      if (RefFrom->Sub==1) {
         return(GeoRef_InterpWD(yin_gdout,yan_gdin,Opt,uuout,vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)]));
      }
   }

   /* To use both Yin and Yang grids in Yin-yang input grid */
   /* Masquer les grilles YY input pour enlever overlap et calculer les X,Y */
   GeoRef_SetCalcYYXY(gset);

   // Interp yinyang to one grid
   if (yyin == 1 && yyout == 0) {
      yincount_yin = gset->yincount_yin;
      yancount_yin = gset->yancount_yin;
      idx=0;
      yin2yin_uuout  = (float*)malloc((4*yincount_yin+4*yancount_yin)*sizeof(float));
      yin2yin_vvout  = &yin2yin_uuout[idx+=yincount_yin];
      yin2yin_spdout = &yin2yin_uuout[idx+=yincount_yin];
      yin2yin_wdout  = &yin2yin_uuout[idx+=yincount_yin];
      yan2yin_uuout  = &yin2yin_uuout[idx+=yincount_yin];
      yan2yin_vvout  = &yin2yin_uuout[idx+=yancount_yin];
      yan2yin_spdout = &yin2yin_uuout[idx+=yancount_yin];
      yan2yin_wdout  = &yin2yin_uuout[idx+=yancount_yin];

      GeoRef_XYUVVal(yin_gdin,Opt,yin2yin_uuout,yin2yin_vvout,uuin,vvin,gset->yin2yin_x,gset->yin2yin_y,gset->yincount_yin);
      GeoRef_UV2WD(yin_gdin,yin2yin_spdout,yin2yin_wdout,yin2yin_uuout,yin2yin_vvout,gset->yin2yin_lat,gset->yin2yin_lon,yincount_yin);

      GeoRef_XYUVVal(yan_gdin,Opt,yan2yin_uuout,yan2yin_vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yin_x,gset->yan2yin_y,yancount_yin);
      GeoRef_UV2WD(yan_gdin,yan2yin_spdout,yan2yin_wdout,yan2yin_uuout, yan2yin_vvout,gset->yan2yin_lat,gset->yan2yin_lon,yancount_yin);
      yincount_yin=0;
      yancount_yin=0;
      for(j=0; j<nj; j++) {
         for (i=0;i<ni; i++) {
            k=(j*ni)+i;
            if (gset->yin_maskout[k] == 1.0) {
               uuout[k]=yan2yin_spdout[yancount_yin];
               vvout[k]=yan2yin_wdout[yancount_yin];
               yancount_yin++;
            } else {
               uuout[k]=yin2yin_spdout[yincount_yin];
               vvout[k]=yin2yin_wdout[yincount_yin];
               yincount_yin++;
            }
         }
      }
      free(yin2yin_uuout);
   }

   // Interp yinyang to yinyang
   if (yyout == 1 && yyin == 1) {
      // Interp input YY grid to YIN
      yincount_yin = gset->yincount_yin;
      yancount_yin = gset->yancount_yin;
      yincount_yan = gset->yincount_yan;
      yancount_yan = gset->yancount_yan;
      idx=0;
      yin2yin_uuout  = (float*)malloc((4*yincount_yin+4*yancount_yin+4*yincount_yan+4*yancount_yan)*sizeof(float));
      yin2yin_vvout  = &yin2yin_uuout[idx+=yincount_yin];
      yin2yin_spdout = &yin2yin_uuout[idx+=yincount_yin];
      yin2yin_wdout  = &yin2yin_uuout[idx+=yincount_yin];
      yan2yin_uuout  = &yin2yin_uuout[idx+=yincount_yin];
      yan2yin_vvout  = &yin2yin_uuout[idx+=yancount_yin];
      yan2yin_spdout = &yin2yin_uuout[idx+=yancount_yin];
      yan2yin_wdout  = &yin2yin_uuout[idx+=yancount_yin];

      yin2yan_uuout  = &yin2yin_uuout[idx+=yancount_yin];
      yin2yan_vvout  = &yin2yin_uuout[idx+=yincount_yan];
      yin2yan_spdout = &yin2yin_uuout[idx+=yincount_yan];
      yin2yan_wdout  = &yin2yin_uuout[idx+=yincount_yan];
      yan2yan_uuout  = &yin2yin_uuout[idx+=yincount_yan];
      yan2yan_vvout  = &yin2yin_uuout[idx+=yancount_yan];
      yan2yan_spdout = &yin2yin_uuout[idx+=yancount_yan];
      yan2yan_wdout  = &yin2yin_uuout[idx+=yancount_yan];

      GeoRef_XYUVVal(yin_gdin,Opt,yin2yin_uuout,yin2yin_vvout,uuin,vvin,gset->yin2yin_x,gset->yin2yin_y,gset->yincount_yin);
      GeoRef_UV2WD(yin_gdin,yin2yin_spdout,yin2yin_wdout,yin2yin_uuout,yin2yin_vvout,gset->yin2yin_lat,gset->yin2yin_lon,gset->yincount_yin);
      GeoRef_XYUVVal(yan_gdin,Opt,yan2yin_uuout,yan2yin_vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yin_x,gset->yan2yin_y,gset->yancount_yin);
      GeoRef_UV2WD(yan_gdin,yan2yin_spdout,yan2yin_wdout,yan2yin_uuout,yan2yin_vvout,gset->yan2yin_lat,gset->yan2yin_lon,yancount_yin);

      GeoRef_XYUVVal(yin_gdin,Opt,yin2yan_uuout,yin2yan_vvout,uuin,vvin,gset->yin2yan_x,gset->yin2yan_y,gset->yincount_yan);
      GeoRef_UV2WD(yin_gdin,yin2yan_spdout,yin2yan_wdout,yin2yan_uuout,yin2yan_vvout,gset->yin2yan_lat,gset->yin2yan_lon,gset->yincount_yan);

      GeoRef_XYUVVal(yan_gdin,Opt,yan2yan_uuout,yan2yan_vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yan_x,gset->yan2yan_y,gset->yancount_yan);
      GeoRef_UV2WD(yan_gdin,yan2yan_spdout,yan2yan_wdout,yan2yan_uuout,yan2yan_vvout,gset->yan2yan_lat,gset->yan2yan_lon,gset->yancount_yan);

      // Build output for YIN output grid
      yincount_yin=0; yancount_yin=0;
      for(j=0; j<nj; j++) {
         for (i=0;i<ni; i++) {
            k=(j*ni)+i;
            if (gset->yin_maskout[k] == 1.0) {
               uuout[k]=yan2yin_spdout[yancount_yin];
               vvout[k]=yan2yin_wdout[yancount_yin];
               yancount_yin++;
            } else {
               uuout[k]=yin2yin_spdout[yincount_yin];
               vvout[k]=yin2yin_wdout[yincount_yin];
               yincount_yin++;
            }
         }
      }

      // Build output for YANG output grid
      yincount_yan=0; yancount_yan=0;
      for(j=0; j<nj; j++) {
         for (i=0;i<ni; i++) {
            k=(j*ni)+i;
            if (gset->yan_maskout[k] == 1.0) {
               uuout[k+(ni*nj)]=yan2yan_spdout[yancount_yan];
               vvout[k+(ni*nj)]=yan2yan_wdout[yancount_yan];
               yancount_yan++;
            } else {
               uuout[k+(ni*nj)]=yin2yan_spdout[yincount_yan];
               vvout[k+(ni*nj)]=yin2yan_wdout[yincount_yan];
               yincount_yan++;
            }
         }
      }
      free(yin2yin_uuout);
   }

   return(TRUE);
}

/**----------------------------------------------------------------------------
 * @brief  Converts meteorological winds (speed/direction) to grid winds (uu/vv).

 *    @param[in]  Ref          geo-reference
 *    @param[out] uuout        Destination UU values
 *    @param[out] vvout        Destination VV values
 *    @param[in]  spdin        Source speed values
 *    @param[in]  wdin         Source direction values
 *    @param[in]  Lat          Latitude of wind value points (NULL, use geo-reference gridpoints)
 *    @param[in]  Lon          Longitude of wind value points (NULL, use geo-reference gridpoints)
 *    @param[in]  Nb           Number of values to convert (size of arrays)

 *    @return                FALSE (0) if operation failed, TRUE (1) otherwise
*/
int32_t GeoRef_WD2UV(TGeoRef *Ref,float *uuout,float *vvout,float *spdin,float *wdin,double *Lat,double *Lon,int32_t Nb) {

   int32_t   ni,nj;
   double *lat_true,*lon_true;

   if (Ref->NbSub > 0 ) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   }

   if (!Lat) {
      GeoRef_CalcLL(Ref);
      Lat=Ref->Lat;
      Lon=Ref->Lon;      
   }

   ni = Nb;
   nj = 1;

   if (uuout != spdin) memcpy(uuout, spdin, Nb*sizeof(float));
   if (vvout != wdin)  memcpy(vvout, wdin, Nb*sizeof(float));

   switch (Ref->GRTYP[0]) {
      case 'E':
         lat_true=(double *)(malloc(2*Nb*sizeof(double)));
         lon_true=&lat_true[Nb];
         GeoRef_RotateXY(lat_true,lon_true,Lon,Lat,ni,Ref->RPNHeadExt.xg1,Ref->RPNHeadExt.xg2,Ref->RPNHeadExt.xg3,Ref->RPNHeadExt.xg4);
         c_ezgfwfllw8(uuout,vvout,Lat,Lon,lat_true,lon_true,&ni,&nj,Ref->GRTYP,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4);
         free(lat_true);
         return(0);
         break;

      case '#':
      case 'Y':
      case 'Z':
         switch(Ref->RPNHeadExt.grref[0]) {
            case 'E':
               lat_true=(double *)(malloc(2*Nb*sizeof(double)));
               lon_true=&lat_true[Nb];
               GeoRef_RotateXY(Lat,Lon,lon_true,lat_true,ni,Ref->RPNHeadExt.xgref1,Ref->RPNHeadExt.xgref2,Ref->RPNHeadExt.xgref3,Ref->RPNHeadExt.xgref4);
               c_ezgfwfllw8(uuout,vvout,Lat,Lon,lat_true,lon_true,&ni,&nj,Ref->RPNHeadExt.grref,&Ref->RPNHeadExt.igref1,&Ref->RPNHeadExt.igref2,&Ref->RPNHeadExt.igref3,&Ref->RPNHeadExt.igref4);
               free(lat_true);
               return(0);
               break;

            default:
               f77name(ez8_gdwfllw)(uuout,vvout,Lon,&ni,&nj,Ref->RPNHeadExt.grref,&Ref->RPNHeadExt.igref1,&Ref->RPNHeadExt.igref2,&Ref->RPNHeadExt.igref3,&Ref->RPNHeadExt.igref4, 1);
               break;
         }

      default:
         f77name(ez8_gdwfllw)(uuout,vvout,Lon,&ni,&nj,Ref->GRTYP,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4, 1);
         break;
   }

   return(0);
}

/**----------------------------------------------------------------------------
 * @brief  Converts grid winds (uu/vv) to meteorological winds (speed/direction).

 *    @param[in]  Ref          geo-reference
 *    @param[out] spdout       Destination speed values
 *    @param[out] wdout        Destination direction values
 *    @param[in]  uuin         Source UU values
 *    @param[in]  vvin         Source VV values
 *    @param[in]  Lat          Latitude of wind value points (NULL, use geo-reference gridpoints)
 *    @param[in]  Lon          Longitude of wind value points (NULL, use geo-reference gridpoints)
 *    @param[in]  Nb           Number of values to convert (size of arrays)

 *    @return                FALSE (0) if operation failed, TRUE (1) otherwise
*/
int32_t GeoRef_UV2WD(TGeoRef *Ref,float *spdout,float *wdout,float *uuin,float *vvin,double *Lat,double *Lon,int32_t Nb) {

   int32_t  ni,nj;
   double  *lat_rot,*lon_rot;

   if (Ref->NbSub > 0 ) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   }

   if (!Lat) {
      GeoRef_CalcLL(Ref);
      Lat=Ref->Lat;
      Lon=Ref->Lon;      
   }

   ni = Nb;
   nj = 1;

   if (spdout != uuin) memcpy(spdout, uuin, Nb*sizeof(float));
   if (wdout != vvin)  memcpy(wdout, vvin, Nb*sizeof(float));

   switch (Ref->GRTYP[0]) {
      case 'E':
         lat_rot=(double *)(malloc(2*Nb*sizeof(double)));
         lon_rot=&lat_rot[Nb];
         GeoRef_RotateXY(Lat,Lon,lon_rot,lat_rot,ni,Ref->RPNHeadExt.xg1,Ref->RPNHeadExt.xg2,Ref->RPNHeadExt.xg3,Ref->RPNHeadExt.xg4);
         c_ezllwfgfw8(spdout,wdout,Lat,Lon,lat_rot,lon_rot,&ni,&nj,Ref->GRTYP,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4);
         free(lat_rot);
         return(0);
         break;

      case '#':
      case 'Y':
      case 'Z':
         switch(Ref->RPNHeadExt.grref[0]) {
	         case 'E':
               lat_rot=(double *)(malloc(2*Nb*sizeof(double)));
               lon_rot=&lat_rot[Nb];
	            GeoRef_RotateXY(Lat,Lon,lon_rot,lat_rot,ni,Ref->RPNHeadExt.xgref1,Ref->RPNHeadExt.xgref2,Ref->RPNHeadExt.xgref3,Ref->RPNHeadExt.xgref4);
	            c_ezllwfgfw8(spdout,wdout,Lat,Lon,lat_rot,lon_rot,&ni,&nj,Ref->RPNHeadExt.grref,&Ref->RPNHeadExt.igref1,&Ref->RPNHeadExt.igref2,&Ref->RPNHeadExt.igref3,&Ref->RPNHeadExt.igref4);
	            free(lat_rot);
	            return(0);
	            break;

	         default:
	            f77name(ez8_llwfgdw)(spdout,wdout,Lon,&ni,&nj,Ref->RPNHeadExt.grref,&Ref->RPNHeadExt.igref1,&Ref->RPNHeadExt.igref2,&Ref->RPNHeadExt.igref3,&Ref->RPNHeadExt.igref4,1);
	            break;
	      }
         break;

      default:
         f77name(ez8_llwfgdw)(spdout,wdout,Lon,&ni,&nj,Ref->GRTYP,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4,1);
         break;
   }

   return(0);
}

/**----------------------------------------------------------------------------
 * @brief  Converts grid winds (uu/vv) to geographical components(uu/vv on EW/NS axis).

 *    @param[in]  Ref          geo-reference
 *    @param[out] uullout       Destination speed values
 *    @param[out] vvllout        Destination direction values
 *    @param[in]  uuin         Source UU values
 *    @param[in]  vvin         Source VV values
 *    @param[in]  Lat          Latitude of wind value points (NULL, use geo-reference gridpoints)
 *    @param[in]  Lon          Longitude of wind value points (NULL, use geo-reference gridpoints)
 *    @param[in]  Nb           Number of values to convert (size of arrays)

 *    @return                FALSE (0) if operation failed, TRUE (1) otherwise
*/
int32_t GeoRef_UV2UV(TGeoRef *Ref,float *uullout,float *vvllout,float *uuin,float *vvin,double *Lat,double *Lon,int32_t Nb) {

   int32_t  ni,nj;
   double  *lat_rot,*lon_rot;

   if (Ref->NbSub>0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(FALSE);
   }

   if (!Lat) {
      GeoRef_CalcLL(Ref);
      Lat=Ref->Lat;
      Lon=Ref->Lon;      
   }

   ni = Nb;
   nj = 1;

   if (uullout != uuin) memcpy(uullout, uuin, Nb*sizeof(float));
   if (vvllout != vvin) memcpy(vvllout, vvin, Nb*sizeof(float));

   switch (Ref->GRTYP[0]) {
      case 'E':
         lat_rot = (double*)malloc(2*Nb*sizeof(double));
         lon_rot = &lat_rot[Nb];
         GeoRef_RotateXY(Lat,Lon,lon_rot,lat_rot,ni,Ref->RPNHeadExt.xg1,Ref->RPNHeadExt.xg2,Ref->RPNHeadExt.xg3,Ref->RPNHeadExt.xg4);
         c_ezllwfgff8(uullout,vvllout,Lat,Lon,lat_rot,lon_rot,&ni,&nj,Ref->GRTYP,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4);
         free(lat_rot);
         return(TRUE);
         break;

      case '#':
      case 'Y':
      case 'Z':
         switch (Ref->RPNHeadExt.grref[0]) {
            case 'E':
               lat_rot = (double*)malloc(2*Nb*sizeof(double));
               lon_rot = &lat_rot[Nb];
	            GeoRef_RotateXY(Lat,Lon,lon_rot,lat_rot,ni,Ref->RPNHeadExt.xgref1,Ref->RPNHeadExt.xgref2,Ref->RPNHeadExt.xgref3,Ref->RPNHeadExt.xgref4);
               c_ezllwfgff8(uullout,vvllout,Lat,Lon,lat_rot,lon_rot,&ni,&nj,Ref->RPNHeadExt.grref,&Ref->RPNHeadExt.igref1,&Ref->RPNHeadExt.igref2,&Ref->RPNHeadExt.igref3,&Ref->RPNHeadExt.igref4);
               free(lat_rot);
               return 0;
               break;

            default:
               break;
         }
         break;

      default:
         break;
    }

   return(TRUE);
}

 /*
    gfwfllw -> GeF Winds From LatLon Winds
    Ce sous-programme effectue la rotation des vents d'un systeme de coordonne
    non tourne a un systeme de coordonnee tourne.
    latin, lonin sont les latlons vraies
    xlatingf, xloningf sont les latlons sur la grille tournee
  */
void c_ezgfwfllw8(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int32_t *ni,int32_t *nj,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4) {

   int32_t zero = 0;
   int32_t npts = *ni * *nj;
   int32_t trois = 3;
   float r[9], ri[9], xlon1, xlat1, xlon2, xlat2;
   double *uvcart, *xyz;
   char grtypl[2];

   uvcart = (double *)malloc(2*3*npts*sizeof(double));
   xyz    = &uvcart[3*npts];

   f77name(cigaxg)(grtyp, &xlat1, &xlon1, &xlat2, &xlon2, ig1, ig2, ig3, ig4,1);
   f77name(ez_crot)(r, ri, &xlon1, &xlat1, &xlon2, &xlat2);
   grtypl[0] = 'L';
   f77name(ez8_gdwfllw)(uullout,vvllout,Lon,ni,nj,grtypl, &zero, &zero, &zero, &zero, 1);
   f77name(ez8_uvacart)(xyz, uullout, vvllout, Lon, Lat, ni, nj);
   f77name(ez8_mxm)(r, &trois, xyz, &trois, uvcart, &npts);
   f77name(ez8_cartauv)(uullout, vvllout, uvcart, xloningf, xlatingf, ni, nj);

   free(uvcart);
}

  /*
    llwfgfw -> LatLon Winds From GeF Winds
    Ce sous-programme effectue la rotation des vents d'un systeme de coordonne
    tourne a un systeme de coordonnee non tourne.
    latin, lonin sont les latlons vraies
    xlatingf, xloningf sont les latlons sur la grille tournee
  */
void c_ezllwfgfw8(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int32_t *ni,int32_t *nj,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4) {

   int32_t zero = 0;
   int32_t npts = *ni * *nj;
   int32_t trois = 3;
   float r[9], ri[9], xlon1, xlat1, xlon2, xlat2;
   double *uvcart, *xyz;
   char grtypl[2];

   uvcart = (double *) malloc(2*3*npts*sizeof(double));
   xyz    = &uvcart[3*npts];

   f77name(cigaxg)(grtyp, &xlat1, &xlon1, &xlat2, &xlon2, ig1, ig2, ig3, ig4,1);
   f77name(ez_crot)(r, ri, &xlon1, &xlat1, &xlon2, &xlat2);
   f77name(ez8_uvacart)(xyz, uullout, vvllout, xloningf, xlatingf, ni, nj);
   f77name(ez8_mxm)(ri, &trois, xyz, &trois, uvcart, &npts);
   f77name(ez8_cartauv)(uullout, vvllout, uvcart, Lon, Lat, ni, nj);
   grtypl[0] = 'L';
   f77name(ez8_llwfgdw)(uullout,vvllout,xloningf,ni,nj,grtypl,&zero,&zero,&zero,&zero,1);

   free(uvcart);
}

  /*
    llwfgff -> LatLon Winds From GeF Winds
    Ce sous-programme effectue la rotation des vents d'un systeme de coordonne
    tourne a un systeme de coordonnee non tourne.
    latin, lonin sont les latlons vraies
    xlatingf, xloningf sont les latlons sur la grille tournee
  */
void c_ezllwfgff8(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int32_t *ni,int32_t *nj,char *grtyp,int32_t *ig1,int32_t *ig2,int32_t *ig3,int32_t *ig4) {

   int32_t npts = *ni * *nj;
   int32_t trois = 3;
   float r[9], ri[9], xlon1, xlat1, xlon2, xlat2;
   double *uvcart, *xyz;

   uvcart = (double *) malloc(2*3*npts*sizeof(double));
   xyz    = &uvcart[3*npts];

   f77name(cigaxg)(grtyp, &xlat1, &xlon1, &xlat2, &xlon2, ig1, ig2, ig3, ig4,1);
   f77name(ez_crot)(r, ri, &xlon1, &xlat1, &xlon2, &xlat2);
   f77name(ez8_uvacart)(xyz, uullout, vvllout, xloningf, xlatingf, ni, nj);
   f77name(ez8_mxm)(ri, &trois, xyz, &trois, uvcart, &npts);
   f77name(ez8_cartauv)(uullout, vvllout, uvcart, Lon, Lat, ni, nj);

   free(uvcart);
}
