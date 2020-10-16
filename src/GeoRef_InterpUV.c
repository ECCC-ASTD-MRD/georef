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

void c_ezgfwfllw(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int *ni,int *nj,char *grtyp,int *ig1,int *ig2,int *ig3,int *ig4);
void c_ezllwfgfw(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int *ni,int *nj,char *grtyp,int *ig1,int *ig2,int *ig3,int *ig4);

int GeoRef_InterpUV(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin) {
   
   int npts,ier,ierc,ierc1,ierc2;
   float *uullout = NULL;
   float *vvllout = NULL;

   GeoRef_SetGet(RefTo,RefFrom);

   if (RefFrom->NbSub > 0 || RefTo->NbSub > 0) {
      ier = GeoRef_InterpYYUV(RefTo,RefFrom,uuout,vvout,uuin,vvin);
   } else {

      ierc = 0;
      ierc1 = 0;
      ierc2 = 0;
   
      npts = RefTo->NX*RefTo->NY;
      ier = GeoRef_CalcLL(RefTo);
   
      RefFrom->Options.VectorMode = TRUE;
      RefFrom->Options.Symmetric = TRUE;

      ierc1=GeoRef_Interp(RefTo,RefFrom,uuout,uuin);
      RefFrom->Options.Symmetric = FALSE;
      ierc2=GeoRef_Interp(RefTo,RefFrom,vvout,vvin);
      RefFrom->Options.Symmetric = TRUE;
      if (ierc1 == 2 || ierc2 == 2) {
         ierc = 2;
      }
   
      if (RefFrom->Options.PolarCorrect == TRUE) {
         ier=GeoRef_CorrectVector(RefTo,RefFrom,uuout,vvout,uuin,vvin);
      }
   
      uullout = (float *)malloc(2*npts*sizeof(float));
      vvllout = &uullout[npts];
   
      GeoRef_UV2WD(RefFrom,uullout,vvllout,uuout,vvout,RefTo->Lat,RefTo->Lon,npts);
      GeoRef_WD2UV(RefTo,uuout,vvout,uullout,vvllout,RefTo->Lat,RefTo->Lon,npts);
   
      RefFrom->Options.VectorMode = FALSE;
      free(uullout);   
   }
   
   return(ierc);
}

int GeoRef_InterpWD(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin) {

   int ier,ierc,ierc1,ierc2;
   float *uullout = NULL;
   float *vvllout = NULL;
   int npts;

   GeoRef_SetGet(RefTo,RefFrom);

   if (RefFrom->NbSub > 0 || RefTo->NbSub > 0) {
      ierc = GeoRef_InterpYYWD(RefTo,RefFrom,uuout,vvout,uuin,vvin);
   } else {
      ierc = 0;
      ierc1 = 0;
      ierc2 = 0;

      npts = RefTo->NX*RefTo->NY;

      RefFrom->Options.VectorMode = TRUE;
      RefFrom->Options.Symmetric = TRUE;

      ierc1=GeoRef_Interp(RefTo,RefFrom,uuout,uuin);
      RefFrom->Options.Symmetric = FALSE;
      ierc2=GeoRef_Interp(RefTo,RefFrom,vvout,vvin);

      if (ierc1 == 2 || ierc2 == 2) {
         ierc = 2;
      }

      RefFrom->Options.Symmetric = TRUE;

      if (RefFrom->Options.PolarCorrect == TRUE) {
         ier = GeoRef_CorrectVector(RefTo,RefFrom,uuout,vvout,uuin,vvin);
      }

      uullout = (float *)malloc(2*npts*sizeof(float));
      vvllout = &uullout[npts];

      /*ezsint does not allocate lat,lon if RefFrom=RefTo*/
      ier = GeoRef_CalcLL(RefTo);

      GeoRef_UV2WD(RefFrom,uullout,vvllout,uuout,vvout,RefTo->Lat,RefTo->Lon,npts);

      memcpy(uuout, uullout, npts*sizeof(float));
      memcpy(vvout, vvllout, npts*sizeof(float));

      RefFrom->Options.VectorMode = FALSE;
      free(uullout);
   }
   return(ierc);
}

int GeoRef_InterpYYUV(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin) {

   TGeoRef *yin_gdin, *yan_gdin, *yin_gdout, *yan_gdout;
   TGridSet *gset=NULL;
   int icode,i,j,k,ierc1,ierc2,ierc;
   int yancount_yin,yincount_yin, yancount_yan,yincount_yan;
   int yyin,yyout;
   int ni, nj,idx;
   float *yin2yin_uuout,*yan2yin_uuout, *yin2yin_vvout,*yan2yin_vvout;
   float *yin2yan_uuout,*yan2yan_uuout, *yin2yan_vvout,*yan2yan_vvout;
   float *yin2yin_spdout,*yan2yin_spdout, *yin2yin_wdout,*yan2yin_wdout;
   float *yin2yan_spdout,*yan2yan_spdout, *yin2yan_wdout,*yan2yan_wdout;
   float *spdout,*wdout;
  
   // Need only access to either yin or Yang info for the lat and lon val 
   
   yyin=0; yyout=0;
   ierc=0;
   ierc1=0;ierc2=0;

   gset=GeoRef_SetGet(RefTo,RefFrom);

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
      ierc1 = GeoRef_InterpUV(yin_gdout,RefFrom,uuout,vvout,uuin,vvin);
      ierc2 = GeoRef_InterpUV(yan_gdout,RefFrom,&uuout[ni*nj],&vvout[ni*nj],uuin,vvin);
      if (ierc1 == 2 || ierc2 == 2) {
         ierc=2;
      }
      return(ierc);
   }

   // Check if one sub grid is identical to one of the sub grids

   if (yin_gdin == RefTo) {
      icode = GeoRef_InterpUV(RefTo,yin_gdin,uuout,vvout,uuin,vvin);
      return(icode);
   }
   if (yan_gdin == RefTo) {
      icode = GeoRef_InterpUV(RefTo,yan_gdin,uuout,vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)]);
      return(icode);
   }

   /* User specifies to use 1 subgrid for interpolation ezsetopt(USE_1SUBGRID) */
   /* User must specify one specific grid ezsetival(SUBGRIDID) */
   /* This is only appropriate if the destination grid is non yin-yang grid */
 
   if (RefFrom->Options.SubGrid) { // User specifies to use 1 grid only
      // Output is a Yin-Yang grid 
      if (yyout == 1) { 
         App_Log(ERROR,"%s: Cannot use subgrid to interpolate to a Yin-Yang grid\n",__func__);
         return(-1);
      }
      // Is specified subgrid within the subgrid list
      if (RefFrom->Options.SubGrid>RefFrom->NbSub) { 
         App_Log(ERROR,"%s: Invalid subgrid: %i\n",__func__,RefFrom->Options.SubGrid);
         return(-1);
      }
      // Use yin input grid
      if (RefFrom->Options.SubGrid==1) { 
         ierc = GeoRef_InterpUV(yin_gdout,yin_gdin,uuout,vvout,uuin,vvin);
         return(ierc);
      }

      // Use yang input grid
      if (RefFrom->Options.SubGrid==2) { 
         ierc = GeoRef_InterpUV(yin_gdout,yan_gdin,uuout,vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)]);
         return(ierc);
      }
   }

   /* To use both Yin and Yang grids in Yin-yang input grid */
   /* Masquer les grilles YY input pour enlever overlap et calculer les X,Y */
   icode = GeoRef_SetCalcYYXY(RefTo,RefFrom);

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

      icode = GeoRef_XYUVVal(yin_gdin,yin2yin_uuout,yin2yin_vvout,uuin,vvin,gset->yin2yin_x,gset->yin2yin_y,gset->yincount_yin);
      icode = GeoRef_UV2WD(yin_gdin,yin2yin_spdout,yin2yin_wdout,yin2yin_uuout,yin2yin_vvout,gset->yin2yin_lat,gset->yin2yin_lon,yincount_yin);
      icode = GeoRef_XYUVVal(yan_gdin,yan2yin_uuout,yan2yin_vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yin_x,gset->yan2yin_y,yancount_yin);
      icode = GeoRef_UV2WD(yan_gdin,yan2yin_spdout,yan2yin_wdout,yan2yin_uuout, yan2yin_vvout,gset->yan2yin_lat,gset->yan2yin_lon,yancount_yin);
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
      icode = GeoRef_WD2UV(RefTo,uuout,vvout,spdout,wdout,gset->yinlat,gset->yinlon,ni*nj);
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

      icode = GeoRef_XYUVVal(yin_gdin,yin2yin_uuout,yin2yin_vvout,uuin,vvin,gset->yin2yin_x,gset->yin2yin_y,yincount_yin);
      icode = GeoRef_UV2WD(yin_gdin,yin2yin_spdout,yin2yin_wdout,yin2yin_uuout,yin2yin_vvout,gset->yin2yin_lat,gset->yin2yin_lon,yincount_yin);
      icode = GeoRef_XYUVVal(yan_gdin,yan2yin_uuout,yan2yin_vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yin_x,gset->yan2yin_y,yancount_yin);
      icode = GeoRef_UV2WD(yan_gdin,yan2yin_spdout,yan2yin_wdout,yan2yin_uuout,yan2yin_vvout,gset->yan2yin_lat,gset->yan2yin_lon,yancount_yin);

      icode = GeoRef_XYUVVal(yin_gdin,yin2yan_uuout,yin2yan_vvout,uuin,vvin,gset->yin2yan_x,gset->yin2yan_y,yincount_yan);
      icode = GeoRef_UV2WD(yin_gdin,yin2yan_spdout,yin2yan_wdout,yin2yan_uuout,yin2yan_vvout,gset->yin2yan_lat,gset->yin2yan_lon,yincount_yan);

      icode = GeoRef_XYUVVal(yan_gdin,yan2yan_uuout,yan2yan_vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yan_x,gset->yan2yan_y,yancount_yan);
      icode = GeoRef_UV2WD(yan_gdin,yan2yan_spdout,yan2yan_wdout,yan2yan_uuout,yan2yan_vvout,gset->yan2yan_lat,gset->yan2yan_lon,yancount_yan);

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
      icode = GeoRef_WD2UV(yin_gdout,uuout,vvout,spdout,wdout,gset->yinlat,gset->yinlon,ni*nj);

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

      icode = GeoRef_WD2UV(yan_gdout,&uuout[ni*nj],&vvout[ni*nj],spdout,wdout,gset->yanlat,gset->yanlon,ni*nj);
      free(spdout);
   }

   return(icode);
}

int GeoRef_InterpYYWD(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin) {

   TGridSet *gset;
   TGeoRef *yin_gdin, *yan_gdin, *yin_gdout, *yan_gdout;
   int icode,i,j,k,ierc1,ierc2,ierc;
   int yancount_yin,yincount_yin, yancount_yan,yincount_yan;
   int yyin,yyout;
   int ni, nj,idx;
   float *yin2yin_uuout,*yan2yin_uuout, *yin2yin_vvout,*yan2yin_vvout;
   float *yin2yan_uuout,*yan2yan_uuout, *yin2yan_vvout,*yan2yan_vvout;
   float *yin2yin_spdout,*yan2yin_spdout, *yin2yin_wdout,*yan2yin_wdout;
   float *yin2yan_spdout,*yan2yan_spdout, *yin2yan_wdout,*yan2yan_wdout;
   float *spdout,*wdout;
  
   // Need only access to either yin or Yang info for the lat and lon val
   
   yyin=0; yyout=0;
   ierc=0;
   ierc1=0;ierc2=0;

   gset=GeoRef_SetGet(RefTo,RefFrom);

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
      ierc1 = GeoRef_InterpWD(yin_gdout,RefFrom,uuout,vvout,uuin,vvin);
      ierc2 = GeoRef_InterpWD(yan_gdout,RefFrom,&uuout[(ni*nj)],&vvout[(ni*nj)],uuin,vvin);
      if (ierc1 == 2 || ierc2 == 2) {
        ierc=2;
      }
      return(ierc);
   }

   // Check if one sub grid is identical to one of the sub grids 
   if (yin_gdin == RefTo) {
      icode = GeoRef_InterpWD(RefTo,yin_gdin,uuout,vvout,uuin,vvin);
      return(icode);
   }
   if (yan_gdin == RefTo) {
      icode = GeoRef_InterpWD(RefTo,yan_gdin,uuout,vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)]);
      return icode;
   }

   /* User specifies to use 1 subgrid for interpolation ezsetopt(USE_1SUBGRID) */
   /* User must specify one specific grid ezsetival(SUBGRIDID) */
   /* This is only appropriate if the destination grid is non yin-yang grid */

   if (RefFrom->Options.SubGrid) { // User specifies to use 1 grid only
      // Output is a Yin-Yang grid 
      if (yyout == 1) { 
         App_Log(ERROR,"%s: Cannot use subgrid to interpolate to a Yin-Yang grid\n",__func__);
         return(-1);
      }
      // Is specified subgrid within the subgrid list
      if (RefFrom->Options.SubGrid>RefFrom->NbSub) { 
         App_Log(ERROR,"%s: Invalid subgrid: %i\n",__func__,RefFrom->Options.SubGrid);
         return(-1);
      }
      // Use yin input grid
      if (RefFrom->Options.SubGrid==1) { 
         ierc = GeoRef_InterpWD(yin_gdout,yin_gdin,uuout,vvout,uuin,vvin);
         return(ierc);
      }

      // Use yang input grid
      if (RefFrom->Options.SubGrid==2) { 
         ierc = GeoRef_InterpWD(yin_gdout,yan_gdin,uuout,vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)]);
         return(ierc);
      }
   }

   /* To use both Yin and Yang grids in Yin-yang input grid */
   /* Masquer les grilles YY input pour enlever overlap et calculer les X,Y */
   icode = GeoRef_SetCalcYYXY(RefTo,RefFrom);

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

      icode = GeoRef_XYUVVal(yin_gdin,yin2yin_uuout,yin2yin_vvout,uuin,vvin,gset->yin2yin_x,gset->yin2yin_y,gset->yincount_yin);
      icode = GeoRef_UV2WD(yin_gdin,yin2yin_spdout,yin2yin_wdout,yin2yin_uuout,yin2yin_vvout,gset->yin2yin_lat,gset->yin2yin_lon,yincount_yin);
  
      icode = GeoRef_XYUVVal(yan_gdin,yan2yin_uuout,yan2yin_vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yin_x,gset->yan2yin_y,yancount_yin);
      icode = GeoRef_UV2WD(yan_gdin,yan2yin_spdout,yan2yin_wdout,yan2yin_uuout, yan2yin_vvout,gset->yan2yin_lat,gset->yan2yin_lon,yancount_yin);
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

      icode = GeoRef_XYUVVal(yin_gdin,yin2yin_uuout,yin2yin_vvout,uuin,vvin,gset->yin2yin_x,gset->yin2yin_y,gset->yincount_yin);
      icode = GeoRef_UV2WD(yin_gdin,yin2yin_spdout,yin2yin_wdout,yin2yin_uuout,yin2yin_vvout,gset->yin2yin_lat,gset->yin2yin_lon,gset->yincount_yin);
      icode = GeoRef_XYUVVal(yan_gdin,yan2yin_uuout,yan2yin_vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yin_x,gset->yan2yin_y,gset->yancount_yin);
      icode = GeoRef_UV2WD(yan_gdin,yan2yin_spdout,yan2yin_wdout,yan2yin_uuout,yan2yin_vvout,gset->yan2yin_lat,gset->yan2yin_lon,yancount_yin);

      icode = GeoRef_XYUVVal(yin_gdin,yin2yan_uuout,yin2yan_vvout,uuin,vvin,gset->yin2yan_x,gset->yin2yan_y,gset->yincount_yan);
      icode = GeoRef_UV2WD(yin_gdin,yin2yan_spdout,yin2yan_wdout,yin2yan_uuout,yin2yan_vvout,gset->yin2yan_lat,gset->yin2yan_lon,gset->yincount_yan);

      icode = GeoRef_XYUVVal(yan_gdin,yan2yan_uuout,yan2yan_vvout,&uuin[(yin_gdin->NX)*(yin_gdin->NY)],&vvin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yan_x,gset->yan2yan_y,gset->yancount_yan);
      icode = GeoRef_UV2WD(yan_gdin,yan2yan_spdout,yan2yan_wdout,yan2yan_uuout,yan2yan_vvout,gset->yan2yan_lat,gset->yan2yan_lon,gset->yancount_yan);

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

   return(icode);
}

int GeoRef_WD2UV(TGeoRef *Ref,float *uugdout,float *vvgdout,float *uullin,float *vvllin,double *Lat,double *Lon,int Nb) {

   int   ni,nj;
   double *lat_true,*lon_true;
  
   if (Ref->NbSub > 0 ) {
      App_Log(ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   }
   ni = Nb;
   nj = 1;

   memcpy(uugdout, uullin, Nb*sizeof(ftnfloat));
   memcpy(vvgdout, vvllin, Nb*sizeof(ftnfloat));
 
//TODO: Convert double
   switch (Ref->GRTYP[0]) {
      case 'E':
         lat_true=(double *)(malloc(2*Nb*sizeof(double)));
         lon_true=&lat_true[Nb];
         f77name(ez8_gfxyfll)(lon_true,lat_true,Lon,Lat,&ni,&Ref->RPNHead.XG[X_LAT1],&Ref->RPNHead.XG[X_LON1],&Ref->RPNHead.XG[X_LAT2],&Ref->RPNHead.XG[X_LON2]);
         c_ezgfwfllw(uugdout,vvgdout,Lat,Lon,lat_true,lon_true,&ni,&nj,Ref->GRTYP,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);
         free(lat_true);
         return(0);
         break;   
        
      case '#':
      case 'Y':
      case 'Z':
         switch(Ref->RPNHead.GRREF[0]) {
            case 'E':
               lat_true=(double *)(malloc(2*Nb*sizeof(double)));
               lon_true=&lat_true[Nb];
               f77name(ez8_gfxyfll)(Lon,Lat,lon_true,lat_true,&ni,&Ref->RPNHead.XGREF[X_LAT1],&Ref->RPNHead.XGREF[X_LON1],&Ref->RPNHead.XGREF[X_LAT2],&Ref->RPNHead.XGREF[X_LON2]);         
               c_ezgfwfllw(uugdout,vvgdout,Lat,Lon,lat_true,lon_true,&ni,&nj,Ref->RPNHead.GRREF,&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4]);
               free(lat_true);
               return(0);
               break;
	    
            default:
               f77name(ez8_gdwfllw)(uugdout,vvgdout,Lon,&ni,&nj,Ref->RPNHead.GRREF,&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4], 1);
               break;
         }
        
      default:
         f77name(ez8_gdwfllw)(uugdout,vvgdout,Lon,&ni,&nj,Ref->GRTYP,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4], 1);
         break;
   }
   
   return(0);
}

int GeoRef_UV2WD(TGeoRef *Ref,float *spd_out,float *wd_out,float *uuin,float *vvin,double *Lat,double *Lon,int Nb) {

   int    ni,nj;
   double *lat_rot,*lon_rot;
   
   if (Ref->NbSub > 0 ) {
      App_Log(ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   }

   ni = Nb;
   nj = 1;

   memcpy(spd_out, uuin, Nb*sizeof(ftnfloat));
   memcpy(wd_out, vvin, Nb*sizeof(ftnfloat));

//TODO: Convert double
   switch (Ref->GRTYP[0]) {
      case 'E':
         lat_rot=(double *)(malloc(2*Nb*sizeof(double)));
         lon_rot=&lat_rot[Nb];
         f77name(ez8_gfxyfll)(Lon,Lat,lon_rot,lat_rot,&ni,&Ref->RPNHead.XG[X_LAT1],&Ref->RPNHead.XG[X_LON1],&Ref->RPNHead.XG[X_LAT2],&Ref->RPNHead.XG[X_LON2]);
         c_ezllwfgfw(spd_out,wd_out,Lat,Lon,lat_rot,lon_rot,&ni,&nj,Ref->GRTYP,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);
         free(lat_rot);
         return(0);  
         break;
       
      case '#':
      case 'Y':
      case 'Z':
         switch(Ref->RPNHead.GRREF[0]) {
	         case 'E':
               lat_rot=(double *)(malloc(2*Nb*sizeof(double)));
               lon_rot=&lat_rot[Nb];
	            f77name(ez8_gfxyfll)(Lon,Lat,lon_rot,lat_rot,&ni,&Ref->RPNHead.XGREF[X_LAT1],&Ref->RPNHead.XGREF[X_LON1],&Ref->RPNHead.XGREF[X_LAT2],&Ref->RPNHead.XGREF[X_LON2]);
	            c_ezllwfgfw(spd_out,wd_out,Lat,Lon,lat_rot,lon_rot,&ni,&nj,Ref->RPNHead.GRREF,&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4]);
	            free(lat_rot);
	            return(0);
	            break;
	   
	         default:
	            f77name(ez8_llwfgdw)(spd_out,wd_out,Lon,&ni,&nj,Ref->RPNHead.GRREF,&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4],1);
	            break;
	      }
         break;
       
      default:
         f77name(ez8_llwfgdw)(spd_out,wd_out,Lon,&ni,&nj,Ref->GRTYP,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4],1);
         break;
   }
      
   return(0);
}

 /*
    gfwfllw -> GeF Winds From LatLon Winds
    Ce sous-programme effectue la rotation des vents d'un systeme de coordonne
    non tourne a un systeme de coordonnee tourne.
    latin, lonin sont les latlons vraies
    xlatingf, xloningf sont les latlons sur la grille tournee
  */

void c_ezgfwfllw(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int *ni,int *nj,char *grtyp,int *ig1,int *ig2,int *ig3,int *ig4) {

   int zero = 0;
   int npts = *ni * *nj;
   int trois = 3;
   float r[9], ri[9], xlon1, xlat1, xlon2, xlat2;
   double *uvcart, *xyz;
   char grtypl[2];

   uvcart = (double *)malloc(2*3*npts*sizeof(double));
   xyz    = &uvcart[3*npts];
  
   f77name(cigaxg)(grtyp, &xlat1, &xlon1, &xlat2, &xlon2, ig1, ig2, ig3, ig4);
   f77name(ez_crot)(r, ri, &xlon1, &xlat1, &xlon2, &xlat2);
   grtypl[0] = 'L';
   f77name(ez8_gdwfllw)(uullout,vvllout,Lon,ni,nj,grtypl, &zero, &zero, &zero, &zero, 1);
   f77name(ez8_uvacart)(xyz, uullout, vvllout, Lon, Lat, ni, nj);
   // TODO Maude: where is this function
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
void c_ezllwfgfw(float *uullout,float *vvllout,double *Lat,double *Lon,double *xlatingf,double *xloningf,int *ni,int *nj,char *grtyp,int *ig1,int *ig2,int *ig3,int *ig4) {

   int zero = 0;
   int npts = *ni * *nj;
   int trois = 3;
   float r[9], ri[9], xlon1, xlat1, xlon2, xlat2;
   double *uvcart, *xyz;
   char grtypl[2];

 
   uvcart = (double *) malloc(2*3*npts*sizeof(double));
   xyz    = &uvcart[3*npts];  
  
   f77name(cigaxg)(grtyp, &xlat1, &xlon1, &xlat2, &xlon2, ig1, ig2, ig3, ig4);
   f77name(ez_crot)(r, ri, &xlon1, &xlat1, &xlon2, &xlat2);
   f77name(ez8_uvacart)(xyz, uullout, vvllout, xloningf, xlatingf, ni, nj); 
   f77name(ez8_mxm)(ri, &trois, xyz, &trois, uvcart, &npts);
   f77name(ez8_cartauv)(uullout, vvllout, uvcart, Lon, Lat, ni, nj); 
   grtypl[0] = 'L';
   f77name(ez8_llwfgdw)(uullout,vvllout,xloningf,ni,nj,grtypl, &zero, &zero, &zero, &zero, 1);

   free(uvcart);
}