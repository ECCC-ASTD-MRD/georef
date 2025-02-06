#include <App.h>
#include "GeoRef.h"
#include "Triangle.h"

int32_t c_gd_isgridrotated(TGeoRef *gr) {

   if (gr->RPNHeadExt.grref[0] == 'E') {
       if (fabs(gr->RPNHeadExt.xgref1-gr->RPNHeadExt.xgref3) < 0.001)
         return(0); // non rotated
      else
         return(1); // rotated
   } else {
      return(0);
   }
   return(0);
}

int32_t c_gdcompatible_grids(TGeoRef *RefFrom, TGeoRef* RefTo) {

   switch(RefTo->GRTYP[0]) {
	   case 'L':
	   case 'A':
	   case 'B':
	   case 'G':
	      return(0);

      // This is a fix from previous version that was never reached
      case 'Z':
         if (RefFrom->RPNHeadExt.grref[0] == 'L')
            return(0);
         else if (RefFrom->RPNHeadExt.grref[0] == 'E')
            if (!c_gd_isgridrotated(RefFrom))
               return(0);
            else
               return(-1);
         break;

	   default:
	      return(-1);
	}

   return 0;
}

int32_t gd_interpm(TGeoRef *Ref,TGeoOptions *Opt,float *Out,float *In,double *X,double *Y,int32_t Nb) {

   Vect3d  b,v;
   int32_t d,n,idx;

   #pragma omp parallel for default(none) private(d,b,idx,n,v) shared(Nb,Ref,Opt,X,Y,Out,In)
   for(d=0;d<Nb;d++) {
      if (X[d]>=0 && Y[d]>=0) {
         b[0]=X[d]-(int)X[d];
         b[1]=Y[d]-(int)Y[d];
         b[2]=1.0-b[0]-b[1];
         idx=(int)X[d]-1;

         if(idx>Ref->NIdx) {
            Out[d]=Opt->NoData;
            continue;
         }

         if (Opt->Interp==IR_NEAREST) {
            n=Bary_Nearest(b);
            Out[d]=In[Ref->Idx[idx+n]];
         } else {
            v[0]=In[Ref->Idx[idx]];
            v[1]=In[Ref->Idx[idx+1]];
            v[2]=In[Ref->Idx[idx+2]];

            Out[d]=Bary_InterpV(b,v);
          }
      }
   }
   return(Nb);
}
         
int32_t GeoRef_InterpFinally(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout,float *zin,double *X,double *Y,int32_t npts,TGeoSet *GSet) {

   int32_t ier, un, j;
   int32_t old_degre_interp;
   double *gdst_lats, tmp, real_un=1.0, real_j, *x,*y;

   // RefTo needed for type 4,5 and Y grid
   // TODO: Check old type 4 and,
   // TODO: need to use new types (do not forget subllinear and subnearest)
   if ((!RefTo && (Opt->Interp==4 || Opt->Interp==5)) || !RefFrom) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid georeference\n",__func__);
      return(-1);
   }

   if (!X || !Y) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Local coordinates not available\n",__func__);
      return(-1);
   }

   if (Opt->CIndex) {
      x=(double*)malloc(npts*sizeof(double)*2);
      y=&x[npts];

      for(j=0;j<npts;j++) {
         x[j]=X[j]+1.0;
         y[j]=Y[j]+1.0;
      }
   } else {
      x=X;
      y=Y;
   }

   old_degre_interp = Opt->Interp;

   switch (Opt->Interp) {
      //TODO: check 4 and 5 (DISTANCE TRIANGLE)
      case IR_AVERAGE:
      case 5:
         ier = c_gdcompatible_grids(RefFrom, RefTo);
         if (ier < 0) {
            Lib_Log(APP_LIBGEOREF,APP_WARNING,"%s: Input and output grids are not compatible for average computation, interpolaton level set to linear\n",__func__);
            Opt->Interp = IR_LINEAR;
         }
         break;

      default:
         break;
   }

   switch(RefFrom->GRTYP[0]) {
      case 'M':
         gd_interpm(RefFrom,Opt,zout,zin,x,y,npts);
         break;
      case '#':
      case 'Z':
      case 'G':
         switch (Opt->Interp) {
            case IR_NEAREST:
               if (!GSet) {
                  f77name(ez8_rgdint_0)(zout,x,y,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&Opt->NoData);
               } else {
                 if (GeoRef_SetEmptyIndex(GSet)) {
                    f77name(ez8_rgd_index_0)(GSet->Index,x,y,&npts,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2);
                 }
                 f77name(ez8_apply_0)(GSet->Index,zout,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&Opt->NoData);
               }
               break;

            case IR_LINEAR:
               if (!GSet) {
                  f77name(ez8_irgdint_1)(zout,x,y,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,RefFrom->AX,RefFrom->AY,&RefFrom->Extension,&Opt->NoData);
               } else {
                  if (GeoRef_SetEmptyIndex(GSet)) {
                     f77name(ez8_irgd_index_1)(GSet->Index,x,y,&npts,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,RefFrom->AX,RefFrom->AY,&RefFrom->Extension);
                  }
                  f77name(ez8_apply_1)(GSet->Index,zout,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&Opt->NoData);
               }
               break;

            case IR_CUBIC:
               f77name(ez8_irgdint_3)(zout,x,y,&npts,zin,&RefFrom->NX,&RefFrom->i1,&RefFrom->i2,&RefFrom->j1,&RefFrom->j2,RefFrom->AX,RefFrom->AY,RefFrom->NCX,RefFrom->NCY,&RefFrom->Extension,&Opt->NoData);
               break;

            case 4:
               f77name(ez_avg)(zout,x,y,&RefTo->NX,&RefTo->NY,zin,&RefFrom->NX,&RefFrom->NY,&RefFrom->Extension);
               break;

            case 5:
               gdst_lats = (double*) malloc(sizeof(double)*npts);
               for (j=0; j < RefTo->NY; j++) {
                  real_j = 1.0 * (j+1);
                  ier = GeoRef_XY2LL(RefTo,&gdst_lats[j],&tmp,&real_un,&real_j,1,TRUE);
               }
               f77name(ez_avg_sph)(zout,x,y,gdst_lats,&RefTo->NX,&RefTo->NY,zin,&RefFrom->NX,&RefFrom->NY,&RefFrom->Extension);
               break;
         }
         break;

      case 'Y':
         un = 1;
         if (RefFrom->NX > 1 && RefFrom->NY > 1 && Opt->Interp==IR_LINEAR) {
            f77name(ez8_rgdint_1)(zout,x,y,&npts,zin,&RefFrom->NX,&un,&RefFrom->NY,0,&Opt->NoData);
         } else {
            if (!GSet) {
               Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: GeoSet not defined\n",__func__);
               return(-1);
            }
            f77name(ez_applywgts)(zout,GSet->wts,GSet->idx,zin,GSet->mask,&RefFrom->NX,&RefFrom->NY,&RefTo->NX,&RefTo->NY,&(GSet->n_wts),&Opt->NoData);
         }
         break;

      default:
         switch (Opt->Interp) {
            case IR_NEAREST:
               if (!GSet) {
                  f77name(ez8_rgdint_0)(zout,x,y,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&Opt->NoData);
               } else {
                  if (GeoRef_SetEmptyIndex(GSet)) {
                     f77name(ez8_rgd_index_0)(GSet->Index,x,y,&npts,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2);
                  }
                  f77name(ez8_apply_0)(GSet->Index,zout,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&Opt->NoData);
               }
               break;

            case IR_LINEAR:
               if (!GSet) {
                  f77name(ez8_rgdint_1)(zout,x,y,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&RefFrom->Extension,&Opt->NoData);
               } else {
                  if (GeoRef_SetEmptyIndex(GSet)) {
                     f77name(ez8_rgd_index_1)(GSet->Index,x,y,&npts,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&RefFrom->Extension);
                  }
                  f77name(ez8_apply_1)(GSet->Index,zout,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&Opt->NoData);
               }
               break;

            case IR_CUBIC:
               if (!GSet) {
                  f77name(ez8_rgdint_3)(zout,x,y,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&RefFrom->Extension,&Opt->NoData);
               } else {
                  if (GeoRef_SetEmptyIndex(GSet)) {
                     f77name(ez8_rgd_index_3)(GSet->Index,x,y,&npts,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&RefFrom->Extension);
                  }
                  f77name(ez8_apply_3)(GSet->Index,zout,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&Opt->NoData);
               }
               break;

            case 4:
               f77name(ez_avg)(zout,x,y,&RefTo->NX,&RefTo->NY,zin,&RefFrom->NX,&RefFrom->NY,&RefFrom->Extension);
               break;

            case 5:
               gdst_lats = (double*) malloc(sizeof(double)*npts);
               for (j=0; j < RefTo->NY; j++) {
                  real_j = 1.0 * (j+1);
                  ier = GeoRef_XY2LL(RefTo,&gdst_lats[j],&tmp,&real_un,&real_j,1,TRUE);
               }
               f77name(ez_avg_sph)(zout,x,y,gdst_lats,&RefTo->NX,&RefTo->NY,zin,&RefFrom->NX,&RefFrom->NY,&RefFrom->Extension);
               free(gdst_lats);
               break;
         }
         break;
   }

   if (Opt->CIndex) {
      free(x);
   }
   Opt->Interp = old_degre_interp;
   return(0);
}

/**----------------------------------------------------------------------------
 * @brief  Interpolates values between 2 georeferences
 *    @param[in]  RefTo      Destination geo-reference
 *    @param[in]  RefFrom    Source geo-reference
 *    @param[in]  Opt        Interpolation options
 *    @param[out] zout       Destination interpolated values
 *    @param[in]  zin        Source values

 *    @return                FALSE (0) if operation failed, TRUE (1) otherwise
*/
int32_t GeoRef_Interp(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout,float *zin) {

   TGeoSet    *gset=NULL;
   TApp_Timer *int_timer = App_TimerCreate();
   float      *lzin,*lxzin;
   int32_t         ok=TRUE;

   lzin  = NULL;
   lxzin = NULL;

   if (!Opt) Opt=&RefTo->Options;
   if (!Opt) Opt=&GeoRef_Options;

   if (!RefFrom || !RefTo) {
      Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Source or target grid undefined\n",__func__);
      return(FALSE);
   }

   if (RefFrom == RefTo) {
      memcpy(zout,zin,RefFrom->NX*RefFrom->NY*sizeof(float));
      return(TRUE);
   }

   App_TimerStart(int_timer);

   if (RefFrom->NbSub > 0 || RefTo->NbSub > 0) {
      // YY mutli grids involved
      return(GeoRef_InterpYY(RefTo,RefFrom,Opt,zout,zin));
   } else {
      gset=GeoRef_SetGet(RefTo,RefFrom,Opt);

      if (RefFrom->Type&GRID_YINVERT) {
         lzin = (float *)malloc(RefFrom->NX * RefFrom->NY * sizeof(float));
         memcpy(lzin,zin,RefFrom->NX*RefFrom->NY*sizeof(float));
         f77name(permut)(lzin,&RefFrom->NX,&RefFrom->NY);
      } else {
         lzin = zin;
      }

      if (RefFrom->Type&GRID_EXPAND) {
         lxzin = (float *)malloc(2*RefFrom->NX*RefFrom->NY*sizeof(float));
         GeoRef_GridGetExpanded(RefFrom,Opt,lxzin,lzin);
      } else {
         lxzin = lzin;
      }

      if (GeoRef_CalcLL(RefTo)) {
         GeoRef_SetCalcXY(gset);
         if (GeoRef_InterpFinally(RefTo,RefFrom,Opt,zout,lxzin,gset->X,gset->Y,RefTo->NX*RefTo->NY,gset)==0) {
            if (Opt->PolarCorrect) {
               GeoRef_SetZoneDefine(gset);
               GeoRef_CorrectValue(gset,zout,lxzin);
            }
         } else {
            ok=FALSE;
         }
      }

      if (lzin && lzin!=zin) {
         free(lzin);
      }

      if (lxzin && lxzin!=lzin && lxzin!=zin) {
         free(lxzin);
      }
   }
   App_TimerStop(int_timer);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Interpolation took \033[1;32m%.3f ms\033[0m\n",__func__,App_TimerTotalTime_ms(int_timer));

   return(ok);
}

int32_t GeoRef_InterpYY(TGeoRef *RefTo, TGeoRef *RefFrom,TGeoOptions *Opt,float *zout,float *zin) {

   TGeoSet *gset=NULL;
   TGeoRef *yin_gdin, *yan_gdin, *yin_gdout, *yan_gdout;
   int32_t k,ni,nj;
   int32_t yancount_yin,yincount_yin, yancount_yan,yincount_yan;
   int32_t yyin,yyout;
   float *yin2yin_zvals,*yan2yin_zvals;
   float *yin2yan_zvals,*yan2yan_zvals;

   yyin=0; yyout=0;

   if (!Opt) Opt=&RefTo->Options;
   if (!Opt) Opt=&GeoRef_Options;

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
   if (!yyin && yyout) {
      if (GeoRef_Interp(yin_gdout,RefFrom,Opt,zout,zin) && GeoRef_Interp(yan_gdout,RefFrom,Opt,&zout[ni*nj],zin)) {
         return(TRUE);
      } else {
         return(FALSE);
      }
   }

   // check if one input sub grid is identical to dest sub grid or dest single grid
   if (yin_gdin == RefTo) {
      return(GeoRef_Interp(RefTo,yin_gdin,Opt,zout,zin));
   }
   if (yan_gdin == RefTo) {
      return(GeoRef_Interp(RefTo,yan_gdin,Opt,zout,&zin[(yin_gdin->NX)*(yin_gdin->NY)]));
   }

   /* User specifies to use 1 subgrid for interpolation ezsetopt(USE_1SUBGRID) */
   /* User must specify the sub grid value in ezsetival(SUBGRIDID) */
   /* This is only appropriate if the destination grid is non yin-yang grid */

   if (RefFrom->Sub>=0) { // User specifies to use 1 grid only
      // Output is a Yin-Yang grid
      if (yyout) {
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
         return(GeoRef_Interp(yin_gdout,yin_gdin,Opt,zout,zin));
      }

      // Use yang input grid
      if (RefFrom->Sub==1) {
         return(GeoRef_Interp(yin_gdout,yan_gdin,Opt,zout,&zin[(yin_gdin->NX)*(yin_gdin->NY)]));
      }
   }

   // To use both Yin and Yang grids in Yin-yang input grid
   // Masquer les grilles YY input pour enlever overlap et calculer les X,Y
   gset=GeoRef_SetGet(RefTo,RefFrom,Opt);
   GeoRef_SetCalcYYXY(gset);

   if (!yyout) {
      // Interp yinyang to one grid

      yincount_yin = gset->yincount_yin;
      yancount_yin = gset->yancount_yin;
      yin2yin_zvals = (float *) malloc(yincount_yin*sizeof(float));
      yan2yin_zvals = (float *) malloc(yancount_yin*sizeof(float));

      GeoRef_XYVal(yin_gdin,Opt,yin2yin_zvals,zin,gset->yin2yin_x,gset->yin2yin_y,gset->yincount_yin);
      GeoRef_XYVal(yan_gdin,Opt,yan2yin_zvals,&zin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yin_x,gset->yan2yin_y,gset->yancount_yin);

      yincount_yin=0;
      yancount_yin=0;
      for(k=0; k<ni*nj; k++) {
         if (gset->yin_maskout[k] == 1.0) {
            zout[k]=yan2yin_zvals[yancount_yin];
            yancount_yin++;
         } else {
            zout[k]=yin2yin_zvals[yincount_yin];
            yincount_yin++;
         }
      }
      free(yin2yin_zvals);
      free(yan2yin_zvals);
      return(TRUE);
   } else {
      // Interp yinyang to yinyang

      yincount_yin = gset->yincount_yin;
      yancount_yin = gset->yancount_yin;
      yincount_yan = gset->yincount_yan;
      yancount_yan = gset->yancount_yan;
      yin2yin_zvals = (float *) malloc(yincount_yin*sizeof(float));
      yan2yin_zvals = (float *) malloc(yancount_yin*sizeof(float));
      yin2yan_zvals = (float *) malloc(yincount_yan*sizeof(float));
      yan2yan_zvals = (float *) malloc(yancount_yan*sizeof(float));

      GeoRef_XYVal(yin_gdin,Opt,yin2yin_zvals,zin,gset->yin2yin_x,gset->yin2yin_y,gset->yincount_yin);
      GeoRef_XYVal(yan_gdin,Opt,yan2yin_zvals,&zin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yin_x,gset->yan2yin_y,gset->yancount_yin);
      GeoRef_XYVal(yin_gdin,Opt,yin2yan_zvals,zin,gset->yin2yan_x,gset->yin2yan_y,gset->yincount_yan);
      GeoRef_XYVal(yan_gdin,Opt,yan2yan_zvals,&zin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yan_x,gset->yan2yan_y,gset->yancount_yan);

      // Interp input YY grid to Yin grid
      yincount_yin=0; yancount_yin=0;
      for(k=0; k<ni*nj; k++) {
         if (gset->yin_maskout[k] == 1.0) {
            zout[k]=yan2yin_zvals[yancount_yin];
            yancount_yin++;
         } else {
            zout[k]=yin2yin_zvals[yincount_yin];
            yincount_yin++;
         }
      }

      // Interp input YY grid to Yang grid
      yincount_yan=0; yancount_yan=0;
      for(k=0; k<ni*nj; k++) {
         if (gset->yan_maskout[k] == 1.0) {
            zout[k+(ni*nj)]=yan2yan_zvals[yancount_yan];
            yancount_yan++;
         } else {
            zout[k+(ni*nj)]=yin2yan_zvals[yincount_yan];
            yincount_yan++;
         }
      }

      free(yin2yin_zvals);
      free(yan2yin_zvals);
      free(yin2yan_zvals);
      free(yan2yan_zvals);
   }

   return(TRUE);
}
