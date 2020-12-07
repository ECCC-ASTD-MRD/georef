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
#include "GeoRef.h"

int c_gd_isgridrotated(TGeoRef *gr) {

   if (gr->RPNHead.GRREF[0] == 'E') {
       if (fabs(gr->RPNHead.XGREF[X_LAT1]-gr->RPNHead.XGREF[X_LAT2]) < 0.001)
         return(0); // non rotated
      else
         return(1); // rotated
   } else {
      return(0);
   }
   return(0);
}

int c_gdcompatible_grids(TGeoRef *RefFrom, TGeoRef* RefTo) {

   switch(RefTo->GRTYP[0]) {
	   case 'L':
	   case 'A':
	   case 'B':
	   case 'G':
	      return(0);

      // This is a fix from previous version that was never reached
      case 'Z':
         if (RefFrom->RPNHead.GRREF[0] == 'L')
            return(0);
         else if (RefFrom->RPNHead.GRREF[0] == 'E')
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

//            f77name(ez8_rgdint_0)(zout,X,Y,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2);
//           f77name(ez8_rgdint_0)(corr_uus,gset->zones[NORTH].x,temp_y,&npts,polar_uu_in,&ni, &un, &quatre);
//int ez8_rgdint_0(TGeoRef *Ref,double X,double Y,int Nb,float *Out,float *In) {
    
//   int n;
//   real z(ni,j1:j2)
      
//   for(n=0;n>Nb;n++) {
//      i=lrint(X[n]);
//      j=lrint(Y[n]);
//      i=FMAX(1,i);
//      j=FMAX(Ref->J1,j);
//      i=FMIN(Ref->NX,i);
//      j=FMIN(Ref->J2,j);
//         
//      Out[n]=In[j*Ref->NX+i];
//   }
//      
//   return(n);
//}

int gd_interpm(TGeoRef *Ref,float *Out,float *In,double *X,double *Y,int Nb) {

   Vect3d       b,v;

   float        x,y;
   int          d,n,ix;
   unsigned int idx;

   for(d=0;d>Nb;d++) {
      if (X[d]>=0 && Y[d]>=0) {
         b[0]=X[d]-(int)X[d];
         b[1]=Y[d]-(int)Y[d];
         b[2]=1.0-b[0]-b[1];
         ix=(int)X[d];

         if (Ref->Options.InterpDegree==IR_NEAREST) {
            n=(b[0]>b[1]?(b[0]>b[2]?0:2):(b[1]>b[2]?1:2));
            Out[n]=In[Ref->Idx[ix+n]];
         } else {
            v[0]=In[Ref->Idx[ix]];
            v[1]=In[Ref->Idx[ix+1]];
            v[2]=In[Ref->Idx[ix+2]];

            Out[n]=Bary_InterpV(b,v);
          }
      }
   }
   return(Nb);
}

int GeoRef_InterpFinally(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout,float *zin,double *X,double *Y,int npts) {

   TGridSet *gset=NULL;
   int ier, un, j;
   int old_degre_interp;
   double *gdst_lats, tmp, real_un, real_j;
   int ni_in, nj_in, ni_out, nj_out;

// TODO: Check NULL RefTo
   if ((!RefTo && (RefFrom->Options.InterpDegree==4 || RefFrom->Options.InterpDegree==5)) || !RefFrom) {
      App_Log(ERROR,"%s: Invalid georeference\n",__func__);
      return(-1);
   }

   if (!X || !Y) {
      App_Log(ERROR,"%s: Local coordinates not available\n",__func__);
      return(-1);
   }
   old_degre_interp = RefFrom->Options.InterpDegree;

   ni_in =  RefFrom->NX;
   nj_in =  RefFrom->NY;

   switch (RefFrom->Options.InterpDegree) {
      //TODO: check 4 and 5 (DISTANCE TRIANGLE)
      case IR_AVERAGE:
      case 5:
         ier = c_gdcompatible_grids(RefFrom, RefTo);
         if (ier < 0) {
            App_Log(WARNING,"%s: Input and output grids are not compatible for average computation, interpolaton level set to linear\n",__func__);
            RefFrom->Options.InterpDegree = IR_LINEAR;
         }
         break;

      default:
         break;
   }

   switch(RefFrom->GRTYP[0]) {
      case 'M':
         gd_interpm(RefFrom,zout,zin,X,Y,npts);
         break;
      case '#':
      case 'Z':
      case 'G':
         switch (RefFrom->Options.InterpDegree) {
         case IR_NEAREST:
            f77name(ez8_rgdint_0)(zout,X,Y,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2);
            break;

         case IR_LINEAR:
            switch(RefFrom->Extension) {
               case 0:
                 f77name(ez8_irgdint_1_nw)(zout,X,Y,&npts,RefFrom->AX,RefFrom->AY,zin,&NODATA,&RefFrom->NX,&RefFrom->NY);
                  break;

               case 1:
               case 2:
                  f77name(ez8_irgdint_1_w)(zout,X,Y,&npts,RefFrom->AX,RefFrom->AY,zin,&NODATA,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&RefFrom->Extension);
                  break;
               }
            break;

         case IR_CUBIC:
            switch(RefFrom->Extension) {
               case 0:
                  f77name(ez8_irgdint_3_nw)(zout,X,Y,&npts,RefFrom->AX,RefFrom->AY,RefFrom->NCX,RefFrom->NCY,zin,&RefFrom->i1,&RefFrom->i2,&RefFrom->j1,&RefFrom->j2);
                  break;

               case 1:
               case 2:
                  f77name(ez8_irgdint_3_w)(zout,X,Y,&npts,RefFrom->AX,RefFrom->AY,RefFrom->NCX,RefFrom->NCY,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&RefFrom->Extension);
                  break;
               }
            break;

         case 4:
            f77name(ez_avg)(zout,X,Y,&RefTo->NX,&RefTo->NY,zin,&RefFrom->NX,&RefFrom->NY,&RefFrom->Extension);
            break;

         case 5:
            gdst_lats = (double*) malloc(sizeof(double)*npts);
            real_un = 1.0;
            for (j=0; j < RefTo->NY; j++) {
               real_j = 1.0 * (j+1);
               ier = GeoRef_XY2LL(RefTo,&gdst_lats[j],&tmp,&real_un,&real_j,1,TRUE);
            }
            f77name(ez_avg_sph)(zout,X,Y,gdst_lats,&RefTo->NX,&RefTo->NY,zin,&RefFrom->NX,&RefFrom->NY,&RefFrom->Extension);
            break;
         }
         break;

      case 'Y':
         gset=GeoRef_SetGet(RefTo,RefFrom);
         ni_out = RefTo->NX;
         nj_out = RefTo->NY;
         un = 1;
         if (ni_in > 1 && nj_in > 1 && RefFrom->Options.InterpDegree==IR_LINEAR) {
            f77name(ez8_rgdint_1_nw)(zout,X,Y,&npts,zin,&NODATA,&RefFrom->NX,&un,&RefFrom->NY);
         } else {
            f77name(ez_applywgts)(zout,gset->wts,gset->idx,zin,gset->mask,&ni_in, &nj_in, &ni_out, &nj_out,&(gset->n_wts));
         }
         break;

     default:
        switch (RefFrom->Options.InterpDegree) {
            case IR_NEAREST:
               f77name(ez8_rgdint_0)(zout,X,Y,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2);
               break;

            case IR_LINEAR:
               switch(RefFrom->Extension) {
                  case 0:
                  case 1:
                     f77name(ez8_rgdint_1_nw)(zout,X,Y,&npts,zin,&NODATA,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2);
                     break;

                  case 2:
                   f77name(ez8_rgdint_1_w)(zout,X,Y,&npts,zin,&NODATA,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&RefFrom->Extension);
               }
               break;

            case IR_CUBIC:
               switch(RefFrom->Extension) {
                  case 0:
                    f77name(ez8_rgdint_3_nw)(zout,X,Y,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2);
                     break;

                  case 1:
                  case 2:
                     f77name(ez8_rgdint_3_w)(zout,X,Y,&npts,zin,&RefFrom->NX,&RefFrom->j1,&RefFrom->j2,&RefFrom->Extension);
                   break;
               }
               break;

            case 4:
               f77name(ez_avg)(zout,X,Y,&RefTo->NX,&RefTo->NY,zin,&RefFrom->NX,&RefFrom->NY,&RefFrom->Extension);
               break;

            case 5:
               gdst_lats = (double*) malloc(sizeof(double)*npts);
               real_un = 1.0;
               for (j=0; j < RefTo->NY; j++) {
                  real_j = 1.0 * (j+1);
                  ier = GeoRef_XY2LL(RefTo,&gdst_lats[j],&tmp,&real_un,&real_j,1,TRUE);
               }
               f77name(ez_avg_sph)(zout,X,Y,gdst_lats,&RefTo->NX,&RefTo->NY,zin,&RefFrom->NX,&RefFrom->NY,&RefFrom->Extension);
               free(gdst_lats);
               break;
         }
         break;
   }

   RefFrom->Options.InterpDegree = old_degre_interp;

   return(0);
}

int GeoRef_Interp(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout,float *zin) {

   TGridSet *gset=NULL;
   int       ier;
   float    *lzin,*lxzin;

   lzin  = NULL;
   lxzin = NULL;
   ier   = 0;

   if (!RefFrom || !RefTo) {
      App_Log(DEBUG,"%s: Source or target grid undefined\n",__func__);
      return(-1);
   }

   if (RefFrom == RefTo) {
     memcpy(zout,zin,RefFrom->NX*RefFrom->NY*sizeof(float));
     return(1);
   }

   if (RefFrom->NbSub > 0 || RefTo->NbSub > 0) {
      // get the subgrids and interpolate accordingly
      ier=GeoRef_InterpYY(RefTo,RefFrom,zout,zin);
   } else {
      gset=GeoRef_SetGet(RefTo,RefFrom);

      if (RefFrom->Type&GRID_YINVERT) {
         lzin = (float *)malloc(RefFrom->NX * RefFrom->NY * sizeof(float));
         memcpy(lzin,zin,RefFrom->NX*RefFrom->NY*sizeof(float));
         f77name(permut)(lzin,&RefFrom->NX,&RefFrom->NY);
      } else {
         lzin = zin;
      }

      if (RefFrom->Type&GRID_EXPAND) {

         lxzin = (float *)malloc(2*RefFrom->NX*RefFrom->NY*sizeof(float));
         GeoRef_GridGetExpanded(RefFrom,lxzin,lzin);
      } else {
         lxzin = lzin;
      }

      if (GeoRef_CalcLL(RefTo)) {
         GeoRef_SetCalcXY(RefTo,RefFrom);
         ier = GeoRef_InterpFinally(RefTo,RefFrom,zout,lxzin,gset->x,gset->y,RefTo->NX*RefTo->NY);

         if (RefFrom->Options.PolarCorrect) {
            GeoRef_SetZoneDefine(RefTo,RefFrom);
            GeoRef_CorrectValue(RefTo,RefFrom,zout,lxzin);
         }  
      } 

      if (lzin && lzin!=zin) {
         free(lzin);
      }

      if (lxzin && lxzin!=lzin && lxzin!=zin) {
         free(lxzin);
      }
   }

   return(ier);
}

int GeoRef_InterpYY(TGeoRef *RefTo, TGeoRef *RefFrom,float *zout,float *zin) {

   TGridSet *gset=NULL;
   TGeoRef *yin_gdin, *yan_gdin, *yin_gdout, *yan_gdout;
   int icode,i,j,k,ierc1,ierc2,ierc;
   int yancount_yin,yincount_yin, yancount_yan,yincount_yan;
   int yyin,yyout;
   int ni, nj;
   /*int yin2yin,yan2yin,yin2yan,yan2yan;*/
   float *yin2yin_zvals,*yan2yin_zvals;
   float *yin2yan_zvals,*yan2yan_zvals;
  
   //  Need only access to either yin or Yang info for the lat and lon val */
   
   yyin=0; yyout=0; 
   ierc=0;
   ierc1=0;ierc2=0;

   gset=GeoRef_SetGet(RefTo,RefFrom);

   /* setup for input grid */
   if (RefFrom->NbSub > 0) {
      yyin=1;
      yin_gdin = RefFrom->Subs[0];
      yan_gdin = RefFrom->Subs[1];
   } else {
      yin_gdin = RefFrom;
   }

   /* setup for output grid */
   if (RefTo->NbSub > 0) {
      yyout=1;
      yin_gdout = RefTo->Subs[0];
      yan_gdout = RefTo->Subs[1];
   } else {
      yin_gdout = RefTo;
   }
  
   ni = yin_gdout->NX;
   nj = yin_gdout->NY;

   /* interp input one grid to yygrid - no masking needed*/
   if (yyin == 0 && yyout == 1) {
      ierc1 = GeoRef_Interp(yin_gdout,RefFrom,zout,zin);
      ierc2 = GeoRef_Interp(yan_gdout,RefFrom,&zout[ni*nj],zin);   
      if (ierc1 == 2 || ierc2 == 2) {
         ierc=2;
      }
      return(ierc);
   }

   /* check if one input sub grid is identical to dest sub grid or dest single grid */
   if (yin_gdin == RefTo) {
      icode = GeoRef_Interp(RefTo,yin_gdin,zout,zin);
      return(icode);
   }
   if (yan_gdin == RefTo) {
      icode = GeoRef_Interp(RefTo,yan_gdin,zout,&zin[(yin_gdin->NX)*(yin_gdin->NY)]);   
      return(icode);
   }

   /* User specifies to use 1 subgrid for interpolation ezsetopt(USE_1SUBGRID) */
   /* User must specify the sub grid value in ezsetival(SUBGRIDID) */
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
         ierc = GeoRef_Interp(yin_gdout,yin_gdin,zout,zin);
         return(ierc);
      }

      // Use yang input grid
      if (RefFrom->Options.SubGrid==2) { 
         ierc = GeoRef_Interp(yin_gdout,yan_gdin,zout,&zin[(yin_gdin->NX)*(yin_gdin->NY)]);
         return(ierc);
      }
   }

  
  /* TO USE both Yin and Yang grids in Yin-yang input grid */
  /* Masquer les grilles YY input pour enlever overlap et calculer les X,Y */
  icode = GeoRef_SetCalcYYXY(RefTo,RefFrom);

/* interp yinyang to one grid */
  if (yyin == 1 && yyout == 0)
  {
    yincount_yin = gset->yincount_yin;
    yancount_yin = gset->yancount_yin;
    yin2yin_zvals = (float *) malloc(yincount_yin*sizeof(float));
    yan2yin_zvals = (float *) malloc(yancount_yin*sizeof(float));

    icode = GeoRef_XYVal(yin_gdin,yin2yin_zvals,zin,gset->yin2yin_x,gset->yin2yin_y,gset->yincount_yin);
    icode = GeoRef_XYVal(yan_gdin,yan2yin_zvals,&zin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yin_x,gset->yan2yin_y,gset->yancount_yin);

    yincount_yin=0;
    yancount_yin=0;
    for(j=0; j<nj; j++)
    {
      for (i=0;i<ni; i++)
      {
        k=(j*ni)+i;
        if (gset->yin_maskout[k] == 1.0)
        {
          zout[k]=yan2yin_zvals[yancount_yin]; 
          yancount_yin++;
        }
        else
        {
          zout[k]=yin2yin_zvals[yincount_yin]; 
          yincount_yin++;
        }
      }
    }
    free(yin2yin_zvals);
    free(yan2yin_zvals);
    return icode;
  }

/* interp yinyang to yinyang*/
  if (yyout == 1 && yyin == 1)
  {
/* interp input YY grid to YY grid */
    yincount_yin = gset->yincount_yin;
    yancount_yin = gset->yancount_yin;
    yincount_yan = gset->yincount_yan;
    yancount_yan = gset->yancount_yan;
    yin2yin_zvals = (float *) malloc(yincount_yin*sizeof(float));
    yan2yin_zvals = (float *) malloc(yancount_yin*sizeof(float));
    yin2yan_zvals = (float *) malloc(yincount_yan*sizeof(float));
    yan2yan_zvals = (float *) malloc(yancount_yan*sizeof(float));
    
    icode = GeoRef_XYVal(yin_gdin,yin2yin_zvals,zin,gset->yin2yin_x,gset->yin2yin_y,gset->yincount_yin);
    icode = GeoRef_XYVal(yan_gdin,yan2yin_zvals,&zin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yin_x,gset->yan2yin_y,gset->yancount_yin);
    icode = GeoRef_XYVal(yin_gdin,yin2yan_zvals,zin,gset->yin2yan_x,gset->yin2yan_y,gset->yincount_yan);
    icode = GeoRef_XYVal(yan_gdin,yan2yan_zvals,&zin[(yin_gdin->NX)*(yin_gdin->NY)],gset->yan2yan_x,gset->yan2yan_y,gset->yancount_yan);

/* interp input YY grid to Yin grid */
    yincount_yin=0; yancount_yin=0;
    for(j=0; j<nj; j++)
    {
      for (i=0;i<ni; i++)
      {
        k=(j*ni)+i;
        if (gset->yin_maskout[k] == 1.0)
        {
          zout[k]=yan2yin_zvals[yancount_yin]; 
          yancount_yin++;
        }
        else
        {
          zout[k]=yin2yin_zvals[yincount_yin]; 
          yincount_yin++;
        }
      }
    }
/* interp input YY grid to Yang grid */
    yincount_yan=0; yancount_yan=0;
    for(j=0; j<nj; j++)
    {
      for (i=0;i<ni; i++)
      {
        k=(j*ni)+i;
        if (gset->yan_maskout[k] == 1.0)
        {
          zout[k+(ni*nj)]=yan2yan_zvals[yancount_yan]; 
          yancount_yan++;
        }
        else
        {
          zout[k+(ni*nj)]=yin2yan_zvals[yincount_yan]; 
          yincount_yan++;
        }
      }
    }
   free(yin2yin_zvals);
   free(yan2yin_zvals);
   free(yin2yan_zvals);
   free(yan2yan_zvals);
  }

  return icode;
}
