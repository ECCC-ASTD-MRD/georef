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

#include "GeoRef.h"

int GeoRef_XYInterp(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout,float *zin,double *X,double *Y,int npts) {

   TGridSet *gset;
   int ier;
   float *lzin, *lxzin;

   lzin  = NULL;
   lxzin = NULL;

   if (RefFrom->Type&GRID_YINVERT) {
      lzin = (float *) malloc(RefFrom->NX*RefFrom->NY*sizeof(float));
      memcpy(lzin, zin, RefFrom->NX*RefFrom->NY*sizeof(float));
      f77name(permut)(lzin, &RefFrom->NX, &RefFrom->NY);
   } else {
     lzin = zin;
   }

   if (RefFrom->Type&GRID_EXPAND) {
     lxzin = (float *) malloc(2*RefFrom->NX*RefFrom->NY*sizeof(float));
     GeoRef_GridGetExpanded(RefFrom,lxzin,lzin);
   } else  {
     lxzin = lzin;
   }

   ier = GeoRef_InterpFinally(RefTo,RefFrom,zout,lxzin,X,Y,npts);

/*   gset.gdin = gdin;
   gset.gdout = gdout;
   gset.x = x;
   gset.y = y;
   Grille[gdout].ni = npts;
   Grille[gdout].nj = 1;
   Grille[gdout].GRTYP = 'L';
*/

/*   if (Opt->PolarCorrect == TRUE && !Opt->VectorMode)
     {
     ier = ez_defzones(&gset);
     ier = GeoRef_CorrectValue(RefTo,RefFrom,zout, lxzin, &gset);
     }*/


/*   if (lzin != zin && lzin != NULL)
      {
      free(lzin);
      }

   if (lxzin != lzin && lxzin != zin && lxzin != NULL)
      {
      free(lxzin);
      }*/

   return(0);
}

int GeoRef_XYVal(TGeoRef *Ref,float *zout,float *zin,double *X,double *Y,int n) {

   TGeoRef *yin_gd, *yan_gd;
   int j, icode, ni, nj;
   float *zoutyin, *zoutyan;
   double *tmpy;

   if (Ref->NbSub > 0) {
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      ni = yin_gd->NX;
      nj = yin_gd->NY;
      tmpy = (double*) malloc(n*sizeof(double));
      zoutyin = (float *) malloc(2*n*sizeof(float));
      zoutyan = &zoutyin[n];
      for (j=0; j< n; j++) {
         if (Y[j] > yin_gd->NY) {
            tmpy[j]=Y[j]-yin_gd->NY;
         } else {
            tmpy[j]=Y[j];
         }
      }

      // TODO: WTF will I do to avoid this or global options
      memcpy(&yin_gd->Options,&Ref->Options,sizeof(TGeoOptions));
      memcpy(&yan_gd->Options,&Ref->Options,sizeof(TGeoOptions));
      
      icode = GeoRef_XYInterp(NULL,yin_gd,zoutyin,zin,X,tmpy,n);
      icode = GeoRef_XYInterp(NULL,yan_gd,zoutyan,&zin[ni*nj],X,tmpy,n);
      for (j=0; j < n; j++) {
         if (Y[j] > yin_gd->NY) {
           zout[j]=zoutyan[j];
         } else {
           zout[j]=zoutyin[j];
         }
      }
      free(tmpy);
      free(zoutyan);
      return(icode);
   } else {
      icode = GeoRef_XYInterp(NULL,Ref,zout,zin,X,Y,n);
      return(icode);
   }
}

int GeoRef_XYUVVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,double *X,double *Y,int n) {

   TGeoRef *yin_gd, *yan_gd;
   int j, icode, ni, nj;
   float *uuyin, *vvyin, *uuyan, *vvyan;
   double *tmpy;

   if (Ref->NbSub > 0) {
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      ni = yin_gd->NX;
      nj = yin_gd->NY;
      tmpy = (double*) malloc(n*sizeof(double));
      uuyin = (float *) malloc(4*n*sizeof(float));
      vvyin = &uuyin[n];
      uuyan = &uuyin[n*2];
      vvyan = &uuyin[n*3];

      for (j=0; j< n; j++) {
         if (Y[j] > yin_gd->NY) {
            tmpy[j]=Y[j]-yin_gd->NY;
         } else {
            tmpy[j]=Y[j];
         }
      }

      // TODO: WTF will I do to avoid this or global options
      memcpy(&yin_gd->Options,&Ref->Options,sizeof(TGeoOptions));
      memcpy(&yan_gd->Options,&Ref->Options,sizeof(TGeoOptions));

      icode = GeoRef_XYUVVal(yin_gd,uuyin,vvyin,uuin,vvin,X,tmpy,n);
      icode = GeoRef_XYUVVal(yan_gd,uuyan,vvyan,&uuin[ni*nj],&vvin[ni*nj],X,tmpy,n);
 
      for (j=0; j < n; j++) {
         if (Y[j] > yin_gd->NY) {
            uuout[j]=uuyan[j];
            vvout[j]=vvyan[j];
         } else {
            uuout[j]=uuyin[j];
            vvout[j]=vvyin[j];
         }
      }
      free(tmpy); 
      free(uuyin);
      return(icode);
   } else {
      Ref->Options.VectorMode = TRUE;
      Ref->Options.Symmetric = TRUE;
      GeoRef_XYInterp(NULL,Ref,uuout,uuin,X,Y,n);
      Ref->Options.Symmetric = FALSE;
      GeoRef_XYInterp(NULL,Ref,vvout,vvin,X,Y,n);
      Ref->Options.Symmetric = TRUE;
      Ref->Options.VectorMode = FALSE;
   }

   return(0);
}

int GeoRef_XYWDVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,double *X,double *Y,int n) {

   TGeoRef *yin_gd, *yan_gd;
   int ier, j, icode, lni, lnj;
   double *tmplat, *tmplon, *tmpy;
   float *uuyin, *vvyin, *uuyan, *vvyan;
   float *tmpuu, *tmpvv;
  
   tmplat = (double*) malloc(n*sizeof(double));
   tmplon = (double*) malloc(n*sizeof(double));
   tmpuu = (float *) malloc(n * sizeof(float));
   tmpvv = (float *) malloc(n * sizeof(float));
  
   if (Ref->NbSub > 0) {
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      lni = yin_gd->NX;
      lnj = yin_gd->NY;
      tmpy = (double *) malloc(n*sizeof(double));
      uuyin = (float *) malloc(n*sizeof(float));
      vvyin = (float *) malloc(n*sizeof(float));
      uuyan = (float *) malloc(n*sizeof(float));
      vvyan = (float *) malloc(n*sizeof(float));
      for (j=0; j< n; j++) {
         if (Y[j] > yin_gd->NY) {
            tmpy[j]=Y[j]-yin_gd->NY;
         } else {
            tmpy[j]=Y[j];
         }
      }
      
      // TODO: WTF will I do to avoid this or global options
      memcpy(&yin_gd->Options,&Ref->Options,sizeof(TGeoOptions));
      memcpy(&yan_gd->Options,&Ref->Options,sizeof(TGeoOptions));

      icode = GeoRef_XYUVVal(yin_gd, tmpuu, tmpvv, uuin, vvin, X, tmpy, n);
      icode = GeoRef_XY2LL(yin_gd,tmplat,tmplon,X,tmpy,n,TRUE);
      icode = GeoRef_UV2WD(yin_gd, uuyin,vvyin,tmpuu,tmpvv,tmplat,tmplon,n);

      icode = GeoRef_XYUVVal(yan_gd, tmpuu, tmpvv, &uuin[(lni*lnj)], &vvin[(lni*lnj)], X, tmpy, n);
      icode = GeoRef_XY2LL(yan_gd,tmplat,tmplon,X,tmpy,n,TRUE);
      icode = GeoRef_UV2WD(yan_gd, uuyan,vvyan,tmpuu,tmpvv,tmplat,tmplon,n);

      for (j=0; j< n; j++) {
         if (Y[j] > yin_gd->NY) {
            uuout[j]=uuyan[j];
            vvout[j]=vvyan[j];
         } else {
            uuout[j]=uuyin[j];
            vvout[j]=vvyin[j];
         }
      }
      free(uuyin); free(vvyin);
      free(uuyan); free(vvyan);
      free(tmpy);
   } else {
      ier = GeoRef_XYUVVal(Ref,tmpuu,tmpvv,uuin,vvin,X,Y,n);
      ier = GeoRef_XY2LL(Ref,tmplat,tmplon,X,Y,n,TRUE);
      ier = GeoRef_UV2WD(Ref,uuout,vvout,tmpuu,tmpvv,tmplat,tmplon,n);
   }

   free(tmplat);
   free(tmplon);
   free(tmpuu);
   free(tmpvv);

   return(0);
}