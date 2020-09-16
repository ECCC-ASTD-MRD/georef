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

int GeoRef_XYInterp(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout,float *zin,float *x,float *y,int npts) {

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
     GeoRef_ExpandGrid(RefFrom, lxzin, lzin);
   } else  {
     lxzin = lzin;
   }

   ier = GeoRef_InterpFinally(RefTo,RefFrom,zout,lxzin,x,y,npts);

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
     ier = ez_corrval(RefTo,RefFrom,zout, lxzin, &gset);
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

int GeoRef_XYVal(TGeoRef *Ref,float *zout,float *zin,float *x,float *y,int n) {

   TGeoRef *yin_gd, *yan_gd;
   int j, icode, ni, nj;
   float *zoutyin, *zoutyan;
   float *tmpy;

   if (Ref->NbSub > 0) {
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      ni = yin_gd->NX;
      nj = yin_gd->NY;
      tmpy = (float *) malloc(n*sizeof(float));
      zoutyin = (float *) malloc(n*sizeof(float));
      zoutyan = (float *) malloc(n*sizeof(float));
      for (j=0; j< n; j++) {
         if (y[j] > yin_gd->NY) {
            tmpy[j]=y[j]-yin_gd->NY;
         } else {
            tmpy[j]=y[j];
         }
      }
      icode = GeoRef_XYInterp(NULL,yin_gd,zoutyin,zin,x,tmpy,n);
      icode = GeoRef_XYInterp(NULL,yan_gd,zoutyan,&zin[ni*nj],x,tmpy,n);
      for (j=0; j < n; j++) {
         if (y[j] > yin_gd->NY) {
           zout[j]=zoutyan[j];
         } else {
           zout[j]=zoutyin[j];
         }
      }
      free(tmpy);
      free(zoutyan); free(zoutyin);
      return(icode);
   } else {
      icode = GeoRef_XYInterp(NULL,Ref,zout,zin,x,y,n);
      return(icode);
   }
}

int GeoRef_XYUVVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,float *x,float *y,int n) {

   TGeoRef *yin_gd, *yan_gd;
   int j, icode, ni, nj;
   float *uuyin, *vvyin, *uuyan, *vvyan, *tmpy;

   if (Ref->NbSub > 0) {
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      ni = yin_gd->NX;
      nj = yin_gd->NY;
      tmpy = (float *) malloc(n*sizeof(float));
      uuyin = (float *) malloc(n*sizeof(float));
      vvyin = (float *) malloc(n*sizeof(float));
      uuyan = (float *) malloc(n*sizeof(float));
      vvyan = (float *) malloc(n*sizeof(float));
      for (j=0; j< n; j++) {
         if (y[j] > yin_gd->NY) {
            tmpy[j]=y[j]-yin_gd->NY;
         } else {
            tmpy[j]=y[j];
         }
      }
      icode = GeoRef_XYUVValN(yin_gd,uuyin,vvyin,uuin,vvin,x,tmpy,n);
      icode = GeoRef_XYUVValN(yan_gd,uuyan,vvyan,&uuin[ni*nj],&vvin[ni*nj],x,tmpy,n);
 
      for (j=0; j < n; j++) {
         if (y[j] > yin_gd->NY) {
            uuout[j]=uuyan[j];
            vvout[j]=vvyan[j];
         } else {
           uuout[j]=uuyin[j];
           vvout[j]=vvyin[j];
         }
      }
      free(tmpy);
      free(uuyan); free(vvyan); 
      free(uuyin); free(vvyin);
      return(icode);
   } else {
      icode = GeoRef_XYUVValN(Ref,uuout,vvout,uuin,vvin,x,y,n);
      return(icode);
   }
}

int GeoRef_XYUVValN(TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin,float *x,float *y,int n) {

   RefFrom->Options.VectorMode = TRUE;
   RefFrom->Options.Symmetric = TRUE;
   GeoRef_XYInterp(NULL,RefFrom,uuout,uuin,x,y,n);
   RefFrom->Options.Symmetric = FALSE;
   GeoRef_XYInterp(NULL,RefFrom,vvout,vvin,x,y,n);
   RefFrom->Options.Symmetric = TRUE;
   RefFrom->Options.VectorMode = FALSE;

   return(0);
}

int GeoRef_XYWDVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,float *x,float *y,int n) {

   TGeoRef *yin_gd, *yan_gd;
   int ier, j, icode, lni, lnj;
   float *tmplat, *tmplon, *tmpy;
   float *uuyin, *vvyin, *uuyan, *vvyan;
   float *tmpuu, *tmpvv;
  
   tmplat = (float *) malloc(n * sizeof(float));
   tmplon = (float *) malloc(n * sizeof(float));
   tmpuu = (float *) malloc(n * sizeof(float));
   tmpvv = (float *) malloc(n * sizeof(float));
  
   if (Ref->NbSub > 0) {
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      lni = yin_gd->NX;
      lnj = yin_gd->NY;
      tmpy = (float *) malloc(n*sizeof(float));
      uuyin = (float *) malloc(n*sizeof(float));
      vvyin = (float *) malloc(n*sizeof(float));
      uuyan = (float *) malloc(n*sizeof(float));
      vvyan = (float *) malloc(n*sizeof(float));
      for (j=0; j< n; j++) {
         if (y[j] > yin_gd->NY) {
            tmpy[j]=y[j]-yin_gd->NY;
         } else {
            tmpy[j]=y[j];
         }
      }
      icode = GeoRef_XYUVValN(yin_gd, tmpuu, tmpvv, uuin, vvin, x, tmpy, n);
      icode = GeoRef_XY2LLN (yin_gd, tmplat, tmplon, x, tmpy, n,FALSE);
      icode = GeoRef_UV2WD (yin_gd, uuyin,vvyin,tmpuu,tmpvv,tmplat,tmplon,n);

      icode = GeoRef_XYUVValN(yan_gd, tmpuu, tmpvv, &uuin[(lni*lnj)], &vvin[(lni*lnj)], x, tmpy, n);
      icode = GeoRef_XY2LLN (yan_gd, tmplat, tmplon, x, tmpy, n,FALSE);
      icode = GeoRef_UV2WD (yan_gd, uuyan,vvyan,tmpuu,tmpvv,tmplat,tmplon,n);

      for (j=0; j< n; j++) {
         if (y[j] > yin_gd->NY) {
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
      ier = GeoRef_XYUVVal(Ref,tmpuu,tmpvv,uuin,vvin,x,y,n);
      ier = GeoRef_XY2LLN(Ref,tmplat,tmplon,x,y,n,FALSE);
      ier = GeoRef_UV2WD(Ref,uuout,vvout,tmpuu,tmpvv,tmplat,tmplon,n);
   }

   free(tmplat);
   free(tmplon);
   free(tmpuu);
   free(tmpvv);

   return(0);
}