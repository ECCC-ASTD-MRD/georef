#include <App.h>
#include "GeoRef.h"

int32_t GeoRef_XYInterp(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout,float *zin,double *X,double *Y,int32_t npts) {

   TGeoSet *gset;
   int32_t ier;
   float *lzin, *lxzin;

   lzin  = NULL;
   lxzin = NULL;
   if (!Opt) Opt=&GeoRef_Options;

   if (RefFrom->Type&GRID_YINVERT) {
      lzin = (float *) malloc(RefFrom->NX*RefFrom->NY*sizeof(float));
      memcpy(lzin, zin, RefFrom->NX*RefFrom->NY*sizeof(float));
      f77name(permut)(lzin, &RefFrom->NX, &RefFrom->NY);
   } else {
      lzin = zin;
   }

   if (RefFrom->Type&GRID_EXPAND) {
      lxzin = (float *) malloc(2*RefFrom->NX*RefFrom->NY*sizeof(float));
      GeoRef_GridGetExpanded(RefFrom,Opt,lxzin,lzin);
   } else  {
      lxzin = lzin;
   }

   ier = GeoRef_InterpFinally(RefTo,RefFrom,Opt,zout,lxzin,X,Y,npts,NULL);

/*   gset.gdin = gdin;
   gset.gdout = gdout;
   gset.x = x;
   gset.y = y;
   Grille[gdout].ni = npts;
   Grille[gdout].nj = 1;
   Grille[gdout].grtyp = 'L';
*/

/*   if (Opt->PolarCorrect == TRUE && !Opt->VectorMode)
     {
     ier = ez_defzones(&gset);
     ier = GeoRef_CorrectValue(set,zout,lxzin,&gset);
     }*/


   if (lzin && zin!=zin) {
      free(lzin);
   }

   if (lxzin && lxzin!=lzin && lxzin!=zin) {
      free(lxzin);
   }

   return(0);
}

int32_t GeoRef_XYVal(TGeoRef *Ref,TGeoOptions *Opt,float *zout,float *zin,double *X,double *Y,int32_t n) {

   TGeoRef *yin_gd, *yan_gd;
   int32_t j, icode, ni, nj;
   float *zoutyin, *zoutyan;
   double *tmpy;

   if (!Opt) Opt=&GeoRef_Options;

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
      
      icode = GeoRef_XYInterp(NULL,yin_gd,Opt,zoutyin,zin,X,tmpy,n);
      icode = GeoRef_XYInterp(NULL,yan_gd,Opt,zoutyan,&zin[ni*nj],X,tmpy,n);
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
      icode = GeoRef_XYInterp(NULL,Ref,Opt,zout,zin,X,Y,n);
      return(icode);
   }
}

int32_t GeoRef_XYUVVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *X,double *Y,int32_t n) {

   TGeoRef *yin_gd, *yan_gd;
   TGeoOptions opt=GeoRef_Options;
   int32_t j, icode, ni, nj;
   float *uuyin, *vvyin, *uuyan, *vvyan;
   double *tmpy;

   if (!Opt) Opt=&GeoRef_Options;

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

      icode = GeoRef_XYUVVal(yin_gd,Opt,uuyin,vvyin,uuin,vvin,X,tmpy,n);
      icode = GeoRef_XYUVVal(yan_gd,Opt,uuyan,vvyan,&uuin[ni*nj],&vvin[ni*nj],X,tmpy,n);
 
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
      opt.VectorMode = TRUE;
      opt.Symmetric = TRUE;
      GeoRef_XYInterp(NULL,Ref,&opt,uuout,uuin,X,Y,n);
      opt.Symmetric = FALSE;
      GeoRef_XYInterp(NULL,Ref,&opt,vvout,vvin,X,Y,n);
   }

   return(0);
}

int32_t GeoRef_XYWDVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *X,double *Y,int32_t n) {

   TGeoRef *yin_gd, *yan_gd;
   int32_t ier, j, icode, lni, lnj;
   double *tmplat, *tmplon, *tmpy;
   float *uuyin, *vvyin, *uuyan, *vvyan;
   float *tmpuu, *tmpvv;
  
   tmplat = (double*) malloc(2*n*sizeof(double));
   tmplon = &tmplat[n];
   tmpuu = (float *) malloc(2*n*sizeof(float));
   tmpvv = &tmpuu[n];
  
   if (!Opt) Opt=&GeoRef_Options;

   if (Ref->NbSub > 0) {
      yin_gd=Ref->Subs[0];
      yan_gd=Ref->Subs[1];
      lni = yin_gd->NX;
      lnj = yin_gd->NY;
      tmpy = (double *) malloc(n*sizeof(double));
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
      
      icode = GeoRef_XYUVVal(yin_gd,Opt,tmpuu,tmpvv,uuin,vvin,X,tmpy,n);
      icode = GeoRef_XY2LL(yin_gd,tmplat,tmplon,X,tmpy,n,TRUE);
      icode = GeoRef_UV2WD(yin_gd, uuyin,vvyin,tmpuu,tmpvv,tmplat,tmplon,n);

      icode = GeoRef_XYUVVal(yan_gd,Opt,tmpuu,tmpvv,&uuin[(lni*lnj)],&vvin[(lni*lnj)],X,tmpy,n);
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
      free(uuyin);
   } else {
      ier = GeoRef_XYUVVal(Ref,Opt,tmpuu,tmpvv,uuin,vvin,X,Y,n);
      ier = GeoRef_XY2LL(Ref,tmplat,tmplon,X,Y,n,TRUE);
      ier = GeoRef_UV2WD(Ref,uuout,vvout,tmpuu,tmpvv,tmplat,tmplon,n);
   }

   free(tmplat);
   free(tmpuu);

   return(0);
}