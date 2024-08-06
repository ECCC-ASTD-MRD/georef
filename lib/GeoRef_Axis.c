#include <math.h>
#include <App.h>
#include "GeoRef.h"

void GeoRef_GridGetExpanded(TGeoRef *Ref,float *zout,float *zin) {   

   switch (Ref->GRTYP[0]) {
      case 'A':
      case 'G':
         f77name(ez_xpngdag2)(zout,zin,&Ref->NX,&Ref->NY,&Ref->j1,&Ref->j2,&Ref->RPNHead.ig1,&Ref->Options.Symmetric);
         break;

      case 'B':
         f77name(ez_xpngdb2)(zout,zin,&Ref->NX,&Ref->NY,&Ref->j1,&Ref->j2,&Ref->RPNHead.ig1,&Ref->Options.Symmetric);
         break;

      default:
         break;
   }
}

void GeoRef_AxisCalcNewtonCoeff(TGeoRef* Ref) {

   int nni,nnj;
 
   if (Ref->GRTYP[0]!='Y' && !Ref->NCX) {

      nni = Ref->NX;
      nnj = Ref->j2 - Ref->j1 + 1;

      Ref->NCX = (float *) malloc(nni*6*sizeof(float));
      Ref->NCY = (float *) malloc(nnj*6*sizeof(float));
      f77name(ez8_nwtncof)(Ref->NCX,Ref->NCY,Ref->AX,Ref->AY,&Ref->NX,&Ref->NY,&Ref->i1,&Ref->i2,&Ref->j1,&Ref->j2,&Ref->Extension);
   }  
}

void GeoRef_AxisCalcExpandCoeff(TGeoRef* Ref) {

   Ref->i1 = 1;
   Ref->i2 = Ref->NX;
   switch (Ref->GRTYP[0]) {
      case '!':
      case 'L':
      case 'N':
      case 'S':
      case 'T':
        Ref->j1 = 1;
        Ref->j2 = Ref->NY;
        Ref->Extension = 0;
        break;

      case 'A':
      case 'G':
	      Ref->Extension = 2;
         switch (Ref->RPNHead.ig1) {
	         case GLOBAL:
	            Ref->j1 = 1;
	            Ref->j2 = Ref->NY;
	            break;

	         case NORTH:
	            Ref->j1 = -Ref->NY+1;
	            Ref->j2 =  Ref->NY;
	            break;

	         case SOUTH:
	            Ref->j1 = 1;
	            Ref->j2 =  2 * Ref->NY;
	            break;
	      }
         break;

      case 'B':
	      Ref->Extension = 1;
         switch (Ref->RPNHead.ig1) {
	         case GLOBAL:
	            Ref->j1 = 1;
	            Ref->j2 = Ref->NY;
	            break;

	         case NORTH:
	            Ref->j1 = -Ref->NY+2;
	            Ref->j2 = Ref->NY;
	            break;

	         case SOUTH:
	            Ref->j1 = 1;
	            Ref->j2 = 2 * Ref->NY - 1;
	            break;
	      }
         break;

      case 'E':
         Ref->j1 = 1;
         Ref->j2 = Ref->NY;
         break;

      case '#':
      case 'Z':
         switch (Ref->RPNHeadExt.grref[0]) {
	         case 'E':
	            Ref->j1 = 1;
	            Ref->j2 = Ref->NY;
	            if ((Ref->AX[Ref->NX-1]-Ref->AX[0]) < 359.0) {
	               Ref->Extension = 0;
	            } else {
	               Ref->Extension = 1;
	            }
               break;

	         default:
	            Ref->j1 = 1;
	            Ref->j2 = Ref->NY;
	            Ref->Extension = 0;
	            break;
         }
	      break;

      default:
	      Ref->j1 = 1;
	      Ref->j2 = Ref->NY;
	      Ref->Extension = 0;
	      break;
   }
}

int GeoRef_AxisGet(TGeoRef *Ref,double *AX,double *AY) {

   int nix,njy;
   
   switch(Ref->GRTYP[0]) {
      case 'Y':
        nix = Ref->NX * Ref->NY;
        njy = nix;
        break;

      default:
        nix = Ref->NX;
        njy = Ref->NY;
        break;
   }
   
   if (Ref->AX) {
      memcpy(AX,Ref->AX, nix*sizeof(double));
      memcpy(AY,Ref->AY, njy*sizeof(double));
   } else {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Descriptor not found\n",__func__);
      return(-1);
   }
   return(0);
}

void GeoRef_AxisDefine(TGeoRef* Ref,double *AX,double *AY) {
  
  int     i,j;
  double *temp, dlon;
  int     zero,deuxnj;

  switch (Ref->GRTYP[0]) {
    case '#':
    case 'Z':
      f77name(cigaxg)(Ref->RPNHeadExt.grref,&Ref->RPNHeadExt.xgref1, &Ref->RPNHeadExt.xgref2, &Ref->RPNHeadExt.xgref3, &Ref->RPNHeadExt.xgref4,
		      &Ref->RPNHeadExt.igref1, &Ref->RPNHeadExt.igref2, &Ref->RPNHeadExt.igref3, &Ref->RPNHeadExt.igref4,1);

      Ref->AX = (double*) malloc(Ref->NX*sizeof(double));
      Ref->AY = (double*) malloc(Ref->NY*sizeof(double));

      memcpy(Ref->AX,AX,Ref->NX*sizeof(double));
      memcpy(Ref->AY,AY,Ref->NY*sizeof(double));
      GeoRef_AxisCalcExpandCoeff(Ref);
      GeoRef_AxisCalcNewtonCoeff(Ref);
      break;

    case 'Y':
      Ref->AX = (double*) malloc(Ref->NX*Ref->NY*sizeof(double));
      Ref->AY = (double*) malloc(Ref->NX*Ref->NY*sizeof(double));
      memcpy(Ref->AX,AX,Ref->NX*Ref->NY*sizeof(double));
      memcpy(Ref->AY,AY,Ref->NX*Ref->NY*sizeof(double));

      GeoRef_AxisCalcExpandCoeff(Ref);
      break;

    case 'G':
      Ref->RPNHeadExt.grref[0] = 'L';
      Ref->RPNHeadExt.igref1 = 0.0;
      Ref->RPNHeadExt.igref2 = 0.0;
      Ref->RPNHeadExt.igref3 = 1.0;
      Ref->RPNHeadExt.igref4 = 1.0;
      f77name(cxgaig)(Ref->RPNHeadExt.grref,&Ref->RPNHeadExt.igref1, &Ref->RPNHeadExt.igref2, &Ref->RPNHeadExt.igref3, &Ref->RPNHeadExt.igref4,
		      &Ref->RPNHeadExt.xgref1, &Ref->RPNHeadExt.xgref2, &Ref->RPNHeadExt.xgref3, &Ref->RPNHeadExt.xgref4,1);

      Ref->AX = (double*)malloc(Ref->NX*sizeof(double));
      dlon = 360.0 / (double) Ref->NX;
      for (i=0; i < Ref->NX; i++) {
	      Ref->AX[i] = (double)i * dlon;
	    }

      zero = 0;
      GeoRef_AxisCalcExpandCoeff(Ref);

      switch (Ref->RPNHead.ig1) {
	      case GLOBAL:
	        Ref->AY = (double*) malloc(Ref->NY*sizeof(double));
	        temp    = (double*) malloc(Ref->NY*sizeof(double));
	        f77name(ez_glat)(Ref->AY,temp,&Ref->NY,&zero);
	        free(temp);
	        break;

	      case NORTH:
	      case SOUTH:
	        deuxnj = 2 * Ref->NY;
	        Ref->AY = (double*) malloc(deuxnj*sizeof(double));
	        temp    = (double*) malloc(deuxnj*sizeof(double));
	        f77name(ez_glat)(Ref->AY,temp,&deuxnj,&zero);
	        free(temp);
	        break;
	      }


      GeoRef_AxisCalcNewtonCoeff(Ref);
      break;

    default:
      GeoRef_AxisCalcExpandCoeff(Ref);
      break;
    }
}

int GeoRef_AxisGetExpanded(TGeoRef* Ref,double *AX,double *AY) {
  
   int nix, njy;
   int istart, jstart;

   if (Ref->NbSub > 0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: This operation is not supported for 'U' grids\n",__func__);
      return(-1);
   }
  
   if (!Ref->AX) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Grid descriptor not found\n",__func__);
      return(-1);
   }

   switch(Ref->GRTYP[0]) {
      case 'Y':
         nix = Ref->NX * Ref->NY;
         memcpy(AX, Ref->AX, nix*sizeof(double));
         memcpy(AY, Ref->AY, nix*sizeof(double));
         break;
      
      default:
         nix = Ref->NX;
         njy = Ref->NY;
         if (Ref->i2 == (nix+1)) istart = 1;
         if (Ref->i2 == (nix+2)) istart = 2;
         if (Ref->i2 == (nix)) istart = 0;

         if (Ref->j2 == (njy+1)) jstart = 1;
         if (Ref->j2 == (njy+2)) jstart = 2;
         if (Ref->j2 == (njy))   jstart = 0;
         memcpy(&AX[istart],Ref->AX, nix*sizeof(double));
         memcpy(&AY[jstart],Ref->AY, njy*sizeof(double));
      
         if (Ref->i2 == (Ref->NX+1)) {
            AX[0] = Ref->AX[nix-2] - 360.0; 
            AX[nix] = AX[2];
         }
      
         if (Ref->i2 == (Ref->NX+2)) {
            AX[0] = Ref->AX[nix-1] - 360.0; 
            AX[nix] = Ref->AX[1]+360.0;
            AX[nix+1] = Ref->AX[2]+360.0;
         }
   }
   return(0);
}
