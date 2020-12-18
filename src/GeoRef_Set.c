/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Fonctions et definitions relatives aux fichiers standards et rmnlib
 * Fichier      : GeoRef_Set.c
 * Creation     : Avril 2006 - J.P. Gauthier
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

void GeoRef_SetZoneFree(TGridSet *GSet) {

  int i;

  for (i=0; i < SET_NZONES; i++) {
     if (GSet->zones[i].npts > 0) {
        free(GSet->zones[i].idx);
        free(GSet->zones[i].x);
        free(GSet->zones[i].y);
        GSet->zones[i].npts = 0;
        GSet->zones[i].idx  = NULL;
        GSet->zones[i].x    = NULL;
        GSet->zones[i].y    = NULL;
      }
   }
}

int GeoRef_SetZoneDefinePole(TGeoRef *RefFrom,TGridSet *GSet,int Zone) {

   TGeoZone *zone=&GSet->zones[Zone];
   double    latpole,lonpole,xpole,ypole;
   int       nhits,i;
   int      *tmpidx;
  
   // On commence par trouver les points au pole
   tmpidx = (int*)malloc(zone->npts*sizeof(int));
  
   nhits = 0;
   if (RefFrom->GRTYP[0] == 'Z' && RefFrom->RPNHead.GRREF[0] == 'E') {
      xpole = 0.5 * RefFrom->NX;
      ypole = (Zone==NORTH)?RefFrom->NY+0.5:0.5;
   } else {
      latpole = (Zone==NORTH)?90.0:-90.0;
      lonpole = 0.0;
      GeoRef_LL2XY(RefFrom,&xpole,&ypole,&latpole,&lonpole,1,TRUE);
   }
  
   for (i=0; i < zone->npts; i++) {
      if (fabs(GSet->Y[i]-ypole) < 1.0e-3) {
         tmpidx[nhits]=i;
         nhits++;
      }
   }
  
   zone->npts = nhits;
   if (nhits > 0) {
      zone->x = (double*)malloc(nhits*sizeof(double));
      zone->y = (double*)malloc(nhits*sizeof(double));
      zone->idx = (int *)malloc(nhits*sizeof(int));
      App_Log(DEBUG,"%s: Number of points at pole: %d\n",__func__,nhits);
      
      for (i=0; i < zone->npts; i++) {
         zone->x[i]   = GSet->X[tmpidx[i]];      
         zone->y[i]   = GSet->Y[tmpidx[i]];
         zone->idx[i] = tmpidx[i];
      }
   }
  
   free(tmpidx);

   return(0);
}

int GeoRef_SetZoneDefineThem(TGeoRef *RefFrom,TGridSet *GSet,int Zone) {

   TGeoZone *zone=&GSet->zones[Zone];
   int      *tmpidx;
   int       nhits,i,jlim;
  
   tmpidx = (int  *) malloc(zone->npts*sizeof(int));
  
   nhits = 0;
   jlim = (Zone==SOUTH)?RefFrom->j1+1:RefFrom->j2-2;
   for (i=0; i < zone->npts; i++) {
      if ((Zone==SOUTH)?((int)GSet->Y[i] < jlim):((int)GSet->Y[i] > jlim)) {
         tmpidx[nhits]=i;
         nhits++;
      }
   }
  
   zone->npts = nhits;
   if (nhits > 0) {
      zone->x = (double*)malloc(nhits*sizeof(double));
      zone->y = (double*)malloc(nhits*sizeof(double));
      zone->idx = (int *) malloc(nhits*sizeof(int));
      App_Log(DEBUG,"%s: Number of points between pole and limit: %d\n",__func__,nhits);
    
      for (i=0; i < zone->npts; i++) {
         zone->x[i]   = GSet->X[tmpidx[i]];      
         zone->y[i]   = GSet->Y[tmpidx[i]];
         zone->idx[i] = tmpidx[i];
      }
   }
  
   free(tmpidx);

   return(0);
}

int GeoRef_SetZoneDefineOut(TGeoRef *RefFrom,TGridSet *GSet,int Zone) {

   TGeoZone *zone=&GSet->zones[Zone];
   int      *tmpidx;
   int       nhits;
   int       i,offsetleft,offsetright,ix,iy;
    
   tmpidx = (int  *) malloc(zone->npts*sizeof(int));
  
   offsetright = 0;
   offsetleft = 0;
/*
   if (groptions.degre_interp == CUBIC) {
     offsetright = 2;
     offsetleft = 1;
   } else {
     offsetright = 0;
     offsetleft = 0;
   }
*/
  
   App_Log(DEBUG,"%s: Offset left: %d Offset right: \n",__func__,offsetleft, offsetright);
   nhits = 0;
   for (i=0; i < zone->npts; i++) {
      ix = (int)(GSet->X[i]+0.5);
      iy = (int)(GSet->Y[i]+0.5);
      if (ix < (1+offsetleft) || iy < (1+offsetleft) || ix > (RefFrom->NX-offsetright) || iy > (RefFrom->NY-offsetright)) {
          tmpidx[nhits]=i;
         nhits++;
      }
   }

   if (nhits > 0) {
      zone->npts = nhits;
      zone->x = (double*)malloc(nhits*sizeof(double));
      zone->y = (double*)malloc(nhits*sizeof(double));
      zone->idx = (int *) malloc(zone->npts*sizeof(int));
      App_Log(DEBUG,"%s: Number opf outside pointst: \n",__func__,offsetleft,zone->npts);
    
      for (i=0; i < zone->npts; i++) {
         zone->x[i]   = GSet->X[tmpidx[i]];      
         zone->y[i]   = GSet->Y[tmpidx[i]];
         zone->idx[i] = tmpidx[i];
      }
   }
  
   free(tmpidx);
  
   return(0);
}

int GeoRef_SetZoneDefine(TGeoRef *RefTo,TGeoRef *RefFrom) {

   TGridSet *gset=NULL;
   int       i;
   int       extrap;

   gset=GeoRef_SetGet(RefTo,RefFrom,NULL);
  
   if (gset->flags & SET_ZONES) {
      return(0);
   }

   extrap = FALSE;
   switch (RefFrom->GRTYP[0]) {
      case 'N':
      case 'S':
      case 'L':
      case '!':
         extrap = TRUE;
         break;
      
      case '#':
      case 'Z':
      case 'Y':
         switch(RefFrom->RPNHead.GRREF[0]) {
	          case 'N':
	          case 'S':
	          case 'L':
	             extrap = TRUE;
	             break;
	  
	          case 'E':
	             if (359.0 > (RefFrom->AX[RefFrom->NX-1] - RefFrom->AX[0])) {
	                extrap = TRUE;
	             }
	             break;
	       }
         break;
   }
  
   for (i=0; i < SET_NZONES; i++) {
      gset->zones[i].npts = 0;
   }
  
    if (extrap) {
       GeoRef_SetZoneDefineOut(RefFrom,gset,OUTSIDE);
    } else {
       GeoRef_SetZoneDefinePole(RefFrom,gset,NORTH_POLE);
       GeoRef_SetZoneDefinePole(RefFrom,gset,SOUTH_POLE);
       GeoRef_SetZoneDefineThem(RefFrom,gset,SOUTH);
       GeoRef_SetZoneDefineThem(RefFrom,gset,NORTH);
    }
  
    gset->flags |= SET_ZONES;
    return(0);
}

int GeoRef_SetCalcXY(TGeoRef *RefTo,TGeoRef *RefFrom) {

   TGridSet *gset=NULL;

   gset=GeoRef_SetGet(RefTo,RefFrom,NULL);

   if (!gset->Index) {
      gset->IndexSize=2*RefTo->NX*RefTo->NY;
      gset->Index = (double*)calloc(gset->IndexSize,sizeof(double));
      gset->X=gset->Index;
      gset->Y=&gset->Index[RefTo->NX*RefTo->NY];

      GeoRef_LL2XY(RefFrom,gset->X,gset->Y,RefTo->Lat,RefTo->Lon,RefTo->NX*RefTo->NY,TRUE);
   }

   return(0);
}

int GeoRef_SetCalcYYXY(TGeoRef *RefTo,TGeoRef *RefFrom) {

   TGridSet *gset=NULL;
   TGeoRef *yin_gdin, *yan_gdin, *yin_gdout, *yan_gdout;
   int icode,nij,i,j,k,ivalue,ni,nj,yni,ynj,yin_mgid;
   int idx_gdin;
   int yancount_yin,yincount_yin, yancount_yan,yincount_yan;
   int yyin,yyout;
   double *yin2yin_lat,*yin2yin_lon,*yan2yin_lat,*yan2yin_lon;
   double *yin2yan_lat,*yin2yan_lon,*yan2yan_lat,*yan2yan_lon;
   float *yin_fld;
    
   //  Need only access to either yin or Yang info for the lat and lon val
      
   yyin=0; yyout=0;

   gset=GeoRef_SetGet(RefTo,RefFrom,NULL);
   /* Mettre du code au cas ou gdx_gdin == -1 */

   if (gset->flags & SET_YYXY) {
      return 0;
   }

   /* Dans un premier temps on calcule la position x-y de tous les points sur la grille */

   /* To be in this routine, input source grid should be Yin-Yang */
   yyin=1;
   yin_gdin = RefFrom->Subs[0];
   yan_gdin = RefFrom->Subs[1];

   /* Check what the destination grid is */
   if (RefTo->NbSub > 0) {
      yyout=1;
      yin_gdout = RefTo->Subs[0];
      yan_gdout = RefTo->Subs[1];
      ni = yin_gdout->NX;
      nj = yin_gdout->NY;
   } else {
      yin_gdout = RefTo;
      ni = RefTo->NX;
      nj = RefTo->NY;
   }

   nij = ni*nj;

   /* Masquer les grilles YY input pour enlever overlap si TRUE */
   k=0;
   yin2yin_lat = (double*)malloc(4*nij*sizeof(double));
   yin2yin_lon = &yin2yin_lat[k+=nij];
   yan2yin_lat = &yin2yin_lat[k+=nij];
   yan2yin_lon = &yin2yin_lat[k+=nij];

   yancount_yin=0;
   yincount_yin=0;

   /* destination grid is one grid */ 
   /* create mask with Yin as a priority choice and store x,y,lat,lon pos */
   gset->yin_maskout = (float *) malloc(ni*nj*sizeof(float));
   gset->yinlat = (double*) malloc(ni*nj*sizeof(double));
   gset->yinlon = (double*) malloc(ni*nj*sizeof(double));
   icode = GeoRef_GetLL(yin_gdout,gset->yinlat,gset->yinlon);
   icode = GeoRef_MaskYYApply(yin_gdout,yin_gdin,ni,nj,gset->yin_maskout,gset->yinlat,gset->yinlon,yin2yin_lat,yin2yin_lon,&yincount_yin,yan2yin_lat,yan2yin_lon,&yancount_yin);
   /* store the lats and lons */
   gset->yincount_yin = yincount_yin;
   gset->yancount_yin = yancount_yin;
   gset->yin2yin_lat = (double*) malloc(yincount_yin*sizeof(double));
   gset->yin2yin_lon = (double*) malloc(yincount_yin*sizeof(double));
   gset->yan2yin_lat = (double*) malloc(yancount_yin*sizeof(double));
   gset->yan2yin_lon = (double*) malloc(yancount_yin*sizeof(double));
   memcpy(gset->yin2yin_lat,yin2yin_lat,yincount_yin*sizeof(double));
   memcpy(gset->yin2yin_lon,yin2yin_lon,yincount_yin*sizeof(double));
   memcpy(gset->yan2yin_lat,yan2yin_lat,yancount_yin*sizeof(double));
   memcpy(gset->yan2yin_lon,yan2yin_lon,yancount_yin*sizeof(double));

   /* store the Xs and Ys */
   gset->yin2yin_x = (double*) malloc(yincount_yin*sizeof(double));
   gset->yin2yin_y = (double*) malloc(yincount_yin*sizeof(double));
   gset->yan2yin_x = (double*) malloc(yancount_yin*sizeof(double));
   gset->yan2yin_y = (double*) malloc(yancount_yin*sizeof(double));
   icode = GeoRef_LL2XY(yin_gdin,gset->yin2yin_x,gset->yin2yin_y,yin2yin_lat,yin2yin_lon,yincount_yin,TRUE);
   icode = GeoRef_LL2XY(yan_gdin,gset->yan2yin_x,gset->yan2yin_y,yan2yin_lat,yan2yin_lon,yancount_yin,TRUE);

   free(yin2yin_lat);

   if (yyout == 1) { 
      k=0;
      yin2yan_lat = (double*)malloc(4*nij*sizeof(double));
      yin2yan_lon = &yin2yan_lat[k+=nij];
      yan2yan_lat = &yin2yan_lat[k+=nij];
      yan2yan_lon = &yin2yan_lat[k+=nij];

      /* create mask (Yin priority) with src Yin,src Yang onto dest Yang and store x,y pos */

      gset->yan_maskout = (float *) malloc(ni*nj*sizeof(float));
      gset->yanlat = (double*) malloc(ni*nj*sizeof(double));
      gset->yanlon = (double*) malloc(ni*nj*sizeof(double));
      icode = GeoRef_GetLL(yan_gdout,gset->yanlat,gset->yanlon);
      icode = GeoRef_MaskYYApply(yan_gdout,yin_gdin,ni,nj,gset->yan_maskout,gset->yanlat,gset->yanlon,yin2yan_lat,yin2yan_lon,&yincount_yan,yan2yan_lat,yan2yan_lon,&yancount_yan);
      gset->yincount_yan = yincount_yan;
      gset->yancount_yan = yancount_yan;
      gset->yin2yan_lat = (double*) malloc(yincount_yan*sizeof(double));
      gset->yin2yan_lon = (double*) malloc(yincount_yan*sizeof(double));
      gset->yan2yan_lat = (double*) malloc(yancount_yan*sizeof(double));
      gset->yan2yan_lon = (double*) malloc(yancount_yan*sizeof(double));
      memcpy(gset->yin2yan_lat,yin2yan_lat,yincount_yan*sizeof(double));
      memcpy(gset->yin2yan_lon,yin2yan_lon,yincount_yan*sizeof(double));
      memcpy(gset->yan2yan_lat,yan2yan_lat,yancount_yan*sizeof(double));
      memcpy(gset->yan2yan_lon,yan2yan_lon,yancount_yan*sizeof(double));
      gset->yin2yan_x = (double*) malloc(yincount_yan*sizeof(double));
      gset->yin2yan_y = (double*) malloc(yincount_yan*sizeof(double));
      gset->yan2yan_x = (double*) malloc(yancount_yan*sizeof(double));
      gset->yan2yan_y = (double*) malloc(yancount_yan*sizeof(double));
      icode = GeoRef_LL2XY(yin_gdin,gset->yin2yan_x,gset->yin2yan_y,yin2yan_lat,yin2yan_lon,yincount_yan,TRUE);
      icode = GeoRef_LL2XY(yan_gdin,gset->yan2yan_x,gset->yan2yan_y,yan2yan_lat,yan2yan_lon,yancount_yan,TRUE);

      free(yin2yan_lat);
   }

   gset->flags |= SET_YYXY;

   return(icode);
}

void GeoRef_SetFree(TGridSet* GSet) {

   int i;

//TODO: Check to free
   if (GSet->RefFrom) {
      GSet->RefFrom == NULL;
   }

   if (GSet->Index) {
      free(GSet->Index);
      GSet->Index=NULL;
      GSet->X=NULL;
      GSet->Y=NULL;
   }

   GeoRef_SetZoneFree(GSet);

 //TODO: to free:
 //   int *mask_in, *mask_out;
 // float *yin_maskout,*yan_maskout;
 // float *yinlat,*yinlon,*yanlat,*yanlon;
 // float *yin2yin_lat,*yin2yin_lon,*yan2yin_lat,*yan2yin_lon;
 // float *yin2yan_lat,*yin2yan_lon,*yan2yan_lat,*yan2yan_lon;
 // float *yin2yin_x,*yin2yin_y,*yan2yin_x,*yan2yin_y;
 // float *yin2yan_x,*yin2yan_y,*yan2yan_x,*yan2yan_y;
}

int GeoRef_SetRead(int FID,TGridSet *GSet,char GFrom,char GTo){

   TRPNHeader h;
   char       typvar[2];

   typvar[0]=GFrom;
   typvar[1]=GTo;

   // Rechercher et lire l'information de l'enregistrement specifie
   if ((h.KEY=cs_fstinf(FID,&h.NI,&h.NJ,&h.NK,-1,"GRIDSET",-1,-1,-1,typvar,"#>>#"))<0) {
      App_Log(ERROR,"%s: Could not find gridset index field (c_fstinf failed)\n",__func__);
      return(FALSE);
   }
   GSet->IndexSize=h.NI;
   GSet->Index=(double*)malloc(GSet->IndexSize*sizeof(double));

   c_fst_data_length(8);
   if (cs_fstluk(GSet->Index,h.KEY,&h.NI,&h.NJ,&h.NK)<0) {
      App_Log(ERROR,"%s: Could not read gridset index field (c_fstlir failed)\n",__func__);
      return(FALSE);
   }
   return(TRUE);
}

int GeoRef_SetWrite(int FID,TGridSet *GSet){


   if (GSet && GSet->Index) {
      c_fst_data_length(8);

      if (cs_fstecr(GSet->Index,-64,FID,0,0,0,GSet->IndexSize,1,1,0,0,0,GSet->G2G,"#>>#","GRIDSET","X",0,0,0,0,5,FALSE)<0) {
         App_Log(ERROR,"%s: Could not write gridset index field (c_fstecr failed)\n",__func__);
         return(FALSE);
      }
   }
   return(TRUE);
}

TGridSet* GeoRef_SetGet(TGeoRef* RefTo,TGeoRef* RefFrom,TGridSet** GSet) {

   TGridSet* gset=NULL;
   int i;

   if (!RefTo || !RefFrom) {
      return(NULL);
   }

   if (!RefTo->Sets) {
      RefTo->Sets=(TGridSet*)calloc(sizeof(TGridSet),SET_MAX);
   }

   // Check for last set (most cases)
   if (RefTo->LastSet && RefTo->LastSet->RefFrom == RefFrom) {
      return(RefTo->LastSet);
   }

   if (RefTo->NbSet>=SET_MAX) {
      for (i=0; i< RefTo->NbSet; i++) {
         if (RefTo->Sets[i].RefFrom!=NULL) {
            GeoRef_SetFree(&RefTo->Sets[i]);
         }
      }
      RefTo->NbSet=0;
   }

   // Otherwise loop on sets
   i=0;
   while (i<RefTo->NbSet) {
      if (RefTo->Sets[i].RefFrom==RefFrom) {
         RefTo->LastSet = &RefTo->Sets[i];
         return(&RefTo->Sets[i]);
      }
      i++;
   }

   // If we get here, we have'nt found any sets, create a new one    
   RefTo->Sets[i].RefFrom = RefFrom;
   RefTo->NbSet++;
   RefTo->Sets[i].G2G[0]=RefFrom->GRTYP[0];
   RefTo->Sets[i].G2G[1]=RefTo->GRTYP[0];
   
   // If a known grid set is passed in, use it
   if (GSet) {
      if (gset=*GSet) {
         if (gset->G2G[0]!=RefTo->Sets[i].G2G[0] || gset->G2G[1]!=RefTo->Sets[i].G2G[1]) {
            App_Log(WARNING,"%s: Invalid grid set %c%c!=%c%c\n",__func__,gset->G2G[0],gset->G2G[1],RefTo->Sets[i].G2G[0],RefTo->Sets[i].G2G[1]);
         } else {
            RefTo->Sets[i].Index=gset->Index;
            RefTo->Sets[i].IndexSize=gset->IndexSize;
         }
      } else {
         *GSet=&RefTo->Sets[i];
      }
   }

   App_Log(DEBUG,"%s: RefFrom : %p RefTo: %p\n",__func__,RefFrom,RefTo);

   return(&RefTo->Sets[i]);
}
