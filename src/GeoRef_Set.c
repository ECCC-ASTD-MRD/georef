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
      GeoRef_LL2XY(RefFrom,&xpole,&ypole,&latpole,&lonpole,1);
   }
  
   for (i=0; i < zone->npts; i++) {
      if (fabs(GSet->y[i]-ypole) < 1.0e-3) {
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
         zone->x[i]   = GSet->x[tmpidx[i]];      
         zone->y[i]   = GSet->y[tmpidx[i]];
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
      if ((Zone==SOUTH)?((int)GSet->y[i] < jlim):((int)GSet->y[i] > jlim)) {
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
         zone->x[i]   = GSet->x[tmpidx[i]];      
         zone->y[i]   = GSet->y[tmpidx[i]];
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
      ix = (int)(GSet->x[i]+0.5);
      iy = (int)(GSet->y[i]+0.5);
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
         zone->x[i]   = GSet->x[tmpidx[i]];      
         zone->y[i]   = GSet->y[tmpidx[i]];
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

   gset=GeoRef_SetGet(RefTo,RefFrom);
  
   if (gset->flags & SET_ZONES) {
      return 0;
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
   int       ninj_out;

   gset=GeoRef_SetGet(RefTo,RefFrom);

   if (gset->x) {
      return(0);
   }

   ninj_out = RefTo->NX*RefTo->NY;

   gset->x = (double*)malloc(ninj_out*sizeof(double));
   gset->y = (double*)malloc(ninj_out*sizeof(double));

   GeoRef_LL2XY(RefFrom,gset->x,gset->y,RefTo->Lat,RefTo->Lon,ninj_out);
   
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

   gset=GeoRef_SetGet(RefTo,RefFrom);
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

   k=0;
   nij = ni*nj;

   /* Masquer les grilles YY input pour enlever overlap si TRUE */
   yin2yin_lat = (double*)malloc(8*nij*sizeof(double));
   yin2yin_lon = &yin2yin_lat[k+=nij];
   yan2yin_lat = &yin2yin_lat[k+=nij];
   yan2yin_lon = &yin2yin_lat[k+=nij];
   yin2yan_lat = &yin2yin_lat[k+=nij];
   yin2yan_lon = &yin2yin_lat[k+=nij];
   yan2yan_lat = &yin2yin_lat[k+=nij];
   yan2yan_lon = &yin2yin_lat[k+=nij];
   yancount_yin=0;
   yincount_yin=0;
   if (yyout == 0) {
      /* destination grid is one grid */ 
      /* create mask with Yin as a priority choice and store x,y,lat,lon pos */
      gset->yin_maskout = (float *) malloc(ni*nj*sizeof(float));
      gset->yinlat = (double*) malloc(ni*nj*sizeof(double));
      gset->yinlon = (double*) malloc(ni*nj*sizeof(double));
      icode = GeoRef_GetLL(yin_gdout,gset->yinlat,gset->yinlon);
      icode = c_ezyymint(yin_gdout,yin_gdin, i,nj,gset->yin_maskout,gset->yinlat,gset->yinlon,yin2yin_lat,yin2yin_lon,&yincount_yin,yan2yin_lat,yan2yin_lon,&yancount_yin);
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
      icode = GeoRef_LL2XY(yin_gdin,gset->yin2yin_x,gset->yin2yin_y,yin2yin_lat,yin2yin_lon,yincount_yin);
      icode = GeoRef_LL2XY(yan_gdin,gset->yan2yin_x,gset->yan2yin_y,yan2yin_lat,yan2yin_lon,yancount_yin);
   }

   if (yyout == 1) { 
      /* destination grid is a U grid*/
      /* create mask (Yin priority)with src Yin,src Yang onto dest Yin and 
                                                            store x,y pos */
      gset->yin_maskout = (float *) malloc(ni*nj*sizeof(float));
      gset->yinlat = (double *) malloc(ni*nj*sizeof(double));
      gset->yinlon = (double *) malloc(ni*nj*sizeof(double));
      icode = GeoRef_GetLL(yin_gdout,gset->yinlat,gset->yinlon);
      icode = c_ezyymint(yin_gdout,yin_gdin,ni,nj,gset->yin_maskout,gset->yinlat,gset->yinlon,yin2yin_lat,yin2yin_lon,&yincount_yin,yan2yin_lat,yan2yin_lon,&yancount_yin);
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
      gset->yin2yin_x = (double*) malloc(yincount_yin*sizeof(double));
      gset->yin2yin_y = (double*) malloc(yincount_yin*sizeof(double));
      gset->yan2yin_x = (double*) malloc(yancount_yin*sizeof(double));
      gset->yan2yin_y = (double*) malloc(yancount_yin*sizeof(double));
      icode = GeoRef_LL2XY(yin_gdin,gset->yin2yin_x,gset->yin2yin_y,yin2yin_lat,yin2yin_lon,yincount_yin);
      icode = GeoRef_LL2XY(yan_gdin,gset->yan2yin_x,gset->yan2yin_y,yan2yin_lat,yan2yin_lon,yancount_yin);
      
      /* create mask (Yin priority) with src Yin,src Yang onto dest Yang and store x,y pos */

      gset->yan_maskout = (float *) malloc(ni*nj*sizeof(float));
      gset->yanlat = (double*) malloc(ni*nj*sizeof(double));
      gset->yanlon = (double*) malloc(ni*nj*sizeof(double));
      icode = GeoRef_GetLL(yan_gdout,gset->yanlat,gset->yanlon);
      icode = c_ezyymint(yan_gdout,yin_gdin,ni,nj,gset->yan_maskout,gset->yanlat,gset->yanlon,yin2yan_lat,yin2yan_lon,&yincount_yan,yan2yan_lat,yan2yan_lon,&yancount_yan);
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
      icode = GeoRef_LL2XY(yin_gdin,gset->yin2yan_x,gset->yin2yan_y,yin2yan_lat,yin2yan_lon,yincount_yan);
      icode = GeoRef_LL2XY(yan_gdin,gset->yan2yan_x,gset->yan2yan_y,yan2yan_lat,yan2yan_lon,yancount_yan);
   }

   free(yin2yin_lat);
   gset->flags |= SET_YYXY;

   return(icode);
}

void GeoRef_SetFree(TGridSet* GSet) {

   int i;

//TODO: Check to free
   if (GSet->RefFrom) {
      GSet->RefFrom == NULL;
   }

   if (GSet->x) {
      free(GSet->x);
      GSet->x  = NULL;
   }

   if (GSet->y) {
      free(GSet->y);
      GSet->y  = NULL;
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

TGridSet* GeoRef_SetGet(TGeoRef* RefTo, TGeoRef* RefFrom) {

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

   App_Log(DEBUG,"%s: RefFrom : %p RefTo: %p\n",__func__,RefFrom,RefTo);

   return(&RefTo->Sets[i]);
}
