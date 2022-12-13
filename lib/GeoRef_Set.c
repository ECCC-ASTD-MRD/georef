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

/*----------------------------------------------------------------------------
 * @brief  Free gridset zone definitions
 * @author Yves Chartier
 * @date   
 *    @param[in]  GSet    Grid set pointer
*/
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

/*----------------------------------------------------------------------------
 * @brief  Finds the points at a pole
 * @author Yves Chartier
 * @date   
 *    @param[in]  GSet    Grid set pointer
 *    @param[in]  Zone    Zone identifier (NORTH_POLE,SOUTH_POLE)
 *    @param[in]  NbPts   Number of points to check
 *
 *    @return             Error code (0=ok)
*/
int GeoRef_SetZoneDefinePole(TGridSet *GSet,int Zone,int NbPts) {

   TGeoZone *zone=&GSet->zones[Zone];
   double    latpole,lonpole,xpole,ypole;
   int      *tmpidx,i;
  
   // On commence par trouver les points au pole
   tmpidx = (int*)malloc(NbPts*sizeof(int));
  
   zone->npts = 0;
   if (GSet->RefFrom->GRTYP[0] == 'Z' && GSet->RefFrom->RPNHead.GRREF[0] == 'E') {
      xpole = 0.5 * GSet->RefFrom->NX;
      ypole = (Zone==NORTH)?GSet->RefFrom->NY+0.5:0.5;
   } else {
      latpole = (Zone==NORTH)?90.0:-90.0;
      lonpole = 0.0;
      GeoRef_LL2XY(GSet->RefFrom,&xpole,&ypole,&latpole,&lonpole,1,TRUE);
   }
  
   for (i=0; i<NbPts; i++) {
      if (fabs(GSet->Y[i]-ypole) < 1.0e-3) {
         tmpidx[zone->npts]=i;
         zone->npts++;
      }
   }
  
   if (zone->npts>0) {
      zone->x = (double*)malloc(zone->npts*sizeof(double));
      zone->y = (double*)malloc(zone->npts*sizeof(double));
      zone->idx = (int *)malloc(zone->npts*sizeof(int));
      Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Number of points at pole: %d\n",__func__,zone->npts);
      
      for (i=0; i<zone->npts; i++) {
         zone->x[i]   = GSet->X[tmpidx[i]];      
         zone->y[i]   = GSet->Y[tmpidx[i]];
         zone->idx[i] = tmpidx[i];
      }
   }
  
   free(tmpidx);

   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Finds the points between the pole and a limit
 * @author Yves Chartier
 * @date   
 *    @param[in]  GSet    Grid set pointer
 *    @param[in]  Zone    Zone identifier (NORTH,SOUTH)
 *    @param[in]  NbPts   Number of points to check
 *
 *    @return             Error code (0=ok)
*/
int GeoRef_SetZoneDefineThem(TGridSet *GSet,int Zone,int NbPts) {

   TGeoZone *zone=&GSet->zones[Zone];
   int      *tmpidx,i,jlim;
  
   tmpidx = (int*) malloc(NbPts*sizeof(int));
  
   zone->npts = 0;
   jlim = (Zone==SOUTH)?GSet->RefFrom->j1+1:GSet->RefFrom->j2-2;
   for (i=0; i<NbPts; i++) {
      if ((Zone==SOUTH)?((int)GSet->Y[i] < jlim):((int)GSet->Y[i] > jlim)) {
         tmpidx[zone->npts]=i;
         zone->npts++;
      }
   }
  
   if (zone->npts>0) {
      zone->x = (double*)malloc(zone->npts*sizeof(double));
      zone->y = (double*)malloc(zone->npts*sizeof(double));
      zone->idx = (int *) malloc(zone->npts*sizeof(int));
      Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Number of points between pole and limit: %d\n",__func__,zone->npts);
    
      for (i=0; i<zone->npts; i++) {
         zone->x[i]   = GSet->X[tmpidx[i]];      
         zone->y[i]   = GSet->Y[tmpidx[i]];
         zone->idx[i] = tmpidx[i];
      }
   }
  
   free(tmpidx);

   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Finds the points outside of the source data
 * @author Yves Chartier
 * @date   
 *    @param[in]  GSet    Grid set pointer
 *    @param[in]  Zone    Zone identifier (OUTSIDE)
 *    @param[in]  NbPts   Number of points to check
 *
 *    @return             Error code (0=ok)
*/
int GeoRef_SetZoneDefineOut(TGridSet *GSet,int Zone,int NbPts) {

   TGeoZone *zone=&GSet->zones[Zone];
   int      *tmpidx,i,offsetleft,offsetright,ix,iy;
    
   tmpidx=(int*)malloc(NbPts*sizeof(int));
  
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

   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: NbPoints %d, Offset left: %d, Offset right: %d\n",__func__,NbPts,offsetleft, offsetright);
   zone->npts=0;
   for (i=0; i<NbPts; i++) {
      ix = (int)(GSet->X[i]+0.5);
      iy = (int)(GSet->Y[i]+0.5);
      if (ix < (1+offsetleft) || iy < (1+offsetleft) || ix > (GSet->RefFrom->NX-offsetright) || iy > (GSet->RefFrom->NY-offsetright)) {
         tmpidx[zone->npts]=i;
         zone->npts++;
      }
   }

   if (zone->npts>0) {
      zone->x = (double*)malloc(zone->npts*sizeof(double));
      zone->y = (double*)malloc(zone->npts*sizeof(double));
      zone->idx = (int *) malloc(zone->npts*sizeof(int));
      Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Number of outside pointst: \n",__func__,offsetleft,zone->npts);
    
      for (i=0; i < zone->npts; i++) {
         zone->x[i]   = GSet->X[tmpidx[i]];      
         zone->y[i]   = GSet->Y[tmpidx[i]];
         zone->idx[i] = tmpidx[i];
      }
   }
  
   free(tmpidx);
  
   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Defines the various zones
 * @author Yves Chartier
 * @date   
 *    @param[in]  GridSet   GridSet
 *
 *    @return             Error code (0=ok)
*/
int GeoRef_SetZoneDefine(TGridSet *GSet) {

   int       i,npts;
   int       extrap;
  
   if (!GSet || (GSet->flags & SET_ZONES)) {
      return(0);
   }

   extrap = FALSE;
   switch (GSet->RefFrom->GRTYP[0]) {
      case 'N':
      case 'S':
      case 'L':
      case '!':
         extrap = TRUE;
         break;
      
      case '#':
      case 'Z':
      case 'Y':
         switch(GSet->RefFrom->RPNHead.GRREF[0]) {
	          case 'N':
	          case 'S':
	          case 'L':
	             extrap = TRUE;
	             break;
	  
	          case 'E':
	             if (359.0 > (GSet->RefFrom->AX[GSet->RefFrom->NX-1] - GSet->RefFrom->AX[0])) {
	                extrap = TRUE;
	             }
	             break;
	       }
         break;
   }
  
   for (i=0; i<SET_NZONES; i++) {
      GSet->zones[i].npts = 0;
   }
   npts = GSet->RefTo->NX * GSet->RefTo->NY;

   if (extrap) {
      GeoRef_SetZoneDefineOut(GSet,OUTSIDE,npts);
   } else {
      GeoRef_SetZoneDefinePole(GSet,NORTH_POLE,npts);
      GeoRef_SetZoneDefinePole(GSet,SOUTH_POLE,npts);
      GeoRef_SetZoneDefineThem(GSet,SOUTH,npts);
      GeoRef_SetZoneDefineThem(GSet,NORTH,npts);
   }
  
   GSet->flags |= SET_ZONES;
   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Calculates XY correspondance of destination points within the source grid
 * @author Yves Chartier
 * @date   
 *    @param[in]  GridSet   GridSet
 *
 *    @return             Error code (0=ok)
*/
int GeoRef_SetCalcXY(TGridSet *GSet) {

   int size=0;

   if (GSet && !GSet->Index) {
      size=GSet->RefTo->NX*GSet->RefTo->NY;
      GSet->IndexDegree=GSet->RefFrom->Options.InterpDegree;
      GSet->Index = (float*)calloc((GSet->IndexDegree==IR_CUBIC?10:(GSet->IndexDegree==IR_LINEAR?6:1))*size,sizeof(float));
      GSet->X = (double*)calloc(size*2,sizeof(double));
      GSet->Y=&GSet->X[size];

      GeoRef_LL2XY(GSet->RefFrom,GSet->X,GSet->Y,GSet->RefTo->Lat,GSet->RefTo->Lon,size,TRUE);
   }

   // Reset index
   if (GSet->IndexDegree!=GSet->RefFrom->Options.InterpDegree) {
      GSet->Index[0]=0;
   }

   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Calculates XY correspondance of destination points within the source grid (for YY grids)
 * @author Yves Chartier
 * @date   
 *    @param[in]  GridSet   GridSet
 *
 *    @return             Error code (0=ok)
*/
int GeoRef_SetCalcYYXY(TGridSet *GSet) {

   TGeoRef *yin_gdin,*yan_gdin,*yin_gdout,*yan_gdout;
   int icode,k,nij,ni,nj;
   int yancount_yin,yincount_yin, yancount_yan,yincount_yan;
   double *yin2yin_lat,*yin2yin_lon,*yan2yin_lat,*yan2yin_lon;
   double *yin2yan_lat,*yin2yan_lon,*yan2yan_lat,*yan2yan_lon;
    
   //  Need only access to either yin or Yang info for the lat and lon val
      
   /* Mettre du code au cas ou gdx_gdin == -1 */

   if (!GSet || (GSet->flags & SET_YYXY)) {
      return 0;
   }

   /* Dans un premier temps on calcule la position x-y de tous les points sur la grille */

   /* To be in this routine, input source grid should be Yin-Yang */
   yin_gdin = GSet->RefFrom->Subs[0];
   yan_gdin = GSet->RefFrom->Subs[1];

   /* Check what the destination grid is */
   if (GSet->RefTo->NbSub > 0) {
      yin_gdout = GSet->RefTo->Subs[0];
      yan_gdout = GSet->RefTo->Subs[1];
      ni = yin_gdout->NX;
      nj = yin_gdout->NY;
   } else {
      yin_gdout = GSet->RefTo;
      ni = GSet->RefTo->NX;
      nj = GSet->RefTo->NY;
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
   GSet->yin_maskout = (float *) malloc(nij*sizeof(float));
   GSet->yinlat = (double*) malloc(nij*sizeof(double));
   GSet->yinlon = (double*) malloc(nij*sizeof(double));
   icode = GeoRef_GetLL(yin_gdout,GSet->yinlat,GSet->yinlon);
   icode = GeoRef_MaskYYApply(yin_gdout,yin_gdin,ni,nj,GSet->yin_maskout,GSet->yinlat,GSet->yinlon,yin2yin_lat,yin2yin_lon,&yincount_yin,yan2yin_lat,yan2yin_lon,&yancount_yin);
   /* store the lats and lons */
   GSet->yincount_yin = yincount_yin;
   GSet->yancount_yin = yancount_yin;
   GSet->yin2yin_lat = (double*) malloc(yincount_yin*sizeof(double));
   GSet->yin2yin_lon = (double*) malloc(yincount_yin*sizeof(double));
   GSet->yan2yin_lat = (double*) malloc(yancount_yin*sizeof(double));
   GSet->yan2yin_lon = (double*) malloc(yancount_yin*sizeof(double));
   memcpy(GSet->yin2yin_lat,yin2yin_lat,yincount_yin*sizeof(double));
   memcpy(GSet->yin2yin_lon,yin2yin_lon,yincount_yin*sizeof(double));
   memcpy(GSet->yan2yin_lat,yan2yin_lat,yancount_yin*sizeof(double));
   memcpy(GSet->yan2yin_lon,yan2yin_lon,yancount_yin*sizeof(double));

   /* store the Xs and Ys */
   GSet->yin2yin_x = (double*) malloc(yincount_yin*sizeof(double));
   GSet->yin2yin_y = (double*) malloc(yincount_yin*sizeof(double));
   GSet->yan2yin_x = (double*) malloc(yancount_yin*sizeof(double));
   GSet->yan2yin_y = (double*) malloc(yancount_yin*sizeof(double));
   icode = GeoRef_LL2XY(yin_gdin,GSet->yin2yin_x,GSet->yin2yin_y,GSet->yin2yin_lat,GSet->yin2yin_lon,GSet->yincount_yin,TRUE);
   icode = GeoRef_LL2XY(yan_gdin,GSet->yan2yin_x,GSet->yan2yin_y,GSet->yan2yin_lat,GSet->yan2yin_lon,GSet->yancount_yin,TRUE);

   free(yin2yin_lat);

   // If destination grid is YY
   if (GSet->RefTo->NbSub > 0) {

      k=0;
      yin2yan_lat = (double*)malloc(4*nij*sizeof(double));
      yin2yan_lon = &yin2yan_lat[k+=nij];
      yan2yan_lat = &yin2yan_lat[k+=nij];
      yan2yan_lon = &yin2yan_lat[k+=nij];

      /* create mask (Yin priority) with src Yin,src Yang onto dest Yang and store x,y pos */

      GSet->yan_maskout = (float *) malloc(nij*sizeof(float));
      GSet->yanlat = (double*) malloc(nij*sizeof(double));
      GSet->yanlon = (double*) malloc(nij*sizeof(double));
      icode = GeoRef_GetLL(yan_gdout,GSet->yanlat,GSet->yanlon);
      icode = GeoRef_MaskYYApply(yan_gdout,yin_gdin,ni,nj,GSet->yan_maskout,GSet->yanlat,GSet->yanlon,yin2yan_lat,yin2yan_lon,&yincount_yan,yan2yan_lat,yan2yan_lon,&yancount_yan);
      GSet->yincount_yan = yincount_yan;
      GSet->yancount_yan = yancount_yan;
      GSet->yin2yan_lat = (double*) malloc(yincount_yan*sizeof(double));
      GSet->yin2yan_lon = (double*) malloc(yincount_yan*sizeof(double));
      GSet->yan2yan_lat = (double*) malloc(yancount_yan*sizeof(double));
      GSet->yan2yan_lon = (double*) malloc(yancount_yan*sizeof(double));
      memcpy(GSet->yin2yan_lat,yin2yan_lat,yincount_yan*sizeof(double));
      memcpy(GSet->yin2yan_lon,yin2yan_lon,yincount_yan*sizeof(double));
      memcpy(GSet->yan2yan_lat,yan2yan_lat,yancount_yan*sizeof(double));
      memcpy(GSet->yan2yan_lon,yan2yan_lon,yancount_yan*sizeof(double));
      GSet->yin2yan_x = (double*) malloc(yincount_yan*sizeof(double));
      GSet->yin2yan_y = (double*) malloc(yincount_yan*sizeof(double));
      GSet->yan2yan_x = (double*) malloc(yancount_yan*sizeof(double));
      GSet->yan2yan_y = (double*) malloc(yancount_yan*sizeof(double));
      icode = GeoRef_LL2XY(yin_gdin,GSet->yin2yan_x,GSet->yin2yan_y,yin2yan_lat,yin2yan_lon,yincount_yan,TRUE);
      icode = GeoRef_LL2XY(yan_gdin,GSet->yan2yan_x,GSet->yan2yan_y,yan2yan_lat,yan2yan_lon,yancount_yan,TRUE);

      free(yin2yan_lat);
   }

   GSet->flags |= SET_YYXY;

   return(icode);
}

/*----------------------------------------------------------------------------
 * @brief  Free grid set structure
 * @author Yves Chartier
 * @date   
 *    @param[in]  GSet      Gridset pointer
 *
*/
void GeoRef_SetFree(TGridSet* GSet) {

   int i;

//TODO: Check to free
   if (GSet->RefFrom) {
      GSet->RefFrom == NULL;
   }

   if (GSet->Index) {
      free(GSet->Index);
      GSet->Index=NULL;
      GSet->IndexDegree=IR_UNDEF;
   }
   if (GSet->X) {
      free(GSet->X);
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

/*----------------------------------------------------------------------------
 * @brief  Reads a gridset definition and index from a file
 * @author Jean-Philippe Gauthier
 * @date   January 2020
 *    @param[in]  FID       FSTD file identifier
 *    @param[in]  GSet      Gridset pointer
 *    @param[in]  Grom      Destination grid type
 *    @param[in]  GTo       Source grid type
 *    @param[in]  Level     Interpolation level (2=bilinear,3=bicubic,-1=any)
 *
 *    @return             Error code (0=ok)
*/
int GeoRef_SetRead(int FID,TGridSet *GSet,char GFrom,char GTo,int Level){

   TRPNHeader h;
   char       typvar[2];

   typvar[0]=GFrom;
   typvar[1]=GTo;

   // Rechercher et lire l'information de l'enregistrement specifie
   if ((h.KEY=cs_fstinf(FID,&h.NI,&h.NJ,&h.NK,-1,"GRIDSET",-1,-1,Level,typvar,"####"))<0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not find gridset index field (c_fstinf failed)\n",__func__);
      return(FALSE);
   }
   GSet->IndexDegree=(TDef_InterpR)Level;
   GSet->Index=(float*)malloc(h.NI*h.NJ*sizeof(double));

   if (cs_fstluk(GSet->Index,h.KEY,&h.NI,&h.NJ,&h.NK)<0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not read gridset index field (c_fstlir failed)\n",__func__);
      return(FALSE);
   }

   if ((h.KEY=cs_fstinf(FID,&h.NI,&h.NJ,&h.NK,-1,"GRIDSET",-1,-1,Level,typvar,"#>>#"))<0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not find gridset longitude field (c_fstinf failed)\n",__func__);
      return(FALSE);
   }
   GSet->X=(double*)malloc(2*h.NI*sizeof(double));
   GSet->Y=&GSet->X[h.NI];

   c_fst_data_length(8);
   if (cs_fstluk(GSet->X,h.KEY,&h.NI,&h.NJ,&h.NK)<0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not read gridset longitude (c_fstlir failed)\n",__func__);
      return(FALSE);
   }
   if ((h.KEY=cs_fstinf(FID,&h.NI,&h.NJ,&h.NK,-1,"GRIDSET",-1,-1,Level,typvar,"#^^#"))<0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not find gridset laitude field (c_fstinf failed)\n",__func__);
      return(FALSE);
   }
   c_fst_data_length(8);
   if (cs_fstluk(GSet->Y,h.KEY,&h.NI,&h.NJ,&h.NK)<0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not read gridset latitude field (c_fstlir failed)\n",__func__);
      return(FALSE);
   }
    return(TRUE);
}

/*----------------------------------------------------------------------------
 * @brief  Writes a gridset definition and index from a file
 * @author Jean-Philippe Gauthier
 * @date   January 2020
 *    @param[in]  GSet      Gridset pointer
 *    @param[in]  FID       FSTD file identifier
 *
 *    @return             Error code (0=ok)
*/
int GeoRef_SetWrite(TGridSet *GSet,int FID){

   int size=0;

   if (GSet && GSet->Index) {
      size=GSet->RefTo->NX*GSet->RefTo->NY;
      Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s:  Writing index (%ix%i)\n",__func__,size,GSet->IndexDegree==IR_CUBIC?10:(GSet->IndexDegree==IR_LINEAR?6:1));

      c_fst_data_length(8);
      if (cs_fstecr(GSet->X,-64,FID,0,0,0,size,1,1,0,0,0,GSet->G2G,"#>>#","GRIDSET","X",0,0,0,0,5,FALSE)<0) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write gridset index field (c_fstecr failed)\n",__func__);
         return(FALSE);
      }
      c_fst_data_length(8);
      if (cs_fstecr(GSet->Y,-64,FID,0,0,0,size,1,1,0,0,0,GSet->G2G,"#^^#","GRIDSET","X",0,0,0,0,5,FALSE)<0) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write gridset index field (c_fstecr failed)\n",__func__);
         return(FALSE);
      }
      if (cs_fstecr(GSet->Index,-32,FID,0,0,0,size,GSet->IndexDegree==IR_CUBIC?10:(GSet->IndexDegree==IR_LINEAR?6:1),1,0,0,GSet->IndexDegree,GSet->G2G,"####","GRIDSET","X",0,0,0,0,5,FALSE)<0) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write gridset index field (c_fstecr failed)\n",__func__);
         return(FALSE);
      }
   }
   return(TRUE);
}

/*----------------------------------------------------------------------------
 * @brief  Find a gridset within the cached list
 * @author Jean-Philippe Gauthier
 * @date   
 *    @param[in]  RefTo     Destination georeference pointer
 *    @param[in]  RefFrom   Source georeference pointer
 *    @param[in]  GSet      Existing gridset to initialise from (optional)
 *
 *    @return               Found gridset pointer, NULL on error
 */
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
      for (i=0; i<RefTo->NbSet; i++) {
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
   RefTo->Sets[i].RefTo = RefTo;
   RefTo->NbSet++;
   RefTo->Sets[i].G2G[0]=RefFrom->GRTYP[0];
   RefTo->Sets[i].G2G[1]=RefTo->GRTYP[0];
   
   // If a known grid set is passed in, use it
   if (GSet) {
      if (gset=*GSet) {
         if (gset->G2G[0]!=RefTo->Sets[i].G2G[0] || gset->G2G[1]!=RefTo->Sets[i].G2G[1]) {
            Lib_Log(APP_LIBGEOREF,APP_WARNING,"%s: Invalid grid set %c%c!=%c%c\n",__func__,gset->G2G[0],gset->G2G[1],RefTo->Sets[i].G2G[0],RefTo->Sets[i].G2G[1]);
         } else {
            RefTo->Sets[i].Index=gset->Index;
            RefTo->Sets[i].IndexDegree=gset->IndexDegree;
            RefTo->Sets[i].X=gset->X;
            RefTo->Sets[i].Y=gset->Y;
         }
      } else {
         *GSet=&RefTo->Sets[i];
      }
   }

   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: RefFrom : %p RefTo: %p\n",__func__,RefFrom,RefTo);

   return(&RefTo->Sets[i]);
}
