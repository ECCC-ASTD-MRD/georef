#include <App.h>
#include "GeoRef.h"

/*----------------------------------------------------------------------------
 * @brief  Free gridset zone definitions
 * @date   
 *    @param[in]  GSet    Grid set pointer
*/
void GeoRef_SetZoneFree(TGeoSet *GSet) {

  int32_t i;

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
 * @date   
 *    @param[in]  GSet    Grid set pointer
 *    @param[in]  Zone    Zone identifier (GRID_NORTH_POLE,GRID_SOUTH_POLE)
 *    @param[in]  NbPts   Number of points to check
 *
 *    @return             Error code (0=ok)
*/
int32_t GeoRef_SetZoneDefinePole(TGeoSet *GSet,int32_t Zone,int32_t NbPts) {

   TGeoZone *zone=&GSet->zones[Zone];
   double    latpole,lonpole,xpole,ypole;
   int32_t      *tmpidx,i;
  
   // On commence par trouver les points au pole
   tmpidx = (int*)malloc(NbPts*sizeof(int));
  
   zone->npts = 0;
   if (GSet->RefFrom->GRTYP[0] == 'Z' && GSet->RefFrom->RPNHeadExt.grref[0] == 'E') {
      xpole = 0.5 * GSet->RefFrom->NX;
      ypole = (Zone==GRID_NORTH)?GSet->RefFrom->NY+0.5:0.5;
   } else {
      latpole = (Zone==GRID_NORTH)?90.0:-90.0;
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
      zone->idx = (int32_t *)malloc(zone->npts*sizeof(int));
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
 * @date   
 *    @param[in]  GSet    Grid set pointer
 *    @param[in]  Zone    Zone identifier (GRID_NORTH,GRID_SOUTH)
 *    @param[in]  NbPts   Number of points to check
 *
 *    @return             Error code (0=ok)
*/
int32_t GeoRef_SetZoneDefineThem(TGeoSet *GSet,int32_t Zone,int32_t NbPts) {

   TGeoZone *zone=&GSet->zones[Zone];
   int32_t      *tmpidx,i,jlim;
  
   tmpidx = (int*) malloc(NbPts*sizeof(int));
  
   zone->npts = 0;
   jlim = (Zone==GRID_SOUTH)?GSet->RefFrom->j1+1:GSet->RefFrom->j2-2;
   for (i=0; i<NbPts; i++) {
      if ((Zone==GRID_SOUTH)?((int)GSet->Y[i] < jlim):((int)GSet->Y[i] > jlim)) {
         tmpidx[zone->npts]=i;
         zone->npts++;
      }
   }
  
   if (zone->npts>0) {
      zone->x = (double*)malloc(zone->npts*sizeof(double));
      zone->y = (double*)malloc(zone->npts*sizeof(double));
      zone->idx = (int32_t *) malloc(zone->npts*sizeof(int));
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
 * @date   
 *    @param[in]  GSet    Grid set pointer
 *    @param[in]  Zone    Zone identifier (GRID_OUTSIDE)
 *    @param[in]  NbPts   Number of points to check
 *
 *    @return             Error code (0=ok)
*/
int32_t GeoRef_SetZoneDefineOut(TGeoSet *GSet,int32_t Zone,int32_t NbPts) {

   TGeoZone *zone=&GSet->zones[Zone];
   int32_t      *tmpidx,i,offsetleft,offsetright,ix,iy;
    
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
      zone->idx = (int32_t *) malloc(zone->npts*sizeof(int));
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
 * @date   
 *    @param[in]  GridSet   GridSet
 *
 *    @return             Error code (0=ok)
*/
int32_t GeoRef_SetZoneDefine(TGeoSet *GSet) {

   int32_t       i,npts;
   int32_t       extrap;
  
   if (!GSet || (GSet->flags & SET_ZONES)) {
      return(0);
   }

   extrap = FALSE;
   switch (GSet->RefFrom->GRTYP[0]) {
      case 'N':
      case 'S':
      case '!':
         extrap = TRUE;
         break;

      case 'L':
         if (GSet->RefTo->Extension == 0) {
               extrap = TRUE;
         }
         break;
      
      case '#':
      case 'Z':
      case 'Y':
         switch(GSet->RefFrom->RPNHeadExt.grref[0]) {
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
      GeoRef_SetZoneDefineOut(GSet,GRID_OUTSIDE,npts);
   } else {
      GeoRef_SetZoneDefinePole(GSet,GRID_NORTH_POLE,npts);
      GeoRef_SetZoneDefinePole(GSet,GRID_SOUTH_POLE,npts);
      GeoRef_SetZoneDefineThem(GSet,GRID_SOUTH,npts);
      GeoRef_SetZoneDefineThem(GSet,GRID_NORTH,npts);
   }
  
   GSet->flags |= SET_ZONES;
   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Calculates XY correspondance of destination points within the source grid
 * @date   
 *    @param[in]  GridSet   GridSet
 *
 *    @return             Error code (0=ok)
*/
int32_t GeoRef_SetCalcXY(TGeoSet *GSet) {

   int32_t size=0,mult=1;
   char *c;

   if (GSet) {
      size=GSet->RefTo->NX*GSet->RefTo->NY;

      if (!GSet->X) {

         GSet->X = (double*)calloc(size*2,sizeof(double));
         GSet->Y=&GSet->X[size];

         GeoRef_LL2XY(GSet->RefFrom,GSet->X,GSet->Y,GSet->RefTo->Lat,GSet->RefTo->Lon,size,TRUE);
      }

      if (!GSet->Index || GSet->IndexDegree!=GSet->Opt.Interp) {   
         GSet->IndexDegree=GSet->Opt.Interp;
         if (GSet->IndexDegree==IR_CONSERVATIVE || GSet->IndexDegree==IR_NORMALIZED_CONSERVATIVE) {
            mult=1024;
            if ((c=getenv("GEOREF_INDEX_SIZE_HINT"))) {
               mult=atoi(c);
            }
         } else {
            mult=GSet->IndexDegree==IR_CUBIC?10:(GSet->IndexDegree==IR_LINEAR?6:1);          
         }
         // Set the index size to the number of items without the multiplicator
         // the index size will be readjusted when in CONSERVATIVE modes otherwise, the multiplicator will be put in nj when writing
         GSet->IndexSize=size;
         GSet->Index = (float*)realloc(GSet->Index,GSet->IndexSize*mult*sizeof(float));
         GSet->Index[0]=REF_INDEX_EMPTY;
      }
   }
   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Calculates XY correspondance of destination points within the source grid (for YY grids)
 * @date   
 *    @param[in]  GridSet   GridSet
 *
 *    @return             Error code (0=ok)
*/
int32_t GeoRef_SetCalcYYXY(TGeoSet *GSet) {

   TGeoRef *yin_gdin,*yan_gdin,*yin_gdout,*yan_gdout;
   int32_t icode,k,nij,ni,nj;
   int32_t yancount_yin,yincount_yin, yancount_yan,yincount_yan;
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
   icode = GeoRef_MaskYYApply(yin_gdout,yin_gdin,&GSet->Opt,ni,nj,GSet->yin_maskout,GSet->yinlat,GSet->yinlon,yin2yin_lat,yin2yin_lon,&yincount_yin,yan2yin_lat,yan2yin_lon,&yancount_yin);
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
      icode = GeoRef_MaskYYApply(yan_gdout,yin_gdin,&GSet->Opt,ni,nj,GSet->yan_maskout,GSet->yanlat,GSet->yanlon,yin2yan_lat,yin2yan_lon,&yincount_yan,yan2yan_lat,yan2yan_lon,&yancount_yan);
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
 * @date   
 *    @param[in]  GSet      Gridset pointer
 *
*/
void GeoRef_SetFree(TGeoSet* GSet) {

   int32_t i;

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
 //   int32_t *mask_in, *mask_out;
 // float *yin_maskout,*yan_maskout;
 // float *yinlat,*yinlon,*yanlat,*yanlon;
 // float *yin2yin_lat,*yin2yin_lon,*yan2yin_lat,*yan2yin_lon;
 // float *yin2yan_lat,*yin2yan_lon,*yan2yan_lat,*yan2yan_lon;
 // float *yin2yin_x,*yin2yin_y,*yan2yin_x,*yan2yin_y;
 // float *yin2yan_x,*yin2yan_y,*yan2yan_x,*yan2yan_y;
}

/*----------------------------------------------------------------------------
 * @brief  Reads a gridset definition and index from a file
 * @date   January 2020
 *    @param[in]  RefTo      Source grid type
 *    @param[in]  RefFrom    Destination grid type
 *    @param[in]  InterpType Interpolation level (2=bilinear,3=bicubic,-1=any)
 *    @param[in]  File       FSTD file pointer
 *
 *    @return             Error code (0=ok)
*/
TGeoSet* GeoRef_SetRead(TGeoRef* RefTo,TGeoRef* RefFrom,int32_t InterpType,fst_file *File) {

   TGeoSet  *gset=NULL;
   fst_record record,crit=default_fst_record;
   char       typvar[2];
   
   if (!(gset=GeoRef_SetGet(RefTo,RefFrom,NULL))) {
      return(NULL);
   }
   
   if (File) {

      typvar[0]=RefTo->GRTYP[0];
      typvar[1]=RefFrom->GRTYP[0];

      // Rechercher et lire l'information de l'enregistrement specifie
      strncpy(crit.etiket,"GRIDSET",FST_ETIKET_LEN);
      strncpy(crit.nomvar,"####",FST_NOMVAR_LEN);
      strncpy(crit.typvar,typvar,FST_TYPVAR_LEN);
      crit.ip3=InterpType;
      if (fst24_read(File,&crit,NULL,&record)!=TRUE) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not find gridset index field (fst24_read failed)\n",__func__);
         return(NULL);
      }
      gset->IndexDegree=(TRef_InterpR)InterpType;
      gset->Index=(float*)record.data;
      record.data=NULL;
   
      strncpy(crit.nomvar,"#>>#",FST_NOMVAR_LEN);
      if (fst24_read(File,&crit,NULL,&record)!=TRUE) {
          Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not find gridset longitude field (fst24_read failed)\n",__func__);
         return(NULL);
      }
      gset->X=record.data;
      record.data=NULL;

      strncpy(crit.nomvar,"#^^#",FST_NOMVAR_LEN);
      if (fst24_read(File,&crit,NULL,&record)!=TRUE) {
          Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not find gridset longitude field (fst24_read failed)\n",__func__);
         return(NULL);
      }
      gset->Y=record.data;
      record.data=NULL;
   }
   return(gset);
}

/*----------------------------------------------------------------------------
 * @brief  Writes a gridset definition and index from a file
 * @date   January 2020
 *    @param[in]  GSet      Gridset pointer
 *    @param[in]  File      FSTD file pointer
 *
 *    @return             Error code (0=ok)
*/
int32_t GeoRef_SetWrite(TGeoSet *GSet,fst_file *File){

   int32_t size=0;
   fst_record record=default_fst_record;

   if (GeoRef_SetHasIndex(GSet)) {
      Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s:  Writing index (%ix%i)\n",__func__,size,GSet->IndexDegree==IR_CUBIC?10:(GSet->IndexDegree==IR_LINEAR?6:1));

      record.data = GSet->X;
      record.pack_bits = 64;
      record.ni   = GSet->RefTo->NX;
      record.nj   = GSet->RefTo->NY;
      record.nk   = 1;
      record.dateo= 0;
      record.deet = 0;
      record.npas = 0;
      record.ip1  = 0;
      record.ip2  = 0;
      record.ip3  = 0;
      strncpy(record.typvar,GSet->G2G,FST_TYPVAR_LEN);
      strncpy(record.nomvar,"#>>#",FST_NOMVAR_LEN);
      strncpy(record.grtyp,"X",FST_GTYP_LEN);
      strncpy(record.etiket,"GRIDSET",FST_ETIKET_LEN);
      record.ig1   = 0;
      record.ig2   = 0;
      record.ig3   = 0;
      record.ig4   = 0;
      record.data_type = FST_TYPE_REAL_IEEE;
      record.data_bits = 64;
      if (fst24_write(File,&record,FST_SKIP)<=0) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write gridset index field (fst24_write failed)\n",__func__);
         return(FALSE);
      }

      record.data = GSet->Y;
      strncpy(record.nomvar,"#^^#",FST_NOMVAR_LEN);
      if (fst24_write(File,&record,FST_SKIP)<=0) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write gridset index field (fst24_write failed)\n",__func__);
         return(FALSE);
      }

      record.data = GSet->Index;
      record.pack_bits = 32;
      record.ni   = GSet->IndexSize;
      record.nj   = GSet->IndexDegree==IR_CUBIC?10:(GSet->IndexDegree==IR_LINEAR?6:1);
      record.ip3  = GSet->IndexDegree;
      strncpy(record.nomvar,"####",FST_NOMVAR_LEN);
      record.data_bits = 32;
      if (fst24_write(File,&record,FST_SKIP)<=0) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write gridset index field (fst24_write failed)\n",__func__);
         return(FALSE);
      }
   }
   return(TRUE);
}

/*----------------------------------------------------------------------------
 * @brief  Find a gridset within the cached list
 * @date   
 *    @param[in]  RefTo     Destination georeference pointer
 *    @param[in]  RefFrom   Source georeference pointer
 *    @param[in]  Opt       Interpolation parameters
 *    @param[in]  GSet      Existing gridset to initialise from (optional)
 *
 *    @return               Found gridset pointer, NULL on error
 */
TGeoSet* GeoRef_SetGet(TGeoRef* RefTo,TGeoRef* RefFrom,TGeoOptions *Opt) {

   TGeoSet* gset=NULL;
   int32_t i;

   if (!RefTo || !RefFrom) {
      return(NULL);
   }

   pthread_mutex_lock(&RefTo->Mutex);
   if (!RefTo->Sets) {
      RefTo->Sets=(TGeoSet*)calloc(sizeof(TGeoSet),SET_MAX);
   }

   // Check for last set (most cases)
   if (RefTo->LastSet && RefTo->LastSet->RefFrom == RefFrom) {
      pthread_mutex_unlock(&RefTo->Mutex);
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
         pthread_mutex_unlock(&RefTo->Mutex);
         return(&RefTo->Sets[i]);
      }
      i++;
   }

   RefTo->NbSet++;
   pthread_mutex_unlock(&RefTo->Mutex);

   // If we get here, we have'nt found any sets, create a new one    
   RefTo->Sets[i].RefFrom = RefFrom;
   RefTo->Sets[i].RefTo = RefTo;
   RefTo->Sets[i].G2G[0]=RefFrom->GRTYP[0];
   RefTo->Sets[i].G2G[1]=RefTo->GRTYP[0];
   if (Opt) memcpy(&RefTo->Sets[i].Opt,Opt,sizeof(TGeoOptions));
   
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: RefFrom : %p RefTo: %p\n",__func__,RefFrom,RefTo);

   return(&RefTo->Sets[i]);
}