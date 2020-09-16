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
#include "Vertex.h"

void GeoRef_SetZoneFree(TGridSet *GSet) {

  int i;

  for (i=0; i < NZONES; i++) {
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
   float *tmpx, *tmpy;
   float latpole, lonpole, xpole, ypole;
   int nhits, i;
   int *tmpidx;
  
   // On commence par trouver les points au pole
   tmpx =   (float *) malloc(zone->npts*sizeof(float));
   tmpy =   (float *) malloc(zone->npts*sizeof(float));
   tmpidx = (int  *) malloc(zone->npts*sizeof(int));
  
   nhits = 0;
   if (RefFrom->GRTYP[0] == 'Z' && RefFrom->RPNHead.GRREF[0] == 'E') {
      xpole = 0.5 * RefFrom->NX;
      ypole = (Zone==NORTH)?RefFrom->NY+0.5:0.5;
   } else {
      latpole = (Zone==NORTH)?90.0:-90.0;
      lonpole = 0.0;
      GeoRef_LL2XYN(RefFrom,&xpole,&ypole,&latpole,&lonpole,1,FALSE);
   }
  
   for (i=0; i < zone->npts; i++) {
      if (fabs(GSet->y[i]-ypole) < 1.0e-3) {
         tmpx[nhits] = GSet->x[i];
         tmpy[nhits] = GSet->y[i];
         tmpidx[nhits]=i;
         nhits++;
      }
   }
  
   zone->npts = nhits;
   if (nhits > 0) {
      zone->x = (float *) malloc(nhits*sizeof(float));
      zone->y = (float *) malloc(nhits*sizeof(float));
      zone->idx = (int *) malloc(nhits*sizeof(int));
      App_Log(DEBUG,"%s: Number of points at pole: %d\n",__func__,nhits);
      
      for (i=0; i < zone->npts; i++) {
         zone->x[i]   = tmpx[i];      
         zone->y[i]   = tmpy[i];     
         zone->idx[i] = tmpidx[i];
      }
   }
  
   free(tmpx);
   free(tmpy);
   free(tmpidx);

   return(0);
}

int GeoRef_SetZoneDefineThem(TGeoRef *RefFrom,TGridSet *GSet,int Zone) {

   TGeoZone *zone=&GSet->zones[Zone];
   float *tmpx, *tmpy;
   int *tmpidx;
   int nhits,i,jlim;
  
   tmpx =   (float *) malloc(zone->npts*sizeof(float));
   tmpy =   (float *) malloc(zone->npts*sizeof(float));
   tmpidx = (int  *) malloc(zone->npts*sizeof(int));
  
   nhits = 0;
   jlim = (Zone==SOUTH)?RefFrom->j1+1:RefFrom->j2-2;
   for (i=0; i < zone->npts; i++) {
      if ((Zone==SOUTH)?((int)GSet->y[i] < jlim):((int)GSet->y[i] > jlim)) {
         tmpx[nhits] = GSet->x[i];
         tmpy[nhits] = GSet->y[i];
         tmpidx[nhits]=i;
         nhits++;
      }
   }
  
   zone->npts = nhits;
   if (nhits > 0) {
      zone->x = (float *) malloc(nhits*sizeof(float));
      zone->y = (float *) malloc(nhits*sizeof(float));
      zone->idx = (int *) malloc(nhits*sizeof(int));
      App_Log(DEBUG,"%s: Number of points between pole and limit: %d\n",__func__,nhits);
    
      for (i=0; i < zone->npts; i++) {
         zone->x[i]   = tmpx[i];      
         zone->y[i]   = tmpy[i];     
         zone->idx[i] = tmpidx[i];
      }
   }
  
   free(tmpx);
   free(tmpy);
   free(tmpidx);

   return(0);
}

int GeoRef_SetZoneDefineOut(TGeoRef *RefFrom,TGridSet *GSet,int Zone) {

   TGeoZone *zone=&GSet->zones[Zone];
   float *tmpx, *tmpy;
   int *tmpidx;
   int nhits;
   int i;
   int offsetleft, offsetright, ix, iy;
    
   tmpx =   (float *) malloc(zone->npts*sizeof(float));
   tmpy =   (float *) malloc(zone->npts*sizeof(float));
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
         tmpx[nhits]  = GSet->x[i];
         tmpy[nhits]  = GSet->y[i];
         tmpidx[nhits]=i;
         nhits++;
      }
   }

   if (nhits > 0) {
      zone->npts = nhits;
      zone->x =   (float *) malloc(zone->npts*sizeof(float));
      zone->y =   (float *) malloc(zone->npts*sizeof(float));
      zone->idx = (int *) malloc(zone->npts*sizeof(int));
      App_Log(DEBUG,"%s: Number opf outside pointst: \n",__func__,offsetleft,zone->npts);
    
      for (i=0; i < zone->npts; i++) {
         zone->x[i] = tmpx[i];      
         zone->y[i] = tmpy[i];     
         zone->idx[i] = tmpidx[i];
      }
   }
  
   free(tmpx);
   free(tmpy);
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
  
   for (i=0; i < NZONES; i++) {
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

int GeoRef_CalcXY_M(TGeoRef *Ref,float *X,float *Y,float *Lat,float *Lon,int Nb) {
 
   TQTree *node;
   Vect3d  b;
   int     n,d,idx;

   for(d=0;d<Nb;d++) {
      if (Ref->QTree) {
         // If there's an index use it
         if ((node=QTree_Find(Ref->QTree,Lon[d],Lat[d])) && node->NbData) {
            
            // Loop on this nodes data payload
            for(n=0;n<node->NbData;n++) {
               idx=(intptr_t)node->Data[n].Ptr-1; // Remove false pointer increment

               if (Bary_Get(b,Ref->Wght?Ref->Wght[idx/3]:0.0,Lon[d],Lat[d],Ref->AX[Ref->Idx[idx]],Ref->AY[Ref->Idx[idx]],
                  Ref->AX[Ref->Idx[idx+1]],Ref->AY[Ref->Idx[idx+1]],Ref->AX[Ref->Idx[idx+2]],Ref->AY[Ref->Idx[idx+2]])) {
                  
                  // Return coordinate as triangle index + barycentric coefficient
                  X[d]=idx+b[0];
                  Y[d]=idx+b[1];
                  return(TRUE);
               }
            }
         }
      } else {
         // Otherwise loop on all
         for(idx=0;idx<Ref->NIdx-3;idx+=3) {
            if (Bary_Get(b,Ref->Wght?Ref->Wght[idx/3]:0.0,Lon[d],Lat[d],Ref->AX[Ref->Idx[idx]],Ref->AY[Ref->Idx[idx]],
               Ref->AX[Ref->Idx[idx+1]],Ref->AY[Ref->Idx[idx+1]],Ref->AX[Ref->Idx[idx+2]],Ref->AY[Ref->Idx[idx+2]])) {

               // Return coordinate as triangle index + barycentric coefficient
               X[d]=idx+b[0];
               Y[d]=idx+b[1];
               return(TRUE);
            }
         }
      }            
   } 
   return(TRUE);
}

int GeoRef_CalcXY_O(TGeoRef *Ref,float *X,float *Y,float *Lat,float *Lon,int Nb) {

   int     out=0;
   int     x,y,n,nd,d,dx,dy,idx,idxs[8];
   double  dists[8],xx,yy;
   Vect2d  pts[4],pt;

   for(d=0;d<Nb;d++) {

      X[d]=-1.0;
      Y[d]=-1.0;

      if ((nd=GeoRef_Nearest(Ref,Lon[d],Lat[d],idxs,dists,8))) {
     
         pt[0]=Lon[d];
         pt[1]=Lat[d];

         // Find which cell includes coordinates
         for(n=0;n<nd;n++) {
            idx=idxs[n];

            // Find within which quad
            dx=-1;dy=-1;
            if (!GeoRef_WithinCell(Ref,pt,pts,idx-Ref->NX-1,idx-1,idx,idx-Ref->NX)) {
            
               dx=0;dy=-1;
               if (!GeoRef_WithinCell(Ref,pt,pts,idx-Ref->NX,idx,idx+1,idx-Ref->NX+1)) {
                  
                  dx=-1;dy=0;
                  if (!GeoRef_WithinCell(Ref,pt,pts,idx-1,idx+Ref->NX-1,idx+Ref->NX,idx)) {
               
                     dx=0;dy=0;
                     if (!GeoRef_WithinCell(Ref,pt,pts,idx,idx+Ref->NX,idx+Ref->NX+1,idx+1)) {
                        idx=-1;
                     }
                  }
               }
            }
            
            // If found, exit loop
            if (idx!=-1) {
               break;
            }
         }

         if (idx!=-1) {
            // Map coordinates to grid
            Vertex_Map(pts,&xx,&yy,Lon[d],Lat[d]);
            X[d]=xx;Y[d]=yy;

            if (!ISNAN(X[d]) && !ISNAN(Y[d])) {
               y=idx/Ref->NX;
               x=idx-y*Ref->NX;
               Y[d]+=y+dy+1;
               X[d]+=x+dx+1; 
            } else {
               App_Log(ERROR,"%s: Invalid coordinate (NAN): ll(%f,%f) xy(%f,%f) %i\n",__func__,Lat[d],Lon[d],X[d],Y[d],idx);
               X[d]=-1,0;
               Y[d]=-1.0;
               out++;
            }
         } else {
 //           App_Log(ERROR,"%s: Point not found: %f %f %i\n",__func__,Lat[d],Lon[d],idx);
            out++;
         }

         // Si on est a l'interieur de la grille
//         if (X[d]>(Ref->X1+0.5) || Y[d]>(Ref->Y1+0.5) || X[d]<(Ref->X0-0.5) || Y[d]<(Ref->Y0-0.5)) {
//            X[d]=-1.0;
//            Y[d]=-1.0;
//            out++;
//         }
      } 
   }
   fprintf(stderr,"---- %i\n",out);
   return(out);
}

int GeoRef_SetCalcXY(TGeoRef *RefTo,TGeoRef *RefFrom) {

   TGridSet *gset=NULL;
   int coordonnee, ni_in, nj_in, ni_out, nj_out, ninj_in, ninj_out;
   int i,j,ier;
   int npts, previous_val_polar_correction;

   _ygrid *ygrid;
   float *gdout_lat, *gdout_lon;

   gset=GeoRef_SetGet(RefTo,RefFrom);

   if (gset->x) {
      return(0);
   }

   // Dans un premier temps on calcule la position x-y de tous les points sur la grille

   ni_in =  RefFrom->NX;
   nj_in =  RefFrom->NY;
   ninj_in = ni_in * nj_in;

   ni_out = RefTo->NX;
   nj_out = RefTo->NY;
   ninj_out = ni_out * nj_out;

   gset->x = (float *)malloc(ninj_out*sizeof(float));
   gset->y = (float *)malloc(ninj_out*sizeof(float));

   switch(RefFrom->GRTYP[0]) {
      case 'A':
      case 'B':
      case 'E':
      case 'L':
      case 'N':
      case 'S':
      case 'T':
      case '!':
         f77name(ez_ll2rgd)(gset->x,gset->y,RefTo->Lat,RefTo->Lon,&ninj_out,&ni_in,&nj_in,&RefFrom->GRTYP,&RefFrom->RPNHead.IG[X_IG1],&RefFrom->RPNHead.IG[X_IG2],&RefFrom->RPNHead.IG[X_IG3],&RefFrom->RPNHead.IG[X_IG4],&RefFrom->Options.Symmetric,RefFrom->AY);
         break;

      case '#':
      case 'Z':
      case 'G':
         coordonnee = RELATIVE;
         f77name(ez_ll2igd)(gset->x,gset->y,RefTo->Lat,RefTo->Lon,&ninj_out,&ni_in,&nj_in,&RefFrom->GRTYP,&RefFrom->RPNHead.GRREF,&RefFrom->RPNHead.IGREF[X_IG1],&RefFrom->RPNHead.IGREF[X_IG2],&RefFrom->RPNHead.IGREF[X_IG3],&RefFrom->RPNHead.IGREF[X_IG4],RefFrom->AX,RefFrom->AY,&coordonnee);
         if (RefFrom->GRTYP[0] == 'G') {
            if (RefFrom->RPNHead.IG[X_IG1] == NORTH) {
               for (j=0; j < ni_out*nj_out; j++) {
                  gset->y[j] -= nj_in;
               }
            }
         }
         break;

      case 'O':
         GeoRef_CalcXY_O(RefFrom,gset->x,gset->y,RefTo->Lat,RefTo->Lon,ninj_out);
         break;

      case 'M':
         GeoRef_CalcXY_M(RefFrom,gset->x,gset->y,RefTo->Lat,RefTo->Lon,ninj_out);
         break;
         
      case 'Y':
         previous_val_polar_correction = RefFrom->Options.PolarCorrect;
         RefFrom->Options.PolarCorrect = FALSE;
         gset->ygrid.n_wts = RefFrom->Options.WeightNum;
         ygrid = &(gset->ygrid);
         ygrid->lat =  (float *) malloc(ninj_in*sizeof(float));
         ygrid->lon =  (float *) malloc(ninj_in*sizeof(float));
         gdout_lat =  (float *) malloc(ninj_out*sizeof(float));
         gdout_lon =  (float *) malloc(ninj_out*sizeof(float));
         ygrid->wts =  (float *) malloc(ninj_out * RefFrom->Options.WeightNum*sizeof(float));
         ygrid->idx =  (int *) malloc(ninj_out * RefFrom->Options.WeightNum*sizeof(int));
         ygrid->mask = (int *) malloc(ninj_out*sizeof(int));
         ier = GeoRef_GetLL(RefFrom, ygrid->lat, ygrid->lon);
         ier = GeoRef_GetLL(RefTo, gdout_lat, gdout_lon);

         if (RefFrom->mask == NULL) {
            f77name(ez_calcxy_y)(ygrid->wts,ygrid->idx,gset->x,gset->y,gdout_lat,gdout_lon,ygrid->lat,ygrid->lon,ygrid->mask,&ni_in,&nj_in,&ni_out,&nj_out,&(RefFrom->Options.WeightNum));
         } else {
            f77name(ez_calcxy_y_m)(ygrid->wts,ygrid->idx,gset->x,gset->y,gdout_lat,gdout_lon,ygrid->mask,ygrid->lat,ygrid->lon,RefFrom->mask,&ni_in,&nj_in,&ni_out,&nj_out,&(RefFrom->Options.WeightNum));
         }

         RefFrom->Options.PolarCorrect = previous_val_polar_correction;
         free(gdout_lat);
         free(gdout_lon);
      break;

      default:
         App_Log(ERROR,"%s: Invalid grid type: %c\n",__func__,RefFrom->GRTYP[0]);
         break;
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
   float *yin2yin_lat,*yin2yin_lon,*yan2yin_lat,*yan2yin_lon;
   float *yin2yan_lat,*yin2yan_lon,*yan2yan_lat,*yan2yan_lon;
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
   yin2yin_lat = (float *) malloc(8*nij*sizeof(float));
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
      gset->yinlat = (float *) malloc(ni*nj*sizeof(float));
      gset->yinlon = (float *) malloc(ni*nj*sizeof(float));
      icode = GeoRef_GetLL(yin_gdout,gset->yinlat,gset->yinlon);
      icode = c_ezyymint(yin_gdout,yin_gdin, i,nj,gset->yin_maskout,gset->yinlat,gset->yinlon,yin2yin_lat,yin2yin_lon,&yincount_yin,yan2yin_lat,yan2yin_lon,&yancount_yin);
      /* store the lats and lons */
      gset->yincount_yin = yincount_yin;
      gset->yancount_yin = yancount_yin;
      gset->yin2yin_lat = (float *) malloc(yincount_yin*sizeof(float));
      gset->yin2yin_lon = (float *) malloc(yincount_yin*sizeof(float));
      gset->yan2yin_lat = (float *) malloc(yancount_yin*sizeof(float));
      gset->yan2yin_lon = (float *) malloc(yancount_yin*sizeof(float));
      memcpy(gset->yin2yin_lat,yin2yin_lat,yincount_yin*sizeof(float));
      memcpy(gset->yin2yin_lon,yin2yin_lon,yincount_yin*sizeof(float));
      memcpy(gset->yan2yin_lat,yan2yin_lat,yancount_yin*sizeof(float));
      memcpy(gset->yan2yin_lon,yan2yin_lon,yancount_yin*sizeof(float));

      /* store the Xs and Ys */
      gset->yin2yin_x = (float *) malloc(yincount_yin*sizeof(float));
      gset->yin2yin_y = (float *) malloc(yincount_yin*sizeof(float));
      gset->yan2yin_x = (float *) malloc(yancount_yin*sizeof(float));
      gset->yan2yin_y = (float *) malloc(yancount_yin*sizeof(float));
      icode = GeoRef_LL2XYN(yin_gdin,gset->yin2yin_x,gset->yin2yin_y,yin2yin_lat,yin2yin_lon,yincount_yin,FALSE);
      icode = GeoRef_LL2XYN(yan_gdin,gset->yan2yin_x,gset->yan2yin_y,yan2yin_lat,yan2yin_lon,yancount_yin,FALSE);
   }

   if (yyout == 1) { 
      /* destination grid is a U grid*/
      /* create mask (Yin priority)with src Yin,src Yang onto dest Yin and 
                                                            store x,y pos */
      gset->yin_maskout = (float *) malloc(ni*nj*sizeof(float));
      gset->yinlat = (float *) malloc(ni*nj*sizeof(float));
      gset->yinlon = (float *) malloc(ni*nj*sizeof(float));
      icode = GeoRef_GetLL(yin_gdout,gset->yinlat,gset->yinlon);
      icode = c_ezyymint(yin_gdout,yin_gdin,ni,nj,gset->yin_maskout,gset->yinlat,gset->yinlon,yin2yin_lat,yin2yin_lon,&yincount_yin,yan2yin_lat,yan2yin_lon,&yancount_yin);
      gset->yincount_yin = yincount_yin;
      gset->yancount_yin = yancount_yin;
      gset->yin2yin_lat = (float *) malloc(yincount_yin*sizeof(float));
      gset->yin2yin_lon = (float *) malloc(yincount_yin*sizeof(float));
      gset->yan2yin_lat = (float *) malloc(yancount_yin*sizeof(float));
      gset->yan2yin_lon = (float *) malloc(yancount_yin*sizeof(float));
      memcpy(gset->yin2yin_lat,yin2yin_lat,yincount_yin*sizeof(float));
      memcpy(gset->yin2yin_lon,yin2yin_lon,yincount_yin*sizeof(float));
      memcpy(gset->yan2yin_lat,yan2yin_lat,yancount_yin*sizeof(float));
      memcpy(gset->yan2yin_lon,yan2yin_lon,yancount_yin*sizeof(float));
      gset->yin2yin_x = (float *) malloc(yincount_yin*sizeof(float));
      gset->yin2yin_y = (float *) malloc(yincount_yin*sizeof(float));
      gset->yan2yin_x = (float *) malloc(yancount_yin*sizeof(float));
      gset->yan2yin_y = (float *) malloc(yancount_yin*sizeof(float));
      icode = GeoRef_LL2XYN(yin_gdin,gset->yin2yin_x,gset->yin2yin_y,yin2yin_lat,yin2yin_lon,yincount_yin,FALSE);
      icode = GeoRef_LL2XYN(yan_gdin,gset->yan2yin_x,gset->yan2yin_y,yan2yin_lat,yan2yin_lon,yancount_yin,FALSE);
      
      /* create mask (Yin priority) with src Yin,src Yang onto dest Yang and store x,y pos */

      gset->yan_maskout = (float *) malloc(ni*nj*sizeof(float));
      gset->yanlat = (float *) malloc(ni*nj*sizeof(float));
      gset->yanlon = (float *) malloc(ni*nj*sizeof(float));
      icode = GeoRef_GetLL(yan_gdout,gset->yanlat,gset->yanlon);
      icode = c_ezyymint(yan_gdout,yin_gdin,ni,nj,gset->yan_maskout,gset->yanlat,gset->yanlon,yin2yan_lat,yin2yan_lon,&yincount_yan,yan2yan_lat,yan2yan_lon,&yancount_yan);
      gset->yincount_yan = yincount_yan;
      gset->yancount_yan = yancount_yan;
      gset->yin2yan_lat = (float *) malloc(yincount_yan*sizeof(float));
      gset->yin2yan_lon = (float *) malloc(yincount_yan*sizeof(float));
      gset->yan2yan_lat = (float *) malloc(yancount_yan*sizeof(float));
      gset->yan2yan_lon = (float *) malloc(yancount_yan*sizeof(float));
      memcpy(gset->yin2yan_lat,yin2yan_lat,yincount_yan*sizeof(float));
      memcpy(gset->yin2yan_lon,yin2yan_lon,yincount_yan*sizeof(float));
      memcpy(gset->yan2yan_lat,yan2yan_lat,yancount_yan*sizeof(float));
      memcpy(gset->yan2yan_lon,yan2yan_lon,yancount_yan*sizeof(float));
      gset->yin2yan_x = (float *) malloc(yincount_yan*sizeof(float));
      gset->yin2yan_y = (float *) malloc(yincount_yan*sizeof(float));
      gset->yan2yan_x = (float *) malloc(yancount_yan*sizeof(float));
      gset->yan2yan_y = (float *) malloc(yancount_yan*sizeof(float));
      icode = GeoRef_LL2XYN(yin_gdin,gset->yin2yan_x,gset->yin2yan_y,yin2yan_lat,yin2yan_lon,yincount_yan,FALSE);
      icode = GeoRef_LL2XYN(yan_gdin,gset->yan2yan_x,gset->yan2yan_y,yan2yan_lat,yan2yan_lon,yancount_yan,FALSE);
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
      RefTo->Sets=(TGridSet*)calloc(sizeof(TGridSet),MAXSETS);
   }

   // Check for last set (most cases)
   if (RefTo->LastSet && RefTo->LastSet->RefFrom == RefFrom) {
      return(RefTo->LastSet);
   }

   if (RefTo->NbSet>=MAXSETS) {
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
