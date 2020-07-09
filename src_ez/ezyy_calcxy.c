/*
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

#include "ezscint.h"
#include "ez_funcdef.h"
#include "../src/GeoRef.h"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint c_ezyy_calcxy(TGeoRef *gdout,TGeoRef *gdin)
{
  wordint icode,nij,i,j,k,ivalue,ni,nj,yni,ynj,yin_mgid;
  int idx_gdin;
  wordint yancount_yin,yincount_yin, yancount_yan,yincount_yan;
  wordint yyin,yyout;
  TGeoRef *yin_gdin, *yan_gdin, *yin_gdout, *yan_gdout;
  wordint yin2yin,yan2yin,yin2yan,yan2yan;
  ftnfloat *yin2yin_lat,*yin2yin_lon,*yan2yin_lat,*yan2yin_lon;
  ftnfloat *yin2yan_lat,*yin2yan_lon,*yan2yan_lat,*yan2yan_lon;
  ftnfloat *yin_fld, global_extrap_value, local_extrap_value;
  char interp_degree[32],extrap_degree[32],extrap_value[32],local_val[32];
  char global_interp_degree[32],global_extrap_degree[32];
  
 /*  need only access to either yin or Yang info for the lat and lon val */
   
  yyin=0; yyout=0;

  idx_gdin = c_find_gdin(gdin, gdout);
  /* Mettre du code au cas ou gdx_gdin == -1 */

  if (gdout->gset[idx_gdin].yyflags & XXX)
      {
      return 0;
      }

   /* Dans un premier temps on calcule la position x-y de tous les points sur la grille */

  /* To be in this routine, input source grid should be Yin-Yang */
  yyin=1;
  yin_gdin = gdin->subgrid[0];
  yan_gdin = gdin->subgrid[1];

  /* Check what the destination grid is */
  if (gdout->nsubgrids > 0)
     {
     yyout=1;
     yin_gdout = gdout->subgrid[0];
     yan_gdout = gdout->subgrid[1];
     ni = yin_gdout->ni;
     nj = yin_gdout->nj;
     }
  else
     {
     yin_gdout = gdout;
     ni = gdout->ni;
     nj = gdout->nj;
     }

  k=0;
  nij = ni*nj;

  /* Masquer les grilles YY input pour enlever overlap si OUI */
  yin2yin_lat = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
  yin2yin_lon = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
  yan2yin_lat = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
  yan2yin_lon = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
  yin2yan_lat = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
  yin2yan_lon = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
  yan2yan_lat = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
  yan2yan_lon = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
  yancount_yin=0;
  yincount_yin=0;
  if (yyout == 0)
    /* destination grid is one grid */ 
    {
    /* create mask with Yin as a priority choice and store x,y,lat,lon pos */
    gdout->gset[idx_gdin].yin_maskout = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yinlat = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yinlon = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
    icode = c_gdll(yin_gdout,gdout->gset[idx_gdin].yinlat,gdout->gset[idx_gdin].yinlon);
    icode = c_ezyymint(yin_gdout,yin_gdin,ni,nj,gdout->gset[idx_gdin].yin_maskout,gdout->gset[idx_gdin].yinlat,gdout->gset[idx_gdin].yinlon,yin2yin_lat,yin2yin_lon,&yincount_yin,yan2yin_lat,yan2yin_lon,&yancount_yin);
    /* store the lats and lons */
    gdout->gset[idx_gdin].yincount_yin = yincount_yin;
    gdout->gset[idx_gdin].yancount_yin = yancount_yin;
    gdout->gset[idx_gdin].yin2yin_lat = (ftnfloat *) malloc(yincount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yin2yin_lon = (ftnfloat *) malloc(yincount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yan2yin_lat = (ftnfloat *) malloc(yancount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yan2yin_lon = (ftnfloat *) malloc(yancount_yin*sizeof(ftnfloat));
    memcpy(gdout->gset[idx_gdin].yin2yin_lat,yin2yin_lat,yincount_yin*sizeof(ftnfloat));
    memcpy(gdout->gset[idx_gdin].yin2yin_lon,yin2yin_lon,yincount_yin*sizeof(ftnfloat));
    memcpy(gdout->gset[idx_gdin].yan2yin_lat,yan2yin_lat,yancount_yin*sizeof(ftnfloat));
    memcpy(gdout->gset[idx_gdin].yan2yin_lon,yan2yin_lon,yancount_yin*sizeof(ftnfloat));

    /* store the Xs and Ys */
    gdout->gset[idx_gdin].yin2yin_x = (ftnfloat *) malloc(yincount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yin2yin_y = (ftnfloat *) malloc(yincount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yan2yin_x = (ftnfloat *) malloc(yancount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yan2yin_y = (ftnfloat *) malloc(yancount_yin*sizeof(ftnfloat));
    icode = c_gdxyfll_orig(yin_gdin,gdout->gset[idx_gdin].yin2yin_x,gdout->gset[idx_gdin].yin2yin_y,yin2yin_lat,yin2yin_lon,yincount_yin);
    icode = c_gdxyfll_orig(yan_gdin,gdout->gset[idx_gdin].yan2yin_x,gdout->gset[idx_gdin].yan2yin_y,yan2yin_lat,yan2yin_lon,yancount_yin);
    }

  if (yyout == 1)
    { /* destination grid is a U grid*/
    /* create mask (Yin priority)with src Yin,src Yang onto dest Yin and 
                                                          store x,y pos */
    gdout->gset[idx_gdin].yin_maskout = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yinlat = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yinlon = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
    icode = c_gdll(yin_gdout,gdout->gset[idx_gdin].yinlat,gdout->gset[idx_gdin].yinlon);
    icode = c_ezyymint(yin_gdout,yin_gdin,ni,nj,gdout->gset[idx_gdin].yin_maskout,gdout->gset[idx_gdin].yinlat,gdout->gset[idx_gdin].yinlon,yin2yin_lat,yin2yin_lon,&yincount_yin,yan2yin_lat,yan2yin_lon,&yancount_yin);
    gdout->gset[idx_gdin].yincount_yin = yincount_yin;
    gdout->gset[idx_gdin].yancount_yin = yancount_yin;
    gdout->gset[idx_gdin].yin2yin_lat = (ftnfloat *) malloc(yincount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yin2yin_lon = (ftnfloat *) malloc(yincount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yan2yin_lat = (ftnfloat *) malloc(yancount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yan2yin_lon = (ftnfloat *) malloc(yancount_yin*sizeof(ftnfloat));
    memcpy(gdout->gset[idx_gdin].yin2yin_lat,yin2yin_lat,yincount_yin*sizeof(ftnfloat));
    memcpy(gdout->gset[idx_gdin].yin2yin_lon,yin2yin_lon,yincount_yin*sizeof(ftnfloat));
    memcpy(gdout->gset[idx_gdin].yan2yin_lat,yan2yin_lat,yancount_yin*sizeof(ftnfloat));
    memcpy(gdout->gset[idx_gdin].yan2yin_lon,yan2yin_lon,yancount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yin2yin_x = (ftnfloat *) malloc(yincount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yin2yin_y = (ftnfloat *) malloc(yincount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yan2yin_x = (ftnfloat *) malloc(yancount_yin*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yan2yin_y = (ftnfloat *) malloc(yancount_yin*sizeof(ftnfloat));
    icode = c_gdxyfll_orig(yin_gdin,gdout->gset[idx_gdin].yin2yin_x,gdout->gset[idx_gdin].yin2yin_y,yin2yin_lat,yin2yin_lon,yincount_yin);
    icode = c_gdxyfll_orig(yan_gdin,gdout->gset[idx_gdin].yan2yin_x,gdout->gset[idx_gdin].yan2yin_y,yan2yin_lat,yan2yin_lon,yancount_yin);
    
    /* create mask (Yin priority)with src Yin,src Yang onto dest Yang and 
                                                          store x,y pos */

    gdout->gset[idx_gdin].yan_maskout = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yanlat = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yanlon = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
    icode = c_gdll(yan_gdout,gdout->gset[idx_gdin].yanlat,gdout->gset[idx_gdin].yanlon);
    icode = c_ezyymint(yan_gdout,yin_gdin,ni,nj,gdout->gset[idx_gdin].yan_maskout,gdout->gset[idx_gdin].yanlat,gdout->gset[idx_gdin].yanlon,yin2yan_lat,yin2yan_lon,&yincount_yan,yan2yan_lat,yan2yan_lon,&yancount_yan);
    gdout->gset[idx_gdin].yincount_yan = yincount_yan;
    gdout->gset[idx_gdin].yancount_yan = yancount_yan;
    gdout->gset[idx_gdin].yin2yan_lat = (ftnfloat *) malloc(yincount_yan*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yin2yan_lon = (ftnfloat *) malloc(yincount_yan*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yan2yan_lat = (ftnfloat *) malloc(yancount_yan*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yan2yan_lon = (ftnfloat *) malloc(yancount_yan*sizeof(ftnfloat));
    memcpy(gdout->gset[idx_gdin].yin2yan_lat,yin2yan_lat,yincount_yan*sizeof(ftnfloat));
    memcpy(gdout->gset[idx_gdin].yin2yan_lon,yin2yan_lon,yincount_yan*sizeof(ftnfloat));
    memcpy(gdout->gset[idx_gdin].yan2yan_lat,yan2yan_lat,yancount_yan*sizeof(ftnfloat));
    memcpy(gdout->gset[idx_gdin].yan2yan_lon,yan2yan_lon,yancount_yan*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yin2yan_x = (ftnfloat *) malloc(yincount_yan*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yin2yan_y = (ftnfloat *) malloc(yincount_yan*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yan2yan_x = (ftnfloat *) malloc(yancount_yan*sizeof(ftnfloat));
    gdout->gset[idx_gdin].yan2yan_y = (ftnfloat *) malloc(yancount_yan*sizeof(ftnfloat));
    icode = c_gdxyfll_orig(yin_gdin,gdout->gset[idx_gdin].yin2yan_x,gdout->gset[idx_gdin].yin2yan_y,yin2yan_lat,yin2yan_lon,yincount_yan);
    icode = c_gdxyfll_orig(yan_gdin,gdout->gset[idx_gdin].yan2yan_x,gdout->gset[idx_gdin].yan2yan_y,yan2yan_lat,yan2yan_lon,yancount_yan);
    }

   free(yin2yin_lat);
   free(yin2yin_lon);
   free(yan2yin_lat);
   free(yan2yin_lon);
   free(yin2yan_lat);
   free(yin2yan_lon);
   free(yan2yan_lat);
   free(yan2yan_lon);
   gdout->gset[idx_gdin].yyflags |= XXX;
   return icode;
}

