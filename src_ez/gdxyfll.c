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

#include "ezscint.h"
#include "ez_funcdef.h"
#include "../src/GeoRef.h"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdxyfll)(PTR_AS_INT GRef, ftnfloat *x, ftnfloat *y, ftnfloat *lat, ftnfloat *lon, wordint *n)
{
  return c_gdxyfll((TGeoRef*)GRef, x, y, lat, lon, *n);
}


wordint c_gdxyfll(TGeoRef *GRef, ftnfloat *x, ftnfloat *y, ftnfloat *lat, ftnfloat *lon, wordint n)
{
  wordint j, icode, maxni,maxnj ;
  TGeoRef *yin_gd, *yan_gd;
  ftnfloat *xyin, *xyan, *yyin, *yyan;

  if (GRef->nsubgrids > 0 )
    {
      yin_gd=GRef->subgrid[0];
      yan_gd=GRef->subgrid[1];
      maxni= yin_gd->ni;
      maxnj= yin_gd->nj;
      xyin = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      xyan = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      yyin = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      yyan = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      icode = c_gdxyfll_orig(yin_gd,xyin,yyin,lat,lon,n);
      icode = c_gdxyfll_orig(yan_gd,xyan,yyan,lat,lon,n);
      for (j=0; j < n; j++)
        {
        if (xyin[j] > maxni || xyin[j] < 1 || yyin[j] > maxnj || yyin[j] < 1)
         {
         /* point is no good, take from YAN eventhough it may not be good*/
         x[j]=xyan[j];
         y[j]=yyan[j]+maxnj;
         }
        else
         {
         x[j]=xyin[j];
         y[j]=yyin[j];
         }
        if (xyan[j] >= yan_gd->mymaskgridi0 &&
            xyan[j] <= yan_gd->mymaskgridi1 &&
            yyan[j] >= yan_gd->mymaskgridj0 &&
            yyan[j] <= yan_gd->mymaskgridj1)
            {
             x[j]=xyan[j];
             y[j]=yyan[j]+maxnj;
            }
        if (xyin[j] >= yin_gd->mymaskgridi0 &&
            xyin[j] <= yin_gd->mymaskgridi1 &&
            yyin[j] >= yin_gd->mymaskgridj0 &&
            yyin[j] <= yin_gd->mymaskgridj1)
            {
             x[j]=xyin[j];
             y[j]=yyin[j];
            }
        }
      free(xyin);free(xyan);free(yyin);free(yyan);
    }
  else
    {
      icode = c_gdxyfll_new(GRef,x,y,lat,lon,n);
    }
  return icode;

}

wordint c_gdxyfll_new(TGeoRef *GRef, ftnfloat *x, ftnfloat *y, ftnfloat *lat, ftnfloat *lon, wordint n)
{
  ftnfloat *tmplons;
  
  wordint j,ni_in, nj_in;
  wordint sym=groptions.symmetrie;
  
  wordint npts;
  wordint coordonnee;

	  npts = n;
	  
	  ni_in =  GRef->ni;
	  nj_in =  GRef->nj;

	  switch(GRef->grtyp[0])
	    {
	    case 'A':
	    case 'B':
	    case 'E':
	    case 'L':
	    case 'N':
	    case 'S':
	    case 'T':
	    case '!':
	      tmplons = (ftnfloat *)malloc(npts * sizeof(ftnfloat));
	      memcpy(tmplons,lon,sizeof(ftnfloat)*npts);
	      
	      f77name(ez_ll2rgd)(x, y,
				 lat, tmplons, &npts,
				 &ni_in, &nj_in, &GRef->grtyp,
				 &GRef->fst.ig[IG1], &GRef->fst.ig[IG2], &GRef->fst.ig[IG3], &GRef->fst.ig[IG4],
				 &sym, GRef->ay);
	      free(tmplons);
	      break;

	    case '#':
	    case 'Z':
	    case 'G':
	      coordonnee = RELATIF;
	      nj_in =  GRef->j2;
	      f77name(ez_ll2igd)(x, y, lat, lon, &npts,
				 &ni_in,&nj_in,&GRef->grtyp, &GRef->grref,
				 &GRef->fst.igref[IG1], &GRef->fst.igref[IG2], 
				 &GRef->fst.igref[IG3], &GRef->fst.igref[IG4],
				 GRef->ax, GRef->ay,&coordonnee);
	      if (GRef->grtyp[0] == 'G' && GRef->fst.ig[IG1] == 1) 
		 {
	          for  (j=0; j < npts; j++)
                       {
	                y[j] = y[j] - nj_in;
	               }
	         }
	      if (GRef->grtyp[0] == 'G' && GRef->fst.ig[IG2] == 1) 
		 {
	          for  (j=0; j < npts; j++)
                       {
	                y[j] = nj_in +1.0 - y[j];
	               }
	         }
      break;
      
      
    default:
      break;
    }
  
  
  return 0;
}
  
wordint c_gdxyfll_orig(TGeoRef *GRef, ftnfloat *x, ftnfloat *y, ftnfloat *lat, ftnfloat *lon, wordint n)
{
  ftnfloat *tmplons;
  
  wordint j,ni_in, nj_in;
  wordint sym=groptions.symmetrie;

  wordint npts;
  wordint coordonnee;
  npts = n;
  
  ni_in =  GRef->ni;
  nj_in =  GRef->nj;

  switch(GRef->grtyp[0])
    {
    case 'A':
    case 'B':
    case 'E':
    case 'L':
    case 'N':
    case 'S':
    case 'T':
    case '!':
      tmplons = (ftnfloat *)malloc(npts * sizeof(ftnfloat));
      memcpy(tmplons,lon,sizeof(ftnfloat)*npts);
      
      f77name(ez_ll2rgd)(x, y,
        lat, tmplons, &npts,
        &ni_in, &nj_in, &GRef->grtyp,
        &GRef->fst.ig[IG1], &GRef->fst.ig[IG2], &GRef->fst.ig[IG3], &GRef->fst.ig[IG4],
        &sym, GRef->ay);
      free(tmplons);
      break;

    case '#':
    case 'Z':
    case 'G':
      coordonnee = RELATIF;
      nj_in =  GRef->j2;
      f77name(ez_ll2igd)(x, y, lat, lon, &npts,
        &ni_in,&nj_in,&GRef->grtyp, &GRef->grref,
        &GRef->fst.igref[IG1], &GRef->fst.igref[IG2], 
        &GRef->fst.igref[IG3], &GRef->fst.igref[IG4],
        GRef->ax, GRef->ay,&coordonnee);
      if (GRef->grtyp[0] == 'G' && GRef->fst.ig[IG1] == 1) 
      {
        for  (j=0; j < npts; j++)
        {
          y[j] = y[j] - nj_in;
        }
      }
      break;
      
      
    default:
      break;
    }
  
  return 0;
}
