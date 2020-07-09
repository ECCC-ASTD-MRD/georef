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
wordint ez_corrval_ausud(ftnfloat *zout, ftnfloat *zin, TGeoRef *gdin, TGeoRef *gdout)
{
  wordint i;
  wordint npts;
  ftnfloat vpolesud;
  ftnfloat *temp, *temp_y, *vals;
  ftnfloat ay[4];
  wordint ni, nj, i1, i2, j1, j2;
  wordint un = 1;
  wordint quatre = 4;

  wordint idx_gdin;
  _gridset *gset;

  idx_gdin = c_find_gdin(gdin, gdout);
  
  gset = &(gdout->gset[idx_gdin]);
  npts = gset->zones[AU_SUD].npts;
  if (npts > 0)
    {
    ni = gdin->ni;

    i1 = 1;
    i2 = ni;
    j1 = gdin->j1 - 1;
    j2 = j1 + 3;

    temp = (ftnfloat *) malloc(4 * ni * sizeof(ftnfloat));
    vals = (ftnfloat *) malloc(npts * sizeof(ftnfloat));
    f77name(ez_calcpoleval)(&vpolesud, zin, &ni, gdin->ax,
			    &gdin->grtyp, &gdin->grref,1,1);
    f77name(ez_fillspole)(temp, zin, &ni, &gdin->j1, &gdin->j2, &vpolesud);

    switch (groptions.degre_interp)
      {
      case CUBIQUE:
	   switch (gdin->grtyp[0])
	     {
	     case 'Z':
	     case 'E':
	     case 'G':
          if  (gdin->ay[gdin->j1-1] == -90.0)
             {
                ay[0] = gdin->ay[0];
                ay[1] = gdin->ay[1];
                ay[2] = gdin->ay[2];
                ay[3] = gdin->ay[3];
                f77name(ez_irgdint_3_wnnc)(vals,gset->zones[AU_SUD].x,
                            gset->zones[AU_SUD].y,&npts,
                            gdin->ax, ay, temp,
                            &ni, &j1, &j2, &gdin->extension);
             }
    else
       {
             ay[0] = -90.0;
             ay[1] = gdin->ay[0];
             ay[2] = gdin->ay[1];
             ay[3] = gdin->ay[2];
             f77name(ez_irgdint_3_wnnc)(vals,gset->zones[AU_SUD].x,
                         gset->zones[AU_SUD].y,&npts,
                         gdin->ax, ay, temp,
                         &ni, &j1, &j2, &gdin->extension);

       }
	       break;

	     default:
	       f77name(ez_rgdint_3_wnnc)(vals,gset->zones[AU_SUD].x,
				         gset->zones[AU_SUD].y,&npts,
				         temp,&ni, &j1, &j2, &gdin->extension);
	       break;
	     }
	break;

      case LINEAIRE:
	   temp_y = (ftnfloat *) malloc(npts*sizeof(ftnfloat));
/*	   for (i=0; i < npts; i++)
	     {
	     temp_y[i] = gset->zones[AU_SUD].y[i] - (1.0*j1);
	     }
	   f77name(ez_rgdint_1_nw)(vals,gset->zones[AU_SUD].x,temp_y,&npts,temp,&ni, &un, &quatre);*/
   	   f77name(ez_rgdint_1_w)(vals,gset->zones[AU_SUD].x,gset->zones[AU_SUD].y,&npts,temp,&ni, &j1, &j2,
&gdin->extension);
	   free(temp_y);
	   break;

      case VOISIN:
  	   temp_y = (ftnfloat *) malloc(npts*sizeof(ftnfloat));
/*	   for (i=0; i < npts; i++)
	     {
	     temp_y[i] = gset->zones[AU_SUD].y[i] - (1.0*j1);
	     }*/
	   f77name(ez_rgdint_0)(vals,gset->zones[AU_SUD].x,gset->zones[AU_SUD].y,&npts,temp,&ni, &j1, &j2);
	   free(temp_y);
	   break;
      }

    for (i=0; i < gset->zones[AU_SUD].npts; i++)
      {
      zout[gset->zones[AU_SUD].idx[i]] = vals[i];
      }

    free(vals);
    free(temp);
    }
  return 0;
}

