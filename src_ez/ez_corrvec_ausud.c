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
wordint ez_corrvec_ausud(ftnfloat *uuout, ftnfloat *vvout,
			 ftnfloat *uuin, ftnfloat *vvin,
			 TGeoRef *gdin, TGeoRef *gdout)
{
  ftnfloat *polar_uu_in, *polar_vv_in, *corr_uus, *corr_vvs, *temp_y, ay[4];
  wordint ni, nj, i1, i2, j1, j2, degree,npts,i;
  wordint idx_gdin;
  _gridset *gset;

  idx_gdin = c_find_gdin(gdin, gdout);

  gset = &(gdout->gset[idx_gdin]);
  npts = gset->zones[AU_SUD].npts;
  ni = gdin->ni;
  nj = gdin->j2 - gdin->j1 + 1;

  i1 = 1;
  i2 = ni;
  j1 = gdin->j1 - 1;
  j2 = j1 + 3;
  degree = 3;

  polar_uu_in = (ftnfloat *) malloc(4 * ni * sizeof(ftnfloat));
  polar_vv_in = (ftnfloat *) malloc(4 * ni * sizeof(ftnfloat));
  corr_uus = (ftnfloat *) malloc(npts * sizeof(ftnfloat));
  corr_vvs = (ftnfloat *) malloc(npts * sizeof(ftnfloat));

  ez_calcspolarwind(polar_uu_in, polar_vv_in, uuin, vvin, ni, nj, gdin);

  switch (groptions.degre_interp)
    {
    case CUBIQUE:
      switch (gdin->grtyp[0])
	{
	case 'Z':
	case 'E':
	case 'G':
	  ay[0] = -90.;
	  ay[1] = gdin->AY[0];
	  ay[2] = gdin->AY[1];
	  ay[3] = gdin->AY[2];
	  f77name(ez_irgdint_3_wnnc)(corr_uus,gset->zones[AU_SUD].x,
				     gset->zones[AU_SUD].y,&npts,
				     gdin->AX, ay, polar_uu_in,
				     &ni, &j1, &j2, &gdin->extension);
	  f77name(ez_irgdint_3_wnnc)(corr_vvs,gset->zones[AU_SUD].x,
				     gset->zones[AU_SUD].y,&npts,
				     gdin->AX, ay, polar_vv_in,
				     &ni, &j1, &j2, &gdin->extension);
	  break;

	default:
	  f77name(ez_rgdint_3_wnnc)(corr_uus,gset->zones[AU_SUD].x,
				    gset->zones[AU_SUD].y,&npts,
				    polar_uu_in,&ni, &j1, &j2, &gdin->extension);
	  f77name(ez_rgdint_3_wnnc)(corr_vvs,gset->zones[AU_SUD].x,
				    gset->zones[AU_SUD].y,&npts,
				    polar_vv_in,&ni, &j1, &j2, &gdin->extension);
	  break;
	}
      break;

    case LINEAIRE:
      f77name(ez_rgdint_1_w)(corr_uus,gset->zones[AU_SUD].x,gset->zones[AU_SUD].y,&npts,polar_uu_in,&ni, &j1, &j2, &gdin->extension);
      f77name(ez_rgdint_1_w)(corr_vvs,gset->zones[AU_SUD].x,gset->zones[AU_SUD].y,&npts,polar_vv_in,&ni, &j1, &j2, &gdin->extension);
      break;

    case VOISIN:
      f77name(ez_rgdint_0)(corr_uus,gset->zones[AU_SUD].x,gset->zones[AU_SUD].y,&npts,polar_uu_in,&ni, &j1, &j2);
      f77name(ez_rgdint_0)(corr_vvs,gset->zones[AU_SUD].x,gset->zones[AU_SUD].y,&npts,polar_vv_in,&ni, &j1, &j2);
      break;
    }


  for (i=0; i < gset->zones[AU_SUD].npts; i++)
    {
    uuout[gset->zones[AU_SUD].idx[i]] = corr_uus[i];
    vvout[gset->zones[AU_SUD].idx[i]] = corr_vvs[i];
    }

  free(polar_uu_in);
  free(polar_vv_in);
  free(corr_uus);
  free(corr_vvs);

  return 0;
}
