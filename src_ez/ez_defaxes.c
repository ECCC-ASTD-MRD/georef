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

void c_ezdefaxes(TGeoRef* GRef, ftnfloat *ax, ftnfloat *ay)
{
  wordint i,j;
  ftnfloat *temp, dlon;
  wordint zero, deuxnj;

  switch (GRef->grtyp[0])
    {
    case '#':
    case 'Z':
      f77name(cigaxg)(&GRef->grref,&GRef->fst.xgref[XLAT1], &GRef->fst.xgref[XLON1], &GRef->fst.xgref[XLAT2], &GRef->fst.xgref[XLON2],
		      &GRef->fst.igref[IG1], &GRef->fst.igref[IG2], &GRef->fst.igref[IG3], &GRef->fst.igref[IG4],1);

      GRef->ax = (ftnfloat *) malloc(GRef->ni*sizeof(ftnfloat));
      GRef->ay = (ftnfloat *) malloc(GRef->nj*sizeof(ftnfloat));

      memcpy(GRef->ax,ax,GRef->ni*sizeof(ftnfloat));
      memcpy(GRef->ay,ay,GRef->nj*sizeof(ftnfloat));
      ez_calcxpncof(GRef);
      ez_calcntncof(GRef);
      break;

    case 'Y':
      GRef->ax = (ftnfloat *) malloc(GRef->ni*GRef->nj*sizeof(ftnfloat));
      GRef->ay = (ftnfloat *) malloc(GRef->ni*GRef->nj*sizeof(ftnfloat));
      memcpy(GRef->ax,ax,GRef->ni*GRef->nj*sizeof(ftnfloat));
      memcpy(GRef->ay,ay,GRef->ni*GRef->nj*sizeof(ftnfloat));

      ez_calcxpncof(GRef);
      break;

    case 'G':
      GRef->grref[0] = 'L';
      GRef->fst.xgref[SWLAT] = 0.0;
      GRef->fst.xgref[SWLON] = 0.0;
      GRef->fst.xgref[DLAT] = 1.0;
      GRef->fst.xgref[DLON] = 1.0;
      f77name(cxgaig)(&GRef->grref,&GRef->fst.igref[IG1], &GRef->fst.igref[IG2], &GRef->fst.igref[IG3], &GRef->fst.igref[IG4],
		      &GRef->fst.xgref[SWLAT], &GRef->fst.xgref[SWLON], &GRef->fst.xgref[DLAT], &GRef->fst.xgref[DLON],1);

      GRef->ax = (ftnfloat *) malloc(GRef->ni*sizeof(ftnfloat));
      dlon = 360. / (ftnfloat) GRef->ni;
      for (i=0; i < GRef->ni; i++)
	      {
	      GRef->ax[i] = (ftnfloat)i * dlon;
	      }

      zero = 0;
      ez_calcxpncof(GRef);

      switch (GRef->fst.ig[IG1])
	      {
	      case GLOBAL:
	        GRef->ay = (ftnfloat *) malloc(GRef->nj*sizeof(ftnfloat));
	        temp    = (ftnfloat *) malloc(GRef->nj*sizeof(ftnfloat));
	        f77name(ez_glat)(GRef->ay,temp,&GRef->nj,&zero);
	        free(temp);
	        break;

	      case NORD:
	      case SUD:
	        deuxnj = 2 * GRef->nj;
	        GRef->ay = (ftnfloat *) malloc(deuxnj*sizeof(ftnfloat));
	        temp    = (ftnfloat *) malloc(deuxnj*sizeof(ftnfloat));
	        f77name(ez_glat)(GRef->ay,temp,&deuxnj,&zero);
	        free(temp);
	        break;
	      }


      ez_calcntncof(GRef);
      GRef->flags |= AX;
      break;

    default:
      ez_calcxpncof(GRef);
      break;
    }


  GRef->flags |= AX;

}
