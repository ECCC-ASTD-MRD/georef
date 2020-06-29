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

void c_ezdefxg(TGeoRef* GRef)
{

/*   TGeoRef *gr;
  wordint gdrow_id, gdcol_id;

  c_gdkey2rowcol(gdid,  &gdrow_id,  &gdcol_id);
  gr = &Grille[gdrow_id][gdcol_id]; */

  switch (GRef->grtyp[0])
    {
    case 'A':
    case 'G':
      GRef->fst.xg[DLON]  = 360. /GRef->ni;
      GRef->fst.xg[SWLON] = 0.0;
      switch (GRef->fst.ig[IG1])
	{
	case 0:
	  GRef->fst.xg[DLAT] = 180./GRef->nj;
	  GRef->fst.xg[SWLAT] = -90. + 0.5*GRef->fst.xg[DLAT];
	  break;

	case 1:
	  GRef->fst.xg[DLAT] = 90./GRef->nj;
	  GRef->fst.xg[SWLAT] = 0.5*GRef->fst.xg[DLAT];
	  GRef->needs_expansion = OUI;
	  break;

	case 2:
	  GRef->fst.xg[DLAT] = 90./GRef->nj;
	  GRef->fst.xg[SWLAT] = -90. + 0.5*GRef->fst.xg[DLAT];
	  GRef->needs_expansion = OUI;
	  break;

	default:
	  fprintf(stderr, "<ez_gdef_fmem> 'A' grid has to be Global/North/South\n");
	  break;
	}

   switch(GRef->fst.ig[IG2])
	   {
	   case 1:
	     GRef->fst.axe_y_inverse = OUI;
	     break;

	   default:
	     break;
	   }

      break;

    case 'B':
      GRef->fst.xg[DLON] = 360. /(GRef->ni-1);
      GRef->fst.xg[SWLON] = 0.0;
      switch (GRef->fst.ig[IG1])
	      {
	      case 0:
	        GRef->fst.xg[DLAT] = 180./(GRef->nj-1);
	        GRef->fst.xg[SWLAT] = -90.;
	        break;

	      case 1:
	        GRef->fst.xg[DLAT] = 90./(GRef->nj-1);
	        GRef->fst.xg[SWLAT] = 0.;
	        GRef->needs_expansion = OUI;
	        break;

	      case 2:
	        GRef->fst.xg[DLAT] = 90./(GRef->nj-1);
	        GRef->fst.xg[SWLAT] = -90.;
	        GRef->needs_expansion = OUI;
	        break;

	      default:
	        fprintf(stderr, "<ezgdef_fmem> 'B' grid has to be Global/North/South\n");
	        break;
	      }

      switch(GRef->fst.ig[IG2])
	      {
	      case 1:
	        GRef->fst.axe_y_inverse = OUI;
	        break;

	      default:
	        break;
	      }
      break;

    case 'E':
      f77name(cigaxg)(&GRef->grtyp,&GRef->fst.xg[XLAT1],&GRef->fst.xg[XLON1],&GRef->fst.xg[XLAT2],&GRef->fst.xg[XLON2],
		      &GRef->fst.ig[IG1],&GRef->fst.ig[IG2],&GRef->fst.ig[IG3],&GRef->fst.ig[IG4],1);
      /*      GRef->fst.xg[DLAT] = 180./GRef->nj;
	      GRef->fst.xg[DLON] = 360./(GRef->ni-1);
	      GRef->fst.xg[SWLON] = 0.0;
	      GRef->fst.xg[SWLAT] = -90. + 0.5*GRef->fst.xg[DLAT];
      */
      break;

    case 'H':
    case 'Y':
    case '!':
      break;

    case '#':
    case 'Z':
      if (GRef->grref[0] == 'N') GRef->fst.hemisphere = 1;
      if (GRef->grref[0] == 'S') GRef->fst.hemisphere = 2;
      if (GRef->grref[0] == 'E')
         {
         f77name(cigaxg)(&GRef->grref,&GRef->fst.xgref[XLAT1], &GRef->fst.xgref[XLON1], &GRef->fst.xgref[XLAT2], &GRef->fst.xgref[XLON2],
            &GRef->fst.igref[IG1], &GRef->fst.igref[IG2], &GRef->fst.igref[IG3], &GRef->fst.igref[IG4],1);
         }

    break;

    case 'L':
      f77name(cigaxg)(&GRef->grtyp,&GRef->fst.xg[SWLAT], &GRef->fst.xg[SWLON], &GRef->fst.xg[DLAT], &GRef->fst.xg[DLON],
		      &GRef->fst.ig[IG1], &GRef->fst.ig[IG2], &GRef->fst.ig[IG3], &GRef->fst.ig[IG4],1);
      break;

    case 'N':
      f77name(cigaxg)(&GRef->grtyp,&GRef->fst.xg[PI], &GRef->fst.xg[PJ], &GRef->fst.xg[D60], &GRef->fst.xg[DGRW],
		      &GRef->fst.ig[IG1], &GRef->fst.ig[IG2], &GRef->fst.ig[IG3], &GRef->fst.ig[IG4],1);
      GRef->fst.hemisphere = 1;
      break;

    case 'S':
      f77name(cigaxg)(&GRef->grtyp,&GRef->fst.xg[PI], &GRef->fst.xg[PJ], &GRef->fst.xg[D60], &GRef->fst.xg[DGRW],
		      &GRef->fst.ig[IG1], &GRef->fst.ig[IG2], &GRef->fst.ig[IG3], &GRef->fst.ig[IG4],1);
      GRef->fst.hemisphere = 2;
      break;

    case 'T':
      f77name(cigaxg)(&GRef->grtyp,&GRef->fst.xg[TD60], &GRef->fst.xg[TDGRW], &GRef->fst.xg[CLAT], &GRef->fst.xg[CLON],
		      &GRef->fst.ig[IG1], &GRef->fst.ig[IG2], &GRef->fst.ig[IG3], &GRef->fst.ig[IG4],1);
      break;

    case 'X':
      fprintf(stderr,"<c_ezgdef> There is no support for grid type 'X'\n");
      return;

    default:
      fprintf(stderr,"<c_ezgdef> Grid type not supported\n");
      return;
    }

  return;

}

