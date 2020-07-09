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

wordint c_gdxysint(ftnfloat *zout, ftnfloat *zin, TGeoRef *gdin, ftnfloat *x, ftnfloat *y, wordint npts);
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdxysint)(ftnfloat *zout, ftnfloat *zin, PTR_AS_INT gdin, ftnfloat *x, ftnfloat *y, wordint *npts)
{
   wordint icode;

   icode = c_gdxysint(zout, zin, (TGeoRef*)gdin, x, y, *npts);
   return icode;
}

wordint c_gdxysint(ftnfloat *zout, ftnfloat *zin, TGeoRef *gdin, ftnfloat *x, ftnfloat *y, wordint npts)
{
   wordint ier;
   ftnfloat *lzin, *lxzin;
  _zone zones[NZONES];
  wordint gdout;
  _gridset tmpset;

   lzin  = NULL;
   lxzin = NULL;

   if (gdin->fst.axe_y_inverse == 1)
      {
      lzin = (ftnfloat *) malloc(gdin->ni*gdin->nj*sizeof(ftnfloat));
      memcpy(lzin, zin, gdin->ni*gdin->nj*sizeof(ftnfloat));
      f77name(permut)(lzin, &gdin->ni, &gdin->nj);
      }
   else
     {
     lzin = zin;
     }

   if (gdin->needs_expansion == OUI)
     {
     lxzin = (ftnfloat *) malloc(2*gdin->ni*gdin->nj*sizeof(ftnfloat));
     ez_xpnsrcgd(gdin, lxzin, lzin);
     }
   else
     {
     lxzin = lzin;
     }

   ier = c_gdinterp(zout, lxzin, gdin, x, y, npts);

   gdout = NMAXGRIDS-1;
/*   tmpset.gdin = gdin;
   tmpset.gdout = gdout;*/
   tmpset.x = x;
   tmpset.y = y;
/*
   Grille[gdout].ni = npts;
   Grille[gdout].nj = 1;
   Grille[gdout].grtyp = 'L';
*/

/*   if (groptions.polar_correction == OUI && groptions.vecteur != VECTEUR)
     {
     ier = ez_defzones(&tmpset);
     ier = ez_corrval(zout, lxzin, &tmpset);
     }*/


/*   if (lzin != zin && lzin != NULL)
      {
      free(lzin);
      }

   if (lxzin != lzin && lxzin != zin && lxzin != NULL)
      {
      free(lxzin);
      }*/

   return 0;
}
