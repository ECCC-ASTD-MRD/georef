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

TGeoRef *c_ezidentify_reg_grid(wordint ni, wordint nj, char *grtyp,
                               wordint ig1, wordint ig2, wordint ig3, wordint ig4)
{
   TGeoRef *GRef, *fref;

   GRef = GeoRef_New();

   GRef->grtyp[0] = grtyp[0];
   GRef->grtyp[1] = '\0';
   GRef->grref[0] = (char)0;
   GRef->ni = ni;
   GRef->nj = nj;
   GRef->fst.ig[IG1] = ig1;
   GRef->fst.ig[IG2] = ig2;
   GRef->fst.ig[IG3] = ig3;
   GRef->fst.ig[IG4] = ig4;
   GRef->fst.igref[IG1] = 0;
   GRef->fst.igref[IG2] = 0;
   GRef->fst.igref[IG3] = 0;
   GRef->fst.igref[IG4] = 0;
   GRef->IG1_JP = ig1;
   GRef->IG2_JP = ig2;
   GRef->IG3_JP = ig3;
   GRef->IG4_JP = ig4;

   GeoRef_Size(GRef, 0, 0, ni - 1, nj - 1, 0);

   // This georef already exists
   if (fref = GeoRef_Find(GRef))
   {
      free(GRef);
      GeoRef_Incr(fref);
      return fref;
   }

   // This is a new georef
   GeoRef_Add(GRef);

   return GRef;
}

TGeoRef *c_ezidentify_irreg_grid(
    wordint ni, wordint nj, char *grtyp, char *grref,
    wordint ig1, wordint ig2, wordint ig3, wordint ig4,
    ftnfloat *ax, ftnfloat *ay)
{
   TGeoRef *GRef, *fref;

   GRef = GeoRef_New();

   GRef->grtyp[0] = grtyp[0];
   GRef->grtyp[1] = '\0';
   GRef->grref[0] = grref[0];
   GRef->grref[1] = '\0';
   GRef->ni = ni;
   GRef->nj = nj;
   GRef->fst.ip1 = 0;
   GRef->fst.ip2 = 0;
   GRef->fst.ip3 = 0;
   GRef->fst.igref[IG1] = ig1;
   GRef->fst.igref[IG2] = ig2;
   GRef->fst.igref[IG3] = ig3;
   GRef->fst.igref[IG4] = ig4;
   GRef->fst.ig[IG1] = ig1;
   GRef->fst.ig[IG2] = ig2;
   GRef->fst.ig[IG3] = ig3;
   GRef->fst.ig[IG4] = ig4;
   GRef->fst.xg[IG1] = 0.0;
   GRef->fst.xg[IG2] = 0.0;
   GRef->fst.xg[IG3] = 0.0;
   GRef->fst.xg[IG4] = 0.0;
   GRef->NbSub = 0;
   strcpy(GRef->fst.nomvarx, "    ");
   strcpy(GRef->fst.nomvary, "    ");
   strcpy(GRef->fst.etiketx, "            ");
   strcpy(GRef->fst.etikety, "            ");
   strcpy(GRef->fst.typvarx, "  ");
   strcpy(GRef->fst.typvary, "  ");
   GRef->fst.deet = 0;
   GRef->fst.npas = 0;
   GRef->fst.nbits = 0;
   GRef->fst.date = 0;
   GRef->i1 = 1;
   GRef->i2 = ni;
   GRef->j1 = 1;
   GRef->j2 = nj;
   GRef->idx_last_gdin = -1;

   switch (grtyp[0])
   {
   case '#':
      GRef->AX = ax;
      GRef->AY = ay;
      break;

   case 'Y':
      GRef->AX = ax;
      GRef->AY = ay;
      break;

   case 'Z':
      f77name(cigaxg)(&(GRef->grref),
                      &GRef->fst.xgref[XLAT1], &GRef->fst.xgref[XLON1], &GRef->fst.xgref[XLAT2], &GRef->fst.xgref[XLON2],
                      &GRef->fst.igref[IG1], &GRef->fst.igref[IG2], &GRef->fst.igref[IG3], &GRef->fst.igref[IG4], 1);
      GRef->AX = ax;
      GRef->AY = ay;
      break;

   case 'G':
      break;

   default:
      fprintf(stderr, "c_ezidentify_irreg_grid : undefined grid type : %c\n", grtyp[0]);
      exit(13);
   }

   GeoRef_Size(GRef, 0, 0, ni - 1, nj - 1, 0);

   // This georef already exists
   if (fref = GeoRef_Find(GRef))
   {
      free(GRef);
      GeoRef_Incr(fref);
      return fref;
   }

   // This is a new georef
   GeoRef_Add(GRef);

   return GRef;
}
