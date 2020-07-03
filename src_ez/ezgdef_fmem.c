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

PTR_AS_INT f77name(ezgdef_fmem)(wordint* ni, wordint* nj, char* grtyp, char* grref,
   wordint* ig1, wordint* ig2, wordint* ig3, wordint* ig4,
   ftnfloat* ax, ftnfloat* ay, F2Cl lengrtyp, F2Cl lengrref)
{
  PTR_AS_INT icode;
  char lgrtyp[2];
  char lgrref[2];

  lgrtyp[0] = grtyp[0];
  lgrtyp[1] = '\0';

  lgrref[0] = grref[0];
  lgrref[1] = '\0';

  icode = (PTR_AS_INT) c_ezgdef_fmem(*ni, *nj, lgrtyp, lgrref, *ig1, *ig2, *ig3, *ig4, ax, ay);
  return icode;
}

//! Insert a grid entry into the list of grids managed by ezscint.  Can be used
//! with regular and irregular ('Y', 'Z') grids, although it is not very useful
//! for regular grids.
//! @param ni Horizontal size of the grid
//! @param nj
//! @param grtyp Grid type ('A', 'B', 'E', 'G', 'L', 'N', 'S','Y', 'Z', '#', '!')
//! @param grref Reference grid type ('E', 'G', 'L', 'N', 'S')
//! @param ig1 ig1 value associated to the reference grid
//! @param ig2 ig2 value associated to the reference grid
//! @param ig3 ig3 value associated to the reference grid
//! @param ig4 ig4 value associated to the reference grid
//! @param ax Positional axis mapped to the '>>' record
//! @param ay Positional axis mapped to the '^^' record
//!
//! If the grid type corresponds to a regular grid type (eg. 'A', 'G', 'N', etc.),
//! then the parameters IG1 through IG4 are taken from an ordinary data record
//! and grref, ax and ay are not used.
//!
//! If grtyp == 'Z' or '#', the dimensions of ax=ni and ay=nj.
//! If grtyp == 'Y', the dimensions of ax=ay=ni*nj. 
TGeoRef* c_ezgdef_fmem(wordint ni, wordint nj, char* grtyp, char* grref,
   wordint ig1, wordint ig2, wordint ig3, wordint ig4, ftnfloat* ax, ftnfloat* ay)
{
   TGeoRef* GRef;

   if (grtyp[0] == '#' || grtyp[0] == 'Y' || grtyp[0] == 'Z' || grtyp[0] == 'G') {
      GRef = c_ezidentify_irreg_grid(ni, nj, grtyp, grref, ig1, ig2, ig3, ig4, ax, ay);
      c_ezdefxg(GRef);
      c_ezdefaxes(GRef, ax, ay);
   } else {
      GRef = c_ezidentify_reg_grid(ni, nj, grtyp, ig1, ig2, ig3, ig4);
      c_ezdefxg(GRef);
   }

   ez_calcxpncof(GRef);

/*    GRef->Grid[0]=grtyp[0];
   GRef->Grid[1]=grtyp[1];
   GRef->Project=GeoRef_RPNProject;
   GRef->UnProject=GeoRef_RPNUnProject;
   GRef->Value=(TGeoRef_Value*)GeoRef_RPNValue;
   GRef->Distance=GeoRef_RPNDistance;
   GRef->Height=NULL; */

   if (groptions.verbose > 0) {
/*       printf("Gdid = \n", gdid); */
      printf("Grille[].grtyp = '%c'\n", GRef->grtyp[0]);
      printf("Grille[].ni    = %d\n",   GRef->ni);
      printf("Grille[].nj    = %d\n",   GRef->nj);
      printf("Grille[].ig[IG1]   = %d\n",   GRef->fst.ig[IG1]);
      printf("Grille[].ig[IG2]   = %d\n",   GRef->fst.ig[IG2]);
      printf("Grille[].ig[IG3]   = %d\n",   GRef->fst.ig[IG3]);
      printf("Grille[].ig[IG4]   = %d\n",   GRef->fst.ig[IG4]);
      printf("Grille[].grref = '%c'\n", GRef->grref[0]);
      printf("Grille[].igref[IG1]= %d\n",   GRef->fst.igref[IG1]);
      printf("Grille[].igref[IG2]= %d\n",   GRef->fst.igref[IG2]);
      printf("Grille[].igref[IG3]= %d\n",   GRef->fst.igref[IG3]);
      printf("Grille[].igref[IG4]= %d\n",   GRef->fst.igref[IG4]);
   }

   return GRef;
}
