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
void ez_calcntncof(TGeoRef* GRef)
{
  wordint nni, nnj, gdcol, gdrow;
  
  if (GRef->flags & NEWTON)
    return;

  nni = GRef->ni;
  nnj = GRef->j2 - GRef->j1 + 1;

  if (GRef->grtyp[0] == (char)'Y') return;
  GRef->ncx = (ftnfloat *) malloc(nni*6*sizeof(ftnfloat));
  GRef->ncy = (ftnfloat *) malloc(nnj*6*sizeof(ftnfloat));
  f77name(ez_nwtncof)(GRef->ncx,GRef->ncy,
		      GRef->ax,GRef->ay,
		      &GRef->ni, &GRef->nj,
		      &GRef->i1, &GRef->i2, 
		      &GRef->j1, &GRef->j2,
		      &GRef->extension);
  
  Grille[gdrow][gdcol].flags |= NEWTON;  
}

