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
wordint ez_interp(ftnfloat *zout, ftnfloat *zin, TGeoRef *gdin, TGeoRef *gdout)
  {
  wordint ni_in, nj_in, ni_out, nj_out, ninj_out;
  
  wordint gdrow_in, gdrow_out, gdcol_in, gdcol_out, npts, cur_gdin, idx_gdin;
  int lcl_ngdin;

  idx_gdin = c_find_gdin(gdin, gdout);
  
  if (gdout->gset[idx_gdin].flags & XXX)
    {
    ni_out = gdout->ni;
    nj_out = gdout->nj;
    ninj_out = ni_out * nj_out;
    
    c_gdinterp(zout, zin, gdin, gdout->gset[idx_gdin].x,
      gdout->gset[idx_gdin].y, ninj_out);
    }
  return 0;   
  }

