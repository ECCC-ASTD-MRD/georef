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
wordint f77name(gdgxpndaxes)(PTR_AS_INT GRef, ftnfloat *ax, ftnfloat *ay)
{
   c_gdgxpndaxes((TGeoRef*)GRef, ax, ay);
   return 0;
}

wordint c_gdgxpndaxes(TGeoRef* GRef, ftnfloat *ax, ftnfloat *ay)
{
  
  wordint nix, njy;
  wordint istart, jstart;

  if (GRef->NbSub > 0)
  {
      fprintf(stderr, "<gdgxpndaxes> This operation is not supported for 'U' grids.\n");
      return -1;
  }
  
  if (!GRef->flags & EZ_AX)
  {
    fprintf(stderr, "(gdgxpndaxes) Erreur! A l'aide! Descripteurs manquants!\n");
    return -1;
  }

  switch(GRef->grtyp[0])
  {
    case 'Y':
      nix = GRef->ni * GRef->nj;
      memcpy(ax, GRef->AX, nix*sizeof(ftnfloat));
      memcpy(ay, GRef->AY, nix*sizeof(ftnfloat));
      break;
      
    default:
      nix = GRef->ni;
      njy = GRef->nj;
      if (GRef->i2 == (nix+1)) istart = 1;
      if (GRef->i2 == (nix+2)) istart = 2;
      if (GRef->i2 == (nix)) istart = 0;

      if (GRef->j2 == (njy+1)) jstart = 1;
      if (GRef->j2 == (njy+2)) jstart = 2;
      if (GRef->j2 == (njy))   jstart = 0;
      memcpy(&ax[istart],GRef->AX, nix*sizeof(ftnfloat));
      memcpy(&ay[jstart],GRef->AY, njy*sizeof(ftnfloat));
      
      if (GRef->i2 == (GRef->ni+1))
      {
        ax[0] = GRef->AX[nix-2] - 360.0; 
        ax[nix] = ax[2];
      }
      
      if (GRef->i2 == (GRef->ni+2))
      {
        ax[0] = GRef->AX[nix-1] - 360.0; 
        ax[nix] = GRef->AX[1]+360.0;
        ax[nix+1] = GRef->AX[2]+360.0;
      }
  }
  return 0;
}
