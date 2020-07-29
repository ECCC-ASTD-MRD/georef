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
wordint f77name(gdgaxes)(PTR_AS_INT GRef, ftnfloat *ax, ftnfloat *ay)
{
   c_gdgaxes((TGeoRef*)GRef, ax, ay);
   return 0;
}

wordint c_gdgaxes(TGeoRef *GRef, ftnfloat *ax, ftnfloat *ay)
{
   wordint nix, njy;
   
   switch(GRef->grtyp[0])
      {
      case 'Y':
        nix = GRef->ni * GRef->nj;
        njy = nix;
        break;

      default:
        nix = GRef->ni;
        njy = GRef->nj;
        break;
      }
   
   if (GRef->flags & EZ_AX)
      {
      memcpy(ax, GRef->AX, nix*sizeof(ftnfloat));
      memcpy(ay, GRef->AY, njy*sizeof(ftnfloat));
      }
   else
      {
      fprintf(stderr, "(gdgaxes) Erreur! A l'aide! Descripteurs manquants!\n");
      return -1;
      }
   return 0;
}
