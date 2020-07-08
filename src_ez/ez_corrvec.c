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
wordint ez_corrvec(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, TGeoRef *gdin, TGeoRef *gdout)
{
  wordint ier;

  wordint idx_gdin;
  _gridset *gset;
  
  idx_gdin = c_find_gdin(gdin, gdout);
  
  gset = &(gdout->gset[idx_gdin]);
  if (gset->zones[AU_NORD].npts > 0)
    {
    ier = ez_corrvec_aunord(uuout,vvout,uuin,vvin, gdin, gdout);
    }
  
  if (gset->zones[AU_SUD].npts > 0)
    {
    ier = ez_corrvec_ausud(uuout,vvout,uuin,vvin, gdin, gdout);
    }
  
  if (gset->zones[POLE_NORD].npts > 0)
    {
    ier = ez_corrvec_aunord(uuout,vvout,uuin,vvin, gdin, gdout);
    }
  
  if (gset->zones[POLE_SUD].npts > 0)
    {
    ier = ez_corrvec_ausud(uuout,vvout,uuin,vvin, gdin, gdout);
    }
      
   return 0;
}
