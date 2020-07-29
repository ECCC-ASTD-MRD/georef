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
wordint ez_defzones(TGeoRef *gdin, TGeoRef *gdout)
{
  wordint i;
  wordint extrap;
  int lcl_ngdin;
  
  wordint npts, idx_gdin;

  idx_gdin = c_find_gdin(gdin, gdout);
  
  if (gdout->gset[idx_gdin].flags & ZONES)
      {
      return 0;
      }
    
  
  npts = gdout->ni * gdout->nj;
  extrap = EZ_NO_EXTRAP;
  switch (gdin->grtyp[0])
    {
    case 'N':
    case 'S':
    case 'L':
    case '!':
      extrap = EZ_EXTRAP;
      break;
      
    case '#':
    case 'Z':
    case 'Y':
      switch(gdin->grref[0])
	{
	case 'N':
	case 'S':
	case 'L':
	  extrap = EZ_EXTRAP;
	  break;
	  
	case 'E':
	  if (359.0 > (gdin->AX[gdin->ni-1] - gdin->AX[0]))
	    {
	    extrap = EZ_EXTRAP;
	    }
	  break;
	}
      break;
    }
  
  for (i=0; i < NZONES; i++)
    {
    gdout->gset[idx_gdin].zones[i].npts = 0;
    }
  
  switch (extrap)
    {
    case EZ_EXTRAP:
      ez_defzone_dehors(gdin, gdout->gset[idx_gdin].x, 
            gdout->gset[idx_gdin].y, npts, 
            &(gdout->gset[idx_gdin].zones[DEHORS]));
      break;
      
    case EZ_NO_EXTRAP:
      ez_defzone_polenord(gdin, gdout->gset[idx_gdin].x, 
            gdout->gset[idx_gdin].y, npts, 
            &(gdout->gset[idx_gdin].zones[POLE_NORD]));
      ez_defzone_polesud(gdin, gdout->gset[idx_gdin].x, 
            gdout->gset[idx_gdin].y, npts, 
            &(gdout->gset[idx_gdin].zones[POLE_SUD]));
      ez_defzone_sud(gdin, gdout->gset[idx_gdin].x, 
            gdout->gset[idx_gdin].y, npts, 
            &(gdout->gset[idx_gdin].zones[AU_SUD]));
      ez_defzone_nord(gdin, gdout->gset[idx_gdin].x, 
            gdout->gset[idx_gdin].y, npts, 
            &(gdout->gset[idx_gdin].zones[AU_NORD]));
    }
  
   gdout->gset[idx_gdin].flags |= ZONES;
  return 0;
}
