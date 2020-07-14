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
wordint f77name(gdllfxyz)(PTR_AS_INT GRef, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, wordint *n)
{
  return c_gdllfxyz((TGeoRef*)GRef, lat, lon, x, y, *n);
}

wordint c_gdllfxyz(TGeoRef* GRef, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, wordint n)
{
  wordint i,npts, hem, un;
  npts = n;
  
  switch(GRef->grtyp[0])
    {
    case 'A':
    case 'B':
    case 'G':
    case 'L':
    case 'N':
    case 'S':
    case 'T':
    case '!':
      c_gdllfxy_orig(GRef, lat, lon, x, y, n);
      break;
      
    case 'Y':
      fprintf(stderr, "********************************************************\n");
      fprintf(stderr, "<gdllfxy>: This operation is not supported for 'Y' grids\n");
      fprintf(stderr, "********************************************************\n");
      break;
      
    case '#':
    case 'Z':
      switch (GRef->grref[0])
  {
  case 'E':
    f77name(ez_gfllfxy)(lon,lat,x,y,&npts,&GRef->fst.xgref[XLAT1],&GRef->fst.xgref[XLON1],
            &GRef->fst.xgref[XLAT2],&GRef->fst.xgref[XLON2]);
    break;
    
  case 'S':
  case 'N':
    if (GRef->grref[0] == 'N') 
      hem = 1;
    else
      hem = 2;

    un = 1;
    f77name(ez_vllfxy)(lat,lon,x,y,&npts,&un,&GRef->fst.xgref[D60],
        &GRef->fst.xgref[DGRW], &GRef->fst.xgref[PI], &GRef->fst.xgref[PJ],&GRef->fst.hemisphere);
    break;
    
  case 'L':
    for (i=0; i < n; i++)
      {
        lat[i] = (y[i])*GRef->fst.xgref[DLAT]+ GRef->fst.xgref[SWLAT];
        lon[i] = (x[i])*GRef->fst.xgref[DLON]+ GRef->fst.xgref[SWLON];
        lon[i] = lon[i] < 0.0 ? lon[i] + 360.0 : lon[i];
      }
    break;
    
  default:
    fprintf(stderr,"<gdllfxy> Errrrrrrrrrrreur!\n");
    break;
  }
      break;
    }
  
  return 0;
  
}
