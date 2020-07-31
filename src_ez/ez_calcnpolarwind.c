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

#include "ez_funcdef.h"
#include "../src/GeoRef.h"
#include "ezscint.h"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint ez_calcnpolarwind(ftnfloat *polar_uu_in, ftnfloat *polar_vv_in, ftnfloat *uuin, ftnfloat *vvin, 
			  wordint ni, wordint nj, TGeoRef *gdin)
{
  wordint k1, k2;
  ftnfloat *polar_wd, *polar_spd,*polar_lat,*polar_lon,*polar_lat_gem, *polar_lon_gem, *polar_x, *polar_y, *polar_uu, *polar_vv;
  char grtypn[2],grtypa[2];
  ftnfloat xlat1, xlat2, xlon1, xlon2;
  wordint ig1n, ig2n, ig3n, ig4n;
  ftnfloat pi, pj, d60, dgrw;
  wordint i,j,ier;
  TGeoRef *gda, *gdps;
  ftnfloat uupole, vvpole;
  ftnfloat quatrevingtdix, zero;

  polar_uu  = (ftnfloat *) malloc(ni*sizeof(ftnfloat));
  polar_vv  = (ftnfloat *) malloc(ni*sizeof(ftnfloat));
  polar_wd  = (ftnfloat *) malloc(ni*sizeof(ftnfloat));
  polar_spd = (ftnfloat *) malloc(ni*sizeof(ftnfloat));
  polar_lat = (ftnfloat *) malloc(ni*sizeof(ftnfloat));
  polar_lon = (ftnfloat *) malloc(ni*sizeof(ftnfloat));
  polar_x   = (ftnfloat *) malloc(ni*sizeof(ftnfloat));
  polar_y   = (ftnfloat *) malloc(ni*sizeof(ftnfloat));

  for (i=0; i < ni; i++)
    {
    polar_x[i] = 1.0 * (i+1);
    polar_y[i] = 1.0 * nj;
    }
  
  c_gdllfxy_orig(gdin, polar_lat, polar_lon, polar_x, polar_y, ni);

  if (gdin->grtyp[0] == 'Z' && gdin->grref[0] == 'E')
    {
    polar_lat_gem   = (ftnfloat *) malloc(ni*sizeof(ftnfloat));
    polar_lon_gem   = (ftnfloat *) malloc(ni*sizeof(ftnfloat));
    
    for (i=0; i < ni; i++)
      {
      polar_lat_gem[i] = polar_lat[i];
      polar_lon_gem[i] = polar_lon[i];
      }
    
    f77name(cigaxg)(gdin->grref, &xlat1, &xlon1, &xlat2, &xlon2, &gdin->fst.igref[IG1], &gdin->fst.igref[IG2], &gdin->fst.igref[IG3], &gdin->fst.igref[IG4],1);
    f77name(ez_gfxyfll)(polar_lon_gem, polar_lat_gem, polar_lon, polar_lat, &ni, &xlat1, &xlon1, &xlat2, &xlon2);
    }

  grtypa[0] = 'A';
  gda = GeoRef_RPNCreate(24,12, grtypa, 0,0,0,0,0);
  c_gdwdfuv(gda, polar_spd, polar_wd,  &uuin[(nj-1)*ni], &vvin[(nj-1)*ni], polar_lat, polar_lon, ni);
  
  pi   = 0.0;
  pj   = 0.0;
  d60  = 1000.0;
  dgrw = 0.0;
  grtypn[0] = 'N';
  f77name(cxgaig)(grtypn, &ig1n, &ig2n, &ig3n, &ig4n, &pi, &pj, &d60, &dgrw,1);
  gdps = GeoRef_RPNCreate(ni, 1, grtypn, ig1n, ig2n, ig3n, ig4n, 0);
  c_gduvfwd(gdps, polar_uu, polar_vv, polar_spd,  polar_wd, polar_lat, polar_lon, ni);

  f77name(ez_calcpoleval)(&uupole, polar_uu, &ni, gdin->AX, &gdin->grtyp, &gdin->grref,1,1);
  f77name(ez_calcpoleval)(&vvpole, polar_vv, &ni, gdin->AX, &gdin->grtyp, &gdin->grref,1,1);

  quatrevingtdix = 90.0;
  zero = 0.0;
  c_gdwdfuv(gdps, polar_spd, polar_wd,  &uupole, &vvpole, &quatrevingtdix, &zero, 1);
 

  polar_lat[0] = 90.0;
  for (i=1; i < ni; i++)
    {
    polar_wd[i]  = polar_wd[0] + polar_lon[i];
    polar_spd[i] = polar_spd[0];
    polar_lat[i] = 90.0;
    }
  polar_wd[0] = polar_wd[0] + polar_lon[0];
  
  c_gduvfwd(gda, polar_uu, polar_vv, polar_spd,  polar_wd, polar_lat, polar_lon, ni);
  
  for (j=0; j < 3; j++)
    {
    for (i=0; i < ni; i++)
      {
      k1 = j * ni + i;
      k2 = (nj-3+j)  * ni + i;
      polar_uu_in[k1] = uuin[k2];
      polar_vv_in[k1] = vvin[k2];
      }
    }
  
  for (i=0; i < ni; i++)
    {
    k1 = 3 * ni+i;
    polar_uu_in[k1] = polar_uu[i];
    polar_vv_in[k1] = polar_vv[i];
    }

  free(polar_y);
  free(polar_x);
  free(polar_lat);
  free(polar_lon);
  free(polar_spd);
  free(polar_wd);
  free(polar_vv);
  free(polar_uu);

  ier = GeoRef_Free(gdps);
  return 0;
}
