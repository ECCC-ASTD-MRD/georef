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
wordint f77name(gduvfwd)(PTR_AS_INT GRef, ftnfloat *uugdout, ftnfloat *vvgdout, 
                         ftnfloat *uullin, ftnfloat *vvllin, ftnfloat *latin, ftnfloat *lonin, wordint *npts)
{
   wordint icode;
   
   icode = c_gduvfwd((TGeoRef*)GRef, uugdout, vvgdout, uullin, vvllin, latin, lonin, *npts);
   return icode;
}

wordint c_gduvfwd(TGeoRef *GRef, ftnfloat *uugdout, ftnfloat *vvgdout, ftnfloat *uullin, ftnfloat *vvllin,
              ftnfloat *latin, ftnfloat *lonin, wordint npts)
  {
  wordint j, icode;

  if (GRef->nsubgrids > 0 )
    {
     fprintf(stderr, "<gduvfwd>: This operation is not supported for 'U' grids\n");
     return -1;
    }
  else
    {
      icode = c_gduvfwd_orig(GRef,uugdout,vvgdout,uullin,vvllin,latin,lonin,npts);
      return icode;
    }
}

wordint c_gduvfwd_orig(TGeoRef *GRef, ftnfloat *uugdout, ftnfloat *vvgdout, ftnfloat *uullin, ftnfloat *vvllin,
              ftnfloat *latin, ftnfloat *lonin, wordint npts)
  {
  ftnfloat *xlatingf, *xloningf, *xlatingf2, *xloningf2, *uvcart, *xyz;
  wordint ni, nj, use_sincos_cache;
  ftnfloat *lat_true,*lon_true;
  
  ni = npts;
  nj = 1;
  
  memcpy(uugdout, uullin, npts*sizeof(ftnfloat));
  memcpy(vvgdout, vvllin, npts*sizeof(ftnfloat));
  
  use_sincos_cache = NON;
  switch (GRef->grtyp[0])
    {
    case 'E':
    lat_true=(ftnfloat *)(malloc(npts*sizeof(ftnfloat)));
    lon_true=(ftnfloat *)(malloc(npts*sizeof(ftnfloat)));
    f77name(ez_gfxyfll)(lon_true,lat_true,lonin,latin,&ni,
                        &GRef->fst.xg[XLAT1],&GRef->fst.xg[XLON1],
                        &GRef->fst.xg[XLAT2],&GRef->fst.xg[XLON2]);
    
    c_ezgfwfllw(uugdout,vvgdout,latin,lonin,lat_true,lon_true,
                &ni,&nj,GRef->grtyp,
                &GRef->fst.ig[IG1],&GRef->fst.ig[IG2],
                &GRef->fst.ig[IG3],&GRef->fst.ig[IG4]);
    free(lat_true);
    free(lon_true);
    return 0;
    break;
        
        
    case '#':
    case 'Y':
    case 'Z':
    switch(GRef->grref[0])
      {
      case 'E':
      lat_true=(ftnfloat *)(malloc(npts*sizeof(ftnfloat)));
      lon_true=(ftnfloat *)(malloc(npts*sizeof(ftnfloat)));
      f77name(ez_gfxyfll)(lonin,latin,lon_true,lat_true,&ni,
                          &GRef->fst.xgref[XLAT1],&GRef->fst.xgref[XLON1],
                          &GRef->fst.xgref[XLAT2],&GRef->fst.xgref[XLON2]);
      
      c_ezgfwfllw(uugdout,vvgdout,latin,lonin,lat_true,lon_true,
                  &ni,&nj,GRef->grref,
                  &GRef->fst.igref[IG1],&GRef->fst.igref[IG2],
                  &GRef->fst.igref[IG3],&GRef->fst.igref[IG4]);
      free(lat_true);
      free(lon_true);
      return 0;
      break;
	    
      default:
      f77name(ez_gdwfllw)(uugdout,vvgdout,lonin,&ni,&nj,&GRef->grref,
                          &GRef->fst.igref[IG1],&GRef->fst.igref[IG2],
                          &GRef->fst.igref[IG3],&GRef->fst.igref[IG4], 1);
      break;
      }
        
    default:
    f77name(ez_gdwfllw)(uugdout,vvgdout,lonin,&ni,&nj,&GRef->grtyp,
                        &GRef->fst.ig[IG1],&GRef->fst.ig[IG2],
                        &GRef->fst.ig[IG3],&GRef->fst.ig[IG4], 1);
    break;
    }
   
   return 0;
}

