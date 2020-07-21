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

#include "../src/GeoRef.h"
#include "ez_funcdef.h"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdllvval)(PTR_AS_INT GRef, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, 
                      ftnfloat *lat, ftnfloat *lon, wordint *n)
{
   wordint icode;
   
   icode = c_gdllvval((TGeoRef*)GRef, uuout,vvout, uuin, vvin, lat, lon, *n);
   return icode;
}

wordint c_gdllvval(TGeoRef *GRef, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, 
               ftnfloat *lat, ftnfloat *lon, wordint n)
{
   ftnfloat *x, *y;
   wordint ier;
   
  if (GRef->NbSub > 0 )
    {
     fprintf(stderr, "<gdllvval>: This operation is not supported for 'U' grids\n");
     return -1;
    }
  else
    {

     x = (ftnfloat *) malloc(n * sizeof(float));
     y = (ftnfloat *) malloc(n * sizeof(float));
   
     ier = c_gdxyfll_orig(GRef, x, y, lat, lon, n);
     ier = c_gdxyvval(GRef, uuout, vvout, uuin, vvin, x, y, n);
   
     free(x);
     free(y);
     return 0;
    }
}
