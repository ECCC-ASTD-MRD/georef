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
wordint f77name(ezwdint)(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, PTR_AS_INT gdout, PTR_AS_INT gdin)
{
   wordint icode;

   icode = c_ezwdint(uuout, vvout, uuin, vvin, (TGeoRef*)gdout, (TGeoRef*)gdin);
   return icode;
}

wordint c_ezwdint(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, TGeoRef *gdout, TGeoRef *gdin)
{
   wordint icode;

   // gdin = iset_gdin;
   // gdout= iset_gdout;
   icode = c_ezdefset(gdout,gdin);

   if (gdin->NbSub > 0 || gdout->NbSub > 0)
      {
      icode = c_ezyywdint(uuout,vvout,uuin,vvin,gdout,gdin);
      // iset_gdin=gdin;
      // iset_gdout=gdout;
      return icode;
      }
   icode = c_ezwdint_orig(uuout, vvout, uuin, vvin,gdout,gdin);
   return icode;
}

wordint c_ezwdint_orig(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, TGeoRef *gdout, TGeoRef *gdin)
{
   wordint ier,ierc,ierc1,ierc2;
   ftnfloat *uullout = NULL;
   ftnfloat *vvllout = NULL;
   wordint npts;

   wordint cur_gdin;
   int lcl_ngdin;

   // gdin = iset_gdin;
   // gdout= iset_gdout;
   icode = c_ezdefset(gdout,gdin);

   ierc = 0;
   ierc1 = 0;
   ierc2 = 0;

   npts = gdout->ni*gdout->nj;

   groptions.vecteur = VECTEUR;

   groptions.symmetrie = SYM;
   ierc1=c_ezsint(uuout,uuin);

   groptions.symmetrie = ANTISYM;
   ierc2=c_ezsint(vvout,vvin);

   if (ierc1 == 2 || ierc2 == 2)
      {
      ierc = 2;
      }

   groptions.symmetrie = SYM;

   if (groptions.polar_correction == OUI)
     {
     ier = ez_corrvec(uuout, vvout, uuin, vvin, gdin, gdout);
     }

   uullout = (ftnfloat *) malloc(npts*sizeof(ftnfloat));
   vvllout = (ftnfloat *) malloc(npts*sizeof(ftnfloat));

   /*ezsint does not allocate lat,lon if gdin=gdout*/
   ier = ez_calclatlon(gdout);

   c_gdwdfuv(gdin, uullout, vvllout, uuout, vvout,
             gdout->lat, gdout->lon, npts);

   memcpy(uuout, uullout, npts*sizeof(ftnfloat));
   memcpy(vvout, vvllout, npts*sizeof(ftnfloat));

   groptions.vecteur = SCALAIRE;
   free(uullout);
   free(vvllout);
   return ierc;
}
