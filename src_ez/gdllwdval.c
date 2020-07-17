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

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdllwdval)(PTR_AS_INT GRef, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, 
                      ftnfloat *lat, ftnfloat *lon, wordint *n)
{
   wordint icode;
   
   icode = c_gdllwdval((TGeoRef*)GRef, uuout,vvout, uuin, vvin, lat, lon, *n);
   return icode;
}

wordint c_gdllwdval(TGeoRef *GRef, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, 
               ftnfloat *lat, ftnfloat *lon, wordint n)
{
   wordint ier,j;
   TGeoRef *yin_gd, *yan_gd;
   ftnfloat *x, *y;
   ftnfloat *uuyin, *vvyin, *uuyan, *vvyan;

   if (GRef->nsubgrids > 0)
      {
      x = (ftnfloat *) malloc(n * sizeof(float));
      y = (ftnfloat *) malloc(n * sizeof(float));
      uuyin = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      vvyin = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      uuyan = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      vvyan = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      ier = c_gdxyfll(GRef, x, y, lat, lon, n);
      ier = c_gdxyvval(GRef, uuout, vvout, uuin, vvin, x, y, n);
      yin_gd=GRef->subgrid[0];
      yan_gd=GRef->subgrid[1];
      ier = c_gdwdfuv_orig(yin_gd,uuyin,vvyin,uuout,vvout,lat,lon,n);
      ier = c_gdwdfuv_orig(yan_gd,uuyan,vvyan,uuout,vvout,lat,lon,n);
      for (j=0; j< n; j++)
        {
          if (y[j] > yin_gd->nj)
             {
             uuout[j]=uuyan[j];
             vvout[j]=vvyan[j];
             }
          else
             {
             uuout[j]=uuyin[j];
             vvout[j]=vvyin[j];
             }
        }
      free(uuyin); free(vvyin);
      free(uuyan); free(vvyan);

      }
   else
      {
      ier = c_gdllvval(GRef, uuout, vvout, uuin, vvin, lat, lon, n);
      ier = c_gdwdfuv(GRef, uuout, vvout, uuout, vvout, lat, lon, n);
      }
   return 0;
}
