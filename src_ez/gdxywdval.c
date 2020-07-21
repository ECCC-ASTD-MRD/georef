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
wordint f77name(gdxywdval)(PTR_AS_INT gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, wordint *n)
{
   wordint icode;
   
   icode = c_gdxywdval((TGeoRef*)gdin, uuout, vvout, uuin, vvin, x, y, *n);
   return icode;

}
wordint c_gdxywdval(TGeoRef *gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, wordint n)
{
  wordint ier, j, icode, lni, lnj;
  TGeoRef *yin_gd, *yan_gd;
  ftnfloat *tmplat, *tmplon, *tmpy;
  ftnfloat *uuyin, *vvyin, *uuyan, *vvyan;
  ftnfloat *tmpuu, *tmpvv;
  
  tmplat = (ftnfloat *) malloc(n * sizeof(ftnfloat));
  tmplon = (ftnfloat *) malloc(n * sizeof(ftnfloat));
  tmpuu = (ftnfloat *) malloc(n * sizeof(ftnfloat));
  tmpvv = (ftnfloat *) malloc(n * sizeof(ftnfloat));
  
  if (gdin->NbSub > 0)
      {
      yin_gd=gdin->Subs[0];
      yan_gd=gdin->Subs[1];
      lni = yin_gd->ni;
      lnj = yin_gd->nj;
      tmpy = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      uuyin = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      vvyin = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      uuyan = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      vvyan = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      for (j=0; j< n; j++)
        {
          if (y[j] > yin_gd->nj)
             {
             tmpy[j]=y[j]-yin_gd->nj;
             }
          else
             {
             tmpy[j]=y[j];
             }
        }
      icode = c_gdxyvval_orig(yin_gd, tmpuu, tmpvv, uuin, vvin, x, tmpy, n);
      icode = c_gdllfxy_orig (yin_gd, tmplat, tmplon, x, tmpy, n);
      icode = c_gdwdfuv_orig (yin_gd, uuyin,vvyin,tmpuu,tmpvv,tmplat,tmplon,n);

      icode = c_gdxyvval_orig(yan_gd, tmpuu, tmpvv, &uuin[(lni*lnj)], &vvin[(lni*lnj)], x, tmpy, n);
      icode = c_gdllfxy_orig (yan_gd, tmplat, tmplon, x, tmpy, n);
      icode = c_gdwdfuv_orig (yan_gd, uuyan,vvyan,tmpuu,tmpvv,tmplat,tmplon,n);
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
      free(tmpy);

      }
  else
      {
      ier = c_gdxyvval(gdin, tmpuu, tmpvv, uuin, vvin, x, y, n);
      ier = c_gdllfxy_orig(gdin, tmplat, tmplon, x, y, n);
      ier = c_gdwdfuv(gdin, uuout, vvout, tmpuu, tmpvv, tmplat, tmplon, n);
      }

  free(tmplat);
  free(tmplon);
  free(tmpuu);
  free(tmpvv);

  return 0;
}

