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

void ez_freezones(_gridset *gdset);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdxysval)(PTR_AS_INT gdin, ftnfloat *zout, ftnfloat *zin, ftnfloat *x, ftnfloat *y, wordint *n)
{
   wordint icode;
   
   icode = c_gdxysval((TGeoRef*)gdin, zout, zin, x, y, *n);
   return icode;
}

wordint c_gdxysval(TGeoRef *gdin, ftnfloat *zout, ftnfloat *zin, ftnfloat *x, ftnfloat *y, wordint n)
{
   wordint j, icode, ni, nj;
   TGeoRef *yin_gd, *yan_gd;
   ftnfloat *zoutyin, *zoutyan;
   ftnfloat *tmpy;

  if (gdin->nsubgrids > 0)
      {
      yin_gd=gdin->subgrid[0];
      yan_gd=gdin->subgrid[1];
      ni = yin_gd->ni;
      nj = yin_gd->nj;
      tmpy = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      zoutyin = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      zoutyan = (ftnfloat *) malloc(n*sizeof(ftnfloat));
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
      icode = c_gdxysval_orig(yin_gd,zoutyin,zin,x,tmpy,n);
      icode = c_gdxysval_orig(yan_gd,zoutyan,&zin[ni*nj],x,tmpy,n);
      for (j=0; j < n; j++)
        {
        if (y[j] > yin_gd->nj)
           {
           zout[j]=zoutyan[j];
           }
        else
           {
           zout[j]=zoutyin[j];
           }
        }
      free(tmpy);
      free(zoutyan); free(zoutyin);
      return icode;
      }
   else
      {
      icode = c_gdxysval_orig(gdin,zout,zin,x,y,n);
      return icode;
      }
   
}

wordint c_gdxysval_orig(TGeoRef *gdin, ftnfloat *zout, ftnfloat *zin, ftnfloat *x, ftnfloat *y, wordint n)
{
   wordint ier;

   ier = c_gdxysint(zout, zin, gdin, x, y, n);

   return ier;
}


void ez_freezones(_gridset *gdset)
{
  wordint i;

  for (i=0; i < NZONES; i++)
    {
    if (gdset->zones[i].npts > 0)
      {
      free(gdset->zones[i].idx);
      free(gdset->zones[i].x);
      free(gdset->zones[i].y);
      gdset->zones[i].npts = 0;
      gdset->zones[i].idx  = NULL;
      gdset->zones[i].x  = NULL;
      gdset->zones[i].y = NULL;
      }
    }

}
