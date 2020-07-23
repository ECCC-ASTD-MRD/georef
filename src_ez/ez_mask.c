/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "ez_funcdef.h"
#include "../src/GeoRef.h"

#define BITPOS(i) (i - ((i >> 5) << 5))
#define GETMSK(fld,i) ((fld[i >> 5]  & (1 << BITPOS(i))) >> BITPOS(i))
#define SETMSK(fld,i) (fld[i >> 5] | (fld[i] << BITPOS(i)))


int f77name(gdsetmask)(PTR_AS_INT gr, int *mask)
{
   return c_gdsetmask((TGeoRef*)gr, mask);
}

int f77name(gdgetmask)(PTR_AS_INT gr, int *mask)
{
   return c_gdgetmask((TGeoRef*)gr, mask);
}

int f77name(ezsint_m)(float *zout, float *zin)
{
   return c_ezsint_m(zout, zin);
}

int f77name(ezuvint_m)(float *uuout, float *vvout, float *uuin, float *vvin)
{
   return c_ezuvint_m(uuout, vvout, uuin, vvin);
}

int f77name(ezsint_mdm)(float *zout, int *mask_out, float *zin, int *mask_in, PTR_AS_INT gdout, PTR_AS_INT gdin)
{
   return c_ezsint_mdm(zout, mask_out, zin, mask_in, (TGeoRef*)gdout, (TGeoRef*)gdin);
}

int f77name(ezuvint_mdm)(float *uuout, float *vvout, int *mask_out, float *uuin, float *vvin, int *mask_in, PTR_AS_INT gdout, PTR_AS_INT gdin)
{
   return c_ezuvint_mdm(uuout, vvout, mask_out, uuin, vvin, mask_in, (TGeoRef*)gdout, (TGeoRef*)gdin);
}

int f77name(ezsint_mask)(int *mask_out, int *mask_in, PTR_AS_INT gdout, PTR_AS_INT gdin)
{
   return c_ezsint_mask(mask_out, mask_in, (TGeoRef*)gdout, (TGeoRef*)gdin);
}

/* ------------------------------------------------------------------------- */

int c_gdsetmask(TGeoRef *gr, int *mask)
{
   int ni, nj;
   if (gr->NbSub > 0)
   {
      fprintf(stderr, "<gdsetmask> This operation is not supported for 'U' grids.\n");
      return -1;
   }
   ni = gr->ni;
   nj = gr->nj;

   if (gr->mask != NULL)
   {
      free(gr->mask);
   }
   gr->mask = (int *)malloc(ni*nj*sizeof(int));
   memcpy(gr->mask, mask, ni*nj*sizeof(int));
   return 0;
}


/* ------------------------------------------------------------------------- */

int c_gdgetmask(TGeoRef *gr, int *mask)
{
   int ni, nj;

   if (gr->NbSub > 0)
   {
      fprintf(stderr, "<gdgetmask> This operation is not supported for 'U' grids.\n");
      return -1;
   }
   ni = gr->ni;
   nj = gr->nj;

   if (gr->mask != NULL)
   {
      memcpy(mask, gr->mask, ni*nj*sizeof(int));
      return 0;
   }
   else
   {
      mask = NULL;
      return -1;
   }
}


/* ------------------------------------------------------------------------- */

int c_ezsint_m(float *zout, float *zin)
{
   fprintf(stderr, "<ezsint_m> This operation is currently not implemented.\n");
   return 0;
} 


/* ------------------------------------------------------------------------- */

int c_ezuvint_m(float *uuout, float *vvout, float *uuin, float *vvin)
{
   fprintf(stderr, "<ezuvint_m> This operation is currently not implemented.\n");
   return 0;
}


/* ------------------------------------------------------------------------- */

int c_ezsint_mdm(float *zout, int *mask_out, float *zin, int *mask_in, TGeoRef *gdout, TGeoRef *gdin)
{
   wordint methode = 2;
   wordint ni_out, nj_out;

   // gdin = c_ezgetgdin();
   // gdout = c_ezgetgdout();

   if (gdout->NbSub > 0 || 
       gdin->NbSub > 0)
      {
       fprintf(stderr, "<ezsint_mdm> This operation is not supported for 'U' grids.\n");
       return -1;
      }

   c_ezdefset(gdout, gdin);
   ni_out = gdout->ni;
   nj_out = gdout->nj;
   c_ezsint(zout, zin, gdout, gdin);
   c_ezsint_mask(mask_out, mask_in, gdout, gdin);
   f77name(lorenzo_mask_fill)(zout, mask_out, &ni_out, &nj_out, &methode);
   return 0;

}


/* ------------------------------------------------------------------------- */

int c_ezuvint_mdm(float *uuout, float *vvout, int *mask_out, float *uuin, float *vvin, int *mask_in, TGeoRef *gdout, TGeoRef *gdin)
{
   wordint methode = 2;
   wordint ni_out, nj_out;

   // gdin = c_ezgetgdin();
   // gdout = c_ezgetgdout();

   if (gdout->NbSub > 0 || 
       gdin->NbSub > 0)
      {
       fprintf(stderr, "<ezuvint_mdm> This operation is not supported for 'U' grids.\n");
       return -1;
      }

   c_ezdefset(gdout, gdin);
   ni_out = gdout->ni;
   nj_out = gdout->nj;
   c_ezsint_mask(mask_out, mask_in, gdout, gdin);
   c_ezuvint(uuout, vvout, uuin, vvin, gdout, gdin);
   f77name(lorenzo_mask_fill)(uuout, mask_out, &ni_out, &nj_out, &methode);
   f77name(lorenzo_mask_fill)(vvout, mask_out, &ni_out, &nj_out, &methode);
   return 0;
}


/* ------------------------------------------------------------------------- */

int c_ezsint_mask(int *mask_out, int *mask_in, TGeoRef *gdout, TGeoRef *gdin)
{
   int i, j, k, npts_in, npts_out, idx_gdin;
   unsigned int bitpos;
   float *fmask_in, *fmask_out, *x, *y;
   char current_option[32], interp_degree[32];
   _ygrid *ygrid;

   // gdin = c_ezgetgdin();
   // gdout = c_ezgetgdout();

   if (gdout->NbSub > 0 || 
       gdin->NbSub > 0)
      {
       fprintf(stderr, "<ezsint_mask> This operation is not supported for 'U' grids.\n");
       return -1;
      }

   c_ezdefset(gdout, gdin);
   idx_gdin = c_find_gdin(gdin, gdout);

   npts_in  = gdin->ni *gdin->nj;
   npts_out = gdout->ni*gdout->nj;

   if (gdin->grtyp[0] == 'Y')
      {
      ygrid = &(gdout->gset[idx_gdin].ygrid);
      memcpy(mask_out, ygrid->mask, gdout->ni*gdout->nj*sizeof(int));
      }
   else
      {
      x = (float *) gdout->gset[idx_gdin].x;
      y = (float *) gdout->gset[idx_gdin].y;
      f77name(qqq_ezsint_mask)(mask_out, x, y, &gdout->ni, &gdout->nj, mask_in, &gdin->ni, &gdin->nj);
      }
   return 0;
}


/* ------------------------------------------------------------------------- */

int f77name(ezget_mask_zones)(int *mask_out, int *mask_in, PTR_AS_INT gdout, PTR_AS_INT gdin)
{
   return c_ezget_mask_zones(mask_out, mask_in, (TGeoRef*)gdout, (TGeoRef*)gdin);
}


/* ------------------------------------------------------------------------- */

int c_ezget_mask_zones(int *mask_out, int *mask_in, TGeoRef *gdout, TGeoRef *gdin)
{
   int i, j, k, npts_in, npts_out, idx_gdin;
   unsigned int bitpos;
   float *x, *y;
   char current_option[32], interp_degree[32];

   strcpy(interp_degree,"interp_degree");
   // gdin = c_ezgetgdin();
   // gdout = c_ezgetgdout();
   if (gdout->NbSub > 0 || 
       gdin->NbSub > 0)
      {
       fprintf(stderr, "<ezget_mask_zones> This operation is not supported for 'U' grids.\n");
       return -1;
      }

   c_ezdefset(gdout, gdin);
   idx_gdin = c_find_gdin(gdin, gdout);

   npts_in  = gdin->ni*gdin->nj;
   npts_out = gdout->ni*gdout->nj;

    x = (float *) gdout->gset[idx_gdin].x;
    y = (float *) gdout->gset[idx_gdin].y;

   f77name(qqq_ezget_mask_zones)(mask_out, x, y, &gdout->ni, &gdout->nj, mask_in, &gdin->ni, &gdin->nj);
    return 0;
}


/* ------------------------------------------------------------------------- */


