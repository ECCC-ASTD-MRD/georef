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

wordint c_gdinterp(ftnfloat *zout, ftnfloat *zin, TGeoRef *gdin, ftnfloat *x, ftnfloat *y, wordint npts)
{
   wordint lnpts;
   wordint cur_gdin, idx_gdin;
   int lcl_ngdin;
   int ier, un, j;
   TGeoRef* gdout;
   wordint old_degre_interp;
   static wordint found = -1;
   static wordint ncalls = 0;
   ftnfloat *gdst_lats, tmp, real_un, real_j;
   char option[32], value[32];

   wordint ni_in, nj_in, ni_out, nj_out, ninj_out;

   lnpts = npts;
   old_degre_interp = groptions.degre_interp;

   ni_in =  gdin->ni;
   nj_in =  gdin->nj;


   switch (groptions.degre_interp)
   {
      case 4:
      case 5:
         gdout = c_ezgetgdout();
         ni_out = gdout->ni;
         nj_out = gdout->nj;
         ninj_out = ni_out * nj_out;
         ier = c_gdcompatible_grids(gdin, gdout);
         if (ier < 0)
         {
            fprintf(stderr, "(gdinterp) input and output grids are not compatible for average computation\n");
            fprintf(stderr, "(gdinterp) interpolaton level set to linear\n");
            groptions.degre_interp = LINEAIRE;
         }
         break;

      default:
         break;
   }


   switch(gdin->grtyp[0])
   {
   case '#':
   case 'Z':
   case 'G':
      switch (groptions.degre_interp)
      {
      case VOISIN:
         f77name(ez_rgdint_0)(zout,x,y,
            &lnpts, zin, &gdin->ni,
            &gdin->j1, &gdin->j2);
         break;

      case LINEAIRE:
         switch(gdin->extension)
         {
         case 0:
            f77name(ez_irgdint_1_nw)(zout,x, y,
               &lnpts, gdin->ax, gdin->ay,
               zin,&gdin->ni, &gdin->nj);
            break;

         case 1:
         case 2:
            f77name(ez_irgdint_1_w)(zout,x, y,
               &lnpts, gdin->ax, gdin->ay,
               zin,&gdin->ni, &gdin->j1, &gdin->j2, &gdin->extension);
            break;
         }
         break;

      case CUBIQUE:
         switch(gdin->extension)
         {
         case 0:
            f77name(ez_irgdint_3_nw)(zout, x, y,
               &lnpts, gdin->ax, gdin->ay,
               gdin->ncx, gdin->ncy, zin,
               &gdin->i1, &gdin->i2,
               &gdin->j1, &gdin->j2);
            break;

         case 1:
         case 2:
            f77name(ez_irgdint_3_w)(zout, x, y,
               &lnpts, gdin->ax, gdin->ay,
               gdin->ncx, gdin->ncy, zin,
               &gdin->ni, &gdin->j1, &gdin->j2,
               &gdin->extension);
            break;
         }
         break;

      case 4:
         gdout = c_ezgetgdout();
         f77name(ez_avg)(zout, x, y, &gdout->ni, &gdout->nj,
         zin, &gdin->ni, &gdin->nj,
         &gdin->extension);
         break;

      case 5:
         gdout = c_ezgetgdout();
         gdst_lats = (float *) malloc(sizeof(ftnfloat)*lnpts);
         real_un = 1.0;
         for (j=0; j < gdout->nj; j++)
         {
            real_j = 1.0 * (j+1);
            ier = c_gdllfxy_orig(gdout, &gdst_lats[j], &tmp, &real_un, &real_j, 1);
         }
         f77name(ez_avg_sph)(zout, x, y, gdst_lats, &gdout->ni, &gdout->nj,
            zin, &gdin->ni, &gdin->nj,
            &gdin->extension);
         break;
      }
      break;

   case 'Y':
      gdout = c_ezgetgdout();
      idx_gdin = c_find_gdin(gdin, gdout);
      ni_out = gdout->ni;
      nj_out = gdout->nj;
      ninj_out = ni_out * nj_out;
      un = 1;
      strcpy(option, "cloud_interp_alg");
      ier = c_ezgetopt(option, value);
      if (ni_in > 1 && nj_in > 1 && 0 == strcmp(value, "linear"))
      {
         f77name(ez_rgdint_1_nw)(zout,x,y,
               &lnpts, zin, &gdin->ni,
               &un, &gdin->nj);
      }
      else
      {
         f77name(ez_applywgts)(zout,
            gdout->gset[idx_gdin].ygrid.wts,
            gdout->gset[idx_gdin].ygrid.idx,
            zin,
            gdout->gset[idx_gdin].ygrid.xx,
            gdout->gset[idx_gdin].ygrid.yy,
            gdout->gset[idx_gdin].ygrid.mask,
            &ni_in, &nj_in, &ni_out, &nj_out,
            &(gdout->gset[idx_gdin].ygrid.n_wts));
      }
      break;

   default:
      switch (groptions.degre_interp)
      {
      case VOISIN:
         f77name(ez_rgdint_0)(zout,x,y,
            &lnpts, zin, &gdin->ni,
            &gdin->j1, &gdin->j2);
         break;

      case LINEAIRE:
         switch(gdin->extension)
         {
         case 0:
         case 1:
            f77name(ez_rgdint_1_nw)(zout,x,y,
                  &lnpts, zin, &gdin->ni,
                  &gdin->j1, &gdin->j2);
            break;

         case 2:
            f77name(ez_rgdint_1_w)(zout,x,y,
                     &lnpts, zin, &gdin->ni,
                     &gdin->j1, &gdin->j2,
                     &gdin->extension);
         }
      break;

      case CUBIQUE:
         switch(gdin->extension)
         {
         case 0:
            f77name(ez_rgdint_3_nw)(zout, x, y,
                     &lnpts, zin, &gdin->ni,
                     &gdin->j1, &gdin->j2);
            break;

         case 1:
         case 2:
            f77name(ez_rgdint_3_w)(zout, x, y,
                     &lnpts, zin,  &gdin->ni,
                     &gdin->j1, &gdin->j2,
                     &gdin->extension);
            break;
         }
      break;

      case 4:
         gdout = c_ezgetgdout();
         f77name(ez_avg)(zout, x, y, &gdout->ni, &gdout->nj,
            zin, &gdin->ni, &gdin->nj,
            &gdin->extension);
         break;

      case 5:
         gdout = c_ezgetgdout();
         gdst_lats = (float *) malloc(sizeof(ftnfloat)*lnpts);
         real_un = 1.0;
         for (j=0; j < gdout->nj; j++)
         {
            real_j = 1.0 * (j+1);
            ier = c_gdllfxy_orig(gdout, &gdst_lats[j], &tmp, &real_un, &real_j, 1);
         }
         f77name(ez_avg_sph)(zout, x, y, gdst_lats, &gdout->ni, &gdout->nj,
            zin, &gdin->ni, &gdin->nj,
            &gdin->extension);
         free(gdst_lats);
         break;
      }
      break;
   }
   groptions.degre_interp = old_degre_interp;
   return 0;

}

/* TODO: it will never reach the 2nd switch statement. is it normal? */
int c_gdcompatible_grids(TGeoRef *gdin, TGeoRef* gdout)
{
   switch(gdout->grtyp[0])
   {
	case 'L':
	case 'A':
	case 'B':
	case 'G':
	   return 0;

	default:
	   return -1;
	}

   switch (gdin->grtyp[0])
   {
   case 'L':
   case 'A':
   case 'B':
   case 'G':
      return 0;
      break;

   case 'Z':
      if (gdin->grref[0] == 'L')
         return 0;
      else if (gdin->grref[0] == 'E')
         if (0 == c_gd_isgridrotated(gdin))
            return 0;
      else
         return -1;
      break;

   default:
      return -1;
   }

   return 0;
}

int c_gd_isgridrotated(TGeoRef *gr)
{
   int ig1, ig2, ig3, ig4;
   float xg1, xg2, xg3, xg4;

   if (gr->grref[0] == 'E')
   {
      xg1 = gr->fst.xgref[XLAT1];
      xg3 = gr->fst.xgref[XLAT2];
/*      fprintf(stderr,"gdid xg1: xg3: %d %f %f\n", gdid, xg1, xg3);*/
      if (fabs(xg1-xg3) < 0.001)
         return 0;/* non rotated*/
      else
         return 1; /*rotated*/
   }
   else
   {
      return 0;
   }
   return 0;
}
