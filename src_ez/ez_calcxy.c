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
wordint ez_calcxy(TGeoRef *gdin, TGeoRef *gdout)
{
   wordint coordonnee, ni_in, nj_in, ni_out, nj_out, ninj_in, ninj_out;
   wordint i,j,ier;
   wordint npts, cur_gdin, previous_val_polar_correction;
   int lcl_ngdin, idx_gdin;
   static wordint found = -1;
   static wordint ncalls = 0;

   _ygrid *ygrid;
   float *gdout_lat, *gdout_lon;

  idx_gdin = c_find_gdin(gdin, gdout);
  /* Mettre du code au cas ou gdx_gdin == -1 */

  if (gdout->gset[idx_gdin].flags & XXX)
      {
      return 0;
      }

   /* Dans un premier temps on calcule la position x-y de tous les points sur la grille */

   ni_in =  gdin->ni;
   nj_in =  gdin->nj;
   ninj_in = ni_in * nj_in;

   ni_out = gdout->ni;
   nj_out = gdout->nj;
   ninj_out = ni_out * nj_out;

   gdout->gset[idx_gdin].x = (ftnfloat *) malloc(ninj_out*sizeof(ftnfloat));
   gdout->gset[idx_gdin].y = (ftnfloat *) malloc(ninj_out*sizeof(ftnfloat));

   switch(gdin->grtyp[0])
      {
      case 'A':
      case 'B':
      case 'E':
      case 'L':
      case 'N':
      case 'S':
      case 'T':
      case '!':
        f77name(ez_ll2rgd)(gdout->gset[idx_gdin].x,
                           gdout->gset[idx_gdin].y,
                           gdout->lat, gdout->lon, &ninj_out,
                           &ni_in, &nj_in, &gdin->grtyp,
                           &gdin->fst.ig[IG1], &gdin->fst.ig[IG2],
                           &gdin->fst.ig[IG3], &gdin->fst.ig[IG4],
                           &groptions.symmetrie, gdin->AY);
        break;


      case '#':
      case 'Z':
      case 'G':
         coordonnee = RELATIF;
         f77name(ez_ll2igd)(gdout->gset[idx_gdin].x,
                            gdout->gset[idx_gdin].y,
                            gdout->lat, gdout->lon, &ninj_out,
                            &ni_in,&nj_in,&gdin->grtyp, &gdin->grref,
                            &gdin->fst.igref[IG1], &gdin->fst.igref[IG2],
                            &gdin->fst.igref[IG3], &gdin->fst.igref[IG4],
                            gdin->AX, gdin->AY,
                            &coordonnee);
         if (gdin->grtyp[0] == 'G')
            {
            if (gdin->fst.ig[IG1] == NORD)
              {
              for (j=0; j < ni_out*nj_out; j++)
                {
                gdout->gset[idx_gdin].y[j] -= nj_in;
                }
               }
            }
         break;

      case 'Y':
         previous_val_polar_correction = groptions.polar_correction;
         groptions.polar_correction = NON;
         gdout->gset[idx_gdin].ygrid.n_wts = groptions.wgt_num;
/*         fprintf(stderr, "(ez_calcxy) %d\n", groptions.wgt_num);*/
         ygrid = &(gdout->gset[idx_gdin].ygrid);
         ygrid->lat =  (float *) malloc(ninj_in*sizeof(ftnfloat));
         ygrid->lon =  (float *) malloc(ninj_in*sizeof(ftnfloat));
         gdout_lat =  (float *) malloc(ninj_out*sizeof(ftnfloat));
         gdout_lon =  (float *) malloc(ninj_out*sizeof(ftnfloat));
         ygrid->wts =  (float *) malloc(ninj_out * groptions.wgt_num*sizeof(ftnfloat));
         ygrid->idx =  (int *) malloc(ninj_out * groptions.wgt_num*sizeof(wordint));
         ygrid->mask = (int *) malloc(ninj_out*sizeof(wordint));
         ier = c_gdll(gdin, ygrid->lat, ygrid->lon);
         ier = c_gdll(gdout, gdout_lat, gdout_lon);

         if (gdin->mask == NULL)
            {
            f77name(ez_calcxy_y)(ygrid->wts, ygrid->idx,
               gdout->gset[idx_gdin].x, gdout->gset[idx_gdin].y, gdout_lat, gdout_lon, ygrid->lat, ygrid->lon, ygrid->mask,
               &ni_in, &nj_in, &ni_out, &nj_out, &(groptions.wgt_num));
            }
         else
            {
            f77name(ez_calcxy_y_m)(ygrid->wts, ygrid->idx,
               gdout->gset[idx_gdin].x, gdout->gset[idx_gdin].y, gdout_lat, gdout_lon, ygrid->mask, ygrid->lat, ygrid->lon, gdin->mask,
               &ni_in, &nj_in, &ni_out, &nj_out, &(groptions.wgt_num));
            }

         groptions.polar_correction = previous_val_polar_correction;
         free(gdout_lat);
         free(gdout_lon);
      break;

      default:
        break;
      }

   gdout->gset[idx_gdin].flags |= XXX;

   return 0;
}

