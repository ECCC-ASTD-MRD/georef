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

void reallocate_gridset_table(TGeoRef* gdid);
void   allocate_gridset_table(TGeoRef* gdid);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezdefset)(PTR_AS_INT gdout, PTR_AS_INT gdin)
{
   wordint icode;

   icode = c_ezdefset((TGeoRef*)gdout, (TGeoRef*)gdin);
   return icode;
}

wordint c_ezdefset(TGeoRef* gdout, TGeoRef* gdin)
{
  /* d'abord trouver si l'ensemble est deja defini */

   wordint i;
   wordint gdrow_in, gdrow_out, gdcol_in, gdcol_out, npts, cur_gdin;
   int idx_gdin;
   static wordint found = -1;
   wordint nsets = 0;

   if (gdout == NULL)
   {
      if (gdin != NULL)
      {
         gdout = gdin;
      }
   }

   if (gdout == gdin)
   {
      iset_gdin = gdin;
      iset_gdout = gdout;
   }

   nsets = gdout->n_gdin;

   if (nsets == 0)
   {
      allocate_gridset_table(gdout);
      gdout->log_chunk_gdin = cur_log_chunk;
   }

   if (nsets >= primes[MAX_LOG_CHUNK])
   {
      for (i=0; i < nsets; i++)
      {
         if (gdout->gset[i].gdin != NULL)
         {
            c_ezfreegridset(gdin, i);
         }
      }
      nsets = 0;
      allocate_gridset_table(gdout);
      gdout->log_chunk_gdin = cur_log_chunk;
   }


   found = -1;

   idx_gdin = gdin % primes[gdout->log_chunk_gdin];
   if (gdout->gset[idx_gdin].gdin == gdin)
   {
      found = 1;
      iset_gdin = gdin;
      iset_gdout = gdout;
      return 1;
   }

   i = idx_gdin;
   while ((found == -1) && (i != idx_gdin-1) && (gdout->gset[i].gdin != -1))
   {
      if (gdin == gdout->gset[i].gdin)
      {
         found = i;
         iset_gdin = gdin;
         iset_gdout = gdout;
         gdout->idx_last_gdin = gdout->gset[i].gdin;
         return 1;
      }
      else
      {
         i++;
         if (0 == (i % (primes[gdout->log_chunk_gdin])))
         {
            i = 0;
         }
      }
   }

   /* si on se rend jusqu'ici alors c'est que le set n'a pas ete trouve */

   /* On initialise le vecteur de gdout pour lequel gdin est utilise en entree
      Ceci sera utile si le vecteur de grilles deborde */

   gdout->gset[i].gdin = gdin;
   cur_gdin = gdin;
   gdout->n_gdin++;

   npts = gdout->ni * gdout->nj;
   gdout->gset[i].x = malloc (sizeof(ftnfloat)*npts);
   gdout->gset[i].y = malloc (sizeof(ftnfloat)*npts);
   gdout->gset[i].use_sincos_cache = NON;

   if (gdout->n_gdin >= (primes[gdout->log_chunk_gdin]/2))
   {
      reallocate_gridset_table(gdout);
   }

   if (gdin->n_gdin_for == 0)
   {
      gdin->gdin_for = malloc(CHUNK *sizeof(TGeoRef*));
      for (i=0; i < CHUNK; i++)
         {
            gdin->gdin_for[i] = NULL;
         }
      gdin->gdin_for[0] = gdout;
      gdin->n_gdin_for++;
   }
   else
   {
      if (0 == (gdin->n_gdin_for % CHUNK))
      {
         gdin->gdin_for = (int *) realloc(gdin->gdin_for, (gdin->n_gdin_for+CHUNK)*sizeof(int));
      }
      gdin->gdin_for[gdin->n_gdin_for] = gdout;
      gdin->n_gdin_for++;
   }

   if (groptions.verbose > 0)
   {
      printf("gdin : %d gdout: %d\n", gdin, gdout);
      printf("cur_gdin                           = %03d\n", cur_gdin);
      printf("n_gdin                             = %03d\n", gdout->n_gdin);
/*       printf("Grille[%03d][%03d].gset[%03d].gdin = %d\n", gdrow_out, gdcol_out, cur_gdin, gdin); */
   }
   iset_gdin = gdin;
   iset_gdout = gdout;
   return 1;
}

void reallocate_gridset_table(TGeoRef* gr)
{
   wordint i;
   wordint newIndex, inserted, curIndex;
   int oldChunkSize, newChunkSize;
   int cur_chunk;
   static wordint found = -1;
   _gridset *gset, *newTable;

   cur_chunk = gr->log_chunk_gdin;
   oldChunkSize = primes[cur_chunk];
   newChunkSize = primes[cur_chunk + 1];
   newTable = (_gridset *) calloc(sizeof(_gridset), newChunkSize);

   for (i=0 ; i < newChunkSize; i++)
   {
      newTable[i].gdin = NULL;
   }

   for (i=0 ; i < oldChunkSize; i++)
   {
      if (gr->gset[i].gdin != NULL)
      {
         newIndex = gr->gset[i].gdin % newChunkSize;
         if (newTable[newIndex].gdin == NULL)
         {
            memcpy(&(newTable[newIndex]), &(gr->gset[i]), sizeof(_gridset));
         }
         else
         {
            newIndex++;
            inserted = -1;
            while (inserted == -1) /** && curIndex != (newIndex-1)) (a reverifier-- Yves-20120228)**/
            {
               fprintf(stderr, "reallocate_gridset_table -- should not be here\n ");
               if (newTable[newIndex].gdin == NULL)
               {
                  memcpy(&(newTable[newIndex]), &(gr->gset[i]), sizeof(_gridset));
                  inserted = 1;
               }
               else
               {
                  newIndex++;
                  if (0 == (newIndex % newChunkSize))
                  {
                     newIndex = 0;
                  }
               }
            }
         }
      }
   }

   free(gr->gset);
   gr->gset = newTable;
   gr->log_chunk_gdin++;
}

void allocate_gridset_table(TGeoRef* gr)
{
   wordint i;
   int chunkSize;

   chunkSize = primes[cur_log_chunk];
   gr->gset = (_gridset *) calloc(sizeof(_gridset), chunkSize);

   for (i=0 ; i < chunkSize; i++)
   {
      gr->gset[i].gdin = NULL;
   }
}
