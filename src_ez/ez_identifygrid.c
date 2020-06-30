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

void c_ez_manageGrillesMemory() {
   int nchunks;
   nchunks = nGrilles / (CHUNK * CHUNK);
   if (nchunks > 0 && 0 == (nGrilles % (CHUNK * CHUNK))) {
      Grille = (TGeoRef **) realloc(Grille, CHUNK * (nchunks+1) * sizeof(TGeoRef *));
   }

   if (0 == (nGrilles % (CHUNK))) {
      Grille[(nGrilles >> LOG2_CHUNK)] = (TGeoRef *) malloc(CHUNK * sizeof(TGeoRef));
   }
}


wordint c_ezidentify_reg_grid(wordint ni, wordint nj, char* grtyp, 
   wordint ig1, wordint ig2, wordint ig3, wordint ig4, TGeoRef* GRef) 
{
   wordint i;
   wordint gdid, gdrow, gdcol, nchunks, newGrid;
   wordint res1, res2, newgrsize, grid_index;
   unsigned int grid_crc;
   char typeGrille;
   TGeoRef* gr;
   TGeoRef newgr;

   gdid = -1;
   newGrid = 0;
   typeGrille = grtyp[0];

/*    if (nGrilles == 0) {
      Grille = (TGeoRef **) calloc(chunks[cur_log_chunk], sizeof(TGeoRef *));
      Grille[0] = (TGeoRef *) calloc(chunks[cur_log_chunk], sizeof(TGeoRef));
      for (i=0; i < chunks[cur_log_chunk]; i++) {
         Grille[0][i].index = -1;
      }
   }

   memset((void *)&newgr, 0, (size_t)sizeof(TGeoRef)); */
   GRef->grtyp[0] = grtyp[0];
   GRef->grref[0] = (char) 0;
   GRef->ni = ni;
   GRef->nj = nj;
   GRef->fst.ig[IG1] = ig1;
   GRef->fst.ig[IG2] = ig2;
   GRef->fst.ig[IG3] = ig3;
   GRef->fst.ig[IG4] = ig4;
   GRef->fst.igref[IG1] = 0;
   GRef->fst.igref[IG2] = 0;
   GRef->fst.igref[IG3] = 0;
   GRef->fst.igref[IG4] = 0;


/*    newgrsize = sizeof(TGeoRef);
   grid_crc = ez_calc_crc((int *)&newgr, &newgrsize, NULL, NULL, 0, 0);
   grid_index = grid_crc % primes_sq[cur_log_chunk]; */

   gdid = c_ez_addgrid(GRef);
/*    if (gr_list[grid_index] == NULL) {
      gdid = c_ez_addgrid(grid_index, &newgr);
      return gdid;
   } else {
      gdid = c_ez_findgrid(grid_index, &newgr);
      if (gdid == -1) {
         gdid = c_ez_addgrid(grid_index, &newgr);
         return gdid;
      } else {
         return gdid;
      }
   } */
   return gdid;
}

wordint c_ezidentify_irreg_grid(
      wordint ni, wordint nj, char* grtyp, char* grref,
      wordint ig1, wordint ig2, wordint ig3, wordint ig4,
      ftnfloat* ax, ftnfloat* ay, TGeoRef* GRef) 
{
   wordint i;
   wordint gdid, gdrow, gdcol, nchunks, newGrid;
   wordint res1, res2, newgrsize, npts, grid_index;
   unsigned int grid_crc;
   char typeGrille;
   TGeoRef *gr, newgr;

   gdid = -1;
   newGrid = 0;
   typeGrille = grtyp[0];

/*    if (nGrilles == 0) {
      Grille = (TGeoRef **) calloc(chunks[cur_log_chunk],sizeof(TGeoRef *));
      Grille[0] = (TGeoRef *) calloc(chunks[cur_log_chunk], sizeof(TGeoRef));
      for (i = 0; i < chunks[cur_log_chunk]; i++) {
         Grille[0][i].index = -1;
      }
   }

   memset((void *)&newgr, (int)0, sizeof(TGeoRef)); */
   GRef->grtyp[0] = grtyp[0];
   GRef->grtyp[1] = '\0';
   GRef->grref[0] = grref[0];
   GRef->grref[1] = '\0';
   GRef->ni = ni;
   GRef->nj = nj;
   GRef->fst.ip1      = 0;
   GRef->fst.ip2      = 0;
   GRef->fst.ip3      = 0;
   GRef->fst.igref[IG1] = ig1;
   GRef->fst.igref[IG2] = ig2;
   GRef->fst.igref[IG3] = ig3;
   GRef->fst.igref[IG4] = ig4;
   GRef->fst.ig[IG1] = ig1;
   GRef->fst.ig[IG2] = ig2;
   GRef->fst.ig[IG3] = ig3;
   GRef->fst.ig[IG4] = ig4;
   GRef->fst.xg[IG1]  = 0.0;
   GRef->fst.xg[IG2]  = 0.0;
   GRef->fst.xg[IG3]  = 0.0;
   GRef->fst.xg[IG4]  = 0.0;
   GRef->nsubgrids = 0;
   strcpy(GRef->fst.nomvarx, "    ");
   strcpy(GRef->fst.nomvary, "    ");
   strcpy(GRef->fst.etiketx, "            ");
   strcpy(GRef->fst.etikety, "            ");
   strcpy(GRef->fst.typvarx, "  ");
   strcpy(GRef->fst.typvary, "  ");
   GRef->fst.deet    = 0;
   GRef->fst.npas    = 0;
   GRef->fst.nbits   = 0;
   GRef->fst.date    = 0;
   GRef->i1=1;
   GRef->i2=ni;
   GRef->j1=1;
   GRef->j2=nj;
   GRef->idx_last_gdin = -1;

/*    newgrsize = sizeof(TGeoRef); */
   switch(typeGrille) {
      case '#':
/*          grid_crc = ez_calc_crc((int *)&newgr, &newgrsize, &(ax[ig3-1]), &(ay[ig4-1]), ni, nj); */
         GRef->ax = ax;
         GRef->ay = ay; 
         break;

      case 'Y':
         npts = ni * nj;
/*          grid_crc = ez_calc_crc((int *)&newgr, &newgrsize, ax, ay, npts, npts); */
         GRef->ax = ax;
         GRef->ay = ay; 
         break;

      case 'Z':
         f77name(cigaxg)(&(GRef->grref),
            &GRef->fst.xgref[XLAT1], &GRef->fst.xgref[XLON1], &GRef->fst.xgref[XLAT2], &GRef->fst.xgref[XLON2],
            &GRef->fst.igref[IG1],   &GRef->fst.igref[IG2],   &GRef->fst.igref[IG3],   &GRef->fst.igref[IG4],1);
/*          grid_crc = ez_calc_crc((int *)&newgr, &newgrsize, ax, ay, ni, nj); */
         GRef->ax = ax;
         GRef->ay = ay; 
         break;

      case 'G':
/*          grid_crc = ez_calc_crc((int *)&newgr, &newgrsize, NULL, NULL, 0, 0); */
         break;

      default :
         fprintf(stderr, "c_ezidentify_irreg_grid : undefined grid type : %c\n", typeGrille);
         exit(13);
   }

/*    grid_index = grid_crc % primes_sq[cur_log_chunk]; */

   gdid = c_ez_addgrid(GRef);
/*    if (gr_list[grid_index] == NULL) {
      gdid = c_ez_addgrid(grid_index, &newgr);
      return gdid;
   } else {
      gdid = c_ez_findgrid(grid_index, &newgr);
      if (gdid == -1) {
         gdid = c_ez_addgrid(grid_index, &newgr);
         return gdid;
      } else {
         return gdid;
      }
   } */
   return gdid;
}
