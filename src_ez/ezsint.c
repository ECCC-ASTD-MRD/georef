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
// Need #include "../src/GeoRef.h" ?

TGeoRef** Grille  = NULL;
wordint nGrilles = 0;
wordint nGrillesMax = CHUNK*CHUNK;
wordint cur_log_chunk = 7;

// These declarations used to have the __thread storage class, but threads
// aren't actually used.  Furthemore, the PGI compiler does not support that
// storage class
wordint iset     = -1;
wordint iset_gdin = -1;
wordint iset_gdout = -1;
_gridset* gridset = NULL;
_groptions groptions = { OUI, CUBIQUE,  MAXIMUM, NON, -1, SYM, SCALAIRE, NON, NON, OUI, 16, 0, DISTANCE, NEAREST, 0.5, 3.0, 0.0  };

wordint log_chunks[]= {0, 1, 2, 3,   4,    5,   6,      7,     8,      9,      10,     11,        12};
wordint primes[]    = {0, 0, 3, 7,  13,   31,   61,   127,   251,    509,    1021,   2039,      4093};
wordint chunks[]    = {0, 0, 4, 8,  16,   32,   64,   128,   256,    512,    1024,   2048,      4096};
wordint primes_sq[] = {0, 0, 3, 61, 251, 1021, 4093, 16381, 65521, 262139, 1048573, 4194301, 16777213};
wordint chunks_sq[] = {0, 0, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216};


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezsint)(ftnfloat *zout, ftnfloat *zin)
{
   wordint icode;
   
   icode = c_ezsint(zout, zin);
   return icode;
}
wordint c_ezsint(ftnfloat *zout, ftnfloat *zin)
{
  wordint icode,gdin,gdout;
  wordint gdrow_in,gdcol_in, gdrow_out,gdcol_out;
   
  if (iset_gdin == UNDEFINED || iset_gdout == UNDEFINED)
    {
    fprintf(stderr,"<c_ezsint> Source or target grid undefined! Aborting...\n");
    return -1;
    }
  
  
  gdin = iset_gdin;
  gdout= iset_gdout;
  
  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
   
  if (iset_gdin == iset_gdout)
    {
    memcpy(zout, zin, Grille[gdrow_in][gdcol_in].ni*Grille[gdrow_in][gdcol_in].nj*sizeof(ftnfloat));
    return 1;
    }


  if (Grille[gdrow_in][gdcol_in].nsubgrids > 0 || Grille[gdrow_out][gdcol_out].nsubgrids > 0)
      {
/* get the subgrids and interpolate accordingly */
      icode = c_ezyysint(zout,zin,gdout,gdin);
      iset_gdin=gdin;
      iset_gdout=gdout;
      return icode;
      }
  icode = c_ezsint_orig(zout, zin);
  return icode;
}

wordint c_ezsint_orig(ftnfloat *zout, ftnfloat *zin)
{
  wordint gdin, gdout;
  wordint ier,ierc;
  wordint npts;
  wordint gdrow_in, gdrow_out, gdcol_in, gdcol_out;
  ftnfloat *lzin, *lxzin;
  
  lzin  = NULL;
  lxzin = NULL;
  ierc  = 0;
  
  if (iset_gdin == UNDEFINED || iset_gdout == UNDEFINED)
    {
    fprintf(stderr,"<c_ezsint_orig> Source or target grid undefined! Aborting...\n");
    return -1;
    }
  
  
  gdin = iset_gdin;
  gdout= iset_gdout;
  
  c_gdkey2rowcol(gdin,  &gdrow_in,  &gdcol_in);
  c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
   
  if (iset_gdin == iset_gdout)
    {
    memcpy(zout, zin, Grille[gdrow_in][gdcol_in].ni*Grille[gdrow_in][gdcol_in].nj*sizeof(ftnfloat));
    return 1;
    }
  
  if (Grille[gdrow_in][gdcol_in].fst.axe_y_inverse == 1)
    {
    lzin = (ftnfloat *) malloc(Grille[gdrow_in][gdcol_in].ni*Grille[gdrow_in][gdcol_in].nj*sizeof(ftnfloat));
    memcpy(lzin, zin, Grille[gdrow_in][gdcol_in].ni*Grille[gdrow_in][gdcol_in].nj*sizeof(ftnfloat));
    f77name(permut)(lzin, &Grille[gdrow_in][gdcol_in].ni, &Grille[gdrow_in][gdcol_in].nj);
    }
  else
    {
    lzin = zin;
    }
  
  if (Grille[gdrow_in][gdcol_in].needs_expansion == OUI)
    {
    lxzin = (ftnfloat *) malloc(2*Grille[gdrow_in][gdcol_in].ni*Grille[gdrow_in][gdcol_in].nj*sizeof(ftnfloat));
    ez_xpnsrcgd(gdin, lxzin, lzin);
    }
  else
    {
    lxzin = lzin;
    }
  
  ier = ez_calclatlon(gdout);
  ier = ez_calcxy(gdin, gdout);
  npts = Grille[gdrow_out][gdcol_out].ni*Grille[gdrow_out][gdcol_out].nj;
  
  ier = ez_interp(zout, lxzin, gdin, gdout);
  
  if (groptions.polar_correction == OUI)
    {
    ier = ez_defzones(gdin, gdout);
    ierc= ez_corrval(zout, lxzin, gdin, gdout);
    }
  
  if (lzin != zin && lzin != NULL)
    {
    free(lzin);
    }
  
  if (lxzin != lzin && lxzin != zin && lxzin != NULL)
    {
    free(lxzin);
    }
  
  return ierc;
}
