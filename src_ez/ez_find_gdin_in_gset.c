 #include "ezscint.h"
 #include "ez_funcdef.h"
 #include "../src/GeoRef.h"
 
 wordint ez_find_gdin_in_gset(TGeoRef* GRef, wordint gdout)
  {
  int i, found, idx_gdin, gdin; 
  int gdrow_out, gdcol_out;
  TGeoRef *gr; 

  c_gdkey2rowcol(gdout, &gdrow_out, &gdcol_out);
  gr = &(Grille[gdrow_out][gdcol_out]);

  /* TODO: remove gdin */
  gdin = GRef->index;
  idx_gdin = gdin % primes[gr->log_chunk_gdin];
  if (gr->gset[idx_gdin].gdin == gdin)
  {
  return idx_gdin;
  }
  
  i = idx_gdin;
  found = -1;
  while ((found == -1) && (i != idx_gdin-1) && (gr->gset[i].gdin != -1))
  {
  if (gdin == gr->gset[i].gdin)
    {
    found = 1;
    idx_gdin = i;
    return idx_gdin;
    }
  else
    {
    i++;
    if (0 == (i % primes[cur_log_chunk]))
      {
      i = 0;
      }
    }
  }
  return -1;
  }  
