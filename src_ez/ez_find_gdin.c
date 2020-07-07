
#include "ezscint.h"
#include "ez_funcdef.h"
#include "../src/GeoRef.h"

int c_find_gdin(TGeoRef* gdin, TGeoRef* gdout)
{
  int i; 
  TGeoRef *idx_gdin;
  _gridset *gset;

  idx_gdin = gdout->idx_last_gdin;
  if (idx_gdin == NULL)
  {
    c_ezdefset(gdout, gdin);
    // idx_gdin = gdin % primes[gdout->log_chunk_gdin];
  
/*   idx_gdin = gdout->idx_last_gdin;*/
  }
    
  gset = gdout->gset;
  // if (gset[idx_gdin].gdin == gdin) 
  // {
  //   return idx_gdin;
  // }
   
  // idx_gdin = gdin % primes[cur_log_chunk];
  i = 0;
  if (gset[i].gdin == gdin) 
  {
    return i;
  }
  
  i++; 
  while (i < gdout->n_gdin)
  {
    if (gset[i].gdin == gdin)
    {
      gdout->idx_last_gdin = gdin;
      return i;
    }
    i++;
  }
  
  return -1;
}
