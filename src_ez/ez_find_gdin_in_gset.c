 #include "ez_funcdef.h"
 #include "../src/GeoRef.h"
 
 wordint ez_find_gdin_in_gset(TGeoRef* gdin, TGeoRef* gdout)
  {
  int idx_gdin;
  
  idx_gdin = 0;
  while (idx_gdin < gdout->n_gdin)
  {
    if (gdin == gdout->gset[idx_gdin].gdin) {
      return idx_gdin;
    }

    idx_gdin++;
  }
  
  return -1;
  }  
