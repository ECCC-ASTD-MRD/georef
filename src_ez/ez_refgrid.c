#include <rpnmacros.h>
#include "ezscint.h"
#include "ez_funcdef.h"


#ifdef MUTEX
// JP
static pthread_mutex_t EZ_MTX=PTHREAD_MUTEX_INITIALIZER;
#endif

int c_ez_refgrid(TGeoRef *GRef)
{
#ifdef MUTEX
// JP
   pthread_mutex_lock(&EZ_MTX);
#endif
// JP
  GRef->access_count++;
#ifdef MUTEX
// JP
   pthread_mutex_unlock(&EZ_MTX);
#endif

   return GRef->access_count;
}