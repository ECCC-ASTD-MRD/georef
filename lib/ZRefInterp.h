#ifndef _ZRefInterp_h
#define _ZRefInterp_h

#include "ZRef.h"

//  Interpolation 
#define ZRNEAREST_NEIGHBOUR 0x001
#define ZRLINEAR            0x002
#define ZRCUBIC_WITH_DERIV  0x004
#define ZRCUBIC_LAGRANGE    0x008

//  Extrapolation 
#define ZRCLAMPED    0x010
#define ZRLAPSERATE  0x020

//  Other options 
#define ZRVERBOSE    0x040

//  check for float exception (will do nothing on SX6) 
#define ZRCHECKFLOAT 0x080

typedef struct TZRefInterp {
   TZRef         *ZRefSrc,*ZRefDest;       // Source and destination vertical references
   int           *Indexes;                 // Interpolation array indices
   int           NIJ;                      // 2D dimensions
   int           Same;                     // Flag indicating source and destination are the same
} TZRefInterp;

int          ZRefInterp_Free(TZRefInterp *Interp);
void         ZRefInterp_Clear(TZRefInterp *Interp);
TZRefInterp *ZRefInterp_Define(TZRef *ZRefDest,TZRef *ZRefSrc,const int NI,const int NJ);
int          ZRefInterp(TZRefInterp *Interp,float *stateOut,float *stateIn,float *derivOut,float *derivIn,float extrapGuideDown,float extrapGuideUp);
int          ZRefInterp_SetOption(const char *Option,const char *Value);
int          ZRefInterp_SetOptioni(const unsigned char option);

#endif
