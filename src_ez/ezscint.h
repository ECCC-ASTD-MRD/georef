#ifndef _EZSCINT

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <malloc.h>
#include "rpnmacros.h"

#define NMAXGRIDS 32
#define NMAXSETS NMAXGRIDS *(NMAXGRIDS - 1)
#define NMAXSUBGRIDS 20

#define LAT 1
#define EZ_AX 2
#define XXX 4
#define NEWTON 256
#define ZONES 2048

#define SCALAIRE 0
#define VECTEUR 1

#define GLOBALE 0
#define LOCALE 1

#define NON 0
#define OUI 1

#define VOISIN 0
#define NEAREST 0
#define LINEAIRE 1
#define LINEAR 1
#define CUBIQUE 3
#define DISTANCE 4
#define TRIANGLE 5
#define LINEAR_AND_NEAREST 6

#define EZ_EXTRAP 1
#define EZ_NO_EXTRAP 0
#define RIEN -1

#define MAXIMUM 4
#define MINIMUM 5
#define VALEUR 6
#define ABORT 13

#define DEHORS 0
#define AU_NORD 1
#define AU_SUD 2
#define POLE_NORD 3
#define POLE_SUD 4

#define GLOBAL 0
#define NORD 1
#define SUD 2

#define SYM 1
#define ANTISYM 0

#define C_TO_FTN(i, j, ni) (wordint)((ni) * (j) + i)

#define ABSOLU 0
#define RELATIF 1

#define SWLAT 0
#define SWLON 1
#define DLAT 2
#define DLON 3

#define TD60 0
#define TDGRW 1
#define CLAT 2
#define CLON 3

#define PI 0
#define PJ 1
#define D60 2
#define DGRW 3

#define IG1 0
#define IG2 1
#define IG3 2
#define IG4 3

#define XLAT1 0
#define XLON1 1
#define XLAT2 2
#define XLON2 3

// TODO: Moved to Hash.h
#define CHUNK 128
#define LOG2_CHUNK 7
#define MAX_LOG_CHUNK 12

#define YES 1
#define NO 0

typedef struct
{
  wordint degre_interp;
  wordint degre_extrap;
  wordint use_1subgrid;
  wordint valeur_1subgrid;
  wordint symmetrie;
  wordint vecteur;
  wordint verbose;
  wordint polar_correction;
  wordint wgt_num;
  wordint msg_pt_tol;
  wordint cld_interp_alg;
  ftnfloat msg_dist_thresh;
  ftnfloat valeur_extrap;
} _groptions;

extern wordint nGrilles;
extern wordint cur_log_chunk;

// These declarations used to have the __thread storage class, but threads
// aren't actually used.  Furthemore, the PGI compiler does not support that
// storage class
extern _groptions groptions;

extern wordint log_chunks[];
extern wordint primes[];
extern wordint chunks[];
extern wordint primes_sq[];
extern wordint chunks_sq[];

#endif

#define _EZSCINT
