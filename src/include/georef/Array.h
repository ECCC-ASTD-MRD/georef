#ifndef _Array_h
#define _Array_h

#include "rmn/Vector.h"

typedef struct T3DArray {
    Vect3d * Data;
    double Value;
    unsigned long Size;
} T3DArray;

T3DArray *T3DArray_Alloc(
    double value,
    unsigned long size
);
void T3DArray_Free(
    T3DArray * const array
);

#endif
