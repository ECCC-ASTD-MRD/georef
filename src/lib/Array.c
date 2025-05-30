//! \file

#include <malloc.h>
#include "georef/Array.h"


//! Allocate 3D coordinate array for a value
T3DArray *T3DArray_Alloc(
    //! [in] 
    double value,
    //! [in] Number of elements in array
    unsigned long size
) {
    //! \return T3DArray pointer
    T3DArray * const array = (T3DArray*)malloc(sizeof(T3DArray));

    if (array) {
        array->Data = (Vect3d*)calloc(size, sizeof(Vect3d));

        if (array->Data) {
            array->Size = size;
            array->Value = value;
        } else {
            free(array);
        }
    }
    return array;
}


//! Free 3D coordinate array
void T3DArray_Free(
    //! [in] Array to free
    T3DArray * const array
) {
    if (array) {
        if (array->Data) free(array->Data);
        free(array);
    }
}
