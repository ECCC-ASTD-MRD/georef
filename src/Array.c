#include <malloc.h>
#include "georef/Array.h"

/*----------------------------------------------------------------------------
 * Nom      : <T3DArray_Alloc>
 * Creation : Aout 1998 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Allocation d'un tableau de coordonnee 3D pour une valeur.
 *
 * Parametres :
 *  <Value>   : Valeur de reference.
 *  <Size>    : Dimension du tableau (Nombre de coordonnee)
 *
 * Retour     :
 *  <TArray*> : Pointeru su le tableau.
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
T3DArray *T3DArray_Alloc(double Value,unsigned long Size) {

   T3DArray *array=(T3DArray*)malloc(sizeof(T3DArray));

   if (array) {
      array->Data=(Vect3d*)calloc(Size,sizeof(Vect3d));

      if (array->Data) {
         array->Size=Size;
         array->Value=Value;
      } else {
         free(array);
         array=NULL;
      }
   }
   return(array);
}

/*----------------------------------------------------------------------------
 * Nom      : <T3DArray_Free>
 * Creation : Aout 1998 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Liberation d'un tableau de coordonnee 3D.
 *
 * Parametres :
 *  <Array>   : Tableau a liberer.
 *
 * Retour     :
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
void T3DArray_Free(T3DArray *Array) {

   if (Array) {
      if (Array->Data)
         free(Array->Data);
      free(Array);
   }
}
