#include <string.h>

#include "GeoRef_Utils.h"

/*----------------------------------------------------------------------------
 * Nom      : <InterpCubic>
 * Creation : Octobre 1998 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Interpolation Cubique.
 *
 * Parametres :
 *  <X0>      :
 *  <X1>      :
 *  <X2>      :
 *  <X3>      :
 *  <F>       :
 *
 * Retour:
 *  <R>  : Valeur interpolee
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
double InterpCubic(double X0,double X1,double X2, double X3,double F) {

   double a0,a1,a2,f2;

   f2=F*F;
   a0=X3-X2-X0+X1;
   a1=X0-X1-a0;
   a2=X2-X0;

   return(a0*F*f2+a1*f2+a2*F+X1);
}

/*----------------------------------------------------------------------------
 * Nom      : <InterpHermite>
 * Creation : Octobre 1998 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Interpolation Hermite.
 *
 * Parametres :
 *  <X0>      :
 *  <X1>      :
 *  <X2>      :
 *  <X3>      :
 *  <F>       :
 *
 * Retour:
 *  <R>  : Valeur interpolee
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
double InterpHermite(double X0,double X1,double X2, double X3,double F,double T,double B) {

  double a0,a1,a2,a3,f0,f1,f2,f3;

   f2=F*F;
   f3=f2*F;

   f0=(X1-X0)*(1+B)*(1-T)/2 + (X2-X1)*(1-B)*(1-T)/2;
   f1=(X2-X1)*(1+B)*(1-T)/2 + (X3-X2)*(1-B)*(1-T)/2;

   a0=2*f3-3*f2+1;
   a1=f3-2*f2+F;
   a2=f3-f2;
   a3=-2*f3+3*f2;

   return(a0*X1+a1*f0+a2*f1+a3*X2);
}

/*----------------------------------------------------------------------------
 * Nom      : <QSort_Double>
 * Creation : Octobre 2005 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Fonction de comparaison pour le tri.
 *
 * Parametres :
 *  <V0>      : Valeur 0
 *  <V1>      : Valeur 1
 *
 * Retour:
 *  <<>=>:      -1< 0= 1>
 *

 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
int QSort_Double(const void *A,const void *B){

   if (*(const double*)A<*(const double*)B) {
      return(-1);
   } else if (*(const double*)A>*(const double*)B) {
      return(1);
   } else {
      return(0);
   }
}

int QSort_Float(const void *A,const void *B){

   if (*(const float*)A<*(const float*)B) {
      return(-1);
   } else if (*(const float*)A>*(const float*)B) {
      return(1);
   } else {
      return(0);
   }
}

int QSort_Int(const void *A, const void *B) {
   return(*(const int*)A)-(*(const int*)B);
}

/*----------------------------------------------------------------------------
 * Nom      : <QSort_Dec*>
 * Creation : Mars 2015 - E. Legault-Ouellet - CMC/CMOE
 *
 * But      : Fonction de comparaison pour le tri decroissant
 *
 * Parametres :
 *  <V0>      : Valeur 0
 *  <V1>      : Valeur 1
 *
 * Retour:
 *      - A negative value if B < A
 *      - A positive value if B > A
 *      - Zero if B == A
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
int QSort_DecDouble(const void *A,const void *B){
   if (*(const double*)B<*(const double*)A) {
      return(-1);
   } else if (*(const double*)B>*(const double*)A) {
      return(1);
   } else {
      return(0);
   }
}

int QSort_DecFloat(const void *A,const void *B){
   if (*(const float*)B<*(const float*)A) {
      return(-1);
   } else if (*(const float*)B>*(const float*)A) {
      return(1);
   } else {
      return(0);
   }
}

int QSort_DecInt(const void *A, const void *B) {
   return(*(const int*)B)-(*(const int*)A);
}

int QSort_StrPtr(const void *A, const void *B) {
   return strcmp(*(const char * const *)A,*(const char * const *)B);
}

/*----------------------------------------------------------------------------
 * Nom      : <Unique>
 * Creation : Mars 2015 - E. Legault-Ouellet - CMC/CMOE
 *
 * But      : Réduit une liste ordonnée de valeurs en une liste ordonnée de
 *            valeurs uniques.
 *
 * Parametres :
 *  <Arr>   : Liste de valeurs (sera modifiée)
 *  <Size>  : Nombre de valeurs (sera ajusté)
 *  <NBytes>: Nombre de bytes par valeur
 *
 * Retour   :
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
void Unique(void *Arr,int* restrict Size,size_t NBytes) {
   int i;
   char *a,*b;

   if( *Size > 0 ) {
      for(i=*Size-1,a=Arr,b=a+NBytes; i; --i,b+=NBytes) {
         if( memcmp(a,b,NBytes) ) {
            a += NBytes;
            if( a != b ) {
               memcpy(a,b,NBytes);
            }
         }
      }

      *Size = (int)((a-(char*)Arr)/NBytes) + 1;
   }
}


double HCentile(double *M,int N,int K) {

   register int    i,j,l,m;
   register double x,d;

   l=0;m=N-1;

   while(l<m) {
      x=M[K];
      i=l;
      j=m;
      do {
         while (M[i]<x) i++;
         while (x<M[j]) j--;
         if (i<=j) {
            d=M[i];M[i]=M[j];M[j]=d;
            i++;
            j--;
         }
      } while (i<=j);
      if (j<K) l=i;
      if (K<i) m=j;
   }
   return(M[K]);
}
