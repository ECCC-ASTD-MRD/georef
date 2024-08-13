#include <malloc.h>
#include <math.h>
#include <strings.h>
#include <stdlib.h>

#include <App.h>
#include "Def.h"
#include "GeoRef.h"
#include "GeoRef_Utils.h"

// Sizes in bytes of the different data types
// TODO: revisit for architecture dependencies
int         TDef_Size[]          = { 0,1,1,1,2,2,4,4,8,8,4,8 };
int         TDef_Mult[]          = { 0,0,0,0,1,1,2,2,3,3,2,3 };
int         TDef_DTYP[]          = { 0,2,4,2,4,2,4,2,4,5,5,-1 };

/*----------------------------------------------------------------------------
 * Nom      : <Def_Clear>
 * Creation : Fevrier 2003- J.P. Gauthier - CMC/CMOE
 *
 * But      : Reinitialiser les valeurs a nodata
 *
 * Parametres :
 *  <Def>     : Structure a reinitialiser
 *
 * Retour:
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
void Def_Clear(TDef *Def){

   int n,i;

   for(n=0;n<Def->NC;n++) {
      for(i=0;i<FSIZE3D(Def);i++) {
         Def_Set(Def,n,i,Def->NoData);
      }
   }
   
   if (Def->Buffer) {
      free(Def->Buffer);
      Def->Buffer=NULL;
   }
   if (Def->Aux) {
      free(Def->Aux);
      Def->Aux=NULL;
   }
   if (Def->Accum) {
      free(Def->Accum);
      Def->Accum=NULL;
   }
   if (Def->Mask) {
      free(Def->Mask);
      Def->Mask=NULL;
   }
   if (Def->Sub) {
      free(Def->Sub);
      Def->Sub=NULL;
   }

   if (Def->Segments) {
      TList_Clear(Def->Segments,(int(*)(void*))T3DArray_Free);
      Def->Segments=NULL;
   }
}

/*----------------------------------------------------------------------------
 * Nom      : <Def_Compat>
 * Creation : Mars 2009- J.P. Gauthier - CMC/CMOE
 *
 * But      : Verifier les dimensiont entre 2 Def et les ajuster en consequence.
 *
 * Parametres :
 *  <ToDef>   : Definition a redimensionner
 *  <DefTFrom>: Definition de laquelle redimensionner
 *
 * Retour:
 *  <Compat>  : Compatibles?
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
int Def_Compat(TDef *ToDef,TDef *FromDef) {

   int ch=1,n,nijk;

   if (FromDef->Idx) {
      return(0);
   }

   if (ToDef->Mode && ToDef->Mode!=ToDef->Data[0]) {
      free(ToDef->Mode);
   }
   ToDef->Mode=NULL;
   ToDef->Dir=NULL;
   
   // Verifier la dimension verticale
   if (ToDef->NK!=FromDef->NK || ToDef->NC!=FromDef->NC) {
      if (ToDef->Idx) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Cannot change the data size for a sub 'U' grid\n",__func__);
         return(0);
      }
      if (ToDef->Data[0]) free(ToDef->Data[0]);
      ToDef->Data[0]=ToDef->Data[1]=ToDef->Data[2]=ToDef->Data[3]=NULL;
      
      ToDef->NK=FromDef->NK;
      ToDef->NC=FromDef->NC;
      nijk=FSIZE3D(ToDef);
      if (!(ToDef->Data[0]=(char*)calloc(nijk*ToDef->NC,TDef_Size[ToDef->Type]))) {
         return(0);
      }
      for(n=1;n<ToDef->NC;n++) {
         ToDef->Data[n]=&ToDef->Data[0][nijk*n*TDef_Size[ToDef->Type]];
      }
      ch=0;      
   }

   // Verifier le masque
   if (FromDef->Mask) {
      if (!ToDef->Mask || ch==0) {
         if (!(ToDef->Mask=(char*)realloc(ToDef->Mask,FSIZE3D(ToDef)))) {
            return(0);
         }
      }
   } else if (ToDef->Mask) {
      free(ToDef->Mask); 
      ToDef->Mask=NULL;
   }

   return(ch);
}

/*----------------------------------------------------------------------------
 * Nom      : <Def_Copy>
 * Creation : Fevrier 2003- J.P. Gauthier - CMC/CMOE
 *
 * But      : Copier une structure TDef.
 *
 * Parametres :
 *  <Def>     : Structure a copier
 *
 * Retour:
 *  <def>     : Pointeur sur le copie de la structure
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
TDef *Def_Copy(TDef *Def){

   int   n,nijk;
   TDef *def=NULL;

   if (!Def->Idx) {
      if (Def && (def=(TDef*)malloc(sizeof(TDef)))) {
         def->Alias=Def->Alias;
         def->CellDim=Def->CellDim;
         def->NI=Def->NI;
         def->NJ=Def->NJ;
         def->NK=Def->NK;
         def->NC=Def->NC;
         def->NIJ=Def->NI*Def->NJ;
         def->NoData=Def->NoData;
         def->Type=Def->Type;
         def->Level=Def->Level;
         def->Idx=Def->Idx;
         def->Buffer=NULL;
         def->Aux=NULL;
         def->Segments=NULL;
         def->Accum=NULL;
         def->Mask=NULL;
         def->Sub=NULL;
         def->Pres=NULL;
         def->PresLS=NULL;
         def->Height=NULL;
         def->Pick=def->Poly=NULL;
         def->Sample=Def->Sample;
         def->SubSample=Def->SubSample;

         memcpy(def->Limits,Def->Limits,6*sizeof(int));
         def->CoordLimits[0][0]=Def->CoordLimits[0][0];
         def->CoordLimits[0][1]=Def->CoordLimits[0][1];
         def->CoordLimits[1][0]=Def->CoordLimits[1][0];
         def->CoordLimits[1][1]=Def->CoordLimits[1][1];

         nijk=FSIZE3D(Def);
         def->Data[0]=def->Data[1]=def->Data[2]=def->Data[3]=NULL;
         if (def->Alias) {
            for(n=0;n<def->NC;n++) def->Data[n]=Def->Data[n];
         } else {
            if (!(def->Data[0]=(char*)malloc(nijk*def->NC*TDef_Size[def->Type]))) {
               Def_Free(def);
               return(NULL);
            }
            memcpy(def->Data[0],Def->Data[0],nijk*def->NC*TDef_Size[Def->Type]);
            for(n=1;n<def->NC;n++) {
               def->Data[n]=&def->Data[0][nijk*n*TDef_Size[Def->Type]];
            }
         }
         def->Mode=def->Data[0];
         def->Dir=NULL;
         
         if (Def->Mask) { 
            if (!(def->Mask=(char*)malloc(nijk))) {
               return(NULL);
            }
            memcpy(def->Mask,Def->Mask,nijk);
         }
      }
   }
   return(def);
}

TDef *Def_CopyPromote(TDef *Def,TDef_Type Type){

   int   n,nijk;
   TDef *def=NULL;

      if (Def && !Def->Idx && (def=(TDef*)malloc(sizeof(TDef)))) {
         def->Alias=0;
         def->CellDim=Def->CellDim;
         def->NI=Def->NI;
         def->NJ=Def->NJ;
         def->NK=Def->NK;
         def->NC=Def->NC;
         def->NIJ=Def->NI*Def->NJ;
         def->NoData=Def->NoData;
         def->Type=Type;
         def->Level=Def->Level;
         def->Idx=Def->Idx;
         def->Buffer=NULL;
         def->Aux=NULL;
         def->Segments=NULL;
         def->Accum=NULL;
         def->Mask=NULL;
         def->Sub=NULL;
         def->Pres=NULL;
         def->PresLS=NULL;
         def->Height=NULL;
         def->Pick=def->Poly=NULL;
         def->Sample=Def->Sample;
         def->SubSample=Def->SubSample;

         memcpy(def->Limits,Def->Limits,6*sizeof(int));
         def->CoordLimits[0][0]=Def->CoordLimits[0][0];
         def->CoordLimits[0][1]=Def->CoordLimits[0][1];
         def->CoordLimits[1][0]=Def->CoordLimits[1][0];
         def->CoordLimits[1][1]=Def->CoordLimits[1][1];

         nijk=FSIZE3D(Def);
         def->Data[0]=def->Data[1]=def->Data[2]=def->Data[3]=NULL;
         if (!(def->Data[0]=(char*)calloc(nijk*def->NC,TDef_Size[def->Type]))) {
            Def_Free(def);
            return(NULL);
         }
         for(n=1;n<def->NC;n++) {
            def->Data[n]=&def->Data[0][nijk*n*TDef_Size[def->Type]];
         }

         def->Mode=def->Data[0];
         def->Dir=NULL;
         
         if (Def->Mask) { 
            if (!(def->Mask=(char*)malloc(nijk))) {
               return(NULL);
            }
            memcpy(def->Mask,Def->Mask,nijk);
         }         
      }

   return(def);
}

TDef *Def_CopyData(TDef *ToDef,TDef *FromDef){

   TDef *def=ToDef;
   int   n,alias=FALSE;
//TODO: define Alias
   if (!def) {
      def=Def_Copy(FromDef);
   }

   if (!FromDef) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid data source\n",__func__);
      return(NULL);
   }
   //TODO:Check Def COMPAT

   for(n=0;n<FromDef->NC;n++) {
      if (FromDef->Data[n]) {
         if (alias) {
            def->Data[n]=FromDef->Data[n];
         } else {
            memcpy(def->Data[n],FromDef->Data[n],FSIZE3D(FromDef)*TDef_Size[FromDef->Type]);
         }
      }
   }
   if (alias)
      def->Mode=FromDef->Mode;

   if (FromDef->Mask) {
      def->Mask=(char*)malloc(FromDef->NIJ);
      memcpy(def->Mask,FromDef->Mask,FromDef->NIJ);
   }
   return (def);
}

/*----------------------------------------------------------------------------
 * Nom      : <Def_Free>
 * Creation : Fevrier 2003- J.P. Gauthier - CMC/CMOE
 *
 * But      : Liberer la structure TDef.
 *
 * Parametres :
 *  <Def>     : Structure a liberer
 *
 * Retour:
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
void Def_Free(TDef *Def){

   if (Def) {
      if (!Def->Alias && !Def->Idx) {
         if (Def->Mode && Def->Mode!=Def->Data[0]) free(Def->Mode);
         if (Def->Data[0])                         free(Def->Data[0]);
         if (Def->Mask)                            free(Def->Mask);
      }

      if (Def->Buffer)             free(Def->Buffer);
      if (Def->Aux)                free(Def->Aux);
      if (Def->Accum)              free(Def->Accum);
      if (Def->Sub)                free(Def->Sub);
      if (Def->Pres>(float*)0x1)   free(Def->Pres);
      if (Def->PresLS>(float*)0x1) free(Def->PresLS);
      if (Def->Height>(float*)0x1) free(Def->Height);
#ifdef HAVE_GDAL
      if (Def->Poly)               OGR_G_DestroyGeometry(Def->Poly);
//Freed by Def->Poly      if (Def->Pick)       OGR_G_DestroyGeometry(Def->Pick);
#endif
      if (Def->Segments)           TList_Clear(Def->Segments,(int(*)(void*))T3DArray_Free);

      free(Def);
   }
}

/*----------------------------------------------------------------------------
 * Nom      : <Def_New>
 * Creation : Fevrier 2003- J.P. Gauthier - CMC/CMOE
 *
 * But      : Initialiser la structure TDef.
 *
 * Parametres :
 *  <NI>      : Dimension du champs a allouer
 *  <NJ>      : Dimension du champs a allouer
 *  <NK>      : Dimension du champs a allouer
 *  <Dim>     : Nombre de composantes
 *  <Type>    : Type de donnnes
 *  <Alias>   : Alias data arrays
 *
 * Retour:
 *  <Def>:      Nouvelle structure
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
TDef *Def_New(int NI,int NJ,int NK,int Dim,TDef_Type Type,int Alias) {

   int   n,nijk;
   TDef *def;

   if (!(def=(TDef*)malloc(sizeof(TDef))))
     return(NULL);

   def->NI=NI;
   def->NJ=NJ;
   def->NK=NK;
   def->NIJ=NI*NJ;
   def->NC=abs(Dim);
   def->Alias=Alias;
   def->CellDim=2;
   def->NoData=nan("NaN");
   def->Level=0;
   def->Idx=0;

   def->Limits[0][0]=0;
   def->Limits[1][0]=0;
   def->Limits[2][0]=0;
   def->Limits[0][1]=NI-1;
   def->Limits[1][1]=NJ-1;
   def->Limits[2][1]=NK-1;

   def->CoordLimits[0][0]=-180;
   def->CoordLimits[0][1]=180;
   def->CoordLimits[1][0]=-90;
   def->CoordLimits[1][1]=90;

   def->SubSample=def->Sample=1;

   def->Type=Type;
   def->Buffer=NULL;
   def->Aux=NULL;
   def->Segments=NULL;
   def->Accum=NULL;
   def->Mask=NULL;
   def->Sub=NULL;
   def->Pres=NULL;
   def->PresLS=NULL;
   def->Height=NULL;
   def->Pick=def->Poly=NULL;

   // Allocate a single buffer and set component pointer to right address
   nijk=FSIZE3D(def);
   def->Data[0]=def->Data[1]=def->Data[2]=def->Data[3]=def->Dir=NULL;
   
   if (!def->Alias) {
      if (!(def->Data[0]=(char*)calloc(nijk*def->NC,TDef_Size[Type]))) {
         Def_Free(def);
         return(NULL);
      }
      
      for(n=1;n<def->NC;n++) {
         def->Data[n]=&def->Data[0][nijk*n*TDef_Size[def->Type]];
      }
   }
   def->Mode=def->Data[0];

   return(def);
}

TDef *Def_Create(int NI,int NJ,int NK,TDef_Type Type,char* Comp0,char* Comp1,char* Mask) {
   
   TDef *def=NULL;
   int dim;

   if (def=(TDef*)malloc(sizeof(TDef))) {
      dim=(Comp0!=NULL)+(Comp1!=NULL);
      def=Def_New(NI,NJ,NK,dim,Type,TRUE);

      def->Data[0]=Comp0;
      def->Data[1]=Comp1;
      def->Mask=Mask;
   }
   return (def);
}
/*
static inline float* Def_ToFloat(TDef *Def,char* Buffer) {

   int n,sz;
   float *buf=NULL;

   sz=FSIZE3D(Def);
   *buf=(float*)malloc(Size*sizeof(float));

   switch(Def->Type) {
      case TD_Unknown:break;
      case TD_Binary: break;
      case TD_UByte:  for(n=0;n<sz;n++) buf[n]=((unsigned char*)Buffer)[n];break;
      case TD_Byte:   for(n=0;n<sz;n++) buf[n]=((char*)Buffer)[n];break;
      case TD_UInt16: for(n=0;n<sz;n++) buf[n]=((unsigned short*)Buffer)[n];break;
      case TD_Int16:  for(n=0;n<sz;n++) buf[n]=((short*)Buffer)[n];break;
      case TD_UInt32: for(n=0;n<sz;n++) buf[n]=((unsigned int*)Buffer)[n];break;
      case TD_Int32:  for(n=0;n<sz;n++) buf[n]=((int*)Buffer)[n];break;
      case TD_UInt64: for(n=0;n<sz;n++) buf[n]=((unsigned long long*)Buffer)[n];break;
      case TD_Int64:  for(n=0;n<sz;n++) buf[n]=((long long*)Buffer)[n];break;
      case TD_Float32: buf=(float*)Buffer;break; 
      case TD_Float64: for(n=0;n<sz;n++) buf[n]=((double*)Buffer)[n];break;
   }

   return(buf);
}
static inline float* Def_FromFloat(TDef *Def) {

   int n;

   sz=FSIZE3D(Def);
   switch(Def->Type) {
      case TD_Unknown:break;
      case TD_Binary: break;
      case TD_UByte:  for(n=0;n<sz;n++) ((unsigned char*)Buffer)[n]=buf[n];break;
      case TD_Byte:   for(n=0;n<sz;n++) ((char*)Buffer)[n]=buf[n];break;
      case TD_UInt16: for(n=0;n<sz;n++) ((unsigned short*)Buffer)[n]=buf[n];break;
      case TD_Int16:  for(n=0;n<sz;n++) ((short*)Buffer)[n]=buf[n];break;
      case TD_UInt32: for(n=0;n<sz;n++) ((unsigned int*)Buffer)[n]=buf[n];break;
      case TD_Int32:  for(n=0;n<sz;n++) ((int*)Buffer)[n]=buf[n];break;
      case TD_UInt64: for(n=0;n<sz;n++) ((unsigned long long*)Buffer)[n]=buf[n];break;
      case TD_Int64:  for(n=0;n<sz;n++) ((long long*)Buffer)[n]=buf[n];break;
      case TD_Float32: (float*)Buffer=buf;break; 
      case TD_Float64: for(n=0;n<sz;n++) ((double*)Buffer)[n]=buf[n];break;
   }

   return(buf);
}
*/

/*----------------------------------------------------------------------------
 * Nom      : <Def_Resize>
 * Creation : Avril 2004- J.P. Gauthier - CMC/CMOE
 *
 * But      : Redimensionner le champs de donnees.
 *
 * Parametres :
 *  <Def>     : Definition a redimensionner
 *  <NI>      : Dimension du champs a allouer
 *  <NJ>      : Dimension du champs a allouer
 *  <NK>      : Dimension du champs a allouer
 *
 * Retour:
 *  <Def>:      Nouvelle structure
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
TDef *Def_Resize(TDef *Def,int NI,int NJ,int NK){

   int n,nijk;

   if (!Def)
     return(NULL);

   if (!Def->Alias && !Def->Idx && (Def->NI!=NI || Def->NJ!=NJ || Def->NK!=NK)) {
      Def->NI=NI;
      Def->NJ=NJ;
      Def->NK=NK;
      Def->NIJ=NI*NJ;

      Def->Limits[0][0]=0;
      Def->Limits[1][0]=0;
      Def->Limits[2][0]=0;
      Def->Limits[0][1]=NI-1;
      Def->Limits[1][1]=NJ-1;
      Def->Limits[2][1]=NK-1;

      Def->CoordLimits[0][0]=-180;
      Def->CoordLimits[0][1]=180;
      Def->CoordLimits[1][0]=-90;
      Def->CoordLimits[1][1]=90;

      Def->SubSample=Def->Sample=1;

      if (Def->Mode && Def->Mode!=Def->Data[0]) {
         free(Def->Mode);
         Def->Mode=NULL;
      }

      // Allocate a single buffer and set component pointer to right address
      Def->Data[1]=Def->Data[2]=Def->Data[3]=NULL;
      nijk=FSIZE3D(Def);
      if (!(Def->Data[0]=(char*)realloc(Def->Data[0],nijk*Def->NC*TDef_Size[Def->Type]))) {
         Def_Free(Def);
         return(NULL);
      }
      for(n=1;n<Def->NC;n++) {
         Def->Data[n]=&Def->Data[0][nijk*n*TDef_Size[Def->Type]];
      }
      Def->Mode=Def->Data[0];
      Def->Dir=NULL;
      
      if (Def->Buffer)             free(Def->Buffer); Def->Buffer=NULL;
      if (Def->Aux)                free(Def->Aux);    Def->Aux=NULL;
      if (Def->Accum)              free(Def->Accum);  Def->Accum=NULL;
      if (Def->Mask)               free(Def->Mask);   Def->Mask=NULL;
      if (Def->Sub)                free(Def->Sub);    Def->Sub=NULL;
      if (Def->Pres>(float*)0x1)   free(Def->Pres);   Def->Pres=NULL;
      if (Def->PresLS>(float*)0x1) free(Def->PresLS); Def->PresLS=NULL;
      if (Def->Height>(float*)0x1) free(Def->Height); Def->Height=NULL;
   }
   return(Def);
}

/*----------------------------------------------------------------------------
 * Nom      : <Def_Paste>
 * Creation : Novembre 2007- J.P. Gauthier - CMC/CMOE
 *
 * But      : Copier les donnees d'un TDef dans un autres (Paste).
 *
 * Parametres :
 *  <ToDef>   : Destination
 *  <DefPaste>: Paste
 *  <X0>      : Point de depart de la tuile en X
 *  <Y0>      : Point de depart de la tuile en Y
 *
 * Retour:
 *  <Code>    : In(1) ou Out(0)
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
int Def_Paste(TDef *ToDef,TDef *DefPaste,int X0, int Y0) {

   int           x,y,dx,dy,x0,y0,x1,y1,c,nc;
   unsigned long idxs,idxd;
   double        a=1.0,src,dst;

   // Check limits
   x0=X0<0?-X0:0;
   y0=Y0<0?-Y0:0;

   x1=DefPaste->NI+X0>ToDef->NI?ToDef->NI:DefPaste->NI;
   y1=DefPaste->NJ+Y0>ToDef->NJ?ToDef->NJ:DefPaste->NJ;

   // If paste is out of destination
   if (x0>ToDef->NI || x1<0 || y0>ToDef->NJ || y1<0) {
      return(0);
   }

   // Maximum number of band to paste
   nc=fmin(ToDef->NC,DefPaste->NC);
   
   dy=Y0;
   for (y=y0;y<y1;y++) {
      dx=X0;
      for (x=x0;x<x1;x++) {
         idxs=FIDX2D(DefPaste,x,y);
         idxd=FIDX2D(ToDef,dx,dy);
         
         if (DefPaste->NC==4) {
            Def_Get(DefPaste,3,idxs,a);
         }

         for(c=0;c<nc;c++) {               
            Def_Get(DefPaste,c,idxs,src);
            
            if (DefPaste->NC==4) {
               Def_Get(ToDef,c,idxd,dst);
               src=src*a+dst*(1.0-a);
            }
            Def_Set(ToDef,c,idxd,src);
         }
         dx++;
      }
      dy++;
   }
   return(1);
}

int Def_GetValue(TGeoRef *Ref,TDef *Def,int C,double X,double Y,double Z,double *Length,double *ThetaXY){

   double       x,y,d,ddir=0.0;
   int          mem,ix,iy;
   float        valf,valdf;
   void        *p0,*p1;
   unsigned int idx;

  *Length=Def->NoData;
   d=Ref->GRTYP[0]=='W'?1.0:0.5;

   //i on est a l'interieur de la grille ou que l'extrapolation est activee
   if (C<Def->NC && X>=(Ref->X0-d) && Y>=(Ref->Y0-d) && Z>=0 && X<=(Ref->X1+d) && Y<=(Ref->Y1+d) && Z<=Def->NK-1) {

      X-=Ref->X0;
      Y-=Ref->Y0;
      DEFCLAMP(Def,X,Y);

      // Index memoire du niveau desire
      mem=Def->NIJ*(int)Z;
      ix=lrint(X);
      iy=lrint(Y);
      idx=iy*Def->NI+ix;

      // Check for mask
      if (Def->Mask && !Def->Mask[idx]) {
         return(FALSE);
      }
      
      if (Def->Data[1] && !C) { 
         Def_Pointer(Def,0,mem,p0);
         Def_Pointer(Def,1,mem,p1);
         GeoRef_XYWDVal(Ref,&valf,&valdf,p0,p1,&x,&y,1);
         *Length=valf;

         // If it's 3D, use the mode for speed since GeoRef_XYWDVal only uses 2D
         if (Def->Data[2])
            GeoRef_XYVal(Ref,&valf,(float*)&Def->Mode[mem],&x,&y,1);
            *Length=valf;
         if (ThetaXY)
            *ThetaXY=valdf;
      } else {            
         Def_Pointer(Def,C,mem,p0);
         GeoRef_XYVal(Ref,&valf,p0,&x,&y,1);
         *Length=valf;
      }
      return(TRUE);
   }
   
   return(FALSE);
}

static void inline Def_SetValue(TDef *Def,int X, int Y,int Z, double Value,TRef_Combine Comb) {

   unsigned long idx;
   double        val;

   if (FIN2D(Def,X,Y)) {
      idx=Z?FIDX3D(Def,X,Y,Z):FIDX2D(Def,X,Y);

      if (Comb==CB_REPLACE) {
         Def_Set(Def,0,idx,Value);
      } else {
         Def_Get(Def,0,idx,val);
         if (!DEFVALID(Def,val)) val=0.0;

         switch(Comb) {
            case CB_MIN    : if (Value<val) Def_Set(Def,0,idx,Value); break;
            case CB_MAX    : if (Value>val) Def_Set(Def,0,idx,Value); break;
            case CB_AVERAGE: Def->Accum[idx]++;
            case CB_SUM    : Value+=val;  Def_Set(Def,0,idx,Value); break;
            case CB_REPLACE: break;
         }
      }
   }
}

/*----------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_Rasterize>
 * Creation     : Juin 2004 J.P. Gauthier - CMC/CMOE
 *
 * But          : Rasteriser des donnees vectorielles dans une bande ou champs
 *
 * Parametres   :
 *   <Def>      : Definition des donnees raster
 *   <Ref>      : Referentiel des donnnes raster
 *   <Geom>     : Donnees vectorielle a rasteriser
 *   <Value>    : Valuer a assigner
 *   <Comb>    : Mode de combinaison des valeurs multiples (CB_REPLACE,CB_MIN,CB_MAX,CB_SUM,CB_AVERAGE)
 *
 * Retour       :
 *
 * Remarques    :
 *
 *---------------------------------------------------------------------------------------------------------------
*/
int GeoRef_Rasterize(TDef *Def,TGeoRef *Ref,OGRGeometryH Geom,double Value) {

#ifdef HAVE_GDAL
   int    i,j,g,ind1,ind2;
   int    x,y,miny,maxy,minx,maxx;
   int    ints,n,ns,np;
   int   *polyInts;
   double dminy,dmaxy,dx1,dy1,dx2,dy2,dy;
   double intersect,tmpd;
   int    horizontal_x1,horizontal_x2;
   int    dnx,dny,x0,x1,y0,y1,fr,sx,sy;

   OGRGeometryH geom;

   if (!Ref || !Def) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid destination\n",__func__);
      return(0);
   }
   if (!Geom) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid source\n",__func__);
      return(0);
   }

   n=0;
   dminy=1e300;
   dmaxy=-1e300;

   for (i=0;i<OGR_G_GetGeometryCount(Geom);i++) {
      geom=OGR_G_GetGeometryRef(Geom,i);
      if (EQUAL(OGR_G_GetGeometryName(geom),"LINEARRING")) {
         n+=ns=OGR_G_GetPointCount(geom);
         for (j=0;j<ns;j++) {
            dy=OGR_G_GetY(geom,j);
            if (dy<dminy)
               dminy=dy;
            if (dy>dmaxy)
               dmaxy=dy;
         }
      } else {
         GeoRef_Rasterize(Def,Ref,geom,Value);
      }
   }
   if (!(n+=OGR_G_GetPointCount(Geom)))
      return(1);

   switch (OGR_G_GetDimension(Geom)) {
      case 0: // Point type
         for (i=0;i<n;i++) {
            dx1=OGR_G_GetX(Geom,i);
            dy1=OGR_G_GetY(Geom,i);
            x=lrint(dx1);
            y=lrint(dy1);
            Def_SetValue(Def,x,y,0,Value,Ref->Options.Combine);
         }
         break;

      case 1: // Line type
         for (i=1;i<n;i++) {
            ind2=i;
            ind1=i==0?n-1:i-1;

            dx1=OGR_G_GetX(Geom,ind1);
            dy1=OGR_G_GetY(Geom,ind1);

            dx2=OGR_G_GetX(Geom,ind2);
            dy2=OGR_G_GetY(Geom,ind2);

            x0=lrint(dx1); y0=lrint(dy1);
            x1=lrint(dx2); y1=lrint(dy2);
            dny=y1-y0;
            dnx=x1-x0;
            if (dny<0) {
               dny=-dny;
               sy=-1;
            } else {
               sy=1;
            }
            if (dnx<0) {
               dnx=-dnx;
               sx=-1;
            } else {
               sx=1;
            }
            dny<<=1;
            dnx<<=1;

            if (FIN2D(Def,x0,y0))
               Def_Set(Def,0,FIDX2D(Def,x0,y0),Value);
            if (dnx>dny) {
               fr=dny-(dnx>>1);
               while(x0!=x1) {
                  if (fr>=0) {
                     y0+=sy;
                     fr-=dnx;
                  }
                  x0+=sx;
                  fr+=dny;
                  Def_SetValue(Def,x0,y0,0,Value,Ref->Options.Combine);
               }
            } else {
               fr=dnx-(dny>>1);
               while(y0!=y1) {
                  if (fr>=0) {
                     x0+=sx;
                     fr-=dny;
                  }
                  y0+=sy;
                  fr+=dnx;
                  Def_SetValue(Def,x0,y0,0,Value,Ref->Options.Combine);
               }
            }
         }
         break;

      case 2: // Polygon type
         miny=(int)(dminy<0?0:dminy);
         maxy=(int)(dmaxy>=Def->NJ?Def->NJ-1:dmaxy);
         minx=0;
         maxx=Def->NI-1;

         polyInts=(int*)malloc(sizeof(int)*n);

         // Fix in 1.3: count a vertex only once
         for (y=miny;y<=maxy;y++) {
            dy=y; // center height of line
            ints=0 ;

            // Initialize polyInts, otherwise it can sometimes causes a segfault
            for (i=0;i<n;i++) {
               polyInts[i]=-1;
            }

            ns=OGR_G_GetGeometryCount(Geom);
            for (g=0;g<(ns==0?1:ns);g++) {
               if (ns) {
                  geom=OGR_G_GetGeometryRef(Geom,g);
               } else {
                  geom=Geom;
               }
               np=OGR_G_GetPointCount(geom);

               for (i=0;i<np;i++) {
                  ind2=i;
                  ind1=i==0?np-1:i-1;

                  dx1=OGR_G_GetX(geom,ind1);
                  dy1=OGR_G_GetY(geom,ind1);

                  dx2=OGR_G_GetX(geom,ind2);
                  dy2=OGR_G_GetY(geom,ind2);

                  if ((dy1<dy && dy2<dy) || (dy1>dy && dy2>dy))
                     continue;

                  if (dy1<dy2) {
                  } else if (dy1>dy2) {
                     tmpd=dy2;
                     dy2=dy1;
                     dy1=tmpd;
                     tmpd=dx2;
                     dx2=dx1;
                     dx1=tmpd;
                  } else { // if (fabs(dy1-dy2)< 1.e-6)
                           /*AE: DO NOT skip bottom horizontal segments
                           -Fill them separately-
                           They are not taken into account twice.*/
                     if (dx1>dx2) {
                        horizontal_x1=lrint(dx2);
                        horizontal_x2=lrint(dx1);
                        if ((horizontal_x1>maxx) || (horizontal_x2<minx))
                           continue;

                        // fill the horizontal segment (separately from the rest)
                        for(x=horizontal_x1;x<horizontal_x2;x++)
                           Def_SetValue(Def,x,y,0,Value,Ref->Options.Combine);
                        continue;
                     } else {
                        // skip top horizontal segments (they are already filled in the regular loop)
                        continue;
                     }
                  }

                  if ((dy<dy2) && (dy>=dy1)) {
                     intersect=(dy-dy1)*(dx2-dx1)/(dy2-dy1)+dx1;
                     polyInts[ints++]=lrint(intersect);
                  }
               }
            }
            qsort(polyInts,ints,sizeof(int),QSort_Int);

            for (i=0;i<ints;i+=2) {
               if (polyInts[i]<=maxx && polyInts[i+1]>=minx) {
                  for(x=polyInts[i];x<=polyInts[i+1];x++)
                     Def_SetValue(Def,x,y,0,Value,Ref->Options.Combine);
               }
            }
         }
         free(polyInts);
         break;
   }
   return(1);
#else
   Lib_Log(APP_LIBGEOREF,APP_ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
   return(0);
#endif
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_Cell2OGR>
 * Creation     : Mars 2006 J.P. Gauthier - CMC/CMOE
 *
 * But          : Projecter une cellule de grille dans le referentiel d'une autre
 *
 * Parametres   :
 *   <Geom>     : Polygone a initialiser
 *   <ToRef>    : GeoReference destination
 *   <FromRef>  : GeoReference source
 *   <I>        : Coordonnee X
 *   <J>        : Coordonnee Y
 *   <Seg>      : Facteur de sectionnement des segments
 *
 * Retour       :
 *   <Nb>       : Nombre de points (negatif si ca passe le wrap)
 *
 * Remarques    :
 *
 *---------------------------------------------------------------------------------------------------------------
*/
int GeoRef_Cell2OGR(OGRGeometryH Geom,TGeoRef *ToRef,TGeoRef *FromRef,int I,int J,int Seg) {

#ifdef HAVE_GDAL
   double n,dn,df;
   double x0,x1,x,y,la,lo;
   int    pt=0;

   if (!Geom || !ToRef || !FromRef) {
      return(0);
   }

   dn=1.0/Seg;
   df=dn*0.5;

   x0=1e32;
   x1=-1e32;

   // Top Left
   for(n=-0.5;n<(0.5+df);n+=dn) {
      x=I-0.5; y=J+n;
      GeoRef_XY2LL(FromRef,&la,&lo,&x,&y,1,TRUE);
      GeoRef_LL2XY(ToRef,&x,&y,&la,&lo,1,TRUE);
      x0=fmin(x0,x);
      x1=fmax(x1,x);
      OGR_G_SetPoint_2D(Geom,pt++,x,y);
   }

   // Top right
   for(n=-0.5;n<(0.5+df);n+=dn) {
      x=I+n; y=J+0.5;
      GeoRef_XY2LL(FromRef,&la,&lo,&x,&y,1,TRUE);
      GeoRef_LL2XY(ToRef,&x,&y,&la,&lo,1,TRUE);
      x0=fmin(x0,x);
      x1=fmax(x1,x);
      OGR_G_SetPoint_2D(Geom,pt++,x,y);
   }

   // Right bottom
   for(n=0.5;n>-(0.5+df);n-=dn) {
      x=I+0.5; y=J+n;
      GeoRef_XY2LL(FromRef,&la,&lo,&x,&y,1,TRUE);
      GeoRef_LL2XY(ToRef,&x,&y,&la,&lo,1,TRUE);
      x0=fmin(x0,x);
      x1=fmax(x1,x);
      OGR_G_SetPoint_2D(Geom,pt++,x,y);
   }

   // Bottom Left
   for(n=0.5;n>-(0.5+df);n-=dn) {
      x=I+n; y=J-0.5;
      GeoRef_XY2LL(FromRef,&la,&lo,&x,&y,1,TRUE);
      GeoRef_LL2XY(ToRef,&x,&y,&la,&lo,1,TRUE);
      x0=fmin(x0,x);
      x1=fmax(x1,x);
      OGR_G_SetPoint_2D(Geom,pt++,x,y);
   }

   // Close the polygon
   x=I-0.5; y=J-0.5;
   GeoRef_XY2LL(FromRef,&la,&lo,&x,&y,1,TRUE);
   GeoRef_LL2XY(ToRef,&x,&y,&la,&lo,1,TRUE);
   OGR_G_SetPoint_2D(Geom,pt++,x,y);

   // If the cell is outside the destination limits
   if ((x0<ToRef->X0 && x1<ToRef->X0) || (x0>ToRef->X1 && x1>ToRef->X1)) {
      return(0);
   }

   // If the size is larger than half the destination, it has to be a wrap around
   if ((x1-x0)>((ToRef->X1-ToRef->X0)>>1)) {
      return(-pt);
   }
   
   return(pt);
#else
   Lib_Log(APP_LIBGEOREF,APP_ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
   return(0);
#endif
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_InterpQuad>
 * Creation     : Novembre 2004 J.P. Gauthier - CMC/CMOE
 *
 * But          : Importer des donnees dans une bande raster
 *
 * Parametres   :
 *   <Def>      : Definition de la donnee
 *   <Ref>      : Referentiel de la donnee
 *   <Geom>     : Geometrie a tester
 *   <Area>     : Aire de la geometrie
 *   <Value>    : Valeur a assigner
 *   <Mode>     : Mode de rasterization
 *   <Type>     : Type de rasterization
 *   <X0>       : Coin inferieur gauche
 *   <Y0>       : Coin inferieur gauche
 *   <X1>       : Coin superieur droit
 *   <Y1>       : Coin superieur droit
 *   <Z>        : Niveau
 *
 * Retour       : Nombre de point de grille affecté
 *
 * Remarques    :
 *
 *---------------------------------------------------------------------------------------------------------------
*/
static int GeoRef_InterpQuad(TDef *Def,TGeoRef *Ref,OGRGeometryH Poly,OGRGeometryH Geom,char Mode,char Type,double Area,double Value,int X0,int Y0,int X1,int Y1,int Z,float **Index) {

#ifdef HAVE_GDAL
   double        dx,dy,dp=1.0,val=0.0;
   int           x,y,n=0,idx2,na;
   OGRGeometryH  inter=NULL,pick;
   OGREnvelope   envg,envp;

   // Setup the intersecting area
   pick=OGR_G_GetGeometryRef(Poly,0);
   dx=(double)X0-0.5;
   dy=(double)Y0-0.5;
   OGR_G_SetPoint(pick,0,dx,dy,0);
   OGR_G_SetPoint(pick,4,dx,dy,0);
   dx=(double)X0-0.5;
   dy=(double)Y1+0.5;
   OGR_G_SetPoint(pick,1,dx,dy,0);
   dx=(double)X1+0.5;
   dy=(double)Y1+0.5;
   OGR_G_SetPoint(pick,2,dx,dy,0);
   dx=(double)X1+0.5;
   dy=(double)Y0-0.5;
   OGR_G_SetPoint(pick,3,dx,dy,0);

   OGR_G_GetEnvelope(pick,&envp);
   OGR_G_GetEnvelope(Geom,&envg);

   na=(Mode=='C' || Mode=='N' || Mode=='A');

   // Test for intersection
   if ((Area>0.0 || !na) && OGM_Intersect(Geom,Poly,&envg,&envp)) {
//   if ((Area>0.0 || !na) && OGR_G_Intersects(Geom,Poly)) {

      // If this is a single pixel
      if (X0==X1 && Y0==Y1) {

         idx2=FIDX2D(Def,X0,Y0);

         // If we are computing areas
         if (Area>0.0) {
           switch(Type) {
               case 'A':  // Area mode
#ifdef HAVE_GPC
                  inter=OGM_GPCOnOGR(GPC_INT,Geom,Poly);
#else                  
                  inter=OGR_G_Intersection(Geom,Poly);
#endif
                  if (Mode=='C' || Mode=='N') {
                     dp=OGR_G_Area(inter)/Area;
                  } else if (Mode=='A') {
                     dp=OGR_G_Area(inter);
                  }
                  break;

               case 'L':  // Length mode
//                  inter=OGR_G_Intersection(Geom,Poly);
                  inter=OGM_Clip(Geom,Poly);
                  if (Mode=='C' || Mode=='N') {
                     dp=OGM_Length(inter)/Area;
                  } else if (Mode=='A') {
                     dp=OGM_Length(inter);
                  }
                  break;

               case 'P':  // Point mode
                  dp=1.0;
                  break;
            }
            val=Value*dp;
            OGR_G_DestroyGeometry(inter);
         } else {
            val=Value;
         }
         // Are we within
         if (Mode!='W' || OGM_Within(Poly,Geom,&envp,&envg)) {
            //Create thread safe region.
            #pragma omp critical
            {
               Def_SetValue(Def,X0,Y0,Z,val,Ref->Options.Combine);

               if (Mode=='N' && Def->Buffer) {
                  Def->Buffer[idx2]+=dp;
               }
            }

            if (*Index) {
               *((*Index)++)=X0;
               *((*Index)++)=Y0;
               *((*Index)++)=dp;
            }
            n=1;
         }
      } else {
         // Otherwise, refine the search by quad dividing
         x=(X1-X0)>>1;
         y=(Y1-Y0)>>1;

         // If within 1 bloc, parse them all
         if (x==0 || y==0) {
            for (x=X0;x<=X1;x++) {
               for (y=Y0;y<=Y1;y++) {
                  n+=GeoRef_InterpQuad(Def,Ref,Poly,Geom,Mode,Type,Area,Value,x,y,x,y,Z,Index);
               }
            }
         } else {
            n+=GeoRef_InterpQuad(Def,Ref,Poly,Geom,Mode,Type,Area,Value,X0,Y0,X0+x,Y0+y,Z,Index);
            n+=GeoRef_InterpQuad(Def,Ref,Poly,Geom,Mode,Type,Area,Value,X0+x+1,Y0,X1,Y0+y,Z,Index);
            n+=GeoRef_InterpQuad(Def,Ref,Poly,Geom,Mode,Type,Area,Value,X0,Y0+y+1,X0+x,Y1,Z,Index);
            n+=GeoRef_InterpQuad(Def,Ref,Poly,Geom,Mode,Type,Area,Value,X0+x+1,Y0+y+1,X1,Y1,Z,Index);
         }
      }
   }

   return(n);
#else
   Lib_Log(APP_LIBGEOREF,APP_ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
   return(0);
#endif
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_InterpOGR>
 * Creation     : Novembre 2004 J.P. Gauthier - CMC/CMOE
 *
 * But          : Importer des donnees vectorielles dans une bande raster
 *
 * Parametres   :
 *   <Def>      : Definition de la donnee
 *   <Ref>      : Referentiel de la donnee
 *   <Layer>    : Couche a importer
 *   <Field>    : Champs de la couche a utiliser
 *   <Value>    : Valeur a assigner
 *   <Final>    : Finalisation de l'operation (Averaging en plusieurs passe)
 *
 * Retour       :
 *   <OK>       : ERROR=0
 *
 * Remarques    :
 *   The index format is a stream of float32 containing for each OGR feature, the list of gridcells and interpolation
 *   factor separated by -1,
 *      [ogr feature index] [feature value area] [grid cell x] [grid cell y] [factor] [grid cell x] [grid cell y] [factor] ... -1.0 ...
 *---------------------------------------------------------------------------------------------------------------
*/
int GeoRef_InterpOGR(TDef *ToDef,TGeoRef *ToRef,OGR_Layer *Layer,TGeoRef *LayerRef,char *Field,double Value,int Final) {

#ifdef HAVE_GDAL
   TGridSet *gset;
   long     f,n=0,nt=0,idx2;
   double   value,val,area,dp,x0,y0;
   int      fld=-1,pi,pj,error=0,isize=0;
   char     mode,type,*c;
   float   *ip=NULL,*lp=NULL,**index=NULL;
   TCoord   co;
   Vect3d   vr;

   OGRSpatialReferenceH          srs=NULL;
   OGRCoordinateTransformationH  tr=NULL;
   OGRGeometryH                  geom=NULL,utmgeom=NULL,hgeom,pick=NULL,poly=NULL;
   OGREnvelope                   env;

   if (!ToRef || !ToDef) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid destination\n",__func__);
      return(0);
   }
   if (!LayerRef || !Layer) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid source\n",__func__);
      return(0);
   }

   if (!Layer->NFeature) {
      return(1);
   }

   // Recuperer la valeur a utiliser dans l'interpolation
   if (Field) {
      if (strcmp(Field,"FEATURE_AREA")==0) {
         fld=-2;
      } else if (strcmp(Field,"FEATURE_AREA_METER")==0) {
         fld=-3;
      } else if (strcmp(Field,"FEATURE_LENGTH")==0) {
         fld=-4;
      } else if (strcmp(Field,"FEATURE_LENGTH_METER")==0) {
         fld=-5;
      } else if (strcmp(Field,"FEATURE_ID")==0) {
         fld=-6;
      } else if (strcmp(Field,"ZCOORD_MIN")==0) {
         fld=-7;
      } else if (strcmp(Field,"ZCOORD_MAX")==0) {
         fld=-8;
      } else if (strcmp(Field,"ZCOORD_AVG")==0) {
         fld=-9;
      } else {
         fld=OGR_FD_GetFieldIndex(Layer->Def,Field);
         if (fld==-1) {
            Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid layer field\n",__func__);
            return(0);
         }
      }
   }

   // In case of average, we need an accumulator
   if (ToRef->Options.Combine==CB_AVERAGE) {
      if (!ToDef->Accum) {
         ToDef->Accum=malloc(FSIZE2D(ToDef)*sizeof(int));
         if (!ToDef->Accum) {
            Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable to allocate accumulation buffer\n",__func__);
            return(0);
         }
      }
      memset(ToDef->Accum,0x0,FSIZE2D(ToDef)*sizeof(int));
   }

   gset=GeoRef_SetGet(ToRef,LayerRef);

   // Do we have and index
   if (gset->Index && gset->Index[0]!=REF_INDEX_EMPTY) {

      // As long as the file or the list is not empty
      ip=gset->Index;
      while(*ip!=REF_INDEX_END) {

         // Get the gridpoint
         f=*(ip++);
         area=*(ip++);
         value=*(ip++);

         if (f>=Layer->NFeature) {
            Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Wrong index, feature index too high (%i)\n",__func__,f);
            return(0);
         }

         // Get value to distribute
         if (fld>=0) {
            value=OGR_F_GetFieldAsDouble(Layer->Feature[f],fld);
         } else if (fld<-9) {
            value=Value;
         }

         // Get the geometry intersections
         while(*ip!=REF_INDEX_SEPARATOR) {
            pi=*(ip++);
            pj=*(ip++);
            dp=*(ip++);

            idx2=FIDX2D(ToDef,pi,pj);

            Def_Get(ToDef,0,idx2,val);
            if (!DEFVALID(ToDef,val)) {
               val=0.0;
            }

            val+=value*dp;

            if (ToRef->Options.Combine==CB_AVERAGE)  {
               ToDef->Accum[idx2]+=1;
            }
            Def_Set(ToDef,0,idx2,val);
         }
         // Skip separator
         ip++;
      }
   } else {

      // Define the max size of the indexes
      isize=1000;
      if ((c=getenv("INTERP_INDEX_SIZE_HINT"))) {
         isize=atoi(c);
      }

      if (gset->Index && gset->Index[0]==REF_INDEX_EMPTY) {
         if (!(index=(float**)malloc(Layer->NFeature*sizeof(float*)))) {
            Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable to allocate local index arrays\n",__func__);
            return(0);
         }
         ip=gset->Index;
      }
      
      // If the request is in meters
      if (fld==-3 || fld==-5) {
         // Create an UTM referential and transform to convert to meters
         x0=LayerRef->X0;
         y0=LayerRef->Y0;
         LayerRef->XY2LL(LayerRef,&x0,&y0,&co.Lat,&co.Lon,1);
         srs=OSRNewSpatialReference(NULL);
         OSRSetUTM(srs,(int)ceil((180+co.Lon)/6),(int)co.Lat);
         tr=OCTNewCoordinateTransformation(LayerRef->Spatial,srs);

         if (!srs || !tr) {
            Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not initiate UTM transormation\n",__func__);
            return(0);
         }
      }

      // Trouve la feature en intersection
      #pragma omp parallel for private(f,geom,hgeom,utmgeom,env,co,value,vr,n,area,mode,type,lp) firstprivate(pick,poly) shared(Layer,LayerRef,ToRef,fld,tr,error,ip,index,isize) reduction(+:nt)
      for(f=0;f<Layer->NFeature;f++) {
         
         if (index) index[f]=NULL;
         if (error) continue;
         n=0;

         if (Layer->Select[f] && Layer->Feature[f]) {

            geom=utmgeom=NULL;
            // Try to access geometry, skipping instead of failing on bad ones
            if (!(hgeom=OGR_F_GetGeometryRef(Layer->Feature[f]))) {
               Lib_Log(APP_LIBGEOREF,APP_WARNING,"%s: Cannot get handle from geometry: %li\n",__func__,f);
               continue;
            }

            // Copie de la geometrie pour transformation
            if (!(geom=OGR_G_Clone(hgeom))) {
               Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not clone the geometry\n",__func__);
               error=1;
               continue;
            }

            if (!pick) {
               pick=OGR_G_CreateGeometry(wkbLinearRing);
               poly=OGR_G_CreateGeometry(wkbPolygon);
               OGR_G_AddGeometryDirectly(poly,pick);
            }

            // If the request is in meters
            if (fld==-3 || fld==-5) {
               if (!(utmgeom=OGR_G_Clone(geom))) {
                  Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not clone the UTM geomtry\n",__func__);
                  OGR_G_DestroyGeometry(geom);
                  error=1;
                  continue;
               }

               // Transform the geom to utm
               OGR_G_Transform(utmgeom,tr);
            }

            // Get value to distribute
            if (fld>=0) {
               value=OGR_F_GetFieldAsDouble(Layer->Feature[f],fld);
            } else if (fld==-2) {
               value=OGR_G_Area(geom);
            } else if (fld==-3) {
               value=OGR_G_Area(utmgeom);
            } else if (fld==-4) {
               value=OGM_Length(geom);
            } else if (fld==-5) {
               value=OGM_Length(utmgeom);
            } else if (fld==-6) {
               value=f;
            } else if (fld==-7) {
               value=OGM_CoordLimit(geom,2,0);
            } else if (fld==-8) {
               value=OGM_CoordLimit(geom,2,1);
            } else if (fld==-9) {
               value=OGM_CoordLimit(geom,2,2);
            } else {
               value=Value;
            }

            // In centroid mode, just project the coordinate into field and set value
            if (ToRef->Options.InterpDegree==IV_CENTROID) {
               OGM_Centroid2D(geom,&vr[0],&vr[1]);
               GeoRef_XY2LL(LayerRef,&co.Lat,&co.Lon,&vr[0],&vr[1],1,TRUE);
               GeoRef_LL2XY(ToRef,&vr[0],&vr[1],&co.Lat,&co.Lon,1,TRUE);
               n=(int)vr[1]*ToDef->NI+(int)vr[0];
               Def_Set(ToDef,0,n,value);
               nt++;
            } else {

               // Transform geometry to field referential
               OGM_OGRProject(geom,LayerRef,ToRef);

               // Use enveloppe limits to initialize the initial lookup range
               OGR_G_GetEnvelope(geom,&env);
               if (!(env.MaxX<(ToRef->X0-0.5) || env.MinX>(ToRef->X1+0.5) || env.MaxY<(ToRef->Y0-0.5) || env.MinY>(ToRef->Y1+0.5))) {

                  if (ToRef->Options.InterpDegree==IV_FAST) {
                     GeoRef_Rasterize(ToDef,ToRef,geom,value);
                  } else {

                     // Get value to split on
                     area=-1.0;
                     switch(ToRef->Options.InterpDegree) {
                        case IV_FAST                           : break;
                        case IV_WITHIN                         : mode='W'; type='A'; break;
                        case IV_INTERSECT                      : mode='I'; type='A'; break;
                        case IV_CENTROID                       : break;
                        case IV_ALIASED                        : mode='A'; type='A'; area=1.0; break;
                        case IV_CONSERVATIVE                   : mode='C'; type='A'; area=OGR_G_Area(geom); break;
                        case IV_NORMALIZED_CONSERVATIVE        : mode='N'; type='A'; area=OGR_G_Area(geom); break;
                        case IV_LENGTH_CONSERVATIVE            : mode='C'; type='L'; area=OGM_Length(geom); break;
                        case IV_LENGTH_NORMALIZED_CONSERVATIVE : mode='N'; type='L'; area=OGM_Length(geom); break;
                        case IV_LENGTH_ALIASED                 : mode='A'; type='L'; area=OGM_Length(geom); break;
                        case IV_POINT_CONSERVATIVE             : mode='C'; type='P'; area=1.0; break;
                     }

                     // If it's nil then nothing to distribute on
                     if (area>0.0 || ToRef->Options.InterpDegree<=IV_CENTROID) {

                        env.MaxX+=0.5;env.MaxY+=0.5;
                        env.MinX=env.MinX<0?0:env.MinX;
                        env.MinY=env.MinY<0?0:env.MinY;
                        env.MaxX=env.MaxX>(ToRef->X1+0.5)?(ToRef->X1+0.5):env.MaxX;
                        env.MaxY=env.MaxY>(ToRef->Y1+0.5)?(ToRef->Y1+0.5):env.MaxY;
                        area=area<0?0.0:area;

                        // Append feature into index
                        lp=NULL;
                        if (ip) {
                           if (!(index[f]=(float*)malloc(isize*sizeof(float)))) {
                              Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable to allocate local index memory (%i)\n",__func__,f);
                              error=1;
                              continue;
                           } 
                           lp=index[f];
                        }
                        if (lp) {
                           *(lp++)=area;
                           *(lp++)=(fld<0 && fld>=-9)?value:-999.0;
                        }

                        nt+=n=GeoRef_InterpQuad(ToDef,ToRef,poly,geom,mode,type,area,value,env.MinX,env.MinY,env.MaxX,env.MaxY,0,&lp);

                        if (lp) {
                           if (n) {
                              *(lp++)=REF_INDEX_SEPARATOR;
                           } else {
                              index[f][0]=REF_INDEX_EMPTY;  // No intersection found, removed previously inserted feature
                           }
                        }
                        Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: %i hits on feature %i of %i (%.0f %.0f x %.0f %.0f)\n",__func__,n,f,Layer->NFeature,env.MinX,env.MinY,env.MaxX,env.MaxY);
                     }
                  }
               }
            }
            if (geom)    OGR_G_DestroyGeometry(geom);
            if (utmgeom) OGR_G_DestroyGeometry(utmgeom);
         }
      }

      Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: %i total hits\n",__func__,nt);

      // Merge indexes
      n=0;
      if (ip && nt && !error) {
         for(f=0;f<Layer->NFeature;f++) {
            if ((lp=index[f]) && *lp!=REF_INDEX_EMPTY) {
               *(ip++)=f;
               while(*lp!=REF_INDEX_SEPARATOR) {
                  *(ip++)=*(lp++);
               }
               *(ip++)=REF_INDEX_SEPARATOR;
            }
            if (index[f]) free(index[f]);
         }
         *(ip++)=REF_INDEX_END;
         gset->IndexSize=(ip-gset->Index)+1;
         free(index);
      }

      if (tr)
         OCTDestroyCoordinateTransformation(tr);

      if (srs)
         OSRDestroySpatialReference(srs);
   }

   if (ToRef->Options.Combine==CB_AVERAGE) {
      for(n=0;n<FSIZE2D(ToDef);n++) {
         if (ToDef->Accum[n]!=0.0) {
            Def_Get(ToDef,0,n,value);
            value/=ToDef->Accum[n];
            Def_Set(ToDef,0,n,value);
         }
      }
   }

   // Return size of index or number of hits, or 1 if nothing found
   nt=gset->Index?gset->IndexSize:nt;
   return((error || nt==0)?1:nt);
#else
   Lib_Log(APP_LIBGEOREF,APP_ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
   return(0);
#endif
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_InterpSub>
 * Creation     : Aout 2013 J.P. Gauthier - CMC/CMOE
 *
 * But          : Interpolate at sub grid resolution
 *
 * Parametres   :
 *   <ToRef>    : Reference du champs destination
 *   <ToDef>    : Description du champs destination
 *   <FromRef>  : Reference du champs source
 *   <FromDef>  : Description du champs source
 *
 * Retour       :
 *  <OK>       : ERROR=0
 *
 * Remarques    :
 *    Interpolates at sub-grid resolution into a temporary "Sub" buffer for use by other functions
 *    like sub grid variance.
 *---------------------------------------------------------------------------------------------------------------
*/
int GeoRef_InterpSub(TGeoRef *ToRef,TDef *ToDef,TGeoRef *FromRef,TDef *FromDef) {

   int  idx,i,j,x,y,x0,x1,y0,y1,s;
   double val,val1,di,dj,d,la,lo;

   if (!ToRef || !ToDef) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid destination\n",__func__);
      return(0);
   }
   if (!FromRef || !FromDef) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid source\n",__func__);
      return(0);
   }

   if (ToDef->SubSample<=1) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Sub sampling is not defined\n",__func__);
      return(0);
   }

   if (!GeoRef_Intersect(FromRef,ToRef,&x0,&y0,&x1,&y1,0)) {
      return(1);
   }

   // Allocate and initialize subgrid
   s=ToDef->SubSample*ToDef->SubSample;
   if (!ToDef->Sub) {
      if ((ToDef->Sub=(float*)malloc(ToDef->NI*ToDef->NJ*s*sizeof(float)))) {
         for(idx=0;idx<ToDef->NI*ToDef->NJ*s;idx++) {
            ToDef->Sub[idx]=ToDef->NoData;
         }
      } else {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable ot allocate subgrid\n",__func__);
         return(0);
      }  
   }

   d=1.0/(ToDef->SubSample-1);

   // Interpolate values on subgrid
   for(y=y0;y<=y1;y++) {
      for(x=x0;x<=x1;x++) {
         idx=(y*ToDef->NI+x)*s;

         for(j=0;j<ToDef->SubSample;j++) {
            for(i=0;i<ToDef->SubSample;i++,idx++) {

               di=(double)x-0.5+d*i;
               dj=(double)y-0.5+d*j;

               GeoRef_XY2LL(FromRef,&la,&lo,&di,&dj,1,TRUE);
               if (GeoRef_LL2XY(ToRef,&di,&dj,&la,&lo,1,FALSE)) {
                  if (Def_GetValue(FromRef,FromDef,0,di,dj,FromDef->Level,&val,&val1) && DEFVALID(FromDef,val)) {
                     ToDef->Sub[idx]=val;
                  }
               }
            }
         }
      }
   }

   return(1);
}

/*----------------------------------------------------------------------------
 * Nom      : <GeoRef_InterpConservative>
 * Creation : Mai 2006 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Effectue l'interpolation conservative
 *
 * Parametres  :
 *  <ToRef>    : Reference du champs destination
 *  <ToDef>    : Description du champs destination
 *  <FromRef>  : Reference du champs source
 *  <FromDef>  : Description du champs source
 *  <Final>    : Finalisation de l'operation (Averaging en plusieurs passe)
 *
 * Retour:
 *  <OK>       : ERROR=0
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
int GeoRef_InterpConservative(TGeoRef *ToRef,TDef *ToDef,TGeoRef *FromRef,TDef *FromDef,int Final) {

#ifdef HAVE_GDAL
   TGridSet    *gset=NULL;
   int          i,j,n,na,nt=0,p=0,pi,pj,idx2,idx3,wrap,k=0,isize,nidx,error=0;
   char        *c;
   double       val0,val1,area,x,y,z,dp;
   float       *ip=NULL,*lp=NULL,**index=NULL;
   OGRGeometryH cell=NULL,ring=NULL,*pick=NULL,*poly=NULL;
   OGREnvelope  env;

   if (!ToRef || !ToDef) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid destination\n",__func__);
      return(0);
   }
   if (!FromRef || !FromDef) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid source\n",__func__);
      return(0);
   }

   // Allocate area buffer if needed
   if (ToRef->Options.InterpDegree==IR_NORMALIZED_CONSERVATIVE && !ToDef->Buffer) {
      ToDef->Buffer=(double*)malloc(FSIZE2D(ToDef)*sizeof(double));
      if (!ToDef->Buffer) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable to allocate area buffer\n",__func__);
         return(0);
      }
   }
   gset=GeoRef_SetGet(ToRef,FromRef);

   // Define the max size of the indexes
   isize=100;
   if ((c=getenv("INTERP_INDEX_SIZE_HINT"))) {
      isize=atoi(c);
   }

   // Process one level at a time
   for (k=0;k<ToDef->NK;k++) {

      if (ToDef->Buffer)
         for(i=0;i<FSIZE2D(ToDef);i++) ToDef->Buffer[i]=0.0;

      // Do we have and index
      if (gset->Index && gset->Index[0]!=REF_INDEX_EMPTY) {

         // As long as the file or the list is not empty
         ip=gset->Index;
         while(*ip!=REF_INDEX_END) {

            // Get the gridpoint
            i=*(ip++);
            j=*(ip++);

            if (!FIN2D(FromDef,i,j)) {
               Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Wrong index, index coordinates (%i,%i)\n",__func__,i,j);
               return(0);
            }

            // Get this gridpoint value
            Def_Get(FromDef,0,FIDX3D(FromDef,i,j,k),val1);
            if (!DEFVALID(FromDef,val1)) {
               continue;
            }

            // Get the geometry intersections
            while(*ip!=REF_INDEX_SEPARATOR) {
               pi=*(ip++);
               pj=*(ip++);
               dp=*(ip++);

               idx2=FIDX2D(ToDef,pi,pj);
               idx3=FIDX3D(ToDef,pi,pj,k);
               Def_Get(ToDef,0,idx3,val0);
               if (!DEFVALID(ToDef,val0))
                  val0=0.0;

               // Assign new value
               val0+=val1*dp;
               Def_Set(ToDef,0,idx3,val0);

               if (ToRef->Options.InterpDegree==IR_NORMALIZED_CONSERVATIVE) {
                  ToDef->Buffer[idx2]+=dp;
               }
            }
            // Skip separator
            ip++;
         }
      } else {
         int cnt, intersect;

         if (gset->Index && gset->Index[0]==REF_INDEX_EMPTY) {
            if (!(index=(float**)malloc(FSIZE2D(FromDef)*sizeof(float*)))) {
               Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable to allocate local index arrays\n",__func__);
               return(0);
            }
            ip=gset->Index;
         }

         // Damn, we dont have the index, do the long run
         #pragma omp parallel for collapse(2) firstprivate(cell,ring,pick,poly) private(i,j,nidx,wrap,intersect,cnt,p,x,y,z,lp,area,val1,env,n,na) shared(k,FromRef,FromDef,ToRef,ToDef,error) reduction(+:nt)
         for(j=0;j<FromDef->NJ;j++) {
            for(i=0;i<FromDef->NI;i++) {

               nidx=j*FromDef->NI+i;
               if (index) index[nidx]=NULL;
               if (error) continue;

               // Create gridcell geometry object
               if (!cell) {
                  cell=OGR_G_CreateGeometry(wkbPolygon);
                  ring=OGR_G_CreateGeometry(wkbLinearRing);
                  OGR_G_AddGeometryDirectly(cell,ring);
               }

               if (!pick) {
                  pick=OGR_G_CreateGeometry(wkbLinearRing);
                  poly=OGR_G_CreateGeometry(wkbPolygon);
                  OGR_G_AddGeometryDirectly(poly,pick);
               }

               // Project the source gridcell into the destination
               wrap=GeoRef_Cell2OGR(ring,ToRef,FromRef,i,j,ToRef->Options.Segment);
               intersect=0;
               cnt = OGR_G_GetPointCount(ring);
               for(p=0;p<cnt;p++) {
                  OGR_G_GetPoint(ring,p,&x,&y,&z);
                  if (x>=(ToRef->X0-0.5) && x<=(ToRef->X1+0.5) && y>=(ToRef->Y0-0.5) && x<=(ToRef->X1+0.5)) {
                     intersect=1;
                     break;
                  }
               }

               if (!wrap || !intersect)
                  continue;

               // Allocate local index
               lp=NULL;
               if (ip) {
                  if (!(index[nidx]=(float*)malloc(isize*sizeof(float)))) {
                     Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable to allocate local index memory\n",__func__);
                     error=1;
                     continue;
                  } 
                  index[nidx][0]=REF_INDEX_EMPTY;
                  lp=index[nidx];
               }
               n=na=0;

               // Are we crossing the wrap around
               if (wrap<0) {
                  // If so, move the wrapped points (assumed greater than NI/2) to the other side
                  for(p=0;p<-wrap;p++) {
                     OGR_G_GetPoint(ring,p,&x,&y,&z);
                     if (x>ToDef->NI>>1) {
                        x-=ToDef->NI;
                        OGR_G_SetPoint_2D(ring,p,x,y);
                     }
                  }

                  // Process the cell
                  area=OGR_G_Area(cell);

                  if (area>0.0) {

                    Def_Get(FromDef,0,FIDX3D(FromDef,i,j,k),val1);
                     // If we are saving the indexes, we have to process even if nodata but use 0.0 so as to not affect results
                     if (!DEFVALID(FromDef,val1)) {
                        if (lp) {
                           val1=0.0;
                        } else {
                           continue;
                        }
                     }

                     // Use enveloppe limits to initialize the initial lookup range
                     OGR_G_GetEnvelope(ring,&env);
                     env.MaxX+=0.5;env.MaxY+=0.5;
                     env.MinX=env.MinX<0?0:env.MinX;
                     env.MinY=env.MinY<0?0:env.MinY;
                     env.MaxX=env.MaxX>ToRef->X1?ToRef->X1:env.MaxX;
                     env.MaxY=env.MaxY>ToRef->Y1?ToRef->Y1:env.MaxY;

                     nt+=na=GeoRef_InterpQuad(ToDef,ToRef,poly,cell,ToRef->Options.InterpDegree==IR_CONSERVATIVE?'C':'N','A',area,val1,env.MinX,env.MinY,env.MaxX,env.MaxY,k,&lp);

                     Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: %i hits on grid point %i %i (%.0f %.0f x %.0f %.0f)\n",__func__,na,i,j,env.MinX,env.MinY,env.MaxX,env.MaxY);
                  }

                  // We have to process the part that was out of the grid limits so translate everything NI points
                  for(p=0;p<-wrap;p++) {
                     OGR_G_GetPoint(ring,p,&x,&y,&z);
                     x+=ToDef->NI;
                     OGR_G_SetPoint_2D(ring,p,x,y);
                  }
               }

              area=OGR_G_Area(cell);

               if (area>0.0) {
                  Def_Get(FromDef,0,FIDX3D(FromDef,i,j,k),val1);
                  if (!DEFVALID(FromDef,val1)) {
                     // If we are saving the indexes, we have to process even if nodata but use 0.0 so as to not affect results
                     if (lp) {
                        val1=0.0;
                     } else {
                        continue;
                     }
                  }

                  // Use enveloppe limits to initialize the initial lookup range
                  OGR_G_GetEnvelope(ring,&env);
                  if (!(env.MaxX<ToRef->X0 || env.MinX>ToRef->X1 || env.MaxY<ToRef->Y0 || env.MinY>ToRef->Y1)) {
                     env.MaxX+=0.5;env.MaxY+=0.5;
                     env.MinX=env.MinX<0?0:env.MinX;
                     env.MinY=env.MinY<0?0:env.MinY;
                     env.MaxX=env.MaxX>ToRef->X1?ToRef->X1:env.MaxX;
                     env.MaxY=env.MaxY>ToRef->Y1?ToRef->Y1:env.MaxY;

                     nt+=n=GeoRef_InterpQuad(ToDef,ToRef,poly,cell,ToRef->Options.InterpDegree=IR_CONSERVATIVE?'C':'N','A',area,val1,env.MinX,env.MinY,env.MaxX,env.MaxY,k,&lp);

                     Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: %i hits on grid point %i %i (%.0f %.0f x %.0f %.0f)\n",__func__,n,i,j,env.MinX,env.MinY,env.MaxX,env.MaxY);
                  }
               }
               if (lp && (n || na)) {
                  *(lp++)=REF_INDEX_SEPARATOR; // End the list for this gridpoint
               }
            }
         }

         // Merge indexes
         n=0;
         if (ip && nt && !error) {
            for(j=0;j<FromDef->NJ;j++) {
               for(i=0;i<FromDef->NI;i++) {
                  nidx=j*FromDef->NI+i;

                  if ((lp=index[nidx]) && *lp!=REF_INDEX_EMPTY) {
                     // Append gridpoint to the index
                     *(ip++)=i;
                     *(ip++)=j;

                     // This gridpoint wraps around
                     while(*lp!=REF_INDEX_SEPARATOR) {
                        *(ip++)=*(lp++);
                     }
                     *(ip++)=REF_INDEX_SEPARATOR;
                  }
                  if (index[nidx]) free(index[nidx]);
               }
            }
            *(ip++)=REF_INDEX_END;
            gset->IndexSize=(ip-gset->Index)+1;
            free(index);
         }
         Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: %i total hits\n",__func__,nt);
      }

      // Finalize and reassign
      idx3=FSIZE2D(ToDef)*k;
      if (Final && GeoRef_Options.InterpDegree==IR_NORMALIZED_CONSERVATIVE) {
         for(n=0;n<FSIZE2D(ToDef);n++) {
            if (ToDef->Buffer[n]!=0.0) {
               Def_Get(ToDef,0,idx3+n,val0);
               val0/=ToDef->Buffer[n];
               Def_Set(ToDef,0,idx3+n,val0);
               ToDef->Buffer[n]=0.0;
            }
         }
      }

   }

   OGR_G_DestroyGeometry(ring);
   OGR_G_DestroyGeometry(cell);

   // Return size of index or number of hits, or 1 if nothing found
   nt=gset->Index?gset->IndexSize:nt;
   return(nt==0?1:nt);
#else
   Lib_Log(APP_LIBGEOREF,APP_ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
   return(0);
#endif

}

/*----------------------------------------------------------------------------
 * Nom      : <GeoRef_InterpAverage>
 * Creation : Mai 2006 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Effectue l'interpolation par moyennage minimum ou maximum
 *
 * Parametres  :
 *  <ToRef>    : Reference du champs destination
 *  <ToDef>    : Description du champs destination
 *  <FromRef>  : Reference du champs source
 *  <FromDef>  : Description du champs source
 *  <Table>    : Table de donnees a verifier
 *  <TmpDef>   : Description du champs de precalcul (ex: pour VARIANCE, moyenne)
 *  <Final>    : Finalisation de l'operation (Averaging en plusieurs passe)
 *
 * Retour:
 *  <OK>       : ERROR=0
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
int GeoRef_InterpAverage(TGeoRef *ToRef,TDef *ToDef,TGeoRef *FromRef,TDef *FromDef,double *Table,TDef **lutDef, int lutSize, TDef *TmpDef,int Final){

   double        val,vx,di[4],dj[4],*fld,*aux,di0,di1,dj0,dj1,ax,ay;
   int          *acc=NULL,x0,x1,y,y0,y1;
   unsigned long idxt,idxk,idxj,n,nijk,nij;
   unsigned int  n2,ndi,ndj,k,t,s,x,dx,dy;
   TGeoScan      gscan;

   if (!ToRef || !ToDef) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid destination\n",__func__);
      return(0);
   }
   if ((!FromRef || !FromDef) && (ToRef->Options.InterpDegree!=IR_NOP && ToRef->Options.InterpDegree!=IR_ACCUM && ToRef->Options.InterpDegree!=IR_BUFFER)) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid source\n",__func__);
      return(0);
   }

   acc=ToDef->Accum;
   fld=ToDef->Buffer;
   aux=ToDef->Aux;
   nij=FSIZE2D(ToDef);
   nijk=FSIZE3D(ToDef);
   val=vx=0.0;

   if (ToRef->Options.InterpDegree!=IR_NOP && ToRef->Options.InterpDegree!=IR_ACCUM && ToRef->Options.InterpDegree!=IR_BUFFER) {
      if (!GeoRef_Intersect(ToRef,FromRef,&x0,&y0,&x1,&y1,0)) {
         return(1);
      }

      // In case of average, we need an accumulator
      if (ToRef->Options.InterpDegree==IR_AVERAGE || ToRef->Options.InterpDegree==IR_VECTOR_AVERAGE || ToRef->Options.InterpDegree==IR_VARIANCE || 
          ToRef->Options.InterpDegree==IR_SQUARE || ToRef->Options.InterpDegree==IR_NORMALIZED_COUNT || ToRef->Options.InterpDegree==IR_COUNT) {
         if (!ToDef->Accum) {
            acc=ToDef->Accum=malloc(nij*sizeof(int));
            if (!ToDef->Accum) {
               Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable to allocate accumulation buffer\n",__func__);
               return(0);
            }
            for(n=0;n<nij;n++) acc[n]=0;
         }
      }

      if (!ToDef->Buffer) {
         fld=ToDef->Buffer=malloc(nijk*sizeof(double));
         if (!ToDef->Buffer) {
            Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable to allocate buffer\n",__func__);
            return(0);
         }
         for(n=0;n<nijk;n++) fld[n]=ToDef->NoData;

        if (ToRef->Options.InterpDegree==IR_VECTOR_AVERAGE) {
            aux=ToDef->Aux=malloc(nijk*sizeof(double));
            if (!ToDef->Aux) {
               Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable to allocate auxiliary buffer\n",__func__);
               return(0);
            }
            for(n=0;n<nijk;n++) aux[n]=0.0;
         }
      }

      if (ToRef->GRTYP[0]=='Y') {
         // Point cloud interpolations
         for(idxt=0;idxt<nij;idxt++) {
            ax=ToRef->AX[idxt];ay=ToRef->AY[idxt];
            if (GeoRef_LL2XY(FromRef,&di0,&dj0,&ax,&ay,1,FALSE)) {
               di0=floor(di0);
               dj0=floor(dj0);
               Def_GetValue(FromRef,FromDef,0,di0,dj0,FromDef->Level,&di[0],&vx);
               Def_GetValue(FromRef,FromDef,0,di0+1.0,dj0,FromDef->Level,&di[1],&vx);
               Def_GetValue(FromRef,FromDef,0,di0,dj0+1.0,FromDef->Level,&di[2],&vx);
               Def_GetValue(FromRef,FromDef,0,di0+1.0,dj0+1.0,FromDef->Level,&di[3],&vx);
            }
            for(s=0;s<4;s++) {
               vx=di[s];
               if (!DEFVALID(FromDef,vx))
                  continue;

               // If the previous value is nodata, initialize the counter
               if (!DEFVALID(FromDef,fld[idxt])) {
                  fld[idxt]=(ToRef->Options.InterpDegree==IR_SUM || ToRef->Options.InterpDegree==IR_AVERAGE|| ToRef->Options.InterpDegree==IR_VECTOR_AVERAGE)?0.0:(ToRef->Options.InterpDegree==IR_MAXIMUM?-HUGE_VAL:HUGE_VAL);
                  if (aux) aux[idxt]=fld[idxt];
               }
               switch(ToRef->Options.InterpDegree) {
                  case IR_MAXIMUM          : if (vx>fld[idxt]) fld[idxt]=vx; break;
                  case IR_MINIMUM          : if (vx<fld[idxt]) fld[idxt]=vx; break;
                  case IR_SUM              : fld[idxt]+=vx;                  break;
                  case IR_AVERAGE          : fld[idxt]+=vx; acc[idxt]++;     break;
                  case IR_VECTOR_AVERAGE   : vx=DEG2RAD(vx); fld[idxt]+=cos(vx); aux[idxt]+=sin(vx); acc[idxt]++; break;
                  default:
                     Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid interpolation type\n",__func__);
                     return(0);
               }
            }
         }
      } else {
         int   *fromClass = NULL;
         int   *toClass = NULL;
         int    nbClass = 0;
         int   *rpnClass = NULL;
         int    i;

         if (lutDef && ToDef->NK > 0) {
            Def_Get(lutDef[0],0,0,val);
            if (val != ToDef->NK) {
               Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid LUT class size (%d) vs Array size (%d)  \n",__func__,(int)val,ToDef->NK);
               return(0);
            }
            // first row is the class count and idno 
            nbClass = lutDef[0]->NI - 1;       
            fromClass = (int *)malloc( sizeof(int)*nbClass );
            for (i = 0; i < nbClass ; i++)
               {
               Def_Get(lutDef[0],0,i+1,val);
               fromClass[i] = val;
               }
            if (lutSize == 2) { // only FROM and TO
               toClass = (int *)malloc( sizeof(int)*nbClass );
               for (i = 0; i < nbClass ; i++) {
                  Def_Get(lutDef[1],0,i+1,val);
                  // make sure to value is within range
                  if ((val >= 1) && (val <= ToDef->NK))
                     toClass[i] = val;
                  else
                     toClass[i] = -1;
               }
            } else {
               double val;
               rpnClass = (int *)malloc( sizeof(int)*lutSize );
               for (i = 1; i < lutSize ; i++) { // skip 1st column : the class id 
                  Def_Get(lutDef[i],0,0,val);
                  rpnClass[i] = val;
                  if ((val < 1)||(val > ToDef->NK)) {
                     Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid LUT index(%d)=%f vs Array size (%d)  \n",__func__,i,val,ToDef->NK);
                     return(0);
                  }
               }
            }
         }

         // grid based interpolations
         GeoScan_Init(&gscan);

         // if > 2048x2048, loop by lines otherwise, do it in one shot
         n2=ToDef->NI>>1;
         dy=((y1-y0)*(x1-x0))>4194304?0:(y1-y0);
         for(y=y0;y<=y1;y+=(dy+1)) {
            if (!(s=GeoScan_Get(&gscan,ToRef,NULL,FromRef,FromDef,x0,y,x1,y+dy,FromDef->CellDim))) {
               Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Unable to allocate coordinate scanning buffer\n",__func__);
               return(0);
            }

            // Loop over source data
            dx=0;
            for(x=0,n=0;x<gscan.N;x++,n++) {

               // Check if we need to skip last x since we change row and last one is end of a cell
               if (s>1 && dx==gscan.DX) {
                  n++;
                  dx=0;
               }
               dx++;

              // Skip if no mask
               if (FromDef->Mask && !FromDef->Mask[gscan.V[x]])
                  continue;

               // Skip if no data
               Def_Get(FromDef,0,gscan.V[x],vx);
               if (!DEFVALID(FromDef,vx) && ToRef->Options.InterpDegree!=IR_COUNT)
                  continue;

               // Figure out ordered coverage
               if (s>1) {
                  di[0]=gscan.X[n];
                  dj[0]=gscan.Y[n];
                  di[1]=gscan.X[n+1]; di0=fmin(di[0],di[1]); di1=fmax(di[0],di[1]);
                  dj[1]=gscan.Y[n+1]; dj0=fmin(dj[0],dj[1]); dj1=fmax(dj[0],dj[1]);

                  di[2]=gscan.X[n+gscan.DX+1]; di0=fmin(di0,di[2]); di1=fmax(di1,di[2]);
                  dj[2]=gscan.Y[n+gscan.DX+1]; dj0=fmin(dj0,dj[2]); dj1=fmax(dj1,dj[2]);
                  di[3]=gscan.X[n+gscan.DX+2]; di0=fmin(di0,di[3]); di1=fmax(di1,di[3]);
                  dj[3]=gscan.Y[n+gscan.DX+2]; dj0=fmin(dj0,dj[3]); dj1=fmax(dj1,dj[3]);

                  di0=ROUND(di0);dj0=ROUND(dj0);
                  di1=ROUND(di1);dj1=ROUND(dj1);
               } else {
                  di0=di1=ROUND(gscan.X[n]);
                  dj0=dj1=ROUND(gscan.Y[n]);
               }

               // Are we within the destination field
               if (di0>=ToDef->NI || dj0>=ToDef->NJ || di1<0 || dj1<0)
                  continue;

               // Test for polar outsidness (Problem we had with yinyang grids)
               if ((di0<0 && di1>ToDef->NI) || (dj0<0 && dj1>ToDef->NJ))
                  continue;

               // Clamp the coordinates
               if (di0<0) di0=0;
               if (dj0<0) dj0=0;
               if (di1>ToDef->NI-1) di1=ToDef->NI-1;
               if (dj1>ToDef->NJ-1) dj1=ToDef->NJ-1;

               // Are we crossing the wrap around
               if (ToRef->Type&GRID_WRAP && di0<n2 && di1>n2 && (di1-di0)>n2) {
                  val=di0;
                  di0=di1;
                  di1=val+ToDef->NI;
               }

               for(ndj=dj0;ndj<=dj1;ndj++) {
                  idxj=ndj*ToDef->NI;
                  for(ndi=di0;ndi<=di1;ndi++) {
                     idxt=idxj+(ndi>=ToDef->NI?ndi-ToDef->NI:ndi);

                     // Skip if no mask
                     if (!ToDef->Mask || ToDef->Mask[idxt]) {

                        // If the previous value is nodata, initialize the counter
                        if (!DEFVALID(ToDef,fld[idxt])) {
                           fld[idxt]=(ToRef->Options.InterpDegree==IR_SUM || ToRef->Options.InterpDegree==IR_AVERAGE  || ToRef->Options.InterpDegree==IR_VECTOR_AVERAGE || ToRef->Options.InterpDegree==IR_VARIANCE || ToRef->Options.InterpDegree==IR_SQUARE || ToRef->Options.InterpDegree==IR_NORMALIZED_COUNT || ToRef->Options.InterpDegree==IR_COUNT)?0.0:(ToRef->Options.InterpDegree==IR_MAXIMUM?-HUGE_VAL:HUGE_VAL);
                           if (aux) aux[idxt]=fld[idxt];
                        }

                        switch(ToRef->Options.InterpDegree) {
                           case IR_MAXIMUM          : if (vx>fld[idxt]) fld[idxt]=vx;
                                                      break;
                           case IR_MINIMUM          : if (vx<fld[idxt]) fld[idxt]=vx;
                                                      break;
                           case IR_SUM              : fld[idxt]+=vx;
                                                      break;
                           case IR_VARIANCE         : acc[idxt]++;
                                                      Def_Get(TmpDef,0,idxt,val);
                                                      fld[idxt]+=(vx-val)*(vx-val);
                                                      break;
                           case IR_SQUARE           : acc[idxt]++;
                                                      fld[idxt]+=vx*vx;
                                                      break;
                           case IR_COUNT            : acc[idxt]++;
                           case IR_AVERAGE          :
                           case IR_NORMALIZED_COUNT : if (Table) {
                                                         t=0;
                                                         while(t<ToDef->NK) {
                                                            if (vx==Table[t]) {
                                                               if (ToRef->Options.InterpDegree!=IR_COUNT) acc[idxt]++;
                                                               fld[t*nij+idxt]+=1.0;
                                                               break;
                                                            }
                                                            t++;
                                                         }
                                                      } else if (lutDef) {
                                                         int hasvalue=0;
                                                         t=0;
                                                         while(t<nbClass) {
                                                            if (vx==fromClass[t]) {
                                                               if (toClass) { // toClass is between 1 and 26 
                                                                  if (toClass[t] > 0)
                                                                     {
                                                                     fld[(toClass[t]-1)*nij+idxt]+=1.0;
                                                                     hasvalue = 1;
                                                                     }
                                                               } else {
                                                                  for (i = 1; i < lutSize ; i++) { // skip 1st column : the class id 
                                                                     Def_Get(lutDef[i],0,t+1,val);
                                                                     if (val > 0) {
                                                                        fld[(rpnClass[i]-1)*nij+idxt]+=val;
                                                                        hasvalue = 1;
                                                                     }
                                                                  }
                                                               }
                                                               break;
                                                            }
                                                            t++;
                                                         }
                                                         // Dont count missing values
                                                         if (hasvalue)
                                                            if (ToRef->Options.InterpDegree!=IR_COUNT) acc[idxt]++; 
                                                      } else {
                                                         if (DEFVALID(FromDef,vx)) {
                                                            fld[idxt]+=vx;
                                                            if (ToRef->Options.InterpDegree!=IR_COUNT) acc[idxt]++;
                                                         }
                                                      }
                                                      break;
                           case IR_VECTOR_AVERAGE   : vx=DEG2RAD(vx); fld[idxt]+=cos(vx); aux[idxt]+=sin(vx); acc[idxt]++; break;
                           default:
                              Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid interpolation type\n",__func__);
                              return(0);
                        }
                     }
                  }
               }
            }
         }

         GeoScan_Clear(&gscan);
         if (fromClass) free(fromClass);
         if (toClass) free(toClass);
         if (rpnClass) free(rpnClass);
      }
   }

   // Finalize and reassign
   if (Final || ToRef->Options.InterpDegree==IR_ACCUM || ToRef->Options.InterpDegree==IR_BUFFER) {
      idxk=0;
      for(k=0;k<ToDef->NK;k++) {
         for(x=0;x<nij;x++,idxk++) {
            vx=ToDef->NoData;

            switch(ToRef->Options.InterpDegree) {
               case IR_ACCUM:
                  if (acc) vx=acc[x];
                  break;

               case IR_BUFFER:
                  if (fld) vx=fld[idxk];
                  break;

               default:
                  if (fld) {
                     if (DEFVALID(ToDef,fld[idxk])) {
                        if (aux) {
                           if (acc && acc[x]!=0) {
                              fld[idxk]/=acc[x];
                              aux[idxk]/=acc[x];
                           }
                           vx=RAD2DEG(atan2(aux[idxk],fld[idxk]));
                           if (vx<0) vx+=360;
                        } else {
                           vx=fld[idxk];
                           if (acc && acc[x]!=0) {
                              vx=fld[idxk]/acc[x];
                           } else {
                              vx=fld[idxk];
                           }
                        }
                     }
                  }
            }

            Def_Set(ToDef,0,idxk,vx);
         }

         // Copy first column to last if it's repeated
         if (ToRef->Type&GRID_REPEAT) {
            idxt=k*nij;
            for(y=ToRef->Y0;y<=ToRef->Y1;y++,idxt+=ToDef->NI) {
               Def_Get(ToDef,0,idxt,vx);
               Def_Set(ToDef,0,idxt+ToDef->NI-1,vx);
            }
         }
      }
   }
   return(1);
}

/*----------------------------------------------------------------------------
 * @brief  Changes a georeference type to W
 * @date   June 2014
 *    @param[in]     Ref          GeoRef pointer
 *    @param[in]     String       WKT description
 *    @param[in]     Transform    Transformation matrix [2][6]
 *    @param[in]     InvTransform Inverse transformation matrix [2][6]
 *    @param[in]     Spatial      Spatial reference (GDAL object)

 *    @return        GeoRef object (NULL=Error)
*/
int GeoRef_InterpDef(TGeoRef *ToRef,TGeoRef *FromRef,TDef *ToDef,TDef *FromDef,int Final) {

   void *pf0,*pt0,*pf1,*pt1;
   int   k,code=FALSE;

   if (!ToRef || !ToDef) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid destination\n",__func__);
      return(FALSE);
   }
   if (!FromRef || !FromDef) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid source\n",__func__);
      return(FALSE);
   }

   // Same GeoRef
   if (GeoRef_Options.InterpDegree<IR_NOP && GeoRef_Equal(ToRef,FromRef)) {
      Def_CopyData(ToDef,FromDef);
      return(TRUE);
   }

   switch(ToRef->Options.InterpDegree) {
      case IR_NEAREST:
      case IR_LINEAR:
      case IR_CUBIC:
         // Loop on vertical levels
         for(k=0;k<ToDef->NK;k++) {
            Def_Pointer(ToDef,0,k*FSIZE2D(ToDef),pt0);
            Def_Pointer(FromDef,0,k*FSIZE2D(FromDef),pf0);

            if (ToDef->Data[1]) {
               // Interpolation vectorielle
               Def_Pointer(ToDef,1,k*FSIZE2D(ToDef),pt1);
               Def_Pointer(FromDef,1,k*FSIZE2D(FromDef),pf1);

               // In case of Y grid, get the speed and dir instead of wind components
               // since grid oriented components dont mean much
               if (ToRef->GRTYP[0]=='Y') {
                  code=GeoRef_InterpWD(ToRef,FromRef,pt0,pt1,pf0,pf1);
               } else {
                  code=GeoRef_InterpUV(ToRef,FromRef,pt0,pt1,pf0,pf1);
               }
            } else{
               // Interpolation scalaire
               code=GeoRef_Interp(ToRef,FromRef,pt0,pf0);
            }

            // Interpolate mask if need be
            if (FromDef->Mask && ToDef->Mask) {
               code=GeoRef_InterpMask(ToRef,FromRef,&ToDef->Mask[k*FSIZE2D(ToDef)],&FromDef->Mask[k*FSIZE2D(FromDef)]);
            }
         }
         break;

      case IR_CONSERVATIVE:
      case IR_NORMALIZED_CONSERVATIVE:
         // Check compatibility between source and destination
//         if (!Def_Compat(field0->Def,field1->Def)) {
//            field0->GRef=GeoRef_Find(GeoRef_Resize(field0->GRef,field0->Def->NI,field0->Def->NJ));
//         }
         code=GeoRef_InterpConservative(ToRef,ToDef,FromRef,FromDef,Final);
         break;

      case IR_MAXIMUM:
      case IR_MINIMUM:
      case IR_SUM:
      case IR_AVERAGE:
      case IR_VARIANCE:
      case IR_SQUARE:
      case IR_NORMALIZED_COUNT:
      case IR_COUNT:
      case IR_VECTOR_AVERAGE:
      case IR_NOP: 
      case IR_ACCUM:
      case IR_BUFFER:
         //ToDo: add double *Table,TDef **lutDef, int lutSize, TDef *TmpDef
         code=GeoRef_InterpAverage(ToRef,ToDef,FromRef,FromDef,NULL,NULL,0,NULL,Final);
     break;

      case IR_SUBNEAREST:
      case IR_SUBLINEAR:
         // Is the internal buffer the right size ?
         if (ToDef->Sub && ToDef->SubSample!=GeoRef_Options.Sampling) {
            free(ToDef->Sub); ToDef->Sub=NULL;
         }
         ToDef->SubSample=GeoRef_Options.Sampling;
         code=GeoRef_InterpSub(ToRef,ToDef,FromRef,FromDef);
         break;
   }

   return(code);
}
