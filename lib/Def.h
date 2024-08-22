#ifndef _Def_h
#define _Def_h

#include <rmn/List.h>
#include <GeoRef.h>
#include <OGR.h>

#define DEFSELECTTYPE(A,B)  (A->Type>B->Type?A:B)
#define DEFSIGNEDTYPE(A)    ((A->Type==TD_UByte || A->Type==TD_UInt16 || A->Type==TD_UInt32 || A->Type==TD_UInt64)?A->Type+1:A->Type)
#define DEFCLAMP(D,X,Y)      X=(X>D->NI-1?D->NI-1:(X<0?0:X));Y=(Y>D->NJ-1?D->NJ-1:(Y<0?0:Y))
#define DEF2DIN(D,I,J)      ((I)>=D->Limits[0][0] && (I)<=D->Limits[0][1] && (J)>=D->Limits[1][0] && (J)<=D->Limits[1][1])
#define DEF3DIN(D,I,J,K)    ((I)>=D->Limits[0][0] && (I)<=D->Limits[0][1] && (J)>=D->Limits[1][0] && (J)<=D->Limits[1][1] && (K)>=D->Limits[2][0] && (K)<=D->Limits[2][1])
#define DEFVALID(Def,val)   (val==val && val!=Def->NoData)

#define Def_Pointer(DEF,COMP,IDX,PTR) PTR=DEF->Data[COMP]+(IDX)*TDef_Size[DEF->Type];
#define Def_PointerMode(DEF,IDX,PTR) PTR=DEF->Mode+(IDX)*TDef_Size[DEF->Type];

#define Def_Set(DEF,COMP,IDX,VAL) {\
switch(DEF->Type) {\
   case TD_Unknown:break;\
   case TD_Binary: break;\
   case TD_UByte:  ((unsigned char*)DEF->Data[COMP])[IDX]=VAL;break;\
   case TD_Byte:   ((char*)DEF->Data[COMP])[IDX]=VAL; break;\
   case TD_UInt16: ((unsigned short*)DEF->Data[COMP])[IDX]=VAL;break;\
   case TD_Int16:  ((short*)DEF->Data[COMP])[IDX]=VAL;break;\
   case TD_UInt32: ((unsigned int*)DEF->Data[COMP])[IDX]=VAL;break;\
   case TD_Int32:  ((int*)DEF->Data[COMP])[IDX]=VAL;break;\
   case TD_UInt64: ((unsigned long long*)DEF->Data[COMP])[IDX]=VAL;break;\
   case TD_Int64:  ((long long*)DEF->Data[COMP])[IDX]=VAL;break;\
   case TD_Float32:((float*)DEF->Data[COMP])[IDX]=VAL;break;\
   case TD_Float64:((double*)DEF->Data[COMP])[IDX]=VAL;break;\
   }\
}

//unsigned char* Def_GetUByte(char* Data,unsigned int Index) {
//   return ((unsigned char*)Data)[Index];
//}
//void Def_SetUByte(char* Data,unsigned int Index,unsigned char Value) {
//   ((unsigned char*)Data)[Index]=Value;
//}

#define Def_GetQuad(DEF,COMP,IDX,VAL) {\
switch(DEF->Type) {\
   case TD_UByte:  VAL[0]=((unsigned char*)DEF->Data[COMP])[IDX[0]];\
                   VAL[1]=((unsigned char*)DEF->Data[COMP])[IDX[1]];\
                   VAL[2]=((unsigned char*)DEF->Data[COMP])[IDX[2]];\
                   VAL[3]=((unsigned char*)DEF->Data[COMP])[IDX[3]];\
                   break;\
   case TD_Byte:   VAL[0]=((char*)DEF->Data[COMP])[IDX[0]];\
                   VAL[1]=((char*)DEF->Data[COMP])[IDX[1]];\
                   VAL[2]=((char*)DEF->Data[COMP])[IDX[2]];\
                   VAL[3]=((char*)DEF->Data[COMP])[IDX[3]];\
                   break;\
   case TD_UInt16: VAL[0]=((unsigned short*)DEF->Data[COMP])[IDX[0]];\
                   VAL[1]=((unsigned short*)DEF->Data[COMP])[IDX[1]];\
                   VAL[2]=((unsigned short*)DEF->Data[COMP])[IDX[2]];\
                   VAL[3]=((unsigned short*)DEF->Data[COMP])[IDX[3]];\
                   break;\
   case TD_Int16:  VAL[0]=((short*)DEF->Data[COMP])[IDX[0]];\
                   VAL[1]=((short*)DEF->Data[COMP])[IDX[1]];\
                   VAL[2]=((short*)DEF->Data[COMP])[IDX[2]];\
                   VAL[3]=((short*)DEF->Data[COMP])[IDX[3]];\
                   break;\
   case TD_UInt32: VAL[0]=((unsigned int*)DEF->Data[COMP])[IDX[0]];\
                   VAL[1]=((unsigned int*)DEF->Data[COMP])[IDX[1]];\
                   VAL[2]=((unsigned int*)DEF->Data[COMP])[IDX[2]];\
                   VAL[3]=((unsigned int*)DEF->Data[COMP])[IDX[3]];\
                   break;\
   case TD_Int32:  VAL[0]=((int*)DEF->Data[COMP])[IDX[0]];\
                   VAL[1]=((int*)DEF->Data[COMP])[IDX[1]];\
                   VAL[2]=((int*)DEF->Data[COMP])[IDX[2]];\
                   VAL[3]=((int*)DEF->Data[COMP])[IDX[3]];\
                   break;\
   case TD_UInt64: VAL[0]=((unsigned long long*)DEF->Data[COMP])[IDX[0]];\
                   VAL[1]=((unsigned long long*)DEF->Data[COMP])[IDX[1]];\
                   VAL[2]=((unsigned long long*)DEF->Data[COMP])[IDX[2]];\
                   VAL[3]=((unsigned long long*)DEF->Data[COMP])[IDX[3]];\
                   break;\
   case TD_Int64:  VAL[0]=((long long*)DEF->Data[COMP])[IDX[0]];\
                   VAL[1]=((long long*)DEF->Data[COMP])[IDX[1]];\
                   VAL[2]=((long long*)DEF->Data[COMP])[IDX[2]];\
                   VAL[3]=((long long*)DEF->Data[COMP])[IDX[3]];\
                   break;\
   case TD_Float32:VAL[0]=((float*)DEF->Data[COMP])[IDX[0]];\
                   VAL[1]=((float*)DEF->Data[COMP])[IDX[1]];\
                   VAL[2]=((float*)DEF->Data[COMP])[IDX[2]];\
                   VAL[3]=((float*)DEF->Data[COMP])[IDX[3]];\
                   break;\
   case TD_Float64:VAL[0]=((double*)DEF->Data[COMP])[IDX[0]];\
                   VAL[1]=((double*)DEF->Data[COMP])[IDX[1]];\
                   VAL[2]=((double*)DEF->Data[COMP])[IDX[2]];\
                   VAL[3]=((double*)DEF->Data[COMP])[IDX[3]];\
                   break;\
   case TD_Unknown:\
   case TD_Binary: \
   default:        VAL[0]=VAL[1]=VAL[2]=VAL[3]=0.0;break;\
   }\
}

#define Def_GetQuadMod(DEF,IDX,VAL) {\
switch(DEF->Type) {\
   case TD_UByte:  VAL[0]=((unsigned char*)DEF->Mode)[IDX[0]];\
                   VAL[1]=((unsigned char*)DEF->Mode)[IDX[1]];\
                   VAL[2]=((unsigned char*)DEF->Mode)[IDX[2]];\
                   VAL[3]=((unsigned char*)DEF->Mode)[IDX[3]];\
                   break;\
   case TD_Byte:   VAL[0]=((char*)DEF->Mode)[IDX[0]];\
                   VAL[1]=((char*)DEF->Mode)[IDX[1]];\
                   VAL[2]=((char*)DEF->Mode)[IDX[2]];\
                   VAL[3]=((char*)DEF->Mode)[IDX[3]];\
                   break;\
   case TD_UInt16: VAL[0]=((unsigned short*)DEF->Mode)[IDX[0]];\
                   VAL[1]=((unsigned short*)DEF->Mode)[IDX[1]];\
                   VAL[2]=((unsigned short*)DEF->Mode)[IDX[2]];\
                   VAL[3]=((unsigned short*)DEF->Mode)[IDX[3]];\
                   break;\
   case TD_Int16:  VAL[0]=((short*)DEF->Mode)[IDX[0]];\
                   VAL[1]=((short*)DEF->Mode)[IDX[1]];\
                   VAL[2]=((short*)DEF->Mode)[IDX[2]];\
                   VAL[3]=((short*)DEF->Mode)[IDX[3]];\
                   break;\
   case TD_UInt32: VAL[0]=((unsigned int*)DEF->Mode)[IDX[0]];\
                   VAL[1]=((unsigned int*)DEF->Mode)[IDX[1]];\
                   VAL[2]=((unsigned int*)DEF->Mode)[IDX[2]];\
                   VAL[3]=((unsigned int*)DEF->Mode)[IDX[3]];\
                   break;\
   case TD_Int32:  VAL[0]=((int*)DEF->Mode)[IDX[0]];\
                   VAL[1]=((int*)DEF->Mode)[IDX[1]];\
                   VAL[2]=((int*)DEF->Mode)[IDX[2]];\
                   VAL[3]=((int*)DEF->Mode)[IDX[3]];\
                   break;\
   case TD_UInt64: VAL[0]=((unsigned long long*)DEF->Mode)[IDX[0]];\
                   VAL[1]=((unsigned long long*)DEF->Mode)[IDX[1]];\
                   VAL[2]=((unsigned long long*)DEF->Mode)[IDX[2]];\
                   VAL[3]=((unsigned long long*)DEF->Mode)[IDX[3]];\
                   break;\
   case TD_Int64:  VAL[0]=((long long*)DEF->Mode)[IDX[0]];\
                   VAL[1]=((long long*)DEF->Mode)[IDX[1]];\
                   VAL[2]=((long long*)DEF->Mode)[IDX[2]];\
                   VAL[3]=((long long*)DEF->Mode)[IDX[3]];\
                   break;\
   case TD_Float32:VAL[0]=((float*)DEF->Mode)[IDX[0]];\
                   VAL[1]=((float*)DEF->Mode)[IDX[1]];\
                   VAL[2]=((float*)DEF->Mode)[IDX[2]];\
                   VAL[3]=((float*)DEF->Mode)[IDX[3]];\
                   break;\
   case TD_Float64:VAL[0]=((double*)DEF->Mode)[IDX[0]];\
                   VAL[1]=((double*)DEF->Mode)[IDX[1]];\
                   VAL[2]=((double*)DEF->Mode)[IDX[2]];\
                   VAL[3]=((double*)DEF->Mode)[IDX[3]];\
                   break;\
   case TD_Unknown:\
   case TD_Binary:\
   default:        VAL[0]=VAL[1]=VAL[2]=VAL[3]=0.0;\
   }\
}

#define Def_Get(DEF,COMP,IDX,VAL) {\
switch(DEF->Type) {\
   case TD_UByte:  VAL=((unsigned char*)DEF->Data[COMP])[IDX];break;\
   case TD_Byte:   VAL=((char*)DEF->Data[COMP])[IDX]; break;\
   case TD_UInt16: VAL=((unsigned short*)DEF->Data[COMP])[IDX];break;\
   case TD_Int16:  VAL=((short*)DEF->Data[COMP])[IDX];break;\
   case TD_UInt32: VAL=((unsigned int*)DEF->Data[COMP])[IDX];break;\
   case TD_Int32:  VAL=((int*)DEF->Data[COMP])[IDX];break;\
   case TD_UInt64: VAL=((unsigned long long*)DEF->Data[COMP])[IDX];break;\
   case TD_Int64:  VAL=((long long*)DEF->Data[COMP])[IDX];break;\
   case TD_Float32:VAL=((float*)DEF->Data[COMP])[IDX];break;\
   case TD_Float64:VAL=((double*)DEF->Data[COMP])[IDX];break;\
   case TD_Unknown:\
   case TD_Binary:\
   default        :VAL=0.0;\
   }\
}

#define Def_GetMod(DEF,IDX,VAL) {\
switch(DEF->Type) {\
   case TD_UByte:  VAL=((unsigned char*)DEF->Mode)[IDX];break;\
   case TD_Byte:   VAL=((char*)DEF->Mode)[IDX]; break;\
   case TD_UInt16: VAL=((unsigned short*)DEF->Mode)[IDX];break;\
   case TD_Int16:  VAL=((short*)DEF->Mode)[IDX];break;\
   case TD_UInt32: VAL=((unsigned int*)DEF->Mode)[IDX];break;\
   case TD_Int32:  VAL=((int*)DEF->Mode)[IDX];break;\
   case TD_UInt64: VAL=((unsigned long long*)DEF->Mode)[IDX];break;\
   case TD_Int64:  VAL=((long long*)DEF->Mode)[IDX];break;\
   case TD_Float32:VAL=((float*)DEF->Mode)[IDX];break;\
   case TD_Float64:VAL=((double*)DEF->Mode)[IDX];break;\
   case TD_Unknown:\
   case TD_Binary:\
   default        :VAL=0.0;\
   }\
}

#define Def_SetMod(DEF,IDX,VAL) {\
switch(DEF->Type) {\
   case TD_Unknown:break;\
   case TD_Binary: break;\
   case TD_UByte:  ((unsigned char*)DEF->Mode)[IDX]=VAL;break;\
   case TD_Byte:   ((char*)DEF->Mode)[IDX]=VAL; break;\
   case TD_UInt16: ((unsigned short*)DEF->Mode)[IDX]=VAL;break;\
   case TD_Int16:  ((short*)DEF->Mode)[IDX]=VAL;break;\
   case TD_UInt32: ((unsigned int*)DEF->Mode)[IDX]=VAL;break;\
   case TD_Int32:  ((int*)DEF->Mode)[IDX]=VAL;break;\
   case TD_UInt64: ((unsigned long long*)DEF->Mode)[IDX]=VAL;break;\
   case TD_Int64:  ((long long*)DEF->Mode)[IDX]=VAL;break;\
   case TD_Float32:((float*)DEF->Mode)[IDX]=VAL;break;\
   case TD_Float64:((double*)DEF->Mode)[IDX]=VAL;break;\
   }\
}

// Value data type
typedef enum {
    TD_Unknown = 0,
    TD_Binary  = 1,
    TD_UByte   = 2,
    TD_Byte    = 3,
    TD_UInt16  = 4,
    TD_Int16   = 5,
    TD_UInt32  = 6,
    TD_Int32   = 7,
    TD_UInt64  = 8,
    TD_Int64   = 9,
    TD_Float32 = 10,
    TD_Float64 = 11,
} TDef_Type;

extern int TDef_Size[];
extern int TDef_DTYP[];

typedef struct TDef {
   double *Buffer,*Aux;          // Buffer temporaire
   int    *Accum;                // Accumulation Buffer temporaire
   char   *Mask;                 // Masque a appliquer au traitement sur le champs
   char   *Data[4];              // Composantes du champs (Pointeurs sur les donnees)
   char   *Mode;                 // Module des champs Data is vectoriel
   char   *Dir;                  // Direction si vectoriel
   float  *Pres,*PresLS,*Height; // Pression au sol
   float  *Sub;                  // Sub grid resolutions values
   OGRGeometryH *Pick,*Poly;     // Geometry used in various interpolation method
   TList  *Segments;             // Liste d'objets de rendue

   double  NoData;               // Valeur de novalue
   TDef_Type Type;               // Type de donnees du champs
   int NI,NJ,NK,NC,NIJ;          // Dimensions du champs
   int Idx;                      // Index displacement into supergrid

   int     CellDim;              // Defined grid point coverage, point=1 or area=2
   double  CoordLimits[2][2];    // Limits of processing in latlon
   int     Limits[3][2];         // Limits of processing in grid points
   int     Level;                // Niveau courant
   int     Sample,SubSample;     // Sample interval in grid points
   char    Alias;                // Alias d'un autre TDef (Pointe sur d'autres donnees)
} TDef;

void  Def_Clear(TDef *Def);
int   Def_Compat(TDef *ToDef,TDef *FromDef);
TDef *Def_Copy(TDef *Def);
TDef *Def_CopyPromote(TDef *Def,TDef_Type Type);
void  Def_Free(TDef *Def);
TDef *Def_New(int NI,int NJ,int NK,int Dim,TDef_Type Type,int Alias);
TDef *Def_Create(int NI,int NJ,int NK,TDef_Type Type,char* Comp0,char* Comp1,char* Mask);
TDef *Def_Resize(TDef *Def,int NI,int NJ,int NK);
int   Def_Paste(TDef *ToDef,TDef *DefPaste,int X0,int Y0);

int   Def_GetValue(TGeoRef *Ref,TDef *Def,TGeoOptions *Opt,int C,double X,double Y,double Z,double *Length,double *ThetaXY);

int   GeoRef_Cell2OGR(OGRGeometryH Geom,TGeoRef *ToRef,TGeoRef *FromRef,int I,int J,int Seg);
int   GeoRef_Rasterize(TGeoRef *ToRef,TDef *ToDef,TGeoOptions *Opt,OGRGeometryH Geom,double Value);
int   GeoRef_InterpDef(TGeoRef *ToRef,TDef *ToDef,TGeoRef *FromRef,TDef *FromDef,TGeoOptions *Opt,int Final);
int   GeoRef_InterpAverage(TGeoRef *ToRef,TDef *ToDef,TGeoRef *FromRef,TDef *FromDef,TGeoOptions *Opt,double *Table,TDef **lutDef, int lutSize, TDef *TmpDef,int Final);
int   GeoRef_InterpConservative(TGeoRef *ToRef,TDef *ToDef,TGeoRef *FromRef,TDef *FromDef,TGeoOptions *Opt,int Final);
int   GeoRef_InterpSub(TGeoRef *ToRef,TDef *ToDef,TGeoRef *FromRef,TDef *FromDef,TGeoOptions *Opt);
int   GeoRef_InterpOGR(TGeoRef *ToRef,TDef *ToDef,TGeoRef *LayerRef,OGR_Layer *Layer,TGeoOptions *Opt,char *Field,double Value,int Final);

void GeoScan_Init(TGeoScan *Scan);
void GeoScan_Clear(TGeoScan *Scan);
int  GeoScan_Get(TGeoScan *Scan,TGeoRef *ToRef,struct TDef *ToDef,TGeoRef *FromRef,struct TDef *FromDef,TGeoOptions *Opt,int X0,int Y0,int X1,int Y1,int Dim);

#endif
