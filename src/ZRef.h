#ifndef _ZRef_h
#define _ZRef_h

#include <rmn.h>

// Level related constants and functions
#define LVL_NIL         -1  //  No conversion
#define LVL_MASL         0  //  Meters above sea level
#define LVL_SIGMA        1  //  P/Ps
#define LVL_PRES         2  //  Pressure mb
#define LVL_UNDEF        3  //  units are user defined
#define LVL_MAGL         4  //  Meters above ground level
#define LVL_HYBRID       5  //  Hybrid levels
#define LVL_THETA        6  //  ?
#define LVL_MBSL         7  //  Meters below sea level
#define LVL_GALCHEN      8  //  Original Gal-Chen -not in convip (JP Defined)
#define LVL_NBR          9  //  Nombre d'elements
#define LVL_HOUR        10  //  Hours
#define LVL_ANGLE       11  //  Radar angles (JP defined)
#define LVL_INT         15  //  Entiers (reserve)
#define LVL_IDX         17  //  Index de matrice
#define LVL_MPRES       21  //  Metres-pression
#define LVL_ETA         32  //  (Pt-P)/(Pt-Ps) -not in convip

#define MB2PA           100.0f           // Constant for converting pressure from [mb] to [Pa]
#define PA2MB           0.01f            // Constant for converting pressure from [Pa] to [mb]

#define PRESS2METER(LVL) (LVL>0?-8409.1*log(LVL)/1200.0:0)
#define SIGMA2METER(LVL) (LVL>0?-8409.1*log(LVL):0)
#define ETA2METER(LVL)   (-8409.1*log(LVL+1e-32))

typedef enum { DEFAULT=1,NEW=2,OLD=3 } TZRef_IP1Mode;

// Vertical referential definition
typedef struct TZRef {
   char  *Name;           // Reference name
   char  *VGD;            // VGrid descriptor
   float *Levels;         // Levels list
   float *A,*B;           // Pressure calculation factors
   float *P0,*P0LS;       // Pressure at surface
   float *PCube;          // 3D Pressure cube
   int32_t    LevelNb;    // Number of Levels
   int32_t    NRef;       // Reference count
   int32_t    Version;    // Version
   int32_t    Type;       // Type of levels
   int32_t    SLEVE;      // Is SLEVE type (boolean 0/1)
   float  POff;           // Pressure offset from level
   float  PTop;           // Pressure at top of atmosphere
   float  PRef;           // Reference pressure
   float  RCoef[2];       // Hybrid level coefficient
   float  ETop;           // Eta coordinate a top

   TZRef_IP1Mode Style;   // IP style
} TZRef;

#define ZRef_Incr(ZREF) __sync_add_and_fetch(&ZREF->NRef,1)
#define ZRef_Decr(ZREF) __sync_sub_and_fetch(&ZREF->NRef,1)

TZRef*  ZRef_New(void);
TZRef*  ZRef_Define(int32_t Type,int32_t NbLevels,float *Levels);
int32_t ZRef_Free(TZRef *ZRef);
int32_t ZRef_Equal(TZRef *Zref0,TZRef *ZRef1);
TZRef*  ZRef_Copy(TZRef *ZRef);
TZRef*  ZRef_HardCopy(TZRef *ZRef);
int32_t ZRef_DecodeRPN(TZRef *ZRef,fst_file* File);
int32_t ZRef_SetRestrictLevels(float *Levels,int32_t NbLevels);
int32_t ZRef_AddRestrictLevel(float Level);
int32_t ZRef_GetLevels(TZRef *ZRef,const fst_record* restrict const H,int32_t Invert);
double  ZRef_K2Pressure(TZRef* restrict const ZRef,double P0,double P0LS,int32_t K);
int32_t ZRef_KCube2Pressure(TZRef* restrict const ZRef,float *P0,float *P0LS,int32_t NIJ,int32_t Log,float *Pres);
int32_t ZRef_KCube2Meter(TZRef* restrict const ZRef,float *GZ,const int32_t NIJ,float *Height);
double  ZRef_Level2Pressure(TZRef* restrict const ZRef,double P0,double P0LS,double Level);
double  ZRef_Pressure2Level(TZRef* restrict const ZRef,double P0,double Pressure);
double  ZRef_IP2Meter(int32_t IP);
double  ZRef_Level2Meter(double Level,int32_t Type);
double  ZRef_IP2Level(int32_t IP,int32_t *Type);
int32_t ZRef_Level2IP(float Level,int32_t Type,TZRef_IP1Mode Mode);
int32_t ZRef_IPFormat(char *Buf,int32_t IP,int32_t Interval);

const char** ZRef_LevelNames();
const char*  ZRef_LevelName(int32_t Type);
const char** ZRef_LevelUnits();
const char*  ZRef_LevelUnit(int32_t Type);

#endif
