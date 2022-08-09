/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Fonctions et definitions relatives aux fichiers standards et rmnlib
 * Fichier      : RPN.h
 * Creation     : Avril 2006 - J.P. Gauthier
 *
 * Description:
 *
 * License:
 *    This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation,
 *    version 2.1 of the License.
 *
 *    This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the
 *    Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 *    Boston, MA 02111-1307, USA.
 *
 *==============================================================================
 */
#ifndef _RPN_h
#define _RPN_h

#define RPNMAX 2048                            ///< Maximum size of fstinl list
#define C_TO_FTN(I,J,NI) (int)((NI)*(J)+(I))   ///< Array index conversion

#define X_PI    0
#define X_PJ    1
#define X_D60   2
#define X_DGRW  3
#define X_IG1   0
#define X_IG2   1
#define X_IG3   2
#define X_IG4   3
#define X_LAT1  0
#define X_LON1  1
#define X_LAT2  2
#define X_LON2  3
#define X_SWLAT 0
#define X_SWLON 1
#define X_DLAT  2
#define X_DLON  3
#define X_TD60  0
#define X_TDGRW 1
#define X_CLAT  2
#define X_CLON  3

struct TGeoRef;

typedef struct TRPNFile {
   char *CId;              ///< Identificateur du fichier
   char *Name;             ///< Path complet du fichier
   int  Open;              ///< Etat du fichier
   unsigned int Id;        ///< Numero d'unite du fichier
   int  NRef;              ///< Nombre de reference
   char Mode;              ///< Mode d'ouverture du fichier (r,w,a)
} TRPNFile;

typedef struct TRPNHeader {
   TRPNFile *File;         ///< Fichier dont provient le champs
   int   FID;               ///< FID dont provient le champs
   int   KEY;               ///< Cle du champs
   int   DATEO;             ///< Date d'origine du champs
   int   DATEV;             ///< Date de validitee du champs
   int   DEET;              ///< Duree d'un pas de temps
   int   NPAS;              ///< Pas de temps
   int   NBITS;             ///< Nombre de bits du champs
   int   DATYP;             ///< Type de donnees
   int   IP1,IP2,IP3;       ///< Specificateur du champs
   int   NI,NJ,NK,NIJ;      ///< Dimensions
   int   IG[4],IGREF[4];    ///< Descripteur de grille
   float XG[4],XGREF[4];    ///< Descripteur de grille
   int   SWA;
   int   LNG;
   int   DLTF;
   int   UBC;
   int   EX1,EX2,EX3;
   char  TYPVAR[3];         ///< Type de variable
   char  NOMVAR[5];         ///< Nom de la variable
   char  ETIKET[13];        ///< Etiquette du champs
   char  GRTYP[2],GRREF[2]; ///< Type de grilles er references
} TRPNHeader;

typedef struct TRPNField {
   TRPNHeader     Head;    ///< Entete du champs
   struct TDef    *Def;    ///< Definition des donnees
   struct TGeoRef *GRef;   ///< Reference geographique horizontale
   struct TZRef   *ZRef;   ///< Reference geographique verticale
} TRPNField;

int  RPN_CopyDesc(int FIdTo,TRPNHeader* const H);
int  RPN_IsDesc(const char* restrict Var);
void RPN_FileLock(void);
void RPN_FileUnlock(void);
void RPN_FieldLock(void);
void RPN_FieldUnlock(void);

TRPNField* RPN_FieldNew();
void       RPN_FieldFree(TRPNField *Fld);
TRPNField* RPN_FieldReadIndex(int FileId,int Index,TRPNField *Fld);
TRPNField* RPN_FieldRead(int FileId,int DateV,char *Eticket,int IP1,int IP2,int IP3,char *TypVar,char *NomVar);
int        RPN_FieldReadComponent(TRPNHeader *Head,float **Ptr,char *Var,int Grid,int Force);
int        RPN_FieldReadLevels(TRPNField *Field);
int        RPN_FieldWrite(int FileId,TRPNField *Field);
struct TGeoRef* RPN_FieldReadGrid(TRPNField *Field);
int        RPN_ReadGrid(struct TGeoRef *Ref);
void       RPN_CopyHead(TRPNHeader *To,TRPNHeader *From);
int        RPN_FieldTile(int FID,struct TDef *Def,TRPNHeader *Head,struct TGeoRef *Ref,struct TZRef *ZRef,int Comp,int NI,int NJ,int Halo,int DATYP,int NPack,int Rewrite,int Compress);

int RPN_GetAllFields(int FID,int DateV,char *Etiket,int Ip1,int Ip2,int Ip3,char *Typvar,char *Nomvar,int **Arr,int *Size);
int RPN_GetAllDates(int *Flds,int NbFlds,int Uniq,int **DateV,int *NbDateV);
int RPN_GetAllIps(int *Flds,int NbFlds,int IpN,int Uniq,int **Ips,int *NbIp);

int RPN_GenerateIG(int *IG1,int *IG2,int *IG3);
int RPN_LinkFiles(char **Files,int N);
int RPN_UnLinkFiles(int FID);
int RPN_LinkPattern(const char* Pattern);

#ifdef HAVE_RMN
#include "rpnmacros.h"

// RPN threadsafe (lock) fstd functions
void cs_fstunlockid(int Unit);
int  cs_fstlockid();
int  cs_fstfrm(int Unit);
int  cs_fstouv(char *Path,char *Mode);
int  cs_fstflush(int Unit);
int  cs_fstinl(int Unit,int *NI,int *NJ,int *NK,int DateO,char *Etiket,int IP1,int IP2,int IP3,char* TypVar,char *NomVar,int *List,int *Nb,int Max);
int  cs_fstinf(int Unit,int *NI,int *NJ,int *NK,int DateO,char *Etiket,int IP1,int IP2,int IP3,char* TypVar,char *NomVar);
int  cs_fstprm(int Unit,int *DateO,int *Deet,int *NPas,int *NI,int *NJ,int *NK,int *NBits,int *Datyp,int *IP1,int *IP2,int *IP3,char* TypVar,char *NomVar,char *Etiket,char *GrTyp,int *IG1,int *IG2,int *IG3,int *IG4,int *Swa,int *Lng,int *DLTF,int *UBC,int *EX1,int *EX2,int *EX3);
int  cs_fstlir(void *Buf,int Unit,int *NI,int *NJ,int *NK,int DateO,char *Etiket,int IP1,int IP2,int IP3,char* TypVar,char *NomVar);
int  cs_fstluk(void *Data,int Idx,int *NI,int *NJ,int *NK);
int  cs_fstsui(int Unit,int *NI,int *NJ,int *NK);
int  cs_fstlukt(void *Data,int Unit,int Idx,char *GRTYP,int *NI,int *NJ,int *NK);
int  cs_fstecr(void *Data,int NPak,int Unit, int DateO,int Deet,int NPas,int NI,int NJ,int NK,int IP1,int IP2,int IP3,char* TypVar,char *NomVar,char *Etiket,char *GrTyp,int IG1,int IG2,int IG3,int IG4,int DaTyp,int Over); 

#endif
#endif
