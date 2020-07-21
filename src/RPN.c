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
#ifdef HAVE_RMN

#include <pthread.h>

#include "App.h"
#include "GeoRef.h"
#include "Def.h"
#include "EZGrid.h"
#include "fnom.h"
#include <glob.h>

static int **LNK_FID = NULL;
static int LNK_NB = 0;

static const char *RPN_Desc[]={ ">>  ","^^  ","^>  ","!!  ","##  ","HY  ","PROJ","MTRX",NULL };

static char FGFDTLock[1000];

static pthread_mutex_t RPNFieldMutex=PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t RPNFileMutex=PTHREAD_MUTEX_INITIALIZER;

void RPN_FileLock() {
   pthread_mutex_lock(&RPNFileMutex);
}
void RPN_FileUnlock() {
   pthread_mutex_unlock(&RPNFileMutex);
}

void RPN_FieldLock() {
   pthread_mutex_lock(&RPNFieldMutex);
}
void RPN_FieldUnlock() {
   pthread_mutex_unlock(&RPNFieldMutex);
}

void cs_fstunlockid(int Unit) {
   pthread_mutex_lock(&RPNFileMutex);
   FGFDTLock[Unit-1]=0;
   pthread_mutex_unlock(&RPNFileMutex);
}

int cs_fstlockid() {

   int id;

   pthread_mutex_lock(&RPNFileMutex);
   
   // Id 5 and 6 are stdout and stderr   
   FGFDTLock[4]=FGFDTLock[5]=1;  

   for (id=0;id<1000;id++) {
      if (FGFDTLock[id]==0) {
         FGFDTLock[id]=1;
         id++;
         break;
      }
   }
   pthread_mutex_unlock(&RPNFileMutex);
   return(id);
}

int cs_fstouv(char *Path,char *Mode) {

   int err=-1,id=-1;
   char mode[32];

   if (Path) {
      id=cs_fstlockid();
      pthread_mutex_lock(&RPNFileMutex);
      if (index(Path,':') && Path[0]!=':') {
         strcpy(mode,Mode);
         strcat(mode,"+REMOTE");
         err=c_fnom(&id,Path,mode,0);
      } else {
         err=c_fnom(&id,Path,Mode,0);
      }
      if (err>=0) {
         err=c_fstouv(id,"RND");
      }
      pthread_mutex_unlock(&RPNFileMutex);
   }
   return(err<0?err:id);
}

int cs_fstflush(int Unit) {
   int err;

   pthread_mutex_lock(&RPNFileMutex);
   err=c_fstfrm(Unit);
   err=c_fstouv(Unit,"RND");
   pthread_mutex_unlock(&RPNFileMutex);

   return(err);
}

int cs_fstfrm(int Unit) {
   int err;

   pthread_mutex_lock(&RPNFileMutex);
   err=c_fstfrm(Unit);
   err=c_fclos(Unit);
   FGFDTLock[Unit-1]=0;
   pthread_mutex_unlock(&RPNFileMutex);

   return(err);
}

int cs_fstlir(void *Buf,int Unit,int *NI,int *NJ,int *NK,int DateV,char *Etiket,int IP1,int IP2,int IP3,char* TypVar,char *NomVar) {

   int err;

   pthread_mutex_lock(&RPNFieldMutex);
   err=c_fstlir(Buf,Unit,NI,NJ,NK,DateV,Etiket,IP1,IP2,IP3,TypVar,NomVar);
   pthread_mutex_unlock(&RPNFieldMutex);

   return(err);
}

int cs_fstinf(int Unit,int *NI,int *NJ,int *NK,int DateO,char *Etiket,int IP1,int IP2,int IP3,char* TypVar,char *NomVar) {

   int err;

   pthread_mutex_lock(&RPNFieldMutex);
   err=c_fstinf(Unit,NI,NJ,NK,DateO,Etiket,IP1,IP2,IP3,TypVar,NomVar);
   pthread_mutex_unlock(&RPNFieldMutex);

   return(err);
}

int cs_fstprm(int Idx,int *DateO,int *Deet,int *NPas,int *NI,int *NJ,int *NK,int *NBits,int *Datyp,int *IP1,int *IP2,int *IP3,char* TypVar,char *NomVar,char *Etiket,char *GrTyp,int *IG1_JP,int *IG2_JP,int *IG3_JP,int *IG4_JP,int *Swa,int *Lng,int *DLTF,int *UBC,int *EX1,int *EX2,int *EX3) {

   int err;

   pthread_mutex_lock(&RPNFieldMutex);
   err=c_fstprm(Idx,DateO,Deet,NPas,NI,NJ,NK,NBits,Datyp,IP1,IP2,IP3,TypVar,NomVar,Etiket,GrTyp,IG1_JP,IG2_JP,IG3_JP,IG4_JP,Swa,Lng,DLTF,UBC,EX1,EX2,EX3);
   pthread_mutex_unlock(&RPNFieldMutex);

   return(err);
}

int cs_fstinl(int Unit,int *NI,int *NJ,int *NK,int DateO,char *Etiket,int IP1,int IP2,int IP3,char* TypVar,char *NomVar,int *List,int *Nb,int Max) {

   int err;

   pthread_mutex_lock(&RPNFieldMutex);
   err=c_fstinl(Unit,NI,NJ,NK,DateO,Etiket,IP1,IP2,IP3,TypVar,NomVar,List,Nb,Max);
   pthread_mutex_unlock(&RPNFieldMutex);

   return(err);
}

int cs_fstluk(void *Data,int Idx,int *NI,int *NJ,int *NK) {

   int err;

   pthread_mutex_lock(&RPNFieldMutex);
   err=c_fstluk(Data,Idx,NI,NJ,NK);
   pthread_mutex_unlock(&RPNFieldMutex);

   return(err);
}

int cs_fstsui(int Unit,int *NI,int *NJ,int *NK) {

   int key;

   pthread_mutex_lock(&RPNFieldMutex);
   key=c_fstsui(Unit,NI,NJ,NK);
   pthread_mutex_unlock(&RPNFieldMutex);

   return(key);
}

int cs_fstlukt(void *Data,int Unit,int Idx,char *GRTYP,int *NI,int *NJ,int *NK) {

   int    err=-1;
   TGrid *grid;

   if (GRTYP[0]=='#') {
      if ((grid=EZGrid_ReadIdx(Unit,Idx,0))) {
         if (EZGrid_TileBurnAll(grid,0,Data)) {
            err=0;
         }
         EZGrid_Free(grid);
      }
      *NI=grid->H.NI;
      *NJ=grid->H.NJ;
      *NK=grid->H.NK;
   } else {
      pthread_mutex_lock(&RPNFieldMutex);
      err=c_fstluk(Data,Idx,NI,NJ,NK);
      pthread_mutex_unlock(&RPNFieldMutex);
   }

   return(err);
}

int cs_fstecr(void *Data,int NPak,int Unit, int DateO,int Deet,int NPas,int NI,int NJ,int NK,int IP1,int IP2,int IP3,char* TypVar,char *NomVar,char *Etiket,char *GrTyp,int IG1_JP,int IG2_JP,int IG3_JP,int IG4_JP,int DaTyp,int Over) {
   int err;

   pthread_mutex_lock(&RPNFieldMutex);
   err=c_fstecr(Data,NULL,NPak,Unit,DateO,Deet,NPas,NI,NJ,NK,IP1,IP2,IP3,TypVar,NomVar,Etiket,GrTyp,IG1_JP,IG2_JP,IG3_JP,IG4_JP,DaTyp,Over);
   pthread_mutex_unlock(&RPNFieldMutex);

   return(err);
}

/*----------------------------------------------------------------------------
 * Nom      : <RPN_IntIdNew>
 * Creation : Janvier 2012 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Wrapper autour de c_ezqkdef pour garder un compte du nombre
 *           d'allocation de chaque grille ezscint
 *
 * Parametres :
 *
 * Retour:
 *  <Id>      : Identificateur de grille ezscint
 *
 * Remarques :
 *----------------------------------------------------------------------------
 */
//TODO: get rid off
int RPN_IntIdNew(int NI,int NJ,char* GRTYP,int IG1_JP,int IG2_JP,int IG3_JP, int IG4_JP,int FID) {

   TGeoRef *GRef;

   if (GRTYP[0]!='M' && GRTYP[0]!='W' && GRTYP[0]!='V') {
//      RPN_IntLock();
      GRef=GeoRef_RPNCreate(NI,NJ,GRTYP,IG1_JP,IG2_JP,IG3_JP,IG4_JP,FID);
//      RPN_IntUnlock();
   }
   
   return GRef;
}

/*----------------------------------------------------------------------------
 * Nom      : <RPN_IntIdFree>
 * Creation : Janvier 2012 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Wrapper autour de c_gdrls pour garder un compte du nombre
 *           d'allocation de chaque grille ezscint
 *
 * Parametres :
 *  <Id>      : Identificateur de grille ezscint
 *
 * Retour:
 *  <n>      : Compte du nobre d'allocation de la grille ezscint
 *
 * Remarques :
 *----------------------------------------------------------------------------
 */
//TODO: get rid off
int RPN_IntIdFree(int Id) {

   int n=-1;

   if (Id<0)
      return(n);
//TODO:Check on grid cache
//  RPN_IntLock();
//         c_gdrls(Id);
//   RPN_IntUnlock();

   return(n);
}


TRPNField* RPN_FieldNew(int NI,int NJ,int NK,int NC,TDef_Type Type) {

   TRPNField  *fld;

   fld=(TRPNField*)calloc(1,sizeof(TRPNField));
   if (!(fld->Def=Def_New(NI,NJ,NK,NC,Type))) {
      App_Log(ERROR,"%s: Could not allocate memory\n",__func__);
      return(NULL);
   }

   fld->GRef=NULL;
   fld->Head.NI=NI;
   fld->Head.NJ=NJ;
   fld->Head.NK=NK;

   return(fld);
}

void RPN_FieldFree(TRPNField *Fld) {

   if (Fld->GRef) GeoRef_Free(Fld->GRef);
   if (Fld->Def)  Def_Free(Fld->Def);

   free(Fld);
}

TRPNField* RPN_FieldReadIndex(int FileId,int Index,TRPNField *Fld) {

   TRPNField  *fld;
   TRPNHeader  h;
   int         ok,type;
   double      nhour;
   float       lvl;

   h.FID=FileId;
   h.KEY=Index;
   h.GRTYP[0]=h.GRTYP[1]='\0';

   strcpy(h.NOMVAR,"    ");
   strcpy(h.TYPVAR,"  ");
   strcpy(h.ETIKET,"            ");
   ok=cs_fstprm(h.KEY,&h.DATEO,&h.DEET,&h.NPAS,&h.NI,&h.NJ,&h.NK,&h.NBITS,
         &h.DATYP,&h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,
         h.GRTYP,&h.IG1_JP,&h.IG2_JP,&h.IG3_JP,&h.IG4_JP,&h.SWA,&h.LNG,&h.DLTF,
         &h.UBC,&h.EX1,&h.EX2,&h.EX3);

   if (ok<0) {
      App_Log(ERROR,"%s: Could not get field information (c_fstprm failed)\n",__func__);
      return(NULL);
   }

   // Calculer la date de validitee du champs
   if (h.DATEO!=0) {
      nhour=((double)h.NPAS*h.DEET)/3600.0;
      f77name(incdatr)(&h.DATEV,&h.DATEO,&nhour);
      if (h.DATEV==101010101) h.DATEV=0;
   } else {
      h.DATEV=0;
   }
   
   // Supprimer les espaces inutiles
   strtrim(h.NOMVAR,' ');
   strtrim(h.TYPVAR,' ');
   strtrim(h.ETIKET,' ');

   // If a TRPNField is passed, fill it instead of new
   if (Fld) {
      fld=Fld;
   } else {
      fld=(TRPNField*)calloc(1,sizeof(TRPNField));
      if (!(fld->Def=Def_New(h.NI,h.NJ,h.NK,1,TD_Float32))) {
         App_Log(ERROR,"%s: Could not allocate memory for fld\n",__func__);
         return(NULL);
      }
   }

   // Recuperer les donnees du champs
   c_fst_data_length(TDef_Size[fld->Def->Type]);
   if ((ok=cs_fstlukt(fld->Def->Data[0],h.FID,h.KEY,h.GRTYP,&h.NI,&h.NJ,&h.NK))<0) {
      App_Log(ERROR,"%s: Could not read field data (c_fstluk failed)\n",__func__);
      return(NULL);
   }

   // If a TRPNField is passed, we assume this has already been done
   if (!Fld) {
      // Recuperer les type de niveaux et forcer ETA pour SIGMA
      lvl=ZRef_IP2Level(h.IP1,&type);
      type=type==LVL_SIGMA?LVL_ETA:type;

      if (h.GRTYP[0]!='W') {
         fld->GRef=GeoRef_RPNCreate(h.NI,h.NJ,h.GRTYP,h.IG1_JP,h.IG2_JP,h.IG3_JP,h.IG4_JP,h.FID);
      }
      fld->ZRef=ZRef_Define(type,h.NK,&lvl);
   //   if (grtyp[0]=='U') {
   //      FSTD_FieldSubBuild(field);
   //   }

      GeoRef_Qualify(fld->GRef);
   }
   memcpy(&fld->Head,&h,sizeof(TRPNHeader));

   return(fld);
}

TRPNField* RPN_FieldRead(int FileId,int DateV,char *Eticket,int IP1,int IP2,int IP3,char *TypVar,char *NomVar) {

   TRPNField  *fld;
   TRPNHeader  h;

   // Rechercher et lire l'information de l'enregistrement specifie
   h.KEY=cs_fstinf(FileId,&h.NI,&h.NJ,&h.NK,DateV,Eticket,IP1,IP2,IP3,TypVar,NomVar);

   if (h.KEY<0) {
      App_Log(ERROR,"%s: Specified field does not exist (c_fstinf failed)\n",__func__);
      return(NULL);
   }

   fld=RPN_FieldReadIndex(FileId,h.KEY,NULL);

   // Recheck for fstsui to be ok
   cs_fstinf(FileId,&h.NI,&h.NJ,&h.NK,DateV,Eticket,IP1,IP2,IP3,TypVar,NomVar);

   return(fld);
}

/*----------------------------------------------------------------------------
 * Nom      : <RPN_FieldReadComponent>
 * Creation : Mars 2006 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Effectue la lecture d'un champ complementaire.
 *
         d=Proj->L*0.5;
         d=Proj->L*0.5;
 * Parametres :
 *  <Head>    : Entete de la donnee
 *  <Ptr>     : Pointeur sur le vecteur a allouer
 *  <Var>     : Variable a lire
 *  <Grid>    : Utiliser les standard grille (1=IP<->IG,0=IP<->IP,-1=IP<->-1);
 *
 * Retour:
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
int RPN_FieldReadComponent(TRPNHeader *Head,float **Ptr,char *Var,int Grid,int Force) {

   int key=0,ni=0,nj=0,nk=0;

   if (Var && (!*Ptr || Force)) {
      if (Grid==1) {
         // Look for corresponding time and if not use any time
         if ((key==cs_fstinf(Head->FID,&ni,&nj,&nk,Head->DATEV,"",Head->IG1_JP,Head->IG2_JP,Head->IG3_JP,"",Var))<=0) {
            key=cs_fstinf(Head->FID,&ni,&nj,&nk,-1,"",Head->IG1_JP,Head->IG2_JP,Head->IG3_JP,"",Var);
         }
      } else if (Grid==0) {
         key=cs_fstinf(Head->FID,&ni,&nj,&nk,Head->DATEV,Head->ETIKET,Head->IP1,Head->IP2,Head->IP3,Head->TYPVAR,Var);
      } else {
         key=cs_fstinf(Head->FID,&ni,&nj,&nk,Head->DATEV,Head->ETIKET,-1,Head->IP2,Head->IP3,Head->TYPVAR,Var);
      }

      if (key<0) {
         // Too many warnings so we catch it later
         // App_Log(WARNING,"%s: Could not find component field %s (c_fstinf failed)\n",__func__,Var);
         return(0);
      } else {
         if (!*Ptr) {
            //TODO: 64 bit descriptorspi++
            if (!(*Ptr=(float*)malloc(ni*nj*nk*sizeof(float)))) {
               App_Log(ERROR,"%s: Not enough memory to read coordinates fields\n",__func__);
               return(0);
            }
         }
         cs_fstluk(*Ptr,key,&ni,&nj,&nk);
      }
   }
   return(ni*nj*nk);
}

/*----------------------------------------------------------------------------
 * Nom      : <FSTD_FieldReadGrid>
 * Creation : Mars 2006 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Effectue les calculs et l'affichage pour le type particule.
 *
 * Parametres :
 *  <Field>   : Adresse des valeurs du champs
 *
 * Retour:
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
int RPN_FieldReadGrid(TRPNField *Field) {

   TRPNHeader *head=&Field->Head;
   int         key,ni,nj,nk;

   // !! -> NI==3 ### Field->Def->NI<4 || 
   if (!Field->GRef || !(Field->GRef->Type&(GRID_SPARSE|GRID_VARIABLE|GRID_VERTICAL)) || (Field->GRef->NY==1 && Field->GRef->Grid[0]!='Y' && Field->GRef->Grid[1]!='Y' && Field->GRef->Grid[0]!='M'))
      return(0);

   if ((!Field->GRef->AY || !Field->GRef->AX) && head->FID>=0) {

      switch(Field->GRef->Grid[0]) {
         case 'M':
            if (!Field->GRef->AY) RPN_FieldReadComponent(head,&Field->GRef->AY,"^^",1,0);
            if (!Field->GRef->AX) RPN_FieldReadComponent(head,&Field->GRef->AX,">>",1,0);

            /* Lire le champs d'indexes*/
            if (!Field->GRef->Idx) {
               key=cs_fstinf(head->FID,&ni,&nj,&nk,-1,"",head->IG1_JP,head->IG2_JP,head->IG3_JP,"","##");
               if (key < 0) {
                  App_Log(ERROR,"%s: Could not find index field %s (c_fstinf failed)",__func__,"##");
                  return(0);
               } else {
                  Field->GRef->NIdx=ni*nj*nk;
                  if (!(Field->GRef->Idx=(unsigned int*)malloc(Field->GRef->NIdx*sizeof(unsigned int)))) {
                     App_Log(ERROR,"%s: Not enough memory to read coordinates fields",__func__);
                     return(0);
                  }
                  cs_fstluk((float*)Field->GRef->Idx,key,&ni,&nj,&nk);
               }
               GeoRef_BuildIndex(Field->GRef);
            }
            break;

         case 'W':
            if (Field->GRef->Grid[1]=='X' || Field->GRef->Grid[1]=='Y' || Field->GRef->Grid[1]=='Z') {
               if (!Field->GRef->AY) RPN_FieldReadComponent(head,&Field->GRef->AY,"^^",1,0);
               if (!Field->GRef->AX) RPN_FieldReadComponent(head,&Field->GRef->AX,">>",1,0);
            }
            
            if (Field->GRef->Grid[1]=='Y') {
               if (!Field->GRef->AY) RPN_FieldReadComponent(head,&Field->GRef->AY,"LA",0,0);
               if (!Field->GRef->AX) RPN_FieldReadComponent(head,&Field->GRef->AX,"LO",0,0);
               if (!Field->GRef->Hgt) RPN_FieldReadComponent(head,&Field->GRef->Hgt,"ZH",0,0);
            }
            break;

         case 'Y':
            if (!Field->GRef->AY) RPN_FieldReadComponent(head,&Field->GRef->AY,"LA",0,0);
            if (!Field->GRef->AY) RPN_FieldReadComponent(head,&Field->GRef->AY,"^^",1,0);
            if (!Field->GRef->AX) RPN_FieldReadComponent(head,&Field->GRef->AX,"LO",0,0);
            if (!Field->GRef->AX) RPN_FieldReadComponent(head,&Field->GRef->AX,">>",1,0);
            if (!Field->GRef->Hgt) RPN_FieldReadComponent(head,&Field->GRef->Hgt,"ZH",0,0);
            break;

         case 'X':
         case 'O':
            if (!Field->GRef->AY) RPN_FieldReadComponent(head,&Field->GRef->AY,"^^",1,0);
            if (!Field->GRef->AX) RPN_FieldReadComponent(head,&Field->GRef->AX,">>",1,0);
            GeoRef_BuildIndex(Field->GRef);          
            break;

          case 'V':
            if (!Field->GRef->AY) RPN_FieldReadComponent(head,&Field->GRef->AY,"^^",1,0);
            if (!Field->GRef->AX) RPN_FieldReadComponent(head,&Field->GRef->AX,">>",1,0);
            RPN_FieldReadLevels(Field);
            break;
      }

      // Need to re-qualify to check AX order
      GeoRef_Qualify(Field->GRef);
   }

   return(Field->GRef->AY && Field->GRef->AX);
}

/*----------------------------------------------------------------------------
 * Nom      : <RPN_FieldReadLevels>
 * Creation : Aout 2014 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Lire la coordonnee verticale pour un profile/xsection.
 *
 * Parametres :
 *  <Field>   : Adresse des valeurs du champs
 *
 * Retour:
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
int RPN_FieldReadLevels(TRPNField *Field) {

   int nj=0;

   if (Field->ZRef->Levels)
      free(Field->ZRef->Levels);

   Field->ZRef->Levels=NULL;
//TODO   Field->ZRef->LevelNb=RPN_FieldReadComponent(&Field->Head,&Field->ZRef->Levels,Field->Spec->Extrude,1,0);

   // If we don't find any level definition, use level index
   if (!Field->ZRef->Levels) {
      if (!(Field->ZRef->Levels=(float*)malloc(Field->Def->NJ*sizeof(float)))) {
         App_Log(ERROR,"%s: Not enough memory to read coordinates fields\n",__func__);
      } else {
         for(nj=0;nj<Field->Def->NJ;nj++) {
            Field->ZRef->Levels[nj]=nj+1;
         }
         Field->ZRef->LevelNb=Field->Def->NJ;
      }
   }
   return(nj);
}

int RPN_FieldWrite(int FileId,TRPNField *Field) {

   int  ok;

   c_fst_data_length(TDef_Size[Field->Def->Type]);
   ok=cs_fstecr(Field->Def->Data[0],-Field->Head.NBITS,FileId,Field->Head.DATEO,Field->Head.DEET,Field->Head.NPAS,
      Field->Head.NI,Field->Head.NJ,Field->Head.NK,Field->Head.IP1,Field->Head.IP2,Field->Head.IP3,Field->Head.TYPVAR,
      Field->Head.NOMVAR,Field->Head.ETIKET,Field->Head.GRTYP,Field->Head.IG1_JP,Field->Head.IG2_JP,Field->Head.IG3_JP,Field->Head.IG4_JP,Field->Head.DATYP,FALSE);

   if (ok<0) {
      App_Log(ERROR,"%s: Could not write field data (c_fstecr failed)\n",__func__);
      return(0);
   }
   return(1);
}

void RPN_CopyHead(TRPNHeader *To,TRPNHeader *From) {

   strncpy(To->NOMVAR,From->NOMVAR,5);
   strncpy(To->TYPVAR,From->TYPVAR,3);
   strncpy(To->ETIKET,From->ETIKET,13);
   To->DATEO=From->DATEO;
   To->DATEV=From->DATEV;
   To->DEET=From->DEET;
   To->NPAS=From->NPAS;
   To->DATYP=From->DATYP;
   To->IP1=From->IP1;
   To->IP2=From->IP2;
   To->IP3=From->IP3;
}

/*----------------------------------------------------------------------------
 * Nom      : <RPN_CopyDesc>
 * Creation : Janvier 2008 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Copier les descripteur de grille d'un fichier dans un autre
 *
 * Parametres :
 *   <FidTo>  : Fichier dans lequel copier
 *   <H>      : RPN Header
 *
 * Retour:
 *  <int>        : Code de reussite (0=erreur, 1=ok)
 *
 * Remarques :
 *----------------------------------------------------------------------------
*/
int RPN_CopyDesc(int FIdTo,TRPNHeader* const H) {

   TRPNHeader  h;
   char       *data=NULL;
   const char *desc;
   int         d=0,ni,nj,nk,sz=0,ip1,ip2;
   int         key;

   if (H->FID>-1) {
      strcpy(h.NOMVAR,"    ");
      strcpy(h.TYPVAR,"  ");
      strcpy(h.ETIKET,"            ");
      strcpy(h.GRTYP," ");

      pthread_mutex_lock(&RPNFieldMutex);
      while((desc=RPN_Desc[d++])) {
         if (strncmp(desc,"HY  ",2)==0) {
            ip1=-1;ip2=-1;
         } else {
            ip1=H->IG1_JP;
            ip2=H->IG2_JP;
         }
         key=c_fstinf(FIdTo,&ni,&nj,&nk,-1,"",ip1,ip2,-1,"",desc);
         if (key<0) {
            // If not already existing in destination
            key=c_fstinf(H->FID,&ni,&nj,&nk,-1,"",ip1,ip2,-1,"",desc);
            if (key>=0) {
               // If existing in source
               if (ni*nj>sz) {
                  // Some descriptors are 64 bit so always allocate double buffer
                  data=(char*)realloc(data,ni*nj*sizeof(double));
                  sz=ni*nj;
               }
               c_fstluk(data,key,&ni,&nj,&nk);

               key=c_fstprm(key,&h.DATEO,&h.DEET,&h.NPAS,&h.NI,&h.NJ,&h.NK,&h.NBITS,&h.DATYP,
                  &h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,&h.IG1_JP,
                  &h.IG2_JP,&h.IG3_JP,&h.IG4_JP,&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
               key=c_fstecr(data,NULL,-h.NBITS,FIdTo,h.DATEO,h.DEET,h.NPAS,h.NI,h.NJ,h.NK,h.IP1,
                  h.IP2,H->GRTYP[0]=='#'?0:h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,h.IG1_JP,h.IG2_JP,h.IG3_JP,h.IG4_JP,h.DATYP,1);
            } else if (key==-29) {
               // Input file not openned
               return(FALSE);
            }
         }
      }

      pthread_mutex_unlock(&RPNFieldMutex);
      if (data) free(data);
   } else {
      return(FALSE);
   }

   return(TRUE);
}

int RPN_IsDesc(const char* restrict Var) {

   const char *desc;
   int         d=0;

   while((desc=RPN_Desc[d++])) {
      if (!strncmp(Var,desc,4)) {
         return(TRUE);
      }
   }

   return(FALSE);
}

/*----------------------------------------------------------------------------
 * Nom      : <RPN_FieldTile>
 * Creation : Janvier 2008 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Save an RPN field as tiles
 *
 * Parametres  :
 *   <FID>     : File id
 *   <Def>     : Data definition
 *   <Head>    : RPN field Header
 *   <GRef>    : Georeference horizontale
 *   <ZRef>    : Georeference verticale
 *   <Comp>    : Component into data array
 *   <NI>      : Horizontal tile size
 *   <NJ>      : Vertical tile size
 *   <Halo>    : Width of the alo
 *   <DATYP>   : Data type (RPN value)
 *   <NPack>   : Facteur de compaction
 *   <Rewrite> : Reecrire le champs ou pas
 *
 * Retour:
 *  <int>        : Code de reussite (0=erreur, 1=ok)
 *
 * Remarques :
 *----------------------------------------------------------------------------
*/
int RPN_FieldTile(int FID,TDef *Def,TRPNHeader *Head,TGeoRef *GRef,TZRef *ZRef,int Comp,int NI,int NJ,int Halo,int DATYP,int NPack,int Rewrite,int Compress) {

   char        *tile=NULL,*data=NULL;;
   int          i,j,k,ip1,ni,nj,di,dj,pj,no,sz,key=0;
   unsigned int idx;

   // Allocate temp tile
   sz=TDef_Size[Def->Type];
   if (!(tile=(char*)malloc((NI+Halo*2)*(NJ+Halo*2)*sz))) {
      return(0);
   }

   pthread_mutex_lock(&RPNFieldMutex);

   for(k=0;k<Def->NK;k++) {
      idx=k*FSIZE2D(Def);

      Def_Pointer(Def,Comp,idx,data);

      // If IP1 is set, use it otherwise, convert it from levels array
      if ((ip1=Head->IP1)==-1 || Def->NK>1) {
         ip1=ZRef_Level2IP(ZRef->Levels[k],ZRef->Type,DEFAULT);
      }

      // Check if tiling asked and if dimensions allow tiling
      if (!NI || !NJ || (Def->NI<NI && Def->NJ<NJ)) {
         c_fst_data_length(TDef_Size[Def->Type]);
         key=c_fstecr(data,NULL,NPack,FID,Head->DATEO,Head->DEET,Head->NPAS,Def->NI,Def->NJ,1,ip1,Head->IP2,Head->IP3,Head->TYPVAR,
            Head->NOMVAR,Head->ETIKET,(GRef?(GRef->Grid[1]!='\0'?&GRef->Grid[1]:GRef->Grid):"X"),Head->IG1_JP,Head->IG2_JP,Head->IG3_JP,Head->IG4_JP,DATYP,Rewrite);
      } else {

         // Build and save the tiles, we adjust the tile size if it is too big
         no=0;
         for(j=0;j<Def->NJ;j+=NJ) {
            nj=((j+NJ>Def->NJ)?(Def->NJ-j):NJ)+Halo*2;
            dj=j-Halo;

            if (dj<0)          { dj+=Halo; nj-=Halo; }
            if (dj+nj>Def->NJ) { nj-=Halo; }

            for(i=0;i<Def->NI;i+=NI) {
               no++;
               ni=((i+NI>Def->NI)?(Def->NI-i):NI)+Halo*2;
               di=i-Halo;

               if (di<0)          { di+=Halo; ni-=Halo; }
               if (di+ni>Def->NI) { ni-=Halo; }

               for(pj=0;pj<nj;pj++) {
                  memcpy(tile+(pj*ni*sz),data+((dj+pj)*Def->NI+di)*sz,ni*sz);
               }
               c_fst_data_length(TDef_Size[Def->Type]);
               key=c_fstecr(tile,NULL,NPack,FID,Head->DATEO,Head->DEET,Head->NPAS,ni,nj,1,ip1,Head->IP2,no,Head->TYPVAR,
                  Head->NOMVAR,Head->ETIKET,"#",Head->IG1_JP,Head->IG2_JP,di+1,dj+1,DATYP,Rewrite);
            }
         }
      }
   }

   pthread_mutex_unlock(&RPNFieldMutex);

   free(tile);

   return(key>=0);
}


/*----------------------------------------------------------------------------
 * Nom      : <RPN_GetAllFields>
 * Creation : Mars 2015 - E. Legault-Ouellet - CMC/CMOE
 *
 * But      : Retourne une liste des champs correspondants aux critères donnés
 *
 * Parametres :
 *  <FID>   : Le file handle vers un fichier standard
 *  <DateV> : La date de validité des champs à chercher
 *  <Etiket>: L'etiket des champs à chercher
 *  <Ip1>   : L'ip1 des champs à chercher
 *  <Ip2>   : L'ip2 des champs à chercher
 *  <Ip3>   : L'ip3 des champs à chercher
 *  <Typvar>: Le typvar des champs à chercher
 *  <Nomvar>: Le nomvar des champs à chercher
 *  <Arr>   : Le pointeur retourné vers la liste des champs
 *  <Size>  : Nombre d'items de la liste
 *
 * Retour   : APP_ERR si erreur, APP_OK si ok.
 *
 * Remarques :
 *  La mémoire est allouée dans la fonction et a besoin d'être libérée par la
 *  fonction appelante. En cas d'erreur, le pointeur retourné sera NULL.
 *
 *----------------------------------------------------------------------------
 */
int RPN_GetAllFields(int FID,int DateV,char *Etiket,int Ip1,int Ip2,int Ip3,char *Typvar,char *Nomvar,int **Arr,int *Size) {
    int ni,nj,nk,err;
    int n=32768; // ~128kiB
    int *arr=NULL,s;

    *Arr    = NULL;
    *Size   = 0;

    do {
        APP_MEM_ASRT(arr,realloc(arr,(n*=2)*sizeof(*arr)));
        if( (err=cs_fstinl(FID,&ni,&nj,&nk,DateV,Etiket,Ip1,Ip2,Ip3,Typvar,Nomvar,arr,&s,n)) && err!=-4762 ) {
            App_Log(ERROR,"(%s) Couldn't get list of fields (cs_fstinl) code=%d\n",__func__,s,n,err);
            APP_FREE(arr);
            return(APP_ERR);
        }
    } while( n == s );

    // Only take the memory we really need
    //APP_MEM_ASRT( arr,realloc(arr,s*sizeof(*arr)) );

    *Size   = s;
    *Arr    = arr;

    return(APP_OK);
}

/*----------------------------------------------------------------------------
 * Nom      : <RPN_GetAllDates>
 * Creation : Mars 2015 - E. Legault-Ouellet - CMC/CMOE
 *
 * But      : Retourne une liste des DateV correspondants aux critères donnés
 *
 * Parametres :
 *  <Flds>  : Les champs dont on veut les dates
 *  <NbFlds>: Le nombre des champs dont on veut les dates
 *  <Uniq>  : Boolean. If true, only returns the sorted list of unique dates.
 *            If false, returns the list of valid dates associated with (and in
 *            the same order as) the list of fields.
 *  <DateV> : [OUT] Les dates valides des champs
 *  <NbDateV: [OUT] Le nombre de dates valides retournées
 *
 * Retour   : APP_ERR si erreur, APP_OK si ok.
 *
 * Remarques :
 *  La mémoire est allouée dans la fonction et a besoin d'être libérée par la
 *  fonction appelante. En cas d'erreur, le pointeur retourné sera NULL.
 *
 *----------------------------------------------------------------------------
 */
int RPN_GetAllDates(int *Flds,int NbFlds,int Uniq,int **DateV,int *NbDateV) {
    TRPNHeader  h;
    int         i,err,*dates;
    double      deltat;

    *DateV = NULL;
    *NbDateV = 0;
    APP_MEM_ASRT(dates,malloc(NbFlds*sizeof(*dates)));

    for(i=0; i<NbFlds; ++i) {
        err=cs_fstprm(Flds[i],&h.DATEO,&h.DEET,&h.NPAS,&h.NI,&h.NJ,&h.NK,&h.NBITS,&h.DATYP,&h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,
                h.GRTYP,&h.IG1_JP,&h.IG2_JP,&h.IG3_JP,&h.IG4_JP,&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
        if( err ) {
            App_Log(ERROR,"(RPN_GetAllDates) Couldn't get info on field (cs_fstprm)\n");
            APP_FREE(dates);
            return(APP_ERR);
        }

        deltat = h.DEET*h.NPAS/3600.0;
        f77name(incdatr)(&dates[i],&h.DATEO,&deltat);
        if( dates[i] == 101010101 ) {
            App_Log(ERROR,"(RPN_GetAllDates) Couldn't get DateV for dateo(%d),deet(%d),npas(%d),deltat(%f) (incdatr)\n",h.DATEO,h.DEET,h.NPAS,deltat);
            APP_FREE(dates);
            return(APP_ERR);
        }
    }

    if( Uniq ) {
        qsort(dates,NbFlds,sizeof(*dates),QSort_Int);
        Unique(dates,&NbFlds,sizeof(*dates));
        // Only take the memory we really need
        //APP_MEM_ASRT( dates,realloc(dates,NbFlds*sizeof(*dates)) );
    }

    *DateV = dates;
    *NbDateV = NbFlds;
    return(APP_OK);
}

/*----------------------------------------------------------------------------
 * Nom      : <RPN_GetAllIps>
 * Creation : Mars 2015 - E. Legault-Ouellet - CMC/CMOE
 *
 * But      : Retourne une liste des IPx correspondants aux critères donnés
 *
 * Parametres :
 *  <Flds>  : Les champs dont on veut les IPx
 *  <NbFlds>: Le nombre des champs dont on veut les IPx
 *  <IpN>   : Le numéro de l'ip voulu (1 for IP1, 2 for IP2, 3 for IP3)
 *  <Uniq>  : Boolean. If true, only returns the sorted list of unique IP1.
 *            If false, returns the list of IP1 associated with (and in
 *            the same order as) the list of fields.
 *  <Ips>   : [OUT] Les IPx des champs
 *  <NbIp>  : [OUT] Le nombre d'IPx retournés
 *
 * Retour   : APP_ERR si erreur, APP_OK si ok.
 *
 * Remarques :
 *  La mémoire est allouée dans la fonction et a besoin d'être libérée par la
 *  fonction appelante. En cas d'erreur, le pointeur retourné sera NULL.
 *
 *----------------------------------------------------------------------------
 */
int RPN_GetAllIps(int *Flds,int NbFlds,int IpN,int Uniq,int **Ips,int *NbIp) {
    TRPNHeader  h;
    int         i,err,*ips;
    double      deltat;

    *Ips = NULL;
    *NbIp = 0;
    APP_MEM_ASRT(ips,malloc(NbFlds*sizeof(*ips)));

    for(i=0; i<NbFlds; ++i) {
        err=cs_fstprm(Flds[i],&h.DATEO,&h.DEET,&h.NPAS,&h.NI,&h.NJ,&h.NK,&h.NBITS,&h.DATYP,&h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,
                h.GRTYP,&h.IG1_JP,&h.IG2_JP,&h.IG3_JP,&h.IG4_JP,&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
        if( err ) {
            App_Log(ERROR,"(RPN_GetAllIps) Couldn't get info on field (cs_fstprm)\n");
            APP_FREE(ips);
            return(APP_ERR);
        }

        switch( IpN ) {
            case 1: ips[i]=h.IP1; break;
            case 2: ips[i]=h.IP2; break;
            case 3: ips[i]=h.IP3; break;
            default:
                App_Log(ERROR,"(RPN_GetAllIps) [%d] is not a valid IP number. Valid numbers are 1,2 and 3.\n",IpN);
                free(ips);
                return(APP_ERR);
        }
    }

    if( Uniq ) {
        qsort(ips,NbFlds,sizeof(*ips),QSort_Int);
        Unique(ips,&NbFlds,sizeof(*ips));
        //APP_MEM_ASRT( ips,realloc(ips,NbFlds*sizeof(*ips)) );
    }

    *Ips = ips;
    *NbIp = NbFlds;
    return(APP_OK);
}

/*----------------------------------------------------------------------------
 * Nom      : <RPN_GenerateIG>
 * Creation : Septembre 2016 - E. Legault-Ouellet - CMC/CMOE
 *
 * But      : Initialise des IG1_JP/2/3 à des valeurs qui se veulent uniques
 *
 * Parametres :
 *  <IG1_JP>   : [OUT] IG1_JP à Initialiser
 *  <IG2_JP>   : [OUT] IG2_JP à Initialiser
 *  <IG3_JP>   : [OUT] IG3_JP à Initialiser
 *
 * Retour   : APP_ERR si erreur, APP_OK si ok.
 *
 * Remarques :
 *  La mémoire est allouée dans la fonction et a besoin d'être libérée par la
 *  fonction appelante. En cas d'erreur, le pointeur retourné sera NULL.
 *
 *----------------------------------------------------------------------------
 */
int RPN_GenerateIG(int *IG1_JP,int *IG2_JP,int *IG3_JP) {
    int64_t bits = (int64_t)time(NULL);

    // IG1_JP, IG2_JP and IG3_JP are limited to 24 bits, but we only have 64 bits to distribute, si IG3_JP will only receive up to 16 bits
    *IG2_JP = bits&0xffff; bits>>=16;
    *IG1_JP = bits&0xffffff; bits>>=24;
    *IG3_JP = bits; //Basically 0. Should turn 1 on January 19th, 2038 (so still way before my retirement)

    return(APP_OK);
}

/*----------------------------------------------------------------------------
 * Nom      : <RPN_LinkFiles>
 * Creation : Octobre 2016 - E. Legault-Ouellet - CMC/CMOE
 *
 * But      : Ouvre et lie plusieurs fichiers standards ensemble
 *
 * Parametres :
 *  <Files> : Liste des fichiers standards à ouvrir et lier
 *  <N>     : Nombre de fichiers à lier (-1 si la liste des fichiers est NULL-terminated)
 *
 * Retour   : Le FID à utiliser ou un nombre négatif si erreur.
 *
 * Remarques : Cette fonction N'EST PAS thread safe
 *  
 *----------------------------------------------------------------------------
 */
int RPN_LinkFiles(char **Files,int N) {
    // Make sure there is at least one file
    if( !Files || N==0 || !*Files )
        return -1;

    // Count the number of files
    if( N<0 ) {
        for(N=1; Files[N]; ++N);
    }

    // Only apply the special treatment if there actually is more than one file
    if( N>1 ) {
        int i,*lst,**ptr;

        // Allocate the memory
        lst = malloc((N+1)*sizeof(*lst));

        // Put the size first
        lst[0] = N;

        // Open all the files
        for(i=0; Files[i]; ++i) {
            if( (lst[i+1]=cs_fstouv(Files[i],"STD+RND+R/O")) < 0 ) {
                App_Log(ERROR,"(%s) Problem opening input file \"%s\"\n",__func__,Files[i]);
                free(lst);
                return -1;
            }
        }

        // Link all files
        if( f77name(fstlnk)(lst+1,lst) != 0 ) {
            App_Log(ERROR,"(%s) Could not link the %d input files together\n",__func__,N);
            free(lst);
            return -1;
        }

        // Add the list of FIDs to the global list
        if( !(ptr=realloc(LNK_FID,(LNK_NB+1)*sizeof(*LNK_FID))) ) {
            App_Log(ERROR,"(%s) Could not allocate memory for the global FIDs array\n",__func__);
            free(lst);
            return -1;
        }
        LNK_FID = ptr;
        LNK_FID[LNK_NB++] = lst;

        // Return the file handle
        return lst[1];
    } else {
        int h;

        // Only one file, no need to link anything
        if( (h=cs_fstouv(Files[0],"STD+RND+R/O")) < 0 ) {
            App_Log(ERROR,"(%s) Problem opening input file %s\n",__func__,Files[0]);
        }

        return h;
    }
}

/*----------------------------------------------------------------------------
 * Nom      : <RPN_UnLinkFiles>
 * Creation : Septembre 2016 - E. Legault-Ouellet - CMC/CMOE
 *
 * But      : Délie et ferme les fichiers standards liés
 *
 * Parametres :
 *  <FID>   : Le FID des fichiers liés à libérer
 *
 * Retour   : APP_ERR si erreur, APP_OK si ok.
 *
 * Remarques : Cette fonction N'EST PAS thread safe
 * 
 *----------------------------------------------------------------------------
 */
int RPN_UnLinkFiles(int FID) {
    // Make sure we have a valid handle
    if( FID >= 0 ) {
        int i,j,**ptr;

        // Find the list containing that FID in the first position
        for(i=0; i<LNK_NB; ++i) {
            if( LNK_FID[i][1] == FID ) {
                // Unlink the files (the number of files is the header of that list)
                f77name(fstunl)(LNK_FID[i]+1,LNK_FID[i]);

                // Close the files
                for(j=1; j<LNK_FID[i][0]; ++j)
                    cs_fstfrm(LNK_FID[i][j]);

                // Free the memory
                free(LNK_FID[i]);

                // Move back the remaining FIDs
                for(j=i+1; j<LNK_NB; )
                    LNK_FID[i++] = LNK_FID[j++];

                // Resize the table
                --LNK_NB;
                if( LNK_NB ) {
                    if( (ptr=realloc(LNK_FID,LNK_NB*sizeof(*LNK_FID))) ) {;
                        LNK_FID = ptr;
                    }
                } else {
                    free(LNK_FID);
                    LNK_FID = NULL;
                }

                return APP_OK;
            }
        }

        // If we are still here, it means there was no list (only one file). Therefore, we just close the file.
        cs_fstfrm(FID);
    }

    return APP_OK;
}

/*----------------------------------------------------------------------------
 * Nom      : <RPN_LinkPattern>
 * Creation : Avril 2017 - E. Legault-Ouellet - CMC/CMOE
 *
 * But      : Ouvre et lie les fichiers correspondant au pattern donné en
 *            mode lecture seulement
 *
 * Parametres :
 *    <Pattern>   : Le pattern des fichiers à ouvrir et lier.
 *
 * Retour   : Le FID à utiliser ou un nombre négatif si erreur.
 *
 * Remarques : Cette fonction N'EST PAS thread safe
 *
 *----------------------------------------------------------------------------
 */
int RPN_LinkPattern(const char* Pattern) {
   // Expand the pattern into a list of files
   glob_t gfiles = (glob_t){0,NULL,0};
   if( glob(Pattern,0,NULL,&gfiles) || !gfiles.gl_pathc ) {
      return -1;
   }

   int fid = RPN_LinkFiles(gfiles.gl_pathv,(int)gfiles.gl_pathc);
   globfree(&gfiles);
   return fid;
}

#endif
