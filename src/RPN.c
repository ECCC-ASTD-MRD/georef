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

int cs_fstprm(int Idx,int *DateO,int *Deet,int *NPas,int *NI,int *NJ,int *NK,int *NBits,int *Datyp,int *IP1,int *IP2,int *IP3,char* TypVar,char *NomVar,char *Etiket,char *GrTyp,int *IG1,int *IG2,int *IG3,int *IG4,int *Swa,int *Lng,int *DLTF,int *UBC,int *EX1,int *EX2,int *EX3) {

   int err;

   pthread_mutex_lock(&RPNFieldMutex);
   err=c_fstprm(Idx,DateO,Deet,NPas,NI,NJ,NK,NBits,Datyp,IP1,IP2,IP3,TypVar,NomVar,Etiket,GrTyp,IG1,IG2,IG3,IG4,Swa,Lng,DLTF,UBC,EX1,EX2,EX3);
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
//   TGrid *grid;
// TODO: reenable this
   if (GRTYP[0]=='#') {
//      if ((grid=EZGrid_ReadIdx(Unit,Idx,0))) {
//         if (EZGrid_TileBurnAll(grid,0,Data)) {
//            err=0;
//         }
//         EZGrid_Free(grid);
//      }
//      *NI=grid->H.NI;
//      *NJ=grid->H.NJ;
//      *NK=grid->H.NK;
   } else {
      pthread_mutex_lock(&RPNFieldMutex);
      err=c_fstluk(Data,Idx,NI,NJ,NK);
      pthread_mutex_unlock(&RPNFieldMutex);
   }

   return(err);
}

int cs_fstecr(void *Data,int NPak,int Unit, int DateO,int Deet,int NPas,int NI,int NJ,int NK,int IP1,int IP2,int IP3,char* TypVar,char *NomVar,char *Etiket,char *GrTyp,int IG1,int IG2,int IG3,int IG4,int DaTyp,int Over) {
   int err;

   pthread_mutex_lock(&RPNFieldMutex);
   err=c_fstecr(Data,NULL,NPak,Unit,DateO,Deet,NPas,NI,NJ,NK,IP1,IP2,IP3,TypVar,NomVar,Etiket,GrTyp,IG1,IG2,IG3,IG4,DaTyp,Over);
   pthread_mutex_unlock(&RPNFieldMutex);

   return(err);
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
   h.GRREF[0]=h.GRREF[1]='\0';

   strcpy(h.NOMVAR,"    ");
   strcpy(h.TYPVAR,"  ");
   strcpy(h.ETIKET,"            ");
   ok=cs_fstprm(h.KEY,&h.DATEO,&h.DEET,&h.NPAS,&h.NI,&h.NJ,&h.NK,&h.NBITS,
         &h.DATYP,&h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,
         h.GRTYP,&h.IG[X_IG1],&h.IG[X_IG2],&h.IG[X_IG3],&h.IG[X_IG4],&h.SWA,&h.LNG,&h.DLTF,
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
   // TODO: Should we do this automatically or only when needed
      // Recuperer les type de niveaux et forcer ETA pour SIGMA
      lvl=ZRef_IP2Level(h.IP1,&type);
      type=type==LVL_SIGMA?LVL_ETA:type;

      fld->GRef=GeoRef_Create(h.NI,h.NJ,h.GRTYP,h.IG[X_IG1],h.IG[X_IG2],h.IG[X_IG3],h.IG[X_IG4],h.FID);
      
      fld->ZRef=ZRef_Define(type,h.NK,&lvl);
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
 *  <GRREF>   : Grid reference returned (Optional)
 * 
 * Retour:
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
*/
int RPN_FieldReadComponent(TRPNHeader *Head,float **Ptr,char *Var,int Grid,int Force) {

   int key=0,ni=0,nj=0,nk=0,i;
   char c[14];

   if (Var && (!*Ptr || Force)) {
      if (Grid==1) {
         // Look for corresponding time and if not use any time
         if ((key==cs_fstinf(Head->FID,&ni,&nj,&nk,Head->DATEV,"",Head->IG[X_IG1],Head->IG[X_IG2],Head->IG[X_IG3],"",Var))<=0) {
            key=cs_fstinf(Head->FID,&ni,&nj,&nk,-1,"",Head->IG[X_IG1],Head->IG[X_IG2],Head->IG[X_IG3],"",Var);
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
 //           cs_fstprm(key,&i,&i,&i,&i,&i,&i,&Head->NBITS,&i,&i,&i,&i,c,c,c,Head->GRREF,&Head->IGREF[X_IG1],&Head->IGREF[X_IG2],&Head->IGREF[X_IG3],&Head->IGREF[X_IG4],&i,&i,&i,&i,&i,&i,&i);
            if (!(*Ptr=(float*)malloc(ni*nj*nk*sizeof(float)))) {
               App_Log(ERROR,"%s: Not enough memory to read coordinates fields\n",__func__);
               return(0);
            }
         }
         cs_fstluk(*Ptr,key,&ni,&nj,&nk);
         if (Force || (Grid && (Head->GRREF[0]==' ' || Head->GRREF[0]=='\0'))) {
            cs_fstprm(key,&i,&i,&i,&i,&i,&i,&i,&i,&i,&i,&i,c,c,c,Head->GRREF,&Head->IGREF[X_IG1],&Head->IGREF[X_IG2],&Head->IGREF[X_IG3],&Head->IGREF[X_IG4],&i,&i,&i,&i,&i,&i,&i);
         }
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
struct TGeoRef* RPN_FieldReadGrid(TRPNField *Field) {

   Field->GRef=GeoRef_Create(Field->Head.NI,Field->Head.NJ,Field->Head.GRTYP,Field->Head.IG[X_IG1],Field->Head.IG[X_IG2],Field->Head.IG[X_IG3],Field->Head.IG[X_IG4],Field->Head.FID);
   
   if (Field->GRef->GRTYP[0]=='V') {
      RPN_FieldReadLevels(Field);
   } 
   return(Field->GRef);
}

int RPN_ReadGrid(struct TGeoRef *GRef) {

   int         key,ni,nj,nk,ig1,ig2,ig3,ig4,idx,s,i,j,offsetx,offsety,dx=0,dy=0;
   float      *ax=NULL,*ay=NULL;
   char        grref[2];
   TRPNHeader *h=&GRef->RPNHead;

   if (GRef->GRTYP[0]=='L' || GRef->GRTYP[0]=='A' || GRef->GRTYP[0]=='B' || GRef->GRTYP[0]=='N' || GRef->GRTYP[0]=='S' || GRef->GRTYP[0]=='G') {
      return(TRUE);
   }

   if ((!GRef->AY || !GRef->AX) && h->FID>=0) {

      switch(GRef->GRTYP[0]) {
         case 'M':
            if (!GRef->AY) dx=RPN_FieldReadComponent(h,&ay,"^^",1,0);
            if (!GRef->AX) dy=RPN_FieldReadComponent(h,&ax,">>",1,0);

            if (ax) { GRef->AX=(double*)malloc(dx*sizeof(double)); for(i=0;i<dx;i++) GRef->AX[i]=ax[i]; }
            if (ay) { GRef->AY=(double*)malloc(dy*sizeof(double)); for(i=0;i<dy;i++) GRef->AY[i]=ay[i]; }

            // Lire le champs d'indexes
            if (!GRef->Idx) {
               key=cs_fstinf(h->FID,&ni,&nj,&nk,-1,"",h->IG[X_IG1],h->IG[X_IG2],h->IG[X_IG3],"","##");
               if (key < 0) {
                  App_Log(ERROR,"%s: Could not find index field %s (c_fstinf failed)",__func__,"##");
                  return(FALSE);
               } else {
                  GRef->NIdx=ni*nj*nk;
                  if (!(GRef->Idx=(unsigned int*)malloc(GRef->NIdx*sizeof(unsigned int)))) {
                     App_Log(ERROR,"%s: Not enough memory to read coordinates fields",__func__);
                     return(FALSE);
                  }
                  cs_fstluk((float*)GRef->Idx,key,&ni,&nj,&nk);
               }
               GeoRef_BuildIndex(GRef);
            }
            break;

         case 'W':
            break;

         case 'Y':
            if (!GRef->AY)  RPN_FieldReadComponent(h,&ay,"LA",0,0);
            if (!ay)        RPN_FieldReadComponent(h,&ay,"^^",1,0);
            if (!GRef->AX)  RPN_FieldReadComponent(h,&ax,"LO",0,0);
            if (!ax)        RPN_FieldReadComponent(h,&ax,">>",1,0);
            if (!GRef->Hgt) RPN_FieldReadComponent(h,&GRef->Hgt,"ZH",0,0);

            if (ax) { GRef->AX=(double*)malloc(GRef->NX*sizeof(double)); for(i=0;i<GRef->NX;i++) GRef->AX[i]=ax[i]; }
            if (ay) { GRef->AY=(double*)malloc(GRef->NY*sizeof(double)); for(i=0;i<GRef->NY;i++) GRef->AY[i]=ay[i]; }
           break;

         case 'X':
         case 'O':
            if (!GRef->AY) dy=RPN_FieldReadComponent(h,&ay,"^^",1,0);
            if (!GRef->AX) dx=RPN_FieldReadComponent(h,&ax,">>",1,0);

            if (ax) { GRef->AX=(double*)malloc(dx*sizeof(double)); for(i=0;i<dx;i++) GRef->AX[i]=ax[i]; }
            if (ay) { GRef->AY=(double*)malloc(dy*sizeof(double)); for(i=0;i<dy;i++) GRef->AY[i]=ay[i]; }
            GeoRef_BuildIndex(GRef);          
            break;

         case 'V':
            if (!GRef->AY) RPN_FieldReadComponent(h,&ay,"^^",1,0);
            if (!GRef->AX) RPN_FieldReadComponent(h,&ax,">>",1,0);

            if (ax) { GRef->AX=(double*)malloc(GRef->NX*sizeof(double)); for(i=0;i<GRef->NX;i++) GRef->AX[i]=ax[i]; }
            if (ay) { GRef->AY=(double*)malloc(GRef->NY*sizeof(double)); for(i=0;i<GRef->NY;i++) GRef->AY[i]=ay[i]; }
//TODO:            RPN_FieldReadLevels(Field);
            break;

         case 'Z':
            if (!GRef->AY) dy=RPN_FieldReadComponent(h,&ay,"^^",1,0);
            if (!GRef->AX) dx=RPN_FieldReadComponent(h,&ax,">>",1,0);

            if (ax) { GRef->AX=(double*)malloc(dx*sizeof(double)); for(i=0;i<dx;i++) GRef->AX[i]=ax[i]; }
            if (ay) { GRef->AY=(double*)malloc(dy*sizeof(double)); for(i=0;i<dy;i++) GRef->AY[i]=ay[i]; }
            break;

         case '#':
            if (!GRef->AY) RPN_FieldReadComponent(h,&ay,"^^",1,0);
            if (!GRef->AX) RPN_FieldReadComponent(h,&ax,">>",1,0);

            if (ax) { GRef->AX=(double*)malloc(GRef->NX*sizeof(double)); for(i=0;i<GRef->NX;i++) GRef->AX[i]=ax[i]; }
            if (ay) { GRef->AY=(double*)malloc(GRef->NY*sizeof(double)); for(i=0;i<GRef->NY;i++) GRef->AY[i]=ay[i]; }

            if (ax && ay) {
               GRef->AX = (double*) malloc(GRef->NX*sizeof(double));
               GRef->AY = (double*) malloc(GRef->NY*sizeof(double));
   // TODO: Check with LireEnrPositionel
   //            offsetx = ip3 - 1;
   //            offsety = ip4 - 1;
   //            for (j=0; j < GRef->NY; j++) GRef->AY[j] = ay[j+offsety];
   //            for (i=0; i < GRef->NX; i++) GRef->AX[i] = ax[i+offsetx];
            }
            break;
         
         case 'U':
            if (!GRef->AX) RPN_FieldReadComponent(h,&ax,"^>",1,0);

            if (ax) {
               GRef->NbSub=(int)ax[2];            // Number of LAM grids (YY=2)
               ni=(int)ax[5];                     // NI size of LAM grid 
               nj=(int)ax[6];                     // NJ size of LAM grid
               GRef->NX = ni;                     // NI size of U grid
               GRef->NY = nj*GRef->NbSub;         // NJ size of U grid
               GRef->AX = (double*)malloc(ni*sizeof(double));
               GRef->AY = (double*)malloc(nj*sizeof(double));
               for(i=0;i<ni;i++) GRef->AX[i]=ax[15+i];
               for(i=0;i<nj;i++) GRef->AY[i]=ax[15+ni+i];

               // Get subgrids
               GRef->Subs = (TGeoRef**)malloc(GRef->NbSub*sizeof(TGeoRef*));
               strcpy(grref,"E");
               idx=11;
               for(s=0;s<GRef->NbSub;s++) {
                  f77name(cxgaig)(grref,&ig1,&ig2,&ig3,&ig4,&ax[idx],&ax[idx+1],&ax[idx+2],&ax[idx+3]);
                  GRef->Subs[s] = GeoRef_CreateInMemory(ni,nj,"Z",grref,ig1,ig2,ig3,ig4,GRef->AX,GRef->AY);
                  //TODO: Do we need to do this here ?
                  GeoRef_MaskYYDefine(GRef->Subs[s]);
                  idx+=ni+nj+10;
               }
               free(ax);
            }
      }

      if (ax) free(ax);
      if (ay) free(ay);

      // Check for 2D AX/AY
      if (dx>GRef->RPNHead.NI) {
         GRef->Type|=GRID_AXY2D;
      }
      // In case of WKT REF
      if (GRef->GRTYP[0]=='W' || h->GRREF[0] == 'W') {
#ifdef HAVE_GDAL
         float  tmpv[6];
         double mtx[6],inv[6],*tm=NULL,*im=NULL;
         char  *proj=NULL;

         if ((key=cs_fstinf(h->FID,&ni,&nj,&nk,-1,"",h->IG[X_IG1],h->IG[X_IG2],h->IG[X_IG3],"","PROJ"))<0) {
            App_Log(ERROR,"%s: Could not find projection description field PROJ (c_fstinf failed)\n",__func__);
            return(FALSE);
         } else {
            c_fst_data_length(1);
            if ((proj=(char*)malloc(ni*nj*4))) {
               cs_fstluk(proj,key,&ni,&nj,&nk);
            }
         }
         if ((key=cs_fstinf(h->FID,&ni,&nj,&nk,-1,"",h->IG[X_IG1],h->IG[X_IG2],h->IG[X_IG3],"","MTRX"))<0) {
            App_Log(ERROR,"%s: Could not find trasform matrix field MTRX (c_fstinf failed)\n",__func__);
            return(FALSE);
         } else {
            cs_fstluk(tmpv,key,&ni,&nj,&nk);
            for(i=0;i<ni;i++) {
               mtx[i]=tmpv[i];
            }
            tm=mtx;
            if (!GDALInvGeoTransform(mtx,inv)) {
               im=NULL;
            } else {
               im=inv;
            }
         }
         GeoRef_SetW(GRef,proj,tm,im,NULL);
         if (proj) free(proj);
#else
   App_Log(ERROR,"W grid support not enabled, needs to be built with GDAL\n",__func__);
   return(FALSE);
#endif
      }
   }

   return(TRUE);
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
//TODO   Field->ZRef->LevelNb=RPN_FieldReadComponent(&Field->Head,&Field->ZRef->Levels,Field->Spec->Extrude,1,0,NULL);

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
      Field->Head.NOMVAR,Field->Head.ETIKET,Field->Head.GRTYP,Field->Head.IG[X_IG1],Field->Head.IG[X_IG2],Field->Head.IG[X_IG3],Field->Head.IG[X_IG4],Field->Head.DATYP,FALSE);

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
      strcpy(h.GRREF," ");

      pthread_mutex_lock(&RPNFieldMutex);
      while((desc=RPN_Desc[d++])) {
         if (strncmp(desc,"HY  ",2)==0) {
            ip1=-1;ip2=-1;
         } else {
            ip1=H->IG[X_IG1];
            ip2=H->IG[X_IG2];
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
                  &h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,&h.IG[X_IG1],
                  &h.IG[X_IG2],&h.IG[X_IG3],&h.IG[X_IG4],&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
               key=c_fstecr(data,NULL,-h.NBITS,FIdTo,h.DATEO,h.DEET,h.NPAS,h.NI,h.NJ,h.NK,h.IP1,
                  h.IP2,H->GRTYP[0]=='#'?0:h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,h.IG[X_IG1],h.IG[X_IG2],h.IG[X_IG3],h.IG[X_IG4],h.DATYP,1);
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
 *   <Ref>    : Georeference horizontale
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
int RPN_FieldTile(int FID,TDef *Def,TRPNHeader *Head,TGeoRef *Ref,TZRef *ZRef,int Comp,int NI,int NJ,int Halo,int DATYP,int NPack,int Rewrite,int Compress) {

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
            Head->NOMVAR,Head->ETIKET,(Ref?(Ref->GRTYP[1]!='\0'?&Ref->GRTYP[1]:Ref->GRTYP):"X"),Head->IG[X_IG1],Head->IG[X_IG2],Head->IG[X_IG3],Head->IG[X_IG4],DATYP,Rewrite);
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
                  Head->NOMVAR,Head->ETIKET,"#",Head->IG[X_IG1],Head->IG[X_IG2],di+1,dj+1,DATYP,Rewrite);
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
                h.GRTYP,&h.IG[X_IG1],&h.IG[X_IG2],&h.IG[X_IG3],&h.IG[X_IG4],&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
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
                h.GRTYP,&h.IG[X_IG1],&h.IG[X_IG2],&h.IG[X_IG3],&h.IG[X_IG4],&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
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
 * But      : Initialise des IG1/2/3 à des valeurs qui se veulent uniques
 *
 * Parametres :
 *  <IG1>   : [OUT] IG1 à Initialiser
 *  <IG2>   : [OUT] IG2 à Initialiser
 *  <IG3>   : [OUT] IG3 à Initialiser
 *
 * Retour   : APP_ERR si erreur, APP_OK si ok.
 *
 * Remarques :
 *  La mémoire est allouée dans la fonction et a besoin d'être libérée par la
 *  fonction appelante. En cas d'erreur, le pointeur retourné sera NULL.
 *
 *----------------------------------------------------------------------------
 */
int RPN_GenerateIG(int *IG1,int *IG2,int *IG3) {
    int64_t bits = (int64_t)time(NULL);

    // IG1, IG2 and IG3 are limited to 24 bits, but we only have 64 bits to distribute, si IG3 will only receive up to 16 bits
    *IG2 = bits&0xffff; bits>>=16;
    *IG1 = bits&0xffffff; bits>>=24;
    *IG3 = bits; //Basically 0. Should turn 1 on January 19th, 2038 (so still way before my retirement)

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
