#include <App.h>
#include "GeoRef_Utils.h"
#include "ZRef.h"

#ifdef HAVE_VGRID
#include "vgrid.h"
#endif

static float       *ZRef_Levels   = NULL;
static uint32_t ZRef_LevelsNb = 0;
static const char  *ZRef_Names[]  = { "MASL","SIGMA","PRESSURE","UNDEFINED","MAGL","HYBRID","THETA","MBSL","GALCHEN","COUNT","HOUR","ANGLE","NIL","NIL","NIL","INT","NIL","IDX","NIL","NIL","NIL","MPRES","NIL","NIL","NIL","NIL","NIL","NIL","NIL","NIL","NIL","NIL","ETA",NULL };
static const char  *ZRef_Units[]  = { "m","sg","mb","-","m","hy","th","m-","m","nb","hr","dg","--","--","--","i","--","x","--","--","--","mp","--","--","--","--","--","--","--","--","--","--","sg",NULL };

int32_t ZREF_IP1MODE=3;

const char **ZRef_LevelNames() {
   return(ZRef_Names);
}

const char *ZRef_LevelName(int32_t Type) {
   return(ZRef_Names[Type<0||Type>21?22:Type]);
}

const char **ZRef_LevelUnits() {
   return(ZRef_Units);
}

const char *ZRef_LevelUnit(int32_t Type) {
   return(ZRef_Units[Type<0||Type>21?22:Type]);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_New>
 * Creation : Octobre 2011 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Initialize a new vertical reference.
 *
 * Parametres  :
 *
 * Retour:
 *  <ZRef>     : Vertical reference (NULL on error)
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
TZRef* ZRef_New(void) {

   TZRef *zref=NULL;

   if (!(zref=(TZRef*)calloc(1,sizeof(TZRef)))) {
      Lib_Log(APP_LIBEER,APP_ERROR,"%s: Unable to allocate memory\n",__func__);
   }
   
   zref->VGD=NULL;
   zref->Levels=NULL;
   zref->Style=NEW;
   zref->Type=LVL_UNDEF;
   zref->LevelNb=0;
   zref->POff=zref->PTop=zref->PRef=zref->ETop=0.0;
   zref->RCoef[0]=zref->RCoef[1]=1.0;
   zref->P0=zref->P0LS=zref->PCube=zref->A=zref->B=NULL;
   zref->Version=-1;
   zref->NRef=1;
   
   return(zref);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_Define>
 * Creation : Octobre 2011 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Initialize a vertical reference.
 *
 * Parametres  :
 *  <ZRef>     : Vertical reference to initialize
 *  <Type>     : Type of levels
 *  <NbLevels> : Number of levels passed in
 *  <Levels>   : List of levels to assign
 *
 * Retour:
 *  <Ok>       : (1=Ok, 0=Bad)
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
TZRef* ZRef_Define(int32_t Type,int32_t NbLevels,float *Levels) {

   TZRef *zref=NULL;
   
   if ((zref=ZRef_New())) {

      zref->Type=Type;
      zref->LevelNb=NbLevels;
      
      if (zref->Levels && zref->LevelNb!=NbLevels) {
         zref->Levels=(float*)realloc(zref->Levels,zref->LevelNb*sizeof(float));
      } else {
         zref->Levels=(float*)calloc(zref->LevelNb,sizeof(float));
      }
      
      if (!zref->Levels) {
         Lib_Log(APP_LIBEER,APP_ERROR,"%s: Unable to allocate memory\n",__func__);
      } else if (Levels) {
         memcpy(zref->Levels,Levels,zref->LevelNb*sizeof(float));
      }
   }
   
   return(zref);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_Free>
 * Creation : Octobre 2011 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Free a vertical reference.
 *
 * Parametres  :
 *  <ZRef>     : Vertical reference to free
 *
 * Retour:
 *  <Ok>       : (1=Ok, 0=Bad)
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
int32_t ZRef_Free(TZRef *ZRef) {

   if (ZRef && !ZRef_Decr(ZRef)) {

#ifdef HAVE_VGRID
      if (ZRef->VGD)    Cvgd_free((vgrid_descriptor**)&ZRef->VGD); ZRef->VGD=NULL;
#endif
      if (ZRef->Levels) free(ZRef->Levels);    ZRef->Levels=NULL;
      if (ZRef->A)      free(ZRef->A);         ZRef->A=NULL;
      if (ZRef->B)      free(ZRef->B);         ZRef->B=NULL;
// P0 is owned by other packages
//      if (ZRef->P0)     free(ZRef->P0);     ZRef->P0=NULL;
//      if (ZRef->P0LS)   free(ZRef->P0LS);   ZRef->P0LS=NULL;
      if (ZRef->PCube)  free(ZRef->PCube);     ZRef->PCube=NULL;

      ZRef->Version=-1;
      ZRef->LevelNb=0;

      free(ZRef);
   }
   return(1);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_Equal>
 * Creation : Octobre 2011 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Check for vertical reference equality
 *
 * Parametres  :
 *  <ZRef0>     : Vertical reference compare
 *  <ZRef1>     : Vertical reference compare
 *
 * Retour:
 *  <Ok>       : (1=Same, 0=not same)
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
int32_t ZRef_Equal(TZRef *ZRef0,TZRef *ZRef1) {

   if (!ZRef0 || !ZRef1)
      return(0);
      
#ifdef HAVE_VGRID
   if (ZRef0->VGD && ZRef1->VGD)
      return(!Cvgd_vgdcmp((vgrid_descriptor*)ZRef0->VGD,(vgrid_descriptor*)ZRef1->VGD));
#endif
   
   if ((ZRef0->LevelNb!=ZRef1->LevelNb) || (ZRef0->Type!=ZRef1->Type) || (ZRef0->Levels && memcmp(ZRef0->Levels,ZRef1->Levels,ZRef0->LevelNb*sizeof(float))!=0))
      return(0);

   return(1);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_Copy>
 * Creation : Octobre 2011 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Copy a vertical reference into another one
 *
 * Parametres  :
 *  <ZRef>     : Vertical reference to copy
 *
 * Retour:
 *  <ZRef>     : Copied vertical reference (Error=NULL)
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
TZRef* ZRef_Copy(TZRef *ZRef) {

   if (ZRef) ZRef_Incr(ZRef);
   
   return(ZRef);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_HardCCopy>
 * Creation : Octobre 2011 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Copy a vertical reference into another one
 *
 * Parametres  :
 *  <ZRef>     : Vertical reference to copy
 *
 * Retour:
 *  <ZRef>     : Copied vertical reference (Error=NULL)
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
TZRef *ZRef_HardCopy(TZRef *ZRef) {

   TZRef *zref=NULL;
   
   if (ZRef && (zref=ZRef_New())) {
   
      zref->LevelNb=ZRef->LevelNb;
      zref->Levels=(float*)malloc(ZRef->LevelNb*sizeof(float));
      memcpy(zref->Levels,ZRef->Levels,ZRef->LevelNb*sizeof(float));
      
      zref->Type=ZRef->Type;
      zref->PTop=ZRef->PTop;
      zref->PRef=ZRef->PRef;
      zref->ETop=ZRef->ETop;
      zref->RCoef[0]=ZRef->RCoef[0];
      zref->RCoef[1]=ZRef->RCoef[1];
      zref->P0=zref->P0LS=zref->A=zref->B=NULL;
      zref->Version=-1;
      zref->NRef=1;
      zref->Style=ZRef->Style;
   }
   return(zref);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_DecodeRPN>
 * Creation : Novembre 2008 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Decoder les parametres des niveaux a partir du champs HY ou !!
 *
 * Parametres   :
 *  <ZRef>      : Vertical referencer
 *  <Unit>      : FSTD file unit
 *
 * Retour:
 *  <Ok>        : (Index du champs <0=erreur).
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
int32_t ZRef_DecodeRPN(TZRef *ZRef,fst_file* File) {

   fst_record h;
   fst_query *queryt,*queryh,*queryp;
   fst_record record = default_fst_record;
   int32_t        cd,key=0,skip,j,k,kind,ip;
   double    *dbuf=NULL;
   float     *fbuf=NULL;
#ifdef HAVE_VGRID
   char       rfls_S[VGD_LEN_RFLS];
#endif
 
   if (!ZRef) {
      return(0);
   }
      
   if (ZRef->Type==LVL_PRES || ZRef->Type==LVL_UNDEF) {
     return(1);
   }

   memset(&h,0,sizeof(fst_record));
   ZRef->SLEVE=0;

   // Check for toctoc (field !!)   
   record.typvar[0]='X';
   record.nomvar[0]='!';record.nomvar[1]='!';

   queryt = fst24_new_query(File, &record, NULL);
   if (fst24_find_next(queryt, &record)) {
      ZRef->Version=record.ig1;
      ZRef->PRef=10.0;
      ZRef->PTop=record.ig2/10.0;
      ZRef->ETop=0.0;
      ZRef->RCoef[0]=0.0f;
      ZRef->RCoef[1]=0.0f;

      if (!ZRef->A) ZRef->A=(float*)malloc(ZRef->LevelNb*sizeof(float));
      if (!ZRef->B) ZRef->B=(float*)malloc(ZRef->LevelNb*sizeof(float));

#ifdef HAVE_VGRID
      if (Cvgd_new_read((vgrid_descriptor**)&ZRef->VGD,fst24_get_unit(File),-1,-1,-1,-1)==VGD_ERROR) {
         Lib_Log(APP_LIBEER,APP_ERROR,"%s: Unable to initialize vgrid descriptor.\n",__func__);
         return(0);
      }
      // Is it SLEVE type?
      Cvgd_get_char((vgrid_descriptor*)ZRef->VGD,"RFLS",rfls_S,1);
      if (strcmp(rfls_S, VGD_NO_REF_NOMVAR)!=0){
         ZRef->SLEVE=1;
      }
#endif                  
      if (fst24_read_record(&record)) {
         dbuf=(double*)record.data;

         // Read in header info
         switch(ZRef->Version) {
            case 1001: ZRef->Type=LVL_SIGMA;  break;
            case 1002: ZRef->Type=LVL_ETA;    ZRef->PTop=dbuf[record.ni]*0.01; break;
            case 2001: ZRef->Type=LVL_PRES;   break;
            case 1003: ZRef->Type=LVL_ETA;    ZRef->PTop=dbuf[record.ni]*0.01; ZRef->PRef=dbuf[record.ni+1]*0.01; ZRef->RCoef[0]=dbuf[record.ni+2]; break;
            case 5001: ZRef->Type=LVL_HYBRID; ZRef->PTop=dbuf[record.ni]*0.01; ZRef->PRef=dbuf[record.ni+1]*0.01; ZRef->RCoef[0]=dbuf[record.ni+2]; break;
            case 5002:
            case 5003: 
            case 5004: 
            case 5005: 
            case 5100: ZRef->Type=LVL_HYBRID; ZRef->PTop=dbuf[record.ni]*0.01; ZRef->PRef=dbuf[record.ni+1]*0.01; ZRef->RCoef[0]=dbuf[record.ni+2]; ZRef->RCoef[1]=dbuf[record.ni+record.ni]; break;
         }
         skip=dbuf[2];

         // Find corresponding level
         for(k=0;k<ZRef->LevelNb;k++) {
            ip=ZRef_Level2IP(ZRef->Levels[k],ZRef->Type,ZRef->Style);
            for(j=skip;j<record.nj;j++) {
               if (dbuf[j*record.ni]==ip) {
                  ZRef->A[k]=dbuf[j*record.ni+1];
                  ZRef->B[k]=dbuf[j*record.ni+2];
                  break;
               }
            }
            if (j==h.nj) {
               Lib_Log(APP_LIBEER,APP_WARNING,"%s: Could not find level definition for %i.\n",__func__,ip);
            }
         }
      } else {
         Lib_Log(APP_LIBEER,APP_WARNING,"%s: Could not read !! field (fst24_read_record).\n",__func__);
      }
   } else {

      // Check fo regular hybrid (field HY)
      record.typvar[0]='X';
      record.nomvar[0]='H';record.nomvar[1]='Y';
      queryh = fst24_new_query(File, &record, NULL);
      if (fst24_find_next(queryh, &record)) {
         ZRef->PTop=ZRef_IP2Level(record.ip1,&kind);
         ZRef->RCoef[0]=record.ig2/1000.0f;
         ZRef->RCoef[1]=0.0f;
         ZRef->PRef=record.ig1;
//       ZRef->Type=LVL_HYBRID;
         // It might be ETA
         if (ZRef->Type==LVL_SIGMA) {
            ZRef->Type=LVL_ETA;
         }
      } else {

         // Try to figure out if it's SIGMA or ETA
         if (ZRef->Type==LVL_SIGMA || (ZRef->Type==LVL_ETA && ZRef->Version==-1)) {
            /*If we find a PT field, we have ETA coordinate otherwise, its'SIGMA*/
            record.typvar[0]=' ';
            record.nomvar[0]='P';record.nomvar[1]='T';
            queryp = fst24_new_query(File, &record, NULL);
            if (fst24_find_next(queryp, &record)) {
               ZRef->PTop=10.0;
               ZRef->Type=LVL_ETA;
               if (fst24_read_record(&record)) {
                  fbuf=(float*)record.data;
                  ZRef->PTop=fbuf[0];
               } else {
                  Lib_Log(APP_LIBEER,APP_WARNING,"%s: Could not read PT field (fst24_read_record), using default PTOP=10.0.\n",__func__);
               }
            }
            fst24_query_free(queryp);
         }
      }
      fst24_query_free(queryh);

      ZRef->Version=0;
   }
   fst24_query_free(queryt);

   if (ZRef->PCube)  free(ZRef->PCube);  ZRef->PCube=NULL;

   return(key>=0);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_SetRestrictLevels>
 * Creation : Janvier 2012 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Definir une liste de restriction de niveaux
 *
 * Parametres   :
 *  <Levels>    : Liste des niveaux
 *  <NbLevels>  : Nombre de niveaux
 *
 * Retour:
 *  <Ok>        : (Index du champs <0=erreur).
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
int32_t ZRef_SetRestrictLevels(float *Levels,int32_t NbLevels) {

   ZRef_LevelsNb=NbLevels;
   ZRef_Levels=(float*)realloc(ZRef_Levels,ZRef_LevelsNb*sizeof(float));

   memcpy(ZRef_Levels,Levels,ZRef_LevelsNb);
   qsort(ZRef_Levels,ZRef_LevelsNb,sizeof(float),QSort_Float);

   return(ZRef_LevelsNb);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_AddRestrictLevel>
 * Creation : Janvier 2012 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Definir une liste de restriction de niveaux
 *
 * Parametres   :
 *  <Level>     : Niveaux
 *
 * Retour:
 *  <Ok>        : (Index du champs <0=erreur).
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
int32_t ZRef_AddRestrictLevel(float Level) {

   ZRef_LevelsNb++;
   ZRef_Levels=(float*)realloc(ZRef_Levels,ZRef_LevelsNb*sizeof(float));
   ZRef_Levels[ZRef_LevelsNb-1]=Level;

   qsort(ZRef_Levels,ZRef_LevelsNb,sizeof(float),QSort_Float);

   return(ZRef_LevelsNb);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_GetLevels>
 * Creation : Janvier 2012 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Lire la liste des niveaux disponibles
 *
 * Parametres   :
 *  <ZRef>      : Vertical referencer
 *  <H>         : RPN header
 *  <Order>     : Ordre de tri des niveaux (IP1) (-1=decroissant, 0=1 seul niveau, 1=croissant)
 *
 * Retour:
 *  <Ok>        : (Index du champs <0=erreur).
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
int32_t ZRef_GetLevels(TZRef *ZRef,const fst_record* restrict const H,int32_t Order) {

   fst_record h,record;
   fst_query  *query;
   int32_t        l,ip1=0;
   int32_t        k;

   if (Order) {
      /*Get the number of levels*/
      /*In case of # grid, set IP3 to 1 to get NK just for the first tile*/
      memcpy(&h,H,sizeof(fst_record));
      h.ip3=H->grtyp[0]=='#'?1:-1;
      h.ip1=-1;

      query = fst24_new_query(h.file, &h, NULL);
      h.nk=fst24_find_count(query);
      if (!(ZRef->Levels=(float*)malloc(h.nk*sizeof(float)))) {
         return(0);
      }

      // Get the levels
      k=0;
      while(fst24_find_next(query,&record)) {
         ZRef->Levels[k]=ZRef_IP2Level(record.ip1,&l);
         if (k==0) ZRef->Type=l;

         // If a list of restrictive levels is defined, check for validity
         if (h.nk>10 && ZRef_Levels) {
            if (!bsearch(&ZRef->Levels[k],ZRef_Levels,ZRef_LevelsNb,sizeof(float),QSort_Float)) {
               continue;
            }
         }

         // Make sure we use a single type of level, the first we get
         if (l==ZRef->Type) {
            k++;
         }
      }
      ZRef->LevelNb=k;
      fst24_query_free(query);

      /*Sort the levels from ground up and remove duplicates*/
      qsort(ZRef->Levels,ZRef->LevelNb,sizeof(float),Order==-1?QSort_DecFloat:QSort_Float);
      Unique(ZRef->Levels,&ZRef->LevelNb,sizeof(ZRef->Levels[0]));
   } else {
      if (!(ZRef->Levels=(float*)malloc(sizeof(float)))) {
         return(0);
      }
      ZRef->LevelNb=1;
      ZRef->Levels[0]=ZRef_IP2Level(H->ip1,&ZRef->Type);
   }

   ZRef->Style=H->ip1>32768?NEW:OLD;
   if (ZRef->PCube)  free(ZRef->PCube);  ZRef->PCube=NULL;

   return(ZRef->LevelNb);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_K2Pressure>
 * Creation : Octobre 2011 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Convert a level into pressure.
 *
 * Parametres  :
 *  <ZRef>     : Vertical reference to free
 *  <P0>       : Pressure at surface in mb
 *  <P0LS>     : Pressure at surface in mb (Smoothed)
 *  <K>        : Level index to convert
 *
 * Retour:
 *  <Pres>     : Pressure in mb
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
double ZRef_K2Pressure(TZRef* restrict const ZRef,double P0,double P0LS, int32_t K) {

  return(ZRef_Level2Pressure(ZRef,P0,P0LS,ZRef->Levels[K]));
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_KCube2Pressure>
 * Creation : Octobre 2011 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Convert each coordinate in a 3D cube into pressure.
 *
 * Parametres  :
 *  <ZRef>     : Vertical reference to free
 *  <P0>       : Pressure at surface in mb
 *  <P0LS>     : Pressure at surface in mb (Smoothed)
 *  <NIJ>      : 2D dimension of grid
 *  <Log>      : Calculate log of pressure
 *  <Pres>     : Output pressure in mb
 *
 * Retour:
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
int32_t ZRef_KCube2Pressure(TZRef* restrict const ZRef,float *P0,float *P0LS,int32_t NIJ,int32_t Log,float *Pres) {

   int32_t   k,idxk=0,ij;
   int32_t   *ips;
   float pref,ptop,*p0,*p0ls;

   if (!P0 && ZRef->Type!=LVL_PRES) {
      Lib_Log(APP_LIBEER,APP_ERROR,"%s: Surface pressure is required\n",__func__);
      return(0);
   }

   pref=ZRef->PRef;
   ptop=ZRef->PTop;

   switch(ZRef->Type) {
      case LVL_PRES:
         for (k=0;k<ZRef->LevelNb;k++,idxk+=NIJ) {
            for (ij=0;ij<NIJ;ij++) {
               Pres[idxk+ij]=ZRef->Levels[k];
            }
         }
         break;

      case LVL_SIGMA:
         for (k=0;k<ZRef->LevelNb;k++,idxk+=NIJ) {
            for (ij=0;ij<NIJ;ij++) {
               Pres[idxk+ij]=P0[ij]*ZRef->Levels[k];
            }
         }
         break;

      case LVL_ETA:
         for (k=0;k<ZRef->LevelNb;k++,idxk+=NIJ) {
            for (ij=0;ij<NIJ;ij++) {
               Pres[idxk+ij]=ptop+(P0[ij]-ptop)*ZRef->Levels[k];
            }
         }
         break;

      case LVL_UNDEF:
         for (k=0;k<ZRef->LevelNb;k++,idxk+=NIJ) {
            for (ij=0;ij<NIJ;ij++) {
               Pres[idxk+ij]=P0[idxk+ij];
            }
         }
         break;
         
      case LVL_HYBRID:
         if (ZRef->Version<=0) {
            ij=1;
            f77name(hyb_to_pres)(Pres,ZRef->Levels,&ptop,&ZRef->RCoef[0],&pref,&ij,P0,&NIJ,&ij,&ZRef->LevelNb);
         } else {
            ips=malloc(ZRef->LevelNb*sizeof(int));
            for (k=0;k<ZRef->LevelNb;k++) {
               ips[k]=ZRef_Level2IP(ZRef->Levels[k],ZRef->Type,ZRef->Style);
            }
            
            // Really not optimal but Cvgd needs pascals
            p0=(float*)malloc(NIJ*sizeof(float));
            for (ij=0;ij<NIJ;ij++) {
               p0[ij]=P0[ij]*MB2PA;
            }
            
#ifdef HAVE_VGRID	    
            if(ZRef->SLEVE){
               p0ls=(float*)malloc(NIJ*sizeof(float));
               for (ij=0;ij<NIJ;ij++) {
                  p0ls[ij]=P0LS[ij]*MB2PA;
               }
               if (Cvgd_levels_2ref((vgrid_descriptor*)ZRef->VGD,NIJ,1,ZRef->LevelNb,ips,Pres,p0,p0ls,0)) {
                  Lib_Log(APP_LIBEER,APP_ERROR,"%s: Problems in Cvgd_levels_2ref\n",__func__);
                  return(0);
               }	      
            } else {
               if (Cvgd_levels((vgrid_descriptor*)ZRef->VGD,NIJ,1,ZRef->LevelNb,ips,Pres,p0,0)) {
                  Lib_Log(APP_LIBEER,APP_ERROR,"%s: Problems in Cvgd_levels\n",__func__);
                  return(0);
               }    
            }
#else
            Lib_Log(APP_LIBEER,APP_ERROR,"%s: Library not built with VGRID\n",__func__);
#endif
            for (ij=0;ij<NIJ*ZRef->LevelNb;ij++) Pres[ij]*=PA2MB;
            free(ips);
            free(p0);            
         }
         break;
         
      default:
         Lib_Log(APP_LIBEER,APP_ERROR,"%s: Invalid level type (%i)\n",__func__,ZRef->Type);
         return(0);
   }
   
   if (ZRef->POff!=0.0) for (ij=0;ij<NIJ*ZRef->LevelNb;ij++) Pres[ij]+=ZRef->POff;
   if (Log)             for (ij=0;ij<NIJ*ZRef->LevelNb;ij++) Pres[ij]=log(Pres[ij]);
   
   return(1);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_KCube2Meter>
 * Creation : Octobre 2011 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Convert each coordinate in a 3D cube into MASL.
 *
 * Parametres  :
 *  <ZRef>     : Vertical reference to free
 *  <GZ>       : Geopotential height 3D
 *  <NIJ>      : 2D dimension of grid
 *  <Height>   : Output height
 *
 * Retour:
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
int32_t ZRef_KCube2Meter(TZRef* restrict const ZRef,float *GZ,const int32_t NIJ,float *Height) {

   uint32_t k,idxk=0,ij;
   float        topo;

   switch (ZRef->Type) {
      case LVL_PRES:
      case LVL_HYBRID:
      case LVL_SIGMA:
      case LVL_ETA:
      case LVL_UNDEF:
      case LVL_THETA:
         for (k=0;k<NIJ*ZRef->LevelNb;k++) {
            Height[k]=GZ[k]*10.0f;
         }
         break;

      case LVL_MASL:
         for (k=0;k<ZRef->LevelNb;k++,idxk+=NIJ) {
            for (ij=0;ij<NIJ;ij++) {
               Height[idxk+ij]=ZRef->Levels[k];
            }
         }
         break;

      case LVL_MBSL:
         for (k=0;k<ZRef->LevelNb;k++,idxk+=NIJ) {
            for (ij=0;ij<NIJ;ij++) {
               Height[idxk+ij]=-ZRef->Levels[k];
            }
         }
         break;

      case LVL_MAGL:
         /*Add the topography to the gz to get the heigth above the sea*/
         for (k=0;k<ZRef->LevelNb;k++,idxk+=NIJ) {
            for (ij=0;ij<NIJ;ij++) {
               Height[idxk+ij]=(GZ[ij]*10.0f)+ZRef->Levels[k];
            }
         }
         break;

      case LVL_GALCHEN:
         /* Height = GALCHEN * (1 - h0/H) + h0
          * Where
          *  - GALCHEN is the level in gal-chen meters
          *  - h0 is the topography
          *  - H is the Height of the top of the atmosphere
          */
         for (k=0;k<ZRef->LevelNb;k++,idxk+=NIJ) {
            for (ij=0;ij<NIJ;ij++) {
               topo=GZ[ij]*10.0f;
               Height[idxk+ij]=ZRef->Levels[k]*(1.0f-topo/ZRef->PTop)+topo;
            }
         }
         break;

      default:
         Lib_Log(APP_LIBEER,APP_ERROR,"%s: Invalid level type (%i)\n",__func__,ZRef->Type);
         return(0);
   }
   return(1);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_Level2Pressure>
 * Creation : Octobre 2011 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Convert a level into pressure.
 *
 * Parametres  :
 *  <ZRef>     : Vertical reference to free
 *  <P0>       : Pressure at surface in mb
 *  <P0LS>     : Pressure at surface in mb (Smoothed)
 *  <K>        : Level index to convert
 *
 * Retour:
 *  <Pres>     : Pressure in mb
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
double ZRef_Level2Pressure(TZRef* restrict const ZRef,double P0,double P0LS,double Level) {

   double pres=-1.0,pref,ptop,rtop;
   int32_t    ip;

   pref=ZRef->PRef;
   ptop=ZRef->PTop;

   switch(ZRef->Type) {
      case LVL_PRES:
         pres=Level;
         break;

      case LVL_SIGMA:
         pres=P0*Level;
         break;

      case LVL_ETA:
         pres=ptop+(P0-ptop)*Level;
         break;

      case LVL_HYBRID:
         if (ZRef->Version<=0) {
            rtop=ptop/pref;
            pres=pref*Level+(P0-pref)*pow((Level-rtop)/(1.0-rtop),ZRef->RCoef[0]);
         } else {
            ip=ZRef_Level2IP(Level,ZRef->Type,ZRef->Style);
            P0*=MB2PA;
#ifdef HAVE_VGRID	    
	      if(ZRef->SLEVE){
	         P0LS*=MB2PA;
	         if (Cvgd_levels_2ref_8((vgrid_descriptor*)ZRef->VGD,1,1,1,&ip,&pres,&P0,&P0LS,0)) {
		         Lib_Log(APP_LIBEER,APP_ERROR,"%s: Problems in Cvgd_levels_2ref_8\n",__func__);
		         return(0);
	         }	      
	      } else {
	         if (Cvgd_levels_8((vgrid_descriptor*)ZRef->VGD,1,1,1,&ip,&pres,&P0,0)) {
		         Lib_Log(APP_LIBEER,APP_ERROR,"%s: Problems in Cvgd_levels_8\n",__func__);
		         return(0);
	         }    
	      }
#else
            Lib_Log(APP_LIBEER,APP_ERROR,"%s: Library not built with VGRID\n",__func__);
#endif
            pres*=PA2MB;
         }
         break;

      default:
         Lib_Log(APP_LIBEER,APP_ERROR,"%s: Invalid level type (%i)\n",__func__,ZRef->Type);
         return(0);
   }
   
   return(pres+ZRef->POff);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_Pressure2Level>
 * Creation : Octobre 2009 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Calculer la pression a un niveaux specifique
 *
 * Parametres :
 *   <Grid>       : Grille
 *   <Pressure>   : Pression du niveau an mb
 *   <P0>         : Pression au sol en mb
 *
 * Retour:
 *   <Level>      : Niveau du modele
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
double ZRef_Pressure2Level(TZRef* restrict const ZRef,double P0,double Pressure) {

   double level=-1.0;
   double a,b,c,d,r,l,err;
   float  pres,pres0=0,pres1;
   int32_t    z;

   switch(ZRef->Type) {
      case LVL_PRES:
         level=Pressure;
         break;

      case LVL_SIGMA:
         level=Pressure/P0;
         break;

      case LVL_ETA:
         level=(Pressure-ZRef->PTop)/(P0-ZRef->PTop);
         break;

      case LVL_HYBRID:
         //TODO need to add analytical function for GEM4 type
         if (ZRef->Version==5002) {
            // Interpolate in log(p)
            pres=log(Pressure*100.0);
            // Find enclosing levels
            for(z=0;z<ZRef->LevelNb;z++) {
               pres1=pres0;
               pres0=ZRef->A[z]+ZRef->B[z]*log(P0/ZRef->PRef);
               if (pres>=pres0 && pres<=pres1) {
                  level=(pres1-pres)/(pres1-pres0);
                  level=ILIN(ZRef->Levels[z],ZRef->Levels[z-1],level);
                  break;
               }
               if (pres>=pres1 && pres<=pres0) {
                  level=(pres0-pres)/(pres0-pres1);
                  level=ILIN(ZRef->Levels[z-1],ZRef->Levels[z],level);
                  break;
               }
            }
         } else {
            a=ZRef->PRef*100.0;
            b=(P0-ZRef->PRef)*100.0;
            c=(ZRef->PTop/ZRef->PRef)*100.0;
            d=Pressure*100.0;
            r=ZRef->RCoef[0];

            /*Use iterative method Newton-Raphson (developped by Alain Malo)*/
            level=0.5;
            err=1.0;
            while(err>0.0001) {
               // TODO: above 12 mb, this code fails with a nan on the pow(level-c,r-1)
               l=level-((a*level+b*pow((level-c)/(1-c),r)-d)/(a+b*r/(pow(1-c,r))*pow(level-c,r-1)));
               err=fabs(l-level)/level;
               level=l;
            }
         }
         break;

      default:
         Lib_Log(APP_LIBEER,APP_ERROR,"%s: Invalid level type (%i)\n",__func__,ZRef->Type);
   }
   return(level);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_Level2Meter>
 * Creation : Avril 2004 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Derminer le niveaux en metre a d'un niveaux d'un autre type
 *
 * Parametres :
 *  <Type>    : Type de niveau.
 *  <Level>   : Valeur du niveau.
 *
 * Retour:
 *  <niveau>  : Niveau en metres
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
double ZRef_Level2Meter(double Level,int32_t Type) {

   double m=0.0;
   
   // Dans le cas d'un niveau 0 et type SIGMA, on force a ETA
   if (Level==0 && Type==LVL_SIGMA) {
      Type=LVL_ETA;
   }

   switch(Type) {
      case LVL_MASL    : m=Level; break;

      case LVL_MBSL    : m=-Level; break;
      
      case LVL_ETA     : m=ETA2METER(Level); break;
      
      case LVL_SIGMA   : m=SIGMA2METER(Level); break;
      
      case LVL_PRES    : if (Level>226.3203) {
                            m=44330.8*(1-pow(Level/1013.25,0.190263));
                         } else if (Level>54.749) {
                            m=11000.0-6341.624*log(Level/226.3202);
                         } else if (Level>8.68) {
                            m=20000.0+216650.0*(pow(Level/54.749,-0.0292173)-1);
                         } else if (Level>0) {
                            m=32000.0+81660.7*(pow(Level/8.680157,-0.0819469)-1);
                         } else {
                            m=0;
                         }
                         break;
                         
      case LVL_UNDEF   : m=Level; break;
      
      case LVL_MAGL    : m=Level; break;
      
      case LVL_HYBRID  : m=ETA2METER(Level); break;
      
      case LVL_THETA   : m=0; break;
      
      case LVL_GALCHEN : m=Level; break;
      
      case LVL_ANGLE   : m=Level; break;
   }

   return(m);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_IP2Meter>
 * Creation : Octobre 2001 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Transfomer une valeur IP1 en elevation en metres.
 *
 * Parametres :
 *  <IP>      : Valeur IP1
 *
 * Retour:
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
double ZRef_IP2Meter(int32_t IP) {

   int32_t   mode=-1,flag=0,kind=LVL_MASL;
   float level=0.0;
   char  format;

   // Si c'est 0 on force a 0 metre car 0 mb=outter space
   if (IP==0)
      return(0);

   // Convertir en niveau reel
   f77name(convip_plus)(&IP,&level,&kind,&mode,&format,&flag,1);

   return(ZRef_Level2Meter(level,kind));
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_IP2Level>
 * Creation : Mai 2003 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Transfomer un IP1 en niveau.
 *
 * Parametres :
 *  <Level>   : Valeur du niveau
 *  <Type>    : Type de niveau (Coordonnees)
 *
 * Retour:
 *  <IP>      : Niveau
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
double ZRef_IP2Level(int32_t IP,int32_t *Type) {

   int32_t    mode=-1,flag=0;
   float  level=0.0;
   char   format;

   // Convertir en niveau reel
   f77name(convip_plus)(&IP,&level,Type,&mode,&format,&flag,1);

   return(level);
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_Level2IP>
 * Creation : Mai 2003 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Transfomer un niveau en IP1.
 *
 * Parametres :
 *  <Level>   : Valeur du niveau
 *  <Type>    : Type de niveau (Coordonnees)
 *  <Mode>    : Mode (NEW,OLD,DEFAULT)
 *
 * Retour:
 *  <IP>      : IP1
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */

int32_t ZRef_Level2IP(float Level,int32_t Type,TZRef_IP1Mode Mode) {

   int32_t    flag=0,ip=0,mode;
   char   format;

   if (Type<0) {
      return(Level);
   } else {
      mode=Mode==DEFAULT?ZREF_IP1MODE:Mode;
         
      // ETA -> SIGMA
      if (Type==LVL_ETA) {
         Type=LVL_SIGMA;
      }

      // GALCHEN -> MAGL
      if (Type==LVL_GALCHEN) {
         // MAGL can't be encoded in old style so set to MASL in this galchen case
         Type=mode==3?LVL_MASL:LVL_MAGL;
      }

      // Hybrid,MAGL and Teta can't be encoded in old style so force new style
      if (Type==LVL_HYBRID || Type==LVL_MAGL || Type==LVL_THETA) {
         mode=2;
      }

      f77name(convip_plus)(&ip,&Level,&Type,&mode,&format,&flag,1);

      return(ip);
   }
}

/*----------------------------------------------------------------------------
 * Nom      : <ZRef_IPFormat>
 * Creation : Mai 2003 - J.P. Gauthier - CMC/CMOE
 *
 * But      : Formatter un IP en niveau et unite
 *
 * Parametres :
 *  <Buf>     : Chaine a remplir
 *  <IP>      : IP a formatter
 *
 * Retour:
 *  <Type>    : Type de niveau
 *
 * Remarques :
 *
 *----------------------------------------------------------------------------
 */
int32_t ZRef_IPFormat(char *Buf,int32_t IP,int32_t Interval) {

   int32_t   type=-1;
   float lvl;
   
   if (Interval && IP<=32000) {
      sprintf(Buf," %8i %-2s",IP,ZRef_Units[LVL_UNDEF]);
   } else {
      lvl=ZRef_IP2Level(IP,&type);

      switch(type) {
         case LVL_SIGMA : 
         case LVL_THETA : sprintf(Buf," %8.4f %-2s",lvl,ZRef_Units[type]); break;
         case LVL_HYBRID: sprintf(Buf," %8.6f %-2s",lvl,ZRef_Units[type]); break;
         case LVL_UNDEF : sprintf(Buf," %8.3f %-2s",lvl,ZRef_Units[type]); break;
         case LVL_MASL  : 
         case LVL_MBSL  : 
         case LVL_MAGL  : 
         case LVL_PRES  : 
         case LVL_HOUR  : 
         case LVL_MPRES :
         default        : sprintf(Buf," %8.1f %-2s",lvl,ZRef_Units[type]); break;
      }    
   }
   return(type);
}
