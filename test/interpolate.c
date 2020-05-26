#include <stdio.h>
#include <strings.h>
#include <malloc.h>

extern int c_fnom();
extern int c_fstouv();
extern int c_fstfrm();
extern int c_fstecr();
extern int c_fstinf();
extern int c_fstprm();
extern int c_fstluk();
extern int c_fstsui();
extern int c_fstinl();
extern int c_fstopc();

typedef struct FSTD_Head {
   int  KEY;              /*Cle du champs*/
   int  DATEO;            /*Date d'origine du champs*/
   int  DATEV;            /*Date de validitee du champs*/
   int  DEET;
   int  NPAS;
   int  NBITS;
   int  DATYP;            /*Type de donnees*/
   int  IP1,IP2,IP3;      /*Specificateur du champs*/
   char TYPVAR[3];        /*Type de variable*/
   char NOMVAR[5];        /*Nom de la variable*/
   char ETIKET[13];       /*Etiquette du champs*/
   char GRTYP[3];
   int  IG1,IG2,IG3,IG4;
   int  SWA;
   int  LNG;
   int  DLTF;
   int  UBC;
   int  EX1,EX2,EX3;
} FSTD_Head;
 
main(int argc,char *argv[]) {

   char *in,*out,*grid,*var;
   int id[3],gid_out,gid_in,err,key,ni,nj,nk;
   float *p[3];
   FSTD_Head h;

   in=argv[1];
   out=argv[2];
   grid=argv[3];
   var=argv[4];

   id[0]=11;
   id[1]=12;
   id[2]=13;
   err=c_fnom(&id[0],in,"STD+RND+R/O");
   if ((err=c_fstouv(id[0],"RND"))<0) {
//      App_Log(ERROR,"Problems opening input file %s\n",In);
      return(0);
   }

   err=c_fnom(&id[1],out,"STD+RND+W");
   if ((err=c_fstouv(id[1],"RND"))<0) {
//      App_Log(ERROR,"Problems opening input file %s\n",In);
      return(0);
   }

   err=c_fnom(&id[2],grid,"STD+RND+R/O");
   if ((err=c_fstouv(id[2],"RND"))<0) {
//      App_Log(ERROR,"Problems opening input file %s\n",In);
      return(0);
   }

   
   h.IP1=h.IP2=h.IP3=-1;
   h.DEET=h.NPAS=h.DATEO=0,h.DATEV=-1;
   strcpy(h.NOMVAR,"    ");
   strcpy(h.TYPVAR,"  ");
   strcpy(h.ETIKET,"            ");
   strcpy(h.GRTYP,"  ");

   key=c_fstinf(id[2],&ni,&nj,&nk,-1,"",-1,-1,-1,"","GRID");
   c_fstprm(key,&h.DATEO,&h.DEET,&h.NPAS,&ni,&nj,&nk,&h.NBITS,&h.DATYP,
               &h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,&h.IG1,
               &h.IG2,&h.IG3,&h.IG4,&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
   p[2]=(float*)calloc(ni*nj*nk,sizeof(float));
   gid_out=c_ezqkdef(ni,nj,h.GRTYP,h.IG1,h.IG2,h.IG3,h.IG4,id[2]);
         
   key=c_fstlir(p[2],id[2],&ni,&nj,&nk,-1,"",-1,-1,-1,"",">>");
   strcpy(h.NOMVAR,"    ");
   strcpy(h.TYPVAR,"  ");
   strcpy(h.ETIKET,"            ");
   strcpy(h.GRTYP,"  ");
   c_fstprm(key,&h.DATEO,&h.DEET,&h.NPAS,&ni,&nj,&nk,&h.NBITS,&h.DATYP,
               &h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,&h.IG1,
               &h.IG2,&h.IG3,&h.IG4,&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
   err=c_fstecr(p[2],NULL,-h.NBITS,id[1],h.DATEO,h.DEET,h.NPAS,ni,nj,nk,h.IP1,h.IP2,h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,h.IG1,h.IG2,h.IG3,h.IG4,h.DATYP,0);
 
   key=c_fstlir(p[2],id[2],&ni,&nj,&nk,-1,"",-1,-1,-1,"","^^");
   strcpy(h.NOMVAR,"    ");
   strcpy(h.TYPVAR,"  ");
   strcpy(h.ETIKET,"            ");
   strcpy(h.GRTYP,"  ");
   c_fstprm(key,&h.DATEO,&h.DEET,&h.NPAS,&ni,&nj,&nk,&h.NBITS,&h.DATYP,
               &h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,&h.IG1,
               &h.IG2,&h.IG3,&h.IG4,&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
   err=c_fstecr(p[2],NULL,-h.NBITS,id[1],h.DATEO,h.DEET,h.NPAS,ni,nj,nk,h.IP1,h.IP2,h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,h.IG1,h.IG2,h.IG3,h.IG4,h.DATYP,0);

   key=c_fstinf(id[0],&ni,&nj,&nk,-1,"",-1,-1,-1,"",var);
   c_fstprm(key,&h.DATEO,&h.DEET,&h.NPAS,&ni,&nj,&nk,&h.NBITS,&h.DATYP,
                  &h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,&h.IG1,
                  &h.IG2,&h.IG3,&h.IG4,&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
   p[0]=(float*)calloc(ni*nj*nk,sizeof(float));
   gid_out=c_ezqkdef(ni,nj,h.GRTYP,h.IG1,h.IG2,h.IG3,h.IG4,id[0]);

   err=c_ezdefset(gid_in,gid_out);
   while(key>=0) {
      
      strcpy(h.NOMVAR,"    ");
      strcpy(h.TYPVAR,"  ");
      strcpy(h.ETIKET,"            ");
      strcpy(h.GRTYP,"  ");
      c_fstprm(key,&h.DATEO,&h.DEET,&h.NPAS,&ni,&nj,&nk,&h.NBITS,&h.DATYP,
                  &h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,&h.IG1,
                  &h.IG2,&h.IG3,&h.IG4,&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
      c_fstluk(p[0],key,&ni,&nj,&nk);
      err=c_ezsint(p[2],p[0]);
      err=c_fstecr(p[2],NULL,-h.NBITS,id[1],h.DATEO,h.DEET,h.NPAS,ni,nj,nk,h.IP1,h.IP2,h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,h.IG1,h.IG2,h.IG3,h.IG4,h.DATYP,0);
   
      key=c_fstsui(id[0],&ni,&nj,&nk);
   }
   
   
   c_fstfrm(id[0]);
   c_fstfrm(id[1]);
   c_fstfrm(id[2]);
   c_fclos(id[0]);
   c_fclos(id[1]);
   c_fclos(id[2]);

   return(0);
}