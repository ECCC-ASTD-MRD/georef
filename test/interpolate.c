#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <malloc.h>
#include "GeoRef.h"
#include "RPN.h"
#include "fnom.h"

extern int c_fstouv();
extern int c_fstfrm();
extern int c_fstecr();
extern int c_fstinf();
extern int c_fstprm();
extern int c_fstluk();
extern int c_fstsui();
extern int c_fstinl();
extern int c_fstopc();
 
main(int argc,char *argv[]) {

   TGeoRef *gref_in,*gref_out;
   char *in,*out,*grid,*var;
   int id[3],err,key,ni,nj,nk,ig;
   float *p[3];
   TRPNHeader h,hg;

   in=argv[1];   // Input data to interpolate
   out=argv[2];  // Output result of interpolation
   grid=argv[3]; // STD file with a GRID field
   var=argv[4];  // Variable to interpolate

   id[0]=21;
   id[1]=22;
   id[2]=23;
   err=c_fnom(&id[0],in,"STD+RND+R/O",0);
   if ((err=c_fstouv(id[0],"RND"))<0) {
//      App_Log(ERROR,"Problems opening input file %s\n",In);
      return(0);
   }

   err=c_fnom(&id[1],out,"STD+RND+W",0);
   if ((err=c_fstouv(id[1],"RND"))<0) {
//      App_Log(ERROR,"Problems opening input file %s\n",In);
      return(0);
   }

   err=c_fnom(&id[2],grid,"STD+RND+R/O",0);
   if ((err=c_fstouv(id[2],"RND"))<0) {
//      App_Log(ERROR,"Problems opening input file %s\n",In);
      return(0);
   }

   h.IP1=h.IP2=h.IP3=-1;
   h.DEET=h.NPAS=h.DATEO=0,h.DATEV=-1;
   strcpy(hg.NOMVAR,"    ");
   strcpy(hg.TYPVAR,"  ");
   strcpy(hg.ETIKET,"            ");
   strcpy(hg.GRTYP,"  ");

   // Read destination grid
   key=c_fstinf(id[2],&ni,&nj,&nk,-1,"",-1,-1,-1,"","GRID");
   c_fstprm(key,&hg.DATEO,&hg.DEET,&hg.NPAS,&hg.NI,&hg.NJ,&hg.NK,&hg.NBITS,&hg.DATYP,
               &hg.IP1,&hg.IP2,&hg.IP3,hg.TYPVAR,hg.NOMVAR,hg.ETIKET,hg.GRTYP,&hg.IG[X_IG1],
               &hg.IG[X_IG2],&hg.IG[X_IG3],&hg.IG[X_IG4],&hg.SWA,&hg.LNG,&hg.DLTF,&hg.UBC,&hg.EX1,&hg.EX2,&hg.EX3);
   p[2]=(float*)calloc(hg.NI*hg.NJ,sizeof(float));
   gref_out=GeoRef_Create(hg.NI,hg.NJ,hg.GRTYP,hg.IG[X_IG1],hg.IG[X_IG2],hg.IG[X_IG3],hg.IG[X_IG4],id[2]);

   // Look for first corresponding field
   key=c_fstinf(id[0],&ni,&nj,&nk,-1,"",-1,-1,-1,"",var);
   c_fstprm(key,&h.DATEO,&h.DEET,&h.NPAS,&ni,&nj,&nk,&h.NBITS,&h.DATYP,
                  &h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,&h.IG[X_IG1],
                  &h.IG[X_IG2],&h.IG[X_IG3],&h.IG[X_IG4],&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
   p[0]=(float*)calloc(ni*nj*nk,sizeof(float));
   gref_in=GeoRef_Create(ni,nj,h.GRTYP,h.IG[X_IG1],h.IG[X_IG2],h.IG[X_IG3],h.IG[X_IG4],id[0]);

   // Loop on all fields
   while(key>=0) {
      
      strcpy(h.NOMVAR,"    ");
      strcpy(h.TYPVAR,"  ");
      strcpy(h.ETIKET,"            ");
      strcpy(h.GRTYP,"  ");
      c_fstprm(key,&h.DATEO,&h.DEET,&h.NPAS,&ni,&nj,&nk,&h.NBITS,&h.DATYP,
                  &h.IP1,&h.IP2,&h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,&h.IG[X_IG1],
                  &h.IG[X_IG2],&h.IG[X_IG3],&h.IG[X_IG4],&h.SWA,&h.LNG,&h.DLTF,&h.UBC,&h.EX1,&h.EX2,&h.EX3);
      c_fstluk(p[0],key,&ni,&nj,&nk);
      err=GeoRef_Interp(gref_out,gref_in,p[2],p[0]);
      err=c_fstecr(p[2],NULL,-h.NBITS,id[1],h.DATEO,h.DEET,h.NPAS,hg.NI,hg.NJ,hg.NK,h.IP1,h.IP2,h.IP3,h.TYPVAR,h.NOMVAR,h.ETIKET,h.GRTYP,hg.IG[X_IG1],hg.IG[X_IG2],hg.IG[X_IG3],hg.IG[X_IG4],h.DATYP,0);
   
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
