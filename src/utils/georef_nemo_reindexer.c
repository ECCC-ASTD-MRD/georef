#include <App.h>
#include <GeoRef.h>
#include <Def.h>
#include <georef_build_info.h>

#define APP_NAME "georef_reindexer"
#define APP_DESC "ECCC/CMC RPN Index conversion tool"

const char *NEMO_Vars[] = { "W001","W002","W003","W004", "I001","I002","I003","I004", "J001","J002","J003","J004", "NAVG","MASK", "ANG"};
typedef enum { W001,W002,W003,W004, I001,I002,I003,I004, J001,J002,J003,J004, NAVG,MASK,ANG } NEMO_Idx;

int ReIndex(char **In,char *Out) {

   fst_record  rec[2][15];
   fst_record  out=default_fst_record,ang=default_fst_record,crit=default_fst_record;
   fst_file   *fin[2],*fout;
   float      *data;
   int         i=0,j=0,n=0,v=0,sz=0,idx,g;

   if (!(fin[0]=fst24_open(In[0],"R/O"))) {
      App_Log(APP_ERROR,"Problems opening input file %s\n",In[0]);
      return(FALSE);
   }

   if (In[1] && !(fin[1]=fst24_open(In[1],"R/O"))) {
      App_Log(APP_ERROR,"Problems opening input file %s\n",In[1]);
      return(FALSE);
   }

   if (!(fout=fst24_open(Out,"R/W"))) {
      App_Log(APP_ERROR,"Problems opening output file %s\n",Out);
      return(FALSE);
   }

   // Create the desination grid
   for(n=0;n<15;n++){
      rec[0][n]=rec[1][n]=default_fst_record;
      strncpy(crit.nomvar,NEMO_Vars[n],FST_NOMVAR_LEN);
      if (!fst24_read(fin[0],&crit,NULL,&rec[0][n])) {
         App_Log(APP_ERROR,"Could not read %s from %s\n",NEMO_Vars[n],In[0]);
         return(FALSE);
      }
      if (In[1] && (!fst24_read(fin[1],&crit,NULL,&rec[1][n]))) {
         App_Log(APP_ERROR,"Could not read %s from %s\n",NEMO_Vars[n],In[1]);
         return(FALSE);
      }
   }
  
   sz=rec[0][0].ni*rec[0][0].nj;
   data=(float*)malloc(sz*2*15*sizeof(float));

   //Yin
   idx=0;
   for(j=0;j<rec[0][0].nj;j++) {
      for(i=0;i<rec[0][0].ni;i++,idx++) {
         if ((g=(((int*)(rec[0][NAVG].data))[idx]))) {
            for(n=0;n<=g-1;n++){
               data[v++]=i;
               data[v++]=j;
               data[v++]=((int*)(rec[0][n+4].data))[idx];
               data[v++]=((int*)(rec[0][n+8].data))[idx];
               data[v++]=((double*)(rec[0][n].data))[idx];
            } 
            data[v++]=REF_INDEX_SEPARATOR;
         }
      }
   }
   // Yang
   if (In[1]) {
      idx=0;
      for(j=0;j<rec[1][0].nj;j++) {
         for(i=0;i<rec[1][0].ni;i++,idx++) {
            if ((g=(((int*)(rec[1][NAVG].data))[idx]))) {
               for(n=0;n<=g-1;n++){
                  data[v++]=i;
                  data[v++]=j+rec[0][0].nj;
                  data[v++]=((int*)(rec[1][n+4].data))[idx];
                  data[v++]=((int*)(rec[1][n+8].data))[idx];
                  data[v++]=((double*)(rec[1][n].data))[idx];
               } 
               data[v++]=REF_INDEX_SEPARATOR;
            }  
         }
      }
   }   
   data[v++]=REF_INDEX_END;

   out.data=data;
   out.data_type = FST_TYPE_REAL_IEEE;
   out.data_bits = 32;
   out.pack_bits = 32;
   out.ni=v;
   out.nj=1;
   out.nk=1;
   out.dateo= 0;
   out.deet = 0;
   out.npas = 0;
   out.ip1  = 0;
   out.ip2  = 0;
   out.ip3  = IR_WEIGHTINDEX;
   out.ig1   = 0;
   out.ig2   = 0;
   out.ig3   = 0;
   out.ig4   = 0;
   out.typvar[0] = 'O';
   out.typvar[1] = In[1]?'U':'Z';
   strncpy(out.etiket, "GRIDSET", FST_ETIKET_LEN);
   strncpy(out.grtyp, "X", FST_GTYP_LEN);
   strncpy(out.nomvar,"####",FST_NOMVAR_LEN);
   fst24_write(fout,&out,FST_SKIP);
 
   for(v=MASK;v<=ANG;v++) {
      //for(i=0;i<sz;i++) {
      //  data[i]=((int*)(rec[0][v].data))[i];
      //}
      memset(data,0x0,2*sz*sizeof(float));
      memcpy(data,rec[0][v].data,sz*rec[0][v].data_bits/8);
      if (In[1]) {
      //  for(i=0;i<sz;i++) {
      //      data[sz+i]=((int*)(rec[1][v].data))[i];
      //  }
         memcpy(&data[sz],rec[1][v].data,sz*rec[1][v].data_bits/8);
      }

      rec[0][v].dateo= 0;
      rec[0][v].deet = 0;
      rec[0][v].npas = 0;
      rec[0][v].ip1  = 0;
      rec[0][v].ip2  = 0;
      rec[0][v].typvar[0] = 'P';
      rec[0][v].typvar[1] = ' ';
      if (In[1]) {
         rec[0][v].grtyp[0] = 'U';
         rec[0][v].nj*=2;
      }
      rec[0][v].data=data;
      fst24_write(fout,&rec[0][v],FST_SKIP);
   }
 
   fst24_close(fin[0]);
   In[1] && fst24_close(fin[0]);
   fst24_close(fout);

   return(TRUE);
}

int main(int argc, char *argv[]) {

   int         ok=0,m=-1,code=0;
   char        *in[2]={ NULL, NULL },*out=NULL;
 
   TApp_Arg appargs[]=
      { { APP_CHAR,  &in,    2,             "i", "input",  "Input file" },
        { APP_CHAR,  &out,   1,             "o", "output", "Output file" },
        { APP_NIL } };

   App_Init(APP_MASTER,APP_NAME,VERSION,APP_DESC,BUILD_TIMESTAMP);

   if (!App_ParseArgs(appargs,argc,argv,APP_NOARGSFAIL|APP_ARGSLOG)) {
      exit(EXIT_FAILURE);
   }

   // Error checking
   if (!in) {
      App_Log(APP_ERROR,"No input standard file specified\n");
      exit(EXIT_FAILURE);
   }
   if (!out) {
      App_Log(APP_ERROR,"No output standard file specified\n");
      exit(EXIT_FAILURE);
   }

   // Launch the app
   App_Start();

   if (code!=EXIT_FAILURE) {
      ok=ReIndex(in,out);
   }
   code=App_End(ok?0:EXIT_FAILURE);
 
   App_Free();

   exit(code);
}
