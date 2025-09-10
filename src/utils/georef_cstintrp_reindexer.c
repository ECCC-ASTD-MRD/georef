#include <App.h>
#include <GeoRef.h>
#include <Def.h>
#include <georef_build_info.h>

#define APP_NAME "georef_reindexer"
#define APP_DESC "ECCC/CMC RPN Index conversion tool"

const char *NEMO_Vars[] = { "W001","W002","W003","W004", "I001","I002","I003","I004", "J001","J002","J003","J004", "NAVG","ANG","MASK" };
typedef enum { W001,W002,W003,W004, I001,I002,I003,I004, J001,J002,J003,J004, NAVG,ANG,MASK } NEMO_Idx;

// Usage example call for CANSIPS
// GEM_to NEMO: georef_cstintrp_reindexer -i /home/smco502/datafiles/constants/cmdn/cansips/atm_ocean//Grille_20240202/weights/weights_gem319x262_to_orca1_default_yin.std /home/smco502/datafiles/constants/cmdn/cansips/atm_ocean//Grille_20240202/weights/weights_gem319x262_to_orca1_default_yang.std -o ./atmos-ocean-grids.fstd -g OU -d 319 131 -b 2
// NEMO_to_GEM: georef_cstintrp_reindexer -i /home/smco502/datafiles/constants/cmdn/cansips/atm_ocean//Grille_20240202/weights/weights_orca1_to_gem319x262_default_yin.std /home/smco502/datafiles/constants/cmdn/cansips/atm_ocean//Grille_20240202/weights/weights_orca1_to_gem319x262_default_yang.std -o ./atmos-ocean-grids.fstd -g UO -d 362 292 --orca 1

#define I16_dest_t int32_t // dty: I 16 gets read into 32 bit
#define I1_dest_t  int32_t // dty: I 1  gets read into 32 bit

int ReIndex(char **In,char *Out,char* FromTo,int *OtherDims,int BDW, int Orca) {

   fst_record  rec[2][15];
   fst_record  out=default_fst_record,ang=default_fst_record,crit=default_fst_record;
   fst_file   *fin[2],*fout;
   float      *data, *angle_data;
   int         i=0,j=0,iy,jy,n=0,v=0,w=0,sz=0,idx,g,in,subgrid;
   float       a;
   I16_dest_t *iy_data[2][4];
   I16_dest_t *jy_data[2][4];
   double     *w_data[2][4];
   I1_dest_t  *mask_data[2];
   I16_dest_t *navg_in;
   float      *ang_in;

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
   for(n=W001;n<=MASK;n++){
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

   // Transfer void* data pointers into pointers of the proper type for
   // indexing into the arrays.
   for(n=W001; n<=W004; n++){
       for(subgrid = 0; subgrid <= (In[1] ? 1 : 0); subgrid++){
            w_data[subgrid][n] = rec[subgrid][n].data;
           iy_data[subgrid][n] = rec[subgrid][I001+n].data;
           jy_data[subgrid][n] = rec[subgrid][J001+n].data;
       }
   }

   mask_data[0] = rec[0][MASK].data;
   if(In[1]){
       mask_data[1] = rec[1][MASK].data;
   }

   sz=rec[0][0].ni*rec[0][0].nj;
   data=(float*)malloc(sz*2*15*sizeof(float));
   angle_data=(float*)malloc(sz*2*3*sizeof(float));

   //Yin
   navg_in = rec[0][NAVG].data;
   ang_in = rec[0][ANG].data;
   idx=0;
   for(j=0;j<rec[0][0].nj;j++) {
      for(i=0;i<rec[0][0].ni;i++,idx++) {
         if ((g=navg_in[idx]) > 0) {
            // Check if inside core grid
            in=TRUE;
            if (BDW) {
               for(n=0;n<g;n++){
                  iy=iy_data[0][n][idx]-1;
                  jy=jy_data[0][n][idx]-1;
                  if (iy<BDW || jy<BDW || iy>=(OtherDims[0]-BDW) || jy>=(OtherDims[1]-BDW)) {
                     in=FALSE;
                     break;
                  }
               }
            }
            if (in) {
               data[v++]=i;
               data[v++]=j;
               a=-ang_in[idx]*M_PI/180;
               angle_data[w++]=cos(a);
               angle_data[w++]=sin(a);
               for(n=0;n<g;n++){
                  iy=iy_data[0][n][idx]-1;
                  data[v++]= (Orca && iy==0) ? OtherDims[0]-2 : ((Orca && iy==OtherDims[0]-1) ? 1 : iy);
                  data[v++]= jy_data[0][n][idx]-1;
                  data[v++]= w_data[0][n][idx];
               } 
               data[v++]=angle_data[w++]=REF_INDEX_SEPARATOR;
            }
         }
      }
   }

   // Yang
   if (In[1]) {
      navg_in = rec[1][NAVG].data;
      ang_in = rec[1][ANG].data;
      idx=0;
      for(j=0;j<rec[1][0].nj;j++) {
         for(i=0;i<rec[1][0].ni;i++,idx++) {
            if ((g=navg_in[idx]) > 0) {
               // Check if inside core grid
               in=TRUE;
               if (BDW) {
                  for(n=0;n<g;n++){
                     iy=iy_data[1][n][idx]-1;
                     jy=jy_data[1][n][idx]-1;
                     if (iy<BDW || jy<BDW || iy>=(OtherDims[0]-BDW) || jy>=(OtherDims[1]-BDW)) {
                        in=FALSE;
                        break;
                     }
                  }
               }
               if (in) {
                  data[v++]=i;
                  data[v++]=(FromTo[0]=='U')?j+rec[0][0].nj:j;
                  a=-ang_in[idx]*M_PI/180;
                  angle_data[w++]=cos(a);
                  angle_data[w++]=sin(a);
                  for(n=0;n<g;n++){
                     iy=iy_data[1][n][idx]-1;
                     data[v++]= (Orca && iy==1) ? OtherDims[0]-2 : ((Orca && iy==OtherDims[0]-1) ? 1 : iy);
                     jy=jy_data[1][n][idx]-1;
                     data[v++]= (FromTo[0]=='O') ? jy+OtherDims[1] : jy;
                     data[v++]= w_data[1][n][idx];
                  }
                  data[v++]=angle_data[w++]=REF_INDEX_SEPARATOR;
               } 
            }  
         }
      }
   }   
   data[v++]=angle_data[w++]=REF_INDEX_END;

   out.data=data;
   ang.data=angle_data;
   out.data_type = ang.data_type = FST_TYPE_REAL_IEEE;
   out.data_bits = ang.data_bits = 32;
   out.pack_bits = ang.data_bits = 32;
   out.ni=v;
   ang.ni=w;
   out.nj=ang.nj=1;
   out.nk=ang.nk=1;
   out.dateo= ang.dateo = 0;
   out.deet = ang.deet  = 0;
   out.npas = ang.npas  = 0;
   out.ip1  = ang.ip1   = 0;
   out.ip2  = ang.ip2   = 0;
   out.ip3  = ang.ip3   = IR_WEIGHTINDEX;
   out.ig1  = ang.ig1   = 0;
   out.ig2  = ang.ig2   = 0;
   out.ig3  = ang.ig3   = 0;
   out.ig4  = ang.ig4   = 0;
   if (FromTo) {
      out.typvar[0] = ang.typvar[0] = FromTo[0];
      out.typvar[1] = ang.typvar[1] = FromTo[1];
   } else {
      // Default to OU or OZ
      out.typvar[0] = ang.typvar[0] = 'O';
      out.typvar[1] = ang.typvar[0] = In[1]?'U':'Z';
   }
   strncpy(out.etiket, "GRIDSET", FST_ETIKET_LEN);
   strncpy(ang.etiket, "GRIDSET", FST_ETIKET_LEN);
   strncpy(out.grtyp, "X", FST_GTYP_LEN);
   strncpy(ang.grtyp, "X", FST_GTYP_LEN);
   strncpy(out.nomvar,"####",FST_NOMVAR_LEN);
   strncpy(ang.nomvar,"#@@#",FST_NOMVAR_LEN);
   fst24_write(fout,&out,FST_YES);
   fst24_write(fout,&ang,FST_YES);

   fst24_close(fin[0]);
   In[1] && fst24_close(fin[1]);
   fst24_close(fout);

   return(TRUE);
}

int main(int argc, char *argv[]) {

   int         ok=0,m=-1,code=0,odim[2]={0,0},bdw=0,orca=0;
   char        *in[2]={ NULL, NULL },*out=NULL,*ft=NULL;
 
   TApp_Arg appargs[]=
      { { APP_CHAR,  &in,    2,             "i", "input",  "Input file" },
        { APP_CHAR,  &out,   1,             "o", "output", "Output file" },
        { APP_CHAR,  &ft,    1,             "g", "grid",   "Grid types from->to (ie: U->O = 'OU')" },
        { APP_INT32, &odim,  2,             "d", "dims",   "NI,NJ dimension of the other (sub)grid " },
        { APP_INT32, &bdw,   1,             "b", "border", "Number of border gridpoint not to be used" },
        { APP_INT32, &orca,  1,             "r", "orca",   "Source grid is NEMO ORCA grid, use real i indices when in halo columns" },
        { APP_NIL } };

   App_Init(APP_MASTER,APP_NAME,VERSION,APP_DESC,BUILD_TIMESTAMP);

   if (!App_ParseArgs(appargs,argc,argv,APP_NOARGSFAIL|APP_ARGSLOG)) {
      exit(EXIT_FAILURE);
   }

   // Error checking
   if (!in[0]) {
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
      ok=ReIndex(in,out,ft,odim,bdw,orca);
   }
   code=App_End(ok?0:EXIT_FAILURE);
 
   App_Free();

   exit(code);
}
