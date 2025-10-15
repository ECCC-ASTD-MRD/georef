// vim: et:ts=3:sw=3:sts=3

#include <App.h>
#include <GeoRef.h>
#include <Def.h>
#include <georef_build_info.h>


#define APP_NAME "georef_reindexer"
#define APP_DESC "ECCC/CMC RPN Index conversion tool"

//const char *NEMO_Vars[] = { "W001","W002","W003","W004", "I001","I002","I003","I004", "J001","J002","J003","J004", "NAVG","ANG","MASK" };
// typedef enum { W001,W002,W003,W004, I001,I002,I003,I004, J001,J002,J003,J004, NAVG,ANG,MASK } NEMO_Idx;

typedef enum {NAVG, ANG, MASK, NOTHER} Other;
const char *OtherStr[NOTHER] = {"NAVG", "ANG", "MASK"};

typedef enum {W, I, J, N} RecType;
const char *FMT[J+1] = {"W%03d", "I%03d", "J%03d"};

// Usage example call for CANSIPS
// GEM_to NEMO: georef_cstintrp_reindexer -i /home/smco502/datafiles/constants/cmdn/cansips/atm_ocean//Grille_20240202/weights/weights_gem319x262_to_orca1_default_yin.std /home/smco502/datafiles/constants/cmdn/cansips/atm_ocean//Grille_20240202/weights/weights_gem319x262_to_orca1_default_yang.std -o ./atmos-ocean-grids.fstd -g OU -d 319 131 -b 2
// NEMO_to_GEM: georef_cstintrp_reindexer -i /home/smco502/datafiles/constants/cmdn/cansips/atm_ocean//Grille_20240202/weights/weights_orca1_to_gem319x262_default_yin.std /home/smco502/datafiles/constants/cmdn/cansips/atm_ocean//Grille_20240202/weights/weights_orca1_to_gem319x262_default_yang.std -o ./atmos-ocean-grids.fstd -g UO -d 362 292 --orca 1

#define I16_dest_t int32_t // dty: I 16 gets read into 32 bit
#define I1_dest_t  int32_t // dty: I 1  gets read into 32 bit

int ReIndex(char **In,char *Out,char* FromTo,int *OtherDims,int BDW, int Orca, int nb_weights) {

   const int  nsubgrid = (In[1] ? 2 : 1);
   fst_record  rec[nsubgrid][nb_weights][N];
   fst_record others[nsubgrid][3];
   fst_record  out=default_fst_record,ang=default_fst_record,crit=default_fst_record;
   fst_file   *fin[nsubgrid],*fout;
   float      *data_out, *angle_data_out;
   int         i=0,j=0,iy,jy,n=0,v=0,w=0,idx,g,in;
   float       a;
   int        glb_ni[nsubgrid], glb_nj[nsubgrid], sz[nsubgrid];

   // Arrays typed pointers for data of W, I, J input records
   I16_dest_t *iy_data[nsubgrid][nb_weights];
   I16_dest_t *jy_data[nsubgrid][nb_weights];
   double     *w_data[nsubgrid][nb_weights];

   // Arrays of typed pointers for NAVG, MASK, ANG input records
   I1_dest_t  *mask_data[nsubgrid];
   I16_dest_t *navg_data[nsubgrid];
   float      *ang_data[nsubgrid];

   // Other vars
   I16_dest_t *navg_in;
   I1_dest_t  *mask_in;
   float      *ang_in;

   if (!(fout=fst24_open(Out,"R/W"))) {
      App_Log(APP_ERROR,"Problems opening output file %s\n",Out);
      return(FALSE);
   }

   for(int sg = 0; sg < nsubgrid; sg++){
      if (!(fin[sg]=fst24_open(In[sg],"R/O"))) {
         App_Log(APP_ERROR,"Problems opening input file %s\n",In[sg]);
         return(FALSE);
      }

      for(n=0; n<nb_weights; n++){
         for(int t = W; t<=J; t++){
            snprintf(crit.nomvar, FST_NOMVAR_LEN, FMT[t], n+1);
            fprintf(stderr, "sg=%d, n=%d, t=%d, crit.nomvar='%s'\n",sg, n,t,crit.nomvar);
            rec[sg][n][t] = default_fst_record;
            if(!fst24_read(fin[sg], &crit, NULL, &rec[sg][n][t])){
               App_Log(APP_ERROR,"Could not read %s from %s\n", crit.nomvar,In[sg]);
               return(FALSE);
            }
         }
         // Cast by assigning
          w_data[sg][n] = rec[sg][n][W].data;
         iy_data[sg][n] = rec[sg][n][I].data;
         jy_data[sg][n] = rec[sg][n][J].data;
      }

      for(n = NAVG; n <= MASK; n++){
         strncpy(crit.nomvar, OtherStr[n],FST_NOMVAR_LEN);
         others[sg][n] = default_fst_record;
         fprintf(stderr, "sg=%d, n=%d, crit.nomvar='%s'\n", sg, n, crit.nomvar);
         if(!fst24_read(fin[sg], &crit, NULL, &others[sg][n])){
            App_Log(APP_ERROR,"Could not read %s from %s\n", crit.nomvar,In[sg]);
            return(FALSE);
         }
      }

      // Cast by assigning
      mask_data[sg] = others[sg][MASK].data;
      navg_data[sg] = others[sg][NAVG].data;
      ang_data[sg]  = others[sg][ANG].data;

      // Use first record of file for sizes
      glb_ni[sg] = rec[sg][0][0].ni;
      glb_nj[sg] = rec[sg][0][0].nj;
      sz[sg] = glb_ni[sg]*glb_nj[sg];

      fst24_close(fin[sg]);
   }

   int data_per_point = 2 + 3 * nb_weights + 1; // <i,j of point>(2) + <i,j,w for each contributing point, up to nbweights of them>(3*nb_weights) + <separator>(1)
   data_out=(float*)malloc(sz[0]*nsubgrid * data_per_point *sizeof(*data_out) + 1); // + 1 for REF_INDEX_END;
   int angle_data_per_point = 3; // cos, sin, separator
   angle_data_out=(float*)malloc(sz[0]*nsubgrid*angle_data_per_point*sizeof(float) + 1); // +1 for REF_INDEX_END;

   //Yin
   navg_in = navg_data[0];
   ang_in = ang_data[0];
   mask_in = mask_data[0];
   idx=0;
   for(j=0;j<glb_nj[0];j++) {
      for(i=0;i<glb_ni[0];i++,idx++) {
         if ((g=navg_in[idx]) > 0) {
            if(g > nb_weights){
               App_Log(APP_ERROR, "NAVG for point (%d,%d) is greater than nb_weights=%d\n", i, j, nb_weights);
               return FALSE;
            }
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
            // Check if masked point
            if (mask_in[idx] == 0) {
               in=FALSE;
            }
            if (in) {
               data_out[v++]=i;
               data_out[v++]=j;
               a=-ang_in[idx]*M_PI/180;
               angle_data_out[w++]=cos(a);
               angle_data_out[w++]=sin(a);
               for(n=0;n<g;n++){
                  iy=iy_data[0][n][idx]-1;
                  data_out[v++]= (Orca && iy==0) ? OtherDims[0]-2 : ((Orca && iy==OtherDims[0]-1) ? 1 : iy);
                  data_out[v++]= jy_data[0][n][idx]-1;
                  data_out[v++]= w_data[0][n][idx];
               } 
               data_out[v++]=angle_data_out[w++]=REF_INDEX_SEPARATOR;
            }
         }
      }
   }

   // Yang
   if (In[1]) {
      navg_in = navg_data[1];
      ang_in = ang_data[1];
      mask_in = mask_data[1];
      idx=0;
      for(j=0;j<glb_nj[1];j++) {
         for(i=0;i<glb_ni[1];i++,idx++) {
            if ((g=navg_in[idx]) > 0) {
               if(g > nb_weights){
                  App_Log(APP_ERROR, "NAVG for point (%d,%d) is greater than nb_weights=%d\n", i, j, nb_weights);
                  return FALSE;
               }
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
               // Check if masked point
               if (mask_in[idx] == 0) {
                  in=FALSE;
               }
               if (in) {
                  data_out[v++]=i;
                  data_out[v++]=(FromTo[0]=='U')?j+glb_nj[0]:j;
                  a=-ang_in[idx]*M_PI/180;
                  angle_data_out[w++]=cos(a);
                  angle_data_out[w++]=sin(a);
                  for(n=0;n<g;n++){
                     iy=iy_data[1][n][idx]-1;
                     data_out[v++]= (Orca && iy==0) ? OtherDims[0]-2 : ((Orca && iy==OtherDims[0]-1) ? 1 : iy);
                     jy=jy_data[1][n][idx]-1;
                     data_out[v++]= (FromTo[0]=='O') ? jy+OtherDims[1] : jy;
                     data_out[v++]= w_data[1][n][idx];
                  }
                  data_out[v++]=angle_data_out[w++]=REF_INDEX_SEPARATOR;
               }
            }
         }
      }
   }
   data_out[v++]=angle_data_out[w++]=REF_INDEX_END;

   out.data=data_out;
   ang.data=angle_data_out;
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

   fst24_close(fout);

   return(TRUE);
}

int main(int argc, char *argv[]) {

   int         ok=0,code=0,odim[2]={0,0},bdw=0,orca=0, nw=-1;
   char        *in[2]={ NULL, NULL },*out=NULL,*ft=NULL;
 
   TApp_Arg appargs[]=
      { { APP_CHAR,   in,    2,             "i", "input",  "Input file" },
        { APP_INT32, &nw,    1,             "n", "nb-weights", "Number of weights" },
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
   if (nw == -1){
      App_Log(APP_ERROR, "No number of weights specified\n");
      exit(EXIT_FAILURE);
   }

   // Launch the app
   App_Start();

   if (code!=EXIT_FAILURE) {
      ok=ReIndex(in,out,ft,odim,bdw,orca, nw);
   }
   code=App_End(ok?0:EXIT_FAILURE);
 
   App_Free();

   exit(code);
}
