#include <App.h>
#include <GeoRef.h>
#include <Def.h>
#include <georef_build_info.h>

#define APP_NAME "georef_interp"
#define APP_DESC "ECCC/CMC RPN fstd interpolation tool"

extern char *TRef_InterpString[];
extern char *TRef_ExtrapString[];

int WriteResults(fst_file *File,fst_record *In,fst_record *Out,char *Etiket) {

   fst24_record_copy_metadata(Out,In,FST_META_TYPE|FST_META_TIME|FST_META_INFO);
   if (Etiket) strncpy(Out->etiket,Etiket,FST_ETIKET_LEN);
   
   Out->data_type = FST_TYPE_REAL_IEEE;
   Out->data_bits = 32;

   if (fst24_write(File,Out,FST_NO) <= 0) {
      App_Log(APP_ERROR, "Unable to write record\n");
      return FALSE;
   }
   return TRUE;
}

int Interpolate(char *In,char *Out,char *Truth,char *Grid,char **Vars,char *Etiket,char *Winds) {

   TGeoRef    *refin=NULL,*refout=NULL;
   TDef       *defin=NULL,*defout=NULL,*defout1=NULL,*defout2=NULL;
   fst_record  crit=default_fst_record,grid=default_fst_record,record=default_fst_record,truth;
   fst_record  critvv=default_fst_record,gridvv=default_fst_record,recordvv=default_fst_record;
   fst_file   *fin,*fout,*fgrid,*ftruth;
   char       *var;
   int         v=0,n=0,vmode=0;

   if (!(fin=fst24_open(In,"R/O"))) {
      App_Log(APP_ERROR,"Problems opening input file %s\n",In);
      return(FALSE);
   }

   if (!(fout=fst24_open(Out,"R/W"))) {
      App_Log(APP_ERROR,"Problems opening output file %s\n",Out);
      return(FALSE);
   }

   if (!(fgrid=fst24_open(Grid,"R/O"))) {
      App_Log(APP_ERROR,"Problems opening grid file %s\n",Grid);
      return(FALSE);
   }

   // Create the desination grid
   strncpy(crit.nomvar,"GRID",FST_NOMVAR_LEN);
// Test OU  strncpy(crit.etiket,"OCEAN",FST_ETIKET_LEN);
// Test UO   strncpy(crit.etiket,"ATMOS",FST_ETIKET_LEN);
   if (fst24_read(fgrid,&crit,NULL,&grid)) {
      refout=GeoRef_CreateFromRecord(&grid);
   } else {
      // Try to for descriptor only
      fst_record tic_tic = default_fst_record;
      fst_record tac_tac = default_fst_record;

      strncpy(crit.nomvar, ">>", FST_NOMVAR_LEN);
      if (!fst24_find_next(fst24_new_query(fgrid,&crit,NULL),&tic_tic)){
         App_Log(APP_ERROR,"Can't find >> field\n");
      }

      strncpy(crit.nomvar, "^^", FST_NOMVAR_LEN);
      if (!fst24_find_next(fst24_new_query(fgrid,&crit,NULL),&tac_tac)){
         App_Log(APP_ERROR,"Can't find ^^ field\n");
      }
      grid.ni=tic_tic.ni;
      grid.nj=tac_tac.nj;
      grid.nk=1;
      grid.ig1=tic_tic.ip1;
      grid.ig2=tic_tic.ip2;
      grid.ig3=tic_tic.ip3;
      grid.ig4=0;
      grid.grtyp[0]='Z';
      grid.grtyp[1]='\0';
      grid.data=(float*)malloc(grid.ni*grid.nj*grid.nk*sizeof(float));
 
      refout=GeoRef_Create(tic_tic.ni,tac_tac.nj,"Z",tic_tic.ip1,tic_tic.ip2,tic_tic.ip3,0,fgrid);
   }

   if (!refout) {
      App_Log(APP_ERROR,"Problems reading grid field\n");
      return(FALSE);
   }

   // Vector component
   fst24_record_copy_metadata(&gridvv,&grid,FST_META_SIZE|FST_META_GRID);
   gridvv.data=(float*)malloc(grid.ni*grid.nj*grid.nk*sizeof(float));

   if (Truth) {
      if (!(ftruth=fst24_open(Truth,"R/O"))) {
         App_Log(APP_ERROR,"Problems opening truth file %s\n",Truth);
         return(FALSE);
      }

      strncpy(crit.nomvar,"GRID",FST_NOMVAR_LEN);
      if (!fst24_read(ftruth,&crit,NULL,&truth)) {
         App_Log(APP_ERROR,"Problems reading truth field\n");
         return(FALSE);
      }
   }

   GeoRef_WriteFST(refout,NULL,-1,-1,-1,-1,fout);
   defout1=Def_Create(grid.ni,grid.nj,grid.nk,TD_Float32,grid.data,NULL,NULL);
   defout1->NoData=GeoRef_Options.NoData;
   defout2=Def_Create(grid.ni,grid.nj,grid.nk,TD_Float32,grid.data,gridvv.data,NULL);
   defout2->NoData=GeoRef_Options.NoData;

   while((var=Vars[v++]) != NULL) {
      crit=default_fst_record;
      strncpy(crit.nomvar,var,FST_NOMVAR_LEN);
// Test UO      strncpy(crit.etiket,"OCEAN",FST_ETIKET_LEN);
      fst_query* query = fst24_new_query(fin,&crit,NULL);
      n=0;

      while(fst24_read_next(query,&record)>0) {

         if(fst24_record_is_descriptor(&record)){
           continue;
         }

         if (!refin) {
            refin=GeoRef_CreateFromRecord(&record);

            if (GeoRef_Options.Interp==IR_WEIGHTINDEX) {
               TGeoSet * gset = GeoRef_SetGet(refout, refin, NULL);
               if (GeoRef_SetReadFST(gset,IR_WEIGHTINDEX,fgrid) == NULL) {
                  App_Log(APP_ERROR, "%s: Unable to read index in grid file %s\n", __func__, fst24_file_name(fgrid));
                  return(FALSE);
               }
            }
         }

         // Clear output buffer
         Def_Clear(defout1);
         Def_Clear(defout2);
         
         // Winds have to be managed with both components
         vmode=FALSE;
         if (strncmp(record.nomvar,"VV",2)==0) {
            // Skip VV, it will be processed with UU
            continue;
         }
         if (strncmp(record.nomvar,"UU",2)==0) {
            // Read VV component
            fst24_record_copy_metadata(&critvv,&record,FST_META_TIME|FST_META_INFO);
            strncpy(critvv.nomvar,"VV",FST_NOMVAR_LEN);
            if (!fst24_read(fin,&critvv,NULL,&recordvv)) {
               App_Log(APP_ERROR,"Unable to find VV component");
               continue;
            }
            vmode=TRUE;
         }

         // Proceed with interpolation
         defin=Def_Create(record.ni,record.nj,record.nk,TD_Float32,record.data,vmode?recordvv.data:NULL,NULL);
         defout=vmode?defout2:defout1;

         if (!GeoRef_InterpDef(refout, defout, refin, defin, &GeoRef_Options,1)) {
   //     if (!GeoRef_Interp(refout,refin,&GeoRef_Options,grid.data,record.data)) {
            App_Log(APP_ERROR,"Interpolation problem");
            return(FALSE);
         }

         // If winds asked geographic mode
         if (vmode && Winds[1]=='E') {
            GeoRef_UV2UV(refout,(float*)defout->Data[0],(float*)defout->Data[1],(float*)defout->Data[0],(float*)defout->Data[1],NULL,NULL,FSIZE2D(defout));
         }

         // Write results
         WriteResults(fout,&record,&grid,Etiket);
         if (vmode) {
            WriteResults(fout,&recordvv,&gridvv,Etiket);
         }
         Def_Free(defin);

         n++;
      }

      if (!n) {
         App_Log(APP_WARNING,"No record found for var %s'\n",var);
      } else {
         App_Log(APP_INFO,"Processed %i x '%s'\n",n,var);
      }
   }

   // Write index
   TGeoSet *gset=GeoRef_SetGet(refout,refin,&GeoRef_Options);
   if (GeoRef_SetHasIndex(gset)) {
      App_Log(APP_INFO,"Saving index containing %i items\n",gset->IndexSize);

      if (!GeoRef_SetWriteFST(gset,fout)){
         return(0);
      }
   }

   fst24_close(fin);
   fst24_close(fout);
   fst24_close(fgrid);

   return(1);
}

int main(int argc, char *argv[]) {

   int         ok=0,m=-1,code=0;
   char        *etiket=NULL,*in=NULL,*out=NULL,*truth=NULL,*grid=NULL,*method=NULL,*extrap=NULL,*wgeo=NULL,*vars[APP_LISTMAX],*ptr;
   char        dmethod[]="LINEAR",dextrap[]="VALUE",dwgeo[]="GRID";
 
   TApp_Arg appargs[]=
      { { APP_CHAR,  &in,    1,             "i", "input",  "Input file" },
        { APP_CHAR,  &out,   1,             "o", "output", "Output file" },
        { APP_CHAR,  &grid,  1,             "g", "grid",   "Grid file" },
        { APP_CHAR,  &truth, 1,             "t", "truth",  "Truth data file to compare with" },
        { APP_CHAR,  &method,1,             "m", "method", "Interpolation method (NEAREST,"APP_COLOR_GREEN"LINEAR"APP_COLOR_RESET",CUBIC,CONSERVATIVE,NORMALIZED_CONSERVATIVE,MAXIMUM,MINIMUM,SUM,AVERAGE,VARIANCE,SQUARE,NORMALIZED_COUNT,COUNT,VECTOR_AVERAGE,SUBNEAREST,SUBLINEAR,WEIGHTINDEX)" },
        { APP_CHAR,  &extrap,1,             "x", "extrap", "Extrapolation method (MAXIMUM,MINIMUM,"APP_COLOR_GREEN"[VALUE]"APP_COLOR_RESET",ABORT)" },
        { APP_CHAR,  &etiket,1,             "e", "etiket", "ETIKET for destination field" },
        { APP_CHAR,  &wgeo,1,               "w", "winds",  "Wind output type ("APP_COLOR_GREEN"GRID"APP_COLOR_RESET",GEO)" },
        { APP_CHAR,  vars,  APP_LISTMAX-1,  "n", "nomvar", "List of variable to process" },
        { APP_NIL } };

   memset(vars,0x0,APP_LISTMAX*sizeof(vars[0]));
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
   if (!grid) {
      App_Log(APP_ERROR,"No grid standard file specified\n");
      exit(EXIT_FAILURE);
   }

   if (!vars[0]) {
      App_Log(APP_ERROR,"No variable specified\n");
      exit(EXIT_FAILURE);
   }

   if (!method) {
      method=dmethod;
   }
   if (!extrap) {
      extrap=dextrap;
   }
   if (!wgeo) {
      wgeo=dwgeo;
   }

   // Launch the app
   App_Start();

   // Validate options
   m=0;
   while (TRef_InterpString[m]) {
      if (strcmp(TRef_InterpString[m],method)==0) {
         GeoRef_Options.Interp=m;
         break;
      }
      m++;
   }
   if (m>IR_WEIGHTINDEX) {
      App_Log(APP_ERROR,"Invalid interpolation method: %s\n",method);
      code=EXIT_FAILURE;
   }

   m=0;
   while (TRef_ExtrapString[m]) {
      if (strcmp(TRef_ExtrapString[m],extrap)==0) {
         GeoRef_Options.Extrap=m;
         break;
      }
      m++;
   }
   if (m>ER_ABORT) {
      GeoRef_Options.NoData=strtof(extrap,&ptr);

      if (ptr==extrap) {
         App_Log(APP_ERROR,"Invalid extrapolation: %s\n",extrap);
         code=EXIT_FAILURE;
      }
   }

   if (code!=EXIT_FAILURE) {
      ok=Interpolate(in,out,truth,grid,vars,etiket,wgeo);
   }
   code=App_End(ok?0:EXIT_FAILURE);
   App_Free();

   exit(code);
}
