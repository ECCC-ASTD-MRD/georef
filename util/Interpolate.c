#include <App.h>
#include <rmn.h>
#include "GeoRef.h"
#include "georef_build_info.h"

#define APP_NAME "Interpolate"
#define APP_DESC "ECCC/CMC RPN fstd interpolation tool."

int Interpolate(char *In,char *Out,char *Truth,char *Grid,char **Vars,char *Etiket,TRef_InterpR Type) {

   TGridSet   *gset=NULL;
   TGeoRef    *refin=NULL,*refout=NULL;
   fst_record  crit=default_fst_record,truth,grid=default_fst_record,record=default_fst_record;
   fst_file   *fin,*fout,*fgrid,*ftruth;
   int         n=0;

  if ((fin=fst24_open(In,"R/O"))<0) {
      App_Log(APP_ERROR,"Problems opening input file %s\n",In);
      return(FALSE);
   }

   if ((fout=fst24_open(Out,"R/W"))<0) {
      App_Log(APP_ERROR,"Problems opening output file %s\n",Out);
      return(FALSE);
   }

   if ((fgrid=fst24_open(Grid,"R/O"))<0) {
      App_Log(APP_ERROR,"Problems opening grid file %s\n",Grid);
      return(FALSE);
   }

   // Create the desination grid
   strncpy(crit.nomvar,"GRID",FST_NOMVAR_LEN);
   if (!fst24_read(fgrid,&crit,NULL,&grid)) {
      App_Log(APP_ERROR,"Problems reading grid field\n");
      return(FALSE);    
   }
   refout=GeoRef_Create(grid.ni,grid.nj,grid.grtyp,grid.ig1,grid.ig2,grid.ig3,grid.ig4,(fst_file*)grid.file);

   if (Truth) {
      if ((ftruth=fst24_open(Truth,"R/O"))<0) {
         App_Log(APP_ERROR,"Problems opening truth file %s\n",Truth);
         return(FALSE);
      }

      strncpy(crit.nomvar,"GRID",FST_NOMVAR_LEN);
      if (!fst24_read(ftruth,&crit,NULL,&truth)) {
         App_Log(APP_ERROR,"Problems reading truth field\n");
         return(FALSE);    
      }
   }
  
//   GeoRef_CopyDesc(fout,&grid);
   GeoRef_Write(refout,fout);   

   strncpy(crit.nomvar,Vars[0],FST_NOMVAR_LEN);
   fst_query* query = fst24_new_query(fin,&crit,NULL);
   while(fst24_read_next(query,&record)>0 && n++<10) {
//      App_Log(APP_INFO,"Processing %s %i\n",Vars[0],record->file_index);

      if (!refin) {
         refin=GeoRef_Create(record.ni,record.nj,record.grtyp,record.ig1,record.ig2,record.ig3,record.ig4,(fst_file*)record.file);
      }

      // Proceed with interpolation
      if (!GeoRef_Interp(refout,refin,grid.data,record.data)) {
         App_Log(APP_ERROR,"Interpolation problem");
         return(FALSE);
      }   
    
      // Write results
	   fst24_record_copy_metadata(&grid,&record,FST_META_TIME|FST_META_INFO);
      if (Etiket) strncpy(grid.etiket,Etiket,FST_ETIKET_LEN);
      if (fst24_write(fout, &grid, FST_NO) < 0) {
         App_Log(APP_ERROR, "Unable to write record\n");
         return FALSE;
      }
   }

   
   // Write index
   gset=GeoRef_SetGet(refout,refin);
   if (GeoRef_SetHasIndex(gset)) {
      App_Log(APP_DEBUG,"Saving index containing %i items\n",gset->IndexSize);
      
      if (!GeoRef_SetWrite(gset,fout)){
         return(0);
      }
   }
    
   fst24_close(fin);
   fst24_close(fout);
   fst24_close(fgrid);

   return(1);
}

int main(int argc, char *argv[]) {

   TRef_InterpR interp=IR_LINEAR;
   int          ok=0,code=EXIT_FAILURE;
   char         *etiket=NULL,*in=NULL,*out=NULL,*truth=NULL,*grid=NULL,*type=NULL,*vars[APP_LISTMAX],dtype[]="LINEAR";

   TApp_Arg appargs[]=
      { { APP_CHAR,  &in,    1,             "i", "input",  "Input file" },
        { APP_CHAR,  &out,   1,             "o", "output", "Output file" },
        { APP_CHAR,  &truth, 1,             "t", "truth",  "Truth data file to compare with" },
        { APP_CHAR,  &grid,  1,             "g", "grid",   "Grid file" },
        { APP_CHAR,  &type,  1,             "m", "type",   "Interpolation type (NEAREST,"APP_COLOR_GREEN"LINEAR"APP_COLOR_RESET",CUBIC)" },
        { APP_CHAR,  &etiket,1,             "e", "etiket", "ETIKET for destination field" },
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
   if (!type) {
      type=dtype;
   }

   // Launch the app
   App_Start();
   switch(type[0]) {
      case 'N': interp=IR_NEAREST; break;
      case 'L': interp=IR_LINEAR; break;
      case 'V': interp=IR_CUBIC; break;
      default:
         App_Log(APP_ERROR,"Invalid interpolation method: %s\n",type);
   }

   GeoRef_Options.InterpDegree=IR_LINEAR;
   GeoRef_Options.NoData=0.0;

   ok=Interpolate(in,out,truth,grid,vars,etiket,interp);
   code=App_End(ok?0:EXIT_FAILURE);
   App_Free();

   exit(code);
}
