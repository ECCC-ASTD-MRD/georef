#include <App.h>
#include "Def.h"
#include "GeoRef.h"
#include "georef_build_info.h"
#include "RPN.h"

#define APP_NAME "ReGrid"
#define APP_DESC "SMC/CMC/EERS RPN fstd re-gridding tool."


int ReGrid(char *In,char *Out,char *Grid,char **Vars) {

   TGridSet   *gset;
   TGeoRef    *refin,*refout;

   fst_record  in,grid,idx,crit=default_fst_record;
   fst_file   *fin,fout,fgrid;

  if ((fin=fst24_open(In,"STD+RND+R/O"))<0) {
      App_Log(APP_ERROR,"Problems opening input file %s\n",In);
      return(0);
   }

   if ((fout=fst24_open(Out,"STD+RND+R/W"))<0) {
      App_Log(APP_ERROR,"Problems opening output file %s\n",Out);
      return(0);
   }

   if ((fgrid=fst24_open(Grid,"STD+RND+R/O"))<0) {
      App_Log(APP_ERROR,"Problems opening grid file %s\n",Grid);
      return(0);
   }

   // Create the desination grid
   strncpy(crit,"P0",FST_NOMVAR_LEN);
   if (!fst24_read(fgrid,&crit,NULL,&recout))) {
      App_Log(APP_ERROR,"Problems reading grid field\n");
      return(0);    
   }
   refout=GeoRef_Create(record.ni,record.nj,record.grtyp,record.ig1,record.ig2,record.ig3,record.ig4,record.file);

   // Create index field
   if (!(idx=RPN_FieldNew(in->Def->NIJ*100,1,1,1,TD_Float32))) {
      return(0);       
   }

   grid->Def->NoData=0.0;
 
   strncpy(crit,Vars[0],FST_NOMVAR_LEN);
   fst_query* query = fst24_new_query(fin,&crit,NULL);
   while(fst24_read_next(query,&record)>0) {
      App_Log(APP_INFO,"Processing %s %i\n",Vars[0],record->file_index);

      if (!(def=Def_New(record.ni,record.nj,record.nk,1,TD_Float32))) {
         Lib_Log(APP_LIBEER,APP_ERROR,"%s: Could not allocate memory for fld\n",__func__);
         return(NULL);
      }
      if (!refin) {
         refin=GeoRef_Create(record.ni,record.nj,record.grtyp,record.ig1,record.ig2,record.ig3,record.ig4,record.file);
      }


      // Reset result grid
      Def_Clear(grid->Def);
     
      // Proceed with interpolation
      if (!(idx->Head.NI=Def_GridInterpConservative(refout,refin,grid->Def,in->Def,TRUE))) {
         return(0);    
      }
      
      // Write results
      RPN_CopyHead(&grid->Head,&in->Head);
      recout
      fst24_write(fout,recout);
   }
   
   // Write index
   gset=GeoRef_SetGet(grid->GRef,in->GRef);
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

   int      ok=0,code=EXIT_FAILURE;
   char     *in=NULL,*out=NULL,*grid=NULL,*type=NULL,*vars[APP_LISTMAX];

   TApp_Arg appargs[]=
      { { APP_CHAR,  &in,   1,             "i", "input",  "Input file" },
        { APP_CHAR,  &out,  1,             "o", "output", "Output file" },
        { APP_CHAR,  &grid, 1,             "g", "grid",   "Grid file" },
        { APP_CHAR,  &type, 1,             "t", "type",   "Interpolation type ("APP_COLOR_GREEN"CONSERVATIVE"APP_COLOR_RESET",NORMALIZED_CONSERVATIVE)" },
        { APP_CHAR,  vars,  APP_LISTMAX-1, "n", "nomvar", "List of variable to process" },
        { APP_NIL } };

   memset(vars,0x0,APP_LISTMAX*sizeof(vars[0]));
   App_Init(APP_MASTER,APP_NAME,VERSION,APP_DESC,BUILD_TIMESTAMP);

   if (!App_ParseArgs(appargs,argc,argv,APP_NOARGSFAIL|APP_ARGSLOG)) {
      exit(EXIT_FAILURE);      
   }
   
   /*Error checking*/
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

   /*Launch the app*/
   App_Start();
   ok=ReGrid(in,out,grid,vars);
   code=App_End(ok?-1:EXIT_FAILURE);
   App_Free();

   exit(code);
}
