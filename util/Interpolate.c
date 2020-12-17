/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Librairie de fonctions utiles
 * Creation     : Janvier 2015
 * Auteur       : Jean-Philippe Gauthier
 *
 * Description: RPN fstd interpolation test
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

#include "App.h"
#include "Def.h"
#include "GeoRef.h"
#include "RPN.h"

#define APP_NAME "Interpolate"
#define APP_DESC "ECCC/CMC RPN fstd interpolation tool."

int Interpolate(char *In,char *Out,char *Grid,char **Vars,TDef_InterpR Type) {

#ifdef HAVE_RMN
   TGridSet  *gset=NULL;
   TRPNField *in,*grid,*idx;
   int  fin,fout,fgrid,n=0;

   if ((fin=cs_fstouv(In,"STD+RND+R/O"))<0) {
      App_Log(ERROR,"Problems opening input file %s\n",In);
      return(0);
   }

   if ((fout=cs_fstouv(Out,"STD+RND+R/W"))<0) {
      App_Log(ERROR,"Problems opening output file %s\n",Out);
      return(0);
   }

   if ((fgrid=cs_fstouv(Grid,"STD+RND+R/O"))<0) {
      App_Log(ERROR,"Problems opening grid file %s\n",Grid);
      return(0);
   }
   if (!(grid=RPN_FieldRead(fgrid,-1,"",-1,-1,-1,"","GRID"))) {
      App_Log(ERROR,"Problems reading grid field\n");
      return(0);    
   }
   //RPN_FieldReadGrid(grid);

   if (!(in=RPN_FieldRead(fin,-1,"",-1,-1,-1,"",Vars[0]))) {
      App_Log(ERROR,"Problems reading input field\n");
      return(0);  
   }
   //RPN_FieldReadGrid(in);

   // Create index field
   if (!(idx=RPN_FieldNew(in->Def->NIJ*100,1,1,1,TD_Float64))) {
      return(0);       
   }

   RPN_CopyDesc(fout,&grid->Head);
  
   memcpy(&idx->Head,&grid->Head,sizeof(TRPNHeader));
   strcpy(idx->Head.NOMVAR,"#%");
   idx->Head.NJ=1;
   idx->Head.GRTYP[0]='X';
   idx->Head.NBITS=32;
   idx->Head.DATYP=5;

   grid->Def->NoData=0.0;
 
   //fprintf(stderr,"...%c...%c...\n",grid->GRef->RPNHead.GRREF[0],in->GRef->RPNHead.GRREF[0]);
   while(in->Head.KEY>0 && n++<10) {
      App_Log(INFO,"Processing %s %i\n",Vars[0],in->Head.KEY);

      // Reset result grid
      Def_Clear(grid->Def);
     
      // Proceed with interpolation
      if (!(Def_GridInterp(grid->GRef,grid->Def,in->GRef,in->Def,Type,ER_VALUE,&gset))) {
         App_Log(ERROR,"%s: EZSCINT interpolation problem",__func__);
         return(0);    
      }
    
      // Write results
      RPN_CopyHead(&grid->Head,&in->Head);
      RPN_FieldWrite(fout,grid);
	
      if ((in->Head.KEY=cs_fstsui(fin,&in->Head.NI,&in->Head.NJ,&in->Head.NK))>0) {     
         if (!RPN_FieldReadIndex(fin,in->Head.KEY,in)) {
            return(0);
         }
      }
   }
   
   // Write index
   if (gset) {
      App_Log(DEBUG,"Saving index containing %i items\n",gset->IndexSize);
      
      if (!GeoRef_SetWrite(fout,gset)) {
         return(0);
      }
   }
   
   cs_fstfrm(fin);
   cs_fstfrm(fout);
   cs_fstfrm(fgrid);
#endif

   return(1);
}

int main(int argc, char *argv[]) {

   TDef_InterpR interp=IR_LINEAR;
   int          ok=0,code=EXIT_FAILURE;
   char        *in=NULL,*out=NULL,*grid=NULL,*type=NULL,*vars[APP_LISTMAX],dtype[]="LINEAR";

   TApp_Arg appargs[]=
      { { APP_CHAR,  &in,   1,             "i", "input",  "Input file" },
        { APP_CHAR,  &out,  1,             "o", "output", "Output file" },
        { APP_CHAR,  &grid, 1,             "g", "grid",   "Grid file" },
        { APP_CHAR,  &type, 1,             "t", "type",   "Interpolation type (NEAREST,"APP_COLOR_GREEN"LINEAR"APP_COLOR_RESET",CUBIC)" },
        { APP_CHAR,  vars,  APP_LISTMAX-1, "n", "nomvar", "List of variable to process" },
        { APP_NIL } };

   memset(vars,0x0,APP_LISTMAX*sizeof(vars[0]));
   App_Init(APP_MASTER,APP_NAME,VERSION,APP_DESC,BUILD_TIMESTAMP);

   if (!App_ParseArgs(appargs,argc,argv,APP_NOARGSFAIL|APP_ARGSLOG)) {
      exit(EXIT_FAILURE);      
   }
   
   /*Error checking*/
   if (!in) {
      App_Log(ERROR,"No input standard file specified\n");
      exit(EXIT_FAILURE);
   }
   if (!out) {
      App_Log(ERROR,"No output standard file specified\n");
      exit(EXIT_FAILURE);
   }
   if (!grid) {
      App_Log(ERROR,"No grid standard file specified\n");
      exit(EXIT_FAILURE);
   }

   if (!vars[0]) {
      App_Log(ERROR,"No variable specified\n");
      exit(EXIT_FAILURE);
   }
   if (!type) {
      type=dtype;
   }

   /*Launch the app*/
   App_Start();
   switch(type[0]) {
      case 'N': interp=IR_NEAREST; break;
      case 'L': interp=IR_LINEAR; break;
      case 'V': interp=IR_CUBIC; break;
      default:
         App_Log(ERROR,"Invalid interpolation method: %s\n",type);
   }

   ok=Interpolate(in,out,grid,vars,interp);
   code=App_End(ok?-1:EXIT_FAILURE);
   App_Free();

   exit(code);
}
