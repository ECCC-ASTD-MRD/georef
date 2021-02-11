/*=========================================================
 * Environnement Canada
 * Centre Meteorologique Canadien
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Lecture et traitements de fichiers raster
 * Fichier      : GeoRef_RPN.c
 * Creation     : Mars 2005 - J.P. Gauthier
 *
 * Description  : Fonctions de manipulations de projections aux standard RPN.
 *
 * Remarques    :
 *
 * License      :
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
 *=========================================================
 */

#include "App.h"
#include "GeoRef.h"
#include "Def.h"
#include "RPN.h"
#include "Vertex.h"

int      GeoRef_RPNProject(TGeoRef *Ref,double X,double Y,double *Lat,double *Lon,int Extrap,int Transform);
int      GeoRef_RPNUnProject(TGeoRef *Ref,double *X,double *Y,double Lat,double Lon,int Extrap,int Transform);

int GeoRef_RPNDefXG(TGeoRef* Ref) {

   switch (Ref->GRTYP[0]) {
      case 'A':
      case 'G':
         Ref->RPNHead.XG[X_DLON]  = 360. /Ref->NX;
         Ref->RPNHead.XG[X_SWLON] = 0.0;
         switch (Ref->RPNHead.IG[X_IG1]) {
				case 0:
				   Ref->RPNHead.XG[X_DLAT] = 180./Ref->NY;
				   Ref->RPNHead.XG[X_SWLAT] = -90. + 0.5*Ref->RPNHead.XG[X_DLAT];
				   break;

				case 1:
				   Ref->RPNHead.XG[X_DLAT] = 90./Ref->NY;
				   Ref->RPNHead.XG[X_SWLAT] = 0.5*Ref->RPNHead.XG[X_DLAT];
				   Ref->Type |= GRID_EXPAND;
				   break;

				case 2:
				   Ref->RPNHead.XG[X_DLAT] = 90./Ref->NY;
				   Ref->RPNHead.XG[X_SWLAT] = -90. + 0.5*Ref->RPNHead.XG[X_DLAT];
				   Ref->Type |= GRID_EXPAND;
				   break;

				default:
			      App_Log(ERROR,"%s: 'A' grid has to be Global/North/South\n",__func__);
               return(-1);
				   break;
			}

         switch(Ref->RPNHead.IG[X_IG2]) {
	         case 1:
	            Ref->Type|=GRID_YINVERT;
	            break;

	         default:
	            break;
	      }
         break;

      case 'B':
         Ref->RPNHead.XG[X_DLON] = 360. /(Ref->NX-1);
         Ref->RPNHead.XG[X_SWLON] = 0.0;
         switch (Ref->RPNHead.IG[X_IG1]) {
	         case 0:
	            Ref->RPNHead.XG[X_DLAT] = 180./(Ref->NY-1);
	            Ref->RPNHead.XG[X_SWLAT] = -90.;
	            break;

	         case 1:
	            Ref->RPNHead.XG[X_DLAT] = 90./(Ref->NY-1);
	            Ref->RPNHead.XG[X_SWLAT] = 0.;
	            Ref->Type |= GRID_EXPAND;
	            break;

	         case 2:
	            Ref->RPNHead.XG[X_DLAT] = 90./(Ref->NY-1);
	            Ref->RPNHead.XG[X_SWLAT] = -90.;
	            Ref->Type |= GRID_EXPAND;
	            break;

	         default:
  			      App_Log(ERROR,"%s: 'B' grid has to be Global/North/South\n",__func__);
	            return(-1);
	      }

         switch(Ref->RPNHead.IG[X_IG2]) {
	         case 1:
	            Ref->Type|=GRID_YINVERT;
	            break;

	         default:
	            break;
	      }
         break;

      case 'E':
         f77name(cigaxg)(Ref->GRTYP,&Ref->RPNHead.XG[X_LAT1],&Ref->RPNHead.XG[X_LON1],&Ref->RPNHead.XG[X_LAT2],&Ref->RPNHead.XG[X_LON2],&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);
      /*      Ref->RPNHead.XG[X_DLAT] = 180./Ref->NY;
	      Ref->RPNHead.XG[X_DLON] = 360./(Ref->NX-1);
	      Ref->RPNHead.XG[X_SWLON] = 0.0;
	      Ref->RPNHead.XG[X_SWLAT] = -90. + 0.5*Ref->RPNHead.XG[X_DLAT];
      */
         break;

      case 'H':
      case 'Y':
      case '!':
         break;

      case '#':
      case 'Z':
         if (Ref->RPNHead.GRREF[0] == 'N') Ref->Hemi = NORTH;
         if (Ref->RPNHead.GRREF[0] == 'S') Ref->Hemi = SOUTH;
         if (Ref->RPNHead.GRREF[0] == 'E' || Ref->RPNHead.GRREF[0]== 'L') {
            f77name(cigaxg)(Ref->RPNHead.GRREF,&Ref->RPNHead.XGREF[X_LAT1], &Ref->RPNHead.XGREF[X_LON1], &Ref->RPNHead.XGREF[X_LAT2], &Ref->RPNHead.XGREF[X_LON2],&Ref->RPNHead.IGREF[X_IG1], &Ref->RPNHead.IGREF[X_IG2], &Ref->RPNHead.IGREF[X_IG3], &Ref->RPNHead.IGREF[X_IG4]);
         }
         break;

      case 'L':
         f77name(cigaxg)(Ref->GRTYP,&Ref->RPNHead.XG[X_SWLAT], &Ref->RPNHead.XG[X_SWLON], &Ref->RPNHead.XG[X_DLAT], &Ref->RPNHead.XG[X_DLON],&Ref->RPNHead.IG[X_IG1], &Ref->RPNHead.IG[X_IG2], &Ref->RPNHead.IG[X_IG3], &Ref->RPNHead.IG[X_IG4]);
         break;

      case 'N':
         f77name(cigaxg)(Ref->GRTYP,&Ref->RPNHead.XG[X_PI], &Ref->RPNHead.XG[X_PJ], &Ref->RPNHead.XG[X_D60], &Ref->RPNHead.XG[X_DGRW],&Ref->RPNHead.IG[X_IG1], &Ref->RPNHead.IG[X_IG2], &Ref->RPNHead.IG[X_IG3], &Ref->RPNHead.IG[X_IG4]);
         Ref->Hemi = NORTH;
         break;

      case 'S':
         f77name(cigaxg)(Ref->GRTYP,&Ref->RPNHead.XG[X_PI], &Ref->RPNHead.XG[X_PJ], &Ref->RPNHead.XG[X_D60], &Ref->RPNHead.XG[X_DGRW],&Ref->RPNHead.IG[X_IG1], &Ref->RPNHead.IG[X_IG2], &Ref->RPNHead.IG[X_IG3], &Ref->RPNHead.IG[X_IG4]);
         Ref->Hemi = SOUTH;
         break;

      case 'T':
		   //TODO: What's T
         f77name(cigaxg)(Ref->GRTYP,&Ref->RPNHead.XG[X_TD60], &Ref->RPNHead.XG[X_TDGRW], &Ref->RPNHead.XG[X_CLAT], &Ref->RPNHead.XG[X_CLON],&Ref->RPNHead.IG[X_IG1], &Ref->RPNHead.IG[X_IG2], &Ref->RPNHead.IG[X_IG3], &Ref->RPNHead.IG[X_IG4]);
         break;

      default:
	      App_Log(DEBUG,"%s: Grid type not supported %c\n",__func__,Ref->GRTYP[0]);
         return(-1);
    }

   return(0);
}

//! Insert a grid entry into the list of grids managed by ezscint.  Can be used
//! with regular and irregular ('Y', 'Z') grids, although it is not very useful
//! for regular grids.
//! @param ni Horizontal size of the grid
//! @param nj
//! @param grtyp Grid type ('A', 'B', 'E', 'G', 'L', 'N', 'S','Y', 'Z', '#', '!')
//! @param grref Reference grid type ('E', 'G', 'L', 'N', 'S')
//! @param ig1 ig1 value associated to the reference grid
//! @param ig2 ig2 value associated to the reference grid
//! @param ig3 ig3 value associated to the reference grid
//! @param ig4 ig4 value associated to the reference grid
//! @param ax Positional axis mapped to the '>>' record
//! @param ay Positional axis mapped to the '^^' record
//!
//! If the grid type corresponds to a regular grid type (eg. 'A', 'G', 'N', etc.),
//! then the parameters IG1 through IG4 are taken from an ordinary data record
//! and grref, ax and ay are not used.
//!
//! If grtyp == 'Z' or '#', the dimensions of ax=ni and ay=nj.
//! If grtyp == 'Y', the dimensions of ax=ay=ni*nj. 

TGeoRef* GeoRef_Set(TGeoRef *Ref,int NI,int NJ,char* GRTYP,char* GRREF,int IG1,int IG2,int IG3,int IG4,double* AX,double* AY) {
   
   Ref->RPNHead.GRTYP[0]=Ref->GRTYP[0] = GRTYP[0];
   Ref->RPNHead.GRTYP[1]=Ref->GRTYP[1] = '\0';
   Ref->RPNHead.GRREF[0] = GRREF?GRREF[0]:'\0';
   Ref->RPNHead.GRREF[1] = '\0';
   Ref->RPNHead.NI=Ref->NX = NI;
   Ref->RPNHead.NJ=Ref->NY = NJ;
   Ref->RPNHead.IG[X_IG1] = IG1;
   Ref->RPNHead.IG[X_IG2] = IG2;
   Ref->RPNHead.IG[X_IG3] = IG3;
   Ref->RPNHead.IG[X_IG4] = IG4;
   Ref->i1 = 1;
   Ref->i2 = NI;
   Ref->j1 = 1;
   Ref->j2 = NJ;

   switch (GRTYP[0]) {
      case 'Z':
         f77name(cigaxg)(Ref->RPNHead.GRREF,&Ref->RPNHead.XGREF[X_LAT1],&Ref->RPNHead.XGREF[X_LON1],&Ref->RPNHead.XGREF[X_LAT2],&Ref->RPNHead.XGREF[X_LON2],&Ref->RPNHead.IGREF[X_IG1],&Ref->RPNHead.IGREF[X_IG2],&Ref->RPNHead.IGREF[X_IG3],&Ref->RPNHead.IGREF[X_IG4]);
      case '#':
      case 'Y':
         Ref->AX = AX;
         Ref->AY = AY;
         break;
   }

   GeoRef_Size(Ref,0,0,NI-1,NJ-1,0);

   return(Ref);
}

TGeoRef* GeoRef_CreateInMemory(int NI,int NJ,char* GRTYP,char* GRREF,int IG1,int IG2,int IG3,int IG4,double* AX,double* AY) {
   
   TGeoRef* ref,*fref;

   ref = GeoRef_New();

   GeoRef_Set(ref,NI,NJ,GRTYP,GRREF,IG1,IG2,IG3,IG4,AX,AY);
 
   // TODO: Would be more efficient to find without creating one
   if (fref = GeoRef_Find(ref)) {
      // This georef already exists
      free(ref);
      GeoRef_Incr(fref);
      return(fref);
   }

   // This is a new georef
   GeoRef_Add(ref);

   GeoRef_RPNDefXG(ref);
   GeoRef_AxisDefine(ref,AX,AY);
   GeoRef_Qualify(ref);

   return(ref);
}

TGeoRef* GeoRef_CreateU(int NI,int NJ,char *GRTYP,char *GRREF,int VerCode,int NbSub,TGeoRef **Subs) {

   int  i;
   TGeoRef *ref,*fref,*sub_gd;
    
   if (NbSub <= 1) {
      App_Log(ERROR,"%s: NbSub given is less than 2\n",__func__);
      return(NULL);
   }
   if (VerCode != 1) {
      App_Log(ERROR,"%s: Invalid VerCode\n",__func__);
      return(NULL);
   }

   ref = GeoRef_New();
  
   if (VerCode == 1) {
      sub_gd = Subs[0];

      ref->RPNHead.GRTYP[0]=ref->GRTYP[0] = GRTYP[0];
      ref->RPNHead.GRREF[0] = GRREF[0];
      ref->NX       = NI;
      ref->NY       = NJ;
         
      // To add more uniqueness to the super-grid index Yin-Yang grid, we also add the rotation of YIN
      ref->RPNHead.IG[X_IG1]  = sub_gd->RPNHead.IG[X_IG1];
      ref->RPNHead.IG[X_IG2]  = sub_gd->RPNHead.IG[X_IG2];
      ref->RPNHead.IG[X_IG3]  = sub_gd->RPNHead.IG[X_IG3];
      ref->RPNHead.IG[X_IG4]  = sub_gd->RPNHead.IG[X_IG4];
      ref->RPNHead.IGREF[X_IG1]=VerCode;
      ref->RPNHead.IGREF[X_IG2]=0;
      ref->RPNHead.IGREF[X_IG3]=0;
      ref->RPNHead.IGREF[X_IG4]=0;
      ref->NbSub= NbSub;
   }
  
   // This georef already exists
   if (fref=GeoRef_Find(ref)) {
      free(ref);
      GeoRef_Incr(fref);
      return(fref);
   }

   // This is a new georef
   GeoRef_Add(ref);

   ref->Subs = (TGeoRef **)malloc(NbSub*sizeof(TGeoRef*));

   for (i=0; i < NbSub; i++) {
      ref->Subs[i] = Subs[i];
      GeoRef_MaskYYDefine(Subs[i]);
      App_Log(DEBUG,"%s: Grille[%p].Subs[%p] has maskgrid=%p\n",__func__,ref,Subs[i],sub_gd->mymaskgrid);
   }

   GeoRef_Qualify(ref);

   App_Log(DEBUG,"%s: grtyp     = '%c'\n",__func__, ref->GRTYP[0]);
   App_Log(DEBUG,"%s: grref     = '%c'\n",__func__, ref->RPNHead.GRREF[0]);
   App_Log(DEBUG,"%s: ni        = %d\n",__func__,ref->NX);
   App_Log(DEBUG,"%s: nj        = %d\n",__func__,ref->NY);
   App_Log(DEBUG,"%s: ig1       = %d\n",__func__,ref->RPNHead.IG[X_IG1]);
   App_Log(DEBUG,"%s: ig2       = %d\n",__func__,ref->RPNHead.IG[X_IG2]);
   App_Log(DEBUG,"%s: ig3       = %d\n",__func__,ref->RPNHead.IG[X_IG3]);
   App_Log(DEBUG,"%s: ig4       = %d\n",__func__,ref->RPNHead.IG[X_IG4]);
   App_Log(DEBUG,"%s: ig1ref    = %d\n",__func__,ref->RPNHead.IGREF[X_IG1]);
   App_Log(DEBUG,"%s: ig2ref    = %d\n",__func__,ref->RPNHead.IGREF[X_IG2]);
   App_Log(DEBUG,"%s: ig3ref    = %d\n",__func__,ref->RPNHead.IGREF[X_IG3]);
   App_Log(DEBUG,"%s: ig4ref    = %d\n",__func__,ref->RPNHead.IGREF[X_IG4]);
   App_Log(DEBUG,"%s: NbSub     = %d\n",__func__,ref->NbSub);
   App_Log(DEBUG,"%s: Subs[0]   = %p\n",__func__,ref->Subs[0]);
   App_Log(DEBUG,"%s: Subs[1]   = %p\n",__func__,ref->Subs[1]);

   App_Log(DEBUG,"%s: RPNHead.XG[1] = %f\n",__func__,ref->RPNHead.XG[1]);
   App_Log(DEBUG,"%s: RPNHead.XG[2] = %f\n",__func__,ref->RPNHead.XG[2]);
   App_Log(DEBUG,"%s: RPNHead.XG[3] = %f\n",__func__,ref->RPNHead.XG[3]);
   App_Log(DEBUG,"%s: RPNHead.XG[4] = %f\n",__func__,ref->RPNHead.XG[4]);

   return(ref);
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_Create>
 * Creation     : Avril 2005 J.P. Gauthier - CMC/CMOE
 *
 * But          : Definir le referetiel de type RPN
 *
 * Parametres   :
 *    <NI>      : Dimension en X
 *    <NJ>      : Dimension en Y
 *    <GRTYP>   : Type de grille
 *    <IG1>     : Descripteur IG1
 *    <IG2>     : Descripteur IG2
 *    <IG3>     : Descripteur IG3
 *    <IG4>     : Descripteur IG4
 *    <FID>     : Identificateur du fichier
 *
 * Retour       :
 *
 * Remarques    :
 *
 *---------------------------------------------------------------------------------------------------------------
*/
TGeoRef* GeoRef_Create(int NI,int NJ,char *GRTYP,int IG1,int IG2,int IG3,int IG4,int FID) {

   TGeoRef *ref,*fref;
   int      id;

   ref=GeoRef_New();

   // If not specified, type is X
   if (GRTYP[0]==' ') GRTYP[0]='X';
   
   GeoRef_Size(ref,0,0,NI-1,NJ-1,0);
   ref->RPNHead.GRTYP[0]=ref->GRTYP[0]=GRTYP[0];
   ref->RPNHead.GRTYP[0]=ref->GRTYP[1]=GRTYP[1];

   if ((NI>1 || NJ>1) && GRTYP[0]!='X' && GRTYP[0]!='P' && GRTYP[0]!='V' && ((GRTYP[0]!='Z' && GRTYP[0]!='Y') || FID!=-1)) {

#ifdef HAVE_RMN
      if (GRTYP[1]=='#') {
         //TODO: CHECK For tiled grids (#) we have to fudge the IG3 ang IG4 to 0 since they're used for tile limit
      }

      if (ref->GRTYP[0]=='L' || ref->GRTYP[0]=='A' || ref->GRTYP[0]=='B' || ref->GRTYP[0]=='N' || ref->GRTYP[0]=='S' || ref->GRTYP[0]=='G') {
         // No need to look for grid descriptors
         return(GeoRef_CreateInMemory(NI,NJ,GRTYP," ",IG1,IG2,IG3,IG4,NULL,NULL));
      }
  
      ref->RPNHead.FID=FID;
      ref->RPNHead.IG[X_IG1] = IG1;
      ref->RPNHead.IG[X_IG2] = IG2;
      ref->RPNHead.IG[X_IG3] = IG3;
      ref->RPNHead.IG[X_IG4] = (GRTYP[0]=='#' || GRTYP[0]=='U')?IG4:0;
 
      // This georef already exists
      if (fref=GeoRef_Find(ref)) {
         free(ref);
         GeoRef_Incr(fref);
         return(fref);
      }

      // This is a new georef
      GeoRef_Add(ref);
      if (!RPN_ReadGrid(ref)) {
         // problems with reading grid descriptors
         return(NULL);
      }

      if (GRTYP[0] != 'U') {
         GeoRef_AxisCalcExpandCoeff(ref);
         ref->i1 = 1;
         ref->i2 = ref->NX;
         ref->j1 = 1;
         ref->j2 = ref->NY;
         if (GRTYP[0]!='Y' && GRTYP[0]!='M' && GRTYP[0]!='O') {
            GeoRef_RPNDefXG(ref);           
            GeoRef_AxisCalcNewtonCoeff(ref);
         } else {
            GeoRef_CalcLL(ref);
         }
      }
#endif
   }

   GeoRef_Qualify(ref);

   return(ref);
}

wordint GeoRef_GridGetParams(TGeoRef *Ref,int *NI,int *NJ,char *GRTYP,int *IG1,int *IG2,int *IG3,int *IG4,char *GRREF,int *IG1REF,int *IG2REF,int *IG3REF,int *IG4REF) {
   
   *NI     = Ref->RPNHead.NI;
   *NJ     = Ref->RPNHead.NJ;

   GRTYP[0]  = Ref->RPNHead.GRTYP[0];
   GRTYP[1]  = '\0';
   GRREF[0]  = Ref->RPNHead.GRREF[0];
   GRREF[1]  = '\0';
  
   *IG1    = Ref->RPNHead.IG[X_IG1];
   *IG2    = Ref->RPNHead.IG[X_IG2];
   *IG3    = Ref->RPNHead.IG[X_IG3];
   *IG4    = Ref->RPNHead.IG[X_IG4];
   *IG1REF = Ref->RPNHead.IGREF[X_IG1];
   *IG2REF = Ref->RPNHead.IGREF[X_IG2];
   *IG3REF = Ref->RPNHead.IGREF[X_IG3];
   *IG4REF = Ref->RPNHead.IGREF[X_IG4];

   return(0);
}

int GEM_grid_param(int *F_bsc_base,int *F_bsc_ext1,int *F_extension ,int F_maxcfl,float *F_lonr,float *F_latr,int *F_ni,int *F_nj,float *F_dx,float *F_dy,double *F_x0_8,double *F_y0_8,double *F_xl_8,double *F_yl_8,int F_overlap,int F_yinyang_L) {

   double delta_8;
   int iref,jref;
  
   // basic global lateral boundary conditions width
   *F_bsc_base = 5;
   if (F_yinyang_L) *F_bsc_base=*F_bsc_base+1;

   // added points for proper de-staggering of u,v at physics interface
   *F_bsc_ext1 = 2;

   // total extension to user specified grid configuration
   *F_extension= F_maxcfl + *F_bsc_base + *F_bsc_ext1;

   if (F_yinyang_L) {

      *F_x0_8 =   45.0 - 3.0*F_overlap;
      *F_xl_8 =  315.0 + 3.0*F_overlap;
      *F_y0_8 = -45.0  -     F_overlap;
      *F_yl_8 =  45.0  +     F_overlap;
      
      delta_8  = ((*F_xl_8)-(*F_x0_8))/(*F_ni-1);
      *F_dx   = delta_8;
      *F_x0_8 = *F_x0_8 - (*F_extension)*delta_8;
      *F_xl_8 = *F_xl_8 + (*F_extension)*delta_8;
      
      delta_8  = ((*F_yl_8)-(*F_y0_8))/(*F_nj-1);
      *F_dy   = delta_8;
      *F_y0_8 = *F_y0_8 - (*F_extension)*delta_8;
      *F_yl_8 = *F_yl_8 + (*F_extension)*delta_8;
      
      *F_ni   = *F_ni + 2  *(*F_extension);
      *F_nj   = *F_nj + 2  *(*F_extension);

   } else {

      iref = *F_ni / 2 + (*F_extension);
      if ((*F_ni)%2==0) {
         *F_lonr = *F_lonr - (*F_dx)/2.0;
      } else {
         iref = iref + 1;
      }
      jref = *F_nj / 2 + (*F_extension);
      if ((*F_nj)%2==0) {
         *F_latr = *F_latr - (*F_dy)/2.0;
      } else {
         jref = *F_nj / 2 + (*F_extension) + 1;
      }
      
      *F_ni   = *F_ni + 2*(*F_extension);
      *F_nj   = *F_nj + 2*(*F_extension);
      *F_x0_8 = *F_lonr - (iref-1) * (*F_dx);
      *F_y0_8 = *F_latr - (jref-1) * (*F_dy);
      *F_xl_8 = *F_x0_8 + (*F_ni  -1) * (*F_dx);
      *F_yl_8 = *F_y0_8 + (*F_nj  -1) * (*F_dy);
      if (*F_x0_8 < 0.) *F_x0_8=*F_x0_8+360.0;
      if (*F_xl_8 < 0.) *F_xl_8=*F_xl_8+360.0;

      if (*F_x0_8 < 0.) {
         fprintf(stderr,"Longitude of WEST %f < 0.0\n",*F_x0_8);
         return(0);
      }
      if (*F_y0_8 < -90.) {
         fprintf(stderr,"Latitude of SOUTH %f < 0.0\n",*F_y0_8);
         return(0);
      }
      if (*F_xl_8 > 360.) {
         fprintf(stderr,"Longitude of EAST %f < 0.0\n",*F_xl_8);
         return(0);
      }
      if (*F_yl_8 > 90.) {
         fprintf(stderr,"Latitude of NORTH %f < 0.0\n",*F_yl_8);
         return(0);
      }
   }
   return(1);
}

void GEM_hgrid4(double *F_xgi_8,double *F_ygi_8,int F_Grd_ni,int F_Grd_nj,float *F_Grd_dx,float *F_Grd_dy,double F_Grd_x0_8,double F_Grd_xl_8,double F_Grd_y0_8,double F_Grd_yl_8, int F_Grd_yinyang_L){

   int i;
   double delta_8;

   delta_8 = (F_Grd_xl_8-F_Grd_x0_8)/(F_Grd_ni-1);
   F_xgi_8[0] = F_Grd_x0_8;
   F_xgi_8[F_Grd_ni-1] = F_Grd_xl_8;
   for(i=1;i<F_Grd_ni-1;i++) F_xgi_8[i]= F_Grd_x0_8 + i*delta_8;

   delta_8 = (F_Grd_yl_8-F_Grd_y0_8)/(F_Grd_nj-1);
   F_ygi_8[0] = F_Grd_y0_8;
   F_ygi_8[F_Grd_nj-1] = F_Grd_yl_8;
   for(i=1;i<F_Grd_nj-1;i++) F_ygi_8[i]= F_Grd_y0_8 + i*delta_8;

   if (F_Grd_yinyang_L) {
      *F_Grd_dx   = fabs(F_xgi_8[1]-F_xgi_8[0]);
      *F_Grd_dy   = fabs(F_ygi_8[1]-F_ygi_8[0]);
   }
}

/*-------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_SetZE>
 * Creation     : Avril 2005 J.P. Gauthier - CMC/CMOE
 *
 * But          : Definir le referentiel de type RPN ZE
 *
 * Parametres   :
 *    <Ref>    : Georef definition
 *    <NI>      : Dimension en X
 *    <NJ>      : Dimension en Y
 *    <DX>      : Resolution en X
 *    <DY>      : Resolution en Y
 *    <LatR>    : 
 *    <LonR>    : 
 *    <MaxCFL>  : 
 *    <XLat1>   : Latitude centrale
 *    <XLon1>   : Longitude centrale
 *    <XLat2>   : Latitude de l'axe de rotation
 *    <XLon2>   : Longitude de l'axe de rotation
 *
 * Retour       :
 *
 * Remarques    :
 *
 *---------------------------------------------------------------------------------------------------------------
*/
TGeoRef* GeoRef_SetZE(TGeoRef *Ref,int NI,int NJ,float DX,float DY,float LatR,float LonR,int MaxCFL,float XLat1,float XLon1,float XLat2,float XLon2) {

#ifdef HAVE_RMN
   int    bsc_base,bsc_ext1,extension;
   double x0,x1,y0,y1;

   if (!Ref) {
      return(NULL);
   }
 
   f77name(cxgaig)("E",&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4],&XLat1,&XLon1,&XLat2,&XLon2);
   f77name(cigaxg)("E",&XLat1,&XLon1,&XLat2,&XLon2,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);
   
   GEM_grid_param(&bsc_base,&bsc_ext1,&extension,MaxCFL,&LonR,&LatR,&NI,&NJ,&DX,&DY,&x0,&y0,&x1,&y1,-1,FALSE);
 
   if (NI!=Ref->NX+1 || NJ!=Ref->NY+1) {
      Ref->AX=realloc(Ref->AX,NI*sizeof(double));
      Ref->AY=realloc(Ref->AY,NJ*sizeof(double));
   }

   // f77name(set_gemhgrid4)(Ref->AX,Ref->AY,&NI,&NJ,&DX,&DY,&x0,&x1,&y0,&y1,FALSE);
   GEM_hgrid4(Ref->AX,Ref->AY,NI,NJ,&DX,&DY,x0,x1,y0,y1,FALSE);
           
   GeoRef_Set(Ref,NI,NJ,"Z","E",Ref->RPNHead.IG[X_IG1],Ref->RPNHead.IG[X_IG2],Ref->RPNHead.IG[X_IG3],Ref->RPNHead.IG[X_IG4],Ref->AX,Ref->AY);
   GeoRef_RPNDefXG(Ref);
   GeoRef_AxisDefine(Ref,Ref->AX,Ref->AY);
   GeoRef_Qualify(Ref);
 #endif
  
   return(Ref);
}
