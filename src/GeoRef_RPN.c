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

double   GeoRef_RPNDistance(TGeoRef *Ref,double X0,double Y0,double X1, double Y1);
int      GeoRef_RPNValue(TGeoRef *Ref,TDef *Def,TDef_InterpR Interp,int C,double X,double Y,double Z,double *Length,double *ThetaXY);
int      GeoRef_RPNProject(TGeoRef *Ref,double X,double Y,double *Lat,double *Lon,int Extrap,int Transform);
int      GeoRef_RPNUnProject(TGeoRef *Ref,double *X,double *Y,double Lat,double Lon,int Extrap,int Transform);

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_RPNDistance>
 * Creation     : Mars 2007 J.P. Gauthier - CMC/CMOE
 *
 * But          : Calculer la distance entre deux points.
 *
 * Parametres    :
 *   <Ref>      : Pointeur sur la reference geographique
 *   <X0>        : coordonnee en X dans la projection/grille
 *   <Y0>        : coordonnee en Y dans la projection/grille
 *   <X0>        : coordonnee en X dans la projection/grille
 *   <Y0>        : coordonnee en Y dans la projection/grille
 *
 * Retour       : Distance
 *
 * Remarques   :
 *
 *---------------------------------------------------------------------------------------------------------------
*/
double GeoRef_RPNDistance(TGeoRef *Ref,double X0,double Y0,double X1, double Y1) {

#ifdef HAVE_RMN
   double i[2],j[2],lat[2],lon[2];

   if (Ref->Type&GRID_EZ) {
      i[0]=X0+1.0;
      j[0]=Y0+1.0;
      i[1]=X1+1.0;
      j[1]=Y1+1.0;

      GeoRef_XY2LL(REFGET(Ref),lat,lon,i,j,2);

      X0=DEG2RAD(lon[0]);
      X1=DEG2RAD(lon[1]);
      Y0=DEG2RAD(lat[0]);
      Y1=DEG2RAD(lat[1]);

      return(DIST(0.0,Y0,X0,Y1,X1));
   }
#endif
   return(hypot(X1-X0,Y1-Y0));
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_RPNValue>
 * Creation     : Mars 2005 J.P. Gauthier - CMC/CMOE
 *
 * But          : Extraire la valeur d'une matrice de donnees.
 *
 * Parametres    :
 *   <Ref>      : Pointeur sur la reference geographique
 *   <Def>       : Pointeur sur la definition de la donnee
 *   <Interp>    : Mode d'interpolation
 *   <C>         : Composante
 *   <X>         : coordonnee en X dans la projection/grille
 *   <Y>         : coordonnee en Y dans la projection/grille
 *   <Z>         : coordonnee en Z dans la projection/grille
 *   <Length>    : Module interpolee
 *   <ThetaXY>   : Direction interpolee
 *
 * Retour       : Inside (1 si a l'interieur du domaine).
 *
 * Remarques   :
 *
 *---------------------------------------------------------------------------------------------------------------
*/

int GeoRef_RPNValue(TGeoRef *Ref,TDef *Def,TDef_InterpR Interp,int C,double X,double Y,double Z,double *Length,double *ThetaXY) {

   Vect3d       b,v;
   double       x,y;
   float        valf,valdf;
   void        *p0,*p1;
   int          mem,ix,iy,n;
   unsigned int idx;

   *Length=Def->NoData;

#ifdef HAVE_RMN
   // In case of triangle meshe
   if (Ref->GRTYP[0]=='M') {
      if (C<Def->NC && X>=0 && Y>=0) {
         b[0]=X-(int)X;
         b[1]=Y-(int)Y;
         b[2]=1.0-b[0]-b[1];
         ix=(int)X;

         if (Interp==IR_NEAREST) {
            n=(b[0]>b[1]?(b[0]>b[2]?0:2):(b[1]>b[2]?1:2));
            Def_Get(Def,C,Ref->Idx[ix+n],x);
         } else {
            Def_Get(Def,C,Ref->Idx[ix],v[0]);
            Def_Get(Def,C,Ref->Idx[ix+1],v[1]);
            Def_Get(Def,C,Ref->Idx[ix+2],v[2]);

            x=Bary_InterpV(b,v);
         }
         *Length=x;
         
         if (Def->Data[1] && !C) {
            if (Interp==IR_NEAREST) {
               n=(b[0]>b[1]?(b[0]>b[2]?0:2):(b[1]>b[2]?1:2));
               Def_Get(Def,1,Ref->Idx[ix+n],y);
            } else {
               Def_Get(Def,1,Ref->Idx[ix],v[0]);
               Def_Get(Def,1,Ref->Idx[ix+1],v[1]);
               Def_Get(Def,1,Ref->Idx[ix+2],v[2]);

               y=Bary_InterpV(b,v);
            }
            *Length=hypot(x,y);
            *ThetaXY=180+RAD2DEG(atan2(x,y));
         }
         return(TRUE);
      } else {
         return(FALSE);
      }
   }

   // Si on est a l'interieur de la grille ou que l'extrapolation est activee
   if (C<Def->NC && X>=(Ref->X0-0.5) && Y>=(Ref->Y0-0.5) && Z>=0 && X<(Ref->X1+0.5) && Y<(Ref->Y1+0.5) && Z<=Def->NK-1) {

      // Index memoire du niveau desire
      mem=Def->NIJ*(int)Z;

      ix=lrint(X);
      iy=lrint(Y);
      idx=iy*Def->NI+ix;
      
      // Check for mask
      if (Def->Mask && !Def->Mask[mem+idx]) {
         return(FALSE);
      }

      // Point cloud
      if (Ref->GRTYP[0]=='Y' || Ref->GRTYP[0]=='P' || Def->NI==1 || Def->NJ==1) {
         mem+=idx;
         Def_Get(Def,C,mem,x);
         if (Def->Data[1] && !C) {
            Def_Get(Def,1,mem,y);
            if (ThetaXY) *ThetaXY=180+RAD2DEG(atan2(x,y));
            *Length=hypot(x,y);
         } else {
            *Length=x;   
         }
         return(TRUE);
      } 

      // Check for nodata in linear interpolation 
      if (Def->Type>=TD_Int64 && Interp!=IR_NEAREST) {
         int idxs[4],dx,dy;
         double vals[4];
         
         dx=X;
         dy=Y;
         idxs[0]=mem+dy*Def->NI+dx;
         idxs[1]=(dx<Ref->X1)?idxs[0]+1:idxs[0];
         idxs[2]=(dy<Ref->Y1)?idxs[0]+Def->NI:idxs[0];
         idxs[3]=(dy<Ref->Y1 && dx<Ref->X1)?idxs[0]+Def->NI+1:idxs[0];
         
         Def_GetQuad(Def,C,idxs,vals);

         // If either value is nodata then interpolation will be nodata as well
         if (!DEFVALID(Def,vals[0]) || !DEFVALID(Def,vals[1]) || !DEFVALID(Def,vals[2]) || !DEFVALID(Def,vals[3])) 
            return(FALSE);        
      }
     
      // XSection
      if (Ref->GRTYP[0]=='V') {
         if (Def->Data[1]) {
            Def_GetMod(Def,FIDX2D(Def,ix,iy),*Length);
         } else {
            *Length=VertexVal(Def,-1,X,Y,0.0);
         }
         return(TRUE);
      }

      // Unstructured or not referenced
      if (Ref->GRTYP[0]=='X' || Ref->GRTYP[0]=='O') {
         if (Def->Type<=TD_Int64 || Interp==IR_NEAREST || (X==ix && Y==iy)) {
            mem+=idx;

           // Pour un champs vectoriel
            if (Def->Data[1] && !C) {
               Def_Get(Def,0,mem,x);
               Def_Get(Def,1,mem,y);
               if (ThetaXY) *ThetaXY=180+RAD2DEG(atan2(x,y)-GeoRef_GeoDir(Ref,X,Y));
               *Length=hypot(x,y);
            } else {
               Def_GetMod(Def,mem,*Length);              
            }
         } else {
            // Pour un champs vectoriel
            if (Def->Data[1] && !C) {
               x=VertexVal(Def,0,X,Y,Z);
               y=VertexVal(Def,1,X,Y,Z);
               if (ThetaXY) *ThetaXY=180+RAD2DEG(atan2(x,y)-GeoRef_GeoDir(Ref,X,Y));
               *Length=hypot(x,y);
            } else {
               *Length=VertexVal(Def,-1,X,Y,Z);               
            }
         }
         return(TRUE);
      }
   
      // RPN grid
      if (Ref->Type&GRID_EZ) {         
         if (Def->Type==TD_Float32 && Def->Data[1] && !C) { 
            x=X+1.0;
            y=Y+1.0;

            if (Interp==IR_NEAREST) {
               x=lrint(x);
               y=lrint(y);
            }
            
            Def_Pointer(Def,0,mem,p0);
            Def_Pointer(Def,1,mem,p1);
            GeoRef_XYWDVal(REFGET(Ref),&valf,&valdf,p0,p1,&x,&y,1);

            // If it's 3D, use the mode for speed since GeoRef_XYWDVal only uses 2D
            if (Def->Data[2])
               GeoRef_XYVal(REFGET(Ref),&valf,(float*)&Def->Mode[mem],&x,&y,1);
               *Length=valf;
            if (ThetaXY)
               *ThetaXY=valdf;
         } else {            
            if (Interp==IR_NEAREST) {
               mem+=idx;
               Def_Get(Def,C,mem,*Length);
            } else {
               *Length=VertexVal(Def,C,X,Y,Z);
//            Def_Pointer(Def,0,mem,p0);
//            x=X+=1.0;y=Y+=1.0;
//                GeoRef_XYVal(REFGET(Ref),&valf,p0,&x,&y,1);
//               fprintf(stderr,"----- %.10e %.10e   ---> %.10e\n",*Length,valf,*Length-valf);
            }
            if (ThetaXY)
               *ThetaXY=0.0;
         }
      }
      return(TRUE);
   }
#endif
   return(FALSE);
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_RPNProject>
 * Creation     : Mars 2005 J.P. Gauthier - CMC/CMOE
 *
 * But          : Projeter une coordonnee de projection en latlon.
 *
 * Parametres    :
 *   <Ref>      : Pointeur sur la reference geographique
 *   <X>         : coordonnee en X dans la projection/grille
 *   <Y>         : coordonnee en Y dans la projection/grille
 *   <Lat>       : Latitude
 *   <Lon>       : Longitude
 *   <Extrap>    : Extrapolation hors grille
 *   <Transform> : Appliquer la transformation
 *
 * Retour       : Inside (1 si a l'interieur du domaine).
 *
 * Remarques   :
 *
 *---------------------------------------------------------------------------------------------------------------
*/
int GeoRef_RPNProject(TGeoRef *Ref,double X,double Y,double *Lat,double *Lon,int Extrap,int Transform) {

#ifdef HAVE_RMN
//TODO: double   GeoRef_XY2LL(REFGET(Ref),Lat,Lon,&X,&Y,1);
#endif

   return(1);
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_RPNUnProject>
 * Creation     : Mars 2005 J.P. Gauthier - CMC/CMOE
 *
 * But          : Projeter une latlon en position grille.
 *
 * Parametres    :
 *   <Ref>      : Pointeur sur la reference geographique
 *   <X>         : coordonnee en X dans la projection/grille
 *   <Y>         : coordonnee en Y dans la projection/grille
 *   <Lat>       : Latitude
 *   <Lon>       : Longitude
 *   <Extrap>    : Extrapolation hors grille
 *   <Transform> : Appliquer la transformation
 *
 * Retour       : Inside (1 si a l'interieur du domaine).
 *
 * Remarques   :
 *
 *---------------------------------------------------------------------------------------------------------------
*/                     
int GeoRef_RPNUnProject(TGeoRef *Ref,double *X,double *Y,double Lat,double Lon,int Extrap,int Transform) {

   int     idx,idxs[8];
   double  dists[8];

   *X=-1.0;
   *Y=-1.0;

   // Invalid coordinates ?
   if (Lat>90.0 || Lat<-90.0 || Lon==-999.0) 
     return(FALSE);

   Lon=GeoRef_Lon(Ref,Lon);
   
#ifdef HAVE_RMN
   if (Ref->Type&GRID_SPARSE) {      
      if (Ref->AX && Ref->AY) {
         if (Ref->GRTYP[0]=='Y') {
            // Get nearest point
            if (GeoRef_Nearest(Ref,Lon,Lat,&idx,dists,1,0.0)) {
               if (dists[0]<1.0) {
                  *Y=(int)(idx/Ref->NX);
                  *X=idx-(*Y)*Ref->NX;
                  return(TRUE);
               }         
            }
         }
      } 
      return(FALSE);
   } else {

      // Extraire la valeur du point de grille
      GeoRef_LL2XY(REFGET(Ref),X,Y,&Lat,&Lon,1);

//TODO: Set CIndex default to 1

      // Fix for G grid 0-360 1/5 gridpoint problem
//TODO: Check to put in fortran ez whatever
      if (Ref->GRTYP[0]=='G' && *X>Ref->X1+0.5) *X-=(Ref->X1+1);
  }
#endif
   return(TRUE);
}

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
         if (Ref->RPNHead.GRREF[0] == 'E') {
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

TGeoRef* GeoRef_RPNCreateInMemory(int NI,int NJ,char* GRTYP,char* GRREF,int IG1,int IG2,int IG3,int IG4,float* AX,float* AY) {
   
   TGeoRef* ref,*fref;

   ref = GeoRef_New();

   ref->RPNHead.GRTYP[0]=ref->GRTYP[0] = GRTYP[0];
   ref->RPNHead.GRTYP[1]=ref->GRTYP[1] = '\0';
   ref->RPNHead.GRREF[0] = GRREF?GRREF[0]:'\0';
   ref->RPNHead.GRREF[1] = '\0';
   ref->RPNHead.NI=ref->NX = NI;
   ref->RPNHead.NJ=ref->NY = NJ;
   ref->RPNHead.IG[X_IG1] = IG1;
   ref->RPNHead.IG[X_IG2] = IG2;
   ref->RPNHead.IG[X_IG3] = IG3;
   ref->RPNHead.IG[X_IG4] = IG4;
   ref->i1 = 1;
   ref->i2 = NI;
   ref->j1 = 1;
   ref->j2 = NJ;

   switch (GRTYP[0]) {
      case 'Z':
         f77name(cigaxg)(ref->RPNHead.GRREF,&ref->RPNHead.XGREF[X_LAT1],&ref->RPNHead.XGREF[X_LON1],&ref->RPNHead.XGREF[X_LAT2],&ref->RPNHead.XGREF[X_LON2],&ref->RPNHead.IGREF[X_IG1],&ref->RPNHead.IGREF[X_IG2],&ref->RPNHead.IGREF[X_IG3],&ref->RPNHead.IGREF[X_IG4]);
      case '#':
      case 'Y':
         ref->AX = AX;
         ref->AY = AY;
         break;
   }

   GeoRef_Size(ref,0,0,NI-1,NJ-1,0);

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

   // TODO: Check for sub-grids (U grids can have sub grids)
   //ref->NbId = GRTYP[0]=='U'? (ref->NbSub==0? 1 : ref->NbSub) : 1;

   ref->Project=GeoRef_RPNProject;
   ref->UnProject=GeoRef_RPNUnProject;
   ref->Value=(TGeoRef_Value*)GeoRef_RPNValue;
   ref->Distance=GeoRef_RPNDistance;

   return(ref);
}

TGeoRef* GeoRef_RPNCreateYY(int NI,int NJ,char *GRTYP,char *GRREF,int VerCode,int NbSub,TGeoRef **Subs) {

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
 * Nom          : <GeoRef_RPNCreate>
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
TGeoRef* GeoRef_RPNCreate(int NI,int NJ,char *GRTYP,int IG1,int IG2,int IG3,int IG4,int FID) {

   TGeoRef *ref,*fref;
   int      id;

   ref=GeoRef_New();

   // If not specified, type is X
   if (GRTYP[0]==' ') GRTYP[0]='X';
   
   GeoRef_Size(ref,0,0,NI-1,NJ-1,0);
   ref->RPNHead.GRTYP[0]=ref->GRTYP[0]=GRTYP[0];
   ref->RPNHead.GRTYP[0]=ref->GRTYP[1]=GRTYP[1];
   ref->Project=GeoRef_RPNProject;
   ref->UnProject=GeoRef_RPNUnProject;
   ref->Value=(TGeoRef_Value*)GeoRef_RPNValue;
   ref->Distance=GeoRef_RPNDistance;
   ref->Height=NULL;

   if ((NI>1 || NJ>1) && GRTYP[0]!='X' && GRTYP[0]!='P' && GRTYP[0]!='M' && GRTYP[0]!='V' && ((GRTYP[0]!='Z' && GRTYP[0]!='Y') || FID!=-1)) {

#ifdef HAVE_RMN
      if (GRTYP[1]=='#') {
         //TODO: CHECK For tiled grids (#) we have to fudge the IG3 ang IG4 to 0 since they're used for tile limit
      }

      if (GRTYP[0]!='#' && GRTYP[0]!='Y' && GRTYP[0]!='Z' && GRTYP[0]!='O' && GRTYP[0]!='U' && GRTYP[0]!=' ') {
         // No need to look for grid descriptors
         return(GeoRef_RPNCreateInMemory(NI,NJ,GRTYP," ",IG1,IG2,IG3,IG4,NULL,NULL));
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
      if (!RPN_ReadGrid(ref,&ref->RPNHead)) {
         // problems with reading grid descriptors
         return(NULL);
      }
      
      if (GRTYP[0] != 'U') {
         GeoRef_AxisCalcExpandCoeff(ref);
         ref->i1 = 1;
         ref->i2 = ref->NX;
         ref->j1 = 1;
         ref->j2 = ref->NY;
         if (GRTYP[0] != 'Y' && GRTYP[0] != 'O') {
            GeoRef_RPNDefXG(ref);
            GeoRef_AxisCalcNewtonCoeff(ref);
         } else {
            GeoRef_CalcLL(ref);
         }
      }
#endif
   }

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

void GEM_hgrid4(float *F_xgi_8,float *F_ygi_8,int F_Grd_ni,int F_Grd_nj,float *F_Grd_dx,float *F_Grd_dy,double F_Grd_x0_8,double F_Grd_xl_8,double F_Grd_y0_8,double F_Grd_yl_8, int F_Grd_yinyang_L){

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
 * Nom          : <GeoRef_RPNGridZE>
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
TGeoRef* GeoRef_RPNGridZE(TGeoRef *Ref,int NI,int NJ,float DX,float DY,float LatR,float LonR,int MaxCFL,float XLat1,float XLon1,float XLat2,float XLon2) {

#ifdef HAVE_RMN
   int    ig1,ig2,ig3,ig4;
   char   gxtyp='E';
   int    bsc_base,bsc_ext1,extension,err;
   int    maxcfl;
   double x0,x1,y0,y1;
   float latr,lonr;

   if (!Ref) {
      return(NULL);
   }
 
   f77name(cxgaig)("E",&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4],&XLat1,&XLon1,&XLat2,&XLon2);
   f77name(cigaxg)("E",&XLat1,&XLon1,&XLat2,&XLon2,&Ref->RPNHead.IG[X_IG1],&Ref->RPNHead.IG[X_IG2],&Ref->RPNHead.IG[X_IG3],&Ref->RPNHead.IG[X_IG4]);
   
   GEM_grid_param(&bsc_base,&bsc_ext1,&extension,MaxCFL,&LonR,&LatR,&NI,&NJ,&DX,&DY,&x0,&y0,&x1,&y1,-1,FALSE);
 
   if (NI!=Ref->NX+1 || NJ!=Ref->NY+1) {
      Ref->AX=realloc(Ref->AX,NI*sizeof(float));
      Ref->AY=realloc(Ref->AY,NJ*sizeof(float));

      GeoRef_Size(Ref,0,0,NI-1,NJ-1,0);
   }

   //   f77name(set_gemhgrid4)(Ref->AX,Ref->AY,&NI,&NJ,&DX,&DY,&x0,&x1,&y0,&y1,FALSE);
   GEM_hgrid4(Ref->AX,Ref->AY,NI,NJ,&DX,&DY,x0,x1,y0,y1,FALSE);
           
 //TODO: Merge with EZ  
 //  Ref->Ids[0]=GeoRef_RPNCreateInMemory(NI,NJ,"Z","E",Ref->IG[X_IG1],Ref->IG[X_IG2],Ref->IG[X_IG3],Ref->IG[X_IG4],Ref->AX,Ref->AY);
   
   Ref->NbSub=1;
   Ref->GRTYP[0]='Z';
   Ref->GRTYP[1]='E';
   Ref->GRTYP[2]='\0';
   Ref->Project=GeoRef_RPNProject;
   Ref->UnProject=GeoRef_RPNUnProject;
   Ref->Value=(TGeoRef_Value*)GeoRef_RPNValue;
   Ref->Distance=GeoRef_RPNDistance;
   Ref->Height=NULL;
 #endif
  
   return(Ref);
}
