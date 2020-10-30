/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Fonctions et definitions relatives aux fichiers standards et rmnlib
 * Fichier      : GeoRef_Type_W.h
 * Creation     : October 2020 - J.P. Gauthier
 *
 * Description:
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
#include "GeoRef.h"
#include "Def.h"
#include "Vertex.h"

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_WKTValue>
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
int GeoRef_WKTValue(TGeoRef *Ref,TDef *Def,TDef_InterpR Interp,int C,double X,double Y,double Z,double *Length,double *ThetaXY){

   double       x,y,d,ddir=0.0;
   int          valid=0,mem,ix,iy;
   unsigned int idx;

  *Length=Def->NoData;
   d=1.0;

   //i on est a l'interieur de la grille ou que l'extrapolation est activee
   if (C<Def->NC && X>=(Ref->X0-d) && Y>=(Ref->Y0-d) && Z>=0 && X<=(Ref->X1+d) && Y<=(Ref->Y1+d) && Z<=Def->NK-1) {

      X-=Ref->X0;
      Y-=Ref->Y0;
      DEFCLAMP(Def,X,Y);

      // Index memoire du niveau desire
      mem=Def->NI*Def->NJ*(int)Z;

      ix=lrint(X);
      iy=lrint(Y);
      idx=iy*Def->NI+ix;

      // Check for mask
      if (Def->Mask && !Def->Mask[idx]) {
         return(valid);
      }
      
      // Reproject vector orientation by adding grid projection's north difference
      if (Def->Data[1] && Ref->Type&GRID_ROTATED) { 
         ddir=GeoRef_GeoDir(Ref,X,Y);
      }

      if (Def->Type<=9 || Interp==IR_NEAREST || (X==ix && Y==iy)) {
         mem+=idx;
         Def_GetMod(Def,mem,*Length);

         // Pour un champs vectoriel
         if (Def->Data[1] && ThetaXY) {
            Def_Get(Def,0,mem,x);
            Def_Get(Def,1,mem,y);
            *ThetaXY=180+RAD2DEG(atan2(x,y)-ddir);
         }
      } else {
         *Length=VertexVal(Def,-1,X,Y,Z);
         // Pour un champs vectoriel
         if (Def->Data[1] && ThetaXY) {
            x=VertexVal(Def,0,X,Y,Z);
            y=VertexVal(Def,1,X,Y,Z);
            *ThetaXY=180+RAD2DEG(atan2(x,y)-ddir);
         }
      }
      valid=DEFVALID(Def,*Length);
   }
   
   return(valid);
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_WKTProject>
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
int GeoRef_WKTProject(TGeoRef *Ref,double X,double Y,double *Lat,double *Lon,int Extrap,int Transform) {

   GeoRef_XY2LL(REFGET(Ref),&X,&Y,Lat,Lon,1);

   return(0);
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_WKTUnProject>
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
int GeoRef_WKTUnProject(TGeoRef *Ref,double *X,double *Y,double Lat,double Lon,int Extrap,int Transform) {

   GeoRef_LL2XY(REFGET(Ref),X,Y,&Lat,&Lon,1);

   return(0);
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_WKTSet>
 * Creation     : Juin 2004 J.P. Gauthier - CMC/CMOE
 *
 * But          : Definir les fonctions de transformations WKT
 *
 * Parametres   :
 *   <Ref>     : Pointeur sur la reference geographique
 *   <String>   : Description de la projection
 *   <Geometry> : Geometrie d'ou extraire la reference spatiale (optionel=NULL)
 *
 * Retour       :
 *
 * Remarques    :
 *
 *---------------------------------------------------------------------------------------------------------------
*/
TGeoRef *GeoRef_SetW(TGeoRef *Ref,char *String,double *Transform,double *InvTransform,OGRSpatialReferenceH Spatial) {

#ifdef HAVE_GDAL
   static OGRSpatialReferenceH llref=NULL;
   char                      *string=NULL;

   if (String && String[0]!='\0') {
      string=strdup(String);
      strtrim(string,' ');
   }

//TODO: Reenable ?   GeoRef_Clear(Ref,0);

   if (Transform || InvTransform) {
      if (!Ref->Transform)
         Ref->Transform=(double*)calloc(6,sizeof(double));
      if (!Ref->InvTransform)
         Ref->InvTransform=(double*)calloc(6,sizeof(double));
   }
   
   if (Transform) {
      memcpy(Ref->Transform,Transform,6*sizeof(double));
   } else {
      if (!InvTransform || !GDALInvGeoTransform(InvTransform,Ref->Transform)) {
         if (Ref->Transform) {
            free(Ref->Transform);
            Ref->Transform=NULL;
         }
      }
   }

   if (InvTransform) {
      memcpy(Ref->InvTransform,InvTransform,6*sizeof(double));
   } else {
      if (!Transform || !GDALInvGeoTransform(Transform,Ref->InvTransform)) {
         if (Ref->InvTransform) {
            free(Ref->InvTransform);
            Ref->InvTransform=NULL;
         }
      }
   }

   if (Spatial) {
      Ref->Spatial=OSRClone(Spatial);
      OSRExportToWkt(Ref->Spatial,&string);
   } else if (string) {
      Ref->Spatial=OSRNewSpatialReference(NULL);
      if (OSRSetFromUserInput(Ref->Spatial,string)==OGRERR_FAILURE) {
        App_Log(WARNING,"%s: Unable to create spatial reference\n",__func__);
        return(NULL);
      }
   } else {
      string=strdup(REFDEFAULT);
      Ref->Spatial=OSRNewSpatialReference(string);
   }

   if (Ref->String)
      free(Ref->String);
   Ref->String=string;

   if (Ref->Spatial) {
      if (!llref) {
         // Create global latlon reference on perfect sphere
         llref=OSRNewSpatialReference(NULL);
         OSRSetFromUserInput(llref,"EPSG:4047");
      }

      if (llref) {
         // Create forward/backward tranformation functions
         Ref->Function=OCTNewCoordinateTransformation(Ref->Spatial,llref);
         Ref->InvFunction=OCTNewCoordinateTransformation(llref,Ref->Spatial);
      }
  } else {
      App_Log(WARNING,"%s: Unable to get spatial reference\n",__func__);
      return(NULL);
   }

   Ref->Project=GeoRef_WKTProject;
   Ref->UnProject=GeoRef_WKTUnProject;
   Ref->Value=(TGeoRef_Value*)GeoRef_WKTValue;
   Ref->Height=NULL;

   return(Ref);
#else
   App_Log(ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
   return(NULL);
#endif
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_WKTCreate>
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
 *
 * Retour       :
 *
 * Remarques    :
 *
 *---------------------------------------------------------------------------------------------------------------
*/
TGeoRef *GeoRef_CreateW(int ni,int nj,char *grtyp,int ig1,int ig2,int ig3,int ig4,char *String,double *Transform,double *InvTransform,OGRSpatialReferenceH Spatial) {

   TGeoRef *ref;

   ref=GeoRef_New();
   GeoRef_Size(ref,0,0,ni>0?ni-1:0,nj>0?nj-1:0,0);
   if (!GeoRef_SetW(ref,String,Transform,InvTransform,Spatial)) {
      return(NULL);
   }
   
   if (grtyp) {
      ref->GRTYP[0]=grtyp[0];
      ref->GRTYP[1]=grtyp[1];
   } else {
      ref->GRTYP[0]='W';
      ref->GRTYP[1]='\0';
   }
   ref->RPNHead.IG[X_IG1]=ig1;
   ref->RPNHead.IG[X_IG2]=ig2;
   ref->RPNHead.IG[X_IG3]=ig3;
   ref->RPNHead.IG[X_IG4]=ig4;
   
   return(ref);
}

#ifdef HAVE_GDAL
static inline int GeoRef_WKTRotate(TRotationTransform *T,double *Lat,double *Lon) {

   double lat,lon,x,y,z,xr,yr,zr;
   
   lon = DEG2RAD(*Lon);
   lat = DEG2RAD(*Lat);

   // Convert from spherical to cartesian coordinates
   x = cos(lon)*cos(lat); 
   y = sin(lon)*cos(lat);
   z = sin(lat);

   xr = T->CosTheta*T->CosPhi*x + T->CosTheta*T->SinPhi*y + T->SinTheta*z;
   yr = -T->SinPhi*x + T->CosPhi*y;
   zr = -T->SinTheta*T->CosPhi*x - T->SinTheta*T->SinPhi*y + T->CosTheta*z;
      
   // Convert cartesian back to spherical coordinates
   *Lon = RAD2DEG(atan2(yr,xr)); 
   *Lat = RAD2DEG(asin(zr));
   
   return(TRUE);
}
#endif

#ifdef HAVE_GDAL
static inline int GeoRef_WKTUnRotate(TRotationTransform *T,double *Lat,double *Lon) {

   double lat,lon,x,y,z,xr,yr,zr;
   
   lon = DEG2RAD(*Lon);
   lat = DEG2RAD(*Lat);

   // Convert from spherical to cartesian coordinates
   x = cos(lon)*cos(lat); 
   y = sin(lon)*cos(lat);
   z = sin(lat);

   xr = T->CosTheta*T->CosPhi*x + -T->SinPhi*y + -T->SinTheta*T->CosPhi*z;
   yr = -T->CosTheta*-T->SinPhi*x + T->CosPhi*y - -T->SinTheta*-T->SinPhi*z;
   zr = T->SinTheta*x + T->CosTheta*z;
      
   // Convert cartesian back to spherical coordinates
   *Lon = RAD2DEG(atan2(yr,xr)); 
   *Lat = RAD2DEG(asin(zr));

   return(TRUE);
}
#endif

int GeoRef_XY2LL_W(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {

#ifdef HAVE_GDAL
   double x,y,z=0.0,d;
   int    n,ok;

   for(n=0;n<Nb;n++) {

      // Transform the point into georeferenced coordinates 
      x=X[n];
      y=Y[n];
      if (Ref->Options.Transform) {
         if (Ref->Transform) {
            x=Ref->Transform[0]+Ref->Transform[1]*X[n]+Ref->Transform[2]*Y[n];
            y=Ref->Transform[3]+Ref->Transform[4]*X[n]+Ref->Transform[5]*Y[n];
         } else if (Ref->GCPTransform) {
            GDALGCPTransform(Ref->GCPTransform,FALSE,1,&x,&y,&z,&ok);
         } else if (Ref->TPSTransform) {
            GDALTPSTransform(Ref->TPSTransform,FALSE,1,&x,&y,&z,&ok);
         } else if (Ref->RPCTransform) {
            GDALRPCTransform(Ref->RPCTransform,FALSE,1,&x,&y,&z,&ok);
         }
      }
      
      // Transform to latlon
      if (Ref->Function) {
         if (!OCTTransform(Ref->Function,1,&x,&y,NULL)) {
            Lon[n]=-999.0;
            Lat[n]=-999.0;
            continue;
         }
      }

      Lon[n]=x;
      Lat[n]=y;
      
      if (Ref->RotTransform) 
         GeoRef_WKTUnRotate(Ref->RotTransform,&Lat[n],&Lon[n]);
   }
#else
   App_Log(ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
#endif

   return(0);
}

int GeoRef_LL2XY_W(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

#ifdef HAVE_GDAL
   double x,y,z=0.0;
   int    n,ok;

   for(n=0;n<Nb;n++) {
      if (Ref->RotTransform) 
         GeoRef_WKTRotate(Ref->RotTransform,&Lat[n],&Lon[n]);

      // Longitude from -180 to 180
      CLAMPLON(Lon[n]);

      x=Lon[n];
      y=Lat[n];

      // Transform from latlon 
      if (Ref->InvFunction) {
         if (!OCTTransform(Ref->InvFunction,1,&x,&y,NULL)) {
            X[n]=-1.0;
            Y[n]=-1.0;
            continue;
         }
      }

      // Transform from georeferenced coordinates
      X[n]=x;
      Y[n]=y;
      if (Ref->Options.Transform) {
         if (Ref->InvTransform) {
            X[n]=Ref->InvTransform[0]+Ref->InvTransform[1]*x+Ref->InvTransform[2]*y;
            Y[n]=Ref->InvTransform[3]+Ref->InvTransform[4]*x+Ref->InvTransform[5]*y;
         } else if (Ref->GCPTransform) {
            GDALGCPTransform(Ref->GCPTransform,TRUE,1,&X[n],&Y[n],&z,&ok);
         } else if (Ref->TPSTransform) {
            GDALTPSTransform(Ref->TPSTransform,TRUE,1,&X[n],&Y[n],&z,&ok);
         } else if (Ref->RPCTransform) {
            GDALRPCTransform(Ref->RPCTransform,TRUE,1,&X[n],&Y[n],&z,&ok);
         }
      }   
   }
#else
   App_Log(ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
#endif
   return(0);
}
