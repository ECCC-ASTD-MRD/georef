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
 * Description  : Fonctions de manipulations de projections aux standard WKT.
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
#include "Vertex.h"

double   GeoRef_WKTDistance(TGeoRef *Ref,double X0,double Y0,double X1, double Y1);
int      GeoRef_WKTValue(TGeoRef *Ref,TDef *Def,TDef_InterpR Interp,int C,double X,double Y,double Z,double *Length,double *ThetaXY);
int      GeoRef_WKTProject(TGeoRef *Ref,double X,double Y,double *Lat,double *Lon,int Extrap,int Transform);
int      GeoRef_WKTUnProject(TGeoRef *Ref,double *X,double *Y,double Lat,double Lon,int Extrap,int Transform);

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_WKTDistance>
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
double GeoRef_WKTDistance(TGeoRef *Ref,double X0,double Y0,double X1, double Y1) {

#ifdef HAVE_GDAL
   double i[2],j[2],lat[2],lon[2],u;
   char *unit,geo;
   
   X0+=Ref->X0;
   X1+=Ref->X0;
   Y0+=Ref->Y0;
   Y1+=Ref->Y0;
   
   // Check for unit type 
   geo=0;
   if (Ref->GRTYP[0]=='Z' || Ref->GRTYP[1]=='Z') {
      geo=1;
   } else {
      if (Ref->Spatial) {
         u=OSRGetLinearUnits(Ref->Spatial,&unit);
         geo=(unit[0]!='M' && unit[0]!='m');
      }
   }

   if (geo) {
      GeoRef_WKTProject(Ref,X0,Y0,&lat[0],&lon[0],1,1);
      GeoRef_WKTProject(Ref,X1,Y1,&lat[1],&lon[1],1,1);
      return(DIST(0.0,DEG2RAD(lat[0]),DEG2RAD(lon[0]),DEG2RAD(lat[1]),DEG2RAD(lon[1])));
   } else {
      if (Ref->Transform) {
         i[0]=Ref->Transform[0]+Ref->Transform[1]*X0+Ref->Transform[2]*Y0;
         j[0]=Ref->Transform[3]+Ref->Transform[4]*X0+Ref->Transform[5]*Y0;
         i[1]=Ref->Transform[0]+Ref->Transform[1]*X1+Ref->Transform[2]*Y1;
         j[1]=Ref->Transform[3]+Ref->Transform[4]*X1+Ref->Transform[5]*Y1;
      } else {
         i[0]=X0;
         j[0]=Y0;
         i[1]=X1;
         j[1]=Y1;
      }
      return(hypot(j[1]-j[0],i[1]-i[0])*u);
   }
#else
   App_Log(ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
   return(0.0);
#endif
}

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

#ifdef HAVE_GDAL
   double d,dx,dy,x,y,z=0.0;
   int    sx,sy,s,ok,gidx;

   d=1.0;

   if( !Extrap && (X>(Ref->X1+d) || Y>(Ref->Y1+d) || X<(Ref->X0-d) || Y<(Ref->Y0-d)) ) {
      *Lon=-999.0;
      *Lat=-999.0;
      return(0);
   }

   // GRTYP cell are corner defined 
   if (Ref->Type&GRID_CORNER) {
      X+=0.5;
      Y+=0.5;
   }

   // Because some grids are defined as WZ and others as ZW, this makes sure we catch the other letter
   gidx = Ref->GRTYP[0]=='W' ? 1 : 0;

   // In case of non-uniform grid, figure out where in the position vector we are
   if (Ref->GRTYP[gidx]=='Z') {
      if (Ref->AX && Ref->AY) {
         // X
         if( X < Ref->X0 )      { sx=Ref->X0; X=Ref->AX[sx]-(Ref->AX[sx]-Ref->AX[sx+1])*(X-sx); }
         else if( X > Ref->X1 ) { sx=Ref->X1; X=Ref->AX[sx]+(Ref->AX[sx]-Ref->AX[sx-1])*(X-sx); }
         else                    { sx=floor(X); X=sx==X?Ref->AX[sx]:ILIN(Ref->AX[sx],Ref->AX[sx+1],X-sx); }

         // Y
         s=Ref->NX;
         if( Y < Ref->Y0 )      { sy=Ref->Y0; Y=Ref->AY[sy*s]-(Ref->AY[sy*s]-Ref->AY[(sy+1)*s])*(Y-sy); }
         else if( Y > Ref->Y1 ) { sy=Ref->Y1; Y=Ref->AY[sy*s]+(Ref->AY[sy*s]-Ref->AY[(sy-1)*s])*(Y-sy); }
         else                    { sy=floor(Y); Y=sy==Y?Ref->AY[sy*s]:ILIN(Ref->AY[sy*s],Ref->AY[(sy+1)*s],Y-sy); }
      }
   } else if (Ref->GRTYP[gidx]=='X' || Ref->GRTYP[gidx]=='Y') {
      if (Ref->AX && Ref->AY) {
         sx=floor(X);sx=CLAMP(sx,Ref->X0,Ref->X1);
         sy=floor(Y);sy=CLAMP(sy,Ref->Y0,Ref->Y1);
         dx=X-sx;;
         dy=Y-sy;

         s=sy*Ref->NX+sx;
         X=Ref->AX[s];
         Y=Ref->AY[s];

         if (++sx<=Ref->X1) {
            s=sy*Ref->NX+sx;
            X+=(Ref->AX[s]-X)*dx;
         }

         if (++sy<=Ref->Y1) {
            s=sy*Ref->NX+(sx-1);
            Y+=(Ref->AY[s]-Y)*dy;
         }
      }
   }

   // Transform the point into georeferenced coordinates 
   x=X;
   y=Y;
   if (Transform) {
      if (Ref->Transform) {
         x=Ref->Transform[0]+Ref->Transform[1]*X+Ref->Transform[2]*Y;
         y=Ref->Transform[3]+Ref->Transform[4]*X+Ref->Transform[5]*Y;
      } else if (Ref->GCPTransform) {
         GDALGCPTransform(Ref->GCPTransform,FALSE,1,&x,&y,&z,&ok);
      } else if (Ref->TPSTransform) {
         GDALGCPTransform(Ref->TPSTransform,FALSE,1,&x,&y,&z,&ok);
      } else if (Ref->RPCTransform) {
         GDALGCPTransform(Ref->RPCTransform,FALSE,1,&x,&y,&z,&ok);
      }
   }
   
   // Transform to latlon
   if (Ref->Function) {
      if (!OCTTransform(Ref->Function,1,&x,&y,NULL)) {
         *Lon=-999.0;
         *Lat=-999.0;
         return(0);
      }
   }

   *Lon=x;
   *Lat=y;
   
   if (Ref->RotTransform) 
      GeoRef_WKTUnRotate(Ref->RotTransform,Lat,Lon);
   
   return(1);
#else
   App_Log(ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
   return(0);
#endif
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

#ifdef HAVE_GDAL
   double x,y,z=0.0,d=1e32,sd;
   int    n,nd,s,dx,dy,ok,idx,idxs[8],gidx;
   double dists[8];
   Vect2d pts[4],pt;
   Vect3d b;
   
   if (Ref->RotTransform) 
      GeoRef_WKTRotate(Ref->RotTransform,&Lat,&Lon);

   if (Lat<=90.0 && Lat>=-90.0 && Lon!=-999.0) {

      // Longitude from -180 to 180
      Lon=Lon>180?Lon-360:Lon;

      x=Lon;
      y=Lat;

      // Transform from latlon 
      if (Ref->InvFunction) {
         if (!OCTTransform(Ref->InvFunction,1,&x,&y,NULL)) {
            *X=-1.0;
            *Y=-1.0;
            return(0);
         }
      }

      // Transform from georeferenced coordinates
      *X=x;
      *Y=y;
      if (Transform) {
         if (Ref->InvTransform) {
            *X=Ref->InvTransform[0]+Ref->InvTransform[1]*x+Ref->InvTransform[2]*y;
            *Y=Ref->InvTransform[3]+Ref->InvTransform[4]*x+Ref->InvTransform[5]*y;
         } else if (Ref->GCPTransform) {
            GDALGCPTransform(Ref->GCPTransform,TRUE,1,X,Y,&z,&ok);
         } else if (Ref->TPSTransform) {
            GDALTPSTransform(Ref->TPSTransform,TRUE,1,X,Y,&z,&ok);
         } else if (Ref->RPCTransform) {
            GDALRPCTransform(Ref->RPCTransform,TRUE,1,X,Y,&z,&ok);
         }

         // Because some grids are defined as WZ and others as ZW, this makes sure we catch the other letter
         gidx = Ref->GRTYP[0]=='W' ? 1 : 0;

         // In case of non-uniform grid, figure out where in the position vector we are 
         if (Ref->GRTYP[gidx]=='Z') {
            if (Ref->AX && Ref->AY) {
               s=Ref->X0;
               // Check if vector is increasing
               if (Ref->AX[s]<Ref->AX[s+1]) {
                  while(s<=Ref->X1 && *X>Ref->AX[s]) s++;
               } else {
                  while(s<=Ref->X1 && *X<Ref->AX[s]) s++;
               }
               if (s>Ref->X0) {
                  // We're in so interpolate postion
                  if (s<=Ref->X1) {
                     *X=(*X-Ref->AX[s-1])/(Ref->AX[s]-Ref->AX[s-1])+s-1;
                  } else {
                     *X=(*X-Ref->AX[Ref->X1])/(Ref->AX[Ref->X1]-Ref->AX[Ref->X1-1])+s-1;
                  }
               } else {
                  // We're out so extrapolate position
                  *X=Ref->X0+(*X-Ref->AX[0])/(Ref->AX[1]-Ref->AX[0]);
               }

               s=Ref->Y0;dx=Ref->NX;
               // Check if vector is increasing
               if (Ref->AY[s*Ref->NX]<Ref->AY[(s+1)*Ref->NX]) {
                  while(s<=Ref->Y1 && *Y>Ref->AY[s*Ref->NX]) s++;
               } else {
                  while(s<=Ref->Y1 && *Y<Ref->AY[s*Ref->NX]) s++;
               }
               if (s>Ref->Y0) {
                  // We're in so interpolate postion
                  if (s<=Ref->Y1) {
                     *Y=(*Y-Ref->AY[(s-1)*Ref->NX])/(Ref->AY[s*Ref->NX]-Ref->AY[(s-1)*Ref->NX])+s-1;
                  } else {
                     *Y=(*Y-Ref->AY[Ref->Y1*Ref->NX])/(Ref->AY[Ref->Y1*Ref->NX]-Ref->AY[(Ref->Y1-1)*Ref->NX])+s-1;
                  }
               } else {
                  // We're out so extrapolate position
                  *Y=Ref->Y0+(*Y-Ref->AY[0])/(Ref->AY[Ref->NX]-Ref->AY[0]);
               }
            }
         } else if (Ref->GRTYP[gidx]=='Y') {
            // Get nearest point
            if (GeoRef_Nearest(Ref,Lon,Lat,&idx,dists,1)) {
               if (dists[0]<1.0) {
                  *Y=(int)(idx/Ref->NX);
                  *X=idx-(*Y)*Ref->NX;
                  return(TRUE);
               }
            }
         } else if (Ref->GRTYP[gidx]=='X') {
            // Get nearest points
            if ((nd=GeoRef_Nearest(Ref,Lon,Lat,idxs,dists,8))) {
               
               pt[0]=Lon;
               pt[1]=Lat;

               // Find which cell includes coordinates
               for(n=0;n<nd;n++) {
                  idx=idxs[n];

                  // Find within which quad
                  dx=-1;dy=-1;
                  if (!GeoRef_WithinCell(Ref,pt,pts,idx-Ref->NX-1,idx-1,idx,idx-Ref->NX)) {
                 
                     dx=0;dy=-1;
                     if (!GeoRef_WithinCell(Ref,pt,pts,idx-Ref->NX,idx,idx+1,idx-Ref->NX+1)) {
                        
                        dx=-1;dy=0;
                        if (!GeoRef_WithinCell(Ref,pt,pts,idx-1,idx+Ref->NX-1,idx+Ref->NX,idx)) {
                     
                           dx=0;dy=0;
                           if (!GeoRef_WithinCell(Ref,pt,pts,idx,idx+Ref->NX,idx+Ref->NX+1,idx+1)) {
                              idx=-1;
                           }
                        }
                     }
                  }
                  
                  // If found, exit loop
                  if (idx!=-1) {
                     break;
                  }
               }
               
               if (idx!=-1) {
                  // Map coordinates to grid
                  Vertex_Map(pts,X,Y,Lon,Lat);
                  
                  if (!ISNAN(*X) && !ISNAN(*Y)) {
                     y=idx/Ref->NX;
                     x=idx-y*Ref->NX;
                     *Y+=y+dy;
                     *X+=x+dx; 
                  } else {
                     *X=-1,0;
                     *Y=-1.0;
                     return(FALSE);
                  }
               }
            }
         }

         // Check the grid limits
         d=1.0;
         if (*X>(Ref->X1+d) || *Y>(Ref->Y1+d) || *X<(Ref->X0-d) || *Y<(Ref->Y0-d)) {
            if (!Extrap) {
               *X=-1.0;
               *Y=-1.0;
            }
            return(0);
         }
      }
      // Grid cell are corner defined 
      if (Ref->Type&GRID_CORNER) {
         *X-=0.5;
         *Y-=0.5;
      }
      
   } else {
      *X=-1.0;
      *Y=-1.0;
      return(0);
   }
   return(1);
#else
   App_Log(ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
   return(0);
#endif
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
int GeoRef_WKTSet(TGeoRef *Ref,char *String,double *Transform,double *InvTransform,OGRSpatialReferenceH Spatial) {

#ifdef HAVE_GDAL
   static OGRSpatialReferenceH llref=NULL;
   char                      *string=NULL;

   if (String && String[0]!='\0') {
      string=strdup(String);
      strtrim(string,' ');
   }

   GeoRef_Clear(Ref,0);

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
        return(0);
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
      return(0);
   }

   Ref->Project=GeoRef_WKTProject;
   Ref->UnProject=GeoRef_WKTUnProject;
   Ref->Value=(TGeoRef_Value*)GeoRef_WKTValue;
   Ref->Distance=GeoRef_WKTDistance;
   Ref->Height=NULL;

   return(1);
#else
   App_Log(ERROR,"Function %s is not available, needs to be built with GDAL\n",__func__);
   return(0);
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
TGeoRef *GeoRef_WKTCreate(int ni,int nj,char *grtyp,int ig1,int ig2,int ig3,int ig4,char *String,double *Transform,double *InvTransform,OGRSpatialReferenceH Spatial) {

   TGeoRef *ref;

   ref=GeoRef_New();
   GeoRef_Size(ref,0,0,ni>0?ni-1:0,nj>0?nj-1:0,0);
   if (!GeoRef_WKTSet(ref,String,Transform,InvTransform,Spatial)) {
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
