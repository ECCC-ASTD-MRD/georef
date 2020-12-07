/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Fonctions et definitions relatives aux fichiers standards et rmnlib
 * Fichier      : GeoRef_Type_R.h
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
 * Nom          : <GeoRef_RDRHeight>
 * Creation     : Decembre 2009 J.P. Gauthier - CMC/CMOE
 *
 * But          : Calculer la hauteur en MAGL d<une coordonnee RADAR.
 *
 * Parametres    :
 *   <Ref>      : Pointeur sur la reference geographique
 *   <ZRef>      : Pointeur sur la reference verticale
 *   <Azimuth>   : coordonnee en X dans la projection/grille
 *   <Bin>       : coordonnee en Y dans la projection/grille
 *   <Sweep>     : coordonnee en Z dans la projection/grille
 *
 * Retour       : Hauteur
 *
 * Remarques   :
 *
 *---------------------------------------------------------------------------------------------------------------
*/
//TODO: vgrid of ZRef ?
double GeoRef_RDRHeight(TGeoRef *Ref,TZRef *ZRef,double Azimuth,double Bin,double Sweep) {

   if (Bin>=0 && Bin<Ref->R && Sweep>=0 && Sweep<ZRef->LevelNb) {
      return(Ref->Loc.Elev+sin(DEG2RAD(ZRef->Levels[(int)Sweep]))*((int)Bin*Ref->ResR));
   } else {
      return(0);
   }
}

/*--------------------------------------------------------------------------------------------------------------
 * Nom          : <GeoRef_RDRCreate>
 * Creation     : Avril 2006 J.P. Gauthier - CMC/CMOE
 *
 * But          : Definir le referetiel de type Radar
 *
 * Parametres   :
 *    <Lat>     : Latitude du centre
 *    <Lon>     : Longitude du centre
 *    <Height>  : Altitude du centre
 *    <NBin>    : Nombre de bin
 *    <ResR>    : Resolution en distance
 *    <ResA>    : Resolution en azimuth
 *
 * Retour       :
 *
 * Remarques    :
 *
 *---------------------------------------------------------------------------------------------------------------
*/
TGeoRef* GeoRef_CreateR(double Lat,double Lon,double Height,int R,double ResR,double ResA) {

   TGeoRef *ref;

   ref=GeoRef_New();

   ref->GRTYP[0]='R';
   ref->Loc.Lat=Lat;
   ref->Loc.Lon=Lon;
   ref->Loc.Elev=Height;
   ref->R=R;
   ref->ResR=ResR;
   ref->ResA=ResA;

   GeoRef_Size(ref,0,0,360/ResA,R-1,0);

   ref->Height=GeoRef_RDRHeight;

   return(ref);
}

int GeoRef_XY2LL_R(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {

   TCoord loc0;
   double x,y,d;
   int    n;

   for(n=0;n<Nb;n++) {

      loc0.Lat=DEG2RAD(Ref->Loc.Lat);
      loc0.Lon=DEG2RAD(Ref->Loc.Lon);

      x=X[n]*Ref->ResA;
      y=Y[n]*Ref->ResR;

      x=DEG2RAD(x);
      d=M2RAD(y*Ref->CTH);

      if (Ref->Options.Transform) {
         Lat[n]=asin(sin(loc0.Lat)*cos(d)+cos(loc0.Lat)*sin(d)*cos(x));
         Lon[n]=fmod(loc0.Lon+(atan2(sin(x)*sin(d)*cos(loc0.Lat),cos(d)-sin(loc0.Lat)*sin(*Lat)))+M_PI,M_2PI)-M_PI;
         Lat[n]=RAD2DEG(*Lat);
         Lon[n]=RAD2DEG(*Lon);
      } else {
         Lat[n]=d;
         Lon[n]=x;
      }
   }

   return(0);
}

int GeoRef_LL2XY_R(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   TCoord  loc0;
   double  x,d,lat,lon;
   int     n;

   for(n=0;n<Nb;n++) {
      loc0.Lat=DEG2RAD(Ref->Loc.Lat);
      loc0.Lon=DEG2RAD(Ref->Loc.Lon);
      lat=DEG2RAD(Lat[n]);
      lon=DEG2RAD(Lon[n]);

      d=fabs(DIST(0.0,loc0.Lat,loc0.Lon,lat,lon));
      x=-RAD2DEG(COURSE(loc0.Lat,loc0.Lon,lat,lon));
      X[n]=x<0.0?x+360.0:x;
      Y[n]=d/Ref->CTH;

      if (Ref->Options.Transform) {
         X[n]/=Ref->ResA;
         Y[n]/=Ref->ResR;
      }
   }

   return(0);
}
