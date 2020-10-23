/*=========================================================
 * Environnement Canada
 * Centre Meteorologique Canadien
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Lecture et traitements de fichiers raster
 * Fichier      : GeoRef.c
 * Creation     : Mars 2005 - J.P. Gauthier
 *
 * Description  : Fonctions de manipulations de projections.
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
#define GEOREF_BUILD

#include <pthread.h>
#include "GeoRef.h"
#include "App.h"
#include "Def.h"
#include "RPN.h"
#include "List.h"

static TList          *GeoRef_List=NULL;                                                                                ///< Global list of known geo references
static pthread_mutex_t GeoRef_Mutex=PTHREAD_MUTEX_INITIALIZER;                                                          ///< Thread lock on geo reference access
static TGeoOptions     GeoRef_Options= { IR_CUBIC, ER_UNDEF, 0.0, 0, TRUE, FALSE, FALSE, 16, TRUE, FALSE, 10.0, 0.0 };  ///< Default options

/**----------------------------------------------------------------------------
 * @brief  Apply thread lock on GeoRef access
 * @author Jean-Philippe Gauthier
 * @date   February 2008
*/
void GeoRef_Lock() {
   pthread_mutex_lock(&GeoRef_Mutex);
}
/**----------------------------------------------------------------------------
 * @brief  Remove thread lock on GeoRef access
 * @author Jean-Philippe Gauthier
 * @date   February 2008
*/
void GeoRef_Unlock() {
   pthread_mutex_unlock(&GeoRef_Mutex);
}

/**----------------------------------------------------------------------------
 * @brief  Re-initialise reprojection buffer
 * @author Jean-Philippe Gauthier
 * @date   February 2008
 *    @param[in]  Scan     Reprojection buffer
*/
void GeoScan_Clear(TGeoScan *Scan) {

   if (Scan) {
      if (Scan->X) free(Scan->X);
      if (Scan->Y) free(Scan->Y);
      if (Scan->V) free(Scan->V);
      if (Scan->D) free(Scan->D);

      Scan->X=Scan->Y=NULL;
      Scan->V=NULL;
      Scan->D=NULL;
      Scan->N=Scan->S=Scan->DX=Scan->DY=0;
   }
}

/**----------------------------------------------------------------------------
 * @brief  Initialiser la structure App
 * @author Jean-Philippe Gauthier
 * @date   Janvier 2017
 *    @param[in]  Type     App type (APP_MASTER=single independent process, APP_THREAD=threaded co-process)
 *    @param[in]  Name     Application name
 *    @param[in]  Version  Application version
 *    @param[in]  Desc     Application description
 *    @param[in]  Stamp    TimeStamp
 *
 *    @return              Parametres de l'application initialisee
*/void GeoScan_Init(TGeoScan *Scan) {

   if (Scan) {
      Scan->X=Scan->Y=NULL;
      Scan->V=NULL;
      Scan->D=NULL;
      Scan->N=Scan->S=Scan->DX=Scan->DY=0;
   }
}

/**----------------------------------------------------------------------------
 * @brief  Reproject a stream of coordinates and extracte values
 * @author Jean-Philippe Gauthier
 * @date   February 2008
 *    @param[in]  Scan     Reprojection buffer
 *    @param[in]  ToRef    Destination geo reference
 *    @param[in]  ToDef    Destination data definition
 *    @param[in]  FromRef  Source geo reference   
 *    @param[in]  FromDef  Destination data definition  
 *    @param[in]  X0       Lower x limit
 *    @param[in]  Y0       Lower y limiy
 *    @param[in]  X1       Higher x limit
 *    @param[in]  Y1       Higher y limit
 *    @param[in]  Dim      Grid cell dimension (1=point, 2=area)
 *    @param[in]  Degree   Interpolation degree 
 *
 *    @return              Size of results
*/
int _GeoScan_Get(TGeoScan *Scan,TGeoRef *ToRef,TDef *ToDef,TGeoRef *FromRef,TDef *FromDef,int X0,int Y0,int X1,int Y1,int Dim,TDef_InterpR Degree) {

   register int idx,x,y,n=0;
   int          d=0,sz,dd;
   double       x0,y0,v;
   
   if (!Scan || !ToRef || !FromRef) {
      return(0);
   }

   // Check limits
   X0=fmax(X0,FromRef->X0);
   Y0=fmax(Y0,FromRef->Y0);
   X1=fmin(X1,FromRef->X1);
   Y1=fmin(Y1,FromRef->Y1);

   // Adjust scan buffer sizes
   Scan->DX=X1-X0+1;
   Scan->DY=Y1-Y0+1;
   dd=(Scan->DX+1)*(Scan->DY+1);
   sz=Scan->DX*Scan->DY;

   if (Scan->S<sz) {
      if (!(Scan->X=(double*)realloc(Scan->X,dd*sizeof(double))))
         return(0);
      if (!(Scan->Y=(double*)realloc(Scan->Y,dd*sizeof(double))))
         return(0);
      if (!(Scan->V=(unsigned int*)realloc(Scan->V,sz*sizeof(unsigned int))))
         return(0);
      if (!(Scan->D=(float*)realloc(Scan->D,sz*sizeof(float))))
         return(0);
      Scan->S=sz;
   }

   dd=Dim-1;
   Scan->N=0;

   for(y=Y0;y<=Y1+dd;y++) {
      idx=(y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
      for(x=X0;x<=X1+dd;x++,idx++,n++) {
         if (x<=X1 && y<=Y1) {
            Scan->V[Scan->N++]=idx;
         }

         x0=dd?x-0.5:x;
         y0=dd?y-0.5:y;
         
         FromRef->Project(FromRef,x0,y0,&Scan->X[n],&Scan->Y[n],0,1);


         if (FromRef->Transform) {
            Scan->X[n]=FromRef->Transform[0]+FromRef->Transform[1]*x0+FromRef->Transform[2]*y0;
            Scan->Y[n]=FromRef->Transform[3]+FromRef->Transform[4]*x0+FromRef->Transform[5]*y0;
         } else {
            Scan->X[n]=x0;
            Scan->Y[n]=y0;
         }
      }
   }

   // WKT grid type
   if (FromRef->GRTYP[0]=='W') {
#ifdef HAVE_GDAL
      for(y=Y0;y<=Y1+dd;y++) {
         idx=(y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
         for(x=X0;x<=X1+dd;x++,idx++,n++) {
            if (x<=X1 && y<=Y1) {
               Scan->V[Scan->N++]=idx;
            }

            x0=dd?x-0.5:x;
            y0=dd?y-0.5:y;
            if (FromRef->Transform) {
               Scan->X[n]=FromRef->Transform[0]+FromRef->Transform[1]*x0+FromRef->Transform[2]*y0;
               Scan->Y[n]=FromRef->Transform[3]+FromRef->Transform[4]*x0+FromRef->Transform[5]*y0;
            } else {
               Scan->X[n]=x0;
               Scan->Y[n]=y0;
            }
         }
      }

      if (FromRef->Function) {
         OCTTransform(FromRef->Function,n,Scan->X,Scan->Y,NULL);
      }
#endif
      d=dd?2:1;
      sz=8;

   // Y GRTYP type
   } else if (FromRef->GRTYP[0]=='Y') {
      for(y=Y0;y<=Y1;y++) {
         idx=(y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
         for(x=X0;x<=X1;x++,idx++,n++) {
            if (x<=X1 && y<=Y1) {
               Scan->V[Scan->N++]=idx;
            }
            ((float*)Scan->X)[n]=FromRef->AX[idx];
            ((float*)Scan->Y)[n]=FromRef->AY[idx];
         }
      }
      d=1;
      sz=4;

   // Other RPN grids
   } else {
#ifdef HAVE_RMN
      for(y=Y0;y<=Y1+dd;y++) {
         idx=(y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
         for(x=X0;x<=X1+dd;x++,idx++,n++) {
            if (x<=X1 && y<=Y1) {
               Scan->V[Scan->N++]=idx;
            }
            Scan->X[n]=dd?x+0.5:x+1.0;
            Scan->Y[n]=dd?y+0.5:y+1.0;
         }
      }
      GeoRef_XY2LL(REFGET(FromRef),Scan->Y,Scan->X,Scan->X,Scan->Y,n);

      d=dd?2:1;
      sz=4;
#else
      App_Log(ERROR,"%s: RMNLIB support not included\n",__func__);
#endif
   }

   // Project to destination grid
   for(x=n-1;x>=0;x--) {
      if (sz==4) {
         x0=(double)((float*)Scan->X)[x];
         y0=(double)((float*)Scan->Y)[x];
      } else {
         x0=Scan->X[x];
         y0=Scan->Y[x];
      }
      
      if (ToDef) {
         Scan->D[x]=ToDef->NoData;
      }

      // If we're inside
      if (ToRef->UnProject(ToRef,&Scan->X[x],&Scan->Y[x],y0,x0,0,1) && ToDef) {
         ToRef->Value(ToRef,ToDef,Degree,0,Scan->X[x],Scan->Y[x],0,&v,NULL);
         Scan->D[x]=v;
      }
   }
   return(d);
}

int GeoScan_Get(TGeoScan *Scan,TGeoRef *ToRef,TDef *ToDef,TGeoRef *FromRef,TDef *FromDef,int X0,int Y0,int X1,int Y1,int Dim,TDef_InterpR Degree) {

   register int idx,x,y,n=0;
   int          d=0,sz,dd;
   double       x0,y0,v;
   int          ix, iy;
   
   if (!Scan || !ToRef || !FromRef) {
      return(0);
   }

   // Check limits
   X0=fmax(X0,FromRef->X0);
   Y0=fmax(Y0,FromRef->Y0);
   X1=fmin(X1,FromRef->X1);
   Y1=fmin(Y1,FromRef->Y1);

   // Adjust scan buffer sizes
   Scan->DX=X1-X0+1;
   Scan->DY=Y1-Y0+1;
   dd=(Scan->DX+1)*(Scan->DY+1);
   sz=Scan->DX*Scan->DY;

   if (Scan->S<sz) {
      if (!(Scan->X=(double*)realloc(Scan->X,dd*sizeof(double))))
         return(0);
      if (!(Scan->Y=(double*)realloc(Scan->Y,dd*sizeof(double))))
         return(0);
      if (!(Scan->V=(unsigned int*)realloc(Scan->V,sz*sizeof(unsigned int))))
         return(0);
      if (!(Scan->D=(float*)realloc(Scan->D,sz*sizeof(float))))
         return(0);
      Scan->S=sz;
   }

   dd=Dim-1;
   Scan->N=0;

   // WKT grid type
   if (FromRef->GRTYP[0]=='W') {
#ifdef HAVE_GDAL
      for(y=Y0;y<=Y1+dd;y++) {
         idx=(y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
         for(x=X0;x<=X1+dd;x++,idx++,n++) {
            if (x<=X1 && y<=Y1) {
               Scan->V[Scan->N++]=idx;
            }

            x0=dd?x-0.5:x;
            y0=dd?y-0.5:y;
            if (FromRef->Transform) {
               Scan->X[n]=FromRef->Transform[0]+FromRef->Transform[1]*x0+FromRef->Transform[2]*y0;
               Scan->Y[n]=FromRef->Transform[3]+FromRef->Transform[4]*x0+FromRef->Transform[5]*y0;
            } else {
               Scan->X[n]=x0;
               Scan->Y[n]=y0;
            }
         }
      }

      if (FromRef->Function) {
         OCTTransform(FromRef->Function,n,Scan->X,Scan->Y,NULL);
      }
#endif
      d=dd?2:1;
      sz=8;

   // Y Grid type
   } else if (FromRef->GRTYP[0]=='Y') {
      for(y=Y0;y<=Y1;y++) {
         idx=(y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
         for(x=X0;x<=X1;x++,idx++,n++) {
            if (x<=X1 && y<=Y1) {
               Scan->V[Scan->N++]=idx;
            }
            Scan->X[n]=FromRef->AX[idx];
            Scan->Y[n]=FromRef->AY[idx];
         }
      }
      d=1;
      sz=4;

   // Other RPN grids
   } else {
#ifdef HAVE_RMN
      for(y=Y0;y<=Y1+dd;y++) {
         idx=(y-FromRef->Y0)*FromDef->NI+(X0-FromRef->X0);
         for(x=X0;x<=X1+dd;x++,idx++,n++) {
            if (x<=X1 && y<=Y1) {
               Scan->V[Scan->N++]=idx;
            }
            Scan->X[n]=dd?x+0.5:x+1.0;
            Scan->Y[n]=dd?y+0.5:y+1.0;
         }
      }
      GeoRef_XY2LL(REFGET(FromRef),Scan->Y,Scan->X,Scan->X,Scan->Y,n);

      d=dd?2:1;
      sz=4;
#else
      App_Log(ERROR,"%s: RMNLIB support not included\n",__func__);
#endif
   }

   // Project to destination grid
   if (ToRef->GRTYP[0]=='W' || ToRef->GRTYP[0]=='M') {
#ifdef HAVE_GDAL
      for(x=n-1;x>=0;x--) {
         x0=Scan->X[x];
         y0=Scan->Y[x];

         if (ToDef) {
            Scan->D[x]=ToDef->NoData;
         }

         if (ToRef->UnProject(ToRef,&Scan->X[x],&Scan->Y[x],y0,x0,0,1)) {
            if (ToDef) {
              ToRef->Value(ToRef,ToDef,Degree,0,Scan->X[x],Scan->Y[x],0,&v,NULL);
              Scan->D[x]=v;
            }
         }
      }

/*
         if (ToRef->Function)
            OCTTransform(ToRef->InvFunction,n,Scan->X,Scan->Y,NULL);

         if (ToRef->InvTransform) {
            for(x=0;x<n;x++) {
               x0=ToRef->InvTransform[0]+ToRef->InvTransform[1]*Scan->X[x]+ToRef->InvTransform[2]*Scan->Y[x];
               y0=ToRef->InvTransform[3]+ToRef->InvTransform[4]*Scan->X[x]+ToRef->InvTransform[5]*Scan->Y[x];
               Scan->X[x]=x0;
               Scan->Y[x]=y0;
            }
         }
*/
#endif
   } else {
#ifdef HAVE_RMN
      GeoRef_LL2XY(REFGET(ToRef),Scan->X,Scan->Y,Scan->Y,Scan->X,n);
//EZFIX
      // If we have the data of source and they're float, get it's values right now
      if (ToDef && ToDef->Type==TD_Float32) {
         if (Degree)
            ToRef->Options.InterpDegree=Degree;
         
         GeoRef_XYVal(REFGET(ToRef),Scan->D,(float*)ToDef->Mode,Scan->X,Scan->Y,n);         
      }

      // Cast back to double (Start from end since type is double, not to overlap values
      for(x=n-1;x>=0;x--) {
         Scan->X[x]=Scan->X[x]-1.0;
         Scan->Y[x]=Scan->Y[x]-1.0;

         if (ToDef) {
            ix = lrint(Scan->X[x]);
            iy = lrint(Scan->Y[x]);
            idx=FIDX2D(ToDef,ix,iy);
            
            if (!FIN2D(ToDef,ix,iy) || (ToDef->Mask && !ToDef->Mask[idx])) {
               // If we're outside, set to nodata
               Scan->D[x]=ToDef->NoData;
            } else if (ToDef->Type<TD_Float32) {
               // Otherwise, set nearest data if not floats
               Def_GetMod(ToDef,idx,Scan->D[x]);
            }
         }
      }
#else
      App_Log(ERROR,"%s: RMNLIB support not included\n",__func__);
#endif
   }
   return(d);
}

/**----------------------------------------------------------------------------
 * @brief  Initialise the georeference grid limits
 * @author Jean-Philippe Gauthier
 * @date   July 2005
 *    @param[in]  Ref   Pointer to geo reference
 *    @param[in]  X0    X lower limit
 *    @param[in]  Y0    Y lower limit
 *    @param[in]  X1    X higher limit
 *    @param[in]  Y0    Y higher limit
 *    @param[in]  BD    Border width
 */
void GeoRef_Size(TGeoRef *Ref,int X0,int Y0,int X1,int Y1,int BD) {

   Ref->X0=X0;
   Ref->X1=X1;
   Ref->Y0=Y0;
   Ref->Y1=Y1;
   Ref->BD=BD;
   Ref->NX=X1-X0+1;
   Ref->NY=Y1-Y0+1;
}

/**----------------------------------------------------------------------------
 * @brief  Free the geo reference resources
 * @author Jean-Philippe Gauthier
 * @date   July 2005
 *    @param[in]  Ref   Pointer to geo reference
 *
 *    @return           Freed code (0=not freed, 1=freed)
 * 
 */
int GeoRef_Free(TGeoRef *Ref) {

  if (!Ref)
      return(0);

   if (__sync_sub_and_fetch(&Ref->NRef,1)!=0) {
      return(0);
   }

   if (Ref->RefFrom) {
//      GeoRef_Free(Ref->RefFrom);
      __sync_sub_and_fetch(&Ref->RefFrom->NRef,1);
   }

   GeoRef_Clear(Ref,1);
   free(Ref);

   // Remove from Georef list
   GeoRef_Lock();
   GeoRef_List=TList_Del(GeoRef_List,(void*)Ref);
   GeoRef_Unlock();

   return(1);
}

/**----------------------------------------------------------------------------
 * @brief  Increment reference count
 * @author Jean-Philippe Gauthier
 * @date   July 2005
 *    @param[in]  Ref   Pointer to geo reference
 *
 *    @return           New reference count
*/
int GeoRef_Incr(TGeoRef *Ref) {

   if (Ref) {
      return(__sync_add_and_fetch(&Ref->NRef,1));
   } else {
      return(0);
   }
}

/**----------------------------------------------------------------------------
 * @brief  Initialiser la structure App
 * @author Jean-Philippe Gauthier
 * @date   July 2005
 *    @param[in]  Ref   Pointer to geo reference
 *    @param[in]  New   Clear the name associated
*/
void GeoRef_Clear(TGeoRef *Ref,int New) {

   int n;
   
   if (Ref) {
      if (Ref->String)       free(Ref->String);       Ref->String=NULL;
      if (Ref->Transform)    free(Ref->Transform);    Ref->Transform=NULL;
      if (Ref->InvTransform) free(Ref->InvTransform); Ref->InvTransform=NULL;
      if (Ref->RotTransform) free(Ref->RotTransform); Ref->RotTransform=NULL;
      if (Ref->Lat)          free(Ref->Lat);          Ref->Lat=NULL;
      if (Ref->Lon)          free(Ref->Lon);          Ref->Lon=NULL;
      if (Ref->Hgt)          free(Ref->Hgt);          Ref->Hgt=NULL;
      if (Ref->Wght)         free(Ref->Wght);         Ref->Wght=NULL;
      if (Ref->Idx)          free(Ref->Idx);          Ref->Idx=NULL; Ref->NIdx=0;
      if (Ref->AX)           free(Ref->AX);           Ref->AX=NULL;
      if (Ref->AY)           free(Ref->AY);           Ref->AY=NULL;
      if (Ref->NCX)          free(Ref->NCX);          Ref->NCX=NULL;
      if (Ref->NCY)          free(Ref->NCY);          Ref->NCY=NULL;
      if (Ref->Subs)         free(Ref->Subs);         Ref->Subs=NULL;

      // Free interpolation sets
      for (n=0;n<Ref->NbSet;n++) {
         GeoRef_SetFree(&Ref->Sets[n]);
      }
      Ref->NbSet=0;

      // Free quadtree index
      if (Ref->QTree) {
         // Check if we actually have a QTree instead of the special kind of index used for more regular grids
         // The hypothesis being that we should at least have children on the head node of any real QTree
         if( Ref->QTree[0].Childs[0] ) {
            QTree_Free(Ref->QTree);
         } else {
            // We have a special kind of index used for more regular grids instead
            for(n=0; n<(GRID_YQTREESIZE+1)*(GRID_YQTREESIZE+1); ++n) {
               QTree_DelData(&Ref->QTree[n]);
            }
            free(Ref->QTree);
         }

         Ref->QTree=NULL;
      }

      memset(&Ref->RPNHead,0x0,sizeof(TRPNHeader));

      if (New) {
         if (Ref->Name)      free(Ref->Name);         Ref->Name=NULL;
      }

#ifdef HAVE_GDAL
      if (Ref->GCPTransform) {
         GDALDestroyGCPTransformer(Ref->GCPTransform);
         Ref->GCPTransform=NULL;
      }
      if (Ref->TPSTransform) {
         GDALDestroyTPSTransformer(Ref->TPSTransform);
         Ref->TPSTransform=NULL;
      }
      if (Ref->RPCTransform) {
         GDALDestroyRPCTransformer(Ref->RPCTransform);
         Ref->RPCTransform=NULL;
      }
      if (Ref->Spatial) {
         OSRDestroySpatialReference(Ref->Spatial);
      }
      
      if (Ref->Function) {
         OCTDestroyCoordinateTransformation(Ref->Function);
         Ref->Function=NULL;
      }
      
      if (Ref->InvFunction) {
         OCTDestroyCoordinateTransformation(Ref->InvFunction);
         Ref->InvFunction=NULL;
      }
#endif

      Ref->RefFrom=NULL;
      Ref->Project=NULL;
      Ref->UnProject=NULL;
      Ref->Value=NULL;
      Ref->Distance=NULL;
      Ref->Height=NULL;
   }
}

/**----------------------------------------------------------------------------
 * @brief  Define qualifying flags for the geo refrence
 * @author Jean-Philippe Gauthier
 * @date   Janvier 2015
 *    @param[in]  Ref   Pointer to geo reference
 */
void GeoRef_Qualify(TGeoRef* __restrict const Ref) {

   TCoord co[2];
   double d[2];
   int    x;

   if (Ref) {
      Ref->Type=GRID_NONE;

      if (Ref->GRTYP[0]!='P' && Ref->GRTYP[0]!='M' && Ref->GRTYP[0]!='V') {
         Ref->Type|=GRID_EZ;
      }

      if (Ref->GRTYP[0]=='M' || Ref->GRTYP[0]=='Y' || Ref->GRTYP[0]=='X' || Ref->GRTYP[0]=='O') {
         Ref->Type|=GRID_SPARSE;
      } else {
         Ref->Type|=GRID_REGULAR;
      }

      if (Ref->GRTYP[0]=='#') {
         Ref->Type|=GRID_TILE;
      }
    
      if (Ref->GRTYP[0]=='X' || Ref->GRTYP[0]=='O') {
         // If grid type is X (ORCA) and a pole in within the grid, mark as wrapping grid
         if (Ref->UnProject(Ref,&d[0],&d[1],89.0,0.0,0,1) || Ref->UnProject(Ref,&d[0],&d[1],-89.0,0.0,0,1)) {
            Ref->Type|=GRID_WRAP;
         }
      }  
     
      if (Ref->GRTYP[0]=='A' || Ref->GRTYP[0]=='B' || Ref->GRTYP[0]=='G') {
         Ref->Type|=GRID_WRAP;
      } else if (Ref->GRTYP[0]!='V' && Ref->X0!=Ref->X1 && Ref->Y0!=Ref->Y1) {
         // Check if north is up by looking at longitude variation on an Y increment at grid limits
         Ref->Project(Ref,Ref->X0,Ref->Y0,&co[0].Lat,&co[0].Lon,1,1);
         Ref->Project(Ref,Ref->X0,Ref->Y0+1,&co[1].Lat,&co[1].Lon,1,1);
         d[0]=co[0].Lon-co[1].Lon;
         Ref->Project(Ref,Ref->X1,Ref->Y1-1,&co[0].Lat,&co[0].Lon,1,1);
         Ref->Project(Ref,Ref->X1,Ref->Y1,&co[1].Lat,&co[1].Lon,1,1);
         d[1]=co[0].Lon-co[1].Lon;

         if (fabs(d[0])>0.0001 || fabs(d[1])>0.0001) {
            Ref->Type|=GRID_ROTATED;
         }
                  
         // Get size of a gridpoint
         Ref->Project(Ref,Ref->X0+(Ref->X1-Ref->X0)/2.0,Ref->Y0+(Ref->Y1-Ref->Y0)/2.0,&co[0].Lat,&co[0].Lon,1,1);
         Ref->Project(Ref,Ref->X0+(Ref->X1-Ref->X0)/2.0+1.0,Ref->Y0+(Ref->Y1-Ref->Y0)/2.0,&co[1].Lat,&co[1].Lon,1,1);
         d[0]=DIST(0.0,DEG2RAD(co[0].Lat),DEG2RAD(co[0].Lon),DEG2RAD(co[1].Lat),DEG2RAD(co[1].Lon));

         // Get distance between first and lat point
         Ref->Project(Ref,Ref->X0,Ref->Y0+(Ref->Y1-Ref->Y0)/2.0,&co[0].Lat,&co[0].Lon,1,1);
         Ref->Project(Ref,Ref->X1,Ref->Y0+(Ref->Y1-Ref->Y0)/2.0,&co[1].Lat,&co[1].Lon,1,1);
         d[1]=DIST(0.0,DEG2RAD(co[0].Lat),DEG2RAD(co[0].Lon),DEG2RAD(co[1].Lat),DEG2RAD(co[1].Lon));

         // If we're within 1.5 grid point, we wrap
         if (d[1]<=(d[0]*1.5)) {
            Ref->Type|=GRID_WRAP;
         }

         // If we're within 0.25 grid point, we repeat
         if (d[1]<=(d[0]*0.25)) {
            Ref->Type|=GRID_REPEAT;
         }
      }

      if (Ref->GRTYP[0]=='V') {
         Ref->Type|=GRID_VERTICAL;
      }

      if (Ref->GRTYP[0]=='R') {
         Ref->Type|=GRID_RADIAL;
      }     
   
      // Check for negative longitude (-180 <-> 180, or 0 <-> 360)
      if (Ref->AX) {
         for(x=0;x<Ref->NX;x++) { 
            if (Ref->AX[x]<0) { 
               Ref->Type|=GRID_NEGLON; 
               break; 
            } 
         } 
      }
   }
}

/**----------------------------------------------------------------------------
 * @brief  Test two GeoRef for equality
 * @author Jean-Philippe Gauthier
 * @date   July 2005
 *    @param[in]  Ref0    First geo reference
 *    @param[in]  Ref1    Second geo reference
 *
 *    @return             Equality (1=True 0=False)
*/
int GeoRef_Equal(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1) {

   if (!Ref0 || !Ref1) {
      return(0);
   }

   if (Ref0->GRTYP[0]!=Ref1->GRTYP[0] || Ref0->GRTYP[1]!=Ref1->GRTYP[1])
      return(0);
   
   if (Ref0->RPNHead.IG[X_IG1]!=Ref1->RPNHead.IG[X_IG1] || Ref0->RPNHead.IG[X_IG2]!=Ref1->RPNHead.IG[X_IG2] || Ref0->RPNHead.IG[X_IG3]!=Ref1->RPNHead.IG[X_IG3] || Ref0->RPNHead.IG[X_IG4]!=Ref1->RPNHead.IG[X_IG4])
     return(0);
   //TOTO: Check AX,AY

   // Cloud point should never be tested as equal
   if (Ref0->GRTYP[0]=='Y' || Ref0->GRTYP[1]=='Y')
      return(0);

   if (Ref1->GRTYP[0]=='Y' || Ref1->GRTYP[1]=='Y')
      return(0);

   // Check on grid limits (exclude U grids since they can be switched internally and they'll be tested earlier anyway)
   if (Ref0->BD!=Ref1->BD || (Ref0->GRTYP[0]!='U' && (Ref0->X0!=Ref1->X0 || Ref0->X1!=Ref1->X1 || Ref0->Y0!=Ref1->Y0 || Ref0->Y1!=Ref1->Y1)))
      return(0);
   
   if (Ref0->Subs && Ref1->Subs && Ref0->Subs[0]!=Ref1->Subs[0])
       return(0);

   if (Ref0->R!=Ref1->R || Ref0->ResR!=Ref1->ResR || Ref0->ResA!=Ref1->ResA || Ref0->Loc.Lat!=Ref1->Loc.Lat || Ref0->Loc.Lon!=Ref1->Loc.Lon || Ref0->Loc.Elev!=Ref1->Loc.Elev)
      return(0);

   if ((Ref0->Spatial && !Ref1->Spatial) || (!Ref0->Spatial && Ref1->Spatial))
      return(0);

#ifdef HAVE_GDAL 
   if (Ref0->Spatial && Ref1->Spatial && !OSRIsSame(Ref0->Spatial,Ref1->Spatial))
      return(0);
#endif

   if ((Ref0->RotTransform && !Ref1->RotTransform) || (!Ref0->RotTransform && Ref1->RotTransform))
      return(0);
   if (Ref0->RotTransform && Ref1->RotTransform)
      if (memcmp(Ref0->RotTransform,Ref1->RotTransform,sizeof(TRotationTransform))!=0)
         return(0);

   if ((Ref0->Transform && !Ref1->Transform) || (!Ref0->Transform && Ref1->Transform))
      return(0);
   if (Ref0->Transform && Ref1->Transform)
      if (memcmp(Ref0->Transform,Ref1->Transform,6*sizeof(double))!=0)
         return(0);

   return(1);
}

/**----------------------------------------------------------------------------
 * @brief  Copy a GeoRef only by incrementing the reference count
 * @author Jean-Philippe Gauthier
 * @date   Janvier 2015
 *    @param[in]  Ref     Pointeur sur la reference
 *
 *    @return             Pointeur sur la copie de la reference
*/
TGeoRef *GeoRef_Copy(TGeoRef* __restrict const Ref) {

   GeoRef_Incr(Ref);
   return(Ref);
}

/**----------------------------------------------------------------------------
 * @brief  Create e new GeoRef but link it to an already existing GeoRef
 * @author Jean-Philippe Gauthier
 * @date   Janvier 2015
 *    @param[in]  Ref     Pointeur sur la reference
 *
 *    @return             Pointeur sur la copie de la reference
*/
TGeoRef *GeoRef_Reference(TGeoRef* __restrict const Ref) {

   TGeoRef *ref;

   ref=GeoRef_New();

   if (Ref) {
      GeoRef_Incr(Ref);
      ref->RefFrom=Ref;
      ref->Project=Ref->Project;
      ref->UnProject=Ref->UnProject;
      ref->Value=Ref->Value;
      ref->Distance=Ref->Distance;

      ref->Loc.Lat=Ref->Loc.Lat;
      ref->Loc.Lon=Ref->Loc.Lon;
      ref->Loc.Elev=Ref->Loc.Elev;
      ref->R=Ref->R;
      ref->ResR=Ref->ResR;
      ref->ResA=Ref->ResA;
   }

   return(ref);
}

/**----------------------------------------------------------------------------
 * @brief  Make a hard(real) copy of the GeoRef structure, not just a reference count increment
 * @author Jean-Philippe Gauthier
 * @date   Janvier 2015
 *    @param[in]  Ref     Pointeur sur la reference
 *
 *    @return             Pointeur sur la copie de la reference
*/
TGeoRef *GeoRef_HardCopy(TGeoRef* __restrict const Ref) {

   TGeoRef *ref;
   int      i;

   ref=GeoRef_New();
   GeoRef_Size(ref,Ref->X0,Ref->Y0,Ref->X1,Ref->Y1,Ref->BD);

   if (Ref) {
      ref->GRTYP[0]=Ref->GRTYP[0];
      ref->GRTYP[1]=Ref->GRTYP[1];
      ref->Project=Ref->Project;
      ref->UnProject=Ref->UnProject;
      ref->Value=Ref->Value;
      ref->Distance=Ref->Distance;
      ref->Type=Ref->Type;
      ref->NbSub=Ref->NbSub;
      ref->QTree=NULL;
      
#ifdef HAVE_RMN
      if (Ref->Subs) {
         ref->Subs=(TGeoRef**)malloc(Ref->NbSub*sizeof(TGeoRef*));
         memcpy(ref->Subs,Ref->Subs,Ref->NbSub*sizeof(TGeoRef*));
         for(i=0;i<ref->NbSub;i++)
            GeoRef_Incr(ref->Subs[i]);
      }
#endif

      memcpy(&ref->RPNHead,&Ref->RPNHead,sizeof(TRPNHeader));
      memcpy(&ref->Options,&Ref->Options,sizeof(TGeoOptions));

      switch(ref->GRTYP[0]) {
         case 'R' :
            ref->Loc.Lat=Ref->Loc.Lat;
            ref->Loc.Lon=Ref->Loc.Lon;
            ref->Loc.Elev=Ref->Loc.Elev;
            ref->R=Ref->R;
            ref->ResR=Ref->ResR;
            ref->ResA=Ref->ResA;
         case 'W' :
            GeoRef_WKTSet(ref,Ref->String,Ref->Transform,Ref->InvTransform,Ref->Spatial);
            if (Ref->RotTransform) {
               ref->RotTransform=(TRotationTransform*)malloc(sizeof(TRotationTransform));
               memcpy(ref->RotTransform,Ref->RotTransform,sizeof(TRotationTransform));
            }
      }
   }
   return(ref);
}

/**----------------------------------------------------------------------------
 * @brief  Resize a geo reference
 * @author Jean-Philippe Gauthier
 * @date   Janvier 2015
 *    @param[in]  Ref     Pointeur sur la reference
 *    @param[in]  NI      New x size
 *    @param[in]  NJ      New y size
 *
 *    @return             Pointer to the new geo reference
*/
TGeoRef *GeoRef_Resize(TGeoRef* __restrict const Ref,int NI,int NJ) {

   TGeoRef *ref;

   if (!Ref) {
      ref=GeoRef_New();
   } else {
      ref=GeoRef_HardCopy(Ref);
   }
   GeoRef_Size(ref,0,0,NI-1,NJ-1,0);

   return(ref);
}

int GeoRef_Project(TGeoRef* __restrict const Ref,double X,double Y,double *Lat,double *Lon,int Extrap,int Transform) {

   if (!Ref) return(0);
   
   if (X>Ref->X1 || Y>Ref->Y1 || X<Ref->X0 || Y<Ref->Y0) {
      if (!Extrap) {
         *Lon=-999.0;
         *Lat=-999.0;
         return(0);
      }
   }
   *Lon=X/(Ref->X1-Ref->X0);
   *Lat=Y/(Ref->Y1-Ref->Y0);

   return(1);
}

int GeoRef_UnProject(TGeoRef* __restrict const Ref,double *X,double *Y,double Lat,double Lon,int Extrap,int Transform) {

   if (!Ref) return(0);
   
   *X=Lon*(Ref->X1-Ref->X0);
   *Y=Lat*(Ref->Y1-Ref->Y0);

   /*Check the grid limits*/
   if (*X>Ref->X1 || *Y>Ref->Y1 || *X<Ref->X0 || *Y<Ref->Y0) {
      if (!Extrap) {
         *X=-1.0;
         *Y=-1.0;
      }
      return(0);
   }
   return(1);
}

/**----------------------------------------------------------------------------
 * @brief  Add a geo reference to the list of known geo reference
 * @author Jean-Philippe Gauthier
 * @date   Janvier 2019
 *    @param[in]  Ref     Pointeur sur la reference
 *
 *    @return             Pointer to the new head of the list
*/
TGeoRef* GeoRef_Add(TGeoRef *Ref) {

   TList *head;

   GeoRef_Lock();
   if (head=TList_Add(GeoRef_List,(void*)Ref)) {
      GeoRef_List=head;
   }
   GeoRef_Unlock();

   return((TGeoRef*)(head?head->Data:NULL));
}

/**----------------------------------------------------------------------------
 * @brief  Find a geo reference
 * @author Jean-Philippe Gauthier
 * @date   Janvier 2015
 *    @param[in]  Ref     Pointeur sur la reference
 *
 *    @return             Pointer to found geo reference (NULL if not found)
*/
TGeoRef* GeoRef_Find(TGeoRef *Ref) {

   TList *item;

   GeoRef_Lock();
   item=TList_Find(GeoRef_List,(TList_CompareProc*)GeoRef_Equal,(void*)Ref);
   GeoRef_Unlock();

   if (item) {
      return((TGeoRef*)item->Data);
   }
   return(NULL);
}

/**----------------------------------------------------------------------------
 * @brief  Create new geo reference
 * @author Jean-Philippe Gauthier
 * @date   Janvier 2015
 *
 *    @return              New geo reference pointer
*/
TGeoRef* GeoRef_New() {

   TGeoRef *ref=malloc(sizeof(TGeoRef));

   GeoRef_Size(ref,0,0,0,0,0);

   // General
   ref->Name=NULL;
   ref->Subs=NULL;
   ref->NbSub=0;
   ref->Type=GRID_NONE;
   ref->NRef=1;
   ref->NIdx=0;
   ref->Lat=NULL;
   ref->Lon=NULL;
   ref->Hgt=NULL;
   ref->Wght=NULL;
   ref->Idx=NULL;
   ref->AX=NULL;
   ref->AY=NULL;
   ref->NCX=NULL;
   ref->NCY=NULL;
   ref->RefFrom=NULL;
   ref->QTree=NULL;
   ref->GRTYP[0]='X';
   ref->GRTYP[1]='\0';
   ref->GRTYP[2]='\0';
   ref->Sets=NULL;
   ref->LastSet=NULL;

   // Assign defualt options
   memcpy(&ref->Options,&GeoRef_Options,sizeof(TGeoOptions));

   // RPN Specific
   memset(&ref->RPNHead,0x0,sizeof(TRPNHeader));

   // WKT Specific
   ref->String=NULL;
   ref->Spatial=NULL;
   ref->Function=NULL;
   ref->InvFunction=NULL;
   ref->Transform=NULL;
   ref->RotTransform=NULL;
   ref->InvTransform=NULL;
   ref->GCPTransform=NULL;
   ref->TPSTransform=NULL;
   ref->RPCTransform=NULL;
   ref->LLExtent.MinX=1e32;
   ref->LLExtent.MinY=1e32;
   ref->LLExtent.MaxX=-1e32;
   ref->LLExtent.MaxY=-1e32;

   // RDR Specific
   ref->Loc.Lat=-999;
   ref->Loc.Lon=-999;
   ref->Loc.Elev=-999;
   ref->R=0;
   ref->CTH=0;
   ref->STH=0;
   ref->ResR=0;
   ref->ResA=0;

   // General functions
   ref->Project=GeoRef_Project;
   ref->UnProject=GeoRef_UnProject;
   ref->Value=NULL;
   ref->Distance=NULL;
   ref->Height=NULL;

   return(ref);
}

/**----------------------------------------------------------------------------
 * @brief  Create spatial index
 * @author Jean-Philippe Gauthier
 * @date   Janvier 2016
 *    @param[in]  Ref     Pointeur sur la reference
 *
 *    @return             Quad tree spatial index
*/
TQTree* GeoRef_BuildIndex(TGeoRef* __restrict const Ref) {

   unsigned int  n,x,y,t;
   double        dx,dy,lat0,lon0,lat1,lon1;
   Vect2d        tr[3],pt;
   
   if (!Ref->AX || !Ref->AY) {
      return(NULL);
   }
            
   // Check data limits   
   lat0=lon0=1e10;
   lat1=lon1=-1e10;
   
   for(n=0;n<Ref->NX*Ref->NY;n++) {
      dy=Ref->AY[n];
      dx=CLAMPLON(Ref->AX[n]);
      lat0=fmin(lat0,dy);
      lon0=fmin(lon0,dx);
      lat1=fmax(lat1,dy);
      lon1=fmax(lon1,dx);
   }
      
   if (Ref->GRTYP[0]=='M') {
      
      // Allocate barycentric weight array if needed
      if (!Ref->Wght) {
         if (!(Ref->Wght=(double*)calloc(Ref->NIdx/3,sizeof(double)))) {
            App_Log(WARNING,"%s: Failed to allocate baricentric weight array\n",__func__);
         }
      }
      
      // Create the tree on the data limits
      if (!(Ref->QTree=QTree_New(lon0,lat0,lon1,lat1,NULL))) {
         App_Log(ERROR,"%s: Failed to create QTree index\n",__func__);
         return(NULL);
      }

      // Loop on triangles
      for(n=0,t=0;n<Ref->NIdx-3;n+=3,t++) {          
         tr[0][0]=Ref->AX[Ref->Idx[n]];     tr[0][1]=Ref->AY[Ref->Idx[n]];
         tr[1][0]=Ref->AX[Ref->Idx[n+1]];   tr[1][1]=Ref->AY[Ref->Idx[n+1]];
         tr[2][0]=Ref->AX[Ref->Idx[n+2]];   tr[2][1]=Ref->AY[Ref->Idx[n+2]];
         
         // Calculate barycentric weight
          if (Ref->Wght)
             Ref->Wght[t]=1.0/((tr[1][0]-tr[0][0])*(tr[2][1]-tr[0][1])-(tr[2][0]-tr[0][0])*(tr[1][1]-tr[0][1]));
 
         // Put it in the quadtree, in any child nodes intersected and set false pointer increment (+1)
         if (!QTree_AddTriangle(Ref->QTree,tr,GRID_MQTREEDEPTH,(void*)(n+1))) {
            App_Log(ERROR,"%s: Failed to add node\n",__func__);
            return(NULL);
         }      
      }
   } else  if (Ref->GRTYP[0]=='Y' || Ref->GRTYP[0]=='X' || Ref->GRTYP[0]=='O' || Ref->GRTYP[1]=='Y' || Ref->GRTYP[1]=='X') {

      // Useless for less than a few thousand points
      if ((Ref->NX*Ref->NY)<500) {
         return(NULL);
      }

      // Create the array on the data limits
      dy=(lat1-lat0)/GRID_YQTREESIZE;
      dx=(lon1-lon0)/GRID_YQTREESIZE;

      if (!(Ref->QTree=(TQTree*)calloc((GRID_YQTREESIZE+1)*(GRID_YQTREESIZE+1),sizeof(TQTree)))) {
         App_Log(ERROR,"%s: Failed to create QTree index\n",__func__);
         return(NULL);
      }

      // Store tree limit on first node
      Ref->QTree[0].BBox[0].X=lon0;
      Ref->QTree[0].BBox[0].Y=lat0;
      Ref->QTree[0].BBox[1].X=lon1;
      Ref->QTree[0].BBox[1].Y=lat1;
   
      // Loop on points
      for(n=0;n<Ref->NX*Ref->NY;n++) {     
         
         pt[0]=Ref->AX[n];
         pt[1]=Ref->AY[n];  
         
         x=(pt[0]-lon0)/dx;
         y=(pt[1]-lat0)/dy;
               
         // Add location and set false pointer increment (+1);
         QTree_AddData(&Ref->QTree[y*GRID_YQTREESIZE+x],pt[0],pt[1],(void*)(n+1));
      }
   }
   
   return(Ref->QTree);
}

/**----------------------------------------------------------------------------
 * @brief  Find closest point(s) to a grid position
 * @author Jean-Philippe Gauthier
 * @date   Janvier 2015
 *    @param[in]  Ref     Pointeur sur la reference
 *    @param[in]  X       X coordinate
 *    @param[in]  Y       Y coordinate
 *    @param[out] Idxs    Pointer to neighbors index found
 *    @param[out] Dists   Squared distances from the neighbors found
 *    @param[in]  NbNear  Number of nearest neighbors to find
 *    @param[in]  MaxDist Maximum distance to find (0.0 = don't care)
 *
 *    @return             Nombre de points trouvé trié du plus près vers le plus loin
*/
int GeoRef_Nearest(TGeoRef* __restrict const Ref,double X,double Y,int *Idxs,double *Dists,int NbNear,double MaxDist) {

   double       dx,dy,l;
   unsigned int n,nn,nr,nnear;
   TQTree      *node;
   int          dxy,x,y,xd,yd,rx;
  
   if (!NbNear || !Idxs || !Dists) return(0);

   Dists[0]=1e32;
   dxy=nnear=0;

   if (Ref->QTree) {     
      
      // Find the closest point(s) by circling larger around cell      
      node=&Ref->QTree[0];

      if (!(Ref->Type&GRID_WRAP) && !FWITHIN(0,node->BBox[0].Y,node->BBox[0].X,node->BBox[1].Y,node->BBox[1].X,Y,X)) {
         return(0);
      }

      dx=(node->BBox[1].X-node->BBox[0].X)/GRID_YQTREESIZE;
      dy=(node->BBox[1].Y-node->BBox[0].Y)/GRID_YQTREESIZE;
      xd=(X-node->BBox[0].X)/dx;
      yd=(Y-node->BBox[0].Y)/dy;
                    
      while(dxy<(GRID_YQTREESIZE>>1)) {

         // Y circling increment
         for(y=yd-dxy;y<=yd+dxy;y++) {
            if (y<0) continue;
            if (y>=GRID_YQTREESIZE) break;
            
            // X Circling increment (avoid revisiting previous cells)
            for(x=xd-dxy;x<=xd+dxy;x+=((!dxy || y==yd-dxy || y==yd+dxy)?1:(dxy+dxy))) {
               rx=x;
               if (x<0) {
                  if (Ref->Type&GRID_WRAP) {
                     rx=x+GRID_YQTREESIZE;
                  } else {
                     continue;                  
                  }
               }
               
               if (x>=GRID_YQTREESIZE) {
                  if (Ref->Type&GRID_WRAP) {
                     rx=x-GRID_YQTREESIZE;
                  } else {
                     break;
                  }
               }
               
               node=&Ref->QTree[y*GRID_YQTREESIZE+rx];

               // Loop on points in this cell and get closest point
               for(n=0;n<node->NbData;n++) {
                      
                  dx=X-node->Data[n].Pos.X;
                  // flip distances for wrap-around
                  if (dx>=180)  dx-=360;
                  if (dx<=-180) dx+=360;
                  
                  dy=Y-node->Data[n].Pos.Y;
                  l=dx*dx+dy*dy;

                  // Loop on number of nearest to find
                  for(nn=0;nn<NbNear;nn++) {
                     // If this is closer
                     if ((MaxDist==0.0 || l<=MaxDist) && l<Dists[nn]) {
                           
                        // Move farther nearest in order
                        for(nr=NbNear-1;nr>nn;nr--) {
                           Dists[nr]=Dists[nr-1];
                           Idxs[nr]=Idxs[nr-1];                       
                        }
                           
                        // Assign found nearest
                        Dists[nn]=l;
                        Idxs[nn]=(intptr_t)node->Data[n].Ptr-1; // Remove false pointer increment
                        nnear++;
                        break;
                     }
                  }
               }
            }            
         }
         
         // If we made at least one cycle and found enough nearest
         if (nnear>=NbNear && dxy) 
            break;
         dxy++;
      }
   } else {
      // Find closest by looping in all points   
      for(n=0;n<Ref->NX*Ref->NY;n++) {
         dx=X-Ref->AX[n];
         dy=Y-Ref->AY[n];
         
         if (Ref->Type&GRID_WRAP) {
            if (dx>180)  dx-=360;
            if (dx<-180) dx+=360;
         }
         
         l=dx*dx+dy*dy;
         
         for(nn=0;nn<NbNear;nn++) {
            if (l<Dists[nn]) {
                  
               // Move farther nearest in order
               for(nr=NbNear-1;nr>nn;nr--) {
                  Dists[nr]=Dists[nr-1];
                  Idxs[nr]=Idxs[nr-1];
               }
               
               // Assign found nearest
               Dists[nn]=l;
               Idxs[nn]=n;
               
               nnear++;
               break;
            }
         }
      }
   }

   // Return found index
   return(nnear>NbNear?NbNear:nnear);
}

/**----------------------------------------------------------------------------
 * @brief  Check for intersection of 2 gep reference
 * @author Jean-Philippe Gauthier
 * @date   Aout 2006 
 *    @param[in]  Ref0     Pointeur sur la reference 0
 *    @param[in]  Ref1     Pointeur sur la reference 1
 *    @param[out] X0       First x corner limit of Ref1 in Ref0
 *    @param[out] Y0       First y corner limit of Ref1 in Ref0
 *    @param[out] X1       Second x corner limit of Ref1 in Ref0
 *    @param[out] Y1       Second y corner limit of Ref1 in Ref0
 *    @param[in]  BD       Include border flag
 *
 *    @return              Intersection (True or False)
*/
int GeoRef_Intersect(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1,int *X0,int *Y0,int *X1,int *Y1,int BD) {

   double lat,lon,di,dj,in=0;
   double x0,y0,x1,y1;
   int    x,y;

   if (!Ref0 || !Ref1) return(0);

   /*Source grid Y*/
   if (Ref1->GRTYP[0]=='Y') {
      *X0=Ref1->X0; *Y0=Ref1->Y0;
      *X1=Ref1->X1; *Y1=Ref1->Y1;
      return(1);
   }

   /*If destination is global*/
   if (Ref0->Type&GRID_WRAP) {
      *X0=Ref1->X0; *Y0=Ref1->Y0;
      *X1=Ref1->X1; *Y1=Ref1->Y1;
      in=1;
   }

   /*Test for limit source inclusion into destination*/
   x0=y0=1e32;
   x1=y1=-1e32;
   x=0;

   if (!in) {
      Ref1->Project(Ref1,Ref1->X0,Ref1->Y0,&lat,&lon,0,1);
      if (Ref0->UnProject(Ref0,&di,&dj,lat,lon,0,1)) {
         x0=Ref1->X0;y0=Ref1->Y0;
         x++;
      }
      Ref1->Project(Ref1,Ref1->X0,Ref1->Y1,&lat,&lon,0,1);
      if (Ref0->UnProject(Ref0,&di,&dj,lat,lon,0,1)) {
         x0=Ref1->X0;y1=Ref1->Y1;
         x++;
      }
      Ref1->Project(Ref1,Ref1->X1,Ref1->Y0,&lat,&lon,0,1);
      if (Ref0->UnProject(Ref0,&di,&dj,lat,lon,0,1)) {
         x1=Ref1->X1;y0=Ref1->Y0;
         x++;
      }
      Ref1->Project(Ref1,Ref1->X1,Ref1->Y1,&lat,&lon,0,1);
      if (Ref0->UnProject(Ref0,&di,&dj,lat,lon,0,1)) {
         x1=Ref1->X1;y1=Ref1->Y1;
      }
      *X0=x0; *Y0=y0;
      *X1=x1; *Y1=y1;

      if (x>=3) {
         in=1;
      }
   }

   if (!in) {

      /*Project Ref0 within Ref1 and get limits*/
      for(x=Ref0->X0;x<=Ref0->X1;x++) {
         Ref0->Project(Ref0,x,Ref0->Y0,&lat,&lon,0,1);
         Ref1->UnProject(Ref1,&di,&dj,lat,lon,1,1);
         x0=fmin(x0,di); y0=fmin(y0,dj);
         x1=fmax(x1,di); y1=fmax(y1,dj);

         Ref0->Project(Ref0,x,Ref0->Y1,&lat,&lon,0,1);
         Ref1->UnProject(Ref1,&di,&dj,lat,lon,1,1);
         x0=fmin(x0,di); y0=fmin(y0,dj);
         x1=fmax(x1,di); y1=fmax(y1,dj);
      }

      for(y=Ref0->Y0;y<=Ref0->Y1;y++) {
         Ref0->Project(Ref0,Ref0->X0,y,&lat,&lon,0,1);
         Ref1->UnProject(Ref1,&di,&dj,lat,lon,1,1);
         x0=fmin(x0,di); y0=fmin(y0,dj);
         x1=fmax(x1,di); y1=fmax(y1,dj);

         Ref0->Project(Ref0,Ref0->X1,y,&lat,&lon,0,1);
         Ref1->UnProject(Ref1,&di,&dj,lat,lon,1,1);
         x0=fmin(x0,di); y0=fmin(y0,dj);
         x1=fmax(x1,di); y1=fmax(y1,dj);
      }

      /*Test for north and south pole including grid*/
      if (Ref0->UnProject(Ref0,&di,&dj,89.9,0.0,0,1) && dj>Ref0->Y0+2 && dj<Ref0->Y1-2 && di>Ref0->X0+2 && di<Ref0->X1-2) {
         Ref1->UnProject(Ref1,&di,&dj,89.9,0.0,1,1);
         x0=fmin(x0,di); y0=fmin(y0,dj);
         x1=fmax(x1,di); y1=fmax(y1,dj);
      }
      if (Ref0->UnProject(Ref0,&di,&dj,-89.9,0.0,0,1) && dj>Ref0->Y0+2 && dj<Ref0->Y1-2 && di>Ref0->X0+2 && di<Ref0->X1-2) {
         Ref1->UnProject(Ref1,&di,&dj,-89.9,0.0,1,1);
         x0=fmin(x0,di); y0=fmin(y0,dj);
         x1=fmax(x1,di); y1=fmax(y1,dj);
      }

      *X0=floor(x0); *Y0=floor(y0);
      *X1=ceil(x1);  *Y1=ceil(y1);

      if (!VOUT(*X0,*X1,Ref1->X0,Ref1->X1) && !VOUT(*Y0,*Y1,Ref1->Y0,Ref1->Y1)) {
         in=1;
      }
   }

   /*Clamp the coordinates*/
   if (BD) {
      REFCLAMP(Ref1,*X0,*Y0,*X1,*Y1);
   } else {
      REFCLAMPBD(Ref1,*X0,*Y0,*X1,*Y1);
   }

   return(in);
}

/**----------------------------------------------------------------------------
 * @brief  Calculates geographical limits of a geo reference
 * @author Jean-Philippe Gauthier
 * @date   Aout 2006
 *    @param[in]  Ref      Pointeur sur la reference
 *    @param[out] Lat0     Latitude of first corner
 *    @param[out] Lon0     Longitude of first corner
 *    @param[out] Lat1     Latitude of second corner
 *    @param[out] Lon1     Longitude of second corner
 *
 *    @return              Error code
*/
int GeoRef_Limits(TGeoRef* __restrict const Ref,double *Lat0,double *Lon0,double *Lat1,double *Lon1) {

   int x,y;
   double di,dj,lat,lon;

   *Lat0=*Lon0=1e32;
   *Lat1=*Lon1=-1e32;

   if (!Ref) return(0);
   
   // Source grid Y
   if (Ref->GRTYP[0]=='Y' || Ref->GRTYP[1]=='Y' || Ref->GRTYP[0]=='X' || Ref->GRTYP[1]=='X' || Ref->GRTYP[0]=='O') {
      for(x=0;x<((Ref->X1-Ref->X0)+1)*((Ref->Y1-Ref->Y0)+1);x++) {
         *Lat0=fmin(*Lat0,Ref->AY[x]); *Lon0=fmin(*Lon0,Ref->AX[x]);
         *Lat1=fmax(*Lat1,Ref->AY[x]); *Lon1=fmax(*Lon1,Ref->AX[x]);
      }
      return(1);
   }

   // If destination is global
   if (Ref->Type&GRID_WRAP) {
      *Lat0=-90.0;  *Lat1=90.0;
      *Lon0=-180.0; *Lon1=180.0;
      return(1);
   }

   // Project Ref0 Border within Ref1 and get limits
   for(x=Ref->X0,y=Ref->Y0;x<=Ref->X1;x++) {
      Ref->Project(Ref,x,y,&lat,&lon,0,1);
      *Lat0=fmin(*Lat0,lat); *Lon0=fmin(*Lon0,lon);
      *Lat1=fmax(*Lat1,lat); *Lon1=fmax(*Lon1,lon);
   }

   for(x=Ref->X0,y=Ref->Y1;x<=Ref->X1;x++) {
      Ref->Project(Ref,x,y,&lat,&lon,0,1);
      *Lat0=fmin(*Lat0,lat); *Lon0=fmin(*Lon0,lon);
      *Lat1=fmax(*Lat1,lat); *Lon1=fmax(*Lon1,lon);
   }

   for(y=Ref->Y0,x=Ref->X0;y<=Ref->Y1;y++) {
      Ref->Project(Ref,x,y,&lat,&lon,0,1);
      *Lat0=fmin(*Lat0,lat); *Lon0=fmin(*Lon0,lon);
      *Lat1=fmax(*Lat1,lat); *Lon1=fmax(*Lon1,lon);
   }

   for(y=Ref->Y0,x=Ref->X1;y<=Ref->Y1;y++) {
      Ref->Project(Ref,x,y,&lat,&lon,0,1);
      *Lat0=fmin(*Lat0,lat); *Lon0=fmin(*Lon0,lon);
      *Lat1=fmax(*Lat1,lat); *Lon1=fmax(*Lon1,lon);
   }

   // Test for north and south pole including grid
   if (Ref->UnProject(Ref,&di,&dj,90.0,0.0,0,1) && dj>Ref->Y0+2 && dj<Ref->Y1-2 && di>Ref->X0+2 && di<Ref->X1-2) {
      *Lat1=90.0;
   }
   if (Ref->UnProject(Ref,&di,&dj,-90.0,0.0,0,1) && dj>Ref->Y0+2 && dj<Ref->Y1-2 && di>Ref->X0+2 && di<Ref->X1-2) {
      *Lat0=-90.0;
   }
   return(1);
}

/**----------------------------------------------------------------------------
 * @brief  Check if a georef is within another
 * @author Jean-Philippe Gauthier
 * @date   February 2009
 *    @param[in]  Ref0     Pointeur sur la reference inclusive
 *    @param[in]  Ref1     Pointeur sur la reference a tester
  *
 *    @return              True or False
*/
int GeoRef_Within(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1) {

   double lat,lon,di,dj;
   int    x,y;

   if (!Ref0 || !Ref1) return(0);
   
   // Project Ref0 Border within Ref1 and get limits
   for(x=Ref0->X0,y=Ref0->Y0;x<=Ref0->X1;x++) {
      Ref0->Project(Ref0,x,y,&lat,&lon,0,1);
      if (!Ref1->UnProject(Ref1,&di,&dj,lat,lon,0,1)) {
         return(0);
      }
   }

   for(x=Ref0->X0,y=Ref0->Y1;x<=Ref0->X1;x++) {
      Ref0->Project(Ref0,x,y,&lat,&lon,0,1);
      if (!Ref1->UnProject(Ref1,&di,&dj,lat,lon,0,1)) {
         return(0);
      }
   }

   for(y=Ref0->Y0,x=Ref0->X0;y<=Ref0->Y1;y++) {
      Ref0->Project(Ref0,x,y,&lat,&lon,0,1);
      if (!Ref1->UnProject(Ref1,&di,&dj,lat,lon,0,1)) {
         return(0);
      };
   }

   for(y=Ref0->Y0,x=Ref0->X1;y<=Ref0->Y1;y++) {
      Ref0->Project(Ref0,x,y,&lat,&lon,0,1);
      if (!Ref1->UnProject(Ref1,&di,&dj,lat,lon,0,1)) {
         return(0);
      }
   }
   return(1);
}
                    
/**----------------------------------------------------------------------------
 * @brief  Check if a geographic range is within a georef
 * @author Jean-Philippe Gauthier
 * @date   February 2009
 *    @param[in]  Ref      Pointeur sur la reference
 *    @param[in]  Lat0     Latitude of first corner
 *    @param[in]  Lon0     Longitude of first corner
 *    @param[in]  Lat1     Latitude of second corner
 *    @param[in]  Lon1     Longitude of second corner
 *    @param[in]  In       Check ofr inside or outside
 * 
 *    @return              True or False
*/
int GeoRef_WithinRange(TGeoRef* __restrict const Ref,double Lat0,double Lon0,double Lat1,double Lon1,int In) {

   double lat[4],lon[4],dl;
   int    d0,d1,d2,d3;

   if (!Ref) return(0);
   
   if (Lat0>Lat1) {
      dl=Lat0;
      Lat0=Lat1;
      Lat1=dl;
   }

   if (Lon0>Lon1) {
      dl=Lon0;
      Lon0=Lon1;
      Lon1=dl;
   }

   if (Lon0*Lon1<0) {
      dl=Lon1-Lon0;
   } else {
      dl=0;
   }
   
   /* Check image within range */
   Ref->Project(Ref,Ref->X0,Ref->Y0,&lat[0],&lon[0],0,1);
   d0=FWITHIN(dl,Lat0,Lon0,Lat1,Lon1,lat[0],lon[0]);
   if (!In && d0) return(1);

   Ref->Project(Ref,Ref->X1,Ref->Y0,&lat[1],&lon[1],0,1);
   d1=FWITHIN(dl,Lat0,Lon0,Lat1,Lon1,lat[1],lon[1]);
   if (!In && d1) return(1);

   Ref->Project(Ref,Ref->X1,Ref->Y1,&lat[2],&lon[2],0,1);
   d2=FWITHIN(dl,Lat0,Lon0,Lat1,Lon1,lat[2],lon[2]);
   if (!In && d2) return(1);

   Ref->Project(Ref,Ref->X0,Ref->Y1,&lat[3],&lon[3],0,1);
   d3=FWITHIN(dl,Lat0,Lon0,Lat1,Lon1,lat[3],lon[3]);
   if (!In && d3) return(1);

   /* Check for all contained */
   if (In) {
      if (d0 && d1 && d2 && d3) {
         return(1);
      } else {
         return(0);
      }
   }

   /* Check range within image */
   lat[0]=fmin(fmin(fmin(lat[0],lat[1]),lat[2]),lat[3]);
   lat[1]=fmax(fmax(fmax(lat[0],lat[1]),lat[2]),lat[3]);
   lon[0]=fmin(fmin(fmin(lon[0],lon[1]),lon[2]),lon[3]);
   lon[1]=fmax(fmax(fmax(lon[0],lon[1]),lon[2]),lon[3]);

   if (FWITHIN(dl,lat[0],lon[0],lat[1],lon[1],Lat0,Lon0)) return(1);
   if (FWITHIN(dl,lat[0],lon[0],lat[1],lon[1],Lat0,Lon1)) return(1);
   if (FWITHIN(dl,lat[0],lon[0],lat[1],lon[1],Lat1,Lon1)) return(1);
   if (FWITHIN(dl,lat[0],lon[0],lat[1],lon[1],Lat1,Lon0)) return(1);

   return(0);
}

int GeoRef_WithinCell(TGeoRef *Ref,Vect2d Pos,Vect2d Pt[4],int Idx0,int Idx1,int Idx2,int Idx3) {
 
   Vect3d b;
   int    t0,t1,sz=Ref->NX*Ref->NY;
   
   if (Idx0<sz && Idx1<sz && Idx2<sz && Idx3<sz && Idx0>=0 && Idx1>=0 && Idx2>=0 && Idx3>=0) {
      
      Pt[0][0]=Ref->AX[Idx0]; Pt[0][1]=Ref->AY[Idx0];
      Pt[1][0]=Ref->AX[Idx1]; Pt[1][1]=Ref->AY[Idx1];
      Pt[2][0]=Ref->AX[Idx2]; Pt[2][1]=Ref->AY[Idx2];
      Pt[3][0]=Ref->AX[Idx3]; Pt[3][1]=Ref->AY[Idx3];
      
      // Make sure all coordinates are on same side of -180/180
      if (Pos[0]>90) {
         if (Pt[0][0]<0) Pt[0][0]+=360;   
         if (Pt[1][0]<0) Pt[1][0]+=360;   
         if (Pt[2][0]<0) Pt[2][0]+=360;   
         if (Pt[3][0]<0) Pt[3][0]+=360;         
      } else if (Pos[0]<-90) {
         if (Pt[0][0]>0) Pt[0][0]-=360;   
         if (Pt[1][0]>0) Pt[1][0]-=360;   
         if (Pt[2][0]>0) Pt[2][0]-=360;   
         if (Pt[3][0]>0) Pt[3][0]-=360;        
      }      
      t0=Bary_Get(b,0.0,Pos[0],Pos[1],Pt[0][0],Pt[0][1],Pt[1][0],Pt[1][1],Pt[2][0],Pt[2][1]);
      t1=Bary_Get(b,0.0,Pos[0],Pos[1],Pt[0][0],Pt[0][1],Pt[2][0],Pt[2][1],Pt[3][0],Pt[3][1]);

      return(t0 || t1);
   }
   return(0);
}
                     
/**----------------------------------------------------------------------------
 * @brief  Project a latlon geographic bounding box into grid coordinates
 * @author Jean-Philippe Gauthier
 * @date   February 2009
 *    @param[in]  Ref      Pointeur sur la reference
 *    @param[in]  Lat0     Latitude of first corner
 *    @param[in]  Lon0     Longitude of first corner
 *    @param[in]  Lat1     Latitude of second corner
 *    @param[in]  Lon1     Longitude of second corner
 *    @param[out] I0       X coordinate of first corner
 *    @param[out] J0       Y coordinate of first corner
 *    @param[out] I1       X coordinate of seconf corner
 *    @param[out] J1       Y coordinate of second corner
 *
 *    @return              Error code
*/
int GeoRef_BoundingBox(TGeoRef* __restrict const Ref,double Lat0,double Lon0,double Lat1,double Lon1,double *I0,double *J0,double *I1,double *J1) {

   double di,dj;

   if (!Ref) return(0);
   
   *I0=*J0=1000000;
   *I1=*J1=-1000000;

   Ref->UnProject(Ref,&di,&dj,Lat0,Lon0,1,1);
   *I0=*I0<di?*I0:di;
   *J0=*J0<dj?*J0:dj;
   *I1=*I1>di?*I1:di;
   *J1=*J1>dj?*J1:dj;

   Ref->UnProject(Ref,&di,&dj,Lat0,Lon1,1,1);
   *I0=*I0<di?*I0:di;
   *J0=*J0<dj?*J0:dj;
   *I1=*I1>di?*I1:di;
   *J1=*J1>dj?*J1:dj;

   Ref->UnProject(Ref,&di,&dj,Lat1,Lon1,1,1);
   *I0=*I0<di?*I0:di;
   *J0=*J0<dj?*J0:dj;
   *I1=*I1>di?*I1:di;
   *J1=*J1>dj?*J1:dj;

   Ref->UnProject(Ref,&di,&dj,Lat1,Lon0,1,1);
   *I0=*I0<di?*I0:di;
   *J0=*J0<dj?*J0:dj;
   *I1=*I1>di?*I1:di;
   *J1=*J1>dj?*J1:dj;

   *I0=*I0<Ref->X0?Ref->X0:*I0;
   *J0=*J0<Ref->Y0?Ref->Y0:*J0;
   *I1=*I1>Ref->X1?Ref->X1:*I1;
   *J1=*J1>Ref->Y1?Ref->Y1:*J1;

   if (fabs(Lon0-Lon1)>=359.99999999) {
      *I0=Ref->X0;
      *I1=Ref->X1;
   }

   if (fabs(Lat0-Lat1)>=179.99999999) {
      *J0=Ref->Y0;
      *J1=Ref->Y1;
   }
   return(1);
}

/**----------------------------------------------------------------------------
 * @brief  Verifier la validite d'une georeference
 * @remark On projete la bounding box et si les latitudes sont en dehors de -90,90 alors c'est pas bon
 * @author Jean-Philippe Gauthier
 * @date   Fevrier 2009
 *    @param[in]  Ref      Pointeur sur la reference
 *
 *    @return              Booleen indiquant la validite
*/
int GeoRef_Valid(TGeoRef* __restrict const Ref) {

   TCoord co[2];

   if (!Ref) return(0);
   
   Ref->Project(Ref,Ref->X0,Ref->Y0,&co[0].Lat,&co[0].Lon,1,1);
   Ref->Project(Ref,Ref->X1,Ref->Y1,&co[1].Lat,&co[1].Lon,1,1);

   if (co[0].Lat<-91 || co[0].Lat>91.0 || co[1].Lat<-91 || co[1].Lat>91.0) {
      return(0);
   }
   return(1);
}

/**----------------------------------------------------------------------------
 * @brief  Assigner des vecteurs de positions X et Y
 * @author Jean-Philippe Gauthier
 * @date   Mars 2010 
 *    @param[in]  Ref      Pointeur sur la reference
 *    @param[in]  XDef     Data definition des positions en X
 *    @param[in]  YDef     Data definition des positions en X
 *
 *    @return              Dimension des resultats
*/
int GeoRef_Positional(TGeoRef *Ref,TDef *XDef,TDef *YDef) {

   int d,dx,dy,nx,ny;
   
   if (!Ref) return(0);

   /* Check the dimensions */
   nx=FSIZE2D(XDef);
   ny=FSIZE2D(YDef);
   dx=(Ref->X1-Ref->X0+1);
   dy=(Ref->Y1-Ref->Y0+1);

   if (Ref->GRTYP[0]=='Y' || Ref->GRTYP[0]=='W') {
      if (nx!=ny || nx!=dx*dy) {
         return(0);
      }
   } else if (Ref->GRTYP[0]=='V') {
      if (nx!=ny || nx!=dx) {
         return(0);
      }
   } else if (!(Ref->X0 & Ref->X1 & Ref->Y0 & Ref->Y1)) {
       GeoRef_Size(Ref,0,0,nx-1,ny-1,Ref->BD);
       dx=nx;
       dy=ny;
   } else if (nx!=dx || ny!=dy) {
      return(0);
   }

   /*Clear arrays*/
   if (Ref->AX) free(Ref->AX);
   if (Ref->AY) free(Ref->AY);

   Ref->AX=(float*)malloc(nx*sizeof(float));
   Ref->AY=(float*)malloc(ny*sizeof(float));

   if (!Ref->AX || !Ref->AY) {
      return(0);
   }

   /*Assign positionals, if size is float, just memcopy otherwise, assign*/
   if (XDef->Type==TD_Float32) {
      memcpy(Ref->AX,XDef->Data[0],nx*sizeof(float));
   } else {
      for(d=0;d<nx;d++) {
         Def_Get(XDef,0,d,Ref->AX[d]);
      }
   }

   if (YDef->Type==TD_Float32) {
      memcpy(Ref->AY,YDef->Data[0],ny*sizeof(float));
   } else {
      for(d=0;d<ny;d++) {
         Def_Get(YDef,0,d,Ref->AY[d]);
      }
   }

#ifdef HAVE_GDAL
   // Get rid of transforms and projection functions if the positionnale are already in latlon (GDAL case)
   if (Ref->GRTYP[1]=='\0') {
      if (Ref->Transform)    free(Ref->Transform);    Ref->Transform=NULL;
      if (Ref->InvTransform) free(Ref->InvTransform); Ref->InvTransform=NULL;
      if (Ref->Function) {
         OCTDestroyCoordinateTransformation(Ref->Function);
         Ref->Function=NULL;
      }
      if (Ref->InvFunction) {
         OCTDestroyCoordinateTransformation(Ref->InvFunction);
         Ref->InvFunction=NULL;
      }
      OSRDestroySpatialReference(Ref->Spatial);
      Ref->Spatial=NULL;
   }
#endif

   /*Set secondary gridtype to Y for the project/unproject functions to work correctly*/
   if (Ref->GRTYP[0]=='W' && Ref->GRTYP[1]=='\0') {
      Ref->GRTYP[1]=nx>dx?'Y':'Z';
   }

   return(1);
}

/**----------------------------------------------------------------------------
 * @brief  Calculer la position latlon de tous les points de grille.
 * @author Jean-Philippe Gauthier
 * @date   June 2015
 *    @param[in]  Ref     Pointeur sur la reference geographique
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 
 *    @return             Number of coordinates
*/
int GeoRef_Coords(TGeoRef *Ref,double *Lat,double *Lon) {

#ifdef HAVE_RMN
   int x,y,nxy;
   double lat,lon;
   
   if (!Ref) return(0);
   
   nxy=Ref->NX*Ref->NY;
//TODO: Check func ez_calclatlon 
   if (!Ref->Lat) {
      Ref->Lat=(double*)malloc(nxy*sizeof(double));
      Ref->Lon=(double*)malloc(nxy*sizeof(double));
   
      nxy=0;
      for(y=Ref->Y0;y<=Ref->Y1;y++) {
         for(x=Ref->X0;x<=Ref->X1;x++) {
            Ref->Project(Ref,x,y,&lat,&lon,FALSE,TRUE);
            Ref->Lat[nxy]=lat;
            Ref->Lon[nxy]=lon;
            nxy++;
         }
      }      
   }
   
   Lat=Ref->Lat;
   Lon=Ref->Lon;
   return(nxy);
#endif
}

/**----------------------------------------------------------------------------
 * @brief  Obtenir les valeurs de distance en X et Y ainsi que l'aire
 *         pour chaque cellule de la grille
 * @remark Si un des tableau est NULL, il ne sera pas remplie
 * @author Jean-Philippe Gauthier
 * @date   Avril 2010 
 *    @param[in]  Grid     Grille
 *    @param[in]  Invert   Invert (1/area)
 *    @param[out] DX       Valeurs de distance en X
 *    @param[out] DY       Valeurs de distance en y
 *    @param[out] DA       Valeurs de l'aire
 *
 *    @return              Code d'erreur (0=erreur, 1=ok)
*/
int GeoRef_CellDims(TGeoRef *Ref,int Invert,float* DX,float* DY,float* DA) {

   unsigned int i,gi,j,gj,nid,pnid,pidx,idx,*tidx;
   int          ig, nx, ny;
   double       di[4],dj[4],dlat[4],dlon[4];
   double       fx,fy,fz,dx[4],dy[4],s,a,b,c;
   char         grtyp[2];
   TGeoRef *gr;

   if (!Ref || Ref->GRTYP[0]=='X' || Ref->GRTYP[0]=='Y') {
      App_Log(WARNING,"%s: DX, DY and DA cannot be calculated on an X or Y grid\n",__func__);        
      return(FALSE);
   } else if (Ref->GRTYP[0]=='M') {
      
      if (DX || DY) {
         App_Log(WARNING,"%s: DX and DY cannot be calculated on an M grid\n",__func__);        
      }
      
      if (DA) {
         a=(EARTHRADIUS*EARTHRADIUS)*0.5;
         tidx=Ref->Idx;
         for(idx=0;idx<Ref->NIdx-3;idx+=3) {
            
            dx[0]=DEG2RAD(Ref->AX[tidx[idx]]);   dy[0]=DEG2RAD(Ref->AY[tidx[idx]]);
            dx[1]=DEG2RAD(Ref->AX[tidx[idx+1]]); dy[1]=DEG2RAD(Ref->AY[tidx[idx+1]]);
            dx[2]=DEG2RAD(Ref->AX[tidx[idx+2]]); dy[2]=DEG2RAD(Ref->AY[tidx[idx+2]]);
            
            s =(dx[1]-dx[0])*(2+sin(dy[0])+sin(dy[1]));
            s+=(dx[2]-dx[1])*(2+sin(dy[1])+sin(dy[2]));
            s+=(dx[0]-dx[2])*(2+sin(dy[2])+sin(dy[0]));          
            s=fabs(s*a);

//             a =(dx[1]-dx[0])*(2+sin(dy[0])+sin(dy[1]));
//             b =(dx[2]-dx[1])*(2+sin(dy[1])+sin(dy[2]));
//             c =(dx[0]-dx[2])*(2+sin(dy[2])+sin(dy[0]));    
//             s = (a+b+c)*0.5;           
//             s = atan(sqrt(tan(s/2.0)*tan((s-a)/2.0)*tan((s-b)/2.0)*tan((s-c)/2.0)))*0.25*EARTHRADIUS;
            
            // Split area over 3 vertices
            s/=3.0;
            DA[tidx[idx]]+=s;
            DA[tidx[idx+1]]+=s;
            DA[tidx[idx+2]]+=s;    
         }

         if (Invert) {
            for(idx=0;idx<Ref->NX;idx++) DA[idx]=1.0/DA[idx];
         }
      }

   } else {
#ifdef HAVE_RMN            
      pnid=Ref->Options.SubGrid;
      pidx=0;
      nx = Ref->NX;
      ny = Ref->NY;
      
      // Loop on the subgrids if needed
/*       for(nid=(pnid?pnid:(Ref->NbId>1?1:0));nid<=(pnid?pnid:(Ref->NbId>1?Ref->NbId:0));nid++) { */
      for(nid=pnid;nid<=(pnid?pnid:(Ref->NbSub>1?(Ref->NbSub-1):0));nid++) {
         if (Ref->NbSub>1 && !pnid) {
/*             c_ezgprm(Ref->Subs[nid],grtyp,&Ref->NX,&Ref->NY,&ig,&ig,&ig,&ig); */ 
            Ref->NX = Ref->Subs[nid]->NX;
            Ref->NY = Ref->Subs[nid]->NY;
         }

         gr = pnid?Ref->Subs[nid-1]:(Ref->NbSub>1?Ref->Subs[nid]:Ref);

         for(j=0,gj=1;j<Ref->NY;j++,gj++) {
            idx=pidx+j*Ref->NX;
            for(i=0,gi=1;i<Ref->NX;i++,idx++,gi++) {
               
               di[0]=gi-0.5; dj[0]=gj;
               di[1]=gi+0.5; dj[1]=gj;
               di[2]=gi;     dj[2]=gj-0.5;
               di[3]=gi;     dj[3]=gj+0.5;

               // Reproject gridpoint length coordinates of segments crossing center of cell
/*                c_gdllfxy(Ref->Subs[nid],dlat,dlon,di,dj,4); */
               GeoRef_XY2LL(gr,dlat,dlon,di,dj,4);
               dx[0]=DEG2RAD(dlon[0]); dy[0]=DEG2RAD(dlat[0]);
               dx[1]=DEG2RAD(dlon[1]); dy[1]=DEG2RAD(dlat[1]);

               dx[2]=DEG2RAD(dlon[2]); dy[2]=DEG2RAD(dlat[2]);
               dx[3]=DEG2RAD(dlon[3]); dy[3]=DEG2RAD(dlat[3]);

               // Get distance in meters
               fx=DIST(0.0,dy[0],dx[0],dy[1],dx[1]);
               fy=DIST(0.0,dy[2],dx[2],dy[3],dx[3]);

               // If x distance is null, we crossed the pole
               if (fx==0.0)
                  fx=(M_PI*fy)/Ref->NX;

               if (DX) DX[idx]=(Invert?1.0/fx:fx);
               if (DY) DY[idx]=(Invert?1.0/fy:fy);
               if (DA) DA[idx]=(Invert?1.0/(fx*fy):(fx*fy));
            }
         }
         pidx+=idx;
      }
      // Set back original grid
      if (Ref->NbSub>1 && !pnid) {
/*          c_ezgprm(Ref->Subs[pnid],grtyp,&Ref->NX,&Ref->NY,&ig,&ig,&ig,&ig); */
         Ref->NX = nx;
         Ref->NY = ny;
      }
#else
      App_Log(ERROR,"%s: RMNLIB support not included\n",__func__);
#endif
   }
   
   return(TRUE);
}
