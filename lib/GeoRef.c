#define GEOREF_BUILD

#include <pthread.h>
#include <App.h>
#include <rmn/List.h>

#include "georef_build_info.h"
#include "GeoRef.h"
#include "Def.h"
#include "Triangle.h"
 
static TList          *GeoRef_List=NULL;                                                                                       ///< Global list of known geo references
static pthread_mutex_t GeoRef_Mutex=PTHREAD_MUTEX_INITIALIZER;                                                                 ///< Thread lock on geo reference access
__thread TGeoOptions   GeoRef_Options= { IR_CUBIC, ER_MAXIMUM, IV_FAST, CB_REPLACE, 0.0, TRUE, FALSE, FALSE, 16, 1, 1, TRUE, FALSE, 10.0, 0.0, 0.0 };  ///< Default options

const char *TRef_InterpVString[] = { "UNDEF","FAST","WITHIN","INTERSECT","CENTROID","ALIASED","CONSERVATIVE","NORMALIZED_CONSERVATIVE","POINT_CONSERVATIVE","LENGTH_CONSERVATIVE","LENGTH_NORMALIZED_CONSERVATIVE","LENGTH_ALIASED",NULL };
const char *TRef_InterpRString[] = { "UNDEF","NEAREST","LINEAR","CUBIC","NORMALIZED_CONSERVATIVE","CONSERVATIVE","MAXIMUM","MINIMUM","SUM","AVERAGE","AVERAGE_VARIANCE","AVERAGE_SQUARE","NORMALIZED_COUNT","COUNT","VECTOR_AVERAGE","NOP","ACCUM","BUFFER","SUBNEAREST","SUBLINEAR",NULL };

__attribute__ ((constructor)) int32_t GeoRef_Init() {
   App_LibRegister(APP_LIBGEOREF,VERSION);
   return(TRUE);
}

/**----------------------------------------------------------------------------
 * @brief  Apply thread lock on GeoRef access
 * @date   February 2008
*/
void GeoRef_Lock() {
   pthread_mutex_lock(&GeoRef_Mutex);
}

/**----------------------------------------------------------------------------
 * @brief  Remove thread lock on GeoRef access
 * @date   February 2008
*/
void GeoRef_Unlock() {
   pthread_mutex_unlock(&GeoRef_Mutex);
}

int32_t GeoRef_Project(struct TGeoRef *Ref,double X,double Y,double *Lat,double *Lon,int32_t Extrap,int32_t Transform) {

   int32_t tr,status;

   if (X<(Ref->X0-0.5) || Y<(Ref->Y0-0.5) || X>(Ref->X1+0.5) || Y>(Ref->Y1+0.5)) {
      if (!Extrap) {
         *Lat=-999.0;
         *Lon=-999.0;
         return(0);
      }
   }

   tr=Ref->Options.Transform;
   Ref->Options.Transform=Transform;
   status=Ref->XY2LL(Ref,Lat,Lon,&X,&Y,1);
   Ref->Options.Transform=tr;

   return(status);
}


int32_t GeoRef_UnProject(struct TGeoRef *Ref,double *X,double *Y,double Lat,double Lon,int32_t Extrap,int32_t Transform) {

   int32_t tr,status;

   if (Lat>90.0 || Lat<-90.0 || Lon==-999.0) 
     return(FALSE);

   tr=Ref->Options.Transform;
   Ref->Options.Transform=Transform;
   status=Ref->LL2XY(Ref,X,Y,&Lat,&Lon,1);
   Ref->Options.Transform=tr;

   if (*X>(Ref->X1+0.5) || *Y>(Ref->Y1+0.5) || *X<(Ref->X0-0.5) || *Y<(Ref->Y0-0.5)) {
      if (!Extrap) {
         *X=-1.0;
         *Y=-1.0;
      }
      return(FALSE);
   }

   return(status);
}

/**----------------------------------------------------------------------------
 * @brief  Initialise the georeference grid limits
 * @date   July 2005
 *    @param[in]  Ref   Pointer to geo reference
 *    @param[in]  X0    X lower limit
 *    @param[in]  Y0    Y lower limit
 *    @param[in]  X1    X higher limit
 *    @param[in]  Y0    Y higher limit
 *    @param[in]  BD    Border width
 */
void GeoRef_Size(TGeoRef *Ref,int32_t X0,int32_t Y0,int32_t X1,int32_t Y1,int32_t BD) {

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
 * @date   July 2005
 *    @param[in]  Ref   Pointer to geo reference
 *
 *    @return           Freed code (0=not freed, 1=freed)
 * 
 */
int32_t GeoRef_Free(TGeoRef *Ref) {

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
 * @date   July 2005
 *    @param[in]  Ref   Pointer to geo reference
 *
 *    @return           New reference count
*/
int32_t GeoRef_Incr(TGeoRef *Ref) {

   if (Ref) {
      return(__sync_add_and_fetch(&Ref->NRef,1));
   } else {
      return(0);
   }
}

/**----------------------------------------------------------------------------
 * @brief  Initialiser la structure App
 * @date   July 2005
 *    @param[in]  Ref   Pointer to geo reference
 *    @param[in]  New   Clear the name associated
*/
void GeoRef_Clear(TGeoRef *Ref,int32_t New) {

   int32_t n;
   
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
      if (Ref->AXY)          free(Ref->AXY);          Ref->AXY=NULL;
      if (Ref->NCX)          free(Ref->NCX);          Ref->NCX=NULL;
      if (Ref->NCY)          free(Ref->NCY);          Ref->NCY=NULL;
      if (Ref->Subs)         free(Ref->Subs);         Ref->Subs=NULL;
      if (Ref->Name)         free(Ref->Name);         Ref->Name=NULL;

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

      memset(&Ref->RPNHead,0x0,sizeof(fst_record_ext));

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
      Ref->XY2LL=NULL;
      Ref->LL2XY=NULL;
      Ref->Height=NULL;
   }
}

/**----------------------------------------------------------------------------
 * @brief  Define qualifying flags for the geo refrence
 * @date   Janvier 2015
 *    @param[in]  Ref   Pointer to geo reference
 */
void GeoRef_Qualify(TGeoRef* __restrict const Ref) {

   TCoord co[2];
   double d[2],lat[2],lon[2],n[2],x[2],y[2];
   int32_t    nx;

   if (Ref) {
      switch(Ref->GRTYP[0]) {
         case 'A': Ref->LL2XY=GeoRef_LL2XY_A; Ref->XY2LL=GeoRef_XY2LL_L; break;
         case 'B': Ref->LL2XY=GeoRef_LL2XY_B; Ref->XY2LL=GeoRef_XY2LL_L; break;
         case 'E': Ref->LL2XY=GeoRef_LL2XY_E; Ref->XY2LL=GeoRef_XY2LL_E; break;
         case 'L': Ref->LL2XY=GeoRef_LL2XY_L; Ref->XY2LL=GeoRef_XY2LL_L; break;
         case 'N': 
         case 'S': Ref->LL2XY=GeoRef_LL2XY_NS; Ref->XY2LL=GeoRef_XY2LL_NS; break;
         case 'T': Ref->LL2XY=GeoRef_LL2XY_T; Ref->XY2LL=GeoRef_XY2LL_T; break;
         case '!': Ref->LL2XY=GeoRef_LL2XY_LAMBERT; Ref->XY2LL=GeoRef_XY2LL_LAMBERT; break;
         case '#':
//TODO: There's something about the G grid         
         case 'G': //GeoRef_LL2XY_G(Ref,X,Y,Lat,Lon,Nb); Ref->XY2LL=GeoRef_XY2LL_G; break;
         case 'Z': Ref->LL2XY=GeoRef_LL2XY_Z; Ref->XY2LL=GeoRef_XY2LL_Z; break;
         case 'O': Ref->LL2XY=GeoRef_LL2XY_O; Ref->XY2LL=GeoRef_XY2LL_O; break;
         case 'R': Ref->LL2XY=GeoRef_LL2XY_R; Ref->XY2LL=GeoRef_XY2LL_R; break;
         case 'M': Ref->LL2XY=GeoRef_LL2XY_M; Ref->XY2LL=GeoRef_XY2LL_M; break;
         case 'W': Ref->LL2XY=GeoRef_LL2XY_W; Ref->XY2LL=GeoRef_XY2LL_W; break;
         case 'Y': Ref->LL2XY=GeoRef_LL2XY_Y; Ref->XY2LL=GeoRef_XY2LL_Y; break;
         default:
            Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid grid type: %c\n",__func__,Ref->GRTYP[0]);
            Ref->Type=GRID_NONE;
            return;
            break;
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
         lat[0]=89.0;lat[1]=-89.0;
         lon[0]=lon[1]=0.0;
         if (GeoRef_LL2XY(Ref,x,y,lat,lon,2,TRUE)) {
            Ref->Type|=GRID_WRAP;
         }
      }  
    
      if (Ref->GRTYP[0]=='A' || Ref->GRTYP[0]=='B' || Ref->GRTYP[0]=='G') {
         Ref->Type|=GRID_WRAP;
      } else if (Ref->GRTYP[0]!='V' && Ref->X0!=Ref->X1 && Ref->Y0!=Ref->Y1) {
         // Check if north is up by looking at longitude variation on an Y increment at grid limits
         x[0]=Ref->X0;x[1]=Ref->X0;
         y[0]=Ref->Y0;y[1]=Ref->Y0+1.0;
         GeoRef_XY2LL(Ref,lat,lon,x,y,2,TRUE);
         d[0]=lon[0]-lon[1];
         x[0]=Ref->X1;x[1]=Ref->X1;
         y[0]=Ref->Y1-1;y[1]=Ref->Y1;
         GeoRef_XY2LL(Ref,lat,lon,x,y,2,TRUE);
         d[1]=lon[0]-lon[1];

         if (fabs(d[0])>0.0001 || fabs(d[1])>0.0001) {
            Ref->Type|=GRID_ROTATED;
         }
                  
         // Get size of a gridpoint
         x[0]=Ref->X0+(Ref->X1-Ref->X0)/2.0;x[1]=x[0]+1.0;
         y[0]=y[1]=Ref->Y0+(Ref->Y1-Ref->Y0)/2.0;
         GeoRef_XY2LL(Ref,lat,lon,x,y,2,TRUE);
         d[0]=DIST(0.0,DEG2RAD(lat[0]),DEG2RAD(lon[0]),DEG2RAD(lat[1]),DEG2RAD(lon[1]));

         // Get distance between first and lat point
         x[0]=Ref->X0;x[1]=Ref->X1;
         y[0]=y[1]=Ref->Y0+(Ref->Y1-Ref->Y0)/2.0;
         GeoRef_XY2LL(Ref,lat,lon,x,y,2,TRUE);
         d[1]=DIST(0.0,DEG2RAD(lat[0]),DEG2RAD(lon[0]),DEG2RAD(lat[1]),DEG2RAD(lon[1]));

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
         for(nx=0;nx<Ref->NX;nx++) { 
            if (Ref->AX[nx]<0) { 
               Ref->Type|=GRID_NEGLON; 
               break; 
            } 
         } 
      }
   }
}         

/**----------------------------------------------------------------------------
 * @brief  Test two GeoRef for equality
 * @date   July 2005
 *    @param[in]  Ref0    First geo reference
 *    @param[in]  Ref1    Second geo reference
 *
 *    @return             Equality (1=True 0=False)
*/
int32_t GeoRef_Equal(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1) {

   if (!Ref0 || !Ref1) {
      return(0);
   }

   if (Ref0->GRTYP[0]!=Ref1->GRTYP[0] || Ref0->GRTYP[1]!=Ref1->GRTYP[1])
      return(0);
   
   if (Ref0->RPNHead.ig1!=Ref1->RPNHead.ig1 || Ref0->RPNHead.ig2!=Ref1->RPNHead.ig2 || Ref0->RPNHead.ig3!=Ref1->RPNHead.ig3 || Ref0->RPNHead.ig4!=Ref1->RPNHead.ig4)
     return(0);
   //TODO: Check AX,AY ?

   // Cloud point32_t should never be tested as equal
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

#ifdef HAVE_GDAL 
   if ((Ref0->Spatial && !Ref1->Spatial) || (!Ref0->Spatial && Ref1->Spatial))
      return(0);

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
 * @brief  Create a new GeoRef but link it to an already existing GeoRef
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
      ref->XY2LL=Ref->XY2LL;
      ref->LL2XY=Ref->LL2XY;
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
 * @date   Janvier 2015
 *    @param[in]  Ref     Pointeur sur la reference
 *
 *    @return             Pointeur sur la copie de la reference
*/
TGeoRef *GeoRef_HardCopy(TGeoRef* __restrict const Ref) {

   TGeoRef *ref;
   int32_t      i;

   ref=GeoRef_New();
   GeoRef_Size(ref,Ref->X0,Ref->Y0,Ref->X1,Ref->Y1,Ref->BD);

   if (Ref) {
      ref->GRTYP[0]=Ref->GRTYP[0];
      ref->GRTYP[1]=Ref->GRTYP[1];
      ref->XY2LL=Ref->XY2LL;
      ref->LL2XY=Ref->LL2XY;
      ref->Type=Ref->Type;
      ref->NbSub=Ref->NbSub;
      ref->QTree=NULL;
      
      if (Ref->Subs) {
         ref->Subs=(TGeoRef**)malloc(Ref->NbSub*sizeof(TGeoRef*));
         memcpy(ref->Subs,Ref->Subs,Ref->NbSub*sizeof(TGeoRef*));
         for(i=0;i<ref->NbSub;i++)
            GeoRef_Incr(ref->Subs[i]);
      }

      memcpy(&ref->RPNHead,&Ref->RPNHead,sizeof(fst_record_ext));
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
            GeoRef_DefineW(ref,Ref->String,Ref->Transform,Ref->InvTransform,Ref->Spatial);
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
 * @date   Janvier 2015
 *    @param[in]  Ref     Pointeur sur la reference
 *    @param[in]  NI      New x size
 *    @param[in]  NJ      New y size
 *
 *    @return             Pointer to the new geo reference
*/
TGeoRef *GeoRef_Resize(TGeoRef* __restrict const Ref,int32_t NI,int32_t NJ) {

   TGeoRef *ref;

   if (!Ref) {
      ref=GeoRef_New();
   } else {
      ref=GeoRef_HardCopy(Ref);
   }
   GeoRef_Size(ref,0,0,NI-1,NJ-1,0);

   return(ref);
}

/**----------------------------------------------------------------------------
 * @brief  Add a geo reference to the list of known geo reference
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
 * @brief  Insert a grid entry into the list of grids managed by ezscint.  Can be used
 *         with regular and irregular ('Y', 'Z') grids, although it is not very useful
 *         for regular grids.
 *         If the grid type corresponds to a regular grid type (eg. 'A', 'G', 'N', etc.),
 *         then the parameters IG1 through IG4 are taken from an ordinary data record
 *         and grref, ax and ay are not used.
 *  
 *         If grtyp == 'Z' or '#', the dimensions of ax=ni and ay=nj.
 *         If grtyp == 'Y', the dimensions of ax=ay=ni*nj. 
 * @date   Avril 2005
 *
 *   @param Ref   Georeference
 *   @param NI    Horizontal size of the grid
 *   @param NJ    Vertical size of the grid
 *   @param GRTYP Grid type ('A', 'B', 'E', 'G', 'L', 'N', 'S','Y', 'Z', '#', '!')
 *   @param grref Reference grid type ('E', 'G', 'L', 'N', 'S')
 *   @param IG1   ig1 value associated to the reference grid
 *   @param IG2   ig2 value associated to the reference grid
 *   @param IG3   ig3 value associated to the reference grid
 *   @param IG4   ig4 value associated to the reference grid
 *   @param AX    Positional axis mapped to the '>>' record
 *   @param AY    Positional axis mapped to the '^^' record
 *  
 *    @return              New geo reference pointer
*/
TGeoRef* GeoRef_Define(TGeoRef *Ref,int32_t NI,int32_t NJ,char* GRTYP,char* grref,int32_t IG1,int32_t IG2,int32_t IG3,int32_t IG4,double* AX,double* AY) {
   
   TGeoRef* ref,*fref;

   ref=Ref?Ref:GeoRef_New();
   if (!ref)
      return(NULL);

   ref->RPNHead.grtyp[0]=ref->GRTYP[0] = GRTYP[0];
   ref->RPNHead.grtyp[1]=ref->GRTYP[1] = '\0';
   ref->RPNHeadExt.grref[0] = grref?grref[0]:'\0';
   ref->RPNHeadExt.grref[1] = '\0';
   ref->RPNHead.ni=ref->NX = NI;
   ref->RPNHead.nj=ref->NY = NJ;
   ref->RPNHead.ig1 = IG1;
   ref->RPNHead.ig2 = IG2;
   ref->RPNHead.ig3 = IG3;
   ref->RPNHead.ig4 = IG4;
   ref->i1 = 1;
   ref->i2 = NI;
   ref->j1 = 1;
   ref->j2 = NJ;
   ref->Type=GRID_NONE;

   switch (GRTYP[0]) {
      case 'Z':
         f77name(cigaxg)(ref->RPNHeadExt.grref,&ref->RPNHeadExt.xgref1,&ref->RPNHeadExt.xgref2,&ref->RPNHeadExt.xgref3,&ref->RPNHeadExt.xgref4,&ref->RPNHeadExt.igref1,&ref->RPNHeadExt.igref2,&ref->RPNHeadExt.igref3,&ref->RPNHeadExt.igref4,1);
      case '#':
      case 'Y':
         ref->AX = AX;
         ref->AY = AY;
         break;
   }

   GeoRef_Size(ref,0,0,NI-1,NJ-1,0);

   if (!Ref) {
      if (fref = GeoRef_Find(ref)) {
         // This georef already exists
         free(ref);
         GeoRef_Incr(fref);
         return(fref);
      }

      // This is a new georef
      GeoRef_Add(ref);
   }

   GeoRef_DefRPNXG(ref);
   GeoRef_AxisDefine(ref,AX,AY);
   GeoRef_Qualify(ref);

   return(ref);
}

int32_t GeoRef_ReadDescriptor(TGeoRef *GRef,void **Ptr,char *Var,int32_t Grid,TApp_Type Type) {

   fst_record *h=&GRef->RPNHead;
   fst_record record=default_fst_record,crit=default_fst_record;
   int32_t ok=0,sz=0,i;

   if (Var && !*Ptr) {
      crit.datev= h->datev;
      crit.ip1  = h->ig1;
      crit.ip2  = h->ig2;
      crit.ip3  = h->ig3;
      strncpy(crit.nomvar,Var,FST_NOMVAR_LEN);

      if (Grid==1) {
         // Look for corresponding time and if not use any time
         if ((ok=fst24_read(h->file,&crit,NULL,&record))<=0) {
            crit.datev= -1;
            ok=fst24_read(h->file,&crit,NULL,&record);
         }
      } else if (Grid==0) {
         strncpy(crit.etiket,h->etiket,FST_ETIKET_LEN);
         strncpy(crit.typvar,h->typvar,FST_TYPVAR_LEN);
         crit.ip1  = h->ip1;
         crit.ip2  = h->ip2;
         crit.ip3  = h->ip3;
         ok=fst24_read(h->file,&crit,NULL,&record);
      } else {
         strncpy(crit.etiket,h->etiket,FST_ETIKET_LEN);
         strncpy(crit.typvar,h->typvar,FST_TYPVAR_LEN);
         crit.ip2  = h->ip2;
         crit.ip3  = h->ip3;
         ok=fst24_read(h->file,&crit,NULL,&record);
      }

      if (ok<=0) {
         return(0);
      } else {
         sz=record.ni*record.nj*record.nk;
         switch(Type) {
            case APP_NIL:
            case APP_FLOAT32:      // Suppose default is float32
               *Ptr=record.data;
               break;

            case APP_FLOAT64:
               if (record.data_bits==64) {
                  *Ptr=record.data;
               } else {
                  if (!(*Ptr=(double*)malloc(sz*sizeof(double)))) {
                     Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Not enough memory to read descriptor field %s\n",__func__,Var);
                     return(0);
                  } 
                  for(i=0;i<sz;i++) ((double*)*Ptr)[i]=((float*)record.data)[i]; 
                  fst24_record_free(&record);
               }
               break;

            case APP_UINT32:     // Suppose default as int
               if (record.data_bits==32 && record.data_type==FST_TYPE_UNSIGNED) {
                  *Ptr=record.data;
               } else {
                  if (!(*Ptr=(unsigned int*)malloc(sz*sizeof(unsigned int)))) {
                     Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Not enough memory to read descriptor field\n",__func__,Var);
                     return(0);
                  } 
                  for(i=0;i<sz;i++) ((unsigned int*)*Ptr)[i]=((int*)record.data)[i];
                  fst24_record_free(&record);
               }
               break;

            case APP_INT32:      // Suppose default as unsigned int
               if (record.data_bits==32 && record.data_type==FST_TYPE_SIGNED) {
                  *Ptr=record.data;
               } else {
                  if (!(*Ptr=(int*)malloc(sz*sizeof(int)))) {
                     Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Not enough memory to read descriptor field\n",__func__,Var);
                     return(0);
                  } 
                  for(i=0;i<sz;i++) ((int*)*Ptr)[i]=((int*)record.data)[i];
                  fst24_record_free(&record);
               }
               break;

            default:
               *Ptr=record.data;
                break;
         }
      }

      // Define reference grid parameters
      strncpy(GRef->RPNHeadExt.grref,record.grtyp,FST_GTYP_LEN);
      GRef->RPNHeadExt.igref1=record.ig1;
      GRef->RPNHeadExt.igref2=record.ig2;
      GRef->RPNHeadExt.igref3=record.ig3;
      GRef->RPNHeadExt.igref4=record.ig4;
   }
   
   return(sz);
}


int32_t GeoRef_Read(struct TGeoRef *GRef) {

   int32_t     key,ni,nj,nk,ig1,ig2,ig3,ig4,idx,s,i,j,offsetx,offsety,sz;
   float      *ax=NULL,*ay=NULL;
   char        grref[2];

   if (!GRef) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid GeoRef object\n",__func__);
      return(FALSE);
   }

   if (GRef->GRTYP[0]=='L' || GRef->GRTYP[0]=='A' || GRef->GRTYP[0]=='B' || GRef->GRTYP[0]=='N' || GRef->GRTYP[0]=='S' || GRef->GRTYP[0]=='G') {
      return(TRUE);
   }

   if ((!GRef->AY || !GRef->AX) && GRef->RPNHead.file) {

      switch(GRef->GRTYP[0]) {
         case 'M':
            if (!GRef->Idx) GRef->NIdx=GeoRef_ReadDescriptor(GRef,(void **)&GRef->Idx,"##",1,APP_UINT32);
            if (!GRef->AY)  GeoRef_ReadDescriptor(GRef,(void **)&GRef->AY,"^^",1,APP_FLOAT64);
            if (!GRef->AX)  sz=GeoRef_ReadDescriptor(GRef,(void **)&GRef->AX,">>",1,APP_FLOAT64);
            GeoRef_BuildIndex(GRef);
            break;

         case 'W':
            break;

         case 'Y':
            if (!GRef->AY) GeoRef_ReadDescriptor(GRef,(void **)&GRef->AY,"LA",0,APP_FLOAT64);
            if (!GRef->AY) GeoRef_ReadDescriptor(GRef,(void **)&GRef->AY,"^^",1,APP_FLOAT64);
            if (!GRef->AX) sz=GeoRef_ReadDescriptor(GRef,(void **)&GRef->AX,"LO",0,APP_FLOAT64);
            if (!GRef->AX) sz=GeoRef_ReadDescriptor(GRef,(void **)&GRef->AX,">>",1,APP_FLOAT64);
            if (!GRef->Hgt) GeoRef_ReadDescriptor(GRef,(void **)&GRef->Hgt,"ZH",0,APP_FLOAT32);

           break;

         case 'X':
         case 'O':
            if (!GRef->AY) GeoRef_ReadDescriptor(GRef,(void **)&GRef->AY,"^^",1,APP_FLOAT64);
            if (!GRef->AX) sz=GeoRef_ReadDescriptor(GRef,(void **)&GRef->AX,">>",1,APP_FLOAT64);
            GeoRef_BuildIndex(GRef);          
            break;

         case 'V':
            if (!GRef->AY) GeoRef_ReadDescriptor(GRef,(void **)&GRef->AY,"^^",1,APP_FLOAT64);
            if (!GRef->AX) sz=GeoRef_ReadDescriptor(GRef,(void **)&GRef->AX,">>",1,APP_FLOAT64);
//TODO:            RPN_FieldReadLevels(Field);
            break;

         case 'Z':
            if (!GRef->AY) GeoRef_ReadDescriptor(GRef,(void **)&GRef->AY,"^^",1,APP_FLOAT64);
            if (!GRef->AX) sz=GeoRef_ReadDescriptor(GRef,(void **)&GRef->AX,">>",1,APP_FLOAT64);
            break;

         case '#':
            if (!GRef->AY) GeoRef_ReadDescriptor(GRef,(void **)&GRef->AY,"^^",1,APP_FLOAT64);
            if (!GRef->AX) sz=GeoRef_ReadDescriptor(GRef,(void **)&GRef->AX,">>",1,APP_FLOAT64);
   // TODO: Check with LireEnrPositionel
   //        if (ax && ay) {
   //            GRef->AX = (double*) malloc(GRef->NX*sizeof(double));
   //            GRef->AY = (double*) malloc(GRef->NY*sizeof(double));
   //            offsetx = ip3 - 1;
   //            offsety = ip4 - 1;
   //            for (j=0; j < GRef->NY; j++) GRef->AY[j] = ay[j+offsety];
   //            for (i=0; i < GRef->NX; i++) GRef->AX[i] = ax[i+offsetx];
   //         }
            break;
         
         case 'U':
            if (!GRef->AXY) GeoRef_ReadDescriptor(GRef,(void **)&GRef->AXY,"^>",1,APP_FLOAT64);

            if (GRef->AXY) {
               GRef->NbSub=(int)ax[2];            // Number of LAM grids (YY=2)
               ni=(int)GRef->AXY[5];              // NI size of LAM grid 
               nj=(int)GRef->AXY[6];              // NJ size of LAM grid
               GRef->NX = ni;                     // NI size of U grid
               GRef->NY = nj*GRef->NbSub;         // NJ size of U grid
               GRef->AX = (double*)malloc(ni*sizeof(double));
               GRef->AY = (double*)malloc(nj*sizeof(double));
               for(i=0;i<ni;i++) GRef->AX[i]=GRef->AXY[15+i];
               for(i=0;i<nj;i++) GRef->AY[i]=GRef->AXY[15+ni+i];

               // Get subgrids
               GRef->Subs = (TGeoRef**)malloc(GRef->NbSub*sizeof(TGeoRef*));
               strcpy(grref,"E");
               idx=11;
               float xg1,xg2,xg3,xg4;

               for(s=0;s<GRef->NbSub;s++) {
                  xg1=GRef->AXY[idx];
                  xg2=GRef->AXY[idx+1];
                  xg3=GRef->AXY[idx+2];
                  xg4=GRef->AXY[idx+3];
                  f77name(cxgaig)(grref,&ig1,&ig2,&ig3,&ig4,&xg1,&xg2,&xg3,&xg4,1);
                  //TODO: chekc AX,AY indexes
                  GRef->Subs[s] = GeoRef_Define(NULL,ni,nj,"Z",grref,ig1,ig2,ig3,ig4,GRef->AX,GRef->AY);
                  //TODO: Do we need to do this here ?
                  GeoRef_MaskYYDefine(GRef->Subs[s]);
                  idx+=ni+nj+10;
               }
            }
      }

      // Check for 2D AX/AY
      if (sz>GRef->NX) {
         GRef->Type|=GRID_AXY2D;
      }
      // In case of WKT REF
      if (GRef->GRTYP[0]=='W' || GRef->RPNHeadExt.grref[0] == 'W') {
#ifdef HAVE_GDAL
         double *mtx=NULL,inv[6],*tm=NULL,*im=NULL;
         char  *proj=NULL;

         if (!GeoRef_ReadDescriptor(GRef,(void **)&proj,"PROJ",1,APP_NIL)) {
            Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not find projection description field PROJ (c_fstinf failed)\n",__func__);
            return(FALSE);
         }
         if (!GeoRef_ReadDescriptor(GRef,(void **)&mtx,"MTRX",1,APP_FLOAT64)) {
            Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not find transform matrix field MTRX (c_fstinf failed)\n",__func__);
            return(FALSE);
         } else {
            tm=mtx;
            if (!GDALInvGeoTransform(mtx,inv)) {
               im=NULL;
            } else {
               im=inv;
            }
         }
         GRef->RPNHeadExt.grref[0]='W';
         GeoRef_DefineW(GRef,proj,tm,im,NULL);
         if (proj) free(proj);
         if (mtx)  free(mtx);

#else
   Lib_Log(APP_LIBGEOREF,APP_ERROR,"W grid support not enabled, needs to be built with GDAL\n",__func__);
   return(FALSE);
#endif
      }
   }

   return(TRUE);
}

/**----------------------------------------------------------------------------
 * @brief  Create new geo reference
 * @date   Avril 2005
 *
 *   @param NI    Horizontal size of the grid
 *   @param NJ    Vertical size of the grid
 *   @param GRTYP Grid type ('A', 'B', 'E', 'G', 'L', 'N', 'S','Y', 'Z', '#', '!')
 *   @param IG1   ig1 value associated to the reference grid
 *   @param IG2   ig2 value associated to the reference grid
 *   @param IG3   ig3 value associated to the reference grid
 *   @param IG4   ig4 value associated to the reference grid
 *   @param FID   FSTD file identifier where to look for grid descriptors
 * 
 *   @return              New geo reference pointer
*/
TGeoRef* GeoRef_Create(int32_t NI,int32_t NJ,char *GRTYP,int32_t IG1,int32_t IG2,int32_t IG3,int32_t IG4,fst_file *File) {

   TGeoRef *ref,*fref;
   int32_t      id;

   ref=GeoRef_New();

   // If not specified, type is X
   if (GRTYP[0]==' ') GRTYP[0]='X';
   
   GeoRef_Size(ref,0,0,NI-1,NJ-1,0);
   ref->RPNHead.grtyp[0]=ref->GRTYP[0]=GRTYP[0];
   ref->RPNHead.grtyp[1]=ref->GRTYP[1]=GRTYP[1];
   ref->RPNHead.ip1 = 0;
   ref->RPNHead.ip2 = 0;
   ref->RPNHead.ip3 = 0;
   ref->Type=GRID_NONE;

   if ((NI>1 || NJ>1) && GRTYP[0]!='X' && GRTYP[0]!='P' && GRTYP[0]!='V' && ((GRTYP[0]!='Z' && GRTYP[0]!='Y') || File)) {

      if (GRTYP[1]=='#') {
         //TODO: CHECK For tiled grids (#) we have to fudge the IG3 ang IG4 to 0 since they're used for tile limit
      }

      if (ref->GRTYP[0]=='L' || ref->GRTYP[0]=='A' || ref->GRTYP[0]=='B' || ref->GRTYP[0]=='N' || ref->GRTYP[0]=='S' || ref->GRTYP[0]=='G') {
         // No need to look for grid descriptors
         return(GeoRef_Define(NULL,NI,NJ,GRTYP," ",IG1,IG2,IG3,IG4,NULL,NULL));
      }
  
      ref->RPNHead.file= File;
      ref->RPNHead.ig1 = IG1;
      ref->RPNHead.ig2 = IG2;
      ref->RPNHead.ig3 = IG3;
//      ref->RPNHead.ig4 = (GRTYP[0]=='#' || GRTYP[0]=='U')?IG4:0;
      ref->RPNHead.ig4 = IG4;

      // This georef already exists
      if (fref=GeoRef_Find(ref)) {
         free(ref);
         GeoRef_Incr(fref);
         return(fref);
      }

      // This is a new georef
      GeoRef_Add(ref);
      if (!GeoRef_Read(ref)) {
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
            GeoRef_DefRPNXG(ref);           
            GeoRef_AxisCalcNewtonCoeff(ref);
         }
      }
   }

   GeoRef_CalcLL(ref);
   GeoRef_Qualify(ref);

   return(ref);
}

/**----------------------------------------------------------------------------
 * @brief  Create new geo reference
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
   ref->Sub=0;
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
   ref->AXY=NULL;
   ref->NCX=NULL;
   ref->NCY=NULL;
   ref->RefFrom=NULL;
   ref->QTree=NULL;
   ref->GRTYP[0]='X';
   ref->GRTYP[1]='\0';
   ref->GRTYP[2]='\0';
   ref->Sets=NULL;
   ref->LastSet=NULL;
   ref->NbSet=0;

   // Assign default options
   GeoRef_Options.NoData=nan("");
   memcpy(&ref->Options,&GeoRef_Options,sizeof(TGeoOptions));

   // RPN Specific
   memset(&ref->RPNHead,0x0,sizeof(fst_record));

   // WKT Specific
   ref->String=NULL;
   ref->Spatial=NULL;
   ref->Function=0;
   ref->InvFunction=0;
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
   ref->XY2LL=NULL;
   ref->LL2XY=NULL;
   ref->Height=NULL;

   pthread_mutex_init(&ref->Mutex, NULL);
   return(ref);
}

/**----------------------------------------------------------------------------
 * @brief  Create spatial index
 * @date   Janvier 2016
 *    @param[in]  Ref     Pointeur sur la reference
 *
 *    @return             Quad tree spatial index
*/
TQTree* GeoRef_BuildIndex(TGeoRef* __restrict const Ref) {

   uint32_t  n,x,y,t;
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
            Lib_Log(APP_LIBGEOREF,APP_WARNING,"%s: Failed to allocate baricentric weight array\n",__func__);
         }
      }
      
      // Create the tree on the data limits
      if (!(Ref->QTree=QTree_New(lon0,lat0,lon1,lat1,NULL))) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Failed to create QTree index\n",__func__);
         return(NULL);
      }

      // Loop on triangles
      for(n=0,t=0;n<Ref->NIdx;n+=3,t++) {          
         tr[0][0]=Ref->AX[Ref->Idx[n]];     tr[0][1]=Ref->AY[Ref->Idx[n]];
         tr[1][0]=Ref->AX[Ref->Idx[n+1]];   tr[1][1]=Ref->AY[Ref->Idx[n+1]];
         tr[2][0]=Ref->AX[Ref->Idx[n+2]];   tr[2][1]=Ref->AY[Ref->Idx[n+2]];
         
         // Calculate barycentric weight
          if (Ref->Wght)
             Ref->Wght[t]=1.0/((tr[1][0]-tr[0][0])*(tr[2][1]-tr[0][1])-(tr[2][0]-tr[0][0])*(tr[1][1]-tr[0][1]));
 
         // Put it in the quadtree, in any child nodes intersected and set false pointer increment (+1)
         if (!QTree_AddTriangle(Ref->QTree,tr,GRID_MQTREEDEPTH,(void*)(n+1))) {
            Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Failed to add node\n",__func__);
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
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Failed to create QTree index\n",__func__);
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
int32_t GeoRef_Nearest(TGeoRef* __restrict const Ref,double X,double Y,int32_t *Idxs,double *Dists,int32_t NbNear,double MaxDist) {

   double       dx,dy,l;
   uint32_t n,nn,nr,nnear;
   TQTree      *node;
   int32_t          dxy,x,y,xd,yd,rx;
  
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
 * @brief  Check for intersection of 2 geo reference
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
int32_t GeoRef_Intersect(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1,int32_t *X0,int32_t *Y0,int32_t *X1,int32_t *Y1,int32_t BD) {

   double la,lo,di,dj,in=0;
   double x0,y0,x1,y1,dx,dy;
   int32_t    x,y;

   if (!Ref0 || !Ref1) return(0);

   // Source grid Y
   if (Ref1->GRTYP[0]=='Y') {
      *X0=Ref1->X0; *Y0=Ref1->Y0;
      *X1=Ref1->X1; *Y1=Ref1->Y1;
      return(1);
   }

   // If destination is global
   if (Ref0->Type&GRID_WRAP) {
      *X0=Ref1->X0; *Y0=Ref1->Y0;
      *X1=Ref1->X1; *Y1=Ref1->Y1;
      in=1;
   }

   // Test for limit source inclusion into destination
   x0=y0=1e32;
   x1=y1=-1e32;
   x=0;

   if (!in) {
      dx=Ref1->X0; dy=Ref1->Y0;
      GeoRef_XY2LL(Ref1,&la,&lo,&dx,&dy,1,TRUE);
      if (GeoRef_LL2XY(Ref0,&di,&dj,&la,&lo,1,FALSE)) {
         x0=Ref1->X0;y0=Ref1->Y0;
         x++;
      }
      dx=Ref1->X0; dy=Ref1->Y1;
      GeoRef_XY2LL(Ref1,&la,&lo,&dx,&dy,1,TRUE);
      if (GeoRef_LL2XY(Ref0,&di,&dj,&la,&lo,1,FALSE)) {
         x0=Ref1->X0;y1=Ref1->Y1;
         x++;
      }
      dx=Ref1->X1; dy=Ref1->Y0;
      GeoRef_XY2LL(Ref1,&la,&lo,&dx,&dy,1,TRUE);
      if (GeoRef_LL2XY(Ref0,&di,&dj,&la,&lo,1,FALSE)) {
         x1=Ref1->X1;y0=Ref1->Y0;
         x++;
      }
      dx=Ref1->X1; dy=Ref1->Y1;
      GeoRef_XY2LL(Ref1,&la,&lo,&dx,&dy,1,TRUE);
      if (GeoRef_LL2XY(Ref0,&di,&dj,&la,&lo,1,FALSE)) {
         x1=Ref1->X1;y1=Ref1->Y1;
      }
      *X0=x0; *Y0=y0;
      *X1=x1; *Y1=y1;

      if (x>=3) {
         in=1;
      }
   }

   if (!in) {

      // Project Ref0 within Ref1 and get limits
      for(x=Ref0->X0;x<=Ref0->X1;x++) {
         dx=x;dy=Ref0->Y0;
         GeoRef_XY2LL(Ref1,&la,&lo,&dx,&dy,1,TRUE);
         GeoRef_LL2XY(Ref0,&di,&dj,&la,&lo,1,FALSE);
         x0=fmin(x0,di); y0=fmin(y0,dj);
         x1=fmax(x1,di); y1=fmax(y1,dj);

         dx=x;dy=Ref0->Y1;
         GeoRef_XY2LL(Ref1,&la,&lo,&dx,&dy,1,TRUE);
         GeoRef_LL2XY(Ref0,&di,&dj,&la,&lo,1,FALSE);
         x0=fmin(x0,di); y0=fmin(y0,dj);
         x1=fmax(x1,di); y1=fmax(y1,dj);
      }

      for(y=Ref0->Y0;y<=Ref0->Y1;y++) {
         dx=Ref0->X0; dy=y;
         GeoRef_XY2LL(Ref1,&la,&lo,&dx,&dy,1,TRUE);
         GeoRef_LL2XY(Ref0,&di,&dj,&la,&lo,1,FALSE);
         x0=fmin(x0,di); y0=fmin(y0,dj);
         x1=fmax(x1,di); y1=fmax(y1,dj);

         dx=Ref0->X1; dy=y;
         GeoRef_XY2LL(Ref1,&la,&lo,&dx,&dy,1,TRUE);
         GeoRef_LL2XY(Ref0,&di,&dj,&la,&lo,1,FALSE);
         x0=fmin(x0,di); y0=fmin(y0,dj);
         x1=fmax(x1,di); y1=fmax(y1,dj);
      }

      // Test for north and south pole including grid
      dx=0.0;dy=89.9;
      if (GeoRef_LL2XY(Ref0,&di,&dj,&dy,&dx,1,FALSE) && dj>Ref0->Y0+2 && dj<Ref0->Y1-2 && di>Ref0->X0+2 && di<Ref0->X1-2) {
         GeoRef_LL2XY(Ref1,&di,&dj,&dy,&dx,1,TRUE);
         x0=fmin(x0,di); y0=fmin(y0,dj);
         x1=fmax(x1,di); y1=fmax(y1,dj);
      }
      dx=0.0;dy=-89.9;
      if (GeoRef_LL2XY(Ref0,&di,&dj,&dy,&dx,1,FALSE) && dj>Ref0->Y0+2 && dj<Ref0->Y1-2 && di>Ref0->X0+2 && di<Ref0->X1-2) {
         GeoRef_LL2XY(Ref1,&di,&dj,&dy,&dx,1,TRUE);
         x0=fmin(x0,di); y0=fmin(y0,dj);
         x1=fmax(x1,di); y1=fmax(y1,dj);
      }

      *X0=floor(x0); *Y0=floor(y0);
      *X1=ceil(x1);  *Y1=ceil(y1);

      if (!VOUT(*X0,*X1,Ref1->X0,Ref1->X1) && !VOUT(*Y0,*Y1,Ref1->Y0,Ref1->Y1)) {
         in=1;
      }
   }

   // Clamp the coordinates
   if (BD) {
      REF_CLAMP(Ref1,*X0,*Y0,*X1,*Y1);
   } else {
      REF_CLAMPBD(Ref1,*X0,*Y0,*X1,*Y1);
   }

   return(in);
}

/**----------------------------------------------------------------------------
 * @brief  Calculates geographical limits of a geo reference
 * @date   Aout 2006
 *    @param[in]  Ref      Pointeur sur la reference
 *    @param[out] Lat0     Latitude of first corner
 *    @param[out] Lon0     Longitude of first corner
 *    @param[out] Lat1     Latitude of second corner
 *    @param[out] Lon1     Longitude of second corner
 *
 *    @return              Error code
*/
int32_t GeoRef_Limits(TGeoRef* __restrict const Ref,double *Lat0,double *Lon0,double *Lat1,double *Lon1) {

   int32_t x,y;
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
      di=x;dj=y;
      GeoRef_XY2LL(Ref,&lat,&lon,&di,&dj,1,FALSE);
      *Lat0=fmin(*Lat0,lat); *Lon0=fmin(*Lon0,lon);
      *Lat1=fmax(*Lat1,lat); *Lon1=fmax(*Lon1,lon);
   }

   for(x=Ref->X0,y=Ref->Y1;x<=Ref->X1;x++) {
      di=x;dj=y;
      GeoRef_XY2LL(Ref,&lat,&lon,&di,&dj,1,FALSE);
      *Lat0=fmin(*Lat0,lat); *Lon0=fmin(*Lon0,lon);
      *Lat1=fmax(*Lat1,lat); *Lon1=fmax(*Lon1,lon);
   }

   for(y=Ref->Y0,x=Ref->X0;y<=Ref->Y1;y++) {
      di=x;dj=y;
      GeoRef_XY2LL(Ref,&lat,&lon,&di,&dj,1,FALSE);
      *Lat0=fmin(*Lat0,lat); *Lon0=fmin(*Lon0,lon);
      *Lat1=fmax(*Lat1,lat); *Lon1=fmax(*Lon1,lon);
   }

   for(y=Ref->Y0,x=Ref->X1;y<=Ref->Y1;y++) {
      di=x;dj=y;
      GeoRef_XY2LL(Ref,&lat,&lon,&di,&dj,1,FALSE);
      *Lat0=fmin(*Lat0,lat); *Lon0=fmin(*Lon0,lon);
      *Lat1=fmax(*Lat1,lat); *Lon1=fmax(*Lon1,lon);
   }

   // Test for north and south pole including grid
   lat=90.0;lon=0.0;
   if (GeoRef_LL2XY(Ref,&di,&dj,&lat,&lon,1,FALSE) && dj>Ref->Y0+2 && dj<Ref->Y1-2 && di>Ref->X0+2 && di<Ref->X1-2) {
      *Lat1=90.0;
   }
   lat=-90.0;lon=0.0;
   if (GeoRef_LL2XY(Ref,&di,&dj,&lat,&lon,1,FALSE) && dj>Ref->Y0+2 && dj<Ref->Y1-2 && di>Ref->X0+2 && di<Ref->X1-2) {
      *Lat0=-90.0;
   }
   return(1);
}

/**----------------------------------------------------------------------------
 * @brief  Check if a georef is within another
 * @date   February 2009
 *    @param[in]  Ref0     Pointeur sur la reference inclusive
 *    @param[in]  Ref1     Pointeur sur la reference a tester
  *
 *    @return              True or False
*/
int32_t GeoRef_Within(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1) {

   double la,lo,x0,y0,x1,y1;
   int32_t    x,y;

   if (!Ref0 || !Ref1) return(0);
   
   // Project Ref0 Border within Ref1 and get limits
   for(x0=Ref0->X0,y0=Ref0->Y0;x0<=Ref0->X1;x0+=1.0) {
      GeoRef_XY2LL(Ref0,&la,&lo,&x0,&y0,1,TRUE);
      if (!GeoRef_LL2XY(Ref1,&x1,&y1,&la,&lo,1,FALSE)) {
         return(0);
      }
   }

   for(x0=Ref0->X0,y0=Ref0->Y1;x0<=Ref0->X1;x0+=1.0) {
      GeoRef_XY2LL(Ref0,&la,&lo,&x0,&y0,1,TRUE);
      if (!GeoRef_LL2XY(Ref1,&x1,&y1,&la,&lo,1,FALSE)) {
         return(0);
      }
   }

   for(y0=Ref0->Y0,x0=Ref0->X0;y0<=Ref0->Y1;y0+=1.0) {
      GeoRef_XY2LL(Ref0,&la,&lo,&x0,&y0,1,TRUE);
      if (!GeoRef_LL2XY(Ref1,&x1,&y1,&la,&lo,1,FALSE)) {
         return(0);
      }
   }

   for(y0=Ref0->Y0,x0=Ref0->X1;y0<=Ref0->Y1;y+=1.0) {
      GeoRef_XY2LL(Ref0,&la,&lo,&x0,&y0,1,TRUE);
      if (!GeoRef_LL2XY(Ref1,&x1,&y1,&la,&lo,1,FALSE)) {
         return(0);
      }
   }
   return(1);
}
                    
/**----------------------------------------------------------------------------
 * @brief  Check if a geographic range is within a georef
 * @date   February 2009
 *    @param[in]  Ref      Pointeur sur la reference
 *    @param[in]  Lat0     Latitude of first corner
 *    @param[in]  Lon0     Longitude of first corner
 *    @param[in]  Lat1     Latitude of second corner
 *    @param[in]  Lon1     Longitude of second corner
 *    @param[in]  In       Check for all inside or outside
 * 
 *    @return              True or False
*/
int32_t GeoRef_WithinRange(TGeoRef* __restrict const Ref,double Lat0,double Lon0,double Lat1,double Lon1,int32_t In) {

   double lat[4],lon[4],x[4],y[4],dl;
   int32_t    d0,d1,d2,d3;

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
   
   // Check image within range
   x[0]=Ref->X0; y[0]=Ref->Y0;
   x[1]=Ref->X1; y[1]=Ref->Y0;
   x[2]=Ref->X1; y[2]=Ref->Y1;
   x[3]=Ref->X0; y[3]=Ref->Y1;
   GeoRef_XY2LL(Ref,lat,lon,x,y,4,FALSE);

   d0=FWITHIN(dl,Lat0,Lon0,Lat1,Lon1,lat[0],lon[0]);
   if (!In && d0) return(1);

   d1=FWITHIN(dl,Lat0,Lon0,Lat1,Lon1,lat[1],lon[1]);
   if (!In && d1) return(1);

   d2=FWITHIN(dl,Lat0,Lon0,Lat1,Lon1,lat[2],lon[2]);
   if (!In && d2) return(1);

   d3=FWITHIN(dl,Lat0,Lon0,Lat1,Lon1,lat[3],lon[3]);
   if (!In && d3) return(1);

   // Check for all contained
   if (In) {
      if (d0 && d1 && d2 && d3) {
         return(1);
      } else {
         return(0);
      }
   }

   // Check range within image
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

int32_t GeoRef_WithinCell(TGeoRef *Ref,Vect2d Pos,Vect2d Pt[4],int32_t Idx0,int32_t Idx1,int32_t Idx2,int32_t Idx3) {
 
   Vect3d b;
   int32_t    t0,t1,sz=Ref->NX*Ref->NY;
   
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
int32_t GeoRef_BoundingBox(TGeoRef* __restrict const Ref,double Lat0,double Lon0,double Lat1,double Lon1,double *I0,double *J0,double *I1,double *J1) {

   double di,dj;

   if (!Ref) return(0);
   
   *I0=*J0=1000000;
   *I1=*J1=-1000000;

   GeoRef_LL2XY(Ref,&di,&dj,&Lat0,&Lon0,1,TRUE);
   *I0=*I0<di?*I0:di;
   *J0=*J0<dj?*J0:dj;
   *I1=*I1>di?*I1:di;
   *J1=*J1>dj?*J1:dj;

   GeoRef_LL2XY(Ref,&di,&dj,&Lat0,&Lon1,1,TRUE);
   *I0=*I0<di?*I0:di;
   *J0=*J0<dj?*J0:dj;
   *I1=*I1>di?*I1:di;
   *J1=*J1>dj?*J1:dj;

   GeoRef_LL2XY(Ref,&di,&dj,&Lat1,&Lon1,1,TRUE);
   *I0=*I0<di?*I0:di;
   *J0=*J0<dj?*J0:dj;
   *I1=*I1>di?*I1:di;
   *J1=*J1>dj?*J1:dj;

   GeoRef_LL2XY(Ref,&di,&dj,&Lat1,&Lon0,1,TRUE);
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
 * @date   Fevrier 2009
 *    @param[in]  Ref      Pointeur sur la reference
 *
 *    @return              Booleen indiquant la validite
*/
int32_t GeoRef_Valid(TGeoRef* __restrict const Ref) {

   double x[2],y[2],lat[2],lon[2];

   if (!Ref) return(0);
   
   x[0]=Ref->X0;y[0]=Ref->Y0;
   x[1]=Ref->X1;y[1]=Ref->Y1;
   GeoRef_XY2LL(Ref,lat,lon,x,y,2,TRUE);

   if (lat[0]<-91 || lat[0]>91.0 || lat[1]<-91 || lat[1]>91.0) {
      return(0);
   }
   return(1);
}

/**----------------------------------------------------------------------------
 * @brief  Assigner des vecteurs de positions X et Y
 * @date   Mars 2010 
 *    @param[in]  Ref      Pointeur sur la reference
 *    @param[in]  XDef     Data definition des positions en X
 *    @param[in]  YDef     Data definition des positions en X
 *
 *    @return              Dimension des resultats
*/
int32_t GeoRef_Positional(TGeoRef *Ref,TDef *XDef,TDef *YDef) {

   int32_t d,dx,dy,nx,ny;
   
   if (!Ref) return(0);

   // Check the dimensions
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

   // Clear arrays
   if (Ref->AX) free(Ref->AX);
   if (Ref->AY) free(Ref->AY);

   Ref->AX=(double*)malloc(nx*sizeof(double));
   Ref->AY=(double*)malloc(ny*sizeof(double));

   if (!Ref->AX || !Ref->AY) {
      return(0);
   }

   // Assign positionals, if size is float, just memcopy otherwise, assign
   if (XDef->Type==TD_Float64) {
      memcpy(Ref->AX,XDef->Data[0],nx*sizeof(float));
   } else {
      for(d=0;d<nx;d++) {
         Def_Get(XDef,0,d,Ref->AX[d]);
      }
   }

   if (YDef->Type==TD_Float64) {
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
 * @brief  Obtenir les valeurs de distance en X et Y ainsi que l'aire
 *         pour chaque cellule de la grille
 * @remark Si un des tableau est NULL, il ne sera pas remplie
 * @date   Avril 2010 
 *    @param[in]  Grid     Grille
 *    @param[in]  Invert   Invert (1/area)
 *    @param[out] DX       Valeurs de distance en X
 *    @param[out] DY       Valeurs de distance en y
 *    @param[out] DA       Valeurs de l'aire
 *
 *    @return              Code d'erreur (0=erreur, 1=ok)
*/
int32_t GeoRef_CellDims(TGeoRef *Ref,int32_t Invert,float* DX,float* DY,float* DA) {

   uint32_t i,gi,j,gj,nid,pnid,pidx,idx,*tidx;
   int32_t          ig, nx, ny;
   double       di[4],dj[4],dlat[4],dlon[4];
   double       fx,fy,fz,dx[4],dy[4],s,a,b,c;
   char         grtyp[2];
   TGeoRef *gr;

   if (!Ref || Ref->GRTYP[0]=='X' || Ref->GRTYP[0]=='Y') {
      Lib_Log(APP_LIBGEOREF,APP_WARNING,"%s: DX, DY and DA cannot be calculated on an X or Y grid\n",__func__);        
      return(FALSE);
   } else if (Ref->GRTYP[0]=='M') {
      
      if (DX || DY) {
         Lib_Log(APP_LIBGEOREF,APP_WARNING,"%s: DX and DY cannot be calculated on an M grid\n",__func__);        
      }
      
      if (DA) {
         a=(EARTHRADIUS*EARTHRADIUS)*0.5;
         tidx=Ref->Idx;
         for(idx=0;idx<Ref->NIdx;idx+=3) {
            
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
      pnid=Ref->Sub;
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

               // Reproject gridpoint32_t length coordinates of segments crossing center of cell
/*                c_gdllfxy(Ref->Subs[nid],dlat,dlon,di,dj,4); */
               GeoRef_XY2LL(gr,dlat,dlon,di,dj,4,TRUE);
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
   }
   
   return(TRUE);
}

int32_t GeoRef_DefRPNXG(TGeoRef* Ref) {

   switch (Ref->GRTYP[0]) {
      case 'A':
      case 'G':
         Ref->RPNHeadExt.xg4  = 360. /Ref->NX;
         Ref->RPNHeadExt.xg2 = 0.0;
         switch (Ref->RPNHead.ig1) {
				case 0:
				   Ref->RPNHeadExt.xg3 = 180./Ref->NY;
				   Ref->RPNHeadExt.xg1 = -90. + 0.5*Ref->RPNHeadExt.xg3;
				   break;

				case 1:
				   Ref->RPNHeadExt.xg3 = 90./Ref->NY;
				   Ref->RPNHeadExt.xg1 = 0.5*Ref->RPNHeadExt.xg3;
				   Ref->Type |= GRID_EXPAND;
				   break;

				case 2:
				   Ref->RPNHeadExt.xg3 = 90./Ref->NY;
				   Ref->RPNHeadExt.xg1 = -90. + 0.5*Ref->RPNHeadExt.xg3;
				   Ref->Type |= GRID_EXPAND;
				   break;

				default:
			      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: 'A' grid has to be Global/North/South\n",__func__);
               return(-1);
				   break;
			}

         switch(Ref->RPNHead.ig2) {
	         case 1:
	            Ref->Type|=GRID_YINVERT;
	            break;

	         default:
	            break;
	      }
         break;

      case 'B':
         Ref->RPNHeadExt.xg4 = 360. /(Ref->NX-1);
         Ref->RPNHeadExt.xg2 = 0.0;
         switch (Ref->RPNHead.ig1) {
	         case 0:
	            Ref->RPNHeadExt.xg3 = 180./(Ref->NY-1);
	            Ref->RPNHeadExt.xg1 = -90.;
	            break;

	         case 1:
	            Ref->RPNHeadExt.xg3 = 90./(Ref->NY-1);
	            Ref->RPNHeadExt.xg1 = 0.;
	            Ref->Type |= GRID_EXPAND;
	            break;

	         case 2:
	            Ref->RPNHeadExt.xg3 = 90./(Ref->NY-1);
	            Ref->RPNHeadExt.xg1 = -90.;
	            Ref->Type |= GRID_EXPAND;
	            break;

	         default:
  			      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: 'B' grid has to be Global/North/South\n",__func__);
	            return(-1);
	      }

         switch(Ref->RPNHead.ig2) {
	         case 1:
	            Ref->Type|=GRID_YINVERT;
	            break;

	         default:
	            break;
	      }
         break;

      case 'E':
         f77name(cigaxg)(Ref->GRTYP,&Ref->RPNHeadExt.xg1,&Ref->RPNHeadExt.xg2,&Ref->RPNHeadExt.xg3,&Ref->RPNHeadExt.xg4,&Ref->RPNHead.ig1,&Ref->RPNHead.ig2,&Ref->RPNHead.ig3,&Ref->RPNHead.ig4,1);
      /*      Ref->RPNHeadExt.xg3 = 180./Ref->NY;
	      Ref->RPNHeadExt.xg4 = 360./(Ref->NX-1);
	      Ref->RPNHeadExt.xg2 = 0.0;
	      Ref->RPNHeadExt.xg1 = -90. + 0.5*Ref->RPNHeadExt.xg3;
      */
         break;

      case 'H':
      case 'Y':
      case '!':
         break;

      case '#':
      case 'Z':
         if (Ref->RPNHeadExt.grref[0] == 'N') Ref->Hemi = GRID_NORTH;
         if (Ref->RPNHeadExt.grref[0] == 'S') Ref->Hemi = GRID_SOUTH;
         if (Ref->RPNHeadExt.grref[0] == 'E' || Ref->RPNHeadExt.grref[0]== 'L') {
            f77name(cigaxg)(Ref->RPNHeadExt.grref,&Ref->RPNHeadExt.xgref1, &Ref->RPNHeadExt.xgref2, &Ref->RPNHeadExt.xgref3, &Ref->RPNHeadExt.xgref4,&Ref->RPNHeadExt.igref1, &Ref->RPNHeadExt.igref2, &Ref->RPNHeadExt.igref3, &Ref->RPNHeadExt.igref4,1);
         }
         break;

      case 'L':
         f77name(cigaxg)(Ref->GRTYP,&Ref->RPNHeadExt.xg1, &Ref->RPNHeadExt.xg2, &Ref->RPNHeadExt.xg3, &Ref->RPNHeadExt.xg4,&Ref->RPNHead.ig1, &Ref->RPNHead.ig2, &Ref->RPNHead.ig3, &Ref->RPNHead.ig4,1);
         break;

      case 'N':
         f77name(cigaxg)(Ref->GRTYP,&Ref->RPNHeadExt.xg1, &Ref->RPNHeadExt.xg2, &Ref->RPNHeadExt.xg3, &Ref->RPNHeadExt.xg4,&Ref->RPNHead.ig1, &Ref->RPNHead.ig2, &Ref->RPNHead.ig3, &Ref->RPNHead.ig4,1);
         Ref->Hemi = GRID_NORTH;
         break;

      case 'S':
         f77name(cigaxg)(Ref->GRTYP,&Ref->RPNHeadExt.xg1, &Ref->RPNHeadExt.xg2, &Ref->RPNHeadExt.xg3, &Ref->RPNHeadExt.xg4,&Ref->RPNHead.ig1, &Ref->RPNHead.ig2, &Ref->RPNHead.ig3, &Ref->RPNHead.ig4,1);
         Ref->Hemi = GRID_SOUTH;
         break;

      case 'T':
		   //TODO: What's T
         f77name(cigaxg)(Ref->GRTYP,&Ref->RPNHeadExt.xg1, &Ref->RPNHeadExt.xg2, &Ref->RPNHeadExt.xg3, &Ref->RPNHeadExt.xg4,&Ref->RPNHead.ig1, &Ref->RPNHead.ig2, &Ref->RPNHead.ig3, &Ref->RPNHead.ig4,1);
         break;

      default:
	      Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Grid type not supported %c\n",__func__,Ref->GRTYP[0]);
         return(-1);
    }

   return(0);
}

int32_t GeoRef_GridGetParams(TGeoRef *Ref,int32_t *NI,int32_t *NJ,char *GRTYP,int32_t *IG1,int32_t *IG2,int32_t *IG3,int32_t *IG4,char *grref,int32_t *IG1REF,int32_t *IG2REF,int32_t *IG3REF,int32_t *IG4REF) {
   
   *NI     = Ref->RPNHead.ni;
   *NJ     = Ref->RPNHead.nj;

   GRTYP[0]  = Ref->RPNHead.grtyp[0];
   GRTYP[1]  = '\0';
   grref[0]  = Ref->RPNHeadExt.grref[0];
   grref[1]  = '\0';
  
   *IG1    = Ref->RPNHead.ig1;
   *IG2    = Ref->RPNHead.ig2;
   *IG3    = Ref->RPNHead.ig3;
   *IG4    = Ref->RPNHead.ig4;
   *IG1REF = Ref->RPNHeadExt.igref1;
   *IG2REF = Ref->RPNHeadExt.igref2;
   *IG3REF = Ref->RPNHeadExt.igref3;
   *IG4REF = Ref->RPNHeadExt.igref4;

   return(0);
}

/*----------------------------------------------------------------------------
 * @brief  Writes a grid definition
 * @date   January 2020
 *    @param[in]  GRef      GeoRef pointer
 *    @param[in]  File      FSTD file pointer
 *
 *    @return             Error code (0=ok)
*/
int32_t GeoRef_Write(TGeoRef *GRef,fst_file *File){

   fst_record record=default_fst_record;
   int32_t i,dbl=FALSE;
   char *c;

   if (!GRef) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid georef (NULL)\n",__func__);
      return(FALSE);
   }
   if (!File) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid file (NULL)\n",__func__);
      return(FALSE);
   }

   if (!GRef->Name) GRef->Name=strdup("Undefined");

   if ((c=getenv("GEOREF_DESCRIPTOR_64"))) {
      dbl=TRUE;
   }
   record.data = (float*)calloc(GRef->NX*GRef->NY,sizeof(float));
   record.ni   = GRef->NX;
   record.nj   = GRef->NY;
   record.dateo = 0;
   record.deet  = 0;
   record.npas  = 0;
   record.ip1  = 0;
   record.ip2  = 0;
   record.ip3  = 0;
   record.ip3  = 0;
   strncpy(record.typvar,"X",FST_TYPVAR_LEN);
   strncpy(record.nomvar,"GRID",FST_NOMVAR_LEN);
   strncpy(record.grtyp,GRef->RPNHead.grtyp,FST_GTYP_LEN);
   strncpy(record.etiket,GRef->Name,FST_ETIKET_LEN);
   record.ig1   = GRef->RPNHead.ig1;
   record.ig2   = GRef->RPNHead.ig2;
   record.ig3   = GRef->RPNHead.ig3;
   record.ig4   = GRef->RPNHead.ig4;
   record.data_type = FST_TYPE_REAL_IEEE;
   record.data_bits = 32;
   record.pack_bits = 32;
   if (fst24_write(File,&record,FST_SKIP)<=0) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write grid record field (fst24_write failed)\n",__func__);
   }

   if (dbl)
      free(record.data);

   if (GRef->AXY) {
      record.data_type = FST_TYPE_REAL_IEEE;
      record.ni   = GRef->NX*GRef->NY*GRef->NbSub+15;
      record.nj   = 1;
      record.nk   = 1;
      record.ip1  = GRef->RPNHead.ig1;
      record.ip2  = GRef->RPNHead.ig2;
      record.ip3  = GRef->RPNHead.ig3;
      strncpy(record.nomvar,"^>",FST_NOMVAR_LEN);
      strncpy(record.grtyp,GRef->RPNHeadExt.grref,FST_GTYP_LEN);
      record.ig1   = GRef->RPNHeadExt.igref1;
      record.ig2   = GRef->RPNHeadExt.igref2;
      record.ig3   = GRef->RPNHeadExt.igref3;
      record.ig4   = GRef->RPNHeadExt.igref4;
      if (dbl) {
         record.data = GRef->AXY;
         record.pack_bits = 64;
         record.data_bits = 64;
      } else {
         for(i=0;i<(record.ni*record.nj);i++) ((float*)record.data)[i]=GRef->AXY[i];
         record.pack_bits = 32;
         record.data_bits = 32;
      }
      if (fst24_write(File,&record,FST_SKIP)<=0) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write >> field (fst24_write failed)\n",__func__);
         return(FALSE);
      } 
   } else if (GRef->AX && GRef->AY) {
      record.data_type = FST_TYPE_REAL_IEEE;
      record.ni   = GRef->NX;
      record.nj   = GRef->Type&GRID_AXY2D?GRef->NY:1;
      record.nk   = 1;
      record.ip1  = GRef->RPNHead.ig1;
      record.ip2  = GRef->RPNHead.ig2;
      record.ip3  = GRef->RPNHead.ig3;
      strncpy(record.nomvar,">>",FST_NOMVAR_LEN);
      strncpy(record.grtyp,GRef->RPNHeadExt.grref,FST_GTYP_LEN);
      record.ig1   = GRef->RPNHeadExt.igref1;
      record.ig2   = GRef->RPNHeadExt.igref2;
      record.ig3   = GRef->RPNHeadExt.igref3;
      record.ig4   = GRef->RPNHeadExt.igref4;
      if (dbl) {
         record.data = GRef->AX;
         record.pack_bits = 64;
         record.data_bits = 64;
      } else {
         for(i=0;i<(GRef->Type&GRID_AXY2D?record.ni*record.nj:record.ni);i++) ((float*)record.data)[i]=GRef->AX[i];
         record.pack_bits = 32;
         record.data_bits = 32;
      }
      if (fst24_write(File,&record,FST_SKIP)<=0) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write >> field (fst24_write failed)\n",__func__);
         return(FALSE);
      }

      strncpy(record.nomvar,"^^",FST_NOMVAR_LEN);
      record.ni   = GRef->Type&GRID_AXY2D?GRef->NX:1;
      record.nj   = GRef->GRTYP[0]=='M'?GRef->NX:GRef->NY;
      if (dbl) {
         record.data = GRef->AY;
      } else {
         for(i=0;i<(GRef->Type&GRID_AXY2D?record.ni*record.nj:record.nj);i++) ((float*)record.data)[i]=GRef->AY[i];
      }
      if (fst24_write(File,&record,FST_SKIP)<=0) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write >> field (fst24_write failed)\n",__func__);
         return(FALSE);
      }
   }

   if (GRef->Idx) {
      record.data = GRef->Idx;
      record.pack_bits = 32;
      record.data_bits = 32;
      record.data_type = FST_TYPE_UNSIGNED;
      record.ni   = GRef->NIdx;
      record.nj   = 1;
      record.nk   = 1;
      record.ip1  = GRef->RPNHead.ig1;
      record.ip2  = GRef->RPNHead.ig2;
      record.ip3  = GRef->RPNHead.ig3;
      strncpy(record.nomvar,"##",FST_NOMVAR_LEN);
      strncpy(record.grtyp,"X",FST_GTYP_LEN);
      record.ig1   = GRef->RPNHeadExt.igref1;
      record.ig2   = GRef->RPNHeadExt.igref2;
      record.ig3   = GRef->RPNHeadExt.igref3;
      record.ig4   = GRef->RPNHeadExt.igref4;
      if (fst24_write(File,&record,FST_SKIP)<=0) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write ## field (fst24_write failed)\n",__func__);
         return(FALSE);
      }
   }

   if (GRef->Transform) {
      record.data = GRef->Transform;
      record.pack_bits = 64;
      record.data_bits = 64;
      record.data_type = FST_TYPE_REAL_IEEE;
      record.ni   = 6;
      record.nj   = 1;
      record.nk   = 1;
      record.ip1  = GRef->RPNHead.ig1;
      record.ip2  = GRef->RPNHead.ig2;
      record.ip3  = GRef->RPNHead.ig3;
      strncpy(record.nomvar,"MTRX",FST_NOMVAR_LEN);
      strncpy(record.grtyp,"X",FST_GTYP_LEN);
      record.ig1   = GRef->RPNHeadExt.igref1;
      record.ig2   = GRef->RPNHeadExt.igref2;
      record.ig3   = GRef->RPNHeadExt.igref3;
      record.ig4   = GRef->RPNHeadExt.igref4;
      if (fst24_write(File,&record,FST_SKIP)<=0) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write MTRX field (fst24_write failed)\n",__func__);
         return(FALSE);
      }
   }

   if (GRef->String) {
      record.data = GRef->String;
      record.pack_bits = 8;
      record.data_bits = 8;
      record.data_type = FST_TYPE_UNSIGNED;
      record.ni   = strlen(GRef->String);
      record.nj   = 1;
      record.nk   = 1;
      record.ip1  = GRef->RPNHead.ig1;
      record.ip2  = GRef->RPNHead.ig2;
      record.ip3  = GRef->RPNHead.ig3;
      strncpy(record.nomvar,"PROJ",FST_NOMVAR_LEN);
      strncpy(record.grtyp,"X",FST_GTYP_LEN);
      record.ig1   = GRef->RPNHeadExt.igref1;
      record.ig2   = GRef->RPNHeadExt.igref2;
      record.ig3   = GRef->RPNHeadExt.igref3;
      record.ig4   = GRef->RPNHeadExt.igref4;
      if (fst24_write(File,&record,FST_SKIP)<=0) {
         Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Could not write PROJ field (fst24_write failed)\n",__func__);
         return(FALSE);
      }
   }
   
   return(TRUE);
}

int32_t GeoRef_CopyDesc(fst_file *FileTo,fst_record* Rec) {

   fst_record  srec = default_fst_record;
   fst_record  rec;
   fst_query  *query;
   char       *data=NULL;
   const char *desc,**descs;
   int32_t         d=0,ni,nj,nk,sz=0,ip1,ip2;
   int32_t         key;

   if (Rec->file) {

      // Loop through possible descriptors NOMVAR
      descs=fst24_record_get_descriptors();
      while((desc=descs[d++])) {
         if (strncmp(desc,"HY  ",FST_NOMVAR_LEN)!=0) {
            srec.ip1=Rec->ig1;
            srec.ip2=Rec->ig2;
         }

         // Does it already exists in the destination file
         strncmp(srec.nomvar,desc,FST_NOMVAR_LEN);
         query = fst24_new_query(FileTo,&srec,NULL);
         if (fst24_find_next(query,&rec)) {
            // If not already existing in destination
            if (fst24_read(Rec->file,&srec,NULL,&rec)) {
               if (!fst24_write(FileTo,&rec,TRUE)) {
                  return(FALSE);
               }
            }
         }
      }

      fst24_record_free(&rec);
   } else {
      return(FALSE);
   }

   return(TRUE);
}