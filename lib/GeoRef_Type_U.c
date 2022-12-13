/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Fonctions et definitions relatives aux fichiers standards et rmnlib
 * Fichier      : GeoRef_Type_Z.h
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

/*----------------------------------------------------------------------------
 * @brief  Define a referential of type RPN U
 * @author Jean-Philippe Gauthier
 * @date   Avril 2015
 *    @param[in] NI      Dimension en X
 *    @param[in] NJ      Dimension en Y
 *    @param[in] GRTYP   Grid type ('A', 'B', 'E', 'G', 'L', 'N', 'S','Y', 'Z', '#', '!')
 *    @param[in] GRREF   Reference grid type ('E', 'G', 'L', 'N', 'S')
 *    @param[in] VerCode Version of the grid    
 *    @param[in] NbSub   Number of sub grids
 *    @param[in] Subs    Array of pointers to subgrid
 *
 *    @return            Georef (NULL=error)
*/
TGeoRef* GeoRef_CreateU(int NI,int NJ,char *GRTYP,char *GRREF,int VerCode,int NbSub,TGeoRef **Subs) {

   int  i;
   TGeoRef *ref,*fref,*sub_gd;
    
   if (NbSub <= 1) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: NbSub given is less than 2\n",__func__);
      return(NULL);
   }
   if (VerCode != 1) {
      Lib_Log(APP_LIBGEOREF,APP_ERROR,"%s: Invalid VerCode\n",__func__);
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
      Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Grille[%p].Subs[%p] has maskgrid=%p\n",__func__,ref,Subs[i],sub_gd->mymaskgrid);
   }

   GeoRef_Qualify(ref);

   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: grtyp     = '%c'\n",__func__, ref->GRTYP[0]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: grref     = '%c'\n",__func__, ref->RPNHead.GRREF[0]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: ni        = %d\n",__func__,ref->NX);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: nj        = %d\n",__func__,ref->NY);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: ig1       = %d\n",__func__,ref->RPNHead.IG[X_IG1]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: ig2       = %d\n",__func__,ref->RPNHead.IG[X_IG2]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: ig3       = %d\n",__func__,ref->RPNHead.IG[X_IG3]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: ig4       = %d\n",__func__,ref->RPNHead.IG[X_IG4]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: ig1ref    = %d\n",__func__,ref->RPNHead.IGREF[X_IG1]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: ig2ref    = %d\n",__func__,ref->RPNHead.IGREF[X_IG2]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: ig3ref    = %d\n",__func__,ref->RPNHead.IGREF[X_IG3]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: ig4ref    = %d\n",__func__,ref->RPNHead.IGREF[X_IG4]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: NbSub     = %d\n",__func__,ref->NbSub);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Subs[0]   = %p\n",__func__,ref->Subs[0]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: Subs[1]   = %p\n",__func__,ref->Subs[1]);

   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: RPNHead.XG[1] = %f\n",__func__,ref->RPNHead.XG[1]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: RPNHead.XG[2] = %f\n",__func__,ref->RPNHead.XG[2]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: RPNHead.XG[3] = %f\n",__func__,ref->RPNHead.XG[3]);
   Lib_Log(APP_LIBGEOREF,APP_DEBUG,"%s: RPNHead.XG[4] = %f\n",__func__,ref->RPNHead.XG[4]);

   return(ref);
}
