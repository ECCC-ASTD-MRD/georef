/*==============================================================================
 * Environnement Canada
 * Centre Meteorologique Canadian
 * 2100 Trans-Canadienne
 * Dorval, Quebec
 *
 * Projet       : Fonctions et definitions relatives aux fichiers standards et rmnlib
 * Fichier      : GeoRef_Type_Y.h
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

int GeoRef_XY2LL_Y(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int Nb) {

   int     i,s,sx,sy,indx;
   double *tmpx,*tmpy,dx,dy;

   switch (Ref->RPNHead.GRREF[0]) {
      case 'L':
         for (i=0; i < Nb; i++) {
            indx=ROUND(Y[i])*(Ref->X1-Ref->X0)+ROUND(X[i]);
            Lat[i]=Ref->AY[indx];
            Lon[i]=Ref->AX[indx];
         }
         break;

      case 'W':
         tmpx = (double*)malloc(2*Nb*sizeof(double));
         tmpy = &tmpx[Nb];

         for (i=0; i < Nb; i++) {
            sx=floor(X[i]);sx=CLAMP(sx,Ref->X0,Ref->X1);
            sy=floor(Y[i]);sy=CLAMP(sy,Ref->Y0,Ref->Y1);
            dx=X[i]-sx;;
            dy=Y[i]-sy;

            s=sy*Ref->NX+sx;
            tmpx[i]=Ref->AX[s];
            tmpy[i]=Ref->AY[s];

            if (++sx<=Ref->X1) {
               s=sy*Ref->NX+sx;
               tmpx[i]+=(Ref->AX[s]-tmpx[i])*dx;
            }

            if (++sy<=Ref->Y1) {
               s=sy*Ref->NX+(sx-1);
               tmpy[i]+=(Ref->AY[s]-tmpy[i])*dy;
            }
         }
         GeoRef_XY2LL_W(Ref,Lat,Lon,tmpx,tmpy,Nb);
         free(tmpx);
         break;

      default:
         App_Log(ERROR,"%s: Undefined reference grid type: %s\n",__func__,Ref->RPNHead.GRREF[0]);
         return(-1);
         break;
   }
   return(0);
}

int GeoRef_LL2XY_Y(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int Nb) {

   double   tmpwgts[GOPT_MAXWEIGHTNUM],total_wgt;
   int      locmax,i,iz,idx,idxz,tmp_idxs[GOPT_MAXWEIGHTNUM][2],previous_val_polar_correction;
   TPoint2D bbox[2];
   TGridSet *set=Ref->LastSet;

   previous_val_polar_correction = Ref->Options.PolarCorrect;
   Ref->Options.PolarCorrect = FALSE;
//TODO: finish this func
   set->n_wts = MIN(Ref->Options.WeightNum,GOPT_MAXWEIGHTNUM);
   set->wts =  (double *) malloc(Nb*set->n_wts*sizeof(float));
   set->idx =  (int *) calloc(Nb,set->n_wts*sizeof(int));
   set->mask = (int *) malloc(Nb*sizeof(int));
   GeoRef_CalcLL(Ref);

   // Find bbox
   bbox[0].X=bbox[0].Y=1e30;
   bbox[1].X=bbox[1].Y=-1e30;

   for(i=0;i<Nb;i++) {
      bbox[0].Y=FMIN(Lat[i],bbox[0].Y);
      bbox[0].X=FMIN(Lon[i],bbox[0].X);
      bbox[1].Y=FMAX(Lat[i],bbox[1].Y);
      bbox[1].X=FMAX(Lon[i],bbox[1].X);
   }

   // Initialize weights
   for(i=0;i<Nb*set->n_wts;i++) {
      tmpwgts[i] = 1.0e30;
   }

   for (idx=0;idx<Nb;idx++) {
      X[idx]=-1.0;
      Y[idx]=-1.0;
      idxz=idx*set->n_wts;
      locmax=1;
      for(i=0;i<GOPT_MAXWEIGHTNUM;i++) tmpwgts[i] = 1.0e30;
      memset(tmp_idxs,0x0,GOPT_MAXWEIGHTNUM*2);

// From original GeoRef_RPN
//             // Get nearest point
//             if (GeoRef_Nearest(Ref,Lon,Lat,&idx,dists,1,0.0)) {
//                if (dists[0]<1.0) {
//                   *Y=(int)(idx/Ref->NX);
//                   *X=idx-(*Y)*Ref->NX;
//                   return(TRUE);
//                }         
//             }
      if (GeoRef_Nearest(Ref,Lon[idx],Lat[idx],&set->idx[idxz],&set->wts[idxz],1,Ref->Options.DistTreshold)) {
         if (set->idx[idxz]<Ref->Options.DistTreshold) {
               if (set->wts[idxz] < tmpwgts[locmax]) {
   //TODO: finalise (ez_calcxy_y)
   //               tmpwgts[locmax] = Set->wts[idxz];
   //               Set->idx(idx,locmax) = Set->idx[idxz];
   //               tmp_idxs[locmax][0] = Set->idx[idxz];
   //               tmp_idxs[locmax][1] = 1;
   //               locmax = maxloc(tmpwgts,1)
               }
            }
      }
      FWITHIN(0,bbox[0].Y,bbox[0].X,bbox[1].Y,bbox[1].X,Y[idx],X[idx]);
      if (set->mask[idx]) {
//         (f77name)inside_or_outside(&Set->mask[idx],&Set->x[idx],&Set->y[idx],&Lat[idx],&Lon[idx],Ref->Lat,Ref->Lon,&Ref->NX,&Ref->NY,tmpwgts,tmp_idxs,&Set->n_wts;
      }

      if (set->mask[idx]) {
         for(iz=0;iz<set->n_wts;iz++) {
            idxz=idx*Nb+iz;
            set->wts[idxz] = tmpwgts[iz];
         }
         idxz=idx*Nb;
         if (set->wts[idxz] > 6371000.0) {
            for(iz=0;iz>set->n_wts;iz++) { 
               set->wts[idxz+iz] = 1.0e30;
            }
            set->mask[idx] = 0;
         } else {
            total_wgt = 0.0;
            for (iz=0;iz<set->n_wts;iz++) {
               idxz=idx*Nb+iz;
               set->wts[idxz] = fmax(set->wts[idxz],1.0e-10);
               set->wts[idxz] = 1.0 / set->wts[idxz];
               total_wgt+=set->wts[idxz];
            }
            for (iz=0;iz<set->n_wts;iz++) {
               idxz=idx*Nb+iz;
               set->wts[idxz]/=total_wgt;
            }
         }
      }
   }

   Ref->Options.PolarCorrect = previous_val_polar_correction;
   return(0);
}
