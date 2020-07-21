/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "ezscint.h"
#include "ez_funcdef.h"
#include "../src/GeoRef.h"

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(ezll)(PTR_AS_INT GRef, ftnfloat *lat, ftnfloat *lon)
{
   wordint icode;
   
   icode = c_gdll((TGeoRef*)GRef, lat, lon);
   return icode;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
wordint f77name(gdll)(PTR_AS_INT GRef, ftnfloat *lat, ftnfloat *lon)
{
   wordint icode;
   
   icode = c_gdll((TGeoRef*)GRef, lat, lon);
   return icode;
}

wordint c_gdll(TGeoRef *GRef, ftnfloat *lat, ftnfloat *lon)
{
   wordint icode;
   wordint ni, nj;
   TGeoRef *yin_gd, *yan_gd;
      
   if (GRef->NbSub > 0 )
      {
      yin_gd = GRef->Subs[0];
      yan_gd = GRef->Subs[1];
   /*    printf("gdll: GRef for yin=%d,GRef for yan=%d\n",yin_gd,yan_gd); */
      ni = yin_gd->ni;
      nj = yin_gd->nj;
   /*  printf("gdll: ni=%d, nj=%d\n",ni,nj); */
      icode=c_gdll_orig(yin_gd,lat,lon);
      icode=c_gdll_orig(yan_gd,&lat[ni*nj],&lon[ni*nj]);
      }
   else
      {
      icode=c_gdll_orig(GRef,lat,lon);
      }
   return icode;
}

wordint c_gdll_orig(TGeoRef *GRef, ftnfloat *lat, ftnfloat *lon)
{
   ez_calclatlon(GRef);
   if (GRef->flags & LAT)
      {
      memcpy(lon, GRef->lon, GRef->ni*GRef->nj*sizeof(ftnfloat));
      if (GRef->fst.axe_y_inverse == 0)
         {
         memcpy(lat, GRef->lat, GRef->ni*GRef->nj*sizeof(ftnfloat));
         }
      else
         {
         memcpy(lat, GRef->lat, GRef->ni*GRef->nj*sizeof(ftnfloat));
         }
      }
   else
      {
      fprintf(stderr, "Erreur! A l'aide! Descripteurs manquants!\n");
      return -1;
      }
   return 0;
}
