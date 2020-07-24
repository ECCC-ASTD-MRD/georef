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

#ifdef MUTEX
//JP
extern pthread_mutex_t EZ_MTX;
#endif

/* TODO: mutex still needed? */
void EliminerGrille(TGeoRef* GRef)
{
   wordint i, index;
    
#ifdef MUTEX
// JP
   pthread_mutex_lock(&EZ_MTX);
#endif

   if (GRef->flags & LAT)
   {
      free(GRef->lat);
      free(GRef->lon);
      GRef->lat = NULL;
      GRef->lon = NULL;
   }

   if (GRef->flags & EZ_AX)
   {
      free(GRef->ax);
      free(GRef->ay);
      GRef->ax = NULL;
      GRef->ay = NULL;
   }

   if (GRef->ncx != NULL)
   {
      free(GRef->ncx);
      free(GRef->ncy);
      GRef->ncx = NULL;
      GRef->ncy = NULL;
   }
   GRef->flags = (int)0;
   
   // Release ezscint sub-grid
   if (Ref->Subs) {
      free(Ref->Subs);  Ref->Subs=NULL;
   }

   for (i=0; i < GRef->n_gdin_for; i++)
   {
      index = ez_find_gdin_in_gset(GRef, GRef->gdin_for[i]);
      c_ezfreegridset(GRef->gdin_for[i], index);
   }

#ifdef MUTEX
// JP
   pthread_mutex_unlock(&EZ_MTX);
#endif

   }
