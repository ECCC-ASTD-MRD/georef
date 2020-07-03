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

#include "ez_funcdef.h"
#include "../src/GeoRef.h"

wordint c_ezfreegridset(TGeoRef* gdid, wordint index)
{
   wordint i;

   if (gdid->gset[index].x != NULL)
      {
      free(gdid->gset[index].x);
      gdid->gset[index].x  = NULL;
      }

   if (gdid->gset[index].y)
      {
      free(gdid->gset[index].y);
      gdid->gset[index].y  = NULL;
      }

   if (gdid->gset[index].gemin.lat_rot)
      {
      free(gdid->gset[index].gemin.lat_rot);
      free(gdid->gset[index].gemin.lon_rot);
      free(gdid->gset[index].gemin.sinlat_rot);
      free(gdid->gset[index].gemin.coslat_rot);
      free(gdid->gset[index].gemin.sinlon_rot);
      free(gdid->gset[index].gemin.coslon_rot);
      free(gdid->gset[index].gemin.sinlat_true);
      free(gdid->gset[index].gemin.coslat_true);
      free(gdid->gset[index].gemin.sinlon_true);
      free(gdid->gset[index].gemin.coslon_true);
      memset(&gdid->gset[index].gemin, (int) 0, sizeof(_gemgrid));
      }

   if (gdid->gset[index].gemout.lat_rot)
      {
      free(gdid->gset[index].gemout.lat_rot);
      free(gdid->gset[index].gemout.lon_rot);
      free(gdid->gset[index].gemout.sinlat_rot);
      free(gdid->gset[index].gemout.coslat_rot);
      free(gdid->gset[index].gemout.sinlon_rot);
      free(gdid->gset[index].gemout.coslon_rot);
      free(gdid->gset[index].gemout.sinlat_true);
      free(gdid->gset[index].gemout.coslat_true);
      free(gdid->gset[index].gemout.sinlon_true);
      free(gdid->gset[index].gemout.coslon_true);
      memset(&gdid->gset[index].gemout, (int) 0, sizeof(_gemgrid));
      }

  for (i=0; i < NZONES; i++)
    {
    if (gdid->gset[index].zones[i].npts > 0)
      {
      free(gdid->gset[index].zones[i].x);
      free(gdid->gset[index].zones[i].y);
      free(gdid->gset[index].zones[i].idx);
      gdid->gset[index].zones[i].x = NULL;
      gdid->gset[index].zones[i].y = NULL;
      gdid->gset[index].zones[i].idx = NULL;
      }
    }
  return 0;
}
