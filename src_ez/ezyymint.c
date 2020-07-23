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

wordint c_ezyymint(TGeoRef *gdout, TGeoRef *gdin, wordint ni, wordint nj, ftnfloat *maskout, ftnfloat *dlat, ftnfloat *dlon, ftnfloat *yinlat, ftnfloat *yinlon, wordint *yyincount, ftnfloat *yanlat, ftnfloat *yanlon, wordint *yyancount)
{
  wordint ivalue,icode,i,j,k;
  TGeoRef *yin_mg;
  wordint yincount,yancount,yni,ynj;
  ftnfloat *yin_fld, global_extrap_value, local_extrap_value;
  char interp_degree[32],extrap_degree[32],extrap_value[32],local_val[32];
  char global_interp_degree[32],global_extrap_degree[32];
  
  yin_mg=gdin->mymaskgrid;
  yni=yin_mg->ni;
  ynj=yin_mg->nj;

  yin_fld = (ftnfloat *) malloc(yni*ynj*sizeof(ftnfloat));
  memset(yin_fld,0.0,yni*ynj*sizeof(ftnfloat));
  /*get original options*/
  strcpy(interp_degree,"interp_degree");
  icode = c_ezgetopt(interp_degree,global_interp_degree);
  strcpy(extrap_degree,"extrap_degree");
  strcpy(extrap_value,"extrap_value");
  icode = c_ezgetopt(extrap_degree,global_extrap_degree);
  ivalue = 0;
  if (0 == strcmp(global_extrap_degree,"value"))
  {
    icode = c_ezgetval(extrap_value,&global_extrap_value);
    ivalue = 1;
  }
  strcpy(local_val,"nearest");
  icode = c_ezsetopt(interp_degree, local_val);

  local_extrap_value = 1.0;
  icode = c_ezsetval(extrap_value,local_extrap_value);
  strcpy(local_val,"value");
  icode = c_ezsetopt(extrap_degree, local_val);
  icode = c_ezsint_orig(maskout,yin_fld,gdout,yin_mg);
  /*masking is done,reset original interp options*/
  icode = c_ezsetopt(interp_degree, global_interp_degree);
  if (ivalue == 1)
  {
    icode = c_ezsetval(extrap_value, global_extrap_value);
  }
  icode = c_ezsetopt(extrap_degree, global_extrap_degree);
  free(yin_fld);

/* now create the destination grids */
  yancount=0;
  yincount=0;
  for (j=0; j<nj; j++)
  {
    for (i=0;i<ni; i++)
    {
      k=(j*ni)+i; 
      if (maskout[k] == 1.0)
      {
        yanlat[yancount]=dlat[k];
        yanlon[yancount]=dlon[k];
        yancount++;
      }
      else
      {
        yinlat[yincount]=dlat[k];
        yinlon[yincount]=dlon[k];
        yincount++;
      }
    }
  }
  *yyincount = yincount;
  *yyancount = yancount;
  return icode;
}

