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
wordint f77name(gdllfxy)(PTR_AS_INT GRef, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, wordint *n)
{
  return c_gdllfxy((TGeoRef*)GRef, lat, lon, x, y, *n);

}

wordint c_gdllfxy(TGeoRef *GRef, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, wordint n)
{
  wordint j, icode;
  TGeoRef *yin_gd, *yan_gd;
  ftnfloat *latyin, *lonyin, *latyan, *lonyan;
  ftnfloat *tmpy;

  if (GRef->NbSub > 0 )
    {
      yin_gd=GRef->Subs[0];
      yan_gd=GRef->Subs[1];
      tmpy = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      latyin = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      lonyin = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      latyan = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      lonyan = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      for (j=0; j< n; j++)
        {
          if (y[j] > yin_gd->nj)
             {
             tmpy[j]=y[j]-yin_gd->nj;
             }
          else
             {
             tmpy[j]=y[j];
             }
        }
      icode = c_gdllfxy_orig(yin_gd,latyin,lonyin,x,tmpy,n);
/*
      for (j=0; j < n; j++)
        {
           printf("gdllfxy yin x %f y %f : lat %f, lon %f \n",x[j],y[j],lat[j],lon[j]);
        }
*/
      icode = c_gdllfxy_orig(yan_gd,latyan,lonyan,x,tmpy,n);
      for (j=0; j < n; j++)
        {
           if (y[j] > yin_gd->nj)
              {
              lat[j]=latyan[j];
              lon[j]=lonyan[j];
/* printf("gdllfxy yan x %f y %f : lat %f, lon %f \n",x[j],y[j],lat[j],lon[j]); */
              }
           else
              {
              lat[j]=latyin[j];
              lon[j]=lonyin[j];
/* printf("gdllfxy yin x %f y %f : lat %f, lon %f \n",x[j],y[j],lat[j],lon[j]); */
              }
        }
        free(tmpy); free(latyin); free(lonyin); free(latyan); free(lonyan);
    }
  else
    {
      icode = c_gdllfxy_new(GRef,lat,lon,x,y,n);
    }
  return icode;
}

wordint c_gdllfxy_new(TGeoRef *GRef, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, wordint n)
{
  ftnfloat xlat1, xlon1, xlat2, xlon2;
  wordint i, npts, un;
  ftnfloat *tmpx, *tmpy, *ytmp=NULL;
  ftnfloat delxx, delyy;
  ftnfloat dlat, dlon, swlat, swlon;
  wordint indx, indy;

  npts = n;
  un = 1;

  switch(GRef->grtyp[0])
    {
    case 'A':
      for (i=0; i < n; i++)
	      {
	      lat[i] = (y[i]-1.0)*GRef->fst.xg[DLAT]+GRef->fst.xg[SWLAT];
	      lon[i] = (x[i]-1.0)*GRef->fst.xg[DLON]+GRef->fst.xg[SWLON];
	      lon[i] = (ftnfloat) (fmod((double) (lon[i] + 360.0), (double) 360.0));
	      }
      break;

    case 'B':
      for (i=0; i < n; i++)
	      {
	      lat[i] = (y[i]-1.0)*GRef->fst.xg[DLAT]+GRef->fst.xg[SWLAT];
	      lon[i] = (x[i]-1.0)*GRef->fst.xg[DLON]+GRef->fst.xg[SWLON];
	      lon[i] = (ftnfloat) (fmod((double) (lon[i] + 360.0), (double) 360.0));
	      }
      break;

    case 'E':
      tmpx = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      tmpy = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      for (i=0; i < n; i++)
	      {
	      dlat  = 180.0 / GRef->nj;
	      dlon  = 360.0 / (GRef->ni - 1);
	      swlat = -90.0 + 0.5 * dlat;
	      swlon = 0.0;
	      tmpx[i] = (x[i]-1.0)*dlon+swlon;
	      tmpy[i] = (y[i]-1.0)*dlat+swlat;
	      }

      f77name(ez_gfllfxy)(lon,lat,tmpx,tmpy,&n,&GRef->fst.xg[XLAT1],&GRef->fst.xg[XLON1],&GRef->fst.xg[XLAT2],&GRef->fst.xg[XLON2]);
      free(tmpx);
      free(tmpy);
      break;


    case 'L':
      for (i=0; i < n; i++)
	      {
	      lon[i] = (x[i]-1.0)*GRef->fst.xg[DLON]+GRef->fst.xg[SWLON];
	      lon[i] = (ftnfloat) (fmod((double) (lon[i] + 360.0), (double) 360.0));
	      lat[i] = (y[i]-1.0)*GRef->fst.xg[DLAT]+GRef->fst.xg[SWLAT];
	      }
      break;

    case 'N':
    case 'S':
      f77name(ez_vllfxy)(lat,lon,x,y,&npts,&un,&GRef->fst.xg[D60],&GRef->fst.xg[DGRW],&GRef->fst.xg[PI],&GRef->fst.xg[PJ],&GRef->fst.hemisphere);
      for (i=0; i < n; i++)
	      {
	      lon[i] = (ftnfloat) (fmod((double) (lon[i] + 360.0), (double) 360.0));
	      }
      break;

    case 'T':
      f77name(ez_vtllfxy)(lat,lon,x,y, &GRef->fst.xg[CLAT], &GRef->fst.xg[CLON], &GRef->fst.xg[TD60], &GRef->fst.xg[TDGRW], &GRef->ni, &GRef->nj, &npts);
      break;

    case '!':
      f77name(ez_llflamb)(lat,lon,x,y,&npts,&GRef->grtyp,&GRef->fst.ig[IG1], &GRef->fst.ig[IG2], &GRef->fst.ig[IG3], &GRef->fst.ig[IG4],1);
      break;


    case 'Y':
       fprintf(stderr, "********************************************************\n");
       fprintf(stderr, "<gdllfxy>: This operation is not supported for 'Y' grids\n");
       fprintf(stderr, "********************************************************\n");
       break;

    case '#':
    case 'Z':
    case 'G':
      tmpx = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      tmpy = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      ytmp = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      for (i=0; i < n; i++)
	      {
	      indx = (int)x[i]-1;
              ytmp[i] = y[i];
              if (GRef->fst.ig[IG2] == 1)
                  {
                    ytmp[i] = GRef->nj +1.0 - y[i];
                  }
	      indy = (int)ytmp[i]-1;

	      indx = indx < 0 ? 0 : indx;
	      indy = indy < 0 ? 0 : indy;
	      indx = indx > GRef->ni-2 ? GRef->ni-2 : indx;
	      indy = indy > GRef->j2-2 ? GRef->j2-2 : indy;
	      delxx = GRef->ax[indx+1]-GRef->ax[indx];
	      tmpx[i] = GRef->ax[indx] + ((x[i]-1.0-indx)*delxx);

	      delyy = GRef->ay[indy+1]-GRef->ay[indy];
	      tmpy[i] = GRef->ay[indy] + ((ytmp[i]-1.0-indy)*delyy);
	      }

      switch (GRef->grref[0])
	      {
	      case 'E':
	        f77name(cigaxg)(&GRef->grref,&xlat1,&xlon1,&xlat2,&xlon2,
			        &GRef->fst.igref[IG1],&GRef->fst.igref[IG2],&GRef->fst.igref[IG3],&GRef->fst.igref[IG4],1);
	        f77name(ez_gfllfxy)(lon, lat, tmpx, tmpy, &npts, &GRef->fst.xgref[XLAT1], &GRef->fst.xgref[XLON1],
			            &GRef->fst.xgref[XLAT2], &GRef->fst.xgref[XLON2]);
	        break;

	      case 'S':
	      case 'N':
	        f77name(ez_vllfxy)(lat,lon,tmpx,tmpy,&npts,&un,&GRef->fst.xgref[D60],
			           &GRef->fst.xgref[DGRW],&GRef->fst.xgref[PI],&GRef->fst.xgref[PJ],&GRef->fst.hemisphere);
	        for (i=0; i < n; i++)
	          {
	          lon[i] = (ftnfloat) (fmod((double) (lon[i] + 360.0), (double) 360.0));
	          }
	        break;

	      case 'L':
	        for (i=0; i < n; i++)
	          {
	          lat[i] = (tmpy[i])*GRef->fst.xgref[DLAT]+GRef->fst.xgref[SWLAT];
	          lon[i] = (tmpx[i])*GRef->fst.xgref[DLON]+GRef->fst.xgref[SWLON];
	          lon[i] = (ftnfloat) (fmod((double) (lon[i] + 360.0), (double) 360.0));
	          }
	        break;

	      default:
	      fprintf(stderr,"<gdllfxy> Errrrrrrrrrrreur!\n");
	      break;
      	}
      free(tmpx);
      free(tmpy);
      free(ytmp);
      break;
    }

  return 0;

}

wordint c_gdllfxy_orig(TGeoRef *GRef, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, wordint n)
{
  ftnfloat xlat1, xlon1, xlat2, xlon2;
  wordint i, npts, un;
  ftnfloat *tmpx, *tmpy;
  ftnfloat delxx, delyy;
  ftnfloat dlat, dlon, swlat, swlon;
  wordint indx, indy;

  npts = n;
  un = 1;

  switch(GRef->grtyp[0])
    {
    case 'A':
      for (i=0; i < n; i++)
	      {
	      lat[i] = (y[i]-1.0)*GRef->fst.xg[DLAT]+GRef->fst.xg[SWLAT];
	      lon[i] = (x[i]-1.0)*GRef->fst.xg[DLON]+GRef->fst.xg[SWLON];
	      lon[i] = (ftnfloat) (fmod((double) (lon[i] + 360.0), (double) 360.0));
	      }
      break;

    case 'B':
      for (i=0; i < n; i++)
	      {
	      lat[i] = (y[i]-1.0)*GRef->fst.xg[DLAT]+GRef->fst.xg[SWLAT];
	      lon[i] = (x[i]-1.0)*GRef->fst.xg[DLON]+GRef->fst.xg[SWLON];
	      lon[i] = (ftnfloat) (fmod((double) (lon[i] + 360.0), (double) 360.0));
	      }
      break;

    case 'E':
      tmpx = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      tmpy = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      for (i=0; i < n; i++)
	      {
	      dlat  = 180.0 / GRef->nj;
	      dlon  = 360.0 / (GRef->ni - 1);
	      swlat = -90.0 + 0.5 * dlat;
	      swlon = 0.0;
	      tmpx[i] = (x[i]-1.0)*dlon+swlon;
	      tmpy[i] = (y[i]-1.0)*dlat+swlat;
	      }

      f77name(ez_gfllfxy)(lon,lat,tmpx,tmpy,&n,&GRef->fst.xg[XLAT1],&GRef->fst.xg[XLON1],&GRef->fst.xg[XLAT2],&GRef->fst.xg[XLON2]);
      free(tmpx);
      free(tmpy);
      break;


    case 'L':
      for (i=0; i < n; i++)
	      {
	      lon[i] = (x[i]-1.0)*GRef->fst.xg[DLON]+GRef->fst.xg[SWLON];
	      lon[i] = (ftnfloat) (fmod((double) (lon[i] + 360.0), (double) 360.0));
	      lat[i] = (y[i]-1.0)*GRef->fst.xg[DLAT]+GRef->fst.xg[SWLAT];
	      }
      break;

    case 'N':
    case 'S':
      f77name(ez_vllfxy)(lat,lon,x,y,&npts,&un,&GRef->fst.xg[D60],&GRef->fst.xg[DGRW],&GRef->fst.xg[PI],&GRef->fst.xg[PJ],&GRef->fst.hemisphere);
      for (i=0; i < n; i++)
	      {
	      lon[i] = (ftnfloat) (fmod((double) (lon[i] + 360.0), (double) 360.0));
	      }
      break;

    case 'T':
      f77name(ez_vtllfxy)(lat,lon,x,y, &GRef->fst.xg[CLAT], &GRef->fst.xg[CLON], &GRef->fst.xg[TD60], &GRef->fst.xg[TDGRW], &GRef->ni, &GRef->nj, &npts);
      break;

    case '!':
      f77name(ez_llflamb)(lat,lon,x,y,&npts,&GRef->grtyp,&GRef->fst.ig[IG1], &GRef->fst.ig[IG2], &GRef->fst.ig[IG3], &GRef->fst.ig[IG4],1);
      break;


    case 'Y':
       fprintf(stderr, "********************************************************\n");
       fprintf(stderr, "<gdllfxy>: This operation is not supported for 'Y' grids\n");
       fprintf(stderr, "********************************************************\n");
       break;

    case '#':
    case 'Z':
    case 'G':
      tmpx = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      tmpy = (ftnfloat *) malloc(n*sizeof(ftnfloat));
      for (i=0; i < n; i++)
	      {
	      indx = (int)x[i]-1;
	      indy = (int)y[i]-1;

	      indx = indx < 0 ? 0 : indx;
	      indy = indy < 0 ? 0 : indy;
	      indx = indx > GRef->ni-2 ? GRef->ni-2 : indx;
	      indy = indy > GRef->j2-2 ? GRef->j2-2 : indy;
	      delxx = GRef->ax[indx+1]-GRef->ax[indx];
	      tmpx[i] = GRef->ax[indx] + ((x[i]-1.0-indx)*delxx);

	      delyy = GRef->ay[indy+1]-GRef->ay[indy];
	      tmpy[i] = GRef->ay[indy] + ((y[i]-1.0-indy)*delyy);
	      }

      switch (GRef->grref[0])
	      {
	      case 'E':
	        f77name(cigaxg)(&GRef->grref,&xlat1,&xlon1,&xlat2,&xlon2,
			        &GRef->fst.igref[IG1],&GRef->fst.igref[IG2],&GRef->fst.igref[IG3],&GRef->fst.igref[IG4],1);
	        f77name(ez_gfllfxy)(lon, lat, tmpx, tmpy, &npts, &GRef->fst.xgref[XLAT1], &GRef->fst.xgref[XLON1],
			            &GRef->fst.xgref[XLAT2], &GRef->fst.xgref[XLON2]);
	        break;

	      case 'S':
	      case 'N':
	        f77name(ez_vllfxy)(lat,lon,tmpx,tmpy,&npts,&un,&GRef->fst.xgref[D60],
			           &GRef->fst.xgref[DGRW],&GRef->fst.xgref[PI],&GRef->fst.xgref[PJ],&GRef->fst.hemisphere);
	        for (i=0; i < n; i++)
	          {
	          lon[i] = (ftnfloat) (fmod((double) (lon[i] + 360.0), (double) 360.0));
	          }
	        break;

	      case 'L':
	        for (i=0; i < n; i++)
	          {
	          lat[i] = (tmpy[i])*GRef->fst.xgref[DLAT]+GRef->fst.xgref[SWLAT];
	          lon[i] = (tmpx[i])*GRef->fst.xgref[DLON]+GRef->fst.xgref[SWLON];
	          lon[i] = (ftnfloat) (fmod((double) (lon[i] + 360.0), (double) 360.0));
	          }
	        break;

	      default:
	      fprintf(stderr,"<gdllfxy> Errrrrrrrrrrreur!\n");
	      break;
      	}
      free(tmpx);
      free(tmpy);
      break;
    }

  return 0;

}
