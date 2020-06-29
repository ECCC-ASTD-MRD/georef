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
wordint ez_calclatlon(TGeoRef* GRef)
   {
   ftnfloat xlat00, xlon00, dlat, dlon;
   wordint i,j,k,ni, nj, npts, hemisphere, gdrow, gdcol;
   ftnfloat *lonp, *latp, *xp, *yp;
   ftnfloat *x, *y;


/*    c_gdkey2rowcol(gdid, &gdrow, &gdcol); */
   if (!(GRef->flags & LAT))
      {
      ni = GRef->ni;
      nj = GRef->nj;
      npts = ni*nj;

      GRef->lat = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
      GRef->lon = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));

      switch(GRef->grtyp[0])
         {
         case 'A':
         case 'B':
         f77name(grll)(GRef->lat,GRef->lon,&ni,&nj,
		    &GRef->fst.xg[SWLAT],&GRef->fst.xg[SWLON], &GRef->fst.xg[DLAT], &GRef->fst.xg[DLON]);
         break;

         case 'E':
         dlon = 360. /(ni-1);
         dlat = 180./(nj);
         xlon00 = 0.0;
         xlat00 = -90. + 0.5*dlat;

         f77name(grll)(GRef->lat,GRef->lon,&ni,&nj,&xlat00,&xlon00,&dlat,&dlon);

         f77name(cigaxg)(&GRef->grtyp, &GRef->fst.xg[XLAT1], &GRef->fst.xg[XLON1],
			&GRef->fst.xg[XLAT2], &GRef->fst.xg[XLON2],
         &GRef->fst.ig[IG1],  &GRef->fst.ig[IG2], &GRef->fst.ig[IG3], &GRef->fst.ig[IG4],1);
         latp = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
         lonp = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
         f77name(ez_gfllfxy)(lonp,latp,GRef->lon,GRef->lat,&npts,
			       &GRef->fst.xg[XLAT1], &GRef->fst.xg[XLON1], &GRef->fst.xg[XLAT2],
			       &GRef->fst.xg[XLON2]);
         memcpy(GRef->lat,latp,npts*sizeof(ftnfloat));
         memcpy(GRef->lon,lonp,npts*sizeof(ftnfloat));
         free(latp);
         free(lonp);
         break;

         case 'L':
         f77name(grll)(GRef->lat,GRef->lon,&ni,&nj,
			 &GRef->fst.xg[SWLAT],&GRef->fst.xg[SWLON],
			 &GRef->fst.xg[DLAT], &GRef->fst.xg[DLON]);
           break;

         case 'N':
         case 'S':
         if (GRef->grtyp[0] == 'N')
	        {
	        hemisphere = 1;
	        }
         else
           {
           hemisphere = 2;
            }
         f77name(grps)(GRef->lat,GRef->lon,&ni,&nj,
			 &GRef->fst.xg[PI],&GRef->fst.xg[PJ],
			 &GRef->fst.xg[D60], &GRef->fst.xg[DGRW], &hemisphere);
         break;

         case 'T':
         npts = ni * nj;
         latp = (ftnfloat *) malloc(npts*sizeof(ftnfloat));
         lonp = (ftnfloat *) malloc(npts*sizeof(ftnfloat));
         xp = (ftnfloat *) malloc(npts*sizeof(ftnfloat));
         yp = (ftnfloat *) malloc(npts*sizeof(ftnfloat));
         for (j=0; j < nj; j++)
            {
            for (i=0; i < ni; i++)
              {
              k = j*ni + i;
              yp[k] = 1.0 * (j+1);
              xp[k] = 1.0 * (i+1);
              }
            }

          f77name(ez_vtllfxy)(latp,lonp, xp,yp,
                           &GRef->fst.xg[CLAT], &GRef->fst.xg[CLON],
                           &GRef->fst.xg[TD60],&GRef->fst.xg[TDGRW],
                           &ni,&nj,&npts);

          memcpy(GRef->lon, lonp, ni*nj*sizeof(ftnfloat));
          memcpy(GRef->lat, latp, ni*nj*sizeof(ftnfloat));
          free(lonp);
          free(latp);
          free(xp);
          free(yp);
          break;

         case 'Y':
         switch (GRef->grref[0])
            {
            case 'N':
            case 'S':
            fprintf(stderr, "<ez_calclatlon> Operation not supported - Y grid on PS Grid\n");
            return -1;
            break;

            case 'L':
            case 'O':
            memcpy(GRef->lon, GRef->ax, GRef->ni*GRef->nj*sizeof(ftnfloat));
	    memcpy(GRef->lat, GRef->ay, GRef->ni*GRef->nj*sizeof(ftnfloat));
	    for (i=0; i < GRef->ni*GRef->nj; i++)
               {
	       if (GRef->lon[i] < 0.0)
                  {
	          GRef->lon[i] = GRef->ax[i] + 360.0;
	          }
	       }
	    break;

	    case 'E':
	    fprintf(stderr, "<ez_calclatlon> Operation not supported - Y grid on E Grid\n");
	    return -1;
            break;

	    }
         break;

         case '#':
         case 'Z':
           case 'G':
             for (j=0; j < nj; j++)
               {
               for (i=0; i < ni; i++)
                 {
                 GRef->lat[C_TO_FTN(i,j,ni)] = GRef->ay[j];
                 GRef->lon[C_TO_FTN(i,j,ni)] = GRef->ax[i];
                 }
               }

	   if (GRef->grtyp[0] == 'G' && GRef->fst.ig[IG1] == NORD)
	     {
	     for (j=0; j < nj; j++)
	       {
	       for (i=0; i < ni; i++)
		 {
		 GRef->lat[C_TO_FTN(i,j,ni)] = GRef->ay[j+nj];
		 }
	       }
	     }

           switch (GRef->grref[0])
	     {
	     case 'N':
	     case 'S':
	       latp = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
	       lonp = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
	       f77name(ez_vllfxy)(latp,lonp,
			       GRef->lon,GRef->lat,&ni,&nj,
			       &GRef->fst.xgref[D60],&GRef->fst.xgref[DGRW],
			       &GRef->fst.xgref[PI], &GRef->fst.xgref[PJ], &GRef->fst.hemisphere);

	       for (i=0; i < ni*nj; i++)
		 {
		 if (lonp[i] < 0.0) lonp[i] += 360.0;
		 }

	       memcpy(GRef->lon, lonp, ni*nj*sizeof(ftnfloat));
	       memcpy(GRef->lat, latp, ni*nj*sizeof(ftnfloat));
	       free(lonp);
	       free(latp);
	       break;

	     case 'L':
	       for (j=0; j < nj; j++)
		 {
		 for (i=0; i < ni; i++)
		   {
		   GRef->lat[C_TO_FTN(i,j,ni)] += 1.0;
		   GRef->lon[C_TO_FTN(i,j,ni)] += 1.0;
		   }
		 }
	       c_llfgr(GRef->lat, GRef->lon, GRef->lon, GRef->lat, ni*nj,
		       GRef->fst.xgref[SWLAT],GRef->fst.xgref[SWLON], GRef->fst.xgref[DLAT], GRef->fst.xgref[DLON]);
	       break;

	     case 'E':
	       latp = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
	       lonp = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
	       f77name(ez_gfllfxy)(lonp,latp,GRef->lon,GRef->lat,&npts,
				   &GRef->fst.xgref[XLAT1],&GRef->fst.xgref[XLON1], &GRef->fst.xgref[XLAT2], &GRef->fst.xgref[XLON2]);
	       memcpy(GRef->lon, lonp, ni*nj*sizeof(ftnfloat));
	       memcpy(GRef->lat, latp, ni*nj*sizeof(ftnfloat));
	       free(lonp);
	       free(latp);
	       break;

	     }
           break;

         case '!':
	   x = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
	   y = (ftnfloat *) malloc(ni*nj*sizeof(ftnfloat));
	   for (j=0; j < nj; j++)
	     {
	     for (i=0; i < ni; i++)
	       {
	       x[C_TO_FTN(i,j,ni)] = (ftnfloat) (i+1.0);
	       y[C_TO_FTN(i,j,ni)] = (ftnfloat) (j+1.0);
	       }
	     }
	   f77name(ez_llflamb)(GRef->lat,GRef->lon,x,y,&npts,
			       &GRef->grtyp, &GRef->fst.ig[IG1],&GRef->fst.ig[IG2],
			       &GRef->fst.ig[IG3],&GRef->fst.ig[IG4],1);
	   for (i=0; i < npts; i++)
	     {
	     if (GRef->lon[i] < 0.0)
	       {
	       GRef->lon[i] += 360.0;
	       }
	     }
	   break;
         }

      switch(GRef->grtyp[0])
	{
	case 'G':
	case 'B':
	case 'A':
	  if (GRef->fst.ig[IG2] == 1)
	    {
	    f77name(permut)(GRef->lat, &GRef->ni, &GRef->nj);
	    }
	  break;

	default:
	  break;
	}

      GRef->flags |= LAT;
      }

   if (groptions.verbose == 2)
     {
     fprintf(stderr, "gdid: %d\n", GRef->index);
     for (j=0; j < nj; j++)
       {
       for (i=0; i < ni; i++)
	 {
	 fprintf(stderr, "%d %d %f %f\n", i,j,GRef->lat[i], GRef->lon[i]);
	 }
       }
     }
   return 0;


}
