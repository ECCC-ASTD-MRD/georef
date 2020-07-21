/*
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
wordint c_ezyysint(ftnfloat *zout, ftnfloat *zin, TGeoRef *gdout, TGeoRef *gdin)
{
  wordint icode,i,j,k,ierc1,ierc2,ierc;
  wordint yancount_yin,yincount_yin, yancount_yan,yincount_yan;
  wordint yyin,yyout;
  wordint ni, nj;
  /*wordint yin2yin,yan2yin,yin2yan,yan2yan;*/
  int idx_gdin;
  ftnfloat *yin2yin_zvals,*yan2yin_zvals;
  ftnfloat *yin2yan_zvals,*yan2yan_zvals;
  
  TGeoRef *yin_gdin, *yan_gdin, *yin_gdout, *yan_gdout;
 /*  need only access to either yin or Yang info for the lat and lon val */
   
  yyin=0; yyout=0; 
  ierc=0;
  ierc1=0;ierc2=0;

  idx_gdin = c_find_gdin(gdin, gdout);

/* setup for input grid */
  if (gdin->NbSub > 0)
     {
     yyin=1;
     yin_gdin = gdin->Subs[0];
     yan_gdin = gdin->Subs[1];
     /*fprintf(stderr,"<c_ezyysint> finds Subs %d and %d in src Yin-Yang grid\n",yin_gdin,yan_gdin);*/
     }
  else
     {
     yin_gdin = gdin;
     }

/* setup for output grid */
  if (gdout->NbSub > 0)
     {
     yyout=1;
     yin_gdout = gdout->Subs[0];
     yan_gdout = gdout->Subs[1];
     }
  else
     {
     yin_gdout = gdout;
     }
  
  ni = yin_gdout->ni;
  nj = yin_gdout->nj;

/* interp input one grid to yygrid - no masking needed*/
  if (yyin == 0 && yyout == 1)
    {
    icode = c_ezdefset(yin_gdout,gdin);
    ierc1 = c_ezsint_orig(zout,zin);
    icode = c_ezdefset(yan_gdout,gdin);
    ierc2 = c_ezsint_orig(&zout[ni*nj],zin);   
    if (ierc1 == 2 || ierc2 == 2)
       {
        ierc=2;
       }
    return ierc;
    }

  /* check if one input sub grid is identical to dest sub grid or dest single grid */
  
  if (yin_gdin == gdout)
     {
     icode = c_ezdefset(gdout,yin_gdin);
     icode = c_ezsint_orig(zout,zin);
     return icode;
     }
  if (yan_gdin == gdout)
     {
     icode = c_ezdefset(gdout,yan_gdin);
     icode = c_ezsint_orig(zout,&zin[(yin_gdin->ni)*(yin_gdin->nj)]);   
     return icode;
     }

  /* User specifies to use 1 subgrid for interpolation ezsetopt(USE_1SUBGRID) */
  /* User must specify the sub grid value in ezsetival(SUBGRIDID) */
  /* This is only appropriate if the destination grid is non yin-yang grid */
  if (groptions.use_1subgrid == 1) /* User specifies to use 1 subgrid only */
     {
     if (yyout == 1) /* output is a Yin-Yang grid */
        {
         fprintf(stderr,"<c_ezyysint> cannot use 1 subgrid to interpolate to a Yin-Yang grid  Aborting...\n");
         return -1;
        }
     if (groptions.valeur_1subgrid != yin_gdin &&  
         groptions.valeur_1subgrid != yan_gdin)  /* chosen subgrid is neither Yin or Yang in source grid */
        {
         fprintf(stderr,"<c_ezyysint> define src subgridid in ezsetival(subgridid)! Aborting...\n");
         return -1;
        }
     if (groptions.valeur_1subgrid == yin_gdin) /*User specifies to use yin input grid */
        {
        icode = c_ezdefset(yin_gdout,groptions.valeur_1subgrid);
        ierc = c_ezsint_orig(zout,zin);
        return ierc;
        }
     if (groptions.valeur_1subgrid == yan_gdin) /* User specifies to use Yang input grid */
        {
        icode = c_ezdefset(yin_gdout,groptions.valeur_1subgrid);
        ierc = c_ezsint_orig(zout,&zin[(yin_gdin->ni)*(yin_gdin->nj)]);
        return ierc;
        }

     }
  /*End of ONE grid option*/   

  /* TO USE both Yin and Yang grids in Yin-yang input grid */
  /* Masquer les grilles YY input pour enlever overlap et calculer les X,Y */
  icode = c_ezyy_calcxy(gdout,gdin);

/* interp yinyang to one grid */
  if (yyin == 1 && yyout == 0)
    {
    yincount_yin = gdout->gset[idx_gdin].yincount_yin;
    yancount_yin = gdout->gset[idx_gdin].yancount_yin;
    yin2yin_zvals = (ftnfloat *) malloc(yincount_yin*sizeof(ftnfloat));
    yan2yin_zvals = (ftnfloat *) malloc(yancount_yin*sizeof(ftnfloat));

    icode = c_gdxysval(yin_gdin,yin2yin_zvals,zin,
           gdout->gset[idx_gdin].yin2yin_x,gdout->gset[idx_gdin].yin2yin_y,
           gdout->gset[idx_gdin].yincount_yin);

    icode = c_gdxysval(yan_gdin,yan2yin_zvals,
            &zin[(yin_gdin->ni)*(yin_gdin->nj)],
    gdout->gset[idx_gdin].yan2yin_x,gdout->gset[idx_gdin].yan2yin_y,
    gdout->gset[idx_gdin].yancount_yin);

    yincount_yin=0;
    yancount_yin=0;
    for(j=0; j<nj; j++)
      {
      for (i=0;i<ni; i++)
        {
        k=(j*ni)+i;
        if (gdout->gset[idx_gdin].yin_maskout[k] == 1.0)
          {
          zout[k]=yan2yin_zvals[yancount_yin]; 
          yancount_yin++;
          }
        else
          {
          zout[k]=yin2yin_zvals[yincount_yin]; 
          yincount_yin++;
          }
        }
      }
    free(yin2yin_zvals);
    free(yan2yin_zvals);
    return icode;
   }

/* interp yinyang to yinyang*/
  if (yyout == 1 && yyin == 1)
    {
/* interp input YY grid to YY grid */
    yincount_yin = gdout->gset[idx_gdin].yincount_yin;
    yancount_yin = gdout->gset[idx_gdin].yancount_yin;
    yincount_yan = gdout->gset[idx_gdin].yincount_yan;
    yancount_yan = gdout->gset[idx_gdin].yancount_yan;
    yin2yin_zvals = (ftnfloat *) malloc(yincount_yin*sizeof(ftnfloat));
    yan2yin_zvals = (ftnfloat *) malloc(yancount_yin*sizeof(ftnfloat));
    yin2yan_zvals = (ftnfloat *) malloc(yincount_yan*sizeof(ftnfloat));
    yan2yan_zvals = (ftnfloat *) malloc(yancount_yan*sizeof(ftnfloat));
    
    icode = c_gdxysval(yin_gdin,yin2yin_zvals,zin,
            gdout->gset[idx_gdin].yin2yin_x,gdout->gset[idx_gdin].yin2yin_y,
            gdout->gset[idx_gdin].yincount_yin);
    icode = c_gdxysval(yan_gdin,yan2yin_zvals,&zin[(yin_gdin->ni)*(yin_gdin->nj)],
            gdout->gset[idx_gdin].yan2yin_x,gdout->gset[idx_gdin].yan2yin_y,
            gdout->gset[idx_gdin].yancount_yin);

    icode = c_gdxysval(yin_gdin,yin2yan_zvals,zin,
            gdout->gset[idx_gdin].yin2yan_x,gdout->gset[idx_gdin].yin2yan_y,
            gdout->gset[idx_gdin].yincount_yan);
    icode = c_gdxysval(yan_gdin,yan2yan_zvals,&zin[(yin_gdin->ni)*(yin_gdin->nj)],
            gdout->gset[idx_gdin].yan2yan_x,gdout->gset[idx_gdin].yan2yan_y,
       gdout->gset[idx_gdin].yancount_yan);

/* interp input YY grid to Yin grid */
    yincount_yin=0; yancount_yin=0;
    for(j=0; j<nj; j++)
      {
      for (i=0;i<ni; i++)
        {
        k=(j*ni)+i;
        if (gdout->gset[idx_gdin].yin_maskout[k] == 1.0)
          {
          zout[k]=yan2yin_zvals[yancount_yin]; 
          yancount_yin++;
          }
        else
          {
          zout[k]=yin2yin_zvals[yincount_yin]; 
          yincount_yin++;
          }
        }
      }
/* interp input YY grid to Yang grid */
    yincount_yan=0; yancount_yan=0;
    for(j=0; j<nj; j++)
      {
      for (i=0;i<ni; i++)
        {
        k=(j*ni)+i;
        if (gdout->gset[idx_gdin].yan_maskout[k] == 1.0)
          {
          zout[k+(ni*nj)]=yan2yan_zvals[yancount_yan]; 
          yancount_yan++;
          }
        else
          {
          zout[k+(ni*nj)]=yin2yan_zvals[yincount_yan]; 
          yincount_yan++;
          }
        }
      }
   free(yin2yin_zvals);
   free(yan2yin_zvals);
   free(yin2yan_zvals);
   free(yan2yan_zvals);
   }

   return icode;
}

