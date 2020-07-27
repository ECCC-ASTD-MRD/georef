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
PTR_AS_INT f77name(ezgdef_supergrid)(wordint *ni, wordint *nj, char *grtyp, char *grref, wordint *vercode, wordint *NbSub, PTR_AS_INT Subs, F2Cl lengrtyp, F2Cl lengrref)
{
  char lgrtyp[2],lgrref[2];

  lgrtyp[0] = grtyp[0];
  lgrtyp[1] = '\0';
   
  lgrref[0] = grref[0];
  lgrref[1] = '\0';
  return (PTR_AS_INT) c_ezgdef_supergrid(*ni, *nj, lgrtyp, lgrref, *vercode, *NbSub,(TGeoRef**)Subs);
}

TGeoRef* c_ezgdef_supergrid(wordint ni, wordint nj, char *grtyp, char *grref, wordint vercode,wordint NbSub, TGeoRef **Subs)
{
  wordint  i;
  ftnfloat *ax,*ay;
  TGeoRef *GRef, *fref, *sub_gd;
    
  if (NbSub <= 1)
  {
    fprintf(stderr,"<c_ezgdef_supergrid> NbSub given is less than 2! Aborting...\n");
    return NULL;
  }
  if (vercode != 1)
  {
    fprintf(stderr,"<c_ezgdef_supergrid> invalid vercode! Aborting...\n");
    return NULL;
  }

  GRef = GeoRef_New();
  strcpy(GRef->fst.etiketx, "            ");
  strcpy(GRef->fst.etikety, "            ");
  strcpy(GRef->fst.typvarx, "  ");
  strcpy(GRef->fst.typvary, "  ");
  strcpy(GRef->fst.nomvarx, "^>  ");
  strcpy(GRef->fst.nomvary, "^>  ");
  
  if (vercode == 1)
  {
    sub_gd = Subs[0];

/*    strcpy(newgr.grtyp, grtyp);
    strcpy(newgr.grref, grref);
*/
    GRef->grtyp[0]=grtyp[0];
    GRef->grref[0]=grref[0];
    RemplirDeBlancs(GRef->fst.nomvarx, 5);
    RemplirDeBlancs(GRef->fst.typvarx, 3);
    RemplirDeBlancs(GRef->fst.etiketx, 13);
    RemplirDeBlancs(GRef->fst.nomvary, 5);
    RemplirDeBlancs(GRef->fst.typvary, 3);
    RemplirDeBlancs(GRef->fst.etikety, 13);
    GRef->ni = ni;
    GRef->nj = nj;
    GRef->idx_last_gdin = -1;
    /* create tictac arrays to add uniqueness in supergrid*/
    ax = (ftnfloat *) malloc(GRef->ni*sizeof(ftnfloat));
    ay = (ftnfloat *) malloc(GRef->nj*sizeof(ftnfloat));
    memcpy(ax,sub_gd->ax,GRef->ni*sizeof(ftnfloat));
    memcpy(ay,sub_gd->ay,GRef->nj*sizeof(ftnfloat));
    GRef->fst.ip1      = sub_gd->fst.ip1;
    GRef->fst.ip2      = sub_gd->fst.ip2;
    GRef->fst.ip3      = sub_gd->fst.ip3;
/*   to add more uniqueness to the super-grid index Yin-Yang grid, we also
add   the rotation of YIN */
    GRef->fst.ig[IG1]  = sub_gd->fst.ig[IG1];
    GRef->fst.ig[IG2]  = sub_gd->fst.ig[IG2];
    GRef->fst.ig[IG3]  = sub_gd->fst.ig[IG3];
    GRef->fst.ig[IG4]  = sub_gd->fst.ig[IG4];
    GRef->fst.xg[IG1]  = 0.0;
    GRef->fst.xg[IG2]  = 0.0;
    GRef->fst.xg[IG3]  = 0.0;
    GRef->fst.xg[IG4]  = 0.0;
    GRef->fst.igref[IG1]=vercode;
    GRef->fst.igref[IG2]=0;
    GRef->fst.igref[IG3]=0;
    GRef->fst.igref[IG4]=0;
    GRef->NbSub= NbSub;
  }

  free(ax); free(ay);
  
  // This georef already exists
  if (fref=GeoRef_Find(GRef)) {
    free(GRef);
    GeoRef_Incr(fref);
    return fref;
  }

  // This is a new georef
  GeoRef_Add(GRef);

  GRef->Subs = (TGeoRef **) malloc(NbSub*sizeof(TGeoRef*));

  for (i=0; i < NbSub; i++)
  {
    GRef->Subs[i] = Subs[i];
    c_ezgdef_yymask(Subs[i]);
    if (groptions.verbose > 0)
    {
      printf("Grille[%p].Subs[%p] has maskgrid=%p\n",GRef,Subs[i],sub_gd->mymaskgrid);
    }
  }

    if (groptions.verbose > 0)
    {
    printf("Grille[%p].nomvarx=%s\n",GRef,GRef->fst.nomvarx);
    printf("Grille[%p].nomvary=%s\n",GRef,GRef->fst.nomvary);
    printf("Grille[%p].etikx=%s\n",GRef,GRef->fst.etiketx);
    printf("Grille[%p].etiky=%s\n",GRef,GRef->fst.etikety);
    printf("Grille[%p].grtyp = '%c'\n", GRef, GRef->grtyp[0]);
    printf("Grille[%p].grref = '%c'\n", GRef, GRef->grref[0]);
    printf("Grille[%p].ni    = %d\n",   GRef, GRef->ni);
    printf("Grille[%p].nj    = %d\n",   GRef, GRef->nj);
    printf("Grille[%p].ip1   = %d\n",   GRef, GRef->fst.ip1);
    printf("Grille[%p].ip2   = %d\n",   GRef, GRef->fst.ip2);
    printf("Grille[%p].ip3   = %d\n",   GRef, GRef->fst.ip3);
    printf("Grille[%p].ig1   = %d\n",   GRef, GRef->fst.ig[IG1]);
    printf("Grille[%p].ig2   = %d\n",   GRef, GRef->fst.ig[IG2]);
    printf("Grille[%p].ig3   = %d\n",   GRef, GRef->fst.ig[IG3]);
    printf("Grille[%p].ig4   = %d\n",   GRef, GRef->fst.ig[IG4]);
    printf("Grille[%p].ig1ref = %d\n",   GRef, GRef->fst.igref[IG1]);
    printf("Grille[%p].ig2ref = %d\n",   GRef, GRef->fst.igref[IG2]);
    printf("Grille[%p].ig3ref = %d\n",   GRef, GRef->fst.igref[IG3]);
    printf("Grille[%p].ig4ref = %d\n",   GRef, GRef->fst.igref[IG4]);
    printf("Grille[%p].NbSub = %d\n", GRef, GRef->NbSub);
    printf("Grille[%p].Subs[0]   = %p\n",   GRef, GRef->Subs[0]);
    printf("Grille[%p].Subs[1]   = %p\n",   GRef, GRef->Subs[1]);

    printf("Grille[%p].fst.xg[1]   = %f\n",   GRef, GRef->fst.xg[1]);
    printf("Grille[%p].fst.xg[2]   = %f\n",   GRef, GRef->fst.xg[2]);
    printf("Grille[%p].fst.xg[3]   = %f\n",   GRef, GRef->fst.xg[3]);
    printf("Grille[%p].fst.xg[4]   = %f\n",   GRef, GRef->fst.xg[4]);
    }
    return GRef;
}

