#include "../src/GeoRef.h"

#ifndef _ezfuncdef
#include "gd_key2rowcol.c"

void EliminerGrille(TGeoRef* GRef);

void f77name(ez_avg)(float* zout, float* x, float* y, int* ni_src, int* nj_src,
            float* zin, int* ni_dst, int* nj_dst, int* extension);

wordint LireEnrPositionnels(TGeoRef* gr, wordint iunit, wordint ip1, wordint ip2, wordint ip3, wordint ip4, wordint read);

void c_llfgr(ftnfloat* lat, ftnfloat* lon, ftnfloat* x, ftnfloat* y, wordint npts,
        ftnfloat latOrigine, ftnfloat lonOrigine, ftnfloat deltaLat, ftnfloat deltaLon);

wordint ez_calclatlon(TGeoRef* GRef);

void ez_calcntncof(TGeoRef* GRef);

wordint ez_calcxpncof(TGeoRef* GRef);

wordint ez_calcxy(TGeoRef* gdin, TGeoRef* gdout);

wordint ez_corrval(ftnfloat* zout, ftnfloat* zin, TGeoRef* gdin, TGeoRef* gdout);
wordint ez_corrvec(ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin, TGeoRef* gdin, TGeoRef* gdout);
wordint ez_corrval_ausud(ftnfloat* zout, ftnfloat* zin,  TGeoRef* gdin, TGeoRef* gdout);
wordint ez_corrval_aunord(ftnfloat* zout, ftnfloat* zin,  TGeoRef* gdin, TGeoRef* gdout);
wordint ez_corrvec_aunord(ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin, TGeoRef* gdin, TGeoRef* gdout);
wordint ez_corrvec_ausud(ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin, TGeoRef* gdin, TGeoRef* gdout);

wordint ez_defzones(TGeoRef* gdin, TGeoRef* gdout);

wordint ez_defzone_dehors(TGeoRef* gdin, ftnfloat* x, ftnfloat* y, wordint npts, _zone* zone);
wordint ez_defzone_polenord(TGeoRef* gdin, ftnfloat* x, ftnfloat* y, wordint npts, _zone* zone);
wordint ez_defzone_polesud(TGeoRef* gdin, ftnfloat* x, ftnfloat* y, wordint npts, _zone* zone);
wordint ez_defzone_nord(TGeoRef* gdin, ftnfloat* x, ftnfloat* y, wordint npts, _zone* zone);
wordint ez_defzone_sud(TGeoRef* gdin, ftnfloat* x, ftnfloat* y, wordint npts, _zone* zone);

wordint ez_interp(ftnfloat* zout, ftnfloat* zin, TGeoRef* gdin, TGeoRef* gdout);

void ez_xpncof(wordint* i1, wordint* i2, wordint* j1, wordint* j2, wordint* couverture,
            wordint ni, wordint nj, char grtyp, char grref,
            wordint ig1, wordint ig2, wordint ig3, wordint ig4, wordint sym,
            ftnfloat* ax, ftnfloat* ay);

void ez_xpnsrcgd(TGeoRef* GRef, ftnfloat* zout, ftnfloat* zin);

wordint c_ezfreegridset(TGeoRef* gr, wordint index);

wordint f77name(ezdefset)(PTR_AS_INT gdout, PTR_AS_INT gdin);
wordint c_ezdefset(TGeoRef* gdout, TGeoRef* gdin);

PTR_AS_INT f77name(ezgdef_supergrid)(wordint* ni, wordint* nj, char* grtyp, char* grref, wordint* vercode, wordint* nsubgrids, PTR_AS_INT subgrid, F2Cl lengrtyp, F2Cl lengrref);
TGeoRef* c_ezgdef_supergrid(wordint ni, wordint nj, char* grtyp, char* grref, wordint vercode, wordint nsubgrids, TGeoRef** subgrid);

wordint c_ezgdef_yymask(TGeoRef* gr);

wordint f77name(ezgenpole)(ftnfloat* vpolnor, ftnfloat* vpolsud, ftnfloat* fld,
                           wordint* ni, wordint* nj, wordint* vecteur,
                           char* grtyp, wordint* hem, F2Cl lengrtyp);
wordint c_ezgenpole(ftnfloat* vpolnor, ftnfloat* vpolsud, ftnfloat* fld,
                           wordint ni, wordint nj, wordint vecteur,
                           char* grtyp, wordint hem);

wordint f77name(ezgetopt)(char* option, char* value, F2Cl lenoption, F2Cl lenvalue);
wordint c_ezgetopt(char* option, char* value);

wordint f77name(ezgetival)(char* option, ftnword* value, F2Cl lenoption);
wordint c_ezgetival(char* option, ftnword* value);

wordint f77name(ezgetval)(char* option, ftnfloat* value, F2Cl lenoption);
wordint c_ezgetval(char* option, ftnfloat* value);

wordint f77name(ezgfstp)(wordint* gdid,
         char* nomvarx, char* typvarx, char* etiketx,
         char* nomvary, char* typvary, char* etikety,
         wordint* ip1, wordint* ip2, wordint* ip3, wordint* dateo,
                     wordint* deet, wordint* npas, wordint* nbits,
         F2Cl lennomvarx, F2Cl lentypvarx, F2Cl lenetiketx,
         F2Cl lennomvary, F2Cl lentypvary, F2Cl lenetikety);
wordint c_ezgfstp(wordint gdid, char* nomvarx, char* typvarx, char* etiketx,
              char* nomvary, char* typvary, char* etikety,
              wordint* ip1, wordint* ip2, wordint* ip3, wordint* dateo, wordint* deet, wordint* npas, wordint* nbits);

wordint f77name(ezgprm)(wordint* gdid, char* grtyp, wordint* ni, wordint* nj,
             wordint* ig1, wordint* ig2, wordint* ig3, wordint* ig4, F2Cl lengrtyp);
wordint   c_ezgprm(wordint gdid, char* grtyp, wordint* ni, wordint* nj, wordint* ig1, wordint* ig2, wordint* ig3, wordint* ig4);

wordint f77name(ezgxprm)(wordint* gdid, wordint* ni, wordint* nj, char* grtyp,
                     wordint* ig1, wordint* ig2, wordint* ig3, wordint* ig4,
                     char* grref, wordint* ig1ref, wordint* ig2ref,
                     wordint* ig3ref, wordint* ig4ref,
                     F2Cl lengrtyp, F2Cl lengrref);
wordint c_ezgxprm(wordint gdid, wordint* ni, wordint* nj,
              char* grtyp, wordint* ig1, wordint* ig2, wordint* ig3, wordint* ig4,
              char* grref, wordint* ig1ref, wordint* ig2ref, wordint* ig3ref, wordint* ig4ref);

wordint f77name(gdll)(PTR_AS_INT GRef, ftnfloat* lat, ftnfloat* lon);
wordint c_gdll(TGeoRef* GRef, ftnfloat* lat, ftnfloat* lon);

wordint f77name(gdrls)(PTR_AS_INT GRef);
wordint c_gdrls(TGeoRef* GRef);

wordint f77name(ezsetopt)(char* option, char* value, F2Cl lenoption, F2Cl lenvalue);
wordint c_ezsetopt(char* option, char* value);

wordint f77name(ezsetival)(char* option, wordint* ivalue, F2Cl lenoption);
wordint c_ezsetival(char* option, wordint ivalue);

wordint f77name(ezsetval)(char* option, ftnfloat* fvalue, F2Cl lenoption);
wordint c_ezsetval(char* option, ftnfloat fvalue);

wordint f77name(ezsint)(ftnfloat* zout, ftnfloat* zin);
wordint c_ezsint(ftnfloat* zout, ftnfloat* zin);

wordint c_find_gdin(TGeoRef* gdin, TGeoRef* gdout);
wordint find_gdin_in_gset(TGeoRef* gdin, TGeoRef* gdout);

wordint f77name(ezuvint)(ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin);
wordint c_ezuvint(ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin);

wordint f77name(ezwdint)(ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin);
wordint c_ezwdint(ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin);

wordint ftnstrclean(char* str, wordint lenstr);

wordint f77name(gdgaxes)(PTR_AS_INT GRef, ftnfloat* ax, ftnfloat* ay);
wordint c_gdgaxes(TGeoRef* GRef, ftnfloat* ax, ftnfloat* ay);

wordint f77name(gdgxpndaxes)(PTR_AS_INT GRef, ftnfloat* ax, ftnfloat* ay);
wordint c_gdgxpndaxes(TGeoRef* GRef, ftnfloat* ax, ftnfloat* ay);

wordint f77name(gdllfxy)(PTR_AS_INT GRef, ftnfloat* lat, ftnfloat* lon, ftnfloat* x, ftnfloat* y, wordint* n);
wordint c_gdllfxy(TGeoRef* GRef, ftnfloat* lat, ftnfloat* lon, ftnfloat* x, ftnfloat* y, wordint n);

wordint f77name(gdllfxyz)(PTR_AS_INT GRef, ftnfloat* lat, ftnfloat* lon, ftnfloat* x, ftnfloat* y, wordint* n);
wordint c_gdllfxyz(TGeoRef* GRef, ftnfloat* lat, ftnfloat* lon, ftnfloat* x, ftnfloat* y, wordint n);

wordint f77name(gdllsval)(PTR_AS_INT GRef, ftnfloat* zout, ftnfloat* zin, ftnfloat* lat, ftnfloat* lon, wordint* n);
wordint c_gdllsval(TGeoRef* GRef, ftnfloat* zout, ftnfloat* zin, ftnfloat* lat, ftnfloat* lon, wordint n);

wordint f77name(gdllvval)(PTR_AS_INT GRef, ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin,
                      ftnfloat* lat, ftnfloat* lon, wordint* n);
wordint c_gdllvval(TGeoRef* GRef, ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin,
               ftnfloat* lat, ftnfloat* lon, wordint n);

wordint f77name(gdllwdval)(PTR_AS_INT GRef, ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin,
                      ftnfloat* lat, ftnfloat* lon, wordint* n);
wordint c_gdllwdval(TGeoRef* GRef, ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin,
               ftnfloat* lat, ftnfloat* lon, wordint n);

wordint f77name(gdxpncf)(wordint* gdin, wordint* i1, wordint* i2, wordint* j1, wordint* j2);
wordint c_gdxpncf(wordint gdin, wordint* i1, wordint* i2, wordint* j1, wordint* j2);

wordint f77name(gdxysval)(PTR_AS_INT gdin, ftnfloat* zout, ftnfloat* zin, ftnfloat* x, ftnfloat* y, wordint* n);
wordint c_gdxysval(TGeoRef* gdin, ftnfloat* zout, ftnfloat* zin, ftnfloat* x, ftnfloat* y, wordint n);

wordint f77name(gdxywdval)(PTR_AS_INT gdin, ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin, ftnfloat* x, ftnfloat* y, wordint* n);
wordint c_gdxywdval(TGeoRef* gdin, ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin, ftnfloat* x, ftnfloat* y, wordint n);

wordint f77name(gdxyvval)(PTR_AS_INT gdin, ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin, ftnfloat* x, ftnfloat* y, wordint* n);
wordint c_gdxyvval(TGeoRef* gdin, ftnfloat* uuout, ftnfloat* vvout, ftnfloat* uuin, ftnfloat* vvin, ftnfloat* x, ftnfloat* y, wordint n);

wordint f77name(gduvfwd)(PTR_AS_INT GRef, ftnfloat* uugdout, ftnfloat* vvgdout,
                     ftnfloat* uullin, ftnfloat* vvllin, ftnfloat* latin, ftnfloat* lonin, wordint* npts);
wordint c_gduvfwd(TGeoRef* GRef,  ftnfloat* uugdout, ftnfloat* vvgdout, ftnfloat* uullin, ftnfloat* vvllin,
              ftnfloat* latin, ftnfloat* lonin, wordint npts);

wordint f77name(gdwdfuv)(PTR_AS_INT GRef, ftnfloat* uullout, ftnfloat* vvllout, ftnfloat* uuin, ftnfloat* vvin,
              ftnfloat* latin, ftnfloat* lonin, wordint* npts);
wordint c_gdwdfuv(TGeoRef* GRef, ftnfloat* uullout, ftnfloat* vvllout, ftnfloat* uuin, ftnfloat* vvin,
              ftnfloat* latin, ftnfloat* lonin, wordint npts);

wordint f77name(gdxyfll)(PTR_AS_INT GRef, ftnfloat* x, ftnfloat* y, ftnfloat* lat, ftnfloat* lon, wordint* n);
wordint c_gdxyfll(TGeoRef* GRef, ftnfloat* x, ftnfloat* y, ftnfloat* lat, ftnfloat* lon, wordint n);

wordint f77name(gdxyzfll)(PTR_AS_INT GRef, ftnfloat* x, ftnfloat* y, ftnfloat* lat, ftnfloat* lon, wordint* n);
wordint c_gdxyzfll(TGeoRef* GRef, ftnfloat* x, ftnfloat* y, ftnfloat* lat, ftnfloat* lon, wordint n);

TGeoRef* c_ezgetgdin();

TGeoRef* c_ezgetgdout();

void c_ezgfllfxy(ftnfloat* lonp, ftnfloat* latp,
                 ftnfloat* lon, ftnfloat* lat,
                 ftnfloat* r, ftnfloat* ri, wordint* npts,
                 ftnfloat* xlat1, ftnfloat* xlon1, ftnfloat* xlat2, ftnfloat* xlon2);

void c_ezgfxyfll(ftnfloat* lonp, ftnfloat* latp,
                 ftnfloat* lon, ftnfloat* lat,
                 ftnfloat* r, ftnfloat* ri, wordint* npts,
                 ftnfloat* xlat1, ftnfloat* xlon1, ftnfloat* xlat2, ftnfloat* xlon2);

void c_ezgfwfllw(ftnfloat* uullout, ftnfloat* vvllout, ftnfloat* latin, ftnfloat* lonin,
                  ftnfloat* xlatingf, ftnfloat* xloningf,
                  wordint* ni, wordint* nj,
                  char* grtyp, wordint* ig1, wordint* ig2, wordint* ig3, wordint* ig4);

void  c_ezllwfgfw(ftnfloat* uullout, ftnfloat* vvllout, ftnfloat* latin, ftnfloat* lonin,
                  ftnfloat* xlatingf, ftnfloat* xloningf,
                 wordint* ni,wordint* nj,
                  char* grtyp,wordint* ig1,wordint* ig2,wordint* ig3,wordint* ig4);

void c_ez_manageGrillesMemory();
int c_ez_refgrid(TGeoRef* GRef);

void c_ezdefxg(TGeoRef* GRef);
void c_ezdefaxes(TGeoRef* GRef, ftnfloat* ax, ftnfloat* ay);
wordint c_gdinterp(ftnfloat* zout, ftnfloat* zin, TGeoRef* gdin, ftnfloat* x, ftnfloat* y, wordint npts);

int f77name(gdsetmask)(PTR_AS_INT gr, int* mask);
int f77name(gdgetmask)(PTR_AS_INT gr, int* mask);
int f77name(ezsint_m)(float* zout, float* zin);
int f77name(ezuvint_m)(float* uuout, float* vvout, float* uuin, float* vvin);
int f77name(ezsint_mdm)(float* zout, int* mask_out, float* zin, int* mask_in);
int f77name(ezuvint_mdm)(float* uuout, float* vvout, int* mask_out, float* uuin, float* vvin, int* mask_in);
int f77name(ezsint_mask)(int* mask_out, int* mask_in);

int c_gdsetmask(TGeoRef* gr, int* mask);
int c_gdgetmask(TGeoRef* gr, int* mask);
int c_ezsint_m(float* zout, float* zin);
int c_ezuvint_m(float* uuout, float* vvout, float* uuin, float* vvin);
int c_ezsint_mdm(float* zout, int* mask_out, float* zin, int* mask_in);
int c_ezuvint_mdm(float* uuout, float* vvout, int* mask_out, float* uuin, float* vvin, int* mask_in);
int c_ezsint_mask(int* mask_out, int* mask_in);

#endif
#define _ezfuncdef
