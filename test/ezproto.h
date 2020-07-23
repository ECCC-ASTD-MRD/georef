#ifndef _EZSCINT_H
#define _EZSCINT_H

#include "rpnmacros.h"
#include "rpn_macros_arch.h"
#include "../src/GeoRef.h"

// RPN external EZscint functions
extern int  c_ezfreegridset(TGeoRef* gr, int index);
extern int  c_ezdefset(TGeoRef* gdout, TGeoRef* gdin);
extern int  c_ezgprm(int gdid, char *grtyp, int *ni, int *nj, int *ig1, int *ig2, int *ig3, int *ig4);
extern int  c_ezgenpole(ftnfloat *vpolnor, ftnfloat *vpolsud, ftnfloat *fld,int ni, int nj, int vecteur,char *grtyp, int hem);
extern int  c_ezgetopt(char *option, char *value);
extern int  c_ezgetval(char *option, ftnfloat *value);
extern int  c_gdll(TGeoRef* GRef, ftnfloat *lat, ftnfloat *lon);

extern int  c_gdrls(TGeoRef* GRef);
extern int  c_ezsetopt(char *option, char *value);
extern int  c_ezsetval(char *option, ftnfloat fvalue);
extern int  c_ezsint(ftnfloat *zout, ftnfloat *zin);
extern int  c_ezuvint(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin);
extern int  c_ezwdint(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin);
extern int  c_gdgaxes(TGeoRef* GRef, ftnfloat *ax, ftnfloat *ay);
extern int  c_gdgxpndaxes(TGeoRef* GRef, ftnfloat *ax, ftnfloat *ay);
extern int  c_gdllfxy(TGeoRef* GRef, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, int n);
extern int  c_gdllfxyz(TGeoRef* GRef, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, int n);
extern int  c_gdllsval(TGeoRef* GRef, ftnfloat *zout, ftnfloat *zin, ftnfloat *lat, ftnfloat *lon, int n);
extern int  c_gdllvval(TGeoRef* GRef, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin,ftnfloat *lat, ftnfloat *lon, int n);
extern int  c_gdllwdval(TGeoRef* GRef, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin,ftnfloat *lat, ftnfloat *lon, int n);
extern int  c_gdxpncf(int gdin, int *i1, int *i2, int *j1, int *j2);
extern int  c_gdxysval(TGeoRef* gdin, ftnfloat *zout, ftnfloat *zin, ftnfloat *x, ftnfloat *y, int n);
extern int  c_gdxywdval(TGeoRef* gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, int n);
extern int  c_gdxyvval(TGeoRef* gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, int n);
extern int  c_gduvfwd(TGeoRef* GRef,  ftnfloat *uugdout, ftnfloat *vvgdout, ftnfloat *uullin, ftnfloat *vvllin,ftnfloat *latin, ftnfloat *lonin, int npts);
extern int  c_gdwdfuv(TGeoRef* GRef, ftnfloat *uullout, ftnfloat *vvllout, ftnfloat *uuin, ftnfloat *vvin,ftnfloat *latin, ftnfloat *lonin, int npts);
extern int  c_gdxyfll(TGeoRef* GRef, ftnfloat *x, ftnfloat *y, ftnfloat *lat, ftnfloat *lon, int n);
extern int  c_gdxyzfll(TGeoRef* GRef, ftnfloat *x, ftnfloat *y, ftnfloat *lat, ftnfloat *lon, int n);
extern void c_ezgfllfxy(ftnfloat *lonp, ftnfloat *latp,ftnfloat *lon, ftnfloat *lat,ftnfloat *r, ftnfloat *ri, int *npts,ftnfloat *xlat1, ftnfloat *xlon1, ftnfloat *xlat2, ftnfloat *xlon2);
extern void c_ezgfxyfll(ftnfloat *lonp, ftnfloat *latp,ftnfloat *lon, ftnfloat *lat,ftnfloat *r, ftnfloat *ri, int *npts,ftnfloat *xlat1, ftnfloat *xlon1, ftnfloat *xlat2, ftnfloat *xlon2);
extern void c_ezgfwfllw(ftnfloat *uullout, ftnfloat *vvllout, ftnfloat *latin, ftnfloat *lonin,ftnfloat *xlatingf, ftnfloat *xloningf,int *ni, int *nj,char *grtyp, int *ig1, int *ig2, int *ig3, int *ig4);
extern void c_ezllwfgfw(ftnfloat *uullout, ftnfloat *vvllout, ftnfloat *latin, ftnfloat *lonin,ftnfloat *xlatingf, ftnfloat *xloningf,int *ni,int *nj,char *grtyp,int *ig1,int *ig2,int *ig3,int *ig4);
extern void c_ezdefxg(TGeoRef* GRef);
extern void c_ezdefaxes(TGeoRef* GRef, ftnfloat *ax, ftnfloat *ay);
extern int  c_gdinterp(ftnfloat *zout, ftnfloat *zin, TGeoRef* gdin, ftnfloat *x, ftnfloat *y, int npts);
extern int  c_gdsetmask(TGeoRef* gr, int *mask);
extern int  c_gdgetmask(TGeoRef* gr, int *mask);
extern int  c_ezsint_m(float *zout, float *zin);
extern int  c_ezuvint_m(float *uuout, float *vvout, float *uuin, float *vvin);
extern int  c_ezsint_mdm(float *zout, int *mask_out, float *zin, int *mask_in);
extern int  c_ezuvint_mdm(float *uuout, float *vvout, int *mask_out, float *uuin, float *vvin, int *mask_in);
extern int  c_ezsint_mask(int *mask_out, int *mask_in);
extern int  c_ez_refgrid(TGeoRef* GRef);
extern int  c_fst_data_length(int size);

#endif