#ifndef _EZSCINT_H
#define _EZSCINT_H

#include "rpnmacros.h"
#include "rpn_macros_arch.h"
#include "../src/GeoRef.h"

// RPN external EZscint functions
extern int  c_ezfreegridset(int gdid, int index);
extern int  c_ezdefset(int gdout, int gdin);
extern int  c_ezgdef(int ni, int nj, char *grtyp, char *grref,int ig1, int ig2, int ig3, int ig4, ftnfloat *ax, ftnfloat *ay);
extern int  c_ezgdef_ffile(int ni, int nj, char *grtyp,int ig1, int ig2, int ig3, int ig4, int iunit);
extern int  c_ezgdef_fll(int ni, int nj,ftnfloat *lat, ftnfloat *lon);
extern TGeoRef*  c_ezgdef_fmem(int ni, int nj, char *grtyp, char *grref,int ig1, int ig2, int ig3, int ig4, ftnfloat *ax, ftnfloat *ay);
extern int  c_ezgprm(int gdid, char *grtyp, int *ni, int *nj, int *ig1, int *ig2, int *ig3, int *ig4);
extern int  c_ezgenpole(ftnfloat *vpolnor, ftnfloat *vpolsud, ftnfloat *fld,int ni, int nj, int vecteur,char *grtyp, int hem);
extern int  c_ezgetopt(char *option, char *value);
extern int  c_ezgetval(char *option, ftnfloat *value);
extern int  c_ezget_nsubgrids(int id);
extern int  c_ezget_subgridids(int id,int *subid);
extern int  c_gdll(int gdid, ftnfloat *lat, ftnfloat *lon);
extern int  c_ezqkdef(int ni, int nj, char *grtyp,int ig1, int ig2, int ig3, int ig4, int iunit);
             
extern int  c_ezquickdef(int ni, int nj, char *grtyp,int ig1, int ig2, int ig3, int ig4, int iunit);
extern int  c_gdrls(TGeoRef* GRef);
extern int  c_ezsetopt(char *option, char *value);
extern int  c_ezsetval(char *option, ftnfloat fvalue);
extern int  c_ezsint(ftnfloat *zout, ftnfloat *zin);
extern int  c_ezuvint(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin);
extern int  c_ezwdint(ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin);
extern int  c_gdgaxes(int gdid, ftnfloat *ax, ftnfloat *ay);
extern int  c_gdgxpndaxes(int gdid, ftnfloat *ax, ftnfloat *ay);
extern int  c_gdllfxy(int gdid, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, int n);
extern int  c_gdllfxyz(int gdid, ftnfloat *lat, ftnfloat *lon, ftnfloat *x, ftnfloat *y, int n);
extern int  c_gdllsval(int gdid, ftnfloat *zout, ftnfloat *zin, ftnfloat *lat, ftnfloat *lon, int n);
extern int  c_gdllvval(int gdid, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin,ftnfloat *lat, ftnfloat *lon, int n);
extern int  c_gdllwdval(int gdid, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin,ftnfloat *lat, ftnfloat *lon, int n);
extern int  c_gdxpncf(int gdin, int *i1, int *i2, int *j1, int *j2);
extern int  c_gdxysval(int gdin, ftnfloat *zout, ftnfloat *zin, ftnfloat *x, ftnfloat *y, int n);
extern int  c_gdxywdval(int gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, int n);
extern int  c_gdxyvval(int gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin, ftnfloat *vvin, ftnfloat *x, ftnfloat *y, int n);
extern int  c_gduvfwd(int gdid,  ftnfloat *uugdout, ftnfloat *vvgdout, ftnfloat *uullin, ftnfloat *vvllin,ftnfloat *latin, ftnfloat *lonin, int npts);
extern int  c_gdwdfuv(int gdid, ftnfloat *uullout, ftnfloat *vvllout, ftnfloat *uuin, ftnfloat *vvin,ftnfloat *latin, ftnfloat *lonin, int npts);
extern int  c_gdxpngd(int gdin, ftnfloat *zxpnded, ftnfloat *zin);
extern int  c_gdxyfll(int gdid, ftnfloat *x, ftnfloat *y, ftnfloat *lat, ftnfloat *lon, int n);
extern int  c_gdxyzfll(int gdid, ftnfloat *x, ftnfloat *y, ftnfloat *lat, ftnfloat *lon, int n);
extern int  c_guval(int gdin, ftnfloat *uuout, ftnfloat *vvout, ftnfloat *uuin,  ftnfloat *vvin, ftnfloat *x, ftnfloat *y, int n);
extern void c_ezgfllfxy(ftnfloat *lonp, ftnfloat *latp,ftnfloat *lon, ftnfloat *lat,ftnfloat *r, ftnfloat *ri, int *npts,ftnfloat *xlat1, ftnfloat *xlon1, ftnfloat *xlat2, ftnfloat *xlon2);
extern void c_ezgfxyfll(ftnfloat *lonp, ftnfloat *latp,ftnfloat *lon, ftnfloat *lat,ftnfloat *r, ftnfloat *ri, int *npts,ftnfloat *xlat1, ftnfloat *xlon1, ftnfloat *xlat2, ftnfloat *xlon2);
extern void c_ezgfwfllw(ftnfloat *uullout, ftnfloat *vvllout, ftnfloat *latin, ftnfloat *lonin,ftnfloat *xlatingf, ftnfloat *xloningf,int *ni, int *nj,char *grtyp, int *ig1, int *ig2, int *ig3, int *ig4);
extern void c_ezllwfgfw(ftnfloat *uullout, ftnfloat *vvllout, ftnfloat *latin, ftnfloat *lonin,ftnfloat *xlatingf, ftnfloat *xloningf,int *ni,int *nj,char *grtyp,int *ig1,int *ig2,int *ig3,int *ig4);
extern void c_ezdefxg(TGeoRef* GRef);
extern void c_ezdefaxes(TGeoRef* GRef, ftnfloat *ax, ftnfloat *ay);
extern int  c_gdinterp(ftnfloat *zout, ftnfloat *zin, int gdin, ftnfloat *x, ftnfloat *y, int npts);
extern int  c_gdsetmask(int gdid, int *mask);
extern int  c_gdgetmask(int gdid, int *mask);
extern int  c_ezsint_m(float *zout, float *zin);
extern int  c_ezuvint_m(float *uuout, float *vvout, float *uuin, float *vvin);
extern int  c_ezsint_mdm(float *zout, int *mask_out, float *zin, int *mask_in);
extern int  c_ezuvint_mdm(float *uuout, float *vvout, int *mask_out, float *uuin, float *vvin, int *mask_in);
extern int  c_ezsint_mask(int *mask_out, int *mask_in);
extern int  c_ez_refgrid(int gdid);
extern int  c_fst_data_length(int size);

#endif