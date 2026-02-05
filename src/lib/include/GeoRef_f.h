#ifndef _GeoRef_f_h
#define _GeoRef_f_h

// Fortran wrapped declaration
void f77name(ez8_rgd_index_0)(float * const index, const double * const px, const double * const py, const int32_t * const npts, const int32_t * const ni, const int32_t * const j1, const int32_t * const j2);
void f77name(ez8_rgd_index_1)(float *index, double *px, double *py, int32_t *npts, int32_t *ni, int32_t *j1, int32_t *j2, int32_t *wrap);
void f77name(ez8_rgd_index_3)(float *index, double *px, double *py, int32_t *npts, int32_t *ni, int32_t *j1, int32_t *j2, int32_t *wrap);
void f77name(ez8_irgd_index_1)(float *index, double *px, double *py, int32_t *npts, int32_t *ni, int32_t *j1, int32_t *j2, double *ax, double *ay, int32_t *wrap);
void f77name(ez8_apply_0)(const float * const index, float * const zo, const int32_t * const npts, const float * const z, const int32_t * const ni, const int32_t * const j1, const int32_t * const j2, const float * const nodata);
void f77name(ez8_apply_1)(float *index, float *zo, int32_t *npts, float *z, int32_t *ni, int32_t *j1, int32_t *j2, float *nodata);
void f77name(ez8_apply_3)(float *index, float *zo, int32_t *npts, float *z, int32_t *ni, int32_t *j1, int32_t *j2, float *nodata);

void f77name(ez8_rgdint_0)      (float * const zo, const double * const px, const double * const py, const int32_t * const npts, const float * const z, const int32_t * const ni, const int32_t * const j1, const int32_t * const j2, const float * const nodata);
void f77name(ez8_rgdint_1)      (float *zo, double *px, double *py, int32_t *npts, float *z, int32_t *ni, int32_t *j1, int32_t *j2, int32_t *wrap, float *nodata);
void f77name(ez8_rgdint_3)      (float *zo, double *px, double *py, int32_t *npts, float *z, int32_t *ni, int32_t *j1, int32_t *j2, int32_t *wrap, float *nodata);
void f77name(ez8_irgdint_1)     (float *zo, double *px, double *py, int32_t *npts, float *z, int32_t *ni, int32_t *j1, int32_t *j2, double *ax, double *ay, int32_t *wrap, float *nodata);
void f77name(ez8_irgdint_3)     (float *zo, double *px, double *py, int32_t *npts, float *z, int32_t *ni, int32_t *i1, int32_t *i2, int32_t *j1, int32_t *j2, double *ax, double *ay, float *cx, float *cy, int32_t *wrap, float *nodata);
void f77name(ez8_irgdint_3_wnnc)(float *zo, double *px, double *py, int32_t *npts, float *z, int32_t *ni, int32_t *j1, int32_t *j2, double *ax, double *ay, int32_t *wrap);

void f77name(ez8_llwfgdw)(float *z1, float *z2, double *xlon, int32_t *li, int32_t *lj, char *grtyp, int32_t *ig1, int32_t *ig2, int32_t *ig3, int32_t *ig4, int32_t sz);
void f77name(ez8_gdwfllw)(float *z1, float *z2, double *xlon, int32_t *li, int32_t *lj, char *grtyp, int32_t *ig1, int32_t *ig2, int32_t *ig3, int32_t *ig4, int32_t sz);
void f77name(ez8_vllfxy)(double *dlat, double *dlon, double *x, double *y, int32_t *ni, int32_t *nj, float *d60, float *dgrw, float *pi, float *pj, int32_t *nhem);
void f77name(ez8_vtllfxy)(double *lat, double *lon, double *x, double *y, float *clat, float *clon, float *d60, float *dgrw, int32_t *ni, int32_t *nj, int32_t *n);
void f77name(ez8_vxyfll)(double *x, double *y, double *dlat, double *dlon, int32_t *nb, float *d60, float *dgrw, float *pi, float *pj, int32_t *nhem);
void f77name(ez8_vtxyfll)(double *x, double *y, double *dlat, double *dlon, float *clat, float *clon, float *d60, float *dgrw, int32_t *ni, int32_t *nj, int32_t *nhem);
void f77name(ez8_mxm)(float *a, int32_t *nar, double *b, int32_t *nac, double *c, int32_t *nbc);
void f77name(ez8_llflamb)(double *xlat, double *xlon, double *x, double *y, int32_t *npts, char *grtyp, int32_t *ig1, int32_t *ig2, int32_t *ig3, int32_t *ig4);
void f77name(ez8_lambfll)(double *x, double *y, double *xlat, double *xlon, int32_t *npts, char *grtyp, int32_t *ig1, int32_t *ig2, int32_t *ig3, int32_t *ig4);
void f77name(ez_lambxyfll99)(double *x, double *y, double *lat, double *lon, int32_t *nb, float *latin1, float *latin2, float *yaxislat, float *yaxislon);
void f77name(ez8_lambllfxy99)(double *lat, double *lon, double *x, double *y, int32_t *nb, float *latin1, float *latin2, float *yaxislat, float *yaxislon);
void f77name(ez8_nwtncof)(float *cx, float *cy, double *ax, double *ay, int32_t *ni, int32_t *nj, int32_t *i1, int32_t *i2, int32_t *j1, int32_t *j2, int32_t *extension);
void f77name(igaxg95)(char *gtypout, float *xg, int32_t *nb, char *grtyp, int32_t *ig1, int32_t *ig2, int32_t *ig3, int32_t *ig4);
void f77name(ez_applywgts)(float*, double*, int*, float*, int*, int*, int*, int*, int*, int*, float*);
void f77name(ez_calcpoleval)(float *poleval, float *z, int32_t *ni, double *ax, char *grtyp, char *grref, int32_t c1, int32_t c2);
void f77name(ez_fillnpole)(float *zout, float *zin, int32_t *ni, int32_t *j1, int32_t *j2, float *valpole);
void f77name(ez_fillspole)(float *zout, float *zin, int32_t *ni, int32_t *j1, int32_t *j2, float *valpole);
void f77name(ez_calcxy_y)();
void f77name(ez_calcxy_y_m)();
void f77name(ez_aminmax)(float*, float*, float*, int*, int*);
void f77name(ez_corrbgd)(float*, int*, int*, int*);
void f77name(ez_glat)(double*, double*, int*, int*);
void f77name(ez_crot)(float *r, float *ri, float *lon1, float *lat1, float *lon2, float *lat2);
void f77name(ez_lac)();
void f77name(ez8_uvacart)(double *XYZ, float *U, float *V, double *LON, double *LAT, int32_t *NI, int32_t *NJ);
void f77name(ez8_cartauv)(float *U, float *V, double *UVCART, double *LON, double *LAT, int32_t *NI, int32_t *NJ);
void f77name(ez_xpngdb2)(float* zout, float* zi, int32_t *ni, int32_t *nj, int32_t *j1, int32_t *j2, int32_t *hem, int32_t *symetrie);
void f77name(ez_xpngdag2)(float* zout, float* zi, int32_t *ni, int32_t *nj, int32_t *j1, int32_t *j2, int32_t *hem, int32_t *symetrie);
void f77name(lorenzo_mask_fill)();
void f77name(qqq_ezget_mask_zones)();
void f77name(qqq_ezsint_mask)(char*, double*, double*, int*, int*, char *, int*, int*, TRef_Interp*);

#endif
