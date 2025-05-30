// EZSCINT merged fonctionnalities
int32_t  GeoRef_DefRPNXG(TGeoRef* Ref);                                                                                                       // c_ezdefxg

int32_t c_gdll(
    //! [in] Grid identifier (returned by ezqkdef or ezgdef)
    int32_t gdid,
    //! [out] Array of grid point latitudes. Must be allocated by client with a dimension corresponding to the grid.
    float *lat,
    //! [out] Array of grid point longitudes. Must be allocated by client with a dimension corresponding to the grid.
    float *lon
) {
    int32_t nb;
    double *lat,*lon;
    nb=Ref->NX*Ref->NY;
    int32_t  GeoRef_GetLL(TGeoRef *Ref,double *Lat,double *Lon);                                                                                  // gdll
}




int32_t  GeoRef_CalcLL(TGeoRef* Ref);                                                                                                         // ez_calclatlon
int32_t  GeoRef_XY2LL(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t N,int32_t Extrap);                                     // c_gdllfxy
int32_t  GeoRef_LL2XY(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t N,int32_t Extrap);                                     // c_gdxyfll

int32_t  GeoRef_XYVal(TGeoRef *Ref,TGeoOptions *Opt,float *zout,float *zin,double *x,double *y,int32_t n);                                    // c_gdxysval
int32_t  GeoRef_XYUVVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *x,double *y,int32_t n);       // c_gdxyvval
int32_t  GeoRef_XYWDVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *x,double *y,int32_t n);       // c_gdxywdval
int32_t  GeoRef_LLVal(TGeoRef *Ref,TGeoOptions *Opt,float *zout,float *zin,double *lat,double *lon,int32_t n);                                // c_gdllsval
int32_t  GeoRef_LLUVVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int32_t n);   // c_gdllvval
int32_t  GeoRef_LLWDVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int32_t n);   // c_gdllwdval
int32_t  GeoRef_Interp(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout, float *zin);                                             // c_ezsint
int32_t  GeoRef_InterpUV(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin);                 // c_ezuvint
int32_t  GeoRef_InterpWD(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin);                 // c_ezwdint
int32_t  GeoRef_InterpYY(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout, float *zin);                                           // c_ezyysint
int32_t  GeoRef_InterpYYUV(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin);               // c_ezyyuvint
int32_t  GeoRef_InterpYYWD(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin);               // c_ezyywdint
int32_t  GeoRef_InterpMask(TGeoRef *RefTo, TGeoRef *RefFrom,TGeoOptions *Opt,char *MaskOut,char *MaskIn);
int32_t  GeoRef_WD2UV(TGeoRef *Ref,float *uugdout,float *vvgdout,float *uullin,float *vvllin,double *Lat,double *Lon,int32_t npts);           // c_gduvfwd
int32_t  GeoRef_UV2WD(TGeoRef *Ref,float *spd_out,float *wd_out,float *uuin,float *vvin,double *Lat,double *Lon,int32_t npts);                // c_gdwdfuv
int32_t  GeoRef_UV2UV(TGeoRef *Ref,float *uullout,float *vvllout,float *uuin,float *vvin,double *Lat,double *Lon,int32_t Nb);                 // c_gdlluvfuv

int32_t  GeoRef_MaskZones(TGeoRef *RefTo,TGeoRef *RefFrom,int32_t *MaskOut,int32_t *MaskIn);
int32_t  GeoRef_MaskYYApply(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,int32_t ni,int32_t nj,float *maskout,double *dlat,double *dlon,double *yinlat,double *yinlon,int32_t *yyincount,double *yanlat,double *yanlon,int32_t *yyancount);
int32_t  GeoRef_MaskYYDefine(TGeoRef *Ref);

void     GeoRef_GridGetExpanded(TGeoRef *Ref,TGeoOptions *Opt,float *zout,float *zin);                                                         // gdxpngd
int32_t  GeoRef_GridGetParams(TGeoRef *Ref,int32_t *NI,int32_t *NJ,char *GRTYP,int32_t *IG1,int32_t *IG2,int32_t *IG3,int32_t *IG4,char *grref,int32_t *IG1REF,int32_t *IG2REF,int32_t *IG3REF,int32_t *IG4REF);  //c_ezgxprm
int32_t  GeoRef_AxisGet(TGeoRef *Ref,double *AX,double *AY);                                                                                   // gdaxes
int32_t  GeoRef_AxisGetExpanded(TGeoRef* Ref,double *AX,double *AY);                                                                           // gdgxpndaxes