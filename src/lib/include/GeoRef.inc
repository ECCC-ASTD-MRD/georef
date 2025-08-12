interface

!  TGeoRef* GeoRef_New();
   function georef_new() result(out) bind(C, name = 'GeoRef_New')
        import :: C_PTR
        implicit none

        type(C_PTR) :: out
    end function georef_new

!   TGeoRef* GeoRef_CreateFromRecord(fst_record *Rec)
    function georef_createfromrecord(rec) result(out) bind(C, name = 'GeoRef_CreateFromRecord')
        import :: C_PTR
        implicit none

        type(C_PTR), intent(in), value  :: rec

        type(C_PTR) :: out
    end function georef_createfromrecord

!   TGeoRef* GeoRef_Create(int32_t NI,int32_t NJ,char *GRTYP,int32_t IG1,int32_t IG2,int32_t IG3,int32_t IG4,fst_file *File);
    function georef_create(ni,nj,grtyp,ig1,ig2,ig3,ig4,file) result(out) bind(C, name = 'GeoRef_Create')
        import :: C_PTR,C_INT32_T,C_CHAR
        implicit none

        integer(C_INT32_T), intent(in) :: ni,nj,ig1,ig2,ig3,ig4
        character(kind = C_CHAR), dimension(*), intent(in) :: grtyp
        type(C_PTR) :: file

        type(C_PTR) :: out
    end function georef_create

!   TGeoRef* GeoRef_CreateU(int32_t NI,int32_t NJ,char *grref,int32_t VerCode,int32_t NbSub,TGeoRef **Subs);
    function georef_createu(ni,nj,grref,vercode,nbsub,subs) result(out) bind(C, name = 'GeoRef_CreateU')
        import :: C_PTR,C_INT32_T,C_CHAR
        implicit none

        integer(C_INT32_T), intent(in) :: ni,nj,nbsub,vercode
        character(kind = C_CHAR), dimension(*), intent(in) :: grref
        type(C_PTR) :: subs

        type(C_PTR) :: out
    end function georef_createu

!   TGeoRef* GeoRef_CreateR(double Lat,double Lon,double Height,int32_t R,double ResR,double ResA);
    function georef_creater(lat,lon,height,r,resr,resa) result(out) bind(C, name = 'GeoRef_CreateR')
        import :: C_PTR,C_DOUBLE
        implicit none

        real(C_DOUBLE), intent(in) :: lat,lon,height,r,resr,resa

        type(C_PTR) :: out
    end function georef_creater

!   TGeoRef* GeoRef_CreateW(int32_t ni,int32_t nj,char *String,double *Transform,double *InvTransform,OGRSpatialReferenceH Spatial);
    function georef_createw(ni,nj,string,transform,invtransform,spatial) result(out) bind(C, name = 'GeoRef_CreateW')
        import :: C_PTR,C_INT32_T,C_CHAR,C_DOUBLE
        implicit none

        integer(C_INT32_T), intent(in) :: ni,nj
        real(C_DOUBLE)    , intent(in), dimension(6) :: transform, invtransform
        type(C_PTR)       , intent(in) :: spatial
        character(kind = C_CHAR), dimension(*), intent(in) :: string

        type(C_PTR) :: out
    end function georef_createw

!   TGeoRef* GeoRef_Define(TGeoRef *Ref,int32_t NI,int32_t NJ,char* GRTYP,char* grref,int32_t IG1,int32_t IG2,int32_t IG3,int32_t IG4,double* AX,double* AY);
    function georef_define(ref,ni,nj,grtyp,grref,ig1,ig2,ig3,ig4,ax,ay) result(out) bind(C, name = 'GeoRef_Define')
        import :: C_PTR,C_INT32_T,C_CHAR,C_DOUBLE
        implicit none

        integer(C_INT32_T), intent(in) :: ni,nj,ig1,ig2,ig3,ig4
        character(kind = C_CHAR), dimension(*), intent(in) :: grtyp, grref
        real(C_DOUBLE), intent(in), dimension(*) ::ax, ay
        type(C_PTR), intent(in), value :: ref

        type(C_PTR) :: out
    end function georef_define

!   TGeoRef* GeoRef_DefineW(TGeoRef *Ref,char *String,double *Transform,double *InvTransform,OGRSpatialReferenceH Spatial);
    function georef_definew(ref,string,transform,invtransform,spatial) result(out) bind(C, name = 'GeoRef_DefineW')
        import :: C_PTR,C_INT32_T,C_CHAR,C_DOUBLE
        implicit none

        real(C_DOUBLE), intent(in), dimension(6) :: transform, invtransform
        type(C_PTR)   , intent(in)    :: spatial
        type(C_PTR)   , intent(in), value :: ref
        character(kind = C_CHAR), dimension(*), intent(in) :: string

        type(C_PTR) :: out
    end function georef_definew

!   TGeoRef* GeoRef_DefineZE(TGeoRef *Ref,int32_t NI,int32_t NJ,float DX,float DY,float LatR,float LonR,int32_t MaxCFL,float XLat1,float XLon1,float XLat2,float XLon2);
    function georef_defineze(ref,ni,nj,dx,dy,latr,lonr,maxcfl,xlat1,xlon1,xlat2,xlon2) result(out) bind(C, name = 'GeoRef_DefineZE')
        import :: C_PTR,C_INT32_T,C_CHAR,C_FLOAT
        implicit none

        integer(C_INT32_T), intent(in) :: ni,nj,maxcfl
        real(C_FLOAT), intent(in) :: dx,dy,latr,lonr,xlat1,xlon1,xlat2,xlon2
        type(C_PTR), intent(in), value :: ref

        type(C_PTR) :: out
    end function georef_defineze

!   TGeoRef* GeoRef_Find(TGeoRef *Ref);
    function georef_find(ref) result(out) bind(C, name = 'GeoRef_Find')
        import :: C_PTR
        implicit none

        type(C_PTR), intent(in), value :: ref

        type(C_PTR) :: out
    end function georef_find

!   TGeoRef* GeoRef_Copy(TGeoRef* __restrict const Ref);
    function georef_copy(ref) result(out) bind(C, name = 'GeoRef_Copy')
        import :: C_PTR
        implicit none

        type(C_PTR), intent(in), value :: ref

        type(C_PTR) :: out
    end function georef_copy

!   TGeoRef* GeoRef_HardCopy(TGeoRef* __restrict const Ref);
    function georef_hardcopy(ref) result(out) bind(C, name = 'GeoRef_HardCopy')
        import :: C_PTR
        implicit none

        type(C_PTR), intent(in), value :: ref

        type(C_PTR) :: out
    end function georef_hardcopy

!   TGeoRef* GeoRef_Reference(TGeoRef* __restrict const Ref);
    function georef_reference(ref) result(out) bind(C, name = 'GeoRef_Reference')
        import :: C_PTR
        implicit none

        type(C_PTR), intent(in), value :: ref

        type(C_PTR) :: out
    end function georef_reference

!   void     GeoRef_Size(TGeoRef *Ref,int32_t X0,int32_t Y0,int32_t X1,int32_t Y1,int32_t BD);
    subroutine georef_size(ref,x0,y0,x1,y1,bd) bind(C, name = 'GeoRef_Size')
        import :: C_PTR, C_INT32_T
        implicit none

        integer(C_INT32_T) :: x0,y0,x1,y1,bd
        type(C_PTR), intent(in), value :: ref

    end subroutine georef_size

!   TGeoRef* GeoRef_Resize(TGeoRef* __restrict const Ref,int32_t NI,int32_t NJ);
    function georef_resize(ref) result(out) bind(C, name = 'GeoRef_Resize')
        import :: C_PTR, C_INT32_T
        implicit none

        integer(C_INT32_T) :: ni,nj
        type(C_PTR), intent(in), value :: ref

        type(C_PTR) :: out
    end function georef_resize

!   int32_t  GeoRef_Free(TGeoRef *Ref);
    function georef_free(ref) result(out) bind(C, name = 'GeoRef_Free')
        import :: C_PTR,C_INT32_T
        implicit none

        type(C_PTR), intent(in), value :: ref

        integer(C_INT32_T) :: out
    end function georef_free

!   void     GeoRef_Clear(TGeoRef *Ref,int32_t New);
    subroutine georef_clear(ref,new) bind(C, name = 'GeoRef_Clear')
        import :: C_PTR, C_INT32_T
        implicit none

        integer(C_INT32_T) :: new
        type(C_PTR), intent(in), value :: ref

    end subroutine georef_clear

!   void     GeoRef_Qualify(TGeoRef* __restrict const Ref);
    subroutine georef_qualify(ref) bind(C, name = 'GeoRef_Qualify')
        import :: C_PTR
        implicit none

        type(C_PTR), intent(in), value :: ref

    end subroutine georef_qualify

!   int32_t  GeoRef_Limits(TGeoRef* __restrict const Ref,double *Lat0,double *Lon0,double *Lat1,double *Lon1);
    function georef_limits(ref,lat0,lon0,lat1,lon1) result(out) bind(C, name = 'GeoRef_Limits')
        import :: C_PTR,C_INT32_T,C_DOUBLE
        implicit none

        real(C_DOUBLE), intent(out) :: lat0,lon0,lat1,lon1
        type(C_PTR)   , intent(in), value :: ref

        integer(C_INT32_T) :: out
    end function georef_limits

!   int32_t  GeoRef_BoundingBox(TGeoRef* __restrict const Ref,double Lat0,double Lon0,double Lat1,double Lon1,double *I0,double *J0,double *I1,double *J1);
    function georef_boundingbox(ref,lat0,lon0,lat1,lon1,i0,j0,i1,j1) result(out) bind(C, name = 'GeoRef_BoundingBox')
        import :: C_PTR,C_INT32_T,C_DOUBLE
        implicit none

        real(C_DOUBLE), intent(in), value :: lat0,lon0,lat1,lon1
        real(C_DOUBLE), intent(out) :: i0,j0,i1,j1
        type(C_PTR)   , intent(in), value :: ref

        integer(C_INT32_T) :: out
    end function georef_boundingbox

!   int32_t GeoRef_Nearest(TGeoRef* __restrict const Ref,double X,double Y,int32_t *Idxs,double *Dists,int32_t NbNear,double MaxDist);
    function georef_nearest(ref,x,y,idxs,dists,nbnear,maxdist) result(out) bind(C, name = 'GeoRef_Nearest')
        import :: C_PTR,C_INT32_T,C_DOUBLE
        implicit none

        real(C_DOUBLE), intent(in), value :: x,y,maxdist
        integer(C_INT32_T), intent(in), dimension(*) :: idxs
        integer(C_INT32_T), intent(in), value :: nbnear
        real(C_DOUBLE), intent(in), dimension(*) :: dists
        type(C_PTR)   , intent(in), value :: ref

        integer(C_INT32_T) :: out
    end function georef_nearest

!   int32_t  GeoRef_Valid(TGeoRef* __restrict const Ref);
    function georef_valid(ref) result(out) bind(C, name = 'GeoRef_Valid')
        import :: C_PTR,C_INT32_T
        implicit none

        type(C_PTR), intent(in), value :: ref

        integer(C_INT32_T) :: out
    end function georef_valid

!   double   GeoRef_XYDistance(TGeoRef *Ref,double X0,double Y0,double X1,double Y1);
    function georef_xydistance(ref,x0,y0,x1,y1) result(out) bind(C, name = 'GeoRef_XYDistance')
        import :: C_PTR,C_DOUBLE
        implicit none

        real(C_DOUBLE), intent(in), value :: x0,y0,x1,y1
        type(C_PTR)   , intent(in), value :: ref

        real(C_DOUBLE) :: out
    end function georef_xydistance

!   double   GeoRef_LLDistance(TGeoRef *Ref,double Lat0,double Lon0,double Lat1,double Lon1);
    function georef_lldistance(ref,lat0,lon0,lat1,lon1) result(out) bind(C, name = 'GeoRef_LLDistance')
        import :: C_PTR,C_DOUBLE
        implicit none

        real(C_DOUBLE), intent(in), value :: lat0,lon0,lat1,lon1
        type(C_PTR)   , intent(in), value :: ref

        real(C_DOUBLE) :: out
    end function georef_lldistance

!   int32_t  GeoRef_WriteFST(TGeoRef *GRef,char *Name,int IG1,int IG2,int IG3,int IG4,fst_file *File);
    function georef_writefst(ref,name,ig1,ig2,ig3,ig4,file) result(out) bind(C, name = 'GeoRef_WriteFST')
        import :: C_PTR,C_CHAR,C_INT32_T
        implicit none

        character(kind = C_CHAR), dimension(*), intent(in) :: name
        type(C_PTR), intent(in), value :: ref,file
        integer(C_INT32_T), intent(in), value :: ig1,ig2,ig3,ig4
        integer(C_INT32_T) :: out

    end function georef_writefst

!   int32_t  GeoRef_Read(struct TGeoRef *GRef);
    function georef_read(ref) result(out) bind(C, name = 'GeoRef_Read')
        import :: C_PTR,C_INT32_T
        implicit none

        type(C_PTR), intent(in), value :: ref

        integer(C_INT32_T) :: out
    end function georef_read

!   int32_t  GeoRef_CopyDesc(fst_file *FileTo,fst_record* Rec);
    function georef_copydesc(file,rec) result(out) bind(C, name = 'GeoRef_CopyDesc')
        import :: C_PTR,C_INT32_T
        implicit none

        type(C_PTR), intent(in) :: file,rec

        integer(C_INT32_T) :: out
    end function georef_copydesc

!   int32_t      GeoRef_GetLL(TGeoRef *Ref,double *Lat,double *Lon);                                                                      // gdll / GeoRef_Coords
    function georef_getll(ref,lat,lon) result(out) bind(C, name = 'GeoRef_GetLL')
        import :: C_PTR,C_INT32_T, C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: ref
        real(C_DOUBLE), intent(in), dimension(*) :: lat,lon

        integer(C_INT32_T) :: out
    end function georef_getll

!   int32_t  GeoRef_CellDims(TGeoRef *Ref,int32_t Invert,float* DX,float* DY,float* DA);
    function georef_celldims(ref,invert,dx,dy,da) result(out) bind(C, name = 'GeoRef_CellDims')
        import :: C_PTR,C_INT32_T,C_FLOAT
        implicit none
        type(C_PTR), intent(in), value :: ref
        integer(C_INT32_T), intent(in) :: invert
        real(C_FLOAT), intent(in), dimension(*) :: dx,dy,da

        integer(C_INT32_T) :: out
    end function georef_celldims

!   int32_t      GeoRef_XY2LL(TGeoRef *Ref,double *Lat,double *Lon,double *X,double *Y,int32_t N,int32_t Extrap);                                 // c_gdllfxy
    function georef_xy2ll(ref,lat,lon,x,y,n,extrap) result(out) bind(C, name = 'GeoRef_XY2LL')
        import :: C_PTR,C_INT32_T, C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: ref
        real(C_DOUBLE), intent(in), dimension(*) :: lat,lon,x,y
        integer(C_INT32_T), value :: n,extrap

        integer(C_INT32_T) :: out
    end function georef_xy2ll

!   int32_t      GeoRef_LL2XY(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int32_t N,int32_t Extrap);                                 // c_gdxyfll
    function georef_ll2xy(ref,x,y,lat,lon,n,extrap) result(out) bind(C, name = 'GeoRef_LL2XY')
        import :: C_PTR,C_INT32_T, C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: ref
        real(C_DOUBLE), intent(in), dimension(*) :: lat,lon,x,y
        integer(C_INT32_T), value :: n,extrap

        integer(C_INT32_T) :: out
    end function georef_ll2xy

!   int32_t      GeoRef_XYVal(TGeoRef *Ref,TGeoOptions *Opt,float *zout,float *zin,double *x,double *y,int32_t n);                                             // c_gdxysval
    function georef_xyval(ref,opt,zout,zin,x,y,n) result(out) bind(C, name = 'GeoRef_XYVal')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: ref
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(in), dimension(*) :: zout
        real(C_FLOAT), intent(in), dimension(*) :: zin
        real(C_DOUBLE), intent(in), dimension(*) :: x,y
        integer(C_INT32_T),value :: n

        integer(C_INT32_T) :: out
    end function georef_xyval

!   int32_t      GeoRef_XYUVVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *x,double *y,int32_t n);                // c_gdxyvval
    function georef_xyuvval(ref,opt,uuout,vvout,uuin,vvin,x,y,n) result(out) bind(C, name = 'GeoRef_XYUVVal')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: ref
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin
        real(C_DOUBLE), intent(out), dimension(*) :: x,y
        integer(C_INT32_T), value :: n

        integer(C_INT32_T) :: out
    end function georef_xyuvval

!   int32_t      GeoRef_XYWDVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *x,double *y,int32_t n);                // c_gdxywdval
    function georef_xywdval(ref,opt,uuout,vvout,uuin,vvin,x,y,n) result(out) bind(C, name = 'GeoRef_XYWDVal')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: ref
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin
        real(C_DOUBLE), intent(out), dimension(*) :: x,y
        integer(C_INT32_T), value :: n

        integer(C_INT32_T) :: out
    end function georef_xywdval

!   int32_t      GeoRef_LLVal(TGeoRef *Ref,TGeoOptions *Opt,float *zout,float *zin,double *lat,double *lon,int32_t n);                                         // c_gdllsval
    function georef_llval(ref,opt,zout,zin,lat,lon,n) result(out) bind(C, name = 'GeoRef_LLVal')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: ref
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(in), dimension(*) :: zout
        real(C_FLOAT), intent(in), dimension(*) :: zin
        real(C_DOUBLE), intent(in), dimension(*) :: lat,lon
        integer(C_INT32_T), value :: n

        integer(C_INT32_T) :: out
    end function georef_llval

!   int32_t      GeoRef_LLUVVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int32_t n);            // c_gdllvval
    function georef_lluvval(ref,opt,uuout,vvout,uuin,vvin,lat,lon,n) result(out) bind(C, name = 'GeoRef_LLUVVal')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: ref
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin
        real(C_DOUBLE), intent(out), dimension(*) :: lat,lon
        integer(C_INT32_T), value :: n

        integer(C_INT32_T) :: out
    end function georef_lluvval
    
!   int32_t      GeoRef_LLWDVal(TGeoRef *Ref,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int32_t n);            // c_gdllwdval
    function georef_llwdval(ref,opt,uuout,vvout,uuin,vvin,lat,lon,n) result(out) bind(C, name = 'GeoRef_LLWDVal')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: ref
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin
        real(C_DOUBLE), intent(out), dimension(*) :: lat,lon
        integer(C_INT32_T), value :: n

        integer(C_INT32_T) :: out
    end function georef_llwdval

!   int32_t  GeoRef_Interp(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout, float *zin);                                                  // c_ezsint
    function georef_interp(refto,reffrom,opt,zout,zin) result(out) bind(C, name = 'GeoRef_Interp')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: refto,reffrom
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(out), dimension(*) :: zout
        real(C_FLOAT), intent(in), dimension(*) :: zin

        integer(C_INT32_T) :: out
    end function georef_interp

!   int32_t GeoRef_InterpDef(TGeoRef *ToRef, TDef *ToDef, TGeoRef *FromRef, TDef *FromDef, TGeoOptions *Opt, int32_t Final);
    function georef_interpdef(refto,reffrom,defto,deffrom,opt,final) result(out) bind(C, name = 'GeoRef_InterpDef')
        import :: C_PTR,C_INT32_T,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: refto,reffrom
        type(C_PTR), intent(in), value :: defto,deffrom
        type(C_PTR), intent(in), value :: opt
        integer(C_INT32_T), value :: final

        integer(C_INT32_T) :: out
    end function georef_interpdef

!   int32_t  GeoRef_InterpUV(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin);                      // c_ezuvint
    function georef_interpuv(refto,reffrom,opt,uuout,vvout,uuin,vvin) result(out) bind(C, name = 'GeoRef_InterpUV')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: refto,reffrom
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin

        integer(C_INT32_T) :: out
    end function georef_interpuv

!   int32_t  GeoRef_InterpWD(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin);                      // c_ezwdint
    function georef_interpwd(refto,reffrom,opt,uuout,vvout,uuin,vvin) result(out) bind(C, name = 'GeoRef_InterpWD')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: refto,reffrom
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin

        integer(C_INT32_T) :: out
    end function georef_interpwd

!   int32_t  GeoRef_InterpYY(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *zout, float *zin);                                                // c_ezyysint
    function georef_interpyy(refto,reffrom,opt,zout,zin) result(out) bind(C, name = 'GeoRef_InterpYY')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: refto,reffrom
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(out), dimension(*) :: zout
        real(C_FLOAT), intent(in), dimension(*) :: zin

        integer(C_INT32_T) :: out
    end function georef_interpyy

!   int32_t  GeoRef_InterpYYUV(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin);                    // c_ezyyuvint
    function georef_interpyyuv(refto,reffrom,opt,uuout,vvout,uuin,vvin) result(out) bind(C, name = 'GeoRef_InterpYYUV')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: refto,reffrom
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin

        integer(C_INT32_T) :: out
    end function georef_interpyyuv

!   int32_t  GeoRef_InterpYYWD(TGeoRef *RefTo,TGeoRef *RefFrom,TGeoOptions *Opt,float *uuout,float *vvout,float *uuin,float *vvin);                    // c_ezyywdint
    function georef_interpyywd(refto,reffrom,opt,uuout,vvout,uuin,vvin) result(out) bind(C, name = 'GeoRef_InterpYYWD')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: refto,reffrom
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin

        integer(C_INT32_T) :: out
    end function georef_interpyywd

!   int32_t  GeoRef_InterpMask(TGeoRef *RefTo, TGeoRef *RefFrom,TGeoOptions *Opt,char *MaskOut,char *MaskIn);
    function georef_interpmask(refto,reffrom,opt,mout,min) result(out) bind(C, name = 'GeoRef_InterpMask')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: refto,reffrom
        type(C_PTR), intent(in), value :: opt
        real(C_FLOAT), intent(out), dimension(*) :: mout
        real(C_FLOAT), intent(in), dimension(*) :: min

        integer(C_INT32_T) :: out
    end function georef_interpmask

!   int32_t  GeoRef_WD2UV(TGeoRef *Ref,float *uugdout,float *vvgdout,float *uullin,float *vvllin,double *Lat,double *Lon,int32_t npts);   // c_gduvfwd
    function georef_wd2uv(ref,uuout,vvout,spdin,dirin,lat,lon,npts) result(out) bind(C, name = 'GeoRef_WD2UV')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: ref
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: spdin,dirin
        real(C_DOUBLE), intent(in), dimension(*) :: lat,lon
        integer(C_INT32_T), value :: npts

        integer(C_INT32_T) :: out
    end function georef_wd2uv

!   int32_t  GeoRef_UV2WD(TGeoRef *Ref,float *spd_out,float *wd_out,float *uuin,float *vvin,double *Lat,double *Lon,int32_t npts);        // c_gdwdfuv
    function georef_uv2wd(ref,spdout,dirout,uuin,vvin,lat,lon,npts) result(out) bind(C, name = 'GeoRef_UV2WD')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: ref
        real(C_FLOAT), intent(out), dimension(*) :: spdout,dirout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin
        real(C_DOUBLE), intent(in), dimension(*) :: lat,lon
        integer(C_INT32_T), value :: npts

        integer(C_INT32_T) :: out
    end function georef_uv2wd

!   int32_t  GeoRef_UV2UV(TGeoRef *Ref,float *uullout,float *vvllout,float *uuin,float *vvin,double *Lat,double *Lon,int32_t Nb);    //c_gdlluvfuv_orig
    function georef_uv2uv(ref,uuout,vvout,uuin,vvin,lat,lon,npts) result(out) bind(C, name = 'GeoRef_UV2UV')
        import :: C_PTR,C_INT32_T, C_FLOAT,C_DOUBLE
        implicit none

        type(C_PTR), intent(in), value :: ref
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin
        real(C_DOUBLE), intent(in), dimension(*) :: lat,lon
        integer(C_INT32_T), value :: npts

        integer(C_INT32_T) :: out
    end function georef_uv2uv

!   int32_t  GeoRef_GridGetParams(TGeoRef *Ref,int32_t *NI,int32_t *NJ,char *GRTYP,int32_t *IG1,int32_t *IG2,int32_t *IG3,int32_t *IG4,char *grref,int32_t *IG1REF,int32_t *IG2REF,int32_t *IG3REF,int32_t *IG4REF);  //c_ezgxprm
    function georef_gridgetparams(ref,ni,nj,grtyp,ig1,ig2,ig3,ig4,grref,ig1ref,ig2ref,ig3ref,ig4ref) result(out) bind(C, name = 'GeoRef_GridGetParams')
        import :: C_PTR,C_INT32_T,C_CHAR,C_DOUBLE
        implicit none

        integer(C_INT32_T), intent(out) :: ni,nj,ig1,ig2,ig3,ig4,ig1ref,ig2ref,ig3ref,ig4ref
        character(kind = C_CHAR), dimension(*), intent(out) :: grtyp, grref
        type(C_PTR), intent(in), value :: ref

        type(C_PTR) :: out
    end function georef_gridgetparams

!    int32_t GeoRef_Within(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1)
    function georef_within(gref0,gref1) result(out) bind(C, name = 'GeoRef_Within')
        import :: C_PTR,C_INT32_T
        implicit none
        type(C_PTR), intent(in), value :: gref0,gref1

        integer(C_INT32_T) :: out
    end function georef_within

!   int32_t  GeoRef_WithinRange(TGeoRef* __restrict const Ref,double Lat0,double Lon0,double Lat1,double Lon1,int32_t In);
    function georef_withinrange(gref0,lat0,lon0,lat1,lon1,in) result(out) bind(C, name = 'GeoRef_WithinRange')
        import :: C_PTR,C_DOUBLE,C_INT32_T
        implicit none
        type(C_PTR), intent(in), value :: gref0
        real(C_DOUBLE), intent(in), value :: lat0,lon0,lat1,lon1
        integer(C_INT32_T), intent(in), value :: in

        integer(C_INT32_T) :: out
    end function georef_withinrange

!   int32_t  GeoRef_Intersect(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1,int32_t *X0,int32_t *Y0,int32_t *X1,int32_t *Y1,int32_t BD);
    function georef_intersect(gref0,gref1,x0,y0,x1,y1,bd) result(out) bind(C, name = 'GeoRef_Intersect')
        import :: C_PTR,C_INT32_T
        implicit none
        type(C_PTR), intent(in), value :: gref0,gref1
        integer(C_INT32_T), intent(out) :: x0,x1,y0,y1
        integer(C_INT32_T), intent(in) :: bd

        integer(C_INT32_T) :: out
    end function georef_intersect

!   int32_t  GeoRef_Equal(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1);
    function georef_equal(gref0,gref1) result(out) bind(C, name = 'GeoRef_Equal')
        import :: C_PTR,C_INT32_T
        implicit none
        type(C_PTR), intent(in), value :: gref0,gref1

        integer(C_INT32_T) :: out
    end function georef_equal

!   int32_t  GeoSet_Free(TGeoSet *Set);
    function geoset_free(set) result(out) bind(C, name = 'GeoSet_Free')
        import :: C_PTR,C_INT32_T
        implicit none

        type(C_PTR), intent(in), value :: set

        integer(C_INT32_T) :: out
    end function geoset_free

!   TGeoSet* GeoRef_SetGet(TGeoRef* RefTo, TGeoRef* RefFrom,TGeoOptions *Opt)
    function georef_setget(refto,reffrom,opt) result(out) bind(C, name = 'GeoRef_SetGet')
        import :: C_PTR
        implicit none
        type(C_PTR), intent(in), value :: refto,reffrom,opt

        type(C_PTR) :: out
    end function georef_setget

!   TGeoSet* GeoRef_SetReadFST(TGeoSet *GSet,const int32_t InterpType,const fst_file * const File)
    function georef_setreadfst(gset,interp,file) result(out) bind(C, name = 'GeoRef_SetReadFST')
        import :: C_PTR, C_INT32_T
        implicit none
        type(C_PTR), intent(in), value :: gset,file
        integer(C_INT32_T), value :: interp

        type(C_PTR) :: out
    end function georef_setreadfst

!   int32_t GeoRef_SetWriteFST(TGeoSet *GSet,fst_file *File){
    function geoset_writefst(gset,file) result(out) bind(C, name = 'GeoRef_SetWriteFST')
        import :: C_PTR,C_INT32_T
        implicit none
        type(C_PTR), intent(in), value :: gset,file

        integer(C_INT32_T) :: out
    end function geoset_writefst

!   TDef *Def_Create(int32_t NI,int32_t NJ,int32_t NK,TDef_Type Type,char* Comp0,char* Comp1,char* Mask);
    function def_create(ni,nj,nk,type,comp0,comp1,mask) result(def) bind(C, name = 'Def_Create')
       import :: C_PTR,C_INT32_T
       implicit none
       type(C_PTR), intent(in), value :: comp0,comp1,mask
       integer(C_INT32_T), value :: ni,nj,nk,type

       type(C_PTR) :: def 
    end function def_create

end interface