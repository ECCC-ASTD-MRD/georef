module georef_mod
    use App
    use f_c_strings_mod
    use rmn_common
    use rmn_fst24
    implicit none

    include 'GeoRef.inc'

    enum, bind(C)
        enumerator :: &
            IR_UNDEF                          = 0, &
            IR_NEAREST                        = 1, &
            IR_LINEAR                         = 2, &
            IR_CUBIC                          = 3, &
            IR_NORMALIZED_CONSERVATIVE        = 4, &
            IR_CONSERVATIVE                   = 5, &
            IR_MAXIMUM                        = 6, &
            IR_MINIMUM                        = 7, &
            IR_SUM                            = 8, &
            IR_AVERAGE                        = 9, &
            IR_VARIANCE                       = 10, &
            IR_SQUARE                         = 11, &
            IR_NORMALIZED_COUNT               = 12, &
            IR_COUNT                          = 13, &
            IR_VECTOR_AVERAGE                 = 14, &
            IR_NOP                            = 15, &
            IR_ACCUM                          = 16, &
            IR_BUFFER                         = 17, &
            IR_SUBNEAREST                     = 18, &
            IR_SUBLINEAR                      = 19
        enumerator :: &
            ER_UNDEF   = 0, &
            ER_MAXIMUM = 1, &
            ER_MINIMUM = 2, &
            ER_VALUE   = 3, &
            ER_ABORT   = 4
        enumerator :: &
            IV_UNDEF                          = 0, &
            IV_FAST                           = 1, &
            IV_WITHIN                         = 2, &
            IV_INTERSECT                      = 3, &
            IV_CENTROID                       = 4, &
            IV_ALIASED                        = 5, &
            IV_CONSERVATIVE                   = 6, &
            IV_NORMALIZED_CONSERVATIVE        = 7, &
            IV_POINT_CONSERVATIVE             = 8, &
            IV_LENGTH_CONSERVATIVE            = 9, &
            IV_LENGTH_NORMALIZED_CONSERVATIVE = 10, &
            IV_LENGTH_ALIASED                 = 11
        enumerator :: &
            CB_REPLACE   = 0, &
            CB_MIN       = 1, &
            CB_MAX       = 2, &
            CB_SUM       = 3, &
            CB_AVERAGE   = 4
    end enum

    type :: georef
        private
        type(C_PTR) :: ptr = c_null_ptr ! Pointer to C control structure
    contains
        procedure, pass   :: init => georef_init_f
        procedure, pass   :: copy => georef_copy_f
        procedure, pass   :: valid => georef_valid_f
        procedure, pass   :: withinrange => georef_withinrange_f
        procedure, pass   :: intersect => georef_intersect_f
        procedure, pass   :: limits => georef_limits_f
        procedure, pass   :: boundingbox => georef_boundingbox_f
        procedure, pass   :: griddistance => georef_griddistance_f
        procedure, pass   :: write => georef_write_f
        procedure, pass   :: read => georef_read_f
        procedure, pass   :: interp => georef_interp_f
        procedure, pass   :: interpuv => georef_interpuv_f
        procedure, pass   :: interpwd => georef_interpwd_f
        procedure, pass   :: ud2wd => georef_uv2wd_f
        procedure, pass   :: wd2uv => georef_wd2uv_f
        procedure, pass   :: uv2uv => georef_uv2uv_f
        procedure, pass   :: llwdval => georef_llwdval_f
        procedure, pass   :: lluvval => georef_lluvval_f
        procedure, pass   :: llval => georef_llval_f
        procedure, pass   :: xywdval => georef_xywdval_f
        procedure, pass   :: xyuvval => georef_xyuvval_f
        procedure, pass   :: xyval => georef_xyval_f
        procedure, pass   :: ll2xy => georef_ll2xy_f
        procedure, pass   :: xy2ll => georef_xy2ll_f
        procedure, pass   :: getll => georef_getll_f
        final :: georef_finalize
    end type georef

    type, bind(C) :: geooptions
        type(C_PTR) :: ptr = c_null_ptr              !< Pointer to C control structure
        integer(C_INT32_T) :: Interp = IR_CUBIC      !< Interpolation degree
        integer(C_INT32_T) :: Extrap = ER_MAXIMUM    !< Extrapolation method
        integer(C_INT32_T) :: InterpVector = IV_FAST !< Vector interpolation method
        integer(C_INT32_T) :: Combine = CB_REPLACE   !< Aggregation type
        integer(C_INT32_T) :: Transform = 1          !< Apply transformation or stay within master referential
        integer(C_INT32_T) :: CIndex = 0             !< C Indexing (starts st 0)
        integer(C_INT32_T) :: Symmetric = 0          !< 
        integer(C_INT32_T) :: Segment = 1            !< How much segmentation (Conservatives/Geometric modes)
        integer(C_INT32_T) :: Sampling = 1           !< Sampling interval
        integer(C_INT8_T)  :: PolarCorrect = 1       !< Apply polar corrections
        integer(C_INT8_T)  :: VectorMode = 0         !< Process data as vector
        real(C_FLOAT)      :: DistTreshold = 10.0    !< Distance treshold for point32_t clouds
        real(C_FLOAT)      :: LonRef = 0.0           !< Longitude referential (-180.0,0.0)
        real(C_FLOAT)      :: NoData = 0.0           !< NoData Value (Default: NaN)
        type(C_PTR) :: Table = c_null_ptr            !< Data table to check of values to check for
        type(C_PTR) :: lutDef = c_null_ptr           !< Lookup table
        integer(C_INT32_T) :: lutSize = 0            !< Number of lookup elements
        integer(C_INT32_T) :: lutDim = 0             !< Dimension of the lookup elements
        type(C_PTR) :: Ancilliary = c_null_ptr       !< PPre calculated field (ex: variance, average,...)
    end type geooptions

contains

    function georef_init_f(this) result(res)
        implicit none
        class(georef), intent(inout) :: this       !< georef instance

        logical :: res                             !< Whether the georef is initialized

        res=.false.
        this%ptr = georef_new()
        if (c_associated(this%ptr)) then
           res=.true.
        endif
    end function georef_init_f

    function georef_valid_f(this) result(res)
        implicit none
        class(georef), intent(inout) :: this       !< georef instance
        logical :: res                             !< Whether the georef is initialized

        integer(C_INT32_T) val

        res=.false.
        val = georef_valid(this%ptr)
        if (val == 1) then
           res=.true.
        endif
    end function georef_valid_f

    function georef_copy_f(this,hard) result(out)
        implicit none
        class(georef), intent(inout) :: this       !< georef instance
        class(georef), allocatable :: out          !< georef instance
        logical, optional :: hard

        if (present(hard)) then
           out%ptr=georef_copy(this%ptr)
        else
           out%ptr=georef_hardcopy(this%ptr)
        endif
    end function georef_copy_f

    subroutine georef_finalize(ref)
        implicit none
        type(georef), intent(inout) :: ref  !< georef instance

        integer(C_INT32_T) val

        val = georef_free(ref%ptr)
        ref%ptr=C_NULL_PTR
    end subroutine georef_finalize

   !> \copybrief georef_within
    function georef_within_f(this,ref) result(res)
        implicit none
        class(georef), intent(in) :: this, ref  !< georef instance

        integer(C_INT32_T) val
        logical :: res                             !< Whether the georef is included within georef objet

        res=.false.;
        val = georef_within(this%ptr,ref%ptr)
        if (val==1) then
           res=.true.
        endif
    end function georef_within_f

   !> \copybrief georef_within
    function georef_withinrange_f(this,lat0,lon0,lat1,lon1,in) result(res)
        class(georef), intent(in) :: this  !< georef instance
        real(C_DOUBLE), intent(in) :: lat0,lon0,lat1,lon1
        logical :: in

        integer(C_INT32_T) val,cin
        logical :: res                             !< Whether the georef is included within range

        res=.false.;
        cin=0;
        if (in) cin=1
        val = georef_withinrange(this%ptr,lat0,lon0,lat1,lon1,cin)
        if (val==1) then
           res=.true.
        endif
    end function georef_withinrange_f

   function georef_intersect_f(this,ref,x0,y0,x1,y1,bd) result(res)
        class(georef), intent(in) :: this, ref  !< georef instance
        integer(C_INT32_T) :: x0,x1,y0,y1,bd

        integer(C_INT32_T) val
        logical :: res                             !< Whether the georef are intersecting

        val=georef_intersect(this%ptr,ref%ptr,x0,y0,x1,y1,bd)
        if (val==1) then
           res=.true.
        endif
      end function georef_intersect_f

    function georef_limits_f(this,lat0,lon0,lat1,lon1) result(res)
        class(georef),  intent(in) :: this  !< georef instance
        real(C_DOUBLE), intent(out) :: lat0,lon0,lat1,lon1
        integer(C_INT32_T) :: val

        logical :: res

        res=.false.;
        val=georef_limits(this%ptr,lat0,lon0,lat1,lon1)
        if (val==1) then
           res=.true.
        endif
    end function georef_limits_f

    function georef_boundingbox_f(this,lat0,lon0,lat1,lon1,i0,j0,i1,j1) result(res)
        class(georef),  intent(in) :: this  !< georef instance
        real(C_DOUBLE), intent(in) :: lat0,lon0,lat1,lon1
        real(C_DOUBLE), intent(out) :: i0,j0,i1,j1
        integer(C_INT32_T) :: val

        logical :: res

        res=.false.;
        val=georef_boundingbox(this%ptr,lat0,lon0,lat1,lon1,i0,j0,i1,j1)
        if (val==1) then
           res=.true.
        endif
    end function georef_boundingbox_f

    function georef_griddistance_f(this,x0,y0,x1,y1) result(out)
        class(georef),  intent(in) :: this  !< georef instance
        real(C_DOUBLE), intent(out) :: x0,y0,x1,y1

        real(C_DOUBLE) :: out

        out=georef_griddistance(this%ptr,x0,y0,x1,y1)
    end function georef_griddistance_f

    function georef_write_f(this,name,file) result(res)
        class(georef),  intent(in) :: this  !< georef instance
        character(len=*), intent(in) :: name
        type(fst_file), intent(in) :: file

        integer(C_INT32_T) :: val
        logical :: res

        res=.false.;
        val=georef_write(this%ptr,trim(name)//C_NULL_CHAR,file%get_c_ptr())
        if (val==1) then
           res=.true.
        endif
    end function georef_write_f

    function georef_read_f(this) result(res)
        class(georef),  intent(in) :: this  !< georef instance

        integer(C_INT32_T) :: val
        logical :: res

        res=.false.;

        val=georef_read(this%ptr)
        if (val==1) then
           res=.true.
        endif
    end function georef_read_f

!   c_gdll
    function georef_getll_f(this,lat,lon) result(out)
        class(georef),  intent(inout) :: this  !< georef instance
        real(C_DOUBLE), intent(out), dimension(*) :: lat,lon

        integer(C_INT32_T) :: out

        out=georef_getll(this%ptr,lat,lon)
    end function georef_getll_f

!   c_gdllfxy
    function georef_xy2ll_f(this,lat,lon,x,y,n,extrap) result(out)
        class(georef),  intent(inout) :: this  !< georef instance
        real(C_DOUBLE), intent(out), dimension(*) :: lat,lon
        real(C_DOUBLE), intent(in), dimension(*) :: x,y
        integer(C_INT32_T) :: n,extrap

        integer(C_INT32_T) :: out

        out=georef_xy2ll(this%ptr,lat,lon,x,y,n,extrap)
    end function georef_xy2ll_f
 
!   c_gdxyfll
    function georef_ll2xy_f(this,x,y,lat,lon,n,extrap) result(out)
        class(georef),  intent(inout) :: this  !< georef instance
        real(C_DOUBLE), intent(in), dimension(*) :: lat,lon
        real(C_DOUBLE), intent(out), dimension(*) :: x,y
        integer(C_INT32_T) :: n,extrap

        integer(C_INT32_T) :: out

        out=georef_ll2xy(this%ptr,x,y,lat,lon,n,extrap)
    end function georef_ll2xy_f
 
!   c_gdxysval
    function georef_xyval_f(this,opt,zout,zin,x,y,n) result(out)
        class(georef),  intent(inout) :: this  !< georef instance
        type(geooptions), intent(in), target :: opt  !< georef instance
        real(C_FLOAT), intent(in), dimension(*) :: zin
        real(C_FLOAT), intent(out), dimension(*) :: zout
        real(C_DOUBLE), intent(out), dimension(*) :: x,y
        integer(C_INT32_T) :: n

        integer(C_INT32_T) :: out

        out=georef_xyval(this%ptr,C_LOC(opt),zout,zin,x,y,n)
    end function georef_xyval_f

!   c_gdxyvval
    function georef_xyuvval_f(this,opt,uuout,vvout,uuin,vvin,x,y,n) result(out)
        class(georef),  intent(inout) :: this  !< georef instance
        type(geooptions),  intent(in), target :: opt  !< georef instance
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin
        real(C_DOUBLE), intent(out), dimension(*) :: x,y
        integer(C_INT32_T) :: n

        integer(C_INT32_T) :: out

        out=georef_xyuvval(this%ptr,C_LOC(opt),uuout,vvout,uuin,vvin,x,y,n)
    end function georef_xyuvval_f

!   c_gdxywdval
    function georef_xywdval_f(this,opt,uuout,vvout,uuin,vvin,x,y,n) result(out)
        class(georef),  intent(inout) :: this  !< georef instance
        type(geooptions),  intent(in), target :: opt  !< georef instance
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin
        real(C_DOUBLE), intent(out), dimension(*) :: x,y
        integer(C_INT32_T) :: n

        integer(C_INT32_T) :: out

        out=georef_xywdval(this%ptr,C_LOC(opt),uuout,vvout,uuin,vvin,x,y,n)
    end function georef_xywdval_f

!   c_gdllsval
    function georef_llval_f(this,opt,zout,zin,lat,lon,n) result(out)
        class(georef),  intent(inout) :: this  !< georef instance
        type(geooptions),  intent(in), target :: opt  !< georef instance
        real(C_FLOAT), intent(out), dimension(*) :: zout
        real(C_FLOAT), intent(in), dimension(*) :: zin
        real(C_DOUBLE), intent(out), dimension(*) :: lat,lon
        integer(C_INT32_T) :: n

        integer(C_INT32_T) :: out

        out=georef_llval(this%ptr,C_LOC(opt),zout,zin,lat,lon,n)
    end function georef_llval_f
  
!   c_gdllvval
    function georef_lluvval_f(this,opt,uuout,vvout,uuin,vvin,lat,lon,n) result(out)
        class(georef),  intent(inout) :: this  !< georef instance
        type(geooptions),  intent(in), target :: opt  !< georef instance
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin
        real(C_DOUBLE), intent(out), dimension(*) :: lat,lon
        integer(C_INT32_T) :: n

        integer(C_INT32_T) :: out

         out=georef_lluvval(this%ptr,C_LOC(opt),uuout,vvout,uuin,vvin,lat,lon,n)
   end function georef_lluvval_f
    
!   c_gdllwdval
    function georef_llwdval_f(this,opt,uuout,vvout,uuin,vvin,lat,lon,n) result(out)
        class(georef),  intent(inout) :: this  !< georef instance
        type(geooptions),  intent(in), target :: opt  !< georef instance
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin
        real(C_DOUBLE), intent(out), dimension(*) :: lat,lon
        integer(C_INT32_T) :: n

        integer(C_INT32_T) :: out

        out=georef_llwdval(this%ptr,C_LOC(opt),uuout,vvout,uuin,vvin,lat,lon,n)
   end function georef_llwdval_f

!   c_ezsint
    function georef_interp_f(this,reffrom,opt,zout,zin) result(out)
        class(georef),  intent(inout) :: this,reffrom  !< georef instance
        type(geooptions),  intent(in), target :: opt  !< georef instance
        real(C_FLOAT), intent(out), dimension(*) :: zout
        real(C_FLOAT), intent(in), dimension(*) :: zin

        integer(C_INT32_T) :: out

        out=georef_interp(this%ptr,reffrom%ptr,C_LOC(opt),zout,zin)
    end function georef_interp_f
 
!   c_ezuvint
    function georef_interpuv_f(this,reffrom,opt,uuout,vvout,uuin,vvin) result(out)
        class(georef),  intent(inout) :: this,reffrom  !< georef instance
        type(geooptions),  intent(in), target :: opt  !< georef instance
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin

        integer(C_INT32_T) :: out

        out=georef_interpuv(this%ptr,reffrom%ptr,C_LOC(opt),uuout,vvout,uuin,vvin)
    end function georef_interpuv_f

!   c_ezwdint
    function georef_interpwd_f(this,reffrom,opt,uuout,vvout,uuin,vvin) result(out)
        class(georef),  intent(inout) :: this,reffrom  !< georef instance
        type(geooptions),  intent(in), target :: opt  !< georef instance
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin

        integer(C_INT32_T) :: out

        out=georef_interpwd(this%ptr,reffrom%ptr,C_LOC(opt),uuout,vvout,uuin,vvin)
    end function georef_interpwd_f

!   c_gduvfwd
    function georef_wd2uv_f(this,uuout,vvout,spdin,dirin,lat,lon,npts) result(out)
        class(georef),  intent(inout) :: this !< georef instance
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: spdin,dirin
        real(C_DOUBLE), intent(in), dimension(*) :: lat,lon
        integer(C_INT32_T) :: npts

        integer(C_INT32_T) :: out

        out=georef_wd2uv(this%ptr,uuout,vvout,spdin,dirin,lat,lon,npts)
    end function georef_wd2uv_f

!   c_gdwdfuv
    function georef_uv2wd_f(this,spdout,dirout,uuin,vvin,lat,lon,npts) result(out)
        class(georef),  intent(inout) :: this !< georef instance
        real(C_FLOAT), intent(out), dimension(*) :: spdout,dirout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin
        real(C_DOUBLE), intent(in), dimension(*) :: lat,lon
        integer(C_INT32_T) :: npts

        integer(C_INT32_T) :: out

        out=georef_uv2wd(this%ptr,spdout,dirout,uuin,vvin,lat,lon,npts)
    end function georef_uv2wd_f

!   c_gdlluvfuv_orig
    function georef_uv2uv_f(this,uuout,vvout,uuin,vvin,lat,lon,npts) result(out)
        class(georef),  intent(inout) :: this !< georef instance
        real(C_FLOAT), intent(out), dimension(*) :: uuout,vvout
        real(C_FLOAT), intent(in), dimension(*) :: uuin,vvin
        real(C_DOUBLE), intent(in), dimension(*) :: lat,lon
        integer(C_INT32_T) :: npts

        integer(C_INT32_T) :: out

        out=georef_uv2uv(this%ptr,uuout,vvout,uuin,vvin,lat,lon,npts)
    end function georef_uv2uv_f
end module georef_mod
