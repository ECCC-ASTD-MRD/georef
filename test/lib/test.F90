program test
    use App
    use rmn_common
    use rmn_fst24
    use georef_mod
    implicit none

    type(fst_file)   :: file1,file2,fileout
    type(fst_record) :: record1,record2
    type(georef)     :: gref1,gref2,grefout

    logical :: success
    integer :: len, status
    character(len=4096) :: argument
    real(C_DOUBLE), dimension(5) :: lat,lon,x,y
    real(C_FLOAT), dimension(5) :: vals
    real(kind = real32), dimension(:), pointer :: data_array1, data_array2
    real(C_DOUBLE) :: lat0,lon0,lat1,lon1,x0,y0,x1,y1,val
    real(C_DOUBLE), dimension(2700000) :: lats,lons
    integer(C_INT32_T) :: ix0,ix1,iy0,iy1

    ! Read first file
    call get_command_argument(1,argument,len,status)
    success = file1%open(trim(argument))
    if (.not. success) then
        call App_log(APP_ERROR, 'Unable to open FST file')
        call exit(-1)
    end if

    ! Read second file
    call get_command_argument(2,argument,len,status)
    success = file2%open(trim(argument))
    if (.not. success) then
        call App_log(APP_ERROR, 'Unable to open FST file')
        call exit(-1)
    end if

    ! Create geo-references
    success=file1%read(record1,NOMVAR="GRID")
    success=gref1%fromrecord(record1)

    success=file2%read(record2,NOMVAR="GRID")
    success=gref2%fromrecord(record2)

    ! Check validity
    write(app_msg,*) 'gref1%valid = ',gref1%valid()
    call App_Log(APP_INFO, app_msg)

    ! Check if equal
    write(app_msg,*) 'gref1%equal = ',gref1%equal(gref2)
    call App_Log(APP_INFO, app_msg)

    ! Test copy
!TODO: compiler error
!    grefout=gref1%copy()

    ! Test write to file
    success = fileout%open('TestF90.fst','R/W')
    if (.not. success) then
        call App_log(APP_ERROR, 'Unable to open FST file')
        call exit(-1)
    end if
!TODO: No sure why this fails
    success=gref1%write(fileout)

    ! Test latlon limits
    success=gref1%limits(lat0,lon0,lat1,lon1)
    write(app_msg,*) 'gref1%limit = ',lat0,',',lon0,' - ',lat1,',',lon1
    call App_Log(APP_INFO, app_msg)
    
    ! Test ij coordinate of latlon bounding box
    lat0=36.0
    lon0=-116.0
    lat1=51.0
    lon1=-78.0
    success=gref1%boundingbox(lat0,lon0,lat1,lon1,x0,y0,x1,y1)
    write(app_msg,*) 'gref1%boundingbox = ',lat0,',',lon0,' - ',lat1,',',lon1,' => ',x0,',',y0,' - ',x1,',',y1
    call App_Log(APP_INFO, app_msg)

    ! Test intersection of 2 georeferences
    success=gref1%intersect(gref2,ix0,iy0,ix1,iy1)
    write(app_msg,*) 'gref1%intersect = ',ix0,',',iy0,' - ',ix1,',',iy1
    call App_Log(APP_INFO, app_msg)

    ! Test if on georeference is within another one
    success=gref1%within(gref2)
    write(app_msg,*) 'gref1%within (in=.FALSE.)= ',success
    call App_Log(APP_INFO, app_msg)

    ! Test if a georeference is within a latlon box
    success=gref1%withinrange(lat0,lon0,lat1,lon1,in=.TRUE.)
    write(app_msg,*) 'gref1%withinrange (in=.TRUE.) = ',success
    call App_Log(APP_INFO, app_msg)

    ! Test meter distance between grid points
    val=gref1%xydistance(x0,y0,x1,y1)
    write(app_msg,*) 'gref1%xydistance = ',val
    call App_Log(APP_INFO, app_msg)

    ! Test meter distancd between latlon points
    val=gref1%lldistance(lat0,lon0,lat1,lon1)
    write(app_msg,*) 'gref1%lldistance = ',val
    call App_Log(APP_INFO, app_msg)

    ! Test transformation from gridpoint to latlon
    x(1)=75.0
    y(1)=30.0
    len=gref1%xy2ll(lat,lon,x,y,1)
    write(app_msg,*) 'gref1%xy2ll = xy=',x(1),',',y(1),' => ll=',lat(1),',',lon(1)
    call App_Log(APP_INFO, app_msg)

    ! Test transformation from latlon to gridpoint
    lat(1)=31.5399963781238
    lon(1)=-127.250003948808
    len=gref1%ll2xy(x,y,lat,lon,1)
    write(app_msg,*) 'gref1%ll2xy = ll=',lat(1),',',lon(1),' => xy=',x(1),',',y(1)
    call App_Log(APP_INFO, app_msg)

    ! Test latlon stream extraction
    len=gref1%getll(lats,lons)
    write(app_msg,*) 'gref1%getll = ll=',len
    call App_Log(APP_INFO, app_msg)

    call record1%get_data_array(data_array1)
    call record2%get_data_array(data_array2)

    ! Test nearest grid point value interpolation
    georef_options%Interp=IR_NEAREST
    len=gref1%xyval(vals,data_array1,x,y,1,opt=georef_options)
    write(app_msg,*) 'gref1%xyval (nearest)= xy=',x(1),',',y(1),' => val=',vals(1)
    call App_Log(APP_INFO, app_msg)

    ! Test linear grid point value interpolation
    georef_options%Interp=IR_LINEAR
    len=gref1%xyval(vals,data_array1,x,y,1,opt=georef_options)
    write(app_msg,*) 'gref1%xyval (linear)= xy=',x(1),',',y(1),' => val=',vals(1)
    call App_Log(APP_INFO, app_msg)

    ! Test cubic grid point value interpolation
    georef_options%Interp=IR_CUBIC
    len=gref1%xyval(vals,data_array1,x,y,1,opt=georef_options)
    write(app_msg,*) 'gref1%xyval (cubic)= xy=',x(1),',',y(1),' => val=',vals(1)
    call App_Log(APP_INFO, app_msg)

    ! Test latlon value extraction
    len=gref1%llval(vals,data_array1,lat,lon,1)
    write(app_msg,*) 'gref1%llval = ll=',lat(1),',',lon(1),' => val=',vals(1)
    call App_Log(APP_INFO, app_msg)

    ! Test interpolation
    len=gref1%interp(gref2,data_array1,data_array2)
    write(app_msg,*) 'gref1%llval = ll=',lat(1),',',lon(1),' => val=',vals(1)
    call App_Log(APP_INFO, app_msg)
    success=fileout%write(record1)

    success = file1%close()
    success = file2%close()
    success = fileout%close()

    if (.not. success) then 
       call App_Log(APP_INFO, 'Test failed');
       call exit(-1)
    end if

    call App_Log(APP_INFO, 'Test successful');
    
end program test
