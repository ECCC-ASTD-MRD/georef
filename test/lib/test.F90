program test
    use App
    use rmn_common
    use rmn_fst24
    use georef_mod
    implicit none

    type(fst_file)   :: file1, file2, fileout
    type(fst_record) :: record1, record2
    type(georef)     :: ref1, ref2, refout
    type(geoset)     :: set1, set2
    type(geodef)     :: def1, def2

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
        call App_log(APP_ERROR, 'Test failed: Unable to open FST file')
        call exit(-1)
    end if

    ! Read second file
    call get_command_argument(2,argument,len,status)
    success = file2%open(trim(argument))
    if (.not. success) then
        call App_log(APP_ERROR, 'Test failed: Unable to open FST file')
        call exit(-1)
    end if

    ! Create geo-references
    success=file1%read(record1,NOMVAR="GRID")
    success=ref1%fromrecord(record1)

    success=file2%read(record2,NOMVAR="GRID")
    success=ref2%fromrecord(record2)

    ! Check validity
    write(app_msg,*) 'ref1%valid = ',ref1%valid()
    call App_Log(APP_INFO, app_msg)

    ! Check if equal
    write(app_msg,*) 'ref1%equal = ',ref1%equal(ref2)
    call App_Log(APP_INFO, app_msg)

    ! Test copy
    refout=ref1%copy(hard=.TRUE.)
    write(app_msg,*) 'refout%valid = ',refout%valid()
    call App_Log(APP_INFO, app_msg)

    ! Test write to file
    success = fileout%open('TestF90.fst','R/W')
    if (.not. success) then
        call App_log(APP_ERROR, 'Test failed: Unable to open FST file')
        call exit(-1)
    end if

    success=ref1%writefst(fileout)
    if (.not. success) then
        call App_log(APP_ERROR, 'Test failed: Unable to write georef to FST file')
        call exit(-1)
    end if

    ! Test latlon limits
    success=ref1%limits(lat0,lon0,lat1,lon1)
    write(app_msg,*) 'ref1%limit = ',lat0,',',lon0,' - ',lat1,',',lon1
    call App_Log(APP_INFO, app_msg)

    ! Test ij coordinate of latlon bounding box
    lat0=36.0
    lon0=-116.0
    lat1=51.0
    lon1=-78.0
    success=ref1%boundingbox(lat0,lon0,lat1,lon1,x0,y0,x1,y1)
    write(app_msg,*) 'ref1%boundingbox = ',lat0,',',lon0,' - ',lat1,',',lon1,' => ',x0,',',y0,' - ',x1,',',y1
    call App_Log(APP_INFO, app_msg)

    ! Test intersection of 2 georeferences
    success=ref1%intersect(ref2,ix0,iy0,ix1,iy1)
    write(app_msg,*) 'ref1%intersect = ',ix0,',',iy0,' - ',ix1,',',iy1
    call App_Log(APP_INFO, app_msg)

    ! Test if on georeference is within another one
    success=ref1%within(ref2)
    write(app_msg,*) 'ref1%within (in=.FALSE.)= ',success
    call App_Log(APP_INFO, app_msg)

    ! Test if a georeference is within a latlon box
    success=ref1%withinrange(lat0,lon0,lat1,lon1,in=.TRUE.)
    write(app_msg,*) 'ref1%withinrange (in=.TRUE.) = ',success
    call App_Log(APP_INFO, app_msg)

    ! Test meter distance between grid points
    val=ref1%xydistance(x0,y0,x1,y1)
    write(app_msg,*) 'ref1%xydistance = ',val
    call App_Log(APP_INFO, app_msg)

    ! Test meter distancd between latlon points
    val=ref1%lldistance(lat0,lon0,lat1,lon1)
    write(app_msg,*) 'ref1%lldistance = ',val
    call App_Log(APP_INFO, app_msg)
    
    ! Test transformation from gridpoint to latlon
    x(1)=75.0
    y(1)=30.0
    len=ref1%xy2ll(lat,lon,x,y,1)
    write(app_msg,*) 'ref1%xy2ll = xy=',x(1),',',y(1),' => ll=',lat(1),',',lon(1)
    call App_Log(APP_INFO, app_msg)

    ! Test transformation from latlon to gridpoint
    lat(1)=31.5399963781238
    lon(1)=-127.250003948808
    len=ref1%ll2xy(x,y,lat,lon,1)
    write(app_msg,*) 'ref1%ll2xy = ll=',lat(1),',',lon(1),' => xy=',x(1),',',y(1)
    call App_Log(APP_INFO, app_msg)

    ! Test latlon stream extraction
    len=ref1%getll(lats,lons)
    write(app_msg,*) 'ref1%getll = ll=',len
    call App_Log(APP_INFO, app_msg)

    call record1%get_data_array(data_array1)
    call record2%get_data_array(data_array2)

    x(1)=75.5
    y(1)=30.5

    ! Test nearest grid point value interpolation
    georef_options%Interp=IR_NEAREST
    len=ref1%xyval(vals,data_array1,x,y,1,opt=georef_options)
    write(app_msg,*) 'ref1%xyval (nearest)= xy=',x(1),',',y(1),' => val=',vals(1)
    call App_Log(APP_INFO, app_msg)

    ! Test linear grid point value interpolation
    georef_options%Interp=IR_LINEAR
    len=ref1%xyval(vals,data_array1,x,y,1,opt=georef_options)
    write(app_msg,*) 'ref1%xyval (linear)= xy=',x(1),',',y(1),' => val=',vals(1)
    call App_Log(APP_INFO, app_msg)

    ! Test cubic grid point value interpolation
    georef_options%Interp=IR_CUBIC
    len=ref1%xyval(vals,data_array1,x,y,1,opt=georef_options)
    write(app_msg,*) 'ref1%xyval (cubic)= xy=',x(1),',',y(1),' => val=',vals(1)
    call App_Log(APP_INFO, app_msg)

    ! Test latlon value extraction
    len=ref1%llval(vals,data_array1,lat,lon,1)
    write(app_msg,*) 'ref1%llval = ll=',lat(1),',',lon(1),' => val=',vals(1)
    call App_Log(APP_INFO, app_msg)

    ! Test interpolation
    len=ref1%interp(ref2,data_array1,data_array2)
    write(app_msg,*) 'ref1%llval = ll=',lat(1),',',lon(1),' => val=',vals(1)
    call App_Log(APP_INFO, app_msg)
    record1%nomvar='data'
    success=fileout%write(record1)

    ! Test set
    set1=ref1%getset(ref2)
    success=set1%writefst(fileout)
    if (.not. success) then
        call App_log(APP_ERROR, 'Test failed: Unable to write geoset to FST file')
        call exit(-1)
    end if

    success=set1%readfst(IR_CUBIC,fileout)
    if (.not. success) then
        call App_log(APP_ERROR, 'Test failed: Unable to read geoset to FST file')
        call exit(-1)
    end if

    ! Test geodef
    data_array1=0.0
    data_array2=22.0
    success=def1%init(record1%ni,record1%nj,record1%nk,TD_Float32,c_loc(data_array1),C_NULL_PTR,C_NULL_PTR)
    success=def2%init(record2%ni,record2%nj,record2%nk,TD_Float32,c_loc(data_array2),C_NULL_PTR,C_NULL_PTR)
    if (.not. success) then
        call App_log(APP_ERROR, 'Test failed: Unable to create geodef')
        call exit(-1)
    end if


    georef_options%interp=IR_CONSERVATIVE
    len=ref1%interpdef(ref2,def1,def2,georef_options,.true.)
    record1%etiket='CONSERV'
    record1%data=c_loc(data_array1)
    success=fileout%write(record1)

    success = file1%close()
    success = file2%close()
    success = fileout%close()

    call App_Log(APP_INFO, 'Test successful');

end program test
