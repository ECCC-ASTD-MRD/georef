MODULE llxy

  !!-----------------------------------------------------------------------
  !!                 ***  MODULE llxy ***
  !!
  !!  **  Purpose: General grid conversion program lat-lon to x,y
  !!      the inverse, F. Roy
  !!
  !!  **  Method: Manage custom structured grids with masks and angles
  !!  **          assuming angle rotation at center of the grid
  !!
  !! history :
  !!     Original code  :   F. Roy (November 2009)
  !!-----------------------------------------------------------------------
  !!--------------------------------------------------------------
  !!

  !! Modules
  USE nemo_proj
  USE triangle_proj
  USE intrp_engine

  IMPLICIT NONE

  PUBLIC ll2xy

  CONTAINS

    SUBROUTINE ll2xy(grds,grdr,loc)

    !! Compute X,Y components and angle A from grid destination
    !! relative to source grid
    !! A is relative to X direction (i==>nx) 
    !! and positive anticlockwise (radians)
    !! Return values of -9. for masked points, M=0

    USE intrp_oce

    IMPLICIT NONE

    ! arguments
    TYPE(type_grid), TARGET, INTENT(inout) :: grds, grdr
    TYPE(type_location), TARGET, INTENT(inout) :: loc
    ! locals
    REAL(KIND=4), DIMENSION(:,:), POINTER :: lat,lon,A
    REAL(KIND=8), DIMENSION(:,:), POINTER :: AA,BB
    INTEGER,      DIMENSION(:,:), POINTER :: II,JJ
    LOGICAL,      DIMENSION(:,:), POINTER :: M
    INTEGER nx, ny, nxs, nys
    REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: x,y
    integer ier, gid
    REAL(KIND=8) :: w1,w2,w3,w4
    REAL(KIND=8) :: as1,as2,as3,as4,asm
    REAL(KIND=8) :: csm,snm,cs1,cs2,cs3,cs4,sn1,sn2,sn3,sn4
    integer i,j,is,js, i1, i2, i3, inext, jnext
    integer, external :: gdxyfll ! rmn functions
    REAL(Kind=4), dimension(3) :: di1, di2, di3, di4, dir, dis
    

    !-----------------------------------------------------------------------------------------
    ! find the fractional position on the source grid of the lat/lon of the destination grid
    !-----------------------------------------------------------------------------------------

    lon => grdr % lld % lont
    lat => grdr % lld % latt
    nx  =  grdr % nx
    ny  =  grdr % ny
    nxs =  grds % std_grd % ni
    nys =  grds % std_grd % nj
    gid =  grds % std_grd % grid_id

    call set_grid(grds)

    allocate( loc % ag(nx,ny) &
            , loc % bg(nx,ny) &
            , loc % ig(nx,ny) &
            , loc % jg(nx,ny) &
            , loc % pg(nx,ny) &
            , loc % mg(nx,ny) )
    if (grds % std_grd % grtyp .eq. 'Y') then
      write(*,*) 'CSTINTRP: This call of ll2xy should not have happened, gdefs="cloud"'
      call my_abort
    endif


    AA => loc % ag
    BB => loc % bg
    II => loc % ig
    JJ => loc % jg
    A  => loc % pg; A(:,:) = 0.
    M  => loc % mg; M(:,:) = .true.


    select case( grds % std_grd % grtyp )
    case('Z','B','U','A','N','S','L','G')
      allocate( x(nx,ny), y(nx,ny) )
      ier = gdxyfll(gid,x,y,lat,lon,nx*ny)
      II(:,:) = int(x(:,:))
      JJ(:,:) = int(y(:,:))
      AA(:,:) = x(:,:) - II(:,:)
      BB(:,:) = y(:,:) - JJ(:,:)
      deallocate(x,y)
      call process_north_pole(grds, grdr, loc)
    case('M')
      call construct_tree(grds)
      call ll2xy_triangle(grds,lat,lon,AA,BB,II,JJ,M,nx,ny)
    case default
      call construct_tree(grds)
      call ll2xy_nemo(grds,lat,lon,AA,BB,II,JJ,M,nx,ny)
    end select

    !-----------------------------------------------------------------------
    ! Compute the difference angle between the destination and source grids
    !-----------------------------------------------------------------------

    do j=1,ny
    do i=1,nx
      if (M(i,j)) then
        dir = grdr % lld % di3d(:,i,j)
        is = ii(i,j)
        js = jj(i,j)
        if ( grds % std_grd % grtyp == 'M') then
          w1 = AA(i,j)
          w2 = BB(i,j)
          w3 = MIN( MAX(zerd, oned - w1 - w2), oned)
          i1 = grds % meshd % in(1,is)
          i2 = grds % meshd % in(2,is)
          i3 = grds % meshd % in(3,is)
          di1 = grds % lld % di3d(:,i1,js)
          di2 = grds % lld % di3d(:,i2,js)
          di3 = grds % lld % di3d(:,i3,js)
          dis = w1 * di1 + w2 * di2 + w3 * di3 ! interpolated vector onto destination point
        else
          di1 = grds % lld % di3d( :, is   , js   )
          call apply_periodic_condition_search(is, js, 1, 0, inext, jnext)
          di2 = grds % lld % di3d( :, inext, jnext)
          call apply_periodic_condition_search(is, js, 1, 1, inext, jnext)
          di3 = grds % lld % di3d( :, inext, jnext)
          call apply_periodic_condition_search(is, js, 0, 1, inext, jnext)
          di4 = grds % lld % di3d( :, inext, jnext)
          w1 = (1._8-AA(i,j)) * ( 1._8-BB(i,j))
          w2 =       AA(i,j)  * ( 1._8-BB(i,j))
          w3 =       AA(i,j)  *        BB(i,j)
          w4 = (1._8-AA(i,j)) *        BB(i,j)
          dis = w1 * di1 + w2 * di2 + w3 * di3 + w4 * di4 ! interpolated vector onto destination point
        endif
        csm = dot_product( dis, dir )
        di1 = grdr % lld % t3d(:,i,j)
        snm = vector_pro( dis, dir, di1 ) ! the cross product is reduced to the normal component to the sphere
        asm = atan2(snm, csm)
        A(i,j) = asm
      endif
    enddo
    enddo

    END SUBROUTINE ll2xy

    SUBROUTINE process_north_pole( grds, grdr, loc )
    ! arguments
    TYPE(type_grid), TARGET, INTENT(inout) :: grds, grdr
    TYPE(type_location), TARGET, INTENT(inout) :: loc
    ! locals
    REAL(KIND=8), DIMENSION(:,:), POINTER :: AA,BB
    INTEGER,      DIMENSION(:,:), POINTER :: II,JJ
    LOGICAL,      DIMENSION(:,:), POINTER :: M
    INTEGER nxs, nys
    INTEGER i, j, inext, jnext
    REAL(KIND=8) x1, y1, x2, y2, x3, y3, a1, a2, d2r

    AA => loc % ag
    BB => loc % bg
    II => loc % ig
    JJ => loc % jg
    M  => loc % mg
    nxs =  grds % std_grd % ni
    nys =  grds % std_grd % nj

    d2r = atan(1.0_8) / 45.0_8  ! degree to radian conversion

    if ( grds % nbdy  == 'N' ) then
       do j=1,size(ii,2)
       do i=1,size(ii,1)
          if ( ii(i,j) > nxs ) then
             call apply_periodic_condition_search(ii(i,j), jj(i,j), 0, 0, inext, jnext, .true.)
             ii(i,j) = inext
             jj(i,j) = jnext
          endif
          if ( jj(i,j) == nys .and. BB(i,j) > zerd ) then ! recompute the linear position based on arc length
             call apply_periodic_condition_search(ii(i,j), jj(i,j), 0, 1, inext, jnext, .true.)
             x1 = grds % lld % lont(ii(i,j), jj(i,j)) * d2r; y1 = grds % lld % latt(ii(i,j), jj(i,j)) * d2r
             x2 = grds % lld % lont(inext  , jnext  ) * d2r; y2 = grds % lld % latt(inext  , jnext  ) * d2r
             x3 = grdr % lld % lont(i,j)              * d2r; y3 = grdr % lld % latt(i,j)              * d2r
             call arc( x1, y1, x2, y2, a1)
             call arc( x1, y1, x3, y3, a2)
             BB(i,j) = min( a2 / a1, oned)
          endif
       enddo
       enddo
    else
       where(ii < 1 .or. ii > nxs .or. jj < 1 .or. jj > nys) M = .false.
    endif

    END SUBROUTINE process_north_pole


END MODULE llxy
