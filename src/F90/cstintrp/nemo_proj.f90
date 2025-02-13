MODULE nemo_proj

  !!-----------------------------------------------------------------------
  !!                 ***  MODULE nemo_proj ***
  !!
  !!  **  Purpose: 1) Read nemo coordinates and reperes
  !!               (getll_nemo_ext,getreps_nemo)
  !!               2) Convert ll to xy and xy to ll
  !!               (ll2xy_nemo,...)
  !!
  !!  **  Method: 1) Read lat-lons
  !!              2) Manage localisation into nemo grid
  !!              using local PS projections
  !!
  !! history :
  !!     Original code : F. Roy (October 2012)
  !!-----------------------------------------------------------------------
  !!--------------------------------------------------------------
  !!

  USE grids
  USE kdtree2_module
  USE mympp
  USE intrp_engine
  USE intrp_options, only: niter_max_weight_solver

  IMPLICIT NONE

  PUBLIC ll2xy_nemo

  Type (kdtree2), pointer :: my_kdtree
  REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: lonred, latred
  INTEGER nxred, nyred

  INTEGER, DIMENSION(:), POINTER :: iil, jjl

  CONTAINS

    !! --------------------------------------------------------------------------------------
    SUBROUTINE ll2xy_nemo(grid,lat,lon,AA,BB,ii,jj,M,nx,ny)

    !! Compute X,Y components from nemo grid
    !! Return values of -9. for masked points, M=0

    !! I/O Variables
    IMPLICIT NONE
    TYPE(type_grid), TARGET :: grid
    INTEGER , INTENT(in) :: nx,ny
    REAL(KIND=4), INTENT(in), DIMENSION (nx,ny) :: lat,lon
    REAL(KIND=8), INTENT(inout), DIMENSION (nx,ny) :: AA,BB
    INTEGER,      INTENT(inout), DIMENSION (nx,ny) :: ii,jj
    LOGICAL,      INTENT(inout), DIMENSION (nx,ny) :: M
    !! Local variables
    LOGICAL :: maskr
    INTEGER :: ie,je    !Grille entree
    REAL(KIND=4), PARAMETER :: r=1.0                    !Projection scale
    REAL(KIND=8) Xw,Yw
    INTEGER ntot, n_series, ipt, ii0, jj0, ilast, jlast, idx, ilast0, jlast0
    REAL(KIND=4) :: start, finish, my_xyz(3)
    INTEGER,      allocatable, dimension(:) :: list_isrc, list_jsrc
    TYPE(kdtree2_result), allocatable, dimension(:), target :: my_res
    INTEGER dxtree, dytree, nxred, nyred

    my_kdtree => grid % kdtree
    dxtree    =  grid % dxtree
    dytree    =  grid % dytree
    nxred     =  grid % nxred
    nyred     =  grid % nyred
    iil       => grid % ired
    jjl       => grid % jred

    call set_grid(grid)

    ! define and allocate search radius
    n_series = 1 ! just give the closest!
    allocate (my_res(n_series))
    n_series = 4
    allocate (list_isrc(n_series), &
              list_jsrc(n_series))

    write(*,*) 'my_kdtree',my_kdtree%n, my_kdtree%dimen
    call cpu_time(start)

    !! Boucle sur la grille d'entree
    DO je=1,ny
    DO ie=1,nx
      ii(ie,je) =-9
      jj(ie,je) =-9
      M (ie,je) = .false.
      AA(ie,je) = 0.
      BB(ie,je) = 0.
      !! ... find closest index from regular grid
      if (mod(ie-1,dxtree)==0) then
         call ez_lac(my_xyz, lon(ie,je), lat(ie,je), 1)
         n_series = 1
         call kdtree2_n_nearest(my_kdtree,my_xyz,n_series,my_res)
         ipt = 1
         idx = my_res(ipt) % idx
         jj0 = (idx-1) / nxred + 1
         ii0 = idx - nxred * (jj0 - 1)
         ilast = iil(ii0)
         jlast = jjl(jj0)
         ilast0 = ilast; jlast0 = jlast
      endif
      n_series = 4
      call refine_search(ilast, jlast, n_series, list_isrc, list_jsrc, lon(ie,je), lat(ie,je))
      if( n_series == 0 ) then
         ! reset the search indices to the last successful ones
         ilast = ilast0; jlast = jlast0
         cycle
      endif
      call calculate_XYM_dst_on_src(AA(ie, je), BB(ie, je), &
                         ii(ie, je), jj(ie, je), M(ie, je), lon(ie,je), lat(ie,je), &
                         n_series, list_isrc, list_jsrc, r)
      
      if (M(ie,je)) then
         ilast = ii(ie,je)
         jlast = jj(ie,je)
         ilast0 = ilast; jlast0 = jlast
      endif
    ENDDO
    ENDDO

    deallocate(my_res, list_isrc, list_jsrc)
    call cpu_time(finish)
    print '("ll2xy Time = ",f6.3," seconds.")',finish-start

    END SUBROUTINE ll2xy_nemo
    !! --------------------------------------------------------------------------------------

    !! --------------------------------------------------------------------------------------
    SUBROUTINE calculate_XYM_dst_on_src(AA, BB, ii, jj, msk, lonw, latw, npt, ipt, jpt, r)
      ! Object: calclualte X,Y and msk of the target grid with respect to the source grid

      implicit none

      ! output coordinates and mask of the target point wrt source grid
      real (kind=8), intent(inout) :: AA, BB
      integer,       intent(inout) :: ii,jj
      logical,       intent(inout) :: msk

      real (kind=4), intent(in)  :: lonw, latw ! coordinates of the target grid point
      !ranges of indices of the source gridpoints to be considered for the given target
      integer, intent (in) :: npt, ipt(npt), jpt(npt)
      REAL(KIND=4), intent (in) :: r ! projection scale
      ! locals
      integer :: k, inext, jnext
      real (kind=8) :: Xw, Yw, w1, w2, w3, w4, rd
      REAL(KIND=8), DIMENSION (4) :: lonin, latin, xpo, ypo  !nemo projection on PS
      logical :: inp ! =true if inpoly
      integer iip
      integer, dimension(4) ::  di = (/ 0, 1, 1, 0 /), dj = (/ 0, 0, 1, 1 /)

      rd = real(r,kind=8) ! conversion to double precision

      DO iip = 1, npt
        !! ... Project four corners on PS grid
        Xw = lonw; Yw = latw ! convert to double precision

        ii = ipt(iip)
        jj = jpt(iip)

        ! take the points anticlockwise and convert to double precision
        do k = 1, 4
            call apply_periodic_condition_search(ii, jj, di(k), dj(k), inext, jnext, extend = .true.)
            lonin(k) = longrd(inext,jnext)
            latin(k) = latgrd(inext,jnext)
            ! all double precision from here on
            call ll_to_polars(rd,Yw,Xw,latin(k),lonin(k),xpo(k),ypo(k))
        enddo

        Xw = zerd; Yw = zerd ! the projection is centered on the interpolated point
        !! ... Project entry grid lat-lon
        !! ... Check if point is inside nemo polygon
        inp=inpoly_complex(xpo,ypo,4,Xw,Yw)

        if (inp) then
          !! ... Project four corners on PS grid
          !! Interpolate x,y coordinates according to grid index ii,jj
          call interp_gen_2D(xpo(1),ypo(1),xpo(2),ypo(2),    &
                             xpo(3),ypo(3),xpo(4),ypo(4),    &
                             Xw,Yw, niter_max_weight_solver, AA,BB)
          msk = .true.
          return ! bail out, no need to continue, the first accceptable polygon is always the best!
        endif ! end condition on inside polygon
      ENDDO

    END SUBROUTINE calculate_XYM_dst_on_src

    !! --------------------------------------------------------------------------------------

END MODULE nemo_proj
