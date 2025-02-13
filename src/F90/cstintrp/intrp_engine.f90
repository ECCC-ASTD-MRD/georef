module intrp_engine
   USE mympp
   USE intrp_const
   USE grids
   
  !!-----------------------------------------------------------------------
  !!                 ***  MODULE intrp_oce  ***
  !!
  !!  **  Purpose: interpolation tools, F. Roy
  !!
  !!  **  Method: To have fun ...
  !!
  !! history :
  !!     Original code  :   F. Roy (March 2009-November 2009)
  !!-----------------------------------------------------------------------
  !!--------------------------------------------------------------
  !!

    IMPLICIT NONE
    ! * Arguments

    PUBLIC

    INTEGER werepeat, nrepeat, dxtree, dytree
    LOGICAL webdy
    CHARACTER(LEN=1) nbdy
    REAL(KIND=4), DIMENSION(:,:), POINTER :: longrd, latgrd
    INTEGER nxgrd, nygrd
    integer ipl, ipr, jpe ! start/end (l=left, r=right, e=end for the top) of the non-periodic domain indices

    CONTAINS

!---------------------------------------------------------------------------------------------
      SUBROUTINE interp_gen_2D(x1,y1,x2,y2,x3,y3,x4,y4,x,y,nitemax,a,b)
!---------------------------------------------------------------------------------------------

! General 2D bilinear interpolation for arbitrary 4 sides polygon
!
!                ------------X3
!      X4---------
!                |
!               |
!              |
!             |
!            |
!           G (point to interpolate to)
!          |
!    X1----------------X2
!
!
      implicit none
! Inputs
      real(kind=8), intent(in) :: x1,y1,x2,y2,x3,y3,x4,y4,x,y ! orientation 1-2-3-4 anticlockwise forms the quadrangle
      real(kind=8), intent(out) :: a, b
      integer, intent(in) :: nitemax ! maximum number of iteration to compute a and b
! Work
      real*8 m124,x1g4,x134,x12g,x123,den,d,e,ax,bx,cx,dx,u1,u2
      real*8, parameter :: maxres = 1e-12_8
      real*8 ap,bp, a0,b0, d0, d1
      integer :: ite, ite2
      real*8, dimension(:), allocatable  :: conv_histo

      allocate(conv_histo(nitemax))



!--------------------------------------------------
! project point and quadrangle onto local axes
!--- vectorial equation
! X1G = a(1-b) X1X2 + (1-a)b X1X4 + ab X1X3
!-apply the right cross product with x X1X4 and divide by X1X2 x X1X4
!--- leads to
! x' = a(1-b) + ab x3'
!--------------------------------------------------

      m124 = (y4-y1) * (x2-x1) - (y2-y1) * (x4-x1)
      x1g4 = (y4-y1) * (x -x1) - (y -y1) * (x4-x1)
      x134 = (y4-y1) * (x3-x1) - (y3-y1) * (x4-x1)
      if ( ABS(m124) < epsilon ) then
         a = -9.; b = -9.
         write(*,*) 'degenerate quadrangle. exiting'
         return
      endif
      x1g4 = x1g4 / m124
      x134 = x134 / m124

!--------------------------------------------------
!-apply the left cross product with X1X2 x and divide by X1X2 x X1X4
!--- leads to
! x" = (1-a)b + ab x3"
!--------------------------------------------------

      x12g = (y -y1) * (x2-x1) - (y2-y1) * (x -x1)
      x123 = (y3-y1) * (x2-x1) - (y2-y1) * (x3-x1)
      x12g = x12g / m124
      x123 = x123 / m124
     
!--------------------------------------------------
! create extra variables
!--------------------------------------------------

     d = x134 - oned; e = x123 - oned

!--------------------------------------------------
! general case
! solve iteratively the system of two non-linear equations using a Newton's method
! x' = a + (x3'-1) ab = a + d ab
! x" = b + (x3"-1) ab = b + e ab
!--------------------------------------------------
! after differentation in terms of a and b:
! x'- a - d * ab = da + d.b.da + d.a.db
! x"- b - e * ab = db + e.b.da + e.a.db
!--------------------------------------------------

      ! first guess
      a = max(min(x1g4, oned), zerd)
      b = max(min(x12g, oned), zerd)
      d1 = abs(a - x1g4 + d * a * b) + abs(b - x12g + e * a * b)
      ite = 1

      do while ( ite <= nitemax .and. d1 > maxres )
         ax = oned + d * b; bx = d * a; cx = e * b; dx = oned + e * a
         den = ax * dx - bx * cx
         u1 = d * a * b + a - x1g4; u2 = e * a * b + b - x12g
         a = a + ( - u1 * dx + u2 * bx ) / den
         b = b + ( - u2 * ax + u1 * cx ) / den
         d1 = abs(a - x1g4 + d * a * b) + abs(b - x12g + e * a * b)
         conv_histo(ite) = d1
         ite = ite + 1
      enddo
      if ( d1 > maxres ) then
            write(*,*) 'unsufficient convergence!', a, b, a+d*a*b-x1g4, b+e*a*b-x12g, x1g4, x12g, u1, u2, d1, ite
            write(*,*) 'x',x1,x2,x3,x4,x
            write(*,*) 'y',y1,y2,y3,y4,y
            write(*,*) 'convergence history'
            do ite2 = 1, ite - 1
               write(*,*) 'iter',ite2,conv_histo(ite2),maxres
            enddo
            call my_abort
      endif

!--------------------------------------------------
! add bound clipping in case of truncation error
! (should be very small with double precision
! and calculates weights for X1,2,3,4
!--------------------------------------------------
      a = max(min(a, oned), zerd)
      b = max(min(b, oned), zerd)

      deallocate(conv_histo)

      return
      END SUBROUTINE interp_gen_2D    


!---------------------------------------------------------------------------------------------
      SUBROUTINE ll_to_polars(r,latc,lonc,lat,lon,xp,yp)
!---------------------------------------------------------------------------------------------
      implicit none

! Lat-Lon to polar orthographic projection with
! center at latc,lonc

! Input/Outputs
      real*8 r,latc,lonc,lat,lon,xp,yp

! ... here latc and lonc and lat and lon are assumed to be degrees
! ... negative west of greenwich and south of equator

! ... latc,lonc are central coordinates
! ... r is the sphere radius
! ... k is the scale factor in both xp,yp directions

! Local variables
      real*8 lambda,lambda0,phi,phi1,k
      real(kind=8), parameter :: pi = atan(1.0_8) * 4.0_8, a=pi/180.0_8

      phi1=latc*a
      phi =lat*a
      lambda0=lonc*a
      lambda =lon*a

!     if stereographic:
!      k=2.d0/( 1.d0+sin(phi1)*sin(phi)    &
!     &             +cos(phi1)*cos(phi)*cos(lambda-lambda0) )

! the orthographic projection is actually easier because k=1
! the equations can be derived by first projecting to cartesian (x,y,z)
! and then applying a rotation that places at the pole the central location
      k = 1.0_8

      xp=r*k*cos(phi)*sin(lambda-lambda0)
      yp=r*k*( cos(phi1)*sin(phi)    &
     &        -sin(phi1)*cos(phi)*cos(lambda-lambda0) )

      return
      END SUBROUTINE ll_to_polars

!-------------------------------------------------------------------------------
      FUNCTION inpoly_complex(x,y,n,xx,yy)
!-------------------------------------------------------------------------------
! Determines if a point xx,yy is inside a complex polygon
! defined by x(n),y(n), the points must be correctly ordered anticloskwise
!-------------------------------------------------------------------------------
      implicit none
! arguments
      integer n
      real*8 x(n),y(n),xx,yy 
      logical inpoly_complex
! Local variables
      real*8, parameter :: epsplus = epsilon * 1e-5_8
      integer ncross,i1,i2
      real*8  dxa, dya, dxb, dyb, vecp, tota

      inpoly_complex=.false.

! Travel around segments
      ncross=0
      tota = zerd
      do i1=1,n        
        i2=mod(i1,n)+1
        dxa = x(i2)-x(i1)
        dya = y(i2)-y(i1)
        dxb = xx   -x(i1)
        dyb = yy   -y(i1)
        vecp = dyb*dxa - dya*dxb
        tota = tota + vecp
        if (vecp > - epsplus ) then
          ncross = ncross + 1
        endif
      enddo

! Now count nodes
      if ( tota > epsilon .and. ncross == n) then
        inpoly_complex=.true.
      endif
   
      END FUNCTION inpoly_complex

    !! --------------------------------------------------------------------------------------
    SUBROUTINE set_grid(grid)

    !! I/O Variables
    IMPLICIT NONE
    TYPE(type_grid), TARGET :: grid

    nxgrd     =  grid % nx
    nygrd     =  grid % ny
    longrd    => grid % lld % lont
    latgrd    => grid % lld % latt
    nbdy      =  grid % nbdy
    nrepeat   =  grid % nrepeat
    webdy     =  grid % webdy
    werepeat  =  grid % werepeat

    select case(werepeat)
    case(0) ! CICE case
      ipl = 1
      ipr = nxgrd
    case(1) ! GEM global grid
      ipl = 1
      ipr = nxgrd - 1
    case(2) ! ORCA grid
      ipl = 2
      ipr = nxgrd - 1
    end select
    jpe = nygrd - nrepeat

    END SUBROUTINE set_grid

    !! --------------------------------------------------------------------------------------
    SUBROUTINE refine_search(ilast, jlast, n_series, list_isrc, list_jsrc, lonw, latw)
      !! Modules

    IMPLICIT NONE
    ! arguments
    INTEGER ilast, jlast, n_series, list_isrc(n_series), list_jsrc(n_series)
    REAL(KIND=4) :: lonw, latw
    ! locals
    REAL(KIND=8) :: lonwd, latwd, x1, y1, x2, y2, xg, yg
    REAL(KIND=8) :: lonin, latin, m012, x0g2, x01g, a, b
    INTEGER inext, jnext, di, dj, d0, d, counter, ilast0, jlast0

    a = 1; b = 1
    inext = -nxgrd ; jnext = -nygrd
    ilast0 = inext; jlast0 = jnext
    d = abs(ilast - inext) + abs(jlast - jnext)
    d0 = d + 1
    counter = 0
   ! apply periodic bc on ref point

    call apply_periodic_condition_search(ilast, jlast, 0, 0, inext, jnext)

    ilast = inext; jlast = jnext
    

    do while ( d < d0 ) ! the next iteration should improve on the previous. if not bail out
      ! all double precision from here on
      ! projection done relative to point ilast, jlast => (0,0)
      lonwd = longrd(ilast,jlast); latwd = latgrd(ilast,jlast)

      ! find projected point +1,0 from ilast,jlast
      call apply_periodic_condition_search(ilast, jlast, 1, 0, inext, jnext, extend=.true.)
      lonin = longrd(inext,jnext); latin = latgrd(inext,jnext)
      call ll_to_polars(oned, latwd, lonwd, latin, lonin, x1, y1)

      ! find projected point 0,+1 from ilast,jlast
      call apply_periodic_condition_search(ilast, jlast, 0, 1, inext, jnext, extend=.true.)

      lonin = longrd(inext,jnext); latin = latgrd(inext,jnext)
      call ll_to_polars(oned, latwd, lonwd, latin, lonin, x2, y2)

      ! the projected target point relative to ilast, jlast
      lonin = lonw; latin = latw
      call ll_to_polars(oned, latwd, lonwd, latin, lonin, xg, yg)

      ! compute strides in x
!--------------------------------------------------
! project point onto local axes (x0,x1,x2)
!--- vectorial equation
! X0G = a X0X1 + b X0X2
!-apply the right cross product with x X0X2 and divide by X0X1 x X0X2
!--- leads to
! x' = a
!--------------------------------------------------

      m012 = y2 * x1 - y1 * x2
      if ( ABS(m012) < epsilon ) then
         n_series = 0
         exit ! likely a grid singularity, we do not want to interpolate close to it, exiting
      endif
      x0g2 = y2 * xg - yg * x2
      a = x0g2 / m012

      ! compute strides in y
!--------------------------------------------------
!-apply the left cross product with X0X1 x and divide by X0X1 x X0X2
!--- leads to
! x" = b
!--------------------------------------------------

      x01g = yg * x1 - y1 * xg
      b = x01g / m012
      di = nint(a); dj = nint(b)
      call apply_periodic_condition_search(ilast,jlast,di,dj,inext,jnext)

      ! preparation for next iteration
      d0 = d
      ilast0 = ilast
      jlast0 = jlast
      d = abs(di) + abs(dj)
      ilast = inext
      jlast = jnext

      counter = counter + 1
    enddo

    ! prepare list of closest cells
    !!! need to rethink this in terms of periodicity !!!
    ilast = ilast0; jlast = jlast0

    ! first search is closest point
    list_isrc(1) = ilast
    list_jsrc(1) = jlast

    ! second search is closest point to the right or left
    if (a<0.5_8) then 
      di = -1
    else
      di =  1
    endif
    call apply_periodic_condition_search(ilast,jlast,di,0,inext,jnext)
    list_isrc(2) = inext
    list_jsrc(2) = jnext

    ! third search is closest point on top or bottom
    if (b<0.5_8) then 
      dj = -1
    else
      dj =  1
    endif
    call apply_periodic_condition_search(ilast,jlast,0,dj,inext,jnext)
    list_isrc(3) = inext
    list_jsrc(3) = jnext

    ! fourth search is closest point in both directions
    call apply_periodic_condition_search(ilast,jlast,di,dj,inext,jnext)
    list_isrc(4) = inext
    list_jsrc(4) = jnext

    END SUBROUTINE refine_search

 
    SUBROUTINE apply_periodic_condition_search(ilast, jlast, di, dj, inext, jnext, extend)
    ! -----------------------------------------------------------------------------
    ! apply periodic conditions on global grids when searching for closest points
    ! -----------------------------------------------------------------------------
    ! arguments
    integer, intent(in) :: ilast, & ! i-index of the starting point
                           jlast, & ! j-index of the starting point
                           di, dj   ! jump stride in i/j
    logical, intent(in), optional :: extend
    integer, intent(out) :: inext, &
                            jnext ! final position after jump
    ! locals
    integer stri, strj, di2, inext2
    logical extendw

    extendw = .false.
    if (present(extend)) extendw = extend

    ! start with west-east
    inext = ilast + di ! apply stride
    if (webdy) then ! w-e periodicity active
      if (inext < ipl ) then
         stri = di - ipl + ilast + 1  ! remaining stride from right border
         inext = ipr + stri
      endif
      if (inext > ipr ) then
         stri = di - ipr + ilast - 1  ! remaining stride from left border
         inext = ipl + stri
      endif
    else ! simple case, no periodicity, force index to stay within grid
      if (inext >= nxgrd ) then
         if (extendw) then
              inext = nxgrd
         else
              inext = nxgrd - 1
         endif
      endif
      if (inext <      1 ) inext = 1
    endif

    ! y-axis north fold
    jnext = jlast + dj
    select case(nbdy)
    case('T') ! T-folding
      if (jnext > jpe) then
         strj = dj - jpe + jlast  ! remaining stride from top border
         jnext = jpe - strj
         inext = nxgrd - inext + 2
      endif
    case('F') ! F-folding
      if (jnext > jpe) then
         strj = dj - jpe + jlast - 1  ! remaining stride from top border
         jnext = jpe - strj
         inext = nxgrd - inext + 1
      endif
    case('N') ! Singular N-Pole
      if (jnext > jpe) then
         strj = dj - jpe + jlast - 1  ! remaining stride from top border
         jnext = jpe - strj
         di2 = (ipr - ipl + 1) / 2 ! shift of 180 degrees on the lat-lon grid
         inext2 = inext + di2
         if (inext2 > ipr ) then
            stri = di2 - ipr + inext - 1  ! remaining stride from left border
            inext = ipl + stri
         else
            inext = inext2
         endif
      endif
    case default
      if (jnext >= nygrd ) then
         if (extendw) then
              jnext = nygrd
         else
              jnext = nygrd - 1
         endif
      endif
    end select

    if (jnext <      1 ) jnext = 1 ! south boundary

    END SUBROUTINE apply_periodic_condition_search

end module intrp_engine
