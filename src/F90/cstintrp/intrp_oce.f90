MODULE intrp_oce
   USE mympp
   use grids
   use intrp_engine

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

    PUBLIC mapospread, cintrp, fillm, mapospread_simple

    
    INTERFACE mapospread_simple ! user in vertical interpolation
       MODULE PROCEDURE mapospread_simple4, mapospread_simple8
    END INTERFACE

    CONTAINS

!---------------------------------------------------------------------------------------------
      SUBROUTINE mapospread( grd, spr, M, nspread, use_local_geo, nearst )
      implicit none
      type(type_grid),  target,  intent(in) :: grd
      type(type_weights), target, intent(inout) :: spr
      real(kind=8), dimension(:,:), intent(in) :: M
      integer, intent(in) :: nspread
      logical, intent(in) :: use_local_geo, nearst

! Spread data d around unmasked area
! If use_local_geo=.true. apply a weight average according to surrounding lat-lon's
! Assumes the mask is binary (0 or 1)
! Loop around nspread times

! nwgts should be set to 4*nspread, safest choice
! Theoretical max for island shape like this around 'x' for spread:
!    1    2      3        4
!  4*1   4*2    4*3      4x4
                          !
                 !       ! !
          !     ! !     !   !
     !   ! !   !   !   !     !
    !x! ! x ! !  x  ! !   x   !
     !   ! !   !   !   !     !
          !     ! !     !   !
                 !       ! !
                          !

! Local variables
      integer nwgts
      real(kind=4), dimension(:,:),   pointer :: lat,lon
      real(kind=8), dimension(:,:,:), pointer :: wgts
      integer,      dimension(:,:,:), pointer :: iwgts,jwgts
      real(kind=8) :: d(4),dtot12,dtot34,dtot   !Distances
      real(kind=8) :: mw
      real(kind=8) :: xp(4),yp(4)   !Projected lat-lon
      real(kind=8), parameter :: r = oned

      real(kind=8), allocatable, dimension(:,:) :: Ms,Mtmp  ! Spread work mask
      integer i,j,is,iw,iiw,nw,nwadd,iwgts_save,jwgts_save,inrs,jnrs
      real(kind=8) wgts_save
      logical first
      real(kind=8) lond,latd,lond1,latd1
      integer isp
      integer, dimension(4) :: ispp = (/-1, 1, 0, 0 /) ! i-array of adjacent sides
      integer, dimension(4) :: jspp = (/ 0, 0,-1, 1 /) ! j-array of adjacent sides
      real(kind=8), allocatable, dimension(:) :: wa
      integer, allocatable, dimension(:) :: iwa, jwa
      integer, dimension(:,:), pointer :: nws
      integer inrs2, jnrs2
      real(kind=8) wgts_save2
      integer nx,ny

      ! allocate the main arrays and some auxiliaries
      nx = grd % nx
      ny = grd % ny
      nwgts = spr % numwgt
      call set_grid(grd)

      allocate(Ms(nx,ny),Mtmp(nx,ny))

      write(*,*) 'CSTINTRP: ETALEMENT DES DONNEES SUR LES REGIONS MASQUEES nsprd:',nspread
      write(*,*) 'MAPOSPREAD: NUMBER OF WEIGHTS ALLOCATED FOR SPREADING=',nwgts
      nw = 4 * nspread * 4
      allocate(wa(nw),iwa(nw),jwa(nw))

      ! default definitions
      Ms=M
      Mtmp=M
       wgts => spr %  wgt
      iwgts => spr % iwgt
      jwgts => spr % jwgt
        nws => spr % nwgt
      lon => grd % lld % lont
      lat => grd % lld % latt

      wgts (1,:,:) = oned
      DO j=1,ny
      DO i=1,nx
        iwgts(1,i,j)=i
        jwgts(1,i,j)=j
        nws(i,j) = 1
        if (Ms(i,j).lt.0.5_8) then
          wgts(1,i,j) = zerd
          nws(i,j) = 0
        endif
      ENDDO
      ENDDO
      !
      ! start spreading iteratively
      !
      DO is=1,nspread
        !... defines the weights for the closest points around i,j using inverse distance
        DO j=1,ny
        DO i=1,nx
          mw = zerd
          do isp = 1,4
            call apply_periodic_condition_search(i, j, ispp(isp), jspp(isp), inrs, jnrs)
            mw = mw + Ms(inrs, jnrs)
          enddo
          if (Ms(i,j).lt.0.5_8.and.mw.gt.0.5_8) then
            if (nearst) then ! why nearst would use a different spreading technique? Moreover the technique is biased towards the last wet point found
              do isp = 1,4
                call apply_periodic_condition_search(i, j, ispp(isp), jspp(isp), inrs, jnrs)
                if (Ms(inrs,jnrs).gt.0.5_8) then
                  wgts (1,i,j) = 1.
                  iwgts(1,i,j) = inrs
                  jwgts(1,i,j) = jnrs
                  nws  (  i,j) = 1
                endif
              enddo
            else
              if (use_local_geo) then
                ! conversion to double precision
                lond = lon(i,j) ; latd = lat(i,j)
                do isp = 1,4
                  call apply_periodic_condition_search(i, j, ispp(isp), jspp(isp), inrs, jnrs)
                ! ... the center lat-lon is at x,y=0,0
                  lond1 = lon(inrs,jnrs); latd1 = lat(inrs,jnrs)
                  call ll_to_polars(r,latd,lond,latd1,lond1,xp(isp),yp(isp))
                  d(isp) = sqrt((xp(isp)**2+yp(isp)**2))
                  d(isp) = max( d(isp), 1e-7_8 ) ! because some points are actually co-located on the ORCA grid
                  d(isp) = r / d(isp) * Ms(inrs,jnrs) ! use the inverse distance metric for weight => the closer, the higher the weight
                enddo
              else
                do isp = 1,4
                  call apply_periodic_condition_search(i, j, ispp(isp), jspp(isp), inrs, jnrs)
                  d(isp) = Ms(inrs,jnrs)
                enddo
              endif
              dtot = SUM( d )
              nw = 0
              do isp = 1,4
                if (d(isp) > 0.0_8) then
                  nw = nw + 1
                  call apply_periodic_condition_search(i, j, ispp(isp), jspp(isp), inrs, jnrs)
                   wgts(nw,i,j) = d(isp) / dtot
                  iwgts(nw,i,j) = inrs
                  jwgts(nw,i,j) = jnrs
                endif
              enddo
              nws(i,j) = nw
            endif
            Mtmp(i,j)=1.
          endif
        ENDDO
        ENDDO
! ... Fix weights if is > 1
! ... Goal here is to trace back points that were precedingly used for spread
! ... so that all weights correspond to wet points
        if (is.gt.1) then
          DO j=1,ny
          DO i=1,nx
            mw = zerd
            do isp = 1,4
              call apply_periodic_condition_search(i, j, ispp(isp), jspp(isp), inrs, jnrs)
              mw = mw + Ms(inrs, jnrs)
            enddo
            if (Ms(i,j).lt.0.5_8.and.mw.gt.0.5_8) then
              nw = nws(i,j)
              iwa(1:nw) = iwgts(1:nw,i,j); iwgts(1:nw,i,j) = 1
              jwa(1:nw) = jwgts(1:nw,i,j); jwgts(1:nw,i,j) = 1
               wa(1:nw) =  wgts(1:nw,i,j);  wgts(1:nw,i,j) = 0.0_8
              DO iw = 1, nws(i,j)
                inrs = iwa(iw)
                jnrs = jwa(iw)
                if (M(inrs,jnrs).lt.0.5_8) then !Dry point weight, need to find a substitute
                  wgts_save = wa(iw)
                  iwa(iw) = 1
                  jwa(iw) = 1
                   wa(iw) = 0.0_8
                  DO iiw = 1, nws(inrs,jnrs)
                    inrs2 = iwgts(iiw,inrs,jnrs)
                    jnrs2 = jwgts(iiw,inrs,jnrs)
                    wgts_save2 = wgts(iiw,inrs,jnrs)
                    nw = nw + 1
                    iwa(nw) = inrs2
                    jwa(nw) = jnrs2
                     wa(nw) = wgts_save2 * wgts_save
                  ENDDO
                endif
              ENDDO
              ! ... eliminate redundancy
              DO iw= 1, nw-1
                inrs = iwa(iw)
                jnrs = jwa(iw)
                DO iiw = iw+1, nw
                  ! find the point if already stored in wa
                  inrs2 = iwa(iiw)
                  jnrs2 = jwa(iiw)
                  if (inrs2 == inrs .and. jnrs2 == jnrs .and. wa(iiw) > 0.0_8) then
                        wa(iw) = wa(iw) + wa(iiw)
                        wa(iiw) = 0.0_8
                        iwa(iiw) = 1
                        jwa(iiw) = 1
                   endif
                ENDDO
              ENDDO
              ! ... compaction of the resulting array
              iw = 1
              DO WHILE (iw <= nw )
                if (wa(iw) == 0.0_8) then
                  iwa(iw:nw-1) = iwa(iw+1:nw)
                  jwa(iw:nw-1) = jwa(iw+1:nw)
                   wa(iw:nw-1) =  wa(iw+1:nw)
                  nw = nw - 1
                else
                  iw = iw + 1
                endif
              ENDDO
              if (nw.gt.nwgts) then
                write(*,*) 'i,j,nw=',i,j,nw,nwgts
                write(*,*) 'nw too large'
                call my_abort
              endif
              ! copy back to main array
              iwgts(1:nw,i,j) = iwa(1:nw)
              jwgts(1:nw,i,j) = jwa(1:nw)
               wgts(1:nw,i,j) =  wa(1:nw)
               nws(i,j) = nw
            endif
          ENDDO
          ENDDO
        endif
        Ms=Mtmp
      ENDDO   !is=1,nspread

      ! reset nwgts based on the compacted array
      nwgts = maxval(nws)
      spr % numwgt = nwgts

      deallocate(Ms,Mtmp)
      deallocate(iwa,jwa,wa)
      write(*,*) 'CSTINTRP: MAX NUMBER OF EXTRA WEIGHTS AFTER SPREADING = ', nwgts

      END SUBROUTINE mapospread

      SUBROUTINE cintrp( loc, intr, bicub, nxsrc, nysrc )
      implicit none
      ! arguments
      type(type_location), intent(in), target :: loc
      type(type_weights), intent(inout), target :: intr
      logical,         intent(in) :: bicub
      integer,         intent(in) :: nxsrc, nysrc
      ! locals
      integer nxdst,nydst
      real(kind=8),    dimension(:,:,:), pointer :: wgt
      integer(kind=4), dimension(:,:,:), pointer :: iwgt,jwgt
      integer,         dimension(:,:)  , pointer :: II,JJ
      real(kind=8),    dimension(:,:)  , pointer :: AA, BB
      logical,         dimension(:,:)  , pointer :: ML
      integer,         dimension(:,:)  , pointer :: nwgt

! knowing that Idst,Jdst correspond to localisation indexes 
! (linking source grid to destination grid)
! If defrm=.true. consider lat-lon grid deformation
! by projecting on a locally centered PS grid
! Exclude computation over masked data (Mdst).
! Bi-Linear, default if bicub=.false.
! Bi-cubique, if bicub=.true.

      integer i,j,k,is,js,inext,jnext
      real(kind=8) :: x1,x2,x3,x4,y1,y2,y3,y4,x,y  !Projected lat-lon
      real(kind=8), parameter :: r=1.0_8
      real(kind=8) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,d4
      real(kind=8) :: one=1.,two=2.,ov6=1./6.
      real(kind=8) lond,latd,lond1,latd1
     
      II => loc % ig; JJ => loc % jg; ML => loc % mg
      AA => loc % ag; BB => loc % bg
       wgt => intr %  wgt
      iwgt => intr % iwgt
      jwgt => intr % jwgt
      nwgt => intr % nwgt
      nxdst = intr % nx
      nydst = intr % ny

      DO j=1,nydst
      DO i=1,nxdst
        if( ML(i,j))  then
          is=II(i,j)
          js=JJ(i,j)
          if (bicub) then
            c3 = AA(i,j) ! a, or the fractional position in i
            d3 = BB(i,j) ! b, or the fractional position in j
            c1 = one - c3
            c0 = - ov6 * c1 * c3
            c2 = c0 * ( two - c3 )
            c4 = c0 * ( one + c3 )
           
            d1 = one - d3
            d0 = - ov6 * d1 * d3
            d2 = d0 * ( two - d3  )
            d4 = d0 * ( one + d3 )
           
            wgt ( 1,i,j) = d1*(c1-2.*c2+c4) - 2.*d2*c1 + d4*c1 ! f_i,j
            iwgt( 1,i,j) = is
            jwgt( 1,i,j) = js
            call apply_periodic_condition_search(is, js, 0, 1, inext, jnext, extend = .true.)
            wgt ( 2,i,j) = d3*(c1-2.*c2+c4) - 2.*d4*c1 + d2*c1 ! f_i,j+1
            iwgt( 2,i,j) = inext
            jwgt( 2,i,j) = jnext 
            call apply_periodic_condition_search(is, js, 1, 1, inext, jnext, extend = .true.)
            wgt ( 3,i,j) = d3*(c3+c2-2.*c4) + d2*c3 - 2.*d4*c3 ! f_i+1,j+1
            iwgt( 3,i,j) = inext
            jwgt( 3,i,j) = jnext
            call apply_periodic_condition_search(is, js, 1, 0, inext, jnext, extend = .true.)
            wgt ( 4,i,j) = d1*(c3+c2-2.*c4) - 2.*d2*c3 + d4*c3 ! f_i+1,j
            iwgt( 4,i,j) = inext
            jwgt( 4,i,j) = jnext
            call apply_periodic_condition_search(is, js, -1, 0, inext, jnext, extend = .true.)
            wgt ( 5,i,j) = d1*c2 ! f_i-1,j
            iwgt( 5,i,j) = inext
            jwgt( 5,i,j) = jnext
            call apply_periodic_condition_search(is, js, 0, -1, inext, jnext, extend = .true.)
            wgt ( 6,i,j) = d2*c1 !f_i,j-1
            iwgt( 6,i,j) = inext
            jwgt( 6,i,j) = jnext
            call apply_periodic_condition_search(is, js, 1, -1, inext, jnext, extend = .true.)
            wgt ( 7,i,j) = c3*d2 ! f_i+1,j-1
            iwgt( 7,i,j) = inext
            jwgt( 7,i,j) = jnext
            call apply_periodic_condition_search(is, js, -1, 1, inext, jnext, extend = .true.)
            wgt ( 8,i,j) = d3*c2 ! f_i-1,j+1
            iwgt( 8,i,j) = inext
            jwgt( 8,i,j) = jnext
            call apply_periodic_condition_search(is, js, 2, 0, inext, jnext, extend = .true.)
            wgt ( 9,i,j) = d1*c4 ! f_i+2,j
            iwgt( 9,i,j) = inext
            jwgt( 9,i,j) = jnext
            call apply_periodic_condition_search(is, js, 2, 1, inext, jnext, extend = .true.)
            wgt (10,i,j) = d3*c4 ! f_i+2,j+1
            iwgt(10,i,j) = inext
            jwgt(10,i,j) = jnext 
            call apply_periodic_condition_search(is, js, 0, 2, inext, jnext, extend = .true.)
            wgt (11,i,j) = d4*c1 ! f_i,j+2
            iwgt(11,i,j) = inext
            jwgt(11,i,j) = jnext
            call apply_periodic_condition_search(is, js, 1, 2, inext, jnext, extend = .true.)
            wgt (12,i,j) = d4*c3 ! f_i+1,j+2
            iwgt(12,i,j) = inext
            jwgt(12,i,j) = jnext
            nwgt(   i,j) = 12
            ! WARNING
          else
            !linear interpolation default (bicub=.false.)
            wgt (1,i,j) = (1._8-AA(i,j)) * (1._8-BB(i,j))
            iwgt(1,i,j) = is
            jwgt(1,i,j) = js
            call apply_periodic_condition_search(is, js, 1, 0, inext, jnext, extend = .true.)
            wgt (2,i,j) =       AA(i,j)  * (1._8-BB(i,j))
            iwgt(2,i,j) = inext
            jwgt(2,i,j) = jnext
            call apply_periodic_condition_search(is, js, 1, 1, inext, jnext, extend = .true.)
            wgt (3,i,j) =       AA(i,j)  *       BB(i,j)
            iwgt(3,i,j) = inext
            jwgt(3,i,j) = jnext
            call apply_periodic_condition_search(is, js, 0, 1, inext, jnext, extend = .true.)
            wgt (4,i,j) = (1._8-AA(i,j)) *       BB(i,j)
            iwgt(4,i,j) = inext
            jwgt(4,i,j) = jnext
            nwgt(  i,j) = 4
            ! interpolate
          endif

        endif
      ENDDO
      ENDDO

      END SUBROUTINE cintrp

      SUBROUTINE caggrt( coli, grds, grdr, intr, defrm, nadis, Ms)
      implicit none
      ! arguments
      type(type_location), intent(in), target :: coli ! reverse location
      type(type_grid)    , intent(in), target :: grds, grdr
      type(type_weights), intent(inout), target :: intr
      real(kind=8),dimension(:,:), intent(in) :: Ms
      logical, intent(in) :: defrm
      integer, intent(in) :: nadis
      ! locals
      real(kind=8),    dimension(:,:,:), pointer :: wgt
      integer(kind=4), dimension(:,:,:), pointer :: iwgt,jwgt
      integer,         dimension(:,:)  , pointer :: II,JJ
      real(kind=8),    dimension(:,:)  , pointer :: AA, BB
      logical,         dimension(:,:)  , pointer :: ML, Mdsta
      integer, dimension(:,:), pointer :: Nadst

! Return mask where no data on destination grid (Mdsta)
! No computation on source mask Msrc

      integer :: naggrmax, naggrtot
      integer i, j, istep, jstep, k, is, js, icl, jcl, icldm, jcldm
      logical, dimension(:,:)  , allocatable :: Msrc, Mtmp
      real(kind=8), allocatable, dimension(:,:)   :: surtot 
      real(kind=8) :: surwrk
      real(kind=8), dimension(:,:), allocatable :: dxs, dys
      integer nmixt, naggr, nwarn, nwgtw
      logical :: print_warning, printw
      integer ntmp
      integer, dimension(:), allocatable :: itmp, jtmp
      integer, dimension(4) :: isi = (/ -1, 1, 1, -1 /)
      integer, dimension(4) :: jsi = (/ -1,-1, 1,  1 /)
      integer, dimension(4) :: ile = (/  1, 0,-1,  0 /)
      integer, dimension(4) :: jle = (/  0, 1, 0, -1 /)
      integer irad, isid, ilen, nsid, nlen

      II => coli % ig; JJ => coli % jg; ML => coli % mg
      AA => coli % ag; BB => coli % bg
       wgt => intr %  wgt; wgt=0.
      iwgt => intr % iwgt; iwgt=1  !To avoid index overflow later
      jwgt => intr % jwgt; jwgt=1
      Mdsta=> intr % mwgt; Mdsta=.false.
      Nadst=> intr % nwgt; Nadst=0
      naggrmax = intr % numwgt

      icldm = -1
      jcldm = -1


      if (defrm) then
         allocate( dxs( grds % nx, grds % ny), dys( grds % nx, grds % ny) )
         call get_cell_dxdy( grds, dxs, dys )
      endif

      ! Action

      allocate(surtot(grdr % nx, grdr % ny)); surtot=0.
      allocate(  Msrc(grds % nx, grds % ny) ); Msrc  = .false.; where( Ms > 0.5_8) Msrc = .true.
      allocate(  Mtmp(grdr % nx, grdr % ny) ); Mtmp(:,:) = .true.

      call remove_periodic_points(grds, Msrc)
      where( .not.ML ) Msrc = .false. ! where the destination grid does not match the source grid

      ! ... First compile closest point in 3D table
      naggrtot=0
      nwgtw=max(nadis-1,0)
      print_warning=.false.
      nwarn=0
      allocate(itmp((2*nwgtw+1)**2))
      allocate(jtmp((2*nwgtw+1)**2))

! Compute weighting surfaces
! from closest point in defined neighbour
      do js = 1, grds % ny
      do is = 1, grds % nx
        if ( Msrc(is,js) ) then
         icl = II(is,js) + nint(AA(is,js))
         jcl = JJ(is,js) + nint(BB(is,js))
         ntmp = 0
         do irad = 0, nwgtw
         nsid = 4
         nlen = 2*irad-1
         if (irad == 0) then
             nsid = 1
             nlen = 0
         endif
         do isid = 1, nsid
         do ilen = 0, nlen
           istep = isi(isid) * irad + ile(isid) * ilen
           jstep = jsi(isid) * irad + jle(isid) * ilen
           call apply_periodic_condition_search(icl, jcl, istep, jstep, i, j, extend = .true.)
           if (Mtmp(i,j)) then ! make sure to go over the following point only once
              Mtmp(i,j) = .false.
              ! follow the points used on the destination grid
              ntmp = ntmp + 1; itmp (ntmp) = i; jtmp(ntmp) = j
              ! Compute weight surface for each point neglecting surface overflows from destination grid point
              if (defrm) then
                 surwrk = dxs(is,js) * dys(is,js)
              else
                 surwrk = oned  !Equal weighting (simple average)
              endif
              k = Nadst(i,j) + 1
              if (k <= naggrmax) then
                Nadst(i,j) = k
                if (k > naggrtot) then
                  icldm = i
                  jcldm = j
                  naggrtot = k
                endif
                surtot(i,j)=surtot(i,j)+surwrk
                wgt (k,i,j)=surwrk
                iwgt(k,i,j)=is
                jwgt(k,i,j)=js
              else
                print_warning=.true. 
                nwarn=nwarn+1
              endif
           endif
         enddo
         enddo
         enddo
         ! reset to true all points used on the destinatin grid
         do k = 1, ntmp
            Mtmp(itmp(k),jtmp(k)) = .true.
         enddo
        endif
      enddo
      enddo
      if (print_warning) then
        write(*,*) 'CSTINTRP: WARNING not enough storage to aggregate all points'
        write(*,*) 'CSTINTRP: FOR THIS NUMBER OF DESTINATION POINTS ==>', nwarn
        write(*,*) 'CSTINTRP: naggrmax=',naggrmax
      endif
   
      if (defrm) then
         deallocate( dxs, dys )
      endif

! Now aggregate over destination grid
      naggr=0
      do j = 1, grdr % ny
      do i = 1, grdr % nx
        if (Nadst(i,j).gt.0) then
          naggr = naggr + 1
          Mdsta(i,j) = .true.
          if (surtot(i,j).gt.0.) then
            do k=1,Nadst(i,j)
              wgt(k,i,j)=wgt(k,i,j)/surtot(i,j)
            enddo
          else !case where all weights are zero (equal weights, no surfaces)
            write(*,*) 'SOMETHING WRONG WITH WEIGHTING AT DESTINATION',i,j
            write(*,*) 'WEIGHTS WERE ALL FIXED TO ZERO, ... HOPING IT IS A DRY POINT'
          endif
        endif
      enddo
      enddo

      write(*,*) 'CSTINTRP: Local max number of points in aggregation, naggrtot=',naggrtot
      write(*,*) 'CSTINTRP: Local max number of points at i,j=',icldm,jcldm
      write(*,*) 'CSTINTRP: Total number of aggregated points=',naggr

      if (.not.defrm) then
        write(*,*) 'CSTINTRP: Local lat-lon deformation neglected'
      endif

      deallocate(surtot, Msrc, Mtmp)
      intr % numwgt = naggrtot

      END SUBROUTINE caggrt


  SUBROUTINE adjust_wgtb( pintr, psprd, Ms )

    implicit none
    ! arguments
    type(type_weights), intent(inout), target :: pintr
    type(type_weights), intent(in)   , target :: psprd
    real(kind=8), intent(in)     :: Ms(:,:)
    ! locals
    integer nxdst,nydst, nxsrc, nysrc
    integer nwgtb, nwgts, nwgtbmax
    real(kind=8), dimension(:,:,:), pointer  :: wgtb, wgts
    integer,      dimension(:,:,:), pointer  :: iwgtb, jwgtb, iwgts, jwgts
    integer, dimension(:,:), pointer :: nws
    integer i,j,nwbase, nw, iw, iiw, nnw
    integer inrs, jnrs, inrs2, jnrs2
    real(kind=8) wgt_save, wgt_save2
    integer ib,jb
    real(kind=8), allocatable, dimension(:) :: wa
    integer, allocatable, dimension(:) :: iwa, jwa

    write(*,*) 'CSTINTRP: ADJUSTING WEIGHTS ACCORDING TO SPREAD'

     wgtb => pintr %  wgt
    iwgtb => pintr % iwgt
    jwgtb => pintr % jwgt
      nws => pintr % nwgt
     wgts => psprd %  wgt
    iwgts => psprd % iwgt
    jwgts => psprd % jwgt

    nwbase = SIZE( pintr % wgt, 1 )
    nwgtb  = pintr % numwgt
    nwgtbmax=0
    nxdst = pintr % nx
    nydst = pintr % ny
    nxsrc = SIZE( Ms, 1)
    nysrc = SIZE( Ms, 2)
    nw = 4 * nwbase
    allocate(wa(nw),iwa(nw),jwa(nw))

    DO j=1,nydst
    DO i=1,nxdst
      nw = nwgtb
      iwa(1:nw) = iwgtb(1:nw,i,j); iwgtb(1:nw,i,j) = 1
      jwa(1:nw) = jwgtb(1:nw,i,j); jwgtb(1:nw,i,j) = 1
       wa(1:nw) =  wgtb(1:nw,i,j);  wgtb(1:nw,i,j) = 0.0_8
      DO iw = 1, pintr % numwgt
        inrs = iwa(iw)
        jnrs = jwa(iw)
        nnw = psprd % nwgt(inrs,jnrs)
        if (wa(iw).ne.0.0_8.and.Ms(inrs,jnrs).lt.0.5_8.and.nnw>0) then !Dry point weight, need to find a substitute
          wgt_save = wa(iw)
          iwa(iw) = 1
          jwa(iw) = 1
           wa(iw) = zerd
          DO iiw =1 , nnw
            inrs2 = iwgts(iiw,inrs,jnrs)
            jnrs2 = jwgts(iiw,inrs,jnrs)
            wgt_save2 = wgts(iiw,inrs,jnrs)
            nw = nw + 1
            iwa(nw) = inrs2
            jwa(nw) = jnrs2
             wa(nw) = wgt_save2 * wgt_save
          ENDDO
        else if (wa(iw).ne.0.0_8.and.Ms(inrs,jnrs).lt.0.5_8.and.nnw==0) then !Dry point weight, no substitute available
          wa(iw) = zerd
        endif
      ENDDO
      ! ... eliminate redundancy
      DO iw = 1, nw-1
        inrs = iwa(iw)
        jnrs = jwa(iw)
        DO iiw = iw+1, nw
          ! find the point if already stored in wa
          inrs2 = iwa(iiw)
          jnrs2 = jwa(iiw)
          if (inrs2 == inrs .and. jnrs2 == jnrs .and. wa(iiw) .ne. 0.0_8) then
                wa(iw) = wa(iw) + wa(iiw)
                wa(iiw) = zerd
                iwa(iiw) = 1
                jwa(iiw) = 1
           endif
        ENDDO
      ENDDO
      ! ... compaction of the resulting array
      iw = 1
      DO WHILE (iw <= nw )
        if (wa(iw) == zerd) then
          iwa(iw:nw-1) = iwa(iw+1:nw)
          jwa(iw:nw-1) = jwa(iw+1:nw)
           wa(iw:nw-1) =  wa(iw+1:nw)
          nw = nw - 1
        else
          iw = iw + 1
        endif
      ENDDO
      if (nw.gt.nwbase) then
        write(*,*) 'i,j,nw=',i,j,nw,nwbase
        write(*,*) 'nw too large'
        call my_abort
      endif
      ! copy back to main array
      iwgtb(1:nw,i,j) = iwa(1:nw)
      jwgtb(1:nw,i,j) = jwa(1:nw)
       wgtb(1:nw,i,j) =  wa(1:nw)
       nws(i,j) = nw
      nwgtbmax = max(nwgtbmax,nw)
    ENDDO
    ENDDO

    ! final verification that sum of weights = 1
    DO j=1,nydst
    DO i=1,nxdst
       wgt_save  = SUM( wgtb(1:nws(i,j),i,j) )
       wgtb(1:nws(i,j),i,j) = wgtb(1:nws(i,j),i,j) / wgt_save
    ENDDO
    ENDDO

    write(*,*) 'CSTINTRP: MAX NUMBER OF WEIGHTS AFTER SPREAD = ', nwgtbmax
    pintr % numwgt = nwgtbmax
    deallocate(wa,iwa,jwa)

  END SUBROUTINE adjust_wgtb


      SUBROUTINE rotate_vector_ninj(u,v,theta,ni,nj)
!--------------------------------------------------------------------------
! Sous-Routine rotate_vector
! Tourne un vecteur (u,v) dans une base orthonormee ayant un angle
! de theta (dans le sens horaire) par rapport a l'ancienne base
!--------------------------------------------------------------------------

! Dans le modele ROM theta=49 degres (il faut transformer en radians)
! si on va de (x,y) vers (est,nord) theta= 49 deg
! si on va de (est,nord) vers (x,y) theta=-49 deg

      implicit none
      integer ni,nj
      real*8 u(ni,nj),v(ni,nj)
      real theta(ni,nj)
      real, allocatable :: utmp(:,:),vtmp(:,:)

      allocate(utmp(ni,nj),vtmp(ni,nj))
      utmp(:,:) = cos(theta(:,:))*u(:,:) - sin(theta(:,:))*v(:,:)
      vtmp(:,:) = sin(theta(:,:))*u(:,:) + cos(theta(:,:))*v(:,:)

      u = utmp
      v = vtmp

      deallocate(utmp,vtmp)
      return
      END SUBROUTINE rotate_vector_ninj

      SUBROUTINE fillm(var,mskd,mskr,nx,ny,npts,npass)
      ! Fill holes in a field using mask (msk)
      ! based on a pseudo-inhouse kriging method (linear)
      implicit none
      integer, intent(in)                      :: nx,ny
      real, intent(inout),    dimension(nx,ny) :: var
      integer, intent(inout), dimension(nx,ny) :: mskd
      integer, intent(in),    dimension(nx,ny) :: mskr
      integer, intent(in), optional            :: npts,npass

      real,     dimension(nx,ny) :: var_tmp
      integer,  dimension(nx,ny) :: msk_tmp

      integer i,j, is, js, ipass, npoints, nmid, nwgt
      integer ileft,iright,jup,jdown,idleft,idright
      logical left,right,up,down,dleft,dright
      real    val,d1,d2

      real,    allocatable, dimension(:,:) :: subxx
      integer, allocatable, dimension(:,:) :: submr,submd

      integer npasse
      logical debug

      npoints=11
      npasse=5
      if (present(npts))  npoints=npts
      if (present(npass)) npasse=npass

      if (mod(npoints,2) == 0) then
        write(*,*) 'fillm: npoints must not be even'
        call my_abort
      endif

      var_tmp(:,:) = var(:,:)
      msk_tmp(:,:) = mskd(:,:)

      !Add 3 extra passes with reduced number of points to deal with boundaries
      do ipass=1,npasse+3
        nmid=int(real(npoints)/2.)
        allocate(subxx(npoints,npoints),submr(npoints,npoints),submd(npoints,npoints))
        do js=1,ny-npoints+1
        do is=1,nx-npoints+1
          subxx(:,:)=var( is:is+npoints-1,js:js+npoints-1)
          submr(:,:)=mskr(is:is+npoints-1,js:js+npoints-1)
          submd(:,:)=mskd(is:is+npoints-1,js:js+npoints-1)
          if ( submr(nmid,nmid) == 1 .and. submd(nmid,nmid) == 0 ) then
            ! Check horizontal 
            left=.false.
            right=.false.
            j=nmid
            do i=1,nmid-1
              if (submd(i,j)==1) then
                left=.true.
                ileft=i
              endif
            enddo
            do i=npoints,nmid+1,-1
              if (submd(i,j)==1) then
                right=.true.
                iright=i
              endif
            enddo
            ! Check vertical
            up=.false.
            down=.false.
            i=nmid
            do j=1,nmid-1
              if (submd(i,j)==1) then
                down=.true.
                jdown=j
              endif
            enddo
            do j=npoints,nmid+1,-1
              if (submd(i,j)==1) then
                up=.true.
                jup=j
              endif
            enddo
            ! Check diagonal
            dleft=.false.
            dright=.false.
            do i=1,nmid-1
              j=i
              if (submd(i,j)==1) then
                dleft=.true.
                idleft=i
              endif
            enddo  
            do i=npoints,nmid+1,-1
              j=i
              if (submd(i,j)==1) then
                dright=.true.
                idright=i
              endif
            enddo
            nwgt=0
            val=0.
            if ( left .and. right ) then
              nwgt=nwgt+1
              d1=nmid-ileft
              d2=iright-nmid
              val=val+( d1 * subxx(iright,nmid) + d2 * subxx(ileft,nmid) ) / (d1+d2)
            endif
            if ( up .and. down ) then
              nwgt=nwgt+1
              d1=nmid-jdown
              d2=jup-nmid
              val=val+( d1 * subxx(nmid,jup)    + d2 * subxx(nmid,jdown)  ) / (d1+d2)
            endif
            if ( dleft .and. dright ) then
              nwgt=nwgt+1
              d1=nmid-idleft
              d2=idright-nmid
              val=val+( d1 * subxx(idright,idright) + d2 * subxx(idleft,idleft) ) / (d1+d2)
            endif
            if ( nwgt > 0 ) then
              val=val/nwgt
              var_tmp(is+nmid-1,js+nmid-1) =  val
              msk_tmp(is+nmid-1,js+nmid-1) =  1
            endif
          endif
        enddo
        enddo
        var(:,:)  = var_tmp(:,:)
        mskd(:,:) = msk_tmp(:,:)
        deallocate(subxx,submr,submd)
        if ( ipass >= npasse ) then
          npoints = npoints / 2
          if (mod(npoints,2) == 0) npoints = npoints + 1
        endif
      enddo
 
      return
      END SUBROUTINE fillm

      SUBROUTINE mapospread_simple4(data4,M,nx,ny,nsprd)
        implicit none
        integer, intent(in) :: nx,ny,nsprd
        real(kind=4),    intent(out),    dimension(nx,ny)      :: data4
        integer,         intent(out),    dimension(nx,ny)      :: M
        
        real(kind=8), dimension(nx, ny) :: data8
        data8 = data4
        call mapospread_simple8(data8, M, nx, ny, nsprd)
        data4 = real(data8, kind=4)
      END SUBROUTINE mapospread_simple4
  
      SUBROUTINE mapospread_simple8(data,M,nx,ny,nsprd)
        implicit none
        integer, intent(in) :: nx,ny,nsprd
        real(kind=8),    intent(out),    dimension(nx,ny)      :: data
        integer,         intent(out),    dimension(nx,ny)      :: M
  
        integer :: is, i, j
        real(kind=8) :: mw, x1, x2, x3, x4, y1, y2, y3, y4, d1, d2, d3 ,d4, &
                        dtot12, dtot34, dta12, a12, b12, a34, b34, dta34,dtot
        real(kind=4), dimension(nx,ny)      :: Ms, Mtmp
  
        Ms(:,:)=M(:,:)
        Mtmp(:,:)=Ms(:,:)
        DO is=1,nsprd
          DO j=2,ny-1
          DO i=2,nx-1
            mw=Ms(i-1,j)+Ms(i+1,j)+Ms(i,j-1)+Ms(i,j+1)
            if (Ms(i,j).lt.0.5.and.mw.gt.0.5) then
              x1=-1.
              y1= 0.
              x2= 1.
              y2= 0.
              x3= 0.
              y3=-1.
              x4= 0.
              y4= 1.
              d1=Ms(i+1,j)
              d2=Ms(i-1,j)
              d3=Ms(i,j+1)
              d4=Ms(i,j-1)
              dtot12=d1+d2
              dtot34=d3+d4
              dta12=0.
              a12=0.
              b12=0.
              dta34=0.
              a34=0.
              b34=0.
              if (dtot12.gt.0.) then
                a12=d1/dtot12
                b12=d2/dtot12
                dta12=a12*data(i+1,j) + &
       &              b12*data(i-1,j)
              endif
              if (dtot34.gt.0.) then
                a34=d3/dtot34
                b34=d4/dtot34
                dta34=a34*data(i,j+1) + &
       &              b34*data(i,j-1)
              endif
              if (dtot12.eq.0.) then
                data(i,j)=dta34
               elseif (dtot34.eq.0.) then
                data(i,j)=dta12
               else     ! at least 2 points unmasked in i and j axis
                        ! dtot12.gt.0 and dtot34.gt.0
                d1=sqrt(1.d0*(x1**2+y1**2))*Ms(i-1,j) ! Weight average 2, 3 or 4 points
                d2=sqrt(1.d0*(x2**2+y2**2))*Ms(i+1,j)
                d3=sqrt(1.d0*(x3**2+y3**2))*Ms(i,j-1)
                d4=sqrt(1.d0*(x4**2+y4**2))*Ms(i,j+1)
                dtot12=d1+d2 
                if (Ms(i-1,j).le.0.5.or.Ms(i+1,j).le.0.5) dtot12=dtot12*4.
                dtot34=d3+d4 
                if (Ms(i,j-1).le.0.5.or.Ms(i,j+1).le.0.5) dtot34=dtot34*4.
                dtot=dtot12+dtot34
                data(i,j)=dtot34/dtot*dta12 + &
       &                  dtot12/dtot*dta34   !BINGO!
              endif
              Mtmp(i,j)=1.
            endif
          ENDDO
          ENDDO
          Ms(:,:)=Mtmp(:,:)
        ENDDO
  
  ! ... Apply zero gradient at corners and extremities if needed
        i=1
        DO j=1,ny
          if (Ms(i,j).lt.0.5) then
            data(i,j)=Ms(i+1,j)*data(i+1,j)
            Mtmp(i,j)=Ms(i+1,j)
          endif
        ENDDO
        Ms=Mtmp
        i=nx
        DO j=1,ny
          if (Ms(i,j).lt.0.5) then
            data(i,j)=Ms(i-1,j)*data(i-1,j)
            Mtmp(i,j)=Ms(i-1,j)
          endif
        ENDDO
        Ms=Mtmp
        j=1
        DO i=1,nx
          if (Ms(i,j).lt.0.5) then
            data(i,j)=Ms(i,j+1)*data(i,j+1)
            Mtmp(i,j)=Ms(i,j+1)
          endif
        ENDDO
        Ms=Mtmp
        j=ny
        DO i=1,nx
          if (Ms(i,j).lt.0.5) then
            data(i,j)=Ms(i,j-1)*data(i,j-1)
            Mtmp(i,j)=Ms(i,j-1)
          endif
        ENDDO
        Ms=Mtmp
        M(:,:)=nint(Ms(:,:))
  
        return
  
      END SUBROUTINE mapospread_simple8

END MODULE intrp_oce
