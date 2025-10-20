MODULE intrp_itfc
  !!---------------------------------------------------------------------------------------------------
  !!                     ***  MODULE  intrp_itfc  ***
  !!
  !!    ** Purpose : Various methods of interpolation, deals with masks
  !!                 output weights, supports pre-mapped NEMO grids
  !!
  !!    ** Method : Receives std file unit numbers (already opened), basic grid
  !!                and std record information (filtered for needs), and cstintrp options.
  !!                From there intepolate fields and write them in the output
  !!                std unit, called by cstintrp
  !!                
  !!
  !!   history:
  !!         Original : F. Roy (2015)
  !!------------------------------------------------------------------------------------------------------

!! Modules:
  USE intrp_oce
  USE intrp_scrip
  USE intrp_engine
  USE intrp_options
  USE llxy
  USE std
  USE utils
  USE grids
  USE mympp


  IMPLICIT NONE

  
  ! ... file units
  integer,save,public       :: units,unitd,unitr  !File units, source, destination , reference

  !... vector treatment
  character(len=4),   &
          save,public       :: nomvar_v
  logical, save,public      :: vector, ln_mask
                    
  !... ezscint variables
  integer ezqkdef,gdll,ier
  external ezqkdef,gdll
  integer,save,public       :: gdids,gdidd   !Grid definition, source and destination

  !... source grid
  type (type_grid), allocatable, dimension(:) :: grid_src
  type (type_grid), allocatable, dimension(:) :: grid_dst
  integer,save,public       :: nxsrc,nysrc,ngsrc
  type (type_data), allocatable, dimension(:), target   :: dtasrc
  type (type_data), allocatable, dimension(:,:), target :: dtadst
  real(kind=8), allocatable, save,public, dimension(:,:,:)   :: Msrc
                             !i,j,angle source grid
                             !Msrc (used as data to interpolate or aggregate src mask)
  logical, save,public      :: Maskin,Maskou    
                             !Logical to manage masks IO

    !... destination grid
  integer,save,public       ::  nxdst,nydst,ngdst
  
  integer, allocatable, save,public, dimension(:,:,:,:) :: &
                Mdstw
                             !Mask destination grid
                             !Mdstw  (used for intermediate work, and the one written in output file at the end)

! Rotation des vecteurs
  real(kind=4), public, parameter  :: pi = acos(-1.)
  character(len=4),   &
   allocatable, save,public, & 
   dimension(:,:)    :: vecttbl
  integer,save,public       :: nvec

! Poids d'interpolation
  integer,save,public       :: nwgts, nwgtbmax, nwgtb
  type(type_weights), allocatable, dimension(:), target, save, public :: sprd
  type(type_weights), allocatable, dimension(:,:), target, save, public :: intr
  type(type_location),allocatable, dimension(:,:), target, save, public :: locr, loci

! Parametres Yin Yang
  type(std_grid), allocatable, save,public, dimension(:), target &
                            :: subgsrc,subgdst

  PRIVATE
  PUBLIC  intrp_itfc_params, intrp_itfc_init, intrp_itfc_main, intrp_itfc_write, &
          intrp_itfc_tri, intrp_itfc_nearst, intrp_itfc_aggr, intrp_itfc_scrip, &
          intrp_itfc_mixtlc, intrp_itfc_bilc, adjust_periodic, apply_weights, reallocate_intrp, &
          mask_tartinage, apply_weights_vector

CONTAINS

  SUBROUTINE intrp_itfc_params(F_grds,F_grdr)
    !! -----------------------------------------------------------------------------------------------
    !!  *** Initialize cstintrp options, dimensions, parameters, validate with 
    !!      source and destination grids
    !! -----------------------------------------------------------------------------------------------
    IMPLICIT NONE

    ! * Arguments
    type (type_grid), target, intent(inout) :: F_grds,F_grdr

    ! * Local variables
    type(std_grid), pointer :: gs, gr
    character(len=40)   :: cig1r,cig2r,cig3r,cig4r,cnir,cnjr, &
                           cig1s,cig2s,cig3s,cig4s,cnis,cnjs
    
    character(len=*), parameter :: int_fmt = '(I0)'

    !----------------------------------------------------------------------------------------
    ! test all cases of interpolatin, adjust parameter or bail out if in error
    !----------------------------------------------------------------------------------------

    gs => F_grds % std_grd
    gr => F_grdr % std_grd

    if (gs % grtyp == 'Y' .and. gr % grtyp == 'Y' ) then
       write(*,*) 'CSTINTRP: CAN NOT TREAT 2 LAT-LON CLOUD GRIDS'
       call my_abort
    endif

    ! FD: why this would be an issue if the grid is the destination?
    if ((gr % grtyp == 'A' .or. gr % grtyp=='B' .or. gr % grtyp == 'G').and.gr % ig2 == 1) then !Special Gaussian
       print*, 'CSTINTRP: no support for inversed grids'
       print*, 'CSTINTRP: F_grdr%grtyp=',gr % grtyp
       print*, 'CSTINTRP: F_grdr%ig2=',gr % ig2
       call my_abort
    endif

    if (gs % grtyp == 'Y' ) then
      if ( intyp == 'nearst' ) then
        write(*,*) 'CSTINTRP: intyp=nearst is not possible in this context'
        write(*,*) 'CSTINTRP: nearst in CSTINTRP is based on closest bilin point'
        call my_abort
      endif
      write(*,*) 'CSTINTRP: WARNING source grid is cloud, forcing aggr to T'
      intyp = 'aggr'
      write(*,*) 'CSTINTRP: WARNING source grid is cloud, forcing sdfrm to F'
      sdfrm=.false.
    endif
    if (gr % grtyp == 'Y' ) then
      select case(intyp)
      case( 'aggr', 'mixt', 'mixtbc' )
        write(*,*) 'CSTINTRP: WARNING destination grid is cloud, forcing bilin to T'
        write(*,*) 'CSTINTRP: WARNING nearst would also be possible'
        write(*,*) 'CSTINTRP: WARNING bicub  would also be possible'
        intyp = 'bilin'
      end select
      if ( nbdyf>0 ) then
        write(*,*) 'CSTINTRP: WARNING destination grid is cloud, no boundary extension possible'
        write(*,*) 'CSTINTRP: WARNING forcing bdyf=0'
        nbdyf = 0
      endif
    endif
    if (gr % grtyp == 'M' ) then
      select case(intyp)
      case( 'aggr', 'mixt', 'mixtbc' )
        write(*,*) 'CSTINTRP: WARNING destination grid is triangular unstructured, forcing bilin to T'
        write(*,*) 'CSTINTRP: WARNING scrip0 would also be possible'
        intyp = 'bilin'
      end select
      if ( nbdyf>0 ) then
        write(*,*) 'CSTINTRP: WARNING destination grid is triangular unstructured, no boundary extension possible'
        write(*,*) 'CSTINTRP: WARNING forcing bdyf=0'
        nbdyf = 0
      endif
    endif

    select case(intyp)
    case('bicub','mixtbc')
      write(*,*) 'CSTINTRP: BI-CUBIC INTERPOLATION'
      nwgtb=12
    case('bilin','mixt')
      write(*,*) 'CSTINTRP: BI-LINEAR INTERPOLATION'
      nwgtb=4
    case('nearst')
      write(*,*) 'CSTINTRP: NEARST INTERPOLATION'
      nwgtb=1
    end select

    !----------------------------------------------------------------------------------------
    ! report the actual parameters that are going to be used
    !----------------------------------------------------------------------------------------

    write(*,*) '=========================================='
    write(*,*) 'CSTINTRP: Using the following parameters'
    write(*,*) 'CSTINTRP:   intyp = ',TRIM(intyp)
    if (intyp=='nearst') &
    write(*,*) 'CSTINTRP:     !!!Note nearst is based on choosing the largest bilin weight'
    write(*,*) 'CSTINTRP:   nadis=',nadis
    write(*,*) 'CSTINTRP:   mtrsh=',mtrsh
    write(*,*) 'CSTINTRP:   nsprd=',nsprd
    write(*,*) 'CSTINTRP:   sdfrm=',sdfrm
    write(*,*) 'CSTINTRP:   nbdyf=',nbdyf
    write(*,*) 'CSTINTRP:   npak=',npak_ou
    write(*,*) 'CSTINTRP:   naggrmax=',naggrmax
    write(*,*) 'CSTINTRP:   owgts=',owgts
    write(*,*) 'CSTINTRP:   uwgdir=',trim(uwgdir)
    write(*,*) 'CSTINTRP:   zland=',zland
    write(*,*) 'CSTINTRP:   lval =',landval
    write(*,*) 'CSTINTRP:   maskd =',maskd
    write(*,*) 'CSTINTRP: periodicity setup'
    write(*,*) 'CSTINTRP:   snbdy=',F_grds % nbdy
    write(*,*) 'CSTINTRP:   snrep=',F_grds % nrepeat
    write(*,*) 'CSTINTRP:   sebdy=',F_grds % webdy
    write(*,*) 'CSTINTRP:   serep=',F_grds % werepeat
    write(*,*) 'CSTINTRP:   dnbdy=',F_grdr % nbdy
    write(*,*) 'CSTINTRP:   dnrep=',F_grdr % nrepeat
    write(*,*) 'CSTINTRP:   debdy=',F_grdr % webdy
    write(*,*) 'CSTINTRP:   derep=',F_grdr % werepeat
    write(*,*) '=========================================='

    ! find the name of weight file if needed
    if (iowgt) then

      write(cig1s, int_fmt) gs%ig1
      write(cig2s, int_fmt) gs%ig2
      write(cig3s, int_fmt) gs%ig3
      write(cig4s, int_fmt) gs%ig4
      write(cnis, int_fmt)  gs%ni
      write(cnjs, int_fmt)  gs%nj

      write(cig1r, int_fmt) gr%ig1
      write(cig2r, int_fmt) gr%ig2
      write(cig3r, int_fmt) gr%ig3
      write(cig4r, int_fmt) gr%ig4
      write(cnir, int_fmt)  gr%ni
      write(cnjr, int_fmt)  gr%nj

      gnams=gs%grtyp// &
         '_'//trim(cig1s)//'_'//trim(cig2s)//'_'//trim(cig3s)//'_'//trim(cig4s)// &
         '_'//trim(cnis)//'_'//trim(cnjs)
      gnamr=gr%grtyp// &
         '_'//trim(cig1r)//'_'//trim(cig2r)//'_'//trim(cig3r)//'_'//trim(cig4r)// &
         '_'//trim(cnir)//'_'//trim(cnjr)
    endif

  END SUBROUTINE intrp_itfc_params

  SUBROUTINE intrp_itfc_init(F_grds,F_grdr,F_nomvar_ref)
    !! -----------------------------------------------------------------------------------------------
    !!  *** Allocate memory and copy grid descriptors to destination file
    !! -----------------------------------------------------------------------------------------------
    IMPLICIT NONE

    ! * Arguments
    type (type_grid), target &
                         :: F_grds,F_grdr
    character(len=4)     :: F_nomvar_ref

    ! * Local variable
    type (std_variable),  &
       dimension(:),      &
       allocatable      :: varr,    &  !All variables from reference file
                           varr_p      !Positional parameters
    type (std_variable) :: varr_pz     !Split YY tmp Z positional parameters
    integer             :: nvr, nvr_p, i,j
    integer             :: ivr, ivr_p, status_st90, ierr, igsrc, igdst
    logical             :: found, Maskouw
    type (std_grid), pointer :: gr
    type (type_grid), pointer :: subs, subr
    character(len=1) :: canum
    type(type_location), pointer :: plocr, ploci
    real(kind=4), allocatable, dimension(:,:,:)     :: wrkr4
    integer,      allocatable, dimension(:,:,:)     :: wrki4

    Maskouw = .false.
    ngsrc = F_grds % ng
    ngdst = F_grdr % ng
    nxsrc = F_grds % nx
    nysrc = F_grds % ny
    nxdst = F_grdr % nx
    nydst = F_grdr % ny
    if ( F_grds % std_grd % grtyp == 'U' ) nysrc = nysrc / 2
    if ( F_grdr % std_grd % grtyp == 'U' ) nydst = nydst / 2

    call allocate_dst
    call allocate_src

! Calcul des xy de la grille source sur la destination

    do igsrc = 1, ngsrc
    do igdst = 1, ngdst

      subs => F_grds % subgrd(igsrc)
      subr => F_grdr % subgrd(igdst)
      subgsrc( igsrc ) = subs % std_grd
      subgdst( igdst ) = subr % std_grd
      plocr => locr(igdst,igsrc)
      ploci => loci(igsrc,igdst)

      ! do the localization of the destination points on the source grid
      select case(intyp)
      case( 'bilin', 'bicub', 'mixt', 'mixtbc', 'nearst')
        write(*,*) 'CSTINTRP: CALCUL SUR LA GRILLE DESTINATION DES XY DU REPERE SOURCE'
        call ll2xy( subs, subr, plocr )
        if ( any( .not. plocr % mg) ) Maskouw = .true.
        if (Maskouw) Maskou=.true.
             write(*,*) 'forward localization'
             write(*,*) 'min/max I',minval( plocr % ig(:,:) ),maxval( plocr % ig(:,:) )
             write(*,*) 'min/max J',minval( plocr % jg(:,:) ),maxval( plocr % jg(:,:) )
             write(*,*) 'min/max A',minval( plocr % ag(:,:) ),maxval( plocr % ag(:,:) )
             write(*,*) 'min/max B',minval( plocr % bg(:,:) ),maxval( plocr % bg(:,:) )
             write(*,*) 'min/max a',minval( plocr % pg(:,:) ),maxval( plocr % pg(:,:) )
      end select

      ! do reverse localization used in aggregation
      select case(intyp)
      case( 'aggr', 'mixt', 'mixtbc')
        write(*,*) 'CSTINTRP: CALCUL SUR LA GRILLE SOURCE DES XY DU REPERE DESTINATION'
        call ll2xy( subr, subs, ploci )
             write(*,*) 'backward localization'
             write(*,*) 'min/max I',minval( ploci % ig(:,:) ),maxval( ploci % ig(:,:) )
             write(*,*) 'min/max J',minval( ploci % jg(:,:) ),maxval( ploci % jg(:,:) )
             write(*,*) 'min/max A',minval( ploci % ag(:,:) ),maxval( ploci % ag(:,:) )
             write(*,*) 'min/max B',minval( ploci % bg(:,:) ),maxval( ploci % bg(:,:) )
             write(*,*) 'min/max a',minval( ploci % pg(:,:) ),maxval( ploci % pg(:,:) )
      end select

    enddo   ! igdst
    enddo   ! igsrc

    ! get all variables from the reference file and isolate the grid variables
    nvr = getstdnvar(unitr,ierr_st90=status_st90)
    if ( nvr <= 0 ) then
      write(*,*) 'CSTINTRP getstdnvar => problem reading reference file variables'
      write(*,*) 'CSTINTRP getstdnvar => nvr,status_st90=',nvr,status_st90
      call my_abort
    endif
    allocate(varr(nvr), varr_p(nvr), stat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTINTRP varr, varr_p allocation problem'
      write(*,*) 'CSTINTRP nvr=',nvr
      call my_abort
    endif
    ierr = fillstdvar(varr,nvr,unitr,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTINTRP fillstdvar => problem filling reference variable informations'
      write(*,*) 'CSTINTRP fillstdvar => ierr,status_st90=',ierr,status_st90
      call my_abort
    endif
    ! recover the grid descriptors (>>, ^^ or more) of the grid by finding variables having the igs as ips.
    gr => F_grdr % std_grd
    ierr = grabgdescstdvar(varr_p,nvr_p, varr, nvr,gr % grtyp, gr % ig1, gr % ig2, gr % ig3, gr % ig4)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTINTRP grabgdescstdvar => problem getting positional paramters'
      write(*,*) 'CSTINTRP grabgdescstdvar => ierr=',ierr
      write(*,*) 'CSTINTRP grabgdescstdvar => F_grdr=',gr
      call my_abort
    endif

    ! copy the grid descriptor from the reference file to the destination file
    do ivr_p=1,nvr_p  !possibly nothing, e.g. PS grid
      if ( varr_p(ivr_p)%nomvar == '^>' .and. splt_yy ) then
        do igdst = 1, ngdst
           subr => F_grdr % subgrd(igdst)
           select case (igdst)
           case(1)
              varr_p(ivr_p) % etiket = 'YinDst'
           case(2)
              varr_p(ivr_p) % etiket = 'YanDst'
           end select
           call write_pos_axis(subr % std_grd, subr % aux_grd, varr_p(ivr_p), '>>', subr % lld % tic, unitd, igdst)
           call write_pos_axis(subr % std_grd, subr % aux_grd, varr_p(ivr_p), '^^', subr % lld % tac, unitd, igdst)
           subgdst(igdst) = subr % std_grd
        enddo
      else
        ierr = xferstdvar(varr_p(ivr_p:ivr_p),1,unitr,unitd)
      endif
    enddo

! Ecriture du champs d'angle et lat-lon (optionel avec les poids d'interpolation)
! Yin Yang:
!   On attend de voir l'interface de couplage YY avant d'assembler
!   Eventuellement extraire les >> ^^ des grilles Z pour le mode separe

    if (owgts) then

      ! Trouve les attributs avec F_nomvar_ref et la grille en question
      found = .false.
      do ivr = 1, nvr 
        if ( trim(varr(ivr)%nomvar) == trim(F_nomvar_ref) .and. &
             varr(ivr)%typvar /= '@@'                     .and. &
             varr(ivr)%grtyp == gr%grtyp              .and. &
             varr(ivr)%ig1   == gr%ig1                .and. &
             varr(ivr)%ig2   == gr%ig2                .and. &
             varr(ivr)%ig3   == gr%ig3                .and. &
             varr(ivr)%ig4   == gr%ig4                .and. &
             varr(ivr)%ni    == gr%ni                 .and. &
             varr(ivr)%nj    == gr%nj) then
          varr_p(1) = varr(ivr)  ! varr_p est une variable temporaire ici
          found = .true.
        endif
      enddo

      if (.not. found) then
        write(*,*) 'CSTINTRP attributes not found in reference file'
        write(*,*) 'CSTINTRP for F_nomvar_ref=',trim(F_nomvar_ref)
        call my_abort
      endif

      ! write angles to destination file
      varr_p(1) % nomvar(1:3)='ANG'
      varr_p(1) % ip1=0
      varr_p(1) % ip2=0
      varr_p(1) % ip3=0
      varr_p(1) % etiket='WEIGHT'

      do igsrc = 1 , ngsrc 
       varr_p(1) % typvar='P@'
       varr_p(1) % nbits=32
       varr_p(1) % datyp=5
       varr_p(1) % xtype=2
       write(canum,'(i1)') igsrc
       if ( splt_yy .and. ngsrc == 2) varr_p(1)%etiket(7:10)='Src'//canum
       varr_p(1) % ip3 = varr_p(1) % ip3 + igsrc - 1
       allocate(wrkr4(nxdst,nydst,ngdst)); wrkr4 = 0.
       do igdst = 1, ngdst
         if (associated(locr(igdst,igsrc) % pg)) & 
               wrkr4(:,:,igdst) = locr(igdst,igsrc) % pg(:,:) * 180./pi
       enddo
       call write_var_std(ngdst, subgdst, varr_p(1), unitd, f2d=wrkr4 )
       deallocate(wrkr4)
       varr_p(1) % typvar='@@'
       varr_p(1) % nbits=1
       varr_p(1) % datyp=2
       varr_p(1) % xtype=1
       allocate(wrki4(nxdst,nydst,ngdst)); wrki4(:,:,:) = 0
       do igdst = 1, ngdst
         if (associated(locr(igdst,igsrc) % mg)) & 
           where( locr(igdst,igsrc) % mg(:,:) ) wrki4(:,:,igdst) = 1
       enddo
       call write_var_std(ngdst, subgdst, varr_p(1), unitd, i2d=wrki4 )
       deallocate(wrki4)
      enddo

      ! write 2D fields of lat/lon
      varr_p(1) % ig3 = gr % ig3 ! might be overwritten
      varr_p(1) % ip3 = 0
      varr_p(1) % etiket='WEIGHT'
      varr_p(1) % typvar='P'
      varr_p(1) % nbits=32
      varr_p(1) % datyp=5
      varr_p(1) % xtype=2
      varr_p(1) % nomvar='LAT'
      allocate(wrkr4(nxdst,nydst,ngdst))
      do igdst = 1, ngdst
         subr => F_grdr % subgrd(igdst)
         wrkr4(:,:,igdst) = subr % lld % latt(:,:)
      enddo
      call write_var_std(ngdst, subgdst, varr_p(1), unitd, f2d=wrkr4 )
      varr_p(1) % ig3 = gr % ig3 ! might be overwritten
      varr_p(1) % ip3 = 0
      varr_p(1) % nomvar='LON'
      do igdst = 1, ngdst
         subr => F_grdr % subgrd(igdst)
         wrkr4(:,:,igdst) = subr % lld % lont(:,:)
      enddo
      call write_var_std(ngdst, subgdst, varr_p(1), unitd, f2d=wrkr4 )
      deallocate(wrkr4)

    endif  ! owgts

  END SUBROUTINE intrp_itfc_init

  SUBROUTINE allocate_dst
    !! -----------------------------------------------------------------------------------------------
    !!  *** allocate destination/ref arrays
    !! -----------------------------------------------------------------------------------------------
    if (allocated(Mdstw)) deallocate(Mdstw)
    allocate(Mdstw(nxdst,nydst,ngdst,ngsrc))
    if (allocated(dtadst)) deallocate(dtadst)
    allocate( dtadst(ngdst,ngsrc))

    Maskou=.false.

    if (allocated(subgdst)) deallocate(subgdst)
    allocate( subgdst(ngdst) )

    if (allocated(intr)) then
       call deallocate_intrp(intr)
       deallocate(intr)
    endif
    allocate( intr(ngdst,ngsrc) )

    if (allocated(locr)) then
       call deallocate_loc(locr)
       deallocate(locr)
    endif
    allocate( locr(ngdst,ngsrc) )

    if (allocated(loci)) then
       call deallocate_loc(loci)
       deallocate(loci)
    endif
    allocate( loci(ngsrc,ngdst) )

  END SUBROUTINE allocate_dst

  SUBROUTINE allocate_src
    !! -----------------------------------------------------------------------------------------------
    !!  *** allocate source arrays
    !! -----------------------------------------------------------------------------------------------
    if (allocated(dtasrc)) deallocate(Msrc, dtasrc)
    allocate(  Msrc(nxsrc,nysrc,ngsrc), dtasrc(ngsrc))

    if (allocated(subgsrc)) deallocate(subgsrc)
    allocate( subgsrc(ngsrc) )

  END SUBROUTINE allocate_src

  SUBROUTINE write_pos_axis(grd, gra, varp, namex, axis1d, unitf, ig )
    !! -----------------------------------------------------------------------------------------------
    !!  *** write to file the tic/tac for separate YY output
    !! -----------------------------------------------------------------------------------------------
    IMPLICIT NONE
    ! arguments
    type(std_grid), intent(inout) :: grd, gra
    type(std_variable), intent(in)  :: varp
    character(len=*) :: namex
    real(kind=4), dimension(*) :: axis1d
    integer unitf
    ! locals
    type(std_variable)  :: varr_pz
    integer ierr, status_st90, ig
    real(kind=4), allocatable, dimension(:,:) :: array2d

    varr_pz = varp
    varr_pz % nomvar = namex
    varr_pz % ig1 = gra % ig1
    varr_pz % ig2 = gra % ig2
    varr_pz % ig3 = gra % ig3
    varr_pz % ig4 = gra % ig4
    varr_pz % grtyp = gra % grtyp
    varr_pz % xtype = 2 !Force E32 type

    varr_pz % ip3 = varr_pz % ip3 + (ig - 1)

    select case (namex)
    case('>>')
      varr_pz % ni = grd % ni
      varr_pz % nj = 1
    case('^^')
      varr_pz % ni = 1
      varr_pz % nj = grd % nj
    end select 

    allocate( array2d(varr_pz % ni, varr_pz % nj) )
    select case (namex)
    case('>>')
      array2d(:,1) = axis1d(1:varr_pz % ni)
    case('^^')
      array2d(1,:) = axis1d(1:varr_pz % nj)
    end select 

    ! ... other parameters do not change compared to ^> ...
    ! ... will add 1 to ip3 for second grid YY grid
    ierr = putstdvar(unitf, array2d, varr_pz % ni, varr_pz % nj, &
                         varr_pz,ierr_st90=status_st90)
    ! make sure the igs of the output grid are correct
    grd % ig1 = varr_pz % ip1
    grd % ig2 = varr_pz % ip2
    grd % ig3 = varr_pz % ip3
    grd % ig4 = varr_pz % deet
    if ( ierr /= 0 ) then
      write(*,*) 'CSTINTRP putstdvar problem => ',namex
      write(*,*) 'CSTINTRP ierr,stdvar=',ierr,varp
      write(*,*) 'CSTINTRP status_st90=',status_st90
      call my_abort
    endif
    deallocate( array2d)

  END SUBROUTINE write_pos_axis

  SUBROUTINE intrp_itfc_main(var, F_grds, F_grdr, F_mask_change)
    !! -----------------------------------------------------------------------------------------------
    !!  *** Performs interpolation using variable and grid attributes
    !! -----------------------------------------------------------------------------------------------
    ! * Arguments

    IMPLICIT NONE
    type (variable), intent(in) :: var
    type (type_grid), intent(in), target &
                         :: F_grds,F_grdr
    logical, intent(in)  :: F_mask_change

    ! * Local variable

    integer              :: igsrc,igdst
    logical              :: fullcal
    type (type_grid), pointer :: subs, subr
    type(type_weights), pointer :: pintr, psprd
    type(type_location), pointer :: plocr, ploci
    character(len=2048) :: fwgt
    character(len=1)   :: cigsrc,cigdst
    character(len=40)  :: cip1
    type(type_data), pointer :: pdst, psrc

!-------------------------------------------------------
! Lecture du champ source et possiblement de son masque
!-------------------------------------------------------

    Maskin=.false. 
    Msrc(:,:,:)=1.
    vector  = var % ln_vector
    ln_mask = var % ln_mask ! found a mask

    call get_data(var, F_grds, units)

!-------------------------------------------------------
! Ajuste la methode de calcul si le masque change
! pour forcer un calcul complet (nouveaux poids)
! Tant que fullcal est a false, les poids restent 
! en memoire
!-------------------------------------------------------

    fullcal=F_mask_change

!-------------------------------------------------------
! Tartinage des regions masquees au besoin
!-------------------------------------------------------

    if (fullcal) then
      call mask_tartinage(F_grds)
    endif

!-------------------------------------------------------
! Interpolation du champs de donnees selon la methode
!-------------------------------------------------------

    do igsrc=1,ngsrc
      psprd => sprd(igsrc)
      subs => F_grds % subgrd(igsrc)
      psrc  => dtasrc(igsrc)
    do igdst=1,ngdst
      fullcal = F_mask_change
      subr => F_grdr % subgrd(igdst)
      pintr => intr(igdst, igsrc)
      plocr => locr(igdst, igsrc)
      ploci => loci(igsrc, igdst)
      pdst  => dtadst(igdst, igsrc)
      allocate( pdst % u (nxdst, nydst) )
      if (vector) &
        allocate( pdst % v (nxdst, nydst) )

      ! generate the name of the weight file if needed
      if (iowgt) then
        write(cigsrc,'(i0)') igsrc
        write(cigdst,'(i0)') igdst
        if ( ln_mask ) then
          write(cip1, '(i0)') var % varstd % ip1
          fwgt=trim(uwgdir)//'/wgt_'//trim(gnams)//'_'//trim(gnamr)//'_'//trim(cip1)// &
                '_'//cigsrc//'_'//cigdst
        else
          fwgt=trim(uwgdir)//'/wgt_'//trim(gnams)//'_'//trim(gnamr)// &
                '_'//cigsrc//'_'//cigdst
        endif
      endif

      ! allocate arrays and load weights if saved on disk
      if (fullcal) call reallocate_intrp( pintr )

      !fullcal is modified only if iowgt=.true. (owgts=.false.) and file is found
      if ( fullcal .and. iowgt ) then
        call intrp_getput_premap_wgt(fullcal,'get', pintr, fwgt)
           !!! Note, if two different masks with same ip1 ... no solution for now
      endif

      if (fullcal) then

        write(*,*) 'CSTINTRP: FIRST FIELD ... FULL COMPUTATION'

        if ( subs % std_grd % grtyp == 'M') then ! unstructured grids
          SELECT CASE ( intyp )
          CASE('nearst', 'bilin') 
            call intrp_itfc_tri(plocr, pintr, subs, Msrc(:,:,igsrc), nearst=intyp == 'nearst')
          CASE('aggr') ! aggr
            call intrp_itfc_aggr(ploci, pintr, subs, subr, Msrc(:,:,igsrc))
          CASE('scrip0')
            call intrp_itfc_scrip(plocr, ploci, pintr, subs, subr, Msrc(:,:,igsrc), .false.)
          END SELECT
        else ! regular grids
          SELECT CASE ( intyp )
          CASE('nearst') 
            call intrp_itfc_nearst(plocr, pintr, psprd, Msrc(:,:,igsrc))
          CASE('bilin','bicub') ! bilin (and possibly nearest) or bicub
            call set_grid(subs)
            call intrp_itfc_bilc(plocr, pintr, psprd, Msrc(:,:,igsrc))
          CASE('aggr') ! aggr
            call intrp_itfc_aggr(ploci, pintr, subs, subr, Msrc(:,:,igsrc))
          CASE('mixt','mixtbc') ! mixt or mixtbc
            call intrp_itfc_mixtlc(plocr, ploci, pintr, psprd, subs, subr, Msrc(:,:,igsrc))
          CASE('scrip0','scrip1')
            call intrp_itfc_scrip(plocr, ploci, pintr, subs, subr, Msrc(:,:,igsrc), intyp=='scrip1')
          END SELECT
        endif

        call adjust_periodic(subs, subr, pintr )

        write(*,*) 'weights', pintr % numwgt
        write(*,*) 'min/max W',minval( pintr %  wgt ),maxval( pintr %  wgt )
        write(*,*) 'min/max I',minval( pintr % iwgt ),maxval( pintr % iwgt )
        write(*,*) 'min/max J',minval( pintr % jwgt ),maxval( pintr % jwgt )
        write(*,*) 'min/max N',minval( pintr % nwgt ),maxval( pintr % nwgt )
        
        if ( iowgt ) then
          call intrp_getput_premap_wgt(fullcal,'put', pintr, fwgt)
        !!! Note, if two different masks with same ip1 ... no solution for now
        endif

      endif ! fullcal


      ! Initialize working mask
      Mdstw(:,:,igdst,igsrc) = 0

      ! mask points masked in the reference field
      if (maskd) then
        where (.not. subr % lld % mask) pintr % mwgt(:,:) = .false.
      endif

      ! ... Mdstw is the one written
      where( pintr % mwgt(:,:) ) Mdstw(:,:,igdst,igsrc) = 1 

      ! flag that would trigger writing the output mask 
      ! if not done already or if no input mask
      if ( any( .not. pintr % mwgt(:,:)) ) Maskou = .true.  
      

!...  Application des poids

      if (vector) then
        write(*,*) 'CSTINTRP: TREATING VECTOR  '
        call apply_weights_vector( psrc, pdst, subs, subr, pintr)
        write(*,*) 'min/max U',minval(pdst % u ),maxval(pdst % u )
        write(*,*) 'min/max V',minval(pdst % v ),maxval(pdst % v )
      else
        call apply_weights( pdst % u, psrc % u, pintr)
        write(*,*) 'min/max F',minval(pdst % u ),maxval(pdst % u )
      endif
    enddo
    enddo

    ! write interpolated and optinoal weights
    call intrp_itfc_write(var,F_grdr,F_mask_change)

  END SUBROUTINE intrp_itfc_main

  SUBROUTINE get_data(var, grid, unitf)
  !-------------------------------------------------------------------
  ! read data from from file for a particular variable, ip1 including mask and v-component
  !-------------------------------------------------------------------
    implicit none
    ! arguments
    type(variable) , target :: var
    type(type_grid), target :: grid
    integer unitf
    ! locals
    type(std_variable), pointer :: F_vars, V_vars, F_varm
    type(std_grid), pointer :: grd
    integer              :: status_st90, ierr
    real(kind=8), allocatable, &
      dimension(:,:)     :: wrkr8
    integer, allocatable, &
      dimension(:,:)     :: wrki4
    integer igsrc

    F_vars => var % varstd
    V_vars => var % varstdv
    F_varm => var % mask

    nomvar_v = V_vars % nomvar

    grd => grid % std_grd
    if (grd % grtyp == 'U') then
      allocate(wrkr8(nxsrc,nysrc*2))
      ierr = getstdvar(wrkr8,nxsrc,nysrc*2,unitf,   &
                       F_vars,ierr_st90=status_st90)
      if ( ierr /= 0 ) then
        write(*,*) 'CSTINTRP getstdvar problem => wrkr8'
        write(*,*) 'CSTINTRP ierr,F_vars=',ierr,F_vars
        write(*,*) 'CSTINTRP status_st90=',status_st90
        call my_abort
      endif
      do igsrc=1,ngsrc
        allocate( dtasrc(igsrc) % u (nxsrc,nysrc) )
        dtasrc(igsrc) % u(:,:) = wrkr8(:,nysrc*(igsrc-1)+1:nysrc*(igsrc-1)+nysrc)
      enddo
      if (vector) then
        ierr = getstdvar(wrkr8,nxsrc,nysrc*2,unitf,V_vars,ierr_st90=status_st90)
        if ( ierr /= 0 ) then
          write(*,*) 'CSTINTRP getstdvar problem => wrkr8'
          write(*,*) 'CSTINTRP ierr,V_vars=',ierr,V_vars
          write(*,*) 'CSTINTRP status_st90=',status_st90
          call my_abort
        endif
        do igsrc=1,ngsrc
          allocate( dtasrc(igsrc) % v (nxsrc,nysrc) )
          dtasrc(igsrc) % v(:,:) = wrkr8(:,nysrc*(igsrc-1)+1:nysrc*(igsrc-1)+nysrc)
        enddo
      endif
      deallocate(wrkr8)
     else
      igsrc = 1
      allocate( dtasrc(igsrc) % u (nxsrc,nysrc) )
      ierr = getstdvar(dtasrc(igsrc) % u, nxsrc,nysrc,unitf,F_vars,ierr_st90=status_st90)
      if ( ierr /= 0 ) then
        write(*,*) 'CSTINTRP getstdvar problem => dtasrc main field'
        write(*,*) 'CSTINTRP ierr,F_vars=',ierr,F_vars
        write(*,*) 'CSTINTRP status_st90=',status_st90
        call my_abort
      endif
      if (vector) then
        allocate( dtasrc(igsrc) % v (nxsrc,nysrc) )
        ierr = getstdvar(dtasrc(igsrc) % v, nxsrc,nysrc,unitf,V_vars,ierr_st90=status_st90)
        if ( ierr /= 0 ) then
          write(*,*) 'CSTINTRP getstdvar problem => dtasrc v field'
          write(*,*) 'CSTINTRP ierr,V_vars=',ierr,V_vars
          write(*,*) 'CSTINTRP status_st90=',status_st90
          call my_abort
        endif
      endif
    endif 

    if (ln_mask) then
      if (grd % grtyp == 'U') then
        allocate(wrki4(nxsrc,nysrc*2))
        ierr = getstdvar(wrki4,nxsrc,nysrc*2,unitf,F_varm,ierr_st90=status_st90)
        if ( ierr /= 0 ) then
          write(*,*) 'CSTINTRP getstdvar problem => wrki4'
          write(*,*) 'CSTINTRP ierr,F_varm=',ierr,F_varm
          write(*,*) 'CSTINTRP status_st90=',status_st90
          call my_abort
        endif
        do igsrc=1,ngsrc
          allocate( dtasrc(igsrc) % m (nxsrc,nysrc) )
          dtasrc(igsrc) % m(:,:) = .false.
          where( wrki4(:,nysrc*(igsrc-1)+1:nysrc*(igsrc-1)+nysrc) == 1 ) dtasrc(igsrc) % m(:,:) = .true.
        enddo
        deallocate(wrki4)
      else
        igsrc=1
          allocate( wrki4(nxsrc,nysrc) )
          ierr = getstdvar(wrki4,nxsrc,nysrc,unitf,F_varm,ierr_st90=status_st90)
          if ( ierr /= 0 ) then
            write(*,*) 'CSTINTRP getstdvar problem => mask source'
            write(*,*) 'CSTINTRP ierr,F_varm=',ierr,F_varm
            write(*,*) 'CSTINTRP status_st90=',status_st90
            call my_abort
          endif
          allocate( dtasrc(igsrc) % m (nxsrc,nysrc) )
          dtasrc(igsrc) % m(:,:) = .false.
          where( wrki4(:,:) ==1 ) dtasrc(igsrc) % m(:,:) = .true.
          deallocate(wrki4)
      endif
      do igsrc=1,ngsrc
        where( .not.dtasrc(igsrc) % m(:,:)) Msrc(:,:,igsrc) = 0
      enddo
      Maskin=.true.
    endif

  END SUBROUTINE get_data

  SUBROUTINE deallocate_loc(loc)
    ! arguments
    type(type_location), dimension(:,:), intent(inout) :: loc
    ! locals
    integer n1, n2, i1, i2
    n1 = SIZE( loc, 1 )
    n2 = SIZE( loc, 2 )
    write(*,*) 'loc dim',n1,n2
    do i2 = 1, n2
       do i1 = 1, n1
          if (associated( loc(i1,i2) % ag ) ) &
             deallocate( loc(i1,i2) % ag &
                       , loc(i1,i2) % bg &
                       , loc(i1,i2) % ig &
                       , loc(i1,i2) % jg &
                       , loc(i1,i2) % pg &
                       , loc(i1,i2) % mg )
       enddo
    enddo
  END SUBROUTINE deallocate_loc

  SUBROUTINE deallocate_intrp(pintr)
    ! arguments
    type(type_weights), dimension(:,:), intent(inout) :: pintr
    ! locals
    integer n1, n2, i1, i2
    n1 = SIZE( pintr, 1 )
    n2 = SIZE( pintr, 2 )
    write(*,*) 'intr dim',n1,n2
    do i2 = 1, n2
       do i1 = 1, n1
          call deallocate_weights( pintr(i1,i2) )
       enddo
    enddo
  END SUBROUTINE deallocate_intrp

  SUBROUTINE reallocate_intrp(pintr)
    ! arguments
    type(type_weights), intent(inout) :: pintr

    ! set max size for bilin, nearst, bicubic

    call deallocate_weights( pintr )

    pintr % nx = nxdst
    pintr % ny = nydst

    select case(intyp)
    case('bilin', 'mixt', 'bicub', 'mixtbc')
      nwgtbmax = max(nwgtb * nwgts, nwgtb)
    case default
      nwgtbmax = 1
    end select

    SELECT CASE(intyp)
    CASE('nearst','bilin','bicub')
      pintr % numwgt = nwgtbmax
    CASE('aggr','mixt','mixtbc','scrip0','scrip1')
      pintr % numwgt = naggrmax
    END SELECT
    call allocate_weights( pintr )

    ! redefine numwgt to its value before spreading
    SELECT CASE(intyp)
    CASE('nearst','bilin','bicub')
      pintr % numwgt = nwgtb
    END SELECT

  END SUBROUTINE reallocate_intrp

  SUBROUTINE reallocate_sprd( ng )
    ! arguments
    integer, intent(in) :: ng
    ! local
    type(type_weights), pointer :: spr
    integer ngdim, ig

    ! if already allocated, destroy existing structure
    if (allocated(sprd)) then
      ngdim = size(sprd)
      do ig = 1, ngdim
        spr => sprd( ig )
        call deallocate_weights( spr )
      enddo
      deallocate( sprd) 
    endif

    ! allocate anew
    allocate( sprd(ng) )
    do ig = 1, ng
      spr => sprd( ig )
      spr % numwgt = nwgts
      spr % nx = nxsrc
      spr % ny = nysrc
      call allocate_weights( spr )
    enddo

  END SUBROUTINE reallocate_sprd

  SUBROUTINE deallocate_weights( pintr )
    ! arguments
    type(type_weights), target, intent(inout) :: pintr

    if (allocated(pintr % wgt) ) &
             deallocate(pintr %  wgt, &
                        pintr % iwgt, &
                        pintr % jwgt, &
                        pintr % nwgt, &
                        pintr % mwgt)

  END SUBROUTINE deallocate_weights

  SUBROUTINE allocate_weights( pintr )
    ! arguments
    type(type_weights), intent(inout) :: pintr
    ! locals
    integer nwgloc, nx, ny

    nx = pintr % nx
    ny = pintr % ny
    nwgloc = pintr % numwgt

    allocate( pintr % nwgt(nx,ny), &
              pintr % mwgt(nx,ny) )

    pintr % nwgt(:,:) = 0
    pintr % mwgt(:,:) = .false.

    allocate(pintr %  wgt(nwgloc,nx,ny), &
             pintr % iwgt(nwgloc,nx,ny), &
             pintr % jwgt(nwgloc,nx,ny))

    pintr %  wgt(:,:,:) = zerd
    pintr % iwgt(:,:,:) = 1
    pintr % jwgt(:,:,:) = 1

  END SUBROUTINE allocate_weights

  SUBROUTINE adjust_periodic(grds, grdr, pintr )
    ! arguments
    type(type_grid) :: grds, grdr
    type(type_weights), target :: pintr
    ! locals
    integer, dimension(:,:,:), pointer :: iwgt, jwgt
    integer, dimension(:,:), pointer :: nwgt
    logical, dimension(:,:), pointer :: mwgt
    real(kind=8), dimension(:,:,:), pointer :: wgt
    integer nx, ny, i, j, im, jm

    iwgt => pintr % iwgt
    jwgt => pintr % jwgt
    nwgt => pintr % nwgt
    mwgt => pintr % mwgt
     wgt => pintr %  wgt
    nx   =  pintr % nx
    ny   =  pintr % ny

    ! still necessary?
      ! Global GEM grid, f(ni)=f(1)
      ! This will avoid array overflow in ocean model
      ! interpolating GEM fields
      if (grds % webdy .and. grds % std_grd % grtyp == 'Z') then
        where (iwgt(:,:,:) == nxsrc) iwgt(:,:,:) = 1
      endif

    ! start with west-east
    if( grdr % webdy ) then ! w-e periodicity active
      select case( grdr % werepeat )
      case(1)
        iwgt(:,nx,:) = iwgt(:,1,:)
        jwgt(:,nx,:) = jwgt(:,1,:)
        nwgt(  nx,:) = nwgt(  1,:)
        mwgt(  nx,:) = mwgt(  1,:)
         wgt(:,nx,:) =  wgt(:,1,:)
      case(2)
        iwgt(:, 1,:) = iwgt(:,nx-1,:)
        jwgt(:, 1,:) = jwgt(:,nx-1,:)
        nwgt(   1,:) = nwgt(  nx-1,:)
        mwgt(   1,:) = mwgt(  nx-1,:)
         wgt(:, 1,:) =  wgt(:,nx-1,:)
        iwgt(:,nx,:) = iwgt(:,   2,:)
        jwgt(:,nx,:) = jwgt(:,   2,:)
        nwgt(  nx,:) = nwgt(     2,:)
        mwgt(  nx,:) = mwgt(     2,:)
         wgt(:,nx,:) =  wgt(:,   2,:)
      end select
    endif

    select case( grdr % nbdy )
    case('F') ! F-folding
      do j = ny, ny - grdr % nrepeat + 1, -1
        jm = j - grdr % nrepeat
        do i = grdr % werepeat, nx
          im = nx - i + 1
          iwgt(:,i,j) = iwgt(:,im,jm)
          jwgt(:,i,j) = jwgt(:,im,jm)
          nwgt(  i,j) = nwgt(  im,jm)
          mwgt(  i,j) = mwgt(  im,jm)
           wgt(:,i,j) =  wgt(:,im,jm)
        enddo
      enddo
    case('T') ! T-folding
      do j = ny, ny - grdr % nrepeat + 1, -1
        jm = j - grdr % nrepeat - 1
        do i = grdr % werepeat, nx
          im = nx - i + 2
          iwgt(:,i,j) = iwgt(:,im,jm)
          jwgt(:,i,j) = jwgt(:,im,jm)
          nwgt(  i,j) = nwgt(  im,jm)
          mwgt(  i,j) = mwgt(  im,jm)
           wgt(:,i,j) =  wgt(:,im,jm)
        enddo
      enddo
      j = ny - grdr % nrepeat
        do i = nx / 2 + 1, nx
          im = nx - i + 2
          iwgt(:,i,j) = iwgt(:,im,j)
          jwgt(:,i,j) = jwgt(:,im,j)
          nwgt(  i,j) = nwgt(  im,j)
          mwgt(  i,j) = mwgt(  im,j)
           wgt(:,i,j) =  wgt(:,im,j)
        enddo
    end select
          
  END SUBROUTINE adjust_periodic

  SUBROUTINE intrp_itfc_nearst(ploc,pintr,psprd, Ms)

    ! Bi-linear interpolation (and possibly nearest) or bicubic

    IMPLICIT NONE
    ! arguments
    integer              :: F_igsrc,F_igdst
    type(type_location)  :: ploc
    type(type_weights) , intent(in), target :: psprd
    type(type_weights) , intent(inout), target :: pintr
    real(kind=8), dimension(nxsrc,nysrc), intent(in) :: Ms
    ! locals
    integer          :: i,j,k,is,js
    real(kind=8) :: as,bs,cs,ds
    real(kind=8), pointer, dimension(:,:) :: AA,BB
    integer     , pointer, dimension(:,:) :: II,JJ, nwgt
    real(kind=8), pointer, dimension(:,:,:) :: wgt
    integer     , pointer, dimension(:,:,:) :: iwgt, jwgt

      AA   => ploc % ag; BB   => ploc % bg
      II   => ploc % ig; JJ   => ploc % jg
       wgt => pintr %  wgt
      iwgt => pintr % iwgt
      jwgt => pintr % jwgt
      nwgt => pintr % nwgt

      k=1
      do j=1,nydst
      do i=1,nxdst
! ... if nearest neighbour corresponds to nint(I+A)
        if( ploc % mg(i,j) ) then
           is = II(i,j) + nint(AA(i,j))
           js = JJ(i,j) + nint(BB(i,j))
            wgt(k,i,j)  = oned
           iwgt(k,i,j)  = is
           jwgt(k,i,j)  = js
           nwgt(  i,j)  = 1
        endif
      enddo
      enddo

      ! adjust the weights according to spreading or mask
      call adjust_wgtb( pintr, psprd, Ms )

! ... Gestion du masque
      pintr % mwgt(:,:) = ploc % mg(:,:)
      if (Maskin) then
        where( pintr % nwgt(:,:) == 0 ) pintr % mwgt(:,:) = .false.
      endif

  END SUBROUTINE intrp_itfc_nearst

  SUBROUTINE intrp_itfc_bilc(ploc, pintr, psprd, Ms)

    ! Bi-linear interpolation (and possibly nearest) or bicubic

    IMPLICIT NONE
    ! arguments
    type(type_location), intent(in), target :: ploc
    type(type_weights) , intent(in), target :: psprd
    type(type_weights) , intent(inout), target :: pintr
    real(kind=8), dimension(nxsrc,nysrc), intent(in) :: Ms

    call cintrp( ploc, pintr, intyp=='bicub', nxsrc, nysrc)

    ! adjust the weights according to spreading or mask
    call adjust_wgtb( pintr, psprd, Ms )
    
! ... Gestion du masque
    pintr % mwgt(:,:) = ploc % mg(:,:)
    if (Maskin) then
        where( pintr % nwgt(:,:) == 0 ) pintr % mwgt(:,:) = .false.
    endif

  END SUBROUTINE intrp_itfc_bilc

  SUBROUTINE intrp_itfc_tri(ploc, pintr, grds, Ms, nearst)

    ! Pure aggregation

    IMPLICIT NONE
    ! arguments
    type(type_location), intent(in), target :: ploc
    type(type_grid)    , intent(in), target  :: grds
    type(type_weights) , intent(inout), target :: pintr
    real(kind=8), dimension(nxsrc,nysrc), intent(in) :: Ms
    logical, optional :: nearst
    ! locals
    integer         :: k,i,j,is,js,kmax(1), icell
    real(kind=8), pointer, dimension(:,:) :: AA,BB
    integer     , pointer, dimension(:,:) :: II,JJ, nwgt
    logical     , pointer, dimension(:,:) :: MI,ML
    real(kind=8), pointer, dimension(:,:,:) :: wgt
    integer     , pointer, dimension(:,:,:) :: iwgt, jwgt
    real(kind=8), allocatable, dimension(:,:) :: fout
    logical nearstw
    real(kind=8), dimension(3) :: wl

      AA   => ploc  % ag; BB   => ploc % bg
      II   => ploc  % ig; JJ   => ploc % jg
      ML   => ploc  % mg
      MI   => pintr % mwgt
       wgt => pintr %  wgt
      iwgt => pintr % iwgt
      jwgt => pintr % jwgt
      nwgt => pintr % nwgt

      nearstw = .false.
      if (present(nearst)) nearstw = nearst

      if (nearstw) then
         pintr % numwgt = 1
      else
         pintr % numwgt = 3
      endif

      write(*,*) 'CSTINTRP: UNSTRUCTURED INTERP'
      do j=1,nydst
      do i=1,nxdst
      if (ML(i,j)) then
         wl(3) = AA(i,j)
         wl(1) = BB(i,j)
         wl(2) = MIN( MAX(zerd, oned - wl(1) - wl(3)), oned)
         icell = ii(i,j)
         if (nearstw) then
           k=1
             wgt(k,i,j) = oned
            kmax = maxloc(wl)
            iwgt(k,i,j) = grds % meshd % in(kmax(k),icell)
            nwgt(i,j)   = k
         else
           do k = 1, grds % meshd % nconnect
             wgt(k,i,j) = wl(k)
            iwgt(k,i,j) = grds % meshd % in(k,icell)
            nwgt(i,j)   = k
           enddo
         endif
      endif
      enddo
      enddo

! ...   Gestion du masque
! ...   Application du seuil sur le masque
      MI(:,:) = ML(:,:)
      if (Maskin) then
        allocate( fout(nxdst, nydst) )
        call apply_weights( fout, Ms, pintr)
        where( fout(:,:) < mtrsh ) MI(:,:) = .false.
        deallocate(fout)
        where( .not.MI(:,:) ) nwgt(:,:) = 0
      endif

  END SUBROUTINE intrp_itfc_tri

  SUBROUTINE intrp_itfc_aggr(ploci, pintr, grds, grdr, Ms)

    ! Pure aggregation

    IMPLICIT NONE
    ! arguments
    type(type_location), intent(in), target :: ploci
    type(type_grid)    , intent(in), target  :: grds, grdr
    type(type_weights) , intent(inout), target :: pintr
    real(kind=8), dimension(nxsrc,nysrc), intent(in) :: Ms

    write(*,*) 'CSTINTRP: AGGREGATION OF THE CLOSEST SQUARED PATCH OF POINTS'
    call set_grid(grdr)
    call caggrt( ploci, grds, grdr, pintr, sdfrm, nadis, Ms ) ! MI also set up here

  END SUBROUTINE intrp_itfc_aggr

  SUBROUTINE intrp_itfc_mixtlc(plocr, ploci, pintr, psprd, grds, grdr, Ms )

    ! Mixt bi-linear (or bi-cubic) and aggregation

    IMPLICIT NONE
    ! arguments
    type(type_location), intent(in), target :: plocr, ploci
    type(type_weights) , intent(in), target :: psprd
    type(type_grid)    , intent(in), target  :: grds, grdr
    type(type_weights) , intent(inout), target :: pintr
    real(kind=8), dimension(nxsrc,nysrc), intent(in) :: Ms
    ! locals
    integer  nmixt, nw, i, j
    type(type_weights) :: bilic

    write(*,*) 'CSTINTRP: MIXT REGRIDDING PART 1 (FIRST FIELD) CAGGR'
    call set_grid(grdr)
    call caggrt( ploci, grds, grdr, pintr, sdfrm, nadis, Ms ) ! MI also set up here

    write(*,*) 'CSTINTRP: MIXT REGRIDDING PART 2 (FIRST FIELD) BILIC'
    bilic % numwgt = nwgtbmax ! set the max number of weights to a conservative number for accomodating spreading later
    bilic % nx = nxdst; bilic % ny = nydst
    call allocate_weights(bilic)
    call set_grid(grds)
    call cintrp( plocr, bilic, intyp=='bicub', nxsrc, nysrc)

    ! adjust the weights according to spreading
    bilic % numwgt = nwgtb ! set it back to a more reasonable number before spreading
    call adjust_wgtb( bilic, psprd, Ms )

    nmixt = nwgtb + 1
    write(*,*) 'Number of points to activate aggregation=',nmixt

    ! mix the aggregation weights with bilic weights
    if (naggrmax < bilic % numwgt) then
      write(*,*) 'increase naggrmax, it must be at least equal to nwgtbmax = ', bilic % numwgt
      call my_abort
    endif
    do j= 1, pintr % ny
    do i= 1, pintr % nx
      if( pintr % nwgt(i,j) < nmixt ) then ! bilinear or bicubic
        nw = bilic % nwgt(i,j)
        pintr % nwgt(i,j) = nw
        pintr %  wgt(1:nw,i,j) = bilic %  wgt(1:nw,i,j) 
        pintr % iwgt(1:nw,i,j) = bilic % iwgt(1:nw,i,j) 
        pintr % jwgt(1:nw,i,j) = bilic % jwgt(1:nw,i,j) 
        pintr %  wgt(nw+1:naggrmax,i,j) = zerd
        pintr % iwgt(nw+1:naggrmax,i,j) = 1
        pintr % jwgt(nw+1:naggrmax,i,j) = 1
      endif
    enddo
    enddo
    pintr % numwgt = maxval( pintr % nwgt )

    
! ...   Gestion du masque
    if (Maskin) then
        pintr % mwgt(:,:) = .true.
        where( pintr % nwgt(:,:) == 0 ) pintr % mwgt(:,:) = .false.
    else
        pintr % mwgt(:,:) = plocr % mg(:,:)
    endif

    ! deallocate working bilin/bicubic weights
    call deallocate_weights( bilic )

  END SUBROUTINE intrp_itfc_mixtlc

  SUBROUTINE intrp_itfc_scrip(plocr, ploci, pintr, grds, grdr, Ms, ln_scrip_grad)

    ! SCRIP (conservative interpolation)

    IMPLICIT NONE
    ! arguments
    type(type_location), intent(in), target :: plocr, ploci
    type(type_grid)    , intent(inout), target  :: grds, grdr
    type(type_weights) , intent(inout), target :: pintr
    real(kind=8), dimension(nxsrc,nysrc), intent(in) :: Ms
    logical,intent(in) :: ln_scrip_grad
    ! locals
    logical, pointer, dimension(:,:) :: mask
    logical, pointer, dimension(:,:) :: MI

    write(*,*) 'CSTINTRP: SCRIP INTERPOLATION'

    ! set the output mask before removal of redundant points
    MI   => pintr % mwgt
    MI(:,:) = grdr % lld % mask(:,:)

    ! define number of points to search
    n_search_scrip = nadis * nadis

    !
    ! create corner point arrays
    !
    if (grds % std_grd % grtyp /= 'M' ) call create_corner_mesh( grds )
    if (grdr % std_grd % grtyp /= 'M' ) call create_corner_mesh( grdr )
    !
    ! transfer variable mask to grid
    !
    allocate( grds % lld % mask( grds % nx, grds % ny ) )
    mask   => grds % lld % mask
    mask(:,:) = .false.
    where( Ms > 0.5_8) mask = .true.
    call scrip( plocr, ploci, grds, grdr, pintr, ln_scrip_grad, naggrmax )
    deallocate( grds % lld % mask )

    ! final adjustement to destination mask
    where( pintr % nwgt == 0 ) MI = .false.

  END SUBROUTINE intrp_itfc_scrip

  SUBROUTINE apply_weights( dtadst, dtasrc, pintr )
  ! apply weights
  implicit none
  ! arguments
  real(kind=8), dimension(:,:), intent(in) :: dtasrc
  real(kind=8), dimension(:,:), intent(out) :: dtadst
  type(type_weights), intent(in) :: pintr
  ! locals
  integer i,j,k,is,js

    dtadst(:,:) = zerd
    do j = 1, pintr % ny
    do i = 1, pintr % nx
    do k = 1, pintr % nwgt(i,j)
        is = pintr % iwgt(k,i,j)
        js = pintr % jwgt(k,i,j)
        dtadst(i,j) = dtadst(i,j) + pintr % wgt(k,i,j) *  dtasrc(is,js)
    enddo
    enddo
    enddo
  END SUBROUTINE apply_weights

  SUBROUTINE apply_weights_vector( dsrc, ddst, gsrc, gdst, pintr )
  ! apply weights
  implicit none
  ! arguments
  type(type_data) :: ddst, dsrc
  type(type_grid) :: gdst, gsrc
  type(type_weights), intent(in) :: pintr
  ! locals
  real(kind=8), dimension(3) :: v1,v2
  real(kind=8) :: vproj
  integer i,j,k,is,js

    do j = 1, pintr % ny
    do i = 1, pintr % nx
    v2(:) = zerd
    do k = 1, pintr % nwgt(i,j)
        is = pintr % iwgt(k,i,j)
        js = pintr % jwgt(k,i,j)
        v1(:) = dsrc % u(is,js) * gsrc % lld % di3d(:,is,js) &
              + dsrc % v(is,js) * gsrc % lld % dj3d(:,is,js)
        v2(:) = v2(:) + pintr % wgt(k,i,j) * v1(:)
    enddo
    ddst % u(i,j) = dot_product( v2, gdst % lld % di3d(:,i,j) )
    ddst % v(i,j) = dot_product( v2, gdst % lld % dj3d(:,i,j) )
    enddo
    enddo

  END SUBROUTINE apply_weights_vector

  SUBROUTINE intrp_getput_premap_wgt(F_fullcal,F_mode,pintr,fwgt)

    ! Read/write pre mapped weights
    ! According to logic defined by interpolation method
    ! WARNING APPLY SAME LOGIC AS IN intrp_itfc

    ! Arguments
    IMPLICIT NONE
    logical, intent(out)         :: F_fullcal
    character(len=3), intent(in) :: F_mode
    type(type_weights), intent(inout) :: pintr
    character(len=*), intent(in) :: fwgt

    ! Local variables
    integer            :: unit
    logical :: file_exists

    inquire(file=fwgt, exist=file_exists)

    if ( F_mode == 'get' .and. file_exists) then
      ! read weight file if it exists
      write(*,*) 'CSTINTRP: TRYING TO OPEN: ', trim(fwgt)
      unit=get_unit()
      if (unit < 0) then
        write(*,*) 'intrp_getput_premap_wgt: could not get unit number to open old '//trim(fwgt)
        call my_abort
      endif
      open(unit,file=fwgt,form='unformatted',status='old')
      write(*,*) 'CSTINTRP: READING wgt FROM: ', trim(fwgt)
      read(unit) pintr % numwgt
      read(unit) pintr % nwgt(:, :)
      read(unit) pintr %  wgt(1:pintr % numwgt,:,:) 
      read(unit) pintr %  iwgt(1:pintr % numwgt,:,:)
      read(unit) pintr %  jwgt(1:pintr % numwgt,:,:)
      read(unit) pintr % mwgt(:,:)
      close(unit)
      F_fullcal=.false.
    endif

    if ( F_mode == 'put' ) then
      ! write weight file if does not exist
      if (.not. file_exists) then
        write(*,*) 'CSTINTRP: TRYING TO OPEN (write): ', trim(fwgt)
        unit=get_unit()
        if (unit < 0) then
          write(*,*) 'intrp_getput_premap_wgt: could not get unit number to open new '//trim(fwgt)
          call my_abort
        endif
        open(unit,file=fwgt,form='unformatted',status='new')
        write(*,*) 'CSTINTRP: WRITING wgt IN: ', trim(fwgt)
        write(unit) pintr % numwgt
        write(unit) pintr % nwgt(:, :)
        write(unit) pintr %  wgt(1:pintr % numwgt,:,:) 
        write(unit) pintr %  iwgt(1:pintr % numwgt,:,:)
        write(unit) pintr %  jwgt(1:pintr % numwgt,:,:)
        write(unit) pintr % mwgt(:,:)
        close(unit)
      else
        write(*,*)  "CSTINTRP: wont't write weights to disk, already exist - ", trim(fwgt)
      endif
    endif

  END SUBROUTINE intrp_getput_premap_wgt

  SUBROUTINE intrp_itfc_write(var,grid,F_mask_change)
    !! -----------------------------------------------------------------------------------------------
    !!  *** Write interpolated records and weights to destination file
    !! -----------------------------------------------------------------------------------------------
    ! * Arguments

    IMPLICIT NONE
    type (variable),  intent(in) :: var
    type (type_grid), intent(in) :: grid
    logical         , intent(in) :: F_mask_change

    ! * Local variable
    type (std_variable)  :: F_vars,F_varm,S_vard,V_vard, MS_vard, MV_vard
    type (std_grid)      :: F_grdr
    integer              :: ierr, status_st90, i,j,k,kk, igsrc, ngsrc_out, igdst
    logical              :: Maskouw, use_Mdstw_o
    integer, allocatable, &
      dimension(:,:,:,:) :: Mdstw_o
    character(len=1)     :: canum
    real(kind=8), dimension(nxdst,nydst,ngdst,ngsrc) :: fdst1, fdst2
    real(kind=8), dimension(nxdst,nydst) :: fdst
    integer     , dimension(nxdst,nydst) :: idst

    F_vars = var % varstd
    V_vard = var % varstdv
    F_varm = var % mask
    F_grdr = grid % std_grd

! ... Define main field attributes

    Maskouw = (Maskin .or. Maskou)
    use_Mdstw_o = .false.

    S_vard = F_vars  !Scalar or U component
    S_vard % grtyp = F_grdr%grtyp
    S_vard % ig1   = F_grdr%ig1
    S_vard % ig2   = F_grdr%ig2
    S_vard % ig3   = F_grdr%ig3
    S_vard % ig4   = F_grdr%ig4
    S_vard % ni    = F_grdr%ni    !nxdst
    S_vard % nj    = F_grdr%nj    !nydst if not YY
    if ( npak_ou /= 999 ) &
            S_vard % nbits = -npak_ou
    if ( Maskouw .and. .not.ln_mask ) &
            S_vard % typvar(2:2)='@'   !source var had no mask

    ngsrc_out=ngsrc
    if ( .not.splt_yy .and. ngsrc == 2 ) then
      ngsrc_out=1
    endif

    ! mask treatment
    if ( Maskouw .and. .not.splt_yy .and. ngsrc == 2 ) then  ! combine src mask
        allocate(Mdstw_o(nxdst,nydst,ngdst,ngsrc))
        Mdstw_o(:,:,:,:) = Mdstw(:,:,:,:)
        use_Mdstw_o = .true.

        where ( Mdstw(:,:,:,1) == 0 .and. Mdstw(:,:,:,2) > 0 )
          Mdstw_o(:,:,:,1) = Mdstw(:,:,:,2)
        endwhere

        Maskouw = .false.
        if ( any( Mdstw_o(:,:,:,1) == 0 )) Maskouw = .true.
        if ( .not. Maskouw ) S_vard%typvar(2:2)=' '
    endif

    ! transfer to flat field
    do igsrc = 1, ngsrc
    do igdst = 1, ngdst
       fdst1(:,:, igdst, igsrc) = dtadst( igdst, igsrc) % u (:,:)
    enddo
    enddo
    if (vector) then
      do igsrc = 1, ngsrc
      do igdst = 1, ngdst
        fdst2(:,:, igdst, igsrc) = dtadst( igdst, igsrc) % v (:,:)
      enddo
      enddo
    endif

    if ( .not.splt_yy .and. ngsrc == 2 ) then
      do igdst = 1, ngdst
        if ( .not. use_Mdstw_o ) then ! no mask
          call combine_src_yy_datad(fdst1(:,:,igdst,:), nxdst, nydst)
          if (vector) &
            call combine_src_yy_datad( &
                                    fdst2(:,:,igdst,:), nxdst, nydst)
        else ! mask
          call combine_src_yy_mask_datad( &
                                    fdst1(:,:,igdst,:), Mdstw(:,:,igdst,:), nxdst, nydst)
          if (vector) &
            call combine_src_yy_mask_datad( &
                                    fdst2(:,:,igdst,:), Mdstw(:,:,igdst,:), nxdst, nydst)
        endif ! mask
      enddo
    endif ! combine source domain

    if ( Maskouw .and. zland ) then
      do igsrc = 1 , ngsrc_out
        do igdst = 1, ngdst
          if ( use_Mdstw_o ) then
              call get_mean_val_and_apply(nxdst, nydst, fdst1(:,:,igdst,igsrc), Mdstw_o(:,:,igdst,igsrc), landval)
            if ( vector ) &
              call get_mean_val_and_apply(nxdst, nydst, fdst2(:,:,igdst,igsrc), Mdstw_o(:,:,igdst,igsrc), landval)
          else
              call get_mean_val_and_apply(nxdst, nydst, fdst1(:,:,igdst,igsrc), Mdstw(  :,:,igdst,igsrc), landval)
            if ( vector ) &
              call get_mean_val_and_apply(nxdst, nydst, fdst2(:,:,igdst,igsrc), Mdstw(  :,:,igdst,igsrc), landval)
          endif
        enddo
      enddo
    endif

    if ( nbdyf > 0 ) then
      do igsrc = 1 , ngsrc_out
        do igdst = 1, ngdst
          do k=nbdyf,1,-1
            i=k+1
            do j=1,nydst
                fdst1(  i-1,j,igdst,igsrc)=fdst1(  i,j,igdst,igsrc)
              if ( vector ) &
                fdst2(  i-1,j,igdst,igsrc)=fdst2(  i,j,igdst,igsrc)
              if ( use_Mdstw_o ) then
                Mdstw_o(i-1,j,igdst,igsrc)=Mdstw_o(i,j,igdst,igsrc)
               else
                Mdstw(  i-1,j,igdst,igsrc)=Mdstw(  i,j,igdst,igsrc)
              endif
            enddo
            i=nxdst-k
            do j=1,nydst
                fdst1(  i+1,j,igdst,igsrc)=fdst1(  i,j,igdst,igsrc)
              if ( vector ) &
                fdst2(  i+1,j,igdst,igsrc)=fdst2(  i,j,igdst,igsrc)
              if ( use_Mdstw_o ) then
                Mdstw_o(i+1,j,igdst,igsrc)=Mdstw_o(i,j,igdst,igsrc)
               else
                Mdstw(  i+1,j,igdst,igsrc)=Mdstw(  i,j,igdst,igsrc)
              endif
            enddo
            j=k+1
            do i=1,nxdst
                fdst1(  i,j-1,igdst,igsrc)=fdst1(  i,j,igdst,igsrc)
              if ( vector ) &
                fdst2(  i,j-1,igdst,igsrc)=fdst2(  i,j,igdst,igsrc)
              if ( use_Mdstw_o ) then
                Mdstw_o(i,j-1,igdst,igsrc)=Mdstw_o(i,j,igdst,igsrc)
               else
                Mdstw(  i,j-1,igdst,igsrc)=Mdstw(  i,j,igdst,igsrc)
              endif
            enddo
            j=nydst-k
            do i=1,nxdst
                fdst1(  i,j+1,igdst,igsrc)=fdst1(  i,j,igdst,igsrc)
              if ( vector ) &
                fdst2(  i,j+1,igdst,igsrc)=fdst2(  i,j,igdst,igsrc)
              if ( use_Mdstw_o ) then
                Mdstw_o(i,j+1,igdst,igsrc)=Mdstw_o(i,j,igdst,igsrc)
               else
                Mdstw(  i,j+1,igdst,igsrc)=Mdstw(  i,j,igdst,igsrc)
              endif
            enddo
          enddo
        enddo
      enddo
    endif

    if ( Maskouw .or. owgts ) then
      if ( ln_mask ) then     !Scalar or U component mask
        MS_vard = F_varm
       else
        MS_vard = F_vars
        MS_vard%typvar='@@'
        MS_vard%datyp=DATYP_UNSIGNED_INT
        MS_vard%nbits=1
        MS_vard%xtype=1
      endif
      MS_vard%grtyp = S_vard%grtyp
      MS_vard%ig1   = F_grdr%ig1
      MS_vard%ig2   = F_grdr%ig2
      MS_vard%ig3   = F_grdr%ig3
      MS_vard%ig4   = F_grdr%ig4
      MS_vard%ni    = S_vard%ni
      MS_vard%nj    = S_vard%nj
      MS_vard%etiket= S_vard%etiket
    endif

    ! main field
    do igsrc = 1 , ngsrc_out
       write(canum,'(i1)') igsrc
       if ( splt_yy .and. ngsrc_out == 2) then
          S_vard % etiket(7:12) = 'Src'//canum//' '
          S_vard % ip3 = igsrc - 1
       endif
       
       ! handle integer case
       if (any(S_vard%datyp == DATYP_INT_LIST)) then
        call write_var_std(ngdst, subgdst, S_vard, unitd, i2d=nint(fdst1(:,:,:,igsrc)))
       else
        call write_var_std(ngdst, subgdst, S_vard, unitd, d2d=fdst1(:,:,:,igsrc))
       endif 

       if ( Maskouw ) then
         if ( splt_yy .and. ngsrc_out == 2) then
            MS_vard % etiket(7:12) = 'Src'//canum//' '
            MS_vard % ip3 = igsrc - 1
         endif
         if ( use_Mdstw_o ) then
            call write_var_std(ngdst, subgdst, MS_vard, unitd, i2d=Mdstw_o(:,:,:,igsrc))
         else
            call write_var_std(ngdst, subgdst, MS_vard, unitd, i2d=Mdstw(:,:,:,igsrc))
         endif
       endif
    enddo

    if ( vector ) then
      V_vard = S_vard  ! V component
      V_vard % nomvar = nomvar_v
      ! main field V
      do igsrc = 1 , ngsrc_out
       S_vard % typvar='P@'
       S_vard % nbits=32
       S_vard % datyp=DATYP_IEEE_FLOATING_POINT
       S_vard % xtype=2
       write(canum,'(i1)') igsrc
       if ( splt_yy .and. ngsrc_out == 2) then
          S_vard % etiket(10:10) = canum
          S_vard % ip3 = igsrc - 1
       endif
            call write_var_std(ngdst, subgdst, V_vard, unitd, d2d=fdst2(:,:,:,igsrc))
       if ( Maskouw ) then
         MV_vard = MS_vard
         MV_vard % nomvar = nomvar_v ! V component mask
         if ( splt_yy .and. ngsrc_out == 2) then
            MS_vard % etiket(10:10) = canum
            MS_vard % ip3 = igsrc - 1
         endif
         if ( use_Mdstw_o ) then
            call write_var_std(ngdst, subgdst, MV_vard, unitd, i2d=Mdstw_o(:,:,:,igsrc))
         else
            call write_var_std(ngdst, subgdst, MV_vard, unitd, i2d=Mdstw(:,:,:,igsrc))
         endif
       endif
      enddo
    endif

! Weight outputs can not be combined for now ...
! Need to figure out how GEM-NEMO communications 
! under YY will be set up

    if (owgts.and.F_mask_change) then 
      ! if owgts, F_mask_change is equivalent to fullcal
      S_vard%nbits=0
      S_vard%typvar='P'
      S_vard%etiket='WEIGHT'
      S_vard%ip1=0
      S_vard%ip2=0
      S_vard%ip3=0
      if (ngdst == 2 .and. .not.splt_yy) S_vard % nj = S_vard % nj / 2
      if (Maskin.or.Maskou) S_vard%typvar(2:2)='@'
      !
      !-----------------------------------------------------------
      ! write weights
      !-----------------------------------------------------------
      ! split subdomains assumed active (looking for trouble?) // FD cannot use new routine because YY can differ in number of weights
        do igsrc = 1 , ngsrc
          if ( ngsrc == 2 ) then
               write(canum,'(i1)') igsrc
               S_vard % etiket(7:10) = 'Src'//canum
          endif
          S_vard % ip3 = igsrc - 1
          do igdst = 1, ngdst 
            S_vard % ig3 = subgdst(igdst) % ig3
            if ( ngdst == 2 ) then
               if ( igdst == 1 ) S_vard%etiket(1:6)='YinDst'
               if ( igdst == 2 ) S_vard%etiket(1:6)='YanDst'
            endif
            ! output number of weights per point
              S_vard%nomvar='NAVG'
              S_vard%typvar='P@'
              S_vard%datyp=DATYP_UNSIGNED_INT
              S_vard%nbits=16
              S_vard%xtype=1
              idst(:,:) = intr(igdst,igsrc) % nwgt(:,:)
              ierr = putstdvar(unitd, idst, nxdst, nydst, S_vard, ierr_st90=status_st90)

            do k = 1, intr(igdst,igsrc) % numwgt
              fdst(:,:) = intr(igdst,igsrc) % wgt(k,:,:)
              write(S_vard%nomvar(2:4),'(i3.3)') k
              S_vard%nomvar(1:1)='W'
              S_vard%typvar='P'
              S_vard%datyp=5
              S_vard%nbits=64
              S_vard%xtype=3
              ierr = putstdvar(unitd, fdst, nxdst, nydst, S_vard, ierr_st90=status_st90)
              if ( ierr /= 0 ) then
                write(*,*) 'CSTINTRP putstdvar problem => wgt(:,:,k,igdst)', k
                write(*,*) 'CSTINTRP ierr,S_vard=',ierr,S_vard
                write(*,*) 'CSTINTRP status_st90=',status_st90
                call my_abort
              endif
              S_vard%datyp=DATYP_UNSIGNED_INT
              S_vard%nbits=16
              S_vard%xtype=1
              S_vard%nomvar(1:1)='I'
              idst(:,:) = intr(igdst,igsrc) % iwgt(k,:,:)
              ierr = putstdvar(unitd, idst, S_vard%ni, S_vard%nj, S_vard, ierr_st90=status_st90)
              if ( ierr /= 0 ) then
                write(*,*) 'CSTINTRP putstdvar problem => iwgt(:,:,k,igdst)', k
                write(*,*) 'CSTINTRP ierr,S_vard=',ierr,S_vard
                write(*,*) 'CSTINTRP status_st90=',status_st90
                call my_abort
              endif
              S_vard%nomvar(1:1)='J'
              idst(:,:) = intr(igdst,igsrc) % jwgt(k,:,:)
              ierr = putstdvar(unitd, idst, S_vard%ni, S_vard%nj, S_vard, ierr_st90=status_st90)
              if ( ierr /= 0 ) then
                write(*,*) 'CSTINTRP putstdvar problem => jwgt(:,:,k,igdst)', k
                write(*,*) 'CSTINTRP ierr,S_vard=',ierr,S_vard
                write(*,*) 'CSTINTRP status_st90=',status_st90
                call my_abort
              endif
            enddo
          enddo
        enddo

      ! write mask
      MS_vard%nomvar='MASK'
      MS_vard%ip1=0
      MS_vard%ip2=0
      MS_vard%ip3=0
      MS_vard%etiket='WEIGHT'
      do igsrc = 1 , ngsrc
        write(canum,'(i1)') igsrc
        if ( splt_yy .and. ngsrc == 2) then
          MS_vard % etiket(7:10) = 'Src'//canum
          MS_vard % ip3 = igsrc - 1
        endif
        MS_vard % ip3 = igsrc - 1
        call write_var_std(ngdst, subgdst, MS_vard, unitd, i2d=Mdstw(:,:,:,igsrc))
      enddo
    endif ! condition on outputting weights info

    if ( use_Mdstw_o ) deallocate(Mdstw_o)
    do igsrc = 1, ngsrc
    do igdst = 1, ngdst
       deallocate( dtadst(igdst,igsrc) % u )
      if (vector) &
       deallocate( dtadst(igdst,igsrc) % v )
    enddo
       deallocate( dtasrc(igsrc) % u )
      if (vector) &
       deallocate( dtasrc(igsrc) % v )
      if (ln_mask) &
       deallocate( dtasrc(igsrc) % m )
    enddo

  END SUBROUTINE intrp_itfc_write

  SUBROUTINE write_var_std(ng, subg, varp, unitf, i2d, f2d, d2d)
    !! -----------------------------------------------------------------------------------------------
    !!  *** write to file data for separate YY output
    !! -----------------------------------------------------------------------------------------------
    IMPLICIT NONE
    ! arguments
    integer, intent(in) :: ng
    type(std_grid), dimension(ng), intent(in) :: subg
    type(std_variable), intent(inout)  :: varp
    integer,      dimension(subg(1) % ni, subg(1) % nj, ng), optional, intent(in) :: i2d
    real(kind=4), dimension(subg(1) % ni, subg(1) % nj, ng), optional, intent(in) :: f2d
    real(kind=8), dimension(subg(1) % ni, subg(1) % nj, ng), optional, intent(in) :: d2d
    integer unitf
    ! locals
    integer ig, ngout
    integer ierr, status_st90
    real(kind=8), allocatable, &
      dimension(:,:)     :: wrk8
    real(kind=4), allocatable, &
      dimension(:,:)     :: wrk4
    integer, allocatable, &
      dimension(:,:)     :: int4
    logical lcomb ! combine the two sub domains in one for output
    integer nx, ny

    if ( ng == 2 .and. splt_yy ) then
      varp % grtyp  = subg(1) % grtyp
      varp % nj     = subg(1) % nj
    endif

    nx = varp % ni
    ny = varp % nj

    lcomb = .false.
    ngout = ng
    if ( ng == 2 .and. .not.splt_yy ) then ! combine the 2 subdomains in one
       lcomb = .true.
       ngout = 1
    endif

    do ig = 1, ngout

    if ( .not.lcomb ) then ! regular output, just copy grid info from subgrid (copy of main if ng=1)
         varp % ig1 = subg(ig) % ig1
         varp % ig2 = subg(ig) % ig2
         varp % ig3 = subg(ig) % ig3
         varp % ig4 = subg(ig) % ig4
      endif

      if ( ng == 2 .and. splt_yy ) then
        select case (ig)
        case(1)
          varp % etiket(1:6) = 'YinDst'
        case(2)
          varp % etiket(1:6) = 'YanDst'
        end select
      endif

      if (present(i2d)) then ! integer case

        if ( lcomb ) then ! combined
           allocate(int4(nx,ny))
           int4(:,1:ny/2     ) = i2d(:,:,1)
           int4(:,  ny/2+1:ny) = i2d(:,:,2)
        else
           allocate(int4(nx,ny))
           int4(:,:) = i2d(:,:,ig)
        endif
        ierr = putstdvar(unitd, int4, nx, ny, varp, ierr_st90=status_st90)
        deallocate( int4 )

      endif

      if (present(f2d)) then ! float case

        if ( lcomb ) then ! combined
           allocate(wrk4(nx,ny))
           wrk4(:,1:ny/2     ) = f2d(:,:,1)
           wrk4(:,  ny/2+1:ny) = f2d(:,:,2)
        else
           allocate(wrk4(nx,ny))
           wrk4(:,:) = f2d(:,:,ig)
        endif

        ierr = putstdvar(unitd, wrk4, nx, ny, varp, ierr_st90=status_st90)
        deallocate( wrk4 )

      endif

      if (present(d2d)) then

        if ( lcomb ) then ! combined
           allocate(wrk8(nx,ny))
           wrk8(:,1:ny/2     ) = d2d(:,:,1)
           wrk8(:,  ny/2+1:ny) = d2d(:,:,2)
        else
           allocate(wrk8(nx,ny))
           wrk8(:,:) = d2d(:,:,ig)
        endif

        ierr = putstdvar(unitd, wrk8, nx, ny, varp, ierr_st90=status_st90)
        deallocate( wrk8 )

      endif

      if ( ierr /= 0 ) then
        write(*,*) 'CSTINTRP putstdvar problem => Adst*180/pi'
        write(*,*) 'CSTINTRP ierr,varr_p(1)=',ierr,varp
        write(*,*) 'CSTINTRP status_st90=',status_st90
        call my_abort
      endif

    enddo

  END SUBROUTINE write_var_std

  SUBROUTINE combine_src_yy_datar(dta,nx,ny)
    IMPLICIT NONE
    ! arguments
    integer nx,ny
    real(Kind=4) :: dta(nx,ny,2)
    ! locals
    integer i,j

    do j=1,nydst
       do i=1,nxdst
          dta(i,j,1)=0.5*(dta(i,j,1)+dta(i,j,2))
       enddo
    enddo

  END SUBROUTINE combine_src_yy_datar

  SUBROUTINE combine_src_yy_datad(dta,nx,ny)
    IMPLICIT NONE
    ! arguments
    integer nx,ny
    real(Kind=8) :: dta(nx,ny,2)

    dta(:,:,1)=0.5_8*(dta(:,:,1)+dta(:,:,2))

  END SUBROUTINE combine_src_yy_datad

  SUBROUTINE combine_src_yy_mask_datad(dta,M,nx,ny)
    IMPLICIT NONE
    ! arguments
    integer nx,ny
    real(Kind=8) :: dta(nx,ny,2)
    integer      :: M  (nx,ny,2)
    ! locals

    where ( M(:,:,1) > 0 .and. M(:,:,2) > 0 ) &
              dta(:,:,1) = 0.5_8 * ( dta(:,:,1) + dta(:,:,2) )

    where ( M(:,:,1)== 0 .and. M(:,:,2) > 0 ) &
              dta(:,:,1) =                        dta(:,:,2)
  END SUBROUTINE combine_src_yy_mask_datad


  SUBROUTINE get_mean_val_and_apply(nx, ny, dta, M, lval)
    IMPLICIT NONE
    ! arguments
    integer nx,ny
    real(kind=8) :: dta(nx,ny)
    integer      :: M  (nx,ny)
    real(kind=4) :: lval
    ! locals
    integer i, j, cnt
    real(kind=8) mval

    cnt   = 0
    mval  = zerd

    if (lval == 999.9) then
       do j=1,ny
          do i=1,nx
             if (M(i,j)==1) then
                cnt = cnt + 1
                mval = mval + dta(i,j)
             endif
          enddo
       enddo
       if ( cnt > 0 ) mval = mval / real( cnt, kind = 8 )
    else
       mval = real( lval, kind=8 )
    endif

    ! apply
    dta(:,:) = dta(:,:) * M(:,:) + mval * (oned - M(:,:) )

  END SUBROUTINE get_mean_val_and_apply



  SUBROUTINE mask_tartinage(F_grds)
    IMPLICIT NONE
    
    type (type_grid), target, intent(in) :: F_grds

    ! local
    integer :: igsrc, nsprdw
    type (type_grid), pointer :: subs
    type(type_weights), pointer :: psprd


    nsprdw = 0
    nwgts = 1
    if ( ln_mask .and. F_grds % std_grd % grtyp /= 'M' ) then
      select case(intyp)
      case( 'nearst', 'bilin', 'mixt', 'bicub', 'mixtbc' )
       if (nsprd > 0) then
         nwgts = 4 * nsprd
         nsprdw = nsprd
       endif
      end select
    endif
    call reallocate_sprd( ngsrc )
    do igsrc=1,ngsrc
      subs => F_grds % subgrd(igsrc)
      psprd => sprd(igsrc)
      psprd % numwgt = nwgts
      if (nsprdw > 0) then
        call mapospread( subs, psprd, Msrc(:,:,igsrc), nsprd, sdfrm, intyp=='nearst' )
        write(*,*) 'spread'
        write(*,*) 'min/max W',minval( psprd %  wgt ),maxval( psprd %  wgt )
        write(*,*) 'min/max I',minval( psprd % iwgt ),maxval( psprd % iwgt )
        write(*,*) 'min/max J',minval( psprd % jwgt ),maxval( psprd % jwgt )
      endif
    enddo

  END SUBROUTINE mask_tartinage

END MODULE intrp_itfc
