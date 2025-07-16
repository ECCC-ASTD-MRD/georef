!> Convert from grid wind components to std meteorological speed and direction
subroutine ez8_llwfgdw(z1, z2, xlon, li, lj, grtyp, ig1, ig2, ig3, ig4)
    use iso_fortran_env, only: real64
    use GeoRef_internal_def, only: rdtodg
    implicit none

    !> X dimension of the grid
    integer, intent(in) :: li
    !> Y dimension of the grid
    integer, intent(in) :: lj
    !> On entry, contains the U grid component
    real, intent(inout) :: z1(li, lj)
    !> On entry, contains the V component
    real, intent(inout) :: z2(li, lj)
    !> Longitudes
    real(kind = real64), intent(in) :: xlon(li, lj)
    !> Grid type
    character, intent(in) :: grtyp
    !> grid descriptors
    integer, intent(in) :: ig1, ig2, ig3, ig4

    !> Converts grid wind components to meteorological wind speed and direction, independent of geographical projection.
    !> On output, a 270 degree wind is a westerly wind, meaning U is +ve and V is zero.

   external cigaxg

    real :: xg1, xg2, xg3, xg4
    integer :: i, j
    real(kind = real64) :: spd0, dir0
    real(kind = real64) :: x1(2 * li * lj), y1(2 * li * lj), lat(2 * li * lj)

!  les #define qui suivent rendent le code plus lisible
#define uu   z1
#define vv   z2
#define spd  z1
#define dir  z2

    if (grtyp .eq. '!') then
        call ez8_lamb_llwfgdw(uu, vv, xlon, li, lj, grtyp, ig1, ig2, ig3, ig4, x1, y1, lat)
        return
    endif

    if (grtyp .eq. 'N')then
        call cigaxg(grtyp, xg1, xg2, xg3, xg4, ig1, ig2, ig3, ig4)
        do j = 1, lj
            do i = 1, li
                spd0 = sqrt(uu(i, j) * uu(i, j) + vv(i, j) * vv(i, j))
                if (spd0 .eq. 0.0) then
                    dir0 = 0.0
                else
                    if (uu(i, j) .eq. 0.0)then
                        if (vv(i, j) .ge. 0.0)then
                            dir0 = xlon(i, j) + xg4 - 90.0
                        else
                            dir0 = xlon(i, j) + xg4 + 90.0
                        endif
                    else
                        dir0 = xlon(i, j) + xg4 - rdtodg * atan2(vv(i, j), uu(i, j))
                    endif
                endif
                dir0 = mod(mod(dir0, 360.0d0) + 360.0, 360.0d0)
                spd(i, j) = real(spd0)
                dir(i, j) = real(dir0)
            enddo
        enddo
        return
    endif

    if (grtyp .eq. 'S') then
        call cigaxg(grtyp, xg1, xg2, xg3, xg4, ig1, ig2, ig3, ig4)
        do j = 1, lj
            do i = 1, li
                spd0 = sqrt(uu(i, j) * uu(i, j) + vv(i, j) * vv(i, j))
                if (spd0 .eq. 0.0) then
                    dir0 = 0.0
                else
                    if (uu(i, j) .eq. 0.0)then
                        if (vv(i, j) .ge. 0.0)then
                            dir0 = 90.0 - xlon(i, j) + xg4
                        else
                            dir0 = 270.0 - xlon(i, j) + xg4
                        endif
                    else
                        dir0 = 180.0 - xlon(i, j) + xg4 - rdtodg * atan2(vv(i, j), uu(i, j))
                    endif
                endif
                dir0 = mod(mod(dir0, 360.0d0) + 360.0, 360.0d0)
                spd(i, j) = real(spd0)
                dir(i, j) = real(dir0)
            enddo
        enddo
        return
    endif

    if (grtyp.eq.'A'.or.grtyp.eq.'B'.or.grtyp.eq.'G'.or.grtyp.eq.'L'.or.grtyp.eq.'M') then
        do j = 1, lj
            do i = 1, li
                spd0 = sqrt(uu(i, j) * uu(i, j) + vv(i, j) * vv(i, j))
                if (spd0 .eq. 0.0) then
                    dir0 = 0.0
                else
                    if (uu(i, j) .eq. 0.0)then
                        if (vv(i, j) .ge. 0.0)then
                            dir0 = 180.0
                        else
                            dir0 = 0.0
                        endif
                    else
                        dir0=270.0 - rdtodg*atan2(vv(i, j), uu(i, j))
                    endif
                endif
                dir0 = mod(mod(dir0, 360.0d0) + 360.0, 360.0d0)
                spd(i, j) = real(spd0)
                dir(i, j) = real(dir0)
            enddo
        enddo
        return
    endif
    write(6, "('0', ' erreur, bad grid type (llwfgdw) - grtyp = ', a11)") grtyp
end
