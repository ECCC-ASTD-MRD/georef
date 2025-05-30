!> Computes the grid co-ordinates of a point
subroutine ez8_vllfxy(dlat, dlon, x, y, ni, nj, d60, dgrw, pi, pj, nhem)
    use iso_fortran_env
    implicit none

    !> \note The companion routine xyfll, which computes the grid co-ordinates given the latitude and longitude, is also available.

    !> First dimension of x, y, dlat, dlon arrays
    integer, intent(in) :: ni
    !> Second dimension of x, y, dlat, dlon arrays
    integer, intent(in) :: nj
    !> Hemisphere : 1 for north, 2 for south
    integer, intent(in) :: nhem
    !> x-co-ordinate of the point32_t as measured with pole as origin
    real(kind = real64), intent(in) :: x(ni, nj)
    !> y-co-ordinate of the point32_t as measured with pole as origin
    real(kind = real64), intent(in) :: y(ni, nj)
    !> Latitude in degrees (-90 to +90, positive n)
    real(kind = real64), intent(out) :: dlat(ni, nj)
    !> Longitude in degrees (-180 to +180, positive e)
    real(kind = real64), intent(out) :: dlon(ni, nj)
    !> Grid length (in metres) of the polar stereographic grid at 60 degrees
    real, intent(in) :: d60
    !> Orientation of greenwich meridian with respect to the grid (in degrees)
    real, intent(in) :: dgrw
    real, intent(in) :: pi
    real, intent(in) :: pj

#include "pi.inc"
#include "qqqpar.inc"

    real(kind = real64) :: x1, y1
    real(kind = real64) :: re, re2, r2, rlat, rlon
    integer :: i, j

    re = 1.866025d0 * 6.371e+6 / d60
    re2 = re**2

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, x1, y1, rlat, rlon, r2) SHARED(ni, nj, x, y, rdtodg, re, re2, nhem, dlat, dlon, pi, pj, dgrw)
    do j = 1, nj
        do i = 1, ni
            x1 = x(i, j) - pi
            y1 = y(i, j) - pj

            if (x1 == 0.  .and. y1 == 0.)  then
                rlat = 90.0d0
                rlon = 0.0d0
            endif

            ! calculate longitude in map coordinates
            if (x1 == 0.0d0) rlon = sign(90.d0, y1)
            if (x1 /= 0.0d0) rlon = atan(y1 / x1) * dble(rdtodg)
            if (x1 < 0.0d0) rlon = rlon + sign(180.d0, y1)

            ! adjust longitude for grid orientation

            rlon = rlon-dgrw
            if (rlon < 0.0d0) rlon = rlon + 3.6d2

            ! calculate latitude
            r2 = x1 * x1 + y1 * y1
            rlat = (re2 - r2) / (re2 + r2)
            rlat = max(-1.0d0, min(rlat, 1.0d0))
            rlat = asin(rlat) * dble(rdtodg)

            ! change signs if in southern hemisphere
            if (nhem == 2) then
                rlat = -rlat
                rlon = -rlon
            endif
            if (rlon < 0.0d0) rlon = rlon + 360.0d0
            dlat(i, j) = rlat
            dlon(i, j) = rlon
        enddo
    enddo
end
