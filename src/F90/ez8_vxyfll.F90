!> Compute the grid co-ordinates measured from the pole of a point, given the latitude and longitude in degrees
subroutine ez8_vxyfll(x, y, dlat, dlon, npts, d60, dgrw, pi, pj, nhem)
    use iso_fortran_env, only: real64
    use GeoRef_internal_def, only: dgtord, nord, sud
    implicit none

    !> Number of elements in the x, y, dlat, dlon arrays
    integer, intent(in) :: npts
    !> Hemisphere : 1 for north, 2 for south
    integer, intent(in) :: nhem
    !> x-coordinate of the point32_t as measured with pole as origin
    real(kind = real64), intent(out) :: x(npts)
    !> y-coordinate of the point32_t as measured with pole as origin
    real(kind = real64), intent(out) :: y(npts)
    !> Latitude in degrees (-90 to +90, positive n)
    real(kind = real64), intent(in) :: dlat(npts)
    !> Longitude in degrees (-180 to +180, positive e)
    real(kind = real64), intent(in) :: dlon(npts)
    !> Grid length (in metres) of the polar stereographic grid at 60 degrees
    real, intent(in) :: d60
    !> Orientation of greenwich meridian with respect to the grid (in degrees)
    real, intent(in) :: dgrw
    !> 
    real, intent(in) :: pi
    !> 
    real, intent(in) :: pj

    !> \note The companion routine llfxy, which computes the latitude and longitude given the grid-coordinates, is also available

    real(kind = real64) :: re, rlon, rlat, sinlat, r
    integer :: i

    re = 1.866025d0 * 6.371d+6 / d60
    if (nhem == nord) then
        !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, rlat, rlon, sinlat, r) SHARED(npts, x, y, dgtord, re, dlat, dlon, pi, pj, dgrw)
        do i = 1, npts
            rlon = dgtord * (dlon(i) + dgrw)
            rlat = dgtord * dlat(i)
            sinlat = sin(rlat)
            r = re * sqrt((1.d0 - sinlat) / (1.d0 + sinlat))
            x(i) = r * cos(rlon) + pi
            y(i) = r * sin(rlon) + pj
        enddo
    elseif (nhem == sud) then
        !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, rlat, rlon, sinlat, r) SHARED(npts, x, y, dgtord, re, dlat, dlon, pi, pj, dgrw)
        do i = 1, npts
            rlon = dlon(i)
            if (rlon > 180.0d0) rlon = rlon - 360.0d0
            rlon = dgtord * (-rlon + dgrw)
            rlat = dgtord * (-dlat(i))
            sinlat = sin(rlat)
            r = re * sqrt((1.d0 - sinlat) / (1.d0 + sinlat))
            x(i) = r * cos(rlon) + pi
            y(i) = r * sin(rlon) + pj
         enddo
    endif
end
