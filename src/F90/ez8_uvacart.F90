!> Compute the winds in the cartesian space from the components
subroutine ez8_uvacart(xyz, u, v, lon, lat, ni, nj)
    use iso_fortran_env, only: real64
    implicit none

    !> E-W dimension of the grid
    integer, intent(in) :: ni
    !> N-S dimension of the grid
    integer, intent(in) :: nj
    !> Unrotated U component of the wind
    real, intent(in) :: u(ni, nj)
    !> Unrotated V component of the wind
    real, intent(in) :: v(ni, nj)
    !> Unrotated winds in cartesian space
    real(kind = real64), intent(out) :: xyz(3, ni * nj)
    !> Grid longitudes in the unrotated system of coordinates
    real(kind = real64), intent(in) :: lon(ni, nj)
    !> Grid latitudes in the unrotated system of coordinates
    real(kind = real64), intent(in) :: lat(ni, nj)

    !> \ingroup ezscint

    integer :: i, j, k
    real(kind = real64) :: a, b, c, d, dar

    dar = acos(-1.0) / 180.0
    k = 0

    !$omp parallel do default(none) private(i, j, a, b, c, d) firstprivate(k) shared(ni, nj, dar, lat, lon, u, v, xyz)
    do j = 1, nj
        do i = 1, ni
            k = k + 1
            a = sin(dar * lon(i, j))
            b = cos(dar * lon(i, j))
            c = sin(dar * lat(i, j))
            d = cos(dar * lat(i, j))
            xyz(1, k) = real( -(u(i, j) * a) - (v(i, j) * b * c) )
            xyz(2, k) = real( (u(i, j) * b) - (v(i, j) * a * c) )
            xyz(3, k) = real( v(i, j) * d )
        end do
    end do
end