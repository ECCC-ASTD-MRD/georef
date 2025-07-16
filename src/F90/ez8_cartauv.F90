!> Compute the components of the winds in the rotated system of coordinates from the winds in the rotated cartesian space
subroutine ez8_cartauv(u, v, uvcart, lon, lat, ni, nj)
    use iso_fortran_env, only: real64
    implicit none

    !> E-W dimension of the grid
    integer, intent(in) :: ni
    !> N-S dimension of the grid
    integer, intent(in) :: nj
    !> Rotated winds in cartesian space
    real(kind = real64), intent(in) :: uvcart(3, ni * nj)
    !> Rotated component u of the wind
    real, intent(out) :: u(ni, nj)
    !> Rotated component v of the wind
    real, intent(out) :: v(ni, nj)
    !> Longitudes of the grid in the rotated system of coordinates
    real(kind = real64), intent(in) :: lon(ni, nj)
    !> Latitudes of the grid in the rotated system of coordinates
    real(kind = real64), intent(in) :: lat(ni, nj)

    !> \ingroup ezscint

    real(kind = real64), parameter :: dar = acos(-1.0) / 180.0

    integer :: i, j, k
    real(kind = real64) :: a, b, c, d, e, f

    k = 0
    !$omp parallel do default(none) private(i, j, a, b, c, d, e, f) firstprivate(k) shared(ni, nj, dar, lat, lon, uvcart, u, v)
    do j = 1, nj
        do i = 1, ni
            k      = k + 1
            a      = cos(dar * lon(i, j))
            b      = sin(dar * lon(i, j))
            e      = cos(dar * lat(i, j))
            f      = sin(dar * lat(i, j))
            u(i, j) = (uvcart(2, k) * a) - (uvcart(1, k) * b)
            c      = (uvcart(1, k) * a) + (uvcart(2, k) * b)
            d      = sqrt(c**2 + uvcart(3, k)**2 )
            v(i, j) = sign(real(d), real((uvcart(3, k) * e) - (c * f)))
        end do
    end do
end