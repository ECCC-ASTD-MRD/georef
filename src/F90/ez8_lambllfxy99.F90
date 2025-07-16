subroutine ez8_lambllfxy99(lat, lon, x, y, n, latin1, latin2, yaxislat, yaxislon)
    use iso_fortran_env, only: real64
    implicit none

    integer, intent(in) :: n
    real(kind = real64), intent(out) :: lat(n), lon(n)
    real(kind = real64), intent(in) :: x(n), y(n)
    real, intent(in) :: latin1, latin2, yaxislat, yaxislon

    real :: pi, pisur4, d2r, rn, rphi1, rphi2, r, f, rtan, rho, theta, rhozero
    integer :: i

    pisur4 = atan(1.0)
    pi     = 4.0 * pisur4
    d2r    = pi / 180.0
    r = 6370997.0

    rphi1 = d2r * latin1
    rphi2 = d2r * latin2

    if (rphi1 == rphi2) then
        rn = sin(rphi1)
    else
        rn = log(cos(rphi1) / cos(rphi2))
        rn = rn / log((tan(pisur4 + 0.5 * rphi2)) / tan(pisur4 + 0.5 * rphi1))
    endif

    rtan = tan(pisur4 + rphi1 * 0.5)
    f = (cos(rphi1) * (rtan ** rn)) / rn
    rhozero = r * f / (tan(pisur4 + yaxislat * d2r * 0.5) ** rn)

    do i = 1, n
        rho = sign(1.0, rn) * sqrt(x(i) * x(i) + ((rhozero - y(i)) ** 2))
        theta = atan(x(i) / (rhozero - y(i)))

        lat(i) = (2.0 * atan((r * f / rho) ** (1.0 / rn)) - 0.5 * pi) / d2r
        lon(i) = theta / (d2r * rn) + yaxislon
    enddo
end