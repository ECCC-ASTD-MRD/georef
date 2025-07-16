subroutine ez8_lambfll(x, y, xlat, xlon, npts, grtyp, ig1, ig2, ig3, ig4)
    use iso_fortran_env, only: real64
    implicit none

    integer, intent(in) :: npts
    real(kind = real64), intent(out) :: x(npts), y(npts)
    real(kind = real64), intent(in) :: xlat(npts), xlon(npts)
    character, intent(in) :: grtyp
    integer, intent(in) :: ig1, ig2, ig3, ig4

    real :: xg(15)
    character :: gtypout
    real :: xlatninj, xlonninj, dx, dy, latin1, latin2
    real :: yaxislat, yaxislon
    real(kind = real64) xlat11, xlon11, x11, y11

    integer i, nxg

    integer lclig1, lclig2, lclig3, lclig4

    lclig1 = ig1
    lclig2 = ig2
    lclig3 = ig3
    lclig4 = ig4

    nxg = 15
    if (grtyp == '!') then
        call igaxg95(gtypout, xg, 15, grtyp, lclig1, lclig2, lclig3, lclig4)
        if (gtypout == 'H') then
            xlat11 =   xg( 1)
            xlon11 =   xg( 2)
            xlatninj = xg(10)
            xlonninj = xg(11)
            yaxislat = real(0.5 * (xlat11 + xlatninj))
            yaxislon = xg(5)
            latin1   = xg(6)
            latin2   = xg(7)
            dx       = xg(3) * 1000.0
            dy       = xg(4) * 1000.0

            call ez8_lambxyfll99(x11, y11, xlat11, xlon11, 1, latin1, latin2, yaxislat, yaxislon)

            call ez8_lambxyfll99(x, y, xlat, xlon, npts, latin1, latin2, yaxislat, yaxislon)
            do i = 1, npts
                x(i) = 1.0 + (x(i) - x11) / dx
                y(i) = 1.0 + (y(i) - y11) / dy
            enddo
        endif
    endif
end
