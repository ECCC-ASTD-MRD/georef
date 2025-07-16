!> Compute the latitudes of gaussian grid
subroutine ez8_glat(latroots, groots, nj, hem)
    use iso_fortran_env
    use rmn_base_const, only: global, north, rdtodg
    implicit none

    integer, intent(in) :: nj
    real(kind = real64), dimension(nj), intent(out) :: latroots
    real(kind = real64), intent(inout) :: groots(*)
    integer, intent(in) :: hem

    external :: dgauss8

    integer :: j, npoly
    real(kind = real64) :: temp

    if (hem /= global) then
        npoly = nj * 2
    else
        npoly = nj
    endif
    call dgauss8(npoly, groots, global)

    do j = 1, npoly / 2
        temp = groots(j)
        groots(j) = groots(npoly + 1 - j)
        groots(npoly + 1 - j) = temp
    enddo

    if (hem /= north) then
        do j = 1, nj
            latroots(j) = 90.0 - rdtodg * acos(groots(j))
        enddo

    endif

    if (hem == north) then
        do j = 1, nj
            latroots(j) = 90.0 - rdtodg * acos(groots(j + nj))
        end do
    endif
end
