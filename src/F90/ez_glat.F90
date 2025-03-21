!> Compute the latitudes of gaussian grid
subroutine ez_glat(latroots, groots, nj, hem)
    use iso_fortran_env
    implicit none

    integer, intent(in) :: nj
    real(kind = real64), dimension(nj), intent(out) :: latroots
    real, intent(inout) :: groots(*)
    integer, intent(in) :: hem

#include "qqqpar.inc"
#include "pi.inc"

    external dgauss

    integer :: j, npoly
    real :: temp

    if (hem /= GLOBAL) then
        npoly = nj * 2
    else
        npoly = nj
    endif
    call dgauss(npoly, groots, global)

    do j = 1, npoly / 2
        temp = groots(j)
        groots(j) = groots(npoly + 1 - j)
        groots(npoly + 1 - j) = temp
    enddo

    if (hem /= nord) then
        do j = 1, nj
            latroots(j) = 90.0 - rdtodg * acos(groots(j))
        enddo

    endif

    if (hem == NORD) then
        do j = 1, nj
            latroots(j) = 90.0 - rdtodg * acos(groots(j + nj))
        end do
    endif
end
