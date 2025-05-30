subroutine ez_xpngdb2(zout, zi, ni, nj, j1, j2, hem, symetrie)
    implicit none
#include "qqqpar.inc"

    integer, intent(in) :: j1
    integer, intent(in) :: j2
    integer, intent(in) :: ni
    integer, intent(in) :: nj
    integer, intent(in) :: hem
    integer, intent(in) :: symetrie

    real, intent(out) :: zout(ni, j1:j2)
    real, intent(in) :: zi(ni, nj)

    integer :: i, j
    real :: sign

    if (symetrie == 0) then
        sign = -1.0
    else
        sign = 1.0
    endif

    if (hem == nord) then
        do j = 1, nj
            do i = 1, ni
                zout(i, j)  = zi(i, j)
            enddo
        enddo

        do j = 2, nj
            do i = 1, ni
                zout(i, 2-j)  = sign * zi(i, j)
            enddo
        enddo
    endif

    if (hem == sud) then
        do j = 1, nj
            do i = 1, ni
                zout(i, j)  = zi(i, j)
            enddo
        enddo

        do j = 2, nj
            do i = 1, ni
                zout(i, nj + j - 1)  = sign * zi(i, nj - j + 1)
            enddo
        enddo
    endif
end
