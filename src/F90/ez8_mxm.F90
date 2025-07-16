subroutine ez8_mxm(a, nar, b, nac, c, nbc)
    use iso_fortran_env, only: real64
    implicit none

    integer, intent(in) :: nar, nac, nbc
    real, intent(in) :: a(nar, 1)
    real(kind = real64), intent(in) :: b(nac, 1)
    real(kind = real64), intent(out) :: c(nar, 1)

    integer :: i, j, k

    do j=1, nbc
        do i = 1, nar
            c(i, j) = 0.0
            do k = 1, nac
                c(i, j) = c(i, j) + a(i, k) * b(k, j)
            end do
        end do
    end do
end
