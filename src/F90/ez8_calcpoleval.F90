subroutine ez8_calcpoleval(poleval, z, ni, ax, grtyp, grref)
    use iso_fortran_env, only: real64
    implicit none

    integer, intent(in) :: ni
    real, intent(out) :: poleval
    real, intent(in) :: z(ni)
    real(kind = real64), intent(in) :: ax(ni)
    character, intent(in) :: grtyp, grref

    integer :: i

    if (grtyp == 'Z'.and.grref == 'E') then
        poleval = 0.0
        do i = 1, ni - 1
            poleval = real(poleval + z(i) * (ax(i + 1) - ax(i)))
        enddo
        poleval = poleval / 360.0
        return
    endif

    poleval = 0.0
    do i = 1, ni
        poleval = poleval + z(i)
    enddo
    poleval = poleval / (1.0 * ni)
end
