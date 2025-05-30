subroutine ez8_rgdint_0(zo, px, py, npts, z, ni, j1, j2, nodata)
    use iso_fortran_env
    implicit none

    integer, intent(in) :: npts, ni, j1, j2
    real, intent(out) :: zo(npts)
    real(kind = real64), intent(in) :: px(npts), py(npts)
    real, intent(in) :: z(ni, j1:j2), nodata

    integer :: i, j, n

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n, i, j) SHARED(npts, j1, j2, ni, px, py, z, zo)
    do n = 1, npts
        i = min(ni, max(1, nint(px(n))))
        j = min(j2, max(j1, nint(py(n))))

        zo(n) = z(i, j)
    enddo
end


subroutine ez8_rgd_index_0(idx, px, py, npts, ni, j1, j2)
    use iso_fortran_env
    implicit none

    integer, intent(in) :: npts, ni, j1, j2
    real(kind = real32), intent(out) :: idx(npts, 2)
    real(kind = real64), intent(in) :: px(npts), py(npts)

    integer :: n

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n) SHARED(npts, idx, ni, j1, j2, px, py)
    do n = 1, npts
        idx(n, 1) = min(ni, max(1, nint(px(n))))
        idx(n, 2) = min(j2, max(j1, nint(py(n))))
    enddo
end


subroutine ez8_apply_0(idx, zo, npts, z, ni, j1, j2, nodata)
    use iso_fortran_env
    implicit none

    integer, intent(in) :: npts, ni, j1, j2
    real(kind = real32), intent(in) :: idx(npts, 2)
    real, intent(in) :: nodata
    real, intent(out) :: zo(npts)
    real, intent(in) :: z(ni, j1:j2)

    integer :: i, j, n

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n, i, j) SHARED(npts, idx, z, zo)
    do n = 1, npts
        i = int(idx(n, 1))
        j = int(idx(n, 2))

        zo(n) = z(i, j)
    enddo
end
