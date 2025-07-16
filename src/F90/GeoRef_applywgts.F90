subroutine GeoRef_applywgts(outfld, wts, idxs, infld, masque, ni_src, nj_src, ni_dst, nj_dst, n_wts, nodata)
    use iso_fortran_env, only: real64
    implicit none

    integer, intent(in) :: ni_src, nj_src, ni_dst, nj_dst, n_wts
    real, intent(out) :: outfld(ni_dst * nj_dst)
    real(kind = real64), intent(in) :: wts(ni_dst, nj_dst, n_wts)
    real, intent(in) :: infld(ni_src * nj_src)
    integer, intent(in) :: idxs(ni_dst, nj_dst, n_wts)
    integer, intent(in) :: masque(ni_dst * nj_dst)
    real, intent(in) :: nodata
    integer :: i, j, k, n

    outfld = nodata

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k, i, j, n) SHARED(ni_dst, nj_dst, masque, outfld, wts, n_wts, infld, idxs)
    do k = 1, ni_dst * nj_dst
        i = mod((k - 1), ni_dst) + 1
        j = 1 + k / ni_dst
        if (masque(k) == 1) then
            outfld(k) = 0.0
            do n = 1, n_wts
                if (idxs(i, j, n) < 1) stop
                outfld(k) = real(outfld(k) + wts(i, j, n) * infld(idxs(i, j, n)))
            enddo
        endif
    enddo
end subroutine
