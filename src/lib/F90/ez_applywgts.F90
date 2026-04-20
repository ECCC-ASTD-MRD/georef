   subroutine ez_applywgts(outfld, wts, idxs, infld, masque, ni_src, nj_src, ni_dst, nj_dst, n_wts, nodata)
   use, intrinsic :: iso_fortran_env
   implicit none
   integer ni_src, nj_src, ni_dst, nj_dst, n_wts, n
   real :: infld(ni_src*nj_src), outfld(ni_dst*nj_dst), nodata
   integer i,j,k

   real(kind=real64) :: wts(ni_dst, nj_dst, n_wts)
   integer :: idxs(ni_dst, nj_dst, n_wts), masque(ni_dst*nj_dst)

   outfld = nodata

   !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,i,j) SHARED(ni_dst,nj_dst,masque,outfld,wts,n_wts,infld,idxs)
   do k=1,ni_dst*nj_dst
     i = mod((k-1),ni_dst) + 1
     j = 1+k/ni_dst
     if (masque(k) == 1) then
        outfld(k) = 0.0
        do n=1,n_wts
           if (idxs(i,j,n) < 1) stop
           outfld(k) = outfld(k)+wts(i,j,n)*infld(idxs(i,j,n))
        enddo
     endif
   enddo

   end subroutine ez_applywgts
