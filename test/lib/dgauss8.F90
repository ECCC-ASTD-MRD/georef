program test_dgauss8
  use iso_fortran_env
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer, parameter :: nroots = 90
  integer :: i
  real(kind = real64), dimension(nroots) :: roots

  call dgauss8(nroots, roots, 0)

!   print *, "Nodes (x):"
!   do i = 1, nroots
!     print *, roots(i)
!   end do

  if (.not. all(ieee_is_finite(roots))) then
    print *, "Error: Bad values found in roots."
    error stop 1
  else
    print *, "All roots computed successfully."
  end if


end program test_dgauss8
