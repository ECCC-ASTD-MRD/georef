! RMNLIB - Library of useful routines for C and FORTRAN programming
! Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
!                          Environnement Canada
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the
! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! Boston, MA 02111-1307, USA.


!> \file

!> Subset of belousovs algorithm used to calculate ordinary legendre polynomials
subroutine ordleg8(sx, coa, ir)
    use iso_fortran_env, only : real64
    implicit none

    !> Legendre polynomial evaluated at coa
    real(kind = real64), intent(out) :: sx
    !> Cosine of colatitude
    real(kind = real64), intent(in) :: coa
    !> Wave number
    integer, intent(in) :: ir

    integer :: irpp, irppm, k, kk, n, n1
    real(kind = real64) :: a, ang, b, c1, c4, delta, fk, fn, fn2, fn2sq, s1, sia, sqr2, theta

    sqr2 = sqrt(2.0)
    irpp = ir + 1
    irppm = irpp - 1
    delta = acos(coa)
    sia = sin(delta)

    theta = delta
    c1 = sqr2

    do n = 1, irppm
      fn = real(n, real64)
      fn2 = 2.0 * fn
      fn2sq = fn2 * fn2
      c1 = c1 * sqrt(1.0 - 1.0 / fn2sq)
    end do

    n = irppm
    ang = fn * theta
    s1 = 0.0
    c4 = 1.0
    a = -1.0
    b = 0.0
    n1 = n + 1

    do kk = 1, n1, 2
        k = kk - 1
        if (k == n) c4 = 0.5 * c4
        s1 = s1 + c4 * cos(ang)
        a = a + 2.0
        b = b + 1.0
        fk = real(k, real64)
        ang = theta * (fn - fk - 2.0)
        c4 = (a * (fn - b + 1.0) / (b * (fn2 - a))) * c4
    end do

    sx = s1 * c1
end


!> Calculates the zeroes of the ordinary legendre polynomial of order n,  i.e. define gaussian grid
subroutine dgauss8(n, roots, kase)
    use iso_fortran_env, only : real64
    use rmn_base_const, only: pie, south
    implicit none

    !> Order of the polynomials
    integer, intent(in) :: n
    !> Zeroes of the ordinary Legendre polynomials
    real(kind = real64), intent(out) :: roots(*)
    !> Area: 0 = global, 1 = North, 2 = South
    integer, intent(in) :: kase

    !> The positive roots are approximated by the best asymptotic formula available to the author, found in
    !> abramowitz and stegun "handbook of mathematical functions".
    !> chapter 22 formula 22.16.6.
    !> Newton's method is used to refine the guess to precision defined by the constant tol.  since the roots are of order
    !> of magnitude unity, absolute precision is adequate, rather than a relative test.
    !> A standard identity is used to determine the derivative of the polynomial in terms of the values of p(n;x), p(n-1;x).
    !> (x**2-1.0)*(dp/dx)=n*(x*p(n;, x)-p(n-1;x)).
    !> See abramowitz and stegun formula 22.8.5
    !> Note that in contrast to other formulas this requires only 2 evaluations of a legendre polynomial per iteration.
    !> Note that the coordinate used is conventionally referred to as mu=cos(theta), running from +1 to -1, for theta from 0 to
    !> pi. The negative roots are  filled by symmetry. for kase=global, all n roots are found, while for
    !> dase=north/south only the +ve/-ve roots are found, (including 0 if n is odd)  i.e. n/2+mod(n, 2) roots.

    real(kind = real64), parameter :: tol = 1.0E-13

    integer :: i, j, l, irt
    real(kind = real64) :: normn, normnm
    real(kind = real64) :: delta, g, gm, pn, pnm, rdpdx, t

    !  ordleg8 returns polynomials normalized to unit integral.
    !  normn, normnmn restore the convention normalization, p(n;1.0)=1.0.
    normn = sqrt(2.0 / (2.0 * n + 1.0))
    normnm = sqrt(2.0 / (2.0 * n - 1.0))
    l = n / 2

    ! calculate asymptotic approximation
    do i = 1, l
        if (kase /= south) j = i
        if (kase == south) j = i + l + mod(n, 2)
        t = (4 * j - 1) * pie / float(4 * n + 2)
        if (kase /= south) irt = i
        if (kase == south) irt = i + mod(n, 2)
        roots(irt) = cos(t + 1.0 / (8.0 * float(n**2) * tan(t)))
    end do

    do i = 1, l
        ! Repeat 1 newton iteration
        delta = huge(0.0)
        do while (abs(delta) > tol)
            call ordleg8(g, roots(i), n)
            call ordleg8(gm, roots(i), n - 1)
            pn = normn * g
            pnm = normnm * gm
            rdpdx = (roots(i) ** 2 - 1.0) / (n * (roots(i) * pn - pnm))
            delta = -pn * rdpdx
            roots(i) = roots(i) + delta
        end do

        roots(n + 1 - i) = -roots(i)
    end do

    if (mod(n, 2) == 0) return
    if (kase /= south) irt = l + 1
    if (kase == south) irt = 1
    roots(irt) = 0.0
end
