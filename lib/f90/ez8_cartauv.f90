      SUBROUTINE ez8_CARTAUV(U, V, UVCART, LON, LAT, NI, NJ)
      implicit none
      INTEGER NI, NJ 
      REAL    U(NI,NJ), V(NI,NJ)
      real(kind=8)  UVCART(3,NI*NJ), LON(NI,nj), LAT(ni,NJ)
!
!author michel roch - april 90
!
!arguments
!    out    U       - rotated component  U of the wind
!           V       - rotated component  V of the wind
!    in     xyz     - rotated winds in cartesian space
!           LON     - longitudes of the grid in the rotated system of coordinates
!           LAT     - latitudes of the grid in the rotated system of coordinates
!           NI      - E-W DIMENSION of the grid
!           NJ      - N-S dimension of the grid
!
!*

      INTEGER I, J, K 
      real(kind=8)    A, B, C, D, E, F, DAR

      DAR = ACOS(-1.)/180.
      K   = 0
 
      !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,A,B,C,D,E,F) FIRSTPRIVATE(K) SHARED(NI,NJ,DAR,LAT,LON,UVCART,U,V)
      DO 20 J=1,NJ
         DO 10 I=1,NI
            K      = K+1
            A      = COS(DAR*LON(I,j))
            B      = SIN(DAR*LON(I,j))
            E      = COS(DAR*LAT(i,J))
            F      = SIN(DAR*LAT(i,J))
            U(I,J) = (UVCART(2,K)*A) - (UVCART(1,K)*B)
            C      = (UVCART(1,K)*A) + (UVCART(2,K)*B)
            D      = SQRT(C**2 + UVCART(3,K)**2 )
            V(I,J) = SIGN(D, (UVCART(3,K)*E)-(C*F))
 10         CONTINUE
 20      CONTINUE

      RETURN
      END 
