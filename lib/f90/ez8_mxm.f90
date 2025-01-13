subroutine ez8_mxm(a,nar,b,nac,c,nbc)

      integer nar,nac,nbc,i,j,k
      real a(nar,1)
      real(kind=8) b(nac,1),c(nar,1)

      do 30 j=1,nbc
            do 20 i = 1,nar
               c(i,j) = 0.0
               do 10 k = 1,nac
                  c(i,j) = c(i,j) + a(i,k)*b(k,j)
      10       continue
      20    continue
      30 continue
      
      return 
      end
