      subroutine ez_xpngdb2(zout,zi,ni,nj,j1,j2,hem,symetrie)
      
      implicit none
#include "qqqpar.cdk"

      external permut
      
      integer j1,j2,ni,nj
      integer hem,symetrie
      integer i,j
      
      
      real zout(ni,j1:j2)
      real zi(ni,nj)
      real sign
      
      if (symetrie.eq.0) then
         sign = -1.0
      else
         sign = 1.0
      endif
      
      if (hem .eq. nord) then
         do j=1,nj
            do i=1,ni
               zout(i,j)  = zi(i,j)
            enddo
         enddo

         do j=2,nj
            do i=1,ni
               zout(i,2-j)  = sign * zi(i,j)
            enddo
         enddo
      endif

      if (hem .eq. sud) then
         do j=1,nj
            do i=1,ni
               zout(i,j)  = zi(i,j)
            enddo
         enddo
         
         do j=2,nj
            do i=1,ni
               zout(i,nj+j-1)  = sign * zi(i,nj-j+1)
            enddo
         enddo
      endif
      
         
      return
      end
      
