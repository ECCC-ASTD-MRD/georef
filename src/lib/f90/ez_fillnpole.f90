      subroutine ez_fillnpole(zout, zin, ni, j1, j2, valpole)
      implicit none

      integer ni,j,i,j1, j2
      real zout(ni, 4), zin(ni,j1:j2), valpole

      do j=1,3
         do i=1,ni
            zout(i,j) = zin(i,j2-3+j)
         enddo
      enddo
      
      do i=1,ni
         zout(i,4) = valpole
      enddo
         
      return
      end
      
