      subroutine ez_fillspole(zout, zin, ni, j1, j2, valpole)
      implicit none

      integer i,j,ni,j1,j2
      real zout(ni, 4),zin(ni,j1:j2),valpole
      
      do j=j1,j1+2
         do i=1,ni
            zout(i,j-j1+2) = zin(i,j)
         enddo
      enddo
      
      do i=1,ni
         zout(i,1) = valpole
      enddo
         

      return
      end
      
