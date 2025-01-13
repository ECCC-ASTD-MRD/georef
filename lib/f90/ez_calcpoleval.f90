      subroutine ez_calcpoleval(poleval, z, ni, ax, grtyp, grref)
      implicit none
      
      integer ni
      real poleval,z(ni)
      real(kind=8) ax(ni)
      character grtyp, grref
      
      integer i
      
      if (grtyp.eq.'Z'.and.grref.eq.'E') then
         poleval = 0.0
         do i=1,ni-1
            poleval = poleval + z(i)*(ax(i+1)-ax(i))
         enddo
         poleval = poleval / 360.0
         return
      endif

      poleval = 0.0
      do i=1,ni
         poleval = poleval + z(i)
      enddo
      poleval = poleval / (1.0 * ni)

      return
      end

