!*********************************************************************
!**s/r ez_glat - calcul des latitudes d'une grille gaussienne
!
!  auteur: Yves Chartier. Mars 1991.
!****
!
      subroutine ez_glat(latroots,groots,nj,hem)
      implicit none

      integer nj, hem

#include "qqqpar.cdk"
#include "pi.cdk"
      
      external dgauss
      real(kind=8) latroots(*)
      real groots(*)
    
      integer j,npoly
      real temp

      if (hem .ne. GLOBAL) then
         npoly = nj * 2
      else
         npoly = nj
      endif
      call dgauss(npoly, groots, global)
      
      do j=1,npoly/2
         temp = groots(j)
         groots(j)=groots(npoly+1-j)
         groots(npoly+1-j)=temp
      enddo
      
      if (hem .ne. nord) then
         do j=1,nj
            latroots(j)= 90. - rdtodg * acos(groots(j))
         enddo
         
      endif
      
      if (hem .eq. NORD) then
         do 20 j=1,nj
            latroots(j)=90.0-rdtodg*acos(groots(j+nj))
 20      continue
      endif

      return
      end
